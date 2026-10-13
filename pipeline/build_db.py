"""Load the cached abstract embeddings into a local LanceDB table.

Does not call an embedding API. Vectors already on disk are normalized the
same way the locked exam normalized them at scoring time. No approximate
index is built, so a search is an exact comparison of the filtered rows.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import numpy as np
import pyarrow as pa

from pipeline.paths import LOCAL, embed_cache

DB_PATH = LOCAL / "ukb" / "pubmed" / "weeks_abstracts.sqlite"
OUT = LOCAL / "ukb" / "final_exam" / "paper_vectors"
MANIFEST = OUT / "manifest.json"
DIM = 1536


def schema() -> pa.Schema:
    return pa.schema(
        [
            pa.field("pmid", pa.int64()),
            pa.field("year", pa.int32()),
            pa.field("month", pa.int32()),
            pa.field("day", pa.int32()),
            pa.field("before_cutoff", pa.int32()),
            pa.field("gene_symbols", pa.list_(pa.string())),
            pa.field("title", pa.string()),
            pa.field("abstract", pa.string()),
            pa.field("vector", pa.list_(pa.float32(), DIM)),
        ]
    )


def load_vectors() -> dict[int, np.ndarray]:
    found: dict[int, np.ndarray] = {}
    paths = sorted(embed_cache().glob("chunk_*.npz"))
    for index, path in enumerate(paths, start=1):
        data = np.load(path)
        pmids = data["pmids"].astype(np.int64)
        vectors = data["vectors"].astype(np.float32)
        for pmid, vector in zip(pmids, vectors):
            found[int(pmid)] = np.array(vector, copy=True)
        if index % 500 == 0 or index == len(paths):
            print(f"CHUNKS {index} of {len(paths)} unique {len(found)}", flush=True)
    return found


def load_meta(pmids: list[int]) -> dict[int, tuple]:
    conn = sqlite3.connect(DB_PATH)
    found = {}
    for start in range(0, len(pmids), 400):
        batch = pmids[start : start + 400]
        marks = ",".join("?" for _ in batch)
        sql = f"""
            SELECT pmid, year, month, day, before_cutoff, title, abstract
            FROM abstracts
            WHERE pmid IN ({marks})
              AND abstract IS NOT NULL
              AND abstract != ''
        """
        for pmid, year, month, day, before, title, abstract in conn.execute(sql, batch):
            found[int(pmid)] = (
                int(year or 0),
                int(month or 0),
                int(day or 0),
                int(before or 0),
                title or "",
                abstract or "",
            )
    genes: dict[int, list[str]] = {pmid: [] for pmid in found}
    for start in range(0, len(pmids), 400):
        batch = [pmid for pmid in pmids[start : start + 400] if pmid in found]
        if not batch:
            continue
        marks = ",".join("?" for _ in batch)
        for pmid, symbol in conn.execute(
            f"SELECT pmid, symbol FROM gene_hits WHERE pmid IN ({marks})",
            batch,
        ):
            bucket = genes.get(int(pmid))
            if bucket is not None:
                bucket.append(symbol)
    conn.close()
    for pmid in genes:
        genes[pmid] = sorted(set(genes[pmid]))
    return {pmid: meta + (genes[pmid],) for pmid, meta in found.items()}


def main() -> None:
    import lancedb

    if MANIFEST.exists():
        saved = json.loads(MANIFEST.read_text(encoding="utf-8"))
        print(f"ALREADY_BUILT rows {saved.get('rows')}", flush=True)
        return
    vectors = load_vectors()
    print(f"VECTORS {len(vectors)}", flush=True)
    meta = load_meta(list(vectors))
    print(f"WITH_TEXT {len(meta)}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    database = lancedb.connect(str(OUT))
    table = None
    batch_pmids = []
    written = 0

    def flush() -> None:
        nonlocal table, batch_pmids, written
        if not batch_pmids:
            return
        years, months, days, befores, titles, abstracts, symbols = [], [], [], [], [], [], []
        matrix = np.vstack([vectors[pmid] for pmid in batch_pmids]).astype(np.float32)
        norms = np.linalg.norm(matrix, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        matrix /= norms
        for pmid in batch_pmids:
            year, month, day, before, title, abstract, gene_list = meta[pmid]
            years.append(year)
            months.append(month)
            days.append(day)
            befores.append(before)
            titles.append(title)
            abstracts.append(abstract)
            symbols.append(gene_list)
        arrow = pa.table(
            {
                "pmid": pa.array(batch_pmids, type=pa.int64()),
                "year": pa.array(years, type=pa.int32()),
                "month": pa.array(months, type=pa.int32()),
                "day": pa.array(days, type=pa.int32()),
                "before_cutoff": pa.array(befores, type=pa.int32()),
                "gene_symbols": pa.array(symbols, type=pa.list_(pa.string())),
                "title": pa.array(titles, type=pa.string()),
                "abstract": pa.array(abstracts, type=pa.string()),
                "vector": pa.FixedSizeListArray.from_arrays(
                    pa.array(matrix.reshape(-1), type=pa.float32()),
                    DIM,
                ),
            },
            schema=schema(),
        )
        if table is None:
            table = database.create_table("papers", arrow, mode="overwrite")
        else:
            table.add(arrow)
        written += len(batch_pmids)
        batch_pmids = []
        print(f"WROTE {written}", flush=True)

    for pmid in meta:
        batch_pmids.append(pmid)
        if len(batch_pmids) >= 2000:
            flush()
    flush()
    MANIFEST.write_text(
        json.dumps({"rows": written, "dim": DIM, "metric": "l2_on_normalized_vectors", "index": "none"}, indent=2),
        encoding="utf-8",
    )
    print(f"DONE {written}", flush=True)


if __name__ == "__main__":
    main()

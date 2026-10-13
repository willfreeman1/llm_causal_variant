"""LanceDB vector search with gene and date prefilters.

Abstracts that name the lead variant are excluded in that same filter, before
the top 25 are chosen. There is no approximate index. Ties break by PubMed id,
ascending. On neighborhood A00407 the same 25 papers come back as in the
locked run, with two swapped, because those scores are equal to 7 decimals.
"""

from __future__ import annotations

import sqlite3

import numpy as np

from pipeline.paths import LOCAL

DB_DIR = LOCAL / "ukb" / "final_exam" / "paper_vectors"
ABSTRACTS = LOCAL / "ukb" / "pubmed" / "weeks_abstracts.sqlite"
TOP_N = 25


def _escape(text: str) -> str:
    return text.replace("'", "''")


def date_sql(cutoff: str, date_mode: str) -> str:
    if date_mode == "october_2019":
        return "before_cutoff = 1"
    year, month, day = (int(part) for part in cutoff.split("-"))
    return (
        "year > 0 AND ("
        f"(month <= 0 AND year < {year}) OR "
        f"(month > 0 AND day <= 0 AND (year < {year} OR (year = {year} AND month < {month}))) OR "
        f"(day > 0 AND (year < {year} OR (year = {year} AND month < {month}) OR "
        f"(year = {year} AND month = {month} AND day < {day})))"
        ")"
    )


def filter_sql(genes: list[str], cutoff: str, rsid: str, position: str, date_mode: str) -> str:
    quoted = ", ".join(f"'{_escape(gene)}'" for gene in genes)
    parts = [
        date_sql(cutoff, date_mode),
        f"array_has_any(gene_symbols, [{quoted}])",
    ]
    rsid = (rsid or "").strip().lower()
    position = str(position or "").strip()
    if rsid:
        parts.append(f"NOT contains(lower(title || ' ' || abstract), '{_escape(rsid)}')")
    if position:
        parts.append(f"NOT contains(lower(title || ' ' || abstract), '{_escape(position)}')")
    return " AND ".join(parts)


_TABLE = None


def _ensure_scalar_indexes(table) -> None:
    have = set()
    for index in table.list_indices():
        for column in list(getattr(index, "columns", []) or []):
            have.add(column)
    wanted = (
        ("gene_symbols", "LABEL_LIST"),
        ("year", "BTREE"),
        ("month", "BTREE"),
        ("day", "BTREE"),
        ("before_cutoff", "BITMAP"),
    )
    for column, kind in wanted:
        if column not in have:
            table.create_scalar_index(column, index_type=kind, replace=False)


def open_table():
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    import lancedb

    database = lancedb.connect(str(DB_DIR))
    table = database.open_table("papers")
    _ensure_scalar_indexes(table)
    for index in table.list_indices():
        kind = str(getattr(index, "index_type", "") or "").upper()
        if "IVF" in kind or "HNSW" in kind or "PQ" in kind:
            raise RuntimeError(f"expected no approximate vector index, found {kind}")
    _TABLE = table
    return table


def search(vector: np.ndarray, genes: list[str], cutoff: str, rsid: str, position: str, date_mode: str) -> list[dict]:
    """LanceDB vector search with gene and date prefilters.

    The lead variant is excluded in the same prefilter, before the top 25.
    Equal distances break by PubMed id, ascending.
    """
    if not genes:
        return []
    query = np.asarray(vector, dtype=np.float32)
    norm = float(np.linalg.norm(query))
    if norm:
        query = query / norm
    sql = filter_sql(genes, cutoff, rsid, position, date_mode)
    table = open_table()

    def run(limit: int) -> list[dict]:
        return (
            table.search(query, vector_column_name="vector")
            .distance_type("dot")
            .where(sql, prefilter=True)
            .select(["pmid", "year", "title", "abstract", "_distance"])
            .limit(limit)
            .to_list()
        )

    rows = run(TOP_N + 1)
    if len(rows) > TOP_N and float(rows[TOP_N - 1]["_distance"]) == float(rows[TOP_N]["_distance"]):
        rows = run(200)
    rows = sorted(rows, key=lambda row: (float(row["_distance"]), int(row["pmid"])))[:TOP_N]
    return [
        {
            "pmid": int(row["pmid"]),
            "year": int(row.get("year") or 0),
            "title": row.get("title") or "",
            "abstract": row.get("abstract") or "",
        }
        for row in rows
    ]


def docs_from_sqlite(pmids: list[int]) -> list[dict]:
    if not pmids:
        return []
    conn = sqlite3.connect(ABSTRACTS)
    marks = ",".join("?" for _ in pmids)
    found = {
        int(pmid): {"pmid": int(pmid), "year": year, "title": title or "", "abstract": abstract or ""}
        for pmid, year, title, abstract in conn.execute(
            f"SELECT pmid, year, title, abstract FROM abstracts WHERE pmid IN ({marks})",
            pmids,
        )
    }
    conn.close()
    return [found[pmid] for pmid in pmids if pmid in found]

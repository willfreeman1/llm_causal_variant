"""Compare database search with the lists the locked runs used.

A question passes when the 25 paper ids are the same set. A different order
is reported and does not fail the check. The exam did not save the paper ids
it sent, so those are recomputed with the old ranking and compared. The older
sealed run saved them.
"""

from __future__ import annotations

import csv
import json
import time

import numpy as np

from pipeline.paths import LOCAL, OLDER, QUESTIONS, RESULTS, embed_cache
from pipeline.questions import load_query_vectors, symbols
from pipeline.retrieve import search

TIE_ORDER_CHANGES: list[int] = []
TIE_AT_CUT: list[int] = []

SEALED = OLDER / "heldout_paper_answers.jsonl"
OUT = LOCAL / "ukb" / "final_exam" / "retrieval_check.json"


def rank_numpy(matrix: np.ndarray, index: dict[int, int], pmids: list[int], query: np.ndarray) -> list[int]:
    usable = [pmid for pmid in pmids if pmid in index]
    if not usable:
        return []
    scores = matrix[[index[pmid] for pmid in usable]] @ query
    order = np.argsort(-scores)
    by_pmid = np.lexsort((np.asarray(usable, dtype=np.int64), -scores))
    if list(order[:25]) != list(by_pmid[:25]):
        TIE_ORDER_CHANGES.append(1)
    if len(order) > 25 and float(scores[order[24]]) == float(scores[order[25]]):
        TIE_AT_CUT.append(1)
    return [usable[int(i)] for i in order[:25]]


def load_needed(needed: set[int]) -> tuple[np.ndarray, dict[int, int]]:
    pmids = []
    vectors = []
    for path in sorted(embed_cache().glob("chunk_*.npz")):
        data = np.load(path)
        arr = data["pmids"].astype(np.int64)
        keep = np.isin(arr, list(needed)) if False else np.fromiter((int(pmid) in needed for pmid in arr), dtype=bool, count=len(arr))
        if not keep.any():
            continue
        pmids.append(arr[keep])
        vectors.append(data["vectors"].astype(np.float32)[keep])
    matrix = np.vstack(vectors)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    matrix /= norms
    index = {}
    for i, pmid in enumerate(np.concatenate(pmids)):
        index[int(pmid)] = i
    return matrix, index


def split_diff(got: list[int], expected: list[int], ident: dict) -> tuple[dict | None, dict | None]:
    if set(got) != set(expected) or len(got) != len(expected):
        ident = dict(ident)
        ident["got"] = got
        ident["expected"] = expected
        return ident, None
    if got != expected:
        ident = dict(ident)
        ident["positions"] = [index for index, (left, right) in enumerate(zip(got, expected)) if left != right]
        return None, ident
    return None, None


def check_sealed() -> dict:
    saved_rows = [json.loads(line) for line in SEALED.read_text(encoding="utf-8").splitlines() if line.strip()]
    phenotypes = {}
    step1 = {}
    with (QUESTIONS / "uk_biobank_questions.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            step1[row["locus_id"]] = row
            phenotypes[row["locus_id"]] = row["phenotype"]
    queries = load_query_vectors(embed_cache() / "query_vectors.npz")
    set_mismatches = []
    order_differences = []
    times = []
    for number, row in enumerate(saved_rows, start=1):
        locus = str(row["locus"])
        meta = step1[locus]
        variant = meta["variant"].strip().strip('"')
        parts = variant.split(":")
        position = parts[1] if len(parts) >= 3 and parts[1].isdigit() else ""
        phenotype = phenotypes[locus]
        started = time.perf_counter()
        docs = search(
            queries[phenotype],
            symbols(row["symbol_gene_string"]),
            "2019-10-01",
            variant,
            position,
            "october_2019",
        )
        times.append(time.perf_counter() - started)
        got = [doc["pmid"] for doc in docs]
        expected = [int(pmid) for pmid in row["pmids"]]
        missing, reordered = split_diff(got, expected, {"id": row["row_number"], "locus": locus})
        if missing:
            set_mismatches.append(missing)
        if reordered:
            order_differences.append(reordered)
        if number % 50 == 0:
            print(f"SEALED {number} set_mismatches {len(set_mismatches)} order {len(order_differences)}", flush=True)
    return {
        "n": len(saved_rows),
        "set_mismatches": len(set_mismatches),
        "order_differences": len(order_differences),
        "set_examples": set_mismatches[:5],
        "order_examples": order_differences[:5],
        "median_seconds": float(np.median(times)) if times else None,
        "p95_seconds": float(np.percentile(times, 95)) if times else None,
    }


def check_exam() -> dict:
    pools = {}
    for line in (LOCAL / "ukb" / "final_exam" / "pools.jsonl").read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            pools[row["id"]] = row
    queries = load_query_vectors(LOCAL / "ukb" / "final_exam" / "query_vectors.npz")

    rows = list(csv.DictReader((QUESTIONS / "new_full.tsv").open(encoding="utf-8"), delimiter="\t"))
    extra = []
    full_keys = {(row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"]) for row in rows}
    for row in csv.DictReader((QUESTIONS / "new_main.tsv").open(encoding="utf-8"), delimiter="\t"):
        key = (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])
        if key not in full_keys:
            extra.append(row)
    rows = rows + extra
    needed = {int(pmid) for pool in pools.values() for pmid in pool["pmids"]}
    print(f"LOAD_MATRIX {len(needed)}", flush=True)
    matrix, index = load_needed(needed)
    set_mismatches = []
    order_differences = []
    times = []
    for number, row in enumerate(rows, start=1):
        query = queries[row["neighborhood_id"]]
        expected = rank_numpy(matrix, index, pools[row["neighborhood_id"]]["pmids"], query)
        started = time.perf_counter()
        docs = search(
            query,
            symbols(row["candidate_genes"]),
            row["earliest_public_date"],
            row.get("noncoding_rsid") or "",
            str(row.get("noncoding_pos") or ""),
            "study",
        )
        times.append(time.perf_counter() - started)
        got = [doc["pmid"] for doc in docs]
        missing, reordered = split_diff(got, expected, {"id": row["neighborhood_id"]})
        if missing:
            set_mismatches.append(missing)
        if reordered:
            order_differences.append(reordered)
        if number % 100 == 0:
            print(f"EXAM {number} set_mismatches {len(set_mismatches)} order {len(order_differences)}", flush=True)
    return {
        "n": len(rows),
        "set_mismatches": len(set_mismatches),
        "order_differences": len(order_differences),
        "set_examples": set_mismatches[:5],
        "order_examples": order_differences[:5],
        "median_seconds": float(np.median(times)) if times else None,
        "p95_seconds": float(np.percentile(times, 95)) if times else None,
    }


def main() -> None:
    print("SEALED", flush=True)
    sealed = check_sealed()
    print(json.dumps(sealed), flush=True)
    print("EXAM", flush=True)
    exam = check_exam()
    print(json.dumps({k: exam[k] for k in ("n", "set_mismatches", "order_differences", "order_examples", "median_seconds", "p95_seconds")}), flush=True)
    payload = {
        "sealed": sealed,
        "exam": exam,
        "numpy_ties_at_cut": len(TIE_AT_CUT),
        "numpy_tie_order_changes": len(TIE_ORDER_CHANGES),
    }
    OUT.write_text(json.dumps(payload), encoding="utf-8")
    if sealed["set_mismatches"] or exam["set_mismatches"]:
        raise SystemExit("retrieval set mismatch")


if __name__ == "__main__":
    main()

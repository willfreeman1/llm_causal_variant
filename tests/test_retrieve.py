"""Retrieval checks on real questions. Needs the local paper database."""

import json
import sqlite3

import pytest

from pipeline.paths import LOCAL, OLDER, QUESTIONS, embed_cache
from pipeline.questions import load_query_vectors, symbols
from pipeline.retrieve import DB_DIR, search

MANIFEST = DB_DIR / "manifest.json"
ABSTRACTS = LOCAL / "ukb" / "pubmed" / "weeks_abstracts.sqlite"


pytestmark = [
    pytest.mark.needs_db,
    pytest.mark.skipif(not MANIFEST.exists(), reason="paper database is not built yet"),
]


def _before(year, month, day, cutoff: str) -> bool:
    cy, cm, cd = (int(part) for part in cutoff.split("-"))
    year = int(year or 0)
    month = int(month or 0)
    day = int(day or 0)
    if year <= 0:
        return False
    if month <= 0:
        return year < cy
    if day <= 0:
        return (year, month) < (cy, cm)
    return (year, month, day) < (cy, cm, cd)


def test_real_questions_respect_the_date_and_the_variant():
    import csv

    rows = list(csv.DictReader((QUESTIONS / "new_full.tsv").open(encoding="utf-8"), delimiter="\t"))
    picked = [row for row in rows if row["noncoding_rsid"] and row["earliest_public_date"]][:3]
    assert len(picked) == 3
    vectors = load_query_vectors(LOCAL / "ukb" / "final_exam" / "query_vectors.npz")
    conn = sqlite3.connect(ABSTRACTS)
    for row in picked:
        docs = search(
            vectors[row["neighborhood_id"]],
            symbols(row["candidate_genes"]),
            row["earliest_public_date"],
            row["noncoding_rsid"],
            row["noncoding_pos"],
            "study",
        )
        assert docs, row["neighborhood_id"]
        pmids = [doc["pmid"] for doc in docs]
        marks = ",".join("?" for _ in pmids)
        for pmid, year, month, day, title, abstract in conn.execute(
            f"SELECT pmid, year, month, day, title, abstract FROM abstracts WHERE pmid IN ({marks})",
            pmids,
        ):
            assert _before(year, month, day, row["earliest_public_date"]), pmid
            text = f"{title or ''} {abstract or ''}".lower()
            assert row["noncoding_rsid"].lower() not in text
            assert str(row["noncoding_pos"]) not in text
    conn.close()


def test_first_saved_paper_list_matches():
    import csv

    saved = json.loads((OLDER / "heldout_paper_answers.jsonl").read_text(encoding="utf-8").splitlines()[0])
    with (QUESTIONS / "uk_biobank_questions.tsv").open(encoding="utf-8", newline="") as handle:
        step1 = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    meta = step1[str(saved["locus"])]
    queries = load_query_vectors(embed_cache() / "query_vectors.npz")
    variant = meta["variant"].strip().strip('"')
    parts = variant.split(":")
    position = parts[1] if len(parts) >= 3 and parts[1].isdigit() else ""
    docs = search(
        queries[meta["phenotype"]],
        symbols(saved["symbol_gene_string"]),
        "2019-10-01",
        variant,
        position,
        "october_2019",
    )
    assert [doc["pmid"] for doc in docs] == [int(pmid) for pmid in saved["pmids"]]

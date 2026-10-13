"""Load the saved questions and answers. No model calls."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from pipeline.paths import QUESTIONS, RESULTS

EXAM = RESULTS
QUESTION_DIR = QUESTIONS
WORKSHOP = {4, 5, 7, 22}


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def answer_key(row: dict) -> tuple:
    if "noncoding_chrom" in row:
        return (
            str(row["noncoding_chrom"]),
            str(row["noncoding_pos"]),
            str(row["noncoding_rsid"]),
            row["label_symbol"],
            row["trait_name"],
        )
    return (str(row["chrom"]), str(row["pos"]), str(row["rsid"]), row["label"], row["trait"])


def load_answers(path: Path) -> dict[tuple, dict]:
    found = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            found[answer_key(row)] = row
    return found


def symbols(gene_string: str) -> list[str]:
    return [part.strip().strip("{}") for part in gene_string.split(",") if part.strip(" {}")]


def load_query_vectors(path: Path) -> dict[str, np.ndarray]:
    data = np.load(path, allow_pickle=True)
    return {key: np.asarray(data[key], dtype=np.float32) for key in data.files}


def exam_rows() -> tuple[list[dict], list[dict], list[dict]]:
    full = read_tsv(QUESTION_DIR / "new_full.tsv")
    main = read_tsv(QUESTION_DIR / "new_main.tsv")
    return full, main, [row for row in full if row["in_exam_b"] == "1"]

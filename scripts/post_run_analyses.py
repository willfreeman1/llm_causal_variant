"""Analyses done after the locked exam. Saved answers only. No model calls.

Reads results_nopaper.jsonl, results_papers.jsonl, and exam_main.tsv.
Writes data/ukb/final_exam/post_run_analyses.json.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
EXAM = ROOT / "data" / "ukb" / "final_exam"
OUT = EXAM / "post_run_analyses.json"
PROTEIN = re.compile(r"^Circulating ([A-Za-z0-9]+) levels\b", re.I)
ALIAS = {"NTPROBNP": "NPPB"}


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


def load_jsonl(path: Path) -> dict[tuple, dict]:
    found = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            found[answer_key(row)] = row
    return found


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def measured_protein(trait: str) -> str:
    match = PROTEIN.match(trait or "")
    if not match:
        return ""
    token = match.group(1).upper()
    return ALIAS.get(token, token)


def mcnemar(only_a: int, only_b: int) -> float:
    total = only_a + only_b
    if total == 0:
        return 1.0
    return float(binomtest(only_a, total, 0.5, alternative="two-sided").pvalue)


def paired_slice(rows: list[dict], nopaper: dict[tuple, dict]) -> dict:
    only_gpt = 0
    only_near = 0
    gpt_right = 0
    near_right = 0
    for row in rows:
        answer = nopaper[answer_key(row)]
        gpt = bool(answer.get("correct"))
        near = row["nearest_gene"].upper() == row["label_symbol"].upper()
        gpt_right += int(gpt)
        near_right += int(near)
        only_gpt += int(gpt and not near)
        only_near += int(near and not gpt)
    return {
        "n": len(rows),
        "gpt_correct": gpt_right,
        "nearest_correct": near_right,
        "only_gpt_right": only_gpt,
        "only_nearest_right": only_near,
        "mcnemar_p": mcnemar(only_gpt, only_near),
    }


def main() -> None:
    nopaper = load_jsonl(EXAM / "results_nopaper.jsonl")
    papers = load_jsonl(EXAM / "results_papers.jsonl")
    main_rows = read_tsv(EXAM / "exam_main.tsv")
    full_rows = read_tsv(EXAM / "exam_a.tsv")
    main_keys = {answer_key(row) for row in main_rows}

    other = [row for row in main_rows if row["trait_group"] not in ("blood_protein", "blood_lipid")]
    protein = [row for row in main_rows if row["trait_group"] == "blood_protein"]
    fresh = [row for row in main_rows if row["previously_reported"] == "0"]
    slices = {
        "all_main": paired_slice(main_rows, nopaper),
        "diseases_and_other": paired_slice(other, nopaper),
        "blood_protein": paired_slice(protein, nopaper),
        "not_previously_reported": paired_slice(fresh, nopaper),
    }

    def hallucination_report(answers: dict[tuple, dict]) -> dict:
        rows = [row for row in answers.values() if row.get("hallucination")]
        named_protein = 0
        blood = 0
        examples = []
        for row in rows:
            token = measured_protein(row.get("trait") or "")
            is_blood = bool(token)
            blood += int(is_blood)
            raw = (row.get("raw_predicted") or "").upper()
            if token and raw == token:
                named_protein += 1
            if len(examples) < 5:
                examples.append({
                    "id": row["id"],
                    "trait": row.get("trait"),
                    "raw_predicted": row.get("raw_predicted"),
                    "measured_protein": token,
                })
        return {
            "n": len(rows),
            "blood_protein_questions": blood,
            "named_the_measured_protein": named_protein,
            "on_main_set": sum(1 for key, row in answers.items() if row.get("hallucination") and key in main_keys),
            "examples": examples,
        }

    full_keys = {answer_key(row) for row in full_rows}
    breaks = []
    for row in full_rows:
        old = nopaper[answer_key(row)]
        new = papers[answer_key(row)]
        if old.get("correct") and not new.get("correct"):
            breaks.append(row["label_symbol"])
    break_counts = Counter(breaks)

    main_fixed = 0
    main_broken = 0
    for row in main_rows:
        old = bool(nopaper[answer_key(row)].get("correct"))
        new = bool(papers[answer_key(row)].get("correct"))
        main_fixed += int((not old) and new)
        main_broken += int(old and (not new))

    payload = {
        "note": "Computed from saved answers after the run. No new model calls. Not locked in advance.",
        "gpt_vs_nearest_no_papers": slices,
        "out_of_list": {
            "no_papers": hallucination_report(nopaper),
            "with_papers": hallucination_report(papers),
            "guardrail": "A gene not in the candidate list was stored as a blank answer and counted wrong.",
        },
        "repeat_collapse": {
            "full_set_n": len(full_rows),
            "full_set_fixed": sum(
                1
                for row in full_rows
                if (not nopaper[answer_key(row)].get("correct"))
                and papers[answer_key(row)].get("correct")
            ),
            "full_set_broken": len(breaks),
            "full_set_broken_mcnemar_p": mcnemar(
                sum(
                    1
                    for row in full_rows
                    if (not nopaper[answer_key(row)].get("correct"))
                    and papers[answer_key(row)].get("correct")
                ),
                len(breaks),
            ),
            "asgr1_breaks": break_counts.get("ASGR1", 0),
            "break_top": break_counts.most_common(8),
            "main_fixed": main_fixed,
            "main_broken": main_broken,
            "main_net": main_fixed - main_broken,
            "main_mcnemar_p": mcnemar(main_fixed, main_broken),
            "answers_not_in_full_set": sorted(row["id"] for key, row in nopaper.items() if key not in full_keys),
        },
        "pooled_papers_caveat": {
            "uk_fixed": 64,
            "uk_broken": 26,
            "main_fixed": main_fixed,
            "main_broken": main_broken,
            "pooled_fixed": 64 + main_fixed,
            "pooled_broken": 26 + main_broken,
            "pooled_mcnemar_p": mcnemar(64 + main_fixed, 26 + main_broken),
            "main_only_mcnemar_p": mcnemar(main_fixed, main_broken),
            "note": "The UK rows added papers only where the two methods disagreed, and those studies were public before the model's training cutoff. The new main set added papers on every question.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

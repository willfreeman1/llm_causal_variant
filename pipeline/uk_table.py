"""Join the older 1,348 questions to the saved model and statistical-model answers."""

from __future__ import annotations

import csv
import json
from collections import defaultdict

from pipeline.paths import OLDER, QUESTIONS, RESULTS


def _norm(gene: str) -> str:
    return (gene or "").strip().upper()


def load_rows() -> list[dict]:
    weeks = {}
    order = []
    with (QUESTIONS / "uk_biobank_questions.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            weeks[row["locus_id"]] = row
            order.append(row["locus_id"])

    by_locus: dict[str, list[dict]] = defaultdict(list)
    with (RESULTS / "uk_caldera_genes.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_locus[row["locus"]].append(row)

    llm_by_row = {}
    with (OLDER / "uk_baseline_answers.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("dataset") == "weeks_et_al":
                llm_by_row[str(row["row_number"])] = row

    nearest = {}
    with (QUESTIONS / "uk_biobank_nearest_gene.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            nearest[row["locus_id"]] = _norm(row["causal_gene"])

    rows = []
    for locus in order:
        label = weeks[locus]
        label_symbol = _norm(label["symbol"])
        label_ensg = label["gene"]
        chrom = int(str(label["chromosome"]).replace("chr", ""))
        llm = llm_by_row[str(label["row_number"])]
        gpt = _norm(llm["predicted"])
        genes = by_locus.get(locus, [])
        if genes:
            genes = sorted(genes, key=lambda item: -float(item["caldera"]))
            top = genes[0]
            cal = _norm(top["gene"])
            cal_hit = cal == label_symbol or top["ensgid"] == label_ensg
            cal_p = float(top["caldera"])
        else:
            cal = ""
            cal_hit = False
            cal_p = 0.0
        near = nearest[locus]
        rows.append({
            "locus": locus,
            "chrom": chrom,
            "gpt": gpt,
            "gpt_hit": gpt == label_symbol,
            "gpt_conf": float(llm["confidence"]),
            "cal": cal,
            "cal_hit": cal_hit,
            "cal_p": cal_p,
            "agree": bool(cal) and gpt == cal,
            "near": near,
            "near_hit": near == label_symbol,
        })
    return rows


def hits(rows: list[dict], choose) -> int:
    total = 0
    for row in rows:
        total += row[choose(row) + "_hit"]
    return total


def choose_fixed(partner: str, threshold: float):
    def choose(row: dict) -> str:
        agreed = row["agree"] if partner == "cal" else row["gpt"] == row["near"]
        if agreed:
            return "gpt"
        if row["gpt_conf"] >= threshold:
            return "gpt"
        return partner

    return choose

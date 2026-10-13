"""Score CALDERA's top gene on the 1,348 neighborhoods against the label,
the language-model pick, and the nearest-gene pick.
"""

import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEEKS = ROOT / "data" / "ukb" / "zenodo_directory_may_2025" / "data" / "benchmark_datasets" / "weeks_et_al_step1.tsv"
CALDERA = ROOT / "data" / "ukb" / "caldera_results" / "weeks_caldera_genes.tsv"
LLM = ROOT / "spikes" / "shringarpure_match" / "results.jsonl"
NEAREST = ROOT / "data" / "ukb" / "predictions" / "weeks_et_al.nearest_gene.csv"
OUT = ROOT / "data" / "ukb" / "caldera_results" / "summary.json"


def norm(gene):
    return (gene or "").strip().upper()


def main():
    weeks = {}
    order = []
    with open(WEEKS, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            weeks[row["locus_id"]] = row
            order.append(row["locus_id"])

    by_locus = defaultdict(list)
    with open(CALDERA, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_locus[row["locus"]].append(row)

    llm_by_row = {}
    with open(LLM, encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["dataset"] == "weeks_et_al":
                llm_by_row[str(row["row_number"])] = row

    nearest = {}
    with open(NEAREST, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            nearest[row["locus_id"]] = row["causal_gene"]

    missing = []
    ties = []
    correct = only_caldera = only_llm = both_right = both_wrong = same_llm = 0
    only_caldera_vs_near = only_near = both_right_near = both_wrong_near = same_near = 0
    near_correct = llm_correct = 0
    all_three_same = 0
    llm_caldera_differ = 0
    any_of_three = 0
    imputed_top = 0

    for locus in order:
        label = weeks[locus]
        label_symbol = norm(label["symbol"])
        label_ensg = label["gene"]
        rows = by_locus.get(locus, [])
        if not rows:
            missing.append(locus)
            pred = ""
            pred_ensg = ""
        else:
            top_score = max(float(row["caldera"]) for row in rows)
            tops = [row for row in rows if float(row["caldera"]) == top_score]
            if len(tops) > 1:
                ties.append(locus)
            chosen = tops[0]
            pred = norm(chosen["gene"])
            pred_ensg = chosen["ensgid"]
            if chosen["impute_pops"] == "TRUE":
                imputed_top += 1
        hit = bool(pred) and (pred == label_symbol or pred_ensg == label_ensg)
        if hit:
            correct += 1

        llm = llm_by_row[str(label["row_number"])]
        llm_hit = bool(llm["correct"])
        llm_gene = norm(llm["predicted"])
        if llm_hit:
            llm_correct += 1
        near_gene = norm(nearest[locus])
        near_hit = near_gene == label_symbol
        if near_hit:
            near_correct += 1

        if pred and pred == llm_gene:
            same_llm += 1
        else:
            llm_caldera_differ += 1
        if hit and llm_hit:
            both_right += 1
        elif hit:
            only_caldera += 1
        elif llm_hit:
            only_llm += 1
        else:
            both_wrong += 1

        if pred and pred == near_gene:
            same_near += 1
        if hit and near_hit:
            both_right_near += 1
        elif hit:
            only_caldera_vs_near += 1
        elif near_hit:
            only_near += 1
        else:
            both_wrong_near += 1

        if pred and pred == llm_gene == near_gene:
            all_three_same += 1
        if hit or llm_hit or near_hit:
            any_of_three += 1

    n = len(order)
    summary = {
        "neighborhoods": n,
        "caldera_correct": correct,
        "caldera_no_gene": len(missing),
        "caldera_ties": len(ties),
        "top_gene_missing_pops_score": imputed_top,
        "gpt54_correct": llm_correct,
        "nearest_correct": near_correct,
        "same_gene_as_gpt54": same_llm,
        "gpt54_and_caldera_differ": llm_caldera_differ,
        "both_right_vs_gpt54": both_right,
        "only_caldera_vs_gpt54": only_caldera,
        "only_gpt54": only_llm,
        "both_wrong_vs_gpt54": both_wrong,
        "same_gene_as_nearest": same_near,
        "both_right_vs_nearest": both_right_near,
        "only_caldera_vs_nearest": only_caldera_vs_near,
        "only_nearest": only_near,
        "both_wrong_vs_nearest": both_wrong_near,
        "all_three_same_gene": all_three_same,
        "at_least_one_right": any_of_three,
        "missing_loci": missing,
    }
    OUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

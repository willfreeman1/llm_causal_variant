"""Score combining rules on the 1,348 UK neighborhoods.

No new model calls. GPT-5.4 answers, CALDERA probabilities, and nearest-gene
picks are already saved. Rules that learn a cutoff are scored by holding out
one chromosome at a time, so the cutoff is fit only on the other chromosomes.
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.tree import DecisionTreeClassifier

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.paths import OLDER, QUESTIONS, RESULTS

WEEKS = QUESTIONS / "uk_biobank_questions.tsv"
CALDERA = RESULTS / "uk_caldera_genes.tsv"
LLM = OLDER / "uk_baseline_answers.jsonl"
NEAREST = QUESTIONS / "uk_biobank_nearest_gene.csv"
OUT = RESULTS / "uk_chooser_summary.json"

THRESHOLDS = [round(x, 2) for x in np.arange(0.50, 1.00, 0.05)]
# CALDERA's probability is spread across the scale, so the grid starts lower than GPT's.
CALDERA_PROB_THRESHOLDS = [round(float(x), 2) for x in np.arange(0.20, 0.95, 0.05)]


def norm(gene):
    return (gene or "").strip().upper()


def load_rows():
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
            nearest[row["locus_id"]] = norm(row["causal_gene"])

    rows = []
    for locus in order:
        label = weeks[locus]
        label_symbol = norm(label["symbol"])
        label_ensg = label["gene"]
        chrom = int(str(label["chromosome"]).replace("chr", ""))
        llm = llm_by_row[str(label["row_number"])]
        gpt = norm(llm["predicted"])
        gpt_hit = gpt == label_symbol
        genes = by_locus.get(locus, [])
        if genes:
            genes = sorted(genes, key=lambda r: -float(r["caldera"]))
            top = genes[0]
            second_p = float(genes[1]["caldera"]) if len(genes) > 1 else 0.0
            cal = norm(top["gene"])
            cal_hit = cal == label_symbol or top["ensgid"] == label_ensg
            best_pops = max(float(g["pops"]) for g in genes)
            cal_p = float(top["caldera"])
            gap = cal_p - second_p
            best_pops_flag = float(top["pops"]) == best_pops
            coding_flag = float(top["coding"]) > 0
            n_genes = int(top["n_genes"])
        else:
            cal = ""
            cal_hit = False
            cal_p = 0.0
            gap = 0.0
            best_pops_flag = False
            coding_flag = False
            n_genes = 0
        near = nearest[locus]
        rows.append({
            "locus": locus,
            "chrom": chrom,
            "gpt": gpt,
            "gpt_hit": gpt_hit,
            "gpt_conf": float(llm["confidence"]),
            "cal": cal,
            "cal_hit": cal_hit,
            "cal_p": cal_p,
            "gap": gap,
            "best_pops": best_pops_flag,
            "coding": coding_flag,
            "n_genes": n_genes,
            "agree": bool(cal) and gpt == cal,
            "near": near,
            "near_hit": near == label_symbol,
        })
    return rows


def hits(rows, choose):
    """choose(row) returns 'gpt', 'cal', or 'near'."""
    n = 0
    for row in rows:
        src = choose(row)
        n += row[src + "_hit"]
    return n


def choose_trait_match(row):
    """Trust CALDERA only when its gene is also the best trait-match gene."""
    if row["cal"] and row["best_pops"]:
        return "cal"
    return "gpt"


def choose_trait_match_then_sureness(threshold):
    """If CALDERA missed the best trait-match gene, take GPT. Otherwise use GPT's sureness line."""

    def choose(row):
        if row["agree"] or not row["best_pops"]:
            return "gpt"
        if row["gpt_conf"] >= threshold:
            return "gpt"
        return "cal"

    return choose


def choose_caldera_probability(threshold):
    """When the names differ, take CALDERA if its probability is at or above the line."""

    def choose(row):
        if row["agree"]:
            return "gpt"
        if row["cal"] and row["cal_p"] >= threshold:
            return "cal"
        return "gpt"

    return choose


def loco_custom(rows, thresholds, choose_from_threshold):
    chroms = sorted({row["chrom"] for row in rows})
    chosen_hits = 0
    fold_thresholds = []
    for chrom in chroms:
        train = [row for row in rows if row["chrom"] != chrom]
        test = [row for row in rows if row["chrom"] == chrom]
        best_t, best_score = None, -1
        for threshold in thresholds:
            score = hits(train, choose_from_threshold(threshold))
            if score > best_score:
                best_score = score
                best_t = threshold
        fold_thresholds.append({"chrom": chrom, "threshold": best_t, "n": len(test)})
        chosen_hits += hits(test, choose_from_threshold(best_t))
    return chosen_hits, fold_thresholds


def disagreement_bins(rows):
    """On neighborhoods where the names differ, how often each method is right."""
    bins = []
    disagree = [row for row in rows if not row["agree"]]
    for lo, hi in [(0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]:
        subset = [row for row in disagree if lo <= row["cal_p"] < hi]
        bins.append({
            "caldera_probability": f"{lo:.1f}-{min(hi, 1.0):.1f}",
            "n": len(subset),
            "caldera_right": sum(row["cal_hit"] for row in subset),
            "gpt_right": sum(row["gpt_hit"] for row in subset),
        })
    return bins


def choose_fixed(partner, threshold):
    def choose(row):
        if row["agree"] if partner == "cal" else row["gpt"] == row["near"]:
            return "gpt"
        if row["gpt_conf"] >= threshold:
            return "gpt"
        return partner
    return choose


def loco_threshold(rows, partner):
    """Pick a confidence cutoff on the other chromosomes. Apply it to the held-out one."""
    chroms = sorted({row["chrom"] for row in rows})
    chosen_hits = 0
    fold_thresholds = []
    for chrom in chroms:
        train = [row for row in rows if row["chrom"] != chrom]
        test = [row for row in rows if row["chrom"] == chrom]
        best_t, best_score = None, -1
        for threshold in THRESHOLDS:
            score = hits(train, choose_fixed(partner, threshold))
            if score > best_score:
                best_score = score
                best_t = threshold
        fold_thresholds.append({"chrom": chrom, "threshold": best_t, "n": len(test)})
        chosen_hits += hits(test, choose_fixed(partner, best_t))
    return chosen_hits, fold_thresholds


def loco_tree(rows, partner, feature_fn, depth):
    chroms = sorted({row["chrom"] for row in rows})
    total = 0
    for chrom in chroms:
        train = [row for row in rows if row["chrom"] != chrom]
        test = [row for row in rows if row["chrom"] == chrom]
        # Fit only where the two methods name different genes and exactly one is right.
        # Both-wrong rows have no preferred source.
        fit_rows = [
            row for row in train
            if row["gpt"] != row[partner] and (row["gpt_hit"] != row[partner + "_hit"])
        ]
        y = np.array([1 if row["gpt_hit"] else 0 for row in fit_rows])
        if len(fit_rows) < 30 or len(set(y.tolist())) < 2:
            choose = choose_fixed(partner, 0.9)
            total += hits(test, choose)
            continue
        model = DecisionTreeClassifier(max_depth=depth, min_samples_leaf=25, random_state=0)
        model.fit(np.array([feature_fn(row) for row in fit_rows]), y)

        def choose(row, model=model):
            same = row["gpt"] == row[partner] and bool(row[partner])
            if same:
                return "gpt"
            use_gpt = int(model.predict(np.array([feature_fn(row)]))[0]) == 1
            return "gpt" if use_gpt else partner

        total += hits(test, choose)
    return total


def main():
    rows = load_rows()
    n = len(rows)
    summary = {
        "neighborhoods": n,
        "gpt_correct": sum(row["gpt_hit"] for row in rows),
        "caldera_correct": sum(row["cal_hit"] for row in rows),
        "nearest_correct": sum(row["near_hit"] for row in rows),
        "gpt_caldera_agree": sum(row["agree"] for row in rows),
        "gpt_nearest_agree": sum(row["gpt"] == row["near"] for row in rows),
    }

    summary["fixed_0.9_gpt_caldera"] = hits(rows, choose_fixed("cal", 0.9))
    summary["fixed_0.9_gpt_nearest"] = hits(rows, choose_fixed("near", 0.9))

    cal_hits, cal_folds = loco_threshold(rows, "cal")
    near_hits, near_folds = loco_threshold(rows, "near")
    summary["loco_threshold_gpt_caldera"] = cal_hits
    summary["loco_threshold_gpt_nearest"] = near_hits
    summary["loco_thresholds_caldera"] = cal_folds
    summary["loco_thresholds_nearest"] = near_folds

    summary["loco_tree_agree_confidence_caldera"] = loco_tree(
        rows, "cal", lambda row: [row["gpt_conf"]], depth=1
    )
    summary["loco_tree_rich_caldera"] = loco_tree(
        rows,
        "cal",
        lambda row: [row["gpt_conf"], row["cal_p"], row["gap"], int(row["best_pops"]), int(row["coding"])],
        depth=2,
    )

    # Descriptive tree fit on every chromosome, for reading the rule. Not the score.
    fit_rows = [
        row for row in rows
        if row["gpt"] != row["cal"] and (row["gpt_hit"] != row["cal_hit"])
    ]
    rich = DecisionTreeClassifier(max_depth=2, min_samples_leaf=25, random_state=0)
    rich.fit(
        np.array([[row["gpt_conf"], row["cal_p"], row["gap"], int(row["best_pops"]), int(row["coding"])] for row in fit_rows]),
        np.array([1 if row["gpt_hit"] else 0 for row in fit_rows]),
    )
    summary["descriptive_rich_tree_text"] = rich_tree_text(rich)

    cal_p_hits, cal_p_folds = loco_custom(rows, CALDERA_PROB_THRESHOLDS, choose_caldera_probability)
    summary["loco_caldera_probability_line"] = cal_p_hits
    summary["loco_caldera_probability_thresholds"] = cal_p_folds
    summary["disagreement_bins_by_caldera_probability"] = disagreement_bins(rows)

    summary["trait_match_rule"] = hits(rows, choose_trait_match)
    trait_hits, trait_folds = loco_custom(rows, THRESHOLDS, choose_trait_match_then_sureness)
    summary["loco_trait_match_then_sureness"] = trait_hits
    summary["loco_trait_match_then_sureness_thresholds"] = trait_folds
    disagree = [row for row in rows if not row["agree"]]
    summary["disagree_caldera_is_best_trait_match"] = {
        "n": sum(row["best_pops"] for row in disagree),
        "caldera_right": sum(row["best_pops"] and row["cal_hit"] for row in disagree),
        "gpt_right": sum(row["best_pops"] and row["gpt_hit"] for row in disagree),
    }
    summary["disagree_caldera_is_not_best_trait_match"] = {
        "n": sum(not row["best_pops"] for row in disagree),
        "caldera_right": sum((not row["best_pops"]) and row["cal_hit"] for row in disagree),
        "gpt_right": sum((not row["best_pops"]) and row["gpt_hit"] for row in disagree),
    }

    OUT.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    hide = {
        "loco_thresholds_caldera",
        "loco_thresholds_nearest",
        "loco_caldera_probability_thresholds",
        "loco_trait_match_then_sureness_thresholds",
        "descriptive_rich_tree_text",
    }
    show = {key: summary[key] for key in summary if key not in hide}
    print(json.dumps(show, indent=2))
    print("thresholds caldera probability", sorted({fold["threshold"] for fold in cal_p_folds}))
    print("thresholds caldera", sorted({fold["threshold"] for fold in cal_folds}))
    print("thresholds nearest", sorted({fold["threshold"] for fold in near_folds}))
    print(summary["descriptive_rich_tree_text"])


def rich_tree_text(model):
    names = ["gpt_conf", "caldera_p", "gap", "best_pops", "coding"]
    lines = []

    def walk(node, depth):
        if model.tree_.feature[node] == -2:
            lines.append("  " * depth + ("take GPT" if model.tree_.value[node][0][1] >= model.tree_.value[node][0][0] else "take CALDERA"))
            return
        name = names[model.tree_.feature[node]]
        cut = model.tree_.threshold[node]
        lines.append("  " * depth + f"if {name} <= {cut:.3f}")
        walk(model.tree_.children_left[node], depth + 1)
        lines.append("  " * depth + f"if {name} > {cut:.3f}")
        walk(model.tree_.children_right[node], depth + 1)

    walk(0, 0)
    return "\n".join(lines)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        LLM = Path(sys.argv[1])
        OUT = Path(sys.argv[2])
    main()

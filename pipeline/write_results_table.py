"""Write results_table.md from the saved scores. No model calls."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from pipeline.decide import LINE, final_decision
from pipeline.paths import EARLY, OLDER, QUESTIONS, RESULTS, ROOT
from pipeline.questions import WORKSHOP, answer_key, load_answers
from pipeline.reproduce import load_sealed, uk_state
from pipeline.uk_table import choose_fixed, hits, load_rows

OUT = ROOT / "results_table.md"
SEED = 20261009
DRAWS = 2000


def interval(flags: list[bool]) -> str:
    arr = np.asarray(flags, dtype=np.float64)
    if len(arr) == 0:
        return ""
    rng = np.random.default_rng(SEED)
    draws = rng.integers(0, len(arr), size=(DRAWS, len(arr)))
    rates = arr[draws].mean(axis=1) * 100
    low, high = np.percentile(rates, [2.5, 97.5])
    return f"{low:.1f} to {high:.1f}"


def money(amount: float) -> str:
    if amount == 0:
        return "$0"
    if amount < 0.01:
        return f"${amount:.4f}"
    return f"${amount:.3f}"


def table_row(name: str, hits: int, n: int, flags: list[bool] | None, cost: float, note: str = "") -> str:
    percent = 100 * hits / n if n else 0
    band = interval(flags) if flags is not None else ""
    label = name if not note else f"{name} {note}"
    per = cost / n if n else 0
    return f"| {label} | {hits} / {n} | {percent:.1f} | {band} | {money(per)} |"


def sum_costs(path: Path, predicate) -> float:
    total = 0.0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        item = json.loads(line)
        if predicate(item):
            total += float(item.get("cost_usd") or 0)
    return total


def main() -> None:
    rows = load_rows()
    summary = json.loads((RESULTS / "uk_chooser_summary.json").read_text(encoding="utf-8"))
    lines = {int(item["chrom"]): float(item["threshold"]) for item in summary["loco_thresholds_caldera"]}
    near_lines = {int(item["chrom"]): float(item["threshold"]) for item in summary["loco_thresholds_nearest"]}
    sealed = load_sealed()
    weeks_cost = {}
    for line in (OLDER / "uk_baseline_answers.jsonl").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        item = json.loads(line)
        if item["dataset"] == "weeks_et_al":
            weeks_cost[str(item["row_number"])] = float(item["cost_usd"])
    # chooser rows have no row_number. Cost is joined through the weeks file.
    weeks_rows = {}
    with (QUESTIONS / "uk_biobank_questions.tsv").open(encoding="utf-8", newline="") as handle:
        for item in csv.DictReader(handle, delimiter="\t"):
            weeks_rows[item["locus_id"]] = item

    def flags_for(subset, line_of, papers: bool) -> list[bool]:
        flags = []
        for row in subset:
            paper = sealed.get(str(row["locus"])) if papers and not row["agree"] else None
            _pick, source, _review = final_decision(uk_state(row, line_of(row), paper))
            if source in ("agreement", "gpt"):
                flags.append(bool(paper.get("correct")) if paper else bool(row["gpt_hit"]))
            else:
                flags.append(bool(row["cal_hit"]))
        return flags

    def cost_of(subset, papers: bool) -> float:
        total = 0.0
        for row in subset:
            total += weeks_cost.get(str(weeks_rows[str(row["locus"])]["row_number"]), 0.0)
            if papers and not row["agree"]:
                paper = sealed.get(str(row["locus"]))
                if paper:
                    total += float(paper.get("cost_usd") or 0)
        return total

    held = [row for row in rows if row["chrom"] not in WORKSHOP]
    loco_flags = flags_for(rows, lambda row: lines[row["chrom"]], False)
    fixed_flags = flags_for(rows, lambda row: LINE, False)
    near_flags = [bool(row["near_hit"]) for row in rows]
    cal_flags = [bool(row["cal_hit"]) for row in rows]
    gpt_flags = [bool(row["gpt_hit"]) for row in rows]
    near_hits = hits(rows, lambda item: choose_fixed("near", near_lines[item["chrom"]])(item))
    near_partner = []
    for item in rows:
        source = choose_fixed("near", near_lines[item["chrom"]])(item)
        near_partner.append(bool(item[source + "_hit"]))

    uk_cost = cost_of(rows, False)
    held_cost = cost_of(held, False)
    held_paper_cost = cost_of(held, True)
    held_no = flags_for(held, lambda row: LINE, False)
    held_yes = flags_for(held, lambda row: LINE, True)

    exam = json.loads((RESULTS / "new_summary.json").read_text(encoding="utf-8"))
    post = json.loads((RESULTS / "post_run_analyses.json").read_text(encoding="utf-8"))
    replay = json.loads((RESULTS / "replay_numbers.json").read_text(encoding="utf-8"))
    nopaper = load_answers(RESULTS / "new_nopaper.jsonl")
    papers = load_answers(RESULTS / "new_papers.jsonl")
    main = list(csv.DictReader((QUESTIONS / "new_main.tsv").open(encoding="utf-8"), delimiter="\t"))
    full = list(csv.DictReader((QUESTIONS / "new_full.tsv").open(encoding="utf-8"), delimiter="\t"))
    exam_b = [row for row in full if row["in_exam_b"] == "1"]

    def exam_cost(subset, path_answers: dict) -> float:
        return sum(float(path_answers[answer_key(row)].get("cost_usd") or 0) for row in subset)

    def band_from_summary(block: dict) -> str:
        low, high = block["question_interval"]
        return f"{low:.1f} to {high:.1f}"

    lines_out = []
    add = lines_out.append
    add("# Results")
    add("")
    add("Each interval redraws whole questions 2,000 times (seed 20261009) and reports the middle 95 percent of those scores. Cost is the model spend for that row, divided by the number of questions.")
    add("")
    add("## Practice, United Kingdom Biobank, 1,348 questions")
    add("")
    add("Practice; possible training exposure. Studies in this table were public before the language model's knowledge cutoff.")
    add("")
    add("| Method | Correct | % | 95% interval | Cost per question |")
    add("|---|---:|---:|---:|---:|")
    add(table_row("Nearest gene", sum(near_flags), 1348, near_flags, 0))
    add(table_row("Statistical model", sum(cal_flags), 1348, cal_flags, 0))
    add(table_row("Language model, no papers", sum(gpt_flags), 1348, gpt_flags, uk_cost))
    add(table_row("Combining rule, per-chromosome line", sum(loco_flags), 1348, loco_flags, uk_cost))
    add(table_row("Combining rule, fixed 0.85 line", sum(fixed_flags), 1348, fixed_flags, uk_cost))
    add(table_row("Rule with the nearest gene as the partner", near_hits, 1348, near_partner, uk_cost))
    add("")
    add("The per-chromosome line is the one fit on the other chromosomes. It is right on 869. A single 0.85 line on all 1,348 is right on "
        f"{sum(fixed_flags)}. The nearest-gene partner, with its own per-chromosome line, is right on {near_hits}.")
    add("")
    add("## 18 chromosomes held out of the paper-setup choice")
    add("")
    add("| Method | Correct | % | 95% interval | Cost per question |")
    add("|---|---:|---:|---:|---:|")
    add(table_row("Combining rule, no papers", sum(held_no), len(held), held_no, held_cost))
    add(table_row("Combining rule, with papers", sum(held_yes), len(held), held_yes, held_paper_cost))
    add("")
    add("Papers were requested only where the two methods disagreed. Chromosomes 4, 5, 7, and 22 are not in this table.")
    add("")
    add("## New studies, main list, 429 questions")
    add("")
    add("This is the primary new-study result. Studies were first published after 31 August 2025.")
    add("")
    add("| Method | Correct | % | 95% interval | Cost per question |")
    add("|---|---:|---:|---:|---:|")
    main_block = exam["slices"]["main"]
    add(f"| Nearest gene | {main_block['nearest']['correct']} / 429 | {main_block['nearest']['percent']:.1f} | {band_from_summary(main_block['nearest'])} | $0 |")
    add(f"| Language model, no papers | {main_block['no_papers']['correct']} / 429 | {main_block['no_papers']['percent']:.1f} | {band_from_summary(main_block['no_papers'])} | {money(exam_cost(main, nopaper) / 429)} |")
    add(f"| Language model, with papers | {main_block['with_papers']['correct']} / 429 | {main_block['with_papers']['percent']:.1f} | {band_from_summary(main_block['with_papers'])} | {money(exam_cost(main, papers) / 429)} |")
    add("")
    add("By kind of trait, language model with no papers against the nearest gene. Counted after the run, not locked beforehand.")
    add("")
    add("| Slice | Questions | Only the language model right | Only the nearest gene right | Paired p |")
    add("|---|---:|---:|---:|---:|")
    labels = {
        "all_main": "All main",
        "diseases_and_other": "Diseases and other traits",
        "blood_protein": "Blood protein levels",
        "not_previously_reported": "Not previously reported",
    }
    for key, label in labels.items():
        item = post["gpt_vs_nearest_no_papers"][key]
        add(f"| {label} | {item['n']} | {item['only_gpt_right']} | {item['only_nearest_right']} | {item['mcnemar_p']:.4f} |")
    add("")
    add("## New studies where both methods have an input, 38 questions")
    add("")
    add("Check that the setup runs; too small to conclude.")
    add("")
    add("| Method | Correct | % | 95% interval | Cost per question |")
    add("|---|---:|---:|---:|---:|")
    b = exam["slices"]
    add(f"| Nearest gene | {b['exam_b_rule']['nearest']['correct']} / 38 | {b['exam_b_rule']['nearest']['percent']:.1f} | {band_from_summary(b['exam_b_rule']['nearest'])} | $0 |")
    add(f"| Language model, no papers | {b['exam_b_gpt']['no_papers']['correct']} / 38 | {b['exam_b_gpt']['no_papers']['percent']:.1f} | {band_from_summary(b['exam_b_gpt']['no_papers'])} | {money(exam_cost(exam_b, nopaper) / 38)} |")
    add(f"| Language model, with papers | {b['exam_b_gpt']['with_papers']['correct']} / 38 | {b['exam_b_gpt']['with_papers']['percent']:.1f} | {band_from_summary(b['exam_b_gpt']['with_papers'])} | {money(exam_cost(exam_b, papers) / 38)} |")
    add(f"| Statistical model alone | {b['exam_b_caldera']['correct']} / 38 | {b['exam_b_caldera']['percent']:.1f} | {band_from_summary(b['exam_b_caldera'])} | $0 |")
    add(f"| Combining rule, no papers | {b['exam_b_rule']['no_papers']['correct']} / 38 | {b['exam_b_rule']['no_papers']['percent']:.1f} | {band_from_summary(b['exam_b_rule']['no_papers'])} | {money(exam_cost(exam_b, nopaper) / 38)} |")
    add(f"| Combining rule, with papers | {b['exam_b_rule']['with_papers']['correct']} / 38 | {b['exam_b_rule']['with_papers']['percent']:.1f} | {band_from_summary(b['exam_b_rule']['with_papers'])} | {money((exam_cost(exam_b, nopaper) + exam_cost(exam_b, papers)) / 38)} |")
    add("")
    add("## Earlier checks")
    add("")
    bias = json.loads((EARLY / "list_order_summary.json").read_text(encoding="utf-8"))
    memo = json.loads((EARLY / "memorization_results.json").read_text(encoding="utf-8"))["summary"]
    add("List order and fame, on 80 neighborhoods asked four times each. The top gene was the same on every order for "
        f"{bias['overall']['loci_same_pick_every_shuffle']} of 80. Putting the true gene in the first third of the list was right on "
        f"{bias['overall']['hit_rate_causal_gene_in_first_third'] * 100:.1f} percent of those trials, and in the last third on "
        f"{bias['overall']['hit_rate_causal_gene_in_last_third'] * 100:.1f} percent. When the answer was wrong, the chosen gene was more studied than the truth "
        f"{bias['overall']['wrong_pick_more_famous_than_truth'] * 100:.0f} percent of the time, against "
        f"{bias['overall']['chance_other_gene_more_famous_than_truth'] * 100:.0f} percent for a random other gene in the same list. "
        f"Spend {money(bias['spent_usd'])}.")
    add("")
    add("Memorization probe. Exact copy of a real Open Targets row: "
        f"{memo['ot_line']['joint_exact']} of {memo['ot_line']['n']}. "
        f"The same question on the iris flower file, which models often have memorized: {memo['iris_line_exact']['exact']} of {memo['iris_line_exact']['n']}. "
        f"Spend {money(memo['spent_usd'])}.")
    add("")
    add("## Guardrails")
    add("")
    notes = json.loads((RESULTS / "retrieval_notes.json").read_text(encoding="utf-8"))
    dropped = int(notes["new_study_variant_drops"])
    sealed_dropped = int(notes["heldout_variant_drops"])
    filt = json.loads((RESULTS / "caldera_filter_report.json").read_text(encoding="utf-8"))
    out_of_list = post["out_of_list"]
    review = replay["review"]
    add(f"Answers outside the candidate list were stored as blank and counted wrong: {out_of_list['no_papers']['n']} with no papers, {out_of_list['with_papers']['n']} with papers. All of them were blood-protein questions, and in every one the model named the measured protein.")
    add(f"Abstracts containing the lead variant's id or position were dropped when the pools were built: {dropped} on the new-study questions and {sealed_dropped} on the held-out chromosomes. Retrieval is LanceDB vector search with gene and date prefilters. On neighborhood A00407 the same 25 papers come back as in the locked run, with two swapped, because those scores are equal to 7 decimals. The locked run did not save a separate count of how many abstracts the date alone removed.")
    add(f"Statistical-model sets dropped for a protein-changing leak before the run: {filt['dropped_either']} (coding probability over 0.5 on {filt['dropped_coding_pip_over_0.5']}, the protein-changing lead inside the set on {filt['dropped_same_signal']}).")
    add(f"Questions flagged for a person to review, without changing the pick: {review['uk_practice_loco']['flagged']['n']} of 1,348 on the practice set "
        f"({review['uk_practice_loco']['flagged']['correct']} correct, {review['uk_practice_loco']['flagged']['percent']} percent) against "
        f"{review['uk_practice_loco']['unflagged']['percent']} percent on the rest. "
        f"On the new main list with no papers, {review['exam_main_no_papers']['flagged']['n']} flagged "
        f"({review['exam_main_no_papers']['flagged']['percent']} percent correct) against {review['exam_main_no_papers']['unflagged']['percent']} percent on the rest.")
    add("One kept input, TREM2, has a protein-changing probability of 0.41. The locked run kept it because the drop line was 0.5, and both methods named TREM2, so it is one of the 24. A new run uses that same line and does not look at the answer, so it keeps the question too.")
    add("")
    spend_file = json.loads((RESULTS / "api_spend.json").read_text(encoding="utf-8"))
    spends = [float(item["usd"]) for item in spend_file["parts"]]
    add(f"Total API spend recorded for the project: {money(sum(spends))}.")
    add("")
    add("On the new main list, adding papers changed the language model by 3 questions, and redrawing the questions puts that change between a loss of 11 and a gain of 18. The same model, with no papers, beats the nearest gene on diseases and other traits and loses to it on blood protein levels. The combining rule's paper change on the 38 new questions that have both inputs was zero. On the 18 held-out chromosomes, papers changed the rule by 20 questions, and those studies may have been in the training data. Counting the same protein gene once, instead of once per protein trait, is what keeps the new result from being dominated by one repeated fact.")
    diagnostic = RESULTS / "diagnostic_summary.json"
    if diagnostic.exists():
        item = json.loads(diagnostic.read_text(encoding="utf-8"))
        protein = item["slices"]["blood_protein"]
        other = item["slices"]["diseases_and_other"]
        lipids = item["slices"]["lipids"]
        main_slice = item["slices"]["main"]
        give = item["slices"]["giveaway"]
        near = item["protein_vs_nearest"]["new"]
        sure = item["sureness_main"]
        add("")
        add("## Diagnostic run after the locked test, on questions already seen")
        add("")
        add("One sentence was added to the instructions, telling the model that the measured protein's gene may be elsewhere and that the question is which listed gene changes the measurement. No papers. The locked scores above are unchanged. Cost "
            f"{money(item['cost_usd'])}.")
        add(f"Blood protein levels, main set: {protein['locked_correct']} to {protein['new_correct']} of {protein['n']}. "
            f"Fixed {protein['fixed']}, broken {protein['broken']}, net {protein['net']:+d}. "
            f"Redrawing the questions puts that change between {protein['net_question_interval'][0]:+d} and {protein['net_question_interval'][1]:+d}.")
        add(f"Diseases and other traits: {other['locked_correct']} to {other['new_correct']} of {other['n']}, net {other['net']:+d}. "
            f"Lipids: {lipids['locked_correct']} of {lipids['n']} either way. "
            f"Main set: {main_slice['locked_correct']} to {main_slice['new_correct']} of {main_slice['n']}, net {main_slice['net']:+d}. "
            f"Giveaway questions, where the protein's own gene is the right answer: {give['locked_correct']} of {give['n']} either way.")
        add(f"Answers outside the list fell from {item['out_of_list']['locked_all']} to {item['out_of_list']['new_all']} overall, "
            f"and from {item['out_of_list']['locked_main_protein']} to {item['out_of_list']['new_main_protein']} on main-set protein questions. "
            f"On those protein questions the model was the only correct method on {near['only_model']} and the nearest gene was the only correct method on {near['only_nearest']}. "
            f"The nearest gene, at {near['nearest_correct']} of {near['n']}, is still ahead of the model, at {near['model_correct']}.")
        add(f"On the main set, the median sureness moved from {sure['locked']['median']:.2f} to {sure['new']['median']:.2f}. "
            f"Answers at or above 0.85, the line the combining rule uses, moved from {sure['locked']['at_or_above_0.85']} to {sure['new']['at_or_above_0.85']} of {sure['locked']['n']}. "
            "The combining rule was not re-run.")
    OUT.write_text("\n".join(lines_out) + "\n", encoding="utf-8")
    print(f"wrote {OUT} spend {sum(spends):.2f}")


if __name__ == "__main__":
    main()

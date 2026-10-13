"""Replay the locked scores from saved answers. No model calls."""

from __future__ import annotations

import json
from pathlib import Path

from pipeline.decide import LINE, final_decision
from pipeline.paths import OLDER, QUESTIONS, RESULTS
from pipeline.questions import WORKSHOP, answer_key, exam_rows, load_answers
from pipeline.uk_table import choose_fixed, hits, load_rows

SEALED = OLDER / "heldout_paper_answers.jsonl"
CHOOSER = RESULTS / "uk_chooser_summary.json"
OUT = RESULTS / "replay_numbers.json"
DIAGNOSTIC = RESULTS / "diagnostic_nopaper.jsonl"


def load_sealed() -> dict[str, dict]:
    found = {}
    for line in SEALED.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            found[str(row["locus"])] = row
    return found


def uk_state(row: dict, line: float, paper: dict | None) -> dict:
    state = {
        "mode": "combined",
        "papers": "on" if paper else "off",
        "gpt_gene": row["gpt"],
        "gpt_confidence": row["gpt_conf"],
        "caldera_gene": row["cal"],
        "caldera_probability": row["cal_p"],
        "sureness_line": line,
    }
    if paper:
        state["paper_gene"] = paper.get("predicted") or ""
        state["paper_confidence"] = float(paper.get("confidence") or 0)
    return state


def correct(gene: str, label: str) -> bool:
    return bool(gene) and gene.upper() == label.upper()


def review_split(flags: list[tuple[bool, bool]]) -> dict:
    """flags are (needs_review, correct)."""
    flagged = [item for item in flags if item[0]]
    clear = [item for item in flags if not item[0]]

    def rate(items: list[tuple[bool, bool]]) -> dict:
        n = len(items)
        hits = sum(1 for _flag, good in items if good)
        return {"n": n, "correct": hits, "percent": round(100 * hits / n, 2) if n else None}

    return {"flagged": rate(flagged), "unflagged": rate(clear)}


def main() -> dict:
    rows = load_rows()
    summary = json.loads(CHOOSER.read_text(encoding="utf-8"))
    lines = {int(item["chrom"]): float(item["threshold"]) for item in summary["loco_thresholds_caldera"]}
    sealed = load_sealed()

    loco_hits = 0
    fixed_hits = 0
    near_loco = hits(rows, choose_fixed("near", 0.85))
    for row in rows:
        pick, source, _review = final_decision(uk_state(row, lines[row["chrom"]], None))
        if source == "agreement" or source == "gpt":
            loco_hits += int(row["gpt_hit"])
        else:
            loco_hits += int(row["cal_hit"])
        pick85, source85, _review = final_decision(uk_state(row, LINE, None))
        if source85 in ("agreement", "gpt"):
            fixed_hits += int(row["gpt_hit"])
        else:
            fixed_hits += int(row["cal_hit"])

    held = [row for row in rows if row["chrom"] not in WORKSHOP]
    held_no = 0
    held_yes = 0
    for row in held:
        _pick, source, _review = final_decision(uk_state(row, LINE, None))
        held_no += int(row["gpt_hit"] if source in ("agreement", "gpt") else row["cal_hit"])
        paper = sealed.get(str(row["locus"])) if not row["agree"] else None
        _pick, source, _review = final_decision(uk_state(row, LINE, paper))
        if source in ("agreement", "gpt"):
            if paper:
                held_yes += int(bool(paper.get("correct")))
            else:
                held_yes += int(row["gpt_hit"])
        else:
            held_yes += int(row["cal_hit"])

    full, main, exam_b = exam_rows()
    nopaper = load_answers(RESULTS / "new_nopaper.jsonl")
    papers = load_answers(RESULTS / "new_papers.jsonl")
    caldera = {}
    with (RESULTS / "new_caldera_genes.tsv").open(encoding="utf-8", newline="") as handle:
        import csv

        best: dict[str, tuple[float, str]] = {}
        for item in csv.DictReader(handle, delimiter="\t"):
            value = float(item["caldera"])
            if item["locus"] not in best or value > best[item["locus"]][0]:
                best[item["locus"]] = (value, item["gene"])
        caldera = {locus: (gene, prob) for locus, (prob, gene) in best.items()}

    def exam_counts(rows_in: list[dict]) -> dict:
        base = papers_n = nearest = 0
        rule_no = rule_yes = 0
        gpt_b = gpt_p = 0
        cal_n = 0
        for row in rows_in:
            key = answer_key(row)
            old = nopaper[key]
            new = papers[key]
            label = row["label_symbol"]
            base += int(bool(old.get("correct")))
            papers_n += int(bool(new.get("correct")))
            nearest += int(row["nearest_gene"].upper() == label.upper())
            gene, prob = caldera.get(row["neighborhood_id"], ("", 0.0))
            # Exam B ids in caldera_genes are the exam_a ids. Main representatives use those ids
            # only when the row itself is the exam_a row. Callers pass exam_a rows for Exam B.
            state = {
                "mode": "combined",
                "papers": "off",
                "gpt_gene": old.get("predicted") or "",
                "gpt_confidence": float(old.get("confidence") or 0),
                "caldera_gene": gene,
                "caldera_probability": prob,
                "sureness_line": LINE,
            }
            pick, _source, _review = final_decision(state)
            rule_no += int(correct(pick, label))
            state["papers"] = "on"
            state["paper_gene"] = new.get("predicted") or ""
            state["paper_confidence"] = float(new.get("confidence") or 0)
            pick, _source, _review = final_decision(state)
            rule_yes += int(correct(pick, label))
            gpt_b += int(bool(old.get("correct")))
            gpt_p += int(bool(new.get("correct")))
            cal_n += int(bool(gene) and gene.upper() == label.upper())
        return {
            "n": len(rows_in),
            "gpt_no_papers": gpt_b,
            "gpt_papers": gpt_p,
            "nearest": nearest,
            "rule_no_papers": rule_no,
            "rule_papers": rule_yes,
            "caldera": cal_n,
        }

    # Exam B rows in exam_a carry the locus ids CALDERA used.
    exam_b_counts = exam_counts(exam_b)
    main_gpt = sum(int(bool(nopaper[answer_key(row)].get("correct"))) for row in main)
    main_papers = sum(int(bool(papers[answer_key(row)].get("correct"))) for row in main)
    main_near = sum(int(row["nearest_gene"].upper() == row["label_symbol"].upper()) for row in main)

    # Review flags. UK practice uses the per-chromosome line. Exams use 0.85 or 0.70.
    uk_flags = []
    for row in rows:
        state = uk_state(row, lines[row["chrom"]], None)
        _pick, source, review = final_decision(state)
        good = row["gpt_hit"] if source in ("agreement", "gpt") else row["cal_hit"]
        uk_flags.append((review, bool(good)))
    held_flags = []
    for row in held:
        paper = sealed.get(str(row["locus"])) if not row["agree"] else None
        state = uk_state(row, LINE, paper)
        _pick, source, review = final_decision(state)
        if source in ("agreement", "gpt"):
            good = bool(paper.get("correct")) if paper else row["gpt_hit"]
        else:
            good = row["cal_hit"]
        held_flags.append((review, bool(good)))

    def llm_flags(rows_in: list[dict], which: str) -> list[tuple[bool, bool]]:
        flags = []
        answers = nopaper if which == "no" else papers
        for row in rows_in:
            answer = answers[answer_key(row)]
            state = {
                "mode": "llm_only",
                "papers": "off" if which == "no" else "on",
                "gpt_gene": answer.get("predicted") or "",
                "gpt_confidence": float(answer.get("confidence") or 0),
                "paper_gene": answer.get("predicted") or "",
                "paper_confidence": float(answer.get("confidence") or 0),
            }
            _pick, _source, review = final_decision(state)
            flags.append((review, bool(answer.get("correct"))))
        return flags

    from pipeline.guardrails import input_refusal

    import csv

    rows_by_locus: dict[str, list[dict]] = {}
    inputs = RESULTS / "caldera_inputs"
    for path in inputs.glob("*_cs.tsv") if inputs.exists() else []:
        for item in csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"):
            rows_by_locus.setdefault(item["locus"], []).append(
                {
                    "locus": item["locus"],
                    "c_gene": item.get("c_gene") or "",
                    "pip": float(item.get("pip") or 0),
                }
            )
    refused = []
    for row in exam_b:
        reason = input_refusal(
            [row["label_symbol"]],
            row["earliest_public_date"],
            rows_by_locus.get(row["neighborhood_id"], []),
            row["neighborhood_id"],
        )
        if reason:
            refused.append({"id": row["neighborhood_id"], "label": row["label_symbol"], "reason": reason})

    payload = {
        "uk_loco_rule": loco_hits,
        "uk_fixed_0.85": fixed_hits,
        "uk_nearest_partner_fixed_0.85": near_loco,
        "uk_nearest_partner_loco": summary["loco_threshold_gpt_nearest"],
        "heldout_rule_no_papers": held_no,
        "heldout_rule_papers": held_yes,
        "heldout_n": len(held),
        "exam_main_gpt_no_papers": main_gpt,
        "exam_main_gpt_papers": main_papers,
        "exam_main_nearest": main_near,
        "exam_main_n": len(main),
        "exam_b": exam_b_counts,
        "full_n": len(full),
        "review": {
            "uk_practice_loco": review_split(uk_flags),
            "heldout_0.85_with_papers": review_split(held_flags),
            "exam_main_no_papers": review_split(llm_flags(main, "no")),
            "exam_main_papers": review_split(llm_flags(main, "papers")),
        },
        "coding_pip_guardrail_refused": refused,
    }
    diagnostic = load_answers(DIAGNOSTIC)
    protein_rows = [row for row in main if row["trait_group"] == "blood_protein"]
    diagnostic_protein = sum(1 for row in protein_rows if diagnostic[answer_key(row)].get("correct"))
    diagnostic_main = sum(1 for row in main if diagnostic[answer_key(row)].get("correct"))

    expected = {
        "uk_loco_rule": 869,
        "uk_heldout_no": 738,
        "uk_heldout_yes": 758,
        "main_no": 200,
        "main_yes": 203,
        "main_near": 186,
        "b_rule": 24,
        "b_gpt_no": 23,
        "b_gpt_yes": 24,
        "b_cal": 12,
        "b_near": 14,
    }
    payload["match"] = {
        "uk_loco_869": loco_hits == expected["uk_loco_rule"],
        "heldout_738": held_no == expected["uk_heldout_no"],
        "heldout_758": held_yes == expected["uk_heldout_yes"],
        "main_200": main_gpt == expected["main_no"],
        "main_203": main_papers == expected["main_yes"],
        "main_186": main_near == expected["main_near"],
        "b_rule_24": exam_b_counts["rule_no_papers"] == expected["b_rule"] and exam_b_counts["rule_papers"] == expected["b_rule"],
        "b_gpt_23_24": exam_b_counts["gpt_no_papers"] == expected["b_gpt_no"] and exam_b_counts["gpt_papers"] == expected["b_gpt_yes"],
        "b_cal_12": exam_b_counts["caldera"] == expected["b_cal"],
        "b_near_14": exam_b_counts["nearest"] == expected["b_near"],
        "diagnostic_protein_83": diagnostic_protein == 83,
        "diagnostic_main_212": diagnostic_main == 212,
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload["match"], indent=2))
    print(
        f"uk {loco_hits} and {fixed_hits} of 1348; nearest partner {summary['loco_threshold_gpt_nearest']}"
    )
    print(f"heldout {held_no} and {held_yes} of {len(held)}")
    print(f"main {main_gpt}, {main_papers}, {main_near} of {len(main)}")
    print(
        "exam_b rule",
        exam_b_counts["rule_no_papers"],
        exam_b_counts["rule_papers"],
        "gpt",
        exam_b_counts["gpt_no_papers"],
        exam_b_counts["gpt_papers"],
        "caldera",
        exam_b_counts["caldera"],
        "nearest",
        exam_b_counts["nearest"],
        "of",
        len(exam_b),
    )
    print("diagnostic_protein", diagnostic_protein, "of 179; diagnostic_main", diagnostic_main, "of 429")
    if not all(payload["match"].values()):
        raise SystemExit("replay numbers do not match")
    return payload


if __name__ == "__main__":
    main()

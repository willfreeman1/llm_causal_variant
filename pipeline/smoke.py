"""Live check of 20 main-list questions, with papers. Cap $5. Does not replace any score."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pipeline.graph import run_question
from pipeline.measurement import measurement
from pipeline.paths import LOCAL, QUESTIONS, RESULTS
from pipeline.questions import answer_key, load_answers, load_query_vectors, symbols

CAP = 5.0
SEED = 20261009
N = 20
OUT = LOCAL / "ukb" / "final_exam" / "smoke.json"


def main() -> None:
    main_rows = []
    import csv

    with (QUESTIONS / "new_main.tsv").open(encoding="utf-8", newline="") as handle:
        main_rows = list(csv.DictReader(handle, delimiter="\t"))
    answers = load_answers(RESULTS / "new_papers.jsonl")
    vectors = load_query_vectors(LOCAL / "ukb" / "final_exam" / "query_vectors.npz")
    by_key = {}
    for key, answer in answers.items():
        by_key[key] = answer
    rng = np.random.default_rng(SEED)
    chosen = [main_rows[int(i)] for i in rng.choice(len(main_rows), size=N, replace=False)]
    trace = LOCAL / "ukb" / "final_exam" / "traces" / "smoke.jsonl"
    if trace.exists():
        trace.unlink()
    matched = 0
    spent = 0.0
    details = []
    for row in chosen:
        answer = by_key[answer_key(row)]
        result = run_question(
            {
                "question_id": row["neighborhood_id"],
                "mode": "llm_only",
                "papers": "on",
                "replay": False,
                "trait": row["trait_name"],
                "gene_string": row["candidate_genes"],
                "cutoff": row["earliest_public_date"],
                "date_mode": "study",
                "lead_rsid": row.get("noncoding_rsid") or "",
                "lead_position": str(row.get("noncoding_pos") or ""),
                "measurement": measurement(row),
                "query_vector": vectors[answer["id"]].tolist(),
                "spend_cap": CAP,
                "trace_path": str(trace),
            }
        )
        live = (result.get("paper_gene") or "").upper()
        saved = (answer.get("predicted") or "").upper()
        same = live == saved
        matched += int(same)
        spent += float(result.get("cost_usd") or 0)
        details.append(
            {
                "id": row["neighborhood_id"],
                "live": live,
                "saved": saved,
                "same": same,
                "cost_usd": result.get("cost_usd") or 0,
            }
        )
        print(f"SMOKE {row['neighborhood_id']} same {same} spent {spent:.3f}", flush=True)
        if spent >= CAP:
            raise SystemExit(f"smoke spend reached ${spent:.2f}")
    payload = {"n": N, "matched_saved": matched, "cost_usd": spent, "cap_usd": CAP, "seed": SEED, "details": details}
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"matched_saved": matched, "n": N, "cost_usd": round(spent, 4)}))


if __name__ == "__main__":
    main()

"""One diagnostic pass of the new-study questions with one extra prompt line.

No papers. Does not change the locked prompt, the locked answer files, or the
locked spend file. Stops at $5. Resumes from the diagnostic answer file.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "spikes" / "shringarpure_match"))

import run_spike
from run_final_exam import rows_to_call, symbols_of

EXAM = ROOT / "data" / "ukb" / "final_exam"
OUT = EXAM / "results_nopaper_diagnostic_prompt.jsonl"
LEDGER = EXAM / "diagnostic_prompt_spend.json"
SUMMARY = EXAM / "diagnostic_prompt_summary.json"
LOCKED = EXAM / "results_nopaper.jsonl"
CAP = 5.0
BOOT = 2000
SEED = 20261009
HINT = (
    "The phenotype may be the level of a specific protein or other molecule. "
    "The gene that encodes that molecule may be elsewhere in the genome and may not be in the list. "
    "The question is which of the listed genes most likely changes the phenotype."
)


def diagnostic_system() -> str:
    marker = "based on literature evidence.\n"
    if run_spike.SYSTEM.count(marker) != 1:
        raise SystemExit("locked system prompt does not contain the expected line")
    system = run_spike.SYSTEM.replace(marker, marker + HINT + "\n", 1)
    if HINT not in system or system == run_spike.SYSTEM:
        raise SystemExit("diagnostic prompt was not inserted")
    return system


def read_spent() -> float:
    if not LEDGER.exists():
        return 0.0
    return float(json.loads(LEDGER.read_text(encoding="utf-8")).get("usd") or 0)


def add_spent(amount: float) -> float:
    total = read_spent() + float(amount)
    LEDGER.write_text(json.dumps({"usd": total}, indent=2), encoding="utf-8")
    return total


def load_key() -> str:
    for line in (ROOT / ".env").read_text(encoding="utf-8").splitlines():
        if line.startswith("OPENAI_API_KEY="):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise SystemExit("OPENAI_API_KEY is missing")


def done_ids() -> set[str]:
    found = set()
    if not OUT.exists():
        return found
    for line in OUT.read_text(encoding="utf-8").splitlines():
        if line.strip():
            found.add(json.loads(line)["id"])
    return found


def reason_of(text: str) -> str:
    raw = text.strip()
    start = raw.find("{")
    end = raw.rfind("}")
    if start < 0 or end <= start:
        return ""
    try:
        parsed = json.loads(raw[start : end + 1])
    except json.JSONDecodeError:
        return ""
    return str(parsed.get("reason") or "").strip()


def chat(api_key: str, system: str, user: str) -> dict:
    body = {
        "model": run_spike.MODEL,
        "temperature": 0,
        "max_completion_tokens": 500,
        "reasoning_effort": "none",
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
    }
    last = ""
    payload = None
    for attempt in range(5):
        request = urllib.request.Request(
            "https://api.openai.com/v1/chat/completions",
            data=json.dumps(body).encode("utf-8"),
            headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as exc:
            last = exc.read().decode("utf-8", errors="replace")[:500]
            if exc.code == 400 and "temperature" in last and "temperature" in body:
                body.pop("temperature", None)
                continue
            if exc.code == 400 and "reasoning_effort" in last and "reasoning_effort" in body:
                body.pop("reasoning_effort", None)
                continue
            if exc.code in (429, 500, 502, 503):
                time.sleep(2 + attempt * 3)
                continue
            raise SystemExit(f"API {exc.code}: {last}") from exc
    if payload is None:
        raise SystemExit(f"API failed after retries: {last}")
    usage = payload.get("usage") or {}
    prompt = int(usage.get("prompt_tokens") or 0)
    completion = int(usage.get("completion_tokens") or 0)
    text = payload["choices"][0]["message"].get("content") or ""
    return {
        "text": text,
        "prompt_tokens": prompt,
        "completion_tokens": completion,
        "cost_usd": prompt * run_spike.INPUT_USD_PER_TOKEN + completion * run_spike.OUTPUT_USD_PER_TOKEN,
    }


def run() -> None:
    system = diagnostic_system()
    print("PROMPT_OK", flush=True)
    api_key = load_key()
    rows = rows_to_call()
    if len(rows) != 1087:
        raise SystemExit(f"expected 1087 questions, found {len(rows)}")
    finished = done_ids()
    todo = [row for row in rows if row["neighborhood_id"] not in finished]
    print(f"TODO {len(todo)} spent {read_spent():.4f}", flush=True)
    with OUT.open("a", encoding="utf-8") as handle:
        for index, row in enumerate(todo, start=1):
            if read_spent() >= CAP:
                raise SystemExit(f"CAP before {row['neighborhood_id']} at {read_spent():.4f}")
            user = (
                "Identify the causal gene.\n"
                f"GWAS phenotype: {row['trait_name']}\n"
                f"Genes in locus: {row['candidate_genes']}\n"
            )
            result = chat(api_key, system, user)
            gene, confidence = run_spike.parse_gene(result["text"])
            allowed = {symbol.upper() for symbol in symbols_of(row)}
            predicted = gene.strip()
            in_list = predicted.upper() in allowed if predicted else False
            correct = bool(predicted) and in_list and predicted.upper() == row["label_symbol"].upper()
            total = add_spent(result["cost_usd"])
            handle.write(json.dumps({
                "id": row["neighborhood_id"],
                "chrom": row["noncoding_chrom"],
                "pos": row["noncoding_pos"],
                "rsid": row["noncoding_rsid"],
                "label": row["label_symbol"],
                "trait": row["trait_name"],
                "predicted": predicted if in_list else "",
                "raw_predicted": predicted,
                "confidence": confidence if confidence is not None else 0.0,
                "correct": correct,
                "hallucination": bool(predicted) and not in_list,
                "parsed": bool(predicted),
                "reason": reason_of(result["text"]),
                "cost_usd": result["cost_usd"],
                "prompt_tokens": result["prompt_tokens"],
                "completion_tokens": result["completion_tokens"],
            }) + "\n")
            handle.flush()
            if index % 50 == 0 or index == len(todo):
                print(f"DONE {index} of {len(todo)} ${total:.4f}", flush=True)
    score()


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
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def mcnemar(fixed: int, broken: int) -> float:
    total = fixed + broken
    if total == 0:
        return 1.0
    return float(binomtest(fixed, total, 0.5, alternative="two-sided").pvalue)


def paired(rows: list[dict], locked: dict, new: dict) -> dict:
    fixed = broken = locked_right = new_right = 0
    nets = []
    for row in rows:
        key = answer_key(row)
        old = bool(locked[key].get("correct"))
        cur = bool(new[key].get("correct"))
        locked_right += int(old)
        new_right += int(cur)
        change = int(cur) - int(old)
        nets.append(change)
        fixed += int((not old) and cur)
        broken += int(old and (not cur))
    draw = np.random.default_rng(SEED)
    chosen = draw.integers(0, len(nets), size=(BOOT, len(nets)))
    samples = np.asarray(nets)[chosen].sum(axis=1)
    low, high = np.percentile(samples, [2.5, 97.5])
    return {
        "n": len(rows),
        "locked_correct": locked_right,
        "new_correct": new_right,
        "fixed": fixed,
        "broken": broken,
        "net": fixed - broken,
        "net_question_interval": [int(low), int(high)],
        "mcnemar_p": mcnemar(fixed, broken),
    }


def against_nearest(rows: list[dict], answers: dict) -> dict:
    only_model = only_near = model_right = near_right = 0
    for row in rows:
        model = bool(answers[answer_key(row)].get("correct"))
        near = row["nearest_gene"].upper() == row["label_symbol"].upper()
        model_right += int(model)
        near_right += int(near)
        only_model += int(model and not near)
        only_near += int(near and not model)
    return {
        "n": len(rows),
        "model_correct": model_right,
        "nearest_correct": near_right,
        "only_model": only_model,
        "only_nearest": only_near,
        "mcnemar_p": mcnemar(only_model, only_near),
    }


def sureness(rows: list[dict], answers: dict) -> dict:
    values = np.asarray([float(answers[answer_key(row)].get("confidence") or 0) for row in rows], dtype=float)
    return {
        "n": int(values.size),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "p10": float(np.percentile(values, 10)),
        "p90": float(np.percentile(values, 90)),
        "at_or_above_0.85": int((values >= 0.85).sum()),
        "below_0.70": int((values < 0.70).sum()),
    }


def out_of_list(rows: list[dict], answers: dict) -> int:
    return sum(1 for row in rows if answers[answer_key(row)].get("hallucination"))


def examples(protein: list[dict], locked: dict, new: dict) -> dict:
    gained = []
    lost = []
    for row in protein:
        key = answer_key(row)
        old = locked[key]
        cur = new[key]
        item = {
            "id": row["neighborhood_id"],
            "trait": row["trait_name"],
            "label": row["label_symbol"],
            "locked": old.get("raw_predicted") or "",
            "new": cur.get("raw_predicted") or "",
            "reason": cur.get("reason") or "",
        }
        if (not old.get("correct")) and cur.get("correct") and len(gained) < 4:
            gained.append(item)
        if old.get("correct") and (not cur.get("correct")) and len(lost) < 4:
            lost.append(item)
    return {"gained": gained, "lost": lost}


def score() -> None:
    locked = load_jsonl(LOCKED)
    new = load_jsonl(OUT)
    main = read_tsv(EXAM / "exam_main.tsv")
    give = read_tsv(EXAM / "exam_giveaway.tsv")
    protein = [row for row in main if row["trait_group"] == "blood_protein"]
    other = [row for row in main if row["trait_group"] not in ("blood_protein", "blood_lipid")]
    lipids = [row for row in main if row["trait_group"] == "blood_lipid"]
    sizes = {
        "protein": len(protein),
        "other": len(other),
        "lipids": len(lipids),
        "main": len(main),
        "giveaway": len(give),
        "answers": len(new),
    }
    if sizes != {"protein": 179, "other": 221, "lipids": 29, "main": 429, "giveaway": 246, "answers": 1087}:
        raise SystemExit(f"unexpected slice sizes {sizes}")
    all_rows = read_tsv(EXAM / "exam_a.tsv")
    called = {answer_key(row) for row in all_rows}
    called.update(answer_key(row) for row in main)
    if len(called) != 1087:
        raise SystemExit(f"call keys {len(called)}")
    locked_hall_main = out_of_list(main, locked)
    payload = {
        "cost_usd": read_spent(),
        "slices": {
            "blood_protein": paired(protein, locked, new),
            "diseases_and_other": paired(other, locked, new),
            "lipids": paired(lipids, locked, new),
            "main": paired(main, locked, new),
            "giveaway": paired(give, locked, new),
        },
        "out_of_list": {
            "locked_all": sum(1 for row in locked.values() if row.get("hallucination")),
            "new_all": sum(1 for row in new.values() if row.get("hallucination")),
            "locked_main": locked_hall_main,
            "new_main": out_of_list(main, new),
            "locked_main_protein": out_of_list(protein, locked),
            "new_main_protein": out_of_list(protein, new),
        },
        "protein_vs_nearest": {
            "locked": against_nearest(protein, locked),
            "new": against_nearest(protein, new),
        },
        "sureness_main": {"locked": sureness(main, locked), "new": sureness(main, new)},
        "sureness_protein": {"locked": sureness(protein, locked), "new": sureness(protein, new)},
        "examples": examples(protein, locked, new),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    if "--show-prompt" in sys.argv:
        print(diagnostic_system())
    elif "--score" in sys.argv:
        score()
    else:
        run()

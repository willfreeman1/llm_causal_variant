"""Check whether the model can recite rows from a saved Open Targets table.

Asks gpt-5.4 to recite file artifacts it could not know from biology:
line number -> exact cells, and a keyed lookup of curator metadata.
Iris line recital is the positive control. Out-of-range lines and fake
variant ids are the negative controls.

Loads OPENAI_API_KEY from the repo .env. Does not print the key.
Stops if estimated spend reaches MAX_USD.
"""

from __future__ import annotations

import csv
import json
import random
import re
import time
import urllib.error
import urllib.request
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENV_PATH = ROOT / ".env"
DATA = Path(os.environ['EARLY_CHECK_DATA'])
OUT = Path(__file__).resolve().parent / "memorization_results.json"

MODEL = "gpt-5.4"
INPUT_USD_PER_TOKEN = 2.50 / 1_000_000
OUTPUT_USD_PER_TOKEN = 15.00 / 1_000_000
MAX_USD = 8.0
SEED = 20261006

GOLD_PATH = DATA / "gold.tsv"
IRIS_PATH = DATA / "iris.data"

LINE_FIELDS = [
    "gold_standard_info.gene_id",
    "gold_standard_info.evidence.class",
    "gold_standard_info.evidence.confidence",
    "gold_standard_info.evidence.curated_by",
    "gold_standard_info.evidence.pubmed_id",
    "metadata.set_label",
    "metadata.date_added",
]
KEYED_FIELDS = [
    "gold_standard_info.evidence.class",
    "gold_standard_info.evidence.confidence",
    "gold_standard_info.evidence.curated_by",
    "gold_standard_info.evidence.pubmed_id",
    "metadata.set_label",
    "metadata.date_added",
]

SYSTEM = (
    "You are being tested on whether a specific file is in your training data. "
    "Reply with exact stored cell values only. If you do not know a value exactly, "
    "write UNKNOWN. Do not guess, do not infer, and do not explain. "
    "Output one field per line as field=value."
)


def load_api_key() -> str:
    for line in ENV_PATH.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == "OPENAI_API_KEY":
            return value.strip().strip('"').strip("'")
    raise SystemExit("OPENAI_API_KEY not found in .env")


def norm(value: str) -> str:
    value = value.strip().strip('"').strip("'")
    value = re.sub(r"\s+", " ", value)
    return value.casefold()


EMPTY_TOKENS = {"", "empty", "blank", "none", "null", "na", "n/a", "(empty)", "nan"}


def is_unknown(value: str) -> bool:
    return norm(value) in {"unknown", "i don't know", "idk", "not in training data"}


def values_match(pred: str, gold: str) -> bool:
    if is_unknown(pred):
        return False
    p, g = norm(pred), norm(gold)
    if g in EMPTY_TOKENS:
        return p in EMPTY_TOKENS
    return p == g


def parse_fields(text: str, fields: list[str]) -> dict[str, str]:
    found = {field: "UNKNOWN" for field in fields}
    for line in text.splitlines():
        if "=" not in line:
            continue
        left, right = line.split("=", 1)
        left = left.strip()
        for field in fields:
            if left == field or left.endswith(field) or field.endswith(left):
                found[field] = right.strip()
    return found


def call_model(api_key: str, user: str) -> dict:
    body = {
        "model": MODEL,
        "temperature": 0,
        "max_completion_tokens": 300,
        "reasoning_effort": "none",
        "messages": [
            {"role": "system", "content": SYSTEM},
            {"role": "user", "content": user},
        ],
    }
    req = urllib.request.Request(
        "https://api.openai.com/v1/chat/completions",
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[:800]
        # Retry once without temperature if the model rejects it.
        if exc.code == 400 and ("temperature" in detail or "reasoning_effort" in detail):
            body.pop("temperature", None)
            if "reasoning_effort" in detail:
                body.pop("reasoning_effort", None)
            req = urllib.request.Request(
                "https://api.openai.com/v1/chat/completions",
                data=json.dumps(body).encode("utf-8"),
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json",
                },
                method="POST",
            )
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    payload = json.loads(resp.read().decode("utf-8"))
            except urllib.error.HTTPError as exc2:
                detail2 = exc2.read().decode("utf-8", errors="replace")[:800]
                raise SystemExit(f"API {exc2.code}: {detail2}") from exc2
        else:
            raise SystemExit(f"API {exc.code}: {detail}") from exc
    usage = payload.get("usage") or {}
    message = payload["choices"][0]["message"]
    text = message.get("content") or ""
    prompt_tokens = usage.get("prompt_tokens") or 0
    completion_tokens = usage.get("completion_tokens") or 0
    cost = prompt_tokens * INPUT_USD_PER_TOKEN + completion_tokens * OUTPUT_USD_PER_TOKEN
    return {
        "text": text,
        "prompt_tokens": prompt_tokens,
        "completion_tokens": completion_tokens,
        "cost_usd": cost,
    }


def field_block(fields: list[str]) -> str:
    return "\n".join(fields)


def line_prompt(n: int) -> str:
    return (
        "File: gold_standards/processed/gwas_gold_standards.191108.tsv\n"
        "Repository: opentargets/genetics-gold-standards\n"
        f"Data row number: {n}\n"
        "Row numbers start at 1 for the first data row after the header.\n"
        "Return exactly these fields as stored in that row:\n"
        f"{field_block(LINE_FIELDS)}"
    )


def keyed_prompt(row: dict) -> str:
    return (
        "File: gold_standards/processed/gwas_gold_standards.191108.tsv\n"
        "Repository: opentargets/genetics-gold-standards\n"
        "Identify the row with these values and return the other fields as stored:\n"
        f"sentinel_variant.rsid={row['sentinel_variant.rsid']}\n"
        f"sentinel_variant.locus_GRCh38.chromosome={row['sentinel_variant.locus_GRCh38.chromosome']}\n"
        f"sentinel_variant.locus_GRCh38.position={row['sentinel_variant.locus_GRCh38.position']}\n"
        f"association_info.gwas_catalog_id={row['association_info.gwas_catalog_id']}\n"
        f"trait_info.standard_trait_name={row['trait_info.standard_trait_name']}\n"
        "Return exactly:\n"
        f"{field_block(KEYED_FIELDS)}"
    )


def iris_prompt(n: int) -> str:
    return (
        "File: the canonical UCI iris.data file (comma-separated, no header, "
        "columns sepal length, sepal width, petal length, petal width, class).\n"
        f"Line number: {n} (1-indexed).\n"
        "Return exactly:\n"
        "line=<the full comma-separated line>"
    )


def majority(rows: list[dict], field: str) -> str:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row[field]] = counts.get(row[field], 0) + 1
    return max(counts, key=counts.get)


def score_fields(pred: dict[str, str], gold: dict, fields: list[str], modal: dict[str, str]) -> dict:
    per = {}
    for field in fields:
        match = values_match(pred.get(field, "UNKNOWN"), gold[field])
        rare = norm(gold[field]) != norm(modal[field])
        per[field] = {
            "match": match,
            "rare_gold": rare,
            "rare_hit": match and rare,
            "unknown": is_unknown(pred.get(field, "UNKNOWN")),
        }
    return per


def summarize(records: list[dict], fields: list[str]) -> dict:
    n = len(records)
    if n == 0:
        return {"n": 0}
    joint = 0
    field_match = {field: 0 for field in fields}
    rare_n = {field: 0 for field in fields}
    rare_hit = {field: 0 for field in fields}
    unknown = {field: 0 for field in fields}
    for rec in records:
        per = rec["per_field"]
        if all(per[field]["match"] for field in fields):
            joint += 1
        for field in fields:
            field_match[field] += int(per[field]["match"])
            rare_n[field] += int(per[field]["rare_gold"])
            rare_hit[field] += int(per[field]["rare_hit"])
            unknown[field] += int(per[field]["unknown"])
    return {
        "n": n,
        "joint_exact": joint,
        "joint_exact_rate": joint / n,
        "field_match_rate": {field: field_match[field] / n for field in fields},
        "rare_value_hit_rate": {
            field: (rare_hit[field] / rare_n[field] if rare_n[field] else None) for field in fields
        },
        "rare_value_counts": rare_n,
        "unknown_rate": {field: unknown[field] / n for field in fields},
    }


def main() -> None:
    api_key = load_api_key()
    rows = list(csv.DictReader(GOLD_PATH.open(encoding="utf-8"), delimiter="\t"))
    iris_lines = [ln.strip() for ln in IRIS_PATH.read_text(encoding="utf-8").splitlines() if ln.strip()]
    rng = random.Random(SEED)

    with_pmid = [i for i, row in enumerate(rows) if row["gold_standard_info.evidence.pubmed_id"].strip()]
    without_pmid = [i for i in range(len(rows)) if i not in set(with_pmid)]
    chosen = rng.sample(with_pmid, 15) + rng.sample(without_pmid, 10)
    rng.shuffle(chosen)

    modal = {field: majority(rows, field) for field in LINE_FIELDS}
    keyed_modal_rate = {}
    for field in KEYED_FIELDS:
        mode = modal[field]
        keyed_modal_rate[field] = sum(1 for row in rows if row[field] == mode) / len(rows)

    tests: list[dict] = []
    for idx in chosen:
        tests.append({"kind": "ot_line", "index": idx, "row_number": idx + 1})
        tests.append({"kind": "ot_keyed", "index": idx, "row_number": idx + 1})
    for n in range(5001, 5011):
        tests.append({"kind": "ot_line_oor", "row_number": n})
    iris_idx = rng.sample(range(len(iris_lines)), 15)
    for i in iris_idx:
        tests.append({"kind": "iris_line", "row_number": i + 1})
    for n in range(401, 406):
        tests.append({"kind": "iris_line_oor", "row_number": n})

    used_rs = {row["sentinel_variant.rsid"] for row in rows}
    for k in range(10):
        fake_rs = f"rs9{k:08d}"
        while fake_rs in used_rs:
            fake_rs = "rs" + str(rng.randint(10**9, 10**10 - 1))
        tests.append(
            {
                "kind": "fake_keyed",
                "row": {
                    "sentinel_variant.rsid": fake_rs,
                    "sentinel_variant.locus_GRCh38.chromosome": "22",
                    "sentinel_variant.locus_GRCh38.position": str(10_000_000 + k),
                    "association_info.gwas_catalog_id": f"GCSTFAKE{k:04d}",
                    "trait_info.standard_trait_name": f"Synthetic phenotype zeta-{k}",
                },
            }
        )

    spent = 0.0
    results = []
    for i, test in enumerate(tests, start=1):
        kind = test["kind"]
        if kind in {"ot_line", "ot_line_oor"}:
            prompt = line_prompt(test["row_number"])
        elif kind == "ot_keyed":
            prompt = keyed_prompt(rows[test["index"]])
        elif kind == "fake_keyed":
            prompt = keyed_prompt(test["row"])
        else:
            prompt = iris_prompt(test["row_number"])

        if spent >= MAX_USD:
            print(f"STOP at ${spent:.4f} before call {i}/{len(tests)}")
            break
        t0 = time.time()
        response = call_model(api_key, prompt)
        spent += response["cost_usd"]
        record = {
            "kind": kind,
            "row_number": test.get("row_number"),
            "response": response["text"],
            "prompt_tokens": response["prompt_tokens"],
            "completion_tokens": response["completion_tokens"],
            "cost_usd": response["cost_usd"],
            "seconds": round(time.time() - t0, 2),
        }
        if kind == "ot_line":
            gold = rows[test["index"]]
            pred = parse_fields(response["text"], LINE_FIELDS)
            record["per_field"] = score_fields(pred, gold, LINE_FIELDS, modal)
            record["pred"] = pred
        elif kind == "ot_keyed":
            gold = rows[test["index"]]
            pred = parse_fields(response["text"], KEYED_FIELDS)
            record["per_field"] = score_fields(pred, gold, KEYED_FIELDS, modal)
            record["pred"] = pred
            record["rsid"] = gold["sentinel_variant.rsid"]
        elif kind == "fake_keyed":
            pred = parse_fields(response["text"], KEYED_FIELDS)
            record["pred"] = pred
            record["all_unknown"] = all(is_unknown(pred[field]) for field in KEYED_FIELDS)
            record["any_specific"] = any(not is_unknown(pred[field]) for field in KEYED_FIELDS)
        elif kind == "iris_line":
            gold_line = iris_lines[test["row_number"] - 1]
            match = re.search(r"line=(.*)", response["text"])
            quoted = match.group(1).strip() if match else response["text"].strip()
            record["match"] = norm(quoted) == norm(gold_line)
            record["unknown"] = is_unknown(quoted) or quoted.upper() == "UNKNOWN"
            record["gold"] = gold_line
            record["pred_line"] = quoted
        elif kind == "iris_line_oor":
            record["unknown"] = "UNKNOWN" in response["text"].upper()
        elif kind == "ot_line_oor":
            pred = parse_fields(response["text"], LINE_FIELDS)
            record["pred"] = pred
            record["all_unknown"] = all(is_unknown(pred[field]) for field in LINE_FIELDS)
        results.append(record)
        flag = ""
        if "per_field" in record:
            hits = sum(1 for field in record["per_field"] if record["per_field"][field]["match"])
            flag = f" hits={hits}"
        elif "match" in record:
            flag = f" match={record['match']}"
        elif "all_unknown" in record:
            flag = f" all_unknown={record['all_unknown']}"
        print(
            f"{i:03d}/{len(tests)} {kind}{flag} "
            f"tok={response['prompt_tokens']}+{response['completion_tokens']} "
            f"call=${response['cost_usd']:.5f} total=${spent:.4f}",
            flush=True,
        )

    ot_line = [r for r in results if r["kind"] == "ot_line"]
    ot_keyed = [r for r in results if r["kind"] == "ot_keyed"]
    summary = {
        "model": MODEL,
        "seed": SEED,
        "n_gold_rows": len(rows),
        "spent_usd": round(spent, 6),
        "calls": len(results),
        "modal_field": {field: modal[field] for field in LINE_FIELDS},
        "modal_rate": keyed_modal_rate,
        "ot_line": summarize(ot_line, LINE_FIELDS),
        "ot_keyed": summarize(ot_keyed, KEYED_FIELDS),
        "ot_line_oor_all_unknown": {
            "n": sum(1 for r in results if r["kind"] == "ot_line_oor"),
            "all_unknown": sum(1 for r in results if r["kind"] == "ot_line_oor" and r.get("all_unknown")),
        },
        "iris_line_exact": {
            "n": sum(1 for r in results if r["kind"] == "iris_line"),
            "exact": sum(1 for r in results if r["kind"] == "iris_line" and r.get("match")),
        },
        "iris_oor_unknown": {
            "n": sum(1 for r in results if r["kind"] == "iris_line_oor"),
            "unknown": sum(1 for r in results if r["kind"] == "iris_line_oor" and r.get("unknown")),
        },
        "fake_keyed": {
            "n": sum(1 for r in results if r["kind"] == "fake_keyed"),
            "all_unknown": sum(1 for r in results if r["kind"] == "fake_keyed" and r.get("all_unknown")),
            "any_specific": sum(1 for r in results if r["kind"] == "fake_keyed" and r.get("any_specific")),
        },
    }
    OUT.write_text(json.dumps({"summary": summary, "results": results}, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

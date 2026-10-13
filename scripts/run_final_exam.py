"""One run of the final exam. Resumes. Stops at $40. Does not change prompts.

Gate checks run before any API call. See final_exam_setup.md, written 9 October 2026.
"""

from __future__ import annotations

import csv
import json
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.embed import BATCH, cache_dir, cap_text, chunk_paths, embed_request, embedded_pmids
from pipeline.measurement import is_blood_protein
from pipeline.model_api import INPUT_USD_PER_TOKEN, MODEL, OUTPUT_USD_PER_TOKEN, SYSTEM, parse_gene

EXAM = ROOT / "data" / "ukb" / "final_exam"
DB = ROOT / "data" / "ukb" / "pubmed" / "weeks_abstracts.sqlite"
SETUP = ROOT / "final_exam_setup.md"
CONTEXT = ROOT / "context.md"
SPEND_PATH = EXAM / "spend.json"
POOLS = EXAM / "pools.jsonl"
NO_PAPER = EXAM / "results_nopaper.jsonl"
PAPERS = EXAM / "results_papers.jsonl"
SUMMARY = EXAM / "results_summary.json"
CAP = 40.0
SHORT_USD = 0.002
PAPER_USD = 0.024
EMBED_USD = 0.02 / 1_000_000
TOP_N = 25
ABSTRACT_CAP = 1800
LINE = 0.85
BOOT = 2000
SEED = 20261009
MODEL_SPEND = {"embed_usd": 0.0, "gpt_usd": 0.0}


def load_key() -> str:
    for line in (ROOT / ".env").read_text(encoding="utf-8").splitlines():
        if line.startswith("OPENAI_API_KEY="):
            return line.split("=", 1)[1].strip()
    raise SystemExit("OPENAI_API_KEY is missing")


def load_spend() -> None:
    if SPEND_PATH.exists():
        saved = json.loads(SPEND_PATH.read_text(encoding="utf-8"))
        MODEL_SPEND.update(saved)


def save_spend() -> None:
    SPEND_PATH.write_text(json.dumps(MODEL_SPEND, indent=2), encoding="utf-8")


def spent() -> float:
    return MODEL_SPEND["embed_usd"] + MODEL_SPEND["gpt_usd"]


def gate() -> None:
    text = SETUP.read_text(encoding="utf-8")
    required = [
        "Written 9 October 2026",
        "## Which numbers count",
        "## How results are reported",
        "## Pooled tests, written down now",
        "The cap is **$40**",
        "TREM2",
    ]
    missing = [item for item in required if item not in text]
    if missing:
        raise SystemExit(f"setup is missing {missing}")
    main = list(csv.DictReader((EXAM / "exam_main.tsv").open(encoding="utf-8"), delimiter="\t"))
    if len(main) < 300:
        raise SystemExit(f"main set {len(main)} is under 300")
    counts: dict[str, int] = defaultdict(int)
    for row in main:
        counts[row["label_symbol"]] += 1
    top = max(counts.values())
    if top / len(main) > 0.05:
        raise SystemExit(f"one gene is {top} of {len(main)}")
    full = list(csv.DictReader((EXAM / "exam_a.tsv").open(encoding="utf-8"), delimiter="\t"))
    exam_b = [row for row in full if row["in_exam_b"] == "1"]

    def key(row: dict) -> tuple:
        return (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])

    have = {key(row) for row in main}
    absent = [row["neighborhood_id"] for row in exam_b if key(row) not in have]
    if absent:
        raise SystemExit(f"Exam B not in the main set: {absent[:5]}")
    if not (EXAM / "caldera_genes.tsv").exists():
        raise SystemExit("CALDERA output is missing")
    estimate = (EXAM / "cost_estimate.log").read_text(encoding="utf-8")
    if "total_usd 31.49" not in estimate:
        raise SystemExit("cost estimate is not the locked $31.49 figure")
    print(f"GATE_OK main {len(main)} max_gene {top} exam_b {len(exam_b)} estimate 31.49", flush=True)


def rows_to_call() -> list[dict]:
    full = list(csv.DictReader((EXAM / "exam_a.tsv").open(encoding="utf-8"), delimiter="\t"))
    main = list(csv.DictReader((EXAM / "exam_main.tsv").open(encoding="utf-8"), delimiter="\t"))

    def key(row: dict) -> tuple:
        return (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])

    have = {key(row) for row in full}
    extra = [row for row in main if key(row) not in have]
    print(f"CALLS full {len(full)} extra {len(extra)}", flush=True)
    return full + extra


def measurement(row: dict) -> str:
    trait = row["trait_name"].strip()
    mapped = (row.get("mapped_trait") or "").strip()
    if not is_blood_protein(trait, mapped) and (row.get("trait_group") == "blood_lipid" or row.get("is_lipid") == "1"):
        name = re.sub(r"\s*\(Model [12]\)", "", trait).strip()
        return f"{name}. Blood lipid measurement."
    left = re.sub(r"[^a-z0-9]+", " ", trait.lower()).strip()
    right = re.sub(r"[^a-z0-9]+", " ", mapped.lower()).strip()
    if not mapped or left == right:
        return f"{trait}."
    return f"{trait}. Catalog term: {mapped}."


def before(year, month, day, cutoff: str) -> bool:
    if not year or not cutoff:
        return False
    cy, cm, cd = (int(part) for part in cutoff.split("-"))
    year = int(year)
    month = int(month or 0)
    day = int(day or 0)
    if month <= 0:
        return year < cy
    if day <= 0:
        return (year, month) < (cy, cm)
    return (year, month, day) < (cy, cm, cd)


def symbols_of(row: dict) -> list[str]:
    return [part.strip().strip("{}") for part in row["candidate_genes"].split(",") if part.strip(" {}")]


def build_pools(rows: list[dict]) -> None:
    if POOLS.exists() and sum(1 for _ in POOLS.open(encoding="utf-8")) >= len(rows):
        print("POOLS_CACHED", flush=True)
        return
    symbols = sorted({gene for row in rows for gene in symbols_of(row)})
    conn = sqlite3.connect(DB)
    hits: dict[str, list[int]] = {symbol: [] for symbol in symbols}
    for start in range(0, len(symbols), 400):
        batch = symbols[start : start + 400]
        marks = ",".join("?" for _ in batch)
        for symbol, pmid in conn.execute(f"SELECT symbol, pmid FROM gene_hits WHERE symbol IN ({marks})", batch):
            hits[symbol].append(int(pmid))
    dated: dict[int, tuple] = {}
    texts: dict[int, str] = {}
    pmids = list({pmid for values in hits.values() for pmid in values})
    for start in range(0, len(pmids), 400):
        batch = pmids[start : start + 400]
        marks = ",".join("?" for _ in batch)
        for pmid, year, month, day, title, abstract in conn.execute(
            f"SELECT pmid, year, month, day, title, abstract FROM abstracts WHERE pmid IN ({marks})",
            batch,
        ):
            if not abstract:
                continue
            dated[int(pmid)] = (year, month, day)
            texts[int(pmid)] = f"{title or ''} {abstract}"
    print(f"POOL_TEXTS {len(texts)}", flush=True)
    with POOLS.open("w", encoding="utf-8") as handle:
        for index, row in enumerate(rows, start=1):
            cutoff = row["earliest_public_date"]
            variant = (row.get("noncoding_rsid") or "").lower()
            position = str(row.get("noncoding_pos") or "")
            kept = []
            dropped = 0
            seen = set()
            for gene in symbols_of(row):
                for pmid in hits.get(gene, []):
                    if pmid in seen:
                        continue
                    stamp = dated.get(pmid)
                    if not stamp or not before(*stamp, cutoff):
                        continue
                    seen.add(pmid)
                    text = texts.get(pmid, "").lower()
                    if (variant and variant in text) or (position and position in text):
                        dropped += 1
                        continue
                    kept.append(pmid)
            record = {
                "id": row["neighborhood_id"],
                "pmids": kept,
                "n_dropped_variant": dropped,
                "measurement": measurement(row),
            }
            handle.write(json.dumps(record) + "\n")
            if index % 100 == 0 or index == len(rows):
                print(f"POOL {index} of {len(rows)}", flush=True)


def load_pools() -> dict[str, dict]:
    found = {}
    for line in POOLS.read_text(encoding="utf-8").splitlines():
        if line.strip():
            record = json.loads(line)
            found[record["id"]] = record
    return found


def embed_missing(api_key: str, pmids: list[int]) -> None:
    done = embedded_pmids()
    missing = [pmid for pmid in pmids if pmid not in done]
    print(f"EMBED_HAVE {len(done)} MISSING {len(missing)} spent {spent():.2f}", flush=True)
    if not missing:
        return
    conn = sqlite3.connect(DB)
    texts = {}
    for start in range(0, len(missing), 400):
        batch = missing[start : start + 400]
        marks = ",".join("?" for _ in batch)
        for pmid, title, abstract in conn.execute(
            f"SELECT pmid, title, abstract FROM abstracts WHERE pmid IN ({marks})",
            batch,
        ):
            texts[int(pmid)] = cap_text(f"{title or ''} {abstract or ''}".strip())
    pending = [(pmid, texts[pmid]) for pmid in missing if texts.get(pmid)]
    blank = len(missing) - len(pending)
    if blank:
        print(f"EMBED_BLANK {blank}", flush=True)
    cache_dir().mkdir(parents=True, exist_ok=True)
    existing = chunk_paths()
    next_chunk = 1 + max((int(path.stem.split("_")[1]) for path in existing), default=-1)
    total = (len(pending) + BATCH - 1) // BATCH
    for batch_index, start in enumerate(range(0, len(pending), BATCH), start=1):
        if spent() >= CAP:
            raise SystemExit(f"CAP embed at {spent():.2f}, batch {batch_index} of {total}")
        batch = pending[start : start + BATCH]
        payload = embed_request(api_key, [text for _, text in batch])
        rows = sorted(payload["data"], key=lambda item: item["index"])
        if len(rows) != len(batch):
            raise SystemExit("embedding response length mismatch")
        vectors = np.asarray([row["embedding"] for row in rows], dtype=np.float32)
        pmid_arr = np.asarray([pmid for pmid, _ in batch], dtype=np.int64)
        out = cache_dir() / f"chunk_{next_chunk:05d}.npz"
        tmp = cache_dir() / f"partial_{next_chunk:05d}.npz"
        np.savez(tmp, pmids=pmid_arr, vectors=vectors)
        tmp.replace(out)
        next_chunk += 1
        usage = payload.get("usage") or {}
        tokens = int(usage.get("total_tokens") or usage.get("prompt_tokens") or 0)
        MODEL_SPEND["embed_usd"] += tokens * EMBED_USD
        save_spend()
        if batch_index == 1 or batch_index % 50 == 0 or batch_index == total:
            print(f"EMBED_BATCH {batch_index} of {total} ${spent():.2f}", flush=True)


def embed_queries(api_key: str, rows: list[dict]) -> dict[str, np.ndarray]:
    path = EXAM / "query_vectors.npz"
    wanted = {row["neighborhood_id"]: measurement(row) for row in rows}
    saved = {}
    if path.exists():
        data = np.load(path, allow_pickle=True)
        saved = {key: data[key] for key in data.files}
    missing_ids = [key for key in wanted if key not in saved]
    print(f"QUERY_MISSING {len(missing_ids)}", flush=True)
    for start in range(0, len(missing_ids), 96):
        if spent() >= CAP:
            raise SystemExit(f"CAP query embed at {spent():.2f}")
        batch_ids = missing_ids[start : start + 96]
        payload = embed_request(api_key, [wanted[key] for key in batch_ids])
        ordered = sorted(payload["data"], key=lambda item: item["index"])
        for key, item in zip(batch_ids, ordered):
            vector = np.asarray(item["embedding"], dtype=np.float32)
            norm = float(np.linalg.norm(vector))
            if norm:
                vector /= norm
            saved[key] = vector
        usage = payload.get("usage") or {}
        tokens = int(usage.get("total_tokens") or usage.get("prompt_tokens") or 0)
        MODEL_SPEND["embed_usd"] += tokens * EMBED_USD
        save_spend()
    np.savez(path, **saved)
    return saved


def load_needed_matrix(needed: set[int]) -> tuple[np.ndarray, dict[int, int]]:
    pmids = []
    vectors = []
    for path in chunk_paths():
        data = np.load(path)
        arr = data["pmids"].astype(np.int64)
        keep = np.array([int(pmid) in needed for pmid in arr])
        if not keep.any():
            continue
        pmids.append(arr[keep])
        vectors.append(data["vectors"].astype(np.float32)[keep])
    pmid_arr = np.concatenate(pmids)
    matrix = np.vstack(vectors)
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    matrix /= norms
    index = {}
    for i, pmid in enumerate(pmid_arr):
        index[int(pmid)] = i
    return matrix, index


def done_ids(path: Path) -> set[str]:
    found = set()
    if not path.exists():
        return found
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            found.add(json.loads(line)["id"])
    return found


def chat(api_key: str, user: str) -> dict:
    body = {
        "model": MODEL,
        "temperature": 0,
        "max_completion_tokens": 500,
        "reasoning_effort": "none",
        "messages": [
            {"role": "system", "content": SYSTEM},
            {"role": "user", "content": user},
        ],
    }
    last = ""
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
    else:
        raise SystemExit(f"API failed after retries: {last}")
    usage = payload.get("usage") or {}
    prompt = int(usage.get("prompt_tokens") or 0)
    completion = int(usage.get("completion_tokens") or 0)
    text = payload["choices"][0]["message"].get("content") or ""
    return {
        "text": text,
        "prompt_tokens": prompt,
        "completion_tokens": completion,
        "cost_usd": prompt * INPUT_USD_PER_TOKEN + completion * OUTPUT_USD_PER_TOKEN,
    }


def record_answer(row: dict, user: str, api_key: str, handle) -> None:
    if spent() >= CAP:
        raise SystemExit(f"CAP before {row['neighborhood_id']} at {spent():.2f}")
    result = chat(api_key, user)
    gene, confidence = parse_gene(result["text"])
    allowed = {symbol.upper() for symbol in symbols_of(row)}
    predicted = gene.strip()
    in_list = predicted.upper() in allowed if predicted else False
    correct = bool(predicted) and in_list and predicted.upper() == row["label_symbol"].upper()
    MODEL_SPEND["gpt_usd"] += result["cost_usd"]
    save_spend()
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
        "cost_usd": result["cost_usd"],
        "prompt_tokens": result["prompt_tokens"],
        "completion_tokens": result["completion_tokens"],
    }) + "\n")
    handle.flush()


def run_nopaper(api_key: str, rows: list[dict]) -> None:
    finished = done_ids(NO_PAPER)
    with NO_PAPER.open("a", encoding="utf-8") as handle:
        todo = [row for row in rows if row["neighborhood_id"] not in finished]
        print(f"NOPAPER todo {len(todo)} spent {spent():.2f}", flush=True)
        for index, row in enumerate(todo, start=1):
            user = f"Identify the causal gene.\nGWAS phenotype: {row['trait_name']}\nGenes in locus: {row['candidate_genes']}\n"
            record_answer(row, user, api_key, handle)
            if index % 50 == 0 or index == len(todo):
                print(f"NOPAPER {index} of {len(todo)} ${spent():.2f}", flush=True)


def paper_user(row: dict, abstracts: list[dict]) -> str:
    lines = [
        "Identify the causal gene.",
        f"GWAS phenotype: {row['trait_name']}",
        f"Measurement: {measurement(row)}",
        f"Genes in locus: {row['candidate_genes']}",
        "",
        "Abstracts:",
    ]
    for number, doc in enumerate(abstracts, start=1):
        abstract = doc["abstract"]
        if len(abstract) > ABSTRACT_CAP:
            abstract = abstract[:ABSTRACT_CAP].rsplit(" ", 1)[0]
        year = doc["year"] if doc["year"] else "unknown year"
        lines.append(f"{number}. ({year}) {doc['title']}. {abstract}")
    return "\n".join(lines) + "\n"


def run_papers(api_key: str, rows: list[dict], pools: dict[str, dict], queries: dict[str, np.ndarray]) -> None:
    needed = {pmid for pool in pools.values() for pmid in pool["pmids"]}
    print(f"LOAD_MATRIX pmids {len(needed)}", flush=True)
    matrix, index = load_needed_matrix(needed)
    conn = sqlite3.connect(DB)
    finished = done_ids(PAPERS)
    todo = [row for row in rows if row["neighborhood_id"] not in finished]
    print(f"PAPERS todo {len(todo)} spent {spent():.2f}", flush=True)
    with PAPERS.open("a", encoding="utf-8") as handle:
        for number, row in enumerate(todo, start=1):
            pool = pools[row["neighborhood_id"]]["pmids"]
            query = queries[row["neighborhood_id"]]
            ranked = []
            usable = [pmid for pmid in pool if pmid in index]
            if usable:
                scores = matrix[[index[pmid] for pmid in usable]] @ query
                order = np.argsort(-scores)[:TOP_N]
                ranked = [usable[int(i)] for i in order]
            docs = []
            if ranked:
                marks = ",".join("?" for _ in ranked)
                found = {
                    int(pmid): {"year": year, "title": title or "", "abstract": abstract or ""}
                    for pmid, year, title, abstract in conn.execute(
                        f"SELECT pmid, year, title, abstract FROM abstracts WHERE pmid IN ({marks})",
                        ranked,
                    )
                }
                for pmid in ranked:
                    if pmid in found:
                        docs.append(found[pmid])
            record_answer(row, paper_user(row, docs), api_key, handle)
            if number % 25 == 0 or number == len(todo):
                print(f"PAPERS {number} of {len(todo)} ${spent():.2f}", flush=True)


def load_results(path: Path) -> dict[tuple, dict]:
    found = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        found[(row["chrom"], row["pos"], row["rsid"], row["label"], row["trait"])] = row
    return found


def correct_of(record: dict | None) -> bool:
    return bool(record and record.get("correct"))


def choose(predicted: str, confidence: float, caldera: str) -> str:
    if predicted and caldera and predicted.upper() == caldera.upper():
        return predicted
    if predicted and confidence >= LINE:
        return predicted
    return caldera or ""


def paired(old: list[bool], new: list[bool]) -> dict:
    fixed = sum(1 for a, b in zip(old, new) if (not a) and b)
    broken = sum(1 for a, b in zip(old, new) if a and (not b))
    net = [b - a for a, b in zip(old, new)]
    return {"fixed": fixed, "broken": broken, "net": fixed - broken, "nets": net}


def bootstrap(flags: list[bool], labels: list[str], nets: list[int] | None = None) -> dict:
    rng = np.random.default_rng(SEED)
    arr = np.asarray(flags, dtype=np.float64)
    n = len(arr)
    draws = rng.integers(0, n, size=(BOOT, n))
    rates = arr[draws].mean(axis=1)
    out = {
        "n": n,
        "correct": int(arr.sum()),
        "percent": round(100 * float(arr.mean()), 2) if n else None,
        "question_interval": [round(100 * float(np.percentile(rates, 2.5)), 2), round(100 * float(np.percentile(rates, 97.5)), 2)] if n else None,
    }
    groups: dict[str, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        groups[label].append(index)
    names = list(groups)
    gene_rates = []
    gene_nets = []
    for _ in range(BOOT):
        picked = rng.choice(len(names), size=len(names), replace=True)
        idxs = [index for name_i in picked for index in groups[names[int(name_i)]]]
        gene_rates.append(float(arr[idxs].mean()) if idxs else 0.0)
        if nets is not None:
            gene_nets.append(int(sum(nets[index] for index in idxs)))
    out["gene_interval"] = [round(100 * float(np.percentile(gene_rates, 2.5)), 2), round(100 * float(np.percentile(gene_rates, 97.5)), 2)] if n else None
    if nets is not None:
        question_nets = np.asarray(nets, dtype=np.int64)[draws].sum(axis=1)
        out["net_question_interval"] = [int(np.percentile(question_nets, 2.5)), int(np.percentile(question_nets, 97.5))]
        out["net_gene_interval"] = [int(np.percentile(gene_nets, 2.5)), int(np.percentile(gene_nets, 97.5))]
    return out


def mcnemar(fixed: int, broken: int) -> float | None:
    from scipy.stats import binomtest

    total = fixed + broken
    if total == 0:
        return 1.0
    return float(binomtest(fixed, total, 0.5, alternative="two-sided").pvalue)


def score() -> None:
    nopaper = load_results(NO_PAPER)
    paper = load_results(PAPERS)
    caldera: dict[str, str] = {}
    with (EXAM / "caldera_genes.tsv").open(encoding="utf-8", newline="") as handle:
        best: dict[str, tuple[float, str]] = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            score_value = float(row["caldera"])
            if row["locus"] not in best or score_value > best[row["locus"]][0]:
                best[row["locus"]] = (score_value, row["gene"])
        caldera = {locus: gene for locus, (_, gene) in best.items()}

    def pack(rows: list[dict], mode: str) -> dict:
        old, new, labels, nearest = [], [], [], []
        for row in rows:
            key = (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])
            base = nopaper.get(key)
            with_papers = paper.get(key)
            nearest.append((row.get("nearest_gene") or "").upper() == row["label_symbol"].upper())
            labels.append(row["label_symbol"])
            if mode == "gpt":
                old.append(correct_of(base))
                new.append(correct_of(with_papers))
            else:
                cal = caldera.get(row["neighborhood_id"], "")
                base_gene = (base or {}).get("predicted") or ""
                base_conf = float((base or {}).get("confidence") or 0)
                paper_gene = (with_papers or {}).get("predicted") or ""
                paper_conf = float((with_papers or {}).get("confidence") or 0)
                agree = bool(base_gene and cal and base_gene.upper() == cal.upper())
                use_gene = base_gene if agree else paper_gene
                use_conf = base_conf if agree else paper_conf
                old.append(choose(base_gene, base_conf, cal).upper() == row["label_symbol"].upper() and bool(choose(base_gene, base_conf, cal)))
                new.append(choose(use_gene, use_conf, cal).upper() == row["label_symbol"].upper() and bool(choose(use_gene, use_conf, cal)))
        pair = paired([bool(item) for item in old], [bool(item) for item in new])
        summary = {
            "no_papers": bootstrap([bool(item) for item in old], labels),
            "with_papers": bootstrap([bool(item) for item in new], labels),
            "nearest": bootstrap(nearest, labels),
            "fixed": pair["fixed"],
            "broken": pair["broken"],
            "net": pair["net"],
            "net_question_interval": bootstrap([bool(item) for item in new], labels, pair["nets"])["net_question_interval"],
            "net_gene_interval": bootstrap([bool(item) for item in new], labels, pair["nets"])["net_gene_interval"],
            "mcnemar_p": mcnemar(pair["fixed"], pair["broken"]),
        }
        return summary

    main = list(csv.DictReader((EXAM / "exam_main.tsv").open(encoding="utf-8"), delimiter="\t"))
    give = list(csv.DictReader((EXAM / "exam_giveaway.tsv").open(encoding="utf-8"), delimiter="\t"))
    full = list(csv.DictReader((EXAM / "exam_a.tsv").open(encoding="utf-8"), delimiter="\t"))
    exam_b = [row for row in full if row["in_exam_b"] == "1"]

    def hypo(rows):
        return [row for row in rows if "hypothyroid" in f"{row['trait_name']} {row['mapped_trait']}".lower()]

    def reported(rows, flag):
        return [row for row in rows if (row["previously_reported"] == "1") is flag]

    slices = {
        "main": pack(main, "gpt"),
        "main_not_previously_reported": pack(reported(main, False), "gpt"),
        "main_lipid": pack([row for row in main if row["trait_group"] == "blood_lipid"], "gpt"),
        "main_protein": pack([row for row in main if row["trait_group"] == "blood_protein"], "gpt"),
        "main_other": pack([row for row in main if row["trait_group"] not in {"blood_lipid", "blood_protein"}], "gpt"),
        "giveaway": pack(give, "gpt"),
        "full": pack(full, "gpt"),
        "hypothyroid": pack(hypo(main), "gpt"),
        "exam_b_gpt": pack(exam_b, "gpt"),
        "exam_b_rule": pack(exam_b, "rule"),
    }
    # CALDERA alone and rule pieces on Exam B
    cal_flags = []
    labels = []
    for row in exam_b:
        cal = caldera.get(row["neighborhood_id"], "")
        cal_flags.append(bool(cal) and cal.upper() == row["label_symbol"].upper())
        labels.append(row["label_symbol"])
    slices["exam_b_caldera"] = bootstrap(cal_flags, labels)
    uk_gpt = {"fixed": 64, "broken": 26}
    uk_rule = {"fixed": 62, "broken": 42}
    pooled = {
        "gpt": {
            "uk": uk_gpt,
            "main": {"fixed": slices["main"]["fixed"], "broken": slices["main"]["broken"]},
            "pooled_fixed": uk_gpt["fixed"] + slices["main"]["fixed"],
            "pooled_broken": uk_gpt["broken"] + slices["main"]["broken"],
        },
        "rule": {
            "uk": uk_rule,
            "exam_b": {"fixed": slices["exam_b_rule"]["fixed"], "broken": slices["exam_b_rule"]["broken"]},
            "pooled_fixed": uk_rule["fixed"] + slices["exam_b_rule"]["fixed"],
            "pooled_broken": uk_rule["broken"] + slices["exam_b_rule"]["broken"],
        },
    }
    pooled["gpt"]["mcnemar_p"] = mcnemar(pooled["gpt"]["pooled_fixed"], pooled["gpt"]["pooled_broken"])
    pooled["rule"]["mcnemar_p"] = mcnemar(pooled["rule"]["pooled_fixed"], pooled["rule"]["pooled_broken"])
    hallucinations = 0
    for path in (NO_PAPER, PAPERS):
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.strip() and json.loads(line).get("hallucination"):
                hallucinations += 1
    payload = {"spend": MODEL_SPEND, "slices": slices, "pooled": pooled, "hallucinations": hallucinations}
    SUMMARY.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_context(payload)
    print("SCORED", flush=True)


def write_context(payload: dict) -> None:
    main = payload["slices"]["main"]
    rule = payload["slices"]["exam_b_rule"]
    text = f"""## Final exam results

Written after the single run on 9 October 2026. The questions and the rules were locked before any call.

The main list has {main['no_papers']['n']} questions, after protein measurements were treated as one trait and after questions that name the measured gene were set aside. With no papers, the language model named the labeled gene on {main['no_papers']['correct']} of them ({main['no_papers']['percent']} percent). With the 25 papers added to every question, it named the labeled gene on {main['with_papers']['correct']} ({main['with_papers']['percent']} percent). Papers fixed {main['fixed']} questions and broke {main['broken']}. The net change is {main['net']} questions. Redrawing whole questions, a plausible range for that net change is {main['net_question_interval'][0]} to {main['net_question_interval'][1]}. Redrawing whole labeled genes, the range is {main['net_gene_interval'][0]} to {main['net_gene_interval'][1]}.

The combined rule on the 38 questions where the statistical model has both inputs was right on {rule['no_papers']['correct']} with no papers and {rule['with_papers']['correct']} with papers. That list is too small to confirm or reject the setup by itself. Net change {rule['net']} (fixed {rule['fixed']}, broke {rule['broken']}).

One kept statistical-model input still gives TREM2, the labeled gene, a protein-changing probability of 0.41. The rule that drops a set only above 0.5 was fixed before the run, so that question stayed in.

The full numbers, the secondary slices, and the pooled counts are in `data/ukb/final_exam/results_summary.json`. Spend on this run: ${payload['spend']['embed_usd']:.2f} to embed abstracts and ${payload['spend']['gpt_usd']:.2f} for the model calls.
"""
    current = CONTEXT.read_text(encoding="utf-8")
    if "## Final exam results" in current:
        start = current.index("## Final exam results")
        rest = current[start + len("## Final exam results"):]
        nxt = rest.find("\n## ")
        end = len(current) if nxt < 0 else start + len("## Final exam results") + nxt
        current = current[:start] + text.strip() + "\n\n" + current[end:].lstrip("\n")
    else:
        marker = "## Files worth knowing"
        current = current.replace(marker, text.strip() + "\n\n" + marker, 1)
    CONTEXT.write_text(current, encoding="utf-8")


def main() -> None:
    load_spend()
    gate()
    rows = rows_to_call()
    api_key = load_key()
    build_pools(rows)
    pools = load_pools()
    needed = sorted({pmid for pool in pools.values() for pmid in pool["pmids"]})
    embed_missing(api_key, needed)
    queries = embed_queries(api_key, rows)
    run_nopaper(api_key, rows)
    run_papers(api_key, rows, pools, queries)
    score()
    print(f"DONE ${spent():.2f}", flush=True)


if __name__ == "__main__":
    main()

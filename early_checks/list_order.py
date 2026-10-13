"""Positional-bias and fame-bias check on gpt-5.4.

Samples high-confidence Open Targets loci, asks the model to rank the
protein-coding genes within 500 kb, and repeats the prompt with the list
shuffled. Fame is ClinVar total submissions for the gene.

Stops at MAX_USD. Reads OPENAI_API_KEY from the repo .env and does not print it.
"""

from __future__ import annotations

import csv
import gzip
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
GOLD_PATH = Path(os.environ['EARLY_CHECK_DATA']) / 'gold.tsv'
CLINVAR_PATH = Path(os.environ['EARLY_CHECK_DATA']) / 'gene_specific_summary.txt'
GENCODE_PATH = Path(os.environ['EARLY_CHECK_DATA']) / 'gencode.v47.annotation.gtf.gz'
CACHE = Path(__file__).resolve().parent / "cache"
LOCI_PATH = CACHE / "loci.json"
RESULTS_PATH = Path(__file__).resolve().parent / "results.jsonl"
SUMMARY_PATH = Path(__file__).resolve().parent / "list_order_summary.json"

MODEL = "gpt-5.4"
INPUT_USD_PER_TOKEN = 2.50 / 1_000_000
OUTPUT_USD_PER_TOKEN = 15.00 / 1_000_000
MAX_USD = 8.0
SEED = 20261006
WINDOW = 500_000
MIN_GENES = 4
PER_ARM = 40
N_SHUFFLES = 4

SYSTEM = (
    "You are an expert in biology and genetics. "
    "Your task is to identify the likely causal gene within a locus for a GWAS phenotype "
    "based on literature evidence. Rank every gene in the list from most likely causal to least likely. "
    "The top gene must be one of the given genes. "
    "Return only JSON with keys causal_gene (string), confidence (number from 0 to 1), "
    "and ranking (array of gene symbols, best first, each given gene exactly once). "
    "No other text."
)


def load_api_key() -> str:
    for line in ENV_PATH.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == "OPENAI_API_KEY":
            return value.strip().strip('"').strip("'")
    raise SystemExit("OPENAI_API_KEY not found")


def http_json(url: str, payload: dict | None = None, timeout: int = 60) -> object:
    data = None if payload is None else json.dumps(payload).encode("utf-8")
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(5):
        req = urllib.request.Request(url, data=data, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            if exc.code in {429, 500, 503} and attempt < 4:
                time.sleep(1.5 * (attempt + 1))
                continue
            detail = exc.read().decode("utf-8", errors="replace")[:400]
            raise SystemExit(f"HTTP {exc.code} {url} {detail}") from exc
    raise SystemExit(f"failed {url}")


def load_clinvar() -> dict[str, int]:
    counts: dict[str, int] = {}
    for line in CLINVAR_PATH.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        try:
            counts[parts[0]] = int(parts[2])
        except ValueError:
            continue
    return counts


def load_loci() -> list[dict]:
    rows = list(csv.DictReader(GOLD_PATH.open(encoding="utf-8"), delimiter="\t"))
    best: dict[tuple, dict] = {}
    for row in rows:
        if row["gold_standard_info.highest_confidence"] != "High":
            continue
        rsid = row["sentinel_variant.rsid"].strip()
        chrom = row["sentinel_variant.locus_GRCh38.chromosome"].strip()
        pos = row["sentinel_variant.locus_GRCh38.position"].strip()
        gene = row["gold_standard_info.gene_id"].split(".")[0].strip()
        trait = row["trait_info.standard_trait_name"].strip() or row["trait_info.reported_trait_name"].strip()
        if not (rsid and chrom and pos and gene and trait):
            continue
        key = (rsid, trait, chrom, pos)
        conf = row["gold_standard_info.evidence.confidence"]
        rank = {"High": 3, "Medium": 2, "Low": 1}.get(conf.split("|")[0], 0)
        current = best.get(key)
        if current is None or rank > current["conf_rank"]:
            best[key] = {
                "rsid": rsid,
                "trait": trait,
                "chrom": chrom,
                "pos": int(pos),
                "gene_id": gene,
                "conf_rank": rank,
            }
    return list(best.values())


def load_gencode() -> dict[str, list[tuple[int, int, str, str]]]:
    by_chrom: dict[str, list[tuple[int, int, str, str]]] = {}
    with gzip.open(GENCODE_PATH, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            cols = line.split("\t")
            if len(cols) < 9 or 'gene_type "protein_coding"' not in cols[8]:
                continue
            chrom = cols[0].removeprefix("chr")
            gid = re.search(r'gene_id "([^"]+)"', cols[8])
            name = re.search(r'gene_name "([^"]+)"', cols[8])
            if not gid or not name:
                continue
            by_chrom.setdefault(chrom, []).append(
                (int(cols[3]), int(cols[4]), gid.group(1).split(".")[0], name.group(1))
            )
    return by_chrom


def genes_in_window(gencode: dict[str, list[tuple[int, int, str, str]]], chrom: str, pos: int) -> list[dict]:
    start = max(1, pos - WINDOW)
    end = pos + WINDOW
    genes = []
    seen = set()
    for gene_start, gene_end, gid, symbol in gencode.get(chrom, []):
        if gene_end < start or gene_start > end or gid in seen or not symbol:
            continue
        seen.add(gid)
        genes.append({"gene_id": gid, "symbol": symbol})
    return genes


def prepare() -> list[dict]:
    CACHE.mkdir(parents=True, exist_ok=True)
    if LOCI_PATH.exists():
        return json.loads(LOCI_PATH.read_text(encoding="utf-8"))
    clinvar = load_clinvar()
    raw = load_loci()
    gencode = load_gencode()
    eligible = []
    for row in raw:
        genes = genes_in_window(gencode, row["chrom"], row["pos"])
        by_id = {gene["gene_id"]: gene["symbol"] for gene in genes}
        if row["gene_id"] not in by_id or len(genes) < MIN_GENES:
            continue
        symbol = by_id[row["gene_id"]]
        eligible.append(
            {
                "rsid": row["rsid"],
                "trait": row["trait"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "gene_id": row["gene_id"],
                "symbol": symbol,
                "clinvar": clinvar.get(symbol, 0),
                "genes": [gene["symbol"] for gene in genes],
                "clinvar_by_gene": {gene["symbol"]: clinvar.get(gene["symbol"], 0) for gene in genes},
            }
        )
    counts = sorted(row["clinvar"] for row in eligible)
    median = counts[len(counts) // 2] if counts else 0
    low = [row for row in eligible if row["clinvar"] <= median]
    high = [row for row in eligible if row["clinvar"] > median]
    rng = random.Random(SEED)
    rng.shuffle(low)
    rng.shuffle(high)
    n = min(PER_ARM, len(low), len(high))
    chosen = low[:n] + high[:n]
    for row in chosen:
        row["fame_arm"] = "obscure" if row["clinvar"] <= median else "famous"
        row["median_clinvar"] = median
    LOCI_PATH.write_text(json.dumps(chosen, indent=2), encoding="utf-8")
    print(
        f"eligible {len(eligible)} of {len(raw)} high-confidence loci "
        f"with at least {MIN_GENES} protein-coding genes; "
        f"median ClinVar submissions {median}; sampled {n} obscure and {n} famous",
        flush=True,
    )
    return chosen


def call_model(api_key: str, trait: str, genes: list[str]) -> dict:
    user = (
        "Identify the causal gene.\n"
        f"GWAS phenotype: {trait}\n"
        "Genes in locus: " + ", ".join(genes)
    )
    body = {
        "model": MODEL,
        "temperature": 0,
        "max_completion_tokens": 800,
        "reasoning_effort": "none",
        "messages": [
            {"role": "system", "content": SYSTEM},
            {"role": "user", "content": user},
        ],
    }
    req = urllib.request.Request(
        "https://api.openai.com/v1/chat/completions",
        data=json.dumps(body).encode("utf-8"),
        headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[:500]
        if exc.code == 400 and "temperature" in detail:
            body.pop("temperature", None)
            req = urllib.request.Request(
                "https://api.openai.com/v1/chat/completions",
                data=json.dumps(body).encode("utf-8"),
                headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=120) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
        else:
            raise SystemExit(f"API {exc.code}: {detail}") from exc
    usage = payload.get("usage") or {}
    text = payload["choices"][0]["message"].get("content") or ""
    prompt_tokens = usage.get("prompt_tokens") or 0
    completion_tokens = usage.get("completion_tokens") or 0
    return {
        "text": text,
        "prompt_tokens": prompt_tokens,
        "completion_tokens": completion_tokens,
        "cost_usd": prompt_tokens * INPUT_USD_PER_TOKEN + completion_tokens * OUTPUT_USD_PER_TOKEN,
    }


def parse_reply(text: str, genes: list[str]) -> dict:
    match = re.search(r"\{.*\}", text, re.S)
    parsed = {}
    if match:
        try:
            parsed = json.loads(match.group(0))
        except json.JSONDecodeError:
            parsed = {}
    by_norm = {gene.casefold(): gene for gene in genes}
    raw_top = str(parsed.get("causal_gene") or "")
    top = by_norm.get(raw_top.casefold(), "")
    ranking = []
    for item in parsed.get("ranking") or []:
        canon = by_norm.get(str(item).casefold())
        if canon and canon not in ranking:
            ranking.append(canon)
    try:
        confidence = float(parsed.get("confidence"))
    except (TypeError, ValueError):
        confidence = None
    return {"top": top, "ranking": ranking, "confidence": confidence}


def rank_of(ranking: list[str], symbol: str, n_genes: int) -> int | None:
    if symbol in ranking:
        return ranking.index(symbol) + 1
    return None


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(xs) != len(ys):
        return None

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = avg
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = sum((a - mx) ** 2 for a in rx) ** 0.5
    deny = sum((b - my) ** 2 for b in ry) ** 0.5
    if denx == 0 or deny == 0:
        return None
    return num / (denx * deny)


def summarize(loci: list[dict], trials: list[dict]) -> dict:
    by_locus: dict[str, list[dict]] = {}
    for trial in trials:
        by_locus.setdefault(trial["locus_id"], []).append(trial)
    stable = 0
    pair_same = 0
    pair_n = 0
    correct = 0
    scored = 0
    pos_x: list[float] = []
    pos_y: list[float] = []
    first_third_hit = [0, 0]
    last_third_hit = [0, 0]
    fame_wrong_higher = 0
    fame_wrong_n = 0
    chance_higher = []
    within_list_rho = []
    for locus in loci:
        rows = by_locus.get(locus["id"], [])
        tops = [row["top"] for row in rows if row["top"]]
        if len(set(tops)) == 1 and len(tops) == N_SHUFFLES:
            stable += 1
        for i in range(len(tops)):
            for j in range(i + 1, len(tops)):
                pair_n += 1
                pair_same += int(tops[i] == tops[j])
        for row in rows:
            if not row["top"]:
                continue
            scored += 1
            correct += int(row["top"] == locus["symbol"])
            if row["output_rank"] is not None and locus["n_genes"] > 1:
                pos_x.append(row["input_index"] / (locus["n_genes"] - 1))
                pos_y.append(row["output_rank"])
            frac = row["input_index"] / max(locus["n_genes"] - 1, 1)
            if frac <= 1 / 3:
                first_third_hit[1] += 1
                first_third_hit[0] += int(row["top"] == locus["symbol"])
            elif frac >= 2 / 3:
                last_third_hit[1] += 1
                last_third_hit[0] += int(row["top"] == locus["symbol"])
            if row["top"] != locus["symbol"]:
                fame_wrong_n += 1
                fame_wrong_higher += int(row["top_clinvar"] > locus["clinvar"])
                others = [g for g in locus["genes"] if g != locus["symbol"]]
                if others:
                    chance_higher.append(
                        sum(1 for g in others if locus["clinvar_by_gene"].get(g, 0) > locus["clinvar"])
                        / len(others)
                    )
            if row["ranking"]:
                xs = [float(locus["clinvar_by_gene"].get(g, 0)) for g in row["ranking"]]
                ys = list(range(1, len(row["ranking"]) + 1))
                rho = spearman(xs, ys)
                if rho is not None:
                    within_list_rho.append(rho)
    def rate(pair: list[int]) -> float | None:
        return pair[0] / pair[1] if pair[1] else None
    return {
        "n_loci": len(loci),
        "n_trials": len(trials),
        "accuracy": correct / scored if scored else None,
        "loci_same_pick_every_shuffle": stable,
        "loci_same_pick_rate": stable / len(loci) if loci else None,
        "shuffle_pair_agreement": pair_same / pair_n if pair_n else None,
        "spearman_input_position_vs_output_rank": spearman(pos_x, pos_y),
        "hit_rate_causal_gene_in_first_third": rate(first_third_hit),
        "hit_rate_causal_gene_in_last_third": rate(last_third_hit),
        "first_third_n": first_third_hit[1],
        "last_third_n": last_third_hit[1],
        "wrong_pick_more_famous_than_truth": (
            fame_wrong_higher / fame_wrong_n if fame_wrong_n else None
        ),
        "wrong_n": fame_wrong_n,
        "chance_other_gene_more_famous_than_truth": (
            sum(chance_higher) / len(chance_higher) if chance_higher else None
        ),
        "mean_within_list_spearman_clinvar_vs_rank": (
            sum(within_list_rho) / len(within_list_rho) if within_list_rho else None
        ),
    }


def main() -> None:
    api_key = load_api_key()
    chosen = prepare()
    done = set()
    if RESULTS_PATH.exists():
        for line in RESULTS_PATH.read_text(encoding="utf-8").splitlines():
            if line.strip():
                rec = json.loads(line)
                done.add((rec["locus_id"], rec["shuffle"]))
    spent = 0.0
    if RESULTS_PATH.exists():
        for line in RESULTS_PATH.read_text(encoding="utf-8").splitlines():
            if line.strip():
                spent += json.loads(line).get("cost_usd") or 0
    rng = random.Random(SEED)
    out = RESULTS_PATH.open("a", encoding="utf-8")
    n_calls = 0
    for locus in chosen:
        locus_id = f"{locus['rsid']}|{locus['trait']}"
        genes = list(locus["genes"])
        for shuffle in range(N_SHUFFLES):
            order = genes[:]
            rng.shuffle(order)
            if (locus_id, shuffle) in done:
                continue
            if spent >= MAX_USD:
                print(f"STOP at ${spent:.4f}", flush=True)
                out.close()
                return
            reply = call_model(api_key, locus["trait"], order)
            spent += reply["cost_usd"]
            n_calls += 1
            parsed = parse_reply(reply["text"], order)
            record = {
                "locus_id": locus_id,
                "shuffle": shuffle,
                "fame_arm": locus["fame_arm"],
                "symbol": locus["symbol"],
                "clinvar": locus["clinvar"],
                "n_genes": len(order),
                "input_index": order.index(locus["symbol"]),
                "order": order,
                "top": parsed["top"],
                "top_clinvar": locus["clinvar_by_gene"].get(parsed["top"], 0) if parsed["top"] else None,
                "output_rank": rank_of(parsed["ranking"], locus["symbol"], len(order)),
                "ranking": parsed["ranking"],
                "confidence": parsed["confidence"],
                "cost_usd": reply["cost_usd"],
                "prompt_tokens": reply["prompt_tokens"],
                "completion_tokens": reply["completion_tokens"],
                "response": reply["text"],
            }
            out.write(json.dumps(record) + "\n")
            out.flush()
            print(
                f"{n_calls:03d} {locus['fame_arm']} {locus['rsid']} "
                f"top={parsed['top'] or 'PARSE?'} truth={locus['symbol']} "
                f"${spent:.4f}",
                flush=True,
            )
    out.close()
    trials = [json.loads(line) for line in RESULTS_PATH.read_text(encoding="utf-8").splitlines() if line.strip()]
    loci = []
    for locus in chosen:
        loci.append(
            {
                "id": f"{locus['rsid']}|{locus['trait']}",
                "symbol": locus["symbol"],
                "clinvar": locus["clinvar"],
                "genes": locus["genes"],
                "clinvar_by_gene": locus["clinvar_by_gene"],
                "n_genes": len(locus["genes"]),
                "fame_arm": locus["fame_arm"],
            }
        )
    by_arm = {}
    for arm in ("obscure", "famous"):
        arm_loci = [row for row in loci if row["fame_arm"] == arm]
        ids = {row["id"] for row in arm_loci}
        by_arm[arm] = summarize(arm_loci, [row for row in trials if row["locus_id"] in ids])
    summary = {
        "model": MODEL,
        "seed": SEED,
        "shuffles": N_SHUFFLES,
        "spent_usd": round(spent, 6),
        "median_clinvar": chosen[0]["median_clinvar"] if chosen else None,
        "overall": summarize(loci, trials),
        "by_arm": by_arm,
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

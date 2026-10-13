"""Build CALDERA inputs for Exam B from the noncoding lead's credible set only.

No scores against labels. A neighborhood whose lead is not in any SuSiE-inf
set is kept on the exam and gets no input file, so CALDERA has no pick there.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAM = ROOT / "data" / "ukb" / "final_exam" / "exam_a.tsv"
from pipeline.paths import trait_scores

PIPS = trait_scores() / "work" / "pips"
EXTRA = ROOT / "data" / "ukb" / "final_exam" / "ot_extra_susie.tsv"
SCORES = trait_scores() / "work" / "scores"
POPS_GZ = ROOT / "data" / "ukb" / "pops" / "PoPS_FullResults.txt.gz"
OUT = ROOT / "data" / "ukb" / "final_exam" / "caldera_inputs"
CACHE = ROOT / "data" / "ukb" / "final_exam" / "variant_consequences.jsonl"
API = "https://api.platform.opentargets.org/api/v4/graphql"
CODING = {
    "missense_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "incomplete_terminal_codon_variant",
}


def post(query: str, variables: dict | None = None) -> dict:
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    request = urllib.request.Request(
        API, data=body, headers={"Content-Type": "application/json", "User-Agent": "exam-consequences"}
    )
    last = ""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode())
            if payload.get("errors"):
                raise RuntimeError(json.dumps(payload["errors"])[:500])
            return payload["data"]
        except urllib.error.HTTPError as exc:
            last = exc.read().decode("utf-8", errors="replace")[:400]
            if exc.code not in (429, 500, 502, 503):
                raise RuntimeError(last) from exc
            time.sleep(2 + attempt * 3)
        except urllib.error.URLError:
            time.sleep(2 + attempt * 3)
    raise RuntimeError(last)


def load_exam() -> list[dict]:
    with EXAM.open(encoding="utf-8", newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t") if row["in_exam_b"] == "1"]


def load_sets(locus_ids: set[str]) -> dict[str, list[dict]]:
    found = defaultdict(list)
    seen = set()
    sources = list(PIPS.glob("GCST*.tsv"))
    if EXTRA.exists():
        sources.append(EXTRA)
    for path in sources:
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                locus = row.get("study_locus_id") or ""
                if locus not in locus_ids:
                    continue
                variant_id = row.get("variant_id") or f"{row['chromosome']}_{row['position']}_{row['ref']}_{row['alt']}"
                key = (locus, variant_id)
                if key in seen:
                    continue
                seen.add(key)
                found[locus].append(
                    {
                        "variant_id": variant_id,
                        "chr": row["chromosome"],
                        "bp": int(row["position"]),
                        "ref": row["ref"].upper(),
                        "alt": row["alt"].upper(),
                        "pip": float(row.get("posterior_probability") or row.get("pip") or 0),
                    }
                )
    return found


def coding_gene(record: dict) -> tuple[str, str]:
    genes = set()
    for transcript in record.get("transcriptConsequences") or []:
        target = transcript.get("target") or {}
        gene_id = target.get("id") or ""
        labels = {item.get("label") for item in (transcript.get("variantConsequences") or [])}
        if gene_id.startswith("ENSG") and labels & CODING:
            genes.add(gene_id)
    if len(genes) == 1:
        return next(iter(genes)), ""
    if len(genes) > 1:
        return "", "multiple"
    return "", ""


def fetch_consequences(variant_ids: list[str]) -> dict[str, dict]:
    known = {}
    if CACHE.exists():
        for line in CACHE.read_text(encoding="utf-8").splitlines():
            if line.strip():
                record = json.loads(line)
                known[record["id"]] = record
    missing = [variant_id for variant_id in variant_ids if variant_id not in known]
    print(f"variants {len(variant_ids)} cached {len(variant_ids) - len(missing)} to_fetch {len(missing)}", flush=True)
    fields = """
      id
      transcriptConsequences {
        target { id }
        variantConsequences { label }
      }
    """
    with CACHE.open("a", encoding="utf-8") as handle:
        for start in range(0, len(missing), 10):
            batch = missing[start:start + 10]
            parts = []
            variables = {}
            for index, variant_id in enumerate(batch):
                key = f"v{index}"
                parts.append(f"{key}: variant(variantId: ${key}) {{ {fields} }}")
                variables[key] = variant_id
            var_defs = " ".join(f"${key}: String!" for key in variables)
            query = f"query({var_defs}) {{ {' '.join(parts)} }}"
            data = post(query, variables)
            for key, variant_id in zip(variables, batch):
                record = data.get(key) or {"id": variant_id, "transcriptConsequences": []}
                record["id"] = variant_id
                known[variant_id] = record
                handle.write(json.dumps(record) + "\n")
            handle.flush()
            print(f"fetched {min(start + 10, len(missing))}/{len(missing)}", flush=True)
    return known


def main() -> None:
    exam = load_exam()
    with_set = [row for row in exam if row["has_susie_noncoding_lead"] == "1" and row["susie_locus_id"]]
    no_set = [row for row in exam if row not in with_set]
    print(f"exam_b {len(exam)} with_set {len(with_set)} no_set {len(no_set)}", flush=True)
    sets = load_sets({row["susie_locus_id"] for row in with_set})
    variant_ids = sorted({item["variant_id"] for items in sets.values() for item in items})
    consequences = fetch_consequences(variant_ids)

    kept = []
    drop_pip = []
    drop_same = []
    label_any = 0
    multi = 0
    for row in with_set:
        annotated = []
        sums = defaultdict(float)
        label_pip = 0.0
        same = False
        for item in sets.get(row["susie_locus_id"], []):
            record = consequences.get(item["variant_id"])
            gene, note = coding_gene(record or {})
            if note == "multiple":
                multi += 1
            copy = dict(item)
            copy["c_gene"] = gene
            copy["locus"] = row["neighborhood_id"]
            annotated.append(copy)
            if gene:
                sums[gene] += item["pip"]
            if gene and gene == row["label_ensembl"] and item["pip"] > 0:
                label_pip += item["pip"]
            allele = (row.get("coding_allele") or "").upper()
            if item["chr"] == row["noncoding_chrom"] and item["bp"] == int(row["coding_pos"] or 0):
                if not allele or allele in {item["ref"], item["alt"]}:
                    same = True
        if label_pip > 0:
            label_any += 1
        pip_drop = any(total > 0.5 for total in sums.values())
        if pip_drop:
            drop_pip.append(row["neighborhood_id"])
        if same:
            drop_same.append(row["neighborhood_id"])
        if pip_drop or same:
            continue
        row["variants"] = annotated
        row["any_c_gene"] = any(item["c_gene"] for item in annotated)
        kept.append(row)

    cgene_rows = sum(1 for row in kept if row["any_c_gene"])
    report = {
        "exam_b": len(exam),
        "with_noncoding_set": len(with_set),
        "no_set_caldera_wrong": len(no_set),
        "dropped_coding_pip_over_0.5": len(drop_pip),
        "dropped_same_signal": len(drop_same),
        "dropped_either": len(drop_pip) + len(drop_same) - len(set(drop_pip) & set(drop_same)),
        "kept": len(kept),
        "kept_with_any_c_gene": cgene_rows,
        "labeled_gene_coding_any_pip": label_any,
        "variants_with_multiple_coding_genes": multi,
        "stop": cgene_rows > 5,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (ROOT / "data" / "ukb" / "final_exam" / "caldera_filter_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2), flush=True)
    if report["stop"]:
        print("STOP more than a handful of kept sets still have a protein-changing gene", flush=True)
        return

    manifest = []
    for row in kept:
        if row["score_source"] == "landstra_file":
            pops_key = row["score_accession"]
        else:
            pops_key = row["pops_trait"]
        manifest.append((pops_key, row["score_source"], row))
    by_pops = defaultdict(list)
    for pops_key, source, row in manifest:
        by_pops[(pops_key, source)].extend(row["variants"])

    for (pops_key, source), variants in by_pops.items():
        cs_path = OUT / f"{pops_key}_cs.tsv"
        with cs_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["locus", "chr", "bp", "pip", "c_gene"], delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for item in variants:
                writer.writerow(
                    {
                        "locus": item["locus"],
                        "chr": item["chr"],
                        "bp": item["bp"],
                        "pip": item["pip"],
                        "c_gene": item["c_gene"],
                    }
                )
        if source == "landstra_file":
            pops_path = OUT / f"{pops_key}_pops.tsv"
            pops_path.write_bytes((SCORES / f"{pops_key}.tsv").read_bytes())
        else:
            write_pops(pops_key, OUT / f"{pops_key}_pops.tsv")
        print(f"wrote {pops_key} variants {len(variants)}", flush=True)
    print(f"groups {len(by_pops)}", flush=True)


def write_pops(trait: str, path: Path) -> None:
    import gzip

    seen = set()
    with gzip.open(POPS_GZ, "rt", encoding="utf-8", newline="") as handle, path.open("w", encoding="utf-8", newline="") as out:
        reader = csv.DictReader(handle, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=["ENSGID", "PoPS_Score"], delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row["trait"] != trait or row["ensgid"] in seen:
                continue
            seen.add(row["ensgid"])
            writer.writerow({"ENSGID": row["ensgid"], "PoPS_Score": row["pops_score"]})
    print(f"pops {trait} genes {len(seen)}", flush=True)


if __name__ == "__main__":
    main()

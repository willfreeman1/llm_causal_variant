"""Build the two files CALDERA reads for the 1,348 UK Biobank neighborhoods.

CALDERA matches a gene score by gene id only, and the same gene has a different
score for each trait. Each trait therefore gets its own pair of files:

- {trait}_cs.tsv with columns locus, chr, bp, pip
- {trait}_pops.tsv with columns ENSGID, PoPS_Score

The source variant table repeats every variant once for each nearby gene, and
its coding column marks the labeled gene. Those repeated gene rows are
collapsed to one row per variant. The coding column is not copied.
"""

import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEEKS = ROOT / "data" / "ukb" / "zenodo_directory_may_2025" / "data" / "benchmark_datasets" / "weeks_et_al_step1.tsv"
VARIANTS = ROOT / "data" / "ukb" / "pops" / "PoPS_UKBB_validation_data.txt.gz"
POPS = ROOT / "data" / "ukb" / "pops" / "PoPS_FullResults.txt.gz"
OUT = ROOT / "data" / "ukb" / "caldera_inputs"
VARIANT_RE = re.compile(r"^chr([0-9]+):([0-9]+):[ACGT]+:[ACGT]+$")


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def main():
    weeks = {}
    with open_text(WEEKS) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            weeks[row["locus_id"]] = row["phenotype"]
    if len(weeks) != 1348:
        raise SystemExit(f"expected 1348 neighborhoods, found {len(weeks)}")

    # One row per variant. The source table is variant by gene.
    seen = {}
    pip_sums = defaultdict(float)
    trait_of = {}
    bad_variants = []
    conflicts = []
    with open_text(VARIANTS) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "coding" not in reader.fieldnames:
            raise SystemExit("expected a coding column so we can be sure it is dropped")
        for row in reader:
            locus = row["locus_id"]
            if locus not in weeks:
                continue
            trait = row["trait"]
            if weeks[locus] != trait:
                raise SystemExit(f"trait mismatch for locus {locus}: {weeks[locus]} vs {trait}")
            trait_of[locus] = trait
            match = VARIANT_RE.match(row["variant"])
            if not match:
                bad_variants.append(row["variant"])
                continue
            chrom, pos = match.group(1), int(match.group(2))
            pip = float(row["pip"])
            key = (locus, row["variant"])
            previous = seen.get(key)
            if previous is None:
                seen[key] = (trait, chrom, pos, pip)
                pip_sums[locus] += pip
            elif previous != (trait, chrom, pos, pip):
                conflicts.append((key, previous, (trait, chrom, pos, pip)))

    if bad_variants:
        raise SystemExit(f"{len(bad_variants)} variants were not chr:pos:ref:alt, e.g. {bad_variants[:3]}")
    if conflicts:
        raise SystemExit(f"{len(conflicts)} variants had conflicting rows, e.g. {conflicts[:2]}")
    missing = sorted(set(weeks) - set(trait_of))
    if missing:
        raise SystemExit(f"{len(missing)} neighborhoods had no variants, e.g. {missing[:5]}")

    by_trait = defaultdict(list)
    for (locus, _variant), (trait, chrom, pos, pip) in seen.items():
        by_trait[trait].append((locus, chrom, pos, pip))

    scores = defaultdict(dict)
    with open_text(POPS) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            trait = row["trait"]
            if trait not in by_trait:
                continue
            gene = row["ensgid"].split(".")[0]
            score = float(row["pops_score"])
            prior = scores[trait].get(gene)
            if prior is None:
                scores[trait][gene] = score
            elif prior != score:
                raise SystemExit(f"{trait} {gene} has two different scores")

    missing_scores = sorted(set(by_trait) - set(scores))
    if missing_scores:
        raise SystemExit(f"no gene scores for {missing_scores}")

    OUT.mkdir(parents=True, exist_ok=True)
    low = []
    summary_traits = []
    for trait in sorted(by_trait):
        if not re.fullmatch(r"[A-Za-z0-9._-]+", trait):
            raise SystemExit(f"unsafe trait name {trait!r}")
        rows = sorted(by_trait[trait], key=lambda item: (item[0], item[1], item[2], -item[3]))
        cs_path = OUT / f"{trait}_cs.tsv"
        with open(cs_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["locus", "chr", "bp", "pip"])
            for locus, chrom, pos, pip in rows:
                writer.writerow([locus, chrom, pos, f"{pip:.8g}"])
        pops_path = OUT / f"{trait}_pops.tsv"
        with open(pops_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["ENSGID", "PoPS_Score"])
            for gene, score in sorted(scores[trait].items()):
                writer.writerow([gene, f"{score:.8g}"])
        loci = {row[0] for row in rows}
        under = sorted(locus for locus in loci if pip_sums[locus] < 0.95)
        low.extend(under)
        summary_traits.append({
            "trait": trait,
            "loci": len(loci),
            "variants": len(rows),
            "genes_with_scores": len(scores[trait]),
            "loci_with_pip_sum_below_0.95": len(under),
        })

    summary = {
        "neighborhoods": len(trait_of),
        "traits": len(by_trait),
        "variants": len(seen),
        "loci_with_pip_sum_below_0.95": len(low),
        "min_pip_sum": min(pip_sums.values()),
        "max_pip_sum": max(pip_sums.values()),
        "coding_column_written": False,
        "traits_detail": summary_traits,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in summary if k != "traits_detail"}, indent=2))


if __name__ == "__main__":
    main()

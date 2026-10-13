"""Regroup blood-protein questions and set aside giveaway questions.

A giveaway is a question whose trait is the level of one gene's product and
that gene is already in the candidate list. No model calls.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import build_final_exam_neighborhoods as base

ROOT = Path(__file__).resolve().parents[1]
EXAM = ROOT / "data" / "ukb" / "final_exam"
RAW = EXAM / "neighborhoods_all.tsv"
FULL = EXAM / "exam_a.tsv"
WINDOW = 500_000
CIRCULATING = re.compile(r"^Circulating\s+([A-Za-z0-9][A-Za-z0-9-]*)\s+levels\b", re.I)
MAPPED_PROTEIN = re.compile(r"\b(level|amount|concentration) of .+ in blood\b", re.I)
ASSAY = re.compile(r"\b(?:OID\d+|SeqId[:\s-]*\d+)\b", re.I)
ALIAS = {"NTPROBNP": "NPPB"}


def is_blood_protein(trait: str, mapped: str) -> bool:
    """One protein measured in blood, not a panel and not a lipid."""
    if CIRCULATING.match(trait or ""):
        return True
    return bool(MAPPED_PROTEIN.search(mapped or ""))


def candidate_symbols(row: dict) -> set[str]:
    found = set()
    for part in (row.get("candidate_genes") or "").split(","):
        symbol = part.strip().strip("{}").upper()
        if symbol:
            found.add(symbol)
    return found


def protein_symbol(trait: str) -> str:
    match = CIRCULATING.match(trait or "")
    if not match:
        return ""
    token = match.group(1).upper()
    return ALIAS.get(token, token)


def is_giveaway(row: dict) -> bool:
    """The measured gene is already sitting in the candidate list."""
    if not is_blood_protein(row["trait_name"], row["mapped_trait"]):
        return False
    candidates = candidate_symbols(row)
    symbol = protein_symbol(row["trait_name"])
    if symbol and symbol in candidates:
        return True
    blob = f"{row['trait_name']} {row['mapped_trait']}"
    for cand in candidates:
        if len(cand) < 3:
            continue
        if re.search(rf"\b{re.escape(cand)}\b", blob, re.I):
            return True
    if ASSAY.search(row["trait_name"] or "") and symbol and symbol in candidates:
        return True
    return False


def protein_genes_already_reported() -> set[str]:
    _, all_traits = base.load_studies()
    protein_keys = set()
    for key, meta in all_traits.items():
        if meta["journal"] >= base.REPORTED_BEFORE:
            continue
        if is_blood_protein(key[2], meta["mapped"]):
            protein_keys.add(key)
    print(f"older_protein_analyses {len(protein_keys)}", flush=True)
    found = set()
    import zipfile

    with zipfile.ZipFile(base.ASSOC_ZIP) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as raw:
            text = (line.decode("utf-8", errors="replace") for line in raw)
            for assoc in csv.DictReader(text, delimiter="\t"):
                day = base.parse_day(assoc.get("DATE") or "")
                if day is None or day >= base.REPORTED_BEFORE:
                    continue
                key = (
                    (assoc.get("PUBMEDID") or "").strip(),
                    (assoc.get("DATE") or "").strip(),
                    (assoc.get("DISEASE/TRAIT") or "").strip(),
                )
                if key not in protein_keys:
                    continue
                found |= base.gene_tokens(assoc.get("MAPPED_GENE") or "")
                found |= base.gene_tokens(assoc.get("REPORTED GENE(S)") or "")
    print(f"older_protein_genes {len(found)}", flush=True)
    return found


def collapse(rows: list[dict]) -> list[dict]:
    buckets = defaultdict(list)
    for index, row in enumerate(rows):
        buckets[(row["label_symbol"], row["trait_group"], row["noncoding_chrom"])].append(index)
    parent = list(range(len(rows)))

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for indexes in buckets.values():
        indexes.sort(key=lambda item: int(rows[item]["noncoding_pos"]))
        start = 0
        for end, index in enumerate(indexes):
            while int(rows[index]["noncoding_pos"]) - int(rows[indexes[start]]["noncoding_pos"]) > WINDOW:
                start += 1
            for other in indexes[start:end]:
                union(index, other)
    groups = defaultdict(list)
    for index in range(len(rows)):
        groups[find(index)].append(index)
    collapsed = []
    for number, indexes in enumerate(groups.values(), start=1):
        members = [rows[index] for index in indexes]

        def sort_key(row: dict) -> tuple:
            lead = row["has_trait_score"] == "1" and row["has_susie_noncoding_lead"] == "1"
            study = row["has_trait_score"] == "1" and row.get("study_has_susie") == "1"
            return (0 if lead else 1 if study else 2, row["study_accession"], int(row["noncoding_pos"]))

        chosen = dict(min(members, key=sort_key))
        chosen["n_collapsed"] = len(members)
        chosen["member_accessions"] = ",".join(sorted({row["study_accession"] for row in members}))
        collapsed.append(chosen)
    return collapsed


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with RAW.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        raw_fields = list(reader.fieldnames or [])
        rows = list(reader)
    known = protein_genes_already_reported()
    for row in rows:
        protein = is_blood_protein(row["trait_name"], row["mapped_trait"])
        row["giveaway"] = "1" if row["single_label"] == "1" and is_giveaway(row) else "0"
        if protein:
            row["trait_group"] = "blood_protein"
            if row["single_label"] == "1":
                row["previously_reported"] = "1" if row["label_symbol"].upper() in known else "0"
    if "giveaway" not in raw_fields:
        raw_fields.append("giveaway")
    write_table(RAW, rows, raw_fields)

    single = [row for row in rows if row["single_label"] == "1"]
    main_rows = collapse([row for row in single if row["giveaway"] != "1"])
    give_rows = collapse([row for row in single if row["giveaway"] == "1"])
    for number, row in enumerate(main_rows, start=1):
        row["neighborhood_id"] = f"M{number:05d}"
    for number, row in enumerate(give_rows, start=1):
        row["neighborhood_id"] = f"G{number:05d}"

    full = list(csv.DictReader(FULL.open(encoding="utf-8"), delimiter="\t"))
    exam_b = [row for row in full if row["in_exam_b"] == "1"]
    exam_b_keys = {
        (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])
        for row in exam_b
    }

    def contains_exam_b(row: dict) -> bool:
        return (
            row["noncoding_chrom"],
            row["noncoding_pos"],
            row["noncoding_rsid"],
            row["label_symbol"],
            row["trait_name"],
        ) in exam_b_keys

    # A collapsed question contains an Exam B row when its representative matches,
    # or when we stored only the representative. Exam B rows are themselves
    # representatives of the old groups. Match the representative fields.
    for row in main_rows:
        row["in_exam_b"] = "1" if contains_exam_b(row) else "0"
    for row in give_rows:
        row["in_exam_b"] = "1" if contains_exam_b(row) else "0"

    fields = raw_fields + ["n_collapsed", "member_accessions", "in_exam_b"]
    # raw rows do not have in_exam_b; collapsed copies do.
    for row in main_rows + give_rows:
        row.setdefault("in_exam_b", "0")
    write_table(EXAM / "exam_main.tsv", main_rows, fields)
    write_table(EXAM / "exam_giveaway.tsv", give_rows, fields)

    def part(items, pred):
        return sum(1 for row in items if pred(row))

    genes = Counter(row["label_symbol"] for row in main_rows)
    top = genes.most_common(10)
    summary = {
        "main": len(main_rows),
        "main_lipid": part(main_rows, lambda row: row["trait_group"] == "blood_lipid"),
        "main_protein": part(main_rows, lambda row: row["trait_group"] == "blood_protein"),
        "main_other": part(main_rows, lambda row: row["trait_group"] not in {"blood_lipid", "blood_protein"}),
        "main_not_previously_reported": part(main_rows, lambda row: row["previously_reported"] != "1"),
        "giveaway": len(give_rows),
        "exam_b_in_main": part(main_rows, lambda row: row["in_exam_b"] == "1"),
        "exam_b_in_giveaway": part(give_rows, lambda row: row["in_exam_b"] == "1"),
        "exam_b_expected": len(exam_b),
        "top_genes": top,
        "max_gene_share": (top[0][1] / len(main_rows)) if main_rows else 1,
    }
    print(summary, flush=True)
    missing = []
    matched = set()
    for row in exam_b:
        key = (row["noncoding_chrom"], row["noncoding_pos"], row["noncoding_rsid"], row["label_symbol"], row["trait_name"])
        hit = contains_exam_b(row) and any(
            (
                item["noncoding_chrom"],
                item["noncoding_pos"],
                item["noncoding_rsid"],
                item["label_symbol"],
                item["trait_name"],
            )
            == key
            for item in main_rows
        )
        if hit:
            matched.add(key)
        else:
            missing.append((row["neighborhood_id"], row["trait_name"][:60], row["label_symbol"]))
    print("exam_b_matched_as_main_representative", len(matched), "missing", missing[:20], flush=True)


if __name__ == "__main__":
    main()

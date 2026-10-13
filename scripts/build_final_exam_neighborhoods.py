"""Rebuild the post-31 August 2025 neighborhood list and collapse repeats.

The clean-neighborhood count uses the same rules as the 6,874 count:
join studies to associations on PubMed id, journal date, and trait name;
a variant is protein-changing when its catalog consequence includes one of
the coding terms below; a neighborhood is clean when exactly one such gene
sits within 500 kb of a lead that is not protein-changing.

No model calls. Open Targets is queried only for studies that already have
an older trait-match score and no local credible-set file.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.request
import zipfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "data" / "ukb" / "gwas_catalog_20261005"
STUDIES = CATALOG / "gwas-catalog-download-studies-v1.0.3.1.txt"
ASSOC_ZIP = CATALOG / "gwas-catalog-associations-full.zip"
GENES = (
    ROOT
    / "third_party"
    / "caldera_run"
    / "data"
    / "gene_locations_gencode_v47_grch38_all_protein_coding.tsv"
)
from pipeline.paths import trait_scores

PIPS = trait_scores() / "work" / "pips"
SCORES = trait_scores() / "work" / "scores"
OUT_DIR = ROOT / "data" / "ukb" / "final_exam"
API = "https://api.platform.opentargets.org/api/v4/graphql"

CUTOFF = datetime(2025, 8, 31)
REPORTED_BEFORE = datetime(2025, 9, 1)
WINDOW = 500_000
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
LIPID_FAMILY = [
    ("triacylglycerol", "triacylglycerol"),
    ("triglyceride", "triacylglycerol"),
    ("diacylglycerol", "diacylglycerol"),
    ("monoacylglycerol", "monoacylglycerol"),
    ("lysophosphatidylcholine", "lysophosphatidylcholine"),
    ("lysophosphatidyl", "lysophosphatidylcholine"),
    ("phosphatidylcholine", "phosphatidylcholine"),
    ("phosphatidylethanolamine", "phosphatidylethanolamine"),
    ("phosphatidylinositol", "phosphatidylinositol"),
    ("hexosylceramide", "hexosylceramide"),
    ("lactosylceramide", "lactosylceramide"),
    ("dihydrosylceramide", "dihydrosylceramide"),
    ("sphingomyelin", "sphingomyelin"),
    ("ceramide", "ceramide"),
    ("cholesteryl", "cholesteryl"),
    ("apolipoprotein", "apolipoprotein"),
    ("lipoprotein", "lipoprotein"),
    ("vldl", "vldl"),
    ("idl", "idl"),
    ("ldl", "ldl"),
    ("hdl", "hdl"),
    ("cholesterol", "cholesterol"),
    ("fatty acid", "fatty acid"),
    ("lipid", "lipid"),
]
TOKEN_SPLIT = re.compile(r"\s+-\s+|,|;|/")


def parse_day(text: str) -> datetime | None:
    text = (text or "").strip()
    if not text or text.upper() == "NA":
        return None
    try:
        return datetime.strptime(text[:10], "%Y-%m-%d")
    except ValueError:
        return None


def consequence_parts(context: str) -> set[str]:
    return {part.strip() for part in context.replace(",", ";").split(";") if part.strip()}


def gene_tokens(text: str) -> set[str]:
    found = set()
    for part in TOKEN_SPLIT.split(text or ""):
        token = part.strip().upper()
        if token and token not in {"NR", "NA", "INTERGENIC"} and " " not in token:
            found.add(token)
    return found


def lipid_family(trait: str, mapped: str) -> str:
    blob = f"{trait} {mapped}".lower()
    for needle, family in LIPID_FAMILY:
        if needle in {"ldl", "hdl", "vldl", "idl", "lipid"}:
            if re.search(rf"\b{needle}\b", blob):
                return family
        elif needle in blob:
            return family
    return ""


def pops_trait(trait: str, mapped: str) -> str:
    blob = f"{trait} {mapped}".lower()
    if "hypothyroid" in blob:
        return "Hypothyroidism"
    if "alzheimer" in blob and "family history" not in blob and "ltfh" not in blob:
        return "Alzheimer"
    if "bone mineral" in blob or "heel bone" in blob:
        return "eBMD"
    if "attention deficit" in blob or re.search(r"\badhd\b", blob):
        return "ADHD"
    if "cystatin" in blob and any(word in blob for word in ("filtr", "egfr", "clearance")):
        return "eGFRcys"
    if any(word in blob for word in ("glomerular", "creatinine clearance", "egfr", "filtration rate")):
        return "eGFR"
    if "fev" in blob and "fvc" in blob:
        return "FEV1FVC"
    return ""


def risk_allele(text: str) -> str:
    text = (text or "").strip()
    if "-" not in text:
        return ""
    allele = text.rsplit("-", 1)[-1].strip().upper()
    if re.fullmatch(r"[ACGT]+", allele):
        return allele
    return ""


def load_genes() -> dict[str, list[tuple[int, int, str, str]]]:
    by_chrom: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    with GENES.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_chrom[row["CHR"].strip()].append(
                (int(row["START"]), int(row["END"]), row["NAME"].strip().upper(), row["ENSGID"].strip())
            )
    return by_chrom


def genes_near(by_chrom, chrom: str, pos: int) -> list[tuple[int, str, str]]:
    found = []
    for start, end, symbol, ensg in by_chrom.get(chrom, []):
        if end < pos - WINDOW or start > pos + WINDOW:
            continue
        if start <= pos <= end:
            distance = 0
        else:
            distance = min(abs(pos - start), abs(pos - end))
        found.append((distance, symbol, ensg))
    return found


def load_studies() -> tuple[dict, dict]:
    """Post-cutoff analyses, plus every study's trait labels for the older-report check."""
    post: dict[tuple, dict] = {}
    all_traits: dict[tuple, dict] = {}
    dup_keys = 0
    with STUDIES.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            journal = parse_day(row["DATE"])
            if journal is None:
                continue
            key = (row["PUBMED ID"].strip(), row["DATE"].strip(), row["DISEASE/TRAIT"].strip())
            mapped = row["MAPPED_TRAIT"].strip()
            trait = row["DISEASE/TRAIT"].strip()
            info = {
                "accession": row["STUDY ACCESSION"].strip(),
                "sumstats": row["FULL SUMMARY STATISTICS"].strip().lower() == "yes",
                "mapped": mapped,
                "submission": (row.get("SUBMISSION DATE") or "").strip(),
                "added": (row.get("DATE ADDED TO CATALOG") or "").strip(),
                "journal": row["DATE"].strip(),
            }
            all_traits[key] = {
                "mapped": mapped,
                "family": lipid_family(trait, mapped),
                "journal": journal,
            }
            if journal <= CUTOFF:
                continue
            if key in post:
                dup_keys += 1
                post[key]["accessions"].append(info["accession"])
                post[key]["sumstats_any"] = post[key]["sumstats_any"] or info["sumstats"]
                for field in ("submission", "added"):
                    old = parse_day(post[key][field])
                    new = parse_day(info[field])
                    if new is not None and (old is None or new < old):
                        post[key][field] = info[field]
            else:
                info["accessions"] = [info["accession"]]
                info["sumstats_first"] = info["sumstats"]
                info["sumstats_any"] = info["sumstats"]
                post[key] = info
    print(f"post_cutoff_studies {len(post)} duplicate_keys {dup_keys}", flush=True)
    return post, all_traits


def load_associations(post: dict) -> tuple[dict, int, int]:
    by_key = defaultdict(list)
    unmatched = 0
    matched = 0
    with zipfile.ZipFile(ASSOC_ZIP) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as raw:
            text = (line.decode("utf-8", errors="replace") for line in raw)
            for row in csv.DictReader(text, delimiter="\t"):
                day = parse_day(row.get("DATE") or "")
                if day is None or day <= CUTOFF:
                    continue
                key = (
                    (row.get("PUBMEDID") or "").strip(),
                    (row.get("DATE") or "").strip(),
                    (row.get("DISEASE/TRAIT") or "").strip(),
                )
                if key not in post:
                    unmatched += 1
                    continue
                matched += 1
                chrom = (row.get("CHR_ID") or "").strip()
                pos_text = (row.get("CHR_POS") or "").strip()
                if not chrom or not pos_text.isdigit():
                    continue
                context = (row.get("CONTEXT") or "").strip()
                is_coding = bool(consequence_parts(context) & CODING)
                by_key[key].append(
                    {
                        "chrom": chrom,
                        "pos": int(pos_text),
                        "coding": is_coding,
                        "mapped": (row.get("MAPPED_GENE") or "").strip(),
                        "rsid": (row.get("SNPS") or "").strip(),
                        "allele": risk_allele(row.get("STRONGEST SNP-RISK ALLELE") or ""),
                        "ensg": (row.get("SNP_GENE_IDS") or "").strip(),
                    }
                )
    print(f"matched_assoc_rows {matched} unmatched_post_cutoff_rows {unmatched}", flush=True)
    return by_key, matched, unmatched


def build_clean(post: dict, by_key: dict) -> list[dict]:
    rows = []
    ambiguous = 0
    multi_symbol = 0
    for key, variants in by_key.items():
        coding = [item for item in variants if item["coding"] and item["mapped"]]
        other = [item for item in variants if not item["coding"]]
        meta = post[key]
        trait = key[2]
        mapped = meta["mapped"]
        family = lipid_family(trait, mapped)
        for lead in other:
            labels = []
            for item in coding:
                if item["chrom"] == lead["chrom"] and abs(item["pos"] - lead["pos"]) <= WINDOW:
                    labels.append(item)
            distinct = {item["mapped"] for item in labels}
            if len(distinct) != 1:
                if len(distinct) > 1:
                    ambiguous += 1
                continue
            coding_lead = min(labels, key=lambda item: (abs(item["pos"] - lead["pos"]), item["pos"]))
            tokens = gene_tokens(coding_lead["mapped"])
            if len(tokens) != 1:
                multi_symbol += 1
                symbol = coding_lead["mapped"].strip().upper()
            else:
                symbol = next(iter(tokens))
            ensg = ""
            ensg_parts = [part.strip() for part in coding_lead["ensg"].split(",") if part.strip().startswith("ENSG")]
            if len(ensg_parts) == 1:
                ensg = ensg_parts[0]
            journal = meta["journal"]
            submission = meta["submission"]
            known = [day for day in (parse_day(journal), parse_day(submission)) if day is not None]
            earliest = min(known).strftime("%Y-%m-%d") if known else ""
            rows.append(
                {
                    "pubmed_id": key[0],
                    "journal_date": journal,
                    "submission_date": submission if parse_day(submission) else "",
                    "date_added": meta["added"] if parse_day(meta["added"]) else "",
                    "earliest_public_date": earliest,
                    "trait_name": trait,
                    "mapped_trait": mapped,
                    "trait_group": "blood_lipid" if family else (mapped or trait),
                    "lipid_family": family,
                    "is_lipid": bool(family),
                    "is_hypothyroid": "hypothyroid" in f"{trait} {mapped}".lower(),
                    "pops_trait": pops_trait(trait, mapped),
                    "study_accessions": ",".join(sorted(set(meta["accessions"]))),
                    "study_accession": sorted(set(meta["accessions"]))[0],
                    "has_sumstats": meta["sumstats_any"],
                    "has_sumstats_first_row": meta["sumstats_first"],
                    "noncoding_rsid": lead["rsid"],
                    "noncoding_chrom": lead["chrom"],
                    "noncoding_pos": lead["pos"],
                    "noncoding_allele": lead["allele"],
                    "coding_rsid": coding_lead["rsid"],
                    "coding_chrom": coding_lead["chrom"],
                    "coding_pos": coding_lead["pos"],
                    "coding_allele": coding_lead["allele"],
                    "label_symbol": symbol,
                    "label_ensembl": ensg,
                    "label_mapped_raw": coding_lead["mapped"],
                    "single_label": len(tokens) == 1,
                }
            )
    print(f"clean_neighborhoods {len(rows)} ambiguous_neighborhoods {ambiguous} multi_symbol_labels {multi_symbol}", flush=True)
    return rows


def index_pips() -> dict[tuple, list[dict]]:
    index = defaultdict(list)
    files = sorted(PIPS.glob("GCST*.tsv"))
    for path in files:
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["finemapping_method"] != "SuSiE-inf":
                    continue
                index[(row["chromosome"].strip(), int(row["position"]))].append(
                    {
                        "accession": row["study_id"],
                        "study_locus_id": row["study_locus_id"],
                        "ref": row["ref"].upper(),
                        "alt": row["alt"].upper(),
                        "pip": float(row["posterior_probability"] or 0),
                        "source": "local",
                    }
                )
    print(f"local_susie_variant_rows {sum(len(v) for v in index.values())} files {len(files)}", flush=True)
    return index


def post_ot(query: str, variables: dict) -> dict:
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        API,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": "exam-credible-sets"},
    )
    last = ""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if payload.get("errors"):
                raise RuntimeError(json.dumps(payload["errors"])[:400])
            return payload["data"]
        except urllib.error.HTTPError as exc:
            last = exc.read().decode("utf-8", errors="replace")[:400]
            if exc.code not in (429, 500, 502, 503):
                raise RuntimeError(f"API {exc.code}: {last}") from exc
            time.sleep(2 + attempt * 3)
        except urllib.error.URLError:
            time.sleep(2 + attempt * 3)
    raise RuntimeError(f"API failed: {last}")


def fetch_ot_susie(accessions: list[str], index: dict) -> None:
    """Add SuSiE-inf variants for studies that are not in the local download."""
    sets_query = """
    query($studyId: String!, $index: Int!) {
      credibleSets(studyIds: [$studyId], page: {index: $index, size: 50}) {
        count
        rows { studyLocusId studyId finemappingMethod }
      }
    }
    """
    locus_query = """
    query($studyLocusId: String!, $index: Int!) {
      credibleSet(studyLocusId: $studyLocusId) {
        locus(page: {index: $index, size: 100}) {
          count
          rows {
            posteriorProbability
            variant { id chromosome position referenceAllele alternateAllele }
          }
        }
      }
    }
    """
    cache = OUT_DIR / "ot_extra_susie.tsv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        with cache.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                index[(row["chromosome"], int(row["position"]))].append(
                    {
                        "accession": row["study_id"],
                        "study_locus_id": row["study_locus_id"],
                        "ref": row["ref"].upper(),
                        "alt": row["alt"].upper(),
                        "pip": float(row["posterior_probability"] or 0),
                        "source": "ot_extra",
                    }
                )
        print(f"ot_extra_cache_rows {sum(1 for _ in cache.open(encoding='utf-8')) - 1}", flush=True)
        return
    written = []
    for number, accession in enumerate(accessions, start=1):
        sets = []
        page = 0
        total = None
        while total is None or len(sets) < total:
            data = post_ot(sets_query, {"studyId": accession, "index": page})
            block = data["credibleSets"]
            total = int(block["count"])
            sets.extend(block["rows"])
            page += 1
            if not block["rows"]:
                break
        susie = [item for item in sets if item["finemappingMethod"] == "SuSiE-inf"]
        print(f"OT {number}/{len(accessions)} {accession} sets {len(sets)} susie {len(susie)}", flush=True)
        for item in susie:
            variants = []
            page = 0
            total = None
            while total is None or len(variants) < total:
                data = post_ot(locus_query, {"studyLocusId": item["studyLocusId"], "index": page})
                locus = data["credibleSet"]["locus"]
                total = int(locus["count"])
                variants.extend(locus["rows"])
                page += 1
                if not locus["rows"]:
                    break
            for variant in variants:
                info = variant["variant"]
                rec = {
                    "study_id": accession,
                    "study_locus_id": item["studyLocusId"],
                    "chromosome": str(info["chromosome"]),
                    "position": str(info["position"]),
                    "ref": info["referenceAllele"].upper(),
                    "alt": info["alternateAllele"].upper(),
                    "posterior_probability": variant["posteriorProbability"],
                }
                written.append(rec)
                index[(rec["chromosome"], int(rec["position"]))].append(
                    {
                        "accession": accession,
                        "study_locus_id": rec["study_locus_id"],
                        "ref": rec["ref"],
                        "alt": rec["alt"],
                        "pip": float(rec["posterior_probability"] or 0),
                        "source": "ot_extra",
                    }
                )
    with cache.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["study_id", "study_locus_id", "chromosome", "position", "ref", "alt", "posterior_probability"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(written)
    print(f"ot_extra_wrote {len(written)}", flush=True)


def match_susie(row: dict, index: dict) -> None:
    accessions = set(row["study_accessions"].split(","))
    hits = []
    for item in index.get((row["noncoding_chrom"], row["noncoding_pos"]), []):
        if item["accession"] not in accessions:
            continue
        allele = row["noncoding_allele"]
        allele_ok = (not allele) or allele in {item["ref"], item["alt"]}
        if allele and not allele_ok:
            continue
        hits.append(item)
    row["has_susie_noncoding_lead"] = bool(hits)
    if not hits:
        row["susie_accession"] = ""
        row["susie_locus_id"] = ""
        row["susie_source"] = ""
        return
    best = max(hits, key=lambda item: item["pip"])
    row["susie_accession"] = best["accession"]
    row["susie_locus_id"] = best["study_locus_id"]
    row["susie_source"] = best["source"]


def mark_scores(rows: list[dict]) -> None:
    have = {path.stem for path in SCORES.glob("GCST*.tsv") if path.stat().st_size > 0}
    for row in rows:
        accessions = set(row["study_accessions"].split(","))
        scored = sorted(accessions & have)
        if scored:
            row["has_trait_score"] = True
            row["score_source"] = "landstra_file"
            row["score_accession"] = scored[0]
        elif row["pops_trait"]:
            row["has_trait_score"] = True
            row["score_source"] = "pops_full"
            row["score_accession"] = ""
        else:
            row["has_trait_score"] = False
            row["score_source"] = ""
            row["score_accession"] = ""


def previously_reported(rows: list[dict], all_traits: dict) -> None:
    symbols = {row["label_symbol"] for row in rows}
    hit_trait = defaultdict(set)
    hit_mapped = defaultdict(set)
    hit_family = defaultdict(set)
    with zipfile.ZipFile(ASSOC_ZIP) as zf:
        name = zf.namelist()[0]
        with zf.open(name) as raw:
            text = (line.decode("utf-8", errors="replace") for line in raw)
            for assoc in csv.DictReader(text, delimiter="\t"):
                day = parse_day(assoc.get("DATE") or "")
                if day is None or day >= REPORTED_BEFORE:
                    continue
                tokens = gene_tokens(assoc.get("MAPPED_GENE") or "") | gene_tokens(assoc.get("REPORTED GENE(S)") or "")
                tokens &= symbols
                if not tokens:
                    continue
                key = (
                    (assoc.get("PUBMEDID") or "").strip(),
                    (assoc.get("DATE") or "").strip(),
                    (assoc.get("DISEASE/TRAIT") or "").strip(),
                )
                meta = all_traits.get(key)
                if meta is None or meta["journal"] >= REPORTED_BEFORE:
                    continue
                trait = key[2]
                for symbol in tokens:
                    hit_trait[symbol].add(trait)
                    if meta["mapped"]:
                        hit_mapped[symbol].add(meta["mapped"])
                    if meta["family"]:
                        hit_family[symbol].add(meta["family"])
    for row in rows:
        symbol = row["label_symbol"]
        if row["is_lipid"]:
            row["previously_reported"] = row["lipid_family"] in hit_family[symbol]
        else:
            row["previously_reported"] = (
                row["trait_name"] in hit_trait[symbol] or row["mapped_trait"] in hit_mapped[symbol]
            )


def add_gene_lists(rows: list[dict], by_chrom) -> None:
    symbol_to_ensg = defaultdict(list)
    for chrom, genes in by_chrom.items():
        for _, _, symbol, ensg in genes:
            symbol_to_ensg[(chrom, symbol)].append(ensg)
    for row in rows:
        near = genes_near(by_chrom, row["noncoding_chrom"], row["noncoding_pos"])
        symbols = sorted({symbol for _, symbol, _ in near})
        row["candidate_genes"] = ",".join(f"{{{symbol}}}" for symbol in symbols)
        row["n_candidates"] = len(symbols)
        if near:
            row["nearest_gene"] = min(near, key=lambda item: (item[0], item[1]))[1]
        else:
            row["nearest_gene"] = ""
        if not row["label_ensembl"]:
            matches = symbol_to_ensg.get((row["noncoding_chrom"], row["label_symbol"]), [])
            if len(matches) == 1:
                row["label_ensembl"] = matches[0]


def collapse(rows: list[dict]) -> list[dict]:
    rows = [row for row in rows if row["single_label"]]
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
        indexes.sort(key=lambda item: rows[item]["noncoding_pos"])
        start = 0
        for end, index in enumerate(indexes):
            while rows[index]["noncoding_pos"] - rows[indexes[start]]["noncoding_pos"] > WINDOW:
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
            lead = row["has_trait_score"] and row["has_susie_noncoding_lead"]
            study = row["has_trait_score"] and row["study_has_susie"]
            return (0 if lead else 1 if study else 2, row["study_accession"], row["noncoding_pos"])

        chosen = min(members, key=sort_key)
        out = dict(chosen)
        out["neighborhood_id"] = f"A{number:05d}"
        out["n_collapsed"] = len(members)
        out["member_accessions"] = ",".join(sorted({row["study_accession"] for row in members}))
        out["in_exam_b"] = bool(chosen["has_trait_score"] and chosen["study_has_susie"] and not chosen["is_hypothyroid"])
        out["hypothyroid_group"] = chosen["is_hypothyroid"]
        out["has_zoodsma"] = any(row["pubmed_id"] == "41044249" for row in members)
        out["has_landstra"] = any(row["pubmed_id"] == "42069741" for row in members)
        out["has_ebmd"] = any(row["pops_trait"] == "eBMD" for row in members)
        out["has_adhd"] = any(row["pops_trait"] == "ADHD" for row in members)
        collapsed.append(out)
    return collapsed


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("has_sumstats", "has_sumstats_first_row", "has_trait_score", "has_susie_noncoding_lead", "study_has_susie", "previously_reported", "is_lipid", "is_hypothyroid", "single_label", "in_exam_b", "hypothyroid_group"):
                if key in out:
                    out[key] = "1" if out[key] else "0"
            writer.writerow(out)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    post, all_traits = load_studies()
    by_key, _, _ = load_associations(post)
    rows = build_clean(post, by_key)
    studies_clean = {row["study_accession"] for row in rows}
    papers = {row["pubmed_id"] for row in rows}
    sumstats_rows = sum(1 for row in rows if row["has_sumstats_first_row"])
    print(f"studies_with_clean_accessions {len(studies_clean)} papers {len(papers)} sumstats_first_row_neighborhoods {sumstats_rows}", flush=True)

    index = index_pips()
    mark_scores(rows)
    need_ot = sorted(
        {
            accession
            for row in rows
            if row["score_source"] == "pops_full"
            for accession in row["study_accessions"].split(",")
            if accession and not (PIPS / f"{accession}.tsv").exists()
        }
    )
    print(f"ot_studies_to_fetch {len(need_ot)}", flush=True)
    if need_ot:
        fetch_ot_susie(need_ot, index)
    susie_studies = {item["accession"] for hits in index.values() for item in hits}
    for row in rows:
        accessions = set(row["study_accessions"].split(","))
        row["study_has_susie"] = bool(accessions & susie_studies)
        match_susie(row, index)
    print("previously_reported pass", flush=True)
    previously_reported(rows, all_traits)
    print("gene lists", flush=True)
    add_gene_lists(rows, load_genes())
    for number, row in enumerate(rows, start=1):
        row["neighborhood_id"] = f"N{number:05d}"

    raw_fields = [
        "neighborhood_id", "study_accession", "study_accessions", "pubmed_id",
        "earliest_public_date", "journal_date", "submission_date", "date_added",
        "trait_name", "mapped_trait", "trait_group", "lipid_family",
        "noncoding_rsid", "noncoding_chrom", "noncoding_pos", "noncoding_allele",
        "coding_rsid", "coding_pos", "coding_allele", "label_symbol", "label_ensembl",
        "candidate_genes", "n_candidates", "nearest_gene",
        "has_sumstats", "has_trait_score", "score_source", "score_accession",
        "has_susie_noncoding_lead", "study_has_susie", "susie_accession", "susie_locus_id", "susie_source",
        "previously_reported", "pops_trait", "is_lipid", "is_hypothyroid", "single_label",
    ]
    write_table(OUT_DIR / "neighborhoods_all.tsv", rows, raw_fields)
    collapsed = collapse(rows)
    exam_fields = raw_fields + [
        "n_collapsed", "member_accessions", "in_exam_b", "hypothyroid_group",
        "has_zoodsma", "has_landstra", "has_ebmd", "has_adhd",
    ]
    write_table(OUT_DIR / "exam_a.tsv", collapsed, exam_fields)

    def count(items, pred):
        return sum(1 for item in items if pred(item))

    exam_b = [row for row in collapsed if row["in_exam_b"]]
    summary = {
        "clean_neighborhoods": len(rows),
        "single_gene_labels": count(rows, lambda row: row["single_label"]),
        "two_gene_mapped_field": count(rows, lambda row: not row["single_label"]),
        "ambiguous_left_out": 1307,
        "exam_a": len(collapsed),
        "exam_a_lipid": count(collapsed, lambda row: row["is_lipid"]),
        "exam_a_non_lipid": count(collapsed, lambda row: not row["is_lipid"]),
        "exam_a_not_previously_reported": count(collapsed, lambda row: not row["previously_reported"]),
        "exam_b": len(exam_b),
        "exam_b_lipid": count(exam_b, lambda row: row["is_lipid"]),
        "exam_b_non_lipid": count(exam_b, lambda row: not row["is_lipid"]),
        "exam_b_not_previously_reported": count(exam_b, lambda row: not row["previously_reported"]),
        "exam_b_lead_in_set": count(exam_b, lambda row: row["has_susie_noncoding_lead"]),
        "exam_b_lead_missing": count(exam_b, lambda row: not row["has_susie_noncoding_lead"]),
        "hypothyroid_group": count(collapsed, lambda row: row["hypothyroid_group"]),
        "hypothyroid_not_previously_reported": count(
            collapsed, lambda row: row["hypothyroid_group"] and not row["previously_reported"]
        ),
        "zoodsma_new_if_scored": count(collapsed, lambda row: row["has_zoodsma"] and not row["in_exam_b"]),
        "bone_new_if_finemapped": count(collapsed, lambda row: row["has_ebmd"] and not row["in_exam_b"]),
        "adhd_new_if_finemapped": count(collapsed, lambda row: row["has_adhd"] and not row["in_exam_b"]),
        "raw_with_susie": count(rows, lambda row: row["has_susie_noncoding_lead"]),
        "raw_with_score": count(rows, lambda row: row["has_trait_score"]),
        "label_ensembl_missing": count(rows, lambda row: not row["label_ensembl"]),
        "nearest_missing": count(rows, lambda row: not row["nearest_gene"]),
        "candidate_empty": count(rows, lambda row: row["n_candidates"] == 0),
    }
    (OUT_DIR / "neighborhood_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

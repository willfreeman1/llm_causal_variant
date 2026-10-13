"""Download Open Targets variant probabilities for the 219 fat-molecule studies.

Each row is one variant in one fine-mapped set. The probability is that
variant's chance of being the causal DNA change in that set. Positions are
the ones Open Targets publishes, on genome build 38. Gene rankings are not
requested.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
from pipeline.paths import trait_scores

JOBS = trait_scores() / "landstra_jobs.tsv"
OUT_DIR = trait_scores() / "work" / "pips"
SUMMARY = OUT_DIR / "summary.json"
API = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY = """
query($studyId: String!, $index: Int!) {
  credibleSets(studyIds: [$studyId], page: {index: $index, size: 50}) {
    count
    rows {
      studyLocusId
      studyId
      finemappingMethod
      confidence
    }
  }
}
"""

LOCUS_QUERY = """
query($studyLocusId: String!, $index: Int!) {
  credibleSet(studyLocusId: $studyLocusId) {
    locus(page: {index: $index, size: 100}) {
      count
      rows {
        posteriorProbability
        is95CredibleSet
        variant { id chromosome position referenceAllele alternateAllele }
      }
    }
  }
}
"""

HEADER = [
    "study_id",
    "study_locus_id",
    "finemapping_method",
    "confidence",
    "variant_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "posterior_probability",
    "in_95_set",
    "assembly",
]


def post(query: str, variables: dict) -> dict:
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    last = ""
    for attempt in range(5):
        request = urllib.request.Request(
            API,
            data=body,
            headers={"Content-Type": "application/json", "User-Agent": "landstra-pip-download"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if payload.get("errors"):
                raise SystemExit(json.dumps(payload["errors"])[:500])
            return payload["data"]
        except urllib.error.HTTPError as exc:
            last = exc.read().decode("utf-8", errors="replace")[:400]
            if exc.code in (429, 500, 502, 503):
                time.sleep(2 + attempt * 3)
                continue
            raise SystemExit(f"API {exc.code}: {last}") from exc
        except urllib.error.URLError:
            time.sleep(2 + attempt * 3)
    raise SystemExit(f"API failed: {last}")


def locus_rows(study_locus_id: str) -> list[dict]:
    rows: list[dict] = []
    index = 0
    total = None
    while total is None or len(rows) < total:
        data = post(LOCUS_QUERY, {"studyLocusId": study_locus_id, "index": index})
        locus = data["credibleSet"]["locus"]
        total = int(locus["count"])
        page = locus["rows"]
        if not page:
            break
        rows.extend(page)
        index += 1
    return rows


def fetch_study(study_id: str) -> list[dict]:
    sets = []
    index = 0
    total = None
    while total is None or len(sets) < total:
        data = post(QUERY, {"studyId": study_id, "index": index})
        block = data["credibleSets"]
        total = int(block["count"])
        page = block["rows"]
        if not page:
            break
        sets.extend(page)
        index += 1
    records = []
    for item in sets:
        variants = locus_rows(item["studyLocusId"])
        for variant in variants:
            info = variant["variant"]
            records.append(
                {
                    "study_id": item["studyId"],
                    "study_locus_id": item["studyLocusId"],
                    "finemapping_method": item.get("finemappingMethod") or "",
                    "confidence": item.get("confidence") or "",
                    "variant_id": info["id"],
                    "chromosome": info["chromosome"],
                    "position": info["position"],
                    "ref": info["referenceAllele"],
                    "alt": info["alternateAllele"],
                    "posterior_probability": variant["posteriorProbability"],
                    "in_95_set": variant["is95CredibleSet"],
                    "assembly": "GRCh38",
                }
            )
    return records


def write_study(study_id: str, records: list[dict]) -> None:
    path = OUT_DIR / f"{study_id}.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    study_ids = []
    for line in JOBS.read_text(encoding="utf-8").splitlines():
        if line.strip():
            study_ids.append(line.split("\t", 1)[0])
    print(f"STUDIES {len(study_ids)}", flush=True)
    done = 0
    for number, study_id in enumerate(study_ids, start=1):
        path = OUT_DIR / f"{study_id}.tsv"
        if path.exists():
            done += 1
            continue
        records = fetch_study(study_id)
        write_study(study_id, records)
        done += 1
        methods = sorted({row["finemapping_method"] or "none" for row in records}) or ["none"]
        print(
            f"{number} {study_id} sets_variants {len(records)} methods {','.join(methods)}",
            flush=True,
        )
        time.sleep(0.15)
    n_variants = 0
    n_empty = 0
    methods: dict[str, int] = {}
    for path in OUT_DIR.glob("GCST*.tsv"):
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if not rows:
            n_empty += 1
            continue
        n_variants += len(rows)
        for row in rows:
            methods[row["finemapping_method"]] = methods.get(row["finemapping_method"], 0) + 1
    summary = {
        "studies": len(study_ids),
        "files": done,
        "empty_studies": n_empty,
        "variant_rows": n_variants,
        "rows_by_method": methods,
        "assembly": "GRCh38",
        "source": API,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

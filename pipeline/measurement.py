"""The measurement sentence pasted above the papers. Same rules as the locked run."""

from __future__ import annotations

import re

CIRCULATING = re.compile(r"^Circulating .+ levels\b", re.I)
MAPPED_PROTEIN = re.compile(r"(level|amount|concentration) of .+ in blood", re.I)


def is_blood_protein(trait: str, mapped: str) -> bool:
    if CIRCULATING.match(trait or ""):
        return True
    return bool(MAPPED_PROTEIN.search(mapped or ""))


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

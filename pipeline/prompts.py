"""Frozen prompts. The system text is the one that produced the locked scores."""

from __future__ import annotations

from pipeline.model_api import SYSTEM
ABSTRACT_CAP = 1800


def no_paper_message(trait: str, gene_string: str) -> str:
    return f"Identify the causal gene.\nGWAS phenotype: {trait}\nGenes in locus: {gene_string}\n"


def paper_message(trait: str, measurement: str, gene_string: str, abstracts: list[dict]) -> str:
    lines = [
        "Identify the causal gene.",
        f"GWAS phenotype: {trait}",
        f"Measurement: {measurement}",
        f"Genes in locus: {gene_string}",
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


def model_payload(state: dict, abstracts: list[dict] | None = None) -> dict:
    """The only fields a model call is allowed to see."""
    abstracts = abstracts or []
    trait = state["trait"]
    gene_string = state["gene_string"]
    if abstracts:
        user = paper_message(trait, state["measurement"], gene_string, abstracts)
    else:
        user = no_paper_message(trait, gene_string)
    return {"system": SYSTEM, "user": user}

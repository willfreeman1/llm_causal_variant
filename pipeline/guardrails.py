"""Checks that run before a model call and after each answer."""

from __future__ import annotations

from collections import defaultdict

from pipeline.model_api import parse_gene

SURENESS_REVIEW = 0.70
CALDERA_REVIEW = 0.30
CODING_PIP_LIMIT = 0.5


def blank_if_outside(gene: str, confidence: float | None, candidates: list[str]) -> tuple[str, float, str]:
    """An answer outside the candidate list becomes blank and is logged."""
    predicted = (gene or "").strip()
    allowed = {item.upper() for item in candidates}
    conf = float(confidence) if confidence is not None else 0.0
    if predicted and predicted.upper() not in allowed:
        return "", conf, predicted
    return predicted, conf, ""


def parse_model_text(text: str, candidates: list[str]) -> tuple[str, float, str]:
    gene, confidence = parse_gene(text or "")
    return blank_if_outside(gene, confidence, candidates)


def input_refusal(
    candidates: list[str],
    cutoff: str,
    caldera_rows: list[dict] | None = None,
    question_locus: str = "",
) -> str:
    """Return a reason to refuse, or an empty string when the inputs are allowed.

    The check uses only the credible-set file. It does not look at the answer.
    """
    if not candidates:
        return "candidate list is empty"
    if not str(cutoff or "").strip():
        return "study date is missing"
    rows = caldera_rows or []
    coding_pip: dict[str, float] = defaultdict(float)
    for row in rows:
        locus = str(row.get("locus") or "")
        if question_locus and locus and locus != question_locus:
            return "credible set is not the noncoding lead's set"
        coded = str(row.get("c_gene") or "").strip()
        if not coded or coded.upper() in ("NA", "NAN"):
            continue
        coding_pip[coded] += float(row.get("pip") or 0)
    if any(total > CODING_PIP_LIMIT for total in coding_pip.values()):
        return "protein-changing probability is above 0.5"
    return ""

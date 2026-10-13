"""The locked combining rule. The cutoffs are fixed. Nothing here is tuned."""

from __future__ import annotations

from pipeline.guardrails import CALDERA_REVIEW, SURENESS_REVIEW

LINE = 0.85


def combining_rule(gpt_gene: str, gpt_confidence: float, caldera_gene: str, line: float = LINE) -> tuple[str, str]:
    """Return the pick and which method it came from: agreement, gpt, or caldera."""
    gpt_gene = (gpt_gene or "").strip()
    caldera_gene = (caldera_gene or "").strip()
    if gpt_gene and caldera_gene and gpt_gene.upper() == caldera_gene.upper():
        return gpt_gene, "agreement"
    if gpt_gene and float(gpt_confidence or 0) >= line:
        return gpt_gene, "gpt"
    return caldera_gene, "caldera"


def final_decision(state: dict) -> tuple[str, str, bool]:
    """Return the pick, the source, and the human-review flag. The flag does not change the pick."""
    mode = state["mode"]
    if mode == "llm_only":
        if state.get("papers") == "on":
            gene = state.get("paper_gene") or ""
            confidence = float(state.get("paper_confidence") or 0)
            source = "gpt_papers"
        else:
            gene = state.get("gpt_gene") or ""
            confidence = float(state.get("gpt_confidence") or 0)
            source = "gpt"
        return gene, source, confidence < SURENESS_REVIEW
    use_gene = state.get("gpt_gene") or ""
    use_conf = float(state.get("gpt_confidence") or 0)
    if state.get("papers") == "on":
        agreed = bool(use_gene and state.get("caldera_gene") and use_gene.upper() == state["caldera_gene"].upper())
        if not agreed:
            use_gene = state.get("paper_gene") or ""
            use_conf = float(state.get("paper_confidence") or 0)
    line = float(state.get("sureness_line") or LINE)
    pick, source = combining_rule(use_gene, use_conf, state.get("caldera_gene") or "", line)
    review = source == "caldera" and float(state.get("caldera_probability") or 0) < CALDERA_REVIEW
    return pick, source, review

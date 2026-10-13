import pytest

from pipeline import spend
from pipeline.graph import node_llm_no_papers


def test_spend_cap_stops_before_the_call(tmp_path, monkeypatch):
    ledger = tmp_path / "spend.json"
    monkeypatch.setattr(spend, "LEDGER", ledger)
    spend.add_spent(1.0, ledger)
    monkeypatch.setattr("pipeline.graph.over_cap", lambda cap, path=None: True)
    state = {
        "replay": False,
        "spend_cap": 0.5,
        "trait": "Albumin",
        "gene_string": "{ALB}",
        "candidates": ["ALB"],
    }
    with pytest.raises(SystemExit, match="spend cap"):
        node_llm_no_papers(state)

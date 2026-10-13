from pipeline.graph import run_question, write_mermaid
from pipeline.prompts import no_paper_message


def test_replay_uses_the_saved_answer_and_blanks_an_outside_gene(tmp_path):
    trace = tmp_path / "trace.jsonl"
    result = run_question(
        {
            "question_id": "demo",
            "mode": "llm_only",
            "papers": "off",
            "replay": True,
            "trait": "Circulating CRTAM levels",
            "gene_string": "{ASGR1},{GNPTAB}",
            "cutoff": "2025-09-01",
            "saved_gpt_raw": "CRTAM",
            "saved_gpt_confidence": 0.91,
            "trace_path": str(trace),
        }
    )
    assert result["gpt_gene"] == ""
    assert result["gpt_raw_logged"] == "CRTAM"
    assert result["final_pick"] == ""
    assert result["nopaper_message"] == no_paper_message(
        "Circulating CRTAM levels", "{ASGR1},{GNPTAB}"
    )
    assert trace.read_text(encoding="utf-8").strip()


def test_graph_drawing_is_written(tmp_path):
    path = tmp_path / "graph.md"
    text = write_mermaid(path)
    assert "guardrail_inputs" in text
    assert "retrieve" in text
    assert path.exists()

import json

from pipeline.guardrails import blank_if_outside, input_refusal, parse_model_text
from pipeline.prompts import model_payload


def test_answer_outside_the_list_becomes_blank():
    gene, confidence, raw = blank_if_outside("CRTAM", 0.8, ["ASGR1", "GNPTAB"])
    assert gene == ""
    assert raw == "CRTAM"
    assert confidence == 0.8


def test_answer_in_the_list_is_kept():
    gene, _confidence, raw = blank_if_outside("asgr1", 0.4, ["ASGR1", "GNPTAB"])
    assert gene == "asgr1"
    assert raw == ""


def test_parse_logs_the_raw_name():
    text = json.dumps({"causal_gene": "CRTAM", "confidence": 0.9, "reason": "the protein"})
    gene, _confidence, raw = parse_model_text(text, ["ASGR1"])
    assert gene == ""
    assert raw == "CRTAM"


def test_refuses_an_extra_credible_set():
    reason = input_refusal(
        ["APOE"],
        "2025-09-01",
        [{"locus": "OTHER", "c_gene": "", "pip": 0.1}, {"locus": "A1", "c_gene": "", "pip": 0.1}],
        "A1",
    )
    assert reason == "credible set is not the noncoding lead's set"


def test_refuses_protein_changing_pip_above_half():
    reason = input_refusal(
        ["APOE"],
        "2025-09-01",
        [{"locus": "A1", "c_gene": "ENSG1", "pip": 0.6}],
        "A1",
    )
    assert reason == "protein-changing probability is above 0.5"


def test_sums_protein_changing_probabilities_for_one_gene():
    reason = input_refusal(
        ["APOE"],
        "2025-09-01",
        [
            {"locus": "A1", "c_gene": "ENSG1", "pip": 0.3},
            {"locus": "A1", "c_gene": "ENSG1", "pip": 0.3},
        ],
        "A1",
    )
    assert reason == "protein-changing probability is above 0.5"


def test_accepts_the_trem2_probability():
    reason = input_refusal(
        ["TREM2"],
        "2025-09-01",
        [{"locus": "A00501", "c_gene": "ENSG00000095970", "pip": 0.41}],
        "A00501",
    )
    assert reason == ""


def test_guardrail_takes_no_label():
    import inspect

    names = set(inspect.signature(input_refusal).parameters)
    assert not any("label" in name for name in names)


def test_label_is_not_in_the_model_payload():
    state = {
        "trait": "Alzheimer's disease",
        "gene_string": "{APOE},{TREM2}",
        "measurement": "Alzheimer's disease. Catalog term: Alzheimer disease.",
        "labeled_ensembl": "ENSG00000095970",
        "label_symbol": "TREM2",
    }
    payload = model_payload(state)
    blob = json.dumps(payload)
    assert "labeled_ensembl" not in payload
    assert "label_symbol" not in payload
    assert "ENSG00000095970" not in blob
    assert payload["user"].startswith("Identify the causal gene.\n")


def test_empty_list_and_missing_date_are_refused():
    assert input_refusal([], "2025-09-01") == "candidate list is empty"
    assert input_refusal(["APOE"], "") == "study date is missing"

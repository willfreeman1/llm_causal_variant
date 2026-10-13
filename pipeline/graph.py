"""LangGraph for the locked method. Replay does not call the API."""

from __future__ import annotations

import json
import os
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import TypedDict

from langgraph.graph import END, START, StateGraph

from pipeline.model_api import (
    INPUT_USD_PER_TOKEN,
    MODEL,
    OUTPUT_USD_PER_TOKEN,
    SYSTEM,
    load_api_key,
    parse_gene,
)

from pipeline.decide import final_decision
from pipeline.paths import LOCAL, ROOT
from pipeline.guardrails import blank_if_outside, input_refusal, parse_model_text
from pipeline.prompts import model_payload
from pipeline.retrieve import search
from pipeline.spend import add_spent, over_cap

TRACE_DIR = LOCAL / "ukb" / "final_exam" / "traces"


class GraphState(TypedDict, total=False):
    question_id: str
    mode: str
    papers: str
    replay: bool
    trait: str
    gene_string: str
    candidates: list[str]
    cutoff: str
    date_mode: str
    lead_rsid: str
    lead_position: str
    measurement: str
    query_vector: list[float]
    sureness_line: float
    question_locus: str
    caldera_cs_path: str
    caldera_pops_path: str
    caldera_rows: list[dict]
    caldera_gene: str
    caldera_probability: float
    saved_gpt_raw: str
    saved_gpt_confidence: float
    saved_paper_raw: str
    saved_paper_confidence: float
    expected_pmids: list[int]
    gpt_gene: str
    gpt_confidence: float
    gpt_raw_logged: str
    paper_gene: str
    paper_confidence: float
    paper_raw_logged: str
    pmids: list[int]
    pmid_mismatch: int
    nopaper_message: str
    paper_message: str
    guardrail_events: list[str]
    path: str
    final_pick: str
    pick_source: str
    needs_review: bool
    refused: bool
    cost_usd: float
    elapsed_s: float
    spend_cap: float
    trace_path: str
    pending_slot: str
    pending_gene: str
    pending_confidence: float
    pending_raw: str
    retrieved_docs: list[dict]
    trace_record: dict
    _started: float


def _events(state: GraphState, extra: str) -> list[str]:
    events = list(state.get("guardrail_events") or [])
    if extra:
        events.append(extra)
    return events


def node_inputs(state: GraphState) -> dict:
    reason = input_refusal(
        state.get("candidates") or [],
        state.get("cutoff") or "",
        state.get("caldera_rows"),
        state.get("question_locus") or "",
    )
    if reason:
        return {"refused": True, "guardrail_events": _events(state, reason), "path": "refused"}
    return {"refused": False, "guardrail_events": list(state.get("guardrail_events") or [])}


def _rscript() -> Path:
    found = shutil.which("Rscript")
    if found:
        return Path(found)
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    matches = sorted((program_files / "R").glob("R-*/bin/Rscript.exe"))
    if matches:
        return matches[-1]
    raise RuntimeError("Rscript was not found. A live statistical-model call needs R.")


def node_caldera(state: GraphState) -> dict:
    if state.get("replay") or state.get("caldera_gene"):
        return {}
    cs_path = state.get("caldera_cs_path") or ""
    pops_path = state.get("caldera_pops_path") or ""
    if not cs_path or not pops_path:
        return {"caldera_gene": "", "caldera_probability": 0.0, "guardrail_events": _events(state, "no caldera output")}
    import subprocess

    out = LOCAL / "ukb" / "final_exam" / "caldera_live_check.tsv"
    command = [
        str(_rscript()),
        str(ROOT / "pipeline" / "caldera_one.R"),
        str(ROOT).replace("\\", "/"),
        str(cs_path).replace("\\", "/"),
        str(pops_path).replace("\\", "/"),
        str(out).replace("\\", "/"),
    ]
    subprocess.run(command, check=True)
    import csv

    best_gene = ""
    best_p = -1.0
    locus = state.get("question_locus") or ""
    with out.open(encoding="utf-8", newline="") as handle:
        for item in csv.DictReader(handle, delimiter="\t"):
            if locus and item.get("locus") != locus:
                continue
            value = float(item["caldera"])
            if value > best_p:
                best_p = value
                best_gene = item["gene"]
    return {"caldera_gene": best_gene, "caldera_probability": best_p if best_p >= 0 else 0.0}


def _call_model(state: GraphState, user: str) -> dict:
    cap = float(state.get("spend_cap") or 0)
    if cap and over_cap(cap):
        raise SystemExit(f"spend cap ${cap:.2f} reached before the call")
    body = {
        "model": MODEL,
        "temperature": 0,
        "max_completion_tokens": 500,
        "reasoning_effort": "none",
        "messages": [
            {"role": "system", "content": SYSTEM},
            {"role": "user", "content": user},
        ],
    }
    key = load_api_key()
    last = ""
    for attempt in range(5):
        request = urllib.request.Request(
            "https://api.openai.com/v1/chat/completions",
            data=json.dumps(body).encode("utf-8"),
            headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            break
        except urllib.error.HTTPError as exc:
            last = exc.read().decode("utf-8", errors="replace")[:500]
            if exc.code == 400 and "temperature" in last and "temperature" in body:
                body.pop("temperature", None)
                continue
            if exc.code == 400 and "reasoning_effort" in last and "reasoning_effort" in body:
                body.pop("reasoning_effort", None)
                continue
            if exc.code in (429, 500, 502, 503):
                time.sleep(2 + attempt * 3)
                continue
            raise SystemExit(f"API {exc.code}: {last}") from exc
    else:
        raise SystemExit(f"API failed after retries: {last}")
    usage = payload.get("usage") or {}
    prompt = int(usage.get("prompt_tokens") or 0)
    completion = int(usage.get("completion_tokens") or 0)
    cost = prompt * INPUT_USD_PER_TOKEN + completion * OUTPUT_USD_PER_TOKEN
    add_spent(cost)
    text = payload["choices"][0]["message"].get("content") or ""
    return {"text": text, "cost_usd": cost}


def node_llm_no_papers(state: GraphState) -> dict:
    user = model_payload(state)["user"]
    if state.get("replay"):
        return {
            "pending_slot": "gpt",
            "pending_gene": state.get("saved_gpt_raw") or "",
            "pending_confidence": float(state.get("saved_gpt_confidence") or 0),
            "pending_raw": state.get("saved_gpt_raw") or "",
            "nopaper_message": user,
        }
    if state.get("spend_cap") and over_cap(float(state["spend_cap"])):
        raise SystemExit("spend cap reached before the no-paper call")
    result = _call_model(state, user)
    gene, confidence = parse_gene(result["text"])
    return {
        "pending_slot": "gpt",
        "pending_gene": gene,
        "pending_confidence": confidence if confidence is not None else 0.0,
        "pending_raw": gene,
        "nopaper_message": user,
        "cost_usd": float(state.get("cost_usd") or 0) + result["cost_usd"],
    }


def node_llm_with_papers(state: GraphState) -> dict:
    docs = state.get("retrieved_docs") or []
    user = model_payload(state, docs)["user"]
    if state.get("replay"):
        return {
            "pending_slot": "paper",
            "pending_gene": state.get("saved_paper_raw") or "",
            "pending_confidence": float(state.get("saved_paper_confidence") or 0),
            "pending_raw": state.get("saved_paper_raw") or "",
            "paper_message": user,
        }
    if state.get("spend_cap") and over_cap(float(state["spend_cap"])):
        raise SystemExit("spend cap reached before the paper call")
    result = _call_model(state, user)
    gene, confidence = parse_gene(result["text"])
    return {
        "pending_slot": "paper",
        "pending_gene": gene,
        "pending_confidence": confidence if confidence is not None else 0.0,
        "pending_raw": gene,
        "paper_message": user,
        "cost_usd": float(state.get("cost_usd") or 0) + result["cost_usd"],
    }


def node_guardrail_output(state: GraphState) -> dict:
    gene, confidence, raw = blank_if_outside(
        state.get("pending_gene") or "",
        state.get("pending_confidence"),
        state.get("candidates") or [],
    )
    events = list(state.get("guardrail_events") or [])
    if raw:
        events.append(f"out of list: {raw}")
    slot = state.get("pending_slot") or "gpt"
    if slot == "paper":
        return {
            "paper_gene": gene,
            "paper_confidence": confidence,
            "paper_raw_logged": raw,
            "guardrail_events": events,
        }
    return {
        "gpt_gene": gene,
        "gpt_confidence": confidence,
        "gpt_raw_logged": raw,
        "guardrail_events": events,
    }


def node_retrieve(state: GraphState) -> dict:
    docs = search(
        state.get("query_vector") or [],
        state.get("candidates") or [],
        state.get("cutoff") or "",
        state.get("lead_rsid") or "",
        state.get("lead_position") or "",
        state.get("date_mode") or "study",
    )
    pmids = [doc["pmid"] for doc in docs]
    expected = state.get("expected_pmids") or []
    mismatch = 0 if not expected else int(pmids != list(expected))
    return {"pmids": pmids, "retrieved_docs": docs, "pmid_mismatch": mismatch, "path": "retrieval"}


def node_decide(state: GraphState) -> dict:
    pick, source, _review = final_decision(state)
    path = state.get("path") or "agreement"
    if path != "retrieval":
        path = "agreement" if source == "agreement" else "no_papers"
    return {"final_pick": pick, "pick_source": source, "path": path}


def node_needs_review(state: GraphState) -> dict:
    _pick, _source, review = final_decision(state)
    return {"needs_review": review, "final_pick": state.get("final_pick") or ""}


def node_trace(state: GraphState) -> dict:
    started = float(state.get("_started") or time.time())
    record = {
        "question_id": state.get("question_id"),
        "mode": state.get("mode"),
        "papers": state.get("papers"),
        "path": state.get("path"),
        "pmids": state.get("pmids") or [],
        "pmid_mismatch": state.get("pmid_mismatch") or 0,
        "guardrail_events": state.get("guardrail_events") or [],
        "needs_review": bool(state.get("needs_review")),
        "final_pick": state.get("final_pick") or "",
        "pick_source": state.get("pick_source") or "",
        "cost_usd": float(state.get("cost_usd") or 0),
        "elapsed_s": round(time.time() - started, 3),
        "refused": bool(state.get("refused")),
        "calls": {
            "no_papers": state.get("nopaper_message") or "",
            "with_papers": state.get("paper_message") or "",
        },
    }
    path = state.get("trace_path") or ""
    if path:
        handle_path = Path(path)
        handle_path.parent.mkdir(parents=True, exist_ok=True)
        with handle_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")
    return {"elapsed_s": record["elapsed_s"], "trace_record": record}


def after_inputs(state: GraphState) -> str:
    if state.get("refused"):
        return "trace"
    if state.get("mode") == "combined":
        return "caldera"
    return "llm_no_papers"


def after_base(state: GraphState) -> str:
    if state.get("mode") == "llm_only":
        return "retrieve" if state.get("papers") == "on" else "decide"
    gpt = (state.get("gpt_gene") or "").upper()
    caldera = (state.get("caldera_gene") or "").upper()
    if gpt and caldera and gpt == caldera:
        return "decide"
    if state.get("papers") == "on":
        return "retrieve"
    return "decide"


def build_graph():
    graph = StateGraph(GraphState)
    graph.add_node("guardrail_inputs", node_inputs)
    graph.add_node("caldera", node_caldera)
    graph.add_node("llm_no_papers", node_llm_no_papers)
    graph.add_node("guardrail_output", node_guardrail_output)
    graph.add_node("retrieve", node_retrieve)
    graph.add_node("llm_with_papers", node_llm_with_papers)
    graph.add_node("guardrail_output_papers", node_guardrail_output)
    graph.add_node("decide", node_decide)
    graph.add_node("needs_review", node_needs_review)
    graph.add_node("trace", node_trace)
    graph.add_edge(START, "guardrail_inputs")
    graph.add_conditional_edges(
        "guardrail_inputs",
        after_inputs,
        {"trace": "trace", "caldera": "caldera", "llm_no_papers": "llm_no_papers"},
    )
    graph.add_edge("caldera", "llm_no_papers")
    graph.add_edge("llm_no_papers", "guardrail_output")
    graph.add_conditional_edges(
        "guardrail_output",
        after_base,
        {"retrieve": "retrieve", "decide": "decide"},
    )
    graph.add_edge("retrieve", "llm_with_papers")
    graph.add_edge("llm_with_papers", "guardrail_output_papers")
    graph.add_edge("guardrail_output_papers", "decide")
    graph.add_edge("decide", "needs_review")
    graph.add_edge("needs_review", "trace")
    graph.add_edge("trace", END)
    return graph.compile()


def write_mermaid(path: Path | None = None) -> str:
    compiled = build_graph()
    text = compiled.get_graph().draw_mermaid()
    destination = path or (Path(__file__).resolve().parent / "graph.md")
    destination.write_text(
        "# Method graph\n\n"
        "Two modes share this graph. `combined` runs the statistical model and the language model. "
        "`llm_only` skips the statistical model. Papers are retrieved only when that mode asks for them, "
        "and in the combined mode only when the two methods disagree.\n\n"
        "```mermaid\n"
        + text
        + "\n```\n",
        encoding="utf-8",
    )
    return text


APP = None


def run_question(state: dict) -> dict:
    global APP
    if APP is None:
        APP = build_graph()
    incoming = dict(state)
    incoming["_started"] = time.time()
    if "candidates" not in incoming and incoming.get("gene_string"):
        incoming["candidates"] = [
            part.strip().strip("{}")
            for part in incoming["gene_string"].split(",")
            if part.strip(" {}")
        ]
    return APP.invoke(incoming)

"""The locked model call. The wording is the published zero-shot prompt."""

from __future__ import annotations

import json
from pathlib import Path

from pipeline.paths import ROOT

MODEL = "gpt-5.4"
INPUT_USD_PER_TOKEN = 2.50 / 1_000_000
OUTPUT_USD_PER_TOKEN = 15.00 / 1_000_000

# Exact string from the earlier study's basic zero-shot prompt, including the leading newline.
SYSTEM = """
You are an expert in biology and genetics.
Your task is to identify likely causal genes within a locus for a given GWAS phenotype based on literature evidence.
From the list, provide the likely causal gene (matching one of the given genes), confidence (0: very unsure to 1: very confident), and a brief reason (50 words or less) for your choice.
Return your response in JSON format, excluding the GWAS phenotype name and gene list in the locus. JSON keys should be 'causal_gene','confidence','reason'. 
Your response must start with '{' and end with '}'.
"""


def load_api_key() -> str:
    env_path = ROOT / ".env"
    if not env_path.exists():
        raise SystemExit("OPENAI_API_KEY is missing. Copy .env.example to .env and add a key.")
    for line in env_path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#") or "=" not in text:
            continue
        key, value = text.split("=", 1)
        if key.strip() == "OPENAI_API_KEY":
            return value.strip().strip('"').strip("'")
    raise SystemExit("OPENAI_API_KEY is missing. Copy .env.example to .env and add a key.")


def parse_gene(text: str) -> tuple[str, float | None]:
    raw = (text or "").strip()
    if raw.startswith("```"):
        raw = raw.strip("`")
        if raw.lower().startswith("json"):
            raw = raw[4:]
        raw = raw.strip()
    start = raw.find("{")
    end = raw.rfind("}")
    parsed = {}
    if start >= 0 and end > start:
        try:
            parsed = json.loads(raw[start : end + 1])
        except json.JSONDecodeError:
            parsed = {}
    gene = str(parsed.get("causal_gene") or "").strip()
    try:
        confidence = float(parsed.get("confidence"))
    except (TypeError, ValueError):
        confidence = None
    return gene, confidence

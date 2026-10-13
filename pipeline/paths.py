"""Folders a fresh copy of the project can see, and the local data it cannot.

Big downloads stay outside the repository. Set CAUSAL_VARIANT_DATA if they
are not in the data folder next to this code. Set EMBED_CACHE for the
abstract-vector cache used to rebuild the paper database.
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QUESTIONS = ROOT / "questions"
RESULTS = ROOT / "results"
OLDER = ROOT / "older_studies"
EARLY = ROOT / "early_checks"
LOCAL = Path(os.environ.get("CAUSAL_VARIANT_DATA", str(ROOT / "data")))


def embed_cache() -> Path:
    raw = os.environ.get("EMBED_CACHE", "")
    if raw:
        return Path(raw)
    legacy = ROOT / "spikes" / "shringarpure_match" / "embed_cache_small"
    if legacy.exists():
        return legacy
    return LOCAL / "embed_cache"


def trait_scores() -> Path:
    raw = os.environ.get("TRAIT_SCORES_DIR", "")
    if not raw:
        raise SystemExit("Set TRAIT_SCORES_DIR to the folder that holds the trait-score work files.")
    return Path(raw)

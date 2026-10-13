"""Spend for this package only. The locked exam's spend file is left alone."""

from __future__ import annotations

import json
from pathlib import Path

from pipeline.paths import RESULTS

LEDGER = RESULTS / "package_spend.json"


def read_spent(path: Path = LEDGER) -> float:
    if not path.exists():
        return 0.0
    return float(json.loads(path.read_text(encoding="utf-8")).get("usd") or 0)


def add_spent(amount: float, path: Path = LEDGER) -> float:
    total = read_spent(path) + float(amount)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"usd": total}, indent=2), encoding="utf-8")
    return total


def over_cap(cap: float, path: Path = LEDGER) -> bool:
    return read_spent(path) >= float(cap)

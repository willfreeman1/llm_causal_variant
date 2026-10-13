"""Turn abstract text into number-lists, and read the lists already saved.

The new-study run calls these. The number-lists themselves stay in the local
cache folder, because they are large.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

from pipeline.paths import embed_cache

EMBED_MODEL = "text-embedding-3-small"
BATCH = 96
TEXT_CAP = 12000


def cache_dir() -> Path:
    return embed_cache()


def embed_request(api_key: str, texts: list[str]) -> dict:
    body = {"model": EMBED_MODEL, "input": texts}
    last_detail = ""
    for attempt in range(6):
        request = urllib.request.Request(
            "https://api.openai.com/v1/embeddings",
            data=json.dumps(body).encode("utf-8"),
            headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last_detail = exc.read().decode("utf-8", errors="replace")[:500]
            if exc.code in (429, 500, 502, 503):
                time.sleep(2 + attempt * 3)
                continue
            raise SystemExit(f"EMBED_API {exc.code}: {last_detail}") from exc
    raise SystemExit(f"EMBED_API failed after retries: {last_detail}")


def cap_text(text: str) -> str:
    text = text or ""
    if len(text) <= TEXT_CAP:
        return text
    return text[:TEXT_CAP].rsplit(" ", 1)[0]


def chunk_paths() -> list[Path]:
    folder = cache_dir()
    if not folder.exists():
        return []
    return sorted(folder.glob("chunk_*.npz"))


def embedded_pmids() -> set[int]:
    done: set[int] = set()
    for path in chunk_paths():
        done.update(int(pmid) for pmid in np.load(path)["pmids"])
    return done

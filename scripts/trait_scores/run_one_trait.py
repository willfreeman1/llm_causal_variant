"""Score the other chromosomes for the first fat-molecule measurement."""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

root = Path(os.environ['TRAIT_SCORES_DIR'])
magma = str(root / "work" / "magma" / "GCST90830355" / "GCST90830355")
out_dir = root / "work" / "one" / "loco"
score_path = root / "work" / "scores" / "GCST90830355.tsv"
chroms = [str(i) for i in range(1, 23)]


def one(chrom):
    out = out_dir / f"chr{chrom}"
    held = Path(str(out) + ".held.tsv")
    if held.exists() and held.stat().st_size > 0:
        print(f"SKIP {chrom}", flush=True)
        return chrom, 0
    env = os.environ.copy()
    env["POPS_THREADS"] = "4"
    proc = subprocess.run(
        [sys.executable, "-u", str(root / "run_chr.py"), chrom, "29", magma, str(out)],
        env=env,
        capture_output=True,
        text=True,
    )
    print(proc.stdout, end="", flush=True)
    if proc.returncode != 0:
        print(proc.stderr, end="", flush=True)
    return chrom, proc.returncode


with ThreadPoolExecutor(max_workers=3) as pool:
    futures = [pool.submit(one, chrom) for chrom in chroms]
    failed = []
    for future in as_completed(futures):
        chrom, code = future.result()
        print(f"DONE {chrom} {code}", flush=True)
        if code != 0:
            failed.append(chrom)
if failed:
    raise SystemExit(f"FAILED {failed}")

parts = [pd.read_csv(out_dir / f"chr{chrom}.held.tsv", sep="\t") for chrom in chroms]
scores = pd.concat(parts, ignore_index=True)[["ENSGID", "PoPS_Score"]]
values = pd.to_numeric(scores.PoPS_Score, errors="coerce")
print(
    f"SCORED {len(scores)} UNIQUE {scores.ENSGID.nunique()} "
    f"FINITE {int(values.notna().sum())} STD {values.std():.4f}",
    flush=True,
)
if scores.ENSGID.nunique() != len(scores) or not values.notna().all() or values.nunique() < 100:
    raise SystemExit("SCORE_FILE_BAD")
score_path.parent.mkdir(parents=True, exist_ok=True)
scores.to_csv(score_path, sep="\t", index=False)
print("TRAIT_OK", score_path, flush=True)

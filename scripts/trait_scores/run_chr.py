"""Score one held-out chromosome for one measurement."""
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

root = Path(os.environ['TRAIT_SCORES_DIR'])
held = sys.argv[1]
chunks = sys.argv[2]
magma = sys.argv[3]
out_prefix = sys.argv[4]

annot_path = root / "work" / "ref" / "gene_annot.txt"
annot = pd.read_csv(annot_path, sep="\t", dtype=str)
chroms = sorted(set(annot.CHR.astype(str)), key=lambda value: (not value.isdigit(), int(value) if value.isdigit() else value))
train = [chrom for chrom in chroms if chrom != held]
threads = os.environ.get("POPS_THREADS", "4")
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = threads

cmd = [
    sys.executable,
    str(root / "pops" / "pops.py"),
    "--gene_annot_path",
    str(annot_path),
    "--feature_mat_prefix",
    str(root / "work" / "features_munged" / "pops_features"),
    "--num_feature_chunks",
    chunks,
    "--magma_prefix",
    magma,
    "--control_features_path",
    str(root / "work" / "ref" / "control_features.txt"),
    "--out_prefix",
    out_prefix,
    "--feature_selection_chromosomes",
    *train,
    "--training_chromosomes",
    *train,
    "--project_out_covariates_chromosomes",
    *train,
]
print("THREADS", threads, "HELD", held, flush=True)
t0 = time.time()
log_path = Path(out_prefix + ".log")
log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w") as log:
    proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
print(f"POPS_EXIT {proc.returncode} SECONDS {time.time() - t0:.1f}", flush=True)
if proc.returncode != 0:
    raise SystemExit(proc.returncode)
preds = pd.read_csv(out_prefix + ".preds", sep="\t")
held_genes = set(annot.loc[annot.CHR.astype(str) == held, "ENSGID"])
kept = preds[preds.ENSGID.isin(held_genes)].copy()
score = pd.to_numeric(kept.PoPS_Score, errors="coerce")
print(
    f"HELD {held} GENES {len(kept)} FINITE {int(score.notna().sum())} NUNIQUE {int(score.nunique())} "
    f"STD {score.std():.4f}",
    flush=True,
)
if len(kept) < 10 or not score.notna().all() or score.nunique() < 10:
    raise SystemExit("score check failed")
kept.to_csv(out_prefix + ".held.tsv", sep="\t", index=False)
print("CHR_OK", held, flush=True)

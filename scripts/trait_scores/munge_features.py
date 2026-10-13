"""Turn the single gene-property table into the chunked files the scorer reads."""
import gzip
import time
import os
from pathlib import Path

import numpy as np

root = Path(os.environ['TRAIT_SCORES_DIR'])
src = root / "dl" / "PoPS.features.txt.gz"
out_dir = root / "work" / "features_munged"
out_dir.mkdir(parents=True, exist_ok=True)
raw_path = out_dir / "matrix.f64"
rows_path = out_dir / "pops_features.rows.txt"
prefix = out_dir / "pops_features"
chunk = 2000

t0 = time.time()
n = 0
with gzip.open(src, "rt") as handle, rows_path.open("w") as rows_out, raw_path.open("wb") as raw_out:
    header = handle.readline().rstrip("\n").split("\t")
    if header[0] != "ENSGID":
        raise SystemExit(f"expected ENSGID, got {header[0]}")
    cols = header[1:]
    n_cols = len(cols)
    (out_dir / "pops_features.columns.txt").write_text("\n".join(cols) + "\n")
    controls = [name for name in cols if name.endswith(".control")]
    (root / "work" / "ref" / "control_features.txt").write_text("\n".join(controls) + "\n")
    print(f"N_FEATURES {n_cols} N_CONTROLS {len(controls)}", flush=True)
    for line in handle:
        tab = line.find("\t")
        gene = line[:tab]
        cleaned = line[tab + 1 :].replace("False", "0").replace("True", "1")
        vals = np.fromstring(cleaned, sep="\t", dtype=np.float64)
        if vals.size != n_cols:
            raise SystemExit(f"row {gene} has {vals.size} values, expected {n_cols}")
        if not np.isfinite(vals).all():
            raise SystemExit(f"row {gene} has a non-finite value")
        rows_out.write(gene + "\n")
        vals.tofile(raw_out)
        n += 1
        if n % 500 == 0:
            print(f"ROWS {n} PER_SEC {n / (time.time() - t0):.1f}", flush=True)

print(f"PARSE_DONE rows {n} seconds {time.time() - t0:.1f}", flush=True)
mm = np.memmap(raw_path, dtype=np.float64, mode="r", shape=(n, n_cols))
n_chunks = 0
for start in range(0, n_cols, chunk):
    stop = min(start + chunk, n_cols)
    block = np.array(mm[:, start:stop], dtype=np.float64)
    np.save(f"{prefix}.mat.{n_chunks}.npy", block)
    (out_dir / f"pops_features.cols.{n_chunks}.txt").write_text("\n".join(cols[start:stop]) + "\n")
    n_chunks += 1
    print(f"CHUNK {n_chunks}", flush=True)
(out_dir / "pops_features.rows.txt").write_text(rows_path.read_text())
(root / "work" / "ref" / "gene_annot.txt").write_text(
    (root / "pops" / "example" / "data" / "utils" / "gene_annot_jun10.txt").read_text()
)
print(f"N_CHUNKS {n_chunks}", flush=True)
print(f"MUNGE_SECONDS {time.time() - t0:.1f}", flush=True)
del mm
raw_path.unlink()
print("MUNGE_DONE", flush=True)

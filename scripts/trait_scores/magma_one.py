"""Turn one full results table into a gene-association file."""
import gzip
import subprocess
import sys
import time
import os
from pathlib import Path

root = Path(os.environ['TRAIT_SCORES_DIR'])
acc = sys.argv[1]
sample_n = sys.argv[2]
gz = Path(sys.argv[3])
out_dir = root / "work" / "magma" / acc
out_dir.mkdir(parents=True, exist_ok=True)
pval = out_dir / f"{acc}.pval"
genes = out_dir / f"{acc}.genes.out"
if genes.exists() and genes.stat().st_size > 0:
    print("MAGMA_EXISTS", acc, flush=True)
    raise SystemExit(0)

if not pval.exists() or pval.stat().st_size == 0:
    kept = 0
    with gzip.open(gz, "rt") as handle, pval.open("w") as dest:
        header = handle.readline().rstrip("\n").split("\t")
        i_rs = header.index("rsid")
        i_p = header.index("p_value")
        best = {}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rs = parts[i_rs]
            p = parts[i_p]
            if not rs.startswith("rs") or p in {"", "NA", "nan"}:
                continue
            try:
                val = float(p)
            except ValueError:
                continue
            prev = best.get(rs)
            if prev is None or val < prev:
                best[rs] = val
        for rs, val in best.items():
            dest.write(f"{rs} {val}\n")
            kept += 1
    print("PVAL_KEPT", kept, flush=True)

bim = next((root / "work" / "ref").rglob("g1000_eur.bim"))
prefix = str(bim.with_suffix(""))
t0 = time.time()
subprocess.run(
    [
        str(root / "work" / "ref" / "magma.exe"),
        "--bfile",
        prefix,
        "--gene-annot",
        str(root / "work" / "ref" / "annot.genes.annot"),
        "--pval",
        str(pval),
        f"N={sample_n}",
        "--gene-model",
        "snp-wise=mean",
        "--out",
        str(out_dir / acc),
    ],
    check=True,
)
windows_out = out_dir / f"{acc}.genes.out.txt"
if windows_out.exists() and not genes.exists():
    genes.write_bytes(windows_out.read_bytes())
print(f"MAGMA_SECONDS {time.time() - t0:.1f}", flush=True)
print("MAGMA_DONE", acc, flush=True)

"""Unpack the reference panel and map variants onto genes. Run once."""
import subprocess
import zipfile
import os
from pathlib import Path

root = Path(os.environ['TRAIT_SCORES_DIR'])
ref = root / "work" / "ref"
zip_path = root / "dl" / "g1000_eur.zip"
magma = ref / "magma.exe"
annot_src = root / "pops" / "example" / "data" / "utils" / "gene_annot_jun10.txt"

ref.mkdir(parents=True, exist_ok=True)
bim = next(ref.rglob("g1000_eur.bim"), None)
if bim is None:
    print("UNZIP", flush=True)
    with zipfile.ZipFile(zip_path) as handle:
        handle.extractall(ref)
    bim = next(ref.rglob("g1000_eur.bim"))
print("BIM", bim, flush=True)
prefix = bim.with_suffix("")
bdir = bim.parent

snploc = ref / "snploc.txt"
if not snploc.exists():
    n = 0
    with bim.open() as src, snploc.open("w") as dest:
        for line in src:
            parts = line.split()
            dest.write(f"{parts[1]} {parts[0]} {parts[3]}\n")
            n += 1
    print("SNPLOC", n, flush=True)

genes = ref / "genes.loc"
if not genes.exists():
    n = 0
    with annot_src.open() as src, genes.open("w") as dest:
        header = src.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in src:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[idx["CHR"]]
            if chrom in {"M", "MT"}:
                continue
            dest.write(f"{parts[idx['ENSGID']]} {chrom} {parts[idx['START']]} {parts[idx['END']]}\n")
            n += 1
    print("GENES", n, flush=True)

annot_out = ref / "annot.genes.annot"
if not annot_out.exists():
    subprocess.run(
        [
            str(magma),
            "--annotate",
            "window=0",
            "--snp-loc",
            str(snploc),
            "--gene-loc",
            str(genes),
            "--out",
            str(ref / "annot"),
        ],
        check=True,
    )
print("BFILE", prefix, flush=True)
print("PREPARE_DONE", flush=True)

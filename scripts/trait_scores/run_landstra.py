"""Score the remaining Landstra fat-molecule measurements on this computer.

Three gene-association jobs and three chromosome fits run at once, at lower
priority than normal programs, so the machine stays usable. A measurement that
already has a score file is skipped.
"""
import os
import queue
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

import pandas as pd

root = Path(os.environ['TRAIT_SCORES_DIR'])
jobs_path = root / "landstra_jobs.tsv"
score_dir = root / "work" / "scores"
dl_dir = root / "dl" / "sumstats"
BELOW_NORMAL = 0x00004000
N_MAGMA = 3
N_POPS = 3

score_dir.mkdir(parents=True, exist_ok=True)
dl_dir.mkdir(parents=True, exist_ok=True)

magma_q = queue.Queue()
chrom_q = queue.Queue()
lock = threading.Lock()
remaining = {}
failed = set()
done_n = 0
t_start = time.time()


def sample_n(text):
    match = re.search(r"[\d,]+", str(text))
    if match is None:
        raise ValueError(f"no sample size in {text!r}")
    value = int(match.group(0).replace(",", ""))
    if not 1000 <= value <= 100000:
        raise ValueError(f"sample size {value} looks wrong")
    return str(value)


def download(acc, url):
    dest = dl_dir / f"{acc}.h.tsv.gz"
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest
    partial = dest.with_suffix(".partial")
    subprocess.run(
        ["curl.exe", "-L", "--fail", "--retry", "5", "--retry-delay", "10", "-o", str(partial), url],
        check=True,
        creationflags=BELOW_NORMAL,
    )
    partial.replace(dest)
    return dest


def run_magma(acc, n_text, gz):
    genes = root / "work" / "magma" / acc / f"{acc}.genes.out"
    if genes.exists() and genes.stat().st_size > 0:
        print(f"MAGMA_EXISTS {acc}", flush=True)
    else:
        subprocess.run(
            [sys.executable, "-u", str(root / "magma_one.py"), acc, sample_n(n_text), str(gz)],
            check=True,
            creationflags=BELOW_NORMAL,
        )
    if gz.exists():
        gz.unlink()


def run_chrom(acc, chrom):
    out = root / "work" / "loco" / acc / f"chr{chrom}"
    held = Path(str(out) + ".held.tsv")
    if held.exists() and held.stat().st_size > 0:
        return
    env = os.environ.copy()
    env["POPS_THREADS"] = "4"
    magma = str(root / "work" / "magma" / acc / acc)
    subprocess.run(
        [sys.executable, "-u", str(root / "run_chr.py"), str(chrom), "29", magma, str(out)],
        check=True,
        env=env,
        creationflags=BELOW_NORMAL,
    )


def write_score(acc):
    loco = root / "work" / "loco" / acc
    parts = [pd.read_csv(loco / f"chr{chrom}.held.tsv", sep="\t") for chrom in range(1, 23)]
    scores = pd.concat(parts, ignore_index=True)[["ENSGID", "PoPS_Score"]]
    values = pd.to_numeric(scores.PoPS_Score, errors="coerce")
    if scores.ENSGID.nunique() != 18383 or len(scores) != 18383 or not values.notna().all():
        raise RuntimeError(
            f"{acc} scored {len(scores)} rows, {scores.ENSGID.nunique()} genes, "
            f"{int(values.notna().sum())} finite"
        )
    path = score_dir / f"{acc}.tsv"
    scores.to_csv(path, sep="\t", index=False)
    for extra in loco.iterdir():
        if extra.suffix != ".tsv":
            extra.unlink()
    return path


def enqueue_chroms(acc):
    missing = []
    for chrom in range(1, 23):
        held = root / "work" / "loco" / acc / f"chr{chrom}.held.tsv"
        if not (held.exists() and held.stat().st_size > 0):
            missing.append(chrom)
    if not missing:
        write_score(acc)
        print(f"TRAIT_OK {acc} FROM_EXISTING", flush=True)
        return
    with lock:
        remaining[acc] = len(missing)
    for chrom in missing:
        chrom_q.put((acc, chrom, 0))


def magma_worker():
    while True:
        item = magma_q.get()
        if item is None:
            magma_q.task_done()
            return
        acc = item.acc
        try:
            print(f"START {acc} {item.trait}", flush=True)
            gz = download(acc, item.url)
            run_magma(acc, item.ntext, gz)
            print(f"MAGMA_READY {acc}", flush=True)
            enqueue_chroms(acc)
        except Exception as exc:
            print(f"MAGMA_FAIL {acc} {exc}", flush=True)
            with lock:
                failed.add(acc)
        finally:
            magma_q.task_done()


def pops_worker():
    global done_n
    while True:
        item = chrom_q.get()
        if item is None:
            chrom_q.task_done()
            return
        acc, chrom, tries = item
        try:
            run_chrom(acc, chrom)
        except Exception as exc:
            print(f"CHR_FAIL {acc} {chrom} try {tries} {exc}", flush=True)
            if tries < 1:
                chrom_q.put((acc, chrom, tries + 1))
            else:
                with lock:
                    failed.add(acc)
            chrom_q.task_done()
            continue
        finish = False
        with lock:
            remaining[acc] -= 1
            if remaining[acc] == 0 and acc not in failed:
                finish = True
        if finish:
            try:
                write_score(acc)
                with lock:
                    done_n += 1
                    finished = done_n
                print(
                    f"TRAIT_OK {acc} FINISHED {finished} ELAPSED_H {(time.time() - t_start) / 3600:.2f}",
                    flush=True,
                )
            except Exception as exc:
                print(f"CONCAT_FAIL {acc} {exc}", flush=True)
                with lock:
                    failed.add(acc)
        chrom_q.task_done()


def main():
    jobs = pd.read_csv(
        jobs_path, sep="\t", header=None, names=["acc", "ntext", "url", "trait"]
    )
    todo = []
    already = 0
    for row in jobs.itertuples(index=False):
        score = score_dir / f"{row.acc}.tsv"
        if score.exists() and score.stat().st_size > 0:
            already += 1
            continue
        todo.append(row)
    print(f"N_JOBS {len(jobs)} ALREADY {already} TODO {len(todo)}", flush=True)
    threads = []
    for _ in range(N_MAGMA):
        thread = threading.Thread(target=magma_worker)
        thread.start()
        threads.append(thread)
    for _ in range(N_POPS):
        thread = threading.Thread(target=pops_worker)
        thread.start()
        threads.append(thread)
    for row in todo:
        genes = root / "work" / "magma" / row.acc / f"{row.acc}.genes.out"
        if genes.exists() and genes.stat().st_size > 0:
            print(f"MAGMA_EXISTS {row.acc}", flush=True)
            enqueue_chroms(row.acc)
        else:
            magma_q.put(row)
    magma_q.join()
    for _ in range(N_MAGMA):
        magma_q.put(None)
    chrom_q.join()
    for _ in range(N_POPS):
        chrom_q.put(None)
    for thread in threads:
        thread.join()
    print(f"BATCH_DONE FAILED {sorted(failed)}", flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

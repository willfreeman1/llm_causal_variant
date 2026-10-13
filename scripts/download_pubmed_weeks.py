"""Save PubMed abstracts for the genes on the UK neighborhood lists.

Each gene symbol is matched to the human gene record in NCBI, then to the
papers NCBI already links to that gene. Title, abstract, and publication
date are stored locally. Papers dated on or after 1 October 2019 are kept
in the file but marked before_cutoff=0, because the UK Biobank fine-mapping
results were posted in October 2019 and a later search must not use them.

The run resumes. Genes with status 'done' in the database are skipped.
The API key is read from .env and is never printed.
"""

import argparse
import csv
import json
import os
import sqlite3
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
WEEKS = ROOT / "questions" / "uk_biobank_gene_lists.tsv"
OUT_DIR = ROOT / "data" / "ukb" / "pubmed"
DB_PATH = OUT_DIR / "weeks_abstracts.sqlite"
STATUS_PATH = OUT_DIR / "download_status.json"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CUTOFF = (2019, 10, 1)
MONTHS = {
    "jan": 1, "january": 1,
    "feb": 2, "february": 2,
    "mar": 3, "march": 3,
    "apr": 4, "april": 4,
    "may": 5,
    "jun": 6, "june": 6,
    "jul": 7, "july": 7,
    "aug": 8, "august": 8,
    "sep": 9, "sept": 9, "september": 9,
    "oct": 10, "october": 10,
    "nov": 11, "november": 11,
    "dec": 12, "december": 12,
}


def load_env():
    for line in (ROOT / ".env").read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.strip().startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        os.environ.setdefault(key.strip(), value.strip())
    return os.environ["NCBI_API_KEY"], os.environ["OPENALEX_MAILTO"]


def gene_symbols():
    symbols = set()
    with WEEKS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            raw = row["symbol_gene_string"].strip().strip("{}")
            for part in raw.split(","):
                symbol = part.strip().strip("{}").strip()
                if symbol:
                    symbols.add(symbol)
    return sorted(symbols)


class Limiter:
    """Start at most about 9 NCBI requests per second, across every worker."""

    def __init__(self, interval=0.11):
        self.interval = interval
        self.lock = threading.Lock()
        self.next_at = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            start_at = max(now, self.next_at)
            self.next_at = start_at + self.interval
        delay = start_at - time.monotonic()
        if delay > 0:
            time.sleep(delay)

    def penalize(self, seconds):
        with self.lock:
            self.next_at = max(self.next_at, time.monotonic() + seconds)


class Ncbi:
    def __init__(self, api_key, email, limiter):
        self.api_key = api_key
        self.email = email
        self.limiter = limiter
        self.session = requests.Session()

    def call(self, endpoint, **params):
        params["api_key"] = self.api_key
        params["tool"] = "llm_causal_variant"
        params["email"] = self.email
        last_error = "no response"
        for attempt in range(6):
            self.limiter.wait()
            try:
                response = self.session.post(BASE + endpoint, data=params, timeout=180)
            except requests.RequestException as exc:
                last_error = type(exc).__name__
                self.limiter.penalize(min(60, 2 ** attempt))
                continue
            if response.status_code == 429 or response.status_code >= 500:
                last_error = f"HTTP {response.status_code}"
                self.limiter.penalize(min(60, 2 ** attempt))
                continue
            if response.status_code != 200:
                body = response.text[:300].replace(self.api_key, "[key]")
                raise RuntimeError(f"{endpoint} HTTP {response.status_code}: {body}")
            if self.api_key in response.text[:200]:
                pass
            return response
        raise RuntimeError(f"{endpoint} failed after retries: {last_error}")

    def json(self, endpoint, **params):
        params.setdefault("retmode", "json")
        # strict=False: an error reply can contain a raw newline inside the message.
        payload = json.loads(self.call(endpoint, **params).text, strict=False)
        if isinstance(payload, dict) and payload.get("error"):
            raise RuntimeError(str(payload["error"])[:300])
        result = payload.get("esearchresult") if isinstance(payload, dict) else None
        if isinstance(result, dict) and result.get("ERROR"):
            raise RuntimeError(str(result["ERROR"])[:300])
        return payload

    def text(self, endpoint, **params):
        return self.call(endpoint, **params).text


def connect():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH, timeout=60, check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS genes (
            symbol TEXT PRIMARY KEY,
            gene_ids TEXT,
            n_pmids INTEGER,
            status TEXT,
            note TEXT
        );
        CREATE TABLE IF NOT EXISTS gene_hits (
            symbol TEXT,
            pmid INTEGER,
            PRIMARY KEY (symbol, pmid)
        );
        CREATE TABLE IF NOT EXISTS abstracts (
            pmid INTEGER PRIMARY KEY,
            pubdate TEXT,
            year INTEGER,
            month INTEGER,
            day INTEGER,
            before_cutoff INTEGER,
            title TEXT,
            abstract TEXT
        );
        """
    )
    return conn


def exact_gene_ids(client, symbol):
    found = client.json(
        "esearch.fcgi",
        db="gene",
        term=f"{symbol}[sym] AND 9606[taxid]",
        retmax=20,
    )
    ids = found.get("esearchresult", {}).get("idlist", [])
    if len(ids) <= 1:
        return ids
    summary = client.json("esummary.fcgi", db="gene", id=",".join(ids))
    result = summary.get("result", {})
    exact = [
        gene_id for gene_id in ids
        if str(result.get(gene_id, {}).get("name", "")).upper() == symbol.upper()
    ]
    return exact or ids[:1]


def linked_pmids(client, gene_id):
    linked = client.json(
        "elink.fcgi",
        dbfrom="gene",
        db="pubmed",
        id=gene_id,
        cmd="neighbor_history",
    )
    linkset = (linked.get("linksets") or [{}])[0]
    histories = linkset.get("linksetdbhistories") or []
    match = next((item for item in histories if item.get("linkname") == "gene_pubmed"), None)
    if not match:
        return []
    webenv = linkset.get("webenv")
    query_key = str(match["querykey"])
    counted = client.json(
        "esearch.fcgi",
        db="pubmed",
        term=f"#{query_key}",
        WebEnv=webenv,
        retmax=0,
        usehistory="y",
    )
    result = counted["esearchresult"]
    count = int(result["count"])
    if count <= 9999:
        return pubmed_ids(client, result["webenv"], result["querykey"], count)
    # One temporary list expires if it is reused for every year. Open a new list per span.
    pmids = []
    for start in range(1950, 2027, 5):
        pmids.extend(pmids_for_span(client, gene_id, start, min(start + 4, 2026)))
    pmids.extend(pmids_for_span(client, gene_id, 1800, 1949))
    pmids.extend(pmids_for_clause(client, gene_id, "NOT 1800:2026[pdat]"))
    return sorted(set(pmids))


def elink_keys(client, gene_id):
    linked = client.json(
        "elink.fcgi",
        dbfrom="gene",
        db="pubmed",
        id=gene_id,
        cmd="neighbor_history",
    )
    linkset = (linked.get("linksets") or [{}])[0]
    histories = linkset.get("linksetdbhistories") or []
    match = next((item for item in histories if item.get("linkname") == "gene_pubmed"), None)
    if not match:
        return None
    return linkset.get("webenv"), str(match["querykey"])


def pubmed_ids(client, webenv, query_key, expected):
    page = client.json(
        "esearch.fcgi",
        db="pubmed",
        WebEnv=webenv,
        query_key=query_key,
        retstart=0,
        retmax=9999,
    )
    ids = [int(pmid) for pmid in page["esearchresult"]["idlist"]]
    if len(ids) < expected:
        raise RuntimeError(f"PubMed returned {len(ids)} of {expected} paper ids")
    return ids


def pmids_for_clause(client, gene_id, clause):
    keys = elink_keys(client, gene_id)
    if keys is None:
        return []
    webenv, query_key = keys
    counted = client.json(
        "esearch.fcgi",
        db="pubmed",
        term=f"#{query_key} AND {clause}" if not clause.startswith("NOT ") else f"#{query_key} {clause}",
        WebEnv=webenv,
        retmax=0,
        usehistory="y",
    )
    result = counted["esearchresult"]
    count = int(result["count"])
    if count == 0:
        return []
    if count > 9999:
        return None
    return pubmed_ids(client, result["webenv"], result["querykey"], count)


def pmids_for_span(client, gene_id, start, end):
    found = pmids_for_clause(client, gene_id, f"{start}:{end}[pdat]")
    if found is not None:
        return found
    if start == end:
        first = pmids_for_clause(client, gene_id, f"{start}/01/01:{start}/06/30[pdat]")
        second = pmids_for_clause(client, gene_id, f"{start}/07/01:{start}/12/31[pdat]")
        if first is None or second is None:
            raise RuntimeError(f"{start} still has more than 9999 papers in one half year")
        return first + second
    mid = (start + end) // 2
    return pmids_for_span(client, gene_id, start, mid) + pmids_for_span(client, gene_id, mid + 1, end)


def parse_date(node):
    if node is None:
        return None
    year = node.findtext("Year")
    month = node.findtext("Month")
    day = node.findtext("Day")
    if year and year.isdigit():
        month_num = None
        if month:
            month_num = MONTHS.get(month.lower())
            if month_num is None and month.isdigit():
                month_num = int(month)
        day_num = int(day) if day and day.isdigit() else None
        return int(year), month_num, day_num
    medline = node.findtext("MedlineDate") or ""
    digits = ""
    for char in medline:
        if char.isdigit():
            digits += char
            if len(digits) == 4:
                break
        elif digits:
            break
    if len(digits) != 4:
        return None
    year_num = int(digits)
    month_num = None
    for name, value in MONTHS.items():
        if name in medline.lower():
            month_num = value
            break
    return year_num, month_num, None


def before_cutoff(parsed):
    """True when the publication date is known to be before October 2019.

    A year with no month is kept only when that year is 2018 or earlier.
    Year-only 2019 is excluded, because that paper may be from October or later.
    """
    if parsed is None:
        return 0
    year, month, _day = parsed
    if year < 2019:
        return 1
    if year > 2019 or month is None or month >= 10:
        return 0
    return 1


def parse_articles(xml_text):
    root = ET.fromstring(xml_text.encode("utf-8") if isinstance(xml_text, str) else xml_text)
    rows = []
    for article in root.iter("PubmedArticle"):
        medline = article.find("MedlineCitation")
        if medline is None or medline.findtext("PMID") is None:
            continue
        pmid = int(medline.findtext("PMID"))
        art = medline.find("Article")
        title = ""
        abstract = ""
        parsed = None
        if art is not None:
            title_node = art.find("ArticleTitle")
            if title_node is not None:
                title = "".join(title_node.itertext()).strip()
            abstract_node = art.find("Abstract")
            parts = []
            if abstract_node is not None:
                for node in abstract_node.findall("AbstractText"):
                    label = node.attrib.get("Label")
                    text = "".join(node.itertext()).strip()
                    if text:
                        parts.append(f"{label}: {text}" if label else text)
            abstract = "\n".join(parts)
            parsed = parse_date(art.find("Journal/JournalIssue/PubDate"))
            if parsed is None:
                parsed = parse_date(art.find("ArticleDate"))
        if parsed is None:
            history = article.find("PubmedData/History")
            if history is not None:
                for stamp in history.findall("PubMedPubDate"):
                    if stamp.attrib.get("PubStatus") == "pubmed":
                        parsed = parse_date(stamp)
                        break
        pubdate = ""
        year = month = day = None
        if parsed is not None:
            year, month, day = parsed
            pubdate = "-".join(str(part) for part in (year, month, day) if part)
        rows.append((pmid, pubdate, year, month, day, before_cutoff(parsed), title, abstract))
    return rows


def store_abstracts(conn, client, pmids, known, locks):
    with locks["known"]:
        missing = [pmid for pmid in pmids if pmid not in known]
        known.update(missing)
    stored = 0
    kept = 0
    finished = 0
    try:
        for start in range(0, len(missing), 200):
            batch = missing[start:start + 200]
            xml_text = client.text(
                "efetch.fcgi",
                db="pubmed",
                id=",".join(str(pmid) for pmid in batch),
                rettype="xml",
                retmode="xml",
            )
            if "<ERROR" in xml_text[:500]:
                raise RuntimeError(xml_text[:300].replace(client.api_key, "[key]"))
            rows = parse_articles(xml_text)
            with locks["db"]:
                conn.executemany(
                    """
                    INSERT OR REPLACE INTO abstracts
                    (pmid, pubdate, year, month, day, before_cutoff, title, abstract)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    rows,
                )
                conn.commit()
            stored += len(rows)
            kept += sum(row[5] for row in rows)
            finished = start + len(batch)
    except Exception:
        with locks["known"]:
            for pmid in missing[finished:]:
                known.discard(pmid)
        raise
    return stored, kept


def write_status(conn, total, status_path=STATUS_PATH):
    genes_done = conn.execute("SELECT COUNT(*) FROM genes WHERE status = 'done'").fetchone()[0]
    abstracts = conn.execute("SELECT COUNT(*) FROM abstracts").fetchone()[0]
    kept = conn.execute("SELECT COUNT(*) FROM abstracts WHERE before_cutoff = 1").fetchone()[0]
    payload = {
        "genes_total": total,
        "genes_done": genes_done,
        "abstracts": abstracts,
        "abstracts_before_cutoff": kept,
    }
    status_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--symbols", default="")
    parser.add_argument("--symbols-file", default="")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    api_key, email = load_env()
    limiter = Limiter()
    if args.symbols_file:
        symbols = [
            line.strip()
            for line in Path(args.symbols_file).read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    else:
        symbols = gene_symbols()
        if args.symbols:
            wanted = {part.strip() for part in args.symbols.split(",") if part.strip()}
            symbols = [symbol for symbol in symbols if symbol in wanted]
    if args.limit:
        symbols = symbols[:args.limit]
    status_path = OUT_DIR / "exam_download_status.json" if args.symbols_file else STATUS_PATH
    conn = connect()
    done = {
        row[0]
        for row in conn.execute("SELECT symbol FROM genes WHERE status = 'done'")
    }
    known = {row[0] for row in conn.execute("SELECT pmid FROM abstracts")}
    locks = {"db": threading.Lock(), "known": threading.Lock()}
    pending = [symbol for symbol in symbols if symbol not in done]
    print(
        f"symbols {len(symbols)} already_done {len(symbols) - len(pending)} abstracts {len(known)} workers {args.workers}",
        flush=True,
    )
    started = time.time()
    finished_count = 0

    def process(symbol):
        client = Ncbi(api_key, email, limiter)
        try:
            gene_ids = exact_gene_ids(client, symbol)
            if not gene_ids:
                with locks["db"]:
                    conn.execute(
                        "INSERT OR REPLACE INTO genes (symbol, gene_ids, n_pmids, status, note) VALUES (?, '', 0, 'done', 'no gene record')",
                        (symbol,),
                    )
                    conn.commit()
                print(f"{symbol} no gene record", flush=True)
                return
            pmids = []
            for gene_id in gene_ids:
                pmids.extend(linked_pmids(client, gene_id))
            pmids = sorted(set(pmids))
            with locks["db"]:
                conn.executemany(
                    "INSERT OR IGNORE INTO gene_hits (symbol, pmid) VALUES (?, ?)",
                    [(symbol, pmid) for pmid in pmids],
                )
                conn.commit()
            stored, kept = store_abstracts(conn, client, pmids, known, locks)
            with locks["db"]:
                conn.execute(
                    "INSERT OR REPLACE INTO genes (symbol, gene_ids, n_pmids, status, note) VALUES (?, ?, ?, 'done', '')",
                    (symbol, ",".join(gene_ids), len(pmids)),
                )
                conn.commit()
            print(
                f"{symbol} genes={','.join(gene_ids)} pmids={len(pmids)} new={stored} kept_in_batch={kept}",
                flush=True,
            )
        except Exception as exc:
            message = str(exc).replace(api_key, "[key]")[:400]
            with locks["db"]:
                conn.execute(
                    "INSERT OR REPLACE INTO genes (symbol, gene_ids, n_pmids, status, note) VALUES (?, '', 0, 'error', ?)",
                    (symbol, message),
                )
                conn.commit()
            print(f"{symbol} ERROR {message}", flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(process, symbol) for symbol in pending]
        for future in as_completed(futures):
            future.result()
            finished_count += 1
            if finished_count % 25 == 0:
                with locks["db"]:
                    status = write_status(conn, len(symbols), status_path)
                print(
                    f"progress genes_done={status['genes_done']} abstracts={status['abstracts']} before_cutoff={status['abstracts_before_cutoff']} elapsed_min={((time.time()-started)/60):.1f}",
                    flush=True,
                )
    status_path = OUT_DIR / "exam_download_status.json" if args.symbols_file else STATUS_PATH
    status = write_status(conn, len(symbols), status_path)
    print("finished", json.dumps(status), flush=True)


if __name__ == "__main__":
    main()

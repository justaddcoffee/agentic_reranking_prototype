"""Fetch PubMed pubdate / epubdate for all PMIDs in pmids.txt using esummary,
batched 200 IDs per call, rate-limited to ~3 req/sec.
Uses curl (with -k) because the local Python urllib hits SSL inspection issues.
Writes pubdates.json:
  {"<pmid>": {"pubdate": "...", "epubdate": "...", "sortpubdate": "...", "earliest_iso": "YYYY-MM-DD"}, ...}
"""
import json
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).parent
PMIDS = (ROOT / "pmids.txt").read_text().split()
OUT = ROOT / "pubdates.json"

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

MONTHS = {"Jan":1,"Feb":2,"Mar":3,"Apr":4,"May":5,"Jun":6,"Jul":7,
          "Aug":8,"Sep":9,"Oct":10,"Nov":11,"Dec":12}

def parse_pubdate(s: str):
    if not s:
        return None
    s = s.strip()
    m = re.match(r"^(\d{4})\b", s)
    if not m:
        return None
    year = int(m.group(1))
    rest = s[m.end():].strip()
    season_month = {"Spring":3,"Summer":6,"Fall":9,"Autumn":9,"Winter":12}
    month = 1
    day = 1
    if rest:
        mm = re.match(r"([A-Za-z]+)", rest)
        if mm:
            tok = mm.group(1)
            if tok in MONTHS:
                month = MONTHS[tok]
                rest2 = rest[mm.end():].strip()
                # day or month range
                if rest2.startswith("-"):
                    pass  # range like Mar-Apr -> earliest is the first month
                else:
                    dm2 = re.match(r"(\d+)", rest2)
                    if dm2:
                        try:
                            day = int(dm2.group(1))
                        except Exception:
                            day = 1
            elif tok in season_month:
                month = season_month[tok]
    try:
        return f"{year:04d}-{month:02d}-{day:02d}"
    except Exception:
        return None

def fetch_batch(ids):
    url = f"{ESUMMARY}?db=pubmed&id={','.join(ids)}&retmode=json"
    cp = subprocess.run(
        ["curl", "-sk", "--max-time", "60",
         "-A", "phenopacket-store-bench/0.1 (justinreese@lbl.gov)",
         url],
        capture_output=True, check=True,
    )
    return json.loads(cp.stdout.decode())

result = {}
B = 200
for i in range(0, len(PMIDS), B):
    chunk = PMIDS[i:i+B]
    last_err = None
    for attempt in range(5):
        try:
            data = fetch_batch(chunk)
            break
        except Exception as e:
            last_err = e
            time.sleep(2 + attempt)
    else:
        print(f"FAILED batch {i}: {last_err}")
        continue
    rs = data.get("result", {})
    uids = rs.get("uids", [])
    for uid in uids:
        rec = rs.get(uid, {}) or {}
        if rec.get("error"):
            result[uid] = {"error": rec.get("error"), "earliest_iso": None}
            continue
        pubdate = rec.get("pubdate", "")
        epubdate = rec.get("epubdate", "")
        sortpubdate = rec.get("sortpubdate", "")
        d_pub = parse_pubdate(pubdate)
        d_epub = parse_pubdate(epubdate) if epubdate else None
        d_sort = None
        if sortpubdate:
            ms = re.match(r"(\d{4})/(\d{2})/(\d{2})", sortpubdate)
            if ms:
                d_sort = f"{ms.group(1)}-{ms.group(2)}-{ms.group(3)}"
        candidates = [x for x in (d_pub, d_epub, d_sort) if x]
        earliest = min(candidates) if candidates else None
        result[uid] = {
            "pubdate": pubdate,
            "epubdate": epubdate,
            "sortpubdate": sortpubdate,
            "earliest_iso": earliest,
        }
    print(f"batch {i//B + 1}/{(len(PMIDS)+B-1)//B} done; got {len(uids)} ids; total resolved {len(result)}")
    time.sleep(0.4)

OUT.write_text(json.dumps(result, indent=2))

unresolved = [p for p in PMIDS if p not in result]
print(f"PMIDs requested: {len(PMIDS)}")
print(f"PMIDs resolved:  {len(result)}")
print(f"Unresolved:      {len(unresolved)}")
if unresolved[:10]:
    print("Sample unresolved:", unresolved[:10])

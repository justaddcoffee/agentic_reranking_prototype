"""Combine phenopackets.csv + pubdates.json: count phenopackets whose earliest
associated PubMed publication date is >= 2026-01-01. Report per-cohort breakdown
of the post-2026-01 subset."""
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
THRESHOLD = "2026-01-01"

pubdates = json.loads((ROOT / "pubdates.json").read_text())

with (ROOT / "phenopackets.csv").open() as f:
    rows = list(csv.DictReader(f))

n_total = len(rows)
n_with_pmid = 0
n_with_resolved_date = 0
n_post_2026 = 0
post_rows = []
unresolved_examples = []
no_date_examples = []
errors = 0

for r in rows:
    pmids = [p for p in r["pmids"].split(";") if p]
    if not pmids:
        continue
    n_with_pmid += 1
    dates = []
    any_resolved = False
    for p in pmids:
        rec = pubdates.get(p)
        if not rec:
            continue
        if rec.get("error"):
            errors += 1
            continue
        d = rec.get("earliest_iso")
        if d:
            dates.append(d)
            any_resolved = True
    if any_resolved:
        n_with_resolved_date += 1
    else:
        if len(no_date_examples) < 5:
            no_date_examples.append((r["file"], pmids))
        continue
    earliest = min(dates)
    r["_earliest"] = earliest
    if earliest >= THRESHOLD:
        n_post_2026 += 1
        post_rows.append(r)

print(f"Total phenopackets:                       {n_total}")
print(f"With at least one PMID listed:            {n_with_pmid}")
print(f"With at least one resolvable PubMed date: {n_with_resolved_date}")
print(f"Resolution errors (PMIDs returning err):  {errors}")
print(f"With earliest pubdate >= {THRESHOLD}:    {n_post_2026}")
print()
if no_date_examples:
    print("Examples with PMIDs but no resolvable pubdate:")
    for ex in no_date_examples:
        print("  ", ex)
    print()

cohort_ct = Counter(r["cohort"] for r in post_rows)
disease_ct = Counter(r["disease_label"] for r in post_rows if r["disease_label"])
print("Top 10 cohorts in post-2026-01-01 subset:")
for c, n in cohort_ct.most_common(10):
    print(f"  {n:5d}  {c}")
print()
print("Top 10 diseases in post-2026-01-01 subset:")
for c, n in disease_ct.most_common(10):
    print(f"  {n:5d}  {c}")
print()

# Distribution by year of earliest pubdate (sanity)
year_ct = Counter()
for r in rows:
    if "_earliest" in r:
        year_ct[r["_earliest"][:4]] += 1
print("Year distribution of earliest pubdate (top 20 years):")
for y, n in sorted(year_ct.items()):
    print(f"  {y}: {n}")

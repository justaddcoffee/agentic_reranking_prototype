"""Walk every phenopacket JSON in extracted/ and collect:
  - file path
  - cohort (immediate parent dir under release dir)
  - disease label(s)
  - PMID list from metaData.externalReferences (id starts with 'PMID:')
  - non-PubMed external reference ids (for caveat reporting)
Writes CSV: phenopackets.csv
"""
import csv
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).parent / "extracted" / "0.1.26"
OUT = Path(__file__).parent / "phenopackets.csv"

PMID_RE = re.compile(r"^PMID[:_]?(\d+)$", re.IGNORECASE)

rows = []
n_total = 0
n_no_extref = 0
n_no_pmid = 0
non_pmid_refs = []

for cohort_dir in sorted(ROOT.iterdir()):
    if not cohort_dir.is_dir():
        continue
    cohort = cohort_dir.name
    for jp in cohort_dir.rglob("*.json"):
        n_total += 1
        try:
            d = json.loads(jp.read_text())
        except Exception as e:
            print(f"parse error {jp}: {e}")
            continue
        meta = d.get("metaData", {}) or {}
        ext_refs = meta.get("externalReferences", []) or []
        pmids = []
        non_pmid = []
        for er in ext_refs:
            rid = (er.get("id") or "").strip()
            m = PMID_RE.match(rid)
            if m:
                pmids.append(m.group(1))
            elif rid:
                non_pmid.append(rid)
        if not ext_refs:
            n_no_extref += 1
        if not pmids:
            n_no_pmid += 1
            for x in non_pmid:
                non_pmid_refs.append(x)
        diseases = d.get("diseases", []) or []
        disease_label = ""
        disease_id = ""
        if diseases:
            t = diseases[0].get("term", {}) or {}
            disease_label = t.get("label", "")
            disease_id = t.get("id", "")
        # Some phenopackets put diagnosis only in interpretations
        if not disease_label:
            interps = d.get("interpretations", []) or []
            if interps:
                diag = interps[0].get("diagnosis", {}) or {}
                t = diag.get("disease", {}) or {}
                disease_label = t.get("label", "")
                disease_id = t.get("id", "")
        rows.append({
            "file": str(jp.relative_to(ROOT)),
            "cohort": cohort,
            "disease_id": disease_id,
            "disease_label": disease_label,
            "pmids": ";".join(pmids),
            "non_pmid_refs": ";".join(non_pmid),
        })

with OUT.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["file", "cohort", "disease_id", "disease_label", "pmids", "non_pmid_refs"])
    w.writeheader()
    w.writerows(rows)

print(f"total phenopackets: {n_total}")
print(f"no externalReferences at all: {n_no_extref}")
print(f"no PMID in externalReferences: {n_no_pmid}")
print(f"distinct non-PubMed ref ids: {len(set(non_pmid_refs))}")
if non_pmid_refs:
    from collections import Counter
    print("Sample non-PubMed:", Counter(non_pmid_refs).most_common(10))

# Distinct PMIDs (for batched lookup)
pmid_set = set()
for r in rows:
    if r["pmids"]:
        pmid_set.update(r["pmids"].split(";"))
print(f"distinct PMIDs: {len(pmid_set)}")
(Path(__file__).parent / "pmids.txt").write_text("\n".join(sorted(pmid_set, key=int)) + "\n")

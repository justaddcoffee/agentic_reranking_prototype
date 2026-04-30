"""Build a per-run cohort manifest from a `CohortSource`.

The manifest is the single source of truth for which cases participate in
this experimental run. Every downstream step (variant spiking, perturbation,
Exomiser, agent, scoring) takes ``--run-id`` and resolves all paths from
``data/derived/runs/<run_id>/``.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO

import yaml

from agentic_reranking.sources.base import CaseBundle, Phenopacket
from agentic_reranking.sources.pps import PhenopacketStoreSource

# Hardcoded for the pilot. Promote to a config table or derive from the
# 2508 OMIM snapshot once Task 9 lands. All four entries below are gene-disease
# links first described in 2025 (i.e., post-OMIM-2508). Netherton (OMIM:256500,
# SPINK5) is the established outlier in the post-2025-09 cohort.
_NOVEL_DISEASE_IDS = {
    "OMIM:621435",  # XP complementation group J (GTF2H4)
    "OMIM:621452",  # Basal ganglia calcification 11, AR (RRP12)
    "OMIM:621386",  # Valence-Farazi cerebellar ataxia syndrome (SKOR2)
    "OMIM:621393",  # Neurodev disorder with brain/craniofacial abnormalities (SNAPIN)
}

MANIFEST_COLUMNS = [
    "case_id",
    "source_id",
    "source_version",
    "source_pmid",
    "source_pubdate",
    "sex",
    "age_years",
    "ancestry",
    "hpo_terms",
    "disease_id",
    "disease_label",
    "truth_gene_symbol",
    "truth_gene_id",
    "truth_hgvs",
    "novelty_class",
]


@dataclass(frozen=True)
class ManifestRow:
    case_id: str
    source_id: str
    source_version: str
    source_pmid: str
    source_pubdate: str
    sex: str
    age_years: str  # blank if unknown; written as int-like or empty string
    ancestry: str
    hpo_terms: str  # ";"-joined HPO IDs, excluding `excluded=True` features
    disease_id: str
    disease_label: str
    truth_gene_symbol: str
    truth_gene_id: str
    truth_hgvs: str
    novelty_class: str  # "novel" | "established"

    def as_dict(self) -> dict[str, str]:
        return {col: getattr(self, col) for col in MANIFEST_COLUMNS}


def _age_years_from_iso(duration: str | None) -> str:
    """Parse ``P5Y`` / ``P5Y6M`` / ``P0Y6M`` → ``"5"`` / ``"5"`` / ``"0"``.

    Returns blank when no parseable years component is present.
    """
    if not duration:
        return ""
    m = re.match(r"^P(\d+)Y", duration)
    return m.group(1) if m else ""


def _sex_from_phenopacket(value: str | None) -> str:
    if not value:
        return ""
    return value.lower() if value in {"FEMALE", "MALE", "OTHER_SEX"} else ""


def _truth_from_phenopacket(pkt: Phenopacket) -> tuple[str, str, str]:
    """Extract (gene_symbol, gene_id, hgvs_c) from the first genomic interpretation."""
    interpretations = pkt.get("interpretations") or []
    if not interpretations:
        return ("", "", "")
    diagnosis = interpretations[0].get("diagnosis") or {}
    gis = diagnosis.get("genomicInterpretations") or []
    if not gis:
        return ("", "", "")
    descriptor = gis[0].get("variantInterpretation", {}).get("variationDescriptor", {})
    gene = descriptor.get("geneContext", {})
    hgvs = ""
    for expr in descriptor.get("expressions", []):
        if expr.get("syntax") == "hgvs.c":
            hgvs = expr.get("value", "")
            break
    return (gene.get("symbol", ""), gene.get("valueId", ""), hgvs)


def _hpo_terms_from_phenopacket(pkt: Phenopacket) -> str:
    """Semicolon-joined HPO IDs, *excluding* phenotypes flagged ``excluded=True``."""
    return ";".join(
        pf["type"]["id"]
        for pf in pkt.get("phenotypicFeatures", [])
        if not pf.get("excluded", False) and "type" in pf and "id" in pf["type"]
    )


def _row_from_bundle(
    bundle: CaseBundle,
    *,
    source_pmid: str,
    source_pubdate: str,
    cohort_disease_id: str,
    cohort_disease_label: str,
) -> ManifestRow:
    pkt = bundle.phenopacket
    subject = pkt.get("subject", {})
    age_iso = (subject.get("timeAtLastEncounter", {}).get("age") or {}).get("iso8601duration")
    truth_symbol, truth_gene_id, truth_hgvs = _truth_from_phenopacket(pkt)
    novelty = "novel" if cohort_disease_id in _NOVEL_DISEASE_IDS else "established"
    return ManifestRow(
        case_id=bundle.case_id,
        source_id=bundle.policy.source_id,
        source_version=bundle.policy.source_version,
        source_pmid=source_pmid,
        source_pubdate=source_pubdate,
        sex=_sex_from_phenopacket(subject.get("sex")),
        age_years=_age_years_from_iso(age_iso),
        ancestry="",  # PPS phenopackets rarely carry ancestry
        hpo_terms=_hpo_terms_from_phenopacket(pkt),
        disease_id=cohort_disease_id,
        disease_label=cohort_disease_label,
        truth_gene_symbol=truth_symbol,
        truth_gene_id=truth_gene_id,
        truth_hgvs=truth_hgvs,
        novelty_class=novelty,
    )


def build_pps_manifest(
    src: PhenopacketStoreSource,
    *,
    min_pub_date: str | None = "2025-09-01",
) -> list[ManifestRow]:
    """Build manifest rows for a PPS source, optionally filtering by publication date.

    PPS-specific because publication date is a PPS concept; GEL/EHR sources
    will have their own builders that filter by encounter or interpretation date.
    """
    rows: list[ManifestRow] = []
    for case_id in src.iter_case_ids():
        pubdate = src.earliest_pub_date(case_id)
        if min_pub_date is not None and (pubdate is None or pubdate < min_pub_date):
            continue
        meta = src.get_row(case_id)
        bundle = src.materialize(case_id)
        rows.append(
            _row_from_bundle(
                bundle,
                source_pmid=meta.pmids[0] if meta.pmids else "",
                source_pubdate=pubdate or "",
                cohort_disease_id=meta.disease_id,
                cohort_disease_label=meta.disease_label,
            )
        )
    return rows


def write_manifest(rows: list[ManifestRow], out: TextIO) -> None:
    writer = csv.DictWriter(
        out, fieldnames=MANIFEST_COLUMNS, dialect="excel-tab", lineterminator="\n"
    )
    writer.writeheader()
    for r in rows:
        writer.writerow(r.as_dict())


def write_run_yaml(
    out: TextIO,
    *,
    source_id: str,
    source_version: str,
    run_id: str,
    filters: dict[str, str],
    n_cases: int,
) -> None:
    yaml.safe_dump(
        {
            "run_id": run_id,
            "created_at": datetime.now(UTC).isoformat(),
            "source": {"id": source_id, "version": source_version},
            "filters": filters,
            "n_cases": n_cases,
        },
        out,
        sort_keys=False,
    )


def _run_dir(run_id: str) -> Path:
    return Path("data/derived/runs") / run_id


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the cohort manifest for a run.")
    parser.add_argument("--source", default="pps", choices=["pps"])
    parser.add_argument("--source-version", required=True)
    parser.add_argument(
        "--source-root",
        default=None,
        help="Override source root path (default: data/raw/<source>/<source-version>)",
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument(
        "--min-pub-date",
        default="2025-09-01",
        help="Filter to PMIDs first published on or after this ISO date (PPS only)",
    )
    args = parser.parse_args()

    if args.source != "pps":
        raise SystemExit(f"unsupported source: {args.source}")

    root = (
        Path(args.source_root) if args.source_root else Path("data/raw/pps") / args.source_version
    )
    src = PhenopacketStoreSource(root)
    rows = build_pps_manifest(src, min_pub_date=args.min_pub_date)

    run_dir = _run_dir(args.run_id)
    run_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "manifest.tsv").open("w") as fh:
        write_manifest(rows, fh)
    with (run_dir / "run.yaml").open("w") as fh:
        write_run_yaml(
            fh,
            source_id="pps",
            source_version=args.source_version,
            run_id=args.run_id,
            filters={"min_pub_date": args.min_pub_date},
            n_cases=len(rows),
        )
    print(f"wrote {len(rows)} rows to {run_dir / 'manifest.tsv'}")


if __name__ == "__main__":
    main()

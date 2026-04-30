"""Unit tests for the run-manifest builder."""

from __future__ import annotations

import io
from pathlib import Path

import pytest
import yaml

from agentic_reranking.cohort.manifest import (
    MANIFEST_COLUMNS,
    build_pps_manifest,
    write_manifest,
    write_run_yaml,
)
from agentic_reranking.sources.pps import PhenopacketStoreSource

FIXTURE_ROOT = Path(__file__).parent.parent / "fixtures" / "pps" / "v0.1.26"
REAL_ROOT = Path("data/raw/pps/v0.1.26")


@pytest.fixture
def src() -> PhenopacketStoreSource:
    return PhenopacketStoreSource(FIXTURE_ROOT)


def test_no_filter_returns_all_fixture_cases(src: PhenopacketStoreSource) -> None:
    rows = build_pps_manifest(src, min_pub_date=None)
    assert len(rows) == 4


def test_post_2025_09_filter_drops_pre_cutoff_cases(src: PhenopacketStoreSource) -> None:
    rows = build_pps_manifest(src, min_pub_date="2025-09-01")
    case_ids = sorted(r.case_id for r in rows)
    assert case_ids == [
        "GTF2H4__PMID_40924475_XP140BR",
        "RRP12__PMID_41059649_Case_C_II_2",
    ]


def test_extracts_demographics_and_truth(src: PhenopacketStoreSource) -> None:
    rows = build_pps_manifest(src, min_pub_date="2025-09-01")
    by_case = {r.case_id: r for r in rows}
    xp = by_case["GTF2H4__PMID_40924475_XP140BR"]
    assert xp.sex == "female"
    assert xp.age_years == "6"
    assert xp.truth_gene_symbol == "GTF2H4"
    assert xp.truth_gene_id == "HGNC:4658"
    assert xp.truth_hgvs == "NM_001517.5:c.138-1G>A"
    assert xp.disease_id == "OMIM:621435"
    assert xp.source_pmid == "40924475"
    assert xp.source_pubdate == "2025-09-09"


def test_excluded_phenotypes_are_dropped(src: PhenopacketStoreSource) -> None:
    """The XP140BR case has a phenotypic feature with `excluded: true` (Anemia).
    It must not appear in the manifest's hpo_terms column."""
    rows = build_pps_manifest(src, min_pub_date="2025-09-01")
    xp = next(r for r in rows if r.case_id == "GTF2H4__PMID_40924475_XP140BR")
    hpo_set = set(xp.hpo_terms.split(";"))
    assert "HP:0001903" not in hpo_set  # Anemia, excluded in this case


def test_novelty_class_uses_disease_id_table(src: PhenopacketStoreSource) -> None:
    rows = build_pps_manifest(src, min_pub_date=None)
    by_case = {r.case_id: r for r in rows}
    assert by_case["GTF2H4__PMID_40924475_XP140BR"].novelty_class == "novel"
    assert by_case["RRP12__PMID_41059649_Case_C_II_2"].novelty_class == "novel"
    assert by_case["SPINK5__PMID_24506793_Case_1"].novelty_class == "established"
    assert by_case["11q_terminal_deletion__PMID_15266616_35"].novelty_class == "established"


def test_write_manifest_emits_tsv_with_correct_columns(src: PhenopacketStoreSource) -> None:
    rows = build_pps_manifest(src, min_pub_date="2025-09-01")
    buf = io.StringIO()
    write_manifest(rows, buf)
    out = buf.getvalue()
    header = out.splitlines()[0].split("\t")
    assert header == MANIFEST_COLUMNS
    assert len(out.splitlines()) == 1 + len(rows)


def test_write_run_yaml_captures_provenance() -> None:
    buf = io.StringIO()
    write_run_yaml(
        buf,
        source_id="pps",
        source_version="v0.1.26",
        run_id="test_run",
        filters={"min_pub_date": "2025-09-01"},
        n_cases=2,
    )
    parsed = yaml.safe_load(buf.getvalue())
    assert parsed["run_id"] == "test_run"
    assert parsed["source"] == {"id": "pps", "version": "v0.1.26"}
    assert parsed["filters"] == {"min_pub_date": "2025-09-01"}
    assert parsed["n_cases"] == 2
    assert "created_at" in parsed


# --- Integration test against the real PPS v0.1.26 release -----------------


@pytest.mark.skipif(
    not (REAL_ROOT / "phenopackets.csv").exists(),
    reason="real PPS data not present",
)
def test_real_cohort_yields_22_cases() -> None:
    src = PhenopacketStoreSource(REAL_ROOT)
    rows = build_pps_manifest(src, min_pub_date="2025-09-01")
    assert len(rows) == 22
    # 5 unique diseases as established earlier
    assert len({r.disease_id for r in rows}) == 5

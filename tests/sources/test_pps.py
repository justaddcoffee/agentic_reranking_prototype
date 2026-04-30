"""Unit tests for `PhenopacketStoreSource` against the vendored fixture."""

from __future__ import annotations

from pathlib import Path

import pytest

from agentic_reranking.sources.pps import PhenopacketStoreSource

FIXTURE_ROOT = Path(__file__).parent.parent / "fixtures" / "pps" / "v0.1.26"
REAL_ROOT = Path("data/raw/pps/v0.1.26")


@pytest.fixture
def src() -> PhenopacketStoreSource:
    return PhenopacketStoreSource(FIXTURE_ROOT)


def test_iter_case_ids_yields_all_fixture_cases(src: PhenopacketStoreSource) -> None:
    case_ids = sorted(src.iter_case_ids())
    assert case_ids == [
        "11q_terminal_deletion__PMID_15266616_35",
        "GTF2H4__PMID_40924475_XP140BR",
        "RRP12__PMID_41059649_Case_C_II_2",
        "SPINK5__PMID_24506793_Case_1",
    ]


def test_materialize_returns_case_bundle_with_correct_policy(
    src: PhenopacketStoreSource,
) -> None:
    bundle = src.materialize("GTF2H4__PMID_40924475_XP140BR")
    assert bundle.case_id == "GTF2H4__PMID_40924475_XP140BR"
    assert bundle.policy.source_id == "pps"
    assert bundle.policy.has_publication_date is True
    assert bundle.policy.has_real_vcf is False
    assert bundle.policy.public_literature_source is True
    assert bundle.policy.has_phi is False


def test_materialize_loads_real_phenopacket_clinical_content(
    src: PhenopacketStoreSource,
) -> None:
    bundle = src.materialize("GTF2H4__PMID_40924475_XP140BR")
    pkt = bundle.phenopacket
    assert pkt["subject"]["id"] == "XP140BR"
    assert pkt["subject"]["sex"] == "FEMALE"
    hpo_ids = [pf["type"]["id"] for pf in pkt["phenotypicFeatures"]]
    assert all(h.startswith("HP:") for h in hpo_ids)
    truth_gene = pkt["interpretations"][0]["diagnosis"]["genomicInterpretations"][0][
        "variantInterpretation"
    ]["variationDescriptor"]["geneContext"]["symbol"]
    assert truth_gene == "GTF2H4"


def test_materialize_unknown_case_id_raises(src: PhenopacketStoreSource) -> None:
    with pytest.raises(KeyError):
        src.materialize("does_not_exist")


def test_earliest_pub_date_returns_iso_date(src: PhenopacketStoreSource) -> None:
    assert src.earliest_pub_date("GTF2H4__PMID_40924475_XP140BR") == "2025-09-09"
    assert src.earliest_pub_date("11q_terminal_deletion__PMID_15266616_35") == "2004-08-15"


# --- Integration tests against the real PPS v0.1.26 release ---------------


@pytest.mark.skipif(
    not (REAL_ROOT / "phenopackets.csv").exists(),
    reason="real PPS data not present (run scripts/fetch_pps.sh)",
)
class TestRealPpsRelease:
    def test_iter_yields_9588_cases(self) -> None:
        src = PhenopacketStoreSource(REAL_ROOT)
        assert sum(1 for _ in src.iter_case_ids()) == 9588

    def test_known_case_materializes(self) -> None:
        src = PhenopacketStoreSource(REAL_ROOT)
        bundle = src.materialize("GTF2H4__PMID_40924475_XP140BR")
        assert bundle.phenopacket["subject"]["id"] == "XP140BR"

    def test_post_2025_09_subset_has_22_cases(self) -> None:
        src = PhenopacketStoreSource(REAL_ROOT)
        post = [
            cid
            for cid in src.iter_case_ids()
            if (d := src.earliest_pub_date(cid)) is not None and d >= "2025-09-01"
        ]
        assert len(post) == 22

from dataclasses import FrozenInstanceError

import pytest

from agentic_reranking.sources.base import CaseBundle, SourcePolicy


def _make_policy() -> SourcePolicy:
    return SourcePolicy(
        source_id="pps",
        source_version="v0.1.26",
        has_real_vcf=False,
        has_publication_date=True,
        has_phi=False,
        restricted_access=False,
        consent_bounded=False,
        public_literature_source=True,
    )


def test_source_policy_is_frozen() -> None:
    p = _make_policy()
    with pytest.raises(FrozenInstanceError):
        p.has_phi = True  # type: ignore[misc]


def test_source_policy_extra_defaults_to_empty_dict() -> None:
    p = _make_policy()
    assert p.extra == {}


def test_source_policy_extra_can_carry_source_specific_tags() -> None:
    p = SourcePolicy(
        source_id="pps",
        source_version="v0.1.26",
        has_real_vcf=False,
        has_publication_date=True,
        has_phi=False,
        restricted_access=False,
        consent_bounded=False,
        public_literature_source=True,
        extra={"phenopacket_store_release": "v0.1.26"},
    )
    assert p.extra["phenopacket_store_release"] == "v0.1.26"


def test_case_bundle_is_frozen() -> None:
    bundle = CaseBundle(case_id="c1", phenopacket={"id": "c1"}, policy=_make_policy())
    with pytest.raises(FrozenInstanceError):
        bundle.case_id = "c2"  # type: ignore[misc]


def test_case_bundle_carries_phenopacket_and_policy() -> None:
    pkt = {"id": "patient1", "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}]}
    bundle = CaseBundle(case_id="case1", phenopacket=pkt, policy=_make_policy())
    assert bundle.phenopacket["id"] == "patient1"
    assert bundle.policy.source_id == "pps"

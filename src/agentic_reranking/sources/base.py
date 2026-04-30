"""Cohort source contract.

A `CohortSource` yields `CaseBundle`s, each pairing a GA4GH Phenopacket
(canonical clinical content, Exomiser-aligned) with a `SourcePolicy` sidecar
(operational metadata: capability flags, consent status, source version).

When a second cohort source materializes (GEL, dark EHR, ...), the Protocol
will grow `leak_controls(cutoff)` and `execution_policy()` methods. Today's
PPS pilot uses hardcoded controls; the Protocol's stub fields document the
intended evolution.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any, Protocol

# A parsed GA4GH Phenopacket v2 JSON document. We don't model the full schema
# here — Exomiser and the rest of the pipeline read the fields they need.
Phenopacket = dict[str, Any]


@dataclass(frozen=True)
class SourcePolicy:
    """Operational metadata for a cohort source.

    Capability flags drive pipeline branches (e.g., variant spiking only
    runs when `has_real_vcf=False`) and gate per-source leak controls.
    """

    source_id: str
    source_version: str
    has_real_vcf: bool
    has_publication_date: bool
    has_phi: bool
    restricted_access: bool
    consent_bounded: bool
    public_literature_source: bool
    extra: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class CaseBundle:
    """One case from a cohort source.

    `phenopacket` is the canonical clinical content (parsed GA4GH JSON,
    consumed directly by Exomiser). `policy` carries operational metadata
    that doesn't fit cleanly in a phenopacket. Don't replace one with the
    other.
    """

    case_id: str
    phenopacket: Phenopacket
    policy: SourcePolicy


class CohortSource(Protocol):
    """A cohort source yields CaseBundles.

    Implementations live under `agentic_reranking.sources.<source_id>`.
    """

    source_id: str
    version: str
    policy: SourcePolicy

    def iter_case_ids(self) -> Iterable[str]:
        """Yield every case identifier this source can materialize."""
        ...

    def materialize(self, case_id: str) -> CaseBundle:
        """Load one case fully into memory."""
        ...

    # Deferred until a second source arrives — see "Multi-source readiness"
    # in docs/plans/2026-04-29-agentic-reranking-prototype.md:
    #
    # def leak_controls(self, cutoff: date) -> Sequence[LeakControl]: ...
    # def execution_policy(self) -> ExecutionPolicy: ...

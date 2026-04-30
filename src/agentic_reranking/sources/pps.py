"""Phenopacket Store cohort source.

Reads from a Phenopacket Store release directory laid out as::

    <root>/
        phenopackets.csv      # one row per case
        pubdates.json         # PMID -> {earliest_iso, pubdate, epubdate, sortpubdate}
        extracted/<release>/<cohort>/<filename>.json
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from agentic_reranking.sources.base import CaseBundle, Phenopacket, SourcePolicy


@dataclass(frozen=True)
class _PpsRow:
    """One row of phenopackets.csv."""

    file: str  # relative to extracted/<release>/, e.g. "GTF2H4/PMID_40924475_XP140BR.json"
    cohort: str
    disease_id: str
    disease_label: str
    pmids: tuple[str, ...]


def _case_id_from_file(file_rel: str) -> str:
    """Stable, unique case id derived from the cohort/filename pair.

    e.g. ``GTF2H4/PMID_40924475_XP140BR.json`` -> ``GTF2H4__PMID_40924475_XP140BR``.
    Using ``__`` as separator avoids collisions with single underscores in patient ids.
    """
    stem = file_rel[:-5] if file_rel.endswith(".json") else file_rel
    return stem.replace("/", "__")


class PhenopacketStoreSource:
    """`CohortSource` for the Monarch Initiative Phenopacket Store."""

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self.source_id = "pps"
        # The release directory under extracted/ is named for the version
        # without the leading "v" (e.g. "0.1.26"). Detect it once.
        self.version = self.root.name  # e.g. "v0.1.26"
        extracted = self.root / "extracted"
        self._release_dir = next(extracted.iterdir()) if extracted.is_dir() else extracted
        self.policy = SourcePolicy(
            source_id="pps",
            source_version=self.version,
            has_real_vcf=False,
            has_publication_date=True,
            has_phi=False,
            restricted_access=False,
            consent_bounded=False,
            public_literature_source=True,
            extra={"phenopacket_store_release": self.version},
        )
        self._index: dict[str, _PpsRow] = self._build_index()
        self._pubdates: dict[str, dict[str, str]] = self._load_pubdates()

    def _build_index(self) -> dict[str, _PpsRow]:
        idx: dict[str, _PpsRow] = {}
        with (self.root / "phenopackets.csv").open() as fh:
            for row in csv.DictReader(fh):
                pmids = tuple(p for p in row["pmids"].split(";") if p)
                idx[_case_id_from_file(row["file"])] = _PpsRow(
                    file=row["file"],
                    cohort=row["cohort"],
                    disease_id=row["disease_id"],
                    disease_label=row["disease_label"],
                    pmids=pmids,
                )
        return idx

    def _load_pubdates(self) -> dict[str, dict[str, str]]:
        pubdates_path = self.root / "pubdates.json"
        if not pubdates_path.exists():
            return {}
        with pubdates_path.open() as fh:
            return cast(dict[str, dict[str, str]], json.load(fh))

    def iter_case_ids(self) -> Iterable[str]:
        return iter(self._index.keys())

    def materialize(self, case_id: str) -> CaseBundle:
        if case_id not in self._index:
            raise KeyError(f"unknown case_id: {case_id}")
        row = self._index[case_id]
        with (self._release_dir / row.file).open() as fh:
            phenopacket = cast(Phenopacket, json.load(fh))
        return CaseBundle(case_id=case_id, phenopacket=phenopacket, policy=self.policy)

    # PPS-specific helpers used by the manifest builder. They are not on the
    # CohortSource Protocol because not every source has them (GEL/EHR don't
    # have publication dates or PMIDs).

    def get_row(self, case_id: str) -> _PpsRow:
        return self._index[case_id]

    def earliest_pub_date(self, case_id: str) -> str | None:
        """Earliest publication ISO date across all PMIDs for this case, or None."""
        row = self._index[case_id]
        dates = [self._pubdates[p]["earliest_iso"] for p in row.pmids if p in self._pubdates]
        return min(dates) if dates else None

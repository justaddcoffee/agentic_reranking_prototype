# Agentic Reranking Prototype — Implementation Plan

> **For Claude:** Use `${SUPERPOWERS_SKILLS_ROOT}/skills/collaboration/executing-plans/SKILL.md` to implement this plan task-by-task.

**Goal:** Build a pilot benchmark that measures whether an agentic AI (Claude Code on Opus 4.6) improves the differential-diagnosis output of Exomiser on a small post-cutoff cohort of rare-disease cases by leveraging literature search, demographic context, and HPO ontology browsing — under a strict anti-leakage protocol.

**Architecture:** Per case: (1) a `CohortSource` plugin (PPS for this pilot) yields a `CaseBundle` = (GA4GH phenopacket + policy sidecar) for each case in the cohort; (2) spike the causal variant into a NA12878 GIAB exome VCF (only when the source's policy says `has_real_vcf=False`); (3) run `pheval.exomiser` on the resulting case to produce a top-50 candidate gene list; (4) hand the case + Exomiser ranking to a Claude Code agent running in a per-case Docker container (shandy `Dockerfile.agent` / `JobContainerRunner` pattern — non-root user, mounted job dir, normal egress) whose **tool surface** is restricted via `claude-agent-sdk` configuration to two MCP servers (HPO ontology lookup, date-filtered PubMed search) — no Bash, WebFetch, WebSearch, or Write/Edit outside the job dir; (5) parse the agent's re-ranked list; (6) score against held-out ground truth. Five ablation conditions isolate the contributions of literature search and demographic context. Knowledge sources are pinned where possible to a 2025-08-31 reference date: HPO ontology and the Exomiser 2508 data bundle use frozen dated snapshots; PubMed is queried live with a strict date filter (a documented approximation — see Task 11). HPO, demographics, and the cohort's gene–disease links are filtered to ensure post-cutoff novelty. **All raw inputs live under `data/raw/<source>/<version>/`; reference snapshots under `data/snapshots/`; per-experiment outputs under `data/derived/runs/<run_id>/`** so additional sources (GEL, dark EHR, future phenopacket repos) can drop in without disturbing the pilot's outputs.

**Tech Stack:**
- Python 3.12, `uv` for env + deps, `pytest` for tests, `ruff` for lint/format, `mypy` for types
- `pheval` and `pheval.exomiser` (Monarch Initiative) for phenopacket handling and Exomiser orchestration
- `pyphetools` / `phenopacket-tools` for phenopacket I/O
- HPO ontology via `pyobo` or direct OBO Foundry releases; `hpo-toolkit` if needed
- NCBI E-utilities (`Bio.Entrez`) for the date-filtered PubMed wrapper
- MCP server stack: `mcp` Python SDK, MCP servers run as stdio child processes inside the per-case agent container
- Docker for per-case filesystem/process isolation (shandy `Dockerfile.agent` + `JobContainerRunner` pattern); container has normal network egress (so the agent reaches `api.anthropic.com` and lit-search reaches NCBI). Leak prevention is at the agent's tool surface, not the network layer.
- `claude-agent-sdk` invoked from inside the container; model locked to `claude-opus-4-6`, temperature=0; built-in tools (Bash, WebFetch, WebSearch, Task, NotebookEdit) disabled; Write/Edit restricted to the job output directory; only the two MCP servers exposed.
- `numpy`, `pandas`, `scipy.stats` for metrics and bootstrap CIs; `matplotlib` for the headline plot

**Scope notes:**
- This is a pilot. N=22 (potentially fewer after the cohort gate), 5 diseases, 2 dominant. The plan deliberately does not scale beyond this — a future plan can swap in a larger cohort once Phenopacket Store releases later 2026 cases.
- No publication or paper-targeted polish in this plan. Goal is reproducible numbers + a results table + one plot.

---

## Workstreams overview

| # | Workstream | Tasks | Approx effort |
|---|---|---|---|
| 1 | Project scaffolding | T1 | 2h |
| 2 | Cohort curation | T2–T4 | 2 days |
| 3 | Background genome + spiking | T5–T6 | 1 day |
| 4 | Phenotype + demographic perturbation | T7–T8 | 1 day |
| 5 | Cohort gate (vanilla Exomiser must rank truth) | T9 | 1 day |
| 6 | Snapshots + lit-search wrapper | T10–T12 | 2 days |
| 7 | MCP servers (HPO, lit-search) | T13–T14 | 2 days |
| 8 | Per-case container + tool-surface restriction (shandy pattern) | T14 | 1 day |
| 9 | Agent harness | T15–T17 | 3 days |
| 10 | Ablation runner | T18 | 1 day |
| 11 | Scoring + reporting (with stratified bootstrap, harm rate) | T19–T21 | 2 days |
| 12 | End-to-end execution + writeup | T22–T24 | 2 days |
| 13 | Reproducibility audit | T25 | 0.5 day |

**Sequencing constraint:** Workstreams 2→3→4→5 are linear (each consumes the previous). Workstream 6 can run in parallel after 1. Workstream 7 needs 6. Workstream 8 needs 7. Workstream 9 needs 8 and 5. Workstream 10 needs 9. Workstream 11 needs 10. The critical path is roughly 1→2→3→4→5→9→10→11→12 (~3 weeks of sequential work).

---

## Task 1: Project scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `README.md`
- Create: `src/agentic_reranking/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/test_smoke.py`
- Create: `Makefile`

**Step 1: Initialize uv project**

```bash
cd /Users/jtr4v/PythonProject/agentic_reranking_prototype
uv init --python 3.12 --package
```

**Step 2: Add dependencies**

```bash
uv add pheval phenopackets pyobo requests pandas numpy scipy matplotlib biopython tenacity pyyaml
uv add --dev pytest pytest-cov ruff mypy types-requests types-PyYAML
```

**Step 3: Configure ruff and mypy**

In `pyproject.toml`, add:

```toml
[tool.ruff]
line-length = 100
target-version = "py312"

[tool.ruff.lint]
select = ["E", "F", "I", "B", "UP", "N", "RUF"]

[tool.mypy]
python_version = "3.12"
strict = true
ignore_missing_imports = true

[tool.pytest.ini_options]
testpaths = ["tests"]
```

**Step 4: Write smoke test**

```python
# tests/test_smoke.py
def test_smoke():
    import agentic_reranking
    assert agentic_reranking is not None
```

**Step 5: Run tests**

`uv run pytest -q` → expect 1 passed.

**Step 6: Add Makefile shortcuts**

```makefile
.PHONY: lint test fmt typecheck
lint: ; uv run ruff check src tests
fmt:  ; uv run ruff format src tests
typecheck: ; uv run mypy src
test: ; uv run pytest -q
all: fmt lint typecheck test
```

**Step 7: Commit**

```bash
git add pyproject.toml uv.lock README.md src tests Makefile .gitignore
git commit -m "scaffold: uv project, ruff, mypy, smoke test"
```

---

## Task 2: Cohort source interface + PPS adapter + run manifest

**Files:**
- Create: `src/agentic_reranking/sources/__init__.py`
- Create: `src/agentic_reranking/sources/base.py` — `CohortSource` Protocol, `CaseBundle`, `SourcePolicy` dataclass
- Create: `src/agentic_reranking/sources/pps.py` — Phenopacket Store adapter (the only source for this pilot)
- Create: `src/agentic_reranking/cohort/__init__.py`
- Create: `src/agentic_reranking/cohort/manifest.py` — generic run-manifest builder that consumes any `CohortSource`
- Create: `data/raw/pps/v0.1.26/.gitkeep` (the actual artifacts move here from `pps/`)
- Output: `data/derived/runs/<run_id>/manifest.tsv`
- Test: `tests/sources/test_base.py`, `tests/sources/test_pps.py`, `tests/cohort/test_manifest.py`

**Goal:** Establish the `CohortSource` Protocol so future sources (GEL, dark EHR, additional phenopacket repos) can drop in without changes to downstream pipeline code, then implement a single PPS adapter for the pilot. The Phenopacket Store adapter consumes `data/raw/pps/v0.1.26/` and yields `CaseBundle = (phenopacket: Phenopacket, policy: SourcePolicy)` per case. The run-manifest builder applies the cohort filter (post-2025-09-01 publication) and writes the per-run manifest.

**Design contract (load-bearing for multi-source readiness — read carefully):**

```python
# src/agentic_reranking/sources/base.py

@dataclass(frozen=True)
class SourcePolicy:
    """Operational metadata that doesn't fit cleanly in a phenopacket.
    Per-source leak controls and execution policy are derived from this."""
    source_id: str                  # "pps", "gel", "ehr_<site>"
    source_version: str             # "v0.1.26", "100k_release_15", "2026-Q1"
    has_real_vcf: bool              # GEL=True; PPS=False (we spike)
    has_publication_date: bool      # PPS=True; GEL/EHR=False
    has_phi: bool                   # EHR=True; PPS/GEL=False
    restricted_access: bool         # GEL=True; PPS=False
    consent_bounded: bool           # GEL/EHR=True
    public_literature_source: bool  # PPS=True
    extra: dict[str, str] = field(default_factory=dict)  # source-specific tags

@dataclass(frozen=True)
class CaseBundle:
    """Glue type: phenopacket holds clinical content (Exomiser-aligned),
    policy holds operational metadata. Don't replace one with the other."""
    phenopacket: dict               # GA4GH Phenopacket v2 as parsed JSON
    policy: SourcePolicy
    case_id: str

class CohortSource(Protocol):
    source_id: str
    version: str
    policy: SourcePolicy

    def iter_case_ids(self) -> Iterable[str]: ...
    def materialize(self, case_id: str) -> CaseBundle: ...
    # Deferred (not implemented in pilot, but documented in contract):
    #   def leak_controls(self, cutoff: date) -> Sequence[LeakControl]: ...
    #   def execution_policy(self) -> ExecutionPolicy: ...
```

For the pilot, the PPS adapter sets `has_publication_date=True`, `public_literature_source=True`, all other flags `False`. The deferred `leak_controls` / `execution_policy` methods are documented in the Protocol but not yet required — they'll fill in when GEL/EHR sources arrive.

**Goal of the run-manifest builder:** Iterate `CohortSource`, materialize each case, apply the (currently PPS-specific) post-2025-09-01 filter, write `data/derived/runs/<run_id>/manifest.tsv`.

**Step 0: Move existing PPS artifacts into the new layout**

```bash
mkdir -p data/raw/pps/v0.1.26
git mv pps/all_phenopackets.zip pps/extracted pps/phenopackets.csv pps/pubdates.json pps/pmids.txt pps/extract_pmids.py pps/fetch_pubdates.py pps/analyze.py data/raw/pps/v0.1.26/
rmdir pps  # only if empty
# update .gitignore to reference the new paths
```

**Step 1: Failing tests** for the Protocol + the PPS adapter.

```python
# tests/sources/test_base.py
from agentic_reranking.sources.base import CohortSource, CaseBundle, SourcePolicy

def test_source_policy_is_immutable():
    p = SourcePolicy(source_id="pps", source_version="v0.1.26", has_real_vcf=False,
                    has_publication_date=True, has_phi=False, restricted_access=False,
                    consent_bounded=False, public_literature_source=True)
    with pytest.raises(FrozenInstanceError):
        p.has_phi = True

# tests/sources/test_pps.py
from agentic_reranking.sources.pps import PhenopacketStoreSource

def test_pps_yields_case_bundles_with_correct_policy():
    src = PhenopacketStoreSource(root="data/raw/pps/v0.1.26")
    case_ids = list(src.iter_case_ids())
    assert len(case_ids) == 9588  # known total in v0.1.26
    bundle = src.materialize(case_ids[0])
    assert bundle.policy.source_id == "pps"
    assert bundle.policy.has_publication_date is True
    assert bundle.policy.has_real_vcf is False
    assert "id" in bundle.phenopacket  # GA4GH phenopacket loaded as dict

def test_pps_extracts_hpo_terms_from_phenopacket():
    bundle = PhenopacketStoreSource("data/raw/pps/v0.1.26").materialize("...")
    hpo_ids = [pf["type"]["id"] for pf in bundle.phenopacket["phenotypicFeatures"]]
    assert all(h.startswith("HP:") for h in hpo_ids)
```

```python
# tests/cohort/test_manifest.py
def test_manifest_filters_by_publication_date():
    src = PhenopacketStoreSource("data/raw/pps/v0.1.26")
    manifest = build_run_manifest(src, run_id="test_run",
                                   filter_pub_date_from=date(2025, 9, 1))
    assert len(manifest) == 22  # 22 post-2025-09 cases
    assert all(row["source_pubdate"] >= "2025-09-01" for row in manifest)
```

**Step 2: Run tests** → fail.

**Step 3: Implement** `sources/base.py`, `sources/pps.py`, `cohort/manifest.py`. The PPS adapter reads from `data/raw/pps/v0.1.26/phenopackets.csv` + `pubdates.json` + `extracted/`. The manifest builder takes any `CohortSource` and a filter; for PPS it applies `source_pubdate >= 2025-09-01`. Manifest output: `data/derived/runs/<run_id>/manifest.tsv` with columns `case_id, source_id, source_version, source_pmid, source_pubdate, sex, age_years, ancestry, hpo_terms, disease_id, disease_label, truth_gene_symbol, truth_gene_id, truth_hgvs, novelty_class`. (`novelty_class` ∈ {`novel`, `established`} — set `novel` for Valence-Farazi and AR BGC-11, else `established`.)

**Step 4: Run tests** → pass.

**Step 5: Run the driver**

```bash
RUN_ID=$(date -u +%Y%m%dT%H%M%SZ)-pilot
uv run python -m agentic_reranking.cohort.manifest \
  --source pps --source-version v0.1.26 \
  --filter min-pub-date=2025-09-01 \
  --run-id "$RUN_ID"
```

Expected: writes `data/derived/runs/$RUN_ID/manifest.tsv` and `data/derived/runs/$RUN_ID/run.yaml` (provenance: source/version/filters/timestamp/git-sha). Sanity-check: 22 rows + header.

**Step 6: Commit**

```bash
git add src/agentic_reranking/sources src/agentic_reranking/cohort tests/sources tests/cohort
git commit -m "sources: CohortSource Protocol + PPS adapter + run manifest"
```

---

## Task 3: HPO term loader and ontology snapshot

**Files:**
- Create: `src/agentic_reranking/hpo/__init__.py`
- Create: `src/agentic_reranking/hpo/snapshot.py`
- Test: `tests/hpo/test_snapshot.py`

**Goal:** Pin HPO ontology to the latest tagged release ≤ 2025-08-31. Load it once; expose lookup functions used by perturbation (Task 7) and the `hpo_terms` MCP server (Task 13).

**Step 1: Identify the right tag.** HPO releases at https://github.com/obophenotype/human-phenotype-ontology/releases — find the release tagged `v2025-XX-XX` with `XX-XX ≤ 08-31`. As of writing, `v2025-08-13` is the candidate. Pin in code.

**Step 2: Failing test**

```python
# tests/hpo/test_snapshot.py
from agentic_reranking.hpo.snapshot import load_snapshot, get_parents, get_ancestors

def test_loaded_snapshot_returns_known_term():
    onto = load_snapshot()
    assert "HP:0001250" in onto
    assert onto["HP:0001250"].label == "Seizure"

def test_parents_of_seizure_includes_neurological_phenotype():
    onto = load_snapshot()
    parents = get_parents(onto, "HP:0001250")
    assert any(p.startswith("HP:") for p in parents)
```

**Step 3: Run test** → fails with import.

**Step 4: Implement** — download the tag's `hp.obo` once into `data/snapshots/hpo/v2025-08-13/hp.obo`, parse with `pyobo` or `obonet`, expose:

```python
def load_snapshot() -> dict[str, HpoTerm]: ...
def get_parents(onto, term_id: str) -> list[str]: ...
def get_ancestors(onto, term_id: str) -> list[str]: ...
def get_children(onto, term_id: str) -> list[str]: ...
def get_label(onto, term_id: str) -> str: ...
```

`HpoTerm` dataclass: `id, label, definition, parents, children`. **Do not** load or expose disease–phenotype annotations (HPOA) — that's the answer key.

**Step 5: Run test** → expect pass.

**Step 6: Commit**

```bash
git add src/agentic_reranking/hpo tests/hpo data/snapshots/hpo/.gitkeep
git commit -m "hpo: pinned 2025-08 ontology snapshot loader"
```

---

## Task 4: Background genome (NA12878 GIAB exome VCF)

**Files:**
- Create: `data/snapshots/giab/README.md`
- Create: `scripts/fetch_giab.sh`
- No test (data download — verified by checksum)

**Goal:** Obtain the GIAB NA12878 (HG001) hg38 whole-exome VCF used by PhEval for variant spiking. We use the same source PhEval uses for reproducibility.

**Step 1: Fetch script**

```bash
#!/usr/bin/env bash
# scripts/fetch_giab.sh
set -euo pipefail
DEST=data/snapshots/giab
mkdir -p "$DEST"
URL="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TBI="$URL.tbi"
curl -L -o "$DEST/HG001.vcf.gz" "$URL"
curl -L -o "$DEST/HG001.vcf.gz.tbi" "$TBI"
echo "Downloaded $(stat -f%z "$DEST/HG001.vcf.gz") bytes"
```

**Step 2: Run** — `bash scripts/fetch_giab.sh`. Expect ~hundreds of MB.

**Step 3: Restrict to exome regions** — using a hg38 exome BED (e.g., refseq_exons.hg38.bed, included via dependency or downloaded), do `bcftools view -R exons.bed -O z -o data/snapshots/giab/HG001.exome.vcf.gz HG001.vcf.gz; tabix -p vcf HG001.exome.vcf.gz`. (If `bcftools` not installed: `brew install bcftools`.)

**Step 4: Document checksums** in `data/snapshots/giab/README.md` so reproducibility is auditable.

**Step 5: Commit script + README only** (gitignore handles the actual VCFs).

```bash
git add scripts/fetch_giab.sh data/snapshots/giab/README.md
git commit -m "data: GIAB NA12878 exome VCF fetch script + README"
```

---

## Task 5: Variant spiking utility

**Files:**
- Create: `src/agentic_reranking/spiking/__init__.py`
- Create: `src/agentic_reranking/spiking/spike.py`
- Test: `tests/spiking/test_spike.py`

**Goal:** Given a background VCF and an HGVS-c notation (from the manifest), produce a per-case VCF with the variant injected. Reuse PhEval's HGVS→VCF normalization where available; otherwise wrap `hgvs` library.

**Step 1: Failing test** — using a minimal fixture VCF (~50 lines) and a known variant.

```python
# tests/spiking/test_spike.py
from agentic_reranking.spiking.spike import spike_variant

def test_spike_adds_variant(tmp_path, fixture_bg_vcf):
    out = tmp_path / "spiked.vcf.gz"
    spike_variant(
        background_vcf=fixture_bg_vcf,
        hgvs_c="NM_000546.5:c.215C>G",
        gene_symbol="TP53",
        output_vcf=out,
    )
    # spiked variant present
    import gzip
    lines = gzip.open(out, "rt").read().splitlines()
    body = [l for l in lines if not l.startswith("#")]
    assert any("17\t" in l and "C\tG" in l for l in body)
    # original variants preserved
    assert len(body) > len([l for l in gzip.open(fixture_bg_vcf, "rt") if not l.startswith("#")])
```

**Step 2: Run test** → fail.

**Step 3: Implement** — uses `pyhgvs` or `hgvs` to resolve HGVS-c → genomic coords on hg38, then appends a VCF record (genotype `0/1` for autosomal dominant or `1/1` for AR — read inheritance from the manifest's `disease_id`, default `0/1`). Sort & re-index with `bcftools sort` + `tabix`.

**Step 4: Run test** → pass.

**Step 5: Edge cases** — add tests for: (a) variant on chrX where sex matters, (b) HGVS-c with version mismatch, (c) variant already present in background (should overwrite or warn).

**Step 6: Commit.**

---

## Task 6: Per-case spiked VCF generator

**Files:**
- Create: `src/agentic_reranking/spiking/build_cohort.py`
- No new test (driver script; correctness covered by Task 5)

**Goal:** Iterate the run manifest, call `spike_variant` for each row whose source policy has `has_real_vcf=False`, write to `data/derived/runs/<run_id>/cases/<source_id>/<case_id>/variants/spiked.vcf.gz`. Sources with `has_real_vcf=True` (future GEL) skip spiking and pass through normalization/QC instead.

**Step 1: Implement driver, dry-run flag, parallelism (4 cores) via `concurrent.futures`.** Driver takes `--run-id` and resolves all input/output paths from it.

**Step 2: Run** — `uv run python -m agentic_reranking.spiking.build_cohort --run-id "$RUN_ID" --background data/snapshots/giab/HG001.exome.vcf.gz`.

**Step 3: Sanity check** — every case has a non-empty spiked VCF.

**Step 4: Commit driver.**

---

## Task 7: Phenotype perturbation

**Files:**
- Create: `src/agentic_reranking/perturbation/__init__.py`
- Create: `src/agentic_reranking/perturbation/phenotype.py`
- Test: `tests/perturbation/test_phenotype.py`

**Goal:** Implement HPO-term perturbation modeled on PhEval's noise-injection patterns. Three operations, each with a configurable rate:
- `drop_rate` (default 0.20): randomly remove terms
- `generalize_rate` (default 0.30): replace term with a parent (or grandparent if no parent)
- `noise_count` (default 1–2): add this many random HPO terms drawn from a clinically-relevant pool (configurable; default = all phenotypic-abnormality descendants)

Seed all random ops from `case_id` so perturbation is deterministic per case.

**Step 1: Failing tests** — one per operation + one for determinism.

```python
def test_drop_reduces_terms_deterministically(): ...
def test_generalize_replaces_with_parent(): ...
def test_noise_adds_terms(): ...
def test_same_case_id_produces_same_perturbation(): ...
```

**Step 2: Run** → fails.

**Step 3: Implement**, using the HPO snapshot from Task 3.

**Step 4: Run** → pass.

**Step 5: Run on real cohort** — generate `data/derived/runs/<run_id>/perturbation/phenotypes.tsv`. Sanity-check: a few cases by eye.

**Step 6: Commit.**

---

## Task 8: Demographic perturbation

**Files:**
- Modify: `src/agentic_reranking/perturbation/__init__.py`
- Create: `src/agentic_reranking/perturbation/demographics.py`
- Test: `tests/perturbation/test_demographics.py`

**Goal:** Shift age by ±1–2 years (clamped ≥0), preserve sex if disease is sex-linked (lookup table for our 5 diseases), randomize ancestry within plausible options. Seeded by case_id.

**Step 1–6:** Same TDD cycle as Task 7. Output `data/derived/runs/<run_id>/perturbation/demographics.tsv`.

---

## Task 9: Cohort gate (vanilla Exomiser must find truth gene)

**Files:**
- Create: `src/agentic_reranking/exomiser/__init__.py`
- Create: `src/agentic_reranking/exomiser/runner.py`
- Create: `src/agentic_reranking/cohort/gate.py`
- Test: `tests/exomiser/test_runner.py` (integration test, marked `@pytest.mark.slow`)

**Goal:** Run `pheval.exomiser` on each perturbed case (perturbed phenotypes + the case's variants — spiked VCF for sources with `has_real_vcf=False`, pass-through real VCF otherwise). Drop any case where the truth gene is outside the top-100 — those cases are unsolvable from the perturbed input and would create noise. Output `data/derived/runs/<run_id>/gated_manifest.tsv` (subset of the run manifest).

**Step 1:** Per the `pheval.exomiser` README, set up its config (Exomiser jar, phenotype.zip 2402, application.properties). **Pin everything explicitly in `config/exomiser.yaml`:**
- Exomiser version (tag)
- Exomiser data bundle: **2508** (released 2025-09-26). Confirmed via [discussion #611](https://github.com/exomiser/Exomiser/discussions/611) that all data sources are pre-2025-08-31: ClinVar 2025-08-12, OMIM 2025-08-12, Orphanet 2025-06-24, HPO 2025-05-06, HPO annotations 2025-05-06. (Exomiser uses YYMM convention: `2508` = August 2025.) Compatible with Exomiser binary 14.0.0, 14.1.0, 15.0.0. Direct download: `2508_hg38.zip`, `2508_phenotype.zip` from the URLs listed in discussion #611.
- `pheval.exomiser` plugin version (commit hash or tag — the plugin is under active development; unpinned versions can change scoring behavior)
- Random seed (Exomiser's older versions use randomized graph traversals in HiPHIVE; set `--seed 42` if the runner accepts it)
- Java version + JVM flags for deterministic ordering (`-Djava.util.Arrays.useLegacyMergeSort=true`, `-XX:+UseSerialGC`)

**Step 2:** Wrap a single-case run as `run_exomiser(case_id, hpo_terms, vcf_path) -> ExomiserResult`. Result: ranked list of `(rank, gene_symbol, score, variant_coords)` tuples.

**Step 3:** Integration test (slow): one case end-to-end produces a non-empty ranking.

**Step 4:** Implement the gate: iterate cohort, call runner, keep cases whose truth gene is in top-100. Log dropped cases.

**Step 5:** Run gate. Expected: ~22 cases survive (some may drop). Commit `gated_manifest.tsv` and the dropped-cases log.

**Step 6:** Commit.

---

## Task 10: ClinVar / OMIM / HPOA decision

ClinVar and OMIM are NOT exposed as agent tools. HPOA is NOT exposed (it's the answer key). No work to do here — this task is just a `docs/decisions/0001-no-clinvar-omim.md` ADR that captures the reasoning. Skipping the TDD cycle.

**Files:**
- Create: `docs/decisions/0001-no-clinvar-omim-hpoa-tools.md`

Content: 1 page summarizing why these tools were dropped (leak vectors, redundancy with Exomiser, lit-search covers the same information space organically).

Commit.

---

## Task 11: PubMed lit-search wrapper

**Files:**
- Create: `src/agentic_reranking/litsearch/__init__.py`
- Create: `src/agentic_reranking/litsearch/wrapper.py`
- Create: `src/agentic_reranking/litsearch/pmid_blocklist.txt` (the 22 source PMIDs)
- Test: `tests/litsearch/test_wrapper.py`

> **⚠ Approximation, not a true snapshot.** This wrapper queries **live** NCBI E-utilities with a `maxdate` filter — it is **not** a frozen 2025-08-31 PubMed dump. NCBI updates retroactively (MeSH terms, corrected publication dates, late-indexing of issued papers, etc.), so the same query may return slightly different results on different dates. This is an explicit, documented compromise per the project owner's preference (a true PubMed snapshot would require a static dump of ~36M records and is out of scope for the pilot). Residual risks: (a) post-2025-08 MeSH re-indexing may surface or hide pre-cutoff papers in unexpected ways; (b) a paper with a journal-issue date after 2025-08-31 but an earlier `epubdate` may slip through depending on which date the wrapper filters on; (c) NCBI database errors or downtime may yield non-deterministic empty results. Mitigations: the source-PMID blocklist (Step 4 below) and the build-time validation step (Step 5) provide defense-in-depth. The PMID blocklist is **not** redundant given live PubMed — it is the load-bearing safeguard against (a) and (b).

**Goal:** A thin wrapper over NCBI E-utilities (`esearch` + `esummary`) with these constraints baked in:
- All queries get `mindate=0001/01/01&maxdate=2025/08/31` appended
- Filter on `epubdate` if present, else `pubdate` — whichever is earliest, must be ≤ 2025-08-31
- Only allow query patterns: `<gene_symbol>[Title/Abstract] AND year`, `<HGVS_string>`, `<gene_symbol> AND <gene_symbol>` (gene–gene). Reject anything else.
- Returned fields: `pmid, title, abstract, authors, pubdate, mesh_terms`. **Drop**: `cited_by`, `linked_records`, `pmcid`, `related_articles`.
- Hard PMID blocklist applied post-fetch: any PMID in `pmid_blocklist.txt` is removed from results before return.

**Step 1: Failing tests**

```python
def test_query_appends_date_filter(): ...
def test_free_text_phenotype_query_rejected(): ...
def test_returned_fields_are_stable_set_only(): ...
def test_post_2025_08_paper_not_returned(): ...  # uses real PubMed; mark @pytest.mark.network
def test_blocklisted_pmids_filtered(): ...
def test_filter_uses_epubdate_when_earlier_than_pubdate(): ...
```

**Step 2:** Implement using `Bio.Entrez`. Validate query with a regex allowlist; reject otherwise with `QueryNotAllowed`. Apply blocklist after fetch.

**Step 3:** Build the blocklist — `scripts/build_pmid_blocklist.py` reads `data/derived/runs/<run_id>/manifest.tsv` and appends any new source PMIDs to `data/snapshots/pubmed/pmid-blocklists/<source_id>.txt` (per-source blocklists; the lit-search wrapper unions all source-specific blocklists at runtime). Re-run whenever the cohort changes. Note: only public-literature sources (where `policy.public_literature_source=True`) contribute PMIDs; private cohorts (GEL/EHR) typically don't have source PMIDs to blocklist.

**Step 4:** Add a build-time validation script `scripts/validate_no_source_paper_findable.py` that, for each cohort case, attempts queries on the truth gene + reasonable variants and confirms the source PMID does not appear (both with and without the blocklist applied — the "without" case tells us whether the date filter alone is doing its job; the "with" case is the production-path check).

**Step 5:** Schedule the validation script to re-run **before every full pipeline execution** (added to `scripts/run_all.sh` as a precondition gate). If a previously-blocked PMID slips through filtering at any later run, fail loudly. This is the load-bearing reproducibility guardrail given the live-PubMed compromise.

**Step 6:** Commit.

---

## Task 12: HPO term info MCP server

**Files:**
- Create: `tools/hpo_terms/server.py`
- Create: `tools/hpo_terms/pyproject.toml`
- Test: `tools/hpo_terms/tests/test_server.py`

**Goal:** Stdio MCP server exposing four tools:
- `lookup_term(hp_id) -> {label, definition, parents[], children[]}`
- `find_term_by_label(label) -> {hp_id, label}` (exact match only)
- `get_ancestors(hp_id) -> [hp_id]`
- `get_descendants(hp_id, max_depth=3) -> [hp_id]`

No disease cross-references in any response.

**Step 1:** Failing tests using the `mcp` Python client.

**Step 2:** Implement using the snapshot loader from Task 3.

**Step 3:** Sanity-check by manually invoking via `mcp inspect`.

**Step 4:** Commit.

---

## Task 13: Lit-search MCP server

**Files:**
- Create: `tools/lit_search/server.py`
- Test: `tools/lit_search/tests/test_server.py`

**Goal:** Stdio MCP server exposing one tool: `search(query, max_results=20) -> [{pmid, title, abstract, authors, pubdate, mesh_terms}]`. Wraps the Task 11 wrapper.

**Step 1–4:** Same pattern.

---

## Task 14: Per-case agent container (shandy pattern)

**Files:**
- Create: `Dockerfile.agent`
- Create: `Dockerfile.base`
- Create: `docker/agent-entrypoint.py`
- Create: `src/agentic_reranking/job_container/__init__.py`
- Create: `src/agentic_reranking/job_container/runner.py`
- Test: `tests/job_container/test_runner.py`

**Reference implementation:** `~/PythonProject/shandy/Dockerfile.agent`, `~/PythonProject/shandy/docker/agent-entrypoint.py`, and `~/PythonProject/shandy/src/openscientist/job_container/runner.py`. The shandy pattern is the model — adapt, don't reinvent.

**Design philosophy (changed from initial plan):** the Docker container provides **filesystem and process isolation, not network isolation.** The agent container has normal egress so it can reach `api.anthropic.com` (for Claude itself) and the lit-search wrapper can reach `eutils.ncbi.nlm.nih.gov`. Leak prevention happens at the **agent's tool surface** — restricted via `claude-agent-sdk` configuration — not at the network layer.

**Container constraints:**
- Non-root user (`agent`, uid 1001)
- Mounts: read-only case input bundle at `/agent/job/input/`, writable output at `/agent/job/output/`, the two MCP server packages at `/agent/tools/` (read-only)
- No host-system mounts beyond the job dir
- Container is short-lived: one run per case per condition, then destroyed
- Resource limits: 4 GB memory, 1 CPU, 10-minute wall-clock cap

**Tool-surface restriction (load-bearing):** `claude-agent-sdk` configured with:
- **Disabled tools:** `Bash`, `WebFetch`, `WebSearch`, `Task`, `NotebookEdit`. `Write` and `Edit` permitted only inside `/agent/job/output/`.
- **Enabled MCP servers:** `hpo_terms` (Task 12), `lit_search` (Task 13). Both stdio-served as child processes inside the container.
- Model: `claude-opus-4-6`. Temperature: 0. Max tokens: 8192. Tool-call budget: 30.

**Step 1: Dockerfile.base** — `python:3.12-slim` + uv + git. Mirrors shandy's `Dockerfile.base` structure.

**Step 2: Dockerfile.agent** — extends base, installs project (`uv pip install --system -e .`), creates `agent` user, copies `agent-entrypoint.py`, sets `WORKDIR /agent/job`. Mirrors shandy's `Dockerfile.agent`.

**Step 3: `docker/agent-entrypoint.py`** — reads `JOB_ID`, `CASE_ID`, `CONDITION` from env, loads input bundle, invokes `claude-agent-sdk` with the configured tool surface, writes output YAML. Mirrors shandy's `agent-entrypoint.py`.

**Step 4: `JobContainerRunner`** — Python class wrapping the Docker SDK to launch one container per (case, condition) pair, mount the right paths, capture exit code, return parsed output. Mirrors shandy's `runner.py` (drop the postgres / database-status pieces — we don't need a job manager).

**Step 5: Failing tests**

```python
def test_runner_launches_container_with_correct_mounts(): ...
def test_runner_passes_case_and_condition_via_env(): ...
def test_runner_captures_output_yaml(): ...
def test_runner_enforces_resource_limits(): ...
def test_runner_kills_container_on_timeout(): ...
def test_disabled_tools_actually_disabled(): ...   # smoke test inside container
```

**Step 6:** Implement runner. Build images once, reuse across cases.

**Step 7: Smoke test** — launch one container with a trivial input, confirm:
- Agent connects to `api.anthropic.com` ✅
- Agent cannot run `Bash` (disabled-tool error) ✅
- Agent CAN call `mcp__hpo_terms__lookup_term` ✅
- Output YAML is written to `/agent/job/output/result.yaml` ✅

**Step 8: Commit.**

**Note on supply-chain risk:** because the container has internet egress, a compromised image could exfiltrate data. We accept this risk — same as shandy, same as any cloud-hosted agent. Mitigated by: pinning base image digest, building only from this repo's Dockerfiles, and not running this on production data.

---

## Task 15: Agent harness — prompt template + invocation

**Files:**
- Create: `src/agentic_reranking/agent/__init__.py`
- Create: `src/agentic_reranking/agent/harness.py`
- Create: `src/agentic_reranking/agent/prompt_template.md`
- Test: `tests/agent/test_harness.py`

**Goal:** Build the per-case agent invocation:
1. Load case from `gated_manifest.tsv` + perturbed phenotypes/demographics + Exomiser top-50
2. Render prompt template
3. Invoke Claude Code (Opus 4.6) via `claude-agent-sdk` inside the per-case container with the two MCP servers attached
4. Capture: tool-call log + final YAML output

**Pinned settings (load-bearing for reproducibility):**
- Model: `claude-opus-4-6`
- Temperature: **0** (deterministic sampling)
- Max tokens: 8192
- Tool-call budget: 30
- System prompt + user prompt seeded from versioned templates; any template change forces a rerun

If we later want to characterize sampling variance, run a single condition × case at temperature=0.7 with N=5 replicates as a one-off study. Default for the headline pipeline is temperature=0.

**Step 1: Failing test** — uses a stubbed Anthropic client returning a known YAML.

**Step 2:** Implement; render prompt with `string.Template` or jinja2.

**Step 3:** Add `--mock` flag for local testing without hitting the API.

**Step 4:** Commit.

---

## Task 16: Output parser

**Files:**
- Create: `src/agentic_reranking/agent/parser.py`
- Test: `tests/agent/test_parser.py`

**Goal:** Parse agent output into `RerankedList = list[RerankedGene]` where `RerankedGene = (rank: int, gene_symbol: str, source: Literal["exomiser_top_50", "lit_search_added"], confidence: Literal["high","medium","low"])`. Validate: ranks are 1..N contiguous, gene symbols look like HGNC symbols (alphanumeric + dashes/underscores), no duplicates.

**Step 1: Failing tests** for: valid input, missing fields, malformed YAML, duplicate gene symbols, non-contiguous ranks.

**Step 2:** Implement with explicit validation; on malformed input, raise `AgentOutputError` carrying the raw output for retry logic.

**Step 3:** Commit.

---

## Task 17: Error handling + retry policy

**Files:**
- Modify: `src/agentic_reranking/agent/harness.py`
- Test: `tests/agent/test_harness_errors.py`

**Goal:** Implement the error-handling policy: malformed YAML → 1 retry with format reminder → fall back to baseline + flag `agent_failed`; tool-call budget exceeded → score with what we have; lit-search transient → exponential backoff (3 attempts).

**Step 1–4:** Add tests for each failure mode using mocked clients.

---

## Task 18: Ablation conditions runner

**Files:**
- Create: `src/agentic_reranking/agent/ablations.py`
- Test: `tests/agent/test_ablations.py`

**Goal:** Five conditions (per Section 4):
1. `exomiser_baseline` — no agent run; emits Exomiser top-50 directly as the "ranking"
2. `agent_no_tools` — agent invoked, but no MCP servers attached. Memorization probe.
3. `agent_hpo_only` — only `hpo_terms` MCP attached; no lit-search.
4. `agent_no_demographics` — both tools, but `demographics` field blanked in input.
5. `agent_full` — both tools, demographics included.

Each condition produces `data/derived/runs/<run_id>/cases/<source_id>/<case_id>/agent/<condition>/result.yaml`.

**Step 1:** Define `Condition` dataclass: tool list, demographics-blanked bool, condition name.

**Step 2:** Tests verify the right tools are wired for each condition.

**Step 3:** Implement runner that iterates conditions × cases. Idempotent (skips already-completed runs).

**Step 4:** Commit.

---

## Task 19: Per-case scorer

**Files:**
- Create: `src/agentic_reranking/scoring/__init__.py`
- Create: `src/agentic_reranking/scoring/per_case.py`
- Test: `tests/scoring/test_per_case.py`

**Goal:** For one case in one condition, compute: `top_1, top_3, top_5, top_10, reciprocal_rank, recall_expanded, harm` and the booleans `agent_failed`.

**Metric definitions:**
- `recall_expanded` (binary): truth gene was outside Exomiser top-50 AND inside agent top-10. Denominator for the reported rate: number of cases where truth gene was outside Exomiser top-50 (i.e., conditional on the recall-expansion opportunity actually existing). Report both numerator and denominator separately.
- `harm` (binary): truth gene's rank moved *down* (numerically larger) from Exomiser baseline to agent output. Captures rank regression — a clinically important failure mode that the prior plan didn't measure.
- Per-condition harm rate = fraction of cases with `harm == True`. Reported alongside top-N improvements; an agent that helps half the cases and hurts the other half is not the same as one that helps half and leaves the rest unchanged.

**Step 1: Failing tests** with hand-crafted ranking lists, including specific tests for: truth gene out of top-50 (recall_expansion eligible), truth gene moved down (harm), truth gene moved up (improvement), no change.

**Step 2:** Implement.

**Step 3:** Commit.

---

## Task 20: Aggregator with stratified bootstrap CIs

**Files:**
- Create: `src/agentic_reranking/scoring/aggregate.py`
- Test: `tests/scoring/test_aggregate.py`

**Goal:** Across cases, compute means + 95% bootstrap CIs (1000 resamples, fixed seed `42`).

**Stratification:** Use a **stratified bootstrap that resamples within disease clusters**, not a naive bootstrap. With 9 cases of disease A and 6 of disease B, a naive bootstrap can produce CIs narrower than the true uncertainty because it doesn't preserve the cluster imbalance. The stratified version resamples each disease's cases independently, then aggregates.

**Reporting:**
- Always pair percentages with raw counts: "Top-1: 5/22 (22.7%)" not "Top-1: 22.7%"
- Stratify by **disease cluster** in addition to `novelty_class` — emit a per-disease breakdown table where each row is one of the 5 diseases with N=9, 6, 5, 1, 1
- Demote MRR from headline metrics; include in supplementary table only

**Outputs:**
- `results/summary.tsv` — one row per condition × overall
- `results/per_disease.tsv` — one row per condition × disease cluster
- `results/per_novelty.tsv` — novel vs established split (kept for completeness)

**Step 1:** Failing tests (synthetic input, known expected means + CI ranges; verify stratified bootstrap produces wider CIs than naive on imbalanced input).

**Step 2:** Implement.

**Step 3:** Commit.

---

## Task 21: Headline plot

**Files:**
- Create: `src/agentic_reranking/scoring/plot.py`
- Output: `results/topn_curves.png`

**Goal:** One plot, top-N curves overlaid for all 5 conditions. Stratified-overlay variant: novel vs established panels side-by-side.

**Step 1:** Implement (matplotlib only, no seaborn).

**Step 2:** Commit.

---

## Task 22: End-to-end execution

**Files:**
- Create: `scripts/run_all.sh`

**Goal:** A single command that runs the entire pipeline for a given source + run-id.

```bash
#!/usr/bin/env bash
set -euo pipefail

SOURCE="${1:-pps}"
SOURCE_VERSION="${2:-v0.1.26}"
RUN_ID="${3:-$(date -u +%Y%m%dT%H%M%SZ)-$SOURCE}"

make all                                            # lint + typecheck + tests
uv run python -m agentic_reranking.cohort.reproducibility_audit --strict   # gate (Task 25)
bash scripts/fetch_giab.sh                          # idempotent
uv run python -m agentic_reranking.cohort.manifest --source "$SOURCE" --source-version "$SOURCE_VERSION" --run-id "$RUN_ID"
uv run python -m agentic_reranking.spiking.build_cohort --run-id "$RUN_ID"
uv run python -m agentic_reranking.perturbation.run_all --run-id "$RUN_ID"
uv run python -m agentic_reranking.cohort.gate --run-id "$RUN_ID"
uv run python -m agentic_reranking.exomiser.run_baseline --run-id "$RUN_ID"
uv run python -m scripts.validate_no_source_paper_findable --run-id "$RUN_ID"   # litsearch leak gate (Task 11 step 5)
uv run python -m agentic_reranking.agent.ablations --run-id "$RUN_ID" --conditions all
uv run python -m agentic_reranking.scoring.aggregate --run-id "$RUN_ID"
uv run python -m agentic_reranking.scoring.plot --run-id "$RUN_ID"
echo "Run complete: data/derived/runs/$RUN_ID/"
```

**Step 1:** Implement.

**Step 2:** Run end-to-end with `--smoke` flag (1 case, 1 condition) to validate plumbing before the full run.

**Step 3:** Full run. Commit results files.

---

## Task 23: Memorization probe analysis

**Files:**
- Create: `docs/results/memorization-probe.md`

**Goal:** Write a short analysis of the `agent_no_tools` condition: per-case top-1 hit rate, how it compares to vanilla Exomiser baseline, what it implies about contamination floor for our 22 cases. Include a per-disease breakdown.

If memorization is high (>30% top-1), explicitly recommend abandoning the 22-case cohort and switching to hand-curated 2026 cases.

---

## Task 24: Results writeup

**Files:**
- Create: `docs/results/2026-04-29-pilot-results.md`

**Goal:** A short results memo: cohort description, ablation results (table + plot), interpretation, threats to validity, what to do next. Not paper-quality; pilot-quality.

**Required reporting elements:**
- All percentages paired with raw counts (e.g., "Top-1: 5/22 (22.7%)")
- Per-disease-cluster breakdown table alongside the aggregate
- Harm-rate alongside improvement metrics
- Memorization-probe result (Task 23) reported as a contamination floor; main metrics framed as "agent value above probe floor," not just "agent value above Exomiser"
- Explicit pilot framing: "N=22 is insufficient for statistical significance on individual metrics; this study aims at signal detection and pipeline validation, not confirmatory inference"
- Stochasticity note: agent runs at temperature=0, but Exomiser scoring is deterministic only if all the Task 9 pins hold; document the test that confirms this

---

## Task 25: Reproducibility audit

**Files:**
- Create: `scripts/reproducibility_audit.py`
- Create: `docs/reproducibility.md`

**Goal:** A single script that audits all the determinism / pinning knobs the pilot depends on, run as the final pre-execution check. If any knob is unpinned or drifted, fail loudly. Output a manifest of versions + hashes that gets stamped into every results file.

**Audited items:**
- Exomiser version, data bundle version, plugin commit hash (Task 9)
- HPO ontology snapshot tag + file SHA-256 (Task 3)
- ClinVar snapshot date + file SHA-256
- `pheval.exomiser` commit hash
- `claude-agent-sdk` version
- Model ID (`claude-opus-4-6`) and temperature (0)
- Bootstrap seed (42)
- Random seeds for perturbation (Tasks 7, 8)
- PMID blocklist size + SHA-256
- The 22 source PMIDs themselves (sanity check that the cohort hasn't drifted)
- Twin-run determinism test: re-run one case in one condition twice and verify identical Exomiser output and identical agent output

**Step 1: Implement the audit script.**

**Step 2:** Run twice — once now (expect: lots of "not yet pinned" warnings as the project scaffolds), once before any results are reported (expect: all green). The audit's "all green" run is a precondition for `scripts/run_all.sh`.

**Step 3:** Commit.

---

## Open questions / known risks

1. **Cohort gate dropout rate.** If perturbation is too aggressive, the gate may drop most of the cohort. Calibrate perturbation strength on a 3-case dry run before applying to all 22.
2. **OMIM/HPOA leak via HPO ontology.** HPO term *labels* contain disease-flavored phrasing (e.g., "Netherton-like ichthyosis"). Audit the perturbed input for any term whose label uniquely identifies a disease; consider redacting term labels and showing only IDs to the agent. Decide before running ablations.
3. **Variant location vs. NA12878 ancestry.** NA12878 is European; for cases with ancestry-specific genetic backgrounds, spiking may produce slightly miscalibrated allele-frequency context. Note in writeup; not blocking.
4. **Live PubMed reproducibility.** Per Task 11, the lit-search wrapper queries live NCBI with a date filter, not a frozen snapshot. NCBI's retroactive MeSH re-indexing and pubdate corrections mean results can drift between runs. Mitigations: PMID blocklist + the validation check at every full-pipeline run. Residual risk: if a source PMID slips past both, we may not detect it. Mitigation 2: re-run the validation check **before every results-publishing run** and store its output alongside results.
5. **Container supply-chain risk.** The agent container has internet egress (so the agent reaches `api.anthropic.com`). A compromised Docker base image could exfiltrate data. Mitigations: pin base image digest, build only from this repo's Dockerfiles, do not run on real patient data. Standard cloud-agent threat model; we accept this risk for a pilot on public phenopackets.
6. **HPO term-label redaction decision.** Decide whether to expose HPO term labels to the agent or only term IDs. Labels are useful (the agent reasons better with them), but for nearly-pathognomonic phenotypes the label itself is a leak vector. Tentative: keep labels for the pilot, audit for pathognomonic terms manually, redact case-by-case if needed. Revisit after the memorization probe.
7. **Naming-probe threshold calibration.** The cohort gate currently fails if >50% of cases are named verbatim by the no-tools probe. The 50% threshold is arbitrary without a null distribution. Plan: run the no-tools probe on a held-out set of post-cutoff diseases unrelated to the cohort, compute the empirical false-positive naming rate, and recalibrate the threshold from that. Capture as a sub-step under Task 23.
8. **Memorization can leak through phenotype association without naming the disease.** The probe currently checks if the agent names the disease. An agent that has memorized HPO→gene associations could re-rank correctly without ever naming the disease. Add a secondary probe: compare top-1 hit rate of the no-tools agent vs. vanilla Exomiser; if the no-tools agent significantly outperforms Exomiser, that itself is evidence of memorization regardless of whether disease names appeared. Add to Task 23.
9. **Inter-rater reliability on perturbation choices.** The HPO perturbation strategy has experimenter degrees of freedom (which terms to drop, generalize, add). For the pilot, document each per-case perturbation deterministically (seeded by `case_id` per Task 7) and include the perturbation log as supplementary material. A real follow-up should add a second annotator.
10. **Qualitative reasoning-trace review.** With N=22, quantitative metrics alone cannot distinguish "agent reasoning well" from "agent got lucky on the dominant disease." Add to Task 24: blinded qualitative review of agent reasoning traces on a random 5 of 22 cases.

---

## Multi-source readiness

The pilot has only one cohort source (Phenopacket Store v0.1.26), but future cohorts will likely include GEL (100K Genomes), dark EHR exports, and possibly other phenopacket repositories. This plan adopts the *minimal* layer needed to make those additions a plugin rather than a refactor, and explicitly defers the deeper machinery until a second source actually arrives.

### What this plan adopts now

- **Directory layout**: `data/raw/<source>/<version>/` for raw inputs, `data/snapshots/` for pinned reference data (HPO ontology, Exomiser bundle, GIAB, PubMed query policy + PMID blocklists), `data/derived/runs/<run_id>/` for per-experiment outputs. Run-scoped derived artifacts let multiple experiments coexist.
- **`CohortSource` Protocol** (Task 2): minimal contract — `iter_case_ids()` and `materialize(case_id) -> CaseBundle`. PPS is the only implementation today.
- **`CaseBundle`** = `(phenopacket: dict, policy: SourcePolicy, case_id: str)`. The phenopacket is canonical for clinical content (Exomiser-aligned); the policy sidecar carries operational metadata (capability flags, consent status, source version) that doesn't fit cleanly in phenopackets.
- **`SourcePolicy`** capability flags: `has_real_vcf`, `has_publication_date`, `has_phi`, `restricted_access`, `consent_bounded`, `public_literature_source`. Pipeline branches on these (e.g., spiking only runs when `has_real_vcf=False`) so adding GEL flips a flag rather than rewriting code.
- **Per-source PMID blocklist files** (`data/snapshots/pubmed/pmid-blocklists/<source_id>.txt`). PPS contributes; GEL/EHR typically don't.
- **Run manifest as the single source of truth per experiment** (`data/derived/runs/<run_id>/manifest.tsv` + `run.yaml` provenance file). Pipeline tools take `--run-id` and resolve all paths from it.

### What this plan deliberately defers (until a second source arrives)

- **`leak_controls(cutoff)` and `execution_policy()` methods on `CohortSource`** — documented in the Protocol comment but not implemented. The PPS pipeline uses today's hardcoded controls (PMID blocklist, post-2025-09-01 filter, etc.). When GEL or EHR is added, these methods get implemented and the shared controls split out from source-specific ones.
- **Generalized `CaseFilter` enum** (`SOURCE_PUBLICATION_DATE` / `ENCOUNTER_DATE` / `VARIANT_INTERPRETATION_DATE` / `INGESTED_AT`). Today's pilot uses a hardcoded `min-pub-date` flag because it only needs that one filter type.
- **Execution profiles for restricted/federated sources** — the current Task 14 container has public egress (fine for PPS, *not* fine for GEL/EHR). When a private cohort lands, add a second container profile with stricter egress and run-on-host-only mode.
- **PHI scanning pre-prompt and small-cell suppression in outputs** — needed for EHR, irrelevant for PPS.
- **Continuous-source versioning** for sources whose state changes between snapshots (EHR exports). Today's PPS adapter pins to v0.1.26 and is done.
- **Genome-build dispatch and VCF normalization policies** — GEL-specific. PPS doesn't ship VCFs.

### Stop-gap: what the next person picking up GEL/EHR should expect

When a second source arrives, the work is:
1. Implement `<NewSource>Source(CohortSource)` in `src/agentic_reranking/sources/<source>.py`.
2. Decide its `SourcePolicy` capability flags.
3. Implement `leak_controls(cutoff)` and `execution_policy()` (these now become required).
4. Generalize `CaseFilter` to multi-field if needed.
5. Add a new container profile under `Dockerfile.agent.<source>` if egress requirements differ.
6. Update `scripts/run_all.sh` to dispatch on source — already takes `--source` as arg 1.

Tasks T3, T5, T7, T8, T9, T15-T20 should remain source-agnostic and not need changes.

---

## What's deliberately NOT in this plan

- Real-time clinical deployment (a "realistic clinical context" framing where the agent has the same internet access as a clinician would) — different cohort, different infrastructure
- Comparison against other prioritization tools (LIRICAL, Phen2Gene) — out of scope for pilot
- Larger cohort, recurring runs, or paper polish
- HPO term-label redaction (deferred to "open questions" — decide based on audit)
- Full citation graph search, related-article expansion, or any tool whose return shape evolves over time
- The deferred multi-source machinery listed above — wait for the second source before building it

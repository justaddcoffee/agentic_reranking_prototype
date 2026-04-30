# agentic-reranking

Pilot benchmark measuring whether an agentic AI (Claude Code on Opus 4.6) improves the differential-diagnosis output of [Exomiser](https://exomiser.readthedocs.io/) on a small post-cutoff cohort of rare-disease cases by leveraging literature search, demographic context, and HPO ontology browsing, with some efforts to address possible leakage.

See [`docs/plans/2026-04-29-agentic-reranking-prototype.md`](docs/plans/2026-04-29-agentic-reranking-prototype.md) for the full implementation plan.

## Development

```bash
uv sync           # install deps + dev deps
make all          # fmt + lint + typecheck + test
```

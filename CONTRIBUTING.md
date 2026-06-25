# Contributing to Locus

Thanks for your interest! Locus is a local-first genome explorer; contributions are welcome.

## Ground rules

- **Never commit personal data.** Everything under `data/` (genomes, reference/annotation
  databases, the DuckDB store, reports, `manifest.json`) is `.gitignore`d and must stay that way.
  Before pushing, check `git status` shows no `*.vcf.gz` / `*.duckdb` / `data/` files.
- **Keep it local.** Don't add code that uploads genetic data anywhere; the MCP server and API
  bind to localhost only.

## Dev setup

```bash
make tools        # bcftools, samtools, htslib, openjdk, uv (Homebrew)
make setup        # uv sync
uv run pytest -q  # full pipeline test on a synthetic genome (no real data needed)
```

Before opening a PR:

```bash
uv run ruff check src tests scripts   # lint (line length 120)
uv run pytest -q                      # tests
cd web && npm run build               # type-check + build the SPA (if you touched web/)
```

The synthetic fixture (`scripts/make_fixture.py`) lets you exercise the whole pipeline without
any real genome. See [CLAUDE.md-style architecture notes in `docs/integration-notes.md`] for the
domain gotchas (gVCF conventions, contig naming, GRCh38).

## Scope

Locus targets Apple Silicon macOS today (native tooling, no Docker). Cross-platform support is
welcome but unverified — call out platform assumptions in your PR.

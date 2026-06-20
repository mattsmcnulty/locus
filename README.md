# Locus

**Explore your genome locally with Claude.**

Locus takes the whole-genome sequencing files from [sequencing.com](https://sequencing.com)
(Illumina DRAGEN output), annotates them with open-source clinical, population, and
pharmacogenomic databases, loads everything into a fast local [DuckDB](https://duckdb.org)
store, and exposes it two ways:

1. **An MCP server** so you can ask **Claude** questions about your genome in plain English.
2. **A local debug SPA** (React) for browsing and visualizing the same data.

> ⚠️ **This is sensitive personal genetic data.** Everything under `data/` is `.gitignore`d
> and never committed. The MCP server and the SPA bind to **localhost only**.

---

## What you feed it

The sequencing.com / DRAGEN bundle:

| File | What it is |
| --- | --- |
| `*.snp-indel.genome.vcf.gz` | Small variants (SNVs + indels) — a genome gVCF |
| `*.cnv.vcf.gz` | Copy-number variants |
| `*.sv.vcf.gz` | Structural variants |
| `*.1.fq.gz` / `*.2.fq.gz` | Raw paired-end reads (optional; not needed to query) |

The variants are already **called**, so the VCFs are the queryable layer — Locus does not
re-process the raw reads.

## Architecture

```
sequencing.com VCFs
      │  ingest    (index, QC, normalize — bcftools)
      ▼
  normalized VCFs
      │  annotate  (ClinVar · dbSNP · gnomAD · VEP/SnpEff · PharmCAT)
      ▼
 annotated variants
      │  load      (cyvcf2 → Arrow → DuckDB)
      ▼
  locus.duckdb  ──►  MCP server   ──►  Claude
              └──►  FastAPI       ──►  React SPA
```

## Quickstart

```bash
# 1. System tools (bcftools, samtools, htslib, java) + uv
make tools

# 2. Python env + deps
make setup

# 3. Put your sequencing.com VCFs in data/genome/ then check the environment
make doctor

# 4. Download & prepare the open-source databases (reference, ClinVar, SnpEff, PharmCAT)
uv run locus download all     # see `locus download --help` for individual targets

# 5. Build the local genome database (ingest → annotate → load)
make pipeline

# 6a. Query with Claude — register the MCP server, then just ask Claude
claude mcp add --scope project --transport stdio locus -- uv run locus-mcp

# 6b. Or browse in the debug SPA
make serve-api      # backend at http://127.0.0.1:8787
make web-dev        # frontend dev server at http://localhost:5173
```

Before you have the real genome, you can exercise the whole pipeline on a tiny
**synthetic** one (no downloads, no real data):

```bash
uv run python scripts/make_fixture.py /tmp/locus-fixture
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus ingest && \
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus load
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus serve api   # then open the SPA
```

### Annotation defaults (and why)

- **Consequences:** default to **SnpEff** (pure-Java, runs natively on Apple Silicon). Ensembl
  **VEP** is supported as the richer option but needs Docker on arm64 — see `docs/integration-notes.md`.
- **gnomAD** allele frequencies are **streamed per-region** from the AWS mirror (no multi-hundred-GB
  download). This runs during `locus annotate --steps gnomad` and can be slow over the network.
- **dbSNP** is skipped by default — DRAGEN already fills rsIDs.
- **PharmCAT** runs via the official `pgkb/pharmcat` Docker image (bundles the required VCF preprocessor).

## Querying with Claude (MCP)

The MCP server exposes read-only tools over your genome:

- `lookup_variant(gene | rsid | region)`
- `clinical_findings(gene?, significance?)` — ClinVar pathogenic / likely-pathogenic
- `pharmacogenomics(drug?, gene?)` — PharmCAT star alleles + CPIC guidance
- `allele_frequency(variant)` — gnomAD population frequency
- `gene_summary(gene)` — all variants in a gene with consequences
- `cnv_sv_overlap(gene | region)` — structural / copy-number hits
- `run_sql(query)` — guarded read-only SQL

Then ask things like *"Do I carry any pathogenic ClinVar variants?"* or
*"What's my CYP2C19 metabolizer status and which drugs does it affect?"*

### Registering the MCP server

**Claude Code** (project-scoped, via the committed `.mcp.json`):

```bash
claude mcp add --scope project --transport stdio locus -- uv run --directory $(pwd) locus-mcp
```

**Claude Desktop (macOS)** — add to `~/Library/Application Support/Claude/claude_desktop_config.json`
under `mcpServers` (use **absolute paths**; Claude Desktop does not inherit your shell `PATH`), then
fully quit and reopen Claude Desktop:

```json
{
  "mcpServers": {
    "locus": {
      "command": "/opt/homebrew/bin/uv",
      "args": ["run", "--directory", "/Users/you/locus", "locus-mcp"]
    }
  }
}
```

## Configuration

Paths are env-overridable (see `.env.example`). To keep the genome and the large
reference/annotation databases on an external drive:

```bash
LOCUS_DATA_DIR=/Volumes/genome/locus-data
```

## Layout

```
src/locus/
  config.py        # env-driven paths (keeps data out of git)
  cli.py           # the `locus` CLI
  ingest.py        # index, QC, normalize
  annotate.py      # ClinVar / dbSNP / gnomAD / VEP / PharmCAT
  load.py          # build the DuckDB store
  db.py            # DuckDB access (read-only for queries)
  mcp_server.py    # the Claude interface
  api.py           # FastAPI backend for the SPA
web/               # Vite + React debug SPA
data/              # (gitignored) genome, reference, annotation DBs, locus.duckdb
```

## Disclaimer

Locus is for personal exploration and education. It is **not** a medical device and its
output is **not** medical advice. Discuss any health-relevant finding with a qualified
clinician or genetic counselor.

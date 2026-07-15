# Locus

**Explore your genome locally with Claude.**

Locus turns the whole-genome sequencing files from [sequencing.com](https://sequencing.com)
(or any Illumina/DRAGEN-style VCFs) into a fast local [DuckDB](https://duckdb.org) store,
annotates them with open-source clinical, population, and pharmacogenomic databases, layers
on deeper interpretation (ancestry, polygenic risk, traits, GWAS), and keeps re-interpreting
your genome as new studies are published — watching PubMed and the GWAS Catalog for papers about
the variants *you* carry — all on your machine. It's exposed two ways:

1. **An MCP server** so you can ask **Claude** questions about your genome in plain English.
2. **A local debug SPA** (React) for browsing and visualizing the same data.

> ⚠️ **This is sensitive personal genetic data.** Everything under `data/` is `.gitignore`d
> and never committed. The MCP server and the SPA bind to **localhost only**, and the refresh
> engine and literature lookups send only *generic* queries outward (release dates, public score
> IDs, rsID lists, and gene symbols) — your genome never leaves the machine.

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
re-process the raw reads. (Read-level analyses like repeat expansions or CNV-aware CYP2D6 would
need the BAM/CRAM and are not yet implemented.)

## What it can tell you

| Area | Capability |
| --- | --- |
| **Clinical** | ClinVar pathogenic / likely-pathogenic variants you carry; **ACMG SF v3.2** secondary findings in the ~81 medically-actionable genes |
| **Predicted impact** | **AlphaMissense** scores for missense variants ClinVar has never classified (the `clnsig = null ≠ benign` gap) |
| **Pharmacogenomics** | **PharmCAT** star-allele diplotypes, metabolizer phenotypes, and CPIC/DPWG drug guidance |
| **Ancestry** | Biogeographic ancestry via PCA + k-NN over 1000 Genomes + HGDP — continental and sub-continental |
| **Polygenic risk** | PGS Catalog scores reported as an **ancestry-matched percentile** (CAD, LDL, T2D, Lp(a), …) |
| **Traits & wellness** | Single-SNP traits (lactose, caffeine, alcohol flush, earwax, eye color, muscle type), the **HLA-B\*57:01** abacavir-hypersensitivity proxy, and your **mtDNA maternal haplogroup** |
| **GWAS breadth** | Which genome-wide-significant (p<5e-8) risk alleles you carry, across the whole GWAS Catalog, queryable by trait |
| **On-demand** | Paste a new paper's rsIDs and get your genotypes live (`ask_about`) |
| **Population frequency** | **gnomAD** allele frequencies (via Ensembl) on the variants where rarity matters — those ClinVar has classified or AlphaMissense calls pathogenic |
| **Literature** | New papers about the **specific variants you carry** (NCBI **LitVar2**, keyed on your clinically-notable rsIDs) and new **PubMed** papers on your notable genes + new **GWAS** associations at your variants, all surfaced into the changelog; ask for the latest research (`literature_for`) or "which variants did this study find that I have?" (`variants_in_study`) |
| **Living updates** | `locus refresh` re-interprets your genome as databases move — **ClinVar reanalysis** (variants you carry newly classified pathogenic), plus new GWAS associations and PubMed papers on your genes |

## Architecture

```
sequencing.com VCFs
      │  ingest    (index, QC, normalize, chr-canonicalize — bcftools)
      ▼
  normalized VCFs
      │  annotate  (ClinVar · SnpEff · AlphaMissense · PharmCAT native; gnomAD AF via Ensembl)
      ▼
 annotated variants
      │  load      (cyvcf2 → Arrow → DuckDB)
      ▼
  locus.duckdb  ◄── ancestry/pgs · traits · gwas   (deeper interpretation; preserved across reloads)
      │       ◄── refresh                          (re-interpret as ClinVar/GWAS/LitVar/PubMed/PGS/CPIC publish)
      ├──►  MCP server   ──►  Claude
      └──►  FastAPI      ──►  React SPA
```

Everything runs **natively on Apple Silicon — no Docker.** PLINK2 (arm64), bcftools/samtools,
and the JVM tools (SnpEff, PharmCAT, Haplogrep2) all run directly.

## Quickstart

> **Just want it working (or setting it up for a non-developer)?** On an Apple Silicon Mac,
> double-click **`setup.command`** (or run `./setup.sh`) — it installs everything, downloads the
> databases, builds your genome store, and registers the MCP server with Claude. See
> [docs/SETUP.md](docs/SETUP.md). The manual steps below are the same thing, broken out.

```bash
# 1. System tools (bcftools, samtools, htslib, java) + uv
make tools

# 2. Python env + deps
make setup

# 3. Put your sequencing.com VCFs in data/genome/, then check the environment
make doctor

# 4. Download & prepare the open-source databases (the table below; `all` fetches everything, ~9 GB)
uv run locus download all

# 5. Build the local genome database (ingest → annotate → load)
make pipeline                # ~5M variants; a few minutes

# 6. Deeper interpretation (each is independent and preserved across variant reloads)
uv run locus ancestry        # ancestry + ancestry-calibrated polygenic risk
uv run locus traits          # single-SNP traits + HLA-B*57:01 proxy + mtDNA haplogroup
uv run locus gwas            # GWAS Catalog risk alleles you carry

# 7a. Query with Claude — register the MCP server (below), then just ask
# 7b. Or browse the SPA
uv run locus serve api       # http://127.0.0.1:8787  (serves the built SPA + the API)

# 8. Keep it current — re-interpret as new studies publish
uv run locus refresh --dry-run     # preview; then `locus refresh`, or schedule it weekly
```

Before you have the real genome, exercise the whole pipeline on a tiny **synthetic** one
(no downloads, no real data):

```bash
uv run python scripts/make_fixture.py /tmp/locus-fixture
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus ingest
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus load
LOCUS_DATA_DIR=/tmp/locus-fixture uv run locus serve api
```

## CLI reference

Run any command with `uv run locus <command>` (or activate the venv / add an alias — see below).
`--help` works on every command and subcommand.

| Command | What it does |
| --- | --- |
| `locus doctor` | Check the toolchain, databases, and input data are in place |
| `locus download <target>` | Download & prepare a database (or `all`) — see targets below |
| `locus ingest` | Index, QC, normalize the VCFs (gVCF → sites; canonicalize contigs to `chr*`) |
| `locus annotate [--steps ...]` | Annotate variants. Steps: `clinvar,gnomad,snpeff,alphamissense,pharmcat` or `all` |
| `locus load` | Build/refresh the DuckDB store (rebuilds variant tables; **preserves** ancestry/PGS/traits/GWAS/watch tables) |
| `locus pipeline` | `ingest → annotate → load`, end to end |
| `locus ancestry` | Biogeographic ancestry (PCA + k-NN over 1000G + HGDP) + ancestry-calibrated polygenic risk |
| `locus traits` | Genotype single-SNP traits/wellness + HLA-B\*57:01 proxy + mtDNA haplogroup |
| `locus gwas` | Genotype GWAS Catalog risk alleles (p<5e-8) and store the ones you carry |
| `locus literature <gene\|rsID\|PMID>` | Recent PubMed papers on a gene/rsID, or (for a PubMed ID) which variants that study reported that you carry |
| `locus refresh [--dry-run] [--force] [--sources clinvar,pgs,cpic,gwas,pubmed,litvar]` | Check tracked sources for new releases and re-interpret what changed (ClinVar reanalysis, new PGS, CPIC updates, new GWAS associations at your variants, new LitVar/PubMed papers on your variants & genes) |
| `locus schedule install [--weekday N] [--hour H]` | Install a weekly macOS launchd job that runs `locus refresh` |
| `locus schedule status` / `locus schedule uninstall` | Show / remove the scheduled job |
| `locus serve mcp` | Start the MCP server (stdio) — what Claude connects to |
| `locus serve api` | Start the FastAPI backend and serve the built SPA at `http://127.0.0.1:8787` |

**Typical first run:** `download all` → `pipeline` → `ancestry` → `traits` → `gwas` → register MCP / `serve api` → `refresh` (and optionally `schedule install`).

**Run `locus` from anywhere** (instead of `uv run locus …`): add a shell alias —

```bash
echo "alias locus='uv run --directory $(pwd) locus'" >> ~/.zshrc && source ~/.zshrc
```

### Download targets

| Target | Contents | Approx size |
| --- | --- | --- |
| `reference` | GRCh38 no-alt analysis-set FASTA (bgzip + faidx) | ~900 MB |
| `clinvar` | ClinVar GRCh38 VCF, chr-renamed | ~180 MB |
| `snpeff` | SnpEff jar + GRCh38 database (consequences) | ~0.5 GB |
| `pharmcat` | PharmCAT pipeline — native jar + Python preprocessor | ~30 MB |
| `alphamissense` | AlphaMissense hg38 (slim, tabix-indexed) | ~640 MB |
| `ancestry` | PLINK2 (arm64) + 1000 Genomes + HGDP reference panels | ~6 GB |
| `haplogrep` | Haplogrep2 jar (mtDNA haplogroups) | ~7 MB |
| `gwas` | NHGRI-EBI GWAS Catalog associations TSV | ~65 MB |

`locus download all` fetches **everything** in the table (~9 GB, idempotent and resumable); you can
also fetch any single target by name (e.g. `locus download ancestry`) if you only want some features.

## Querying with Claude (MCP)

The MCP server exposes **read-only, typed** tools over your genome (object-typed output so strict
clients dispatch them reliably):

**Lookups** — `genome_overview`, `lookup_variant_by_rsid`, `lookup_variants_in_gene`,
`lookup_variants_in_region`, `allele_frequency`, `run_sql` (guarded read-only SQL)
**Clinical** — `clinical_findings` (ClinVar P/LP), `predicted_damaging` (rare AlphaMissense-pathogenic),
`secondary_findings` (ACMG SF), `structural_variants` (CNV/SV)
**Pharmacogenomics** — `pharmacogenomics` (PharmCAT diplotypes + CPIC/DPWG guidance)
**Ancestry & risk** — `ancestry`, `polygenic_risk`
**Traits & breadth** — `traits`, `gwas_associations`, `ask_about` (paste rsIDs or a trait)
**Literature** — `literature_for` (recent PubMed on a gene/rsID/topic), `variants_in_study` (which of a paper's variants you carry)
**Living updates** — `whats_new` (the ranked changelog from `locus refresh`)

Then ask things like:
- *"Do I carry any pathogenic ClinVar variants? Any ACMG secondary findings?"*
- *"What's my CYP2C19 metabolizer status and which drugs does it affect?"*
- *"What's my CAD polygenic risk percentile, and where do I fall in ancestry?"*
- *"What's the latest research on my BRCA2, and does this new paper (PMID …) report variants I carry?"*
- *"What's my mtDNA haplogroup? Am I lactose-persistent? Can I take abacavir?"*
- *"What GWAS associations do I carry for type 2 diabetes?"*
- *"Check rs429358 and rs7412"* (APOE) — or paste any rsIDs from a new paper.
- *"What's new in my genome this month?"*

### Registering the MCP server

**Claude Code** (project-scoped, via the committed `.mcp.json`):

```bash
claude mcp add --scope project --transport stdio locus -- uv run --directory $(pwd) locus-mcp
```

**Claude Desktop (macOS)** — add to `~/Library/Application Support/Claude/claude_desktop_config.json`
under `mcpServers` (use **absolute paths**; Claude Desktop does not inherit your shell `PATH`), then
**fully quit and reopen** Claude Desktop:

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

## The living genome (refresh & scheduling)

Your genome is sequenced once and static; the *interpretation databases* move every month.
`locus refresh` polls each source, detects new releases against `data/manifest.json`, re-interprets
only what changed, and writes a ranked "what's new since last run" changelog (the `whats_new` MCP
tool, the SPA **Changelog** tab, and `data/reports/whats_new.md`).

- **ClinVar reanalysis** — snapshots your current classifications, fetches the new ClinVar,
  re-annotates + reloads, and surfaces variants you carry that became (or stopped being) pathogenic,
  tiered by ClinVar review status.
- **LitVar2 (variant-level literature)** — the highest-signal source: new papers about the *specific
  variants you carry* (every clinically-notable rsID — pathogenic, VUS/conflicting, risk/drug-response),
  via NCBI's variant→literature index. *"New paper about rs…, which you carry"* beats a gene-level hit.
- **PubMed (gene-level literature)** — recent papers on your notable genes (ones you carry pathogenic,
  pharmacogenomic, or AlphaMissense-predicted-damaging variants in). Both literature watchers send only
  rsIDs / gene symbols, and are deduped so a paper never surfaces twice.
- **GWAS Catalog** — when a new catalog release lands, re-scans the variants you carry and flags
  newly-published genome-wide-significant associations (weak single hits, not a calibrated score).
  Fully local — nothing variant-specific leaves the machine.
- **PGS Catalog** — reports newly published scores (suggest-only; you add relevant IDs to track).
- **CPIC** — flags pharmacogenomic guideline updates touching genes you carry.

Schedule it weekly (native macOS launchd):

```bash
uv run locus refresh --dry-run     # preview what would change
uv run locus refresh               # run it
uv run locus schedule install      # weekly autopilot (Sun 03:00 by default)
```

## The debug SPA

`locus serve api` serves both the API and the built React app at `http://127.0.0.1:8787`. Tabs:
**Search** (rsID / gene / region), **Clinical**, **Pharmacogenomics**, **Ancestry** (continental +
sub-continental bars + a population PCA scatter), **Risk** (polygenic percentiles), **Traits**,
**GWAS**, **Changelog**, and **SQL**.

To work on the UI, run the Vite dev server for hot-reload (`make web-dev` on :5173, already
CORS-allowed to call the API); otherwise `cd web && npm run build` and refresh.

### Launch as a Mac app

Prefer a Dock icon over a command? Build a native `Locus.app`:

```bash
scripts/build_macos_app.sh        # installs /Applications/Locus.app (DNA icon)
```

It's a tiny native Cocoa app (Swift) that starts `locus serve api`, opens the browser, shows in the
Dock, and **stops the server when you Quit** (Cmd-Q / Dock → Quit). Double-click it, or drag it to the
Dock. (`scripts/make_app_icon.py` renders the icon; `scripts/locus_app.swift` is the source; the `.app`
itself is machine-local, not committed.) Prefer a terminal window with live logs instead? Use
`scripts/locus-serve.command`.

## Configuration

Paths are env-overridable (see `.env.example`, prefix `LOCUS_`). To keep the genome and the large
reference/annotation databases on an external drive:

```bash
LOCUS_DATA_DIR=/Volumes/genome/locus-data
```

## Project layout

```
src/locus/
  config.py        # env-driven paths (keeps data out of git)
  cli.py           # the `locus` CLI
  ingest.py        # index, QC, normalize, chr-canonicalize
  annotate.py      # ClinVar / gnomAD / SnpEff / AlphaMissense / PharmCAT
  load.py          # build the DuckDB store (+ writers for ancestry/traits/gwas/watch)
  db.py            # DuckDB access (read-only for queries)
  queries.py       # shared, typed query layer (used by BOTH the MCP server and the API)
  mcp_server.py    # the Claude interface (typed MCP tools)
  api.py           # FastAPI backend for the SPA
  ancestry.py      # PCA + k-NN ancestry; markers_genotypes() (hom-ref-aware genotyping primitive)
  pgs.py           # PGS Catalog scoring + ancestry-matched calibration
  panels.py        # ACMG SF genes + curated trait/wellness tag SNPs
  traits.py        # single-SNP traits + mtDNA haplogroup
  gwas.py          # GWAS Catalog "associations you carry" + on-demand rsID genotyping
  refresh.py       # the living-genome refresh engine (ClinVar reanalysis, PGS/CPIC watchers)
  manifest.py      # data/manifest.json — source version tracking
  schedule.py      # macOS launchd scheduling for `locus refresh`
  vcfutils.py      # contig/gVCF helpers; shell.py / artifacts.py — shell-outs & canonical paths
web/               # Vite + React debug SPA
scripts/           # make_fixture.py — synthetic genome for tests
data/              # (gitignored) genome, reference, annotation DBs, locus.duckdb, manifest.json
```

## Development

```bash
uv run pytest -q                       # full pipeline test on the synthetic fixture (no real data)
uv run ruff check src tests scripts    # lint (line length 120)
cd web && npm run build                # type-check + build the SPA
```

## Privacy & data

Your genetic data **stays on your Mac**. Everything under `data/` is `.gitignore`d and never
committed; the MCP server and the web app bind to **localhost only**; nothing is uploaded. The
`locus refresh` updater and the literature/`ask_about` lookups send only *generic* queries to
public databases — release dates, public database/score IDs, **rsID lists, and gene symbols**
(to NCBI PubMed/LitVar2, the GWAS Catalog, and Ensembl). Your **genotypes never leave the machine**;
all matching against your genome is done locally.

## Acknowledgments

Locus stands on excellent open tools, databases, and public APIs. It **never redistributes them** —
please respect each source's license:

- **Tools:** [bcftools/samtools/htslib](https://www.htslib.org), [PLINK2](https://www.cog-genomics.org/plink/2.0/)
  (GPLv3), [SnpEff](https://pcingola.github.io/SnpEff/), [PharmCAT](https://pharmcat.org),
  [Haplogrep](https://haplogrep.i-med.ac.at), [DuckDB](https://duckdb.org).
- **Databases** (downloaded to your machine, queried locally):
  [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) (NCBI, public domain),
  the [1000 Genomes Project](https://www.internationalgenome.org)
  & [HGDP](https://www.internationalgenome.org/data-portal/data-collection/hgdp) reference panels,
  the [PGS Catalog](https://www.pgscatalog.org), the [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/),
  and [Phylotree](https://www.phylotree.org).
- **Public APIs** (queried on demand; only rsIDs, gene symbols, and release dates are ever sent):
  [Ensembl REST](https://rest.ensembl.org) — variant coordinates and
  [gnomAD](https://gnomad.broadinstitute.org) allele frequencies;
  [NCBI PubMed + LitVar2](https://www.ncbi.nlm.nih.gov/research/litvar2/) — literature;
  [CPIC](https://cpicpgx.org) — pharmacogenomic guidelines.
- **AlphaMissense** (Google DeepMind): the predictions are licensed **CC-BY-NC 4.0 (non-commercial)**.
  Locus uses them for personal/research interpretation only and does not redistribute them.

Locus is an independent project — **not affiliated with or endorsed by** sequencing.com, Illumina,
Anthropic, or any data provider.

## License

Locus's source code is released under the **MIT License** — see [LICENSE](LICENSE). The third-party
tools and databases it downloads at runtime carry their own licenses (see Acknowledgments above).

## Disclaimer

Locus is for personal exploration and education. It is **not** a medical device and its output is
**not** medical advice. Polygenic percentiles are research-grade and only valid within an
ancestry-matched reference; GWAS single-hit associations are weak and must not be summed; the
HLA-B\*57:01 result is a screening proxy. Discuss any health-relevant finding with a qualified
clinician or genetic counselor.

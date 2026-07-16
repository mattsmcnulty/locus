"""``locus`` command-line interface.

    locus setup       # one-command guided full install (for friends/family)
    locus doctor      # check toolchain + data are in place
    locus mcp install # register the MCP server with Claude (Desktop + Code)
    locus ingest      # index, QC, normalize the sequencing.com VCFs
    locus annotate    # ClinVar / dbSNP / gnomAD / VEP / PharmCAT
    locus load        # build the DuckDB store
    locus pipeline    # ingest -> annotate -> load, end to end
    locus serve mcp   # start the MCP server (query with Claude)
    locus serve api   # start the FastAPI backend (debug SPA)
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from .config import settings

app = typer.Typer(
    name="locus",
    help="Explore your genome locally with Claude.",
    no_args_is_help=True,
    add_completion=False,
)
serve_app = typer.Typer(help="Run the query interfaces.", no_args_is_help=True)
app.add_typer(serve_app, name="serve")
schedule_app = typer.Typer(help="Schedule periodic `locus refresh` (macOS launchd).", no_args_is_help=True)
app.add_typer(schedule_app, name="schedule")
mcp_app = typer.Typer(help="Register the MCP server with Claude.", no_args_is_help=True)
app.add_typer(mcp_app, name="mcp")

console = Console()

# External tools the pipeline shells out to.
REQUIRED_TOOLS = ["bcftools", "samtools", "tabix", "bgzip"]
OPTIONAL_TOOLS = ["java"]

# On macOS, `which java` resolves to a stub that errors unless a JDK is installed.
# The Homebrew openjdk is keg-only, so check its real location too.
_JAVA_CANDIDATES = [
    "/opt/homebrew/opt/openjdk/bin/java",
    "/usr/local/opt/openjdk/bin/java",
]


def _resolve_java() -> str | None:
    for cand in [shutil.which("java"), *_JAVA_CANDIDATES]:
        if not cand or not Path(cand).exists():
            continue
        try:
            out = subprocess.run([cand, "-version"], capture_output=True, text=True, timeout=10)
            text = out.stderr or out.stdout
            if "version" in text.lower() and "unable to locate" not in text.lower():
                return text.strip().splitlines()[0]
        except (subprocess.SubprocessError, OSError):
            continue
    return None


def _tool_version(tool: str) -> str | None:
    if tool == "java":
        return _resolve_java()
    path = shutil.which(tool)
    if not path:
        return None
    for flag in ("--version", "-version", "version"):
        try:
            out = subprocess.run(
                [tool, flag], capture_output=True, text=True, timeout=10
            )
            line = (out.stdout or out.stderr).strip().splitlines()
            if line:
                return line[0]
        except (subprocess.SubprocessError, OSError):
            continue
    return path


@app.command()
def doctor() -> None:
    """Check that the toolchain and input data are ready."""
    table = Table(title="Locus environment", show_lines=False)
    table.add_column("Check")
    table.add_column("Status")
    table.add_column("Detail")

    for tool in REQUIRED_TOOLS:
        v = _tool_version(tool)
        table.add_row(tool, "[green]ok[/]" if v else "[red]MISSING[/]", v or "brew install htslib bcftools samtools")
    for tool in OPTIONAL_TOOLS:
        v = _tool_version(tool)
        table.add_row(f"{tool} (optional)", "[green]ok[/]" if v else "[yellow]absent[/]", v or "—")

    from . import artifacts, download, mcp_install, vcfutils

    # Genome inputs + GRCh38 build sanity (the #1 silent failure is a non-GRCh38 file).
    inputs = artifacts.classify_inputs(settings.genome_dir) if settings.genome_dir.exists() else None
    if inputs and inputs.small_variants:
        try:
            build = vcfutils.detect_build(vcfutils.read_info(inputs.small_variants))
        except Exception:  # noqa: BLE001 - best-effort
            build = "unknown"
        grch37 = build.startswith("GRCh37")
        table.add_row(
            "genome VCFs",
            "[red]GRCh37?[/]" if grch37 else "[green]ok[/]",
            f"{inputs.small_variants.name} (build {build})"
            + ("  ← Locus needs GRCh38; re-export from sequencing.com" if grch37 else ""),
        )
    else:
        table.add_row("genome VCFs", "[yellow]none[/]",
                      f"drop your sequencing.com files into {settings.genome_dir}")

    # Java-backed steps self-skip silently without a JDK — surface that instead of a false "ok".
    java_ok = _resolve_java() is not None

    def _db_row(label: str, present: bool, hint: str, needs_java: bool = False) -> None:
        if present and needs_java and not java_ok:
            table.add_row(label, "[yellow]needs java[/]",
                          "installed, but no Java runtime → this step SKIPS. Fix: brew install openjdk")
            return
        table.add_row(label, "[green]ok[/]" if present else "[yellow]absent[/]", hint)

    ann = settings.annotations_dir
    _db_row("reference FASTA", artifacts.find_reference() is not None, "locus download reference")
    _db_row("ClinVar", (ann / download.CLINVAR_CHR_VCF).exists(), "locus download clinvar")
    _db_row("SnpEff", (ann / "snpEff" / "snpEff.jar").exists(), "locus download snpeff", needs_java=True)
    _db_row("AlphaMissense", (ann / "alphamissense" / "AlphaMissense_hg38.slim.tsv.bgz").exists(),
            "locus download alphamissense")
    _db_row("PharmCAT (native)", (artifacts.pharmcat_install_dir() / "pharmcat_pipeline").exists(),
            "locus download pharmcat", needs_java=True)
    _db_row("ancestry panel",
            (ann / "ancestry" / "all_hg38.pgen").exists() and (settings.data_dir / "tools" / "plink2").exists(),
            "locus download ancestry")
    _db_row("GWAS Catalog", (ann / "gwas" / "gwas-catalog-associations.tsv").exists(), "locus download gwas")
    _db_row("Haplogrep (mtDNA)", (ann / "haplogrep" / "haplogrep.jar").exists(), "locus download haplogrep",
            needs_java=True)

    # MCP registration with Claude.
    reg = mcp_install.is_registered()
    both = reg["desktop"] and reg["desktop_matches"] and reg["code"] and reg["code_matches"]
    table.add_row("MCP (Claude)", "[green]ok[/]" if both else "[yellow]not registered[/]",
                  "Desktop + Code" if both else "run `locus mcp install`")

    table.add_row(
        "DuckDB store",
        "[green]ok[/]" if settings.db_path.exists() else "[yellow]not built[/]",
        str(settings.db_path) if settings.db_path.exists() else "run `locus setup` (or `locus pipeline`)",
    )

    # Coverage, not just file existence. Every annotation step self-skips when its inputs are
    # missing, so the store can be fully built and still have a column that is 100% NULL — a
    # feature that silently answers nothing. Checking files alone reports that as green.
    if settings.db_path.exists():
        for (check, status), detail in _coverage_rows():
            table.add_row(check, status, detail)
    console.print(table)


# What each variant column powers, and which step fills it — so a gap names its own fix.
# `scoped` marks columns only ever filled for a deliberate subset, where a low % is correct
# and only zero is a bug (ClinVar/AlphaMissense/gnomAD annotate a slice of the genome by design).
_COVERAGE_CHECKS = (
    ("consequence", "SnpEff", "consequences + gene names", "locus annotate --steps snpeff", False),
    ("gene", "SnpEff", "gene lookups, literature watch", "locus annotate --steps snpeff", False),
    ("clnsig", "ClinVar", "clinical + secondary findings", "locus annotate --steps clinvar", True),
    ("am_class", "AlphaMissense", "predicted_damaging", "locus annotate --steps alphamissense", True),
    ("gnomad_af", "gnomAD/Ensembl", "allele_frequency, rarity filter", "locus annotate --steps gnomad", True),
)


def _coverage_rows() -> list[tuple[tuple[str, str], str]]:
    """Report how much of the store each annotation actually populated."""
    from .db import connect

    rows: list[tuple[tuple[str, str], str]] = []
    try:
        with connect(read_only=True) as con:
            total = con.execute("SELECT count(*) FROM variants").fetchone()[0]
            if not total:
                return [(("variant coverage", "[yellow]empty[/]"), "no variants loaded — run `locus load`")]
            for col, step, powers, fix, scoped in _COVERAGE_CHECKS:
                try:
                    n = con.execute(f"SELECT count({col}) FROM variants").fetchone()[0]
                except Exception:  # noqa: BLE001 - column absent on an older store
                    n = 0
                if n == 0:
                    status, detail = "[red]NOT APPLIED[/]", f"{step} produced nothing → {powers} dead. Fix: {fix}"
                elif scoped:
                    status, detail = "[green]ok[/]", f"{n:,} variants — {step} (scoped by design)"
                else:
                    status, detail = "[green]ok[/]", f"{n:,}/{total:,} ({100 * n / total:.1f}%) — {step}"
                rows.append(((f"  ↳ {col}", status), detail))
            # Deep-interpretation tables are written by their own steps and preserved across reloads.
            for tbl, fix in (("pgs_scores", "locus ancestry"), ("traits", "locus traits"),
                             ("associations", "locus gwas")):
                try:
                    n = con.execute(f"SELECT count(*) FROM {tbl}").fetchone()[0]
                except Exception:  # noqa: BLE001
                    n = 0
                rows.append(((f"  ↳ {tbl}", "[green]ok[/]" if n else "[yellow]empty[/]"),
                             f"{n:,} rows" if n else f"not populated — run `{fix}`"))
    except Exception as e:  # noqa: BLE001 - never let doctor itself blow up
        return [(("variant coverage", "[yellow]unreadable[/]"), str(e)[:60])]
    return rows


@app.command()
def ingest(
    vcf_dir: Path = typer.Option(None, help="Directory of input VCFs (default: data/genome)."),
    normalize: bool = typer.Option(True, help="Left-align/split multiallelics (needs reference FASTA)."),
) -> None:
    """Index, QC, and normalize the sequencing.com VCFs."""
    from . import ingest as _ingest

    settings.ensure_dirs()
    _ingest.run(vcf_dir=vcf_dir or settings.genome_dir, normalize=normalize)


@app.command()
def download(
    target: str = typer.Argument(
        "all", help="reference|clinvar|snpeff|pharmcat|alphamissense|ancestry|haplogrep|gwas|all"),
) -> None:
    """Download & prepare reference and annotation databases."""
    from . import download as _download

    settings.ensure_dirs()
    _download.run(target)


@app.command()
def annotate(
    steps: str = typer.Option("all", help="Comma list: clinvar,snpeff,alphamissense,gnomad,pharmcat or 'all'."),
    force: bool = typer.Option(False, "--force",
                               help="Allow overwriting the store with FEWER annotations than it already has."),
) -> None:
    """Annotate variants against open-source databases."""
    from . import annotate as _annotate

    settings.ensure_dirs()
    _annotate.run(steps=steps, force=force)


@app.command()
def load() -> None:
    """Load annotated variants into the DuckDB store."""
    from . import load as _load

    settings.ensure_dirs()
    _load.run()


@app.command()
def ancestry() -> None:
    """Estimate biogeographic ancestry and ancestry-calibrated polygenic risk scores."""
    from . import ancestry as _ancestry
    from . import load as _load
    from . import pgs as _pgs

    settings.ensure_dirs()
    anc = _ancestry.run()
    scores = _pgs.run(nearest_superpop=anc.nearest)
    _load.write_ancestry(anc, scores)
    console.print("[green]Ancestry + polygenic scores written to the database.[/]")


@app.command()
def traits() -> None:
    """Genotype single-SNP traits/wellness (and the HLA-B*57:01 proxy) into the database."""
    from . import traits as _traits

    settings.ensure_dirs()
    _traits.run()


@app.command()
def gwas() -> None:
    """Genotype GWAS Catalog risk alleles (p<5e-8) and store the ones this genome carries."""
    from . import gwas as _gwas

    settings.ensure_dirs()
    _gwas.run()


@app.command()
def setup(
    skip_app: bool = typer.Option(False, "--skip-app", help="Don't build the macOS Dock app."),
    yes: bool = typer.Option(False, "--yes", "-y", help="Don't pause at the genome-file step."),
    skip_download: bool = typer.Option(False, "--skip-download", hidden=True),
) -> None:
    """Guided full install: download everything, build your genome store, register with Claude."""
    from . import setup as _setup

    settings.ensure_dirs()
    _setup.run(skip_app=skip_app, assume_yes=yes, skip_download=skip_download)


@app.command()
def pipeline(
    normalize: bool = typer.Option(True, help="Normalize during ingest."),
    steps: str = typer.Option("all", help="Annotation steps to run."),
) -> None:
    """Run ingest -> annotate -> load end to end."""
    from . import annotate as _annotate
    from . import ingest as _ingest
    from . import load as _load

    settings.ensure_dirs()
    _ingest.run(vcf_dir=settings.genome_dir, normalize=normalize)
    _annotate.run(steps=steps)
    _load.run()
    console.print("[green]Pipeline complete.[/] Try `locus serve mcp` or `locus serve api`.")


@app.command()
def refresh(
    sources: str = typer.Option("all", help="Comma list: clinvar,pgs,cpic,gwas,pubmed,litvar or 'all'."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Probe + report what would change; write nothing."),
    force: bool = typer.Option(False, "--force", help="Run the per-source work even if versions look unchanged."),
) -> None:
    """Check tracked sources for new releases and re-interpret what changed (ClinVar reanalysis,
    new GWAS associations at your variants, and new PubMed papers on your genes)."""
    from . import refresh as _refresh

    settings.ensure_dirs()
    _refresh.run(sources=sources, dry_run=dry_run, force=force)


@app.command()
def report(
    open_it: bool = typer.Option(True, "--open/--no-open", help="Open the report in your browser."),
) -> None:
    """Write a self-contained HTML summary of your genome (offline, shareable, no scripts)."""
    from . import report as _report

    settings.ensure_dirs()
    out = _report.build()
    if open_it:
        subprocess.run(["open", str(out)], check=False)


@app.command()
def literature(
    query: str = typer.Argument(..., help="A gene ('BRCA2'), an rsID ('rs7903146'), or a PubMed ID."),
    since: str = typer.Option("", help="Only papers on/after this ISO date (e.g. 2026-01-01)."),
) -> None:
    """Look up recent research: a gene/rsID → recent PubMed papers; a PubMed ID → which variants
    that study reported that THIS genome carries. Sends only gene symbols / rsIDs — never genotypes."""
    from . import literature as _lit

    settings.ensure_dirs()
    if query.strip().isdigit():  # a PubMed ID → paper → your variants
        res = _lit.study_variants(query.strip())
        console.print(f"[bold]PMID {res['pmid']}[/] — {res['total']} variant(s) reported, "
                      f"you carry {res['carried']}.")
        for m in res["markers"]:
            mark = "●" if _lit._is_carried(m.get("genotype"), m.get("ref")) else "○"
            console.print(f"  {mark} {m['rsid']} {m.get('genotype', '—')}"
                          + (f"  ({m['gene']})" if m.get("gene") else ""))
        console.print(f"[dim]{res['note']}[/]")
        return
    hits = _lit.literature_for(query, since=since or None)
    console.print(f"[bold]{len(hits)} recent paper(s)[/] for '{query}':")
    for h in hits:
        console.print(f"  • [bold]{h.title}[/] — {h.journal} {h.year}  [dim]{h.url}[/]")


@schedule_app.command("install")
def schedule_install(
    weekday: int = typer.Option(0, help="Day of week (0=Sun … 6=Sat)."),
    hour: int = typer.Option(3, help="Hour of day (0-23), local time."),
) -> None:
    """Install a weekly launchd job that runs `locus refresh`."""
    from . import schedule as _schedule

    settings.ensure_dirs()
    _schedule.install(weekday=weekday, hour=hour)


@schedule_app.command("uninstall")
def schedule_uninstall() -> None:
    """Remove the scheduled refresh job."""
    from . import schedule as _schedule

    _schedule.uninstall()


@schedule_app.command("status")
def schedule_status() -> None:
    """Show whether the refresh job is scheduled and loaded."""
    from . import schedule as _schedule

    _schedule.status()


@mcp_app.command("install")
def mcp_install_cmd() -> None:
    """Register the Locus MCP server with Claude Desktop and Claude Code (safe-merge + backup)."""
    from . import mcp_install

    settings.ensure_dirs()
    mcp_install.run()


@mcp_app.command("status")
def mcp_status_cmd() -> None:
    """Show whether the MCP server is registered with Claude (and points at this repo)."""
    from . import mcp_install

    mcp_install.status()


@serve_app.command("mcp")
def serve_mcp() -> None:
    """Start the MCP server so Claude can query your genome (stdio)."""
    from .mcp_server import main as mcp_main

    mcp_main()


@serve_app.command("api")
def serve_api(
    host: str = typer.Option(None, help="Bind host (default from config: localhost)."),
    port: int = typer.Option(None, help="Bind port."),
) -> None:
    """Start the FastAPI backend for the debug SPA."""
    import uvicorn

    uvicorn.run(
        "locus.api:app",
        host=host or settings.api_host,
        port=port or settings.api_port,
        reload=False,
    )


if __name__ == "__main__":
    app()

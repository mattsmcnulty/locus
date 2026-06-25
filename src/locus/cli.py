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

    def _db_row(label: str, present: bool, hint: str) -> None:
        table.add_row(label, "[green]ok[/]" if present else "[yellow]absent[/]", hint)

    ann = settings.annotations_dir
    _db_row("reference FASTA", artifacts.find_reference() is not None, "locus download reference")
    _db_row("ClinVar", (ann / download.CLINVAR_CHR_VCF).exists(), "locus download clinvar")
    _db_row("SnpEff", (ann / "snpEff" / "snpEff.jar").exists(), "locus download snpeff")
    _db_row("AlphaMissense", (ann / "alphamissense" / "AlphaMissense_hg38.slim.tsv.bgz").exists(),
            "locus download alphamissense")
    _db_row("PharmCAT (native)", (artifacts.pharmcat_install_dir() / "pharmcat_pipeline").exists(),
            "locus download pharmcat")
    _db_row("ancestry panel",
            (ann / "ancestry" / "all_hg38.pgen").exists() and (settings.data_dir / "tools" / "plink2").exists(),
            "locus download ancestry")
    _db_row("GWAS Catalog", (ann / "gwas" / "gwas-catalog-associations.tsv").exists(), "locus download gwas")
    _db_row("Haplogrep (mtDNA)", (ann / "haplogrep" / "haplogrep.jar").exists(), "locus download haplogrep")

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
    console.print(table)


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
    steps: str = typer.Option("all", help="Comma list: clinvar,gnomad,snpeff,pharmcat or 'all'."),
) -> None:
    """Annotate variants against open-source databases."""
    from . import annotate as _annotate

    settings.ensure_dirs()
    _annotate.run(steps=steps)


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
    sources: str = typer.Option("all", help="Comma list: clinvar,pgs or 'all'."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Probe + report what would change; write nothing."),
    force: bool = typer.Option(False, "--force", help="Run the per-source work even if versions look unchanged."),
) -> None:
    """Check tracked sources for new releases and re-interpret what changed (ClinVar reanalysis)."""
    from . import refresh as _refresh

    settings.ensure_dirs()
    _refresh.run(sources=sources, dry_run=dry_run, force=force)


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

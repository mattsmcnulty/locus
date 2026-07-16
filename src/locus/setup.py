"""Guided full install: download everything → locate the genome → build → interpret →
register with Claude → optional Dock app.

Built for non-engineers: every step is wrapped so a failure prints a plain-English
"paused, safe to re-run" panel instead of a Python traceback. All underlying steps
self-skip on re-run, so the whole thing is resumable.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel

from . import artifacts, shell
from .config import settings

console = Console()
REPO_ROOT = Path(__file__).resolve().parents[2]

EXAMPLE_QUESTIONS = [
    "Using my locus genome, give me an overview of what's loaded.",
    "What's my biogeographic ancestry, and where do my ancestors come from?",
    "Do I have any ACMG secondary findings I should know about?",
    "What does my genome say about caffeine, lactose, and alcohol?",
    "What are my pharmacogenomic results for clopidogrel and warfarin?",
    "What's my polygenic risk for coronary artery disease, as a percentile?",
]


def _step(title: str, fn: Callable[[], None]) -> None:
    console.rule(f"[bold]{title}")
    try:
        fn()
    except typer.Exit:
        raise
    except Exception as e:  # noqa: BLE001 - friendly, non-engineer-facing surface
        console.print(Panel(
            f"[red]'{title}' didn't finish.[/]\n\n{e}\n\n"
            "This is safe to retry — run [bold]./setup.command[/] again (or "
            "[bold]uv run locus setup[/]). Finished steps are skipped automatically.",
            title="Setup paused", border_style="red"))
        raise typer.Exit(code=1) from e


def _verify_deps() -> None:
    missing = [t for t in ("bcftools", "samtools", "tabix", "bgzip") if not shell.have(t)]
    if missing:
        raise RuntimeError(
            f"Missing required tools: {', '.join(missing)}.\n"
            "Run ./setup.command first — it installs them with Homebrew.")
    if shell.resolve_java() is None:
        console.print("[yellow]Java not found — SnpEff/PharmCAT/mtDNA haplogroup will be skipped. "
                      "Install with `brew install openjdk` and re-run to enable them.[/]")
    console.print("[green]Required tools present.[/]")


def _download_all() -> None:
    from . import download

    console.print("Downloading reference + databases (~9 GB; resumable if interrupted)…")
    download.run("all")


def _genome_gate(assume_yes: bool) -> None:
    from . import vcfutils

    while True:
        inputs = artifacts.classify_inputs(settings.genome_dir)
        if inputs.small_variants is not None:
            break
        console.print(Panel(
            "Export your genome from sequencing.com (log in → your 30× WGS order → Download) "
            f"and save the files into:\n\n  [bold]{settings.genome_dir}[/]\n\n"
            "Required: the file with [bold]snp-indel.genome.vcf.gz[/] in its name.\n"
            "Recommended too: the [bold]cnv.vcf.gz[/] and [bold]sv.vcf.gz[/] files.",
            title="Add your genome files", border_style="yellow"))
        if assume_yes or not sys.stdin.isatty():
            raise typer.Exit(code=2)
        if not typer.confirm("I've added the files — check again?", default=True):
            raise typer.Exit(code=2)
    # GRCh38 build check — a GRCh37/hg19 file runs end-to-end but mis-annotates everything.
    build = vcfutils.detect_build(vcfutils.read_info(inputs.small_variants))
    if build.startswith("GRCh37"):
        raise RuntimeError(
            f"{inputs.small_variants.name} looks like {build} (GRCh37/hg19). Locus only "
            "understands GRCh38 — please re-export the GRCh38 file from sequencing.com.")
    console.print(f"[green]Found your genome[/] {inputs.small_variants.name}  (build {build})")


def _build_db() -> None:
    from . import annotate, ingest, load

    ingest.run(vcf_dir=settings.genome_dir, normalize=True)
    annotate.run(steps="all")
    load.run()


def _ancestry() -> None:
    from . import ancestry, load, pgs

    anc = ancestry.run()
    load.write_ancestry(anc, pgs.run(nearest_superpop=anc.nearest))


def _traits() -> None:
    from . import traits

    traits.run()


def _gwas() -> None:
    from . import gwas

    gwas.run()


def _register_mcp() -> None:
    from . import mcp_install

    mcp_install.run()


def _final() -> None:
    qs = "\n".join(f"  • {q}" for q in EXAMPLE_QUESTIONS)
    console.print(Panel(
        "[bold green]Setup complete![/]\n\n"
        "1. Fully quit Claude Desktop (Cmd-Q, not just close the window) and reopen it.\n"
        "2. Then ask Claude (Desktop or Code) things like:\n"
        f"{qs}\n\n"
        "Prefer a visual UI? Open the [bold]Locus[/] app (Dock) or run [bold]locus serve api[/].\n\n"
        "[dim]Not medical advice — for research/education. Confirm anything health-relevant with a "
        "clinician or genetic counselor. Polygenic percentiles are estimates, valid only within a "
        "matched ancestry.[/]",
        title="🧬 Locus is ready", border_style="green"))


def run(*, skip_app: bool = False, assume_yes: bool = False, skip_download: bool = False) -> None:
    console.print(Panel(
        "[bold]Locus — full setup[/]\nDownloads ~9 GB of reference databases and builds your local "
        "genome store, then registers it with Claude. ~20–40 minutes, one-time. Your genome and your "
        "genotypes stay on your Mac; lookups send only public IDs (rsIDs, gene names).",
        border_style="cyan"))
    settings.ensure_dirs()

    _step("1/8  Checking tools", _verify_deps)
    if not skip_download:
        _step("2/8  Downloading databases (~9 GB)", _download_all)
    _step("3/8  Locating your genome", lambda: _genome_gate(assume_yes))
    _step("4/8  Building the genome database", _build_db)
    _step("5/8  Ancestry + polygenic risk", _ancestry)
    _step("6/8  Traits + maternal haplogroup", _traits)
    _step("7/8  GWAS associations you carry", _gwas)
    _step("8/8  Registering with Claude", _register_mcp)

    if not skip_app:
        console.rule("[bold]Optional: macOS Dock app")
        try:
            shell.run(["bash", str(REPO_ROOT / "scripts" / "build_macos_app.sh")])
        except Exception as e:  # noqa: BLE001 - the app is a nicety, never fatal
            console.print(f"[yellow]Dock app couldn't build ({e}). Skipping — everything else is set; "
                          "you can build it later with scripts/build_macos_app.sh.[/]")

    _final()

"""Local ancestry / chromosome painting (Phase 3) via Gnomix.

Produces per-segment ancestry along each chromosome (the 23andMe karyogram view).
This is the heaviest Locus capability and has two important caveats:

1. **Reference build.** Gnomix's *pretrained* models (AI-sandbox) are GRCh37, but
   Locus runs on GRCh38 — so we use Gnomix in *train* mode against the GRCh38
   HGDP+1KG phased panel (`locus download localancestry`). Training + inference is
   per-chromosome and takes hours for a whole genome (run it in the background).

2. **Value depends on admixture.** Local ancestry is informative for *admixed*
   genomes (the mosaic of segments). For a non-admixed individual it is essentially
   one ancestry throughout — a monochrome painting. Check `locus ancestry` first.

The surfacing (DuckDB ``ancestry_segments``, the ``ancestry_painting`` query/MCP
tool, and the SPA karyogram) works regardless; only this compute step is heavy.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rich.console import Console

from . import ancestry, shell
from .config import settings

console = Console()

# Autosomes only (local ancestry painting is conventionally autosomal).
CHROMS = [str(i) for i in range(1, 23)]


def _la_dir() -> Path:
    return settings.annotations_dir / "localancestry"


def _gnomix() -> Path:
    g = _la_dir() / "gnomix" / "gnomix.py"
    if not g.exists():
        raise FileNotFoundError("Gnomix not installed. Run `locus download localancestry`.")
    return g


def _panel_vcf(chrom: str) -> Path:
    return _la_dir() / "panel" / f"hgdp1kgp_chr{chrom}.shapeit5_phased.filter1_SNP_maf005.rechr.vcf.gz"


def _genetic_map(chrom: str) -> Path:
    return _la_dir() / "maps" / f"plink.chr{chrom}.GRCh38.map"


def _sample_map() -> Path:
    """Sample -> population map for the HGDP+1KG panel (built from the panel metadata)."""
    return _la_dir() / "hgdp1kgp_sample_map.tsv"


def _query_vcf(chrom: str) -> Path:
    """Matt's genotypes at the panel's chr positions (chr-prefixed, hom-ref aware)."""
    dest = settings.work_dir / f"query_chr{chrom}.vcf.gz"
    bed = settings.work_dir / f"panel_chr{chrom}.bed"
    shell.sh(
        f"bcftools query -f '%CHROM\\t%POS0\\t%END\\n' {_panel_vcf(chrom)} > {bed}"
    )
    ancestry.markers_genotypes(bed, dest)
    return dest


def run_chromosome(chrom: str) -> list[tuple]:
    """Train Gnomix on the GRCh38 panel for one chromosome and infer the sample's segments."""
    out_dir = settings.reports_dir / "localancestry" / f"chr{chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)
    query = _query_vcf(chrom)

    # python gnomix.py <query> <out> <chr> <phase=True> <genetic_map> <reference> <sample_map>
    import os

    env = dict(os.environ)
    java = shell.resolve_java()
    if java:
        env["PATH"] = f"{Path(java).parent}:{env.get('PATH', '')}"
    shell.run_env(
        [sys.executable, str(_gnomix()), str(query), str(out_dir), chrom, "True",
         str(_genetic_map(chrom)), str(_panel_vcf(chrom)), str(_sample_map())],
        env=env,
    )
    return _parse_msp(next(out_dir.glob("*.msp"), None) or next(out_dir.glob("*.msp.tsv")), chrom)


def _parse_msp(msp: Path | None, chrom: str) -> list[tuple]:
    """Parse a Gnomix .msp into (haplotype, chrom, start, end, ancestry, posterior) rows."""
    if msp is None or not msp.exists():
        return []
    rows: list[tuple] = []
    labels: dict[int, str] = {}
    with open(msp) as fh:
        for line in fh:
            if line.startswith("#Subpopulation") or "=" in line and line.startswith("#"):
                # header maps numeric code -> population name
                for tok in line.replace("#Subpopulation order/codes:", "").split():
                    if "=" in tok:
                        name, code = tok.split("=")
                        labels[int(code)] = name
                continue
            if line.startswith("#") or line.startswith("chm"):
                continue
            c = line.rstrip("\n").split("\t")
            spos, epos = int(c[1]), int(c[2])
            # columns 6+ are per-haplotype ancestry codes
            for hap, code in enumerate(c[6:8]):
                anc = labels.get(int(code), code)
                rows.append((hap, f"chr{chrom}", spos, epos, anc, None))
    return rows


def run(chroms: list[str] | None = None) -> int:
    """Run local ancestry across chromosomes and write the painting to the DB. Heavy — see module docs."""
    from .load import write_segments

    if not _gnomix().exists():  # raises if missing
        return 0
    targets = chroms or CHROMS
    console.rule("[bold]Local ancestry (chromosome painting)")
    console.print(f"[yellow]Heavy: per-chromosome Gnomix training on the GRCh38 panel ({len(targets)} chroms).[/]")
    segments: list[tuple] = []
    for chrom in targets:
        if not _panel_vcf(chrom).exists():
            console.print(f"[yellow]chr{chrom}: panel missing — run `locus download localancestry`.[/]")
            continue
        console.print(f"chr{chrom}…")
        segments.extend(run_chromosome(chrom))
    write_segments(segments)
    console.print(f"[green]Painting written[/] — {len(segments)} segments across {len(targets)} chromosomes.")
    return len(segments)

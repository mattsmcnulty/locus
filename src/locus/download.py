"""Download & prepare reference and annotation databases.

Each function is idempotent (skips work already done) and encodes the exact
sources / prep steps verified in docs/integration-notes.md. Big downloads (VEP
cache ~25 GB, dbSNP ~28 GB, gnomAD ~hundreds of GB) are gated behind explicit
opt-in; the defaults here are the small, high-value ones (reference, ClinVar)
plus the native-friendly SnpEff and the PharmCAT Docker image.
"""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from . import shell
from .config import settings
from .vcfutils import plain_to_chr_map, write_rename_file

console = Console()

REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
    "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
SNPEFF_URL = "https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"
SNPEFF_DB = "GRCh38.mane.1.0.ensembl"

# PharmCAT runs natively (no Docker): a Java jar + a Python preprocessor.
PHARMCAT_VERSION = "3.2.0"
PHARMCAT_PIPELINE_URL = (
    f"https://github.com/PharmGKB/PharmCAT/releases/download/v{PHARMCAT_VERSION}/"
    f"pharmcat-pipeline-{PHARMCAT_VERSION}.tar.gz"
)
PHARMCAT_DEPS = ["colorama>=0.4.6", "pandas>=2.1.3", "packaging~=24.1"]

REFERENCE_FASTA = "GRCh38_no_alt.fa.gz"
CLINVAR_CHR_VCF = "clinvar.chr.vcf.gz"


def _curl(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    # -C - resumes partial downloads; -L follows redirects; --fail surfaces HTTP errors.
    shell.run(["curl", "-fL", "-C", "-", "-o", str(dest), url])


def download_reference() -> Path:
    """GRCh38 no-alt analysis set → bgzip → faidx. ~873 MB download."""
    out = settings.reference_dir / REFERENCE_FASTA
    if out.exists() and Path(str(out) + ".fai").exists():
        console.print(f"[green]reference present[/] → {out}")
        return out
    settings.reference_dir.mkdir(parents=True, exist_ok=True)
    raw = settings.reference_dir / "GRCh38_no_alt.fna.gz"
    console.print("[bold]Downloading GRCh38 no-alt analysis set (~873 MB)…[/]")
    _curl(REFERENCE_URL, raw)
    # NCBI ships PLAIN gzip — must re-bgzip before faidx (see integration-notes).
    console.print("Re-compressing as bgzip and indexing…")
    plain_fa = settings.reference_dir / "GRCh38_no_alt.fa"
    shell.sh(f"gzip -dc {raw} > {plain_fa}")
    shell.run(["bgzip", "-f", str(plain_fa)])  # -> GRCh38_no_alt.fa.gz
    shell.run(["samtools", "faidx", str(out)])
    raw.unlink(missing_ok=True)
    console.print(f"[green]reference ready[/] → {out}")
    return out


def download_clinvar() -> Path:
    """ClinVar GRCh38 VCF (~183 MB), md5-checked, renamed to chr-prefixed contigs."""
    out = settings.annotations_dir / CLINVAR_CHR_VCF
    if out.exists():
        console.print(f"[green]ClinVar present[/] → {out}")
        return out
    ann = settings.annotations_dir
    ann.mkdir(parents=True, exist_ok=True)
    raw = ann / "clinvar.vcf.gz"
    console.print("[bold]Downloading ClinVar GRCh38 VCF (~183 MB)…[/]")
    _curl(f"{CLINVAR_BASE}/clinvar.vcf.gz", raw)
    _curl(f"{CLINVAR_BASE}/clinvar.vcf.gz.tbi", ann / "clinvar.vcf.gz.tbi")
    _curl(f"{CLINVAR_BASE}/clinvar.vcf.gz.md5", ann / "clinvar.vcf.gz.md5")
    # Verify md5 (file format: "<md5>  clinvar.vcf.gz").
    expected = (ann / "clinvar.vcf.gz.md5").read_text().split()[0]
    actual = shell.capture(["bash", "-c", f"md5 -q {raw} 2>/dev/null || md5sum {raw} | cut -d' ' -f1"]).strip()
    if expected and actual and expected != actual:
        raise shell.ToolError(f"ClinVar md5 mismatch: expected {expected}, got {actual}")
    console.print("md5 OK. Renaming contigs to chr-prefixed (1→chr1, MT→chrM)…")
    rename = write_rename_file(plain_to_chr_map(), ann / "clinvar2chr.txt")
    shell.run(["bcftools", "annotate", "--rename-chrs", str(rename), str(raw), "-Oz", "-o", str(out)])
    shell.run(["bcftools", "index", "-t", str(out)])
    raw.unlink(missing_ok=True)
    console.print(f"[green]ClinVar ready[/] → {out}")
    return out


def setup_snpeff() -> Path:
    """Download the SnpEff jar (pure Java, native arm64) + the GRCh38 database."""
    snpeff_dir = settings.annotations_dir / "snpEff"
    jar = snpeff_dir / "snpEff.jar"
    if jar.exists():
        console.print(f"[green]SnpEff present[/] → {jar}")
    else:
        settings.annotations_dir.mkdir(parents=True, exist_ok=True)
        zip_path = settings.annotations_dir / "snpEff_latest_core.zip"
        console.print("[bold]Downloading SnpEff core (~50 MB)…[/]")
        _curl(SNPEFF_URL, zip_path)
        shell.run(["unzip", "-o", "-q", str(zip_path), "-d", str(settings.annotations_dir)])
        zip_path.unlink(missing_ok=True)
    # Download the GRCh38 database into snpEff/data.
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — install it with `brew install openjdk`, then re-run.[/]")
        return jar
    console.print(f"Fetching SnpEff database {SNPEFF_DB} (~0.5 GB, one-time)…")
    shell.run(shell.java_cmd(["-jar", str(jar), "download", "-v", SNPEFF_DB]))
    console.print(f"[green]SnpEff ready[/] (db {SNPEFF_DB})")
    return jar


def setup_pharmcat() -> Path:
    """Install PharmCAT natively (no Docker): download the pipeline tarball + Python deps.

    The tarball bundles the jar, the ``pharmcat_pipeline`` wrapper, the
    ``pharmcat_vcf_preprocessor``, and the PGx positions/regions files.
    """
    install_dir = settings.annotations_dir / "pharmcat"
    pipeline = install_dir / "pharmcat_pipeline"
    if not pipeline.exists():
        install_dir.mkdir(parents=True, exist_ok=True)
        tarball = settings.annotations_dir / "pharmcat-pipeline.tar.gz"
        console.print(f"[bold]Downloading PharmCAT {PHARMCAT_VERSION} pipeline (~28 MB)…[/]")
        _curl(PHARMCAT_PIPELINE_URL, tarball)
        shell.run(["tar", "-xzf", str(tarball), "-C", str(install_dir)])
        tarball.unlink(missing_ok=True)

    # Install the preprocessor's Python deps into the active environment.
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — PharmCAT needs Java 17+. `brew install openjdk`.[/]")
    console.print("Installing PharmCAT preprocessor Python deps (colorama, pandas, packaging)…")
    shell.run(["uv", "pip", "install", "-q", *PHARMCAT_DEPS])
    console.print(f"[green]PharmCAT ready[/] (native, v{PHARMCAT_VERSION}) → {install_dir}")
    return pipeline


# Targets that require an explicit, deliberate opt-in (very large).
def guidance_large() -> None:
    console.print(
        "\n[bold]Large optional databases[/] (opt-in, see docs/integration-notes.md):\n"
        "  • gnomAD v4.1 — streamed per-region during `locus annotate` (no full download).\n"
        "  • dbSNP 157 (~28 GB) — usually unnecessary; DRAGEN already fills rsIDs.\n"
        "  • VEP offline cache (~25 GB) — richer consequences; run via Docker (ensemblorg/ensembl-vep).\n"
    )


TARGETS = {
    "reference": download_reference,
    "clinvar": download_clinvar,
    "snpeff": setup_snpeff,
    "pharmcat": setup_pharmcat,
}


def run(target: str) -> None:
    if target == "all":
        # Resilient: a failure in one step shouldn't abort the others.
        for name, fn in TARGETS.items():
            try:
                fn()
            except Exception as e:  # noqa: BLE001 - report and continue
                console.print(f"[red]'{name}' failed:[/] {e}\n[dim]Re-run later with `locus download {name}`.[/]")
        guidance_large()
        return
    fn = TARGETS.get(target)
    if not fn:
        raise SystemExit(f"Unknown target '{target}'. Choose from: {', '.join(TARGETS)}, all.")
    fn()

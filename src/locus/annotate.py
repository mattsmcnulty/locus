"""Annotate phase: layer open-source databases onto the sites VCF.

Each step is independent and self-skips when its database isn't present, so you
can run a subset (``locus annotate --steps clinvar,pharmcat``) or ``all``. The
small-variant chain is successive ``bcftools annotate`` calls; PharmCAT runs
separately (it needs the gVCF, not the sites VCF) and writes a JSON report that
the load phase ingests.
"""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from . import artifacts, download, shell
from .config import settings
from .vcfutils import read_info

console = Console()

# ClinVar fields to transfer (curated; not the whole INFO block).
CLINVAR_FIELDS = (
    "INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC,INFO/CLNDISDB,"
    "INFO/CLNHGVS,INFO/MC,INFO/ALLELEID,INFO/GENEINFO"
)
# gnomAD AF fields, prefixed so they don't collide with any existing AF.
GNOMAD_TRANSFER = (
    "INFO/gnomAD_AF:=INFO/AF,INFO/gnomAD_AF_grpmax:=INFO/AF_grpmax,"
    "INFO/gnomAD_grpmax:=INFO/grpmax,INFO/gnomAD_AC:=INFO/AC,INFO/gnomAD_AN:=INFO/AN"
)
GNOMAD_GENOMES = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/"
    "gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz"
)

ALL_STEPS = ["clinvar", "gnomad", "snpeff", "alphamissense", "pharmcat"]


def _index(vcf: Path) -> None:
    shell.run(["bcftools", "index", "-f", "-t", str(vcf)])


def _present_info(vcf: Path) -> set[str]:
    """INFO field IDs defined in a VCF header (so we don't request absent fields)."""
    header = shell.capture(["bcftools", "view", "-h", str(vcf)])
    return {
        ln.split("ID=", 1)[1].split(",", 1)[0]
        for ln in header.splitlines()
        if ln.startswith("##INFO=<ID=")
    }


def annotate_clinvar(src: Path, dest: Path) -> Path:
    clinvar = settings.annotations_dir / download.CLINVAR_CHR_VCF
    if not clinvar.exists():
        console.print("[yellow]ClinVar DB missing — skipping. Run `locus download clinvar`.[/]")
        return src
    # Only transfer fields actually present in this ClinVar build.
    present = _present_info(clinvar)
    fields = [f for f in CLINVAR_FIELDS.split(",") if f.split("/")[-1] in present]
    if not fields:
        console.print("[yellow]No expected ClinVar INFO fields found — skipping.[/]")
        return src
    console.print("Annotating ClinVar clinical significance…")
    shell.run([
        "bcftools", "annotate", "-a", str(clinvar), "-c", ",".join(fields),
        str(src), "-Oz", "-o", str(dest),
    ])
    _index(dest)
    # Validate the join actually matched something (0 ⇒ contig mismatch — see notes).
    n = shell.capture(["bash", "-o", "pipefail", "-c", f"bcftools view -H {dest} | grep -c CLNSIG || true"]).strip()
    console.print(f"  ClinVar-annotated records: {n}")
    return dest


def annotate_gnomad(src: Path, dest: Path) -> Path:
    """Stream gnomAD v4.1 genome AF per-chromosome (no full download) and transfer AF fields."""
    console.print("Annotating gnomAD v4.1 allele frequencies (streamed; can be slow over network)…")
    info = read_info(src)
    chroms = [c for c in info.contigs if c.startswith("chr")] or sorted(
        {ln.split("\t")[0] for ln in shell.capture(
            ["bash", "-o", "pipefail", "-c", f"bcftools view -H {src} | cut -f1 | sort -u"]).splitlines()}
    )
    work = settings.work_dir
    per_chrom: list[Path] = []
    try:
        for chrom in chroms:
            regions = work / f"{chrom}.regions.bed"
            shell.sh(f"bcftools query -f '%CHROM\\t%POS0\\t%END\\n' -r {chrom} {src} > {regions}")
            if regions.stat().st_size == 0:
                continue
            slim = work / f"gnomad.{chrom}.slim.vcf.gz"
            url = GNOMAD_GENOMES.format(chrom=chrom)
            shell.sh(
                f"bcftools view -R {regions} '{url}' "
                f"| bcftools annotate -x '^INFO/AF,INFO/AC,INFO/AN,INFO/AF_grpmax,INFO/grpmax' "
                f"-Oz -o {slim}"
            )
            _index(slim)
            out_c = work / f"{chrom}.gnomad.vcf.gz"
            shell.run([
                "bcftools", "annotate", "-a", str(slim), "-c", GNOMAD_TRANSFER,
                "-r", chrom, str(src), "-Oz", "-o", str(out_c),
            ])
            _index(out_c)
            per_chrom.append(out_c)
        if per_chrom:
            shell.run(["bcftools", "concat", "-Oz", "-o", str(dest), *map(str, per_chrom)])
            _index(dest)
            return dest
    except shell.ToolError as e:
        console.print(f"[yellow]gnomAD streaming failed ({e}); leaving variants without AF.[/]")
    return src


def annotate_snpeff(src: Path, dest: Path) -> Path:
    jar = settings.annotations_dir / "snpEff" / "snpEff.jar"
    if not jar.exists():
        console.print("[yellow]SnpEff missing — skipping consequences. Run `locus download snpeff`.[/]")
        return src
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — skipping SnpEff. `brew install openjdk`.[/]")
        return src
    console.print("Annotating functional consequences (SnpEff)…")
    # snpEff writes uncompressed VCF to stdout; bgzip it.
    cmd = " ".join(shell.java_cmd(["-Xmx6g", "-jar", str(jar), "-noStats", download.SNPEFF_DB, str(src)]))
    shell.sh(f"{cmd} | bgzip -c > {dest}")
    _index(dest)
    return dest


def annotate_alphamissense(src: Path, dest: Path) -> Path:
    """Annotate missense pathogenicity from AlphaMissense (calibrated score for ~every missense).

    Fills the gap where ClinVar is silent: a high am_pathogenicity on a variant ClinVar has never
    seen is real signal, not 'nothing'.
    """
    am = settings.annotations_dir / "alphamissense" / "AlphaMissense_hg38.slim.tsv.bgz"
    if not am.exists():
        console.print("[yellow]AlphaMissense missing — skipping. Run `locus download alphamissense`.[/]")
        return src
    header = settings.work_dir / "am.header.txt"
    header.write_text(
        '##INFO=<ID=am_pathogenicity,Number=1,Type=Float,Description="AlphaMissense pathogenicity (0-1)">\n'
        '##INFO=<ID=am_class,Number=1,Type=String,Description="AlphaMissense class (benign/ambiguous/pathogenic)">\n'
    )
    console.print("Annotating AlphaMissense missense pathogenicity…")
    shell.run([
        "bcftools", "annotate", "-a", str(am), "-h", str(header),
        "-c", "CHROM,POS,REF,ALT,am_pathogenicity,am_class",
        str(src), "-Oz", "-o", str(dest),
    ])
    _index(dest)
    n = shell.capture(["bash", "-c", f"bcftools view -H {dest} 2>/dev/null | grep -c am_pathogenicity || true"]).strip()
    console.print(f"  AlphaMissense-annotated records: {n}")
    return dest


def _pharmcat_input(inputs, reference: Path) -> Path:
    """Build PharmCAT's input VCF: PGx regions only, chr-prefixed, hom-ref blocks expanded.

    Restrict the gVCF to PharmCAT's regions first (small subset), rename contigs to
    chr-prefixed, then ``bcftools convert --gvcf2vcf`` to expand ``ALT="."`` hom-ref
    blocks into explicit per-position ``0/0`` calls — so PGx sites that are reference
    get genotyped instead of becoming no-calls.
    """
    install = artifacts.pharmcat_install_dir()
    regions = install / "pharmcat_regions.bed"
    pgx_input = artifacts.pharmcat_input_vcf()
    info = read_info(inputs.small_variants)

    if info.chr_prefixed:
        region_src = f"bcftools view -R {regions} {inputs.small_variants}"
    else:
        # Match the gVCF's non-chr contigs (chrM→MT), restrict, then rename to chr.
        nochr = settings.work_dir / "pharmcat_regions.nochr.bed"
        shell.sh(f"sed -e 's/^chrM\\t/MT\\t/' -e 's/^chr//' {regions} > {nochr}")
        rename = settings.work_dir / "contigs2chr.txt"
        if not rename.exists():
            from .vcfutils import chr_rename_map, write_rename_file
            write_rename_file(chr_rename_map(info.contigs), rename)
        region_src = f"bcftools view -R {nochr} {inputs.small_variants} | bcftools annotate --rename-chrs {rename}"

    # GT-only: PharmCAT only needs genotypes, and keeping AD/DP makes its internal
    # `bcftools norm` fail to merge per-allele FORMAT tags ("could not merge AD").
    shell.sh(
        f"{region_src} | bcftools convert --gvcf2vcf -f {reference} -Ou "
        f"| bcftools annotate -x '^FORMAT/GT' -Oz -o {pgx_input}"
    )
    _index(pgx_input)
    return pgx_input


def annotate_pharmcat() -> Path | None:
    """Run PharmCAT natively (jar + Python preprocessor) on the PGx-restricted VCF."""
    pipeline = artifacts.pharmcat_install_dir() / "pharmcat_pipeline"
    if not pipeline.exists():
        console.print("[yellow]PharmCAT not installed — run `locus download pharmcat`.[/]")
        return None
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — PharmCAT needs Java 17+. `brew install openjdk`.[/]")
        return None
    inputs = artifacts.classify_inputs(settings.genome_dir)
    if not inputs.small_variants:
        console.print("[yellow]No small-variant gVCF for PharmCAT.[/]")
        return None
    reference = artifacts.find_reference()
    if reference is None:
        console.print("[yellow]Reference FASTA needed for PharmCAT. Run `locus download reference`.[/]")
        return None

    out_dir = artifacts.pharmcat_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    console.print("Preparing PharmCAT input (PGx regions, chr-prefixed, hom-ref expanded)…")
    pgx_input = _pharmcat_input(inputs, reference)

    console.print("Running PharmCAT (native)…")
    # PharmCAT calls `java`; put the working JDK first on PATH for the subprocess.
    import os
    import sys

    env = dict(os.environ)
    java = shell.resolve_java()
    if java:
        env["PATH"] = f"{Path(java).parent}:{env.get('PATH', '')}"
        env["JAVA_HOME"] = str(Path(java).parent.parent)
    shell.run_env(
        [sys.executable, str(pipeline), str(pgx_input), "-o", str(out_dir),
         "-reporterJson", "-reporterCallsOnlyTsv"],
        env=env,
    )
    report = next(out_dir.glob("*.report.json"), None)
    console.print(f"[green]PharmCAT report:[/] {report}" if report else "[yellow]No PharmCAT report produced.[/]")
    return report


def run(steps: str = "all") -> Path:
    """Run the requested annotation steps, producing the annotated sites VCF."""
    src = artifacts.sites_vcf()
    if not src.exists():
        raise FileNotFoundError(f"No sites VCF ({src}). Run `locus ingest` first.")

    requested = ALL_STEPS if steps in ("all", "") else [s.strip() for s in steps.split(",")]
    console.rule(f"[bold]Annotate ({', '.join(requested)})")

    # bcftools annotate needs the input indexed — self-heal in case ingest's index is missing.
    from .ingest import ensure_index as _ensure_index

    _ensure_index(src)

    work = settings.work_dir
    cur = src
    if "clinvar" in requested:
        cur = annotate_clinvar(cur, work / f"{settings.sample_id}.clinvar.vcf.gz")
    if "gnomad" in requested:
        cur = annotate_gnomad(cur, work / f"{settings.sample_id}.gnomad.vcf.gz")
    if "snpeff" in requested or "vep" in requested:
        cur = annotate_snpeff(cur, work / f"{settings.sample_id}.snpeff.vcf.gz")
    if "alphamissense" in requested:
        cur = annotate_alphamissense(cur, work / f"{settings.sample_id}.am.vcf.gz")

    # Finalize the small-variant annotated VCF.
    dest = artifacts.annotated_vcf()
    if cur != src:
        shell.run(["bcftools", "view", str(cur), "-Oz", "-o", str(dest)])
        _index(dest)
        console.print(f"[green]Annotated VCF ready[/] → {dest}")
    else:
        console.print("[yellow]No small-variant annotations applied (no DBs present).[/]")

    if "pharmcat" in requested:
        annotate_pharmcat()

    return dest if dest.exists() else src

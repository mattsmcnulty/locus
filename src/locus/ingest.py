"""Ingest phase: index, QC, and normalize the sequencing.com / DRAGEN VCFs.

The headline transform is converting the DRAGEN small-variant **gVCF**
(``*.snp-indel.genome.vcf.gz``, which carries ``<NON_REF>`` blocks) into a
normalized **sites** VCF that the rest of the pipeline can treat as ordinary
variants. See docs/integration-notes.md for the why.
"""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from . import artifacts, shell
from .config import settings
from .shell import ToolError
from .vcfutils import VcfInfo, chr_rename_map, detect_build, read_info, write_rename_file

console = Console()


def ensure_index(vcf: Path) -> None:
    """Make sure a bgzipped VCF has a tabix index (verifies it was actually created)."""
    if Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists():
        return
    shell.run(["bcftools", "index", "-t", str(vcf)])
    if not (Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists()):
        raise ToolError(f"Failed to create a tabix index for {vcf}")


def qc_summary(vcf: Path) -> None:
    """Print a quick `bcftools stats` sanity check (counts, ts/tv, het/hom)."""
    text = shell.capture(["bash", "-o", "pipefail", "-c", f"bcftools stats {vcf} | grep -E '^SN'"])
    console.print(f"[bold]QC ({vcf.name}):[/]")
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:
            console.print(f"  {parts[2]:<35} {parts[3]}")


def gvcf_to_sites(
    gvcf: Path, reference: Path | None, dest: Path, *, normalize: bool, rename_file: Path | None = None
) -> Path:
    """Convert a gVCF to a normalized, chr-prefixed sites VCF.

    Handles both gVCF conventions: DRAGEN/GATK append a symbolic ``<NON_REF>``
    allele (real variants look like ``A  G,<NON_REF>``); bcftools/mpileup gVCFs
    use ``ALT="."`` + ``INFO/END`` hom-ref blocks. The pipeline:

    0. (optional) ``bcftools annotate --rename-chrs`` — canonicalize contigs to
       UCSC ``chr`` names so they match the reference / ClinVar / gnomAD.
    1. ``bcftools norm -m -any`` — split multiallelics; with ``-f`` also
       left-aligns indels so positions/alleles match the annotation DBs.
    2. ``bcftools view -e 'ALT="." || ALT="<NON_REF>"'`` — drop hom-ref block
       rows and split-off symbolic alleles, keeping only real variants.
    3. ``bcftools annotate -x INFO/END`` — strip the gVCF block-end tag.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    source = f"bcftools annotate --rename-chrs {rename_file} {gvcf}" if rename_file else f"bcftools view {gvcf}"
    if normalize and reference is not None:
        norm = f"bcftools norm -f {reference} -m -any -"
    else:
        if normalize and reference is None:
            console.print(
                "[yellow]No reference FASTA found — splitting multiallelics but NOT left-aligning. "
                "Run `locus download reference` for full normalization.[/]"
            )
        norm = "bcftools norm -m -any -"
    drop_blocks = "bcftools view -e 'ALT=\".\" || ALT=\"<NON_REF>\"'"
    strip_end = f"bcftools annotate -x INFO/END -Oz -o {dest}"
    shell.sh(f"{source} | {norm} | {drop_blocks} | {strip_end}")
    ensure_index(dest)
    return dest


def run(vcf_dir: Path, *, normalize: bool = True) -> VcfInfo:
    """Classify inputs, index them, QC + convert the small-variant gVCF to sites."""
    inputs = artifacts.classify_inputs(vcf_dir)
    if inputs.small_variants is None:
        raise ToolError(
            f"No small-variant VCF found in {vcf_dir}. Expected a *.snp-indel.genome.vcf.gz file."
        )

    console.rule("[bold]Ingest")
    console.print(f"small variants : {inputs.small_variants.name}")
    console.print(f"CNV            : {inputs.cnv.name if inputs.cnv else '—'}")
    console.print(f"SV             : {inputs.sv.name if inputs.sv else '—'}")

    # Index everything present.
    for vcf in filter(None, [inputs.small_variants, inputs.cnv, inputs.sv, *inputs.others]):
        ensure_index(vcf)

    info = read_info(inputs.small_variants)
    build = detect_build(info)
    console.print(
        f"\nbuild={build}  chr_prefixed={info.chr_prefixed}  gVCF={info.is_gvcf}  rsIDs={info.has_rsids}"
    )
    if build.startswith("GRCh37"):
        console.print("[red]This looks like GRCh37/hg19. Locus targets GRCh38 — annotations will be wrong.[/]")

    # Canonicalize contigs to chr-prefixed so they match the reference / ClinVar / gnomAD.
    rename_file = None
    if not info.chr_prefixed:
        mapping = chr_rename_map(info.contigs)
        rename_file = write_rename_file(mapping, settings.work_dir / "contigs2chr.txt")
        console.print(
            f"[cyan]Contigs are non-chr-prefixed (e.g. '1', 'MT') — renaming {len(mapping)} contigs "
            f"to chr-prefixed (1→chr1, MT→chrM) to match the reference & annotation DBs.[/]"
        )

    qc_summary(inputs.small_variants)

    reference = artifacts.find_reference()
    if reference is None and normalize:
        console.print("[yellow]Reference FASTA not present; see `locus download reference`.[/]")

    dest = artifacts.sites_vcf()
    console.print(f"\nConverting gVCF → normalized sites VCF: {dest}")
    gvcf_to_sites(inputs.small_variants, reference, dest, normalize=normalize, rename_file=rename_file)

    n = _sites_count(dest)
    console.print(f"[green]Sites VCF ready[/] — {n:,} variant records → {dest}")
    return info


def _sites_count(vcf: Path) -> int:
    return int(
        shell.capture(["bash", "-o", "pipefail", "-c", f"bcftools view -H {vcf} | wc -l"]).strip() or "0"
    )

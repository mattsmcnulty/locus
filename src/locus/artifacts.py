"""Canonical paths for pipeline inputs and intermediates.

Keeping these in one place lets ingest -> annotate -> load agree on filenames
without threading paths through every call.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .config import settings


@dataclass(frozen=True)
class GenomeInputs:
    """The classified sequencing.com / DRAGEN files found in the genome dir."""

    small_variants: Path | None  # *.snp-indel.genome.vcf.gz (a gVCF)
    cnv: Path | None             # *.cnv.vcf.gz
    sv: Path | None              # *.sv.vcf.gz
    others: list[Path]


def classify_inputs(vcf_dir: Path) -> GenomeInputs:
    """Sort the VCFs in a directory into small-variant / CNV / SV by DRAGEN naming."""
    small = cnv = sv = None
    others: list[Path] = []
    for p in sorted(vcf_dir.glob("*.vcf.gz")):
        name = p.name.lower()
        if "cnv" in name:
            cnv = p
        elif name.endswith(".sv.vcf.gz") or ".sv." in name:
            sv = p
        elif "snp-indel" in name or "genome.vcf" in name or "hard-filtered" in name:
            small = small or p
        else:
            others.append(p)
    # If nothing matched the small-variant pattern, take the first non-CNV/SV file.
    if small is None:
        for p in others:
            small = p
            others = [o for o in others if o != p]
            break
    return GenomeInputs(small_variants=small, cnv=cnv, sv=sv, others=others)


def find_reference() -> Path | None:
    """Locate the prepared GRCh38 reference FASTA (bgzipped, faidx-indexed)."""
    rd = settings.reference_dir
    if not rd.exists():
        return None
    # Prefer a bgzipped, indexed FASTA.
    for fa in sorted(rd.glob("*.fa.gz")) + sorted(rd.glob("*.fna.gz")) + sorted(rd.glob("*.fasta.gz")):
        if fa.with_suffix(fa.suffix + ".fai").exists() or Path(str(fa) + ".fai").exists():
            return fa
    # Fall back to any plain FASTA with a .fai.
    for fa in sorted(rd.glob("*.fa")) + sorted(rd.glob("*.fna")) + sorted(rd.glob("*.fasta")):
        if Path(str(fa) + ".fai").exists():
            return fa
    return None


# Intermediate artifacts (under data/work).
def sites_vcf() -> Path:
    """Normalized sites VCF derived from the DRAGEN gVCF (no annotations yet)."""
    return settings.work_dir / f"{settings.sample_id}.sites.vcf.gz"


def annotated_vcf() -> Path:
    """Sites VCF after ClinVar / dbSNP / gnomAD / VEP annotation."""
    return settings.work_dir / f"{settings.sample_id}.annotated.vcf.gz"


def pharmcat_input_vcf() -> Path:
    """gVCF expanded to a regular VCF (REF blocks kept) for PharmCAT."""
    return settings.work_dir / f"{settings.sample_id}.pharmcat-input.vcf.gz"


def pharmcat_dir() -> Path:
    """Where PharmCAT reports are written."""
    return settings.reports_dir / "pharmcat"


def pharmcat_install_dir() -> Path:
    """Where the native PharmCAT pipeline (jar + preprocessor) is installed."""
    return settings.annotations_dir / "pharmcat"

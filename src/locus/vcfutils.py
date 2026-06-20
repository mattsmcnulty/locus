"""VCF introspection helpers.

These guard the two failure modes the research flagged as *silent*:
1. contig-naming mismatches (chr-prefixed vs not) that make annotation match zero records;
2. mistaking a DRAGEN gVCF (with ``<NON_REF>`` blocks) for a normal sites VCF.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .shell import capture

# Canonical chromosome stems used to build rename maps (1..22, X, Y, MT/M).
_AUTOSOMES = [str(i) for i in range(1, 23)]


@dataclass
class VcfInfo:
    path: Path
    contigs: list[str]
    chr_prefixed: bool          # do data contigs look like "chr1"?
    is_gvcf: bool               # contains <NON_REF> symbolic allele?
    reference: str | None       # ##reference header value, if any
    has_rsids: bool             # any non-"." ID seen in the sampled records?


def _header(path: Path) -> str:
    return capture(["bcftools", "view", "-h", str(path)])


def read_info(path: Path, sample_records: int = 2000) -> VcfInfo:
    """Introspect a VCF/gVCF without loading it fully."""
    header = _header(path)

    contigs = [
        line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">")
        for line in header.splitlines()
        if line.startswith("##contig=")
    ]
    # ##contig lines can be absent (e.g. ClinVar) — fall back to scanning CHROM later.
    # gVCF block conventions vary: DRAGEN/GATK use a <NON_REF> symbolic allele;
    # bcftools/mpileup gVCFs use ALT="." + INFO/END hom-ref blocks. Detect both.
    is_gvcf = "##ALT=<ID=NON_REF" in header or "ID=<NON_REF" in header

    reference = None
    for line in header.splitlines():
        if line.startswith("##reference="):
            reference = line.split("=", 1)[1].strip()
            break

    # Sample a few data rows for contig style, rsID presence, and hom-ref blocks.
    chr_prefixed = any(c.startswith("chr") for c in contigs)
    has_rsids = False
    # No pipefail here: `head` closing the pipe SIGPIPEs bcftools by design — expected, not an error.
    rows = capture(
        ["bash", "-c", f"bcftools view -H {path} 2>/dev/null | head -n {sample_records}"]
    )
    seen_chroms: set[str] = set()
    for row in rows.splitlines():
        cols = row.split("\t")
        if len(cols) < 8:
            continue
        seen_chroms.add(cols[0])
        if cols[2] not in (".", ""):
            has_rsids = True
        if cols[4] in (".", "<NON_REF>") or "END=" in cols[7]:
            is_gvcf = True  # reference-block row → this is a gVCF
    if seen_chroms:
        chr_prefixed = any(c.startswith("chr") for c in seen_chroms)
        if not contigs:
            contigs = sorted(seen_chroms)

    return VcfInfo(
        path=Path(path),
        contigs=contigs,
        chr_prefixed=chr_prefixed,
        is_gvcf=is_gvcf,
        reference=reference,
        has_rsids=has_rsids,
    )


def plain_to_chr_map() -> dict[str, str]:
    """Map Ensembl/ClinVar-style names (1, X, MT) -> UCSC chr names (chr1, chrX, chrM)."""
    m = {c: f"chr{c}" for c in _AUTOSOMES}
    m.update({"X": "chrX", "Y": "chrY", "MT": "chrM"})
    return m


def canonical_chrom(contig: str) -> str:
    """Canonicalize a single contig to its UCSC chr-prefixed name.

    Prepend ``chr`` to a non-prefixed contig, with the mitochondrion ``MT`` ->
    ``chrM``; scaffolds get the prefix too (``1_KI270706v1_random`` ->
    ``chr1_KI270706v1_random``). Already-``chr`` names are returned unchanged, so
    this is idempotent.
    """
    if contig.startswith("chr"):
        return contig
    return "chrM" if contig == "MT" else f"chr{contig}"


def chr_rename_map(contigs: list[str]) -> dict[str, str]:
    """Rename map to canonicalize *any* contig list to UCSC chr-prefixed names.

    The rule that matches the GRCh38 no-alt analysis set (and ClinVar/gnomAD via
    chr-prefixing) is :func:`canonical_chrom`. Handles unplaced/unlocalized
    scaffolds too (``Un_KI270302v1`` -> ``chrUn_KI270302v1``). Already-``chr``
    names are left as-is. Only contigs that change are returned.
    """
    out: dict[str, str] = {}
    for c in contigs:
        canon = canonical_chrom(c)
        if canon != c:
            out[c] = canon
    return out


def chr_to_plain_map() -> dict[str, str]:
    """Inverse of :func:`plain_to_chr_map`."""
    return {v: k for k, v in plain_to_chr_map().items()}


def write_rename_file(mapping: dict[str, str], dest: Path) -> Path:
    """Write a `bcftools annotate --rename-chrs` two-column file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("".join(f"{k}\t{v}\n" for k, v in mapping.items()))
    return dest


def detect_build(info: VcfInfo) -> str:
    """Best-effort genome build label from header/contig evidence (GRCh38 expected for DRAGEN)."""
    ref = (info.reference or "").lower()
    if "38" in ref or "hg38" in ref:
        return "GRCh38"
    if "37" in ref or "hg19" in ref:
        return "GRCh37"
    # Fall back to a chromosome length signal would require contig lengths; default to GRCh38.
    return "GRCh38?"

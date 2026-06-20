#!/usr/bin/env python3
"""Generate a tiny, self-consistent synthetic DRAGEN-like genome for testing.

Writes into <base>/reference and <base>/genome so that `LOCUS_DATA_DIR=<base>`
makes the whole pipeline runnable with no downloads and no real genetic data.

The gVCF deliberately exercises the tricky bits:
  - hom-ref blocks (pure ``<NON_REF>`` rows with INFO/END) that must be dropped
  - a real SNV carrying the appended ``<NON_REF>`` allele
  - a multiallelic SNV (must split into two records)
  - a deletion inside a homopolymer run (must left-align against the reference)

Usage: python scripts/make_fixture.py /tmp/locus-fixture
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# A fixed 60 bp "chr21" with a homopolymer A-run at 9-12 to test left-alignment.
SEQ = "ACGTACGTAAAACGTAGGGTTTTACGTACGTACGTACGTTTACGTACGTACGTACGTACGT"
CONTIG = "chr21"


def ref_at(pos: int, length: int = 1) -> str:
    """1-based reference allele of given length."""
    return SEQ[pos - 1 : pos - 1 + length]


def build_vcf() -> str:
    L = len(SEQ)
    header = f"""##fileformat=VCFv4.2
##reference=file:///synthetic/GRCh38_no_alt.fa
##contig=<ID={CONTIG},length={L}>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowDepth,Description="DP<=1">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"""

    def row(pos, id_, ref, alt, info, gt, qual="50"):
        return f"{CONTIG}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\tPASS\t{info}\tGT:DP:GQ\t{gt}:30:50"

    rows = [
        # hom-ref block 1..4 -> dropped
        row(1, ".", ref_at(1), "<NON_REF>", "END=4", "0/0"),
        # real SNV with appended <NON_REF>, carries an rsID
        row(5, "rs5000", ref_at(5), "G,<NON_REF>", "DP=30", "0/1"),
        # multiallelic SNV -> splits into two records (G>T, G>C)
        row(7, ".", ref_at(7), "T,C,<NON_REF>", "DP=28", "1/2"),
        # deletion inside homopolymer A-run (9-12) -> must left-align
        row(11, ".", ref_at(11, 2), "A,<NON_REF>", "DP=25", "0/1"),
        # hom-ref block 30..39 -> dropped
        row(30, ".", ref_at(30), "<NON_REF>", "END=39", "0/0"),
    ]
    return header + "\n" + "\n".join(rows) + "\n"


# CNV/SV are emitted with a NON-chr-prefixed contig ("21") on purpose: that is how
# sequencing.com/DRAGEN ships these files, and the loader must canonicalize them to
# "chr21" so region queries (which are always chr-prefixed) can match.
def build_cnv_vcf() -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=21,length=60>\n"
        '##ALT=<ID=DEL,Description="Deletion">\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "21\t10\t.\tN\t<DEL>\t50\tPASS\tEND=30;SVTYPE=DEL;SVLEN=-20\tGT:CN\t0/1:1\n"
    )


def build_sv_vcf() -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=21,length=60>\n"
        '##ALT=<ID=INS,Description="Insertion">\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
        '##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate id">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Paired-read support">\n'
        '##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split-read support">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "21\t15\t.\tN\t<INS>\t50\tPASS\tSVTYPE=INS;SVLEN=58\tGT:PR:SR\t0/1:10,5:3,2\n"
    )


def _write_indexed_vcf(text: str, path: Path) -> Path:
    path.write_text(text)
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    gz = Path(str(path) + ".gz")
    subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
    return gz


def main(base: Path) -> None:
    ref_dir = base / "reference"
    genome_dir = base / "genome"
    ref_dir.mkdir(parents=True, exist_ok=True)
    genome_dir.mkdir(parents=True, exist_ok=True)

    # Reference FASTA -> bgzip -> faidx
    fa = ref_dir / "synthetic_GRCh38.fa"
    fa.write_text(f">{CONTIG}\n{SEQ}\n")
    subprocess.run(["bgzip", "-f", str(fa)], check=True)
    fa_gz = ref_dir / "synthetic_GRCh38.fa.gz"
    subprocess.run(["samtools", "faidx", str(fa_gz)], check=True)

    # gVCF -> bgzip -> tabix
    vcf = genome_dir / "Synthetic-SQF-30x-WGS.snp-indel.genome.vcf"
    vcf.write_text(build_vcf())
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    vcf_gz = genome_dir / "Synthetic-SQF-30x-WGS.snp-indel.genome.vcf.gz"
    subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)

    # CNV + SV (non-chr-prefixed contigs, like the real DRAGEN outputs).
    cnv_gz = _write_indexed_vcf(build_cnv_vcf(), genome_dir / "Synthetic-SQF-30x-WGS.cnv.vcf")
    sv_gz = _write_indexed_vcf(build_sv_vcf(), genome_dir / "Synthetic-SQF-30x-WGS.sv.vcf")

    print(f"reference : {fa_gz}")
    print(f"gVCF      : {vcf_gz}")
    print(f"CNV       : {cnv_gz}")
    print(f"SV        : {sv_gz}")
    print(f"\nRun:  LOCUS_DATA_DIR={base} uv run locus ingest")


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/locus-fixture").resolve()
    main(out)

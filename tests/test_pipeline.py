"""End-to-end pipeline test on a synthetic genome (no downloads, no real data).

Builds the fixture, runs ingest → ClinVar annotate → load, and asserts the
tricky transforms: NON_REF blocks dropped, multiallelics split, the homopolymer
indel left-aligned, and ClinVar clinical filtering.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _have(tool: str) -> bool:
    from shutil import which

    return which(tool) is not None


pytestmark = pytest.mark.skipif(
    not (_have("bcftools") and _have("samtools") and _have("bgzip")),
    reason="requires bcftools/samtools/htslib on PATH",
)


@pytest.fixture()
def genome(tmp_path: Path):
    """Build the synthetic fixture under tmp and point Locus at it.

    Mutates the shared ``settings`` singleton in place (rather than env+reload) so
    every module that did ``from .config import settings`` sees the same paths.
    """
    base = tmp_path / "data"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_fixture.py"), str(base)], check=True)
    from locus.config import settings

    settings.data_dir = base
    settings.db_path = base / "locus.duckdb"
    return base


def _synthetic_clinvar(base: Path) -> None:
    ann = base / "annotations"
    ann.mkdir(parents=True, exist_ok=True)
    vcf = ann / "clinvar.chr.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="">\n'
        '##INFO=<ID=CLNDN,Number=.,Type=String,Description="">\n'
        '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr21\t5\t12345\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNDN=Test_disease;GENEINFO=TESTGENE:1\n"
        "chr21\t7\t12346\tG\tT\t.\t.\tCLNSIG=Benign;GENEINFO=TESTGENE:1\n"
    )
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    subprocess.run(["bcftools", "index", "-f", "-t", str(ann / "clinvar.chr.vcf.gz")], check=True)


def test_ingest_load_query(genome, monkeypatch):
    from locus import annotate, ingest, load, queries
    from locus.config import settings

    # Ingest: gVCF -> normalized sites VCF.
    ingest.run(settings.genome_dir, normalize=True)

    # ClinVar annotate (synthetic DB).
    _synthetic_clinvar(genome)
    annotate.run(steps="clinvar")

    # Load into DuckDB.
    load.run()

    # 4 variant records: SNV(5), split multiallelic(7 x2), left-aligned indel(8).
    ov = queries.overview()
    assert ov["variants"] == 4
    assert ov["clinvar_annotated"] == 2

    # Indel was left-aligned from pos 11 (AA>A) to pos 8 (TA>T).
    page = queries.lookup_by_region("chr21:8")
    assert any(v.ref == "TA" and v.alt == "T" for v in page.hits)

    # rsID preserved.
    assert queries.lookup_by_rsid("rs5000").total == 1

    # Clinical findings default to pathogenic only (Benign excluded).
    cf = queries.clinical_findings()
    assert cf.total == 1
    assert cf.hits[0].clnsig == "Pathogenic"
    assert cf.hits[0].gene == "TESTGENE"

    # Gene lookup picks up both TESTGENE rows.
    assert queries.lookup_by_gene("TESTGENE").total == 2


def test_sql_guard_blocks_mutations(genome):
    from locus import annotate, ingest, load, queries
    from locus.config import settings

    ingest.run(settings.genome_dir, normalize=True)
    _synthetic_clinvar(genome)
    annotate.run(steps="clinvar")
    load.run()

    with pytest.raises(ValueError):
        queries.run_sql("DROP TABLE variants")
    with pytest.raises(ValueError):
        queries.run_sql("UPDATE variants SET pos = 0")
    # A legitimate read works.
    assert queries.run_sql("SELECT count(*) FROM variants")["rows"][0][0] == 4


def test_noncprefixed_bcftools_gvcf(tmp_path):
    """A sequencing.com-style gVCF: non-chr contigs + ALT='.' hom-ref blocks.

    Ingest must canonicalize contigs to chr-prefixed and drop the hom-ref blocks.
    """
    from locus import artifacts, ingest
    from locus.config import settings

    base = tmp_path / "data"
    genome = base / "genome"
    genome.mkdir(parents=True)
    vcf = genome / "Sample.snp-indel.genome.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh38.p13\n"
        "##contig=<ID=21,length=300>\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "21\t100\t.\tA\t.\t50\tPASS\tEND=110\tGT\t0/0\n"   # hom-ref block -> dropped
        "21\t120\trs1\tC\tG\t50\tPASS\t.\tGT\t0/1\n"        # real SNV -> chr21
        "21\t130\t.\tG\tT,A\t50\tPASS\t.\tGT\t1/2\n"        # multiallelic -> splits to 2
    )
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(genome / "Sample.snp-indel.genome.vcf.gz")], check=True)

    settings.data_dir = base
    settings.db_path = base / "locus.duckdb"

    info = ingest.run(genome, normalize=False)  # no reference needed to test rename + block drop
    assert info.chr_prefixed is False  # input detected as non-prefixed
    assert info.is_gvcf is True        # ALT='.' + END blocks detected as gVCF

    out = subprocess.run(
        ["bcftools", "view", "-H", str(artifacts.sites_vcf())], capture_output=True, text=True
    ).stdout.splitlines()
    assert len(out) == 3  # SNV + two split multiallelic alleles
    assert all(r.startswith("chr21\t") for r in out)  # contigs canonicalized
    assert not any(r.split("\t")[4] == "." for r in out)  # no hom-ref blocks left


def test_region_parsing():
    from locus.queries import parse_region

    assert parse_region("chr7:117,480,000-117,670,000") == ("chr7", 117480000, 117670000)
    assert parse_region("7:55050000") == ("chr7", 55050000, 55050000)
    with pytest.raises(ValueError):
        parse_region("not-a-region")

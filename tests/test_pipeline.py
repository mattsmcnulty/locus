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


def ingest_and_clinvar(base: Path) -> None:
    """Ingest the fixture and stage a synthetic ClinVar so annotate has something to join."""
    from locus import ingest
    from locus.config import settings

    ingest.run(settings.genome_dir, normalize=True)
    _synthetic_clinvar(base)


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


def test_structural_overlap_chr_canonicalized(genome):
    """CNV/SV ship non-chr-prefixed ('21'); the loader must canonicalize to 'chr21'
    so region queries (always chr-prefixed) can match. Regression: this silently
    returned 0 for every region before the fix."""
    from locus import ingest, load, queries
    from locus.config import settings

    ingest.run(settings.genome_dir, normalize=True)
    load.run()

    ov = queries.overview()
    assert ov["cnv"] >= 1 and ov["sv"] >= 1

    hits = queries.structural_overlap("chr21:1-60")
    assert len(hits) >= 2, "structural_overlap must find the fixture CNV + SV"
    assert all(h.chrom == "chr21" for h in hits)  # contigs canonicalized at load
    assert {"cnv", "sv"} <= {h.kind for h in hits}


def test_mcp_stdio_tools_respond(genome):
    """End-to-end over the real stdio MCP path a client (Claude) uses: list tools,
    then call the two handlers that were reported timing out. Guards against both
    the structural chr-prefix bug and any future serialization/hang regression."""
    import asyncio
    import os

    from locus import ingest, load
    from locus.config import settings

    ingest.run(settings.genome_dir, normalize=True)
    load.run()

    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    env = {**os.environ, "LOCUS_DATA_DIR": str(genome)}
    params = StdioServerParameters(command=sys.executable, args=["-m", "locus.mcp_server"], env=env)

    async def _run():
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await asyncio.wait_for(session.initialize(), timeout=30)
                tools = await asyncio.wait_for(session.list_tools(), timeout=30)
                by_name = {t.name: t for t in tools.tools}
                assert {"structural_variants", "polygenic_risk", "genome_overview"} <= set(by_name)

                # These two must expose a named object output schema, NOT FastMCP's
                # generic {"result": [...]} list wrapper, which strict clients won't
                # dispatch a call against (the original hang).
                assert set(by_name["structural_variants"].outputSchema["properties"]) == {"total", "hits"}
                assert list(by_name["polygenic_risk"].outputSchema["properties"]) == ["scores"]

                # The handler that was effectively broken: must return >0 hits now.
                sv = await asyncio.wait_for(
                    session.call_tool("structural_variants", {"region": "chr21:1-60"}), timeout=30
                )
                assert not sv.isError
                assert sv.structuredContent["total"] >= 2
                assert len(sv.structuredContent["hits"]) >= 2

                # Must respond promptly without hanging (empty is fine on the fixture,
                # which has no ancestry/PGS step).
                pr = await asyncio.wait_for(session.call_tool("polygenic_risk", {}), timeout=30)
                assert not pr.isError
                assert isinstance(pr.structuredContent["scores"], list)

    asyncio.run(_run())


def test_region_parsing():
    from locus.queries import parse_region

    assert parse_region("chr7:117,480,000-117,670,000") == ("chr7", 117480000, 117670000)
    assert parse_region("7:55050000") == ("chr7", 55050000, 55050000)
    with pytest.raises(ValueError):
        parse_region("not-a-region")


def test_classify_clinvar_delta():
    """ClinVar reanalysis diff: newly-pathogenic (tiered by review status), reclassified,
    de-pathogenized, and withdrawn. Pure function — no I/O."""
    from locus.refresh import classify_clinvar_delta

    def v(clnsig, rev="criteria_provided,_multiple_submitters,_no_conflicts", gene="G"):
        return {"clnsig": clnsig, "clnrevstat": rev, "gene": gene, "rsid": "rs1", "clndn": "Some_disease"}

    prev = {
        ("chr1", 100, "A", "G"): v("Benign"),                  # -> becomes pathogenic
        ("chr1", 200, "C", "T"): v("Pathogenic"),              # -> loses pathogenic (depathogenized)
        ("chr1", 250, "C", "A"): v("Benign"),                  # -> Benign->VUS (true reclassified)
        ("chr1", 300, "G", "A"): v("Likely_pathogenic"),       # -> de-pathogenized
        ("chr1", 400, "T", "C"): v("Pathogenic"),              # -> withdrawn (gone in cur)
        ("chr1", 500, "A", "T"): v("Benign"),                  # -> unchanged, no finding
    }
    cur = {
        ("chr1", 100, "A", "G"): v("Pathogenic", rev="criteria_provided,_single_submitter"),  # 1-star
        ("chr1", 200, "C", "T"): v("Uncertain_significance"),
        ("chr1", 250, "C", "A"): v("Uncertain_significance"),
        ("chr1", 300, "G", "A"): v("Benign"),
        ("chr1", 500, "A", "T"): v("Benign"),
    }
    out = {(f.kind, f.pos): f for f in classify_clinvar_delta(prev, cur)}
    assert out[("newly_pathogenic", 100)].tier == "moderate"   # single-submitter -> 1 star
    assert out[("depathogenized", 200)].tier == "moderate"     # lost pathogenic status
    assert out[("reclassified", 250)].tier == "weak"           # neither side pathogenic
    assert out[("depathogenized", 300)].tier == "moderate"
    assert out[("withdrawn", 400)].tier == "moderate"
    assert not any(f.pos == 500 for f in classify_clinvar_delta(prev, cur))  # unchanged -> nothing

    # Multi-submitter newly-pathogenic is the strong tier.
    strong = classify_clinvar_delta({("chr2", 9, "A", "T"): v("Benign")},
                                     {("chr2", 9, "A", "T"): v("Pathogenic")})
    assert strong[0].kind == "newly_pathogenic" and strong[0].tier == "strong"


def test_whats_new_query(genome):
    """whats_new ranks strongest-first, filters by tier/since, and is empty before refresh."""
    from locus import ingest, load, queries
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()

    # No findings yet.
    assert queries.whats_new().total == 0

    rows = [
        ("2026-06-01T00:00:00", "clinvar", "newly_pathogenic", "strong", "chr1", 100,
         "A", "G", "rs1", "BRCA1", "strong one", "", "Benign", "Pathogenic", None, None),
        ("2026-06-10T00:00:00", "pubmed", "new_study", "info", None, None,
         None, None, None, "BRCA1", "info one", "abstract…", None, None, "2026",
         "https://pubmed.ncbi.nlm.nih.gov/999/"),
    ]
    with connect(read_only=False) as con:
        load.append_findings(con, rows)

    wn = queries.whats_new()
    assert wn.total == 2
    assert wn.findings[0].tier == "strong"            # strongest first
    assert wn.counts_by_tier == {"strong": 1, "info": 1}
    assert queries.whats_new(tier="strong").total == 1
    assert queries.whats_new(since="2026-06-05").total == 1  # only the info one
    # v4 `url` column round-trips (clickable citation on PubMed/GWAS findings).
    info = queries.whats_new(tier="info").findings[0]
    assert info.url == "https://pubmed.ncbi.nlm.nih.gov/999/"


def test_panels_integrity():
    """Tag-SNP data sanity: each has 0/1/2 interpretations, chr-prefixed coords, distinct
    single-base effect/other alleles; ACMG gene set is non-empty + uppercase."""
    from locus.panels import ACMG_SF_GENES, TAG_SNPS

    assert len(ACMG_SF_GENES) >= 70
    assert all(g == g.upper() and g.isalnum() for g in ACMG_SF_GENES)
    assert {"BRCA1", "BRCA2", "LDLR", "KCNQ1"} <= ACMG_SF_GENES
    seen = set()
    for s in TAG_SNPS:
        assert set(s.interp) == {0, 1, 2}, f"{s.rsid} needs 0/1/2 interpretations"
        assert s.chrom.startswith("chr")
        assert s.effect_allele in "ACGT" and s.other_allele in "ACGT"
        assert s.effect_allele != s.other_allele
        assert s.rsid not in seen
        seen.add(s.rsid)
    # The clinically-important HLA-B*57:01 proxy must be present and flagged pharmacogenomic.
    hla = next(s for s in TAG_SNPS if s.rsid == "rs2395029")
    assert hla.category == "pharmacogenomic"


def test_secondary_findings_acmg_filter(genome):
    """ACMG SF returns pathogenic variants only in actionable genes."""
    from locus import ingest, load, queries
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    with connect(read_only=False) as con:
        con.executemany(
            "INSERT INTO variants (chrom, pos, ref, alt, gene, clnsig) VALUES (?,?,?,?,?,?)",
            [("chr17", 43000000, "A", "G", "BRCA1", "Pathogenic"),       # ACMG gene -> included
             ("chr1", 12345, "C", "T", "MADEUPGENE", "Pathogenic")],     # not ACMG -> excluded
        )
    sf = queries.secondary_findings()
    genes = {h.gene for h in sf.hits}
    assert "BRCA1" in genes
    assert "MADEUPGENE" not in genes


def test_traits_compute_runs(genome):
    """traits.compute() runs end-to-end via markers_genotypes; on the chr21-only fixture every
    tag SNP is off-contig, so all come back not-callable (exercises the wiring without crashing)."""
    from locus import ingest, traits
    from locus.config import settings
    from locus.panels import TAG_SNPS

    ingest.run(settings.genome_dir, normalize=True)
    results = traits.compute()
    assert len(results) == len(TAG_SNPS)
    assert all(r.dosage is None and r.genotype == "—" for r in results)  # off-contig fixture


def test_gwas_parse_filters(tmp_path):
    """GWAS parse keeps only genome-wide-significant single lead SNPs with a real risk allele."""
    from locus.gwas import parse

    def row(**kw):
        f = [""] * 36
        f[1], f[7] = kw.get("pmid", "1"), kw.get("trait", "T")
        f[11], f[12] = kw.get("chr", "4"), kw.get("pos", "100")
        f[20], f[21] = kw.get("strongest", "rs1-A"), kw.get("snps", "rs1")
        f[27], f[30], f[34] = kw.get("p", "1e-10"), kw.get("orb", "1.2"), kw.get("mapped", "trait X")
        return "\t".join(f)

    tsv = tmp_path / "g.tsv"
    tsv.write_text("\n".join([
        "header",
        row(snps="rs1", strongest="rs1-A", chr="4", pos="100", p="1e-10"),    # keep
        row(snps="rs2", strongest="rs2-T", chr="7", pos="200", p="1e-3"),      # drop: not significant
        row(snps="rs3", strongest="rs3-?", chr="7", pos="300", p="1e-20"),     # drop: no risk allele
        row(snps="rs4; rs5", strongest="rs4-A", chr="7", pos="400", p="1e-20"),  # drop: multi-SNP
        row(snps="rs6", strongest="rs6-G", chr="7 x 8", pos="500", p="1e-20"),   # drop: multi-locus
    ]) + "\n")
    out = parse(tsv)
    assert {a.rsid for a in out} == {"rs1"}
    a = out[0]
    assert a.chrom == "chr4" and a.pos == 100 and a.risk_allele == "A"


def test_risk_dosage():
    """Risk-allele dosage: hom-ref-at-non-risk is 0 (not excluded); strand flip handled;
    risk allele absent from site -> None."""
    from locus.gwas import _risk_dosage

    g = lambda a, b: [a, b, False]  # noqa: E731 - cyvcf2-style genotype
    assert _risk_dosage("A", ["G"], g(0, 0), "G") == 0    # hom-ref, risk=alt not present
    assert _risk_dosage("A", ["G"], g(0, 1), "G") == 1    # het
    assert _risk_dosage("A", ["G"], g(1, 1), "G") == 2    # hom-alt
    assert _risk_dosage("A", ["G"], g(0, 0), "A") == 2    # hom-ref, risk == ref
    assert _risk_dosage("A", ["G"], g(0, 1), "C") == 1    # strand flip C->G
    assert _risk_dosage("A", ["T"], g(0, 1), "C") is None  # risk neither at site nor its complement
    assert _risk_dosage("A", ["G"], g(-1, -1), "G") is None  # no-call


def test_refresh_dry_run_no_writes(genome, monkeypatch):
    """Orchestrator (monkeypatched checkers, no network): dry-run surfaces findings but
    writes neither the manifest nor the DB."""
    from locus import manifest, refresh
    from locus.config import settings

    monkeypatch.setattr(refresh, "check_clinvar",
                        lambda m: {"version": "newmd5", "checksum": "newmd5", "url": "u", "changed": True})
    monkeypatch.setattr(refresh, "check_pgs",
                        lambda m: {"version": "2099-01-01", "url": "u", "n_new": 3, "ids": ["PGS1"], "changed": True})
    # Keep the test hermetic: stub the network-touching checkers (CPIC/GWAS probe live APIs).
    monkeypatch.setattr(refresh, "check_cpic", lambda m: None)
    monkeypatch.setattr(refresh, "check_gwas",
                        lambda m: {"version": "2099-01-01", "url": "u", "changed": True})

    findings = refresh.run(dry_run=True)
    kinds = {f.kind for f in findings}
    assert "release" in kinds                 # PGS summary
    assert "update_available" in kinds        # ClinVar + GWAS updates flagged (work NOT run in dry-run)
    assert not manifest.manifest_path().exists()   # nothing written
    assert not (settings.data_dir / "reports" / "whats_new.md").exists()


def test_classify_gwas_delta():
    """GWAS re-analysis diff surfaces only newly-carried (rsID, trait) associations, as weak hits
    with a PubMed citation. Pure function — no I/O."""
    from locus.refresh import classify_gwas_delta

    def a(dosage, pmid, orb="1.2"):
        return {"dosage": dosage, "pmid": pmid, "chrom": "chr7", "pos": 9, "or_beta": orb, "pval": 1e-12}

    prev = {("rs1", "height"): a(1, "100")}
    cur = {
        ("rs1", "height"): a(1, "100"),                 # unchanged -> nothing
        ("rs2", "type 2 diabetes"): a(2, "200"),        # newly carried -> one finding
    }
    out = classify_gwas_delta(prev, cur)
    assert len(out) == 1
    f = out[0]
    assert f.rsid == "rs2" and f.kind == "new_association" and f.tier == "weak"
    assert f.url == "https://pubmed.ncbi.nlm.nih.gov/200/"
    assert "2 risk allele" in f.title
    assert classify_gwas_delta(cur, cur) == []          # nothing new vs itself


def test_is_carried_from_genotype():
    """Carrier detection needs the reference allele to tell hom-alt from hom-ref."""
    from locus.literature import _is_carried

    assert _is_carried("A/G", "A") is True     # het
    assert _is_carried("G/G", "A") is True     # hom-alt
    assert _is_carried("A/A", "A") is False    # hom-ref
    assert _is_carried("—", "A") is False      # not callable
    assert _is_carried("A/G", None) is True    # het still detectable without ref
    assert _is_carried("G/G", None) is False   # hom-alt indistinguishable from hom-ref w/o ref


def test_literature_term_for():
    from locus.literature import _term_for

    assert _term_for("rs334") == "rs334"                      # rsID stays literal
    assert _term_for("BRCA2").startswith("BRCA2[Gene]")       # bare symbol gets [Gene]
    assert _term_for("breast cancer risk") == "breast cancer risk"  # free text passes through


def test_pubmed_search_parses(monkeypatch):
    """pubmed_search: esearch → PMIDs, efetch XML → hydrated PubMedHit (title/abstract/year)."""
    from locus import literature

    xml = (
        '<?xml version="1.0"?><PubmedArticleSet><PubmedArticle><MedlineCitation>'
        "<PMID>111</PMID><Article>"
        "<ArticleTitle>BRCA1 and cancer risk</ArticleTitle>"
        "<Abstract><AbstractText>We found a variant.</AbstractText></Abstract>"
        "<Journal><Title>J Test</Title><JournalIssue><PubDate><Year>2026</Year>"
        "</PubDate></JournalIssue></Journal></Article></MedlineCitation></PubmedArticle>"
        "</PubmedArticleSet>"
    )

    def fake_eutils(endpoint, params, *, want_json):
        return {"esearchresult": {"idlist": ["111"]}} if endpoint == "esearch.fcgi" else xml

    monkeypatch.setattr(literature, "_eutils", fake_eutils)
    hits = literature.pubmed_search("BRCA1[Gene]", gene="BRCA1")
    assert len(hits) == 1
    h = hits[0]
    assert h.pmid == "111" and "cancer risk" in h.title
    assert h.abstract == "We found a variant." and h.year == "2026" and h.gene == "BRCA1"
    assert h.url == "https://pubmed.ncbi.nlm.nih.gov/111/"


def test_study_rsids(monkeypatch):
    """A PubMed ID → GWAS Catalog study → its reported rsIDs (regex-harvested, resilient)."""
    from locus import literature

    def fake_get(path_or_url):
        if "findByPublicationIdPubmedId" in path_or_url:
            return {"_embedded": {"studies": [
                {"accessionId": "GCST1", "_links": {"associations": {"href": "http://x/assoc"}}}]}}
        return {"_embedded": {"associations": [
            {"loci": [{"strongestRiskAlleles": [{"riskAlleleName": "rs7903146-T"}]}]}]}}

    monkeypatch.setattr(literature, "_gwas_get", fake_get)
    assert literature.study_rsids("12345") == ["rs7903146"]


def test_pubmed_findings_dedup(genome, monkeypatch):
    """_pubmed_findings skips already-seen PMIDs, and on first run seeds the baseline silently."""
    from locus import ingest, literature, load, refresh
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    with connect(read_only=False) as con:  # a notable gene + one already-seen PMID
        con.execute("INSERT INTO variants (chrom, pos, ref, alt, gene, clnsig) "
                    "VALUES ('chr17', 43000000, 'A', 'G', 'BRCA1', 'Pathogenic')")
        con.execute("CREATE TABLE IF NOT EXISTS watch_seen_ids(source VARCHAR, external_id VARCHAR)")
        con.execute("INSERT INTO watch_seen_ids VALUES ('pubmed', 'SEEN1')")

    hits = [literature.PubMedHit(pmid="SEEN1", title="old", gene="BRCA1"),
            literature.PubMedHit(pmid="NEW1", title="new paper", abstract="A", gene="BRCA1")]
    monkeypatch.setattr(literature, "pubmed_search", lambda *a, **k: hits)

    check = {"since": "2026-01-01", "first_run": False}
    findings = refresh._pubmed_findings(check)
    assert [f.title for f in findings] == ["New paper on BRCA1: new paper"]
    assert check["new_pmids"] == ["NEW1"]          # SEEN1 excluded, not re-marked

    check2 = {"since": "2026-01-01", "first_run": True}
    assert refresh._pubmed_findings(check2) == []  # baseline: emit nothing…
    assert check2["new_pmids"] == ["NEW1"]         # …but still record the unseen PMID


def test_annotate_refuses_to_shrink_the_store(genome, monkeypatch):
    """annotate always rebuilds from the sites VCF, so a step that is skipped (or simply not
    requested) is DROPPED from the result. That is how gnomAD AF and SnpEff consequences each
    silently vanished from a working store. Re-running a subset must refuse, not quietly shrink."""
    from locus import annotate

    ingest_and_clinvar(genome)
    annotate.run(steps="clinvar")
    # applied_steps reads the VCF's own INFO headers — it cannot drift from what's really there.
    assert annotate.applied_steps() == {"clinvar"}

    # Pretend the store also carries snpeff+gnomad; re-running clinvar-only would drop them.
    monkeypatch.setattr(annotate, "applied_steps",
                        lambda vcf=None: {"clinvar", "snpeff", "gnomad"})
    with pytest.raises(Exception) as e:
        annotate.run(steps="clinvar")
    msg = str(e.value)
    assert "gnomad" in msg and "snpeff" in msg, "the error must name what would be lost"
    assert "--force" in msg, "and how to override it deliberately"

    # --force is the explicit escape hatch.
    annotate.run(steps="clinvar", force=True)


def test_annotate_hard_fails_on_zero_join(genome):
    """A step that runs without error but matches NOTHING is the failure this codebase keeps
    hitting: the annotation is absent, every tool answers null, and the pipeline reports success.
    The count was already being computed and thrown away; now it must raise."""
    from locus import annotate, artifacts, ingest, shell
    from locus.config import settings

    ingest.run(settings.genome_dir, normalize=True)
    sites = artifacts.sites_vcf()

    # A tag no record carries is exactly what a contig mismatch produces: the command exits 0,
    # writes a valid VCF, and joins nothing.
    with pytest.raises(shell.ToolError) as e:
        annotate._require_join(sites, "NO_SUCH_TAG", "ClinVar", "fix hint here")
    assert "annotated 0 records" in str(e.value)
    assert "contig" in str(e.value).lower(), "the error must name the usual cause"
    assert "fix hint here" in str(e.value)

    # A tag that IS present returns the count instead of raising.
    assert annotate._require_join(sites, "chr21", "Test", "") == 4


def test_download_all_raises_on_failure(monkeypatch):
    """`download all` kept going after a failure AND returned normally, so setup marched on and
    printed 'Locus is ready' over a half-built genome. Resilient must not mean silent."""
    from locus import download, shell

    monkeypatch.setattr(download, "TARGETS", {
        "good": lambda: None,
        "bad": lambda: (_ for _ in ()).throw(RuntimeError("network down")),
    })
    monkeypatch.setattr(download, "guidance_large", lambda: None)
    with pytest.raises(shell.ToolError) as e:
        download.run("all")
    assert "bad" in str(e.value)          # names what failed
    assert "1 of 2 downloads failed" in str(e.value)


def test_reanalyze_clinvar_aborts_on_annotation_collapse(monkeypatch):
    """If a re-annotate loses ClinVar, `cur` is empty and diffing would emit a fabricated
    'pathogenic record withdrawn — re-check' for every P/LP variant carried."""
    from locus import refresh, shell

    prev = {("chr1", i, "A", "G"): {"clnsig": "Pathogenic", "clnrevstat": "", "gene": "G",
                                    "rsid": "rs1", "clndn": "d"} for i in range(100)}
    snaps = iter([prev, {}])   # before: 100 annotated; after: collapsed to 0
    monkeypatch.setattr(refresh, "_snapshot_clinvar", lambda: next(snaps))
    monkeypatch.setattr(refresh.download, "download_clinvar", lambda: None)
    monkeypatch.setattr(refresh.annotate, "run", lambda **k: None)
    monkeypatch.setattr(refresh.load, "run", lambda: None)
    with pytest.raises(shell.ToolError) as e:
        refresh._reanalyze_clinvar()
    assert "collapsed" in str(e.value)

    # Sanity: the pure differ WOULD have produced the fiction we're preventing.
    assert len(refresh.classify_clinvar_delta(prev, {})) == 100


def test_study_variants_distinguishes_failure_from_zero(monkeypatch):
    """'We couldn't look it up' and 'you carry none' are opposite answers."""
    from locus import gwas, literature

    monkeypatch.setattr(literature, "study_rsids", lambda pmid: ["rs1", "rs2"])
    monkeypatch.setattr(gwas, "ask_markers", lambda rsids: [])   # Ensembl down
    res = literature.study_variants("123")
    assert res["carried"] == 0
    assert "NOT" in res["note"] and "unknown" in res["note"].lower()

    monkeypatch.setattr(gwas, "ask_markers",
                        lambda rsids: [{"rsid": "rs1", "genotype": "A/G", "ref": "A"}])
    res = literature.study_variants("123")
    assert res["carried"] == 1
    assert "could not be resolved" in res["note"]   # 1 of 2 unresolved is surfaced, not hidden


def test_carrier_status_zygosity_and_honesty(genome):
    """Carrier = one copy (silent for you, matters for your children). Two copies = possibly
    affected. And the report must ALWAYS name what it cannot assess: SMN1/FMR1 are among the most
    important carrier screens and a VCF cannot answer them, so an empty result that implied
    otherwise would be the most dangerous output this tool could produce."""
    from locus import ingest, load, queries
    from locus.config import settings
    from locus.db import connect
    from locus.panels import CARRIER_PANEL

    ingest.run(settings.genome_dir, normalize=True)
    load.run()

    # Even with zero hits, the "cannot assess" list must be present and non-empty.
    empty = queries.carrier_status()
    assert empty.total == 0
    assert {n.gene for n in empty.not_assessed} >= {"SMN1", "FMR1"}
    assert "not a negative screen" in empty.note.lower() or "NOT a clinical carrier screen" in empty.note

    with connect(read_only=False) as con:
        con.executemany(
            "INSERT INTO variants (chrom,pos,ref,alt,gene,clnsig,gt,rsid) VALUES (?,?,?,?,?,?,?,?)",
            [("chr12", 102917130, "C", "T", "PAH", "Pathogenic", "0/1", "rsA"),      # 1 copy -> carrier
             ("chr7", 117559590, "A", "G", "CFTR", "Pathogenic", "1/1", "rsB"),      # hom -> affected
             ("chr15", 72346580, "G", "A", "HEXA", "Pathogenic", "0/1", "rsC"),      # 2 P/LP in one
             ("chr15", 72346590, "T", "C", "HEXA", "Likely_pathogenic", "0/1", "rsD"),  # gene -> affected
             ("chr1", 100, "A", "G", "PAH", "Benign", "0/1", "rsE")],                # benign -> ignored
        )
    by_gene = {h.gene: h for h in queries.carrier_status().hits}
    assert by_gene["PAH"].status == "carrier" and by_gene["PAH"].n_variants == 1
    assert by_gene["CFTR"].status == "likely_affected", "homozygous P/LP is not merely a carrier"
    assert by_gene["HEXA"].status == "likely_affected", "two P/LP in one gene = presumed compound het"
    assert by_gene["PAH"].condition == "Phenylketonuria (PKU)"
    # Benign never contributes.
    assert by_gene["PAH"].n_variants == 1

    # likely_affected sorts first — it's the part that isn't just about family planning.
    assert queries.carrier_status().hits[0].status == "likely_affected"
    assert len(CARRIER_PANEL) == queries.carrier_status().panel_size


def test_carrier_panel_integrity():
    from locus.panels import CARRIER_PANEL, CARRIER_UNASSESSABLE

    genes = [c.gene for c in CARRIER_PANEL]
    assert len(genes) == len(set(genes)), "no duplicate genes"
    assert all(c.inheritance in ("AR", "XL") for c in CARRIER_PANEL)
    assert all(c.condition for c in CARRIER_PANEL)
    # The unassessable genes must NOT be in the panel — listing them as screened is the lie.
    unassessable = {g for g, _, _ in CARRIER_UNASSESSABLE}
    assert not (unassessable & set(genes)), "a gene we cannot assess must not appear as screened"
    assert all(why for _, _, why in CARRIER_UNASSESSABLE), "each must explain why"


def test_acmg_panel_matches_v33():
    """The ACMG SF panel is the whole basis of `secondary_findings`, whose empty result is
    reported to the user as reassuring. A gene missing here is a screen that silently never
    happens. Pinned to ACMG SF v3.3 (84 genes; Lee K, et al. Genet Med. 2025;27(8):101454)."""
    from locus.panels import ACMG_SF_GENES, ACMG_SF_RECESSIVE, ACMG_SF_VERSION

    assert ACMG_SF_VERSION == "v3.3"
    assert len(ACMG_SF_GENES) == 84, "v3.3 is exactly 84 genes"
    # Regression: these four (calmodulinopathy / CPVT — sudden cardiac death) were absent,
    # so a pathogenic variant in them returned "no findings" and Claude reassured the user.
    for g in ("CALM1", "CALM2", "CALM3", "TRDN"):
        assert g in ACMG_SF_GENES, f"{g} is on ACMG SF v3.3"
    for g in ("ABCD1", "CYP27A1", "PLN"):   # added in v3.3
        assert g in ACMG_SF_GENES
    # CDH1 is not an ACMG SF gene at any version — carrying it made the panel not-ACMG-SF.
    assert "CDH1" not in ACMG_SF_GENES
    assert ACMG_SF_RECESSIVE <= ACMG_SF_GENES, "recessive genes must be part of the panel"
    assert {"MUTYH", "HFE", "ATP7B"} <= ACMG_SF_RECESSIVE


def test_predicted_damaging_excludes_clinvar_classified(genome):
    """The tool's whole premise is the 'ClinVar is silent' set. Without a clnsig IS NULL filter it
    returned ClinVar-BENIGN variants as damaging findings ClinVar had missed — a false alarm on
    health data, stated confidently to Claude."""
    from locus import ingest, load, queries
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    with connect(read_only=False) as con:
        con.executemany(
            "INSERT INTO variants (chrom,pos,ref,alt,gene,clnsig,am_class) VALUES (?,?,?,?,?,?,?)",
            [("chr1", 900001, "A", "G", "GENEA", None, "likely_pathogenic"),       # silent -> include
             ("chr1", 900002, "A", "G", "GENEB", "Benign", "likely_pathogenic"),   # benign  -> EXCLUDE
             ("chr1", 900003, "A", "G", "GENEC", "Pathogenic", "likely_pathogenic")],  # known -> exclude
        )
    hits = queries.predicted_damaging(limit=50).hits
    genes = {h.gene for h in hits}
    assert "GENEA" in genes
    assert "GENEB" not in genes, "ClinVar-benign must never be reported as predicted-damaging"
    assert "GENEC" not in genes, "ClinVar already classified it — not the 'silent' set"
    assert all(h.clnsig is None for h in hits)


def test_secondary_findings_recessive_gating(genome):
    """ACMG reports its recessive genes only when two P/LP variants are present. A single het
    carrier (HFE p.C282Y is ~10% of Europeans) is not a secondary finding."""
    from locus import ingest, load, queries
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    with connect(read_only=False) as con:
        con.executemany(
            "INSERT INTO variants (chrom,pos,ref,alt,gene,clnsig,gt) VALUES (?,?,?,?,?,?,?)",
            [("chr6", 26090951, "G", "A", "HFE", "Pathogenic", "0/1"),      # lone het carrier -> excluded
             ("chr11", 108100000, "C", "T", "ATM", "Pathogenic", "0/1"),    # not on ACMG SF -> excluded
             ("chr17", 43000000, "A", "G", "BRCA1", "Pathogenic", "0/1")],  # dominant het -> reported
        )
    genes = {h.gene for h in queries.secondary_findings().hits}
    assert "BRCA1" in genes, "dominant gene: a single P/LP het is a real finding"
    assert "HFE" not in genes, "recessive gene: a lone het carrier is not an ACMG secondary finding"
    assert "ATM" not in genes

    # Homozygous in a recessive gene IS biallelic -> report it.
    with connect(read_only=False) as con:
        con.execute("UPDATE variants SET gt='1/1' WHERE gene='HFE'")
    assert "HFE" in {h.gene for h in queries.secondary_findings().hits}


def test_write_associations_refuses_to_wipe(genome):
    """DROP was unconditional while INSERT was guarded, so an empty compute (a bcftools/network
    hiccup) silently erased a populated table — and the next refresh diff then saw prev->0 and
    reported nothing wrong."""
    from locus import gwas, ingest, load
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    carried = [gwas.Carried("rs1", "chr1", 1, "A", 1, "heterozygous", "T", "t", 1e-9, "1.2", "123")]
    load.write_associations(carried)
    with connect(read_only=True) as con:
        assert con.execute("SELECT count(*) FROM associations").fetchone()[0] == 1

    load.write_associations([])   # must NOT erase the good table
    with connect(read_only=True) as con:
        assert con.execute("SELECT count(*) FROM associations").fetchone()[0] == 1, \
            "an empty result means upstream failed, not that you carry zero associations"


def test_gnomad_scope_expression(tmp_path, monkeypatch):
    """AF is only fetched where it can change an answer (AlphaMissense-pathogenic / non-benign
    ClinVar). Scoping is what keeps this ~1k rsIDs against a rate-limited public API instead of
    5.1M. The expression must use only INFO fields actually present (a bogus tag makes bcftools
    error), and be None when there's nothing to scope by so the caller skips entirely."""
    from locus import annotate

    def with_info(*ids):
        monkeypatch.setattr(annotate, "_present_info", lambda _p: set(ids))
        return annotate._gnomad_scope(tmp_path / "x.vcf.gz")

    # Nothing upstream ran -> no scoping possible -> skip (never "fetch everything").
    assert with_info() is None
    assert with_info("DP", "MQ") is None

    # AlphaMissense only: pathogenic (incl. likely_pathogenic) — not every scored variant.
    assert with_info("am_class") == 'INFO/am_class~"pathogenic"'

    # ClinVar only: present AND not one of the benign calls.
    e = with_info("CLNSIG")
    assert 'INFO/CLNSIG!="."' in e
    for b in ("Benign", "Benign/Likely_benign", "Likely_benign"):
        assert f'INFO/CLNSIG!="{b}"' in e
    # Regression: bcftools' `!~` silently fails to filter CLNSIG (3,792 in -> 3,792 out), which
    # made the "tight" scope select 43,460 variants instead of ~1k. Exact != is required.
    assert "!~" not in e, "must not use the regex operator to exclude benign — it does not filter"

    # Absent fields must never leak into the expression.
    assert "am_class" not in with_info("CLNSIG")
    assert "CLNSIG" not in with_info("am_class")

    # Both present -> OR'd together, shell-safe.
    e = with_info("CLNSIG", "am_class")
    assert 'INFO/am_class~"pathogenic"' in e and 'INFO/CLNSIG!="."' in e and "||" in e
    assert "'" not in e, "expression is single-quoted into a shell command"


def test_ensembl_af_batching_is_pinned():
    """`pops=1` returns ~120 population records per variant, so the RESPONSE is the constraint,
    not the request count. Measured against the live API: 25 ids = ~220KB in 20-36s; 50 ids
    reliably times out. A 180-id batch with a 30s timeout failed 100% of the time and looked
    exactly like rate-limiting — it wasn't. Keep batch small and the timeout well above 36s."""
    import inspect

    from locus import annotate

    sig = inspect.signature(annotate._ensembl_gnomad_af)
    assert sig.parameters["batch"].default <= 25, "batches >25 time out against Ensembl"
    src = inspect.getsource(annotate._ensembl_gnomad_af)
    assert "timeout=90" in src, "timeout must exceed the observed 20-36s response time"


def test_gnomad_runs_last_in_chain():
    """gnomAD must run after clinvar/snpeff/alphamissense — it scopes by their annotations, so
    running it earlier would leave nothing to scope by and force a full 5.1M-position fetch."""
    import inspect

    from locus import annotate, refresh

    body = inspect.getsource(annotate.run)
    order = [body.index(f'"{s}" in requested') for s in ("clinvar", "alphamissense", "gnomad")]
    assert order == sorted(order), "gnomad must be applied last in annotate.run()"
    # gnomAD is deliberately OUT of the weekly refresh chain: streaming its AF is latency-bound
    # and would burn ~25min per run for nothing. Every other local-DB step must stay in, or a
    # refresh silently drops that annotation from the store (the original AF-disappeared bug).
    for step in ("clinvar", "snpeff", "alphamissense"):
        assert step in refresh._REANNOTATE_STEPS


def test_doctor_flags_missing_java(genome, monkeypatch, capsys):
    """SnpEff/PharmCAT/Haplogrep silently self-skip without a JDK. doctor used to report a green
    'ok' just because the jar existed — that false-green hid SnpEff never running. Guard it."""
    from locus import cli
    from locus.config import settings

    ann = settings.annotations_dir
    (ann / "snpEff").mkdir(parents=True, exist_ok=True)
    (ann / "snpEff" / "snpEff.jar").write_text("")       # jar present…
    monkeypatch.setattr(cli, "_resolve_java", lambda: None)  # …but no working Java
    cli.doctor()
    out = capsys.readouterr().out
    assert "needs java" in out, "doctor must flag a present-but-unusable Java-backed tool"

    # With Java working, the same row is a legitimate ok.
    monkeypatch.setattr(cli, "_resolve_java", lambda: "/usr/bin/java")
    cli.doctor()
    assert "needs java" not in capsys.readouterr().out


def test_litvar_pmids(monkeypatch):
    """LitVar2 rsID→PMIDs: keep the newest (highest) PMIDs, capped; reject non-rs ids."""
    from locus import literature

    monkeypatch.setattr(literature, "_litvar_get",
                        lambda path: {"pmids": [100, 300, 200, 50], "pmids_count": 4})
    assert literature.litvar_pmids("rs1800896", recent=2) == ["300", "200"]  # newest first, capped
    assert literature.litvar_pmids("notanrs") == []                          # non-rs id rejected


def test_litvar_findings_dedup(genome, monkeypatch):
    """Variant-level watcher: new papers about a specific rsID you carry, deduped; baseline silent."""
    from locus import ingest, literature, load, refresh
    from locus.config import settings
    from locus.db import connect

    ingest.run(settings.genome_dir, normalize=True)
    load.run()
    with connect(read_only=False) as con:  # a clinically-notable rsID + one already-seen PMID
        con.execute("INSERT INTO variants (chrom,pos,ref,alt,rsid,gene,clnsig) "
                    "VALUES ('chr2', 47000000, 'A','G','rs777','MSH2','Pathogenic')")
        con.execute("CREATE TABLE IF NOT EXISTS watch_seen_ids(source VARCHAR, external_id VARCHAR)")
        con.execute("INSERT INTO watch_seen_ids VALUES ('litvar','SEEN9')")

    monkeypatch.setattr(literature, "litvar_pmids", lambda rsid, **k: ["SEEN9", "NEW9"])
    monkeypatch.setattr(literature, "fetch_pubmed_meta",
                        lambda pmids: {"NEW9": literature.PubMedHit(pmid="NEW9", title="fresh finding", abstract="B")})

    check = {"first_run": False}
    f = refresh._litvar_findings(check)
    assert [x.rsid for x in f] == ["rs777"]
    assert "MSH2" in f[0].title and f[0].url.endswith("/NEW9/")
    assert check["new_pmids"] == ["NEW9"]           # SEEN9 excluded

    check2 = {"first_run": True}
    assert refresh._litvar_findings(check2) == []   # baseline: silent
    assert check2["new_pmids"] == ["NEW9"]

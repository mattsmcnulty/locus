"""Locus MCP server — lets Claude query your genome.

A stdio MCP server exposing read-only tools over the DuckDB store. Returns
typed Pydantic models so the SDK emits a precise output schema, and paginates
everything (Claude Code truncates MCP output past ~10k tokens).

Run directly (``python -m locus.mcp_server``) or via ``locus serve mcp``.
Register with Claude Code:

    claude mcp add --scope project --transport stdio locus -- uv run locus-mcp
"""

from __future__ import annotations

from mcp.server.fastmcp import FastMCP
from pydantic import BaseModel, Field

from . import queries
from .db import db_exists

mcp = FastMCP("locus")


# Tools must return an object-typed structured output. A bare ``list[...]`` makes
# FastMCP emit a generic ``{"result": [...]}`` wrapper schema that strict MCP
# clients (e.g. Claude Desktop) refuse to dispatch a call against — the request
# never even reaches the server. Wrap list returns in a named model so the output
# schema matches the working tools (like VariantPage).
class PolygenicRiskReport(BaseModel):
    scores: list[queries.PgsResult] = Field(description="One entry per curated polygenic score")


class StructuralVariantsResult(BaseModel):
    total: int = Field(description="Number of CNV/SV records overlapping the region")
    hits: list[queries.StructuralHit]


class LiteratureHit(BaseModel):
    pmid: str
    title: str
    abstract: str = ""
    journal: str = ""
    year: str = ""
    gene: str = ""
    url: str = Field(description="PubMed link")


class LiteratureResult(BaseModel):
    query: str
    total: int
    hits: list[LiteratureHit]
    note: str


class StudyVariantsResult(BaseModel):
    pmid: str
    total: int = Field(description="Variants the study reported")
    carried: int = Field(description="How many this genome carries a non-reference allele at")
    markers: list[queries.MarkerGenotype]
    note: str


def _require_db() -> str | None:
    if not db_exists():
        return "No Locus database found. Build it first with `locus pipeline`."
    return None


@mcp.tool()
def genome_overview() -> dict:
    """Summarize the loaded genome: variant counts, genome build, and how much is
    annotated (ClinVar, gnomAD, pharmacogenomics). Call this first to see what's available."""
    err = _require_db()
    if err:
        return {"error": err}
    return queries.overview()


@mcp.tool()
def lookup_variant_by_rsid(rsid: str) -> queries.VariantPage:
    """Look up variant(s) by dbSNP rsID (e.g. 'rs1799853'). Returns the genotype and
    any ClinVar/gnomAD annotations carried for that variant."""
    return queries.lookup_by_rsid(rsid)


@mcp.tool()
def lookup_variants_in_gene(gene: str, limit: int = 100, offset: int = 0) -> queries.VariantPage:
    """List variants in a gene by HGNC symbol (e.g. 'BRCA1'). Paginated — pass offset to page.
    `total` is the full match count. Requires the SnpEff annotation step to have run."""
    return queries.lookup_by_gene(gene, limit=limit, offset=offset)


@mcp.tool()
def lookup_variants_in_region(region: str, limit: int = 200, offset: int = 0) -> queries.VariantPage:
    """List variants in a genomic region. `region` is 'chr7:117480000-117670000' (1-based, chr-prefixed)
    or a single position 'chr7:117559590'. Paginated."""
    return queries.lookup_by_region(region, limit=limit, offset=offset)


@mcp.tool()
def clinical_findings(gene: str = "", significance: str = "", limit: int = 100, offset: int = 0) -> queries.VariantPage:
    """ClinVar clinical findings in this genome. With no args, returns pathogenic / likely-pathogenic
    variants. Optionally filter by `gene` (HGNC symbol) and/or `significance` substring
    (e.g. 'pathogenic', 'risk_factor', 'drug_response'). Always confirm health-relevant hits clinically."""
    return queries.clinical_findings(
        gene=gene or None, significance=significance or None, limit=limit, offset=offset
    )


@mcp.tool()
def pharmacogenomics(gene: str = "", drug: str = "") -> queries.PgxResult:
    """Pharmacogenomic results from PharmCAT: star-allele diplotypes, metabolizer phenotypes, and
    CPIC/DPWG drug guidance. Filter by `gene` (e.g. 'CYP2C19') or `drug` (e.g. 'clopidogrel').
    Some genes report 'Unknown'/'No Result' — notably CYP2D6, whose calling needs copy-number data
    from raw reads that Locus does not process. Treat those as "not determined", not "normal".
    Requires the PharmCAT annotation step to have run."""
    return queries.pharmacogenomics(gene=gene or None, drug=drug or None)


@mcp.tool()
def allele_frequency(region: str) -> queries.VariantPage:
    """How common is a variant? Returns gnomAD allele frequencies for variants at a position/region
    ('chr1:55052000' or 'chr1:55050000-55060000'). Low gnomad_af = rare. Frequencies are carried for
    the variants where rarity is informative — those ClinVar has classified or AlphaMissense calls
    pathogenic; other variants return a null gnomad_af rather than a frequency."""
    return queries.allele_frequency(region)


@mcp.tool()
def structural_variants(region: str, limit: int = 100) -> StructuralVariantsResult:
    """Copy-number (CNV) and structural (SV) events overlapping a region ('chr1:1000000-2000000').
    For CNVs, `cn` is the estimated copy number (2 = normal diploid)."""
    hits = queries.structural_overlap(region, limit=limit)
    return StructuralVariantsResult(total=len(hits), hits=hits)


@mcp.tool()
def predicted_damaging(gene: str = "", limit: int = 100) -> queries.VariantPage:
    """Rare, predicted-damaging missense variants from AlphaMissense — the 'ClinVar is silent' set:
    variants ClinVar has NOT classified but AlphaMissense scores as likely pathogenic. Optionally
    filter by gene. Use this when ClinVar returns nothing for a gene of interest. These are
    computational predictions, not clinical assertions — far weaker evidence than a ClinVar
    pathogenic call, and most are benign in reality. Rarity filtering is best-effort: gnomAD
    frequencies are only carried for some variants, so a hit is not guaranteed to be rare."""
    err = _require_db()
    if err:
        return {"error": err}  # type: ignore[return-value]
    return queries.predicted_damaging(gene=gene or None, limit=limit)


@mcp.tool()
def secondary_findings(limit: int = 100) -> queries.VariantPage:
    """ACMG SF v3.3 secondary (incidental) findings: pathogenic / likely-pathogenic ClinVar variants
    in the 84 medically-actionable genes ACMG recommends reporting (cancer, cardiac, metabolic).
    ACMG's recessive genes are reported only when two P/LP variants are present, so single carriers
    are excluded by design. Empty is the common, reassuring result — but it means "no P/LP ClinVar
    variant in these 84 genes", NOT a clean bill of health: it cannot see variants ClinVar hasn't
    classified, non-coding or structural variants, or any gene off this list. Confirm hits clinically."""
    err = _require_db()
    if err:
        return queries.VariantPage(total=0, limit=limit, offset=0, hits=[])
    return queries.secondary_findings(limit=limit)


@mcp.tool()
def carrier_status(limit: int = 100) -> queries.CarrierReport:
    """Carrier status for common recessive conditions — which ones this person carries ONE
    pathogenic copy of. This is about their children, not their own health: a carrier is
    typically unaffected. It matters when both partners carry the same condition (1-in-4 risk per
    pregnancy), so frame results around family planning and a genetic counselor, never as a
    personal diagnosis. 'likely_affected' = two pathogenic copies (phase unknown from a VCF, so a
    compound-het is presumed, not proven).

    IMPORTANT: always report the `not_assessed` list alongside any result. This is a curated panel
    of VCF-assessable conditions, not a clinical carrier screen — and several of the most
    important screens (SMN1/spinal muscular atrophy, FMR1/Fragile X) CANNOT be answered from this
    data at all. An empty `hits` list is NOT a negative carrier screen; do not present it as one."""
    err = _require_db()
    if err:
        return queries.CarrierReport(total=0, hits=[], not_assessed=[], panel_size=0, note=err)
    return queries.carrier_status(limit=limit)


@mcp.tool()
def traits(category: str = "") -> queries.TraitsReport:
    """Single-SNP traits/wellness (lactose, caffeine, alcohol flush, earwax, eye color, muscle type),
    the HLA-B*57:01 abacavir-hypersensitivity screening proxy, and the mtDNA maternal-lineage
    haplogroup. Filter by `category`: 'wellness' | 'pharmacogenomic' | 'maternal lineage' (the
    haplogroup lives under the last one). Call with no filter to see everything. Informational,
    not diagnostic. Requires `locus traits` to have run."""
    err = _require_db()
    if err:
        return queries.TraitsReport(total=0, traits=[])
    return queries.traits(category=category or None)


@mcp.tool()
def ancestry() -> queries.AncestrySummary:
    """Estimated biogeographic ancestry: continental proportions (k-NN over 1000 Genomes + HGDP) and the
    PCA placement. Continental ancestry is reliable; finer/admixed breakdowns are approximate.
    Requires `locus ancestry` to have run."""
    err = _require_db()
    if err:
        return {"error": err}  # type: ignore[return-value]
    return queries.ancestry()


@mcp.tool()
def polygenic_risk() -> PolygenicRiskReport:
    """Polygenic (aggregate) risk scores for the curated set (coronary artery disease, LDL
    cholesterol, type 2 diabetes, Lp(a)) — call it to see what is actually scored, reported as
    an ancestry-matched percentile where available. Percentiles are only meaningful within the matched
    ancestry; these are research-grade estimates, not diagnoses. Requires `locus ancestry` to have run."""
    err = _require_db()
    if err:
        return PolygenicRiskReport(scores=[])
    return PolygenicRiskReport(scores=queries.polygenic_risk())


@mcp.tool()
def gwas_associations(trait: str = "", limit: int = 100) -> queries.AssociationPage:
    """GWAS Catalog risk alleles THIS genome carries (genome-wide significant, p<5e-8), optionally
    filtered by `trait` substring (e.g. 'type 2 diabetes', 'height'). WEAK / EXPLORATORY: single hits
    with tiny effects — do NOT sum them or read them like a calibrated score (use polygenic_risk for
    that). Ordered by significance. Requires `locus gwas` to have run."""
    err = _require_db()
    if err:
        return queries.AssociationPage(total=0, limit=limit, offset=0, hits=[])
    return queries.gwas_associations(trait=trait or None, limit=limit)


@mcp.tool()
def whats_new(since: str = "", tier: str = "") -> queries.WhatsNew:
    """What changed about THIS genome since the last `locus refresh` — the deterministic
    changelog (e.g. ClinVar variants you carry that were newly classified pathogenic or
    reclassified). Ranked strongest-first. Optionally filter by `since` (ISO date like
    '2026-06-01') or `tier` ('strong'|'moderate'|'weak'|'info'). Empty if refresh hasn't run."""
    err = _require_db()
    if err:
        return queries.WhatsNew(total=0, counts_by_tier={}, findings=[])
    return queries.whats_new(since=since or None, tier=tier or None)


@mcp.tool()
def ask_about(query: str, limit: int = 50) -> queries.AskResult:
    """On-demand: given a condition/trait name OR a paper's rsIDs (space/comma-separated), report what
    THIS genome carries. rsIDs → live genotype lookup at those positions (hom-ref-aware) plus any
    ClinVar/AlphaMissense we already have; a trait name → your carried GWAS associations for it. This is
    the 'keep up with new studies' tool: paste rsIDs from a new paper to see your genotypes. Weak
    single-variant evidence — informational, not diagnostic."""
    import re as _re

    err = _require_db()
    if err:
        return queries.AskResult(query=query, mode="error", markers=[], associations=[], note=err)
    rsids = _re.findall(r"rs\d+", query, _re.IGNORECASE)
    if rsids:
        import contextlib
        import sys

        from . import gwas

        # markers_genotypes shells out and logs to stdout; redirect to stderr so it can't
        # corrupt the stdio JSON-RPC stream (safe — FastMCP isn't writing during this sync call).
        with contextlib.redirect_stdout(sys.stderr):
            rows = gwas.ask_markers(rsids)
        markers = [queries.MarkerGenotype(**r) for r in rows]
        note = ("Live genotype lookup at the requested SNPs (Ensembl GRCh38 coords, hom-ref-aware). "
                "Single-variant evidence — confirm context before interpreting."
                if markers else "Could not resolve those rsIDs (Ensembl lookup failed or unknown IDs).")
        return queries.AskResult(query=query, mode="rsids", markers=markers, associations=[], note=note)
    page = queries.gwas_associations(trait=query, limit=limit)
    note = page.note if page.hits else ("No carried GWAS associations for that trait. Run `locus gwas` "
                                        "first, or pass rsIDs to genotype specific variants.")
    return queries.AskResult(query=query, mode="trait", markers=[], associations=page.hits, note=note)


@mcp.tool()
def literature_for(query: str, since: str = "") -> LiteratureResult:
    """Latest PubMed research on a gene ('BRCA2'), an rsID ('rs7903146'), or a topic, to read in
    the context of THIS genome. Returns raw paper titles + abstracts + PubMed links for you to
    summarize and relate to the person's variants — it does NOT itself filter by genotype. Optionally
    pass `since` (ISO date like '2026-01-01') to limit to recent papers. Informational, not diagnostic."""
    import contextlib
    import sys

    from . import literature

    # pubmed_search logs to stdout; keep it off the stdio JSON-RPC stream.
    with contextlib.redirect_stdout(sys.stderr):
        raw = literature.literature_for(query, since=since or None)
    hits = [LiteratureHit(pmid=h.pmid, title=h.title, abstract=h.abstract, journal=h.journal,
                          year=h.year, gene=h.gene, url=h.url) for h in raw]
    note = ("Recent PubMed papers matching the query. Summarize them against the genome's other "
            "findings; single papers are not clinical guidance." if hits
            else "No PubMed papers found (or the lookup failed). Try a broader query or drop `since`.")
    return LiteratureResult(query=query, total=len(hits), hits=hits, note=note)


@mcp.tool()
def variants_in_study(pmid: str) -> StudyVariantsResult:
    """Given the PubMed ID of a GWAS study, pull the variants it reported (via the GWAS Catalog) and
    show which ones THIS genome carries — live, hom-ref-aware genotyping. Answers 'this new paper
    found N variants, which do I have?'. Weak single-variant evidence — informational, not diagnostic."""
    err = _require_db()
    if err:
        return StudyVariantsResult(pmid=pmid, total=0, carried=0, markers=[], note=err)
    import contextlib
    import sys

    from . import literature

    with contextlib.redirect_stdout(sys.stderr):
        res = literature.study_variants(pmid)
    markers = [queries.MarkerGenotype(**m) for m in res["markers"]]
    return StudyVariantsResult(pmid=pmid, total=res["total"], carried=res["carried"],
                               markers=markers, note=res["note"])


@mcp.tool()
def run_sql(query: str) -> dict:
    """Run a read-only SELECT against the genome database for power queries. Tables: variants
    (chrom,pos,ref,alt,rsid,gt,filter,gene,consequence,clnsig,clndn,clnrevstat,gnomad_af,gnomad_af_grpmax,
    am_pathogenicity,am_class), cnv, sv, pgx_genes, pgx_drugs, meta, traits, associations (carried GWAS
    hits), pgs_scores, ancestry_global, ancestry_pca, watch_findings (the refresh changelog),
    watch_seen_ids, sources. Mutating statements are rejected; results capped at 200 rows."""
    err = _require_db()
    if err:
        return {"error": err}
    try:
        return queries.run_sql(query)
    except ValueError as e:
        return {"error": str(e)}


def main() -> None:
    mcp.run()  # stdio transport (default)


if __name__ == "__main__":
    main()

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
    `total` is the full match count. Requires gene annotation (VEP/SnpEff) to be loaded."""
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
    Requires the PharmCAT annotation step to have run."""
    return queries.pharmacogenomics(gene=gene or None, drug=drug or None)


@mcp.tool()
def allele_frequency(region: str) -> queries.VariantPage:
    """How common is a variant? Returns gnomAD population allele frequencies for variants at a
    position/region ('chr1:55052000' or 'chr1:55050000-55060000'). Low gnomad_af = rare."""
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
    variants ClinVar never classified but AlphaMissense scores as likely pathogenic (and rare, AF<1%).
    Optionally filter by gene. Use this when ClinVar returns nothing for a gene of interest."""
    err = _require_db()
    if err:
        return {"error": err}  # type: ignore[return-value]
    return queries.predicted_damaging(gene=gene or None, limit=limit)


@mcp.tool()
def ancestry() -> queries.AncestrySummary:
    """Estimated biogeographic ancestry: continental proportions (k-NN over 1000 Genomes) and the
    PCA placement. Continental ancestry is reliable; finer/admixed breakdowns are approximate.
    Requires `locus ancestry` to have run."""
    err = _require_db()
    if err:
        return {"error": err}  # type: ignore[return-value]
    return queries.ancestry()


@mcp.tool()
def polygenic_risk() -> PolygenicRiskReport:
    """Polygenic (aggregate) risk scores for common traits (CAD, LDL, T2D, AFib, Lp(a)), reported as
    an ancestry-matched percentile where available. Percentiles are only meaningful within the matched
    ancestry; these are research-grade estimates, not diagnoses. Requires `locus ancestry` to have run."""
    err = _require_db()
    if err:
        return PolygenicRiskReport(scores=[])
    return PolygenicRiskReport(scores=queries.polygenic_risk())


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
def run_sql(query: str) -> dict:
    """Run a read-only SELECT against the genome database for power queries. Tables: variants
    (chrom,pos,ref,alt,rsid,gt,filter,gene,consequence,clnsig,clndn,clnrevstat,gnomad_af,gnomad_af_grpmax),
    cnv, sv, pgx_genes, pgx_drugs, meta. Mutating statements are rejected; results capped at 200 rows."""
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

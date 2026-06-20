"""Shared read-only query layer over the Locus DuckDB store.

Both the MCP server (``mcp_server.py``) and the SPA backend (``api.py``) call
these functions, so query logic lives in exactly one place. Everything here is
read-only and returns Pydantic models, which gives the MCP tools a precise
output schema and the API clean JSON.
"""

from __future__ import annotations

import re

from pydantic import BaseModel, Field

from .db import connect

# Columns surfaced for a variant, in a fixed order.
_VAR_COLS = (
    "chrom", "pos", "ref", "alt", "rsid", "gt", "filter", "gene", "consequence",
    "clnsig", "clndn", "clnrevstat", "gnomad_af", "gnomad_af_grpmax",
    "am_pathogenicity", "am_class",
)


class Variant(BaseModel):
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str | None = None
    gt: str | None = Field(default=None, description="Genotype, e.g. 0/1 (het) or 1/1 (hom-alt)")
    filter: str | None = None
    gene: str | None = None
    consequence: str | None = None
    clnsig: str | None = Field(default=None, description="ClinVar clinical significance")
    clndn: str | None = Field(default=None, description="ClinVar disease name(s)")
    clnrevstat: str | None = Field(default=None, description="ClinVar review status (star rating)")
    gnomad_af: float | None = Field(default=None, description="gnomAD global allele frequency")
    gnomad_af_grpmax: float | None = Field(default=None, description="gnomAD max per-ancestry AF")
    am_pathogenicity: float | None = Field(default=None, description="AlphaMissense pathogenicity (0-1)")
    am_class: str | None = Field(default=None, description="AlphaMissense class (benign/ambiguous/pathogenic)")


class VariantPage(BaseModel):
    total: int = Field(description="Total matching rows (before limit/offset)")
    limit: int
    offset: int
    hits: list[Variant]


class PgxGene(BaseModel):
    gene: str
    diplotype: str | None = None
    phenotype: str | None = None
    activity_score: str | None = None


class PgxDrug(BaseModel):
    drug: str
    gene: str | None = None
    source: str | None = None
    recommendation: str | None = None


class PgxResult(BaseModel):
    genes: list[PgxGene]
    drugs: list[PgxDrug]


class StructuralHit(BaseModel):
    kind: str  # "cnv" | "sv"
    chrom: str
    pos: int
    end: int | None = None
    svtype: str | None = None
    cn: int | None = None
    svlen: int | None = None
    filter: str | None = None
    genes: str | None = None


def _rows_to_variants(rows) -> list[Variant]:
    return [Variant(**dict(zip(_VAR_COLS, r, strict=True))) for r in rows]


def _var_select(where: str, params: list, limit: int, offset: int) -> VariantPage:
    cols = ", ".join(_VAR_COLS)
    with connect(read_only=True) as con:
        total = con.execute(f"SELECT count(*) FROM variants WHERE {where}", params).fetchone()[0]
        rows = con.execute(
            f"SELECT {cols} FROM variants WHERE {where} ORDER BY chrom, pos LIMIT ? OFFSET ?",
            [*params, limit, offset],
        ).fetchall()
    return VariantPage(total=total, limit=limit, offset=offset, hits=_rows_to_variants(rows))


# ── Region parsing ────────────────────────────────────────────────────────────
_REGION_RE = re.compile(r"^(chr[\w]+|[\w]+):([\d,]+)(?:-([\d,]+))?$", re.IGNORECASE)


def parse_region(region: str) -> tuple[str, int, int]:
    """Parse 'chr7:117,480,000-117,670,000' or 'chr7:117559590' into (chrom, start, end)."""
    m = _REGION_RE.match(region.strip())
    if not m:
        raise ValueError(f"Unrecognized region '{region}'. Use chr:start-end or chr:pos.")
    chrom = m.group(1)
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", "")) if m.group(3) else start
    return chrom, start, end


# ── Public query functions ─────────────────────────────────────────────────────
def lookup_by_rsid(rsid: str, limit: int = 50, offset: int = 0) -> VariantPage:
    rsid = rsid.strip()
    if not rsid.lower().startswith("rs"):
        rsid = "rs" + rsid
    return _var_select("rsid = ?", [rsid], limit, offset)


def lookup_by_gene(gene: str, limit: int = 100, offset: int = 0) -> VariantPage:
    return _var_select("upper(gene) = upper(?)", [gene.strip()], limit, offset)


def lookup_by_region(region: str, limit: int = 200, offset: int = 0) -> VariantPage:
    chrom, start, end = parse_region(region)
    return _var_select("chrom = ? AND pos BETWEEN ? AND ?", [chrom, start, end], limit, offset)


def clinical_findings(
    gene: str | None = None, significance: str | None = None, limit: int = 100, offset: int = 0
) -> VariantPage:
    """ClinVar-annotated variants. Defaults to pathogenic / likely-pathogenic."""
    where = ["clnsig IS NOT NULL"]
    params: list = []
    if significance:
        where.append("lower(clnsig) LIKE ?")
        params.append(f"%{significance.lower()}%")
    else:
        # Default to clean pathogenic/likely-pathogenic: exclude benign and the
        # lower-confidence "Conflicting_classifications_of_pathogenicity".
        where.append(
            "(lower(clnsig) LIKE '%pathogenic%' "
            "AND lower(clnsig) NOT LIKE '%benign%' "
            "AND lower(clnsig) NOT LIKE '%conflicting%')"
        )
    if gene:
        where.append("upper(gene) = upper(?)")
        params.append(gene.strip())
    return _var_select(" AND ".join(where), params, limit, offset)


def predicted_damaging(
    gene: str | None = None, max_af: float = 0.01, limit: int = 100, offset: int = 0
) -> VariantPage:
    """Rare, predicted-damaging missense variants (AlphaMissense pathogenic + rare/unannotated).

    This is the 'ClinVar is silent' set — variants ClinVar has never classified but AlphaMissense
    scores as likely damaging. Defaults to allele frequency < 1%.
    """
    where = ["am_class LIKE '%pathogenic%'", "(gnomad_af IS NULL OR gnomad_af < ?)"]
    params: list = [max_af]
    if gene:
        where.append("upper(gene) = upper(?)")
        params.append(gene.strip())
    return _var_select(" AND ".join(where), params, limit, offset)


def allele_frequency(region_or_variant: str) -> VariantPage:
    """How common is a variant? Looks up gnomAD AF for a position/region."""
    return lookup_by_region(region_or_variant, limit=50, offset=0)


def pharmacogenomics(gene: str | None = None, drug: str | None = None) -> PgxResult:
    with connect(read_only=True) as con:
        gwhere, gparams = ("TRUE", [])
        if gene:
            gwhere, gparams = ("upper(gene) = upper(?)", [gene.strip()])
        genes = con.execute(
            f"SELECT gene, diplotype, phenotype, activity_score FROM pgx_genes WHERE {gwhere} ORDER BY gene",
            gparams,
        ).fetchall()
        dwhere, dparams = ("TRUE", [])
        conds = []
        if drug:
            conds.append("lower(drug) LIKE ?")
            dparams.append(f"%{drug.lower()}%")
        if gene:
            conds.append("upper(gene) LIKE upper(?)")
            dparams.append(f"%{gene.strip()}%")
        if conds:
            dwhere = " AND ".join(conds)
        drugs = con.execute(
            f"SELECT drug, gene, source, recommendation FROM pgx_drugs WHERE {dwhere} ORDER BY drug",
            dparams,
        ).fetchall()
    return PgxResult(
        genes=[PgxGene(gene=g[0], diplotype=g[1], phenotype=g[2], activity_score=g[3]) for g in genes],
        drugs=[PgxDrug(drug=d[0], gene=d[1], source=d[2], recommendation=d[3]) for d in drugs],
    )


def structural_overlap(region: str, limit: int = 100) -> list[StructuralHit]:
    """CNV/SV records overlapping a region.

    Matches the chromosome with *or* without a ``chr`` prefix, so it is correct
    whether the cnv/sv tables were loaded with canonicalized contigs (the current
    loader) or the older non-prefixed form.
    """
    chrom, start, end = parse_region(region)
    stem = chrom[3:] if chrom.startswith("chr") else chrom
    chrom_opts = [chrom, stem]
    if stem == "M":
        chrom_opts.append("MT")  # mitochondrion: chrM <-> MT
    placeholders = ", ".join("?" for _ in chrom_opts)
    hits: list[StructuralHit] = []
    with connect(read_only=True) as con:
        for kind, tbl in (("cnv", "cnv"), ("sv", "sv")):
            cn_col = "cn" if kind == "cnv" else "NULL AS cn"
            rows = con.execute(
                f"""SELECT chrom, pos, "end", svtype, {cn_col}, svlen, filter, genes
                    FROM {tbl}
                    WHERE chrom IN ({placeholders}) AND pos <= ? AND "end" >= ? LIMIT ?""",
                [*chrom_opts, end, start, limit],
            ).fetchall()
            for r in rows:
                hits.append(StructuralHit(
                    kind=kind, chrom=r[0], pos=r[1], end=r[2], svtype=r[3],
                    cn=r[4], svlen=r[5], filter=r[6], genes=r[7],
                ))
    return hits


_FORBIDDEN_SQL = re.compile(
    r"\b(insert|update|delete|drop|alter|create|attach|copy|pragma|install|load|export)\b", re.IGNORECASE
)


def run_sql(sql: str, max_rows: int = 200) -> dict:
    """Run a read-only SELECT against the store. Rejects any mutating statement."""
    s = sql.strip().rstrip(";")
    if not re.match(r"^\s*(select|with)\b", s, re.IGNORECASE):
        raise ValueError("Only SELECT/WITH queries are allowed.")
    if _FORBIDDEN_SQL.search(s):
        raise ValueError("Query contains a disallowed keyword.")
    with connect(read_only=True) as con:
        cur = con.execute(s)
        names = [d[0] for d in cur.description]
        rows = cur.fetchmany(max_rows)
    return {"columns": names, "rows": [list(r) for r in rows], "truncated_to": max_rows}


class AncestryComponent(BaseModel):
    code: str
    name: str
    proportion: float = Field(description="Fraction (0-1); k-NN estimate over 1000 Genomes")


class PcaPoint(BaseModel):
    label: str
    pc1: float
    pc2: float
    is_sample: bool


class AncestrySummary(BaseModel):
    components: list[AncestryComponent]      # continental rollup
    populations: list[AncestryComponent]     # sub-continental (fine 1000 Genomes populations)
    pca: list[PcaPoint]
    note: str = (
        "k-NN placement among 1000 Genomes populations. Continental ancestry is robust; the "
        "sub-continental breakdown is 'genetically closest to' (sensitive to reference panel sizes), "
        "not a calibrated admixture percentage — and far coarser than 23andMe's proprietary panels."
    )


class PgsResult(BaseModel):
    pgs_id: str
    trait: str
    raw: float
    percentile: float | None = Field(default=None, description="Within the ancestry-matched reference")
    ancestry: str | None = None
    n_used: int
    coverage: float


def ancestry() -> AncestrySummary:
    with connect(read_only=True) as con:
        cont = con.execute(
            "SELECT code, name, proportion FROM ancestry_global WHERE level='continental' "
            "ORDER BY proportion DESC"
        ).fetchall()
        pops = con.execute(
            "SELECT code, name, proportion FROM ancestry_global WHERE level='population' "
            "ORDER BY proportion DESC"
        ).fetchall()
        pca = con.execute("SELECT label, pc1, pc2, is_sample FROM ancestry_pca").fetchall()
    return AncestrySummary(
        components=[AncestryComponent(code=c[0], name=c[1], proportion=c[2]) for c in cont],
        populations=[AncestryComponent(code=c[0], name=c[1], proportion=c[2]) for c in pops],
        pca=[PcaPoint(label=p[0], pc1=p[1], pc2=p[2], is_sample=p[3]) for p in pca],
    )


def polygenic_risk() -> list[PgsResult]:
    with connect(read_only=True) as con:
        rows = con.execute(
            "SELECT pgs_id, trait, raw, percentile, ancestry, n_used, coverage "
            "FROM pgs_scores ORDER BY trait"
        ).fetchall()
    return [
        PgsResult(pgs_id=r[0], trait=r[1], raw=r[2], percentile=r[3], ancestry=r[4],
                  n_used=r[5], coverage=r[6])
        for r in rows
    ]


def overview() -> dict:
    """Summary stats about the loaded genome (counts, build, annotation coverage)."""
    with connect(read_only=True) as con:
        meta = dict(con.execute("SELECT key, value FROM meta").fetchall())
        n = con.execute("SELECT count(*) FROM variants").fetchone()[0]
        annotated = con.execute("SELECT count(*) FROM variants WHERE clnsig IS NOT NULL").fetchone()[0]
        with_af = con.execute("SELECT count(*) FROM variants WHERE gnomad_af IS NOT NULL").fetchone()[0]
        pgx = con.execute("SELECT count(*) FROM pgx_genes").fetchone()[0]
        cnv = con.execute("SELECT count(*) FROM cnv").fetchone()[0]
        sv = con.execute("SELECT count(*) FROM sv").fetchone()[0]
    return {
        "meta": meta,
        "variants": n,
        "clinvar_annotated": annotated,
        "gnomad_annotated": with_af,
        "pgx_genes": pgx,
        "cnv": cnv,
        "sv": sv,
    }

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

    The 'ClinVar is silent' set — variants ClinVar has never classified but AlphaMissense scores as
    likely damaging. ``clnsig IS NULL`` is what makes that true and is not optional: without it a
    third of the results are variants ClinVar *has* seen, including ones it calls **benign**, and
    they get reported as damaging findings ClinVar supposedly missed.

    Rarity (default AF < 1%) is best-effort: gnomAD AF is only carried for variants where rarity is
    informative, so `gnomad_af IS NULL` is admitted rather than dropping unscored variants entirely.
    """
    where = ["am_class LIKE '%pathogenic%'", "clnsig IS NULL", "(gnomad_af IS NULL OR gnomad_af < ?)"]
    params: list = [max_af]
    if gene:
        where.append("upper(gene) = upper(?)")
        params.append(gene.strip())
    return _var_select(" AND ".join(where), params, limit, offset)


def allele_frequency(region_or_variant: str) -> VariantPage:
    """How common is a variant? Looks up gnomAD AF for a position/region."""
    return lookup_by_region(region_or_variant, limit=50, offset=0)


def secondary_findings(limit: int = 100, offset: int = 0) -> VariantPage:
    """ACMG SF secondary findings: pathogenic / likely-pathogenic ClinVar variants in the
    medically-actionable gene set. 'No findings' is a defensible, reassuring result.

    ACMG reports its recessive genes only when TWO P/LP variants are present, so a lone
    heterozygous carrier is excluded — otherwise common carrier states (HFE p.C282Y sits in ~10%
    of Europeans) would surface as actionable findings and generate false alarms.
    """
    from .panels import ACMG_SF_GENES, ACMG_SF_RECESSIVE

    genes = sorted(ACMG_SF_GENES)
    placeholders = ", ".join("?" for _ in genes)
    plp = ("(lower(clnsig) LIKE '%pathogenic%' AND lower(clnsig) NOT LIKE '%benign%' "
           "AND lower(clnsig) NOT LIKE '%conflicting%')")
    rec = sorted(ACMG_SF_RECESSIVE)
    rec_ph = ", ".join("?" for _ in rec)
    # Dominant genes: any P/LP. Recessive genes: only if biallelic — homozygous (gt not 0/x),
    # or ≥2 P/LP variants in that same gene (presumed compound het; phase is unknown from a VCF).
    where = (
        f"{plp} AND upper(gene) IN ({placeholders}) AND ("
        f"  upper(gene) NOT IN ({rec_ph})"
        f"  OR gt IN ('1/1', '1|1')"
        f"  OR (SELECT count(*) FROM variants v2 WHERE upper(v2.gene) = upper(variants.gene)"
        f"      AND {plp.replace('clnsig', 'v2.clnsig')}) >= 2"
        f")"
    )
    params = [g.upper() for g in genes] + [g.upper() for g in rec]
    return _var_select(where, params, limit, offset)


class CarrierHit(BaseModel):
    gene: str
    condition: str
    inheritance: str = Field(description="AR (autosomal recessive) | XL (X-linked)")
    status: str = Field(description="carrier (one copy) | likely_affected (two copies)")
    n_variants: int = Field(description="Pathogenic/likely-pathogenic variants found in this gene")
    zygosity: str | None = Field(default=None, description="Genotype of the top variant, e.g. 0/1")
    rsid: str | None = None
    clnsig: str | None = None
    clndn: str | None = Field(default=None, description="ClinVar's disease name for the variant")


class NotAssessed(BaseModel):
    gene: str
    condition: str
    why: str = Field(description="Why this genome's data cannot answer it")


class CarrierReport(BaseModel):
    total: int
    hits: list[CarrierHit]
    not_assessed: list[NotAssessed] = Field(
        description="Conditions this data CANNOT speak to — absence here is not a negative result")
    panel_size: int
    note: str


def carrier_status(limit: int = 100) -> CarrierReport:
    """Carrier status for common recessive conditions: which you carry one pathogenic copy of.

    The complement of ``secondary_findings``, which drops lone heterozygous carriers because they
    aren't findings *for you*. They are exactly what matters for family planning: two carriers of
    the same condition have a 1-in-4 risk per pregnancy.

    Always returns ``not_assessed`` alongside the hits. Several of the most important carrier
    tests (SMN1, FMR1) are simply not answerable from a VCF, and an empty result that quietly
    implied they had been checked would be the most dangerous thing this function could do.
    """
    from .panels import CARRIER_PANEL, CARRIER_UNASSESSABLE

    by_gene = {c.gene: c for c in CARRIER_PANEL}
    genes = sorted(by_gene)
    ph = ", ".join("?" for _ in genes)
    plp = ("lower(clnsig) LIKE '%pathogenic%' AND lower(clnsig) NOT LIKE '%benign%' "
           "AND lower(clnsig) NOT LIKE '%conflicting%'")
    with connect(read_only=True) as con:
        rows = con.execute(
            f"SELECT upper(gene), count(*), max(CASE WHEN gt IN ('1/1','1|1') THEN 1 ELSE 0 END), "
            f"       any_value(gt), any_value(rsid), any_value(clnsig), any_value(clndn) "
            f"FROM variants WHERE {plp} AND upper(gene) IN ({ph}) GROUP BY upper(gene)",
            [g.upper() for g in genes],
        ).fetchall()

    hits: list[CarrierHit] = []
    for gene, n, has_hom, gt, rsid, clnsig, clndn in rows:
        info = by_gene[gene]
        # Two copies (homozygous, or ≥2 P/LP variants — phase unknown from a VCF, so compound-het
        # is presumed rather than proven) means possibly affected, not merely a carrier.
        affected = bool(has_hom) or n >= 2
        hits.append(CarrierHit(
            gene=gene, condition=info.condition, inheritance=info.inheritance,
            status="likely_affected" if affected else "carrier",
            n_variants=n, zygosity=gt, rsid=rsid, clnsig=clnsig,
            clndn=(clndn or "").replace("_", " ") or None))
    hits.sort(key=lambda h: (h.status != "likely_affected", h.gene))

    note = (
        "Carrier status is about your children, not your health: one copy of a recessive variant "
        "is typically silent for you. It matters when both partners carry the same condition "
        "(1-in-4 risk per pregnancy), so interpret it with a partner's results and a genetic "
        "counselor. 'likely_affected' means two pathogenic copies were seen — phase is unknown "
        "from a VCF, so a compound-het is presumed, not proven; confirm clinically. "
        "This is a curated panel of common, VCF-assessable conditions — NOT a clinical carrier "
        "screen (ACMG's Tier 3 panel is 113 genes). An empty result means nothing was found in "
        "these genes; it is not a negative screen, and it says nothing about `not_assessed`."
    )
    return CarrierReport(
        total=len(hits), hits=hits[:limit], panel_size=len(CARRIER_PANEL), note=note,
        not_assessed=[NotAssessed(gene=g, condition=c, why=w) for g, c, w in CARRIER_UNASSESSABLE],
    )


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
    proportion: float = Field(description="Fraction (0-1); k-NN estimate over 1000 Genomes + HGDP")


class PcaPoint(BaseModel):
    label: str
    pc1: float
    pc2: float
    is_sample: bool
    group: str | None = Field(default=None, description="Superpopulation, for coloring")


class AncestrySummary(BaseModel):
    components: list[AncestryComponent]      # continental rollup
    populations: list[AncestryComponent]     # sub-continental (fine 1000 Genomes + HGDP populations)
    pca: list[PcaPoint]
    note: str = (
        "k-NN placement among 1000 Genomes + HGDP reference populations (so fine labels include HGDP "
        "groups such as French, Orcadian, Sardinian, Basque). Continental ancestry is robust; the "
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
        pca = con.execute('SELECT label, pc1, pc2, is_sample, "group" FROM ancestry_pca').fetchall()
    return AncestrySummary(
        components=[AncestryComponent(code=c[0], name=c[1], proportion=c[2]) for c in cont],
        populations=[AncestryComponent(code=c[0], name=c[1], proportion=c[2]) for c in pops],
        pca=[PcaPoint(label=p[0], pc1=p[1], pc2=p[2], is_sample=p[3], group=p[4]) for p in pca],
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


class WatchFinding(BaseModel):
    ts: str = Field(description="When this finding was recorded (ISO8601)")
    source: str
    kind: str = Field(description="newly_pathogenic | reclassified | withdrawn | depathogenized | release")
    tier: str = Field(description="strong | moderate | weak | info — confidence/actionability")
    title: str
    detail: str | None = None
    chrom: str | None = None
    pos: int | None = None
    gene: str | None = None
    rsid: str | None = None
    old_value: str | None = Field(default=None, description="Prior classification/value")
    new_value: str | None = Field(default=None, description="New classification/value")
    release: str | None = None
    url: str | None = Field(default=None, description="Citation / source link (e.g. PubMed, GWAS Catalog)")


class WhatsNew(BaseModel):
    total: int
    since: str | None = None
    counts_by_tier: dict[str, int]
    findings: list[WatchFinding]
    note: str = (
        "Deterministic changelog from `locus refresh`. 'strong' = high-confidence ClinVar "
        "reanalysis (multi-submitter/expert-panel); always confirm health-relevant hits clinically."
    )


_WATCH_COLS = ("ts", "source", "kind", "tier", "chrom", "pos", "rsid", "gene",
               "title", "detail", "old_value", "new_value", "release", "url")


def whats_new(since: str | None = None, tier: str | None = None, limit: int = 200) -> WhatsNew:
    """Ranked 'what changed about your genome' findings written by `locus refresh`.

    Optionally filter by `since` (ISO date/datetime) and/or `tier`. Ordered strongest-first.
    """
    with connect(read_only=True) as con:
        exists = con.execute(
            "SELECT count(*) FROM information_schema.tables WHERE table_name = 'watch_findings'"
        ).fetchone()[0]
        if not exists:
            return WhatsNew(total=0, since=since, counts_by_tier={}, findings=[])
        # `url` may be absent on a pre-v4 store we can't ALTER from a read-only connection.
        have = {c for (c,) in con.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name = 'watch_findings'"
        ).fetchall()}
        cols = tuple(c for c in _WATCH_COLS if c in have)
        where, params = ["TRUE"], []
        if since:
            where.append("ts >= ?")
            params.append(since)
        if tier:
            where.append("tier = ?")
            params.append(tier)
        clause = " AND ".join(where)
        total = con.execute(f"SELECT count(*) FROM watch_findings WHERE {clause}", params).fetchone()[0]
        counts = dict(con.execute(
            f"SELECT tier, count(*) FROM watch_findings WHERE {clause} GROUP BY tier", params
        ).fetchall())
        rank = "CASE tier WHEN 'strong' THEN 0 WHEN 'moderate' THEN 1 WHEN 'weak' THEN 2 ELSE 3 END"
        rows = con.execute(
            f"SELECT {', '.join(cols)} FROM watch_findings WHERE {clause} ORDER BY {rank}, ts DESC LIMIT ?",
            [*params, limit],
        ).fetchall()
    findings = [WatchFinding(**dict(zip(cols, r, strict=True))) for r in rows]
    return WhatsNew(total=total, since=since, counts_by_tier=counts, findings=findings)


class TraitResult(BaseModel):
    rsid: str
    category: str = Field(description="wellness | pharmacogenomic | maternal lineage (the mtDNA haplogroup)")
    trait: str
    genotype: str = Field(description="Observed genotype, e.g. 'A/G'; '—' if not callable")
    dosage: int | None = Field(default=None, description="Effect-allele copies (0/1/2)")
    effect_allele: str
    interpretation: str
    note: str | None = None


class TraitsReport(BaseModel):
    total: int
    traits: list[TraitResult]
    note: str = (
        "Single well-characterized SNPs — informational, not diagnostic. The HLA-B*57:01 entry "
        "is a European-validated screening proxy (confirm with HLA typing before acting)."
    )


def traits(category: str | None = None) -> TraitsReport:
    """Single-SNP trait/wellness results (and the HLA-B*57:01 proxy) from `locus traits`."""
    with connect(read_only=True) as con:
        exists = con.execute(
            "SELECT count(*) FROM information_schema.tables WHERE table_name = 'traits'"
        ).fetchone()[0]
        if not exists:
            return TraitsReport(total=0, traits=[])
        where, params = "TRUE", []
        if category:
            where, params = "category = ?", [category]
        rows = con.execute(
            f"SELECT rsid, category, trait, genotype, dosage, effect_allele, interpretation, note "
            f"FROM traits WHERE {where} ORDER BY category, trait", params
        ).fetchall()
    items = [TraitResult(rsid=r[0], category=r[1], trait=r[2], genotype=r[3], dosage=r[4],
                         effect_allele=r[5], interpretation=r[6], note=r[7]) for r in rows]
    return TraitsReport(total=len(items), traits=items)


class Association(BaseModel):
    rsid: str
    chrom: str
    pos: int
    risk_allele: str
    dosage: int = Field(description="Copies of the risk allele carried (1 or 2)")
    zygosity: str
    trait: str
    mapped_trait: str
    pval: float
    or_beta: str | None = None
    pmid: str | None = None


class AssociationPage(BaseModel):
    total: int
    limit: int
    offset: int
    trait: str | None = None
    hits: list[Association]
    note: str = (
        "GWAS Catalog risk alleles you carry — WEAK / EXPLORATORY single hits (genome-wide "
        "significant, p<5e-8), each with a tiny effect. Not a calibrated score (see polygenic_risk) "
        "and never sum these ORs. Ordered by significance."
    )


def gwas_associations(trait: str | None = None, limit: int = 100, offset: int = 0) -> AssociationPage:
    """GWAS Catalog risk alleles the genome carries, optionally filtered by trait substring."""
    with connect(read_only=True) as con:
        exists = con.execute(
            "SELECT count(*) FROM information_schema.tables WHERE table_name = 'associations'"
        ).fetchone()[0]
        if not exists:
            return AssociationPage(total=0, limit=limit, offset=offset, trait=trait, hits=[])
        where, params = ["TRUE"], []
        if trait:
            where.append("(lower(trait) LIKE ? OR lower(mapped_trait) LIKE ?)")
            params += [f"%{trait.lower()}%", f"%{trait.lower()}%"]
        clause = " AND ".join(where)
        total = con.execute(f"SELECT count(*) FROM associations WHERE {clause}", params).fetchone()[0]
        rows = con.execute(
            f"SELECT rsid, chrom, pos, risk_allele, dosage, zygosity, trait, mapped_trait, pval, "
            f"or_beta, pmid FROM associations WHERE {clause} ORDER BY pval ASC LIMIT ? OFFSET ?",
            [*params, limit, offset],
        ).fetchall()
    hits = [Association(rsid=r[0], chrom=r[1], pos=r[2], risk_allele=r[3], dosage=r[4],
                        zygosity=r[5], trait=r[6], mapped_trait=r[7], pval=r[8],
                        or_beta=r[9], pmid=r[10]) for r in rows]
    return AssociationPage(total=total, limit=limit, offset=offset, trait=trait, hits=hits)


class MarkerGenotype(BaseModel):
    rsid: str
    chrom: str
    pos: int
    genotype: str = Field(description="Observed genotype, e.g. 'A/G'; '—' if not callable")
    ref: str | None = Field(default=None, description="Reference allele (to tell hom-alt from hom-ref)")
    gene: str | None = None
    clnsig: str | None = None
    am_class: str | None = None


class AskResult(BaseModel):
    query: str
    mode: str = Field(description="rsids | trait | error")
    markers: list[MarkerGenotype]
    associations: list[Association]
    note: str


class VariantDossier(BaseModel):
    rsid: str
    found: bool = Field(description="False = no record in this genome's variant table")
    gene: str | None = None
    genotype: str | None = None
    zygosity: str | None = Field(default=None, description="heterozygous | homozygous-alternate | homozygous-reference")
    consequence: str | None = None
    clinvar: str | None = Field(default=None,
                                description="CLINICAL evidence — the strongest signal here")
    clinvar_disease: str | None = None
    clinvar_review: str | None = Field(default=None,
                                       description="ClinVar review status (stars) — how far to trust the call")
    alphamissense: str | None = Field(default=None,
                                      description="A COMPUTATIONAL prediction, far weaker than a ClinVar call")
    alphamissense_score: float | None = None
    gnomad_af: float | None = Field(default=None,
                                    description="Global allele frequency. Common = unlikely to be severe")
    gnomad_af_grpmax: float | None = None
    acmg_sf_gene: bool = Field(default=False,
                               description="Gene is on the ACMG SF v3.3 medically-actionable list")
    carrier_condition: str | None = Field(default=None,
                                          description="Recessive condition, if the gene is on the carrier panel")
    gwas_associations: list[Association] = Field(default_factory=list)
    trait: TraitResult | None = None
    literature_url: str
    note: str


def _zygosity(gt: str | None) -> str | None:
    if not gt:
        return None
    a = [x for x in re.split(r"[/|]", gt) if x not in ("", ".")]
    if len(a) < 2:
        return None
    if a[0] == a[1]:
        return "homozygous-reference" if a[0] == "0" else "homozygous-alternate"
    return "heterozygous"


def variant_dossier(rsid: str) -> VariantDossier:
    """Everything this genome knows about one variant, in a single call.

    Exists because the evidence types disagree and only mean something together: an
    AlphaMissense-"pathogenic" call looks alarming until you see ClinVar says benign and 40% of
    people carry it. Assembling that picture from separate tools invites a confident wrong answer,
    so the dossier puts the contradicting fields side by side and states the precedence.
    """
    from .panels import ACMG_SF_GENES, CARRIER_PANEL

    rsid = rsid.strip()
    if not rsid.lower().startswith("rs"):
        rsid = "rs" + rsid
    litvar = f"https://www.ncbi.nlm.nih.gov/research/litvar2/docsum?text={rsid}"

    page = lookup_by_rsid(rsid, limit=1)
    if not page.hits:
        return VariantDossier(
            rsid=rsid, found=False, literature_url=litvar,
            note=(f"{rsid} has no record in this genome's variant table. That almost always means "
                  f"homozygous-reference (the store keeps non-reference sites, so hom-ref calls "
                  f"simply aren't rows) — it does NOT mean the variant was ruled out or that the "
                  f"lookup failed. Use `ask_about` with this rsID to genotype the position live "
                  f"and be certain before telling anyone they 'don't have' it."))

    v = page.hits[0]
    gene = (v.gene or "").upper() or None
    assoc = gwas_associations(limit=25)
    carried = [a for a in assoc.hits if a.rsid == rsid]
    tr = next((t for t in traits().traits if t.rsid == rsid), None)
    cond = next((c.condition for c in CARRIER_PANEL if gene and c.gene == gene), None)

    note = (
        "Weigh these together — they routinely disagree, and the order matters. ClinVar is a "
        "clinical assertion and outranks everything else here; check its review status (a "
        "no-assertion call is weak). AlphaMissense is a computational prediction that is wrong "
        "often and NEVER overrides a ClinVar benign call. A high gnomAD frequency is strong "
        "evidence against a severe effect no matter what the prediction says: a variant carried "
        "by a large fraction of people is not causing a severe disease in all of them. GWAS hits "
        "are weak single associations and must not be summed. If the fields conflict, say so "
        "rather than picking the scariest one. Not diagnostic — confirm anything health-relevant "
        "with a clinician or genetic counselor."
    )
    return VariantDossier(
        rsid=rsid, found=True, gene=v.gene, genotype=v.gt, zygosity=_zygosity(v.gt),
        consequence=v.consequence, clinvar=v.clnsig,
        clinvar_disease=(v.clndn or "").replace("_", " ") or None, clinvar_review=v.clnrevstat,
        alphamissense=v.am_class, alphamissense_score=v.am_pathogenicity,
        gnomad_af=v.gnomad_af, gnomad_af_grpmax=v.gnomad_af_grpmax,
        acmg_sf_gene=bool(gene and gene in ACMG_SF_GENES), carrier_condition=cond,
        gwas_associations=carried, trait=tr, literature_url=litvar, note=note,
    )


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

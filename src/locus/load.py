"""Load phase: stream the (annotated) sites VCF into DuckDB.

Reads the annotated VCF if present, else the plain sites VCF, via cyvcf2 →
pyarrow batches → DuckDB (zero-copy Arrow). Also loads CNV/SV records and the
PharmCAT report when those exist. The result is a single read-only-queryable
``locus.duckdb`` consumed by the MCP server and the SPA.
"""

from __future__ import annotations

import datetime as _dt
import json
from pathlib import Path

import duckdb
import pyarrow as pa
from cyvcf2 import VCF
from rich.console import Console

from . import artifacts
from .config import settings
from .db import SCHEMA_VERSION

console = Console()

BATCH = 50_000

# INFO fields we promote to typed columns when the annotation step has added them.
_CLINVAR_FIELDS = ("CLNSIG", "CLNDN", "CLNREVSTAT", "CLNVC", "CLNDISDB", "CLNHGVS", "MC", "ALLELEID")
_GNOMAD_FIELDS = ("gnomAD_AF", "gnomAD_AF_grpmax", "gnomAD_grpmax", "gnomAD_AC", "gnomAD_AN")


def _info_ids(vcf: VCF) -> set[str]:
    ids = set()
    for line in vcf.raw_header.split("\n"):
        if line.startswith("##INFO=<ID="):
            ids.add(line.split("ID=", 1)[1].split(",", 1)[0])
    return ids


def _as_str(val) -> str | None:
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        return ",".join(str(v) for v in val)
    return str(val)


def _gene_from_info(v, info_ids: set[str]) -> str | None:
    """Best-effort gene symbol: ClinVar GENEINFO ('BRCA1:672|...') or VEP CSQ SYMBOL."""
    if "GENEINFO" in info_ids:
        gi = v.INFO.get("GENEINFO")
        if gi:
            return str(gi).split(":", 1)[0].split("|", 1)[0]
    if "ANN" in info_ids:  # SnpEff ANN: allele|effect|impact|gene|...
        ann = v.INFO.get("ANN")
        if ann:
            parts = str(ann).split("|")
            if len(parts) > 3 and parts[3]:
                return parts[3]
    return None


def _consequence_from_info(v, info_ids: set[str]) -> str | None:
    if "ANN" in info_ids:
        ann = v.INFO.get("ANN")
        if ann:
            parts = str(ann).split("|")
            if len(parts) > 1:
                return parts[1]
    return None


def _empty_cols() -> dict[str, list]:
    cols = {
        "chrom": [], "pos": [], "ref": [], "alt": [], "rsid": [],
        "qual": [], "filter": [], "gt": [], "gene": [], "consequence": [],
    }
    for f in _CLINVAR_FIELDS:
        cols[f.lower()] = []
    for f in _GNOMAD_FIELDS:
        cols[f.lower()] = []
    cols["am_pathogenicity"] = []
    cols["am_class"] = []
    return cols


def _variant_batches(path: Path):
    vcf = VCF(str(path))
    info_ids = _info_ids(vcf)
    cols = _empty_cols()
    n = 0
    for v in vcf:
        for alt in (v.ALT or [""]):
            cols["chrom"].append(v.CHROM)
            cols["pos"].append(int(v.POS))
            cols["ref"].append(v.REF)
            cols["alt"].append(alt)
            cols["rsid"].append(None if v.ID in (None, ".") else v.ID)
            cols["qual"].append(float(v.QUAL) if v.QUAL is not None else None)
            cols["filter"].append("PASS" if v.FILTER is None else v.FILTER)
            cols["gt"].append(_fmt_gt(v))
            cols["gene"].append(_gene_from_info(v, info_ids))
            cols["consequence"].append(_consequence_from_info(v, info_ids))
            for f in _CLINVAR_FIELDS:
                cols[f.lower()].append(_as_str(v.INFO.get(f)) if f in info_ids else None)
            for f in _GNOMAD_FIELDS:
                val = v.INFO.get(f) if f in info_ids else None
                cols[f.lower()].append(float(val) if isinstance(val, (int, float)) else _as_str(val))
            amp = v.INFO.get("am_pathogenicity") if "am_pathogenicity" in info_ids else None
            cols["am_pathogenicity"].append(float(amp) if isinstance(amp, (int, float)) else None)
            cols["am_class"].append(_as_str(v.INFO.get("am_class")) if "am_class" in info_ids else None)
            n += 1
            if len(cols["pos"]) >= BATCH:
                yield pa.table(cols)
                cols = _empty_cols()
    if cols["pos"]:
        yield pa.table(cols)
    vcf.close()


def _fmt_gt(v) -> str | None:
    try:
        a, b, *_ = v.genotypes[0]
        sep = "|" if v.genotypes[0][2] else "/"
        return f"{'.' if a < 0 else a}{sep}{'.' if b < 0 else b}"
    except (IndexError, TypeError):
        return None


def _create_schema(con: duckdb.DuckDBPyConnection) -> None:
    con.execute("""
        CREATE TABLE IF NOT EXISTS meta(
            key VARCHAR, value VARCHAR
        );
        CREATE TABLE IF NOT EXISTS variants(
            chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR,
            rsid VARCHAR, qual DOUBLE, filter VARCHAR, gt VARCHAR,
            gene VARCHAR, consequence VARCHAR,
            clnsig VARCHAR, clndn VARCHAR, clnrevstat VARCHAR, clnvc VARCHAR,
            clndisdb VARCHAR, clnhgvs VARCHAR, mc VARCHAR, alleleid VARCHAR,
            gnomad_af DOUBLE, gnomad_af_grpmax DOUBLE, gnomad_grpmax VARCHAR,
            gnomad_ac VARCHAR, gnomad_an VARCHAR,
            am_pathogenicity DOUBLE, am_class VARCHAR
        );
        CREATE TABLE IF NOT EXISTS cnv(
            chrom VARCHAR, pos BIGINT, "end" BIGINT, svtype VARCHAR, alt VARCHAR,
            cn INTEGER, svlen BIGINT, filter VARCHAR, genes VARCHAR
        );
        CREATE TABLE IF NOT EXISTS sv(
            chrom VARCHAR, pos BIGINT, "end" BIGINT, svtype VARCHAR, alt VARCHAR,
            svlen BIGINT, mateid VARCHAR, filter VARCHAR, pr VARCHAR, sr VARCHAR, genes VARCHAR
        );
        CREATE TABLE IF NOT EXISTS pgx_genes(
            gene VARCHAR, diplotype VARCHAR, phenotype VARCHAR, activity_score VARCHAR
        );
        CREATE TABLE IF NOT EXISTS pgx_drugs(
            drug VARCHAR, gene VARCHAR, source VARCHAR, recommendation VARCHAR, classification VARCHAR
        );
        CREATE TABLE IF NOT EXISTS ancestry_global(
            superpop VARCHAR, name VARCHAR, proportion DOUBLE
        );
        CREATE TABLE IF NOT EXISTS ancestry_pca(
            label VARCHAR, pc1 DOUBLE, pc2 DOUBLE, is_sample BOOLEAN
        );
        CREATE TABLE IF NOT EXISTS pgs_scores(
            pgs_id VARCHAR, trait VARCHAR, raw DOUBLE, percentile DOUBLE,
            ancestry VARCHAR, n_used INTEGER, coverage DOUBLE
        );
        CREATE TABLE IF NOT EXISTS ancestry_segments(
            haplotype INTEGER, chrom VARCHAR, start BIGINT, "end" BIGINT,
            ancestry VARCHAR, posterior DOUBLE
        );
    """)


def write_segments(segments: list) -> None:
    """Write local-ancestry segments (chromosome painting) into the DuckDB store."""
    import duckdb as _d

    from .config import settings

    con = _d.connect(str(settings.db_path))
    try:
        _create_schema(con)
        con.execute("DELETE FROM ancestry_segments")
        if segments:
            con.executemany("INSERT INTO ancestry_segments VALUES (?,?,?,?,?,?)", segments)
    finally:
        con.close()


def write_ancestry(ancestry_result, pgs_scores: list) -> None:
    """Write ancestry + polygenic-score results into the existing DuckDB store.

    Runs as a standalone step (``locus ancestry``) — replaces just these tables,
    leaving the variant tables intact. Requires no other process holding the DB.
    """
    import duckdb as _d

    from .config import settings

    con = _d.connect(str(settings.db_path))
    try:
        _create_schema(con)
        con.execute("DELETE FROM ancestry_global; DELETE FROM ancestry_pca; DELETE FROM pgs_scores")
        if ancestry_result is not None:
            from .ancestry import SUPERPOPS

            con.executemany("INSERT INTO ancestry_global VALUES (?,?,?)", [
                (sp, SUPERPOPS.get(sp, sp), frac)
                for sp, frac in ancestry_result.proportions.items() if frac > 0
            ])
            rows = [("you", ancestry_result.sample_pcs[0], ancestry_result.sample_pcs[1], True)]
            for sp, (p1, p2) in ancestry_result.ref_centroids.items():
                rows.append((sp, p1, p2, False))
            con.executemany("INSERT INTO ancestry_pca VALUES (?,?,?,?)", rows)
        if pgs_scores:
            con.executemany("INSERT INTO pgs_scores VALUES (?,?,?,?,?,?,?)", [
                (s.pgs_id, s.trait, s.raw, s.percentile, s.ancestry, s.n_used, s.coverage)
                for s in pgs_scores
            ])
    finally:
        con.close()


def _finalize(con: duckdb.DuckDBPyConnection) -> None:
    # Physically order by coordinate (zonemap pruning) and add point-lookup indexes.
    con.execute("CREATE TABLE variants_sorted AS SELECT * FROM variants ORDER BY chrom, pos")
    con.execute("DROP TABLE variants")
    con.execute("ALTER TABLE variants_sorted RENAME TO variants")
    con.execute("CREATE INDEX idx_var_pos ON variants(chrom, pos)")
    con.execute("CREATE INDEX idx_var_rsid ON variants(rsid)")
    con.execute("CREATE INDEX idx_var_gene ON variants(gene)")


def load_cnv(con: duckdb.DuckDBPyConnection, cnv_vcf: Path) -> int:
    vcf = VCF(str(cnv_vcf))
    rows = []
    for v in vcf:
        cn = v.format("CN")
        cn_val = int(cn[0][0]) if cn is not None and len(cn) else None
        rows.append((
            v.CHROM, int(v.POS), int(v.INFO.get("END") or v.end),
            _as_str(v.INFO.get("SVTYPE")), _as_str(v.ALT),
            cn_val, _coerce_int(v.INFO.get("SVLEN")),
            "PASS" if v.FILTER is None else v.FILTER, None,
        ))
    vcf.close()
    if rows:
        con.executemany("INSERT INTO cnv VALUES (?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def load_sv(con: duckdb.DuckDBPyConnection, sv_vcf: Path) -> int:
    vcf = VCF(str(sv_vcf))
    rows = []
    for v in vcf:
        pr = v.format("PR")
        sr = v.format("SR")
        rows.append((
            v.CHROM, int(v.POS), int(v.INFO.get("END") or v.end or v.POS),
            _as_str(v.INFO.get("SVTYPE")), _as_str(v.ALT),
            _coerce_int(v.INFO.get("SVLEN")), _as_str(v.INFO.get("MATEID")),
            "PASS" if v.FILTER is None else v.FILTER,
            _as_str(pr[0]) if pr is not None and len(pr) else None,
            _as_str(sr[0]) if sr is not None and len(sr) else None, None,
        ))
    vcf.close()
    if rows:
        con.executemany("INSERT INTO sv VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def _coerce_int(val):
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        val = val[0]
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _diplotype_str(dip: dict) -> tuple[str | None, str | None, str | None]:
    """From a PharmCAT diplotype object -> (diplotype, phenotype, activity_score)."""
    a1 = (dip.get("allele1") or {}).get("name")
    a2 = (dip.get("allele2") or {}).get("name")
    if a1 and a2:
        diplo = f"{a1}/{a2}"
    else:
        diplo = a1 or a2
    phenos = dip.get("phenotypes") or []
    return diplo, ("; ".join(phenos) if phenos else None), _as_str(dip.get("activityScore"))


def load_pharmcat(con: duckdb.DuckDBPyConnection, report_json: Path) -> tuple[int, int]:
    """Parse a PharmCAT report.json. genes -> {symbol: {...sourceDiplotypes}}, drugs -> {source: {drug: {...}}}."""
    data = json.loads(report_json.read_text())

    gene_rows = []
    for gene, payload in (data.get("genes") or {}).items():
        if not isinstance(payload, dict):
            continue
        dips = payload.get("sourceDiplotypes") or payload.get("recommendationDiplotypes") or []
        if dips:
            diplo, pheno, activity = _diplotype_str(dips[0])
        else:
            diplo = pheno = activity = None
        gene_rows.append((gene, diplo, pheno, activity))

    drug_rows = []
    seen: set = set()
    for source_name, drugs in (data.get("drugs") or {}).items():
        if not isinstance(drugs, dict):
            continue
        for drug, payload in drugs.items():
            if not isinstance(payload, dict):
                continue
            for gl in payload.get("guidelines") or []:
                for ann in gl.get("annotations") or []:
                    rec = ann.get("drugRecommendation")
                    if not rec:
                        continue
                    gene = None
                    gts = ann.get("genotypes") or []
                    if gts:
                        gdips = gts[0].get("diplotypes") or []
                        if gdips:
                            gene = (gdips[0].get("allele1") or {}).get("gene")
                    key = (drug, gene, rec)
                    if key in seen:
                        continue
                    seen.add(key)
                    drug_rows.append((drug, gene, gl.get("source") or source_name, rec, ann.get("classification")))

    if gene_rows:
        con.executemany("INSERT INTO pgx_genes VALUES (?,?,?,?)", gene_rows)
    if drug_rows:
        con.executemany("INSERT INTO pgx_drugs VALUES (?,?,?,?,?)", drug_rows)
    return len(gene_rows), len(drug_rows)


def run() -> Path:
    """Build the DuckDB store from the pipeline artifacts."""
    src = artifacts.annotated_vcf() if artifacts.annotated_vcf().exists() else artifacts.sites_vcf()
    if not src.exists():
        raise FileNotFoundError(f"No sites/annotated VCF found ({src}). Run `locus ingest` first.")

    db_path = settings.db_path
    db_path.parent.mkdir(parents=True, exist_ok=True)

    console.rule("[bold]Load → DuckDB")
    console.print(f"source : {src.name}")
    con = duckdb.connect(str(db_path))
    try:
        # Rebuild only the variant-derived tables; preserve ancestry/PGS (written by `locus ancestry`).
        for t in ("variants", "cnv", "sv", "pgx_genes", "pgx_drugs", "meta"):
            con.execute(f"DROP TABLE IF EXISTS {t}")
        _create_schema(con)
        total = 0
        for tbl in _variant_batches(src):
            con.register("batch_arrow", tbl)
            con.execute("INSERT INTO variants SELECT * FROM batch_arrow")
            con.unregister("batch_arrow")
            total += tbl.num_rows
        console.print(f"variants : {total:,}")

        inputs = artifacts.classify_inputs(settings.genome_dir)
        if inputs.cnv and inputs.cnv.exists():
            console.print(f"cnv      : {load_cnv(con, inputs.cnv):,}")
        if inputs.sv and inputs.sv.exists():
            console.print(f"sv       : {load_sv(con, inputs.sv):,}")

        pcat_dir = artifacts.pharmcat_dir()
        pharmcat_report = next(pcat_dir.glob("*.report.json"), None) if pcat_dir.exists() else None
        if pharmcat_report:
            g, d = load_pharmcat(con, pharmcat_report)
            console.print(f"pgx      : {g} genes, {d} drugs")

        _finalize(con)
        con.execute("DELETE FROM meta")
        con.executemany("INSERT INTO meta VALUES (?, ?)", [
            ("schema_version", str(SCHEMA_VERSION)),
            ("sample_id", settings.sample_id),
            ("source_vcf", src.name),
            ("created_at", _dt.datetime.now().isoformat(timespec="seconds")),
            ("variant_count", str(total)),
        ])
    finally:
        con.close()
    console.print(f"[green]DuckDB store ready[/] → {db_path}")
    return db_path

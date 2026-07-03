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
from .vcfutils import canonical_chrom

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
            level VARCHAR, code VARCHAR, name VARCHAR, proportion DOUBLE
        );
        CREATE TABLE IF NOT EXISTS ancestry_pca(
            label VARCHAR, pc1 DOUBLE, pc2 DOUBLE, is_sample BOOLEAN, "group" VARCHAR
        );
        CREATE TABLE IF NOT EXISTS pgs_scores(
            pgs_id VARCHAR, trait VARCHAR, raw DOUBLE, percentile DOUBLE,
            ancestry VARCHAR, n_used INTEGER, coverage DOUBLE
        );
        -- The "living" refresh spine: which source versions we last saw, IDs already
        -- processed, and the ranked "what's new since last run" findings. All three
        -- are PRESERVED across a variant reload (not in run()'s drop list).
        CREATE TABLE IF NOT EXISTS sources(
            name VARCHAR, version VARCHAR, url VARCHAR, checksum VARCHAR,
            license VARCHAR, last_checked VARCHAR, last_updated VARCHAR
        );
        CREATE TABLE IF NOT EXISTS watch_seen_ids(
            source VARCHAR, external_id VARCHAR
        );
        CREATE TABLE IF NOT EXISTS watch_findings(
            ts VARCHAR, source VARCHAR, kind VARCHAR, tier VARCHAR,
            chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR, rsid VARCHAR, gene VARCHAR,
            title VARCHAR, detail VARCHAR, old_value VARCHAR, new_value VARCHAR,
            release VARCHAR, url VARCHAR
        );
        -- Single-SNP traits / wellness + HLA proxy (preserved across a variant reload).
        CREATE TABLE IF NOT EXISTS traits(
            rsid VARCHAR, category VARCHAR, trait VARCHAR, genotype VARCHAR,
            dosage INTEGER, effect_allele VARCHAR, interpretation VARCHAR, note VARCHAR
        );
        -- GWAS Catalog risk alleles the genome carries (weak/exploratory; preserved).
        CREATE TABLE IF NOT EXISTS associations(
            rsid VARCHAR, chrom VARCHAR, pos BIGINT, risk_allele VARCHAR, dosage INTEGER,
            zygosity VARCHAR, trait VARCHAR, mapped_trait VARCHAR, pval DOUBLE,
            or_beta VARCHAR, pmid VARCHAR
        );
    """)


def write_ancestry(ancestry_result, pgs_scores: list) -> None:
    """Write ancestry + polygenic-score results into the existing DuckDB store.

    Runs as a standalone step (``locus ancestry``) — replaces just these tables,
    leaving the variant tables intact. Requires no other process holding the DB.
    """
    import duckdb as _d

    from .config import settings

    con = _d.connect(str(settings.db_path))
    try:
        # Drop+recreate (tolerates schema changes), then repopulate.
        for t in ("ancestry_global", "ancestry_pca", "pgs_scores"):
            con.execute(f"DROP TABLE IF EXISTS {t}")
        _create_schema(con)
        if ancestry_result is not None:
            from .ancestry import POPULATIONS, SUPERPOPS

            rows = [("continental", sp, SUPERPOPS.get(sp, sp), frac)
                    for sp, frac in ancestry_result.proportions.items() if frac > 0]
            rows += [("population", pop, POPULATIONS.get(pop, pop), frac)
                     for pop, frac in ancestry_result.populations.items() if frac > 0]
            con.executemany("INSERT INTO ancestry_global VALUES (?,?,?,?)", rows)
            # PCA scatter: the sample + the fine-population reference cloud (colored by
            # continent via `group`). Falls back to superpop centroids if pop_points absent.
            rows = [("you", ancestry_result.sample_pcs[0], ancestry_result.sample_pcs[1],
                     True, ancestry_result.nearest)]
            pop_points = getattr(ancestry_result, "pop_points", None)
            if pop_points:
                rows += [(label, p1, p2, False, sp) for label, p1, p2, sp in pop_points]
            else:
                rows += [(sp, p1, p2, False, sp) for sp, (p1, p2) in ancestry_result.ref_centroids.items()]
            con.executemany("INSERT INTO ancestry_pca VALUES (?,?,?,?,?)", rows)
        if pgs_scores:
            con.executemany("INSERT INTO pgs_scores VALUES (?,?,?,?,?,?,?)", [
                (s.pgs_id, s.trait, s.raw, s.percentile, s.ancestry, s.n_used, s.coverage)
                for s in pgs_scores
            ])
    finally:
        con.close()


def write_traits(results: list) -> None:
    """Replace the ``traits`` table with freshly-computed tag-SNP results.

    Standalone step (``locus traits``) — like write_ancestry, it touches only its own
    table and leaves the variant tables intact.
    """
    import duckdb as _d

    con = _d.connect(str(settings.db_path))
    try:
        con.execute("DROP TABLE IF EXISTS traits")
        _create_schema(con)
        if results:
            con.executemany(
                "INSERT INTO traits VALUES (?,?,?,?,?,?,?,?)",
                [(r.rsid, r.category, r.trait, r.genotype, r.dosage, r.effect_allele,
                  r.interpretation, r.note) for r in results],
            )
    finally:
        con.close()


def write_associations(carried: list) -> None:
    """Replace the ``associations`` table with carried GWAS risk alleles. Standalone step
    (``locus gwas``); leaves variant tables intact."""
    import duckdb as _d

    con = _d.connect(str(settings.db_path))
    try:
        con.execute("DROP TABLE IF EXISTS associations")
        _create_schema(con)
        if carried:
            con.executemany(
                "INSERT INTO associations VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                [(c.rsid, c.chrom, c.pos, c.risk_allele, c.dosage, c.zygosity, c.trait,
                  c.mapped_trait, c.pval, c.or_beta, c.pmid) for c in carried],
            )
    finally:
        con.close()


def upsert_source(con: duckdb.DuckDBPyConnection, name: str, *, version: str, url: str = "",
                  checksum: str = "", license: str = "", last_checked: str = "",
                  last_updated: str = "") -> None:
    """Record/replace the last-seen state of an external source (projection of the manifest)."""
    con.execute("CREATE TABLE IF NOT EXISTS sources(name VARCHAR, version VARCHAR, url VARCHAR, "
                "checksum VARCHAR, license VARCHAR, last_checked VARCHAR, last_updated VARCHAR)")
    con.execute("DELETE FROM sources WHERE name = ?", [name])
    con.execute("INSERT INTO sources VALUES (?,?,?,?,?,?,?)",
                [name, version, url, checksum, license, last_checked, last_updated])


def append_findings(con: duckdb.DuckDBPyConnection, rows: list[tuple]) -> int:
    """Append ``watch_findings`` rows. Each row matches the watch_findings column order:
    (ts, source, kind, tier, chrom, pos, ref, alt, rsid, gene, title, detail,
    old_value, new_value, release, url)."""
    if not rows:
        return 0
    con.execute("CREATE TABLE IF NOT EXISTS watch_findings(ts VARCHAR, source VARCHAR, kind VARCHAR, "
                "tier VARCHAR, chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR, rsid VARCHAR, "
                "gene VARCHAR, title VARCHAR, detail VARCHAR, old_value VARCHAR, new_value VARCHAR, "
                "release VARCHAR, url VARCHAR)")
    # Migrate a pre-v4 store in place so old changelog history is preserved.
    con.execute("ALTER TABLE watch_findings ADD COLUMN IF NOT EXISTS url VARCHAR")
    con.executemany("INSERT INTO watch_findings VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    return len(rows)


def mark_seen(con: duckdb.DuckDBPyConnection, source: str, ids: list[str]) -> None:
    """Record external IDs already processed for a source (e.g. PGS score IDs)."""
    if not ids:
        return
    con.execute("CREATE TABLE IF NOT EXISTS watch_seen_ids(source VARCHAR, external_id VARCHAR)")
    con.executemany("INSERT INTO watch_seen_ids VALUES (?, ?)", [(source, i) for i in ids])


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
            canonical_chrom(v.CHROM), int(v.POS), int(v.INFO.get("END") or v.end),
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
            canonical_chrom(v.CHROM), int(v.POS), int(v.INFO.get("END") or v.end or v.POS),
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

"""The "living genome" refresh engine.

The genome is sequenced once and static; the *interpretation* databases move. This
module polls each tracked source, detects new releases against ``data/manifest.json``,
re-interprets only what changed, and writes a ranked "what's new since last run"
changelog into the DuckDB ``watch_findings`` table (surfaced by ``queries.whats_new``,
the ``whats_new`` MCP tool, and the SPA Changelog tab).

Headline capability: **ClinVar reanalysis** — when a new ClinVar release lands, diff
its classifications against the ones the genome currently carries and surface variants
that became (or stopped being) pathogenic. The diff is naturally restricted to carried
positions because it compares the genome's own ``variants.clnsig`` before vs after a
re-annotate + reload — reusing the existing pipeline rather than re-implementing it.

Privacy: every outbound call sends only generic queries (release dates, public md5s,
score IDs). The genome never leaves the machine; all matching is local.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rich.console import Console

from . import annotate, download, load, manifest
from .config import settings
from .db import connect

console = Console()

# Re-annotate the full local-DB chain on a ClinVar refresh so gnomAD AF + SnpEff genes +
# AlphaMissense all survive. annotate.run() always rebuilds from the sites VCF, so any step left
# out here is silently *dropped* from the store — that's the bug that quietly removed gnomAD AF.
# Keep every local-DB step listed. (gnomAD is affordable here now that AF comes from Ensembl in
# seconds rather than streaming hundreds of GB of gnomAD VCFs.)
_REANNOTATE_STEPS = "clinvar,snpeff,alphamissense,gnomad"

# Tiers, strongest first — controls ranking so signal isn't buried under noise.
TIER_ORDER = ["strong", "moderate", "weak", "info"]


@dataclass
class Finding:
    source: str
    kind: str            # e.g. "newly_pathogenic", "reclassified", "withdrawn", "release"
    tier: str            # strong | moderate | weak | info
    title: str
    detail: str = ""
    chrom: str | None = None
    pos: int | None = None
    ref: str | None = None
    alt: str | None = None
    rsid: str | None = None
    gene: str | None = None
    old_value: str | None = None
    new_value: str | None = None
    release: str | None = None
    url: str | None = None       # citation / source link (PubMed, GWAS Catalog)

    def row(self, ts: str) -> tuple:
        return (ts, self.source, self.kind, self.tier, self.chrom, self.pos, self.ref,
                self.alt, self.rsid, self.gene, self.title, self.detail,
                self.old_value, self.new_value, self.release, self.url)


# ── HTTP probes (small, generic queries only) ──────────────────────────────────
def _http_text(url: str, timeout: int = 30) -> str | None:
    import httpx
    try:
        r = httpx.get(url, timeout=timeout, follow_redirects=True)
        r.raise_for_status()
        return r.text
    except Exception as e:  # noqa: BLE001 - network is best-effort; report & continue
        console.print(f"[yellow]probe failed[/] {url}: {e}")
        return None


def _http_json(url: str, timeout: int = 30):
    import httpx
    try:
        r = httpx.get(url, timeout=timeout, follow_redirects=True)
        r.raise_for_status()
        return r.json()
    except Exception as e:  # noqa: BLE001
        console.print(f"[yellow]probe failed[/] {url}: {e}")
        return None


def _http_stamp(url: str, timeout: int = 30) -> str | None:
    """A cheap change-detector for a bulk file: its Last-Modified (or ETag) header."""
    import httpx
    try:
        r = httpx.head(url, timeout=timeout, follow_redirects=True)
        r.raise_for_status()
        return r.headers.get("last-modified") or r.headers.get("etag")
    except Exception as e:  # noqa: BLE001
        console.print(f"[yellow]probe failed[/] {url}: {e}")
        return None


def _days_ago_iso(days: int) -> str:
    import datetime as _dt
    return (_dt.datetime.now() - _dt.timedelta(days=days)).date().isoformat()


def _seen_ids(source: str) -> set[str]:
    """External IDs already processed for a source (dedup, e.g. PubMed PMIDs)."""
    with connect(read_only=True) as con:
        try:
            return {i for (i,) in con.execute(
                "SELECT external_id FROM watch_seen_ids WHERE source = ?", [source]).fetchall()}
        except Exception:  # noqa: BLE001 - table may be absent
            return set()


# ── ClinVar classification helpers (pure; unit-testable) ────────────────────────
def is_pathogenic(clnsig: str | None) -> bool:
    """Pathogenic / likely-pathogenic, excluding benign and conflicting (matches
    queries.clinical_findings)."""
    if not clnsig:
        return False
    s = clnsig.lower()
    return "pathogenic" in s and "benign" not in s and "conflicting" not in s


def review_strength(clnrevstat: str | None) -> str:
    """Map a ClinVar review status (star rating) to a confidence band."""
    s = (clnrevstat or "").lower()
    if "practice_guideline" in s or "reviewed_by_expert_panel" in s or "multiple_submitters" in s:
        return "strong"   # >=2 stars
    if "single_submitter" in s:
        return "moderate"  # 1 star
    return "weak"          # 0 stars / no assertion


def classify_clinvar_delta(prev: dict, cur: dict) -> list[Finding]:
    """Diff two snapshots of the genome's ClinVar state.

    ``prev``/``cur`` map (chrom, pos, ref, alt) -> dict(clnsig, clnrevstat, gene, rsid, clndn).
    Surfaces variants that *became* pathogenic (the actionable case, tiered by review
    status), variants reclassified to/from pathogenic, and pathogenic classifications
    that were withdrawn. Pure function — no I/O.
    """
    findings: list[Finding] = []
    for key, c in cur.items():
        p = prev.get(key)
        old_sig = (p or {}).get("clnsig")
        new_sig = c.get("clnsig")
        if old_sig == new_sig:
            continue
        chrom, pos, ref, alt = key
        common = dict(chrom=chrom, pos=pos, ref=ref, alt=alt,
                      rsid=c.get("rsid"), gene=c.get("gene"),
                      old_value=old_sig, new_value=new_sig)
        now_path, was_path = is_pathogenic(new_sig), is_pathogenic(old_sig)
        if now_path and not was_path:
            band = review_strength(c.get("clnrevstat"))
            findings.append(Finding(
                source="clinvar", kind="newly_pathogenic", tier=band,
                title=f"{c.get('gene') or chrom}:{pos} now {new_sig}",
                detail=(c.get("clndn") or "").replace("_", " ")[:300], **common))
        elif was_path and not now_path:
            findings.append(Finding(
                source="clinvar", kind="depathogenized", tier="moderate",
                title=f"{(p or {}).get('gene') or chrom}:{pos} no longer pathogenic ({old_sig}→{new_sig})",
                detail="ClinVar reclassified away from pathogenic — re-check.", **common))
        else:
            findings.append(Finding(
                source="clinvar", kind="reclassified", tier="weak",
                title=f"{c.get('gene') or chrom}:{pos} reclassified ({old_sig}→{new_sig})",
                detail=(c.get("clndn") or "").replace("_", " ")[:300], **common))
    # Pathogenic positions that disappeared entirely from ClinVar (withdrawn record).
    for key, p in prev.items():
        if key not in cur and is_pathogenic(p.get("clnsig")):
            chrom, pos, ref, alt = key
            findings.append(Finding(
                source="clinvar", kind="withdrawn", tier="moderate",
                title=f"{p.get('gene') or chrom}:{pos} pathogenic record withdrawn",
                detail="A previously-pathogenic ClinVar record is gone — re-check.",
                chrom=chrom, pos=pos, ref=ref, alt=alt, rsid=p.get("rsid"),
                gene=p.get("gene"), old_value=p.get("clnsig"), new_value=None))
    return findings


def _snapshot_clinvar() -> dict:
    """Current annotated-ClinVar state of the genome, keyed by (chrom,pos,ref,alt)."""
    snap: dict = {}
    with connect(read_only=True) as con:
        rows = con.execute(
            "SELECT chrom, pos, ref, alt, clnsig, clnrevstat, gene, rsid, clndn "
            "FROM variants WHERE clnsig IS NOT NULL"
        ).fetchall()
    for chrom, pos, ref, alt, clnsig, clnrevstat, gene, rsid, clndn in rows:
        snap[(chrom, pos, ref, alt)] = dict(
            clnsig=clnsig, clnrevstat=clnrevstat, gene=gene, rsid=rsid, clndn=clndn)
    return snap


# ── Source checkers ─────────────────────────────────────────────────────────────
def check_clinvar(man: dict) -> dict | None:
    """Probe ClinVar's rolling md5. Returns {version, url, changed} or None on failure."""
    txt = _http_text(f"{download.CLINVAR_BASE}/clinvar.vcf.gz.md5")
    if not txt:
        return None
    remote_md5 = txt.split()[0] if txt.split() else ""
    if not remote_md5:
        return None
    prev = manifest.get(man, "clinvar")
    current = prev.get("checksum")
    if not current:
        # Bootstrap from the locally-stored md5 (written by download_clinvar) if present.
        local = settings.annotations_dir / "clinvar.vcf.gz.md5"
        if local.exists():
            parts = local.read_text().split()
            current = parts[0] if parts else None
    return {"version": remote_md5, "checksum": remote_md5,
            "url": f"{download.CLINVAR_BASE}/clinvar.vcf.gz", "changed": remote_md5 != current}


def check_cpic(man: dict) -> dict | None:
    """Probe CPIC guideline versions. Returns {version, versions, names, prev_versions, changed}."""
    import hashlib

    data = _http_json("https://api.cpicpgx.org/v1/guideline?select=id,name,version")
    if not isinstance(data, list) or not data:
        return None
    versions = {str(g["id"]): g.get("version") for g in data if "id" in g}
    names = {str(g["id"]): g.get("name", "") for g in data}
    sig = ";".join(f"{k}:{versions[k]}" for k in sorted(versions))
    version = hashlib.md5(sig.encode()).hexdigest()[:12]
    prev = manifest.get(man, "cpic")
    prev_versions = prev.get("versions") if isinstance(prev.get("versions"), dict) else {}
    return {"version": version, "versions": versions, "names": names,
            "prev_versions": prev_versions, "url": "https://api.cpicpgx.org/v1/guideline",
            "changed": version != prev.get("version")}


def _cpic_findings(check: dict) -> list[Finding]:
    """Surface CPIC guideline version bumps that touch a gene the genome has a PGx call for."""
    prev, cur, names = check["prev_versions"], check["versions"], check["names"]
    if not prev:  # first run — just record the baseline, don't flood
        return []
    with connect(read_only=True) as con:
        try:
            genes = {g for (g,) in con.execute("SELECT DISTINCT gene FROM pgx_genes").fetchall() if g}
        except Exception:  # noqa: BLE001 - table may be absent
            genes = set()
    out: list[Finding] = []
    for gid, ver in cur.items():
        if prev.get(gid) == ver:
            continue
        name = names.get(gid, "")
        if any(g and g in name for g in genes):
            out.append(Finding(
                source="cpic", kind="guideline_update", tier="moderate",
                title=f"CPIC guideline updated: {name}",
                detail="A pharmacogenomic guideline for a gene you carry a result for changed — "
                       "re-check the dosing guidance via `pharmacogenomics`.",
                new_value=str(ver), release=str(ver)))
    return out


def check_pgs(man: dict) -> dict | None:
    """Probe the PGS Catalog current release. Returns {version, n_new, ids, changed} or None."""
    data = _http_json("https://www.pgscatalog.org/rest/release/current")
    if not isinstance(data, dict) or "date" not in data:
        return None
    date = data.get("date")
    ids = data.get("released_score_ids") or []
    prev = manifest.get(man, "pgs_catalog")
    return {"version": date, "url": "https://www.pgscatalog.org/rest/release/current",
            "n_new": len(ids), "ids": ids, "changed": date != prev.get("version")}


# ── Per-source work ──────────────────────────────────────────────────────────────
def _reanalyze_clinvar() -> list[Finding]:
    """Snapshot → fetch the new ClinVar → re-annotate + reload → diff."""
    console.print("[bold]ClinVar reanalysis:[/] snapshotting current classifications…")
    prev = _snapshot_clinvar()

    # Force a fresh ClinVar download by removing the cached outputs, then re-prepare.
    ann = settings.annotations_dir
    for f in (download.CLINVAR_CHR_VCF, download.CLINVAR_CHR_VCF + ".tbi",
              "clinvar.vcf.gz.md5"):
        (ann / f).unlink(missing_ok=True)
    download.download_clinvar()

    console.print(f"Re-annotating ({_REANNOTATE_STEPS}) + reloading against the new ClinVar…")
    annotate.run(steps=_REANNOTATE_STEPS)
    load.run()

    cur = _snapshot_clinvar()
    findings = classify_clinvar_delta(prev, cur)
    console.print(f"ClinVar reanalysis: {len(findings)} change(s) at carried positions "
                  f"(prev {len(prev):,} → now {len(cur):,} annotated).")
    return findings


def _summarize_pgs(check: dict) -> list[Finding]:
    """Suggest-only PGS watcher: report new scores; never auto-add (ancestry-mismatch +
    runtime risk). The user adds relevant IDs to pgs.CURATED_PGS by hand."""
    n = check.get("n_new", 0)
    if not n:
        return []
    ids = check.get("ids", [])
    sample = ", ".join(ids[:8]) + (" …" if len(ids) > 8 else "")
    return [Finding(source="pgs_catalog", kind="release", tier="info",
                    title=f"{n} new polygenic score(s) published — review to track",
                    detail=f"Suggest-only: not auto-added. Add relevant IDs to pgs.CURATED_PGS, "
                           f"then `locus ancestry`. New: {sample}",
                    release=check.get("version"))]


# ── GWAS Catalog re-analysis (fully local; mirrors ClinVar reanalysis) ───────────
def _snapshot_gwas() -> dict:
    """Carried GWAS associations keyed by (rsid, mapped_trait)."""
    snap: dict = {}
    with connect(read_only=True) as con:
        try:
            rows = con.execute(
                "SELECT rsid, mapped_trait, dosage, or_beta, pval, pmid, chrom, pos FROM associations"
            ).fetchall()
        except Exception:  # noqa: BLE001 - table absent until `locus gwas`/first re-analysis
            rows = []
    for rsid, mapped_trait, dosage, or_beta, pval, pmid, chrom, pos in rows:
        snap[(rsid, mapped_trait)] = dict(dosage=dosage, or_beta=or_beta, pval=pval,
                                          pmid=pmid, chrom=chrom, pos=pos)
    return snap


def classify_gwas_delta(prev: dict, cur: dict) -> list[Finding]:
    """Newly-carried genome-wide-significant associations (keys in ``cur`` not in ``prev``).

    Pure function — no I/O. Each is a WEAK single hit (surfaced as such), never to be summed.
    """
    findings: list[Finding] = []
    for key, c in cur.items():
        if key in prev:
            continue
        rsid, trait = key
        pmid = c.get("pmid")
        findings.append(Finding(
            source="gwas", kind="new_association", tier="weak",
            title=f"{rsid} ↔ {trait} — you carry {c.get('dosage')} risk allele(s)",
            detail="New genome-wide-significant GWAS association at a variant you carry. WEAK single "
                   "hit — do not sum with others or read it like a calibrated risk score.",
            rsid=rsid, chrom=c.get("chrom"), pos=c.get("pos"),
            new_value=c.get("or_beta"), release=pmid,
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None))
    return findings


def check_gwas(man: dict) -> dict | None:
    """Probe the GWAS Catalog bulk file's Last-Modified/ETag. Returns {version, url, changed}."""
    stamp = _http_stamp(download.GWAS_ASSOC_ZIP)
    if not stamp:
        return None
    prev = manifest.get(man, "gwas_catalog")
    return {"version": stamp, "url": download.GWAS_ASSOC_ZIP, "changed": stamp != prev.get("version")}


def _reanalyze_gwas() -> list[Finding]:
    """Snapshot carried associations → re-fetch the catalog → recompute → diff (new carried hits)."""
    from . import gwas

    console.print("[bold]GWAS re-analysis:[/] snapshotting carried associations…")
    prev = _snapshot_gwas()
    # Force a fresh catalog download (setup_gwas skips when the TSV is present).
    (settings.annotations_dir / "gwas" / "gwas-catalog-associations.tsv").unlink(missing_ok=True)
    tsv = download.setup_gwas()
    console.print("Re-parsing GWAS Catalog (p < 5e-8 lead SNPs) + re-genotyping carried alleles…")
    carried = gwas.compute(gwas.parse(tsv))
    load.write_associations(carried)
    cur = _snapshot_gwas()
    if not prev:  # first time associations exist — record the baseline, don't flood
        console.print(f"GWAS baseline recorded ({len(cur):,} carried associations).")
        return []
    findings = classify_gwas_delta(prev, cur)
    console.print(f"GWAS re-analysis: {len(findings)} newly-carried association(s) "
                  f"(prev {len(prev):,} → now {len(cur):,}).")
    return findings


# ── PubMed literature watcher (sends only gene symbols; genotypes never leave) ────
def _watch_genes(cap: int = 200) -> list[str]:
    """The genome's notable genes for the gene-level PubMed sweep, most-actionable first:
    ClinVar P/LP genes, then PharmCAT PGx genes, then any gene with an AlphaMissense-predicted-
    damaging variant. Priority-ordered so the cap (a safety bound) drops the weakest signal first."""
    with connect(read_only=True) as con:
        def q(sql: str) -> list[str]:
            try:
                return [g for (g,) in con.execute(sql).fetchall() if g]
            except Exception:  # noqa: BLE001
                return []
        plp = q("SELECT DISTINCT gene FROM variants WHERE gene IS NOT NULL "
                "AND lower(clnsig) LIKE '%pathogenic%' AND lower(clnsig) NOT LIKE '%benign%' "
                "AND lower(clnsig) NOT LIKE '%conflicting%'")
        pgx = q("SELECT DISTINCT gene FROM pgx_genes")
        amd = q("SELECT DISTINCT gene FROM variants WHERE am_class LIKE '%pathogenic%'")
    ordered: list[str] = []
    seen: set[str] = set()
    for group in (sorted(set(plp)), sorted(set(pgx)), sorted(set(amd))):
        for g in group:
            if g not in seen:
                seen.add(g)
                ordered.append(g)
    return ordered[:cap]


def check_pubmed(man: dict) -> dict | None:
    """PubMed has no 'release' — we always look for papers published since we last checked."""
    prev = manifest.get(man, "pubmed")
    return {"version": manifest.now_iso(), "since": prev.get("last_checked"),
            "first_run": not prev, "changed": True, "url": "https://pubmed.ncbi.nlm.nih.gov/"}


def _pubmed_findings(check: dict) -> list[Finding]:
    """Search PubMed for recent papers on the genome's notable genes; dedup against seen PMIDs.

    On the first run we seed the seen-set from a bounded recent window but emit nothing (baseline),
    so later runs surface only genuinely new papers. Mutates ``check['new_pmids']`` for mark_seen.
    """
    from . import literature

    genes = _watch_genes()
    check["new_pmids"] = []
    if not genes:
        return []
    since = check.get("since") or _days_ago_iso(60)
    seen = _seen_ids("pubmed")
    first = check.get("first_run")
    findings: list[Finding] = []
    new_ids: list[str] = []
    for gene in genes:
        for h in literature.pubmed_search(
                f"{gene}[Gene] AND (variant OR mutation OR polymorphism OR pathogenic OR clinical)",
                mindate=since, retmax=5, gene=gene):
            if h.pmid in seen or h.pmid in new_ids:
                continue
            new_ids.append(h.pmid)
            if first:
                continue  # baseline only — don't flood with a backlog on first run
            snippet = (h.abstract[:300] + "…") if len(h.abstract) > 300 else h.abstract
            findings.append(Finding(
                source="pubmed", kind="new_study", tier="info", gene=gene,
                title=f"New paper on {gene}: {h.title[:160]}",
                detail=snippet or f"{h.journal} {h.year}".strip(),
                release=h.year or None, url=h.url))
    check["new_pmids"] = new_ids
    return findings


# ── LitVar2 variant watcher (papers about the specific variants you carry) ────────
def _litvar_variants(cap: int = 800) -> list[tuple[str, str | None]]:
    """The clinically-notable rsIDs to watch for new papers: anything ClinVar hasn't called plainly
    benign (P/LP, VUS/conflicting, risk/drug-response/protective/association), most-actionable first.
    Variant-centric — far higher signal than gene-level free-text, and it scales like ClinVar/GWAS."""
    with connect(read_only=True) as con:
        try:
            rows = con.execute(
                "SELECT DISTINCT rsid, gene, clnsig FROM variants "
                "WHERE rsid IS NOT NULL AND clnsig IS NOT NULL AND lower(clnsig) NOT LIKE '%benign%' "
                "AND (lower(clnsig) LIKE '%pathogenic%' OR lower(clnsig) LIKE '%uncertain%' "
                "  OR lower(clnsig) LIKE '%conflicting%' OR lower(clnsig) LIKE '%risk%' "
                "  OR lower(clnsig) LIKE '%drug_response%' OR lower(clnsig) LIKE '%protective%' "
                "  OR lower(clnsig) LIKE '%association%')"
            ).fetchall()
        except Exception:  # noqa: BLE001
            rows = []

    def rank(clnsig: str | None) -> int:
        s = (clnsig or "").lower()
        if "pathogenic" in s and "conflicting" not in s:
            return 0
        if "conflicting" in s or "uncertain" in s:
            return 1
        return 2  # risk_factor / drug_response / protective / association

    rows.sort(key=lambda r: rank(r[2]))
    out: list[tuple[str, str | None]] = []
    seen: set[str] = set()
    for rsid, gene, _ in rows:
        if rsid not in seen:
            seen.add(rsid)
            out.append((rsid, gene))
    return out[:cap]


def check_litvar(man: dict) -> dict:
    """LitVar has no 'release' — novelty is handled by the seen-PMID dedup, like PubMed."""
    prev = manifest.get(man, "litvar")
    return {"version": manifest.now_iso(), "first_run": not prev, "changed": True,
            "url": "https://www.ncbi.nlm.nih.gov/research/litvar2/"}


def _litvar_findings(check: dict) -> list[Finding]:
    """New papers about the specific variants you carry (LitVar2 rsID→literature), deduped by PMID.

    First run seeds the baseline silently. Mutates ``check['new_pmids']`` for mark_seen."""
    from . import literature

    variants = _litvar_variants()
    check["new_pmids"] = []
    if not variants:
        return []
    seen = _seen_ids("litvar")
    first = check.get("first_run")
    new: dict[str, tuple[str, str | None]] = {}  # pmid -> (rsid, gene)
    for rsid, gene in variants:
        for pmid in literature.litvar_pmids(rsid):
            if pmid in seen or pmid in new:
                continue
            new[pmid] = (rsid, gene)
    check["new_pmids"] = list(new)
    if first or not new:
        return []  # baseline: record the seen-set, emit nothing
    meta = literature.fetch_pubmed_meta(list(new))
    findings: list[Finding] = []
    for pmid, (rsid, gene) in new.items():
        h = meta.get(pmid)
        title = (h.title if h else "") or "(title unavailable)"
        snippet = ((h.abstract[:300] + "…") if h and len(h.abstract) > 300 else (h.abstract if h else "")) \
            or f"New paper mentioning {rsid}."
        label = rsid + (f" ({gene})" if gene else "")
        findings.append(Finding(
            source="litvar", kind="new_study", tier="info", gene=gene, rsid=rsid,
            title=f"New paper on {label}: {title[:150]}", detail=snippet,
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"))
    return findings


# ── Orchestrator ─────────────────────────────────────────────────────────────────
ALL_SOURCES = ["clinvar", "pgs", "cpic", "gwas", "pubmed", "litvar"]


def _write_report(findings: list[Finding], ts: str) -> Path | None:
    if not findings:
        return None
    settings.reports_dir.mkdir(parents=True, exist_ok=True)
    out = settings.reports_dir / "whats_new.md"
    lines = [f"# Locus — what's new ({ts})", ""]
    for tier in TIER_ORDER:
        tf = [f for f in findings if f.tier == tier]
        if not tf:
            continue
        lines.append(f"## {tier.title()} ({len(tf)})")
        for f in tf:
            loc = f" `{f.chrom}:{f.pos}`" if f.chrom else ""
            lines.append(f"- **{f.title}**{loc} — {f.detail}".rstrip(" —"))
        lines.append("")
    out.write_text("\n".join(lines))
    return out


def run(sources: str = "all", *, dry_run: bool = False, force: bool = False) -> list[Finding]:
    """Check each source, re-interpret what changed, and persist ranked findings.

    sources: comma-separated subset of {clinvar, pgs, cpic, gwas, pubmed, litvar} or "all".
    dry_run: probe + report what *would* change; write nothing.
    force:   run the per-source work even if the version looks unchanged.
    """
    requested = ALL_SOURCES if sources in ("all", "") else [s.strip() for s in sources.split(",")]
    man = manifest.load()
    ts = manifest.now_iso()
    console.rule("[bold]locus refresh" + (" (dry-run)" if dry_run else ""))

    findings: list[Finding] = []
    source_updates: list[tuple] = []  # (name, version, url, checksum, changed)
    pgs_new_ids: list[str] = []       # new PGS score IDs to mark seen
    pubmed_new_ids: list[str] = []    # PubMed PMIDs to mark seen (dedup)
    litvar_new_ids: list[str] = []    # LitVar PMIDs to mark seen (dedup)
    cpic_versions: dict | None = None  # per-guideline versions to persist

    if "clinvar" in requested:
        chk = check_clinvar(man)
        if chk is None:
            console.print("[yellow]ClinVar: probe failed — skipping.[/]")
        else:
            changed = chk["changed"] or force
            state = "NEW" if chk["changed"] else ("unchanged, forced" if force else "unchanged")
            console.print(f"ClinVar: remote {chk['version'][:12]}… ({state})")
            if changed and not dry_run:
                findings += _reanalyze_clinvar()
            elif changed and dry_run:
                findings.append(Finding(
                    source="clinvar", kind="update_available", tier="info",
                    title="ClinVar update available — run `locus refresh` to reanalyze",
                    new_value=chk["version"][:12]))
            source_updates.append(("clinvar", chk["version"], chk["url"], chk["checksum"], chk["changed"]))

    if "pgs" in requested:
        chk = check_pgs(man)
        if chk is None:
            console.print("[yellow]PGS Catalog: probe failed — skipping.[/]")
        else:
            console.print(f"PGS Catalog: release {chk['version']} "
                          f"({'NEW, ' + str(chk['n_new']) + ' scores' if chk['changed'] else 'unchanged'})")
            if chk["changed"]:
                findings += _summarize_pgs(chk)  # pure/cheap — fine to show in dry-run too
                pgs_new_ids = chk["ids"]
            source_updates.append(("pgs_catalog", chk["version"], chk["url"], "", chk["changed"]))

    if "cpic" in requested:
        chk = check_cpic(man)
        if chk is None:
            console.print("[yellow]CPIC: probe failed — skipping.[/]")
        else:
            console.print(f"CPIC: {len(chk['versions'])} guidelines "
                          f"({'NEW' if chk['changed'] else 'unchanged'})")
            if chk["changed"]:
                findings += _cpic_findings(chk)
            cpic_versions = chk["versions"]
            source_updates.append(("cpic", chk["version"], chk["url"], "", chk["changed"]))

    if "gwas" in requested:
        chk = check_gwas(man)
        if chk is None:
            console.print("[yellow]GWAS Catalog: probe failed — skipping.[/]")
        else:
            changed = chk["changed"] or force
            state = "NEW" if chk["changed"] else ("unchanged, forced" if force else "unchanged")
            console.print(f"GWAS Catalog: {str(chk['version'])[:32]} ({state})")
            if changed and not dry_run:
                findings += _reanalyze_gwas()
            elif changed and dry_run:
                findings.append(Finding(
                    source="gwas", kind="update_available", tier="info",
                    title="GWAS Catalog update available — run `locus refresh` to re-scan carried variants",
                    new_value=str(chk["version"])[:32]))
            source_updates.append(("gwas_catalog", chk["version"], chk["url"], "", chk["changed"]))

    if "pubmed" in requested:
        chk = check_pubmed(man)
        if dry_run:
            console.print("PubMed: [dim]dry-run — would search recent papers for your notable genes.[/]")
        else:
            pf = _pubmed_findings(chk)
            findings += pf
            pubmed_new_ids = chk.get("new_pmids", [])
            tag = " (baseline seeded)" if chk.get("first_run") and pubmed_new_ids else ""
            console.print(f"PubMed: {len(pf)} new paper finding(s){tag}")
        source_updates.append(("pubmed", chk["version"], chk["url"], "", True))

    if "litvar" in requested:
        chk = check_litvar(man)
        if dry_run:
            console.print("LitVar: [dim]dry-run — would scan the literature for the variants you carry.[/]")
        else:
            lf = _litvar_findings(chk)
            findings += lf
            litvar_new_ids = chk.get("new_pmids", [])
            tag = " (baseline seeded)" if chk.get("first_run") and litvar_new_ids else ""
            console.print(f"LitVar: {len(lf)} new variant-paper finding(s){tag}")
        source_updates.append(("litvar", chk["version"], chk["url"], "", True))

    # Rank strongest-first for the printed digest.
    findings.sort(key=lambda f: TIER_ORDER.index(f.tier) if f.tier in TIER_ORDER else len(TIER_ORDER))
    _print_digest(findings, dry_run)

    if dry_run:
        console.print("[dim]dry-run: no changes written.[/]")
        return findings

    # Persist: manifest (source of truth) + sources table + findings + report.
    for name, version, url, checksum, changed in source_updates:
        manifest.record(man, name, version=version, url=url, checksum=checksum, changed=changed)
    if cpic_versions is not None:  # persist the per-guideline version map for precise diffs
        man.setdefault("cpic", {})["versions"] = cpic_versions
    manifest.save(man)
    with connect(read_only=False) as con:
        for name, version, url, checksum, _changed in source_updates:
            entry = manifest.get(man, name)
            load.upsert_source(con, name, version=version, url=url, checksum=checksum,
                               license=entry.get("license", ""),
                               last_checked=entry.get("last_checked", ts),
                               last_updated=entry.get("last_updated", ts))
        load.append_findings(con, [f.row(ts) for f in findings])
        # Mark new PGS IDs seen so Phase 2's suggest-only watcher can diff later.
        load.mark_seen(con, "pgs_catalog", pgs_new_ids)
        # Mark PubMed + LitVar PMIDs seen so we never re-surface the same paper.
        load.mark_seen(con, "pubmed", pubmed_new_ids)
        load.mark_seen(con, "litvar", litvar_new_ids)
    report = _write_report(findings, ts)
    if report:
        console.print(f"[green]Wrote changelog[/] → {report}")
    console.print(f"[green]Refresh complete[/] — {len(findings)} new finding(s).")
    return findings


def _print_digest(findings: list[Finding], dry_run: bool) -> None:
    if not findings:
        console.print("[green]Nothing new.[/]")
        return
    verb = "Would surface" if dry_run else "New"
    console.print(f"\n[bold]{verb} ({len(findings)}):[/]")
    for f in findings:
        tag = {"strong": "[red]", "moderate": "[yellow]", "weak": "[dim]", "info": "[cyan]"}.get(f.tier, "")
        console.print(f"  {tag}[{f.tier}][/] {f.title}")

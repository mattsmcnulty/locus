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

# Re-annotate the full local-DB chain on a ClinVar refresh so SnpEff genes +
# AlphaMissense survive (each step copies the previous; running clinvar alone would
# strip them). gnomAD is deferred/streamed, so it's left out of the refresh chain.
_REANNOTATE_STEPS = "clinvar,snpeff,alphamissense"

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

    def row(self, ts: str) -> tuple:
        return (ts, self.source, self.kind, self.tier, self.chrom, self.pos, self.ref,
                self.alt, self.rsid, self.gene, self.title, self.detail,
                self.old_value, self.new_value, self.release)


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


# ── Orchestrator ─────────────────────────────────────────────────────────────────
ALL_SOURCES = ["clinvar", "pgs"]


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

    sources: comma-separated subset of {clinvar, pgs} or "all".
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

    # Rank strongest-first for the printed digest.
    findings.sort(key=lambda f: TIER_ORDER.index(f.tier) if f.tier in TIER_ORDER else len(TIER_ORDER))
    _print_digest(findings, dry_run)

    if dry_run:
        console.print("[dim]dry-run: no changes written.[/]")
        return findings

    # Persist: manifest (source of truth) + sources table + findings + report.
    for name, version, url, checksum, changed in source_updates:
        manifest.record(man, name, version=version, url=url, checksum=checksum, changed=changed)
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

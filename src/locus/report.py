"""A single self-contained HTML summary of the genome.

For the people who won't run a CLI or ask an LLM anything — the report is the whole product for
them. It is deliberately plain HTML with inline CSS and no scripts or external requests, so it
opens from a double-click, works offline forever, and can't phone anywhere.

It reuses ``queries`` rather than re-deriving anything, so it cannot drift from what the MCP tools
and the SPA say. It also repeats their caveats instead of trimming them: a static document gets
read without a conversation around it, so every place a reader could infer "nothing found = all
clear" has to say otherwise in the document itself.

The output contains personal genetic data. It is written under ``reports_dir`` (gitignored) and
should be shared with the same care as the VCF it came from.
"""

from __future__ import annotations

import datetime as _dt
import html
from pathlib import Path

from rich.console import Console

from . import queries
from .config import settings

console = Console()

_CSS = """
:root { color-scheme: light dark; }
* { box-sizing: border-box; }
body { font: 15px/1.55 -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
       margin: 0; padding: 2rem 1rem 4rem; background: #fbfbfd; color: #1d1d1f; }
@media (prefers-color-scheme: dark) { body { background: #131316; color: #e8e8ed; }
  .card { background: #1c1c21 !important; border-color: #2e2e36 !important; }
  th { background: #232329 !important; } code { background: #2a2a31 !important; } }
main { max-width: 860px; margin: 0 auto; }
h1 { font-size: 1.6rem; margin: 0 0 .25rem; } h2 { font-size: 1.15rem; margin: 2rem 0 .6rem; }
.sub { color: #6e6e73; margin: 0 0 1.5rem; font-size: .9rem; }
.card { background: #fff; border: 1px solid #e3e3e8; border-radius: 10px; padding: 1rem 1.15rem;
        margin: 0 0 1rem; }
table { border-collapse: collapse; width: 100%; font-size: .9rem; }
th, td { text-align: left; padding: .45rem .6rem; border-bottom: 1px solid #e3e3e880; vertical-align: top; }
th { background: #f5f5f7; font-weight: 600; }
.hint { color: #6e6e73; font-size: .85rem; }
.pill { display: inline-block; padding: .1rem .5rem; border-radius: 99px; font-size: .75rem; font-weight: 600; }
.ok { background: #e6f4ea; color: #137333; } .warn { background: #fef7e0; color: #8a6116; }
.info { background: #e8f0fe; color: #1a56b3; }
.note { border-left: 3px solid #c7c7cc; padding-left: .8rem; color: #6e6e73; font-size: .85rem; margin: .6rem 0 0; }
.disclaimer { border: 1px solid #f0c36d; background: #fffaf0; border-radius: 10px; padding: 1rem 1.15rem; }
@media (prefers-color-scheme: dark) { .disclaimer { background: #2a2318; border-color: #6b5a2e; } }
"""


def _e(v) -> str:
    return html.escape(str(v)) if v is not None else "—"


def _card(title: str, body: str) -> str:
    return f'<h2>{title}</h2>\n<div class="card">{body}</div>'


def _table(headers: list[str], rows: list[list[str]], empty: str) -> str:
    if not rows:
        return f'<p class="hint">{empty}</p>'
    head = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _overview() -> str:
    o = queries.overview()
    m = o.get("meta", {})
    steps = m.get("annotate_steps", "unknown")
    rows = [
        ["Variants", f"{o.get('variants', 0):,}"],
        ["ClinVar-annotated", f"{o.get('clinvar_annotated', 0):,}"],
        ["Pharmacogenomic genes", f"{o.get('pgx_genes', 0):,}"],
        ["Copy-number / structural", f"{o.get('cnv', 0):,} / {o.get('sv', 0):,}"],
        ["Annotations applied", _e(steps)],
        ["Built", _e(m.get("created_at"))],
    ]
    return _table(["", ""], rows, "No genome loaded.")


def _secondary() -> str:
    sf = queries.secondary_findings(limit=50)
    rows = [[_e(h.gene), _e(h.clnsig), _e((h.clndn or "").replace("_", " ")[:70]), _e(h.gt)]
            for h in sf.hits]
    t = _table(["Gene", "ClinVar", "Condition", "Genotype"], rows,
               '<span class="pill ok">No findings</span> — no pathogenic ClinVar variant in the '
               "84 ACMG-actionable genes. This is the common result.")
    return t + (
        '<p class="note"><strong>What this does not mean.</strong> "No findings" covers only '
        "pathogenic/likely-pathogenic ClinVar variants in 84 specific genes. It is not a clean "
        "bill of health: it cannot see variants ClinVar has never classified, non-coding or "
        "structural variants, repeat expansions, or any gene off that list.</p>")


def _carrier() -> str:
    c = queries.carrier_status()
    rows = [[_e(h.gene), _e(h.condition), f'<span class="pill {"warn" if h.status == "likely_affected" else "info"}">'
             f'{h.status.replace("_", " ")}</span>', _e(h.zygosity)] for h in c.hits]
    t = _table(["Gene", "Condition", "Status", "Zygosity"], rows,
               "No carrier findings in this panel.")
    na = _table(["Condition", "Why this data cannot answer it"],
                [[f"<strong>{_e(n.gene)}</strong> — {_e(n.condition)}", _e(n.why)] for n in c.not_assessed],
                "—")
    return (t + '<p class="note">Carrier status is about children, not your own health: one copy is '
            "typically silent. It matters when both partners carry the same condition (1-in-4 risk "
            "per pregnancy). This is a curated panel of common conditions, <strong>not a clinical "
            "carrier screen</strong> (ACMG's is 113 genes).</p>"
            f"<h3 style='font-size:.95rem;margin:1.2rem 0 .4rem'>Not assessed</h3>{na}"
            '<p class="note">These were <strong>not tested</strong> — their absence above is not a '
            "negative result. Two of them are among the most commonly offered carrier screens.</p>")


def _pgx() -> str:
    p = queries.pharmacogenomics()
    rows = [[_e(g.gene), _e(g.diplotype), _e(g.phenotype)] for g in p.genes]
    return _table(["Gene", "Diplotype", "Phenotype"], rows, "PharmCAT has not run.") + (
        '<p class="note">Genes reporting <code>Unknown</code> / <code>No Result</code> were '
        "<strong>not determined</strong> — not normal. CYP2D6 in particular needs copy-number data "
        "from raw reads that Locus does not process. Share this with a prescriber; do not change "
        "any medication on your own.</p>")


def _ancestry() -> str:
    try:
        a = queries.ancestry()
    except Exception:  # noqa: BLE001 - ancestry step may not have run
        return '<p class="hint">Ancestry has not been computed (run <code>locus ancestry</code>).</p>'
    rows = [[_e(c.name), f"{c.proportion * 100:.1f}%"] for c in a.components]
    return _table(["Population", "Proportion"], rows, "Not computed.") + \
        f'<p class="note">{_e(a.note)}</p>'


def _risk() -> str:
    rows = [[_e(s.trait), f"{s.percentile:.0f}th" if s.percentile is not None else "—",
             _e(s.ancestry), _e(s.pgs_id)] for s in queries.polygenic_risk()]
    return _table(["Trait", "Percentile", "Ancestry-matched", "Score"], rows,
                  "No polygenic scores (run <code>locus ancestry</code>).") + (
        '<p class="note">Research-grade estimates, meaningful only within a matched ancestry. A high '
        "percentile is not a diagnosis and a low one is not protection — lifestyle and family "
        "history usually matter more.</p>")


def _traits() -> str:
    rows = [[_e(t.trait), _e(t.genotype), _e(t.interpretation)] for t in queries.traits().traits]
    return _table(["Trait", "Genotype", "Interpretation"], rows, "Run <code>locus traits</code>.")


def _whats_new() -> str:
    wn = queries.whats_new(limit=15)
    rows = [[f'<span class="pill {"warn" if f.tier in ("strong", "moderate") else "info"}">{_e(f.tier)}</span>',
             _e(f.title[:90]), _e(f.source)] for f in wn.findings]
    return _table(["Tier", "Finding", "Source"], rows, "Nothing yet — run <code>locus refresh</code>.")


def build(dest: Path | None = None) -> Path:
    """Render the report. Returns the path written."""
    out = dest or (settings.reports_dir / "genome-report.html")
    out.parent.mkdir(parents=True, exist_ok=True)
    now = _dt.datetime.now().strftime("%d %B %Y, %H:%M")
    sections = [
        _card("Overview", _overview()),
        _card("Ancestry", _ancestry()),
        _card("ACMG secondary findings", _secondary()),
        _card("Carrier status", _carrier()),
        _card("Pharmacogenomics", _pgx()),
        _card("Polygenic risk", _risk()),
        _card("Traits &amp; wellness", _traits()),
        _card("What's new", _whats_new()),
    ]
    disclaimer = (
        '<div class="disclaimer"><strong>Please read.</strong> This is for personal exploration and '
        "education. It is <strong>not a medical device and not medical advice</strong>, and it is "
        "not a diagnostic test. Nothing here has been clinically confirmed. Discuss anything "
        "health-relevant with a qualified clinician or genetic counselor before acting on it — "
        "and note that an empty section means <em>nothing was found by these specific checks</em>, "
        "never that a condition has been ruled out.<br><br>"
        "This file contains your personal genetic data. It was generated on your own machine and "
        "sends nothing anywhere, but treat sharing it as you would sharing the genome itself.</div>"
    )
    body = "\n".join(sections)
    out.write_text(
        f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        f"<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>Locus — genome report</title><style>{_CSS}</style></head><body><main>"
        f"<h1>Genome report</h1><p class='sub'>Generated {now} · locally, from your own data</p>"
        f"{disclaimer}{body}"
        f"<p class='hint' style='margin-top:2rem'>Generated by Locus. Ask Claude for anything this "
        f"doesn't cover.</p></main></body></html>"
    )
    console.print(f"[green]Report written[/] → {out}")
    return out

"""Literature watcher — keep the genome current with new PubMed papers + GWAS studies.

Two capabilities, both scoped to the variants *this* genome carries:

  1. A weekly watcher (driven by ``locus refresh``): new PubMed papers mentioning genes you
     carry notable variants in, surfaced into the ``watch_findings`` changelog. (The GWAS
     re-analysis watcher lives in ``refresh.py`` next to the ClinVar one.)
  2. Interactive lookups (MCP ``literature_for`` / ``variants_in_study``): "the latest research
     on my BRCA2", or "this new paper (PMID) reported N variants — which do I carry?".

Privacy: only generic identifiers leave the machine — gene symbols and rsIDs to NCBI
E-utilities and the GWAS Catalog REST API. Genotypes never leave; all matching is local.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass

from rich.console import Console

console = Console()

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_GWAS_REST = "https://www.ebi.ac.uk/gwas/rest/api"
_TOOL = "locus"  # NCBI etiquette: identify the tool. No personal email is sent.
_MIN_INTERVAL = 0.34  # ≤3 requests/second without an API key.
_last_call = 0.0


@dataclass
class PubMedHit:
    pmid: str
    title: str
    abstract: str = ""
    journal: str = ""
    year: str = ""
    gene: str = ""

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"


# ── HTTP helpers (generic queries only; best-effort) ────────────────────────────
def _throttle() -> None:
    global _last_call
    import time as _t

    wait = _MIN_INTERVAL - (_t.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = _t.monotonic()


def _eutils(endpoint: str, params: dict, *, want_json: bool):
    import httpx

    _throttle()
    params = {**params, "tool": _TOOL}
    try:
        r = httpx.get(f"{_EUTILS}/{endpoint}", params=params, timeout=30, follow_redirects=True)
        r.raise_for_status()
        return r.json() if want_json else r.text
    except Exception as e:  # noqa: BLE001 - network is best-effort; report & continue
        console.print(f"[yellow]PubMed probe failed[/] ({endpoint}): {e}")
        return None


def _gwas_get(path_or_url: str):
    import httpx

    url = path_or_url if path_or_url.startswith("http") else f"{_GWAS_REST}{path_or_url}"
    try:
        r = httpx.get(url, timeout=30, follow_redirects=True,
                      headers={"Accept": "application/json"})
        r.raise_for_status()
        return r.json()
    except Exception as e:  # noqa: BLE001
        console.print(f"[yellow]GWAS REST probe failed[/] {url}: {e}")
        return None


# ── PubMed ──────────────────────────────────────────────────────────────────────
def _fmt_date(iso: str | None) -> str | None:
    """ISO8601 / 'YYYY-MM-DD' → PubMed 'YYYY/MM/DD' (esearch date format)."""
    if not iso:
        return None
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", iso)
    return f"{m.group(1)}/{m.group(2)}/{m.group(3)}" if m else None


def _parse_pubmed_xml(xml_text: str, gene: str = "") -> list[PubMedHit]:
    """Parse an efetch PubmedArticleSet into PubMedHits (title + abstract + journal + year)."""
    import xml.etree.ElementTree as ET

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        console.print(f"[yellow]PubMed XML parse failed:[/] {e}")
        return []
    hits: list[PubMedHit] = []
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID") or ""
        if not pmid:
            continue
        title = "".join(art.find(".//Article/ArticleTitle").itertext()).strip() \
            if art.find(".//Article/ArticleTitle") is not None else ""
        abstract = " ".join(
            "".join(node.itertext()).strip()
            for node in art.findall(".//Article/Abstract/AbstractText")
        ).strip()
        journal = art.findtext(".//Article/Journal/Title") or ""
        year = (art.findtext(".//Article/Journal/JournalIssue/PubDate/Year")
                or art.findtext(".//Article/Journal/JournalIssue/PubDate/MedlineDate") or "")[:4]
        hits.append(PubMedHit(pmid=pmid, title=title, abstract=abstract,
                              journal=journal, year=year, gene=gene))
    return hits


def pubmed_search(term: str, *, mindate: str | None = None, maxdate: str | None = None,
                  retmax: int = 20, gene: str = "") -> list[PubMedHit]:
    """Search PubMed and return hydrated hits (title + abstract). ``mindate``/``maxdate`` are
    ISO dates that scope by publication date (used by the weekly watcher)."""
    params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": retmax, "sort": "date"}
    lo, hi = _fmt_date(mindate), _fmt_date(maxdate)
    if lo:
        params.update({"datetype": "pdat", "mindate": lo, "maxdate": hi or "3000/01/01"})
    res = _eutils("esearch.fcgi", params, want_json=True)
    ids = (((res or {}).get("esearchresult") or {}).get("idlist")) or []
    if not ids:
        return []
    xml_text = _eutils("efetch.fcgi", {"db": "pubmed", "id": ",".join(ids),
                                       "retmode": "xml", "rettype": "abstract"}, want_json=False)
    if not xml_text:
        return []
    hits = _parse_pubmed_xml(xml_text, gene=gene)
    order = {pmid: i for i, pmid in enumerate(ids)}
    hits.sort(key=lambda h: order.get(h.pmid, len(ids)))
    return hits


# A gene must co-occur with an actual genetics/variant word in the title or abstract — not just be
# gene2pubmed-linked. Deliberately [tiab] (title/abstract), NOT [MeSH]: the weekly watcher queries
# the last ~week, and brand-new papers aren't MeSH-indexed yet, so a MeSH filter would silently drop
# most recent hits (measured: ~⅔ of the last fortnight's BRCA2 papers). The old catch-all clause
# `(variant OR mutation OR ... OR clinical)` matched ~everything ("clinical"/"variant" are ubiquitous)
# — this cut a noisy gene like MET from ~11k papers to a few hundred while keeping BRCA2's real signal.
_GENE_GENETICS_TIAB = (
    '(mutation[tiab] OR variant[tiab] OR variants[tiab] OR polymorphism[tiab] '
    'OR genotype[tiab] OR allele[tiab] OR "loss of function"[tiab] OR pathogenic[tiab] '
    'OR "genetic variant"[tiab])'
)


def _gene_term(gene: str) -> str:
    """High-precision PubMed term for a gene symbol: the gene entity AND a genetics word in the
    title/abstract. Used by both the weekly watcher and the interactive lookup so they stay in sync."""
    return f"{gene}[Gene] AND {_GENE_GENETICS_TIAB}"


def _mentions_gene(hit: PubMedHit, gene: str) -> bool:
    """Backstop for symbol-collision / pure-linkage noise: keep a hit only if the gene symbol
    actually appears (whole word) in its title or abstract. Catches papers that [Gene] links but
    never name — including recent ones that arrive before MeSH indexing could have filtered them.

    Case-SENSITIVE on purpose: HGNC symbols are uppercase in text (MET, BRCA2, C9orf72), whereas the
    English words they collide with are lowercase ("met", "set"). Measured on a noisy gene, this drops
    the false hits (MET: 5→1) with no loss on real ones (BRCA2/RIPK1 unchanged). A human-genome tool
    can accept missing the rare rodent-cased ortholog ("Lgals3") for that precision."""
    pattern = re.compile(rf"\b{re.escape(gene)}\b")
    return bool(pattern.search(hit.title or "") or pattern.search(hit.abstract or ""))


def _term_for(query: str) -> str:
    """Build a PubMed term from a free query: an rsID stays literal, a bare gene symbol gets the
    high-precision gene term, anything else passes through."""
    q = query.strip()
    if re.fullmatch(r"rs\d+", q, re.IGNORECASE):
        return q
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,9}", q):
        return _gene_term(q)
    return q


def literature_for(query: str, *, since: str | None = None, retmax: int = 20) -> list[PubMedHit]:
    """On-demand PubMed lookup for a gene / rsID / free-text query (interactive MCP tool)."""
    gene = query.strip() if re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,9}", query.strip()) else ""
    return pubmed_search(_term_for(query), mindate=since, retmax=retmax, gene=gene)


def fetch_pubmed_meta(pmids: list[str]) -> dict[str, PubMedHit]:
    """Hydrate a batch of PMIDs into PubMedHits (title + abstract) via efetch → {pmid: hit}."""
    pmids = [p for p in pmids if p]
    if not pmids:
        return {}
    xml_text = _eutils("efetch.fcgi", {"db": "pubmed", "id": ",".join(pmids),
                                       "retmode": "xml", "rettype": "abstract"}, want_json=False)
    return {h.pmid: h for h in (_parse_pubmed_xml(xml_text) if xml_text else [])}


# ── LitVar2 — variant-centric literature (papers about a specific rsID) ───────────
_LITVAR = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"


def _litvar_get(path: str):
    import httpx

    _throttle()
    try:
        r = httpx.get(f"{_LITVAR}{path}", timeout=30, follow_redirects=True,
                      headers={"Accept": "application/json"})
        if r.status_code in (400, 404):
            return None  # rsID not in LitVar's index (no papers cite it) — expected, stay quiet
        r.raise_for_status()
        return r.json()
    except Exception as e:  # noqa: BLE001 - network is best-effort; report & continue
        console.print(f"[yellow]LitVar probe failed[/] ({path}): {e}")
        return None


def litvar_pmids(rsid: str, recent: int = 25) -> list[str]:
    """The most-recent PubMed IDs mentioning a specific rsID, via NCBI LitVar2.

    PMIDs are ~chronological, so we keep the highest (newest) as a bound — enough to catch papers
    published between weekly runs without pulling a heavily-studied variant's entire back-catalogue.
    """
    if not rsid.lower().startswith("rs"):
        return []
    vid = f"litvar%40{rsid}%23%23"  # LitVar variant id form: litvar@rs...##
    data = _litvar_get(f"/variant/get/{vid}/publications")
    pmids = [str(p) for p in ((data or {}).get("pmids") or [])]
    pmids.sort(key=lambda p: int(p) if p.isdigit() else 0, reverse=True)
    return pmids[:recent]


# ── GWAS Catalog study → variants ────────────────────────────────────────────────
def study_rsids(pmid: str) -> list[str]:
    """Every rsID reported by the GWAS Catalog study/studies for a PubMed ID."""
    data = _gwas_get(f"/studies/search/findByPublicationIdPubmedId?pubmedId={pmid}")
    studies = ((data or {}).get("_embedded") or {}).get("studies") or []
    rsids: set[str] = set()
    for st in studies:
        href = (((st.get("_links") or {}).get("associations") or {}).get("href"))
        acc = st.get("accessionId")
        assoc = _gwas_get(href) if href else (_gwas_get(f"/studies/{acc}/associations") if acc else None)
        if assoc:
            rsids.update(re.findall(r"rs\d+", json.dumps(assoc)))
    return sorted(rsids, key=lambda r: int(r[2:]))


def study_variants(pmid: str) -> dict:
    """Given a PubMed ID, pull the study's reported rsIDs and genotype *this* genome at them."""
    from . import gwas

    rsids = study_rsids(pmid)
    if not rsids:
        return {"pmid": pmid, "total": 0, "carried": 0, "markers": [],
                "note": "No GWAS Catalog study/variants found for that PubMed ID."}
    markers = gwas.ask_markers(rsids)
    if not markers:
        # The rsIDs came from the GWAS Catalog, so a genotyping failure (Ensembl down, bcftools
        # error) yields an empty list while `len(rsids)` stays large — which would otherwise be
        # reported as a confident "you carry 0 of them". "We couldn't look it up" and "you carry
        # none" are opposite answers; never let the first masquerade as the second.
        return {"pmid": pmid, "total": len(rsids), "carried": 0, "markers": [],
                "note": (f"This study reported {len(rsids)} variant(s), but genotyping them failed "
                         f"(the Ensembl lookup returned nothing). This is NOT a result of 'you carry "
                         f"none' — it is unknown. Try again shortly.")}
    carried = sum(1 for m in markers if _is_carried(m.get("genotype"), m.get("ref")))
    unresolved = len(rsids) - len(markers)
    note = (f"{len(rsids)} variant(s) reported by this study; genotyped {len(markers)}; you carry a "
            f"non-reference allele at {carried}. Hom-ref-aware live genotyping (Ensembl GRCh38). "
            + (f"{unresolved} could not be resolved and are unknown, not absent. " if unresolved else "")
            + "Single-variant evidence — informational, not diagnostic.")
    return {"pmid": pmid, "total": len(rsids), "carried": carried, "markers": markers, "note": note}


def _is_carried(genotype: str | None, ref: str | None) -> bool:
    """True if any called allele differs from the reference (heterozygous or homozygous-alt).

    Needs ``ref`` to tell a homozygous-alt call apart from a homozygous-reference one; without it
    only heterozygosity is detectable."""
    if not genotype or genotype in ("—", "."):
        return False
    alleles = [a for a in re.split(r"[/|]", genotype) if a and a != "."]
    if not alleles:
        return False
    if ref:
        return any(a != ref for a in alleles)
    return len(set(alleles)) > 1

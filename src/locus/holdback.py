"""Held-back topics: results the person has chosen not to see (e.g. Alzheimer's / APOE).

Opt-in and empty by default. Once a topic is held back, every surface that reads the store hides
it: each read-only connection shadows the user-facing tables with filtered TEMP views (so the
MCP tools, the API/SPA, the HTML report and ``run_sql`` all see the same filtered store), and the
paths that never touch the store — live genotyping, PubMed, the weekly report — check
:func:`is_held` directly.

Hidden is not deleted: the pipeline and the weekly refresh still write everything, so removing
a topic later brings it all back. This is a "don't show me" preference that protects against
*stumbling onto* a result, not a security boundary against someone who has the files.

Held rows are hidden whether the match is by gene, by rsID, or by a trait/disease keyword. A
gene's rsIDs are hidden for EVERY trait: APOE rs429358 is also a top LDL-cholesterol hit, and
showing that genotype would reveal the very thing held back.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from .config import settings


@dataclass(frozen=True)
class Topic:
    name: str
    label: str
    genes: tuple[str, ...] = ()
    rsids: tuple[str, ...] = ()
    keywords: tuple[str, ...] = ()  # matched case-insensitively as substrings of trait/disease text


# Curated presets. Conservative: genes whose results mostly ARE the topic, and keywords specific
# enough not to swallow unrelated traits (no bare "als" — it matches "also").
PRESETS: dict[str, Topic] = {t.name: t for t in (
    Topic("alzheimers", "Alzheimer's disease & dementia (incl. APOE)",
          genes=("APOE", "APP", "PSEN1", "PSEN2", "TREM2", "SORL1"),
          rsids=("rs429358", "rs7412"),  # the two SNPs that define APOE e2/e3/e4
          keywords=("alzheimer", "dementia", "cognitive decline")),
    Topic("huntington", "Huntington's disease", genes=("HTT",), keywords=("huntington",)),
    Topic("parkinsons", "Parkinson's disease",
          genes=("LRRK2", "SNCA", "PRKN", "PINK1", "PARK7", "VPS35"), keywords=("parkinson",)),
    Topic("als", "ALS / motor neuron disease",
          genes=("SOD1", "C9ORF72", "TARDBP", "FUS"),
          keywords=("amyotrophic lateral sclerosis", "motor neuron disease", "frontotemporal")),
    Topic("hereditary-cancer", "Hereditary cancer risk (BRCA, Lynch, …)",
          genes=("BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2", "EPCAM", "APC", "TP53",
                 "CDH1", "STK11", "PTEN", "CHEK2", "ATM"),
          keywords=("hereditary breast", "ovarian cancer", "lynch syndrome", "li-fraumeni")),
    Topic("psychiatric", "Psychiatric conditions",
          keywords=("schizophrenia", "bipolar", "major depressive", "depression", "autism", "adhd",
                    "attention deficit", "anorexia nervosa", "obsessive-compulsive", "suicid")),
)}

_GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9._-]*$")
_RSID_RE = re.compile(r"^rs\d+$")
_KW_RE = re.compile(r"^[a-z0-9][a-z0-9 '\-]*$")


@dataclass
class HoldbackSpec:
    """The merged set of everything currently held back."""
    topics: list[Topic] = field(default_factory=list)

    @property
    def genes(self) -> set[str]:
        return {g for t in self.topics for g in t.genes}

    @property
    def rsids(self) -> set[str]:
        return {r for t in self.topics for r in t.rsids}

    @property
    def keywords(self) -> set[str]:
        return {k for t in self.topics for k in t.keywords}

    def __bool__(self) -> bool:
        return bool(self.topics)


# ── Persistence (data/holdback.json — under data/, so gitignored and per-machine) ──────────────
def config_path() -> Path:
    return settings.data_dir / "holdback.json"


def _read() -> dict:
    p = config_path()
    if not p.exists():
        return {"topics": [], "custom": {}}
    try:
        data = json.loads(p.read_text() or "{}")
    except json.JSONDecodeError as e:
        # Fail CLOSED-ish and loudly: a corrupt file must not silently un-hide everything.
        raise RuntimeError(f"{p} is not valid JSON ({e}); fix or delete it.") from e
    return {"topics": list(data.get("topics", [])), "custom": dict(data.get("custom", {}))}


def _write(data: dict) -> None:
    p = config_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, indent=2) + "\n")


def load() -> HoldbackSpec:
    data = _read()
    topics = [PRESETS[n] for n in data["topics"] if n in PRESETS]
    for name, c in data["custom"].items():
        topics.append(Topic(name=name, label=c.get("label") or name, genes=tuple(c.get("genes", ())),
                            rsids=tuple(c.get("rsids", ())), keywords=tuple(c.get("keywords", ()))))
    return HoldbackSpec(topics)


def add_preset(name: str) -> Topic:
    if name not in PRESETS:
        raise ValueError(f"Unknown topic '{name}'. Presets: {', '.join(PRESETS)} (or use --gene/--rsid/--keyword).")
    data = _read()
    if name not in data["topics"]:
        data["topics"].append(name)
        _write(data)
    return PRESETS[name]


def add_custom(name: str, *, genes=(), rsids=(), keywords=(), label: str | None = None) -> Topic:
    genes = tuple(g.strip().upper() for g in genes if g.strip())
    rsids = tuple(r.strip().lower() for r in rsids if r.strip())
    keywords = tuple(k.strip().lower() for k in keywords if k.strip())
    bad = ([g for g in genes if not _GENE_RE.match(g)] + [r for r in rsids if not _RSID_RE.match(r)]
           + [k for k in keywords if not _KW_RE.match(k)])
    if bad:
        raise ValueError(f"Not a valid gene / rsID / keyword: {', '.join(bad)}")
    if not (genes or rsids or keywords):
        raise ValueError("A custom topic needs at least one --gene, --rsid or --keyword.")
    if name in PRESETS:
        raise ValueError(f"'{name}' is a preset name; pick another name for a custom topic.")
    data = _read()
    data["custom"][name] = {"label": label or name, "genes": list(genes), "rsids": list(rsids),
                            "keywords": list(keywords)}
    _write(data)
    return Topic(name, label or name, genes, rsids, keywords)


def remove(name: str) -> bool:
    data = _read()
    found = name in data["topics"] or name in data["custom"]
    data["topics"] = [t for t in data["topics"] if t != name]
    data["custom"].pop(name, None)
    if found:
        _write(data)
    return found


# ── Gene footprints (genotype-INDEPENDENT positions) ──────────────────────────────────────────
def gene_spans(spec: HoldbackSpec | None = None) -> dict[str, tuple[str, int, int]]:
    """gene -> (chrom, start, end) for held genes, from the local ClinVar file's GENEINFO.

    Hiding by position catches variants SnpEff labels with a neighbouring gene and live-genotyped
    markers that have no store row. It must not depend on the person's genotype — otherwise
    "hidden" vs "shown" would itself reveal whether they carry something — so the spans come from
    ClinVar (the same for everyone), never from this genome's variant table. Cached per gene set.
    """
    import subprocess

    spec = load() if spec is None else spec
    genes = spec.genes
    clinvar = settings.annotations_dir / "clinvar.chr.vcf.gz"
    if not genes or not clinvar.exists():
        return {}
    cache = settings.work_dir / "holdback_spans.json"
    key = {"genes": sorted(genes), "clinvar_mtime": clinvar.stat().st_mtime}
    if cache.exists():
        try:
            cached = json.loads(cache.read_text())
            if cached.get("key") == key:
                return {g: tuple(v) for g, v in cached["spans"].items()}
        except (json.JSONDecodeError, KeyError):
            pass
    hits: dict[str, dict[str, list[int]]] = {}
    proc = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO/GENEINFO\n", str(clinvar)],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return {}
    for line in proc.stdout.splitlines():
        chrom, pos, info = line.split("\t", 2)
        for entry in info.split("|"):
            g = entry.split(":", 1)[0].upper()
            if g in genes:
                hits.setdefault(g, {}).setdefault(chrom, []).append(int(pos))
    spans = {}
    for g, by_chrom in hits.items():
        chrom, ps = max(by_chrom.items(), key=lambda kv: len(kv[1]))
        spans[g] = (chrom, min(ps), max(ps))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"key": key, "spans": spans}))
    return spans


# ── Matching (Python side — for paths that don't read the store) ──────────────────────────────
def _kw_regex(spec: HoldbackSpec) -> re.Pattern | None:
    kws = sorted(spec.keywords)
    return re.compile("|".join(re.escape(k) for k in kws), re.IGNORECASE) if kws else None


def is_held(*, gene: str | None = None, rsid: str | None = None, text: str | None = None,
            chrom: str | None = None, pos: int | None = None, spec: HoldbackSpec | None = None) -> bool:
    """Does this gene / rsID / position / free text (trait, disease, title) touch a held-back topic?"""
    spec = load() if spec is None else spec
    if not spec:
        return False
    if chrom and pos is not None:
        c = chrom if chrom.startswith("chr") else ("chrM" if chrom == "MT" else f"chr{chrom}")
        if any(c == sc and a <= int(pos) <= b for sc, a, b in gene_spans(spec).values()):
            return True
    if gene and any(g.strip().upper() in spec.genes for g in re.split(r"[,;&/ ]+", gene) if g.strip()):
        return True
    if rsid and rsid.strip().lower() in spec.rsids:
        return True
    rx = _kw_regex(spec)
    return bool(text and rx and rx.search(text))


def query_is_held(query: str, spec: HoldbackSpec | None = None) -> bool:
    """Does a free-form request (gene name, rsIDs, or a topic) directly ask for held-back results?"""
    spec = load() if spec is None else spec
    if not spec:
        return False
    tokens = re.findall(r"[A-Za-z0-9._-]+", query)
    return (any(is_held(gene=t, rsid=t, spec=spec) for t in tokens) or is_held(text=query, spec=spec))


def note(spec: HoldbackSpec | None = None) -> str | None:
    spec = load() if spec is None else spec
    if not spec:
        return None
    labels = "; ".join(t.label for t in spec.topics)
    return (f"The person has chosen to hold back results about: {labels}. Every Locus tool already "
            "filters these out. Do not try to reconstruct them from other results, and do not treat "
            "their absence as a negative result. If asked about one, say it is held back by their own "
            "setting and that `locus holdback remove <topic>` shows it again.")


def held_message(what: str) -> str:
    return (f"{what} is held back by this person's own Locus setting, so no result is shown — this is "
            "NOT a negative or homozygous-reference result. `locus holdback list` shows what's held back; "
            "`locus holdback remove <topic>` shows it again.")


# ── SQL side: filtered TEMP views shadowing the user-facing tables ────────────────────────────
def _sql_list(values: set[str]) -> str:
    return ", ".join("'" + v.replace("'", "''") + "'" for v in sorted(values)) or "NULL"


def _held_predicates(spec: HoldbackSpec) -> dict[str, str]:
    """table -> SQL boolean expression that is TRUE for held rows (NULL-safe)."""
    genes, rsids = _sql_list(spec.genes), _sql_list(spec.rsids)
    kws = sorted(spec.keywords)
    rx = "|".join(re.escape(k) for k in kws).replace("'", "''")

    def gene_in(col: str) -> str:
        return f"coalesce(upper({col}) IN ({genes}), false)"

    def gene_list_has(col: str) -> str:  # comma/space-separated gene lists (CNV/SV, PGx)
        return (f"coalesce(len(list_intersect(string_split_regex(upper({col}), '[^A-Z0-9._-]+'), "
                f"[{genes}])) > 0, false)" if spec.genes else "false")

    def rsid_in(col: str) -> str:
        return f"coalesce(lower({col}) IN ({rsids}), false)"

    def kw(*cols: str) -> str:
        if not kws:
            return "false"
        return " OR ".join(f"coalesce(regexp_matches(lower({c}), '{rx}'), false)" for c in cols)

    spans = list(gene_spans(spec).values())

    def in_span(chrom: str = "chrom", pos: str = "pos", end: str | None = None) -> str:
        if not spans:
            return "false"
        # Store contigs are chr-prefixed except CNV/SV, which keep the caller's naming — match both.
        terms = [f"(({chrom} = '{c}' OR {chrom} = '{c.removeprefix('chr')}') AND "
                 + (f"{pos} <= {b} AND coalesce({end}, {pos}) >= {a})" if end else f"{pos} BETWEEN {a} AND {b})")
                 for c, a, b in spans]
        return "coalesce(" + " OR ".join(terms) + ", false)"

    sv_span = in_span(end='"end"')
    return {
        "variants": f"{gene_in('gene')} OR {rsid_in('rsid')} OR {kw('clndn')} OR {in_span()}",
        "associations": f"{rsid_in('rsid')} OR {kw('trait', 'mapped_trait')} OR {in_span()}",
        "traits": f"{rsid_in('rsid')} OR {kw('trait', 'interpretation')}",
        "pgs_scores": kw("trait"),
        "personality_pgi": kw("trait"),
        "cnv": f"{gene_list_has('genes')} OR {sv_span}",
        "sv": f"{gene_list_has('genes')} OR {sv_span}",
        "pgx_genes": gene_in("gene"),
        "pgx_drugs": f"{gene_list_has('gene')} OR {kw('drug', 'recommendation')}",
        "watch_findings": (f"{gene_list_has('gene')} OR {rsid_in('rsid')} OR {kw('title', 'detail')} "
                           f"OR {in_span()}"),
    }


def install_views(con, spec: HoldbackSpec | None = None) -> None:
    """Shadow each user-facing table with a TEMP view that drops held rows (read-only connections)."""
    spec = load() if spec is None else spec
    if not spec:
        return
    db = con.execute("SELECT current_database()").fetchone()[0]
    present = {r[0] for r in con.execute(
        "SELECT table_name FROM duckdb_tables() WHERE database_name = current_database() "
        "AND schema_name = 'main'").fetchall()}
    for table, held in _held_predicates(spec).items():
        if table in present:
            con.execute(f'CREATE OR REPLACE TEMP VIEW {table} AS '
                        f'SELECT * FROM "{db}".main.{table} WHERE NOT ({held})')


def check_sql(sql: str, db_name: str = "locus", spec: HoldbackSpec | None = None) -> None:
    """Refuse raw SQL that would step around the filtered views."""
    spec = load() if spec is None else spec
    if not spec:
        return
    if re.search(rf'(\b{re.escape(db_name)}|"{re.escape(db_name)}")\s*\.', sql, re.IGNORECASE):
        raise ValueError("Catalog-qualified table names aren't allowed while topics are held back — "
                         "use plain table names (e.g. `variants`).")
    if re.search(r"\b(query|query_table)\s*\(", sql, re.IGNORECASE):
        raise ValueError("query()/query_table() aren't allowed while topics are held back.")

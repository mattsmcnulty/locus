"""ClinGen Gene-Disease Validity — is a gene *itself* established as causing a disease?

This answers the gap the variant-level watchers can't: a new study makes gene X notable — a new
gene-disease association, or an existing one strengthened/refuted — *without* reclassifying any
variant you carry. ClinVar tracks variants; GWAS tracks associations; LitVar/PubMed track papers.
None of them fire on "the gene you carry a variant in was just established as disease-causing".

ClinGen curates exactly that: expert panels assign each gene-disease pair a validity classification
(Definitive → Strong → Moderate → Limited → Disputed → Refuted). It's a small, freely-downloadable
CSV (no API key), and the assertions are structured and rare-to-change — so watching it is high
signal, unlike a raw literature sweep.

The watcher (in refresh.py) scopes to genes where you carry a *non-benign* variant, because that's
the only place a gene-disease change means something for you: a gene gaining a definitive disease
link matters if you carry an uncertain variant in it, and is noise if you carry only a benign one.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .config import settings

# Strongest → weakest. Used to detect upgrades/downgrades and to sort interactive lookups.
CLASSIFICATION_RANK = {
    "Definitive": 6,
    "Strong": 5,
    "Moderate": 4,
    "Limited": 3,
    "Animal Model Only": 2,
    "Disputed": 1,
    "Refuted": 0,
    "No Known Disease Relationship": 0,
}
# Classifications strong enough that a *new* one is worth surfacing.
ESTABLISHED = frozenset({"Definitive", "Strong", "Moderate"})
# A gene-disease link being actively doubted is itself notable.
REFUTING = frozenset({"Disputed", "Refuted"})


@dataclass(frozen=True)
class Assertion:
    gene: str
    disease: str
    mondo: str            # MONDO disease id — the stable key, labels get reworded
    moi: str              # mode of inheritance: AD | AR | XL | …
    classification: str
    date: str             # classification date (ISO)
    url: str              # the online curation report
    panel: str = ""       # gene curation expert panel

    @property
    def key(self) -> tuple[str, str]:
        return (self.gene, self.mondo)

    @property
    def rank(self) -> int:
        return CLASSIFICATION_RANK.get(self.classification, 0)


def _csv_path() -> Path:
    return settings.annotations_dir / "clingen" / "gene-disease-summary.csv"


def parse(path: Path) -> list[Assertion]:
    """Parse a ClinGen Gene-Disease-Summary CSV.

    The file has a few banner rows, then a ``GENE SYMBOL,…`` header, then a ``+++`` separator, then
    data. We key on the header row rather than a fixed line number so a change in banner text
    doesn't silently shift the columns.
    """
    if not path.exists():
        return []
    def cell(row: list[str], cols: dict[str, int], name: str) -> str:
        i = cols.get(name)
        return row[i].strip() if i is not None and i < len(row) else ""

    out: list[Assertion] = []
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        cols: dict[str, int] | None = None
        for row in csv.reader(fh):
            if not row or not row[0] or row[0].startswith("+++"):
                continue
            if cols is None:
                if row[0].strip().upper() == "GENE SYMBOL":
                    cols = {name.strip().upper(): i for i, name in enumerate(row)}
                continue
            gene, cls = cell(row, cols, "GENE SYMBOL"), cell(row, cols, "CLASSIFICATION")
            if not gene or not cls:
                continue
            out.append(Assertion(
                gene=gene, disease=cell(row, cols, "DISEASE LABEL"),
                mondo=cell(row, cols, "DISEASE ID (MONDO)"), moi=cell(row, cols, "MOI"),
                classification=cls, date=cell(row, cols, "CLASSIFICATION DATE"),
                url=cell(row, cols, "ONLINE REPORT"), panel=cell(row, cols, "GCEP"),
            ))
    return out


def load_assertions(path: Path | None = None) -> list[Assertion]:
    """All assertions from the downloaded ClinGen CSV (empty if not downloaded yet)."""
    return parse(path or _csv_path())


def for_gene(gene: str) -> list[Assertion]:
    """ClinGen's curated disease associations for one gene, strongest-validity first.

    Reads the local CSV — no network, nothing leaves the machine (the gene name never even goes
    out). Used by the interactive `gene_disease_validity` MCP tool.
    """
    g = gene.strip().upper()
    hits = [a for a in load_assertions() if a.gene.upper() == g]
    hits.sort(key=lambda a: (-a.rank, a.disease))
    return hits

"""Polygenic scores (PGS Catalog) — Phase 1.

The sample's raw score for a PGS is simply ``sum(effect_allele_count × weight)``
over the score's markers. We compute that directly in Python from the
gVCF-expanded genotypes (so hom-ref positions count as 0 effect alleles, not as
missing) and report exact coverage. The reference distribution (for the
ancestry-matched percentile) is scored separately with PLINK2 over the panel.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from cyvcf2 import VCF
from rich.console import Console

from . import ancestry
from .config import settings

console = Console()

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass
class ScoreVariant:
    chrom: str       # chr-prefixed
    pos: int
    effect_allele: str
    other_allele: str
    weight: float


@dataclass
class SampleScore:
    pgs_id: str
    raw: float
    n_total: int      # variants in the score
    n_used: int       # variants with a confident sample genotype
    coverage: float   # n_used / n_total
    # (chrom_nochr, pos) -> (effect_allele, weight), for the variants actually used (to score the panel).
    used: dict[tuple[str, int], tuple[str, float]] = None  # type: ignore[assignment]
    # (chrom_nochr, pos) -> sample effect-allele count, so the raw score can be re-summed over a subset.
    counts: dict[tuple[str, int], int] = None  # type: ignore[assignment]

    def raw_over(self, keys) -> float:
        """The sample's score summed over just ``keys`` (a subset of ``used``)."""
        return sum(self.counts[k] * self.used[k][1] for k in keys)


@dataclass
class RefComparison:
    """Where the sample sits in its ancestry-matched reference — computed over ONE variant set.

    PLINK drops some of the sample's variants when scoring the panel (multiallelic / non-SNP /
    duplicate positions, or absent from the panel). Comparing a sample sum over N variants against
    panel sums over N−k is apples-to-oranges, so ``raw``/``n_used`` are re-derived over exactly the
    variants the panel kept.
    """
    raw: float
    n_used: int
    percentile: float | None
    z: float | None


# A small curated, high-value starter set (GRCh38-harmonized PGS Catalog IDs).
CURATED_PGS = [
    ("PGS000018", "Coronary artery disease"),   # genome-wide metaGRS (~1.7M variants; slow)
    ("PGS000065", "LDL cholesterol"),
    ("PGS000868", "Type 2 diabetes"),
    ("PGS000192", "Cholesterol (LDL/HDL, 9-SNP)"),
]


@dataclass
class CalibratedScore:
    pgs_id: str
    trait: str
    raw: float
    percentile: float | None     # within the ancestry-matched reference (0-100)
    ancestry: str | None         # superpop the percentile is relative to
    n_used: int
    coverage: float


def parse_scoring_file(path: Path) -> tuple[list[ScoreVariant], dict]:
    """Parse a PGS Catalog harmonized (hmPOS_GRCh38) scoring file."""
    meta: dict = {}
    variants: list[ScoreVariant] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header: list[str] | None = None
        for line in fh:
            if line.startswith("#"):
                if "=" in line:
                    k, _, v = line[1:].strip().partition("=")
                    meta[k] = v
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            row = dict(zip(header, cols, strict=False))
            chrom = row.get("hm_chr") or row.get("chr_name")
            pos = row.get("hm_pos") or row.get("chr_position")
            ea = row.get("effect_allele")
            wt = row.get("effect_weight")
            if not (chrom and pos and ea and wt):
                continue
            try:
                variants.append(ScoreVariant(
                    chrom=chrom if chrom.startswith("chr") else f"chr{chrom}",
                    pos=int(pos), effect_allele=ea.upper(),
                    other_allele=(row.get("other_allele") or "").upper(),
                    weight=float(wt),
                ))
            except ValueError:
                continue
    return variants, meta


def _effect_allele_count(ref: str, alts: list[str], gt, ea: str, oa: str) -> int | None:
    """Count copies of the effect allele in a sample genotype, handling strand flip."""
    a, b = gt[0], gt[1]
    if a < 0 or b < 0:
        return None  # no-call
    alleles = [ref] + list(alts)

    def allele(i: int) -> str:
        return alleles[i] if 0 <= i < len(alleles) else ""

    called = [allele(a), allele(b)]
    # Direct match.
    if ea in called or oa in called or ea == ref or oa == ref:
        return sum(1 for c in called if c == ea)
    # Strand flip (effect/other are on the opposite strand).
    eaf, oaf = _COMPLEMENT.get(ea, ea), _COMPLEMENT.get(oa, oa)
    if eaf in called or oaf in called or eaf == ref:
        return sum(1 for c in called if c == eaf)
    return None  # alleles don't reconcile — exclude from the score


def score_sample(scoring_file: Path, pgs_id: str) -> SampleScore:
    """Compute the sample's raw PGS = sum(effect_allele_count × weight)."""
    variants, _ = parse_scoring_file(scoring_file)
    settings.work_dir.mkdir(parents=True, exist_ok=True)
    bed = settings.work_dir / f"{pgs_id}.markers.bed"
    bed.write_text("".join(f"{v.chrom}\t{v.pos - 1}\t{v.pos}\n" for v in variants))

    geno_vcf = settings.work_dir / f"{pgs_id}.geno.vcf.gz"
    ancestry.markers_genotypes(bed, geno_vcf)

    # Index sample genotypes by (chrom, pos).
    calls: dict[tuple[str, int], tuple] = {}
    vcf = VCF(str(geno_vcf))
    for rec in vcf:
        calls[(rec.CHROM, rec.POS)] = (rec.REF, rec.ALT, rec.genotypes[0])
    vcf.close()

    raw = 0.0
    used: dict[tuple[str, int], tuple[str, float]] = {}
    counts: dict[tuple[str, int], int] = {}
    for v in variants:
        hit = calls.get((v.chrom, v.pos))
        if hit is None:
            continue
        ref, alts, gt = hit
        n = _effect_allele_count(ref, alts, gt, v.effect_allele, v.other_allele)
        if n is None:
            continue
        raw += n * v.weight
        # Key by the panel's non-chr naming; keep effect allele + weight to score the panel identically.
        key = (v.chrom.replace("chr", ""), v.pos)
        used[key] = (v.effect_allele, v.weight)
        counts[key] = n

    return SampleScore(
        pgs_id=pgs_id, raw=raw, n_total=len(variants), n_used=len(used),
        coverage=(len(used) / len(variants) if variants else 0.0), used=used, counts=counts,
    )


def _ref_labels() -> tuple[list[str], dict[str, str]]:
    """Reference panel sample IIDs (psam order) + IID->superpop map."""
    from . import ancestry

    psam = ancestry._panel()[2]
    iids, superpop = [], {}
    with open(psam) as fh:
        header = next(fh).lstrip("#").split()
        si = header.index("SuperPop")
        for line in fh:
            f = line.split()
            iids.append(f[0])
            superpop[f[0]] = f[si]
    return iids, superpop


def score_reference(score: SampleScore) -> tuple[np.ndarray, set[tuple[str, int]]]:
    """Score the reference panel on the sample's variants (PLINK2 --score).

    Returns the per-reference-sample sums AND the set of variants PLINK actually scored, which can
    be smaller than ``score.used`` — callers must compare the sample over that same set.
    """
    from . import ancestry, shell

    pgen, pvar, psam = ancestry._panel()
    wd = settings.work_dir
    sf = wd / f"{score.pgs_id}.refscore.txt"
    extract = wd / f"{score.pgs_id}.refscore.ids"
    sf.write_text("ID\tA1\tW\n" + "".join(f"{c}:{p}\t{ea}\t{w}\n" for (c, p), (ea, w) in score.used.items()))
    extract.write_text("".join(f"{c}:{p}\n" for (c, p) in score.used))

    out = wd / f"{score.pgs_id}.refscore"
    vars_file = Path(str(out) + ".sscore.vars")
    vars_file.unlink(missing_ok=True)
    shell.run([
        str(ancestry.plink2()), "--pgen", str(pgen), "--pvar", str(pvar), "--psam", str(psam),
        # Biallelic SNPs + dedup so chrom:pos IDs are unique (multiallelics collide otherwise).
        "--max-alleles", "2", "--snps-only", "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all",
        "--extract", str(extract),
        "--score", str(sf), "1", "2", "3", "header-read", "cols=+scoresums", "list-variants",
        "--out", str(out),
    ])
    # .sscore: per-sample; SCORE1_SUM column = sum(dosage*weight).
    sscore = {}
    with open(str(out) + ".sscore") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = next((i for i, h in enumerate(header) if h.endswith("SCORE1_SUM") or h == "SCORE1_SUM"), -1)
        iid_col = header.index("IID")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sscore[f[iid_col]] = float(f[col])
    iids, _ = _ref_labels()
    return np.array([sscore.get(i, np.nan) for i in iids]), _scored_variants(vars_file, score)


def _scored_variants(vars_file: Path, score: SampleScore) -> set[tuple[str, int]]:
    """Parse PLINK's ``list-variants`` output (``chrom:pos`` per line) into ``used`` keys."""
    if not vars_file.exists():
        raise RuntimeError(f"PLINK did not write {vars_file.name}; can't align sample and panel variant sets")
    kept: set[tuple[str, int]] = set()
    for line in vars_file.read_text().splitlines():
        c, _, p = line.strip().rpartition(":")
        if c and p.isdigit() and (c, int(p)) in score.used:
            kept.add((c, int(p)))
    return kept


def compare_to_reference(score: SampleScore, superpop: str | None) -> RefComparison | None:
    """Percentile + z of the sample within its ancestry-matched panel, over ONE shared variant set."""
    if not (superpop and score.used):
        return None
    ref, kept = score_reference(score)
    iids, sp = _ref_labels()
    mask = np.array([sp.get(i) == superpop for i in iids]) & ~np.isnan(ref)
    if mask.sum() < 20 or not kept:
        return None
    raw = score.raw_over(kept)
    dist = ref[mask]
    sd = float(dist.std())
    return RefComparison(
        raw=raw, n_used=len(kept),
        percentile=round(float((dist < raw).mean() * 100), 1),
        z=round(float((raw - dist.mean()) / sd), 2) if sd > 0 else None,
    )


def calibrate(score: SampleScore, trait: str, nearest_superpop: str | None) -> CalibratedScore:
    """Percentile of the sample's raw score within its ancestry-matched reference distribution."""
    cmp = compare_to_reference(score, nearest_superpop)
    raw, n_used = (cmp.raw, cmp.n_used) if cmp else (score.raw, score.n_used)
    return CalibratedScore(
        pgs_id=score.pgs_id, trait=trait, raw=round(raw, 4), percentile=cmp.percentile if cmp else None,
        ancestry=nearest_superpop, n_used=n_used,
        coverage=round(n_used / score.n_total, 3) if score.n_total else 0.0,
    )


def download_score(pgs_id: str) -> Path:
    """Fetch a PGS Catalog harmonized GRCh38 scoring file (cached)."""
    import httpx

    d = settings.annotations_dir / "ancestry" / "pgs"
    d.mkdir(parents=True, exist_ok=True)
    dest = d / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    if dest.exists():
        return dest
    meta = httpx.get(f"https://www.pgscatalog.org/rest/score/{pgs_id}", timeout=60).json()
    url = meta["ftp_harmonized_scoring_files"]["GRCh38"]["positions"]
    with httpx.stream("GET", url, timeout=120, follow_redirects=True) as r:
        r.raise_for_status()
        with open(dest, "wb") as fh:
            for chunk in r.iter_bytes():
                fh.write(chunk)
    return dest


def run(nearest_superpop: str | None = None) -> list[CalibratedScore]:
    """Compute the curated polygenic scores, ancestry-calibrated when a superpop is given."""
    console.rule("[bold]Polygenic scores")
    results = []
    for pgs_id, trait in CURATED_PGS:
        try:
            sf = download_score(pgs_id)
            ss = score_sample(sf, pgs_id)
            cal = calibrate(ss, trait, nearest_superpop)
            pct = f"{cal.percentile:.0f}th pct ({cal.ancestry})" if cal.percentile is not None else "raw only"
            console.print(f"  {trait:28} {pct:22} coverage {cal.coverage:.0%}")
            results.append(cal)
        except Exception as e:  # noqa: BLE001 - report and continue
            console.print(f"  [yellow]{trait}: {e}[/]")
    return results

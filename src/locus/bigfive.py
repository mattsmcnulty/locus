"""Big Five personality polygenic indices — a SEPARATE spot from health polygenic risk.

Scores the five Big Five traits (extraversion, agreeableness, conscientiousness, neuroticism,
openness) with the polygenic indices (PGIs) published alongside Schwaba et al., *Nature* 2026
("Robust inference and correlates from genetic associations with personality"; 1,260 lead variants,
N up to 1.1M). The weights are the paper's own SBayesR/HapMap3 EUR PGIs, public on its OSF project,
and are scored with the same engine as the curated health scores (``pgs.score_sample`` + the
ancestry-matched 1000 Genomes percentile).

Why its own module, table and tool rather than a row in ``pgs_scores``: a health polygenic score reads
as clinical liability, whereas a personality PGI is a weak behavioural *tilt* — these indices explain
~1–4% of trait variance — and the honest caveat differs in kind. Keeping them apart means a reader never
weighs "99th percentile extraversion" like "99th percentile coronary risk".

The weights ship on GRCh37; ``convert`` lifts them to GRCh38 by coordinate (UCSC chain via pyliftover)
into Locus's ``hmPOS_GRCh38`` scoring format. That is sufficient because ``pgs.score_reference`` keys the
reference panel on GRCh38 chrom:pos, so rsIDs are never needed. Validated on the first run: ≥99.95% of
SNPs lift, the lifted neuroticism PGI correlates r≈0.67 with independent PGS Catalog neuroticism scores
across the EUR panel, and the five PGIs reproduce the known Big Five cross-trait genetic correlations.

Privacy: the weight files and the chain are public downloads; scoring is entirely local.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console

from . import download, pgs, shell

console = Console()

PAPER = "Schwaba et al., Nature 2026 — doi 10.1038/s41586-026-10992-9 (weights: OSF project HGNSM)"

CAVEAT = (
    "A personality polygenic index is a tilt, not a prediction and not a personality test: these indices "
    "explain only ~1–4% of trait variance (extraversion the most, agreeableness the least), so even an "
    "extreme percentile shifts the expected trait by a fraction of a standard deviation — see "
    "expected_shift_sd. Percentiles are relative to the ancestry-matched 1000 Genomes reference. The "
    "paper's within-family analyses show these indices capture direct genetic effects, minimally "
    "confounded by family environment. Kept separate from polygenic_risk on purpose."
)

# The paper's EUR validation betas (SBayesR, HapMap3, population level; averaged across its five
# validation cohorts). expected_shift_sd = z × beta: what a given z on the index means in trait SDs.
VALIDATION_BETA = {"ext": 0.19, "neu": 0.17, "ope": 0.17, "con": 0.16, "agr": 0.10}

# Below this lift rate something is wrong (wrong chain, wrong file) and the score would silently be
# computed on a fraction of its variants — refuse rather than report a confident-looking number.
MIN_LIFT_RATE = 0.95

_COMP = str.maketrans("ACGT", "TGCA")


@dataclass
class BigFiveScore:
    code: str            # ext | agr | con | neu | ope
    trait: str
    raw: float
    percentile: float | None   # within the ancestry-matched reference (0-100)
    z: float | None            # standardized against that same reference distribution
    ancestry: str | None
    n_used: int
    coverage: float
    beta: float                # the paper's validation beta for this trait

    @property
    def expected_shift_sd(self) -> float | None:
        return None if self.z is None else round(self.z * self.beta, 3)


# ── GRCh37 → GRCh38 conversion ───────────────────────────────────────────────────
def parse_weights(path: Path) -> Iterator[tuple[str, int, str, str, str]]:
    """Yield (chrom, pos, effect_allele, other_allele, weight) from an SBayesR weights file
    (space-separated; columns Chrom, Position, A1 = effect allele, A2, A1Effect)."""
    with open(path) as fh:
        header = fh.readline().split()
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.split()
            if not f:
                continue
            yield f[ci["Chrom"]], int(f[ci["Position"]]), f[ci["A1"]], f[ci["A2"]], f[ci["A1Effect"]]


def lift_rows(rows: Iterable[tuple[str, int, str, str, str]], lifter) -> tuple[list, dict]:
    """Pure: lift GRCh37 rows to GRCh38 via ``lifter.convert_coordinate`` (pyliftover semantics:
    0-based in, list of (chrom, 0-based pos, strand, score) out). Complements both alleles on a
    minus-strand lift; drops lifts that fail or land on another chromosome. Returns (rows, stats)."""
    out: list[tuple[str, int, str, str, str]] = []
    st = {"in": 0, "lifted": 0, "fail": 0, "cross_chrom": 0, "minus_strand": 0}
    for chrom, pos, a1, a2, w in rows:
        st["in"] += 1
        hit = lifter.convert_coordinate(f"chr{chrom}", pos - 1)
        if not hit:
            st["fail"] += 1
            continue
        nchr, npos0, strand, _ = hit[0]
        if nchr != f"chr{chrom}":
            st["cross_chrom"] += 1
            continue
        if strand == "-":
            st["minus_strand"] += 1
            a1, a2 = a1.translate(_COMP), a2.translate(_COMP)
        out.append((chrom, npos0 + 1, a1, a2, w))
        st["lifted"] += 1
    return out, st


def scoring_file(code: str) -> Path:
    return download.bigfive_dir() / f"ReGPC_{code}_hmPOS_GRCh38.txt.gz"


def convert(code: str, *, force: bool = False) -> Path:
    """Lift one trait's GRCh37 weights into a Locus ``hmPOS_GRCh38`` scoring file (cached)."""
    dst = scoring_file(code)
    if dst.exists() and not force:
        return dst
    from pyliftover import LiftOver

    src = download.bigfive_weights_path(code)
    lifter = LiftOver(str(download.liftover_chain_path()))
    rows, st = lift_rows(parse_weights(src), lifter)
    rate = st["lifted"] / st["in"] if st["in"] else 0.0
    if rate < MIN_LIFT_RATE:
        raise shell.ToolError(
            f"Big Five {code}: only {rate:.1%} of {st['in']:,} variants lifted GRCh37→GRCh38 "
            f"(fail={st['fail']:,}, cross-chrom={st['cross_chrom']:,}). Wrong chain or weights file? "
            f"Re-run `locus download bigfive` — refusing to score a partial index.")
    label = download.BIGFIVE_WEIGHTS[code][0]
    with gzip.open(dst, "wt") as out:
        out.write(f"#pgs_id=ReGPC_{code}\n#trait_reported={label} (Big Five; {PAPER})\n"
                  "#genome_build=GRCh38 (lifted from GRCh37 via UCSC hg19ToHg38 chain)\n"
                  "#weight_type=beta (trait SD per effect allele)\n"
                  "hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\n")
        for chrom, pos, a1, a2, w in rows:
            out.write(f"{chrom}\t{pos}\t{a1}\t{a2}\t{w}\n")
    console.print(f"  {label:24} lifted {st['lifted']:,}/{st['in']:,} ({rate:.2%}; "
                  f"{st['minus_strand']} minus-strand complemented)")
    return dst


# ── Scoring ───────────────────────────────────────────────────────────────────────
def _calibrate(ss: pgs.SampleScore, superpop: str | None) -> tuple[float | None, float | None]:
    """Percentile AND z of the sample within the ancestry-matched reference (one PLINK2 pass)."""
    if not (superpop and ss.used):
        return None, None
    ref = pgs.score_reference(ss)
    iids, sp = pgs._ref_labels()
    mask = np.array([sp.get(i) == superpop for i in iids]) & ~np.isnan(ref)
    if mask.sum() < 20:
        return None, None
    dist = ref[mask]
    pct = round(float((dist < ss.raw).mean() * 100), 1)
    sd = float(dist.std())
    z = round(float((ss.raw - dist.mean()) / sd), 2) if sd > 0 else None
    return pct, z


def stored_superpop() -> str | None:
    """The nearest continental ancestry `locus ancestry` recorded for this genome (None if not run)."""
    from .db import connect, db_exists

    if not db_exists():
        return None
    try:
        with connect(read_only=True) as con:
            row = con.execute('SELECT "group" FROM ancestry_pca WHERE is_sample LIMIT 1').fetchone()
        return row[0] if row else None
    except Exception:  # noqa: BLE001 - table absent = ancestry not run
        return None


def weights_present() -> bool:
    return download.liftover_chain_path().exists() and all(
        download.bigfive_weights_path(c).exists() for c in download.BIGFIVE_WEIGHTS)


def run(nearest_superpop: str | None = None) -> list[BigFiveScore]:
    """Score the five Big Five PGIs, ancestry-calibrated. Self-skips (returns []) without weights."""
    console.rule("[bold]Big Five personality polygenic indices")
    if not weights_present():
        console.print("[yellow]Big Five weights not found — run `locus download bigfive` first. Skipping.[/]")
        return []
    superpop = nearest_superpop or stored_superpop()
    if superpop is None:
        console.print("[yellow]No ancestry on record (run `locus ancestry`) — reporting raw scores only.[/]")
    results: list[BigFiveScore] = []
    for code, (label, _url) in download.BIGFIVE_WEIGHTS.items():
        try:
            sf = convert(code)
            ss = pgs.score_sample(sf, f"ReGPC_{code}")
            pct, z = _calibrate(ss, superpop)
            s = BigFiveScore(code=code, trait=label, raw=round(ss.raw, 4), percentile=pct, z=z,
                             ancestry=superpop if pct is not None else None, n_used=ss.n_used,
                             coverage=round(ss.coverage, 3), beta=VALIDATION_BETA[code])
            shown = f"{pct:.0f}th pct ({superpop}), z {z:+.2f}" if pct is not None else "raw only"
            console.print(f"  {label:24} {shown:30} coverage {s.coverage:.0%}")
            results.append(s)
        except Exception as e:  # noqa: BLE001 - report and continue, like pgs.run
            console.print(f"  [yellow]{label}: {e}[/]")
    return results

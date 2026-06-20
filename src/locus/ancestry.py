"""Ancestry: PC projection + global-ancestry assignment, and the shared genotype
primitive used by polygenic scoring.

Pipeline (all native, no Docker/ADMIXTURE):
  1. build_reference_model() — one-time: LD-prune the 1000 Genomes panel, run
     PLINK2 PCA (with allele weights + freqs), cache reference PCs + population labels.
  2. harmonize_sample() — get the sample's genotypes (incl. hom-ref 0/0) at the pruned
     SNPs, reconciled to the panel's REF/ALT, so PLINK2 can project consistently.
  3. project_sample() — PLINK2 --score (variance-standardize) projects the sample onto
     the reference PCs.
  4. assign_ancestry() — k-NN over reference PCs → global ancestry proportions + the
     nearest superpopulation (used to pick the ancestry-matched PRS reference).

The marker primitive markers_genotypes() is also reused by polygenic scoring (pgs.py).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rich.console import Console

from . import artifacts, shell
from .config import settings
from .vcfutils import chr_rename_map, read_info, write_rename_file

console = Console()

# Superpopulation labels and display names (1000 Genomes).
SUPERPOPS = {
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
}

# Fine population display names (1000 Genomes 26 populations + HGDP groups if present).
POPULATIONS = {
    # EUR
    "GBR": "British (England/Scotland)", "CEU": "Northwest European", "IBS": "Iberian (Spain)",
    "TSI": "Tuscan (Italy)", "FIN": "Finnish",
    # SAS
    "GIH": "Gujarati", "PJL": "Punjabi", "BEB": "Bengali", "STU": "Sri Lankan Tamil", "ITU": "Telugu",
    # EAS
    "CHB": "Han Chinese (Beijing)", "CHS": "Han Chinese (South)", "JPT": "Japanese",
    "CDX": "Dai Chinese", "KHV": "Kinh (Vietnam)",
    # AFR
    "YRI": "Yoruba (Nigeria)", "LWK": "Luhya (Kenya)", "GWD": "Gambian", "MSL": "Mende (Sierra Leone)",
    "ESN": "Esan (Nigeria)", "ASW": "African-American (SW US)", "ACB": "African-Caribbean (Barbados)",
    # AMR
    "MXL": "Mexican-American", "PUR": "Puerto Rican", "CLM": "Colombian", "PEL": "Peruvian",
    # HGDP European groups (projected in for finer within-Europe resolution)
    "French": "French", "Orcadian": "Orcadian (Orkney, Scotland)", "Basque": "Basque (Spain/France)",
    "Sardinian": "Sardinian", "BergamoItalian": "Bergamo Italian", "Tuscan": "Tuscan (HGDP)",
    "Russian": "Russian", "Adygei": "Adygei (Caucasus)",
}

# HGDP region -> 1000G-style superpopulation (for the continental rollup).
HGDP_REGION_SUPERPOP = {
    "EUROPE": "EUR", "EAST_ASIA": "EAS", "AFRICA": "AFR", "AMERICA": "AMR",
    "CENTRAL_SOUTH_ASIA": "SAS", "MIDDLE_EAST": "MID", "OCEANIA": "OCE",
}
N_PCS = 10
KNN_K = 50  # nearest reference neighbours for the fine-population breakdown


def _ancestry_dir() -> Path:
    return settings.annotations_dir / "ancestry"


def _model_dir() -> Path:
    return _ancestry_dir() / "model"


def plink2() -> Path:
    p = settings.data_dir / "tools" / "plink2"
    if not p.exists():
        raise FileNotFoundError("PLINK2 not installed. Run `locus download ancestry`.")
    return p


def _panel() -> tuple[Path, Path, Path]:
    d = _ancestry_dir()
    pgen, pvar, psam = d / "all_hg38.pgen", d / "all_hg38.pvar.zst", d / "hg38_corrected.psam"
    if not (pgen.exists() and pvar.exists() and psam.exists()):
        raise FileNotFoundError("Reference panel missing. Run `locus download ancestry`.")
    return pgen, pvar, psam


# ── Shared primitive ────────────────────────────────────────────────────────────
def markers_genotypes(regions_bed: Path, dest: Path) -> Path:
    """Sample genotypes (incl. hom-ref 0/0) at the chr-prefixed positions in ``regions_bed``."""
    inputs = artifacts.classify_inputs(settings.genome_dir)
    if not inputs.small_variants:
        raise FileNotFoundError("No small-variant gVCF found in the genome dir.")
    reference = artifacts.find_reference()
    if reference is None:
        raise FileNotFoundError("Reference FASTA needed. Run `locus download reference`.")

    info = read_info(inputs.small_variants)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if info.chr_prefixed:
        region_src = f"bcftools view -R {regions_bed} {inputs.small_variants}"
    else:
        nochr = settings.work_dir / f"{dest.stem}.markers.nochr.bed"
        shell.sh(f"sed -e 's/^chrM\\t/MT\\t/' -e 's/^chr//' {regions_bed} > {nochr}")
        rename = settings.work_dir / "contigs2chr.txt"
        if not rename.exists():
            write_rename_file(chr_rename_map(info.contigs), rename)
        region_src = f"bcftools view -R {nochr} {inputs.small_variants} | bcftools annotate --rename-chrs {rename}"

    # gvcf2vcf expands whole hom-ref blocks; clip back to exactly the marker positions and
    # sort (block-boundary expansion can emit slightly out-of-order records at scale).
    shell.sh(
        f"{region_src} | bcftools convert --gvcf2vcf -f {reference} -Ou "
        f"| bcftools view -T {regions_bed} -Ou | bcftools sort -Oz -o {dest}"
    )
    shell.run(["bcftools", "index", "-f", "-t", str(dest)])
    return dest


# ── Reference model (one-time build) ────────────────────────────────────────────
def build_reference_model(maf: float = 0.05) -> Path:
    """LD-prune the panel, run PCA (allele weights + freqs), cache PCs + labels.

    Variant IDs are set to ``@:#`` (chrom:pos, panel's non-chr naming) so the sample
    can be matched by position regardless of rsID availability.
    """
    md = _model_dir()
    if (md / "pca.eigenvec.allele").exists() and (md / "ref_pcs.npz").exists():
        return md
    md.mkdir(parents=True, exist_ok=True)
    pgen, pvar, psam = _panel()
    pk = str(plink2())
    base = ["--pgen", str(pgen), "--pvar", str(pvar), "--psam", str(psam)]

    # 1. Clean to common biallelic autosomal SNPs with chrom:pos IDs (deduped). This makes
    #    every later step (prune, PCA, sample matching) key on position, not rsID.
    console.print("Building ancestry model — cleaning panel to common SNPs (chrom:pos IDs)…")
    shell.run([pk, *base, "--autosome", "--maf", str(maf), "--snps-only", "--max-alleles", "2",
               "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all", "--new-id-max-allele-len", "100",
               "--make-pgen", "--out", str(md / "clean")])
    # 2. LD-prune.
    console.print("LD-pruning…")
    shell.run([pk, "--pfile", str(md / "clean"), "--indep-pairwise", "1000", "50", "0.2",
               "--out", str(md / "prune")])
    # 3. Pruned panel + PCA (allele weights + freqs for projection).
    console.print("Extracting pruned panel…")
    shell.run([pk, "--pfile", str(md / "clean"), "--extract", str(md / "prune.prune.in"),
               "--make-pgen", "--out", str(md / "panel_pruned")])
    console.print(f"Running PCA ({N_PCS} PCs)…")
    shell.run([pk, "--pfile", str(md / "panel_pruned"), "--freq",
               "--pca", str(N_PCS), "allele-wts", "--out", str(md / "pca")])
    # Project the reference samples onto their own PCs via --score, so reference and
    # sample PCs are computed identically (same scale) for distance comparison.
    console.print("Projecting reference samples…")
    shell.run([pk, "--pfile", str(md / "panel_pruned"), "--read-freq", str(md / "pca.afreq"),
               "--score", str(md / "pca.eigenvec.allele"), "2", "5",
               "header-read", "variance-standardize", "--score-col-nums", f"6-{6 + N_PCS - 1}",
               "--out", str(md / "refproj")])

    _cache_reference_pcs(md)
    console.print(f"[green]Ancestry model ready[/] → {md}")
    return md


def _parse_sscore_pcs(path: Path) -> tuple[list[str], np.ndarray, list[str], list[str]]:
    """Parse a PLINK2 .sscore — (IIDs, PC matrix from PC*_AVG cols, superpops, fine populations)."""
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        iidx = header.index("IID")
        avg = [header.index(f"PC{k}_AVG") for k in range(1, N_PCS + 1)]
        sp_idx = header.index("SuperPop") if "SuperPop" in header else None
        pop_idx = header.index("Population") if "Population" in header else None
        iids, pcs, sps, pops = [], [], [], []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            iids.append(c[iidx])
            pcs.append([float(c[i]) for i in avg])
            sps.append(c[sp_idx] if sp_idx is not None else "NA")
            pops.append(c[pop_idx] if pop_idx is not None else "NA")
    return iids, np.array(pcs, dtype=float), sps, pops


def _cache_reference_pcs(md: Path) -> None:
    """Cache reference projected PCs + superpopulation + fine-population labels (from refproj sscore)."""
    iids, pcs, sps, pops = _parse_sscore_pcs(md / "refproj.sscore")
    np.savez(md / "ref_pcs.npz", iids=np.array(iids), pcs=pcs,
             superpop=np.array(sps), population=np.array(pops))


def add_hgdp_reference() -> int:
    """Project HGDP samples onto the 1000G PC space and append them to the cached reference.

    No panel merge: HGDP is projected with the same allele weights as the 1000G reference, then
    appended to ref_pcs.npz so the k-NN sees finer European populations (French, Orcadian, Basque…).
    Idempotent (guarded by a marker file).
    """
    md = _model_dir()
    if (md / "hgdp_added").exists():
        return 0
    hgdp = _ancestry_dir() / "hgdp"
    pgen, pvar, psam = hgdp / "hgdp_all.pgen", hgdp / "hgdp_all.pvar.zst", hgdp / "hgdp.psam"
    if not pgen.exists():
        console.print("[yellow]HGDP not present — skipping (download hgdp_all.pgen first).[/]")
        return 0
    out = settings.work_dir / "hgdp_proj"
    console.print("Projecting HGDP onto the 1000G PC space…")
    shell.run([str(plink2()), "--pgen", str(pgen), "--pvar", str(pvar), "--psam", str(psam),
               "--max-alleles", "2", "--snps-only", "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all",
               "--extract", str(md / "prune.prune.in"), "--read-freq", str(md / "pca.afreq"),
               "--score", str(md / "pca.eigenvec.allele"), "2", "5", "header-read", "variance-standardize",
               "--score-col-nums", f"6-{6 + N_PCS - 1}", "--out", str(out)])
    iids, pcs, _, _ = _parse_sscore_pcs(Path(str(out) + ".sscore"))
    # Join population + region from the HGDP psam (IID is column 0).
    pop_of, region_of = {}, {}
    with open(psam) as fh:
        h = next(fh).lstrip("#").rstrip("\n").split("\t")
        ip, ir = h.index("population"), h.index("region")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            pop_of[c[0]], region_of[c[0]] = c[ip], c[ir]
    pops = [pop_of.get(i, "NA") for i in iids]
    sps = [HGDP_REGION_SUPERPOP.get(region_of.get(i, ""), "NA") for i in iids]

    data = np.load(md / "ref_pcs.npz", allow_pickle=True)
    np.savez(md / "ref_pcs.npz",
             iids=np.concatenate([data["iids"], np.array(iids)]),
             pcs=np.concatenate([data["pcs"], pcs]),
             superpop=np.concatenate([data["superpop"], np.array(sps)]),
             population=np.concatenate([data["population"], np.array(pops)]))
    (md / "hgdp_added").write_text("1")
    console.print(f"[green]Added {len(iids)} HGDP reference samples[/] (finer European resolution).")
    return len(iids)


def pruned_sites() -> list[tuple[str, int, str, str]]:
    """(chrom[non-chr], pos, ref, alt) for the pruned panel SNPs."""
    pvar = _model_dir() / "panel_pruned.pvar"
    sites = []
    with open(pvar) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            # CHROM POS ID REF ALT
            sites.append((c[0], int(c[1]), c[3], c[4].rstrip("\n")))
    return sites


# ── Sample harmonization + projection ───────────────────────────────────────────
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def _gt_vs_panel(sref: str, salts: list[str], gt, pref: str, palt: str) -> str:
    """Sample genotype expressed against panel REF/ALT — copies of panel ALT (strand-aware)."""
    a, b = gt[0], gt[1]
    if a < 0 or b < 0:
        return "./."
    alleles = [sref, *salts]
    called = [alleles[i] if 0 <= i < len(alleles) else "N" for i in (a, b)]

    def count(alt: str, ref: str) -> str | None:
        if all(c in (alt, ref) for c in called):
            return f"{sum(1 for c in called if c == alt)}"
        return None

    for pr, pa in ((pref, palt), (_COMPLEMENT.get(pref, pref), _COMPLEMENT.get(palt, palt))):
        n = count(pa, pr)
        if n is not None:
            return {"0": "0/0", "1": "0/1", "2": "1/1"}[n]
    return "./."


def harmonize_sample(dest: Path) -> Path:
    """Sample genotypes at the pruned panel SNPs, written against the panel's REF/ALT."""
    from cyvcf2 import VCF

    sites = pruned_sites()
    bed = settings.work_dir / "ancestry.markers.bed"
    bed.write_text("".join(f"chr{c}\t{p - 1}\t{p}\n" for c, p, _, _ in sites))

    geno = settings.work_dir / "ancestry.geno.vcf.gz"
    markers_genotypes(bed, geno)
    calls: dict[tuple[str, int], tuple] = {}
    vcf = VCF(str(geno))
    for rec in vcf:
        calls[(rec.CHROM.replace("chr", ""), rec.POS)] = (rec.REF, rec.ALT, rec.genotypes[0])
    vcf.close()

    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    for c, p, ref, alt in sites:
        hit = calls.get((c, p))
        gt = _gt_vs_panel(hit[0], hit[1], hit[2], ref, alt) if hit else "./."
        lines.append(f"{c}\t{p}\t{c}:{p}\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}")
    raw = settings.work_dir / "ancestry.harmonized.vcf"
    raw.write_text("\n".join(lines) + "\n")
    shell.sh(f"bgzip -f {raw}")
    return Path(str(raw) + ".gz")


def project_sample(harmonized_vcf: Path) -> np.ndarray:
    """Project the harmonized sample onto the reference PCs via PLINK2 --score (same recipe as the ref)."""
    md = _model_dir()
    pk = str(plink2())
    out = settings.work_dir / "sample_proj"
    shell.run([pk, "--vcf", str(harmonized_vcf), "--set-all-var-ids", "@:#",
               "--read-freq", str(md / "pca.afreq"),
               "--score", str(md / "pca.eigenvec.allele"), "2", "5",
               "header-read", "variance-standardize",
               "--score-col-nums", f"6-{6 + N_PCS - 1}", "--out", str(out)])
    _, pcs, _, _ = _parse_sscore_pcs(Path(str(out) + ".sscore"))
    return pcs[0]


# ── Ancestry assignment ─────────────────────────────────────────────────────────
@dataclass
class AncestryResult:
    proportions: dict[str, float]       # superpop -> fraction (continental rollup, sums to 1)
    populations: dict[str, float]       # fine population code -> fraction (sub-continental, sums to 1)
    nearest: str                        # nearest superpopulation code (for PRS calibration)
    nearest_population: str             # nearest fine population code
    sample_pcs: list[float]
    ref_centroids: dict[str, list[float]]  # superpop -> [PC1, PC2] for plotting
    pop_points: list[list]                 # [population, PC1, PC2, superpop] per fine population


def assign_ancestry(sample_pcs: np.ndarray, k: int = KNN_K) -> AncestryResult:
    """Soft ancestry from the k nearest reference neighbours, at fine-population and continental levels.

    Fine-population fractions are a k-NN estimate (sensitive to reference panel sizes — read as
    'genetically closest to', not a calibrated admixture %). The continental rollup is robust.
    """
    md = _model_dir()
    data = np.load(md / "ref_pcs.npz", allow_pickle=True)
    ref_pcs, ref_sp, ref_pop = data["pcs"], data["superpop"], data["population"]
    d = np.linalg.norm(ref_pcs - sample_pcs[None, :], axis=1)
    nn = np.argsort(d)[:k]

    populations: dict[str, float] = {}
    for pop, c in zip(*np.unique(ref_pop[nn], return_counts=True), strict=False):
        populations[str(pop)] = round(int(c) / k, 4)
    populations = dict(sorted(populations.items(), key=lambda x: -x[1]))

    proportions = {sp: 0.0 for sp in SUPERPOPS}
    for sp, c in zip(*np.unique(ref_sp[nn], return_counts=True), strict=False):
        if sp in proportions:
            proportions[sp] = round(int(c) / k, 4)

    nearest = max(proportions, key=proportions.get)
    nearest_pop = next(iter(populations), "NA")
    centroids = {}
    for sp in SUPERPOPS:
        m = ref_sp == sp
        if m.any():
            centroids[sp] = [round(float(ref_pcs[m, 0].mean()), 4), round(float(ref_pcs[m, 1].mean()), 4)]
    # Per-fine-population centroids — the dense, colored "world populations" cloud the SPA plots.
    pop_points: list[list] = []
    for pop in sorted({str(p) for p in ref_pop}):
        m = ref_pop == pop
        sp = str(ref_sp[m][0]) if m.any() else "NA"
        pop_points.append([pop, round(float(ref_pcs[m, 0].mean()), 4),
                           round(float(ref_pcs[m, 1].mean()), 4), sp])
    return AncestryResult(proportions, populations, nearest, nearest_pop,
                          [round(float(x), 4) for x in sample_pcs], centroids, pop_points)


def run() -> AncestryResult:
    """End-to-end ancestry: build model (cached), add HGDP if present, project the sample, assign."""
    build_reference_model()
    add_hgdp_reference()  # finer European resolution if HGDP was downloaded
    console.rule("[bold]Ancestry")
    harmonized = harmonize_sample(settings.work_dir / "ancestry.harmonized.vcf.gz")
    pcs = project_sample(harmonized)
    res = assign_ancestry(pcs)
    console.print("Continental ancestry:")
    for sp, frac in sorted(res.proportions.items(), key=lambda x: -x[1]):
        if frac > 0:
            console.print(f"  {SUPERPOPS[sp]:18} {frac:5.0%}")
    console.print("Sub-continental (closest 1000 Genomes populations):")
    for pop, frac in res.populations.items():
        if frac > 0:
            console.print(f"  {POPULATIONS.get(pop, pop):28} {frac:5.0%}")
    return res

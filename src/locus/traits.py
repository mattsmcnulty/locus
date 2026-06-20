"""Single-SNP traits / wellness + the HLA-B*57:01 screening proxy.

Genotypes the curated ``panels.TAG_SNPS`` against the genome using the same
hom-ref-aware primitives as polygenic scoring: ``ancestry.markers_genotypes`` (gVCF
expansion at the marker positions) + ``pgs._effect_allele_count`` (orientation-safe
effect-allele dosage against the genome's own REF/ALT). Results land in the ``traits``
table (preserved across a variant reload, like ancestry/PGS).
"""

from __future__ import annotations

from dataclasses import dataclass

from cyvcf2 import VCF
from rich.console import Console

from . import ancestry, panels
from .config import settings
from .pgs import _effect_allele_count

console = Console()


@dataclass
class TraitResult:
    rsid: str
    category: str          # wellness | pharmacogenomic
    trait: str
    genotype: str          # e.g. "A/G", or "—" if not callable
    dosage: int | None     # effect-allele copies (0/1/2), None if no-call/unreconciled
    effect_allele: str
    interpretation: str
    note: str


def compute() -> list[TraitResult]:
    """Genotype every tag SNP and map effect-allele dosage to its phenotype text."""
    snps = panels.TAG_SNPS
    settings.work_dir.mkdir(parents=True, exist_ok=True)
    bed = settings.work_dir / "tag_snps.bed"
    bed.write_text(panels.tag_snps_bed())
    geno = settings.work_dir / "tag_snps.geno.vcf.gz"
    ancestry.markers_genotypes(bed, geno)

    calls: dict[tuple[str, int], tuple] = {}
    vcf = VCF(str(geno))
    for rec in vcf:
        calls[(rec.CHROM, rec.POS)] = (rec.REF, rec.ALT, rec.genotypes[0])
    vcf.close()

    results: list[TraitResult] = []
    for s in snps:
        hit = calls.get((s.chrom, s.pos))
        if hit is None:
            results.append(TraitResult(s.rsid, s.category, s.trait, "—", None, s.effect_allele,
                                       "Not callable in this genome.", s.note))
            continue
        ref, alts, gt = hit
        alleles = [ref, *list(alts)]

        def _base(i: int, alleles: list = alleles) -> str:  # bind per-iteration alleles
            return alleles[i] if 0 <= i < len(alleles) else "."

        gtstr = f"{_base(gt[0])}/{_base(gt[1])}"
        dosage = _effect_allele_count(ref, alts, gt, s.effect_allele, s.other_allele)
        interp = (s.interp.get(dosage, "") if dosage is not None
                  else "Genotype could not be reconciled with the expected alleles.")
        results.append(TraitResult(s.rsid, s.category, s.trait, gtstr, dosage,
                                    s.effect_allele, interp, s.note))
    return results


def haplogroup() -> TraitResult | None:
    """mtDNA maternal-lineage haplogroup via Haplogrep2 (self-skips if the jar, Java, or
    chrM variants are absent)."""
    from . import artifacts, shell

    jar = settings.annotations_dir / "haplogrep" / "haplogrep.jar"
    if not jar.exists() or shell.resolve_java() is None:
        return None
    src = artifacts.annotated_vcf() if artifacts.annotated_vcf().exists() else artifacts.sites_vcf()
    if not src.exists():
        return None
    work = settings.work_dir
    chrm, out = work / "chrM.vcf.gz", work / "haplogroup.txt"
    try:
        shell.run(["bcftools", "view", "-r", "chrM", str(src), "-Oz", "-o", str(chrm)], quiet=True)
        shell.run(["bcftools", "index", "-f", "-t", str(chrm)], quiet=True)
        n = shell.capture(["bash", "-c", f"bcftools view -H {chrm} 2>/dev/null | wc -l"]).strip()
        if n in ("", "0"):
            return None
        shell.run(shell.java_cmd(["-jar", str(jar), "classify", "--in", str(chrm),
                                  "--format", "vcf", "--out", str(out)]), quiet=True)
    except shell.ToolError:
        return None
    lines = out.read_text().splitlines() if out.exists() else []
    if len(lines) < 2:
        return None
    f = [c.strip().strip('"') for c in lines[1].split("\t")]
    hg, qual = f[1], (f[3] if len(f) > 3 else "?")
    return TraitResult(
        rsid="mtDNA", category="maternal lineage", trait="mtDNA haplogroup",
        genotype=hg, dosage=None, effect_allele="-",
        interpretation=f"Maternal-line haplogroup {hg} (Haplogrep quality {qual}).",
        note="Deep maternal ancestry from mitochondrial DNA (Phylotree 17, rCRS).")


def run() -> list[TraitResult]:
    from . import load

    console.rule("[bold]Traits")
    results = compute()
    hg = haplogroup()
    if hg is not None:
        results.append(hg)
    load.write_traits(results)
    tags = [r for r in results if r.rsid != "mtDNA"]
    n_called = sum(1 for r in tags if r.dosage is not None)
    extra = " + mtDNA haplogroup" if hg is not None else ""
    console.print(f"[green]Traits written[/] — {n_called}/{len(tags)} tag SNP(s) callable{extra}.")
    for r in results:
        console.print(f"  {r.trait}: {r.genotype} → {r.interpretation}")
    return results

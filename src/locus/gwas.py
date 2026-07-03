"""GWAS Catalog "associations you carry".

Breadth, heavily caveated: which published genome-wide-significant risk alleles does
this genome carry? Downloads the NHGRI-EBI GWAS Catalog bulk associations, keeps lead
SNPs at p < 5e-8 (deduped per rsID×trait), genotypes them hom-ref-aware via
``ancestry.markers_genotypes``, and stores the ones whose *risk allele* is present.

These are SINGLE hits with tiny, winner's-curse-inflated effects — they are surfaced as
weak/exploratory, shown with OR/beta + p, and must never be summed or styled like a
calibrated PGS (that's what ``polygenic_risk`` is for).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from cyvcf2 import VCF
from rich.console import Console

from . import ancestry
from .config import settings
from .pgs import _COMPLEMENT
from .vcfutils import canonical_chrom

console = Console()

GENOME_WIDE_SIG = 5e-8

# 0-based column indices in the ontology-annotated "alt-full" TSV.
_C_PMID, _C_TRAIT, _C_CHR, _C_POS = 1, 7, 11, 12
_C_STRONGEST, _C_SNPS, _C_PVAL, _C_ORBETA, _C_MAPPED = 20, 21, 27, 30, 34


@dataclass
class Assoc:
    rsid: str
    chrom: str
    pos: int
    risk_allele: str
    trait: str
    mapped_trait: str
    pval: float
    or_beta: str
    pmid: str


def _risk_allele(strongest: str) -> str | None:
    """'rs7903146-T' -> 'T'; 'rs..-?' / intergenic / multi -> None."""
    if "-" not in strongest:
        return None
    a = strongest.rsplit("-", 1)[1].strip().upper()
    return a if a in ("A", "C", "G", "T") else None


def parse(tsv: Path, pmax: float = GENOME_WIDE_SIG) -> list[Assoc]:
    """Stream the catalog → genome-wide-significant single lead SNPs, deduped per
    (rsID, mapped-trait) keeping the most significant."""
    best: dict[tuple[str, str], Assoc] = {}
    with open(tsv, encoding="utf-8", errors="replace") as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= _C_MAPPED:
                continue
            rsid, chr_id, chr_pos = f[_C_SNPS].strip(), f[_C_CHR].strip(), f[_C_POS].strip()
            # single SNP only (skip haplotype / multi-SNP rows)
            if not rsid.startswith("rs") or any(c in rsid for c in " x;,") or not chr_id or not chr_pos:
                continue
            if any(c in chr_id for c in " x;,") or any(c in chr_pos for c in " x;,"):
                continue
            ra = _risk_allele(f[_C_STRONGEST])
            if not ra:
                continue
            try:
                pos, pval = int(chr_pos), float(f[_C_PVAL])
            except ValueError:
                continue
            if not (0 < pval < pmax):
                continue
            mapped = f[_C_MAPPED].strip() or f[_C_TRAIT].strip()
            key = (rsid, mapped)
            if key not in best or pval < best[key].pval:
                best[key] = Assoc(rsid, canonical_chrom(chr_id), pos, ra, f[_C_TRAIT].strip(),
                                  mapped, pval, f[_C_ORBETA].strip(), f[_C_PMID].strip())
    return list(best.values())


def _risk_dosage(ref: str, alts: list[str], gt, risk: str) -> int | None:
    """Copies of the risk allele (0/1/2). Unlike pgs._effect_allele_count this returns 0
    (not None) when the person is hom-ref and the risk allele simply isn't present —
    'carries zero risk alleles' is a real answer. Returns None only when the risk allele
    isn't represented at the site at all (can't assess)."""
    a, b = gt[0], gt[1]
    if a < 0 or b < 0:
        return None
    alleles = [ref, *list(alts)]
    called = [alleles[i] if 0 <= i < len(alleles) else "" for i in (a, b)]
    site = set(alleles)
    if risk in site:
        return sum(1 for c in called if c == risk)
    flip = _COMPLEMENT.get(risk, risk)
    if flip in site:
        return sum(1 for c in called if c == flip)
    return None  # risk allele not at this site (e.g. multiallelic mismatch) — exclude


@dataclass
class Carried:
    rsid: str
    chrom: str
    pos: int
    risk_allele: str
    dosage: int
    zygosity: str
    trait: str
    mapped_trait: str
    pval: float
    or_beta: str
    pmid: str


def compute(assocs: list[Assoc]) -> list[Carried]:
    """Genotype the lead SNPs once and keep associations whose risk allele is carried."""
    positions = sorted({(a.chrom, a.pos) for a in assocs})
    settings.work_dir.mkdir(parents=True, exist_ok=True)
    bed = settings.work_dir / "gwas.markers.bed"
    bed.write_text("".join(f"{c}\t{p - 1}\t{p}\n" for c, p in positions))
    geno = settings.work_dir / "gwas.geno.vcf.gz"
    console.print(f"Genotyping {len(positions):,} lead-SNP positions…")
    ancestry.markers_genotypes(bed, geno)

    calls: dict[tuple[str, int], tuple] = {}
    vcf = VCF(str(geno))
    for rec in vcf:
        calls[(rec.CHROM, rec.POS)] = (rec.REF, rec.ALT, rec.genotypes[0])
    vcf.close()

    carried: list[Carried] = []
    for a in assocs:
        hit = calls.get((a.chrom, a.pos))
        if hit is None:
            continue
        ref, alts, gt = hit
        d = _risk_dosage(ref, alts, gt, a.risk_allele)
        if not d:  # None (unassessable) or 0 (risk allele absent)
            continue
        carried.append(Carried(a.rsid, a.chrom, a.pos, a.risk_allele, d,
                               "homozygous" if d == 2 else "heterozygous",
                               a.trait, a.mapped_trait, a.pval, a.or_beta, a.pmid))
    return carried


def _resolve_rsids(rsids: list[str]) -> dict[str, tuple[str, int]]:
    """rsID -> (chrom, GRCh38 pos) via Ensembl (generic query; only rsIDs leave the machine)."""
    import httpx

    out: dict[str, tuple[str, int]] = {}
    try:
        r = httpx.post("https://rest.ensembl.org/variation/homo_sapiens",
                       headers={"Content-Type": "application/json", "Accept": "application/json"},
                       json={"ids": rsids}, timeout=60)
        r.raise_for_status()
        data = r.json()
    except Exception as e:  # noqa: BLE001
        console.print(f"[yellow]Ensembl lookup failed:[/] {e}")
        return out
    for rsid, info in (data or {}).items():
        for m in info.get("mappings") or []:
            sr = str(m.get("seq_region_name", ""))
            if m.get("assembly_name") == "GRCh38" and "_" not in sr:
                out[rsid] = (canonical_chrom(sr), int(m["start"]))
                break
    return out


def ask_markers(rsids: list[str]) -> list[dict]:
    """On-demand genotype lookup at arbitrary rsIDs (e.g. from a new paper), hom-ref-aware,
    joined with any local ClinVar/AlphaMissense annotations we carry."""
    from .db import connect

    rsids = [r for r in rsids if r.lower().startswith("rs")]
    pos = _resolve_rsids(rsids)
    if not pos:
        return []
    settings.work_dir.mkdir(parents=True, exist_ok=True)
    bed = settings.work_dir / "ask.markers.bed"
    bed.write_text("".join(f"{c}\t{p - 1}\t{p}\n" for c, p in pos.values()))
    geno = settings.work_dir / "ask.geno.vcf.gz"
    ancestry.markers_genotypes(bed, geno)

    calls: dict[tuple[str, int], tuple] = {}
    vcf = VCF(str(geno))
    for rec in vcf:
        calls[(rec.CHROM, rec.POS)] = (rec.REF, rec.ALT, rec.genotypes[0])
    vcf.close()

    ann: dict[str, tuple] = {}
    with connect(read_only=True) as con:
        if rsids:
            ph = ", ".join("?" for _ in rsids)
            for r in con.execute(
                f"SELECT rsid, gene, clnsig, am_class FROM variants WHERE rsid IN ({ph})", rsids
            ).fetchall():
                ann[r[0]] = (r[1], r[2], r[3])

    results = []
    for rsid, (chrom, p) in pos.items():
        hit = calls.get((chrom, p))
        if hit is None:
            gtstr, ref = "—", None
        else:
            ref, alts, gt = hit
            alleles = [ref, *list(alts)]
            gtstr = "/".join(alleles[i] if 0 <= i < len(alleles) else "." for i in gt[:2])
        gene, clnsig, am = ann.get(rsid, (None, None, None))
        results.append({"rsid": rsid, "chrom": chrom, "pos": p, "genotype": gtstr, "ref": ref,
                        "gene": gene, "clnsig": clnsig, "am_class": am})
    return results


def run() -> list[Carried]:
    from . import download, load

    console.rule("[bold]GWAS associations you carry")
    tsv = download.setup_gwas()
    console.print("Parsing GWAS Catalog (p < 5e-8 lead SNPs)…")
    assocs = parse(tsv)
    console.print(f"{len(assocs):,} genome-wide-significant lead associations after dedup.")
    carried = compute(assocs)
    load.write_associations(carried)
    n_snps = len({c.rsid for c in carried})
    console.print(f"[green]Carried[/] {len(carried):,} associations across {n_snps:,} risk SNPs "
                  f"(weak/exploratory — single hits, not a calibrated score).")
    return carried

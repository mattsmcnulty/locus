"""Curated gene/variant panels for VCF-only interpretation breadth.

Two things live here:

1. ``ACMG_SF_GENES`` — the ACMG SF v3.3 medically-actionable "secondary findings"
   genes. ``queries.secondary_findings()`` reports pathogenic / likely-pathogenic
   ClinVar variants the genome carries in these genes; "nothing found" is a defensible,
   reassuring result — which is exactly why the list must be complete and the
   recessive genes must be zygosity-gated. Both are load-bearing for that reassurance.

2. ``TAG_SNPS`` — well-characterized single-SNP traits (the 23andMe-style wellness
   layer) plus the HLA-B*57:01 screening proxy. Coordinates are GRCh38 and alleles are
   **forward-strand, verified against Ensembl** — interpretation orientation is encoded
   per SNP by *effect-allele dosage* (0/1/2 copies), and dosage is computed against the
   genome's own REF/ALT via ``pgs._effect_allele_count`` so strand never silently flips.
"""

from __future__ import annotations

from dataclasses import dataclass

# ACMG SF v3.3 (84 genes) — Lee K, et al. Genet Med. 2025;27(8):101454.
# A filter over ClinVar P/LP — confirm any hit clinically. Keep this list exact: a gene missing
# here is a screen that silently never happens, reported to the user as a reassuring "no findings".
ACMG_SF_GENES = frozenset({
    # Cancer risk (28)
    "APC", "BMPR1A", "BRCA1", "BRCA2", "MAX", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH",
    "NF2", "PALB2", "PMS2", "PTEN", "RB1", "RET", "SDHAF2", "SDHB", "SDHC", "SDHD",
    "SMAD4", "STK11", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WT1",
    # Cardiovascular disease (41)
    "ACTA2", "ACTC1", "APOB", "BAG3", "CALM1", "CALM2", "CALM3", "CASQ2", "COL3A1",
    "DES", "DSC2", "DSG2", "DSP", "FBN1", "FLNC", "KCNH2", "KCNQ1", "LDLR", "LMNA",
    "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3", "PCSK9", "PKP2", "PLN", "PRKAG2",
    "RBM20", "RYR2", "SCN5A", "SMAD3", "TGFBR1", "TGFBR2", "TMEM43", "TNNC1", "TNNI3",
    "TNNT2", "TPM1", "TRDN", "TTN",
    # Inborn errors of metabolism (6)
    "ABCD1", "BTD", "CYP27A1", "GAA", "GLA", "OTC",
    # Other genetic disease (9)
    "ACVRL1", "ATP7B", "CACNA1S", "ENG", "HFE", "HNF1A", "RPE65", "RYR1", "TTR",
})
ACMG_SF_VERSION = "v3.3"

# ACMG reports these only when TWO P/LP variants are present (biallelic) — a single heterozygous
# carrier is not a secondary finding. Without this gate, common carrier states (HFE p.C282Y is
# carried by ~10% of Europeans) surface as actionable "findings", which is a false alarm.
ACMG_SF_RECESSIVE = frozenset({
    "MUTYH", "CASQ2", "TRDN", "BTD", "CYP27A1", "GAA", "ATP7B", "HFE", "HNF1A", "RPE65",
})


# ── Carrier screening ────────────────────────────────────────────────────────────
# Recessive conditions where carrying ONE pathogenic copy is usually silent for you but matters
# for family planning (two carriers → 1-in-4 risk per pregnancy). This is the complement of
# `secondary_findings`, which deliberately drops lone heterozygous carriers.
#
# Scope, stated honestly: this is a curated panel of common, **VCF-assessable** conditions. It is
# NOT a clinical carrier screen — ACMG's Tier 3 panel is 113 genes (Gregg et al, Genet Med 2021).
# Nothing here should be read as "you are not a carrier" for anything off this list.
@dataclass(frozen=True)
class CarrierGene:
    gene: str
    condition: str
    inheritance: str   # "AR" (autosomal recessive) | "XL" (X-linked)


CARRIER_PANEL: tuple[CarrierGene, ...] = (
    CarrierGene("CFTR", "Cystic fibrosis", "AR"),
    CarrierGene("HBB", "Sickle cell disease / beta-thalassemia", "AR"),
    CarrierGene("HEXA", "Tay-Sachs disease", "AR"),
    CarrierGene("GBA", "Gaucher disease", "AR"),
    CarrierGene("PAH", "Phenylketonuria (PKU)", "AR"),
    CarrierGene("GALT", "Classic galactosemia", "AR"),
    CarrierGene("ATP7B", "Wilson disease", "AR"),
    CarrierGene("SERPINA1", "Alpha-1 antitrypsin deficiency", "AR"),
    CarrierGene("BTD", "Biotinidase deficiency", "AR"),
    CarrierGene("ACADM", "MCAD deficiency", "AR"),
    CarrierGene("ASPA", "Canavan disease", "AR"),
    CarrierGene("ELP1", "Familial dysautonomia", "AR"),
    CarrierGene("FANCC", "Fanconi anemia, group C", "AR"),
    CarrierGene("BLM", "Bloom syndrome", "AR"),
    CarrierGene("MCOLN1", "Mucolipidosis IV", "AR"),
    CarrierGene("NPC1", "Niemann-Pick disease, type C1", "AR"),
    CarrierGene("SLC26A4", "Pendred syndrome / hearing loss", "AR"),
    CarrierGene("GJB2", "Nonsyndromic hearing loss (connexin 26)", "AR"),
    CarrierGene("SLC22A5", "Primary carnitine deficiency", "AR"),
    CarrierGene("CLN3", "Batten disease (juvenile NCL)", "AR"),
    CarrierGene("USH2A", "Usher syndrome type 2A", "AR"),
    CarrierGene("PCDH15", "Usher syndrome type 1F", "AR"),
    CarrierGene("DMD", "Duchenne/Becker muscular dystrophy", "XL"),
    CarrierGene("F8", "Hemophilia A", "XL"),
    CarrierGene("F9", "Hemophilia B", "XL"),
    CarrierGene("G6PD", "G6PD deficiency", "XL"),
)

# Conditions this data CANNOT speak to. Reported explicitly as "not assessed" rather than being
# quietly absent — the whole failure mode we keep hitting is a screen that never ran being read as
# a clean result. Several of the most important carrier tests are exactly these.
CARRIER_UNASSESSABLE: tuple[tuple[str, str, str], ...] = (
    ("SMN1", "Spinal muscular atrophy",
     "carrier status is a copy-number call (SMN1 exon-7 dosage, and SMN2 paralog interference); "
     "short-read VCFs cannot determine it. This is one of the two universally-offered screens."),
    ("FMR1", "Fragile X syndrome",
     "caused by a CGG repeat expansion, which a VCF does not represent at all — needs the reads."),
    ("HBA1/HBA2", "Alpha-thalassemia",
     "most carriers have whole-gene deletions (copy-number), not SNVs."),
    ("CYP21A2", "Congenital adrenal hyperplasia (21-hydroxylase deficiency)",
     "near-identical CYP21A1P pseudogene makes short-read calls unreliable."),
)


@dataclass(frozen=True)
class TagSnp:
    rsid: str
    chrom: str            # chr-prefixed, GRCh38
    pos: int              # GRCh38 1-based
    effect_allele: str    # forward strand (Ensembl-verified)
    other_allele: str     # forward strand
    category: str         # "wellness" | "pharmacogenomic"
    trait: str
    interp: dict[int, str]  # effect-allele dosage (0/1/2) -> phenotype text
    note: str = ""
    recessive: bool = False  # phenotype expressed only with 2 effect alleles (for display)


# Forward-strand, GRCh38, Ensembl-verified (rs-id → chrom/pos/alleles checked 2026-06-20).
TAG_SNPS: list[TagSnp] = [
    TagSnp("rs4988235", "chr2", 135851076, "A", "G", "wellness", "Lactase persistence (LCT/MCM6)",
           {0: "Likely lactose intolerant as an adult (no persistence allele).",
            1: "Likely lactase-persistent — can digest lactose.",
            2: "Lactase-persistent — can digest lactose."}),
    TagSnp("rs671", "chr12", 111803962, "A", "G", "wellness", "Alcohol flush (ALDH2)",
           {0: "Normal ALDH2 — no alcohol-flush reaction.",
            1: "Reduced ALDH2 — alcohol flush and higher sensitivity.",
            2: "ALDH2-deficient — strong alcohol flush; elevated alcohol-related risk."},
           note="The flush (A) allele is rare in Europeans."),
    TagSnp("rs17822931", "chr16", 48224287, "T", "C", "wellness", "Earwax type & body odor (ABCC11)",
           {0: "Wet earwax; typical body odor.",
            1: "Wet earwax (carrier of the dry allele).",
            2: "Dry earwax; markedly reduced body odor."}, recessive=True),
    TagSnp("rs762551", "chr15", 74749576, "A", "C", "wellness", "Caffeine metabolism (CYP1A2)",
           {0: "Slow caffeine metabolizer (CYP1A2 *1A/*1A).",
            1: "Intermediate caffeine metabolizer.",
            2: "Fast caffeine metabolizer (CYP1A2 *1F/*1F)."}),
    TagSnp("rs12913832", "chr15", 28120472, "A", "G", "wellness", "Eye color (HERC2/OCA2)",
           {0: "Brown eyes likely (GG).",
            1: "Intermediate — green/hazel possible.",
            2: "Blue eyes likely (AA)."}),
    TagSnp("rs1815739", "chr11", 66560624, "T", "C", "wellness", "Muscle fiber type (ACTN3 R577X)",
           {0: "Two functional ACTN3 copies (CC) — power/sprint-associated.",
            1: "One functional ACTN3 copy (CT) — mixed.",
            2: "No functional ACTN3 (TT) — endurance-associated."}),
    TagSnp("rs2395029", "chr6", 31464003, "G", "T", "pharmacogenomic",
           "HLA-B*57:01 proxy — abacavir hypersensitivity (HCP5)",
           {0: "No HLA-B*57:01 tag — standard abacavir risk.",
            1: "Carries the HLA-B*57:01 tag — likely B*57:01 positive: abacavir hypersensitivity "
               "risk (CPIC: do not prescribe abacavir).",
            2: "Homozygous HLA-B*57:01 tag — likely B*57:01 positive: abacavir hypersensitivity risk."},
           note="Tag SNP in near-perfect LD with HLA-B*57:01 in Europeans (r²≈1) — a screening "
                "proxy, not a clinical HLA type. Confirm with HLA typing before acting."),
]


def tag_snps_bed() -> str:
    """BED text (0-based start) for all tag-SNP positions, for markers_genotypes()."""
    return "".join(f"{s.chrom}\t{s.pos - 1}\t{s.pos}\n" for s in TAG_SNPS)

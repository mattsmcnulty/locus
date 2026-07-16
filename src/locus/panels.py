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

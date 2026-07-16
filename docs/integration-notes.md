# Locus integration notes

Verified 2026-06-19 via a multi-agent research + adversarial-verification pass. These are the
external facts the pipeline depends on. **No genetic data here — safe to commit.** Re-verify URLs
and versions periodically (NCBI/Ensembl/gnomAD refresh on schedules noted below).

---

## Reference genome (GRCh38)

Use the **NCBI "no-alt analysis set", chr-prefixed** — this is the assembly DRAGEN's hg38 reference
is built from, so coordinates match the sequencing.com VCFs exactly.

- URL: `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`
- ~873 MB gzipped / ~3.1 GB uncompressed. Contigs: `chr1..chr22, chrX, chrY, chrM`, EBV, scaffolds.
- **Gotcha:** NCBI ships **plain gzip, not bgzip** → `samtools faidx` fails ("not a BGZF file").
  Must `gzip -d` then `bgzip` then `samtools faidx` (regenerate the `.fai`/`.gzi`).
- Do **not** use the `_plus_hs38d1` decoy set — decoy/ALT/HLA contigs never carry reportable variants.

```bash
gzip -d GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38_no_alt.fa
bgzip -@ 4 GRCh38_no_alt.fa            # -> GRCh38_no_alt.fa.gz
samtools faidx GRCh38_no_alt.fa.gz    # -> .fai + .gzi
```

## ⚠️ Observed sequencing.com format (2026-06, this account)

The research below assumed DRAGEN, but the **actual** sequencing.com 30× WGS file is different —
verify per file with `bcftools view -h`:

- **Caller is bcftools/mpileup, not DRAGEN** (`##source=Sequencing.com (30x WGS)`; INFO tags
  `INDEL/IDV/IMF/VDB/RPBZ/MQBZ/BQBZ/SGB/MQ0F`). No DRAGEN-specific FILTERs — just `PASS`/`FAIL`.
- **It IS a gVCF**, but hom-ref blocks are `ALT="."` + `INFO/END` rows (~half the records), **not**
  the `<NON_REF>` convention. The ingest drop-filter therefore matches `ALT="." OR ALT="<NON_REF>"`.
- **Contigs are non-chr-prefixed**: `1..22, X, Y, MT` + scaffolds (`1_KI270706v1_random`,
  `Un_KI270302v1`). Ingest renames all to `chr*` (`MT`→`chrM`) to match the chr-prefixed reference/DBs.
- **rsIDs are already populated** in the ID column → the dbSNP step is unnecessary.
- Reference is `GRCh38.p13`; the NCBI no-alt analysis set is concordant (0 REF mismatches observed).

Everything below is still the correct integration approach; `ingest` is caller-agnostic once blocks
are dropped and contigs canonicalized.

## DRAGEN VCF format (sequencing.com output)

- `*.snp-indel.genome.vcf.gz` is a **gVCF**: symbolic `<NON_REF>` allele + hom-ref blocks with `INFO/END`.
  Loading it raw massively over-counts "reference" rows. Detect via `##ALT=<ID=NON_REF` / `ALT=<NON_REF>`.
  Convert to a normalized sites VCF:
  ```bash
  bcftools view -e 'ALT="<NON_REF>"' in.genome.vcf.gz \
    | bcftools norm -f GRCh38_no_alt.fa.gz -m -any --keep-sum AD - \
    | bcftools annotate -x INFO/END -Oz -o sample.sites.vcf.gz
  tabix -p vcf sample.sites.vcf.gz
  ```
  For **PharmCAT** specifically, instead expand blocks so hom-ref PGx positions survive:
  `bcftools convert --gvcf2vcf -f <ref> ...` (see PharmCAT below).
- Contigs are **chr-prefixed**. chrM is rCRS and gets separate mito filtering (not the genome-wide QUAL filters).
- DRAGEN-only tags (`DRAGstrInfo`, `DRAGstrParams`) break strict parsers — strip with `bcftools annotate -x` if needed.
- Small-variant default FILTERs: `DRAGENSnpHardQUAL` (QUAL<10.41), `DRAGENIndelHardQUAL` (QUAL<7.83),
  `LowDepth` (DP<=1). Per-sample status in `FORMAT/FT`; germline confidence in `FORMAT/GQ`.
- `*.cnv.vcf.gz`: `SVTYPE=CNV`, ALT `<DEL>`/`<DUP>`, `INFO/END,SVLEN,REFLEN`, **`FORMAT/CN`** (int copy number;
  2 = normal diploid) + `SM` (copy-ratio). Classify on **CN**, not GT. FILTERs `cnvQual`, `cnvCopyRatio`.
- `*.sv.vcf.gz`: Manta-style. `SVTYPE` DEL/DUP/INS/INV/**BND**. BND = paired breakends linked by `INFO/MATEID`
  with bracket ALT notation — a naive POS..END BED is **wrong** for BND. `FORMAT/PR` (paired-read) + `SR`
  (split-read) support. For gene overlap prefer **AnnotSV** (handles BND, adds OMIM/ACMG); bedtools as fallback.

## ClinVar (clinical significance)

- URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` (+ `.tbi`, `.md5`). ~183 MB.
  Weekly refresh; monthly first-Thursday is the permanently archived one (pin to that for reproducibility,
  or just store our own dated copy + md5). Self-describing (embeds INFO headers → no `-h` needed).
- **Gotcha:** ClinVar uses **non-prefixed contigs (1, X, MT) with NO `##contig` lines**. Against a
  chr-prefixed DRAGEN VCF, `bcftools annotate` matches **ZERO** records silently. Rename ClinVar →
  chr-prefixed (note `MT`→`chrM`, not `chrMT`).
- Curated field set to transfer: `CLNSIG, CLNDN, CLNREVSTAT, CLNVC, GENEINFO, CLNDISDB, CLNHGVS, MC, ALLELEID`.
- Validate after: `bcftools view -H out.vcf.gz | grep -c CLNSIG` must be > 0 (0 ⇒ contig mismatch).

```bash
bcftools annotate --rename-chrs clinvar2chr.txt clinvar.vcf.gz -Oz -o clinvar.chr.vcf.gz
bcftools index -t clinvar.chr.vcf.gz
bcftools annotate -a clinvar.chr.vcf.gz \
  -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC,INFO/GENEINFO,INFO/CLNDISDB,INFO/CLNHGVS,INFO/MC,INFO/ALLELEID \
  sample.sites.vcf.gz -Oz -o sample.clinvar.vcf.gz
```

## dbSNP (rsIDs) — usually skippable

- DRAGEN typically already fills the ID column with rsIDs. Only run dbSNP if IDs are missing/`.`.
- If needed: build 157, `https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz` (~28 GB).
  CHROM uses **RefSeq accessions** (`NC_000001.11`) → rename via the GRCh38.p14 assembly_report
  (col7 RefSeq-Accn → col10 UCSC name; drop `na` rows). Gate the 28 GB download behind the "IDs missing" check.

## gnomAD v4.1 (population allele frequency)

- Latest sites VCFs: `release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr{1..22,X,Y}.vcf.bgz`
  (a v4.1.1 re-release also exists with same AFs / updated VRS+VEP annotations — either works).
- Mirrors (avoid GCP egress from outside Google Cloud — prefer AWS no-sign-request):
  - GCS: `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/<file>`
  - AWS: `s3://gnomad-public-us-east-1/release/4.1/vcf/genomes/<file>` (`--no-sign-request`, range reads OK)
- **Genomes** set (~563 GB total) is right for WGS. **Do not bulk download** — `tabix`/`bcftools view -R`
  stream only the regions overlapping the sample's variants over HTTPS range reads.
- Contigs are **chr-prefixed**. v4 fields: `AF, AC, AN, AF_grpmax, grpmax`, per-ancestry
  `AF_{afr,amr,asj,eas,fin,mid,nfe,sas}`. **`AF_popmax` was renamed `AF_grpmax`** in v4 (gone).
  Consider `bcftools view -f PASS` first (many rows are AC0/AS_VQSR). Prefix transferred fields `gnomAD_*`.

## VEP (gene + functional consequence) — default annotator

- Release 116 (June 2026). On Apple Silicon, **Docker is most reliable** (`ensemblorg/ensembl-vep`);
  bioconda `ensembl-vep` deps have no native osx-arm64 build.
- Offline cache: `https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz`
  (~22–26 GB). Reference FASTA (required for HGVS): `.../release-116/fasta/homo_sapiens/dna_index/Homo_sapiens.GRCh38.dna.toplevel.fa.gz`.
- Run `vep --offline --cache --dir_cache ... --fasta ... --vcf --everything --fork N --synonyms chr_synonyms.txt`.
  Output is a pipe-delimited `CSQ` INFO field; parse field names from the `##INFO=<ID=CSQ` header.
- **SnpEff** is the lighter fallback (pure-Java jar, runs native on arm64): `snpEff_latest_core.zip`,
  DB `GRCh38.mane.1.0.ensembl`, emits `ANN`. Pair with SnpSift for filtering.

## PharmCAT (pharmacogenomics) — Docker-first

- **v3.2.0** (Feb 2026). Requires Java 17+. **Cannot take a raw DRAGEN VCF** — requires its VCF Preprocessor,
  which rejects gVCFs.
- **Use the `pgkb/pharmcat` Docker image** (note org is `pgkb`, not `pharmgkb`): bundles JRE 17, bcftools,
  bgzip, the preprocessor, and the versioned `pharmcat_positions.vcf.bgz` — sidesteps the bare-metal
  reference-FASTA/positions-file versioning issues.
- Flow: convert the DRAGEN **gVCF → regular VCF with REF blocks expanded** (`bcftools convert --gvcf2vcf -f <ref>`)
  so hom-ref PGx positions survive (PharmCAT treats absent positions as **no-call**, not reference), then
  `pharmcat_pipeline sample.regular.vcf.gz -o results`.
- Output: parse `*.report.json` (per-gene `sourceDiplotype`, phenotype/activity; drug section with CPIC/DPWG
  recommendations) and/or the calls-only TSV. **Breaking rename:** `wildtypeAllele` → `referenceAllele` in
  recent JSON. Do **not** use `--missing-to-ref` for a real report (research-only, can manufacture calls).
- Caveat: CYP2D6 copy-number is imperfect from short-read WGS — flag clinically-actionable calls for confirmation.

## MCP Python SDK + DuckDB

- Package `mcp` 1.28.x. **Pin `mcp>=1.27,<2` in pyproject** — the upper bound is load-bearing, not
  decorative. Import FastMCP from **`mcp.server.fastmcp`** (the vendored 1.x), NOT the standalone `fastmcp` 3.x.
- **v2 migration — assessed 2026-07, deliberately deferred.** Measured against `mcp==2.0.0b2` in a throwaway
  venv rather than assumed:
  - v2 is **pre-release only** (2.0.0a1–a3, b1, b2); newest stable is 1.28.x. Migrating now means building
    against a beta whose API can still move before 2.0.0 — i.e. migrating twice.
  - `mcp.server.fastmcp` is **removed outright** in v2 (ModuleNotFoundError). FastMCP is not renamed or
    relocated: `mcp.FastMCP` and `mcp.server.FastMCP` are both absent. The whole server would need porting.
  - The replacement is **`from mcp.server import MCPServer`**, which keeps the same decorator shape —
    `.tool()`, `.prompt()`, `.resource()`, `.run()` all exist — so the port looks close to
    `FastMCP("locus")` → `MCPServer("locus")` plus a re-verify of the output schemas.
  - The **client** API our stdio test uses (`ClientSession`, `StdioServerParameters`, `mcp.client.stdio`)
    survives v2 unchanged.
  - **Trigger to revisit:** 2.0.0 final on PyPI. Then: bump the bound, swap the import/constructor, and
    re-run `test_mcp_stdio_tools_respond` — it exercises the real stdio path and asserts the object-typed
    output schemas, which is exactly what a framework swap is most likely to break.
- Tools: `@mcp.tool()` with a **return type annotation** (Pydantic model / `list[...]`) → auto output schema +
  validated structured JSON. Without an annotation you only get unstructured text.
- Claude Code truncates MCP output > ~10k tokens → **return concise, paginated results** (`limit`/`offset`,
  include `total`). `MAX_MCP_OUTPUT_TOKENS` raises the cap.
- Register: `claude mcp add --scope project --transport stdio locus -- uv run locus-mcp` (writes `.mcp.json`).
  Flags go **before** the name; everything after `--` is the server command. Claude Desktop: edit
  `~/Library/Application Support/Claude/claude_desktop_config.json` with the same `mcpServers` shape.
- **DuckDB concurrency:** never mix read-write + read-only across processes. Ingest writes with nothing else
  attached; then MCP **and** API open with `read_only=True`. Re-ingest = write a new file and swap.
- Load path: `cyvcf2` stream → `pyarrow` batches → `con.register(view, table); INSERT INTO ... SELECT`.
  Physically `ORDER BY chrom,pos` for zonemap pruning; ART indexes on `(chrom,pos)`, `rsid`, `gene` for
  selective point lookups. cyvcf2: `v.ALT` is a **list**, `v.FILTER is None` means PASS, `v.INFO.get(...)`.

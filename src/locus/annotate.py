"""Annotate phase: layer open-source databases onto the sites VCF.

Each step is independent and self-skips when its database isn't present, so you
can run a subset (``locus annotate --steps clinvar,pharmcat``) or ``all``. The
small-variant chain is successive ``bcftools annotate`` calls; PharmCAT runs
separately (it needs the gVCF, not the sites VCF) and writes a JSON report that
the load phase ingests.
"""

from __future__ import annotations

from pathlib import Path

from rich.console import Console

from . import artifacts, download, shell
from .config import settings
from .vcfutils import read_info

console = Console()

# ClinVar fields to transfer (curated; not the whole INFO block).
CLINVAR_FIELDS = (
    "INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC,INFO/CLNDISDB,"
    "INFO/CLNHGVS,INFO/MC,INFO/ALLELEID,INFO/GENEINFO"
)
# gnomAD allele frequencies come from Ensembl's REST mirror, keyed by rsID.
#
# Why not gnomAD's own VCFs: they're only published per-chromosome (chr21 alone is 7.2 GB, the
# set is ~300 GB) with no AF-only slim file. Streaming them is latency-bound — ~0.3s per position
# plus a large per-chromosome index fetch — so even a scoped ~52k-position set took hours and died
# mid-run. Ensembl answers a batch of ~180 rsIDs in seconds, keeps setup lightweight, and reuses
# the rsID-only egress posture we already have (see gwas._resolve_rsids).
ENSEMBL_VARIATION = "https://rest.ensembl.org/variation/homo_sapiens"
# 'ALL' is the global AF; 'remaining' isn't a real population, so it's excluded from grpmax.
_GNOMAD_GRPMAX_SKIP = {"ALL", "remaining"}

# Declaration order only; run() applies the real order (gnomad last — see run()).
ALL_STEPS = ["clinvar", "snpeff", "alphamissense", "gnomad", "pharmcat"]


def _index(vcf: Path) -> None:
    shell.run(["bcftools", "index", "-f", "-t", str(vcf)])


def _present_info(vcf: Path) -> set[str]:
    """INFO field IDs defined in a VCF header (so we don't request absent fields)."""
    header = shell.capture(["bcftools", "view", "-h", str(vcf)])
    return {
        ln.split("ID=", 1)[1].split(",", 1)[0]
        for ln in header.splitlines()
        if ln.startswith("##INFO=<ID=")
    }


def annotate_clinvar(src: Path, dest: Path) -> Path:
    clinvar = settings.annotations_dir / download.CLINVAR_CHR_VCF
    if not clinvar.exists():
        console.print("[yellow]ClinVar DB missing — skipping. Run `locus download clinvar`.[/]")
        return src
    # Only transfer fields actually present in this ClinVar build.
    present = _present_info(clinvar)
    fields = [f for f in CLINVAR_FIELDS.split(",") if f.split("/")[-1] in present]
    if not fields:
        console.print("[yellow]No expected ClinVar INFO fields found — skipping.[/]")
        return src
    console.print("Annotating ClinVar clinical significance…")
    shell.run([
        "bcftools", "annotate", "-a", str(clinvar), "-c", ",".join(fields),
        str(src), "-Oz", "-o", str(dest),
    ])
    _index(dest)
    # Validate the join actually matched something (0 ⇒ contig mismatch — see notes).
    n = shell.capture(["bash", "-o", "pipefail", "-c", f"bcftools view -H {dest} | grep -c CLNSIG || true"]).strip()
    console.print(f"  ClinVar-annotated records: {n}")
    return dest


# ClinVar's benign calls, excluded from the AF scope. Matched exactly, NOT by regex: bcftools'
# `!~` silently fails to filter this field (verified: 3,792 rows in, 3,792 out), while `!=` works.
_CLNSIG_BENIGN = ("Benign", "Benign/Likely_benign", "Likely_benign")


def _gnomad_scope(src: Path) -> str | None:
    """A bcftools ``-i`` expression selecting the few positions where AF changes an answer:
    AlphaMissense-pathogenic, or ClinVar-classified as anything other than benign.

    Deliberately tight (~1k variants here, not 52k and certainly not 5.1M). AF is meaningless for
    the ~3M intergenic / ~1.7M intronic variants a WGS carries, and it only *changes* a conclusion
    for variants already suspicious — chiefly `predicted_damaging`'s rarity filter. Every extra
    variant is another rsID batched to a rate-limited public API, so breadth here is paid for in
    setup time and in getting throttled. Ad-hoc "how common is X" is better answered live.

    Returns None when nothing upstream ran — the caller then skips rather than fetching everything.
    """
    have = _present_info(src)
    terms: list[str] = []
    if "am_class" in have:
        terms.append('INFO/am_class~"pathogenic"')
    if "CLNSIG" in have:
        not_benign = " && ".join(f'INFO/CLNSIG!="{b}"' for b in _CLNSIG_BENIGN)
        terms.append(f'(INFO/CLNSIG!="." && {not_benign})')
    return " || ".join(terms) if terms else None


def _af_cache_path() -> Path:
    return settings.work_dir / "gnomad.af.cache.json"


def _load_af_cache() -> dict[str, dict]:
    """rsID -> {allele: [af, grpmax, pop]} persisted across runs.

    gnomAD frequencies are a fixed release — refetching thousands of rsIDs on every weekly
    refresh would be pure waste, so results are cached and only unseen rsIDs are looked up.
    """
    import json

    p = _af_cache_path()
    if not p.exists():
        return {}
    try:
        data = json.loads(p.read_text())
        return data if isinstance(data, dict) else {}
    except (json.JSONDecodeError, OSError):
        return {}


def _save_af_cache(cache: dict[str, dict]) -> None:
    import json

    p = _af_cache_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(cache))
    tmp.replace(p)  # atomic: a crash can't truncate the cache


def _ensembl_gnomad_af(rsids: list[str], batch: int = 25) -> dict[str, dict]:
    """rsID -> {alt_allele: (af_all, grpmax_af, grpmax_pop)} from Ensembl's gnomAD frequencies.

    Batched POSTs; best-effort with retries — a failed batch costs those variants their AF, not
    the whole step. Only rsIDs leave the machine (same posture as gwas._resolve_rsids). Results
    are cached on disk so re-runs (every weekly refresh) fetch only genuinely new rsIDs.

    Batch size and timeout are empirically pinned, not arbitrary. `pops=1` returns ~120 population
    records per variant, so the response — not the request count — is the constraint: 25 ids is
    ~220 KB and answers in 20-36s, while 50 ids reliably times out. The timeout must stay well
    above that 36s or it rejects responses that were about to succeed. Don't raise `batch`.
    """
    import time

    import httpx

    cache = _load_af_cache()
    out: dict[str, dict] = {r: {a: tuple(v) for a, v in cache[r].items()} for r in rsids if r in cache}
    todo = [r for r in rsids if r not in cache]
    if out:
        console.print(f"  gnomAD AF: {len(out):,} rsIDs from cache; {len(todo):,} to fetch")
    if not todo:
        return out

    fetched, misses = 0, 0
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        data = None
        for attempt in range(2):
            try:
                r = httpx.post(ENSEMBL_VARIATION,
                               headers={"Content-Type": "application/json", "Accept": "application/json"},
                               json={"ids": chunk}, params={"pops": "1"}, timeout=90)
                if r.status_code in (429, 500, 502, 503, 504):  # throttled or unwell — do not hammer
                    raise RuntimeError(f"HTTP {r.status_code}")
                r.raise_for_status()
                data = r.json()
                break
            except Exception as e:  # noqa: BLE001 - network is best-effort
                console.print(f"[yellow]  Ensembl AF batch failed (try {attempt + 1}/2): {e}[/]")
                time.sleep(3 * (attempt + 1))
        if not data:
            # Bail out rather than grind: a rate-limited/unwell API won't recover inside this run,
            # and hammering it is how you get blocked. AF is optional — skip it and move on.
            misses += 1
            if misses >= 2:
                console.print("[yellow]  Ensembl unavailable — skipping gnomAD AF for this run "
                              "(cached results kept; re-run later).[/]")
                break
            continue
        misses = 0
        for rsid, info in data.items():
            per: dict[str, tuple] = {}
            for p in info.get("populations") or []:
                name = str(p.get("population", ""))
                if not name.startswith("gnomADg:"):
                    continue
                grp = name.split(":", 1)[1]
                allele, freq = p.get("allele"), p.get("frequency")
                if allele is None or freq is None:
                    continue
                af_all, gmax, gpop = per.get(allele, (None, None, None))
                if grp == "ALL":
                    af_all = freq
                elif grp not in _GNOMAD_GRPMAX_SKIP and (gmax is None or freq > gmax):
                    gmax, gpop = freq, grp
                per[allele] = (af_all, gmax, gpop)
            out[rsid] = per
            cache[rsid] = {a: list(v) for a, v in per.items()}  # cache misses too (empty = no AF)
        fetched += len(chunk)
        console.print(f"  Ensembl AF: {fetched:,}/{len(todo):,} rsIDs fetched")
        _save_af_cache(cache)  # checkpoint: an interrupted run keeps its progress
    return out


def annotate_gnomad(src: Path, dest: Path) -> Path:
    """Transfer gnomAD genome allele frequencies, sourced from Ensembl REST (keyed by rsID).

    Runs *last* in the chain so it can scope lookups to positions where AF actually means
    something (see ``_gnomad_scope``) — ~50k rather than ~5.1M. Variants without an rsID, and
    indels whose Ensembl allele representation doesn't match the VCF ALT, simply get no AF.
    """
    expr = _gnomad_scope(src)
    filt = f"-i '{expr}' " if expr else ""
    if not expr:
        console.print("[yellow]gnomAD: nothing to scope by — annotating every position's rsID.[/]")
    rows = shell.capture(["bash", "-o", "pipefail", "-c",
                          f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\n' {filt}{src}"]).splitlines()
    variants = []
    for ln in rows:
        f = ln.split("\t")
        if len(f) != 5:
            continue
        rsid = f[4].split(";")[0]
        if rsid.startswith("rs"):
            variants.append((f[0], f[1], f[2], f[3], rsid))
    if not variants:
        console.print("[yellow]gnomAD: no scoped variants carry an rsID — skipping.[/]")
        return src

    rsids = sorted({v[4] for v in variants})
    console.print(f"Fetching gnomAD allele frequencies from Ensembl ({len(rsids):,} variants)…")
    afs = _ensembl_gnomad_af(rsids)
    if not afs:
        console.print("[yellow]gnomAD: Ensembl returned no frequencies; leaving variants without AF.[/]")
        return src

    work = settings.work_dir
    tsv = work / "gnomad.af.tsv"
    n = 0
    with tsv.open("w") as fh:
        for chrom, pos, ref, alt, rsid in variants:
            hit = (afs.get(rsid) or {}).get(alt)
            if not hit:
                continue
            af_all, gmax, gpop = hit
            if af_all is None and gmax is None:
                continue
            fh.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{'.' if af_all is None else af_all}\t"
                     f"{'.' if gmax is None else gmax}\t{gpop or '.'}\n")
            n += 1
    if n == 0:
        console.print("[yellow]gnomAD: no frequencies matched the scoped ALT alleles — skipping.[/]")
        return src

    shell.sh(f"sort -k1,1 -k2,2n {tsv} | bgzip -f > {tsv}.gz")
    shell.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", f"{tsv}.gz"])
    header = work / "gnomad.header.txt"
    _d = "via Ensembl"
    header.write_text(
        f'##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD genomes AF ({_d})">\n'
        f'##INFO=<ID=gnomAD_AF_grpmax,Number=1,Type=Float,Description="Max gnomAD AF across populations ({_d})">\n'
        f'##INFO=<ID=gnomAD_grpmax,Number=1,Type=String,Description="gnomAD population with the max AF ({_d})">\n'
    )
    shell.run(["bcftools", "annotate", "-a", f"{tsv}.gz", "-h", str(header),
               "-c", "CHROM,POS,REF,ALT,gnomAD_AF,gnomAD_AF_grpmax,gnomAD_grpmax",
               str(src), "-Oz", "-o", str(dest)])
    _index(dest)
    console.print(f"  gnomAD-annotated records: {n:,}")
    return dest


def annotate_snpeff(src: Path, dest: Path) -> Path:
    jar = settings.annotations_dir / "snpEff" / "snpEff.jar"
    if not jar.exists():
        console.print("[yellow]SnpEff missing — skipping consequences. Run `locus download snpeff`.[/]")
        return src
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — skipping SnpEff. `brew install openjdk`.[/]")
        return src
    console.print("Annotating functional consequences (SnpEff)…")
    # snpEff writes uncompressed VCF to stdout; bgzip it.
    cmd = " ".join(shell.java_cmd(["-Xmx6g", "-jar", str(jar), "-noStats", download.SNPEFF_DB, str(src)]))
    shell.sh(f"{cmd} | bgzip -c > {dest}")
    _index(dest)
    return dest


def annotate_alphamissense(src: Path, dest: Path) -> Path:
    """Annotate missense pathogenicity from AlphaMissense (calibrated score for ~every missense).

    Fills the gap where ClinVar is silent: a high am_pathogenicity on a variant ClinVar has never
    seen is real signal, not 'nothing'.
    """
    am = settings.annotations_dir / "alphamissense" / "AlphaMissense_hg38.slim.tsv.bgz"
    if not am.exists():
        console.print("[yellow]AlphaMissense missing — skipping. Run `locus download alphamissense`.[/]")
        return src
    header = settings.work_dir / "am.header.txt"
    header.write_text(
        '##INFO=<ID=am_pathogenicity,Number=1,Type=Float,Description="AlphaMissense pathogenicity (0-1)">\n'
        '##INFO=<ID=am_class,Number=1,Type=String,Description="AlphaMissense class (benign/ambiguous/pathogenic)">\n'
    )
    console.print("Annotating AlphaMissense missense pathogenicity…")
    shell.run([
        "bcftools", "annotate", "-a", str(am), "-h", str(header),
        "-c", "CHROM,POS,REF,ALT,am_pathogenicity,am_class",
        str(src), "-Oz", "-o", str(dest),
    ])
    _index(dest)
    n = shell.capture(["bash", "-c", f"bcftools view -H {dest} 2>/dev/null | grep -c am_pathogenicity || true"]).strip()
    console.print(f"  AlphaMissense-annotated records: {n}")
    return dest


def _pharmcat_input(inputs, reference: Path) -> Path:
    """Build PharmCAT's input VCF: PGx regions only, chr-prefixed, hom-ref blocks expanded.

    Restrict the gVCF to PharmCAT's regions first (small subset), rename contigs to
    chr-prefixed, then ``bcftools convert --gvcf2vcf`` to expand ``ALT="."`` hom-ref
    blocks into explicit per-position ``0/0`` calls — so PGx sites that are reference
    get genotyped instead of becoming no-calls.
    """
    install = artifacts.pharmcat_install_dir()
    regions = install / "pharmcat_regions.bed"
    pgx_input = artifacts.pharmcat_input_vcf()
    info = read_info(inputs.small_variants)

    if info.chr_prefixed:
        region_src = f"bcftools view -R {regions} {inputs.small_variants}"
    else:
        # Match the gVCF's non-chr contigs (chrM→MT), restrict, then rename to chr.
        nochr = settings.work_dir / "pharmcat_regions.nochr.bed"
        shell.sh(f"sed -e 's/^chrM\\t/MT\\t/' -e 's/^chr//' {regions} > {nochr}")
        rename = settings.work_dir / "contigs2chr.txt"
        if not rename.exists():
            from .vcfutils import chr_rename_map, write_rename_file
            write_rename_file(chr_rename_map(info.contigs), rename)
        region_src = f"bcftools view -R {nochr} {inputs.small_variants} | bcftools annotate --rename-chrs {rename}"

    # GT-only: PharmCAT only needs genotypes, and keeping AD/DP makes its internal
    # `bcftools norm` fail to merge per-allele FORMAT tags ("could not merge AD").
    shell.sh(
        f"{region_src} | bcftools convert --gvcf2vcf -f {reference} -Ou "
        f"| bcftools annotate -x '^FORMAT/GT' -Oz -o {pgx_input}"
    )
    _index(pgx_input)
    return pgx_input


def annotate_pharmcat() -> Path | None:
    """Run PharmCAT natively (jar + Python preprocessor) on the PGx-restricted VCF."""
    pipeline = artifacts.pharmcat_install_dir() / "pharmcat_pipeline"
    if not pipeline.exists():
        console.print("[yellow]PharmCAT not installed — run `locus download pharmcat`.[/]")
        return None
    if shell.resolve_java() is None:
        console.print("[yellow]No working Java found — PharmCAT needs Java 17+. `brew install openjdk`.[/]")
        return None
    inputs = artifacts.classify_inputs(settings.genome_dir)
    if not inputs.small_variants:
        console.print("[yellow]No small-variant gVCF for PharmCAT.[/]")
        return None
    reference = artifacts.find_reference()
    if reference is None:
        console.print("[yellow]Reference FASTA needed for PharmCAT. Run `locus download reference`.[/]")
        return None

    out_dir = artifacts.pharmcat_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    console.print("Preparing PharmCAT input (PGx regions, chr-prefixed, hom-ref expanded)…")
    pgx_input = _pharmcat_input(inputs, reference)

    console.print("Running PharmCAT (native)…")
    # PharmCAT calls `java`; put the working JDK first on PATH for the subprocess.
    import os
    import sys

    env = dict(os.environ)
    java = shell.resolve_java()
    if java:
        env["PATH"] = f"{Path(java).parent}:{env.get('PATH', '')}"
        env["JAVA_HOME"] = str(Path(java).parent.parent)
    shell.run_env(
        [sys.executable, str(pipeline), str(pgx_input), "-o", str(out_dir),
         "-reporterJson", "-reporterCallsOnlyTsv"],
        env=env,
    )
    report = next(out_dir.glob("*.report.json"), None)
    console.print(f"[green]PharmCAT report:[/] {report}" if report else "[yellow]No PharmCAT report produced.[/]")
    return report


def run(steps: str = "all") -> Path:
    """Run the requested annotation steps, producing the annotated sites VCF."""
    src = artifacts.sites_vcf()
    if not src.exists():
        raise FileNotFoundError(f"No sites VCF ({src}). Run `locus ingest` first.")

    requested = ALL_STEPS if steps in ("all", "") else [s.strip() for s in steps.split(",")]
    console.rule(f"[bold]Annotate ({', '.join(requested)})")

    # bcftools annotate needs the input indexed — self-heal in case ingest's index is missing.
    from .ingest import ensure_index as _ensure_index

    _ensure_index(src)

    work = settings.work_dir
    cur = src
    if "clinvar" in requested:
        cur = annotate_clinvar(cur, work / f"{settings.sample_id}.clinvar.vcf.gz")
    if "snpeff" in requested or "vep" in requested:
        cur = annotate_snpeff(cur, work / f"{settings.sample_id}.snpeff.vcf.gz")
    if "alphamissense" in requested:
        cur = annotate_alphamissense(cur, work / f"{settings.sample_id}.am.vcf.gz")
    # gnomAD LAST: it scopes its (expensive, streamed) AF lookups to the positions the steps above
    # marked interesting — ClinVar-annotated, AlphaMissense-scored, or coding. Running it earlier
    # would leave it nothing to scope by and force a ~5.1M-position fetch.
    if "gnomad" in requested:
        cur = annotate_gnomad(cur, work / f"{settings.sample_id}.gnomad.vcf.gz")

    # Finalize the small-variant annotated VCF.
    dest = artifacts.annotated_vcf()
    if cur != src:
        shell.run(["bcftools", "view", str(cur), "-Oz", "-o", str(dest)])
        _index(dest)
        console.print(f"[green]Annotated VCF ready[/] → {dest}")
    else:
        console.print("[yellow]No small-variant annotations applied (no DBs present).[/]")

    if "pharmcat" in requested:
        annotate_pharmcat()

    return dest if dest.exists() else src

# Setting up Locus (for non-developers)

This guide gets Locus running on your own Mac so you can ask **Claude** about your
[sequencing.com](https://sequencing.com) whole-genome data. It's a one-time, ~20–40 minute
setup, and **everything stays on your computer** — nothing is uploaded anywhere.

**You need:** an **Apple Silicon Mac** (M1/M2/M3/M4), **Claude Desktop** (or Claude Code)
installed, ~25 GB of free disk, and your sequencing.com 30× WGS files.

---

## 1. Get the Locus folder

Open **Terminal** and clone the repo:

```bash
git clone https://github.com/mattsmcnulty/locus.git
cd locus
```

*(First time using `git`? macOS may pop up "Install Command Line Developer Tools" — click
**Install**, let it finish, then run those two commands again.)*

## 2. Download your genome from sequencing.com

Log in → open your **30× Whole Genome Sequencing** order → **Download**, and save these files
into the `data/genome/` folder inside the Locus folder:

- **`…snp-indel.genome.vcf.gz`** — required
- `…cnv.vcf.gz` and `…sv.vcf.gz` — recommended (copy-number & structural variants)

Keep the filenames exactly as downloaded. (The big `.fastq`/`.fq.gz` read files are **not**
needed — skip them.)

## 3. Run the setup

In Finder, **double-click `setup.command`** in the Locus folder. A Terminal window opens and
it does the rest:

- installs the tools it needs (Homebrew, etc. — it may ask for your Mac password once, and may
  pop a "Install Command Line Tools" dialog — click Install, then double-click setup again);
- downloads ~9 GB of reference databases (one time);
- builds your private genome database and runs all the interpretation;
- registers Locus with Claude.

It's **safe to close and re-run** — it picks up where it left off. Total time is mostly the
download: ~20–40 minutes depending on your connection.

## 4. Restart Claude and ask away

**Fully quit Claude Desktop (Cmd-Q, not just close the window) and reopen it.** Then ask things
like:

- "Using my locus genome, give me an overview of what's loaded."
- "What's my biogeographic ancestry, and where do my ancestors come from?"
- "Do I have any ACMG secondary findings I should know about?"
- "What does my genome say about caffeine, lactose, and alcohol?"
- "What are my pharmacogenomic results for clopidogrel and warfarin?"
- "What's my polygenic risk for coronary artery disease, as a percentile?"

Prefer clicking around? Open the **Locus** app (it was added to your Applications), or run
`locus serve api` and visit http://127.0.0.1:8787.

---

## If something looks off

Run `uv run locus doctor` in the Locus folder — it shows what's installed, whether your genome
is detected (and is GRCh38), and whether Claude is registered.

## Important — please read

Locus is for **personal exploration and education. It is not medical advice and not a medical
device.** Polygenic risk percentiles are research-grade estimates and only meaningful within a
matched ancestry; single-variant ("GWAS") associations are weak; the HLA-B\*57:01 result is a
screening proxy. **Discuss anything health-relevant with a qualified clinician or genetic
counselor.** Your genetic data never leaves your Mac.

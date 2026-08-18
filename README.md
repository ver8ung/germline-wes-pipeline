# NGS germline short-variant pipeline (WES, GRCh38)

End-to-end, reproducible variant calling from **raw FASTQ reads → annotated VCF**,
following [GATK germline short-variant best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932).
Orchestrated with **Snakemake**, every tool pinned in its own **Conda** env, and
the whole thing wrapped in a single **Docker** image so it runs the same on
Windows, macOS, Linux, or an HPC node.

- **Analysis:** germline SNPs + small indels
- **Assay:** whole-exome (WES), restricted to your capture kit's targets
- **Reference:** GRCh38 (Broad/GATK bundle, `chr`-named contigs)
- **Cohort size:** any `n ≥ 1` (joint genotyping + hard-filtering)

---

## Pipeline stages

```
                 ┌─ FastQC (raw QC) ─────────────────────────────┐
 FASTQ (R1/R2) ──┤                                                │
   per unit      └─ fastp (trim/adapter) ─ bwa-mem ─ sort ──┐     │
                                                            │     │
              merge units → MarkDuplicates → BQSR  ◄────────┘     │
                                  │                               │
              HaplotypeCaller (GVCF, per sample, on targets)      │
                                  │                               │
              CombineGVCFs → GenotypeGVCFs   (cohort joint call)  │
                                  │                               │
              SelectVariants → VariantFiltration → MergeVcfs      │  (QC metrics
                                  │                               │   from every
              SnpEff  (+ optional SnpSift / ClinVar)              │   stage)
                                  │                               │
        results/annotated/all.annotated.vcf.gz                    │
                                                                  ▼
                                              MultiQC ─ results/qc/multiqc_report.html
```

| Stage | Tool | Key output |
|-------|------|-----------|
| Raw QC | FastQC | `results/qc/fastqc/` |
| Trim | fastp | `results/trimmed/`, `results/qc/fastp/` |
| Align | bwa-mem + samtools | `results/mapped/` (temp) |
| Dedup | GATK MarkDuplicates | `results/dedup/` (temp) + metrics |
| Recalibrate | GATK BQSR | `results/bqsr/*.recal.bam` |
| Call | GATK HaplotypeCaller (GVCF) | `results/called/*.g.vcf.gz` |
| Joint genotype | CombineGVCFs + GenotypeGVCFs | `results/genotyped/cohort.vcf.gz` |
| Filter | GATK hard-filter | `results/filtered/all.filtered.vcf.gz` (flagged) |
| Select PASS | bcftools | `results/filtered/all.pass.vcf.gz` |
| Annotate | SnpEff (+ SnpSift) | **`results/annotated/all.annotated.vcf.gz`** (PASS-only) |
| Report | MultiQC | `results/qc/multiqc_report.html` |

---

## Layout

```
ngs-germline-wes/
├── config/
│   ├── config.yaml             # all knobs: reference URLs, intervals, filters, annotation
│   ├── samples.tsv             # one row per biological sample
│   ├── units.tsv               # one row per lane/library (FASTQ pair) of a sample
│   ├── twist_onso.yaml         # overlay: the higher-depth Twist-Onso run
│   ├── units.twist_onso.tsv    #   its unit sheet (data not publicly available)
│   ├── giab_trio.yaml          # overlay: the GIAB Ashkenazim trio (n=3, joint calling)
│   ├── samples.giab_trio.tsv   #   its sample sheet, with pedigree columns
│   └── units.giab_trio.tsv     #   its unit sheet (reads reverted from BAM)
├── benchmarks/                 # committed accuracy evidence, one dir per run
├── workflow/
│   ├── Snakefile               # includes the rule modules, defines `all` + `setup_reference`
│   ├── rules/*.smk             # common, resources, capture_bed, qc, trim, map,
│   │                           #   dedup_bqsr, calling, filtering, annotation,
│   │                           #   stats, benchmark, smoke, bam_revert
│   └── envs/*.yaml             # pinned per-tool Conda environments
├── docker/
│   ├── Dockerfile              # Snakemake + conda/mamba orchestrator image
│   └── base-image.txt          # base image tag, passed in as a build-arg
├── .tests/smoke/               # tiny simulated-read fixture for `smoke`
├── dev/                        # author analysis tooling; NOT part of the pipeline
├── data/                       # ← put your FASTQ here (gitignored)
├── resources/                  # ← downloaded reference/known-sites/db (gitignored)
├── results/                    # ← pipeline outputs (gitignored)
├── run.ps1                     # Windows launcher
├── run.cmd                     # cmd.exe wrapper around run.ps1
└── run.sh                      # Linux/macOS/WSL launcher
```

---

## Prerequisites

- **Docker Desktop** (Windows/macOS) or Docker Engine (Linux). Nothing else is
  installed on the host — Snakemake and every bioinformatics tool live in the image
  / per-rule Conda envs.
- **Memory:** give Docker **≥ 16 GB** RAM (Docker Desktop → Settings → Resources).
  `bwa index` of GRCh38 and GATK each need several GB. The launchers request 16 GB
  by default, so provisioning less than that means asking Docker for more than the
  VM has. Override with `WES_DOCKER_MEMORY` (see *Environment variables* below).
- **Disk:** budget **~70 GB** for the documented GIAB walkthrough, not 30:
  - ~11 GB reference bundle + known sites + SnpEff DB + bwa index (`resources/`)
  - ~16 GB input FASTQ for the default Garvan run (`data/`)
  - ~35 GB working space and outputs (`results/`) — the recalibrated BAM alone is
    the largest single file.

  A run on your own smaller exome needs correspondingly less, but the reference
  bundle is a fixed ~11 GB floor.
- **Network:** the first run downloads the reference, known-sites VCFs, and the
  SnpEff database.

---

## Quick start

**1. Point the sample sheets at your data.** Edit `config/units.tsv`
(tab-separated — one row per FASTQ pair):

```tsv
sample	unit	platform	fq1	fq2
NA12878	L001	ILLUMINA	data/NA12878_L001_R1.fastq.gz	data/NA12878_L001_R2.fastq.gz
```

and `config/samples.tsv` (one row per sample). A sample may span several `unit`
rows (multiple lanes/libraries) — they're aligned separately with distinct read
groups and merged before duplicate marking.

**2. Supply your exome capture targets.** The default `config/config.yaml` →
`intervals.bed` is `resources/intervals/nextera_expandedexome.GRCh38.bed`, which
`.\run get_capture_bed` builds for you (it is the kit that produced the default
GIAB dataset, so kit and BED match out of the box). To use a different kit, point
`intervals.bed` at your own GRCh38, `chr`-named target BED from the vendor
(Agilent SureSelect, IDT xGen, Twist, …); three other providers are already wired
up — see the `intervals` block in `config/config.yaml`.

**3. Run it.**

```powershell
# Windows, from the project root. Note the leading .\ — PowerShell does not
# resolve commands from the current directory, so a bare `run` fails with
# CommandNotFoundException. (In cmd.exe, bare `run` does work.)
.\run smoke -n                # dry-run the smoke test (no downloads)
.\run setup_reference         # download + index reference (once, ~1.5–2 h)
.\run --cores 8               # full pipeline → results/annotated/all.annotated.vcf.gz
```

> `setup_reference` is dominated by `bwa index` on GRCh38, which is single-threaded
> and takes roughly 45–60 minutes by itself on top of the ~4.9 GB of downloads
> (`logs/resources/bwa_index.log` recorded 2,762 s = 46 min here). It is a
> one-time cost — the index is reused by every later run.

> **Why not just type `run.ps1`?** Double-clicking a `.ps1` or typing its bare
> name opens it in **Notepad** — Windows associates `.ps1` with *edit*, not *run*.
> Use `run.cmd` (above), or invoke the script explicitly in PowerShell with the
> `.\` prefix: `.\run.ps1 --cores 8`. If PowerShell blocks it with an
> execution-policy error, either run via `run.cmd`, or allow local scripts once:
> `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

```bash
# Linux / macOS / WSL:
./run.sh -n
./run.sh setup_reference
./run.sh --cores 8
```

The launchers build the Docker image on first use, bind-mount the project at
`/workflow`, and run `snakemake --use-conda` (per-rule envs are created once and
cached under `.snakemake/conda`). Any extra arguments pass straight through to
Snakemake, e.g. `./run.sh --cores 8 -p` or a specific target file.

> The image **build** is the slow part — it bakes all tool envs in (budget
> ~15–30 min, once). After that the first **run** only needs the one-time
> reference download; the per-rule envs are already inside the image.

---

## Image versioning & rebuilds

The launchers keep each image in lockstep with the code that defines it:

- **Content-tagged images.** The tag is `src-<hash>`, a SHA256 over exactly the
  files that determine what the image contains: `workflow/envs/*.yaml`, the
  Dockerfile, and `docker/base-image.txt`. Change any of them and you get a new
  tag and a rebuild; change anything else and the existing image is reused.
  `:latest` also points at the most recent build.

  It is deliberately **not** the git commit SHA. Tagging by commit meant every
  commit invalidated the tag and forced a needless 15–30 min rebuild, while two
  commits with identical image inputs pointlessly built the same image twice.

  Both launchers compute this hash identically (same file order, one hash over the
  concatenated bytes, lowercase, first 12 chars), so a tree checked out on Windows
  and on Linux resolves to the same image rather than building two.
- **Pinnable base.** The base image is read from
  [docker/base-image.txt](docker/base-image.txt) and passed as a build-arg. It
  currently holds a floating tag, so two people building from the same commit can
  get different base layers — pin it to a digest there if you need byte-level
  reproducible rebuilds (instructions in the file).
- **Staleness guard.** Before every run the launcher runs the image's baked
  `check-bake` script, which re-hashes the bind-mounted `workflow/envs/*.yaml` and
  compares to the digest baked at build time. If they differ you get a loud
  WARNING (otherwise the per-rule envs would be silently re-created every run).
  Re-bake with `docker rmi <image>` then re-run the launcher.

**Needs a rebuild:** changes to `workflow/envs/*.yaml` (tool versions) or the
Dockerfile. **Does *not*:** rule logic, config, sample sheets, params, intervals —
those are bind-mounted and take effect on the next run.

## Configuration highlights (`config/config.yaml`)

- **`ref`** — GRCh38 FASTA URL (Broad bundle). Downloaded + indexed
  (`.fai`, `.dict`, bwa index) automatically.
- **`known_sites`** — dbSNP, Mills gold-standard indels, 1000G known indels, for BQSR.
- **`intervals`** — capture BED + `padding` (default 100 bp each side) applied at
  HaplotypeCaller/BQSR time.
- **`filtering`** — GATK hard-filter thresholds, split for SNPs vs indels. Edit to taste.
- **`annotation.snpeff.db`** — defaults to **`hg38`** (UCSC/RefSeq, `chr`-named) so it
  matches the Broad reference contigs with no renaming. See *Design notes*.
- **`annotation.clinvar.activate`** — set `true` to overlay ClinVar IDs (downloaded
  and contig-renamed to match). Off by default.
- **`resources`** — per-step threads / memory; tune for laptop vs server.

---

## Outputs

- **`results/annotated/all.annotated.vcf.gz`** — the headline: joint-called,
  hard-filtered, functionally annotated multi-sample VCF (+ `.tbi`). This is
  **PASS-only**: records that failed a hard filter are not in it. Use it directly;
  no further filtering is needed.
- **`results/filtered/all.filtered.vcf.gz`** — the same callset *before* PASS
  selection, with every record present and its `FILTER` tag intact. Reach for this
  when you want to see what was rejected and why — computing a filter-failure rate,
  or hunting rescue candidates near a threshold. `results/filtered/all.pass.vcf.gz`
  is the PASS subset that annotation consumes.
- **`results/qc/multiqc_report.html`** — one report aggregating FastQC, fastp,
  MarkDuplicates, samtools, bcftools, and SnpEff metrics across all samples.

  > Two scopes appear in the outputs, deliberately. `bcftools stats` (and therefore
  > the MultiQC variant counts) describes the **PASS** callset, matching the
  > deliverable. The hap.py benchmark scores the **flagged** VCF, because hap.py
  > reads the `FILTER` column itself and reports ALL and PASS rows separately — it
  > needs the rejected records to do that. So the two will not agree, and should
  > not.
- Intermediates under `results/` (recalibrated BAMs, per-sample GVCFs, the raw
  genotyped VCF). Per-rule logs under `logs/`.

---

## Design notes & gotchas

- **Contig naming.** The Broad GRCh38 bundle, its known-sites VCFs, and SnpEff
  `hg38` all use `chr1…chr22,chrX,chrY,chrM` — kept consistent end-to-end. If you
  switch to the Ensembl SnpEff DB (`GRCh38.105`, no `chr`), you must rename contigs
  around annotation (`bcftools annotate --rename-chrs`). ClinVar (Ensembl-named) is
  already handled: when enabled it's downloaded and renamed to `chr` once.
- **Hard-filtering, not VQSR.** VQSR needs a large cohort to train; hard-filtering
  works for any size including a single exome. For large cohorts (≳ 30 WES) switch
  to VQSR (`VariantRecalibrator` / `ApplyVQSR`) for better sensitivity/specificity.
- **CombineGVCFs vs GenomicsDBImport.** `CombineGVCFs` is used for portability and
  robustness inside containers. For large cohorts, `GenomicsDBImport` scales better
  — swap it into `calling.smk` and feed its workspace to `GenotypeGVCFs`.
- **bwa (classic), not bwa-mem2.** Classic `bwa` index is ~6 GB resident; bwa-mem2 is
  faster but needs ~30 GB RAM for the GRCh38 index — impractical on a laptop. Swap in
  `map.smk`/`envs/bwa.yaml` if you have the memory.
- **No interval scatter.** HaplotypeCaller runs over the whole target list per
  sample. For WES that's fine; for WGS or speed, scatter by interval and gather.
- **Pre-baked Conda envs.** Each rule's tools come from a pinned Conda env
  ([workflow/envs/](workflow/envs/)), created at **image build time** into
  `/opt/snakemake-envs` and reused at run time via `--conda-prefix` (the launchers
  pass it), so runs are offline with no first-run env creation. Snakemake matches a
  baked env by its YAML *content* hash — see **Image versioning & rebuilds** above
  for how the staleness guard catches drift. On Linux/HPC you can instead attach
  per-rule `container:` directives and run `snakemake --sdm apptainer` for true
  per-tool images.

---

## Handy commands

```bash
./run.sh -n                                  # dry-run (plan only)
./run.sh --cores 8 -p                        # full run, print shell commands
./run.sh setup_reference                     # only fetch + index reference data
./run.sh --cores 4 results/called/NA12878.g.vcf.gz   # build one target
./run.sh --dag | dot -Tsvg > dag.svg         # render the DAG (needs graphviz)
./run.sh --report report.html                # Snakemake provenance report (after a run)
docker run --rm -it -v "${PWD}:/workflow" -w /workflow ngs-germline-wes bash  # poke around
```

Inside that last shell, `smoke` is **not** a Snakemake target — it is a shortcut the
launchers expand to `--configfile .tests/smoke/config.yaml --cores 4`. Use the full
form when you are driving `snakemake` directly.

## Environment variables

Both launchers read these; they are optional.

| Variable | Default | Effect |
|---|---|---|
| `WES_DOCKER_MEMORY` | `16g` | Memory limit passed to `docker run`, and the basis for Snakemake's `--resources mem_mb` budget (90% of it). Must fit inside the Docker Desktop / WSL2 VM limit. Raise it for large multi-unit runs, e.g. `32g`. |
| `WES_SCRATCH_BAMS` | unset | When set to `1`, mounts Docker **named volumes** over `results/{mapped,dedup,bqsr}` so heavy BAM I/O hits VM-native storage instead of the (slow) Windows bind mount. |

`WES_SCRATCH_BAMS` has consequences worth knowing before you enable it:

- The recalibrated BAM then lives **on a Docker volume, not on your disk**, so host
  tools cannot see `results/bqsr/`. To read it (e.g. for `dev/diag_capture.sh`), add
  `--mount type=volume,source=wes-bqsr,target=/workflow/results/bqsr` to your
  `docker run`.
- The volumes persist between runs. Reset them with
  `docker volume rm wes-mapped wes-dedup wes-bqsr`.
- Toggling it mid-project remaps those paths, so use it for fresh runs rather than
  switching partway through.

```powershell
$env:WES_DOCKER_MEMORY = '32g'      # PowerShell
```
```bash
WES_DOCKER_MEMORY=32g ./run.sh --cores 8    # bash
```

> **There is a floor.** The launchers derive Snakemake's `--resources mem_mb`
> budget from this value (90% of it), and refuse to start if the budget is smaller
> than the largest single rule's reservation — `bwa_mem` and `benchmark` each
> reserve 12000 MiB, so anything below 16g is rejected up front with a message
> naming the offending rule. That check exists because the failure it replaces was
> miserable: Snakemake only raises the error once the oversized job becomes
> schedulable, i.e. after FastQC and trimming have already run, and the message
> blames pipes.
>
> This is why `resources.benchmark.mem_mb` is set *equal* to `bwa_mem`'s rather
> than higher: the floor is the maximum over all rules, so raising it above 14745
> would break the default 16g run for everyone, smoke test included.

> **RTG sizes its JVM heap from the host, not the container.** `rtg format` and
> `vcfeval` default to 90% of *detected physical* RAM, and inside Docker
> `/proc/meminfo` still reports the **host's** memory — so the JVM happily
> overshoots `--memory` and the container is OOM-killed with no Java stack trace at
> all. Both rules that invoke rtg export `RTG_MEM` derived from their `mem_mb`
> reservation. The documented `rtg RTG_MEM=8g <cmd>` first-argument form is not
> reachable here because hap.py invokes rtg itself; the environment variable is the
> only lever. If you add another rtg rule, export it there too.

## Smoke test (built in)

A self-contained smoke test exercises the **entire DAG** without your own data or
a full exome. It simulates ~40x paired reads (wgsim) over a 100 kb chr20 window
and runs them all the way to an annotated VCF. It uses a config overlay
([.tests/smoke/](.tests/smoke/)) that repoints only the inputs — reference, known
sites, filters, and SnpEff DB are inherited unchanged.

```powershell
.\run.ps1 smoke -n     # dry-run: resolves the whole plan, downloads nothing
.\run.ps1 smoke        # full run: simulate -> align -> call -> filter -> annotate
```
```bash
./run.sh smoke -n
./run.sh smoke
```

Equivalent to `snakemake --configfile .tests/smoke/config.yaml --cores 4 [...]`.

- **`smoke -n`** is the fast wiring check: it confirms every rule connects with
  **zero downloads** and zero compute — do this first.
- **`smoke`** (no `-n`) still needs the GRCh38 reference (download + BWA index is
  unavoidable for real alignment/BQSR/calling), but after that the run finishes in
  minutes instead of hours — a true end-to-end validation before you commit to a
  full exome.

To smoke-test with real data instead, grab a small public exome (e.g. a Genome in
a Bottle NA12878 subset), point `units.tsv` at the FASTQs, and use a small BED.

## Evaluating against GIAB

The repo's **default configuration is this benchmark run**, so no config edits are
needed: `config/samples.tsv` + `config/units.tsv` describe one `NA12878` sample
with 4 units (two libraries × two lanes), and `intervals.bed` already points at
the Nextera BED — the kit that produced these reads. Kit and evaluation region
therefore match, which is what makes the resulting numbers meaningful.

1. **Download the WES reads** (≈16 GB, 8 files) into `data\garvan\` (where
   `config/units.tsv` expects them):
   ```powershell
   $u='https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/'
   $dst='data\garvan'; New-Item -ItemType Directory -Force $dst | Out-Null
   'NIST7035_TAAGGCGA_L001_R1_001','NIST7035_TAAGGCGA_L001_R2_001',
   'NIST7035_TAAGGCGA_L002_R1_001','NIST7035_TAAGGCGA_L002_R2_001',
   'NIST7086_CGTACTAG_L001_R1_001','NIST7086_CGTACTAG_L001_R2_001',
   'NIST7086_CGTACTAG_L002_R1_001','NIST7086_CGTACTAG_L002_R2_001' |
     ForEach-Object { curl.exe -L -o "$dst\$_.fastq.gz" "$u$_.fastq.gz" }
   ```

2. **Build the GRCh38 capture BED** (Illumina ships only hg19 → lifted over):
   ```powershell
   .\run get_capture_bed
   ```
   → `resources/intervals/nextera_expandedexome.GRCh38.bed` (what `intervals.bed` points at).

3. **Run the pipeline:** `.\run --cores 8` → `results/annotated/all.annotated.vcf.gz`

4. **Benchmark** vs. GIAB high-confidence calls (hap.py, scored on *confident
   regions ∩ capture targets*; truth files auto-download):
   ```powershell
   .\run benchmark
   ```
   → `results/benchmark/happy.NA12878.summary.csv` (precision / recall / F1 for
   SNPs & indels), plus `happy.NA12878.extended.csv` with the same metrics broken
   down by genomic context.

### How scoring works

**Benchmarking is per sample.** Truth sets are configured by sample name under
`benchmark.truth` in `config/config.yaml` (HG001/NA12878 and the GIAB Ashkenazim
trio ship configured), and a sample with no entry is simply not benchmarked — so
the `benchmark` target does the right thing for a cohort containing a mix.

Each sample's calls are extracted from the cohort callset before scoring, because
**hap.py cannot score a multi-sample VCF and has no flag to select a column.**
Under `vcfeval` it aborts; under the older `xcmp` engine it did not — it silently
scored the first sample column and reported the result under whichever truth set
you named. The extraction therefore runs at every cohort size, including n=1.

**Comparison engine** is `benchmark.engine`, default `vcfeval` (RTG). On the
`twist_onso` benchmark the two engines returned identical SNP recall, precision and F1
and identical TP/FN/FP/UNK counts, differing only in how they attribute false
positives between genotype and allele errors (`FP.gt`/`FP.al`); on indels exactly one
variant crossed from FP to TP, worth ~0.001 F1 — hap.py already left-shifts and
decomposes both callsets before comparing, so most variant-representation differences
are normalised before either engine sees them. Note what that control is: the two
engines scored independently generated callsets whose counts matched, not literally
the same VCF. vcfeval is the default because it is
the stricter and more widely cited engine, **not** because it improves the
numbers. One label = one engine: switching re-runs and replaces that label's
output rather than sitting alongside it.

vcfeval needs the reference in RTG's SDF format, built once by `rtg_format_sdf`
(~1 min, ~1.2 GB). It is deliberately **not** part of `setup_reference`, so users
who never benchmark do not pay for it.

**Stratification** (`benchmark.stratification`) scores 14 curated
[GIAB genome-stratification](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/)
subsets alongside the headline, turning one F1 into a map of *where* calling
fails — homopolymers, tandem repeats, segmental duplications, low mappability, GC
extremes, MHC, RefSeq CDS and their complements. This is the single most
informative thing in the benchmark output; see `happy.extended.csv`, where
`Subset=*` is the headline and the named subsets are diagnostic. The list is
deliberately curated rather than GIAB's full 298-entry manifest, which would take
`happy.extended.csv` from ~66 rows to ~6,600; widen it by copying more lines from
`<base_url>GRCh38-all-stratifications.tsv` into `stratification.regions`.

> **Migration note.** Truth files are now stored as
> `resources/benchmark/<sample>.truth.{vcf.gz,bed}` rather than under their GIAB
> filenames, because remote and local names cannot be assumed equal — HG001
> publishes `..._benchmark.bed` while the whole Ashkenazim trio publishes
> `..._benchmark_noinconsistent.bed`. If you already have the HG001 files, rename
> rather than re-download them:
> ```
> cd resources/benchmark
> mv HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz     NA12878.truth.vcf.gz
> mv HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi NA12878.truth.vcf.gz.tbi
> mv HG001_GRCh38_1_22_v4.2.1_benchmark.bed        NA12878.truth.bed
> ```

### Measured accuracy

Committed under [`benchmarks/`](benchmarks/), one directory per run, each with a
`PROVENANCE.md` recording the commit, config overlay, capture BED, scoring engine and
truth set the numbers came from. Both runs below are GIAB HG001/NA12878 GRCh38 v4.2.1
scored with hap.py 0.3.15.

| Run | Type | Recall (ALL / PASS) | Precision (ALL / PASS) | F1 (ALL / PASS) |
|---|---|---|---|---|
| `default` (Garvan/Nextera) | SNP | 0.8924 / 0.8210 | 0.9737 / 0.9859 | 0.9313 / 0.8959 |
| | INDEL | 0.7643 / 0.7633 | 0.6686 / 0.7402 | 0.7133 / 0.7515 |
| `twist_onso` | SNP | 0.9406 / 0.9076 | 0.9859 / 0.9979 | 0.9627 / 0.9506 |
| | INDEL | 0.8822 / 0.8810 | 0.7120 / 0.8269 | 0.7880 / 0.8531 |

**Do not compare the two runs' recall or F1 to each other.** The evaluation region is
the capture BED intersected with GIAB high confidence, so different kits give
different denominators — 23,713 truth SNPs for `twist_onso` against 46,032 for
`default`. They answer different questions.

What *is* comparable, and what both runs agree on, is the **ALL-to-PASS behaviour:
the default hard filters help indels and cost more than they buy on SNPs.**

| | true variants lost | false positives removed | F1 |
|---|---|---|---|
| SNP, `default` | 3,287 | 568 | 0.9313 → 0.8959 |
| SNP, `twist_onso` | 782 | 273 | 0.9627 → 0.9506 |
| INDEL, `default` | 5 | 550 | 0.7133 → **0.7515** |
| INDEL, `twist_onso` | 1 | 142 | 0.7880 → **0.8531** |

The two effects come from two different filters that never touch the same variant
type: `SOR3` is defined only under `filtering.snvs`, so it cannot fire on an indel,
while `QD2` does the indel work. Measured on the retained `twist_onso` callset, `SOR3`
carries 74.0% of all filtered records and 85.9% of filtered SNPs, and `QD2` carries
90.3% of filtered indels. **The thresholds are
left at GATK best practice and documented, not tuned** — tuning them against the one
sample we benchmark would fit the filter to the test set. If your application is SNP
recall-sensitive, `SOR3` in `config/config.yaml` is the first threshold to
reconsider, and `results/filtered/all.filtered.vcf.gz` keeps every failing record
with its FILTER tag intact so you can re-decide without re-running the pipeline.

Each run's `NOTES.md` records dataset caveats and stage-by-stage runtimes.

### Running a second dataset without destroying the first

`benchmark.label` namespaces the hap.py output: unset (the default) writes
`results/benchmark/happy.*`, and setting it writes `happy.<label>.*`. The
`config/twist_onso.yaml` overlay uses it:

```powershell
.\run --configfile config/twist_onso.yaml --cores 8 all benchmark
```

Put `--cores` immediately after the configfile value — Snakemake 8's `--configfile`
is greedy and will otherwise swallow the target name.

> **Only the hap.py output is namespaced.** `results/genotyped/`,
> `results/filtered/` and `results/annotated/` are cohort-level and share fixed
> paths across configs, and both configs here call the same sample name, so
> `results/bqsr/NA12878.recal.bam` collides too. **Run different configs
> sequentially, and copy `results/annotated` and `results/filtered` aside between
> them if you want to keep both.** Running them concurrently is not supported and
> Snakemake's workdir lock will refuse it anyway.

> **Switching to a config with DIFFERENT sample names needs the cohort results
> cleared**, and this is not optional:
> ```
> rm -rf results/genotyped results/filtered results/annotated
> ```
> Snakemake decides a target is up to date by comparing it against its *direct*
> inputs. When the whole cohort-level chain already exists it never descends as far
> as `combine_gvcfs`, so it never notices the sample set changed — `all` would
> report success while `results/annotated/` still described the previous cohort,
> and the benchmark would then try to pull this cohort's samples out of it.
> `common.smk` raises a `WorkflowError` at DAG-build time rather than letting that
> happen, and the message names the command above.
>
> Only the cohort level needs clearing. `results/called/*.g.vcf.gz` and
> `results/bqsr/*.recal.bam` are not `temp()`, so samples already processed are not
> realigned and returning to a previous config costs the cohort steps only.

### Benchmarking a trio (GIAB Ashkenazim)

`config/giab_trio.yaml` runs HG002 (son), HG003 (father) and HG004 (mother) as one
cohort — the first configuration here with more than one sample, and therefore the
only one that exercises joint genotyping at n>1 or produces a Mendelian violation
rate.

GIAB publishes **no FASTQ** for these exomes. The only Illumina exome data for the
trio is a position-sorted, duplicate-marked BAM per sample, so download those by
hand (≈27.5 GB total; `curl -C -` resumes a broken transfer, which matters at this
size):

```powershell
$b='https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/'
$dst='data\giab_trio'; New-Item -ItemType Directory -Force $dst | Out-Null
@{ HG002='HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam';
   HG003='HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam';
   HG004='HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam' }.GetEnumerator() |
  ForEach-Object { curl.exe -L -C - -o "$dst\$($_.Key).exome.bam" "$b$($_.Value)" }
```

`bam_revert.smk` then reverts each BAM to paired FASTQ under `results/reverted/`,
which is what `config/units.giab_trio.tsv` points at. The source alignment is
discarded entirely, so this pipeline's own mapping, duplicate marking and
recalibration all apply as they would to any other input.

```powershell
rm -r -force results\genotyped, results\filtered, results\annotated   # see above
.\run --configfile config/giab_trio.yaml --cores 8 all benchmark
```

> **The trio is scored against a PROXY capture BED, and its absolute recall is
> depressed as a result.** This data was captured with Agilent SureSelect V5, whose
> target BEDs are available only through a SureDesign account and cannot be
> redistributed, so there is no public equivalent to ship. Twist Exome v2 is used
> instead: a core-coding design that mostly sits inside V5's larger footprint, but
> any region Twist targets that V5 did not enrich becomes an uncovered false
> negative for reasons that have nothing to do with the caller.
>
> That is acceptable for what this run is for. HG002 vs HG003 vs HG004 share a kit,
> a BED and their denominators, so the comparison between them is internally valid,
> and the evaluation region is identical to the `twist_onso` run's. **Do not
> compare these absolute recalls to the default (Nextera) run.** Quantify the
> mismatch rather than asserting it — `SAMPLE=HG002 bash dev/diag_capture.sh` splits
> the false negatives into covered and uncovered, and the `refseq_cds`
> stratification subset gives a design-independent slice for free.

**Mendelian violation rate** is emitted to `results/benchmark/mendelian.giab_trio.txt`
whenever `samples.tsv` carries `sex` / `paternal_id` / `maternal_id` columns and some
sample names both parents. It asks whether the child's genotypes could have been
inherited, which needs **no truth set at all** — so unlike every hap.py number it
covers the whole callset rather than the high-confidence regions, and cannot be
flattered by anything specific to the benchmark. Restricted to chr1–22, matching the
truth sets' own scope.

> **Rebuild the image first.** These features add tool envs (`ucsc-liftover`,
> `hap.py` + `rtg-tools`), so a pre-baked image built before them is stale — the
> launcher's staleness guard will warn. Re-bake with `docker rmi <image>` and then
> re-run the launcher, which
> rebuilds under the tag the launcher actually looks for. Do **not** build by hand
> with `docker build -t ngs-germline-wes ...`: that produces only the `:latest` tag,
> the launcher still won't find the tag it wants, and you pay the 15–30 min build
> twice. A manual build also skips `--build-arg BASE=`, silently ignoring
> `docker/base-image.txt`.

---

## Troubleshooting

- **`Missing input files … twist_exome_v2.GRCh38.bed`** — the default capture BED
  hasn't been built yet; run `run get_capture_bed` (or point `intervals.bed` at your
  own kit-specific BED).
- **OOM / killed JVM** — raise Docker's memory and/or lower the per-step `mem_mb`
  in `config.yaml` so `-Xmx` fits.
- **SnpEff "database not found"** — the `snpeff_download` rule needs network; or set
  `annotation.snpeff.db` to one shown by `snpEff databases | grep -i grch38`.
- **VariantFiltration JEXL warnings about missing annotations** — expected when an
  annotation (e.g. `MQRankSum`) is absent at a site; that site simply isn't filtered
  by that expression.
- **`cohort.vcf.gz exists but was built WITHOUT these samples`** — you switched to a
  config with different sample names. Clear the cohort-level results as the message
  says; per-sample BAMs and GVCFs are kept, so nothing is realigned.
- **Container OOM-killed during `benchmark` or `rtg_format_sdf`, no Java stack
  trace** — RTG sized its heap from the *host's* RAM rather than the container's.
  Both rules export `RTG_MEM`; if you added an rtg rule of your own, export it there
  too.
- **hap.py aborts on a multi-sample VCF** — score the per-sample query
  (`results/benchmark/<sample>.query.vcf.gz`), which `benchmark_query_vcf` produces.
  hap.py has no sample-selection flag.
- **`revert_bam_to_fastq` produces almost only singletons** — the input BAM was not
  collated, so mates were never adjacent. The rule collates first; if you adapted it,
  do not feed `samtools fastq` a position-sorted BAM directly, and never drop `-s`
  (without it singletons are written into the R1/R2 streams and destroy pairing).

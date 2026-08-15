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
│   └── units.twist_onso.tsv    #   its unit sheet (data not publicly available)
├── workflow/
│   ├── Snakefile               # includes the rule modules, defines `all` + `setup_reference`
│   ├── rules/*.smk             # common, resources, capture_bed, qc, trim, map,
│   │                           #   dedup_bqsr, calling, filtering, annotation,
│   │                           #   stats, benchmark, smoke
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
> and takes 60–90 minutes by itself on top of the ~4.9 GB of downloads. It is a
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
> than the largest single rule's reservation — `bwa_mem` reserves 12000 MiB, so
> anything below 16g is rejected up front with a message naming the offending
> rule. That check exists because the failure it replaces was miserable: Snakemake
> only raises the error once the oversized job becomes schedulable, i.e. after
> FastQC and trimming have already run, and the message blames pipes.

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

## Evaluating against GIAB (NA12878 / HG001)

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
   → `results/benchmark/happy.summary.csv` (precision / recall / F1 for SNPs & indels).

### Measured accuracy

Committed under [`benchmarks/`](benchmarks/), one directory per run, each with a
`PROVENANCE.md` recording the commit, config overlay, capture BED, scoring engine and
truth set the numbers came from. Both runs below are GIAB HG001/NA12878 GRCh38 v4.2.1
scored with hap.py 0.3.15.

| Run | Type | Recall (ALL / PASS) | Precision (ALL / PASS) | F1 (ALL / PASS) |
|---|---|---|---|---|
| `default` (Garvan/Nextera) | SNP | 0.8925 / 0.8201 | 0.9715 / 0.9839 | 0.9303 / 0.8946 |
| | INDEL | 0.7643 / 0.7629 | 0.6642 / 0.7365 | 0.7107 / 0.7495 |
| `twist_onso` | SNP | 0.9406 / 0.9076 | 0.9859 / 0.9979 | 0.9627 / 0.9506 |
| | INDEL | 0.8810 / 0.8797 | 0.7111 / 0.8257 | 0.7869 / 0.8519 |

**Do not compare the two runs' recall or F1 to each other.** The evaluation region is
the capture BED intersected with GIAB high confidence, so different kits give
different denominators — 23,713 truth SNPs for `twist_onso` against 46,032 for
`default`. They answer different questions.

What *is* comparable, and what both runs agree on, is the **ALL-to-PASS behaviour:
the default hard filters help indels and cost more than they buy on SNPs.**

| | true variants lost | false positives removed | F1 |
|---|---|---|---|
| SNP, `default` | 3,333 | 588 | 0.9303 → 0.8946 |
| SNP, `twist_onso` | 782 | 273 | 0.9627 → 0.9506 |
| INDEL, `twist_onso` | 1 | 142 | 0.7869 → **0.8519** |

`SOR3` accounts for ~74% of all filtered records in both runs. **The thresholds are
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

> **Rebuild the image first.** These features add two tool envs (`ucsc-liftover`,
> `hap.py`), so the pre-baked image is now stale — the launcher's staleness guard
> will warn. Re-bake with `docker rmi <image>` and then re-run the launcher, which
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

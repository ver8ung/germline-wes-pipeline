# Dataset and interpretation notes -- twist_onso run

Hand-written, unlike `PROVENANCE.md`, which `dev/collect_benchmark_evidence.sh`
regenerates. Keep observations that the collector cannot derive in this file so a
re-collection does not silently drop them.

## What a full re-run showed

This run realigned from FASTQ -- all 8 `fastp` and all 8 `bwa_mem` jobs re-executed
on 2026-08-16, finishing and scoring on 2026-08-18 -- against an archived run scored
2026-08-15 with hap.py's `xcmp` engine. It is an independent regeneration of the
callset, not a re-score of the retained one, and the SNP metrics came back unchanged:

| | archived (xcmp) | this run (vcfeval) |
|---|---|---|
| SNP ALL | TP 22,304 / FP 318 / F1 0.962728 | identical |
| SNP PASS | TP 21,522 / FP 45 / F1 0.950618 | identical |
| INDEL ALL | TP 703 / FP 295 / F1 0.786946 | TP 704 / FP 294 / F1 0.788046 |
| INDEL PASS | TP 702 / FP 153 / F1 0.851866 | TP 703 / FP 152 / F1 0.853060 |

**Scope that precisely.** `TRUTH.TOTAL/TP/FN`, `QUERY.TOTAL/FP/UNK` and recall,
precision and F1 are identical on both SNP rows to every printed digit. What did move
is the false-positive *cause* split: `FP.gt`/`FP.al` go 43/3 -> 44/18 at ALL and
17/1 -> 18/1 at PASS, and on indels 20/63 -> 48/32 and 20/54 -> 48/24. The engine
relabels how the same 318 and 45 false positives are attributed between genotype and
allele errors without changing which records they are. It is not "bit-identical" --
one diff falsifies that.

Net, exactly one indel moves from FP to TP in each filter tier. That is the entire
xcmp -> vcfeval difference on this callset, and it is the control the `default` run
could not provide: there the `@RG LB:` duplicate-marking fix and the engine switch
landed together, so the engine's contribution was inferred. Here the LB fix does not
apply -- these 8 units genuinely are separate libraries, and MarkDuplicates confirms
it by reporting 8 library rows -- so the engine is the only variable, and it moves one
variant.

MarkDuplicates metrics are identical apart from the `Started on:` header, and the
samtools `stats` and `flagstat` files are byte-identical. That demonstrates the scored
callset and the alignment summary statistics matched. It does **not** demonstrate
byte-identical BAMs, which was never tested.

## Library characteristics

From `qc/NA12878.metrics.txt` (whole library, before interval restriction):

| Metric | Value |
|---|---|
| Libraries | 8, one per sequencing unit |
| Read pairs examined | 231,794,186 (463,589,704 reads, 2 x 100 bp) |
| Unmapped | 1,584 (0.00034%) |
| Duplicate pairs | 6,464,426 (2.79%) |

From `qc/NA12878.stats.txt`, measured on the **target-restricted** recal BAM, so its
denominators are smaller by construction: 402,208,343 reads, 40,160,162,277 bases
mapped (cigar), 2.93% duplicates, 99.85% properly paired, 0.108% error rate, insert
222.5. Capture target 36,458,262 bp across 283,942 intervals.

**Depth is an upper bound, and so is the comparison to `default`.** 40.16 Gb mapped
over a 36.46 Mb target gives ~1,101x, reached only if every read were on target; the
`default` equivalent is 10.91 Gb over 62.05 Mb = ~176x. The ratio of those two bounds
is 6.3x, but each assumes 100% on-target and neither QC file records an on-target
fraction, so **the true depth ratio between the runs is unmeasured.** Earlier notes in
this repo quoted "~8x"; that figure was never measured either. `dev/diag_capture.sh`
computes the on-target fraction if the exact number ever matters.

**Duplicate rate is low for the depth** -- 2.79% of pairs, with an estimated library
size in the hundreds of millions. The depth is real molecular coverage rather than PCR
echo, which is the favourable case for call quality.

## The input is not raw sequencer output

Only 1,584 reads out of 463.6 million are unmapped -- 99.99966% map to GRCh38. Real
exome FASTQs do not behave this way. The dataset was almost certainly extracted from
an already-aligned BAM, with unmapped reads discarded upstream.

The recal BAM reports 1,050 unmapped, but it is interval-restricted, so that figure is
a lower bound on the library and should not be quoted as the library total. The
MarkDuplicates count above is the one measured before `-L` applies.

This is not a pipeline defect, but it makes the input cleaner than what a stranger
running their own FASTQs will feed in. Do not read these numbers as a prediction of
performance on raw data: the pipeline never had to cope with unmappable reads here.

## Where SNP calling fails (stratified, SNP PASS)

The headline F1 of 0.9506 averages regimes that differ by 59 recall points.

| Subset | truth N | recall | F1 |
|---|---|---|---|
| notinalldifficultregions | 16,639 | 0.9392 | 0.9685 |
| refseq_cds | 19,639 | 0.9080 | 0.9508 |
| **(headline)** | 23,713 | 0.9076 | 0.9506 |
| alldifficultregions | 7,074 | 0.8332 | 0.9063 |
| segdups | 1,387 | 0.5443 | 0.6965 |
| lowmappabilityall | 763 | 0.3486 | 0.5067 |
| **MHC** | **553** | **0.0000** | -- |

`alldifficultregions` and its complement partition the truth denominator exactly
(7,074 + 16,639 = 23,713) and the false positives exactly too: 39 of the 45 surviving
PASS false positives are in difficult regions, 6 outside. Outside difficult regions
this run emits **6 false-positive SNPs against 15,628 true positives**.

**The MHC fails completely, and it replicates across capture kits.** 553 confident
truth SNPs here, zero recovered, and -- unlike `default`, which still emits 16 PASS
records (4 FP, 12 UNK) against its 1,049 truth SNPs -- twist_onso emits **no PASS
records at all** in the MHC.

The cause is mapping, not filtering, and the ALL rows show it: before any filter is
applied the MHC already has 4 true positives out of 553, from 11 query records. So
549 of 553 were lost at the calling stage; hard filtering only converts the last 4
into 0. The MHC is the most polymorphic region of the human genome and this pipeline
is **not ALT-aware**: reads from a divergent haplotype have nowhere unique to go
against GRCh38's single reference haplotype. Anyone who needs MHC calls needs an
ALT-aware or graph-based alignment strategy; no amount of depth or threshold tuning
will help, and this run has several times the depth of `default` to prove it.

> **Do not add the MHC's share of failures to the `alldifficultregions` share.** They
> overlap only partially: 3,937,174 of the MHC's 4,970,557 bp (79.2%) lie *outside*
> `alldifficultregions`. `segdups` and `lowmappabilityall` are genuine subsets of it;
> the MHC is not.

## Where indel calling fails (stratified, INDEL PASS)

| Subset | truth N | query emitted | recall | F1 |
|---|---|---|---|---|
| notinalldifficultregions | 279 | 330 | 0.9713 | 0.9731 |
| notinAllTandemRepeatsandHomopolymers | 412 | 504 | 0.9466 | 0.9559 |
| SimpleRepeat_homopolymer_7to11 | 106 | 131 | 0.8962 | 0.9018 |
| SimpleRepeat_homopolymer_4to6 | 196 | 324 | 0.8980 | 0.8892 |
| **(headline)** | 798 | 1,399 | 0.8810 | 0.8531 |
| alldifficultregions | 519 | 1,067 | 0.8324 | 0.7943 |
| **SimpleRepeat_homopolymer_ge12** | **76** | **299** | 0.6316 | **0.4760** |

The same shape as `default`: excellent outside repeats, collapsing inside long
homopolymers, and the collapse driven by **over-calling** rather than by missed
calls -- 299 records emitted where 76 truth variants exist.

Resist the obvious cross-kit comparison. The raw over-call ratio is 3.93x here against
3.33x on `default`, which invites "deeper sequencing over-calls more". It does not
survive scrutiny: `QUERY.TOTAL` includes UNK records -- ones falling outside the GIAB
confident BED, which hap.py cannot score at all -- and the UNK share inside this band
differs by kit (155 of 299 = 51.8% here, 1,094 of 2,414 = 45.3% on `default`).
Excluding them the ratios are 1.89x and 1.82x, i.e. indistinguishable. What the two
runs agree on is the *existence* of the homopolymer collapse, not its magnitude.

Note also that `alldifficultregions` partitions the truth denominator (519 + 279 =
798) but **not** the query side (1,067 + 330 = 1,397 against 1,399), so the two
subsets are not a clean split of the callset.

## Hard filtering: helps indels, hurts SNPs

| | true variants lost | false positives removed | ratio | F1 |
|---|---|---|---|---|
| SNP | 782 | 273 | 2.9 : 1 against | 0.9627 -> 0.9506 |
| INDEL | 1 | 142 | 142 : 1 for | 0.7880 -> **0.8531** |

The `default` run, on a different kit at a fraction of the depth:

| | true variants lost | false positives removed | ratio | F1 |
|---|---|---|---|---|
| SNP | 3,287 | 568 | 5.8 : 1 against | 0.9313 -> 0.8959 |
| INDEL | 5 | 550 | 110 : 1 for | 0.7133 -> 0.7515 |

**The direction reproduces; the magnitude does not.** Both runs lose far more true
SNPs than they remove false ones, and both gain on indels. But 2.9:1 against and 5.8:1
against differ by a factor of two, and the runs differ in capture kit, target size,
region-difficulty composition, duplicate rate and the LB fix as well as in depth.
Nothing here isolates depth as the cause.

**The two effects come from two different filters that never touch the same variant
type.** `SOR3` is defined only under `filtering.snvs` in `config/config.yaml`, so it
cannot fire on an indel by construction. Of the 17,193 filtered records in this
callset, 12,727 carry `SOR3` (74.0%); restricted to the 14,819 filtered SNP-class
records that is 85.9%. On the indel side `QD2` tags 2,143 of 2,374 (90.3%), and inside
the evaluation region it is the sole failing tag on every filtered indel. So the SNP
recall loss is a `SOR3` story and the indel gain is a `QD2` story, and they are
independent knobs.

The `default` figure previously recorded here -- 13,312 of 18,077, 73.6% -- was
measured before that run was repeated with the LB fix, and its callset is not
retained, so it cannot be re-derived from the committed evidence. Treat the 74.0%
above as a twist_onso measurement only.

**The defaults are deliberately left at GATK best practice and documented, not
tuned.** Tuning them against GIAB on the only genome we benchmark would be fitting the
filter to the test set. The observation is recorded so a user with a recall-sensitive
application knows which threshold to reconsider first, and knows it is SNP recall --
not indel precision -- that they would be buying.

## These numbers are NOT comparable to the `default` run

The evaluation region is the capture BED intersected with GIAB high confidence, so the
denominators differ between kits. This run scores 23,713 truth SNPs against 46,032 in
`default` -- roughly half the region. Twist's higher scores describe a different,
smaller question, not a better pipeline. Compare ALL-to-PASS behaviour and stratified
*shape* across runs; do not compare recall or F1 directly.

## Runtime, for capacity planning

Single 8-core container, 16 GB cap, Docker Desktop bind mount on a Windows host. Stage
durations are from this run's Snakemake logs.

| Stage | Wall clock |
|---|---|
| fastp (8 units) | 6 min 15 s |
| bwa_mem (8 units) | 93 min 36 s |
| merge_sample_bams | 8 min 38 s |
| mark_duplicates | 52 min 42 s |
| base_recalibrator | 69 min 31 s |
| apply_bqsr | 79 min 36 s |
| haplotype_caller | 115 min 33 s |
| genotype -> filter | ~15 s |
| benchmark_happy (vcfeval + 14 strata) | 18 min 41 s |
| snpeff | 2 min 45 s |

Stages sum to **~7 h 27 min**. Nothing ran concurrently -- SnpEff starts after hap.py
finishes, not alongside it. The run itself was deliberately split across three
sessions with long idle gaps, so elapsed calendar time is not a capacity figure; the
sum above is.

MarkDuplicates is the memory peak: JVM heap reached 4.73 GB against its 7.13 GB
`-Xmx`, comfortably inside the 16 GB container cap. Container RSS was not captured, so
that is a heap figure, not a process figure.

`BaseRecalibrator` writes its ~9 MB recalibration report with many small writes, which
is pathologically slow over a Docker Desktop bind mount (~14 min of its 69 min here,
versus seconds on native storage). Large sequential BAM writes on the same mount
sustain ~40 MB/s, so this affects only the report. Native-storage or `WES_SCRATCH_BAMS`
users will not see it.

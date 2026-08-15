# Dataset and interpretation notes -- twist_onso run

Hand-written, unlike `PROVENANCE.md`, which `dev/collect_benchmark_evidence.sh`
regenerates. Keep observations that the collector cannot derive in this file so a
re-collection does not silently drop them.

## Library characteristics

| Metric | Value |
|---|---|
| Reads | 402,208,343 (2 x 100 bp, insert 222.5) |
| Bases mapped (cigar) | 40,160,162,277 |
| Capture target | 36,458,262 bp across 283,942 intervals |
| Duplicates | 11,767,063 (2.93%) |
| Properly paired | 99.85% |
| Error rate | 0.108% |

**Depth.** 40.16 Gb mapped over a 36.46 Mb target is an upper bound of ~1101x,
which would only be reached if every read were on target. At the 60-80% on-target
typical of exome capture the mean on-target depth is ~660-880x. Quote it as a
range: this run never measured the on-target fraction directly. `dev/diag_capture.sh`
computes it if the exact figure ever matters.

**Duplicate rate is low for the depth** -- 2.93%, with an estimated library size of
532 million. The depth is real molecular coverage rather than PCR echo, which is the
favourable case for call quality. It also means MarkDuplicates removed almost
nothing, so HaplotypeCaller saw essentially the full depth.

## The input is not raw sequencer output

Only 1,050 reads out of 402 million are unmapped -- 99.9997% map to GRCh38. Real
exome FASTQs do not behave this way. The dataset was almost certainly extracted from
an already-aligned BAM, with unmapped reads discarded upstream.

This is not a pipeline defect, but it makes the input cleaner than what a stranger
running their own FASTQs will feed in. Do not read these numbers as a prediction of
performance on raw data: the pipeline never had to cope with unmappable reads here.

## Hard filtering: helps indels, hurts SNPs

Derived from the ALL and PASS rows in `PROVENANCE.md`.

| | true variants lost | false positives removed | ratio | F1 |
|---|---|---|---|---|
| SNP | 782 | 273 | 2.9 : 1 against | 0.9627 -> 0.9506 |
| INDEL | 1 | 142 | 142 : 1 for | 0.7869 -> 0.8519 |

The `default` (Garvan/Nextera) run shows the same asymmetry: SNPs 3,333 lost against
588 removed (5.7:1 against, F1 0.9303 -> 0.8946), indels net positive. Two capture
kits, different targets, an ~8x depth difference, same sign both times.

`SOR3` dominates the filtered set in both runs -- 12,727 of 17,193 filtered records
here (74.0%), 13,312 of 18,077 (73.6%) in `default`. It is the single knob driving
the SNP recall loss.

**The defaults are deliberately left at GATK best practice and documented, not
tuned.** Tuning them against GIAB on the only sample we benchmark would be fitting
the filter to the test set. The observation is recorded so a user with a
recall-sensitive application knows which threshold to reconsider first, and knows it
is SNP recall -- not indel precision -- that they would be buying.

## These numbers are NOT comparable to the `default` run

The evaluation region is the capture BED intersected with GIAB high confidence, so
the denominators differ between kits. This run scores 23,713 truth SNPs against
46,032 in `default` -- roughly half the region. Twist's higher scores describe a
different, smaller question, not a better pipeline. Compare ALL-to-PASS behaviour
across runs; do not compare recall or F1 directly.

## Runtime, for capacity planning

Single 8-core container, 16 GB cap, Docker Desktop bind mount on a Windows host.

| Stage | Wall clock |
|---|---|
| bwa_mem (8 units, run on 14.08) | ~2 h |
| mark_duplicates | 52 min |
| base_recalibrator | 66 min (51.6 min traversal + ~14 min report write) |
| apply_bqsr | ~76 min |
| haplotype_caller | 113.4 min, 1.5 GB peak heap |
| genotype -> annotate | ~5 min |
| benchmark_happy | ~23 min |

Peak container memory was ~4.5 GB during MarkDuplicates, well inside the 16 GB cap,
so this fits comfortably on a 16 GB machine.

`BaseRecalibrator` writes its ~9 MB recalibration report with many small writes,
which is pathologically slow over a Docker Desktop bind mount (~14 min here versus
seconds on native storage). Large sequential BAM writes on the same mount sustain
~40 MB/s, so this affects only the report. Native-storage or `WES_SCRATCH_BAMS`
users will not see it.

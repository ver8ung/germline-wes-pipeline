# Dataset and interpretation notes -- default run (Garvan/Nextera)

Hand-written, unlike `PROVENANCE.md`, which `dev/collect_benchmark_evidence.sh`
regenerates. Keep observations that the collector cannot derive in this file so a
re-collection does not silently drop them.

## Hard filtering: helps indels, hurts SNPs

Derived from the ALL and PASS rows in `PROVENANCE.md`.

| | true variants lost | false positives removed | ratio | F1 |
|---|---|---|---|---|
| SNP | 3,333 | 588 | 5.7 : 1 against | 0.9303 -> 0.8946 |
| INDEL | -6 (recovered) | 562 | net positive | 0.7107 -> 0.7495 |

`SOR3` tags 13,312 of the 18,077 filtered records (73.6%). The `twist_onso` run
reproduces both the asymmetry and the `SOR3` share (74.0%) on a different capture
kit at roughly 8x the depth -- see `benchmarks/twist_onso/NOTES.md`.

Thresholds are deliberately left at GATK best practice and documented, not tuned;
tuning against the only sample we benchmark would fit the filter to the test set.

## These numbers are NOT comparable to the `twist_onso` run

The evaluation region is the capture BED intersected with GIAB high confidence, so
the denominators differ between kits. This run scores 46,032 truth SNPs against
23,713 in `twist_onso`. Compare ALL-to-PASS behaviour across runs; do not compare
recall or F1 directly.

## What this run actually established

The correctness fixes that preceded it (GenotypeGVCFs interval padding, MIXED/MNP
recovery, `bwa -K` determinism, the PASS-only deliverable) **did not move the
accuracy numbers**. Against the archived pre-fix run from June, SNP ALL is identical
to four decimal places (4,947 false negatives both times), SNP PASS differs by two,
and INDEL ALL recovered six variants.

That is the honest result: this was correctness and reproducibility work, not
accuracy work. What the run does establish is **determinism** -- the same numbers
reproduced across a code change, an engine pin, a host restart and a mid-run pause.

## Runtime, for capacity planning

Single 8-core container, 16 GB cap, Docker Desktop bind mount on a Windows host.
Reconstructed from per-rule log timestamps, so gaps include Snakemake scheduling.

| Stage | Wall clock |
|---|---|
| bwa_mem (4 units) | ~28 min |
| base_recalibrator -> apply_bqsr | ~21 min for apply |
| haplotype_caller | 34.8 min, 3.0 GB peak heap |
| genotype -> annotate | ~14 min |
| benchmark_happy | ~17 min |

Peak HaplotypeCaller heap here (3.0 GB) is *higher* than the far deeper `twist_onso`
run (1.5 GB), so heap does not scale with depth the way one might assume.

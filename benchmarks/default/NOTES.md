# Dataset and interpretation notes -- default run (Garvan/Nextera)

Hand-written, unlike `PROVENANCE.md`, which `dev/collect_benchmark_evidence.sh`
regenerates. Keep observations that the collector cannot derive in this file so a
re-collection does not silently drop them.

## What changed since the previous (archived) numbers

Two things changed at once, and they can be separated cleanly.

| | archived (xcmp, LB per lane) | this run (vcfeval, LB fixed) |
|---|---|---|
| SNP ALL F1 | 0.9303 | 0.9313 |
| SNP PASS F1 | 0.8946 | 0.8959 |
| INDEL ALL F1 | 0.7107 | 0.7133 |
| INDEL PASS F1 | 0.7495 | 0.7515 |

F1 barely moves, but that understates what happened. **The change is in false
positives, not recall:**

| | FP before | FP after | change |
|---|---|---|---|
| SNP ALL | 1,207 | 1,108 | -8.2% |
| SNP PASS | 619 | 540 | -12.8% |
| INDEL ALL | 1,914 | 1,877 | -1.9% |
| INDEL PASS | 1,352 | 1,327 | -1.8% |

True positives are static (SNP ALL 41,085 -> 41,080; SNP PASS 37,752 -> 37,793).

**The SNP gain is entirely the `LB:` fix, not the engine.** The `twist_onso` run --
which does *not* receive the LB fix, because its 8 units genuinely are separate
libraries -- was later re-run from FASTQ and scored with vcfeval. Its SNP
`TRUTH.TOTAL/TP/FN`, `QUERY.TOTAL/FP/UNK`, recall, precision and F1 came back
identical to the archived xcmp values to every printed digit; only the false-positive
cause split (`FP.gt`/`FP.al`) moved, and on indels exactly one variant crossed from FP
to TP. So the engine is worth about one indel and nothing at all on SNPs, and the
whole improvement here is attributable to duplicate marking. See
`benchmarks/twist_onso/NOTES.md` for that comparison in full.

`@RG LB:` used to be fabricated as `{sample}.{unit}`, making each lane look like its
own library. This sheet is two libraries over two lanes each, and MarkDuplicates
keys on (library, 5' position), so cross-lane duplicates were never marked:

| | libraries | duplicate pairs | rate |
|---|---|---|---|
| before | 4 (one per lane) | 4,749,702 | ~6.0% |
| after | 2 (NIST7035, NIST7086) | 8,364,400 | ~10.5% |

76% more duplicates detected. Mechanistically the precision gain is what you would
predict: unmarked PCR duplicates supply spurious *independent* support for an
artifact, inflating confidence in a wrong call, while true variants already have
ample non-duplicate support. Hence fewer false positives at unchanged recall.

## Where calling actually fails (stratified, SNP PASS)

The headline F1 of 0.8959 is an average over wildly different regimes. From
`happy.extended.csv`:

| Subset | truth N | recall | F1 |
|---|---|---|---|
| refseq_cds | 18,626 | 0.8625 | 0.9214 |
| notinalldifficultregions | 32,627 | 0.8478 | 0.9135 |
| **(headline)** | 46,032 | 0.8210 | 0.8959 |
| alldifficultregions | 13,405 | 0.7559 | 0.8512 |
| segdups | 3,859 | 0.5955 | 0.7291 |
| lowmappabilityall | 1,940 | 0.3969 | 0.5488 |
| **MHC** | **1,049** | **0.0000** | **--** |

**The MHC is a total failure, and it replicates across capture kits.** 1,049
confident truth SNPs in the evaluation region, **zero** recovered, 16 query records
emitted. The `twist_onso` run shows the same thing on a completely different capture
design and at several times the depth -- 553 truth SNPs, zero recovered, and not one
PASS record emitted -- so this is a property of the pipeline rather than of either
kit. (Its 11 MHC records are the ALL row; all 11 fail a hard filter. Do not quote the
ALL record count next to a PASS true-positive count, as an earlier version of this
sentence did.)

The cause is mapping, not capture or filtering. Measured on the twist_onso BAM,
**98.2% of reads in chr6:28.5-33.5 Mb carry MAPQ 0**, against 0.07% in a control
region on the same chromosome -- so the reads are captured and present, but bwa
cannot place them uniquely and HaplotypeCaller discards them. The MHC is the most
polymorphic region of the human genome and this pipeline is **not ALT-aware**: reads
from a divergent haplotype have nowhere unique to go against GRCh38's single
reference haplotype. Anyone who needs MHC calls needs an ALT-aware or graph-based
alignment strategy; no amount of depth or filter tuning will help.

## Where indel calling fails (stratified, INDEL PASS)

| Subset | truth N | query emitted | recall | F1 |
|---|---|---|---|---|
| notinalldifficultregions | 1,805 | 1,875 | 0.9091 | 0.9464 |
| notinAllTandemRepeatsandHomopolymers | 2,395 | 2,524 | 0.8894 | 0.9311 |
| **(headline)** | 4,871 | 7,389 | 0.7633 | 0.7515 |
| alldifficultregions | 3,075 | 5,524 | 0.6780 | 0.6488 |
| AllTandemRepeats_le50bp | 530 | 674 | 0.5868 | 0.6411 |
| **SimpleRepeat_homopolymer_ge12** | **724** | **2,414** | 0.4599 | **0.3454** |

Indel calling is excellent outside repeats (F1 0.946) and collapses inside long
homopolymers, where it emits **3.3x more records than there are truth variants**.
That over-calling, not missed calls, is what drags the headline indel F1 down: the
`*` row's 7,389 query records against 4,871 truth variants is mostly repeat-region
noise. This is a known, expected weakness of short reads in homopolymers rather than
anything specific to this pipeline -- but it is now measured rather than assumed,
and it says clearly which indels from this callset to distrust.

## Hard filtering: helps indels, hurts SNPs

| | true variants lost | false positives removed | ratio | F1 |
|---|---|---|---|---|
| SNP | 3,287 | 568 | 5.8 : 1 against | 0.9313 -> 0.8959 |
| INDEL | 5 | 550 | 110 : 1 for | 0.7133 -> 0.7515 |

Unchanged in direction from the previous run and reproduced by `twist_onso` on a
different kit at substantially greater depth -- though the depth ratio between the two
runs is an unmeasured quantity, not the "~8x" this file used to assert: both figures
are mapped-bases-over-target-size upper bounds that assume a 100% on-target rate
neither run measured. The direction reproduces; the magnitude does not (5.8:1 against
here, 2.9:1 there). `SOR3` is the dominant tag in this run's filtered set, and in
`twist_onso` it carries 74.0% of all filtered records and 85.9% of filtered SNPs.
Thresholds are deliberately left at GATK best practice and documented, not tuned:
tuning against the only sample we benchmark would be fitting the filter to the test
set.

Note the stratified view sharpens this. Low-mappability SNP recall falls from 0.756
(ALL) to 0.397 (PASS) -- the filters' cost is concentrated exactly where mapping is
already marginal, which is also where their benefit is largest. A recall-sensitive
application should revisit `SOR3` first, and can now see precisely which regions it
would be trading.

## These numbers are NOT comparable to the other runs

The evaluation region is the capture BED intersected with GIAB high confidence, so
denominators differ between kits and truth sets. This run scores 46,032 truth SNPs
against 23,713 for `twist_onso`. Compare ALL-to-PASS behaviour and stratified
*shape* across runs; do not compare absolute recall or F1.

## Runtime, for capacity planning

Single 8-core container, 16 GB cap, Docker Desktop bind mount on a Windows host.
Wall clock 12:37 -> 15:24, **2 h 47 min** end to end including scoring.

| Stage | Wall clock |
|---|---|
| fastp (4 units) | ~2 min |
| bwa_mem (4 units) | ~10 min |
| merge + mark_duplicates | ~45 min |
| base_recalibrator | ~26 min |
| apply_bqsr | ~21 min |
| haplotype_caller | ~35 min |
| genotype -> filter -> annotate | ~10 min |
| benchmark_happy (vcfeval + 14 strata) | ~23 min |

`rtg format` (the vcfeval SDF, ~1.2 GB) is a one-off ~1 min and is shared by every
subsequent run. Peak container memory stayed around 4.5 GB against the 16 GB cap.

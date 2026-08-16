#!/usr/bin/env bash
# Callable-region diagnosis: is a recall deficit capture-driven (truth variants
# sitting in target regions the capture never covered) or pipeline-driven
# (covered regions that were mis-called, under-called or over-filtered)?
#
# It derives the actually-callable region from read depth in the recalibrated
# BAM, then attributes every hap.py false negative to "uncovered" or "covered".
# A high uncovered fraction means the caller is at its ceiling and more depth or
# a different BED will not help; a high covered fraction points at the pipeline.
#
# THIS IS NOT PART OF THE PIPELINE. It is an ad-hoc analysis script that runs
# inside the project's container, against results a run has already produced.
#
# Usage (from the repo root, after a run + `benchmark` have completed):
#
#   docker run --rm -v "$PWD:/workflow" -w /workflow ngs-germline-wes:latest \
#       bash dev/diag_capture.sh
#
#   # for a labelled run (see config/twist_onso.yaml):
#   docker run --rm -v "$PWD:/workflow" -w /workflow ngs-germline-wes:latest \
#       bash -c 'LABEL=twist_onso BED=resources/intervals/twist_exome_v2.GRCh38.bed \
#                bash dev/diag_capture.sh'
#
#   # for one sample of a multi-sample run (see config/giab_trio.yaml):
#   docker run --rm -v "$PWD:/workflow" -w /workflow ngs-germline-wes:latest \
#       bash -c 'LABEL=giab_trio SAMPLE=HG002 \
#                BED=resources/intervals/twist_exome_v2.GRCh38.bed \
#                bash dev/diag_capture.sh'
#
# The trio is the case this matters most for: it is scored against a PROXY capture
# BED (Twist v2 standing in for the unobtainable Agilent SureSelect V5 targets), so
# some of its false negatives are certainly regions the real kit never enriched.
# FRAC_CALLABLE_10x below is what turns that caveat into a number.
#
# If the run used the WES_SCRATCH_BAMS scratch volumes, results/bqsr lives on a
# Docker volume rather than the bind mount, so add:
#   --mount type=volume,source=wes-bqsr,target=/workflow/results/bqsr
set -euo pipefail

# LABEL selects which hap.py output to read. Empty (the default) = the
# unlabelled output of the default config; otherwise happy.<LABEL>.*
LABEL="${LABEL:-}"
# SAMPLE picks which sample of the cohort to diagnose. hap.py output, the recal
# BAM and the truth BED are ALL per-sample now, so deriving the three of them from
# one variable is what keeps them describing the same sample -- overriding only
# BAM and leaving the truth BED pointing at another genome would silently compare
# one sample's coverage against another's confident regions.
SAMPLE="${SAMPLE:-NA12878}"
BED="${BED:-resources/intervals/nextera_expandedexome.GRCh38.bed}"
BAM="${BAM:-results/bqsr/${SAMPLE}.recal.bam}"
CONF="${CONF:-resources/benchmark/${SAMPLE}.truth.bed}"

if [ -n "$LABEL" ]; then
  HAPPY="results/benchmark/happy.${LABEL}.${SAMPLE}.vcf.gz"
else
  HAPPY="results/benchmark/happy.${SAMPLE}.vcf.gz"
fi

SAMTOOLS=$(ls /opt/snakemake-envs/*/bin/samtools 2>/dev/null | head -1 || true)
BEDTOOLS=$(ls /opt/snakemake-envs/*/bin/bedtools 2>/dev/null | head -1 || true)
echo "LABEL=${LABEL:-<none>}"
echo "SAMPLE=$SAMPLE"
echo "BED=$BED"
echo "BAM=$BAM"
echo "CONF=$CONF"
echo "HAPPY=$HAPPY"
echo "TOOL_samtools=$SAMTOOLS"
echo "TOOL_bedtools=$BEDTOOLS"
[ -n "$SAMTOOLS" ] || { echo "FATAL: samtools not in baked envs"; exit 3; }
[ -n "$BEDTOOLS" ] || { echo "FATAL: bedtools not in baked envs"; exit 3; }
[ -f "$BED"   ] || { echo "FATAL: $BED not found (run get_capture_bed)"; exit 3; }
[ -f "$BAM"   ] || { echo "FATAL: $BAM not found (if WES_SCRATCH_BAMS was used, mount the wes-bqsr volume at /workflow/results/bqsr)"; exit 3; }
[ -f "$HAPPY" ] || { echo "FATAL: $HAPPY not found (run the 'benchmark' target first; set LABEL= if the run was labelled)"; exit 3; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
bp() { awk '{s+=$3-$2} END{printf "%d", s+0}' "$1"; }

# --- eval region = capture targets ∩ GIAB confident (what hap.py scored over) ---
sort -k1,1 -k2,2n "$BED"  > "$TMP/bed.sorted"
sort -k1,1 -k2,2n "$CONF" > "$TMP/conf.sorted"
"$BEDTOOLS" intersect -a "$TMP/bed.sorted" -b "$TMP/conf.sorted" \
  | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -i - > "$TMP/eval.bed"
BED_BP=$(bp "$TMP/bed.sorted"); EVAL_BP=$(bp "$TMP/eval.bed")
echo "TARGET_BP=$BED_BP"
echo "EVAL_BP=$EVAL_BP"

# --- callable bp within eval from actual read depth (MQ>=20, exclude dups) ---
"$SAMTOOLS" depth -a -b "$TMP/eval.bed" -Q 20 -G 1024 "$BAM" \
  | awk -v B="$BEDTOOLS" -v T="$TMP" '
      { if($3>=10){ print $1"\t"($2-1)"\t"$2 | (B" merge -i - > "T"/callable10.bed") }
        if($3>=20){ print $1"\t"($2-1)"\t"$2 | (B" merge -i - > "T"/callable20.bed") } }'
CALL10_BP=$(bp "$TMP/callable10.bed"); CALL20_BP=$(bp "$TMP/callable20.bed")
echo "CALLABLE10_BP=$CALL10_BP"
echo "CALLABLE20_BP=$CALL20_BP"
awk -v e="$EVAL_BP" -v c="$CALL10_BP" 'BEGIN{printf "FRAC_CALLABLE_10x=%.4f\n",(e>0?c/e:0)}'
awk -v e="$EVAL_BP" -v c="$CALL20_BP" 'BEGIN{printf "FRAC_CALLABLE_20x=%.4f\n",(e>0?c/e:0)}'

# --- FN truth variants from hap.py per-variant VCF, split by type ---
echo "HAPPY_SAMPLES: $(zcat "$HAPPY" | grep -m1 '^#CHROM' | cut -f10-)"
zcat "$HAPPY" | awk -F'\t' '
  /^#CHROM/ { for(i=10;i<=NF;i++){ if($i=="TRUTH")tc=i }; next }
  /^#/ { next }
  { n=split($9,fmt,":"); bd=0; bvt=0;
    for(i=1;i<=n;i++){ if(fmt[i]=="BD")bd=i; if(fmt[i]=="BVT")bvt=i }
    if(tc==0) next
    split($tc,Tr,":"); dec=(bd>0?Tr[bd]:"."); typ=(bvt>0?Tr[bvt]:".")
    if(dec=="FN"){ print $1"\t"($2-1)"\t"$2"\t"typ } }' \
  | sort -k1,1 -k2,2n > "$TMP/fn.bed"

awk '$4=="SNP"{print $1"\t"$2"\t"$3}'  "$TMP/fn.bed" | sort -k1,1 -k2,2n > "$TMP/fn_snv.bed"
awk '$4!="SNP"{print $1"\t"$2"\t"$3}'  "$TMP/fn.bed" | sort -k1,1 -k2,2n > "$TMP/fn_indel.bed"
FN_TOTAL=$(wc -l < "$TMP/fn.bed"); FN_SNV=$(wc -l < "$TMP/fn_snv.bed"); FN_INDEL=$(wc -l < "$TMP/fn_indel.bed")
echo "FN_TOTAL=$FN_TOTAL"
echo "FN_SNV=$FN_SNV"
echo "FN_INDEL=$FN_INDEL"

isect_u() { "$BEDTOOLS" intersect -a "$1" -b "$2" -u | wc -l; }
FN_SNV_C10=$(isect_u "$TMP/fn_snv.bed" "$TMP/callable10.bed")
FN_SNV_C20=$(isect_u "$TMP/fn_snv.bed" "$TMP/callable20.bed")
FN_IND_C10=$(isect_u "$TMP/fn_indel.bed" "$TMP/callable10.bed")
echo "FN_SNV_IN_CALLABLE10=$FN_SNV_C10"
echo "FN_SNV_IN_CALLABLE20=$FN_SNV_C20"
echo "FN_SNV_UNCOVERED10=$((FN_SNV-FN_SNV_C10))"
echo "FN_INDEL_IN_CALLABLE10=$FN_IND_C10"
echo "FN_INDEL_UNCOVERED10=$((FN_INDEL-FN_IND_C10))"
echo "DONE"

#!/usr/bin/env bash
# Ad-hoc diagnosis: is the NA12878 (Twist-Onso) recall loss capture-driven
# (truth variants in uncovered target regions) or pipeline-driven (covered
# regions where we mis-/under-called or over-filtered)?
set -euo pipefail
cd /workflow

BAM=results/bqsr/NA12878.recal.bam
HAPPY=results/benchmark/happy.vcf.gz
TWIST=resources/intervals/twist_exome_v2.GRCh38.bed
CONF=resources/benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed

SAMTOOLS=$(ls /opt/snakemake-envs/*/bin/samtools 2>/dev/null | head -1 || true)
BEDTOOLS=$(ls /opt/snakemake-envs/*/bin/bedtools 2>/dev/null | head -1 || true)
echo "TOOL_samtools=$SAMTOOLS"
echo "TOOL_bedtools=$BEDTOOLS"
[ -n "$SAMTOOLS" ] || { echo "FATAL: samtools not in baked envs"; exit 3; }
[ -n "$BEDTOOLS" ] || { echo "FATAL: bedtools not in baked envs"; exit 3; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
bp() { awk '{s+=$3-$2} END{printf "%d", s+0}' "$1"; }

# --- eval region = Twist targets ∩ GIAB confident (what hap.py scored over) ---
sort -k1,1 -k2,2n "$TWIST" > "$TMP/twist.sorted"
sort -k1,1 -k2,2n "$CONF"  > "$TMP/conf.sorted"
"$BEDTOOLS" intersect -a "$TMP/twist.sorted" -b "$TMP/conf.sorted" \
  | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -i - > "$TMP/eval.bed"
TWIST_BP=$(bp "$TMP/twist.sorted"); EVAL_BP=$(bp "$TMP/eval.bed")
echo "TWIST_BP=$TWIST_BP"
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

#!/bin/bash

set -euo pipefail

# Function to check command success
check_success() {
  if [ $? -ne 0 ]; then
    echo "❌ Error in step: $1"
    exit 1
  else
    echo "✅ Completed: $1"
  fi
}

# Set input and reference files
FASTQ1="sample_R1.fastq.gz"
FASTQ2="sample_R2.fastq.gz"
REF="hg38.fa"
KNOWN_SITES="dbsnp.vcf.gz"
SAMPLE="sample"

# Step 1: Quality Control
fastqc $FASTQ1 $FASTQ2 -o qc_reports/
check_success "Quality Control (FastQC)"

# Step 2: Trimming
trimmomatic PE -threads 4 \
  $FASTQ1 $FASTQ2 \
  ${SAMPLE}_R1_paired.fq.gz ${SAMPLE}_R1_unpaired.fq.gz \
  ${SAMPLE}_R2_paired.fq.gz ${SAMPLE}_R2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
check_success "Trimming (Trimmomatic)"

# Step 3: Alignment
bwa mem -t 4 $REF ${SAMPLE}_R1_paired.fq.gz ${SAMPLE}_R2_paired.fq.gz > ${SAMPLE}.sam
check_success "Alignment (BWA-MEM)"

# Step 4: Convert, Sort, and Index
samtools view -Sb ${SAMPLE}.sam | samtools sort -o ${SAMPLE}_sorted.bam
check_success "SAM to sorted BAM (SAMtools)"
samtools index ${SAMPLE}_sorted.bam
check_success "Indexing BAM (SAMtools)"

# Step 5: Mark Duplicates
picard MarkDuplicates I=${SAMPLE}_sorted.bam O=${SAMPLE}_dedup.bam M=${SAMPLE}_metrics.txt
check_success "Mark Duplicates (Picard)"
samtools index ${SAMPLE}_dedup.bam
check_success "Indexing Deduplicated BAM"

# Step 6: Base Quality Score Recalibration
gatk BaseRecalibrator \
  -I ${SAMPLE}_dedup.bam \
  -R $REF \
  --known-sites $KNOWN_SITES \
  -O ${SAMPLE}_recal_data.table
check_success "BaseRecalibrator (GATK)"

gatk ApplyBQSR \
  -R $REF \
  -I ${SAMPLE}_dedup.bam \
  --bqsr-recal-file ${SAMPLE}_recal_data.table \
  -O ${SAMPLE}_recal.bam
check_success "ApplyBQSR (GATK)"

# Step 7: Variant Calling
gatk HaplotypeCaller \
  -R $REF \
  -I ${SAMPLE}_recal.bam \
  -O ${SAMPLE}_raw.vcf.gz
check_success "Variant Calling (GATK HaplotypeCaller)"

# Step 8: Variant Filtering
gatk VariantFiltration \
  -R $REF \
  -V ${SAMPLE}_raw.vcf.gz \
  -O ${SAMPLE}_filtered.vcf.gz \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_snp_filter"
check_success "Variant Filtering (GATK)"

# Step 9: Annotation
table_annovar.pl ${SAMPLE}_filtered.vcf.gz humandb/ \
  -buildver hg38 \
  -out ${SAMPLE}_annotated \
  -remove \
  -protocol refGene,dbnsfp42a \
  -operation g,f \
  -nastring . \
  -vcfinput
check_success "Variant Annotation (ANNOVAR)"

echo "🎉 Pipeline completed successfully!"

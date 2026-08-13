# =============================================================================
# stats.smk — alignment + variant QC metrics, aggregated into one MultiQC report
# that spans every stage (FastQC, fastp, MarkDuplicates, samtools, bcftools, SnpEff).
# =============================================================================


rule samtools_stats:
    input:
        bam="results/bqsr/{sample}.recal.bam",
        bai="results/bqsr/{sample}.recal.bam.bai",
    output:
        "results/qc/samtools/{sample}.stats.txt",
    log:
        "logs/stats/samtools_stats_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools stats {input.bam} > {output} 2> {log}"


rule samtools_flagstat:
    input:
        bam="results/bqsr/{sample}.recal.bam",
        bai="results/bqsr/{sample}.recal.bam.bai",
    output:
        "results/qc/samtools/{sample}.flagstat.txt",
    log:
        "logs/stats/samtools_flagstat_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools flagstat {input.bam} > {output} 2> {log}"


rule bcftools_stats:
    input:
        # The PASS callset, matching what annotation ships. Run against the
        # flagged VCF this counted records the pipeline had marked as failing,
        # so the MultiQC variant totals disagreed with the actual deliverable.
        vcf="results/filtered/all.pass.vcf.gz",
        tbi="results/filtered/all.pass.vcf.gz.tbi",
    output:
        "results/qc/bcftools/all.stats.txt",
    log:
        "logs/stats/bcftools_stats.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "bcftools stats {input.vcf} > {output} 2> {log}"


def multiqc_inputs(wildcards):
    """Every per-stage QC artifact MultiQC should aggregate."""
    files = []
    for s, u in units.index:
        files += [
            f"results/qc/fastqc/{s}.{u}.R1_fastqc.zip",
            f"results/qc/fastqc/{s}.{u}.R2_fastqc.zip",
            f"results/qc/fastp/{s}.{u}.fastp.json",
        ]
    for s in samples["sample"]:
        files += [
            f"results/qc/dedup/{s}.metrics.txt",
            f"results/qc/samtools/{s}.stats.txt",
            f"results/qc/samtools/{s}.flagstat.txt",
        ]
    files += [
        "results/qc/bcftools/all.stats.txt",
        "results/qc/snpeff/all.snpeff.csv",
    ]
    return files


rule multiqc:
    input:
        multiqc_inputs,
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/stats/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        # Scan the qc + logs trees; --force overwrites a stale report.
        "multiqc --force -o results/qc -n multiqc_report.html results/qc logs 2> {log}"

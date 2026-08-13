# =============================================================================
# dedup_bqsr.smk — mark PCR/optical duplicates, then recalibrate base qualities
# (BQSR) against known variation. Output: analysis-ready, per-sample BAMs.
# =============================================================================


rule mark_duplicates:
    input:
        "results/mapped/{sample}.merged.bam",
    output:
        bam=temp("results/dedup/{sample}.dedup.bam"),
        bai=temp("results/dedup/{sample}.dedup.bam.bai"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    params:
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["markdup"]["mem_mb"],
    log:
        "logs/markdup/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        ( gatk --java-options "-Xmx{params.xmx}m" MarkDuplicates \
              -I {input} -O {output.bam} -M {output.metrics}
          samtools index {output.bam} ) 2> {log}
        """


rule base_recalibrator:
    input:
        bam="results/dedup/{sample}.dedup.bam",
        bai="results/dedup/{sample}.dedup.bam.bai",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
        known=known_sites_files(),
        intervals=INTERVALS,
    output:
        "results/bqsr/{sample}.recal.table",
    params:
        known_args=known_sites_args,
        padding=PADDING,
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["bqsr"]["mem_mb"],
    log:
        "logs/bqsr/{sample}.recal_table.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "-Xmx{params.xmx}m" BaseRecalibrator \
            -I {input.bam} -R {input.fasta} \
            {params.known_args} \
            -L {input.intervals} --interval-padding {params.padding} \
            -O {output} 2> {log}
        """


rule apply_bqsr:
    input:
        bam="results/dedup/{sample}.dedup.bam",
        bai="results/dedup/{sample}.dedup.bam.bai",
        recal="results/bqsr/{sample}.recal.table",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
        intervals=INTERVALS,
    output:
        bam="results/bqsr/{sample}.recal.bam",
        bai="results/bqsr/{sample}.recal.bam.bai",
    params:
        padding=PADDING,
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["bqsr"]["mem_mb"],
    log:
        "logs/bqsr/{sample}.apply.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        # -L restricts the WRITTEN BAM to capture targets (+padding), matching the
        # downstream HaplotypeCaller window. Off-target reads are dropped here
        # instead of being recalibrated and written genome-wide — this is the big
        # I/O win: the persistent recal.bam shrinks from a full-genome BAM to just
        # the on/near-target reads HC actually consumes (no calls change).
        r"""
        ( gatk --java-options "-Xmx{params.xmx}m" ApplyBQSR \
              -I {input.bam} -R {input.fasta} \
              --bqsr-recal-file {input.recal} \
              -L {input.intervals} --interval-padding {params.padding} \
              -O {output.bam}
          samtools index {output.bam} ) 2> {log}
        """

# =============================================================================
# calling.smk — per-sample variant discovery in GVCF mode, then cohort-level
# joint genotyping. CombineGVCFs scales to any n (incl. n=1); for very large
# cohorts swap in GenomicsDBImport (see README).
# =============================================================================


rule haplotype_caller:
    input:
        bam="results/bqsr/{sample}.recal.bam",
        bai="results/bqsr/{sample}.recal.bam.bai",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
        intervals=INTERVALS,
    output:
        gvcf="results/called/{sample}.g.vcf.gz",
        tbi="results/called/{sample}.g.vcf.gz.tbi",
    params:
        padding=PADDING,
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    threads: config["resources"]["haplotypecaller"]["threads"]
    resources:
        mem_mb=config["resources"]["haplotypecaller"]["mem_mb"],
    log:
        "logs/haplotypecaller/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "-Xmx{params.xmx}m" HaplotypeCaller \
            -I {input.bam} -R {input.fasta} \
            -ERC GVCF \
            -L {input.intervals} --interval-padding {params.padding} \
            --native-pair-hmm-threads {threads} \
            -O {output.gvcf} 2> {log}
        """


rule combine_gvcfs:
    input:
        gvcfs=all_sample_gvcfs(),
        tbis=[g + ".tbi" for g in all_sample_gvcfs()],
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
    output:
        gvcf="results/genotyped/cohort.g.vcf.gz",
        tbi="results/genotyped/cohort.g.vcf.gz.tbi",
    params:
        variant_args=gvcf_combine_args,
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["genotype"]["mem_mb"],
    log:
        "logs/genotype/combine.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "-Xmx{params.xmx}m" CombineGVCFs \
            -R {input.fasta} {params.variant_args} \
            -O {output.gvcf} 2> {log}
        """


rule genotype_gvcfs:
    input:
        gvcf="results/genotyped/cohort.g.vcf.gz",
        tbi="results/genotyped/cohort.g.vcf.gz.tbi",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
        intervals=INTERVALS,
    output:
        vcf="results/genotyped/cohort.vcf.gz",
        tbi="results/genotyped/cohort.vcf.gz.tbi",
    params:
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["genotype"]["mem_mb"],
    log:
        "logs/genotype/genotype.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "-Xmx{params.xmx}m" GenotypeGVCFs \
            -R {input.fasta} -V {input.gvcf} \
            -L {input.intervals} \
            -O {output.vcf} 2> {log}
        """

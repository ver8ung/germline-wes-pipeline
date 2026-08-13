# =============================================================================
# filtering.smk — GATK hard-filtering. SNPs and indels use different thresholds,
# so split by type, flag each with VariantFiltration, then merge back. Records
# that fail get a FILTER tag (kept, not dropped) so nothing is silently lost.
# =============================================================================


def filter_args(vartype):
    """Build `--filter-name N --filter-expression "E"` pairs from config.

    Expressions contain '<'/'>' and spaces, so they stay double-quoted to
    survive the bash shell that Snakemake invokes.
    """
    spec = config["filtering"][vartype]
    return " ".join(
        f'--filter-name "{name}" --filter-expression "{expr}"'
        for name, expr in spec.items()
    )


rule select_variants:
    input:
        vcf="results/genotyped/cohort.vcf.gz",
        tbi="results/genotyped/cohort.vcf.gz.tbi",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
    output:
        vcf=temp("results/filtered/cohort.{vartype}.vcf.gz"),
        tbi=temp("results/filtered/cohort.{vartype}.vcf.gz.tbi"),
    wildcard_constraints:
        vartype="snvs|indels",
    params:
        # Every variant class GenotypeGVCFs can emit must land in exactly one of
        # these two branches, or it silently disappears at merge_filtered: a class
        # matching neither selection is written to neither temp VCF, both rules
        # still succeed, and nothing warns. MIXED (multiallelic sites carrying both
        # a SNP and an indel allele) and MNP were being lost that way.
        # MIXED goes with indels because indel thresholds are the safe ones to
        # apply to a record that contains an indel allele; MNP goes with SNPs.
        select=lambda wc: (
            "--select-type-to-include SNP --select-type-to-include MNP"
            if wc.vartype == "snvs"
            else "--select-type-to-include INDEL --select-type-to-include MIXED"
        ),
    log:
        "logs/filter/select_{vartype}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk SelectVariants \
            -R {input.fasta} -V {input.vcf} \
            {params.select} \
            -O {output.vcf} 2> {log}
        """


rule hard_filter:
    input:
        vcf="results/filtered/cohort.{vartype}.vcf.gz",
        tbi="results/filtered/cohort.{vartype}.vcf.gz.tbi",
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        dict=ref_dict(),
    output:
        vcf=temp("results/filtered/cohort.{vartype}.filtered.vcf.gz"),
        tbi=temp("results/filtered/cohort.{vartype}.filtered.vcf.gz.tbi"),
    wildcard_constraints:
        vartype="snvs|indels",
    params:
        filters=lambda wc: filter_args(wc.vartype),
    log:
        "logs/filter/hardfilter_{vartype}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk VariantFiltration \
            -R {input.fasta} -V {input.vcf} \
            {params.filters} \
            -O {output.vcf} 2> {log}
        """


rule merge_filtered:
    input:
        snvs="results/filtered/cohort.snvs.filtered.vcf.gz",
        snvs_tbi="results/filtered/cohort.snvs.filtered.vcf.gz.tbi",
        indels="results/filtered/cohort.indels.filtered.vcf.gz",
        indels_tbi="results/filtered/cohort.indels.filtered.vcf.gz.tbi",
    output:
        vcf="results/filtered/all.filtered.vcf.gz",
        tbi="results/filtered/all.filtered.vcf.gz.tbi",
    log:
        "logs/filter/merge.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk MergeVcfs \
            -I {input.snvs} -I {input.indels} \
            -O {output.vcf} 2> {log}
        """

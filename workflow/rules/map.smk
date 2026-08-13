# =============================================================================
# map.smk — align each unit with bwa-mem (read group attached), coordinate-sort,
# then merge a sample's units into a single BAM for duplicate marking.
# =============================================================================


rule bwa_mem:
    input:
        r1="results/trimmed/{sample}.{unit}.R1.fastq.gz",
        r2="results/trimmed/{sample}.{unit}.R2.fastq.gz",
        fasta=config["ref"]["fasta"],
        idx=multiext(config["ref"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("results/mapped/{sample}.{unit}.sorted.bam"),
    params:
        # '\t' escapes are interpreted by bwa itself; keep them literal here.
        rg=get_read_group,
    threads: config["resources"]["bwa_mem"]["threads"]
    resources:
        mem_mb=config["resources"]["bwa_mem"]["mem_mb"],
    log:
        "logs/bwa_mem/{sample}.{unit}.log",
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
        ( bwa mem -M -t {threads} -R '{params.rg}' \
              {input.fasta} {input.r1} {input.r2} \
          | samtools sort -@ {threads} -o {output} - ) 2> {log}
        """


rule merge_sample_bams:
    input:
        get_sample_unit_bams,
    output:
        temp("results/mapped/{sample}.merged.bam"),
    threads: config["resources"]["sort"]["threads"]
    log:
        "logs/merge/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        # Single-unit samples are the common case; skip merge (some samtools
        # builds reject a 1-input merge) and just copy. Multi-unit -> real merge.
        r"""
        ( n=$(echo {input} | wc -w)
          if [ "$n" -eq 1 ]; then
              cp {input} {output}
          else
              samtools merge -@ {threads} -f {output} {input}
          fi ) 2> {log}
        """

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
        # -K fixes the batch size bwa processes at once. Without it the batch is
        # 10000000 * nThreads, so the per-batch insert-size model -- and therefore
        # MAPQ and pair rescue near batch boundaries -- changes with --cores, and
        # the same input yields different calls. Reproducibility requires it.
        # -Y (soft-clip supplementary) rather than -M (mark them secondary): GATK4
        # handles proper supplementary records natively; -M is the legacy Picard
        # compatibility mode and hides them from duplicate marking.
        # samtools sort: -m caps per-thread memory (default 768M is multiplied by
        # {threads} on top of the resident bwa index) and -T keeps spill files off
        # the slow bind mount.
        ( bwa mem -Y -K 100000000 -t {threads} -R '{params.rg}' \
              {input.fasta} {input.r1} {input.r2} \
          | samtools sort -@ {threads} -m 1G -T {resources.tmpdir}/{wildcards.sample}.{wildcards.unit}.sort \
                -o {output} - ) 2> {log}
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

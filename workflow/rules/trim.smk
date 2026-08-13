# =============================================================================
# trim.smk — fastp adapter + quality trimming, per sequencing unit.
# Auto-detects PE adapters; emits a JSON report consumed by MultiQC.
# =============================================================================


rule fastp:
    input:
        unpack(get_unit_fastqs),
    output:
        # temp(): trimmed FASTQ are a pure intermediate (bwa-mem is the sole
        # consumer) — auto-delete once mapped, like results/mapped + results/dedup.
        # Regenerable from data/ raw FASTQ; keeps results/trimmed/ from ballooning.
        r1=temp("results/trimmed/{sample}.{unit}.R1.fastq.gz"),
        r2=temp("results/trimmed/{sample}.{unit}.R2.fastq.gz"),
        json="results/qc/fastp/{sample}.{unit}.fastp.json",
        html="results/qc/fastp/{sample}.{unit}.fastp.html",
    threads: 4
    log:
        "logs/fastp/{sample}.{unit}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe \
          --thread {threads} \
          --json {output.json} --html {output.html} \
          > {log} 2>&1
        """

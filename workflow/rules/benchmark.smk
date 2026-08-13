# =============================================================================
# benchmark.smk — score the pipeline against the GIAB HG001/NA12878 truth set
# with hap.py (GA4GH benchmarking). Compares results/filtered/all.filtered.vcf.gz
# to the GIAB high-confidence calls, restricted to (confident regions ∩ capture
# targets). Opt-in: `snakemake benchmark` (or `run benchmark`).
# =============================================================================
import os

TRUTH_VCF = config["benchmark"]["truth_vcf"]
TRUTH_BED = config["benchmark"]["truth_bed"]
_BENCH_BASE = config["benchmark"]["base_url"]

# Optional run label so multiple experiments (e.g. different capture kits) don't
# overwrite each other's hap.py output. Unset -> "results/benchmark/happy.*"
# (the Docker pre-bake target); set via an overlay (config/twist_onso.yaml).
_BENCH_LABEL = config["benchmark"].get("label", "")
_BENCH_PREFIX = "results/benchmark/happy" + (f".{_BENCH_LABEL}" if _BENCH_LABEL else "")

localrules:
    download_truth_vcf,
    download_truth_bed,
    benchmark,


rule download_truth_vcf:
    output:
        vcf=TRUTH_VCF,
        tbi=TRUTH_VCF + ".tbi",
    params:
        url=_BENCH_BASE + os.path.basename(TRUTH_VCF),
    log:
        "logs/benchmark/download_truth_vcf.log",
    conda:
        "../envs/download.yaml"
    shell:
        r"""
        ( wget --tries=3 -q -O {output.vcf} {params.url}
          wget --tries=3 -q -O {output.tbi} {params.url}.tbi ) 2> {log}
        """


rule download_truth_bed:
    output:
        TRUTH_BED,
    params:
        url=_BENCH_BASE + os.path.basename(TRUTH_BED),
    log:
        "logs/benchmark/download_truth_bed.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule benchmark_happy:
    input:
        query="results/filtered/all.filtered.vcf.gz",
        query_tbi="results/filtered/all.filtered.vcf.gz.tbi",
        truth=TRUTH_VCF,
        truth_tbi=TRUTH_VCF + ".tbi",
        truth_bed=TRUTH_BED,
        targets=config["intervals"]["bed"],
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
    output:
        summary=_BENCH_PREFIX + ".summary.csv",
        extended=_BENCH_PREFIX + ".extended.csv",
    params:
        prefix=_BENCH_PREFIX,
    threads: 4
    log:
        "logs/benchmark/happy.log",
    conda:
        "../envs/happy.yaml"
    shell:
        # -f = GIAB confident regions, -T = capture targets -> scored on the
        # intersection. Default xcmp engine (no RTG/vcfeval dependency).
        r"""
        hap.py {input.truth} {input.query} \
            -f {input.truth_bed} \
            -T {input.targets} \
            -r {input.fasta} \
            -o {params.prefix} \
            --threads {threads} > {log} 2>&1
        """


# Convenience target — run the benchmark.
rule benchmark:
    input:
        _BENCH_PREFIX + ".summary.csv",

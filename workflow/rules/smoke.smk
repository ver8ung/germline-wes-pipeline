# =============================================================================
# smoke.smk — tiny end-to-end smoke test. Simulates paired reads from a small
# region of the (downloaded) reference with wgsim, so the FULL DAG runs in
# minutes against real GRCh38 instead of a whole exome.
#
# Activated by running with the smoke config overlay, which repoints samples/
# units/intervals at .tests/smoke/ (see README, or `run.ps1 smoke` / `run.sh smoke`):
#
#     snakemake --configfile .tests/smoke/config.yaml --cores 4
#
# Under the DEFAULT config this rule is defined but never triggered (nothing
# requests results/smoke/*.fastq.gz), so it is harmless to always include.
# =============================================================================


rule smoke_simulate_reads:
    input:
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        bed=config["intervals"]["bed"],
    output:
        r1="results/smoke/smoke.R1.fastq.gz",
        r2="results/smoke/smoke.R2.fastq.gz",
    params:
        nreads=20000,   # ~40x over a 100 kb window at 2x100 bp — fast but callable
        readlen=100,
        seed=11,        # fixed seed -> reproducible smoke reads
    log:
        "logs/smoke/simulate.log",
    conda:
        "../envs/wgsim.yaml"
    shell:
        # Extract the first BED interval (1-based for samtools), simulate, bgzip.
        r"""
        ( reg=$(awk 'NR==1{{print $1":"$2+1"-"$3}}' {input.bed})
          echo "simulating $reg" >&2
          samtools faidx {input.fasta} "$reg" > results/smoke/region.fa
          wgsim -N {params.nreads} -1 {params.readlen} -2 {params.readlen} -S {params.seed} \
                results/smoke/region.fa \
                results/smoke/smoke.R1.fastq results/smoke/smoke.R2.fastq
          bgzip -f results/smoke/smoke.R1.fastq
          bgzip -f results/smoke/smoke.R2.fastq
          rm -f results/smoke/region.fa results/smoke/region.fa.fai ) 2> {log}
        """

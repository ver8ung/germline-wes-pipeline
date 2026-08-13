# =============================================================================
# qc.smk — FastQC on the raw reads (one report per read end, per unit).
# fastp (trim.smk) and the alignment/variant metrics also feed MultiQC.
# =============================================================================


def _fastqc_input(wildcards):
    """Raw FASTQ for one read end of a unit (R1 -> fq1, R2 -> fq2)."""
    col = "fq1" if wildcards.read == "R1" else "fq2"
    return units.loc[(wildcards.sample, wildcards.unit), col]


rule fastqc:
    input:
        _fastqc_input,
    output:
        html="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.html",
        zip="results/qc/fastqc/{sample}.{unit}.{read}_fastqc.zip",
    wildcard_constraints:
        read="R1|R2",
    log:
        "logs/fastqc/{sample}.{unit}.{read}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        # FastQC names outputs after the input basename, so run into a temp dir
        # and move the single html/zip to deterministic names.
        r"""
        tmp=$(mktemp -d)
        fastqc -q -o "$tmp" {input} > {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.html}
        mv "$tmp"/*_fastqc.zip {output.zip}
        rm -rf "$tmp"
        """

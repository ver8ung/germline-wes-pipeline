# =============================================================================
# common.smk — sample sheets, helper functions, wildcard constraints, targets.
# Included first by the Snakefile so every other rule file can use these.
# =============================================================================
import os
import re
import pandas as pd


# --- Load sample / unit sheets ----------------------------------------------
samples = (
    pd.read_table(config["samples"], dtype=str, comment="#")
    .set_index("sample", drop=False)
    .sort_index()
)

units = (
    pd.read_table(config["units"], dtype=str, comment="#")
    .set_index(["sample", "unit"], drop=False)
    .sort_index()
)

# Fail loudly on the most common sample-sheet mistakes.
assert not samples.empty, "config/samples.tsv has no rows"
assert not units.empty, "config/units.tsv has no rows"
_missing = set(units["sample"]) - set(samples["sample"])
assert not _missing, f"units.tsv references samples absent from samples.tsv: {_missing}"
assert units["fq2"].notna().all(), (
    "this pipeline expects paired-end reads — every unit needs fq1 AND fq2"
)

# Check the FASTQs exist up front. Without this the first sign of a bad path in
# units.tsv is an opaque MissingInputException naming only a fastp/fastqc output,
# possibly after the reference/index rules have already run for many minutes.
# A warning rather than an error, so `setup_reference` still works before the
# sequencing data has landed. Paths under results/ are skipped: those are
# produced by the workflow itself (the smoke test simulates its reads there).
_missing_fastqs = sorted(
    {
        str(p)
        for p in pd.concat([units["fq1"], units["fq2"]]).dropna()
        if not str(p).replace("\\", "/").startswith("results/")
        and not os.path.exists(str(p))
    }
)
if _missing_fastqs:
    logger.warning(
        "\nunits.tsv references FASTQ files that do not exist:\n  "
        + "\n  ".join(_missing_fastqs)
        + f"\n  Relative paths resolve against the repo root ({os.getcwd()}).\n"
        f"  Fix the fq1/fq2 columns in {config['units']} or mount the data\n"
        "  directory. Any target needing these reads will fail until then.\n"
    )


# --- Wildcard constraints ---------------------------------------------------
# Restrict wildcards to declared values so '{sample}.{unit}' filenames split
# unambiguously even when one name is a prefix of another (regex backtracks).
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples["sample"].unique()),
    unit="|".join(re.escape(u) for u in units["unit"].unique()),


# --- Shared constants -------------------------------------------------------
INTERVALS = "resources/intervals/targets.interval_list"   # padded exome targets
PADDING = int(config["intervals"]["padding"])             # bp flanking each target


# --- Reference helpers ------------------------------------------------------
def ref_dict():
    """GATK sequence dictionary path: <ref-prefix>.dict (drops the FASTA ext)."""
    fasta = config["ref"]["fasta"]
    for ext in (".fasta", ".fa", ".fna"):
        if fasta.endswith(ext):
            return fasta[: -len(ext)] + ".dict"
    return fasta + ".dict"


def known_sites_files():
    """Every known-sites VCF plus its .tbi — used as BQSR rule inputs."""
    out = []
    for ks in config["known_sites"].values():
        out += [ks["path"], ks["path"] + ".tbi"]
    return out


def known_sites_args(wildcards):
    """`--known-sites a.vcf.gz --known-sites b.vcf.gz ...` for BaseRecalibrator."""
    return " ".join(f"--known-sites {ks['path']}" for ks in config["known_sites"].values())


# --- Per-unit / per-sample helpers ------------------------------------------
def get_unit_fastqs(wildcards):
    """Raw R1/R2 for one sequencing unit (consumed by fastp and fastqc)."""
    u = units.loc[(wildcards.sample, wildcards.unit)]
    return {"r1": u.fq1, "r2": u.fq2}


def get_read_group(wildcards):
    """SAM @RG header for bwa -R. ID per-unit, SM per-sample (BQSR needs both)."""
    platform = units.loc[(wildcards.sample, wildcards.unit), "platform"]
    return (
        r"@RG\tID:{s}.{u}\tSM:{s}\tLB:{s}.{u}\tPU:{s}.{u}\tPL:{pl}"
    ).format(s=wildcards.sample, u=wildcards.unit, pl=platform)


def sample_units(sample):
    """All unit ids belonging to a sample (list, even when there is just one)."""
    return units.loc[[sample]]["unit"].tolist()


def get_sample_unit_bams(wildcards):
    """Per-unit sorted BAMs to merge into one per-sample BAM before MarkDuplicates."""
    return expand(
        "results/mapped/{sample}.{unit}.sorted.bam",
        sample=wildcards.sample,
        unit=sample_units(wildcards.sample),
    )


def all_sample_gvcfs():
    """Every per-sample GVCF — the cohort fan-in for CombineGVCFs."""
    return expand("results/called/{sample}.g.vcf.gz", sample=list(samples["sample"]))


def gvcf_combine_args(wildcards):
    """`--variant s1.g.vcf.gz --variant s2.g.vcf.gz ...` for CombineGVCFs."""
    return " ".join(f"--variant {g}" for g in all_sample_gvcfs())


# --- Annotation source switch -----------------------------------------------
def annotation_source():
    """Final annotation stage feeding the canonical output (ClinVar optional)."""
    if config["annotation"]["clinvar"]["activate"]:
        return "results/annotated/all.clinvar.vcf.gz"
    return "results/annotated/all.snpeff.vcf.gz"


# --- Top-level targets ------------------------------------------------------
def get_final_output():
    return [
        "results/annotated/all.annotated.vcf.gz",
        "results/annotated/all.annotated.vcf.gz.tbi",
        "results/qc/multiqc_report.html",
    ]

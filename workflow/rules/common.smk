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


# --- Guard: a cohort callset left over from a DIFFERENT set of samples --------
# results/genotyped/, results/filtered/ and results/annotated/ are cohort-level
# and share fixed paths across configs. Snakemake decides a target is up to date
# by comparing it with its DIRECT inputs, so when every cohort-level output
# already exists it never descends as far as combine_gvcfs and never notices that
# the sample set changed. Switching from one cohort to another therefore does NOT
# re-run anything: `all` reports success while results/annotated/ still describes
# the PREVIOUS cohort, and the benchmark then tries to pull this cohort's samples
# out of the previous cohort's VCF.
#
# The per-sample GVCFs are the exact test, and they are not temp(): a cohort VCF
# cannot possibly contain a sample whose GVCF was never produced.
_cohort_vcf = "results/genotyped/cohort.vcf.gz"
if os.path.exists(_cohort_vcf):
    _no_gvcf = sorted(
        s
        for s in samples["sample"]
        if not os.path.exists(f"results/called/{s}.g.vcf.gz")
    )
    if _no_gvcf:
        raise WorkflowError(
            f"\n{_cohort_vcf} exists but was built WITHOUT these samples:\n  "
            + "\n  ".join(_no_gvcf)
            + "\n\nThe cohort-level outputs on disk belong to a different cohort, and"
            "\nSnakemake will treat them as up to date rather than rebuilding them."
            "\nMove or delete the stale cohort results before running this config:"
            "\n  rm -rf results/genotyped results/filtered results/annotated"
            "\n(Per-sample results/called/*.g.vcf.gz and results/bqsr/*.recal.bam are"
            "\nkept, so samples already processed are not realigned.)\n"
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
    """SAM @RG header for bwa -R. ID per-unit, SM per-sample, LB per LIBRARY.

    LB used to be fabricated as `{sample}.{unit}`, i.e. one library per unit. That
    is only right when every unit really is its own library, and it is wrong for
    the default sheet: NIST7035 and NIST7086 are two libraries sequenced over two
    lanes each. MarkDuplicates keys duplicate detection on (LB, 5' position), so
    calling each lane its own library meant cross-lane PCR duplicates of the same
    library were NEVER marked -- under-counting duplicates and overstating
    effective depth.

    The `library` column is optional: without it the old fabricated value is used,
    so sheets that never had one behave exactly as before (correctly so for
    units.twist_onso.tsv, whose 8 units genuinely are 8 separate libraries).
    """
    u = units.loc[(wildcards.sample, wildcards.unit)]
    if "library" in units.columns and pd.notna(u["library"]):
        library = u["library"]
    else:
        library = f"{wildcards.sample}.{wildcards.unit}"
    return (
        r"@RG\tID:{s}.{u}\tSM:{s}\tLB:{lb}\tPU:{s}.{u}\tPL:{pl}"
    ).format(s=wildcards.sample, u=wildcards.unit, lb=library, pl=u["platform"])


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

# =============================================================================
# benchmark.smk — score the pipeline against the GIAB truth sets with hap.py
# (GA4GH benchmarking), PER SAMPLE.
#
# For each sample that has a `benchmark.truth` entry, one sample's calls are
# extracted from the cohort callset and compared to that sample's GIAB
# high-confidence calls, restricted to (confident regions ∩ capture targets).
# Opt-in: `snakemake benchmark` (or `run benchmark`).
#
# Two things here are easy to get wrong and expensive to discover late:
#   * hap.py has NO sample-selection flag. A multi-sample query is not a
#     supported input, so the per-sample extraction below is mandatory -- see
#     rule benchmark_query_vcf.
#   * vcfeval is a JVM tool that sizes its own heap from the HOST's RAM even
#     inside a container -- see RTG_MEM in the two rules that invoke rtg.
# =============================================================================
import os
import re

_TRUTH = config["benchmark"]["truth"]

# Optional run label so multiple experiments (e.g. different capture kits) don't
# overwrite each other's hap.py output. Unset -> "results/benchmark/happy.<sample>.*"
# (the Docker pre-bake target); set via an overlay (config/twist_onso.yaml).
_BENCH_LABEL = config["benchmark"].get("label", "")
_LBL = f".{_BENCH_LABEL}" if _BENCH_LABEL else ""

# --- Comparison engine -------------------------------------------------------
# Validate here rather than letting a typo reach hap.py: a misspelling would
# otherwise surface as an argparse error deep inside a rule that only runs after
# the whole pipeline has completed.
_BENCH_ENGINE = config["benchmark"].get("engine", "vcfeval")
if _BENCH_ENGINE not in ("xcmp", "vcfeval"):
    raise WorkflowError(
        f"benchmark.engine must be 'xcmp' or 'vcfeval', got {_BENCH_ENGINE!r}"
    )

# RTG's Sequence Data File — the indexed reference vcfeval requires. Kept beside
# the other derived reference data, but deliberately NOT part of `setup_reference`
# (see rule rtg_format_sdf).
RTG_SDF = "resources/rtg/" + os.path.basename(config["ref"]["fasta"]) + ".sdf"

_ENGINE_ARGS = (
    f"--engine vcfeval --engine-vcfeval-template {RTG_SDF}"
    if _BENCH_ENGINE == "vcfeval"
    else "--engine xcmp"
)

# --- Stratifications ---------------------------------------------------------
_STRAT = config["benchmark"].get("stratification", {})
_STRAT_ON = bool(_STRAT.get("activate", False))
_STRAT_DIR = "resources/stratifications/" + str(_STRAT.get("version", "v3.6"))
_STRAT_TSV = f"{_STRAT_DIR}/manifest.tsv"
_STRAT_REGIONS = _STRAT.get("regions", {}) if _STRAT_ON else {}
_STRAT_BEDS = [f"{_STRAT_DIR}/{name}.bed.gz" for name in _STRAT_REGIONS]


def bench_prefix(sample):
    """hap.py -o prefix for one sample of the current (optionally labelled) run."""
    return f"results/benchmark/happy{_LBL}.{sample}"


def benchmarkable_samples():
    """Samples in THIS cohort that have a truth set configured.

    The smoke sample has none, so `benchmark` under the smoke overlay resolves to
    an empty input list rather than to a MissingInputException.
    """
    return sorted(set(samples["sample"]) & set(_TRUTH))


localrules:
    download_truth_vcf,
    download_truth_bed,
    download_stratification_bed,
    stratification_manifest,
    benchmark,


# --- Truth data --------------------------------------------------------------
# Local paths are FIXED (<sample>.truth.*) and the remote filename comes from
# config. Deriving the URL from the local basename -- as this file used to --
# silently requires remote and local names to match, which GIAB does not honour:
# HG001 publishes ..._benchmark.bed, the Ashkenazim trio ..._benchmark_noinconsistent.bed.
rule download_truth_vcf:
    output:
        vcf="resources/benchmark/{sample}.truth.vcf.gz",
        tbi="resources/benchmark/{sample}.truth.vcf.gz.tbi",
    params:
        url=lambda wc: _TRUTH[wc.sample]["base_url"] + _TRUTH[wc.sample]["vcf"],
    wildcard_constraints:
        # Narrower than the cohort-wide constraint in common.smk: an unconfigured
        # sample must be "no rule to produce this" rather than a KeyError raised
        # from inside the params lambda while the DAG is being built.
        sample="|".join(re.escape(s) for s in _TRUTH),
    log:
        "logs/benchmark/download_truth_vcf.{sample}.log",
    conda:
        "../envs/download.yaml"
    shell:
        r"""
        ( wget --tries=3 -q -O {output.vcf} {params.url}
          wget --tries=3 -q -O {output.tbi} {params.url}.tbi ) 2> {log}
        """


rule download_truth_bed:
    output:
        "resources/benchmark/{sample}.truth.bed",
    params:
        url=lambda wc: _TRUTH[wc.sample]["base_url"] + _TRUTH[wc.sample]["bed"],
    wildcard_constraints:
        sample="|".join(re.escape(s) for s in _TRUTH),
    log:
        "logs/benchmark/download_truth_bed.{sample}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


# --- Stratification regions --------------------------------------------------
rule download_stratification_bed:
    output:
        f"{_STRAT_DIR}/{{name}}.bed.gz",
    params:
        # NOT base_url + basename(output): the local name is the short manifest
        # key while the remote path carries a Category/ prefix and a GRCh38_ file
        # prefix, so the two genuinely differ.
        url=lambda wc: _STRAT["base_url"] + _STRAT_REGIONS[wc.name],
    wildcard_constraints:
        name="|".join(re.escape(n) for n in _STRAT_REGIONS) or "^$",
    log:
        "logs/benchmark/strat_{name}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule stratification_manifest:
    """The TSV hap.py's --stratification consumes: `name <TAB> path`.

    hap.py resolves each path RELATIVE TO THE TSV'S OWN DIRECTORY, so writing
    bare filenames next to the BEDs keeps the manifest independent of where the
    repo is checked out or mounted.
    """
    input:
        _STRAT_BEDS,
    output:
        _STRAT_TSV,
    run:
        with open(output[0], "w") as fh:
            for name in _STRAT_REGIONS:
                fh.write(f"{name}\t{name}.bed.gz\n")


# --- vcfeval reference index -------------------------------------------------
rule rtg_format_sdf:
    """RTG Sequence Data File — the reference index vcfeval needs.

    Deliberately NOT part of `rule setup_reference`: that target is the run-once
    bootstrap every user pays, and this is a ~2.5 GB / ~20 min index only people
    running the benchmark need.

    The output is a DIRECTORY, and Snakemake's completeness check for a directory
    output is only "does the path exist". An SDF is ~10 binary files with no
    single member whose presence proves the set is complete, so a container
    OOM-kill or Ctrl-C partway through `rtg format` would leave a half-written SDF
    that the next run treats as finished -- and vcfeval then fails deep inside
    hap.py complaining about the template rather than about the SDF. Building into
    a sibling .tmp and renaming means the declared path appears only after rtg
    exits 0.
    """
    input:
        fasta=config["ref"]["fasta"],
    output:
        sdf=directory(RTG_SDF),
    params:
        rtg_mem=lambda wildcards, resources: f"{int(resources.mem_mb * 0.66)}m",
    resources:
        mem_mb=config["resources"]["benchmark"]["mem_mb"],
    log:
        "logs/benchmark/rtg_format.log",
    conda:
        "../envs/happy.yaml"
    shell:
        r"""
        # RTG sizes its JVM heap at 90% of DETECTED PHYSICAL RAM. Inside a Docker
        # container /proc/meminfo still reports the HOST's memory, so that default
        # is 90% of the host rather than of --memory: the JVM overshoots the cgroup
        # and the container is OOM-killed with no Java stack trace at all. RTG_MEM
        # is the only override reachable here -- the documented `rtg RTG_MEM=8g
        # <cmd>` first-argument form is unavailable because hap.py invokes rtg
        # itself, as plain `rtg vcfeval ...`.
        export RTG_MEM={params.rtg_mem}
        # Snakemake does not pre-remove directory outputs, and `rtg format` aborts
        # if its target exists, so a retry needs the rm.
        ( rm -rf {output.sdf} {output.sdf}.tmp
          rtg format -f fasta -o {output.sdf}.tmp {input.fasta}
          mv {output.sdf}.tmp {output.sdf} ) > {log} 2>&1
        """


# --- Per-sample query --------------------------------------------------------
rule benchmark_query_vcf:
    """One sample's calls, extracted from the cohort callset for hap.py.

    hap.py cannot score a multi-sample VCF and exposes no flag to pick a column
    (verified against 0.3.15's own --help: there is no such option). Under
    --engine vcfeval, rtg aborts outright. Under the old xcmp engine it did NOT
    abort -- it silently scored the FIRST sample column and reported the result
    under whichever truth set was named, which at n=1 is right by accident and at
    n>1 is three plausible wrong numbers.

    This therefore runs unconditionally, including at n=1: a rule that only exists
    for n>1 is a rule nobody ever tests.
    """
    input:
        vcf="results/filtered/all.filtered.vcf.gz",
        tbi="results/filtered/all.filtered.vcf.gz.tbi",
    output:
        # Lives under results/benchmark/, never results/filtered/, so it cannot
        # collide with all.filtered.vcf.gz or cohort.{vartype}.* even if someone
        # names a sample "all" or "cohort".
        vcf="results/benchmark/{sample}.query.vcf.gz",
        tbi="results/benchmark/{sample}.query.vcf.gz.tbi",
    log:
        "logs/benchmark/query_{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        # -a trims ALT alleles this sample does not carry and recomputes AC/AN;
        # GT="alt" then drops sites where the sample is hom-ref or no-call, which a
        # cohort VCF is full of and which hap.py would otherwise weigh as query
        # records.
        # FILTER is left exactly as the cohort computed it. The hard filters ARE
        # cohort-level, and hap.py's ALL-vs-PASS split has to reflect the filters
        # the pipeline actually applied, not a per-sample recomputation.
        ( bcftools view -s {wildcards.sample} -a -Ou {input.vcf} \
          | bcftools view -i 'GT="alt"' -Oz -o {output.vcf}
          tabix -p vcf {output.vcf} ) 2> {log}
        """


# --- Scoring -----------------------------------------------------------------
rule benchmark_happy:
    input:
        query="results/benchmark/{sample}.query.vcf.gz",
        query_tbi="results/benchmark/{sample}.query.vcf.gz.tbi",
        truth="resources/benchmark/{sample}.truth.vcf.gz",
        truth_tbi="resources/benchmark/{sample}.truth.vcf.gz.tbi",
        truth_bed="resources/benchmark/{sample}.truth.bed",
        # Decoupled from intervals.bed so the same callset can be re-scored
        # against a narrower region without re-running the pipeline.
        targets=config["benchmark"].get("targets_bed", config["intervals"]["bed"]),
        fasta=config["ref"]["fasta"],
        fai=config["ref"]["fasta"] + ".fai",
        sdf=[RTG_SDF] if _BENCH_ENGINE == "vcfeval" else [],
        strat_tsv=[_STRAT_TSV] if _STRAT_ON else [],
        strat_beds=_STRAT_BEDS,
    output:
        summary=f"results/benchmark/happy{_LBL}.{{sample}}.summary.csv",
        extended=f"results/benchmark/happy{_LBL}.{{sample}}.extended.csv",
        # Previously undeclared, though dev/collect_benchmark_evidence.sh reads the
        # runinfo and dev/diag_capture.sh reads the VCF -- the DAG was lying about
        # two files the repo's own tooling depends on.
        # .metrics.json.gz and the five .roc.*.csv.gz remain undeclared on purpose:
        # they are genuine hap.py side-products and nothing here consumes them.
        runinfo=f"results/benchmark/happy{_LBL}.{{sample}}.runinfo.json",
        vcf=f"results/benchmark/happy{_LBL}.{{sample}}.vcf.gz",
        vcf_tbi=f"results/benchmark/happy{_LBL}.{{sample}}.vcf.gz.tbi",
    params:
        prefix=f"results/benchmark/happy{_LBL}.{{sample}}",
        # In params, not just inlined into the shell string: changing an engine or
        # stratification setting alters no input file's mtime, so without this
        # Snakemake would leave a stale result in place and call the DAG complete.
        # Snakemake's default rerun-triggers include `params`, so this makes the
        # config knob actually re-score.
        engine_args=_ENGINE_ARGS,
        strat_args=(f"--stratification {_STRAT_TSV}" if _STRAT_ON else ""),
        rtg_mem=lambda wildcards, resources: f"{int(resources.mem_mb * 0.66)}m",
    threads: 4
    resources:
        mem_mb=config["resources"]["benchmark"]["mem_mb"],
    log:
        # Was a hardcoded logs/benchmark/happy.log for EVERY run, so labelled and
        # unlabelled runs -- and now every sample -- overwrote each other's log,
        # and the log never described the run you were looking at.
        f"logs/benchmark/happy{_LBL}.{{sample}}.log",
    conda:
        "../envs/happy.yaml"
    shell:
        # -f = that sample's GIAB confident regions, -T = capture targets, so
        # scoring happens on the intersection.
        r"""
        export RTG_MEM={params.rtg_mem}
        hap.py {input.truth} {input.query} \
            -f {input.truth_bed} \
            -T {input.targets} \
            -r {input.fasta} \
            -o {params.prefix} \
            {params.engine_args} \
            {params.strat_args} \
            --scratch-prefix {resources.tmpdir} \
            --threads {threads} > {log} 2>&1
        """


# --- Mendelian consistency ---------------------------------------------------
# An accuracy signal that needs NO truth set: it asks whether the child's
# genotypes could have been inherited from the parents'. That makes it completely
# independent of GIAB -- it covers the whole callset rather than just the
# high-confidence regions, and it cannot be flattered by anything specific to the
# benchmark regions. Only meaningful when the cohort actually contains a trio.
_PED_COLS = ("sex", "paternal_id", "maternal_id")


def has_pedigree():
    """True when samples.tsv carries pedigree columns AND some sample names a parent."""
    if not all(c in samples.columns for c in _PED_COLS):
        return False
    known = set(samples["sample"])
    return any(
        str(samples.loc[s, "paternal_id"]) in known
        and str(samples.loc[s, "maternal_id"]) in known
        for s in samples["sample"]
    )


rule write_pedigree:
    """PLINK PED derived from samples.tsv.

    Kept derived rather than shipped as a standalone .ped so there is one source
    of truth: a separate file would be a second sample list to keep in sync with
    the first, and the failure mode of them disagreeing is silent.

    bcftools reads columns <ignored>,proband,father,mother,sex; the sixth
    (phenotype) column is written for standard PLINK compatibility and ignored.
    """
    output:
        "results/pedigree/cohort.ped",
    log:
        "logs/benchmark/pedigree.log",
    run:
        with open(output[0], "w") as fh:
            for s in samples["sample"]:
                r = samples.loc[s]
                fh.write(f"FAM\t{s}\t{r.paternal_id}\t{r.maternal_id}\t{r.sex}\t0\n")


rule mendelian_check:
    """Trio consistency on the PASS callset -- what the pipeline actually ships."""
    input:
        vcf="results/filtered/all.pass.vcf.gz",
        tbi="results/filtered/all.pass.vcf.gz.tbi",
        ped="results/pedigree/cohort.ped",
    output:
        f"results/benchmark/mendelian{_LBL}.txt",
    params:
        # Autosomes only. The sex chromosomes follow different inheritance rules
        # (bcftools has --rules for that) and the GIAB truth sets this run is
        # scored against are chr1-22 anyway, so including them would mix a
        # different question into the same number.
        regions=",".join(f"chr{i}" for i in range(1, 23)),
    log:
        f"logs/benchmark/mendelian{_LBL}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        bcftools +mendelian2 -P {input.ped} -m c \
            -r {params.regions} {input.vcf} > {output} 2> {log}
        """


# Convenience target — score every sample in this cohort that has a truth set,
# plus the Mendelian check when the cohort is a trio.
rule benchmark:
    input:
        expand(
            f"results/benchmark/happy{_LBL}.{{sample}}.summary.csv",
            sample=benchmarkable_samples(),
        ),
        [f"results/benchmark/mendelian{_LBL}.txt"] if has_pedigree() else [],

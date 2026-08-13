# =============================================================================
# Self-contained NGS orchestrator image — tools PRE-BAKED.
#
# Every per-rule Conda environment (workflow/envs/*.yaml) is created at BUILD
# time into /opt/snakemake-envs, which lives OUTSIDE the runtime bind-mount.
# The container therefore runs fully offline with no slow first-run env setup.
#
# Snakemake identifies a conda env by the hash of its YAML *content* (+ the
# conda prefix), not its path — so the bind-mounted workflow/envs/*.yaml at run
# time map onto these baked envs as long as the YAML content is unchanged and
# the SAME --conda-prefix is passed (the run.ps1 / run.sh launchers do this).
#
# BUILD CONTEXT MUST BE THE PROJECT ROOT (the COPYs below read workflow/ and
# config/). Build with:  docker build -t ngs-germline-wes -f docker/Dockerfile .
#
# BASE is pinnable by digest for reproducible rebuilds — the launchers pass the
# value from docker/base-image.txt as --build-arg BASE=...; default stays on the
# floating tag so a bare `docker build` still works.
# =============================================================================
ARG BASE=condaforge/miniforge3:latest
FROM ${BASE}

# Stamped by the launchers (git sha / content hash) for traceability:
#   docker inspect --format '{{ index .Config.Labels "org.opencontainers.image.version" }}' <image>
ARG IMAGE_TAG=unknown
LABEL org.opencontainers.image.title="ngs-germline-wes" \
      org.opencontainers.image.description="GATK germline WES pipeline (FASTQ -> annotated VCF), per-rule tools pre-baked" \
      org.opencontainers.image.version="${IMAGE_TAG}"

# Must match the --conda-prefix passed by the launchers at run time.
ENV SMK_CONDA_PREFIX=/opt/snakemake-envs

# 1) Snakemake in its OWN env with a pinned Python, so it can never conflict
#    with whatever Python the base image ships. Snakemake 8 supports modern
#    Python and is the version this workflow was validated against. Putting the
#    env on PATH makes `snakemake`/`python` resolve here; conda+mamba remain on
#    PATH from /opt/conda/bin for env creation.
RUN mamba create -y -n smk -c conda-forge -c bioconda \
        python=3.12 \
        "snakemake-minimal>=8.20,<9" \
        pandas=2.2.* \
        tabulate \
    && mamba clean -afy
ENV PATH=/opt/conda/envs/smk/bin:$PATH

# 2) Pre-create every per-rule conda env.
#    A tiny build-only config is generated in-image (no dependency on any host
#    test assets being present in the build context) purely so the DAG resolves
#    for `--conda-create-envs-only`. It pulls in the exact same envs the real
#    pipeline uses: gatk, bwa, samtools, fastp, fastqc, snpeff, multiqc,
#    download (wget/htslib/bcftools), and wgsim.
WORKDIR /build
COPY workflow/ workflow/
COPY config/ config/
RUN mkdir -p _buildcfg \
    && printf 'chr20\t30000000\t30100000\n' > _buildcfg/smoke.bed \
    && printf 'sample\tcondition\nsmoke\tnormal\n' > _buildcfg/samples.tsv \
    && printf 'sample\tunit\tplatform\tfq1\tfq2\nsmoke\tL001\tILLUMINA\tresults/smoke/smoke.R1.fastq.gz\tresults/smoke/smoke.R2.fastq.gz\n' > _buildcfg/units.tsv \
    && printf 'samples: _buildcfg/samples.tsv\nunits: _buildcfg/units.tsv\nintervals:\n  bed: _buildcfg/smoke.bed\n  padding: 100\n' > _buildcfg/config.yaml \
    && snakemake \
        --snakefile workflow/Snakefile \
        --configfile _buildcfg/config.yaml \
        --use-conda \
        --conda-prefix "$SMK_CONDA_PREFIX" \
        --conda-create-envs-only --cores 1 \
        all \
        resources/intervals/twist_exome_v2.GRCh38.bed \
        resources/intervals/illumina_exome_v2.5.GRCh38.bed \
        resources/intervals/idt_xgen_exome_v2.GRCh38.bed \
        resources/intervals/nextera_expandedexome.GRCh38.bed \
        results/benchmark/happy.summary.csv \
    && mamba clean -afy \
    && conda clean -afy

# Record a single digest over the baked env YAMLs; the launchers recompute this
# from the bind-mounted files at run time and warn on drift (see check-bake).
RUN cat $(ls workflow/envs/*.yaml | sort) | sha256sum | cut -d' ' -f1 \
        > "$SMK_CONDA_PREFIX/.bake-digest"

# Staleness guard (option 5): `check-bake` re-hashes the bind-mounted env files
# (project mounted at /workflow) and compares to the baked digest. Exit 0 = in
# sync, non-zero = the image is stale and should be rebuilt. The launchers run
# this and warn; keeping it as a baked script avoids shell-quoting pitfalls.
RUN printf '%s\n' \
    '#!/usr/bin/env bash' \
    'cur=$(cat $(ls /workflow/workflow/envs/*.yaml | sort) | sha256sum | cut -d" " -f1)' \
    'baked=$(cat /opt/snakemake-envs/.bake-digest 2>/dev/null)' \
    '[ "$cur" = "$baked" ]' \
    > /usr/local/bin/check-bake \
    && chmod +x /usr/local/bin/check-bake

# 3) Runtime workdir — the launchers bind-mount the project here.
WORKDIR /workflow
CMD ["bash"]

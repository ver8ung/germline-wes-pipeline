# =============================================================================
# annotation.smk — functional annotation with SnpEff, an optional ClinVar
# overlay with SnpSift, then a finalize step that produces the canonical
# results/annotated/all.annotated.vcf.gz regardless of which path ran.
# =============================================================================

SNPEFF_DB = config["annotation"]["snpeff"]["db"]
SNPEFF_DATADIR = config["annotation"]["snpeff"]["datadir"]

# Contig-naming guard. The Broad GRCh38 reference is chr-named, and so is SnpEff's
# `hg38` (UCSC/RefSeq) database. An Ensembl-style database (`GRCh38.*`) uses bare
# contig names, and SnpEff does not error on the mismatch: it simply matches
# nothing and emits a VCF where every variant is un-annotated. Warn loudly rather
# than fail, since renaming contigs to suit an Ensembl DB is a legitimate setup.
if SNPEFF_DB.upper().startswith("GRCH"):
    logger.warning(
        f"\nSnpEff database {SNPEFF_DB!r} is Ensembl-style (no 'chr' prefix), but the\n"
        "  Broad GRCh38 reference used here is chr-named. SnpEff will NOT report an\n"
        "  error on this mismatch — it will silently annotate nothing. Use the\n"
        "  chr-named 'hg38' database, or rename contigs before annotation.\n"
    )

localrules:
    snpeff_download,
    finalize_annotation,


rule snpeff_download:
    output:
        directory(f"{SNPEFF_DATADIR}/{SNPEFF_DB}"),
    log:
        "logs/annotate/snpeff_download.log",
    conda:
        "../envs/snpeff.yaml"
    shell:
        # Absolute, QUOTED -dataDir so SnpEff neither falls back to its install
        # location nor word-splits a path containing spaces.
        'snpEff download -dataDir "$(pwd)/{SNPEFF_DATADIR}" {SNPEFF_DB} 2> {log}'


rule snpeff:
    input:
        # PASS-only, not the flagged VCF -- see rule select_pass in filtering.smk.
        vcf="results/filtered/all.pass.vcf.gz",
        tbi="results/filtered/all.pass.vcf.gz.tbi",
        db=f"{SNPEFF_DATADIR}/{SNPEFF_DB}",
    output:
        vcf="results/annotated/all.snpeff.vcf.gz",
        tbi="results/annotated/all.snpeff.vcf.gz.tbi",
        csv="results/qc/snpeff/all.snpeff.csv",          # MultiQC-readable stats
        html="results/qc/snpeff/all.snpeff.summary.html",
    params:
        xmx=lambda wildcards, resources: int(resources.mem_mb * 0.85),
    resources:
        mem_mb=config["resources"]["annotate"]["mem_mb"],
    log:
        "logs/annotate/snpeff.log",
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        ( snpEff -Xmx{params.xmx}m ann \
              -dataDir "$(pwd)/{SNPEFF_DATADIR}" \
              -nodownload \
              -csvStats {output.csv} \
              -stats {output.html} \
              {SNPEFF_DB} {input.vcf} \
          | bgzip -c > {output.vcf}
          tabix -p vcf {output.vcf} ) 2> {log}
        """


rule snpsift_clinvar:
    input:
        vcf="results/annotated/all.snpeff.vcf.gz",
        tbi="results/annotated/all.snpeff.vcf.gz.tbi",
        clinvar=config["annotation"]["clinvar"]["path"],
        clinvar_tbi=config["annotation"]["clinvar"]["path"] + ".tbi",
    output:
        vcf="results/annotated/all.clinvar.vcf.gz",
        tbi="results/annotated/all.clinvar.vcf.gz.tbi",
    log:
        "logs/annotate/snpsift_clinvar.log",
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        ( SnpSift annotate -tabix \
              {input.clinvar} {input.vcf} \
          | bgzip -c > {output.vcf}
          tabix -p vcf {output.vcf} ) 2> {log}
        """


# Canonical final output. annotation_source() resolves to the SnpSift output
# when ClinVar is active, otherwise the SnpEff output.
rule finalize_annotation:
    input:
        vcf=annotation_source(),
        tbi=annotation_source() + ".tbi",
    output:
        vcf="results/annotated/all.annotated.vcf.gz",
        tbi="results/annotated/all.annotated.vcf.gz.tbi",
    shell:
        "cp {input.vcf} {output.vcf} && cp {input.tbi} {output.tbi}"

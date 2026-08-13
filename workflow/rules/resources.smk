# =============================================================================
# resources.smk — fetch + index the GRCh38 reference, known-sites VCFs, and
# build the padded exome interval list. Run once via `--until setup_reference`,
# or let these rules trigger on demand. All are localrules (no cluster submit).
# =============================================================================

# output-path -> source-URL lookup for the known-sites downloader
KNOWN_BY_PATH = {ks["path"]: ks["url"] for ks in config["known_sites"].values()}

localrules:
    download_reference,
    download_known_site,
    bed_to_intervals,
    make_contig_map,
    download_clinvar,


# --- Reference FASTA + its indexes ------------------------------------------
rule download_reference:
    output:
        config["ref"]["fasta"],
    params:
        url=config["ref"]["fasta_url"],
    log:
        "logs/resources/download_reference.log",
    conda:
        "../envs/download.yaml"
    shell:
        # No --continue: combining resume with -O can yield a corrupt file.
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule samtools_faidx:
    input:
        config["ref"]["fasta"],
    output:
        config["ref"]["fasta"] + ".fai",
    log:
        "logs/resources/faidx.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} 2> {log}"


rule sequence_dict:
    input:
        config["ref"]["fasta"],
    output:
        ref_dict(),
    log:
        "logs/resources/create_dict.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} 2> {log}"


rule bwa_index:
    input:
        config["ref"]["fasta"],
    output:
        multiext(config["ref"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/resources/bwa_index.log",
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input} 2> {log}"


# --- Known-sites VCFs (for BQSR) --------------------------------------------
# Generic downloader: matches any configured known-sites path under resources/known/.
rule download_known_site:
    output:
        "resources/known/{name}.vcf.gz",
    params:
        url=lambda wc: KNOWN_BY_PATH["resources/known/" + wc.name + ".vcf.gz"],
    log:
        "logs/resources/download_{name}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule tabix_known_site:
    input:
        "resources/known/{name}.vcf.gz",
    output:
        "resources/known/{name}.vcf.gz.tbi",
    log:
        "logs/resources/tabix_{name}.log",
    conda:
        "../envs/download.yaml"
    shell:
        "tabix -p vcf {input} 2> {log}"


# --- Exome capture intervals ------------------------------------------------
# BED -> GATK interval_list. Padding is applied where intervals are *consumed*
# (HaplotypeCaller / BaseRecalibrator use --interval-padding), so this stays raw.
rule bed_to_intervals:
    input:
        bed=config["intervals"]["bed"],
        dict=ref_dict(),
    output:
        INTERVALS,
    log:
        "logs/resources/bed_to_intervals.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        "gatk BedToIntervalList -I {input.bed} -O {output} -SD {input.dict} 2> {log}"


# --- Optional ClinVar overlay (only built when annotation.clinvar.activate) --
# NCBI ClinVar is Ensembl-named (1, 2, ... MT). Our VCF is chr-named, so we
# generate an ens->chr map and rewrite ClinVar's contigs once, here.
rule make_contig_map:
    output:
        "resources/annotation/ens2chr.txt",
    log:
        "logs/resources/contig_map.log",
    shell:
        r"""
        ( for c in $(seq 1 22) X Y; do echo -e "${{c}}\tchr${{c}}"; done
          echo -e "MT\tchrM" ) > {output} 2> {log}
        """


rule download_clinvar:
    input:
        chrmap="resources/annotation/ens2chr.txt",
    output:
        vcf=config["annotation"]["clinvar"]["path"],
        tbi=config["annotation"]["clinvar"]["path"] + ".tbi",
    params:
        url=config["annotation"]["clinvar"]["url"],
    log:
        "logs/resources/clinvar.log",
    conda:
        "../envs/download.yaml"
    shell:
        r"""
        ( wget --tries=3 -q -O resources/annotation/clinvar.raw.vcf.gz {params.url}
          bcftools annotate --rename-chrs {input.chrmap} \
              -Oz -o {output.vcf} resources/annotation/clinvar.raw.vcf.gz
          tabix -p vcf {output.vcf}
          rm -f resources/annotation/clinvar.raw.vcf.gz ) 2> {log}
        """

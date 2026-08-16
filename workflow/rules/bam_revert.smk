# =============================================================================
# bam_revert.smk — turn an aligned BAM back into paired FASTQ so it can enter the
# pipeline at the top.
#
# Needed for the GIAB Ashkenazim trio: the only Illumina exome data GIAB
# publishes for HG002/HG003/HG004 is a position-sorted, duplicate-marked BAM
# (OsloUniversityHospital_Exome, 2015, Agilent SureSelect V5), not FASTQ. The
# source BAM's own alignment is discarded entirely -- it was made against an older
# reference with older tools, and reverting to reads means this pipeline's
# alignment, duplicate marking and recalibration all apply as they would to any
# other input.
#
# Same shape as smoke.smk: it PRODUCES reads into results/, and the units sheet
# points at them. common.smk deliberately skips the FASTQ-existence check for
# paths under results/ for exactly this reason.
#
# Opt-in: only samples listed under `revert_bams` in the config get this rule, so
# it is inert for every config that supplies real FASTQ.
# =============================================================================
import re

_REVERT_BAMS = config.get("revert_bams", {}) or {}


rule revert_bam_to_fastq:
    input:
        bam=lambda wc: _REVERT_BAMS[wc.sample],
    output:
        r1="results/reverted/{sample}.R1.fastq.gz",
        r2="results/reverted/{sample}.R2.fastq.gz",
        # Declared as real outputs rather than sent to /dev/null. A large
        # singleton count is the signature of a revert that lost pairing, and
        # discarding it silently is how that goes unnoticed until the insert-size
        # distribution looks strange three stages later.
        singleton="results/reverted/{sample}.singleton.fastq.gz",
        other="results/reverted/{sample}.other.fastq.gz",
    wildcard_constraints:
        # Empty config -> a constraint that matches nothing, so the rule simply
        # never applies instead of raising a KeyError from the input lambda.
        sample="|".join(re.escape(s) for s in _REVERT_BAMS) or "^$",
    threads: config["resources"]["sort"]["threads"]
    resources:
        mem_mb=config["resources"]["sort"]["mem_mb"],
    log:
        "logs/reverted/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
        ( set -e
          # Fail on a truncated download NOW rather than ~40 minutes into collate.
          samtools quickcheck -v {input.bam}

          # `collate`, NOT `sort -n`: it groups mates into adjacency without a full
          # sort. This is load-bearing -- the input is POSITION-sorted, so mates sit
          # megabases apart, and `samtools fastq` pairs by ADJACENCY rather than by
          # name. Feeding it a position-sorted BAM directly would emit essentially
          # every read as a singleton and yield zero usable pairs, while exiting 0.
          samtools collate -u -O -@ {threads} \
              -T {resources.tmpdir}/collate.{wildcards.sample} \
              {input.bam} \
          | samtools fastq -@ {threads} \
                -F 0x900 \
                -n \
                -1 {output.r1} -2 {output.r2} \
                -0 {output.other} -s {output.singleton} \
                -
          # -F 0x900 (secondary 0x100 + supplementary 0x800) is samtools' own
          #   default, restated because it matters: supplementary records are
          #   hard-clipped, so reverting them would emit truncated reads carrying
          #   duplicate names and break pairing. 0x400 is deliberately NOT added --
          #   duplicates must survive the revert so MarkDuplicates sees the library's
          #   true complexity rather than a pre-filtered view of it.
          # -s is mandatory, not tidiness: without it samtools writes singletons
          #   INTO the -1/-2 streams, which silently destroys pairing.
          # -n drops the /1 /2 name suffixes; bwa strips them anyway.

          # Counts into the log so a broken revert is visible immediately.
          echo "R1 pairs:       $(( $(zcat {output.r1} | wc -l) / 4 ))"
          echo "R2 pairs:       $(( $(zcat {output.r2} | wc -l) / 4 ))"
          echo "singletons:     $(( $(zcat {output.singleton} | wc -l) / 4 ))"
          echo "unpaired/other: $(( $(zcat {output.other} | wc -l) / 4 ))"
        ) > {log} 2>&1
        """

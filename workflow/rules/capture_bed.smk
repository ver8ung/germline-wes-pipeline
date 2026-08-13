# =============================================================================
# capture_bed.smk — build a GRCh38 capture-target BED. Four providers:
#   (1) Twist Exome 2.0    — native GRCh38; download + trim + sort + merge.
#   (2) Illumina Exome 2.5 — Twist-built, native GRCh38; same recipe as (1).
#   (3) IDT xGen Exome v2  — native GRCh38 (independent design); same recipe.
#   (4) Nextera Expanded   — Illumina ships hg19 only, so we download that BED +
#       the UCSC hg19->hg38 chain, lift it over, then primary-filter+sort+merge.
#
# Pick one via config["intervals"]["bed"]; `get_capture_bed` builds that path.
# Opt-in: `snakemake get_capture_bed` (or `run get_capture_bed`). For your OWN
# kit, just place its GRCh38 BED at the configured path instead.
# =============================================================================

# Fixed output paths (NOT config-driven): the producers must keep the same
# output regardless of which config/overlay is active, so the Docker build can
# pre-bake their envs by requesting these paths directly. config/config.yaml's
# intervals.bed points at one of them.
CAPTURE_BED_TWIST    = "resources/intervals/twist_exome_v2.GRCh38.bed"
CAPTURE_BED_ILLUMINA = "resources/intervals/illumina_exome_v2.5.GRCh38.bed"
CAPTURE_BED_IDT      = "resources/intervals/idt_xgen_exome_v2.GRCh38.bed"
CAPTURE_BED_GRCH38   = "resources/intervals/nextera_expandedexome.GRCh38.bed"

localrules:
    capture_bed_twist,
    capture_bed_illumina,
    capture_bed_idt,
    download_liftover_chain,
    download_nextera_hg19_bed,
    liftover_capture_bed,
    get_capture_bed,


# Shared recipe for native-GRCh38 exome BEDs (Twist Exome 2.0, the Twist-built
# Illumina Exome 2.5, IDT xGen Exome v2). Each ships chr-named with extra columns
# and some (IDT) carry a few overlapping targets. Trim to 3 cols, keep primary
# contigs (chr1..22,X,Y,M) so it matches the reference .dict, sort, and merge
# overlapping/abutting intervals (portable awk == `bedtools merge -d 0`), so the
# stored BED is clean regardless of the vendor's formatting.
_NATIVE_BED_SHELL = r"""
( tmp=$(mktemp -d)
  wget --tries=3 -q -O "$tmp/raw.bed" {params.url}
  cut -f1-3 "$tmp/raw.bed" \
    | awk '$1 ~ /^chr([0-9]+|X|Y|M)$/' \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{{OFS="\t"}} {{ if($1==c && $2<=e){{ if($3>e) e=$3 }} else {{ if(c!="") print c,s,e; c=$1; s=$2; e=$3 }} }} END{{ if(c!="") print c,s,e }}' \
    > {output}
  echo "raw=$(wc -l < "$tmp/raw.bed")  merged_targets=$(wc -l < {output})  covered_bp=$(awk '{{s+=$3-$2}} END{{print s}}' {output})" >&2
  rm -rf "$tmp" ) 2> {log}
"""


# --- Provider (1): Twist Exome 2.0 (native GRCh38, no liftover) ---------------
rule capture_bed_twist:
    output:
        CAPTURE_BED_TWIST,
    params:
        url=config["capture_bed"]["twist_url"],
    log:
        "logs/capture_bed/twist.log",
    conda:
        "../envs/download.yaml"
    shell:
        _NATIVE_BED_SHELL


# --- Provider (2): Illumina Exome 2.5 (Twist-built; native GRCh38) ------------
# ~99% identical to Twist 2.0 (same design); 2.5 is newer / ~1 Mb larger.
rule capture_bed_illumina:
    output:
        CAPTURE_BED_ILLUMINA,
    params:
        url=config["capture_bed"]["illumina_url"],
    log:
        "logs/capture_bed/illumina.log",
    conda:
        "../envs/download.yaml"
    shell:
        _NATIVE_BED_SHELL


# --- Provider (3): IDT xGen Exome Research Panel v2 (native GRCh38) -----------
# Independent design (~87% interval overlap with Twist, not 99%) — the useful
# "different vendor" option for cross-kit validation. Public IDT CDN BED6.
rule capture_bed_idt:
    output:
        CAPTURE_BED_IDT,
    params:
        url=config["capture_bed"]["idt_url"],
    log:
        "logs/capture_bed/idt.log",
    conda:
        "../envs/download.yaml"
    shell:
        _NATIVE_BED_SHELL


# --- Provider (4): Nextera Expanded Exome (hg19 BED -> GRCh38 liftover) -------
rule download_liftover_chain:
    output:
        "resources/liftover/hg19ToHg38.over.chain.gz",
    params:
        url=config["capture_bed"]["chain_url"],
    log:
        "logs/capture_bed/chain.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule download_nextera_hg19_bed:
    output:
        "resources/intervals/nextera_expandedexome.hg19.bed.gz",
    params:
        url=config["capture_bed"]["hg19_url"],
    log:
        "logs/capture_bed/nextera_hg19.log",
    conda:
        "../envs/download.yaml"
    shell:
        "wget --tries=3 -q -O {output} {params.url} 2> {log}"


rule liftover_capture_bed:
    input:
        bed="resources/intervals/nextera_expandedexome.hg19.bed.gz",
        chain="resources/liftover/hg19ToHg38.over.chain.gz",
    output:
        CAPTURE_BED_GRCH38,
    log:
        "logs/capture_bed/liftover.log",
    conda:
        "../envs/liftover.yaml"
    shell:
        # Keep only chr-prefixed data lines, lift, restrict to primary contigs
        # (chr1..22,X,Y,M) so it matches the reference .dict, then merge so
        # split-lifts don't leave overlapping/adjacent intervals.
        r"""
        ( tmp=$(mktemp -d)
          gunzip -c {input.bed} | grep '^chr' > "$tmp/hg19.bed"
          gunzip -c {input.chain} > "$tmp/chain"
          liftOver "$tmp/hg19.bed" "$tmp/chain" "$tmp/hg38.bed" "$tmp/unmapped.bed"
          awk '$1 ~ /^chr([0-9]+|X|Y|M)$/' "$tmp/hg38.bed" | sort -k1,1 -k2,2n > "$tmp/hg38.sorted.bed"
          bedtools merge -i "$tmp/hg38.sorted.bed" > {output}
          echo "hg19=$(wc -l < "$tmp/hg19.bed")  lifted_primary=$(wc -l < "$tmp/hg38.sorted.bed")  merged=$(wc -l < {output})" >&2
          rm -rf "$tmp" ) 2> {log}
        """


# Convenience target — build just the ACTIVE capture BED (whichever provider
# config["intervals"]["bed"] points at: Twist by default, or the Nextera lift).
rule get_capture_bed:
    input:
        config["intervals"]["bed"],

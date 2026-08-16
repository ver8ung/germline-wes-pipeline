#!/usr/bin/env bash
# Collect a benchmark run's results AND the provenance needed to interpret them
# into benchmarks/<label>/, ready to commit.
#
# WHY THIS EXISTS. This repo previously carried hap.py numbers that could not be
# tied to anything: the scoring engine was unpinned and its version unrecoverable
# (hap.py's own runinfo.json leaves the version field empty), the capture BED and
# unit sheet in use at the time were not recorded, and the code had moved on by
# the time anyone looked. A benchmark number without that context is not evidence,
# it is decoration. This script captures the context at the moment of collection,
# so the same thing cannot happen again.
#
# Run it on the HOST (it needs git), from the repo root, after a run + benchmark:
#
#   bash dev/collect_benchmark_evidence.sh                    # default run, all samples
#   bash dev/collect_benchmark_evidence.sh twist_onso         # a labelled run
#   bash dev/collect_benchmark_evidence.sh giab_trio HG002    # one sample of a run
#
# With one sample the layout is benchmarks/<label>/; with several it is
# benchmarks/<label>/<sample>/ plus a shared PROVENANCE.md index.
set -euo pipefail

LABEL="${1:-}"
shift || true
REQUESTED_SAMPLES=("$@")

if [ -n "$LABEL" ]; then
  PREFIX_BASE="results/benchmark/happy.${LABEL}"
  OUTROOT="benchmarks/${LABEL}"
  CONFIGFILE="config/${LABEL}.yaml"
else
  PREFIX_BASE="results/benchmark/happy"
  OUTROOT="benchmarks/default"
  CONFIGFILE=""
fi

# --- Which samples? ----------------------------------------------------------
# hap.py output is now per sample (happy[.label].<sample>.summary.csv), so with no
# explicit argument, collect every sample this run actually scored.
if [ ${#REQUESTED_SAMPLES[@]} -eq 0 ]; then
  mapfile -t REQUESTED_SAMPLES < <(
    for f in "${PREFIX_BASE}".*.summary.csv; do
      [ -f "$f" ] || continue
      s="${f#"${PREFIX_BASE}".}"
      echo "${s%.summary.csv}"
    done | sort -u
  )
fi

if [ ${#REQUESTED_SAMPLES[@]} -eq 0 ]; then
  echo "FATAL: no ${PREFIX_BASE}.<sample>.summary.csv found." >&2
  echo "       Run the pipeline and the 'benchmark' target first." >&2
  exit 1
fi

# --- Capture git state BEFORE writing anything ------------------------------
# Order matters. benchmarks/ is tracked (the whole point is that it is committed),
# so creating OUTROOT below makes `git status --porcelain` non-empty. Reading it
# afterwards stamped EVERY provenance record as dirty, including runs from a
# pristine checkout -- turning the one field that certifies the numbers into a
# constant, and making a genuinely dirty run indistinguishable from a clean one.
SHA="$(git rev-parse HEAD 2>/dev/null || echo 'not a git checkout')"
DIRTY=""
if [ -n "$(git status --porcelain -- . ':(exclude)benchmarks' 2>/dev/null)" ]; then
  DIRTY="  **(working tree was DIRTY -- these numbers do not describe a clean commit)**"
fi

IMG="ngs-germline-wes:latest"

# --- Metadata, parsed properly rather than grepped --------------------------
# Everything below is read from the artifacts of the run that ACTUALLY happened
# (hap.py's runinfo.json) rather than re-derived from config, so a config edit
# between running and collecting cannot silently rewrite history. The one
# exception is the truth-set version: local truth paths are now
# resources/benchmark/<sample>.truth.* by design, so the version string lives only
# in the remote filename in config.
#
# Parsed with real yaml/json inside the image rather than grep/awk on the host:
# the truth block is nested three levels deep, and a `grep -E '^\s+vcf:'` would
# happily return another sample's value after any reformatting.
read_meta() {   # $1 = sample
  MSYS_NO_PATHCONV=1 docker run --rm -i -v "$(pwd):/workflow" -w /workflow "$IMG" \
    python - "$1" "$PREFIX_BASE" "$CONFIGFILE" <<'PY' 2>/dev/null || true
import json, os, sys, yaml

sample, prefix_base, overlay = sys.argv[1], sys.argv[2], sys.argv[3]
out = {}

runinfo = f"{prefix_base}.{sample}.runinfo.json"
if os.path.exists(runinfo):
    with open(runinfo) as fh:
        ri = json.load(fh)
    fa = ri.get("final_args", {}) or {}
    out["ENGINE"] = fa.get("engine") or ""
    out["SDF"] = fa.get("engine_vcfeval_template") or ""
    out["STRAT"] = fa.get("strat_tsv") or ""
    out["TRUTH_BED"] = fa.get("fp_bedfile") or ""
    # hap.py records neither the truth VCF nor -T in final_args; both are only in
    # the recorded command line.
    cmd = ""
    for entry in ri.get("runInfo", []) or []:
        v = entry.get("value", "")
        if "hap.py" in v:
            cmd = v
            break
    out["CMDLINE"] = " ".join(cmd.split())

cfg = {}
with open("config/config.yaml") as fh:
    cfg = yaml.safe_load(fh) or {}
if overlay and os.path.exists(overlay):
    with open(overlay) as fh:
        ov = yaml.safe_load(fh) or {}
    def merge(a, b):
        for k, v in b.items():
            if isinstance(v, dict) and isinstance(a.get(k), dict):
                merge(a[k], v)
            else:
                a[k] = v
    merge(cfg, ov)

out["UNITS"] = cfg.get("units", "")
out["BED"] = (cfg.get("intervals", {}) or {}).get("bed", "")
bench = cfg.get("benchmark", {}) or {}
out["TARGETS"] = bench.get("targets_bed", out["BED"])
truth = (bench.get("truth", {}) or {}).get(sample, {}) or {}
# The version-bearing name, e.g. HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
out["TRUTH_REMOTE"] = truth.get("vcf", "")
out["TRUTH_URL"] = (truth.get("base_url", "") or "") + (truth.get("vcf", "") or "")
strat = bench.get("stratification", {}) or {}
out["STRAT_VERSION"] = str(strat.get("version", "")) if strat.get("activate") else ""
out["STRAT_N"] = str(len(strat.get("regions", {}) or {})) if strat.get("activate") else "0"

for k, v in out.items():
    print(f"{k}={v}")
PY
}

# Tool versions, read from the baked conda envs rather than from runinfo.json,
# whose version field hap.py leaves empty. This is the field whose absence made
# the previous numbers unreproducible.
pkg_version() {  # $1 = conda package name prefix
  MSYS_NO_PATHCONV=1 docker run --rm "$IMG" bash -c \
    "ls /opt/snakemake-envs/*/conda-meta/$1-*.json 2>/dev/null | head -1 | xargs -r basename" \
    2>/dev/null | sed 's/\.json$//' || true
}
HAPPY_VER="$(pkg_version 'hap.py')"
RTG_VER="$(pkg_version 'rtg-tools')"
[ -n "$HAPPY_VER" ] || HAPPY_VER="UNKNOWN (could not read the baked env)"

MULTI=0
[ ${#REQUESTED_SAMPLES[@]} -gt 1 ] && MULTI=1

collect_one() {   # $1 = sample, $2 = output dir
  local sample="$1" outdir="$2" prefix="${PREFIX_BASE}.${1}"

  [ -f "${prefix}.summary.csv" ] || {
    echo "FATAL: ${prefix}.summary.csv not found." >&2
    exit 1
  }

  mkdir -p "$outdir/qc"
  cp "${prefix}.summary.csv" "$outdir/happy.summary.csv"
  [ -f "${prefix}.extended.csv" ] && cp "${prefix}.extended.csv" "$outdir/happy.extended.csv"
  [ -f "${prefix}.runinfo.json" ] && cp "${prefix}.runinfo.json" "$outdir/happy.runinfo.json"

  # QC: the small per-sample metrics. Scoped to THIS sample rather than globbed --
  # results/ is shared across configs, so a glob would sweep in another cohort's
  # leftovers and attribute them to this run.
  for f in "results/qc/samtools/${sample}.stats.txt" \
           "results/qc/samtools/${sample}.flagstat.txt" \
           "results/qc/dedup/${sample}.metrics.txt"; do
    [ -f "$f" ] && cp "$f" "$outdir/qc/" || true
  done

  # Cohort-level artifacts live at the run root, not under a sample's qc/ -- they
  # describe the whole callset, and in a multi-sample run copying them into every
  # sample directory would duplicate them. multiqc_report.html was already here;
  # bcftools stats joins it. Never multiqc_data.json (~10 MB).
  [ -f results/qc/bcftools/all.stats.txt ] && cp results/qc/bcftools/all.stats.txt "$OUTROOT/bcftools_all.stats.txt" || true
  [ -f results/qc/multiqc_report.html ] && cp results/qc/multiqc_report.html "$OUTROOT/" || true

  # --- Provenance -----------------------------------------------------------
  local meta; meta="$(read_meta "$sample")"
  get() { echo "$meta" | grep -m1 "^$1=" | cut -d= -f2- || true; }

  local engine sdf strat truth_bed truth_remote truth_url units bed targets strat_ver strat_n cmdline
  engine="$(get ENGINE)";           sdf="$(get SDF)"
  strat="$(get STRAT)";             truth_bed="$(get TRUTH_BED)"
  truth_remote="$(get TRUTH_REMOTE)"; truth_url="$(get TRUTH_URL)"
  units="$(get UNITS)";             bed="$(get BED)"
  targets="$(get TARGETS)";         strat_ver="$(get STRAT_VERSION)"
  strat_n="$(get STRAT_N)";         cmdline="$(get CMDLINE)"

  [ -n "$units" ] || { echo "FATAL: could not read 'units' from the merged config." >&2; exit 1; }
  [ -n "$bed" ]   || { echo "FATAL: could not read 'intervals.bed' from the merged config." >&2; exit 1; }

  local engine_note="${engine:-UNKNOWN}"
  local strat_row="none"
  # $strat_ver already carries its own leading 'v' (v3.6), so do not add another.
  [ -n "$strat_ver" ] && strat_row="GIAB genome stratifications ${strat_ver}, ${strat_n} curated subsets"

  {
    echo "# Benchmark provenance -- ${LABEL:-default} run, sample ${sample}"
    echo
    echo "Generated by \`dev/collect_benchmark_evidence.sh\`."
    echo
    echo "| Field | Value |"
    echo "|---|---|"
    echo "| Commit | \`${SHA}\`${DIRTY} |"
    echo "| Sample | \`${sample}\` |"
    echo "| Config overlay | ${CONFIGFILE:-none (default config)} |"
    echo "| Units sheet | \`${units}\` |"
    echo "| Capture BED (calling) | \`${bed}\` |"
    echo "| Evaluation BED (hap.py -T) | \`${targets}\` |"
    echo "| Truth set | \`${truth_remote:-UNKNOWN}\` |"
    echo "| Truth confident regions | \`${truth_bed:-UNKNOWN}\` |"
    echo "| Truth source | ${truth_url:-UNKNOWN} |"
    echo "| Scoring tool | \`${HAPPY_VER}\` |"
    echo "| Comparison engine | \`${engine_note}\`${RTG_VER:+ (\`${RTG_VER}\`)} |"
    [ -n "$sdf" ] && echo "| vcfeval template | \`${sdf}\` |"
    echo "| Stratification | ${strat_row} |"
    echo "| Evaluation region | capture targets ∩ GIAB high-confidence |"
    echo
    echo "## Results"
    echo
    echo '```'
    cat "$outdir/happy.summary.csv"
    echo '```'
    echo
    echo "## Reading these numbers"
    echo
    echo "- **PASS** rows score only records that passed hard filtering; **ALL** rows"
    echo "  score every record. Both matter: a large ALL-to-PASS recall drop means the"
    echo "  filters are discarding true variants."
    echo "- Recall here is bounded by capture coverage, not only by the caller. Truth"
    echo "  variants in target regions the capture never covered count as false"
    echo "  negatives however good the calling is. Use \`dev/diag_capture.sh\` to split"
    echo "  the false negatives into covered and uncovered before drawing conclusions."
    echo "- The evaluation region is the capture BED intersected with the GIAB"
    echo "  high-confidence regions, so numbers are NOT comparable across runs that"
    echo "  used different capture kits or different truth sets -- the denominators"
    echo "  differ. Compare ALL-to-PASS BEHAVIOUR across runs, not absolute recall."
    if [ -n "$strat_ver" ]; then
      echo "- With stratification on, \`Subset=*\` is the headline row and the named"
      echo "  subsets are diagnostic: they say WHERE calling fails, not just how often."
      echo "  See \`happy.extended.csv\`."
    fi
    if [ "$engine_note" = "vcfeval" ]; then
      echo "- Scored with the **vcfeval** (RTG) comparison engine rather than hap.py's"
      echo "  default \`xcmp\`. Measured on identical input, the two agree exactly on SNPs"
      echo "  and to within ~0.001 F1 on indels: hap.py already left-shifts and decomposes"
      echo "  both callsets before comparing (\`preprocessing_leftshift\`), so most"
      echo "  variant-representation differences are normalised before either engine sees"
      echo "  them. vcfeval is used because it is the stricter, more widely cited engine,"
      echo "  NOT because it materially changes these numbers."
    fi
    if [ -n "$cmdline" ]; then
      echo
      echo "## Exact command"
      echo
      echo '```'
      echo "$cmdline"
      echo '```'
    fi
  } > "$outdir/PROVENANCE.md"
}

for s in "${REQUESTED_SAMPLES[@]}"; do
  if [ "$MULTI" -eq 1 ]; then
    collect_one "$s" "$OUTROOT/$s"
  else
    collect_one "$s" "$OUTROOT"
  fi
done

# Mendelian violation rate, when the run produced one (trio configs only).
MENDEL="results/benchmark/mendelian${LABEL:+.$LABEL}.txt"
[ -f "$MENDEL" ] && cp "$MENDEL" "$OUTROOT/mendelian.txt" || true

echo "Collected into $OUTROOT:"
find "$OUTROOT" -type f | sort | sed 's/^/  /'
echo
echo "samples:       ${REQUESTED_SAMPLES[*]}"
echo "scoring tool:  $HAPPY_VER"
[ -n "$RTG_VER" ] && echo "rtg-tools:     $RTG_VER"
[ -n "$DIRTY" ] && echo "WARNING: working tree is dirty; commit before trusting these numbers." >&2
exit 0

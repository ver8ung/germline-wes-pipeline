#!/usr/bin/env bash
# Run the NGS germline WES pipeline inside the Docker orchestrator (Linux/Mac/WSL).
#
# Builds a commit-tagged image on first use, then forwards all args to Snakemake.
# Image tag = git short sha (+ -dirty), or a content hash of the image inputs
# (Dockerfile + workflow/envs/*.yaml) when this isn't a git checkout. The base
# image is pinned via docker/base-image.txt. A staleness guard warns before each
# run if the env files changed since the image was built.
#
#   ./run.sh -n                  # dry-run the whole DAG
#   ./run.sh smoke -n            # dry-run the tiny smoke test (no downloads)
#   ./run.sh smoke               # full smoke test: simulate reads -> annotated VCF
#   ./run.sh setup_reference     # only download + index the reference
#   ./run.sh --cores 8           # full run, FASTQ -> annotated VCF
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="ngs-germline-wes"

# --- Image tag: content hash of exactly the files that determine image content ---
# Deliberately NOT the git commit SHA. The image is built from the Dockerfile, the
# per-rule env YAMLs and the base image tag, and nothing else; tagging by commit
# meant every commit invalidated the tag and triggered a needless 15-30 min
# rebuild, while two commits with identical image inputs pointlessly built twice.
#
# MUST match Get-ContentTag in run.ps1 byte for byte, or the same tree builds two
# different images on Windows vs Linux. Shared contract: sort inputs by relative
# POSIX path, cat their raw bytes in that order, one sha256, lowercase, first 12.
# base-image.txt is included because it is fed to the build as --build-arg BASE;
# leaving it out meant re-pinning the base silently reused a stale image.
tag="src-$(
  {
    find "$PROJ/workflow/envs" -name '*.yaml' -print0 | sort -z |
      xargs -0 cat
    cat "$PROJ/docker/Dockerfile" "$PROJ/docker/base-image.txt"
  } | sha256sum | cut -c1-12
)"
IMAGE="$REPO:$tag"

# --- Build this exact tag if it isn't present yet ---
if [[ -z "$(docker images -q "$IMAGE")" ]]; then
  base="$(grep -vE '^\s*#' "$PROJ/docker/base-image.txt" | grep -E '\S' | head -1 | tr -d '[:space:]')"
  echo "==> Building $IMAGE   (base: $base)"
  docker build -t "$IMAGE" -t "$REPO:latest" \
    --build-arg BASE="$base" --build-arg IMAGE_TAG="$tag" \
    -f "$PROJ/docker/Dockerfile" "$PROJ"
fi

# --- Staleness guard: warn if env files drifted from what the image baked ---
if ! docker run --rm -v "$PROJ:/workflow" "$IMAGE" check-bake >/dev/null 2>&1; then
  echo "WARNING: workflow/envs/*.yaml changed since image '$IMAGE' was built." >&2
  echo "         Per-rule envs would be rebuilt every run (slow, needs network)." >&2
  echo "         Re-bake them:  docker rmi $IMAGE ; then re-run this script." >&2
fi

# --- `smoke` shortcut: apply the smoke config overlay, forward the rest (e.g. -n) ---
if [[ "${1:-}" == "smoke" ]]; then
  shift
  set -- --configfile .tests/smoke/config.yaml --cores 4 "$@"
fi

# --- Default to a full run on 8 cores if nothing was passed ---
if [[ $# -eq 0 ]]; then
  set -- --cores 8
fi

# --- Only allocate a TTY when stdin actually is one. Hardcoding -it breaks any
# --- redirected use, e.g. the documented `./run.sh --dag | dot -Tsvg > dag.svg`,
# --- which Docker aborts with "the input device is not a TTY".
tty_args=()
if [[ -t 0 ]]; then
  tty_args=(-it)
fi

# Memory limit, overridable. Must fit inside the Docker Desktop / WSL2 VM limit.
DOCKER_MEM="${WES_DOCKER_MEMORY:-16g}"

# --- Translate the container limit into a Snakemake mem_mb budget -------------
# The per-rule `resources: mem_mb` values in config.yaml are ADVISORY unless the
# run also passes --resources mem_mb=<budget>. Without it Snakemake schedules
# purely on cores, so two GATK jobs that each reserved 8 GB can be co-scheduled
# inside a 16 GB container alongside the resident bwa index, and the OOM killer
# decides the outcome. 90% leaves headroom for Snakemake itself.
case "$DOCKER_MEM" in
  *[gG]) mem_mib=$(( ${DOCKER_MEM%[gG]} * 1024 )) ;;
  *[mM]) mem_mib=${DOCKER_MEM%[mM]} ;;
  *)     mem_mib=$(( 16 * 1024 )) ;;
esac
mem_budget=$(( mem_mib * 9 / 10 ))

# Preflight: a per-rule mem_mb larger than the whole budget is a hard Snakemake
# error that surfaces only once that job becomes schedulable, with a message that
# blames pipes. Catch it here, while we can still say what to do about it.
max_rule_mem=$(grep -oE '^[[:space:]]+mem_mb:[[:space:]]*[0-9]+' "$PROJ/config/config.yaml" |
               grep -oE '[0-9]+$' | sort -n | tail -1)
if [[ -n "${max_rule_mem:-}" ]] && (( mem_budget < max_rule_mem )); then
  need_gib=$(( (max_rule_mem * 10 / 9 + 1023) / 1024 ))
  {
    echo "ERROR: container memory '$DOCKER_MEM' is too small for this config."
    echo "       Budget would be ${mem_budget} MiB (90% of ${mem_mib} MiB), but the largest"
    echo "       rule in config/config.yaml reserves ${max_rule_mem} MiB."
    echo "       Set WES_DOCKER_MEMORY to at least ${need_gib}g, or lower resources.*.mem_mb."
  } >&2
  exit 1
fi

# --resources must go LAST: it takes a variable-length list, so placing it before
# the caller's arguments makes it swallow a positional target name.
extra_args=()
if [[ " $* " != *" --resources "* ]]; then
  extra_args=(--resources "mem_mb=$mem_budget")
fi

# Banner goes to stderr: on stdout it would be piped into whatever consumes the
# run's output (again, `--dag | dot`) and corrupt it.
echo "==> [$IMAGE] snakemake --use-conda $* ${extra_args[*]}" >&2
exec docker run --rm "${tty_args[@]}" --init \
  -v "$PROJ:/workflow" \
  -w /workflow \
  "--memory=$DOCKER_MEM" \
  "$IMAGE" \
  snakemake --use-conda --conda-prefix /opt/snakemake-envs "$@" "${extra_args[@]}"

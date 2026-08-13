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

# --- Image tag: git short sha (+ -dirty), else content hash of image inputs ---
tag="$(git -C "$PROJ" rev-parse --short=12 HEAD 2>/dev/null || true)"
if [[ -n "$tag" ]]; then
  [[ -n "$(git -C "$PROJ" status --porcelain 2>/dev/null)" ]] && tag="${tag}-dirty"
  tag="git-${tag}"
else
  # shellcheck disable=SC2012
  tag="src-$(cat $(ls "$PROJ"/workflow/envs/*.yaml | sort) "$PROJ/docker/Dockerfile" | sha256sum | cut -c1-12)"
fi
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

echo "==> [$IMAGE] snakemake --use-conda $*"
exec docker run --rm -it --init \
  -v "$PROJ:/workflow" \
  -w /workflow \
  --memory=16g \
  "$IMAGE" \
  snakemake --use-conda --conda-prefix /opt/snakemake-envs "$@"

#requires -Version 5
<#
.SYNOPSIS
  Run the NGS germline WES pipeline inside the Docker orchestrator (Windows).

.DESCRIPTION
  Builds a commit-tagged image on first use, then runs Snakemake with --use-conda
  inside a container that has the project bind-mounted at /workflow. Any arguments
  you pass are forwarded to Snakemake.

  Image tag = git short sha (+ -dirty), or a content hash of the image inputs
  (Dockerfile + workflow/envs/*.yaml) when this isn't a git checkout. Each distinct
  tag builds its own image, so image <-> code stay in lockstep. The base image is
  pinned via docker/base-image.txt. Before each run, a staleness guard warns if the
  env files changed since the image was built.

.EXAMPLE
  .\run.ps1 -n                       # dry-run the whole DAG (nothing executes)
  .\run.ps1 smoke -n                 # dry-run the tiny smoke test (no downloads)
  .\run.ps1 smoke                    # full smoke test: simulate reads -> annotated VCF
  .\run.ps1 setup_reference          # only download + index the reference
  .\run.ps1 --cores 8                # full run, FASTQ -> annotated VCF
#>
[CmdletBinding()]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]] $SnakemakeArgs
)

$ErrorActionPreference = 'Stop'
$proj = Split-Path -Parent $MyInvocation.MyCommand.Path
$repo = 'ngs-germline-wes'

# --- Image tag: git short sha (+ -dirty), else content hash of image inputs ---
function Get-ContentTag {
    $files = @(Join-Path $proj 'docker\Dockerfile')
    $files += Get-ChildItem (Join-Path $proj 'workflow\envs\*.yaml') |
        Sort-Object FullName | Select-Object -ExpandProperty FullName
    $acc = ($files | ForEach-Object { (Get-FileHash -Algorithm SHA256 -Path $_).Hash }) -join ''
    $bytes = [System.Security.Cryptography.SHA256]::Create().ComputeHash([Text.Encoding]::ASCII.GetBytes($acc))
    return 'src-' + ((($bytes | ForEach-Object { $_.ToString('x2') }) -join '').Substring(0, 12))
}

$tag = $null
try {
    $sha = (& git -C $proj rev-parse --short=12 HEAD 2>$null)
    if ($LASTEXITCODE -eq 0 -and $sha) {
        # Capture porcelain output first (a plain assignment); testing the native
        # command's `2>$null` redirection directly inside `if (...)` trips
        # PSScriptAnalyzer's PSPossibleIncorrectUsageOfRedirectionOperator (a
        # false positive - `>` here redirects stderr, it is not a comparison).
        $porcelain = (& git -C $proj status --porcelain 2>$null)
        $dirty = if ($porcelain) { '-dirty' } else { '' }
        $tag = "git-$($sha.Trim())$dirty"
    }
} catch {}
if (-not $tag) { $tag = Get-ContentTag }
$image = "${repo}:$tag"

# --- Build this exact tag if it isn't present yet ---
if (-not (docker images -q $image)) {
    $base = (Get-Content (Join-Path $proj 'docker\base-image.txt') |
        Where-Object { $_ -and ($_ -notmatch '^\s*#') } | Select-Object -First 1).Trim()
    Write-Host "==> Building $image   (base: $base)" -ForegroundColor Cyan
    docker build -t $image -t "${repo}:latest" `
        --build-arg BASE=$base --build-arg IMAGE_TAG=$tag `
        -f (Join-Path $proj 'docker\Dockerfile') $proj
    if ($LASTEXITCODE -ne 0) { throw "docker build failed" }
}

# --- Staleness guard: warn if env files drifted from what the image baked ---
docker run --rm -v "${proj}:/workflow" $image check-bake 2>$null | Out-Null
if ($LASTEXITCODE -ne 0) {
    Write-Host "WARNING: workflow/envs/*.yaml changed since image '$image' was built." -ForegroundColor Yellow
    Write-Host "         Per-rule envs would be rebuilt every run (slow, needs network)." -ForegroundColor Yellow
    Write-Host "         Re-bake them:  docker rmi $image ; then re-run this script." -ForegroundColor Yellow
}

# --- `smoke` shortcut: apply the smoke config overlay, forward the rest (e.g. -n) ---
if ($SnakemakeArgs.Count -ge 1 -and $SnakemakeArgs[0] -eq 'smoke') {
    $rest = @()
    if ($SnakemakeArgs.Count -gt 1) { $rest = $SnakemakeArgs[1..($SnakemakeArgs.Count - 1)] }
    $SnakemakeArgs = @('--configfile', '.tests/smoke/config.yaml', '--cores', '4') + $rest
}

# --- Default to a full run on 8 cores if the caller passed nothing ---
if (-not $SnakemakeArgs -or $SnakemakeArgs.Count -eq 0) {
    $SnakemakeArgs = @('--cores', '8')
}

Write-Host "==> [$image] snakemake --use-conda $($SnakemakeArgs -join ' ')" -ForegroundColor Cyan
# Allocate an interactive TTY only when STDIN is a real console. In a
# non-interactive context (CI, background job) `-t` fails with
# "the input device is not a TTY", so fall back to a plain batch run.
$ttyArgs = if ([Console]::IsInputRedirected) { @() } else { @('-it') }
# Container memory cap. Default 16g; override for big multi-unit runs that align
# several libraries in parallel, e.g.  $env:WES_DOCKER_MEMORY = '40g'
$dockerMem = if ($env:WES_DOCKER_MEMORY) { $env:WES_DOCKER_MEMORY } else { '16g' }

# Optional fast scratch for the heavy intermediate BAMs. Docker Desktop bind
# mounts (the E: drive) write at only ~6 MB/s, and the genome-scale BAMs in
# results/{mapped,dedup,bqsr} dominate wall-clock. Mounting Docker-managed named
# volumes over just those three subdirs moves their I/O onto fast VM-native
# storage, while the rest of results/ stays on the bind mount so the deliverables
# (VCFs, benchmark, QC/MultiQC) stay visible on the host. Enable with:
#     $env:WES_SCRATCH_BAMS = '1'
# Those dirs hold only regenerable intermediates (recal.bam + temp mapped/dedup
# BAMs), never final outputs. Toggling this changes where Snakemake sees those
# files, so it re-runs mapping->BQSR, so use it for fresh runs. Tools that read
# results/bqsr on the host (e.g. diag_capture.sh) must mount wes-bqsr too.
# Reset/reclaim:  docker volume rm wes-mapped wes-dedup wes-bqsr
$scratchArgs = @()
if ($env:WES_SCRATCH_BAMS) {
    foreach ($d in 'mapped', 'dedup', 'bqsr') {
        $scratchArgs += @('--mount', "type=volume,source=wes-$d,target=/workflow/results/$d,volume-nocopy")
    }
    Write-Host "==> scratch volumes ON: results/{mapped,dedup,bqsr} -> wes-{mapped,dedup,bqsr}" -ForegroundColor Cyan
}

docker run --rm @ttyArgs --init `
    -v "${proj}:/workflow" `
    @scratchArgs `
    -w /workflow `
    --memory=$dockerMem `
    $image `
    snakemake --use-conda --conda-prefix /opt/snakemake-envs @SnakemakeArgs

exit $LASTEXITCODE

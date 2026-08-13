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
# MUST match the content-tag branch of run.sh byte for byte, or the same tree
# builds two different images on Windows vs Linux. The previous version could
# never agree with run.sh: it hashed each file separately, joined the UPPERCASE
# hex digests and hashed that string, while run.sh cat'd the raw bytes and hashed
# once -- mathematically unrelated functions -- and it also ordered the Dockerfile
# first where run.sh put it last.
#
# Shared contract: sort inputs by relative POSIX path, concatenate their raw
# bytes in that order, one SHA256 over the result, lowercase, first 12 chars.
# base-image.txt is included because it is fed to the build as --build-arg BASE;
# omitting it meant re-pinning the base silently reused a stale image.
function Get-ContentTag {
    $files = @()
    $files += Get-ChildItem (Join-Path $proj 'workflow\envs\*.yaml') |
        Sort-Object Name | Select-Object -ExpandProperty FullName
    $files += (Join-Path $proj 'docker\Dockerfile')
    $files += (Join-Path $proj 'docker\base-image.txt')

    $stream = New-Object System.IO.MemoryStream
    foreach ($f in $files) {
        $bytes = [System.IO.File]::ReadAllBytes($f)
        $stream.Write($bytes, 0, $bytes.Length)
    }
    $stream.Position = 0
    $hash = (Get-FileHash -InputStream $stream -Algorithm SHA256).Hash.ToLower()
    $stream.Dispose()
    return 'src-' + $hash.Substring(0, 12)
}

# Deliberately NOT the git commit SHA -- see the comment on Get-ContentTag.
$tag = Get-ContentTag
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
# $ErrorActionPreference='Stop' plus `2>$null` on a NATIVE exe is a trap in
# PS 5.1: each stderr line comes back as a NativeCommandError that terminates the
# script even when the command exited 0. Any incidental docker chatter ("Unable
# to find image ... locally", a WSL2 integration notice) would kill the launcher
# before Snakemake ever started, and the warning branch below would be
# unreachable. Drop to Continue for the probe, then restore.
$prevEAP = $ErrorActionPreference
$ErrorActionPreference = 'Continue'
docker run --rm -v "${proj}:/workflow" $image check-bake 2>$null | Out-Null
$bakeExit = $LASTEXITCODE
$ErrorActionPreference = $prevEAP
if ($bakeExit -ne 0) {
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

# --- Translate the container limit into a Snakemake mem_mb budget -------------
# The per-rule `resources: mem_mb` values in config.yaml are ADVISORY unless the
# run also passes --resources mem_mb=<budget>. Without it Snakemake schedules
# purely on cores, so two GATK jobs that each reserved 8 GB can be co-scheduled
# inside a 16 GB container alongside the resident bwa index, and the OOM killer
# decides the outcome. 90% leaves headroom for Snakemake itself.
# Validate rather than silently defaulting. Docker accepts spellings this parser
# did not recognise ('8gb', '1.5g'); falling back to 16 GiB then told Snakemake it
# could co-schedule ~14.7 GB inside an 8 GB container -- the exact OOM the budget
# exists to prevent, with the preflight below reporting all clear.
if ($dockerMem -notmatch '^(\d+)\s*([gGmM])[bB]?$') {
    Write-Host "ERROR: WES_DOCKER_MEMORY='$dockerMem' is not a form this launcher can size." -ForegroundColor Red
    Write-Host "       Use an integer with a g or m suffix, e.g. 16g, 32g, 24000m." -ForegroundColor Red
    exit 1
}
# Read both capture groups out of $Matches BEFORE any other -match runs. The
# -match operator overwrites $Matches on success, so testing the unit with a
# second -match clobbers the number: '16g' silently became a 0 MiB budget, while
# '24000m' survived only because a FAILED match leaves $Matches untouched.
$memNum = [int]$Matches[1]
$memUnit = $Matches[2].ToLower()
$memMib = if ($memUnit -eq 'g') { $memNum * 1024 } else { $memNum }
# Floor, not [int] rounding: run.sh computes mem_mib*9/10 in integer arithmetic,
# which truncates. [int] rounds, so 16g gave 14746 here against 14745 there.
$memBudget = [int][math]::Floor($memMib * 0.9)

# Preflight: a per-rule mem_mb larger than the whole budget is a hard Snakemake
# error, and it surfaces only once that job becomes schedulable -- i.e. after
# FastQC and trimming have already run -- with a message that blames pipes.
# Catch it here, while we can still say what to do about it.
$maxRuleMem = 0
Get-Content (Join-Path $proj 'config\config.yaml') |
    Select-String -Pattern '^\s+mem_mb:\s*(\d+)' |
    ForEach-Object { $v = [int]$_.Matches[0].Groups[1].Value; if ($v -gt $maxRuleMem) { $maxRuleMem = $v } }
if ($maxRuleMem -gt 0 -and $memBudget -lt $maxRuleMem) {
    $needGib = [math]::Ceiling($maxRuleMem / 0.9 / 1024)
    Write-Host "ERROR: container memory '$dockerMem' is too small for this config." -ForegroundColor Red
    Write-Host "       Budget would be ${memBudget} MiB (90% of $memMib MiB), but the largest" -ForegroundColor Red
    Write-Host "       rule in config/config.yaml reserves ${maxRuleMem} MiB." -ForegroundColor Red
    Write-Host "       Set WES_DOCKER_MEMORY to at least ${needGib}g, or lower resources.*.mem_mb." -ForegroundColor Red
    exit 1
}

# --resources must go LAST: it takes a variable-length list, so placing it before
# the caller's arguments makes it swallow a positional target name.
$extraArgs = @()
if (-not ($SnakemakeArgs -contains '--resources')) {
    $extraArgs += @('--resources', "mem_mb=$memBudget")
}

Write-Host "==> [$image] snakemake --use-conda $($SnakemakeArgs -join ' ') $($extraArgs -join ' ')" -ForegroundColor Cyan

docker run --rm @ttyArgs --init `
    -v "${proj}:/workflow" `
    @scratchArgs `
    -w /workflow `
    --memory=$dockerMem `
    $image `
    snakemake --use-conda --conda-prefix /opt/snakemake-envs @SnakemakeArgs @extraArgs

exit $LASTEXITCODE

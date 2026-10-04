# Windows version of bench_gto_gpu.sh: builds mdlib (Vulkan backend, Release) in
# build-bench\ and compares the GPU GTO kernels:
#   * electron density: the reference kernel, the tiled kernel and the two-pass GEMM
#     path in its configurations,
#   * molecular orbitals (1 orbital psi and psi^2, 32 orbitals sum of psi^2): the
#     reference kernel, the shell kernel in its configurations and the GEMM path.
# Every result is checked against the reference kernel and a CPU evaluation.
#
#   powershell -ExecutionPolicy Bypass -File .\bench_gto_gpu.ps1            full run
#   powershell -ExecutionPolicy Bypass -File .\bench_gto_gpu.ps1 --quick    smoke test
#
# Requirements: Visual Studio (C toolchain), CMake, Python 3 and git on PATH, network
# access on the first configure (Slang, the Vulkan headers and volk are downloaded), and a
# GPU driver with Vulkan support. No Vulkan SDK and no HDF5 needed. Run it from a
# "Developer PowerShell for VS" if CMake does not find Visual Studio by itself.
#
# Extra arguments are passed to md_bench_gto_gpu (--case NAME, --algos LIST, --dim N,
# --iters N, --seconds S, --scratch-mb MB); see benchmark\bench_gto_gpu.c.
#
# Choosing the GPU (e.g. the integrated Intel GPU instead of a discrete NVIDIA one):
#   .\bench_gto_gpu.ps1 -ListDevices                 lists the Vulkan adapters and exits
#   .\bench_gto_gpu.ps1 -Device intel                part of the adapter name ...
#   .\bench_gto_gpu.ps1 -Device 1                    ... or its number from -ListDevices
#   .\bench_gto_gpu.ps1 -Prefer low-power            integrated before discrete
# Without these, $env:MD_GPU_DEVICE is used when set, else the discrete GPU. A -Device that
# matches nothing is an error (the run does not fall back to another GPU).
#
# Results go to bench_gto_gpu_<host>.txt, or bench_gto_gpu_<host>_<device>.txt with
# -Device / -Prefer / MD_GPU_DEVICE, so runs on different GPUs keep separate logs. The log
# starts with a description of the system (OS, CPU, memory, compiler, mdlib revision, the
# adapters, the GPU used and its driver).

[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Device = "",
    [ValidateSet("", "high-performance", "low-power")][string]$Prefer = "",
    [switch]$ListDevices,
    [Parameter(ValueFromRemainingArguments = $true)][string[]]$BenchArgs = @()
)

# "Continue", not "Stop": Windows PowerShell 5.1 turns every stderr line of a redirected
# native command (cmake warnings, the benchmark's log) into an error record, which "Stop"
# would make fatal. Exit codes are checked explicitly instead.
$ErrorActionPreference = "Continue"
Set-Location -Path $PSScriptRoot

$build = "build-bench"

if (-not $BenchArgs) { $BenchArgs = @() }
if ($Device) { $BenchArgs = @("--device", $Device) + $BenchArgs }
if ($Prefer) { $BenchArgs = @("--prefer", $Prefer) + $BenchArgs }

# Log name: one per host, and per GPU selection.
$tag = ""
$sel = if ($Device) { $Device } elseif ($env:MD_GPU_DEVICE) { $env:MD_GPU_DEVICE } else { "" }
if ($sel)    { $tag += "_" + $sel }
if ($Prefer) { $tag += "_" + $Prefer }
$tag = $tag -replace "[^A-Za-z0-9_.-]", "-"
$log = Join-Path $PSScriptRoot "bench_gto_gpu_$($env:COMPUTERNAME)$tag.txt"

# Shows lines on the console and appends them to the log as UTF-8 (Tee-Object would
# write UTF-16 on Windows PowerShell 5.1 and has no -Encoding there).
function Write-Log {
    param([Parameter(ValueFromPipeline = $true)][AllowEmptyString()][string]$Line)
    process {
        Write-Host $Line
        Out-File -FilePath $log -Append -Encoding utf8 -InputObject $Line
    }
}

$cmakeArgs = @("-S", ".", "-B", $build, "-DMD_ENABLE_GPU=ON", "-DMD_ENABLE_HDF5=OFF", "-DMD_UNITTEST=OFF", "-DMD_BENCHMARK=ON")
# Reuse an already downloaded slang instead of fetching it again.
if (Test-Path "build\third_party\slang") {
    $cmakeArgs += "-DSLANG_CACHE_DIR=$(Join-Path $PSScriptRoot 'build\third_party\slang')"
}

Write-Host "Configuring $build ..."
& cmake @cmakeArgs *> "$build.configure.log"
if ($LASTEXITCODE -ne 0) {
    Get-Content "$build.configure.log" -Tail 40
    Write-Host "Configure failed, see $build.configure.log"
    exit 1
}

Write-Host "Building md_bench_gto_gpu ..."
& cmake --build $build --config Release --target md_bench_gto_gpu --parallel *> "$build.build.log"
if ($LASTEXITCODE -ne 0) {
    Select-String -Path "$build.build.log" -Pattern "error" | Select-Object -First 40
    Write-Host "Build failed, see $build.build.log"
    exit 1
}

# Prefer the Release binary (multi-config generators put each configuration in its own folder).
$bin = Get-ChildItem -Path $build -Recurse -Filter "md_bench_gto_gpu.exe" |
    Sort-Object @{ Expression = { if ($_.FullName -match "\\Release\\") { 0 } else { 1 } } }, @{ Expression = { $_.LastWriteTime }; Descending = $true } |
    Select-Object -First 1
if (-not $bin) { Write-Host "md_bench_gto_gpu.exe not found under $build"; exit 1 }

# Runs the benchmark (streaming its output), drops debug log lines, optionally keeps only
# the lines matching $Filter.
function Invoke-Bench([string[]]$Arguments, [string]$Filter = "") {
    & $bin.FullName @Arguments 2>&1 |
        ForEach-Object { "$_" } |
        Where-Object { $_ -notmatch "\[debug\]" -and ($Filter -eq "" -or $_ -match $Filter) } |
        ForEach-Object { $_ -replace "^.*profile:", "  GEMM profile:" } |
        Write-Log
}

if ($ListDevices) {
    & $bin.FullName --list-devices 2>&1 | ForEach-Object { "$_" } | Where-Object { $_ -notmatch "\[debug\]" }
    exit 0
}

Write-Host "Running (results go to $log) ..."
# Extra lines next to the system description that md_bench_gto_gpu prints itself.
Set-Content -Path $log -Value $null
"# bench_gto_gpu.ps1, $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz'), PowerShell $($PSVersionTable.PSVersion)" | Write-Log
try {
    Get-CimInstance Win32_VideoController | ForEach-Object {
        "# display adapter: $($_.Name), driver $($_.DriverVersion), $($_.DriverDate)"
    } | Write-Log
} catch { }
if (Get-Command nvidia-smi -ErrorAction SilentlyContinue) {
    & nvidia-smi --query-gpu=name,driver_version,memory.total,clocks.max.sm --format=csv,noheader 2>$null |
        ForEach-Object { "# nvidia-smi: $_" } | Write-Log
}

Invoke-Bench $BenchArgs
if ($LASTEXITCODE -ne 0) {
    "md_bench_gto_gpu failed (exit code $LASTEXITCODE)" | Write-Log
    exit 1
}

$profileFilter = "=== case|-- grid|^  \[|profile:"
"" | Write-Log
"# GEMM path phase breakdown (MD_GTO_GEMM_PROFILE=1, synchronises between passes)" | Write-Log
$env:MD_GTO_GEMM_PROFILE = "1"
Invoke-Bench (@("--iters", "2", "--seconds", "0", "--case", "mol", "--case", "c60f", "--case", "c240") + $BenchArgs + @("--algos", "gemm-v1,gemm")) $profileFilter
"" | Write-Log
"# GEMM path phase breakdown, orbitals (mo-gemm)" | Write-Log
Invoke-Bench (@("--iters", "2", "--seconds", "0", "--case", "mol", "--case", "c240") + $BenchArgs + @("--algos", "mo-gemm")) $profileFilter
Remove-Item Env:\MD_GTO_GEMM_PROFILE

Write-Host ""
Write-Host "Results written to $log"

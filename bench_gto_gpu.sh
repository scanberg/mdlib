#!/usr/bin/env bash
# Builds mdlib (GPU backend, Release) in build-bench/ and compares the GPU GTO kernels:
#   * electron density: the reference kernel, the tiled kernel and the two-pass GEMM
#     path in its configurations,
#   * molecular orbitals (1 orbital psi and psi^2, 32 orbitals sum of psi^2): the
#     reference kernel, the shell kernel in its configurations and the GEMM path.
# Every result is checked against the reference kernel and a CPU evaluation.
# Linux (Vulkan) and macOS (Metal). Windows: bench_gto_gpu.ps1.
#
#   ./bench_gto_gpu.sh            full run (5-30 minutes depending on the GPU)
#   ./bench_gto_gpu.sh --quick    smoke test (64^3 grids, 2 iterations)
#
# Requirements: CMake, a C compiler, Python 3, git and network access on the first
# configure (Slang and, for Vulkan, the Vulkan headers and volk are downloaded), and a
# GPU driver with Vulkan (Linux) or Metal (macOS). HDF5 is optional.
#
# Extra arguments are passed to md_bench_gto_gpu (--case NAME, --algos LIST, --dim N,
# --iters N, --seconds S, --scratch-mb MB); see benchmark/bench_gto_gpu.c.
#
# Choosing the GPU on machines with several:
#   ./bench_gto_gpu.sh --list-devices          lists the adapters and exits
#   ./bench_gto_gpu.sh --device intel          part of the adapter name, or its number
#   ./bench_gto_gpu.sh --prefer low-power      integrated before discrete
# Without these, MD_GPU_DEVICE is used when set, else the discrete GPU. A --device that
# matches nothing is an error (no silent fallback to another GPU).
#
# Results go to bench_gto_gpu_<host>.txt, or bench_gto_gpu_<host>_<device>.txt when a GPU
# is selected, so runs on different GPUs keep separate logs. The log starts with a
# description of the system
# (OS, CPU, memory, compiler, mdlib revision, GPU and driver).
set -euo pipefail
cd "$(dirname "$0")"

BUILD=build-bench
HOST="$(hostname -s 2>/dev/null || hostname)"
# Log name: one per host and GPU selection.
SEL="${MD_GPU_DEVICE:-}"; PREF=""; LIST=0
ARGS=("$@")
for ((i = 0; i < ${#ARGS[@]}; ++i)); do
    case "${ARGS[$i]}" in
        --device) SEL="${ARGS[$((i + 1))]:-}" ;;
        --prefer) PREF="${ARGS[$((i + 1))]:-}" ;;
        --list-devices) LIST=1 ;;
    esac
done
TAG=""
[ -n "$SEL" ]  && TAG+="_$SEL"
[ -n "$PREF" ] && TAG+="_$PREF"
TAG="$(printf '%s' "$TAG" | tr -c 'A-Za-z0-9_.-' '-')"
LOG="bench_gto_gpu_${HOST}${TAG}.txt"

CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release -DMD_ENABLE_GPU=ON -DMD_UNITTEST=OFF -DMD_BENCHMARK=ON)
if [ -d /opt/homebrew/lib/cmake/hdf5 ]; then
    CMAKE_ARGS+=(-DHDF5_DIR=/opt/homebrew/lib/cmake/hdf5)
fi
# Reuse an already downloaded slang instead of fetching it again.
if [ -d build/third_party/slang ]; then
    CMAKE_ARGS+=(-DSLANG_CACHE_DIR="$PWD/build/third_party/slang")
fi

echo "Configuring $BUILD ..."
# With HDF5 if it can be found (the 'mol' case then uses a real SCF density), else without.
if ! cmake -S . -B "$BUILD" "${CMAKE_ARGS[@]}" -DMD_ENABLE_HDF5=ON > "$BUILD.configure.log" 2>&1; then
    echo "  (HDF5 not usable, configuring without it)"
    rm -f "$BUILD/CMakeCache.txt"
    if ! cmake -S . -B "$BUILD" "${CMAKE_ARGS[@]}" -DMD_ENABLE_HDF5=OFF > "$BUILD.configure.log" 2>&1; then
        tail -40 "$BUILD.configure.log" | tee "$LOG"
        echo "Configure failed, see $BUILD.configure.log" | tee -a "$LOG"
        exit 1
    fi
fi
JOBS="$( (nproc || sysctl -n hw.ncpu) 2>/dev/null || echo 4)"
echo "Building md_bench_gto_gpu ..."
if ! cmake --build "$BUILD" --config Release --target md_bench_gto_gpu -j "$JOBS" > "$BUILD.build.log" 2>&1; then
    grep -E "error|Error" "$BUILD.build.log" | head -40 | tee "$LOG"
    echo "Build failed, see $BUILD.build.log" | tee -a "$LOG"
    exit 1
fi

BIN="$BUILD/bin/md_bench_gto_gpu"
[ -x "$BIN" ] || BIN="$(find "$BUILD" -name md_bench_gto_gpu -type f -perm -u+x | head -1)"

if [ "$LIST" = 1 ]; then
    "$BIN" --list-devices 2>&1 | grep -v "\[debug\]"
    exit 0
fi

echo "Running (results go to $LOG) ..."
{
    echo "# bench_gto_gpu.sh, $(date)"
    if command -v nvidia-smi >/dev/null; then
        nvidia-smi --query-gpu=name,driver_version,memory.total,clocks.max.sm --format=csv,noheader 2>/dev/null \
            | sed 's/^/# nvidia-smi: /' || true
    fi
    if command -v system_profiler >/dev/null; then
        system_profiler SPDisplaysDataType 2>/dev/null | grep -E "Chipset Model|Total Number of Cores|Metal" \
            | sed -E 's/^ +/# display: /' || true
    fi
    if ! "$BIN" "$@"; then
        echo "md_bench_gto_gpu failed"
        exit 1
    fi
    echo
    echo "# GEMM path phase breakdown (MD_GTO_GEMM_PROFILE=1, synchronises between passes)"
    MD_GTO_GEMM_PROFILE=1 "$BIN" --iters 2 --seconds 0 --case mol --case c60f --case c240 "$@" --algos gemm-v1,gemm \
        2>&1 | grep -E "=== case|-- grid|^  \[|profile:" | sed -E 's/^.*profile:/  GEMM profile:/' || true
    echo
    echo "# GEMM path phase breakdown, orbitals (mo-gemm)"
    MD_GTO_GEMM_PROFILE=1 "$BIN" --iters 2 --seconds 0 --case mol --case c240 "$@" --algos mo-gemm \
        2>&1 | grep -E "=== case|-- grid|^  \[|profile:" | sed -E 's/^.*profile:/  GEMM profile:/' || true
} 2>&1 | grep -v "\[debug\]" | tee "$LOG"

echo
echo "Results written to $(pwd)/$LOG"

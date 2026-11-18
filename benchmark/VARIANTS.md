# Kernel variants: the measurements behind the pruning (temporary)

What each variant is measured against, per GPU, before it goes (`md_gto_gemm_bench --variants`,
`GTO_ISO_VARIANTS=1 md_unittest --filter=gto_iso_gpu_lib.variants`). The rule: a variant goes where what replaces it
is within 5% on all four GPUs, on the workloads that would take it - so one GPU where it is more than 5% ahead keeps
it, and the others need not measure it. This file goes with the switches once decided.

Ratios are B / A, the median of 7 rounds in alternating order; below 1, B is faster. Two runs on the M4 Pro
(gemm-bench r5 and r6, isosurface r41 and r42): both shown where both measured it. The GTX 1060 (Windows, Vulkan,
driver 582.66; gemm-bench r6): a few of its lines were lost to the console (-). The HD 530 (the same machine, Vulkan
1.3.215; gemm-bench r6).

## Grid kernels (md_gto)

### Block shape: a kernel per shape (A) against one kernel, the shape an argument (B)

| pass             | GPU     | mol vol     | mol plane   | protein vol     | protein plane |
|------------------|---------|-------------|-------------|-----------------|---------------|
| density GEMM     | M4 Pro  | 1.006 1.000 | 0.961 1.031 | 0.995 0.999     | 1.032 1.017   |
| ELF GEMM         | M4 Pro  | 0.994 0.987 | 1.012 0.951 | 1.001 1.001     | 1.006 0.990   |
| orbitals shell 1 | M4 Pro  | 1.063 1.071 | 1.011 1.016 | **1.051 1.051** | 0.996 1.007   |
| orbitals shell 8 | M4 Pro  | 0.961 1.012 | 0.992 0.984 | 0.998 1.000     | 1.002 1.091   |
| orbitals GEMM 32 | M4 Pro  | 1.015 0.989 | 1.029 1.007 | 0.998 1.003     | 0.974 1.129   |
| density GEMM     | 1060    | 1.009       | 1.027       | 1.000           | 0.999         |
| ELF GEMM         | 1060    | 0.990       | 1.006       | 0.997           | 1.000         |
| orbitals shell 1 | 1060    | 0.995       | -           | 1.002           | -             |
| orbitals shell 8 | 1060    | 1.088       | 0.994       | **1.176**       | **1.189**     |
| orbitals GEMM 32 | 1060    | 1.022       | 1.065       | 0.992           | 1.009         |
| density GEMM     | HD 530  | 1.000       | 0.997       | 1.000           | 0.996         |
| ELF GEMM         | HD 530  | 1.000       | 1.000       | 1.000           | 1.002         |
| orbitals shell 1 | HD 530  | 0.998       | 1.002       | 0.970           | 0.967         |
| orbitals shell 8 | HD 530  | 0.996       | 0.997       | 0.956           | 0.964         |
| orbitals GEMM 32 | HD 530  | 1.001       | 1.031       | 1.000           | 0.996         |

(The planes' runs are 0.2-1.6 ms: their spread is noise - the protein plane's 1.091 and 1.129 had ranges of
0.93-1.25 and 0.98-1.18.) Every pair the same bits, but the shell kernel on the HD 530 (2e-7 to 7e-7 of the largest
value: its compiler contracts the runtime shape's arithmetic otherwise).

### AO tile: 32 (A) against 64 with 32 for blocks of <= 32 AOs (B); the vendor's choice today: Apple 32, else 64

| pass    | GPU    | mol vol     | mol plane   | protein vol | protein plane |
|---------|--------|-------------|-------------|-------------|---------------|
| density | M4 Pro | 1.112 1.112 | 0.950 1.000 | 1.043 1.046 | 1.072 1.042   |
| ELF     | M4 Pro | 0.985 0.998 | 1.043 1.105 | 0.970 0.963 | 0.958 0.957   |
| density | 1060   | 1.209       | 1.063       | 1.137       | 1.134         |
| ELF     | 1060   | 0.983       | 0.999       | **0.935**   | **0.936**     |
| density | HD 530 | **0.854**   | **0.817**   | **0.766**   | **0.733**     |
| ELF     | HD 530 | 0.951       | **0.898**   | **0.823**   | **0.784**     |

### Orbitals: the shell kernel, 8 a pass (A), against the GEMM path (B), volume

| orbitals | M4 Pro mol  | M4 Pro protein | 1060 mol | 1060 protein | HD 530 mol | HD 530 protein |
|----------|-------------|----------------|----------|--------------|------------|----------------|
| 2        | 3.238 3.375 | 4.805 4.752    | 4.980    | 2.030        | 5.166      | 2.382          |
| 4        | 2.707 2.750 | 3.691 3.713    | 3.975    | 1.556        | 3.640      | 1.941          |
| 8        | 2.400 2.421 | 2.780 2.750    | 3.643    | 1.091        | 2.768      | 1.457          |
| 16       | 1.304       | 1.396          | 1.894    | 0.545        | 1.390      | 0.738          |
| 24       | 0.903       | 0.949          | 1.234    | 0.365        | 0.938      | 0.507          |
| 32       | 0.671       | 0.728          | 0.912    | -            | 0.714      | 0.392          |

(The GEMM path's time is flat from 2 to 32 orbitals: on the protein 124 ms on the M4 Pro, 194 ms on the 1060, 2.4 s
on the HD 530; the shell kernel per pass of 8: 45, 178 ms and 1.66 s. The crossover - GEMM over the shell kernel's
passes - is not the vendor's alone: on mol at 19-29 orbitals on all three, on the protein at 22 on the M4 Pro, 9 on
the 1060, 12 on the HD 530. Past 8 today: at most 40% lost at 9-16 orbitals on mol, against up to 84% that 16 would
lose on the 1060's protein. Left as it is.)

## Isosurface GEMMs (md_gto_iso_gpu)

Densities of the N highest occupied orbitals at 0.002: the extraction's total (the GEMMs' stages: corners, project,
relax). mol.h5 at cell 0.15, protein.h5 at 0.2. A wide GEMM of fewer columns than N runs in passes.

| case                          | GPU    | B/A total   | B/A gemm    |
|-------------------------------|--------|-------------|-------------|
| tiled N 8, 16 vs 32           | M4 Pro | 1.049       | 1.065       |
| tiled N 16, 16 vs 32          | M4 Pro | 1.048       | 1.078       |
| tiled N 24, 32 vs 64          | M4 Pro | 1.349       | 1.502       |
| tiled N 32, 32 vs 64          | M4 Pro | 1.370       | 1.516       |
| N 96, tiled 64 vs w128        | M4 Pro | 0.796 0.797 | 0.756 0.757 |
| N 96, w128 vs y256            | M4 Pro | 1.570 1.574 | 1.717 1.720 |
| N 128, w128 vs y256           | M4 Pro | 1.577 1.577 | 1.711 1.707 |
| N 160, w192 vs y256           | M4 Pro | 1.243 1.244 | 1.280 1.280 |
| N 192, w192 vs y256           | M4 Pro | 1.238 1.241 | 1.273 1.273 |
| N 160, w128 (2 passes) vs w192 | M4 Pro | 0.716      | 0.687       |
| N 192, w128 (2 passes) vs w192 | M4 Pro | 0.718      | 0.691       |
| N 256, w128 (2) vs y256       | M4 Pro | 0.897       | 0.887       |
| N 256, w192 (2) vs y256       | M4 Pro | 0.694       | 0.672       |
| N 512, w128 (4) vs y256 (2)   | M4 Pro | 0.884       | 0.878       |
| N 512, w192 (3) vs y256 (2)   | M4 Pro | 0.857       | 0.851       |
| N 1024, w128 (8) vs y256 (4)  | M4 Pro | 0.889       | 0.886       |
| N 1024, w192 (6) vs y256 (4)  | M4 Pro | 0.855       | 0.852       |
| N 96, tiled 64 vs w128        | 1060   | 0.986       | 0.979       |
| N 96, w128 vs y256            | 1060   | 1.495       | 1.644       |
| N 128, w128 vs y256           | 1060   | 1.502       | 1.643       |
| N 160, w192 vs y256           | 1060   | 1.159       | 1.186       |
| N 192, w192 vs y256           | 1060   | 1.160       | 1.187       |
| N 160, w128 (2 passes) vs w192 | 1060  | 0.750       | 0.719       |
| N 192, w128 (2 passes) vs w192 | 1060  | 0.739       | 0.709       |
| N 256, w128 (2) vs y256       | 1060   | 0.855       | 0.842       |
| N 256, w192 (2) vs y256       | 1060   | 0.644       | 0.615       |

(The 1060 lost the device at N 512 with w128 in four passes: md_gto_iso_gpu_extract made the corners one launch,
over Windows' 2 s watchdog. Isosurface r43 chunks by the GEMMs' columns, the blocking extraction too.)

## Decided (one GPU more than 5% ahead keeps a variant)

- Wide 128, 192 and 256 stay, each the fastest in its band on the M4 Pro: at 160-192 w192 is 28% ahead of w128 in
  two passes and 24% ahead of y256 (whose idle quarter of columns that is); past 192 y256 is 10-15% ahead of passes
  of w128 or w192.
- Tiled 32 stays (64: 35-37% behind at 24-32 functions); tiled 64 is the only one for sets of more than 32 functions
  that are not canonical.
- The 8-orbital shell kernel stays (the GEMM path 2.4-4.8x behind at 2-8 orbitals).
- The shell kernel keeps a kernel per shape, at 1 orbital (one kernel 5.1% behind on the M4 Pro's volume, twice) and
  at 8 (9-19% behind on the 1060).
- The density's GEMM keeps both tiles: 32 (64 4-11% behind on the M4 Pro, 6-21% on the 1060) and 64 (15-27% ahead
  of 32 on the HD 530). NVIDIA's 64 today is the wrong choice for the density on the 1060.
- ELF/LOL's GEMM keeps both: 64 ahead of 32 everywhere (3-4% on the M4 Pro, 6.5% on the 1060's protein, 5-22% on
  the HD 530), 32 for the small blocks.
- So the tile becomes a choice per pass and vendor, not per vendor: ELF/LOL 64; the density 32 on Apple and (the
  1060; the 3080 to see) NVIDIA, 64 on Intel.

## Open (the 1060, the Intel iGPU, the 3080)

- One kernel for both shapes: the counting and Phi passes, the density's GEMM pass, ELF/LOL's finish and Phi with
  gradients, the orbitals' GEMM pass (neutral on the M4 Pro, the 1060 and the HD 530): the 3080 to see.
- Tiled 16 against 32: 4.8-4.9% on the M4 Pro's totals (6.5-7.8% on the GEMMs' stages); the 1060's, the HD 530's
  and the 3080's to run (GTO_ISO_VARIANTS_ONLY=tiled).

# Kernel variants: the measurements behind the pruning (temporary)

What each variant is measured against, per GPU, before it goes (`md_gto_gemm_bench --variants`,
`GTO_ISO_VARIANTS=1 md_unittest --filter=gto_iso_gpu_lib.variants`). The rule: a variant goes where what replaces it
is within 5% on all four GPUs, on the workloads that would take it. This file goes with the switches once decided.

Ratios are B / A, the median of 7 rounds in alternating order; below 1, B is faster.

## Grid kernels (md_gto)

### Block shape: a kernel per shape (A) against one kernel, the shape an argument (B)

| pass             | GPU     | mol vol | mol plane | protein vol | protein plane |
|------------------|---------|---------|-----------|-------------|---------------|
| density GEMM     | M4 Pro  | 1.006   | 0.961     | 0.995       | 1.032         |
| ELF GEMM         | M4 Pro  | 0.994   | 1.012     | 1.001       | 1.006         |
| orbitals shell 1 | M4 Pro  | 1.063   | 1.011     | **1.051**   | 0.996         |
| orbitals shell 8 | M4 Pro  | 0.961   | 0.992     | 0.998       | 1.002         |
| orbitals GEMM 32 | M4 Pro  | 1.015   | 1.029     | 0.998       | 0.974         |

(mol's plane runs are 0.2-1 ms: noise.) Every pair the same bits.

### AO tile: 32 (A) against 64 with 32 for blocks of <= 32 AOs (B); the vendor's choice today: Apple 32, else 64

| pass    | GPU    | mol vol | mol plane | protein vol | protein plane |
|---------|--------|---------|-----------|-------------|---------------|
| density | M4 Pro | 1.112   | 0.950     | 1.043       | 1.072         |
| ELF     | M4 Pro | 0.985   | 1.043     | 0.970       | 0.958         |

### Orbitals: the shell kernel, 8 a pass (A), against the GEMM path (B), volume

| orbitals | GPU    | mol   | protein |
|----------|--------|-------|---------|
| 2        | M4 Pro | 3.238 | 4.805   |
| 4        | M4 Pro | 2.707 | 3.691   |
| 8        | M4 Pro | 2.400 | 2.780   |

(protein: the GEMM path 124 ms for 2 to 32 orbitals; the shell kernel 45 ms for 8.)

## Isosurface GEMMs (md_gto_iso_gpu)

Densities of the N highest occupied orbitals at 0.002: the extraction's total (the GEMMs' stages: corners, project,
relax). mol.h5 at cell 0.15, protein.h5 at 0.2.

| case                     | GPU    | B/A total | B/A gemm |
|--------------------------|--------|-----------|----------|
| tiled N 8, 16 vs 32      | M4 Pro | 1.049     | 1.065    |
| tiled N 16, 16 vs 32     | M4 Pro | 1.048     | 1.078    |
| tiled N 24, 32 vs 64     | M4 Pro | 1.349     | 1.502    |
| tiled N 32, 32 vs 64     | M4 Pro | 1.370     | 1.516    |
| wide N 96, tiled 64 vs w128 | M4 Pro | 0.796  | 0.756    |
| wide N 96, w128 vs y256  | M4 Pro | 1.570     | 1.717    |
| wide N 128, w128 vs y256 | M4 Pro | 1.577     | 1.711    |
| wide N 160, w192 vs y256 | M4 Pro | 1.243     | 1.280    |
| wide N 192, w192 vs y256 | M4 Pro | 1.238     | 1.273    |

## Standing (M4 Pro only)

- One kernel for both block shapes: neutral for the GEMM passes and the 8-orbital shell kernel; the 1-orbital shell
  kernel on a volume 5-6% slower: it keeps a kernel per shape unless that changes.
- AO tile: Apple's 32 confirmed for the density (64: 4-11% slower); ELF 3-4% faster at 64.
- The 8-orbital shell kernel stays (the GEMM path 2.4-4.8x slower at 2-8 orbitals); the crossover is past 8 here.
- Tiled 16: within 5% on the totals (5-8% on the GEMMs' stages) - open. Tiled 32 stays (64: 35-37% slower).
- Wide 128 stays (tiled 64: 26% slower at 96; y256 57%). w192 against y256: y256 24% slower even full - y256 looks
  the weak one: open whether passes of w128 or w192 replace it (193 to 1024 functions).

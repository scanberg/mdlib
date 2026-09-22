#pragma once

#include <stddef.h>
#include <stdbool.h>

// Thin wrapper around PFFFT (ext/pffft) for the transforms mdlib needs.
//
// Conventions
// - Forward transforms use the kernel exp(-2*pi*i*j*k/N) and are NOT normalized.
// - Complex data is stored interleaved (re, im).
// - All buffers passed to the transform functions must be 64 byte aligned. Use md_fft_alloc / md_fft_free
//   or any other aligned allocation.
// - Plans are immutable after creation and can be shared between threads, scratch buffers can not.
//
// Size restrictions
// PFFFT supports lengths N = N_min * 2^a * 3^b * 5^c, where N_min depends on the SIMD width and transform type.
// Use md_fft_valid_size to find the closest supported length.

#ifdef __cplusplus
extern "C" {
#endif

// Returns the smallest supported transform length >= n.
// real: true for real to complex transforms (stricter requirement), false for complex to complex.
int md_fft_valid_size(int n, bool real);

// 64 byte aligned heap allocations suitable for the transform buffers.
void* md_fft_alloc(size_t bytes);
void  md_fft_free(void* ptr);

typedef struct md_fft_2d_t md_fft_2d_t;

// 2D real to complex transform of a nx * ny real field, stored row major with x as the fastest dimension.
// nx must be a valid real size and ny a valid complex size (both are satisfied by md_fft_valid_size(n, true)).
md_fft_2d_t* md_fft_2d_create(int nx, int ny);
void         md_fft_2d_destroy(md_fft_2d_t* fft);

// Number of floats required for the scratch buffer passed to md_fft_2d_r2c
size_t md_fft_2d_scratch_size(const md_fft_2d_t* fft);

// Forward 2D real to complex transform
// in:      nx * ny real values (row major, x fastest). Not modified.
// out:     (nx/2 + 1) * ny complex values (interleaved), row major with kx fastest:
//          out[2 * (ky * (nx/2+1) + kx) + 0] = Re, [.. + 1] = Im
//          ky follows standard FFT ordering (ky > ny/2 represent negative frequencies ky - ny).
// kx_count: Only columns kx < kx_count are transformed along y (the rest are left undefined).
//          Pass 0 or (nx/2+1) to compute all columns. Useful when only a low frequency disk is required.
// scratch: md_fft_2d_scratch_size(fft) floats, 64 byte aligned.
void md_fft_2d_r2c(const md_fft_2d_t* fft, const float* in, float* out, int kx_count, float* scratch);

#ifdef __cplusplus
}
#endif

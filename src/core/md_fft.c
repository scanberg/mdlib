#include <core/md_fft.h>
#include <core/md_common.h>

#include <pffft.h>

#include <string.h>

// Number of columns gathered and transformed together in the column pass, to improve cache utilization
#define COL_BLOCK 8

struct md_fft_2d_t {
    int nx, ny;
    PFFFT_Setup* row;   // Real, length nx
    PFFFT_Setup* col;   // Complex, length ny
};

int md_fft_valid_size(int n, bool real) {
    if (n < 1) n = 1;
    // A size which is valid for real transforms of length n is also required to be valid for complex
    // transforms when used as the second dimension in 2D, so we search for a size which satisfies both.
    int N = pffft_nearest_transform_size(n, real ? PFFFT_REAL : PFFFT_COMPLEX, 1);
    if (real) {
        while (!pffft_is_valid_size(N, PFFFT_COMPLEX) || !pffft_is_valid_size(N, PFFFT_REAL)) {
            N = pffft_nearest_transform_size(N + 1, PFFFT_REAL, 1);
        }
    }
    return N;
}

void* md_fft_alloc(size_t bytes) {
    return pffft_aligned_malloc(bytes);
}

void md_fft_free(void* ptr) {
    pffft_aligned_free(ptr);
}

md_fft_2d_t* md_fft_2d_create(int nx, int ny) {
    if (nx <= 0 || ny <= 0) return NULL;
    if (!pffft_is_valid_size(nx, PFFFT_REAL) || !pffft_is_valid_size(ny, PFFFT_COMPLEX)) {
        return NULL;
    }
    md_fft_2d_t* fft = (md_fft_2d_t*)pffft_aligned_malloc(sizeof(md_fft_2d_t));
    if (!fft) return NULL;
    fft->nx = nx;
    fft->ny = ny;
    fft->row = pffft_new_setup(nx, PFFFT_REAL);
    fft->col = pffft_new_setup(ny, PFFFT_COMPLEX);
    if (!fft->row || !fft->col) {
        md_fft_2d_destroy(fft);
        return NULL;
    }
    return fft;
}

void md_fft_2d_destroy(md_fft_2d_t* fft) {
    if (!fft) return;
    if (fft->row) pffft_destroy_setup(fft->row);
    if (fft->col) pffft_destroy_setup(fft->col);
    pffft_aligned_free(fft);
}

size_t md_fft_2d_scratch_size(const md_fft_2d_t* fft) {
    if (!fft) return 0;
    const size_t nx = (size_t)fft->nx;
    const size_t ny = (size_t)fft->ny;
    // row: nx (tmp) + nx (work)
    // col: COL_BLOCK * 2*ny (gathered columns) + 2*ny (tmp) + 2*ny (work)
    const size_t row = 2 * nx;
    const size_t col = (COL_BLOCK + 2) * 2 * ny;
    return MAX(row, col) + 16;
}

void md_fft_2d_r2c(const md_fft_2d_t* fft, const float* in, float* out, int kx_count, float* scratch) {
    ASSERT(fft);
    ASSERT(in);
    ASSERT(out);
    ASSERT(scratch);

    const int nx = fft->nx;
    const int ny = fft->ny;
    const int nkx = nx / 2 + 1;
    if (kx_count <= 0 || kx_count > nkx) kx_count = nkx;

    // Row pass (real -> complex along x)
    {
        float* tmp  = scratch;
        float* work = scratch + nx;
        for (int y = 0; y < ny; ++y) {
            pffft_transform_ordered(fft->row, in + (size_t)y * nx, tmp, work, PFFFT_FORWARD);
            float* dst = out + (size_t)y * nkx * 2;
            // Unpack PFFFT's packed format: tmp[0] = DC, tmp[1] = Nyquist, then (re,im) for k = 1 .. nx/2-1
            const int n = MIN(kx_count, nkx);
            dst[0] = tmp[0];
            dst[1] = 0.0f;
            const int k_end = MIN(n, nx / 2);
            if (k_end > 1) {
                memcpy(dst + 2, tmp + 2, sizeof(float) * 2 * (size_t)(k_end - 1));
            }
            if (n == nkx) {
                dst[2 * (nx / 2) + 0] = tmp[1];
                dst[2 * (nx / 2) + 1] = 0.0f;
            }
        }
    }

    // Column pass (complex -> complex along y), in blocks of columns
    {
        float* block = scratch;                          // COL_BLOCK * 2 * ny
        float* tmp   = scratch + (size_t)COL_BLOCK * 2 * ny;  // 2 * ny
        float* work  = tmp + (size_t)2 * ny;               // 2 * ny
        const size_t stride = (size_t)nkx * 2;

        for (int kx0 = 0; kx0 < kx_count; kx0 += COL_BLOCK) {
            const int nb = MIN(COL_BLOCK, kx_count - kx0);
            // Gather
            for (int y = 0; y < ny; ++y) {
                const float* src = out + y * stride + (size_t)kx0 * 2;
                for (int b = 0; b < nb; ++b) {
                    block[(size_t)b * 2 * ny + 2 * y + 0] = src[2 * b + 0];
                    block[(size_t)b * 2 * ny + 2 * y + 1] = src[2 * b + 1];
                }
            }
            // Transform
            for (int b = 0; b < nb; ++b) {
                float* col = block + (size_t)b * 2 * ny;
                pffft_transform_ordered(fft->col, col, tmp, work, PFFFT_FORWARD);
                memcpy(col, tmp, sizeof(float) * 2 * (size_t)ny);
            }
            // Scatter
            for (int y = 0; y < ny; ++y) {
                float* dst = out + y * stride + (size_t)kx0 * 2;
                for (int b = 0; b < nb; ++b) {
                    dst[2 * b + 0] = block[(size_t)b * 2 * ny + 2 * y + 0];
                    dst[2 * b + 1] = block[(size_t)b * 2 * ny + 2 * y + 1];
                }
            }
        }
    }
}

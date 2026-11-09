#include "utest.h"

#include <core/md_linalg.h>
#include <core/md_executor.h>
#include <core/md_os.h>

#include <stdio.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// A reproducible symmetric matrix with a spread spectrum: an LCG, so no platform rand() differences.
static void sym_random(double* a, size_t n, unsigned seed) {
    unsigned s = seed;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            s = s * 1664525u + 1013904223u;
            const double v = (double)(s >> 8) / (double)(1u << 24) - 0.5;
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
        a[i * n + i] += (double)i * 0.1;   // spread the spectrum
    }
}

// max |A v_j - w_j v_j| and max |V^T V - I|
static void residuals(double* out_eig, double* out_orth, const double* A, const double* V, const double* w, size_t n) {
    double re = 0.0, ro = 0.0;
    for (size_t j = 0; j < n; ++j) {
        const double* v = V + j * n;
        for (size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (size_t k = 0; k < n; ++k) s += A[i * n + k] * v[k];
            re = fmax(re, fabs(s - w[j] * v[i]));
        }
        for (size_t k = 0; k < n; ++k) {
            double d = 0.0;
            for (size_t i = 0; i < n; ++i) d += v[i] * V[k * n + i];
            ro = fmax(ro, fabs(d - (j == k ? 1.0 : 0.0)));
        }
    }
    *out_eig = re;
    *out_orth = ro;
}

UTEST(linalg, sym_eigen_random) {
    const size_t sizes[] = { 1, 2, 3, 7, 16, 61, 128 };
    for (size_t t = 0; t < sizeof(sizes) / sizeof(sizes[0]); ++t) {
        const size_t n = sizes[t];
        double* A = malloc(sizeof(double) * n * n);
        double* V = malloc(sizeof(double) * n * n);
        double* w = malloc(sizeof(double) * n);
        sym_random(A, n, 17u + (unsigned)n);
        memcpy(V, A, sizeof(double) * n * n);

        ASSERT_TRUE(md_linalg_sym_eigen(V, w, n));
        for (size_t j = 1; j < n; ++j) {
            EXPECT_LE(w[j - 1], w[j]);
        }
        double re, ro;
        residuals(&re, &ro, A, V, w, n);
        EXPECT_LT(re, 1.0e-10 * (double)n);
        EXPECT_LT(ro, 1.0e-10 * (double)n);

        free(A); free(V); free(w);
    }
}

// Degenerate and already diagonal inputs: the QL iteration must not stall on them.
UTEST(linalg, sym_eigen_degenerate) {
    enum { N = 9 };
    double A[N * N] = {0};
    double V[N * N];
    double w[N];
    for (int i = 0; i < N; ++i) A[i * N + i] = (i % 3 == 0) ? 2.0 : -1.0;
    memcpy(V, A, sizeof(V));
    ASSERT_TRUE(md_linalg_sym_eigen(V, w, N));
    for (int i = 0; i < 6; ++i) EXPECT_NEAR(w[i], -1.0, 1.0e-14);
    for (int i = 6; i < N; ++i) EXPECT_NEAR(w[i],  2.0, 1.0e-14);
    double re, ro;
    residuals(&re, &ro, A, V, w, N);
    EXPECT_LT(re, 1.0e-13);
    EXPECT_LT(ro, 1.0e-13);

    // zero matrix
    memset(V, 0, sizeof(V));
    ASSERT_TRUE(md_linalg_sym_eigen(V, w, N));
    for (int i = 0; i < N; ++i) EXPECT_EQ(w[i], 0.0);
}

UTEST(linalg, sym_eigen_rejects) {
    double w[1];
    EXPECT_FALSE(md_linalg_sym_eigen(NULL, w, 1));
    double a[1] = { 3.0 };
    EXPECT_FALSE(md_linalg_sym_eigen(a, w, 0));
    ASSERT_TRUE(md_linalg_sym_eigen(a, w, 1));
    EXPECT_EQ(w[0], 3.0);
    EXPECT_EQ(a[0], 1.0);
}

// A matrix with a prescribed spectrum: diag(lam) turned by three Householder reflections
static void sym_with_spectrum(double* a, size_t n, const double* lam, unsigned seed) {
    memset(a, 0, sizeof(double) * n * n);
    for (size_t i = 0; i < n; ++i) a[i * n + i] = lam[i];
    double* v = malloc(sizeof(double) * n);
    double* y = malloc(sizeof(double) * n);
    unsigned s = seed;
    for (int h = 0; h < 3; ++h) {
        double nn = 0.0;
        for (size_t i = 0; i < n; ++i) {
            s = s * 1664525u + 1013904223u;
            v[i] = (double)(s >> 8) / (double)(1u << 24) - 0.5;
            nn += v[i] * v[i];
        }
        for (size_t i = 0; i < n; ++i) v[i] /= sqrt(nn);
        for (size_t i = 0; i < n; ++i) {
            double t = 0.0;
            for (size_t k = 0; k < n; ++k) t += a[i * n + k] * v[k];
            y[i] = t;
        }
        double vy = 0.0;
        for (size_t i = 0; i < n; ++i) vy += v[i] * y[i];
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) a[i * n + j] += -2.0 * v[i] * y[j] - 2.0 * y[i] * v[j] + 4.0 * vy * v[i] * v[j];
        }
    }
    free(v);
    free(y);
}

// The spectra that divide and conquer finds hard: many equal eigenvalues (deflation by rotation), a
// spread over 24 orders of magnitude, a cluster within 1e-13, rank one. Large enough to recurse.
UTEST(linalg, sym_eigen_hard_spectra) {
    const size_t n = 157;
    double* lam = malloc(sizeof(double) * n);
    double* A = malloc(sizeof(double) * n * n);
    double* V = malloc(sizeof(double) * n * n);
    double* w = malloc(sizeof(double) * n);
    for (int kind = 0; kind < 5; ++kind) {
        double scale = 1.0;
        for (size_t i = 0; i < n; ++i) {
            switch (kind) {
            case 0: lam[i] = (double)(i % 4); break;
            case 1: lam[i] = pow(10.0, -12.0 + 24.0 * (double)i / (double)(n - 1)); scale = 1.0e12; break;
            case 2: lam[i] = 1.0 + 1.0e-13 * sin((double)i); break;
            case 3: lam[i] = i == n - 1 ? (double)n : 0.0; scale = (double)n; break;
            case 4: lam[i] = 0.0; break;
            }
        }
        sym_with_spectrum(A, n, lam, 99u + (unsigned)kind);
        memcpy(V, A, sizeof(double) * n * n);
        ASSERT_TRUE(md_linalg_sym_eigen(V, w, n));
        for (size_t j = 1; j < n; ++j) EXPECT_LE(w[j - 1], w[j]);
        double re, ro;
        residuals(&re, &ro, A, V, w, n);
        EXPECT_LT(re, 1.0e-13 * scale);
        EXPECT_LT(ro, 1.0e-13);
        if (kind == 0) {
            // 157 = 40 + 39 + 39 + 39: 0 x40, 1 x39, 2 x39, 3 x39
            for (size_t j = 0; j < n; ++j) EXPECT_NEAR(w[j], j < 40 ? 0.0 : j < 79 ? 1.0 : j < 118 ? 2.0 : 3.0, 1.0e-13);
        }
    }
    free(lam); free(A); free(V); free(w);
}

// An executor that runs the pieces in reverse, each on a different thread_idx out of four: a result
// that depended on the order or on the thread would differ from the one without an executor.
static bool reverse_parallel_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    (void)user;
    uint32_t t = 0;
    uint32_t end = count;
    while (end > 0) {
        const uint32_t beg = end > grain ? end - grain : 0;
        body(arg, beg, end, t);
        t = (t + 1) % 4;
        end = beg;
    }
    return true;
}

UTEST(linalg, sym_eigen_executor_is_bitwise) {
    const size_t sizes[] = { 30, 260, 700 };
    md_executor_t ex = { reverse_parallel_for, NULL, 4 };
    for (size_t t = 0; t < sizeof(sizes) / sizeof(sizes[0]); ++t) {
        const size_t n = sizes[t];
        double* A = malloc(sizeof(double) * n * n);
        double* B = malloc(sizeof(double) * n * n);
        double* wa = malloc(sizeof(double) * n);
        double* wb = malloc(sizeof(double) * n);
        sym_random(A, n, 5u + (unsigned)n);
        memcpy(B, A, sizeof(double) * n * n);
        ASSERT_TRUE(md_linalg_sym_eigen_ex(A, wa, n, NULL));
        ASSERT_TRUE(md_linalg_sym_eigen_ex(B, wb, n, &ex));
        EXPECT_EQ(0, memcmp(A, B, sizeof(double) * n * n));
        EXPECT_EQ(0, memcmp(wa, wb, sizeof(double) * n));
        free(A); free(B); free(wa); free(wb);
    }
}

UTEST(linalg, sym_eigen_rejects_nonfinite) {
    const size_t n = 40;
    double* A = malloc(sizeof(double) * n * n);
    double* w = malloc(sizeof(double) * n);
    sym_random(A, n, 3u);
    A[7 * n + 11] = A[11 * n + 7] = NAN;
    EXPECT_FALSE(md_linalg_sym_eigen(A, w, n));
    sym_random(A, n, 3u);
    A[0] = INFINITY;
    EXPECT_FALSE(md_linalg_sym_eigen(A, w, n));
    free(A); free(w);
}

// Timings on this machine (LINALG_REPORT=1), one thread. The revision says which build produced them.
#define LINALG_REPORT_REVISION "linalg-dc 2026-11-09 r1"

UTEST(linalg, report) {
    const char* env = getenv("LINALG_REPORT");
    if (!env || env[0] != '1') UTEST_SKIP("set LINALG_REPORT=1 (about 15 s)");
    printf("  revision %s\n", LINALG_REPORT_REVISION);

    const size_t gn = 1000;
    double* A = malloc(sizeof(double) * gn * gn);
    double* B = malloc(sizeof(double) * gn * gn);
    double* C = malloc(sizeof(double) * gn * gn);
    sym_random(A, gn, 1u);
    sym_random(B, gn, 2u);
    double best = 1.0e30;
    for (int r = 0; r < 3; ++r) {
        const md_tick_t t0 = md_tick_now();
        ASSERT_TRUE(md_linalg_gemm(false, false, gn, gn, gn, 1.0, A, gn, B, gn, 0.0, C, gn, NULL));
        best = fmin(best, md_tick_to_seconds(md_tick_now() - t0));
    }
    printf("  md_linalg_gemm %zu^3: %.3f s, %.1f GFLOP/s\n", gn, best, 2.0 * (double)gn * gn * gn / best * 1.0e-9);
    free(A); free(B); free(C);

    const size_t sizes[] = { 250, 500, 1000, 2000 };
    for (size_t t = 0; t < sizeof(sizes) / sizeof(sizes[0]); ++t) {
        const size_t n = sizes[t];
        double* M = malloc(sizeof(double) * n * n);
        double* V = malloc(sizeof(double) * n * n);
        double* w = malloc(sizeof(double) * n);
        sym_random(M, n, 7u + (unsigned)n);
        best = 1.0e30;
        for (int r = 0; r < (n >= 2000 ? 1 : 3); ++r) {
            memcpy(V, M, sizeof(double) * n * n);
            const md_tick_t t0 = md_tick_now();
            ASSERT_TRUE(md_linalg_sym_eigen(V, w, n));
            best = fmin(best, md_tick_to_seconds(md_tick_now() - t0));
        }
        // spot check: a few eigenpairs
        double re = 0.0;
        for (size_t j = 0; j < n; j += n / 7) {
            for (size_t i = 0; i < n; ++i) {
                double s = 0.0;
                for (size_t k = 0; k < n; ++k) s += M[i * n + k] * V[j * n + k];
                re = fmax(re, fabs(s - w[j] * V[j * n + i]));
            }
        }
        printf("  md_linalg_sym_eigen n = %4zu: %.3f s   (residual %.1e)\n", n, best, re);
        EXPECT_LT(re, 1.0e-10 * (double)n);
        free(M); free(V); free(w);
    }
}

#include "utest.h"

#include <core/md_linalg.h>

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

#pragma once

#include <stdbool.h>
#include <stddef.h>

// Dense linear algebra that more than one part of mdlib needs, in ONE place.
//
// The procedures here are deliberately shaped like their LAPACK counterparts (column major, in
// place, eigenvalues ascending), so that a build which has a LAPACK can route them there without
// any caller noticing. Until then they are portable C (with md_simd for the inner kernels), fast enough
// that the LAPACK route is an option rather than a necessity.

#ifdef __cplusplus
extern "C" {
#endif

// Eigendecomposition of a real symmetric n x n matrix, in place (the LAPACK dsyev contract).
//
// a:      n * n doubles. In: the symmetric matrix (as it is symmetric, its storage order does not
//         matter). Out: the eigenvectors, eigenvector j CONTIGUOUS at a[j * n .. j * n + n) - column
//         j in column major terms - orthonormal.
// eigval: n doubles. Out: the eigenvalues in ascending order, eigval[j] that of eigenvector j.
//
// Returns false if the input held a NaN or an Inf, if the iteration did not converge (not seen in
// practice), if the executor was cancelled, or if the workspace could not be allocated; the contents
// of 'a' and 'eigval' are then undefined.
//
// LAPACK dsyevd's method: Householder tridiagonalization, divide and conquer on the tridiagonal
// problem, and the reflectors applied in blocks. Almost all of the work is matrix products and one
// pass over the matrix per reduction step. It needs 3 n^2 doubles of workspace besides 'a', and about
// 2 MB per thread.
// About 0.27 s for n = 1000 and 1.9 s for n = 2000 on one core of a 2.1 GHz Xeon (AVX2), within
// 1.5x of OpenBLAS dsyevd; the textbook method this replaced (EISPACK tred2/tql2) took 1.4 s and 12 s.
bool md_linalg_sym_eigen(double* a, double* eigval, size_t n);

struct md_executor_t;

// The same, its loops spread over an executor (md_executor.h; NULL runs them on the calling thread).
// The result does not depend on the executor or its number of threads: it is the same, bit for bit.
bool md_linalg_sym_eigen_ex(double* a, double* eigval, size_t n, const struct md_executor_t* executor);

// C = alpha op(A) op(B) + beta C, column major (the BLAS dgemm contract): op(A) = A or A^T is m x k,
// op(B) = B or B^T is k x n, C is m x n; lda, ldb, ldc the distances between columns. beta = 0 ignores
// what C held (NaN included). Runs on 'executor' (NULL: the calling thread), and the result does not
// depend on it. Row major operands are their own transposes: for row major C = A B, compute the
// column major C^T = B^T A^T, which is the same memory. False only if the scratch (a few MB per thread)
// could not be allocated, or the executor was cancelled.
bool md_linalg_gemm(bool trans_a, bool trans_b, size_t m, size_t n, size_t k, double alpha, const double* A, size_t lda,
                    const double* B, size_t ldb, double beta, double* C, size_t ldc, const struct md_executor_t* executor);


#ifdef __cplusplus
}
#endif

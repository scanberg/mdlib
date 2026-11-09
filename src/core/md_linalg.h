#pragma once

#include <stdbool.h>
#include <stddef.h>

// Dense linear algebra that more than one part of mdlib needs, in ONE place.
//
// The procedures here are deliberately shaped like their LAPACK counterparts (column major, in
// place, eigenvalues ascending), so that a build which has a LAPACK can route them there without
// any caller noticing. Until then they are plain C reference implementations: correct, robust and
// single threaded, and the place to make them faster when that is needed.

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
// Returns false if the iteration did not converge, which in practice means the input held a NaN or
// an Inf; the contents of 'a' and 'eigval' are then undefined.
//
// Cost is O(n^3) with a modest constant: Householder tridiagonalization and the implicit QL
// algorithm (EISPACK tred2/tql2, as in JAMA). About 0.6 s for n = 1000 and 16 s for n = 2750 on a
// single core; nothing is parallel yet.
bool md_linalg_sym_eigen(double* a, double* eigval, size_t n);

#ifdef __cplusplus
}
#endif

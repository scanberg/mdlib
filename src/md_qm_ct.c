#include <md_qm_ct.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_executor.h>
#include <core/md_linalg.h>
#include <core/md_log.h>

#include <float.h>
#include <math.h>

// ---------------------------------------------------------------------------
// C = A^T B, all row major: A is {k,m}, B is {k,n}, C is {m,n}.
//
// The one product kernel here. Every product this file needs fits it once its operands are read
// the right way round, and in this form both operands are read along their rows: the inner loop
// is an outer product of four entries of a row of A with eight of a row of B, which a compiler
// keeps in registers. Blocked over n and k so the slab of B in use stays in cache; parallel over
// tiles of the rows of C, each written by exactly one piece.
// ---------------------------------------------------------------------------

enum { GEMM_MR = 4, GEMM_NR = 8, GEMM_KB = 256, GEMM_NB = 512 };

typedef struct gemm_tn_args_t {
    double* C;
    const double* A;
    const double* B;
    size_t m, n, k;
} gemm_tn_args_t;

static void gemm_tn_rows(void* user, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const gemm_tn_args_t* g = (const gemm_tn_args_t*)user;
    const size_t m = g->m, n = g->n, k = g->k;
    const size_t i_beg = (size_t)beg * GEMM_MR;
    const size_t i_end = MD_MIN(m, (size_t)end * GEMM_MR);

    for (size_t i = i_beg; i < i_end; ++i) {
        md_memset(g->C + i * n, 0, sizeof(double) * n);
    }

    for (size_t j0 = 0; j0 < n; j0 += GEMM_NB) {
        const size_t j1 = MD_MIN(n, j0 + GEMM_NB);
        for (size_t k0 = 0; k0 < k; k0 += GEMM_KB) {
            const size_t k1 = MD_MIN(k, k0 + GEMM_KB);
            for (size_t i = i_beg; i < i_end; i += GEMM_MR) {
                const size_t mr = MD_MIN((size_t)GEMM_MR, i_end - i);
                for (size_t j = j0; j < j1; j += GEMM_NR) {
                    const size_t nr = MD_MIN((size_t)GEMM_NR, j1 - j);
                    double acc[GEMM_MR][GEMM_NR] = {{0}};
                    if (mr == GEMM_MR && nr == GEMM_NR) {
                        for (size_t p = k0; p < k1; ++p) {
                            const double* a = g->A + p * m + i;
                            const double* b = g->B + p * n + j;
                            for (int r = 0; r < GEMM_MR; ++r) {
                                for (int c = 0; c < GEMM_NR; ++c) {
                                    acc[r][c] += a[r] * b[c];
                                }
                            }
                        }
                    } else {
                        for (size_t p = k0; p < k1; ++p) {
                            const double* a = g->A + p * m + i;
                            const double* b = g->B + p * n + j;
                            for (size_t r = 0; r < mr; ++r) {
                                for (size_t c = 0; c < nr; ++c) {
                                    acc[r][c] += a[r] * b[c];
                                }
                            }
                        }
                    }
                    for (size_t r = 0; r < mr; ++r) {
                        double* dst = g->C + (i + r) * n + j;
                        for (size_t c = 0; c < nr; ++c) dst[c] += acc[r][c];
                    }
                }
            }
        }
    }
}

static bool gemm_tn(double* C, const double* A, const double* B, size_t m, size_t n, size_t k, const md_executor_t* ex) {
    gemm_tn_args_t args = { C, A, B, m, n, k };
    const uint32_t tiles = (uint32_t)((m + GEMM_MR - 1) / GEMM_MR);
    return md_parallel_for(ex, tiles, 4, gemm_tn_rows, &args);
}

static void transpose(double* out, const double* in, size_t rows, size_t cols) {
    enum { TB = 32 };
    for (size_t r0 = 0; r0 < rows; r0 += TB) {
        const size_t r1 = MD_MIN(rows, r0 + TB);
        for (size_t c0 = 0; c0 < cols; c0 += TB) {
            const size_t c1 = MD_MIN(cols, c0 + TB);
            for (size_t r = r0; r < r1; ++r) {
                for (size_t c = c0; c < c1; ++c) {
                    out[c * rows + r] = in[r * cols + c];
                }
            }
        }
    }
}

static double dot(const double* a, const double* b, size_t n) {
    double s[4] = { 0, 0, 0, 0 };
    size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        for (int l = 0; l < 4; ++l) s[l] += a[i + l] * b[i + l];
    }
    double r = (s[0] + s[1]) + (s[2] + s[3]);
    for (; i < n; ++i) r += a[i] * b[i];
    return r;
}

// ---------------------------------------------------------------------------
// Lowdin orbitals
// ---------------------------------------------------------------------------

bool md_qm_lowdin_orbitals(double* out, const double* C, size_t num_mo, size_t num_ao, const md_executor_t* ex) {
    if (!out || !C || num_mo == 0 || num_ao == 0) return false;
    if (num_mo > num_ao) {
        MD_LOG_ERROR("Lowdin orbitals: %zu orbitals over %zu AOs cannot be linearly independent", num_mo, num_ao);
        return false;
    }
    const size_t M = num_mo, A = num_ao;

    md_temp_scope_t temp = md_temp_begin();
    bool ok = false;

    // G = C C^T over the rows, {M,M}: the Gram matrix of the coefficient vectors, which is C^T C in
    // column terms. As A^T B with A = B = C^T.
    double* G  = md_temp_alloc_array(temp, double, M * M);
    double* w  = md_temp_alloc_array(temp, double, M);
    double* Vs = md_temp_alloc_array(temp, double, M * M);
    if (!G || !w || !Vs) {
        MD_LOG_ERROR("Lowdin orbitals: out of memory (%zu orbitals)", M);
        goto done;
    }
    {
        const md_temp_mark_t mark = md_temp_mark(temp);
        double* Ct = md_temp_alloc_array(temp, double, A * M);
        if (!Ct) {
            MD_LOG_ERROR("Lowdin orbitals: out of memory (%zu orbitals)", M);
            goto done;
        }
        transpose(Ct, C, M, A);
        if (!gemm_tn(G, Ct, Ct, M, M, A, ex)) goto done;
        md_temp_reset(temp, mark);
    }

    // G = V diag(w) V^T, eigenvector j as row j of G afterwards
    if (!md_linalg_sym_eigen(G, w, M)) {
        MD_LOG_ERROR("Lowdin orbitals: the eigendecomposition did not converge");
        goto done;
    }
    {
        const double w_max = w[M - 1];
        if (!(w[0] > w_max * (double)M * DBL_EPSILON)) {
            MD_LOG_ERROR("Lowdin orbitals: the orbitals are linearly dependent (Gram eigenvalues %g .. %g)", w[0], w_max);
            goto done;
        }
    }

    // W = G^(-1/2) = V^T diag(w^-1/2) V in row terms. As A^T B with A = V, B = rows of V scaled.
    for (size_t j = 0; j < M; ++j) {
        const double f = 1.0 / sqrt(w[j]);
        const double* v = G + j * M;
        double* s = Vs + j * M;
        for (size_t i = 0; i < M; ++i) s[i] = f * v[i];
    }
    {
        double* W = md_temp_alloc_array(temp, double, M * M);
        if (!W) {
            MD_LOG_ERROR("Lowdin orbitals: out of memory (%zu orbitals)", M);
            goto done;
        }
        if (!gemm_tn(W, G, Vs, M, M, M, ex)) goto done;
        // out = W C, and W is symmetric, so W^T C
        if (!gemm_tn(out, W, C, M, A, M, ex)) goto done;
    }
    ok = true;

done:
    md_temp_end(temp);
    return ok;
}

void md_qm_lowdin_transform(double* out, const double* x, size_t count, const double* S,
                            const double* C, const double* C_lowdin, size_t num_mo, size_t num_ao) {
    if (!out || !x || !S || !C || !C_lowdin || count == 0) return;
    const size_t M = num_mo, A = num_ao;
    md_temp_scope_t temp = md_temp_begin();
    double* y = md_temp_alloc_array(temp, double, A);
    double* z = md_temp_alloc_array(temp, double, M);
    for (size_t v = 0; v < count; ++v) {
        const double* xv = x + v * A;
        double* ov = out + v * A;
        // y = S x (S symmetric: rows), z = C y over the orbitals, out = sum_m z_m C~_m
        for (size_t mu = 0; mu < A; ++mu) y[mu] = dot(S + mu * A, xv, A);
        for (size_t m = 0; m < M; ++m) z[m] = dot(C + m * A, y, A);
        md_memset(ov, 0, sizeof(double) * A);
        for (size_t m = 0; m < M; ++m) {
            const double f = z[m];
            const double* row = C_lowdin + m * A;
            for (size_t mu = 0; mu < A; ++mu) ov[mu] += f * row[mu];
        }
    }
    md_temp_end(temp);
}

void md_qm_lowdin_inverse_transform(double* out, const double* x_lowdin, size_t count,
                                    const double* C, const double* C_lowdin, size_t num_mo, size_t num_ao) {
    if (!out || !x_lowdin || !C || !C_lowdin || count == 0) return;
    const size_t M = num_mo, A = num_ao;
    md_temp_scope_t temp = md_temp_begin();
    double* z = md_temp_alloc_array(temp, double, M);
    for (size_t v = 0; v < count; ++v) {
        const double* xv = x_lowdin + v * A;
        double* ov = out + v * A;
        for (size_t m = 0; m < M; ++m) z[m] = dot(C_lowdin + m * A, xv, A);
        md_memset(ov, 0, sizeof(double) * A);
        for (size_t m = 0; m < M; ++m) {
            const double f = z[m];
            const double* row = C + m * A;
            for (size_t mu = 0; mu < A; ++mu) ov[mu] += f * row[mu];
        }
    }
    md_temp_end(temp);
}

// ---------------------------------------------------------------------------
// Group blocks and Omega
// ---------------------------------------------------------------------------

void md_qm_ct_group_blocks(double* out, const double* x, size_t K, size_t num_ao, const uint32_t* ao_group, size_t num_groups) {
    if (!out || num_groups == 0 || K == 0) return;
    md_memset(out, 0, sizeof(double) * num_groups * K * K);
    if (!x || !ao_group) return;

    for (size_t mu = 0; mu < num_ao; ++mu) {
        const uint32_t g = ao_group[mu];
        if (g >= num_groups) continue;
        double* blk = out + (size_t)g * K * K;
        for (size_t k = 0; k < K; ++k) {
            const double xk = x[k * num_ao + mu];
            if (xk == 0.0) continue;
            double* row = blk + k * K;
            for (size_t l = 0; l <= k; ++l) {
                row[l] += xk * x[l * num_ao + mu];
            }
        }
    }
    // mirror the lower triangle
    for (size_t g = 0; g < num_groups; ++g) {
        double* blk = out + g * K * K;
        for (size_t k = 0; k < K; ++k) {
            for (size_t l = k + 1; l < K; ++l) blk[k * K + l] = blk[l * K + k];
        }
    }
}

void md_qm_ct_omega(double* out, const double* lambda, size_t K,
                    const double* hole_blocks, size_t num_hole_groups,
                    const double* part_blocks, size_t num_part_groups) {
    if (!out || !lambda || !hole_blocks || !part_blocks || K == 0) return;
    const size_t KK = K * K;

    md_temp_scope_t temp = md_temp_begin();
    // H'[a] = (s s^T) o H[a], so that Omega_ab = <H'[a], P[b]>: one dot product of K^2 per entry.
    double* Hs = md_temp_alloc_array(temp, double, num_hole_groups * KK);
    double* s  = md_temp_alloc_array(temp, double, K);
    for (size_t k = 0; k < K; ++k) s[k] = lambda[k] > 0.0 ? sqrt(lambda[k]) : 0.0;
    for (size_t a = 0; a < num_hole_groups; ++a) {
        const double* h = hole_blocks + a * KK;
        double* hs = Hs + a * KK;
        for (size_t k = 0; k < K; ++k) {
            for (size_t l = 0; l < K; ++l) hs[k * K + l] = s[k] * s[l] * h[k * K + l];
        }
    }
    for (size_t a = 0; a < num_hole_groups; ++a) {
        const double* hs = Hs + a * KK;
        double* row = out + a * num_part_groups;
        for (size_t b = 0; b < num_part_groups; ++b) {
            // a sum of squares in exact arithmetic; never let round off say otherwise
            row[b] = MD_MAX(0.0, dot(hs, part_blocks + b * KK, KK));
        }
    }
    md_temp_end(temp);
}

// ---------------------------------------------------------------------------
// Block NTOs
// ---------------------------------------------------------------------------
//
// X = H_a D P_b^T, with H_a the hole NTOs restricted to group a (columns), D = diag(s) and P_b the
// particle NTOs restricted to b. Only their K x K blocks are known: Hb = H_a^T H_a, Pb = P_b^T P_b.
//
//   Hb = Q diag(e) Q^T,  R = diag(e)^(1/2) Q^T   (rows of R for e above round off only)
//   M  = R D Pb D R^T  = Y diag(sigma^2) Y^T     (r x r, symmetric)
//
// The nonzero singular values of X are sigma; the left vector is u = H_a alpha with
// alpha = Q diag(e)^(-1/2) y, and the right v = P_b beta with beta = D R^T y / sigma. Both of unit norm,
// and tr M = tr(D Pb D Hb) = Omega_ab.

size_t md_qm_ct_block_nto(double* out_weight, double* out_hole_coef, double* out_part_coef,
                          const double* lambda, size_t K,
                          const double* hole_block_a, const double* part_block_b) {
    if (!out_weight || !lambda || !hole_block_a || !part_block_b || K == 0) return 0;
    for (size_t k = 0; k < K; ++k) out_weight[k] = 0.0;

    md_temp_scope_t temp = md_temp_begin();
    size_t count = 0;

    double* s  = md_temp_alloc_array(temp, double, K);
    double* Q  = md_temp_alloc_array(temp, double, K * K);
    double* e  = md_temp_alloc_array(temp, double, K);
    double* R  = md_temp_alloc_array(temp, double, K * K);   // {r,K}
    double* T  = md_temp_alloc_array(temp, double, K * K);   // D Pb D R^T, {K,r}
    double* Mm = md_temp_alloc_array(temp, double, K * K);   // {r,r}
    double* sig2 = md_temp_alloc_array(temp, double, K);

    for (size_t k = 0; k < K; ++k) s[k] = lambda[k] > 0.0 ? sqrt(lambda[k]) : 0.0;

    md_memcpy(Q, hole_block_a, sizeof(double) * K * K);
    if (!md_linalg_sym_eigen(Q, e, K)) goto done;

    // the hole block's range: eigenvalues above round off of its largest
    const double e_max = e[K - 1];
    if (!(e_max > 0.0)) goto done;
    const double e_tol = e_max * (double)K * 1.0e3 * DBL_EPSILON;
    size_t r = 0;
    size_t* idx = md_temp_alloc_array(temp, size_t, K);
    for (size_t j = K; j-- > 0;) {
        if (e[j] > e_tol) idx[r++] = j;
    }
    for (size_t q = 0; q < r; ++q) {
        const double f = sqrt(e[idx[q]]);
        const double* v = Q + idx[q] * K;          // eigenvector idx[q], contiguous
        for (size_t k = 0; k < K; ++k) R[q * K + k] = f * v[k];
    }

    // T[k][q] = s_k sum_l Pb[k][l] s_l R[q][l]
    for (size_t k = 0; k < K; ++k) {
        const double* pk = part_block_b + k * K;
        for (size_t q = 0; q < r; ++q) {
            const double* rq = R + q * K;
            double acc = 0.0;
            for (size_t l = 0; l < K; ++l) acc += pk[l] * s[l] * rq[l];
            T[k * r + q] = s[k] * acc;
        }
    }
    // M[p][q] = sum_k R[p][k] T[k][q]
    for (size_t p = 0; p < r; ++p) {
        for (size_t q = 0; q <= p; ++q) {
            double acc = 0.0;
            for (size_t k = 0; k < K; ++k) acc += R[p * K + k] * T[k * r + q];
            Mm[p * r + q] = acc;
            Mm[q * r + p] = acc;
        }
    }
    if (!md_linalg_sym_eigen(Mm, sig2, r)) goto done;

    double total = 0.0;
    for (size_t j = 0; j < r; ++j) total += MD_MAX(0.0, sig2[j]);
    const double w_tol = MD_MAX(total, 0.0) * (double)K * 1.0e3 * DBL_EPSILON;

    for (size_t jj = r; jj-- > 0;) {        // descending
        const double w2 = sig2[jj];
        if (!(w2 > w_tol)) break;
        out_weight[count] = w2;
        const double* y = Mm + jj * r;      // eigenvector jj, contiguous
        if (out_hole_coef) {
            double* alpha = out_hole_coef + count * K;
            for (size_t k = 0; k < K; ++k) alpha[k] = 0.0;
            for (size_t q = 0; q < r; ++q) {
                const double f = y[q] / sqrt(e[idx[q]]);
                const double* v = Q + idx[q] * K;
                for (size_t k = 0; k < K; ++k) alpha[k] += f * v[k];
            }
        }
        if (out_part_coef) {
            double* beta = out_part_coef + count * K;
            const double inv = 1.0 / sqrt(w2);
            for (size_t k = 0; k < K; ++k) {
                double acc = 0.0;
                for (size_t q = 0; q < r; ++q) acc += R[q * K + k] * y[q];
                beta[k] = s[k] * acc * inv;
            }
        }
        count++;
    }

done:
    md_temp_end(temp);
    return count;
}

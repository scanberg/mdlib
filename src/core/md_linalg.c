#include <core/md_linalg.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_simd.h>
#include <core/md_executor.h>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

// The symmetric eigensolver has LAPACK dsyevd's structure, which is what makes it fast: almost all of
// its O(n^3) work is matrix products, not the rotations and rank-one updates of the textbook method.
//
//   1. Householder reduction to tridiagonal form, T = Q^T A Q, Q kept as its reflectors (dsytd2).
//   2. The eigenvectors Z of T by divide and conquer (Cuppen; dstedc/dlaed*): T split in two halves,
//      each solved recursively, and the halves joined by solving the rank-one update of their
//      eigenvalues. The new eigenvectors are the halves' eigenvectors times a small dense matrix,
//      which is where the time goes: one product per join. Deflation removes the eigenpairs a join
//      does not change (common, and the reason this is often far cheaper than its bound), and the
//      Gu-Eisenstat recomputation of the update vector keeps the eigenvectors orthogonal to working
//      precision however close the eigenvalues are. The leaves are solved by implicit QL (tql2).
//   3. The eigenvectors of A, Q Z, by applying the reflectors in blocks (compact WY), again as
//      matrix products.
//
// Column major throughout: A(r, c) = a[c * n + r], so eigenvector j is contiguous at a[j * n].

#define LA_EPS (0.5 * DBL_EPSILON)   // unit roundoff, LAPACK's dlamch('E')

// ---------------------------------------------------------------------------------------------------
// The one product kernel: C = alpha op(A) B + beta C, column major, op(A) = A or A^T. Packed panels
// and an 8 x 6 register tile (GotoBLAS structure) so that it runs near the core's FMA rate whatever
// the operand layout; md_simd maps it to AVX2/AVX-512 or NEON.
// ---------------------------------------------------------------------------------------------------

enum { LA_MR = 8, LA_NR = 6, LA_MC = 96, LA_KC = 256, LA_NC = 768 };

static void la_pack_a(double* dst, bool ta, const double* A, size_t lda, size_t i0, size_t mc, size_t p0, size_t kc) {
    for (size_t ir = 0; ir < mc; ir += LA_MR) {
        const size_t mr = MD_MIN((size_t)LA_MR, mc - ir);
        double* d = dst + ir * kc;
        if (!ta) {
            for (size_t p = 0; p < kc; ++p) {
                const double* s = A + (p0 + p) * lda + i0 + ir;
                size_t r = 0;
                for (; r < mr; ++r) d[p * LA_MR + r] = s[r];
                for (; r < LA_MR; ++r) d[p * LA_MR + r] = 0.0;
            }
        } else {
            for (size_t r = 0; r < LA_MR; ++r) {
                if (r < mr) {
                    const double* s = A + (i0 + ir + r) * lda + p0;
                    for (size_t p = 0; p < kc; ++p) d[p * LA_MR + r] = s[p];
                } else {
                    for (size_t p = 0; p < kc; ++p) d[p * LA_MR + r] = 0.0;
                }
            }
        }
    }
}

static void la_pack_b(double* dst, bool tb, const double* B, size_t ldb, size_t p0, size_t kc, size_t j0, size_t nc) {
    for (size_t jr = 0; jr < nc; jr += LA_NR) {
        const size_t nr = MD_MIN((size_t)LA_NR, nc - jr);
        double* d = dst + jr * kc;
        if (tb) {
            for (size_t p = 0; p < kc; ++p) {
                const double* s = B + (p0 + p) * ldb + j0 + jr;
                size_t c = 0;
                for (; c < nr; ++c) d[p * LA_NR + c] = s[c];
                for (; c < LA_NR; ++c) d[p * LA_NR + c] = 0.0;
            }
            continue;
        }
        for (size_t c = 0; c < LA_NR; ++c) {
            if (c < nr) {
                const double* s = B + (j0 + jr + c) * ldb + p0;
                for (size_t p = 0; p < kc; ++p) d[p * LA_NR + c] = s[p];
            } else {
                for (size_t p = 0; p < kc; ++p) d[p * LA_NR + c] = 0.0;
            }
        }
    }
}

static void la_kernel(size_t kc, const double* Ap, const double* Bp, double* C, size_t ldc, size_t mr, size_t nr, double alpha) {
    md_256d c00 = md_mm256_setzero_pd(), c01 = md_mm256_setzero_pd();
    md_256d c10 = md_mm256_setzero_pd(), c11 = md_mm256_setzero_pd();
    md_256d c20 = md_mm256_setzero_pd(), c21 = md_mm256_setzero_pd();
    md_256d c30 = md_mm256_setzero_pd(), c31 = md_mm256_setzero_pd();
    md_256d c40 = md_mm256_setzero_pd(), c41 = md_mm256_setzero_pd();
    md_256d c50 = md_mm256_setzero_pd(), c51 = md_mm256_setzero_pd();
    for (size_t p = 0; p < kc; ++p) {
        const md_256d a0 = md_mm256_loadu_pd(Ap);
        const md_256d a1 = md_mm256_loadu_pd(Ap + 4);
        md_256d b;
        b = md_mm256_set1_pd(Bp[0]); c00 = md_mm256_fmadd_pd(a0, b, c00); c01 = md_mm256_fmadd_pd(a1, b, c01);
        b = md_mm256_set1_pd(Bp[1]); c10 = md_mm256_fmadd_pd(a0, b, c10); c11 = md_mm256_fmadd_pd(a1, b, c11);
        b = md_mm256_set1_pd(Bp[2]); c20 = md_mm256_fmadd_pd(a0, b, c20); c21 = md_mm256_fmadd_pd(a1, b, c21);
        b = md_mm256_set1_pd(Bp[3]); c30 = md_mm256_fmadd_pd(a0, b, c30); c31 = md_mm256_fmadd_pd(a1, b, c31);
        b = md_mm256_set1_pd(Bp[4]); c40 = md_mm256_fmadd_pd(a0, b, c40); c41 = md_mm256_fmadd_pd(a1, b, c41);
        b = md_mm256_set1_pd(Bp[5]); c50 = md_mm256_fmadd_pd(a0, b, c50); c51 = md_mm256_fmadd_pd(a1, b, c51);
        Ap += LA_MR;
        Bp += LA_NR;
    }
    const md_256d al = md_mm256_set1_pd(alpha);
    if (mr == LA_MR && nr == LA_NR) {
#define LA_STORE(c, lo, hi) { double* dst = C + (c) * ldc; \
        md_mm256_storeu_pd(dst,     md_mm256_fmadd_pd(lo, al, md_mm256_loadu_pd(dst))); \
        md_mm256_storeu_pd(dst + 4, md_mm256_fmadd_pd(hi, al, md_mm256_loadu_pd(dst + 4))); }
        LA_STORE(0, c00, c01) LA_STORE(1, c10, c11) LA_STORE(2, c20, c21)
        LA_STORE(3, c30, c31) LA_STORE(4, c40, c41) LA_STORE(5, c50, c51)
#undef LA_STORE
    } else {
        double t[LA_NR][LA_MR];
        md_mm256_storeu_pd(t[0], c00); md_mm256_storeu_pd(t[0] + 4, c01);
        md_mm256_storeu_pd(t[1], c10); md_mm256_storeu_pd(t[1] + 4, c11);
        md_mm256_storeu_pd(t[2], c20); md_mm256_storeu_pd(t[2] + 4, c21);
        md_mm256_storeu_pd(t[3], c30); md_mm256_storeu_pd(t[3] + 4, c31);
        md_mm256_storeu_pd(t[4], c40); md_mm256_storeu_pd(t[4] + 4, c41);
        md_mm256_storeu_pd(t[5], c50); md_mm256_storeu_pd(t[5] + 4, c51);
        for (size_t c = 0; c < nr; ++c) {
            for (size_t r = 0; r < mr; ++r) C[c * ldc + r] += alpha * t[c][r];
        }
    }
}

// Doubles of packing scratch la_gemm_serial needs, at most, for products of at most n in each dimension
static size_t la_pack_size(size_t n) {
    const size_t mc = MD_MIN((size_t)LA_MC, (n + LA_MR - 1) / LA_MR * LA_MR);
    const size_t nc = MD_MIN((size_t)LA_NC, (n + LA_NR - 1) / LA_NR * LA_NR);
    const size_t kc = MD_MIN((size_t)LA_KC, n);
    return mc * kc + kc * nc;
}

// C (m x n) = alpha op(A) op(B) + beta C; op(A) is m x k (A is k x m when ta), op(B) is k x n (B is n x k
// when tb). pack: la_pack_size doubles. Every entry of C is summed in the same order however the product
// is split into blocks.
static void la_gemm_serial(bool ta, bool tb, size_t m, size_t n, size_t k, double alpha, const double* A, size_t lda,
                           const double* B, size_t ldb, double beta, double* C, size_t ldc, double* pack) {
    if (m == 0 || n == 0) return;
    if (beta != 1.0) {
        for (size_t j = 0; j < n; ++j) {
            double* c = C + j * ldc;
            if (beta == 0.0) memset(c, 0, sizeof(double) * m);
            else for (size_t i = 0; i < m; ++i) c[i] *= beta;
        }
    }
    if (k == 0 || alpha == 0.0) return;

    const size_t mc_max = MD_MIN((size_t)LA_MC, (m + LA_MR - 1) / LA_MR * LA_MR);
    const size_t kc_max = MD_MIN((size_t)LA_KC, k);
    double* pa = pack;
    double* pb = pack + mc_max * kc_max;

    for (size_t jc = 0; jc < n; jc += LA_NC) {
        const size_t nc = MD_MIN((size_t)LA_NC, n - jc);
        for (size_t pc = 0; pc < k; pc += LA_KC) {
            const size_t kc = MD_MIN((size_t)LA_KC, k - pc);
            la_pack_b(pb, tb, B, ldb, pc, kc, jc, nc);
            for (size_t ic = 0; ic < m; ic += LA_MC) {
                const size_t mc = MD_MIN((size_t)LA_MC, m - ic);
                la_pack_a(pa, ta, A, lda, ic, mc, pc, kc);
                for (size_t jr = 0; jr < nc; jr += LA_NR) {
                    for (size_t ir = 0; ir < mc; ir += LA_MR) {
                        la_kernel(kc, pa + ir * kc, pb + jr * kc, C + (jc + jr) * ldc + ic + ir, ldc,
                                  MD_MIN((size_t)LA_MR, mc - ir), MD_MIN((size_t)LA_NR, nc - jr), alpha);
                    }
                }
            }
        }
    }
}

// The scratch every parallel piece of the solver may use: per thread_idx
typedef struct la_par_t {
    const md_executor_t* ex;
    uint32_t num_threads;
    double*  pack;           // [num_threads * pack_size]
    size_t   pack_size;
    double*  vec;            // [num_threads * vec_size] one vector per thread
    size_t   vec_size;
} la_par_t;

typedef struct la_gemm_args_t {
    bool ta, tb, split_rows;
    size_t m, n, k, chunk;
    double alpha, beta;
    const double* A; size_t lda;
    const double* B; size_t ldb;
    double* C; size_t ldc;
    const la_par_t* par;
} la_gemm_args_t;

static void la_gemm_piece(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    const la_gemm_args_t* g = (const la_gemm_args_t*)arg;
    double* pack = g->par->pack + (size_t)thread_idx * g->par->pack_size;
    if (g->split_rows) {
        const size_t i0 = (size_t)beg * g->chunk, i1 = MD_MIN(g->m, (size_t)end * g->chunk);
        const double* A = g->ta ? g->A + i0 * g->lda : g->A + i0;
        la_gemm_serial(g->ta, g->tb, i1 - i0, g->n, g->k, g->alpha, A, g->lda, g->B, g->ldb, g->beta, g->C + i0, g->ldc, pack);
    } else {
        const size_t j0 = (size_t)beg * g->chunk, j1 = MD_MIN(g->n, (size_t)end * g->chunk);
        const double* B = g->tb ? g->B + j0 : g->B + j0 * g->ldb;
        la_gemm_serial(g->ta, g->tb, g->m, j1 - j0, g->k, g->alpha, g->A, g->lda, B, g->ldb, g->beta, g->C + j0 * g->ldc, g->ldc, pack);
    }
}

// la_gemm_serial over the executor, C split by rows or columns, whichever is longer
static bool la_gemm_tt(const la_par_t* par, bool ta, bool tb, size_t m, size_t n, size_t k, double alpha, const double* A, size_t lda,
                       const double* B, size_t ldb, double beta, double* C, size_t ldc) {
    if (par->num_threads == 1 || (double)m * (double)n * (double)k < 1.0e6) {
        la_gemm_serial(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, par->pack);
        return true;
    }
    la_gemm_args_t g = { ta, tb, m >= n, m, n, k, 0, alpha, beta, A, lda, B, ldb, C, ldc, par };
    const size_t len = g.split_rows ? m : n;
    const size_t unit = g.split_rows ? LA_MC : 8 * LA_NR;
    const size_t want = (len + 4 * par->num_threads - 1) / (4 * par->num_threads);   // ~4 pieces per thread
    g.chunk = MD_MAX(unit, (want + unit - 1) / unit * unit);
    return md_parallel_for(par->ex, (uint32_t)((len + g.chunk - 1) / g.chunk), 1, la_gemm_piece, &g);
}

static inline bool la_gemm(const la_par_t* par, bool ta, size_t m, size_t n, size_t k, double alpha, const double* A, size_t lda,
                           const double* B, size_t ldb, double beta, double* C, size_t ldc) {
    return la_gemm_tt(par, ta, false, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

bool md_linalg_gemm(bool trans_a, bool trans_b, size_t m, size_t n, size_t k, double alpha, const double* A, size_t lda,
                    const double* B, size_t ldb, double beta, double* C, size_t ldc, const struct md_executor_t* executor) {
    if (m == 0 || n == 0) return true;
    if (!A || !B || !C) return false;
    const uint32_t nt = md_executor_num_threads(executor);
    const size_t pack = la_pack_size(MD_MAX(m, MD_MAX(n, k)));
    const size_t bytes = sizeof(double) * (size_t)nt * pack;
    md_allocator_i* heap = md_get_heap_allocator();
    double* mem = (double*)md_alloc(heap, bytes);
    if (!mem) return false;
    la_par_t par = { executor, nt, mem, pack, NULL, 0 };
    const bool ok = la_gemm_tt(&par, trans_a, trans_b, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    md_free(heap, mem, bytes);
    return ok;
}

// ---------------------------------------------------------------------------------------------------
// 1. Householder tridiagonalization of the lower triangle (dsytd2). Reflector i is H_i = I - tau_i v v^T
// with v[i+1] = 1 and v[i+2..n) stored in A(i+2..n, i); Q = H_0 H_1 ... H_{n-2}. d gets the diagonal,
// e[i] the coupling of i and i+1. The trailing update of step i and the matrix-vector product of step
// i+1 share one pass over the matrix, which is what this costs: memory traffic, not arithmetic.
// ---------------------------------------------------------------------------------------------------

static double la_nrm2(const double* x, size_t n) {
    double amax = 0.0;
    for (size_t i = 0; i < n; ++i) amax = MD_MAX(amax, fabs(x[i]));
    if (amax == 0.0 || !isfinite(amax)) return amax;
    const double s = 1.0 / amax;
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double t = x[i] * s;
        sum += t * t;
    }
    return amax * sqrt(sum);
}

static inline double la_hsum(md_256d v) {
    double t[4];
    md_mm256_storeu_pd(t, v);
    return (t[0] + t[1]) + (t[2] + t[3]);
}

// Turns x (length m) into a reflector: x[0] = beta, x[1..m) the tail of v. Returns tau.
static double la_householder(double* x, size_t m) {
    const double alpha = x[0];
    const double xnorm = m > 1 ? la_nrm2(x + 1, m - 1) : 0.0;
    if (xnorm == 0.0) return 0.0;
    const double beta = -copysign(hypot(alpha, xnorm), alpha);
    const double s = 1.0 / (alpha - beta);
    for (size_t i = 1; i < m; ++i) x[i] *= s;
    x[0] = beta;
    return (beta - alpha) / beta;
}

// One pass of a step over the trailing block A(i+1.., i+1..) of size m, in pieces of columns: finish the
// pending update of each column (pv, pw; none if NULL), then add its part of y = A v. Each piece sums
// into its own partial y, so the pieces are independent and the sum is the same whatever ran where.
enum { LA_TRD_PIECES = 16, LA_TRD_COLS = 192 };

typedef struct la_trd_args_t {
    double*         T;      // &A(i+1, i+1)
    size_t          ld;
    size_t          m;
    const double*   v;
    const double*   pv;
    const double*   pw;
    double*         ypart;  // piece p's partial at ypart + p * ld, rows [cb[p], m)
    const uint32_t* cb;     // [P + 1] the pieces' first columns
} la_trd_args_t;

static void la_trd_piece(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const la_trd_args_t* t = (const la_trd_args_t*)arg;
    const size_t m = t->m;
    const double* v = t->v;
    const double* pv = t->pv;
    const double* pw = t->pw;
    for (uint32_t p = beg; p < end; ++p) {
        const size_t c0 = t->cb[p], c1 = t->cb[p + 1];
        double* yv = t->ypart + (size_t)p * t->ld;
        memset(yv + c0, 0, sizeof(double) * (m - c0));
        for (size_t c = c0; c < c1; ++c) {
            double* col = t->T + c * t->ld;     // col[r] = A(i+1+r, i+1+c), r >= c
            const double vc = v[c];
            md_256d sacc = md_mm256_setzero_pd();
            const md_256d bv = md_mm256_set1_pd(vc);
            size_t r = c;
            double s;
            if (pv) {
                const double uw = pw[c], uv = pv[c];
                col[r] -= pv[r] * uw + pw[r] * uv;
                s = col[r] * vc;
                ++r;
                const md_256d buw = md_mm256_set1_pd(uw), buv = md_mm256_set1_pd(uv);
                for (; r + 4 <= m; r += 4) {
                    md_256d cc = md_mm256_loadu_pd(col + r);
                    cc = md_mm256_fnmadd_pd(md_mm256_loadu_pd(pv + r), buw, cc);
                    cc = md_mm256_fnmadd_pd(md_mm256_loadu_pd(pw + r), buv, cc);
                    md_mm256_storeu_pd(col + r, cc);
                    sacc = md_mm256_fmadd_pd(cc, md_mm256_loadu_pd(v + r), sacc);
                    md_mm256_storeu_pd(yv + r, md_mm256_fmadd_pd(cc, bv, md_mm256_loadu_pd(yv + r)));
                }
                for (; r < m; ++r) {
                    col[r] -= pv[r] * uw + pw[r] * uv;
                    s += col[r] * v[r];
                    yv[r] += col[r] * vc;
                }
            } else {
                s = col[r] * vc;
                ++r;
                for (; r + 4 <= m; r += 4) {
                    const md_256d cc = md_mm256_loadu_pd(col + r);
                    sacc = md_mm256_fmadd_pd(cc, md_mm256_loadu_pd(v + r), sacc);
                    md_mm256_storeu_pd(yv + r, md_mm256_fmadd_pd(cc, bv, md_mm256_loadu_pd(yv + r)));
                }
                for (; r < m; ++r) {
                    s += col[r] * v[r];
                    yv[r] += col[r] * vc;
                }
            }
            yv[c] += s + la_hsum(sacc);
        }
    }
}

// Pending update of the whole trailing block without a product (a step whose reflector is the identity)
static void la_trd_update_only(double* T, size_t ld, size_t m, const double* pv, const double* pw) {
    for (size_t c = 0; c < m; ++c) {
        double* col = T + c * ld;
        const double vc = pv[c], wc = pw[c];
        for (size_t r = c; r < m; ++r) col[r] -= pv[r] * wc + pw[r] * vc;
    }
}

// y, w: [2n] each; ypart: [LA_TRD_PIECES * n]
static bool la_tridiag(double* a, size_t n, double* d, double* e, double* tau, double* y, double* w, double* ypart, const la_par_t* par) {
#define A_(r, c) a[(size_t)(c) * n + (size_t)(r)]
    // Pending rank-2 update of the trailing matrix from the previous step: A -= pv pw^T + pw pv^T over
    // rows/columns [i, n). Applied lazily, column by column, inside the next step's pass.
    const double* pv = NULL;
    const double* pw = NULL;
    // v and w of this and the previous step, alternating
    double* vbuf[2] = { y, y + n };
    double* wbuf[2] = { w, w + n };
    uint32_t cb[LA_TRD_PIECES + 1];

    for (size_t i = 0; i + 1 < n; ++i) {
        const size_t m = n - i - 1;           // the reflector's length, rows i+1 .. n-1
        // the previous step's update of column i (rows i..n), the column this step reduces
        if (pv) {
            const double vi = pv[0], wi = pw[0];
            double* col = &A_(i, i);
            for (size_t r = 0; r <= m; ++r) col[r] -= pv[r] * wi + pw[r] * vi;
            pv += 1;      // the update's vectors now start at row i+1
            pw += 1;
        }
        d[i] = A_(i, i);
        double* x = &A_(i + 1, i);
        const double t = la_householder(x, m);
        e[i] = x[0];
        tau[i] = t;
        if (t == 0.0) {
            if (pv) la_trd_update_only(&A_(i + 1, i + 1), n, m, pv, pw);
            pv = pw = NULL;
            continue;
        }

        // v, explicit with its leading 1
        double* v = vbuf[i & 1];
        v[0] = 1.0;
        memcpy(v + 1, x + 1, sizeof(double) * (m - 1));

        // the pieces: equal areas of the lower triangle, fixed by m alone
        const size_t P = MD_CLAMP(m / LA_TRD_COLS, (size_t)1, (size_t)LA_TRD_PIECES);
        cb[0] = 0;
        for (size_t p = 1; p < P; ++p) {
            const size_t c = m - (size_t)((double)m * sqrt(1.0 - (double)p / (double)P) + 0.5);
            cb[p] = (uint32_t)MD_CLAMP(c, (size_t)cb[p - 1], m);
        }
        cb[P] = (uint32_t)m;
        la_trd_args_t args = { &A_(i + 1, i + 1), n, m, v, pv, pw, ypart, cb };
        if (P == 1) la_trd_piece(&args, 0, 1, 0);
        else if (!md_parallel_for(par->ex, (uint32_t)P, 1, la_trd_piece, &args)) return false;

        // y = the pieces' partials, in order
        double* yv = wbuf[i & 1];
        memcpy(yv, ypart, sizeof(double) * m);
        for (size_t p = 1; p < P; ++p) {
            const double* yp = ypart + p * n;
            for (size_t r = cb[p]; r < m; ++r) yv[r] += yp[r];
        }

        // w = tau y - (tau^2 / 2) (y . v) v, and the update A -= v w^T + w v^T is pending
        double yv_dot = 0.0;
        for (size_t r = 0; r < m; ++r) {
            yv[r] *= t;
            yv_dot += yv[r] * v[r];
        }
        const double al = -0.5 * t * yv_dot;
        for (size_t r = 0; r < m; ++r) yv[r] += al * v[r];
        pv = v;
        pw = yv;
    }
    d[n - 1] = A_(n - 1, n - 1);
    e[n - 1] = 0.0;
    return true;
#undef A_
}

// ---------------------------------------------------------------------------------------------------
// 3. Z <- Q Z, the reflectors applied in blocks of LA_NB as I - Y T Y^T (dlarft forward/columnwise,
// dlarfb), from the last block to the first. A block is thin (LA_NB columns of Y), so as products its
// cost would be moving Z through memory, several times; instead it is applied to a few columns of Z
// at a time while they are in cache: w = Y^T z, w = T w, z -= Y w. One read and one write of Z per
// block, and the columns are independent.
// ---------------------------------------------------------------------------------------------------

enum { LA_NB = 32, LA_QG = 4 };   // reflectors per block, columns of Z per group

// Zg -= Y W^T for g <= LA_QG columns of Zg; Wt holds W transposed (row c of Wt is column c of W, ld ldw)
static void la_q_update_cols(double* Zg, size_t ldz, const double* Y, size_t nr, size_t kb, const double* Wt, size_t ldw, size_t g) {
    double* z[LA_QG];
    for (size_t c = 0; c < LA_QG; ++c) z[c] = Zg + MD_MIN(c, g - 1) * ldz;
    double w[LA_QG][LA_NB];
    for (size_t c = 0; c < LA_QG; ++c) {
        for (size_t q = 0; q < kb; ++q) w[c][q] = Wt[MD_MIN(c, g - 1) + q * ldw];
    }
    size_t r = 0;
    for (; r + 8 <= nr; r += 8) {
        md_256d z00 = md_mm256_loadu_pd(z[0] + r), z01 = md_mm256_loadu_pd(z[0] + r + 4);
        md_256d z10 = md_mm256_loadu_pd(z[1] + r), z11 = md_mm256_loadu_pd(z[1] + r + 4);
        md_256d z20 = md_mm256_loadu_pd(z[2] + r), z21 = md_mm256_loadu_pd(z[2] + r + 4);
        md_256d z30 = md_mm256_loadu_pd(z[3] + r), z31 = md_mm256_loadu_pd(z[3] + r + 4);
        // rows r .. r+8 of Y[:, q] are zero for q > r + 7: those columns add nothing
        const size_t qe = MD_MIN(kb, r + 8);
        for (size_t q = 0; q < qe; ++q) {
            const md_256d y0 = md_mm256_loadu_pd(Y + q * nr + r), y1 = md_mm256_loadu_pd(Y + q * nr + r + 4);
            md_256d b;
            b = md_mm256_set1_pd(w[0][q]); z00 = md_mm256_fnmadd_pd(y0, b, z00); z01 = md_mm256_fnmadd_pd(y1, b, z01);
            b = md_mm256_set1_pd(w[1][q]); z10 = md_mm256_fnmadd_pd(y0, b, z10); z11 = md_mm256_fnmadd_pd(y1, b, z11);
            b = md_mm256_set1_pd(w[2][q]); z20 = md_mm256_fnmadd_pd(y0, b, z20); z21 = md_mm256_fnmadd_pd(y1, b, z21);
            b = md_mm256_set1_pd(w[3][q]); z30 = md_mm256_fnmadd_pd(y0, b, z30); z31 = md_mm256_fnmadd_pd(y1, b, z31);
        }
        // a column repeated to fill the group computes the same values: storing it twice is harmless
        md_mm256_storeu_pd(z[3] + r, z30); md_mm256_storeu_pd(z[3] + r + 4, z31);
        md_mm256_storeu_pd(z[2] + r, z20); md_mm256_storeu_pd(z[2] + r + 4, z21);
        md_mm256_storeu_pd(z[1] + r, z10); md_mm256_storeu_pd(z[1] + r + 4, z11);
        md_mm256_storeu_pd(z[0] + r, z00); md_mm256_storeu_pd(z[0] + r + 4, z01);
    }
    for (; r < nr; ++r) {
        const size_t qe = MD_MIN(kb, r + 1);
        for (size_t c = 0; c < g; ++c) {
            double s = 0.0;
            for (size_t q = 0; q < qe; ++q) s += Y[q * nr + r] * w[c][q];
            z[c][r] -= s;
        }
    }
}

typedef struct la_aq_args_t {
    double* Z;
    size_t n, r0, nr, kb;
    const double* Y;
    const double* T;
    double* W;
} la_aq_args_t;

// Groups of LA_QG columns: their rows of W^T times T^T, then z -= Y w
static void la_aq_piece(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const la_aq_args_t* q = (const la_aq_args_t*)arg;
    const size_t n = q->n, kb = q->kb;
    for (size_t grp = beg; grp < end; ++grp) {
        const size_t j0 = grp * LA_QG;
        const size_t g = MD_MIN((size_t)LA_QG, n - j0);
        for (size_t j = j0; j < j0 + g; ++j) {
            double wj[LA_NB];
            for (size_t c = 0; c < kb; ++c) wj[c] = q->W[j + c * n];
            for (size_t c = 0; c < kb; ++c) {
                double s = 0.0;
                for (size_t p = c; p < kb; ++p) s += q->T[c + p * kb] * wj[p];
                q->W[j + c * n] = s;
            }
        }
        la_q_update_cols(q->Z + j0 * n + q->r0, n, q->Y, q->nr, kb, q->W + j0, n, g);
    }
}

static bool la_apply_q(const double* a, size_t n, const double* tau, double* Z, double* Y, double* T, double* W, const la_par_t* par) {
    if (n < 3) return true;     // reflector n-2 and beyond have length 1: none for n < 3
    const size_t nrefl = n - 2;
    const size_t nblocks = (nrefl + LA_NB - 1) / LA_NB;
    for (size_t b = nblocks; b-- > 0;) {
        const size_t i0 = b * LA_NB;
        const size_t kb = MD_MIN((size_t)LA_NB, nrefl - i0);
        const size_t r0 = i0 + 1;
        const size_t nr = n - r0;

        for (size_t c = 0; c < kb; ++c) {
            const size_t i = i0 + c;
            double* yc = Y + c * nr;
            const double* src = a + i * n;
            for (size_t r = 0; r < nr; ++r) {
                const size_t row = r0 + r;
                yc[r] = row < i + 1 ? 0.0 : (row == i + 1 ? 1.0 : src[row]);
            }
        }
        // T upper triangular: T[0:c, c] = -tau_c T[0:c, 0:c] Y[:, 0:c]^T Y[:, c]
        for (size_t c = 0; c < kb; ++c) {
            const double tc = tau[i0 + c];
            T[c + c * kb] = tc;
            if (c == 0) continue;
            double* tcol = T + c * kb;
            const double* yc = Y + c * nr;
            for (size_t q = 0; q < c; ++q) {
                const double* yq = Y + q * nr;
                double s = 0.0;
                for (size_t r = c; r < nr; ++r) s += yq[r] * yc[r];   // Y[:, c] is zero above row c
                tcol[q] = -tc * s;
            }
            for (size_t q = 0; q < c; ++q) {
                double s = 0.0;
                for (size_t p = q; p < c; ++p) s += T[q + p * kb] * tcol[p];
                tcol[q] = s;
            }
        }
        // W^T = Z^T Y (n x kb, ld n): as a product with Y the small packed operand, Z read once
        if (!la_gemm(par, true, n, kb, nr, 1.0, Z + r0, n, Y, nr, 0.0, W, n)) return false;
        la_aq_args_t args = { Z, n, r0, nr, kb, Y, T, W };
        if (!md_parallel_for(par->ex, (uint32_t)((n + LA_QG - 1) / LA_QG), 16, la_aq_piece, &args)) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------------------------------
// 2. Divide and conquer for the symmetric tridiagonal eigenproblem.
// ---------------------------------------------------------------------------------------------------

enum { LA_LEAF = 25 };

// Implicit QL with the eigenvectors (EISPACK tql2): d the diagonal, e[i] the coupling of i and i+1
// (e[m-1] = 0, used as scratch), z an m x m block (leading dimension ldz) to rotate. Sorted ascending.
static bool la_tql2(double* d, double* e, double* z, size_t m, size_t ldz) {
#define Z_(r, c) z[(size_t)(c) * ldz + (size_t)(r)]
    double f = 0.0, tst1 = 0.0;
    const double eps = DBL_EPSILON;
    e[m - 1] = 0.0;
    for (size_t l = 0; l < m; ++l) {
        tst1 = MD_MAX(tst1, fabs(d[l]) + fabs(e[l]));
        size_t mm = l;
        while (mm < m && fabs(e[mm]) > eps * tst1) ++mm;
        if (mm == m) mm = m - 1;
        if (mm > l) {
            int iter = 0;
            do {
                if (++iter > 64) return false;
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (size_t i = l + 2; i < m; ++i) d[i] -= h;
                f += h;

                p = d[mm];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (size_t ii = mm; ii-- > l;) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[ii];
                    h = c * p;
                    r = hypot(p, e[ii]);
                    e[ii + 1] = s * r;
                    s = e[ii] / r;
                    c = p / r;
                    p = c * d[ii] - s * g;
                    d[ii + 1] = h + s * (c * g + s * d[ii]);
                    double* x = &Z_(0, ii);
                    double* y = &Z_(0, ii + 1);
                    for (size_t k = 0; k < m; ++k) {
                        const double t = y[k];
                        y[k] = s * x[k] + c * t;
                        x[k] = c * x[k] - s * t;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
            } while (fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }
    for (size_t i = 0; i + 1 < m; ++i) {
        size_t k = i;
        double p = d[i];
        for (size_t j = i + 1; j < m; ++j) {
            if (d[j] < p) { k = j; p = d[j]; }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            double* x = &Z_(0, i);
            double* y = &Z_(0, k);
            for (size_t j = 0; j < m; ++j) { const double t = x[j]; x[j] = y[j]; y[j] = t; }
        }
    }
#undef Z_
    return true;
}

// The j-th root of the secular equation 1/rho + sum_i z_i^2 / (dl_i - lambda) = 0, dl ascending and
// rho > 0, so that dl_j < lambda_j < dl_{j+1} (dl_K-1 + rho |z|^2 for the last). Returns lambda_j and
// delta_i = dl_i - lambda_j, computed against the nearer pole so that the small ones are accurate:
// the eigenvectors are built from them (dlaed4's approach). The iteration is the two-pole rational
// model of Bunch, Nielsen and Sorensen, kept inside a shrinking bracket by bisection, so it always
// terminates.
static double la_secular(size_t K, size_t j, const double* dl, const double* z, double rho, double* delta) {
    if (K == 1) {
        const double t = rho * z[0] * z[0];
        delta[0] = -t;
        return dl[0] + t;
    }
    const double rhoinv = 1.0 / rho;
    size_t org;
    double lo, hi;
    if (j + 1 < K) {
        const double mid = 0.5 * (dl[j + 1] - dl[j]);
        double w = rhoinv;
        for (size_t i = 0; i < K; ++i) w += z[i] * z[i] / ((dl[i] - dl[j]) - mid);
        if (w >= 0.0) { org = j;     lo = 0.0;  hi = mid; }
        else          { org = j + 1; lo = -mid; hi = 0.0; }
    } else {
        double zz = 0.0;
        for (size_t i = 0; i < K; ++i) zz += z[i] * z[i];
        org = j; lo = 0.0; hi = rho * zz;
    }
    // the model's two poles: the nearest on either side, or the last two for the last root
    const size_t pa = (j + 1 < K) ? j : K - 2;
    const size_t pb = pa + 1;
    const double dorg = dl[org];
    double tau = 0.5 * (lo + hi);

    for (int it = 0; it < 200; ++it) {
        double psi = 0.0, dpsi = 0.0, phi = 0.0, dphi = 0.0, err = 0.0;
        for (size_t i = 0; i < K; ++i) {
            const double di = (dl[i] - dorg) - tau;
            delta[i] = di;
            const double t = z[i] / di;
            const double term = z[i] * t;
            if (i <= pa) { psi += term; dpsi += t * t; }
            else         { phi += term; dphi += t * t; }
            err += fabs(term);
        }
        const double w = rhoinv + psi + phi;
        const double dw = dpsi + dphi;
        if (fabs(w) <= LA_EPS * (8.0 * (err + rhoinv) + fabs(tau) * dw)) break;
        if (w < 0.0) lo = tau; else hi = tau;
        if (!(hi - lo > 2.0 * LA_EPS * MD_MAX(fabs(lo), fabs(hi)))) break;

        const double Da = delta[pa], Db = delta[pb];
        const double C = w - Da * dpsi - Db * dphi;
        const double A = (Da + Db) * w - Da * Db * (dpsi + dphi);
        const double B = Da * Db * w;
        double eta;
        if (C == 0.0) {
            eta = B / A;
        } else {
            const double disc = sqrt(fabs(A * A - 4.0 * B * C));
            eta = (A <= 0.0) ? (A - disc) / (2.0 * C) : (2.0 * B) / (A + disc);
        }
        if (!(w * eta < 0.0)) eta = -w / dw;         // the model stepped the wrong way (or NaN): Newton
        double tn = tau + eta;
        if (!(tn > lo && tn < hi)) tn = 0.5 * (lo + hi);
        if (tn == tau) break;
        tau = tn;
    }
    for (size_t i = 0; i < K; ++i) delta[i] = (dl[i] - dorg) - tau;
    return dorg + tau;
}

typedef struct la_dc_t {
    const la_par_t* par;
    size_t    ldz;
    double*   e2;       // [n] leaf scratch
    double*   z;        // [n] the update vector
    double*   dl;       // [n] the undeflated poles
    double*   zz;       // [n] their weights
    double*   lam;      // [n] the new eigenvalues
    double*   gw;       // [n] Gu-Eisenstat products, then zhat
    double*   dfv;      // [n] the deflated eigenvalues
    double*   dv;       // [n] the values in the final order
    double*   col;      // [n] one column, for the permutation
    double*   mat;      // [2 n^2] the deflated columns, the gathered columns, the update's eigenvectors
    uint32_t* idx;      // [n] the sort of d
    uint32_t* nd;       // [n] undeflated (block columns)
    uint32_t* df;       // [n] deflated (block columns)
    uint32_t* ord;      // [n] the undeflated by column type, as positions into nd
    uint32_t* pos;      // [n] inverse of ord
    uint32_t* perm;     // [n] the final order
    uint8_t*  type;     // [n] 1: rows of the top half only, 3: bottom only, 2: both
    uint8_t*  seen;     // [n]
} la_dc_t;

// The eigenproblem of the update, D + rho z z^T over the K undeflated poles, in three parallel steps
typedef struct la_up_args_t {
    size_t K;
    double rho;
    const double* dl;
    const double* zz;
    double* U;              // K x K: the deltas, then the eigenvectors (rows in type order)
    double* lam;
    double* gw;
    const uint32_t* pos;
    const la_par_t* par;
} la_up_args_t;

static void la_up_roots(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const la_up_args_t* u = (const la_up_args_t*)arg;
    for (size_t j = beg; j < end; ++j) u->lam[j] = la_secular(u->K, j, u->dl, u->zz, u->rho, u->U + j * u->K);
}

static void la_up_weights(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const la_up_args_t* u = (const la_up_args_t*)arg;
    const size_t K = u->K;
    for (size_t i = beg; i < end; ++i) {
        double g = u->U[i * K + i];
        for (size_t j = 0; j < K; ++j) {
            if (j != i) g *= u->U[j * K + i] / (u->dl[i] - u->dl[j]);
        }
        u->gw[i] = copysign(sqrt(MD_MAX(0.0, -g)), u->zz[i]);
    }
}

static void la_up_vectors(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    const la_up_args_t* u = (const la_up_args_t*)arg;
    const size_t K = u->K;
    double* tmp = u->par->vec + (size_t)thread_idx * u->par->vec_size;
    for (size_t j = beg; j < end; ++j) {
        double* uj = u->U + j * K;
        double nrm = 0.0;
        for (size_t i = 0; i < K; ++i) { tmp[i] = u->gw[i] / uj[i]; nrm += tmp[i] * tmp[i]; }
        nrm = 1.0 / sqrt(nrm);
        for (size_t i = 0; i < K; ++i) uj[u->pos[i]] = tmp[i] * nrm;
    }
}

// Joins two solved halves of size k and m - k into the m x m block Z (ldz), eigenvalues d.
static bool la_dc_merge(double* d, double* Z, size_t m, size_t k, double beta, la_dc_t* w) {
    const size_t ldz = w->ldz;
    double rho = fabs(beta);
    const double sgn = beta < 0.0 ? -1.0 : 1.0;
    double* z = w->z;
    // z = Q^T v, v = e_{k-1} + sgn e_k: the last row of the top half, the first of the bottom one
    for (size_t j = 0; j < k; ++j) z[j] = Z[j * ldz + (k - 1)];
    for (size_t j = k; j < m; ++j) z[j] = sgn * Z[j * ldz + k];
    // |z|^2 = 2: normalize, and rho with it
    const double s2 = 1.0 / sqrt(2.0);
    for (size_t j = 0; j < m; ++j) z[j] *= s2;
    rho *= 2.0;

    // merge the two sorted halves
    uint32_t* idx = w->idx;
    {
        size_t a = 0, b = k, o = 0;
        while (a < k && b < m) idx[o++] = (uint32_t)(d[b] < d[a] ? b++ : a++);
        while (a < k) idx[o++] = (uint32_t)a++;
        while (b < m) idx[o++] = (uint32_t)b++;
    }
    double dmax = 0.0, zmax = 0.0;
    for (size_t j = 0; j < m; ++j) { dmax = MD_MAX(dmax, fabs(d[j])); zmax = MD_MAX(zmax, fabs(z[j])); }
    const double tol = 8.0 * LA_EPS * MD_MAX(dmax, zmax);

    uint8_t* type = w->type;
    for (size_t j = 0; j < m; ++j) type[j] = j < k ? 1 : 3;

    // Deflation (dlaed2): a negligible weight leaves its pair as it is; two close poles are rotated so
    // that one of them has no weight, and leaves too.
    uint32_t* nd = w->nd;
    uint32_t* df = w->df;
    size_t K = 0, ndf = 0;
    if (rho * zmax <= tol) {
        for (size_t t = 0; t < m; ++t) df[ndf++] = idx[t];
    } else {
        int64_t pj = -1;
        for (size_t t = 0; t < m; ++t) {
            const uint32_t nj = idx[t];
            if (rho * fabs(z[nj]) <= tol) { df[ndf++] = nj; continue; }
            if (pj < 0) { pj = nj; continue; }
            double s = z[pj], c = z[nj];
            const double tau = hypot(c, s);
            const double tt = d[nj] - d[pj];
            c /= tau;
            s = -s / tau;
            if (fabs(tt * c * s) <= tol) {
                z[nj] = tau;
                z[pj] = 0.0;
                if (type[nj] != type[pj]) type[nj] = type[pj] = 2;
                double* x = Z + (size_t)pj * ldz;
                double* y = Z + (size_t)nj * ldz;
                for (size_t r = 0; r < m; ++r) {
                    const double xr = x[r], yr = y[r];
                    x[r] = c * xr + s * yr;
                    y[r] = c * yr - s * xr;
                }
                const double t2 = d[pj] * c * c + d[nj] * s * s;
                d[nj] = d[pj] * s * s + d[nj] * c * c;
                d[pj] = t2;
                df[ndf++] = (uint32_t)pj;
            } else {
                nd[K++] = (uint32_t)pj;
            }
            pj = nj;
        }
        if (pj >= 0) nd[K++] = (uint32_t)pj;
    }
    // both lists by value (they are nearly so already; the rotations move values by at most tol)
    for (size_t i = 1; i < K; ++i) {
        const uint32_t v = nd[i];
        size_t j = i;
        while (j > 0 && d[nd[j - 1]] > d[v]) { nd[j] = nd[j - 1]; --j; }
        nd[j] = v;
    }
    for (size_t i = 1; i < ndf; ++i) {
        const uint32_t v = df[i];
        size_t j = i;
        while (j > 0 && d[df[j - 1]] > d[v]) { df[j] = df[j - 1]; --j; }
        df[j] = v;
    }

    // The deflated columns, saved before the block is overwritten
    double* DF = w->mat;
    for (size_t i = 0; i < ndf; ++i) {
        memcpy(DF + i * m, Z + (size_t)df[i] * ldz, sizeof(double) * m);
        w->dfv[i] = d[df[i]];
    }

    if (K > 0) {
        // the undeflated by type: rows [0,k) are nonzero only in types 1 and 2, rows [k,m) in 2 and 3
        uint32_t* ord = w->ord;
        uint32_t* pos = w->pos;
        size_t n1 = 0, n2 = 0, n3 = 0, o = 0;
        for (size_t i = 0; i < K; ++i) n1 += type[nd[i]] == 1, n2 += type[nd[i]] == 2, n3 += type[nd[i]] == 3;
        for (uint8_t ty = 1; ty <= 3; ++ty) {
            for (size_t i = 0; i < K; ++i) if (type[nd[i]] == ty) { pos[i] = (uint32_t)o; ord[o++] = (uint32_t)i; }
        }
        double* Qt = DF + ndf * m;                 // k x (n1 + n2)
        double* Qb = Qt + k * (n1 + n2);           // (m - k) x (n2 + n3)
        double* U  = Qb + (m - k) * (n2 + n3);     // K x K, rows in type order
        for (size_t q = 0; q < n1 + n2; ++q) memcpy(Qt + q * k, Z + (size_t)nd[ord[q]] * ldz, sizeof(double) * k);
        for (size_t q = n1; q < K; ++q) memcpy(Qb + (q - n1) * (m - k), Z + (size_t)nd[ord[q]] * ldz + k, sizeof(double) * (m - k));

        double* dl = w->dl;
        double* zz = w->zz;
        for (size_t i = 0; i < K; ++i) { dl[i] = d[nd[i]]; zz[i] = z[nd[i]]; }

        la_up_args_t up = { K, rho, dl, zz, U, w->lam, w->gw, pos, w->par };
        const la_par_t* par = w->par;
        // the roots, delta into U's columns (pole order for now)
        if (!md_parallel_for(par->ex, (uint32_t)K, 16, la_up_roots, &up)) return false;
        // Gu-Eisenstat: the weights for which the computed roots are exact
        if (!md_parallel_for(par->ex, (uint32_t)K, 64, la_up_weights, &up)) return false;
        // eigenvector j of the update: zhat_i / delta_ij, normalized, rows into type order
        if (!md_parallel_for(par->ex, (uint32_t)K, 16, la_up_vectors, &up)) return false;

        if (!la_gemm(par, false, k,     K, n1 + n2, 1.0, Qt, k,     U,      K, 0.0, Z,     ldz)) return false;
        if (!la_gemm(par, false, m - k, K, n2 + n3, 1.0, Qb, m - k, U + n1, K, 0.0, Z + k, ldz)) return false;
    }
    for (size_t i = 0; i < ndf; ++i) memcpy(Z + (K + i) * ldz, DF + i * m, sizeof(double) * m);

    // Into ascending order: two sorted runs, lam then the deflated, merged by a permutation of the
    // block's columns followed cycle by cycle
    double* dv = w->dv;
    uint32_t* perm = w->perm;      // perm[p] = the column that goes to p
    {
        size_t a = 0, b = 0, o = 0;
        while (a < K && b < ndf) perm[o++] = (uint32_t)(w->dfv[b] < w->lam[a] ? K + b++ : a++);
        while (a < K)   perm[o++] = (uint32_t)a++;
        while (b < ndf) perm[o++] = (uint32_t)(K + b++);
        for (size_t p = 0; p < m; ++p) dv[p] = perm[p] < K ? w->lam[perm[p]] : w->dfv[perm[p] - K];
    }
    memcpy(d, dv, sizeof(double) * m);
    uint8_t* seen = w->seen;
    memset(seen, 0, m);
    for (size_t s = 0; s < m; ++s) {
        if (seen[s] || perm[s] == s) { seen[s] = 1; continue; }
        memcpy(w->col, Z + s * ldz, sizeof(double) * m);
        size_t p = s;
        for (;;) {
            seen[p] = 1;
            const size_t src = perm[p];
            if (src == s) { memcpy(Z + p * ldz, w->col, sizeof(double) * m); break; }
            memcpy(Z + p * ldz, Z + src * ldz, sizeof(double) * m);
            p = src;
        }
    }
    return true;
}

static bool la_dc_rec(double* d, const double* e, double* Z, size_t m, la_dc_t* w) {
    const size_t ldz = w->ldz;
    if (m <= LA_LEAF) {
        for (size_t c = 0; c < m; ++c) {
            memset(Z + c * ldz, 0, sizeof(double) * m);
            Z[c * ldz + c] = 1.0;
        }
        if (m > 1) memcpy(w->e2, e, sizeof(double) * (m - 1));
        return la_tql2(d, w->e2, Z, m, ldz);
    }
    const size_t k = m / 2;
    const double beta = e[k - 1];
    d[k - 1] -= fabs(beta);
    d[k]     -= fabs(beta);
    if (!la_dc_rec(d, e, Z, k, w)) return false;
    if (!la_dc_rec(d + k, e + k, Z + k * ldz + k, m - k, w)) return false;
    return la_dc_merge(d, Z, m, k, beta, w);
}

bool md_linalg_sym_eigen_ex(double* a, double* eigval, size_t n, const struct md_executor_t* executor) {
    if (!a || !eigval || n == 0) return false;
    if (n == 1) {
        eigval[0] = a[0];
        a[0] = 1.0;
        return isfinite(eigval[0]);
    }

    const uint32_t nt = md_executor_num_threads(executor);
    const size_t nn = n * n;
    const size_t pack = la_pack_size(n);
    // one allocation: Z, the merge's matrices, the vectors, the per thread scratch, the index arrays
    const size_t num_dbl = nn + 2 * nn + 16 * n + LA_TRD_PIECES * n + LA_NB * LA_NB + 2 * LA_NB * n + (size_t)nt * (pack + n);
    const size_t bytes = sizeof(double) * num_dbl + sizeof(uint32_t) * 6 * n + 2 * n;
    md_allocator_i* heap = md_get_heap_allocator();
    double* mem = (double*)md_alloc(heap, bytes);
    if (!mem) return false;

    double* v = mem;
    double* Z = v; v += nn;
    la_par_t par = { executor, nt, NULL, pack, NULL, n };
    la_dc_t w = {0};
    w.par = &par;
    w.ldz = n;
    w.mat = v; v += 2 * nn;
    double* e     = v; v += n;
    double* tau   = v; v += n;
    double* ty    = v; v += 2 * n;
    double* tw    = v; v += 2 * n;
    double* ypart = v; v += LA_TRD_PIECES * n;
    w.e2  = v; v += n;
    w.z   = v; v += n;
    w.dl  = v; v += n;
    w.zz  = v; v += n;
    w.lam = v; v += n;
    w.gw  = v; v += n;
    w.dfv = v; v += n;
    w.dv  = v; v += n;
    w.col = v; v += n;
    v += n;   // spare, keeps the count at 16 n
    double* bT = v; v += LA_NB * LA_NB;
    double* bY = v; v += LA_NB * n;
    double* bW = v; v += LA_NB * n;
    par.pack = v; v += (size_t)nt * pack;
    par.vec  = v; v += (size_t)nt * n;
    uint32_t* iv = (uint32_t*)v;
    w.idx  = iv; iv += n;
    w.nd   = iv; iv += n;
    w.df   = iv; iv += n;
    w.ord  = iv; iv += n;
    w.pos  = iv; iv += n;
    w.perm = iv; iv += n;
    w.type = (uint8_t*)iv;
    w.seen = w.type + n;

    bool ok = false;
    if (!la_tridiag(a, n, eigval, e, tau, ty, tw, ypart, &par)) goto done;

    // scaled to norm 1, as dstedc does, so that the tolerances are relative
    // (and a NaN or Inf anywhere in the input has reached d or e by now: refused here, explicitly, as
    // a maximum would pass a NaN over)
    double nrm = 0.0;
    bool finite = true;
    for (size_t i = 0; i < n; ++i) {
        finite &= isfinite(eigval[i]) && isfinite(e[i]);
        nrm = MD_MAX(nrm, MD_MAX(fabs(eigval[i]), fabs(e[i])));
    }
    if (!finite) goto done;
    memset(Z, 0, sizeof(double) * nn);
    if (nrm == 0.0) {
        for (size_t i = 0; i < n; ++i) Z[i * n + i] = 1.0;
    } else {
        const double s = 1.0 / nrm;
        for (size_t i = 0; i < n; ++i) { eigval[i] *= s; e[i] *= s; }
        if (!la_dc_rec(eigval, e, Z, n, &w)) goto done;
        for (size_t i = 0; i < n; ++i) eigval[i] *= nrm;
    }
    if (!la_apply_q(a, n, tau, Z, bY, bT, bW, &par)) goto done;
    memcpy(a, Z, sizeof(double) * nn);
    ok = true;
done:
    md_free(heap, mem, bytes);
    return ok;
}

bool md_linalg_sym_eigen(double* a, double* eigval, size_t n) {
    return md_linalg_sym_eigen_ex(a, eigval, n, NULL);
}

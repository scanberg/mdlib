#include "utest.h"

#include <md_qm_ct.h>
#include <core/md_linalg.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Self contained: a random AO metric, orbitals orthonormal in it, and a transition built from
// random NTO pairs. Every quantity md_qm_ct computes from the NTOs is checked against the
// definition, evaluated the expensive way through an explicit S^(1/2) and the full A x A
// transition density.

enum { A = 40, M = 40, NOCC = 12, K = 5, G = 7 };

static unsigned lcg(unsigned* s) { *s = *s * 1664525u + 1013904223u; return *s >> 8; }
static double urand(unsigned* s) { return (double)lcg(s) / (double)(1u << 24) - 0.5; }

// f(S) for a symmetric S through its eigendecomposition: out = V diag(f(w)) V^T
static void sym_func(double* out, const double* S, size_t n, double power) {
    double* V = malloc(sizeof(double) * n * n);
    double* w = malloc(sizeof(double) * n);
    memcpy(V, S, sizeof(double) * n * n);
    md_linalg_sym_eigen(V, w, n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double acc = 0.0;
            for (size_t k = 0; k < n; ++k) acc += V[k * n + i] * pow(w[k], power) * V[k * n + j];
            out[i * n + j] = acc;
        }
    }
    free(V); free(w);
}

// A random orthogonal n x n matrix: the eigenvectors of a random symmetric one. Rows orthonormal.
static void random_orthogonal(double* Q, size_t n, unsigned seed) {
    double* w = malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j <= i; ++j) {
        Q[i * n + j] = Q[j * n + i] = urand(&seed);
    }
    md_linalg_sym_eigen(Q, w, n);
    free(w);
}

typedef struct problem_t {
    double S[A * A];
    double Sh[A * A];       // S^(1/2)
    double C[M * A];        // orbitals as rows, C^T S C = I in column terms
    double hole[K * A];     // NTOs as rows, AO basis
    double part[K * A];
    double lambda[K];
    uint32_t ao_group[A];
} problem_t;

static void make_problem(problem_t* p) {
    unsigned seed = 4711;
    // S = B^T B / A + I/2: symmetric positive definite, eigenvalues spread over a decade or so
    double* B = malloc(sizeof(double) * A * A);
    for (size_t i = 0; i < A * A; ++i) B[i] = urand(&seed) * 2.0;
    for (size_t i = 0; i < A; ++i) for (size_t j = 0; j < A; ++j) {
        double acc = 0.0;
        for (size_t k = 0; k < A; ++k) acc += B[k * A + i] * B[k * A + j];
        p->S[i * A + j] = acc / A + (i == j ? 0.5 : 0.0);
    }
    free(B);
    sym_func(p->Sh, p->S, A, 0.5);

    // C (columns) = S^(-1/2) Q, so C^T S C = Q^T Q = I. As rows: C_row[m] = (S^(-1/2) Q)[:, m]
    double* Sih = malloc(sizeof(double) * A * A);
    double* Q = malloc(sizeof(double) * A * A);
    sym_func(Sih, p->S, A, -0.5);
    random_orthogonal(Q, A, 99);
    for (size_t m = 0; m < M; ++m) for (size_t mu = 0; mu < A; ++mu) {
        double acc = 0.0;
        for (size_t nu = 0; nu < A; ++nu) acc += Sih[mu * A + nu] * Q[m * A + nu];   // Q's row m as column m
        p->C[m * A + mu] = acc;
    }

    // NTOs: orthonormal combinations of the occupied (hole) and the virtual (particle) orbitals
    double* U = malloc(sizeof(double) * NOCC * NOCC);
    double* V = malloc(sizeof(double) * (M - NOCC) * (M - NOCC));
    random_orthogonal(U, NOCC, 5);
    random_orthogonal(V, M - NOCC, 6);
    for (size_t k = 0; k < K; ++k) {
        for (size_t mu = 0; mu < A; ++mu) {
            double h = 0.0, q = 0.0;
            for (size_t i = 0; i < NOCC; ++i)     h += U[k * NOCC + i] * p->C[i * A + mu];
            for (size_t a = 0; a < M - NOCC; ++a) q += V[k * (M - NOCC) + a] * p->C[(NOCC + a) * A + mu];
            p->hole[k * A + mu] = h;
            p->part[k * A + mu] = q;
        }
    }
    // comparable weights, so the cross terms between pairs are large
    const double lam[K] = { 0.35, 0.25, 0.2, 0.15, 0.05 };
    memcpy(p->lambda, lam, sizeof(lam));

    for (size_t mu = 0; mu < A; ++mu) p->ao_group[mu] = (uint32_t)((mu * 3 + 1) % G);   // groups interleave

    free(Sih); free(Q); free(U); free(V);
}

// Omega the expensive way: T~ = S^(1/2) (sum_k s_k h_k p_k^T) S^(1/2), then sums of squares by group
static void omega_reference(double* out, const problem_t* p) {
    double* T = calloc(A * A, sizeof(double));
    double* Tt = calloc(A * A, sizeof(double));
    double* tmp = calloc(A * A, sizeof(double));
    for (size_t k = 0; k < K; ++k) {
        const double s = sqrt(p->lambda[k]);
        for (size_t i = 0; i < A; ++i) for (size_t j = 0; j < A; ++j) T[i * A + j] += s * p->hole[k * A + i] * p->part[k * A + j];
    }
    for (size_t i = 0; i < A; ++i) for (size_t j = 0; j < A; ++j) {
        double acc = 0.0;
        for (size_t k = 0; k < A; ++k) acc += p->Sh[i * A + k] * T[k * A + j];
        tmp[i * A + j] = acc;
    }
    for (size_t i = 0; i < A; ++i) for (size_t j = 0; j < A; ++j) {
        double acc = 0.0;
        for (size_t k = 0; k < A; ++k) acc += tmp[i * A + k] * p->Sh[k * A + j];
        Tt[i * A + j] = acc;
    }
    memset(out, 0, sizeof(double) * G * G);
    for (size_t i = 0; i < A; ++i) for (size_t j = 0; j < A; ++j) {
        out[p->ao_group[i] * G + p->ao_group[j]] += Tt[i * A + j] * Tt[i * A + j];
    }
    free(T); free(Tt); free(tmp);
}

UTEST(qm_ct, lowdin_orbitals_are_s_half_c) {
    problem_t* p = malloc(sizeof(problem_t));
    make_problem(p);

    double* Ct = malloc(sizeof(double) * M * A);
    ASSERT_TRUE(md_qm_lowdin_orbitals(Ct, p->C, M, A, NULL));

    double err = 0.0, orth = 0.0;
    for (size_t m = 0; m < M; ++m) {
        for (size_t mu = 0; mu < A; ++mu) {
            double ref = 0.0;
            for (size_t nu = 0; nu < A; ++nu) ref += p->Sh[mu * A + nu] * p->C[m * A + nu];
            err = fmax(err, fabs(Ct[m * A + mu] - ref));
        }
        for (size_t n = 0; n < M; ++n) {
            double d = 0.0;
            for (size_t mu = 0; mu < A; ++mu) d += Ct[m * A + mu] * Ct[n * A + mu];
            orth = fmax(orth, fabs(d - (m == n ? 1.0 : 0.0)));
        }
    }
    EXPECT_LT(err, 1.0e-11);
    EXPECT_LT(orth, 1.0e-11);

    // the transforms, both ways
    double* ht = malloc(sizeof(double) * K * A);
    double* back = malloc(sizeof(double) * K * A);
    md_qm_lowdin_transform(ht, p->hole, K, p->S, p->C, Ct, M, A);
    md_qm_lowdin_inverse_transform(back, ht, K, p->C, Ct, M, A);
    double terr = 0.0, berr = 0.0;
    for (size_t k = 0; k < K; ++k) for (size_t mu = 0; mu < A; ++mu) {
        double ref = 0.0;
        for (size_t nu = 0; nu < A; ++nu) ref += p->Sh[mu * A + nu] * p->hole[k * A + nu];
        terr = fmax(terr, fabs(ht[k * A + mu] - ref));
        berr = fmax(berr, fabs(back[k * A + mu] - p->hole[k * A + mu]));
    }
    EXPECT_LT(terr, 1.0e-11);
    EXPECT_LT(berr, 1.0e-11);

    // more orbitals than AOs cannot be independent
    EXPECT_FALSE(md_qm_lowdin_orbitals(Ct, p->C, A + 1, A, NULL));

    free(Ct); free(ht); free(back); free(p);
}

UTEST(qm_ct, omega_is_exact_with_interference) {
    problem_t* p = malloc(sizeof(problem_t));
    make_problem(p);

    double* Ct = malloc(sizeof(double) * M * A);
    double* ht = malloc(sizeof(double) * K * A);
    double* pt = malloc(sizeof(double) * K * A);
    ASSERT_TRUE(md_qm_lowdin_orbitals(Ct, p->C, M, A, NULL));
    md_qm_lowdin_transform(ht, p->hole, K, p->S, p->C, Ct, M, A);
    md_qm_lowdin_transform(pt, p->part, K, p->S, p->C, Ct, M, A);

    double Hb[G * K * K], Pb[G * K * K], omega[G * G], ref[G * G];
    md_qm_ct_group_blocks(Hb, ht, K, A, p->ao_group, G);
    md_qm_ct_group_blocks(Pb, pt, K, A, p->ao_group, G);
    md_qm_ct_omega(omega, p->lambda, K, Hb, G, Pb, G);
    omega_reference(ref, p);

    double err = 0.0, total = 0.0, incoherent_err = 0.0;
    for (size_t a = 0; a < G; ++a) for (size_t b = 0; b < G; ++b) {
        err = fmax(err, fabs(omega[a * G + b] - ref[a * G + b]));
        EXPECT_GE(omega[a * G + b], 0.0);
        total += omega[a * G + b];
        // one product per pair, the cross terms dropped
        double inc = 0.0;
        for (size_t k = 0; k < K; ++k) inc += p->lambda[k] * Hb[a * K * K + k * K + k] * Pb[b * K * K + k * K + k];
        incoherent_err = fmax(incoherent_err, fabs(inc - ref[a * G + b]));
    }
    EXPECT_LT(err, 1.0e-12);
    EXPECT_NEAR(total, 1.0, 1.0e-12);   // the lambdas sum to one, and the groups cover every AO
    // the test would prove nothing if the cross terms happened to be negligible here
    EXPECT_GT(incoherent_err, 1.0e-3);

    // the hole and particle per group are the row and column sums, and they DO decompose per pair
    for (size_t a = 0; a < G; ++a) {
        double row = 0.0, diag = 0.0;
        for (size_t b = 0; b < G; ++b) row += omega[a * G + b];
        for (size_t k = 0; k < K; ++k) diag += p->lambda[k] * Hb[a * K * K + k * K + k];
        EXPECT_NEAR(row, diag, 1.0e-12);
    }

    // blocks add: the block of groups {0,1} is the sum of theirs
    uint32_t merged[A];
    for (size_t mu = 0; mu < A; ++mu) merged[mu] = p->ao_group[mu] <= 1 ? 0 : p->ao_group[mu] - 1;
    double Hm[(G - 1) * K * K];
    md_qm_ct_group_blocks(Hm, ht, K, A, merged, G - 1);
    double add_err = 0.0;
    for (size_t i = 0; i < K * K; ++i) add_err = fmax(add_err, fabs(Hm[i] - (Hb[i] + Hb[K * K + i])));
    EXPECT_LT(add_err, 1.0e-14);

    free(Ct); free(ht); free(pt); free(p);
}

UTEST(qm_ct, block_nto) {
    problem_t* p = malloc(sizeof(problem_t));
    make_problem(p);

    double* Ct = malloc(sizeof(double) * M * A);
    double* ht = malloc(sizeof(double) * K * A);
    double* pt = malloc(sizeof(double) * K * A);
    ASSERT_TRUE(md_qm_lowdin_orbitals(Ct, p->C, M, A, NULL));
    md_qm_lowdin_transform(ht, p->hole, K, p->S, p->C, Ct, M, A);
    md_qm_lowdin_transform(pt, p->part, K, p->S, p->C, Ct, M, A);

    double Hb[G * K * K], Pb[G * K * K], omega[G * G];
    md_qm_ct_group_blocks(Hb, ht, K, A, p->ao_group, G);
    md_qm_ct_group_blocks(Pb, pt, K, A, p->ao_group, G);
    md_qm_ct_omega(omega, p->lambda, K, Hb, G, Pb, G);

    for (size_t a = 0; a < G; a += 2) {
        for (size_t b = 1; b < G; b += 3) {
            double w[K], hc[K * K], pc[K * K];
            const size_t n = md_qm_ct_block_nto(w, hc, pc, p->lambda, K, Hb + a * K * K, Pb + b * K * K);
            ASSERT_GT(n, 0u);
            double sum = 0.0;
            for (size_t j = 0; j < n; ++j) {
                sum += w[j];
                if (j > 0) EXPECT_LE(w[j], w[j - 1]);
            }
            EXPECT_NEAR(sum, omega[a * G + b], 1.0e-12);

            // u_j = sum_k hc[j][k] h~_k on group a, v_j likewise: unit norm, and u_j^T T~_ab v_j = sigma_j
            for (size_t j = 0; j < n; ++j) {
                double u[A] = {0}, v[A] = {0};
                for (size_t mu = 0; mu < A; ++mu) {
                    for (size_t k = 0; k < K; ++k) {
                        if (p->ao_group[mu] == a) u[mu] += hc[j * K + k] * ht[k * A + mu];
                        if (p->ao_group[mu] == b) v[mu] += pc[j * K + k] * pt[k * A + mu];
                    }
                }
                double nu = 0.0, nv = 0.0, uTv = 0.0;
                for (size_t mu = 0; mu < A; ++mu) { nu += u[mu] * u[mu]; nv += v[mu] * v[mu]; }
                for (size_t k = 0; k < K; ++k) {
                    double hu = 0.0, pv = 0.0;
                    for (size_t mu = 0; mu < A; ++mu) { hu += ht[k * A + mu] * u[mu]; pv += pt[k * A + mu] * v[mu]; }
                    uTv += sqrt(p->lambda[k]) * hu * pv;
                }
                EXPECT_NEAR(nu, 1.0, 1.0e-10);
                EXPECT_NEAR(nv, 1.0, 1.0e-10);
                EXPECT_NEAR(uTv, sqrt(w[j]), 1.0e-10);
            }
        }
    }
    free(Ct); free(ht); free(pt); free(p);
}

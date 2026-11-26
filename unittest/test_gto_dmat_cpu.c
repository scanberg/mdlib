// md_gto_dmat: a density matrix in tiles of the AO pairs that meet, in double on the CPU (no GPU).
//
//   gto_dmat_cpu.matches_file_matrix   h2o, amide and mol: the tiles built from the orbitals against the density
//                                      matrix VeloxChem wrote (its D, read raw from the file), every kept entry
//   gto_dmat_cpu.pattern_keeps_what_meets  mol, a sparse pattern (cutoff 1e-2): every AO pair whose shells meet is
//                                      kept and equal to the matrix; any other pair is 0 or kept in a tile with one
//   gto_dmat_cpu.signed_weights        a set of orbitals and its negation: D - D is 0 to rounding; -w is -D exactly
//   gto_dmat_cpu.executor_bitwise      the build and the evaluation on threads and shuffled: bitwise the calling
//                                      thread's
//   gto_dmat_cpu.evaluate_matches_matrix  mol: rho, its gradient and tau at points from the tiles against the whole
//                                      matrix with every AO (qm_test_ao_at), and ELF / LOL against md_qm_elf_evaluate
//                                      (at a negligible cutoff within its rounding; at 1e-7 within what the cutoff
//                                      leaves out)

#include "utest.h"
#include "vlx_test_util.h"
#include "test_executor.h"

#include <md_gto.h>
#include <md_gto_dmat.h>
#include <md_qm.h>
#include <core/md_allocator.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct dc_src_t {
    vlx_test_t     t;
    md_gto_basis_t basis;
    float*         xyz;             // bohr, per atom
    size_t         num_atoms, num_ao, num_mo;
    double*        rows;            // the SCF density's orbitals: num_mo x num_ao, the basis order
    double*        w;
    double*        D;               // the same as a matrix, num_ao x num_ao
} dc_src_t;

static bool dc_load(dc_src_t* s, const char* file) {
    memset(s, 0, sizeof(*s));
    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    if (!vlx_test_load(&s->t, md_str_from_cstr(path), MD_MEGABYTES(256))) return false;
    if (!qm_test_basis(&s->basis, &s->t)) return false;
    s->num_atoms = s->t.sys.atom.count;
    s->xyz = (float*)malloc(sizeof(float) * 3 * s->num_atoms);
    if (vlx_test_atom_xyz_bohr(s->xyz, 3 * s->num_atoms, &s->t) != s->num_atoms) return false;
    s->num_mo = md_qm_density_orbitals(&s->rows, &s->w, &s->num_ao, &s->t.sys, MD_STR_LIT("orbital/total/density"), md_get_heap_allocator());
    if (s->num_mo == 0 || s->num_ao != md_gto_basis_num_ao(&s->basis)) return false;
    s->D = (double*)malloc(sizeof(double) * s->num_ao * s->num_ao);
    return s->D && md_gto_density_from_orbitals(s->D, s->rows, s->num_ao, 1, s->w, s->num_mo, s->num_ao, NULL);
}

static void dc_free(dc_src_t* s) {
    if (s->rows) md_free(md_get_heap_allocator(), s->rows, sizeof(double) * s->num_mo * s->num_ao);
    if (s->w) md_free(md_get_heap_allocator(), s->w, sizeof(double) * s->num_mo);
    free(s->D);
    free(s->xyz);
    qm_test_free(&s->t);
}

static bool dc_init(md_gto_dmat_t* dm, const dc_src_t* s, double cutoff, double margin) {
    const md_gto_dmat_desc_t desc = { .basis = &s->basis, .atom_xyz = s->xyz, .cutoff = cutoff, .margin = margin };
    return md_gto_dmat_init(dm, &desc, md_get_heap_allocator());
}

UTEST(gto_dmat_cpu, matches_file_matrix) {
    static const char* files[] = { "h2o.h5", "amide.h5", "mol.h5" };
    for (int fi = 0; fi < 3; ++fi) {
        dc_src_t s;
        ASSERT_TRUE(dc_load(&s, files[fi]));
        size_t n = 0;
        // What VeloxChem wrote: its D, read from the file (symmetrised, in shell order and Cartesian)
        const double* F = qm_test_matrix(&s.t, MD_STR_LIT("orbital/total/density"), &n);
        ASSERT_TRUE(F != NULL);
        ASSERT_EQ(n, s.num_ao);
        md_gto_dmat_t dm;
        ASSERT_TRUE(dc_init(&dm, &s, 1.0e-7, 1.0));
        ASSERT_TRUE(md_gto_dmat_build(&dm, s.rows, s.num_ao, 1, s.w, s.num_mo, NULL));
        double d = 0.0, m = 0.0;
        size_t kept = 0;
        for (uint32_t i = 0; i < s.num_ao; ++i) {
            for (uint32_t j = 0; j < s.num_ao; ++j) {
                const double v = md_gto_dmat_at(&dm, i, j);
                m = fmax(m, fabs(F[(size_t)i * n + j]));
                if (v == 0.0) continue;
                ++kept;
                d = fmax(d, fabs(v - F[(size_t)i * n + j]));
            }
        }
        printf("  %s: %zu AOs, %llu tile pairs; the %zu entries kept against the file's D: max |diff| %.2e (max |D| %.3g)\n",
               files[fi], s.num_ao, (unsigned long long)dm.num_tiles, kept, d, m);
        EXPECT_LT(d, 1.0e-9 * m);
        md_gto_dmat_free(&dm);
        dc_free(&s);
    }
}

// The shells of AO i (by the basis' AO order) as an index, and whether shells meet as the pattern defines it
static void dc_ao_shells(uint32_t* ao_shell, const md_gto_basis_t* b) {
    uint32_t n = 0;
    for (uint32_t s = 0; s < b->num_shells; ++s) for (uint32_t c = 0; c < md_gto_num_cart_ao(b->shells[s].l); ++c) ao_shell[n++] = s;
}

UTEST(gto_dmat_cpu, pattern_keeps_what_meets) {
    dc_src_t s;
    ASSERT_TRUE(dc_load(&s, "mol.h5"));
    const double cutoff = 1.0e-2, margin = 0.0;
    md_gto_dmat_t dm;
    ASSERT_TRUE(dc_init(&dm, &s, cutoff, margin));
    ASSERT_TRUE(md_gto_dmat_build(&dm, s.rows, s.num_ao, 1, s.w, s.num_mo, NULL));
    const uint64_t total = (uint64_t)dm.nt * (dm.nt + 1) / 2;
    printf("  mol at %.0e: %llu of %llu tile pairs\n", cutoff, (unsigned long long)dm.num_tiles, (unsigned long long)total);
    EXPECT_LT(dm.num_tiles, total);     // sparse, or the test tests nothing
    uint32_t* ao_shell = (uint32_t*)malloc(sizeof(uint32_t) * s.num_ao);
    float* r = (float*)malloc(sizeof(float) * s.basis.num_shells);
    dc_ao_shells(ao_shell, &s.basis);
    md_gto_basis_shell_radii(r, &s.basis, cutoff);
    size_t meet = 0, wrong = 0;
    for (uint32_t i = 0; i < s.num_ao; ++i) {
        for (uint32_t j = 0; j < s.num_ao; ++j) {
            const uint32_t a = ao_shell[i], b = ao_shell[j];
            const float* p = s.xyz + 3 * s.basis.shells[a].atom_idx;
            const float* q = s.xyz + 3 * s.basis.shells[b].atom_idx;
            const double dx = p[0] - q[0], dy = p[1] - q[1], dz = p[2] - q[2];
            const double lim = (double)r[a] + (double)r[b] + margin;
            const bool meets = r[a] > 0.0f && r[b] > 0.0f && dx * dx + dy * dy + dz * dz <= lim * lim;
            const double v = md_gto_dmat_at(&dm, i, j), ref = s.D[(size_t)i * s.num_ao + j];
            meet += meets;
            // kept: the matrix to rounding (the tile holds its whole block); not kept: 0
            if (meets ? fabs(v - ref) > 1.0e-12 : (v != 0.0 && fabs(v - ref) > 1.0e-12)) ++wrong;
        }
    }
    printf("  %zu of %zu AO pairs meet; %zu wrong\n", meet, s.num_ao * s.num_ao, wrong);
    EXPECT_EQ(wrong, (size_t)0);
    free(ao_shell);
    free(r);
    md_gto_dmat_free(&dm);
    dc_free(&s);
}

UTEST(gto_dmat_cpu, signed_weights) {
    dc_src_t s;
    ASSERT_TRUE(dc_load(&s, "mol.h5"));
    md_gto_dmat_t a, b;
    ASSERT_TRUE(dc_init(&a, &s, 1.0e-7, 1.0));
    ASSERT_TRUE(dc_init(&b, &s, 1.0e-7, 1.0));
    // the orbitals twice, with w and -w: a difference of a density with itself
    double* rows2 = (double*)malloc(sizeof(double) * 2 * s.num_mo * s.num_ao);
    double* w2 = (double*)malloc(sizeof(double) * 2 * s.num_mo);
    memcpy(rows2, s.rows, sizeof(double) * s.num_mo * s.num_ao);
    memcpy(rows2 + s.num_mo * s.num_ao, s.rows, sizeof(double) * s.num_mo * s.num_ao);
    for (size_t k = 0; k < s.num_mo; ++k) { w2[k] = s.w[k]; w2[s.num_mo + k] = -s.w[k]; }
    ASSERT_TRUE(md_gto_dmat_build(&a, rows2, s.num_ao, 1, w2, 2 * s.num_mo, NULL));
    ASSERT_TRUE(md_gto_dmat_build(&b, s.rows, s.num_ao, 1, s.w, s.num_mo, NULL));
    double zero = 0.0, m = 0.0;
    for (uint64_t e = 0; e < 256 * a.num_tiles; ++e) { zero = fmax(zero, fabs(a.tiles[e])); m = fmax(m, fabs(b.tiles[e])); }
    printf("  D - D: max |entry| %.2e (max |D| %.3g)\n", zero, m);
    EXPECT_LT(zero, 1.0e-13 * m);
    for (size_t k = 0; k < s.num_mo; ++k) w2[k] = -s.w[k];
    ASSERT_TRUE(md_gto_dmat_build(&a, s.rows, s.num_ao, 1, w2, s.num_mo, NULL));
    bool neg = true;
    for (uint64_t e = 0; e < 256 * a.num_tiles; ++e) neg &= a.tiles[e] == -b.tiles[e];
    EXPECT_TRUE(neg);
    free(rows2);
    free(w2);
    md_gto_dmat_free(&a);
    md_gto_dmat_free(&b);
    dc_free(&s);
}

// Points around the molecule: around atoms taken in a stride, at a few distances, on a Fibonacci sphere each
static size_t dc_points(float* pts, size_t want, const dc_src_t* s) {
    static const float dist[] = { 0.3f, 0.9f, 1.6f, 2.5f, 4.0f };
    size_t n = 0;
    for (size_t k = 0; n < want; ++k) {
        const size_t a = (k * 7) % s->num_atoms;
        const float r = dist[k % 5];
        const double z = 1.0 - 2.0 * ((double)(k % 61) + 0.5) / 61.0, ph = 2.399963229728653 * (double)k;
        const double rr = sqrt(fmax(0.0, 1.0 - z * z));
        pts[3 * n + 0] = s->xyz[3 * a + 0] + r * (float)(rr * cos(ph));
        pts[3 * n + 1] = s->xyz[3 * a + 1] + r * (float)(rr * sin(ph));
        pts[3 * n + 2] = s->xyz[3 * a + 2] + r * (float)z;
        ++n;
    }
    return n;
}

UTEST(gto_dmat_cpu, executor_bitwise) {
    dc_src_t s;
    ASSERT_TRUE(dc_load(&s, "mol.h5"));
    enum { NP = 512 };
    float* pts = (float*)malloc(sizeof(float) * 3 * NP);
    dc_points(pts, NP, &s);
    test_executor_t ex[2] = { test_executor_threads(4), test_executor_shuffled(3) };
    md_gto_dmat_t ref;
    ASSERT_TRUE(dc_init(&ref, &s, 1.0e-7, 1.0));
    ASSERT_TRUE(md_gto_dmat_build(&ref, s.rows, s.num_ao, 1, s.w, s.num_mo, NULL));
    double* out[2][3];
    for (int k = 0; k < 2; ++k) {
        out[k][0] = (double*)malloc(sizeof(double) * NP);
        out[k][1] = (double*)malloc(sizeof(double) * 3 * NP);
        out[k][2] = (double*)malloc(sizeof(double) * NP);
    }
    ASSERT_EQ(md_gto_dmat_evaluate(out[0][0], out[0][1], out[0][2], &ref, &s.basis, s.xyz, 0, pts, NP, 0, NULL), (size_t)NP);
    for (int x = 0; x < 2; ++x) {
        md_gto_dmat_t dm;
        ASSERT_TRUE(dc_init(&dm, &s, 1.0e-7, 1.0));
        const md_executor_t* exe = test_executor(&ex[x]);
        ASSERT_TRUE(md_gto_dmat_build(&dm, s.rows, s.num_ao, 1, s.w, s.num_mo, exe));
        EXPECT_EQ(memcmp(dm.tiles, ref.tiles, sizeof(double) * 256 * ref.num_tiles), 0);
        ASSERT_EQ(md_gto_dmat_evaluate(out[1][0], out[1][1], out[1][2], &dm, &s.basis, s.xyz, 0, pts, NP, 0, exe), (size_t)NP);
        EXPECT_EQ(memcmp(out[0][0], out[1][0], sizeof(double) * NP), 0);
        EXPECT_EQ(memcmp(out[0][1], out[1][1], sizeof(double) * 3 * NP), 0);
        EXPECT_EQ(memcmp(out[0][2], out[1][2], sizeof(double) * NP), 0);
        md_gto_dmat_free(&dm);
    }
    for (int k = 0; k < 2; ++k) for (int j = 0; j < 3; ++j) free(out[k][j]);
    free(pts);
    md_gto_dmat_free(&ref);
    dc_free(&s);
}

UTEST(gto_dmat_cpu, evaluate_matches_matrix) {
    dc_src_t s;
    ASSERT_TRUE(dc_load(&s, "mol.h5"));
    enum { NP = 400 };
    float* pts = (float*)malloc(sizeof(float) * 3 * NP);
    dc_points(pts, NP, &s);
    const size_t n = s.num_ao;
    double* ao = (double*)malloc(sizeof(double) * n);
    double* aog = (double*)malloc(sizeof(double) * 3 * n);
    double* rho = (double*)malloc(sizeof(double) * NP);
    double* grad = (double*)malloc(sizeof(double) * 3 * NP);
    double* tau = (double*)malloc(sizeof(double) * NP);
    // a cutoff far below what matters, so that the screening leaves out nothing measurable: the tiles' own error
    md_gto_dmat_t dm;
    ASSERT_TRUE(dc_init(&dm, &s, 1.0e-14, 0.0));
    ASSERT_TRUE(md_gto_dmat_build(&dm, s.rows, n, 1, s.w, s.num_mo, NULL));
    ASSERT_EQ(md_gto_dmat_evaluate(rho, grad, tau, &dm, &s.basis, s.xyz, 0, pts, NP, 0, NULL), (size_t)NP);
    double e_rho = 0.0, e_g = 0.0, e_tau = 0.0, m_rho = 0.0, m_g = 0.0, m_tau = 0.0;
    for (size_t p = 0; p < NP; ++p) {
        const double x[3] = { pts[3 * p], pts[3 * p + 1], pts[3 * p + 2] };
        qm_test_ao_at(ao, aog, &s.basis, s.xyz, x);
        double r = 0.0, g[3] = { 0, 0, 0 }, t = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double u = 0.0, ug[3] = { 0, 0, 0 };
            for (size_t j = 0; j < n; ++j) {
                const double d = s.D[i * n + j];
                u += d * ao[j];
                for (int q = 0; q < 3; ++q) ug[q] += d * aog[3 * j + q];
            }
            r += ao[i] * u;
            for (int q = 0; q < 3; ++q) { g[q] += 2.0 * aog[3 * i + q] * u; t += 0.5 * aog[3 * i + q] * ug[q]; }
        }
        e_rho = fmax(e_rho, fabs(rho[p] - r));
        e_tau = fmax(e_tau, fabs(tau[p] - t));
        m_rho = fmax(m_rho, fabs(r));
        m_tau = fmax(m_tau, fabs(t));
        for (int q = 0; q < 3; ++q) { e_g = fmax(e_g, fabs(grad[3 * p + q] - g[q])); m_g = fmax(m_g, fabs(g[q])); }
    }
    printf("  mol, %d points, against the whole matrix: max |diff| rho %.2e (of %.3g), grad %.2e (of %.3g), tau %.2e (of %.3g)\n",
           NP, e_rho, m_rho, e_g, m_g, e_tau, m_tau);
    EXPECT_LT(e_rho, 1.0e-11 * m_rho);
    EXPECT_LT(e_g, 1.0e-11 * m_g);
    EXPECT_LT(e_tau, 1.0e-11 * m_tau);
    md_gto_dmat_free(&dm);

    // ELF and LOL from the tiles against the orbitals' (md_qm_elf_evaluate, every AO): at a negligible cutoff within
    // that reference's own rounding (its tau is over the derivative basis, whose coefficients a_p c_p are floats:
    // 6e-8 relative); at 1e-7, viamd's for ELF and LOL, within what that cutoff leaves out
    double* elf = (double*)malloc(sizeof(double) * NP);
    double* lol = (double*)malloc(sizeof(double) * NP);
    ASSERT_EQ(md_qm_elf_evaluate(elf, lol, pts, NP, 0, &s.t.sys, s.xyz, 0, 0.0, MD_QM_ELF_D_MIN), (size_t)NP);
    static const double cutoffs[2] = { 1.0e-14, 1.0e-7 }, tol_elf[2] = { 1.0e-6, 2.0e-5 }, tol_lol[2] = { 1.0e-7, 2.0e-6 };
    for (int c = 0; c < 2; ++c) {
        ASSERT_TRUE(dc_init(&dm, &s, cutoffs[c], 1.0));
        ASSERT_TRUE(md_gto_dmat_build(&dm, s.rows, n, 1, s.w, s.num_mo, NULL));
        ASSERT_EQ(md_gto_dmat_evaluate(rho, grad, tau, &dm, &s.basis, s.xyz, 0, pts, NP, 0, NULL), (size_t)NP);
        double e_elf = 0.0, e_lol = 0.0;
        size_t counted = 0;
        for (size_t p = 0; p < NP; ++p) {
            if (!(rho[p] > 1.0e-3)) continue;
            const double g2 = grad[3 * p] * grad[3 * p] + grad[3 * p + 1] * grad[3 * p + 1] + grad[3 * p + 2] * grad[3 * p + 2];
            e_elf = fmax(e_elf, fabs(md_qm_elf_closed_shell(rho[p], g2, tau[p], 0.0, MD_QM_ELF_D_MIN) - elf[p]));
            e_lol = fmax(e_lol, fabs(md_qm_lol_closed_shell(rho[p], tau[p], 0.0) - lol[p]));
            ++counted;
        }
        printf("  ELF and LOL at %.0e where rho > 1e-3 (%zu points): max |diff| %.2e, %.2e\n", cutoffs[c], counted, e_elf, e_lol);
        EXPECT_GT(counted, (size_t)NP / 4);
        EXPECT_LT(e_elf, tol_elf[c]);
        EXPECT_LT(e_lol, tol_lol[c]);
        md_gto_dmat_free(&dm);
    }
    free(elf);
    free(lol);
    free(ao); free(aog); free(rho); free(grad); free(tau); free(pts);
    dc_free(&s);
}

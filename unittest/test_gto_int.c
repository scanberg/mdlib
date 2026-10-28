#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include "qm_test_util.h"
#ifdef MD_HDF5
#include "vlx_test_util.h"
#endif

#include <md_gto.h>
#include <md_gto_int.h>
#include <md_molden.h>

#include <float.h>
#include <math.h>

// Tests for md_gto_int: the electrostatic potential of an AO density plus point charges.
//
// The reference is PySCF (int1e_grids), which shares no code with this library: the same molden
// file loaded by pyscf.tools.molden, its density D = C occ C^T over PySCF's own spherical basis,
// and V = sum_A Z_A/|C-A| - sum D (mu|1/|r-C||nu). Agreement therefore covers the reader, the
// spherical to Cartesian conversion of the coefficients, the AO convention and the integrals in
// one go.

#define INT_TEST_PI 3.14159265358979323846

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Everything the potential of a loaded QM system needs, read the way a consumer reads it.
typedef struct {
    md_gto_basis_t basis;
    float*         atom_xyz;    // bohr
    double*        atom_z;
    size_t         num_atoms;
    double*        D;
    size_t         num_ao;
} int_test_qm_t;

static bool int_test_extract(int_test_qm_t* out, const qm_test_t* t) {
    MEMSET(out, 0, sizeof(*out));
    if (!qm_test_basis(&out->basis, t)) return false;
    out->num_ao = md_gto_basis_num_ao(&out->basis);

    const size_t cap = 3 * 1024;
    out->atom_xyz = (float*)md_alloc(t->alloc, sizeof(float) * cap);
    {
        double* xyz = (double*)md_alloc(t->alloc, sizeof(double) * cap);
        const size_t nv = qm_test_series(xyz, cap, t, STR_LIT("qm/atom/coordinate"));   // Angstrom, 3 per atom
        for (size_t i = 0; i < nv; ++i) out->atom_xyz[i] = (float)(xyz[i] * QM_TEST_ANGSTROM_TO_BOHR);
        out->num_atoms = nv / 3;
    }
    if (out->num_atoms == 0) return false;

    // The charge the electrons were solved against, which is not the atomic number under an ECP
    out->atom_z = (double*)md_alloc(t->alloc, sizeof(double) * out->num_atoms);
    if (qm_test_series(out->atom_z, out->num_atoms, t, STR_LIT("qm/atom/nuclear_charge")) != out->num_atoms) return false;

    size_t dim = 0;
    out->D = qm_test_matrix(t, STR_LIT("orbital/total/density"), &dim);
    if (!out->D) out->D = qm_test_matrix(t, STR_LIT("orbital/alpha/density"), &dim);
    return out->D && dim == out->num_ao;
}

static bool int_test_charges(md_gto_int_charges_t* q, const int_test_qm_t* m, double threshold, bool nuclei, md_allocator_i* alloc) {
    md_gto_int_charges_desc_t desc = {
        .basis          = &m->basis,
        .atom_xyz       = m->atom_xyz,
        .density_matrix = m->D,
        .density_scale  = -1.0,
        .point_xyz      = nuclei ? m->atom_xyz : NULL,
        .point_charge   = nuclei ? m->atom_z : NULL,
        .num_points     = nuclei ? m->num_atoms : 0,
        .threshold      = threshold,
    };
    return md_gto_int_charges_init(q, &desc, alloc);
}

// ---------------------------------------------------------------------------
// Analytic cases
// ---------------------------------------------------------------------------

// A single normalised s function phi, density phi^2: a spherical Gaussian of exponent 2a and unit
// charge, whose potential is erf(sqrt(2a) r) / r. Two primitives so that the cross term (a merged
// gaussian of a different exponent) is exercised too, checked against the same closed form summed
// over the four products.
UTEST(gto_int, s_shell_closed_form) {
    float alpha[2] = { 1.3f, 0.35f };
    float coeff[2] = { 0.6f, 0.5f };
    md_gto_shell_t shell = { .atom_idx = 0, .primitive_offset = 0, .num_primitives = 2, .l = 0 };
    md_gto_basis_t basis = { .num_shells = 1, .num_primitives = 2, .shells = &shell, .alpha = alpha, .coeff = coeff };
    const float  atom[3] = { 0.25f, -0.5f, 1.0f };
    const double D[1] = { 2.0 };

    md_gto_int_charges_t q = {0};
    md_gto_int_charges_desc_t desc = { .basis = &basis, .atom_xyz = atom, .density_matrix = D, .density_scale = 1.0 };
    ASSERT_TRUE(md_gto_int_charges_init(&q, &desc, md_get_heap_allocator()));
    EXPECT_EQ(3u, q.num_gaussians);   // (a,a), (b,b) and the merged (a,b) = (b,a)
    EXPECT_EQ(0u, q.max_order);

    const float pts[4][3] = { {0.25f, -0.5f, 1.0f}, {1.0f, 0.0f, 1.5f}, {-2.0f, 3.0f, 0.0f}, {10.0f, 10.0f, -10.0f} };
    double V[4];
    md_gto_int_potential_xyz(V, NULL, &pts[0][0], 4, 0, &q);

    for (int k = 0; k < 4; ++k) {
        const double dx = pts[k][0] - atom[0], dy = pts[k][1] - atom[1], dz = pts[k][2] - atom[2];
        const double r = sqrt(dx * dx + dy * dy + dz * dz);
        double ref = 0.0;
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
            const double p = (double)alpha[i] + (double)alpha[j];
            const double charge = D[0] * coeff[i] * coeff[j] * pow(INT_TEST_PI / p, 1.5);
            ref += charge * (r < 1e-12 ? 2.0 * sqrt(p / INT_TEST_PI) : erf(sqrt(p) * r) / r);
        }
        EXPECT_NEAR(ref, V[k], 1e-12 * fmax(1.0, fabs(ref)));
    }

    const md_gto_int_moments_t m = md_gto_int_charges_moments(&q, NULL);
    double charge = 0.0;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        const double p = (double)alpha[i] + (double)alpha[j];
        charge += D[0] * coeff[i] * coeff[j] * pow(INT_TEST_PI / p, 1.5);
    }
    EXPECT_NEAR(charge, m.charge, 1e-12);
    // A spherical distribution's dipole about the coordinate origin is its charge times its centre,
    // and each Gaussian of exponent p adds q / 2p to the diagonal second moments.
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(charge * atom[k], m.dipole[k], 1e-12);
    double xx = 0.0;
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) {
        const double p = (double)alpha[i] + (double)alpha[j];
        xx += D[0] * coeff[i] * coeff[j] * pow(INT_TEST_PI / p, 1.5) * ((double)atom[0] * atom[0] + 0.5 / p);
    }
    EXPECT_NEAR(xx, m.second[0], 1e-12);
    EXPECT_NEAR(charge * atom[0] * atom[1], m.second[1], 1e-12);
    EXPECT_NEAR(charge * atom[1] * atom[2], m.second[4], 1e-12);

    md_gto_int_charges_free(&q, md_get_heap_allocator());
}

// Spherical gaussian charges, q (p/pi)^(3/2) exp(-p r^2): each a smeared point charge with the
// potential q erf(sqrt(p) r) / r, 2 q sqrt(p/pi) at its centre and q / r far out. Alone, and merged
// into a distribution that already has order 0 and higher gaussians from a density and point
// charges, where the potential has to be the sum of the parts - the merge shifts every gaussian
// and coefficient after the order 0 run.
UTEST(gto_int, gaussian_charges_closed_form) {
    md_allocator_i* heap = md_get_heap_allocator();
    const double gxyz[3][3] = { {0.0, 0.0, 0.0}, {1.2, -0.4, 0.3}, {-0.8, 2.0, -1.1} };
    const double gexp[3]    = { 86.9, 27.1, 3.0 };
    const double gq[3]      = { -1.2e-3, 4.9e-3, 0.7 };

    md_gto_int_charges_t alone = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&alone, &(md_gto_int_charges_desc_t){
        .gaussian_xyz = &gxyz[0][0], .gaussian_exponent = gexp, .gaussian_charge = gq, .num_gaussian_charges = 3 }, heap));
    EXPECT_EQ(3u, alone.num_gaussians);
    EXPECT_EQ(0u, alone.max_order);

    const float pts[5][3] = { {0, 0, 0}, {1.2f, -0.4f, 0.3f}, {0.05f, 0.02f, -0.03f}, {3.0f, 1.0f, -2.0f}, {40.0f, -30.0f, 20.0f} };
    double V[5], E[15];
    md_gto_int_potential_xyz(V, E, &pts[0][0], 5, 0, &alone);
    for (int i = 0; i < 5; ++i) {
        double ref = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double dx = pts[i][0] - gxyz[k][0], dy = pts[i][1] - gxyz[k][1], dz = pts[i][2] - gxyz[k][2];
            const double r = sqrt(dx * dx + dy * dy + dz * dz);
            ref += gq[k] * (r < 1e-7 ? 2.0 * sqrt(gexp[k] / INT_TEST_PI) : erf(sqrt(gexp[k]) * r) / r);
        }
        EXPECT_NEAR(ref, V[i], 1e-11 * fmax(1.0, fabs(ref)));
    }
    const md_gto_int_moments_t m = md_gto_int_charges_moments(&alone, NULL);
    EXPECT_NEAR(gq[0] + gq[1] + gq[2], m.charge, 1e-14);

    // Merged with an s and a p shell's density and two point charges
    float alpha[2] = { 1.3f, 0.35f };
    float coeff[2] = { 0.6f, 0.5f };
    md_gto_shell_t shells[2] = {
        { .atom_idx = 0, .primitive_offset = 0, .num_primitives = 2, .l = 0 },
        { .atom_idx = 1, .primitive_offset = 0, .num_primitives = 2, .l = 1 },
    };
    md_gto_basis_t basis = { .num_shells = 2, .num_primitives = 2, .shells = shells, .alpha = alpha, .coeff = coeff };
    const float atoms[2][3] = { {0.25f, -0.5f, 1.0f}, {-1.0f, 0.5f, 0.0f} };
    double D[16] = {0};
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) D[i * 4 + j] = (i == j) ? 1.0 : 0.1 * (i + j);
    const double z[2] = { 1.0, 2.0 };
    const md_gto_int_charges_desc_t base = { .basis = &basis, .atom_xyz = &atoms[0][0], .density_matrix = D,
                                             .point_xyz = &atoms[0][0], .point_charge = z, .num_points = 2 };
    md_gto_int_charges_desc_t both = base;
    both.gaussian_xyz = &gxyz[0][0];
    both.gaussian_exponent = gexp;
    both.gaussian_charge = gq;
    both.num_gaussian_charges = 3;

    md_gto_int_charges_t q_base = {0}, q_both = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&q_base, &base, heap));
    ASSERT_TRUE(md_gto_int_charges_init(&q_both, &both, heap));
    EXPECT_EQ(q_base.num_gaussians + 3, q_both.num_gaussians);
    EXPECT_EQ(q_base.max_order, q_both.max_order);
    EXPECT_TRUE(q_base.max_order >= 1);   // the merge has gaussians after the order 0 run to move

    double Vb[5], Eb[15], Vt[5], Et[15];
    md_gto_int_potential_xyz(Vb, Eb, &pts[0][0], 5, 0, &q_base);
    md_gto_int_potential_xyz(Vt, Et, &pts[0][0], 5, 0, &q_both);
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(Vb[i] + V[i], Vt[i], 1e-11 * fmax(1.0, fabs(Vt[i])));
        for (int c = 0; c < 3; ++c) EXPECT_NEAR(Eb[3 * i + c] + E[3 * i + c], Et[3 * i + c], 1e-10 * fmax(1.0, fabs(Et[3 * i + c])));
    }

    // And through the far field, which builds its tree over all the gaussians
    both.far_field_theta = MD_GTO_INT_FAR_FIELD_THETA;
    md_gto_int_charges_t q_far = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&q_far, &both, heap));
    double Vf[5];
    md_gto_int_potential_xyz(Vf, NULL, &pts[0][0], 5, 0, &q_far);
    for (int i = 0; i < 5; ++i) EXPECT_NEAR(Vt[i], Vf[i], 1e-4 * fmax(1.0, fabs(Vt[i])));

    // A non-positive exponent is no gaussian
    const double bad_exp[1] = { 0.0 };
    md_gto_int_charges_t q_bad = {0};
    EXPECT_FALSE(md_gto_int_charges_init(&q_bad, &(md_gto_int_charges_desc_t){
        .gaussian_xyz = &gxyz[0][0], .gaussian_exponent = bad_exp, .gaussian_charge = gq, .num_gaussian_charges = 1 }, heap));

    md_gto_int_charges_free(&alone, heap);
    md_gto_int_charges_free(&q_base, heap);
    md_gto_int_charges_free(&q_both, heap);
    md_gto_int_charges_free(&q_far, heap);
}

// Point charges only: Coulomb's law, and the field of each.
UTEST(gto_int, point_charges) {
    const float  xyz[2][3] = { {0, 0, 0}, {1.5f, 0, 0} };
    const double z[2] = { 1.0, -2.0 };
    md_gto_int_charges_t q = {0};
    md_gto_int_charges_desc_t desc = { .point_xyz = &xyz[0][0], .point_charge = z, .num_points = 2 };
    ASSERT_TRUE(md_gto_int_charges_init(&q, &desc, md_get_heap_allocator()));
    EXPECT_EQ(0u, q.num_gaussians);

    const float c[3] = { 0.5f, 1.0f, 0.0f };
    double V, E[3];
    md_gto_int_potential_xyz(&V, E, c, 1, 0, &q);
    const double r0 = sqrt(0.25 + 1.0), r1 = sqrt(1.0 + 1.0);
    EXPECT_NEAR(1.0 / r0 - 2.0 / r1, V, 1e-14);
    EXPECT_NEAR(0.5 / (r0 * r0 * r0) - 2.0 * (-1.0) / (r1 * r1 * r1), E[0], 1e-14);
    EXPECT_NEAR(1.0 / (r0 * r0 * r0) - 2.0 * 1.0 / (r1 * r1 * r1), E[1], 1e-14);

    md_gto_int_charges_free(&q, md_get_heap_allocator());
}

// A cluster of point charges a little around A, against the point multipoles it has there: its
// charge, dipole and second moment about A reproduce its potential and field up to the octupole,
// and its moments about any origin exactly. The cluster shares no code with the multipole terms, so
// this pins the convention - the signs, the 1/2, the trace - and not only self consistency.
UTEST(gto_int, point_multipoles_match_a_charge_cluster) {
    md_allocator_i* heap = md_get_heap_allocator();
    const float  A[3] = { 0.3f, -0.2f, 0.5f };
    const double h = 3.0e-2;
    const double dir[5][3] = { { 0.7, -0.1, 0.3 }, { -0.4, 0.8, -0.2 }, { 0.1, -0.6, -0.9 }, { -0.5, -0.3, 0.6 }, { 0.2, 0.5, 0.1 } };
    const double cq[5] = { 1.0, -0.6, 0.8, -1.4, 0.5 };

    // The cluster, and its moments about A from the positions as they are stored (float)
    float  pos[5][3];
    double charge = 0.0, mu[3] = {0}, Q[6] = {0};
    for (int i = 0; i < 5; ++i) {
        double sd[3];
        for (int k = 0; k < 3; ++k) {
            pos[i][k] = (float)(A[k] + h * dir[i][k]);
            sd[k] = (double)pos[i][k] - (double)A[k];
        }
        charge += cq[i];
        for (int k = 0; k < 3; ++k) mu[k] += cq[i] * sd[k];
        Q[0] += cq[i] * sd[0] * sd[0]; Q[1] += cq[i] * sd[0] * sd[1]; Q[2] += cq[i] * sd[0] * sd[2];
        Q[3] += cq[i] * sd[1] * sd[1]; Q[4] += cq[i] * sd[1] * sd[2]; Q[5] += cq[i] * sd[2] * sd[2];
    }

    md_gto_int_charges_t cluster = {0}, full = {0}, no_quad = {0}, no_dip = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&cluster, &(md_gto_int_charges_desc_t){ .point_xyz = &pos[0][0], .point_charge = cq, .num_points = 5 }, heap));
    ASSERT_TRUE(md_gto_int_charges_init(&full,    &(md_gto_int_charges_desc_t){ .point_xyz = A, .point_charge = &charge, .num_points = 1, .point_dipole = mu, .point_quadrupole = Q }, heap));
    ASSERT_TRUE(md_gto_int_charges_init(&no_quad, &(md_gto_int_charges_desc_t){ .point_xyz = A, .point_charge = &charge, .num_points = 1, .point_dipole = mu }, heap));
    ASSERT_TRUE(md_gto_int_charges_init(&no_dip,  &(md_gto_int_charges_desc_t){ .point_xyz = A, .point_charge = &charge, .num_points = 1, .point_quadrupole = Q }, heap));
    ASSERT_TRUE(full.point_dipole != NULL && full.point_quadrupole != NULL);
    EXPECT_TRUE(no_quad.point_quadrupole == NULL);
    EXPECT_TRUE(no_dip.point_dipole == NULL);

    // What remains is the octupole, ~ h/r of the quadrupole term, which is itself small. Summed over
    // the points, as at any one of them a term can be near a node of its angular dependence
    const float pts[6][3] = { { 3.3f, -0.2f, 0.5f }, { 0.3f, 3.0f, 0.9f }, { -2.1f, -1.9f, 1.7f }, { 1.4f, 0.8f, -2.6f }, { -0.9f, 2.2f, -1.8f }, { 2.0f, -2.5f, 2.4f } };
    double err = 0.0, errq = 0.0, errd = 0.0, eerr = 0.0, eerrq = 0.0;
    for (int i = 0; i < 6; ++i) {
        double Vc, Ec[3], Vf, Ef[3], Vq, Eq[3], Vd, Ed[3];
        md_gto_int_potential_xyz(&Vc, Ec, pts[i], 1, 0, &cluster);
        md_gto_int_potential_xyz(&Vf, Ef, pts[i], 1, 0, &full);
        md_gto_int_potential_xyz(&Vq, Eq, pts[i], 1, 0, &no_quad);
        md_gto_int_potential_xyz(&Vd, Ed, pts[i], 1, 0, &no_dip);
        EXPECT_LT(fabs(Vf - Vc), 1e-6);
        err  += fabs(Vf - Vc);
        errq += fabs(Vq - Vc);   // the quadrupole term, to within the octupole
        errd += fabs(Vd - Vc);   // the dipole term, likewise
        for (int k = 0; k < 3; ++k) {
            EXPECT_LT(fabs(Ef[k] - Ec[k]), 1e-6);
            eerr  += fabs(Ef[k] - Ec[k]);
            eerrq += fabs(Eq[k] - Ec[k]);
        }
    }
    EXPECT_LT(err,  0.02 * errq);
    EXPECT_LT(err,  0.001 * errd);
    EXPECT_LT(eerr, 0.05 * eerrq);

    // The field of the multipoles is minus their gradient
    const float step = 1.0f / 256.0f;
    for (int i = 0; i < 6; ++i) {
        double V, E[3];
        md_gto_int_potential_xyz(&V, E, pts[i], 1, 0, &full);
        for (int k = 0; k < 3; ++k) {
            float st[4][3];
            const float off[4] = { -2 * step, -step, step, 2 * step };
            for (int j = 0; j < 4; ++j) { MEMCPY(st[j], pts[i], sizeof(st[j])); st[j][k] += off[j]; }
            double Vs[4];
            md_gto_int_potential_xyz(Vs, NULL, &st[0][0], 4, 0, &full);
            const double dVdx = (Vs[0] - 8.0 * Vs[1] + 8.0 * Vs[2] - Vs[3]) / (12.0 * (double)step);
            EXPECT_NEAR(-dVdx, E[k], 1e-8);
        }
    }

    // Moments about an origin elsewhere: the multipoles carry the cluster's exactly
    const double O[3] = { -1.0, 2.0, 0.5 };
    const md_gto_int_moments_t mc = md_gto_int_charges_moments(&cluster, O);
    const md_gto_int_moments_t mf = md_gto_int_charges_moments(&full, O);
    EXPECT_NEAR(mc.charge, mf.charge, 1e-12);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(mc.dipole[k], mf.dipole[k], 1e-12);
    for (int k = 0; k < 6; ++k) EXPECT_NEAR(mc.second[k], mf.second[k], 1e-12);

    md_gto_int_charges_free(&cluster, heap);
    md_gto_int_charges_free(&full, heap);
    md_gto_int_charges_free(&no_quad, heap);
    md_gto_int_charges_free(&no_dip, heap);
}

// The trace of a quadrupole has no potential and no field: (1/2) Q:T with T traceless
UTEST(gto_int, point_quadrupole_trace_is_silent) {
    md_allocator_i* heap = md_get_heap_allocator();
    const float  A[3] = { 0, 0, 0 };
    const double zero = 0.0;
    const double iso[6] = { 2.5, 0, 0, 2.5, 0, 2.5 };
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&q, &(md_gto_int_charges_desc_t){ .point_xyz = A, .point_charge = &zero, .num_points = 1, .point_quadrupole = iso }, heap));
    const float c[3] = { 1.1f, -0.7f, 2.3f };
    double V, E[3];
    md_gto_int_potential_xyz(&V, E, c, 1, 0, &q);
    EXPECT_NEAR(0.0, V, 1e-14);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(0.0, E[k], 1e-14);
    md_gto_int_charges_free(&q, heap);
}

// ---------------------------------------------------------------------------
// Against PySCF
// ---------------------------------------------------------------------------

// Water, cc-pVDZ, from test_data/molden/h2o_ccpvdz.molden (spherical d). Points in bohr, potential
// in hartree/e, computed by PySCF 2.14 from the same file (see the comment at the top).
static const double h2o_ccpvdz_esp_ref[][4] = {
    {  9.748630, 13.606017, 14.258239,  9.632954183973e+00 },   // 0.4 bohr from O
    {  7.687530, 13.320320, 14.393436,  3.500177680931e+00 },   // 0.2 bohr from H
    { 10.648630, 13.106017, 14.958239,  1.914824635800e-01 },
    {  9.448630, 13.506017, 16.958239, -7.541033851893e-02 },
    {  6.448630, 14.506017, 14.958239,  5.162453965006e-02 },
    { 11.448630, 16.006017, 12.958239, -2.035581314166e-02 },
    {  9.448630,  9.506017, 16.458239, -1.602105971253e-02 },
    { 14.448630, 13.506017, 14.458239, -1.052795924192e-02 },
    {  3.448630,  7.506017, 17.458239,  3.527456977674e-03 },
    {  9.448630, 13.506017,  2.458239,  4.809274375542e-03 },
    { 29.448630, 23.506017,  9.458239, -5.785403846836e-04 },   // 20 bohr out
    {  8.667174, 12.533370, 11.865224,  8.632135811114e-02 },
};
// PySCF's dipole moment of the same density and nuclei, e bohr.
static const double h2o_ccpvdz_dipole_ref[3] = { -0.39479203722594036, -0.3250355388393018, -0.6656448852176311 };

UTEST(gto_int, h2o_molden_matches_pyscf) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(16));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));

    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    ASSERT_EQ(3u, m.num_atoms);

    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 0.0, true, t.alloc));

    const size_t n = ARRAY_SIZE(h2o_ccpvdz_esp_ref);
    float  pts[ARRAY_SIZE(h2o_ccpvdz_esp_ref)][3];
    double V[ARRAY_SIZE(h2o_ccpvdz_esp_ref)];
    for (size_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) pts[i][k] = (float)h2o_ccpvdz_esp_ref[i][k];
    md_gto_int_potential_xyz(V, NULL, &pts[0][0], n, 0, &q);

    // Not the integrals but the inputs limit the agreement: the basis stores exponents and
    // coefficients as float, and atoms and points are float positions, ~5e-7 bohr off at 10 bohr
    // from the origin. Near a nucleus the potential changes by Z/r^2 per bohr, so the first point,
    // 0.4 bohr from oxygen, is ~3e-5 au (3e-6 relative) away from PySCF for that reason alone.
    double worst = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double ref = h2o_ccpvdz_esp_ref[i][3];
        const double err = fabs(V[i] - ref);
        worst = MAX(worst, err / MAX(1.0, fabs(ref)));
        EXPECT_NEAR(ref, V[i], 2e-6 + 5e-6 * fabs(ref));
    }
    printf("h2o cc-pVDZ vs PySCF: max |dV| / max(1, |V|) = %.2e\n", worst);

    const md_gto_int_moments_t mom = md_gto_int_charges_moments(&q, NULL);
    EXPECT_NEAR(0.0, mom.charge, 2e-6);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(h2o_ccpvdz_dipole_ref[k], mom.dipole[k], 2e-5);

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}

// Screening is bounded: every dropped gaussian is below the threshold everywhere, so the error is
// at most the threshold times the number dropped, and in practice far less.
UTEST(gto_int, screening_error_is_bounded) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(16));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));

    md_gto_int_charges_t full = {0}, screened = {0};
    ASSERT_TRUE(int_test_charges(&full, &m, 0.0, true, t.alloc));
    const double threshold = 1e-6;
    ASSERT_TRUE(int_test_charges(&screened, &m, threshold, true, t.alloc));
    EXPECT_LT(screened.num_gaussians, full.num_gaussians);

    // The gaussians that survive carry identical coefficients: screening only removes.
    size_t dropped = full.num_gaussians - screened.num_gaussians;
    printf("h2o cc-pVDZ: %u gaussians, %u above %.0e\n", full.num_gaussians, screened.num_gaussians, threshold);

    const size_t n = ARRAY_SIZE(h2o_ccpvdz_esp_ref);
    float  pts[ARRAY_SIZE(h2o_ccpvdz_esp_ref)][3];
    double Va[ARRAY_SIZE(h2o_ccpvdz_esp_ref)], Vb[ARRAY_SIZE(h2o_ccpvdz_esp_ref)];
    for (size_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) pts[i][k] = (float)h2o_ccpvdz_esp_ref[i][k];
    md_gto_int_potential_xyz(Va, NULL, &pts[0][0], n, 0, &full);
    md_gto_int_potential_xyz(Vb, NULL, &pts[0][0], n, 0, &screened);
    for (size_t i = 0; i < n; ++i) EXPECT_LT(fabs(Va[i] - Vb[i]), threshold * (double)dropped);

    md_gto_int_charges_free(&full, t.alloc);
    md_gto_int_charges_free(&screened, t.alloc);
    qm_test_free(&t);
}

// Only the symmetric part of an AO matrix has a density, and that is all the build uses: adding an
// antisymmetric matrix changes nothing. (A transition density goes in as it is.)
UTEST(gto_int, antisymmetric_part_is_ignored) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(16));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));

    md_gto_int_charges_t a = {0}, b = {0};
    ASSERT_TRUE(int_test_charges(&a, &m, 0.0, false, t.alloc));
    const size_t N = m.num_ao;
    for (size_t i = 0; i < N; ++i) for (size_t j = 0; j < i; ++j) {
        const double x = 0.01 * (double)((i * 7 + j * 13) % 17) - 0.08;
        m.D[i * N + j] += x;
        m.D[j * N + i] -= x;
    }
    ASSERT_TRUE(int_test_charges(&b, &m, 0.0, false, t.alloc));
    EXPECT_EQ(a.num_gaussians, b.num_gaussians);

    const float pts[3][3] = { {9.0f, 13.0f, 14.0f}, {12.0f, 10.0f, 15.0f}, {0.0f, 0.0f, 0.0f} };
    double Va[3], Vb[3];
    md_gto_int_potential_xyz(Va, NULL, &pts[0][0], 3, 0, &a);
    md_gto_int_potential_xyz(Vb, NULL, &pts[0][0], 3, 0, &b);
    for (int i = 0; i < 3; ++i) EXPECT_NEAR(Va[i], Vb[i], 1e-12 * MAX(1.0, fabs(Va[i])));

    md_gto_int_charges_free(&a, t.alloc);
    md_gto_int_charges_free(&b, t.alloc);
    qm_test_free(&t);
}

// The field against central differences of the potential.
UTEST(gto_int, field_is_minus_gradient) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(16));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 0.0, true, t.alloc));

    const float h = 1.0f / 256.0f;   // exact in float, so the stencil points are exactly symmetric
    double worst = 0.0;
    for (size_t i = 1; i < ARRAY_SIZE(h2o_ccpvdz_esp_ref); ++i) {
        float c[3] = { (float)h2o_ccpvdz_esp_ref[i][0], (float)h2o_ccpvdz_esp_ref[i][1], (float)h2o_ccpvdz_esp_ref[i][2] };
        double V, E[3];
        md_gto_int_potential_xyz(&V, E, c, 1, 0, &q);
        for (int k = 0; k < 3; ++k) {
            float st[4][3];
            const float off[4] = { -2 * h, -h, h, 2 * h };
            for (int s = 0; s < 4; ++s) { MEMCPY(st[s], c, sizeof(c)); st[s][k] += off[s]; }
            double Vs[4];
            md_gto_int_potential_xyz(Vs, NULL, &st[0][0], 4, 0, &q);
            // fourth order central difference
            const double dVdx = (Vs[0] - 8.0 * Vs[1] + 8.0 * Vs[2] - Vs[3]) / (12.0 * (double)h);
            const double err = fabs(-dVdx - E[k]);
            worst = MAX(worst, err / MAX(1.0, fabs(E[k])));
            EXPECT_NEAR(-dVdx, E[k], 1e-6 * MAX(1.0, fabs(E[k])));
        }
    }
    printf("field vs -grad V (finite differences): max rel err %.2e\n", worst);

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}

// A cluster of copies of one molecule: n^3 of them on a lattice of the given spacing (bohr), each
// turned a little differently, with a block diagonal density. Big enough for a far field to matter.
typedef struct {
    md_gto_basis_t basis;
    float*  atom_xyz;
    double* atom_z;
    size_t  num_atoms;
    double* D;
    size_t  num_ao;
} int_test_cluster_t;

static void int_test_cluster(int_test_cluster_t* c, const int_test_qm_t* m, int n, float spacing, md_allocator_i* alloc) {
    const size_t copies = (size_t)n * n * n;
    MEMSET(c, 0, sizeof(*c));
    c->num_atoms = m->num_atoms * copies;
    c->num_ao    = m->num_ao * copies;
    c->atom_xyz  = (float*)md_alloc(alloc, sizeof(float) * 3 * c->num_atoms);
    c->atom_z    = (double*)md_alloc(alloc, sizeof(double) * c->num_atoms);
    c->D         = (double*)md_alloc(alloc, sizeof(double) * c->num_ao * c->num_ao);
    MEMSET(c->D, 0, sizeof(double) * c->num_ao * c->num_ao);
    c->basis.num_shells     = (uint32_t)(m->basis.num_shells * copies);
    c->basis.num_primitives = m->basis.num_primitives;
    c->basis.alpha          = m->basis.alpha;
    c->basis.coeff          = m->basis.coeff;
    c->basis.shells         = (md_gto_shell_t*)md_alloc(alloc, sizeof(md_gto_shell_t) * c->basis.num_shells);

    // the molecule about its first atom, so that a turn keeps it in its cell
    const float* o = m->atom_xyz;
    size_t ci = 0;
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x, ++ci) {
        const float a = 0.7f * (float)ci, ca = cosf(a), sa = sinf(a);
        for (size_t i = 0; i < m->num_atoms; ++i) {
            const float d[3] = { m->atom_xyz[3 * i] - o[0], m->atom_xyz[3 * i + 1] - o[1], m->atom_xyz[3 * i + 2] - o[2] };
            float* dst = c->atom_xyz + 3 * (ci * m->num_atoms + i);
            dst[0] = spacing * x + ca * d[0] - sa * d[1];
            dst[1] = spacing * y + sa * d[0] + ca * d[1];
            dst[2] = spacing * z + d[2];
            c->atom_z[ci * m->num_atoms + i] = m->atom_z[i];
        }
        for (uint32_t s = 0; s < m->basis.num_shells; ++s) {
            md_gto_shell_t sh = m->basis.shells[s];
            sh.atom_idx += (uint32_t)(ci * m->num_atoms);
            c->basis.shells[ci * m->basis.num_shells + s] = sh;
        }
        // A turn about z mixes the x and y AOs of each shell, so the copy's density is not the
        // original's in its AOs: the far field is exercised on a cluster of molecules with
        // densities that are not quite theirs. Harmless here - every check is against the exact
        // sum over the same distribution.
        const size_t a0 = ci * m->num_ao;
        for (size_t i = 0; i < m->num_ao; ++i) {
            MEMCPY(c->D + (a0 + i) * c->num_ao + a0, m->D + i * m->num_ao, sizeof(double) * m->num_ao);
        }
    }
}

// The far field against the exact sum on a cluster of 4^3 waters, at points among and around
// them: the potential and the field to the accuracy md_gto_int_charges_desc_t states, and exactly
// the same distribution otherwise (moments, gaussians).
UTEST(gto_int, far_field_matches_exact) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(256));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    int_test_cluster_t c;
    int_test_cluster(&c, &m, 4, 5.6f, t.alloc);

    md_gto_int_charges_desc_t desc = {
        .basis = &c.basis, .atom_xyz = c.atom_xyz, .density_matrix = c.D, .density_scale = -1.0,
        .point_xyz = c.atom_xyz, .point_charge = c.atom_z, .num_points = c.num_atoms, .threshold = 1e-8,
    };
    md_gto_int_charges_t exact = {0}, far = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&exact, &desc, t.alloc));
    desc.far_field_theta = MD_GTO_INT_FAR_FIELD_THETA;
    ASSERT_TRUE(md_gto_int_charges_init(&far, &desc, t.alloc));
    EXPECT_EQ(0u, exact.num_nodes);
    EXPECT_GT(far.num_nodes, 1u);
    EXPECT_EQ(exact.num_gaussians, far.num_gaussians);
    EXPECT_GE(far.far_field_order, far.max_order);

    // every member of the root, once
    {
        const uint32_t n = far.num_gaussians + far.num_points;
        uint8_t* seen = (uint8_t*)md_alloc(t.alloc, n);
        MEMSET(seen, 0, n);
        for (uint32_t i = 0; i < n; ++i) { ASSERT_LT(far.item[i], n); seen[far.item[i]]++; }
        for (uint32_t i = 0; i < n; ++i) EXPECT_EQ(1, seen[i]);
        EXPECT_EQ(0u, far.node[0].item_begin);
        EXPECT_EQ(n, far.node[0].item_end);
    }

    // Among the molecules (a lattice offset by half a cell), next to them, and far out
    enum { NP = 120 };
    float pts[NP][3];
    uint32_t seed = 99;
    for (int i = 0; i < NP; ++i) {
        for (int k = 0; k < 3; ++k) {
            seed = seed * 1664525u + 1013904223u;
            const float u = (float)(seed >> 8) / (float)(1u << 24);
            pts[i][k] = i < 80 ? -2.0f + u * 4.0f * 5.6f + 2.0f : (i < 100 ? -6.0f - u * 4.0f : 40.0f + 30.0f * u);
        }
    }
    double Ve[NP], Vf[NP], Ee[NP * 3], Ef[NP * 3];
    md_gto_int_potential_xyz(Ve, Ee, &pts[0][0], NP, 0, &exact);
    md_gto_int_potential_xyz(Vf, Ef, &pts[0][0], NP, 0, &far);
    double ev = 0.0, ee = 0.0;
    for (int i = 0; i < NP; ++i) {
        ev = MAX(ev, fabs(Ve[i] - Vf[i]));
        for (int k = 0; k < 3; ++k) ee = MAX(ee, fabs(Ee[3 * i + k] - Ef[3 * i + k]));
    }
    // and the potential alone, through the written out low orders
    double Vp[NP];
    md_gto_int_potential_xyz(Vp, NULL, &pts[0][0], NP, 0, &far);
    double ep = 0.0;
    for (int i = 0; i < NP; ++i) ep = MAX(ep, fabs(Vp[i] - Vf[i]));
    printf("64 waters, %u gaussians, %u nodes: far field max |dV| %.2e au, max |dE| %.2e au; work per point %llu of %llu\n",
           far.num_gaussians, far.num_nodes, ev, ee,
           (unsigned long long)md_gto_int_charges_work_per_point(&far), (unsigned long long)md_gto_int_charges_work_per_point(&exact));
    EXPECT_LT(ev, 1e-4);
    EXPECT_LT(ee, 1e-4);
    EXPECT_LT(ep, 1e-12);
    EXPECT_LT(md_gto_int_charges_work_per_point(&far), md_gto_int_charges_work_per_point(&exact));

    md_gto_int_charges_free(&exact, t.alloc);
    md_gto_int_charges_free(&far, t.alloc);
    qm_test_free(&t);
}

// The grid evaluator places its points exactly where the documentation says.
UTEST(gto_int, grid_matches_points) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(16));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 1e-10, true, t.alloc));

    // Rotated, anisotropic, not a multiple of anything.
    const float c = cosf(0.4f), s = sinf(0.4f);
    md_grid_t grid = {
        .orientation = { .elem = { {c, s, 0}, {-s, c, 0}, {0, 0, 1} } },
        .origin  = vec3_set(5.0f, 9.0f, 10.5f),
        .spacing = vec3_set(0.7f, 0.6f, 0.8f),
        .dim     = { 7, 9, 6 },
    };
    const float so[3] = { 0.5f, 0.5f, 0.5f };
    const size_t nv = md_grid_num_points(&grid);
    float* vol = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    md_gto_int_potential_grid(vol, &grid, so, &q);

    double worst = 0.0;
    for (int k = 0; k < grid.dim[2]; ++k) for (int j = 0; j < grid.dim[1]; ++j) for (int i = 0; i < grid.dim[0]; ++i) {
        const float l[3] = { (i + 0.5f) * 0.7f, (j + 0.5f) * 0.6f, (k + 0.5f) * 0.8f };
        const float p[3] = { 5.0f + c * l[0] - s * l[1], 9.0f + s * l[0] + c * l[1], 10.5f + l[2] };
        double V;
        md_gto_int_potential_xyz(&V, NULL, p, 1, 0, &q);
        const double v = vol[((size_t)k * grid.dim[1] + j) * grid.dim[0] + i];
        worst = MAX(worst, fabs(v - V) / MAX(1.0, fabs(V)));
    }
    EXPECT_LT(worst, 1e-5);   // float output, and float vs double point placement

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}

// The evaluators use a tabulated Boys function, and a written out potential for orders 0..2; a
// distribution without the table (boys_table NULL) is evaluated by the series and the general
// recursion. The two must agree to rounding for every order, and over every regime of T = p r^2:
// the Taylor expansion of the table, the asymptote with exp(-T), and the asymptote without it.
UTEST(gto_int, tabulated_boys_matches_series) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));

    // One normalised primitive per shell, l = 0..4 on each of three atoms: pairs up to g x g (L = 8)
    enum { NA = 3, NL = 5, NS = NA * NL };
    const float atom_xyz[NA][3] = { {0.0f, 0.0f, 0.0f}, {1.6f, 0.4f, -0.3f}, {-0.7f, 1.3f, 0.9f} };
    const float expo[NL] = { 9.5f, 1.3f, 0.9f, 0.7f, 0.15f };
    md_gto_shell_t shells[NS];
    float alpha[NS], coeff[NS];
    for (int a = 0; a < NA; ++a) for (int l = 0; l < NL; ++l) {
        const int s = a * NL + l;
        alpha[s] = expo[l] * (1.0f + 0.1f * a);
        coeff[s] = (float)md_qm_primitive_norm_factor((uint32_t)l, alpha[s]);
        shells[s] = (md_gto_shell_t){ .atom_idx = (uint32_t)a, .primitive_offset = (uint32_t)s, .num_primitives = 1, .l = (uint32_t)l };
    }
    md_gto_basis_t basis = { .num_shells = NS, .num_primitives = NS, .shells = shells, .alpha = alpha, .coeff = coeff };
    const size_t N = md_gto_basis_num_ao(&basis);
    double* D = (double*)md_alloc(alloc, sizeof(double) * N * N);
    uint32_t rng = 4711u;
    for (size_t i = 0; i < N; ++i) for (size_t j = 0; j <= i; ++j) {
        rng = rng * 1664525u + 1013904223u;
        const double v = ((double)(rng >> 8) / (double)(1u << 24) - 0.5) * (i == j ? 1.0 : 0.2);
        D[i * N + j] = D[j * N + i] = v;
    }
    const double Z[NA] = { 3.0, 1.0, 2.0 };

    md_gto_int_charges_t q = {0};
    md_gto_int_charges_desc_t desc = {
        .basis = &basis, .atom_xyz = &atom_xyz[0][0], .density_matrix = D, .density_scale = -1.0,
        .point_xyz = &atom_xyz[0][0], .point_charge = Z, .num_points = NA,
    };
    ASSERT_TRUE(md_gto_int_charges_init(&q, &desc, alloc));
    ASSERT_TRUE(q.boys_table != NULL);
    EXPECT_EQ(8u, q.max_order);

    // Points on rays out from near the first atom to 60 bohr, geometrically spaced: with exponents
    // from 0.3 to 21 that puts T anywhere from ~0 to ~1e5
    enum { NR = 7, NSTEP = 40, NP = NR * NSTEP };
    float pts[NP][3];
    for (int r = 0; r < NR; ++r) {
        const float dir[3] = { cosf(0.9f * r) * sinf(0.5f + 0.4f * r), sinf(0.9f * r) * sinf(0.5f + 0.4f * r), cosf(0.5f + 0.4f * r) };
        for (int s = 0; s < NSTEP; ++s) {
            const float d = 0.02f * powf(3000.0f, (float)s / (NSTEP - 1));
            for (int k = 0; k < 3; ++k) pts[r * NSTEP + s][k] = 0.05f + d * dir[k];
        }
    }

    double V_tab[NP], E_tab[NP * 3], V_ser[NP], E_ser[NP * 3], V_tab_only[NP];
    md_gto_int_potential_xyz(V_tab, E_tab, &pts[0][0], NP, 0, &q);
    md_gto_int_potential_xyz(V_tab_only, NULL, &pts[0][0], NP, 0, &q);   // the written out orders
    double* table = q.boys_table;
    q.boys_table = NULL;
    md_gto_int_potential_xyz(V_ser, E_ser, &pts[0][0], NP, 0, &q);
    q.boys_table = table;

    double worst_V = 0.0, worst_E = 0.0;
    for (int i = 0; i < NP; ++i) {
        worst_V = MAX(worst_V, fabs(V_tab[i] - V_ser[i]) / MAX(1.0, fabs(V_ser[i])));
        worst_V = MAX(worst_V, fabs(V_tab_only[i] - V_ser[i]) / MAX(1.0, fabs(V_ser[i])));
        for (int k = 0; k < 3; ++k) worst_E = MAX(worst_E, fabs(E_tab[3 * i + k] - E_ser[3 * i + k]) / MAX(1.0, fabs(E_ser[3 * i + k])));
    }
    printf("tabulated vs series Boys, s..g: max rel diff V %.2e, E %.2e\n", worst_V, worst_E);
    EXPECT_LT(worst_V, 1e-13);
    EXPECT_LT(worst_E, 1e-13);

    md_arena_allocator_destroy(alloc);
}

#ifdef MD_HDF5
// VeloxChem water: no external potential to compare with, but the reader publishes the ground state
// dipole VeloxChem computed, and the far field of a correct density and nuclei IS that dipole. This
// covers the spherical to Cartesian conversion of the density (T^T D T) end to end.
UTEST(gto_int, vlx_h2o_dipole_matches_published) {
    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/h2o.h5"), MEGABYTES(32)));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));

    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 0.0, true, t.alloc));
    const md_gto_int_moments_t mom = md_gto_int_charges_moments(&q, NULL);
    EXPECT_NEAR(0.0, mom.charge, 1e-5);

    const md_attribute_t* a = md_attributes_find(&t.sys.attributes, STR_LIT("dipole/ground_state/vector"));
    if (!a) {
        UTEST_SKIP("h2o.h5 carries no ground state dipole");
    }
    double dip[3];
    ASSERT_EQ(3u, md_attribute_extract_f64(dip, 3, a, md_attribute_slice_all(), md_unit_elementary_charge_bohr()));
    printf("vlx h2o dipole: published (%.6f %.6f %.6f), from the density (%.6f %.6f %.6f) e bohr\n",
           dip[0], dip[1], dip[2], mom.dipole[0], mom.dipole[1], mom.dipole[2]);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(dip[k], mom.dipole[k], 1e-4);

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
}
#endif

// ---------------------------------------------------------------------------
// GPU against the CPU reference
// ---------------------------------------------------------------------------
#if MD_ENABLE_GPU

static const char* int_test_no_device_reason(void) {
    const char* err = md_gpu_last_error();
    return (err && err[0]) ? err : "No GPU device available";
}

// Evaluates `q` on `grid` on the GPU and reads the result back into out.
static bool int_test_gpu_grid(float* out, md_gpu_stream_t stream, const md_gto_int_charges_t* q, const md_grid_t* grid,
                              const float so[3], uint64_t work_per_dispatch) {
    const size_t nv = md_grid_num_points(grid);
    md_gpu_texture_t tex = md_gpu_texture_create(stream, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid->dim[0], .height = (uint32_t)grid->dim[1], .depth_or_layers = (uint32_t)grid->dim[2],
    });
    md_gpu_mem_t rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * nv);
    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, q);
    if (!tex || !rb.cpu || !gq) return false;

    md_gto_int_gpu_potential_desc_t desc = {
        .charges = gq, .out_tex = tex, .grid = grid,
        .sample_offset = { so[0], so[1], so[2] }, .op = MD_GTO_OP_SET,
        .work_per_dispatch = work_per_dispatch,
    };
    if (!md_gto_int_gpu_potential_launch(stream, &desc)) return false;
    md_gpu_copy_from_texture(stream, rb.gpu, tex, NULL);
    md_gpu_stream_sync(stream);
    MEMCPY(out, rb.cpu, sizeof(float) * nv);

    md_gto_int_gpu_charges_destroy(stream, gq);
    md_gpu_free(stream, rb.gpu);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_sync(stream);
    return true;
}

// max |gpu - cpu| / max(1, |cpu|): the potential spans ~10 au near a nucleus to ~1e-3 far out,
// and float terms are good to a few parts in 1e7 of the larger of the two.
static double int_test_compare(const float* gpu, const float* cpu, size_t n) {
    double worst = 0.0;
    for (size_t i = 0; i < n; ++i) {
        uint32_t u;
        MEMCPY(&u, &gpu[i], sizeof(u));
        if ((u & 0x7F800000u) == 0x7F800000u) return DBL_MAX;   // not finite (fast-math safe test)
        worst = MAX(worst, fabs((double)gpu[i] - (double)cpu[i]) / MAX(1.0, fabs((double)cpu[i])));
    }
    return worst;
}

UTEST(gto_int, h2o_gpu_matches_cpu) {
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        UTEST_SKIP(int_test_no_device_reason());
    }
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);

    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(32));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 1e-10, true, t.alloc));

    // Around the molecule, not a multiple of the 4x4x4 group on any axis.
    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin  = vec3_set(4.0f, 8.5f, 8.0f),
        .spacing = vec3_set(0.37f, 0.41f, 0.43f),
        .dim     = { 30, 27, 31 },
    };
    const float so[3] = { 0.5f, 0.5f, 0.5f };
    const size_t nv = md_grid_num_points(&grid);
    float* cpu = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    float* gpu = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    float* gpu_split = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    md_gto_int_potential_grid(cpu, &grid, so, &q);

    ASSERT_TRUE(int_test_gpu_grid(gpu, stream, &q, &grid, so, 0));
    const double rel = int_test_compare(gpu, cpu, nv);
    // Within a bohr of a nucleus the potential changes by up to ~100 au per bohr, and the float
    // voxel position (~2e-7 bohr) is what limits the agreement there. Further out the terms are.
    double rel_far = 0.0;
    for (int k = 0; k < grid.dim[2]; ++k) for (int j = 0; j < grid.dim[1]; ++j) for (int i = 0; i < grid.dim[0]; ++i) {
        const double c[3] = { 4.0 + (i + 0.5) * 0.37f, 8.5 + (j + 0.5) * 0.41f, 8.0 + (k + 0.5) * 0.43f };
        double dmin = HUGE_VAL;
        for (size_t a = 0; a < m.num_atoms; ++a) {
            const double dx = c[0] - m.atom_xyz[a * 3 + 0], dy = c[1] - m.atom_xyz[a * 3 + 1], dz = c[2] - m.atom_xyz[a * 3 + 2];
            dmin = MIN(dmin, sqrt(dx * dx + dy * dy + dz * dz));
        }
        if (dmin < 1.0) continue;
        const size_t v = ((size_t)k * grid.dim[1] + j) * grid.dim[0] + i;
        rel_far = MAX(rel_far, fabs((double)gpu[v] - (double)cpu[v]) / MAX(1.0, fabs((double)cpu[v])));
    }
    printf("h2o cc-pVDZ potential, GPU vs CPU: max |diff| / max(1,|V|) = %.3e, beyond 1 bohr of a nucleus %.3e\n", rel, rel_far);
    EXPECT_LT(rel, 5e-6);
    EXPECT_LT(rel_far, 2e-6);

    // What is left once the position is accounted for: |diff| - |E| * 1e-6 bohr, with E the field.
    // That is the error of the float terms themselves, and it is small everywhere.
    {
        float*  pts = (float*)md_alloc(t.alloc, sizeof(float) * 3 * nv);
        double* V   = (double*)md_alloc(t.alloc, sizeof(double) * nv);
        double* E   = (double*)md_alloc(t.alloc, sizeof(double) * 3 * nv);
        for (int k = 0; k < grid.dim[2]; ++k) for (int j = 0; j < grid.dim[1]; ++j) for (int i = 0; i < grid.dim[0]; ++i) {
            const size_t v = ((size_t)k * grid.dim[1] + j) * grid.dim[0] + i;
            pts[v * 3 + 0] = 4.0f + (i + 0.5f) * 0.37f;
            pts[v * 3 + 1] = 8.5f + (j + 0.5f) * 0.41f;
            pts[v * 3 + 2] = 8.0f + (k + 0.5f) * 0.43f;
        }
        md_gto_int_potential_xyz(V, E, pts, nv, 0, &q);
        double excess = 0.0;
        for (size_t v = 0; v < nv; ++v) {
            const double e = sqrt(E[v * 3] * E[v * 3] + E[v * 3 + 1] * E[v * 3 + 1] + E[v * 3 + 2] * E[v * 3 + 2]);
            excess = MAX(excess, (fabs((double)gpu[v] - V[v]) - e * 1e-6) / MAX(1.0, fabs(V[v])));
        }
        printf("  beyond what 1e-6 bohr of position explains: %.3e\n", excess);
        EXPECT_LT(excess, 1e-6);
    }

    // Split into many slabs and runs of gaussians: the fixed point sum makes it bit identical.
    ASSERT_TRUE(int_test_gpu_grid(gpu_split, stream, &q, &grid, so, 20000));
    EXPECT_EQ(0, memcmp(gpu, gpu_split, sizeof(float) * nv));

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
}

// Evaluates `q` at the points `xyz` (stride bytes apart, float) on the GPU and reads the result back.
static bool int_test_gpu_xyz(float* out, md_gpu_stream_t stream, const md_gto_int_charges_t* q, const float* xyz, size_t n,
                             uint32_t stride, uint64_t work_per_dispatch) {
    const size_t bytes = (size_t)(stride ? stride : 12) * n;
    md_gpu_mem_t pts = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, bytes);
    md_gpu_mem_t res = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * n);
    md_gpu_mem_t rb  = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * n);
    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, q);
    if (!pts.gpu || !res.gpu || !rb.cpu || !gq) return false;
    void* dst = md_gpu_upload_begin(stream, pts.gpu, bytes);
    MEMCPY(dst, xyz, bytes);
    md_gpu_upload_end(stream);

    md_gto_int_gpu_potential_xyz_desc_t desc = {
        .charges = gq, .xyz = pts.gpu, .xyz_stride = stride, .num_xyz = (uint32_t)n,
        .out_potential = res.gpu, .work_per_dispatch = work_per_dispatch,
    };
    if (!md_gto_int_gpu_potential_xyz_launch(stream, &desc)) return false;
    md_gpu_copy(stream, rb.gpu, res.gpu, sizeof(float) * n);
    md_gpu_stream_sync(stream);
    MEMCPY(out, rb.cpu, sizeof(float) * n);

    md_gto_int_gpu_charges_destroy(stream, gq);
    md_gpu_free(stream, pts.gpu);
    md_gpu_free(stream, res.gpu);
    md_gpu_free(stream, rb.gpu);
    md_gpu_stream_sync(stream);
    return true;
}

#ifdef MD_HDF5
// md_gto_int_gpu_potential_xyz_launch against md_gto_int_potential_xyz at scattered points, packed and
// with the 16 byte stride of a float4 vertex buffer (md_gto_iso_gpu_mesh_t::position)
UTEST(gto_int, vlx_amide_gpu_xyz_matches_cpu) {
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        UTEST_SKIP(int_test_no_device_reason());
    }
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);

    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.h5"), MEGABYTES(64)));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    md_gto_int_charges_t q = {0};
    ASSERT_TRUE(int_test_charges(&q, &m, 1e-10, true, t.alloc));

    // Points 0.5 - 4 bohr from the atoms, in every direction: where isosurfaces are
    const size_t n = 10007;     // not a multiple of the group size
    float*  p3  = (float*)md_alloc(t.alloc, sizeof(float) * 3 * n);
    float*  p4  = (float*)md_alloc(t.alloc, sizeof(float) * 4 * n);
    double* V   = (double*)md_alloc(t.alloc, sizeof(double) * n);
    double* E   = (double*)md_alloc(t.alloc, sizeof(double) * 3 * n);
    float*  gpu = (float*)md_alloc(t.alloc, sizeof(float) * n);
    float*  gpu4 = (float*)md_alloc(t.alloc, sizeof(float) * n);
    float*  gpu_split = (float*)md_alloc(t.alloc, sizeof(float) * n);
    uint32_t rng = 12345u;
    for (size_t i = 0; i < n; ++i) {
        const size_t a = i % m.num_atoms;
        float d[3], len2 = 0.0f;
        do {
            len2 = 0.0f;
            for (int k = 0; k < 3; ++k) {
                rng = rng * 1664525u + 1013904223u;
                d[k] = (float)(rng >> 8) * (2.0f / 16777216.0f) - 1.0f;
                len2 += d[k] * d[k];
            }
        } while (len2 < 1e-4f || len2 > 1.0f);
        rng = rng * 1664525u + 1013904223u;
        const float r = 0.5f + 3.5f * (float)(rng >> 8) / 16777216.0f;
        for (int k = 0; k < 3; ++k) {
            p3[i * 3 + k] = (float)m.atom_xyz[a * 3 + k] + r * d[k] / sqrtf(len2);
            p4[i * 4 + k] = p3[i * 3 + k];
        }
        p4[i * 4 + 3] = 1.0f;
    }
    md_gto_int_potential_xyz(V, E, p3, n, 0, &q);

    ASSERT_TRUE(int_test_gpu_xyz(gpu, stream, &q, p3, n, 0, 0));
    ASSERT_TRUE(int_test_gpu_xyz(gpu4, stream, &q, p4, n, 16, 0));
    EXPECT_EQ(0, memcmp(gpu, gpu4, sizeof(float) * n));

    // As for the grid: what is left beyond what the position explains. The points are absolute floats,
    // moved into the frame of the distribution in float on the device (the grid's are formed in double):
    // a rounding of their own, so 2e-6 bohr rather than the grid's 1e-6.
    // What is left is the rounding of the float terms, which goes with the size of the terms, not of V:
    // a few bohr from the molecule the nuclei's potential and the electrons' cancel to a fraction of
    // either (V 0.16 from 12.9 and -12.7), so it is measured against |V_nuc| + |V_el|, the size of what
    // cancels, in float epsilons. Backends round differently (Metal compiles with fast math, Vulkan
    // allows its exp, rsqrt and division a few ulps): a couple of epsilons, not a fraction of one.
    double rel = 0.0, excess = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double e = sqrt(E[i * 3] * E[i * 3] + E[i * 3 + 1] * E[i * 3 + 1] + E[i * 3 + 2] * E[i * 3 + 2]);
        const double d = fabs((double)gpu[i] - V[i]);
        double vn = 0.0;
        for (uint32_t k = 0; k < q.num_points; ++k) {
            const double dx = p3[i * 3 + 0] - q.point_xyz[k * 3 + 0], dy = p3[i * 3 + 1] - q.point_xyz[k * 3 + 1], dz = p3[i * 3 + 2] - q.point_xyz[k * 3 + 2];
            vn += q.point_charge[k] / sqrt(dx * dx + dy * dy + dz * dz);
        }
        rel    = MAX(rel, d / MAX(1.0, fabs(V[i])));
        excess = MAX(excess, (d - e * 2e-6) / ((fabs(vn) + fabs(V[i] - vn)) * FLT_EPSILON));
    }
    printf("amide def2-SVP potential at %zu points, GPU vs CPU: max |diff| / max(1,|V|) = %.3e, beyond 2e-6 bohr of position %.2f eps of |V_nuc| + |V_el|\n", n, rel, excess);
    EXPECT_LT(rel, 5e-6);
    EXPECT_LT(excess, 2.0);

    // Split into many batches and runs of gaussians: bit identical
    ASSERT_TRUE(int_test_gpu_xyz(gpu_split, stream, &q, p3, n, 0, 20000));
    EXPECT_EQ(0, memcmp(gpu, gpu_split, sizeof(float) * n));

    md_gto_int_charges_free(&q, t.alloc);
    qm_test_free(&t);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
}
#endif

// Every Hermite order the GPU has a kernel for: s, p, d, f and g shells on three centres, so the
// pairs reach g x g (L = 8), with a random symmetric AO matrix.
UTEST(gto_int, all_orders_gpu_match_cpu) {
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        UTEST_SKIP(int_test_no_device_reason());
    }
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));

    // One normalised primitive per shell, l = 0..4 on each of three atoms.
    enum { NA = 3, NL = 5, NS = NA * NL };
    const float atom_xyz[NA][3] = { {0.0f, 0.0f, 0.0f}, {1.6f, 0.4f, -0.3f}, {-0.7f, 1.3f, 0.9f} };
    const float expo[NL] = { 2.1f, 1.3f, 0.9f, 0.7f, 0.55f };
    md_gto_shell_t shells[NS];
    float alpha[NS], coeff[NS];
    for (int a = 0; a < NA; ++a) for (int l = 0; l < NL; ++l) {
        const int s = a * NL + l;
        alpha[s] = expo[l] * (1.0f + 0.1f * a);
        coeff[s] = (float)md_qm_primitive_norm_factor((uint32_t)l, alpha[s]);
        shells[s] = (md_gto_shell_t){ .atom_idx = (uint32_t)a, .primitive_offset = (uint32_t)s, .num_primitives = 1, .l = (uint32_t)l };
    }
    md_gto_basis_t basis = { .num_shells = NS, .num_primitives = NS, .shells = shells, .alpha = alpha, .coeff = coeff };
    const size_t N = md_gto_basis_num_ao(&basis);
    double* D = (double*)md_alloc(alloc, sizeof(double) * N * N);
    uint32_t rng = 12345u;
    for (size_t i = 0; i < N; ++i) for (size_t j = 0; j <= i; ++j) {
        rng = rng * 1664525u + 1013904223u;
        const double v = ((double)(rng >> 8) / (double)(1u << 24) - 0.5) * (i == j ? 1.0 : 0.2);
        D[i * N + j] = D[j * N + i] = v;
    }
    const double Z[NA] = { 3.0, 1.0, 2.0 };

    md_gto_int_charges_t q = {0};
    md_gto_int_charges_desc_t desc = {
        .basis = &basis, .atom_xyz = &atom_xyz[0][0], .density_matrix = D, .density_scale = -1.0,
        .point_xyz = &atom_xyz[0][0], .point_charge = Z, .num_points = NA,
    };
    ASSERT_TRUE(md_gto_int_charges_init(&q, &desc, alloc));
    EXPECT_EQ(8u, q.max_order);
    for (uint32_t L = 0; L <= 8; ++L) EXPECT_GT(q.order_offset[L + 1], q.order_offset[L]);

    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin  = vec3_set(-4.0f, -3.5f, -3.8f),
        .spacing = vec3_set(0.31f, 0.33f, 0.29f),
        .dim     = { 26, 25, 27 },
    };
    const float so[3] = { 0.5f, 0.5f, 0.5f };
    const size_t nv = md_grid_num_points(&grid);
    float* cpu = (float*)md_alloc(alloc, sizeof(float) * nv);
    float* gpu = (float*)md_alloc(alloc, sizeof(float) * nv);
    md_gto_int_potential_grid(cpu, &grid, so, &q);
    ASSERT_TRUE(int_test_gpu_grid(gpu, stream, &q, &grid, so, 0));
    const double rel = int_test_compare(gpu, cpu, nv);
    printf("s..g on three centres, GPU vs CPU: max |diff| / max(1,|V|) = %.3e\n", rel);
    EXPECT_LT(rel, 5e-6);

    md_arena_allocator_destroy(alloc);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
}

// The far field on the device against the exact sum on the CPU, on a cluster of 3^3 waters: a grid
// through and around it, whole; then again into the same texture with md_gto_op ADD for a mask
// (only the bricks holding a masked voxel change); and at points scattered over the same box, the
// worst case for a group's walk.
UTEST(gto_int, far_field_gpu_matches_exact) {
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        UTEST_SKIP(int_test_no_device_reason());
    }
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);

    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(128));
    ASSERT_TRUE(md_molden_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    int_test_qm_t m;
    ASSERT_TRUE(int_test_extract(&m, &t));
    int_test_cluster_t c;
    int_test_cluster(&c, &m, 3, 5.6f, t.alloc);
    md_gto_int_charges_desc_t cd = {
        .basis = &c.basis, .atom_xyz = c.atom_xyz, .density_matrix = c.D, .density_scale = -1.0,
        .point_xyz = c.atom_xyz, .point_charge = c.atom_z, .num_points = c.num_atoms, .threshold = 1e-8,
    };
    md_gto_int_charges_t exact = {0}, far = {0};
    ASSERT_TRUE(md_gto_int_charges_init(&exact, &cd, t.alloc));
    cd.far_field_theta = MD_GTO_INT_FAR_FIELD_THETA;
    ASSERT_TRUE(md_gto_int_charges_init(&far, &cd, t.alloc));

    md_grid_t grid = {
        .orientation = mat3_ident(),
        .origin  = vec3_set(-4.0f, -3.7f, -3.9f),
        .spacing = vec3_set(0.71f, 0.69f, 0.73f),
        .dim     = { 22, 21, 23 },   // not a multiple of the 4x4x4 brick on any axis
    };
    const float so[3] = { 0.5f, 0.5f, 0.5f };
    const size_t nv = md_grid_num_points(&grid);
    float* ref = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    float* gpu = (float*)md_alloc(t.alloc, sizeof(float) * nv);
    md_gto_int_potential_grid(ref, &grid, so, &exact);

    md_gpu_texture_t tex = md_gpu_texture_create(stream, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid.dim[0], .height = (uint32_t)grid.dim[1], .depth_or_layers = (uint32_t)grid.dim[2],
    });
    md_gpu_mem_t rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * nv);
    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, &far);
    ASSERT_TRUE(tex && rb.cpu && gq);

    md_gto_int_gpu_potential_desc_t desc = {
        .charges = gq, .out_tex = tex, .grid = &grid, .sample_offset = { so[0], so[1], so[2] }, .op = MD_GTO_OP_SET,
    };
    ASSERT_TRUE(md_gto_int_gpu_potential_launch(stream, &desc));
    md_gpu_copy_from_texture(stream, rb.gpu, tex, NULL);
    md_gpu_stream_sync(stream);
    MEMCPY(gpu, rb.cpu, sizeof(float) * nv);
    // relative to max(1, |V|), as the other GPU tests: next to a nucleus the potential is tens of au,
    // and float is what limits it there
    const double worst = int_test_compare(gpu, ref, nv);
    printf("27 waters, far field on the GPU vs exact: grid max |dV| / max(1, |V|) %.2e\n", worst);
    EXPECT_LT(worst, 2e-4);

    // A mask of every other brick's one voxel: those bricks gain the potential once more, whole
    uint8_t* mask = (uint8_t*)md_alloc(t.alloc, nv);
    MEMSET(mask, 0, nv);
    for (int z = 0; z < grid.dim[2]; ++z) for (int y = 0; y < grid.dim[1]; ++y) for (int x = 0; x < grid.dim[0]; ++x) {
        if (((x / 4 + y / 4 + z / 4) & 1) == 0 && x % 4 == 1 && y % 4 == 2 && z % 4 == 3) mask[((size_t)z * grid.dim[1] + y) * grid.dim[0] + x] = 1;
    }
    desc.op = MD_GTO_OP_ADD;
    desc.voxel_mask = mask;
    ASSERT_TRUE(md_gto_int_gpu_potential_launch(stream, &desc));
    md_gpu_copy_from_texture(stream, rb.gpu, tex, NULL);
    md_gpu_stream_sync(stream);
    double worst_in = 0.0, worst_out = 0.0;
    for (int z = 0; z < grid.dim[2]; ++z) for (int y = 0; y < grid.dim[1]; ++y) for (int x = 0; x < grid.dim[0]; ++x) {
        const size_t i = ((size_t)z * grid.dim[1] + y) * grid.dim[0] + x;
        // a brick has its masked voxel unless the grid ends before it
        const int bx = x / 4, by = y / 4, bz = z / 4;
        const bool brick_masked = ((bx + by + bz) & 1) == 0 && bx * 4 + 1 < grid.dim[0] && by * 4 + 2 < grid.dim[1] && bz * 4 + 3 < grid.dim[2];
        const float v = ((const float*)rb.cpu)[i];
        if (brick_masked) worst_in = MAX(worst_in, fabs((double)v - 2.0 * (double)gpu[i]));
        else              worst_out = MAX(worst_out, fabs((double)v - (double)gpu[i]));
    }
    EXPECT_LT(worst_in, 1e-5);
    EXPECT_EQ(0.0, worst_out);

    // Scattered points
    enum { NP = 300 };
    float pts[NP][3];
    uint32_t seed = 4242;
    for (int i = 0; i < NP; ++i) for (int k = 0; k < 3; ++k) {
        seed = seed * 1664525u + 1013904223u;
        pts[i][k] = -4.0f + 19.0f * (float)(seed >> 8) / (float)(1u << 24);
    }
    double Ve[NP];
    float  Vg[NP];
    md_gto_int_potential_xyz(Ve, NULL, &pts[0][0], NP, 0, &exact);
    ASSERT_TRUE(int_test_gpu_xyz(Vg, stream, &far, &pts[0][0], NP, 0, 0));
    double worst_p = 0.0;
    for (int i = 0; i < NP; ++i) worst_p = MAX(worst_p, fabs((double)Vg[i] - Ve[i]) / MAX(1.0, fabs(Ve[i])));
    printf("27 waters, far field on the GPU vs exact: points max |dV| / max(1, |V|) %.2e\n", worst_p);
    EXPECT_LT(worst_p, 2e-4);

    // The same points shuffled, in a float4 layout: the device sorts them, so each gets the same value
    // to the bit, whatever the order they came in
    uint32_t perm[NP];
    for (int i = 0; i < NP; ++i) perm[i] = (uint32_t)i;
    for (int i = NP - 1; i > 0; --i) {
        seed = seed * 1664525u + 1013904223u;
        const int j = (int)((seed >> 8) % (uint32_t)(i + 1));
        const uint32_t tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    float shuffled[NP][4];
    for (int i = 0; i < NP; ++i) {
        for (int k = 0; k < 3; ++k) shuffled[i][k] = pts[perm[i]][k];
        shuffled[i][3] = 1.0f;
    }
    float Vs[NP];
    ASSERT_TRUE(int_test_gpu_xyz(Vs, stream, &far, &shuffled[0][0], NP, 16, 0));
    int differ = 0;
    for (int i = 0; i < NP; ++i) differ += Vs[i] != Vg[perm[i]];
    EXPECT_EQ(0, differ);

    // As jobs, in small pieces polled without waiting: the same, to the bit
    {
        md_gto_int_gpu_potential_desc_t jd = {
            .charges = gq, .out_tex = tex, .grid = &grid, .sample_offset = { so[0], so[1], so[2] }, .op = MD_GTO_OP_SET,
            .work_per_dispatch = 1u << 20,
        };
        md_gto_int_gpu_job_t* job = md_gto_int_gpu_potential_job_begin(stream, &jd);
        ASSERT_TRUE(job != NULL);
        int polls = 0;
        float last = 0.0f;
        bool monotone = true;
        md_gto_int_gpu_job_status_t st;
        while ((st = md_gto_int_gpu_job_poll(job, 0.0)) == MD_GTO_INT_GPU_JOB_RUNNING) {
            const float p = md_gto_int_gpu_job_progress(job);
            monotone = monotone && p >= last && p <= 1.0f;
            last = p;
            polls++;
        }
        EXPECT_EQ(MD_GTO_INT_GPU_JOB_DONE, st);
        EXPECT_EQ(1.0f, md_gto_int_gpu_job_progress(job));
        EXPECT_TRUE(md_gto_int_gpu_job_end(job));
        EXPECT_TRUE(monotone);
        md_gpu_copy_from_texture(stream, rb.gpu, tex, NULL);
        md_gpu_stream_sync(stream);
        int diff = 0;
        for (size_t i = 0; i < nv; ++i) diff += ((const float*)rb.cpu)[i] != gpu[i];
        EXPECT_EQ(0, diff);
        printf("27 waters, the grid as a job: %d polls while it ran, %d voxels differ from the launch\n", polls, diff);

        // ended at once: cancelled
        job = md_gto_int_gpu_potential_job_begin(stream, &jd);
        ASSERT_TRUE(job != NULL);
        EXPECT_FALSE(md_gto_int_gpu_job_end(job));
        md_gpu_stream_sync(stream);

        // the points
        md_gpu_mem_t pbuf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(pts));
        md_gpu_mem_t vbuf = md_gpu_malloc(stream, MD_GPU_MEM_DEVICE, sizeof(float) * NP);
        md_gpu_mem_t vrb  = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * NP);
        ASSERT_TRUE(pbuf.gpu && vbuf.gpu && vrb.cpu);
        MEMCPY(md_gpu_upload_begin(stream, pbuf.gpu, sizeof(pts)), pts, sizeof(pts));
        md_gpu_upload_end(stream);
        md_gto_int_gpu_potential_xyz_desc_t xd = {
            .charges = gq, .xyz = pbuf.gpu, .num_xyz = NP, .out_potential = vbuf.gpu, .work_per_dispatch = 1u << 20,
        };
        job = md_gto_int_gpu_potential_xyz_job_begin(stream, &xd);
        ASSERT_TRUE(job != NULL);
        EXPECT_EQ(MD_GTO_INT_GPU_JOB_DONE, md_gto_int_gpu_job_poll(job, -1.0));
        EXPECT_TRUE(md_gto_int_gpu_job_end(job));
        md_gpu_copy(stream, vrb.gpu, vbuf.gpu, sizeof(float) * NP);
        md_gpu_stream_sync(stream);
        diff = 0;
        for (int i = 0; i < NP; ++i) diff += ((const float*)vrb.cpu)[i] != Vg[i];
        EXPECT_EQ(0, diff);
        md_gpu_free(stream, pbuf.gpu);
        md_gpu_free(stream, vbuf.gpu);
        md_gpu_free(stream, vrb.gpu);
    }

    md_gto_int_gpu_charges_destroy(stream, gq);
    md_gpu_free(stream, rb.gpu);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_sync(stream);
    md_gto_int_charges_free(&exact, t.alloc);
    md_gto_int_charges_free(&far, t.alloc);
    qm_test_free(&t);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
}

#ifdef MD_HDF5
// Whether the device has not been lost (a lost device makes every later call return at once, and a
// timing of nothing)
static bool int_test_gpu_alive(void) {
    const char* e = md_gpu_last_error();
    if (e && strstr(e, "DEVICE_LOST")) {
        printf("  the device was lost: %s\n", e);
        return false;
    }
    return true;
}

// The electrostatic potential of a large calculation as viamd evaluates it, timed: the distribution
// with a far field, the band of a surface on the GPU (voxels 2.5 - 4.5 bohr from the nearest atom,
// in 0.4 bohr voxels, as the bricks holding them), the same voxels as points, and both against the
// exact sum at some of them. GTO_INT_FILE=<file in test_data/vlx>, e.g. protein.h5.
UTEST(gto_int, benchmark) {
    const char* file = getenv("GTO_INT_FILE");
    if (!file || !file[0]) UTEST_SKIP("set GTO_INT_FILE=<file in test_data/vlx>");
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) UTEST_SKIP(int_test_no_device_reason());
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);
    md_gpu_device_info_t dinfo;
    md_gpu_device_info(device, &dinfo);

    char path[512];
    snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/vlx/%s", file);
    vlx_test_t t = {0};
    md_tick_t t0 = md_tick_now();
    ASSERT_TRUE(vlx_test_load(&t, str_from_cstr(path), MEGABYTES(1024)));
    md_gto_basis_t basis = {0};
    ASSERT_TRUE(qm_test_basis(&basis, &t));
    const size_t na = t.sys.atom.count;
    float* xyz = (float*)malloc(sizeof(float) * 3 * na);
    double* Z = (double*)malloc(sizeof(double) * na);
    ASSERT_EQ(na, vlx_test_atom_xyz_bohr(xyz, 3 * na, &t));
    ASSERT_EQ(na, qm_test_series(Z, na, &t, STR_LIT("qm/atom/nuclear_charge")));
    double* rows = NULL; double* w = NULL; size_t num_ao = 0;
    const size_t num_mo = md_qm_density_orbitals(&rows, &w, &num_ao, &t.sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
    ASSERT_TRUE(num_mo > 0 && num_ao == md_gto_basis_num_ao(&basis));
    double* D = (double*)malloc(sizeof(double) * num_ao * num_ao);
    ASSERT_TRUE(md_gto_density_from_orbitals(D, rows, num_ao, 1, w, num_mo, num_ao));
    printf("\n  md_gto_int far field %s\n  device: %s | %s: %zu atoms, %zu AOs (loaded, and D, in %.1f s)\n",
           MD_GTO_INT_FAR_FIELD_REVISION, dinfo.name, file, na, num_ao, md_tick_to_seconds(md_tick_now() - t0));

    md_allocator_i* heap = md_get_heap_allocator();
    md_gto_int_charges_desc_t cd = {
        .basis = &basis, .atom_xyz = xyz, .density_matrix = D, .density_scale = -1.0,
        .point_xyz = xyz, .point_charge = Z, .num_points = na, .threshold = 1e-6,
    };
    md_gto_int_charges_t exact = {0}, far = {0};
    t0 = md_tick_now();
    ASSERT_TRUE(md_gto_int_charges_init(&exact, &cd, heap));
    const double t_exact = md_tick_to_seconds(md_tick_now() - t0);
    cd.far_field_theta = MD_GTO_INT_FAR_FIELD_THETA;
    t0 = md_tick_now();
    ASSERT_TRUE(md_gto_int_charges_init(&far, &cd, heap));
    const double t_far = md_tick_to_seconds(md_tick_now() - t0);
    printf("  distribution: %u gaussians, built in %.2f s; with the far field (%u nodes) in %.2f s; work per point %llu, exactly %llu\n",
           far.num_gaussians, t_exact, far.num_nodes, t_far,
           (unsigned long long)md_gto_int_charges_work_per_point(&far), (unsigned long long)md_gto_int_charges_work_per_point(&exact));

    // the grid: the atoms' box and 6 bohr, 0.4 bohr voxels; the band: 2.5 - 4.5 bohr from the nearest atom
    float lo[3] = { 1e30f, 1e30f, 1e30f }, hi[3] = { -1e30f, -1e30f, -1e30f };
    for (size_t i = 0; i < na; ++i) for (int k = 0; k < 3; ++k) { lo[k] = MIN(lo[k], xyz[3 * i + k]); hi[k] = MAX(hi[k], xyz[3 * i + k]); }
    const float h = 0.4f;
    md_grid_t grid = { .orientation = mat3_ident(), .origin = vec3_set(lo[0] - 6.0f, lo[1] - 6.0f, lo[2] - 6.0f), .spacing = vec3_set(h, h, h) };
    for (int k = 0; k < 3; ++k) grid.dim[k] = (int)ceilf((hi[k] - lo[k] + 12.0f) / h);
    const size_t nv = md_grid_num_points(&grid);
    float* dmin = (float*)malloc(sizeof(float) * nv);
    uint8_t* mask = (uint8_t*)malloc(nv);
    for (size_t i = 0; i < nv; ++i) dmin[i] = 1e30f;
    const int rr = (int)ceilf(4.5f / h) + 1;
    for (size_t a = 0; a < na; ++a) {
        int c[3];
        for (int k = 0; k < 3; ++k) c[k] = (int)((xyz[3 * a + k] - grid.origin.elem[k]) / h - 0.5f);
        for (int z = MAX(0, c[2] - rr); z <= MIN(grid.dim[2] - 1, c[2] + rr); ++z)
        for (int y = MAX(0, c[1] - rr); y <= MIN(grid.dim[1] - 1, c[1] + rr); ++y)
        for (int x = MAX(0, c[0] - rr); x <= MIN(grid.dim[0] - 1, c[0] + rr); ++x) {
            const float p[3] = { grid.origin.x + (x + 0.5f) * h, grid.origin.y + (y + 0.5f) * h, grid.origin.z + (z + 0.5f) * h };
            const float d2 = (p[0] - xyz[3 * a]) * (p[0] - xyz[3 * a]) + (p[1] - xyz[3 * a + 1]) * (p[1] - xyz[3 * a + 1]) + (p[2] - xyz[3 * a + 2]) * (p[2] - xyz[3 * a + 2]);
            float* dm = dmin + ((size_t)z * grid.dim[1] + y) * grid.dim[0] + x;
            *dm = MIN(*dm, d2);
        }
    }
    size_t num_band = 0;
    for (size_t i = 0; i < nv; ++i) { mask[i] = (dmin[i] >= 2.5f * 2.5f && dmin[i] <= 4.5f * 4.5f) ? 1 : 0; num_band += mask[i]; }
    printf("  grid %d x %d x %d (%zu voxels), band %zu voxels\n", grid.dim[0], grid.dim[1], grid.dim[2], nv, num_band);

    md_gpu_texture_t tex = md_gpu_texture_create(stream, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)grid.dim[0], .height = (uint32_t)grid.dim[1], .depth_or_layers = (uint32_t)grid.dim[2],
    });
    md_gpu_mem_t rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, sizeof(float) * nv);
    ASSERT_TRUE(tex && rb.cpu);
    t0 = md_tick_now();
    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, &far);
    md_gpu_stream_sync(stream);
    ASSERT_TRUE(gq != NULL);
    printf("  on the device in %.2f s\n", md_tick_to_seconds(md_tick_now() - t0));

    // First a sixteenth of the band (every 16th layer of bricks): it makes the kernels, and gives a
    // rate even when the whole band does not go through
    uint8_t* part = (uint8_t*)malloc(nv);
    size_t num_part = 0;
    for (size_t i = 0; i < nv; ++i) {
        const size_t z = i / ((size_t)grid.dim[0] * grid.dim[1]);
        part[i] = mask[i] && (z / 4) % 16 == 0;
        num_part += part[i];
    }
    md_gto_int_gpu_potential_desc_t pd = {
        .charges = gq, .out_tex = tex, .grid = &grid, .sample_offset = { 0.5f, 0.5f, 0.5f }, .op = MD_GTO_OP_SET, .voxel_mask = part,
    };
    for (int r = 0; r < 2; ++r) {
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_int_gpu_potential_launch(stream, &pd));
        md_gpu_stream_sync(stream);
        const double ts = md_tick_to_seconds(md_tick_now() - t0);
        ASSERT_TRUE(int_test_gpu_alive());
        printf("  a sixteenth of the band on the GPU (%zu voxels): %.2f s (%.1f us per voxel)%s\n", num_part, ts, ts * 1e6 / (double)MAX(num_part, 1),
               r == 0 ? ", with the kernels made" : "");
        if (ts > 10.0) break;
    }
    free(part);

    pd.voxel_mask = mask;
    double best = 1e30;
    for (int r = 0; r < 2; ++r) {   // a second only where it is quick
        t0 = md_tick_now();
        ASSERT_TRUE(md_gto_int_gpu_potential_launch(stream, &pd));
        md_gpu_stream_sync(stream);
        best = MIN(best, md_tick_to_seconds(md_tick_now() - t0));
        ASSERT_TRUE(int_test_gpu_alive());
        if (best > 10.0) break;
    }
    md_gpu_copy_from_texture(stream, rb.gpu, tex, NULL);
    md_gpu_stream_sync(stream);
    printf("  band on the GPU: %.2f s (%.1f us per band voxel)\n", best, best * 1e6 / (double)MAX(num_band, 1));

    // The band as a job, as viamd runs it: polled for 2 ms of every 16 ms frame
    if (best < 10.0) {
        t0 = md_tick_now();
        md_gto_int_gpu_job_t* job = md_gto_int_gpu_potential_job_begin(stream, &pd);
        ASSERT_TRUE(job != NULL);
        int frames = 0;
        double longest_ms = 0.0;
        md_gto_int_gpu_job_status_t st;
        for (;;) {
            const md_tick_t f0 = md_tick_now();
            st = md_gto_int_gpu_job_poll(job, 2.0);
            const double ms = md_tick_to_milliseconds(md_tick_now() - f0);
            longest_ms = MAX(longest_ms, ms);
            frames++;
            if (st != MD_GTO_INT_GPU_JOB_RUNNING) break;
            if (ms < 16.0) md_thread_sleep((size_t)(16.0 - ms));
        }
        EXPECT_EQ(MD_GTO_INT_GPU_JOB_DONE, st);
        md_gto_int_gpu_job_end(job);
        ASSERT_TRUE(int_test_gpu_alive());
        printf("  band as a job, polled 2 ms a 16 ms frame: %.2f s over %d frames, the longest poll %.1f ms\n",
               md_tick_to_seconds(md_tick_now() - t0), frames, longest_ms);
    }

    // the band's voxels as points, shuffled: the device sorts them, so this is the time of any order
    size_t* idx = (size_t*)malloc(sizeof(size_t) * num_band);
    size_t np = 0;
    for (size_t i = 0; i < nv; ++i) if (mask[i]) idx[np++] = i;
    uint64_t seed = 0x9E3779B97F4A7C15ull;
    for (size_t i = np; i > 1; --i) {
        seed = seed * 6364136223846793005ull + 1442695040888963407ull;
        const size_t j = (size_t)((seed >> 33) % i);
        const size_t tmp = idx[i - 1]; idx[i - 1] = idx[j]; idx[j] = tmp;
    }
    float* pts = (float*)malloc(sizeof(float) * 3 * np);
    for (size_t p = 0; p < np; ++p) {
        const size_t i = idx[p];
        const size_t x = i % grid.dim[0], y = (i / grid.dim[0]) % grid.dim[1], z = i / ((size_t)grid.dim[0] * grid.dim[1]);
        pts[3 * p + 0] = grid.origin.x + (x + 0.5f) * h;
        pts[3 * p + 1] = grid.origin.y + (y + 0.5f) * h;
        pts[3 * p + 2] = grid.origin.z + (z + 0.5f) * h;
    }
    float* Vp = (float*)malloc(sizeof(float) * np);
    // a sixteenth of them first, which makes the kernels
    t0 = md_tick_now();
    ASSERT_TRUE(int_test_gpu_xyz(Vp, stream, &far, pts, np / 16, 0, 0));
    printf("  a sixteenth of them as points, with the kernels made: %.2f s\n", md_tick_to_seconds(md_tick_now() - t0));
    t0 = md_tick_now();
    ASSERT_TRUE(int_test_gpu_xyz(Vp, stream, &far, pts, np, 0, 0));
    ASSERT_TRUE(int_test_gpu_alive());
    printf("  the same voxels as points, shuffled, on the GPU: %.2f s (with the upload and the distribution again)\n", md_tick_to_seconds(md_tick_now() - t0));

    // against the exact sum at some of them
    enum { NS = 200 };
    float  sp[NS][3];
    double Ve[NS], Vc[NS];
    size_t pick[NS];
    for (int i = 0; i < NS; ++i) {
        pick[i] = (size_t)((uint64_t)i * np / NS);
        for (int k = 0; k < 3; ++k) sp[i][k] = pts[3 * pick[i] + k];
    }
    t0 = md_tick_now();
    md_gto_int_potential_xyz(Ve, NULL, &sp[0][0], NS, 0, &exact);
    const double t_ce = md_tick_to_seconds(md_tick_now() - t0) / NS;
    t0 = md_tick_now();
    md_gto_int_potential_xyz(Vc, NULL, &sp[0][0], NS, 0, &far);
    const double t_cf = md_tick_to_seconds(md_tick_now() - t0) / NS;
    double eg = 0.0, ep = 0.0, ec = 0.0, span_lo = 1e30, span_hi = -1e30;
    for (int i = 0; i < NS; ++i) {
        eg = MAX(eg, fabs((double)((const float*)rb.cpu)[idx[pick[i]]] - Ve[i]));
        ep = MAX(ep, fabs((double)Vp[pick[i]] - Ve[i]));
        ec = MAX(ec, fabs(Vc[i] - Ve[i]));
        span_lo = MIN(span_lo, Ve[i]); span_hi = MAX(span_hi, Ve[i]);
    }
    printf("  against the exact sum at %d band voxels (V from %.3f to %.3f au): grid %.2e, points %.2e, CPU far field %.2e au\n",
           NS, span_lo, span_hi, eg, ep, ec);
    printf("  CPU per point: exact %.2f ms, far field %.2f ms\n", t_ce * 1e3, t_cf * 1e3);
    EXPECT_LT(eg, 5e-4);
    EXPECT_LT(ep, 5e-4);

    free(Vp); free(idx); free(pts); free(mask); free(dmin);
    md_gto_int_gpu_charges_destroy(stream, gq);
    md_gpu_free(stream, rb.gpu);
    md_gpu_texture_destroy(tex);
    md_gpu_stream_sync(stream);
    md_gto_int_charges_free(&exact, heap);
    md_gto_int_charges_free(&far, heap);
    md_free(heap, rows, sizeof(double) * num_mo * num_ao);
    md_free(heap, w, sizeof(double) * num_mo);
    free(D); free(Z); free(xyz);
    qm_test_free(&t);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
}
#endif

#endif

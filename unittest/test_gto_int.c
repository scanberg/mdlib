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

    out->atom_z = (double*)md_alloc(t->alloc, sizeof(double) * out->num_atoms);
    if (qm_test_series(out->atom_z, out->num_atoms, t, STR_LIT("qm/atom/atomic_number")) != out->num_atoms) return false;

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
    md_gto_int_gpu_potential_launch(stream, &desc);
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

#endif

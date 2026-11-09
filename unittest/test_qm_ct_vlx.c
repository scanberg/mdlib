#include "utest.h"

#include "vlx_test_util.h"

#include <md_qm_ct.h>
#include <core/md_allocator.h>

#include <math.h>

// The charge transfer analysis on real VeloxChem output, read the way a consumer reads it: the
// attribute table's orbitals, overlap and NTOs, and the AO to atom map of its basis. The reference
// numbers come from an independent evaluation of Omega_AB = sum T~^2 from the response vector in
// the file's own spherical basis (numpy). The table publishes the Cartesian embedding of that basis,
// whose Lowdin populations differ at the level of a few 1e-3, hence the tolerance.

typedef struct ct_vlx_t {
    size_t num_ao, num_mo, num_atoms, K;
    double* omega;      // {N,N} atoms, normalized to sum 1
    double  lambda_sum;
    double  omega_sum;
} ct_vlx_t;

static bool ct_vlx_compute(ct_vlx_t* out, vlx_test_t* t, uint32_t state) {
    md_allocator_i* alloc = t->alloc;

    size_t num_ao = 0;
    const double* S = qm_test_matrix(t, MD_STR_LIT("basis/overlap"), &num_ao);
    const md_attribute_t* coeff  = qm_test_attr(t, MD_STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* lam_a  = qm_test_attr(t, MD_STR_LIT("vlx/rsp/nto/lambda"));
    const md_attribute_t* hole_a = qm_test_attr(t, MD_STR_LIT("vlx/rsp/nto/hole/coefficient"));
    const md_attribute_t* part_a = qm_test_attr(t, MD_STR_LIT("vlx/rsp/nto/particle/coefficient"));
    if (!S || !coeff || !lam_a || !hole_a || !part_a) return false;

    const size_t num_mo = coeff->format.shape[0];
    if (coeff->format.shape[1] != num_ao) return false;
    double* C = md_alloc(alloc, sizeof(double) * num_mo * num_ao);
    if (md_attribute_extract_f64(C, num_mo * num_ao, coeff, md_attribute_slice_all(), md_unit_none()) != num_mo * num_ao) return false;

    const md_attribute_slice_t slice = md_attribute_slice_1(state);
    const size_t K = lam_a->format.shape[1];
    double* lambda = md_alloc(alloc, sizeof(double) * K);
    double* hole   = md_alloc(alloc, sizeof(double) * K * num_ao);
    double* part   = md_alloc(alloc, sizeof(double) * K * num_ao);
    if (md_attribute_extract_f64(lambda, K, lam_a, slice, md_unit_none()) != K) return false;
    if (md_attribute_extract_f64(hole, K * num_ao, hole_a, slice, md_unit_none()) != K * num_ao) return false;
    if (md_attribute_extract_f64(part, K * num_ao, part_a, slice, md_unit_none()) != K * num_ao) return false;

    md_gto_basis_t basis = {0};
    if (!qm_test_basis(&basis, t) || md_gto_basis_num_ao(&basis) != num_ao) return false;
    uint32_t* ao_atom = md_alloc(alloc, sizeof(uint32_t) * num_ao);
    if (md_gto_basis_ao_to_atom(ao_atom, &basis) != num_ao) return false;
    size_t num_atoms = 0;
    for (size_t i = 0; i < num_ao; ++i) num_atoms = ao_atom[i] + 1 > num_atoms ? ao_atom[i] + 1 : num_atoms;

    double* Ct = md_alloc(alloc, sizeof(double) * num_mo * num_ao);
    if (!md_qm_lowdin_orbitals(Ct, C, num_mo, num_ao, NULL)) return false;

    double* ht = md_alloc(alloc, sizeof(double) * K * num_ao);
    double* pt = md_alloc(alloc, sizeof(double) * K * num_ao);
    md_qm_lowdin_transform(ht, hole, K, S, C, Ct, num_mo, num_ao);
    md_qm_lowdin_transform(pt, part, K, S, C, Ct, num_mo, num_ao);

    double* Hb = md_alloc(alloc, sizeof(double) * num_atoms * K * K);
    double* Pb = md_alloc(alloc, sizeof(double) * num_atoms * K * K);
    md_qm_ct_group_blocks(Hb, ht, K, num_ao, ao_atom, num_atoms);
    md_qm_ct_group_blocks(Pb, pt, K, num_ao, ao_atom, num_atoms);

    double* omega = md_alloc(alloc, sizeof(double) * num_atoms * num_atoms);
    md_qm_ct_omega(omega, lambda, K, Hb, num_atoms, Pb, num_atoms);

    double lsum = 0.0, osum = 0.0;
    for (size_t k = 0; k < K; ++k) lsum += lambda[k];
    for (size_t i = 0; i < num_atoms * num_atoms; ++i) osum += omega[i];
    for (size_t i = 0; i < num_atoms * num_atoms; ++i) omega[i] /= osum;

    *out = (ct_vlx_t){ num_ao, num_mo, num_atoms, K, omega, lsum, osum };
    return true;
}

static double ct_group_sum(const ct_vlx_t* ct, const int* from, size_t num_from, const int* to, size_t num_to) {
    double s = 0.0;
    for (size_t i = 0; i < num_from; ++i) for (size_t j = 0; j < num_to; ++j) s += ct->omega[from[i] * ct->num_atoms + to[j]];
    return s;
}

// N,N-dimethylformamide, S1: an n -> pi* excitation that stays on the formamide group.
UTEST(qm_ct_vlx, amide_s1_is_local) {
    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/amide.h5"), MD_MEGABYTES(64)));
    ct_vlx_t ct = {0};
    ASSERT_TRUE(ct_vlx_compute(&ct, &t, 0));
    ASSERT_EQ(ct.num_atoms, 12u);

    // every pair kept: Omega carries exactly the norm of the excitation
    EXPECT_NEAR(ct.omega_sum, ct.lambda_sum, 1.0e-9 * ct.lambda_sum);
    for (size_t i = 0; i < ct.num_atoms * ct.num_atoms; ++i) EXPECT_GE(ct.omega[i], 0.0);

    const int formamide[] = { 0, 1, 2, 5 };
    const int methyl_a[]  = { 3, 6, 7, 8 };
    const int methyl_b[]  = { 4, 9, 10, 11 };
    EXPECT_NEAR(ct_group_sum(&ct, formamide, 4, formamide, 4), 0.901, 0.01);
    EXPECT_NEAR(ct_group_sum(&ct, formamide, 4, methyl_a, 4),  0.037, 0.01);
    EXPECT_NEAR(ct_group_sum(&ct, formamide, 4, methyl_b, 4),  0.030, 0.01);

    // between ATOMS most of it is off the diagonal - which is why atoms are the wrong level for
    // charge transfer, and fragments the right one
    double diag = 0.0;
    for (size_t a = 0; a < ct.num_atoms; ++a) diag += ct.omega[a * ct.num_atoms + a];
    EXPECT_NEAR(1.0 - diag, 0.741, 0.02);

    qm_test_free(&t);
}

// A secondary amine on a phenylpropyl chain, S1: the amine lone pair to the ring.
UTEST(qm_ct_vlx, mol_s1_amine_to_phenyl) {
    vlx_test_t t = {0};
    ASSERT_TRUE(vlx_test_load(&t, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.h5"), MD_MEGABYTES(128)));
    ct_vlx_t ct = {0};
    ASSERT_TRUE(ct_vlx_compute(&ct, &t, 0));
    ASSERT_EQ(ct.num_atoms, 26u);
    EXPECT_NEAR(ct.omega_sum, ct.lambda_sum, 1.0e-9 * ct.lambda_sum);

    const int amine[]  = { 0, 7, 17, 20, 21, 22 };
    const int linker[] = { 1, 2, 4, 11, 12, 13, 14, 15, 16 };
    const int phenyl[] = { 3, 5, 6, 8, 9, 10, 18, 19, 23, 24, 25 };
    EXPECT_NEAR(ct_group_sum(&ct, amine,  6, phenyl, 11), 0.787, 0.01);
    EXPECT_NEAR(ct_group_sum(&ct, linker, 9, phenyl, 11), 0.137, 0.01);
    EXPECT_NEAR(ct_group_sum(&ct, phenyl, 11, phenyl, 11), 0.048, 0.01);

    qm_test_free(&t);
}

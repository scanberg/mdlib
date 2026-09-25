#include "utest.h"

#include "qm_test_util.h"

#include <md_trexio.h>
#include <md_molden.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>

// Everything a TREXIO file carries reaches a consumer as an md_system_t and its ATTRIBUTE TABLE:
// md_trexio.h has two entry points and no reader object to ask questions of. So these tests read
// what a real consumer reads, by the same paths.
//
// See test_data/trexio/README.md for where the file came from.

static bool trexio_load(qm_test_t* t, const char* path, size_t arena_bytes) {
    qm_test_init(t, arena_bytes);
    return md_trexio_system_init_from_file(&t->sys, &t->state, str_from_cstr(path));
}

UTEST(trexio, parse) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MEGABYTES(16)));

    EXPECT_EQ(3u, t.sys.atom.count);
    EXPECT_EQ(3u, t.state.num_atoms);

    // TREXIO stores coordinates in bohr; a system's are in Angstrom, and so is what is published.
    const md_attribute_t* coord = qm_test_attr(&t, STR_LIT("qm/atom/coordinate"));
    ASSERT_TRUE(coord != NULL);
    ASSERT_EQ(md_attribute_components(&coord->format), 3u);
    ASSERT_EQ(md_attribute_value_count(&coord->format), 3u);

    double xyz[9] = {0};
    ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_unit_none()), ARRAY_SIZE(xyz));
    EXPECT_NEAR(4.999999640091311, xyz[0], 1.0e-9);
    EXPECT_NEAR(7.147076485540978, xyz[1], 1.0e-9);
    EXPECT_NEAR(7.650970449269813, xyz[2], 1.0e-9);

    // ...and they are the same three numbers the Molden file beside it was written from.
    EXPECT_NEAR(4.068065707173540, xyz[3], 1.0e-9);

    EXPECT_NEAR(xyz[0], (double)t.state.xyz[0].x, 1.0e-4);

    // nucleus_label is authoritative for the element; nucleus_charge is what the electrons see and
    // is not the atomic number for an effective core potential.
    double z[3] = {0};
    const md_attribute_t* za = qm_test_attr(&t, STR_LIT("qm/atom/atomic_number"));
    ASSERT_TRUE(za != NULL);
    ASSERT_EQ(md_attribute_extract_f64(z, 3, za, md_unit_none()), 3u);
    EXPECT_EQ(8.0, z[0]);
    EXPECT_EQ(1.0, z[1]);
    EXPECT_EQ(1.0, z[2]);

    EXPECT_EQ(10.0, qm_test_scalar(&t, STR_LIT("trexio/electron_count/total"), -1.0));
    EXPECT_EQ( 5.0, qm_test_scalar(&t, STR_LIT("trexio/electron_count/up"),    -1.0));
    EXPECT_EQ( 5.0, qm_test_scalar(&t, STR_LIT("trexio/electron_count/down"),  -1.0));

    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("trexio/basis_type")),      STR_LIT("Gaussian")));
    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("trexio/mo_type")),         STR_LIT("MO")));
    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("trexio/ao_convention")),   STR_LIT("spherical")));
    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("trexio/metadata/package_version")), STR_LIT("2.6.1")));

    // This file carries no mo_energy and no nucleus_repulsion. An absent path is how a consumer
    // learns a block is missing; a column of zeros would read as a legitimate set of energies.
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("orbital/alpha/energy")));
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("trexio/nuclear_repulsion_energy")));

    md_gto_basis_t basis = {0};
    ASSERT_TRUE(qm_test_basis(&basis, &t));
    EXPECT_EQ(12u, basis.num_shells);
    EXPECT_EQ(22u, basis.num_primitives);
    EXPECT_EQ(25u, (uint32_t)md_gto_basis_num_ao(&basis));      // Cartesian, whatever the file stored
    EXPECT_EQ(24u, (uint32_t)md_gto_basis_num_sph_ao(&basis));  // which is what the file did store

    qm_test_free(&t);
}

// TREXIO states one occupation per orbital over BOTH spins, and says whether that is what it means
// through mo_spin. What the orbital/ tree wants is per channel, so a restricted file is split - and
// the file's own numbers stay reachable under trexio/mo/, which is the part that would otherwise be
// lost.
UTEST(trexio, restricted_occupations_are_split_per_spin) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MEGABYTES(16)));

    double file_occ[24] = {0};
    ASSERT_EQ(qm_test_series(file_occ, ARRAY_SIZE(file_occ), &t, STR_LIT("trexio/mo/occupation")), ARRAY_SIZE(file_occ));
    EXPECT_NEAR(2.0, file_occ[0], 1.0e-12);
    EXPECT_NEAR(2.0, file_occ[4], 1.0e-12);
    EXPECT_NEAR(0.0, file_occ[5], 1.0e-12);

    double occ[24] = {0};
    ASSERT_EQ(qm_test_series(occ, ARRAY_SIZE(occ), &t, STR_LIT("orbital/alpha/occupation")), ARRAY_SIZE(occ));
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR(1.0, occ[i], 1.0e-12);
    }
    EXPECT_NEAR(0.0, occ[5], 1.0e-12);

    // Beta is a second NAME for alpha here, not a second copy.
    const md_attribute_t* a = qm_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* b = qm_test_attr(&t, STR_LIT("orbital/beta/coefficient"));
    ASSERT_TRUE(a != NULL);
    ASSERT_TRUE(b != NULL);
    EXPECT_TRUE(md_attribute_same_data(a, b));

    qm_test_free(&t);
}

// The one test that pins every convention at once - see the note on the Molden equivalent. TREXIO
// orders its spherical functions m = 0, +1, -1, +2, -2 and md_gto orders them ascending in m, so an
// unpermuted read of this file gives orbitals that are visibly not orthonormal.
UTEST(trexio, orbitals_are_orthonormal) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MEGABYTES(32)));
    EXPECT_LT(qm_test_orthonormality(&t, STR_LIT("orbital/alpha/coefficient")), 1.0e-5);
    qm_test_free(&t);
}

UTEST(trexio, total_density_holds_ten_electrons) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MEGABYTES(32)));
    EXPECT_NEAR(10.0, qm_test_electron_count(&t, STR_LIT("orbital/total/density")), 1.0e-5);
    EXPECT_NEAR( 5.0, qm_test_electron_count(&t, STR_LIT("orbital/alpha/density")), 1.0e-5);
    EXPECT_NEAR( 0.0, qm_test_electron_count(&t, STR_LIT("orbital/difference/density")), 1.0e-9);
    qm_test_free(&t);
}

// THE cross format test, and the reason the Molden file in test_data was converted from this one
// rather than found somewhere: the two readers have nothing in common but the contract, they take
// completely different routes to it - one parses text and the other reads HDF5, and they disagree
// about primitive normalisation, about atomic orbital order and about where the occupations live -
// and the system that comes out the far end must be the same system.
//
// If either reader drifts, this separates immediately, and it says which one by which side moved.
UTEST(trexio, agrees_with_the_same_calculation_in_molden) {
    qm_test_t tx = {0};
    qm_test_t md = {0};
    ASSERT_TRUE(trexio_load(&tx, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MEGABYTES(32)));
    qm_test_init(&md, MEGABYTES(32));
    ASSERT_TRUE(md_molden_system_init_from_file(&md.sys, &md.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));

    ASSERT_EQ(tx.sys.atom.count, md.sys.atom.count);

    md_gto_basis_t a = {0};
    md_gto_basis_t b = {0};
    ASSERT_TRUE(qm_test_basis(&a, &tx));
    ASSERT_TRUE(qm_test_basis(&b, &md));
    ASSERT_EQ(a.num_shells,     b.num_shells);
    ASSERT_EQ(a.num_primitives, b.num_primitives);

    // Both files state a spherical basis, so both contractions are normalised against the spherical
    // atomic orbital and the published exponents and coefficients must agree outright.
    for (uint32_t s = 0; s < a.num_shells; ++s) {
        EXPECT_EQ(a.shells[s].l,              b.shells[s].l);
        EXPECT_EQ(a.shells[s].atom_idx,       b.shells[s].atom_idx);
        EXPECT_EQ(a.shells[s].num_primitives, b.shells[s].num_primitives);
    }
    for (uint32_t p = 0; p < a.num_primitives; ++p) {
        EXPECT_NEAR(a.alpha[p], b.alpha[p], 1.0e-4f * fabsf(a.alpha[p]) + 1.0e-6f);
        EXPECT_NEAR(a.coeff[p], b.coeff[p], 1.0e-4f * fabsf(a.coeff[p]) + 1.0e-6f);
    }

    // ...and so must the molecular orbitals, in the Cartesian order both publish them in.
    size_t na = 0, nb = 0;
    double* ca = qm_test_matrix(&tx, STR_LIT("orbital/alpha/coefficient"), &na);
    double* cb = qm_test_matrix(&md, STR_LIT("orbital/alpha/coefficient"), &nb);
    ASSERT_TRUE(ca != NULL);
    ASSERT_TRUE(cb != NULL);
    ASSERT_EQ(na, nb);
    ASSERT_EQ(25u, (uint32_t)na);

    double worst = 0.0;
    for (size_t i = 0; i < 24 * na; ++i) {
        const double e = fabs(ca[i] - cb[i]);
        worst = (e > worst) ? e : worst;
    }
    EXPECT_LT(worst, 1.0e-9);

    // The occupations too - one file stated them over both spins and the other wrote Occup= 2, and
    // both land on the same per channel numbers.
    double occ_a[24] = {0};
    double occ_b[24] = {0};
    ASSERT_EQ(qm_test_series(occ_a, ARRAY_SIZE(occ_a), &tx, STR_LIT("orbital/alpha/occupation")), ARRAY_SIZE(occ_a));
    ASSERT_EQ(qm_test_series(occ_b, ARRAY_SIZE(occ_b), &md, STR_LIT("orbital/alpha/occupation")), ARRAY_SIZE(occ_b));
    for (size_t i = 0; i < ARRAY_SIZE(occ_a); ++i) {
        EXPECT_NEAR(occ_a[i], occ_b[i], 1.0e-12);
    }

    qm_test_free(&tx);
    qm_test_free(&md);
}

// What the sniffer is for: .h5 is shared with several other formats, so the extension does not
// settle it and a caller has to look inside.
UTEST(trexio, file_recognition) {
    EXPECT_TRUE(md_trexio_file_is_trexio(STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio")));
    EXPECT_FALSE(md_trexio_file_is_trexio(STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb")));
    EXPECT_FALSE(md_trexio_file_is_trexio(STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/no-such-file.trexio")));
}

UTEST(trexio, a_file_that_is_not_hdf5_is_declined) {
    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(1));
    EXPECT_FALSE(md_trexio_system_init_from_file(&t.sys, &t.state, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    qm_test_free(&t);
}

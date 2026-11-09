#include "utest.h"

#include "qm_test_util.h"
#include "trexio_pyscf_fixtures.h"

#include <md_trexio.h>
#include <md_molden.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <stdio.h>

// Everything a TREXIO file carries reaches a consumer as an md_system_t and its ATTRIBUTE TABLE:
// md_trexio.h has two entry points and no reader object to ask questions of. So these tests read
// what a real consumer reads, by the same paths.
//
// See test_data/trexio/README.md for where the file came from.

static bool trexio_load(qm_test_t* t, const char* path, size_t arena_bytes) {
    qm_test_init(t, arena_bytes);
    return md_trexio_system_init_from_file(&t->sys, &t->state, md_str_from_cstr(path));
}

UTEST(trexio, parse) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MD_MEGABYTES(16)));

    EXPECT_EQ(3u, t.sys.atom.count);
    EXPECT_EQ(3u, t.state.num_atoms);

    // TREXIO stores coordinates in bohr; a system's are in Angstrom, and so is what is published.
    const md_attribute_t* coord = qm_test_attr(&t, MD_STR_LIT("qm/atom/coordinate"));
    ASSERT_TRUE(coord != NULL);
    ASSERT_EQ(coord->format.components, 3u);
    ASSERT_EQ(md_attribute_value_count(&coord->format), 3u);

    double xyz[9] = {0};
    ASSERT_EQ(md_attribute_extract_f64(xyz, MD_COUNTOF(xyz), coord, md_attribute_slice_all(), md_unit_none()), MD_COUNTOF(xyz));
    EXPECT_NEAR(4.999999640091311, xyz[0], 1.0e-9);
    EXPECT_NEAR(7.147076485540978, xyz[1], 1.0e-9);
    EXPECT_NEAR(7.650970449269813, xyz[2], 1.0e-9);

    // ...and they are the same three numbers the Molden file beside it was written from.
    EXPECT_NEAR(4.068065707173540, xyz[3], 1.0e-9);

    EXPECT_NEAR(xyz[0], (double)t.state.xyz[0].x, 1.0e-4);

    // nucleus_label is authoritative for the element; nucleus_charge is what the electrons see and
    // is not the atomic number for an effective core potential.
    double z[3] = {0};
    const md_attribute_t* za = qm_test_attr(&t, MD_STR_LIT("qm/atom/atomic_number"));
    ASSERT_TRUE(za != NULL);
    ASSERT_EQ(md_attribute_extract_f64(z, 3, za, md_attribute_slice_all(), md_unit_none()), 3u);
    EXPECT_EQ(8.0, z[0]);
    EXPECT_EQ(1.0, z[1]);
    EXPECT_EQ(1.0, z[2]);

    // ...and nucleus_charge is published beside it as the nuclear charge. All electron here, so the two agree.
    double q[3] = {0};
    const md_attribute_t* qa = qm_test_attr(&t, MD_STR_LIT("qm/atom/nuclear_charge"));
    ASSERT_TRUE(qa != NULL);
    ASSERT_EQ(md_attribute_extract_f64(q, 3, qa, md_attribute_slice_all(), md_unit_elementary_charge()), 3u);
    EXPECT_EQ(8.0, q[0]);
    EXPECT_EQ(1.0, q[1]);
    EXPECT_EQ(1.0, q[2]);

    EXPECT_EQ(10.0, qm_test_scalar(&t, MD_STR_LIT("trexio/electron_count/total"), -1.0));
    EXPECT_EQ( 5.0, qm_test_scalar(&t, MD_STR_LIT("trexio/electron_count/up"),    -1.0));
    EXPECT_EQ( 5.0, qm_test_scalar(&t, MD_STR_LIT("trexio/electron_count/down"),  -1.0));

    EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("trexio/basis_type")),      MD_STR_LIT("Gaussian")));
    EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("trexio/mo_type")),         MD_STR_LIT("MO")));
    EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("trexio/ao_convention")),   MD_STR_LIT("spherical")));
    EXPECT_TRUE(md_str_eq(qm_test_string(&t, MD_STR_LIT("trexio/metadata/package_version")), MD_STR_LIT("2.6.1")));

    // This file carries no mo_energy and no nucleus_repulsion. An absent path is how a consumer
    // learns a block is missing; a column of zeros would read as a legitimate set of energies.
    EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("orbital/alpha/energy")));
    EXPECT_FALSE(qm_test_has(&t, MD_STR_LIT("trexio/nuclear_repulsion_energy")));

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
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MD_MEGABYTES(16)));

    double file_occ[24] = {0};
    ASSERT_EQ(qm_test_series(file_occ, MD_COUNTOF(file_occ), &t, MD_STR_LIT("trexio/mo/occupation")), MD_COUNTOF(file_occ));
    EXPECT_NEAR(2.0, file_occ[0], 1.0e-12);
    EXPECT_NEAR(2.0, file_occ[4], 1.0e-12);
    EXPECT_NEAR(0.0, file_occ[5], 1.0e-12);

    double occ[24] = {0};
    ASSERT_EQ(qm_test_series(occ, MD_COUNTOF(occ), &t, MD_STR_LIT("orbital/alpha/occupation")), MD_COUNTOF(occ));
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR(1.0, occ[i], 1.0e-12);
    }
    EXPECT_NEAR(0.0, occ[5], 1.0e-12);

    // Beta is a second NAME for alpha here, not a second copy.
    const md_attribute_t* a = qm_test_attr(&t, MD_STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* b = qm_test_attr(&t, MD_STR_LIT("orbital/beta/coefficient"));
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
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MD_MEGABYTES(32)));
    EXPECT_LT(qm_test_orthonormality(&t, MD_STR_LIT("orbital/alpha/coefficient")), 1.0e-5);
    qm_test_free(&t);
}

UTEST(trexio, total_density_holds_ten_electrons) {
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MD_MEGABYTES(32)));
    EXPECT_NEAR(10.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/total/density")), 1.0e-5);
    EXPECT_NEAR( 5.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/alpha/density")), 1.0e-5);
    EXPECT_NEAR( 0.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/difference/density")), 1.0e-9);
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
    ASSERT_TRUE(trexio_load(&tx, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio", MD_MEGABYTES(32)));
    qm_test_init(&md, MD_MEGABYTES(32));
    ASSERT_TRUE(md_molden_system_init_from_file(&md.sys, &md.state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));

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
    double* ca = qm_test_matrix(&tx, MD_STR_LIT("orbital/alpha/coefficient"), &na);
    double* cb = qm_test_matrix(&md, MD_STR_LIT("orbital/alpha/coefficient"), &nb);
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
    ASSERT_EQ(qm_test_series(occ_a, MD_COUNTOF(occ_a), &tx, MD_STR_LIT("orbital/alpha/occupation")), MD_COUNTOF(occ_a));
    ASSERT_EQ(qm_test_series(occ_b, MD_COUNTOF(occ_b), &md, MD_STR_LIT("orbital/alpha/occupation")), MD_COUNTOF(occ_b));
    for (size_t i = 0; i < MD_COUNTOF(occ_a); ++i) {
        EXPECT_NEAR(occ_a[i], occ_b[i], 1.0e-12);
    }

    qm_test_free(&tx);
    qm_test_free(&md);
}

// What the sniffer is for: .h5 is shared with several other formats, so the extension does not
// settle it and a caller has to look inside.
UTEST(trexio, file_recognition) {
    EXPECT_TRUE(md_trexio_file_is_trexio(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio")));
    EXPECT_FALSE(md_trexio_file_is_trexio(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb")));
    EXPECT_FALSE(md_trexio_file_is_trexio(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/no-such-file.trexio")));
}

UTEST(trexio, a_file_that_is_not_hdf5_is_declined) {
    qm_test_t t = {0};
    qm_test_init(&t, MD_MEGABYTES(1));
    EXPECT_FALSE(md_trexio_system_init_from_file(&t.sys, &t.state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// Files as PySCF writes them (trexio_pyscf_fixtures.h)
// ---------------------------------------------------------------------------

static bool tx_write(const char* path, const tx_item_t* items, size_t count) {
    const hid_t file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    bool ok = true;
    for (size_t i = 0; ok && i < count; ++i) {
        const tx_item_t* it = &items[i];
        const hid_t g = H5Lexists(file, it->group, H5P_DEFAULT) > 0 ? H5Gopen2(file, it->group, H5P_DEFAULT) : H5Gcreate2(file, it->group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (g < 0) { ok = false; break; }
        if (it->kind == TX_ATTR_I64) {
            const hid_t space = H5Screate(H5S_SCALAR);
            const hid_t a = H5Acreate2(g, it->name, H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT);
            ok = a >= 0 && H5Awrite(a, H5T_NATIVE_INT64, &it->value) >= 0;
            if (a >= 0) H5Aclose(a);
            H5Sclose(space);
        } else {
            const bool f64 = it->kind == TX_DSET_F64;
            const hid_t space = H5Screate_simple(it->rank, it->dims, NULL);
            const hid_t d = H5Dcreate2(g, it->name, f64 ? H5T_IEEE_F64LE : H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ok = d >= 0 && H5Dwrite(d, f64 ? H5T_NATIVE_DOUBLE : H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, it->data) >= 0;
            if (d >= 0) H5Dclose(d);
            H5Sclose(space);
        }
        H5Gclose(g);
    }
    H5Fclose(file);
    return ok;
}

// The AOs as the formula states them. Each AO's norm is ao_normalization x its radial part x its angular one, which
// for PySCF's xx is sqrt(4 pi / 5), and for its s functions 1 only with ao_normalization's 1/sqrt(4 pi) in it: read
// as unit normalised with ao_normalization on top (as before), the orbitals of this file were off by up to a factor
// sqrt(4 pi), and its density held the wrong number of electrons.
UTEST(trexio, pyscf_cartesian_functions_as_stated) {
    const char* path = "md_unittest_pyscf_h2o_cart.trexio";
    ASSERT_TRUE(tx_write(path, h2o_631gs_cart, MD_COUNTOF(h2o_631gs_cart)));
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, path, MD_MEGABYTES(32)));
    EXPECT_LT(qm_test_orthonormality(&t, MD_STR_LIT("orbital/alpha/coefficient")), 1.0e-5);
    EXPECT_NEAR(10.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/total/density")), 1.0e-5);
    qm_test_free(&t);
    remove(path);
}

// A restricted open shell: one set of orbitals, occupations 2 2 2 2 1 0 and no mo_spin to say more. The singly
// occupied orbital holds an alpha electron - halved into both channels, the radical had no spin density at all.
UTEST(trexio, restricted_open_shell) {
    const char* path = "md_unittest_pyscf_oh_rohf.trexio";
    ASSERT_TRUE(tx_write(path, oh_sto3g_rohf, MD_COUNTOF(oh_sto3g_rohf)));
    qm_test_t t = {0};
    ASSERT_TRUE(trexio_load(&t, path, MD_MEGABYTES(32)));
    EXPECT_LT(qm_test_orthonormality(&t, MD_STR_LIT("orbital/alpha/coefficient")), 1.0e-5);
    double occ_a[6] = {0}, occ_b[6] = {0};
    ASSERT_EQ(6u, (uint32_t)qm_test_series(occ_a, 6, &t, MD_STR_LIT("orbital/alpha/occupation")));
    ASSERT_EQ(6u, (uint32_t)qm_test_series(occ_b, 6, &t, MD_STR_LIT("orbital/beta/occupation")));
    EXPECT_EQ(1.0, occ_a[4]);
    EXPECT_EQ(0.0, occ_b[4]);
    EXPECT_EQ(1.0, occ_b[3]);
    EXPECT_NEAR(9.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/total/density")), 1.0e-5);
    EXPECT_NEAR(5.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/alpha/density")), 1.0e-5);
    EXPECT_NEAR(1.0, qm_test_electron_count(&t, MD_STR_LIT("orbital/difference/density")), 1.0e-5);
    qm_test_free(&t);
    remove(path);
}

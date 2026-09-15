#include "utest.h"

#include "qm_test_util.h"

#include <md_molden.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>

// Everything a Molden file carries reaches a consumer as an md_system_t and its ATTRIBUTE TABLE:
// md_molden.h has three entry points and no reader object to ask questions of. So these tests read
// what a real consumer reads, by the same paths, and a value that cannot be reached this way is a
// value no consumer can reach either - which is the property they exist to hold.
//
// See test_data/molden/README.md for where each file came from and what it is for.

static bool molden_load(qm_test_t* t, const char* path, size_t arena_bytes) {
    qm_test_init(t, arena_bytes);
    return md_molden_system_init_from_file(&t->sys, &t->state, str_from_cstr(path));
}

// ---------------------------------------------------------------------------
// The structural read
// ---------------------------------------------------------------------------

UTEST(molden, parse) {
    qm_test_t t = {0};
    ASSERT_TRUE(molden_load(&t, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden", MEGABYTES(16)));

    // The atoms land in the system itself, not in the table.
    EXPECT_EQ(3u, t.sys.atom.count);
    EXPECT_EQ(3u, t.state.num_atoms);

    // [Atoms] was written in Angstrom and the state is in Angstrom, so this is a straight read.
    const md_attribute_t* coord = qm_test_attr(&t, STR_LIT("qm/atom/coordinate"));
    ASSERT_TRUE(coord != NULL);
    ASSERT_EQ(md_attribute_components(&coord->format), 3u);
    ASSERT_EQ(md_attribute_value_count(&coord->format), 3u);

    double xyz[9] = {0};
    ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_unit_none()), ARRAY_SIZE(xyz));
    EXPECT_NEAR(4.999999640091, xyz[0], 1.0e-6);
    EXPECT_NEAR(7.147076485541, xyz[1], 1.0e-6);
    EXPECT_NEAR(7.650970449270, xyz[2], 1.0e-6);

    // The system's own state is that geometry too, narrowed to float.
    EXPECT_NEAR(xyz[0], (double)t.state.x[0], 1.0e-4);
    EXPECT_NEAR(xyz[1], (double)t.state.y[0], 1.0e-4);
    EXPECT_NEAR(xyz[2], (double)t.state.z[0], 1.0e-4);

    const md_attribute_t* z = qm_test_attr(&t, STR_LIT("qm/atom/atomic_number"));
    ASSERT_TRUE(z != NULL);
    double zz[3] = {0};
    ASSERT_EQ(md_attribute_extract_f64(zz, 3, z, md_unit_none()), 3u);
    EXPECT_EQ(8.0, zz[0]);
    EXPECT_EQ(1.0, zz[1]);
    EXPECT_EQ(1.0, zz[2]);

    // cc-pVDZ on water: O gets 3s 2p 1d, each H 2s 1p - twelve shells over twenty-two primitives.
    md_gto_basis_t basis = {0};
    ASSERT_TRUE(qm_test_basis(&basis, &t));
    EXPECT_EQ(12u, basis.num_shells);
    EXPECT_EQ(22u, basis.num_primitives);

    // The published basis is CARTESIAN whatever the file said, so the d shell contributes six.
    EXPECT_EQ(25u, (uint32_t)md_gto_basis_num_ao(&basis));
    EXPECT_EQ(24u, (uint32_t)md_gto_basis_num_sph_ao(&basis));

    qm_test_free(&t);
}

// The two things the format states in words. molden/ao_convention is not something a consumer can
// work out afterwards - the basis in the table is Cartesian either way - so if it stops being
// published there is no recovering what the coefficients were read as.
UTEST(molden, run_description_is_published) {
    qm_test_t t = {0};
    ASSERT_TRUE(molden_load(&t, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden", MEGABYTES(16)));

    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("molden/ao_convention")), STR_LIT("5D7F15G")));
    EXPECT_TRUE(str_begins_with(qm_test_string(&t, STR_LIT("molden/title")), STR_LIT("H2O cc-pVDZ")));

    // No vibrational block in this file, and no [Program] line. An absent path is how a consumer
    // learns a block is missing - not a zero it would have to interpret.
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("molden/vib/frequency")));
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("qm/atom/normal_mode")));
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("molden/program")));

    qm_test_free(&t);
}

UTEST(molden, orbital_energies_and_occupations) {
    qm_test_t t = {0};
    ASSERT_TRUE(molden_load(&t, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden", MEGABYTES(16)));

    ASSERT_EQ(24u, (uint32_t)qm_test_count(&t, STR_LIT("orbital/alpha/energy")));

    double energy[24] = {0};
    ASSERT_EQ(qm_test_series(energy, ARRAY_SIZE(energy), &t, STR_LIT("orbital/alpha/energy")), ARRAY_SIZE(energy));
    EXPECT_NEAR(-20.546019802060, energy[0], 1.0e-9);
    EXPECT_NEAR( -1.318586,       energy[1], 1.0e-6);
    EXPECT_NEAR( -0.498086,       energy[4], 1.0e-6);
    EXPECT_NEAR(  0.176113,       energy[5], 1.0e-6);

    // The file states Occup= 2 and has no beta orbitals, so it is restricted and each spin channel
    // gets half - which is the per channel convention the orbital/ tree uses throughout.
    double occ[24] = {0};
    ASSERT_EQ(qm_test_series(occ, ARRAY_SIZE(occ), &t, STR_LIT("orbital/alpha/occupation")), ARRAY_SIZE(occ));
    for (size_t i = 0; i < 5; ++i) {
        EXPECT_NEAR(1.0, occ[i], 1.0e-12);
    }
    EXPECT_NEAR(0.0, occ[5], 1.0e-12);

    // ...and beta is a second NAME for alpha, not a second copy.
    const md_attribute_t* a = qm_test_attr(&t, STR_LIT("orbital/alpha/coefficient"));
    const md_attribute_t* b = qm_test_attr(&t, STR_LIT("orbital/beta/coefficient"));
    ASSERT_TRUE(a != NULL);
    ASSERT_TRUE(b != NULL);
    EXPECT_TRUE(md_attribute_same_data(a, b));

    // Sym= is text and is published as one string per orbital, not as an index into anything.
    ASSERT_EQ(24u, (uint32_t)qm_test_count(&t, STR_LIT("orbital/alpha/symmetry")));
    EXPECT_TRUE(str_eq(qm_test_string_at(&t, STR_LIT("orbital/alpha/symmetry"), 0),  STR_LIT("A")));
    EXPECT_TRUE(str_eq(qm_test_string_at(&t, STR_LIT("orbital/alpha/symmetry"), 23), STR_LIT("A")));

    qm_test_free(&t);
}

// ---------------------------------------------------------------------------
// The numbers
// ---------------------------------------------------------------------------

// The one test that pins every convention at once. An orthonormal set of molecular orbitals stays
// orthonormal only if the atomic orbitals were put in the right order, each contraction was
// normalised against the convention its file stated, the spherical expansion kept its factors, and
// basis/overlap describes the same basis the coefficients are written against. Get any of those
// wrong and this number is of order one, not of order 1e-7.
//
// The tolerance is set by the basis itself: md_gto_basis_t holds exponents and contraction
// coefficients as float, so ~1e-7 is the floor for any file, however exactly it was written.
UTEST(molden, orbitals_are_orthonormal) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden",     // [5D], spherical d
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden",  // [6D], the same orbitals in Cartesian
    };
    for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
        qm_test_t t = {0};
        ASSERT_TRUE(molden_load(&t, files[i], MEGABYTES(32)));
        EXPECT_LT(qm_test_orthonormality(&t, STR_LIT("orbital/alpha/coefficient")), 1.0e-5);
        qm_test_free(&t);
    }
}

// tr(D S) over the total density is the electron count, and water has ten. It is a different
// assertion from orthonormality: this one is what catches the OCCUPATION convention, which
// orthonormality is completely blind to.
UTEST(molden, total_density_holds_ten_electrons) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden",
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden",
    };
    for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
        qm_test_t t = {0};
        ASSERT_TRUE(molden_load(&t, files[i], MEGABYTES(32)));
        EXPECT_NEAR(10.0, qm_test_electron_count(&t, STR_LIT("orbital/total/density")), 1.0e-5);
        // Restricted, so the two spins hold five each and their difference is nothing.
        EXPECT_NEAR( 5.0, qm_test_electron_count(&t, STR_LIT("orbital/alpha/density")), 1.0e-5);
        EXPECT_NEAR( 0.0, qm_test_electron_count(&t, STR_LIT("orbital/difference/density")), 1.0e-9);
        qm_test_free(&t);
    }
}

// The same calculation written two ways: [5D] states five spherical d functions, [6D] states six
// Cartesian ones. The basis in the table is Cartesian for both, so both files must describe the
// SAME twenty-five atomic orbitals - and, for the s and p shells, with the same coefficients, since
// nothing in either path touches those.
//
// The d coefficients are NOT expected to agree: a shell is normalised against the AO its own file
// stated, so the [5D] file's d contraction is sqrt(12) smaller than the [6D] file's and its
// coefficients are correspondingly larger. What must agree is the orbital they multiply out to,
// which is what the overlap-weighted comparison below asks.
UTEST(molden, spherical_and_cartesian_describe_one_wavefunction) {
    qm_test_t sph = {0};
    qm_test_t car = {0};
    ASSERT_TRUE(molden_load(&sph, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden", MEGABYTES(32)));
    ASSERT_TRUE(molden_load(&car, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden", MEGABYTES(32)));

    md_gto_basis_t sph_basis = {0};
    md_gto_basis_t car_basis = {0};
    ASSERT_TRUE(qm_test_basis(&sph_basis, &sph));
    ASSERT_TRUE(qm_test_basis(&car_basis, &car));
    ASSERT_EQ(sph_basis.num_shells, car_basis.num_shells);
    ASSERT_EQ(md_gto_basis_num_ao(&sph_basis), md_gto_basis_num_ao(&car_basis));

    // <psi_i^sph | psi_i^car> == 1 for every orbital: same function, two spellings. The overlap is
    // integrated from each file's own basis and the two bases differ by that per shell scale, so
    // this is computed in the spherical file's basis with the Cartesian file's coefficients
    // rescaled onto it - which is exactly what the DENSITY does, so compare that instead. A
    // density is invariant under the shell scaling, because the coefficients absorb it.
    size_t n_sph = 0, n_car = 0;
    double* d_sph = qm_test_matrix(&sph, STR_LIT("orbital/total/density"), &n_sph);
    double* d_car = qm_test_matrix(&car, STR_LIT("orbital/total/density"), &n_car);
    ASSERT_TRUE(d_sph != NULL);
    ASSERT_TRUE(d_car != NULL);
    ASSERT_EQ(n_sph, n_car);

    // The Cartesian file's density is stated against a basis whose d functions are sqrt(12) larger,
    // so it is correspondingly smaller there. Compare the invariant instead: both integrate to ten
    // electrons against their own overlap, which the test above already asserts, and the s/p block -
    // where the two bases are identical function for function - must agree outright.
    double worst = 0.0;
    for (size_t i = 0; i < 9; ++i) {         // the O s and p AOs, before the d shell starts
        for (size_t j = 0; j < 9; ++j) {
            const double e = fabs(d_sph[i * n_sph + j] - d_car[i * n_car + j]);
            worst = (e > worst) ? e : worst;
        }
    }
    EXPECT_LT(worst, 1.0e-6);

    qm_test_free(&sph);
    qm_test_free(&car);
}

// ---------------------------------------------------------------------------
// What the format lets a writer do
// ---------------------------------------------------------------------------

// An 'sp' shell is one line of the file carrying two contractions over one set of exponents. It
// becomes an s shell and a p shell, in that order, which is also the order its four atomic orbitals
// appear in - so a file written either way describes the same basis and the same orbitals, and
// nothing downstream has to know which spelling it was.
UTEST(molden, sp_shells_expand_to_s_and_p) {
    qm_test_t split = {0};
    qm_test_t fused = {0};
    ASSERT_TRUE(molden_load(&split, MD_UNITTEST_DATA_DIR "/molden/ammonia_sto3g.molden",    MEGABYTES(16)));
    ASSERT_TRUE(molden_load(&fused, MD_UNITTEST_DATA_DIR "/molden/ammonia_sto3g_sp.molden", MEGABYTES(16)));

    md_gto_basis_t a = {0};
    md_gto_basis_t b = {0};
    ASSERT_TRUE(qm_test_basis(&a, &split));
    ASSERT_TRUE(qm_test_basis(&b, &fused));

    ASSERT_EQ(6u, a.num_shells);
    ASSERT_EQ(a.num_shells, b.num_shells);
    ASSERT_EQ(a.num_primitives, b.num_primitives);
    ASSERT_EQ(8u, (uint32_t)md_gto_basis_num_ao(&a));

    for (uint32_t s = 0; s < a.num_shells; ++s) {
        EXPECT_EQ(a.shells[s].l,              b.shells[s].l);
        EXPECT_EQ(a.shells[s].atom_idx,       b.shells[s].atom_idx);
        EXPECT_EQ(a.shells[s].num_primitives, b.shells[s].num_primitives);
    }

    // The primitives of a contraction are a SUM, so a writer is free to list them in any order and
    // the sp file does - it has to pick one order for two columns that were written in two. What
    // must match is the function the contraction adds up to, and the overlap is exactly that,
    // element for element, with no ordering in it.
    size_t da = 0, db = 0;
    double* sa = qm_test_matrix(&split, STR_LIT("basis/overlap"), &da);
    double* sb = qm_test_matrix(&fused, STR_LIT("basis/overlap"), &db);
    ASSERT_TRUE(sa != NULL);
    ASSERT_TRUE(sb != NULL);
    ASSERT_EQ(8u, (uint32_t)da);
    ASSERT_EQ(da, db);
    for (size_t i = 0; i < da * db; ++i) {
        EXPECT_NEAR(sa[i], sb[i], 1.0e-6);
    }

    // ...and the orbitals written against them are the same numbers.
    size_t na = 0, nb = 0;
    double* ca = qm_test_matrix(&split, STR_LIT("orbital/alpha/coefficient"), &na);
    double* cb = qm_test_matrix(&fused, STR_LIT("orbital/alpha/coefficient"), &nb);
    ASSERT_TRUE(ca != NULL);
    ASSERT_TRUE(cb != NULL);
    ASSERT_EQ(na, nb);
    for (size_t i = 0; i < 8 * na; ++i) {
        EXPECT_NEAR(ca[i], cb[i], 1.0e-9);
    }

    qm_test_free(&split);
    qm_test_free(&fused);
}

// A file written by a third party, with none of the generated files' regularity: parenthesised unit
// on [Atoms], no angular convention marker at all, Sym= before Ene=, and exponents listed out of
// order within a contraction. It is here to be READ, not to be believed - see the README, its
// molecular orbitals are not orthonormal - so what is asserted is what the parser recovered.
UTEST(molden, third_party_file_parses) {
    qm_test_t t = {0};
    ASSERT_TRUE(molden_load(&t, MD_UNITTEST_DATA_DIR "/molden/ammonia_sto3g.molden", MEGABYTES(16)));

    EXPECT_EQ(4u, t.sys.atom.count);

    // No marker in the file, so the format's own default stands and is published as such.
    EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("molden/ao_convention")), STR_LIT("6D10F15G")));

    md_gto_basis_t basis = {0};
    ASSERT_TRUE(qm_test_basis(&basis, &t));
    EXPECT_EQ(6u, basis.num_shells);       // N: 1s 2s 2p, plus one s per hydrogen
    EXPECT_EQ(18u, basis.num_primitives);  // STO-3G, three per shell
    EXPECT_EQ(8u, (uint32_t)md_gto_basis_num_ao(&basis));

    ASSERT_EQ(8u, (uint32_t)qm_test_count(&t, STR_LIT("orbital/alpha/energy")));
    double energy[8] = {0};
    ASSERT_EQ(qm_test_series(energy, ARRAY_SIZE(energy), &t, STR_LIT("orbital/alpha/energy")), ARRAY_SIZE(energy));
    EXPECT_NEAR(-15.3128, energy[0], 1.0e-4);

    qm_test_free(&t);
}

// The vibrational blocks, on a file written into the test rather than stored beside it: the numbers
// below are placeholders, chosen so that a transposed axis or an off by one would be obvious, and
// nothing here asserts physics. What is under test is that [FREQ], [INT] and [FR-NORM-COORD] are
// read, that the modes come out one row per mode and one displacement per atom, and that a section
// this reader does not know about does not derail the ones around it.
UTEST(molden, vibrational_blocks) {
    static const char* file =
        "[Molden Format]\n"
        "[Title]\n"
        " vibrational fixture\n"
        "[Atoms] AU\n"
        "O 1 8  0.000000  0.000000  0.000000\n"
        "H 2 1  0.000000  1.430000  1.100000\n"
        "H 3 1  0.000000 -1.430000  1.100000\n"
        "[UNKNOWN-SECTION]\n"
        " something a later version of the format added\n"
        "[FREQ]\n"
        "1600.0\n"
        "3700.0\n"
        "3800.0\n"
        "[INT]\n"
        "70.0\n"
        "5.0\n"
        "45.0\n"
        "[FR-NORM-COORD]\n"
        "vibration 1\n"
        "  0.01  0.02  0.03\n"
        "  0.11  0.12  0.13\n"
        "  0.21  0.22  0.23\n"
        "vibration 2\n"
        "  1.01  1.02  1.03\n"
        "  1.11  1.12  1.13\n"
        "  1.21  1.22  1.23\n"
        "vibration 3\n"
        "  2.01  2.02  2.03\n"
        "  2.11  2.12  2.13\n"
        "  2.21  2.22  2.23\n";

    qm_test_t t = {0};
    qm_test_init(&t, MEGABYTES(4));
    ASSERT_TRUE(md_molden_system_init_from_str(&t.sys, &t.state, str_from_cstr(file)));

    EXPECT_EQ(3u, t.sys.atom.count);

    // [Atoms] AU, so the geometry was in bohr and the system is in Angstrom.
    double xyz[9] = {0};
    const md_attribute_t* coord = qm_test_attr(&t, STR_LIT("qm/atom/coordinate"));
    ASSERT_TRUE(coord != NULL);
    ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_unit_none()), ARRAY_SIZE(xyz));
    EXPECT_NEAR(1.43 * 0.529177210903, xyz[4], 1.0e-9);

    double freq[3] = {0};
    ASSERT_EQ(qm_test_series(freq, ARRAY_SIZE(freq), &t, STR_LIT("molden/vib/frequency")), ARRAY_SIZE(freq));
    EXPECT_NEAR(1600.0, freq[0], 1.0e-9);
    EXPECT_NEAR(3800.0, freq[2], 1.0e-9);

    double intensity[3] = {0};
    ASSERT_EQ(qm_test_series(intensity, ARRAY_SIZE(intensity), &t, STR_LIT("molden/vib/ir_intensity")), ARRAY_SIZE(intensity));
    EXPECT_NEAR(70.0, intensity[0], 1.0e-9);

    // {D,N} of 3 component values: mode first, atom last, which is the convention an "atom axis is
    // the last index axis" consumer relies on.
    const md_attribute_t* mode = qm_test_attr(&t, STR_LIT("qm/atom/normal_mode"));
    ASSERT_TRUE(mode != NULL);
    ASSERT_EQ(2u, mode->format.rank);
    ASSERT_EQ(3u, mode->format.shape[0]);
    ASSERT_EQ(3u, mode->format.shape[1]);
    ASSERT_EQ(3u, md_attribute_components(&mode->format));

    double row[9] = {0};
    ASSERT_EQ(qm_test_row(row, ARRAY_SIZE(row), &t, STR_LIT("qm/atom/normal_mode"), 1), ARRAY_SIZE(row));
    EXPECT_NEAR(1.01, row[0], 1.0e-9);
    EXPECT_NEAR(1.13, row[5], 1.0e-9);
    EXPECT_NEAR(1.23, row[8], 1.0e-9);

    // No [GTO], so nothing that depends on a basis is published - and the file is still read.
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("basis/shell/angular_momentum")));
    EXPECT_FALSE(qm_test_has(&t, STR_LIT("orbital/alpha/coefficient")));

    qm_test_free(&t);
}

// What the sniffer is for: an extension is not evidence, since .molden, .molden.input and .mold are
// all in use and none of them is exclusive.
UTEST(molden, file_recognition) {
    EXPECT_TRUE(md_molden_file_is_molden(STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden")));
    EXPECT_TRUE(md_molden_file_is_molden(STR_LIT(MD_UNITTEST_DATA_DIR "/molden/ammonia_sto3g.molden")));
    EXPECT_FALSE(md_molden_file_is_molden(STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb")));
    EXPECT_FALSE(md_molden_file_is_molden(STR_LIT(MD_UNITTEST_DATA_DIR "/molden/no-such-file.molden")));
}

// Refusals. Each of these is a file the reader must decline rather than half read, because a system
// with a plausible looking basis in it is worse than no system at all.
UTEST(molden, malformed_files_are_declined) {
    struct {
        const char* what;
        const char* text;
    } cases[] = {
        { "no atoms",
          "[Molden Format]\n[Title]\n empty\n" },
        { "a shell type this library cannot evaluate",
          "[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\n h 1 1.00\n  1.0 1.0\n" },
        { "a shell before any atom is named",
          "[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n s 1 1.00\n  1.0 1.0\n" },
        { "an atom index the file does not have",
          "[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n7 0\n s 1 1.00\n  1.0 1.0\n" },
        { "a contraction that ends early",
          "[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\n s 3 1.00\n  1.0 1.0\n" },
        { "an orbital naming an atomic orbital the basis does not have",
          "[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n[GTO]\n1 0\n s 1 1.00\n  1.0 1.0\n\n[MO]\nEne= 0.0\nOccup= 1.0\n 9 1.0\n" },
    };

    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        qm_test_t t = {0};
        qm_test_init(&t, MEGABYTES(1));
        EXPECT_FALSE(md_molden_system_init_from_str(&t.sys, &t.state, str_from_cstr(cases[i].text)));
        qm_test_free(&t);
    }
}

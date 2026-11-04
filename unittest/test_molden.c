#include "utest.h"

#include "qm_test_util.h"

#include <md_molden.h>
#include <md_qm.h>
#include <md_gto.h>
#include <md_gto_iso.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

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
    // A QM calculation alone: every atom is in its QM region
    for (size_t i = 0; i < t.sys.atom.count; ++i) {
        EXPECT_TRUE((md_atom_flags(&t.sys.atom, i) & MD_ATOM_FLAG_QM) != 0);
    }

    // [Atoms] was written in Angstrom and the state is in Angstrom, so this is a straight read.
    const md_attribute_t* coord = qm_test_attr(&t, STR_LIT("qm/atom/coordinate"));
    ASSERT_TRUE(coord != NULL);
    ASSERT_EQ(coord->format.components, 3u);
    ASSERT_EQ(md_attribute_value_count(&coord->format), 3u);

    double xyz[9] = {0};
    ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_attribute_slice_all(), md_unit_none()), ARRAY_SIZE(xyz));
    EXPECT_NEAR(4.999999640091, xyz[0], 1.0e-6);
    EXPECT_NEAR(7.147076485541, xyz[1], 1.0e-6);
    EXPECT_NEAR(7.650970449270, xyz[2], 1.0e-6);

    // The system's own state is that geometry too, narrowed to float.
    EXPECT_NEAR(xyz[0], (double)t.state.xyz[0].x, 1.0e-4);
    EXPECT_NEAR(xyz[1], (double)t.state.xyz[0].y, 1.0e-4);
    EXPECT_NEAR(xyz[2], (double)t.state.xyz[0].z, 1.0e-4);

    const md_attribute_t* z = qm_test_attr(&t, STR_LIT("qm/atom/atomic_number"));
    ASSERT_TRUE(z != NULL);
    double zz[3] = {0};
    ASSERT_EQ(md_attribute_extract_f64(zz, 3, z, md_attribute_slice_all(), md_unit_none()), 3u);
    EXPECT_EQ(8.0, zz[0]);
    EXPECT_EQ(1.0, zz[1]);
    EXPECT_EQ(1.0, zz[2]);

    // Molden states no effective core charges, so the nuclear charges are the atomic numbers - and
    // they are published all the same, so that a consumer never has to choose between the two
    const md_attribute_t* nq = qm_test_attr(&t, STR_LIT("qm/atom/nuclear_charge"));
    ASSERT_TRUE(nq != NULL);
    double qq[3] = {0};
    ASSERT_EQ(md_attribute_extract_f64(qq, 3, nq, md_attribute_slice_all(), md_unit_elementary_charge()), 3u);
    EXPECT_EQ(8.0, qq[0]);
    EXPECT_EQ(1.0, qq[1]);
    EXPECT_EQ(1.0, qq[2]);

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
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart.molden",// [6D10F15G] from PySCF: d, f and g
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d_shell.molden",   // the two above with their Cartesian
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart_shell.molden", // functions normalised as a shell
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
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart.molden",
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d_shell.molden",
        MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart_shell.molden",
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

// The electron density of a loaded file at n points (bohr), from its own basis and orbitals. 'atom' is the
// water's three atoms in bohr, read from the first file a comparison loads.
static bool molden_density_at(double* rho, const float* p, size_t n, qm_test_t* t, const float atom[9]) {
    md_gto_basis_t basis = {0};
    if (!qm_test_basis(&basis, t)) return false;
    double *rows = NULL, *w = NULL;
    size_t na = 0;
    const size_t k = md_qm_density_orbitals(&rows, &w, &na, &t->sys, STR_LIT("orbital/total/density"), md_get_heap_allocator());
    md_gto_iso_field_t* f = k ? md_gto_iso_field_create_products(&basis, atom, 0, rows, NULL, w, k, md_get_heap_allocator()) : NULL;
    if (f) {
        md_gto_iso_field_evaluate(f, p, n, 0, rho, NULL);
        md_gto_iso_field_destroy(f);
    }
    if (rows) md_free(md_get_heap_allocator(), rows, sizeof(double) * k * na);
    if (w)    md_free(md_get_heap_allocator(), w, sizeof(double) * k);
    return f != NULL;
}

// The largest relative difference between two files' densities at points around the water (> 1e-6 only)
static double molden_density_difference(const char* file_a, const char* file_b) {
    qm_test_t a = {0}, b = {0};
    double worst = INFINITY;
    if (molden_load(&a, file_a, MEGABYTES(32)) && molden_load(&b, file_b, MEGABYTES(32))) {
        const md_attribute_t* ca = qm_test_attr(&a, STR_LIT("qm/atom/coordinate"));
        double xyz[9] = {0};
        if (ca && md_attribute_extract_f64(xyz, 9, ca, md_attribute_slice_all(), md_unit_none()) == 9) {
            float atom[9];
            for (int i = 0; i < 9; ++i) atom[i] = (float)(xyz[i] * 1.8897261246257702);
            enum { N = 2000 };
            float* p = (float*)malloc(sizeof(float) * 3 * N);
            double* ra = (double*)malloc(sizeof(double) * N);
            double* rb = (double*)malloc(sizeof(double) * N);
            uint64_t r = 12345;
            for (int i = 0; i < 3 * N; ++i) {
                r = r * 6364136223846793005ull + 1442695040888963407ull;
                p[i] = atom[i % 3] + 6.0f * ((float)(r >> 40) / (float)(1ull << 24) - 0.5f);
            }
            if (molden_density_at(ra, p, N, &a, atom) && molden_density_at(rb, p, N, &b, atom)) {
                worst = 0.0;
                for (int i = 0; i < N; ++i) {
                    if (ra[i] > 1e-6) worst = fmax(worst, fabs(rb[i] - ra[i]) / ra[i]);
                }
            }
            free(p); free(ra); free(rb);
        }
    }
    qm_test_free(&a);
    qm_test_free(&b);
    return worst;
}

// The same calculation written two ways: [5D] states five spherical d functions, [6D] six Cartesian ones, each
// normalised on its own as the format has them. The basis in the table is Cartesian for both, so both files
// describe the same twenty-five atomic orbitals - and the same density everywhere, d functions included.
//
// The coefficients are not compared: a shell is normalised against the AO its own file stated, so the two
// d contractions differ by a constant the coefficients absorb. The density does not depend on that, so it
// is compared at points: the electron density of each file, evaluated from its own basis and orbitals. Over
// the s and p shells alone this would have passed with the Cartesian coefficients misread, as it did while
// the reader took them against functions sharing the xx normalisation (PySCF: off by 3.6% in the density).
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
    qm_test_free(&sph);
    qm_test_free(&car);

    const double worst = molden_density_difference(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden", MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden");
    printf("  [5D] and [6D] water: the densities at 2000 points agree to %.1e (relative)\n", worst);
    EXPECT_LT(worst, 1.0e-5);
}

// Two normalisations of Cartesian functions are written: each on its own (the format's: Molden, PySCF), and all
// of a shell with its (l,0,0) monomial's constant (Turbomole, CFOUR 2.1, PSI4 <= 1.3.2), the files not saying
// which. The reader tells them apart by the orbitals, which are normalised in only one of them, and says which it
// took (molden/cartesian_normalisation). Each _shell file is its counterpart with its coefficients restated in
// the shell convention - the cc-pVDZ water over d functions, the cc-pVQZ one over d, f and g - so read right, the
// two give one density.
UTEST(molden, cartesian_normalisation_is_read_from_the_orbitals) {
    static const struct { const char* file; const char* expect; } cases[] = {
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden",            NULL },        // spherical d: nothing to decide
        { MD_UNITTEST_DATA_DIR "/molden/ammonia_sto3g.molden",         NULL },        // no d at all
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden",         "function" },
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d_shell.molden",   "shell" },
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart.molden",       "function" },
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart_shell.molden", "shell" },
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        qm_test_t t = {0};
        ASSERT_TRUE(molden_load(&t, cases[i].file, MEGABYTES(32)));
        const bool has = qm_test_has(&t, STR_LIT("molden/cartesian_normalisation"));
        if (cases[i].expect) {
            EXPECT_TRUE(has);
            EXPECT_TRUE(str_eq(qm_test_string(&t, STR_LIT("molden/cartesian_normalisation")), str_from_cstr(cases[i].expect)));
        } else {
            EXPECT_FALSE(has);
        }
        qm_test_free(&t);
    }
    static const char* pairs[][2] = {
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden",   MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d_shell.molden" },
        { MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart.molden", MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvqz_cart_shell.molden" },
    };
    for (size_t i = 0; i < ARRAY_SIZE(pairs); ++i) {
        const double worst = molden_density_difference(pairs[i][0], pairs[i][1]);
        printf("  %s: each function and shell normalised read to one density, %.1e (relative)\n", pairs[i][0] + strlen(MD_UNITTEST_DATA_DIR "/molden/"), worst);
        EXPECT_LT(worst, 1.0e-9);
    }
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
    ASSERT_EQ(md_attribute_extract_f64(xyz, ARRAY_SIZE(xyz), coord, md_attribute_slice_all(), md_unit_none()), ARRAY_SIZE(xyz));
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
    ASSERT_EQ(3u, mode->format.components);

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

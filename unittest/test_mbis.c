#include "utest.h"

#include "qm_test_util.h"
#include "test_executor.h"

#include <md_mbis.h>
#include <md_molden.h>
#include <md_qm.h>
#include <md_gto.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

// MBIS against an INDEPENDENT implementation: the reference values below come from a numpy MBIS (the fixed
// point of Verstraelen et al. 2016, same initial guess) on PySCF's own molecular grid (Treutler radial,
// Lebedev angular, Becke partition) and PySCF's own density evaluation from the same Molden files, converged
// in the grid (levels 3, 5 and 7 agree to 1e-6): tools/mbis_reference.py. It shares nothing with md_mbis.c but the
// definition.

typedef struct mbis_test_input_t {
    qm_test_t       t;
    md_gto_basis_t  basis;
    dvec3_t         xyz[16];
    uint8_t         z[16];
    double          q[16];
    size_t          num_atoms;
    double*         D;
} mbis_test_input_t;

static bool mbis_test_load(mbis_test_input_t* in, const char* path) {
    memset(in, 0, sizeof(*in));
    qm_test_init(&in->t, MD_MEGABYTES(16));
    if (!md_molden_system_init_from_file(&in->t.sys, &in->t.state, md_str_from_cstr(path))) return false;
    if (!qm_test_basis(&in->basis, &in->t)) return false;

    const md_attribute_t* ca = qm_test_attr(&in->t, MD_STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* za = qm_test_attr(&in->t, MD_STR_LIT("qm/atom/atomic_number"));
    const md_attribute_t* qa = qm_test_attr(&in->t, MD_STR_LIT("qm/atom/nuclear_charge"));
    if (!ca || !za || !qa) return false;
    in->num_atoms = md_attribute_value_count(&ca->format);
    if (in->num_atoms > MD_COUNTOF(in->xyz)) return false;

    double buf[3 * 16];
    if (md_attribute_extract_f64(buf, 3 * in->num_atoms, ca, md_attribute_slice_all(), md_unit_bohr_radius()) != 3 * in->num_atoms) return false;
    for (size_t i = 0; i < in->num_atoms; ++i) in->xyz[i] = (dvec3_t){{{buf[3 * i], buf[3 * i + 1], buf[3 * i + 2]}}};
    if (md_attribute_extract_f64(buf, in->num_atoms, za, md_attribute_slice_all(), md_unit_none()) != in->num_atoms) return false;
    for (size_t i = 0; i < in->num_atoms; ++i) in->z[i] = (uint8_t)buf[i];
    if (md_attribute_extract_f64(in->q, in->num_atoms, qa, md_attribute_slice_all(), md_unit_elementary_charge()) != in->num_atoms) return false;

    size_t dim = 0;
    in->D = qm_test_matrix(&in->t, MD_STR_LIT("orbital/total/density"), &dim);
    return in->D && dim == md_gto_basis_num_ao(&in->basis);
}

static md_mbis_desc_t mbis_test_desc(const mbis_test_input_t* in) {
    return (md_mbis_desc_t){
        .basis = &in->basis, .atom_coord = in->xyz, .atomic_number = in->z, .nuclear_charge = in->q,
        .num_atoms = in->num_atoms, .density_matrix = in->D,
    };
}

UTEST(mbis, water_against_an_independent_implementation) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    const md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t r;
    ASSERT_TRUE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));

    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(10.0, r.electron_count, 1.0e-5);
    // Default grid: within 2e-5 e of the converged reference
    EXPECT_NEAR(-0.869937, r.charge[0], 2.0e-5);
    EXPECT_NEAR( 0.434968, r.charge[1], 2.0e-5);
    EXPECT_NEAR( 0.434968, r.charge[2], 2.0e-5);
    EXPECT_NEAR(r.charge[1], r.charge[2], 1.0e-6);              // the two H are equivalent
    EXPECT_NEAR(0.0, r.charge[0] + r.charge[1] + r.charge[2], 1.0e-5);

    // The pro-atoms: O has a core and a valence shell, each H one
    ASSERT_EQ(4u, r.num_shells);
    ASSERT_EQ(0u, r.shell_offset[0]);
    ASSERT_EQ(2u, r.shell_offset[1]);
    ASSERT_EQ(3u, r.shell_offset[2]);
    ASSERT_EQ(4u, r.shell_offset[3]);
    EXPECT_NEAR(1.616555, r.shell_population[0], 1.0e-4);
    EXPECT_NEAR(0.056098, r.shell_width[0],      1.0e-5);
    EXPECT_NEAR(7.253382, r.shell_population[1], 1.0e-4);
    EXPECT_NEAR(0.403099, r.shell_width[1],      1.0e-5);
    EXPECT_NEAR(0.565032, r.shell_population[2], 1.0e-4);
    EXPECT_NEAR(0.346771, r.shell_width[2],      1.0e-5);

    md_mbis_result_free(&r, md_get_heap_allocator());
    qm_test_free(&in.t);
}

// A finer grid converges ONTO the reference, rather than merely landing near it with the defaults.
UTEST(mbis, converges_with_the_grid) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    md_mbis_desc_t desc = mbis_test_desc(&in);
    desc.radial_points = 100;
    desc.angular_order = 47;
    desc.threshold = 1.0e-9;   // the grid is what is tested: the iterations converged well beyond it
    md_mbis_result_t r;
    ASSERT_TRUE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    EXPECT_NEAR(10.0, r.electron_count, 1.0e-6);
    EXPECT_NEAR(-0.869937, r.charge[0], 2.0e-6);
    EXPECT_NEAR( 0.434968, r.charge[1], 2.0e-6);
    EXPECT_NEAR( 0.434968, r.charge[2], 2.0e-6);
    md_mbis_result_free(&r, md_get_heap_allocator());
    qm_test_free(&in.t);
}

// The same calculation stored with Cartesian (6D) and with spherical (5D) d functions: one density, so
// one partition. A partition of the BASIS (Mulliken) would agree too; one of space must.
UTEST(mbis, independent_of_the_angular_convention_of_the_file) {
    mbis_test_input_t a, b;
    ASSERT_TRUE(mbis_test_load(&a, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    ASSERT_TRUE(mbis_test_load(&b, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz_6d.molden"));
    const md_mbis_desc_t da = mbis_test_desc(&a), db = mbis_test_desc(&b);
    md_mbis_result_t ra, rb;
    ASSERT_TRUE(md_mbis_compute(&ra, &da, md_get_heap_allocator()));
    ASSERT_TRUE(md_mbis_compute(&rb, &db, md_get_heap_allocator()));
    for (size_t i = 0; i < 3; ++i) EXPECT_NEAR(ra.charge[i], rb.charge[i], 1.0e-6);
    md_mbis_result_free(&ra, md_get_heap_allocator());
    md_mbis_result_free(&rb, md_get_heap_allocator());
    qm_test_free(&a.t);
    qm_test_free(&b.t);
}

// An unrestricted radical, partitioned through its TOTAL density (alpha + beta).
UTEST(mbis, open_shell) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/oh_631gs_uhf.molden"));
    const md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t r;
    ASSERT_TRUE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    EXPECT_TRUE(r.converged);
    EXPECT_NEAR(9.0, r.electron_count, 1.0e-5);
    EXPECT_NEAR(-0.413525, r.charge[0], 2.0e-5);
    EXPECT_NEAR( 0.413525, r.charge[1], 2.0e-5);
    md_mbis_result_free(&r, md_get_heap_allocator());
    qm_test_free(&in.t);
}

// A lone atom owns all of its density, and its one shell's width is a third of the mean distance of its
// electrons from the nucleus.
UTEST(mbis, lone_atom) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/he_ccpvdz.molden"));
    const md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t r;
    ASSERT_TRUE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    EXPECT_TRUE(r.converged);
    ASSERT_EQ(1u, r.num_shells);
    EXPECT_NEAR(0.0, r.charge[0], 1.0e-6);
    EXPECT_NEAR(0.307786, r.shell_width[0], 1.0e-5);
    md_mbis_result_free(&r, md_get_heap_allocator());
    qm_test_free(&in.t);
}

// Threads, and pieces in any order, give the same bits.
UTEST(mbis, executor_does_not_change_the_result) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t ref;
    ASSERT_TRUE(md_mbis_compute(&ref, &desc, md_get_heap_allocator()));

    test_executor_t threads = test_executor_threads(4);
    test_executor_t shuffled = test_executor_shuffled(3);
    const md_executor_t* executors[2] = {test_executor(&threads), test_executor(&shuffled)};
    for (size_t e = 0; e < 2; ++e) {
        desc.executor = executors[e];
        md_mbis_result_t r;
        ASSERT_TRUE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
        EXPECT_EQ(ref.iterations, r.iterations);
        EXPECT_EQ(ref.num_points, r.num_points);
        EXPECT_TRUE(memcmp(ref.charge, r.charge, sizeof(double) * ref.num_atoms) == 0);
        EXPECT_TRUE(memcmp(ref.shell_width, r.shell_width, sizeof(double) * ref.num_shells) == 0);
        md_mbis_result_free(&r, md_get_heap_allocator());
    }
    md_mbis_result_free(&ref, md_get_heap_allocator());
    qm_test_free(&in.t);
}

// Started from its own answer, a run is converged at once and lands where it started; the pro-atoms of other atoms
// are refused rather than used.
UTEST(mbis, warm_start) {
    mbis_test_input_t in, other;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    ASSERT_TRUE(mbis_test_load(&other, MD_UNITTEST_DATA_DIR "/molden/oh_631gs_uhf.molden"));
    md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t cold, warm, foreign;
    ASSERT_TRUE(md_mbis_compute(&cold, &desc, md_get_heap_allocator()));
    desc.initial = &cold;
    ASSERT_TRUE(md_mbis_compute(&warm, &desc, md_get_heap_allocator()));
    EXPECT_TRUE(warm.converged);
    EXPECT_LE(warm.iterations, 3u);
    EXPECT_LT(warm.iterations, cold.iterations);
    for (size_t i = 0; i < 3; ++i) EXPECT_NEAR(cold.charge[i], warm.charge[i], 2.0e-6);

    const md_mbis_desc_t od = mbis_test_desc(&other);
    ASSERT_TRUE(md_mbis_compute(&foreign, &od, md_get_heap_allocator()));
    desc.initial = &foreign;   // OH: two atoms, three shells
    md_mbis_result_t r;
    EXPECT_FALSE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));

    md_mbis_result_free(&cold, md_get_heap_allocator());
    md_mbis_result_free(&warm, md_get_heap_allocator());
    md_mbis_result_free(&foreign, md_get_heap_allocator());
    qm_test_free(&in.t);
    qm_test_free(&other.t);
}

UTEST(mbis, refuses_incomplete_input_and_stops_when_cancelled) {
    mbis_test_input_t in;
    ASSERT_TRUE(mbis_test_load(&in, MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
    md_mbis_desc_t desc = mbis_test_desc(&in);
    md_mbis_result_t r;
    desc.density_matrix = NULL;
    EXPECT_FALSE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    desc = mbis_test_desc(&in);
    desc.num_atoms = 2;     // the basis has shells on a third atom
    EXPECT_FALSE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    desc = mbis_test_desc(&in);
    volatile int32_t cancel = 1;
    desc.cancel = &cancel;
    EXPECT_FALSE(md_mbis_compute(&r, &desc, md_get_heap_allocator()));
    EXPECT_TRUE(r.charge == NULL);
    qm_test_free(&in.t);
}

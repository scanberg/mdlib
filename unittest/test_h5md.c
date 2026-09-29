#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_h5md.h>
#include <md_tpr.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_str.h>

#include "run_check.h"

// See test_data/h5md/README.md for where the files come from. The GROMACS ones are mdrun's own
// output beside the run input it was given, so a system read from the trajectory can be held
// against one read from the tpr. The spec_*.h5md files are written by make_spec_files.py, and every
// value in them follows from one formula, repeated here:
//
//     position[f][i][k] = 0.1 f + 0.01 i + 0.001 k
#define H5MD_DIR MD_UNITTEST_DATA_DIR "/h5md/"

static double spec_position(size_t f, size_t i, size_t k) {
    return 0.1 * (double)f + 0.01 * (double)i + 0.001 * (double)k;
}

typedef struct h5md_test_t {
    md_allocator_i*   arena;
    md_system_t       sys;
    md_system_state_t state;
} h5md_test_t;

static bool h5md_test_load(h5md_test_t* t, const char* path) {
    MEMSET(t, 0, sizeof(*t));
    t->arena = md_vm_arena_create(GIGABYTES(1));
    t->sys.alloc = t->arena;
    t->state.alloc = t->arena;
    return md_h5md_system_init_from_file(&t->sys, &t->state, str_from_cstr(path));
}

static void h5md_test_free(h5md_test_t* t) {
    md_vm_arena_destroy(t->arena);
}

static const md_attribute_t* attr(const md_system_t* sys, const char* path) {
    return md_attributes_find(&sys->attributes, str_from_cstr(path));
}

// The values of a resident or virtual attribute's slice as doubles, in its own unit
static size_t slice_f64(double* dst, size_t cap, const md_attribute_t* a, md_attribute_slice_t slice) {
    return a ? md_attribute_extract_slice_f64(dst, cap, a, &slice, md_unit_none()) : 0;
}

static bool has_bond(const md_system_t* sys, int a, int b) {
    return md_bond_find(&sys->bond, a, b) != -1;
}

// ### GROMACS ###

// The frames of peptide_tip3p.h5md, read with h5py: sums over all atoms, the first and last atom,
// and the cell, Angstrom. The pressure coupling moves the box from frame to frame.
static const run_ref_t tip3p_refs[] = {
    { 0, {19338.455118, 19465.512666, 19690.742607}, {10.58959, 11.21027, 13.18069}, {6.07, 10.97, 22.5},           {25.237999, 25.237999, 25.237999, 0, 0, 0} },
    { 2, {19324.492904, 19449.853706, 19686.644035}, {10.59674, 11.18981, 13.15262}, {6.044874, 10.96786, 22.45431}, {25.2230239, 25.2230239, 25.2230239, 0, 0, 0} },
    { 4, {19325.086883, 19444.218226, 19687.764425}, {10.59322, 11.13981, 13.16189}, {6.018226, 10.9787, 22.38792}, {25.2180123, 25.2180123, 25.2180123, 0, 0, 0} },
};

static const run_ref_t tip4p_refs[] = {
    { 0, {26161.935300, 25148.294790, 24941.768123}, {10.43946, 11.24066, 13.07944}, {20.25, 17.52, 1.84},           {25.237999, 25.237999, 25.237999, 0, 0, 0} },
    { 2, {26102.907444, 25082.000287, 24881.991101}, {10.43022, 11.2131, 13.00079},  {20.21934, 17.48667, 1.828081}, {25.1761699, 25.1761699, 25.1761699, 0, 0, 0} },
    { 4, {26155.918094, 25121.442569, 24927.620054}, {10.45516, 11.21008, 12.97113}, {20.24551, 17.5236, 1.844893},  {25.2195954, 25.2195954, 25.2195954, 0, 0, 0} },
};

// The same topology read from the trajectory and from the run input it was simulated from has to
// be the same system: the atoms, their types, the residues and their numbers, the charges and the
// bonds. What the GROMACS module leaves out (force field types, non-bonded parameters) is the only
// difference allowed.
static void compare_with_tpr(int* utest_result, const char* h5md_path, const char* tpr_path, bool vsites) {
    h5md_test_t h = {0};
    ASSERT_TRUE(h5md_test_load(&h, h5md_path));

    md_system_t tpr = { .alloc = h.arena };
    md_system_state_t tpr_state = { .alloc = h.arena };
    ASSERT_TRUE(md_tpr_system_init_from_file(&tpr, &tpr_state, str_from_cstr(tpr_path)));

    ASSERT_EQ(tpr.atom.count, h.sys.atom.count);
    EXPECT_TRUE(str_eq(tpr.description, h.sys.description));
    for (size_t i = 0; i < tpr.atom.count; ++i) {
        EXPECT_TRUE(str_eq(md_atom_name(&tpr.atom, i), md_atom_name(&h.sys.atom, i)));
        EXPECT_EQ(md_atom_atomic_number(&tpr.atom, i), md_atom_atomic_number(&h.sys.atom, i));
        EXPECT_EQ(md_atom_mass(&tpr.atom, i), md_atom_mass(&h.sys.atom, i));
        EXPECT_EQ(md_atom_flags(&tpr.atom, i) & MD_FLAG_COARSE_GRAINED, md_atom_flags(&h.sys.atom, i) & MD_FLAG_COARSE_GRAINED);
    }

    ASSERT_EQ(tpr.component.count, h.sys.component.count);
    for (size_t c = 0; c < tpr.component.count; ++c) {
        EXPECT_TRUE(str_eq(md_component_name(&tpr.component, c), md_component_name(&h.sys.component, c)));
        EXPECT_EQ(md_component_seq_id(&tpr.component, c), md_component_seq_id(&h.sys.component, c));
        EXPECT_EQ(md_component_atom_range(&tpr.component, c).beg, md_component_atom_range(&h.sys.component, c).beg);
    }
    // The dipeptide keeps the numbers it had in the PDB it was built from
    EXPECT_EQ(221, md_component_seq_id(&h.sys.component, 0));
    EXPECT_EQ(222, md_component_seq_id(&h.sys.component, 1));

    const md_attribute_t* qa = attr(&tpr, "atom/charge");
    const md_attribute_t* qb = attr(&h.sys, "atom/charge");
    ASSERT_TRUE(qa && qb);
    ASSERT_EQ(md_attribute_element_count(&qa->format), md_attribute_element_count(&qb->format));
    EXPECT_EQ(0, memcmp(qa->data, qb->data, md_attribute_byte_size(&qa->format)));

    // GROMACS writes the chemical bonds, constraints and SETTLE. The tpr reader also ties each
    // virtual site to its first constructing atom, which the trajectory does not have.
    size_t missing = 0;
    for (size_t b = 0; b < h.sys.bond.count; ++b) {
        EXPECT_TRUE(h.sys.bond.flags[b] & MD_BOND_FLAG_TOPOLOGY);
        missing += md_bond_find(&tpr.bond, h.sys.bond.pairs[b].idx[0], h.sys.bond.pairs[b].idx[1]) == -1;
    }
    EXPECT_EQ(0u, missing);
    if (vsites) {
        EXPECT_LT(h.sys.bond.count, tpr.bond.count);
    } else {
        EXPECT_EQ(tpr.bond.count, h.sys.bond.count);
    }

    // The first frame is where the run started, which is what the tpr holds - after mdrun has
    // constrained it and put the atoms in the box, so the same up to a box vector and a little
    ASSERT_EQ(tpr_state.num_atoms, h.state.num_atoms);
    const double L = tpr_state.unitcell.x;
    double max_diff = 0;
    for (size_t i = 0; i < h.state.num_atoms; ++i) {
        const vec3_t d = vec3_sub(tpr_state.xyz[i], h.state.xyz[i]);
        const double dx = d.x - L * round(d.x / L), dy = d.y - L * round(d.y / L), dz = d.z - L * round(d.z / L);
        max_diff = MAX(max_diff, sqrt(dx * dx + dy * dy + dz * dz));
    }
    EXPECT_LT(max_diff, 0.05);
    EXPECT_NEAR(tpr_state.unitcell.x, h.state.unitcell.x, 1.0e-4);
    EXPECT_TRUE(md_unitcell_flags(&h.state.unitcell) & MD_UNITCELL_PBC_Z);

    const md_attribute_t* creator = attr(&h.sys, "h5md/creator/name");
    ASSERT_TRUE(creator != NULL);
    EXPECT_TRUE(str_eq(md_attribute_str(&h.sys.attributes, creator, 0), STR_LIT("GROMACS")));

    h5md_test_free(&h);
}

UTEST(h5md, gromacs_tip3p_system_matches_tpr) {
    compare_with_tpr(utest_result, H5MD_DIR "peptide_tip3p.h5md", H5MD_DIR "peptide_tip3p.tpr", false);
}

UTEST(h5md, gromacs_tip4p_system_matches_tpr) {
    compare_with_tpr(utest_result, H5MD_DIR "peptide_tip4p.h5md", H5MD_DIR "peptide_tip4p.tpr", true);
}

// Positions every 5 steps, velocities every 10 and forces every 20: the positions make the run, and
// velocities and forces each get a frame axis of their own, so each value keeps the step it was
// written at.
static void check_gromacs_run(int* utest_result, const char* path, const run_ref_t* refs, size_t num_refs, size_t num_atoms, const double vel_sum[3][3]) {
    h5md_test_t t = {0};
    ASSERT_TRUE(h5md_test_load(&t, path));
    const str_t run = STR_LIT("run/gmx");
    ASSERT_TRUE(md_h5md_system_publish_run(&t.sys, str_from_cstr(path), run, MD_RUN_FLAG_NONE));
    run_check_refs(utest_result, &t.sys, run, 5, num_atoms, refs, num_refs);

    const md_attribute_t* time = attr(&t.sys, "run/gmx/time");
    const md_attribute_t* step = attr(&t.sys, "run/gmx/step");
    ASSERT_TRUE(time && step);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_NEAR(0.04, ((const double*)time->data)[4], 1.0e-12);
    EXPECT_EQ(15, ((const int64_t*)step->data)[3]);

    // Every frame is one uncompressed chunk, so each has a place in the file
    const md_attribute_t* offset = attr(&t.sys, "run/gmx/source/offset");
    ASSERT_TRUE(offset != NULL);
    for (size_t f = 0; f < 5; ++f) EXPECT_GT(((const int64_t*)offset->data)[f], 0);

    EXPECT_TRUE(attr(&t.sys, "run/gmx/atom/velocity") == NULL);
    EXPECT_TRUE(attr(&t.sys, "run/gmx/atom/force") == NULL);
    const md_attribute_t* vel = attr(&t.sys, "run/gmx/h5md/particles/velocity/atom/velocity");
    const md_attribute_t* frc = attr(&t.sys, "run/gmx/h5md/particles/force/atom/force");
    ASSERT_TRUE(vel && frc);
    EXPECT_EQ(3u, vel->format.shape[0]);
    EXPECT_EQ(2u, frc->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(vel->unit, md_unit_div(md_unit_nanometer(), md_unit_picosecond())));
    EXPECT_TRUE(md_unit_equal(frc->unit, md_unit_div(md_unit_div(md_unit_scl(md_unit_joule(), 1.0e3), md_unit_mole()), md_unit_nanometer())));
    EXPECT_EQ(attr(&t.sys, "run/gmx/h5md/particles/velocity/time"), md_attributes_axis(&t.sys.attributes, vel));

    // Through extraction: velocities at frames 0 and 2 (steps 0 and 10), none at frame 1 (step 5)
    md_system_state_t st = { .alloc = t.arena };
    md_system_state_init(&st, num_atoms);
    const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("h5md/particles/velocity/atom/velocity") };
    md_system_extract_t* ex = md_system_extract_begin(&t.sys, run, paths, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    for (int64_t f = 0; f < 5; ++f) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &st));
        const md_attribute_t* v = md_attributes_find(&st.attributes, STR_LIT("h5md/particles/velocity/atom/velocity"));
        if (f % 2) {
            EXPECT_TRUE(v == NULL);
            continue;
        }
        ASSERT_TRUE(v && v->data);
        const float* d = (const float*)v->data;
        double sum[3] = {0};
        for (size_t i = 0; i < num_atoms; ++i) for (int k = 0; k < 3; ++k) sum[k] += d[i * 3 + k];
        for (int k = 0; k < 3; ++k) EXPECT_NEAR(vel_sum[f / 2][k], sum[k], 1.0e-3);
    }
    md_system_extract_end(ex);

    // Removing the run removes all of it
    md_attributes_remove_prefix(&t.sys.attributes, run);
    EXPECT_EQ(0u, md_attributes_query(NULL, 0, &t.sys.attributes, run));
    h5md_test_free(&t);
}

UTEST(h5md, gromacs_tip3p_run) {
    const double vel_sum[3][3] = {
        {-56.37175536, -27.86827456, 57.61431187},
        { 29.34552022, -21.48351553, 26.86924008},
        { -2.28961452,   8.71451763, 13.51708333},
    };
    check_gromacs_run(utest_result, H5MD_DIR "peptide_tip3p.h5md", tip3p_refs, ARRAY_SIZE(tip3p_refs), 1577, vel_sum);
}

UTEST(h5md, gromacs_tip4p_run) {
    const double vel_sum[3][3] = {
        {-10.72599577, -31.5178012,  15.92132909},
        { 39.64799054, -14.61398636,  6.68716171},
        { -2.16551102, -26.56838597,  4.46162295},
    };
    check_gromacs_run(utest_result, H5MD_DIR "peptide_tip4p.h5md", tip4p_refs, ARRAY_SIZE(tip4p_refs), 2074, vel_sum);
}

// ### THE SPECIFICATION ###

UTEST(h5md, spec_system) {
    h5md_test_t t = {0};
    ASSERT_TRUE(h5md_test_load(&t, H5MD_DIR "spec_fixed.h5md"));

    // The larger of the two particle groups
    ASSERT_EQ(4u, t.sys.atom.count);
    ASSERT_EQ(4u, t.state.num_atoms);

    // An enumerated species names the particles, and the name gives the element
    const char* names[] = { "C", "O", "H", "H" };
    const md_atomic_number_t z[] = { 6, 8, 1, 1 };
    const float mass[] = { 12.011f, 15.999f, 1.008f, 1.008f };
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_TRUE(str_eq_cstr(md_atom_name(&t.sys.atom, i), names[i]));
        EXPECT_EQ(z[i], md_atom_atomic_number(&t.sys.atom, i));
        EXPECT_NEAR(mass[i], md_atom_mass(&t.sys.atom, i), 1.0e-5);
    }
    EXPECT_EQ(md_atom_type_idx(&t.sys.atom, 2), md_atom_type_idx(&t.sys.atom, 3));

    // Pairs name particles by id. The fill value, a pair naming no particle, pairs of the other
    // group and triples are not bonds
    EXPECT_EQ(3u, t.sys.bond.count);
    EXPECT_TRUE(has_bond(&t.sys, 0, 1));
    EXPECT_TRUE(has_bond(&t.sys, 1, 2));
    EXPECT_TRUE(has_bond(&t.sys, 1, 3));

    double q[4];
    ASSERT_EQ(4u, slice_f64(q, 4, attr(&t.sys, "atom/charge"), md_attribute_slice_all()));
    EXPECT_NEAR(-0.8, q[1], 1.0e-6);

    // The first frame, nm to Angstrom
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_NEAR(10.0 * spec_position(0, i, 0), t.state.xyz[i].x, 1.0e-5);
        EXPECT_NEAR(10.0 * spec_position(0, i, 2), t.state.xyz[i].z, 1.0e-5);
    }
    // A cuboid box, and z without a boundary is not periodic
    EXPECT_NEAR(20.0, t.state.unitcell.x, 1.0e-9);
    EXPECT_NEAR(30.0, t.state.unitcell.y, 1.0e-9);
    EXPECT_NEAR(0.0,  t.state.unitcell.z, 1.0e-9);
    EXPECT_TRUE(md_unitcell_flags(&t.state.unitcell) & MD_UNITCELL_PBC_X);
    EXPECT_FALSE(md_unitcell_flags(&t.state.unitcell) & MD_UNITCELL_PBC_Z);

    const md_attribute_t* version = attr(&t.sys, "h5md/version");
    ASSERT_TRUE(version != NULL);
    EXPECT_EQ(1, ((const int32_t*)version->data)[0]);
    EXPECT_TRUE(str_eq(md_attribute_str(&t.sys.attributes, attr(&t.sys, "h5md/author/name"), 0), STR_LIT("mdlib unittest")));
    h5md_test_free(&t);
}

UTEST(h5md, spec_run) {
    h5md_test_t t = {0};
    ASSERT_TRUE(h5md_test_load(&t, H5MD_DIR "spec_fixed.h5md"));
    const str_t run = STR_LIT("run/spec");
    ASSERT_TRUE(md_h5md_system_publish_run(&t.sys, STR_LIT(H5MD_DIR "spec_fixed.h5md"), run, MD_RUN_FLAG_NONE));

    // A fixed step and time: an increment and an offset
    const md_attribute_t* time = attr(&t.sys, "run/spec/time");
    const md_attribute_t* step = attr(&t.sys, "run/spec/step");
    ASSERT_TRUE(time && step);
    ASSERT_EQ(6u, time->format.shape[0]);
    for (size_t f = 0; f < 6; ++f) {
        EXPECT_NEAR(0.2 + 0.02 * (double)f, ((const double*)time->data)[f], 1.0e-12);
        EXPECT_EQ(100 + 10 * (int64_t)f, ((const int64_t*)step->data)[f]);
    }

    // Big endian doubles, read straight from the file; then the same through a state
    const md_attribute_t* pos = attr(&t.sys, "run/spec/atom/position");
    ASSERT_TRUE(pos != NULL);
    EXPECT_GT(((const int64_t*)attr(&t.sys, "run/spec/source/offset")->data)[0], 0);
    double x[12];
    for (size_t f = 0; f < 6; ++f) {
        ASSERT_EQ(12u, slice_f64(x, 12, pos, md_attribute_slice_1((uint32_t)f)));
        for (size_t i = 0; i < 4; ++i) for (size_t k = 0; k < 3; ++k) {
            EXPECT_NEAR(10.0 * spec_position(f, i, k), x[i * 3 + k], 1.0e-5);
        }
    }
    ASSERT_EQ(3u, slice_f64(x, 3, pos, md_attribute_slice_2(5, 3)));
    EXPECT_NEAR(10.0 * spec_position(5, 3, 1), x[1], 1.0e-5);

    md_system_state_t st = { .alloc = t.arena };
    md_system_state_init(&st, 4);
    ASSERT_TRUE(run_extract_one(&st, &t.sys, run, 4));
    EXPECT_NEAR(10.0 * spec_position(4, 2, 1), st.xyz[2].y, 1.0e-5);
    EXPECT_NEAR(20.0, st.unitcell.x, 1.0e-6);
    EXPECT_NEAR(0.0,  st.unitcell.z, 1.0e-6);

    // Force at every position step: beside the positions. It is compressed, so HDF5 reads it
    const md_attribute_t* force = attr(&t.sys, "run/spec/atom/force");
    ASSERT_TRUE(force != NULL);
    EXPECT_EQ(-1, ((const int64_t*)attr(&t.sys, "run/spec/h5md/particles/force/source/offset")->data)[0]);
    ASSERT_EQ(12u, slice_f64(x, 12, force, md_attribute_slice_1(3)));
    EXPECT_NEAR(-spec_position(3, 1, 2), x[5], 1.0e-6);

    // Velocity at steps 100, 120, 140 and no time of its own: the run's time step gives it one.
    // Two frames to a chunk, without a filter, so still read straight from the file
    const md_attribute_t* vel = attr(&t.sys, "run/spec/h5md/particles/velocity/atom/velocity");
    const md_attribute_t* vel_time = attr(&t.sys, "run/spec/h5md/particles/velocity/time");
    ASSERT_TRUE(vel && vel_time);
    EXPECT_TRUE(md_unit_equal(vel_time->unit, md_unit_picosecond()));
    EXPECT_NEAR(0.24, ((const double*)vel_time->data)[1], 1.0e-12);
    EXPECT_GT(((const int64_t*)attr(&t.sys, "run/spec/h5md/particles/velocity/source/offset")->data)[2], 0);
    ASSERT_EQ(12u, slice_f64(x, 12, vel, md_attribute_slice_1(2)));
    EXPECT_NEAR(2.0 * spec_position(2, 3, 0), x[9], 1.0e-6);
    size_t row = 0;
    EXPECT_TRUE(md_attribute_axis_map(&row, time, 4, vel_time));
    EXPECT_EQ(2u, row);
    EXPECT_FALSE(md_attribute_axis_map(&row, time, 3, vel_time));

    // Observables, as the file has them
    const md_attribute_t* energy = attr(&t.sys, "run/spec/h5md/observables/potential_energy/value");
    ASSERT_TRUE(energy != NULL);
    EXPECT_TRUE(energy->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
    double e[3];
    ASSERT_EQ(2u, slice_f64(e, 2, energy, md_attribute_slice_all()));
    EXPECT_EQ(-11.25, e[1]);
    EXPECT_TRUE(md_attributes_axis(&t.sys.attributes, energy) == attr(&t.sys, "run/spec/h5md/observables/potential_energy/time"));

    const md_attribute_t* com = attr(&t.sys, "run/spec/h5md/observables/center/of_mass/value");
    ASSERT_TRUE(com != NULL);
    EXPECT_EQ(2u, com->format.rank);
    ASSERT_EQ(3u, slice_f64(e, 3, com, md_attribute_slice_1(5)));
    EXPECT_EQ(15.0, e[0]);
    EXPECT_NEAR(0.3, ((const double*)attr(&t.sys, "run/spec/h5md/observables/center/of_mass/time")->data)[5], 1.0e-12);

    const md_attribute_t* volume = attr(&t.sys, "run/spec/h5md/observables/volume/value");
    ASSERT_TRUE(volume != NULL);
    EXPECT_EQ(0u, volume->format.rank);
    double vol = 0;
    ASSERT_EQ(1u, md_attribute_extract_f64(&vol, 1, volume, md_unit_pow(md_unit_angstrom(), 3)));
    EXPECT_NEAR(24000.0, vol, 1.0e-6);
    EXPECT_TRUE(attr(&t.sys, "run/spec/h5md/observables/label/value") == NULL);

    h5md_test_free(&t);
}

// No units, no time, integer species, compressed positions and a box per frame that misses one
UTEST(h5md, spec_bare) {
    h5md_test_t t = {0};
    ASSERT_TRUE(h5md_test_load(&t, H5MD_DIR "spec_bare.h5md"));
    ASSERT_EQ(4u, t.sys.atom.count);
    EXPECT_TRUE(str_eq_cstr(md_atom_name(&t.sys.atom, 2), "8"));
    EXPECT_EQ(8, md_atom_atomic_number(&t.sys.atom, 2));    // from its mass
    EXPECT_EQ(1, md_atom_atomic_number(&t.sys.atom, 3));
    EXPECT_EQ(0u, t.sys.bond.count);
    // Without a unit, a length is Angstrom
    EXPECT_NEAR(spec_position(0, 3, 2), t.state.xyz[3].z, 1.0e-6);
    EXPECT_TRUE(md_unitcell_is_triclinic(&t.state.unitcell));
    EXPECT_NEAR(5.0, t.state.unitcell.xy, 1.0e-6);

    const str_t run = STR_LIT("run/bare");
    ASSERT_TRUE(md_h5md_system_publish_run(&t.sys, STR_LIT(H5MD_DIR "spec_bare.h5md"), run, MD_RUN_FLAG_NONE));
    const md_attribute_t* time = attr(&t.sys, "run/bare/time");
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(3.0, ((const double*)time->data)[3]);
    EXPECT_EQ(-1, ((const int64_t*)attr(&t.sys, "run/bare/source/offset")->data)[0]);

    md_system_state_t st = { .alloc = t.arena };
    md_system_state_init(&st, 4);
    for (int64_t f = 0; f < 6; ++f) {
        ASSERT_TRUE(run_extract_one(&st, &t.sys, run, f));
        EXPECT_NEAR(spec_position((size_t)f, 1, 1), st.xyz[1].y, 1.0e-6);
        if (f == 2) {
            EXPECT_EQ(0u, md_unitcell_flags(&st.unitcell));
        } else {
            EXPECT_NEAR(20.0 * (1.0 + 0.1 * (double)f), st.unitcell.x, 1.0e-4);
        }
    }

    // Integer images, beside the positions they share their step with
    const md_attribute_t* image = attr(&t.sys, "run/bare/atom/image");
    ASSERT_TRUE(image != NULL);
    double img[12];
    ASSERT_EQ(12u, slice_f64(img, 12, image, md_attribute_slice_1(4)));
    EXPECT_EQ(1.0, img[7]);
    h5md_test_free(&t);
}

UTEST(h5md, declines) {
    h5md_test_t t = {0};
    // HDF5, but not H5MD
    EXPECT_FALSE(h5md_test_load(&t, MD_UNITTEST_DATA_DIR "/trexio/h2o.trexio"));
    h5md_test_free(&t);
    EXPECT_FALSE(h5md_test_load(&t, H5MD_DIR "does_not_exist.h5md"));
    h5md_test_free(&t);
    EXPECT_FALSE(h5md_test_load(&t, MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    h5md_test_free(&t);

    // A run for another system leaves nothing behind
    ASSERT_TRUE(h5md_test_load(&t, H5MD_DIR "spec_fixed.h5md"));
    EXPECT_FALSE(md_h5md_system_publish_run(&t.sys, STR_LIT(H5MD_DIR "peptide_tip3p.h5md"), STR_LIT("run/x"), MD_RUN_FLAG_NONE));
    EXPECT_EQ(0u, md_attributes_query(NULL, 0, &t.sys.attributes, STR_LIT("run/x")));
    h5md_test_free(&t);
}

// Frames read from several threads at once, each with its own extraction context, are the frames
// read from one: through the io straight from the file (the GROMACS positions), and through HDF5
// under the reader's lock (the compressed positions of spec_bare).
typedef struct h5md_thread_job_t {
    const md_system_t* sys;
    str_t   run;
    size_t  num_atoms;
    size_t  num_frames;
    size_t  rounds;
    double  sum;        // over every coordinate of every frame, in the first round
    size_t  mismatch;   // later rounds that summed to something else
    bool    ok;
} h5md_thread_job_t;

static void h5md_thread_extract(void* data) {
    h5md_thread_job_t* job = (h5md_thread_job_t*)data;
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_state_t st = { .alloc = arena };
    md_system_state_init(&st, job->num_atoms);
    const str_t paths[] = { STR_INIT("atom/position") };
    md_system_extract_t* ex = md_system_extract_begin(job->sys, job->run, paths, 1, arena);
    job->ok = ex != NULL;
    for (size_t r = 0; job->ok && r < job->rounds; ++r) {
        double sum = 0;
        for (size_t f = 0; job->ok && f < job->num_frames; ++f) {
            job->ok = md_system_extract_frame(ex, (int64_t)f, &st);
            for (size_t i = 0; job->ok && i < job->num_atoms; ++i) sum += st.xyz[i].x + st.xyz[i].y + st.xyz[i].z;
        }
        if (r == 0) job->sum = sum;
        else job->mismatch += sum != job->sum;
    }
    if (ex) md_system_extract_end(ex);
    md_vm_arena_destroy(arena);
}

static void check_threads(int* utest_result, const char* path) {
    h5md_test_t t = {0};
    ASSERT_TRUE(h5md_test_load(&t, path));
    const str_t run = STR_LIT("run/threads");
    ASSERT_TRUE(md_h5md_system_publish_run(&t.sys, str_from_cstr(path), run, MD_RUN_FLAG_NONE));

    h5md_thread_job_t one = { .sys = &t.sys, .run = run, .num_atoms = t.sys.atom.count, .num_frames = run_num_frames(&t.sys, run), .rounds = 1 };
    h5md_thread_extract(&one);
    ASSERT_TRUE(one.ok);

    enum { NUM_THREADS = 4, ROUNDS = 8 };
    h5md_thread_job_t jobs[NUM_THREADS];
    md_thread_t* threads[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; ++i) {
        jobs[i] = one;
        jobs[i].rounds = ROUNDS;
        jobs[i].sum = 0;
        threads[i] = md_thread_create(h5md_thread_extract, &jobs[i]);
    }
    for (int i = 0; i < NUM_THREADS; ++i) {
        md_thread_join(threads[i]);
        EXPECT_TRUE(jobs[i].ok);
        EXPECT_EQ(one.sum, jobs[i].sum);
        EXPECT_EQ(0u, jobs[i].mismatch);
    }
    h5md_test_free(&t);
}

UTEST(h5md, threads_read_from_the_file) {
    check_threads(utest_result, H5MD_DIR "peptide_tip3p.h5md");
}

UTEST(h5md, threads_read_through_hdf5) {
    check_threads(utest_result, H5MD_DIR "spec_bare.h5md");
}

#include "utest.h"
#include <string.h>
#include <stdio.h>

#include <md_trr.h>
#include <md_gro.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include "run_check.h"

#define TRR_RUN   STR_LIT("run/tryptophan")
#define TRR_FILE  STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr")

// Recorded from the trajectory reader this replaced
static const run_ref_t trr_refs[] = {
    { 0,   {127924.769, 127848.061, 128821.903}, {20.8102722, 26.3608513, 20.7103443}, {38.0198212, 1.62372589, 37.3288422}, {40, 40, 40, 0, 0, 0} },
    { 50,  {130043.601, 129768.828, 129645.515}, {19.9605198, 22.9012718, 22.1681995}, {37.9994278, 1.77006245, 1.10822511}, {40, 40, 40, 0, 0, 0} },
    { 100, {130022.237, 129490.678, 129577.395}, {14.8264885, 19.565136, 17.9483681},  {3.74613476, 35.8030586, 38.5597649}, {40, 40, 40, 0, 0, 0} },
};

// Velocities of the same frames: sums over all atoms and the first atom, Angstrom/ps
static const struct { size_t frame; double sum[3]; double first[3]; } trr_vel_refs[] = {
    { 0,   {636.961223, 435.094655, 237.041484},    {-5.90699244, 6.62543249, 7.37030125} },
    { 50,  {-187.007794, -866.681276, -932.603337}, {-0.289574832, -7.29495096, 17.1043167} },
    { 100, {-451.388558, -949.039574, -189.758381}, {-4.07094908, -6.2315712, -14.9569607} },
};

static bool trr_load(md_system_t* sys, md_allocator_i* arena) {
    sys->alloc = arena;
    md_system_state_t sys_state = {.alloc = arena};
    return md_gro_system_init_from_file(sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro")) &&
        md_trr_system_publish_run(sys, TRR_FILE, TRR_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE);
}

UTEST(trr, run_matches_reference) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {0};
    ASSERT_TRUE(trr_load(&sys, arena));
    run_check_refs(utest_result, &sys, TRR_RUN, 101, 6495, trr_refs, ARRAY_SIZE(trr_refs));

    const md_attributes_t* t = &sys.attributes;
    const md_attribute_t* time = md_attributes_find(t, STR_LIT("run/tryptophan/time"));
    const md_attribute_t* step = md_attributes_find(t, STR_LIT("run/tryptophan/step"));
    const md_attribute_t* vel  = md_attributes_find(t, STR_LIT("run/tryptophan/atom/velocity"));
    ASSERT_TRUE(time && step && vel);
    // No forces in the file, and velocities in every frame, so no group of their own
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/tryptophan/atom/force")) == NULL);
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/tryptophan/trr/velocity/time")) == NULL);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_EQ(time, md_attributes_axis(t, vel));
    EXPECT_EQ(1.0, ((const double*)time->data)[1]);

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, sys.atom.count);
    const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("atom/velocity") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, TRR_RUN, paths, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    for (size_t r = 0; r < ARRAY_SIZE(trr_vel_refs); ++r) {
        ASSERT_TRUE(md_system_extract_frame(ex, (int64_t)trr_vel_refs[r].frame, &st));
        const md_attribute_t* v = md_attributes_find(&st.attributes, STR_LIT("atom/velocity"));
        ASSERT_TRUE(v && v->data);
        const float* d = (const float*)v->data;
        double sum[3] = {0};
        for (size_t i = 0; i < sys.atom.count; ++i) for (int k = 0; k < 3; ++k) sum[k] += d[i * 3 + k];
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(trr_vel_refs[r].sum[k], sum[k], 1.0e-2);
            EXPECT_NEAR(trr_vel_refs[r].first[k], d[k], 1.0e-5);
        }
        EXPECT_TRUE(md_unit_equal(v->unit, md_unit_div(md_unit_angstrom(), md_unit_picosecond())));
    }
    md_system_extract_end(ex);

    // One atom of one frame, without a context: only its three values are read.
    float full[3], one_xyz[3];
    ASSERT_TRUE(run_extract_one(&st, &sys, TRR_RUN, 7));
    full[0] = st.xyz[11].x; full[1] = st.xyz[11].y; full[2] = st.xyz[11].z;
    const md_attribute_t* pos = md_attributes_find(t, STR_LIT("run/tryptophan/atom/position"));
    md_attribute_slice_t one = md_attribute_slice_2(7, 11);
    ASSERT_EQ(3u, md_attribute_extract_slice_f32(one_xyz, 3, pos, &one, md_unit_none()));
    EXPECT_EQ(0, MEMCMP(full, one_xyz, sizeof(full)));

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(trr, nonexistent_file) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    EXPECT_FALSE(md_trr_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.trr"), TRR_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(0u, run_num_frames(&sys, TRR_RUN));
    md_vm_arena_destroy(arena);
}

UTEST(trr, frame_boundary_conditions) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {0};
    ASSERT_TRUE(trr_load(&sys, arena));
    const int64_t F = (int64_t)run_num_frames(&sys, TRR_RUN);

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, sys.atom.count);
    EXPECT_FALSE(run_extract_one(&st, &sys, TRR_RUN, -1));
    EXPECT_FALSE(run_extract_one(&st, &sys, TRR_RUN, F));
    EXPECT_TRUE(run_extract_one(&st, &sys, TRR_RUN, 0));
    EXPECT_TRUE(run_extract_one(&st, &sys, TRR_RUN, F - 1));

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(trr, run_refuses_a_different_system) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro")));
    EXPECT_FALSE(md_trr_system_publish_run(&sys, TRR_FILE, TRR_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(0u, run_num_frames(&sys, TRR_RUN));
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A TRR written by hand, big endian as XDR is, so each section can be present or not per frame.
static void trr_put_i32(md_file_t file, int32_t v) {
    const uint32_t u = (uint32_t)v;
    const uint8_t b[4] = { (uint8_t)(u >> 24), (uint8_t)(u >> 16), (uint8_t)(u >> 8), (uint8_t)u };
    md_file_write(file, b, 4);
}

static void trr_put_f32(md_file_t file, float v) {
    int32_t i;
    MEMCPY(&i, &v, 4);
    trr_put_i32(file, i);
}

static void trr_put_frame(md_file_t file, int natoms, int step, float time, const float box[9], const float* x, const float* v, const float* f) {
    const int vec = natoms * 3 * 4;
    trr_put_i32(file, 1993);
    trr_put_i32(file, 13);
    trr_put_i32(file, 12);
    md_file_write(file, "GMX_trn_file", 12);
    const int32_t fields[13] = { 0, 0, box ? 36 : 0, 0, 0, 0, 0, x ? vec : 0, v ? vec : 0, f ? vec : 0, natoms, step, 0 };
    for (int i = 0; i < 13; ++i) trr_put_i32(file, fields[i]);
    trr_put_f32(file, time);
    trr_put_f32(file, 0.0f);
    if (box) for (int i = 0; i < 9; ++i) trr_put_f32(file, box[i]);
    if (x) for (int i = 0; i < natoms * 3; ++i) trr_put_f32(file, x[i]);
    if (v) for (int i = 0; i < natoms * 3; ++i) trr_put_f32(file, v[i]);
    if (f) for (int i = 0; i < natoms * 3; ++i) trr_put_f32(file, f[i]);
}

// Velocities written at every other frame and in one frame without coordinates, forces in every
// frame with coordinates, a triclinic box, and a last frame cut short as a run still being written
// leaves it. The run's frames are the complete frames with coordinates; the forces sit beside the
// positions and the velocities get a group of their own, each at the time it was written.
UTEST(trr, run_publishes_sections_at_their_own_times) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_LIT("md_unittest_trr_sections.trr");
    enum { NA = 5 };

    // nm, row i box vector i, tilted as GROMACS writes it: below the diagonal
    const float box[9] = { 3.0f, 0.0f, 0.0f,  1.0f, 2.5f, 0.0f,  -0.5f, 0.75f, 2.0f };
    float x[NA * 3], v[NA * 3], f[NA * 3];
    for (int i = 0; i < NA * 3; ++i) { v[i] = 0.5f * i; f[i] = -2.0f * i; }

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    for (int fr = 0; fr < 4; ++fr) {
        for (int i = 0; i < NA * 3; ++i) x[i] = 0.1f * fr + 0.01f * i;
        trr_put_frame(file, NA, fr * 100, 0.2f * fr, box, x, (fr % 2 == 0) ? v : NULL, f);
    }
    // Velocities alone, between the coordinates
    for (int i = 0; i < NA * 3; ++i) v[i] = 7.0f;
    trr_put_frame(file, NA, 350, 0.7f, box, NULL, v, NULL);
    // Half a frame: the header and the box, then nothing.
    trr_put_frame(file, NA, 400, 0.8f, box, NULL, NULL, NULL);
    const int64_t full = (int64_t)md_file_size(file);
    md_file_close(&file);
    {
        // Rewrite the last header to promise coordinates it does not have.
        md_file_t again = {0};
        ASSERT_TRUE(md_file_open(&again, path, MD_FILE_READ | MD_FILE_WRITE));
        const int64_t frame_bytes = 4 + 4 + 4 + 12 + 13 * 4 + 8 + 36;
        uint8_t xsize[4] = { 0, 0, 0, NA * 12 };
        md_file_write_at(again, full - frame_bytes + 4 + 4 + 4 + 12 + 7 * 4, xsize, 4);
        md_file_close(&again);
    }

    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_trr_system_publish_run(&sys, path, TRR_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attributes_t* t = &sys.attributes;

    const md_attribute_t* time = md_attributes_find(t, STR_LIT("run/tryptophan/time"));
    const md_attribute_t* step = md_attributes_find(t, STR_LIT("run/tryptophan/step"));
    ASSERT_TRUE(time && step);
    ASSERT_EQ(4u, time->format.shape[0]);
    EXPECT_NEAR(0.6, ((const double*)time->data)[3], 1.0e-6);
    EXPECT_EQ(300, ((const int64_t*)step->data)[3]);

    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/tryptophan/atom/position")) != NULL);
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/tryptophan/atom/force"))    != NULL);
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/tryptophan/atom/velocity")) == NULL);

    // The velocities: frames 0 and 2, and the one without coordinates
    const md_attribute_t* vtime = md_attributes_find(t, STR_LIT("run/tryptophan/trr/velocity/time"));
    const md_attribute_t* vel   = md_attributes_find(t, STR_LIT("run/tryptophan/trr/velocity/atom/velocity"));
    ASSERT_TRUE(vtime && vel);
    ASSERT_EQ(3u, vtime->format.shape[0]);
    EXPECT_EQ(vtime, md_attributes_axis(t, vel));
    EXPECT_NEAR(0.4, ((const double*)vtime->data)[1], 1.0e-6);
    EXPECT_NEAR(0.7, ((const double*)vtime->data)[2], 1.0e-6);
    float vv[NA * 3];
    md_attribute_slice_t s2 = md_attribute_slice_1(2);
    ASSERT_EQ((size_t)(NA * 3), md_attribute_extract_slice_f32(vv, NA * 3, vel, &s2, md_unit_none()));
    EXPECT_EQ(70.0f, vv[4]);   // 7 nm/ps
    md_attribute_slice_t s1 = md_attribute_slice_1(1);
    ASSERT_EQ((size_t)(NA * 3), md_attribute_extract_slice_f32(vv, NA * 3, vel, &s1, md_unit_none()));
    EXPECT_EQ(0.5f * 4 * 10.0f, vv[4]);

    const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("unitcell"), STR_LIT("atom/force") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, TRR_RUN, paths, ARRAY_SIZE(paths), md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, NA);
    ASSERT_TRUE(md_system_extract_frame(ex, 3, &st));
    for (int i = 0; i < NA; ++i) {
        EXPECT_NEAR((0.3f + 0.01f * (i * 3 + 0)) * 10.0f, st.xyz[i].x, 1.0e-5f);
        EXPECT_NEAR((0.3f + 0.01f * (i * 3 + 2)) * 10.0f, st.xyz[i].z, 1.0e-5f);
    }
    const md_attribute_t* force = md_attributes_find(&st.attributes, STR_LIT("atom/force"));
    ASSERT_TRUE(force != NULL);
    EXPECT_EQ(-2.0f * 7, ((const float*)force->data)[7]);

    // The whole box, tilt included, in Angstrom.
    EXPECT_NEAR(30.0, st.unitcell.x,  1.0e-4);
    EXPECT_NEAR(10.0, st.unitcell.xy, 1.0e-4);
    EXPECT_NEAR(25.0, st.unitcell.y,  1.0e-4);
    EXPECT_NEAR(-5.0, st.unitcell.xz, 1.0e-4);
    EXPECT_NEAR(7.5,  st.unitcell.yz, 1.0e-4);
    EXPECT_NEAR(20.0, st.unitcell.z,  1.0e-4);
    md_system_extract_end(ex);

    // Along the run, the velocities are there at the frames that have them and absent at the
    // others: one state reused through all of them never shows the last frame's in their place.
    const str_t with_velocity[] = { STR_LIT("atom/position"), STR_LIT("trr/velocity/atom/velocity") };
    ex = md_system_extract_begin(&sys, TRR_RUN, with_velocity, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    for (int64_t fr = 0; fr < 4; ++fr) {
        ASSERT_TRUE(md_system_extract_frame(ex, fr, &st));
        const md_attribute_t* v = md_attributes_find(&st.attributes, STR_LIT("trr/velocity/atom/velocity"));
        if (fr % 2 == 0) {
            ASSERT_TRUE(v != NULL);
            EXPECT_EQ(0.5f * 4 * 10.0f, ((const float*)v->data)[4]);
        } else {
            EXPECT_TRUE(v == NULL);
        }
    }
    md_system_extract_end(ex);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

#pragma once

// Checking a published run against values recorded from the trajectory readers before they were
// removed. Each reference is one frame: the coordinate sums over all atoms, the first and last atom,
// and the cell. The tolerances leave room for the last bit a float parse may differ in between
// platforms, and none for a wrong frame, a wrong atom order or a wrong unit.

#include "utest.h"
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <math.h>

typedef struct run_ref_t {
    size_t frame;
    double sum[3];
    double first[3];
    double last[3];
    double cell[6];     // x, y, z, xy, xz, yz
} run_ref_t;

static inline size_t run_num_frames(const md_system_t* sys, str_t run) {
    char buf[512];
    const md_attribute_t* time = md_attributes_find(&sys->attributes, md_run_path(buf, sizeof(buf), run, STR_LIT("time")));
    return time ? time->format.shape[0] : 0;
}

static inline size_t run_num_atoms(const md_system_t* sys, str_t run) {
    char buf[512];
    const md_attribute_t* pos = md_attributes_find(&sys->attributes, md_run_path(buf, sizeof(buf), run, STR_LIT("atom/position")));
    return pos ? pos->format.shape[1] : 0;
}

// One frame of positions and cell into state, which has room for the run's atoms.
static inline bool run_extract_one(md_system_state_t* state, const md_system_t* sys, str_t run, int64_t frame) {
    const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 2, md_get_heap_allocator());
    if (!ex) return false;
    const bool ok = md_system_extract_frame(ex, frame, state);
    md_system_extract_end(ex);
    return ok;
}

static inline bool run_near(double ref, double got, double rel, double abs_tol) {
    return fabs(ref - got) <= rel * fabs(ref) + abs_tol;
}

static inline void run_check_refs(int* utest_result, const md_system_t* sys, str_t run, size_t expected_frames, size_t expected_atoms, const run_ref_t* refs, size_t num_refs) {
    ASSERT_EQ(expected_frames, run_num_frames(sys, run));
    const size_t N = run_num_atoms(sys, run);
    ASSERT_EQ(expected_atoms, N);

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, N);

    const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 2, md_get_heap_allocator());
    EXPECT_TRUE(ex != NULL);
    for (size_t r = 0; ex && r < num_refs; ++r) {
        const run_ref_t* ref = &refs[r];
        if (!md_system_extract_frame(ex, (int64_t)ref->frame, &st)) {
            *utest_result = UTEST_TEST_FAILURE;
            UTEST_PRINTF("run_check: frame %zu failed to extract\n", ref->frame);
            continue;
        }
        double sum[3] = {0};
        for (size_t i = 0; i < N; ++i) {
            sum[0] += st.xyz[i].x;
            sum[1] += st.xyz[i].y;
            sum[2] += st.xyz[i].z;
        }
        const double first[3] = { st.xyz[0].x, st.xyz[0].y, st.xyz[0].z };
        const double last[3]  = { st.xyz[N - 1].x, st.xyz[N - 1].y, st.xyz[N - 1].z };
        const double cell[6]  = { st.unitcell.x, st.unitcell.y, st.unitcell.z, st.unitcell.xy, st.unitcell.xz, st.unitcell.yz };
        for (int k = 0; k < 3; ++k) {
            // A sum of N values each good to a few ulps
            if (!run_near(ref->sum[k], sum[k], 1.0e-6, 1.0e-4 * (double)N)) {
                *utest_result = UTEST_TEST_FAILURE;
                UTEST_PRINTF("run_check: frame %zu sum[%d] %.9g, expected %.9g\n", ref->frame, k, sum[k], ref->sum[k]);
            }
            if (!run_near(ref->first[k], first[k], 1.0e-6, 1.0e-5) || !run_near(ref->last[k], last[k], 1.0e-6, 1.0e-5)) {
                *utest_result = UTEST_TEST_FAILURE;
                UTEST_PRINTF("run_check: frame %zu atom coordinate %d differs\n", ref->frame, k);
            }
        }
        for (int k = 0; k < 6; ++k) {
            if (!run_near(ref->cell[k], cell[k], 1.0e-6, 1.0e-4)) {
                *utest_result = UTEST_TEST_FAILURE;
                UTEST_PRINTF("run_check: frame %zu cell[%d] %.9g, expected %.9g\n", ref->frame, k, cell[k], ref->cell[k]);
            }
        }
        EXPECT_EQ((double)ref->frame, st.frame);
    }
    md_system_extract_end(ex);
    md_vm_arena_destroy(arena);
}

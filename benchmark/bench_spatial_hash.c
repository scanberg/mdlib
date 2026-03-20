#include "ubench.h"

#include <md_gro.h>
#include <md_system.h>
#include <core/md_intrinsics.h>
#include <core/md_spatial_acc.h>
#include <core/md_arena_allocator.h>

#define RADIUS 10.0

void spatial_acc_pair_count_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx;
    (void)j_idx;
    (void)ij_dist2;
    uint32_t* count = (uint32_t*)user_param;
    *count += (uint32_t)num_pairs;
}

// This makes little sense logically
// There are specialized routines for internal pair checks within a radius which should be favored over this.
UBENCH_EX(spatial_acc, query_ext_vs_int_pair) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    md_system_t sys = { .alloc = arena };
    if (!md_gro_system_init_from_file(&sys, STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro"))) {
        fprintf(stderr, "Failed to load system for spatial acc benchmark\n");
        md_vm_arena_destroy(arena);
        return;
    }

    md_coord_stream_t stream = md_coord_stream_create_soa(sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count);
    md_spatial_acc_t acc = { .alloc = arena };
    md_spatial_acc_init(&acc, &stream, RADIUS, &sys.unitcell, 0);

    uint32_t count = 0;
    UBENCH_DO_BENCHMARK() {
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, &stream, RADIUS, spatial_acc_pair_count_callback, &count, 0);
        UBENCH_DO_NOTHING(&count);
    }

    md_vm_arena_destroy(arena);
}

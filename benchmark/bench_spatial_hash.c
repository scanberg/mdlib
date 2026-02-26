#include "ubench.h"

#include <md_gro.h>
#include <md_system.h>
#include <core/md_intrinsics.h>
#include <core/md_spatial_acc.h>
#include <core/md_arena_allocator.h>

#define RADIUS 10.0

void spatial_acc_pair_count_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    uint32_t* count = (uint32_t*)user_param;
    *count += (uint32_t)num_pairs;
}

UBENCH_EX(spatial_hash, query_n2_centered) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    md_system_t sys = { 0 };
    bool result = (md_gro_system_loader()->init_from_file(&sys, STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro"), NULL, arena));

    if (result) {
        md_spatial_acc_t acc = { .alloc = arena };
        md_spatial_acc_init(&acc, sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, 10.0, &sys.unitcell);

        uint32_t count = 0;
        UBENCH_DO_BENCHMARK() {
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, sys.atom.x, sys.atom.y, sys.atom.z, NULL, sys.atom.count, RADIUS, spatial_acc_pair_count_callback, &count);
            UBENCH_DO_NOTHING(&count);
        }
    }

    md_vm_arena_destroy(arena);
}

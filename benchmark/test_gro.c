#include "ubench.h"

#include <md_gro.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

UBENCH_EX(gro, load) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, MEGABYTES(1));

    UBENCH_DO_BENCHMARK() {
        md_gro_data_t gro = {0};
        md_gro_data_parse_file(&gro, STR(MD_BENCHMARK_DATA_DIR"/centered.gro"), alloc);
    }

    md_arena_allocator_destroy(alloc);
}
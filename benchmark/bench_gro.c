#include "ubench.h"

#include <md_gro.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UBENCH_EX(gro, load) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_file_o* file = md_file_open(path, MD_FILE_READ);
    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_gro_data_t gro = {0};
        md_gro_data_parse_file(&gro, path, alloc);
    }

    md_arena_allocator_destroy(alloc);
}
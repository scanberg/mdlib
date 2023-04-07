#include "ubench.h"

#include <md_xyz.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UBENCH_EX(xyz, xmol) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/40-40-2-ddba-dyna.xmol");

    md_file_o* file = md_file_open(path, MD_FILE_READ);
    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_xyz_data_t xyz = {0};
        md_xyz_data_parse_file(&xyz, path, alloc);
    }

    md_arena_allocator_destroy(alloc);
}

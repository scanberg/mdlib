#include "ubench.h"

#include <md_gro.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_linear_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <md_util.h>

UBENCH_EX(gro, load) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    str_t path = STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_file_o* file = md_file_open(path, MD_FILE_READ);
    if (!file) {
        MD_LOG_ERROR("Could not open file '%.*s'", path.len, path.ptr);
        return;
    }
    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_gro_data_t gro = {0};
        md_gro_data_parse_file(&gro, path, alloc);
    }

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(gro, postprocess) {
    size_t capacity = MEGABYTES(16);
    void* buffer = md_alloc(md_get_heap_allocator(), capacity);
    md_allocator_i* alloc = md_linear_allocator_create(buffer, capacity);
    str_t path = STR_LIT(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_system_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, 0, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    size_t reset_pos = md_linear_allocator_get_pos(alloc);
    UBENCH_DO_BENCHMARK() {
        md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_BOND_BIT);
        md_linear_allocator_set_pos_back(alloc, reset_pos);
    }

    md_free(md_get_heap_allocator(), buffer, capacity);
}
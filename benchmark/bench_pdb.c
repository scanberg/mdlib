#include "ubench.h"

#include <md_pdb.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UBENCH_EX(pdb, dppc64) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/dppc64.pdb");

    str_t text = load_textfile(path, md_heap_allocator);

    UBENCH_SET_BYTES(text.len);

    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_pdb_data_t pdb = {0};
        //md_pdb_data_parse_file(&pdb, STR(MD_BENCHMARK_DATA_DIR"/dppc64.pdb"), alloc);
        md_pdb_data_parse_str(&pdb, text, alloc);
    }

    md_arena_allocator_destroy(alloc);
    str_free(text, md_heap_allocator);
}

UBENCH_EX(pdb, _100frames) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/100frames.pdb");
    
    md_file_o* file = md_file_open(path, MD_FILE_READ);
    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);
    
    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_pdb_data_t pdb = {0};
        md_pdb_data_parse_file(&pdb, path, alloc);
    }

    md_arena_allocator_destroy(alloc);
}

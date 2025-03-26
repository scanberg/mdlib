#include "ubench.h"

#include <md_pdb.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

UBENCH_EX(pdb, dppc64) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    str_t path = STR_LIT(MD_BENCHMARK_DATA_DIR "/dppc64.pdb");

    str_t text = load_textfile(path, md_get_heap_allocator());
    if (str_empty(text)) {
    	MD_LOG_ERROR("Failed to load file: %.*s\n", STR_ARG(path));
    	return;
    }  
    UBENCH_SET_BYTES(text.len);

    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_pdb_data_t pdb = {0};
        //md_pdb_data_parse_file(&pdb, STR_LIT(MD_BENCHMARK_DATA_DIR"/dppc64.pdb"), alloc);
        md_pdb_data_parse_str(&pdb, text, alloc);
    }

    md_arena_allocator_destroy(alloc);
    str_free(text, md_get_heap_allocator());
}

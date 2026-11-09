#include "ubench.h"

#include <md_pdb.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

UBENCH_EX(pdb, dppc64) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    str_t path = STR_INIT(MD_BENCHMARK_DATA_DIR "/dppc64.pdb");

    str_t text = load_textfile(path, md_get_heap_allocator());
    if (str_empty(text)) {
    	MD_LOG_ERROR("Failed to load file: %.*s\n", STR_ARG(path));
    	return;
    }  
    UBENCH_SET_BYTES(text.len);

    // The system of the file: what loading it is
    UBENCH_DO_BENCHMARK() {
        md_arena_allocator_reset(alloc);
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        md_pdb_system_init_from_str(&sys, &st, text, MD_PDB_OPTION_NONE, NULL);
    }

    md_arena_allocator_destroy(alloc);
    str_free(text, md_get_heap_allocator());
}

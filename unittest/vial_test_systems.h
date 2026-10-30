#pragma once

// The reference's two structures, T and P, loaded once per test file: the VIAL tests compile against them.
//
//   T  test_data/tubulin-A-B.pdb    tubulin with vinblastine, chains A, B and C
//   P  test_data/1ALA-560ns.pdb     an alanine peptide with no chain ids

#include <core/md_arena_allocator.h>
#include <md_system.h>
#include <md_pdb.h>
#include <md_util.h>

#ifndef MD_UNITTEST_DATA_DIR
#error "MD_UNITTEST_DATA_DIR names the test data"
#endif

// 'T' or 'P'; NULL for another name, or when the file does not load
static inline const md_system_t* vial_test_system(char which) {
    static md_system_t sys[2];
    static bool loaded[2];
    if (which != 'T' && which != 'P') return NULL;
    const int i = which == 'T' ? 0 : 1;
    if (!loaded[i]) {
        loaded[i] = true;
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
        sys[i] = (md_system_t){.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        const char* path = i == 0 ? MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb" : MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb";
        if (md_pdb_system_init_from_file(&sys[i], &state, str_from_cstr(path), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE)) {
            md_util_system_infer(&sys[i], &state, MD_UTIL_INFER_ALL);
        }
    }
    return sys[i].atom.count ? &sys[i] : NULL;
}

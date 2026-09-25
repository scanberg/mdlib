#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include "run_check.h"

UTEST(pdb, parse_ordinary) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/1k4r.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 9084);
    EXPECT_EQ(pdb_data.num_connections, 36);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 24);
    EXPECT_EQ(pdb_data.num_sheets, 102);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, tryptophan) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/tryptophan.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 28);

    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_TRUE(md_pdb_system_init_from_data(&sys, &sys_state, &pdb_data, MD_PDB_OPTION_NONE));

    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, unmatched_model_entry) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/dppc64.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 14738);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, parse_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 38);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 5814);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (int64_t i = 0; i < pdb_data.num_models; ++i) {
        char data[6] = {0};
        md_file_seek(file, pdb_data.models[i].byte_offset, MD_FILE_BEG);
        md_file_read(file, data, 5);
        EXPECT_EQ(strncmp(data, "MODEL", 5), 0);
    }
    md_file_close(&file);

    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, create_system) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, arena));
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = { .alloc = arena };
    EXPECT_TRUE(md_pdb_system_init_from_data(&sys, &sys_state, &pdb_data, MD_PDB_OPTION_NONE));
    ASSERT_EQ(sys.atom.count, pdb_data.num_atom_coordinates);

    EXPECT_EQ(1185, sys.component.count);
    EXPECT_EQ(3, sys.instance.count);
    EXPECT_EQ(1, sys.entity.count);

    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, pdb_data.atom_coordinates[i].x);
        EXPECT_EQ(sys_state.xyz[i].y, pdb_data.atom_coordinates[i].y);
        EXPECT_EQ(sys_state.xyz[i].z, pdb_data.atom_coordinates[i].z);
    }

    md_system_free(&sys);
    
    md_pdb_data_free(&pdb_data, arena);
    md_arena_allocator_destroy(arena);
}

UTEST(pdb, parse_nonexistent_file) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/nonexistent.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_FALSE(result);
    
    // Should be safe to free even when parsing failed
    md_pdb_data_free(&pdb_data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, parse_empty_path) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = {0};  // Empty string
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, md_get_heap_allocator());
    EXPECT_FALSE(result);
    
    md_pdb_data_free(&pdb_data, md_get_heap_allocator());
    md_arena_allocator_destroy(alloc);
}
// ### RUN ###

#define PDB_RUN STR_LIT("run/1ala")

// Recorded from the trajectory reader this replaced
static const run_ref_t pdb_refs[] = {
    { 0,  {3568.34002, 7394.96399, 3699.66098}, {23.7040005, 23, 21.3549995},         {25.0020008, 75.5550003, 23.5319996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 19, {3297.43299, 6833.70401, 3418.89901}, {22.4950008, 43.8180008, 11.6870003}, {7.09700012, 53.9360008, 23.4740009}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
    { 37, {3298.315, 6835.46697, 3419.76601},   {28.4710007, 57.6910019, 25.2150002}, {8.61200047, 35.6669998, 14.9329996}, {46.6450005, 96.6660004, 48.3619995, 0, 0, 0} },
};

// One frame per model, time as model ordinals, and the file's one CRYST1 cell at every frame.
UTEST(pdb, run_matches_reference) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_system_t sys = {.alloc = arena};
    md_system_state_t sys_state = {.alloc = arena};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &sys_state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, path, PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, PDB_RUN, 38, 153, pdb_refs, ARRAY_SIZE(pdb_refs));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/1ala/time"));
    const md_attribute_t* pos  = md_attributes_find(&sys.attributes, STR_LIT("run/1ala/atom/position"));
    ASSERT_TRUE(time && pos);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(9.0, ((const double*)time->data)[9]);

    // The first model is the structure's own coordinates.
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, sys.atom.count);
    ASSERT_TRUE(run_extract_one(&got, &sys, PDB_RUN, 0));
    EXPECT_EQ(0, MEMCMP(sys_state.xyz, got.xyz, sys.atom.count * sizeof(vec3_t)));

    // One atom of one model: the text is read up to that atom and no further.
    float xyz[3];
    ASSERT_TRUE(run_extract_one(&got, &sys, PDB_RUN, 9));
    md_attribute_slice_t one = md_attribute_slice_2(9, 100);
    ASSERT_EQ(3u, md_attribute_extract_slice_f32(xyz, 3, pos, &one, md_unit_none()));
    EXPECT_EQ(got.xyz[100].x, xyz[0]);
    EXPECT_EQ(got.xyz[100].z, xyz[2]);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A file of one model is a structure, not a run.
UTEST(pdb, run_needs_several_models) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    EXPECT_FALSE(md_pdb_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), PDB_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/1ala/time")) == NULL);
    md_vm_arena_destroy(arena);
}

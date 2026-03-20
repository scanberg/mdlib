#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_trajectory.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UTEST(pdb, parse_ordinary) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1k4r.pdb");
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

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/tryptophan.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 28);

    md_system_t sys = {.alloc = alloc};
    EXPECT_TRUE(md_pdb_system_init(&sys, &pdb_data, MD_PDB_OPTION_NONE));

    md_arena_allocator_destroy(alloc);
}

UTEST(pdb, unmatched_model_entry) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/dppc64.pdb");
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
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
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

UTEST(pdb, trajectory_i) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));

    ASSERT_TRUE(sys.trajectory);

    EXPECT_EQ(md_trajectory_num_atoms(sys.trajectory), 153);
    EXPECT_EQ(md_trajectory_num_frames(sys.trajectory), 38);

    const size_t mem_size = md_trajectory_num_atoms(sys.trajectory) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(arena, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(sys.trajectory) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(sys.trajectory) * 2;
    for (int64_t i = 0; i < md_trajectory_num_frames(sys.trajectory); ++i) {
        md_trajectory_frame_header_t header = {0};
        EXPECT_TRUE(md_trajectory_load_frame(sys.trajectory, i, &header, x, y, z));
    }

    md_arena_allocator_destroy(arena);
}

UTEST(pdb, trajectory_reader_i) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_trajectory_i* traj = md_pdb_trajectory_create(path, arena, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(arena, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_reader_i reader = {0};
    ASSERT_TRUE(md_trajectory_reader_init(&reader, traj));

    md_trajectory_frame_header_t header = {0};
    EXPECT_TRUE(md_trajectory_reader_load_frame(reader, 0, &header, x, y, z));
    EXPECT_EQ(153, header.num_atoms);
    EXPECT_TRUE(md_trajectory_reader_load_frame(reader, md_trajectory_num_frames(traj) - 1, &header, x, y, z));
    EXPECT_EQ(153, header.num_atoms);

    md_trajectory_reader_free(&reader);
    md_pdb_trajectory_free(traj);
    md_arena_allocator_destroy(arena);
}

UTEST(pdb, create_system) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, arena));
    md_system_t sys = {.alloc = arena};
    EXPECT_TRUE(md_pdb_system_init_from_data(&sys, &pdb_data, MD_PDB_OPTION_NONE));
    ASSERT_EQ(sys.atom.count, pdb_data.num_atom_coordinates);

    EXPECT_EQ(1185, sys.component.count);
    EXPECT_EQ(3, sys.instance.count);
    EXPECT_EQ(1, sys.entity.count);

    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys.atom.x[i], pdb_data.atom_coordinates[i].x);
        EXPECT_EQ(sys.atom.y[i], pdb_data.atom_coordinates[i].y);
        EXPECT_EQ(sys.atom.z[i], pdb_data.atom_coordinates[i].z);
    }

    md_system_free(&sys);
    
    md_pdb_data_free(&pdb_data, arena);
    md_arena_allocator_destroy(arena);
}

UTEST(pdb, parse_nonexistent_file) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.pdb");
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
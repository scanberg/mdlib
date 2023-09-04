#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

UTEST(pdb, parse_ordinary) {
    str_t path = STR(MD_UNITTEST_DATA_DIR"/1k4r.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, md_heap_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 9084);
    EXPECT_EQ(pdb_data.num_connections, 36);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 24);
    EXPECT_EQ(pdb_data.num_sheets, 102);

    md_pdb_data_free(&pdb_data, md_heap_allocator);
}

UTEST(pdb, unmatched_model_entry) {
    str_t path = STR(MD_UNITTEST_DATA_DIR"/dppc64.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, md_heap_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 0);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 14738);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_pdb_data_free(&pdb_data, md_heap_allocator);
}

UTEST(pdb, parse_trajectory) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_pdb_data_t pdb_data = {0};
    bool result = md_pdb_data_parse_file(&pdb_data, path, md_heap_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(pdb_data.num_models, 38);
    EXPECT_EQ(pdb_data.num_atom_coordinates, 5814);
    EXPECT_EQ(pdb_data.num_cryst1, 1);
    EXPECT_EQ(pdb_data.num_connections, 0);
    EXPECT_EQ(pdb_data.num_helices, 0);
    EXPECT_EQ(pdb_data.num_sheets, 0);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (int64_t i = 0; i < pdb_data.num_models; ++i) {
        char data[6] = {0};
        md_file_seek(file, pdb_data.models[i].byte_offset, MD_FILE_BEG);
        md_file_read(file, data, 5);
        EXPECT_EQ(strncmp(data, "MODEL", 5), 0);
    }
    md_file_close(file);

    md_pdb_data_free(&pdb_data, md_heap_allocator);
}

UTEST(pdb, trajectory_i) {
    const str_t path = STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    md_trajectory_i* traj = md_pdb_trajectory_create(path, md_heap_allocator);
    ASSERT_TRUE(traj);

    EXPECT_EQ(md_trajectory_num_atoms(traj), 153);
    EXPECT_EQ(md_trajectory_num_frames(traj), 38);

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_temp_allocator, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
    }

    md_free(md_temp_allocator, mem_ptr, mem_size);
    md_pdb_trajectory_free(traj);
}

UTEST(pdb, create_molecule) {
    md_allocator_i* alloc = md_heap_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR "/1k4r.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, alloc));

    md_molecule_t mol = {0};
    EXPECT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));
    ASSERT_EQ(mol.atom.count, pdb_data.num_atom_coordinates);

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], pdb_data.atom_coordinates[i].x);
        EXPECT_EQ(mol.atom.y[i], pdb_data.atom_coordinates[i].y);
        EXPECT_EQ(mol.atom.z[i], pdb_data.atom_coordinates[i].z);
    }

    md_molecule_free(&mol, alloc);
    
    /*
    
    EXPECT_TRUE(md_pdb_molecule_loader()->init_from_file(&mol, path, alloc));

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], pdb_data.atom_coordinates[i].x);
        EXPECT_EQ(mol.atom.y[i], pdb_data.atom_coordinates[i].y);
        EXPECT_EQ(mol.atom.z[i], pdb_data.atom_coordinates[i].z);
    }

    EXPECT_TRUE(md_pdb_molecule_loader()->free(&mol, alloc));
    */
    md_pdb_data_free(&pdb_data, alloc);
}
#include "utest.h"
#include <string.h>

#include <md_xyz.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_file.h>

UTEST(xyz, xyz_standard) {
    str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 10);
    EXPECT_EQ(data.num_coordinates, 2280 * 10);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "2280", 4), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, xyz_xmol) {
    str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 50);
    EXPECT_EQ(data.num_coordinates, 50 * 540);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "   540", 6), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, xyz_tinker) {
    str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/full.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 1);
    EXPECT_EQ(data.num_coordinates, 404);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "         404  molden generated tinker .xyz (mm3 param.)", 54), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, xyz_tinker_arc) {
    str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/full.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 10);
    EXPECT_EQ(data.num_coordinates, 10 * 404);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[64];
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "   404  molden generated tinker .xyz (mm3 param.)", 48), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, create_molecule) {
    md_allocator_i* alloc = default_allocator;
    str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");

    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_file(&data, path, alloc));

    md_molecule_t mol = {0};
    EXPECT_TRUE(md_xyz_molecule_init(&mol, &data, alloc));
    ASSERT_GT(data.num_models, 0);
    ASSERT_EQ(mol.atom.count, data.models[0].end_coord_index - data.models[0].beg_coord_index);

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.coordinates[i].x);
        EXPECT_EQ(mol.atom.y[i], data.coordinates[i].y);
        EXPECT_EQ(mol.atom.z[i], data.coordinates[i].z);
    }

    md_molecule_free(&mol, alloc);

    md_xyz_data_free(&data, alloc);
}

UTEST(xyz, trajectory_i) {
    const str_t path = MAKE_STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
    md_trajectory_i* traj = md_xyz_trajectory_create(path, default_allocator);
    ASSERT_TRUE(traj);

    EXPECT_EQ(md_trajectory_num_atoms(traj), 2280);
    EXPECT_EQ(md_trajectory_num_frames(traj), 10);

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(default_temp_allocator, mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(header.num_atoms, 2280);
    }

    md_free(default_temp_allocator, mem_ptr, mem_size);
    md_xyz_trajectory_free(traj);
}

#include "utest.h"
#include <string.h>

#include <md_xyz.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_array.h>

UTEST(xyz, xyz_standard) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 10);
    EXPECT_EQ(data.num_coordinates, 2280 * 10);

    EXPECT_NEAR(data.coordinates[0].x, 57.834f, 1.0e-5f);
    EXPECT_NEAR(data.coordinates[0].y, 36.568f, 1.0e-5f);
    EXPECT_NEAR(data.coordinates[0].z, 62.491f, 1.0e-5f);

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
    str_t path = STR(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 50);
    EXPECT_EQ(data.num_coordinates, 50 * 540);

    EXPECT_STREQ(data.coordinates[0].element_symbol, "Au");
    EXPECT_NEAR( data.coordinates[0].x, -2.264467f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].y,  1.246472f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].z,  3.629187f, 1.0e-5f);

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
    str_t path = STR(MD_UNITTEST_DATA_DIR "/full.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 1);
    EXPECT_EQ(data.num_coordinates, 404);

    EXPECT_EQ(   data.coordinates[0].atom_index, 1);
    EXPECT_STREQ(data.coordinates[0].element_symbol, "C");
    EXPECT_NEAR( data.coordinates[0].x, 18.673994f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].y,  1.292906f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].z,  0.733642f, 1.0e-5f);
    EXPECT_EQ(   data.coordinates[0].atom_type, 2);
    EXPECT_EQ(   data.coordinates[0].connectivity[0], 2);
    EXPECT_EQ(   data.coordinates[0].connectivity[1], 20);
    EXPECT_EQ(   data.coordinates[0].connectivity[2], 270);
    
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
    str_t path = STR(MD_UNITTEST_DATA_DIR "/full.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 10);
    EXPECT_EQ(data.num_coordinates, 10 * 404);

    EXPECT_EQ(   data.coordinates[0].atom_index, 1);
    EXPECT_STREQ(data.coordinates[0].element_symbol, "C");
    EXPECT_NEAR( data.coordinates[0].x, 17.935708f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].y,  1.368677f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].z,  0.600876f, 1.0e-5f);
    EXPECT_EQ(   data.coordinates[0].atom_type, 2);
    EXPECT_EQ(   data.coordinates[0].connectivity[0], 2);
    EXPECT_EQ(   data.coordinates[0].connectivity[1], 20);
    EXPECT_EQ(   data.coordinates[0].connectivity[2], 270);

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

UTEST(xyz, o2_arc) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/o2.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 2000);
    EXPECT_EQ(data.num_coordinates, 2000 * 2);

    EXPECT_STREQ(data.coordinates[0].element_symbol, "O");
    EXPECT_NEAR( data.coordinates[0].x, -1.977261f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].y, -1.149547f, 1.0e-5f);
    EXPECT_NEAR( data.coordinates[0].z, -0.000606f, 1.0e-5f);
    EXPECT_EQ(   data.coordinates[0].atom_type, 7);
    EXPECT_EQ(   data.coordinates[0].connectivity[0], 2);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[64];
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "     2", 6), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, h2o_arc) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/h2o.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 2000);
    EXPECT_EQ(data.num_coordinates, 2000 * 3);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[64];
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "     3  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, ch4_arc) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/ch4.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, default_allocator);
    EXPECT_TRUE(result);
    EXPECT_EQ(data.num_models, 2000);
    EXPECT_EQ(data.num_coordinates, 2000 * 5);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_NE((FILE*)file, NULL);
    for (int64_t i = 0; i < data.num_models; ++i) {
        char str[64];
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        fgets(str, sizeof(str), (FILE*)file);
        EXPECT_EQ(strncmp(str, "     5  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, default_allocator);
}

UTEST(xyz, create_molecule) {
    md_allocator_i* alloc = default_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");

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
    const str_t path = STR(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
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

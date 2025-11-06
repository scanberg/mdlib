#include "utest.h"
#include <string.h>

#include <md_xyz.h>
#include <md_trajectory.h>
#include <md_system.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_array.h>

UTEST(xyz, xyz_standard) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(10, data.num_models);
    EXPECT_EQ(2280 * 10, data.num_coordinates);

    EXPECT_NEAR(57.834f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(36.568f, data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(62.491f, data.coordinates[0].z, 1.0e-5f);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "2280", 4), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, c60) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/C60-Ih.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(1, data.num_models);
    EXPECT_EQ(60, data.num_coordinates);

    EXPECT_EQ('C', data.coordinates->element_symbol[0]);
    EXPECT_NEAR(2.16650f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(0.59060f, data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(2.58740f, data.coordinates[0].z, 1.0e-5f);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, xyz_xmol) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(50, data.num_models);
    EXPECT_EQ(50 * 540, data.num_coordinates);

    EXPECT_STREQ("Au", data.coordinates[0].element_symbol);
    EXPECT_NEAR(-2.264467f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR( 1.246472f, data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR( 3.629187f, data.coordinates[0].z, 1.0e-5f);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "   540", 6), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, xyz_tinker) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/full.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(1, data.num_models);
    EXPECT_EQ(404, data.num_coordinates);

    EXPECT_EQ(1,            data.coordinates[0].atom_index);
    EXPECT_STREQ("C",       data.coordinates[0].element_symbol);
    EXPECT_NEAR(18.673994f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(1.292906f,  data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(0.733642f,  data.coordinates[0].z, 1.0e-5f);
    EXPECT_EQ(2,            data.coordinates[0].atom_type);
    EXPECT_EQ(2,            data.coordinates[0].connectivity[0]);
    EXPECT_EQ(20,           data.coordinates[0].connectivity[1]);
    EXPECT_EQ(270,          data.coordinates[0].connectivity[2]);
    
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "         404  molden generated tinker .xyz (mm3 param.)", 54), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, xyz_tinker_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/full.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(10, data.num_models);
    EXPECT_EQ(10 * 404, data.num_coordinates);

    EXPECT_EQ(1,            data.coordinates[0].atom_index);
    EXPECT_STREQ("C",       data.coordinates[0].element_symbol);
    EXPECT_NEAR(17.935708f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(1.368677f,  data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(0.600876f,  data.coordinates[0].z, 1.0e-5f);
    EXPECT_EQ(2,            data.coordinates[0].atom_type);
    EXPECT_EQ(2,            data.coordinates[0].connectivity[0]);
    EXPECT_EQ(20,           data.coordinates[0].connectivity[1]);
    EXPECT_EQ(270,          data.coordinates[0].connectivity[2]);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "   404  molden generated tinker .xyz (mm3 param.)", 48), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, o2_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/o2.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 2, data.num_coordinates);

    EXPECT_STREQ("O", data.coordinates[0].element_symbol);
    EXPECT_NEAR(-1.977261f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(-1.149547f, data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(-0.000606f, data.coordinates[0].z, 1.0e-5f);
    EXPECT_EQ(7,            data.coordinates[0].atom_type);
    EXPECT_EQ(2,            data.coordinates[0].connectivity[0]);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     2", 6), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, h2o_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/h2o.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 3, data.num_coordinates);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     3  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, ch4_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/ch4.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 5, data.num_coordinates);

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     5  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, extended_xyz) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/extended.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());

    ASSERT_TRUE(result);
    EXPECT_EQ(1,   data.num_models);
    EXPECT_EQ(456, data.num_coordinates);

    if (data.num_coordinates > 0) {
        EXPECT_STREQ("Zr", data.coordinates[0].element_symbol);
        EXPECT_NEAR(18.53562587f, data.coordinates[0].x, 1.0e-5f);
        EXPECT_NEAR(10.57149039f, data.coordinates[0].y, 1.0e-5f);
        EXPECT_NEAR(10.42623774f, data.coordinates[0].z, 1.0e-5f);
    }

    // Lattice="
    // 20.94815017098275 -3.412045517350664e-05 -2.2269710615728675e-05
    // -3.431827648979917e-05 20.947967304256764 -1.0983820669246559e-05
    // -2.245407280127725e-05 -1.1267829933312672e-05 20.94797217314631" 
    if (data.num_models > 0) {
        EXPECT_NEAR(data.models[0].cell[0][0], 20.94815017098275,       1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][1], -3.412045517350664e-05,  1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][2], -2.2269710615728675e-05, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[1][0], -3.431827648979917e-05,  1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][1], 20.947967304256764,      1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][2], -1.0983820669246559e-05, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[2][0], -2.245407280127725e-05,  1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][1], -1.1267829933312672e-05, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][2], 20.94797217314631,       1.0e-5f);
    }

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, extended1_xyz) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/extended1.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());

    ASSERT_TRUE(result);
    EXPECT_EQ(1,   data.num_models);
    EXPECT_EQ(192, data.num_coordinates);

    if (data.num_coordinates > 0) {
        EXPECT_STREQ("Si", data.coordinates[0].element_symbol);
        EXPECT_NEAR(0, data.coordinates[0].x, 1.0e-5f);
        EXPECT_NEAR(0, data.coordinates[0].y, 1.0e-5f);
        EXPECT_NEAR(0, data.coordinates[0].z, 1.0e-5f);
    }

    // Lattice="
    // 20.94815017098275 -3.412045517350664e-05 -2.2269710615728675e-05
    // -3.431827648979917e-05 20.947967304256764 -1.0983820669246559e-05
    // -2.245407280127725e-05 -1.1267829933312672e-05 20.94797217314631" 
    if (data.num_models > 0) {
        EXPECT_NEAR(data.models[0].cell[0][0], 14.24,   1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][1], 0,       1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][2], 0,       1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[1][0], 0,       1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][1], 14.24,   1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][2], 0,       1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[2][0], 0,       1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][1], 0,       1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][2], 14.24,   1.0e-5f);
    }

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, extended2_xyz) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/extended2.xyz");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());

    ASSERT_TRUE(result);
    EXPECT_EQ(1,  data.num_models);
    EXPECT_EQ(22, data.num_coordinates);

    if (data.num_coordinates > 0) {
        // 12.01320886      16.80849934      -2.98503290
        EXPECT_STREQ("H", data.coordinates[0].element_symbol);
        EXPECT_NEAR(12.01320886, data.coordinates[0].x, 1.0e-5f);
        EXPECT_NEAR(16.80849934, data.coordinates[0].y, 1.0e-5f);
        EXPECT_NEAR(-2.98503290, data.coordinates[0].z, 1.0e-5f);
    }

    if (data.num_models > 0) {
        EXPECT_NEAR(data.models[0].cell[0][0], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][1], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][2], 0, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[1][0], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][1], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][2], 0, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[2][0], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][1], 0, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][2], 0, 1.0e-5f);
    }

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, extended_xyz_lattice_braced) {
    str_t input = STR_LIT(
        "4\n"
        "Lattice={1 2 3 4 5 6 7 8 9} "
        "Properties=species:S:1:pos:R:3:forces:R:3:energies:R:1\n"
        "Zr      18.53562587      10.57149039      10.42623774       0.00024027       0.00000174       0.00000768      -2.78259732\n"
        "Zr      18.53563066       0.09750089      20.90023074       0.00021729      -0.00005997      -0.00005202      -2.78258954\n"
        "Zr       8.06153933      10.57150349      20.90023721       0.00023774       0.00003751      -0.00003056      -2.78259378\n"
        "Zr       8.06156708       0.09752582      10.42625768       0.00025192      -0.00002267      -0.00007721      -2.78258456\n"
	);

    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_str(&data, input, md_get_heap_allocator());

    ASSERT_TRUE(result);
    EXPECT_EQ(1, data.num_models);
    EXPECT_EQ(4, data.num_coordinates);

    if (data.num_coordinates > 0) {
        EXPECT_STREQ("Zr", data.coordinates[0].element_symbol);
        EXPECT_NEAR(18.53562587f, data.coordinates[0].x, 1.0e-5f);
        EXPECT_NEAR(10.57149039f, data.coordinates[0].y, 1.0e-5f);
        EXPECT_NEAR(10.42623774f, data.coordinates[0].z, 1.0e-5f);
    }

    if (data.num_models > 0) {
        EXPECT_NEAR(data.models[0].cell[0][0], 1, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][1], 2, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][2], 3, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[1][0], 4, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][1], 5, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][2], 6, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[2][0], 7, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][1], 8, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][2], 9, 1.0e-5f);
    }

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, extended_xyz_lattice_array) {
    str_t input = STR_LIT(
        "4\n"
        "Lattice=[[1,2, 3], [4,5,6], [ 7 , 8, 9]] "
        "Properties=species:S:1:pos:R:3:forces:R:3:energies:R:1\n"
        "Zr      18.53562587      10.57149039      10.42623774       0.00024027       0.00000174       0.00000768      -2.78259732\n"
        "Zr      18.53563066       0.09750089      20.90023074       0.00021729      -0.00005997      -0.00005202      -2.78258954\n"
        "Zr       8.06153933      10.57150349      20.90023721       0.00023774       0.00003751      -0.00003056      -2.78259378\n"
        "Zr       8.06156708       0.09752582      10.42625768       0.00025192      -0.00002267      -0.00007721      -2.78258456\n"
    );

    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_str(&data, input, md_get_heap_allocator());

    ASSERT_TRUE(result);
    EXPECT_EQ(1, data.num_models);
    EXPECT_EQ(4, data.num_coordinates);

    if (data.num_coordinates > 0) {
        EXPECT_STREQ("Zr", data.coordinates[0].element_symbol);
        EXPECT_NEAR(18.53562587f, data.coordinates[0].x, 1.0e-5f);
        EXPECT_NEAR(10.57149039f, data.coordinates[0].y, 1.0e-5f);
        EXPECT_NEAR(10.42623774f, data.coordinates[0].z, 1.0e-5f);
    }

    if (data.num_models > 0) {
        EXPECT_NEAR(data.models[0].cell[0][0], 1, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][1], 2, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[0][2], 3, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[1][0], 4, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][1], 5, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[1][2], 6, 1.0e-5f);

        EXPECT_NEAR(data.models[0].cell[2][0], 7, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][1], 8, 1.0e-5f);
        EXPECT_NEAR(data.models[0].cell[2][2], 9, 1.0e-5f);
    }

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, create_molecule) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");

    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_file(&data, path, alloc));

    md_system_t mol = {0};
    EXPECT_TRUE(md_xyz_molecule_init(&mol, &data, alloc));
    ASSERT_GT(data.num_models, 0);
    ASSERT_EQ(mol.atom.count, data.models[0].end_coord_index - data.models[0].beg_coord_index);

    for (size_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.coordinates[i].x);
        EXPECT_EQ(mol.atom.y[i], data.coordinates[i].y);
        EXPECT_EQ(mol.atom.z[i], data.coordinates[i].z);
    }

    md_system_free(&mol, alloc);

    md_xyz_data_free(&data, alloc);
}

UTEST(xyz, trajectory_i) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
    md_trajectory_i* traj = md_xyz_trajectory_create(path, md_get_heap_allocator(), MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    EXPECT_EQ(2280, md_trajectory_num_atoms(traj));
    EXPECT_EQ(10, md_trajectory_num_frames(traj));

    const size_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_temp_allocator(), mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float *z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (size_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(2280, header.num_atoms);
    }

    md_free(md_get_temp_allocator(), mem_ptr, mem_size);
    md_xyz_trajectory_free(traj);
}


UTEST(xyz, comprehensive_c720) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/c720.xyz");
    
    md_system_t mol = {0};
    bool result = md_xyz_system_loader()->init_from_file(&mol, path, NULL, alloc);
    ASSERT_TRUE(result);
    
    // C720 should have exactly 720 carbon atoms
    EXPECT_EQ(mol.atom.count, 720);
    
    // All atoms should be carbon
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(md_atom_atomic_number(&mol.atom, i), 6); // Carbon atomic number
    }
    
    // Check that coordinates are not all the same (should be a 3D structure)
    bool has_variation_x = false, has_variation_y = false, has_variation_z = false;
    float first_x = mol.atom.x[0], first_y = mol.atom.y[0], first_z = mol.atom.z[0];
    for (int64_t i = 1; i < mol.atom.count; ++i) {
        if (fabsf(mol.atom.x[i] - first_x) > 0.01f) has_variation_x = true;
        if (fabsf(mol.atom.y[i] - first_y) > 0.01f) has_variation_y = true;
        if (fabsf(mol.atom.z[i] - first_z) > 0.01f) has_variation_z = true;
    }
    EXPECT_TRUE(has_variation_x);
    EXPECT_TRUE(has_variation_y);
    EXPECT_TRUE(has_variation_z);
    
    md_system_free(&mol, alloc);
}

UTEST(xyz, error_handling) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test nonexistent file
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.xyz");
    md_system_t mol = {0};
    bool result = md_xyz_system_loader()->init_from_file(&mol, path, NULL, alloc);
    EXPECT_FALSE(result);
    md_system_free(&mol, alloc);
    
    // Test empty path
    str_t empty_path = {0,0};
    result = md_xyz_system_loader()->init_from_file(&mol, empty_path, NULL, alloc);
    EXPECT_FALSE(result);
    md_system_free(&mol, alloc);
}

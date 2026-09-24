#include "utest.h"
#include <string.h>
#include <stdio.h>

#include "run_check.h"

#include <md_xyz.h>
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

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "2280", 4), 0);
    }
    md_file_close(&file);

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

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[8] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "   540", 6), 0);
    }
    md_file_close(&file);

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
    
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "         404  molden generated tinker .xyz (mm3 param.)", 54), 0);
    }
    md_file_close(&file);

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

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "   404  molden generated tinker .xyz (mm3 param.)", 48), 0);
    }
    md_file_close(&file);

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

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     2", 6), 0);
    }
    md_file_close(&file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, h2o_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/h2o.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 3, data.num_coordinates);

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     3  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(&file);

    md_xyz_data_free(&data, md_get_heap_allocator());
}

UTEST(xyz, ch4_arc) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/ch4.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 5, data.num_coordinates);

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    for (size_t i = 0; i < data.num_models; ++i) {
        char str[64] = {0};
        md_file_seek(file, data.models[i].byte_offset, MD_FILE_BEG);
        EXPECT_NE(md_file_read(file, str, sizeof(str)), 0);
        EXPECT_EQ(strncmp(str, "     5  molden generated tinker .xyz (mm3 param.)", 49), 0);
    }
    md_file_close(&file);

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

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_TRUE(md_xyz_system_init_from_data(&sys, &sys_state, &data, MD_XYZ_OPTION_DISABLE_CACHE_WRITE));
    ASSERT_GT(data.num_models, 0);
    ASSERT_EQ(sys.atom.count, data.models[0].end_coord_index - data.models[0].beg_coord_index);

    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, data.coordinates[i].x);
        EXPECT_EQ(sys_state.xyz[i].y, data.coordinates[i].y);
        EXPECT_EQ(sys_state.xyz[i].z, data.coordinates[i].z);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
    md_xyz_data_free(&data, alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST(xyz, comprehensive_c720) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/c720.xyz");
    
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    bool result = md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_NONE);
    ASSERT_TRUE(result);
    
    // C720 should have exactly 720 carbon atoms
    EXPECT_EQ(sys.atom.count, 720);
    
    // All atoms should be carbon
    for (int64_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(md_atom_atomic_number(&sys.atom, i), 6); // Carbon atomic number
    }
    
    // Check that coordinates are not all the same (should be a 3D structure)
    bool has_variation_x = false, has_variation_y = false, has_variation_z = false;
    float first_x = sys_state.xyz[0].x, first_y = sys_state.xyz[0].y, first_z = sys_state.xyz[0].z;
    for (int64_t i = 1; i < sys.atom.count; ++i) {
        if (fabsf(sys_state.xyz[i].x - first_x) > 0.01f) has_variation_x = true;
        if (fabsf(sys_state.xyz[i].y - first_y) > 0.01f) has_variation_y = true;
        if (fabsf(sys_state.xyz[i].z - first_z) > 0.01f) has_variation_z = true;
    }
    EXPECT_TRUE(has_variation_x);
    EXPECT_TRUE(has_variation_y);
    EXPECT_TRUE(has_variation_z);
    
    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

UTEST(xyz, error_handling) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test nonexistent file
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/nonexistent.xyz");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    bool result = md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_NONE);
    EXPECT_FALSE(result);

    md_system_reset(&sys);
    
    // Test empty path
    str_t empty_path = {0};
    result = md_xyz_system_init_from_file(&sys, &sys_state, empty_path, MD_XYZ_OPTION_NONE);
    EXPECT_FALSE(result);
    md_system_free(&sys);
}

// ### RUN ###

#define XYZ_RUN STR_LIT("run/xyz")

// Recorded from the trajectory reader this replaced
static const run_ref_t xyz_refs_plain[] = {
    { 0, {134882.57, 90638.2101, 143800.936}, {57.8339996, 36.5680008, 62.4910011}, {60.1839981, 35.7579994, 94.0210037}, {0, 0, 0, 0, 0, 0} },
    { 5, {136689.2, 89723.1401, 139868.596},  {59.2939987, 35.118, 60.5960007},     {61.223999, 37.7280006, 94.1460037},  {0, 0, 0, 0, 0, 0} },
    { 9, {132960.48, 90949.11, 142797.449},   {57.7290001, 35.8530006, 61.6160011}, {59.769001, 38.0229988, 96.6760025},  {0, 0, 0, 0, 0, 0} },
};
static const run_ref_t xyz_refs_arc[] = {
    { 0, {-33.6560587, 64.7285536, 8.83717018}, {17.9357071, 1.36867702, 0.600875974}, {-12.2290363, -5.84506416, -6.32934189}, {150, 150, 150, 0, 0, 0} },
    { 5, {-36.4959451, 65.8312044, 10.811501},  {18.0747089, 1.47082901, 0.44283101},  {-12.1220617, -5.82219696, -6.06302977}, {150, 150, 150, 0, 0, 0} },
    { 9, {-34.4564086, 69.4655409, 11.9528428}, {17.9201984, 1.24945498, 0.381440997}, {-12.1343231, -5.83449793, -5.91649723}, {150, 150, 150, 0, 0, 0} },
};
static const run_ref_t xyz_refs_o2[] = {
    { 0,    {-2.49199992, -2.296, 0},       {-1.97726095, -1.14954698, -0.000606000016}, {-0.514738977, -1.14645302, 0.000606000016}, {27.9360008, 27.9360008, 27.9360008, 0, 0, 0} },
    { 1000, {-2.49199998, -2.29599994, 0},  {-1.38184595, -1.71667695, -0.222642004},    {-1.11015403, -0.579322994, 0.222642004},    {27.9360008, 27.9360008, 27.9360008, 0, 0, 0} },
    { 1999, {-2.49200004, -2.296, 0},       {-0.668142021, -1.44052899, -0.114527002},   {-1.82385802, -0.855471015, 0.114527002},    {27.9360008, 27.9360008, 27.9360008, 0, 0, 0} },
};
static const run_ref_t xyz_refs_xmol[] = {
    { 0,  {37.6413734, 848.181407, 1665.6992}, {-2.264467, 1.246472, 3.62918711},   {-2.5212729, 11.655159, 7.37770414}, {0, 0, 0, 0, 0, 0} },
    { 25, {363.554169, 515.712128, 1647.947},  {-1.89301896, 1.145679, 3.36118007}, {2.51823092, 21.0426464, 10.6155357}, {0, 0, 0, 0, 0, 0} },
    { 49, {449.09689, 219.117665, 1670.5264},  {-2.252141, 1.22666299, 3.30314994}, {4.11364412, 14.8061962, 11.1623373}, {0, 0, 0, 0, 0, 0} },
};

// A file of several frames as a run: its first frame is the structure's own coordinates, and the
// frames agree with what the reader it replaced gave. A single frame is refused. The number of frames
// the run holds goes to out_frames when given.
//
// @NOTE: void, because the EXPECT and ASSERT macros expand to a bare return.
static void xyz_run_check(int* utest_result, str_t path, const run_ref_t* refs, size_t num_refs, size_t expected_frames, size_t* out_frames) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    size_t frames = 0;
    md_system_t sys = { .alloc = arena };
    md_system_state_t sys_state = { .alloc = arena };
    if (!md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_DISABLE_CACHE_WRITE)) {
        *utest_result = UTEST_TEST_FAILURE;
        goto done;
    }
    if (!md_xyz_system_publish_run(&sys, path, XYZ_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        EXPECT_EQ(0u, expected_frames);
        goto done;
    }
    frames = run_num_frames(&sys, XYZ_RUN);
    EXPECT_EQ(expected_frames, frames);
    if (refs) {
        run_check_refs(utest_result, &sys, XYZ_RUN, expected_frames, sys.atom.count, refs, num_refs);
    }

    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, sys.atom.count);
    EXPECT_TRUE(run_extract_one(&got, &sys, XYZ_RUN, 0));
    EXPECT_EQ(0, MEMCMP(sys_state.xyz, got.xyz, sys.atom.count * sizeof(vec3_t)));

done:
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    if (out_frames) *out_frames = frames;
}

UTEST(xyz, run_matches_reference) {
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz"), xyz_refs_plain, ARRAY_SIZE(xyz_refs_plain), 10, NULL);
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/full.arc"), xyz_refs_arc, ARRAY_SIZE(xyz_refs_arc), 10, NULL);
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/o2.arc"), xyz_refs_o2, ARRAY_SIZE(xyz_refs_o2), 2000, NULL);
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol"), xyz_refs_xmol, ARRAY_SIZE(xyz_refs_xmol), 50, NULL);
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/h2o.arc"), NULL, 0, 2000, NULL);
    xyz_run_check(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/ch4.arc"), NULL, 0, 2000, NULL);
    // Single frames: structures, not runs
    const char* single[] = { "full.xyz", "extended.xyz", "extended1.xyz", "extended2.xyz", "c720.xyz" };
    for (size_t i = 0; i < ARRAY_SIZE(single); ++i) {
        char buf[1024];
        const int len = snprintf(buf, sizeof(buf), MD_UNITTEST_DATA_DIR "/%s", single[i]);
        xyz_run_check(utest_result, (str_t){buf, (size_t)len}, NULL, 0, 0, NULL);
    }
}

// Extended XYZ with a tilted Lattice that changes from frame to frame: each frame's own cell.
UTEST(xyz, run_extended_cell_per_frame) {
    const str_t path = STR_LIT("md_unittest_run_extended.xyz");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    for (int f = 0; f < 3; ++f) {
        md_file_printf(file, "3\n");
        md_file_printf(file, "Lattice=\"%.1f 0.0 0.0 1.5 9.0 0.0 -0.5 0.75 8.0\" Properties=species:S:1:pos:R:3 frame=%d\n", 10.0 + f, f);
        md_file_printf(file, "O %.3f 1.000 2.000\n", 0.5 * f);
        md_file_printf(file, "H 1.000 %.3f 2.500\n", 1.0 + 0.25 * f);
        md_file_printf(file, "H 3.000 1.000 %.3f\n", 2.0 - 0.125 * f);
    }
    md_file_close(&file);

    size_t frames = 0;
    xyz_run_check(utest_result, path, NULL, 0, 3, &frames);
    EXPECT_EQ(3u, frames);

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = arena };
    md_system_state_t sys_state = { .alloc = arena };
    ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_DISABLE_CACHE_WRITE));
    ASSERT_TRUE(md_xyz_system_publish_run(&sys, path, XYZ_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/xyz/unitcell"));
    ASSERT_TRUE(cell && cell->data);
    const float* box = (const float*)cell->data;
    EXPECT_NEAR(12.0f, box[2 * 9 + 0], 1.0e-5f);   // frame 2, a
    EXPECT_NEAR(1.5f,  box[2 * 9 + 3], 1.0e-5f);   // b tilts along x
    EXPECT_NEAR(0.75f, box[2 * 9 + 7], 1.0e-5f);   // c tilts along y
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

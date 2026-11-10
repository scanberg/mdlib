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
#include <md_unitcell.h>
#include <math.h>

UTEST(xyz, xyz_standard) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/C60-Ih.xyz");
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

// The names of the atoms are the symbols as the file writes them, without the zeros which pad them in the parsed data
UTEST(xyz, atom_names) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* files[] = { MD_UNITTEST_DATA_DIR "/C60-Ih.xyz", MD_UNITTEST_DATA_DIR "/tryptophan.xyz", MD_UNITTEST_DATA_DIR "/extended.xyz", MD_UNITTEST_DATA_DIR "/h2o.arc" };
    for (size_t f = 0; f < MD_COUNTOF(files); ++f) {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &state, md_str_from_cstr(files[f]), MD_XYZ_OPTION_NONE));
        ASSERT_GT(sys.atom.count, 0u);
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const md_str_t name = md_atom_name(&sys.atom, i);
            EXPECT_GT(name.len, 0u);
            for (size_t k = 0; k < name.len; ++k) EXPECT_NE(0, (int)name.ptr[k]);
        }
    }
    md_system_t c60 = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_xyz_system_init_from_file(&c60, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/C60-Ih.xyz"), MD_XYZ_OPTION_NONE));
    EXPECT_TRUE(md_str_eq(md_atom_name(&c60.atom, 0), MD_STR_LIT("C")));
    md_arena_allocator_destroy(alloc);
}

UTEST(xyz, xyz_xmol) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/full.xyz");
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
    EXPECT_EQ(2,            data.coordinates[0].particle_type);
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/full.arc");
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
    EXPECT_EQ(2,            data.coordinates[0].particle_type);
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/o2.arc");
    md_xyz_data_t data = {0};
    bool result = md_xyz_data_parse_file(&data, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    EXPECT_EQ(2000, data.num_models);
    EXPECT_EQ(2000 * 2, data.num_coordinates);

    EXPECT_STREQ("O", data.coordinates[0].element_symbol);
    EXPECT_NEAR(-1.977261f, data.coordinates[0].x, 1.0e-5f);
    EXPECT_NEAR(-1.149547f, data.coordinates[0].y, 1.0e-5f);
    EXPECT_NEAR(-0.000606f, data.coordinates[0].z, 1.0e-5f);
    EXPECT_EQ(7,            data.coordinates[0].particle_type);
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/h2o.arc");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/ch4.arc");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/extended.xyz");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/extended1.xyz");
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/extended2.xyz");
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
    md_str_t input = MD_STR_INIT(
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
    md_str_t input = MD_STR_INIT(
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz");

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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/c720.xyz");
    
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
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/nonexistent.xyz");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    bool result = md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_NONE);
    EXPECT_FALSE(result);

    md_system_reset(&sys);
    
    // Test empty path
    md_str_t empty_path = {0};
    result = md_xyz_system_init_from_file(&sys, &sys_state, empty_path, MD_XYZ_OPTION_NONE);
    EXPECT_FALSE(result);
    md_system_free(&sys);
}

// ### RUN ###

#define XYZ_RUN MD_STR_LIT("run/xyz")

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
// @NOTE: void, because the EXPECT and MD_ASSERT macros expand to a bare return.
static void xyz_run_check(int* utest_result, md_str_t path, const run_ref_t* refs, size_t num_refs, size_t expected_frames, size_t* out_frames) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
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
    EXPECT_EQ(0, md_memcmp(sys_state.xyz, got.xyz, sys.atom.count * sizeof(vec3_t)));

done:
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    if (out_frames) *out_frames = frames;
}

UTEST(xyz, run_matches_reference) {
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz"), xyz_refs_plain, MD_COUNTOF(xyz_refs_plain), 10, NULL);
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/full.arc"), xyz_refs_arc, MD_COUNTOF(xyz_refs_arc), 10, NULL);
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/o2.arc"), xyz_refs_o2, MD_COUNTOF(xyz_refs_o2), 2000, NULL);
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol"), xyz_refs_xmol, MD_COUNTOF(xyz_refs_xmol), 50, NULL);
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/h2o.arc"), NULL, 0, 2000, NULL);
    xyz_run_check(utest_result, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/ch4.arc"), NULL, 0, 2000, NULL);
    // Single frames: structures, not runs
    const char* single[] = { "full.xyz", "extended.xyz", "extended1.xyz", "extended2.xyz", "c720.xyz" };
    for (size_t i = 0; i < MD_COUNTOF(single); ++i) {
        char buf[1024];
        const int len = snprintf(buf, sizeof(buf), MD_UNITTEST_DATA_DIR "/%s", single[i]);
        xyz_run_check(utest_result, (md_str_t){buf, (size_t)len}, NULL, 0, 0, NULL);
    }
}

// Extended XYZ with a tilted Lattice that changes from frame to frame: each frame's own cell.
UTEST(xyz, run_extended_cell_per_frame) {
    const md_str_t path = MD_STR_INIT("md_unittest_run_extended.xyz");
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

    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = arena };
    md_system_state_t sys_state = { .alloc = arena };
    ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &sys_state, path, MD_XYZ_OPTION_DISABLE_CACHE_WRITE));
    ASSERT_TRUE(md_xyz_system_publish_run(&sys, path, XYZ_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* cell = md_attributes_find(&sys.attributes, MD_STR_LIT("run/xyz/unitcell"));
    ASSERT_TRUE(cell && cell->data);
    const float* box = (const float*)cell->data;
    EXPECT_NEAR(12.0f, box[2 * 9 + 0], 1.0e-5f);   // frame 2, a
    EXPECT_NEAR(1.5f,  box[2 * 9 + 3], 1.0e-5f);   // b tilts along x
    EXPECT_NEAR(0.75f, box[2 * 9 + 7], 1.0e-5f);   // c tilts along y
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// ### DIALECTS ###

static bool xyz_load_str(md_system_t* sys, md_system_state_t* state, const char* text) {
    return md_xyz_system_init_from_str(sys, state, md_str_from_cstr(text), MD_XYZ_OPTION_NONE);
}

// Plain XYZ with more columns than the element and coordinates: they are not read (nor taken for Tinker's)
UTEST(xyz, plain_extra_columns) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(xyz_load_str(&sys, &st, "2\nwith charges and forces\nC 0.0 0.0 0.0 -0.5 1.0 2.0\nO 1.2 0.0 0.0 0.5 -1.0 0.0\n"));
    ASSERT_EQ(2u, sys.atom.count);
    EXPECT_EQ(MD_Z_C, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, 1));
    EXPECT_EQ(1.2f, st.xyz[1].x);
    EXPECT_EQ(0u, sys.bond.count);

    // Fortran exponents are numbers; what is no number is an error, not a zero
    md_system_t sys2 = {.alloc = arena};
    ASSERT_TRUE(xyz_load_str(&sys2, &st, "1\nfortran\nC 0.15D+01 -2.5d-1 0\n"));
    EXPECT_EQ(1.5f, st.xyz[0].x);
    EXPECT_EQ(-0.25f, st.xyz[0].y);
    md_system_t sys3 = {.alloc = arena};
    EXPECT_FALSE(xyz_load_str(&sys3, &st, "1\nbroken\nC 0.0 abc 0.0\n"));
    md_arena_allocator_destroy(arena);
}

// The element as an atomic number, or a symbol in another case or with a number after it: a name and an element
UTEST(xyz, element_forms) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(xyz_load_str(&sys, &st, "4\nlabels\n8 0 0 0\nC1 1 0 0\nCL2 2 0 0\ncu 3 0 0\n"));
    ASSERT_EQ(4u, sys.atom.count);
    EXPECT_EQ(MD_Z_O,  md_atom_atomic_number(&sys.atom, 0));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 0), MD_STR_LIT("O")));
    EXPECT_EQ(MD_Z_C,  md_atom_atomic_number(&sys.atom, 1));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 1), MD_STR_LIT("C1")));
    EXPECT_EQ(MD_Z_Cl, md_atom_atomic_number(&sys.atom, 2));
    EXPECT_EQ(MD_Z_Cu, md_atom_atomic_number(&sys.atom, 3));
    md_arena_allocator_destroy(arena);
}

// Extended XYZ: the columns are those the Properties say, in any order, the element a symbol or Z
UTEST(xyz, extended_properties_order) {
    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_str(&data, MD_STR_LIT("2\nProperties=pos:R:3:species:S:1 energy=-1.5\n1.0 2.0 3.0 O\n1.5 2.0 3.0 H\n"), md_get_heap_allocator()));
    ASSERT_EQ(2u, data.num_coordinates);
    EXPECT_STREQ("O", data.coordinates[0].element_symbol);
    EXPECT_EQ(8, data.coordinates[0].atomic_number);
    EXPECT_EQ(1.0f, data.coordinates[0].x);
    EXPECT_EQ(3.0f, data.coordinates[0].z);
    EXPECT_STREQ("H", data.coordinates[1].element_symbol);
    md_xyz_data_free(&data, md_get_heap_allocator());

    ASSERT_TRUE(md_xyz_data_parse_str(&data, MD_STR_LIT("2\nProperties=Z:I:1:forces:R:3:pos:R:3\n8 9 9 9 1.0 2.0 3.0\n1 9 9 9 1.5 2.0 3.0\n"), md_get_heap_allocator()));
    ASSERT_EQ(2u, data.num_coordinates);
    EXPECT_EQ(8, data.coordinates[0].atomic_number);
    EXPECT_STREQ("O", data.coordinates[0].element_symbol);
    EXPECT_EQ(1.0f, data.coordinates[0].x);
    EXPECT_EQ(1.5f, data.coordinates[1].x);
    md_xyz_data_free(&data, md_get_heap_allocator());

    // Positions are required
    EXPECT_FALSE(md_xyz_data_parse_str(&data, MD_STR_LIT("1\nProperties=species:S:1:forces:R:3\nO 0 0 0\n"), md_get_heap_allocator()));
    md_xyz_data_free(&data, md_get_heap_allocator());
}

// A Lattice in another orientation than a along x and b in the xy plane (an fcc primitive cell, as ASE writes it): the
// cell is turned into that orientation, and the atoms with it, so that what they are within it is kept
UTEST(xyz, extended_general_lattice) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(xyz_load_str(&sys, &st, "2\nLattice=\"0 1.8 1.8 1.8 0 1.8 1.8 1.8 0\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\nCu 0 0 0\nCu 0.9 0.9 0\n"));
    ASSERT_EQ(2u, sys.atom.count);
    const md_unitcell_t* cell = &st.unitcell;
    double a, b, c, alpha, beta, gamma;
    md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, cell);
    const double len = 1.8 * sqrt(2.0);
    EXPECT_NEAR(len, a, 1e-5);
    EXPECT_NEAR(len, b, 1e-5);
    EXPECT_NEAR(len, c, 1e-5);
    EXPECT_NEAR(60.0, alpha, 1e-4);
    EXPECT_NEAR(60.0, beta, 1e-4);
    EXPECT_NEAR(60.0, gamma, 1e-4);
    // (0.9, 0.9, 0) is a half of c: so it is still
    double A[3][3];
    md_unitcell_A_extract_double(A, cell);
    for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(0.5 * A[2][k], st.xyz[1].elem[k], 1e-5);
    }
    md_arena_allocator_destroy(arena);
}

// Tinker: atom names of a force field, a first line of only the number of atoms, and the connectivity as the bonds
UTEST(xyz, tinker_names_and_bonds) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    const char* files[] = {
        "3  water, tinker\n     1  OW     0.000000    0.000000    0.000000    63     2     3\n     2  HW     0.957200    0.000000    0.000000    64     1\n     3  HW    -0.239988    0.926627    0.000000    64     1\n",
        "3\n     1  O1     0.000000    0.000000    0.000000    63     2     3\n     2  H1     0.957200    0.000000    0.000000    64     1\n     3  H2    -0.239988    0.926627    0.000000    64     1\n",
    };
    for (size_t f = 0; f < MD_COUNTOF(files); ++f) {
        md_system_t sys = {.alloc = arena};
        md_system_state_t st = {.alloc = arena};
        ASSERT_TRUE(xyz_load_str(&sys, &st, files[f]));
        ASSERT_EQ(3u, sys.atom.count);
        EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&sys.atom, 1));
        EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&sys.atom, 2));
        ASSERT_EQ(2u, sys.bond.count);
        EXPECT_EQ(0u, sys.bond.pairs[0].idx[0]);
        EXPECT_EQ(1u, sys.bond.pairs[0].idx[1]);
        EXPECT_EQ(0u, sys.bond.pairs[1].idx[0]);
        EXPECT_EQ(2u, sys.bond.pairs[1].idx[1]);
    }

    md_xyz_data_t data = {0};
    ASSERT_TRUE(md_xyz_data_parse_str(&data, md_str_from_cstr(files[0]), md_get_heap_allocator()));
    EXPECT_TRUE(md_str_eq(data.models[0].comment, MD_STR_LIT("water, tinker")));
    md_xyz_data_free(&data, md_get_heap_allocator());
    md_arena_allocator_destroy(arena);
}

UTEST(xyz, invalid_headers) {
    md_xyz_data_t data = {0};
    EXPECT_FALSE(md_xyz_data_parse_str(&data, MD_STR_LIT("0\nnone\nH 0 0 0\n"), md_get_heap_allocator()));
    md_xyz_data_free(&data, md_get_heap_allocator());
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    EXPECT_FALSE(xyz_load_str(&sys, &st, "0\nnone\nH 0 0 0\n"));
    md_arena_allocator_destroy(arena);
}

// A trajectory of more columns than the coordinates (forces): each frame's positions are those of the columns of pos
UTEST(xyz, run_extra_columns) {
    const md_str_t paths[2] = { MD_STR_INIT("md_unittest_run_forces.xyz"), MD_STR_INIT("md_unittest_run_plain_extra.xyz") };
    for (int p = 0; p < 2; ++p) {
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, paths[p], MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
        for (int f = 0; f < 3; ++f) {
            md_file_printf(file, "2\n");
            if (p == 0) {
                md_file_printf(file, "Lattice=\"10 0 0 0 10 0 0 0 10\" Properties=species:S:1:pos:R:3:forces:R:3\n");
                md_file_printf(file, "O %d.0 1.0 2.0 0.1 0.2 0.3\nH 1.0 %d.5 2.0 -0.1 -0.2 -0.3\n", f, f);
            } else {
                md_file_printf(file, "frame %d\n", f);
                md_file_printf(file, "O %d.0 1.0 2.0 -0.1\nH 1.0 %d.5 2.0 0.2\n", f, f);
            }
        }
        md_file_close(&file);

        size_t frames = 0;
        xyz_run_check(utest_result, paths[p], NULL, 0, 3, &frames);
        md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
        md_system_t sys = { .alloc = arena };
        md_system_state_t st = { .alloc = arena };
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &st, paths[p], MD_XYZ_OPTION_DISABLE_CACHE_WRITE));
        ASSERT_TRUE(md_xyz_system_publish_run(&sys, paths[p], XYZ_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
        md_system_state_t got = {.alloc = arena};
        md_system_state_init(&got, sys.atom.count);
        for (int f = 0; f < 3; ++f) {
            ASSERT_TRUE(run_extract_one(&got, &sys, XYZ_RUN, f));
            EXPECT_EQ((float)f, got.xyz[0].x);
            EXPECT_EQ(1.0f, got.xyz[0].y);
            EXPECT_EQ(2.0f, got.xyz[0].z);
            EXPECT_EQ(1.0f, got.xyz[1].x);
            EXPECT_EQ((float)f + 0.5f, got.xyz[1].y);
        }
        md_system_free(&sys);
        md_vm_arena_destroy(arena);
        remove(paths[p].ptr);
    }
}

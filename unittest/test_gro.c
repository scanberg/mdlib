#include "utest.h"
#include <string.h>

#include <md_gro.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

UTEST(gro, parse_small) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/catalyst.gro");
    md_gro_data_t gro_data = {0};
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(gro_data.num_atoms, 1336);

    md_molecule_t mol = {0};

    md_gro_molecule_init(&mol, &gro_data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }
    md_molecule_free(&mol, alloc);

    EXPECT_TRUE(md_gro_molecule_api()->init_from_file(&mol, path, NULL, alloc));
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }

    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, parse_big) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/centered.gro");
    md_gro_data_t gro_data = { 0 };
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(gro_data.num_atoms, 161742);

    md_molecule_t mol = { 0 };

    md_gro_molecule_init(&mol, &gro_data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }
    md_molecule_free(&mol, alloc);

    EXPECT_TRUE(md_gro_molecule_api()->init_from_file(&mol, path, NULL, alloc));
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }

    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, parse_small_water) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/water.gro");
    md_gro_data_t gro_data = {0};
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(gro_data.num_atoms, 12165);

    // Check that atoms have reasonable coordinates (should be in a 5x5x5 box)
    for (int64_t i = 0; i < 100 && i < gro_data.num_atoms; ++i) { // Test first 100 atoms
        EXPECT_GT(gro_data.atom_data[i].x, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].x, 6.0f);
        EXPECT_GT(gro_data.atom_data[i].y, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].y, 6.0f);
        EXPECT_GT(gro_data.atom_data[i].z, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].z, 6.0f);
    }

    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/nonexistent.gro");
    md_gro_data_t gro_data = {0};
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_FALSE(result);
    
    md_gro_data_free(&gro_data, alloc);
}

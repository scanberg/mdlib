#include "utest.h"
#include <string.h>

#include <md_lammps.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_log.h>

UTEST(lammps, parse_small) {
    md_allocator_i* alloc = md_heap_allocator;

    str_t path = STR(MD_UNITTEST_DATA_DIR"/Water_Ethane_Cubic_Init.data");
    md_lammps_data_t lammps_data = {0};
    data_format_t *formatPtr, format;
    formatPtr = &format;
    get_data_format(formatPtr, full);

    bool result = md_lammps_data_parse_file(&lammps_data, path, alloc, formatPtr);
    EXPECT_TRUE(result);
    EXPECT_EQ(lammps_data.num_atoms, 7800);

    /*

    md_molecule_t mol = {0};

    md_gro_molecule_init(&mol, &gro_data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }
    md_molecule_free(&mol, alloc);

    EXPECT_TRUE(md_gro_molecule_api()->init_from_file(&mol, path, alloc));
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }

    */

    md_lammps_data_free(&lammps_data, alloc);
}

/*
UTEST(gro, parse_big) {
    md_allocator_i* alloc = md_heap_allocator;

    str_t path = STR(MD_UNITTEST_DATA_DIR"/centered.gro");
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

    EXPECT_TRUE(md_gro_molecule_api()->init_from_file(&mol, path, alloc));
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], gro_data.atom_data[i].x * 10.0f);
        EXPECT_EQ(mol.atom.y[i], gro_data.atom_data[i].y * 10.0f);
        EXPECT_EQ(mol.atom.z[i], gro_data.atom_data[i].z * 10.0f);
    }

    md_gro_data_free(&gro_data, alloc);
}

*/
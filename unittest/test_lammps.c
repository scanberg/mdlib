#include "utest.h"
#include <string.h>

#include <md_lammps.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_log.h>

UTEST(lammps, water_ethane_cubic) {
    md_allocator_i* alloc = md_heap_allocator;

    str_t path = STR(MD_UNITTEST_DATA_DIR"/Water_Ethane_Cubic_Init.data");
    md_lammps_data_t data = {0};

    bool result = md_lammps_data_parse_file(&data, path, MD_LAMMPS_ATOM_FORMAT_FULL, alloc);
    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7800);
    EXPECT_EQ(data.num_atom_types, 4);
    EXPECT_EQ(data.num_bonds, 6200);
    EXPECT_EQ(data.num_bond_types, 3);
    EXPECT_EQ(data.num_angles, 8200);
    EXPECT_EQ(data.num_angle_types, 2);
    EXPECT_EQ(data.num_dihedrals, 5400);
    EXPECT_EQ(data.num_dihedral_types, 1);

    EXPECT_EQ(data.atoms[0].id,    1);
    EXPECT_EQ(data.atoms[0].resid, 1);
    EXPECT_EQ(data.atoms[0].type,  2);
    EXPECT_NEAR(data.atoms[0].charge, -0.06824945f,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].x, 0.7171309449132868, 1.0e-5f);
    EXPECT_NEAR(data.atoms[0].y, 20.30016060370651,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].z, 16.45385018655536,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].mass, 12.011f, 1.0e-5f);

    EXPECT_EQ(data.bonds[0].id,   1);
    EXPECT_EQ(data.bonds[0].type, 2);
    EXPECT_EQ(data.bonds[0].atom_id[0], 1);
    EXPECT_EQ(data.bonds[0].atom_id[1], 2);

    EXPECT_EQ(data.cell.xlo, 0);
    EXPECT_EQ(data.cell.ylo, 0);
    EXPECT_EQ(data.cell.zlo, 0);
    EXPECT_NEAR(data.cell.xhi, 39.121263316592f, 1.0e-5);
    EXPECT_NEAR(data.cell.yhi, 39.121263316592f, 1.0e-5);
    EXPECT_NEAR(data.cell.zhi, 39.121263316592f, 1.0e-5);
    EXPECT_EQ(data.cell.xy, 0);
    EXPECT_EQ(data.cell.xz, 0);
    EXPECT_EQ(data.cell.yz, 0);

    md_molecule_t mol = {0};

    md_lammps_molecule_init(&mol, &data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.atoms[i].x);
        EXPECT_EQ(mol.atom.y[i], data.atoms[i].y);
        EXPECT_EQ(mol.atom.z[i], data.atoms[i].z);
        EXPECT_NE(mol.atom.mass[i], 0.0f);
        EXPECT_NE(mol.atom.element[i], 0);
    }

    md_molecule_free(&mol, alloc);

    md_lammps_data_free(&data, alloc);
}

UTEST(lammps, water_ethane_triclinic) {
    md_allocator_i* alloc = md_heap_allocator;

    str_t path = STR(MD_UNITTEST_DATA_DIR"/Water_Ethane_Triclinic_Init.data");
    md_lammps_data_t data = {0};

    bool result = md_lammps_data_parse_file(&data, path, MD_LAMMPS_ATOM_FORMAT_FULL, alloc);
    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7722);
    EXPECT_EQ(data.num_atom_types, 4);
    EXPECT_EQ(data.num_bonds, 6138);
    EXPECT_EQ(data.num_bond_types, 3);
    EXPECT_EQ(data.num_angles, 8118);
    EXPECT_EQ(data.num_angle_types, 2);
    EXPECT_EQ(data.num_dihedrals, 5346);
    EXPECT_EQ(data.num_dihedral_types, 1);

    EXPECT_EQ(data.atoms[0].id,    1);
    EXPECT_EQ(data.atoms[0].resid, 1);
    EXPECT_EQ(data.atoms[0].type,  2);
    EXPECT_NEAR(data.atoms[0].charge, -0.06824945f,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].x, 12.231602469911088, 1.0e-5f);
    EXPECT_NEAR(data.atoms[0].y, 29.301075267966024,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].z, 23.680986746474638,  1.0e-5f);
    EXPECT_NEAR(data.atoms[0].mass, 12.011f, 1.0e-5f);

    EXPECT_EQ(data.bonds[0].id,   1);
    EXPECT_EQ(data.bonds[0].type, 2);
    EXPECT_EQ(data.bonds[0].atom_id[0], 1);
    EXPECT_EQ(data.bonds[0].atom_id[1], 2);

    EXPECT_EQ(data.cell.xlo, 0);
    EXPECT_EQ(data.cell.ylo, 0);
    EXPECT_EQ(data.cell.zlo, 0);
    EXPECT_NEAR(data.cell.xhi, 39.12f, 1.0e-5f);
    EXPECT_NEAR(data.cell.yhi, 35.78331355545549f,  1.0e-5f);
    EXPECT_NEAR(data.cell.zhi, 42.35032810895f,     1.0e-5f);
    EXPECT_NEAR(data.cell.xy,  3.13063427949586f,   1.0e-5f);
    EXPECT_NEAR(data.cell.xz, -7.487709420998051f,  1.0e-5f);
    EXPECT_NEAR(data.cell.yz, -3.1174214811242806f, 1.0e-5f);

    md_molecule_t mol = {0};

    md_lammps_molecule_init(&mol, &data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.atoms[i].x);
        EXPECT_EQ(mol.atom.y[i], data.atoms[i].y);
        EXPECT_EQ(mol.atom.z[i], data.atoms[i].z);
        EXPECT_EQ(mol.atom.mass[i], data.atoms[i].mass);
        EXPECT_EQ(mol.atom.resid[i], data.atoms[i].resid);
        EXPECT_NE(mol.atom.element[i], 0);
    }

    md_molecule_free(&mol, alloc);

    md_lammps_data_free(&data, alloc);
}

UTEST(lammps, read_standardASCII_lammpstrj_cubic) {
    md_allocator_i* alloc = md_heap_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR"/cubic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc);
    ASSERT_TRUE(traj);

    EXPECT_EQ(7800, md_trajectory_num_atoms(traj));
    EXPECT_EQ(10, md_trajectory_num_frames(traj));

    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_temp_allocator, mem_size);
    float* x = (float*)mem_ptr;
    float* y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float* z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;
    md_trajectory_frame_header_t header;
    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(7800, header.num_atoms);
    }
    EXPECT_TRUE(md_trajectory_load_frame(traj, 0, &header, x, y, z));
    EXPECT_NE(x, 0);
    EXPECT_NEAR(header.unit_cell.basis.col[0].x, 39.121262, 0.0001);

    md_free(md_temp_allocator, mem_ptr, mem_size);
    md_lammps_trajectory_free(traj);
}

UTEST(lammps, read_standardASCII_lammpstrj_triclinic) {
    md_allocator_i* alloc = md_heap_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR"/triclinic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc);
    ASSERT_TRUE(traj);

    EXPECT_EQ(7722, md_trajectory_num_atoms(traj));
    EXPECT_EQ(10, md_trajectory_num_frames(traj));
    const int64_t mem_size = md_trajectory_num_atoms(traj) * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_temp_allocator, mem_size);
    float* x = (float*)mem_ptr;
    float* y = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 1;
    float* z = (float*)mem_ptr + md_trajectory_num_atoms(traj) * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < md_trajectory_num_frames(traj); ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(7722, header.num_atoms);
    }

    EXPECT_TRUE(md_trajectory_load_frame(traj, 0, &header, x, y, z));
    EXPECT_NE(x, 0);
    MD_LOG_DEBUG("X value is %f", x);
    MD_LOG_DEBUG("X value is %f", header.unit_cell.basis.elem[0][0]);
    //The coordinates are scaled by a vector, don't know what the correct number should be

    md_free(md_temp_allocator, mem_ptr, mem_size);
    md_lammps_trajectory_free(traj);

}

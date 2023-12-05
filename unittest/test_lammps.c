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
    ASSERT_TRUE(result);
    EXPECT_EQ(lammps_data.num_atoms, 7800);
    EXPECT_EQ(lammps_data.num_atom_types, 4);
    EXPECT_EQ(lammps_data.atom_type_mass[0].mass, 1.008f);
    EXPECT_EQ(lammps_data.bonds[0].second_atom_idx, 2);
    EXPECT_EQ(lammps_data.cell_def.xy, 0);

    md_molecule_t mol = {0};

    md_lammps_molecule_init(&mol, &lammps_data, alloc);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], lammps_data.atom_data[i].x);
        EXPECT_EQ(mol.atom.y[i], lammps_data.atom_data[i].y);
        EXPECT_EQ(mol.atom.z[i], lammps_data.atom_data[i].z);
        EXPECT_NE(mol.atom.mass[i], 0.0f);
        EXPECT_NE(mol.atom.element[i], 0);
    }

    md_molecule_free(&mol, alloc);

    md_lammps_data_free(&lammps_data, alloc);
}

UTEST(lammps, read_frames_triclinic) {
    md_allocator_i* alloc = md_heap_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR"/triclinic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc);
    ASSERT_TRUE(traj);

    EXPECT_EQ(7722, md_trajectory_num_atoms(traj));
    EXPECT_EQ(121, md_trajectory_num_frames(traj));

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
    //The coordinates are scaled by a vector, don't know what the correct number should be

    md_free(md_temp_allocator, mem_ptr, mem_size);
    md_lammps_trajectory_free(traj);

}

UTEST(lammps, read_frames_cubic) {
    md_allocator_i* alloc = md_heap_allocator;
    str_t path = STR(MD_UNITTEST_DATA_DIR"/cubic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc);
    ASSERT_TRUE(traj);

    EXPECT_EQ(7800, md_trajectory_num_atoms(traj));
    EXPECT_EQ(120, md_trajectory_num_frames(traj));

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
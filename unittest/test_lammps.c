#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_lammps.h>
#include <md_trajectory.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_log.h>

#define MAX_VALIDATION_SAMPLES 100

UTEST(lammps, water_ethane_cubic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Cubic_Init.data");
    md_lammps_data_t data = {0};

    const char** atom_formats = md_lammps_atom_format_strings();
    ASSERT_TRUE(atom_formats);

    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    ASSERT_TRUE(atom_format);

    bool result = md_lammps_data_parse_file(&data, path, atom_format, alloc);

    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7800);
    EXPECT_EQ(data.num_atom_types, 4);
    EXPECT_EQ(data.num_bonds, 6200);
    EXPECT_EQ(data.num_bond_types, 3);
    EXPECT_EQ(data.num_angles, 8200);
    EXPECT_EQ(data.num_angle_types, 2);
    EXPECT_EQ(data.num_dihedrals, 5400);
    EXPECT_EQ(data.num_dihedral_types, 1);

    EXPECT_EQ(data.cell.xlo, 0);
    EXPECT_EQ(data.cell.ylo, 0);
    EXPECT_EQ(data.cell.zlo, 0);
    EXPECT_NEAR(data.cell.xhi, 39.121263316592f, 1.0e-5);
    EXPECT_NEAR(data.cell.yhi, 39.121263316592f, 1.0e-5);
    EXPECT_NEAR(data.cell.zhi, 39.121263316592f, 1.0e-5);
    EXPECT_EQ(data.cell.xy, 0);
    EXPECT_EQ(data.cell.xz, 0);
    EXPECT_EQ(data.cell.yz, 0);

    const md_lammps_atom_t ref_atoms[] = {
        {1, 1, 2, -0.06824945, 0.7171309449132868,  20.30016060370651,  16.45385018655536 },
        {2, 1, 2, -0.06824945, 38.92628229188471,   19.676107297642947, 17.524030273606662 },
        {3, 1, 1,  0.02274982, 0.35600433814655547, 21.295980033518028, 16.190304644421158 },
        {4, 1, 1,  0.02274982, 0.7222350352475837,  19.67644034589673,  15.558042260533268 },
        {5, 1, 1,  0.02274982, 1.736920897842471,   20.382201091257382, 16.834406173978778 },
        {6, 1, 1,  0.02274982, 37.906493379107445,  19.594063602913423, 17.143473777787637 },
    };

    for (size_t i = 0; i < ARRAY_SIZE(ref_atoms); ++i) {
        EXPECT_EQ(  ref_atoms[i].id,        data.atoms[i].id);
        EXPECT_EQ(  ref_atoms[i].resid,     data.atoms[i].resid);
        EXPECT_EQ(  ref_atoms[i].type,      data.atoms[i].type);
        EXPECT_NEAR(ref_atoms[i].charge,    data.atoms[i].charge,   1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].x,         data.atoms[i].x,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].y,         data.atoms[i].y,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].z,         data.atoms[i].z,        1.0e-5f);
    }

    const md_lammps_atom_type_t ref_atom_types[] = {
        {1, 1.008f,    1.1f},
        {2, 12.011f,   1.7f},
        {3, 15.9994f,  1.52f},
        {4, 1.008f,    1.1f},
    };

    for (size_t i = 0; i < ARRAY_SIZE(ref_atom_types); ++i) {
        EXPECT_EQ(ref_atom_types[i].id,     data.atom_types[i].id);
        EXPECT_NEAR(ref_atom_types[i].mass, data.atom_types[i].mass,   1.0e-5f);
        EXPECT_NEAR(ref_atom_types[i].radius, data.atom_types[i].radius, 1.0e-5f);
    }

    const md_lammps_bond_t ref_bonds[] = {
        {1, 2, {1, 2}},
        {2, 1, {1, 3}},
        {3, 1, {1, 4}},
        {4, 1, {1, 5}},
        {5, 1, {2, 6}},
        {6, 1, {2, 7}},
    };

    for (size_t i = 0; i < ARRAY_SIZE(ref_bonds); ++i) {
        EXPECT_EQ(ref_bonds[i].id, data.bonds[i].id);
        EXPECT_EQ(ref_bonds[i].type, data.bonds[i].type);
        EXPECT_EQ(ref_bonds[i].atom_id[0], data.bonds[i].atom_id[0]);
        EXPECT_EQ(ref_bonds[i].atom_id[1], data.bonds[i].atom_id[1]);
    }

    const md_lammps_angle_t ref_angles[] = {
        {1, 1, {2, 1, 3}},
        {2, 1, {2, 1, 4}},
        {3, 1, {2, 1, 5}},
        {4, 1, {3, 1, 4}},
        {5, 1, {3, 1, 5}},
        {6, 1, {4, 1, 5}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(ref_angles); ++i) {
        EXPECT_EQ(ref_angles[i].id, data.angles[i].id);
        EXPECT_EQ(ref_angles[i].type, data.angles[i].type);
        EXPECT_EQ(ref_angles[i].atom_id[0], data.angles[i].atom_id[0]);
        EXPECT_EQ(ref_angles[i].atom_id[1], data.angles[i].atom_id[1]);
        EXPECT_EQ(ref_angles[i].atom_id[2], data.angles[i].atom_id[2]);
    }

    const md_lammps_dihedral_t ref_dihedrals[] = {
        {1, 1, {3, 1, 2, 6}},
        {2, 1, {3, 1, 2, 7}},
        {3, 1, {3, 1, 2, 8}},
        {4, 1, {4, 1, 2, 6}},
        {5, 1, {4, 1, 2, 7}},
        {6, 1, {4, 1, 2, 8}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(ref_dihedrals); ++i) {
        EXPECT_EQ(ref_dihedrals[i].id,          data.dihedrals[i].id);
        EXPECT_EQ(ref_dihedrals[i].type,        data.dihedrals[i].type);
        EXPECT_EQ(ref_dihedrals[i].atom_id[0],  data.dihedrals[i].atom_id[0]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[1],  data.dihedrals[i].atom_id[1]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[2],  data.dihedrals[i].atom_id[2]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[3],  data.dihedrals[i].atom_id[3]);
    }

    md_system_t mol = {0};

    md_lammps_molecule_init(&mol, &data, alloc);
    for (size_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.atoms[i].x);
        EXPECT_EQ(mol.atom.y[i], data.atoms[i].y);
        EXPECT_EQ(mol.atom.z[i], data.atoms[i].z);
    }

    // Skip first atom type == Unknown
    for (size_t i = 1; i < mol.atom.type.count; ++i) {
        str_t type_id = md_atom_type_name(&mol.atom.type, i);
        float mass    = md_atom_type_mass(&mol.atom.type, i);
        float radius  = md_atom_type_radius(&mol.atom.type, i);
        md_atomic_number_t z = md_atom_type_atomic_number(&mol.atom.type, i);

        bool found = false;
        for (size_t j = 0; j < data.num_atom_types; ++j) {
            char buf[8];
            int len = snprintf(buf, sizeof(buf), "type_%i", data.atom_types[j].id);
            str_t ref_type_id = {buf, len};
            if (str_eq(type_id, ref_type_id)) {
                float ref_mass = data.atom_types[j].mass;
                float ref_radius = data.atom_types[j].radius;
                md_atomic_number_t ref_z = md_atomic_number_infer_from_mass(ref_mass);

                found = true;
                EXPECT_NEAR(mass, ref_mass, 1.0e-5f);
                EXPECT_NEAR(radius, ref_radius, 1.0e-5f);
                EXPECT_EQ(z, ref_z);
                break;
            }
        }
        EXPECT_TRUE(found);
    }

    md_system_free(&mol, alloc);

    md_lammps_data_free(&data, alloc);
}

UTEST(lammps, water_ethane_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Triclinic_Init.data");
    md_lammps_data_t data = {0};

    const char** atom_formats = md_lammps_atom_format_strings();
    ASSERT_TRUE(atom_formats);

    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    ASSERT_TRUE(atom_format);

    bool result = md_lammps_data_parse_file(&data, path, atom_format, alloc);

    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7722);
    EXPECT_EQ(data.num_atom_types, 4);
    EXPECT_EQ(data.num_bonds, 6138);
    EXPECT_EQ(data.num_bond_types, 3);
    EXPECT_EQ(data.num_angles, 8118);
    EXPECT_EQ(data.num_angle_types, 2);
    EXPECT_EQ(data.num_dihedrals, 5346);
    EXPECT_EQ(data.num_dihedral_types, 1);

    EXPECT_EQ(data.cell.xlo, 0);
    EXPECT_EQ(data.cell.ylo, 0);
    EXPECT_EQ(data.cell.zlo, 0);
    EXPECT_NEAR(data.cell.xhi, 39.12f, 1.0e-5f);
    EXPECT_NEAR(data.cell.yhi, 35.78331355545549f,  1.0e-5f);
    EXPECT_NEAR(data.cell.zhi, 42.35032810895f,     1.0e-5f);
    EXPECT_NEAR(data.cell.xy,  3.13063427949586f,   1.0e-5f);
    EXPECT_NEAR(data.cell.xz, -7.487709420998051f,  1.0e-5f);
    EXPECT_NEAR(data.cell.yz, -3.1174214811242806f, 1.0e-5f);


    const md_lammps_atom_t ref_atoms[] = {
        {1, 1, 2, -0.06824945, 12.231602469911088, 29.301075267966024, 23.680986746474638},
        {2, 1, 2, -0.06824945, 10.905906375367397, 28.52065203980296,  23.693054990351918},
        {3, 1, 1,  0.02274982, 12.296277457009811, 29.94048021774919,  24.56331518414414 },
        {4, 1, 1,  0.02274982, 12.289658698028006, 29.924248444294474, 22.78666604407197 },
        {5, 1, 1,  0.02274982, 13.073614466503752, 28.606431083336908, 23.684196670413105},
        {6, 1, 1,  0.02274982, 10.063894441320446, 29.21529541522476,  23.689841754669413},
    };

    for (size_t i = 0; i < ARRAY_SIZE(ref_atoms); ++i) {
        EXPECT_EQ(  ref_atoms[i].id,        data.atoms[i].id);
        EXPECT_EQ(  ref_atoms[i].resid,     data.atoms[i].resid);
        EXPECT_EQ(  ref_atoms[i].type,      data.atoms[i].type);
        EXPECT_NEAR(ref_atoms[i].charge,    data.atoms[i].charge,   1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].x,         data.atoms[i].x,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].y,         data.atoms[i].y,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].z,         data.atoms[i].z,        1.0e-5f);
    }

    const md_lammps_bond_t ref_bonds[] = {
        {1, 2, {1, 2}},
        {2, 1, {1, 3}},
        {3, 1, {1, 4}},
        {4, 1, {1, 5}},
        {5, 1, {2, 6}},
        {6, 1, {2, 7}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(ref_bonds); ++i) {
        EXPECT_EQ(ref_bonds[i].id, data.bonds[i].id);
        EXPECT_EQ(ref_bonds[i].type, data.bonds[i].type);
        EXPECT_EQ(ref_bonds[i].atom_id[0], data.bonds[i].atom_id[0]);
        EXPECT_EQ(ref_bonds[i].atom_id[1], data.bonds[i].atom_id[1]);
    }

    const md_lammps_angle_t ref_angles[] = {
        {1, 1, {2, 1, 3}},
        {2, 1, {2, 1, 4}},
        {3, 1, {2, 1, 5}},
        {4, 1, {3, 1, 4}},
        {5, 1, {3, 1, 5}},
        {6, 1, {4, 1, 5}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(ref_angles); ++i) {
        EXPECT_EQ(ref_angles[i].id, data.angles[i].id);
        EXPECT_EQ(ref_angles[i].type, data.angles[i].type);
        EXPECT_EQ(ref_angles[i].atom_id[0], data.angles[i].atom_id[0]);
        EXPECT_EQ(ref_angles[i].atom_id[1], data.angles[i].atom_id[1]);
        EXPECT_EQ(ref_angles[i].atom_id[2], data.angles[i].atom_id[2]);
    }

    const md_lammps_dihedral_t ref_dihedrals[] = {
        {1, 1, {3, 1, 2, 6}},
        {2, 1, {3, 1, 2, 7}},
        {3, 1, {3, 1, 2, 8}},
        {4, 1, {4, 1, 2, 6}},
        {5, 1, {4, 1, 2, 7}},
        {6, 1, {4, 1, 2, 8}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(ref_dihedrals); ++i) {
        EXPECT_EQ(ref_dihedrals[i].id,          data.dihedrals[i].id);
        EXPECT_EQ(ref_dihedrals[i].type,        data.dihedrals[i].type);
        EXPECT_EQ(ref_dihedrals[i].atom_id[0],  data.dihedrals[i].atom_id[0]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[1],  data.dihedrals[i].atom_id[1]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[2],  data.dihedrals[i].atom_id[2]);
        EXPECT_EQ(ref_dihedrals[i].atom_id[3],  data.dihedrals[i].atom_id[3]);
    }

    md_system_t mol = {0};

    md_lammps_molecule_init(&mol, &data, alloc);
    for (size_t i = 0; i < mol.atom.count; ++i) {
        EXPECT_EQ(mol.atom.x[i], data.atoms[i].x);
        EXPECT_EQ(mol.atom.y[i], data.atoms[i].y);
        EXPECT_EQ(mol.atom.z[i], data.atoms[i].z);
        EXPECT_NE(mol.atom.type_idx[i], 0);
    }

    md_system_free(&mol, alloc);

    md_lammps_data_free(&data, alloc);
}

UTEST(lammps, read_standardASCII_lammpstrj_cubic) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/cubic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    size_t num_atoms = md_trajectory_num_atoms(traj);
    size_t num_frames = md_trajectory_num_frames(traj);
    EXPECT_EQ(7800, num_atoms);
    EXPECT_EQ(10, num_frames);
    size_t stride = ALIGN_TO(num_atoms, 16);
    const size_t bytes = stride * 3 * sizeof(float);
    void* mem = md_alloc(md_get_temp_allocator(), bytes);
    float* x = (float*)mem + stride * 0;
    float* y = (float*)mem + stride * 1;
    float* z = (float*)mem + stride * 2;

    md_trajectory_frame_header_t header;
    for (size_t i = 0; i < num_frames; ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(7800, header.num_atoms);
    }
    EXPECT_TRUE(md_trajectory_load_frame(traj, 0, &header, x, y, z));
    EXPECT_NE(x[0], 0);
    

    EXPECT_NEAR(x[0], 0.018331 * 39.121262, 0.0001); //Should be about 0.018331 of cell
    EXPECT_NEAR(header.unitcell.x, 39.121262, 0.0001);

    EXPECT_NEAR(y[0], 0.518904 * 39.121262, 0.0001); //Should be about 0.518904 of cell
    EXPECT_NEAR(header.unitcell.y, 39.121262, 0.0001);

    EXPECT_NEAR(z[0], 0.420586 * 39.121262, 0.0001); //Should be about 0.420586 of cell
    EXPECT_NEAR(header.unitcell.z, 39.121262, 0.0001);

    md_free(md_get_temp_allocator(), mem, bytes);
    md_lammps_trajectory_free(traj);
}

UTEST(lammps, read_standardASCII_lammpstrj_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/triclinic_standardASCII.lammpstrj");
    //md_trajectory_loader_i* traj_load = md_lammps_trajectory_loader();
    md_trajectory_i* traj = md_lammps_trajectory_create(path, alloc, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    size_t num_atoms = md_trajectory_num_atoms(traj);
    size_t num_frames = md_trajectory_num_frames(traj);
    EXPECT_EQ(7722, num_atoms);
    EXPECT_EQ(10, num_frames);
    size_t stride = ALIGN_TO(num_atoms, 16);
    const size_t bytes = stride * 3 * sizeof(float);
    void* mem = md_alloc(md_get_temp_allocator(), bytes);
    float* x = (float*)mem + stride * 0;
    float* y = (float*)mem + stride * 1;
    float* z = (float*)mem + stride * 2;

    md_trajectory_frame_header_t header;

    for (size_t i = 0; i < num_frames; ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
        EXPECT_EQ(7722, header.num_atoms);
    }

    EXPECT_TRUE(md_trajectory_load_frame(traj, 0, &header, x, y, z));

    EXPECT_NEAR(12.2316074, x[0], 0.0001); //Should be about 0.35 of cell
    EXPECT_NEAR(39.1199989, header.unitcell.x, 0.0001);

    EXPECT_NEAR(29.3010769, y[0], 0.0001); //Should be about 0.87 of cell
    EXPECT_NEAR(35.7833138, header.unitcell.y, 0.0001);

    EXPECT_NEAR(23.6809902, z[0], 0.0001); //Should be about 0.56 of cell
    EXPECT_NEAR(42.3503265, header.unitcell.z, 0.0001);

    md_free(md_get_temp_allocator(), mem, bytes);
    md_lammps_trajectory_free(traj);

}


UTEST(lammps, comprehensive_data_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test both cubic and triclinic data files comprehensively
    str_t paths[] = {
        STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"),
        STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Triclinic_Init.data")
    };
    
    const char** atom_formats = md_lammps_atom_format_strings();
    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    md_lammps_molecule_loader_arg_t args = md_lammps_molecule_loader_arg(atom_format);
    
    for (int p = 0; p < 2; ++p) {
        md_system_t mol = {0};
        bool result = md_lammps_system_loader()->init_from_file(&mol, paths[p], &args, alloc);
        ASSERT_TRUE(result);
        
        EXPECT_GT(mol.atom.count, 0);
        
        // Validate coordinates are finite for a sample of atoms
        for (int64_t i = 0; i < MIN(MAX_VALIDATION_SAMPLES, mol.atom.count); ++i) {
            EXPECT_FALSE(isnan(mol.atom.x[i]));
            EXPECT_FALSE(isnan(mol.atom.y[i]));
            EXPECT_FALSE(isnan(mol.atom.z[i]));
            EXPECT_FALSE(isinf(mol.atom.x[i]));
            EXPECT_FALSE(isinf(mol.atom.y[i]));
            EXPECT_FALSE(isinf(mol.atom.z[i]));
        }
        
        md_system_free(&mol, alloc);
    }
}

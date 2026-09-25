#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_lammps.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>

#include "run_check.h"

#define MAX_VALIDATION_SAMPLES 100

UTEST(lammps, water_ethane_cubic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Cubic_Init.data");
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

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_lammps_system_init_from_data(&sys, &sys_state, &data);
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, data.atoms[i].x);
        EXPECT_EQ(sys_state.xyz[i].y, data.atoms[i].y);
        EXPECT_EQ(sys_state.xyz[i].z, data.atoms[i].z);
    }

    // Skip first atom type == Unknown
    for (size_t i = 1; i < sys.atom.type.count; ++i) {
        str_t type_id = md_atom_type_name(&sys.atom.type, i);
        float mass    = md_atom_type_mass(&sys.atom.type, i);
        float radius  = md_atom_type_radius(&sys.atom.type, i);
        md_atomic_number_t z = md_atom_type_atomic_number(&sys.atom.type, i);

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

    md_system_free(&sys);
    md_system_state_free(&sys_state);

    md_lammps_data_free(&data, alloc);
}

UTEST(lammps, water_ethane_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Triclinic_Init.data");
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

    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_lammps_system_init_from_data(&sys, &sys_state, &data));
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, data.atoms[i].x);
        EXPECT_EQ(sys_state.xyz[i].y, data.atoms[i].y);
        EXPECT_EQ(sys_state.xyz[i].z, data.atoms[i].z);
        EXPECT_NE(sys.atom.type_idx[i], 0);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);

    md_lammps_data_free(&data, alloc);
}

#define LAMMPS_RUN STR_LIT("run/dump")

// Recorded from the trajectory reader this replaced
static const run_ref_t lammps_refs_cubic[] = {
    { 0, {151603.819, 153255.739, 151588.291}, {0.717131853, 20.3001785, 16.4538536}, {18.2944317, 28.1543579, 25.7364693}, {39.1212633, 39.1212633, 39.1212633, 0, 0, 0} },
    { 5, {152777.457, 152629.797, 152331.6},   {38.5581512, 19.8276329, 16.9380589},  {19.7969627, 29.9873066, 27.0132313}, {39.1212633, 39.1212633, 39.1212633, 0, 0, 0} },
    { 9, {152191.176, 152829.139, 152926.806}, {38.6221504, 19.5720139, 16.8795338},  {19.8374538, 29.1452999, 27.2056694}, {39.1212633, 39.1212633, 39.1212633, 0, 0, 0} },
};
static const run_ref_t lammps_refs_triclinic[] = {
    { 0, {133787.04, 125832.176, 163041.841},  {12.2316074, 29.3010769, 23.6809902}, {17.5067368, 1.37374055, 25.4883747}, {39.12, 35.7833136, 42.3503281, 3.13063428, -7.48770942, -3.11742148} },
    { 5, {134817.816, 126200.858, 163380.641}, {12.568203, 28.2144489, 24.4033165},  {19.1436806, 2.89150667, 24.8202133}, {39.12, 35.7833136, 42.3503281, 3.13063428, -7.48770942, -3.11742148} },
    { 9, {134242.501, 126280.058, 163813.305}, {11.6662769, 28.0301437, 25.3092327}, {19.9502373, 2.64502144, 25.6623058}, {39.12, 35.7833136, 42.3503281, 3.13063428, -7.48770942, -3.11742148} },
};

UTEST(lammps, run_cubic) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, LAMMPS_RUN, 10, 7800, lammps_refs_cubic, ARRAY_SIZE(lammps_refs_cubic));

    // Scaled coordinates in the file: the first atom sits at these fractions of the cell
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 7800);
    ASSERT_TRUE(run_extract_one(&st, &sys, LAMMPS_RUN, 0));
    EXPECT_NEAR(st.xyz[0].x, 0.018331 * 39.121262, 0.0001);
    EXPECT_NEAR(st.xyz[0].y, 0.518904 * 39.121262, 0.0001);
    EXPECT_NEAR(st.xyz[0].z, 0.420586 * 39.121262, 0.0001);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A LAMMPS dump records the TIMESTEP but not dt, so real time is not derivable. The steps are
// reported exactly and the times fall back to ordinals with no unit, rather than the steps being
// handed out as femtoseconds.
UTEST(lammps, frame_steps_and_time_fallback) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/dump/time"));
    const md_attribute_t* step = md_attributes_find(&sys.attributes, STR_LIT("run/dump/step"));
    ASSERT_TRUE(time && step);

    // No unit: the times are ordinals standing in for time, not picoseconds or femtoseconds
    EXPECT_TRUE(md_unit_is_none(time->unit));
    const size_t F = time->format.shape[0];
    for (size_t i = 0; i < F; ++i) {
        EXPECT_NEAR(((const double*)time->data)[i], (double)i, 1.0e-9);
    }
    // Steps are whatever the file said, and non decreasing
    for (size_t i = 1; i < F; ++i) {
        EXPECT_GE(((const int64_t*)step->data)[i], ((const int64_t*)step->data)[i - 1]);
    }

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(lammps, run_triclinic) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/triclinic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, LAMMPS_RUN, 10, 7722, lammps_refs_triclinic, ARRAY_SIZE(lammps_refs_triclinic));

    // One atom: the whole frame is parsed, the atom lines being in no order, and one kept.
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 7722);
    ASSERT_TRUE(run_extract_one(&st, &sys, LAMMPS_RUN, 9));
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, STR_LIT("run/dump/atom/position"));
    float xyz[3];
    md_attribute_slice_t one = md_attribute_slice_2(9, 17);
    ASSERT_EQ(3u, md_attribute_extract_slice_f32(xyz, 3, pos, &one, md_unit_none()));
    EXPECT_EQ(st.xyz[17].x, xyz[0]);
    EXPECT_EQ(st.xyz[17].z, xyz[2]);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(lammps, comprehensive_data_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test both cubic and triclinic data files comprehensively
    str_t paths[] = {
        STR_INIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"),
        STR_INIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Triclinic_Init.data")
    };
    
    const char** atom_formats = md_lammps_atom_format_strings();
    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    
    for (int p = 0; p < 2; ++p) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, paths[p], atom_format));
        
        EXPECT_GT(sys.atom.count, 0);
        
        // Validate coordinates are finite for a sample of atoms
        for (int64_t i = 0; i < MIN(MAX_VALIDATION_SAMPLES, sys.atom.count); ++i) {
            EXPECT_FALSE(isnan(sys_state.xyz[i].x));
            EXPECT_FALSE(isnan(sys_state.xyz[i].y));
            EXPECT_FALSE(isnan(sys_state.xyz[i].z));
            EXPECT_FALSE(isinf(sys_state.xyz[i].x));
            EXPECT_FALSE(isinf(sys_state.xyz[i].y));
            EXPECT_FALSE(isinf(sys_state.xyz[i].z));
        }
        
        md_system_free(&sys);
        md_system_state_free(&sys_state);
    }
}

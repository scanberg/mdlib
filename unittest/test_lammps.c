#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_lammps.h>
#include <md_system.h>
#include <md_run.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <md_unitcell.h>
#include <stdio.h>

#include "run_check.h"

#define MAX_VALIDATION_SAMPLES 100

UTEST(lammps, water_ethane_cubic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Cubic_Init.data");
    md_lammps_data_t data = {0};

    const char** atom_formats = md_lammps_atom_format_strings();
    ASSERT_TRUE(atom_formats);

    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    ASSERT_TRUE(atom_format);

    bool result = md_lammps_data_parse_file(&data, path, atom_format, alloc);

    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7800);
    EXPECT_EQ(data.num_particle_types, 4);
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

    for (size_t i = 0; i < MD_COUNTOF(ref_atoms); ++i) {
        EXPECT_EQ(  ref_atoms[i].id,        data.atoms[i].id);
        EXPECT_EQ(  ref_atoms[i].resid,     data.atoms[i].resid);
        EXPECT_EQ(  ref_atoms[i].type,      data.atoms[i].type);
        EXPECT_NEAR(ref_atoms[i].charge,    data.atoms[i].charge,   1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].x,         data.atoms[i].x,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].y,         data.atoms[i].y,        1.0e-5f);
        EXPECT_NEAR(ref_atoms[i].z,         data.atoms[i].z,        1.0e-5f);
    }

    const md_lammps_atom_type_t ref_particle_types[] = {
        {1, 1.008f,    1.1f},
        {2, 12.011f,   1.7f},
        {3, 15.9994f,  1.52f},
        {4, 1.008f,    1.1f},
    };

    for (size_t i = 0; i < MD_COUNTOF(ref_particle_types); ++i) {
        EXPECT_EQ(ref_particle_types[i].id, data.particle_types[i].id);
        EXPECT_NEAR(ref_particle_types[i].mass, data.particle_types[i].mass, 1.0e-5f);
        EXPECT_NEAR(ref_particle_types[i].radius, data.particle_types[i].radius, 1.0e-5f);
    }

    const md_lammps_bond_t ref_bonds[] = {
        {1, 2, {1, 2}},
        {2, 1, {1, 3}},
        {3, 1, {1, 4}},
        {4, 1, {1, 5}},
        {5, 1, {2, 6}},
        {6, 1, {2, 7}},
    };

    for (size_t i = 0; i < MD_COUNTOF(ref_bonds); ++i) {
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
    for (size_t i = 0; i < MD_COUNTOF(ref_angles); ++i) {
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
    for (size_t i = 0; i < MD_COUNTOF(ref_dihedrals); ++i) {
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

    // Skip first particle type == Unknown
    for (size_t i = 1; i < sys.atom.type.count; ++i) {
        md_str_t type_id = md_particle_type_name(&sys.atom.type, i);
        float mass    = md_particle_type_mass(&sys.atom.type, i);
        float radius  = md_particle_type_radius(&sys.atom.type, i);
        md_atomic_number_t z = md_particle_type_atomic_number(&sys.atom.type, i);

        bool found = false;
        for (size_t j = 0; j < data.num_particle_types; ++j) {
            char buf[8];
            int len = snprintf(buf, sizeof(buf), "type_%i", data.particle_types[j].id);
            md_str_t ref_type_id = {buf, len};
            if (md_str_eq(type_id, ref_type_id)) {
                float ref_mass = data.particle_types[j].mass;
                float ref_radius = data.particle_types[j].radius;
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

    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/Water_Ethane_Triclinic_Init.data");
    md_lammps_data_t data = {0};

    const char** atom_formats = md_lammps_atom_format_strings();
    ASSERT_TRUE(atom_formats);

    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    ASSERT_TRUE(atom_format);

    bool result = md_lammps_data_parse_file(&data, path, atom_format, alloc);

    ASSERT_TRUE(result);
    EXPECT_EQ(data.num_atoms, 7722);
    EXPECT_EQ(data.num_particle_types, 4);
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

    for (size_t i = 0; i < MD_COUNTOF(ref_atoms); ++i) {
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
    for (size_t i = 0; i < MD_COUNTOF(ref_bonds); ++i) {
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
    for (size_t i = 0; i < MD_COUNTOF(ref_angles); ++i) {
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
    for (size_t i = 0; i < MD_COUNTOF(ref_dihedrals); ++i) {
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

#define LAMMPS_RUN MD_STR_LIT("run/dump")

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
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, LAMMPS_RUN, 10, 7800, lammps_refs_cubic, MD_COUNTOF(lammps_refs_cubic));

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
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, MD_STR_LIT("run/dump/time"));
    const md_attribute_t* step = md_attributes_find(&sys.attributes, MD_STR_LIT("run/dump/step"));
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
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/triclinic_standardASCII.lammpstrj"), LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    run_check_refs(utest_result, &sys, LAMMPS_RUN, 10, 7722, lammps_refs_triclinic, MD_COUNTOF(lammps_refs_triclinic));

    // One atom: the whole frame is parsed, the atom lines being in no order, and one kept.
    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, 7722);
    ASSERT_TRUE(run_extract_one(&st, &sys, LAMMPS_RUN, 9));
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, MD_STR_LIT("run/dump/atom/position"));
    float xyz[3];
    md_attribute_slice_t one = md_attribute_slice_2(9, 17);
    ASSERT_EQ(3u, md_attribute_extract_f32(xyz, 3, pos, one, md_unit_none()));
    EXPECT_EQ(st.xyz[17].x, xyz[0]);
    EXPECT_EQ(st.xyz[17].z, xyz[2]);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

UTEST(lammps, comprehensive_data_validation) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Test both cubic and triclinic data files comprehensively
    md_str_t paths[] = {
        MD_STR_INIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"),
        MD_STR_INIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Triclinic_Init.data")
    };
    
    const char** atom_formats = md_lammps_atom_format_strings();
    const char* atom_format = atom_formats[MD_LAMMPS_ATOM_FORMAT_FULL];
    
    for (int p = 0; p < 2; ++p) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, paths[p], atom_format));
        
        EXPECT_GT(sys.atom.count, 0);
        
        // Validate coordinates are finite for a sample of atoms
        for (int64_t i = 0; i < MD_MIN(MAX_VALIDATION_SAMPLES, sys.atom.count); ++i) {
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

// ### DATA FILES ###

static bool lammps_load_str(md_system_t* sys, md_system_state_t* state, const char* text) {
    return md_lammps_system_init_from_str(sys, state, md_str_from_cstr(text), NULL);
}

// Comments after the counts and the keywords, as LAMMPS allows them anywhere after a blank
UTEST(lammps, data_comments) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st,
        "LAMMPS data file # title\n\n"
        "3 atoms # three\n2 bonds\n2 atom types\n1 bond types\n\n"
        "0.0 10.0 xlo xhi # box\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\n"
        "Masses # masses\n\n1 15.9994 # O\n2 1.008 # H\n\n"
        "Atoms # full\n\n1 1 1 -0.8 1.0 1.0 1.0 # O\n2 1 2 0.4 1.9 1.0 1.0\n3 1 2 0.4 0.8 1.9 1.0\n\n"
        "Bonds\n\n1 1 1 2\n2 1 1 3\n"));
    EXPECT_EQ(3u, sys.atom.count);
    EXPECT_EQ(2u, sys.bond.count);
    EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&sys.atom, 1));
    md_arena_allocator_destroy(arena);
}

// Masses is optional: the types are there without it, of no element
UTEST(lammps, data_without_masses) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st,
        "no masses\n\n3 atoms\n2 atom types\n\n0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\n"
        "Atoms # atomic\n\n1 1 1.0 1.0 1.0\n2 2 1.9 1.0 1.0\n3 2 0.8 1.9 1.0\n"));
    ASSERT_EQ(3u, sys.atom.count);
    EXPECT_NE(sys.atom.type_idx[0], sys.atom.type_idx[1]);
    EXPECT_EQ(sys.atom.type_idx[1], sys.atom.type_idx[2]);
    md_arena_allocator_destroy(arena);
}

// Atom ids need not be 1 to N (atoms deleted), nor Masses in the order of the types: bonds are by the ids
UTEST(lammps, data_sparse_ids) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st,
        "sparse\n\n3 atoms\n2 bonds\n2 atom types\n\n0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\n"
        "Masses\n\n2 1.008\n1 15.9994\n\n"
        "Atoms # full\n\n40 7 2 0.4 0.8 1.9 1.0\n10 7 1 -0.8 1.0 1.0 1.0\n25 7 2 0.4 1.9 1.0 1.0\n\n"
        "Bonds\n\n1 1 10 25\n2 1 40 10\n"));
    ASSERT_EQ(3u, sys.atom.count);
    EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(1.0f, st.xyz[0].x);
    EXPECT_EQ(0.8f, st.xyz[2].x);
    ASSERT_EQ(2u, sys.bond.count);
    EXPECT_EQ(0u, sys.bond.pairs[0].idx[0]);
    EXPECT_EQ(1u, sys.bond.pairs[0].idx[1]);
    EXPECT_EQ(0u, sys.bond.pairs[1].idx[0]);
    EXPECT_EQ(2u, sys.bond.pairs[1].idx[1]);
    md_arena_allocator_destroy(arena);
}

// Type labels: in Masses and Atoms in place of the numbers, and the names of the types
UTEST(lammps, data_type_labels) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st,
        "labels\n\n3 atoms\n2 atom types\n\n0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\n"
        "Atom Type Labels\n\n1 OW\n2 HW\n\n"
        "Masses\n\nOW 15.9994\nHW 1.008\n\n"
        "Atoms # full\n\n1 1 OW -0.8 1.0 1.0 1.0\n2 1 HW 0.4 1.9 1.0 1.0\n3 1 2 0.4 0.8 1.9 1.0\n"));
    ASSERT_EQ(3u, sys.atom.count);
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 0), MD_STR_LIT("OW")));
    EXPECT_TRUE(md_str_eq(md_atom_name(&sys.atom, 1), MD_STR_LIT("HW")));
    EXPECT_EQ(sys.atom.type_idx[1], sys.atom.type_idx[2]);
    EXPECT_EQ(MD_Z_O, md_atom_atomic_number(&sys.atom, 0));
    EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&sys.atom, 2));
    md_arena_allocator_destroy(arena);
}

// More than 99 types are as many types (their names are cut, not merged), and molecules past 99 molecules
UTEST(lammps, data_many_types) {
    char* text = md_alloc(md_get_heap_allocator(), 64 * 1024);
    int len = snprintf(text, 64 * 1024, "many\n\n120 atoms\n120 atom types\n\n0 100 xlo xhi\n0 100 ylo yhi\n0 100 zlo zhi\n\nMasses\n\n");
    for (int t = 1; t <= 120; ++t) len += snprintf(text + len, 64 * 1024 - len, "%d %s\n", t, t % 2 ? "12.011" : "1.008");
    len += snprintf(text + len, 64 * 1024 - len, "\nAtoms # molecular\n\n");
    for (int i = 1; i <= 120; ++i) len += snprintf(text + len, 64 * 1024 - len, "%d %d %d %g 1.0 1.0\n", i, i, i, i * 0.5);
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st, text));
    EXPECT_EQ(120u, sys.atom.count);
    EXPECT_EQ(121u, sys.atom.type.count);
    EXPECT_EQ(120u, sys.residue.count);
    for (size_t i = 1; i < sys.atom.count; ++i) EXPECT_NE(sys.atom.type_idx[i - 1], sys.atom.type_idx[i]);
    EXPECT_EQ(105, sys.residue.seq_id[104]);
    md_arena_allocator_destroy(arena);
    md_free(md_get_heap_allocator(), text, 64 * 1024);
}

// Without the style in the Atoms keyword line, the lines tell it: full, molecular, charge or atomic, with or without
// image flags
UTEST(lammps, data_detect_style) {
    const char* head = "detect\n\n2 atoms\n2 atom types\n\n0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\nMasses\n\n1 15.9994\n2 1.008\n\nAtoms\n\n";
    const struct { const char* atoms; md_lammps_atom_format_t format; } cases[] = {
        { "1 1 1 -0.8 1.0 2.0 3.0 0 0 0\n2 1 2 0.4 1.5 2.0 3.0 0 0 0\n", MD_LAMMPS_ATOM_FORMAT_FULL },
        { "1 1 1 -0.8 1.0 2.0 3.0\n2 1 2 0.4 1.5 2.0 3.0\n",             MD_LAMMPS_ATOM_FORMAT_FULL },
        { "1 1 1 1.0 2.0 3.0\n2 1 2 1.5 2.0 3.0\n",                       MD_LAMMPS_ATOM_FORMAT_MOLECULAR },
        { "1 1 -0.8 1.0 2.0 3.0\n2 2 0.4 1.5 2.0 3.0\n",                  MD_LAMMPS_ATOM_FORMAT_CHARGE },
        { "1 1 1.0 2.0 3.0 1 0 -1\n2 2 1.5 2.0 3.0 0 0 0\n",              MD_LAMMPS_ATOM_FORMAT_ATOMIC },
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        char buf[1024];
        const int len = snprintf(buf, sizeof(buf), "%s%s", head, cases[i].atoms);
        EXPECT_EQ((int)cases[i].format, (int)md_lammps_atom_format_from_str((md_str_t){buf, (size_t)len}));
        md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
        md_system_t sys = {.alloc = arena};
        md_system_state_t st = {.alloc = arena};
        EXPECT_TRUE(lammps_load_str(&sys, &st, buf));
        EXPECT_EQ(1.5f, st.xyz[1].x);
        EXPECT_EQ(MD_Z_H, md_atom_atomic_number(&sys.atom, 1));
        md_arena_allocator_destroy(arena);
    }
    // Named styles of the newer releases
    EXPECT_EQ((int)MD_LAMMPS_ATOM_FORMAT_HYBRID, (int)md_lammps_atom_format_from_str(MD_STR_LIT("t\n\n1 atoms\n1 atom types\n\nAtoms # hybrid sphere dipole\n\n1 1 1.0 2.0 3.0 1.0 1.0 0 0 0 0\n")));
}

// Argon and potassium told apart by their masses (Ar is the heavier, out of the order of the elements)
UTEST(lammps, element_from_mass) {
    EXPECT_EQ(MD_Z_Ar, md_atomic_number_infer_from_mass(39.948f));
    EXPECT_EQ(MD_Z_K,  md_atomic_number_infer_from_mass(39.0983f));
    EXPECT_EQ(MD_Z_Co, md_atomic_number_infer_from_mass(58.933f));
    EXPECT_EQ(MD_Z_Ni, md_atomic_number_infer_from_mass(58.693f));
    EXPECT_NEAR(39.948f, md_atomic_number_mass(MD_Z_Ar), 1e-3f);
    EXPECT_NEAR(39.098f, md_atomic_number_mass(MD_Z_K), 1e-3f);
}

// A general triclinic box: the cell turned so that a lies along x, and the atoms with it
UTEST(lammps, data_general_triclinic) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st,
        "general\n\n2 atoms\n1 atom types\n\n0 4 4 avec\n4 0 4 bvec\n4 4 0 cvec\n1 1 1 abc origin\n\nMasses\n\n1 63.546\n\n"
        "Atoms # atomic\n\n1 1 1.0 1.0 1.0\n2 1 5.0 5.0 1.0\n"));
    double a, b, c, alpha, beta, gamma;
    md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &st.unitcell);
    EXPECT_NEAR(4.0 * sqrt(2.0), a, 1e-5);
    EXPECT_NEAR(60.0, gamma, 1e-4);
    EXPECT_NEAR(0.0, st.xyz[0].x, 1e-6);
    // Atom 2 is c from the origin
    double A[3][3];
    md_unitcell_A_extract_double(A, &st.unitcell);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(A[2][k], st.xyz[1].elem[k], 1e-5);
    md_arena_allocator_destroy(arena);
}

// ### DUMPS ###

static bool write_text(md_str_t path, const char* text) {
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    const size_t len = strlen(text);
    const bool ok = md_file_write(file, text, len) == len;
    md_file_close(&file);
    return ok;
}

static const char* lammps_two_atoms = "two\n\n2 atoms\n1 atom types\n\n0.0 10.0 xlo xhi\n0.0 10.0 ylo yhi\n0.0 10.0 zlo zhi\n\nMasses\n\n1 39.948\n\nAtoms # atomic\n\n1 1 1.0 1.0 1.0\n2 1 5.0 5.0 5.0\n";

// The cell of the bounds of a dump's triclinic box: those of the box with its tilts, xy + xz among them
UTEST(lammps, dump_triclinic_bounds) {
    const md_str_t path = MD_STR_LIT("md_unittest_dump_tri.lammpstrj");
    ASSERT_TRUE(write_text(path,
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0.0 13.5 2.0\n0.0 10.0 1.5\n0.0 10.0 0.0\n"
        "ITEM: ATOMS id type x y z\n2 1 5.0 5.0 5.0\n1 1 1.0 1.0 1.0\n"));
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st, lammps_two_atoms));
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, path, LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, 2);
    ASSERT_TRUE(run_extract_one(&got, &sys, LAMMPS_RUN, 0));
    EXPECT_NEAR(10.0, got.unitcell.x, 1e-6);
    EXPECT_NEAR(2.0, got.unitcell.xy, 1e-6);
    EXPECT_NEAR(1.5, got.unitcell.xz, 1e-6);
    EXPECT_EQ(1.0f, got.xyz[0].x);
    EXPECT_EQ(5.0f, got.xyz[1].x);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// UNITS and TIME before the TIMESTEP, boundaries other than periodic, scaled coordinates unwrapped by their image flags
UTEST(lammps, dump_time_boundaries_images) {
    const md_str_t path = MD_STR_LIT("md_unittest_dump_fixed.lammpstrj");
    ASSERT_TRUE(write_text(path,
        "ITEM: UNITS\nreal\nITEM: TIME\n0.0\nITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS ff pp fm\n0.0 10.0\n0.0 10.0\n0.0 10.0\n"
        "ITEM: ATOMS id type xs ys zs ix iy iz\n1 1 0.1 0.1 0.1 1 0 0\n2 1 0.5 0.5 0.5 0 0 0\n"
        "ITEM: TIME\n1.0\nITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS ff pp fm\n0.0 10.0\n0.0 10.0\n0.0 10.0\n"
        "ITEM: ATOMS id type xs ys zs ix iy iz\n1 1 0.2 0.1 0.1 1 0 0\n2 1 0.5 0.6 0.5 0 0 0\n"));
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st, lammps_two_atoms));
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, path, LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    EXPECT_EQ(2u, run_num_frames(&sys, LAMMPS_RUN));
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, 2);
    ASSERT_TRUE(run_extract_one(&got, &sys, LAMMPS_RUN, 1));
    EXPECT_NEAR(12.0f, got.xyz[0].x, 1e-5f);
    EXPECT_NEAR(6.0f, got.xyz[1].y, 1e-5f);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// Unwrapped coordinates are not unwrapped again by the image flags beside them
UTEST(lammps, dump_unwrapped_with_images) {
    const md_str_t path = MD_STR_LIT("md_unittest_dump_unwrapped.lammpstrj");
    ASSERT_TRUE(write_text(path,
        "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0.0 10.0\n0.0 10.0\n0.0 10.0\n"
        "ITEM: ATOMS id type xu yu zu ix iy iz\n1 1 11.0 1.0 1.0 1 0 0\n2 1 5.0 5.0 5.0 0 0 0\n"));
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    ASSERT_TRUE(lammps_load_str(&sys, &st, lammps_two_atoms));
    ASSERT_TRUE(md_lammps_system_publish_run(&sys, path, LAMMPS_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, 2);
    ASSERT_TRUE(run_extract_one(&got, &sys, LAMMPS_RUN, 0));
    EXPECT_EQ(11.0f, got.xyz[0].x);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

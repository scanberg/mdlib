#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_itp.h>
#include <md_gro.h>
#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_unit.h>

// Mesoscale cellulose (Mehandzhiyski et al.) in the layout of cellulose_msc_50.itp: every slice is a
// CC hub bonded to two OC and four IC beads, and the CC hubs are chained. The whole molecule is one
// residue. Written with the file's own clutter: parameters, angles and POSRES behind an #ifdef.
static str_t make_msc_itp(md_allocator_i* alloc, const char* molname, int slices) {
    md_strb_t sb = md_strb_create(alloc);
    md_strb_fmt(&sb, "\n[ defaults ]\n; nbfunc comb-rule\n    1  3  no\n");
    md_strb_fmt(&sb, "[ atomtypes ]\n; name mass charge ptype sigma epsilon\n    IC 3237.780 0.000 A 0.0 0.0\n    OC 3237.780 0.000 A 0.0 0.0\n    CC 3885.336 0.000 A 0.0 0.0\n");
    md_strb_fmt(&sb, "[ nonbond_params ]\n  IC IC 1 1.95000E-00 25.00000E-00\n");
    md_strb_fmt(&sb, "[ moleculetype ]\n; molname nrexcl\n%s\t    4\n\n[ atoms ]\n", molname);
    static const char* names[7] = {"CC", "OC", "OC", "IC", "IC", "IC", "IC"};
    int nr = 1;
    for (int s = 0; s < slices; ++s) {
        for (int b = 0; b < 7; ++b, ++nr) {
            // The first slice is heavier in the real file, which is what the type mass test below relies on
            const double mass = (b == 0) ? (s == 0 ? 3981.330 : 3885.336) : (s == 0 ? 3317.775 : 3237.780);
            md_strb_fmt(&sb, "%6d %5s %6d  MSC %5s %6d  0.00000E+00  %10.3f\n", nr, names[b], 1, names[b], nr, mass);
        }
    }
    md_strb_fmt(&sb, "\n[ bonds ]\n;   ai     aj funct   table          k\n");
    for (int s = 0; s < slices; ++s) {
        const int c = 1 + 7 * s;
        md_strb_fmt(&sb, "%5d %5d 1 1.745 15000.0\n%5d %5d 1 1.745 15000.0\n", c, c + 4, c, c + 6);
        md_strb_fmt(&sb, "%5d %5d 1 1.227 21000.0\n%5d %5d 1 1.227 21000.0\n", c, c + 3, c, c + 5);
        md_strb_fmt(&sb, "%5d %5d 1 1.268 22000.0\n%5d %5d 1 1.268 22000.0\n", c, c + 1, c, c + 2);
        if (s + 1 < slices) md_strb_fmt(&sb, "%5d %5d 1 2.068 500000.0\n", c, c + 7);
    }
    md_strb_fmt(&sb, "[ angles ]\n    4    1    5    1   41.32   65000.000\n");
    md_strb_fmt(&sb, "#ifdef POSRES\n[ position_restraints ]\n    1     1    1000.0 1000.0 1000.0\n#endif\n");
    return md_strb_to_str(sb);
}

// One gro residue per molecule. Slices advance 2.19 nm along x from each molecule's origin, the six
// satellites sit around the hub, and everything is wrapped into the box.
static str_t make_msc_gro(md_allocator_i* alloc, int num_mol, int slices, const float origin[][3], float box[3], int trailing_water) {
    md_strb_t sb = md_strb_create(alloc);
    const int n = num_mol * slices * 7 + trailing_water;
    md_strb_fmt(&sb, "Mesoscale cellulose\n%d\n", n);
    static const char* names[7] = {"CC", "OC", "OC", "IC", "IC", "IC", "IC"};
    static const float off[7][3] = {{0,0,0}, {0,0.4f,-1.2f}, {0,-0.4f,1.2f}, {0,-0.8f,-0.9f}, {0,0.8f,-0.9f}, {0,-0.8f,0.9f}, {0,0.8f,0.9f}};
    int nr = 1;
    for (int m = 0; m < num_mol; ++m) {
        for (int s = 0; s < slices; ++s) {
            for (int b = 0; b < 7; ++b, ++nr) {
                float p[3] = { origin[m][0] + 2.19f * s + off[b][0], origin[m][1] + off[b][1], origin[m][2] + off[b][2] };
                for (int k = 0; k < 3; ++k) p[k] = fmodf(fmodf(p[k], box[k]) + box[k], box[k]);
                md_strb_fmt(&sb, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", m + 1, "MSC", names[b], nr, p[0], p[1], p[2]);
            }
        }
    }
    for (int w = 0; w < trailing_water; ++w, ++nr) {
        md_strb_fmt(&sb, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", num_mol + 1 + w, "W", "W", nr, 1.0f, 1.0f, 1.0f + w);
    }
    md_strb_fmt(&sb, "%10.5f%10.5f%10.5f\n", box[0], box[1], box[2]);
    return md_strb_to_str(sb);
}

static size_t count_flagged(const md_system_t* sys, md_bond_flags_t flag) {
    size_t n = 0;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (sys->bond.flags[i] & flag) n += 1;
    }
    return n;
}

UTEST(itp, parse_sections_and_preprocessor) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    str_t src = STR_LIT(
        "; comment line\n"
        "#define FLEXIBLE\n"
        "[ atomtypes ]\n"
        "  OW   8  15.9994  -0.834  A  0.315 0.636\n"            // name at.num mass charge ptype
        "  HW  HW  1  1.008  0.417  A  0.0 0.0\n"                // name bond_type at.num mass charge ptype
        "  X   12.0  0.0  A  0.0 0.0\n"                         // name mass charge ptype
        "\n"
        "[ moleculetype ]\n"
        "SOL  2\n"
        "[ atoms ]\n"
        "  1  OW  1  SOL  OW   1  -0.834\n"                     // no mass: atomtype fallback
        "  2  HW  1  SOL  HW1  1  +0.417  1.008\n"                 // explicit plus sign
        "  3  HW  1  SOL  HW2  1   0.417  1.008\n"
        "#ifndef FLEXIBLE\n"
        "[ settles ]\n"
        "  1  1  0.09572  0.15139\n"
        "#else\n"
        "[ bonds ]\n"
        "  1  2  1  0.09572 \\\n"                                 // continuation
        "        502416.0\n"
        "  1  3  1  0.09572 502416.0 ; trailing comment\n"
        "#ifdef NOT_DEFINED\n"
        "  2  3  1  0.15139 502416.0\n"
        "#endif\n"
        "#endif\n"
        "[ constraints ]\n"
        "  1  2  1  0.09572\n"                                   // duplicate of a bond
        "\n"
        "[ moleculetype ]\n"
        "VS 1\n"
        "[ atoms ]\n"
        "  1  X  1  VS  A  1  0.0  12.0\n"
        "  2  X  1  VS  B  1  0.0  12.0\n"
        "  3  X  1  VS  M  1  0.0  0.0\n"
        "  4  X  1  VS  N  1  0.0  0.0\n"
        "[ bonds ]\n"
        "  1  2  1\n"
        "[ virtual_sites2 ]\n"
        "  3  1  2  1  0.5\n"
        "[ virtual_sitesn ]\n"
        "  4  2  1  2\n"
        "[ intermolecular_interactions ]\n"
        "[ bonds ]\n"
        "  1  4  1\n"                                            // global indices, must not land in VS
        "[ system ]\n"
        "Water and sites\n"
        "[ molecules ]\n"
        "SOL  2\n"
        "VS   1\n");

    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, src, STR_LIT(""), alloc));

    ASSERT_EQ(md_array_size(data.atomtypes), 3u);
    EXPECT_EQ(data.atomtypes[0].atomic_number, 8);
    EXPECT_NEAR(data.atomtypes[0].mass, 15.9994f, 1e-4f);
    EXPECT_NEAR(data.atomtypes[0].charge, -0.834f, 1e-4f);
    EXPECT_EQ(data.atomtypes[1].atomic_number, 1);
    EXPECT_NEAR(data.atomtypes[1].mass, 1.008f, 1e-4f);
    EXPECT_EQ(data.atomtypes[2].atomic_number, -1);
    EXPECT_NEAR(data.atomtypes[2].mass, 12.0f, 1e-4f);

    ASSERT_EQ(md_array_size(data.moleculetypes), 2u);
    const md_itp_moleculetype_t* sol = &data.moleculetypes[0];
    EXPECT_TRUE(str_eq_cstr(sol->name, "SOL"));
    ASSERT_EQ(md_array_size(sol->atoms), 3u);
    EXPECT_TRUE(str_eq_cstr(sol->atoms[1].name, "HW1"));
    EXPECT_TRUE(sol->atoms[0].has_charge);
    EXPECT_FALSE(sol->atoms[0].has_mass);
    EXPECT_NEAR(sol->atoms[1].charge, 0.417f, 1e-5f);
    // FLEXIBLE is defined: the bonds branch, not settles, and NOT_DEFINED stays out. The constraint repeats 1-2.
    ASSERT_EQ(md_array_size(sol->bonds), 2u);
    EXPECT_EQ(sol->bonds[0].idx[0], 0);
    EXPECT_EQ(sol->bonds[0].idx[1], 1);
    EXPECT_EQ(sol->bonds[1].idx[0], 0);
    EXPECT_EQ(sol->bonds[1].idx[1], 2);

    const md_itp_moleculetype_t* vs = &data.moleculetypes[1];
    // 1-2, site 3 to its first constructing atom 1, site 4 (funct 2, from 1 2) to 1; nothing from intermolecular
    ASSERT_EQ(md_array_size(vs->bonds), 3u);
    EXPECT_EQ(vs->bonds[0].idx[0], 0); EXPECT_EQ(vs->bonds[0].idx[1], 1);
    EXPECT_EQ(vs->bonds[1].idx[0], 0); EXPECT_EQ(vs->bonds[1].idx[1], 2);
    EXPECT_EQ(vs->bonds[2].idx[0], 0); EXPECT_EQ(vs->bonds[2].idx[1], 3);

    EXPECT_TRUE(str_eq_cstr(data.system_name, "Water and sites"));
    ASSERT_EQ(md_array_size(data.molecules), 2u);
    EXPECT_EQ(data.molecules[0].count, 2);

    {
        // settles when FLEXIBLE is not defined
        str_t rigid = STR_LIT("[ moleculetype ]\nSOL 2\n[ atoms ]\n1 OW 1 SOL OW 1\n2 HW 1 SOL HW1 1\n3 HW 1 SOL HW2 1\n"
                              "#ifndef FLEXIBLE\n[ settles ]\n1 1 0.1 0.16\n#else\n[ bonds ]\n2 3 1\n#endif\n");
        md_itp_data_t d = {0};
        ASSERT_TRUE(md_itp_data_parse_str(&d, rigid, STR_LIT(""), alloc));
        ASSERT_EQ(md_array_size(d.moleculetypes[0].bonds), 2u);
        EXPECT_EQ(d.moleculetypes[0].bonds[1].idx[1], 2);
    }

    {
        str_t broken = STR_LIT("#ifdef A\n[ moleculetype ]\nX 1\n");
        md_itp_data_t d = {0};
        EXPECT_FALSE(md_itp_data_parse_str(&d, broken, STR_LIT(""), alloc));
    }

    md_vm_arena_destroy(alloc);
}

// The case this exists for: a coarse grained gro has no bonds and no bond inference that can find them.
UTEST(itp, supplement_coarse_grained_cellulose) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    enum { SLICES = 20, MOLS = 3 };
    float box[3] = {30.0f, 30.0f, 20.0f};
    const float origin[MOLS][3] = {{5, 5, 5}, {12, 15, 10}, {20, 25, 15}};  // 20 slices = 41.6 nm, longer than the box

    str_t gro = make_msc_gro(alloc, MOLS, SLICES, origin, box, 2);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, gro));
    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL);

    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, make_msc_itp(alloc, "MSC", SLICES), STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));

    const size_t bonds_per_mol = SLICES * 7 - 1;
    EXPECT_EQ(count_flagged(&sys, MD_BOND_FLAG_TOPOLOGY), MOLS * bonds_per_mol);
    EXPECT_EQ(sys.bond.count, MOLS * bonds_per_mol);
    // Three fibrils and the two water beads the topology says nothing about
    EXPECT_EQ(md_structure_count(&sys.structure), (size_t)(MOLS + 2));

    // Per atom mass is published, and it varies (first slice heavier)
    const md_attribute_t* mass_attr = md_attributes_find(&sys.attributes, STR_LIT("atom/mass"));
    ASSERT_TRUE(mass_attr != NULL);
    float masses[SLICES * 7 * MOLS + 2];
    ASSERT_EQ(md_attribute_extract_f32(masses, ARRAY_SIZE(masses), mass_attr, md_unit_dalton()), ARRAY_SIZE(masses));
    EXPECT_NEAR(masses[0], 3981.33f, 0.01f);
    EXPECT_NEAR(masses[7], 3885.336f, 0.01f);

    // Unwrap now makes each fibril whole: consecutive hubs 21.9 A apart along x
    md_util_unwrap_system(&state, &sys);
    for (int m = 0; m < MOLS; ++m) {
        for (int s = 0; s + 1 < SLICES; ++s) {
            const size_t a = (size_t)(m * SLICES * 7 + s * 7);
            EXPECT_NEAR(state.xyz[a + 7].x - state.xyz[a].x, 21.9f, 0.05f);
            EXPECT_NEAR(state.xyz[a + 7].y - state.xyz[a].y, 0.0f, 0.05f);
        }
    }

    // Re-inferring bonds keeps the topology's and adds nothing on top of them
    md_util_infer_covalent_bonds(&sys.bond, &state, &sys, sys.alloc);
    md_bond_build_connectivity(&sys.bond, sys.atom.count, sys.alloc);
    EXPECT_EQ(count_flagged(&sys, MD_BOND_FLAG_TOPOLOGY), MOLS * bonds_per_mol);

    md_vm_arena_destroy(alloc);
}

// Atom names repeat along a polymer, so two consecutive 5-slice molecules carry exactly the names of one
// 10-slice molecule. Residue boundaries have to decide, and scanning longest first must not glue them.
UTEST(itp, residue_boundaries_disambiguate_polymers) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    float box[3] = {100.0f, 100.0f, 100.0f};
    const float origin[2][3] = {{10, 10, 10}, {10, 60, 10}};
    str_t gro = make_msc_gro(alloc, 2, 5, origin, box, 0);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, gro));

    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, make_msc_itp(alloc, "MSC10", 10), STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_data_parse_str(&data, make_msc_itp(alloc, "MSC5", 5), STR_LIT(""), alloc));
    ASSERT_EQ(md_array_size(data.moleculetypes), 2u);

    md_itp_instance_t* inst = 0;
    ASSERT_EQ(md_itp_match_system(&inst, &data, &sys, alloc), 2u);
    EXPECT_EQ(inst[0].type, 1u);
    EXPECT_EQ(inst[0].atom_offset, 0u);
    EXPECT_EQ(inst[1].type, 1u);
    EXPECT_EQ(inst[1].atom_offset, 35u);

    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));
    EXPECT_EQ(md_structure_count(&sys.structure), 2u);

    md_vm_arena_destroy(alloc);
}

// A .top lays molecules down in order. When its [ molecules ] does not fit, placement falls back to scanning.
UTEST(itp, molecules_section_sequential_and_fallback) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    float box[3] = {100.0f, 100.0f, 100.0f};
    const float origin[2][3] = {{10, 10, 10}, {10, 60, 10}};
    str_t gro = make_msc_gro(alloc, 2, 3, origin, box, 1);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, gro));

    str_t water = STR_LIT("[ moleculetype ]\nW 1\n[ atoms ]\n1 W 1 W W 1 0.0 72.0\n");
    {
        md_strb_t top = md_strb_create(alloc);
        md_strb_push_str(&top, make_msc_itp(alloc, "MSC", 3));
        md_strb_push_str(&top, water);
        md_strb_push_cstr(&top, "[ system ]\ntest\n[ molecules ]\nMSC 2\nW 1\n");
        md_itp_data_t data = {0};
        ASSERT_TRUE(md_itp_data_parse_str(&data, md_strb_to_str(top), STR_LIT(""), alloc));
        md_itp_instance_t* inst = 0;
        ASSERT_EQ(md_itp_match_system(&inst, &data, &sys, alloc), 3u);
        EXPECT_EQ(inst[2].type, 1u);
        EXPECT_EQ(inst[2].atom_offset, 42u);
    }
    {
        // Wrong order: W first cannot match atom 0, so the scan places them anyway
        md_strb_t top = md_strb_create(alloc);
        md_strb_push_str(&top, make_msc_itp(alloc, "MSC", 3));
        md_strb_push_str(&top, water);
        md_strb_push_cstr(&top, "[ molecules ]\nW 1\nMSC 2\n");
        md_itp_data_t data = {0};
        ASSERT_TRUE(md_itp_data_parse_str(&data, md_strb_to_str(top), STR_LIT(""), alloc));
        md_itp_instance_t* inst = 0;
        EXPECT_EQ(md_itp_match_system(&inst, &data, &sys, alloc), 3u);
    }
    {
        // Nothing matches
        md_itp_data_t data = {0};
        ASSERT_TRUE(md_itp_data_parse_str(&data, STR_LIT("[ moleculetype ]\nFOO 1\n[ atoms ]\n1 X 1 FOO X 1\n"), STR_LIT(""), alloc));
        EXPECT_FALSE(md_itp_system_supplement(&sys, &data));
    }

    md_vm_arena_destroy(alloc);
}

// Bonds outside matched molecules and user defined bonds survive; inferred bonds inside are replaced.
UTEST(itp, supplement_replaces_only_covered_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    float box[3] = {100.0f, 100.0f, 100.0f};
    const float origin[1][3] = {{10, 10, 10}};
    str_t gro = make_msc_gro(alloc, 1, 2, origin, box, 2);  // 14 beads, then W at 14 and 15
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, gro));

    md_bond_data_clear(&sys.bond);
    md_system_bond_insert(&sys, 1, 2, MD_BOND_FLAG_COVALENT);       // inside: replaced
    md_system_bond_insert(&sys, 14, 15, MD_BOND_FLAG_COVALENT);     // outside: kept
    md_system_bond_insert(&sys, 3, 4, MD_BOND_FLAG_USER_DEFINED);   // user: kept, and last
    md_bond_build_connectivity(&sys.bond, sys.atom.count, sys.alloc);

    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, make_msc_itp(alloc, "MSC", 2), STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));

    EXPECT_EQ(sys.bond.count, 1u + 13u + 1u);
    EXPECT_EQ(sys.bond.pairs[0].idx[0], 14);
    EXPECT_EQ((int)sys.bond.flags[sys.bond.count - 1], (int)MD_BOND_FLAG_USER_DEFINED);
    EXPECT_EQ(md_bond_find(&sys.bond, 1, 2), -1);
    EXPECT_NE(md_bond_find(&sys.bond, 0, 7), -1);

    md_vm_arena_destroy(alloc);
}

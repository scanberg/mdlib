#include "utest.h"

#include <md_system.h>
#include <md_topology_build.h>
#include <md_attributes.h>
#include <md_pdb.h>
#include <md_xyz.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_unit.h>

// Joining two systems of the same atoms (md_system_match, md_system_merge): the rules of md_io.h, LOAD, kind by kind,
// on systems built in memory. Two waters: a topology of them, a PDB file of them, an XYZ file of them, and their
// coordinates alone (as a restart has them).

typedef struct merge_test_t {
    md_allocator_i* arena;
} merge_test_t;

// The heap rather than an arena: what a join frees, it frees for real, so that a sanitizer sees a double free or a leak
UTEST_F_SETUP(merge_test_t) {
    utest_fixture->arena = md_get_heap_allocator();
}

UTEST_F_TEARDOWN(merge_test_t) {
    (void)utest_fixture;
}

static const md_atom_pair_t water_bonds[4] = { {{0, 1}}, {{0, 2}}, {{3, 4}}, {{3, 5}} };

// The topology of two waters: atoms without coordinates, their bonds complete
static bool water_topology(md_system_t* sys, md_system_state_t* state) {
    static const md_str_t name[6] = { MD_STR_INIT("OW"), MD_STR_INIT("HW1"), MD_STR_INIT("HW2"), MD_STR_INIT("OW"), MD_STR_INIT("HW1"), MD_STR_INIT("HW2") };
    static const md_str_t res[2] = { MD_STR_INIT("SOL"), MD_STR_INIT("SOL") };
    static const md_str_t chain[2] = { MD_STR_INIT("W"), MD_STR_INIT("W") };
    static const uint32_t res_offset[3] = { 0, 3, 6 };
    static const double charge[6] = { -0.834, 0.417, 0.417, -0.834, 0.417, 0.417 };
    static const double mass[6] = { 15.999, 1.008, 1.008, 15.999, 1.008, 1.008 };
    const md_topology_desc_t desc = {
        .format = "test topology",
        .num_atoms = 6, .atom_name = name, .charge = charge, .mass = mass,
        .num_residues = 2, .res_name = res, .res_offset = res_offset, .res_chain = chain,
        .num_bonds = 4, .bonds = water_bonds,
    };
    return md_topology_system_load(sys, state, &desc, MD_STR_LIT("water.top")) == MD_LOAD_OK;
}

// The PDB file of them: names cut short as some writers do, the water named HOH, one bond given by CONECT
static const char* water_pdb =
    "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
    "ATOM      1  OW  HOH A   1       1.000   1.000   1.000  1.00  0.00           O\n"
    "ATOM      2  HW1 HOH A   1       1.957   1.000   1.000  1.00  0.00           H\n"
    "ATOM      3  HW2 HOH A   1       0.760   1.927   1.000  1.00  0.00           H\n"
    "ATOM      4  OW  HOH A   2       5.000   5.000   5.000  1.00  0.00           O\n"
    "ATOM      5  HW1 HOH A   2       5.957   5.000   5.000  1.00  0.00           H\n"
    "ATOM      6  HW2 HOH A   2       4.760   5.927   5.000  1.00  0.00           H\n"
    "CONECT    4    5\n"
    "END\n";

static const char* water_xyz =
    "6\n"
    "two waters\n"
    "O  1.000  1.000  1.000\n"
    "H  1.957  1.000  1.000\n"
    "H  0.760  1.927  1.000\n"
    "O  5.000  5.000  5.000\n"
    "H  5.957  5.000  5.000\n"
    "H  4.760  5.927  5.000\n";

// The coordinates alone, as a restart has them: a system known by its count. Spread far apart unless near.
static void water_points(md_system_t* sys, md_system_state_t* state, bool near) {
    md_system_reset(sys);
    sys->atom.count = 6;
    md_system_state_init(state, 6);
    for (int i = 0; i < 6; ++i) {
        const float s = near ? 1.0f : 15.0f;
        state->xyz[i] = (md_float3_t){ { 1.0f + s * (float)i, 2.0f, 3.0f } };
    }
}

UTEST_F(merge_test_t, onto_no_atoms_is_a_copy) {
    md_system_t a = { .alloc = utest_fixture->arena }, b = { .alloc = utest_fixture->arena };
    md_system_state_t as = { .alloc = utest_fixture->arena }, bs = { .alloc = utest_fixture->arena };
    ASSERT_TRUE(water_topology(&b, &bs));
    uint32_t override = ~0u;
    EXPECT_TRUE(md_system_match(&override, &a, &as, &b, &bs));
    EXPECT_EQ(0u, override);
    ASSERT_TRUE(md_system_merge(&a, &as, &b, &bs));
    EXPECT_EQ(6u, a.atom.count);
    EXPECT_EQ(4u, md_bond_list_count(&a.given_bond));
    EXPECT_TRUE(md_attributes_find(&a.attributes, MD_STR_LIT("atom/charge")) != NULL);
    md_system_free(&a); md_system_state_free(&as);
    md_system_free(&b); md_system_state_free(&bs);
}

// ATOMS: elements and names where both give them
UTEST_F(merge_test_t, same_atoms) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t top = { .alloc = alloc }, pdb = { .alloc = alloc }, xyz = { .alloc = alloc };
    md_system_state_t top_s = { .alloc = alloc }, pdb_s = { .alloc = alloc }, xyz_s = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&top, &top_s));
    ASSERT_TRUE(md_pdb_system_init_from_str(&pdb, &pdb_s, md_str_from_cstr(water_pdb), NULL));
    ASSERT_TRUE(md_xyz_system_init_from_str(&xyz, &xyz_s, md_str_from_cstr(water_xyz), NULL));

    EXPECT_TRUE(md_system_match(NULL, &top, &top_s, &pdb, &pdb_s));
    EXPECT_TRUE(md_system_match(NULL, &pdb, &pdb_s, &top, &top_s));
    EXPECT_TRUE(md_system_match(NULL, &top, &top_s, &xyz, &xyz_s));    // by their elements

    // An oxygen where the topology has a hydrogen
    md_system_t swapped = { .alloc = alloc };
    md_system_state_t swapped_s = { .alloc = alloc };
    const char* text = "6\n\nH 1 1 1\nO 2 1 1\nH 1 2 1\nO 5 5 5\nH 6 5 5\nH 5 6 5\n";
    ASSERT_TRUE(md_xyz_system_init_from_str(&swapped, &swapped_s, md_str_from_cstr(text), NULL));
    EXPECT_FALSE(md_system_match(NULL, &top, &top_s, &swapped, &swapped_s));

    // Fewer atoms
    md_system_t three = { .alloc = alloc };
    md_system_state_t three_s = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&three, &three_s, MD_STR_LIT("3\n\nO 1 1 1\nH 2 1 1\nH 1 2 1\n"), NULL));
    EXPECT_FALSE(md_system_match(NULL, &top, &top_s, &three, &three_s));

    // Not joined: nothing changed
    const uint64_t version = top.topology_version;
    EXPECT_FALSE(md_system_merge(&top, &top_s, &swapped, &swapped_s));
    EXPECT_EQ(version, top.topology_version);
    EXPECT_FALSE(md_system_state_has_coords(&top_s));

    md_system_free(&top); md_system_state_free(&top_s); md_system_free(&pdb); md_system_state_free(&pdb_s); md_system_free(&xyz); md_system_state_free(&xyz_s); md_system_free(&swapped); md_system_state_free(&swapped_s); md_system_free(&three); md_system_state_free(&three_s);
}

// ATOMS by their count alone: the topology's bonds in the coordinates
UTEST_F(merge_test_t, by_count_the_bonds_decide) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t top = { .alloc = alloc }, pts = { .alloc = alloc };
    md_system_state_t top_s = { .alloc = alloc }, pts_s = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&top, &top_s));

    water_points(&pts, &pts_s, false);
    EXPECT_FALSE(md_system_match(NULL, &top, &top_s, &pts, &pts_s));

    water_points(&pts, &pts_s, true);
    uint32_t override = ~0u;
    EXPECT_TRUE(md_system_match(&override, &top, &top_s, &pts, &pts_s));
    EXPECT_EQ(0u, override);
    ASSERT_TRUE(md_system_merge(&top, &top_s, &pts, &pts_s));
    // COORDINATES filled; the atoms are the topology's
    ASSERT_TRUE(md_system_state_has_coords(&top_s));
    EXPECT_EQ(4.0f, top_s.xyz[3].x);
    EXPECT_EQ(2u, top.residue.count);
    EXPECT_EQ(4u, md_bond_list_count(&top.given_bond));

    md_system_free(&top); md_system_state_free(&top_s); md_system_free(&pts); md_system_state_free(&pts_s);
}

// The points first: the topology names their atoms, and what the points had stays theirs
UTEST_F(merge_test_t, points_first) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t pts = { .alloc = alloc }, top = { .alloc = alloc };
    md_system_state_t pts_s = { .alloc = alloc }, top_s = { .alloc = alloc };
    water_points(&pts, &pts_s, true);
    const float value[6] = { 1, 2, 3, 4, 5, 6 };
    md_attributes_publish_atom_column(&pts.attributes, MD_STR_LIT("points/value"), md_unit_none(), 1, value, 6);
    ASSERT_TRUE(md_system_bond_insert(&pts, 0, 3, MD_BOND_FLAG_NONE));
    ASSERT_EQ(1u, md_bond_list_count(&pts.bond_edit));
    ASSERT_TRUE(water_topology(&top, &top_s));

    ASSERT_TRUE(md_system_merge(&pts, &pts_s, &top, &top_s));
    EXPECT_EQ(6u, pts.atom.count);
    EXPECT_EQ(2u, pts.residue.count);
    EXPECT_TRUE(md_str_eq_cstr(md_atom_name(&pts.atom, 1), "HW1"));
    EXPECT_EQ(4u, md_bond_list_count(&pts.given_bond));
    EXPECT_EQ(1u, md_bond_list_count(&pts.bond_edit));       // the user's
    EXPECT_EQ(4.0f, pts_s.xyz[3].x);                         // the points' coordinates
    EXPECT_TRUE(md_attributes_find(&pts.attributes, MD_STR_LIT("points/value")) != NULL);
    EXPECT_TRUE(md_attributes_find(&pts.attributes, MD_STR_LIT("atom/charge")) != NULL);

    md_system_free(&pts); md_system_state_free(&pts_s); md_system_free(&top); md_system_state_free(&top_s);
}

// BONDS: complete over complete replaces, given in part only between atoms whose bonds are not complete
UTEST_F(merge_test_t, bonds) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t top = { .alloc = alloc }, pdb = { .alloc = alloc };
    md_system_state_t top_s = { .alloc = alloc }, pdb_s = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&pdb, &pdb_s, md_str_from_cstr(water_pdb), NULL));
    ASSERT_EQ(1u, md_bond_list_count(&pdb.given_bond));
    ASSERT_TRUE(water_topology(&top, &top_s));

    // The structure first: its CONECT bond, of atoms the topology completes, gives way to the topology's
    ASSERT_TRUE(md_system_merge(&pdb, &pdb_s, &top, &top_s));
    EXPECT_EQ(4u, md_bond_list_count(&pdb.given_bond));

    // The topology first: the CONECT bond is of atoms whose bonds are complete, and is not added
    md_system_t top2 = { .alloc = alloc }, pdb2 = { .alloc = alloc };
    md_system_state_t top2_s = { .alloc = alloc }, pdb2_s = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&top2, &top2_s));
    ASSERT_TRUE(md_pdb_system_init_from_str(&pdb2, &pdb2_s, md_str_from_cstr(water_pdb), NULL));
    uint32_t override = ~0u;
    EXPECT_TRUE(md_system_match(&override, &top2, &top2_s, &pdb2, &pdb2_s));
    EXPECT_EQ(0u, override);
    ASSERT_TRUE(md_system_merge(&top2, &top2_s, &pdb2, &pdb2_s));
    EXPECT_EQ(4u, md_bond_list_count(&top2.given_bond));

    // Complete over complete: replaced, which the match says
    md_system_t top3 = { .alloc = alloc };
    md_system_state_t top3_s = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&top3, &top3_s));
    EXPECT_TRUE(md_system_match(&override, &top2, &top2_s, &top3, &top3_s));
    EXPECT_TRUE(override & MD_SYSTEM_KIND_BONDS);
    EXPECT_TRUE(override & MD_SYSTEM_KIND_ATTRIBUTES);
    EXPECT_FALSE(override & MD_SYSTEM_KIND_COORDINATES);

    md_system_free(&top); md_system_state_free(&top_s); md_system_free(&pdb); md_system_state_free(&pdb_s); md_system_free(&top2); md_system_state_free(&top2_s); md_system_free(&pdb2); md_system_state_free(&pdb2_s); md_system_free(&top3); md_system_state_free(&top3_s);
}

// COORDINATES and LABELS: the structure's coordinates and cell, the topology's chains where the structure gives none;
// a second structure replaces the coordinates
UTEST_F(merge_test_t, coordinates_and_labels) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t xyz = { .alloc = alloc }, pdb = { .alloc = alloc }, top = { .alloc = alloc };
    md_system_state_t xyz_s = { .alloc = alloc }, pdb_s = { .alloc = alloc }, top_s = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&xyz, &xyz_s, md_str_from_cstr(water_xyz), NULL));
    ASSERT_TRUE(md_pdb_system_init_from_str(&pdb, &pdb_s, md_str_from_cstr(water_pdb), NULL));
    ASSERT_TRUE(water_topology(&top, &top_s));

    // The XYZ file knows its atoms by their elements: the topology names them, and gives its chains
    ASSERT_EQ(0u, xyz.instance.count);
    ASSERT_TRUE(md_system_merge(&xyz, &xyz_s, &top, &top_s));
    EXPECT_EQ(2u, xyz.residue.count);
    EXPECT_EQ(1u, xyz.instance.count);
    EXPECT_EQ(1.957f, xyz_s.xyz[1].x);

    // The PDB file onto it: its coordinates and cell replace those there, which the match says
    uint32_t override = 0;
    EXPECT_TRUE(md_system_match(&override, &xyz, &xyz_s, &pdb, &pdb_s));
    EXPECT_TRUE(override & MD_SYSTEM_KIND_COORDINATES);
    ASSERT_TRUE(md_system_merge(&xyz, &xyz_s, &pdb, &pdb_s));
    EXPECT_NEAR(20.0f, xyz_s.unitcell.x, 1e-4f);
    // The names stay those it had
    EXPECT_TRUE(md_str_eq_cstr(md_system_residue_name(&xyz, 0), "SOL"));

    md_system_free(&xyz); md_system_state_free(&xyz_s); md_system_free(&pdb); md_system_state_free(&pdb_s); md_system_free(&top); md_system_state_free(&top_s);
}

static size_t constant_provider(void* dst, size_t cap, const md_attributes_t* table, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)table; (void)attr; (void)slice; (void)io;
    float* out = (float*)dst;
    for (size_t i = 0; i < cap; ++i) out[i] = *(const float*)user_data;
    return cap;
}

// ATTRIBUTES: of every kind, copied into the system's own table; the same path replaced
UTEST_F(merge_test_t, attributes) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t a = { .alloc = alloc }, b = { .alloc = alloc };
    md_system_state_t as = { .alloc = alloc }, bs = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&a, &as));
    ASSERT_TRUE(water_topology(&b, &bs));

    const md_str_t text[2] = { MD_STR_INIT("first"), MD_STR_INIT("second") };
    md_attributes_create(&b.attributes, &(md_attribute_desc_t){
        .path = MD_STR_LIT("test/text"), .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = {2} },
        .data = text, .byte_size = sizeof(text) });
    float* constant = md_attributes_alloc_user_data(&b.attributes, sizeof(float));
    *constant = 7.0f;
    const md_attribute_virtual_t virt = { .provider = constant_provider, .user_data = constant, .user_data_size = sizeof(float) };
    const md_attribute_id_t computed = md_attributes_create(&b.attributes, &(md_attribute_desc_t){
        .path = MD_STR_LIT("test/computed"), .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = {6} },
        .virt = &virt });
    md_attributes_alias(&b.attributes, computed, MD_STR_LIT("test/alias"), (md_str_t){0}, (md_str_t){0});

    ASSERT_TRUE(md_system_merge(&a, &as, &b, &bs));
    md_system_free(&b); md_system_state_free(&bs);     // what a has is its own

    const md_attribute_t* t = md_attributes_find(&a.attributes, MD_STR_LIT("test/text"));
    ASSERT_TRUE(t != NULL);
    EXPECT_TRUE(md_str_eq_cstr(md_attribute_str(&a.attributes, t, 1), "second"));
    const md_attribute_t* c = md_attributes_find(&a.attributes, MD_STR_LIT("test/alias"));
    ASSERT_TRUE(c != NULL);
    EXPECT_TRUE(md_attribute_is_alias(c));
    float values[6] = { 0 };
    EXPECT_EQ(6u, md_attribute_extract_f32(values, 6, &a.attributes, c, md_attribute_slice_all(), md_unit_none()));
    EXPECT_EQ(7.0f, values[5]);

    md_system_free(&a); md_system_state_free(&as);
}

// GROUPS: added, with the roles the system gives none of
UTEST_F(merge_test_t, groups) {
    md_allocator_i* alloc = utest_fixture->arena;
    md_system_t a = { .alloc = alloc }, b = { .alloc = alloc };
    md_system_state_t as = { .alloc = alloc }, bs = { .alloc = alloc };
    ASSERT_TRUE(water_topology(&a, &as));
    ASSERT_TRUE(water_topology(&b, &bs));
    md_system_group_add(&a, MD_STR_LIT("first"), md_bitset_range(0, 3), MD_GROUP_ORIGIN_GIVEN);
    const uint32_t g[2] = { (uint32_t)md_system_group_add(&b, MD_STR_LIT("first"), md_bitset_range(0, 3), MD_GROUP_ORIGIN_GIVEN),
                            (uint32_t)md_system_group_add(&b, MD_STR_LIT("second"), md_bitset_range(3, 6), MD_GROUP_ORIGIN_GIVEN) };
    md_system_group_set_role(&b, MD_GROUP_ROLE_ENERGY, g, 2);
    md_system_group_add(&b, MD_STR_LIT("mine"), md_bitset_range(0, 1), MD_GROUP_ORIGIN_USER);
    const uint64_t groups = a.group.version;

    ASSERT_TRUE(md_system_merge(&a, &as, &b, &bs));
    ASSERT_EQ(3u, a.group.count);
    EXPECT_NE(groups, a.group.version);
    EXPECT_EQ(MD_GROUP_ORIGIN_USER, (md_group_origin_t)a.group.origin[2]);
    size_t count = 0;
    const uint32_t* order = md_system_group_role_order(&a.group, MD_GROUP_ROLE_ENERGY, &count);
    ASSERT_EQ(2u, count);
    EXPECT_TRUE(md_str_eq_cstr(a.group.name[order[1]], "second"));

    md_system_free(&a); md_system_state_free(&as); md_system_free(&b); md_system_state_free(&bs);
}

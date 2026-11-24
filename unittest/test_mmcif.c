#include "utest.h"
#include <math.h>

#include <md_mmcif.h>
#include <md_system.h>
#include <md_infer.h>
#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_arena_allocator.h>
#include <core/md_vec_math.h>

#define MAX_VALIDATION_SAMPLES 100

// The atom of an instance (by its label asym id), a residue of it (by its seq id) and a name, or -1
static int conn_atom(const md_system_t* sys, const char* asym, int seq, const char* name) {
    for (size_t i = 0; i < sys->instance.count; ++i) {
        if (!md_str_eq_cstr(MD_LABEL_STR(sys->instance.id[i]), asym)) continue;
        for (uint32_t r = sys->instance.res_offset[i]; r < sys->instance.res_offset[i + 1]; ++r) {
            if (sys->residue.seq_id[r] != seq) continue;
            for (uint32_t a = sys->residue.atom_offset[r]; a < sys->residue.atom_offset[r + 1]; ++a) {
                if (md_str_eq_cstr(md_atom_name(&sys->atom, a), name)) return (int)a;
            }
        }
    }
    return -1;
}

// The flags of the bond between two atoms, of the graph or an interaction; 0xFFFF for none
static uint32_t conn_bond(const md_system_t* sys, int a, int b) {
    if (a < 0 || b < 0) return 0xFFFF;
    md_bond_idx_t idx = md_system_bond_find(sys, a, b);
    if (idx < 0) idx = md_bond_find_interaction(&sys->bond, a, b);
    return idx < 0 ? 0xFFFF : (uint32_t)md_system_bond_flags(sys, idx);
}

// The kind and the origin of a bond, of its flags, as the file gives one
static uint32_t conn_given(md_bond_kind_t kind) {
    return (uint32_t)md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, kind), MD_BOND_ORIGIN_FILE);
}
#define CONN_FIELDS ((uint32_t)(MD_BOND_FLAG_KIND_MASK | MD_BOND_FLAG_ORIGIN_MASK | MD_BOND_FLAG_KIND_GIVEN))

UTEST(mmcif, 1fez) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/1fez.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_load_from_file(&sys, &sys_state, path, NULL) == MD_LOAD_OK;
    EXPECT_TRUE(result);

    if (result) {
        EXPECT_EQ(4097, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", md_str_ptr(md_atom_name(&sys.atom, 0)));
        //EXPECT_STREQ("LYS", mol.atom.resname[0].buf);
        //EXPECT_EQ(1, mol.atom.resid[0]);
        //EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(52.489, sys_state.xyz[0].x, 0.001);
        EXPECT_NEAR(21.292, sys_state.xyz[0].y, 0.001);
        EXPECT_NEAR(84.339, sys_state.xyz[0].z, 0.001);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

UTEST(mmcif, 2or2) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/2or2.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_load_from_file(&sys, &sys_state, path, NULL) == MD_LOAD_OK;
    EXPECT_TRUE(result);
    //md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_SYSTEM_INFER_ALL);

    if (result) {
        EXPECT_EQ(5382, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", md_str_ptr(md_atom_name(&sys.atom, 0)));
        //EXPECT_STREQ("ALA", sys.atom.resname[0].buf);
        //EXPECT_EQ(1, sys.atom.resid[0]);
        //EXPECT_STREQ("A", sys.atom.chainid[0].buf);
        EXPECT_NEAR(58.157, sys_state.xyz[0].x, 0.001);
        EXPECT_NEAR(49.822, sys_state.xyz[0].y, 0.001);
        EXPECT_NEAR(80.569, sys_state.xyz[0].z, 0.001);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

UTEST(mmcif, 8g7u) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/8g7u.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_load_from_file(&sys, &sys_state, path, NULL) == MD_LOAD_OK;
    EXPECT_TRUE(result);
    md_system_infer(&sys, &sys_state, NULL);

    if (result) {
        EXPECT_EQ(14229, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", md_str_ptr(md_atom_name(&sys.atom, 0)));
        //EXPECT_STREQ("PHE", sys.atom.resname[0].buf);
        //EXPECT_EQ(241, sys.atom.resid[0]);
        //EXPECT_STREQ("A", sys.atom.chainid[0].buf);
        EXPECT_NEAR(77.862,  sys_state.xyz[0].x, 0.001);
        EXPECT_NEAR(105.453, sys_state.xyz[0].y, 0.001);
        EXPECT_NEAR(80.951,  sys_state.xyz[0].z, 0.001);

        // The bonds of _struct_conn, of the kinds it says: two disulfides, two covalent links of the RNA, five zinc
        // coordinations and 46 hydrogen bonds of its base pairs
        EXPECT_EQ(55, (int)md_bond_list_count(&sys.given_bond));
        EXPECT_EQ(46, (int)sys.bond.interaction_count);
        EXPECT_EQ(conn_given(MD_BOND_KIND_COVALENT),   CONN_FIELDS & conn_bond(&sys, conn_atom(&sys, "B", 347, "SG"), conn_atom(&sys, "B", 402, "SG")));
        EXPECT_EQ(conn_given(MD_BOND_KIND_COVALENT),   CONN_FIELDS & conn_bond(&sys, conn_atom(&sys, "E", 1, "O3'"), conn_atom(&sys, "E", 2, "P")));
        EXPECT_EQ(conn_given(MD_BOND_KIND_COORDINATE), CONN_FIELDS & conn_bond(&sys, conn_atom(&sys, "A", 813, "SG"), conn_atom(&sys, "G", 1001, "ZN")));
        EXPECT_EQ(conn_given(MD_BOND_KIND_HYDROGEN),   CONN_FIELDS & conn_bond(&sys, conn_atom(&sys, "E", 1, "O6"), conn_atom(&sys, "F", 24, "N4")));
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
}


// STRUCT_CONN. A lysine and an aspartate of chain A, whose OD2 has two alternate locations, a ligand of a C=O and a zinc.
// The connections: a salt bridge, the same to a symmetry mate, a double bond of the ligand, a coordination of the second
// alternate location, a hydrogen bond of a mismatched pair, and a connection of a type of no bond.
static const char conn_cif[] =
"data_CONN\n"
"loop_\n"
"_struct_conn.id\n"
"_struct_conn.conn_type_id\n"
"_struct_conn.ptnr1_label_asym_id\n"
"_struct_conn.ptnr1_label_comp_id\n"
"_struct_conn.ptnr1_label_seq_id\n"
"_struct_conn.ptnr1_label_atom_id\n"
"_struct_conn.pdbx_ptnr1_label_alt_id\n"
"_struct_conn.ptnr1_auth_seq_id\n"
"_struct_conn.ptnr1_symmetry\n"
"_struct_conn.ptnr2_label_asym_id\n"
"_struct_conn.ptnr2_label_comp_id\n"
"_struct_conn.ptnr2_label_seq_id\n"
"_struct_conn.ptnr2_label_atom_id\n"
"_struct_conn.pdbx_ptnr2_label_alt_id\n"
"_struct_conn.ptnr2_auth_seq_id\n"
"_struct_conn.ptnr2_symmetry\n"
"_struct_conn.pdbx_value_order\n"
"saltbr1 saltbr A LYS 1 NZ ? 1 1_555 A ASP 2 OD1 ? 2 1_555 ?\n"
"saltbr2 saltbr A LYS 1 NZ ? 1 1_555 A ASP 2 OD2 ? 2 2_555 ?\n"
"covale1 covale B LIG . C1 ? 101 1_555 B LIG . O1 ? 101 1_555 doub\n"
"metalc1 metalc C ZN . ZN ? 201 1_555 A ASP 2 OD2 B 2 1_555 ?\n"
"mismat1 MISMAT A LYS 1 N ? 1 1_555 A ASP 2 O ? 2 1_555 ?\n"
"other1 unknown A LYS 1 N ? 1 1_555 A ASP 2 OD1 ? 2 1_555 ?\n"
"#\n"
"loop_\n"
"_atom_site.group_PDB\n"
"_atom_site.id\n"
"_atom_site.type_symbol\n"
"_atom_site.label_atom_id\n"
"_atom_site.label_alt_id\n"
"_atom_site.label_comp_id\n"
"_atom_site.label_asym_id\n"
"_atom_site.label_entity_id\n"
"_atom_site.label_seq_id\n"
"_atom_site.Cartn_x\n"
"_atom_site.Cartn_y\n"
"_atom_site.Cartn_z\n"
"_atom_site.auth_seq_id\n"
"ATOM 1 N N . LYS A 1 1 0.000 0.000 0.000 1\n"
"ATOM 2 N NZ . LYS A 1 1 5.000 0.000 0.000 1\n"
"ATOM 3 O O . ASP A 1 2 0.000 2.900 0.000 2\n"
"ATOM 4 O OD1 . ASP A 1 2 7.900 0.000 0.000 2\n"
"ATOM 5 O OD2 A ASP A 1 2 5.000 2.900 0.000 2\n"
"ATOM 6 O OD2 B ASP A 1 2 5.000 -2.900 0.000 2\n"
"HETATM 7 C C1 . LIG B 2 . 20.000 0.000 0.000 101\n"
"HETATM 8 O O1 . LIG B 2 . 21.220 0.000 0.000 101\n"
"HETATM 9 ZN ZN . ZN C 3 . 30.000 0.000 0.000 201\n"
"#\n";

UTEST(mmcif, struct_conn) {
    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t state = { .alloc = md_get_heap_allocator() };
    ASSERT_TRUE(md_mmcif_system_load_from_str(&sys, &state, md_str_from_cstr(conn_cif), NULL) == MD_LOAD_OK);
    ASSERT_EQ(8, (int)sys.atom.count);     // The first alternate location of OD2
    // The salt bridge, the double bond and the hydrogen bond; not to a symmetry mate, nor to an atom not loaded, nor of
    // a type of no bond
    EXPECT_EQ(3, (int)md_bond_list_count(&sys.given_bond));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));

    const int nz = conn_atom(&sys, "A", 1, "NZ"), od1 = conn_atom(&sys, "A", 2, "OD1"), od2 = conn_atom(&sys, "A", 2, "OD2");
    EXPECT_EQ(conn_given(MD_BOND_KIND_IONIC), CONN_FIELDS & conn_bond(&sys, nz, od1));
    EXPECT_EQ(0xFFFFu, conn_bond(&sys, nz, od2));
    EXPECT_EQ(0xFFFFu, conn_bond(&sys, conn_atom(&sys, "C", 201, "ZN"), od2));
    EXPECT_EQ(conn_given(MD_BOND_KIND_HYDROGEN), CONN_FIELDS & conn_bond(&sys, conn_atom(&sys, "A", 1, "N"), conn_atom(&sys, "A", 2, "O")));
    EXPECT_EQ(2, (int)sys.bond.interaction_count);
    // The ligand's bond as the file gives it, of its order: inference finds it too, and keeps the file's
    const uint32_t c_o = conn_bond(&sys, conn_atom(&sys, "B", 101, "C1"), conn_atom(&sys, "B", 101, "O1"));
    EXPECT_EQ(conn_given(MD_BOND_KIND_COVALENT), CONN_FIELDS & c_o);
    EXPECT_EQ(2, md_bond_order((md_bond_flags_t)c_o));

    md_system_free(&sys);
    md_system_state_free(&state);
}

// ASSEMBLIES. A chain A of two alanines, a chain B of one and two waters in asym C, and assemblies which take them by
// rows of their own, by a product of operators and by ranges.
static const char assembly_cif[] =
"data_TEST\n"
"loop_\n"
"_entity.id\n"
"_entity.type\n"
"_entity.pdbx_description\n"
"1 polymer 'Peptide'\n"
"2 water 'water'\n"
"#\n"
"_entity_poly.entity_id 1\n"
"_entity_poly.type 'polypeptide(L)'\n"
"#\n"
"loop_\n"
"_pdbx_struct_assembly_gen.assembly_id\n"
"_pdbx_struct_assembly_gen.oper_expression\n"
"_pdbx_struct_assembly_gen.asym_id_list\n"
"1 1 A,C\n"
"1 2 B\n"
"2 '(1,2)(3)' A\n"
"3 1 A\n"
"3 1 C\n"
"4 '1-2' B\n"
"5 1 Q\n"
"#\n"
"loop_\n"
"_pdbx_struct_oper_list.id\n"
"_pdbx_struct_oper_list.type\n"
"_pdbx_struct_oper_list.matrix[1][1]\n"
"_pdbx_struct_oper_list.matrix[1][2]\n"
"_pdbx_struct_oper_list.matrix[1][3]\n"
"_pdbx_struct_oper_list.vector[1]\n"
"_pdbx_struct_oper_list.matrix[2][1]\n"
"_pdbx_struct_oper_list.matrix[2][2]\n"
"_pdbx_struct_oper_list.matrix[2][3]\n"
"_pdbx_struct_oper_list.vector[2]\n"
"_pdbx_struct_oper_list.matrix[3][1]\n"
"_pdbx_struct_oper_list.matrix[3][2]\n"
"_pdbx_struct_oper_list.matrix[3][3]\n"
"_pdbx_struct_oper_list.vector[3]\n"
"1 'identity operation' 1 0 0 0 0 1 0 0 0 0 1 0\n"
"2 'point symmetry operation' -1 0 0 10 0 -1 0 0 0 0 1 0\n"
"3 'point symmetry operation' 0 -1 0 0 1 0 0 0 0 0 1 5\n"
"#\n"
"loop_\n"
"_atom_site.group_PDB\n"
"_atom_site.id\n"
"_atom_site.type_symbol\n"
"_atom_site.label_atom_id\n"
"_atom_site.label_alt_id\n"
"_atom_site.label_comp_id\n"
"_atom_site.label_asym_id\n"
"_atom_site.label_entity_id\n"
"_atom_site.label_seq_id\n"
"_atom_site.Cartn_x\n"
"_atom_site.Cartn_y\n"
"_atom_site.Cartn_z\n"
"_atom_site.auth_asym_id\n"
"_atom_site.auth_seq_id\n"
"ATOM 1 N N . ALA A 1 1 0.000 0.000 0.000 A 1\n"
"ATOM 2 C CA . ALA A 1 1 1.458 0.000 0.000 A 1\n"
"ATOM 3 C C . ALA A 1 1 2.009 1.420 0.000 A 1\n"
"ATOM 4 O O . ALA A 1 1 1.251 2.390 0.000 A 1\n"
"ATOM 5 N N . ALA A 1 2 3.332 1.536 0.000 A 2\n"
"ATOM 6 C CA . ALA A 1 2 3.988 2.839 0.000 A 2\n"
"ATOM 7 C C . ALA A 1 2 5.504 2.693 0.000 A 2\n"
"ATOM 8 O O . ALA A 1 2 6.030 1.580 0.000 A 2\n"
"ATOM 9 N N . ALA B 1 1 0.000 10.000 0.000 B 1\n"
"ATOM 10 C CA . ALA B 1 1 1.458 10.000 0.000 B 1\n"
"ATOM 11 C C . ALA B 1 1 2.009 11.420 0.000 B 1\n"
"ATOM 12 O O . ALA B 1 1 1.251 12.390 0.000 B 1\n"
"HETATM 13 O O . HOH C 2 . 0.000 0.000 6.000 A 101\n"
"HETATM 14 O O . HOH C 2 . 0.000 0.000 9.000 A 102\n"
"#\n";

static bool part_is(const md_assembly_data_t* a, uint32_t part, const md_instance_idx_t* inst, size_t inst_count, const uint32_t* oper, size_t oper_count) {
    const md_urange_t ir = md_assembly_group_inst_range(a, part);
    const md_urange_t orr = md_assembly_group_oper_range(a, part);
    if (ir.end - ir.beg != inst_count || orr.end - orr.beg != oper_count) return false;
    for (size_t i = 0; i < inst_count; ++i) if (a->inst[ir.beg + i] != inst[i]) return false;
    for (size_t i = 0; i < oper_count; ++i) if (a->oper[orr.beg + i] != oper[i]) return false;
    return true;
}

UTEST(mmcif, assemblies) {
    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t state = { .alloc = md_get_heap_allocator() };
    ASSERT_TRUE(md_mmcif_system_load_from_str(&sys, &state, (md_str_t){assembly_cif, sizeof(assembly_cif) - 1}, NULL) == MD_LOAD_OK);
    ASSERT_EQ(14, sys.atom.count);
    // A, B, and the waters of C each an instance of their own
    ASSERT_EQ(4, sys.instance.count);

    const md_assembly_data_t* a = &sys.assembly;
    // 5 names an asym which has no atoms: it has nothing to show and is left out
    ASSERT_EQ(4, a->count);
    EXPECT_STREQ("1", md_str_ptr(md_assembly_id(a, 0)));
    EXPECT_STREQ("4", md_str_ptr(md_assembly_id(a, 3)));
    // The three operators, and the two of the product (1)(3), (2)(3)
    ASSERT_EQ(5, a->oper_count);
    EXPECT_STREQ("1x3", md_str_ptr(MD_LABEL_STR(a->oper_id[3])));
    EXPECT_STREQ("2x3", md_str_ptr(MD_LABEL_STR(a->oper_id[4])));

    // 1: A and the waters by the identity, B by the other
    {
        const md_urange_t parts = md_assembly_group_range(a, 0);
        ASSERT_EQ(2, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {0, 2, 3};
        const uint32_t o0[] = {0};
        const md_instance_idx_t i1[] = {1};
        const uint32_t o1[] = {1};
        EXPECT_TRUE(part_is(a, parts.beg + 0, i0, 3, o0, 1));
        EXPECT_TRUE(part_is(a, parts.beg + 1, i1, 1, o1, 1));
        EXPECT_EQ(2, md_assembly_copy_count(a, 0));
    }
    // 2: A by the products, the factor of the last group applied first: 2 (half a turn about z and 10 along x) after 3
    // (a quarter turn about z and 5 along z), the matrices of the file by rows
    {
        const md_urange_t parts = md_assembly_group_range(a, 1);
        ASSERT_EQ(1, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {0};
        const uint32_t o0[] = {3, 4};
        EXPECT_TRUE(part_is(a, parts.beg, i0, 1, o0, 2));
        const md_float4_t p = {1, 2, 3, 1};
        const md_float4_t q = md_float4x4_mul_float4(a->oper_transform[4], p);
        EXPECT_NEAR(12.0f, q.x, 1e-5f);
        EXPECT_NEAR(-1.0f, q.y, 1e-5f);
        EXPECT_NEAR(8.0f, q.z, 1e-5f);
    }
    // 3: two rows of the same operator are one part
    {
        const md_urange_t parts = md_assembly_group_range(a, 2);
        ASSERT_EQ(1, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {0, 2, 3};
        const uint32_t o0[] = {0};
        EXPECT_TRUE(part_is(a, parts.beg, i0, 3, o0, 1));
    }
    // 4: a range of operators
    {
        const md_urange_t parts = md_assembly_group_range(a, 3);
        ASSERT_EQ(1, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {1};
        const uint32_t o0[] = {0, 1};
        EXPECT_TRUE(part_is(a, parts.beg, i0, 1, o0, 2));
    }

    md_system_free(&sys);
    md_system_state_free(&state);
}

// The assembly of a file in the key-value form of mmCIF (a single row of each category): all asyms by the identity
UTEST(mmcif, assembly_key_value) {
    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t state = { .alloc = md_get_heap_allocator() };
    ASSERT_TRUE(md_mmcif_system_load_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"), NULL) == MD_LOAD_OK);
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    ASSERT_EQ(1, a->group_count);
    ASSERT_EQ(1, a->oper_count);
    const md_urange_t ir = md_assembly_group_inst_range(a, 0);
    EXPECT_EQ(sys.instance.count, ir.end - ir.beg);
    const md_float4x4_t I = md_float4x4_ident();
    EXPECT_EQ(0, md_memcmp(&I, &a->oper_transform[0], sizeof(md_float4x4_t)));
    md_system_free(&sys);
    md_system_state_free(&state);
}

#include <md_mmcif.c>

UTEST(mmcif, tokenizer) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    {
        md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/1fez.cif");
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
        md_file_reader_t reader;
        ASSERT_TRUE(md_file_reader_init(&reader, file, MD_MEGABYTES(1), alloc));

        mmcif_parse_state_t state = {
            .reader = &reader.r,
            .sb = md_str_writer_create(alloc),
        };

        md_str_t tok;
        
        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("data_XXXX")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("loop_")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_author.name")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_author.pdbx_ordinal")));

        ASSERT_TRUE(mmcif_peek_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("Morais, M.C.")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("Morais, M.C.")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("1")));

        for (size_t i = 0; i < 1000; ++i) {
            ASSERT_TRUE(mmcif_next_token(&tok, &state));
        }

        md_file_close(&file);
    }

    {
        md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/8g7u.cif");
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
        md_file_reader_t reader;
        ASSERT_TRUE(md_file_reader_init(&reader, file, MD_MEGABYTES(1), alloc));

        mmcif_parse_state_t state = {
            .reader = &reader.r,
            .sb = md_str_writer_create(alloc),
        };

        md_str_t tok;

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("data_8G7U")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_entry.id")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("8G7U")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_conform.dict_name")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("mmcif_pdbx.dic")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_conform.dict_version")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("5.381")));

        for (size_t i = 0; i < 1000; ++i) {
            ASSERT_TRUE(mmcif_next_token(&tok, &state));
        }
        
        md_file_close(&file);
    }


    // TODO Test parsing of the above section

    md_temp_end(temp);
}

UTEST(mmcif, parse_section) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_str_t section = MD_STR_INIT(
        "_entity_poly.entity_id                      1 \n"
        "_entity_poly.type                           \"polypeptide(L)\" \n"
        "_entity_poly.nstd_linkage                   no \n"
        "_entity_poly.nstd_monomer                   no \n"
        "_entity_poly.pdbx_seq_one_letter_code       \n"
        ";KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
        "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
        "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
        "TIETMQELESVMEHIE\n"
        ";\n"
        "_entity_poly.pdbx_seq_one_letter_code_can   \n"
        ";KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
        "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
        "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
        "TIETMQELESVMEHIE\n"
        ";\n"
        "_entity_poly.pdbx_strand_id                 A,B \n"
        "_entity_poly.pdbx_target_identifier         ? \n"
        "# \n"
    );

    {
        md_reader_i reader = md_reader_from_str(section);
        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, alloc));
        EXPECT_EQ(8, sec.num_fields);
        EXPECT_EQ(1, sec.num_rows);

        EXPECT_STREQ("_entity_poly.entity_id", md_str_ptr(sec.headers[0]));
        EXPECT_STREQ("1", md_str_ptr(sec.values[0]));

        EXPECT_STREQ("_entity_poly.type", md_str_ptr(sec.headers[1]));
        EXPECT_STREQ("polypeptide(L)", md_str_ptr(sec.values[1]));

        EXPECT_STREQ("_entity_poly.nstd_linkage", md_str_ptr(sec.headers[2]));
        EXPECT_STREQ("no", md_str_ptr(sec.values[2]));

        EXPECT_STREQ("_entity_poly.nstd_monomer", md_str_ptr(sec.headers[3]));
        EXPECT_STREQ("no", md_str_ptr(sec.values[3]));

        EXPECT_STREQ("_entity_poly.pdbx_seq_one_letter_code", md_str_ptr(sec.headers[4]));
        EXPECT_STREQ(
            "KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
            "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
            "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
            "TIETMQELESVMEHIE",
            md_str_ptr(sec.values[4]));

        EXPECT_STREQ("_entity_poly.pdbx_seq_one_letter_code_can", md_str_ptr(sec.headers[5]));
        EXPECT_STREQ(
            "KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
            "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
            "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
            "TIETMQELESVMEHIE",
            md_str_ptr(sec.values[5]));

        EXPECT_STREQ("_entity_poly.pdbx_strand_id", md_str_ptr(sec.headers[6]));
        EXPECT_STREQ("A,B", md_str_ptr(sec.values[6]));

        EXPECT_STREQ("_entity_poly.pdbx_target_identifier", md_str_ptr(sec.headers[7]));
        EXPECT_STREQ("?", md_str_ptr(sec.values[7]));
    }

    md_str_t item_looped = MD_STR_INIT(
        "_entity.id \n"
        "_entity.type \n"
        "_entity.src_method \n"
        "_entity.pdbx_description \n"
        "_entity.formula_weight \n"
        "_entity.pdbx_number_of_molecules \n"
        "_entity.pdbx_ec \n"
        "_entity.pdbx_mutation \n"
        "_entity.pdbx_fragment \n"
        "_entity.details \n"
        "1 polymer     man 'Antiviral innate immune response receptor RIG-I' 106740.555 2 3.6.4.13 ? ? ? \n"
        "2 polymer     man 'E3 ubiquitin-protein ligase RNF135'              47946.305  2 2.3.2.27 ? ? ? \n"
        "3 polymer     man p3dsRNA24a                                        7885.581   1 ?        ? ? ? \n"
        "4 polymer     man p3dsRNA24b                                        7788.551   1 ?        ? ? ? \n"
        "5 non-polymer syn 'ZINC ION'                                        65.409     2 ?        ? ? ? \n"
        "# \n"
    );

    {
        md_reader_i reader = md_reader_from_str(item_looped);
        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
            .in_loop = true,
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, alloc));

        EXPECT_EQ(10, sec.num_fields);
        EXPECT_EQ(5, sec.num_rows);

        EXPECT_STREQ("_entity.id", md_str_ptr(sec.headers[0]));
        EXPECT_STREQ("1", md_str_ptr(mmcif_section_value(&sec, 0, 0)));
        EXPECT_STREQ("2", md_str_ptr(mmcif_section_value(&sec, 1, 0)));
        EXPECT_STREQ("3", md_str_ptr(mmcif_section_value(&sec, 2, 0)));
        EXPECT_STREQ("4", md_str_ptr(mmcif_section_value(&sec, 3, 0)));
        EXPECT_STREQ("5", md_str_ptr(mmcif_section_value(&sec, 4, 0)));

        EXPECT_STREQ("_entity.type", md_str_ptr(sec.headers[1]));
        EXPECT_STREQ("polymer", md_str_ptr(mmcif_section_value(&sec, 0, 1)));
        EXPECT_STREQ("polymer", md_str_ptr(mmcif_section_value(&sec, 1, 1)));
        EXPECT_STREQ("polymer", md_str_ptr(mmcif_section_value(&sec, 2, 1)));
        EXPECT_STREQ("polymer", md_str_ptr(mmcif_section_value(&sec, 3, 1)));
        EXPECT_STREQ("non-polymer", md_str_ptr(mmcif_section_value(&sec, 4, 1)));

        EXPECT_STREQ("_entity.src_method", md_str_ptr(sec.headers[2]));
        EXPECT_STREQ("man", md_str_ptr(mmcif_section_value(&sec, 0, 2)));
        EXPECT_STREQ("man", md_str_ptr(mmcif_section_value(&sec, 1, 2)));
        EXPECT_STREQ("man", md_str_ptr(mmcif_section_value(&sec, 2, 2)));
        EXPECT_STREQ("man", md_str_ptr(mmcif_section_value(&sec, 3, 2)));
        EXPECT_STREQ("syn", md_str_ptr(mmcif_section_value(&sec, 4, 2)));

        EXPECT_STREQ("_entity.pdbx_description", md_str_ptr(sec.headers[3]));
        EXPECT_STREQ("Antiviral innate immune response receptor RIG-I", md_str_ptr(mmcif_section_value(&sec, 0, 3)));
        EXPECT_STREQ("E3 ubiquitin-protein ligase RNF135", md_str_ptr(mmcif_section_value(&sec, 1, 3)));
        EXPECT_STREQ("p3dsRNA24a", md_str_ptr(mmcif_section_value(&sec, 2, 3)));
        EXPECT_STREQ("p3dsRNA24b", md_str_ptr(mmcif_section_value(&sec, 3, 3)));
        EXPECT_STREQ("ZINC ION", md_str_ptr(mmcif_section_value(&sec, 4, 3)));

        EXPECT_STREQ("_entity.formula_weight", md_str_ptr(sec.headers[4]));
        EXPECT_STREQ("106740.555", md_str_ptr(mmcif_section_value(&sec, 0, 4)));
        EXPECT_STREQ("47946.305", md_str_ptr(mmcif_section_value(&sec, 1, 4)));
        EXPECT_STREQ("7885.581", md_str_ptr(mmcif_section_value(&sec, 2, 4)));
        EXPECT_STREQ("7788.551", md_str_ptr(mmcif_section_value(&sec, 3, 4)));
        EXPECT_STREQ("65.409", md_str_ptr(mmcif_section_value(&sec, 4, 4)));

        EXPECT_STREQ("_entity.pdbx_number_of_molecules", md_str_ptr(sec.headers[5]));
        EXPECT_STREQ("2", md_str_ptr(mmcif_section_value(&sec, 0, 5)));
        EXPECT_STREQ("2", md_str_ptr(mmcif_section_value(&sec, 1, 5)));
        EXPECT_STREQ("1", md_str_ptr(mmcif_section_value(&sec, 2, 5)));
        EXPECT_STREQ("1", md_str_ptr(mmcif_section_value(&sec, 3, 5)));
        EXPECT_STREQ("2", md_str_ptr(mmcif_section_value(&sec, 4, 5)));
    }

    md_temp_end(temp);
}

UTEST(mmcif, parse_2or2_comprehensive) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/2or2.cif");
    
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    bool result = md_mmcif_system_load_from_file(&sys, &sys_state, path, NULL) == MD_LOAD_OK;
    ASSERT_TRUE(result);
    
    // Check basic structure properties
    EXPECT_GT(sys.atom.count, 0);
    EXPECT_GT(sys.residue.count, 0);
    EXPECT_GT(sys.instance.count, 0);
    
    // Check that coordinates are reasonable (not all zeros or infinities)
    bool has_nonzero_coord = false;
    for (int64_t i = 0; i < sys.atom.count && i < MAX_VALIDATION_SAMPLES; ++i) {
        EXPECT_FALSE(isnan(sys_state.xyz[i].x));
        EXPECT_FALSE(isnan(sys_state.xyz[i].y));
        EXPECT_FALSE(isnan(sys_state.xyz[i].z));
        EXPECT_FALSE(isinf(sys_state.xyz[i].x));
        EXPECT_FALSE(isinf(sys_state.xyz[i].y));
        EXPECT_FALSE(isinf(sys_state.xyz[i].z));
        
        if (sys_state.xyz[i].x != 0.0f || sys_state.xyz[i].y != 0.0f || sys_state.xyz[i].z != 0.0f) {
            has_nonzero_coord = true;
        }
    }
    EXPECT_TRUE(has_nonzero_coord);
    
    md_temp_end(temp);
}

UTEST(mmcif, nonexistent_file) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/nonexistent.cif");
    
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_FALSE(md_mmcif_system_load_from_file(&sys, &sys_state, path, NULL) == MD_LOAD_OK);
    
    // Should be safe to free even when init failed
    md_system_free(&sys);
    md_system_state_free(&sys_state);
    md_temp_end(temp);
}

UTEST(mmcif, advance_to_next_control) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/1fez.cif");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    md_file_reader_t reader;
    ASSERT_TRUE(md_file_reader_init(&reader, file, MD_MEGABYTES(1), alloc));

    mmcif_parse_state_t state = {
        .reader = &reader.r,
        .sb = md_str_writer_create(alloc),
    };

    md_str_t tok;

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("data_XXXX")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("loop_")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_author.name")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_audit_author.pdbx_ordinal")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(md_str_eq(tok, MD_STR_LIT("_struct.entry_id")));

    md_temp_end(temp);
}

UTEST(mmcif, atom_site_many_columns) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    // More _atom_site columns than a fixed row buffer of 64 would hold, with the ones the reader uses after them
    const int num_extra = 100;
    md_str_writer_t sb = md_str_writer_create(alloc);
    md_str_writer_push_cstr(&sb, "data_test\nloop_\n");
    for (int i = 0; i < num_extra; ++i) {
        md_str_writer_fmt(&sb, "_atom_site.extra_%d\n", i);
    }
    md_str_writer_push_cstr(&sb,
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_asym_id\n");
    const char* rows[] = {
        "ATOM 1 N N  . ALA A 1 1 ? 1.000 2.000 3.000 1 A\n",
        "ATOM 2 C CA . ALA A 1 1 ? 4.000 5.000 6.000 1 A\n",
        "ATOM 3 C C  . ALA A 1 1 ? 7.000 8.000 9.000 1 A\n",
    };
    for (size_t r = 0; r < MD_COUNTOF(rows); ++r) {
        for (int i = 0; i < num_extra; ++i) {
            md_str_writer_push_cstr(&sb, "x ");
        }
        md_str_writer_push_cstr(&sb, rows[r]);
    }
    md_str_writer_push_cstr(&sb, "#\n");

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_load_from_str(&sys, &sys_state, md_str_writer_to_str(&sb), NULL) == MD_LOAD_OK);
    ASSERT_EQ(3, sys.atom.count);
    EXPECT_STREQ("CA", md_str_ptr(md_atom_name(&sys.atom, 1)));
    EXPECT_NEAR(7.0f, sys_state.xyz[2].x, 1.0e-5f);
    EXPECT_NEAR(8.0f, sys_state.xyz[2].y, 1.0e-5f);
    EXPECT_NEAR(9.0f, sys_state.xyz[2].z, 1.0e-5f);

    md_temp_end(temp);
}

// 9gfu misses residues within both of its chains (_pdbx_unobs_or_zero_occ_residues): label_seq_id 40-41, 67-72 and
// 259-260 of chain A, and 62-71 of chain B. Each leaves a break between the residues on either side, whose numbers (the
// residues of the sequence) count the missing ones exactly. The residues missing at the ends of the chains leave none.
UTEST(mmcif, chain_breaks) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_mmcif_system_load_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/9gfu.cif"), NULL) == MD_LOAD_OK);
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));

    const struct {
        const char* chain;
        int before;
        int after;
    } expected[] = {
        {"A", 39, 42}, {"A", 66, 73}, {"A", 258, 261}, {"B", 61, 72},
    };
    size_t found = 0;
    for (size_t i = 0; i + 1 < sys.residue.count; ++i) {
        if (!(md_system_residue_flags(&sys, i) & MD_RESIDUE_FLAG_BREAK_AFTER)) continue;
        if (found < MD_COUNTOF(expected)) {
            const md_instance_idx_t inst = md_system_instance_find_by_residue_idx(&sys, i);
            EXPECT_TRUE(md_str_eq_cstr(md_system_instance_id(&sys, (size_t)inst), expected[found].chain));
            EXPECT_EQ(expected[found].before, md_system_residue_seq_id(&sys, i));
            EXPECT_EQ(expected[found].after,  md_system_residue_seq_id(&sys, i + 1));
            EXPECT_TRUE(md_system_residue_flags(&sys, i + 1) & MD_RESIDUE_FLAG_BREAK_BEFORE);
            EXPECT_EQ(inst, md_system_instance_find_by_residue_idx(&sys, i + 1));
            md_atom_idx_t atoms[2];
            EXPECT_TRUE(md_system_chain_break_atoms(atoms, &sys, i));
        }
        ++found;
    }
    EXPECT_EQ(MD_COUNTOF(expected), found);

    md_vm_arena_destroy(alloc);
}

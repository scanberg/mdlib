#include "utest.h"
#include <math.h>

#include <md_mmcif.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#define MAX_VALIDATION_SAMPLES 100

UTEST(mmcif, 1fez) {
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/1fez.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_init_from_file(&sys, &sys_state, path);
    EXPECT_TRUE(result);

    if (result) {
        EXPECT_EQ(4097, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&sys.atom, 0)));
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
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/2or2.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_init_from_file(&sys, &sys_state, path);
    EXPECT_TRUE(result);
    //md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_UTIL_INFER_ALL);

    if (result) {
        EXPECT_EQ(5382, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&sys.atom, 0)));
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
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/8g7u.cif");

    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t sys_state = { .alloc = md_get_heap_allocator() };
    bool result = md_mmcif_system_init_from_file(&sys, &sys_state, path);
    EXPECT_TRUE(result);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    if (result) {
        EXPECT_EQ(14229, sys.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&sys.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&sys.atom, 0)));
        //EXPECT_STREQ("PHE", sys.atom.resname[0].buf);
        //EXPECT_EQ(241, sys.atom.resid[0]);
        //EXPECT_STREQ("A", sys.atom.chainid[0].buf);
        EXPECT_NEAR(77.862,  sys_state.xyz[0].x, 0.001);
        EXPECT_NEAR(105.453, sys_state.xyz[0].y, 0.001);
        EXPECT_NEAR(80.951,  sys_state.xyz[0].z, 0.001);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
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
    const md_urange_t ir = md_assembly_part_inst_range(a, part);
    const md_urange_t orr = md_assembly_part_oper_range(a, part);
    if (ir.end - ir.beg != inst_count || orr.end - orr.beg != oper_count) return false;
    for (size_t i = 0; i < inst_count; ++i) if (a->inst[ir.beg + i] != inst[i]) return false;
    for (size_t i = 0; i < oper_count; ++i) if (a->oper[orr.beg + i] != oper[i]) return false;
    return true;
}

UTEST(mmcif, assemblies) {
    md_system_t sys = { .alloc = md_get_heap_allocator() };
    md_system_state_t state = { .alloc = md_get_heap_allocator() };
    ASSERT_TRUE(md_mmcif_system_init_from_str(&sys, &state, (str_t){assembly_cif, sizeof(assembly_cif) - 1}));
    ASSERT_EQ(14, sys.atom.count);
    // A, B, and the waters of C each an instance of their own
    ASSERT_EQ(4, sys.instance.count);

    const md_assembly_data_t* a = &sys.assembly;
    // 5 names an asym which has no atoms: it has nothing to show and is left out
    ASSERT_EQ(4, a->count);
    EXPECT_STREQ("1", str_ptr(md_assembly_id(a, 0)));
    EXPECT_STREQ("4", str_ptr(md_assembly_id(a, 3)));
    // The three operators, and the two of the product (1)(3), (2)(3)
    ASSERT_EQ(5, a->oper_count);
    EXPECT_STREQ("1x3", str_ptr(LBL_TO_STR(a->oper_id[3])));
    EXPECT_STREQ("2x3", str_ptr(LBL_TO_STR(a->oper_id[4])));

    // 1: A and the waters by the identity, B by the other
    {
        const md_urange_t parts = md_assembly_part_range(a, 0);
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
        const md_urange_t parts = md_assembly_part_range(a, 1);
        ASSERT_EQ(1, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {0};
        const uint32_t o0[] = {3, 4};
        EXPECT_TRUE(part_is(a, parts.beg, i0, 1, o0, 2));
        const vec4_t p = {1, 2, 3, 1};
        const vec4_t q = mat4_mul_vec4(a->oper_transform[4], p);
        EXPECT_NEAR(12.0f, q.x, 1e-5f);
        EXPECT_NEAR(-1.0f, q.y, 1e-5f);
        EXPECT_NEAR(8.0f, q.z, 1e-5f);
    }
    // 3: two rows of the same operator are one part
    {
        const md_urange_t parts = md_assembly_part_range(a, 2);
        ASSERT_EQ(1, parts.end - parts.beg);
        const md_instance_idx_t i0[] = {0, 2, 3};
        const uint32_t o0[] = {0};
        EXPECT_TRUE(part_is(a, parts.beg, i0, 3, o0, 1));
    }
    // 4: a range of operators
    {
        const md_urange_t parts = md_assembly_part_range(a, 3);
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
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
    const md_assembly_data_t* a = &sys.assembly;
    ASSERT_EQ(1, a->count);
    ASSERT_EQ(1, a->part_count);
    ASSERT_EQ(1, a->oper_count);
    const md_urange_t ir = md_assembly_part_inst_range(a, 0);
    EXPECT_EQ(sys.instance.count, ir.end - ir.beg);
    const mat4_t I = mat4_ident();
    EXPECT_EQ(0, MEMCMP(&I, &a->oper_transform[0], sizeof(mat4_t)));
    md_system_free(&sys);
    md_system_state_free(&state);
}

#include <md_mmcif.c>

UTEST(mmcif, tokenizer) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    {
        str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1fez.cif");
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
        char* buf = md_temp_alloc(temp, MEGABYTES(1));
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
        };

        str_t tok;
        
        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("data_XXXX")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("loop_")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.name")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.pdbx_ordinal")));

        ASSERT_TRUE(mmcif_peek_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("Morais, M.C.")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("Morais, M.C.")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("1")));

        for (size_t i = 0; i < 1000; ++i) {
            ASSERT_TRUE(mmcif_next_token(&tok, &state));
        }

        md_file_close(&file);
    }

    {
        str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/8g7u.cif");
        md_file_t file = {0};
        ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
        char* buf = md_temp_alloc(temp, MEGABYTES(1));
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
        };

        str_t tok;

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("data_8G7U")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("_entry.id")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("8G7U")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_conform.dict_name")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("mmcif_pdbx.dic")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_conform.dict_version")));

        ASSERT_TRUE(mmcif_next_token(&tok, &state));
        EXPECT_TRUE(str_eq(tok, STR_LIT("5.381")));

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

    str_t section = STR_INIT(
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
        md_buffered_reader_t reader = md_buffered_reader_from_str(section);
        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, alloc));
        EXPECT_EQ(8, sec.num_fields);
        EXPECT_EQ(1, sec.num_rows);

        EXPECT_STREQ("_entity_poly.entity_id", str_ptr(sec.headers[0]));
        EXPECT_STREQ("1", str_ptr(sec.values[0]));

        EXPECT_STREQ("_entity_poly.type", str_ptr(sec.headers[1]));
        EXPECT_STREQ("polypeptide(L)", str_ptr(sec.values[1]));

        EXPECT_STREQ("_entity_poly.nstd_linkage", str_ptr(sec.headers[2]));
        EXPECT_STREQ("no", str_ptr(sec.values[2]));

        EXPECT_STREQ("_entity_poly.nstd_monomer", str_ptr(sec.headers[3]));
        EXPECT_STREQ("no", str_ptr(sec.values[3]));

        EXPECT_STREQ("_entity_poly.pdbx_seq_one_letter_code", str_ptr(sec.headers[4]));
        EXPECT_STREQ(
            "KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
            "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
            "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
            "TIETMQELESVMEHIE",
            str_ptr(sec.values[4]));

        EXPECT_STREQ("_entity_poly.pdbx_seq_one_letter_code_can", str_ptr(sec.headers[5]));
        EXPECT_STREQ(
            "KIEAVIFDWAGTTVDYGCFAPLEVFMEIFHKRGVAITAEEARKPMGLLKIDHVRALTEMPRIASEWNRVFRQLPTEADIQ\n"
            "EMYEEFEEILFAILPRYASPINAVKEVIASLRERGIKIGSTTGYTREMMDIVAKEAALQGYKPDFLVTPDDVPAGRPYPW\n"
            "MCYKNAMELGVYPMNHMIKVGDTVSDMKEGRNAGMWTVGVILGSSELGLTEEEVENMDSVELREKIEVVRNRFVENGAHF\n"
            "TIETMQELESVMEHIE",
            str_ptr(sec.values[5]));

        EXPECT_STREQ("_entity_poly.pdbx_strand_id", str_ptr(sec.headers[6]));
        EXPECT_STREQ("A,B", str_ptr(sec.values[6]));

        EXPECT_STREQ("_entity_poly.pdbx_target_identifier", str_ptr(sec.headers[7]));
        EXPECT_STREQ("?", str_ptr(sec.values[7]));
    }

    str_t item_looped = STR_INIT(
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
        md_buffered_reader_t reader = md_buffered_reader_from_str(item_looped);
        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_str_writer_create(alloc),
            .in_loop = true,
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, alloc));

        EXPECT_EQ(10, sec.num_fields);
        EXPECT_EQ(5, sec.num_rows);

        EXPECT_STREQ("_entity.id", str_ptr(sec.headers[0]));
        EXPECT_STREQ("1", str_ptr(mmcif_section_value(&sec, 0, 0)));
        EXPECT_STREQ("2", str_ptr(mmcif_section_value(&sec, 1, 0)));
        EXPECT_STREQ("3", str_ptr(mmcif_section_value(&sec, 2, 0)));
        EXPECT_STREQ("4", str_ptr(mmcif_section_value(&sec, 3, 0)));
        EXPECT_STREQ("5", str_ptr(mmcif_section_value(&sec, 4, 0)));

        EXPECT_STREQ("_entity.type", str_ptr(sec.headers[1]));
        EXPECT_STREQ("polymer", str_ptr(mmcif_section_value(&sec, 0, 1)));
        EXPECT_STREQ("polymer", str_ptr(mmcif_section_value(&sec, 1, 1)));
        EXPECT_STREQ("polymer", str_ptr(mmcif_section_value(&sec, 2, 1)));
        EXPECT_STREQ("polymer", str_ptr(mmcif_section_value(&sec, 3, 1)));
        EXPECT_STREQ("non-polymer", str_ptr(mmcif_section_value(&sec, 4, 1)));

        EXPECT_STREQ("_entity.src_method", str_ptr(sec.headers[2]));
        EXPECT_STREQ("man", str_ptr(mmcif_section_value(&sec, 0, 2)));
        EXPECT_STREQ("man", str_ptr(mmcif_section_value(&sec, 1, 2)));
        EXPECT_STREQ("man", str_ptr(mmcif_section_value(&sec, 2, 2)));
        EXPECT_STREQ("man", str_ptr(mmcif_section_value(&sec, 3, 2)));
        EXPECT_STREQ("syn", str_ptr(mmcif_section_value(&sec, 4, 2)));

        EXPECT_STREQ("_entity.pdbx_description", str_ptr(sec.headers[3]));
        EXPECT_STREQ("Antiviral innate immune response receptor RIG-I", str_ptr(mmcif_section_value(&sec, 0, 3)));
        EXPECT_STREQ("E3 ubiquitin-protein ligase RNF135", str_ptr(mmcif_section_value(&sec, 1, 3)));
        EXPECT_STREQ("p3dsRNA24a", str_ptr(mmcif_section_value(&sec, 2, 3)));
        EXPECT_STREQ("p3dsRNA24b", str_ptr(mmcif_section_value(&sec, 3, 3)));
        EXPECT_STREQ("ZINC ION", str_ptr(mmcif_section_value(&sec, 4, 3)));

        EXPECT_STREQ("_entity.formula_weight", str_ptr(sec.headers[4]));
        EXPECT_STREQ("106740.555", str_ptr(mmcif_section_value(&sec, 0, 4)));
        EXPECT_STREQ("47946.305", str_ptr(mmcif_section_value(&sec, 1, 4)));
        EXPECT_STREQ("7885.581", str_ptr(mmcif_section_value(&sec, 2, 4)));
        EXPECT_STREQ("7788.551", str_ptr(mmcif_section_value(&sec, 3, 4)));
        EXPECT_STREQ("65.409", str_ptr(mmcif_section_value(&sec, 4, 4)));

        EXPECT_STREQ("_entity.pdbx_number_of_molecules", str_ptr(sec.headers[5]));
        EXPECT_STREQ("2", str_ptr(mmcif_section_value(&sec, 0, 5)));
        EXPECT_STREQ("2", str_ptr(mmcif_section_value(&sec, 1, 5)));
        EXPECT_STREQ("1", str_ptr(mmcif_section_value(&sec, 2, 5)));
        EXPECT_STREQ("1", str_ptr(mmcif_section_value(&sec, 3, 5)));
        EXPECT_STREQ("2", str_ptr(mmcif_section_value(&sec, 4, 5)));
    }

    md_temp_end(temp);
}

UTEST(mmcif, parse_2or2_comprehensive) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/2or2.cif");
    
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    bool result = md_mmcif_system_init_from_file(&sys, &sys_state, path);
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
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/nonexistent.cif");
    
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_FALSE(md_mmcif_system_init_from_file(&sys, &sys_state, path));
    
    // Should be safe to free even when init failed
    md_system_free(&sys);
    md_system_state_free(&sys_state);
    md_temp_end(temp);
}

UTEST(mmcif, advance_to_next_control) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR "/1fez.cif");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));
    char* buf = md_temp_alloc(temp, MEGABYTES(1));
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

    mmcif_parse_state_t state = {
        .reader = &reader,
        .sb = md_str_writer_create(alloc),
    };

    str_t tok;

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("data_XXXX")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("loop_")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.name")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.pdbx_ordinal")));

    EXPECT_TRUE(mmcif_advance_to_next_control_token(&state));
    EXPECT_TRUE(mmcif_next_token(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("_struct.entry_id")));

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
    for (size_t r = 0; r < ARRAY_SIZE(rows); ++r) {
        for (int i = 0; i < num_extra; ++i) {
            md_str_writer_push_cstr(&sb, "x ");
        }
        md_str_writer_push_cstr(&sb, rows[r]);
    }
    md_str_writer_push_cstr(&sb, "#\n");

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_mmcif_system_init_from_str(&sys, &sys_state, md_str_writer_to_str(&sb)));
    ASSERT_EQ(3, sys.atom.count);
    EXPECT_STREQ("CA", str_ptr(md_atom_name(&sys.atom, 1)));
    EXPECT_NEAR(7.0f, sys_state.xyz[2].x, 1.0e-5f);
    EXPECT_NEAR(8.0f, sys_state.xyz[2].y, 1.0e-5f);
    EXPECT_NEAR(9.0f, sys_state.xyz[2].z, 1.0e-5f);

    md_temp_end(temp);
}

// 9gfu misses residues within both of its chains (_pdbx_unobs_or_zero_occ_residues): label_seq_id 40-41, 67-72 and
// 259-260 of chain A, and 62-71 of chain B. Each leaves a break between the residues on either side, whose numbers (the
// residues of the sequence) count the missing ones exactly. The residues missing at the ends of the chains leave none.
UTEST(mmcif, chain_breaks) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/9gfu.cif")));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

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
        if (found < ARRAY_SIZE(expected)) {
            const md_instance_idx_t inst = md_system_instance_find_by_residue_idx(&sys, i);
            EXPECT_TRUE(str_eq_cstr(md_system_instance_id(&sys, (size_t)inst), expected[found].chain));
            EXPECT_EQ(expected[found].before, md_system_residue_seq_id(&sys, i));
            EXPECT_EQ(expected[found].after,  md_system_residue_seq_id(&sys, i + 1));
            EXPECT_TRUE(md_system_residue_flags(&sys, i + 1) & MD_RESIDUE_FLAG_BREAK_BEFORE);
            EXPECT_EQ(inst, md_system_instance_find_by_residue_idx(&sys, i + 1));
            md_atom_idx_t atoms[2];
            EXPECT_TRUE(md_util_chain_break_atoms(atoms, &sys, i));
        }
        ++found;
    }
    EXPECT_EQ(ARRAY_SIZE(expected), found);

    md_vm_arena_destroy(alloc);
}

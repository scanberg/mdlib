#include "utest.h"
#include <math.h>

#include <md_mmcif.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_allocator.h>

#define MAX_VALIDATION_SAMPLES 100

UTEST(mmcif, 1fez) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1fez.cif");

    md_system_t mol;
    bool result = md_mmcif_system_loader()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);

    if (result) {
        EXPECT_EQ(4097, mol.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&mol.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&mol.atom, 0)));
        //EXPECT_STREQ("LYS", mol.atom.resname[0].buf);
        //EXPECT_EQ(1, mol.atom.resid[0]);
        //EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(52.489, mol.atom.x[0], 0.001);
        EXPECT_NEAR(21.292, mol.atom.y[0], 0.001);
        EXPECT_NEAR(84.339, mol.atom.z[0], 0.001);
    }

    md_system_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 2or2) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/2or2.cif");

    md_system_t mol;
    bool result = md_mmcif_system_loader()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);
    //md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_UTIL_POSTPROCESS_ALL);

    if (result) {
        EXPECT_EQ(5382, mol.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&mol.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&mol.atom, 0)));
        //EXPECT_STREQ("ALA", mol.atom.resname[0].buf);
        //EXPECT_EQ(1, mol.atom.resid[0]);
        //EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(58.157, mol.atom.x[0], 0.001);
        EXPECT_NEAR(49.822, mol.atom.y[0], 0.001);
        EXPECT_NEAR(80.569, mol.atom.z[0], 0.001);
    }

    md_system_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 8g7u) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/8g7u.cif");

    md_system_t mol;
    bool result = md_mmcif_system_loader()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
    EXPECT_TRUE(result);
    md_util_molecule_postprocess(&mol, md_get_heap_allocator(), MD_UTIL_POSTPROCESS_ALL);

    if (result) {
        EXPECT_EQ(14229, mol.atom.count);

        EXPECT_EQ(7, md_atom_atomic_number(&mol.atom, 0));
        EXPECT_STREQ("N", str_ptr(md_atom_name(&mol.atom, 0)));
        //EXPECT_STREQ("PHE", mol.atom.resname[0].buf);
        //EXPECT_EQ(241, mol.atom.resid[0]);
        //EXPECT_STREQ("A", mol.atom.chainid[0].buf);

        EXPECT_NEAR(77.862,  mol.atom.x[0], 0.001);
        EXPECT_NEAR(105.453, mol.atom.y[0], 0.001);
        EXPECT_NEAR(80.951,  mol.atom.z[0], 0.001);
    }

    md_system_free(&mol, md_get_heap_allocator());
}

#include <md_mmcif.c>

UTEST(mmcif, tokenizer) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif");
        md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
        char* buf = md_vm_arena_push(alloc, MEGABYTES(1));
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_strb_create(alloc),
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

        md_file_close(file);
    }

    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif");
        md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
        char* buf = md_vm_arena_push(alloc, MEGABYTES(1));
        md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_strb_create(alloc),
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
        
        md_file_close(file);
    }


    // TODO Test parsing of the above section

    md_vm_arena_destroy(alloc);
}

UTEST(mmcif, parse_section) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    str_t section = STR_LIT(
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
            .sb = md_strb_create(alloc),
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, false, alloc));
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

    str_t section_looped = STR_LIT(
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
        md_buffered_reader_t reader = md_buffered_reader_from_str(section_looped);
        mmcif_parse_state_t state = {
            .reader = &reader,
            .sb = md_strb_create(alloc),
        };

        mmcif_section_t sec = {0};
        ASSERT_TRUE(mmcif_parse_section(&sec, &state, true, alloc));

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

    md_vm_arena_destroy(alloc);
}

UTEST(mmcif, parse_2or2_comprehensive) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/2or2.cif");
    
    md_system_t mol = {0};
    bool result = md_mmcif_system_loader()->init_from_file(&mol, path, NULL, alloc);
    ASSERT_TRUE(result);
    
    // Check basic structure properties
    EXPECT_GT(mol.atom.count, 0);
    EXPECT_GT(mol.comp.count, 0);
    EXPECT_GT(mol.inst.count, 0);
    
    // Check that coordinates are reasonable (not all zeros or infinities)
    bool has_nonzero_coord = false;
    for (int64_t i = 0; i < mol.atom.count && i < MAX_VALIDATION_SAMPLES; ++i) {
        EXPECT_FALSE(isnan(mol.atom.x[i]));
        EXPECT_FALSE(isnan(mol.atom.y[i]));
        EXPECT_FALSE(isnan(mol.atom.z[i]));
        EXPECT_FALSE(isinf(mol.atom.x[i]));
        EXPECT_FALSE(isinf(mol.atom.y[i]));
        EXPECT_FALSE(isinf(mol.atom.z[i]));
        
        if (mol.atom.x[i] != 0.0f || mol.atom.y[i] != 0.0f || mol.atom.z[i] != 0.0f) {
            has_nonzero_coord = true;
        }
    }
    EXPECT_TRUE(has_nonzero_coord);
    
    md_system_free(&mol, alloc);
}

UTEST(mmcif, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/nonexistent.cif");
    
    md_system_t mol = {0};
    bool result = md_mmcif_system_loader()->init_from_file(&mol, path, NULL, alloc);
    EXPECT_FALSE(result);
    
    // Should be safe to free even when init failed
    md_system_free(&mol, alloc);
}

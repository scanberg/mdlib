#include "utest.h"

#include <md_mmcif.h>
#include <md_molecule.h>
#include <md_util.h>
#include <core/md_allocator.h>

UTEST(mmcif, 1fez) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/1fez.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
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

    md_molecule_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 2or2) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/2or2.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
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

    md_molecule_free(&mol, md_get_heap_allocator());
}

UTEST(mmcif, 8g7u) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR"/8g7u.cif");

    md_molecule_t mol;
    bool result = md_mmcif_molecule_api()->init_from_file(&mol, path, NULL, md_get_heap_allocator());
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

    md_molecule_free(&mol, md_get_heap_allocator());
}

#include <core/md_parse.h>
#include <core/md_str_builder.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>

typedef struct {
    md_buffered_reader_t* reader;
    str_t cur_line;
    md_strb_t sb;
} mmcif_parse_state_t;

// Returns true if token was extracted,
// False if EOF / stream
static bool mmcif_token_next(str_t* out_tok, mmcif_parse_state_t* s) {
    str_t tok = {0};
    while (true) {
        if (str_empty(s->cur_line)) {
            if (!md_buffered_reader_extract_line(&s->cur_line, s->reader)) {
                return false;
            }
        }
        if (!extract_token(&tok, &s->cur_line) || tok.ptr[0] == '#') {
            s->cur_line = (str_t){0};
            continue;
        }

        if (tok.ptr[0] == '\'' || tok.ptr[0] == '\"') {
            char delim = tok.ptr[0];
            if (tok.len > 1 && tok.ptr[tok.len-1] == delim) {
                tok.ptr +=1;
                tok.len -=2;
            } else {
                size_t loc = 0;
                if (!str_find_char(&loc, s->cur_line, delim)) {
                    MD_LOG_ERROR("Unmatched string");
                    return false;
                }
                tok.ptr  += 1;
                tok.len  += loc;
                s->cur_line.ptr += (loc+1);
                s->cur_line.len -= (loc+1);
            }
        } else if (tok.ptr[0] == ';') {
            md_strb_reset(&s->sb);
            md_strb_push_str(&s->sb, str_substr(tok, 1, SIZE_MAX));
            md_strb_push_str(&s->sb, s->cur_line);
            md_strb_push_char(&s->sb, '\n');

            bool closed = false;
            while (md_buffered_reader_extract_line(&s->cur_line, s->reader)) {
                if (s->cur_line.len > 0 && s->cur_line.ptr[0] == ';') {
                    md_strb_pop(&s->sb, 1); // remove newline char
                    closed = true;
                    s->cur_line = (str_t) {0};
                    break;
                }
                md_strb_push_str (&s->sb, s->cur_line);
                md_strb_push_char(&s->sb, '\n');
            }
            if (!closed) {
                MD_LOG_ERROR("Unterminated multiline string");
                return false;
            }
            tok = md_strb_to_str(s->sb);
        }
        break;
    }

    *out_tok = tok;
    return true;
}

UTEST(mmcif, tokenizer) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif");

    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    char* buf = md_vm_arena_push(alloc, MEGABYTES(1));
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, MEGABYTES(1), file);

    mmcif_parse_state_t state = {
        .reader = &reader,
        .cur_line = {0},
        .sb = md_strb_create(alloc),
    };

    str_t tok;
    
    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("data_XXXX")));

    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("loop_")));

    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.name")));

    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("_audit_author.pdbx_ordinal")));

    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("Morais, M.C.")));

    ASSERT_TRUE(mmcif_token_next(&tok, &state));
    EXPECT_TRUE(str_eq(tok, STR_LIT("1")));

    for (size_t i = 0; i < 1000; ++i) {
        ASSERT_TRUE(mmcif_token_next(&tok, &state));
    }
}

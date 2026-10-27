#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <stdio.h>

#include "vial_spec.h"

#define VIAL_REFERENCE_PATH MD_UNITTEST_SOURCE_DIR "/../docs/vial_reference.md"
#define VIAL_PRELUDE_PATH   MD_UNITTEST_SOURCE_DIR "/../src/vial/prelude.vial"

struct vial_spec {
    md_allocator_i* alloc;
    str_t text;
    vial_spec_t spec;
};

UTEST_F_SETUP(vial_spec) {
    utest_fixture->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    utest_fixture->text  = load_textfile(STR_LIT(VIAL_REFERENCE_PATH), utest_fixture->alloc);
    ASSERT_FALSE(str_empty(utest_fixture->text));
    vial_spec_read(&utest_fixture->spec, utest_fixture->text, utest_fixture->alloc);
}

UTEST_F_TEARDOWN(vial_spec) {
    md_arena_allocator_destroy(utest_fixture->alloc);
}

// The reference reads as a test suite: blocks are found, and every comment that starts like an expectation follows
// the grammar of "How to read this reference".
UTEST_F(vial_spec, reads) {
    const vial_spec_t* spec = &utest_fixture->spec;
    const size_t num_blocks = md_array_size(spec->blocks);
    EXPECT_GE(num_blocks, 20);
    EXPECT_EQ(0, spec->num_unparsed);

    size_t num_result = 0, num_error = 0, num_warning = 0, num_filter_blocks = 0, num_cont = 0;
    size_t by_structure[3] = {0};  // T, P, none
    for (size_t i = 0; i < num_blocks; ++i) {
        const vial_block_t* b = &spec->blocks[i];
        num_filter_blocks += b->mode == VIAL_BLOCK_FILTERS;
        num_cont += b->continues;
        by_structure[b->structure == 'T' ? 0 : b->structure == 'P' ? 1 : 2] += 1;
        for (size_t j = 0; j < md_array_size(b->expects); ++j) {
            switch (b->expects[j].kind) {
            case VIAL_EXPECT_RESULT:  num_result  += 1; break;
            case VIAL_EXPECT_ERROR:   num_error   += 1; break;
            case VIAL_EXPECT_WARNING: num_warning += 1; break;
            default: break;
            }
        }
    }
    printf("vial reference: %zu blocks (%zu on T, %zu on P, %zu other; %zu as filters, %zu continuing the one before)\n",
        num_blocks, by_structure[0], by_structure[1], by_structure[2], num_filter_blocks, num_cont);
    printf("vial reference: %zu expectations: %zu results, %zu errors, %zu with warnings\n",
        num_result + num_error + num_warning, num_result, num_error, num_warning);
    printf("vial reference: %zu codes in the diagnostics table\n", md_array_size(spec->diag_codes));
    EXPECT_GT(num_result, 100);
    EXPECT_GT(num_error, 10);
}

// Every code an example expects is in the diagnostics table
UTEST_F(vial_spec, codes_are_listed) {
    const vial_spec_t* spec = &utest_fixture->spec;
    EXPECT_GT(md_array_size(spec->diag_codes), 30);
    for (size_t i = 0; i < md_array_size(spec->blocks); ++i) {
        const vial_block_t* b = &spec->blocks[i];
        for (size_t j = 0; j < md_array_size(b->expects); ++j) {
            const vial_expect_t* e = &b->expects[j];
            for (int k = 0; k < e->num_codes; ++k) {
                const bool listed = vial_spec_has_code(spec, e->codes[k]);
                if (!listed) printf("line %d: %s is not in the diagnostics table\n", e->line, e->codes[k]);
                EXPECT_TRUE(listed);
            }
        }
    }
}

// src/vial/prelude.vial is the prelude the reference prints, character for character
UTEST_F(vial_spec, prelude_matches) {
    const vial_spec_t* spec = &utest_fixture->spec;
    const vial_block_t* prelude = NULL;
    for (size_t i = 0; i < md_array_size(spec->blocks); ++i) {
        if (str_eq(spec->blocks[i].section, STR_LIT("The prelude"))) { prelude = &spec->blocks[i]; break; }
    }
    ASSERT_TRUE(prelude != NULL);
    str_t file = load_textfile(STR_LIT(VIAL_PRELUDE_PATH), utest_fixture->alloc);
    ASSERT_FALSE(str_empty(file));
    const bool same = str_eq(str_trim(file), str_trim(prelude->src));
    if (!same) printf("src/vial/prelude.vial differs from the prelude in the reference (line %d)\n", prelude->line);
    EXPECT_TRUE(same);
}

// The grammar of expectations, on the forms the reference uses
UTEST(vial_spec, grammar) {
    vial_expect_t e;

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("-> float[433][3] Å, frame: a comment")));
    EXPECT_EQ(VIAL_EXPECT_RESULT, e.kind);
    EXPECT_TRUE(str_eq(e.result, STR_LIT("float[433][3] Å")));
    EXPECT_TRUE(str_eq(e.note, STR_LIT("frame: a comment")));

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("-> 'MET1/A'")));
    EXPECT_TRUE(str_eq(e.result, STR_LIT("'MET1/A'")));

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("error E0301: [?] and [433] do not broadcast")));
    EXPECT_EQ(VIAL_EXPECT_ERROR, e.kind);
    EXPECT_EQ(1, e.num_codes);
    EXPECT_STREQ("E0301", e.codes[0]);

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("warning W0601, W0405 -> float Å | missing: empty on T")));
    EXPECT_EQ(VIAL_EXPECT_WARNING, e.kind);
    EXPECT_EQ(2, e.num_codes);
    EXPECT_STREQ("W0405", e.codes[1]);
    EXPECT_TRUE(str_eq(e.result, STR_LIT("float Å | missing")));

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("warning W0405; R0101 and a missing value")));
    EXPECT_EQ(1, e.num_codes);
    EXPECT_TRUE(str_eq(e.note, STR_LIT("R0101 and a missing value")));

    e = (vial_expect_t){0};
    EXPECT_TRUE(vial_spec__parse_expect(&e, STR_LIT("the CAs of the group's twenty amino acids")));
    EXPECT_EQ(VIAL_EXPECT_NONE, e.kind);

    e = (vial_expect_t){0};
    EXPECT_FALSE(vial_spec__parse_expect(&e, STR_LIT("error: no code")));
    e = (vial_expect_t){0};
    EXPECT_FALSE(vial_spec__parse_expect(&e, STR_LIT("error E0301 -> float")));
}

// ### PHASE 1 ###

#include <vial/vial_internal.h>

// Whether a code is one the front end reports
static bool vial_spec__is_syntax_code(const char* code) {
    static const char* codes[] = {"E0101", "E0102", "E0103", "E0104", "E0701", "E0801", "E0803"};
    for (size_t i = 0; i < ARRAY_SIZE(codes); ++i) if (strncmp(code, codes[i], 5) == 0) return true;
    return false;
}

static const vial_expect_t* vial_spec__expect_at(const vial_block_t* b, int line_in_block) {
    for (size_t i = 0; i < md_array_size(b->expects); ++i) {
        if (b->expects[i].line_in_block == line_in_block) return &b->expects[i];
    }
    return NULL;
}

static bool vial_spec__expects_code(const vial_expect_t* e, const char* code) {
    if (!e) return false;
    for (int k = 0; k < e->num_codes; ++k) if (strncmp(e->codes[k], code, 5) == 0) return true;
    return false;
}

// Every block lexes, and the lexer reports only what a line expects
UTEST_F(vial_spec, lexes) {
    const vial_spec_t* spec = &utest_fixture->spec;
    size_t num_tokens = 0;
    for (size_t i = 0; i < md_array_size(spec->blocks); ++i) {
        const vial_block_t* b = &spec->blocks[i];
        vial_diags_t diags;
        vial_diags_init(&diags, b->src, utest_fixture->alloc);
        vial_token_t* tokens = vial_lex(b->src, &diags, utest_fixture->alloc);
        num_tokens += md_array_size(tokens);
        for (size_t j = 0; j < md_array_size(diags.items); ++j) {
            const md_vial_diagnostic_t* d = &diags.items[j];
            const int line = (int)vial_line_of(&diags, d->span.beg);
            const bool expected = vial_spec__expects_code(vial_spec__expect_at(b, line), d->code);
            if (!expected) {
                printf("reference line %d: unexpected %s %.*s\n", b->line + line - 1, d->code, (int)d->message.len, d->message.ptr);
            }
            EXPECT_TRUE(expected);
        }
    }
    printf("vial reference: %zu tokens\n", num_tokens);
}

#include <vial/ast.h>

// The phase 1 gate: every block parses, each line gives exactly the syntax codes it expects, and the prelude parses
UTEST_F(vial_spec, parses) {
    const vial_spec_t* spec = &utest_fixture->spec;
    size_t num_checked = 0, num_nodes = 0;
    for (size_t i = 0; i < md_array_size(spec->blocks); ++i) {
        const vial_block_t* b = &spec->blocks[i];
        vial_diags_t diags;
        vial_diags_init(&diags, b->src, utest_fixture->alloc);
        vial_token_t* tokens = vial_lex(b->src, &diags, utest_fixture->alloc);
        vial_tree_t tree = vial_parse(tokens, b->src, &diags, utest_fixture->alloc);
        num_nodes += md_array_size(tree.nodes);

        // Every diagnostic is one its line expects
        for (size_t j = 0; j < md_array_size(diags.items); ++j) {
            const md_vial_diagnostic_t* d = &diags.items[j];
            const int line = (int)vial_line_of(&diags, d->span.beg);
            const bool expected = vial_spec__expects_code(vial_spec__expect_at(b, line), d->code);
            if (!expected) {
                printf("reference line %d: unexpected %s %.*s\n", b->line + line - 1, d->code, (int)d->message.len, d->message.ptr);
            }
            EXPECT_TRUE(expected);
        }
        // Every syntax code a line expects is given on it
        for (size_t j = 0; j < md_array_size(b->expects); ++j) {
            const vial_expect_t* e = &b->expects[j];
            for (int k = 0; k < e->num_codes; ++k) {
                if (!vial_spec__is_syntax_code(e->codes[k])) continue;
                bool found = false;
                for (size_t q = 0; q < md_array_size(diags.items); ++q) {
                    const md_vial_diagnostic_t* d = &diags.items[q];
                    if ((int)vial_line_of(&diags, d->span.beg) == e->line_in_block && strncmp(d->code, e->codes[k], 5) == 0) found = true;
                }
                if (!found) printf("reference line %d: expected %s, not given\n", e->line, e->codes[k]);
                EXPECT_TRUE(found);
                num_checked += 1;
            }
        }
    }

    // The prelude, as the compiler will read it
    str_t prelude = load_textfile(STR_LIT(VIAL_PRELUDE_PATH), utest_fixture->alloc);
    vial_diags_t diags;
    vial_diags_init(&diags, prelude, utest_fixture->alloc);
    vial_token_t* tokens = vial_lex(prelude, &diags, utest_fixture->alloc);
    vial_parse(tokens, prelude, &diags, utest_fixture->alloc);
    for (size_t j = 0; j < md_array_size(diags.items); ++j) {
        printf("prelude.vial line %u: %s %.*s\n", vial_line_of(&diags, diags.items[j].span.beg), diags.items[j].code,
            (int)diags.items[j].message.len, diags.items[j].message.ptr);
    }
    EXPECT_EQ(0, md_array_size(diags.items));

    printf("vial reference: every block parses into %zu nodes; %zu expected syntax codes given\n", num_nodes, num_checked);
}

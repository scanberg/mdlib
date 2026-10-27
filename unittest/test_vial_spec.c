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

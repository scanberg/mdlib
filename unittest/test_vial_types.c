// The types of the checker, and the harness that compares them with the reference's results

#include "utest.h"

#include <core/md_arena_allocator.h>
#include <core/md_common.h>

#include <stdio.h>

#include "vial_spec.h"
#include <vial/types.h>

#define VIAL_TYPES_REFERENCE MD_UNITTEST_SOURCE_DIR "/../docs/vial_reference.md"

static md_unit_t unit_of(const char* s) {
    md_unit_t u;
    md_unit_parse(&u, str_from_cstr(s));
    return u;
}

static const char* print(const vial_type_t* t) {
    static char buf[256];
    vial_type_print(buf, sizeof(buf), t);
    return buf;
}

// The printer writes the reference's notation
UTEST(vial_types, print) {
    vial_type_t t = vial_type_scalar(VIAL_DT_FLOAT);
    EXPECT_STREQ("float", print(&t));

    t.unit = unit_of(u8"Å");
    EXPECT_STREQ(u8"float Å", print(&t));
    vial_type_push_dim(&t, vial_extent_static(433));
    vial_type_push_dim(&t, vial_extent_static(3));
    EXPECT_STREQ(u8"float[433][3] Å", print(&t));

    vial_type_t sq = vial_type_scalar(VIAL_DT_FLOAT);
    sq.unit = md_unit_mul(unit_of(u8"Å"), unit_of(u8"Å"));
    vial_type_push_dim(&sq, vial_extent_static(45));
    EXPECT_STREQ(u8"float[45] Å²", print(&sq));

    vial_type_t s = vial_type_scalar(VIAL_DT_SET);
    vial_type_push_dim(&s, vial_extent_dynamic(7));
    EXPECT_STREQ("set[?]", print(&s));
    vial_type_t s2 = vial_type_scalar(VIAL_DT_SET);
    vial_type_push_dim(&s2, vial_extent_static(14));
    vial_type_push_dim(&s2, vial_extent_static(4));
    EXPECT_STREQ("set[14][4]", print(&s2));

    vial_type_t m = vial_type_scalar(VIAL_DT_FLOAT);
    m.unit = unit_of(u8"Å");
    m.flags |= VIAL_TYPE_MISSING;
    EXPECT_STREQ(u8"float Å | missing", print(&m));

    vial_type_t r = vial_type_scalar(VIAL_DT_FLOAT);
    r.unit = unit_of(u8"Å");
    vial_type_push_dim(&r, (vial_extent_t){.kind = VIAL_EXT_RANGED, .n = 4, .m = 14});
    vial_type_push_dim(&r, vial_extent_static(3));
    EXPECT_STREQ(u8"float[4..14][3] Å", print(&r));

    vial_type_t q = vial_type_scalar(VIAL_DT_FLOAT);
    q.unit_kind = VIAL_UNIT_UNKNOWN;
    EXPECT_STREQ("float ?", print(&q));

    vial_type_t deg = vial_type_scalar(VIAL_DT_FLOAT);
    deg.unit = unit_of("deg");
    EXPECT_STREQ(u8"float °", print(&deg));

    vial_type_t e = vial_type_scalar(VIAL_DT_FLOAT);
    e.unit = unit_of("kJ/mol");
    EXPECT_STREQ("float kJ/mol", print(&e));

    vial_type_t t64 = vial_type_scalar(VIAL_DT_FLOAT64);
    t64.unit = unit_of("ps");
    EXPECT_STREQ("float64 ps", print(&t64));

    // A tuple, as closest_pair gives
    vial_type_t el[3] = {vial_type_scalar(VIAL_DT_FLOAT), vial_type_scalar(VIAL_DT_INT), vial_type_scalar(VIAL_DT_INT)};
    el[0].unit = unit_of(u8"Å");
    vial_type_t tup = vial_type_scalar(VIAL_DT_TUPLE);
    tup.elems = el;
    tup.num_elems = 3;
    EXPECT_STREQ(u8"(float Å, int, int)", print(&tup));

    // A seventh dimension does not fit
    vial_type_t big = vial_type_scalar(VIAL_DT_FLOAT);
    for (int i = 0; i < VIAL_MAX_DIMS; ++i) EXPECT_TRUE(vial_type_push_dim(&big, vial_extent_static(2)));
    EXPECT_FALSE(vial_type_push_dim(&big, vial_extent_static(2)));

    vial_type_t tr = vial_type_scalar(VIAL_DT_FLOAT);
    tr.stage = VIAL_STAGE_TRAJECTORY;
    EXPECT_STREQ("trajectory", vial_stage_str(&tr));
    tr.flags |= VIAL_TYPE_PER_FRAME;
    EXPECT_STREQ("trajectory, per frame", vial_stage_str(&tr));
}

UTEST(vial_types, matches) {
    EXPECT_TRUE(vial_spec_type_matches(STR_LIT(u8"float[433][3] Å"), STR_LIT(u8"float[433][3] Å")));
    EXPECT_TRUE(vial_spec_type_matches(STR_LIT(u8"float[n][3] Å"), STR_LIT(u8"float[6784][3] Å")));
    EXPECT_TRUE(vial_spec_type_matches(STR_LIT(u8"float[45] Å²"), STR_LIT(u8"float[45] Å^2")));
    EXPECT_FALSE(vial_spec_type_matches(STR_LIT(u8"float[433] Å"), STR_LIT(u8"float[432] Å")));
    EXPECT_FALSE(vial_spec_type_matches(STR_LIT("set[?]"), STR_LIT("set[433]")));
    EXPECT_FALSE(vial_spec_type_matches(STR_LIT(u8"float Å"), STR_LIT(u8"float Å | missing")));
}

// Every result in the reference reads as types or a value, with its stage and pass where written
UTEST(vial_types, results_read) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    str_t text = load_textfile(STR_LIT(VIAL_TYPES_REFERENCE), alloc);
    ASSERT_FALSE(str_empty(text));
    vial_spec_t spec;
    vial_spec_read(&spec, text, alloc);
    int num_types = 0, num_values = 0, num_stages = 0, num_passes = 0;
    for (size_t i = 0; i < md_array_size(spec.blocks); ++i) {
        const vial_block_t* b = &spec.blocks[i];
        for (size_t j = 0; j < md_array_size(b->expects); ++j) {
            const vial_expect_t* e = &b->expects[j];
            if (str_empty(e->result)) continue;
            vial_result_t r;
            const bool ok = vial_spec_read_result(&r, e);
            if (!ok) printf("reference line %d: cannot read the result '%.*s'\n", e->line, (int)e->result.len, e->result.ptr);
            EXPECT_TRUE(ok);
            if (r.is_value) num_values += 1; else num_types += r.num_types;
            num_stages += r.has_stage;
            num_passes += r.pass > 0;
        }
    }
    printf("vial reference: %d types, %d values, %d stages and %d passes to check\n", num_types, num_values, num_stages, num_passes);
    EXPECT_GT(num_types, 150);
    md_arena_allocator_destroy(alloc);
}

#include "utest.h"
#include <string.h>

#include <md_gro.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_array.h>

#define NM_TO_ANGSTROM 10.0f

UTEST(gro, parse_small) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/catalyst.gro");
    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, path, alloc));
    EXPECT_EQ(gro_data.num_atoms, 1336);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_gro_system_init_from_data(&sys, &sys_state, &gro_data);
    for (int64_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, gro_data.atom_data[i].x * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].y, gro_data.atom_data[i].y * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].z, gro_data.atom_data[i].z * NM_TO_ANGSTROM);
    }
    md_system_reset(&sys);

    EXPECT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, path));
    // A format without force field types still keeps the column aligned, with empty entries
    EXPECT_EQ(sys.atom.type.count, md_array_size(sys.atom.type.ff_type));
    for (size_t t = 0; t < sys.atom.type.count; ++t) {
        EXPECT_TRUE(str_empty(md_atom_type_ff_type(&sys.atom.type, t)));
    }
    for (int64_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, gro_data.atom_data[i].x * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].y, gro_data.atom_data[i].y * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].z, gro_data.atom_data[i].z * NM_TO_ANGSTROM);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, parse_big) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/centered.gro");
    md_gro_data_t gro_data = { 0 };
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, path, alloc));
    EXPECT_EQ(gro_data.num_atoms, 161742);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_gro_system_init_from_data(&sys, &sys_state, &gro_data);
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, gro_data.atom_data[i].x * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].y, gro_data.atom_data[i].y * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].z, gro_data.atom_data[i].z * NM_TO_ANGSTROM);
    }
    md_system_reset(&sys);

    EXPECT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, path));
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(sys_state.xyz[i].x, gro_data.atom_data[i].x * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].y, gro_data.atom_data[i].y * NM_TO_ANGSTROM);
        EXPECT_EQ(sys_state.xyz[i].z, gro_data.atom_data[i].z * NM_TO_ANGSTROM);
    }

    md_system_free(&sys);
    md_system_state_free(&sys_state);
    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, parse_small_water) {
    md_allocator_i* alloc = md_get_heap_allocator();

    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/water.gro");
    md_gro_data_t gro_data = {0};
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_TRUE(result);
    EXPECT_EQ(gro_data.num_atoms, 12165);

    // Check that atoms have reasonable coordinates (should be in a 5x5x5 box)
    for (size_t i = 0; i < 100 && i < gro_data.num_atoms; ++i) { // Test first 100 atoms
        EXPECT_GT(gro_data.atom_data[i].x, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].x, 6.0f);
        EXPECT_GT(gro_data.atom_data[i].y, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].y, 6.0f);
        EXPECT_GT(gro_data.atom_data[i].z, -1.0f);
        EXPECT_LT(gro_data.atom_data[i].z, 6.0f);
    }

    md_gro_data_free(&gro_data, alloc);
}

UTEST(gro, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t path = STR_INIT(MD_UNITTEST_DATA_DIR"/nonexistent.gro");
    md_gro_data_t gro_data = {0};
    bool result = md_gro_data_parse_file(&gro_data, path, alloc);
    EXPECT_FALSE(result);
    
    md_gro_data_free(&gro_data, alloc);
}

// ---------------------------------------------------------------------------
// The coordinate columns
//
// A line that keeps to the columns of the first atom line is read field by field; one that does not
// is read by its tokens. Either way the values are the ones parse_float gives the numbers written.
// ---------------------------------------------------------------------------

#include <core/md_parse.h>
#include <math.h>

static bool gro_test_same(float a, float b) {
    uint32_t ua, ub;
    memcpy(&ua, &a, 4);
    memcpy(&ub, &b, 4);
    return ua == ub;
}

static bool gro_test_absent(float v) {
    uint32_t u;
    memcpy(&u, &v, 4);
    return (u & 0x7FFFFFFFu) > 0x7F800000u;
}

UTEST(gro, coordinate_columns) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const str_t text = STR_LIT(
        "columns, and lines that do not keep to them\n"
        "    7\n"
        "    1SOL     OW    1   0.126   1.624   1.679  0.1427 -0.4840  0.0544\n"
        "    1SOL    HW1    2  -1.234-123.456 999.999\n"
        "    1SOL    HW2    3    0.5     1.25    2.0\n"
        "    2SOL     OW    4   0.126   1.624   1.679  0.1427 -0.4840\n"
        "    2SOL    HW1    5   0.126   1.624   1.679  \n"
        "    2SOL    HW2    6 -1234.567   1.000   2.000\n"
        "    3SOL     OW    7   1.000   2.000   3.000  0.1000  0.2000  0.3000  extra\n"
        "   1.86206   1.86206   1.86206\n");

    md_gro_data_t gro = {0};
    ASSERT_TRUE(md_gro_data_parse_str(&gro, text, alloc));
    ASSERT_EQ(7, gro.num_atoms);
    const md_gro_atom_t* a = gro.atom_data;

    // The columns, with velocities one decimal finer
    EXPECT_TRUE(gro_test_same(a[0].x, 0.126f));
    EXPECT_TRUE(gro_test_same(a[0].z, 1.679f));
    EXPECT_TRUE(gro_test_same(a[0].vx, 0.1427f));
    EXPECT_TRUE(gro_test_same(a[0].vy, -0.4840f));
    // Fields that touch are still the columns
    EXPECT_TRUE(gro_test_same(a[1].x, -1.234f));
    EXPECT_TRUE(gro_test_same(a[1].y, -123.456f));
    EXPECT_TRUE(gro_test_same(a[1].z, 999.999f));
    EXPECT_TRUE(gro_test_absent(a[1].vx));
    // Not the columns: read by tokens
    EXPECT_TRUE(gro_test_same(a[2].x, 0.5f));
    EXPECT_TRUE(gro_test_same(a[2].y, 1.25f));
    EXPECT_TRUE(gro_test_same(a[2].z, 2.0f));
    // Two velocities are none, trailing spaces are nothing
    EXPECT_TRUE(gro_test_same(a[3].x, 0.126f));
    EXPECT_TRUE(gro_test_absent(a[3].vx));
    EXPECT_TRUE(gro_test_absent(a[3].vz));
    EXPECT_TRUE(gro_test_same(a[4].z, 1.679f));
    EXPECT_TRUE(gro_test_absent(a[4].vx));
    // A value wider than its field pushes the others out of their columns: tokens again
    EXPECT_TRUE(gro_test_same(a[5].x, -1234.567f));
    EXPECT_TRUE(gro_test_same(a[5].y, 1.0f));
    EXPECT_TRUE(gro_test_same(a[5].z, 2.0f));
    // What follows the velocities is not read
    EXPECT_TRUE(gro_test_same(a[6].z, 3.0f));
    EXPECT_TRUE(gro_test_same(a[6].vz, 0.3f));

    EXPECT_TRUE(gro_test_same(gro.box[0][0], 1.86206f));
    md_gro_data_free(&gro, alloc);
}

// gmx editconf -ndec: wider columns with more decimals, the velocities one more again
UTEST(gro, coordinate_columns_with_more_decimals) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const str_t text = STR_LIT(
        "five decimals\n"
        "    2\n"
        "    1SOL     OW    1   0.12600   1.62400  -1.67900  0.142700 -0.484000  0.054400\n"
        "    1SOL    HW1    2 -12.34567 123.45678   0.00001\n"
        "   1.86206   1.86206   1.86206\n");

    md_gro_data_t gro = {0};
    ASSERT_TRUE(md_gro_data_parse_str(&gro, text, alloc));
    ASSERT_EQ(2, gro.num_atoms);
    const md_gro_atom_t* a = gro.atom_data;
    EXPECT_TRUE(gro_test_same(a[0].x, 0.126f));
    EXPECT_TRUE(gro_test_same(a[0].z, -1.679f));
    EXPECT_TRUE(gro_test_same(a[0].vx, 0.1427f));
    EXPECT_TRUE(gro_test_same(a[0].vz, 0.0544f));
    EXPECT_TRUE(gro_test_same(a[1].x, -12.34567f));
    EXPECT_TRUE(gro_test_same(a[1].y, 123.45678f));
    EXPECT_TRUE(gro_test_same(a[1].z, 0.00001f));
    EXPECT_TRUE(gro_test_absent(a[1].vx));
    md_gro_data_free(&gro, alloc);
}

// The reference files, every coordinate against parse_float on its token
static size_t gro_test_compare_with_tokens(str_t path) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_gro_data_t gro = {0};
    if (!md_gro_data_parse_file(&gro, path, alloc)) return SIZE_MAX;
    str_t text = load_textfile(path, alloc);
    str_t rest = text, line;
    str_extract_line(&line, &rest);
    str_extract_line(&line, &rest);
    size_t mismatches = 0;
    for (size_t i = 0; i < gro.num_atoms && str_extract_line(&line, &rest); ++i) {
        str_t coords = str_substr(line, 20, SIZE_MAX);
        str_t tok[6];
        const size_t n = extract_tokens(tok, 6, &coords);
        const md_gro_atom_t* a = &gro.atom_data[i];
        const float* v[6] = { &a->x, &a->y, &a->z, &a->vx, &a->vy, &a->vz };
        for (size_t k = 0; k < n; ++k) {
            if (!gro_test_same(*v[k], (float)parse_float(tok[k]))) mismatches += 1;
        }
    }
    str_free(text, alloc);
    md_gro_data_free(&gro, alloc);
    return mismatches;
}

UTEST(gro, coordinate_columns_match_the_tokens) {
    EXPECT_EQ(0, gro_test_compare_with_tokens(STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro")));
    EXPECT_EQ(0, gro_test_compare_with_tokens(STR_LIT(MD_UNITTEST_DATA_DIR "/water.gro")));
    EXPECT_EQ(0, gro_test_compare_with_tokens(STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
}

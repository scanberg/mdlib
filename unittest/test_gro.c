#include "utest.h"
#include <string.h>

#include <md_gro.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_file.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <md_attributes.h>

#include "test_executor.h"

#define NM_TO_ANGSTROM 10.0f

// A coordinate of the file in nm as the loader keeps it, in Angstrom
static inline float nm(float x) {
    return x * NM_TO_ANGSTROM;
}

// The velocities of a loaded system, in nm/ps as the file has them: NULL where the file has none
static const float* gro_test_velocities(const md_system_t* sys) {
    return (const float*)md_attribute_view(md_attributes_find(&sys->attributes, MD_STR_LIT("atom/velocity")), MD_ATTRIBUTE_TYPE_F32, 3, 1);
}

// The reference files, every coordinate (and velocity) of the system against md_parse_float on its token, defined later
static size_t gro_test_compare_with_tokens(md_str_t path, const md_executor_t* ex);

UTEST(gro, parse_small) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/catalyst.gro");

    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, path, NULL));
    EXPECT_EQ((size_t)1336, sys.atom.count);
    EXPECT_EQ((size_t)1336, sys_state.num_atoms);
    // A format without force field types still keeps the column aligned, with empty entries
    EXPECT_EQ(sys.atom.type.count, md_array_size(sys.atom.type.ff_type));
    for (size_t t = 0; t < sys.atom.type.count; ++t) {
        EXPECT_TRUE(md_str_empty(md_particle_type_ff_type(&sys.atom.type, t)));
    }
    EXPECT_EQ((size_t)0, gro_test_compare_with_tokens(path, NULL));

    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

UTEST(gro, parse_big) {
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/centered.gro");
    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, path, NULL));
    EXPECT_EQ((size_t)161742, sys.atom.count);
    md_system_free(&sys);
    md_system_state_free(&sys_state);
    EXPECT_EQ((size_t)0, gro_test_compare_with_tokens(path, NULL));
}

UTEST(gro, parse_small_water) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR"/water.gro");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, path, NULL));
    EXPECT_EQ((size_t)12165, sys.atom.count);
    EXPECT_EQ((size_t)4055, sys.residue.count);

    // Reasonable coordinates: in a box of 5 nm
    for (size_t i = 0; i < 100 && i < sys.atom.count; ++i) {
        EXPECT_GT(sys_state.xyz[i].x, -10.0f);
        EXPECT_LT(sys_state.xyz[i].x, 60.0f);
        EXPECT_GT(sys_state.xyz[i].y, -10.0f);
        EXPECT_LT(sys_state.xyz[i].y, 60.0f);
        EXPECT_GT(sys_state.xyz[i].z, -10.0f);
        EXPECT_LT(sys_state.xyz[i].z, 60.0f);
    }
    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

UTEST(gro, nonexistent_file) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    EXPECT_FALSE(md_gro_system_init_from_file(&sys, &sys_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR"/nonexistent.gro"), NULL));
    EXPECT_EQ((size_t)0, sys.atom.count);
    md_system_free(&sys);
    md_system_state_free(&sys_state);
}

// What is not a GRO file, or ends before it says it would: no system, and the state of none
UTEST(gro, malformed) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const char* texts[] = {
        "",
        "title only\n",
        "title\nabc\n    1SOL     OW    1   0.126   1.624   1.679\n   1.0   1.0   1.0\n",
        "title\n    3\n    1SOL     OW    1   0.126   1.624   1.679\n   1.0   1.0   1.0\n",
        "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n",
    };
    for (size_t i = 0; i < MD_COUNTOF(texts); ++i) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        EXPECT_FALSE(md_gro_system_init_from_str(&sys, &st, md_str_from_cstr(texts[i]), NULL));
        EXPECT_EQ((size_t)0, sys.atom.count);
        EXPECT_EQ((size_t)0, st.num_atoms);
    }
    md_vm_arena_destroy(alloc);
}

// A file of several frames: the first is the system
UTEST(gro, frames) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    const md_str_t text = MD_STR_LIT(
        "frame 1\n"
        "    2\n"
        "    1SOL     OW    1   0.100   0.200   0.300\n"
        "    1SOL    HW1    2   0.400   0.500   0.600\n"
        "   1.00000   2.00000   3.00000\n"
        "frame 2\n"
        "    2\n"
        "    1SOL     OW    1   9.100   9.200   9.300\n"
        "    1SOL    HW1    2   9.400   9.500   9.600\n"
        "   4.00000   5.00000   6.00000\n");
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, text, NULL));
    ASSERT_EQ((size_t)2, sys.atom.count);
    EXPECT_EQ((size_t)1, sys.residue.count);
    EXPECT_NEAR(4.0f, st.xyz[1].x, 1e-5f);
    EXPECT_NEAR(10.0, st.unitcell.x, 1e-5);
    EXPECT_NEAR(30.0, st.unitcell.z, 1e-5);
    md_vm_arena_destroy(alloc);
}

// ---------------------------------------------------------------------------
// The coordinate columns
//
// A line that keeps to the columns of the first atom line is read field by field; one that does not
// is read by its tokens. Either way the values are the ones md_parse_float gives the numbers written.
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
    const md_str_t text = MD_STR_LIT(
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

    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, text, NULL));
    ASSERT_EQ((size_t)7, sys.atom.count);
    const md_float3_t* x = st.xyz;
    const float* v = gro_test_velocities(&sys);
    ASSERT_TRUE(v != NULL);

    // The columns, with velocities one decimal finer
    EXPECT_TRUE(gro_test_same(x[0].x, nm(0.126f)));
    EXPECT_TRUE(gro_test_same(x[0].z, nm(1.679f)));
    EXPECT_TRUE(gro_test_same(v[0 * 3 + 0], 0.1427f));
    EXPECT_TRUE(gro_test_same(v[0 * 3 + 1], -0.4840f));
    // Fields that touch are still the columns
    EXPECT_TRUE(gro_test_same(x[1].x, nm(-1.234f)));
    EXPECT_TRUE(gro_test_same(x[1].y, nm(-123.456f)));
    EXPECT_TRUE(gro_test_same(x[1].z, nm(999.999f)));
    EXPECT_TRUE(gro_test_absent(v[1 * 3 + 0]));
    // Not the columns: read by tokens
    EXPECT_TRUE(gro_test_same(x[2].x, nm(0.5f)));
    EXPECT_TRUE(gro_test_same(x[2].y, nm(1.25f)));
    EXPECT_TRUE(gro_test_same(x[2].z, nm(2.0f)));
    // Two velocities are none, trailing spaces are nothing
    EXPECT_TRUE(gro_test_same(x[3].x, nm(0.126f)));
    EXPECT_TRUE(gro_test_absent(v[3 * 3 + 0]));
    EXPECT_TRUE(gro_test_absent(v[3 * 3 + 2]));
    EXPECT_TRUE(gro_test_same(x[4].z, nm(1.679f)));
    EXPECT_TRUE(gro_test_absent(v[4 * 3 + 0]));
    // A value wider than its field pushes the others out of their columns: tokens again
    EXPECT_TRUE(gro_test_same(x[5].x, nm(-1234.567f)));
    EXPECT_TRUE(gro_test_same(x[5].y, nm(1.0f)));
    EXPECT_TRUE(gro_test_same(x[5].z, nm(2.0f)));
    // What follows the velocities is not read
    EXPECT_TRUE(gro_test_same(x[6].z, nm(3.0f)));
    EXPECT_TRUE(gro_test_same(v[6 * 3 + 2], 0.3f));

    EXPECT_NEAR(18.6206, st.unitcell.x, 1e-4);
    md_system_free(&sys);
    md_system_state_free(&st);
}

// gmx editconf -ndec: wider columns with more decimals, the velocities one more again
UTEST(gro, coordinate_columns_with_more_decimals) {
    md_allocator_i* alloc = md_get_heap_allocator();
    const md_str_t text = MD_STR_LIT(
        "five decimals\n"
        "    2\n"
        "    1SOL     OW    1   0.12600   1.62400  -1.67900  0.142700 -0.484000  0.054400\n"
        "    1SOL    HW1    2 -12.34567 123.45678   0.00001\n"
        "   1.86206   1.86206   1.86206\n");

    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, text, NULL));
    ASSERT_EQ((size_t)2, sys.atom.count);
    const md_float3_t* x = st.xyz;
    const float* v = gro_test_velocities(&sys);
    ASSERT_TRUE(v != NULL);
    EXPECT_TRUE(gro_test_same(x[0].x, nm(0.126f)));
    EXPECT_TRUE(gro_test_same(x[0].z, nm(-1.679f)));
    EXPECT_TRUE(gro_test_same(v[0], 0.1427f));
    EXPECT_TRUE(gro_test_same(v[2], 0.0544f));
    EXPECT_TRUE(gro_test_same(x[1].x, nm(-12.34567f)));
    EXPECT_TRUE(gro_test_same(x[1].y, nm(123.45678f)));
    EXPECT_TRUE(gro_test_same(x[1].z, nm(0.00001f)));
    EXPECT_TRUE(gro_test_absent(v[3]));
    md_system_free(&sys);
    md_system_state_free(&st);
}

// The reference files, every coordinate against md_parse_float on its token (in Angstrom, as the loader makes it), every
// velocity against its own
static size_t gro_test_compare_with_tokens(md_str_t path, const md_executor_t* ex) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    size_t mismatches = SIZE_MAX;
    if (md_gro_system_init_from_file(&sys, &st, path, ex)) {
        mismatches = 0;
        const float* vel = gro_test_velocities(&sys);
        md_str_t text = md_load_textfile(path, alloc);
        md_str_t rest = text, line;
        md_parse_next_line(&line, &rest);
        md_parse_next_line(&line, &rest);
        for (size_t i = 0; i < sys.atom.count && md_parse_next_line(&line, &rest); ++i) {
            md_str_t coords = md_str_substr(line, 20, SIZE_MAX);
            md_str_t tok[6];
            const size_t n = md_parse_next_tokens(tok, 6, &coords);
            const float xyz[3] = { st.xyz[i].x, st.xyz[i].y, st.xyz[i].z };
            for (size_t k = 0; k < n && k < 3; ++k) {
                if (!gro_test_same(xyz[k], nm((float)md_parse_float(tok[k])))) mismatches += 1;
            }
            for (size_t k = 3; vel && k < n; ++k) {
                if (!gro_test_same(vel[i * 3 + k - 3], (float)md_parse_float(tok[k]))) mismatches += 1;
            }
        }
    }
    md_vm_arena_destroy(alloc);
    return mismatches;
}

UTEST(gro, coordinate_columns_match_the_tokens) {
    test_executor_t threads = test_executor_threads(3);
    const char* files[] = { MD_UNITTEST_DATA_DIR "/catalyst.gro", MD_UNITTEST_DATA_DIR "/water.gro", MD_UNITTEST_DATA_DIR "/centered.gro", MD_UNITTEST_DATA_DIR "/npt.gro" };
    for (size_t f = 0; f < MD_COUNTOF(files); ++f) {
        EXPECT_EQ((size_t)0, gro_test_compare_with_tokens(md_str_from_cstr(files[f]), NULL));
        EXPECT_EQ((size_t)0, gro_test_compare_with_tokens(md_str_from_cstr(files[f]), test_executor(&threads)));
    }
}

// ---------------------------------------------------------------------------
// Straight from the text, in blocks and pieces
//
// A file of more than one block (md_reader.h) whose residues straddle the blocks and the pieces, with velocities on
// some lines and none on others: the same system read on the calling thread, on threads, and with the pieces out of
// order; from the file and from a string.
// ---------------------------------------------------------------------------

#include "test_load_blocks.h"
#include <core/md_unit.h>
#include <stdio.h>

UTEST(gro, blocks_and_pieces) {
    const size_t num_waters = 300000;   // 900000 lines of about 68 characters: two blocks
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(4));
    md_array(char) text = 0;
    char line[128];
    int n = snprintf(line, sizeof(line), "waters\n%5zu\n", num_waters * 3);
    md_array_push_array(text, line, (size_t)n, alloc);
    static const char* names[3] = { "OW", "HW1", "HW2" };
    for (size_t w = 0; w < num_waters; ++w) {
        for (int a = 0; a < 3; ++a) {
            const size_t i = w * 3 + a;
            const float x = (float)(i % 997) * 0.01f, y = (float)(i % 991) * 0.01f, z = (float)(i % 983) * 0.01f;
            // GROMACS numbers residues and atoms modulo 100000
            if (w % 7 == 3) {
                n = snprintf(line, sizeof(line), "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", (int)((w + 1) % 100000), "SOL", names[a], (int)((i + 1) % 100000), x, y, z);
            } else {
                n = snprintf(line, sizeof(line), "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", (int)((w + 1) % 100000), "SOL", names[a], (int)((i + 1) % 100000), x, y, z, 0.1f * a, -0.2f, 0.3f);
            }
            md_array_push_array(text, line, (size_t)n, alloc);
        }
    }
    n = snprintf(line, sizeof(line), "   9.97000   9.91000   9.83000\n");
    md_array_push_array(text, line, (size_t)n, alloc);
    const md_str_t str = { text, md_array_size(text) };
    ASSERT_GT(str.len, (size_t)MD_MEGABYTES(33));

    const char* path = "md_unittest_blocks.gro";
    ASSERT_TRUE(load_blocks_write(path, str));

    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&ref, &ref_st, str, NULL));
    EXPECT_EQ(num_waters * 3, ref.atom.count);
    EXPECT_EQ(num_waters, ref.residue.count);

    test_executor_t threads  = test_executor_threads(4);
    test_executor_t shuffled = test_executor_shuffled(3);
    const md_executor_t* executors[] = { NULL, test_executor(&threads), test_executor(&shuffled) };
    for (size_t e = 0; e < MD_COUNTOF(executors); ++e) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, md_str_from_cstr(path), executors[e]));
        EXPECT_EQ((size_t)0, load_blocks_diff(&ref, &ref_st, &sys, &st));
        md_system_t sys2 = { .alloc = alloc };
        md_system_state_t st2 = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_str(&sys2, &st2, str, executors[e]));
        EXPECT_EQ((size_t)0, load_blocks_diff(&ref, &ref_st, &sys2, &st2));
    }

    // The velocities: of every line which has them, absent where a line has none
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, md_str_from_cstr(path), test_executor(&threads)));
    const md_attribute_t* vel = md_attributes_find(&sys.attributes, MD_STR_LIT("atom/velocity"));
    ASSERT_TRUE(vel != NULL);
    const md_unit_t nm_per_ps = md_unit_div(md_unit_nanometer(), md_unit_picosecond());
    float v[3] = {0};
    // Water 3 has none, water 4 has them
    ASSERT_EQ((size_t)3, md_attribute_extract_f32(v, 3, vel, md_attribute_slice_1(3 * 3 + 1), nm_per_ps));
    EXPECT_TRUE(gro_test_absent(v[0]));
    ASSERT_EQ((size_t)3, md_attribute_extract_f32(v, 3, vel, md_attribute_slice_1(4 * 3 + 1), nm_per_ps));
    EXPECT_TRUE(gro_test_same(v[0], 0.1f));
    EXPECT_TRUE(gro_test_same(v[1], -0.2f));

    remove(path);
    md_vm_arena_destroy(alloc);
}

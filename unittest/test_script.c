#include "utest.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <core/md_sync.h>
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_gro.h>
#include <md_pdb.h>

#include <md_script.c>

// Create molecule for evaulation
#define ATOM_COUNT 16
static float x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
static float y[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
static float z[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
static float r[] = {1,2,3,4,4,4,5,1,1,2,3,4,4,4,5,1};
static float m[] = {1,2,2,2,2,4,4,4,1,2,2,2,2,4,4,4};
static uint8_t e[] = {1,8,1,2,6,8,6,8,1,8,1,2,6,7,6,8};
static md_label_t n[] = {
    MAKE_LABEL("H"),
    MAKE_LABEL("O"),
    MAKE_LABEL("H"),
    MAKE_LABEL("He"),
    MAKE_LABEL("C"),
    MAKE_LABEL("N"),
    MAKE_LABEL("CA"),
    MAKE_LABEL("O"),
    MAKE_LABEL("H"),
    MAKE_LABEL("O"),
    MAKE_LABEL("H"),
    MAKE_LABEL("He"),
    MAKE_LABEL("C"),
    MAKE_LABEL("N"),
    MAKE_LABEL("CA"),
    MAKE_LABEL("O")
};
static md_residue_idx_t r_idx[] = {0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3};
static md_chain_idx_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

#define RES_COUNT 4
static md_label_t r_name[] = {MAKE_LABEL("SOL"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT"), MAKE_LABEL("PFT")};
static md_residue_id_t r_id[] = {1, 2, 3, 4};
static md_range_t r_range[] = {{0, 3}, {3, 8}, {8,12}, {12,16}};

#define CHAIN_COUNT 1
static md_label_t c_id[] = {MAKE_LABEL("A")};
static md_range_t c_arange[] = {0,16};
static md_range_t c_rrange[] = {0,4};

md_molecule_t test_mol = {
    .atom = {
        .count = ATOM_COUNT,
        .x = x,
        .y = y,
        .z = z,
        .radius = r,
        .mass = m,
        .element = e,
        .name = n,
        .residue_idx = r_idx,
        .chain_idx = c_idx
},
.residue = {
        .count = RES_COUNT,
        .name = r_name,
        .id = r_id,
        .atom_range = r_range
},
.chain = {
        .count = CHAIN_COUNT,
        .id = c_id,
        .atom_range = c_arange,
        .residue_range = c_rrange
}
};

static bool eval_selection(md_bitfield_t* bitfield, str_t expr, md_molecule_t* mol) {
    ASSERT(bitfield);
    ASSERT(mol);
    data_t data = {0};
    if (eval_expression(&data, expr, mol, default_temp_allocator)) {
        if (data.type.base_type == TYPE_BITFIELD) {
            md_bitfield_t* res = (md_bitfield_t*)data.ptr;
            const int64_t len = type_info_array_len(data.type);
            for (int64_t i = 0; i < len; ++i) {
                md_bitfield_or_inplace(bitfield, &res[i]);
            }
            return true;
        }
    }
    return false;
}

static uint64_t make_bits(const char* bit_str) {    
    uint64_t bits = 0;
    const uint64_t num_bits = strlen(bit_str);
    ASSERT(num_bits < 64);
    for (uint64_t i = 0; i < num_bits; i++) {
        const uint64_t bit_mask = (1LLU << i);
        if (bit_str[i] != '0')
            bits |= bit_mask;
        else
            bits &= ~bit_mask;
    }
    return bits;
}

static void print_bits(uint64_t* bits, uint64_t num_bits) {
    for (uint64_t i = 0; i < num_bits; ++i) {
        const uint64_t blk_idx = i / 64;
        const uint64_t bit_mask = (1LLU << i);
        printf("%i", bits[blk_idx] & bit_mask ? 1 : 0);
    }
}

UTEST(script, basic_expressions) {
    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, MAKE_STR("'this is a string'"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_STRING);
        str_t str = as_string(data);
        EXPECT_STREQ(str.ptr, "this is a string");
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, MAKE_STR("2 + 5"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(as_int(data), 7);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, MAKE_STR("2 + 5.0"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
        EXPECT_EQ(as_float(data), 7.0);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, MAKE_STR("{2,1} + {1,8}"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(data.type.dim[0], 2);
        EXPECT_EQ(as_int_arr(data)[0], 3);
        EXPECT_EQ(as_int_arr(data)[1], 9);
    }
}

static bool test_selection(const char* expr, const char* ref_bit_str) {
    bool result = false;
    uint64_t ref = make_bits(ref_bit_str);
    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, default_temp_allocator);

    if (!eval_selection(&bf, str_from_cstr(expr), &test_mol)) {
        printf("Failed evaluation of expression: '%s'\n", expr);
        goto done;
    }
    bool cmp_res = bit_cmp(bf.bits, &ref, 0, ATOM_COUNT);
    if (!cmp_res) {
        printf("Got:\n");
        print_bits(bf.bits, ATOM_COUNT);
        printf("\nExpected\n");
        print_bits(&ref, ATOM_COUNT);
        printf("\n");
        goto done;
    }
    result = true;
done:
    md_bitfield_free(&bf);
    return result;
}

UTEST(script, selection) {
    EXPECT_TRUE(test_selection("all",               "1111111111111111"));
    EXPECT_TRUE(test_selection("resname('SOL')",    "1110000000000000"));
    EXPECT_TRUE(test_selection("element('C')",      "0000101000001010"));
    EXPECT_TRUE(test_selection("label('CA')",       "0000001000000010"));
    EXPECT_TRUE(test_selection("atom(1) in resname('PFT')", "0000000010001000"));
    //TEST_SELECTION("atom(1:2) or element('O') in residue(:)", "1001000010001000");
}

UTEST(script, compile_script) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));
    const str_t gro_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/centered.gro");

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, gro_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    str_t script_src = load_textfile(MAKE_STR(MD_UNITTEST_DATA_DIR "/script.txt"), alloc);
    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_source(ir, script_src, &mol, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));

    md_arena_allocator_destroy(alloc);
}

UTEST(script, semantic) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));

    md_molecule_t mol = { 0 };
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, MAKE_STR(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc));

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        md_script_ir_compile_source(ir, MAKE_STR("p1 = resname('ALA') resname('GLY');"), &mol, NULL);
        EXPECT_FALSE(md_script_ir_valid(ir));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST(script, selection_big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));

    const str_t gro_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/centered.gro");

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro_data, gro_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);

    bool result = eval_selection(&bf, MAKE_STR("atom(1:20) and element('O') in chain(:)"), &mol);
    EXPECT_TRUE(result);

    md_arena_allocator_destroy(alloc);
}

UTEST(script, dynamic_length) {
    md_allocator_i* alloc = default_allocator;
    const str_t pdb_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_trajectory_i* traj = md_pdb_trajectory_create(pdb_file, alloc);
    ASSERT_TRUE(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = MAKE_STR("sel1 = residue(within_z(1:50));");
        md_script_ir_compile_source(ir, src, &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    md_script_ir_free(ir);
    md_pdb_trajectory_free(traj);
    md_molecule_free(&mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}

UTEST(script, property_compute) {
    md_allocator_i* alloc = default_allocator;
    const str_t pdb_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_trajectory_i* traj = md_pdb_trajectory_create(pdb_file, alloc);
    ASSERT_TRUE(traj);
    uint32_t num_frames = (uint32_t)md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = MAKE_STR("prop1 = distance_pair(com(resname(\"ALA\")), 1);");
        md_script_ir_compile_source(ir, src, &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_compile_source(ir, MAKE_STR("d1 = 1:5 in residue(1:3);"), &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_compile_source(ir, MAKE_STR("prop1 = rdf(element('C'), element('O'), 20.0);"), &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_compile_source(ir, MAKE_STR("sel = within_x(0:100);\np1  = distance(com(sel), 100);"), &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        str_t src = MAKE_STR(
            "s1 = residue(1:10);\n"
            "s2 = residue(11:15);\n"
            "s = {s1, s2};"
        );
        md_script_ir_compile_source(ir, src, &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        str_t src = MAKE_STR("s1 = count(within(10, residue(:)));");
        md_script_ir_compile_source(ir, src, &mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &mol, traj, 0, num_frames));
        md_script_eval_free(eval);
    }

    md_script_ir_free(ir);
    md_pdb_trajectory_free(traj);
    md_molecule_free(&mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}

#define NUM_THREADS 8

typedef struct thread_data_t {
    const md_script_ir_t* ir;
    const md_molecule_t* mol;
    const md_trajectory_i* traj;
    const md_script_eval_t* ref_eval;
    md_script_eval_t* eval;
    int num_corrupt_values;
} thread_data_t;

void func(void* user_data) {
    thread_data_t* data = (thread_data_t*)user_data;
    const int64_t num_cur_prop = md_script_eval_num_properties(data->eval);
    const md_script_property_t* cur_prop = md_script_eval_properties(data->eval);

    const int64_t num_ref_prop = md_script_eval_num_properties(data->ref_eval);
    const md_script_property_t* ref_prop = md_script_eval_properties(data->ref_eval);

    ASSERT(num_cur_prop == num_ref_prop);
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(data->traj);
    if (md_script_eval_frame_range(data->eval, data->ir, data->mol, data->traj, 0, num_frames)) {
        for (int64_t p_idx = 0; p_idx < num_cur_prop; ++p_idx) {
            ASSERT(cur_prop[p_idx].data.num_values == ref_prop[p_idx].data.num_values);
            for (int64_t i = 0; i < cur_prop[p_idx].data.num_values; ++i) {
                if (cur_prop[p_idx].data.values[i] != ref_prop[p_idx].data.values[i]) {
                    data->num_corrupt_values += 1;
                    fprintf(stderr, "Corruption occured in thread %"PRIu64" at frame %i, expected: '%g', got: '%g'\n", md_thread_id(), (int)i, ref_prop[p_idx].data.values[i], cur_prop[p_idx].data.values[i]);
                }
            }
        }
    }
}

UTEST(script, parallel_evaluation) {
    md_allocator_i* alloc = default_allocator;
    const str_t pdb_file = MAKE_STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
    const str_t script = MAKE_STR("p1 = distance(1,10);");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, pdb_file, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_trajectory_i* traj = md_pdb_trajectory_create(pdb_file, alloc);
    ASSERT_TRUE(traj);

    md_script_eval_t* eval[NUM_THREADS] = {0};
    md_thread_t* threads[NUM_THREADS] = {0};
    thread_data_t thread_data[NUM_THREADS] = {0};

    int64_t num_frames = md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_source(ir, script, &mol, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));

    md_script_eval_t* ref_eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(ref_eval, ir, &mol, traj, 0, (uint32_t)num_frames));

    ASSERT_EQ(1, md_script_eval_num_properties(ref_eval));
    ASSERT_EQ(num_frames, md_script_eval_properties(ref_eval)[0].data.num_values);

#if 0
    printf("Ref Values:\n");
    for (int64_t i = 0; i < ref_eval.properties[0].data.num_values; ++i) {
        printf("[%lli]: %g\n", i, ref_eval.properties[0].data.values[i]);
    }
#endif

    for (int pass = 0; pass < 10; ++pass) {
        for (int i = 0; i < NUM_THREADS; ++i) {
            eval[i] = md_script_eval_create(num_frames, ir, alloc);
            EXPECT_NE(NULL, eval[i]);
        }

        for (int i = 0; i < NUM_THREADS; ++i) {
            thread_data[i] = (thread_data_t) {
                .ir = ir,
                .mol = &mol,
                .traj = traj,
                .ref_eval = ref_eval,
                .eval = eval[i],
                .num_corrupt_values = 0
            };
            threads[i] = md_thread_create(func, &thread_data[i]);
        }

        int total_corrupt_values = 0;
        for (int i = 0; i < NUM_THREADS; ++i) {
            md_thread_join(threads[i]);
            total_corrupt_values += thread_data[i].num_corrupt_values;
        }
        EXPECT_EQ(0, total_corrupt_values);
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        md_script_eval_free(eval[i]);
    }

    md_script_ir_free(ir);
    md_script_eval_free(ref_eval);

    md_pdb_trajectory_free(traj);
    md_molecule_free(&mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}
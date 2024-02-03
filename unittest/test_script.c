#include "utest.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
#include <core/md_bitop.inl>
#include <core/md_os.h>
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
static md_label_t t[] = {
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
static md_label_t rname[] = {
    MAKE_LABEL("SOL"),
    MAKE_LABEL("SOL"),
    MAKE_LABEL("SOL"),
    MAKE_LABEL("LYS"),
    MAKE_LABEL("LYS"),
    MAKE_LABEL("LYS"),
    MAKE_LABEL("LYS"),
    MAKE_LABEL("LYS"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT"),
    MAKE_LABEL("PFT")
};
static md_chain_idx_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

#define RES_COUNT 4
static md_label_t r_name[] = {MAKE_LABEL("SOL"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT"), MAKE_LABEL("PFT")};
static md_residue_id_t r_id[] = {1, 2, 3, 4};
static uint32_t r_off[] = {0, 3, 8, 12, 16};

#define CHAIN_COUNT 1
static md_label_t c_id[] = {MAKE_LABEL("A")};
static uint32_t c_aoff[] = {0,16};
static uint32_t c_roff[] = {0,4};

md_molecule_t test_mol = {
    .atom = {
        .count = ATOM_COUNT,
        .x = x,
        .y = y,
        .z = z,
        .radius = r,
        .mass = m,
        .element = e,
        .type = t,
        .resname = rname,
},
.residue = {
        .count = RES_COUNT,
        .name = r_name,
        .id = r_id,
        .atom_offset = r_off
},
.chain = {
        .count = CHAIN_COUNT,
        .id = c_id,
        .atom_offset = c_aoff,
        .res_offset = c_roff
}
};

struct script {
    md_vm_arena_t arena;
    md_allocator_i alloc;
    md_molecule_t amy;
    md_molecule_t ala;
    md_trajectory_i* ala_traj;
};

UTEST_F_SETUP(script) {
    md_vm_arena_init(&utest_fixture->arena, GIGABYTES(4));

    utest_fixture->alloc = md_vm_arena_create_interface(&utest_fixture->arena);

    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&utest_fixture->amy, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"),   NULL, &utest_fixture->alloc));
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&utest_fixture->ala, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, &utest_fixture->alloc));

    md_util_molecule_postprocess(&utest_fixture->amy, &utest_fixture->alloc, MD_UTIL_POSTPROCESS_ALL);
    md_util_molecule_postprocess(&utest_fixture->ala, &utest_fixture->alloc, MD_UTIL_POSTPROCESS_ALL);

    utest_fixture->ala_traj = md_pdb_trajectory_create(STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), &utest_fixture->alloc);
}

UTEST_F_TEARDOWN(script) {
    md_vm_arena_free(&utest_fixture->arena);
}

static bool eval_selection(md_bitfield_t* bitfield, str_t expr, md_molecule_t* mol) {
    ASSERT(bitfield);
    ASSERT(mol);
    data_t data = {0};
    if (eval_expression(&data, expr, mol, md_temp_allocator)) {
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

UTEST(script, type) {
    {
        type_info_t a = {.base_type = TYPE_INT, .dim = {1}, .len_dim = 0};
        type_info_t b = {.base_type = TYPE_INT, .dim = {1}, .len_dim = 0};
        EXPECT_TRUE(type_info_equal(a,b));
    }
    {
        type_info_t a = {.base_type = TYPE_INT, .dim = {1,1}, .len_dim = 0};
        type_info_t b = {.base_type = TYPE_INT, .dim = {1,0}, .len_dim = 0};
        EXPECT_TRUE(type_info_equal(a,b));
    }
    {
        type_info_t a = {.base_type = TYPE_INT, .dim = {1,0}, .len_dim = 0};
        type_info_t b = {.base_type = TYPE_INT, .dim = {2,0}, .len_dim = 0};
        EXPECT_FALSE(type_info_equal(a,b));
    }
}

UTEST(script, basic_expressions) {
    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("'this is a string'"), &test_mol, md_temp_allocator);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_STRING);
            str_t str = as_string(data);
            EXPECT_STREQ(str.ptr, "this is a string");
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("2 + 5"), &test_mol, md_temp_allocator);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_INT);
            EXPECT_EQ(as_int(data), 7);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("2 + 5.0"), &test_mol, md_temp_allocator);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
            EXPECT_EQ(as_float(data), 7.0);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("{2,1} + {1,8}"), &test_mol, md_temp_allocator);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_INT);
            EXPECT_EQ(data.type.dim[0], 2);
            EXPECT_EQ(as_int_arr(data)[0], 3);
            EXPECT_EQ(as_int_arr(data)[1], 9);
        }
    }
}

ast_node_t* parse_and_type_check_expression(str_t expr, md_script_ir_t* ir, md_molecule_t* mol, md_vm_arena_t* temp_arena) {
    // @HACK: We use alloc here: If the data type is a str_t, then it gets a shallow copy
    // Which means that the actual string data is contained within the ir->arena => temp_alloc
    ir->str = str_copy(expr, ir->arena);

    md_allocator_i temp_alloc = md_vm_arena_create_interface(temp_arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    bool result = false;

    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = &temp_alloc});
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .mol = mol,
            .temp_arena = temp_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
        };

        if (static_check_node(node, &ctx)) {
            return node;
        }
    }

    if (ir->errors) {
        for (int64_t i = 0; i < md_array_size(ir->errors); ++i) {
            MD_LOG_ERROR("%.*s", ir->errors[i].text.len, ir->errors[i].text.ptr);
        }
    }
    return NULL;
}

UTEST(script, assignment) {
    SETUP_TEMP_ALLOC(GIGABYTES(4));
    md_script_ir_t* ir = create_ir(&temp_alloc);

    {
        md_script_ir_clear(ir);

        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,2}"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("a"));
            EXPECT_NE(NULL, ident);
            EXPECT_NE(NULL, ident->data);
        }
    }

    {
        md_script_ir_clear(ir);

        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,distance(1,2)}"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* a = get_identifier(ir, STR_LIT("a"));
            EXPECT_NE(NULL, a);
            EXPECT_NE(NULL, a->data);

            identifier_t* b = get_identifier(ir, STR_LIT("b"));
            EXPECT_NE(NULL, b);
            EXPECT_NE(NULL, b->data);
        }
    }

    FREE_TEMP_ALLOC();
}

static bool test_selection(const char* expr, const char* ref_bit_str) {
    bool result = false;
    uint64_t ref = make_bits(ref_bit_str);
    md_bitfield_t bf = {0};
    md_bitfield_init(&bf, md_temp_allocator);

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

UTEST_F(script, compile_script) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    str_t script_src = load_textfile(STR_LIT(MD_UNITTEST_DATA_DIR "/script.txt"), alloc);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, script_src, &utest_fixture->amy, NULL, NULL));
    EXPECT_TRUE(md_script_ir_valid(ir));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, semantic) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->amy;

    md_script_ir_t* ir = md_script_ir_create(alloc);

    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("p1 = resname('ALA') resname('GLY');"), mol, NULL, NULL));
    
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "ads = residue({10629,10633,10635,10637,10653,10659,10661,10665,10678:10679,10684,10686});"
        "ads_end = residue({10628:10629,10633,10635:10638,10643,10647,10650,10653:10654,10656,10659,10661,10663:10665,10669:10671,10678:10679,10682,10684,10686});"
        "not_ads = residue(not ads_end and not protein);"
        "super_ads = {residue(10633),residue(10637),residue(10665),residue(10678),residue(10679)};"
        "dih1_pft = dihedral(22, 20, 1, 2) in resname('PFT');"
        "dih2_pft = dihedral(2, 3, 6, 10) in resname('PFT');"
        "dih3_pft = dihedral(29, 27, 9, 10) in resname('PFT');"
        "dih4_pft = dihedral(35, 33, 31 , 29) in resname('PFT');"
        "dih_center_pft={dih2_pft,dih3_pft};"
        "dih1_not_ads = dihedral(22, 20, 1, 2) in not_ads;"
        "dih2_not_ads = dihedral(2, 3, 6, 10) in not_ads;"
        "dih3_not_ads = dihedral(29, 27, 9, 10) in not_ads;"
        "dih4_not_ads = dihedral(35, 33, 31 , 29) in not_ads;"
        "dih_center_not_ads = {dih2_not_ads,dih3_not_ads};"
    ), mol, NULL, NULL));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, selection_big) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->amy;

    md_bitfield_t bf = md_bitfield_create(alloc);
    EXPECT_TRUE(eval_selection(&bf, STR_LIT("atom(1:20) and element('O') in chain(:)"), mol));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, dynamic_length) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->ala;
    md_trajectory_i* traj = utest_fixture->ala_traj;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = STR_LIT("sel1 = residue(within_z(1:50));");
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, property_compute) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->ala;
    md_trajectory_i* traj = utest_fixture->ala_traj;
    uint32_t num_frames = (uint32_t)md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("num = count(residue(resname('ALA') and within(3.0, protein)));");
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("prop1 = distance_pair(com(resname(\"ALA\")), 1);");
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("d1 = 1:5 in residue(1:3);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("V = sdf(residue(1), element('H'), 5.0) * 1;"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("prop1 = rdf(element('C'), element('O'), 20.0);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("sel = within_x(0:100);\np1  = distance(com(sel), 100);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(md_script_eval_num_properties(eval), 1);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        str_t src = STR_LIT(
            "s1 = residue(1:10);\n"
            "s2 = residue(11:15);\n"
            "s = {s1, s2};"
        );

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        str_t src = STR_LIT("s1 = count(within(10, residue(:)));");

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));
        md_script_eval_free(eval);
    }

    md_arena_allocator_destroy(alloc);
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

UTEST_F(script, parallel_evaluation) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->ala;
    md_trajectory_i* traj = utest_fixture->ala_traj;

    const str_t script = STR_LIT("p1 = distance(1,10);");

    md_script_eval_t* eval[NUM_THREADS] = {0};
    md_thread_t* threads[NUM_THREADS] = {0};
    thread_data_t thread_data[NUM_THREADS] = {0};

    int64_t num_frames = md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, script, mol, traj, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));

    md_script_eval_t* ref_eval = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
    ASSERT_TRUE(md_script_eval_frame_range(ref_eval, ir, mol, traj, 0, (uint32_t)num_frames));

    ASSERT_EQ(1, md_script_eval_num_properties(ref_eval));
    ASSERT_EQ(num_frames, md_script_eval_properties(ref_eval)[0].data.num_values);

#if 0
    printf("Ref Values:\n");
    const md_script_property_t* ref_prop = md_script_eval_properties(ref_eval);
    for (int64_t i = 0; i < ref_prop->data.num_values; ++i) {
        printf("[%lli]: %g\n", i, ref_prop->data.values[i]);
    }
#endif

    for (int pass = 0; pass < 10; ++pass) {
        for (int i = 0; i < NUM_THREADS; ++i) {
            eval[i] = md_script_eval_create(num_frames, ir, STR_LIT(""), alloc);
            EXPECT_NE(NULL, eval[i]);
        }

        for (int i = 0; i < NUM_THREADS; ++i) {
            thread_data[i] = (thread_data_t) {
                .ir = ir,
                .mol = mol,
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

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, parse_unary_binary) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sqrt(2) * -4;"), mol, NULL, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5-4;"), mol, NULL, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = -4;"), mol, NULL, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = (-4);"), mol, NULL, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5 * (-4);"), mol, NULL, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5 * -4;"), mol, NULL, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = (5) - 4;"), mol, NULL, NULL));
    }

    md_arena_allocator_destroy(alloc);
}
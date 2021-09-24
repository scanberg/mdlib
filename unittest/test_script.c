#include "utest.h"
#include <md_script.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>
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
static const char* n[] = {"H", "O", "H", "He", "C", "N", "CA", "O", "H", "O", "H", "He", "C", "N", "CA", "O"};
static md_residue_idx_t r_idx[] = {0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3};
static md_chain_idx_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

#define RES_COUNT 4
static const char* r_name[] = {"SOL", "LYS", "PFT", "PFT"};
static md_residue_id_t r_id[] = {1, 2, 3, 4};
static md_range_t r_range[] = {{0, 3}, {3, 8}, {8,12}, {12,16}};

#define CHAIN_COUNT 1
static const char* c_id[] = {"A"};
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

static bool eval_selection(md_exp_bitfield_t* bitfield, str_t expr, md_molecule_t* mol) {
    ASSERT(bitfield);
    ASSERT(mol);
    data_t data = {0};
    if (eval_expression(&data, expr, mol, default_temp_allocator)) {
        if (data.type.base_type == TYPE_BITFIELD) {
            md_exp_bitfield_t* res = (md_exp_bitfield_t*)data.ptr;
            const int64_t len = type_info_array_len(data.type);
            for (int64_t i = 0; i < len; ++i) {
                md_bitfield_or(bitfield, bitfield, &res[i]);
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
        EXPECT_TRUE(eval_expression(&data, make_cstr("'this is a string'"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_STRING);
        EXPECT_STREQ(as_string(data).ptr, "this is a string");
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("2 + 5"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(as_int(data), 7);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("2 + 5.0"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
        EXPECT_EQ(as_float(data), 7.0);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("{2,1} + {1,8}"), &test_mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(data.type.dim[0], 2);
        EXPECT_EQ(as_int_arr(data)[0], 3);
        EXPECT_EQ(as_int_arr(data)[1], 9);
    }
}

#define TEST_SELECTION(expr, ref_bit_str) \
{ \
uint64_t ref = make_bits(ref_bit_str); \
md_exp_bitfield_t bf = {0}; \
md_bitfield_init(&bf, default_temp_allocator); \
ASSERT_TRUE(eval_selection(&bf, make_cstr(expr), &test_mol)); \
bool cmp_res = bit_cmp(bf.bits, &ref, 0, ATOM_COUNT); \
EXPECT_TRUE(cmp_res); \
if (!cmp_res) { \
    printf("Got:\n"); \
    print_bits(bf.bits, ATOM_COUNT); \
    printf("\nExpected\n"); \
    print_bits(&ref, ATOM_COUNT); \
    printf("\n"); \
} \
md_bitfield_free(&bf); \
}

UTEST(script, selection) {
    TEST_SELECTION("all",               "1111111111111111");
    TEST_SELECTION("resname('SOL')",    "1110000000000000");
    TEST_SELECTION("element('C')",      "0000101000001010");
    TEST_SELECTION("label('CA')",       "0000001000000010");
    TEST_SELECTION("atom(1) in resname('PFT')", "0000000010001000");
    //TEST_SELECTION("atom(1:2) or element('O') in residue(:)", "1001000010001000");
}

UTEST(script, compile_script) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));
    const str_t gro_file = make_cstr(MD_UNITTEST_DATA_DIR "/centered.gro");

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(gro_file, &gro_data, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    str_t script_src = load_textfile(make_cstr(MD_UNITTEST_DATA_DIR "/script.txt"), alloc);
    md_script_ir_compile_args_t args = {
        .src = script_src,
        .mol = &mol,
        .alloc = default_allocator
    };
    md_script_ir_t ir = {0};
    EXPECT_TRUE(md_script_ir_compile(&ir, args));

    md_arena_allocator_destroy(alloc);
}

UTEST(script, selection_big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));

    const str_t gro_file = make_cstr(MD_UNITTEST_DATA_DIR "/centered.gro");

    md_gro_data_t gro_data = {0};
    ASSERT_TRUE(md_gro_data_parse_file(gro_file, &gro_data, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro_data, alloc));

    md_exp_bitfield_t bf = {0};
    md_bitfield_init(&bf, alloc);

    bool result = eval_selection(&bf, make_cstr("atom(1:20) and element('O') in chain(:)"), &mol);
    EXPECT_TRUE(result);

    md_arena_allocator_destroy(alloc);
}

UTEST(script, property_compute) {
    md_allocator_i* alloc = default_allocator;
    const str_t pdb_file = make_cstr(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(pdb_file, &pdb_data, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_trajectory_i traj = {0};
    ASSERT_TRUE(md_pdb_trajectory_open(&traj, pdb_file, alloc));

    md_script_ir_t ir = {0};
    md_script_eval_t eval = {0};
    {
        md_script_ir_compile_args_t compile_args = {
            .src = make_cstr("prop1 = rdf(element('C'), element('O'), 20.0);"),
            .mol = &mol,
            .alloc = alloc
        };
        EXPECT_TRUE(md_script_ir_compile(&ir, compile_args));

        EXPECT_TRUE(md_script_eval_alloc(&eval, traj.num_frames, &ir, alloc));

        md_script_eval_args_t eval_args = {
            .ir = &ir,
            .mol = &mol,
            .traj = &traj,
        };
        ASSERT_TRUE(md_script_eval_compute(&eval, eval_args));

        EXPECT_EQ(eval.num_properties, 1);
    }

    {
        md_script_ir_compile_args_t compile_args = {
            .src = make_cstr(
                "sel = x(0:10);\n"
                "p1  = distance(com(sel), 100);"
                ),
            .mol = &mol,
            .alloc = alloc
        };
        EXPECT_TRUE(md_script_ir_compile(&ir, compile_args));

        EXPECT_TRUE(md_script_eval_alloc(&eval, traj.num_frames, &ir, alloc));

        md_script_eval_args_t eval_args = {
            .ir = &ir,
            .mol = &mol,
            .traj = &traj,
        };
        ASSERT_TRUE(md_script_eval_compute(&eval, eval_args));

        EXPECT_EQ(eval.num_properties, 1);
        const md_script_property_t* props = eval.properties;

        //EXPECT_EQ(props->data.num_values, traj_header.num_frames);
    }

    {
        md_script_ir_compile_args_t compile_args = {
            .src = make_cstr(
                "d1 = distance(10:2, 100);"
            ),
            .mol = &mol,
            .alloc = alloc
        };
        EXPECT_FALSE(md_script_ir_compile(&ir, compile_args));
    }

    {
        md_script_ir_compile_args_t compile_args = {
            .src = make_cstr(
                "s1 = residue(1:10);\n"
                "s2 = residue(11:15);\n"
                "s = {s1, s2};"
            ),
            .mol = &mol,
            .alloc = alloc
        };
        EXPECT_TRUE(md_script_ir_compile(&ir, compile_args));
    }

    md_script_ir_free(&ir);
    md_script_eval_free(&eval);

    md_pdb_trajectory_close(&traj);
    md_pdb_molecule_free(&mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}
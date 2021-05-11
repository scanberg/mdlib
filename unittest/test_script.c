#include "utest.h"
#include <md_script.h>
#include <md_molecule.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitfield.h>

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
static uint32_t r_idx[] = {0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3};
static uint32_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

#define RES_COUNT 4
static const char* r_name[] = {"SOL", "LYS", "PFT", "PFT"};
static uint32_t r_id[] = {1, 2, 3, 4};
static md_range r_range[] = {{0, 3}, {3, 8}, {8,12}, {12,16}};

#define CHAIN_COUNT 1
static const char* c_id[] = {"A"};
static md_range c_arange[] = {0,16};
static md_range c_rrange[] = {0,4};

static md_molecule mol = {
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



static bool eval_selection(md_bitfield_t* bitfield, str_t expr, md_molecule* mol) {
    data_t data = {0};
    if (eval_expression(&data, expr, mol, default_temp_allocator)) {
        if (is_scalar(data.type) && data.type.base_type == TYPE_BITFIELD) {
            md_bitfield_t* res = (md_bitfield_t*)data.ptr;
            bitfield->bits = res->bits;
            bitfield->num_bits = res->num_bits;
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
        const uint64_t blk_idx = i / 64;
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
        EXPECT_TRUE(eval_expression(&data, make_cstr("'this is a string'"), &mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_STRING);
        EXPECT_STREQ(as_string(data).ptr, "this is a string");
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("2 + 5"), &mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(as_int(data), 7);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("2 + 5.0"), &mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
        EXPECT_EQ(as_float(data), 7.0);
    }

    {
        data_t data = {0};
        EXPECT_TRUE(eval_expression(&data, make_cstr("{2,1} + {1,8}"), &mol, default_temp_allocator));
        EXPECT_EQ(data.type.base_type, TYPE_INT);
        EXPECT_EQ(data.type.dim[0], 2);
        EXPECT_EQ(as_int_arr(data)[0], 3);
        EXPECT_EQ(as_int_arr(data)[1], 9);
    }
}

#define TEST_SELECTION(expr, ref_bit_str) \
{ \
uint64_t ref = make_bits(ref_bit_str); \
md_bitfield_t bf = {0}; \
ASSERT_TRUE(eval_selection(&bf, make_cstr(expr), &mol)); \
bool cmp_res = bit_cmp(bf.bits, &ref, 0, ATOM_COUNT); \
EXPECT_TRUE(cmp_res); \
if (!cmp_res) { \
    printf("Got:\n"); \
    print_bits(bf.bits, bf.num_bits); \
    printf("\nExpected\n"); \
    print_bits(&ref, bf.num_bits); \
    printf("\n"); \
} \
}

UTEST(script, selection) {
    TEST_SELECTION("all",               "1111111111111111");
    TEST_SELECTION("resname('SOL')",    "1110000000000000");
    TEST_SELECTION("element('C')",      "0000101000001010");
    TEST_SELECTION("label('CA')",       "0000001000000010");
    TEST_SELECTION("atom(1) within resname('PFT')", "0000000010001000");
}

UTEST(script, compile_script) {

    const char script_file_raw[] = MD_UNITTEST_DATA_DIR "/script.txt";
    const str_t script_file = {script_file_raw, ARRAY_SIZE(script_file_raw)};

    str_t script_src = load_textfile(script_file, default_allocator);

    struct md_script_o* ir = md_script_compile(script_src, &mol, default_allocator);

    free_str(script_src, default_allocator);
}

#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_pdb.h>

UTEST(script, property_compute) {
    md_allocator_i* alloc = default_allocator;
    const char pdb_file_raw[] = MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb";
    const str_t pdb_file = {pdb_file_raw, ARRAY_SIZE(pdb_file_raw)};

    struct md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(pdb_file, &pdb_data, alloc));

    struct md_pdb_molecule_t pdb_mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&pdb_mol, &pdb_data, alloc));

    struct md_trajectory_i* traj = md_pdb_trajectory_open(pdb_file, alloc);
    ASSERT_NE(traj, NULL);

    md_trajectory_header_t traj_header = {0};
    EXPECT_TRUE(md_trajectory_extract_header(traj, &traj_header));

    str_t src = make_cstr("prop1 = rdf(element('C'), element('O'), 20.0);");

    struct md_script_o* script = md_script_compile(src, &pdb_mol.mol, alloc);
    ASSERT_TRUE(md_script_evaluate(script, &pdb_mol.mol, traj, alloc));
    EXPECT_EQ(md_script_get_num_properties(script), 1);
    const md_script_property_t* props = md_script_get_properties(script);

    EXPECT_EQ(props->data.num_values, traj_header.num_frames);

    md_script_free(script);

    md_pdb_trajectory_close(traj);
    md_pdb_molecule_free(&pdb_mol, alloc);
    md_pdb_data_free(&pdb_data, alloc);
}
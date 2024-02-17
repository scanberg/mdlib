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
static float mol_x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
static float mol_y[] = {4,3,2,1,4,3,2,1,1,1,1,1,1,1,1,1};
static float mol_z[] = {3,2,1,4,3,2,1,2,2,2,2,2,2,2,2,2};
static float mol_r[] = {1,2,3,4,4,4,5,1,1,2,3,4,4,4,5,1};
static float mol_m[] = {1,2,2,2,2,4,4,4,1,2,2,2,2,4,4,4};
static uint8_t mol_e[] = {1,8,1,2,6,8,6,8,1,8,1,2,6,7,6,8};
static md_label_t mol_t[] = {
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
        .x = mol_x,
        .y = mol_y,
        .z = mol_z,
        .radius = mol_r,
        .mass = mol_m,
        .element = mol_e,
        .type = mol_t,
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
    bool initialized;
    md_vm_arena_t arena;
    md_allocator_i alloc;
    md_molecule_t amy;
    md_molecule_t ala;
    md_trajectory_i* ala_traj;
};

static md_molecule_t* amy = 0;
static md_molecule_t* ala = 0;
static md_trajectory_i* ala_traj = 0;

UTEST_F_SETUP(script) {
    md_vm_arena_init(&utest_fixture->arena, GIGABYTES(4));
    utest_fixture->alloc = md_vm_arena_create_interface(&utest_fixture->arena);

    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&utest_fixture->amy, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"),   NULL, &utest_fixture->alloc));
    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&utest_fixture->ala, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, &utest_fixture->alloc));

    md_util_molecule_postprocess(&utest_fixture->amy, &utest_fixture->alloc, MD_UTIL_POSTPROCESS_ALL);
    md_util_molecule_postprocess(&utest_fixture->ala, &utest_fixture->alloc, MD_UTIL_POSTPROCESS_ALL);

    utest_fixture->ala_traj = md_pdb_trajectory_create(STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), &utest_fixture->alloc, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
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

UTEST(script, type_equal) {
    {
        type_info_t a = {.base_type = TYPE_INT, .dim = {1}};
        type_info_t b = {.base_type = TYPE_INT, .dim = {1}};
        EXPECT_TRUE(type_info_equal(a,b));
    }
    {
        /*
        type_info_t a = {.base_type = TYPE_INT, .dim = {1,1}};
        type_info_t b = {.base_type = TYPE_INT, .dim = {1,0}};
        EXPECT_TRUE(type_info_equal(a,b));
        */
    }
    {
        type_info_t a = {.base_type = TYPE_INT, .dim = {1,0}};
        type_info_t b = {.base_type = TYPE_INT, .dim = {2,0}};
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

#define SETUP_EVAL_CTX(ir, mol, arena) \
	md_allocator_i temp_alloc = md_vm_arena_create_interface(&arena); \
	eval_context_t ctx = { \
		.ir = ir, \
		.mol = mol, \
		.temp_arena = &temp_arena, \
		.temp_alloc = &temp_alloc, \
		.alloc = &temp_alloc, \
	}

ast_node_t* parse_and_type_check_expression(str_t expr, md_script_ir_t* ir, md_molecule_t* mol, md_vm_arena_t* temp_arena) {
    // @HACK: We use alloc here: If the data type is a str_t, then it gets a shallow copy
    // Which means that the actual string data is contained within the ir->arena => temp_alloc
    ir->str = str_copy(expr, ir->arena);

    md_allocator_i temp_alloc = md_vm_arena_create_interface(temp_arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    bool result = false;

    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = &temp_alloc});
    node = prune_expressions(node);
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
    SETUP_TEMP_ALLOC(GIGABYTES(1));
    md_script_ir_t* ir = create_ir(&temp_alloc);

    {
        md_script_ir_clear(ir);

        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,2}"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* a = get_identifier(ir, STR_LIT("a"));
            ASSERT_NE(NULL, a);
            ASSERT_NE(NULL, a->data);
            EXPECT_EQ(TYPE_INT, a->data->type.base_type);
            EXPECT_EQ(1, a->data->type.dim[0]);

            identifier_t* b = get_identifier(ir, STR_LIT("b"));
            ASSERT_NE(NULL, b);
            ASSERT_NE(NULL, b->data);
            EXPECT_EQ(TYPE_INT, b->data->type.base_type);
            EXPECT_EQ(1, b->data->type.dim[0]);
        }
    }

    {
        md_script_ir_clear(ir);

        //@NOTE: Implicit conversion will kick in in the array composition converting 1 to 1.0f since that is the common 'compatible' type in array arguments
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,distance(1,2)}"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* a = get_identifier(ir, STR_LIT("a"));
            ASSERT_NE(NULL, a);
            ASSERT_NE(NULL, a->data);
            EXPECT_EQ(TYPE_FLOAT, a->data->type.base_type);
            EXPECT_EQ(1, a->data->type.dim[0]);

            identifier_t* b = get_identifier(ir, STR_LIT("b"));
            ASSERT_NE(NULL, b);
            ASSERT_NE(NULL, b->data);
            EXPECT_EQ(TYPE_FLOAT, b->data->type.base_type);
            EXPECT_EQ(1, b->data->type.dim[0]);
        }
    }

    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{x,y,z} = coord(1)"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* x = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, x);
            ASSERT_NE(NULL, x->data);
            EXPECT_EQ(TYPE_FLOAT, x->data->type.base_type);
            EXPECT_EQ(1, x->data->type.dim[0]);

            identifier_t* y = get_identifier(ir, STR_LIT("y"));
            ASSERT_NE(NULL, y);
            ASSERT_NE(NULL, y->data);
            EXPECT_EQ(TYPE_FLOAT, y->data->type.base_type);
            EXPECT_EQ(1, y->data->type.dim[0]);

            identifier_t* z = get_identifier(ir, STR_LIT("z"));
            ASSERT_NE(NULL, z);
            ASSERT_NE(NULL, z->data);
            EXPECT_EQ(TYPE_FLOAT, z->data->type.base_type);
            EXPECT_EQ(1, z->data->type.dim[0]);
        }
    }

    {
        md_script_ir_clear(ir);
        // @NOTE: LHS contains 2 arguments, RHS contains 3 => No match in assignment
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{x,y} = coord(1)"), ir, &test_mol, &vm_arena);
        EXPECT_FALSE(node);
    }

    FREE_TEMP_ALLOC();
}

UTEST(script, array) {
    SETUP_TEMP_ALLOC(GIGABYTES(1));
    md_script_ir_t* ir = create_ir(&temp_alloc);

    {
        str_t src = STR_LIT(
            "s1 = residue(1:2);\n"
            "s2 = residue(2:4);\n"
            "s = {s1, s2};"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL, NULL);
        EXPECT_TRUE(result);
        identifier_t* s1 = get_identifier(ir, STR_LIT("s1"));
        identifier_t* s2 = get_identifier(ir, STR_LIT("s2"));
        identifier_t* s  = get_identifier(ir, STR_LIT("s"));

        if (s1) {
            EXPECT_EQ(s1->data->type.base_type, TYPE_BITFIELD);
            EXPECT_EQ(s1->data->type.dim[0], 2);
        }
        if (s2) {
            EXPECT_EQ(s2->data->type.base_type, TYPE_BITFIELD);
            EXPECT_EQ(s2->data->type.dim[0], 3);
        }
        if (s) {
            EXPECT_EQ(s->data->type.base_type, TYPE_BITFIELD);
            EXPECT_EQ(s->data->type.dim[0], 5);
        }
    }

    {
        str_t src = STR_LIT(
            "d1 = distance(1,2) in residue(1:2);\n"
            "d2 = distance(1,2) in residue(2:4);\n"
            "d = {d1, d2};"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL, NULL);
        EXPECT_TRUE(result);
        identifier_t* d1 = get_identifier(ir, STR_LIT("d1"));
        identifier_t* d2 = get_identifier(ir, STR_LIT("d2"));
        identifier_t* d  = get_identifier(ir, STR_LIT("d"));

        if (d1) {
            EXPECT_EQ(d1->data->type.base_type, TYPE_FLOAT);
            EXPECT_EQ(d1->data->type.dim[0], 2);
        }
        if (d2) {
            EXPECT_EQ(d2->data->type.base_type, TYPE_FLOAT);
            EXPECT_EQ(d2->data->type.dim[0], 3);
        }
        if (d) {
            EXPECT_EQ(d->data->type.base_type,  TYPE_FLOAT);
            EXPECT_EQ(d->data->type.dim[0], 5);
        }
    }

    {
        str_t src = STR_LIT(
            "c1 = coord(1:2);\n"
            "c2 = coord(4:6);\n"
            "c = {c1, c2};"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL, NULL);
        EXPECT_TRUE(result);
        identifier_t* c1 = get_identifier(ir, STR_LIT("c1"));
        identifier_t* c2 = get_identifier(ir, STR_LIT("c2"));
        identifier_t* c  = get_identifier(ir, STR_LIT("c"));

        if (c1) {
            EXPECT_EQ(c1->data->type.base_type, TYPE_FLOAT);
            EXPECT_EQ(c1->data->type.dim[0], 2);
            EXPECT_EQ(c1->data->type.dim[1], 3);
        }
        if (c2) {
            EXPECT_EQ(c2->data->type.base_type, TYPE_FLOAT);
            EXPECT_EQ(c2->data->type.dim[0], 3);
            EXPECT_EQ(c1->data->type.dim[1], 3);
        }
        if (c) {
            EXPECT_EQ(c->data->type.base_type,  TYPE_FLOAT);
            EXPECT_EQ(c->data->type.dim[0], 5);
            EXPECT_EQ(c->data->type.dim[1], 3);
        }
    }

    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = coord(1) in residue(:)"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(4, ident->data->type.dim[0]);
            EXPECT_EQ(3, ident->data->type.dim[1]);
        }
    }
    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = angle(1,2,3) in residue(:)"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(4, ident->data->type.dim[0]);
        }
    }
    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = {angle(1,2,3), angle(3,2,1)} in residue(:)"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(4, ident->data->type.dim[0]);
            EXPECT_EQ(2, ident->data->type.dim[1]);
        }
    }
    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = {angle(1,2,3) in residue(:), angle(3,2,1) in residue(:)}"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(8, ident->data->type.dim[0]);
        }
    }

    FREE_TEMP_ALLOC();
}

UTEST(script, array_subscript) {
    SETUP_TEMP_ALLOC(GIGABYTES(1));
    md_script_ir_t* ir = create_ir(&temp_alloc);

    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[1]"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(3, ident->data->type.dim[0]);
            EXPECT_EQ(0, ident->data->type.dim[1]);
        }
    }
    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[1,:]"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(3, ident->data->type.dim[0]);
            EXPECT_EQ(0, ident->data->type.dim[1]);
        }
    }
    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[:,1]"), ir, &test_mol, &vm_arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(4, ident->data->type.dim[0]);
            EXPECT_EQ(0, ident->data->type.dim[1]);
        }
    }
    {
        str_t src = STR_LIT(
            "xyz  = coord(:);"
            "xyz1 = xyz[1,:];"
            "xyz2 = xyz[2,2:];"
            "xyz3 = xyz[2];"
            "x    = xyz[:,1];"
            "y    = xyz[:,2];"
            "z    = xyz[:,3];"
        );
        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL, NULL);

        eval_context_t ctx = {
            .ir = ir,
            .mol = &test_mol,
            .temp_arena = &vm_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
        };

        identifier_t* xyz = get_identifier(ir, STR_LIT("xyz"));
        EXPECT_TRUE(xyz);
        if (xyz) {
            EXPECT_EQ(16, xyz->data->type.dim[0]);
            EXPECT_EQ(3,  xyz->data->type.dim[1]);
        }

        identifier_t* xyz1 = get_identifier(ir, STR_LIT("xyz1"));
        EXPECT_TRUE(xyz1);
        if (xyz1) {
            EXPECT_TRUE(xyz1->data);
            EXPECT_EQ(3, xyz1->data->type.dim[0]);
            EXPECT_EQ(0, xyz1->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, xyz1->data->type, &temp_alloc);
            evaluate_node(&data, xyz1->node, &ctx);
            const float* coord = (const float*)data.ptr;
            EXPECT_NEAR(mol_x[0], coord[0], 1.0e-6f);
            EXPECT_NEAR(mol_y[0], coord[1], 1.0e-6f);
            EXPECT_NEAR(mol_z[0], coord[2], 1.0e-6f);
        }

        identifier_t* xyz2 = get_identifier(ir, STR_LIT("xyz2"));
        EXPECT_TRUE(xyz2);
        if (xyz2) {
            EXPECT_TRUE(xyz2->data);
            EXPECT_EQ(2, xyz2->data->type.dim[0]);
            EXPECT_EQ(0, xyz2->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, xyz2->data->type, &temp_alloc);
            evaluate_node(&data, xyz2->node, &ctx);
            const float* coord = (const float*)data.ptr;
            EXPECT_NEAR(mol_y[1], coord[0], 1.0e-6f);
            EXPECT_NEAR(mol_z[1], coord[1], 1.0e-6f);
        }

        identifier_t* xyz3 = get_identifier(ir, STR_LIT("xyz3"));
        EXPECT_TRUE(xyz3);
        if (xyz3) {
            EXPECT_TRUE(xyz1->data);
            EXPECT_EQ(3, xyz3->data->type.dim[0]);
            EXPECT_EQ(0, xyz3->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, xyz3->data->type, &temp_alloc);
            evaluate_node(&data, xyz3->node, &ctx);
            const float* coord = (const float*)data.ptr;
            EXPECT_NEAR(mol_x[1], coord[0], 1.0e-6f);
            EXPECT_NEAR(mol_y[1], coord[1], 1.0e-6f);
            EXPECT_NEAR(mol_z[1], coord[2], 1.0e-6f);
        }

        identifier_t* x = get_identifier(ir, STR_LIT("x"));
        EXPECT_TRUE(x);
        if (x) {
            EXPECT_TRUE(x->data);
            EXPECT_EQ(16, x->data->type.dim[0]);
            EXPECT_EQ(0,  x->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, x->data->type, &temp_alloc);
            evaluate_node(&data, x->node, &ctx);
            const float* coord = (const float*)data.ptr;
            for (int i = 0; i < ATOM_COUNT; ++i) {
                EXPECT_NEAR(mol_x[i], coord[i], 1.0e-6f);
            }
        }

        identifier_t* y = get_identifier(ir, STR_LIT("y"));
        EXPECT_TRUE(y);
        if (y) {
            EXPECT_TRUE(y->data);
            EXPECT_EQ(16, y->data->type.dim[0]);
            EXPECT_EQ(0,  y->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, y->data->type, &temp_alloc);
            evaluate_node(&data, y->node, &ctx);
            const float* coord = (const float*)data.ptr;
            for (int i = 0; i < ATOM_COUNT; ++i) {
                EXPECT_NEAR(mol_y[i], coord[i], 1.0e-6f);
            }
        }

        identifier_t* z = get_identifier(ir, STR_LIT("z"));
        EXPECT_TRUE(z);
        if (z) {
            EXPECT_TRUE(z->data);
            EXPECT_EQ(16, z->data->type.dim[0]);
            EXPECT_EQ(0,  z->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, z->data->type, &temp_alloc);
            evaluate_node(&data, z->node, &ctx);
            const float* coord = (const float*)data.ptr;
            for (int i = 0; i < ATOM_COUNT; ++i) {
                EXPECT_NEAR(mol_z[i], coord[i], 1.0e-6f);
            }
        }
    }

    FREE_TEMP_ALLOC();
}

UTEST(script, dim_op) {
    SETUP_TEMP_ALLOC(GIGABYTES(1));
    md_script_ir_t* ir = create_ir(&temp_alloc);

    {
        str_t src = STR_LIT(
            "xyz  = coord(1:15);\n"
            "flat = flatten(xyz);\n"
            "xyz_t = transpose(xyz);\n"
            "sw = shape_weights(residue(1:4));\n"
            "{lin,plan,iso} = transpose(shape_weights(residue(1:4)));"
        );
        md_script_ir_compile_from_source(ir, src, &test_mol, NULL, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        eval_context_t ctx = {
            .ir = ir,
            .mol = &test_mol,
            .temp_arena = &vm_arena,
            .temp_alloc = &temp_alloc,
            .alloc = &temp_alloc,
        };

        {
            identifier_t* xyz = get_identifier(ir, STR_LIT("xyz"));
            EXPECT_TRUE(xyz);
            if (xyz) {
                EXPECT_EQ(TYPE_FLOAT, xyz->data->type.base_type);
                EXPECT_EQ(15, xyz->data->type.dim[0]);
                EXPECT_EQ(3,  xyz->data->type.dim[1]);
            }

            data_t data = {0};
            allocate_data(&data, xyz->data->type, &temp_alloc);
            evaluate_node(&data, xyz->node, &ctx);
            const vec3_t* coord = (const vec3_t*)data.ptr;
            for (int i = 0; i < xyz->data->type.dim[0]; ++i) {
                EXPECT_NEAR(mol_x[i], coord[i].x, 1.0e-6f);
                EXPECT_NEAR(mol_y[i], coord[i].y, 1.0e-6f);
                EXPECT_NEAR(mol_z[i], coord[i].z, 1.0e-6f);
            }
        }

        {
            identifier_t* flat = get_identifier(ir, STR_LIT("flat"));
            EXPECT_TRUE(flat);
            if (flat) {
                EXPECT_EQ(TYPE_FLOAT, flat->data->type.base_type);
                EXPECT_EQ(45, flat->data->type.dim[0]);
                EXPECT_EQ(0,  flat->data->type.dim[1]);
            }

            data_t data = {0};
            allocate_data(&data, flat->data->type, &temp_alloc);
            evaluate_node(&data, flat->node, &ctx);
            const float* coord = (const float*)data.ptr;
            for (int i = 0; i < flat->data->type.dim[0]; ++i) {
                switch(i % 3) {
                case 0: EXPECT_NEAR(mol_x[i / 3], coord[i], 1.0e-6f); break;
                case 1: EXPECT_NEAR(mol_y[i / 3], coord[i], 1.0e-6f); break;
                case 2: EXPECT_NEAR(mol_z[i / 3], coord[i], 1.0e-6f); break;
                }
            }
        }

        {
            identifier_t* xyz_t = get_identifier(ir, STR_LIT("xyz_t"));
            EXPECT_TRUE(xyz_t);
            if (xyz_t) {
                EXPECT_EQ(TYPE_FLOAT, xyz_t->data->type.base_type);
                EXPECT_EQ(3,  xyz_t->data->type.dim[0]);
                EXPECT_EQ(15, xyz_t->data->type.dim[1]);

                data_t data = {0};
                allocate_data(&data, xyz_t->data->type, &temp_alloc);
                evaluate_node(&data, xyz_t->node, &ctx);
                const float* coords = (const float*)data.ptr;
                const float* x = coords + 0 * xyz_t->data->type.dim[1];
                const float* y = coords + 1 * xyz_t->data->type.dim[1];
                const float* z = coords + 2 * xyz_t->data->type.dim[1];
                for (int i = 0; i < xyz_t->data->type.dim[1]; ++i) {
                    EXPECT_NEAR(mol_x[i], x[i], 1.0e-6f);
                    EXPECT_NEAR(mol_y[i], y[i], 1.0e-6f);
                    EXPECT_NEAR(mol_z[i], z[i], 1.0e-6f);
                }
            }
        }

        {
            identifier_t* sw = get_identifier(ir, STR_LIT("sw"));
            EXPECT_TRUE(sw);
            if (sw) {
                EXPECT_EQ(TYPE_FLOAT, sw->data->type.base_type);
                EXPECT_EQ(4, sw->data->type.dim[0]);
                EXPECT_EQ(3, sw->data->type.dim[1]);
            }
        }

        {
            identifier_t* lin  = get_identifier(ir, STR_LIT("lin"));
            identifier_t* plan = get_identifier(ir, STR_LIT("plan"));
            identifier_t* iso  = get_identifier(ir, STR_LIT("iso"));

            EXPECT_TRUE(lin);
            EXPECT_TRUE(plan);
            EXPECT_TRUE(iso);

            if (lin) {
                EXPECT_EQ(TYPE_FLOAT, lin->data->type.base_type);
                EXPECT_EQ(4, lin->data->type.dim[0]);
                EXPECT_EQ(0, lin->data->type.dim[1]);
            }
            if (plan) {
                EXPECT_EQ(TYPE_FLOAT, plan->data->type.base_type);
                EXPECT_EQ(4, plan->data->type.dim[0]);
                EXPECT_EQ(0, plan->data->type.dim[1]);
            }
            if (iso) {
                EXPECT_EQ(TYPE_FLOAT, iso->data->type.base_type);
                EXPECT_EQ(4, iso->data->type.dim[0]);
                EXPECT_EQ(0, iso->data->type.dim[1]);
            }
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

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, implicit_conversion) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->amy;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sel = residue({1,2,3,4});"), mol, NULL, NULL));

    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("v = sdf(chain(:), resname(\"PFT\"), 50);"), mol, NULL, NULL));

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

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("{lin, plan, iso} = shape_weights(:);");
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        identifier_t* lin  = get_identifier(ir, STR_LIT("lin"));
        identifier_t* plan = get_identifier(ir, STR_LIT("plan"));
        identifier_t* iso  = get_identifier(ir, STR_LIT("iso"));

        EXPECT_TRUE(lin);
        if (lin) {
            EXPECT_EQ(TYPE_FLOAT, lin->data->type.base_type);
            EXPECT_EQ(1, lin->data->type.dim[0]);
        }

        EXPECT_TRUE(plan);
        if (plan) {
            EXPECT_EQ(TYPE_FLOAT, plan->data->type.base_type);
            EXPECT_EQ(1, plan->data->type.dim[0]);
        }

        EXPECT_TRUE(iso);
        if (iso) {
            EXPECT_EQ(TYPE_FLOAT, iso->data->type.base_type);
            EXPECT_EQ(1, iso->data->type.dim[0]);
        }
        EXPECT_EQ(md_script_ir_property_count(ir), 3);
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        if (eval) {
            EXPECT_TRUE(md_script_eval_property_data(eval, STR_LIT("lin")));
            EXPECT_TRUE(md_script_eval_property_data(eval, STR_LIT("plan")));
            EXPECT_TRUE(md_script_eval_property_data(eval, STR_LIT("iso")));
        }
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));
        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("prop1 = distance_pair(com(resname(\"ALA\")), 1);");
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        ASSERT_TRUE(md_script_ir_valid(ir));

        identifier_t* prop1 = get_identifier(ir, STR_LIT("prop1"));
        ASSERT_TRUE(prop1);
        EXPECT_EQ(TYPE_FLOAT, prop1->data->type.base_type);
        EXPECT_EQ(1, prop1->data->type.dim[0]);

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("d1 = 1:5 in residue(1:3);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("V = sdf(residue(1), element('H'), 5.0) * 2;"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("prop1 = rdf(element('C'), element('O'), 20.0);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("sel = within_x(0:100);\np1  = distance(com(sel), 100);"), mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, traj, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        str_t src = STR_LIT("s1 = count(within(10, residue(:)));");

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, src, mol, traj, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
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
    const size_t num_props = md_script_ir_property_count(data->ir);
    const str_t* props = md_script_ir_property_names(data->ir);
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(data->traj);
    if (md_script_eval_frame_range(data->eval, data->ir, data->mol, data->traj, 0, num_frames)) {
        for (size_t p_idx = 0; p_idx < num_props; ++p_idx) {
            const md_script_property_data_t* cur_data = md_script_eval_property_data(data->eval,     props[p_idx]);
            const md_script_property_data_t* ref_data = md_script_eval_property_data(data->ref_eval, props[p_idx]);
            ASSERT(cur_data->num_values == ref_data->num_values);
            for (size_t i = 0; i < cur_data->num_values; ++i) {
                if (cur_data->values[i] != ref_data->values[i]) {
                    data->num_corrupt_values += 1;
                    fprintf(stderr, "Corruption occured in thread %"PRIu64" at frame %i, expected: '%g', got: '%g'\n", md_thread_id(), (int)i, ref_data->values[i], cur_data->values[i]);
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

    size_t num_frames = md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, script, mol, traj, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));
    ASSERT_EQ(1, md_script_ir_property_count(ir));

    md_script_eval_t* ref_eval = md_script_eval_create(num_frames, ir, alloc);
    EXPECT_EQ(num_frames, md_script_eval_frame_count(ref_eval));
    ASSERT_TRUE(md_script_eval_frame_range(ref_eval, ir, mol, traj, 0, (uint32_t)num_frames));

    const md_script_property_data_t* data = md_script_eval_property_data(ref_eval, STR_LIT("p1"));
    ASSERT_TRUE(data);
    ASSERT_EQ(num_frames, data->num_values);

#if 0
    printf("Ref Values:\n");
    const md_script_property_t* ref_prop = md_script_eval_properties(ref_eval);
    for (int64_t i = 0; i < ref_prop->data.num_values; ++i) {
        printf("[%lli]: %g\n", i, ref_prop->data.values[i]);
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

UTEST_F(script, visualize) {
    md_allocator_i* alloc = md_arena_allocator_create(&utest_fixture->alloc, MEGABYTES(1));
    md_molecule_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = STR_LIT(
        "w = shape_weights(residue(:));"
        "{x,y,z} = {w[:,1], w[:,2], w[:,3]};"
        "sx = w[5,1];"
        );
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, src, mol, NULL, NULL));

        identifier_t*  x = get_identifier(ir, STR_LIT("x"));
        identifier_t* sx = get_identifier(ir, STR_LIT("sx"));

        ASSERT_TRUE(x);
        ASSERT_TRUE(sx);

        md_script_vis_t vis = {0};
        md_script_vis_init(&vis, alloc);
        md_script_vis_ctx_t ctx = {
            .ir = ir,
            .mol = mol,
        };
        
        EXPECT_TRUE(md_script_vis_eval_payload(&vis, (const md_script_vis_payload_o*)x->node, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        EXPECT_EQ(mol->atom.count, md_bitfield_popcount(&vis.atom_mask));
        
        md_script_vis_clear(&vis);
        EXPECT_TRUE(md_script_vis_eval_payload(&vis, (const md_script_vis_payload_o*)sx->node, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        size_t res_atom_count = md_residue_atom_count(mol->residue, 4);
        size_t pop_count = md_bitfield_popcount(&vis.atom_mask);
        EXPECT_EQ(res_atom_count, pop_count);
    }

    md_arena_allocator_destroy(alloc);
}
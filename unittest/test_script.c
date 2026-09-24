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
#include <md_system.h>
#include <md_trajectory.h>
#include <md_gro.h>
#include <md_pdb.h>

// Test only instrumentation of the evaluator, see test_hook_proc_eval
#define MD_SCRIPT_TEST_HOOKS
#include <md_script.c>

// Create molecule for evaulation
#define ATOM_COUNT 16
static float mol_x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
static float mol_y[] = {4,3,2,1,4,3,2,1,1,1,1,1,1,1,1,1};
static float mol_z[] = {3,2,1,4,3,2,1,2,2,2,2,2,2,2,2,2};
static md_atom_type_idx_t mol_ti[] = {0, 4, 0, 1, 2, 4, 2, 4, 0, 4, 0, 1, 2, 3, 2, 4};

#define ATOM_TYPE_COUNT 5
static md_atomic_number_t atom_type_z[] = {1, 2, 6, 7, 8};
static md_label_t atom_type_id[] = {
    MAKE_LABEL("H1"),
    MAKE_LABEL("HE"),
    MAKE_LABEL("CA"),
    MAKE_LABEL("N"),
    MAKE_LABEL("O")
};
static float atom_type_mass[] = {
    1.008f,
    4.002602f,
    12.011f,
    14.007f,
    15.999f
};
static float atom_type_radius[] = {
    1.2f,
    1.4f,
    1.7f,
    1.55f,
    1.52f
};

#define COMP_COUNT 4
static md_label_t r_name[] = {MAKE_LABEL("SOL"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT"), MAKE_LABEL("PFT")};
static md_sequence_id_t  r_id[] = {1, 2, 3, 4};
static uint32_t    r_off[] = {0, 3, 8, 12, 16};

#define ENT_COUNT 3
static md_label_t    e_id[] = {MAKE_LABEL("WAT"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT")};
static md_flags_t e_flags[] = {MD_FLAG_WATER, MD_FLAG_AMINO_ACID, MD_FLAG_HETERO};
static str_t       e_desc[] = {MAKE_LABEL("water"), MAKE_LABEL("polypeptide"), MAKE_LABEL("non-polymer")};

#define INST_COUNT 4
static md_label_t        i_id[] = {MAKE_LABEL("A"), MAKE_LABEL("B"), MAKE_LABEL("C"), MAKE_LABEL("D")};
static uint32_t         i_off[] = {0, 3, 8, 12, 16};
static md_entity_idx_t i_eidx[] = {0, 1, 2, 2};

md_system_t test_mol = {
    // A hand built system still needs a reference state: it is what evaluation reads coordinates
    // from when there is no trajectory frame.
    .reference = {
        .num_atoms = ATOM_COUNT,
        .x = mol_x,
        .y = mol_y,
        .z = mol_z,
    },
    .atom = {
        .count = ATOM_COUNT,
        .type_idx = mol_ti,
        .type = {
            .count = ATOM_TYPE_COUNT,
            .name = atom_type_id,
            .z = atom_type_z,
            .mass = atom_type_mass,
            .radius = atom_type_radius,
        },
    },
    .component = {
        .count = COMP_COUNT,
        .name = r_name,
        .seq_id = r_id,
        .atom_offset = r_off
    },
    .instance = {
        .count = INST_COUNT,
        .id = i_id,
        .comp_offset = i_off,
    },
    .entity = {
        .count = ENT_COUNT,
        .id = e_id,
        .flags = e_flags
    },
};

struct script {
    bool initialized;
    md_allocator_i* arena;
    md_system_t amy;
    md_system_t ala;
    md_system_t npt;   // Triclinic (rhombic dodecahedron) cell - see the within_* tests below
};

UTEST_F_SETUP(script) {
    utest_fixture->arena = md_vm_arena_create(GIGABYTES(4));

    utest_fixture->amy.alloc = utest_fixture->arena;
    md_system_state_t amy_state = { .alloc = utest_fixture->arena };
    ASSERT_TRUE(md_gro_system_init_from_file(&utest_fixture->amy, &amy_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
    md_util_system_infer(&utest_fixture->amy, &amy_state, MD_UTIL_INFER_ALL);

    utest_fixture->ala.alloc = utest_fixture->arena;
    md_system_state_t ala_state = { .alloc = utest_fixture->arena };
    ASSERT_TRUE(md_pdb_system_init_from_file(&utest_fixture->ala, &ala_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    md_util_system_infer(&utest_fixture->ala, &ala_state, MD_UTIL_INFER_ALL);

    // A triclinic cell, not just a box. Every other system in this fixture is orthorhombic, which
    // is exactly the shape that hides mistakes in the periodic neighbour walk.
    utest_fixture->npt.alloc = utest_fixture->arena;
    md_system_state_t npt_state = { .alloc = utest_fixture->arena };
    ASSERT_TRUE(md_gro_system_init_from_file(&utest_fixture->npt, &npt_state, STR_LIT(MD_UNITTEST_DATA_DIR "/npt.gro")));
    md_util_system_infer(&utest_fixture->npt, &npt_state, MD_UTIL_INFER_ALL);
}

UTEST_F_TEARDOWN(script) {
    md_vm_arena_destroy(utest_fixture->arena);
}

static bool eval_selection(md_bitfield_t* bitfield, str_t expr, md_system_t* mol) {
    ASSERT(bitfield);
    ASSERT(mol);

    md_temp_scope_t temp = md_temp_begin_avoid(bitfield->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    bool success = false;

    data_t data = {0};
    if (eval_expression(&data, expr, mol, temp_arena)) {
        if (data.type.base_type == TYPE_BITFIELD) {
            md_bitfield_t* res = (md_bitfield_t*)data.ptr;
            const int64_t len = type_info_array_len(data.type);
            for (int64_t i = 0; i < len; ++i) {
                md_bitfield_or_inplace(bitfield, &res[i]);
            }
            success = true;
        }
    }
    md_temp_end(temp);
    return success;
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

static bool compare_selection(const char* exprA, const char* exprB) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    bool result = false;
    md_bitfield_t a = {0};
    md_bitfield_t b = {0};
    md_bitfield_init(&a, temp_arena);
    md_bitfield_init(&b, temp_arena);

    if (!eval_selection(&a, str_from_cstr(exprA), &test_mol)) {
        printf("Failed evaluation of expression: '%s'\n", exprA);
        goto done;
    }
    if (!eval_selection(&b, str_from_cstr(exprB), &test_mol)) {
        printf("Failed evaluation of expression: '%s'\n", exprB);
        goto done;
    }

    bool cmp_res = bit_cmp((uint64_t*)a.bits, (uint64_t*)b.bits, 0, ATOM_COUNT);
    if (!cmp_res) {
        printf("Got (A):\n");
        print_bits((uint64_t*)a.bits, ATOM_COUNT);
        printf("\nExpected (B)\n");
        print_bits((uint64_t*)b.bits, ATOM_COUNT);
        printf("\n");
        goto done;
    }
    result = true;
done:
    md_temp_end(temp);
    return result;
}

#define SETUP_EVAL_CTX(ir, mol, arena) \
	md_allocator_i temp_alloc = md_vm_arena_create_interface(&arena); \
	eval_context_t ctx = { \
		.ir = ir, \
		.mol = mol, \
		.state = &(mol)->reference, \
		.reference = &(mol)->reference, \
		.temp_arena = &temp_arena, \
		.temp_alloc = &temp_alloc, \
		.alloc = &temp_alloc, \
	}

ast_node_t* parse_and_type_check_expression(str_t expr, md_script_ir_t* ir, md_system_t* sys, md_allocator_i* arena) {
    // @HACK: We use alloc here: If the data type is a str_t, then it gets a shallow copy
    // Which means that the actual string data is contained within the ir->arena => temp_alloc
    ir->str = str_copy(expr, ir->arena);

    tokenizer_t tokenizer = tokenizer_init(ir->str);
    bool result = false;

    ast_node_t* node = parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = arena});
    node = prune_expressions(node);
    if (node) {
        eval_context_t ctx = {
            .ir = ir,
            .sys = sys,
            .cur_state = &sys->reference,
            .ref_state = &sys->reference,
            .temp_alloc = arena,
            .alloc = arena,
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

UTEST_F(script, common_subexpression_elimination) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);
    str_t src = STR_LIT(
        "x = element('H') in resname('ALA');\n"
        "y = element('O') in resname('ALA');\n"
    );

    bool result = md_script_ir_compile_from_source(ir, src, &utest_fixture->ala, NULL);
    ASSERT_TRUE(result);

    identifier_t* x = get_identifier(ir, STR_LIT("x"));
    identifier_t* y = get_identifier(ir, STR_LIT("y"));

    uint64_t hx = hash_node(x->node, 0);
    uint64_t hy = hash_node(y->node, 0);

    EXPECT_NE(hx, hy);

    md_temp_end(temp_scope);
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
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);
    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("'this is a string'"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_STRING);
            str_t str = as_string(data);
            EXPECT_STREQ(str.ptr, "this is a string");
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("2 + 5"), &test_mol, arena);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_INT);
            EXPECT_EQ(as_int(data), 7);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("2 + 5.0"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
            EXPECT_EQ(as_float(data), 7.0);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("{2,1} + {1,8}"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_INT);
            EXPECT_EQ(data.type.dim[0], 2);
            EXPECT_EQ(as_int_arr(data)[0], 3);
            EXPECT_EQ(as_int_arr(data)[1], 9);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("{2.5,1.5} + {1.0,8.0}"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
            EXPECT_EQ(data.type.dim[0], 2);
            EXPECT_EQ(as_float_arr(data)[0], 3.5f);
            EXPECT_EQ(as_float_arr(data)[1], 9.5f);
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("7 - 4 + 2 - 3"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_INT);
            EXPECT_EQ(data.type.dim[0], 1);
            int val = as_int(data);
            EXPECT_EQ(val, (7-4+2-3));
        }
    }

    {
        data_t data = {0};
        bool result = eval_expression(&data, STR_LIT("2.5 + 1.0 - 2.5 + 1.0"), &test_mol, arena);
        EXPECT_TRUE(result);
        if (result) {
            EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
            EXPECT_EQ(data.type.dim[0], 1);
            float val = as_float(data);
            EXPECT_EQ(val, (2.5f + 1.0f - 2.5f + 1.0f));
        }
    }
    md_temp_end(temp);
}

UTEST(script, type_compatability) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    {
        md_script_ir_clear(ir);

        ast_node_t* node = parse_and_type_check_expression(STR_LIT("v = normalize(coord(4));"), ir, &test_mol, arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* v = get_identifier(ir, STR_LIT("v"));
            ASSERT_NE(NULL, v);
            ASSERT_NE(NULL, v->data);
            EXPECT_EQ(TYPE_FLOAT, v->data->type.base_type);
        }
    }
    md_temp_end(temp_scope);
}

#if 0
UTEST(script, vector_operations) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);
    data_t data = {0};
    bool result = eval_expression(&data, STR_LIT("vec3(2.5, 1.5, 1.0) + vec3(1.0, 8.0, 2.0)"), &test_mol, arena);
    EXPECT_TRUE(result);
    if (result) {
        EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
        EXPECT_EQ(data.type.dim[0], 3);
        EXPECT_EQ(as_float_arr(data)[0], 3.5f);
        EXPECT_EQ(as_float_arr(data)[1], 9.5f);
        EXPECT_EQ(as_float_arr(data)[2], 3.0f);
    }
    md_temp_end(temp);
}

UTEST(script, array_type_compatability) {    
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);
    data_t data = {0};
    bool result = eval_expression(&data, STR_LIT("coord(1) + coord(4)"), &test_mol, arena);
    EXPECT_TRUE(result);
    if (result) {
        EXPECT_EQ(data.type.base_type, TYPE_FLOAT);
        EXPECT_EQ(data.type.dim[0], 3);
    }
    md_temp_end(temp);
}
#endif

UTEST(script, assignment) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    {
        md_script_ir_clear(ir);

        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,2}"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{a,b} = {1,distance(1,2)}"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{x,y,z} = coord(1)"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("{x,y} = coord(1)"), ir, &test_mol, arena);
        EXPECT_FALSE(node);
    }

    md_temp_end(temp_scope);
}

UTEST(script, array) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    {
        str_t src = STR_LIT(
            "s1 = residue(1:2);\n"
            "s2 = residue(2:4);\n"
            "s = {s1, s2};"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
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
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
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
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
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
        // Stride on selection ranges: residue(1:2:3) -> residues {1,3}
        str_t src = STR_LIT(
            "s = residue(1:2:3);"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
        EXPECT_TRUE(result);
        identifier_t* s = get_identifier(ir, STR_LIT("s"));
        EXPECT_TRUE(s);
        if (s) {
            EXPECT_EQ(s->data->type.base_type, TYPE_BITFIELD);
            EXPECT_EQ(s->data->type.dim[0], 2);
        }
    }

    {
        // Stride on index ranges for coordinate selection: coord(1:2:5) -> atoms {1,3,5}
        str_t src = STR_LIT(
            "c = coord(1:2:5);"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
        EXPECT_TRUE(result);
        identifier_t* c  = get_identifier(ir, STR_LIT("c"));
        EXPECT_TRUE(c);
        if (c) {
            EXPECT_EQ(c->data->type.base_type,  TYPE_FLOAT);
            EXPECT_EQ(c->data->type.dim[0], 3);
            EXPECT_EQ(c->data->type.dim[1], 3);

            // Evaluate and verify values match mol_x/y/z at indices 0,2,4
            eval_context_t ctx = {
                .ir = ir,
                .sys = &test_mol,
                .cur_state = &test_mol.reference,
                .ref_state = &test_mol.reference,
                .temp_alloc = arena,
                .alloc = arena,
            };
            data_t data = {0};
            allocate_data(&data, c->data->type, arena);
            evaluate_node(&data, c->node, &ctx);
            const float* coords = (const float*)data.ptr;
            for (int i = 0; i < 3; ++i) {
                const int idx = i * 3;
                const int src_idx = i * 2; // 0,2,4
                EXPECT_NEAR(mol_x[src_idx], coords[idx + 0], 1.0e-6f);
                EXPECT_NEAR(mol_y[src_idx], coords[idx + 1], 1.0e-6f);
                EXPECT_NEAR(mol_z[src_idx], coords[idx + 2], 1.0e-6f);
            }
        }
    }

    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = coord(1) in residue(:)"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = angle(1,2,3) in residue(:)"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = {angle(1,2,3), angle(3,2,1)} in residue(:)"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = {angle(1,2,3) in residue(:), angle(3,2,1) in residue(:)}"), ir, &test_mol, arena);
        EXPECT_TRUE(node != NULL);
        if (node) {
            identifier_t* ident = get_identifier(ir, STR_LIT("x"));
            ASSERT_NE(NULL, ident);
            ASSERT_NE(NULL, ident->data);
            EXPECT_EQ(8, ident->data->type.dim[0]);
        }
    }

    md_temp_end(temp_scope);
}

UTEST(script, array_subscript) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    {
        md_script_ir_clear(ir);
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[1]"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[1,:]"), ir, &test_mol, arena);
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
        ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = (coord(1) in residue(:))[:,1]"), ir, &test_mol, arena);
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
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);

        eval_context_t ctx = {
            .ir = ir,
            .sys = &test_mol,
            .cur_state = &test_mol.reference,
            .ref_state = &test_mol.reference,
            .temp_alloc = arena,
            .alloc = arena,
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
            allocate_data(&data, xyz1->data->type, arena);
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
            allocate_data(&data, xyz2->data->type, arena);
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
            allocate_data(&data, xyz3->data->type, arena);
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
            allocate_data(&data, x->data->type, arena);
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
            allocate_data(&data, y->data->type, arena);
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
            allocate_data(&data, z->data->type, arena);
            evaluate_node(&data, z->node, &ctx);
            const float* coord = (const float*)data.ptr;
            for (int i = 0; i < ATOM_COUNT; ++i) {
                EXPECT_NEAR(mol_z[i], coord[i], 1.0e-6f);
            }
        }
    }

    {
        // Strided subscript over first dimension: pick x-component for atoms {1,3,5,7}
        str_t src = STR_LIT(
            "xyz = coord(1:8);"
            "xs  = xyz[1:2:7,1];"
        );
        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
        EXPECT_TRUE(result);

        eval_context_t ctx = {
            .ir = ir,
            .sys = &test_mol,
            .cur_state = &test_mol.reference,
            .ref_state = &test_mol.reference,
            .temp_alloc = arena,
            .alloc = arena,
        };

        identifier_t* xs = get_identifier(ir, STR_LIT("xs"));
        EXPECT_TRUE(xs);
        if (xs) {
            EXPECT_TRUE(xs->data);
            EXPECT_EQ(4, xs->data->type.dim[0]);
            EXPECT_EQ(0, xs->data->type.dim[1]);

            data_t data = {0};
            allocate_data(&data, xs->data->type, arena);
            evaluate_node(&data, xs->node, &ctx);
            const float* x = (const float*)data.ptr;
            EXPECT_NEAR(mol_x[0], x[0], 1.0e-6f);
            EXPECT_NEAR(mol_x[2], x[1], 1.0e-6f);
            EXPECT_NEAR(mol_x[4], x[2], 1.0e-6f);
            EXPECT_NEAR(mol_x[6], x[3], 1.0e-6f);
        }
    }

    md_temp_end(temp_scope);
}

UTEST(script, dim_op) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    {
        str_t src = STR_LIT(
            "xyz  = coord(1:15);\n"
            "flat = flatten(xyz);\n"
            "xyz_t = transpose(xyz);\n"
            "sw = shape_weights(residue(1:4));\n"
            "{lin,plan,iso} = transpose(shape_weights(residue(1:4)));"
        );
        md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        eval_context_t ctx = {
            .ir = ir,
            .sys = &test_mol,
            .cur_state = &test_mol.reference,
            .ref_state = &test_mol.reference,
            .temp_alloc = arena,
            .alloc = arena,
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
            allocate_data(&data, xyz->data->type, arena);
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
            allocate_data(&data, flat->data->type, arena);
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
                allocate_data(&data, xyz_t->data->type, arena);
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

    md_temp_end(temp_scope);
}

static bool test_selection(const char* expr, const char* ref_bit_str) {
    bool result = false;
    uint64_t ref = make_bits(ref_bit_str);
    md_bitfield_t bf = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);
    md_bitfield_init(&bf, arena);

    if (!eval_selection(&bf, str_from_cstr(expr), &test_mol)) {
        printf("Failed evaluation of expression: '%s'\n", expr);
        goto done;
    }
    bool cmp_res = bit_cmp((uint64_t*)bf.bits, &ref, 0, ATOM_COUNT);
    if (!cmp_res) {
        printf("Got:\n");
        print_bits((uint64_t*)bf.bits, ATOM_COUNT);
        printf("\nExpected\n");
        print_bits(&ref, ATOM_COUNT);
        printf("\n");
        goto done;
    }
    result = true;
done:
    md_temp_end(temp);
    return result;
}

UTEST(script, selection) {
    EXPECT_TRUE(test_selection("all",               "1111111111111111"));
    EXPECT_TRUE(test_selection("resname('SOL')",    "1110000000000000"));
    EXPECT_TRUE(test_selection("resname('LYS')",    "0001111100000000"));
    EXPECT_TRUE(test_selection("element('C')",      "0000101000001010"));
    EXPECT_TRUE(test_selection("label('CA')",       "0000101000001010"));
    EXPECT_TRUE(test_selection("atom(1) in resname('PFT')", "0000000010001000"));
    //TEST_SELECTION("atom(1:2) or element('O') in residue(:)", "1001000010001000");
}

UTEST(script, stride_procedures_selectors) {
    // Stride over residues (components) equals explicit union of picks
    EXPECT_TRUE(compare_selection("residue(1:2:4)", "residue(1) or residue(3)"));

    // Stride over instances equals explicit union
    EXPECT_TRUE(compare_selection("instance(1:2:4)", "instance(1) or instance(3)"));

    // Note: test dataset has a single chain; skip chain stride validation here.

    // Stride over atoms
    EXPECT_TRUE(test_selection("atom(1:2:8)", "1010101000000000"));
    EXPECT_TRUE(test_selection("atom(2:3:16)", "0100100100100100"));

    // Stride over elements by atomic number ranges (odd Z up to 10 -> H(1), N(7))
    EXPECT_TRUE(test_selection("element(1:2:10)", "1010000010100100"));
}

UTEST_F(script, compile_script) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    str_t script_src = load_textfile(STR_LIT(MD_UNITTEST_DATA_DIR "/script.txt"), alloc);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, script_src, &utest_fixture->amy, NULL));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, implicit_conversion) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->amy;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sel = residue({1,2,3,4});"), mol, NULL));

    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("v = sdf(chain(:), resname(\"PFT\"), 50);"), mol, NULL));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, semantic) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->amy;

    md_script_ir_t* ir = md_script_ir_create(alloc);

    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "dih1_pft = dihedral(22, 20, 1, 2) in resname('PFT');"
    ), mol, NULL));

    md_script_ir_clear(ir);
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("p1 = resname('ALA') resname('GLY');"), mol, NULL));
    
    md_script_ir_clear(ir);
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
    ), mol, NULL));

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, selection_big) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->amy;

    md_bitfield_t bf = md_bitfield_create(alloc);
    EXPECT_TRUE(eval_selection(&bf, STR_LIT("atom(1:20) and element('O') in chain(:)"), mol));

    md_arena_allocator_destroy(alloc);
}

// ================================================================================================
// within() with a search radius that reaches or exceeds the periodic cell.
//
// Reported: "set the range argument of within(...) in a selection in a representation larger than
// the cutoff in my trajectories and viamd crashes immediately with a segfault".
//
// within() builds an md_spatial_acc grid sized from the radius and then walks a +-ncell neighbourhood
// of cells around every query point. Two things change character once the radius approaches the box:
// cell_dim is CLAMPed to at least 1, so an axis collapses to a single cell and its cell size becomes
// the whole box rather than the requested extent; and the neighbourhood then has to reach across more
// than one periodic image, which the wrap logic may or may not be able to express.
//
// These sweep the whole radius range - well inside the box, around half the box, past the box, and
// absurd - against an orthorhombic box, a very elongated box, and a triclinic one.

// Ascending, and deliberately straddling the smallest box dimension of every system in the fixture.
static const double WITHIN_RADII[] = { 2.0, 10.0, 24.0, 47.0, 60.0, 120.0, 700.0, 5000.0 };

// within(r, atom(1)) is every atom within r of atom 1, minus atom 1 itself (within excludes its own
// input). Two properties have to hold for any r, and neither depends on the geometry:
//   - monotonicity: a larger radius can only ever select more atoms
//   - saturation:   once r exceeds the cell diagonal, every other atom is inside it
// Monotonicity is the interesting one. A guard that bails out and returns an empty set on a large
// radius does not crash, but it does silently answer "nothing is nearby" - and this catches that too.
// NOTE: utest's EXPECT_* macros write through a variable named utest_result that UTEST_F declares
// in the test body, so a shared helper has to take it along explicitly.
static void within_radius_sweep(int* utest_result, md_system_t* mol, const char* tag, md_allocator_i* alloc) {
    size_t prev_count = 0;
    for (size_t i = 0; i < ARRAY_SIZE(WITHIN_RADII); ++i) {
        char expr[128];
        snprintf(expr, sizeof(expr), "within(%.1f, atom(1))", WITHIN_RADII[i]);

        md_bitfield_t bf = md_bitfield_create(alloc);
        const bool ok = eval_selection(&bf, str_from_cstr(expr), mol);
        EXPECT_TRUE(ok);
        if (!ok) {
            printf("  [%s] r=%.1f: evaluation failed outright\n", tag, WITHIN_RADII[i]);
            md_bitfield_free(&bf);
            continue;
        }

        const size_t count = md_bitfield_popcount(&bf);
        if (count < prev_count) {
            printf("  [%s] r=%.1f selected %zu atoms, but r=%.1f already selected %zu\n",
                   tag, WITHIN_RADII[i], count, WITHIN_RADII[i - 1], prev_count);
        }
        EXPECT_GE(count, prev_count);
        prev_count = count;
        md_bitfield_free(&bf);
    }

    // The last radius in the sweep dwarfs every cell in the fixture, so nothing can be outside it.
    if (prev_count != mol->atom.count - 1) {
        printf("  [%s] at r=%.1f expected all %zu atoms bar the query atom, got %zu\n",
               tag, WITHIN_RADII[ARRAY_SIZE(WITHIN_RADII) - 1], mol->atom.count - 1, prev_count);
    }
    EXPECT_EQ(mol->atom.count - 1, prev_count);
}

UTEST_F(script, within_radius_sweep_ortho) {
    // 1ALA-560ns: orthorhombic, 46.6 x 96.7 x 48.4 A. Radii past ~46 A exceed the shortest axis.
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    within_radius_sweep(utest_result, &utest_fixture->ala, "1ALA ortho 47x97x48", alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, within_radius_sweep_elongated) {
    // centered.gro: 216 x 216 x 642 A. Very anisotropic, so the cell grid collapses along x and y
    // long before it does along z.
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    within_radius_sweep(utest_result, &utest_fixture->amy, "centered ortho 216x216x642", alloc);
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, within_radius_sweep_triclinic) {
    // npt.gro: rhombic dodecahedron, a = b = 47.8 A, c = 33.8 A, third vector leaning half a cell in
    // x and y. The triclinic neighbour walk is a separate code path from the orthorhombic one.
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    within_radius_sweep(utest_result, &utest_fixture->npt, "npt triclinic dodecahedron", alloc);
    md_arena_allocator_destroy(alloc);
}

// The frange overload takes a different path into the same machinery, and it is the one the report
// names ("the range argument"). An inner bound of zero has to agree with the scalar overload.
UTEST_F(script, within_radius_sweep_frange) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;

    for (size_t i = 0; i < ARRAY_SIZE(WITHIN_RADII); ++i) {
        char expr_rng[128];
        char expr_flt[128];
        snprintf(expr_rng, sizeof(expr_rng), "within(0:%.1f, atom(1))", WITHIN_RADII[i]);
        snprintf(expr_flt, sizeof(expr_flt), "within(%.1f, atom(1))",   WITHIN_RADII[i]);

        md_bitfield_t bf_rng = md_bitfield_create(alloc);
        md_bitfield_t bf_flt = md_bitfield_create(alloc);
        EXPECT_TRUE(eval_selection(&bf_rng, str_from_cstr(expr_rng), mol));
        EXPECT_TRUE(eval_selection(&bf_flt, str_from_cstr(expr_flt), mol));

        const size_t n_rng = md_bitfield_popcount(&bf_rng);
        const size_t n_flt = md_bitfield_popcount(&bf_flt);
        if (n_rng != n_flt) {
            printf("  r=%.1f: within(0:r) selected %zu, within(r) selected %zu\n", WITHIN_RADII[i], n_rng, n_flt);
        }
        EXPECT_EQ(n_flt, n_rng);

        md_bitfield_free(&bf_rng);
        md_bitfield_free(&bf_flt);
    }
    md_arena_allocator_destroy(alloc);
}

// Closest to what the report actually did: a large radius evaluated per frame over a trajectory,
// which is what a representation does. The spatial acc is rebuilt from each frame's coordinates and
// unit cell, so a frame dependent failure only shows up here.
UTEST_F(script, within_radius_over_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);
    ASSERT_GT(num_frames, 0u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    for (size_t i = 0; i < ARRAY_SIZE(WITHIN_RADII); ++i) {
        char src_buf[160];
        snprintf(src_buf, sizeof(src_buf), "n%zu = count(within(%.1f, atom(1)));", i, WITHIN_RADII[i]);

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, str_from_cstr(src_buf), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        if (!md_script_ir_valid(ir)) continue;

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));
        md_script_eval_free(eval);
    }
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, dynamic_length) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = STR_LIT("sel1 = residue(within_z(1:50));");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, property_compute) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);

    md_script_ir_t* ir = md_script_ir_create(alloc);

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("num = count(residue(resname('ALA') and within(3.0, protein)));");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("{lin, plan, iso} = shape_weights(:);");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
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
            const md_attributes_t* attributes = md_script_eval_attributes(eval);
            EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/lin")));
            EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/plan")));
            EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/iso")));
        }
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));
        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_LIT("prop1 = distance_pair(com(resname(\"ALA\")), 1);");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        ASSERT_TRUE(md_script_ir_valid(ir));

        identifier_t* prop1 = get_identifier(ir, STR_LIT("prop1"));
        ASSERT_TRUE(prop1);
        EXPECT_EQ(TYPE_FLOAT, prop1->data->type.base_type);
        EXPECT_EQ(1, prop1->data->type.dim[0]);

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("d1 = 1:5 in residue(1:3);"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        // Strided integer range constant in an expression
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("d2 = 1:2:5 in residue(1:3);"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("V = sdf(residue(1), element('H'), 5.0) * 2;"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("prop1 = rdf(element('C'), element('O'), 20.0);"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("sel = within_x(0:100);\np1  = distance(com(sel), 100);"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        str_t src = STR_LIT("s1 = count(within(10, residue(:)));");

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));
        md_script_eval_free(eval);
    }

    md_arena_allocator_destroy(alloc);
}

// Every property an evaluation computes is also published as an attribute under 'script/'. What
// this checks is that the attribute IS the property - the same storage and the shape the property
// kind implies - rather than a second copy of it that could drift.
UTEST_F(script, property_attributes) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);

    md_script_ir_t* ir = md_script_ir_create(alloc);

    // A temporal property with a single value per frame.
    {
        md_script_ir_clear(ir);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("d = distance(1, 2);"), mol, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);

        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        ASSERT_NE(NULL, attributes);
        EXPECT_EQ(num_frames, attributes->num_frames);

        const md_attribute_t* attr = md_attributes_find(attributes, STR_LIT("script/d"));
        ASSERT_TRUE(attr != NULL);
        EXPECT_EQ(MD_ATTRIBUTE_FLAG_TEMPORAL, attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
        EXPECT_EQ(MD_ATTRIBUTE_TYPE_F32, attr->format.type);
        EXPECT_EQ(1u, attr->format.components);
        EXPECT_EQ(2u, attr->format.rank);
        EXPECT_EQ(num_frames, attr->format.shape[0]);
        EXPECT_EQ(1u, attr->format.shape[1]);

        const md_attribute_t* range = md_attributes_find(attributes, STR_LIT("script/d/range"));
        ASSERT_TRUE(range != NULL);
        EXPECT_EQ(0u, range->format.rank);
        EXPECT_EQ(2u, range->format.components);
        EXPECT_TRUE(md_unit_equal(attr->unit, range->unit));

        // Nothing to summarise over: there is no population axis worth the name.
        EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/d/mean")) == NULL);

        const uint64_t version = md_attributes_version(attributes, attr->id);
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));
        EXPECT_TRUE(md_attributes_version(attributes, attr->id) > version);

        // One frame out of the middle, through the attribute rather than the pointer.
        float value = 0;
        md_attribute_slice_t slice = md_attribute_slice_1(num_frames / 2);
        EXPECT_EQ(1, md_attribute_extract_slice_f32(&value, 1, attr, &slice, md_unit_none()));
        EXPECT_EQ(((const float*)attr->data)[num_frames / 2], value);

        // The range covers every value that was evaluated.
        float minmax[2] = {0, 0};
        EXPECT_EQ(2, md_attribute_extract_f32(minmax, 2, range, md_unit_none()));
        for (uint32_t i = 0; i < num_frames; ++i) {
            const float v = ((const float*)attr->data)[i];
            EXPECT_TRUE(minmax[0] <= v && v <= minmax[1]);
        }

        md_script_eval_free(eval);
    }

    // A temporal property with a population: the per frame summary over it is published beside it.
    {
        md_script_ir_clear(ir);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("g = distance_pair(residue(:), 1);"), mol, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);

        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        ASSERT_TRUE(attributes != NULL);

        const md_attribute_t* attr = md_attributes_find(attributes, STR_LIT("script/g"));
        ASSERT_TRUE(attr != NULL);
        EXPECT_EQ(2u, attr->format.rank);
        EXPECT_EQ(num_frames, attr->format.shape[0]);
        ASSERT_TRUE(attr->format.shape[1] > 1);

        const md_attribute_t* mean = md_attributes_find(attributes, STR_LIT("script/g/mean"));
        const md_attribute_t* var  = md_attributes_find(attributes, STR_LIT("script/g/variance"));
        const md_attribute_t* ext  = md_attributes_find(attributes, STR_LIT("script/g/extent"));
        ASSERT_TRUE(mean != NULL);
        ASSERT_TRUE(var  != NULL);
        ASSERT_TRUE(ext  != NULL);

        EXPECT_EQ(1u, mean->format.rank);
        EXPECT_EQ(num_frames, mean->format.shape[0]);
        EXPECT_EQ(MD_ATTRIBUTE_FLAG_TEMPORAL, mean->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);

        // min and max are two COMPONENTS of one value, not two positions along an axis.
        EXPECT_EQ(2u, ext->format.components);
        EXPECT_EQ(1u, ext->format.rank);
        EXPECT_EQ(num_frames, ext->format.shape[0]);

        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        float minmax[2] = {0, 0};
        md_attribute_slice_t slice = md_attribute_slice_1(0);
        EXPECT_EQ(2, md_attribute_extract_slice_f32(minmax, 2, ext, &slice, md_unit_none()));
        EXPECT_TRUE(minmax[0] <= minmax[1]);

        // The per frame summary agrees with the population it summarises.
        const uint32_t pop = attr->format.shape[1];
        const float* values = (const float*)attr->data;
        float lo = FLT_MAX, hi = -FLT_MAX;
        for (uint32_t k = 0; k < pop; ++k) {
            lo = MIN(lo, values[k]);
            hi = MAX(hi, values[k]);
        }
        EXPECT_EQ(lo, minmax[0]);
        EXPECT_EQ(hi, minmax[1]);

        md_script_eval_free(eval);
    }

    // A distribution: values, per bin weights, and the bin coordinates as a sibling.
    {
        md_script_ir_clear(ir);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("h = rdf(element('C'), element('O'), 20.0);"), mol, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);

        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        ASSERT_TRUE(attributes != NULL);

        const md_attribute_t* attr   = md_attributes_find(attributes, STR_LIT("script/h"));
        const md_attribute_t* weight = md_attributes_find(attributes, STR_LIT("script/h/weight"));
        const md_attribute_t* bin    = md_attributes_find(attributes, STR_LIT("script/h/bin"));
        const md_attribute_t* range  = md_attributes_find(attributes, STR_LIT("script/h/range"));
        ASSERT_TRUE(attr   != NULL);
        ASSERT_TRUE(weight != NULL);
        ASSERT_TRUE(bin    != NULL);
        ASSERT_TRUE(range  != NULL);

        const uint32_t num_bins = attr->format.shape[0];
        ASSERT_TRUE(num_bins > 1);
        // The range is over the bin axis, so it is in the bins' unit rather than the values'.
        EXPECT_TRUE(md_unit_equal(bin->unit, range->unit));

        // A distribution is not indexed by frame, so it carries no temporal claim.
        EXPECT_EQ(0, attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
        EXPECT_EQ(1u, attr->format.rank);
        EXPECT_EQ(num_bins, attr->format.shape[0]);
        EXPECT_EQ(num_bins, weight->format.shape[0]);
        EXPECT_EQ(num_bins, bin->format.shape[0]);

        // The values and the weights are separate buffers; they used to be one.
        EXPECT_TRUE(attr->data != weight->data);

        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, 0, num_frames));

        // The bin coordinates are computed from the range the evaluation settled on.
        md_array(float) coords = md_array_create(float, num_bins, alloc);
        EXPECT_EQ(num_bins, md_attribute_extract_f32(coords, num_bins, bin, md_unit_none()));
        float x_range[2] = {0, 0};
        EXPECT_EQ(2, md_attribute_extract_f32(x_range, 2, range, md_unit_none()));
        const float x_min = x_range[0];
        const float x_max = x_range[1];
        EXPECT_NEAR(0.0f,  x_min, 1.0e-4f);
        EXPECT_NEAR(20.0f, x_max, 1.0e-4f);
        const float x_scl = (x_max - x_min) / (float)num_bins;
        EXPECT_NEAR(x_min + 0.5f * x_scl, coords[0], 1.0e-4f);
        EXPECT_NEAR(x_max - 0.5f * x_scl, coords[num_bins - 1], 1.0e-4f);
        for (uint32_t i = 1; i < num_bins; ++i) {
            ASSERT_TRUE(coords[i] > coords[i - 1]);
        }

        md_script_eval_free(eval);
    }

    // A volume: three index axes and nothing else to say about it.
    {
        md_script_ir_clear(ir);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("V = sdf(residue(1), element('H'), 5.0);"), mol, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);

        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        ASSERT_TRUE(attributes != NULL);

        const md_attribute_t* attr = md_attributes_find(attributes, STR_LIT("script/V"));
        ASSERT_TRUE(attr != NULL);
        EXPECT_EQ(3u, attr->format.rank);
        EXPECT_TRUE(attr->format.shape[0] > 1);
        EXPECT_EQ(attr->format.shape[0], attr->format.shape[1]);
        EXPECT_EQ(attr->format.shape[0], attr->format.shape[2]);
        EXPECT_EQ(0, attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
        EXPECT_TRUE(attr->data != NULL);
        // A volume has no range to speak of.
        EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/V/range")) == NULL);

        md_script_eval_free(eval);
    }
}

UTEST_F(script, identifier_names) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("first = 1;\nsecond = 2;"), &utest_fixture->ala, NULL));

    // Every identifier is listed under its own name
    const size_t num_ident = md_script_ir_num_identifiers(ir);
    const str_t* ident = md_script_ir_identifiers(ir);
    bool found[2] = {false, false};
    for (size_t i = 0; i < num_ident; ++i) {
        if (str_eq(ident[i], STR_LIT("first")))  found[0] = true;
        if (str_eq(ident[i], STR_LIT("second"))) found[1] = true;
    }
    EXPECT_TRUE(found[0]);
    EXPECT_TRUE(found[1]);
    md_script_ir_free(ir);
}

UTEST(script, stride_invalid_ranges) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);

    // Step zero should fail static check
    md_script_ir_clear(ir);
    ast_node_t* node = parse_and_type_check_expression(STR_LIT("x = residue(1:0:5)"), ir, &test_mol, arena);
    EXPECT_FALSE(node);

    // Descending range with positive step should fail (only ascending ranges supported)
    md_script_ir_clear(ir);
    node = parse_and_type_check_expression(STR_LIT("x = residue(5:2:1)"), ir, &test_mol, arena);
    EXPECT_FALSE(node);

    md_temp_end(temp_scope);
}

#define NUM_THREADS 8

typedef struct thread_data_t {
    const md_script_ir_t* ir;
    const md_system_t* sys;
    const md_script_eval_t* ref_eval;
    md_script_eval_t* eval;
    int num_corrupt_values;
} thread_data_t;

void func(void* user_data) {
    thread_data_t* data = (thread_data_t*)user_data;
    const size_t num_props = md_script_ir_property_count(data->ir);
    const str_t* props = md_script_ir_property_names(data->ir);
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(data->sys->trajectory);
    if (md_script_eval_frame_range(data->eval, data->ir, data->sys, 0, num_frames)) {
        for (size_t p_idx = 0; p_idx < num_props; ++p_idx) {
            char path[256];
            snprintf(path, sizeof(path), "script/%.*s", (int)props[p_idx].len, props[p_idx].ptr);
            const md_attribute_t* cur = md_attributes_find(md_script_eval_attributes(data->eval),     str_from_cstr(path));
            const md_attribute_t* ref = md_attributes_find(md_script_eval_attributes(data->ref_eval), str_from_cstr(path));
            ASSERT(cur && ref);
            const size_t num_values = md_attribute_element_count(&cur->format);
            ASSERT(num_values == md_attribute_element_count(&ref->format));
            const float* cur_values = (const float*)cur->data;
            const float* ref_values = (const float*)ref->data;
            for (size_t i = 0; i < num_values; ++i) {
                if (cur_values[i] != ref_values[i]) {
                    data->num_corrupt_values += 1;
                    fprintf(stderr, "Corruption occured in thread %"PRIu64" at frame %i, expected: '%g', got: '%g'\n", md_thread_id(), (int)i, ref_values[i], cur_values[i]);
                }
            }
        }
    }
}

UTEST_F(script, parallel_evaluation) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    md_trajectory_i* traj = utest_fixture->ala.trajectory;

    const str_t script = STR_LIT("p1 = distance(1,10);");

    md_script_eval_t* eval[NUM_THREADS] = {0};
    md_thread_t* threads[NUM_THREADS] = {0};
    thread_data_t thread_data[NUM_THREADS] = {0};

    size_t num_frames = md_trajectory_num_frames(traj);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, script, mol, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));
    ASSERT_EQ(1, md_script_ir_property_count(ir));

    md_script_eval_t* ref_eval = md_script_eval_create(num_frames, ir, alloc);
    EXPECT_EQ(num_frames, md_script_eval_frame_count(ref_eval));
    ASSERT_TRUE(md_script_eval_frame_range(ref_eval, ir, mol, 0, (uint32_t)num_frames));

    const md_attribute_t* p1 = md_attributes_find(md_script_eval_attributes(ref_eval), STR_LIT("script/p1"));
    ASSERT_TRUE(p1);
    ASSERT_EQ(num_frames, md_attribute_element_count(&p1->format));


    for (int pass = 0; pass < 10; ++pass) {
        for (int i = 0; i < NUM_THREADS; ++i) {
            eval[i] = md_script_eval_create(num_frames, ir, alloc);
            EXPECT_NE(NULL, eval[i]);
        }

        for (int i = 0; i < NUM_THREADS; ++i) {
            thread_data[i] = (thread_data_t) {
                .ir = ir,
                .sys = mol,
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
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sqrt(2) * -4;"), mol, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5-4;"), mol, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = -4;"), mol, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = (-4);"), mol, NULL));
        
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5 * (-4);"), mol, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = 5 * -4;"), mol, NULL));

        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("x = (5) - 4;"), mol, NULL));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, visualize) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = STR_LIT(
        "w = shape_weights(residue(:));"
        "{x,y,z} = {w[:,1], w[:,2], w[:,3]};"
        "sx = w[5,1];"
        );
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, src, mol, NULL));

        identifier_t*  x = get_identifier(ir, STR_LIT("x"));
        identifier_t* sx = get_identifier(ir, STR_LIT("sx"));

        ASSERT_TRUE(x);
        ASSERT_TRUE(sx);

        md_script_vis_t vis = {0};
        md_script_vis_init(&vis, alloc);
        md_script_vis_ctx_t ctx = {
            .ir = ir,
            .sys = mol,
            .state = &mol->reference,
        };
        
        EXPECT_TRUE(md_script_vis_eval_payload(&vis, (const md_script_vis_payload_o*)x->node, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        EXPECT_EQ(mol->atom.count, md_bitfield_popcount(&vis.atom_mask));
        
        md_script_vis_clear(&vis);
        EXPECT_TRUE(md_script_vis_eval_payload(&vis, (const md_script_vis_payload_o*)sx->node, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        size_t res_atom_count = md_component_atom_count(&mol->component, 4);
        size_t pop_count = md_bitfield_popcount(&vis.atom_mask);
        EXPECT_EQ(res_atom_count, pop_count);
    }

    md_arena_allocator_destroy(alloc);
}

// ### NAMED ARGUMENTS ###

UTEST(script, named_args_signature_table) {
    const str_t keywords[] = {
        STR_LIT("in"), STR_LIT("of"), STR_LIT("out"), STR_LIT("and"), STR_LIT("or"), STR_LIT("xor"), STR_LIT("not"),
    };

    for (size_t s = 0; s < ARRAY_SIZE(signatures); ++s) {
        const proc_sig_t* sig = &signatures[s];
        char msg[256];
        snprintf(msg, sizeof(msg), "signature of '%.*s'", STR_ARG(sig->proc));

        EXPECT_TRUE_MSG(sig->num_params > 0 && sig->num_params <= MAX_SUPPORTED_PROC_ARGS, msg);

        // One signature per procedure name
        for (size_t t = s + 1; t < ARRAY_SIZE(signatures); ++t) {
            EXPECT_FALSE_MSG(str_eq(sig->proc, signatures[t].proc), msg);
        }

        bool seen_omittable = false;     // A positional parameter which can be left out
        bool seen_optional = false;
        for (size_t p = 0; p < sig->num_params; ++p) {
            const param_sig_t* param = &sig->param[p];
            snprintf(msg, sizeof(msg), "parameter '%.*s' of '%.*s'", STR_ARG(param->name), STR_ARG(sig->proc));

            EXPECT_TRUE_MSG(md_script_identifier_name_valid(param->name), msg);
            for (size_t k = 0; k < ARRAY_SIZE(keywords); ++k) {
                EXPECT_FALSE_MSG(str_eq(param->name, keywords[k]), msg);
            }
            for (size_t q = p + 1; q < sig->num_params; ++q) {
                EXPECT_FALSE_MSG(str_eq(param->name, sig->param[q].name), msg);
            }

            const uint32_t omit = param->flags & (PARAM_OPTIONAL | PARAM_DEFAULT | PARAM_NULLABLE);
            const bool kw_only = param->flags & PARAM_KW_ONLY;
            // At most one way of being left out
            EXPECT_TRUE_MSG(omit == 0 || omit == PARAM_OPTIONAL || omit == PARAM_DEFAULT || omit == PARAM_NULLABLE, msg);
            // As in Python: a positional parameter which has to be given does not follow one which can be left out
            EXPECT_FALSE_MSG(!omit && !kw_only && seen_omittable, msg);
            // Optional parameters are trailing: leaving them out selects an overload taking fewer arguments
            EXPECT_FALSE_MSG(seen_optional && !(param->flags & PARAM_OPTIONAL), msg);
            if (param->flags & PARAM_DEFAULT) {
                EXPECT_TRUE_MSG(param->def_type.base_type != TYPE_UNDEFINED && is_scalar(param->def_type), msg);
            }
            seen_omittable |= omit && !kw_only;
            seen_optional  |= (param->flags & PARAM_OPTIONAL) != 0;
        }

        // Every overload fits the signature
        size_t num_overloads = 0;
        size_t min_arity = SIZE_MAX;
        for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
            const procedure_t* proc = &procedures[i];
            if (!str_eq(proc->name, sig->proc)) continue;
            num_overloads += 1;
            min_arity = MIN(min_arity, proc->num_args);

            snprintf(msg, sizeof(msg), "overload of '%.*s' taking %i arguments", STR_ARG(sig->proc), (int)proc->num_args);
            EXPECT_FALSE_MSG(proc->flags & FLAG_SYMMETRIC_ARGS, msg);
            EXPECT_LE_MSG(proc->num_args, sig->num_params, msg);
            // Anything an overload does not take must be possible to leave out
            for (size_t p = proc->num_args; p < sig->num_params; ++p) {
                EXPECT_TRUE_MSG(sig->param[p].flags & PARAM_OPTIONAL, msg);
            }
            // A nullable parameter is passed as absent, so every overload takes it
            for (size_t p = 0; p < sig->num_params; ++p) {
                if (sig->param[p].flags & PARAM_NULLABLE) {
                    EXPECT_LT_MSG(p, proc->num_args, msg);
                }
            }
        }
        snprintf(msg, sizeof(msg), "procedure '%.*s' of signature", STR_ARG(sig->proc));
        EXPECT_GT_MSG(num_overloads, (size_t)0, msg);

        // An optional parameter is only meaningful if some overload can be reached without it
        for (size_t p = 0; p < sig->num_params; ++p) {
            if (sig->param[p].flags & PARAM_OPTIONAL) {
                snprintf(msg, sizeof(msg), "optional parameter '%.*s' of '%.*s'", STR_ARG(sig->param[p].name), STR_ARG(sig->proc));
                EXPECT_LE_MSG(min_arity, p, msg);
            }
        }
    }
}

static bool named_args_same_data(const data_t* x, const data_t* y, md_allocator_i* alloc) {
    if (!type_info_equal(x->type, y->type)) return false;
    if (x->type.base_type == TYPE_BITFIELD) {
        const int64_t len = type_info_array_len(x->type);
        const md_bitfield_t* bx = (const md_bitfield_t*)x->ptr;
        const md_bitfield_t* by = (const md_bitfield_t*)y->ptr;
        for (int64_t i = 0; i < len; ++i) {
            const size_t pop = md_bitfield_popcount(&bx[i]);
            if (pop != md_bitfield_popcount(&by[i])) return false;
            md_bitfield_t both = md_bitfield_create(alloc);
            md_bitfield_and(&both, &bx[i], &by[i]);
            if (md_bitfield_popcount(&both) != pop) return false;
        }
        return true;
    }
    return x->size == y->size && memcmp(x->ptr, y->ptr, x->size) == 0;
}

// Evaluates a positional reference call and a variant of it with named arguments and requires identical results
static bool named_args_equivalent(const char* positional, const char* named, md_system_t* sys, md_allocator_i* alloc) {
    data_t a = {0};
    data_t b = {0};
    if (!eval_expression(&a, str_from_cstr(positional), sys, alloc)) {
        printf("Failed to evaluate '%s'\n", positional);
        return false;
    }
    if (!eval_expression(&b, str_from_cstr(named), sys, alloc)) {
        printf("Failed to evaluate '%s'\n", named);
        return false;
    }
    if (!named_args_same_data(&a, &b, alloc)) {
        printf("Results differ between '%s' and '%s'\n", positional, named);
        return false;
    }
    return true;
}

UTEST_F(script, named_args_equivalence) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;

    static const char* cases[][2] = {
        {"distance(1, 10)",                         "distance(a=1, b=10)"},
        {"distance(1, 10)",                         "distance(b=10, a=1)"},
        {"distance(1, 10)",                         "distance(1, b=10)"},
        {"distance(com(residue(1)), com(residue(3)))", "distance(b=com(residue(3)), a=com(residue(1)))"},
        {"distance_min(residue(1), residue(3))",    "distance_min(b=residue(3), a=residue(1))"},
        {"distance_max(residue(1), residue(3))",    "distance_max(b=residue(3), a=residue(1))"},
        // Order matters here: the result is laid out [a][b]
        {"distance_pair(1:2, 4:6)",                 "distance_pair(b=4:6, a=1:2)"},
        {"angle(1, 2, 3)",                          "angle(c=3, a=1, b=2)"},
        {"dihedral(1, 2, 3, 4)",                    "dihedral(1, 2, d=4, c=3)"},
        {"within(4.0, residue(1))",                 "within(around=residue(1), radius=4.0)"},
        {"within(2.0:4.0, residue(1))",             "within(radius=2.0:4.0, around=residue(1))"},
        {"within(4.0) in residue(1)",               "within(radius=4.0) in residue(1)"},
        {"within_xyz(0:10, 0:12, 0:14)",            "within_xyz(z=0:14, y=0:12, x=0:10)"},
        {"count(protein)",                          "count(sel=protein)"},
        {"count(protein, 'residue')",               "count(unit='residue', sel=protein)"},
        {"split(protein, 3)",                       "split(parts=3, sel=protein)"},
        {"rdf(element('C'), element('H'), 10.0)",   "rdf(element('C'), element('H'), cutoff=10.0)"},
        {"contact_count(residue(:), residue(:), 3.0)", "contact_count(b=residue(:), a=residue(:), cutoff=3.0)"},
        // Named calls nested in arrays and as arguments of other calls
        {"{distance(1, 10), distance(2, 10)}",      "{distance(a=1, b=10), distance(b=10, a=2)}"},
        {"count(within(4.0, residue(1)))",          "count(sel=within(around=residue(1), radius=4.0))"},
    };

    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        EXPECT_TRUE(named_args_equivalent(cases[i][0], cases[i][1], sys, alloc));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, named_args_compile) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* sys = &utest_fixture->ala;

    // The default script of viamd, and the same script with its arguments given by name in a different order.
    // Binding happens before anything else looks at the calls, so the two must compile to the same thing.
    str_t positional = STR_LIT(
        "s1 = resname(\"ALA\")[2:8];\n"
        "d1 = distance(10,30);\n"
        "a1 = angle(2,1,3) in resname(\"ALA\");\n"
        "r = rdf(element('C'), element('H'), 10.0);\n"
        "v = sdf(s1, element('H'), 10.0);\n"
        "{lin,plan,iso} = shape_weights(all);\n");
    str_t named = STR_LIT(
        "s1 = resname(\"ALA\")[2:8];\n"
        "d1 = distance(b=30, a=10);\n"
        "a1 = angle(c=3, a=2, b=1) in resname(\"ALA\");\n"
        "r = rdf(element('C'), element('H'), cutoff=10.0);\n"
        "v = sdf(target=element('H'), structures=s1, extent=10.0);\n"
        "{lin,plan,iso} = shape_weights(all);\n");

    md_script_ir_t* ir_pos = md_script_ir_create(alloc);
    md_script_ir_t* ir_named = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir_pos, positional, sys, NULL));
    ASSERT_TRUE(md_script_ir_compile_from_source(ir_named, named, sys, NULL));
    EXPECT_EQ(md_script_ir_fingerprint(ir_pos), md_script_ir_fingerprint(ir_named));

    // Whitespace around '=' and names which are also identifiers in the script
    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "a = residue(1);\n"
        "b = residue(3);\n"
        "d = distance_min( b = b , a = a );\n"), sys, NULL));

    md_arena_allocator_destroy(alloc);
}

static bool named_args_compile_fails_with(md_script_ir_t* ir, md_system_t* sys, const char* src, const char* expected) {
    md_script_ir_clear(ir);
    if (md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL)) {
        printf("Expected compilation of '%s' to fail\n", src);
        return false;
    }
    const size_t num_errors = md_script_ir_num_errors(ir);
    const md_log_token_t* errors = md_script_ir_errors(ir);
    for (size_t i = 0; i < num_errors; ++i) {
        char buf[512];
        snprintf(buf, sizeof(buf), "%.*s", STR_ARG(errors[i].text));
        if (strstr(buf, expected)) return true;
    }
    printf("Compilation of '%s' failed, but without the expected error '%s'. Got:\n", src, expected);
    for (size_t i = 0; i < num_errors; ++i) {
        printf("  %.*s\n", STR_ARG(errors[i].text));
    }
    return false;
}

UTEST_F(script, named_args_errors) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* sys = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(a=1, c=2);",         "has no parameter named 'c'"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(a=1, a=2);",         "given more than once"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(a=1, 2);",           "positional argument cannot follow"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(1, a=2);",           "already given by position"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(a=1);",              "Missing argument 'b'"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(1);",                "Missing argument 'b'"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = distance(1, 2, 3, b=4);",     "Too many arguments"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = within(around=residue(1));",  "Missing argument 'radius'"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = sqrt(x=2.0);",                "does not accept named arguments"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = flatten(a=residue(1:2));",    "does not accept named arguments"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = {a=1, 2};",                   "only valid in procedure calls"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = residue(1:3)[a=1];",          "only valid in procedure calls"));
    // A name that parses fine but binds to an argument of the wrong type still goes through overload resolution
    EXPECT_TRUE(named_args_compile_fails_with(ir, sys, "x = split(sel=protein, parts='three');", "Could not find matching procedure"));

    md_arena_allocator_destroy(alloc);
}

static ast_node_t* named_args_parse_only(md_script_ir_t* ir, const char* expr, md_allocator_i* alloc) {
    ir->str = str_copy(str_from_cstr(expr), ir->arena);
    tokenizer_t tokenizer = tokenizer_init(ir->str);
    return prune_expressions(parse_expression(&(parse_context_t){ .ir = ir, .tokenizer = &tokenizer, .temp_alloc = alloc }));
}

// The binder on its own, against signatures which exercise what the procedures ported so far do not:
// defaults, keyword only parameters, omitted optional parameters and parameter names shadowing procedures.
UTEST_F(script, named_args_binding) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_script_ir_t* ir = create_ir(alloc);
    eval_context_t ctx = {
        .ir = ir,
        .sys = &utest_fixture->ala,
        .temp_alloc = alloc,
        .alloc = alloc,
    };

    const proc_sig_t sig_default = {
        .proc = STR_LIT("angle"), .num_params = 3,
        .param = {
            {STR_LIT("a"), PARAM_REQUIRED},
            {STR_LIT("b"), PARAM_REQUIRED},
            {STR_LIT("c"), PARAM_DEFAULT, .def_type = TI_INT, .def = {._int = 7}},
        },
    };
    const proc_sig_t sig_optional = {
        .proc = STR_LIT("angle"), .num_params = 3,
        .param = {
            {STR_LIT("a"), PARAM_REQUIRED},
            {STR_LIT("b"), PARAM_OPTIONAL},
            {STR_LIT("c"), PARAM_OPTIONAL},
        },
    };
    const proc_sig_t sig_kw_only = {
        .proc = STR_LIT("angle"), .num_params = 3,
        .param = {
            {STR_LIT("a"), PARAM_REQUIRED},
            {STR_LIT("b"), PARAM_REQUIRED},
            {STR_LIT("c"), PARAM_KW_ONLY},
        },
    };
    const proc_sig_t sig_shadow = {
        .proc = STR_LIT("angle"), .num_params = 3,
        .param = {
            {STR_LIT("count"),   PARAM_REQUIRED},
            {STR_LIT("residue"), PARAM_REQUIRED},
            {STR_LIT("c"),       PARAM_REQUIRED},
        },
    };

    ast_node_t* node = 0;

    // Default supplied for a positional call
    node = named_args_parse_only(ir, "angle(1, 2)", alloc);
    ASSERT_TRUE(node && node->type == AST_PROC_CALL);
    ASSERT_TRUE(bind_arguments(node, &sig_default, &ctx));
    ASSERT_EQ(md_array_size(node->children), (size_t)3);
    EXPECT_EQ(node->children[2]->type, AST_CONSTANT_VALUE);
    EXPECT_EQ(node->children[2]->value._int, 7);

    // Default overridden by name, arguments reordered
    node = named_args_parse_only(ir, "angle(c=3, b=2, a=1)", alloc);
    ASSERT_TRUE(node);
    ASSERT_TRUE(node->named_args);
    ASSERT_TRUE(bind_arguments(node, &sig_default, &ctx));
    EXPECT_FALSE(node->named_args);
    ASSERT_EQ(md_array_size(node->children), (size_t)3);
    EXPECT_EQ(node->children[0]->value._int, 1);
    EXPECT_EQ(node->children[1]->value._int, 2);
    EXPECT_EQ(node->children[2]->value._int, 3);

    // Binding is idempotent
    ASSERT_TRUE(bind_arguments(node, &sig_default, &ctx));
    ASSERT_EQ(md_array_size(node->children), (size_t)3);
    EXPECT_EQ(node->children[2]->value._int, 3);

    // Trailing optional parameters may be left out, but not ones before a given argument
    node = named_args_parse_only(ir, "angle(a=1)", alloc);
    ASSERT_TRUE(bind_arguments(node, &sig_optional, &ctx));
    EXPECT_EQ(md_array_size(node->children), (size_t)1);
    node = named_args_parse_only(ir, "angle(1, c=3)", alloc);
    EXPECT_FALSE(bind_arguments(node, &sig_optional, &ctx));

    // Keyword only
    node = named_args_parse_only(ir, "angle(1, 2, 3)", alloc);
    EXPECT_FALSE(bind_arguments(node, &sig_kw_only, &ctx));
    node = named_args_parse_only(ir, "angle(1, 2, c=3)", alloc);
    EXPECT_TRUE(bind_arguments(node, &sig_kw_only, &ctx));

    // Parameter names which are also procedure names are recognised as names, not parsed as calls
    node = named_args_parse_only(ir, "angle(residue=2, c=3, count=1)", alloc);
    ASSERT_TRUE(node);
    ASSERT_TRUE(bind_arguments(node, &sig_shadow, &ctx));
    ASSERT_EQ(md_array_size(node->children), (size_t)3);
    EXPECT_EQ(node->children[0]->value._int, 1);
    EXPECT_EQ(node->children[1]->value._int, 2);
    EXPECT_EQ(node->children[2]->value._int, 3);

    // '==' is a comparison, not a named argument
    node = named_args_parse_only(ir, "count(x == 1)", alloc);
    EXPECT_TRUE(node == NULL || node->named_args == NULL);

    md_arena_allocator_destroy(alloc);
}

// ### IDENTIFIER MEMOIZATION ###
// A dynamic identifier is evaluated once per frame; references to it reuse that value.

// Evaluates src over the whole trajectory and returns the number of times proc_ptr was evaluated
static size_t memo_count_evaluations(md_script_eval_t** out_eval, str_t src, int (*proc_ptr)(data_t*, data_t[], eval_context_t*), md_system_t* sys, md_allocator_i* alloc) {
    md_script_ir_t* ir = md_script_ir_create(alloc);
    if (!md_script_ir_compile_from_source(ir, src, sys, NULL)) {
        printf("Failed to compile '%.*s'\n", STR_ARG(src));
        return SIZE_MAX;
    }
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);

    test_hook_proc_eval.proc_ptr = proc_ptr;
    test_hook_proc_eval.count = 0;
    const bool ok = md_script_eval_frame_range(eval, ir, sys, 0, (uint32_t)num_frames);
    const size_t count = test_hook_proc_eval.count;
    test_hook_proc_eval.proc_ptr = NULL;
    test_hook_proc_eval.count = 0;

    if (!ok) {
        printf("Failed to evaluate '%.*s'\n", STR_ARG(src));
        return SIZE_MAX;
    }
    if (out_eval) *out_eval = eval;
    return count;
}

// Requires the temporal property 'name' to hold the same values in both evaluations
static bool memo_same_property(const md_script_eval_t* a, const md_script_eval_t* b, str_t name_a, str_t name_b) {
    const md_attribute_t* pa = md_attributes_find(md_script_eval_attributes(a), name_a);
    const md_attribute_t* pb = md_attributes_find(md_script_eval_attributes(b), name_b);
    if (!pa || !pb) {
        printf("Missing property '%.*s' or '%.*s'\n", STR_ARG(name_a), STR_ARG(name_b));
        return false;
    }
    const size_t na = md_attribute_element_count(&pa->format);
    const size_t nb = md_attribute_element_count(&pb->format);
    if (na != nb || na == 0) {
        printf("Property sizes differ: %zu vs %zu\n", na, nb);
        return false;
    }
    const float* va = (const float*)pa->data;
    const float* vb = (const float*)pb->data;
    for (size_t i = 0; i < na; ++i) {
        if (va[i] != vb[i]) {
            printf("Property values differ at %zu: %g vs %g\n", i, va[i], vb[i]);
            return false;
        }
    }
    return true;
}

UTEST_F(script, identifier_memoization) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);
    ASSERT_GT(num_frames, (size_t)1);

    // Referenced three times, directly and through other identifiers: evaluated once per frame
    {
        md_script_eval_t* memo = 0;
        md_script_eval_t* inline_eval = 0;
        EXPECT_EQ(num_frames, memo_count_evaluations(&memo, STR_LIT(
            "d = distance(1, 10);"
            "e = d * 2;"
            "f = d + e;"), _distance, sys, alloc));

        // Same values as spelling the expressions out
        EXPECT_EQ(4 * num_frames, memo_count_evaluations(&inline_eval, STR_LIT(
            "d = distance(1, 10);"
            "e = distance(1, 10) * 2;"
            "f = distance(1, 10) + distance(1, 10) * 2;"), _distance, sys, alloc));
        ASSERT_TRUE(memo && inline_eval);
        EXPECT_TRUE(memo_same_property(memo, inline_eval, STR_LIT("script/e"), STR_LIT("script/e")));
        EXPECT_TRUE(memo_same_property(memo, inline_eval, STR_LIT("script/f"), STR_LIT("script/f")));
    }

    // Variable length: the length of r is only known per frame. Neither sizing nor evaluating the
    // references may run the spatial query again.
    {
        md_script_eval_t* memo = 0;
        md_script_eval_t* inline_eval = 0;
        EXPECT_EQ(num_frames, memo_count_evaluations(&memo, STR_LIT(
            "s = within(3.0, residue(1));"
            "r = residue(s);"
            "n = count(r);"
            "m = count(r, 'residue') + count(s);"), _within_expl_flt, sys, alloc));

        memo_count_evaluations(&inline_eval, STR_LIT(
            "n = count(residue(within(3.0, residue(1))));"
            "m = count(residue(within(3.0, residue(1))), 'residue') + count(within(3.0, residue(1)));"), _within_expl_flt, sys, alloc);
        ASSERT_TRUE(memo && inline_eval);
        EXPECT_TRUE(memo_same_property(memo, inline_eval, STR_LIT("script/n"), STR_LIT("script/n")));
        EXPECT_TRUE(memo_same_property(memo, inline_eval, STR_LIT("script/m"), STR_LIT("script/m")));
    }

    // Destructured identifiers are reused as well
    {
        EXPECT_EQ(num_frames, memo_count_evaluations(NULL, STR_LIT(
            "{x, y, z} = com(residue(1));"
            "a = x + y + z;"
            "b = x * y;"), _com, sys, alloc));
    }

    md_arena_allocator_destroy(alloc);
}

// Referenced within a context, the expression is evaluated per context and must not be replaced by the
// value computed outside of it.
UTEST_F(script, identifier_memoization_context) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    const bool compiled = md_script_ir_compile_from_source(ir, STR_LIT(
        "d = distance(1, 2);"
        "x = d in residue(1:3);"), sys, NULL);

    if (compiled) {
        md_script_eval_t* memo = 0;
        md_script_eval_t* inline_eval = 0;
        memo_count_evaluations(&memo, STR_LIT(
            "d = distance(1, 2);"
            "x = d in residue(1:3);"), _distance, sys, alloc);
        memo_count_evaluations(&inline_eval, STR_LIT(
            "d = distance(1, 2);"
            "x = distance(1, 2) in residue(1:3);"), _distance, sys, alloc);
        ASSERT_TRUE(memo && inline_eval);
        EXPECT_TRUE(memo_same_property(memo, inline_eval, STR_LIT("script/x"), STR_LIT("script/x")));
    }

    md_arena_allocator_destroy(alloc);
}


// Identifiers declared by destructuring used to be evaluated by evaluating the whole expression into the
// storage of a single element, which asserted. Now they are memoized, and where that cannot apply (within a
// context) the element is extracted from an evaluation of the whole expression.
UTEST_F(script, identifier_destructured_reference) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);

    static const char* srcs[] = {
        "{x, y, z} = com(residue(1)); a = x * 2;",
        "{x, y, z} = coord(1); a = x + y + z;",
        "{x, y, z} = com(residue(1)); a = x in residue(1:3);",
    };
    for (size_t i = 0; i < ARRAY_SIZE(srcs); ++i) {
        md_script_ir_t* ir = md_script_ir_create(alloc);
        EXPECT_TRUE_MSG(md_script_ir_compile_from_source(ir, str_from_cstr(srcs[i]), sys, NULL), srcs[i]);
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_TRUE_MSG(md_script_eval_frame_range(eval, ir, sys, 0, (uint32_t)num_frames), srcs[i]);
    }

    // The value of a destructured identifier is its element of the expression
    {
        md_script_eval_t* a = 0;
        md_script_eval_t* b = 0;
        memo_count_evaluations(&a, STR_LIT("{x, y, z} = com(residue(1)); v = y * 2;"), NULL, sys, alloc);
        memo_count_evaluations(&b, STR_LIT("c = com(residue(1)); v = c[2] * 2;"), NULL, sys, alloc);
        ASSERT_TRUE(a && b);
        EXPECT_TRUE(memo_same_property(a, b, STR_LIT("script/v"), STR_LIT("script/v")));
    }

    md_arena_allocator_destroy(alloc);
}


// ### VARIABLE LENGTH EVALUATION ###
// A value whose length is only known per frame is evaluated once: its arguments are evaluated a single time
// and the length follows from them. This used to take one evaluation to find the length and another to fill
// in the result, doubling the work for every level of nesting.

UTEST_F(script, variable_length_single_evaluation) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);

    // Each level of residue() has a length which depends on the frame
    md_script_eval_t* e1 = 0;
    md_script_eval_t* e2 = 0;
    md_script_eval_t* e3 = 0;
    EXPECT_EQ(num_frames, memo_count_evaluations(&e1, STR_LIT("n = count(residue(within(3.0, residue(1))));"), _within_expl_flt, sys, alloc));
    EXPECT_EQ(num_frames, memo_count_evaluations(&e2, STR_LIT("n = count(residue(residue(within(3.0, residue(1)))));"), _within_expl_flt, sys, alloc));
    EXPECT_EQ(num_frames, memo_count_evaluations(&e3, STR_LIT("n = count(residue(residue(residue(within(3.0, residue(1))))));"), _within_expl_flt, sys, alloc));
    ASSERT_TRUE(e1 && e2 && e3);
    // residue() of whole residues is the same residues
    EXPECT_TRUE(memo_same_property(e1, e2, STR_LIT("script/n"), STR_LIT("script/n")));
    EXPECT_TRUE(memo_same_property(e1, e3, STR_LIT("script/n"), STR_LIT("script/n")));

    // The atoms of the residues within 3 Å, counted two ways
    md_script_eval_t* by_residue = 0;
    md_script_eval_t* by_atom = 0;
    memo_count_evaluations(&by_residue, STR_LIT("n = count(residue(within(3.0, residue(1))));"), NULL, sys, alloc);
    memo_count_evaluations(&by_atom, STR_LIT("n = count(within(3.0, residue(1)) , 'residue');"), NULL, sys, alloc);
    ASSERT_TRUE(by_residue && by_atom);
    {
        const md_attribute_t* a = md_attributes_find(md_script_eval_attributes(by_atom), STR_LIT("script/n"));
        const md_attribute_t* r = md_attributes_find(md_script_eval_attributes(e1), STR_LIT("script/n"));
        ASSERT_TRUE(a && r);
        // count(x) counts atoms: the atoms of the residues is at least the number of residues
        const float* av = (const float*)a->data;
        const float* rv = (const float*)r->data;
        for (size_t i = 0; i < num_frames; ++i) {
            EXPECT_GE(rv[i], av[i]);
            EXPECT_GT(av[i], 0.0f);
        }
    }

    // A length deduced from an argument whose length is only known per frame
    {
        md_script_eval_t* a = 0;
        md_script_eval_t* b = 0;
        EXPECT_EQ(num_frames, memo_count_evaluations(&a, STR_LIT("m = max(coord_x(within(3.0, residue(1))) * 2);"), _within_expl_flt, sys, alloc));
        memo_count_evaluations(&b, STR_LIT("m = max(coord_x(within(3.0, residue(1)))) * 2;"), NULL, sys, alloc);
        ASSERT_TRUE(a && b);
        EXPECT_TRUE(memo_same_property(a, b, STR_LIT("script/m"), STR_LIT("script/m")));
    }

    // Arguments whose lengths are only known per frame, and which must match
    {
        md_script_eval_t* a = 0;
        md_script_eval_t* b = 0;
        memo_count_evaluations(&a, STR_LIT("s = within(3.0, residue(1)); m = max(coord_x(s) + coord_y(s));"), NULL, sys, alloc);
        memo_count_evaluations(&b, STR_LIT("m = max(coord_x(within(3.0, residue(1))) + coord_y(within(3.0, residue(1))));"), NULL, sys, alloc);
        ASSERT_TRUE(a && b);
        EXPECT_TRUE(memo_same_property(a, b, STR_LIT("script/m"), STR_LIT("script/m")));

        // Different lengths are an evaluation error, not a read past the end of the shorter one
        md_script_ir_t* ir = md_script_ir_create(alloc);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("m = max(coord_x(within(3.0, residue(1))) + coord_x(within(8.0, residue(1))));"), sys, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_FALSE(md_script_eval_frame_range(eval, ir, sys, 0, (uint32_t)num_frames));
    }

    // Through a declared identifier whose length is only known per frame
    {
        md_script_eval_t* a = 0;
        EXPECT_EQ(num_frames, memo_count_evaluations(&a, STR_LIT(
            "r = residue(within(3.0, residue(1)));"
            "n = count(residue(r));"), _within_expl_flt, sys, alloc));
        ASSERT_TRUE(a);
        EXPECT_TRUE(memo_same_property(a, e1, STR_LIT("script/n"), STR_LIT("script/n")));
    }

    md_arena_allocator_destroy(alloc);
}


// ### CONTACTS ###

#include <md_contact.h>

static bool contacts_compile_fails_with(md_system_t* sys, md_allocator_i* alloc, const char* src, const char* expected) {
    md_script_ir_t* ir = md_script_ir_create(alloc);
    return named_args_compile_fails_with(ir, sys, src, expected);
}

UTEST_F(script, contacts_compile) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* sys = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "c = contacts(residue(:), cutoff=4.5);"
        "n = count(c);"
        "na = count(c, 'atom');"
        "d = degree(c);"
        "c2 = contacts(residue(1:5), residue(6:15), cutoff=4.0, exclude_bonds=0);"
        "n2 = count(c2);"
        "c3 = contacts(chunks(residue(:), 5), cutoff=5, min_separation=3, parent=residue(:), exclude_within=chain(:));"
        "d3 = degree(c3);"), sys, NULL));

    // d is a temporal property with one value per residue
    EXPECT_EQ(md_script_ir_property_flags(ir, STR_LIT("d")), MD_SCRIPT_PROPERTY_FLAG_TEMPORAL);

    // The cutoff has to be named: after b, which may be left out, positions are ambiguous
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "c = contacts(residue(:), residue(:), 4.5);", "can only be given by name"));
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "c = contacts(residue(:));", "Missing argument 'cutoff'"));
    // Groups which change from frame to frame have no fixed pairs to speak of
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "c = contacts(residue(within(3.0, residue(1))), cutoff=4.0);", "cannot depend on the frame"));
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "c = contacts(residue(:), cutoff=0);", "has to be positive"));
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "c = contacts(residue(:), cutoff=4.0); n = count(c, 'bead');", "'group'"));
    EXPECT_TRUE(contacts_compile_fails_with(sys, alloc, "s = chunks(residue(:), 0);", "has to be positive"));

    md_arena_allocator_destroy(alloc);
}

// The same groups for md_contact, built without the script: one per residue
static md_bitfield_t* contacts_residue_groups(size_t* out_num, const md_system_t* sys, md_allocator_i* alloc) {
    md_bitfield_t* groups = md_alloc(alloc, sizeof(md_bitfield_t) * sys->component.count);
    for (size_t c = 0; c < sys->component.count; ++c) {
        const md_urange_t range = md_system_component_atom_range(sys, c);
        groups[c] = md_bitfield_create(alloc);
        md_bitfield_set_range(&groups[c], range.beg, range.end);
    }
    *out_num = sys->component.count;
    return groups;
}

// Per frame, what the script publishes against md_contact on the same frame
UTEST_F(script, contacts_matches_kernel) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);

    md_script_eval_t* eval = 0;
    // Evaluated once per frame, although four properties refer to it
    EXPECT_EQ(num_frames, memo_count_evaluations(&eval, STR_LIT(
        "c = contacts(residue(:), cutoff=6.0, min_separation=2);"
        "n = count(c);"
        "g = count(c, 'group');"
        "na = count(c, 'atom');"
        "d = degree(c);"), _contacts, sys, alloc));
    ASSERT_TRUE(eval);

    const md_attributes_t* attr = md_script_eval_attributes(eval);
    const md_attribute_t* n  = md_attributes_find(attr, STR_LIT("script/n"));
    const md_attribute_t* g  = md_attributes_find(attr, STR_LIT("script/g"));
    const md_attribute_t* na = md_attributes_find(attr, STR_LIT("script/na"));
    const md_attribute_t* d  = md_attributes_find(attr, STR_LIT("script/d"));
    ASSERT_TRUE(n && g && na && d);

    size_t num = 0;
    md_bitfield_t* res = contacts_residue_groups(&num, sys, alloc);
    ASSERT_EQ(md_attribute_element_count(&d->format), num_frames * num);

    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 6.0, .exclude_bonds = 3, .min_separation = 2 };
    md_contact_query_t q;
    ASSERT_TRUE(md_contact_query_init(&q, &desc, sys, alloc));
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_system_state_init(&state, sys->atom.count));

    float* deg = md_alloc(alloc, sizeof(float) * num);
    size_t total = 0;
    for (size_t f = 0; f < num_frames; ++f) {
        ASSERT_TRUE(md_trajectory_load_frame(sys->trajectory, f, &state));
        md_contact_set_t set;
        ASSERT_TRUE(md_contact_query_eval(&set, &q, &state, alloc));
        size_t atom_pairs = 0;
        memset(deg, 0, sizeof(float) * num);
        for (size_t k = 0; k < set.count; ++k) {
            atom_pairs += set.atom_pairs[k];
            deg[set.i[k]] += 1;
            deg[set.j[k]] += 1;
        }
        EXPECT_EQ((float)set.count, ((const float*)n->data)[f]);
        EXPECT_EQ((float)set.count, ((const float*)g->data)[f]);
        EXPECT_EQ((float)atom_pairs, ((const float*)na->data)[f]);
        for (size_t r = 0; r < num; ++r) {
            EXPECT_EQ(deg[r], ((const float*)d->data)[f * num + r]);
        }
        total += set.count;
        md_contact_set_free(&set);
    }
    EXPECT_GT(total, (size_t)0);
    md_contact_query_free(&q);
    md_arena_allocator_destroy(alloc);
}

// chunks: consecutive runs of the given size, per element, in index order
UTEST_F(script, chunks) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* sys = &utest_fixture->ala;

    // The evaluation's transient allocator is a temp arena, as in every evaluation
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    data_t data = {0};
    md_script_ir_t* ir = create_ir(alloc);
    ast_node_t* node = parse_and_type_check_expression(STR_LIT("chunks(residue(:), 7)"), ir, sys, temp_alloc);
    ASSERT_TRUE(node);
    eval_context_t ctx = { .ir = ir, .sys = sys, .temp_alloc = temp_alloc, .alloc = temp_alloc, .cur_state = &sys->reference, .ref_state = &sys->reference };
    ASSERT_TRUE(evaluate_node_alloc(&data, node, &ctx, temp_alloc));
    const md_bitfield_t* chunk = (const md_bitfield_t*)data.ptr;
    const size_t num_chunks = element_count(data);

    size_t expected = 0;
    for (size_t c = 0; c < sys->component.count; ++c) {
        const md_urange_t range = md_system_component_atom_range(sys, c);
        for (uint32_t beg = range.beg; beg < range.end; beg += 7) {
            const uint32_t end = MIN(beg + 7, range.end);
            ASSERT_LT(expected, num_chunks);
            EXPECT_EQ((size_t)(end - beg), md_bitfield_popcount(&chunk[expected]));
            EXPECT_EQ((size_t)(end - beg), md_bitfield_popcount_range(&chunk[expected], beg, end));
            expected += 1;
        }
    }
    EXPECT_EQ(expected, num_chunks);
    md_temp_end(temp);
    md_arena_allocator_destroy(alloc);
}

// parent and exclude_within against md_contact given the same parents and labels, derived here without the script
UTEST_F(script, contacts_parent_and_exclude_within) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;

    size_t num_res = 0;
    md_bitfield_t* res = contacts_residue_groups(&num_res, sys, alloc);

    // Chunks of 5 particles per residue, and for each the residue it is part of
    md_array(md_bitfield_t) chunk = 0;
    md_array(uint32_t) parent = 0;
    uint32_t* label = md_alloc(alloc, sizeof(uint32_t) * sys->atom.count);
    for (size_t c = 0; c < num_res; ++c) {
        const md_urange_t range = md_system_component_atom_range(sys, c);
        for (uint32_t k = range.beg; k < range.end; ++k) label[k] = (uint32_t)c;
        for (uint32_t beg = range.beg; beg < range.end; beg += 5) {
            md_bitfield_t bf = md_bitfield_create(alloc);
            md_bitfield_set_range(&bf, beg, MIN(beg + 5, range.end));
            md_array_push(chunk, bf, alloc);
            md_array_push(parent, (uint32_t)c, alloc);
        }
    }
    const size_t num_chunks = md_array_size(chunk);

    struct { const char* src; md_contact_desc_t desc; } cases[] = {
        { "c = contacts(chunks(residue(:), 5), cutoff=5.0, min_separation=3, parent=residue(:));",
          { .group_a = chunk, .num_a = num_chunks, .cutoff = 5.0, .exclude_bonds = 3, .min_separation = 3, .group_parent = parent } },
        { "c = contacts(chunks(residue(:), 5), cutoff=5.0, exclude_within=residue(:));",
          { .group_a = chunk, .num_a = num_chunks, .cutoff = 5.0, .exclude_bonds = 3, .particle_label = label } },
    };

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_script_ir_t* ir = create_ir(alloc);
        ast_node_t* node = parse_and_type_check_expression(str_from_cstr(cases[i].src), ir, sys, temp_alloc);
        ASSERT_TRUE(node);
        eval_context_t ctx = { .ir = ir, .sys = sys, .temp_alloc = temp_alloc, .alloc = temp_alloc, .cur_state = &sys->reference, .ref_state = &sys->reference };
        data_t data = {0};
        ASSERT_TRUE(evaluate_node_alloc(&data, node, &ctx, temp_alloc));
        const md_contact_set_t* got = (const md_contact_set_t*)data.ptr;

        md_contact_set_t want;
        ASSERT_TRUE(md_contact_compute(&want, &cases[i].desc, sys, &sys->reference, alloc));
        EXPECT_GT(want.count, (size_t)0);
        ASSERT_EQ(want.count, got->count);
        for (size_t k = 0; k < want.count; ++k) {
            EXPECT_EQ(want.i[k], got->i[k]);
            EXPECT_EQ(want.j[k], got->j[k]);
            EXPECT_EQ(want.atom_pairs[k], got->atom_pairs[k]);
        }
        // Nothing within one residue
        if (cases[i].desc.particle_label) {
            for (size_t k = 0; k < got->count; ++k) {
                EXPECT_NE(parent[got->i[k]], parent[got->j[k]]);
            }
        }
    }
    md_temp_end(temp);
    md_arena_allocator_destroy(alloc);
}

// Neither compiling nor visualizing searches for contacts: in a large system either would stall the application.
// Only the evaluation of frames does.
UTEST_F(script, contacts_no_search_outside_evaluation) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);

    test_hook_contact_searches = 0;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "c = contacts(residue(:), cutoff=6.0);"
        "d = degree(c);"
        "n = count(c, 'atom');"), sys, NULL));
    EXPECT_EQ((size_t)0, test_hook_contact_searches);
    // The shape was still known: one value per residue
    EXPECT_EQ(md_script_ir_property_flags(ir, STR_LIT("d")), MD_SCRIPT_PROPERTY_FLAG_TEMPORAL);

    // Hovering the expressions
    const str_t names[] = { STR_LIT("c"), STR_LIT("d"), STR_LIT("n") };
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        identifier_t* ident = get_identifier(ir, names[i]);
        ASSERT_TRUE(ident);
        md_script_vis_t vis = {0};
        md_script_vis_init(&vis, alloc);
        md_script_vis_ctx_t vctx = { .ir = ir, .sys = sys, .state = &sys->reference };
        EXPECT_TRUE(md_script_vis_eval_payload(&vis, (const md_script_vis_payload_o*)ident->node, -1, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        // What is shown are the groups
        EXPECT_EQ(sys->atom.count, md_bitfield_popcount(&vis.atom_mask));
        md_script_vis_free(&vis);
    }
    EXPECT_EQ((size_t)0, test_hook_contact_searches);

    // Evaluating the trajectory searches once per frame
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, sys, 0, (uint32_t)num_frames));
    EXPECT_EQ(num_frames, test_hook_contact_searches);

    md_arena_allocator_destroy(alloc);
}


// ### CONTACTS: THE SCRIPT AGAINST THE KERNEL ###
// The cases of test_contact.c, written as scripts. Each case states its groups twice: in C, exactly as the kernel
// test does, and as script. The script's groups are first checked to be the same groups, so both ask the same
// question. The contact set the script produces is then held to md_contact_compute given the kernel description
// (bit for bit), and to the brute force reference of the kernel tests (contact_reference.h). The script is
// evaluated the way a frame evaluation does it: compiled once, the contacts found through the query prepared for
// the call site. Criteria the script does not expose (radii, type pairs) and the particle pair layer have no
// counterpart here.

#include "contact_reference.h"

// One group per component whose flags intersect the mask (or every component if the mask is 0), as test_contact.c
static md_bitfield_t* mirror_residue_groups(size_t* out_count, const md_system_t* sys, md_flags_t mask, md_allocator_i* alloc) {
    md_bitfield_t* groups = md_alloc(alloc, sizeof(md_bitfield_t) * MAX(sys->component.count, 1));
    size_t n = 0;
    for (size_t c = 0; c < sys->component.count; ++c) {
        if (mask && !(md_system_component_flags(sys, c) & mask)) continue;
        const md_urange_t range = md_system_component_atom_range(sys, c);
        groups[n] = md_bitfield_create(alloc);
        md_bitfield_set_range(&groups[n], range.beg, range.end);
        ++n;
    }
    *out_count = n;
    return groups;
}

static bool mirror_bitfield_equal(const md_bitfield_t* a, const md_bitfield_t* b, md_allocator_i* alloc) {
    const size_t pa = md_bitfield_popcount(a);
    if (pa != md_bitfield_popcount(b)) return false;
    md_bitfield_t both = md_bitfield_create(alloc);
    md_bitfield_and(&both, a, b);
    return md_bitfield_popcount(&both) == pa;
}

// The value of an identifier of a compiled script, evaluated on a state as the evaluation of a frame does it
static bool mirror_eval(data_t* out, md_script_ir_t* ir, const char* name, md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    identifier_t* ident = get_identifier(ir, str_from_cstr(name));
    if (!ident || !ident->node) {
        printf("  no identifier '%s'\n", name);
        return false;
    }
    eval_context_t ctx = { .ir = ir, .sys = sys, .temp_alloc = alloc, .alloc = alloc, .cur_state = state, .ref_state = &sys->reference };
    return evaluate_node_alloc(out, ident->node, &ctx, alloc);
}

static bool mirror_same_groups(const data_t* d, const md_bitfield_t* groups, size_t num, const char* name, md_allocator_i* alloc) {
    if (d->type.base_type != TYPE_BITFIELD) {
        printf("  '%s' is not a selection\n", name);
        return false;
    }
    const size_t n = element_count(*d);
    if (n != num) {
        printf("  '%s' has %zu groups, the kernel case %zu\n", name, n, num);
        return false;
    }
    const md_bitfield_t* bf = as_bitfield(*d);
    for (size_t k = 0; k < n; ++k) {
        if (!mirror_bitfield_equal(&bf[k], &groups[k], alloc)) {
            printf("  group %zu of '%s' differs from the kernel case (%zu vs %zu particles)\n", k, name, md_bitfield_popcount(&bf[k]), md_bitfield_popcount(&groups[k]));
            return false;
        }
    }
    return true;
}

static bool mirror_same_sets(const md_contact_set_t* s, const md_contact_set_t* k) {
    if (s->count != k->count || s->num_a != k->num_a || s->num_b != k->num_b || s->flags != k->flags) {
        printf("  script: %zu pairs of %u x %u groups (flags %u), kernel: %zu pairs of %u x %u groups (flags %u)\n",
            s->count, s->num_a, s->num_b, s->flags, k->count, k->num_a, k->num_b, k->flags);
        return false;
    }
    for (size_t x = 0; x < s->count; ++x) {
        if (s->i[x] != k->i[x] || s->j[x] != k->j[x] || s->atom_pairs[x] != k->atom_pairs[x] || s->d_min[x] != k->d_min[x]) {
            printf("  pair %zu: script (%u, %u) %u pairs d_min %f, kernel (%u, %u) %u pairs d_min %f\n", x,
                s->i[x], s->j[x], s->atom_pairs[x], s->d_min[x], k->i[x], k->j[x], k->atom_pairs[x], k->d_min[x]);
            return false;
        }
    }
    return true;
}

// src declares the groups as 'ga' (and 'gb' if desc has a B set) and the contacts as 'c'. The groups of desc are
// those of the kernel case. On success the script's set is returned in out, if given.
static bool mirror_check(md_system_t* sys, const char* src, const md_contact_desc_t* desc, bool expect_contacts, md_contact_set_t* out, md_allocator_i* alloc) {
    md_script_ir_t* ir = md_script_ir_create(alloc);
    if (!md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL)) {
        printf("  failed to compile '%s'\n", src);
        for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
            printf("  " STR_FMT "\n", STR_ARG(md_script_ir_errors(ir)[i].text));
        }
        return false;
    }
    const md_system_state_t* state = &sys->reference;

    data_t ga = {0}, gb = {0}, c = {0};
    if (!mirror_eval(&ga, ir, "ga", sys, state, alloc) || !mirror_same_groups(&ga, desc->group_a, desc->num_a, "ga", alloc)) return false;
    if (desc->group_b) {
        if (!mirror_eval(&gb, ir, "gb", sys, state, alloc) || !mirror_same_groups(&gb, desc->group_b, desc->num_b, "gb", alloc)) return false;
    }

    test_hook_contact_local_prepares = 0;
    if (!mirror_eval(&c, ir, "c", sys, state, alloc) || c.type.base_type != TYPE_CONTACT) {
        printf("  failed to evaluate the contacts of '%s'\n", src);
        return false;
    }
    if (test_hook_contact_local_prepares != 0) {
        printf("  the contacts were not found through the query prepared for the call site\n");
        return false;
    }
    const md_contact_set_t* set = (const md_contact_set_t*)c.ptr;

    md_contact_set_t kernel = {0};
    if (!md_contact_compute(&kernel, desc, sys, state, alloc)) {
        printf("  md_contact_compute failed\n");
        return false;
    }
    bool ok = mirror_same_sets(set, &kernel);

    md_array(ref_pair_t) ref = ref_contacts(desc, sys, state, alloc);
    ok = same_contacts(set, ref, md_array_size(ref)) && ok;
    if (expect_contacts && set->count == 0) {
        printf("  no contacts, the case is vacuous\n");
        ok = false;
    }
    if (out) *out = *set;
    return ok;
}

// "atom(b:e)" of the particles of components [beg, end), one based and inclusive
static int mirror_atom_range(char* buf, size_t cap, const md_system_t* sys, size_t beg, size_t end) {
    const uint32_t a = md_system_component_atom_range(sys, beg).beg;
    const uint32_t b = md_system_component_atom_range(sys, end - 1).end;
    return snprintf(buf, cap, "atom(%u:%u)", a + 1, b);
}

UTEST_F(script, contacts_mirror_self_distance) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* res = mirror_residue_groups(&num, sys, 0, alloc);
    ASSERT_GT(num, (size_t)10);
    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 3.5 };
    EXPECT_TRUE(mirror_check(sys, "ga = residue(:); c = contacts(ga, cutoff=3.5, exclude_bonds=0);", &desc, true, NULL, alloc));
    // Written inline, as one would
    EXPECT_TRUE(mirror_check(sys, "ga = residue(:); c = contacts(residue(:), cutoff=3.5, exclude_bonds=0);", &desc, true, NULL, alloc));
    md_temp_end(temp);
}

UTEST_F(script, contacts_mirror_between_sets) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* res = mirror_residue_groups(&num, sys, 0, alloc);
    ASSERT_GT(num, (size_t)10);

    // Overlapping sets: the first two thirds of the residues against the last two thirds
    const size_t third = num / 3;
    md_contact_desc_t desc = { .group_a = res, .num_a = num - third, .group_b = res + third, .num_b = num - third, .cutoff = 4.0 };
    char src[512];
    snprintf(src, sizeof(src), "ga = residue(1:%zu); gb = residue(%zu:%zu); c = contacts(ga, gb, cutoff=4.0, exclude_bonds=0);", num - third, third + 1, num);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));

    // Protein against everything
    size_t num_prot = 0;
    md_bitfield_t* prot = mirror_residue_groups(&num_prot, sys, MD_FLAG_AMINO_ACID, alloc);
    ASSERT_GT(num_prot, (size_t)0);
    desc = (md_contact_desc_t){ .group_a = prot, .num_a = num_prot, .group_b = res, .num_b = num, .cutoff = 3.0 };
    EXPECT_TRUE(mirror_check(sys, "ga = residue(protein); gb = residue(:); c = contacts(ga, gb, cutoff=3.0, exclude_bonds=0);", &desc, true, NULL, alloc));
    md_temp_end(temp);
}

UTEST_F(script, contacts_mirror_bond_exclusion_and_separation) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* prot = mirror_residue_groups(&num, sys, MD_FLAG_AMINO_ACID, alloc);
    ASSERT_GT(num, (size_t)3);

    // Three bonds is the default of the script
    md_contact_desc_t desc = { .group_a = prot, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3 };
    md_contact_set_t with = {0};
    EXPECT_TRUE(mirror_check(sys, "ga = residue(protein); c = contacts(ga, cutoff=4.5);", &desc, true, &with, alloc));

    // Neighbouring residues are always in contact, the separation leaves the ones further along the chain
    desc.cutoff = 6.0;
    desc.min_separation = 2;
    md_contact_set_t sep = {0};
    EXPECT_TRUE(mirror_check(sys, "ga = residue(protein); c = contacts(ga, cutoff=6.0, min_separation=2);", &desc, true, &sep, alloc));
    for (size_t k = 0; k < sep.count; ++k) {
        EXPECT_GE(sep.j[k] - sep.i[k], 2u);
    }

    // Exclusion removes atom pairs across the peptide bonds, and never adds any
    desc.cutoff = 4.5;
    desc.min_separation = 0;
    desc.exclude_bonds = 0;
    md_contact_set_t without = {0};
    EXPECT_TRUE(mirror_check(sys, "ga = residue(protein); c = contacts(ga, cutoff=4.5, exclude_bonds=0);", &desc, true, &without, alloc));
    size_t sum_with = 0, sum_without = 0;
    for (size_t k = 0; k < with.count; ++k) sum_with += with.atom_pairs[k];
    for (size_t k = 0; k < without.count; ++k) sum_without += without.atom_pairs[k];
    EXPECT_LT(sum_with, sum_without);
    EXPECT_LE(with.count, without.count);
    md_temp_end(temp);
}

UTEST_F(script, contacts_mirror_triclinic) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->npt;
    size_t num_all = 0;
    md_bitfield_t* all = mirror_residue_groups(&num_all, sys, 0, alloc);
    ASSERT_GT(num_all, (size_t)10);
    // Every third residue: still spread over the whole cell, at a ninth of the cost of the reference
    md_bitfield_t* res = md_alloc(alloc, sizeof(md_bitfield_t) * num_all);
    size_t num = 0;
    for (size_t k = 0; k < num_all; k += 3) res[num++] = all[k];

    char src[512];
    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 5.0 };
    snprintf(src, sizeof(src), "ga = residue(1:3:%zu); c = contacts(ga, cutoff=5.0, exclude_bonds=0);", num_all);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));

    const size_t split = num / 3;
    desc = (md_contact_desc_t){ .group_a = res, .num_a = split, .group_b = res + split, .num_b = num - split, .cutoff = 6.0 };
    snprintf(src, sizeof(src), "ga = residue(1:3:%zu); gb = residue(%zu:3:%zu); c = contacts(ga, gb, cutoff=6.0, exclude_bonds=0);",
        1 + 3 * (split - 1), 1 + 3 * split, num_all);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));
    md_temp_end(temp);
}

// Groups sharing atoms, and atoms belonging to no group
UTEST_F(script, contacts_mirror_overlapping_groups) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* res = mirror_residue_groups(&num, sys, MD_FLAG_AMINO_ACID, alloc);
    ASSERT_GT(num, (size_t)6);
    // The amino acids are the leading components, so their indices are those of the components
    for (size_t k = 0; k < num; ++k) ASSERT_TRUE(md_system_component_flags(sys, k) & MD_FLAG_AMINO_ACID);

    // Windows of three consecutive residues, each overlapping the next by two. In the script: an array of selections
    const size_t num_win = num - 2;
    md_bitfield_t* win = md_alloc(alloc, sizeof(md_bitfield_t) * num_win);
    char windows[4096];
    size_t len = snprintf(windows, sizeof(windows), "{");
    for (size_t w = 0; w < num_win; ++w) {
        win[w] = md_bitfield_create(alloc);
        md_bitfield_or(&win[w], &res[w], &res[w + 1]);
        md_bitfield_or_inplace(&win[w], &res[w + 2]);
        if (w) len += snprintf(windows + len, sizeof(windows) - len, ", ");
        len += mirror_atom_range(windows + len, sizeof(windows) - len, sys, w, w + 3);
    }
    len += snprintf(windows + len, sizeof(windows) - len, "}");
    ASSERT_LT(len, sizeof(windows));

    char src[8192];
    md_contact_desc_t desc = { .group_a = win, .num_a = num_win, .cutoff = 4.0, .exclude_bonds = 3 };
    snprintf(src, sizeof(src), "ga = %s; c = contacts(ga, cutoff=4.0);", windows);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));

    desc = (md_contact_desc_t){ .group_a = win, .num_a = num_win, .group_b = res, .num_b = num, .cutoff = 4.0 };
    snprintf(src, sizeof(src), "ga = %s; gb = residue(protein); c = contacts(ga, gb, cutoff=4.0, exclude_bonds=0);", windows);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));
    md_temp_end(temp);
}

// With a parent per group, only groups of the same parent are neighbours
UTEST_F(script, contacts_mirror_min_separation_parent) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* res = mirror_residue_groups(&num, sys, MD_FLAG_AMINO_ACID, alloc);
    ASSERT_GT(num, (size_t)8);
    for (size_t k = 0; k < num; ++k) ASSERT_TRUE(md_system_component_flags(sys, k) & MD_FLAG_AMINO_ACID);

    // Two 'chains': residues [0, 7) and [7, num). In the script the parents are selections
    uint32_t* parent = md_alloc(alloc, sizeof(uint32_t) * num);
    for (size_t k = 0; k < num; ++k) parent[k] = k < 7 ? 0 : 1;
    char first[64], second[64], src[512];
    mirror_atom_range(first, sizeof(first), sys, 0, 7);
    mirror_atom_range(second, sizeof(second), sys, 7, num);

    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3, .min_separation = 3, .group_parent = parent };
    md_contact_set_t with_parent = {0}, without_parent = {0};
    snprintf(src, sizeof(src), "ga = residue(protein); c = contacts(ga, cutoff=4.5, min_separation=3, parent={%s, %s});", first, second);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, &with_parent, alloc));

    // Without parents the chain is one: this frame happens to have no contacts three residues apart
    desc.group_parent = NULL;
    EXPECT_TRUE(mirror_check(sys, "ga = residue(protein); c = contacts(ga, cutoff=4.5, min_separation=3);", &desc, false, &without_parent, alloc));

    // Residues 6 and 7 are consecutive, so in contact, but belong to different parents
    bool found_with = false, found_without = false;
    for (size_t k = 0; k < with_parent.count; ++k)    found_with    |= (with_parent.i[k] == 6 && with_parent.j[k] == 7);
    for (size_t k = 0; k < without_parent.count; ++k) found_without |= (without_parent.i[k] == 6 && without_parent.j[k] == 7);
    EXPECT_TRUE(found_with);
    EXPECT_FALSE(found_without);
    md_temp_end(temp);
}

UTEST_F(script, contacts_mirror_empty) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;

    // Groups too far apart
    md_bitfield_t g[2] = { md_bitfield_create(alloc), md_bitfield_create(alloc) };
    md_bitfield_set_bit(&g[0], 0);
    md_bitfield_set_bit(&g[1], 1);
    md_contact_desc_t desc = { .group_a = g, .num_a = 2, .cutoff = 0.01 };
    md_contact_set_t set = {0};
    EXPECT_TRUE(mirror_check(sys, "ga = {atom(1), atom(2)}; c = contacts(ga, cutoff=0.01, exclude_bonds=0);", &desc, false, &set, alloc));
    EXPECT_EQ(set.count, (size_t)0);

    // An empty B is an empty set: no contacts. Not the contacts within A, which is what no B means
    size_t num = 0;
    md_bitfield_t* res = mirror_residue_groups(&num, sys, 0, alloc);
    static const md_bitfield_t no_groups = {0};
    desc = (md_contact_desc_t){ .group_a = res, .num_a = num, .group_b = &no_groups, .num_b = 0, .cutoff = 5.0 };
    EXPECT_TRUE(mirror_check(sys, "ga = residue(:); gb = residue(atom(1) and not atom(1)); c = contacts(ga, gb, cutoff=5.0);", &desc, false, &set, alloc));
    EXPECT_EQ(set.count, (size_t)0);
    EXPECT_EQ(set.num_b, 0u);
    EXPECT_EQ(set.flags, (uint32_t)MD_CONTACT_FLAG_NONE);
    md_temp_end(temp);
}

// A plain selection is one group. Legal, but it is the likeliest mistake, and one which looks like a broken
// procedure: a single group has no contacts with itself, and two single groups are in contact or not.
UTEST_F(script, contacts_single_group_warnings) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* sys = &utest_fixture->ala;
    const struct { const char* src; const char* warning; } cases[] = {
        { "c = contacts(all, cutoff=5.0);", "never in contact with itself" },
        { "c = contacts(residue(1), residue(5), cutoff=6.0);", "0 or 1" },
        { "c = contacts(residue(:), cutoff=6.0);", NULL },
        { "c = contacts(residue(1:3), residue(5), cutoff=6.0);", NULL },
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_script_ir_t* ir = md_script_ir_create(alloc);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(cases[i].src), sys, NULL));
        const size_t num = md_script_ir_num_warnings(ir);
        const md_log_token_t* warnings = md_script_ir_warnings(ir);
        bool found = false;
        for (size_t w = 0; w < num; ++w) {
            if (cases[i].warning && str_find_str(NULL, warnings[w].text, str_from_cstr(cases[i].warning))) found = true;
        }
        if (cases[i].warning) {
            EXPECT_TRUE(found);
        } else {
            EXPECT_EQ((size_t)0, num);
        }
    }
    md_arena_allocator_destroy(alloc);
}

// The whole path of an application: compiled once, evaluated over the trajectory by several threads at once, the
// properties read afterwards. Every frame against the kernel query on the same frame, and the kernel against the
// reference.
typedef struct mirror_job_t {
    md_script_eval_t* eval;
    md_script_ir_t* ir;
    md_system_t* sys;
    uint32_t beg, end;
    bool ok;
} mirror_job_t;

static void mirror_job(void* data) {
    mirror_job_t* job = (mirror_job_t*)data;
    job->ok = md_script_eval_frame_range(job->eval, job->ir, job->sys, job->beg, job->end);
}

UTEST_F(script, contacts_mirror_query_over_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    ASSERT_TRUE(sys->trajectory);
    const size_t num_frames = md_trajectory_num_frames(sys->trajectory);
    ASSERT_GT(num_frames, (size_t)3);

    size_t num = 0;
    md_bitfield_t* prot = mirror_residue_groups(&num, sys, MD_FLAG_AMINO_ACID, alloc);
    md_contact_desc_t desc = { .group_a = prot, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3, .min_separation = 3 };

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "c = contacts(residue(protein), cutoff=4.5, min_separation=3);"
        "n = count(c);"
        "na = count(c, 'atom');"
        "d = degree(c);"), sys, NULL));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(eval);

    // Four ranges at once, as a pool of workers would
    test_hook_contact_local_prepares = 0;
    enum { NUM_JOBS = 4 };
    mirror_job_t jobs[NUM_JOBS];
    md_thread_t* threads[NUM_JOBS];
    for (uint32_t t = 0; t < NUM_JOBS; ++t) {
        jobs[t] = (mirror_job_t){ eval, ir, sys, (uint32_t)(num_frames * t / NUM_JOBS), (uint32_t)(num_frames * (t + 1) / NUM_JOBS), false };
        threads[t] = md_thread_create(mirror_job, &jobs[t]);
        ASSERT_TRUE(threads[t]);
    }
    for (uint32_t t = 0; t < NUM_JOBS; ++t) {
        md_thread_join(threads[t]);
        EXPECT_TRUE(jobs[t].ok);
    }
    EXPECT_EQ((size_t)0, test_hook_contact_local_prepares);

    const md_attributes_t* attr = md_script_eval_attributes(eval);
    const md_attribute_t* n  = md_attributes_find(attr, STR_LIT("script/n"));
    const md_attribute_t* na = md_attributes_find(attr, STR_LIT("script/na"));
    const md_attribute_t* d  = md_attributes_find(attr, STR_LIT("script/d"));
    ASSERT_TRUE(n && na && d);
    ASSERT_EQ(md_attribute_element_count(&d->format), num_frames * num);

    md_contact_query_t q;
    ASSERT_TRUE(md_contact_query_init(&q, &desc, sys, alloc));
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_system_state_init(&state, sys->atom.count));
    float* deg = md_alloc(alloc, sizeof(float) * num);

    size_t frames_differing = 0;
    for (size_t f = 0; f < num_frames; ++f) {
        ASSERT_TRUE(md_trajectory_load_frame(sys->trajectory, f, &state));
        md_contact_set_t set;
        ASSERT_TRUE(md_contact_query_eval(&set, &q, &state, alloc));
        md_array(ref_pair_t) ref = ref_contacts(&desc, sys, &state, alloc);
        EXPECT_TRUE(same_contacts(&set, ref, md_array_size(ref)));

        size_t atom_pairs = 0;
        memset(deg, 0, sizeof(float) * num);
        for (size_t k = 0; k < set.count; ++k) {
            atom_pairs += set.atom_pairs[k];
            deg[set.i[k]] += 1;
            deg[set.j[k]] += 1;
        }
        EXPECT_EQ((float)set.count, ((const float*)n->data)[f]);
        EXPECT_EQ((float)atom_pairs, ((const float*)na->data)[f]);
        for (size_t r = 0; r < num; ++r) {
            EXPECT_EQ(deg[r], ((const float*)d->data)[f * num + r]);
        }
        if (f && ((const float*)n->data)[f] != ((const float*)n->data)[f - 1]) frames_differing += 1;
        md_contact_set_free(&set);
    }
    // The contacts do change over the trajectory
    EXPECT_GT(frames_differing, (size_t)0);

    md_contact_query_free(&q);
    md_arena_allocator_destroy(alloc);
}

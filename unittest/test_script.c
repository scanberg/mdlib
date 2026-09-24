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
    static const str_t keywords[] = {
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

        bool seen_non_required = false;
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

            const bool required = (param->flags & (PARAM_OPTIONAL | PARAM_DEFAULT)) == 0;
            EXPECT_FALSE_MSG(required && seen_non_required, msg);              // required parameters come first
            EXPECT_FALSE_MSG((param->flags & PARAM_DEFAULT) && seen_optional, msg);  // no default after optional
            EXPECT_FALSE_MSG((param->flags & PARAM_OPTIONAL) && (param->flags & PARAM_DEFAULT), msg);
            if (param->flags & PARAM_DEFAULT) {
                EXPECT_TRUE_MSG(param->def_type.base_type != TYPE_UNDEFINED && is_scalar(param->def_type), msg);
            }
            seen_non_required |= !required;
            seen_optional     |= (param->flags & PARAM_OPTIONAL) != 0;
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

// ---------------------------------------------------------------------------------------------------------------
// Regressions
// ---------------------------------------------------------------------------------------------------------------

static const float* eval_property_data(md_script_eval_t* eval, const char* name) {
    char path[128];
    snprintf(path, sizeof(path), "script/%s", name);
    const md_attribute_t* attr = md_attributes_find(md_script_eval_attributes(eval), str_from_cstr(path));
    return attr ? (const float*)attr->data : NULL;
}

// distance_max used to call the routine for the minimum distance.
UTEST_F(script, distance_max_is_the_largest_distance) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(compiles(ir, "lo = distance_min(residue(1), residue(2)); hi = distance_max(residue(1), residue(2));", mol));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(eval != NULL);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, (str_t){0}, 0, num_frames));

    const float* lo = eval_property_data(eval, "lo");
    const float* hi = eval_property_data(eval, "hi");
    ASSERT_TRUE(lo && hi);
    for (uint32_t i = 0; i < num_frames; ++i) {
        EXPECT_TRUE(hi[i] > lo[i]);
    }
    md_script_eval_free(eval);
}

// A bit that is set in only one operand survives an xor, wherever it is, so the result cannot be limited to the
// range that the two operands share.
UTEST(script, xor_of_selections_covers_the_union_of_the_ranges) {
    EXPECT_TRUE(test_selection("atom(1:5) xor atom(3:8)",     "1100011100000000"));
    EXPECT_TRUE(test_selection("atom(1:3) xor atom(9:11)",    "1110000011100000"));
    EXPECT_TRUE(test_selection("atom(1:8) xor atom(1:8)",     "0000000000000000"));
    EXPECT_TRUE(test_selection("resname('SOL') xor all",      "0001111111111111"));
}

// contact_count started each element with the total of the previous ones.
UTEST_F(script, contact_count_is_counted_per_element) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(compiles(ir,
        "all4 = contact_count(residue(1:4), residue(5:8), 5.0);"
        "c1 = contact_count(residue(1), residue(5:8), 5.0);"
        "c2 = contact_count(residue(2), residue(5:8), 5.0);"
        "c3 = contact_count(residue(3), residue(5:8), 5.0);"
        "c4 = contact_count(residue(4), residue(5:8), 5.0);", mol));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(eval != NULL);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, (str_t){0}, 0, num_frames));

    const float* all4 = eval_property_data(eval, "all4");
    const float* c[4] = {
        eval_property_data(eval, "c1"), eval_property_data(eval, "c2"),
        eval_property_data(eval, "c3"), eval_property_data(eval, "c4"),
    };
    ASSERT_TRUE(all4 && c[0] && c[1] && c[2] && c[3]);
    for (uint32_t f = 0; f < num_frames; ++f) {
        for (int k = 0; k < 4; ++k) {
            EXPECT_EQ(c[k][f], all4[f * 4 + k]);
        }
    }
    md_script_eval_free(eval);
}

// An identifier that unpacks a right hand side which is known at compile time used to hold a null pointer.
UTEST_F(script, destructuring_a_constant) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    // Selections: the identifiers hold one residue each, and can be used in later statements
    ASSERT_TRUE(compiles(ir, "{a, b} = residue(1:2); na = count(a); nb = count(b) + 0.0 * distance(1, 2);", mol));
    {
        md_script_eval_t* eval = md_script_eval_create((uint32_t)md_trajectory_num_frames(mol->trajectory), ir, alloc);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, (str_t){0}, 0, 1));
        const float* nb = eval_property_data(eval, "nb");
        ASSERT_TRUE(nb != NULL);
        EXPECT_EQ(10.0f, nb[0]); // residue 2 of the peptide has ten atoms
        md_script_eval_free(eval);
    }

    // Vectors
    ASSERT_TRUE(compiles(ir, "{x, y} = vec2(1, 2); z = x + y + 0.0 * distance(1, 2);", mol));
    {
        md_script_eval_t* eval = md_script_eval_create((uint32_t)md_trajectory_num_frames(mol->trajectory), ir, alloc);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, (str_t){0}, 0, 1));
        const float* z = eval_property_data(eval, "z");
        ASSERT_TRUE(z != NULL);
        EXPECT_EQ(3.0f, z[0]);
        md_script_eval_free(eval);
    }
}

// Comparisons used to be rejected in every form: the array overload came first and swallowed a scalar operand, and
// a comparison could not be followed by an identifier.
UTEST_F(script, comparison_operators_compile) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)md_trajectory_num_frames(mol->trajectory);
    md_script_ir_t* ir = md_script_ir_create(alloc);

    const char* ok[] = {
        "d = distance(1, 2); b = d < 2.0;",
        "d = distance(1, 2); b = 2.0 > d;",
        "d = distance(1, 2); e = distance(1, 3); b = d <= e; c = d >= e; f = d == e;",
        "a = distance(1, 2) in residue(1:3); b = a < 1.01; c = 1.01 > a; d = a < a;",
        "a = distance(1, 2) in residue(1:3); b = (a < 1.01) and not (a > 1.005); c = (a >= 1.0) xor (a <= 1.5);",
        "d = distance(1, 2); b = (d < 2.0) and (d > 0.5);",
    };
    for (size_t i = 0; i < ARRAY_SIZE(ok); ++i) {
        EXPECT_TRUE(compiles(ir, ok[i], mol));
        if (md_script_ir_valid(ir)) {
            md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
            ASSERT_TRUE(eval != NULL);
            EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, (str_t){0}, 0, num_frames));
            md_script_eval_free(eval);
        }
    }

    // A comparison of arrays of different lengths is still an error
    EXPECT_FALSE(compiles(ir, "a = distance(1, 2) in residue(1:3); b = distance(1, 2) in residue(1:4); c = a < b;", mol));
}
}

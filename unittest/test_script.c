#include "utest.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <inttypes.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_bitset.h>
#include <core/md_bitop.inl>
#include <core/md_os.h>
#include <md_script.h>
#include <md_filter.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_trr.h>
#include "run_check.h"

// Test only instrumentation of the evaluator, see test_hook_proc_eval
#define MD_SCRIPT_TEST_HOOKS
#include <md_script.c>
#include <md_chem.h>
#include "component_util.h"
#include <md_pucker.h>

// Create molecule for evaulation
#define ATOM_COUNT 16
static float mol_x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
static float mol_y[] = {4,3,2,1,4,3,2,1,1,1,1,1,1,1,1,1};
static float mol_z[] = {3,2,1,4,3,2,1,2,2,2,2,2,2,2,2,2};
// The same coordinates packed, as a state holds them
static vec3_t mol_xyz[] = {
    {1,4,3}, {2,3,2}, {3,2,1}, {4,1,4}, {5,4,3}, {6,3,2}, {7,2,1}, {8,1,2},
    {1,1,2}, {2,1,2}, {3,1,2}, {4,1,2}, {5,1,2}, {6,1,2}, {7,1,2}, {8,1,2},
};
static md_particle_type_idx_t mol_ti[] = {0, 4, 0, 1, 2, 4, 2, 4, 0, 4, 0, 1, 2, 3, 2, 4};

#define PARTICLE_TYPE_COUNT 5
static md_atomic_number_t particle_type_z[] = {1, 2, 6, 7, 8};
static md_label_t particle_type_id[] = {
    MAKE_LABEL("H1"),
    MAKE_LABEL("HE"),
    MAKE_LABEL("CA"),
    MAKE_LABEL("N"),
    MAKE_LABEL("O")
};
static float particle_type_mass[] = {
    1.008f,
    4.002602f,
    12.011f,
    14.007f,
    15.999f
};
static float particle_type_radius[] = {
    1.2f,
    1.4f,
    1.7f,
    1.55f,
    1.52f
};

#define RES_COUNT 4
static md_label_t r_name[] = {MAKE_LABEL("SOL"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT"), MAKE_LABEL("PFT")};
static md_sequence_id_t  r_id[] = {1, 2, 3, 4};
static uint32_t    r_off[] = {0, 3, 8, 12, 16};

#define ENT_COUNT 3
static md_label_t    e_id[] = {MAKE_LABEL("WAT"), MAKE_LABEL("LYS"), MAKE_LABEL("PFT")};
// The kind is the low bits of the entity flags
static md_part_type_flags_t e_flags[] = {(md_part_type_flags_t)MD_PART_KIND_WATER, (md_part_type_flags_t)MD_PART_KIND_PEPTIDE, (md_part_type_flags_t)MD_PART_KIND_NON_POLYMER};
static str_t       e_desc[] = {MAKE_LABEL("water"), MAKE_LABEL("polypeptide"), MAKE_LABEL("non-polymer")};

#define INST_COUNT 4
static md_label_t        i_id[] = {MAKE_LABEL("A"), MAKE_LABEL("B"), MAKE_LABEL("C"), MAKE_LABEL("D")};
static uint32_t         i_off[] = {0, 3, 8, 12, 16};
static md_part_type_idx_t i_eidx[] = {0, 1, 2, 2};

md_system_t test_mol = {
    // A hand built system still needs a reference state: it is what evaluation reads coordinates
    // from when there is no trajectory frame.
    .reference = {
        .num_atoms = ATOM_COUNT,
        .xyz = mol_xyz,
    },
    .atom = {
        .count = ATOM_COUNT,
        .type_idx = mol_ti,
        .type = {
            .count = PARTICLE_TYPE_COUNT,
            .name = particle_type_id,
            .z = particle_type_z,
            .mass = particle_type_mass,
            .radius = particle_type_radius,
        },
    },
    .residue = {
        .count = RES_COUNT,
        .name = r_name,
        .seq_id = r_id,
        .atom_offset = r_off
    },
    .instance = {
        .count = INST_COUNT,
        .id = i_id,
        .res_offset = i_off,
    },
    .part = {
        .type = {
            .count = ENT_COUNT,
            .id = e_id,
            .flags = e_flags
        },
    },
};

// The run the fixture's trajectory is published as, and its number of frames
#define SCRIPT_RUN STR_LIT("run/ala")


static uint32_t script_frames(const md_system_t* sys) {
    char buf[64];
    const md_attribute_t* time = md_attributes_find_in(&sys->attributes, SCRIPT_RUN, STR_LIT("time"));
    return time ? time->format.shape[0] : 0;
}

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
    // Its models are the frames the evaluations below step through.
    ASSERT_TRUE(md_pdb_system_publish_run(&utest_fixture->ala, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), SCRIPT_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
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

// The union of what expr evaluates to, allocated with alloc
static bool eval_selection(md_bitset_t* out, str_t expr, md_system_t* mol, md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(mol);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    bool success = false;
    *out = md_bitset_empty_set();

    data_t data = {0};
    if (eval_expression(&data, expr, mol, temp_arena)) {
        if (data.type.base_type == TYPE_ATOMSET) {
            md_bitset_t* res = (md_bitset_t*)data.ptr;
            const int64_t len = type_info_array_len(data.type);
            md_bitmask_t b = md_bitmask_create(0, (uint32_t)mol->atom.count, temp_arena);
            for (int64_t i = 0; i < len; ++i) {
                md_bitmask_or(&b, res[i]);
            }
            *out = md_bitmask_extract(&b, alloc);
            success = true;
        }
    }
    md_temp_end(temp);
    return success;
}

// The first 64 members of a set as bits
static uint64_t bitset_bits64(md_bitset_t s) {
    return md_bitset_word(s, 0);
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
    md_bitset_t a = {0};
    md_bitset_t b = {0};

    if (!eval_selection(&a, str_from_cstr(exprA), &test_mol, temp_arena)) {
        printf("Failed evaluation of expression: '%s'\n", exprA);
        goto done;
    }
    if (!eval_selection(&b, str_from_cstr(exprB), &test_mol, temp_arena)) {
        printf("Failed evaluation of expression: '%s'\n", exprB);
        goto done;
    }

    uint64_t a_bits = bitset_bits64(a);
    uint64_t b_bits = bitset_bits64(b);
    bool cmp_res = bit_cmp(&a_bits, &b_bits, 0, ATOM_COUNT);
    if (!cmp_res) {
        printf("Got (A):\n");
        print_bits(&a_bits, ATOM_COUNT);
        printf("\nExpected (B)\n");
        print_bits(&b_bits, ATOM_COUNT);
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
            .compile_ir = ir,
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
    str_t src = STR_INIT(
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

// Two constant expressions which differ only in their operator are not the same expression: the second must
// not be given the value of the first.
UTEST_F(script, common_subexpression_elimination_operators) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp_scope);
    md_script_ir_t* ir = create_ir(arena);
    str_t src = STR_INIT(
        "a = 1.0 + 2.0;\n"
        "b = 1.0 - 2.0;\n"
        "c = 4.0 * 2.0;\n"
        "d = 4.0 / 2.0;\n"
    );

    ASSERT_TRUE(md_script_ir_compile_from_source(ir, src, &utest_fixture->ala, NULL));

    const char* names[4] = {"a", "b", "c", "d"};
    const float expected[4] = {3.0f, -1.0f, 8.0f, 2.0f};
    for (int i = 0; i < 4; ++i) {
        identifier_t* ident = get_identifier(ir, str_from_cstr(names[i]));
        ASSERT_TRUE(ident && ident->data && ident->data->ptr);
        EXPECT_NEAR(expected[i], *(const float*)ident->data->ptr, 1e-6f);
    }

    md_temp_end(temp_scope);
}

// The fingerprint of property 'name' once src is compiled into ir, 0 if it does not compile
static uint64_t compiled_fingerprint(md_script_ir_t* ir, const char* src, const char* name, md_system_t* sys) {
    md_script_ir_clear(ir);
    if (!md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL)) {
        return 0;
    }
    return md_script_ir_property_fingerprint(ir, str_from_cstr(name));
}

// A property's fingerprint is what it is computed from: not its name, its spelling, or the rest of the script
UTEST_F(script, property_fingerprint) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* sys = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    const uint64_t base = compiled_fingerprint(ir, "s = residue(1:2);\nd = distance(com(s), 10);", "d", sys);
    ASSERT_NE(0ull, base);

    // Spacing, comments and the names of the property and of what it uses
    EXPECT_EQ(base, compiled_fingerprint(ir, "# the first two residues\ns  =  residue(1:2);\n\n  d=distance(com(s),10); # to atom 10", "d", sys));
    EXPECT_EQ(base, compiled_fingerprint(ir, "s = residue(1:2);\nlen = distance(com(s), 10);", "len", sys));
    EXPECT_EQ(base, compiled_fingerprint(ir, "sel = residue(1:2);\nd = distance(com(sel), 10);", "d", sys));

    // Another statement before it changes nothing about it
    EXPECT_EQ(base, compiled_fingerprint(ir, "a = angle(1, 2, 3);\ns = residue(1:2);\nd = distance(com(s), 10);", "d", sys));

    // What it uses, however indirectly, and how it combines it
    EXPECT_NE(base, compiled_fingerprint(ir, "s = residue(1:3);\nd = distance(com(s), 10);", "d", sys));
    EXPECT_NE(base, compiled_fingerprint(ir, "s = residue(1:2);\nd = distance(com(s), 11);", "d", sys));
    EXPECT_NE(base, compiled_fingerprint(ir, "s = residue(1:2);\nd = distance(com(s), 10) * 2.0;", "d", sys));

    const uint64_t plus  = compiled_fingerprint(ir, "d = distance(1, 10) + 1.0;", "d", sys);
    const uint64_t minus = compiled_fingerprint(ir, "d = distance(1, 10) - 1.0;", "d", sys);
    EXPECT_NE(0ull, plus);
    EXPECT_NE(plus, minus);

    // The elements of a destructured value are different values
    md_script_ir_clear(ir);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("{x, y, z} = com(residue(1));"), sys, NULL));
    const uint64_t fx = md_script_ir_property_fingerprint(ir, STR_LIT("x"));
    const uint64_t fy = md_script_ir_property_fingerprint(ir, STR_LIT("y"));
    EXPECT_NE(0ull, fx);
    EXPECT_NE(0ull, fy);
    EXPECT_NE(fx, fy);
    EXPECT_EQ(0ull, md_script_ir_property_fingerprint(ir, STR_LIT("w")));

    md_arena_allocator_destroy(alloc);
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
        str_t src = STR_INIT(
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
            EXPECT_EQ(s1->data->type.base_type, TYPE_ATOMSET);
            EXPECT_EQ(s1->data->type.dim[0], 2);
        }
        if (s2) {
            EXPECT_EQ(s2->data->type.base_type, TYPE_ATOMSET);
            EXPECT_EQ(s2->data->type.dim[0], 3);
        }
        if (s) {
            EXPECT_EQ(s->data->type.base_type, TYPE_ATOMSET);
            EXPECT_EQ(s->data->type.dim[0], 5);
        }
    }

    {
        str_t src = STR_INIT(
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
        str_t src = STR_INIT(
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
        str_t src = STR_INIT(
            "s = residue(1:2:3);"
        );

        md_script_ir_clear(ir);
        bool result = md_script_ir_compile_from_source(ir, src, &test_mol, NULL);
        EXPECT_TRUE(result);
        identifier_t* s = get_identifier(ir, STR_LIT("s"));
        EXPECT_TRUE(s);
        if (s) {
            EXPECT_EQ(s->data->type.base_type, TYPE_ATOMSET);
            EXPECT_EQ(s->data->type.dim[0], 2);
        }
    }

    {
        // Stride on index ranges for coordinate selection: coord(1:2:5) -> atoms {1,3,5}
        str_t src = STR_INIT(
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
        str_t src = STR_INIT(
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
        str_t src = STR_INIT(
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
        str_t src = STR_INIT(
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
    md_bitset_t bf = {0};
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);

    if (!eval_selection(&bf, str_from_cstr(expr), &test_mol, arena)) {
        printf("Failed evaluation of expression: '%s'\n", expr);
        goto done;
    }
    uint64_t bits = bitset_bits64(bf);
    bool cmp_res = bit_cmp(&bits, &ref, 0, ATOM_COUNT);
    if (!cmp_res) {
        printf("Got:\n");
        print_bits(&bits, ATOM_COUNT);
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
    // Stride over residues equals explicit union of picks
    EXPECT_TRUE(compare_selection("residue(1:2:4)", "residue(1) or residue(3)"));


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
    // The script was written when the strands of a system without chain ids were chains found by position; they are
    // polymers now (chain() is the data's ids only). Read it as such until the file says polymer().
    {
        md_str_writer_t w = md_str_writer_create(alloc);
        str_t rest = script_src;
        size_t pos = 0;
        while (str_find_str(&pos, rest, STR_LIT("chain("))) {
            md_str_writer_push_str(&w, str_substr(rest, 0, pos));
            md_str_writer_push_str(&w, STR_LIT("polymer("));
            rest = str_substr(rest, pos + 6, SIZE_MAX);
        }
        md_str_writer_push_str(&w, rest);
        script_src = md_str_writer_to_str(&w);
    }

    md_script_ir_t* ir = md_script_ir_create(alloc);
    bool success = md_script_ir_compile_from_source(ir, script_src, &utest_fixture->amy, NULL);
    EXPECT_TRUE(success);
    if (!success) {
        md_script_ir_write_diagnostics(md_writer_stderr(), ir);
    }

    md_arena_allocator_destroy(alloc);
}

// sum() of a constant array is folded at compile time, ints converted to floats
UTEST_F(script, sum) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    data_t data = {0};
    ASSERT_TRUE(eval_expression(&data, STR_LIT("sum({1.0, 2.0, 3.5})"), &test_mol, alloc));
    EXPECT_EQ(6.5f, as_float(data));
    ASSERT_TRUE(eval_expression(&data, STR_LIT("sum({1, 2})"), &test_mol, alloc));
    EXPECT_EQ(3.0f, as_float(data));
    ASSERT_TRUE(eval_expression(&data, STR_LIT("sum({-1.0, 1.0}) + sum({2.0})"), &test_mol, alloc));
    EXPECT_EQ(2.0f, as_float(data));
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, implicit_conversion) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->amy;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sel = residue({1,2,3,4});"), mol, NULL));

    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("v = sdf(polymer(:), resname(\"PFT\"), 50);"), mol, NULL));

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

    md_bitset_t bf = {0};
    EXPECT_TRUE(eval_selection(&bf, STR_LIT("atom(1:20) and element('O') in polymer(:)"), mol, alloc));

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

        md_bitset_t bf = {0};
        const bool ok = eval_selection(&bf, str_from_cstr(expr), mol, alloc);
        EXPECT_TRUE(ok);
        if (!ok) {
            printf("  [%s] r=%.1f: evaluation failed outright\n", tag, WITHIN_RADII[i]);
            md_bitset_free(&bf, alloc);
            continue;
        }

        const size_t count = md_bitset_count(bf);
        if (count < prev_count) {
            printf("  [%s] r=%.1f selected %zu atoms, but r=%.1f already selected %zu\n",
                   tag, WITHIN_RADII[i], count, WITHIN_RADII[i - 1], prev_count);
        }
        EXPECT_GE(count, prev_count);
        prev_count = count;
        md_bitset_free(&bf, alloc);
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

        md_bitset_t bf_rng = {0};
        md_bitset_t bf_flt = {0};
        EXPECT_TRUE(eval_selection(&bf_rng, str_from_cstr(expr_rng), mol, alloc));
        EXPECT_TRUE(eval_selection(&bf_flt, str_from_cstr(expr_flt), mol, alloc));

        const size_t n_rng = md_bitset_count(bf_rng);
        const size_t n_flt = md_bitset_count(bf_flt);
        if (n_rng != n_flt) {
            printf("  r=%.1f: within(0:r) selected %zu, within(r) selected %zu\n", WITHIN_RADII[i], n_rng, n_flt);
        }
        EXPECT_EQ(n_flt, n_rng);

        md_bitset_free(&bf_rng, alloc);
        md_bitset_free(&bf_flt, alloc);
    }
    md_arena_allocator_destroy(alloc);
}

// Closest to what the report actually did: a large radius evaluated per frame over a trajectory,
// which is what a representation does. The spatial acc is rebuilt from each frame's coordinates and
// unit cell, so a frame dependent failure only shows up here.
UTEST_F(script, within_radius_over_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)script_frames(mol);
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
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
        md_script_eval_free(eval);
    }
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, dynamic_length) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    {
        str_t src = STR_INIT("sel1 = residue(within_z(1:50));");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
    }

    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, property_compute) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    uint32_t num_frames = (uint32_t)script_frames(mol);

    md_script_ir_t* ir = md_script_ir_create(alloc);

    {
        md_script_ir_clear(ir);
        str_t src = STR_INIT("num = count(residue(resname('ALA') and within(3.0, protein)));");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_INIT("{lin, plan, iso} = shape_weights(:);");
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
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        str_t src = STR_INIT("prop1 = distance_pair(com(resname(\"ALA\")), 1);");
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        ASSERT_TRUE(md_script_ir_valid(ir));

        identifier_t* prop1 = get_identifier(ir, STR_LIT("prop1"));
        ASSERT_TRUE(prop1);
        EXPECT_EQ(TYPE_FLOAT, prop1->data->type.base_type);
        EXPECT_EQ(1, prop1->data->type.dim[0]);

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

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
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, STR_LIT("sel = within_x(0:100);\np1  = distance(com(sel), 100);"), mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));

        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_NE(NULL, eval);
        EXPECT_EQ(1, md_script_eval_property_count(eval));
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

        md_script_eval_free(eval);
    }

    {
        str_t src = STR_INIT("s1 = count(within(10, residue(:)));");

        md_script_ir_clear(ir);
        md_script_ir_compile_from_source(ir, src, mol, NULL);
        EXPECT_TRUE(md_script_ir_valid(ir));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
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
    uint32_t num_frames = (uint32_t)script_frames(mol);

    md_script_ir_t* ir = md_script_ir_create(alloc);

    // A temporal property with a single value per frame.
    {
        md_script_ir_clear(ir);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("d = distance(1, 2);"), mol, NULL));
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_NE(NULL, eval);

        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        ASSERT_NE(NULL, attributes);
        // The evaluation's own frame axis, which every temporal property is checked against.
        const md_attribute_t* axis = md_attributes_find(attributes, STR_LIT("time"));
        ASSERT_TRUE(axis != NULL);
        EXPECT_EQ(num_frames, axis->format.shape[0]);

        const md_attribute_t* attr = md_attributes_find(attributes, STR_LIT("script/d"));
        ASSERT_TRUE(attr != NULL);
        EXPECT_EQ(axis, md_attributes_axis(attributes, attr));
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
        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
        EXPECT_TRUE(md_attributes_version(attributes, attr->id) > version);

        // One frame out of the middle, through the attribute rather than the pointer.
        float value = 0;
        md_attribute_slice_t slice = md_attribute_slice_1(num_frames / 2);
        EXPECT_EQ(1, md_attribute_extract_f32(&value, 1, attr, slice, md_unit_none()));
        EXPECT_EQ(((const float*)attr->data)[num_frames / 2], value);

        // The range covers every value that was evaluated.
        float minmax[2] = {0, 0};
        EXPECT_EQ(2, md_attribute_extract_f32(minmax, 2, range, md_attribute_slice_all(), md_unit_none()));
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

        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

        float minmax[2] = {0, 0};
        md_attribute_slice_t slice = md_attribute_slice_1(0);
        EXPECT_EQ(2, md_attribute_extract_f32(minmax, 2, ext, slice, md_unit_none()));
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

        EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

        // The bin coordinates are computed from the range the evaluation settled on.
        md_array(float) coords = md_array_create(float, num_bins, alloc);
        EXPECT_EQ(num_bins, md_attribute_extract_f32(coords, num_bins, bin, md_attribute_slice_all(), md_unit_none()));
        float x_range[2] = {0, 0};
        EXPECT_EQ(2, md_attribute_extract_f32(x_range, 2, range, md_attribute_slice_all(), md_unit_none()));
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
    const uint32_t num_frames = (uint32_t)script_frames(data->sys);
    if (md_script_eval_frame_range(data->eval, data->ir, data->sys, SCRIPT_RUN, 0, num_frames)) {
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

    const str_t script = STR_INIT("p1 = distance(1,10);");

    md_script_eval_t* eval[NUM_THREADS] = {0};
    md_thread_t* threads[NUM_THREADS] = {0};
    thread_data_t thread_data[NUM_THREADS] = {0};

    size_t num_frames = script_frames(mol);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, script, mol, NULL);
    EXPECT_TRUE(md_script_ir_valid(ir));
    ASSERT_EQ(1, md_script_ir_property_count(ir));

    md_script_eval_t* ref_eval = md_script_eval_create(num_frames, ir, alloc);
    EXPECT_EQ(num_frames, md_script_eval_frame_count(ref_eval));
    ASSERT_TRUE(md_script_eval_frame_range(ref_eval, ir, mol, SCRIPT_RUN, 0, (uint32_t)num_frames));

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
        str_t src = STR_INIT(
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
        
        EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, x->node), -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        EXPECT_EQ(mol->atom.count, md_bitset_count(vis.atom_mask));
        
        md_script_vis_clear(&vis);
        EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, sx->node), -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        size_t res_atom_count = md_system_residue_atom_count(mol, 4);
        size_t pop_count = md_bitset_count(vis.atom_mask);
        EXPECT_EQ(res_atom_count, pop_count);
    }

    md_arena_allocator_destroy(alloc);
}

// ### NAMED ARGUMENTS ###

UTEST(script, named_args_signature_table) {
    const str_t keywords[] = {
        STR_INIT("in"), STR_INIT("of"), STR_INIT("out"), STR_INIT("and"), STR_INIT("or"), STR_INIT("xor"), STR_INIT("not"),
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
    if (x->type.base_type == TYPE_ATOMSET) {
        const int64_t len = type_info_array_len(x->type);
        const md_bitset_t* bx = (const md_bitset_t*)x->ptr;
        const md_bitset_t* by = (const md_bitset_t*)y->ptr;
        for (int64_t i = 0; i < len; ++i) {
            const size_t pop = md_bitset_count(bx[i]);
            if (pop != md_bitset_count(by[i])) return false;
            md_bitset_t both = md_bitset_and(bx[i], by[i], alloc);
            if (md_bitset_count(both) != pop) return false;
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
    str_t positional = STR_INIT(
        "s1 = resname(\"ALA\")[2:8];\n"
        "d1 = distance(10,30);\n"
        "a1 = angle(2,1,3) in resname(\"ALA\");\n"
        "r = rdf(element('C'), element('H'), 10.0);\n"
        "v = sdf(s1, element('H'), 10.0);\n"
        "{lin,plan,iso} = shape_weights(all);\n");
    str_t named = STR_INIT(
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

// Index selectors (residue, chain, instance; ring in test_ring.c) count what has an atom in their context, in order,
// each cut to the context. A position past the last selects nothing within a context, and is an error without one.
UTEST(script, index_selectors_in_context) {
    // SOL 1-3, LYS 4-8, PFT 9-12, PFT 13-16
    EXPECT_TRUE(test_selection("residue(2) in union(resname('PFT'))",   "0000000000001111"));
    EXPECT_TRUE(test_selection("residue(1:2) in union(resname('PFT'))", "0000000011111111"));
    // atom(3:9) has atoms of SOL, LYS and the first PFT
    EXPECT_TRUE(test_selection("residue(2) in atom(3:9)",  "0001111100000000"));
    EXPECT_TRUE(test_selection("residue(3) in atom(3:9)",  "0000000010000000"));    // Cut to the context
    EXPECT_TRUE(test_selection("residue(:) in atom(3:9)",  "0011111110000000"));
    EXPECT_TRUE(test_selection("residue(4) in atom(3:9)",  "0000000000000000"));    // Past the last: nothing
    EXPECT_TRUE(test_selection("residue(1) in atom(4:16)", "0001111100000000"));
    // One context per residue of resname('PFT')
    EXPECT_TRUE(test_selection("residue(1) in resname('PFT')", "0000000011111111"));
    EXPECT_TRUE(test_selection("residue(2) in resname('PFT')", "0000000000000000"));

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_script_ir_t* ir = md_script_ir_create(alloc);
    EXPECT_TRUE(named_args_compile_fails_with(ir, &test_mol, "x = residue(5);",  "past the last residue (4)"));
    EXPECT_TRUE(named_args_compile_fails_with(ir, &test_mol, "x = residue(0);",  "before the first residue"));
    // test_mol gives no chain ids, so it has no chains
    EXPECT_TRUE(named_args_compile_fails_with(ir, &test_mol, "x = chain(1);",    "no chain ids"));
    // Instances are not a word of the script
    EXPECT_TRUE(named_args_compile_fails_with(ir, &test_mol, "x = instance(1);", "instance"));
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
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
        .compile_ir = ir,     // Binding arguments is part of a compilation
        .sys = &utest_fixture->ala,
        .temp_alloc = alloc,
        .alloc = alloc,
    };

    const proc_sig_t sig_default = {
        .proc = STR_INIT("angle"), .num_params = 3,
        .param = {
            {STR_INIT("a"), PARAM_REQUIRED},
            {STR_INIT("b"), PARAM_REQUIRED},
            {STR_INIT("c"), PARAM_DEFAULT, .def_type = TI_INT, .def = {._int = 7}},
        },
    };
    const proc_sig_t sig_optional = {
        .proc = STR_INIT("angle"), .num_params = 3,
        .param = {
            {STR_INIT("a"), PARAM_REQUIRED},
            {STR_INIT("b"), PARAM_OPTIONAL},
            {STR_INIT("c"), PARAM_OPTIONAL},
        },
    };
    const proc_sig_t sig_kw_only = {
        .proc = STR_INIT("angle"), .num_params = 3,
        .param = {
            {STR_INIT("a"), PARAM_REQUIRED},
            {STR_INIT("b"), PARAM_REQUIRED},
            {STR_INIT("c"), PARAM_KW_ONLY},
        },
    };
    const proc_sig_t sig_shadow = {
        .proc = STR_INIT("angle"), .num_params = 3,
        .param = {
            {STR_INIT("count"),   PARAM_REQUIRED},
            {STR_INIT("residue"), PARAM_REQUIRED},
            {STR_INIT("c"),       PARAM_REQUIRED},
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

// ### attr() ###

static md_attribute_id_t publish_f64(md_attributes_t* t, const char* path, md_attribute_format_t fmt, md_attribute_flags_t flags, md_unit_t unit, const double* data) {
    return md_attributes_create(t, &(md_attribute_desc_t){
        .path = str_from_cstr(path), .format = fmt, .flags = flags, .unit = unit,
        .data = data, .byte_size = md_attribute_element_count(&fmt) * sizeof(double)});
}

static md_attribute_format_t series_fmt(uint32_t n, uint32_t components) {
    return (md_attribute_format_t){.type = MD_ATTRIBUTE_TYPE_F64, .components = components, .rank = 1, .shape = {n}};
}

static bool compiles(md_script_ir_t* ir, const char* src, const md_system_t* sys) {
    md_script_ir_clear(ir);
    md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL);
    return md_script_ir_valid(ir);
}

// attr() reads a temporal attribute of the system at the frame being evaluated. The attribute keeps
// its own axis - here one sampled twice as often as the trajectory - and the row is found by time.
UTEST(script, attr_reads_a_temporal_attribute_at_the_frame) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {.alloc = arena};
    md_system_state_t state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro")));
    ASSERT_TRUE(md_trr_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), STR_LIT("run/t"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));

    md_attributes_t* t = &sys.attributes;
    const md_attribute_t* run_time = md_attributes_find(t, STR_LIT("run/t/time"));
    ASSERT_TRUE(run_time != NULL);
    const uint32_t F = run_time->format.shape[0];
    ASSERT_GT(F, 2u);
    const double* frame_times = (const double*)run_time->data;
    const md_unit_t time_unit = run_time->unit;
    ASSERT_TRUE(md_unit_equal(time_unit, md_unit_picosecond()));

    // Twice the rate: every frame time and the midpoints between them.
    const uint32_t R = 2 * F - 1;
    double* obs_time = md_alloc(arena, R * sizeof(double));
    double* obs_val  = md_alloc(arena, R * sizeof(double));
    double* obs_vec  = md_alloc(arena, R * 3 * sizeof(double));
    double* obs_ten  = md_alloc(arena, R * 9 * sizeof(double));
    for (uint32_t r = 0; r < R; ++r) {
        obs_time[r] = (r % 2 == 0) ? frame_times[r / 2] : 0.5 * (frame_times[r / 2] + frame_times[r / 2 + 1]);
        obs_val[r]  = 2.0 * obs_time[r];
        for (int k = 0; k < 3; ++k) obs_vec[r * 3 + k] = obs_time[r] + k;
        for (int k = 0; k < 9; ++k) obs_ten[r * 9 + k] = k;
    }
    ASSERT_NE(publish_f64(t, "run/t/obs/time",  series_fmt(R, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_picosecond(), obs_time), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/t/obs/value", series_fmt(R, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_kelvin(), obs_val), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/t/obs/vec",   series_fmt(R, 3), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), obs_vec), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/t/obs/tensor",
        (md_attribute_format_t){.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3, .shape = {R, 3, 3}},
        MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), obs_ten), MD_ATTRIBUTE_INVALID);
    const double constant[3] = {1, 2, 3};
    ASSERT_NE(publish_f64(t, "run/t/static", series_fmt(3, 1), MD_ATTRIBUTE_FLAG_NONE, md_unit_none(), constant), MD_ATTRIBUTE_INVALID);

    md_script_ir_t* ir = md_script_ir_create(arena);

    // Relative to the run, evaluated over every frame: the value at frame f is the one at the
    // frame's own time, not at row f of the finer axis.
    ASSERT_TRUE(compiles(ir, "e = attr(\"obs/value\");", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        const md_attribute_t* e = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/e"));
        ASSERT_TRUE(e != NULL);
        EXPECT_TRUE(md_unit_equal(e->unit, md_unit_kelvin()));
        for (uint32_t f = 0; f < F; ++f) {
            EXPECT_NEAR((float)(2.0 * frame_times[f]), ((const float*)e->data)[f], 1.0e-3f);
        }
        md_script_eval_free(eval);
    }

    // A vector is float[3], a tensor float[3][3] - which a property cannot hold as it is, so it is
    // flattened - and a path may also be given in full.
    ASSERT_TRUE(compiles(ir, "v = attr(\"obs/vec\");", &sys));
    ASSERT_TRUE(compiles(ir, "m = flatten(attr(\"run/t/obs/tensor\"));", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        const md_attribute_t* m = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/m"));
        ASSERT_TRUE(m != NULL);
        EXPECT_EQ(9u, m->format.shape[1]);
        EXPECT_EQ(5.0f, ((const float*)m->data)[9 * 3 + 5]);
        md_script_eval_free(eval);
    }

    // What attr() refuses, each for its own reason.
    EXPECT_FALSE(compiles(ir, "x = attr(\"obs/missing\");", &sys));
    EXPECT_FALSE(compiles(ir, "x = attr(\"static\");", &sys));
    EXPECT_FALSE(compiles(ir, "x = attr(\"run/t/time\" + 1);", &sys));

    // A second run holding the same relative path makes the short form ambiguous, and it is an
    // error rather than a choice. Naming the run in full still works.
    ASSERT_NE(publish_f64(t, "run/u/time", series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, time_unit, frame_times), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/u/obs/time",  series_fmt(R, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_picosecond(), obs_time), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/u/obs/value", series_fmt(R, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), obs_val), MD_ATTRIBUTE_INVALID);
    EXPECT_FALSE(compiles(ir, "e = attr(\"obs/value\");", &sys));
    EXPECT_TRUE(compiles(ir, "e = attr(\"run/u/obs/value\");", &sys));

    // import() is gone: files are loaded into the attributes, and read with attr().
    EXPECT_FALSE(compiles(ir, "t = import(\"" MD_UNITTEST_DATA_DIR "/ener.edr\");", &sys));
    EXPECT_FALSE(compiles(ir, "t = import(\"" MD_UNITTEST_DATA_DIR "/energy.xvg\");", &sys));

    md_script_ir_free(ir);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// A temporal f64 series along the run, under a label of its own when one is given
static md_attribute_id_t publish_series(md_attributes_t* t, const char* path, const char* label, uint32_t n, uint32_t components, md_unit_t unit, const double* data) {
    return md_attributes_create(t, &(md_attribute_desc_t){
        .path = str_from_cstr(path), .label = label ? str_from_cstr(label) : (str_t){0},
        .format = series_fmt(n, components), .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = unit,
        .data = data, .byte_size = (size_t)n * components * sizeof(double)});
}

// Whether src fails to compile with an error that says what because says
static bool fails_because(md_script_ir_t* ir, const char* src, const md_system_t* sys, const char* because) {
    if (compiles(ir, src, sys)) return false;
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        if (str_find_str(NULL, md_script_ir_errors(ir)[i].text, str_from_cstr(because))) return true;
    }
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        printf("  %s: " STR_FMT "\n", src, STR_ARG(md_script_ir_errors(ir)[i].text));
    }
    return false;
}

static str_t eval_label(const md_script_eval_t* eval, const char* path, size_t i) {
    const md_attributes_t* table = md_script_eval_attributes(eval);
    const md_attribute_t* attr = md_attributes_find(table, str_from_cstr(path));
    if (!attr || attr->format.type != MD_ATTRIBUTE_TYPE_STR || i >= attr->format.shape[0]) return (str_t){0};
    return md_attribute_str(table, attr, i);
}

// attr("edr/coul_sr_*"): every attribute the pattern matches, an array over them in path order, whose members are named
// in the evaluation's table (script/<ident>/labels) wherever what a property holds is still that array.
UTEST(script, attr_patterns) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {.alloc = arena};
    md_system_state_t state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro")));
    ASSERT_TRUE(md_trr_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), STR_LIT("run/t"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));

    md_attributes_t* t = &sys.attributes;
    const md_attribute_t* run_time = md_attributes_find(t, STR_LIT("run/t/time"));
    ASSERT_TRUE(run_time != NULL);
    const uint32_t F = run_time->format.shape[0];
    const double* frame_times = (const double*)run_time->data;
    const md_unit_t energy = md_unit_kelvin();  // any unit will do, as long as the members agree on it

    // Member k is k + 1 times the frame number. Published out of order: the array is in path order.
    double* val[4];
    for (int k = 0; k < 4; ++k) {
        val[k] = md_alloc(arena, F * sizeof(double));
        for (uint32_t f = 0; f < F; ++f) val[k][f] = (double)(k + 1) * f;
    }
    double* vec = md_alloc(arena, F * 3 * sizeof(double));
    for (uint32_t i = 0; i < F * 3; ++i) vec[i] = i;
    ASSERT_NE(publish_series(t, "run/t/edr/coul_sr_c", "Coul-SR:C", F, 1, energy, val[2]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/coul_sr_a", "Coul-SR:A", F, 1, energy, val[0]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/coul_sr_b", "Coul-SR:B", F, 1, energy, val[1]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/lj_sr_a",   NULL, F, 1, energy, val[0]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/lj_sr_b",   NULL, F, 1, energy, val[1]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/unit_a",    NULL, F, 1, energy, val[0]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/unit_b",    NULL, F, 1, md_unit_none(), val[1]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/shape_a",   NULL, F, 1, energy, val[0]), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/shape_b",   NULL, F, 3, energy, vec), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/box_a",     NULL, F, 3, energy, vec), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/box_b",     NULL, F, 3, energy, vec), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/t/edr/deeper/coul_sr_z", NULL, F, 1, energy, val[3]), MD_ATTRIBUTE_INVALID);
    const double constant[3] = {1, 2, 3};
    ASSERT_NE(publish_f64(t, "run/t/edr/const_a", series_fmt(3, 1), MD_ATTRIBUTE_FLAG_NONE, energy, constant), MD_ATTRIBUTE_INVALID);

    md_script_ir_t* ir = md_script_ir_create(arena);

    // In path order, a '*' within one segment (not edr/deeper/coul_sr_z), named by their labels
    ASSERT_TRUE(compiles(ir, "e = attr(\"edr/coul_sr_*\");", &sys));
    {
        const ast_node_t* node = labelled_pattern(ir->property_nodes[0], 0);
        ASSERT_TRUE(node != NULL);
        ASSERT_EQ(3u, (uint32_t)md_array_size(node->attr_members));
        EXPECT_EQ(3, node->data.type.dim[0]);
        EXPECT_EQ(0, node->data.type.dim[1]);
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        const md_attribute_t* e = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/e"));
        ASSERT_TRUE(e != NULL);
        EXPECT_EQ(3u, e->format.shape[1]);
        EXPECT_TRUE(md_unit_equal(e->unit, energy));
        for (uint32_t f = 0; f < F; ++f) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_EQ((float)((k + 1) * f), ((const float*)e->data)[f * 3 + k]);
            }
        }
        EXPECT_TRUE(str_eq(eval_label(eval, "script/e/labels", 0), STR_LIT("Coul-SR:A")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/e/labels", 1), STR_LIT("Coul-SR:B")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/e/labels", 2), STR_LIT("Coul-SR:C")));
        md_script_eval_free(eval);
    }

    // Without labels of their own, by the segment the '*' stood in. Absolute, and '*' in more than one place.
    ASSERT_TRUE(compiles(ir, "l = attr(\"run/t/edr/lj_*_*\");", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/l/labels", 0), STR_LIT("lj_sr_a")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/l/labels", 1), STR_LIT("lj_sr_b")));
        md_script_eval_free(eval);
    }
    ASSERT_TRUE(compiles(ir, "l = attr(\"*/lj_sr_*\");", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/l/labels", 0), STR_LIT("edr/lj_sr_a")));
        md_script_eval_free(eval);
    }

    // Through names and arithmetic that keeps the members' axis, and not past what takes it apart
    ASSERT_TRUE(compiles(ir,
        "x = attr(\"edr/coul_sr_*\");\n"
        "y = x * 2 - x;\n"
        "z = -y;\n"
        "w = x + attr(\"edr/lj_sr_*\")[1];\n"
        "s = x[1] + x[2];\n"
        "u = x + attr(\"edr/coul_sr_*\");\n"
        "v = x[1:2];\n"
        "m = x[1:2] + attr(\"edr/lj_sr_*\");\n"
        "k = attr(\"edr/lj_sr_*\") + attr(\"edr/coul_sr_*\")[2:3];\n", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/y/labels", 2), STR_LIT("Coul-SR:C")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/z/labels", 0), STR_LIT("Coul-SR:A")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/w/labels", 1), STR_LIT("Coul-SR:B")));
        EXPECT_TRUE(str_eq(eval_label(eval, "script/u/labels", 1), STR_LIT("Coul-SR:B")));
        EXPECT_TRUE(md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/s/labels")) == NULL);
        EXPECT_TRUE(md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/v/labels")) == NULL);
        // Each of these is as much one array's members as the other's
        EXPECT_TRUE(md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/m/labels")) == NULL);
        EXPECT_TRUE(md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/k/labels")) == NULL);
        const md_attribute_t* z = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/z"));
        ASSERT_TRUE(z != NULL);
        EXPECT_EQ(-(float)(3 * 2), ((const float*)z->data)[2 * 3 + 2]);
        md_script_eval_free(eval);
    }

    // Summed per frame, in the unit of the members; min and max keep it too
    ASSERT_TRUE(compiles(ir, "s = sum(attr(\"edr/coul_sr_*\"));", &sys));
    EXPECT_TRUE(md_unit_equal(ir->property_nodes[0]->data.unit[1], energy));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        const md_attribute_t* sum = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/s"));
        ASSERT_TRUE(sum != NULL);
        EXPECT_TRUE(md_unit_equal(sum->unit, energy));
        for (uint32_t f = 0; f < F; ++f) {
            EXPECT_EQ((float)(6 * f), ((const float*)sum->data)[f]);   // 1 + 2 + 3 times the frame
        }
        EXPECT_TRUE(md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/s/labels")) == NULL);
        md_script_eval_free(eval);
    }
    ASSERT_TRUE(compiles(ir, "m = max(attr(\"edr/coul_sr_*\"));", &sys));
    EXPECT_TRUE(md_unit_equal(ir->property_nodes[0]->data.unit[1], energy));

    // One match is an array of one, vectors an array of vectors
    ASSERT_TRUE(compiles(ir, "o = attr(\"edr/coul_sr_c*\");", &sys));
    EXPECT_EQ(1, ir->property_nodes[0]->data.type.dim[0]);
    ASSERT_TRUE(compiles(ir, "b = attr(\"edr/box_*\");", &sys));
    ASSERT_EQ(1u, (uint32_t)md_array_size(ir->type_checked_expressions));
    EXPECT_EQ(2, ir->type_checked_expressions[0]->data.type.dim[0]);
    EXPECT_EQ(3, ir->type_checked_expressions[0]->data.type.dim[1]);
    ASSERT_TRUE(compiles(ir, "b = flatten(attr(\"edr/box_*\"));", &sys));
    {
        md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, &sys, STR_LIT("run/t"), 0, F));
        const md_attribute_t* b = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/b"));
        ASSERT_TRUE(b != NULL);
        ASSERT_EQ(6u, b->format.shape[1]);
        // Frame 1 of each: the second vector of the series, box_a then box_b
        for (int k = 0; k < 6; ++k) EXPECT_EQ((float)(3 + k % 3), ((const float*)b->data)[6 + k]);
        md_script_eval_free(eval);
    }

    // Refused: nothing matches, only what attr() does not read, a frame axis, members that disagree on what a
    // value is or on its unit
    EXPECT_TRUE(fails_because(ir, "x = attr(\"edr/none_*\");", &sys, "nothing in any run matches"));
    EXPECT_TRUE(fails_because(ir, "x = attr(\"edr/const_*\");", &sys, "none that attr() reads"));
    EXPECT_TRUE(fails_because(ir, "x = attr(\"tim*\");", &sys, "only frame axes"));
    EXPECT_TRUE(fails_because(ir, "x = attr(\"edr/unit_*\");", &sys, "does not make one array"));
    EXPECT_TRUE(fails_because(ir, "x = attr(\"edr/shape_*\");", &sys, "does not make one array"));

    // What it matched is part of what a property is computed from: one more match is another property
    ASSERT_TRUE(compiles(ir, "e = attr(\"edr/coul_sr_*\");", &sys));
    const uint64_t before = md_script_ir_property_fingerprint(ir, STR_LIT("e"));
    ASSERT_NE(publish_series(t, "run/t/edr/coul_sr_d", "Coul-SR:D", F, 1, energy, val[3]), MD_ATTRIBUTE_INVALID);
    ASSERT_TRUE(compiles(ir, "e = attr(\"edr/coul_sr_*\");", &sys));
    EXPECT_NE(before, md_script_ir_property_fingerprint(ir, STR_LIT("e")));
    EXPECT_EQ(4, ir->property_nodes[0]->data.type.dim[0]);

    // Matches in two runs: the run has to be named, as for a path
    ASSERT_NE(publish_f64(t, "run/u/time", series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, run_time->unit, frame_times), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_series(t, "run/u/edr/coul_sr_a", "Coul-SR:A", F, 1, energy, val[0]), MD_ATTRIBUTE_INVALID);
    EXPECT_TRUE(fails_because(ir, "e = attr(\"edr/coul_sr_*\");", &sys, "'run/u/edr/coul_sr_*'"));
    EXPECT_TRUE(fails_because(ir, "e = attr(\"run/*/edr/coul_sr_a\");", &sys, "'run/t/edr/coul_sr_a'"));
    EXPECT_TRUE(compiles(ir, "e = attr(\"run/t/edr/coul_sr_*\");", &sys));
    EXPECT_TRUE(compiles(ir, "e = attr(\"run/u/edr/coul_sr_*\");", &sys));

    md_script_ir_free(ir);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// '*' stands for any run of characters, anywhere and as many times as it is written
UTEST(script, glob_match) {
    EXPECT_TRUE (glob_match(STR_LIT("CA"),    STR_LIT("CA")));
    EXPECT_FALSE(glob_match(STR_LIT("CA"),    STR_LIT("CAX")));
    EXPECT_TRUE (glob_match(STR_LIT("C*"),    STR_LIT("C")));
    EXPECT_TRUE (glob_match(STR_LIT("C*"),    STR_LIT("CA")));
    EXPECT_FALSE(glob_match(STR_LIT("C*"),    STR_LIT("AC")));
    EXPECT_TRUE (glob_match(STR_LIT("*A"),    STR_LIT("CA")));
    EXPECT_FALSE(glob_match(STR_LIT("*A"),    STR_LIT("AC")));
    // Used to match anything at all
    EXPECT_TRUE (glob_match(STR_LIT("*H*"),   STR_LIT("OHO")));
    EXPECT_TRUE (glob_match(STR_LIT("*H*"),   STR_LIT("H")));
    EXPECT_FALSE(glob_match(STR_LIT("*H*"),   STR_LIT("OC")));
    // Used to let the parts before and after the '*' overlap
    EXPECT_FALSE(glob_match(STR_LIT("AB*BC"), STR_LIT("ABC")));
    EXPECT_TRUE (glob_match(STR_LIT("AB*BC"), STR_LIT("ABBC")));
    EXPECT_TRUE (glob_match(STR_LIT("AB*BC"), STR_LIT("ABxyBC")));
    // Used to be refused: '*' in the middle, and more than two
    EXPECT_TRUE (glob_match(STR_LIT("C*1*"),  STR_LIT("C21x")));
    EXPECT_FALSE(glob_match(STR_LIT("C*1*"),  STR_LIT("C2")));
    EXPECT_TRUE (glob_match(STR_LIT("A*B*C"), STR_LIT("AxxBxxC")));
    EXPECT_FALSE(glob_match(STR_LIT("A*B*C"), STR_LIT("AxxCxxB")));
    EXPECT_TRUE (glob_match(STR_LIT("**A"),   STR_LIT("A")));
    EXPECT_TRUE (glob_match(STR_LIT("*"),     STR_LIT("anything")));
    // What a '*' takes has to grow past a false start
    EXPECT_TRUE (glob_match(STR_LIT("*AB"),   STR_LIT("AAB")));
    EXPECT_TRUE (glob_match(STR_LIT("*ABC*ABD"), STR_LIT("ABCABCABD")));
    EXPECT_FALSE(glob_match(STR_LIT("*ABC*ABD"), STR_LIT("ABCABCABE")));
    // Case is part of the name
    EXPECT_FALSE(glob_match(STR_LIT("c*"),    STR_LIT("CA")));
    // The selectors: nothing matches an empty name, and an empty pattern is no pattern
    EXPECT_FALSE(match_query(STR_LIT("*"), STR_LIT("")));
    EXPECT_FALSE(match_query(STR_LIT(""),  STR_LIT("CA")));
    EXPECT_FALSE(validate_query(STR_LIT("")));
    EXPECT_TRUE (validate_query(STR_LIT("*A*B*")));

    // In paths, a '*' does not reach across a '/'
    EXPECT_TRUE (glob_match_path(STR_LIT("edr/*"),        STR_LIT("edr/potential")));
    EXPECT_FALSE(glob_match_path(STR_LIT("edr/*"),        STR_LIT("edr/sub/potential")));
    EXPECT_FALSE(glob_match_path(STR_LIT("edr*"),         STR_LIT("edr/potential")));
    EXPECT_TRUE (glob_match_path(STR_LIT("*/potential"),  STR_LIT("edr/potential")));
    EXPECT_TRUE (glob_match_path(STR_LIT("xvg/*/coul_*"), STR_LIT("xvg/energy/coul_sr")));
    EXPECT_FALSE(glob_match_path(STR_LIT("xvg/*/coul_*"), STR_LIT("xvg/energy/lj_sr")));

    md_allocator_i* arena = md_vm_arena_create(MEGABYTES(1));
    EXPECT_TRUE(str_eq(glob_varying_segments(STR_LIT("edr/coul_*"),   STR_LIT("edr/coul_sol"),   arena), STR_LIT("coul_sol")));
    EXPECT_TRUE(str_eq(glob_varying_segments(STR_LIT("xvg/*/coul_*"), STR_LIT("xvg/e/coul_x"),   arena), STR_LIT("e/coul_x")));
    md_vm_arena_destroy(arena);
}

// name(), resname() and instance() select by the same patterns
UTEST(script, name_patterns) {
    md_allocator_i* arena = md_vm_arena_create(MEGABYTES(64));
    const struct { const char* pattern; int kind; char a; char b; } cases[] = {
        {"*A*", 0, 'A', 0},     // contains A
        {"*E",  1, 'E', 0},     // ends with E
        {"H*",  2, 'H', 0},     // begins with H
        {"C*A", 3, 'C', 'A'},   // begins with C, ends with A
        {"*",   4, 0, 0},       // anything
    };
    for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
        char src[64];
        snprintf(src, sizeof(src), "name(\"%s\")", cases[c].pattern);
        md_bitset_t sel = {0};
        ASSERT_TRUE(eval_selection(&sel, str_from_cstr(src), &test_mol, arena));
        for (size_t i = 0; i < test_mol.atom.count; ++i) {
            const str_t n = md_atom_name(&test_mol.atom, i);
            bool expected = false;
            switch (cases[c].kind) {
            case 0: expected = str_find_char(NULL, n, cases[c].a); break;
            case 1: expected = n.len && n.ptr[n.len - 1] == cases[c].a; break;
            case 2: expected = n.len && n.ptr[0] == cases[c].a; break;
            case 3: expected = n.len >= 2 && n.ptr[0] == cases[c].a && n.ptr[n.len - 1] == cases[c].b; break;
            case 4: expected = n.len > 0; break;
            }
            char msg[96];
            snprintf(msg, sizeof(msg), "%s, atom %zu '" STR_FMT "'", cases[c].pattern, i, STR_ARG(n));
            EXPECT_TRUE_MSG(expected == md_bitset_test(sel, (uint32_t)i), msg);
        }
    }
    md_vm_arena_destroy(arena);
}

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
    const uint32_t num_frames = (uint32_t)script_frames(mol);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(compiles(ir, "lo = distance_min(residue(1), residue(2)); hi = distance_max(residue(1), residue(2));", mol));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(eval != NULL);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

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
    const uint32_t num_frames = (uint32_t)script_frames(mol);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(compiles(ir,
        "all4 = contact_count(residue(1:4), residue(5:8), 5.0);"
        "c1 = contact_count(residue(1), residue(5:8), 5.0);"
        "c2 = contact_count(residue(2), residue(5:8), 5.0);"
        "c3 = contact_count(residue(3), residue(5:8), 5.0);"
        "c4 = contact_count(residue(4), residue(5:8), 5.0);", mol));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(eval != NULL);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

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
        md_script_eval_t* eval = md_script_eval_create((uint32_t)script_frames(mol), ir, alloc);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 1));
        const float* nb = eval_property_data(eval, "nb");
        ASSERT_TRUE(nb != NULL);
        EXPECT_EQ(10.0f, nb[0]); // residue 2 of the peptide has ten atoms
        md_script_eval_free(eval);
    }

    // Vectors
    ASSERT_TRUE(compiles(ir, "{x, y} = vec2(1, 2); z = x + y + 0.0 * distance(1, 2);", mol));
    {
        md_script_eval_t* eval = md_script_eval_create((uint32_t)script_frames(mol), ir, alloc);
        ASSERT_TRUE(eval != NULL);
        ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 1));
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
    const uint32_t num_frames = (uint32_t)script_frames(mol);
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
            EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
            md_script_eval_free(eval);
        }
    }

    // A comparison of arrays of different lengths is still an error
    EXPECT_FALSE(compiles(ir, "a = distance(1, 2) in residue(1:3); b = distance(1, 2) in residue(1:4); c = a < b;", mol));
}

// The script reference (docs/script_reference.md) gives every procedure an entry with a machine readable comment,
// `<!-- proc name=... aliases=a,b ... -->`, which is what the in-app help navigates by. Keep the two in step: every
// procedure the parser knows is documented, and every name the reference documents is one the parser knows.
static str_t reference_meta_value(str_t line, str_t key) {
    for (size_t p = 0; p + key.len < line.len; ++p) {
        if ((p == 0 || line.ptr[p - 1] == ' ') && str_eq(str_substr(line, p, key.len), key) && line.ptr[p + key.len] == '=') {
            const size_t beg = p + key.len + 1;
            size_t end = beg;
            while (end < line.len && line.ptr[end] != ' ') ++end;
            return str_substr(line, beg, end - beg);
        }
    }
    return (str_t){0};
}

static bool reference_contains(const str_t* names, size_t count, str_t name) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq(names[i], name)) return true;
    }
    return false;
}

static bool procedure_table_contains(const procedure_t* table, size_t count, str_t name) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq(table[i].name, name)) return true;
    }
    return false;
}

UTEST(script, reference_documents_every_procedure) {
    md_allocator_i* alloc = md_get_heap_allocator();
    str_t doc = load_textfile(STR_LIT(MD_UNITTEST_DATA_DIR "/../docs/script_reference.md"), alloc);
    ASSERT_FALSE(str_empty(doc));

    // Procedures that the parser handles itself rather than through the procedure table
    static const str_t intrinsics[] = { STR_INIT("attr"), STR_INIT("flatten"), STR_INIT("transpose") };

    str_t names[1024];  // documented names and aliases
    size_t num_names = 0;
    char msg[256];

    str_t rest = doc;
    str_t line;
    while (str_extract_line(&line, &rest)) {
        if (!str_begins_with(line, STR_LIT("<!-- proc"))) continue;
        const str_t name = reference_meta_value(line, STR_LIT("name"));
        snprintf(msg, sizeof(msg), "procedure comment without a name: '%.*s'", STR_ARG(line));
        EXPECT_FALSE_MSG(str_empty(name), msg);
        ASSERT_LT(num_names, ARRAY_SIZE(names));
        names[num_names++] = name;

        str_t list = reference_meta_value(line, STR_LIT("aliases"));
        while (!str_empty(list)) {
            size_t comma = list.len;
            str_find_char(&comma, list, ',');
            const str_t alias = str_trim(str_substr(list, 0, comma));
            if (!str_empty(alias)) {
                ASSERT_LT(num_names, ARRAY_SIZE(names));
                names[num_names++] = alias;
            }
            list = str_substr(list, MIN(comma + 1, list.len), SIZE_MAX);
        }
    }
    EXPECT_GT(num_names, (size_t)0);

    for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        if (procedure_table_contains(procedures, i, procedures[i].name)) continue;  // another overload, already checked
        snprintf(msg, sizeof(msg), "procedure '%.*s' has no entry in docs/script_reference.md", STR_ARG(procedures[i].name));
        EXPECT_TRUE_MSG(reference_contains(names, num_names, procedures[i].name), msg);
    }
    for (size_t i = 0; i < ARRAY_SIZE(intrinsics); ++i) {
        snprintf(msg, sizeof(msg), "procedure '%.*s' has no entry in docs/script_reference.md", STR_ARG(intrinsics[i]));
        EXPECT_TRUE_MSG(reference_contains(names, num_names, intrinsics[i]), msg);
    }
    for (size_t i = 0; i < num_names; ++i) {
        const bool known = procedure_table_contains(procedures, ARRAY_SIZE(procedures), names[i]) ||
                           procedure_table_contains(operators,  ARRAY_SIZE(operators),  names[i]) ||
                           procedure_table_contains(casts,      ARRAY_SIZE(casts),      names[i]) ||
                           reference_contains(intrinsics, ARRAY_SIZE(intrinsics), names[i]);
        snprintf(msg, sizeof(msg), "docs/script_reference.md documents '%.*s', which the parser does not know", STR_ARG(names[i]));
        EXPECT_TRUE_MSG(known, msg);
    }

    str_free(doc, alloc);
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
    const size_t num_frames = script_frames(sys);
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);

    test_hook_proc_eval.proc_ptr = proc_ptr;
    test_hook_proc_eval.count = 0;
    const bool ok = md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, (uint32_t)num_frames);
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
    const size_t num_frames = script_frames(sys);
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
    const size_t num_frames = script_frames(sys);

    static const char* srcs[] = {
        "{x, y, z} = com(residue(1)); a = x * 2;",
        "{x, y, z} = coord(1); a = x + y + z;",
        "{x, y, z} = com(residue(1)); a = x in residue(1:3);",
    };
    for (size_t i = 0; i < ARRAY_SIZE(srcs); ++i) {
        md_script_ir_t* ir = md_script_ir_create(alloc);
        EXPECT_TRUE_MSG(md_script_ir_compile_from_source(ir, str_from_cstr(srcs[i]), sys, NULL), srcs[i]);
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        EXPECT_TRUE_MSG(md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, (uint32_t)num_frames), srcs[i]);
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
    const size_t num_frames = script_frames(sys);

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
        EXPECT_FALSE(md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, (uint32_t)num_frames));
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
        "c3 = contacts(chunks(residue(:), 5), cutoff=5, min_separation=3, parent=residue(:), exclude_within=polymer(:));"
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
static md_bitset_t* contacts_residue_groups(size_t* out_num, const md_system_t* sys, md_allocator_i* alloc) {
    md_bitset_t* groups = md_alloc(alloc, sizeof(md_bitset_t) * sys->residue.count);
    for (size_t c = 0; c < sys->residue.count; ++c) {
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        groups[c] = md_bitset_range(range.beg, range.end);
    }
    *out_num = sys->residue.count;
    return groups;
}

// Per frame, what the script publishes against md_contact on the same frame
UTEST_F(script, contacts_matches_kernel) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    const size_t num_frames = script_frames(sys);

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
    md_bitset_t* res = contacts_residue_groups(&num, sys, alloc);
    ASSERT_EQ(md_attribute_element_count(&d->format), num_frames * num);

    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 6.0, .exclude_bonds = 3, .min_separation = 2 };
    md_contact_query_t q;
    ASSERT_TRUE(md_contact_query_init(&q, &desc, sys, alloc));
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_system_state_init(&state, sys->atom.count));

    float* deg = md_alloc(alloc, sizeof(float) * num);
    size_t total = 0;
    for (size_t f = 0; f < num_frames; ++f) {
        ASSERT_TRUE(run_extract_one(&state, sys, SCRIPT_RUN, (int64_t)f));
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
    const md_bitset_t* chunk = (const md_bitset_t*)data.ptr;
    const size_t num_chunks = element_count(data);

    size_t expected = 0;
    for (size_t c = 0; c < sys->residue.count; ++c) {
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        for (uint32_t beg = range.beg; beg < range.end; beg += 7) {
            const uint32_t end = MIN(beg + 7, range.end);
            ASSERT_LT(expected, num_chunks);
            EXPECT_EQ((size_t)(end - beg), md_bitset_count(chunk[expected]));
            EXPECT_EQ((size_t)(end - beg), md_bitset_count_range(chunk[expected], beg, end));
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
    md_bitset_t* res = contacts_residue_groups(&num_res, sys, alloc);

    // Chunks of 5 particles per residue, and for each the residue it is part of
    md_array(md_bitset_t) chunk = 0;
    md_array(uint32_t) parent = 0;
    uint32_t* label = md_alloc(alloc, sizeof(uint32_t) * sys->atom.count);
    for (size_t c = 0; c < num_res; ++c) {
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        for (uint32_t k = range.beg; k < range.end; ++k) label[k] = (uint32_t)c;
        for (uint32_t beg = range.beg; beg < range.end; beg += 5) {
            md_bitset_t bf = md_bitset_range(beg, MIN(beg + 5, range.end));
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
    const size_t num_frames = script_frames(sys);

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
    const str_t names[] = { STR_INIT("c"), STR_INIT("d"), STR_INIT("n") };
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        identifier_t* ident = get_identifier(ir, names[i]);
        ASSERT_TRUE(ident);
        md_script_vis_t vis = {0};
        md_script_vis_init(&vis, alloc);
        md_script_vis_ctx_t vctx = { .ir = ir, .sys = sys, .state = &sys->reference };
        EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, ident->node), -1, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT));
        // What is shown are the groups
        EXPECT_EQ(sys->atom.count, md_bitset_count(vis.atom_mask));
        md_script_vis_free(&vis);
    }
    EXPECT_EQ((size_t)0, test_hook_contact_searches);

    // Evaluating the trajectory searches once per frame
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, (uint32_t)num_frames));
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

#define ANY_KIND -1

// One group per component of the kind (or every component for ANY_KIND), as test_contact.c
static md_bitset_t* mirror_residue_groups(size_t* out_count, const md_system_t* sys, int kind, md_allocator_i* alloc) {
    md_bitset_t* groups = md_alloc(alloc, sizeof(md_bitset_t) * MAX(sys->residue.count, 1));
    size_t n = 0;
    for (size_t c = 0; c < sys->residue.count; ++c) {
        if (kind != ANY_KIND && (int)md_system_residue_kind(sys, c) != kind) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        groups[n] = md_bitset_range(range.beg, range.end);
        ++n;
    }
    *out_count = n;
    return groups;
}

static bool mirror_bitfield_equal(const md_bitset_t* a, const md_bitset_t* b, md_allocator_i* alloc) {
    (void)alloc;
    return md_bitset_equal(*a, *b);
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

static bool mirror_same_groups(const data_t* d, const md_bitset_t* groups, size_t num, const char* name, md_allocator_i* alloc) {
    if (d->type.base_type != TYPE_ATOMSET) {
        printf("  '%s' is not a selection\n", name);
        return false;
    }
    const size_t n = element_count(*d);
    if (n != num) {
        printf("  '%s' has %zu groups, the kernel case %zu\n", name, n, num);
        return false;
    }
    const md_bitset_t* bf = as_bitset(*d);
    for (size_t k = 0; k < n; ++k) {
        if (!mirror_bitfield_equal(&bf[k], &groups[k], alloc)) {
            printf("  group %zu of '%s' differs from the kernel case (%zu vs %zu particles)\n", k, name, md_bitset_count(bf[k]), md_bitset_count(groups[k]));
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

// "atom(b:e)" of the particles of residues [beg, end), one based and inclusive
static int mirror_atom_range(char* buf, size_t cap, const md_system_t* sys, size_t beg, size_t end) {
    const uint32_t a = md_system_residue_atom_range(sys, beg).beg;
    const uint32_t b = md_system_residue_atom_range(sys, end - 1).end;
    return snprintf(buf, cap, "atom(%u:%u)", a + 1, b);
}

UTEST_F(script, contacts_mirror_self_distance) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitset_t* res = mirror_residue_groups(&num, sys, ANY_KIND, alloc);
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
    md_bitset_t* res = mirror_residue_groups(&num, sys, ANY_KIND, alloc);
    ASSERT_GT(num, (size_t)10);

    // Overlapping sets: the first two thirds of the residues against the last two thirds
    const size_t third = num / 3;
    md_contact_desc_t desc = { .group_a = res, .num_a = num - third, .group_b = res + third, .num_b = num - third, .cutoff = 4.0 };
    char src[512];
    snprintf(src, sizeof(src), "ga = residue(1:%zu); gb = residue(%zu:%zu); c = contacts(ga, gb, cutoff=4.0, exclude_bonds=0);", num - third, third + 1, num);
    EXPECT_TRUE(mirror_check(sys, src, &desc, true, NULL, alloc));

    // Protein against everything
    size_t num_prot = 0;
    md_bitset_t* prot = mirror_residue_groups(&num_prot, sys, MD_RESIDUE_KIND_AMINO_ACID, alloc);
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
    md_bitset_t* prot = mirror_residue_groups(&num, sys, MD_RESIDUE_KIND_AMINO_ACID, alloc);
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
    md_bitset_t* all = mirror_residue_groups(&num_all, sys, ANY_KIND, alloc);
    ASSERT_GT(num_all, (size_t)10);
    // Every third residue: still spread over the whole cell, at a ninth of the cost of the reference
    md_bitset_t* res = md_alloc(alloc, sizeof(md_bitset_t) * num_all);
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
    md_bitset_t* res = mirror_residue_groups(&num, sys, MD_RESIDUE_KIND_AMINO_ACID, alloc);
    ASSERT_GT(num, (size_t)6);
    // The amino acids are the leading residues, so their indices are those of the residues
    for (size_t k = 0; k < num; ++k) ASSERT_EQ(MD_RESIDUE_KIND_AMINO_ACID, md_system_residue_kind(sys, k));

    // Windows of three consecutive residues, each overlapping the next by two. In the script: an array of selections
    const size_t num_win = num - 2;
    md_bitset_t* win = md_alloc(alloc, sizeof(md_bitset_t) * num_win);
    char windows[4096];
    size_t len = snprintf(windows, sizeof(windows), "{");
    for (size_t w = 0; w < num_win; ++w) {
        win[w] = md_bitset_or(md_bitset_or(res[w], res[w + 1], alloc), res[w + 2], alloc);
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
    md_bitset_t* res = mirror_residue_groups(&num, sys, MD_RESIDUE_KIND_AMINO_ACID, alloc);
    ASSERT_GT(num, (size_t)8);
    for (size_t k = 0; k < num; ++k) ASSERT_EQ(MD_RESIDUE_KIND_AMINO_ACID, md_system_residue_kind(sys, k));

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
    md_bitset_t g[2] = { md_bitset_range(0, 1), md_bitset_range(1, 2) };
    md_contact_desc_t desc = { .group_a = g, .num_a = 2, .cutoff = 0.01 };
    md_contact_set_t set = {0};
    EXPECT_TRUE(mirror_check(sys, "ga = {atom(1), atom(2)}; c = contacts(ga, cutoff=0.01, exclude_bonds=0);", &desc, false, &set, alloc));
    EXPECT_EQ(set.count, (size_t)0);

    // An empty B is an empty set: no contacts. Not the contacts within A, which is what no B means
    size_t num = 0;
    md_bitset_t* res = mirror_residue_groups(&num, sys, ANY_KIND, alloc);
    static const md_bitset_t no_groups = {0};
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
    job->ok = md_script_eval_frame_range(job->eval, job->ir, job->sys, SCRIPT_RUN, job->beg, job->end);
}

UTEST_F(script, contacts_mirror_query_over_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* sys = &utest_fixture->ala;
    ASSERT_GT(script_frames(sys), 0u);
    const size_t num_frames = script_frames(sys);
    ASSERT_GT(num_frames, (size_t)3);

    size_t num = 0;
    md_bitset_t* prot = mirror_residue_groups(&num, sys, MD_RESIDUE_KIND_AMINO_ACID, alloc);
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
        ASSERT_TRUE(run_extract_one(&state, sys, SCRIPT_RUN, (int64_t)f));
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

// The names that editors highlight: every keyword is one to the tokenizer, and the built-in identifiers are exactly
// the procedures, the parser's own procedures and the constants, each once.
UTEST(script, keywords_and_builtin_identifiers) {
    const size_t num_keywords = md_script_num_keywords();
    const str_t* keywords = md_script_keywords();
    ASSERT_GT(num_keywords, (size_t)0);
    char msg[128];
    for (size_t i = 0; i < num_keywords; ++i) {
        tokenizer_t tokenizer = tokenizer_init(keywords[i]);
        const token_t token = tokenizer_consume_next(&tokenizer);
        snprintf(msg, sizeof(msg), "keyword '%.*s'", STR_ARG(keywords[i]));
        EXPECT_TRUE_MSG(token.type != TOKEN_IDENT && str_eq(token.str, keywords[i]), msg);
    }

    const size_t count = md_script_builtin_identifiers(NULL, 0);
    ASSERT_GT(count, (size_t)0);
    str_t names[512];
    ASSERT_LE(count, ARRAY_SIZE(names));
    EXPECT_EQ(md_script_builtin_identifiers(names, ARRAY_SIZE(names)), count);

    // A smaller buffer gets the first names, and the same count
    str_t few[4];
    EXPECT_EQ(md_script_builtin_identifiers(few, ARRAY_SIZE(few)), count);
    for (size_t i = 0; i < ARRAY_SIZE(few); ++i) EXPECT_TRUE(str_eq(few[i], names[i]));

    bool has_attr = false, has_transpose = false, has_pi = false;
    for (size_t i = 0; i < count; ++i) {
        snprintf(msg, sizeof(msg), "builtin identifier '%.*s'", STR_ARG(names[i]));
        EXPECT_TRUE_MSG(md_script_identifier_name_valid(names[i]), msg);
        for (size_t j = i + 1; j < count; ++j) EXPECT_FALSE_MSG(str_eq(names[i], names[j]), msg);
        for (size_t k = 0; k < num_keywords; ++k) EXPECT_FALSE_MSG(str_eq(names[i], keywords[k]), msg);
        has_attr      |= str_eq(names[i], STR_LIT("attr"));
        has_transpose |= str_eq(names[i], STR_LIT("transpose"));
        has_pi        |= str_eq(names[i], STR_LIT("PI"));
    }
    EXPECT_TRUE(has_attr);
    EXPECT_TRUE(has_transpose);
    EXPECT_TRUE(has_pi);

    for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        bool found = false;
        for (size_t j = 0; j < count && !found; ++j) found = str_eq(procedures[i].name, names[j]);
        snprintf(msg, sizeof(msg), "procedure '%.*s'", STR_ARG(procedures[i].name));
        EXPECT_TRUE_MSG(found, msg);
    }
}

// PERIODIC INVARIANCE
//
// Translating every atom by the same vector and wrapping each atom back into the cell does not change
// the configuration, so it must not change anything a script measures. The translations are chosen so
// that the probed structure ends up straddling each face of the cell in turn, which is exactly the
// case where a procedure that treats the raw coordinates as contiguous goes wrong.

typedef struct pbc_probe_t {
    const char* name;
    int         kind;   // 0: scalar, 1: float[N], 2: volume
    float       tol;    // absolute for scalars, relative L1 for volumes
} pbc_probe_t;

static void pbc_translate(md_system_state_t* dst, const md_system_state_t* src, vec3_t shift, bool wrap) {
    for (size_t i = 0; i < src->num_atoms; ++i) {
        dst->xyz[i] = vec3_add(src->xyz[i], shift);
    }
    if (wrap) {
        md_util_pbc(dst->xyz, NULL, dst->num_atoms, &dst->unitcell);
    }
}

// As an evaluation of a frame does it, including the per atom masses and radii
static bool pbc_eval(data_t* out, md_script_ir_t* ir, const char* name, md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    identifier_t* ident = get_identifier(ir, str_from_cstr(name));
    if (!ident || !ident->node) return false;
    float* mass   = md_alloc(alloc, sizeof(float) * ALIGN_TO(sys->atom.count, 16));
    float* radius = md_alloc(alloc, sizeof(float) * ALIGN_TO(sys->atom.count, 16));
    md_atom_extract_masses(mass,   0, sys->atom.count, &sys->atom);
    md_atom_extract_radii (radius, 0, sys->atom.count, &sys->atom);
    eval_context_t ctx = { .ir = ir, .sys = sys, .atom_mass = mass, .atom_radius = radius, .temp_alloc = alloc, .alloc = alloc, .cur_state = state, .ref_state = &sys->reference };
    return evaluate_node_alloc(out, ident->node, &ctx, alloc);
}

static void pbc_invariance(int* utest_result, md_system_t* sys, const char* src, const char* anchor, const pbc_probe_t* probes, size_t num_probes, md_allocator_i* alloc) {
    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL);
    if (!md_script_ir_valid(ir)) {
        for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
            str_t err = md_script_ir_errors(ir)[i].text;
            printf("  %.*s\n", (int)err.len, err.ptr);
        }
    }
    ASSERT_TRUE(md_script_ir_valid(ir));

    const md_system_state_t* ref = &sys->reference;

    // The structure the translations are placed around: its plain centroid in the reference state
    md_bitset_t anchor_bf = {0};
    ASSERT_TRUE(eval_selection(&anchor_bf, str_from_cstr(anchor), sys, alloc));
    vec3_t c = {0};
    size_t n = 0;
    md_bitset_iter_t it = md_bitset_iter(anchor_bf);
    while (md_bitset_iter_next(&it)) {
        c = vec3_add(c, ref->xyz[it.idx]);
        n += 1;
    }
    ASSERT_GT(n, (size_t)0);
    c = vec3_div1(c, (float)n);

    mat3_t A = {0}, I = {0};
    md_unitcell_A_extract_float(A.elem, &ref->unitcell);
    md_unitcell_I_extract_float(I.elem, &ref->unitcell);
    const vec3_t fc = mat3_mul_vec3(I, c);

    // Values in the reference state
    data_t* want = md_alloc(alloc, sizeof(data_t) * num_probes);
    for (size_t p = 0; p < num_probes; ++p) {
        ASSERT_TRUE(pbc_eval(&want[p], ir, probes[p].name, sys, ref, alloc));
    }

    md_system_state_t state = *ref;
    state.xyz = md_alloc(alloc, sizeof(vec3_t) * ref->num_atoms);

    // Wrapped per atom: the centroid onto each face (fractional 0 along one axis), onto a corner,
    // and a generic offset. Not wrapped: the whole system carried out of the cell, as a trajectory
    // with whole molecules or without jumps has it.
    const struct { vec3_t f; bool wrap; } targets[] = {
        { { 0.0f,   fc.y,  fc.z }, true },
        { { fc.x,   0.0f,  fc.z }, true },
        { { fc.x,   fc.y,  0.0f }, true },
        { { 0.0f,   0.0f,  0.0f }, true },
        { { 1.0f,   0.02f, 0.97f }, true },
        { { fc.x + 0.37f, fc.y - 0.61f, fc.z + 0.29f }, true },
        { { fc.x + 2.0f,  fc.y, fc.z }, false },
        { { fc.x - 0.9f,  fc.y + 1.3f, fc.z - 2.2f }, false },
        { { 1.0f,   1.0f,  -1.0f }, false },
    };

    for (size_t t = 0; t < ARRAY_SIZE(targets); ++t) {
        const vec3_t shift = mat3_mul_vec3(A, vec3_sub(targets[t].f, fc));
        pbc_translate(&state, ref, shift, targets[t].wrap);

        for (size_t p = 0; p < num_probes; ++p) {
            data_t got = {0};
            ASSERT_TRUE(pbc_eval(&got, ir, probes[p].name, sys, &state, alloc));
            const float* a = (const float*)want[p].ptr;
            const float* b = (const float*)got.ptr;
            // Every probe is float data; volumes are one element of 128^3 floats, so count floats by size
            const size_t len = got.size / sizeof(float);
            ASSERT_EQ(want[p].size, got.size);

            if (probes[p].kind == 2) {
                double l1 = 0, sum = 0;
                for (size_t i = 0; i < len; ++i) {
                    l1  += fabs((double)a[i] - (double)b[i]);
                    sum += fabs((double)a[i]);
                }
                const double rel = sum > 0 ? l1 / sum : l1;
                if (!(rel <= probes[p].tol)) printf("  %s: target %zu, relative L1 difference %.4f (sum %.0f)\n", probes[p].name, t, rel, sum);
                EXPECT_GT(sum, 0.0);
                EXPECT_LE(rel, (double)probes[p].tol);
            } else {
                for (size_t i = 0; i < len; ++i) {
                    if (!(fabsf(a[i] - b[i]) <= probes[p].tol)) printf("  %s[%zu]: target %zu, reference %f, translated %f\n", probes[p].name, i, t, a[i], b[i]);
                    EXPECT_NEAR(a[i], b[i], probes[p].tol);
                }
            }
        }
    }
}

static const pbc_probe_t PBC_PROBES[] = {
    { "r",  0, 1.0e-3f },
    { "rf", 0, 1.0e-3f },
    { "a",  0, 1.0e-3f },
    { "d",  0, 1.0e-3f },
    { "x",  0, 1.0e-2f },
    { "w",  1, 1.0e-3f },
    { "v",  2, 2.0e-2f },
};

UTEST_F(script, pbc_invariance_ortho) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    pbc_invariance(utest_result, &utest_fixture->ala,
        "r = rmsd(residue(1:8));\n"
        "rf = rmsd(union(residue(1:8)));\n"
        "a = angle(atom(1), atom(45), atom(90));\n"
        "d = dihedral(atom(1), atom(30), atom(60), atom(90));\n"
        "x = distance(residue(1), residue(8));\n"
        "w = shape_weights(residue(1:8));\n"
        "v = sdf(residue(4), element('C'), 8.0);\n",
        "residue(4)", PBC_PROBES, ARRAY_SIZE(PBC_PROBES), alloc);
    md_vm_arena_destroy(alloc);
}

UTEST_F(script, pbc_invariance_triclinic) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    pbc_invariance(utest_result, &utest_fixture->npt,
        "r = rmsd(residue(1:8));\n"
        "rf = rmsd(union(residue(1:8)));\n"
        "a = angle(atom(1), atom(45), atom(90));\n"
        "d = dihedral(atom(1), atom(30), atom(60), atom(90));\n"
        "x = distance(residue(1), residue(8));\n"
        "w = shape_weights(residue(1:8));\n"
        "v = sdf(residue(4), element('O'), 8.0);\n",
        "residue(4)", PBC_PROBES, ARRAY_SIZE(PBC_PROBES), alloc);
    md_vm_arena_destroy(alloc);
}

// rmsd takes each bitfield of its argument as a structure of its own: one fit and one value per bitfield.
// union() pools them into a single structure.
UTEST_F(script, rmsd_per_structure) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    md_system_t* sys = &utest_fixture->ala;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, STR_LIT(
        "rr = rmsd(residue(1:8));\n"
        "rf = rmsd(union(residue(1:8)));\n"
        "r1 = rmsd(residue(1));\n"
        "r3 = rmsd(residue(3));\n"
        "r8 = rmsd(residue(8));\n"
        "rp = rmsd(residue(1:8) or residue(1:8));\n"
        "rc = rmsd(all()) in residue(1:8);\n"
        "rx = rmsd(residue(1:2)) in union(residue(3:4));\n"
        "ry = rmsd(residue(1:8)) in residue(3);\n"), sys, NULL);
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        str_t err = md_script_ir_errors(ir)[i].text;
        printf("  %.*s\n", (int)err.len, err.ptr);
    }
    ASSERT_TRUE(md_script_ir_valid(ir));

    // Every atom displaced by a pseudo random amount, so that each residue deviates by something different
    const md_system_state_t* ref = &sys->reference;
    md_system_state_t state = *ref;
    state.xyz = md_alloc(alloc, sizeof(vec3_t) * ref->num_atoms);
    srand(1234);
    for (size_t i = 0; i < ref->num_atoms; ++i) {
        const float s = 0.2f + 0.05f * (float)(i % 13);
        const vec3_t d = { s * ((float)rand() / RAND_MAX - 0.5f), s * ((float)rand() / RAND_MAX - 0.5f), s * ((float)rand() / RAND_MAX - 0.5f) };
        state.xyz[i] = vec3_add(ref->xyz[i], d);
    }

    data_t rr = {0}, rf = {0}, r1 = {0}, r3 = {0}, r8 = {0}, rc = {0}, rx = {0};
    ASSERT_TRUE(pbc_eval(&rr, ir, "rr", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&rf, ir, "rf", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&r1, ir, "r1", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&r3, ir, "r3", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&r8, ir, "r8", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&rc, ir, "rc", sys, &state, alloc));
    ASSERT_TRUE(pbc_eval(&rx, ir, "rx", sys, &state, alloc));

    ASSERT_EQ((size_t)8, rr.size / sizeof(float));
    ASSERT_EQ((size_t)1, rf.size / sizeof(float));
    ASSERT_EQ((size_t)1, r1.size / sizeof(float));

    const float* v = (const float*)rr.ptr;
    EXPECT_GT(v[0], 0.01f);
    EXPECT_NEAR(v[0], as_float(r1), 1.0e-5f);
    EXPECT_NEAR(v[2], as_float(r3), 1.0e-5f);
    EXPECT_NEAR(v[7], as_float(r8), 1.0e-5f);

    // The pooled fit is a different quantity from any of the separate ones: one rigid fit of all eight
    // residues cannot do better than eight separate fits, so it deviates at least as much as their mean square
    double ms = 0;
    for (int i = 0; i < 8; ++i) ms += (double)v[i] * v[i];
    EXPECT_GE(as_float(rf) + 1.0e-4f, (float)sqrt(ms / 8.0) * 0.9f);
    EXPECT_NE(as_float(rf), v[0]);

    // A selection combined with 'or' is still one bitfield
    data_t rp = {0};
    ASSERT_TRUE(pbc_eval(&rp, ir, "rp", sys, &state, alloc));
    EXPECT_EQ((size_t)1, rp.size / sizeof(float));
    EXPECT_NEAR(as_float(rp), as_float(rf), 1.0e-5f);

    // In a context the atoms are restricted to the context: all() in residue(i) is residue(i)
    ASSERT_EQ((size_t)8, rc.size / sizeof(float));
    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(((const float*)rc.ptr)[i], v[i], 1.0e-5f);
    }
    // An array argument inside a context gives [context][N]: residue indices count the residues of the context,
    // so residue(1:2) in union(residue(3:4)) is residues 3 and 4
    EXPECT_EQ((size_t)2, rx.size / sizeof(float));
    EXPECT_NEAR(((const float*)rx.ptr)[0], v[2], 1.0e-5f);
    EXPECT_NEAR(((const float*)rx.ptr)[1], v[3], 1.0e-5f);
    // and residue(3) holds one residue: residue(2:8) are past its last and select nothing
    data_t ry = {0};
    ASSERT_TRUE(pbc_eval(&ry, ir, "ry", sys, &state, alloc));
    EXPECT_EQ((size_t)1, ry.size / sizeof(float));
    EXPECT_NEAR(((const float*)ry.ptr)[0], v[2], 1.0e-5f);

    // Over the trajectory, as a property: one value per residue per frame
    const uint32_t num_frames = script_frames(sys);
    ASSERT_GT(num_frames, 0u);
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    EXPECT_TRUE(md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, num_frames));
    md_script_eval_free(eval);

    md_vm_arena_destroy(alloc);
}

static bool same_atoms(const md_bitset_t* a, const md_bitset_t* b, md_allocator_i* alloc) {
    (void)alloc;
    return md_bitset_equal(*a, *b);
}

// The closest (farthest) approach of the atoms of selection a to the atoms of selection b, under the minimum image
static float extent_reference(md_system_t* sys, const md_system_state_t* state, const char* a, const char* b, bool largest, md_allocator_i* alloc) {
    md_bitset_t bf[2] = {0};
    if (!eval_selection(&bf[0], str_from_cstr(a), sys, alloc) || !eval_selection(&bf[1], str_from_cstr(b), sys, alloc)) return -1.0f;
    float ext = largest ? 0.0f : FLT_MAX;
    md_bitset_iter_t it_a = md_bitset_iter(bf[0]);
    while (md_bitset_iter_next(&it_a)) {
        md_bitset_iter_t it_b = md_bitset_iter(bf[1]);
        while (md_bitset_iter_next(&it_b)) {
            vec3_t dx = vec3_sub(state->xyz[it_a.idx], state->xyz[it_b.idx]);
            md_util_min_image_vec3(&dx, 1, &state->unitcell);
            const float d = vec3_length(dx);
            ext = largest ? MAX(ext, d) : MIN(ext, d);
        }
    }
    return ext;
}

// distance_min and distance_max give one value per element of their first argument: the closest (farthest) approach of
// the atoms of that element to everything in the second argument, which is taken as a whole.
UTEST_F(script, distance_min_max_per_element_of_the_first_argument) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(4));
    md_system_t* sys = &utest_fixture->ala;
    const md_system_state_t* state = &sys->reference;

    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_script_ir_compile_from_source(ir, STR_LIT(
        "mn = distance_min(residue(1:4), residue(6));\n"
        "mx = distance_max(residue(1:4), residue(6));\n"
        "mn3 = distance_min(residue(3), residue(6));\n"
        "mb = distance_min(residue(1), residue(5:8));\n"
        "mi = distance_min({1, 5, 9}, residue(6));\n"
        "mr = distance_max(1:9, residue(6));\n"
        "mp = distance_min(vec3(0, 0, 0), residue(6));\n"
        "mq = distance_min(residue(1:4), com(residue(6)));\n"
        "mc = distance_min(all(), out com(residue(6))) in residue(1:4);\n"
        "dyn = distance_min(within(5.0, residue(6)) and not residue(6), residue(6));\n"
        "dl = distance_min(residue(within(3.0, residue(6)) and not residue(6)), residue(6));\n"), sys, NULL);
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        str_t err = md_script_ir_errors(ir)[i].text;
        printf("  %.*s\n", (int)err.len, err.ptr);
    }
    ASSERT_TRUE(md_script_ir_valid(ir));

    data_t mn = {0}, mx = {0}, mn3 = {0}, mb = {0}, mi = {0}, mr = {0}, mp = {0}, mq = {0}, mc = {0}, dyn = {0}, dl = {0};
    ASSERT_TRUE(pbc_eval(&mn,  ir, "mn",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mx,  ir, "mx",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mn3, ir, "mn3", sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mb,  ir, "mb",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mi,  ir, "mi",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mr,  ir, "mr",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mp,  ir, "mp",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mq,  ir, "mq",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&mc,  ir, "mc",  sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&dyn, ir, "dyn", sys, state, alloc));
    ASSERT_TRUE(pbc_eval(&dl,  ir, "dl",  sys, state, alloc));

    // One value per residue of the first argument, measured from its atoms
    ASSERT_EQ((size_t)4, mn.size / sizeof(float));
    ASSERT_EQ((size_t)4, mx.size / sizeof(float));
    for (int k = 0; k < 4; ++k) {
        char a[32];
        snprintf(a, sizeof(a), "residue(%i)", k + 1);
        EXPECT_NEAR(((const float*)mn.ptr)[k], extent_reference(sys, state, a, "residue(6)", false, alloc), 1.0e-4f);
        EXPECT_NEAR(((const float*)mx.ptr)[k], extent_reference(sys, state, a, "residue(6)", true,  alloc), 1.0e-4f);
        EXPECT_GT(((const float*)mx.ptr)[k], ((const float*)mn.ptr)[k]);
    }
    // A single selection gives a single value
    ASSERT_EQ((size_t)1, mn3.size / sizeof(float));
    EXPECT_EQ(((const float*)mn.ptr)[2], as_float(mn3));

    // The second argument is one set of atoms, not one centre of mass per residue
    ASSERT_EQ((size_t)1, mb.size / sizeof(float));
    EXPECT_NEAR(as_float(mb), extent_reference(sys, state, "residue(1)", "residue(5:8)", false, alloc), 1.0e-4f);

    // Atom indices are one element each, a range is one element
    ASSERT_EQ((size_t)3, mi.size / sizeof(float));
    EXPECT_NEAR(((const float*)mi.ptr)[0], extent_reference(sys, state, "atom(1)", "residue(6)", false, alloc), 1.0e-4f);
    EXPECT_NEAR(((const float*)mi.ptr)[1], extent_reference(sys, state, "atom(5)", "residue(6)", false, alloc), 1.0e-4f);
    EXPECT_NEAR(((const float*)mi.ptr)[2], extent_reference(sys, state, "atom(9)", "residue(6)", false, alloc), 1.0e-4f);
    ASSERT_EQ((size_t)1, mr.size / sizeof(float));
    EXPECT_NEAR(as_float(mr), extent_reference(sys, state, "atom(1:9)", "residue(6)", true, alloc), 1.0e-4f);
    ASSERT_EQ((size_t)1, mp.size / sizeof(float));
    EXPECT_GT(as_float(mp), 0.0f);

    // In a context, one value per context
    ASSERT_EQ((size_t)4, mq.size / sizeof(float));
    ASSERT_EQ((size_t)4, mc.size / sizeof(float));
    for (int k = 0; k < 4; ++k) {
        EXPECT_GT(((const float*)mq.ptr)[k], 0.0f);
        EXPECT_NEAR(((const float*)mc.ptr)[k], ((const float*)mq.ptr)[k], 1.0e-5f);
    }

    // A selection which changes size from frame to frame is still one element
    ASSERT_EQ((size_t)1, dyn.size / sizeof(float));
    EXPECT_GT(as_float(dyn), 0.0f);
    EXPECT_LE(as_float(dyn), 5.0f + 1.0e-4f);

    // While a number of elements which changes from frame to frame is a result which does, as for rmsd
    md_bitset_t near = {0};
    ASSERT_TRUE(eval_selection(&near, STR_LIT("within(3.0, residue(6)) and not residue(6)"), sys, alloc));
    size_t num_near = 0;
    for (size_t i = 0; i < sys->residue.count; ++i) {
        md_urange_t range = md_system_residue_atom_range(sys, i);
        num_near += md_bitset_test_any_range(near, range.beg, range.end);
    }
    EXPECT_GT(num_near, (size_t)0);
    EXPECT_EQ(num_near, dl.size / sizeof(float));
    EXPECT_EQ(0u, (unsigned)md_script_ir_property_flags(ir, STR_LIT("dl")));

    // Over the trajectory, as properties: four values per frame, and one for the dynamic selection
    const uint32_t num_frames = script_frames(sys);
    ASSERT_GT(num_frames, 0u);
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, num_frames));
    const float* p_mn  = eval_property_data(eval, "mn");
    const float* p_mn3 = eval_property_data(eval, "mn3");
    const float* p_dyn = eval_property_data(eval, "dyn");
    ASSERT_TRUE(p_mn && p_mn3 && p_dyn);
    for (uint32_t f = 0; f < num_frames; ++f) {
        EXPECT_EQ(p_mn3[f], p_mn[f * 4 + 2]);
    }
    md_script_eval_free(eval);

    // Shown: the closest approach of each element, or of the hovered one only
    md_bitset_t res[2] = {0};
    ASSERT_TRUE(eval_selection(&res[0], STR_LIT("residue(1:4) or residue(6)"), sys, alloc));
    ASSERT_TRUE(eval_selection(&res[1], STR_LIT("residue(3) or residue(6)"), sys, alloc));
    identifier_t* ident = get_identifier(ir, STR_LIT("mn"));
    ASSERT_TRUE(ident);
    md_script_vis_ctx_t vctx = { .ir = ir, .sys = sys, .state = state };
    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, ident->node), -1, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_EQ((size_t)4, md_array_size(vis.text));
    EXPECT_TRUE(same_atoms(&vis.atom_mask, &res[0], alloc));
    md_script_vis_clear(&vis);
    EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, ident->node), 2, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    ASSERT_EQ((size_t)1, md_array_size(vis.text));
    EXPECT_NEAR((float)vis.text[0].value, ((const float*)mn.ptr)[2], 1.0e-5f);
    EXPECT_TRUE(same_atoms(&vis.atom_mask, &res[1], alloc));
    md_script_vis_free(&vis);

    md_vm_arena_destroy(alloc);
}

// ### COMPLETION ###

// Completions at the '|' in src, which is taken out
static md_script_completions_t complete_at(const char* src_with_cursor, const md_system_t* sys, md_allocator_i* alloc) {
    const char* bar = strchr(src_with_cursor, '|');
    const int cursor = bar ? (int)(bar - src_with_cursor) : (int)strlen(src_with_cursor);
    char buf[512];
    snprintf(buf, sizeof(buf), "%.*s%s", cursor, src_with_cursor, bar ? bar + 1 : "");
    return md_script_complete(str_copy_cstr(buf, alloc), cursor, sys, alloc);
}

static int count_kind(md_script_completions_t c, md_script_completion_kind_t kind) {
    int n = 0;
    for (size_t i = 0; i < c.count; ++i) n += c.items[i].kind == kind;
    return n;
}

static bool has_item(md_script_completions_t c, md_script_completion_kind_t kind, const char* text) {
    for (size_t i = 0; i < c.count; ++i) {
        if (c.items[i].kind == kind && str_eq(c.items[i].text, str_from_cstr(text))) return true;
    }
    return false;
}

UTEST(script, completion_values_in_strings) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
    md_script_completions_t c;

    // In a string that is not closed yet: the names, and the quote that closes it
    c = complete_at("x = resname(\"P|", &test_mol, alloc);
    EXPECT_EQ(3, (int)c.count);
    EXPECT_EQ(3, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "PFT\""));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "SOL\""));
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("P")));
    EXPECT_EQ(13, c.range.beg);
    EXPECT_EQ(14, c.range.end);

    // In a closed string: the names alone, replacing its contents
    c = complete_at("x = resname(\"L|YZ\");", &test_mol, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "LYS"));
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("L")));
    EXPECT_EQ(13, c.range.beg);
    EXPECT_EQ(16, c.range.end);

    // Single quotes, the other name selectors
    c = complete_at("x = element('|')", &test_mol, alloc);
    EXPECT_EQ(5, (int)c.count);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "He"));
    c = complete_at("x = name('|", &test_mol, alloc);
    EXPECT_EQ(5, (int)c.count);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "CA'"));

    // An element of an array literal
    c = complete_at("x = resname({\"SOL\", \"|", &test_mol, alloc);
    EXPECT_EQ(3, (int)c.count);

    // A parameter that is not the first, given by position and by name
    c = complete_at("x = count(all, \"|\")", &test_mol, alloc);
    EXPECT_EQ(5, (int)c.count);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "residue"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "instance"));
    c = complete_at("x = count(all, unit=\"|\")", &test_mol, alloc);
    EXPECT_EQ(5, (int)c.count);

    // Values that are not from the system are there without one, the others are not
    c = complete_at("x = count(all, \"|\")", NULL, alloc);
    EXPECT_EQ(5, (int)c.count);
    c = complete_at("x = resname(\"|\")", NULL, alloc);
    EXPECT_EQ(0, (int)c.count);

    // Strings that take no known values, comments and numbers get nothing
    c = complete_at("x = attr(\"edr/|", &test_mol, alloc);  // No such attributes
    EXPECT_EQ(0, (int)c.count);
    c = complete_at("x = resname(\"ALA\"); # res|", &test_mol, alloc);
    EXPECT_EQ(0, (int)c.count);
    c = complete_at("x = 12|", &test_mol, alloc);
    EXPECT_EQ(0, (int)c.count);
    c = complete_at("x = resname(\"ALA\"|", &test_mol, alloc);
    EXPECT_EQ(0, (int)c.count);

    md_vm_arena_destroy(alloc);
}

UTEST(script, completion_identifiers) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
    md_script_completions_t c;

    // In an expression: variables, procedures, constants and keywords, and the identifier is replaced as a whole
    c = complete_at("d1 = distance(1,2);\n{a, b} = shape_weights(all);\ny = di|st", &test_mol, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "d1"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "a"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "b"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "dist"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "transpose"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_CONSTANT, "PI"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "in"));
    EXPECT_EQ(3, count_kind(c, MD_SCRIPT_COMPLETION_VARIABLE));
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "y"));   // Not in its own definition
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("di")));
    EXPECT_EQ((int)(c.range.end - c.range.beg), 4);

    // Every procedure name once
    const size_t num_builtin = md_script_builtin_identifiers(NULL, 0);
    EXPECT_EQ((int)(num_builtin - ARRAY_SIZE(constants)), count_kind(c, MD_SCRIPT_COMPLETION_PROCEDURE));

    // A keyword being typed is an identifier being typed
    c = complete_at("x = resname(\"ALA\") in|", &test_mol, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "molecule"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "instance"));

    // Where a variable is named, only the names that are taken
    c = complete_at("d1 = distance(1,2);\nd|", &test_mol, alloc);
    EXPECT_EQ(1, (int)c.count);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "d1"));
    c = complete_at("abc = 1;\n{ab|c, d} = shape_weights(all);", &test_mol, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "abc"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "d"));  // Being written, in the same statement
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_PROCEDURE));

    md_vm_arena_destroy(alloc);
}

UTEST(script, completion_arguments) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
    md_script_completions_t c;

    // At the start of an argument: the values it takes (quoted), next to everything else
    c = complete_at("x = resname(|", &test_mol, alloc);
    EXPECT_EQ(3, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "\"PFT\""));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    c = complete_at("x = resname(P|", &test_mol, alloc);
    EXPECT_EQ(3, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));
    // Not further into the argument
    c = complete_at("x = resname(y + |", &test_mol, alloc);
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));

    // The parameters that can still be given by name
    c = complete_at("x = count(|", &test_mol, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "sel="));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "unit="));
    c = complete_at("x = count(all, u|", &test_mol, alloc);
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "sel="));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "unit="));
    EXPECT_EQ(5, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));
    c = complete_at("x = contacts(a, cutoff=3, |", &test_mol, alloc);
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "a="));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "cutoff="));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "b="));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "exclude_bonds="));
    // Not for the value of an argument given by name, or in a nested call
    c = complete_at("x = count(sel=|", &test_mol, alloc);
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_PARAMETER));
    c = complete_at("x = count(residue(|", &test_mol, alloc);
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_PARAMETER));
    EXPECT_EQ(3, count_kind(c, MD_SCRIPT_COMPLETION_VALUE));

    md_vm_arena_destroy(alloc);
}

// The completions of a filter at the '|' in expr, which is taken out
static md_script_completions_t complete_filter_at(const char* expr_with_cursor, const md_script_ir_t* ir, md_allocator_i* alloc) {
    const char* bar = strchr(expr_with_cursor, '|');
    const int cursor = bar ? (int)(bar - expr_with_cursor) : (int)strlen(expr_with_cursor);
    char buf[512];
    snprintf(buf, sizeof(buf), "%.*s%s", cursor, expr_with_cursor, bar ? bar + 1 : "");
    const md_filter_context_t ctx = {&test_mol, NULL, ir};
    return md_filter_complete(str_copy_cstr(buf, alloc), cursor, &ctx, alloc);
}

static const md_script_completion_t* find_item(md_script_completions_t c, md_script_completion_kind_t kind, const char* text) {
    for (size_t i = 0; i < c.count; ++i) {
        if (c.items[i].kind == kind && str_eq(c.items[i].text, str_from_cstr(text))) return &c.items[i];
    }
    return NULL;
}

// Offered, and of the type expected at the cursor
static bool fits_item(md_script_completions_t c, md_script_completion_kind_t kind, const char* text) {
    const md_script_completion_t* item = find_item(c, kind, text);
    return item && item->fits;
}

// In a call, what has the type of the parameter the argument is for fits; at the top of a script nothing in particular
UTEST(script, completion_fits_the_parameter) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
    md_script_completions_t c;

    c = complete_at("x = within(5, |", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "protein"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));   // Offered, a script is not a filter
    EXPECT_FALSE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "not"));

    c = complete_at("x = within(|", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_CONSTANT, "PI"));
    EXPECT_FALSE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));

    // Given by name, and through an array literal
    c = complete_at("x = count(sel=|", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "residue"));
    c = complete_at("x = within(5, {resname(\"SOL\"), |", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "element"));

    // The values an argument takes and the names of parameters always do
    c = complete_at("x = resname(|", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_VALUE, "\"PFT\""));
    c = complete_at("x = count(|", &test_mol, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PARAMETER, "unit="));

    // At the top of a script anything goes
    c = complete_at("x = |", &test_mol, alloc);
    ASSERT_GT((int)c.count, 0);
    int num_fit = 0;
    for (size_t i = 0; i < c.count; ++i) num_fit += c.items[i].fits;
    EXPECT_EQ(0, num_fit);

    md_vm_arena_destroy(alloc);
}

UTEST(script, completion_of_filters) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(64));
    md_script_completions_t c;

    // Outside of a call: what picks atoms, and the keyword that negates, but nothing that joins two
    c = complete_filter_at("|", NULL, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "within"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "not"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_CONSTANT, "PI"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "and"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "of"));
    EXPECT_EQ(0, count_kind(c, MD_SCRIPT_COMPLETION_VARIABLE));
    c = complete_filter_at("protein and |", NULL, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "water"));

    // After an operand, only what joins it to another, also while it is typed
    c = complete_filter_at("resname('PFT') |", NULL, alloc);
    EXPECT_EQ(4, (int)c.count);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "and"));
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "or"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_KEYWORD, "in"));
    c = complete_filter_at("protein a|", NULL, alloc);
    EXPECT_EQ(4, (int)c.count);
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("a")));

    // In a call everything, with what fits the parameter marked
    c = complete_filter_at("within(|", NULL, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));
    EXPECT_FALSE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));
    c = complete_filter_at("within(5, |", NULL, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "resname"));
    EXPECT_FALSE(fits_item(c, MD_SCRIPT_COMPLETION_PROCEDURE, "distance"));

    // The values of the system in strings, as in a script
    c = complete_filter_at("resname('P|", NULL, alloc);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "PFT'"));

    // The names of a compiled script, those which pick atoms where only that goes
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("sel = resname(\"PFT\");\nd = distance(1, 2);"), &test_mol, NULL));
    c = complete_filter_at("|", ir, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "sel"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "d"));
    c = complete_filter_at("within(|", ir, alloc);
    EXPECT_TRUE(fits_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "d"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "sel"));
    EXPECT_FALSE(fits_item(c, MD_SCRIPT_COMPLETION_VARIABLE, "sel"));
    md_script_ir_free(ir);

    md_vm_arena_destroy(alloc);
}

// Every procedure that takes a string says what the string is, and every entry is for a string parameter
UTEST(script, completion_domain_table) {
    char msg[128];
    for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
        for (size_t p = 0; p < procedures[i].num_args; ++p) {
            if (procedures[i].arg_type[p].base_type != TYPE_STRING) continue;
            snprintf(msg, sizeof(msg), "procedure '%.*s' parameter %d", STR_ARG(procedures[i].name), (int)p);
            EXPECT_TRUE_MSG(find_param_domain(procedures[i].name, (int)p) != DOMAIN_NONE, msg);
        }
    }
    for (size_t i = 0; i < ARRAY_SIZE(param_domains); ++i) {
        bool found = false;
        for (size_t j = 0; j < ARRAY_SIZE(script_intrinsics) && !found; ++j) {
            found = str_eq(script_intrinsics[j], param_domains[i].proc);  // Parsed, not in the procedure table
        }
        for (size_t j = 0; j < ARRAY_SIZE(procedures) && !found; ++j) {
            found = str_eq(procedures[j].name, param_domains[i].proc) && param_domains[i].param < procedures[j].num_args &&
                    procedures[j].arg_type[param_domains[i].param].base_type == TYPE_STRING;
        }
        snprintf(msg, sizeof(msg), "domain of '%.*s' parameter %d", STR_ARG(param_domains[i].proc), (int)param_domains[i].param);
        EXPECT_TRUE_MSG(found, msg);
    }
}

// The tokenizer stops at the end of what it is given, also when that is in the middle of a larger buffer: a comment
// that runs to the end, and a carriage return at the end of an unterminated string, used to look one past it.
UTEST(script, tokenizer_stays_within_its_string) {
    const char buf[] = "x = 1; # comment\nabc";
    const str_t str = {buf, 16};  // Ends inside the comment, the buffer goes on with a newline and an identifier
    tokenizer_t tokenizer = tokenizer_init(str);
    int types[8] = {0};
    int count = 0;
    token_t token;
    while (token = tokenizer_consume_next(&tokenizer), token.type != TOKEN_END && count < 8) {
        types[count++] = token.type;
    }
    EXPECT_EQ(4, count);
    EXPECT_EQ(TOKEN_IDENT, types[0]);
    EXPECT_EQ(';', types[3]);
    EXPECT_EQ(TOKEN_END, (int)token.type);
    EXPECT_EQ(16, token.beg);
    EXPECT_EQ(TOKEN_END, (int)tokenizer_peek_next(&tokenizer).type);

    const char buf2[] = "\"abc\r\n";
    tokenizer = tokenizer_init((str_t){buf2, 5});  // Ends at the carriage return
    token = tokenizer_consume_next(&tokenizer);
    EXPECT_EQ(TOKEN_UNDEF, (int)token.type);  // An unterminated string
    EXPECT_EQ(5, token.end);

    // A script that ends in a comment compiles
    md_allocator_i* arena = md_vm_arena_create(MEGABYTES(64));
    md_script_ir_t* ir = md_script_ir_create(arena);
    EXPECT_TRUE(compiles(ir, "x = 1; # comment", &test_mol));
    md_script_ir_free(ir);
    md_vm_arena_destroy(arena);
}

UTEST(script, completion_attribute_paths) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {.alloc = arena};
    md_system_state_t state = {.alloc = arena};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro")));
    ASSERT_TRUE(md_trr_system_publish_run(&sys, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), STR_LIT("run/t"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_attributes_t* t = &sys.attributes;
    const md_attribute_t* run_time = md_attributes_find(t, STR_LIT("run/t/time"));
    ASSERT_TRUE(run_time != NULL);
    const uint32_t F = run_time->format.shape[0];
    const double* frame_times = (const double*)run_time->data;
    double* values = md_alloc(arena, F * sizeof(double));
    for (uint32_t f = 0; f < F; ++f) values[f] = (double)f;
    ASSERT_NE(publish_f64(t, "run/t/obs/value", series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), values), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/t/static",    series_fmt(3, 1), MD_ATTRIBUTE_FLAG_NONE, md_unit_none(), values), MD_ATTRIBUTE_INVALID);

    md_script_ir_t* ir = md_script_ir_create(arena);
    char msg[256];
    char src[256];

    // Every path offered is one that attr() reads: each compiles, as it is offered
    md_script_completions_t c = complete_at("x = attr(\"|", &sys, arena);
    ASSERT_GT(c.count, (size_t)0);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "obs/value\""));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "time\""));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "static\""));
    for (size_t i = 0; i < c.count; ++i) {
        snprintf(src, sizeof(src), "x = attr(\"" STR_FMT ");", STR_ARG(c.items[i].text));
        snprintf(msg, sizeof(msg), "offered '%.*s'", STR_ARG(c.items[i].label));
        EXPECT_TRUE_MSG(compiles(ir, src, &sys), msg);
    }

    // Full paths, once that is what is being typed
    c = complete_at("x = attr(\"run|\")", &sys, arena);
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "run/t/obs/value"));
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "obs/value"));

    // md_script_attr_path gives the form completion offers, and nothing for what attr() cannot read
    EXPECT_TRUE(str_eq(md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/t/obs/value"))), STR_LIT("obs/value")));
    EXPECT_TRUE(str_eq(md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/t/time"))), STR_LIT("time")));
    EXPECT_TRUE(str_empty(md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/t/static")))));
    EXPECT_TRUE(str_empty(md_script_attr_path(t, NULL)));

    // In full where the relative path names something outside the runs, which attr() would find first
    ASSERT_NE(publish_f64(t, "run/t/obs/shadowed", series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), values), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "obs/shadowed",       series_fmt(3, 1), MD_ATTRIBUTE_FLAG_NONE,     md_unit_none(), values), MD_ATTRIBUTE_INVALID);
    {
        const str_t path = md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/t/obs/shadowed")));
        EXPECT_TRUE(str_eq(path, STR_LIT("run/t/obs/shadowed")));
        snprintf(src, sizeof(src), "x = attr(\"" STR_FMT "\");", STR_ARG(path));
        EXPECT_TRUE(compiles(ir, src, &sys));
        EXPECT_FALSE(compiles(ir, "x = attr(\"obs/shadowed\");", &sys));
    }

    // A second run with the same path makes the short form ambiguous: both are offered in full, and still compile
    ASSERT_NE(publish_f64(t, "run/u/time",      series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, run_time->unit, frame_times), MD_ATTRIBUTE_INVALID);
    ASSERT_NE(publish_f64(t, "run/u/obs/value", series_fmt(F, 1), MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), values), MD_ATTRIBUTE_INVALID);
    c = complete_at("x = attr(\"|\")", &sys, arena);
    EXPECT_FALSE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "obs/value"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "run/t/obs/value"));
    EXPECT_TRUE(has_item(c, MD_SCRIPT_COMPLETION_VALUE, "run/u/obs/value"));
    EXPECT_TRUE(str_eq(md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/t/obs/value"))), STR_LIT("run/t/obs/value")));
    EXPECT_TRUE(str_eq(md_script_attr_path(t, md_attributes_find(t, STR_LIT("run/u/obs/value"))), STR_LIT("run/u/obs/value")));
    for (size_t i = 0; i < c.count; ++i) {
        snprintf(src, sizeof(src), "x = attr(\"" STR_FMT "\");", STR_ARG(c.items[i].text));
        snprintf(msg, sizeof(msg), "offered '%.*s'", STR_ARG(c.items[i].label));
        EXPECT_TRUE_MSG(compiles(ir, src, &sys), msg);
    }

    md_script_ir_free(ir);
    md_system_free(&sys);
    md_vm_arena_destroy(arena);
}

// ---------------------------------------------------------------------------------------------------------------
// A compiled IR is read-only to evaluations and visualizations: one IR is shared by the evaluations of frame
// ranges running on worker threads, and by visualizations requested while they run.

// Everything an evaluation could write to in an IR: its own arrays, and the AST nodes and what they hold
static uint64_t ir_snapshot(const md_script_ir_t* ir) {
    uint64_t h = 0xcbf29ce484222325ULL;
#define SNAP_MIX(v) (h = (h ^ (uint64_t)(v)) * 0x100000001b3ULL)
    SNAP_MIX(md_array_size(ir->errors));
    SNAP_MIX(md_array_size(ir->warnings));
    SNAP_MIX(md_array_size(ir->identifiers));
    SNAP_MIX(md_array_size(ir->contact_queries));
    SNAP_MIX(md_array_size(ir->nodes));
    SNAP_MIX(ir->flags);
    for (size_t i = 0; i < md_array_size(ir->identifiers); ++i) {
        SNAP_MIX((uintptr_t)ir->identifiers[i].node);
        SNAP_MIX((uintptr_t)ir->identifiers[i].data);
    }
    for (size_t i = 0; i < md_array_size(ir->nodes); ++i) {
        const ast_node_t* node = ir->nodes[i];
        SNAP_MIX(node->type);
        SNAP_MIX(node->flags);
        SNAP_MIX(md_array_size(node->children));
        SNAP_MIX((uintptr_t)node->children);
        SNAP_MIX(node->data.type.base_type);
        for (int d = 0; d < MAX_NUM_DIMS; ++d) SNAP_MIX((uint32_t)node->data.type.dim[d]);
        SNAP_MIX((uintptr_t)node->data.ptr);
        SNAP_MIX(node->data.size);
        if (node->data.ptr && node->data.size) {
            if (node->data.type.base_type == TYPE_ATOMSET) {
                const md_bitset_t* bf = (const md_bitset_t*)node->data.ptr;
                for (size_t k = 0; k < node->data.size / sizeof(md_bitset_t); ++k) {
                    SNAP_MIX(md_bitset_hash64(bf[k], 0));
                }
            } else {
                const uint8_t* bytes = (const uint8_t*)node->data.ptr;
                for (size_t k = 0; k < node->data.size; ++k) SNAP_MIX(bytes[k]);
            }
        }
    }
#undef SNAP_MIX
    return h;
}

UTEST_F(script, evaluation_does_not_write_ir) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)script_frames(mol);
    ASSERT_GT(num_frames, 0u);

    // Contexts, subscripts of selections, identifiers and prepared contact queries
    str_t src = STR_LIT(
        "sel = residue(:);"
        "sub = sel[2:4];"
        "d   = distance(com(sub[1]), com(sub[2]));"
        "x   = coord(1) in residue(1:3);"
        "c   = contacts(residue(:), cutoff=4.5);"
        "n   = count(c);"
    );
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, src, mol, NULL));
    ASSERT_TRUE(md_script_ir_valid(ir));
    const uint64_t compiled = ir_snapshot(ir);

    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    EXPECT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
    EXPECT_EQ(compiled, ir_snapshot(ir));

    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    const md_script_vis_ctx_t vis_ctx = { .ir = ir, .sys = mol, .state = &mol->reference };
    const md_script_vis_token_t* tokens = md_script_ir_vis_tokens(ir);
    for (size_t i = 0; i < md_script_ir_num_vis_tokens(ir); ++i) {
        for (int subidx = -1; subidx < 2; ++subidx) {
            md_script_vis_clear(&vis);
            md_script_vis_eval_ref(&vis, tokens[i].ref, subidx, &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT);
        }
    }
    EXPECT_EQ(compiled, ir_snapshot(ir));

    // A visualization of an expression compiles one of its own, against the IR as context
    md_script_vis_clear(&vis);
    EXPECT_TRUE(md_script_vis_eval_string(&vis, STR_LIT("dd = distance(com(sub[1]), com(sub[2]));"), &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_EQ(compiled, ir_snapshot(ir));

    md_script_vis_free(&vis);
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// The diagnostics of an evaluation go to its own log, never into the IR. Only a compilation records them there.
UTEST_F(script, evaluation_diagnostics_go_to_its_log) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_script_ir_t* ir = create_ir(alloc);
    const token_t tok = {0};

    eval_log_t log = { .alloc = alloc };
    eval_context_t ctx = { .ir = ir, .log = &log };
    EVAL_LOG_WARNING(&ctx, tok, "a warning %d", 1);
    EVAL_LOG_ERROR(&ctx, tok, "first error");
    EVAL_LOG_ERROR(&ctx, tok, "second error");   // As in a compilation, only the first error is kept
    EXPECT_EQ(0, (int)md_array_size(ir->errors));
    EXPECT_EQ(0, (int)md_array_size(ir->warnings));
    ASSERT_EQ(1, (int)md_array_size(log.warnings));
    ASSERT_EQ(1, (int)md_array_size(log.errors));
    EXPECT_TRUE(str_eq(log.warnings[0].text, STR_LIT("a warning 1")));
    EXPECT_TRUE(str_eq(log.errors[0].text, STR_LIT("first error")));

    // Without a log they go to md_log, and still not into the IR
    ctx.log = NULL;
    EVAL_LOG_ERROR(&ctx, tok, "unlogged");
    EXPECT_EQ(0, (int)md_array_size(ir->errors));

    // A compilation records them in the IR, which is how they reach the editor
    ir->record_log = true;
    ctx.compile_ir = ir;
    EVAL_LOG_ERROR(&ctx, tok, "compile error");
    ASSERT_EQ(1, (int)md_array_size(ir->errors));
    EXPECT_TRUE(str_eq(ir->errors[0].text, STR_LIT("compile error")));

    md_arena_allocator_destroy(alloc);
}

// A visualization of an expression reports whether it produced one, and marks the atoms involved
UTEST_F(script, vis_eval_string_reports_success) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    const md_script_vis_ctx_t vis_ctx = { .sys = mol, .state = &mol->reference };

    EXPECT_TRUE(md_script_vis_eval_string(&vis, STR_LIT("d = distance(1, 2);"), &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_EQ(2, (int)md_bitset_count(vis.atom_mask));

    md_script_vis_clear(&vis);
    EXPECT_FALSE(md_script_vis_eval_string(&vis, STR_LIT("d = distance(1,;"), &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT));

    md_script_vis_free(&vis);
    md_arena_allocator_destroy(alloc);
}

// A visualization reference names a part of one compilation. It resolves in that IR only: not in another, not
// in one compiled later from the same source, and not in its own IR once that has been recompiled.
UTEST_F(script, vis_ref_resolves_in_its_own_compilation_only) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    const str_t src = STR_LIT(
        "sel = residue(1:4);"
        "d   = distance(com(sel[1]), com(sel[2]));"
        "x   = coord(1) in residue(1:3);"
    );

    md_script_ir_t* a = md_script_ir_create(alloc);
    md_script_ir_t* b = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(a, src, mol, NULL));
    ASSERT_TRUE(md_script_ir_compile_from_source(b, src, mol, NULL));

    const md_script_vis_ref_t none = {0};
    EXPECT_TRUE(md_script_vis_ref_empty(none));
    EXPECT_FALSE(md_script_vis_ref_valid(a, none));
    EXPECT_TRUE(md_script_vis_ref_empty(md_script_ir_property_vis_ref(a, STR_LIT("no_such_property"))));
    EXPECT_TRUE(md_script_vis_ref_empty(md_script_ir_property_vis_ref(NULL, STR_LIT("d"))));

    const md_script_vis_ref_t ref_a = md_script_ir_property_vis_ref(a, STR_LIT("d"));
    const md_script_vis_ref_t ref_b = md_script_ir_property_vis_ref(b, STR_LIT("d"));
    ASSERT_FALSE(md_script_vis_ref_empty(ref_a));
    ASSERT_FALSE(md_script_vis_ref_empty(ref_b));
    EXPECT_FALSE(md_script_vis_ref_equal(ref_a, ref_b));
    EXPECT_TRUE (md_script_vis_ref_valid(a, ref_a));
    EXPECT_FALSE(md_script_vis_ref_valid(b, ref_a));
    EXPECT_TRUE (md_script_vis_ref_valid(b, ref_b));
    EXPECT_FALSE(md_script_vis_ref_valid(a, ref_b));
    EXPECT_FALSE(md_script_vis_ref_valid(NULL, ref_a));

    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    md_script_vis_ctx_t ctx = { .ir = a, .sys = mol, .state = &mol->reference };
    // A distance is shown as geometry
    EXPECT_TRUE(md_script_vis_eval_ref(&vis, ref_a, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_GT(md_array_size(vis.points) + md_array_size(vis.lines), (size_t)0);
    md_script_vis_clear(&vis);
    EXPECT_FALSE(md_script_vis_eval_ref(&vis, ref_b, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_EQ((size_t)0, md_array_size(vis.points) + md_array_size(vis.lines));
    EXPECT_EQ((size_t)0, md_bitset_count(vis.atom_mask));

    // The editor's tokens resolve in the compilation they came from
    const md_script_vis_token_t* tokens = md_script_ir_vis_tokens(a);
    ASSERT_GT(md_script_ir_num_vis_tokens(a), (size_t)0);
    for (size_t i = 0; i < md_script_ir_num_vis_tokens(a); ++i) {
        EXPECT_TRUE (md_script_vis_ref_valid(a, tokens[i].ref));
        EXPECT_FALSE(md_script_vis_ref_valid(b, tokens[i].ref));
    }

    // Recompiled, its own IR does not take an old reference for one of its own
    md_script_ir_clear(a);
    ASSERT_TRUE(md_script_ir_compile_from_source(a, src, mol, NULL));
    EXPECT_FALSE(md_script_vis_ref_valid(a, ref_a));
    EXPECT_TRUE(str_empty(md_script_vis_ref_ident(a, ref_a)));
    EXPECT_EQ(0, md_script_vis_ref_dim(a, ref_a));
    EXPECT_TRUE(md_script_vis_ref_valid(a, md_script_ir_property_vis_ref(a, STR_LIT("d"))));

    // The array 'x' in residue(1:3) has one element per context
    const md_script_vis_ref_t ref_x = md_script_ir_property_vis_ref(a, STR_LIT("x"));
    if (!md_script_vis_ref_empty(ref_x)) {
        EXPECT_EQ(3, md_script_vis_ref_dim(a, ref_x));
    }

    md_script_vis_free(&vis);
    md_script_ir_free(a);
    md_script_ir_free(b);
    md_arena_allocator_destroy(alloc);
}

// A node is referred to by its position in the IR: that has to hold for every node, also those the parser and
// the static check make by copying another (subscripts, casts)
UTEST_F(script, node_index_is_position_in_ir) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "sel = residue(:);"
        "sub = sel[2:4];"
        "xyz = coord(1:8);"
        "xs  = xyz[1:2:7,1];"
        "f   = 1 + 2.5;"
        "d   = distance(com(sub[1]), com(sub[2])) * 2;"
    ), mol, NULL));
    for (size_t i = 0; i < md_array_size(ir->nodes); ++i) {
        EXPECT_EQ((uint32_t)i, ir->nodes[i]->index);
    }
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

static const eval_property_t* find_eval_property(const md_script_eval_t* eval, str_t ident) {
    for (size_t i = 0; i < md_array_size(eval->props); ++i) {
        if (str_eq(eval->props[i].ident, ident)) return &eval->props[i];
    }
    return NULL;
}

// A distribution is shown as values / weights, also with bins merged (sum of values / sum of weights). Over several
// frames that has to be sum(values) / sum(weights) over the frames evaluated so far: the weights of an rdf differ
// between frames (it normalizes by the pair count of its frame), so those of any one frame will not do.
UTEST_F(script, distribution_weights_accumulate) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)script_frames(mol);
    ASSERT_GT(num_frames, 2u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("r = rdf(element('C'), element('O'), 6.0);"), mol, NULL));
    ASSERT_EQ(MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION, md_script_ir_property_flags(ir, STR_LIT("r")));

    // Each frame on its own: its values (counts) and weights
    size_t num_bins = 0;
    double* sum_v = NULL;
    double* sum_w = NULL;
    float* first_w = NULL;
    bool weights_differ = false;
    for (uint32_t f = 0; f < num_frames; ++f) {
        md_script_eval_t* single = md_script_eval_create(num_frames, ir, alloc);
        ASSERT_TRUE(md_script_eval_frame_range(single, ir, mol, SCRIPT_RUN, f, f + 1));
        const eval_property_t* p = find_eval_property(single, STR_LIT("r"));
        ASSERT_TRUE(p);
        if (!sum_v) {
            num_bins = p->count;
            sum_v   = md_alloc(alloc, num_bins * sizeof(double));
            sum_w   = md_alloc(alloc, num_bins * sizeof(double));
            first_w = md_alloc(alloc, num_bins * sizeof(float));
            MEMSET(sum_v, 0, num_bins * sizeof(double));
            MEMSET(sum_w, 0, num_bins * sizeof(double));
            MEMCPY(first_w, p->weights, num_bins * sizeof(float));
        }
        for (size_t i = 0; i < num_bins; ++i) {
            sum_v[i] += p->values[i];
            sum_w[i] += p->weights[i];
            weights_differ |= fabsf(p->weights[i] - first_w[i]) > 1.0e-6f * MAX(1.0f, fabsf(first_w[i]));
        }
        md_script_eval_free(single);
    }
    // Otherwise the test shows nothing
    ASSERT_TRUE(weights_differ);

    // All frames, in two ranges: evaluated in either order, and as the frames come in
    for (int order = 0; order < 2; ++order) {
        md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
        const uint32_t mid = num_frames / 2;
        if (order == 0) {
            ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, mid));
            ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, mid, num_frames));
        } else {
            ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, mid, num_frames));
            ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, mid));
        }
        const eval_property_t* p = find_eval_property(eval, STR_LIT("r"));
        ASSERT_TRUE(p);
        ASSERT_EQ(num_bins, (size_t)p->count);

        int mismatches = 0;
        for (size_t i = 0; i < num_bins; ++i) {
            if (sum_w[i] <= 0.0) continue;
            const double expected = sum_v[i] / sum_w[i];
            const double shown    = p->weights[i] > 0.0f ? (double)p->values[i] / (double)p->weights[i] : 0.0;
            if (fabs(shown - expected) > 1.0e-4 * MAX(1.0, fabs(expected))) {
                mismatches += 1;
            }
        }
        EXPECT_EQ(0, mismatches);

        // Merged into a quarter of the bins, as the viewer resamples: still the ratio of the sums
        const size_t factor = 4;
        int merged_mismatches = 0;
        for (size_t d = 0; d < num_bins / factor; ++d) {
            double v = 0, w = 0, ev = 0, ew = 0;
            for (size_t k = 0; k < factor; ++k) {
                const size_t i = d * factor + k;
                v += p->values[i]; w += p->weights[i];
                ev += sum_v[i];    ew += sum_w[i];
            }
            if (w <= 0.0 || ew <= 0.0) continue;
            if (fabs(v / w - ev / ew) > 1.0e-4 * MAX(1.0, fabs(ev / ew))) {
                merged_mismatches += 1;
            }
        }
        EXPECT_EQ(0, merged_mismatches);

        // Cleared, nothing is shown, and the weights are back to their default
        md_script_eval_clear_data(eval);
        for (size_t i = 0; i < num_bins; ++i) {
            if (p->values[i] != 0.0f || p->weights[i] != 1.0f) { EXPECT_TRUE(false); break; }
        }
        md_script_eval_free(eval);
    }

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

static bool bits_equal(md_bitset_t a, md_bitset_t b) {
    return md_bitset_equal(a, b);
}

// The ratio values / weights of every bin of a distribution, as shown
static bool dist_same(const eval_property_t* a, const eval_property_t* b) {
    if (!a || !b || a->count != b->count) return false;
    for (size_t i = 0; i < a->count; ++i) {
        const double ra = a->weights[i] > 0 ? (double)a->values[i] / a->weights[i] : 0.0;
        const double rb = b->weights[i] > 0 ? (double)b->values[i] / b->weights[i] : 0.0;
        if (fabs(ra - rb) > 1.0e-4 * MAX(1.0, fabs(rb))) return false;
    }
    return true;
}

static bool floats_near(const float* a, const float* b, size_t n, float rel);

// The property evaluated on its own, over the given frames only (every frame for NULL): what it publishes is the
// aggregate over the frames evaluated
static md_script_eval_t* eval_over(const md_script_ir_t* ir, md_system_t* sys, uint32_t num_frames, const md_bitset_t* frames, md_allocator_i* alloc) {
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    if (!eval) return NULL;
    if (!frames) {
        md_script_eval_frame_range(eval, ir, sys, SCRIPT_RUN, 0, num_frames);
        return eval;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    const size_t count = md_bitset_count(*frames);
    uint32_t* indices = md_temp_alloc_array(temp, uint32_t, MAX(count, 1));
    md_bitset_extract_indices(indices, count, *frames);
    md_script_eval_frames(eval, ir, sys, SCRIPT_RUN, indices, count);
    md_temp_end(temp);
    return eval;
}

// What a slot holds, extracted into a property of the same shape, so that it compares like one (dist_same)
static eval_property_t extract_slot(md_script_eval_slot_info_t* out_info, const md_script_eval_t* eval, uint64_t key, str_t ident, md_allocator_i* alloc) {
    eval_property_t out = {0};
    md_script_eval_slot_info_t info = {0};
    if (!md_script_eval_slot_info(&info, eval, key, ident)) return out;
    out.count      = (uint32_t)info.num_values;
    out.num_values = info.num_values;
    out.values  = md_alloc(alloc, info.num_values * sizeof(float));
    out.weights = md_alloc(alloc, info.num_values * sizeof(float));
    if (!md_script_eval_slot_extract(out_info, out.values, out.weights, info.num_values, eval, key, ident)) {
        out.values = NULL;
    }
    return out;
}

static md_script_eval_slot_info_t slot_info_of(const md_script_eval_t* eval, uint64_t key, str_t ident) {
    md_script_eval_slot_info_t info = {0};
    md_script_eval_slot_info(&info, eval, key, ident);
    return info;
}

// A distribution kept per frame: a slot over any frames is summed from what is kept, and nothing is evaluated for it
UTEST_F(script, slots_of_a_distribution_kept_per_frame) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "d = distance(1, 10);"
        "r = rdf(element('C'), element('O'), 6.0);"
    ), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(n, ir, alloc);
    ASSERT_TRUE(eval);
    ASSERT_TRUE(find_eval_property(eval, STR_LIT("r"))->frame_values != NULL);

    // Only aggregates have slots
    const md_bitset_t sub = md_bitset_range(2, n - 2);
    ASSERT_TRUE (md_script_eval_slot_request(eval, 1, STR_LIT("r"), &sub));
    EXPECT_FALSE(md_script_eval_slot_request(eval, 1, STR_LIT("d"), &sub));
    EXPECT_FALSE(md_script_eval_slot_request(eval, 1, STR_LIT("none"), &sub));
    EXPECT_FALSE(md_script_eval_slot_info(&(md_script_eval_slot_info_t){0}, eval, 2, STR_LIT("r")));

    // Asked for before anything is evaluated, it holds nothing and waits for its frames
    md_script_eval_slot_info_t info = slot_info_of(eval, 1, STR_LIT("r"));
    EXPECT_NE((uint64_t)0, info.version);
    EXPECT_EQ(0u, info.num_frames);
    EXPECT_EQ(n - 4, info.num_pending);
    EXPECT_EQ((size_t)n, md_bitset_count(md_script_eval_pending_frames(eval, alloc)));

    // Each frame evaluated once, for the published aggregate and the slot at once
    test_hook_proc_eval.proc_ptr = _rdf_flt;
    test_hook_proc_eval.count = 0;
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n));
    EXPECT_EQ((size_t)n, test_hook_proc_eval.count);
    EXPECT_EQ((size_t)0, md_bitset_count(md_script_eval_pending_frames(eval, alloc)));

    eval_property_t held = extract_slot(&info, eval, 1, STR_LIT("r"), alloc);
    ASSERT_TRUE(held.values != NULL);
    EXPECT_EQ(n - 4, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);

    // The distribution over its frames alone, not over all of them
    md_script_eval_t* all  = eval_over(ir, mol, n, NULL, alloc);
    md_script_eval_t* part = eval_over(ir, mol, n, &sub, alloc);
    ASSERT_TRUE(all && part);
    EXPECT_TRUE (dist_same(&held, find_eval_property(part, STR_LIT("r"))));
    EXPECT_FALSE(dist_same(&held, find_eval_property(all,  STR_LIT("r"))));
    // and the published one is over every frame, as before
    EXPECT_TRUE (dist_same(find_eval_property(eval, STR_LIT("r")), find_eval_property(all, STR_LIT("r"))));

    // Moved: nothing is evaluated, it holds its new frames at once, under a new version
    const uint64_t before = info.version;
    const md_bitset_t moved = md_bitset_range(0, 3);
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &moved));
    EXPECT_EQ((size_t)0, md_bitset_count(md_script_eval_pending_frames(eval, alloc)));
    held = extract_slot(&info, eval, 1, STR_LIT("r"), alloc);
    EXPECT_EQ(3u, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);
    EXPECT_NE(before, info.version);
    md_script_eval_t* moved_alone = eval_over(ir, mol, n, &moved, alloc);
    EXPECT_TRUE(dist_same(&held, find_eval_property(moved_alone, STR_LIT("r"))));

    // Asked for the same frames again, nothing changes
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &moved));
    EXPECT_EQ(info.version, slot_info_of(eval, 1, STR_LIT("r")).version);

    // Another key is another slot, and one let go is gone
    ASSERT_TRUE(md_script_eval_slot_request(eval, 2, STR_LIT("r"), &sub));
    EXPECT_EQ(n - 4, slot_info_of(eval, 2, STR_LIT("r")).num_frames);
    md_script_eval_slot_release(eval, 2, STR_LIT("r"));
    EXPECT_FALSE(md_script_eval_slot_info(&info, eval, 2, STR_LIT("r")));
    EXPECT_EQ(3u, slot_info_of(eval, 1, STR_LIT("r")).num_frames);

    test_hook_proc_eval.proc_ptr = NULL;
    test_hook_proc_eval.count = 0;
    md_script_eval_free(moved_alone);
    md_script_eval_free(part);
    md_script_eval_free(all);
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// A distribution not kept per frame is added up as frames are evaluated: a slot moved evaluates the frames it no
// longer wants and those it now wants, nothing else, and holds exactly what it would over its frames alone
UTEST_F(script, slots_added_up_as_frames_are_evaluated) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("r = rdf(element('C'), element('O'), 6.0);"), mol, NULL));
    const md_script_eval_desc_t desc = { .num_frames = n, .frame_budget = 1 };
    md_script_eval_t* eval = md_script_eval_create_desc(ir, &desc, alloc);
    ASSERT_TRUE(eval);
    ASSERT_TRUE(find_eval_property(eval, STR_LIT("r"))->frame_values == NULL);

    const md_bitset_t sub = md_bitset_range(2, n - 2);
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &sub));

    test_hook_proc_eval.proc_ptr = _rdf_flt;
    test_hook_proc_eval.count = 0;
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n));
    EXPECT_EQ((size_t)n, test_hook_proc_eval.count);

    md_script_eval_slot_info_t info = {0};
    eval_property_t held = extract_slot(&info, eval, 1, STR_LIT("r"), alloc);
    EXPECT_EQ(n - 4, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);
    md_script_eval_t* part = eval_over(ir, mol, n, &sub, alloc);
    EXPECT_TRUE(dist_same(&held, find_eval_property(part, STR_LIT("r"))));

    // Moved: what remains is the frames to take out and those to add
    const md_bitset_t moved = md_bitset_range(0, 3);
    const md_bitset_t moving = md_bitset_xor(sub, moved, alloc);
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &moved));
    EXPECT_TRUE(bits_equal(moving, md_script_eval_pending_frames(eval, alloc)));
    info = slot_info_of(eval, 1, STR_LIT("r"));
    EXPECT_EQ(n - 4, info.num_frames);  // until they are evaluated, it holds what it held
    EXPECT_EQ((uint32_t)md_bitset_count(moving), info.num_pending);

    test_hook_proc_eval.count = 0;
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n));
    EXPECT_EQ((size_t)md_bitset_count(moving), test_hook_proc_eval.count);
    held = extract_slot(&info, eval, 1, STR_LIT("r"), alloc);
    EXPECT_EQ(3u, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);
    md_script_eval_t* moved_alone = eval_over(ir, mol, n, &moved, alloc);
    EXPECT_TRUE(dist_same(&held, find_eval_property(moved_alone, STR_LIT("r"))));

    // Asked for nothing, it holds nothing: values of zero and weights of one, exactly
    const md_bitset_t none = md_bitset_empty_set();
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &none));
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n));
    held = extract_slot(&info, eval, 1, STR_LIT("r"), alloc);
    EXPECT_EQ(0u, info.num_frames);
    bool empty = held.values != NULL;
    for (size_t i = 0; i < held.num_values && empty; ++i) {
        empty = held.values[i] == 0.0f && held.weights[i] == 1.0f;
    }
    EXPECT_TRUE(empty);

    // Interrupted, nothing is lost: what remains stays to be done
    ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &sub));
    md_script_eval_interrupt(eval);
    md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n);
    EXPECT_EQ(n - 4, slot_info_of(eval, 1, STR_LIT("r")).num_pending);
    md_script_eval_reset_interrupt(eval);

    test_hook_proc_eval.proc_ptr = NULL;
    test_hook_proc_eval.count = 0;
    md_script_eval_free(moved_alone);
    md_script_eval_free(part);
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

UTEST_F(script, slots_of_a_volume) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(64));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("v = sdf(residue(1), element('H'), 5.0);"), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(n, ir, alloc);
    ASSERT_TRUE(eval);

    const md_bitset_t a = md_bitset_range(0, 4);
    const md_bitset_t b = md_bitset_range(2, 6);
    ASSERT_TRUE(md_script_eval_slot_request(eval, 7, STR_LIT("v"), &a));
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 6));

    const size_t num_voxels = find_eval_property(eval, STR_LIT("v"))->num_values;
    float* values = md_alloc(alloc, num_voxels * sizeof(float));
    md_script_eval_slot_info_t info = {0};
    EXPECT_FALSE(md_script_eval_slot_extract(&info, values, NULL, num_voxels - 1, eval, 7, STR_LIT("v")));
    EXPECT_EQ(num_voxels, info.num_values);
    ASSERT_TRUE(md_script_eval_slot_extract(&info, values, NULL, num_voxels, eval, 7, STR_LIT("v")));
    EXPECT_EQ(4u, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);
    md_script_eval_t* over_a = eval_over(ir, mol, n, &a, alloc);
    EXPECT_TRUE(floats_near(values, find_eval_property(over_a, STR_LIT("v"))->values, num_voxels, 1e-4f));
    md_script_eval_free(over_a);

    // Moved: two frames out, two in
    ASSERT_TRUE(md_script_eval_slot_request(eval, 7, STR_LIT("v"), &b));
    EXPECT_EQ(4u, slot_info_of(eval, 7, STR_LIT("v")).num_pending);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 6));
    ASSERT_TRUE(md_script_eval_slot_extract(&info, values, NULL, num_voxels, eval, 7, STR_LIT("v")));
    EXPECT_EQ(4u, info.num_frames);
    EXPECT_EQ(0u, info.num_pending);
    md_script_eval_t* over_b = eval_over(ir, mol, n, &b, alloc);
    EXPECT_TRUE(floats_near(values, find_eval_property(over_b, STR_LIT("v"))->values, num_voxels, 1e-4f));
    md_script_eval_free(over_b);

    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// Recompiled, the slots of what did not change keep what they hold; those of what changed are kept as asked for
UTEST_F(script, slots_carry_over_to_the_next_evaluation) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);
    const md_bitset_t sub = md_bitset_range(2, n - 2);

    for (int kept = 0; kept < 2; ++kept) {
        const md_script_eval_desc_t desc = { .num_frames = n, .frame_budget = kept ? 0 : 1 };
        md_script_ir_t* ir1 = md_script_ir_create(alloc);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir1, STR_LIT(
            "r = rdf(element('C'), element('O'), 6.0);"
            "s = rdf(element('N'), element('O'), 6.0);"), mol, NULL));
        md_script_eval_t* eval1 = md_script_eval_create_desc(ir1, &desc, alloc);
        ASSERT_TRUE(md_script_eval_slot_request(eval1, 1, STR_LIT("r"), &sub));
        ASSERT_TRUE(md_script_eval_slot_request(eval1, 1, STR_LIT("s"), &sub));
        ASSERT_TRUE(md_script_eval_frame_range(eval1, ir1, mol, SCRIPT_RUN, 0, n));
        const uint64_t version = slot_info_of(eval1, 1, STR_LIT("r")).version;

        // r renamed, s computed from something else
        md_script_ir_t* ir2 = md_script_ir_create(alloc);
        ASSERT_TRUE(md_script_ir_compile_from_source(ir2, STR_LIT(
            "g = rdf(element('C'), element('O'), 6.0);"
            "s = rdf(element('N'), element('O'), 5.0);"), mol, NULL));
        md_script_eval_t* eval2 = md_script_eval_create_from(eval1, ir2, &desc, alloc);
        ASSERT_TRUE(eval2);

        md_script_eval_slot_info_t g = {0}, s = {0};
        ASSERT_TRUE(md_script_eval_slot_info(&g, eval2, 1, STR_LIT("g")));
        ASSERT_TRUE(md_script_eval_slot_info(&s, eval2, 1, STR_LIT("s")));
        EXPECT_EQ(n - 4, g.num_frames);
        EXPECT_EQ(0u, g.num_pending);
        EXPECT_GT(g.version, version);
        EXPECT_EQ(0u, s.num_frames);
        EXPECT_EQ(n - 4, s.num_pending);

        // What g holds is what r held
        eval_property_t r_held = extract_slot(NULL, eval1, 1, STR_LIT("r"), alloc);
        eval_property_t g_held = extract_slot(NULL, eval2, 1, STR_LIT("g"), alloc);
        EXPECT_TRUE(dist_same(&g_held, &r_held));

        md_script_eval_free(eval2);
        md_script_eval_free(eval1);
        md_script_ir_free(ir2);
        md_script_ir_free(ir1);
    }
    md_arena_allocator_destroy(alloc);
}

typedef struct eval_job_t {
    md_script_eval_t* eval;
    const md_script_ir_t* ir;
    md_system_t* sys;
    uint32_t num_frames;
    volatile uint64_t done;
} eval_job_t;

static void eval_job_run(void* data) {
    eval_job_t* job = (eval_job_t*)data;
    md_script_eval_frame_range(job->eval, job->ir, job->sys, SCRIPT_RUN, 0, job->num_frames);
    md_atomic_fetch_add64(&job->done, 1);
}

// Read while another thread evaluates: whatever it holds when it is read, it is exactly the aggregate over the frames
// it holds. Evaluated in order by one thread, those are the first of its frames.
UTEST_F(script, slots_read_while_evaluating) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("r = rdf(element('C'), element('O'), 6.0);"), mol, NULL));

    for (int kept = 0; kept < 2; ++kept) {
        // What each frame holds on its own
        md_script_eval_t* ref = md_script_eval_create(n, ir, alloc);
        ASSERT_TRUE(md_script_eval_frame_range(ref, ir, mol, SCRIPT_RUN, 0, n));
        const eval_property_t* per_frame = find_eval_property(ref, STR_LIT("r"));
        ASSERT_TRUE(per_frame->frame_values != NULL);
        const size_t bins = per_frame->num_values;

        const md_script_eval_desc_t desc = { .num_frames = n, .frame_budget = kept ? 0 : 1 };
        md_script_eval_t* eval = md_script_eval_create_desc(ir, &desc, alloc);
        const md_bitset_t all = md_bitset_range(0, n);
        ASSERT_TRUE(md_script_eval_slot_request(eval, 1, STR_LIT("r"), &all));

        eval_job_t job = { .eval = eval, .ir = ir, .sys = mol, .num_frames = n };
        md_thread_t* thread = md_thread_create(eval_job_run, &job);
        ASSERT_TRUE(thread != NULL);

        float* values  = md_alloc(alloc, bins * sizeof(float));
        float* weights = md_alloc(alloc, bins * sizeof(float));
        double* expect = md_alloc(alloc, bins * sizeof(double));
        uint32_t reads = 0, mismatches = 0, last = 0;
        bool backwards = false;
        for (;;) {
            const bool finished = md_atomic_load64(&job.done) != 0;
            md_script_eval_slot_info_t info = {0};
            if (!md_script_eval_slot_extract(&info, values, weights, bins, eval, 1, STR_LIT("r"))) break;
            reads += 1;
            backwards |= info.num_frames < last;
            last = info.num_frames;
            if (info.num_frames > 0) {
                // The first num_frames frames, values over weights summed bin by bin
                for (size_t i = 0; i < bins; ++i) expect[i] = 0;
                for (uint32_t f = 0; f < info.num_frames; ++f) {
                    for (size_t i = 0; i < bins; ++i) expect[i] += per_frame->frame_values[f * bins + i];
                }
                for (size_t i = 0; i < bins; ++i) {
                    const double e = expect[i] / info.num_frames;
                    if (fabs(values[i] - e) > 1e-4 * MAX(1.0, fabs(e))) { mismatches += 1; break; }
                }
            }
            if (finished) break;
        }
        md_thread_join(thread);
        EXPECT_GT(reads, 0u);
        EXPECT_FALSE(backwards);
        EXPECT_EQ(0u, mismatches);
        EXPECT_EQ(n, slot_info_of(eval, 1, STR_LIT("r")).num_frames);

        md_script_eval_free(eval);
        md_script_eval_free(ref);
    }
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// What an edit changes is known before anything is evaluated, and it is what md_script_eval_create_from then does
UTEST_F(script, plan_is_what_create_from_does) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);
    const md_script_eval_desc_t desc = { .num_frames = n, .frame_budget = 1 };

    md_script_ir_t* ir1 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir1, STR_LIT(
        "d = distance(1, 2);"
        "r = rdf(element('C'), element('O'), 6.0);"
        "s = rdf(element('N'), element('O'), 6.0);"), mol, NULL));
    md_script_eval_t* eval1 = md_script_eval_create_desc(ir1, &desc, alloc);
    ASSERT_TRUE(eval1);

    // Against no evaluation, everything is new
    md_script_eval_plan_t plans[8] = {0};
    md_bitset_t pending = {0};
    ASSERT_TRUE(md_script_eval_plan(plans, ARRAY_SIZE(plans), &pending, NULL, ir1, n, alloc));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_NEW, plans[0].status);
    EXPECT_EQ(n, plans[0].num_pending);
    EXPECT_EQ((size_t)n, md_bitset_count(pending));

    // Half of it evaluated, and s over frames of its own
    const md_bitset_t sub = md_bitset_range(2, n - 2);
    ASSERT_TRUE(md_script_eval_slot_request(eval1, 1, STR_LIT("s"), &sub));
    ASSERT_TRUE(md_script_eval_frame_range(eval1, ir1, mol, SCRIPT_RUN, 0, n / 2));
    EXPECT_EQ(md_script_ir_property_fingerprint(ir1, STR_LIT("r")), md_script_eval_property_fingerprint(eval1, STR_LIT("r")));
    EXPECT_EQ((uint64_t)0, md_script_eval_property_fingerprint(eval1, STR_LIT("none")));

    // Against itself, what remains is what remains
    ASSERT_TRUE(md_script_eval_plan(plans, ARRAY_SIZE(plans), &pending, eval1, ir1, n, alloc));
    EXPECT_TRUE(bits_equal(pending, md_script_eval_pending_frames(eval1, alloc)));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_HELD, plans[0].status);
    EXPECT_EQ(n / 2, plans[0].num_frames);
    EXPECT_EQ(n - n / 2, plans[0].num_pending);

    // d renamed, r changed, s as it was, a new
    md_script_ir_t* ir2 = md_script_ir_create(alloc);
    const str_t src2 = STR_LIT(
        "dist = distance(1, 2);\n"
        "r = rdf(element('C'), element('O'), 5.0);\n"
        "s = rdf(element('N'), element('O'), 6.0);\n"
        "a = angle(1, 2, 3);\n"
        "{cx, cy, cz} = coord(1);");
    ASSERT_TRUE(md_script_ir_compile_from_source(ir2, src2, mol, NULL));
    {
        // Where each is named
        const char* named[3] = { "dist", "s", "cy" };
        for (int i = 0; i < 3; ++i) {
            const md_script_range_marker_t range = md_script_ir_property_range(ir2, str_from_cstr(named[i]));
            ASSERT_LT(range.beg, range.end);
            EXPECT_TRUE(str_eq(str_substr(src2, range.beg, range.end - range.beg), str_from_cstr(named[i])));
        }
        const md_script_range_marker_t none = md_script_ir_property_range(ir2, STR_LIT("none"));
        EXPECT_EQ(none.beg, none.end);
    }
    ASSERT_TRUE(md_script_eval_plan(plans, ARRAY_SIZE(plans), &pending, eval1, ir2, n, alloc));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_HELD, plans[0].status);
    EXPECT_TRUE(str_eq(plans[0].from, STR_LIT("d")));
    EXPECT_EQ(n / 2, plans[0].num_frames);
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_CHANGED, plans[1].status);
    EXPECT_EQ(0u, plans[1].num_frames);
    EXPECT_EQ(n, plans[1].num_pending);
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_HELD, plans[2].status);
    EXPECT_TRUE(str_eq(plans[2].from, STR_LIT("s")));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_NEW, plans[3].status);

    // and so it goes
    md_script_eval_t* eval2 = md_script_eval_create_from(eval1, ir2, &desc, alloc);
    ASSERT_TRUE(eval2);
    EXPECT_TRUE(bits_equal(pending, md_script_eval_pending_frames(eval2, alloc)));
    const char* names[4] = { "dist", "r", "s", "a" };
    for (int i = 0; i < 4; ++i) {
        const str_t name = str_from_cstr(names[i]);
        EXPECT_EQ((uint64_t)plans[i].num_frames, md_bitset_count(md_script_eval_property_frames(eval2, name, alloc)));
    }

    // With s alone: what it has not evaluated, and what its slot has to add, which is within that
    md_script_ir_t* ir3 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir3, STR_LIT("s = rdf(element('N'), element('O'), 6.0);"), mol, NULL));
    ASSERT_TRUE(md_script_eval_plan(plans, ARRAY_SIZE(plans), &pending, eval1, ir3, n, alloc));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_HELD, plans[0].status);
    md_script_eval_t* eval3 = md_script_eval_create_from(eval1, ir3, &desc, alloc);
    EXPECT_TRUE(bits_equal(pending, md_script_eval_pending_frames(eval3, alloc)));
    EXPECT_EQ((size_t)(n - n / 2), md_bitset_count(pending));
    EXPECT_EQ(n - n / 2, plans[0].num_pending);

    // Too little room: nothing is written
    EXPECT_FALSE(md_script_eval_plan(plans, 3, &pending, eval1, ir2, n, alloc));

    // Over another number of frames, nothing is taken over
    ASSERT_TRUE(md_script_eval_plan(plans, ARRAY_SIZE(plans), &pending, eval1, ir3, n - 1, alloc));
    EXPECT_EQ(MD_SCRIPT_EVAL_STATUS_CHANGED, plans[0].status);

    md_script_eval_free(eval3);
    md_script_eval_free(eval2);
    md_script_eval_free(eval1);
    md_script_ir_free(ir3);
    md_script_ir_free(ir2);
    md_script_ir_free(ir1);
    md_arena_allocator_destroy(alloc);
}

// A distribution taken over is kept as it was, per frame or not, whatever the budget of the new evaluation: what it
// holds comes along, its slots too
UTEST_F(script, taken_over_distribution_is_kept_as_it_was) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 6u);
    const md_bitset_t sub = md_bitset_range(2, n - 2);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("r = rdf(element('C'), element('O'), 6.0);"), mol, NULL));
    for (int kept = 0; kept < 2; ++kept) {
        const md_script_eval_desc_t before = { .num_frames = n, .frame_budget = kept ? 0 : 1 };
        const md_script_eval_desc_t after  = { .num_frames = n, .frame_budget = kept ? 1 : 0 };
        md_script_eval_t* eval1 = md_script_eval_create_desc(ir, &before, alloc);
        ASSERT_TRUE(md_script_eval_slot_request(eval1, 1, STR_LIT("r"), &sub));
        ASSERT_TRUE(md_script_eval_frame_range(eval1, ir, mol, SCRIPT_RUN, 0, n));

        md_script_eval_t* eval2 = md_script_eval_create_from(eval1, ir, &after, alloc);
        EXPECT_EQ(kept != 0, find_eval_property(eval2, STR_LIT("r"))->frame_values != NULL);
        EXPECT_EQ((size_t)0, md_bitset_count(md_script_eval_pending_frames(eval2, alloc)));
        EXPECT_EQ(n - 4, slot_info_of(eval2, 1, STR_LIT("r")).num_frames);

        // and a new one is kept as the budget says
        md_script_eval_t* eval3 = md_script_eval_create_desc(ir, &after, alloc);
        EXPECT_EQ(kept == 0, find_eval_property(eval3, STR_LIT("r"))->frame_values != NULL);

        md_script_eval_free(eval3);
        md_script_eval_free(eval2);
        md_script_eval_free(eval1);
    }
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// ### DISTRIBUTIONS ###
// A distribution is its own type, [n][bins], with the weight of every bin beside its value rather than in the type

static bool floats_near(const float* a, const float* b, size_t n, float rel) {
    for (size_t i = 0; i < n; ++i) {
        if (fabsf(a[i] - b[i]) > rel * MAX(1.0f, fabsf(b[i]))) return false;
    }
    return true;
}

UTEST_F(script, distribution_type) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    #define RDF_CO "rdf(element('C'), element('O'), 6.0)"
    // Scaled, and pooled, it is a distribution
    EXPECT_TRUE(compiles(ir, "g = " RDF_CO ";", mol));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION, md_script_ir_property_flags(ir, STR_LIT("g")));
    EXPECT_TRUE(compiles(ir, "g = 2.0 * " RDF_CO ";", mol));
    EXPECT_TRUE(compiles(ir, "g = " RDF_CO " * 2.0;", mol));
    EXPECT_TRUE(compiles(ir, "g = " RDF_CO " / 2.0;", mol));
    EXPECT_TRUE(compiles(ir, "g = pool(" RDF_CO ", " RDF_CO ");", mol));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION, md_script_ir_property_flags(ir, STR_LIT("g")));

    // Nothing else: per frame it would not be what it is over the frames
    EXPECT_FALSE(compiles(ir, "g = " RDF_CO " + " RDF_CO ";", mol));
    EXPECT_FALSE(compiles(ir, "g = " RDF_CO " - " RDF_CO ";", mol));
    EXPECT_FALSE(compiles(ir, "g = " RDF_CO " * " RDF_CO ";", mol));
    EXPECT_FALSE(compiles(ir, "g = " RDF_CO " / " RDF_CO ";", mol));
    EXPECT_FALSE(compiles(ir, "g = " RDF_CO " + 1.0;", mol));
    EXPECT_FALSE(compiles(ir, "g = -" RDF_CO ";", mol));
    EXPECT_FALSE(compiles(ir, "g = 2.0 / " RDF_CO ";", mol));

    // An array of them is one too
    EXPECT_TRUE(compiles(ir, "g = " RDF_CO " in residue(1:3);", mol));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION, md_script_ir_property_flags(ir, STR_LIT("g")));
    EXPECT_TRUE(compiles(ir, "g = pool(" RDF_CO " in residue(1:3), " RDF_CO " in residue(4:6));", mol));
    EXPECT_FALSE(compiles(ir, "g = pool(" RDF_CO " in residue(1:3), " RDF_CO " in residue(4:5));", mol));
    #undef RDF_CO

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// Scaling scales the values and leaves the weights; pooling adds both; an array is per element, published [n][bins]
UTEST_F(script, distribution_ops_and_arrays) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(32));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t n = (uint32_t)script_frames(mol);
    ASSERT_GT(n, 0u);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "a  = rdf(element('C'), element('O'), 6.0);"
        "b  = rdf(element('N'), element('O'), 6.0);"
        "s  = 2.0 * a;"
        "h  = a / 4.0;"
        "p  = pool(a, b);"
        "ga = rdf(element('C'), element('O'), 6.0) in residue(1:3);"
        "g2 = ga[2];"
        "r2 = rdf(element('C'), element('O'), 6.0) in residue(2);"
    ), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(n, ir, alloc);
    ASSERT_TRUE(eval);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, n));

    const eval_property_t* a  = find_eval_property(eval, STR_LIT("a"));
    const eval_property_t* bb = find_eval_property(eval, STR_LIT("b"));
    const eval_property_t* s  = find_eval_property(eval, STR_LIT("s"));
    const eval_property_t* h  = find_eval_property(eval, STR_LIT("h"));
    const eval_property_t* p  = find_eval_property(eval, STR_LIT("p"));
    const eval_property_t* ga = find_eval_property(eval, STR_LIT("ga"));
    const eval_property_t* g2 = find_eval_property(eval, STR_LIT("g2"));
    const eval_property_t* r2 = find_eval_property(eval, STR_LIT("r2"));
    ASSERT_TRUE(a && bb && s && h && p && ga && g2 && r2);
    const size_t bins = a->count;
    ASSERT_EQ((size_t)a->num_values, bins);

    float* tmp = md_alloc(alloc, bins * sizeof(float));

    // Scaled: values times 2, weights as they were
    for (size_t i = 0; i < bins; ++i) tmp[i] = 2.0f * a->values[i];
    EXPECT_TRUE(floats_near(s->values,  tmp, bins, 1e-5f));
    EXPECT_TRUE(floats_near(s->weights, a->weights, bins, 1e-5f));
    for (size_t i = 0; i < bins; ++i) tmp[i] = a->values[i] / 4.0f;
    EXPECT_TRUE(floats_near(h->values,  tmp, bins, 1e-5f));
    EXPECT_TRUE(floats_near(h->weights, a->weights, bins, 1e-5f));

    // Pooled: both added
    for (size_t i = 0; i < bins; ++i) tmp[i] = a->values[i] + bb->values[i];
    EXPECT_TRUE(floats_near(p->values, tmp, bins, 1e-4f));
    for (size_t i = 0; i < bins; ++i) tmp[i] = a->weights[i] + bb->weights[i];
    EXPECT_TRUE(floats_near(p->weights, tmp, bins, 1e-4f));

    // An array: one distribution per element, each with its own weights, published as [3][bins]
    EXPECT_EQ((size_t)3 * bins, (size_t)ga->num_values);
    const md_attributes_t* table = md_script_eval_attributes(eval);
    const md_attribute_t* attr = md_attributes_find(table, STR_LIT("script/ga"));
    const md_attribute_t* wttr = md_attributes_find(table, STR_LIT("script/ga/weight"));
    ASSERT_TRUE(attr && wttr);
    EXPECT_EQ(2u, (uint32_t)attr->format.rank);
    EXPECT_EQ(3u, attr->format.shape[0]);
    EXPECT_EQ((uint32_t)bins, attr->format.shape[1]);
    EXPECT_EQ(2u, (uint32_t)wttr->format.rank);
    const md_attribute_t* single = md_attributes_find(table, STR_LIT("script/a"));
    ASSERT_TRUE(single);
    EXPECT_EQ(1u, (uint32_t)single->format.rank);

    // Its second element is the distribution over residue 2 alone, whether taken by subscript or evaluated so
    EXPECT_TRUE(floats_near(ga->values  + bins, r2->values,  bins, 1e-5f));
    EXPECT_TRUE(floats_near(ga->weights + bins, r2->weights, bins, 1e-5f));
    EXPECT_TRUE(floats_near(g2->values,  r2->values,  bins, 1e-5f));
    EXPECT_TRUE(floats_near(g2->weights, r2->weights, bins, 1e-5f));
    // and not the same as the first
    EXPECT_FALSE(floats_near(ga->weights, ga->weights + bins, bins, 1e-5f) && floats_near(ga->values, ga->values + bins, bins, 1e-5f));

    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// ---------------------------------------------------------------------------------------------------------------
// Review fixes

// An element past the end of what a reference refers to resolves to nothing, and an empty selection has no closest
// pair to draw: neither may reach the procedures, which index unchecked
UTEST_F(script, vis_ref_out_of_range_and_empty) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "d = distance_pair(residue(1), residue(2));"
        "e = distance_min(residue(1), within_x(1000.0:2000.0));"), mol, NULL));

    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    const md_script_vis_ctx_t ctx = { .ir = ir, .sys = mol, .state = &mol->reference };
    const md_script_vis_ref_t d = md_script_ir_property_vis_ref(ir, STR_LIT("d"));
    const int dim = md_script_vis_ref_dim(ir, d);
    ASSERT_GT(dim, 0);
    EXPECT_TRUE (md_script_vis_eval_ref(&vis, d, dim - 1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_FALSE(md_script_vis_eval_ref(&vis, d, dim,     &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));

    const md_script_vis_ref_t e = md_script_ir_property_vis_ref(ir, STR_LIT("e"));
    md_script_vis_clear(&vis);
    EXPECT_TRUE(md_script_vis_eval_ref(&vis, e, -1, &ctx, MD_SCRIPT_VISUALIZE_DEFAULT));

    md_script_vis_free(&vis);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

static double dist_shown_sum(const eval_property_t* p, size_t off) {
    double sum = 0;
    for (size_t i = 0; i < MD_DIST_BINS; ++i) {
        const float w = p->weights[off + i];
        sum += w != 0.0f ? p->values[off + i] / w : 0.0;
    }
    return sum;
}

// density() is the three axis profiles: each is the one of density_x/y/z over the longest extent, so the same mass
// per slab volume, whose cross section is the cell's
UTEST_F(script, density_axes_agree) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(8));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT(
        "all3 = density(all()); dx = density_x(all()); dy = density_y(all()); dz = density_z(all());"), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(1, ir, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 1));

    const eval_property_t* all3 = find_eval_property(eval, STR_LIT("all3"));
    const char* names[3] = {"dx", "dy", "dz"};
    // Over its own extent an axis has bins of its own thickness: the mass per bin times the bin thickness, summed,
    // is the total mass over the cross section whichever the extent. Compare that.
    const double L[3] = { mol->reference.unitcell.x, mol->reference.unitcell.y, mol->reference.unitcell.z };
    const double Lmax = MAX(L[0], MAX(L[1], L[2]));
    for (int a = 0; a < 3; ++a) {
        const eval_property_t* single = find_eval_property(eval, str_from_cstr(names[a]));
        ASSERT_TRUE(single != NULL);
        const double s_all = dist_shown_sum(all3, (size_t)a * MD_DIST_BINS) * (Lmax / MD_DIST_BINS);
        const double s_one = dist_shown_sum(single, 0) * (L[a] / MD_DIST_BINS);
        EXPECT_NEAR(s_one, s_all, 1e-3 * s_one);
    }
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// An rdf's weights are the pairs a uniform distribution would put in each shell: over the window they add up to
// the pairs counted, and the first shell is a shell, also when the window starts away from 0
UTEST_F(script, rdf_shell_weights) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(8));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("g = rdf(element('C'), element('O'), 3.0:8.0);"), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(1, ir, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, 1));
    const eval_property_t* g = find_eval_property(eval, STR_LIT("g"));
    ASSERT_TRUE(g != NULL);
    double sum_v = 0, sum_w = 0;
    for (size_t i = 0; i < MD_DIST_BINS; ++i) { sum_v += g->values[i]; sum_w += g->weights[i]; }
    ASSERT_GT(sum_v, 0.0);
    EXPECT_NEAR(sum_v, sum_w, 1e-3 * sum_v);
    EXPECT_LT(g->weights[0], 1.01f * g->weights[1]);
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}


// ### PLANE ###
// plane(): one least squares plane per group or one through positions, oriented, and what centroid(), normal(),
// radius(), deviation() and angle() read of it. Checked against references computed here, on dictionary components
// (no test data needed).

typedef struct plane_sys_t {
    md_allocator_i*   alloc;
    md_system_t       sys;
    md_system_state_t st;
} plane_sys_t;

static bool plane_sys_init(plane_sys_t* p, const char* id) {
    p->alloc = md_vm_arena_create(GIGABYTES(1));
    md_array(ref_comp_t) residues = read_components(STR_LIT(MD_UNITTEST_SOURCE_DIR "/ccd_chem.txt"), p->alloc);
    const ref_comp_t* c = find_component(residues, id);
    if (!c) return false;
    p->sys = (md_system_t){ .alloc = p->alloc };
    p->st  = (md_system_state_t){ .alloc = p->alloc };
    int* map = md_alloc(p->alloc, sizeof(int) * (md_array_size(c->atom) + 1));
    build_component(&p->sys, &p->st, map, c, false, false);
    md_chem_perceive(&p->sys, &p->st, MD_CHEM_FLAG_NONE);
    // The structure the topology is of, as a file gives it
    p->sys.reference.alloc = p->alloc;
    md_system_state_copy(&p->sys.reference, &p->st);
    return true;
}

static bool plane_eval(data_t* out, md_script_ir_t* ir, const char* name, plane_sys_t* p, const md_system_state_t* state) {
    identifier_t* ident = get_identifier(ir, str_from_cstr(name));
    if (!ident || !ident->node) return false;
    eval_context_t ctx = { .ir = ir, .sys = &p->sys, .temp_alloc = p->alloc, .alloc = p->alloc, .cur_state = state, .ref_state = &p->sys.reference };
    return evaluate_node_alloc(out, ident->node, &ctx, p->alloc);
}

static bool plane_compile(md_script_ir_t* ir, const char* src, plane_sys_t* p) {
    if (md_script_ir_compile_from_source(ir, str_from_cstr(src), &p->sys, NULL)) return true;
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        const md_log_token_t* e = md_script_ir_errors(ir) + i;
        printf("  '%s': %.*s\n", src, (int)e->text.len, e->text.ptr);
    }
    return false;
}

static bool plane_compile_fails_with(plane_sys_t* p, const char* src, const char* expected) {
    md_script_ir_t* ir = md_script_ir_create(p->alloc);
    const bool ok = md_script_ir_compile_from_source(ir, str_from_cstr(src), &p->sys, NULL);
    bool found = false;
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) {
        const md_log_token_t* e = md_script_ir_errors(ir) + i;
        found |= str_find_str(NULL, e->text, str_from_cstr(expected));
    }
    if (ok || !found) printf("  '%s' was expected to fail with '%s'\n", src, expected);
    md_script_ir_free(ir);
    return !ok && found;
}

UTEST(script, plane) {
    plane_sys_t trp = {0};
    ASSERT_TRUE(plane_sys_init(&trp, "TRP"));
    md_system_t* sys = &trp.sys;
    ASSERT_EQ(2, (int)md_ring_data_count(&sys->ring));

    md_script_ir_t* ir = md_script_ir_create(trp.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = plane(ring()); c = centroid(p); n = normal(p); r = radius(p); dv = deviation(p); a = angle(p[1], p[2]); s = plane(residue(1));", &trp));

    data_t p = {0}, c = {0}, n = {0}, r = {0}, dv = {0}, a = {0};
    ASSERT_TRUE(plane_eval(&p,  ir, "p",  &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&c,  ir, "c",  &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&n,  ir, "n",  &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&r,  ir, "r",  &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&dv, ir, "dv", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&a,  ir, "a",  &trp, &trp.st));
    ASSERT_EQ(TYPE_PLANE, p.type.base_type);
    ASSERT_EQ(2, (int)element_count(p));
    ASSERT_EQ(2, (int)element_count(c));
    ASSERT_EQ(2, (int)element_count(r));
    ASSERT_EQ(1, (int)element_count(a));

    for (int i = 0; i < 2; ++i) {
        // The plane of ring i: the references from its atoms, in the order of the ring
        const plane_t* pl = &((const plane_t*)p.ptr)[i];
        const int32_t* ring = md_ring_atoms(&sys->ring, i);
        const int nr = (int)md_ring_size(&sys->ring, i);
        EXPECT_EQ(nr, pl->count);
        vec3_t mean = {0};
        for (int k = 0; k < nr; ++k) mean = vec3_add(mean, trp.st.xyz[ring[k]]);
        mean = vec3_mul1(mean, 1.0f / nr);
        const vec3_t ci = ((const vec3_t*)c.ptr)[i];
        const vec3_t ni = ((const vec3_t*)n.ptr)[i];
        EXPECT_NEAR(mean.x, ci.x, 1e-4f);
        EXPECT_NEAR(mean.y, ci.y, 1e-4f);
        EXPECT_NEAR(mean.z, ci.z, 1e-4f);
        EXPECT_NEAR(1.0f, vec3_length(ni), 1e-4f);

        float rad = 0, dev = 0;
        vec3_t newell = {0};
        for (int k = 0; k < nr; ++k) {
            const vec3_t d = vec3_sub(trp.st.xyz[ring[k]], mean);
            rad += vec3_length(d);
            dev += vec3_dot(d, ni) * vec3_dot(d, ni);
            newell = vec3_add(newell, vec3_cross(d, vec3_sub(trp.st.xyz[ring[(k + 1) % nr]], mean)));
        }
        EXPECT_NEAR(rad / nr, ((const float*)r.ptr)[i], 1e-4f);
        EXPECT_NEAR(sqrtf(dev / nr), ((const float*)dv.ptr)[i], 1e-4f);
        // Aromatic rings are flat: the normal is that of the polygon, on the side the ring goes around anticlockwise
        EXPECT_LT(((const float*)dv.ptr)[i], 0.02f);
        EXPECT_GT(vec3_dot(ni, vec3_normalize(newell)), 0.999f);
    }
    // The two rings of indole lie in one plane
    EXPECT_LT(RAD_TO_DEG(((const float*)a.ptr)[0]), 3.0f);

    // Turned over (180 degrees about an axis in the plane of the first ring): the planes turn with the atoms, and the
    // normal of a group which is not a ring keeps to its atoms as well
    data_t s0 = {0}, s1 = {0}, p1 = {0};
    ASSERT_TRUE(plane_eval(&s0, ir, "s", &trp, &trp.st));
    const plane_t ring0 = ((const plane_t*)p.ptr)[0];
    const plane_t res0  = ((const plane_t*)s0.ptr)[0];
    const vec3_t axis = ring0.axis;
    md_system_state_t rot = { .alloc = trp.alloc };
    md_system_state_init(&rot, trp.st.num_atoms);
    rot.num_atoms = trp.st.num_atoms;
    for (size_t k = 0; k < trp.st.num_atoms; ++k) {
        const vec3_t x = trp.st.xyz[k];
        // 180 degrees about axis through the origin: 2 (x.a) a - x
        rot.xyz[k] = vec3_sub(vec3_mul1(axis, 2.0f * vec3_dot(x, axis)), x);
    }
    ASSERT_TRUE(plane_eval(&p1, ir, "p", &trp, &rot));
    ASSERT_TRUE(plane_eval(&s1, ir, "s", &trp, &rot));
    const vec3_t want_ring = vec3_sub(vec3_mul1(axis, 2.0f * vec3_dot(ring0.normal, axis)), ring0.normal);
    const vec3_t want_res  = vec3_sub(vec3_mul1(axis, 2.0f * vec3_dot(res0.normal, axis)), res0.normal);
    EXPECT_GT(vec3_dot(((const plane_t*)p1.ptr)[0].normal, want_ring), 0.999f);
    EXPECT_GT(vec3_dot(((const plane_t*)s1.ptr)[0].normal, want_res), 0.999f);

    // Errors before evaluation
    EXPECT_TRUE(plane_compile_fails_with(&trp, "x = plane(element('N'));", "a plane needs at least 3"));
    EXPECT_TRUE(plane_compile_fails_with(&trp, "x = plane(within(3, element('N')));", "cannot depend on the frame"));
    EXPECT_TRUE(plane_compile_fails_with(&trp, "x = angle(plane(ring()), plane({residue(1), residue(1), residue(1)}));", "as many or one of them a single plane"));
    // Positions: the centres of the two rings, and a single one
    EXPECT_TRUE(plane_compile_fails_with(&trp, "x = plane(coord(ring()));", "2 positions, a plane needs at least 3"));
    EXPECT_TRUE(plane_compile_fails_with(&trp, "x = plane(com(ring()));", "1 position, a plane needs at least 3"));

    md_script_ir_free(ir);
    md_vm_arena_destroy(trp.alloc);
}

// The number of warnings of ir which contain text
static int warnings_with(const md_script_ir_t* ir, const char* text) {
    int n = 0;
    for (size_t i = 0; i < md_script_ir_num_warnings(ir); ++i) {
        n += str_find_str(NULL, md_script_ir_warnings(ir)[i].text, str_from_cstr(text)) ? 1 : 0;
    }
    return n;
}

// plane(union(groups)): one plane through all of the atoms of the groups. plane(positions): one through the positions,
// atom indices or float[3], its normal on the side from which they go around anticlockwise in the order given.
UTEST(script, plane_union_positions) {
    plane_sys_t trp = {0};
    ASSERT_TRUE(plane_sys_init(&trp, "TRP"));
    md_system_t* sys = &trp.sys;

    md_script_ir_t* ir = md_script_ir_create(trp.alloc);
    ASSERT_TRUE(plane_compile(ir,
        "p = plane(ring());\n"
        "u = plane(union(ring()));\n"
        "t = plane({1, 2, 3});\n"
        "b = plane({3, 2, 1});\n"
        "f = plane(coord(1:3));\n"
        "i = plane(1:4);\n"
        "k = plane(coord(within(3.0, atom(1))));\n", &trp));
    EXPECT_EQ(0, (int)md_script_ir_num_warnings(ir));

    data_t p = {0}, u = {0}, t = {0}, b = {0}, f = {0}, i = {0};
    ASSERT_TRUE(plane_eval(&p, ir, "p", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&u, ir, "u", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&t, ir, "t", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&b, ir, "b", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&f, ir, "f", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&i, ir, "i", &trp, &trp.st));
    ASSERT_EQ(2, (int)element_count(p));
    ASSERT_EQ(1, (int)element_count(u));
    ASSERT_EQ(1, (int)element_count(t));
    ASSERT_EQ(1, (int)element_count(i));

    // Both rings of indole as one: their 9 atoms (two shared), in the plane of either
    const plane_t* pu = (const plane_t*)u.ptr;
    const plane_t* p0 = (const plane_t*)p.ptr;
    EXPECT_EQ(9, pu->count);
    EXPECT_GT(fabsf(vec3_dot(pu->normal, p0->normal)), 0.999f);
    EXPECT_LT(pu->deviation, 0.05f);

    // Three atoms in order, and the other way around: the normal is on the side they go around anticlockwise. The
    // same positions as float[3] give the same plane.
    const plane_t* pt = (const plane_t*)t.ptr;
    const plane_t* pb = (const plane_t*)b.ptr;
    const plane_t* pf = (const plane_t*)f.ptr;
    EXPECT_EQ(3, pt->count);
    const vec3_t x0 = trp.st.xyz[0], x1 = trp.st.xyz[1], x2 = trp.st.xyz[2];
    const vec3_t turn = vec3_normalize(vec3_cross(vec3_sub(x1, x0), vec3_sub(x2, x0)));
    EXPECT_GT(vec3_dot(pt->normal, turn), 0.999f);
    EXPECT_GT(vec3_dot(pb->normal, turn), -1.001f);
    EXPECT_LT(vec3_dot(pb->normal, turn), -0.999f);
    EXPECT_GT(vec3_dot(pf->normal, pt->normal), 0.999f);
    const vec3_t c3 = vec3_mul1(vec3_add(vec3_add(x0, x1), x2), 1.0f / 3.0f);
    EXPECT_NEAR(c3.x, pt->centroid.x, 1e-4f);
    EXPECT_NEAR(c3.y, pt->centroid.y, 1e-4f);
    EXPECT_NEAR(c3.z, pt->centroid.z, 1e-4f);
    EXPECT_NEAR(0.0f, pt->deviation, 1e-4f);

    // Positions which depend on the frame: the atoms around atom 1
    data_t k = {0};
    ASSERT_TRUE(plane_eval(&k, ir, "k", &trp, &trp.st));
    ASSERT_EQ(1, (int)element_count(k));
    EXPECT_GE(((const plane_t*)k.ptr)->count, 3);

    // Atom indices: the atoms
    const plane_t* pi = (const plane_t*)i.ptr;
    EXPECT_EQ(4, pi->count);
    vec3_t mean = {0};
    for (int k = 0; k < 4; ++k) mean = vec3_add(mean, trp.st.xyz[k]);
    mean = vec3_mul1(mean, 0.25f);
    EXPECT_NEAR(mean.x, pi->centroid.x, 1e-4f);
    EXPECT_NEAR(mean.y, pi->centroid.y, 1e-4f);
    EXPECT_NEAR(mean.z, pi->centroid.z, 1e-4f);
    (void)sys;
    md_script_ir_free(ir);

    // The old names of union compile, with one warning each
    ir = md_script_ir_create(trp.alloc);
    ASSERT_TRUE(plane_compile(ir, "x = plane(join(ring()));\ny = plane(flatten(ring()));\nz = count(flatten(ring()));\n", &trp));
    EXPECT_EQ(1, warnings_with(ir, "join() is deprecated, use union()"));
    EXPECT_EQ(2, warnings_with(ir, "use union()") - warnings_with(ir, "join() is deprecated"));
    data_t x = {0}, y = {0};
    ASSERT_TRUE(plane_eval(&x, ir, "x", &trp, &trp.st));
    ASSERT_TRUE(plane_eval(&y, ir, "y", &trp, &trp.st));
    EXPECT_EQ(9, ((const plane_t*)x.ptr)->count);
    EXPECT_EQ(9, ((const plane_t*)y.ptr)->count);
    md_script_ir_free(ir);

    md_vm_arena_destroy(trp.alloc);
}

// A sugar ring is in order from its ring O to the anomeric C, and its normal follows that order
UTEST(script, plane_sugar) {
    plane_sys_t ade = {0};
    ASSERT_TRUE(plane_sys_init(&ade, "A"));
    md_script_ir_t* ir = md_script_ir_create(ade.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = plane(ring(kind='sugar'));", &ade));
    data_t p = {0};
    ASSERT_TRUE(plane_eval(&p, ir, "p", &ade, &ade.st));
    ASSERT_EQ(1, (int)element_count(p));
    const plane_t* pl = (const plane_t*)p.ptr;
    // The ring: O4', C1', C2', C3', C4'
    int64_t r = -1;
    for (size_t i = 0; i < md_ring_data_count(&ade.sys.ring); ++i) {
        if (md_system_ring_flags(&ade.sys, i) & MD_RING_FLAG_SUGAR) r = (int64_t)i;
    }
    ASSERT_GE(r, 0);
    const int32_t* ring = md_ring_atoms(&ade.sys.ring, (size_t)r);
    EXPECT_TRUE(str_eq(md_atom_name(&ade.sys.atom, ring[0]), STR_LIT("O4'")));
    vec3_t newell = {0};
    for (int k = 0; k < 5; ++k) {
        newell = vec3_add(newell, vec3_cross(vec3_sub(ade.st.xyz[ring[k]], pl->centroid), vec3_sub(ade.st.xyz[ring[(k + 1) % 5]], pl->centroid)));
    }
    EXPECT_GT(vec3_dot(pl->normal, vec3_normalize(newell)), 0.99f);
    // A furanose is puckered: some atoms out of the plane
    EXPECT_GT(pl->deviation, 0.05f);
    md_script_ir_free(ir);
    md_vm_arena_destroy(ade.alloc);
}

// ### HOVER ###
// What hovering shows: the type, and what the value is. Calls with defaults or absent arguments show their own type
// (not that of an argument the call did not write), selections their atoms, and md_script_vis_value_text the value at
// the given state, for contacts and planes as well.

static const md_script_vis_token_t* hover_token(const md_script_ir_t* ir, const char* src, const char* text) {
    // The outermost token which spans exactly text
    const size_t len = strlen(text);
    const char* at = strstr(src, text);
    if (!at) return NULL;
    const int beg = (int)(at - src);
    const md_script_vis_token_t* tok = md_script_ir_vis_tokens(ir);
    const md_script_vis_token_t* best = NULL;
    for (size_t i = 0; i < md_script_ir_num_vis_tokens(ir); ++i) {
        if (tok[i].range.beg == beg && tok[i].range.end == beg + (int)len && (!best || tok[i].depth < best->depth)) best = &tok[i];
    }
    return best;
}

UTEST(script, hover_text) {
    plane_sys_t trp = {0};
    ASSERT_TRUE(plane_sys_init(&trp, "TRP"));
    const char* src = "k = contacts(ring(), cutoff=5.0);\nm = smiles('cc', level='residue');\nr = ring();\np = plane(ring());\nd = distance(1, 2);";
    md_script_ir_t* ir = md_script_ir_create(trp.alloc);
    ASSERT_TRUE(plane_compile(ir, src, &trp));

    for (size_t i = 0; i < md_script_ir_num_vis_tokens(ir); ++i) {
        const str_t t = md_script_ir_vis_tokens(ir)[i].text;
        EXPECT_FALSE(str_find_str(NULL, t, STR_LIT("NULL")));
        EXPECT_FALSE(str_find_str(NULL, t, STR_LIT("undefined")));
    }

    const md_script_vis_token_t* t;
    t = hover_token(ir, src, "contacts(ring(), cutoff=5.0)");
    ASSERT_TRUE(t != NULL);
    EXPECT_TRUE(str_begins_with(t->text, STR_LIT("contact")));
    t = hover_token(ir, src, "smiles('cc', level='residue')");
    ASSERT_TRUE(t != NULL);
    EXPECT_TRUE(str_begins_with(t->text, STR_LIT("atomset")));
    t = hover_token(ir, src, "ring()");
    ASSERT_TRUE(t != NULL);
    EXPECT_TRUE(str_find_str(NULL, t->text, STR_LIT("2 selections, 9 atoms in all")));

    // At the state given
    const md_script_vis_ctx_t ctx = { .ir = ir, .sys = &trp.sys, .state = &trp.st };
    char buf[1024];
    t = hover_token(ir, src, "contacts(ring(), cutoff=5.0)");
    EXPECT_GT(md_script_vis_value_text(buf, sizeof(buf), t->ref, -1, &ctx), 0);
    EXPECT_TRUE(strstr(buf, "of groups in contact") != NULL);
    t = hover_token(ir, src, "plane(ring())");
    EXPECT_GT(md_script_vis_value_text(buf, sizeof(buf), t->ref, -1, &ctx), 0);
    EXPECT_TRUE(strstr(buf, "2 planes") != NULL);
    EXPECT_GT(md_script_vis_value_text(buf, sizeof(buf), t->ref, 1, &ctx), 0);
    EXPECT_TRUE(strstr(buf, "plane 2 of 2") != NULL && strstr(buf, "normal") != NULL);
    t = hover_token(ir, src, "distance(1, 2)");
    EXPECT_GT(md_script_vis_value_text(buf, sizeof(buf), t->ref, -1, &ctx), 0);
    const float d = vec3_distance(trp.st.xyz[0], trp.st.xyz[1]);
    char want[32];
    snprintf(want, sizeof(want), "%.2f", d);
    EXPECT_TRUE(strstr(buf, want) != NULL);

    md_script_ir_free(ir);
    md_vm_arena_destroy(trp.alloc);
}

// ### PUCKER ###
// pucker(): the Cremer-Pople coordinates of rings (md_pucker), in the order the conventions are written for, and the
// name of the conformation. The numbers themselves are checked against an independent implementation outside of the
// unit tests; here what the script makes of them.

static bool pucker_eval(data_t* out, md_script_ir_t* ir, const char* name, plane_sys_t* p, const md_system_state_t* state) {
    return plane_eval(out, ir, name, p, state);
}

UTEST(script, pucker_pyranose) {
    plane_sys_t nag = {0};
    ASSERT_TRUE(plane_sys_init(&nag, "NAG"));
    md_script_ir_t* ir = md_script_ir_create(nag.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = pucker(ring(kind='sugar')); q = amplitude(p); t = theta(p); f = phase(p);", &nag));
    data_t p = {0}, q = {0}, t = {0};
    ASSERT_TRUE(pucker_eval(&p, ir, "p", &nag, &nag.st));
    ASSERT_TRUE(pucker_eval(&q, ir, "q", &nag, &nag.st));
    ASSERT_TRUE(pucker_eval(&t, ir, "t", &nag, &nag.st));
    ASSERT_EQ(1, (int)element_count(p));
    const pucker_rec_t* r = (const pucker_rec_t*)p.ptr;
    // beta-D-GlcNAc is a 4C1 chair: theta near 0 in the order O5, C1, ..., C5
    EXPECT_STREQ("4C1 chair", r->name);
    EXPECT_LT(RAD_TO_DEG(((const float*)t.ptr)[0]), 10.0f);
    EXPECT_NEAR(0.59f, ((const float*)q.ptr)[0], 0.05f);
    EXPECT_TRUE(r->anchored);

    // The atoms given the other way around: the same ring seen from the other side, theta goes to 180 - theta
    // O5 C1 C2 C3 C4 C5 by their indices, then reversed
    const int32_t* ring = NULL;
    for (size_t i = 0; i < md_ring_data_count(&nag.sys.ring); ++i) {
        if (md_system_ring_flags(&nag.sys, i) & MD_RING_FLAG_SUGAR) ring = md_ring_atoms(&nag.sys.ring, i);
    }
    ASSERT_TRUE(ring != NULL);
    char src[256];
    snprintf(src, sizeof(src), "u = theta(pucker({%d, %d, %d, %d, %d, %d}));", ring[0] + 1, ring[5] + 1, ring[4] + 1, ring[3] + 1, ring[2] + 1, ring[1] + 1);
    md_script_ir_t* ir2 = md_script_ir_create(nag.alloc);
    ASSERT_TRUE(plane_compile(ir2, src, &nag));
    data_t u = {0};
    ASSERT_TRUE(pucker_eval(&u, ir2, "u", &nag, &nag.st));
    EXPECT_NEAR(180.0f - RAD_TO_DEG(((const float*)t.ptr)[0]), RAD_TO_DEG(((const float*)u.ptr)[0]), 0.01f);

    // Errors before evaluation
    EXPECT_TRUE(plane_compile_fails_with(&nag, "x = pseudorotation(pucker(ring(kind='sugar')));", "is not a furanose"));
    EXPECT_TRUE(plane_compile_fails_with(&nag, "x = pucker(residue(1));", "a ring of 5 to 7 is needed"));
    snprintf(src, sizeof(src), "x = pucker({%d, %d, %d, %d, %d, %d});", ring[0] + 1, ring[2] + 1, ring[1] + 1, ring[3] + 1, ring[4] + 1, ring[5] + 1);
    EXPECT_TRUE(plane_compile_fails_with(&nag, src, "are not bonded"));
    EXPECT_TRUE(plane_compile_fails_with(&nag, "x = pucker({1, 2, 3});", "a ring of 5 to 7 atoms in order"));

    md_script_ir_free(ir2);
    md_script_ir_free(ir);
    md_vm_arena_destroy(nag.alloc);
}

// A furanose: the pseudorotation of nucleic acids, and endo and exo by the side of C5'. In a mirror image (the L-sugar)
// the IUPAC name turns over and P moves by 180 degrees, but endo stays on the side of C5'.
UTEST(script, pucker_furanose) {
    plane_sys_t ade = {0};
    ASSERT_TRUE(plane_sys_init(&ade, "A"));
    md_script_ir_t* ir = md_script_ir_create(ade.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = pucker(ring(kind='sugar')); P = pseudorotation(p); f = phase(p);", &ade));
    data_t p = {0}, P = {0}, f = {0};
    ASSERT_TRUE(pucker_eval(&p, ir, "p", &ade, &ade.st));
    ASSERT_TRUE(pucker_eval(&P, ir, "P", &ade, &ade.st));
    ASSERT_TRUE(pucker_eval(&f, ir, "f", &ade, &ade.st));
    const pucker_rec_t* r = (const pucker_rec_t*)p.ptr;
    // Altona-Sundaralingam P is Cremer-Pople phi2 plus 90 degrees, to within a degree
    const float d = fmodf(RAD_TO_DEG(((const float*)P.ptr)[0]) - RAD_TO_DEG(((const float*)f.ptr)[0]) - 90.0f + 720.0f, 360.0f);
    EXPECT_LT(MIN(d, 360.0f - d), 1.5f);
    // The name agrees with P: an envelope at P = 18 + 36 k, a twist between
    const char* phrase = strchr(r->name, ',');
    ASSERT_TRUE(phrase != NULL);
    char before[64];
    snprintf(before, sizeof(before), "%s", r->name);
    const float P0 = r->pseudo_phase;
    // A D-ribose, so the convention of nucleic acids holds: P 72 to 108 is O4'-endo, 0 to 36 C3'-endo, 144 to 180 C2'-endo
    if (P0 > 72.0f && P0 < 108.0f)       EXPECT_TRUE(strstr(r->name, "O4'-endo") != NULL);
    else if (P0 >= 0.0f && P0 < 36.0f)   EXPECT_TRUE(strstr(r->name, "C3'-endo") != NULL);
    else if (P0 > 144.0f && P0 < 180.0f) EXPECT_TRUE(strstr(r->name, "C2'-endo") != NULL);
    else                                 EXPECT_TRUE(strstr(r->name, "-endo") != NULL);

    // The mirror image
    md_system_state_t mir = { .alloc = ade.alloc };
    md_system_state_init(&mir, ade.st.num_atoms);
    mir.num_atoms = ade.st.num_atoms;
    for (size_t k = 0; k < ade.st.num_atoms; ++k) mir.xyz[k] = vec3_set(-ade.st.xyz[k].x, ade.st.xyz[k].y, ade.st.xyz[k].z);
    data_t pm = {0};
    ASSERT_TRUE(pucker_eval(&pm, ir, "p", &ade, &mir));
    const pucker_rec_t* rm = (const pucker_rec_t*)pm.ptr;
    const float dP = fmodf(rm->pseudo_phase - P0 + 720.0f, 360.0f);
    EXPECT_NEAR(180.0f, dP, 1.0f);
    EXPECT_STRNE(before, rm->name);
    EXPECT_STREQ(strchr(before, ','), strchr(rm->name, ','));

    md_script_ir_free(ir);
    md_vm_arena_destroy(ade.alloc);
}

// Proline from N towards CA; rings alike all around get the kind of their conformation only, and warn of a phase
UTEST(script, pucker_order) {
    plane_sys_t pro = {0}, phe = {0};
    ASSERT_TRUE(plane_sys_init(&pro, "PRO"));
    ASSERT_TRUE(plane_sys_init(&phe, "PHE"));
    int32_t order[MD_RING_MAX_SIZE];
    bool anchored = false;
    ASSERT_EQ(5, (int)md_pucker_ring_order(order, &anchored, &pro.sys, 0));
    EXPECT_TRUE(anchored);
    EXPECT_TRUE(str_eq(md_atom_name(&pro.sys.atom, order[0]), STR_LIT("N")));
    EXPECT_TRUE(str_eq(md_atom_name(&pro.sys.atom, order[1]), STR_LIT("CA")));
    EXPECT_TRUE(str_eq(md_atom_name(&pro.sys.atom, order[2]), STR_LIT("CB")));
    EXPECT_TRUE(str_eq(md_atom_name(&pro.sys.atom, order[4]), STR_LIT("CD")));

    ASSERT_EQ(6, (int)md_pucker_ring_order(order, &anchored, &phe.sys, 0));
    EXPECT_TRUE(str_eq(md_atom_name(&phe.sys.atom, order[0]), STR_LIT("CG")));      // the C bonded to CB
    EXPECT_FALSE(anchored);                                                          // but either way round

    md_script_ir_t* ir = md_script_ir_create(phe.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = pucker(ring()); f = phase(p);", &phe));
    EXPECT_GT((int)md_script_ir_num_warnings(ir), 0);
    data_t p = {0};
    ASSERT_TRUE(pucker_eval(&p, ir, "p", &phe, &phe.st));
    EXPECT_STREQ("planar", ((const pucker_rec_t*)p.ptr)->name);

    md_script_ir_free(ir);
    md_vm_arena_destroy(pro.alloc);
    md_vm_arena_destroy(phe.alloc);
}

// The table of canonical conformations names itself back, and reads as the conventions do
UTEST(script, pucker_canonical) {
    md_pucker_canonical_t c6[64], c5[64];
    ASSERT_EQ(38, (int)md_pucker_canonical(c6, ARRAY_SIZE(c6), 6));
    ASSERT_EQ(20, (int)md_pucker_canonical(c5, ARRAY_SIZE(c5), 5));
    EXPECT_EQ(0, (int)md_pucker_canonical(NULL, 0, 7));
    EXPECT_EQ(38, (int)md_pucker_canonical(NULL, 0, 6));

    const char* loc6[] = {"O", "1", "2", "3", "4", "5"};
    const char* loc5[] = {"O", "1", "2", "3", "4"};
    char buf[32];
    for (int i = 0; i < 38; ++i) {
        md_pucker_t p = {6, 0.6f, 0, c6[i].phi2, 0, 0, c6[i].theta};
        const md_pucker_conformer_t c = md_pucker_conformer(&p);
        EXPECT_EQ((int)c6[i].conf.kind, (int)c.kind);
        EXPECT_EQ(c6[i].conf.above[0], c.above[0]);
        EXPECT_EQ(c6[i].conf.below[0], c.below[0]);
    }
    for (int i = 0; i < 20; ++i) {
        md_pucker_t p = {5, 0.4f, 0, c5[i].phi2, 0, 0, 0};
        const md_pucker_conformer_t c = md_pucker_conformer(&p);
        EXPECT_EQ((int)c5[i].conf.kind, (int)c.kind);
        EXPECT_EQ(c5[i].conf.above[0], c.above[0]);
        EXPECT_EQ(c5[i].conf.below[0], c.below[0]);
    }
    md_pucker_conformer_name(buf, sizeof(buf), c6[0].conf, loc6);
    EXPECT_STREQ("4C1", buf);
    md_pucker_locants(buf, sizeof(buf), c6[2].conf.above, loc6);     // the boat at phi 0
    EXPECT_STREQ("3,O", buf);
    md_pucker_locants(buf, sizeof(buf), c6[2].conf.below, loc6);
    EXPECT_STREQ("", buf);
    // 3E at P = phi2 + 90 = 18, C3'-endo
    for (int i = 0; i < 20; ++i) {
        if (c5[i].phi2 == 288.0f) {
            md_pucker_conformer_name(buf, sizeof(buf), c5[i].conf, loc5);
            EXPECT_STREQ("3E", buf);
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------
// Taking over from a previous evaluation (md_script_eval_create_from)

// What is computed from the same thing as before is taken over, values and frames, and not evaluated again; what is
// new, or computed from something else now, starts over. The versions published carry on from the previous ones.
UTEST_F(script, eval_create_from_keeps_what_did_not_change) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)script_frames(mol);
    ASSERT_GE(num_frames, 4u);
    const md_script_eval_desc_t desc = { .num_frames = num_frames };

    md_script_ir_t* ir1 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir1, STR_LIT("d = distance(1, 2);\ng = rdf(element('C'), element('O'), 6.0);"), mol, NULL));
    md_script_eval_t* eval1 = md_script_eval_create_desc(ir1, &desc, alloc);
    ASSERT_TRUE(eval1);
    ASSERT_TRUE(md_script_eval_frame_range(eval1, ir1, mol, SCRIPT_RUN, 0, num_frames));
    const eval_property_t* d1 = find_eval_property(eval1, STR_LIT("d"));
    const eval_property_t* g1 = find_eval_property(eval1, STR_LIT("g"));
    ASSERT_TRUE(d1 && g1);

    uint64_t max_version = 0;
    for (md_attribute_iter_t it = md_attributes_iter(md_script_eval_attributes(eval1), (str_t){0}); md_attributes_next(&it);) {
        max_version = MAX(max_version, it.attr->version);
    }

    // d and g as they were, under other names and after a line of its own; a is new
    md_script_ir_t* ir2 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir2, STR_LIT("a = angle(1, 2, 3);\ndist = distance(1, 2);\ngr = rdf(element('C'), element('O'), 6.0);"), mol, NULL));
    md_script_eval_t* eval2 = md_script_eval_create_from(eval1, ir2, &desc, alloc);
    ASSERT_TRUE(eval2);

    eval_property_t* dist = (eval_property_t*)find_eval_property(eval2, STR_LIT("dist"));
    const eval_property_t* gr = find_eval_property(eval2, STR_LIT("gr"));
    const eval_property_t* a  = find_eval_property(eval2, STR_LIT("a"));
    ASSERT_TRUE(dist && gr && a);

    EXPECT_EQ((uint64_t)num_frames, md_bitset_count(md_script_eval_property_frames(eval2, STR_LIT("dist"), alloc)));
    EXPECT_EQ((uint64_t)num_frames, md_bitset_count(md_script_eval_property_frames(eval2, STR_LIT("gr"), alloc)));
    EXPECT_EQ((uint64_t)0,          md_bitset_count(md_script_eval_property_frames(eval2, STR_LIT("a"), alloc)));
    EXPECT_TRUE(floats_near(dist->values, d1->values, num_frames, 0.0f));
    EXPECT_TRUE(floats_near(gr->values, g1->values, gr->num_values, 0.0f));
    EXPECT_EQ(g1->accum_count, gr->accum_count);

    // Every frame remains to be evaluated, for a alone
    EXPECT_EQ((uint64_t)num_frames, md_bitset_count(md_script_eval_pending_frames(eval2, alloc)));

    // Nothing published now has a version of before
    for (md_attribute_iter_t it = md_attributes_iter(md_script_eval_attributes(eval2), (str_t){0}); md_attributes_next(&it);) {
        EXPECT_GT(it.attr->version, max_version);
    }

    // Evaluating does not touch what was taken over: a value planted in it stays
    const float planted = 12345.0f;
    dist->values[0] = planted;
    ASSERT_TRUE(md_script_eval_frame_range(eval2, ir2, mol, SCRIPT_RUN, 0, num_frames));
    EXPECT_EQ(planted, dist->values[0]);
    EXPECT_EQ(g1->accum_count, gr->accum_count);
    EXPECT_EQ((uint64_t)0, md_bitset_count(md_script_eval_pending_frames(eval2, alloc)));

    // and what is new comes out as it would from scratch
    md_script_eval_t* fresh = md_script_eval_create_desc(ir2, &desc, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(fresh, ir2, mol, SCRIPT_RUN, 0, num_frames));
    const eval_property_t* a_fresh = find_eval_property(fresh, STR_LIT("a"));
    ASSERT_TRUE(a_fresh);
    EXPECT_TRUE(floats_near(a->values, a_fresh->values, num_frames, 1e-6f));

    // Computed from something else now: starts over, whatever it is called
    md_script_ir_t* ir3 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir3, STR_LIT("dist = distance(1, 3);"), mol, NULL));
    md_script_eval_t* eval3 = md_script_eval_create_from(eval2, ir3, &desc, alloc);
    ASSERT_TRUE(eval3);
    EXPECT_EQ((uint64_t)0, md_bitset_count(md_script_eval_property_frames(eval3, STR_LIT("dist"), alloc)));

    // Over another number of frames, another run: nothing is taken over, yet the versions still carry on
    for (md_attribute_iter_t it = md_attributes_iter(md_script_eval_attributes(eval2), (str_t){0}); md_attributes_next(&it);) {
        max_version = MAX(max_version, it.attr->version);
    }
    const md_script_eval_desc_t fewer = { .num_frames = num_frames - 1 };
    md_script_eval_t* eval4 = md_script_eval_create_from(eval2, ir2, &fewer, alloc);
    ASSERT_TRUE(eval4);
    EXPECT_EQ((uint64_t)0, md_bitset_count(md_script_eval_property_frames(eval4, STR_LIT("dist"), alloc)));
    EXPECT_EQ((uint64_t)(num_frames - 1), md_bitset_count(md_script_eval_pending_frames(eval4, alloc)));
    for (md_attribute_iter_t it = md_attributes_iter(md_script_eval_attributes(eval4), (str_t){0}); md_attributes_next(&it);) {
        EXPECT_GT(it.attr->version, max_version);
    }

    md_script_eval_free(eval4);
    md_script_eval_free(eval3);
    md_script_eval_free(fresh);
    md_script_eval_free(eval2);
    md_script_eval_free(eval1);
    md_script_ir_free(ir3);
    md_script_ir_free(ir2);
    md_script_ir_free(ir1);
    md_arena_allocator_destroy(alloc);
}

// A frame only some properties still need evaluates what those need: the variables they use, not the others
UTEST_F(script, eval_evaluates_only_what_is_needed) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)script_frames(mol);
    ASSERT_GE(num_frames, 4u);
    const md_script_eval_desc_t desc = { .num_frames = num_frames };

    md_script_ir_t* ir1 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir1, STR_LIT("s = distance(1, 2);\nt = s * 2.0;\nu = angle(1, 2, 3);"), mol, NULL));
    md_script_eval_t* eval1 = md_script_eval_create_desc(ir1, &desc, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(eval1, ir1, mol, SCRIPT_RUN, 0, num_frames));

    // u changes, s and t do not: t needs s, u needs neither
    md_script_ir_t* ir2 = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir2, STR_LIT("s = distance(1, 2);\nt = s * 2.0;\nu = angle(1, 2, 4) + 0.0 * s;"), mol, NULL));
    md_script_eval_t* eval2 = md_script_eval_create_from(eval1, ir2, &desc, alloc);
    ASSERT_TRUE(eval2);

    eval_property_t* s = (eval_property_t*)find_eval_property(eval2, STR_LIT("s"));
    eval_property_t* t = (eval_property_t*)find_eval_property(eval2, STR_LIT("t"));
    const eval_property_t* u = find_eval_property(eval2, STR_LIT("u"));
    ASSERT_TRUE(s && t && u);
    EXPECT_EQ((uint64_t)num_frames, md_bitset_count(md_script_eval_property_frames(eval2, STR_LIT("t"), alloc)));

    const float planted = 12345.0f;
    s->values[1] = planted;
    t->values[1] = planted;
    ASSERT_TRUE(md_script_eval_frame_range(eval2, ir2, mol, SCRIPT_RUN, 0, num_frames));
    EXPECT_EQ(planted, s->values[1]);
    EXPECT_EQ(planted, t->values[1]);

    // u used s as it was evaluated for the frame, not as it was planted
    md_script_eval_t* fresh = md_script_eval_create_desc(ir2, &desc, alloc);
    ASSERT_TRUE(md_script_eval_frame_range(fresh, ir2, mol, SCRIPT_RUN, 0, num_frames));
    EXPECT_TRUE(floats_near(u->values, find_eval_property(fresh, STR_LIT("u"))->values, num_frames, 1e-6f));

    md_script_eval_free(fresh);
    md_script_eval_free(eval2);
    md_script_eval_free(eval1);
    md_script_ir_free(ir2);
    md_script_ir_free(ir1);
    md_arena_allocator_destroy(alloc);
}

// The diagnostics of the last compilation of ir, as md_script_ir_write_diagnostics writes them
static str_t diagnostics_text(md_str_writer_t* sw, const md_script_ir_t* ir) {
    md_str_writer_reset(sw);
    md_script_ir_write_diagnostics(&sw->w, ir);
    return md_str_writer_to_str(sw);
}

// The diagnostic written for one log token in src
static const char* diagnostics_of(md_str_writer_t* sw, const char* src, int beg, int end, bool is_error, const char* msg) {
    md_str_writer_reset(sw);
    const md_log_token_t tok = { .range = { beg, end }, .text = str_from_cstr(msg) };
    diag_append(&sw->w, str_from_cstr(src), &tok, is_error);
    return md_str_writer_to_cstr(sw);
}

static bool diagnostics_contain(str_t txt, const char* expected) {
    if (str_find_str(NULL, txt, str_from_cstr(expected))) return true;
    printf("Expected the diagnostics to contain:\n%s\nGot:\n" STR_FMT, expected, STR_ARG(txt));
    return false;
}

UTEST(script, diagnostics) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_str_writer_t sw = md_str_writer_create(alloc);

    // Line and column count characters: the tab and the two byte 'Å' are one each. The echoed line expands the tab
    // and the marker beneath it follows, so the caret is under the 'y'
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("a = 1;\n\tb = \"\xC3\x85\" + y;"), &test_mol, NULL));
    str_t txt = diagnostics_text(&sw, ir);
    EXPECT_TRUE(diagnostics_contain(txt, "<script>:2:12: error: "));
    EXPECT_TRUE(diagnostics_contain(txt,
        "    2 |     b = \"\xC3\x85\" + y;\n"
        "      |               ^\n"
        "1 error generated.\n"));

    // A warning of a script which compiles, its range underlined
    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("c = contacts(all, cutoff=5.0);"), &test_mol, NULL));
    txt = diagnostics_text(&sw, ir);
    EXPECT_TRUE(diagnostics_contain(txt, "<script>:1:14: warning: "));
    EXPECT_TRUE(diagnostics_contain(txt,
        "    1 | c = contacts(all, cutoff=5.0);\n"
        "      |              ^~~\n"
        "1 warning generated.\n"));

    // Nothing to report
    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("a = 1;"), &test_mol, NULL));
    EXPECT_TRUE(str_empty(diagnostics_text(&sw, ir)));

    // A range over several lines is underlined on each, indentation left out
    EXPECT_STREQ(
        "<script>:1:5: error: E\n"
        "    1 | x = f(a,\n"
        "      |     ^~~~\n"
        "    2 |   b);\n"
        "      |   ~~~\n",
        diagnostics_of(&sw, "x = f(a,\n  b);", 4, 15, true, "E"));

    // Over more than five lines, the first and last two are shown
    EXPECT_STREQ(
        "<script>:1:5: error: E\n"
        "    1 | x = f(a,\n"
        "      |     ^~~~\n"
        "    2 |   b,\n"
        "      |   ~~\n"
        "      | ...\n"
        "    5 |   e,\n"
        "      |   ~~\n"
        "    6 |   g);\n"
        "      |   ~~~\n",
        diagnostics_of(&sw, "x = f(a,\n  b,\n  c,\n  d,\n  e,\n  g);", 4, 34, true, "E"));

    // CRLF line endings: the carriage return is not echoed
    EXPECT_STREQ(
        "<script>:2:5: warning: W\n"
        "    2 | b = y;\n"
        "      |     ^\n",
        diagnostics_of(&sw, "a = 1;\r\nb = y;\r\n", 12, 13, false, "W"));

    // A range outside the source is clamped to it: an empty range at its end
    EXPECT_STREQ(
        "<script>:1:7: error: E\n"
        "    1 | a = 1;\n"
        "      |       ^\n",
        diagnostics_of(&sw, "a = 1;", 100, 200, true, "E"));

    // The gutter widens for line numbers of more than five digits
    md_str_writer_t big = md_str_writer_create(alloc);
    md_write_fill(&big.w, '\n', 100000);
    md_str_writer_push_cstr(&big, "y;");
    md_str_writer_t sw_big = md_str_writer_create(alloc);
    EXPECT_STREQ(
        "<script>:100001:1: error: E\n"
        "100001 | y;\n"
        "       | ^\n",
        diagnostics_of(&sw_big, md_str_writer_to_cstr(&big), 100000, 100001, true, "E"));

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// An operator where an operand belongs is an error. It used to be accepted by the parser as a node with a NULL
// operand, which the static check then crashed on
UTEST(script, stray_mult) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_str_writer_t sw = md_str_writer_create(alloc);

    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("a = distance(1,*2);"), &test_mol, NULL));
    str_t txt = diagnostics_text(&sw, ir);
    EXPECT_TRUE(diagnostics_contain(txt, "<script>:1:16: error: Invalid arithmetic expression, left operand of '*' is missing."));
    EXPECT_TRUE(diagnostics_contain(txt,
        "    1 | a = distance(1,*2);\n"
        "      |                ^\n"
        "1 error generated.\n"));

    // Only '-' has a unary form: every other binary operator needs its left operand
    const char* missing_left[] = {
        "a = distance(1,+2);",
        "a = distance(1,/2);",
        "a = distance(1,//2);",
        "a = distance(1,<2);",
        "a = distance(1,==2);",
        "a = *2;",
        "a = 1 * * 2;",
        "a = (*2);",
        "a = not *2;",
        "*2;",
    };
    for (size_t i = 0; i < ARRAY_SIZE(missing_left); ++i) {
        md_script_ir_clear(ir);
        EXPECT_FALSE(md_script_ir_compile_from_source(ir, str_from_cstr(missing_left[i]), &test_mol, NULL));
        EXPECT_TRUE(diagnostics_contain(diagnostics_text(&sw, ir), "left operand of"));
    }

    // A comparison without its right operand says so, where it is
    md_script_ir_clear(ir);
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("a = 1 < ;"), &test_mol, NULL));
    EXPECT_TRUE(diagnostics_contain(diagnostics_text(&sw, ir), "<script>:1:7: error: Invalid comparison, right operand of '<' is missing."));

    // A statement which logged an error is not kept, whatever came back from it. The argument list which failed to
    // parse made this a call to all() without arguments, which compiled
    md_script_ir_clear(ir);
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, STR_LIT("a = all(<2);"), &test_mol, NULL));
    EXPECT_EQ((size_t)1, md_script_ir_num_errors(ir));

    // The unary minus is unaffected
    md_script_ir_clear(ir);
    EXPECT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("a = -2 * 3;\nb = 2 * -3;\nc = 1 - -a;"), &test_mol, NULL));

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// What is wrong with a script is an error of its compilation: never an assertion, nor a crash in a release build
// (where the assertion is gone), nor a script which compiles into something other than what was written
UTEST_F(script, malformed_scripts_are_errors) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);
    md_str_writer_t sw = md_str_writer_create(alloc);

    struct {
        const char* src;
        const char* error;
    } cases[] = {
        // Expressions next to each other, without an operator between them: the second took the place of the first
        {"a = 3 4;",                                "<script>:1:5: error: Invalid syntax, expressions must be separated by an operator"},
        {"a = distance(1 2, 3);",                   "<script>:1:14: error: Invalid syntax, expressions must be separated by an operator"},
        {"a = residue(1) not residue(2);",          "error: Invalid syntax, expressions must be separated by an operator"},
        // An assignment within an expression declared its identifier as often as the expression was checked
        {"a = distance(1, (b = 2));",               "<script>:1:18: error: An assignment is a statement of its own"},
        {"c = contacts(residue(:), cutoff = 5 x = 3);", "error: "},
        // More arguments than the evaluation has room for
        {"v = vec4(1, 2, 3, 4, 5, 6, 7, 8, 9);",    "<script>:1:5: error: Too many arguments (9) in call to 'vec4'"},
        // A context which is not known at compile time was reported, and the statement kept without a type
        {"x = distance(1, 2) in within(3, residue(1));", "error: Right hand side of 'in' must be known at compile time."},
        // Integer division traps on a zero divisor
        {"a = -1 // 0;",                            "<script>:1:5: error: Integer division by zero"},
        {"a = {4, 2} // {2, 0};",                   "error: Integer division by zero"},
        // The int[2] could neither be taken as it was nor converted, and did not fill its part of the array
        {"f = {{1, 2}, {3, 4:5}};",                 "<script>:1:6: error: Incompatible types wihin array construct"},
        // A selection which is empty (the system has no water) where one atom is needed
        {"m = match(water);",                       "error: match: the reference is empty"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_script_ir_clear(ir);
        EXPECT_FALSE(md_script_ir_compile_from_source(ir, str_from_cstr(cases[i].src), mol, NULL));
        EXPECT_FALSE(md_script_ir_valid(ir));
        EXPECT_TRUE(diagnostics_contain(diagnostics_text(&sw, ir), cases[i].error));
    }

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// Visualizing what an editor hovers: the procedures are called without a destination, which only some of them take
UTEST_F(script, visualize_without_destination) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    const char* src =
        "arr = {vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0)};\n"
        "cx = coord_x(residue(1)) in residue(1:2);\n"
        "c = coord(residue(1)) in residue(1:2);\n"
        "r = rmsd(residue(1:2));\n";
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(src), mol, NULL));

    md_script_vis_ctx_t vctx = { .ir = ir, .sys = mol, .state = &mol->reference };
    const char* names[] = {"arr", "cx", "c", "r"};
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        identifier_t* ident = get_identifier(ir, str_from_cstr(names[i]));
        ASSERT_TRUE(ident);
        // An element of an array: the others are skipped, they used to be called as if for the length of their result
        for (int subidx = -1; subidx < 2; ++subidx) {
            md_script_vis_t vis = {0};
            md_script_vis_init(&vis, alloc);
            EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, ident->node), subidx, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT));
            md_script_vis_free(&vis);
        }
        char buf[256];
        md_script_vis_value_text(buf, sizeof(buf), make_vis_ref(ir, ident->node), -1, &vctx);
    }

    // In a context each part has the length fixed at compile time. All of a part is written, which coord_x did not
    // (it wrote one value of each, and asserted), and what it has room for beyond the values of its context is zero
    data_t res1 = {0}, res2 = {0}, both = {0};
    ASSERT_TRUE(eval_expression(&res1, STR_LIT("coord_x(residue(1))"), mol, alloc));
    ASSERT_TRUE(eval_expression(&res2, STR_LIT("coord_x(residue(2))"), mol, alloc));
    ASSERT_TRUE(eval_expression(&both, STR_LIT("coord_x(residue(1)) in residue(1:2)"), mol, alloc));
    const size_t n1 = res1.size / sizeof(float);
    const size_t n2 = res2.size / sizeof(float);
    ASSERT_GT(n1, n2);
    ASSERT_EQ(2 * n1 * sizeof(float), both.size);
    const float* x = (const float*)both.ptr;
    for (size_t i = 0; i < n1; ++i) {
        EXPECT_EQ(((const float*)res1.ptr)[i], x[i]);
        EXPECT_EQ(i < n2 ? ((const float*)res2.ptr)[i] : 0.0f, x[n1 + i]);
    }

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// Scripts which compile, and crashed when they were evaluated or visualized
UTEST_F(script, evaluate_edge_cases) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(4));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)MIN(script_frames(mol), 4);
    md_script_ir_t* ir = md_script_ir_create(alloc);

    const char* srcs[] = {
        // The center of mass of a range, within a context
        "c = com(1:3) in polymer(:); d = distance(com(1:3), 1) in residue(1:2);",
        // Element wise operations on arrays which can be empty: the dims of an empty array made for one element
        "s = within(3.0, residue(1:)); m = max(coord_x(s) + coord_y(s));",
        // The residues of nothing
        "n = count(residue(residue(within(0:0, residue(1)))));",
        // Parts of a context which are empty
        "x = coord(not residue(1)) in residue(1:2); y = coord_x(not residue(1)) in residue(1:2);",
        // An empty array is typed as a scalar (its dims are all 0) and matched with the scalar operations, which have
        // nothing to read: the evaluation fails, with an error
        "e = -coord_z(flatten(water)); f = sqrt(coord_z(flatten(water)));",
        "g = -coord_x(not residue(1)) in residue(1:2);",
    };
    md_script_vis_ctx_t vctx = { .ir = ir, .sys = mol, .state = &mol->reference };
    for (size_t i = 0; i < ARRAY_SIZE(srcs); ++i) {
        md_script_ir_clear(ir);
        EXPECT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(srcs[i]), mol, NULL));
        if (!md_script_ir_valid(ir)) {
            md_script_ir_write_diagnostics(md_writer_stderr(), ir);
            continue;
        }
        if (md_script_ir_property_count(ir) > 0) {
            md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
            ASSERT_NE(NULL, eval);
            md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames);
            md_script_eval_free(eval);
        }
        // What hovering each of its tokens does
        const md_script_vis_token_t* toks = md_script_ir_vis_tokens(ir);
        for (size_t t = 0; t < md_script_ir_num_vis_tokens(ir); ++t) {
            for (int subidx = -1; subidx < 1; ++subidx) {
                md_script_vis_t vis = {0};
                md_script_vis_init(&vis, alloc);
                md_script_vis_eval_ref(&vis, toks[t].ref, subidx, &vctx, MD_SCRIPT_VISUALIZE_DEFAULT);
                md_script_vis_free(&vis);
            }
            char buf[256];
            md_script_vis_value_text(buf, sizeof(buf), toks[t].ref, -1, &vctx);
        }
    }

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// The operands of an operator are taken in the order they are written. FLAG_SYMMETRIC_ARGS lets a match swap them,
// and was set on '-', '/' and '//' with an array and a scalar: '2.0 - {1.0, 2.0}' was {1.0, 2.0} - 2.0
UTEST_F(script, operand_order) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* mol = &utest_fixture->ala;

    struct {
        const char* expr;
        float expected[2];
    } cases[] = {
        {"2.0 - {1.0, 2.0}",    {1.0f, 0.0f}},
        {"{1.0, 2.0} - 2.0",    {-1.0f, 0.0f}},
        {"10.0 / {2.0, 5.0}",   {5.0f, 2.0f}},
        {"{2.0, 5.0} / 10.0",   {0.2f, 0.5f}},
        {"2.0 + {1.0, 2.0}",    {3.0f, 4.0f}},
        {"3.0 * {1.0, 2.0}",    {3.0f, 6.0f}},
        // Integers are taken as floats by '-' and '/', and stay integers for '//'
        {"2 - {1, 2}",          {1.0f, 0.0f}},
        {"{1, 2} - 2",          {-1.0f, 0.0f}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        data_t d = {0};
        ASSERT_TRUE(eval_expression(&d, str_from_cstr(cases[i].expr), mol, alloc));
        ASSERT_EQ((size_t)2, (size_t)type_info_array_len(d.type));
        for (int j = 0; j < 2; ++j) {
            const float v = d.type.base_type == TYPE_INT ? (float)((int*)d.ptr)[j] : ((float*)d.ptr)[j];
            EXPECT_NEAR(cases[i].expected[j], v, 1e-6f);
        }
    }

    struct {
        const char* expr;
        int expected[2];
    } int_cases[] = {
        {"10 // {2, 5}",        {5, 2}},
        {"{10, 20} // 5",       {2, 4}},
    };
    for (size_t i = 0; i < ARRAY_SIZE(int_cases); ++i) {
        data_t d = {0};
        ASSERT_TRUE(eval_expression(&d, str_from_cstr(int_cases[i].expr), mol, alloc));
        ASSERT_EQ(TYPE_INT, d.type.base_type);
        EXPECT_EQ(int_cases[i].expected[0], ((int*)d.ptr)[0]);
        EXPECT_EQ(int_cases[i].expected[1], ((int*)d.ptr)[1]);
    }

    // Volumes: s - v and s / v against v
    data_t v = {0}, sv = {0}, dv = {0};
    ASSERT_TRUE(eval_expression(&v,  STR_LIT("sdf(residue(1), element('H'), 5.0) + 1.0"), mol, alloc));
    ASSERT_TRUE(eval_expression(&sv, STR_LIT("3.0 - (sdf(residue(1), element('H'), 5.0) + 1.0)"), mol, alloc));
    ASSERT_TRUE(eval_expression(&dv, STR_LIT("3.0 / (sdf(residue(1), element('H'), 5.0) + 1.0)"), mol, alloc));
    ASSERT_EQ(v.size, sv.size);
    ASSERT_EQ(v.size, dv.size);
    const float* fv = (const float*)v.ptr;
    for (size_t i = 0; i < v.size / sizeof(float); i += 997) {
        EXPECT_NEAR(3.0f - fv[i], ((const float*)sv.ptr)[i], 1e-5f);
        EXPECT_NEAR(3.0f / fv[i], ((const float*)dv.ptr)[i], 1e-5f);
    }

    // Only what is commutative swaps its operands
    const procedure_t* tables[] = {operators, procedures};
    const size_t table_len[] = {ARRAY_SIZE(operators), ARRAY_SIZE(procedures)};
    const str_t commutative[] = {STR_LIT("+"), STR_LIT("*"), STR_LIT("and"), STR_LIT("or"), STR_LIT("xor"), STR_LIT("min"), STR_LIT("max")};
    for (size_t t = 0; t < ARRAY_SIZE(tables); ++t) {
        for (size_t i = 0; i < table_len[t]; ++i) {
            if (!(tables[t][i].flags & FLAG_SYMMETRIC_ARGS)) continue;
            bool ok = false;
            for (size_t c = 0; c < ARRAY_SIZE(commutative); ++c) {
                ok |= str_eq(tables[t][i].name, commutative[c]);
            }
            if (!ok) printf("FLAG_SYMMETRIC_ARGS on '" STR_FMT "', which is not commutative\n", STR_ARG(tables[t][i].name));
            EXPECT_TRUE(ok);
        }
    }

    md_arena_allocator_destroy(alloc);
}

// An array of volumes is a property of its own: published {n, x, y, z}, its i-th volume what the i-th expression of
// it is as a volume alone. A context makes one, and so does an array of volumes, which stacks them: '{sdf(a), sdf(b)}'
// was one volume twice as long along x, and not a property.
UTEST_F(script, volume_arrays) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(16));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = (uint32_t)MIN(script_frames(mol), 2);
    md_script_ir_t* ir = md_script_ir_create(alloc);

    const char* src =
        // The density of the hydrogens of the system ('out' of the context) around each residue
        "ctx = sdf(residue(1), out element('H'), 5.0) in residue(1:2);\n"
        "ctx1 = sdf(residue(1), out element('H'), 5.0) in residue(1);\n"
        "ctx2 = sdf(residue(1), out element('H'), 5.0) in residue(2);\n"
        "arr = {sdf(residue(1), element('H'), 5.0), sdf(residue(2), element('H'), 5.0)};\n"
        "arr1 = sdf(residue(1), element('H'), 5.0);\n"
        "arr2 = sdf(residue(2), element('H'), 5.0);\n";
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(src), mol, NULL));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_VOLUME, md_script_ir_property_flags(ir, STR_LIT("ctx")) & MD_SCRIPT_PROPERTY_FLAG_VOLUME);
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_VOLUME, md_script_ir_property_flags(ir, STR_LIT("arr")) & MD_SCRIPT_PROPERTY_FLAG_VOLUME);

    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));
    const md_attributes_t* table = md_script_eval_attributes(eval);

    struct { const char* array; const char* member[2]; } cases[] = {
        {"script/ctx", {"script/ctx1", "script/ctx2"}},
        {"script/arr", {"script/arr1", "script/arr2"}},
    };
    const size_t vol_size = (size_t)MD_VOL_DIM * MD_VOL_DIM * MD_VOL_DIM;
    for (size_t c = 0; c < ARRAY_SIZE(cases); ++c) {
        const md_attribute_t* a = md_attributes_find(table, str_from_cstr(cases[c].array));
        ASSERT_TRUE(a);
        ASSERT_EQ(4u, a->format.rank);
        EXPECT_EQ(2u, a->format.shape[0]);
        EXPECT_EQ((uint32_t)MD_VOL_DIM, a->format.shape[1]);
        EXPECT_EQ((uint32_t)MD_VOL_DIM, a->format.shape[3]);
        for (int m = 0; m < 2; ++m) {
            const md_attribute_t* v = md_attributes_find(table, str_from_cstr(cases[c].member[m]));
            ASSERT_TRUE(v);
            ASSERT_EQ(3u, v->format.rank);
            const float* x = (const float*)a->data + m * vol_size;
            const float* y = (const float*)v->data;
            size_t num_diff = 0;
            double sum = 0;
            for (size_t i = 0; i < vol_size; ++i) {
                num_diff += x[i] != y[i];
                sum += y[i];
            }
            EXPECT_EQ((size_t)0, num_diff);
            EXPECT_GT(sum, 0.0);
        }
    }
    md_script_eval_free(eval);

    // Visualized, element m of the array is the structure of its own volume: residue m + 1
    md_script_vis_ctx_t vctx = { .ir = ir, .sys = mol, .state = &mol->reference };
    const char* arrays[] = {"ctx", "arr"};
    for (size_t a = 0; a < ARRAY_SIZE(arrays); ++a) {
        identifier_t* ident = get_identifier(ir, str_from_cstr(arrays[a]));
        ASSERT_TRUE(ident);
        for (int m = 0; m < 2; ++m) {
            md_script_vis_t vis = {0};
            md_script_vis_init(&vis, alloc);
            EXPECT_TRUE(md_script_vis_eval_ref(&vis, make_vis_ref(ir, ident->node), m, &vctx, MD_SCRIPT_VISUALIZE_SDF));
            ASSERT_EQ((size_t)1, md_array_size(vis.sdf.structures));
            const md_urange_t res = md_system_residue_atom_range(mol, (size_t)m);
            EXPECT_EQ((size_t)(res.end - res.beg), md_bitset_count(vis.sdf.structures[0]));
            EXPECT_TRUE(md_bitset_test(vis.sdf.structures[0], res.beg));
            md_script_vis_free(&vis);
        }
    }

    // No more than so many volumes as one property
    md_script_ir_clear(ir);
    char big[1024];
    int len = snprintf(big, sizeof(big), "v = sdf(residue(1), element('H'), 5.0);\nw = {v");
    for (int i = 1; i < MAX_PROPERTY_VOLUMES + 1; ++i) len += snprintf(big + len, sizeof(big) - len, ", v");
    snprintf(big + len, sizeof(big) - len, "};");
    EXPECT_FALSE(md_script_ir_compile_from_source(ir, str_from_cstr(big), mol, NULL));
    md_str_writer_t sw = md_str_writer_create(alloc);
    EXPECT_TRUE(diagnostics_contain(diagnostics_text(&sw, ir), "An array of 33 volumes"));

    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}

// ### BONDS ###

// Acetate with Mg at 2.08 Å from one oxygen, a water, and two Cu 2.45 Å apart (as in test_bonds.c)
static const char bonds_site_pdb[] =
    "HETATM    1  C1  ACT A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    "HETATM    2  C2  ACT A   1       1.520   0.000   0.000  1.00  0.00           C\n"
    "HETATM    3  O1  ACT A   1       2.150   1.080   0.000  1.00  0.00           O\n"
    "HETATM    4  O2  ACT A   1       2.150  -1.080   0.000  1.00  0.00           O\n"
    "HETATM    5  MG  MG  A   2       3.198   2.877   0.000  1.00  0.00          MG\n"
    "HETATM    6  O   HOH A   3       0.000   6.000   0.000  1.00  0.00           O\n"
    "HETATM    7  H1  HOH A   3       0.760   6.590   0.000  1.00  0.00           H\n"
    "HETATM    8  H2  HOH A   3      -0.760   6.590   0.000  1.00  0.00           H\n"
    "HETATM    9  CU  CU  A   4      20.000   0.000   0.000  1.00  0.00          CU\n"
    "HETATM   10  CU  CU  A   5      22.450   0.000   0.000  1.00  0.00          CU\n"
    "END\n";

// The number of bonds of a bond set identifier, -1 if it is not one
static int bondset_count_of(md_script_ir_t* ir, const char* name, md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    data_t d = {0};
    if (!mirror_eval(&d, ir, name, sys, state, alloc)) return -1;
    if (d.type.base_type != TYPE_BONDSET) {
        printf("  '%s' is not a bond set\n", name);
        return -1;
    }
    int n = 0;
    const bondset_t* b = d.ptr;
    for (size_t i = 0; i < element_count(d); ++i) n += (int)md_bitset_count(b[i].bonds) + (int)b[i].num_extra;
    return n;
}

static float float_of(md_script_ir_t* ir, const char* name, md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    data_t d = {0};
    if (!mirror_eval(&d, ir, name, sys, state, alloc) || d.type.base_type != TYPE_FLOAT) return -1.0f;
    return as_float(d);
}

UTEST(script, bondsets) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, str_from_cstr(bonds_site_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    // C1-C2, C2-O1, C2-O2, Mg-O1, O-H1, O-H2, Cu-Cu
    ASSERT_EQ(7, (int)sys.bond.count);
    // A hydrogen bond a user marks, from a water hydrogen to O2
    md_system_bond_insert(&sys, 6, 3, md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_HYDROGEN), MD_BOND_ORIGIN_USER));
    ASSERT_EQ(1, (int)sys.bond.noncovalent_count);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    const char* src =
        "first_r = bond(residue(:))[1];"
        "two_r   = bond(residue(:))[1:2];"
        "pick    = {bond(element('Mg') or within(3, element('Mg')) or element('Cu')), bond()}[1];"
        "all_b   = bond();"
        "cov     = bond(kind='covalent');"
        "coord_b = bond(kind='coordinate');"
        "hb      = bond(kind='hydrogen');"
        "nc      = bond(kind='noncovalent');"
        "inter_b = bond(kind='interaction');"
        "graph   = bond(kind={'covalent', 'coordinate'});"
        "user    = bond(origin='user');"
        "ordered = bond(order=0:4);"
        "water_b = bond(resname('HOH'));"
        "c_o     = bond(element('C'), element('O'));"
        "at_mg   = bond(b=element('Mg'));"
        "per_res = bond(residue(:));"
        "none_b  = cov and coord_b;"
        "either  = cov or coord_b;"
        "other   = not cov;"
        "odd     = graph xor cov;"
        "mg_side = atom(coord_b);"
        "ends_r  = atom(bond(residue(:)));"
        "n_all   = count(all_b);"
        "l       = length(coord_b);"
        "in_ctx  = bond() in resname('ACT');"
        "n_first = count(first_r);"
        "by_atom = bond(atom(2), atom(3));"
        "by_idx  = bond(2, 3);"
        "by_rev  = bond(3, 2);"
        "by_arr  = bond(2, {3, 4});";
    const bool compiled = md_script_ir_compile_from_source(ir, str_from_cstr(src), &sys, NULL);
    for (size_t i = 0; i < md_script_ir_num_errors(ir); ++i) printf("%.*s\n", STR_ARG(md_script_ir_errors(ir)[i].text));
    ASSERT_TRUE(compiled);
    ASSERT_TRUE(md_script_ir_valid(ir));

    EXPECT_EQ(8, bondset_count_of(ir, "all_b", &sys, &st, alloc));
    EXPECT_EQ(6, bondset_count_of(ir, "cov", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "coord_b", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "hb", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "nc", &sys, &st, alloc));
    EXPECT_EQ(0, bondset_count_of(ir, "inter_b", &sys, &st, alloc));
    EXPECT_EQ(7, bondset_count_of(ir, "graph", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "user", &sys, &st, alloc));
    // The order is of the bonds of the graph: the hydrogen bond has none
    EXPECT_EQ(7, bondset_count_of(ir, "ordered", &sys, &st, alloc));
    EXPECT_EQ(2, bondset_count_of(ir, "water_b", &sys, &st, alloc));
    EXPECT_EQ(2, bondset_count_of(ir, "c_o", &sys, &st, alloc));
    // A bond by the indices of its atoms, as the context menu writes it: in either order, the ints taken as atoms
    EXPECT_EQ(1, bondset_count_of(ir, "by_atom", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "by_idx", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "by_rev", &sys, &st, alloc));
    EXPECT_EQ(2, bondset_count_of(ir, "by_arr", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "at_mg", &sys, &st, alloc));
    // One set per residue: the bonds within each, so neither Mg-O1 nor the hydrogen bond (between residues) nor Cu-Cu
    data_t per_res = {0};
    ASSERT_TRUE(mirror_eval(&per_res, ir, "per_res", &sys, &st, alloc));
    EXPECT_EQ((int)sys.residue.count, (int)element_count(per_res));
    EXPECT_EQ(5, bondset_count_of(ir, "per_res", &sys, &st, alloc));
    EXPECT_EQ(0, bondset_count_of(ir, "none_b", &sys, &st, alloc));
    EXPECT_EQ(7, bondset_count_of(ir, "either", &sys, &st, alloc));
    EXPECT_EQ(2, bondset_count_of(ir, "other", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "odd", &sys, &st, alloc));
    EXPECT_EQ(3, bondset_count_of(ir, "in_ctx", &sys, &st, alloc));
    // A subscript is a set of its own, not a view into the array it was taken from: evaluated when the script is
    // compiled, the array is gone by the time it is read
    EXPECT_EQ(3, bondset_count_of(ir, "first_r", &sys, &st, alloc));
    EXPECT_EQ(3, bondset_count_of(ir, "two_r", &sys, &st, alloc));
    EXPECT_EQ(3.0f, float_of(ir, "n_first", &sys, &st, alloc));
    {
        // Evaluated with scratch memory of its own, which is written over once it is done: nothing of the result may
        // be in it. Of the frame (within), so that it is evaluated here rather than when the script is compiled, and so is
        // the set it is taken from.
        md_allocator_i* scratch = md_vm_arena_create(GIGABYTES(1));
        identifier_t* ident = get_identifier(ir, STR_LIT("pick"));
        ASSERT_TRUE(ident && ident->node);
        eval_context_t ctx = { .ir = ir, .sys = &sys, .temp_alloc = scratch, .alloc = alloc, .cur_state = &st, .ref_state = &sys.reference };
        data_t d = {0};
        ASSERT_TRUE(evaluate_node_alloc(&d, ident->node, &ctx, alloc));
        md_vm_arena_reset(scratch);
        MEMSET(md_vm_arena_push(scratch, MEGABYTES(4)), 0xFF, MEGABYTES(4));
        ASSERT_EQ(TYPE_BONDSET, d.type.base_type);
        ASSERT_EQ(1, (int)element_count(d));
        // Two bonds apart from one another, Mg-O1 and Cu-Cu: words of their own
        const bondset_t* b = (const bondset_t*)d.ptr;
        EXPECT_EQ(2, (int)md_bitset_count(b->bonds));
        EXPECT_TRUE(md_bitset_test(b->bonds, 3));
        EXPECT_TRUE(md_bitset_test(b->bonds, 6));
        EXPECT_EQ(0, (int)b->num_extra);
        md_vm_arena_destroy(scratch);
    }
    {
        // As a filter names it: the value the compilation stored
        const md_filter_context_t fctx = { .sys = &sys, .state = &st, .ir = ir };
        md_array(md_atom_pair_t) pairs = 0;
        EXPECT_TRUE(md_filter_bonds(&pairs, STR_LIT("first_r"), &fctx, alloc).success);
        EXPECT_EQ(3, (int)md_array_size(pairs));
        for (size_t i = 0; i < md_array_size(pairs); ++i) {
            EXPECT_TRUE(0 <= pairs[i].idx[0] && pairs[i].idx[0] < 4 && 0 <= pairs[i].idx[1] && pairs[i].idx[1] < 4);
        }
    }

    data_t atoms = {0};
    ASSERT_TRUE(mirror_eval(&atoms, ir, "mg_side", &sys, &st, alloc));
    ASSERT_EQ(TYPE_ATOMSET, atoms.type.base_type);
    EXPECT_EQ(2, (int)md_bitset_count(*as_bitset(atoms)));
    EXPECT_TRUE(md_bitset_test(*as_bitset(atoms), 2));
    EXPECT_TRUE(md_bitset_test(*as_bitset(atoms), 4));
    // atom(b) of several sets: one set of atoms each, the atoms of the bonds within each residue
    data_t ends = {0};
    ASSERT_TRUE(mirror_eval(&ends, ir, "ends_r", &sys, &st, alloc));
    ASSERT_EQ(TYPE_ATOMSET, ends.type.base_type);
    ASSERT_EQ((int)sys.residue.count, (int)element_count(ends));
    EXPECT_EQ(4, (int)md_bitset_count(as_bitset(ends)[0]));
    EXPECT_EQ(0, (int)md_bitset_count(as_bitset(ends)[1]));
    EXPECT_EQ(3, (int)md_bitset_count(as_bitset(ends)[2]));

    // Bonds and atoms do not mix within an expression; the error says how to write it
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = resname('HOH') and bond();", "atom(b)"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = within(3.0, bond());", "atom(b)"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = resname('HOH') in bond();", "atom(b)"));
    EXPECT_EQ(8.0f, float_of(ir, "n_all", &sys, &st, alloc));
    data_t len = {0};
    ASSERT_TRUE(mirror_eval(&len, ir, "l", &sys, &st, alloc));
    ASSERT_EQ(1, (int)element_count(len));
    EXPECT_NEAR(2.08f, as_float_arr(len)[0], 0.01f);

    // Shown as the atoms of its bonds and a line along each
    md_script_vis_t vis = {0};
    md_script_vis_init(&vis, alloc);
    const md_script_vis_ctx_t vis_ctx = { .sys = &sys, .state = &st };
    EXPECT_TRUE(md_script_vis_eval_string(&vis, STR_LIT("v = bond(kind='coordinate');"), &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    EXPECT_EQ(2, (int)md_bitset_count(vis.atom_mask));
    EXPECT_EQ(2, (int)md_array_size(vis.lines));
    md_script_vis_clear(&vis);
    EXPECT_TRUE(md_script_vis_eval_string(&vis, STR_LIT("v = formed(bond_delta());"), &vis_ctx, MD_SCRIPT_VISUALIZE_DEFAULT));
    md_script_vis_free(&vis);

    // Hovering says what it is
    char text[256];
    data_t hv = {0};
    ASSERT_TRUE(mirror_eval(&hv, ir, "all_b", &sys, &st, alloc));
    md_str_writer_t sw = md_str_writer_create(alloc);
    summarize_data(&sw.w, hv, -1);
    snprintf(text, sizeof(text), "%.*s", STR_ARG(md_str_writer_to_str(&sw)));
    EXPECT_STREQ("8 bonds", text);

    // A bond set is no selection: it has to be made one
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = bond() and element('C');", "Could not find supporting operator"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = bond(kind='weak');", "Unknown bond kind 'weak'"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = bond(origin='hand');", "Unknown bond origin 'hand'"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "x = bond(resonance='aromatic_ish');", "Unknown bond resonance"));
    md_vm_arena_destroy(alloc);
}

// Two waters, a hydrogen of the first pointing at the oxygen of the second (as in test_bonds.c)
static const char bonds_water_pair_pdb[] =
    "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
    "HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H\n"
    "HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H\n"
    "HETATM    4  O   HOH A   2       2.800   0.000   0.000  1.00  0.00           O\n"
    "HETATM    5  H1  HOH A   2       3.040   0.930   0.000  1.00  0.00           H\n"
    "HETATM    6  H2  HOH A   2       3.040  -0.930   0.000  1.00  0.00           H\n"
    "END\n";

UTEST(script, bond_delta) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t ref = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &ref, str_from_cstr(bonds_water_pair_pdb), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &ref, MD_UTIL_INFER_ALL));
    ASSERT_EQ(4, (int)sys.bond.count);

    md_script_ir_t* ir = md_script_ir_create(alloc);
    const char* src =
        "d       = bond_delta();"
        "f       = formed(d);"
        "b       = broken(d);"
        "e       = bond(d);"
        "eb      = bond(d, b=resid(1));"
        "eo      = bond(d, order=1);"
        "ek      = bond(d, kind='covalent', origin='inferred');"
        "nf      = count(f);"
        "nb      = count(b);"
        "d2      = bond_delta(resid(2));"
        "f2      = formed(d2);"
        "e2      = bond(d2);"
        "pinned  = broken(bond_delta(breakable='user'));"
        "narrow  = formed(bond_delta(tolerance=0.3));";
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(src), &sys, NULL));
    ASSERT_TRUE(md_script_ir_valid(ir));
    EXPECT_EQ(md_script_ir_property_flags(ir, STR_LIT("nf")), MD_SCRIPT_PROPERTY_FLAG_TEMPORAL);

    // In the reference nothing changes
    EXPECT_EQ(0, bondset_count_of(ir, "f", &sys, &ref, alloc));
    EXPECT_EQ(0, bondset_count_of(ir, "b", &sys, &ref, alloc));
    EXPECT_EQ(4, bondset_count_of(ir, "e", &sys, &ref, alloc));
    EXPECT_EQ(2, bondset_count_of(ir, "eb", &sys, &ref, alloc));
    EXPECT_EQ(4, bondset_count_of(ir, "eo", &sys, &ref, alloc));
    EXPECT_EQ(4, bondset_count_of(ir, "ek", &sys, &ref, alloc));

    // The proton moved to the second water
    md_system_state_t st = ref;
    st.xyz = md_alloc(alloc, sizeof(vec3_t) * ref.num_atoms);
    MEMCPY(st.xyz, ref.xyz, sizeof(vec3_t) * ref.num_atoms);
    st.xyz[1] = vec3_set(1.80f, 0.0f, 0.0f);
    EXPECT_EQ(1, bondset_count_of(ir, "f", &sys, &st, alloc));
    EXPECT_EQ(1, bondset_count_of(ir, "b", &sys, &st, alloc));
    EXPECT_EQ(4, bondset_count_of(ir, "e", &sys, &st, alloc));
    // With an atom of the first water: its O-H2, and the bond its H1 formed to the second
    EXPECT_EQ(2, bondset_count_of(ir, "eb", &sys, &st, alloc));
    // The bond formed is of no known order, the others single
    EXPECT_EQ(3, bondset_count_of(ir, "eo", &sys, &st, alloc));
    // Inferred and covalent, as the bonds of the waters
    EXPECT_EQ(4, bondset_count_of(ir, "ek", &sys, &st, alloc));
    EXPECT_EQ(1.0f, float_of(ir, "nf", &sys, &st, alloc));
    EXPECT_EQ(1.0f, float_of(ir, "nb", &sys, &st, alloc));
    // The scope of the second water: its two bonds, and the one formed to it
    EXPECT_EQ(1, bondset_count_of(ir, "f2", &sys, &st, alloc));
    EXPECT_EQ(3, bondset_count_of(ir, "e2", &sys, &st, alloc));
    // Inferred bonds do not break when only a user's may
    EXPECT_EQ(0, bondset_count_of(ir, "pinned", &sys, &st, alloc));
    // A wider band: the proton at 1.0 Å from O2 is not yet within 0.7 of its cutoff
    EXPECT_EQ(0, bondset_count_of(ir, "narrow", &sys, &st, alloc));

    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "d = bond_delta(tolerance=0.5);", "tolerance"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "d = bond_delta(breakable='hand');", "Unknown bond origin 'hand'"));
    EXPECT_TRUE(contacts_compile_fails_with(&sys, alloc, "e = bond(bond_delta(), kind='ionic');", "Unknown bond kind 'ionic'"));
    md_vm_arena_destroy(alloc);
}

// ### PERIODIC VALUES ###
// What the compilation says of the period of an expression (see PERIODIC VALUES in md_script.c), what the evaluation
// publishes of it, and the circular summary over a population of periodic values.

static period_t ident_period(md_script_ir_t* ir, const char* name) {
    const identifier_t* ident = get_identifier(ir, str_from_cstr(name));
    return (ident && ident->data) ? ident->data->period : (period_t){0};
}

#define EXPECT_PERIOD(ir, name, b, e, s) do { \
    const period_t _p = ident_period(ir, name); \
    EXPECT_TRUE(period_valid(_p)); \
    EXPECT_NEAR((float)(b), _p.beg, 1.0e-4f); \
    EXPECT_NEAR((float)(e), _p.end, 1.0e-4f); \
    EXPECT_NEAR((float)(s), _p.seam, 1.0e-4f); \
} while (0)

#define EXPECT_NOT_PERIODIC(ir, name) EXPECT_FALSE(period_valid(ident_period(ir, name)))

UTEST_F(script, periodic_propagation) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    md_script_ir_t* ir = md_script_ir_create(alloc);

    // phi of residue 2 (C1 N2 CA2 C2), psi of residue 2 (N2 CA2 C2 N3) and phi of residue 3
    const char* src =
        "d    = dihedral(11, 13, 15, 21);\n"
        "deg  = dihedral(11, 13, 15, 21) * 180 / PI;\n"
        "half = dihedral(11, 13, 15, 21) / 2;\n"
        "sh   = dihedral(11, 13, 15, 21) + 1.0;\n"
        "lsh  = 1.0 + dihedral(11, 13, 15, 21);\n"
        "rsh  = 1.0 - dihedral(11, 13, 15, 21);\n"
        "neg  = -dihedral(11, 13, 15, 21);\n"
        "diff = dihedral(11, 13, 15, 21) - dihedral(13, 15, 21, 23);\n"
        "sum2 = dihedral(11, 13, 15, 21) + dihedral(13, 15, 21, 23);\n"
        "arr  = {dihedral(11, 13, 15, 21), dihedral(21, 23, 25, 31)};\n"
        "elem = arr[2];\n"
        "lo   = min(arr);\n"
        "tot  = sum(arr);\n"
        "ref  = d;\n"
        "at   = atan2(distance(1, 5), distance(5, 11));\n"
        "mix  = {dihedral(11, 13, 15, 21), angle(11, 13, 15)};\n"
        "ab   = abs(dihedral(11, 13, 15, 21));\n"
        "co   = cos(dihedral(11, 13, 15, 21));\n"
        "sq   = dihedral(11, 13, 15, 21) * dihedral(13, 15, 21, 23);\n"
        "inv  = 1.0 / dihedral(11, 13, 15, 21);\n"
        "an   = angle(11, 13, 15);\n"
        "dist = distance(1, 5);\n";
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(src), mol, NULL));

    // A dihedral is given in [-180, 180[ and cut at 0
    EXPECT_PERIOD(ir, "d",    -PI, PI, 0);
    // Scaled: degrees, and half angles
    EXPECT_PERIOD(ir, "deg",  -180, 180, 0);
    EXPECT_PERIOD(ir, "half", -PI / 2, PI / 2, 0);
    // Moved along by a constant, on either side
    EXPECT_PERIOD(ir, "sh",   -PI + 1, PI + 1, 1);
    EXPECT_PERIOD(ir, "lsh",  -PI + 1, PI + 1, 1);
    EXPECT_PERIOD(ir, "rsh",  1 - PI, 1 + PI, 1);
    EXPECT_PERIOD(ir, "neg",  -PI, PI, 0);
    // Two angles: about the difference (sum) of the middles, cut half a period from it
    EXPECT_PERIOD(ir, "diff", -PI, PI, -PI);
    EXPECT_PERIOD(ir, "sum2", -PI, PI, -PI);
    // Arrays of them, their elements, the least of them, references
    EXPECT_PERIOD(ir, "arr",  -PI, PI, 0);
    EXPECT_PERIOD(ir, "elem", -PI, PI, 0);
    EXPECT_PERIOD(ir, "lo",   -PI, PI, 0);
    EXPECT_PERIOD(ir, "tot",  -PI, PI, 0);
    EXPECT_PERIOD(ir, "ref",  -PI, PI, 0);
    // atan2 is cut where it cuts
    EXPECT_PERIOD(ir, "at",   -PI, PI, -PI);

    // What no longer lies on a circle
    EXPECT_NOT_PERIODIC(ir, "mix");
    EXPECT_NOT_PERIODIC(ir, "ab");
    EXPECT_NOT_PERIODIC(ir, "co");
    EXPECT_NOT_PERIODIC(ir, "sq");
    EXPECT_NOT_PERIODIC(ir, "inv");
    EXPECT_NOT_PERIODIC(ir, "an");
    EXPECT_NOT_PERIODIC(ir, "dist");

    // -x is in the unit of x
    const identifier_t* neg = get_identifier(ir, STR_LIT("neg"));
    ASSERT_TRUE(neg && neg->data);
    EXPECT_TRUE(md_unit_equal(md_unit_radian(), neg->data->unit[1]));

    // The flags of the properties
    const md_script_property_flags_t temporal_periodic = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_PERIODIC;
    EXPECT_EQ(temporal_periodic, md_script_ir_property_flags(ir, STR_LIT("d")));
    EXPECT_EQ(temporal_periodic, md_script_ir_property_flags(ir, STR_LIT("arr")));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_TEMPORAL, md_script_ir_property_flags(ir, STR_LIT("ab")));
    EXPECT_EQ(MD_SCRIPT_PROPERTY_FLAG_TEMPORAL, md_script_ir_property_flags(ir, STR_LIT("an")));

    // Published beside the values, in their unit; nothing for what is not periodic
    const uint32_t num_frames = script_frames(mol);
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    const md_attributes_t* attributes = md_script_eval_attributes(eval);
    const md_attribute_t* values = md_attributes_find(attributes, STR_LIT("script/deg"));
    const md_attribute_t* period = md_attributes_find(attributes, STR_LIT("script/deg/period"));
    const md_attribute_t* seam   = md_attributes_find(attributes, STR_LIT("script/deg/seam"));
    ASSERT_TRUE(values && period && seam);
    EXPECT_EQ(0u, period->format.rank);
    EXPECT_EQ(2u, period->format.components);
    EXPECT_EQ(0u, seam->format.rank);
    EXPECT_EQ(1u, seam->format.components);
    EXPECT_TRUE(md_unit_equal(values->unit, period->unit));
    EXPECT_TRUE(md_unit_equal(values->unit, seam->unit));
    float pv[2] = {0}, sv = 1;
    EXPECT_EQ(2, md_attribute_extract_f32(pv, 2, period, md_attribute_slice_all(), md_unit_none()));
    EXPECT_EQ(1, md_attribute_extract_f32(&sv, 1, seam, md_attribute_slice_all(), md_unit_none()));
    EXPECT_NEAR(-180.0f, pv[0], 1.0e-3f);
    EXPECT_NEAR( 180.0f, pv[1], 1.0e-3f);
    EXPECT_NEAR(   0.0f, sv,    1.0e-3f);
    EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/ab/period")) == NULL);
    EXPECT_TRUE(md_attributes_find(attributes, STR_LIT("script/ab/seam")) == NULL);

    md_script_eval_free(eval);
    md_script_ir_free(ir);
}

// The phases of puckers: Cremer-Pople phi2 cut at 0, the pseudorotation of a furanose at 270 degrees
UTEST(script, periodic_pucker_phases) {
    plane_sys_t ade = {0};
    ASSERT_TRUE(plane_sys_init(&ade, "A"));
    md_script_ir_t* ir = md_script_ir_create(ade.alloc);
    ASSERT_TRUE(plane_compile(ir, "p = pucker(ring(kind='sugar')); P = pseudorotation(p); f = phase(p); t = theta(pucker(ring(size=6))); q = amplitude(p);", &ade));
    EXPECT_PERIOD(ir, "P", 0, TWO_PI, 1.5 * PI);
    EXPECT_PERIOD(ir, "f", 0, TWO_PI, 0);
    EXPECT_NOT_PERIODIC(ir, "t");
    EXPECT_NOT_PERIODIC(ir, "q");
    md_script_ir_free(ir);
    md_vm_arena_destroy(ade.alloc);
}

static double circ_diff(double a, double b, double len) {
    const double d = period_wrap(a - b, -0.5 * len, len);
    return fabs(d);
}

UTEST(script, circular_summary) {
    float scratch[16];
    float mn, mx, mean, var;
    const period_t deg = period_make(-180, 180, 0);

    // Across the end of the period: about 180, not 0
    {
        const float v[] = {179, -179};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 2, deg, scratch);
        EXPECT_LT(circ_diff(mean, 180, 360), 1.0e-3);
        EXPECT_NEAR(1.0f, sqrtf(var), 1.0e-3f);     // two values 2 degrees apart, as a standard deviation would have it
        EXPECT_NEAR(179.0f, mn, 1.0e-4f);
        EXPECT_NEAR(181.0f, mx, 1.0e-4f);           // past the end: the arc crosses it
    }
    // Within it: as the least and the greatest, and the mean as the linear one
    {
        const float v[] = {10, 30, 20};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 3, deg, scratch);
        EXPECT_NEAR(20.0f, mean, 1.0e-3f);
        EXPECT_NEAR(10.0f, mn, 1.0e-4f);
        EXPECT_NEAR(30.0f, mx, 1.0e-4f);
        // Close to the linear standard deviation (8.165 degrees) when the values are close together
        const double R = (1.0 + 2.0 * cos(DEG_TO_RAD(10.0))) / 3.0;
        EXPECT_NEAR((float)RAD_TO_DEG(sqrt(-2.0 * log(R))), sqrtf(var), 1.0e-3f);
        EXPECT_NEAR(8.165f, sqrtf(var), 0.05f);
    }
    // Given outside of the period, a whole number of periods away
    {
        const float v[] = {370, -350, 10};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 3, deg, scratch);
        EXPECT_NEAR(10.0f, mean, 1.0e-3f);
        EXPECT_NEAR(10.0f, mn, 1.0e-3f);
        EXPECT_NEAR(10.0f, mx, 1.0e-3f);
        EXPECT_NEAR(0.0f, var, 1.0e-3f);
    }
    // All around the circle: the deviation covers the whole period, and so does the arc, but for the widest gap
    {
        const float v[] = {0, 90, 180, -90};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 4, deg, scratch);
        EXPECT_NEAR(180.0f, sqrtf(var), 1.0e-3f);
        EXPECT_NEAR(270.0f, mx - mn, 1.0e-3f);
        EXPECT_TRUE(-180.0f <= mn && mn < 180.0f);
    }
    // One value
    {
        const float v[] = {-42};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 1, deg, scratch);
        EXPECT_NEAR(-42.0f, mean, 1.0e-3f);
        EXPECT_EQ(-42.0f, mn);
        EXPECT_EQ(-42.0f, mx);
        EXPECT_NEAR(0.0f, var, 1.0e-6f);
    }
    // A period of another length and place: [0, 2 PI[
    {
        const period_t rad = period_make(0, TWO_PI, 0);
        const float v[] = {(float)DEG_TO_RAD(350.0), (float)DEG_TO_RAD(20.0)};
        compute_circular_summary(&mn, &mx, &mean, &var, v, 2, rad, scratch);
        EXPECT_LT(circ_diff(mean, DEG_TO_RAD(5.0), TWO_PI), 1.0e-4);
        EXPECT_NEAR((float)DEG_TO_RAD(350.0), mn, 1.0e-4f);
        EXPECT_NEAR((float)DEG_TO_RAD(380.0), mx, 1.0e-4f);
    }
}

// The summary over the population of a periodic property, frame by frame, against circular statistics done here
UTEST_F(script, periodic_summary_over_trajectory) {
    md_allocator_i* alloc = md_arena_allocator_create(utest_fixture->arena, MEGABYTES(1));
    md_system_t* mol = &utest_fixture->ala;
    const uint32_t num_frames = script_frames(mol);
    ASSERT_GT(num_frames, 1u);
    md_script_ir_t* ir = md_script_ir_create(alloc);

    // psi of residues 2 to 6, and the same in degrees
    const char* src =
        "psi = {dihedral(13, 15, 21, 23), dihedral(23, 25, 31, 33), dihedral(33, 35, 41, 43), dihedral(43, 45, 51, 53), dihedral(53, 55, 61, 63)};\n"
        "psi_deg = psi * 180 / PI;\n";
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, str_from_cstr(src), mol, NULL));
    md_script_eval_t* eval = md_script_eval_create(num_frames, ir, alloc);
    ASSERT_NE(NULL, eval);
    ASSERT_TRUE(md_script_eval_frame_range(eval, ir, mol, SCRIPT_RUN, 0, num_frames));

    const char* names[2] = {"psi", "psi_deg"};
    const double lens[2] = {TWO_PI, 360.0};
    for (int n = 0; n < 2; ++n) {
        char path[64];
        const md_attributes_t* attributes = md_script_eval_attributes(eval);
        snprintf(path, sizeof(path), "script/%s", names[n]);
        const md_attribute_t* attr = md_attributes_find(attributes, str_from_cstr(path));
        snprintf(path, sizeof(path), "script/%s/mean", names[n]);
        const md_attribute_t* mean = md_attributes_find(attributes, str_from_cstr(path));
        snprintf(path, sizeof(path), "script/%s/variance", names[n]);
        const md_attribute_t* var  = md_attributes_find(attributes, str_from_cstr(path));
        snprintf(path, sizeof(path), "script/%s/extent", names[n]);
        const md_attribute_t* ext  = md_attributes_find(attributes, str_from_cstr(path));
        ASSERT_TRUE(attr && mean && var && ext);
        const uint32_t pop = attr->format.shape[1];
        ASSERT_EQ(5u, pop);

        const double len = lens[n];
        const double beg = -0.5 * len;
        const double to_angle = TWO_PI / len;
        for (uint32_t f = 0; f < num_frames; ++f) {
            const float* v = (const float*)attr->data + (size_t)f * pop;
            double c = 0, s = 0;
            for (uint32_t k = 0; k < pop; ++k) {
                c += cos(v[k] * to_angle);
                s += sin(v[k] * to_angle);
            }
            c /= pop;
            s /= pop;
            const double R = sqrt(c * c + s * s);
            const double ref_mean = atan2(s, c) / to_angle;
            const float m = ((const float*)mean->data)[f];
            EXPECT_TRUE(beg <= m && m < beg + len);
            EXPECT_LT(circ_diff(m, ref_mean, len), 1.0e-3 * len);
            const double ref_sd = sqrt(MIN(-2.0 * log(R), PI * PI)) / to_angle;
            EXPECT_NEAR((float)ref_sd, sqrtf(((const float*)var->data)[f]), (float)(1.0e-3 * len));

            // Every value lies on the arc, which is no longer than it has to be: it begins and ends at values
            const float* e = (const float*)ext->data + 2 * (size_t)f;
            EXPECT_LE(e[0], e[1]);
            EXPECT_LT(e[1] - e[0], len);
            bool at_beg = false, at_end = false;
            for (uint32_t k = 0; k < pop; ++k) {
                const double w = period_wrap(v[k], e[0], len);
                EXPECT_LE(w, e[1] + 1.0e-4 * len);
                at_beg |= circ_diff(v[k], e[0], len) < 1.0e-4 * len;
                at_end |= circ_diff(v[k], e[1], len) < 1.0e-4 * len;
            }
            EXPECT_TRUE(at_beg && at_end);
        }
    }

    md_script_eval_free(eval);
    md_script_ir_free(ir);
}

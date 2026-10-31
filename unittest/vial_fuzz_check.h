#pragma once

// What the VIAL front end promises for any input at all, as one check. The built-in fuzzer (test_vial_fuzz.c) and
// the libFuzzer target (fuzz_vial.c) both run it; it returns NULL when every promise holds and otherwise what broke.
//
// Tokens       spans in bounds and in order, the last one EOF; every byte no token covers is layout
// Tree         node 0 NONE, the root a MODULE over the whole source, every child id valid and used once, children
//              inside their parent and in source order, the node count linear in the token count
// Diagnostics  spans in bounds, in source order, a known code, a message
// The API      md_vial_compile gives the same tokens and diagnostics as the stages run by hand, twice over
// Lowering     the array IR verifies, whatever was lowered of it
// Evaluation   with a system and its coordinates, the program evaluates for one frame without crashing, within a limit
//              on the size of values, and every value it computes has the dtype and rank of its type

#include <md_vial.h>
#include <vial/vial_internal.h>
#include <vial/ast.h>
#include <vial/ir.h>
#include <vial/interp.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>

#include <stdio.h>
#include <string.h>

// The codes the front end may give
static inline bool vial_fuzz__front_end_code(const char* code) {
    static const char* codes[] = {"E0101", "E0102", "E0103", "E0104", "E0502", "E0701", "E0801", "E0803"};
    for (size_t i = 0; i < sizeof(codes) / sizeof(codes[0]); ++i) {
        if (memcmp(code, codes[i], 6) == 0) return true;
    }
    return false;
}

// The codes the checker gives so far
static inline bool vial_fuzz__checker_code(const char* code) {
    static const char* codes[] = {"E0201", "W0201", "E0305", "E0312", "E0701", "E0801", "E0802", "E0803",
                                  "E0301", "E0303", "E0304", "E0308", "E0319", "W0501",
                                  "E0202", "E0203", "E0204", "E0205", "E0206", "E0207", "E0313", "E0316", "W0601", "W0602",
                                  "E0302", "E0318", "W0608", "E0307", "E0309", "E0310", "E0311", "W0405", "W0603", "E0314", "W0607", "E0315", "E0317", "E0401", "E0403", "E0404", "E0405", "W0404", "W0406", "E0303", "E0306", "W0605"};
    for (size_t i = 0; i < sizeof(codes) / sizeof(codes[0]); ++i) {
        if (memcmp(code, codes[i], 6) == 0) return true;
    }
    return false;
}

static inline bool vial_fuzz__layout_byte(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\\';
}

static inline const char* vial_fuzz__check_tree(const vial_tree_t* t, uint32_t len, size_t num_tokens, md_allocator_i* alloc) {
    const size_t num_nodes = md_array_size(t->nodes);
    const size_t num_children = md_array_size(t->children);
    if (num_nodes < 2) return "the tree has no root";
    if (t->nodes[0].kind != VIAL_N_NONE || t->nodes[0].child_count != 0) return "node 0 is not NONE";
    if (t->root == 0 || t->root >= num_nodes) return "the root is not a node";
    if (t->nodes[t->root].kind != VIAL_N_MODULE) return "the root is not a MODULE";
    if (t->nodes[t->root].span.beg != 0 || t->nodes[t->root].span.end != len) return "the MODULE does not span the source";
    if (num_nodes > 8 * num_tokens + 16) return "more nodes than the tokens account for";

    uint8_t* seen = md_alloc(alloc, num_nodes);
    memset(seen, 0, num_nodes);
    seen[t->root] = 1;
    for (size_t id = 1; id < num_nodes; ++id) {
        const vial_node_t* n = &t->nodes[id];
        if (n->kind == VIAL_N_NONE || n->kind >= VIAL_N_COUNT) return "a node of no kind";
        if (n->span.beg > n->span.end || n->span.end > len) return "a node's span is out of bounds";
        if ((size_t)n->child_beg + n->child_count > num_children) return "a node's children are out of bounds";
        uint32_t prev_end = n->span.beg;
        for (uint32_t i = 0; i < n->child_count; ++i) {
            const vial_node_id c = t->children[n->child_beg + i];
            if (c == 0) continue;
            if (c >= num_nodes) return "a child id is not a node";
            if (c == id) return "a node is its own child";
            if (seen[c]) return "a node has two parents";
            seen[c] = 1;
            const md_vial_span_t s = t->nodes[c].span;
            if (s.beg < n->span.beg || s.end > n->span.end) return "a child lies outside its parent";
            if (s.beg < prev_end) return "children are out of source order";
            prev_end = s.end;
        }
    }
    for (size_t id = 1; id < num_nodes; ++id) {
        if (!seen[id]) return "a node is in no tree";
    }
    return NULL;
}

static inline const char* vial_fuzz__check_diags(const md_vial_diagnostic_t* d, size_t n, uint32_t len, bool checked) {
    for (size_t i = 0; i < n; ++i) {
        if (d[i].span.beg > d[i].span.end || d[i].span.end > len) return "a diagnostic's span is out of bounds";
        if (i > 0 && d[i - 1].span.beg > d[i].span.beg) return "diagnostics are out of source order";
        if (d[i].code[5] != '\0' || !(vial_fuzz__front_end_code(d[i].code) || (checked && vial_fuzz__checker_code(d[i].code)))) return "a diagnostic with a code the compiler does not give";
        if (str_empty(d[i].message)) return "a diagnostic without a message";
        if (memchr(d[i].message.ptr, '\n', d[i].message.len) || memchr(d[i].message.ptr, '\r', d[i].message.len)) return "a message is more than one line";
    }
    if (n > len + 1) return "more diagnostics than bytes";
    return NULL;
}

// src need not end in a NUL; the stages must not read past src.len. alloc should be an arena: nothing is freed.
static inline bool vial_fuzz__frame(void* user, int64_t index, vial_frame_t* out) {
    const md_system_state_t* st = user;
    if (index != 0 || !st) return false;
    *out = (vial_frame_t){.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell, .index = 0};
    return true;
}

// One evaluation over one frame, its interpreter's memory from eval_alloc
static inline bool vial_fuzz__run(vial_interp_t* in, const vial_ir_t* ir, vial_sys_t* vs, const struct md_system_t* sys,
                                  const md_system_state_t* state, md_allocator_i* eval_alloc) {
    if (!vial_interp_init(in, ir, vs, sys, eval_alloc)) return false;
    in->max_elements = (size_t)1 << 20;
    in->load_frame = vial_fuzz__frame;
    in->load_user = (void*)state;
    vial_interp_run(in, 1);
    return true;
}

// The program over one frame: NULL, or which promise broke. It is evaluated twice, the second time in an arena whose
// memory is dirty and at other addresses, and must give the same bits.
static inline const char* vial_fuzz__evaluate(const md_vial_program_t* prog, const struct md_system_t* sys, const md_system_state_t* state,
                                              md_allocator_i* alloc) {
    const vial_ir_t* ir = md_vial__run_ir(prog);
    if (!ir || !sys || !state) return NULL;
    vial_sys_t vs;
    if (!vial_sys_init(&vs, sys, alloc)) return NULL;
    vial_interp_t in, again;
    if (!vial_fuzz__run(&in, ir, &vs, sys, state, md_get_heap_allocator())) return NULL;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
    void* junk = md_alloc(arena, KILOBYTES(32));
    memset(junk, 0xA5, KILOBYTES(32));
    md_free(arena, junk, KILOBYTES(32));
    const bool ran = vial_fuzz__run(&again, ir, &vs, sys, state, arena);
    static char buf[200];
    const char* why = NULL;
    for (size_t k = 0; k < md_array_size(ir->outputs) && !why; ++k) {
        const vial_val_t* v = vial_interp_at(&in, ir->outputs[k].value, 0);
        const vial_type_t* t = &ir->values[ir->outputs[k].value].type;
        const int nl = (int)(ir->outputs[k].name.len < 40 ? ir->outputs[k].name.len : 40);
        if (ran) {
            const vial_val_t* w = vial_interp_at(&again, ir->outputs[k].value, 0);
            if ((!v) != (!w) || (v && !vial_val_equal(v, w))) {
                snprintf(buf, sizeof(buf), "evaluating twice gives different values: %.*s", nl, ir->outputs[k].name.ptr);
                why = buf;
                break;
            }
        }
        if (!v || t->dtype == VIAL_DT_TUPLE || t->dtype == VIAL_DT_ERROR || t->dtype == VIAL_DT_DEF) continue;
        const uint8_t nd = (v->flags & VIAL_VAL_FRAMES) ? v->ndim - 1 : v->ndim;
        if (v->dtype != t->dtype || nd != t->ndim) {
            snprintf(buf, sizeof(buf), "an evaluated value differs from its type: %.*s is dtype %d rank %u, typed %d rank %u",
                nl, ir->outputs[k].name.ptr, v->dtype, nd, t->dtype, t->ndim);
            why = buf;
        }
    }
    if (ran) vial_interp_free(&again);
    vial_interp_free(&in);
    md_arena_allocator_destroy(arena);
    return why;
}

static inline const char* vial_fuzz_check(str_t src, const struct md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    if (src.len > UINT32_MAX / 2) return NULL;
    const uint32_t len = (uint32_t)src.len;

    vial_diags_t diags;
    vial_diags_init(&diags, src, alloc);
    vial_token_t* tokens = vial_lex(src, &diags, alloc);
    const size_t num_tokens = md_array_size(tokens);
    if (num_tokens == 0 || tokens[num_tokens - 1].kind != MD_VIAL_TOKEN_EOF) return "the tokens do not end in EOF";

    // Tokens: in bounds, in order, and nothing but layout between them
    uint32_t covered = 0;
    for (size_t i = 0; i < num_tokens; ++i) {
        const md_vial_span_t s = tokens[i].span;
        if (s.beg > s.end || s.end > len) return "a token's span is out of bounds";
        if (i > 0 && tokens[i - 1].span.beg > s.beg) return "tokens are out of order";
        if (s.beg == s.end) continue;   // INDENT, DEDENT, EOF and the NEWLINE that ends the text
        if (s.beg < covered) return "two tokens overlap";
        for (uint32_t j = covered; j < s.beg; ++j) {
            if (!vial_fuzz__layout_byte(src.ptr[j])) return "a byte that is not layout belongs to no token";
        }
        covered = s.end;
    }
    for (uint32_t j = covered; j < len; ++j) {
        if (!vial_fuzz__layout_byte(src.ptr[j])) return "a byte that is not layout belongs to no token";
    }

    const size_t num_lex_diags = md_array_size(diags.items);
    vial_tree_t tree = vial_parse(tokens, src, &diags, alloc);
    const char* err = vial_fuzz__check_tree(&tree, len, num_tokens, alloc);
    if (err) return err;
    const size_t num_diags = md_array_size(diags.items);
    if (num_diags < num_lex_diags) return "the parser dropped a diagnostic of the lexer";
    err = vial_fuzz__check_diags(diags.items, num_diags, len, false);
    if (err) return err;

    str_t dump = vial_tree_dump(&tree, src, alloc);
    if (str_empty(dump)) return "the tree dump is empty";

    // The public API, twice without a system and, given one, twice with it: the same tokens and diagnostics each time,
    // so compiling is deterministic
    size_t first_count = 0;
    const int rounds = sys ? 4 : 2;
    for (int round = 0; round < rounds; ++round) {
        const md_vial_compile_desc_t desc = {.sys = round >= 2 ? sys : NULL};
        md_vial_program_t* prog = md_vial_compile(src, &desc, alloc);
        if (!prog) return "md_vial_compile returned NULL";
        const md_vial_token_t* pt = NULL;
        const md_vial_diagnostic_t* pd = NULL;
        const size_t npt = md_vial_tokens(prog, &pt);
        const size_t npd = md_vial_diagnostics(prog, &pd);
        const char* why = NULL;
        if (npt != num_tokens) why = "md_vial_tokens differs from the lexer";
        for (size_t i = 0; !why && i < npt; ++i) {
            if (pt[i].kind != tokens[i].kind || pt[i].span.beg != tokens[i].span.beg || pt[i].span.end != tokens[i].span.end) {
                why = "md_vial_tokens differs from the lexer";
            }
        }
        // The front end's diagnostics are all there, as the stages gave them, among the checker's; and the whole list
        // keeps the promises, in source order. The two rounds give the same list.
        if (!why) why = vial_fuzz__check_diags(pd, npd, len, true);
        size_t k = 0;
        for (size_t i = 0; !why && i < npd && k < num_diags; ++i) {
            if (memcmp(pd[i].code, diags.items[k].code, 6) == 0 && pd[i].span.beg == diags.items[k].span.beg &&
                pd[i].span.end == diags.items[k].span.end && str_eq(pd[i].message, diags.items[k].message) &&
                str_eq(pd[i].fix, diags.items[k].fix)) {
                k += 1;
            }
        }
        if (!why && k != num_diags) why = "md_vial_diagnostics lost a diagnostic of the front end";
        if (!why && (round & 1) == 0) first_count = npd;
        if (!why && (round & 1) == 1 && npd != first_count) why = "compiling twice gives different diagnostics";
        if (!why && !vial_ir_verify(md_vial__ir(prog), NULL, 0)) why = "the array IR does not verify";
        if (!why && md_vial__run_ir(prog) && !vial_ir_verify(md_vial__run_ir(prog), NULL, 0)) why = "the inlined IR does not verify";
        if (!why && round == 3) why = vial_fuzz__evaluate(prog, sys, state, alloc);
        md_vial_program_free(prog);
        if (why) return why;
    }
    return NULL;
}

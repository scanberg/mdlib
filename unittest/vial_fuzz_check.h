#pragma once

// What the VIAL front end promises for any input at all, as one check. The built-in fuzzer (test_vial_fuzz.c) and
// the libFuzzer target (fuzz_vial.c) both run it; it returns NULL when every promise holds and otherwise what broke.
//
// Tokens       spans in bounds and in order, the last one EOF; every byte no token covers is layout
// Tree         node 0 NONE, the root a MODULE over the whole source, every child id valid and used once, children
//              inside their parent and in source order, the node count linear in the token count
// Diagnostics  spans in bounds, in source order, a known code, a message
// The API      md_vial_compile gives the same tokens and diagnostics as the stages run by hand, twice over

#include <md_vial.h>
#include <vial/vial_internal.h>
#include <vial/ast.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

#include <string.h>

// The codes the front end may give
static inline bool vial_fuzz__front_end_code(const char* code) {
    static const char* codes[] = {"E0101", "E0102", "E0103", "E0104", "E0502", "E0701", "E0801", "E0803"};
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

static inline const char* vial_fuzz__check_diags(const md_vial_diagnostic_t* d, size_t n, uint32_t len) {
    for (size_t i = 0; i < n; ++i) {
        if (d[i].span.beg > d[i].span.end || d[i].span.end > len) return "a diagnostic's span is out of bounds";
        if (i > 0 && d[i - 1].span.beg > d[i].span.beg) return "diagnostics are out of source order";
        if (d[i].code[5] != '\0' || !vial_fuzz__front_end_code(d[i].code)) return "a diagnostic with a code the front end does not give";
        if (str_empty(d[i].message)) return "a diagnostic without a message";
        if (d[i].message.ptr[d[i].message.len - 1] == '\n') return "a message ends in a newline";
    }
    if (n > len + 1) return "more diagnostics than bytes";
    return NULL;
}

// src need not end in a NUL; the stages must not read past src.len. alloc should be an arena: nothing is freed.
static inline const char* vial_fuzz_check(str_t src, md_allocator_i* alloc) {
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
    err = vial_fuzz__check_diags(diags.items, num_diags, len);
    if (err) return err;

    str_t dump = vial_tree_dump(&tree, src, alloc);
    if (str_empty(dump)) return "the tree dump is empty";

    // The public API, twice: the same tokens and diagnostics, so compiling is deterministic
    for (int round = 0; round < 2; ++round) {
        md_vial_program_t* prog = md_vial_compile(src, NULL, alloc);
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
        if (!why && npd != num_diags) why = "md_vial_diagnostics differs from the stages";
        for (size_t i = 0; !why && i < npd; ++i) {
            if (memcmp(pd[i].code, diags.items[i].code, 6) != 0 || pd[i].span.beg != diags.items[i].span.beg ||
                pd[i].span.end != diags.items[i].span.end || !str_eq(pd[i].message, diags.items[i].message) ||
                !str_eq(pd[i].fix, diags.items[i].fix)) {
                why = "md_vial_diagnostics differs from the stages";
            }
        }
        md_vial_program_free(prog);
        if (why) return why;
    }
    return NULL;
}

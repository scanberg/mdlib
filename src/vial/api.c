#include "check.h"
#include "infer.h"

#include <core/md_arena_allocator.h>
#include <core/md_common.h>

#include <string.h>
#include <stdio.h>

struct md_vial_program_t {
    md_allocator_i*  arena;
    str_t            src;
    vial_token_t*    tokens;
    md_vial_token_t* pub_tokens;
    vial_diags_t     diags;
    vial_tree_t      tree;
    vial_check_t     check;
    vial_types_t     types;
    // The prelude, compiled with every script
    str_t            prelude_src;
    vial_diags_t     prelude_diags;
    vial_tree_t      prelude;
};

#include <vial_prelude.inl>   // prelude_vial[], prelude_vial_size: src/vial/prelude.vial, embedded when configuring

md_vial_program_t* md_vial_compile(str_t src, const md_vial_compile_desc_t* desc, md_allocator_i* alloc) {
    (void)desc;
    md_allocator_i* arena = md_arena_allocator_create(alloc, KILOBYTES(64));
    if (!arena) return NULL;
    md_vial_program_t* prog = md_alloc(arena, sizeof(md_vial_program_t));
    memset(prog, 0, sizeof(*prog));
    prog->arena = arena;
    prog->src = str_copy(src, arena);

    vial_diags_init(&prog->diags, prog->src, arena);
    prog->tokens = vial_lex(prog->src, &prog->diags, arena);
    prog->tree = vial_parse(prog->tokens, prog->src, &prog->diags, arena);

    const bool with_prelude = !(desc && (desc->flags & MD_VIAL_COMPILE_NO_PRELUDE));
    if (with_prelude) {
        prog->prelude_src = (str_t){(const char*)prelude_vial, prelude_vial_size - 1};
        vial_diags_init(&prog->prelude_diags, prog->prelude_src, arena);
        vial_token_t* prelude_tokens = vial_lex(prog->prelude_src, &prog->prelude_diags, arena);
        prog->prelude = vial_parse(prelude_tokens, prog->prelude_src, &prog->prelude_diags, arena);
    }

    prog->check = (vial_check_t){.tree = &prog->tree, .src = prog->src, .diags = &prog->diags,
        .prelude = with_prelude ? &prog->prelude : NULL, .prelude_src = prog->prelude_src,
        .prelude_diags = &prog->prelude_diags, .sys = desc ? desc->sys : NULL, .alloc = arena};
    vial_check_names(&prog->check);
    vial_check_types(&prog->check, &prog->types);

    const size_t n = md_array_size(prog->tokens);
    md_array_resize(prog->pub_tokens, n, arena);
    for (size_t i = 0; i < n; ++i) {
        prog->pub_tokens[i] = (md_vial_token_t){.kind = prog->tokens[i].kind, .span = prog->tokens[i].span};
    }
    return prog;
}

void md_vial_program_free(md_vial_program_t* prog) {
    if (prog) md_arena_allocator_destroy(prog->arena);
}

size_t md_vial_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->diags.items;
    return md_array_size(prog->diags.items);
}

size_t md_vial_tokens(const md_vial_program_t* prog, const md_vial_token_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->pub_tokens;
    return md_array_size(prog->pub_tokens);
}

size_t md_vial__prelude_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->prelude_diags.items;
    return md_array_size(prog->prelude_diags.items);
}

// The names a statement binds, flattened in order
static void target_names(const vial_tree_t* t, vial_node_id target, vial_node_id** out, md_allocator_i* alloc) {
    const vial_node_kind_t k = (vial_node_kind_t)t->nodes[target].kind;
    if (k == VIAL_N_NAME) { md_array_push(*out, target, alloc); return; }
    if (k == VIAL_N_TUPLE || k == VIAL_N_LIST) {
        for (uint32_t i = 0; i < t->nodes[target].child_count; ++i) target_names(t, vial_child(t, target, i), out, alloc);
    }
}

bool md_vial__line_result(const md_vial_program_t* prog, uint32_t line, md_vial__result_t* r) {
    memset(r, 0, sizeof(*r));
    const vial_tree_t* t = &prog->tree;
    const vial_node_id root = t->root;
    for (uint32_t i = 0; i < t->nodes[root].child_count; ++i) {
        const vial_node_id s = vial_child(t, root, i);
        if (vial_line_of(&prog->diags, t->nodes[s].span.beg) != line) continue;
        vial_node_id* nodes = NULL;
        if (t->nodes[s].kind == VIAL_N_ASSIGN) {
            const uint32_t n = t->nodes[s].child_count;
            if (n >= 2) target_names(t, vial_child(t, s, 0), &nodes, prog->arena);
        } else if (t->nodes[s].kind == VIAL_N_EXPR_STMT) {
            md_array_push(nodes, vial_child(t, s, 0), prog->arena);
        }
        if (!md_array_size(nodes)) return false;
        r->num_types = 0;
        r->typed = true;
        for (size_t k = 0; k < md_array_size(nodes) && k < 4; ++k) {
            const vial_type_t* ty = &prog->types.types[nodes[k]];
            if (vial_type_is_error(ty)) r->typed = false;
            vial_type_print(r->types[k], sizeof(r->types[k]), ty);
            r->num_types += 1;
        }
        const vial_type_t* first = &prog->types.types[nodes[0]];
        snprintf(r->stage, sizeof(r->stage), "%s", vial_stage_str(first));
        r->pass = first->pass;
        const vial_value_t* v = &prog->types.values[nodes[0]];
        if (vial_value_is_scalar(v)) { r->has_value = true; vial_value_print(r->value, sizeof(r->value), v); }
        return true;
    }
    return false;
}

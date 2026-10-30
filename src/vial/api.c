#include "check.h"

#include <md_system.h>
#include <md_attributes.h>
#include "infer.h"
#include "ir.h"

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
    vial_ir_t        ir;        // the script lowered
    vial_ir_t        run_ir;    // and inlined, for evaluation; empty when it does not inline
    bool             has_run_ir;
    md_vial_output_t* outputs;  // md_array
    // The prelude, compiled with every script
    str_t            prelude_src;
    vial_diags_t     prelude_diags;
    vial_tree_t      prelude;
};

#include <vial_prelude.inl>   // prelude_vial[], prelude_vial_size: src/vial/prelude.vial, embedded when configuring

static md_vial_stage_t public_stage(const vial_type_t* t) {
    switch (t->stage) {
    case VIAL_STAGE_FRAME: return MD_VIAL_STAGE_FRAME;
    case VIAL_STAGE_TRAJECTORY: return (t->flags & VIAL_TYPE_PER_FRAME) ? MD_VIAL_STAGE_TRAJECTORY_PER_FRAME : MD_VIAL_STAGE_TRAJECTORY;
    default: return MD_VIAL_STAGE_STATIC;
    }
}

static uint32_t public_pass(const vial_type_t* t) { return t->stage != VIAL_STAGE_STATIC && t->pass == 0 ? 1 : t->pass; }

static void target_names(const vial_tree_t* t, vial_node_id target, vial_node_id** out, md_allocator_i* alloc);

// The published values: top-level names, bound once each, to what can be published
static void collect_outputs(md_vial_program_t* prog) {
    const vial_tree_t* t = &prog->tree;
    const vial_node_id root = t->root;
    for (uint32_t i = 0; i < t->nodes[root].child_count; ++i) {
        const vial_node_id s = vial_child(t, root, i);
        if (t->nodes[s].kind != VIAL_N_ASSIGN || t->nodes[s].child_count < 2) continue;
        vial_node_id* names = NULL;
        for (uint32_t k = 0; k + 1 < t->nodes[s].child_count; ++k) target_names(t, vial_child(t, s, k), &names, prog->arena);
        for (size_t k = 0; k < md_array_size(names); ++k) {
            const vial_node_id nn = names[k];
            const md_vial_span_t sp = t->nodes[nn].span;
            const str_t name = {prog->src.ptr + sp.beg, sp.end - sp.beg};
            if (!name.len || name.ptr[0] == '_') continue;
            const vial_type_t* ty = &prog->types.types[nn];
            const bool number = ty->dtype == VIAL_DT_BOOL || ty->dtype == VIAL_DT_INT || ty->dtype == VIAL_DT_FLOAT || ty->dtype == VIAL_DT_FLOAT64;
            const bool set = ty->dtype == VIAL_DT_SET;
            if (!(number && ty->ndim <= 2) && !set) continue;
            md_vial_output_t o = {.name = name, .span = sp, .stage = public_stage(ty), .pass = public_pass(ty), .selection = set};
            vial_type_print(o.type, sizeof(o.type), ty);
            md_array_push(prog->outputs, o, prog->arena);
        }
    }
}

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
        .prelude_diags = &prog->prelude_diags, .sys = desc ? desc->sys : NULL,
        .num_frames = desc ? desc->num_frames : 0, .run = desc ? str_copy(desc->run, arena) : (str_t){0}, .alloc = arena};
    if (!prog->check.num_frames && prog->check.sys && !str_empty(prog->check.run)) {
        // The run's frame axis
        char path[256];
        snprintf(path, sizeof(path), "%.*s/time", (int)prog->check.run.len, prog->check.run.ptr);
        const md_attribute_t* time = md_attributes_find(&prog->check.sys->attributes, str_from_cstr(path));
        if (time && time->format.rank >= 1) prog->check.num_frames = time->format.shape[0];
    }
    vial_ir_init(&prog->ir, arena);
    prog->check.ir = &prog->ir;
    vial_check_names(&prog->check);
    vial_check_types(&prog->check, &prog->types);
    if (!md_array_size(prog->ir.failures)) prog->has_run_ir = vial_ir_inline(&prog->ir, &prog->run_ir, arena);

    collect_outputs(prog);

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
        r->pass = first->stage != VIAL_STAGE_STATIC && first->pass == 0 ? 1 : first->pass;    // frame values are pass 1
        const vial_value_t* v = &prog->types.values[nodes[0]];
        if (vial_value_is_scalar(v)) { r->has_value = true; vial_value_print(r->value, sizeof(r->value), v); }
        return true;
    }
    return false;
}

size_t md_vial_outputs(const md_vial_program_t* prog, const md_vial_output_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->outputs;
    return md_array_size(prog->outputs);
}

bool md_vial_hover(const md_vial_program_t* prog, uint32_t offset, md_vial_hover_t* out) {
    memset(out, 0, sizeof(*out));
    if (!prog) return false;
    // The smallest expression with a type whose span holds offset. Expressions inside a def's body are typed per call
    // and have no one type: hover shows them at the calls.
    const vial_tree_t* t = &prog->tree;
    vial_node_id best = 0;
    uint32_t best_len = UINT32_MAX;
    for (vial_node_id i = 1; i < md_array_size(t->nodes); ++i) {
        const vial_node_t* n = &t->nodes[i];
        if (n->kind < VIAL_N_NAME || n->kind == VIAL_N_CMPOP || n->kind == VIAL_N_KWARG || n->kind == VIAL_N_SLICE || n->kind == VIAL_N_UNIT) continue;
        if (offset < n->span.beg || offset >= n->span.end) continue;
        if (!prog->types.typed[i] || vial_type_is_error(&prog->types.types[i])) continue;
        const uint32_t len = n->span.end - n->span.beg;
        if (len < best_len) { best = i; best_len = len; }
    }
    if (!best) return false;
    const vial_type_t* ty = &prog->types.types[best];
    out->span = t->nodes[best].span;
    vial_type_print(out->type, sizeof(out->type), ty);
    out->stage = public_stage(ty);
    out->pass = public_pass(ty);
    const vial_value_t* v = &prog->types.values[best];
    if (vial_value_is_scalar(v)) vial_value_print(out->value, sizeof(out->value), v);
    return true;
}

const struct vial_ir_t* md_vial__ir(const md_vial_program_t* prog) { return prog ? &prog->ir : NULL; }
const struct vial_ir_t* md_vial__run_ir(const md_vial_program_t* prog) { return prog && prog->has_run_ir ? &prog->run_ir : NULL; }

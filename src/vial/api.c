#include "check.h"

#include <md_system.h>
#include <md_attributes.h>
#include "infer.h"
#include "ir.h"
#include "builtins.h"
#include "infer_internal.h"
#include "select.h"
#include "volume.h"

#include <md_util.h>
#include <core/md_array.h>
#include <core/md_intrinsics.h>
#include <core/md_os.h>
#include <stdlib.h>

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
    uint32_t*        run_of;    // by value of ir: its value in run_ir, for the top level
    md_vial_output_t* outputs;  // md_array
    md_vial_mark_t*  marks;     // md_array, made when compiling
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
static void make_marks(md_vial_program_t* prog);

// A name's value in the inlined IR, 0 when it has none
static vial_ssa_t run_output(const md_vial_program_t* prog, str_t name) {
    if (!prog->has_run_ir) return 0;
    for (size_t i = 0; i < md_array_size(prog->run_ir.outputs); ++i) {
        if (str_eq(prog->run_ir.outputs[i].name, name)) return prog->run_ir.outputs[i].value;
    }
    return 0;
}

// Whether v is what histogram() gives, or an element of it
static bool from_histogram(const vial_ir_t* ir, vial_ssa_t v) {
    uint32_t op = ir->values[v].op;
    if (op != VIAL_IR_NO_OP && ir->ops[op].kind == VIAL_IR_ELEM && ir->ops[op].num_args) op = ir->values[ir->args[ir->ops[op].arg_beg]].op;
    return op != VIAL_IR_NO_OP && ir->ops[op].kind == VIAL_IR_BUILTIN && !strcmp(vial_builtin_at(ir->ops[op].aux)->name, "histogram");
}

// The volume op v is made by, sdf or density; NULL for another value
static const char* volume_of(const vial_ir_t* ir, vial_ssa_t v) {
    const uint32_t op = ir->values[v].op;
    if (op == VIAL_IR_NO_OP || ir->ops[op].kind != VIAL_IR_BUILTIN) return NULL;
    const char* name = vial_builtin_at(ir->ops[op].aux)->name;
    return vial_volume_builtin(name) ? name : NULL;
}

// A volume is added up over the frames at the top level: one inside a map, as a comprehension or a def mapped over
// groups gives, would be an array of volumes, which is not supported yet (E0406). The program is not evaluated.
static bool volumes_at_top(md_vial_program_t* prog) {
    const vial_ir_t* ir = &prog->run_ir;
    const uint32_t top = ir->funcs[0].body;
    bool ok = true;
    for (size_t i = 0; i < md_array_size(ir->ops); ++i) {
        const vial_ir_op_t* op = &ir->ops[i];
        if (op->kind != VIAL_IR_BUILTIN || op->parent == top || !vial_volume_builtin(vial_builtin_at(op->aux)->name)) continue;
        if (!vial_diag_at(&prog->diags, op->span.beg)) {
            vial_diag(&prog->diags, "E0406", op->span, "`%s` is added up over the frames at the top level; an array of volumes, from a comprehension or a def mapped over groups, is not supported yet",
                vial_builtin_at(op->aux)->name);
        }
        ok = false;
    }
    if (!ok) vial_diags_sort(&prog->diags);
    return ok;
}

// The published values: top-level names, bound once each, to what can be published
static void collect_outputs(md_vial_program_t* prog) {
    const vial_tree_t* t = &prog->tree;
    // Fingerprints share what they hash
    const size_t nv = prog->has_run_ir ? md_array_size(prog->run_ir.values) : 0;
    uint64_t* memo = nv ? md_alloc(prog->arena, nv * sizeof(uint64_t)) : NULL;
    if (memo) memset(memo, 0, nv * sizeof(uint64_t));
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
            if (!vial_types_published(&prog->types, name)) continue;
            const vial_type_t* ty = &prog->types.types[nn];
            const bool number = ty->dtype == VIAL_DT_BOOL || ty->dtype == VIAL_DT_INT || ty->dtype == VIAL_DT_FLOAT || ty->dtype == VIAL_DT_FLOAT64;
            const bool set = ty->dtype == VIAL_DT_SET;
            const vial_ssa_t rv = run_output(prog, name);
            const char* vol = rv ? volume_of(&prog->run_ir, rv) : NULL;
            if (!(number && (ty->ndim <= 2 || (vol && ty->ndim == 3))) && !set) continue;
            md_vial_output_t o = {.name = name, .span = sp, .stage = public_stage(ty), .pass = public_pass(ty), .selection = set};
            vial_type_print(o.type, sizeof(o.type), ty);
            if (o.stage == MD_VIAL_STAGE_FRAME || o.stage == MD_VIAL_STAGE_TRAJECTORY_PER_FRAME) o.flags |= MD_VIAL_OUTPUT_TEMPORAL;
            if (ty->period_end > ty->period_beg) o.flags |= MD_VIAL_OUTPUT_PERIODIC;
            if (vol) o.flags |= MD_VIAL_OUTPUT_VOLUME | (!strcmp(vol, "sdf") ? MD_VIAL_OUTPUT_SDF : 0);
            if (rv) {
                o.fingerprint = vial_ir_fingerprint(&prog->run_ir, rv, prog->check.sys, prog->check.run, memo);
                if (from_histogram(&prog->run_ir, rv)) o.flags |= MD_VIAL_OUTPUT_DISTRIBUTION;
            }
            md_array_push(prog->outputs, o, prog->arena);
        }
    }
}

// The prelude, lexed, parsed and resolved once for every program: it depends on no script and no system, and it was
// most of what compiling a short script, or a filter, cost. Read only once made; never freed.
static struct {
    volatile uint64_t   ticket;     // the first to take one builds it
    volatile uint64_t   ready;
    md_allocator_i*     arena;
    vial_diags_t        diags;
    vial_tree_t         tree;
    vial_names_t        names;
} shared_prelude;

static void prelude_build(void) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    const str_t src = {(const char*)prelude_vial, prelude_vial_size - 1};
    vial_diags_init(&shared_prelude.diags, src, arena);
    vial_token_t* tokens = vial_lex(src, &shared_prelude.diags, arena);
    shared_prelude.tree = vial_parse(tokens, src, &shared_prelude.diags, arena);
    vial_check_prelude(&shared_prelude.tree, src, &shared_prelude.diags, &shared_prelude.names, arena);
    shared_prelude.arena = arena;
}

static void prelude_ensure(void) {
    if (md_atomic_load64(&shared_prelude.ready)) return;
    if (md_atomic_fetch_add64(&shared_prelude.ticket, 1) == 0) {
        prelude_build();
        md_atomic_fetch_add64(&shared_prelude.ready, 1);
        return;
    }
    while (!md_atomic_load64(&shared_prelude.ready)) md_thread_sleep(0);
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
        prelude_ensure();
        prog->prelude_src = (str_t){(const char*)prelude_vial, prelude_vial_size - 1};
        prog->prelude_diags = shared_prelude.diags;     // read only: what the prelude itself says, normally nothing
        prog->prelude = shared_prelude.tree;
    }

    prog->check = (vial_check_t){.tree = &prog->tree, .src = prog->src, .diags = &prog->diags,
        .prelude = with_prelude ? &prog->prelude : NULL, .prelude_resolved = with_prelude ? &shared_prelude.names : NULL,
        .prelude_src = prog->prelude_src,
        .prelude_diags = &prog->prelude_diags, .sys = desc ? desc->sys : NULL,
        .num_frames = desc ? desc->num_frames : 0, .run = desc ? str_copy(desc->run, arena) : (str_t){0}, .alloc = arena};
    if (desc && desc->num_sets && desc->sets) {
        // The stored selections, copied: those with a name a script may bind
        md_vial_named_set_t* st = md_alloc(arena, desc->num_sets * sizeof(md_vial_named_set_t));
        size_t n = 0;
        for (size_t i = 0; i < desc->num_sets; ++i) {
            const md_vial_named_set_t* s = &desc->sets[i];
            if (!md_vial_name_valid(s->name)) continue;
            if (s->ndim > 2 || (s->ndim && !s->offsets)) continue;
            uint32_t* at = md_alloc(arena, (s->num_atoms ? s->num_atoms : 1) * sizeof(uint32_t));
            if (s->num_atoms) memcpy(at, s->atoms, s->num_atoms * sizeof(uint32_t));
            md_vial_named_set_t c = {.name = str_copy(s->name, arena), .atoms = at, .num_atoms = s->num_atoms, .ndim = s->ndim};
            if (s->ndim) {
                const size_t count = (size_t)s->dims[0] * (s->ndim > 1 ? s->dims[1] : 1);
                if (s->offsets[count] > s->num_atoms) continue;
                uint32_t* off = md_alloc(arena, (count + 1) * sizeof(uint32_t));
                memcpy(off, s->offsets, (count + 1) * sizeof(uint32_t));
                c.dims[0] = s->dims[0];
                c.dims[1] = s->ndim > 1 ? s->dims[1] : 0;
                c.offsets = off;
            }
            st[n++] = c;
        }
        prog->check.stored = st;
        prog->check.num_stored = n;
    }
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
    if (!md_array_size(prog->ir.failures)) prog->has_run_ir = vial_ir_inline_map(&prog->ir, &prog->run_ir, arena, &prog->run_of);
    if (prog->has_run_ir && !volumes_at_top(prog)) prog->has_run_ir = false;
    if (prog->has_run_ir) {
        vial_ir_pairs(&prog->run_ir, arena);
        vial_ir_contract(&prog->run_ir, arena);
    }

    collect_outputs(prog);
    make_marks(prog);

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

// The smallest expression with a type whose span holds offset, or 0. Expressions inside a def's body are typed per call
// and have no one type: hover shows them at the calls.
// An expression with a type: what hover shows, the editor marks and drawing takes
static bool is_expr(const md_vial_program_t* prog, vial_node_id i) {
    const vial_tree_t* t = &prog->tree;
    if (!i || i >= md_array_size(t->nodes)) return false;
    const vial_node_t* n = &t->nodes[i];
    if (n->kind < VIAL_N_NAME || n->kind == VIAL_N_CMPOP || n->kind == VIAL_N_KWARG || n->kind == VIAL_N_SLICE || n->kind == VIAL_N_UNIT) return false;
    return i < md_array_size(prog->types.typed) && prog->types.typed[i] && !vial_type_is_error(&prog->types.types[i]);
}

static vial_node_id hover_node(const md_vial_program_t* prog, uint32_t offset) {
    const vial_tree_t* t = &prog->tree;
    vial_node_id best = 0;
    uint32_t best_len = UINT32_MAX;
    for (vial_node_id i = 1; i < md_array_size(t->nodes); ++i) {
        const vial_node_t* n = &t->nodes[i];
        if (offset < n->span.beg || offset >= n->span.end) continue;
        if (!is_expr(prog, i)) continue;
        const uint32_t len = n->span.end - n->span.beg;
        if (len < best_len) { best = i; best_len = len; }
    }
    return best;
}

uint32_t md_vial_expr_at(const md_vial_program_t* prog, uint32_t offset) { return prog ? hover_node(prog, offset) : 0; }

bool md_vial_hover(const md_vial_program_t* prog, uint32_t offset, md_vial_hover_t* out) {
    memset(out, 0, sizeof(*out));
    return prog && md_vial_expr_info(prog, hover_node(prog, offset), out);
}

size_t md_vial_marks(const md_vial_program_t* prog, const md_vial_mark_t** out) {
    if (out) *out = prog ? prog->marks : NULL;
    return prog ? md_array_size(prog->marks) : 0;
}

// Every expression with a type, and how many of them hold it
static void make_marks(md_vial_program_t* prog) {
    const vial_tree_t* t = &prog->tree;
    const vial_names_t* nm = &prog->check.names;
    for (vial_node_id i = 1; i < md_array_size(t->nodes); ++i) {
        if (!is_expr(prog, i)) continue;
        uint32_t depth = 0;
        for (vial_node_id p = i < md_array_size(nm->parent) ? nm->parent[i] : 0; p; p = p < md_array_size(nm->parent) ? nm->parent[p] : 0) {
            if (is_expr(prog, p)) ++depth;
        }
        md_vial_mark_t m = {.expr = i, .span = t->nodes[i].span, .depth = depth};
        md_array_push(prog->marks, m, prog->arena);
    }
    // In source order: by start, the outer first
    const size_t n = md_array_size(prog->marks);
    for (size_t i = 1; i < n; ++i) {
        const md_vial_mark_t m = prog->marks[i];
        size_t j = i;
        while (j && (prog->marks[j - 1].span.beg > m.span.beg || (prog->marks[j - 1].span.beg == m.span.beg && prog->marks[j - 1].span.end < m.span.end))) {
            prog->marks[j] = prog->marks[j - 1];
            --j;
        }
        prog->marks[j] = m;
    }
}

uint32_t md_vial__expr_value(const md_vial_program_t* prog, uint32_t expr) {
    if (!prog || !is_expr(prog, expr) || expr - prog->types.base >= md_array_size(prog->types.ssa)) return 0;
    return prog->types.ssa[expr - prog->types.base];
}

bool md_vial_expr_info(const md_vial_program_t* prog, uint32_t expr, md_vial_hover_t* out) {
    memset(out, 0, sizeof(*out));
    if (!prog || !is_expr(prog, expr)) return false;
    const vial_node_id best = expr;
    const vial_type_t* ty = &prog->types.types[best];
    out->span = prog->tree.nodes[best].span;
    vial_type_print(out->type, sizeof(out->type), ty);
    out->stage = public_stage(ty);
    out->pass = public_pass(ty);
    const vial_value_t* v = &prog->types.values[best];
    if (vial_value_is_scalar(v)) vial_value_print(out->value, sizeof(out->value), v);
    if (ty->dtype != VIAL_DT_TUPLE && ty->unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(ty->unit)) {
        // The type prints its dtype and extents, then a space and the unit: neither of the first has a space
        char u[128];
        const size_t un = vial_unit_print(u, sizeof(u), ty->unit);
        const char* sp = strchr(out->type, ' ');
        if (sp && un && strncmp(sp + 1, u, un) == 0) {
            out->unit = ty->unit;
            out->unit_at = (uint32_t)(sp + 1 - out->type);
            out->unit_len = (uint32_t)un;
            if (v->kind == VIAL_VAL_FLOAT) { out->has_number = true; out->number = v->f; }
            else if (v->kind == VIAL_VAL_INT) { out->has_number = true; out->number = (double)v->i; }
        }
    }
    return true;
}

uint32_t md_vial__hover_value(const md_vial_program_t* prog, uint32_t offset, md_vial_span_t* span) {
    if (!prog) return 0;
    const vial_node_id n = hover_node(prog, offset);
    if (!n || n - prog->types.base >= md_array_size(prog->types.ssa)) return 0;
    if (span) *span = prog->tree.nodes[n].span;
    return prog->types.ssa[n - prog->types.base];
}

const struct vial_ir_t* md_vial__ir(const md_vial_program_t* prog) { return prog ? &prog->ir : NULL; }

uint32_t md_vial__run_value(const md_vial_program_t* prog, uint32_t v) {
    if (!prog || !prog->has_run_ir || !prog->run_of || v >= md_array_size(prog->ir.values)) return 0;
    return prog->run_of[v];
}
const struct vial_ir_t* md_vial__run_ir(const md_vial_program_t* prog) { return prog && prog->has_run_ir ? &prog->run_ir : NULL; }

// ### EDITOR ###

size_t md_vial_keywords(str_t* out, size_t cap) {
    size_t n = 0;
    for (int k = VIAL_KW_NONE + 1; k < VIAL_KW_COUNT; ++k, ++n) {
        if (out && n < cap) out[n] = str_from_cstr(vial_kw_str((vial_kw_t)k));
    }
    return n;
}

static bool ident_char(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'; }

// The prelude's top-level defs a script calls: those not starting with _, read from its source
static size_t prelude_defs(str_t* out, size_t cap) {
    const str_t src = {(const char*)prelude_vial, prelude_vial_size - 1};
    size_t n = 0;
    for (size_t i = 0; i + 4 < src.len; ++i) {
        if ((i && src.ptr[i - 1] != '\n') || strncmp(src.ptr + i, "def ", 4)) continue;
        size_t e = i + 4;
        while (e < src.len && ident_char(src.ptr[e])) ++e;
        if (e == i + 4 || src.ptr[i + 4] == '_') continue;
        if (out && n < cap) out[n] = (str_t){src.ptr + i + 4, e - i - 4};
        ++n;
    }
    return n;
}

static bool is_prelude_def(str_t name) {
    str_t defs[64];
    const size_t n = MIN(prelude_defs(defs, 64), 64);
    for (size_t i = 0; i < n; ++i) if (str_eq(defs[i], name)) return true;
    return false;
}

size_t md_vial_builtin_names(str_t* out, size_t cap) {
    size_t n = 0;
    for (size_t i = 0; i < vial_builtin_count(); ++i) {
        const char* nm = vial_builtin_at(i)->name;
        if (nm[0] == '_') continue;
        if (out && n < cap) out[n] = str_from_cstr(nm);
        ++n;
    }
    str_t defs[64];
    const size_t nd = MIN(prelude_defs(defs, 64), 64);
    for (size_t i = 0; i < nd; ++i) {
        if (vial_builtin_find(defs[i])) continue;
        if (out && n < cap) out[n] = defs[i];
        ++n;
    }
    return n;
}

static bool is_keyword(str_t name) {
    for (int k = VIAL_KW_NONE + 1; k < VIAL_KW_COUNT; ++k) if (str_eq(name, str_from_cstr(vial_kw_str((vial_kw_t)k)))) return true;
    return false;
}

bool md_vial_name_valid(str_t name) {
    if (!name.len || (name.ptr[0] >= '0' && name.ptr[0] <= '9')) return false;
    for (size_t i = 0; i < name.len; ++i) if (!ident_char(name.ptr[i])) return false;
    return !is_keyword(name) && !vial_builtin_find(name) && !is_prelude_def(name);
}

str_t md_vial_attr_path(str_t run, str_t name) {
    if (!run.len || name.len <= run.len + 1 || strncmp(name.ptr, run.ptr, run.len) || name.ptr[run.len] != '/') return (str_t){0};
    return (str_t){name.ptr + run.len + 1, name.len - run.len - 1};
}

// Completion

typedef struct cmpl_t {
    md_vial_completion_t*   items;  // md_array
    str_t                   prefix;
    md_allocator_i*         alloc;
} cmpl_t;

static void cmpl_add(cmpl_t* c, str_t label, str_t text, md_vial_completion_kind_t kind) {
    if (!label.len || label.len < c->prefix.len || (c->prefix.len && strncmp(label.ptr, c->prefix.ptr, c->prefix.len))) return;
    for (size_t i = 0; i < md_array_size(c->items); ++i) if (str_eq(c->items[i].label, label)) return;
    md_vial_completion_t it = {.label = str_copy(label, c->alloc), .text = str_copy(text, c->alloc), .kind = kind};
    md_array_push(c->items, it, c->alloc);
}

static bool tok_op(const vial_token_t* t, vial_op_t op) {
    return (t->kind == MD_VIAL_TOKEN_OPERATOR || t->kind == MD_VIAL_TOKEN_PUNCT) && t->sub == op;
}

static str_t tok_text(const md_vial_program_t* prog, const vial_token_t* t) {
    return (str_t){prog->src.ptr + t->span.beg, t->span.end - t->span.beg};
}

// The innermost call whose parentheses hold offset, and the name it calls; 0 when none
static vial_node_id call_at(const md_vial_program_t* prog, uint32_t offset, str_t callee) {
    const vial_tree_t* t = &prog->tree;
    vial_node_id best = 0;
    uint32_t best_len = UINT32_MAX;
    for (vial_node_id i = 1; i < md_array_size(t->nodes); ++i) {
        const vial_node_t* n = &t->nodes[i];
        if (n->kind != VIAL_N_CALL || offset <= n->span.beg || offset > n->span.end) continue;
        const vial_node_id f = vial_child(t, i, 0);
        if (!f) continue;
        // The callee's name: f, or the name after the dot of x.f
        vial_node_id nm = f;
        if (t->nodes[f].kind == VIAL_N_ATTR) nm = vial_child(t, f, 1);
        if (!nm || t->nodes[nm].kind != VIAL_N_NAME) continue;
        const md_vial_span_t sp = t->nodes[nm].span;
        if (!str_eq((str_t){prog->src.ptr + sp.beg, sp.end - sp.beg}, callee) || offset <= t->nodes[f].span.end) continue;
        const uint32_t len = n->span.end - n->span.beg;
        if (len < best_len) { best = i; best_len = len; }
    }
    return best;
}

// A parameter as name=, unless the call gives it already
static void param_add(cmpl_t* c, str_t name, const str_t* given, size_t num_given) {
    if (!name.len || name.ptr[0] == '_') return;
    for (size_t g = 0; g < num_given; ++g) if (str_eq(given[g], name)) return;
    char buf[128];
    const int len = snprintf(buf, sizeof(buf), "%.*s=", (int)name.len, name.ptr);
    if (len > 0 && len < (int)sizeof(buf)) cmpl_add(c, name, (str_t){buf, (size_t)len}, MD_VIAL_COMPLETION_PARAMETER);
}

// The parameters of the call around the cursor: a built-in's as the checker bound it, a def's from its PARAM nodes
static void call_params(cmpl_t* c, const md_vial_program_t* prog, vial_node_id call) {
    const vial_tree_t* t = &prog->tree;
    str_t given[32];
    size_t ng = 0;
    for (uint32_t i = 1; i < t->nodes[call].child_count && ng < 32; ++i) {
        const vial_node_id a = vial_child(t, call, i);
        if (!a || t->nodes[a].kind != VIAL_N_KWARG) continue;
        const vial_node_id kn = vial_child(t, a, 0);
        if (kn) given[ng++] = (str_t){prog->src.ptr + t->nodes[kn].span.beg, t->nodes[kn].span.end - t->nodes[kn].span.beg};
    }
    if (call < md_array_size(prog->types.sigs) && prog->types.sigs[call].params) {
        const vx_sig_t* sg = &prog->types.sigs[call];
        for (uint32_t k = 0; k < sg->np; ++k) param_add(c, str_from_cstr(sg->params[k].name), given, ng);
        return;
    }
    vial_node_id f = vial_child(t, call, 0);
    if (f && t->nodes[f].kind == VIAL_N_ATTR) f = vial_child(t, f, 1);
    // A selector: the keywords of its level
    if (f && t->nodes[f].kind == VIAL_N_NAME) {
        const str_t fname = {prog->src.ptr + t->nodes[f].span.beg, t->nodes[f].span.end - t->nodes[f].span.beg};
        for (uint8_t lv = VIAL_LEVEL_PARTICLE; lv <= VIAL_LEVEL_RING; ++lv) {
            const bool alias = lv == VIAL_LEVEL_PARTICLE && (str_eq(fname, STR_LIT("atom")) || str_eq(fname, STR_LIT("bead")));
            if (!alias && !str_eq(fname, str_from_cstr(vial_select_name(lv)))) continue;
            static const char* words[] = {"name", "element", "kind", "resid", "resname", "chain", "id", "size", "scope"};  // not serial, until it is supported
            for (size_t w = 0; w < ARRAY_SIZE(words); ++w) {
                if (alias && lv == VIAL_LEVEL_PARTICLE && str_eq(str_from_cstr(words[w]), STR_LIT("kind"))) continue;    // atom() and bead() are a kind
                if (vial_select_takes(lv, vial_select_key(str_from_cstr(words[w])))) param_add(c, str_from_cstr(words[w]), given, ng);
            }
            return;
        }
    }
    // A def: the script's, or the prelude's
    const vial_symbol_t* sym = f ? vial_names_symbol(&prog->check.names, f) : NULL;
    if (!sym || (sym->kind != VIAL_SYM_DEF && sym->kind != VIAL_SYM_PRELUDE)) return;
    const vial_tree_t* dt = sym->in_prelude || sym->kind == VIAL_SYM_PRELUDE ? &prog->prelude : t;
    const str_t dsrc = dt == t ? prog->src : prog->prelude_src;
    // The symbol's node is the def's, or its name's
    const vial_names_t* dn = dt == t ? &prog->check.names : &prog->check.prelude_names;
    vial_node_id def = sym->node;
    if (def && def < md_array_size(dt->nodes) && dt->nodes[def].kind != VIAL_N_DEF && def < md_array_size(dn->parent)) def = dn->parent[def];
    if (!def || def >= md_array_size(dt->nodes) || dt->nodes[def].kind != VIAL_N_DEF) return;
    const vial_node_id ps = vial_child(dt, def, 1);
    for (uint32_t i = 0; ps && i < dt->nodes[ps].child_count; ++i) {
        const vial_node_id pn = vial_child(dt, ps, i);
        const vial_node_id nm = pn ? vial_child(dt, pn, 0) : 0;
        if (nm) param_add(c, (str_t){dsrc.ptr + dt->nodes[nm].span.beg, dt->nodes[nm].span.end - dt->nodes[nm].span.beg}, given, ng);
    }
}

// Whether the scope a symbol is bound in holds offset: the module always, a def, a comprehension or a loop where its
// node's span does
static bool scope_holds(const md_vial_program_t* prog, uint32_t scope, uint32_t offset) {
    const vial_names_t* nm = &prog->check.names;
    while (scope) {
        const vial_node_id n = nm->scopes[scope].node;
        if (n && (offset < prog->tree.nodes[n].span.beg || offset > prog->tree.nodes[n].span.end)) return false;
        const uint32_t up = nm->scopes[scope].parent;
        if (up == scope) break;
        scope = up;
    }
    return true;
}

// The values a selector keyword matches, in the system
static void system_values(cmpl_t* c, const md_system_t* sys, str_t callee, str_t kw) {
    #define IS(a, b) (str_eq(a, STR_LIT(b)))
    const bool particles = IS(callee, "particle") || IS(callee, "atom") || IS(callee, "bead");
    const bool residues = IS(callee, "residue");
    if (particles && IS(kw, "kind")) {
        static const char* k[] = {"atom", "bead", "virtual"};
        for (size_t i = 0; i < ARRAY_SIZE(k); ++i) cmpl_add(c, str_from_cstr(k[i]), str_from_cstr(k[i]), MD_VIAL_COMPLETION_VALUE);
        return;
    }
    if (residues && IS(kw, "kind")) {
        static const char* k[] = {"amino_acid", "nucleotide", "water", "ion", "lipid", "other"};
        for (size_t i = 0; i < ARRAY_SIZE(k); ++i) cmpl_add(c, str_from_cstr(k[i]), str_from_cstr(k[i]), MD_VIAL_COMPLETION_VALUE);
        return;
    }
    if (IS(callee, "part") && IS(kw, "kind")) {
        static const char* k[] = {"polymer", "peptide", "nucleic", "dna", "rna", "glycan", "ligand", "lipid", "ion", "water"};
        for (size_t i = 0; i < ARRAY_SIZE(k); ++i) cmpl_add(c, str_from_cstr(k[i]), str_from_cstr(k[i]), MD_VIAL_COMPLETION_VALUE);
        return;
    }
    if (!sys) return;
    if (particles && (IS(kw, "name") || !kw.len)) {
        for (size_t t = 0; t < sys->atom.type.count; ++t) {
            const str_t n = md_particle_type_name(&sys->atom.type, t);
            cmpl_add(c, n, n, MD_VIAL_COMPLETION_VALUE);
        }
    } else if (particles && IS(kw, "element")) {
        for (size_t t = 0; t < sys->atom.type.count && sys->atom.type.z; ++t) {
            const str_t n = md_util_element_symbol(sys->atom.type.z[t]);
            cmpl_add(c, n, n, MD_VIAL_COMPLETION_VALUE);
        }
    } else if ((particles && IS(kw, "resname")) || (residues && (IS(kw, "name") || !kw.len))) {
        for (size_t r = 0; r < md_system_residue_count(sys); ++r) {
            const str_t n = md_system_residue_name(sys, r);
            cmpl_add(c, n, n, MD_VIAL_COMPLETION_VALUE);
        }
    } else if (IS(kw, "chain") || (IS(callee, "chain") && (IS(kw, "id") || !kw.len))) {
        for (size_t i = 0; i < md_system_chain_count(sys); ++i) {
            const str_t n = md_system_chain_id(sys, i);
            cmpl_add(c, n, n, MD_VIAL_COMPLETION_VALUE);
        }
    } else if (IS(callee, "molecule") && (IS(kw, "name") || !kw.len)) {
        for (size_t i = 0; i < md_system_molecule_type_count(sys); ++i) {
            const str_t n = md_system_molecule_type_name(sys, i);
            cmpl_add(c, n, n, MD_VIAL_COMPLETION_VALUE);
        }
    }
    #undef IS
}

static int cmpl_rank(md_vial_completion_kind_t k) {
    switch (k) {
    case MD_VIAL_COMPLETION_PARAMETER: return 0;
    case MD_VIAL_COMPLETION_VALUE:     return 0;
    case MD_VIAL_COMPLETION_PROPERTY:  return 1;
    case MD_VIAL_COMPLETION_VARIABLE:  return 1;
    case MD_VIAL_COMPLETION_FUNCTION:  return 2;
    case MD_VIAL_COMPLETION_CONSTANT:  return 2;
    default:                           return 3;
    }
}

static int cmpl_cmp(const void* a, const void* b) {
    const md_vial_completion_t* x = a;
    const md_vial_completion_t* y = b;
    const int r = cmpl_rank(x->kind) - cmpl_rank(y->kind);
    if (r) return r;
    const size_t n = MIN(x->label.len, y->label.len);
    const int c = strncmp(x->label.ptr, y->label.ptr, n);
    return c ? c : (int)x->label.len - (int)y->label.len;
}

bool md_vial_complete(md_vial_completions_t* out, const md_vial_program_t* prog, uint32_t cursor, md_allocator_i* alloc) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    if (!prog || !alloc || cursor > prog->src.len) return false;
    const vial_token_t* toks = prog->tokens;
    const size_t nt = md_array_size(toks);
    cmpl_t c = {.alloc = alloc};
    // The token the cursor is in or at the end of
    size_t ti = SIZE_MAX;
    for (size_t i = 0; i < nt; ++i) {
        const vial_token_t* t = &toks[i];
        if (t->kind == MD_VIAL_TOKEN_EOF) break;
        if (t->span.beg < cursor && cursor <= t->span.end &&
            (t->kind == MD_VIAL_TOKEN_NAME || t->kind == MD_VIAL_TOKEN_KEYWORD || t->kind == MD_VIAL_TOKEN_STRING)) { ti = i; break; }
        if (t->span.beg >= cursor) break;
    }
    md_vial_span_t span = {cursor, cursor};
    bool in_string = false;
    if (ti != SIZE_MAX && toks[ti].kind == MD_VIAL_TOKEN_STRING) {
        // Inside the quotes: the string's contents are replaced
        const vial_token_t* t = &toks[ti];
        const str_t text = tok_text(prog, t);
        size_t q = 0;
        while (q < text.len && text.ptr[q] != '\'' && text.ptr[q] != '"') ++q;
        const size_t ql = (t->flags & VIAL_STR_TRIPLE) ? 3 : 1;
        const uint32_t beg = t->span.beg + (uint32_t)(q + ql);
        uint32_t end = t->span.end;
        if (!(t->flags & VIAL_STR_UNCLOSED) && end >= beg + ql) end -= (uint32_t)ql;
        if (cursor < beg || cursor > end) return true;     // on a quote
        span = (md_vial_span_t){beg, end};
        in_string = true;
    } else if (ti != SIZE_MAX) {
        span = toks[ti].span;
    }
    c.prefix = prog->src.ptr ? (str_t){prog->src.ptr + span.beg, cursor - span.beg} : (str_t){0};
    // The tokens before what is replaced, back to the call around it
    size_t prev = SIZE_MAX;
    for (size_t i = 0; i < nt && toks[i].span.end <= span.beg && toks[i].kind != MD_VIAL_TOKEN_EOF; ++i) {
        if (toks[i].kind == MD_VIAL_TOKEN_COMMENT) continue;
        if (ti != SIZE_MAX && i == ti) break;
        prev = i;
    }
    str_t kw = {0}, callee = {0};
    bool subscript = false;     // inside a subscript's brackets, s['|']: a particle's name
    bool after_dot = prev != SIZE_MAX && tok_op(&toks[prev], VIAL_OP_DOT);
    if (prev != SIZE_MAX && prev >= 1 && tok_op(&toks[prev], VIAL_OP_ASSIGN) && toks[prev - 1].kind == MD_VIAL_TOKEN_NAME) {
        kw = tok_text(prog, &toks[prev - 1]);
    }
    {
        int depth = 0;
        for (size_t j = prev; j != SIZE_MAX; --j) {
            const vial_token_t* t = &toks[j];
            if (t->kind == MD_VIAL_TOKEN_NEWLINE || t->kind == MD_VIAL_TOKEN_INDENT || t->kind == MD_VIAL_TOKEN_DEDENT) { if (!depth) break; continue; }
            if (tok_op(t, VIAL_OP_RPAREN) || tok_op(t, VIAL_OP_RBRACKET) || tok_op(t, VIAL_OP_RBRACE)) { ++depth; continue; }
            if (tok_op(t, VIAL_OP_LPAREN) || tok_op(t, VIAL_OP_LBRACKET) || tok_op(t, VIAL_OP_LBRACE)) {
                if (depth) { --depth; continue; }
                if (tok_op(t, VIAL_OP_LPAREN) && j >= 1 && toks[j - 1].kind == MD_VIAL_TOKEN_NAME) callee = tok_text(prog, &toks[j - 1]);
                if (tok_op(t, VIAL_OP_LBRACKET) && j >= 1 && (toks[j - 1].kind == MD_VIAL_TOKEN_NAME || tok_op(&toks[j - 1], VIAL_OP_RBRACKET) ||
                                                             tok_op(&toks[j - 1], VIAL_OP_RPAREN))) subscript = true;
                break;
            }
            if (j == 0) break;
        }
    }
    // A kwarg's value is not the call's own argument when an inner bracket opened since
    if (kw.len && !callee.len) kw = (str_t){0};

    if (in_string) {
        if (str_eq(callee, STR_LIT("attr")) && prog->check.sys && prog->check.run.len) {
            // attr()'s paths: the run's attributes
            const md_attributes_t* at = &prog->check.sys->attributes;
            for (size_t i = 0; i < md_array_size(at->attr); ++i) {
                const str_t path = md_vial_attr_path(prog->check.run, at->attr[i].path);
                cmpl_add(&c, path, path, MD_VIAL_COMPLETION_VALUE);
            }
        } else if (callee.len) {
            system_values(&c, prog->check.sys, callee, kw);
        } else if (subscript) {
            system_values(&c, prog->check.sys, STR_LIT("particle"), STR_LIT("name"));
        }
    } else if (after_dot) {
        for (size_t i = 0; i < vial_builtin_count(); ++i) {
            const vial_builtin_t* b = vial_builtin_at(i);
            if (b->name[0] == '_') continue;
            if (b->flags & VIAL_BUILTIN_PROPERTY) cmpl_add(&c, str_from_cstr(b->name), str_from_cstr(b->name), MD_VIAL_COMPLETION_PROPERTY);
        }
        // The selectors within each set, for .residue[i] and .particle(...)
        static const char* views[] = {"particle", "atom", "bead", "residue", "ring", "part", "molecule", "head", "tail"};
        for (size_t i = 0; i < ARRAY_SIZE(views); ++i) cmpl_add(&c, str_from_cstr(views[i]), str_from_cstr(views[i]), MD_VIAL_COMPLETION_PROPERTY);
    } else {
        // The call's parameters first
        if (callee.len && !kw.len) {
            const vial_node_id call = call_at(prog, cursor, callee);
                    if (call) call_params(&c, prog, call);
        }
        // Names the script binds where the cursor is: top-level names bound before the line it is on, and the
        // parameters, locals and loop variables of what holds it
        const vial_names_t* nm = &prog->check.names;
        // The top-level statement the cursor is in: names it and later ones bind are not bound yet
        uint32_t here = 0;
        {
            const vial_tree_t* t = &prog->tree;
            const vial_node_id root = t->root;
            for (uint32_t k = 0; root && k < t->nodes[root].child_count; ++k) {
                const vial_node_id st = vial_child(t, root, k);
                if (st && t->nodes[st].span.beg < span.beg) here = k;
                else if (st && t->nodes[st].span.beg >= span.beg) break;
            }
            // Past the end of the last statement and on a line of its own: a new statement
            if (root && t->nodes[root].child_count) {
                const vial_node_id st = vial_child(t, root, here);
                if (st && span.beg > t->nodes[st].span.end) {
                    bool newline = false;
                    for (uint32_t q = t->nodes[st].span.end; q < span.beg; ++q) newline |= prog->src.ptr[q] == '\n';
                    if (newline) here += 1;
                }
            }
        }
        for (size_t i = 1; i < md_array_size(nm->symbols); ++i) {
            const vial_symbol_t* s = &nm->symbols[i];
            if (s->in_prelude || s->kind == VIAL_SYM_BUILTIN || s->kind == VIAL_SYM_PRELUDE || !s->name.len) continue;
            if (!scope_holds(prog, s->scope, cursor)) continue;
            if ((s->kind == VIAL_SYM_GLOBAL || s->kind == VIAL_SYM_DEF) && s->stmt >= here) continue;
            cmpl_add(&c, s->name, s->name, s->kind == VIAL_SYM_DEF ? MD_VIAL_COMPLETION_FUNCTION : MD_VIAL_COMPLETION_VARIABLE);
        }
        for (size_t i = 0; i < vial_builtin_count(); ++i) {
            const vial_builtin_t* b = vial_builtin_at(i);
            if (b->name[0] == '_') continue;
            if (b->flags & (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_SELECTOR)) cmpl_add(&c, str_from_cstr(b->name), str_from_cstr(b->name), MD_VIAL_COMPLETION_FUNCTION);
            else if (b->flags & VIAL_BUILTIN_VALUE) cmpl_add(&c, str_from_cstr(b->name), str_from_cstr(b->name), MD_VIAL_COMPLETION_CONSTANT);
        }
        // Stored selections
        for (size_t i = 0; i < prog->check.num_stored; ++i) cmpl_add(&c, prog->check.stored[i].name, prog->check.stored[i].name, MD_VIAL_COMPLETION_CONSTANT);
        str_t defs[64];
        const size_t nd = MIN(prelude_defs(defs, 64), 64);
        for (size_t i = 0; i < nd; ++i) cmpl_add(&c, defs[i], defs[i], MD_VIAL_COMPLETION_FUNCTION);
        // The keywords a script writes; the others are reserved and E08xx
        static const vial_kw_t kws[] = {VIAL_KW_FALSE, VIAL_KW_NONE_, VIAL_KW_TRUE, VIAL_KW_AND, VIAL_KW_ASSERT, VIAL_KW_DEF, VIAL_KW_ELIF,
                                        VIAL_KW_ELSE, VIAL_KW_FOR, VIAL_KW_IF, VIAL_KW_IN, VIAL_KW_IS, VIAL_KW_NOT, VIAL_KW_OR, VIAL_KW_RETURN};
        for (size_t i = 0; i < ARRAY_SIZE(kws); ++i) {
            const str_t w = str_from_cstr(vial_kw_str(kws[i]));
            cmpl_add(&c, w, w, MD_VIAL_COMPLETION_KEYWORD);
        }
    }
    out->span = span;
    out->prefix = str_copy(c.prefix, alloc);
    out->count = md_array_size(c.items);
    if (out->count) qsort(c.items, out->count, sizeof(*c.items), cmpl_cmp);
    out->items = c.items;
    return true;
}

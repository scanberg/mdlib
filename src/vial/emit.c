// Lowering: the check emits the array IR as it types each node (see infer_internal.h). The checker has already made
// every decision lowering needs, so this file only translates them: which built-in a call is and how its arguments
// bound, which branch a static if takes, which specialization of a def a call uses. A node whose value the check
// knows becomes a constant, so a selection is evaluated once, when checking.

#include "infer_internal.h"
#include "ir.h"
#include "builtins.h"

#include <core/md_allocator.h>
#include <core/md_array.h>

#include <stdio.h>
#include <string.h>

typedef vial_type_t T;

static inline uint16_t sub_of(const infer_t* x, vial_node_id id) { return x->t->nodes[id].sub; }

// ### PLUMBING ###

void vx_sizes_ssa(infer_t* x, vial_types_t* out, size_t num_nodes, size_t num_syms) {
    (void)x;
    md_array_resize(out->ssa, num_nodes, x->alloc);
    md_array_resize(out->bind_of, num_nodes, x->alloc);
    memset(out->ssa, 0, num_nodes * sizeof(*out->ssa));
    memset(out->bind_of, 0, num_nodes * sizeof(*out->bind_of));
    if (num_syms) {
        md_array_resize(out->sym_ssa, num_syms, x->alloc);
        memset(out->sym_ssa, 0, num_syms * sizeof(*out->sym_ssa));
    }
}

vial_ssa_t vx_ssa(const infer_t* x, vial_node_id id) {
    if (!id || !x->out->ssa) return 0;
    return x->out->ssa[id - x->out->base];
}

void vx_set_ssa(infer_t* x, vial_node_id id, vial_ssa_t v) {
    if (id && x->out->ssa) x->out->ssa[id - x->out->base] = v;
}

md_vial_span_t vx_span(const infer_t* x, vial_node_id id) {
    // A def of the prelude has no span in the script: its ops have the span of the call that led there
    if (x->t != x->c->tree) {
        const vx_frame_t* f = x->frame;
        while (f && !f->site) f = f->parent;
        return f && f->site ? x->c->tree->nodes[f->site].span : (md_vial_span_t){0};
    }
    return id ? x->t->nodes[id].span : (md_vial_span_t){0};
}

void vx_lower_fail(infer_t* x, vial_node_id id, const char* what) {
    if (!x->ir) return;
    for (vx_frame_t* f = x->frame; f; f = f->parent) f->lower_failed = true;
    const md_vial_span_t sp = vx_span(x, id);
    vial_ir_failure_t f = {.span = sp};
    if (x->t != x->c->tree && id) {
        // In the prelude: which part of it, as the span is the script's call
        const md_vial_span_t ps = x->t->nodes[id].span;
        const int len = (int)(ps.end - ps.beg < 40 ? ps.end - ps.beg : 40);
        snprintf(f.what, sizeof(f.what), "%s, in the prelude's '%.*s'", what, len, x->c->prelude_src.ptr + ps.beg);
    } else {
        snprintf(f.what, sizeof(f.what), "%s", what);
    }
    for (size_t i = 0; i < md_array_size(x->ir->failures); ++i) {
        const vial_ir_failure_t* g = &x->ir->failures[i];
        if (g->span.beg == sp.beg && g->span.end == sp.end && strcmp(g->what, f.what) == 0) return;
    }
    md_array_push(x->ir->failures, f, x->ir->alloc);
}

vial_ssa_t vx_sym_ssa(const infer_t* x, uint32_t sym) {
    for (const vx_frame_t* f = x->frame; f; f = f->parent) {
        if (f->t != x->t || !f->syms.num_buckets) continue;
        const uint32_t* s = md_hashmap_get((md_hashmap32_t*)&f->syms, (uint64_t)sym);
        if (s) return *s < md_array_size(f->sym_ssa) ? f->sym_ssa[*s] : 0;
    }
    if (x->t == x->c->tree && x->top->sym_ssa && sym < md_array_size(x->top->sym_ssa)) return x->top->sym_ssa[sym];
    return 0;
}

void vx_frame_set_ssa(infer_t* x, uint32_t sym, vial_ssa_t v) {
    vx_frame_t* f = x->frame;
    if (!f) {
        if (x->top->sym_ssa && sym < md_array_size(x->top->sym_ssa)) x->top->sym_ssa[sym] = v;
        return;
    }
    const uint32_t* s = f->syms.num_buckets ? md_hashmap_get(&f->syms, (uint64_t)sym) : NULL;
    if (!s) return;
    while (md_array_size(f->sym_ssa) <= *s) md_array_push(f->sym_ssa, 0, x->alloc);
    f->sym_ssa[*s] = v;
}

void vx_record_sig(infer_t* x, vial_node_id call, const vx_param_t* params, uint32_t np) {
    if (!x->top || !x->c || x->t != x->c->tree) return;
    vial_types_t* t = x->top;
    const size_t nn = md_array_size(x->t->nodes);
    if (call >= nn) return;
    const size_t had = md_array_size(t->sigs);
    if (had < nn) {
        md_array_resize(t->sigs, nn, x->alloc);
        memset(t->sigs + had, 0, (nn - had) * sizeof(vx_sig_t));
    }
    t->sigs[call] = (vx_sig_t){.params = params, .np = np};
}

void vx_record_bind(infer_t* x, vial_node_id call, const vx_param_t* params, uint32_t np, const vial_node_id* args) {
    if (!x->ir || !x->out->bind_of || np > 12) return;
    vx_bound_t b = {.params = params, .np = np};
    memcpy(b.args, args, np * sizeof(vial_node_id));
    const uint32_t at = x->out->bind_of[call - x->out->base];
    if (at) { x->out->binds[at - 1] = b; return; }
    md_array_push(x->out->binds, b, x->alloc);
    x->out->bind_of[call - x->out->base] = (uint32_t)md_array_size(x->out->binds);
}

// ### CONSTANTS ###

// The check's value of a node as a runtime value of type t; false when it has none, or none that is complete
static bool const_of(infer_t* x, const vial_value_t* v, const T* t, vial_val_t* out) {
    md_allocator_i* a = x->ir->alloc;
    if (t->stage != VIAL_STAGE_STATIC) return false;
    switch (v->kind) {
    case VIAL_VAL_BOOL: case VIAL_VAL_INT: case VIAL_VAL_FLOAT: {
        if (t->ndim) return false;
        const double d = v->kind == VIAL_VAL_FLOAT ? v->f : v->kind == VIAL_VAL_INT ? (double)v->i : (double)v->b;
        switch (t->dtype) {
        case VIAL_DT_BOOL:    *out = vial_val_bool(d != 0.0, a); return true;
        case VIAL_DT_INT:     if (d < INT32_MIN || d > INT32_MAX) return false; *out = vial_val_int((int32_t)d, a); return true;
        case VIAL_DT_FLOAT:   *out = vial_val_float((float)d, a); return true;
        case VIAL_DT_FLOAT64: *out = vial_val_float64(d, a); return true;
        default: return false;
        }
    }
    case VIAL_VAL_STR:
        if (t->dtype != VIAL_DT_STR || t->ndim) return false;
        *out = vial_val_str(v->s, a);
        return true;
    case VIAL_VAL_NULL:
        *out = vial_val_none();
        return true;
    case VIAL_VAL_SETS: {
        if (t->dtype != VIAL_DT_SET || !v->sets) return false;
        *out = vial_val_sets(v->sets);
        // An array of sets of several dimensions, set[k][4], keeps them: the sets are listed in row-major order
        size_t n = 1;
        bool known = t->ndim > 1;
        for (uint8_t d = 0; d < t->ndim && known; ++d) {
            if (t->dims[d].kind != VIAL_EXT_STATIC) known = false;
            else n *= t->dims[d].n;
        }
        if (known && !v->sets->single && n == v->sets->count) {
            out->ndim = t->ndim;
            for (uint8_t d = 0; d < t->ndim; ++d) out->dims[d] = t->dims[d].n;
        }
        return true;
    }
    case VIAL_VAL_STRS: {
        const uint32_t n = v->n;
        *out = vial_val_make(VIAL_DT_STR, 1, &n, a);
        for (uint32_t i = 0; i < n; ++i) ((str_t*)out->data)[i] = str_copy(v->strs[i], a);
        return true;
    }
    case VIAL_VAL_LABELS: {
        if (!x->sys || !v->sets || !v->sets->unit) return false;
        const uint32_t n = v->sets->count;
        *out = vial_val_make(VIAL_DT_STR, 1, &n, a);
        for (uint32_t i = 0; i < n; ++i) {
            char buf[64];
            const size_t len = vial_unit_label(buf, sizeof(buf), x->sys, v->sets->level, v->sets->unit[i]);
            ((str_t*)out->data)[i] = str_copy((str_t){buf, len}, a);
        }
        return true;
    }
    case VIAL_VAL_SHAPE: {
        // A tuple of the extents, when every one is known
        if (!v->type) return false;
        const uint32_t n = v->type->ndim;
        vial_val_t* el = md_alloc(a, (n ? n : 1) * sizeof(vial_val_t));
        for (uint32_t i = 0; i < n; ++i) {
            if (v->type->dims[i].kind != VIAL_EXT_STATIC) return false;
            el[i] = vial_val_int((int32_t)v->type->dims[i].n, a);
        }
        *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = n};
        return true;
    }
    default:
        return false;
    }
}

static str_t node_text(const infer_t* x, vial_node_id id) {
    const md_vial_span_t s = x->t->nodes[id].span;
    const str_t src = x->t == x->c->tree ? x->c->src : x->c->prelude_src;
    return (str_t){src.ptr + s.beg, s.end - s.beg};
}

// ### NODES ###

static vial_ssa_t op1(infer_t* x, vial_node_id id, uint16_t kind, uint16_t sub, uint32_t aux, const vial_ssa_t* args, uint32_t n, const T* t) {
    for (uint32_t i = 0; i < n; ++i) if (!args[i]) return 0;
    return vial_ir_op1(x->ir, kind, sub, aux, args, n, t, vx_span(x, id));
}

static const T* type_at(const infer_t* x, vial_node_id id) { return &x->out->types[id - x->out->base]; }

// b converted to a's unit for + - % and comparisons, as the check folds it: compatible units convert to the left's
static vial_ssa_t in_unit_of(infer_t* x, vial_node_id at, const T* a, vial_ssa_t bv, const T* b) {
    if (a->unit_kind != VIAL_UNIT_KNOWN || b->unit_kind != VIAL_UNIT_KNOWN) return bv;
    double f;
    if (!md_unit_conversion_factor(&f, b->unit, a->unit) || f == 1.0) return bv;
    T ft = vial_type_scalar(b->dtype == VIAL_DT_FLOAT64 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT);
    vial_val_t fv = ft.dtype == VIAL_DT_FLOAT64 ? vial_val_float64(f, x->ir->alloc) : vial_val_float((float)f, x->ir->alloc);
    const vial_ssa_t c = vial_ir_const(x->ir, &fv, &ft, (str_t){0}, vx_span(x, at));
    T r = *b;
    if (r.dtype != VIAL_DT_FLOAT64) r.dtype = VIAL_DT_FLOAT;
    r.unit = a->unit;
    return op1(x, at, VIAL_IR_BINARY, VIAL_OP_STAR, 0, (vial_ssa_t[]){bv, c}, 2, &r);
}

// An operand evaluated as a float64 beside one (vx_operand_dtype): a bare float number is emitted again as a float64
// constant of the value the check knows, so 0.1 beside a float64 is 0.1 and not the float nearest it
static vial_ssa_t in_operand_dtype(infer_t* x, vial_node_id n, vial_ssa_t v, const T* self, const T* other) {
    if (!v || vx_operand_dtype(self, other) == self->dtype) return v;
    const vial_value_t val = vx_value_of(x, n);
    double d;
    if (val.kind == VIAL_VAL_FLOAT) d = val.f;
    else if (val.kind == VIAL_VAL_INT) d = (double)val.i;
    else return v;
    T ft = *self;
    ft.dtype = VIAL_DT_FLOAT64;
    const vial_val_t fv = vial_val_float64(d, x->ir->alloc);
    return vial_ir_const(x->ir, &fv, &ft, (str_t){0}, vx_span(x, n));
}

static vial_ssa_t emit_binary(infer_t* x, vial_node_id id, const T* t) {
    const vial_node_id an = vx_child(x, id, 0), bn = vx_child(x, id, 1);
    const uint16_t op = sub_of(x, id);
    vial_ssa_t a = vx_ssa(x, an), b = vx_ssa(x, bn);
    a = in_operand_dtype(x, an, a, type_at(x, an), type_at(x, bn));
    b = in_operand_dtype(x, bn, b, type_at(x, bn), type_at(x, an));
    if (op == VIAL_OP_PLUS || op == VIAL_OP_MINUS || op == VIAL_OP_PERCENT) b = in_unit_of(x, id, type_at(x, an), b, type_at(x, bn));
    return op1(x, id, VIAL_IR_BINARY, op, 0, (vial_ssa_t[]){a, b}, 2, t);
}

static vial_ssa_t emit_compare(infer_t* x, vial_node_id id, const T* t) {
    // A chain is the and of its comparisons: a < b < c is a < b and b < c, each operand evaluated once
    const uint32_t n = vx_num_children(x, id);
    vial_ssa_t acc = 0;
    for (uint32_t i = 0; i + 2 < n; i += 2) {
        const vial_node_id an = vx_child(x, id, i), on = vx_child(x, id, i + 1), bn = vx_child(x, id, i + 2);
        const uint16_t cmp = x->t->nodes[on].sub;
        const T* at = type_at(x, an);
        const T* bt = type_at(x, bn);
        if (cmp == VIAL_CMP_IS || cmp == VIAL_CMP_IS_NOT || bt->dtype == VIAL_DT_TUPLE || at->dtype == VIAL_DT_NONE || bt->dtype == VIAL_DT_NONE) {
            vx_lower_fail(x, on, "a comparison whose value is known only when compiling");
            return 0;
        }
        T pt = *t;
        pt.ndim = 0;
        if (!vx_broadcast(at, bt, &pt) && cmp != VIAL_CMP_IN && cmp != VIAL_CMP_NOT_IN) pt = *t;
        if (cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) { pt.ndim = at->ndim; memcpy(pt.dims, at->dims, sizeof(pt.dims)); }
        pt.dtype = VIAL_DT_BOOL;
        pt.unit_kind = VIAL_UNIT_KNOWN;
        pt.unit = md_unit_none();
        vial_ssa_t a = vx_ssa(x, an), b = vx_ssa(x, bn);
        a = in_operand_dtype(x, an, a, at, bt);
        b = in_operand_dtype(x, bn, b, bt, at);
        if (bt->dtype != VIAL_DT_RANGE) b = in_unit_of(x, on, at, b, bt);
        const vial_ssa_t c = op1(x, on, VIAL_IR_COMPARE, cmp, 0, (vial_ssa_t[]){a, b}, 2, i + 3 >= n && !acc ? t : &pt);
        if (!c) return 0;
        acc = acc ? op1(x, id, VIAL_IR_LOGIC, VIAL_KW_AND, 0, (vial_ssa_t[]){acc, c}, 2, i + 3 >= n ? t : &pt) : c;
    }
    return acc;
}

static vial_ssa_t scope_of(infer_t* x, vial_node_id id) {
    // The particles `not` takes from: every particle of the system, the scope of a script
    if (!x->sys) { vx_lower_fail(x, id, "`not` of a set without a system"); return 0; }
    if (!x->all) x->all = vial_sets_all(x->sys->num_atoms, x->alloc);
    const vial_val_t v = vial_val_sets(x->all);
    T t = vial_type_scalar(VIAL_DT_SET);
    return vial_ir_const(x->ir, &v, &t, STR_LIT("system"), vx_span(x, id));
}

static vial_ssa_t emit_boolop(infer_t* x, vial_node_id id, const T* t) {
    const uint32_t n = vx_num_children(x, id);
    const uint16_t kw = sub_of(x, id);
    if (t->dtype != VIAL_DT_SET) {
        vial_ssa_t acc = 0;
        for (uint32_t i = 0; i < n; ++i) {
            const vial_ssa_t c = vx_ssa(x, vx_child(x, id, i));
            if (!c) return 0;
            acc = acc ? op1(x, id, VIAL_IR_LOGIC, kw, 0, (vial_ssa_t[]){acc, c}, 2, t) : c;
        }
        return acc;
    }
    // Sets: under an `and`, `not b` subtracts b, as the check computes it
    vial_ssa_t acc = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const vial_node_id cn = vx_child(x, id, i);
        if (kw == VIAL_KW_AND && vx_kind(x, cn) == VIAL_N_NOT) continue;
        const vial_ssa_t c = vx_ssa(x, cn);
        if (!c) return 0;
        acc = acc ? op1(x, id, VIAL_IR_SETOP, kw, 0, (vial_ssa_t[]){acc, c}, 2, t) : c;
    }
    if (kw == VIAL_KW_AND) {
        for (uint32_t i = 0; i < n; ++i) {
            const vial_node_id cn = vx_child(x, id, i);
            if (vx_kind(x, cn) != VIAL_N_NOT) continue;
            if (!acc) acc = scope_of(x, id);
            const vial_ssa_t b = vx_ssa(x, vx_child(x, cn, 0));
            if (!acc || !b) return 0;
            acc = op1(x, id, VIAL_IR_SETOP, VIAL_OP_MINUS, 0, (vial_ssa_t[]){acc, b}, 2, t);
        }
    }
    return acc;
}

static vial_ssa_t emit_not(infer_t* x, vial_node_id id, const T* t) {
    const vial_ssa_t a = vx_ssa(x, vx_child(x, id, 0));
    if (t->dtype != VIAL_DT_SET) return op1(x, id, VIAL_IR_NOT, 0, 0, &a, 1, t);
    const vial_node_id up = x->n->parent[id];
    if (up && vx_kind(x, up) == VIAL_N_BOOLOP && x->t->nodes[up].sub == VIAL_KW_AND) return a;   // the `and` subtracts it
    const vial_ssa_t s = scope_of(x, id);
    return op1(x, id, VIAL_IR_SETOP, VIAL_KW_NOT, 0, (vial_ssa_t[]){s, a}, 2, t);
}

static vial_ssa_t emit_items(infer_t* x, vial_node_id id, uint16_t kind, const T* t) {
    const uint32_t n = vx_num_children(x, id);
    vial_ssa_t args[64];
    if (n > 64) { vx_lower_fail(x, id, "a list of more than 64 items written out"); return 0; }
    for (uint32_t i = 0; i < n; ++i) {
        args[i] = vx_ssa(x, vx_child(x, id, i));
        if (!args[i]) return 0;
    }
    return vial_ir_op1(x->ir, kind, 0, 0, args, n, t, vx_span(x, id));
}

static vial_ssa_t emit_subscript(infer_t* x, vial_node_id id, const T* t) {
    const vial_node_id on = vx_child(x, id, 0), in = vx_child(x, id, 1);
    const vial_ssa_t a = vx_ssa(x, on);
    if (!a || !in) return 0;
    const T* at = type_at(x, on);
    if (at->dtype == VIAL_DT_TUPLE) {
        const vial_value_t k = vx_value_of(x, in);
        if (k.kind != VIAL_VAL_INT) return 0;
        const int64_t i = k.i < 0 ? k.i + at->num_elems : k.i;
        return op1(x, id, VIAL_IR_ELEM, 0, (uint32_t)i, &a, 1, t);
    }
    vial_node_id items[16];
    uint32_t num = 0;
    if (vx_kind(x, in) == VIAL_N_TUPLE) {
        for (uint32_t i = 0; i < vx_num_children(x, in) && num < 16; ++i) items[num++] = vx_child(x, in, i);
    } else {
        items[num++] = in;
    }
    vial_ssa_t args[1 + 3 * 16];
    uint32_t na = 0;
    args[na++] = a;
    const uint32_t first = (uint32_t)md_array_size(x->ir->index_items);
    for (uint32_t i = 0; i < num; ++i) {
        const vial_node_id it = items[i];
        const vial_node_kind_t k = vx_kind(x, it);
        vial_index_item_t item = {0};
        if (k == VIAL_N_ELLIPSIS) item.kind = VIAL_INDEX_ELLIPSIS;
        else if (k == VIAL_N_CONST && x->t->nodes[it].sub == VIAL_KW_NONE_) item.kind = VIAL_INDEX_NEWAXIS;
        else if (k == VIAL_N_SLICE) {
            item.kind = VIAL_INDEX_SLICE;
            const vial_node_id parts[3] = {vx_child(x, it, 0), vx_child(x, it, 1), vx_child(x, it, 2)};
            vial_ssa_t* dst[3] = {&item.a, &item.b, &item.c};
            for (int q = 0; q < 3; ++q) {
                if (!parts[q]) continue;
                *dst[q] = vx_ssa(x, parts[q]);
                if (!*dst[q]) return 0;
                args[na++] = *dst[q];
            }
        } else {
            const T* it_t = type_at(x, it);
            item.kind = it_t->ndim == 0 && it_t->dtype == VIAL_DT_INT ? VIAL_INDEX_SCALAR : VIAL_INDEX_ARRAY;
            if (it_t->dtype == VIAL_DT_SET && !(at->dtype == VIAL_DT_SET)) {
                // a[s]: the values of s's particles along a's dimension over particles, which runs over known particles
                int pd = -1;
                for (uint8_t q = 0; q < at->ndim && pd < 0; ++q) {
                    const uint32_t dsc = at->dims[q].desc;
                    if (dsc && dsc < md_array_size(*x->descs) && (*x->descs)[dsc].level == VIAL_LEVEL_PARTICLE) pd = q;
                }
                const vial_sets_t* over = pd >= 0 ? (*x->descs)[at->dims[pd].desc].sets : NULL;
                if (!over || !over->single) {
                    vx_lower_fail(x, it, "a set as a position, along particles not known when compiling");
                    return 0;
                }
                const vial_val_t ov = vial_val_sets(over);
                const T ot = vial_type_scalar(VIAL_DT_SET);
                const vial_val_t dv = vial_val_int(pd, x->ir->alloc);
                const T dt = vial_type_scalar(VIAL_DT_INT);
                item.kind = VIAL_INDEX_SET;
                item.b = vial_ir_const(x->ir, &ov, &ot, (str_t){0}, vx_span(x, it));
                item.c = vial_ir_const(x->ir, &dv, &dt, (str_t){0}, vx_span(x, it));
                item.a = vx_ssa(x, it);
                if (!item.a) return 0;
                args[na++] = item.a;
                args[na++] = item.b;
                args[na++] = item.c;
                md_array_push(x->ir->index_items, item, x->ir->alloc);
                continue;
            }
            item.a = vx_ssa(x, it);
            if (!item.a) return 0;
            args[na++] = item.a;
        }
        md_array_push(x->ir->index_items, item, x->ir->alloc);
    }
    const uint32_t op = vial_ir_emit(x->ir, &(vial_ir_emit_t){.kind = VIAL_IR_INDEX, .sub = (uint16_t)num, .aux = first, .args = args,
                                                              .num_args = na, .types = t, .num_results = 1, .span = vx_span(x, id)});
    return x->ir->ops[op].result;
}

static vial_ssa_t builtin(infer_t* x, vial_node_id id, const char* name, const vial_ssa_t* args, const str_t* names, uint32_t n, const T* t) {
    const vial_builtin_t* b = vial_builtin_find(str_from_cstr(name));
    if (!b) return 0;
    const uint32_t op = vial_ir_emit(x->ir, &(vial_ir_emit_t){.kind = VIAL_IR_BUILTIN, .aux = vial_builtin_index(b), .args = args,
                                                              .names = names, .num_args = n, .types = t, .num_results = 1, .span = vx_span(x, id)});
    return x->ir->ops[op].result;
}

static vial_ssa_t emit_attr(infer_t* x, vial_node_id id, const T* t) {
    const vial_node_id on = vx_child(x, id, 0), nn = vx_child(x, id, 1);
    const vial_ssa_t a = vx_ssa(x, on);
    if (!a || !nn) return 0;
    const str_t name = node_text(x, nn);
    char buf[16];
    snprintf(buf, sizeof(buf), "%.*s", (int)(name.len < 15 ? name.len : 15), name.ptr);
    // .xyz of a set is its positions; of coordinates, itself
    const T* at = type_at(x, on);
    if (str_eq(name, STR_LIT("xyz")) && at->dtype != VIAL_DT_SET) return a;
    if (str_eq(name, STR_LIT("xyz"))) return builtin(x, id, "positions", &a, NULL, 1, t);
    return builtin(x, id, buf, &a, NULL, 1, t);
}

// A tuple or list written as an argument whose items the built-in typed one by one, as concatenate's: a tuple of them
static vial_ssa_t items_as_tuple(infer_t* x, vial_node_id an) {
    if (vx_kind(x, an) != VIAL_N_TUPLE && vx_kind(x, an) != VIAL_N_LIST) return 0;
    const uint32_t n = vx_num_children(x, an);
    if (n > 64) return 0;
    vial_ssa_t items[64];
    T* elems = md_alloc(x->alloc, (n ? n : 1) * sizeof(T));
    for (uint32_t i = 0; i < n; ++i) {
        const vial_node_id c = vx_child(x, an, i);
        items[i] = vx_ssa(x, c);
        if (!items[i]) return 0;
        elems[i] = *type_at(x, c);
    }
    T tt = vial_type_scalar(VIAL_DT_TUPLE);
    for (uint32_t i = 0; i < n; ++i) vx_join_stage(&tt, &tt, &elems[i]);
    tt.elems = elems;
    tt.num_elems = n;
    return vial_ir_op1(x->ir, VIAL_IR_TUPLE, 0, 0, items, n, &tt, vx_span(x, an));
}

static vial_ssa_t emit_call(infer_t* x, vial_node_id id, const T* t) {
    const vial_node_id callee = vx_child(x, id, 0);
    const vial_symbol_t* sym = callee && vx_kind(x, callee) == VIAL_N_NAME ? vial_names_symbol(x->n, callee) : NULL;
    if (!sym || sym->kind != VIAL_SYM_BUILTIN || !sym->builtin) {
        vx_lower_fail(x, id, "a call that is not lowered yet");
        return 0;
    }
    const uint32_t bi = x->out->bind_of ? x->out->bind_of[id - x->out->base] : 0;
    vial_ssa_t args[16];
    str_t names[16];
    uint32_t n = 0;
    if (bi) {
        const vx_bound_t* b = &x->out->binds[bi - 1];
        for (uint32_t k = 0; k < b->np && n < 16; ++k, ++n) {
            names[n] = str_from_cstr(b->params[k].name);
            args[n] = b->args[k] ? vx_ssa(x, b->args[k]) : 0;
            if (b->args[k] && !args[n]) return 0;
        }
    } else {
        // As written: positional arguments, then keywords by name
        for (uint32_t i = 1; i < vx_num_children(x, id) && n < 16; ++i) {
            const vial_node_id an = vx_child(x, id, i);
            if (!an) continue;
            if (vx_kind(x, an) == VIAL_N_KWARG) {
                const vial_node_id kn = vx_child(x, an, 0), vn = vx_child(x, an, 1);
                names[n] = kn ? node_text(x, kn) : (str_t){0};
                args[n] = vx_ssa(x, vn);
            } else {
                names[n] = (str_t){0};
                args[n] = vx_ssa(x, an);
                if (!args[n]) args[n] = items_as_tuple(x, an);
            }
            if (!args[n]) return 0;
            ++n;
        }
    }
    const uint32_t op = vial_ir_emit(x->ir, &(vial_ir_emit_t){.kind = VIAL_IR_BUILTIN, .aux = vial_builtin_index(sym->builtin), .args = args,
                                                              .names = names, .num_args = n, .types = t, .num_results = 1, .span = vx_span(x, id)});
    return x->ir->ops[op].result;
}

// Whether something inside span already failed to lower, so a failure here would only repeat it
static bool failed_within(const infer_t* x, md_vial_span_t span, bool prelude) {
    if (prelude) return false;      // every part of a prelude def has the call's span
    for (size_t i = 0; i < md_array_size(x->ir->failures); ++i) {
        const md_vial_span_t f = x->ir->failures[i].span;
        if (f.beg >= span.beg && f.end <= span.end) return true;
    }
    return false;
}

void vx_emit(infer_t* x, vial_node_id id, const T* t) {
    if (!x->ir || !id || vial_type_is_error(t)) return;
    const vial_node_kind_t k = vx_kind(x, id);
    // What the check knows is a constant, even when a handler lowered it already
    const vial_value_t v = vx_value_of(x, id);
    vial_val_t c;
    // A static tuple of arrays of sets whose dimensions describe exactly their sets, as common() gives: a constant
    if (t->dtype == VIAL_DT_TUPLE && t->stage == VIAL_STAGE_STATIC && k == VIAL_N_CALL && t->num_elems && t->num_elems <= 8) {
        vial_val_t* el = md_alloc(x->ir->alloc, t->num_elems * sizeof(vial_val_t));
        bool all = true;
        for (uint32_t i = 0; i < t->num_elems && all; ++i) {
            const T* e = &t->elems[i];
            const uint32_t dsc = e->ndim == 1 ? e->dims[0].desc : 0;
            const vial_sets_t* ss = dsc && dsc < md_array_size(*x->descs) ? (*x->descs)[dsc].sets : NULL;
            all = e->dtype == VIAL_DT_SET && ss && !ss->single && e->dims[0].kind == VIAL_EXT_STATIC && ss->count == e->dims[0].n;
            if (all) el[i] = vial_val_sets(ss);
        }
        if (all) {
            const vial_val_t tv = {.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = t->num_elems};
            vx_set_ssa(x, id, vial_ir_const(x->ir, &tv, t, node_text(x, id), vx_span(x, id)));
            return;
        }
    }
    const vial_symbol_t* nsym = k == VIAL_N_NAME ? vial_names_symbol(x->n, id) : NULL;
    const bool user_name = k == VIAL_N_NAME && !(nsym && nsym->kind == VIAL_SYM_BUILTIN);
    if (v.kind != VIAL_VAL_NONE && !user_name && const_of(x, &v, t, &c)) {
        const str_t text = k == VIAL_N_NAME || k == VIAL_N_CALL || k == VIAL_N_BOOLOP || k == VIAL_N_NOT || k == VIAL_N_SUBSCRIPT ? node_text(x, id) : (str_t){0};
        vx_set_ssa(x, id, vial_ir_const(x->ir, &c, t, text, vx_span(x, id)));
        return;
    }
    if (vx_ssa(x, id)) return;
    vial_ssa_t r = 0;
    switch (k) {
    case VIAL_N_NAME: {
        const uint32_t sym = x->n->ref[id];
        const vial_symbol_t* s = vial_names_symbol(x->n, id);
        if (s && s->kind == VIAL_SYM_BUILTIN && s->builtin) {
            // A selector written bare, all its units, or a predefined set, selected when evaluated
            r = vial_ir_op1(x->ir, VIAL_IR_BUILTIN, 0, vial_builtin_index(s->builtin), NULL, 0, t, vx_span(x, id));
            break;
        }
        if (t->dtype == VIAL_DT_DEF) {
            // A def passed as a value: each one passed has its own specialization, so nothing is passed when evaluating
            const vial_val_t none = vial_val_none();
            r = vial_ir_const(x->ir, &none, t, node_text(x, id), vx_span(x, id));
            break;
        }
        r = sym ? vx_sym_ssa(x, sym) : 0;
        if (!r) vx_lower_fail(x, id, "a name whose value is not lowered");
        break;
    }
    case VIAL_N_RANGE: {
        // a..b, both included: a constant of two numbers
        const vial_value_t a = vx_value_of(x, vx_child(x, id, 0)), b = vx_value_of(x, vx_child(x, id, 1));
        double lo, hi;
        if ((a.kind == VIAL_VAL_INT || a.kind == VIAL_VAL_FLOAT) && (b.kind == VIAL_VAL_INT || b.kind == VIAL_VAL_FLOAT)) {
            lo = a.kind == VIAL_VAL_INT ? (double)a.i : a.f;
            hi = b.kind == VIAL_VAL_INT ? (double)b.i : b.f;
            vial_val_t rv = vial_val_make(VIAL_DT_RANGE, 0, NULL, x->ir->alloc);
            ((double*)rv.data)[0] = lo;
            ((double*)rv.data)[1] = hi;
            r = vial_ir_const(x->ir, &rv, t, node_text(x, id), vx_span(x, id));
        } else {
            vx_lower_fail(x, id, "a range whose ends change per frame");
        }
        break;
    }
    case VIAL_N_QUANTITY: r = vx_ssa(x, vx_child(x, id, 0)); break;
    case VIAL_N_UNARY: {
        const vial_ssa_t a = vx_ssa(x, vx_child(x, id, 0));
        r = sub_of(x, id) == VIAL_OP_PLUS ? a : op1(x, id, VIAL_IR_UNARY, sub_of(x, id), 0, &a, 1, t);
        break;
    }
    case VIAL_N_NOT:      r = emit_not(x, id, t); break;
    case VIAL_N_BINARY:   r = emit_binary(x, id, t); break;
    case VIAL_N_COMPARE:  r = emit_compare(x, id, t); break;
    case VIAL_N_BOOLOP:   r = emit_boolop(x, id, t); break;
    case VIAL_N_LIST:     r = emit_items(x, id, VIAL_IR_LIST, t); break;
    case VIAL_N_TUPLE:    r = emit_items(x, id, VIAL_IR_TUPLE, t); break;
    case VIAL_N_SUBSCRIPT: r = emit_subscript(x, id, t); break;
    case VIAL_N_ATTR:     r = emit_attr(x, id, t); break;
    case VIAL_N_CALL:     r = emit_call(x, id, t); break;
    case VIAL_N_TERNARY: {
        // A static condition takes one branch; the check typed only that one
        const vial_value_t cv = vx_value_of(x, vx_child(x, id, 1));
        if (cv.kind == VIAL_VAL_BOOL) r = vx_ssa(x, vx_child(x, id, cv.b ? 0 : 2));
        else vx_lower_fail(x, id, "a conditional per frame whose branches differ in dtype");
        break;
    }
    default:
        break;
    }
    if (!r && !failed_within(x, vx_span(x, id), x->t != x->c->tree)) {
        char what[64];
        snprintf(what, sizeof(what), "%s is not lowered yet", vial_node_kind_str(k));
        vx_lower_fail(x, id, what);
    }
    vx_set_ssa(x, id, r);
}

vial_ssa_t vx_unpack(infer_t* x, vial_node_id target, vial_ssa_t value, const T* t, uint32_t i, const T* et) {
    if (!x->ir || !value) return 0;
    if (t->dtype == VIAL_DT_TUPLE) return op1(x, target, VIAL_IR_ELEM, 0, i, &value, 1, et);
    // a[i] along the first dimension
    T it = vial_type_scalar(VIAL_DT_INT);
    vial_val_t iv = vial_val_int((int32_t)i, x->ir->alloc);
    const vial_ssa_t c = vial_ir_const(x->ir, &iv, &it, (str_t){0}, vx_span(x, target));
    const uint32_t first = (uint32_t)md_array_size(x->ir->index_items);
    md_array_push(x->ir->index_items, ((vial_index_item_t){.kind = VIAL_INDEX_SCALAR, .a = c}), x->ir->alloc);
    const uint32_t op = vial_ir_emit(x->ir, &(vial_ir_emit_t){.kind = VIAL_IR_INDEX, .sub = 1, .aux = first, .args = (vial_ssa_t[]){value, c},
                                                              .num_args = 2, .types = et, .num_results = 1, .span = vx_span(x, target)});
    return x->ir->ops[op].result;
}

static void target_names(const vial_tree_t* t, vial_node_id target, vial_node_id** out, md_allocator_i* alloc) {
    const vial_node_kind_t k = (vial_node_kind_t)t->nodes[target].kind;
    if (k == VIAL_N_NAME) { md_array_push(*out, target, alloc); return; }
    if (k == VIAL_N_TUPLE || k == VIAL_N_LIST) {
        for (uint32_t i = 0; i < t->nodes[target].child_count; ++i) target_names(t, vial_child(t, target, i), out, alloc);
    }
}

void vx_emit_finish(infer_t* x) {
    vial_ir_t* ir = x->ir;
    ir->descs = *x->descs;
    vial_ir_yield(ir, NULL, 0, (md_vial_span_t){0});
    const vial_tree_t* t = x->c->tree;
    const vial_node_id root = t->root;
    for (uint32_t i = 0; i < t->nodes[root].child_count; ++i) {
        const vial_node_id s = vial_child(t, root, i);
        if (t->nodes[s].kind != VIAL_N_ASSIGN || t->nodes[s].child_count < 2) continue;
        vial_node_id* names = NULL;
        for (uint32_t k = 0; k + 1 < t->nodes[s].child_count; ++k) target_names(t, vial_child(t, s, k), &names, x->alloc);
        for (size_t k = 0; k < md_array_size(names); ++k) {
            const vial_node_id nn = names[k];
            const uint32_t sym = x->n->ref[nn];
            const vial_ssa_t v = sym ? x->top->sym_ssa[sym] : 0;
            const md_vial_span_t sp = t->nodes[nn].span;
            const str_t name = {x->c->src.ptr + sp.beg, sp.end - sp.beg};
            if (!v || !name.len) continue;
            md_array_push(ir->outputs, ((vial_ir_output_t){.name = name, .value = v, .span = sp}), ir->alloc);
        }
    }
}

void vx_emit_return(infer_t* x, vial_node_id value, vial_node_id at) {
    vial_ssa_t r = value ? vx_ssa(x, value) : 0;
    if (!value) {
        const vial_val_t none = vial_val_none();
        const T nt = vial_type_scalar(VIAL_DT_NONE);
        r = vial_ir_const(x->ir, &none, &nt, (str_t){0}, vx_span(x, at));
    }
    if (!r) {
        // The body still ends, so the IR keeps its shape; the function is not called. A value that failed to check
        // was reported already.
        if (!vial_type_is_error(&x->out->types[value - x->out->base])) vx_lower_fail(x, at, "a return whose value is not lowered");
        vial_ir_yield(x->ir, NULL, 0, vx_span(x, at));
        return;
    }
    vial_ir_yield(x->ir, &r, 1, vx_span(x, at));
    vx_close_pending(x);
}

// An argument as the parameter's kind takes it: a set's positions, a set of one particle's position, a point as a
// group of one. Identity otherwise.
static vial_ssa_t to_kind(infer_t* x, vial_node_id at, vial_ssa_t v, const T* vt, uint8_t kind, const T* core) {
    const bool from_set = vt->dtype == VIAL_DT_SET && (kind == VX_KIND_POINT || kind == VX_KIND_POINTS || kind == VX_KIND_GROUP);
    const bool one_point = vt->dtype != VIAL_DT_SET && (kind == VX_KIND_POINTS || kind == VX_KIND_GROUP) && vt->ndim == 1;
    if (!from_set && !one_point) return v;
    return op1(x, at, VIAL_IR_KIND, kind, 0, &v, 1, core);
}

// The call itself: of a function, or applying a nested def's body
static vial_ssa_t call_op(infer_t* x, uint32_t fn, uint32_t body, const vial_ssa_t* args, uint32_t np, const T* t, md_vial_span_t sp) {
    vial_ir_t* ir = x->ir;
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = fn ? VIAL_IR_CALL : VIAL_IR_APPLY, .aux = fn, .args = args, .num_args = np,
                                                      .types = t, .num_results = 1, .span = sp});
    if (!fn) {
        if (ir->regions[body].owner != VIAL_IR_NO_OP) return 0;     // a body is applied once
        vial_ir_attach(ir, op, 0, body);
    }
    return ir->ops[op].result;
}

vial_ssa_t vx_emit_call(infer_t* x, vial_node_id call, uint32_t fn, uint32_t body, const uint8_t* kinds, const uint32_t* argv, const T* atypes,
                        const T* loops, const T* cores, uint32_t np, const T* result) {
    vial_ir_t* ir = x->ir;
    bool mapped = false;
    for (uint32_t k = 0; k < np; ++k) {
        if (!argv[k]) return 0;
        mapped = mapped || loops[k].ndim > 0;
    }
    const md_vial_span_t sp = vx_span(x, call);
    vial_ssa_t conv[32];
    if (!mapped) {
        for (uint32_t k = 0; k < np; ++k) conv[k] = to_kind(x, call, argv[k], &atypes[k], kinds[k], &cores[k]);
        return call_op(x, fn, body, conv, np, result, sp);
    }
    // A map over the dimensions the kinds leave: each call gets one element of each mapped argument
    T elems[32];
    uint32_t mapping[32];
    for (uint32_t k = 0; k < np; ++k) {
        elems[k] = atypes[k];
        mapping[k] = loops[k].ndim;
        if (mapping[k]) {
            elems[k].ndim -= (uint8_t)mapping[k];
            memmove(elems[k].dims, elems[k].dims + mapping[k], elems[k].ndim * sizeof(elems[k].dims[0]));
        }
    }
    const uint32_t region = vial_ir_region(ir, VIAL_IR_NO_OP, ir->regions[ir->insert].func, elems, np);
    const uint32_t saved = ir->insert;
    ir->insert = region;
    for (uint32_t k = 0; k < np; ++k) conv[k] = to_kind(x, call, vial_ir_param(ir, region, k), &elems[k], kinds[k], &cores[k]);
    T one = *result;
    {
        // What one call gives: the result without the mapped dimensions
        uint8_t lnd = 0;
        for (uint32_t k = 0; k < np; ++k) if (loops[k].ndim > lnd) lnd = loops[k].ndim;
        one.ndim -= lnd;
        memmove(one.dims, one.dims + lnd, one.ndim * sizeof(one.dims[0]));
    }
    const vial_ssa_t cr = call_op(x, fn, body, conv, np, &one, sp);
    vial_ir_yield(ir, &cr, cr ? 1 : 0, sp);
    ir->insert = saved;
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_MAP, .args = argv, .num_args = np, .types = result,
                                                      .num_results = 1, .attr = mapping, .num_attr = np, .span = sp});
    vial_ir_attach(ir, op, 0, region);
    return ir->ops[op].result;
}

// ### COMPREHENSIONS ###

static bool callee_named(const infer_t* x, vial_node_id node, const char* name) {
    if (!node || vx_kind(x, node) != VIAL_N_CALL) return false;
    const vial_node_id c = vx_child(x, node, 0);
    if (!c || vx_kind(x, c) != VIAL_N_NAME) return false;
    const vial_symbol_t* s = vial_names_symbol(x->n, c);
    return s && s->kind == VIAL_SYM_BUILTIN && s->builtin && strcmp(s->builtin->name, name) == 0;
}

// What one iterable maps over, emitted in the current region: range() as its ints, a set as its particles, an array
// as itself. Every operand maps over its first dimension.
static vial_ssa_t iter_operand(infer_t* x, vial_node_id it, const T* elem, T* param) {
    *param = *elem;
    if (callee_named(x, it, "range")) {
        vial_ssa_t args[3];
        uint32_t n = 0;
        for (uint32_t i = 1; i < vx_num_children(x, it) && n < 3; ++i) {
            const vial_node_id c = vx_child(x, it, i);
            if (!c) continue;
            args[n] = vx_ssa(x, c);
            if (!args[n]) return 0;
            ++n;
        }
        T t = vial_type_scalar(VIAL_DT_INT);
        t.stage = elem->stage;
        vial_type_push_dim(&t, (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = (uint32_t)it});
        return builtin(x, it, "range", args, NULL, n, &t);
    }
    const vial_ssa_t a = vx_ssa(x, it);
    if (!a) return 0;
    const T* at = type_at(x, it);
    if (at->dtype == VIAL_DT_SET && at->ndim == 0) {
        T t = *elem;
        vial_type_push_dim(&t, (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = (uint32_t)it});
        return op1(x, it, VIAL_IR_KIND, VIAL_KIND_EACH, 0, &a, 1, &t);
    }
    return a;
}

bool vx_comp_begin(infer_t* x, vial_node_id it, vial_node_id target, const T* elem, vx_comp_t* c) {
    memset(c, 0, sizeof(*c));
    T params[8];
    if (callee_named(x, it, "zip")) {
        if (elem->dtype != VIAL_DT_TUPLE) return false;
        uint32_t k = 0;
        for (uint32_t i = 1; i < vx_num_children(x, it); ++i) {
            const vial_node_id a = vx_child(x, it, i);
            if (!a) continue;
            if (k >= 8 || k >= elem->num_elems) return false;
            c->operands[k] = iter_operand(x, a, &elem->elems[k], &params[k]);
            if (!c->operands[k]) return false;
            ++k;
        }
        c->num = k;
    } else {
        c->operands[0] = iter_operand(x, it, elem, &params[0]);
        if (!c->operands[0]) return false;
        c->num = 1;
    }
    vial_ir_t* ir = x->ir;
    c->region = vial_ir_region(ir, VIAL_IR_NO_OP, ir->regions[ir->insert].func, params, c->num);
    c->saved = ir->insert;
    ir->insert = c->region;
    // The target: one name for one operand, or a tuple of them; names in a tuple, one operand each
    if (vx_kind(x, target) == VIAL_N_NAME) {
        vial_ssa_t v = vial_ir_param(ir, c->region, 0);
        if (c->num > 1) {
            vial_ssa_t ps[8];
            for (uint32_t i = 0; i < c->num; ++i) ps[i] = vial_ir_param(ir, c->region, i);
            v = vial_ir_op1(ir, VIAL_IR_TUPLE, 0, 0, ps, c->num, elem, vx_span(x, target));
        }
        const uint32_t sym = x->n->ref[target];
        if (sym) vx_frame_set_ssa(x, sym, v);
        vx_set_ssa(x, target, v);
    } else if (vx_kind(x, target) == VIAL_N_TUPLE || vx_kind(x, target) == VIAL_N_LIST) {
        if (vx_num_children(x, target) != c->num) { ir->insert = c->saved; return false; }
        for (uint32_t i = 0; i < c->num; ++i) {
            const vial_node_id tn = vx_child(x, target, i);
            if (vx_kind(x, tn) != VIAL_N_NAME) { ir->insert = c->saved; return false; }
            const uint32_t sym = x->n->ref[tn];
            if (sym) vx_frame_set_ssa(x, sym, vial_ir_param(ir, c->region, i));
            vx_set_ssa(x, tn, vial_ir_param(ir, c->region, i));
        }
    } else {
        ir->insert = c->saved;
        return false;
    }
    c->ok = true;
    return true;
}

void vx_comp_end(infer_t* x, vial_node_id id, vx_comp_t* c, vial_ssa_t value, vial_ssa_t keep, const T* t) {
    vial_ir_t* ir = x->ir;
    if (!c->ok) return;
    const md_vial_span_t sp = vx_span(x, id);
    if (!value) {
        vial_ir_yield(ir, NULL, 0, sp);
        ir->insert = c->saved;
        return;
    }
    const vial_ssa_t ys[2] = {value, keep};
    vial_ir_yield(ir, ys, keep ? 2 : 1, sp);
    ir->insert = c->saved;
    uint32_t mapping[8];
    for (uint32_t i = 0; i < c->num; ++i) mapping[i] = 1;
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_MAP, .sub = keep ? VIAL_MAP_FILTERED : 0, .args = c->operands,
                                                      .num_args = c->num, .types = t, .num_results = 1, .attr = mapping, .num_attr = c->num, .span = sp});
    vial_ir_attach(ir, op, 0, c->region);
    vx_set_ssa(x, id, ir->ops[op].result);
}

void vx_emit_unreachable(infer_t* x) {
    for (vx_frame_t* f = x->frame; f; f = f->parent) f->lower_failed = true;
    if (x->frame) vial_ir_yield(x->ir, NULL, 0, (md_vial_span_t){0});
}

void vx_emit_unreachable_mark(infer_t* x) {
    for (vx_frame_t* f = x->frame; f; f = f->parent) f->lower_failed = true;
}

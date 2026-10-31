// Passes over the array IR. Inlining: every call of a def's specialization, and every nested def's applied body, is
// replaced by a copy of the body, so the program is one function and what a def computes over the trajectory (a fold,
// a mean over frames) is at the top level, where the trajectory's sweeps see it.

#include "ir.h"
#include "ast.h"
#include "builtins.h"

#include <core/md_allocator.h>
#include <core/md_array.h>

#include <string.h>

typedef struct inl_t {
    const vial_ir_t*    src;
    vial_ir_t*          dst;
    vial_ssa_t*         map;        // by src ssa: its value in dst
    uint32_t            depth;
} inl_t;

static vial_ssa_t remap(const inl_t* c, vial_ssa_t v) { return v ? c->map[v] : 0; }

static bool copy_ops(inl_t* c, uint32_t region, vial_ssa_t* yields, uint32_t max_yields, uint32_t* num_yields);

// A region of an op, copied as a region of the op's copy
static uint32_t copy_region(inl_t* c, uint32_t src_region) {
    const vial_ir_t* src = c->src;
    vial_ir_t* dst = c->dst;
    const vial_ir_region_t* r = &src->regions[src_region];
    vial_type_t params[64];
    for (uint32_t i = 0; i < r->num_params && i < 64; ++i) params[i] = src->values[r->param_beg + i].type;
    const uint32_t nr = vial_ir_region(dst, VIAL_IR_NO_OP, 0, params, r->num_params);
    for (uint32_t i = 0; i < r->num_params; ++i) c->map[r->param_beg + i] = vial_ir_param(dst, nr, i);
    const uint32_t saved = dst->insert;
    dst->insert = nr;
    vial_ssa_t ys[64];
    uint32_t ny = 0;
    copy_ops(c, src_region, ys, 64, &ny);
    vial_ir_yield(dst, ys, ny, (md_vial_span_t){0});
    dst->insert = saved;
    return nr;
}

// A body inlined with its parameters bound to args: what it yields
static bool inline_body(inl_t* c, uint32_t body, const vial_ssa_t* src_args, uint32_t n, vial_ssa_t* result) {
    const vial_ir_region_t* r = &c->src->regions[body];
    if (r->num_params != n || c->depth > 64) return false;
    for (uint32_t i = 0; i < n; ++i) c->map[r->param_beg + i] = remap(c, src_args[i]);
    c->depth += 1;
    vial_ssa_t ys[1];
    uint32_t ny = 0;
    const bool ok = copy_ops(c, body, ys, 1, &ny);
    c->depth -= 1;
    if (!ok || ny != 1) return false;
    *result = ys[0];
    return true;
}

static bool copy_ops(inl_t* c, uint32_t region, vial_ssa_t* yields, uint32_t max_yields, uint32_t* num_yields) {
    const vial_ir_t* src = c->src;
    vial_ir_t* dst = c->dst;
    const vial_ir_region_t* r = &src->regions[region];
    *num_yields = 0;
    for (size_t i = 0; i < md_array_size(r->ops); ++i) {
        const vial_ir_op_t* op = &src->ops[r->ops[i]];
        const vial_ssa_t* args = op->num_args ? src->args + op->arg_beg : NULL;
        if (op->kind == VIAL_IR_YIELD) {
            for (uint32_t k = 0; k < op->num_args && k < max_yields; ++k) yields[k] = remap(c, args[k]);
            *num_yields = op->num_args < max_yields ? op->num_args : max_yields;
            return true;
        }
        if (op->kind == VIAL_IR_CALL || op->kind == VIAL_IR_APPLY) {
            const uint32_t body = op->kind == VIAL_IR_CALL ? src->funcs[op->aux].body : op->region[0];
            vial_ssa_t res;
            if (!inline_body(c, body, args, op->num_args, &res)) return false;
            c->map[op->result] = res;
            continue;
        }
        // shown(x, d) is x: the drawing is for hover, which reads the IR as lowered, so evaluation never makes it
        if (op->kind == VIAL_IR_BUILTIN && op->num_args >= 1 && !strcmp(vial_builtin_at(op->aux)->name, "shown")) {
            c->map[op->result] = remap(c, args[0]);
            continue;
        }
        vial_ssa_t a[256];
        if (op->num_args > 256) return false;
        for (uint32_t k = 0; k < op->num_args; ++k) a[k] = remap(c, args[k]);
        uint32_t aux = op->aux;
        if (op->kind == VIAL_IR_CONST) {
            aux = (uint32_t)md_array_size(dst->consts);
            md_array_push(dst->consts, src->consts[op->aux], dst->alloc);
            md_array_push(dst->const_text, src->const_text[op->aux], dst->alloc);
        } else if (op->kind == VIAL_IR_INDEX) {
            aux = (uint32_t)md_array_size(dst->index_items);
            for (uint32_t k = 0; k < op->sub; ++k) {
                vial_index_item_t it = src->index_items[op->aux + k];
                it.a = remap(c, it.a);
                it.b = remap(c, it.b);
                it.c = remap(c, it.c);
                md_array_push(dst->index_items, it, dst->alloc);
            }
        }
        vial_type_t types[8];
        for (uint32_t k = 0; k < op->num_results && k < 8; ++k) types[k] = src->values[op->result + k].type;
        // Attributes in ints: as many as the kind keeps there
        const uint32_t nattr = op->kind == VIAL_IR_MAP ? op->num_args : op->kind == VIAL_IR_PAIRS ? 2 : 0;
        const uint32_t nop = vial_ir_emit(dst, &(vial_ir_emit_t){.kind = op->kind, .sub = op->sub, .aux = aux, .args = a,
            .names = op->num_args ? src->names + op->arg_beg : NULL, .num_args = op->num_args, .types = types, .num_results = op->num_results,
            .attr = nattr ? src->ints + op->attr : NULL, .num_attr = nattr, .span = op->span});
        for (int k = 0; k < 2; ++k) {
            if (!op->region[k]) continue;
            const uint32_t nr = copy_region(c, op->region[k]);
            vial_ir_attach(dst, nop, (uint32_t)k, nr);
        }
        for (uint32_t k = 0; k < op->num_results; ++k) c->map[op->result + k] = dst->ops[nop].result + k;
    }
    return true;
}

bool vial_ir_inline(const vial_ir_t* src, vial_ir_t* dst, md_allocator_i* alloc) {
    vial_ir_init(dst, alloc);
    dst->descs = src->descs;
    const size_t nv = md_array_size(src->values);
    inl_t c = {.src = src, .dst = dst};
    c.map = md_alloc(alloc, (nv ? nv : 1) * sizeof(vial_ssa_t));
    memset(c.map, 0, nv * sizeof(vial_ssa_t));
    vial_ssa_t ys[1];
    uint32_t ny = 0;
    if (!copy_ops(&c, src->funcs[0].body, ys, 1, &ny)) return false;
    vial_ir_yield(dst, NULL, 0, (md_vial_span_t){0});
    for (size_t i = 0; i < md_array_size(src->outputs); ++i) {
        vial_ir_output_t o = src->outputs[i];
        o.value = remap(&c, o.value);
        md_array_push(dst->outputs, o, alloc);
    }
    for (size_t i = 0; i < md_array_size(src->failures); ++i) md_array_push(dst->failures, src->failures[i], alloc);
    md_free(alloc, c.map, (nv ? nv : 1) * sizeof(vial_ssa_t));
    return true;
}

// ### CONTRACTION ###
// a * b + c becomes fma(a, b, c), rounded once, where the product has no other use and every value is a float of one
// dtype. Decided here, once, so every backend executes exactly these fusions and none of its own.

static void count_uses(const vial_ir_t* ir, uint32_t* uses) {
    for (size_t i = 0; i < md_array_size(ir->ops); ++i) {
        const vial_ir_op_t* op = &ir->ops[i];
        for (uint32_t k = 0; k < op->num_args; ++k) uses[ir->args[op->arg_beg + k]] += 1;
    }
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) uses[ir->outputs[i].value] += 1;
    uses[0] = 0;
}

static bool is_float_dtype(uint8_t dt) { return dt == VIAL_DT_FLOAT || dt == VIAL_DT_FLOAT64; }

uint32_t vial_ir_contract(vial_ir_t* ir, md_allocator_i* alloc) {
    const size_t nv = md_array_size(ir->values);
    uint32_t* uses = md_alloc(alloc, (nv ? nv : 1) * sizeof(uint32_t));
    memset(uses, 0, nv * sizeof(uint32_t));
    count_uses(ir, uses);
    uint32_t fused = 0;
    for (size_t i = 0; i < md_array_size(ir->ops); ++i) {
        vial_ir_op_t* op = &ir->ops[i];
        if (op->kind != VIAL_IR_BINARY || op->sub != VIAL_OP_PLUS || op->num_args != 2) continue;
        const uint8_t dt = ir->values[op->result].type.dtype;
        if (!is_float_dtype(dt)) continue;
        for (int side = 0; side < 2; ++side) {
            const vial_ssa_t p = ir->args[op->arg_beg + side], c = ir->args[op->arg_beg + 1 - side];
            const uint32_t pop = ir->values[p].op;
            if (pop == VIAL_IR_NO_OP || uses[p] != 1) continue;
            const vial_ir_op_t* mul = &ir->ops[pop];
            if (mul->kind != VIAL_IR_BINARY || mul->sub != VIAL_OP_STAR || mul->parent != op->parent) continue;
            const vial_ssa_t a = ir->args[mul->arg_beg], b = ir->args[mul->arg_beg + 1];
            if (ir->values[p].type.dtype != dt || ir->values[a].type.dtype != dt || ir->values[b].type.dtype != dt || ir->values[c].type.dtype != dt) continue;
            // The add becomes the fma, its arguments appended: a, b, c
            op->kind = VIAL_IR_FMA;
            op->sub = 0;
            op->arg_beg = (uint32_t)md_array_size(ir->args);
            op->num_args = 3;
            const vial_ssa_t abc[3] = {a, b, c};
            for (int k = 0; k < 3; ++k) {
                md_array_push(ir->args, abc[k], ir->alloc);
                md_array_push(ir->names, (str_t){0}, ir->alloc);
            }
            uses[p] = 0;
            fused += 1;
            break;
        }
    }
    md_free(alloc, uses, (nv ? nv : 1) * sizeof(uint32_t));
    return fused;
}

// ### REDUCTIONS OVER EVERY PAIR ###
// red(norm(displacement(x, y))) and red(norm(displacement(x, y)) < r), with x of shape [n][1] and y [m] (points or sets
// of one particle each), become a VIAL_IR_PAIRS op of x, y and r. What the reduction reads stays in the IR; with no
// other reader it is not live, and the matrix is never built.

static const char* builtin_name(const vial_ir_t* ir, const vial_ir_op_t* op) {
    return op->kind == VIAL_IR_BUILTIN ? vial_builtin_at(op->aux)->name : "";
}

// The op's argument for a parameter, by name, or 0
static vial_ssa_t arg_named(const vial_ir_t* ir, const vial_ir_op_t* op, const char* name) {
    for (uint32_t k = 0; k < op->num_args; ++k) {
        const str_t nm = ir->names[op->arg_beg + k];
        if (nm.len && str_eq(nm, str_from_cstr(name))) return ir->args[op->arg_beg + k];
    }
    return 0;
}

static uint32_t def_op(const vial_ir_t* ir, vial_ssa_t v) { return v ? ir->values[v].op : VIAL_IR_NO_OP; }

// The shape of x or y's points: a float type's dimensions but its last of 3, or a set type's
static bool lead_shape(const vial_type_t* t, uint32_t* dims, uint8_t* nd) {
    if (t->dtype == VIAL_DT_SET) { *nd = t->ndim; }
    else if (t->dtype == VIAL_DT_FLOAT && t->ndim >= 1 && t->dims[t->ndim - 1].kind == VIAL_EXT_STATIC && t->dims[t->ndim - 1].n == 3) { *nd = t->ndim - 1; }
    else return false;
    for (uint8_t d = 0; d < *nd; ++d) dims[d] = t->dims[d].kind == VIAL_EXT_STATIC ? t->dims[d].n : 0;
    return true;
}

uint32_t vial_ir_pairs(vial_ir_t* ir, md_allocator_i* alloc) {
    (void)alloc;
    uint32_t lowered = 0;
    for (size_t i = 0; i < md_array_size(ir->ops); ++i) {
        vial_ir_op_t* op = &ir->ops[i];
        if (op->kind != VIAL_IR_BUILTIN) continue;
        const char* name = builtin_name(ir, op);
        uint32_t red = 0;
        if (!strcmp(name, "min")) red = VIAL_PAIRS_MIN;
        else if (!strcmp(name, "max")) red = VIAL_PAIRS_MAX;
        else if (!strcmp(name, "argmin")) red = VIAL_PAIRS_ARGMIN;
        else if (!strcmp(name, "argmax")) red = VIAL_PAIRS_ARGMAX;
        else if (!strcmp(name, "sum")) red = VIAL_PAIRS_COUNT;
        else if (!strcmp(name, "any")) red = VIAL_PAIRS_ANY;
        else continue;
        // Only a reduction of a, by axis=: not min(a, b), not by=
        if (arg_named(ir, op, "b") || arg_named(ir, op, "by")) continue;
        const vial_ssa_t a = arg_named(ir, op, "a");
        const vial_ssa_t axis_v = arg_named(ir, op, "axis");
        uint32_t axis = VIAL_PAIRS_ALL;
        if (axis_v) {
            const uint32_t ao = def_op(ir, axis_v);
            if (ao == VIAL_IR_NO_OP || ir->ops[ao].kind != VIAL_IR_CONST) continue;
            const vial_val_t* c = &ir->consts[ir->ops[ao].aux];
            if (c->dtype != VIAL_DT_INT || c->ndim || !c->data) continue;
            const int32_t ax = *(const int32_t*)c->data;
            if (ax == 1 || ax == -1) axis = 1;
            else if (ax == 0 || ax == -2) axis = 0;
            else continue;
        } else if (red == VIAL_PAIRS_ARGMIN || red == VIAL_PAIRS_ARGMAX) {
            continue;   // a flat index over every pair: not one of the rule's
        }
        // The distances, or their test against r
        uint32_t test = 0;
        vial_ssa_t dist = a, r = 0;
        if (red == VIAL_PAIRS_COUNT || red == VIAL_PAIRS_ANY) {
            const uint32_t co = def_op(ir, a);
            if (co == VIAL_IR_NO_OP || ir->ops[co].kind != VIAL_IR_COMPARE) continue;
            const vial_ir_op_t* c = &ir->ops[co];
            if (c->sub != VIAL_CMP_LT && c->sub != VIAL_CMP_LE) continue;
            dist = ir->args[c->arg_beg];
            r = ir->args[c->arg_beg + 1];
            if (!r || ir->values[r].type.ndim != 0) continue;
            test = c->sub;
        }
        const uint32_t no = def_op(ir, dist);
        if (no == VIAL_IR_NO_OP || strcmp(builtin_name(ir, &ir->ops[no]), "norm") || ir->ops[no].num_args < 1) continue;
        const vial_ir_op_t* nop = &ir->ops[no];
        if (ir->values[dist].type.dtype != VIAL_DT_FLOAT || ir->values[dist].type.ndim != 2) continue;
        const vial_ssa_t disp = ir->args[nop->arg_beg];
        const uint32_t dop = def_op(ir, disp);
        if (dop == VIAL_IR_NO_OP || strcmp(builtin_name(ir, &ir->ops[dop]), "displacement")) continue;
        const vial_ir_op_t* dp = &ir->ops[dop];
        const vial_ssa_t x = arg_named(ir, dp, "a") ? arg_named(ir, dp, "a") : ir->args[dp->arg_beg];
        const vial_ssa_t y = arg_named(ir, dp, "b");
        if (!x || !y) continue;
        // x is [n][1] and y [m]
        uint32_t xd[VIAL_MAX_DIMS], yd[VIAL_MAX_DIMS];
        uint8_t xn, yn;
        if (!lead_shape(&ir->values[x].type, xd, &xn) || !lead_shape(&ir->values[y].type, yd, &yn)) continue;
        if (xn != 2 || xd[1] != 1 || yn != 1) continue;
        // x and y are visible to the reduction: displacement read them, and the reduction reads what it gave
        // The reduction becomes the pairs op: its arguments appended, its attributes in ints
        const vial_ssa_t args[3] = {x, y, r};
        op->kind = VIAL_IR_PAIRS;
        op->sub = (uint16_t)red;
        op->aux = 0;
        op->arg_beg = (uint32_t)md_array_size(ir->args);
        op->num_args = 3;
        for (int k = 0; k < 3; ++k) {
            md_array_push(ir->args, args[k], ir->alloc);
            md_array_push(ir->names, (str_t){0}, ir->alloc);
        }
        op->attr = (uint32_t)md_array_size(ir->ints);
        md_array_push(ir->ints, axis, ir->alloc);
        md_array_push(ir->ints, test, ir->alloc);
        lowered += 1;
    }
    return lowered;
}

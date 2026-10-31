// Passes over the array IR. Inlining: every call of a def's specialization, and every nested def's applied body, is
// replaced by a copy of the body, so the program is one function and what a def computes over the trajectory (a fold,
// a mean over frames) is at the top level, where the trajectory's sweeps see it.

#include "ir.h"

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
    vial_ssa_t ys[8];
    uint32_t ny = 0;
    copy_ops(c, src_region, ys, 8, &ny);
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
        const uint32_t nattr = op->kind == VIAL_IR_MAP ? op->num_args : 0;
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

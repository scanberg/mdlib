#include "interp.h"
#include "builtins.h"
#include "ast.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

bool vial_rt_fail(vial_kctx_t* k, const char code[6], const char* fmt, ...) {
    vial_interp_t* in = k->in;
    vial_rt_error_t e = {.span = k->op ? k->op->span : (md_vial_span_t){0}};
    memcpy(e.code, code, 5);
    va_list args;
    va_start(args, fmt);
    vsnprintf(e.message, sizeof(e.message), fmt, args);
    va_end(args);
    if (!in->failed) in->error = e;
    in->failed = true;
    // Once each: the same failure at the same place in every frame is one
    for (size_t i = 0; i < md_array_size(in->errors); ++i) {
        const vial_rt_error_t* g = &in->errors[i];
        if (g->span.beg == e.span.beg && g->span.end == e.span.end && !strcmp(g->message, e.message)) return false;
    }
    md_array_push(in->errors, e, in->alloc);
    return false;
}

bool vial_interp_init(vial_interp_t* in, const vial_ir_t* ir, vial_sys_t* sys, const struct md_system_t* msys, md_allocator_i* alloc) {
    memset(in, 0, sizeof(*in));
    in->ir = ir;
    in->sys = sys;
    in->msys = msys;
    in->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));
    in->temp = md_arena_allocator_create(alloc, MEGABYTES(1));
    if (!in->alloc || !in->temp) return false;
    const size_t nv = md_array_size(ir->values);
    md_array_resize(in->vals, nv, in->alloc);
    md_array_resize(in->state, nv, in->alloc);
    memset(in->vals, 0, nv * sizeof(*in->vals));
    memset(in->state, 0, nv);
    // What the outputs need: the ops that define them, their arguments' ops, and everything in a live op's regions
    const size_t no = md_array_size(ir->ops);
    md_array_resize(in->live, no, in->alloc);
    memset(in->live, 0, no);
    uint32_t* stack = NULL;
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        const uint32_t op = ir->values[ir->outputs[i].value].op;
        if (op != VIAL_IR_NO_OP) md_array_push(stack, op, in->alloc);
    }
    while (md_array_size(stack)) {
        const uint32_t op = *md_array_last(stack);
        md_array_pop(stack);
        if (in->live[op]) continue;
        in->live[op] = 1;
        const vial_ir_op_t* o = &ir->ops[op];
        for (uint32_t k = 0; k < o->num_args; ++k) {
            const vial_ssa_t a = ir->args[o->arg_beg + k];
            if (a && ir->values[a].op != VIAL_IR_NO_OP) md_array_push(stack, ir->values[a].op, in->alloc);
        }
        if (o->kind == VIAL_IR_CALL) {
            const vial_ir_region_t* r = &ir->regions[ir->funcs[o->aux].body];
            for (size_t q = 0; q < md_array_size(r->ops); ++q) md_array_push(stack, r->ops[q], in->alloc);
        }
        for (int k = 0; k < 2; ++k) {
            if (!o->region[k]) continue;
            const vial_ir_region_t* r = &ir->regions[o->region[k]];
            for (size_t q = 0; q < md_array_size(r->ops); ++q) md_array_push(stack, r->ops[q], in->alloc);
        }
    }
    return true;
}

// at(x, frame=i): x as it is at frame i, computed in a step of its own before frames are evaluated
static bool eval_at(vial_interp_t* in, uint32_t opi) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_op_t* op = &ir->ops[opi];
    vial_kctx_t k = {.in = in, .op = op};
    const vial_ssa_t xv = op->num_args > 0 ? ir->args[op->arg_beg] : 0;
    const vial_ssa_t fv = op->num_args > 1 ? ir->args[op->arg_beg + 1] : 0;
    if (!xv || !fv || !in->state[fv]) return vial_rt_fail(&k, "R0103", "at() needs a value and a frame");
    const vial_val_t* f = &in->vals[fv];
    if (f->dtype != VIAL_DT_INT || !f->data) return vial_rt_fail(&k, "R0103", "at()'s frame is an int");
    const int64_t index = *(const int32_t*)f->data;
    if (!in->load_frame) return vial_rt_fail(&k, "R0103", "at() needs the run's frames");
    vial_frame_t fr;
    if (!in->load_frame(in->load_user, index, &fr)) return vial_rt_fail(&k, "R0104", "at() reads frame %lld, which the run does not have", (long long)index);
    vial_interp_t child;
    if (!vial_interp_init(&child, ir, in->sys, in->msys, md_get_heap_allocator())) return vial_rt_fail(&k, "R0103", "out of memory");
    // The child computes only what the value needs: the live ops are those of x
    memset(child.live, 0, md_array_size(child.live));
    {
        uint32_t* stack = NULL;
        md_array_push(stack, ir->values[xv].op, child.alloc);
        while (md_array_size(stack)) {
            const uint32_t q = *md_array_last(stack);
            md_array_pop(stack);
            if (q == VIAL_IR_NO_OP || child.live[q]) continue;
            child.live[q] = 1;
            const vial_ir_op_t* o = &ir->ops[q];
            for (uint32_t a = 0; a < o->num_args; ++a) {
                const vial_ssa_t s = ir->args[o->arg_beg + a];
                if (s) md_array_push(stack, ir->values[s].op, child.alloc);
            }
            if (o->kind == VIAL_IR_CALL) {
                const vial_ir_region_t* r = &ir->regions[ir->funcs[o->aux].body];
                for (size_t z = 0; z < md_array_size(r->ops); ++z) md_array_push(stack, r->ops[z], child.alloc);
            }
            for (int z = 0; z < 2; ++z) {
                if (!o->region[z]) continue;
                const vial_ir_region_t* r = &ir->regions[o->region[z]];
                for (size_t w = 0; w < md_array_size(r->ops); ++w) md_array_push(stack, r->ops[w], child.alloc);
            }
        }
    }
    child.load_frame = in->load_frame;
    child.load_user = in->load_user;
    const bool ok = vial_interp_frame(&child, &fr);
    const vial_val_t* v = ok ? vial_interp_value(&child, xv) : NULL;
    bool done = false;
    if (v) {
        in->vals[op->result] = vial_val_copy(v, in->alloc);
        done = true;
    } else {
        vial_rt_fail(&k, child.failed ? child.error.code : "R0103", "at(): %s", child.failed ? child.error.message : "the value is not computed");
    }
    vial_interp_free(&child);
    return done;
}

void vial_interp_free(vial_interp_t* in) {
    if (in->temp) md_arena_allocator_destroy(in->temp);
    if (in->alloc) md_arena_allocator_destroy(in->alloc);
    memset(in, 0, sizeof(*in));
}

const vial_val_t* vial_interp_value(const vial_interp_t* in, vial_ssa_t v) {
    return v && v < md_array_size(in->state) && in->state[v] ? &in->vals[v] : NULL;
}

// ### REGIONS ###

static bool eval_region(vial_interp_t* in, uint32_t region, const vial_val_t* params, vial_val_t* yields, uint32_t max_yields,
                        uint32_t* num_yields, md_allocator_i* alloc);

// One element of a value along its leading dimensions: a view of numbers, or the picked set
static vial_val_t element(const vial_val_t* a, uint32_t lead, size_t i, md_allocator_i* alloc) {
    vial_val_t e = *a;
    e.ndim = a->ndim - (uint8_t)lead;
    memmove(e.dims, a->dims + lead, e.ndim * sizeof(uint32_t));
    size_t inner = 1;
    for (uint8_t k = 0; k < e.ndim; ++k) inner *= e.dims[k];
    if (a->dtype == VIAL_DT_SET) {
        if (e.ndim == 0) {
            uint32_t pos = (uint32_t)i;
            vial_sets_t* s = vial_sets_pick(a->sets, &pos, 1, alloc);
            s->single = true;
            e.sets = s;
        } else {
            uint32_t* pos = md_alloc(alloc, (inner ? inner : 1) * sizeof(uint32_t));
            for (size_t k = 0; k < inner; ++k) pos[k] = (uint32_t)(i * inner + k);
            e.sets = vial_sets_pick(a->sets, pos, (uint32_t)inner, alloc);
        }
    } else if (a->data) {
        e.data = (uint8_t*)a->data + i * inner * vial_dtype_size(a->dtype);
    }
    e.missing = NULL;
    if (a->missing) {
        for (size_t k = 0; k < inner; ++k) if (vial_val_is_missing(a, i * inner + k)) vial_val_set_missing(&e, k, alloc);
    }
    return e;
}

// Writes element i of a stacked result from one call's value
static bool place(vial_kctx_t* k, vial_val_t* out, size_t i, const vial_val_t* e, uint32_t** set_pos, const vial_val_t** set_src) {
    size_t inner = 1;
    for (uint8_t q = out->ndim - e->ndim; q < out->ndim; ++q) inner *= out->dims[q];
    if (vial_val_count(e) != inner) return vial_rt_fail(k, "R0103", "the calls of a map give values of different shapes");
    if (out->dtype == VIAL_DT_SET) {
        // Sets are gathered at the end: remember which element of which call's sets
        (void)set_pos; (void)set_src;
        return true;
    }
    if (e->dtype != out->dtype) {
        const vial_val_t c = vial_val_as(e, out->dtype, k->alloc);
        return place(k, out, i, &c, set_pos, set_src);
    }
    const size_t es = vial_dtype_size(out->dtype);
    if (es && e->data) memcpy((uint8_t*)out->data + i * inner * es, e->data, inner * es);
    if (e->missing) {
        for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(e, q)) vial_val_set_missing(out, i * inner + q, k->alloc);
    }
    return true;
}

static bool eval_map(vial_interp_t* in, const vial_ir_op_t* op, const vial_val_t* const* av, vial_val_t* out, md_allocator_i* alloc) {
    vial_kctx_t k = {.in = in, .op = op, .type = &in->ir->values[op->result].type, .alloc = alloc};
    const uint32_t* mapping = in->ir->ints + op->attr;
    // The mapped dimensions of every operand, broadcast as NumPy's are, aligned on the right
    uint8_t lnd = 0;
    uint32_t ldims[VIAL_MAX_DIMS] = {0};
    for (uint32_t i = 0; i < op->num_args; ++i) if (mapping[i] > lnd) lnd = (uint8_t)mapping[i];
    for (uint8_t d = 0; d < lnd; ++d) ldims[d] = 1;
    for (uint32_t i = 0; i < op->num_args; ++i) {
        const uint32_t m = mapping[i];
        for (uint32_t d = 0; d < m; ++d) {
            const uint32_t e = av[i]->dims[d];
            uint32_t* dst = &ldims[lnd - m + d];
            if (*dst == 1) *dst = e;
            else if (e != 1 && e != *dst) return vial_rt_fail(&k, "R0103", "a map over extents %u and %u, which do not broadcast", *dst, e);
        }
    }
    size_t count = 1;
    for (uint8_t d = 0; d < lnd; ++d) count *= ldims[d];
    const vial_ir_region_t* r = &in->ir->regions[op->region[0]];
    vial_val_t* params = md_alloc(alloc, (r->num_params ? r->num_params : 1) * sizeof(vial_val_t));
    const bool filtered = op->sub == VIAL_MAP_FILTERED;
    const vial_type_t* rt = k.type;
    vial_val_t* results = md_alloc(alloc, (count ? count : 1) * sizeof(vial_val_t));
    uint8_t* keep = filtered ? md_alloc(alloc, count ? count : 1) : NULL;
    size_t kept = 0;
    for (size_t i = 0; i < count; ++i) {
        // The element of each operand at loop position i, by its own broadcast position
        size_t rem = i;
        uint32_t idx[VIAL_MAX_DIMS];
        for (int d = lnd - 1; d >= 0; --d) { idx[d] = (uint32_t)(rem % ldims[d]); rem /= ldims[d]; }
        for (uint32_t p = 0; p < op->num_args; ++p) {
            const uint32_t m = mapping[p];
            if (!m) { params[p] = *av[p]; continue; }
            size_t lin = 0;
            for (uint32_t d = 0; d < m; ++d) {
                const uint32_t e = av[p]->dims[d];
                const uint32_t j = e == 1 ? 0 : idx[lnd - m + d];
                lin = lin * e + j;
            }
            params[p] = element(av[p], m, lin, alloc);
        }
        vial_val_t ys[2];
        uint32_t ny = 0;
        if (!eval_region(in, op->region[0], params, ys, 2, &ny, alloc)) return false;
        if (ny < 1) return vial_rt_fail(&k, "R0103", "a map's call gave no value");
        if (filtered) {
            if (ny < 2 || ys[1].dtype != VIAL_DT_BOOL || !ys[1].data) return vial_rt_fail(&k, "R0103", "a filter gave no bool");
            keep[i] = *(const uint8_t*)ys[1].data && !vial_val_is_missing(&ys[1], 0);
            if (!keep[i]) continue;
            kept += 1;
        }
        results[i] = ys[0];
    }
    // Stacked: the loop's dimensions (one, when filtered), then one call's
    vial_val_t o = {.dtype = rt->dtype};
    if (filtered) { o.ndim = 1; o.dims[0] = (uint32_t)kept; }
    else { o.ndim = lnd; memcpy(o.dims, ldims, lnd * sizeof(uint32_t)); }
    const vial_val_t* first = NULL;
    for (size_t i = 0; i < count && !first; ++i) if (!filtered || keep[i]) first = &results[i];
    vial_val_t inner_shape = first ? *first : (vial_val_t){0};
    if (!first) {
        // No call: the element's shape is the type's, where it is static
        const uint8_t ond = filtered ? 1 : lnd;
        inner_shape.ndim = rt->ndim > ond ? rt->ndim - ond : 0;
        for (uint8_t q = 0; q < inner_shape.ndim; ++q) inner_shape.dims[q] = rt->dims[ond + q].kind == VIAL_EXT_STATIC ? rt->dims[ond + q].n : 0;
    }
    for (uint8_t q = 0; q < inner_shape.ndim; ++q) o.dims[o.ndim++] = inner_shape.dims[q];
    if (o.dtype == VIAL_DT_SET) {
        // An array of sets: each call's sets, in order
        uint32_t total = 0, nsets = 0;
        for (size_t i = 0; i < count; ++i) {
            if (filtered && !keep[i]) continue;
            const vial_sets_t* s = results[i].sets;
            nsets += s->count;
            total += s->offset[s->count];
        }
        vial_sets_t* s = md_alloc(alloc, sizeof(vial_sets_t));
        memset(s, 0, sizeof(*s));
        s->count = nsets;
        s->offset = md_alloc(alloc, (nsets + 1) * sizeof(uint32_t));
        s->idx = md_alloc(alloc, (total ? total : 1) * sizeof(uint32_t));
        uint32_t at = 0, n = 0;
        s->offset[0] = 0;
        for (size_t i = 0; i < count; ++i) {
            if (filtered && !keep[i]) continue;
            const vial_sets_t* e = results[i].sets;
            for (uint32_t q = 0; q < e->count; ++q) {
                const uint32_t sz = vial_sets_size(e, q);
                memcpy(s->idx + at, vial_sets_items(e, q), sz * sizeof(uint32_t));
                at += sz;
                s->offset[++n] = at;
            }
        }
        if (o.ndim == 0) s->single = true;
        o.sets = s;
        *out = o;
        return true;
    }
    if (o.dtype == VIAL_DT_TUPLE || o.dtype == VIAL_DT_NONE || o.dtype == VIAL_DT_DRAWING) {
        *out = (vial_val_t){.dtype = o.dtype, .ndim = o.ndim};
        memcpy(out->dims, o.dims, sizeof(o.dims));
        return true;
    }
    const vial_val_t made = vial_val_make(o.dtype, o.ndim, o.dims, alloc);
    o.data = made.data;
    size_t j = 0;
    for (size_t i = 0; i < count; ++i) {
        if (filtered && !keep[i]) continue;
        if (!place(&k, &o, j++, &results[i], NULL, NULL)) return false;
    }
    *out = o;
    return true;
}

static bool eval_op(vial_interp_t* in, uint32_t opi, md_allocator_i* alloc) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_op_t* op = &ir->ops[opi];
    vial_kctx_t k = {.in = in, .op = op, .type = op->num_results ? &ir->values[op->result].type : NULL, .alloc = alloc};
    const vial_val_t* av[64];
    if (op->num_args > 64) return vial_rt_fail(&k, "R0103", "an op with more than 64 arguments");
    for (uint32_t i = 0; i < op->num_args; ++i) {
        const vial_ssa_t a = ir->args[op->arg_beg + i];
        av[i] = a ? &in->vals[a] : NULL;
    }
    vial_val_t* out = op->num_results ? &in->vals[op->result] : NULL;
    switch (op->kind) {
    case VIAL_IR_CONST:   *out = ir->consts[op->aux]; return true;
    case VIAL_IR_UNARY:   return vial_k_unary(&k, op->sub, av[0], out);
    case VIAL_IR_NOT:     return vial_k_not(&k, av[0], out);
    case VIAL_IR_BINARY:  return vial_k_binary(&k, op->sub, av[0], av[1], out);
    case VIAL_IR_FMA:     return vial_k_fma(&k, av[0], av[1], av[2], out);
    case VIAL_IR_COMPARE: return vial_k_compare(&k, op->sub, av[0], av[1], out);
    case VIAL_IR_LOGIC:   return vial_k_logic(&k, op->sub, av[0], av[1], out);
    case VIAL_IR_SETOP:   return vial_k_setop(&k, op->sub, av[0], av[1], out);
    case VIAL_IR_CONVERT: return vial_k_convert(&k, av[0], k.type->dtype, out);
    case VIAL_IR_KIND:    return vial_k_kind(&k, op->sub, av[0], out);
    case VIAL_IR_INDEX:   return vial_k_index(&k, av[0], ir->index_items + op->aux, op->sub, out);
    case VIAL_IR_LIST: {
        vial_val_t items[64];
        for (uint32_t i = 0; i < op->num_args; ++i) items[i] = *av[i];
        return vial_k_stack(&k, items, op->num_args, out);
    }
    case VIAL_IR_TUPLE: {
        vial_val_t* el = md_alloc(alloc, (op->num_args ? op->num_args : 1) * sizeof(vial_val_t));
        for (uint32_t i = 0; i < op->num_args; ++i) el[i] = *av[i];
        *out = (vial_val_t){.dtype = VIAL_DT_TUPLE, .elems = el, .num_elems = op->num_args};
        return true;
    }
    case VIAL_IR_ELEM:
        if (av[0]->dtype != VIAL_DT_TUPLE || op->aux >= av[0]->num_elems) return vial_rt_fail(&k, "R0103", "an element of what is not a tuple");
        *out = av[0]->elems[op->aux];
        return true;
    case VIAL_IR_BUILTIN: {
        const vial_builtin_t* b = vial_builtin_at(op->aux);
        if (!b) return vial_rt_fail(&k, "R0103", "an unknown built-in");
        return vial_k_builtin(&k, b->name, av, ir->names + op->arg_beg, op->num_args, out);
    }
    case VIAL_IR_CALL: {
        const uint32_t body = ir->funcs[op->aux].body;
        vial_val_t params[32];
        for (uint32_t i = 0; i < op->num_args && i < 32; ++i) params[i] = *av[i];
        uint32_t ny = 0;
        if (!eval_region(in, body, params, out, 1, &ny, alloc)) return false;
        if (ny != 1) return vial_rt_fail(&k, "R0103", "a def gave no value");
        return true;
    }
    case VIAL_IR_APPLY: {
        vial_val_t params[32];
        for (uint32_t i = 0; i < op->num_args && i < 32; ++i) params[i] = *av[i];
        uint32_t ny = 0;
        if (!eval_region(in, op->region[0], params, out, 1, &ny, alloc)) return false;
        if (ny != 1) return vial_rt_fail(&k, "R0103", "a def gave no value");
        return true;
    }
    case VIAL_IR_MAP:
        return eval_map(in, op, av, out, alloc);
    case VIAL_IR_IF: {
        if (av[0]->dtype != VIAL_DT_BOOL || !av[0]->data) return vial_rt_fail(&k, "R0103", "an if without a bool");
        const bool c = *(const uint8_t*)av[0]->data;
        uint32_t ny = 0;
        return eval_region(in, op->region[c ? 0 : 1], NULL, out, op->num_results, &ny, alloc);
    }
    case VIAL_IR_FOLD:
        return vial_rt_fail(&k, "R0103", "a fold is evaluated over the trajectory, not per frame");
    default:
        return vial_rt_fail(&k, "R0103", "%s is not evaluated yet", vial_ir_kind_str((vial_ir_kind_t)op->kind));
    }
}

static bool eval_region(vial_interp_t* in, uint32_t region, const vial_val_t* params, vial_val_t* yields, uint32_t max_yields,
                        uint32_t* num_yields, md_allocator_i* alloc) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_region_t* r = &ir->regions[region];
    for (uint32_t i = 0; i < r->num_params; ++i) {
        in->vals[r->param_beg + i] = params[i];
        in->state[r->param_beg + i] = 1;
    }
    *num_yields = 0;
    for (size_t i = 0; i < md_array_size(r->ops); ++i) {
        const uint32_t opi = r->ops[i];
        const vial_ir_op_t* op = &ir->ops[opi];
        if (op->kind == VIAL_IR_YIELD) {
            for (uint32_t k = 0; k < op->num_args && k < max_yields; ++k) yields[k] = in->vals[ir->args[op->arg_beg + k]];
            *num_yields = op->num_args;
            return true;
        }
        if (!eval_op(in, opi, alloc)) return false;
        for (uint32_t k = 0; k < op->num_results; ++k) in->state[op->result + k] = 1;
    }
    return true;
}

// ### THE TOP LEVEL ###

static bool is_static(const vial_ir_t* ir, const vial_ir_op_t* op) {
    return op->num_results && ir->values[op->result].type.stage == VIAL_STAGE_STATIC;
}

static bool is_trajectory(const vial_ir_t* ir, const vial_ir_op_t* op) {
    return op->num_results && ir->values[op->result].type.stage == VIAL_STAGE_TRAJECTORY;
}

static bool eval_top(vial_interp_t* in, bool static_only) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_region_t* r = &ir->regions[ir->funcs[0].body];
    bool ok = true;
    for (size_t i = 0; i < md_array_size(r->ops); ++i) {
        const uint32_t opi = r->ops[i];
        const vial_ir_op_t* op = &ir->ops[opi];
        if (op->kind == VIAL_IR_YIELD) break;
        if (!op->num_results || !in->live[opi]) continue;
        const bool st = is_static(ir, op);
        if (st && in->state[op->result] == 2) continue;
        if (static_only && !st) continue;
        if (!st && is_trajectory(ir, op)) continue;
        if (st && op->kind == VIAL_IR_BUILTIN && !strcmp(vial_builtin_at(op->aux)->name, "at")) {
            if (eval_at(in, opi)) in->state[op->result] = 2;
            else ok = false;
            continue;
        }
        // A value computed from one that failed is not computed either
        bool ready = true;
        for (uint32_t k = 0; k < op->num_args; ++k) {
            const vial_ssa_t a = ir->args[op->arg_beg + k];
            if (a && !in->state[a]) ready = false;
        }
        if (!ready) { ok = false; continue; }
        if (!eval_op(in, opi, st ? in->alloc : in->temp)) { ok = false; continue; }
        for (uint32_t k = 0; k < op->num_results; ++k) in->state[op->result + k] = st ? 2 : 1;
    }
    return ok;
}

bool vial_interp_static(vial_interp_t* in) {
    in->frame = NULL;
    return eval_top(in, true);
}

bool vial_interp_frame(vial_interp_t* in, const vial_frame_t* frame) {
    // What the last frame computed goes
    const size_t nv = md_array_size(in->state);
    for (size_t i = 0; i < nv; ++i) if (in->state[i] == 1) in->state[i] = 0;
    md_arena_allocator_reset(in->temp);
    const bool st = eval_top(in, true);
    in->frame = frame;
    return eval_top(in, false) && st;
}

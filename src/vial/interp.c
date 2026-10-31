// No contraction: where a multiply and an add fuse is the IR's decision (VIAL_IR_FMA). GCC and Clang get
// -ffp-contract=off for this file from CMake; MSVC's /fp:contract is taken back here, which it allows per file.
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fp_contract(off)
#endif

#include "interp.h"
#include "builtins.h"
#include "ast.h"
#include "infer_internal.h"

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_intrinsics.h>

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

bool vial_rt_fail(vial_kctx_t* k, const char code[6], const char* fmt, ...) {
    vial_interp_t* in = k->in;
    vial_rt_error_t e = {.span = k->op ? k->op->span : (md_vial_span_t){0}, .sweep = in->sweep,
                         .frame = in->frame ? in->frame->index : INT64_MAX, .seq = in->err_seq++};
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

bool vial_k_room(vial_kctx_t* k, size_t count) {
    const size_t limit = k->in->max_elements ? k->in->max_elements : VIAL_MAX_ELEMENTS;
    if (count <= limit) return true;
    return vial_rt_fail(k, "R0105", "a value of %zu elements; the limit is %zu", count, limit);
}

bool vial_k_make(vial_kctx_t* k, uint8_t dtype, uint8_t ndim, const uint32_t* dims, vial_val_t* out) {
    // The count saturates rather than wrapping: six extents of 2^32 overflow 64 bits
    size_t n = 1;
    for (uint8_t i = 0; i < ndim; ++i) {
        if (dims[i] && n > SIZE_MAX / dims[i]) { n = SIZE_MAX; break; }
        n *= dims[i];
    }
    if (!vial_k_room(k, n)) return false;
    // Not zeroed: every kernel writes every element it makes, and what it marks missing is zeroed by the mask. A
    // debug build fills it with ones bits instead (NaN, -1, True), so a kernel that reads what it did not write
    // shows in the tests' values.
    *out = vial_val_make_uninit(dtype, ndim, dims, k->alloc);
#ifndef NDEBUG
    if (out->data) memset(out->data, 0xFF, vial_val_count(out) * vial_dtype_size(dtype));
#endif
    return true;
}

bool vial_k_make_zero(vial_kctx_t* k, uint8_t dtype, uint8_t ndim, const uint32_t* dims, vial_val_t* out) {
    if (!vial_k_make(k, dtype, ndim, dims, out)) return false;
    if (out->data) memset(out->data, 0, vial_val_count(out) * vial_dtype_size(dtype));
    return true;
}

// The live ops: those that define the roots, their arguments' ops, and in a live op's regions, the yield and what
// it needs. Every op is pure, so an op nothing live reads is not evaluated, in a region as at the top: a map's calls
// skip the witnesses of a built-in that only a drawing, dropped when inlining, would have read.
void vial_interp_only(vial_interp_t* in, const vial_ssa_t* roots, size_t n) {
    const vial_ir_t* ir = in->ir;
    memset(in->live, 0, md_array_size(in->live));
    uint32_t* stack = NULL;
    for (size_t i = 0; i < n; ++i) if (roots[i] && ir->values[roots[i]].op != VIAL_IR_NO_OP) md_array_push(stack, ir->values[roots[i]].op, in->alloc);
    const uint32_t* yields[3];
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
        uint32_t nr = 0;
        if (o->kind == VIAL_IR_CALL) yields[nr++] = &ir->funcs[o->aux].body;
        for (int k = 0; k < 2; ++k) if (o->region[k]) yields[nr++] = &o->region[k];
        for (uint32_t k = 0; k < nr; ++k) {
            const vial_ir_region_t* r = &ir->regions[*yields[k]];
            // A region ends in its yield
            const size_t m = md_array_size(r->ops);
            if (m) md_array_push(stack, r->ops[m - 1], in->alloc);
        }
    }
}

bool vial_interp_init(vial_interp_t* in, const vial_ir_t* ir, vial_sys_t* sys, const struct md_system_t* msys, md_allocator_i* alloc) {
    memset(in, 0, sizeof(*in));
    in->ir = ir;
    in->sys = sys;
    in->msys = msys;
    in->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));
    // A frame's values: a virtual memory arena, reset per frame and keeping what it committed, so a frame allocates
    // no pages and faults none in after the first
    in->temp = md_vm_arena_create(GIGABYTES(4));
    if (!in->alloc || !in->temp) return false;
    const size_t nv = md_array_size(ir->values);
    md_array_resize(in->vals, nv, in->alloc);
    md_array_resize(in->state, nv, in->alloc);
    memset(in->vals, 0, nv * sizeof(*in->vals));
    memset(in->state, 0, nv);
    // What the outputs need
    const size_t no = md_array_size(ir->ops);
    md_array_resize(in->live, no, in->alloc);
    md_array_resize(in->cell_ext, no, in->alloc);
    if (no) memset(in->cell_ext, 0, no * sizeof(double));
    uint32_t* roots = NULL;
    // What the outputs need, but a hidden name's (_x): nothing reads it but hover, which evaluates for itself, so a
    // matrix bound to _dm and only reduced is never built
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        const str_t nm = ir->outputs[i].name;
        if (nm.len && nm.ptr[0] == '_') continue;
        md_array_push(roots, ir->outputs[i].value, in->alloc);
    }
    vial_interp_only(in, roots, md_array_size(roots));
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
    vial_interp_only(&child, &xv, 1);
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
    for (int d = 0; d < VIAL_MAP_SCRATCH_DEPTH; ++d) if (in->map_scratch[d]) md_vm_arena_destroy(in->map_scratch[d]);
    if (in->temp) md_vm_arena_destroy(in->temp);
    if (in->alloc) md_arena_allocator_destroy(in->alloc);
    memset(in, 0, sizeof(*in));
}

const vial_val_t* vial_interp_value(const vial_interp_t* in, vial_ssa_t v) {
    return v && v < md_array_size(in->state) && in->state[v] ? &in->vals[v] : NULL;
}

// ### REGIONS ###

static bool eval_region(vial_interp_t* in, uint32_t region, const vial_val_t* params, vial_val_t* yields, uint32_t max_yields,
                        uint32_t* num_yields, md_allocator_i* alloc);
static bool eval_op(vial_interp_t* in, uint32_t opi, md_allocator_i* alloc);
static bool eval_region_pub(vial_interp_t* in, uint32_t region, const vial_val_t* params, vial_val_t* yields, uint32_t max_yields,
                            uint32_t* num_yields, md_allocator_i* alloc) {
    return eval_region(in, region, params, yields, max_yields, num_yields, alloc);
}
static bool eval_op_pub(vial_interp_t* in, uint32_t opi, md_allocator_i* alloc) { return eval_op(in, opi, alloc); }

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

// Writes element i of a stacked result from one call's value; a conversion of its dtype is made in scratch
static bool place(vial_kctx_t* k, vial_val_t* out, size_t i, const vial_val_t* e, md_allocator_i* scratch) {
    size_t inner = 1;
    for (uint8_t q = out->ndim - e->ndim; q < out->ndim; ++q) inner *= out->dims[q];
    if (vial_val_count(e) != inner) return vial_rt_fail(k, "R0103", "the calls of a map give values of different shapes");
    if (e->dtype != out->dtype) {
        const vial_val_t c = vial_val_as(e, out->dtype, scratch);
        if (c.dtype != out->dtype) return vial_rt_fail(k, "R0103", "the calls of a map give values of different kinds");
        return place(k, out, i, &c, scratch);
    }
    const size_t es = vial_dtype_size(out->dtype);
    if (es && e->data) memcpy((uint8_t*)out->data + i * inner * es, e->data, inner * es);
    if (e->missing) {
        for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(e, q)) vial_val_set_missing(out, i * inner + q, k->alloc);
    }
    return true;
}

// The scratch arena of a map at the current depth, made on first use. A call of a map is evaluated in it and it is
// reset once the call's value is placed, so a map holds its result and one call's intermediates, not every call's.
// A map in a call is one deeper and has its own, so resetting it leaves the outer call's values alone.
static md_allocator_i* map_scratch(vial_interp_t* in) {
    if (in->map_depth >= VIAL_MAP_SCRATCH_DEPTH) return NULL;
    md_allocator_i** a = &in->map_scratch[in->map_depth];
    if (!*a) *a = md_vm_arena_create(GIGABYTES(4));
    return *a;
}

// The parameter of map region r that value v is, or -1
static int param_of(const vial_ir_t* ir, uint32_t r, vial_ssa_t v) {
    if (!v || ir->values[v].op != VIAL_IR_NO_OP || ir->values[v].region != r) return -1;
    return (int)ir->values[v].index;
}

// A map whose calls are only the distance of closest_pair or farthest_pair between two of its parameters, each one set
// taken as a group, as _between's map over arrays of groups is once inlined: the built-in gets every call at once
// (vial_k_pair_batch), which gives the same bits as the calls one by one without a call's overhead each. -1 when the
// map is not one; otherwise whether it succeeded.
static int eval_map_pairs(vial_interp_t* in, const vial_ir_op_t* op, const vial_val_t* const* av, const uint32_t* mapping,
                          const uint32_t* ldims, uint8_t lnd, size_t count, vial_val_t* out, md_allocator_i* alloc) {
    const vial_ir_t* ir = in->ir;
    if (in->plain || op->sub == VIAL_MAP_FILTERED || in->ir->values[op->result].type.dtype != VIAL_DT_FLOAT) return -1;
    const uint32_t ri = op->region[0];
    const vial_ir_region_t* r = &ir->regions[ri];
    // Its live ops: a kind of each parameter it reads, the built-in, the distance out of its tuple, the yield
    uint32_t kinds[2], nk = 0, bop = VIAL_IR_NO_OP, eop = VIAL_IR_NO_OP, yop = VIAL_IR_NO_OP;
    for (size_t i = 0; i < md_array_size(r->ops); ++i) {
        const uint32_t opi = r->ops[i];
        const vial_ir_op_t* o = &ir->ops[opi];
        if (o->kind != VIAL_IR_YIELD && !in->live[opi]) continue;
        switch (o->kind) {
        case VIAL_IR_KIND:    if (nk == 2) return -1; kinds[nk++] = opi; break;
        case VIAL_IR_BUILTIN: if (bop != VIAL_IR_NO_OP) return -1; bop = opi; break;
        case VIAL_IR_ELEM:    if (eop != VIAL_IR_NO_OP) return -1; eop = opi; break;
        case VIAL_IR_YIELD:   yop = opi; break;
        default:              return -1;
        }
        if (yop != VIAL_IR_NO_OP) break;
    }
    if (bop == VIAL_IR_NO_OP || eop == VIAL_IR_NO_OP || yop == VIAL_IR_NO_OP) return -1;
    const vial_ir_op_t* b = &ir->ops[bop];
    const char* bname = vial_builtin_at(b->aux)->name;
    const bool closest = !strcmp(bname, "closest_pair");
    if (!closest && strcmp(bname, "farthest_pair")) return -1;
    const vial_ir_op_t* e = &ir->ops[eop];
    const vial_ir_op_t* y = &ir->ops[yop];
    if (e->aux != 0 || e->num_args != 1 || ir->args[e->arg_beg] != b->result) return -1;
    if (y->num_args != 1 || ir->args[y->arg_beg] != e->result) return -1;
    // p and q: kinds as groups or points of mapped parameters, each element one set; nothing else given
    int par[2] = {-1, -1};
    for (uint32_t a = 0; a < b->num_args; ++a) {
        const vial_ssa_t v = ir->args[b->arg_beg + a];
        if (!v) continue;
        const str_t nm = ir->names[b->arg_beg + a];
        const int which = str_eq(nm, STR_LIT("p")) ? 0 : str_eq(nm, STR_LIT("q")) ? 1 : -1;
        if (which < 0) return -1;
        const uint32_t ko = ir->values[v].op;
        if (ko == VIAL_IR_NO_OP || ir->ops[ko].kind != VIAL_IR_KIND || (ko != kinds[0] && (nk < 2 || ko != kinds[1]))) return -1;
        if (ir->ops[ko].sub != VX_KIND_GROUP && ir->ops[ko].sub != VX_KIND_POINTS) return -1;
        par[which] = param_of(ir, ri, ir->args[ir->ops[ko].arg_beg]);
    }
    for (int w = 0; w < 2; ++w) {
        if (par[w] < 0 || (uint32_t)par[w] >= op->num_args) return -1;
        const vial_val_t* g = av[par[w]];
        if (g->dtype != VIAL_DT_SET || !g->sets || !mapping[par[w]] || g->ndim != mapping[par[w]] || g->sets->count != vial_val_count(g)) return -1;
    }
    // Which group of each the calls pair, by the map's broadcast
    vial_kctx_t k = {.in = in, .op = b, .type = &ir->values[op->result].type, .alloc = alloc};
    uint32_t* ix[2];
    for (int w = 0; w < 2; ++w) ix[w] = md_alloc(alloc, (count ? count : 1) * sizeof(uint32_t));
    for (size_t i = 0; i < count; ++i) {
        size_t rem = i;
        uint32_t idx[VIAL_MAX_DIMS];
        for (int d = lnd - 1; d >= 0; --d) { idx[d] = (uint32_t)(rem % ldims[d]); rem /= ldims[d]; }
        for (int w = 0; w < 2; ++w) {
            const vial_val_t* g = av[par[w]];
            const uint32_t m = mapping[par[w]];
            size_t lin = 0;
            for (uint32_t d = 0; d < m; ++d) {
                const uint32_t ext = g->dims[d];
                lin = lin * ext + (ext == 1 ? 0 : idx[lnd - m + d]);
            }
            ix[w][i] = (uint32_t)lin;
        }
    }
    vial_val_t o;
    if (!vial_k_make(&k, VIAL_DT_FLOAT, lnd, ldims, &o)) return 0;
    if (!vial_k_pair_batch(&k, closest, av[par[0]]->sets, av[par[1]]->sets, ix[0], ix[1], count, &o)) return 0;
    *out = o;
    return 1;
}

static bool eval_map(vial_interp_t* in, const vial_ir_op_t* op, const vial_val_t* const* av, vial_val_t* out, md_allocator_i* alloc) {
    vial_kctx_t k = {.in = in, .op = op, .type = &in->ir->values[op->result].type, .alloc = alloc};
    const uint32_t* mapping = op->num_args ? in->ir->ints + op->attr : NULL;
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
    if (!vial_k_room(&k, count)) return false;
    const int batched = eval_map_pairs(in, op, av, mapping, ldims, lnd, count, out, alloc);
    if (batched >= 0) return batched == 1;
    const vial_ir_region_t* r = &in->ir->regions[op->region[0]];
    const bool filtered = op->sub == VIAL_MAP_FILTERED;
    const vial_type_t* rt = k.type;
    md_allocator_i* scratch = map_scratch(in);
    md_allocator_i* ea = scratch ? scratch : alloc;
    const size_t entry = scratch ? md_vm_arena_get_pos(scratch) : 0;
    in->map_depth += 1;

    // The result: the loop's dimensions (one, when filtered), then one call's, made when the first call gives its shape
    vial_val_t o = {.dtype = rt->dtype};
    // Filtered: room for every call while they run, of which the kept are the first; cut to them at the end
    if (filtered) { o.ndim = 1; o.dims[0] = (uint32_t)count; }
    else { o.ndim = lnd; memcpy(o.dims, ldims, lnd * sizeof(uint32_t)); }
    const uint8_t outer_nd = o.ndim;
    const bool has_data = o.dtype != VIAL_DT_SET && o.dtype != VIAL_DT_TUPLE && o.dtype != VIAL_DT_NONE && o.dtype != VIAL_DT_DRAWING;
    bool shaped = false;
    size_t kept = 0;
    uint32_t* set_off = NULL;   // md_array in alloc: an array of sets, gathered call by call
    uint32_t* set_idx = NULL;
    if (o.dtype == VIAL_DT_SET) md_array_push(set_off, 0, alloc);
    vial_val_t* params = md_alloc(alloc, (r->num_params ? r->num_params : 1) * sizeof(vial_val_t));
    bool ok = true;
    for (size_t i = 0; i < count && ok; ++i) {
        // Interruptible within a frame: a map can be most of a frame's work
        if ((i & 255) == 0 && in->stop && md_atomic_load64(in->stop)) { in->interrupted = true; ok = false; break; }
        const size_t mark = scratch ? md_vm_arena_get_pos(scratch) : 0;
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
            params[p] = element(av[p], m, lin, ea);
        }
        vial_val_t ys[2];
        uint32_t ny = 0;
        if (!eval_region(in, op->region[0], params, ys, 2, &ny, ea)) { ok = false; break; }
        if (ny < 1) { ok = vial_rt_fail(&k, "R0103", "a map's call gave no value"); break; }
        bool keep = true;
        if (filtered) {
            if (ny < 2 || ys[1].dtype != VIAL_DT_BOOL || !ys[1].data) { ok = vial_rt_fail(&k, "R0103", "a filter gave no bool"); break; }
            keep = *(const uint8_t*)ys[1].data && !vial_val_is_missing(&ys[1], 0);
        }
        if (keep) {
            const vial_val_t* y = &ys[0];
            if (!shaped) {
                shaped = true;
                for (uint8_t q = 0; q < y->ndim && o.ndim < VIAL_MAX_DIMS; ++q) o.dims[o.ndim++] = y->dims[q];
                if (has_data) {
                    vial_val_t made;
                    if (!vial_k_make(&k, o.dtype, o.ndim, o.dims, &made)) { ok = false; break; }
                    o.data = made.data;
                }
            }
            if (o.dtype == VIAL_DT_SET) {
                // Every call's sets fill the shape the first call gave: an array of sets is not ragged
                size_t inner = 1;
                for (uint8_t q = outer_nd; q < o.ndim; ++q) inner *= o.dims[q];
                const vial_sets_t* e = y->sets;
                if (!e || e->count != inner || vial_val_count(y) != inner) {
                    ok = vial_rt_fail(&k, "R0103", "the calls of a map give values of different shapes");
                    break;
                }
                for (uint32_t q = 0; q < e->count; ++q) {
                    const uint32_t sz = vial_sets_size(e, q);
                    md_array_push_array(set_idx, vial_sets_items(e, q), sz, alloc);
                    md_array_push(set_off, (uint32_t)md_array_size(set_idx), alloc);
                }
            } else if (has_data) {
                if (!place(&k, &o, kept, y, ea)) { ok = false; break; }
            }
            kept += 1;
        }
        if (scratch) md_vm_arena_set_pos_back(scratch, mark);
    }
    in->map_depth -= 1;
    if (scratch) md_vm_arena_set_pos_back(scratch, entry);
    if (!ok) return false;
    if (filtered) o.dims[0] = (uint32_t)kept;
    if (!shaped) {
        // No call: the element's shape is the type's, where it is static
        const uint8_t nd = rt->ndim > outer_nd ? rt->ndim - outer_nd : 0;
        for (uint8_t q = 0; q < nd; ++q) o.dims[o.ndim++] = rt->dims[outer_nd + q].kind == VIAL_EXT_STATIC ? rt->dims[outer_nd + q].n : 0;
        if (has_data) {
            vial_val_t made;
            if (!vial_k_make(&k, o.dtype, o.ndim, o.dims, &made)) return false;
            o.data = made.data;
        }
    }
    if (o.dtype == VIAL_DT_SET) {
        // An array of sets: each call's sets, in order
        vial_sets_t* s = md_alloc(alloc, sizeof(vial_sets_t));
        memset(s, 0, sizeof(*s));
        s->count = (uint32_t)md_array_size(set_off) - 1;
        s->offset = set_off;
        s->idx = set_idx ? set_idx : md_alloc(alloc, sizeof(uint32_t));
        if (o.ndim == 0) s->single = true;
        o.sets = s;
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
        return vial_k_builtin(&k, b->name, av, op->num_args ? ir->names + op->arg_beg : NULL, op->num_args, out);
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
    case VIAL_IR_PAIRS:
        return vial_k_pairs(&k, op->sub, ir->ints[op->attr], ir->ints[op->attr + 1], av[0], av[1], av[2], out);
    case VIAL_IR_IF: {
        if (av[0]->dtype != VIAL_DT_BOOL || !av[0]->data) return vial_rt_fail(&k, "R0103", "an if without a bool");
        const bool c = *(const uint8_t*)av[0]->data;
        uint32_t ny = 0;
        return eval_region(in, op->region[c ? 0 : 1], NULL, out, op->num_results, &ny, alloc);
    }
    case VIAL_IR_FOLD:
        return vial_rt_fail(&k, "R0103", "a fold is evaluated over the trajectory (vial_interp_run), not one frame alone");
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
        if (!in->live[opi]) continue;
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
    md_vm_arena_reset(in->temp);
    const bool st = eval_top(in, true);
    in->frame = frame;
    return eval_top(in, false) && st;
}

// ### THE TRAJECTORY ###
// A static op is computed once. A per-frame op runs for each frame, in the sweep its inputs allow: a value that reads
// a whole series, such as ends - mean(ends, axis='frame'), runs in the sweep after the series ends. A once op runs
// after a sweep: one that reads per-frame values as series (mean over frames, smooth, histogram), and one computed from
// whole-trajectory values alone (eigh(C)). A fold is a per-frame op whose state carries from frame to frame.

static bool reads_frames(const vial_ir_t* ir, const vial_ir_op_t* op) {
    for (uint32_t k = 0; k < op->num_args; ++k) {
        const vial_ssa_t a = ir->args[op->arg_beg + k];
        if (!a) continue;
        const vial_type_t* t = &ir->values[a].type;
        if (t->stage == VIAL_STAGE_FRAME || (t->stage == VIAL_STAGE_TRAJECTORY && (t->flags & VIAL_TYPE_PER_FRAME))) return true;
    }
    return false;
}

static void classify(vial_interp_t* in) {
    const vial_ir_t* ir = in->ir;
    const size_t no = md_array_size(ir->ops);
    md_array_resize(in->op_class, no, in->alloc);
    md_array_resize(in->op_sweep, no, in->alloc);
    memset(in->op_class, 0, no);
    memset(in->op_sweep, 0, no);
    const vial_ir_region_t* r = &ir->regions[ir->funcs[0].body];
    for (size_t i = 0; i < md_array_size(r->ops); ++i) {
        const uint32_t opi = r->ops[i];
        const vial_ir_op_t* op = &ir->ops[opi];
        if (!op->num_results || op->kind == VIAL_IR_YIELD) continue;
        const vial_type_t* t = &ir->values[op->result].type;
        uint8_t cls = VIAL_OP_STATIC;
        if (t->stage == VIAL_STAGE_FRAME) cls = VIAL_OP_PER_FRAME;
        else if (t->stage == VIAL_STAGE_TRAJECTORY) {
            // A window over frames (smooth, derivative) keeps the frame but needs the whole series; every other value
            // that keeps the frame is computed frame by frame, a fold in order
            const char* bname = op->kind == VIAL_IR_BUILTIN ? vial_builtin_at(op->aux)->name : "";
            const bool window = (!strcmp(bname, "smooth") || !strcmp(bname, "derivative")) && reads_frames(ir, op);
            cls = op->kind == VIAL_IR_FOLD ? VIAL_OP_PER_FRAME :
                  window ? VIAL_OP_ONCE :
                  (t->flags & VIAL_TYPE_PER_FRAME) ? VIAL_OP_PER_FRAME : VIAL_OP_ONCE;
        }
        in->op_class[opi] = cls;
        // The sweep: a per-frame op after its per-frame inputs' sweeps and after the sweeps once inputs follow
        uint8_t sweep = cls == VIAL_OP_STATIC ? 0 : 1;
        for (uint32_t k = 0; k < op->num_args; ++k) {
            const vial_ssa_t a = ir->args[op->arg_beg + k];
            if (!a || ir->values[a].op == VIAL_IR_NO_OP) continue;
            const uint32_t src = ir->values[a].op;
            const uint8_t sc = in->op_class[src], ss = in->op_sweep[src];
            uint8_t need = 0;
            if (sc == VIAL_OP_PER_FRAME) need = ss;
            else if (sc == VIAL_OP_ONCE) need = cls == VIAL_OP_PER_FRAME ? ss + 1 : ss;
            if (need > sweep) sweep = need;
        }
        in->op_sweep[opi] = sweep;
    }
}

const vial_val_t* vial_interp_at(const vial_interp_t* in, vial_ssa_t v, int64_t f) {
    if (!v || v >= md_array_size(in->vals)) return NULL;
    const uint32_t op = in->ir->values[v].op;
    if (op != VIAL_IR_NO_OP && in->op_class && in->op_class[op] == VIAL_OP_PER_FRAME) {
        if (!in->series || !in->series[v] || f < 0 || f >= in->num_frames) return NULL;
        return in->series[v][f].dtype || in->series[v][f].ndim ? &in->series[v][f] : NULL;
    }
    return vial_interp_value(in, v);
}

// One frame of a series, made by a once op such as smooth: the value at frame f
static vial_val_t frame_of(const vial_val_t* s, int64_t f, md_allocator_i* alloc) {
    vial_val_t e = *s;
    e.flags &= ~VIAL_VAL_FRAMES;
    e.ndim = s->ndim - 1;
    memmove(e.dims, s->dims + 1, e.ndim * sizeof(uint32_t));
    size_t inner = 1;
    for (uint8_t k = 0; k < e.ndim; ++k) inner *= e.dims[k];
    if (s->data) e.data = (uint8_t*)s->data + (size_t)f * inner * vial_dtype_size(s->dtype);
    e.missing = NULL;
    for (size_t k = 0; k < inner; ++k) if (vial_val_is_missing(s, (size_t)f * inner + k)) vial_val_set_missing(&e, k, alloc);
    return e;
}

// A per-frame value's series, stacked: [frames] and its shape, flagged as frames
static bool stacked(vial_interp_t* in, vial_ssa_t a, vial_val_t* out) {
    const vial_val_t* s = in->series ? in->series[a] : NULL;
    if (!s) return false;
    const int64_t nf = in->num_frames;
    const vial_val_t* first = &s[0];
    vial_val_t o = {.dtype = first->dtype, .ndim = first->ndim + 1, .flags = VIAL_VAL_FRAMES};
    o.dims[0] = (uint32_t)nf;
    memcpy(o.dims + 1, first->dims, first->ndim * sizeof(uint32_t));
    const size_t inner = vial_val_count(first), es = vial_dtype_size(first->dtype);
    const vial_val_t made = vial_val_make(o.dtype, o.ndim, o.dims, in->alloc);
    o.data = made.data;
    for (int64_t f = 0; f < nf; ++f) {
        const vial_val_t* v = &s[f];
        if (vial_val_count(v) != inner || v->dtype != first->dtype) {
            // A frame that was not computed: missing throughout
            for (size_t k = 0; k < inner; ++k) vial_val_set_missing(&o, (size_t)f * inner + k, in->alloc);
            continue;
        }
        if (es && v->data) memcpy((uint8_t*)o.data + (size_t)f * inner * es, v->data, inner * es);
        for (size_t k = 0; k < inner; ++k) if (vial_val_is_missing(v, k)) vial_val_set_missing(&o, (size_t)f * inner + k, in->alloc);
    }
    *out = o;
    return true;
}

// A fold's step for this frame: the state after it, or the state carried over when x is missing and skipped
static bool eval_fold(vial_interp_t* in, uint32_t opi) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_op_t* op = &ir->ops[opi];
    vial_kctx_t k = {.in = in, .op = op, .map_index = -1};
    const vial_ssa_t r = op->result;
    const vial_val_t* x = &in->vals[ir->args[op->arg_beg]];
    const vial_val_t* init = &in->vals[ir->args[op->arg_beg + 1]];
    if (md_array_size(in->fold_state) <= r) {
        const size_t old = md_array_size(in->fold_state);
        md_array_resize(in->fold_state, r + 1, in->alloc);
        md_array_resize(in->fold_has, r + 1, in->alloc);
        memset(in->fold_state + old, 0, (r + 1 - old) * sizeof(vial_val_t));
        memset(in->fold_has + old, 0, r + 1 - old);
    }
    bool x_missing = false;
    for (size_t i = 0; i < vial_val_count(x); ++i) x_missing = x_missing || vial_val_is_missing(x, i);
    const bool skip = op->sub != 0;
    if (!in->fold_has[r]) {
        // The initial state: init, taken from the first frame where it is present
        bool init_missing = false;
        for (size_t i = 0; i < vial_val_count(init); ++i) init_missing = init_missing || vial_val_is_missing(init, i);
        if (init_missing) { in->vals[r] = *init; return true; }
        in->fold_state[r] = vial_val_copy(init, in->alloc);
        in->fold_has[r] = 1;
    }
    if (!(skip && x_missing)) {
        vial_val_t params[2] = {in->fold_state[r], *x};
        vial_val_t next;
        uint32_t ny = 0;
        if (!eval_region_pub(in, op->region[0], params, &next, 1, &ny, in->temp)) return false;
        if (ny != 1) return vial_rt_fail(&k, "R0103", "a fold's step gave no state");
        in->fold_state[r] = vial_val_copy(&next, in->alloc);
    }
    in->vals[r] = in->fold_state[r];
    return true;
}

static bool needs_series(const vial_interp_t* in, vial_ssa_t v) {
    const vial_ir_t* ir = in->ir;
    // Outputs, and what a once op or a later sweep reads
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (ir->outputs[i].value == v) return true;
    return false;
}

enum { FRAME_OK = 0, FRAME_FAILED = 1, FRAME_STOP = 2 };

// One frame of a sweep, in `in`, which is main or one of its workers: the per-frame ops of this sweep, and what earlier
// sweeps computed read at this frame. FRAME_FAILED when an op failed, which leaves the rest of the frame computed;
// FRAME_STOP when the frame does not load or the evaluation is interrupted.
static int run_frame(vial_interp_t* in, vial_interp_t* main, uint8_t sweep, int64_t f) {
    const vial_ir_t* ir = in->ir;
    const vial_ir_region_t* top = &ir->regions[ir->funcs[0].body];
    const size_t nv = md_array_size(ir->values);
    if (main->stop && md_atomic_load64(main->stop)) return FRAME_STOP;
    md_atomic_fetch_add64(&main->frames_done, 1);
    in->sweep = sweep;
    in->frame = NULL;
    in->err_seq = 0;
    vial_frame_t fr;
    if (!in->load_frame || !in->load_frame(in->load_user, f, &fr)) {
        vial_frame_t at = {.index = f};
        in->frame = &at;
        vial_kctx_t k = {.in = in, .map_index = -1};
        vial_rt_fail(&k, "R0104", "frame %lld does not load", (long long)f);
        in->frame = NULL;
        return FRAME_STOP;
    }
    main->times[f] = fr.time;
    in->frame = &fr;
    for (size_t v = 0; v < nv; ++v) if (in->state[v] == 1) in->state[v] = 0;
    md_vm_arena_reset(in->temp);
    int result = FRAME_OK;
    for (size_t i = 0; i < md_array_size(top->ops); ++i) {
        const uint32_t opi = top->ops[i];
        const vial_ir_op_t* op = &ir->ops[opi];
        if (!op->num_results || !in->live[opi]) continue;
        // Static ops were computed before the sweeps; once ops after an earlier sweep, their series read below
        if (in->op_class[opi] != VIAL_OP_PER_FRAME) continue;
        const uint8_t sw = in->op_sweep[opi];
        if (sw > sweep) continue;
        if (sw < sweep) {
            // Computed in an earlier sweep: its value at this frame, when kept
            if (in->series[op->result] && (in->series[op->result][f].dtype || in->series[op->result][f].ndim)) {
                in->vals[op->result] = in->series[op->result][f];
                in->state[op->result] = 1;
            }
            continue;
        }
        // Inputs that are series made after an earlier sweep are read at this frame
        bool ready = true;
        const vial_val_t* saved_vals[64];
        vial_ssa_t saved_ssa[64];
        uint32_t ns = 0;
        for (uint32_t k = 0; k < op->num_args; ++k) {
            const vial_ssa_t a = ir->args[op->arg_beg + k];
            if (!a) continue;
            if (!in->state[a]) { ready = false; break; }
            if (in->state[a] == 3 && (in->vals[a].flags & VIAL_VAL_FRAMES) && ns < 64) {
                saved_vals[ns] = md_alloc(in->temp, sizeof(vial_val_t));
                *(vial_val_t*)saved_vals[ns] = in->vals[a];
                saved_ssa[ns++] = a;
                in->vals[a] = frame_of(&in->vals[a], f, in->temp);
            }
        }
        bool done = false;
        if (ready) done = op->kind == VIAL_IR_FOLD ? eval_fold(in, opi) : eval_op_pub(in, opi, in->temp);
        for (uint32_t q = 0; q < ns; ++q) in->vals[saved_ssa[q]] = *saved_vals[q];
        if (!done) { result = FRAME_FAILED; continue; }
        in->state[op->result] = 1;
        if (in->series[op->result]) in->series[op->result][f] = vial_val_copy(&in->vals[op->result], in->alloc);
    }
    in->frame = NULL;
    return result;
}

// A sweep's frames in chunks, each taken by a free worker
typedef struct sweep_job_t {
    vial_interp_t*      main;
    uint8_t             sweep;
    int64_t             num_frames;
    int64_t             chunk;
    volatile uint64_t   failed;
    volatile uint64_t   stop;
} sweep_job_t;

static void sweep_chunk(void* arg, size_t c) {
    sweep_job_t* job = arg;
    vial_interp_t* main = job->main;
    // A free worker: there are as many as calls run at once, so one is always free soon
    vial_interp_t* w = NULL;
    while (!w) {
        for (size_t q = 0; q < main->num_workers && !w; ++q) {
            vial_interp_t* cand = &main->workers[q];
            if (md_atomic_fetch_add64(&cand->busy, 1) == 0) w = cand;
            else md_atomic_fetch_add64(&cand->busy, UINT64_MAX);
        }
    }
    const int64_t beg = (int64_t)c * job->chunk, end = MIN(beg + job->chunk, job->num_frames);
    for (int64_t f = beg; f < end && !md_atomic_load64(&job->stop); ++f) {
        const int r = run_frame(w, main, job->sweep, f);
        if (r == FRAME_FAILED) md_atomic_fetch_add64(&job->failed, 1);
        if (r == FRAME_STOP) { md_atomic_fetch_add64(&job->failed, 1); md_atomic_fetch_add64(&job->stop, 1); }
    }
    if (w->load_release) w->load_release(w->load_user);
    md_atomic_fetch_add64(&w->busy, UINT64_MAX);
}

// A worker ready for a sweep: what main holds that is not per frame, static values and what earlier sweeps finished
static void worker_sync(vial_interp_t* w, const vial_interp_t* main) {
    const size_t nv = md_array_size(main->state);
    for (size_t v = 0; v < nv; ++v) {
        if (main->state[v] >= 2) { w->vals[v] = main->vals[v]; w->state[v] = main->state[v]; }
        else w->state[v] = 0;
    }
    w->op_class = main->op_class;
    w->op_sweep = main->op_sweep;
    w->series = main->series;
    w->num_frames = main->num_frames;
    w->times = main->times;
    w->stop = main->stop;
    w->max_elements = main->max_elements;
    w->plain = main->plain;
    memcpy(w->live, main->live, md_array_size(main->live));
}

static int cmp_error(const void* a, const void* b) {
    const vial_rt_error_t* x = a;
    const vial_rt_error_t* y = b;
    if (x->sweep != y->sweep) return x->sweep < y->sweep ? -1 : 1;
    if (x->frame != y->frame) return x->frame < y->frame ? -1 : 1;
    return x->seq < y->seq ? -1 : x->seq > y->seq;
}

// Every interpreter's errors as one: in the order one thread meets them, each failure once
static void merge_errors(vial_interp_t* in) {
    vial_rt_error_t* all = NULL;
    for (size_t i = 0; i < md_array_size(in->errors); ++i) md_array_push(all, in->errors[i], in->alloc);
    for (size_t w = 0; w < in->num_workers; ++w) {
        const vial_interp_t* wk = &in->workers[w];
        for (size_t i = 0; i < md_array_size(wk->errors); ++i) md_array_push(all, wk->errors[i], in->alloc);
    }
    const size_t n = md_array_size(all);
    if (n > 1) qsort(all, n, sizeof(vial_rt_error_t), cmp_error);
    md_array_shrink(in->errors, 0);
    for (size_t i = 0; i < n; ++i) {
        bool seen = false;
        for (size_t j = 0; j < md_array_size(in->errors) && !seen; ++j) {
            const vial_rt_error_t* g = &in->errors[j];
            seen = g->span.beg == all[i].span.beg && g->span.end == all[i].span.end && !strcmp(g->message, all[i].message);
        }
        if (!seen) md_array_push(in->errors, all[i], in->alloc);
    }
    if (md_array_size(in->errors)) { in->failed = true; in->error = in->errors[0]; }
}

bool vial_interp_run(vial_interp_t* in, int64_t num_frames) {
    const vial_ir_t* ir = in->ir;
    classify(in);
    in->num_frames = num_frames;
    in->times = md_alloc(in->alloc, (num_frames ? num_frames : 1) * sizeof(double));
    in->sweep = 0;
    bool ok = vial_interp_static(in);
    // Which per-frame values are kept: outputs, and those read across sweeps or as series
    const size_t nv = md_array_size(ir->values);
    uint8_t* keep = md_alloc(in->alloc, nv ? nv : 1);
    memset(keep, 0, nv);
    uint8_t max_sweep = 0;
    const vial_ir_region_t* top = &ir->regions[ir->funcs[0].body];
    for (size_t i = 0; i < md_array_size(top->ops); ++i) {
        const uint32_t opi = top->ops[i];
        const vial_ir_op_t* op = &ir->ops[opi];
        if (!op->num_results || !in->live[opi]) continue;
        if (in->op_sweep[opi] > max_sweep) max_sweep = in->op_sweep[opi];
        if (in->op_class[opi] == VIAL_OP_PER_FRAME && needs_series(in, op->result)) keep[op->result] = 1;
        for (uint32_t k = 0; k < op->num_args; ++k) {
            const vial_ssa_t a = ir->args[op->arg_beg + k];
            if (!a || ir->values[a].op == VIAL_IR_NO_OP) continue;
            const uint32_t src = ir->values[a].op;
            if (in->op_class[src] != VIAL_OP_PER_FRAME) continue;
            if (in->op_class[opi] == VIAL_OP_ONCE || in->op_sweep[opi] > in->op_sweep[src]) keep[a] = 1;
        }
    }
    md_array_resize(in->series, nv, in->alloc);
    memset(in->series, 0, nv * sizeof(*in->series));
    for (size_t v = 0; v < nv; ++v) {
        if (!keep[v]) continue;
        in->series[v] = md_alloc(in->alloc, (num_frames ? num_frames : 1) * sizeof(vial_val_t));
        memset(in->series[v], 0, (num_frames ? num_frames : 1) * sizeof(vial_val_t));
    }
    // The sweeps with per-frame ops, and those of them with a fold, whose frames go in order
    uint8_t has_frames[256] = {0}, has_fold[256] = {0};
    for (size_t i = 0; i < md_array_size(top->ops); ++i) {
        const uint32_t opi = top->ops[i];
        if (!in->live[opi] || in->op_class[opi] != VIAL_OP_PER_FRAME) continue;
        has_frames[in->op_sweep[opi]] = 1;
        if (ir->ops[opi].kind == VIAL_IR_FOLD) has_fold[in->op_sweep[opi]] = 1;
    }
    int64_t sweeps_with_frames = 0;
    for (uint32_t sweep = 1; sweep <= max_sweep; ++sweep) sweeps_with_frames += has_frames[sweep];
    in->frames_done = 0;
    in->frames_total = (uint64_t)(sweeps_with_frames * num_frames);
    in->interrupted = false;
    bool stopped = false;
    for (uint32_t sweep = 1; sweep <= max_sweep && !stopped; ++sweep) {
        in->sweep = (uint8_t)sweep;
        if (has_frames[sweep] && num_frames > 0) {
            if (in->num_workers && in->parallel_for && !has_fold[sweep]) {
                for (size_t w = 0; w < in->num_workers; ++w) worker_sync(&in->workers[w], in);
                // Chunks of a few frames: enough of them that the workers finish together
                const int64_t chunk = MAX((int64_t)1, num_frames / (int64_t)(in->num_workers * 8));
                sweep_job_t job = {.main = in, .sweep = (uint8_t)sweep, .num_frames = num_frames, .chunk = chunk};
                in->parallel_for(in->parallel_user, (size_t)((num_frames + chunk - 1) / chunk), sweep_chunk, &job);
                if (md_atomic_load64(&job.failed)) ok = false;
                if (md_atomic_load64(&job.stop)) stopped = true;
            } else {
                for (int64_t f = 0; f < num_frames; ++f) {
                    const int r = run_frame(in, in, (uint8_t)sweep, f);
                    if (r != FRAME_OK) ok = false;
                    if (r == FRAME_STOP) { stopped = true; break; }
                }
            }
        }
        if (stopped) break;
        // After the sweep: the once ops it completes, in program order
        in->frame = NULL;
        in->err_seq = 0;
        for (size_t i = 0; i < md_array_size(top->ops); ++i) {
            const uint32_t opi = top->ops[i];
            const vial_ir_op_t* op = &ir->ops[opi];
            if (!op->num_results || !in->live[opi] || in->op_class[opi] != VIAL_OP_ONCE || in->op_sweep[opi] != sweep) continue;
            vial_ssa_t replaced[64];
            vial_val_t old[64];
            uint32_t nr = 0;
            bool ready = true;
            for (uint32_t k = 0; k < op->num_args && ready; ++k) {
                const vial_ssa_t a = ir->args[op->arg_beg + k];
                if (!a || ir->values[a].op == VIAL_IR_NO_OP) continue;
                const uint32_t src = ir->values[a].op;
                if (in->op_class[src] == VIAL_OP_PER_FRAME) {
                    vial_val_t st;
                    if (!stacked(in, a, &st)) { ready = false; break; }
                    if (nr < 64) { replaced[nr] = a; old[nr++] = in->vals[a]; }
                    in->vals[a] = st;
                    in->state[a] = 1;
                } else if (!in->state[a]) {
                    ready = false;
                }
            }
            bool done = false;
            if (ready) done = eval_op_pub(in, opi, in->alloc);
            for (uint32_t q = 0; q < nr; ++q) in->vals[replaced[q]] = old[q];
            if (!done) { ok = false; continue; }
            in->state[op->result] = 3;
        }
    }
    in->frame = NULL;
    if (in->load_release) in->load_release(in->load_user);
    merge_errors(in);
    if (stopped && in->stop && md_atomic_load64(in->stop)) { in->interrupted = true; return false; }
    return ok && !stopped;
}

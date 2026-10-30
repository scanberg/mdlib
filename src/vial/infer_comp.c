// The checker, step 7: comprehensions and loops.
//
// A comprehension iterates over an array (its first dimension), a set (its particles, one at a time), range() or
// zip(). Its result has the iterable's extent, with its description, so `[mean(rmsf[r]) for r in residues]` is labelled
// by residue; an `if` gives a new extent. When the iterable's value is known and the element is a set, the
// comprehension is evaluated when compiling, element by element: `[particles(r, name='CA') for r in res_a]` is known to
// hold one particle per set, so it is points.

#include "infer_internal.h"

#include <core/md_common.h>
#include <core/md_array.h>

#include <stdio.h>
#include <string.h>

typedef vial_type_t T;

static T err(void) { return vial_type_error(); }

#define MAX_EVALUATED 4096      // elements a comprehension is evaluated for when compiling

// What iterating gives: the element's type, the extent, and per element its value when known
typedef struct iter_t {
    T               elem;
    vial_extent_t   extent;
    uint8_t         stage;
    // Known values: sets (an array's sets, or a set's particles), or ints from range()
    const vial_sets_t* sets;
    bool            particles;  // sets is one set, iterated particle by particle
    bool            ints;
    int64_t         lo, step;
    uint32_t        count;      // when known
    bool            known;
} iter_t;

static bool callee_is(infer_t* x, vial_node_id node, const char* name) {
    if (!node || vx_kind(x, node) != VIAL_N_CALL) return false;
    const vial_node_id c = vx_child(x, node, 0);
    if (!c || vx_kind(x, c) != VIAL_N_NAME) return false;
    const vial_symbol_t* s = vial_names_symbol(x->n, c);
    return s && s->kind == VIAL_SYM_BUILTIN && s->builtin && strcmp(s->builtin->name, name) == 0;
}

static bool int_value(infer_t* x, vial_node_id n, int64_t* out) {
    const vial_value_t v = vx_value_of(x, n);
    if (v.kind != VIAL_VAL_INT) return false;
    *out = v.i;
    return true;
}

static bool iterate(infer_t* x, vial_node_id it, iter_t* out);

// range(n), range(a, b), range(a, b, step): ints, their count known when the bounds are
static bool iterate_range(infer_t* x, vial_node_id call, iter_t* out) {
    vial_node_id a[3] = {0};
    uint32_t n = 0;
    bool ok = true;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id c = vx_child(x, call, i);
        if (!c) continue;
        if (vx_kind(x, c) == VIAL_N_KWARG || n == 3) {
            vx_diag_fix(x, "E0204", c, NULL, "`range` takes 1 to 3 ints");
            ok = false;
            continue;
        }
        a[n++] = c;
        const T t = vx_infer(x, c);
        if (vial_type_is_error(&t)) { ok = false; continue; }
        if (t.dtype != VIAL_DT_INT || t.ndim) { vx_diag_fix(x, "E0313", c, NULL, "`range` takes ints; this is %s", vx_tstr(&t)); ok = false; }
        if (t.stage > out->stage) out->stage = t.stage;
    }
    if (!ok) return false;
    if (n == 0) { vx_diag_fix(x, "E0204", call, NULL, "`range` takes 1 to 3 ints"); return false; }
    out->elem = vial_type_scalar(VIAL_DT_INT);
    out->elem.stage = out->stage;
    int64_t lo = 0, hi = 0, step = 1;
    const bool known = (n == 1 ? int_value(x, a[0], &hi) : int_value(x, a[0], &lo) && int_value(x, a[1], &hi)) &&
                       (n < 3 || int_value(x, a[2], &step)) && step != 0;
    if (known && out->stage == VIAL_STAGE_STATIC) {
        int64_t c = step > 0 ? (hi - lo + step - 1) / step : (lo - hi - step - 1) / (-step);
        if (c < 0) c = 0;
        out->extent = vial_extent_static((uint32_t)c);
        out->ints = true; out->lo = lo; out->step = step; out->count = (uint32_t)c; out->known = true;
    } else {
        out->extent = (vial_extent_t){.kind = out->stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call};
    }
    return true;
}

// zip(a, b, ...): a tuple of their elements, which pair in order, so their extents must match
static bool iterate_zip(infer_t* x, vial_node_id call, iter_t* out) {
    uint32_t n = 0;
    T* elems = NULL;
    bool ok = true;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id c = vx_child(x, call, i);
        if (!c) continue;
        if (vx_kind(x, c) == VIAL_N_KWARG) { vx_diag_fix(x, "E0202", c, NULL, "`zip` takes no keywords"); ok = false; continue; }
        iter_t one = {0};
        if (!iterate(x, c, &one)) { ok = false; continue; }
        if (n == 0) {
            out->extent = one.extent;
        } else {
            vial_extent_t e;
            if (!vx_extent_join(out->extent, one.extent, &e) || (out->extent.kind == VIAL_EXT_STATIC && one.extent.kind == VIAL_EXT_STATIC &&
                out->extent.n != one.extent.n)) {
                vx_diag_fix(x, "E0301", c, NULL, "`zip` pairs elements in order, and these have %u and %u", out->extent.n, one.extent.n);
                ok = false;
                continue;
            }
            out->extent = e;
        }
        if (one.stage > out->stage) out->stage = one.stage;
        md_array_push(elems, one.elem, x->alloc);
        n += 1;
    }
    if (!ok) return false;
    if (n < 2) { vx_diag_fix(x, "E0204", call, NULL, "`zip` pairs two or more arrays"); return false; }
    out->elem = vial_type_scalar(VIAL_DT_TUPLE);
    out->elem.elems = elems;
    out->elem.num_elems = n;
    return true;
}

static bool iterate(infer_t* x, vial_node_id it, iter_t* out) {
    memset(out, 0, sizeof(*out));
    if (callee_is(x, it, "range")) return iterate_range(x, it, out);
    if (callee_is(x, it, "zip")) return iterate_zip(x, it, out);
    const T a = vx_infer(x, it);
    if (vial_type_is_error(&a)) return false;
    out->stage = a.stage;
    const vial_value_t v = vx_value_of(x, it);
    if (a.dtype == VIAL_DT_SET && a.ndim == 0) {
        // A set's particles, one at a time
        out->elem = vial_type_scalar(VIAL_DT_SET);
        out->elem.flags |= VIAL_TYPE_ONE_PARTICLE;
        out->elem.stage = a.stage;
        if (v.kind == VIAL_VAL_SETS && a.stage == VIAL_STAGE_STATIC) {
            out->count = vial_sets_size(v.sets, 0);
            out->extent = vial_extent_static(out->count);
            out->extent.desc = vx_desc(x, VIAL_LEVEL_PARTICLE, v.sets);
            out->sets = v.sets; out->particles = true; out->known = true;
        } else {
            out->extent = (vial_extent_t){.kind = a.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = a.ident ? a.ident : (uint32_t)it};
        }
        return true;
    }
    if (a.ndim >= 1 && a.dtype != VIAL_DT_TUPLE && a.dtype != VIAL_DT_DEF) {
        out->elem = a;
        out->elem.ndim -= 1;
        memmove(out->elem.dims, a.dims + 1, out->elem.ndim * sizeof(a.dims[0]));
        if (a.dtype == VIAL_DT_SET) {
            out->elem.flags &= ~VIAL_TYPE_EMPTY;
        }
        out->extent = a.dims[0];
        if (a.dtype == VIAL_DT_SET && a.ndim == 1 && v.kind == VIAL_VAL_SETS && a.stage == VIAL_STAGE_STATIC && !v.sets->single) {
            out->sets = v.sets; out->count = v.sets->count; out->known = true;
        }
        return true;
    }
    vx_diag_fix(x, "E0313", it, NULL, "a comprehension iterates over an array, a set, range() or zip(); this is %s", vx_tstr(&a));
    return false;
}

// Binds a comprehension's or loop's target in the frame
static bool bind(infer_t* x, vial_node_id target, const T* t, const vial_value_t* v) {
    if (!target) return false;
    const vial_node_kind_t k = vx_kind(x, target);
    if (k == VIAL_N_NAME) {
        x->out->types[(target) - x->out->base] = *t;
        x->out->typed[(target) - x->out->base] = 1;
        x->out->values[(target) - x->out->base] = v ? *v : (vial_value_t){0};
        const uint32_t sym = x->n->ref[target];
        const vial_value_t none = {0};
        if (sym) vx_frame_set(x, sym, t, v ? v : &none, 0);
        return true;
    }
    if (k == VIAL_N_TUPLE || k == VIAL_N_LIST) {
        const uint32_t n = vx_num_children(x, target);
        if (t->dtype != VIAL_DT_TUPLE || t->num_elems != n) {
            vx_diag_fix(x, "E0319", target, NULL, "%u names, and %s does not unpack into %u", n, vx_tstr(t), n);
            return false;
        }
        bool ok = true;
        for (uint32_t i = 0; i < n; ++i) ok = bind(x, vx_child(x, target, i), &t->elems[i], NULL) && ok;
        return ok;
    }
    return false;
}

// Every node of a subtree, to check it again for another element
static void subtree(infer_t* x, vial_node_id root, vial_node_id** out) {
    vial_node_id* stack = NULL;
    md_array_push(stack, root, x->alloc);
    while (md_array_size(stack)) {
        const vial_node_id n = *md_array_last(stack);
        md_array_pop(stack);
        md_array_push(*out, n, x->alloc);
        for (uint32_t i = 0; i < vx_num_children(x, n); ++i) {
            const vial_node_id c = vx_child(x, n, i);
            if (c) md_array_push(stack, c, x->alloc);
        }
    }
}

static vx_frame_t child_frame(const infer_t* x, vial_node_id at) {
    vx_frame_t f = {.parent = x->frame, .t = x->t, .def = at, .depth = x->frame ? x->frame->depth : 0};
    f.prelude = x->frame && x->frame->prelude;
    f.site = x->frame ? x->frame->site : 0;
    f.syms.allocator = x->alloc;
    return f;
}

T vx_listcomp(infer_t* x, vial_node_id id) {
    const uint32_t nc = vx_num_children(x, id);
    const vial_node_id elem = nc ? vx_child(x, id, 0) : 0;
    if (!elem) return err();
    vx_frame_t f = child_frame(x, id);
    infer_t y = *x;
    y.frame = &f;

    T dims = vial_type_scalar(VIAL_DT_INT);
    uint8_t stage = VIAL_STAGE_STATIC;
    iter_t first = {0};
    uint32_t num_for = 0;
    bool filtered = false, ok = true;
    for (uint32_t i = 1; i < nc; ++i) {
        const vial_node_id cf = vx_child(x, id, i);
        if (!cf || vx_kind(x, cf) != VIAL_N_COMP_FOR) continue;
        const vial_node_id target = vx_child(x, cf, 0), it = vx_child(x, cf, 1);
        iter_t r;
        if (!it || !iterate(&y, it, &r)) { ok = false; continue; }
        if (num_for == 0) first = r;
        num_for += 1;
        if (r.stage > stage) stage = r.stage;
        if (!bind(&y, target, &r.elem, NULL)) ok = false;
        // A nested for runs inside the one before: the extents multiply
        if (dims.ndim == 0) vial_type_push_dim(&dims, r.extent);
        else if (dims.dims[0].kind == VIAL_EXT_STATIC && r.extent.kind == VIAL_EXT_STATIC) { dims.dims[0].n *= r.extent.n; dims.dims[0].desc = 0; }
        else dims.dims[0] = (vial_extent_t){.kind = stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)id};
        for (uint32_t k = 2; k < vx_num_children(x, cf); ++k) {
            const vial_node_id c = vx_child(x, cf, k);
            if (!c) continue;
            const T ct = vx_infer(&y, c);
            if (vial_type_is_error(&ct)) { ok = false; continue; }
            if (ct.dtype != VIAL_DT_BOOL || ct.ndim) {
                vx_diag_fix(&y, "E0303", c, NULL, "a comprehension's if needs one bool; this is `%s`", vx_tstr(&ct));
                ok = false;
                continue;
            }
            if (ct.stage > stage) stage = ct.stage;
            filtered = true;
        }
    }
    if (!ok || num_for == 0) { vx_infer(&y, elem); return err(); }
    if (filtered) dims.dims[0] = (vial_extent_t){.kind = stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)id};

    const T e = vx_infer(&y, elem);
    if (vial_type_is_error(&e)) return err();
    if (e.dtype == VIAL_DT_TUPLE || e.dtype == VIAL_DT_DEF || e.dtype == VIAL_DT_NONE || e.dtype == VIAL_DT_RANGE) {
        vx_diag_fix(x, "E0319", elem, NULL, "a list holds numbers, strings, sets or drawings; this holds %s", vx_tstr(&e));
        return err();
    }
    T t = e;
    t.ndim = 0;
    vial_type_push_dim(&t, dims.dims[0]);
    for (uint8_t k = 0; k < e.ndim; ++k) {
        if (!vial_type_push_dim(&t, e.dims[k])) {
            vx_diag_fix(x, "E0308", id, NULL, "the result would have 7 dimensions; the limit is 6");
            return err();
        }
    }
    T o = t;
    vx_join_stage(&t, &o, &(T){.stage = stage});
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    if (e.dtype == VIAL_DT_SET && (e.flags & VIAL_TYPE_ONE_PARTICLE)) t.flags |= VIAL_TYPE_ONE_PARTICLE;

    // Evaluated element by element, when everything about it is known: an array of sets
    if (num_for == 1 && !filtered && first.known && first.count <= MAX_EVALUATED && e.dtype == VIAL_DT_SET && e.ndim == 0 &&
        e.stage == VIAL_STAGE_STATIC && t.stage == VIAL_STAGE_STATIC) {
        const vial_node_id cf = vx_child(x, id, 1);
        const vial_node_id target = vx_child(x, cf, 0);
        vial_node_id* nodes = NULL;
        subtree(&y, elem, &nodes);
        uint32_t* offset = md_alloc(x->alloc, (first.count + 1) * sizeof(uint32_t));
        uint32_t* idx = NULL;
        offset[0] = 0;
        bool all_one = first.count > 0, known = true;
        for (uint32_t k = 0; k < first.count && known; ++k) {
            const vial_sets_t* one;
            if (first.particles) one = vial_sets_one(vial_sets_items(first.sets, 0) + k, 1, x->alloc);
            else one = vial_sets_pick(first.sets, &k, 1, x->alloc), ((vial_sets_t*)one)->single = true;
            const vial_value_t v = {.kind = VIAL_VAL_SETS, .sets = one};
            T et = first.elem;
            et.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
            if (vial_sets_size(one, 0) == 1) et.flags |= VIAL_TYPE_ONE_PARTICLE;
            if (vx_kind(&y, target) != VIAL_N_NAME) { known = false; break; }
            bind(&y, target, &et, &v);
            for (size_t q = 0; q < md_array_size(nodes); ++q) y.out->typed[nodes[q] - y.out->base] = 0;
            vx_infer(&y, elem);
            const vial_value_t r = vx_value_of(&y, elem);
            if (r.kind != VIAL_VAL_SETS || !r.sets->single) { known = false; break; }
            const uint32_t sz = vial_sets_size(r.sets, 0);
            all_one = all_one && sz == 1;
            md_array_push_array(idx, vial_sets_items(r.sets, 0), sz, x->alloc);
            offset[k + 1] = (uint32_t)md_array_size(idx);
        }
        if (known) {
            vial_sets_t* s = md_alloc(x->alloc, sizeof(vial_sets_t));
            memset(s, 0, sizeof(*s));
            s->count = first.count;
            s->offset = offset;
            s->idx = idx ? idx : md_alloc(x->alloc, sizeof(uint32_t));
            // The units the iterable's sets are, when it has them: [particles(r, name='CA') for r in res_a] runs over res_a's residues
            if (!first.particles && first.sets->unit) { s->level = first.sets->level; s->unit = first.sets->unit; }
            vx_set_value(x, id, (vial_value_t){.kind = VIAL_VAL_SETS, .sets = s});
            if (all_one) t.flags |= VIAL_TYPE_ONE_PARTICLE;
        }
        // The element's nodes keep the types of one element; the target, the element type
        bind(&y, target, &first.elem, NULL);
    }
    return t;
}

bool vx_check_for(infer_t* x, vial_node_id s) {
    const vial_node_id target = vx_child(x, s, 0), it = vx_child(x, s, 1), body = vx_child(x, s, 2);
    vx_frame_t f = child_frame(x, s);
    infer_t y = *x;
    if (!x->frame) y.frame = &f;    // at top level a loop has its own names; in a def they are the def's
    iter_t r;
    if (!it || !iterate(&y, it, &r)) return false;
    bind(&y, target, &r.elem, NULL);
    vx_check_block(&y, body);
    return false;   // a loop may run no times
}

// ### FRAMES: fold, frame(), time() ###

static T call_fold(infer_t* x, vial_node_id call) {
    static const vx_param_t params[] = {{"step", VX_KIND_WHOLE, true}, {"x", VX_KIND_WHOLE, true}, {"init", VX_KIND_WHOLE, true},
                                         {"skip", VX_KIND_WHOLE, false}};
    vial_node_id a[4];
    if (!vx_bind(x, call, "fold", params, 4, a)) return err();
    // The step is a def, named
    const vial_node_id sn = a[0];
    const vial_symbol_t* sym = vx_kind(x, sn) == VIAL_N_NAME ? vial_names_symbol(x->n, sn) : NULL;
    vial_node_id def = sym && sym->kind == VIAL_SYM_DEF ? sym->node : 0;
    if (def && vx_kind(x, def) == VIAL_N_NAME) def = x->n->parent[def];
    if (!def || vx_kind(x, def) != VIAL_N_DEF) {
        const T st = x->out->types[(sn) - x->out->base];
        if (!vial_type_is_error(&st)) vx_diag_fix(x, "E0313", sn, NULL, "`fold`'s step is a def, step(state, value); this is %s", vx_tstr(&st));
        return err();
    }
    const T v = x->out->types[(a[1]) - x->out->base], init = x->out->types[(a[2]) - x->out->base];
    if (vial_type_is_error(&v) || vial_type_is_error(&init)) return err();
    if (v.stage == VIAL_STAGE_STATIC) {
        vx_diag_fix(x, "E0313", a[1], NULL, "`fold` steps over frames, and %s is static", vx_tstr(&v));
        return err();
    }
    bool skip = true;
    if (a[3]) {
        const vial_value_t sv = vx_value_of(x, a[3]);
        if (sv.kind != VIAL_VAL_BOOL) { vx_diag_fix(x, "E0313", a[3], NULL, "`skip` is True or False, known when compiling"); return err(); }
        skip = sv.b;
    }
    // One state and one frame's value: the state as init has it, any frame
    T s = init;
    s.stage = VIAL_STAGE_STATIC; s.pass = 0; s.flags &= ~VIAL_TYPE_PER_FRAME;
    T one = v;
    one.flags &= ~VIAL_TYPE_PER_FRAME;
    if (!skip) { s.flags |= v.flags & VIAL_TYPE_MISSING; }
    else one.flags &= ~VIAL_TYPE_MISSING;
    const T types[2] = {s, one};
    const T r = vx_call_def_types(x, call, def, types, NULL, 2);
    if (vial_type_is_error(&r)) return err();
    bool same = r.dtype == init.dtype && r.ndim == init.ndim && r.unit_kind == init.unit_kind &&
                (r.unit_kind != VIAL_UNIT_KNOWN || md_unit_base_equal(r.unit, init.unit) || (md_unit_is_none(r.unit) && md_unit_is_none(init.unit)));
    for (uint8_t k = 0; same && k < r.ndim; ++k) {
        vial_extent_t e;
        same = vx_extent_join(r.dims[k], init.dims[k], &e) && !(r.dims[k].kind == VIAL_EXT_STATIC && init.dims[k].kind == VIAL_EXT_STATIC &&
               r.dims[k].n != init.dims[k].n);
    }
    if (!same && !(r.dtype == VIAL_DT_INT && init.dtype == VIAL_DT_INT)) {
        T rr = r, ii = init;
        rr.flags = 0; ii.flags = 0;
        vx_diag_fix(x, "E0315", call, NULL, "the step returns %s; init is %s, and a fold's state keeps its type", vx_tstr(&rr), vx_tstr(&ii));
        return err();
    }
    T t = init;
    t.stage = VIAL_STAGE_TRAJECTORY;
    t.flags |= VIAL_TYPE_PER_FRAME;
    t.flags |= r.flags & VIAL_TYPE_MISSING;
    const uint8_t pv = vx_pass(&v), pi = vx_pass(&init);
    t.pass = pv > pi ? pv : pi;
    return t;
}

bool vx_call_frames2(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out);

bool vx_call_frames(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    if (str_eq(name, STR_LIT("fold"))) { *out = call_fold(x, call); return true; }
    if (vx_call_frames2(x, call, b, out)) return true;
    if (str_eq(name, STR_LIT("frame")) || str_eq(name, STR_LIT("time"))) {
        const bool time = str_eq(name, STR_LIT("time"));
        if (vx_num_children(x, call) > 1) {
            for (uint32_t i = 1; i < vx_num_children(x, call); ++i) if (vx_child(x, call, i)) vx_infer(x, vx_child(x, call, i));
            vx_diag_fix(x, "E0204", call, NULL, "`%s` takes no arguments", b->name);
            *out = err();
            return true;
        }
        T t = vial_type_scalar(time ? VIAL_DT_FLOAT64 : VIAL_DT_INT);
        if (time) vial_unit_name(STR_LIT("ps"), &t.unit);
        t.stage = VIAL_STAGE_FRAME;
        *out = t;
        return true;
    }
    return false;
}

// ### OVER FRAMES: smooth, derivative, histogram, autocorrelation, at ###

static bool static_int(infer_t* x, vial_node_id n, const char* what, int64_t lo, int64_t* out) {
    const vial_value_t v = vx_value_of(x, n);
    if (v.kind != VIAL_VAL_INT || v.i < lo) {
        vx_diag_fix(x, "E0313", n, NULL, "%s is an int of at least %lld, known when compiling", what, (long long)lo);
        return false;
    }
    *out = v.i;
    return true;
}

static bool needs_frames(infer_t* x, vial_node_id n, const T* t, const char* fname) {
    if (t->stage != VIAL_STAGE_STATIC) return true;
    vx_diag_fix(x, "E0313", n, NULL, "`%s` runs over frames, and %s is static", fname, vx_tstr(t));
    return false;
}

// Does a subtree call at()
static bool calls_at(infer_t* x, vial_node_id root) {
    vial_node_id* stack = NULL;
    md_array_push(stack, root, x->alloc);
    while (md_array_size(stack)) {
        const vial_node_id n = *md_array_last(stack);
        md_array_pop(stack);
        if (callee_is(x, n, "at")) return true;
        for (uint32_t i = 0; i < vx_num_children(x, n); ++i) if (vx_child(x, n, i)) md_array_push(stack, vx_child(x, n, i), x->alloc);
    }
    return false;
}

static T call_over_frames(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    if (str_eq(name, STR_LIT("smooth")) || str_eq(name, STR_LIT("derivative"))) {
        const bool smooth = str_eq(name, STR_LIT("smooth"));
        static const vx_param_t ps[] = {{"x", VX_KIND_WHOLE, true}, {"n", VX_KIND_WHOLE, true}};
        static const vx_param_t pd[] = {{"x", VX_KIND_WHOLE, true}};
        vial_node_id a[2] = {0};
        if (!vx_bind(x, call, fname, smooth ? ps : pd, smooth ? 2 : 1, a)) return err();
        const T v = x->out->types[(a[0]) - x->out->base];
        if (vial_type_is_error(&v) || !needs_frames(x, a[0], &v, fname)) return err();
        if (smooth) {
            int64_t n;
            if (!static_int(x, a[1], "the window", 1, &n)) return err();
            if (n % 2 == 0) { vx_diag_fix(x, "E0313", a[1], NULL, "a centred window has an odd number of frames; this is %lld", (long long)n); return err(); }
        }
        T t = v;
        t.stage = VIAL_STAGE_TRAJECTORY;
        t.flags |= VIAL_TYPE_PER_FRAME;
        t.pass = vx_pass(&v);
        if (!smooth) {
            md_unit_t ps_unit;
            vial_unit_name(STR_LIT("ps"), &ps_unit);
            if (v.unit_kind == VIAL_UNIT_KNOWN) t.unit = md_unit_div(v.unit, ps_unit);
            t.flags &= ~VIAL_TYPE_PERIODIC;
            if (t.dtype != VIAL_DT_FLOAT64) t.dtype = VIAL_DT_FLOAT;
        }
        return t;
    }
    if (str_eq(name, STR_LIT("histogram"))) {
        static const vx_param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"axis", VX_KIND_WHOLE, false}, {"bins", VX_KIND_WHOLE, true},
                                       {"range", VX_KIND_WHOLE, true}, {"density", VX_KIND_WHOLE, false}};
        vial_node_id a[5];
        if (!vx_bind(x, call, fname, p, 5, a)) return err();
        const T v = x->out->types[(a[0]) - x->out->base];
        if (vial_type_is_error(&v)) return err();
        if (a[1]) {
            const vial_value_t av = vx_value_of(x, a[1]);
            if (av.kind != VIAL_VAL_STR || !str_eq(av.s, STR_LIT("frame"))) {
                vx_diag_fix(x, "E0316", a[1], NULL, "`histogram` bins over frames, axis='frame'");
                return err();
            }
        }
        if (!needs_frames(x, a[0], &v, fname)) return err();
        int64_t bins;
        if (!static_int(x, a[2], "bins", 1, &bins)) return err();
        const T r = x->out->types[(a[3]) - x->out->base];
        if (vial_type_is_error(&r)) return err();
        if (vx_kind(x, a[3]) != VIAL_N_TUPLE || vx_num_children(x, a[3]) != 2) {
            vx_diag_fix(x, "E0313", a[3], NULL, "range= is a tuple of two numbers, (lo, hi)");
            return err();
        }
        T t = vial_type_scalar(VIAL_DT_FLOAT);
        vial_type_push_dim(&t, vial_extent_static((uint32_t)bins));
        t.stage = VIAL_STAGE_TRAJECTORY;
        t.pass = vx_pass(&v);
        return t;
    }
    if (str_eq(name, STR_LIT("autocorrelation"))) {
        static const vx_param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"max_lag", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a)) return err();
        const T v = x->out->types[(a[0]) - x->out->base];
        if (vial_type_is_error(&v) || !needs_frames(x, a[0], &v, fname)) return err();
        int64_t lag;
        if (!static_int(x, a[1], "max_lag", 0, &lag)) return err();
        T t = v;
        t.dtype = VIAL_DT_FLOAT;
        t.unit_kind = VIAL_UNIT_KNOWN; t.unit = md_unit_none();
        t.flags &= ~(VIAL_TYPE_PERIODIC | VIAL_TYPE_PER_FRAME | VIAL_TYPE_MISSING);
        vial_type_push_dim(&t, vial_extent_static((uint32_t)lag + 1));
        t.stage = VIAL_STAGE_TRAJECTORY;
        t.pass = vx_pass(&v);
        return t;
    }
    if (str_eq(name, STR_LIT("at"))) {
        static const vx_param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"frame", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a)) return err();
        const T v = x->out->types[(a[0]) - x->out->base];
        if (vial_type_is_error(&v)) return err();
        if (v.stage == VIAL_STAGE_TRAJECTORY || calls_at(x, a[0])) {
            vx_diag_fix(x, "E0404", a[0], NULL, v.stage == VIAL_STAGE_TRAJECTORY ? "at() reads one frame, and %s depends on other frames"
                : "at() reads one frame, and this reads another with at()", vx_tstr(&v));
            return err();
        }
        int64_t k;
        if (!static_int(x, a[1], "frame=", 0, &k)) return err();
        if (x->c->num_frames && (uint64_t)k >= x->c->num_frames) {
            vx_diag_fix(x, "E0405", a[1], NULL, "frame %lld of a run with %llu frames", (long long)k, (unsigned long long)x->c->num_frames);
            return err();
        }
        T t = v;
        t.stage = VIAL_STAGE_STATIC;
        t.pass = 0;
        t.flags &= ~VIAL_TYPE_PER_FRAME;
        for (uint8_t q = 0; q < t.ndim; ++q) {
            // An extent that changes per frame is the one frame k has: known when compiling
            if (t.dims[q].kind == VIAL_EXT_DYNAMIC) t.dims[q] = (vial_extent_t){.kind = VIAL_EXT_UNKNOWN, .id = (uint32_t)call, .desc = t.dims[q].desc};
        }
        if (t.dtype == VIAL_DT_SET) t.ident = 0;
        return t;
    }
    return err();
}

bool vx_call_frames2(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    static const char* names[] = {"smooth", "derivative", "histogram", "autocorrelation", "at"};
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        if (strcmp(b->name, names[i]) == 0) { *out = call_over_frames(x, call, str_from_cstr(b->name), b->name); return true; }
    }
    return false;
}

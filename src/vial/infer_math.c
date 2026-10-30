// The checker, step 6: reductions, math, and building, reshaping and joining arrays.
//
// Reductions remove dimensions, never the frame unless axis='frame' names it; by= reduces in groups of particles.
// Math works elementwise with the reference's unit rules. Shapes are values: shape(a) remembers a's extents and what
// each runs over, so zeros(shape(res_a)) carries the residue labels.

#include "infer_internal.h"

#include <md_system.h>

#include <core/md_common.h>
#include <core/md_array.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef vial_type_t T;
typedef vx_param_t param_t;

static T err(void) { return vial_type_error(); }

static bool is_float(uint8_t d) { return d == VIAL_DT_FLOAT || d == VIAL_DT_FLOAT64; }
static bool is_number(const T* t) { return t->dtype == VIAL_DT_INT || is_float(t->dtype) || t->dtype == VIAL_DT_BOOL; }
static bool has_unit(const T* t) { return t->unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(t->unit); }

static void join_stage(T* out, const T* a) {
    if (a->stage > out->stage) out->stage = a->stage;
    if (a->pass > out->pass) out->pass = a->pass;
    out->flags |= a->flags & (VIAL_TYPE_MISSING | VIAL_TYPE_PER_FRAME);
}

static void set_unit(T* t, md_unit_t u) { t->unit_kind = VIAL_UNIT_KNOWN; t->unit = u; }
static void no_unit(T* t) { set_unit(t, md_unit_none()); }

static md_unit_t unit_named(const char* name) {
    md_unit_t u = md_unit_none();
    vial_unit_name(str_from_cstr(name), &u);
    return u;
}

static vial_value_t value(infer_t* x, vial_node_id n) { return n ? vx_value_of(x, n) : (vial_value_t){0}; }
static T type_of(infer_t* x, vial_node_id n) { return n ? x->out->types[n] : err(); }

static const char* tstr(const T* t) { return vx_tstr(t); }

// The square root of a unit, when every exponent is even
static bool unit_sqrt(md_unit_t u, md_unit_t* out) {
    md_unit_t r = u;
    for (int i = 0; i < 8; ++i) {
        if (u.base.arr[i] % 2) return false;
        r.base.arr[i] = (int8_t)(u.base.arr[i] / 2);
    }
    r.mult = sqrt(md_unit_scale(u));
    *out = r;
    return true;
}

static bool is_angle_unit(const T* t) {
    if (t->unit_kind != VIAL_UNIT_KNOWN) return t->unit_kind == VIAL_UNIT_OPEN;
    if (md_unit_is_none(t->unit)) return true;      // a bare number is in radians, as in Python
    md_unit_base_t b = t->unit.base;
    const int8_t angle = b.dim.angle;
    b.dim.angle = 0;
    return angle == 1 && b.raw_bits == 0;
}

// ### AXES ###

#define AXIS_FRAME 0x40     // bit set: the frame; low bits: dimensions

typedef struct axes_t {
    bool     all;           // no axis= given: every dimension
    bool     frame;
    uint8_t  mask;          // dimensions named
} axes_t;

// axis= as an int, negative from the end, 'frame', or a tuple of them; E0316 for one the array does not have
static bool read_axes(infer_t* x, vial_node_id node, const T* a, axes_t* out) {
    memset(out, 0, sizeof(*out));
    if (!node) { out->all = true; return true; }
    vial_node_id items[8];
    uint32_t n = 0;
    if (vx_kind(x, node) == VIAL_N_TUPLE) {
        for (uint32_t i = 0; i < vx_num_children(x, node) && n < 8; ++i) items[n++] = vx_child(x, node, i);
    } else {
        items[n++] = node;
    }
    for (uint32_t i = 0; i < n; ++i) {
        const vial_value_t v = value(x, items[i]);
        if (v.kind == VIAL_VAL_STR && str_eq(v.s, STR_LIT("frame"))) { out->frame = true; continue; }
        if (v.kind == VIAL_VAL_INT) {
            int64_t k = v.i < 0 ? v.i + a->ndim : v.i;
            if (k >= 0 && k < a->ndim) { out->mask |= (uint8_t)(1u << k); continue; }
            char fix[32] = "";
            if (a->ndim) snprintf(fix, sizeof(fix), "axis=%u", a->ndim - 1);
            vx_diag_fix(x, "E0316", items[i], fix[0] ? fix : NULL, "axis=%lld on `%s`, which has %u dimension%s", (long long)v.i, tstr(a), a->ndim,
                a->ndim == 1 ? "" : "s");
            return false;
        }
        const T it = type_of(x, items[i]);
        if (vial_type_is_error(&it)) return false;
        vx_diag_fix(x, "E0316", items[i], NULL, "an axis is an int known when compiling, or 'frame'; this is %s", tstr(&it));
        return false;
    }
    return true;
}

// The dimensions an axes_t keeps
static void remove_axes(T* t, const T* a, const axes_t* ax) {
    t->ndim = 0;
    if (ax->all) return;
    for (uint8_t k = 0; k < a->ndim; ++k) if (!((ax->mask >> k) & 1)) vial_type_push_dim(t, a->dims[k]);
}

// ### REDUCTIONS ###

typedef enum red_t { RED_SUM, RED_MEAN, RED_MIN, RED_MAX, RED_STD, RED_VAR, RED_ANY, RED_ALL, RED_ARGMIN, RED_ARGMAX } red_t;

static int reduction_of(str_t name) {
    static const char* names[] = {"sum", "mean", "min", "max", "std", "var", "any", "all", "argmin", "argmax"};
    for (int i = 0; i < (int)ARRAY_SIZE(names); ++i) if (str_eq(name, str_from_cstr(names[i]))) return i;
    return -1;
}

// any and all of sets: the union and the intersection of an array of sets
static T reduce_sets(infer_t* x, vial_node_id call, vial_node_id an, const T* a, bool any) {
    T t = vial_type_scalar(VIAL_DT_SET);
    join_stage(&t, a);
    t.ident = (uint32_t)call;
    if (a->ndim == 0) { vx_set_value(x, call, value(x, an)); t.flags = a->flags; return t; }
    const vial_value_t v = value(x, an);
    if (v.kind == VIAL_VAL_SETS && a->stage == VIAL_STAGE_STATIC && v.sets->count) {
        const vial_sets_t* s = v.sets;
        vial_sets_t* r;
        if (any) {
            r = vial_sets_merge(s, x->alloc);
        } else {
            r = vial_sets_one(vial_sets_items(s, 0), vial_sets_size(s, 0), x->alloc);
            for (uint32_t i = 1; i < s->count; ++i) r = vial_sets_and(r, 0, s, i, x->alloc);
        }
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_SETS, .sets = r});
        const uint32_t n = vial_sets_size(r, 0);
        if (n == 1) t.flags |= VIAL_TYPE_ONE_PARTICLE;
        if (n == 0) t.flags |= VIAL_TYPE_EMPTY;
    }
    return t;
}

// The first dimension of a that describes particles, for by= and set indices; -1 when there is none
static int particle_dim(infer_t* x, const T* a, const vial_desc_t** desc) {
    for (uint8_t k = 0; k < a->ndim; ++k) {
        const uint32_t d = a->dims[k].desc;
        if (d && d < md_array_size(x->out->descs) && x->out->descs[d].level == VIAL_LEVEL_PARTICLE) {
            if (desc) *desc = &x->out->descs[d];
            return k;
        }
    }
    return -1;
}

static T reduce_by(infer_t* x, vial_node_id call, int red, vial_node_id an, const T* a, vial_node_id gn) {
    const T g = type_of(x, gn);
    if (vial_type_is_error(&g)) return err();
    if (g.dtype != VIAL_DT_SET || g.ndim != 1) {
        vx_diag_fix(x, "E0313", gn, NULL, "`by` takes an array of sets; this is %s", tstr(&g));
        return err();
    }
    const vial_desc_t* d = NULL;
    const int pd = particle_dim(x, a, &d);
    if (pd < 0) {
        vx_diag_fix(x, "E0311", an, NULL, "by= groups particles, and `%s` has no dimension over particles", tstr(a));
        return err();
    }
    T t = *a;
    t.ndim = 0;
    vial_type_push_dim(&t, g.dims[0]);
    join_stage(&t, &g);
    // Groups with none of a's particles give a missing value
    const vial_value_t gv = value(x, gn);
    if (d && d->sets && d->sets->single && gv.kind == VIAL_VAL_SETS && !gv.sets->single && red != RED_SUM && red != RED_ANY && red != RED_ALL) {
        uint32_t empty = 0;
        for (uint32_t i = 0; i < gv.sets->count; ++i) {
            const vial_sets_t* both = vial_sets_and(d->sets, 0, gv.sets, i, x->alloc);
            empty += vial_sets_size(both, 0) == 0;
        }
        if (empty) {
            char ex[VIAL_EXCERPT_MAX + 4];
            const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, an));
            vx_diag_fix(x, "W0603", gn, NULL, "%u of %u groups have no values in `%.*s`", empty, gv.sets->count, (int)e.len, e.ptr);
            vx_diag_fix(x, "W0405", call, NULL, "the groups with no values are missing");
            t.flags |= VIAL_TYPE_MISSING;
        }
    }
    return t;
}

static T call_reduce(infer_t* x, vial_node_id call, int red, const char* fname) {
    // min(a, b) and max(a, b) with two arguments work elementwise
    uint32_t num_pos = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id c = vx_child(x, call, i);
        num_pos += c && vx_kind(x, c) != VIAL_N_KWARG;
    }
    if ((red == RED_MIN || red == RED_MAX) && num_pos == 2) {
        static const param_t two[] = {{"a", VX_KIND_WHOLE, true}, {"b", VX_KIND_WHOLE, true}};
        vial_node_id args[2];
        if (!vx_bind(x, call, fname, two, 2, args)) return err();
        const T a = type_of(x, args[0]), b = type_of(x, args[1]);
        if (vial_type_is_error(&a) || vial_type_is_error(&b)) return err();
        if (!is_number(&a) || !is_number(&b)) {
            vx_diag_fix(x, "E0313", call, NULL, "`%s` of two values takes numbers; this is %s and %s", fname, tstr(&a), tstr(&b));
            return err();
        }
        T t = vial_type_scalar(vx_promote(&a, &b));
        if (!vx_broadcast(&a, &b, &t)) {
            vx_diag_fix(x, "E0301", call, NULL, "`%s` and `%s` do not broadcast", tstr(&a), tstr(&b));
            return err();
        }
        vx_unit_join(x, call, "compares", &a, &b, &t);
        join_stage(&t, &a);
        join_stage(&t, &b);
        return t;
    }
    static const param_t params[] = {{"a", VX_KIND_WHOLE, true}, {"axis", VX_KIND_WHOLE, false}, {"by", VX_KIND_WHOLE, false}};
    vial_node_id args[3];
    if (!vx_bind(x, call, fname, params, 3, args)) return err();
    const T a = type_of(x, args[0]);
    if (vial_type_is_error(&a)) return err();
    if ((red == RED_ANY || red == RED_ALL) && a.dtype == VIAL_DT_SET && !args[1] && !args[2]) return reduce_sets(x, call, args[0], &a, red == RED_ANY);
    if (!is_number(&a)) {
        vx_diag_fix(x, "E0313", args[0], NULL, "`%s` takes numbers or bools; this is %s", fname, tstr(&a));
        return err();
    }
    T t;
    if (args[2]) {
        if (args[1]) {
            vx_diag_fix(x, "E0316", args[1], NULL, "`%s` takes axis= or by=, not both", fname);
            return err();
        }
        t = reduce_by(x, call, red, args[0], &a, args[2]);
        if (vial_type_is_error(&t)) return t;
    } else {
        axes_t ax;
        if (args[1] && vial_type_is_error(&x->out->types[args[1]])) return err();
        if (!read_axes(x, args[1], &a, &ax)) return err();
        t = a;
        remove_axes(&t, &a, &ax);
        if (ax.frame) {
            if (a.stage == VIAL_STAGE_STATIC) {
                vx_diag_fix(x, "E0316", args[1], NULL, "axis='frame' on a static value, which has no frames");
                return err();
            }
            t.stage = VIAL_STAGE_TRAJECTORY;
            t.pass = (uint8_t)(a.pass + 1);
            t.flags &= ~VIAL_TYPE_PER_FRAME;
        }
        // A reduction over an extent that is empty in some frame has no value there, where it has no identity
        bool dynamic = false;
        for (uint8_t k = 0; k < a.ndim; ++k) if ((ax.all || ((ax.mask >> k) & 1)) && a.dims[k].kind == VIAL_EXT_DYNAMIC) dynamic = true;
        if (dynamic && red != RED_SUM && red != RED_ANY && red != RED_ALL) {
            vx_diag_fix(x, "W0405", call, NULL, "`%s` over an extent that may be empty is missing in the frames where it is", fname);
            t.flags |= VIAL_TYPE_MISSING;
        }
    }
    // Dtype and unit
    const bool periodic = a.flags & VIAL_TYPE_PERIODIC;
    t.flags &= ~VIAL_TYPE_PERIODIC;
    switch (red) {
    case RED_SUM:
        t.dtype = a.dtype == VIAL_DT_BOOL ? VIAL_DT_INT : a.dtype;
        break;
    case RED_MEAN:
        t.dtype = a.dtype == VIAL_DT_FLOAT64 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT;
        if (periodic) { t.flags |= VIAL_TYPE_PERIODIC; t.period_beg = a.period_beg; t.period_end = a.period_end; t.seam = a.seam; }
        break;
    case RED_MIN: case RED_MAX:
        if (periodic) { t.flags |= VIAL_TYPE_PERIODIC; t.period_beg = a.period_beg; t.period_end = a.period_end; t.seam = a.seam; }
        break;
    case RED_STD:
        t.dtype = a.dtype == VIAL_DT_FLOAT64 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT;
        break;
    case RED_VAR:
        t.dtype = a.dtype == VIAL_DT_FLOAT64 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT;
        if (has_unit(&a)) set_unit(&t, md_unit_mul(a.unit, a.unit));
        break;
    case RED_ANY: case RED_ALL:
        t.dtype = VIAL_DT_BOOL;
        no_unit(&t);
        break;
    case RED_ARGMIN: case RED_ARGMAX:
        t.dtype = VIAL_DT_INT;
        no_unit(&t);
        break;
    default: break;
    }
    return t;
}

// ### MATH ###

static T call_unary(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    static const param_t params[] = {{"x", VX_KIND_WHOLE, true}};
    vial_node_id arg;
    if (!vx_bind(x, call, fname, params, 1, &arg)) return err();
    const T a = type_of(x, arg);
    if (vial_type_is_error(&a)) return err();
    if (!is_number(&a)) {
        vx_diag_fix(x, "E0313", arg, NULL, "`%s` takes numbers; this is %s", fname, tstr(&a));
        return err();
    }
    T t = a;
    t.flags &= ~VIAL_TYPE_PERIODIC;
    if (str_eq(name, STR_LIT("abs"))) {
        if (t.dtype == VIAL_DT_BOOL) t.dtype = VIAL_DT_INT;
        return t;
    }
    if (str_eq(name, STR_LIT("floor")) || str_eq(name, STR_LIT("ceil")) || str_eq(name, STR_LIT("round"))) return t;
    if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT;
    if (str_eq(name, STR_LIT("sqrt"))) {
        if (has_unit(&a)) {
            md_unit_t u;
            if (unit_sqrt(a.unit, &u)) set_unit(&t, u);
            else {
                vx_diag_fix(x, "W0501", call, "unitless(...)", "the square root of %s has no unit VIAL can write; the result's unit becomes unknown",
                    tstr(&a));
                t.unit_kind = VIAL_UNIT_UNKNOWN;
            }
        }
        return t;
    }
    const bool trig = str_eq(name, STR_LIT("sin")) || str_eq(name, STR_LIT("cos")) || str_eq(name, STR_LIT("tan"));
    if (trig) {
        if (a.unit_kind == VIAL_UNIT_KNOWN && !is_angle_unit(&a)) {
            vx_diag_fix(x, "W0501", call, "unitless(...)", "`%s` takes an angle; this is %s", fname, tstr(&a));
        }
        no_unit(&t);
        return t;
    }
    const bool inverse = str_eq(name, STR_LIT("asin")) || str_eq(name, STR_LIT("acos")) || str_eq(name, STR_LIT("atan"));
    if (has_unit(&a)) {
        vx_diag_fix(x, "W0501", call, "unitless(...)", "`%s` takes a value without a unit; this is %s", fname, tstr(&a));
    }
    if (inverse) set_unit(&t, unit_named("deg"));
    else no_unit(&t);
    return t;
}

// Several values broadcast together, the first's unit for the rest
static bool broadcast_all(infer_t* x, vial_node_id call, const vial_node_id* args, uint32_t n, T* out) {
    T t = type_of(x, args[0]);
    for (uint32_t i = 1; i < n; ++i) {
        if (!args[i]) continue;
        const T b = type_of(x, args[i]);
        T m = t;
        if (!vx_broadcast(&t, &b, &m)) {
            vx_diag_fix(x, "E0301", call, NULL, "`%s` and `%s` do not broadcast", tstr(&t), tstr(&b));
            return false;
        }
        t.ndim = m.ndim;
        memcpy(t.dims, m.dims, sizeof(m.dims));
        join_stage(&t, &b);
    }
    *out = t;
    return true;
}

static bool all_typed(infer_t* x, const vial_node_id* args, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) if (args[i] && vial_type_is_error(&x->out->types[args[i]])) return false;
    return true;
}

static bool all_numbers(infer_t* x, const char* fname, const vial_node_id* args, uint32_t n) {
    for (uint32_t i = 0; i < n; ++i) {
        if (!args[i]) continue;
        const T a = type_of(x, args[i]);
        if (!is_number(&a)) {
            vx_diag_fix(x, "E0313", args[i], NULL, "`%s` takes numbers; this is %s", fname, tstr(&a));
            return false;
        }
    }
    return true;
}

static bool last_is(const T* t, uint32_t n) {
    return t->ndim >= 1 && t->dims[t->ndim - 1].kind == VIAL_EXT_STATIC && t->dims[t->ndim - 1].n == n;
}

static T call_math(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    if (str_eq(name, STR_LIT("clip"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"lo", VX_KIND_WHOLE, true}, {"hi", VX_KIND_WHOLE, true}};
        vial_node_id a[3];
        if (!vx_bind(x, call, fname, p, 3, a) || !all_typed(x, a, 3) || !all_numbers(x, fname, a, 3)) return err();
        T t;
        if (!broadcast_all(x, call, a, 3, &t)) return err();
        const T v = type_of(x, a[0]), lo = type_of(x, a[1]), hi = type_of(x, a[2]);
        T u = t;
        vx_unit_join(x, call, "clips", &v, &lo, &u);
        vx_unit_join(x, call, "clips", &v, &hi, &u);
        t.unit_kind = u.unit_kind; t.unit = u.unit;
        return t;
    }
    if (str_eq(name, STR_LIT("atan2"))) {
        static const param_t p[] = {{"y", VX_KIND_WHOLE, true}, {"x", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2) || !all_numbers(x, fname, a, 2)) return err();
        T t;
        if (!broadcast_all(x, call, a, 2, &t)) return err();
        const T y = type_of(x, a[0]), xx = type_of(x, a[1]);
        T u = t;
        vx_unit_join(x, call, "compares", &y, &xx, &u);
        if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT;
        set_unit(&t, unit_named("deg"));
        return t;
    }
    if (str_eq(name, STR_LIT("dot")) || str_eq(name, STR_LIT("cross"))) {
        const bool cross = str_eq(name, STR_LIT("cross"));
        static const param_t p[] = {{"u", VX_KIND_WHOLE, true}, {"v", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2) || !all_numbers(x, fname, a, 2)) return err();
        const T u = type_of(x, a[0]), v = type_of(x, a[1]);
        if (u.ndim == 0 || v.ndim == 0 || (cross && (!last_is(&u, 3) || !last_is(&v, 3)))) {
            vx_diag_fix(x, "E0313", call, NULL, cross ? "`cross` takes vectors, a last dimension of 3; this is %s and %s" :
                "`dot` takes vectors, arrays of one dimension or more; this is %s and %s", tstr(&u), tstr(&v));
            return err();
        }
        T t;
        if (!broadcast_all(x, call, a, 2, &t)) return err();
        t.dtype = vx_promote(&u, &v);
        if (!cross) t.ndim -= 1;
        // Units multiply
        if (u.unit_kind == VIAL_UNIT_UNKNOWN || v.unit_kind == VIAL_UNIT_UNKNOWN) t.unit_kind = VIAL_UNIT_UNKNOWN;
        else set_unit(&t, md_unit_mul(u.unit_kind == VIAL_UNIT_KNOWN ? u.unit : md_unit_none(), v.unit_kind == VIAL_UNIT_KNOWN ? v.unit : md_unit_none()));
        t.flags &= ~VIAL_TYPE_PERIODIC;
        return t;
    }
    if (str_eq(name, STR_LIT("norm")) || str_eq(name, STR_LIT("normalize"))) {
        const bool normalize = str_eq(name, STR_LIT("normalize"));
        static const param_t p[] = {{"v", VX_KIND_WHOLE, true}, {"axis", VX_KIND_WHOLE, false}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, normalize ? 1 : 2, a) || !all_typed(x, a, normalize ? 1 : 2) || !all_numbers(x, fname, a, 1)) return err();
        const T v = type_of(x, a[0]);
        if (v.ndim == 0) {
            vx_diag_fix(x, "E0313", a[0], NULL, "`%s` takes a vector, an array; this is %s", fname, tstr(&v));
            return err();
        }
        T t = v;
        if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT;
        t.flags &= ~VIAL_TYPE_PERIODIC;
        if (normalize) { no_unit(&t); return t; }
        int64_t axis = -1;
        if (a[1]) {
            const vial_value_t av = value(x, a[1]);
            if (av.kind != VIAL_VAL_INT) {
                vx_diag_fix(x, "E0316", a[1], NULL, "`norm` takes one axis, an int known when compiling");
                return err();
            }
            axis = av.i;
        }
        const int64_t k = axis < 0 ? axis + v.ndim : axis;
        if (k < 0 || k >= v.ndim) {
            vx_diag_fix(x, "E0316", a[1] ? a[1] : call, NULL, "axis=%lld on `%s`, which has %u dimensions", (long long)axis, tstr(&v), v.ndim);
            return err();
        }
        t.ndim = 0;
        for (uint8_t q = 0; q < v.ndim; ++q) if (q != k) vial_type_push_dim(&t, v.dims[q]);
        return t;
    }
    if (str_eq(name, STR_LIT("where"))) {
        static const param_t p[] = {{"c", VX_KIND_WHOLE, true}, {"a", VX_KIND_WHOLE, true}, {"b", VX_KIND_WHOLE, true}};
        vial_node_id a[3];
        if (!vx_bind(x, call, fname, p, 3, a) || !all_typed(x, a, 3)) return err();
        const T c = type_of(x, a[0]), va = type_of(x, a[1]), vb = type_of(x, a[2]);
        if (c.dtype != VIAL_DT_BOOL) {
            vx_diag_fix(x, "E0313", a[0], NULL, "`where` takes a bool condition; this is %s", tstr(&c));
            return err();
        }
        if (!all_numbers(x, fname, a + 1, 2)) return err();
        T t;
        if (!broadcast_all(x, call, a, 3, &t)) return err();
        t.dtype = vx_promote(&va, &vb);
        T u = t;
        vx_unit_join(x, call, "joins", &va, &vb, &u);
        t.unit_kind = u.unit_kind; t.unit = u.unit;
        t.flags &= ~VIAL_TYPE_PERIODIC;
        return t;
    }
    if (str_eq(name, STR_LIT("unitless")) || str_eq(name, STR_LIT("float")) || str_eq(name, STR_LIT("int"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1) || !all_numbers(x, fname, &a, 1)) return err();
        T t = type_of(x, a);
        if (str_eq(name, STR_LIT("unitless"))) { no_unit(&t); t.flags &= ~VIAL_TYPE_PERIODIC; }
        else if (str_eq(name, STR_LIT("float"))) { if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT; }
        else t.dtype = VIAL_DT_INT;
        const vial_value_t v = value(x, a);
        if (str_eq(name, STR_LIT("int")) && v.kind == VIAL_VAL_FLOAT) vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_INT, .i = (int64_t)v.f});
        else if (str_eq(name, STR_LIT("float")) && v.kind == VIAL_VAL_INT) vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_FLOAT, .f = (double)v.i});
        else vx_set_value(x, call, v);
        return t;
    }
    return call_unary(x, call, name, fname);
}

// ### SHAPES AND BUILDING ###

// A shape argument: an int, a tuple of ints, or shape(a); its extents into t. E0313 otherwise.
static bool shape_from(infer_t* x, vial_node_id node, T* t) {
    const T s = type_of(x, node);
    const vial_value_t v = value(x, node);
    t->ndim = 0;
    if (v.kind == VIAL_VAL_SHAPE && v.type) {
        for (uint8_t k = 0; k < v.type->ndim; ++k) vial_type_push_dim(t, v.type->dims[k]);
        join_stage(t, &s);
        return true;
    }
    if (s.dtype == VIAL_DT_INT && s.ndim == 0) {
        if (v.kind != VIAL_VAL_INT || v.i < 0) {
            vx_diag_fix(x, "E0313", node, NULL, "a shape is known when compiling, an int or a tuple of ints");
            return false;
        }
        vial_type_push_dim(t, vial_extent_static((uint32_t)v.i));
        return true;
    }
    if (vx_kind(x, node) == VIAL_N_TUPLE) {
        for (uint32_t i = 0; i < vx_num_children(x, node); ++i) {
            const vial_value_t c = value(x, vx_child(x, node, i));
            if (c.kind != VIAL_VAL_INT || c.i < 0) {
                vx_diag_fix(x, "E0313", node, NULL, "a shape is known when compiling, an int or a tuple of ints");
                return false;
            }
            if (!vial_type_push_dim(t, vial_extent_static((uint32_t)c.i))) {
                vx_diag_fix(x, "E0308", node, NULL, "the result would have 7 dimensions; the limit is 6");
                return false;
            }
        }
        return true;
    }
    if (!vial_type_is_error(&s)) vx_diag_fix(x, "E0313", node, NULL, "a shape is an int, a tuple of ints or shape(a); this is %s", tstr(&s));
    return false;
}

// The number of elements of a type, when every extent is static
static bool elements(const T* t, uint64_t* n) {
    *n = 1;
    for (uint8_t k = 0; k < t->ndim; ++k) {
        if (t->dims[k].kind != VIAL_EXT_STATIC) return false;
        *n *= t->dims[k].n;
    }
    return true;
}

static T call_build(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    if (str_eq(name, STR_LIT("shape"))) {
        static const param_t p[] = {{"a", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1)) return err();
        const T at = type_of(x, a);
        T* keep = md_alloc(x->alloc, sizeof(T));
        *keep = at;
        T* elems = md_alloc(x->alloc, (at.ndim + 1) * sizeof(T));
        for (uint8_t k = 0; k < at.ndim; ++k) {
            elems[k] = vial_type_scalar(VIAL_DT_INT);
            if (at.dims[k].kind == VIAL_EXT_DYNAMIC) elems[k].stage = at.stage;
        }
        T t = vial_type_scalar(VIAL_DT_TUPLE);
        t.elems = elems;
        t.num_elems = at.ndim;
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_SHAPE, .type = keep});
        return t;
    }
    if (str_eq(name, STR_LIT("zeros")) || str_eq(name, STR_LIT("ones")) || str_eq(name, STR_LIT("full"))) {
        const bool full = str_eq(name, STR_LIT("full"));
        static const param_t pz[] = {{"shape", VX_KIND_WHOLE, true}, {"unit", VX_KIND_WHOLE, false}};
        static const param_t pf[] = {{"shape", VX_KIND_WHOLE, true}, {"value", VX_KIND_WHOLE, true}, {"unit", VX_KIND_WHOLE, false}};
        vial_node_id a[3] = {0};
        if (!vx_bind(x, call, fname, full ? pf : pz, full ? 3 : 2, a) || !all_typed(x, a, full ? 3 : 2)) return err();
        T t = vial_type_scalar(VIAL_DT_FLOAT);
        t.unit_kind = VIAL_UNIT_OPEN;
        if (full) {
            const T v = type_of(x, a[1]);
            if (!is_number(&v) || v.ndim) {
                vx_diag_fix(x, "E0313", a[1], NULL, "`full` fills with one number; this is %s", tstr(&v));
                return err();
            }
            t.unit_kind = v.unit_kind; t.unit = v.unit;
            if (is_float(v.dtype)) t.dtype = v.dtype;
            join_stage(&t, &v);
        }
        if (!shape_from(x, a[0], &t)) return err();
        const vial_node_id un = a[full ? 2 : 1];
        if (un) {
            const vial_value_t u = value(x, un);
            md_unit_t unit;
            if (u.kind != VIAL_VAL_STR || !vial_unit_string(u.s, &unit)) {
                vx_diag_fix(x, "E0502", un, NULL, "unit= takes a unit string, such as 'nm' or 'kJ/mol'");
                return err();
            }
            set_unit(&t, unit);
        }
        return t;
    }
    if (str_eq(name, STR_LIT("zeros_like")) || str_eq(name, STR_LIT("ones_like"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1) || !all_numbers(x, fname, &a, 1)) return err();
        T t = type_of(x, a);
        t.flags &= ~(VIAL_TYPE_MISSING | VIAL_TYPE_PERIODIC);
        return t;
    }
    if (str_eq(name, STR_LIT("linspace")) || str_eq(name, STR_LIT("arange"))) {
        const bool lin = str_eq(name, STR_LIT("linspace"));
        static const param_t pl[] = {{"lo", VX_KIND_WHOLE, true}, {"hi", VX_KIND_WHOLE, true}, {"n", VX_KIND_WHOLE, true}};
        static const param_t pa[] = {{"lo", VX_KIND_WHOLE, true}, {"hi", VX_KIND_WHOLE, true}, {"step", VX_KIND_WHOLE, false}};
        vial_node_id a[3];
        if (!vx_bind(x, call, fname, lin ? pl : pa, 3, a) || !all_typed(x, a, 3) || !all_numbers(x, fname, a, 3)) return err();
        const T lo = type_of(x, a[0]), hi = type_of(x, a[1]);
        T t = vial_type_scalar(VIAL_DT_FLOAT);
        vx_unit_join(x, call, "joins", &lo, &hi, &t);
        uint32_t n = 0;
        bool known = false;
        if (lin) {
            const vial_value_t v = value(x, a[2]);
            if (v.kind == VIAL_VAL_INT && v.i >= 0) { n = (uint32_t)v.i; known = true; }
        } else {
            const vial_value_t l = value(x, a[0]), h = value(x, a[1]), s = a[2] ? value(x, a[2]) : (vial_value_t){.kind = VIAL_VAL_INT, .i = 1};
            double fl, fh, fs;
            if ((l.kind == VIAL_VAL_INT || l.kind == VIAL_VAL_FLOAT) && (h.kind == VIAL_VAL_INT || h.kind == VIAL_VAL_FLOAT) &&
                (s.kind == VIAL_VAL_INT || s.kind == VIAL_VAL_FLOAT)) {
                fl = l.kind == VIAL_VAL_INT ? (double)l.i : l.f;
                fh = h.kind == VIAL_VAL_INT ? (double)h.i : h.f;
                fs = s.kind == VIAL_VAL_INT ? (double)s.i : s.f;
                if (fs != 0.0) { const double c = ceil((fh - fl) / fs); n = c > 0 ? (uint32_t)c : 0; known = true; }
            }
        }
        if (!known) {
            vx_diag_fix(x, "E0313", call, NULL, "`%s` makes an array whose length is known when compiling", fname);
            return err();
        }
        vial_type_push_dim(&t, vial_extent_static(n));
        return t;
    }
    if (str_eq(name, STR_LIT("eye"))) {
        static const param_t p[] = {{"n", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1)) return err();
        const vial_value_t v = value(x, a);
        if (v.kind != VIAL_VAL_INT || v.i < 0) {
            vx_diag_fix(x, "E0313", a, NULL, "`eye` takes an int known when compiling");
            return err();
        }
        T t = vial_type_scalar(VIAL_DT_FLOAT);
        vial_type_push_dim(&t, vial_extent_static((uint32_t)v.i));
        vial_type_push_dim(&t, vial_extent_static((uint32_t)v.i));
        return t;
    }
    if (str_eq(name, STR_LIT("reshape"))) {
        static const param_t p[] = {{"a", VX_KIND_WHOLE, true}, {"shape", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2)) return err();
        const T src = type_of(x, a[0]);
        T t = src;
        // -1 infers one extent
        T want = vial_type_scalar(VIAL_DT_INT);
        int infer_at = -1;
        if (vx_kind(x, a[1]) == VIAL_N_TUPLE) {
            for (uint32_t i = 0; i < vx_num_children(x, a[1]); ++i) {
                const vial_value_t c = value(x, vx_child(x, a[1], i));
                if (c.kind == VIAL_VAL_INT && c.i == -1 && infer_at < 0) { infer_at = (int)i; vial_type_push_dim(&want, vial_extent_static(1)); continue; }
                if (c.kind != VIAL_VAL_INT || c.i < 0) {
                    vx_diag_fix(x, "E0313", a[1], NULL, "a shape is known when compiling, an int or a tuple of ints");
                    return err();
                }
                vial_type_push_dim(&want, vial_extent_static((uint32_t)c.i));
            }
        } else if (!shape_from(x, a[1], &want)) {
            return err();
        }
        uint64_t have, need;
        if (elements(&src, &have) && elements(&want, &need)) {
            if (infer_at >= 0 && need && have % need == 0) { want.dims[infer_at].n = (uint32_t)(have / need); need = have; }
            if (have != need) {
                vx_diag_fix(x, "E0307", call, NULL, "%llu values do not fill %llu", (unsigned long long)have, (unsigned long long)need);
                return err();
            }
        }
        t.ndim = want.ndim;
        memcpy(t.dims, want.dims, sizeof(want.dims));
        return t;
    }
    if (str_eq(name, STR_LIT("flatten"))) {
        static const param_t p[] = {{"a", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1)) return err();
        const T src = type_of(x, a);
        T t = src;
        t.ndim = 0;
        uint64_t n;
        vial_type_push_dim(&t, elements(&src, &n) && n <= UINT32_MAX ? vial_extent_static((uint32_t)n) :
            (vial_extent_t){.kind = src.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
        return t;
    }
    if (str_eq(name, STR_LIT("diff"))) {
        static const param_t p[] = {{"a", VX_KIND_WHOLE, true}, {"axis", VX_KIND_WHOLE, false}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2) || !all_numbers(x, fname, a, 1)) return err();
        const T src = type_of(x, a[0]);
        int64_t axis = 0;
        if (a[1]) {
            const vial_value_t v = value(x, a[1]);
            if (v.kind != VIAL_VAL_INT) { vx_diag_fix(x, "E0316", a[1], NULL, "`diff` takes one axis, an int known when compiling"); return err(); }
            axis = v.i;
        }
        const int64_t k = axis < 0 ? axis + src.ndim : axis;
        if (k < 0 || k >= src.ndim) {
            vx_diag_fix(x, "E0316", a[1] ? a[1] : a[0], NULL, "axis=%lld on `%s`, which has %u dimensions", (long long)axis, tstr(&src), src.ndim);
            return err();
        }
        T t = src;
        vial_extent_t* e = &t.dims[k];
        if (e->kind == VIAL_EXT_STATIC) e->n = e->n ? e->n - 1 : 0;
        else if (e->kind == VIAL_EXT_DYNAMIC) e->offset += 1;
        if (t.dtype == VIAL_DT_BOOL) t.dtype = VIAL_DT_INT;
        // The difference of two periodic values is cut half a period from 0
        if (t.flags & VIAL_TYPE_PERIODIC) {
            const float w = t.period_end - t.period_beg;
            t.period_beg = -w / 2; t.period_end = w / 2; t.seam = w / 2;
        }
        return t;
    }
    return err();
}

// ### JOINING AND LINEAR ALGEBRA ###

static T call_concatenate(infer_t* x, vial_node_id call) {
    // The list or tuple written in the call is not typed as a list, which would stack its items: its items are
    vial_node_id a[2] = {0, 0};
    bool ok = true;
    uint32_t num_pos = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id c = vx_child(x, call, i);
        if (!c) continue;
        if (vx_kind(x, c) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, c, 0), vn = vx_child(x, c, 1);
            if (vn) vx_infer(x, vn);
            if (kn && str_eq(vx_text(x, kn), STR_LIT("axis")) && !a[1]) { a[1] = vn; continue; }
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            vx_diag_fix(x, "E0202", c, NULL, "`concatenate` takes no keyword `%.*s`; it takes arrays and axis", (int)kw.len, kw.ptr);
            ok = false;
            continue;
        }
        num_pos += 1;
        if (num_pos == 1) {
            a[0] = c;
            if (vx_kind(x, c) == VIAL_N_LIST || vx_kind(x, c) == VIAL_N_TUPLE) {
                for (uint32_t q = 0; q < vx_num_children(x, c); ++q) if (vx_child(x, c, q)) vx_infer(x, vx_child(x, c, q));
            } else {
                vx_infer(x, c);
            }
        } else if (num_pos == 2 && !a[1]) {
            a[1] = c;
            vx_infer(x, c);
        } else {
            vx_infer(x, c);
            vx_diag_fix(x, "E0204", c, NULL, "`concatenate` takes at most 2 arguments");
            ok = false;
        }
    }
    if (!ok) return err();
    if (!a[0]) { vx_diag_fix(x, "E0204", call, NULL, "`concatenate` needs arrays"); return err(); }
    const vial_node_kind_t k = vx_kind(x, a[0]);
    if (k != VIAL_N_LIST && k != VIAL_N_TUPLE) {
        const T it = type_of(x, a[0]);
        if (!vial_type_is_error(&it)) vx_diag_fix(x, "E0313", a[0], NULL, "`concatenate` takes a tuple or a list of arrays written in the call");
        return err();
    }
    int64_t axis = 0;
    if (a[1]) {
        const vial_value_t v = value(x, a[1]);
        if (v.kind != VIAL_VAL_INT) { vx_diag_fix(x, "E0316", a[1], NULL, "`concatenate` takes one axis, an int known when compiling"); return err(); }
        axis = v.i;
    }
    const uint32_t n = vx_num_children(x, a[0]);
    if (n == 0) { vx_diag_fix(x, "E0313", a[0], NULL, "`concatenate` needs at least one array"); return err(); }
    T t = type_of(x, vx_child(x, a[0], 0));
    if (vial_type_is_error(&t)) return err();
    const int64_t ax = axis < 0 ? axis + t.ndim : axis;
    if (t.ndim == 0 || ax < 0 || ax >= t.ndim) {
        vx_diag_fix(x, "E0316", a[1] ? a[1] : a[0], NULL, "`concatenate` joins along a dimension the arrays have; `%s` has %u", tstr(&t), t.ndim);
        return err();
    }
    bool static_n = t.dims[ax].kind == VIAL_EXT_STATIC;
    uint32_t total = t.dims[ax].n;
    vial_sets_t* joined = NULL;
    const vial_value_t v0 = value(x, vx_child(x, a[0], 0));
    bool sets_known = t.dtype == VIAL_DT_SET && t.ndim == 1 && v0.kind == VIAL_VAL_SETS;
    for (uint32_t i = 1; i < n; ++i) {
        const vial_node_id in = vx_child(x, a[0], i);
        const T b = type_of(x, in);
        if (vial_type_is_error(&b)) return err();
        if (b.dtype != t.dtype && !(is_number(&b) && is_number(&t))) {
            vx_diag_fix(x, "E0319", in, NULL, "`concatenate` joins arrays of one kind; this is %s and %s", tstr(&t), tstr(&b));
            return err();
        }
        if (b.ndim != t.ndim) {
            vx_diag_fix(x, "E0309", in, NULL, "`%s` and `%s` have different numbers of dimensions", tstr(&t), tstr(&b));
            return err();
        }
        for (uint8_t q = 0; q < t.ndim; ++q) {
            if (q == ax) continue;
            vial_extent_t e;
            if (!vx_extent_join(t.dims[q], b.dims[q], &e) || (t.dims[q].kind == VIAL_EXT_STATIC && b.dims[q].kind == VIAL_EXT_STATIC && t.dims[q].n != b.dims[q].n)) {
                vx_diag_fix(x, "E0309", in, NULL, "`%s` and `%s` differ in dimension %u", tstr(&t), tstr(&b), q);
                return err();
            }
        }
        static_n = static_n && b.dims[ax].kind == VIAL_EXT_STATIC;
        total += b.dims[ax].n;
        if (is_number(&b)) { t.dtype = vx_promote(&t, &b); T u = t; vx_unit_join(x, in, "joins", &t, &b, &u); t.unit_kind = u.unit_kind; t.unit = u.unit; }
        join_stage(&t, &b);
        const vial_value_t bv = value(x, in);
        sets_known = sets_known && bv.kind == VIAL_VAL_SETS;
    }
    t.dims[ax] = static_n ? vial_extent_static(total) :
        (vial_extent_t){.kind = t.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call};
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    if (sets_known && static_n && t.stage == VIAL_STAGE_STATIC) {
        // The sets one after the other: a new array, its units of no one level
        uint32_t w = 0;
        joined = md_alloc(x->alloc, sizeof(vial_sets_t));
        memset(joined, 0, sizeof(*joined));
        uint32_t idx_total = 0;
        for (uint32_t i = 0; i < n; ++i) { const vial_sets_t* s = value(x, vx_child(x, a[0], i)).sets; idx_total += s->offset[s->count]; }
        joined->offset = md_alloc(x->alloc, (total + 1) * sizeof(uint32_t));
        joined->idx = md_alloc(x->alloc, (idx_total + 1) * sizeof(uint32_t));
        joined->offset[0] = 0;
        uint32_t c = 0;
        for (uint32_t i = 0; i < n; ++i) {
            const vial_sets_t* s = value(x, vx_child(x, a[0], i)).sets;
            for (uint32_t q = 0; q < s->count; ++q) {
                const uint32_t sz = vial_sets_size(s, q);
                memcpy(joined->idx + w, vial_sets_items(s, q), sz * sizeof(uint32_t));
                w += sz;
                joined->offset[++c] = w;
            }
        }
        joined->count = c;
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_SETS, .sets = joined});
        t.dims[ax].desc = vx_desc(x, VIAL_LEVEL_NONE, joined);
    }
    return t;
}

static bool square(infer_t* x, vial_node_id at, const char* fname, const T* m) {
    const bool sq = m->ndim >= 2 && m->dims[m->ndim - 1].kind == m->dims[m->ndim - 2].kind &&
                    (m->dims[m->ndim - 1].kind != VIAL_EXT_STATIC || m->dims[m->ndim - 1].n == m->dims[m->ndim - 2].n);
    if (!sq) vx_diag_fix(x, "E0310", at, NULL, "`%s` needs a square matrix; this is `%s`", fname, tstr(m));
    return sq;
}

static T call_linalg(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    if (str_eq(name, STR_LIT("covariance"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"axis", VX_KIND_WHOLE, false}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2) || !all_numbers(x, fname, a, 1)) return err();
        const T v = type_of(x, a[0]);
        if (a[1]) {
            const vial_value_t av = value(x, a[1]);
            if (av.kind != VIAL_VAL_STR || !str_eq(av.s, STR_LIT("frame"))) {
                vx_diag_fix(x, "E0316", a[1], NULL, "`covariance` is over frames, axis='frame'");
                return err();
            }
        }
        if (v.ndim != 1) {
            vx_diag_fix(x, "E0313", a[0], NULL, "`covariance` takes values of one dimension, [n]; this is %s", tstr(&v));
            return err();
        }
        if (v.stage == VIAL_STAGE_STATIC) {
            vx_diag_fix(x, "E0316", a[1] ? a[1] : a[0], NULL, "covariance over frames of a static value, which has no frames");
            return err();
        }
        T t = v;
        t.dtype = VIAL_DT_FLOAT;
        vial_type_push_dim(&t, v.dims[0]);
        if (has_unit(&v)) set_unit(&t, md_unit_mul(v.unit, v.unit));
        t.stage = VIAL_STAGE_TRAJECTORY;
        t.pass = (uint8_t)(v.pass + 1);
        t.flags &= ~(VIAL_TYPE_PERIODIC | VIAL_TYPE_PER_FRAME);
        return t;
    }
    if (str_eq(name, STR_LIT("eigh")) || str_eq(name, STR_LIT("trace")) || str_eq(name, STR_LIT("det")) || str_eq(name, STR_LIT("inv"))) {
        static const param_t p[] = {{"m", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1) || !all_numbers(x, fname, &a, 1)) return err();
        const T m = type_of(x, a);
        if (!square(x, call, fname, &m)) return err();
        T t = m;
        if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT;
        t.flags &= ~VIAL_TYPE_PERIODIC;
        if (str_eq(name, STR_LIT("trace"))) { t.ndim -= 2; return t; }
        if (str_eq(name, STR_LIT("det"))) {
            t.ndim -= 2;
            if (has_unit(&m)) {
                if (m.dims[m.ndim - 1].kind == VIAL_EXT_STATIC) set_unit(&t, md_unit_pow(m.unit, (int)m.dims[m.ndim - 1].n));
                else t.unit_kind = VIAL_UNIT_UNKNOWN;
            }
            return t;
        }
        if (str_eq(name, STR_LIT("inv"))) { if (has_unit(&m)) set_unit(&t, md_unit_pow(m.unit, -1)); return t; }
        // eigh: the eigenvalues in m's unit, the eigenvectors without one
        T w = t;
        w.ndim -= 1;
        T v = t;
        no_unit(&v);
        T* elems = md_alloc(x->alloc, 2 * sizeof(T));
        elems[0] = w;
        elems[1] = v;
        T r = vial_type_scalar(VIAL_DT_TUPLE);
        r.elems = elems;
        r.num_elems = 2;
        join_stage(&r, &m);
        return r;
    }
    if (str_eq(name, STR_LIT("solve"))) {
        static const param_t p[] = {{"m", VX_KIND_WHOLE, true}, {"b", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2) || !all_numbers(x, fname, a, 2)) return err();
        const T m = type_of(x, a[0]), b = type_of(x, a[1]);
        if (!square(x, call, fname, &m)) return err();
        T t = b;
        if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT;
        if (has_unit(&m) || has_unit(&b)) {
            set_unit(&t, md_unit_div(b.unit_kind == VIAL_UNIT_KNOWN ? b.unit : md_unit_none(), m.unit_kind == VIAL_UNIT_KNOWN ? m.unit : md_unit_none()));
        }
        join_stage(&t, &m);
        return t;
    }
    return err();
}

// ### MISSING VALUES ###

static T call_missing(infer_t* x, vial_node_id call, str_t name, const char* fname) {
    if (str_eq(name, STR_LIT("missing_where"))) {
        static const param_t p[] = {{"c", VX_KIND_WHOLE, true}, {"x", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2)) return err();
        const T c = type_of(x, a[0]), v = type_of(x, a[1]);
        if (c.dtype != VIAL_DT_BOOL) { vx_diag_fix(x, "E0313", a[0], NULL, "`missing_where` takes a bool condition; this is %s", tstr(&c)); return err(); }
        T t;
        if (!broadcast_all(x, call, (vial_node_id[]){a[1], a[0]}, 2, &t)) return err();
        t.dtype = v.dtype;
        t.unit_kind = v.unit_kind; t.unit = v.unit;
        t.flags |= VIAL_TYPE_MISSING;
        return t;
    }
    if (str_eq(name, STR_LIT("valid"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}};
        vial_node_id a;
        if (!vx_bind(x, call, fname, p, 1, &a) || !all_typed(x, &a, 1)) return err();
        T t = type_of(x, a);
        t.dtype = VIAL_DT_BOOL;
        no_unit(&t);
        t.flags &= ~(VIAL_TYPE_MISSING | VIAL_TYPE_PERIODIC);
        return t;
    }
    if (str_eq(name, STR_LIT("filled"))) {
        static const param_t p[] = {{"x", VX_KIND_WHOLE, true}, {"v", VX_KIND_WHOLE, true}};
        vial_node_id a[2];
        if (!vx_bind(x, call, fname, p, 2, a) || !all_typed(x, a, 2)) return err();
        const T v = type_of(x, a[0]), f = type_of(x, a[1]);
        T t = v;
        T u = t;
        if (is_number(&v) && is_number(&f)) vx_unit_join(x, call, "fills", &v, &f, &u);
        t.unit_kind = u.unit_kind; t.unit = u.unit;
        t.flags &= ~VIAL_TYPE_MISSING;
        return t;
    }
    return err();
}

// ### SETS: common ###

static bool matches_by(const vial_sys_t* s, uint8_t level, uint32_t ua, uint32_t ub, bool label) {
    const struct md_system_t* sys = s->sys;
    (void)sys;
    if (level != VIAL_LEVEL_RESIDUE) return false;
    if (md_system_residue_seq_id(s->sys, ua) != md_system_residue_seq_id(s->sys, ub)) return false;
    return !label || str_eq(md_system_residue_name(s->sys, ua), md_system_residue_name(s->sys, ub));
}

static T call_common(infer_t* x, vial_node_id call) {
    static const param_t p[] = {{"a", VX_KIND_WHOLE, true}, {"b", VX_KIND_WHOLE, true}, {"by", VX_KIND_WHOLE, false}};
    vial_node_id a[3];
    if (!vx_bind(x, call, "common", p, 3, a) || !all_typed(x, a, 3)) return err();
    const T ta = type_of(x, a[0]), tb = type_of(x, a[1]);
    bool label = false;
    if (a[2]) {
        const vial_value_t v = value(x, a[2]);
        if (v.kind != VIAL_VAL_STR || !(str_eq(v.s, STR_LIT("id")) || str_eq(v.s, STR_LIT("label")))) {
            vx_diag_fix(x, "E0313", a[2], NULL, "`by` takes 'id' or 'label'");
            return err();
        }
        label = str_eq(v.s, STR_LIT("label"));
    }
    if (ta.dtype != VIAL_DT_SET || tb.dtype != VIAL_DT_SET || ta.ndim != 1 || tb.ndim != 1) {
        vx_diag_fix(x, "E0311", call, NULL, "`common` pairs arrays of units, such as residues; this is %s and %s", tstr(&ta), tstr(&tb));
        return err();
    }
    T r = ta;
    r.dims[0] = (vial_extent_t){.kind = ta.stage == VIAL_STAGE_STATIC && tb.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call};
    join_stage(&r, &tb);
    r.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    T* elems = md_alloc(x->alloc, 2 * sizeof(T));
    elems[0] = r;
    elems[1] = r;
    const vial_value_t va = value(x, a[0]), vb = value(x, a[1]);
    if (x->sys && va.kind == VIAL_VAL_SETS && vb.kind == VIAL_VAL_SETS && va.sets->unit && vb.sets->unit && va.sets->level == vb.sets->level &&
        va.sets->level == VIAL_LEVEL_RESIDUE && r.stage == VIAL_STAGE_STATIC) {
        // Each of a's units with the first of b's that matches, in a's order; one desc for both, as common paired them
        uint32_t* pa = NULL;
        uint32_t* pb = NULL;
        bool* used = md_alloc(x->alloc, vb.sets->count + 1);
        memset(used, 0, vb.sets->count + 1);
        for (uint32_t i = 0; i < va.sets->count; ++i) {
            for (uint32_t j = 0; j < vb.sets->count; ++j) {
                if (used[j] || !matches_by(x->sys, VIAL_LEVEL_RESIDUE, va.sets->unit[i], vb.sets->unit[j], label)) continue;
                used[j] = true;
                md_array_push(pa, i, x->alloc);
                md_array_push(pb, j, x->alloc);
                break;
            }
        }
        const uint32_t n = (uint32_t)md_array_size(pa);
        const vial_sets_t* sa = vial_sets_pick(va.sets, pa ? pa : (uint32_t[]){0}, n, x->alloc);
        const vial_sets_t* sb = vial_sets_pick(vb.sets, pb ? pb : (uint32_t[]){0}, n, x->alloc);
        elems[0].dims[0] = vial_extent_static(n);
        elems[1].dims[0] = vial_extent_static(n);
        elems[0].dims[0].desc = vx_desc(x, VIAL_LEVEL_RESIDUE, sa);
        elems[1].dims[0].desc = vx_desc(x, VIAL_LEVEL_RESIDUE, sb);
    }
    T t = vial_type_scalar(VIAL_DT_TUPLE);
    t.elems = elems;
    t.num_elems = 2;
    join_stage(&t, &r);
    return t;
}

// ### DISPATCH ###

bool vx_call_math(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    const char* fname = b->name;
    const int red = reduction_of(name);
    if (red >= 0) { *out = call_reduce(x, call, red, fname); return true; }
    static const char* math[] = {"abs", "floor", "ceil", "round", "clip", "sqrt", "exp", "log", "log10", "sin", "cos", "tan", "asin", "acos",
                                 "atan", "atan2", "dot", "cross", "norm", "normalize", "where", "unitless", "float", "int"};
    for (size_t i = 0; i < ARRAY_SIZE(math); ++i) if (str_eq(name, str_from_cstr(math[i]))) { *out = call_math(x, call, name, fname); return true; }
    static const char* build[] = {"shape", "zeros", "ones", "full", "zeros_like", "ones_like", "linspace", "arange", "eye", "reshape", "flatten", "diff"};
    for (size_t i = 0; i < ARRAY_SIZE(build); ++i) if (str_eq(name, str_from_cstr(build[i]))) { *out = call_build(x, call, name, fname); return true; }
    static const char* linalg[] = {"covariance", "eigh", "trace", "det", "inv", "solve"};
    for (size_t i = 0; i < ARRAY_SIZE(linalg); ++i) if (str_eq(name, str_from_cstr(linalg[i]))) { *out = call_linalg(x, call, name, fname); return true; }
    static const char* missing[] = {"missing_where", "valid", "filled"};
    for (size_t i = 0; i < ARRAY_SIZE(missing); ++i) if (str_eq(name, str_from_cstr(missing[i]))) { *out = call_missing(x, call, name, fname); return true; }
    if (str_eq(name, STR_LIT("concatenate"))) { *out = call_concatenate(x, call); return true; }
    if (str_eq(name, STR_LIT("common")))      { *out = call_common(x, call); return true; }
    return false;
}

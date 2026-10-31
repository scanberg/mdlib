// The checker, step 5: kinds, and the geometry built-ins.
//
// A parameter's kind says how much of an argument one call takes; what is left over on the left is mapped over, and
// the mapped dimensions of all arguments broadcast, as NumPy's generalized ufuncs do. The same split serves defs with
// annotated parameters (step 7). Positions are where a selection turns into coordinates: .xyz, a point or a points
// parameter. A set of several particles never quietly becomes one point (E0302).

#include "infer_internal.h"

#include <core/md_common.h>
#include <core/md_array.h>

#include <stdio.h>
#include <string.h>

typedef vial_type_t T;

static T err(void) { return vial_type_error(); }

static md_unit_t unit_named(const char* name) {
    md_unit_t u = md_unit_none();
    vial_unit_name(str_from_cstr(name), &u);
    return u;
}

static void join_stage(T* out, const T* a) { const T o = *out; vx_join_stage(out, &o, a); }

static bool is_float(const T* t) { return t->dtype == VIAL_DT_FLOAT || t->dtype == VIAL_DT_FLOAT64; }
static bool is_number(const T* t) { return t->dtype == VIAL_DT_INT || is_float(t) || t->dtype == VIAL_DT_BOOL; }

static bool last_is(const T* t, uint32_t n) {
    return t->ndim >= 1 && t->dims[t->ndim - 1].kind == VIAL_EXT_STATIC && t->dims[t->ndim - 1].n == n;
}

static const char* excerpt(char* buf, size_t cap, infer_t* x, vial_node_id node) {
    const str_t e = vial_excerpt(buf, cap, vx_text(x, node));
    (void)e;
    return buf;
}

// ### POSITIONS ###

// The smallest and largest set of an array of sets, when its value is known
static bool set_sizes(infer_t* x, vial_node_id node, uint32_t* lo, uint32_t* hi) {
    const vial_value_t v = vx_value_of(x, node);
    if (v.kind != VIAL_VAL_SETS || v.sets->single || v.sets->count == 0) return false;
    *lo = UINT32_MAX; *hi = 0;
    for (uint32_t i = 0; i < v.sets->count; ++i) {
        const uint32_t n = vial_sets_size(v.sets, i);
        if (n < *lo) *lo = n;
        if (n > *hi) *hi = n;
    }
    return true;
}

// The extent of a set's particles: static when its value is known, ? per frame when it changes
static vial_extent_t particle_extent(infer_t* x, const T* s, vial_node_id node) {
    const vial_value_t v = vx_value_of(x, node);
    if (s->stage != VIAL_STAGE_STATIC) return (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = s->ident ? s->ident : (uint32_t)node};
    if (v.kind == VIAL_VAL_SETS && v.sets->single) {
        vial_extent_t e = vial_extent_static(vial_sets_size(v.sets, 0));
        e.desc = vx_desc(x, VIAL_LEVEL_PARTICLE, v.sets);
        return e;
    }
    return (vial_extent_t){.kind = VIAL_EXT_UNKNOWN, .id = (uint32_t)node};
}

// A length in Å, static: a result's stage comes from its arguments
static T length(void) {
    T t = vial_type_scalar(VIAL_DT_FLOAT);
    t.unit = unit_named("angstrom");
    return t;
}

// Positions of particles: Å, frame
static T coords(void) {
    T t = length();
    t.stage = VIAL_STAGE_FRAME;
    return t;
}

// A fix that writes the whole of a node's text inside pre and post; NULL when it would be too long to offer
static const char* fix_around(char* buf, size_t cap, infer_t* x, vial_node_id node, const char* pre, const char* post) {
    const str_t t = vx_text(x, node);
    const int n = snprintf(buf, cap, "%s%.*s%s", pre, (int)t.len, t.ptr, post);
    return n > 0 && (size_t)n < cap ? buf : NULL;
}

// .xyz of a set: float[n][3] Å, frame; of an array of one-particle sets, one position each; of other arrays, E0318
static bool positions_of(infer_t* x, vial_node_id node, const T* s, vial_node_id at, T* out) {
    T t = coords();
    join_stage(&t, s);
    if (s->ndim == 0) {
        vial_type_push_dim(&t, particle_extent(x, s, node));
        vial_type_push_dim(&t, vial_extent_static(3));
        *out = t;
        return true;
    }
    if (s->flags & VIAL_TYPE_ONE_PARTICLE) {
        for (uint8_t k = 0; k < s->ndim; ++k) vial_type_push_dim(&t, s->dims[k]);
        if (!vial_type_push_dim(&t, vial_extent_static(3))) {
            vx_diag_fix(x, "E0308", at, NULL, "the result would have 7 dimensions; the limit is 6");
            return false;
        }
        *out = t;
        return true;
    }
    char ex[VIAL_EXCERPT_MAX + 4], buf[256];
    excerpt(ex, sizeof(ex), x, node);
    const char* fix = fix_around(buf, sizeof(buf), x, node, "particles(", ").xyz");
    uint32_t lo, hi;
    if (set_sizes(x, node, &lo, &hi) && lo != hi) {
        vx_diag_fix(x, "E0318", at, fix, "`%s` holds %u sets of %u to %u particles; their positions do not form one array", ex,
            s->dims[0].n, lo, hi);
    } else {
        vx_diag_fix(x, "E0318", at, fix, "`%s` holds sets whose sizes may differ; their positions do not form one array", ex);
    }
    return false;
}

// ### KINDS ###

enum {
    VK_WHOLE = VX_KIND_WHOLE,   // unannotated: taken whole
    VK_SET = VX_KIND_SET,
    VK_GROUP = VX_KIND_GROUP,
    VK_POINT = VX_KIND_POINT,
    VK_POINTS = VX_KIND_POINTS,
    VK_FLOAT = VX_KIND_FLOAT,
    VK_INT = VX_KIND_INT,
    VK_BOOL = VX_KIND_BOOL,
};

// What one call takes of an argument, and what is mapped over
typedef struct split_t {
    T       core;
    T       loop;       // only its dims count
} split_t;

static void set_loop(split_t* s, const T* a, uint8_t ndim) {
    s->loop = vial_type_scalar(VIAL_DT_INT);
    for (uint8_t k = 0; k < ndim; ++k) vial_type_push_dim(&s->loop, a->dims[k]);
}

// The core of an array of sets under points or group: each set's positions, of a size that may change from set to set
static T ragged_positions(infer_t* x, vial_node_id node, const T* a) {
    T t = coords();
    join_stage(&t, a);
    uint32_t lo, hi;
    if (a->flags & VIAL_TYPE_ONE_PARTICLE) vial_type_push_dim(&t, vial_extent_static(1));
    else if (a->stage == VIAL_STAGE_STATIC && set_sizes(x, node, &lo, &hi)) {
        vial_type_push_dim(&t, lo == hi ? vial_extent_static(lo) : (vial_extent_t){.kind = VIAL_EXT_RANGED, .n = lo, .m = hi});
    } else {
        vial_type_push_dim(&t, (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = (uint32_t)node});
    }
    // One of the sets' particles, a different set per call: described as particles, which are not known here
    t.dims[0].desc = vx_desc(x, VIAL_LEVEL_PARTICLE, NULL);
    vial_type_push_dim(&t, vial_extent_static(3));
    return t;
}

static const char* kind_name(uint8_t k) {
    static const char* names[] = {"anything", "set", "group", "point", "points", "float", "int", "bool"};
    return k < ARRAY_SIZE(names) ? names[k] : "?";
}

static bool split(infer_t* x, vial_node_id node, const T* a, uint8_t kind, const char* fname, const char* pname, split_t* out) {
    memset(out, 0, sizeof(*out));
    out->loop = vial_type_scalar(VIAL_DT_INT);
    char ex[VIAL_EXCERPT_MAX + 4], buf[256];
    const char* fix = NULL;
    switch (kind) {
    case VK_WHOLE:
        out->core = *a;
        return true;
    case VK_POINT:
        if (is_float(a) && last_is(a, 3)) {
            out->core = *a;
            out->core.ndim = 1;
            out->core.dims[0] = a->dims[a->ndim - 1];
            if (out->core.unit_kind == VIAL_UNIT_OPEN) { out->core.unit_kind = VIAL_UNIT_KNOWN; out->core.unit = unit_named("angstrom"); }
            set_loop(out, a, a->ndim - 1);
            return true;
        }
        if (a->dtype == VIAL_DT_SET && (a->flags & VIAL_TYPE_ONE_PARTICLE)) {
            out->core = coords();
            join_stage(&out->core, a);
            vial_type_push_dim(&out->core, vial_extent_static(3));
            set_loop(out, a, a->ndim);
            return true;
        }
        if (a->dtype == VIAL_DT_SET && a->ndim == 0 && (a->flags & VIAL_TYPE_EMPTY)) {
            // A set known to be empty (W0601 where it was selected) has no position in any frame
            out->core = coords();
            out->core.flags |= VIAL_TYPE_MISSING;
            vial_type_push_dim(&out->core, vial_extent_static(3));
            vx_diag_fix(x, "W0405", node, NULL, "`%s` is empty, so what `%s` gives is missing in every frame", excerpt(ex, sizeof(ex), x, node), fname);
            return true;
        }
        excerpt(ex, sizeof(ex), x, node);
        if (a->dtype == VIAL_DT_SET) {
            fix = fix_around(buf, sizeof(buf), x, node, "", ".com");
            vx_diag_fix(x, "E0302", node, fix, "`%s` needs a point for %s; `%s` %s, so write .com for its centre or between='closest'", fname,
                pname, ex, a->ndim ? "holds sets that may hold several particles" : "may hold several particles");
            return false;
        }
        if (a->dtype == VIAL_DT_INT && a->ndim == 0) {
            fix = fix_around(buf, sizeof(buf), x, node, "particles(serial=", ")");
            vx_diag_fix(x, "E0302", node, fix, "`%s` needs a point for %s; an int is not one", fname, pname);
            return false;
        }
        vx_diag_fix(x, "E0313", node, NULL, "`%s` takes a point for %s, a last dimension of 3 or a set of one particle; this is %s", fname,
            pname, vx_tstr(a));
        return false;
    case VK_POINTS:
    case VK_GROUP:
        if (is_float(a) && last_is(a, 3)) {
            if (a->ndim == 1 || kind == VK_GROUP) {
                // One point: a group of one location. Under group, every point is a group of its own.
                out->core = *a;
                out->core.ndim = 2;
                out->core.dims[0] = vial_extent_static(1);
                out->core.dims[1] = a->dims[a->ndim - 1];
                set_loop(out, a, a->ndim - 1);
            } else {
                out->core = *a;
                out->core.ndim = 2;
                out->core.dims[0] = a->dims[a->ndim - 2];
                out->core.dims[1] = a->dims[a->ndim - 1];
                set_loop(out, a, a->ndim - 2);
            }
            if (out->core.unit_kind == VIAL_UNIT_OPEN) { out->core.unit_kind = VIAL_UNIT_KNOWN; out->core.unit = unit_named("angstrom"); }
            return true;
        }
        if (a->dtype == VIAL_DT_SET) {
            if (a->ndim == 0) {
                T t = coords();
                join_stage(&t, a);
                vial_type_push_dim(&t, particle_extent(x, a, node));
                vial_type_push_dim(&t, vial_extent_static(3));
                out->core = t;
                return true;
            }
            out->core = ragged_positions(x, node, a);
            set_loop(out, a, a->ndim);
            return true;
        }
        if (a->dtype == VIAL_DT_INT && a->ndim == 0) {
            fix = fix_around(buf, sizeof(buf), x, node, "particles(serial=", ")");
            vx_diag_fix(x, "E0302", node, fix, "`%s` needs %s for %s; an int is not one", fname, kind_name(kind), pname);
            return false;
        }
        vx_diag_fix(x, "E0313", node, NULL, "`%s` takes %s for %s, positions or a set; this is %s", fname, kind_name(kind), pname, vx_tstr(a));
        return false;
    case VK_SET:
        if (a->dtype == VIAL_DT_SET) {
            out->core = *a;
            out->core.ndim = 0;
            set_loop(out, a, a->ndim);
            return true;
        }
        vx_diag_fix(x, "E0313", node, NULL, "`%s` takes a set for %s; this is %s", fname, pname, vx_tstr(a));
        return false;
    case VK_FLOAT:
    case VK_INT:
    case VK_BOOL: {
        const bool fits = kind == VK_BOOL ? a->dtype == VIAL_DT_BOOL : kind == VK_INT ? (a->dtype == VIAL_DT_INT || a->dtype == VIAL_DT_BOOL) : is_number(a);
        if (!fits) {
            vx_diag_fix(x, "E0313", node, NULL, "`%s` takes %s for %s; this is %s", fname, kind_name(kind), pname, vx_tstr(a));
            return false;
        }
        out->core = *a;
        out->core.ndim = 0;
        set_loop(out, a, a->ndim);
        return true;
    }
    default:
        return false;
    }
}

bool vx_split(infer_t* x, vial_node_id node, const T* a, uint8_t kind, const char* fname, const char* pname, T* core, T* loop) {
    split_t s;
    if (!split(x, node, a, kind, fname, pname, &s)) return false;
    *core = s.core;
    *loop = s.loop;
    return true;
}

// ### ARGUMENTS ###

typedef vx_param_t param_t;

// Binds a call's arguments to parameters by position, then by keyword; reports E0202 and E0204. Every argument is typed.
static bool bind(infer_t* x, vial_node_id call, const char* fname, const param_t* params, uint32_t np, vial_node_id* args) {
    vx_record_sig(x, call, params, np);
    bool ok = true;
    uint32_t pos = 0;
    memset(args, 0, np * sizeof(vial_node_id));
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, a, 0), vn = vx_child(x, a, 1);
            if (vn) vx_infer(x, vn);
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            uint32_t k = 0;
            while (k < np && !str_eq(kw, str_from_cstr(params[k].name))) ++k;
            if (k == np) {
                char names[160] = "";
                size_t n = 0;
                for (uint32_t q = 0; q < np && n < sizeof(names); ++q) n += snprintf(names + n, sizeof(names) - n, "%s%s", q ? ", " : "", params[q].name);
                vx_diag_fix(x, "E0202", a, NULL, "`%s` takes no keyword `%.*s`; it takes %s", fname, (int)kw.len, kw.ptr, names);
                ok = false;
                continue;
            }
            if (args[k]) {
                vx_diag_fix(x, "E0204", a, NULL, "`%s` is given %s twice", fname, params[k].name);
                ok = false;
                continue;
            }
            args[k] = vn;
            continue;
        }
        vx_infer(x, a);
        if (pos >= np) {
            vx_diag_fix(x, "E0204", a, NULL, "`%s` takes at most %u arguments", fname, np);
            ok = false;
            continue;
        }
        args[pos++] = a;
    }
    for (uint32_t k = 0; ok && k < np; ++k) {
        if (params[k].required && !args[k]) {
            vx_diag_fix(x, "E0204", call, NULL, "`%s` needs %s", fname, params[k].name);
            ok = false;
        }
    }
    if (ok && x->ir) vx_record_bind(x, call, params, np, args);
    return ok;
}

bool vx_bind(infer_t* x, vial_node_id call, const char* fname, const vx_param_t* params, uint32_t np, vial_node_id* args) {
    return bind(x, call, fname, params, np, args);
}

// The types of bound arguments, split by kind, and their mapped dimensions broadcast. false when an argument failed
// or the mapped dimensions do not broadcast (E0301).
static bool map_args(infer_t* x, vial_node_id call, const char* fname, const param_t* params, uint32_t np, const vial_node_id* args,
                     split_t* splits, T* loop) {
    *loop = vial_type_scalar(VIAL_DT_INT);
    bool ok = true;
    for (uint32_t k = 0; k < np; ++k) {
        if (!args[k]) { memset(&splits[k], 0, sizeof(splits[k])); continue; }
        const T a = x->out->types[(args[k]) - x->out->base];
        if (vial_type_is_error(&a)) { ok = false; continue; }
        if (!split(x, args[k], &a, params[k].kind, fname, params[k].name, &splits[k])) { ok = false; continue; }
        T m;
        vx_check_pairing(x, call, loop, &splits[k].loop);
        if (!vx_broadcast(loop, &splits[k].loop, &m)) {
            vx_diag_fix(x, "E0301", call, NULL, "the dimensions `%s` maps over do not broadcast: `%s` and `%s`", fname, vx_tstr(loop),
                vx_tstr(&splits[k].loop));
            return false;
        }
        loop->ndim = m.ndim;
        memcpy(loop->dims, m.dims, sizeof(m.dims));
        join_stage(loop, &splits[k].core);
        join_stage(loop, &a);
    }
    return ok;
}

// The mapped dimensions, then the core's
static bool with_loop(infer_t* x, vial_node_id call, const T* loop, T core, T* out) {
    T t = core;
    t.ndim = 0;
    for (uint8_t k = 0; k < loop->ndim; ++k) vial_type_push_dim(&t, loop->dims[k]);
    for (uint8_t k = 0; k < core.ndim; ++k) {
        if (!vial_type_push_dim(&t, core.dims[k])) {
            vx_diag_fix(x, "E0308", call, NULL, "the result would have 7 dimensions; the limit is 6");
            return false;
        }
    }
    join_stage(&t, loop);
    *out = t;
    return true;
}

static vial_value_t value(infer_t* x, vial_node_id n) { return n ? vx_value_of(x, n) : (vial_value_t){0}; }

// ### THE BUILT-INS ###

static T call_com(infer_t* x, vial_node_id call, const char* fname, vial_node_id s_node, vial_node_id weights_node) {
    static const param_t params[] = {{"s", VK_POINTS, true}, {"weights", VK_WHOLE, false}};
    vial_node_id args[2] = {s_node, weights_node};
    if (!s_node && !bind(x, call, fname, params, 2, args)) return err();
    if (args[1]) {
        const vial_value_t w = value(x, args[1]);
        const T wt = x->out->types[(args[1]) - x->out->base];
        const bool ok = (w.kind == VIAL_VAL_STR && str_eq(w.s, STR_LIT("mass"))) || w.kind == VIAL_VAL_NULL || vial_type_is_error(&wt);
        if (!ok) {
            vx_diag_fix(x, "E0313", args[1], NULL, "`weights` takes 'mass' or None");
            return err();
        }
    }
    split_t sp[2];
    T loop;
    if (!map_args(x, call, fname, params, 2, args, sp, &loop)) return err();
    T core = length();
    vial_type_push_dim(&core, vial_extent_static(3));
    T t;
    return with_loop(x, call, &loop, core, &t) ? t : err();
}

// mass, charge, element, name, bfactor: one value per particle of a set
static T call_property(infer_t* x, vial_node_id call, str_t name) {
    static const param_t params[] = {{"s", VK_WHOLE, true}};
    vial_node_id args[1];
    char fname[16];
    snprintf(fname, sizeof(fname), "%.*s", (int)name.len, name.ptr);
    if (!bind(x, call, fname, params, 1, args)) return err();
    const T s = x->out->types[(args[0]) - x->out->base];
    if (vial_type_is_error(&s)) return err();
    if (s.dtype != VIAL_DT_SET) {
        vx_diag_fix(x, "E0313", args[0], NULL, "`%s` takes a set; this is %s", fname, vx_tstr(&s));
        return err();
    }
    if (s.ndim) {
        char ex[VIAL_EXCERPT_MAX + 4], buf[256], pre[32];
        excerpt(ex, sizeof(ex), x, args[0]);
        snprintf(pre, sizeof(pre), "%s(particles(", fname);
        const char* fix = fix_around(buf, sizeof(buf), x, args[0], pre, "))");
        vx_diag_fix(x, "E0318", args[0], fix, "`%s` holds sets whose sizes may differ; their particles' %ss do not form one array", ex, fname);
        return err();
    }
    const bool str = str_eq(name, STR_LIT("element")) || str_eq(name, STR_LIT("name"));
    T t = vial_type_scalar(str ? VIAL_DT_STR : VIAL_DT_FLOAT);
    if (str_eq(name, STR_LIT("mass")))    t.unit = unit_named("u");
    if (str_eq(name, STR_LIT("charge")))  t.unit = unit_named("e");
    if (str_eq(name, STR_LIT("bfactor"))) t.unit = md_unit_mul(unit_named("angstrom"), unit_named("angstrom"));
    t.stage = s.stage;
    vial_type_push_dim(&t, particle_extent(x, &s, args[0]));
    return t;
}

static T call_displacement(infer_t* x, vial_node_id call) {
    static const param_t params[] = {{"a", VK_POINT, true}, {"b", VK_POINT, true}};
    vial_node_id args[2];
    split_t sp[2];
    T loop;
    if (!bind(x, call, "displacement", params, 2, args) || !map_args(x, call, "displacement", params, 2, args, sp, &loop)) return err();
    T core = length();
    vial_type_push_dim(&core, vial_extent_static(3));
    T t;
    return with_loop(x, call, &loop, core, &t) ? t : err();
}

// angle(a, b, c) and dihedral(a, b, c, d) on points, or on one argument whose last two dimensions are [n][3], or an
// array of one-particle sets whose last dimension is n
static T call_angle(infer_t* x, vial_node_id call, bool dihedral) {
    const char* fname = dihedral ? "dihedral" : "angle";
    const uint32_t n = dihedral ? 4 : 3;
    static const param_t params4[] = {{"a", VK_POINT, true}, {"b", VK_POINT, true}, {"c", VK_POINT, true}, {"d", VK_POINT, true}};
    vx_record_sig(x, call, params4, n);     // what completion offers while the call is being written
    uint32_t num_pos = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) num_pos += vx_child(x, call, i) && vx_kind(x, vx_child(x, call, i)) != VIAL_N_KWARG;
    T loop;
    if (num_pos == 1) {
        static const param_t one[] = {{"points", VK_WHOLE, true}};
        vial_node_id arg;
        if (!bind(x, call, fname, one, 1, &arg)) return err();
        const T a = x->out->types[(arg) - x->out->base];
        if (vial_type_is_error(&a)) return err();
        loop = vial_type_scalar(VIAL_DT_INT);
        const bool rows = is_float(&a) && a.ndim >= 2 && last_is(&a, 3) && a.dims[a.ndim - 2].kind == VIAL_EXT_STATIC && a.dims[a.ndim - 2].n == n;
        const bool sets = a.dtype == VIAL_DT_SET && (a.flags & VIAL_TYPE_ONE_PARTICLE) && last_is(&a, n);
        if (!rows && !sets) {
            vx_diag_fix(x, "E0313", arg, NULL, "`%s` takes %u points, or one array of them, float[%u][3] or set[%u] of single particles; this is %s",
                fname, n, n, n, vx_tstr(&a));
            return err();
        }
        for (uint8_t k = 0; k + (rows ? 2 : 1) < a.ndim + 0u; ++k) vial_type_push_dim(&loop, a.dims[k]);
        join_stage(&loop, &a);
        if (sets) loop.stage = loop.stage < VIAL_STAGE_FRAME ? VIAL_STAGE_FRAME : loop.stage;
    } else {
        if (num_pos != n && num_pos != 0) {
            for (uint32_t i = 1; i < vx_num_children(x, call); ++i) if (vx_child(x, call, i)) vx_infer(x, vx_child(x, call, i));
            vx_diag_fix(x, "E0204", call, NULL, "`%s` takes %u points; it got %u", fname, n, num_pos);
            return err();
        }
        vial_node_id args[4];
        split_t sp[4];
        if (!bind(x, call, fname, params4, n, args) || !map_args(x, call, fname, params4, n, args, sp, &loop)) return err();
    }
    T core = vial_type_scalar(VIAL_DT_FLOAT);
    core.unit = unit_named("deg");
    if (dihedral) {
        core.flags |= VIAL_TYPE_PERIODIC;
        core.period_beg = -180.0f; core.period_end = 180.0f; core.seam = 180.0f;
    }
    T t;
    return with_loop(x, call, &loop, core, &t) ? t : err();
}

// The descriptions of a points core's point dimension and of the mapped ones, for W0608
static bool shares_units(const split_t* p, const split_t* q) {
    uint32_t dp[VIAL_MAX_DIMS + 1], nq = 0;
    uint32_t n = 0;
    for (uint8_t k = 0; k < p->loop.ndim; ++k) if (p->loop.dims[k].desc) dp[n++] = p->loop.dims[k].desc;
    if (p->core.ndim >= 2 && p->core.dims[0].desc) dp[n++] = p->core.dims[0].desc;
    const uint32_t want = q->core.ndim >= 2 ? q->core.dims[0].desc : 0;
    (void)nq;
    for (uint32_t i = 0; i < n && want; ++i) if (dp[i] == want) return true;
    return false;
}

static T call_pair(infer_t* x, vial_node_id call, bool closest) {
    const char* fname = closest ? "closest_pair" : "farthest_pair";
    static const param_t params[] = {{"p", VK_POINTS, true}, {"q", VK_POINTS, true}, {"distinct", VK_WHOLE, false}};
    vial_node_id args[3];
    split_t sp[3];
    T loop;
    if (!bind(x, call, fname, params, 3, args) || !map_args(x, call, fname, params, 3, args, sp, &loop)) return err();
    bool distinct = false;
    if (args[2]) {
        const vial_value_t v = value(x, args[2]);
        if (v.kind != VIAL_VAL_BOOL) {
            vx_diag_fix(x, "E0313", args[2], NULL, "`distinct` is True or False, known when compiling");
            return err();
        }
        distinct = v.b;
    }
    if (closest && !distinct && shares_units(&sp[0], &sp[1])) {
        const str_t t = vx_text(x, call);
        char fix[160] = "";
        if (t.len > 1 && t.len < 120 && t.ptr[t.len - 1] == ')') snprintf(fix, sizeof(fix), "%.*s, distinct=True)", (int)(t.len - 1), t.ptr);
        vx_diag_fix(x, "W0608", call, fix[0] ? fix : NULL, "`%s` pairs every point with itself, so the shortest distance is zero; distinct=True leaves those pairs out", fname);
    }
    T d = length();
    T i = vial_type_scalar(VIAL_DT_INT);
    T* elems = md_alloc(x->alloc, 3 * sizeof(T));
    if (!with_loop(x, call, &loop, d, &elems[0]) || !with_loop(x, call, &loop, i, &elems[1])) return err();
    elems[2] = elems[1];
    T t = vial_type_scalar(VIAL_DT_TUPLE);
    t.elems = elems;
    t.num_elems = 3;
    join_stage(&t, &elems[0]);
    return t;
}

static T call_fit(infer_t* x, vial_node_id call) {
    static const param_t params[] = {{"p", VK_POINTS, true}, {"ref", VK_POINTS, true}, {"on", VK_WHOLE, false}, {"weights", VK_WHOLE, false}};
    vial_node_id args[4];
    split_t sp[4];
    T loop;
    if (!bind(x, call, "fit", params, 4, args) || !map_args(x, call, "fit", params, 4, args, sp, &loop)) return err();
    if (args[2]) {
        const T on = x->out->types[(args[2]) - x->out->base];
        if (on.dtype != VIAL_DT_SET && on.dtype != VIAL_DT_NONE) {
            vx_diag_fix(x, "E0313", args[2], NULL, "`on` takes a set or None; this is %s", vx_tstr(&on));
            return err();
        }
    }
    T core = sp[0].core;
    core.stage = core.stage < VIAL_STAGE_FRAME ? VIAL_STAGE_FRAME : core.stage;
    T t;
    if (!with_loop(x, call, &loop, core, &t)) return err();
    join_stage(&t, &sp[1].core);
    return t;
}

bool vx_call_geom(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    if (str_eq(name, STR_LIT("com")))          { *out = call_com(x, call, "com", 0, 0); return true; }
    if (str_eq(name, STR_LIT("positions"))) {
        static const param_t params[] = {{"s", VK_WHOLE, true}};
        vial_node_id arg;
        *out = err();
        if (!bind(x, call, "positions", params, 1, &arg)) return true;
        const T a = x->out->types[(arg) - x->out->base];
        if (vial_type_is_error(&a)) return true;
        if (a.dtype == VIAL_DT_SET) positions_of(x, arg, &a, call, out);
        else if (is_float(&a) && last_is(&a, 3)) *out = a;
        else vx_diag_fix(x, "E0313", arg, NULL, "`positions` takes a set or coordinates; this is %s", vx_tstr(&a));
        return true;
    }
    if (str_eq(name, STR_LIT("mass")) || str_eq(name, STR_LIT("charge")) || str_eq(name, STR_LIT("element")) ||
        str_eq(name, STR_LIT("name")) || str_eq(name, STR_LIT("bfactor"))) { *out = call_property(x, call, name); return true; }
    if (str_eq(name, STR_LIT("displacement")))  { *out = call_displacement(x, call); return true; }
    if (str_eq(name, STR_LIT("angle")))         { *out = call_angle(x, call, false); return true; }
    if (str_eq(name, STR_LIT("dihedral")))      { *out = call_angle(x, call, true); return true; }
    if (str_eq(name, STR_LIT("closest_pair")))  { *out = call_pair(x, call, true); return true; }
    if (str_eq(name, STR_LIT("farthest_pair"))) { *out = call_pair(x, call, false); return true; }
    if (str_eq(name, STR_LIT("fit")))           { *out = call_fit(x, call); return true; }
    return false;
}

bool vx_attr_geom(infer_t* x, vial_node_id attr, vial_node_id on, const T* a, str_t name, T* out) {
    if (str_eq(name, STR_LIT("com"))) {
        // s.com is com(s)
        vx_infer(x, on);
        *out = call_com(x, attr, "com", on, 0);
        return true;
    }
    if (a->dtype != VIAL_DT_SET) return false;
    const bool comp = str_eq(name, STR_LIT("x")) || str_eq(name, STR_LIT("y")) || str_eq(name, STR_LIT("z"));
    if (!comp && !str_eq(name, STR_LIT("xyz"))) return false;     // an unknown name is E0201, from name resolution
    T t;
    if (!positions_of(x, on, a, attr, &t)) { *out = err(); return true; }
    if (comp) t.ndim -= 1;
    *out = t;
    return true;
}

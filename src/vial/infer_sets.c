// The checker, step 4: the selectors, the predefined sets, set algebra and indexing, and len, count and labels.
//
// With a system, every static selector is evaluated when compiling: its sets are known, so their sizes, counts and
// labels are too, and W0601 (matches nothing), W0602 (a resid in several chains) and E0207 (no chains) are given.
// Without one, the same calls are typed with extents not known. A selector with a frame-stage argument, such as a
// spatial set or a scope that changes per frame, is frame stage, its length known only per frame.

#include "infer_internal.h"

#include <md_system.h>

#include <core/md_common.h>
#include <core/md_array.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef vial_type_t T;

static T err(void) { return vial_type_error(); }

static void join_stage(T* out, const T* a) { const T o = *out; vx_join_stage(out, &o, a); }

static vial_value_t sets_value(const vial_sets_t* s) { return (vial_value_t){.kind = VIAL_VAL_SETS, .sets = s}; }

T vx_sets_type(infer_t* x, const vial_sets_t* s, vial_node_id at) {
    (void)at;
    T t = vial_type_scalar(VIAL_DT_SET);
    if (s->single) {
        const uint32_t n = vial_sets_size(s, 0);
        if (n == 1) t.flags |= VIAL_TYPE_ONE_PARTICLE;
        if (n == 0) t.flags |= VIAL_TYPE_EMPTY;
        return t;
    }
    vial_extent_t e = vial_extent_static(s->count);
    e.desc = vx_desc(x, s->level, s);
    vial_type_push_dim(&t, e);
    // An array whose every set holds one particle is an array of points, such as the waters of an X-ray structure
    bool one = s->count > 0;
    for (uint32_t i = 0; one && i < s->count; ++i) one = vial_sets_size(s, i) == 1;
    if (one) t.flags |= VIAL_TYPE_ONE_PARTICLE;
    return t;
}

// ### SELECTORS ###

static uint8_t selector_level(str_t name, uint8_t* particle_kind) {
    *particle_kind = 0;
    if (str_eq(name, STR_LIT("particles"))) return VIAL_LEVEL_PARTICLE;
    if (str_eq(name, STR_LIT("atoms"))) { *particle_kind = 1; return VIAL_LEVEL_PARTICLE; }
    if (str_eq(name, STR_LIT("beads"))) { *particle_kind = 2; return VIAL_LEVEL_PARTICLE; }
    if (str_eq(name, STR_LIT("residues"))) return VIAL_LEVEL_RESIDUE;
    if (str_eq(name, STR_LIT("parts"))) return VIAL_LEVEL_PART;
    if (str_eq(name, STR_LIT("molecules"))) return VIAL_LEVEL_MOLECULE;
    if (str_eq(name, STR_LIT("chains"))) return VIAL_LEVEL_CHAIN;
    return VIAL_LEVEL_NONE;
}

// The type of a selector whose sets are not known: one set, or an array of a length not known
static T selector_unknown(uint8_t level, uint8_t stage, vial_node_id at) {
    T t = vial_type_scalar(VIAL_DT_SET);
    t.stage = stage;
    t.ident = (uint32_t)at;
    if (level == VIAL_LEVEL_PARTICLE) return t;
    vial_extent_t e = {.kind = stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)at};
    vial_type_push_dim(&t, e);
    return t;
}

static bool no_chains(const infer_t* x) { return x->sys && x->sys->num_chains == 0; }

static void report_no_chains(infer_t* x, vial_node_id at) {
    vx_diag_fix(x, "E0207", at, "parts(kind='polymer')",
        "the data gives no chain ids, and chains are never inferred; a strand is found as a part");
}

// Whether a keyword takes strings or ints
static bool key_takes_str(uint8_t level, uint8_t key) {
    switch (key) {
    case VIAL_KEY_NAME: case VIAL_KEY_ELEMENT: case VIAL_KEY_KIND: case VIAL_KEY_RESNAME: case VIAL_KEY_CHAIN: return true;
    case VIAL_KEY_ID: return level == VIAL_LEVEL_CHAIN;
    default: return false;
    }
}

static const char* key_word(uint8_t key) {
    static const char* words[] = {"", "name", "serial", "element", "kind", "resid", "resname", "chain", "id", "scope"};
    return key < ARRAY_SIZE(words) ? words[key] : "";
}

static bool int_of(const vial_value_t* v, int64_t* out) {
    if (v->kind == VIAL_VAL_INT) { *out = v->i; return true; }
    if (v->kind == VIAL_VAL_FLOAT && v->f >= -9.0e15 && v->f <= 9.0e15 && v->f == floor(v->f)) { *out = (int64_t)v->f; return true; }
    return false;
}

// One value of a keyword, from a node already typed: a string, an int or a range. 1 when known, 0 when not known
// when compiling, -1 when it is a value the keyword does not take.
static int match_of(infer_t* x, vial_node_id node, bool wants_str, vial_match_t* m) {
    const T t = x->out->types[(node) - x->out->base];
    const vial_value_t v = vx_value_of(x, node);
    memset(m, 0, sizeof(*m));
    if (wants_str) {
        if (t.dtype != VIAL_DT_STR || t.ndim) return -1;
        if (v.kind != VIAL_VAL_STR) return 0;
        m->kind = 0; m->s = v.s;
        return 1;
    }
    if (t.dtype == VIAL_DT_RANGE && vx_kind(x, node) == VIAL_N_RANGE) {
        const vial_value_t lo = vx_value_of(x, vx_child(x, node, 0)), hi = vx_value_of(x, vx_child(x, node, 1));
        if (!int_of(&lo, &m->lo) || !int_of(&hi, &m->hi)) return 0;
        m->kind = 2;
        return 1;
    }
    if (t.dtype == VIAL_DT_RANGE) return 0;
    if (t.dtype != VIAL_DT_INT || t.ndim) return -1;
    if (!int_of(&v, &m->lo)) return 0;
    m->kind = 1;
    return 1;
}

// A keyword's value: one value, or a list of them. Typed already. 1 known, 0 not known, -1 refused (reported).
static int filter_of(infer_t* x, uint8_t level, uint8_t key, vial_node_id node, vial_filter_t* f) {
    const bool wants_str = key_takes_str(level, key);
    const T t = x->out->types[(node) - x->out->base];
    f->key = key;
    f->any = NULL;
    int known = 1;
    if (vx_kind(x, node) == VIAL_N_LIST && t.ndim == 1) {
        for (uint32_t i = 0; i < vx_num_children(x, node); ++i) {
            vial_match_t m;
            const int r = match_of(x, vx_child(x, node, i), wants_str, &m);
            if (r < 0) { known = -1; break; }
            if (r == 0) known = 0;
            else md_array_push(f->any, m, x->alloc);
        }
    } else if (t.ndim == 1 && t.dtype == (wants_str ? VIAL_DT_STR : VIAL_DT_INT)) {
        known = 0;      // a list from elsewhere: typed, not evaluated here
    } else {
        vial_match_t m;
        known = match_of(x, node, wants_str, &m);
        if (known == 1) md_array_push(f->any, m, x->alloc);
    }
    if (known < 0) {
        vx_diag_fix(x, "E0206", node, NULL, "`%s` takes %s; this is %s", key_word(key),
            wants_str ? "a string or a list of strings" : "an int, a range a..b, or a list of ints", vx_tstr(&t));
        return -1;
    }
    if (key == VIAL_KEY_KIND && known == 1) {
        for (size_t i = 0; i < md_array_size(f->any); ++i) {
            if (vial_select_kind_known(level, f->any[i].s)) continue;
            // A kind that starts the same, such as amino for amino_acid, is offered as the fix
            const char* kinds = vial_select_kinds(level);
            char fix[64] = "";
            const str_t w = f->any[i].s;
            for (const char* p = kinds; *p && w.len;) {
                const char* e = p;
                while (*e && *e != ',' && *e != ' ') ++e;
                if ((size_t)(e - p) > w.len && strncmp(p, w.ptr, w.len) == 0) { snprintf(fix, sizeof(fix), "'%.*s'", (int)(e - p), p); break; }
                p = e;
                while (*p == ',' || *p == ' ' || strncmp(p, "and ", 4) == 0) p += *p == 'a' ? 4 : 1;
            }
            char ex[VIAL_EXCERPT_MAX + 4];
            const str_t e = vial_excerpt(ex, sizeof(ex), w);
            vx_diag_fix(x, "E0206", node, fix[0] ? fix : NULL, "'%.*s' is not a %s kind; the kinds are %s", (int)e.len, e.ptr,
                level == VIAL_LEVEL_PARTICLE ? "particle" : level == VIAL_LEVEL_RESIDUE ? "residue" : "part", kinds);
            return -1;
        }
    }
    return known;
}

// The keyword a selector may take instead of one it does not, as E0202's fix
static uint8_t key_instead(uint8_t level, uint8_t key) {
    if (level == VIAL_LEVEL_CHAIN && key == VIAL_KEY_NAME) return VIAL_KEY_ID;
    if (level == VIAL_LEVEL_RESIDUE && key == VIAL_KEY_RESID) return VIAL_KEY_ID;
    if (level == VIAL_LEVEL_RESIDUE && key == VIAL_KEY_RESNAME) return VIAL_KEY_NAME;
    if (level == VIAL_LEVEL_PARTICLE && key == VIAL_KEY_ID) return VIAL_KEY_RESID;
    return VIAL_KEY_NONE;
}

static const vial_filter_t* find_filter(const vial_select_t* q, uint8_t key) {
    for (size_t i = 0; i < md_array_size(q->filters); ++i) if (q->filters[i].key == key) return &q->filters[i];
    return NULL;
}

static T selector_call(infer_t* x, vial_node_id call, str_t name) {
    uint8_t pkind = 0;
    const uint8_t level = selector_level(name, &pkind);
    vial_select_t q = {.level = level, .particle_kind = pkind};
    const char* sel = vial_select_name(level);
    T stage = vial_type_scalar(VIAL_DT_SET);
    bool failed = false, known = x->sys != NULL;
    uint32_t num_pos = 0;
    vial_node_id resid_arg = 0;
    const uint32_t n = vx_num_children(x, call);

    if (level == VIAL_LEVEL_CHAIN && no_chains(x)) { report_no_chains(x, call); failed = true; }

    for (uint32_t i = 1; i < n; ++i) {
        const vial_node_id arg = vx_child(x, call, i);
        if (!arg) continue;
        if (vx_kind(x, arg) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, arg, 0), vn = vx_child(x, arg, 1);
            const T vt = vn ? vx_infer(x, vn) : err();
            if (!kn || !vn || vial_type_is_error(&vt)) { failed = true; continue; }
            const str_t kw = vx_text(x, kn);
            const uint8_t key = vial_select_key(kw);
            if (!key || !vial_select_takes(level, key)) {
                const uint8_t instead = key_instead(level, key);
                char fix[96] = "";
                if (instead) {
                    char ex[VIAL_EXCERPT_MAX + 4];
                    const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, vn));
                    snprintf(fix, sizeof(fix), "%s=%.*s", key_word(instead), (int)e.len, e.ptr);
                }
                vx_diag_fix(x, "E0202", arg, fix[0] ? fix : NULL, "`%s` takes no keyword `%.*s`; it takes %s", sel, (int)kw.len, kw.ptr,
                    vial_select_keywords(level));
                failed = true;
                continue;
            }
            join_stage(&stage, &vt);
            if (key == VIAL_KEY_SCOPE) {
                if (vt.dtype != VIAL_DT_SET) {
                    vx_diag_fix(x, "E0206", vn, NULL, "`scope` takes a set or system; this is %s", vx_tstr(&vt));
                    failed = true;
                    continue;
                }
                const vial_value_t v = vx_value_of(x, vn);
                if (v.kind == VIAL_VAL_SETS) q.scope = v.sets->single ? v.sets : vial_sets_merge(v.sets, x->alloc);
                else known = false;
                continue;
            }
            if (key == VIAL_KEY_CHAIN && no_chains(x)) { report_no_chains(x, arg); failed = true; continue; }
            vial_filter_t f;
            const int r = filter_of(x, level, key, vn, &f);
            if (r < 0) { failed = true; continue; }
            if (r == 0) known = false;
            if (key == VIAL_KEY_SERIAL) known = false;      // mdlib keeps no serial numbers yet
            if (key == VIAL_KEY_RESID) resid_arg = arg;
            md_array_push(q.filters, f, x->alloc);
            continue;
        }
        // A positional argument: at most one
        const T vt = vx_infer(x, arg);
        num_pos += 1;
        if (num_pos > 1) {
            char fix[96] = "";
            char ex[VIAL_EXCERPT_MAX + 4];
            const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, arg));
            if (vt.dtype == VIAL_DT_STR && vial_select_takes(level, level == VIAL_LEVEL_CHAIN ? VIAL_KEY_ID : VIAL_KEY_NAME)) {
                snprintf(fix, sizeof(fix), "%s=%.*s", level == VIAL_LEVEL_CHAIN ? "id" : "name", (int)e.len, e.ptr);
            } else if ((vt.dtype == VIAL_DT_INT || vt.dtype == VIAL_DT_RANGE) && (level == VIAL_LEVEL_RESIDUE || level == VIAL_LEVEL_PARTICLE)) {
                snprintf(fix, sizeof(fix), "%s=%.*s", level == VIAL_LEVEL_RESIDUE ? "id" : "serial", (int)e.len, e.ptr);
            }
            vx_diag_fix(x, "E0203", arg, fix[0] ? fix : NULL, "`%s` takes one positional argument; `%.*s` is the second", sel, (int)e.len, e.ptr);
            failed = true;
            continue;
        }
        if (vial_type_is_error(&vt)) { failed = true; continue; }
        if (!(level == VIAL_LEVEL_PARTICLE && vt.ndim >= 1 && vt.dtype != VIAL_DT_STR && vt.dtype != VIAL_DT_INT && vt.dtype != VIAL_DT_SET)) {
            join_stage(&stage, &vt);    // particles(x) on an array reads only what x runs over, which is static
        }
        if (vt.dtype == VIAL_DT_SET) {
            const vial_value_t v = vx_value_of(x, arg);
            if (v.kind == VIAL_VAL_SETS) q.within = v.sets;
            else known = false;
            continue;
        }
        if (level == VIAL_LEVEL_PARTICLE && vt.ndim >= 1 && vt.dtype != VIAL_DT_STR && vt.dtype != VIAL_DT_INT) {
            // particles(x) on an array: the particles its first particle dimension describes
            int pd = -1;
            for (uint8_t k = 0; k < vt.ndim && pd < 0; ++k) {
                const uint32_t d = vt.dims[k].desc;
                if (d && d < md_array_size((*x->descs)) && (*x->descs)[d].level == VIAL_LEVEL_PARTICLE) pd = k;
            }
            if (pd < 0) {
                vx_diag_fix(x, "E0311", arg, NULL, "particles(x) takes the particles an array runs over, and `%s` has no dimension over particles",
                    vx_tstr(&vt));
                failed = true;
                continue;
            }
            const vial_desc_t* d = &(*x->descs)[vt.dims[pd].desc];
            if (d->sets && vt.dims[pd].kind == VIAL_EXT_STATIC) q.within = d->sets;
            else known = false;
            continue;
        }
        const bool is_str = vt.dtype == VIAL_DT_STR;
        const bool is_int = vt.dtype == VIAL_DT_INT || vt.dtype == VIAL_DT_RANGE;
        uint8_t key = VIAL_KEY_NONE;
        if (is_str && level == VIAL_LEVEL_CHAIN) key = VIAL_KEY_ID;
        else if (is_str && vial_select_takes(level, VIAL_KEY_NAME)) key = VIAL_KEY_NAME;
        else if (is_int && level == VIAL_LEVEL_RESIDUE) key = VIAL_KEY_ID;
        else if (is_int && level == VIAL_LEVEL_PARTICLE) key = VIAL_KEY_SERIAL;
        if (!key) {
            char fix[96] = "";
            const char* what = level == VIAL_LEVEL_CHAIN ? "an id string" : level == VIAL_LEVEL_PART ? "a set; parts have no name or id" :
                               level == VIAL_LEVEL_MOLECULE ? "a set or a name" : "a set, a name or an id";
            if (is_int && vt.dtype == VIAL_DT_INT && vx_num_children(x, call) == 2) {
                char ex[VIAL_EXCERPT_MAX + 4];
                const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, arg));
                snprintf(fix, sizeof(fix), "%.*s[%.*s]", (int)name.len, name.ptr, (int)e.len, e.ptr);
            }
            vx_diag_fix(x, "E0205", fix[0] ? call : arg, fix[0] ? fix : NULL, "`%s` takes %s, not %s%s", sel, what, vx_tstr(&vt),
                is_int ? "; positions are subscripts" : "");
            failed = true;
            continue;
        }
        vial_filter_t f;
        const int r = filter_of(x, level, key, arg, &f);
        if (r < 0) { failed = true; continue; }
        if (r == 0 || key == VIAL_KEY_SERIAL) known = false;
        md_array_push(q.filters, f, x->alloc);
    }
    if (failed) return err();
    if (stage.stage != VIAL_STAGE_STATIC) return selector_unknown(level, stage.stage, call);
    if (!known) return selector_unknown(level, VIAL_STAGE_STATIC, call);

    const vial_sets_t* s = vial_select(x->sys, &q, x->alloc);
    const bool empty = s->single ? vial_sets_size(s, 0) == 0 : s->count == 0;
    if (empty) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, call));
        vx_diag_fix(x, "W0601", call, NULL, "`%.*s` matches nothing in its scope", (int)e.len, e.ptr);
    } else if (level == VIAL_LEVEL_PARTICLE && resid_arg && !find_filter(&q, VIAL_KEY_CHAIN) && x->sys->num_chains > 1) {
        int32_t chains[2];
        if (vial_select_resid_chains(x->sys, find_filter(&q, VIAL_KEY_RESID), q.scope, chains) == 2) {
            const str_t a = vial_sys_chain_id(x->sys, chains[0]), b = vial_sys_chain_id(x->sys, chains[1]);
            const str_t at = vx_text(x, resid_arg);
            char fix[128];
            snprintf(fix, sizeof(fix), "%.*s, chain='%.*s'", (int)at.len, at.ptr, (int)a.len, a.ptr);
            char ex[VIAL_EXCERPT_MAX + 4];
            const str_t e = vial_excerpt(ex, sizeof(ex), at);
            vx_diag_fix(x, "W0602", resid_arg, at.len < 96 ? fix : NULL, "`%.*s` matches residues in chains %.*s and %.*s; the merged set hides which",
                (int)e.len, e.ptr, (int)a.len, a.ptr, (int)b.len, b.ptr);
        }
    }
    T t = vx_sets_type(x, s, call);
    vx_set_value(x, call, sets_value(s));
    return t;
}

bool vx_name_sets(infer_t* x, vial_node_id id, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    if (b->flags & VIAL_BUILTIN_SELECTOR) {
        uint8_t pkind = 0;
        const uint8_t level = selector_level(name, &pkind);
        if (level == VIAL_LEVEL_CHAIN && no_chains(x)) { report_no_chains(x, id); *out = err(); return true; }
        if (!x->sys) { *out = selector_unknown(level, VIAL_STAGE_STATIC, id); return true; }
        const vial_select_t q = {.level = level, .particle_kind = pkind};
        const vial_sets_t* s = vial_select(x->sys, &q, x->alloc);
        *out = vx_sets_type(x, s, id);
        vx_set_value(x, id, sets_value(s));
        return true;
    }
    if (b->flags & VIAL_BUILTIN_VALUE) {
        if (!x->sys) { *out = vial_type_scalar(VIAL_DT_SET); return true; }
        const vial_sets_t* s = NULL;
        size_t k = 0;
        for (; k < ARRAY_SIZE(x->predefined) && x->predefined[k].sets; ++k) {
            if (str_eq(x->predefined[k].name, name)) { s = x->predefined[k].sets; break; }
        }
        if (!s) {
            s = vial_select_predefined(x->sys, name, x->alloc);
            if (!s) return false;
            if (k < ARRAY_SIZE(x->predefined)) { x->predefined[k].name = name; x->predefined[k].sets = s; }
        }
        *out = vx_sets_type(x, s, id);
        vx_set_value(x, id, sets_value(s));
        return true;
    }
    return false;
}

// ### SET ALGEBRA ###

vial_value_t vx_sets_boolop(infer_t* x, bool is_and, const vial_value_t* a, const vial_value_t* b) {
    if (a->kind != VIAL_VAL_SETS || b->kind != VIAL_VAL_SETS || !a->sets->single || !b->sets->single) return (vial_value_t){0};
    return sets_value(is_and ? vial_sets_and(a->sets, 0, b->sets, 0, x->alloc) : vial_sets_or(a->sets, 0, b->sets, 0, x->alloc));
}

vial_value_t vx_sets_not(infer_t* x, const vial_value_t* a) {
    if (a->kind != VIAL_VAL_SETS || !a->sets->single || !x->sys) return (vial_value_t){0};
    if (!x->all) x->all = vial_sets_all(x->sys->num_atoms, x->alloc);
    return sets_value(vial_sets_not(x->all, a->sets, 0, x->alloc));
}

// ### INDEXING ###

// Python's slice over n: the positions it takes, in order
static uint32_t* slice_positions(infer_t* x, vial_node_id sl, uint32_t n, uint32_t* count, bool* ok) {
    int64_t lo = 0, hi = 0, step = 1;
    const vial_node_id ln = vx_child(x, sl, 0), hn = vx_child(x, sl, 1), sn = vx_child(x, sl, 2);
    *ok = true;
    if (ln) { const vial_value_t v = vx_value_of(x, ln); if (v.kind != VIAL_VAL_INT) *ok = false; else lo = v.i; }
    if (hn) { const vial_value_t v = vx_value_of(x, hn); if (v.kind != VIAL_VAL_INT) *ok = false; else hi = v.i; }
    if (sn) { const vial_value_t v = vx_value_of(x, sn); if (v.kind != VIAL_VAL_INT || v.i == 0) *ok = false; else step = v.i; }
    *count = 0;
    if (!*ok) return NULL;
    const int64_t N = n;
    int64_t b, e;
    if (step > 0) {
        b = !ln ? 0 : (lo < 0 ? (lo + N < 0 ? 0 : lo + N) : (lo > N ? N : lo));
        e = !hn ? N : (hi < 0 ? (hi + N < 0 ? 0 : hi + N) : (hi > N ? N : hi));
    } else {
        b = !ln ? N - 1 : (lo < 0 ? (lo + N < -1 ? -1 : lo + N) : (lo >= N ? N - 1 : lo));
        e = !hn ? -1 : (hi < 0 ? (hi + N < -1 ? -1 : hi + N) : (hi >= N ? N - 1 : hi));
    }
    uint32_t* pos = NULL;
    for (int64_t k = b; step > 0 ? k < e : k > e; k += step) md_array_push(pos, (uint32_t)k, x->alloc);
    *count = (uint32_t)md_array_size(pos);
    return pos;
}

// The positions an index item takes along an extent of n: one int, a slice, or a list of ints. single says the item
// removes the dimension.
static uint32_t* positions(infer_t* x, vial_node_id item, uint32_t n, uint32_t* count, bool* single) {
    *count = 0;
    *single = false;
    if (vx_kind(x, item) == VIAL_N_SLICE) {
        bool ok;
        uint32_t* p = slice_positions(x, item, n, count, &ok);
        return ok ? p : NULL;
    }
    const vial_value_t v = vx_value_of(x, item);
    uint32_t* pos = NULL;
    if (v.kind == VIAL_VAL_INT) {
        const int64_t k = v.i < 0 ? v.i + n : v.i;
        if (k < 0 || k >= n) return NULL;
        md_array_push(pos, (uint32_t)k, x->alloc);
        *single = true;
        *count = 1;
        return pos;
    }
    if (vx_kind(x, item) == VIAL_N_LIST) {
        for (uint32_t i = 0; i < vx_num_children(x, item); ++i) {
            const vial_value_t c = vx_value_of(x, vx_child(x, item, i));
            if (c.kind != VIAL_VAL_INT) return NULL;
            const int64_t k = c.i < 0 ? c.i + n : c.i;
            if (k < 0 || k >= n) return NULL;
            md_array_push(pos, (uint32_t)k, x->alloc);
        }
        *count = (uint32_t)md_array_size(pos);
        return pos ? pos : md_alloc(x->alloc, sizeof(uint32_t));
    }
    return NULL;
}

vial_value_t vx_sets_subscript(infer_t* x, const vial_value_t* a, vial_node_id item) {
    vial_value_t none = {0};
    if (!item) return none;
    if (a->kind == VIAL_VAL_LABELS) {
        // One label is written when it is taken; a slice keeps them unwritten
        uint32_t cnt; bool single;
        const uint32_t* pos = positions(x, item, a->sets->count, &cnt, &single);
        if (!pos || !x->sys) return none;
        if (single) {
            char buf[64];
            const size_t len = vial_unit_label(buf, sizeof(buf), x->sys, a->sets->level, a->sets->unit[pos[0]]);
            return (vial_value_t){.kind = VIAL_VAL_STR, .s = str_copy((str_t){buf, len}, x->alloc)};
        }
        return (vial_value_t){.kind = VIAL_VAL_LABELS, .sets = vial_sets_pick(a->sets, pos, cnt, x->alloc)};
    }
    if (a->kind == VIAL_VAL_STRS) {
        uint32_t cnt; bool single;
        const uint32_t* pos = positions(x, item, a->n, &cnt, &single);
        if (!pos) return none;
        if (single) return (vial_value_t){.kind = VIAL_VAL_STR, .s = a->strs[pos[0]]};
        str_t* strs = md_alloc(x->alloc, (cnt + 1) * sizeof(str_t));
        for (uint32_t k = 0; k < cnt; ++k) strs[k] = a->strs[pos[k]];
        return (vial_value_t){.kind = VIAL_VAL_STRS, .strs = strs, .n = cnt};
    }
    if (a->kind != VIAL_VAL_SETS) return none;
    const vial_sets_t* s = a->sets;
    uint32_t cnt; bool single;
    if (s->single) {
        const uint32_t* pos = positions(x, item, vial_sets_size(s, 0), &cnt, &single);
        if (!pos) return none;
        return sets_value(vial_sets_pick_particles(s, 0, pos, cnt, x->alloc));
    }
    const uint32_t* pos = positions(x, item, s->count, &cnt, &single);
    if (!pos) return none;
    vial_sets_t* r = vial_sets_pick(s, pos, cnt, x->alloc);
    r->single = single;
    return sets_value(r);
}

// ### len, count, labels ###

static vial_node_id only_positional(infer_t* x, vial_node_id call, const char* fname, vial_node_id* kw_axis) {
    vial_node_id arg = 0;
    uint32_t num = 0;
    if (kw_axis) *kw_axis = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, a, 0);
            if (kw_axis && kn && str_eq(vx_text(x, kn), STR_LIT("axis"))) { *kw_axis = vx_child(x, a, 1); continue; }
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            if (vx_child(x, a, 1)) vx_infer(x, vx_child(x, a, 1));
            vx_diag_fix(x, "E0202", a, NULL, "`%s` takes no keyword `%.*s`", fname, (int)kw.len, kw.ptr);
            return 0;
        }
        num += 1;
        if (num == 1) arg = a;
    }
    if (num != 1) {
        for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
            const vial_node_id a = vx_child(x, call, i);
            if (a && vx_kind(x, a) != VIAL_N_KWARG) vx_infer(x, a);
        }
        vx_diag_fix(x, "E0204", call, NULL, "`%s` takes one argument; it got %u", fname, num);
        return 0;
    }
    return arg;
}

static T call_len(infer_t* x, vial_node_id call) {
    const vial_node_id arg = only_positional(x, call, "len", NULL);
    if (!arg) return err();
    const T a = vx_infer(x, arg);
    if (vial_type_is_error(&a)) return err();
    if (a.ndim == 0 && a.dtype != VIAL_DT_SET && a.dtype != VIAL_DT_TUPLE) {
        vx_diag_fix(x, "E0313", arg, NULL, "`len` takes an array or a set; this is %s", vx_tstr(&a));
        return err();
    }
    T t = vial_type_scalar(VIAL_DT_INT);
    t.stage = a.stage; t.pass = a.pass;
    if (a.dtype == VIAL_DT_TUPLE) {
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_INT, .i = a.num_elems});
    } else if (a.ndim == 0) {
        const vial_value_t v = vx_value_of(x, arg);
        if (v.kind == VIAL_VAL_SETS && a.stage == VIAL_STAGE_STATIC) vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_INT, .i = vial_sets_size(v.sets, 0)});
    } else if (a.dims[0].kind == VIAL_EXT_STATIC) {
        t.stage = VIAL_STAGE_STATIC; t.pass = 0;
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_INT, .i = a.dims[0].n});
    } else if (a.dims[0].kind == VIAL_EXT_UNKNOWN) {
        t.stage = VIAL_STAGE_STATIC; t.pass = 0;
    }
    return t;
}

static T call_count(infer_t* x, vial_node_id call) {
    const vial_node_id arg = only_positional(x, call, "count", NULL);
    if (!arg) return err();
    const T a = vx_infer(x, arg);
    if (vial_type_is_error(&a)) return err();
    if (a.dtype != VIAL_DT_SET) {
        vx_diag_fix(x, "E0313", arg, a.dtype == VIAL_DT_BOOL ? "sum(...)" : NULL,
            "`count` takes a set or an array of sets; this is %s%s", vx_tstr(&a), a.dtype == VIAL_DT_BOOL ? ", and sum(mask) counts the true values" : "");
        return err();
    }
    T t = a;
    t.dtype = VIAL_DT_INT;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY | VIAL_TYPE_OVERLAP);
    t.unit = md_unit_none();
    t.unit_kind = VIAL_UNIT_KNOWN;
    const vial_value_t v = vx_value_of(x, arg);
    if (a.ndim == 0 && v.kind == VIAL_VAL_SETS && a.stage == VIAL_STAGE_STATIC) {
        vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_INT, .i = vial_sets_size(v.sets, 0)});
    }
    return t;
}

static T call_labels(infer_t* x, vial_node_id call) {
    vial_node_id axis_node = 0;
    const vial_node_id arg = only_positional(x, call, "labels", &axis_node);
    if (!arg) return err();
    const T a = vx_infer(x, arg);
    int64_t axis = 0;
    if (axis_node) {
        const T at = vx_infer(x, axis_node);
        if (vial_type_is_error(&at)) return err();
        const vial_value_t av = vx_value_of(x, axis_node);
        if (av.kind != VIAL_VAL_INT) {
            vx_diag_fix(x, "E0316", axis_node, NULL, "`axis` is an int known when compiling; this is %s", vx_tstr(&at));
            return err();
        }
        axis = av.i;
    }
    if (vial_type_is_error(&a)) return err();
    if (axis < 0) axis += a.ndim;
    if (a.ndim == 0 || axis < 0 || axis >= a.ndim) {
        vx_diag_fix(x, "E0316", axis_node ? axis_node : arg, NULL, "`labels` takes a dimension of %s, which has %u", vx_tstr(&a), a.ndim);
        return err();
    }
    const vial_extent_t e = a.dims[axis];
    T t = vial_type_scalar(VIAL_DT_STR);
    vial_type_push_dim(&t, e);
    t.stage = e.kind == VIAL_EXT_DYNAMIC ? a.stage : VIAL_STAGE_STATIC;
    // The labels themselves, when the dimension describes known units
    if (e.desc && e.desc < md_array_size((*x->descs)) && x->sys) {
        const vial_desc_t* d = &(*x->descs)[e.desc];
        if (d->sets && d->level != VIAL_LEVEL_NONE && d->sets->unit && e.kind == VIAL_EXT_STATIC && d->sets->count == e.n) {
            vx_set_value(x, call, (vial_value_t){.kind = VIAL_VAL_LABELS, .sets = d->sets});
        }
    }
    return t;
}

// within(r, of=x): the scope's particles within r of x, per frame. of= is always written.
static T call_within(infer_t* x, vial_node_id call) {
    vial_node_id r = 0, of = 0;
    uint32_t num_pos = 0;
    bool failed = false;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, a, 0), vn = vx_child(x, a, 1);
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            if (str_eq(kw, STR_LIT("of"))) { of = vn; continue; }
            if (str_eq(kw, STR_LIT("r")) && !r) { r = vn; continue; }
            if (vn) vx_infer(x, vn);
            vx_diag_fix(x, "E0202", a, NULL, "`within` takes no keyword `%.*s`; it takes r and of", (int)kw.len, kw.ptr);
            failed = true;
            continue;
        }
        num_pos += 1;
        if (num_pos == 1) { r = a; continue; }
        vx_infer(x, a);
        char fix[96] = "";
        char ex[VIAL_EXCERPT_MAX + 4];
        const str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, a));
        if (num_pos == 2 && !of) snprintf(fix, sizeof(fix), "of=%.*s", (int)e.len, e.ptr);
        vx_diag_fix(x, "E0203", a, fix[0] ? fix : NULL, "`within` takes the radius as its one positional argument; of= is always written, of=%.*s",
            (int)e.len, e.ptr);
        failed = true;
    }
    const T rt = r ? vx_infer(x, r) : err();
    const T ot = of ? vx_infer(x, of) : err();
    if (failed) return err();
    if (!r || !of) {
        vx_diag_fix(x, "E0204", call, NULL, "`within` takes a radius and of=, such as within(5.0, of=protein)");
        return err();
    }
    if (vial_type_is_error(&rt) || vial_type_is_error(&ot)) return err();
    md_unit_t length;
    vial_unit_name(STR_LIT("angstrom"), &length);
    const bool number = rt.dtype == VIAL_DT_INT || rt.dtype == VIAL_DT_FLOAT || rt.dtype == VIAL_DT_FLOAT64;
    if (!number || rt.ndim || (rt.unit_kind == VIAL_UNIT_KNOWN && !md_unit_base_equal(rt.unit, length))) {
        vx_diag_fix(x, "E0313", r, NULL, "`within` takes a radius, one length; this is %s", vx_tstr(&rt));
        return err();
    }
    const bool points = (ot.dtype == VIAL_DT_FLOAT || ot.dtype == VIAL_DT_FLOAT64) && ot.ndim >= 1 &&
                        ot.dims[ot.ndim - 1].kind == VIAL_EXT_STATIC && ot.dims[ot.ndim - 1].n == 3;
    if (ot.dtype != VIAL_DT_SET && !points) {
        vx_diag_fix(x, "E0313", of, NULL, "`of` takes a set or points, a last dimension of 3; this is %s", vx_tstr(&ot));
        return err();
    }
    T t = vial_type_scalar(VIAL_DT_SET);
    t.stage = VIAL_STAGE_FRAME;
    t.ident = (uint32_t)call;
    if (rt.stage > t.stage) t.stage = rt.stage;
    if (ot.stage > t.stage) t.stage = ot.stage;
    return t;
}

// bonds(s=None): per bond with both ends in s, its two particles as one-particle sets, set[n][2]
static T call_bonds(infer_t* x, vial_node_id call) {
    static const vx_param_t p[] = {{"s", VX_KIND_WHOLE, false}, {"kind", VX_KIND_WHOLE, false}};
    vial_node_id a[2];
    if (!vx_bind(x, call, "bonds", p, 2, a)) return err();
    T st = a[0] ? x->out->types[(a[0]) - x->out->base] : vial_type_scalar(VIAL_DT_NONE);
    if (vial_type_is_error(&st)) return err();
    if (st.dtype != VIAL_DT_SET && st.dtype != VIAL_DT_NONE) {
        vx_diag_fix(x, "E0313", a[0], NULL, "`bonds` takes a set; this is %s", vx_tstr(&st));
        return err();
    }
    if (st.dtype == VIAL_DT_SET && st.ndim) {
        vx_diag_fix(x, "E0313", a[0], NULL, "`bonds` takes one set; particles(...) merges an array of them");
        return err();
    }
    T t = vial_type_scalar(VIAL_DT_SET);
    t.flags |= VIAL_TYPE_ONE_PARTICLE;
    t.stage = st.stage;
    const vial_value_t sv = a[0] ? vx_value_of(x, a[0]) : (vial_value_t){0};
    const bool all = !a[0] || st.dtype == VIAL_DT_NONE;
    if (a[1] || !x->sys || st.stage != VIAL_STAGE_STATIC || (!all && sv.kind != VIAL_VAL_SETS)) {
        // kind= arrives with md_bond's kinds; a set per frame, or no system, gives a count not known here
        vial_type_push_dim(&t, (vial_extent_t){.kind = st.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
        vial_type_push_dim(&t, vial_extent_static(2));
        return t;
    }
    const md_system_t* sys = x->sys->sys;
    const size_t nb = sys->bond.count;
    const uint32_t* in = all ? NULL : vial_sets_items(sv.sets, 0);
    const uint32_t nin = all ? 0 : vial_sets_size(sv.sets, 0);
    uint32_t* idx = NULL;
    for (size_t i = 0; i < nb; ++i) {
        const md_atom_pair_t pr = sys->bond.pairs[i];
        if (pr.idx[0] < 0 || pr.idx[1] < 0) continue;
        if (!all) {
            // Both ends in s, by binary search in its ascending list
            bool both = true;
            for (int e = 0; e < 2 && both; ++e) {
                uint32_t lo = 0, hi = nin;
                const uint32_t want = (uint32_t)pr.idx[e];
                while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (in[mid] < want) lo = mid + 1; else hi = mid; }
                both = lo < nin && in[lo] == want;
            }
            if (!both) continue;
        }
        md_array_push(idx, (uint32_t)pr.idx[0], x->alloc);
        md_array_push(idx, (uint32_t)pr.idx[1], x->alloc);
    }
    const uint32_t n2 = (uint32_t)md_array_size(idx);
    vial_sets_t* r = md_alloc(x->alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = n2;
    r->idx = idx ? idx : md_alloc(x->alloc, sizeof(uint32_t));
    r->offset = md_alloc(x->alloc, (n2 + 1) * sizeof(uint32_t));
    for (uint32_t k = 0; k <= n2; ++k) r->offset[k] = k;
    vial_type_push_dim(&t, vial_extent_static(n2 / 2));
    vial_type_push_dim(&t, vial_extent_static(2));
    vx_set_value(x, call, sets_value(r));
    return t;
}

// spread(values, groups, combine=None): each particle of a group gets that group's value, the reverse of by=
static T call_spread(infer_t* x, vial_node_id call) {
    static const vx_param_t p[] = {{"values", VX_KIND_WHOLE, true}, {"groups", VX_KIND_WHOLE, true}, {"combine", VX_KIND_WHOLE, false}};
    vial_node_id a[3];
    if (!vx_bind(x, call, "spread", p, 3, a)) return err();
    const T v = x->out->types[(a[0]) - x->out->base], g = x->out->types[(a[1]) - x->out->base];
    if (vial_type_is_error(&v) || vial_type_is_error(&g)) return err();
    if (g.dtype != VIAL_DT_SET || g.ndim != 1) {
        vx_diag_fix(x, "E0313", a[1], NULL, "`spread` takes groups, an array of sets; this is %s", vx_tstr(&g));
        return err();
    }
    if (v.ndim == 0) {
        vx_diag_fix(x, "E0313", a[0], NULL, "`spread` takes one value per group; this is %s", vx_tstr(&v));
        return err();
    }
    vial_extent_t e;
    if (!vx_extent_join(v.dims[0], g.dims[0], &e) || (v.dims[0].kind == VIAL_EXT_STATIC && g.dims[0].kind == VIAL_EXT_STATIC && v.dims[0].n != g.dims[0].n)) {
        vx_diag_fix(x, "E0301", call, NULL, "`spread` gives each group one value, and these are %s for %s", vx_tstr(&v), vx_tstr(&g));
        return err();
    }
    bool combine = false;
    if (a[2]) {
        const vial_value_t c = vx_value_of(x, a[2]);
        if (c.kind != VIAL_VAL_STR || !(str_eq(c.s, STR_LIT("mean")) || str_eq(c.s, STR_LIT("max")) || str_eq(c.s, STR_LIT("first")))) {
            vx_diag_fix(x, "E0313", a[2], NULL, "`combine` takes 'mean', 'max' or 'first'");
            return err();
        }
        combine = true;
    }
    T t = v;
    vx_join_stage(&t, &v, &g);
    const vial_value_t gv = vx_value_of(x, a[1]);
    if (gv.kind == VIAL_VAL_SETS && !gv.sets->single && g.stage == VIAL_STAGE_STATIC) {
        const vial_sets_t* merged = vial_sets_merge(gv.sets, x->alloc);
        const uint32_t n = vial_sets_size(merged, 0);
        if (n != gv.sets->offset[gv.sets->count] && !combine) {
            vx_diag_fix(x, "E0306", call, "combine='mean'", "the groups share particles, so each of those would get several values; combine= says how");
            return err();
        }
        t.dims[0] = vial_extent_static(n);
        t.dims[0].desc = vx_desc(x, VIAL_LEVEL_PARTICLE, merged);
    } else {
        if ((g.flags & VIAL_TYPE_OVERLAP) && !combine) {
            vx_diag_fix(x, "E0306", call, "combine='mean'", "the groups may share particles, so each of those would get several values; combine= says how");
            return err();
        }
        t.dims[0] = (vial_extent_t){.kind = g.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call,
                                    .desc = vx_desc(x, VIAL_LEVEL_PARTICLE, NULL)};
    }
    return t;
}

bool vx_call_sets(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    if ((b->flags & VIAL_BUILTIN_SELECTOR)) { *out = selector_call(x, call, name); return true; }
    if (str_eq(name, STR_LIT("len")))    { *out = call_len(x, call); return true; }
    if (str_eq(name, STR_LIT("count")))  { *out = call_count(x, call); return true; }
    if (str_eq(name, STR_LIT("labels"))) { *out = call_labels(x, call); return true; }
    if (str_eq(name, STR_LIT("within"))) { *out = call_within(x, call); return true; }
    if (str_eq(name, STR_LIT("bonds")))  { *out = call_bonds(x, call); return true; }
    if (str_eq(name, STR_LIT("spread"))) { *out = call_spread(x, call); return true; }
    return false;
}

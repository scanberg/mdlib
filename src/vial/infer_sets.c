// The checker, step 4: the selectors, the predefined sets, set algebra and indexing, and len, count and labels.
//
// With a system, every static selector is evaluated when compiling: its sets are known, so their sizes, counts and
// labels are too, and W0601 (matches nothing), W0602 (a resid in several chains), E0207 (no chains) and E0208 (no
// force-field types) are given.
// Without one, the same calls are typed with extents not known. A selector with a frame-stage argument, such as a
// spatial set or a scope that changes per frame, is frame stage, its length known only per frame.

#include "infer_internal.h"

#include <md_system.h>
#include <md_match.h>
#include <md_attributes.h>
#include <md_hbond.h>
#include <md_smiles.h>
#include <md_element.h>
#include <core/md_bitset.h>
#include <core/md_allocator.h>
#include <core/md_hash.h>

#include <core/md_common.h>
#include <core/md_array.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

static uint8_t selector_level(md_str_t name, uint8_t* particle_kind) {
    *particle_kind = 0;
    if (md_str_eq(name, MD_STR_LIT("particle"))) return VIAL_LEVEL_PARTICLE;
    if (md_str_eq(name, MD_STR_LIT("atom"))) { *particle_kind = 1; return VIAL_LEVEL_PARTICLE; }
    if (md_str_eq(name, MD_STR_LIT("bead"))) { *particle_kind = 2; return VIAL_LEVEL_PARTICLE; }
    if (md_str_eq(name, MD_STR_LIT("residue"))) return VIAL_LEVEL_RESIDUE;
    if (md_str_eq(name, MD_STR_LIT("part"))) return VIAL_LEVEL_PART;
    if (md_str_eq(name, MD_STR_LIT("polymer"))) { *particle_kind = 1; return VIAL_LEVEL_PART; }  // the strands: parts of a polymer kind
    if (md_str_eq(name, MD_STR_LIT("molecule"))) return VIAL_LEVEL_MOLECULE;
    if (md_str_eq(name, MD_STR_LIT("chain"))) return VIAL_LEVEL_CHAIN;
    if (md_str_eq(name, MD_STR_LIT("tail"))) return VIAL_LEVEL_TAIL;
    if (md_str_eq(name, MD_STR_LIT("head"))) return VIAL_LEVEL_HEAD;
    if (md_str_eq(name, MD_STR_LIT("ring"))) return VIAL_LEVEL_RING;
    if (md_str_eq(name, MD_STR_LIT("group"))) return VIAL_LEVEL_GROUP;
    if (md_str_eq(name, MD_STR_LIT("fragment"))) return VIAL_LEVEL_FRAGMENT;
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
    vx_diag_fix(x, "E0207", at, "part(kind='polymer')",
        "the data gives no chain ids, and chains are never inferred; a strand is found as a part");
}

// Whether the data gives no force-field types: only a topology does, and what a loader types from the names alone has
// none (md_particle_type_ff_type)
bool vx_no_ff_types(const infer_t* x) {
    if (!x->sys) return false;
    const md_particle_type_data_t* td = &x->sys->sys->atom.type;
    for (size_t t = 0; t < td->count; ++t) if (!md_str_empty(md_particle_type_ff_type(td, t))) return false;
    return true;
}

void vx_report_no_ff_types(infer_t* x, vial_node_id at) {
    vx_diag_fix(x, "E0208", at, NULL, "the data gives no force-field types; they come from a topology (.tpr, .top, a LAMMPS data "
        "file) and are never inferred: select by name= and element=, or by hybridization=, aromatic= and hydrogens=");
}

// mdlib keeps no serial numbers yet, and a serial is not its position where the file skips records, so a particle is
// selected by its position until it does, rather than by a number that would match the wrong one
static void report_no_serials(infer_t* x, vial_node_id at) {
    vx_diag_fix(x, "E0801", at, NULL, "selecting particles by serial is not supported yet; select them by position, particle[i]");
}

// What a keyword's values are: strings, ints or bools
static uint8_t key_dtype(uint8_t level, uint8_t key) {
    switch (key) {
    case VIAL_KEY_NAME: case VIAL_KEY_ELEMENT: case VIAL_KEY_KIND: case VIAL_KEY_TYPE: case VIAL_KEY_RESNAME: case VIAL_KEY_CHAIN:
    case VIAL_KEY_ROLE: case VIAL_KEY_HYBRIDIZATION: case VIAL_KEY_GAP: return VIAL_DT_STR;
    case VIAL_KEY_AROMATIC: case VIAL_KEY_METALLIC: case VIAL_KEY_RESOLVED: return VIAL_DT_BOOL;
    case VIAL_KEY_ID: return level == VIAL_LEVEL_CHAIN ? VIAL_DT_STR : VIAL_DT_INT;
    default: return VIAL_DT_INT;
    }
}

static bool int_of(const vial_value_t* v, int64_t* out) {
    if (v->kind == VIAL_VAL_INT) { *out = v->i; return true; }
    if (v->kind == VIAL_VAL_FLOAT && v->f >= -9.0e15 && v->f <= 9.0e15 && v->f == floor(v->f)) { *out = (int64_t)v->f; return true; }
    return false;
}

// One value of a keyword, from a node already typed: a string, an int or a range. 1 when known, 0 when not known
// when compiling, -1 when it is a value the keyword does not take.
// A want of match_of and expr_of besides a dtype: an int, a range or a string, as a bond's order= takes
#define WANT_INT_OR_STR 0xF0

static int match_of(infer_t* x, vial_node_id node, uint8_t want, vial_match_t* m) {
    const T t = x->out->types[(node) - x->out->base];
    const vial_value_t v = vx_value_of(x, node);
    memset(m, 0, sizeof(*m));
    if (want == WANT_INT_OR_STR) want = t.dtype == VIAL_DT_STR ? VIAL_DT_STR : VIAL_DT_INT;
    if (want == VIAL_DT_BOOL) {
        if (t.dtype != VIAL_DT_BOOL || t.ndim) return -1;
        if (v.kind != VIAL_VAL_BOOL) return 0;
        m->kind = VIAL_MATCH_BOOL; m->lo = v.b;
        return 1;
    }
    if (want == VIAL_DT_STR) {
        if (t.dtype != VIAL_DT_STR || t.ndim) return -1;
        if (v.kind != VIAL_VAL_STR) return 0;
        m->kind = VIAL_MATCH_STR; m->s = v.s;
        return 1;
    }
    if (t.dtype == VIAL_DT_RANGE && vx_kind(x, node) == VIAL_N_RANGE) {
        const vial_value_t lo = vx_value_of(x, vx_child(x, node, 0)), hi = vx_value_of(x, vx_child(x, node, 1));
        if (!int_of(&lo, &m->lo) || !int_of(&hi, &m->hi)) return 0;
        m->kind = VIAL_MATCH_RANGE;
        return 1;
    }
    if (t.dtype == VIAL_DT_RANGE) return 0;
    if (t.dtype != VIAL_DT_INT || t.ndim) return -1;
    if (!int_of(&v, &m->lo)) return 0;
    m->kind = VIAL_MATCH_INT;
    return 1;
}

static bool is_logic(const infer_t* x, vial_node_id n) {
    const vial_node_kind_t k = vx_kind(x, n);
    return k == VIAL_N_BOOLOP || k == VIAL_N_NOT;
}

// The leaves under and/or/not: false when one is no value a keyword takes (a set, a float), or failed
static bool logic_leaves_are_values(infer_t* x, vial_node_id n, bool* any_str, uint8_t* stage) {
    if (is_logic(x, n)) {
        bool ok = true;
        for (uint32_t i = 0; i < vx_num_children(x, n); ++i) ok = logic_leaves_are_values(x, vx_child(x, n, i), any_str, stage) && ok;
        return ok;
    }
    const T t = vx_infer(x, n);
    if (vial_type_is_error(&t)) return false;
    const bool value = (t.ndim == 0 && (t.dtype == VIAL_DT_STR || t.dtype == VIAL_DT_INT || t.dtype == VIAL_DT_RANGE)) ||
                       (t.ndim == 1 && vx_kind(x, n) == VIAL_N_LIST && (t.dtype == VIAL_DT_STR || t.dtype == VIAL_DT_INT));
    *any_str = *any_str || t.dtype == VIAL_DT_STR;
    if (t.stage > *stage) *stage = t.stage;
    return value;
}

static void type_logic(infer_t* x, vial_node_id n, const T* t) {
    if (!is_logic(x, n)) return;
    x->out->types[(n) - x->out->base] = *t;
    x->out->typed[(n) - x->out->base] = 1;
    for (uint32_t i = 0; i < vx_num_children(x, n); ++i) type_logic(x, vx_child(x, n, i), t);
}

// Whether a keyword's value is and/or/not over the values it takes, as name='H*' and not 'HA', and if so its type. Its
// leaves are typed as any value is; the and/or/not over them stand for the values they accept, so they are typed as
// their leaves, a str or an int, and are no expression of their own (strings have no and/or, E0319): nothing of them
// is lowered. and/or over sets stays their intersection and union, which is no keyword's value.
static bool value_logic(infer_t* x, vial_node_id n, T* out) {
    if (!is_logic(x, n)) return false;
    bool any_str = false;
    uint8_t stage = VIAL_STAGE_STATIC;
    if (!logic_leaves_are_values(x, n, &any_str, &stage)) return false;
    T t = vial_type_scalar(any_str ? VIAL_DT_STR : VIAL_DT_INT);
    t.stage = stage;
    type_logic(x, n, &t);
    *out = t;
    return true;
}

// A keyword's value as its expression, appended to *e in prefix order (vial_filter_t): one value, a list (an or of its
// values), or and/or/not over them. 1 known, 0 not known when compiling, -1 a value the keyword does not take, whose
// node *bad receives.
static int expr_of(infer_t* x, vial_node_id node, uint8_t want, vial_match_t** e, vial_node_id* bad) {
    const T t = x->out->types[(node) - x->out->base];
    const vial_node_kind_t k = vx_kind(x, node);
    const bool list = k == VIAL_N_LIST && t.ndim == 1;
    if (is_logic(x, node) || list) {
        const uint32_t n = vx_num_children(x, node);
        const uint8_t mk = list ? VIAL_MATCH_OR : k == VIAL_N_NOT ? VIAL_MATCH_NOT : x->t->nodes[node].sub == VIAL_KW_AND ? VIAL_MATCH_AND : VIAL_MATCH_OR;
        md_array_push(*e, ((vial_match_t){.kind = mk, .num = n}), x->alloc);
        int known = 1;
        for (uint32_t i = 0; i < n; ++i) {
            const int r = expr_of(x, vx_child(x, node, i), want, e, bad);
            if (r < 0) return -1;
            if (r == 0) known = 0;
        }
        return known;
    }
    if (t.ndim == 1 && (t.dtype == want || (want == WANT_INT_OR_STR && (t.dtype == VIAL_DT_STR || t.dtype == VIAL_DT_INT)))) {
        return 0;   // a list from elsewhere: typed, not evaluated here
    }
    vial_match_t m;
    const int r = match_of(x, node, want, &m);
    if (r < 0) { *bad = node; return -1; }
    md_array_push(*e, m, x->alloc);
    return r;
}

// A keyword's value: one value, a list of them, or and/or/not over them. Typed already. 1 known, 0 not known, -1
// refused (reported).
static int filter_of(infer_t* x, uint8_t level, uint8_t key, vial_node_id node, vial_filter_t* f) {
    const uint8_t want = key_dtype(level, key);
    f->key = key;
    f->expr = NULL;
    vial_node_id bad = node;
    const int known = expr_of(x, node, want, &f->expr, &bad);
    if (known < 0) {
        const T bt = x->out->types[(bad) - x->out->base];
        vx_diag_fix(x, "E0206", bad, NULL, "`%s` takes %s; this is %s", vial_select_key_word(key),
            want == VIAL_DT_STR ? "a string, a list of strings, or and/or/not over strings" :
            want == VIAL_DT_BOOL ? "True or False" : "an int, a range a..b, a list of ints, or and/or/not over them",
            vx_tstr(&bt));
        return -1;
    }
    if (key == VIAL_KEY_HYBRIDIZATION && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_hybridization_known(f->expr[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), f->expr[i].s);
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is not a hybridization; the hybridizations are sp, sp2 and sp3", (int)e.len, e.ptr);
            return -1;
        }
    }
    if (key == VIAL_KEY_KIND && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_kind_known(level, f->expr[i].s)) continue;
            // A kind that starts the same, such as amino for amino_acid, is offered as the fix
            const char* kinds = vial_select_kinds(level);
            char fix[64] = "";
            const md_str_t w = f->expr[i].s;
            for (const char* p = kinds; *p && w.len;) {
                const char* e = p;
                while (*e && *e != ',' && *e != ' ') ++e;
                if ((size_t)(e - p) > w.len && strncmp(p, w.ptr, w.len) == 0) { snprintf(fix, sizeof(fix), "'%.*s'", (int)(e - p), p); break; }
                p = e;
                while (*p == ',' || *p == ' ' || strncmp(p, "and ", 4) == 0) p += *p == 'a' ? 4 : 1;
            }
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), w);
            vx_diag_fix(x, "E0206", node, fix[0] ? fix : NULL, "'%.*s' is not a %s kind; the kinds are %s", (int)e.len, e.ptr,
                level == VIAL_LEVEL_PARTICLE ? "particle" : level == VIAL_LEVEL_RESIDUE ? "residue" : level == VIAL_LEVEL_RING ? "ring" : "part", kinds);
            return -1;
        }
    }
    if (key == VIAL_KEY_ROLE && level == VIAL_LEVEL_PARTICLE && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_particle_role_known(f->expr[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), f->expr[i].s);
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is not a role of particles; the roles are backbone, main, side, base, sugar, "
                "terminus, n_terminus, c_terminus, 5_terminus, 3_terminus and qm", (int)e.len, e.ptr);
            return -1;
        }
    }
    if (key == VIAL_KEY_ROLE && level == VIAL_LEVEL_RESIDUE && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_residue_role_known(f->expr[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), f->expr[i].s);
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is not a role of residues; the roles are terminus, n_terminus, c_terminus, 5_terminus "
                "and 3_terminus", (int)e.len, e.ptr);
            return -1;
        }
    }
    if (key == VIAL_KEY_GAP && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_gap_known(f->expr[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), f->expr[i].s);
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is no gap of a residue; a gap is before or after it", (int)e.len, e.ptr);
            return -1;
        }
    }
    if (key == VIAL_KEY_ROLE && level == VIAL_LEVEL_GROUP && known == 1) {
        for (size_t i = 0; i < md_array_size(f->expr); ++i) {
            if (f->expr[i].kind != VIAL_MATCH_STR || vial_select_role_known(f->expr[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), f->expr[i].s);
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is not a role of groups; the roles are tc, energy, acc, freeze, user1, user2, comm, "
                "compressed-x, orire-fit and qmmm", (int)e.len, e.ptr);
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

static T selector_call(infer_t* x, vial_node_id call, md_str_t name) {
    uint8_t pkind = 0;
    const uint8_t level = selector_level(name, &pkind);
    vial_select_t q = {.level = level, .particle_kind = pkind};
    const bool polymers = level == VIAL_LEVEL_PART && pkind;
    const char* sel = polymers ? "polymer" : vial_select_name(level);
    T stage = vial_type_scalar(VIAL_DT_SET);
    bool failed = false, known = x->sys != NULL;
    uint32_t num_pos = 0;
    vial_node_id resid_arg = 0;
    vial_node_id logic_at = 0;      // the first keyword value of and/or/not, which only a selection made when compiling takes
    const uint32_t n = vx_num_children(x, call);

    if (level == VIAL_LEVEL_CHAIN && no_chains(x)) { report_no_chains(x, call); failed = true; }

    for (uint32_t i = 1; i < n; ++i) {
        const vial_node_id arg = vx_child(x, call, i);
        if (!arg) continue;
        if (vx_kind(x, arg) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, arg, 0), vn = vx_child(x, arg, 1);
            T vt;
            if (vn && value_logic(x, vn, &vt)) { if (!logic_at) logic_at = vn; }
            else vt = vn ? vx_infer(x, vn) : err();
            if (!kn || !vn || vial_type_is_error(&vt)) { failed = true; continue; }
            const md_str_t kw = vx_text(x, kn);
            const uint8_t key = vial_select_key(kw);
            if (!key || !vial_select_takes(level, key)) {
                const uint8_t instead = key_instead(level, key);
                char fix[96] = "";
                if (instead) {
                    char ex[VIAL_EXCERPT_MAX + 4];
                    const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, vn));
                    snprintf(fix, sizeof(fix), "%s=%.*s", vial_select_key_word(instead), (int)e.len, e.ptr);
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
            if (key == VIAL_KEY_TYPE && vx_no_ff_types(x)) { vx_report_no_ff_types(x, arg); failed = true; continue; }
            if (key == VIAL_KEY_SERIAL) { report_no_serials(x, arg); failed = true; continue; }
            vial_filter_t f;
            const int r = filter_of(x, level, key, vn, &f);
            if (r < 0) { failed = true; continue; }
            if (r == 0) known = false;
            if (key == VIAL_KEY_RESID) resid_arg = arg;
            md_array_push(q.filters, f, x->alloc);
            continue;
        }
        // A positional argument: at most one
        T vt;
        if (value_logic(x, arg, &vt)) { if (!logic_at) logic_at = arg; }
        else vt = vx_infer(x, arg);
        num_pos += 1;
        if (num_pos > 1) {
            char fix[96] = "";
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, arg));
            if (vt.dtype == VIAL_DT_STR && vial_select_takes(level, level == VIAL_LEVEL_CHAIN ? VIAL_KEY_ID : VIAL_KEY_NAME)) {
                snprintf(fix, sizeof(fix), "%s=%.*s", level == VIAL_LEVEL_CHAIN ? "id" : "name", (int)e.len, e.ptr);
            } else if ((vt.dtype == VIAL_DT_INT || vt.dtype == VIAL_DT_RANGE) && level == VIAL_LEVEL_RESIDUE) {
                snprintf(fix, sizeof(fix), "id=%.*s", (int)e.len, e.ptr);
            }
            vx_diag_fix(x, "E0203", arg, fix[0] ? fix : NULL, "`%s` takes one positional argument; `%.*s` is the second", sel, (int)e.len, e.ptr);
            failed = true;
            continue;
        }
        if (vial_type_is_error(&vt)) { failed = true; continue; }
        if (!(level == VIAL_LEVEL_PARTICLE && vt.ndim >= 1 && vt.dtype != VIAL_DT_STR && vt.dtype != VIAL_DT_INT && vt.dtype != VIAL_DT_SET)) {
            join_stage(&stage, &vt);    // particle(x) on an array reads only what x runs over, which is static
        }
        if (vt.dtype == VIAL_DT_SET) {
            const vial_value_t v = vx_value_of(x, arg);
            if (v.kind == VIAL_VAL_SETS) q.within = v.sets;
            else known = false;
            continue;
        }
        if (level == VIAL_LEVEL_PARTICLE && vt.ndim >= 1 && vt.dtype != VIAL_DT_STR && vt.dtype != VIAL_DT_INT) {
            // particle(x) on an array: the particles its first particle dimension describes
            int pd = -1;
            for (uint8_t k = 0; k < vt.ndim && pd < 0; ++k) {
                const uint32_t d = vt.dims[k].desc;
                if (d && d < md_array_size((*x->descs)) && (*x->descs)[d].level == VIAL_LEVEL_PARTICLE) pd = k;
            }
            if (pd < 0) {
                vx_diag_fix(x, "E0311", arg, NULL, "particle(x) takes the particles an array runs over, and `%s` has no dimension over particles",
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
        if (is_int && (level == VIAL_LEVEL_RESIDUE || level == VIAL_LEVEL_PARTICLE)) {
            // A number alone says neither which count it is in nor where it starts: v1's residue(5) is the fifth
            // residue, by position from 1, and the file's residue 5 is another wherever numbering restarts per chain
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, arg));
            const char* key_name = level == VIAL_LEVEL_RESIDUE ? "id" : "serial";
            char pos[64] = "i";
            int64_t v = 0;
            bool lit = e.len > 0 && e.len < 12;
            for (size_t c = 0; c < e.len && lit; ++c) { lit = e.ptr[c] >= '0' && e.ptr[c] <= '9'; v = v * 10 + (e.ptr[c] - '0'); }
            if (lit && v > 0) snprintf(pos, sizeof(pos), "%lld", (long long)(v - 1));
            const int nl = (int)name.len;
            vx_diag_fix(x, "E0205", arg, NULL, "`%.*s` takes a set or a name; a number is the file's, %.*s(%s=%.*s), or a position, "
                "%.*s[%s], counted from 0, which is v1's %s(%.*s)", nl, name.ptr, nl, name.ptr, key_name, (int)e.len, e.ptr, nl, name.ptr,
                pos, level == VIAL_LEVEL_RESIDUE ? "residue" : "atom", (int)e.len, e.ptr);
            failed = true;
            continue;
        }
        if (!key) {
            char fix[96] = "";
            const char* what = level == VIAL_LEVEL_CHAIN ? "an id string" :
                               level == VIAL_LEVEL_PART || level == VIAL_LEVEL_MOLECULE ? "a set or a name" : "a set, a name or an id";
            if (is_int && vt.dtype == VIAL_DT_INT && vx_num_children(x, call) == 2) {
                char ex[VIAL_EXCERPT_MAX + 4];
                const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, arg));
                snprintf(fix, sizeof(fix), "%.*s[%.*s]", (int)name.len, name.ptr, (int)e.len, e.ptr);
            }
            vx_diag_fix(x, "E0205", fix[0] ? call : arg, fix[0] ? fix : NULL, "`%s` takes %s, not %s%s", sel, what, vx_tstr(&vt),
                is_int ? "; positions are subscripts" : "");
            failed = true;
            continue;
        }
        if (key == VIAL_KEY_SERIAL) { report_no_serials(x, arg); failed = true; continue; }
        vial_filter_t f;
        const int r = filter_of(x, level, key, arg, &f);
        if (r < 0) { failed = true; continue; }
        if (r == 0) known = false;
        md_array_push(q.filters, f, x->alloc);
    }
    // The kinds of a strand: polymer(kind='ligand') would pick nothing by construction
    for (size_t i = 0; polymers && !failed && i < md_array_size(q.filters); ++i) {
        const vial_filter_t* f = &q.filters[i];
        if (f->key != VIAL_KEY_KIND) continue;
        for (size_t j = 0; j < md_array_size(f->expr); ++j) {
            static const char* strands[] = {"polymer", "peptide", "nucleic", "dna", "rna"};
            if (f->expr[j].kind != VIAL_MATCH_STR) continue;
            bool ok = false;
            for (size_t k = 0; k < MD_COUNTOF(strands) && !ok; ++k) ok = md_str_eq(f->expr[j].s, md_str_from_cstr(strands[k]));
            if (ok) continue;
            vx_diag_fix(x, "E0206", call, NULL, "'%.*s' is not a kind of polymer; the kinds are polymer, peptide, nucleic, dna and rna",
                (int)f->expr[j].s.len, f->expr[j].s.ptr);
            failed = true;
            break;
        }
    }
    if (failed) return err();
    if (logic_at && (stage.stage != VIAL_STAGE_STATIC || !known)) {
        // Evaluation takes a keyword's values, and and/or/not over strings is no value
        vx_diag_fix(x, "E0801", logic_at, NULL, "and/or/not in a keyword's value is resolved when compiling; a selection %s does not take it yet",
            stage.stage != VIAL_STAGE_STATIC ? "that changes per frame" : "of values not known when compiling");
        return err();
    }
    if (stage.stage != VIAL_STAGE_STATIC) return selector_unknown(level, stage.stage, call);
    if (!known) return selector_unknown(level, VIAL_STAGE_STATIC, call);

    const vial_sets_t* s = vial_select(x->sys, &q, x->alloc);
    const bool empty = s->single ? vial_sets_size(s, 0) == 0 : s->count == 0;
    if (empty) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, call));
        vx_diag_fix(x, "W0601", call, NULL, "`%.*s` matches nothing in its scope", (int)e.len, e.ptr);
    } else if (level == VIAL_LEVEL_PARTICLE && resid_arg && !find_filter(&q, VIAL_KEY_CHAIN) && x->sys->num_chains > 1) {
        int32_t chains[2];
        if (vial_select_resid_chains(x->sys, find_filter(&q, VIAL_KEY_RESID), q.scope, chains) == 2) {
            const md_str_t a = vial_sys_chain_id(x->sys, chains[0]), b = vial_sys_chain_id(x->sys, chains[1]);
            const md_str_t at = vx_text(x, resid_arg);
            char fix[128];
            snprintf(fix, sizeof(fix), "%.*s, chain='%.*s'", (int)at.len, at.ptr, (int)a.len, a.ptr);
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), at);
            vx_diag_fix(x, "W0602", resid_arg, at.len < 96 ? fix : NULL, "`%.*s` matches residues in chains %.*s and %.*s; the merged set hides which",
                (int)e.len, e.ptr, (int)a.len, a.ptr, (int)b.len, b.ptr);
        }
    }
    T t = vx_sets_type(x, s, call);
    vx_set_value(x, call, sets_value(s));
    return t;
}

// The particles of a set of a stored selection, ascending, each once, those the system has, written at out: how many
static uint32_t stored_set(const uint32_t* atoms, size_t n_in, uint32_t na, uint32_t* idx) {
    uint32_t n = 0;
    bool sorted = true;
    for (size_t i = 0; i < n_in; ++i) {
        if (atoms[i] >= na) continue;
        if (n && atoms[i] <= idx[n - 1]) sorted = false;
        idx[n++] = atoms[i];
    }
    if (!sorted) {
        // Insertion into place, then each once: stored selections are rarely out of order
        for (uint32_t i = 1; i < n; ++i) {
            const uint32_t v = idx[i];
            uint32_t j = i;
            while (j && idx[j - 1] > v) { idx[j] = idx[j - 1]; --j; }
            idx[j] = v;
        }
        uint32_t m = 0;
        for (uint32_t i = 0; i < n; ++i) if (!m || idx[i] != idx[m - 1]) idx[m++] = idx[i];
        n = m;
    }
    return n;
}

bool vx_name_stored(infer_t* x, vial_node_id id, uint32_t stored, T* out) {
    if (!x->c || stored >= x->c->num_stored) return false;
    const md_vial_named_set_t* ns = &x->c->stored[stored];
    const uint32_t na = x->sys ? x->sys->num_atoms : UINT32_MAX;
    const uint32_t count = ns->ndim ? ns->dims[0] * (ns->ndim > 1 ? ns->dims[1] : 1) : 1;
    vial_sets_t* s = md_alloc(x->alloc, sizeof(vial_sets_t));
    memset(s, 0, sizeof(*s));
    uint32_t* off = md_alloc(x->alloc, (count + 1) * sizeof(uint32_t));
    uint32_t* idx = md_alloc(x->alloc, (ns->num_atoms ? ns->num_atoms : 1) * sizeof(uint32_t));
    off[0] = 0;
    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t b = ns->ndim ? ns->offsets[i] : 0;
        const uint32_t e = ns->ndim ? ns->offsets[i + 1] : (uint32_t)ns->num_atoms;
        off[i + 1] = off[i] + (e > b ? stored_set(ns->atoms + b, e - b, na, idx + off[i]) : 0);
    }
    s->count = count;
    s->offset = off;
    s->idx = idx;
    s->single = ns->ndim == 0;
    if (ns->ndim < 2) {
        *out = vx_sets_type(x, s, id);
    } else {
        T t = vial_type_scalar(VIAL_DT_SET);
        vial_type_push_dim(&t, vial_extent_static(ns->dims[0]));
        vial_type_push_dim(&t, vial_extent_static(ns->dims[1]));
        bool one = count > 0;
        for (uint32_t i = 0; one && i < count; ++i) one = vial_sets_size(s, i) == 1;
        if (one) t.flags |= VIAL_TYPE_ONE_PARTICLE;
        *out = t;
    }
    vx_set_value(x, id, sets_value(s));
    return true;
}

bool vx_name_sets(infer_t* x, vial_node_id id, const vial_builtin_t* b, T* out) {
    const md_str_t name = md_str_from_cstr(b->name);
    if (b->flags & VIAL_BUILTIN_SELECTOR) {
        uint8_t pkind = 0;
        const uint8_t level = selector_level(name, &pkind);
        if (level == VIAL_LEVEL_CHAIN && no_chains(x)) { report_no_chains(x, id); *out = err(); return true; }
        if (!x->sys) { *out = selector_unknown(level, VIAL_STAGE_STATIC, id); return true; }
        // All of a level's units, once per compile, as the predefined sets
        const vial_sets_t* s = NULL;
        size_t k = 0;
        for (; k < MD_COUNTOF(x->predefined) && x->predefined[k].sets; ++k) {
            if (md_str_eq(x->predefined[k].name, name)) { s = x->predefined[k].sets; break; }
        }
        if (!s) {
            const vial_select_t q = {.level = level, .particle_kind = pkind};
            s = vial_select(x->sys, &q, x->alloc);
            if (k < MD_COUNTOF(x->predefined)) { x->predefined[k].name = name; x->predefined[k].sets = s; }
        }
        *out = vx_sets_type(x, s, id);
        vx_set_value(x, id, sets_value(s));
        return true;
    }
    if (md_str_eq(name, MD_STR_LIT("pi"))) {
        // A float, as VIAL's numbers are: beside a float64 a bare float is read as one (in_operand_dtype)
        *out = vial_type_scalar(VIAL_DT_FLOAT);
        vx_set_value(x, id, (vial_value_t){.kind = VIAL_VAL_FLOAT, .f = 3.14159265358979323846f});
        return true;
    }
    if (b->flags & VIAL_BUILTIN_VALUE) {
        if (!x->sys) { *out = vial_type_scalar(VIAL_DT_SET); return true; }
        const vial_sets_t* s = NULL;
        size_t k = 0;
        for (; k < MD_COUNTOF(x->predefined) && x->predefined[k].sets; ++k) {
            if (md_str_eq(x->predefined[k].name, name)) { s = x->predefined[k].sets; break; }
        }
        if (!s) {
            s = vial_select_predefined(x->sys, name, x->alloc);
            if (!s) return false;
            if (k < MD_COUNTOF(x->predefined)) { x->predefined[k].name = name; x->predefined[k].sets = s; }
        }
        *out = vx_sets_type(x, s, id);
        vx_set_value(x, id, sets_value(s));
        return true;
    }
    return false;
}

// ### SET ALGEBRA ###

vial_value_t vx_sets_boolop(infer_t* x, bool is_and, const vial_value_t* a, const vial_value_t* b) {
    if (a->kind != VIAL_VAL_SETS || b->kind != VIAL_VAL_SETS) return (vial_value_t){0};
    // An array with one set broadcasts: residue('ALA') and particle(name='CA') is known set by set, so that each is
    // known to hold one particle where it does
    if (a->sets->single != b->sets->single) {
        const vial_sets_t* arr = a->sets->single ? b->sets : a->sets;
        const vial_sets_t* one = a->sets->single ? a->sets : b->sets;
        return sets_value(vial_sets_boolop_each(arr, one, is_and, x->alloc));
    }
    if (!a->sets->single) return (vial_value_t){0};
    return sets_value(is_and ? vial_sets_and(a->sets, 0, b->sets, 0, x->alloc) : vial_sets_or(a->sets, 0, b->sets, 0, x->alloc));
}

vial_value_t vx_sets_not(infer_t* x, const vial_value_t* a) {
    if (a->kind != VIAL_VAL_SETS || !a->sets->single || !x->sys) return (vial_value_t){0};
    if (!x->all) x->all = vial_sets_all(x->sys->num_atoms, x->alloc);
    return sets_value(vial_sets_not(x->all, a->sets, 0, x->alloc));
}

// ### NAMES ###

bool vx_name_is_glob(md_str_t name) {
    return vial_is_glob(name);
}

static vial_sets_t* select_named(infer_t* x, uint8_t level, md_str_t name) {
    vial_select_t q = {.level = level};
    vial_filter_t f = {.key = VIAL_KEY_NAME};
    md_array_push(f.expr, ((vial_match_t){.kind = VIAL_MATCH_STR, .s = name}), x->alloc);
    md_array_push(q.filters, f, x->alloc);
    return vial_select(x->sys, &q, x->alloc);
}

const vial_sets_t* vx_named_particles(infer_t* x, md_str_t name) {
    return x->sys ? select_named(x, VIAL_LEVEL_PARTICLE, name) : NULL;
}

T vx_named(infer_t* x, vial_node_id id, vial_node_id on, vial_node_id in, const T* a) {
    const T st = vx_infer(x, in);
    if (vial_type_is_error(&st)) return st;
    const vial_value_t nv = vx_value_of(x, in);
    if (nv.kind != VIAL_VAL_STR) {
        vx_diag_fix(x, "E0319", in, NULL, "a set's particles are named by a constant string");
        return vial_type_error();
    }
    const md_str_t name = nv.s;
    const bool glob = vx_name_is_glob(name);
    const bool said = vial_diag_at(x->diags, x->t->nodes[in].span.beg);    // an f-string, refused already
    T t = *a;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    if (!glob) t.flags |= VIAL_TYPE_ONE_PARTICLE;   // one particle per set by what the subscript means
    if (!x->sys) return t;
    const vial_sets_t* named = vx_named_particles(x, name);
    if (said) return vial_type_error();
    if (vial_sets_size(named, 0) == 0) {
        // residue['ALA'] names particles: a residue name is residue('ALA')
        const vial_sets_t* res = select_named(x, VIAL_LEVEL_RESIDUE, name);
        if (res->count) {
            char fix[96];
            snprintf(fix, sizeof(fix), "residue('%.*s')", (int)MD_MIN(name.len, 60), name.ptr);
            vx_diag_fix(x, "W0601", in, fix, "no particle is named '%.*s'; a subscript names particles, and residues of that name are %s",
                (int)name.len, name.ptr, fix);
        } else {
            vx_diag_fix(x, "W0601", in, NULL, "no particle is named '%.*s'", (int)name.len, name.ptr);
        }
    }
    const vial_value_t av = vx_value_of(x, on);
    if (av.kind != VIAL_VAL_SETS || t.stage != VIAL_STAGE_STATIC) return t;
    vial_sets_t* r = av.sets->single ? vial_sets_and(av.sets, 0, named, 0, x->alloc) : vial_sets_boolop_each(av.sets, named, true, x->alloc);
    if (av.sets->single) r->single = true;
    if (!glob) {
        // The sets without the name, or with it twice, are emptied: their point is missing
        uint32_t none = 0, several = 0, first_size = r->count ? vial_sets_size(r, 0) : 0;
        uint32_t* keep = md_alloc(x->alloc, (r->count + 1) * sizeof(uint32_t));
        uint32_t at = 0;
        keep[0] = 0;
        for (uint32_t i = 0; i < r->count; ++i) {
            const uint32_t sz = vial_sets_size(r, i);
            none += sz == 0;
            several += sz > 1;
            if (sz == 1) r->idx[at++] = vial_sets_items(r, i)[0];
            keep[i + 1] = at;
        }
        r->offset = keep;
        if ((none || several) && vial_sets_size(named, 0) == 0) {
            t.flags |= VIAL_TYPE_MISSING;      // no particle has the name at all: W0601 said so
        } else if (none || several) {
            t.flags |= VIAL_TYPE_MISSING;
            char what[96];
            if (none && several) snprintf(what, sizeof(what), "%u %s none and %u several", none, none == 1 ? "has" : "have", several);
            else if (none) snprintf(what, sizeof(what), "%u %s none", none, none == 1 ? "has" : "have");
            else snprintf(what, sizeof(what), "%u %s several", several, several == 1 ? "has" : "have");
            if (av.sets->single) {
                if (none) {
                    vx_diag_fix(x, "W0609", id, NULL, "the set has no particle named '%.*s', so its point is missing", (int)name.len, name.ptr);
                } else {
                    vx_diag_fix(x, "W0609", id, NULL, "the set has %u particles named '%.*s', and a name picks one per set: subscript its units, "
                        "such as residue(s)['%.*s']", first_size, (int)name.len, name.ptr, (int)name.len, name.ptr);
                }
            } else {
                vx_diag_fix(x, "W0609", id, NULL, "a particle named '%.*s' is one per set in %u of %u sets; %s, and %s missing",
                    (int)name.len, name.ptr, r->count - none - several, r->count, what, none + several == 1 ? "its entry is" : "their entries are");
            }
        }
    }
    vx_set_value(x, id, sets_value(r));
    return t;
}

// ### THE PARTICLES OF EACH SET ###

uint8_t vx_view_level(const infer_t* x, vial_node_id n, uint8_t* particle_kind, md_str_t* name) {
    *particle_kind = 0;
    if (!n || vx_kind(x, n) != VIAL_N_ATTR) return VIAL_LEVEL_NONE;
    const vial_node_id nn = vx_child(x, n, 1);
    if (!nn) return VIAL_LEVEL_NONE;
    const md_str_t t = vx_text(x, nn);
    if (name) *name = t;
    return selector_level(t, particle_kind);
}

bool vx_is_particles_view(const infer_t* x, vial_node_id n) {
    uint8_t pk;
    return vx_view_level(x, n, &pk, NULL) != VIAL_LEVEL_NONE;
}

const vial_sets_t* vx_view_units(infer_t* x, md_str_t name) {
    if (!x->sys) return NULL;
    for (size_t k = 0; k < MD_COUNTOF(x->predefined) && x->predefined[k].sets; ++k) {
        if (md_str_eq(x->predefined[k].name, name)) return x->predefined[k].sets;
    }
    uint8_t pkind = 0;
    const uint8_t level = selector_level(name, &pkind);
    const vial_select_t q = {.level = level, .particle_kind = pkind};
    const vial_sets_t* s = vial_select(x->sys, &q, x->alloc);
    for (size_t k = 0; k < MD_COUNTOF(x->predefined); ++k) {
        if (!x->predefined[k].sets) { x->predefined[k].name = name; x->predefined[k].sets = s; break; }
    }
    return s;
}

static T view_at(infer_t* x, vial_node_id id, vial_node_id in, const T* a, vial_value_t av, const vial_sets_t* units, md_str_t what);

// Whether a view stands alone, as a value: not the start of a subscript or a call that picks within each set. On one
// set a slice slices that value, chain[0].residue[0:2], since within one set there is nothing else it could mean.
bool vx_view_bare(const infer_t* x, vial_node_id id) {
    const vial_node_id p = x->n->parent[id];
    if (!p || vx_child(x, p, 0) != id) return true;
    const vial_node_kind_t pk = vx_kind(x, p);
    if (pk == VIAL_N_CALL) return false;
    if (pk != VIAL_N_SUBSCRIPT) return true;
    const vial_node_id on = vx_child(x, id, 0), in = vx_child(x, p, 1);
    const T* a = on ? &x->out->types[on - x->out->base] : NULL;
    return a && a->ndim == 0 && in && vx_kind(x, in) == VIAL_N_SLICE;
}

// s.residue alone, on one set: every residue within s, each with the particles it shares with s, which is
// residue(scope=s). On an array of sets the result would be ragged, so that is refused by the caller.
T vx_view_all(infer_t* x, vial_node_id id, vial_node_id on, const T* a) {
    uint8_t pkind = 0;
    md_str_t name = {0};
    const uint8_t level = vx_view_level(x, id, &pkind, &name);
    const vial_value_t av = vx_value_of(x, on);
    T st = vial_type_scalar(VIAL_DT_SET);
    join_stage(&st, a);
    if (st.stage != VIAL_STAGE_STATIC || av.kind != VIAL_VAL_SETS || !x->sys) return selector_unknown(level, st.stage, id);
    const vial_select_t q = {.level = level, .particle_kind = pkind, .scope = av.sets->single ? av.sets : vial_sets_merge(av.sets, x->alloc)};
    const vial_sets_t* s = vial_select(x->sys, &q, x->alloc);
    if (!s) return selector_unknown(level, VIAL_STAGE_STATIC, id);
    const bool empty = s->single ? vial_sets_size(s, 0) == 0 : s->count == 0;
    if (empty) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, on));
        vx_diag_fix(x, "W0601", id, NULL, "`%.*s` holds no %.*s", (int)e.len, e.ptr, (int)name.len, name.ptr);
    }
    T t = vx_sets_type(x, s, id);
    vx_set_value(x, id, sets_value(s));
    return t;
}

T vx_particle_at(infer_t* x, vial_node_id id, vial_node_id obj, vial_node_id in, const T* a) {
    return view_at(x, id, in, a, vx_value_of(x, obj), NULL, MD_STR_LIT("particle"));
}

T vx_view_at(infer_t* x, vial_node_id id, vial_node_id view, vial_node_id in, const T* a) {
    uint8_t pkind;
    md_str_t name = {0};
    const uint8_t level = vx_view_level(x, view, &pkind, &name);
    vial_value_t av = vx_value_of(x, vx_child(x, view, 0));
    if (level == VIAL_LEVEL_PARTICLE) {
        // atom, bead: the particles of the kind, then positions among them
        if (pkind && av.kind == VIAL_VAL_SETS) {
            const vial_sets_t* kind = vx_view_units(x, name);
            av = kind ? sets_value(av.sets->single ? vial_sets_and(av.sets, 0, kind, 0, x->alloc) : vial_sets_boolop_each(av.sets, kind, true, x->alloc))
                      : (vial_value_t){0};
        }
        // s.particle['CA'] and s.particle[['N', 'CA']] are s['CA'] and s[['N', 'CA']] written in full
        if (in && vx_kind(x, in) == VIAL_N_STRING) { vx_set_value(x, view, av); return vx_named(x, id, view, in, a); }
        if (in && vx_kind(x, in) == VIAL_N_LIST && vx_num_children(x, in) && vx_kind(x, vx_child(x, in, 0)) == VIAL_N_STRING) {
            vx_set_value(x, view, av);
            return vx_named_list(x, id, view, in, a);
        }
        return view_at(x, id, in, a, av, NULL, name);
    }
    const vial_sets_t* units = vx_view_units(x, name);
    if (!units && x->sys) return err();
    return view_at(x, id, in, a, av, units ? units : NULL, name);
}

static T view_at(infer_t* x, vial_node_id id, vial_node_id in, const T* a, vial_value_t av, const vial_sets_t* units, md_str_t what) {
    const bool grouping = !md_str_eq(what, MD_STR_LIT("particle")) && !md_str_eq(what, MD_STR_LIT("atom")) && !md_str_eq(what, MD_STR_LIT("bead"));
    const int wl = (int)what.len;
    if (!in || vx_kind(x, in) == VIAL_N_SLICE || vx_kind(x, in) == VIAL_N_TUPLE) {
        if (in) vx_infer(x, in);
        vx_diag_fix(x, "E0319", id, NULL, "`.%.*s[i]` takes a position, or a list of them, the same within each set", wl, what.ptr);
        return err();
    }
    const T ix = vx_infer(x, in);
    if (vial_type_is_error(&ix)) return err();
    const bool list = ix.dtype == VIAL_DT_INT && ix.ndim == 1;
    if (!(ix.dtype == VIAL_DT_INT && ix.ndim == 0) && !list) {
        vx_diag_fix(x, "E0319", in, NULL, "`.%.*s[i]` takes an int, or a list of them, positions within each set; this is %s", wl, what.ptr,
            vx_tstr(&ix));
        return err();
    }
    if (list && a->ndim >= VIAL_MAX_DIMS) {
        vx_diag_fix(x, "E0308", id, NULL, "the result would have 7 dimensions; the limit is 6");
        return err();
    }
    T t = *a;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    if (!grouping) t.flags |= VIAL_TYPE_ONE_PARTICLE;
    if (list) {
        vial_extent_t e = ix.dims[0];
        e.desc = 0;
        vial_type_push_dim(&t, e);
    }
    join_stage(&t, &ix);
    // The positions, when known
    int64_t one;
    int64_t* pos = &one;
    uint32_t m = 1;
    if (list) {
        if (vx_kind(x, in) != VIAL_N_LIST) return t;
        m = vx_num_children(x, in);
        pos = md_alloc(x->alloc, (m ? m : 1) * sizeof(int64_t));
        for (uint32_t j = 0; j < m; ++j) {
            const vial_value_t c = vx_value_of(x, vx_child(x, in, j));
            if (c.kind != VIAL_VAL_INT) return t;
            pos[j] = c.i;
        }
    } else {
        const vial_value_t iv = vx_value_of(x, in);
        if (iv.kind != VIAL_VAL_INT) return t;
        one = iv.i;
    }
    if (av.kind != VIAL_VAL_SETS || t.stage != VIAL_STAGE_STATIC) return t;
    const vial_sets_t* s = av.sets;
    if (grouping) {
        if (!units) return t;
        uint32_t outside = 0;
        vial_sets_t* r = vial_sets_units_at(s, units, pos, NULL, m, &outside, x->alloc);
        r->single = s->single && !list;
        if (outside) {
            t.flags |= VIAL_TYPE_MISSING;
            if (s->single) vx_diag_fix(x, "W0609", id, NULL, "a position is beyond the %.*ss the set holds, so its entry is missing", wl, what.ptr);
            else vx_diag_fix(x, "W0609", id, NULL, "a position is outside %u of %u sets, which hold fewer %.*ss; %s missing", outside, s->count,
                wl, what.ptr, outside == 1 ? "its entry is" : "their entries are");
        }
        const T rt = vx_sets_type(x, r, id);
        t.flags |= rt.flags & VIAL_TYPE_EMPTY;
        vx_set_value(x, id, sets_value(r));
        return t;
    }
    const uint32_t total = s->count * m;
    vial_sets_t* r = md_alloc(x->alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = total;
    r->single = s->single && !list;
    if (!list) { r->level = s->level; r->unit = s->unit; }
    r->offset = md_alloc(x->alloc, (total + 1) * sizeof(uint32_t));
    r->idx = md_alloc(x->alloc, (total ? total : 1) * sizeof(uint32_t));
    uint32_t at = 0, outside = 0, c = 0, worst = 0;
    int64_t worst_pos = 0;
    r->offset[0] = 0;
    for (uint32_t i = 0; i < s->count; ++i) {
        const int64_t n = vial_sets_size(s, i);
        bool short_set = false;
        for (uint32_t j = 0; j < m; ++j) {
            const int64_t q = pos[j] < 0 ? pos[j] + n : pos[j];
            if (q >= 0 && q < n) r->idx[at++] = vial_sets_items(s, i)[q];
            else { short_set = true; if (!worst) worst_pos = pos[j]; worst = 1; }
            r->offset[++c] = at;
        }
        outside += short_set;
    }
    if (outside) {
        t.flags |= VIAL_TYPE_MISSING;
        if (s->single) {
            vx_diag_fix(x, "W0609", id, NULL, "position %lld is outside the set of %u particles, so its point is missing", (long long)worst_pos,
                vial_sets_size(s, 0));
        } else {
            vx_diag_fix(x, "W0609", id, NULL, "position %lld is outside %u of %u sets, which have fewer particles; %s missing", (long long)worst_pos,
                outside, s->count, outside == 1 ? "its entry is" : "their entries are");
        }
    }
    vx_set_value(x, id, sets_value(r));
    return t;
}

T vx_named_list(infer_t* x, vial_node_id id, vial_node_id on, vial_node_id in, const T* a) {
    const uint32_t m = vx_num_children(x, in);
    md_str_t* names = md_alloc(x->alloc, (m ? m : 1) * sizeof(md_str_t));
    for (uint32_t j = 0; j < m; ++j) {
        const vial_node_id c = vx_child(x, in, j);
        const T ct = vx_infer(x, c);
        if (vial_type_is_error(&ct)) return err();
        const vial_value_t v = vx_value_of(x, c);
        if (v.kind != VIAL_VAL_STR || vx_name_is_glob(v.s)) {
            vx_diag_fix(x, "E0319", c, NULL, "a list names one particle each, as constant strings without * or ?");
            return err();
        }
        names[j] = v.s;
    }
    if (!m) { vx_diag_fix(x, "E0319", in, NULL, "a list of names names at least one"); return err(); }
    if (a->ndim >= VIAL_MAX_DIMS) { vx_diag_fix(x, "E0308", id, NULL, "the result would have 7 dimensions; the limit is 6"); return err(); }
    T t = *a;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    t.flags |= VIAL_TYPE_ONE_PARTICLE;
    vial_type_push_dim(&t, vial_extent_static(m));
    if (!x->sys) return t;
    const vial_sets_t** named = md_alloc(x->alloc, m * sizeof(vial_sets_t*));
    for (uint32_t j = 0; j < m; ++j) {
        named[j] = vx_named_particles(x, names[j]);
        if (vial_sets_size(named[j], 0) == 0) vx_diag_fix(x, "W0601", vx_child(x, in, j), NULL, "no particle is named '%.*s'", (int)names[j].len, names[j].ptr);
    }
    const vial_value_t av = vx_value_of(x, on);
    if (av.kind != VIAL_VAL_SETS || t.stage != VIAL_STAGE_STATIC) return t;
    const vial_sets_t* s = av.sets;
    const uint32_t total = s->count * m;
    vial_sets_t* r = md_alloc(x->alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = total;
    r->offset = md_alloc(x->alloc, (total + 1) * sizeof(uint32_t));
    r->idx = md_alloc(x->alloc, (total ? total : 1) * sizeof(uint32_t));
    uint32_t at = 0, c = 0, lacking = 0;
    int first_bad = -1;
    r->offset[0] = 0;
    for (uint32_t i = 0; i < s->count; ++i) {
        bool bad = false;
        for (uint32_t j = 0; j < m; ++j) {
            const vial_sets_t* both = vial_sets_and(s, i, named[j], 0, x->alloc);
            if (vial_sets_size(both, 0) == 1) r->idx[at++] = vial_sets_items(both, 0)[0];
            else { bad = true; if (first_bad < 0) first_bad = (int)j; }
            r->offset[++c] = at;
        }
        lacking += bad;
    }
    if (lacking) {
        t.flags |= VIAL_TYPE_MISSING;
        const md_str_t n0 = names[first_bad];
        if (s->single) {
            vx_diag_fix(x, "W0609", id, NULL, "the set has no particle named '%.*s', or several, so its point is missing", (int)n0.len, n0.ptr);
        } else {
            vx_diag_fix(x, "W0609", id, NULL, "%u of %u sets lack a particle of a name, or have it twice, such as '%.*s'; %s missing",
                lacking, s->count, (int)n0.len, n0.ptr, lacking == 1 ? "its entry is" : "their entries are");
        }
    }
    vx_set_value(x, id, sets_value(r));
    return t;
}

T vx_particles_call(infer_t* x, vial_node_id call, vial_node_id view, const T* a) {
    uint8_t pkind;
    md_str_t vname = MD_STR_LIT("particle");
    vx_view_level(x, view, &pkind, &vname);
    // Its arguments as the selector takes them: residue('HIS') within each set is s.residue('HIS')
    // The selector with these keywords, then each set of a and the particles of what it gives
    const T sel = selector_call(x, call, vname);
    if (vial_type_is_error(&sel)) return err();
    const vial_value_t sv = vx_value_of(x, call);
    vx_set_value(x, call, (vial_value_t){0});
    vx_set_value(x, view, sv);     // what emit reads: the selection, which the keywords make known when compiling
    T t = *a;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    join_stage(&t, &sel);
    const vial_value_t av = vx_value_of(x, vx_child(x, view, 0));
    if (av.kind == VIAL_VAL_SETS && sv.kind == VIAL_VAL_SETS && t.stage == VIAL_STAGE_STATIC) {
        const vial_sets_t* one = sv.sets->single ? sv.sets : vial_sets_merge(sv.sets, x->alloc);
        const vial_value_t r = av.sets->single ? sets_value(vial_sets_and(av.sets, 0, one, 0, x->alloc))
                                               : sets_value(vial_sets_boolop_each(av.sets, one, true, x->alloc));
        const T st = vx_sets_type(x, r.sets, call);
        t.flags |= st.flags & (VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
        vx_set_value(x, call, r);
    }
    return t;
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

vial_value_t vx_sets_subscript(infer_t* x, const vial_value_t* a, vial_node_id item, uint32_t inner) {
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
            return (vial_value_t){.kind = VIAL_VAL_STR, .s = md_alloc_str(x->alloc, (md_str_t){buf, len})};
        }
        return (vial_value_t){.kind = VIAL_VAL_LABELS, .sets = vial_sets_pick(a->sets, pos, cnt, x->alloc)};
    }
    if (a->kind == VIAL_VAL_STRS) {
        uint32_t cnt; bool single;
        const uint32_t* pos = positions(x, item, a->n, &cnt, &single);
        if (!pos) return none;
        if (single) return (vial_value_t){.kind = VIAL_VAL_STR, .s = a->strs[pos[0]]};
        md_str_t* strs = md_alloc(x->alloc, (cnt + 1) * sizeof(md_str_t));
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
    // An array of more dimensions than one: positions along the first, each a block of the rest
    if (inner == 0 || s->count % inner) return none;
    const uint32_t* pos = positions(x, item, s->count / inner, &cnt, &single);
    if (!pos) return none;
    if (inner > 1) {
        uint32_t* all = md_alloc(x->alloc, ((size_t)cnt * inner + 1) * sizeof(uint32_t));
        for (uint32_t k = 0; k < cnt; ++k) for (uint32_t j = 0; j < inner; ++j) all[k * inner + j] = pos[k] * inner + j;
        vial_sets_t* r = vial_sets_pick(s, all, cnt * inner, x->alloc);
        r->single = false;
        return sets_value(r);
    }
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
            if (kw_axis && kn && md_str_eq(vx_text(x, kn), MD_STR_LIT("axis"))) { *kw_axis = vx_child(x, a, 1); continue; }
            const md_str_t kw = kn ? vx_text(x, kn) : MD_STR_LIT("");
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
    const vial_node_id arg = only_positional(x, call, "label", &axis_node);
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
        vx_diag_fix(x, "E0316", axis_node ? axis_node : arg, NULL, "`label` takes a dimension of %s, which has %u", vx_tstr(&a), a.ndim);
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
            const md_str_t kw = kn ? vx_text(x, kn) : MD_STR_LIT("");
            if (md_str_eq(kw, MD_STR_LIT("of"))) { of = vn; continue; }
            if (md_str_eq(kw, MD_STR_LIT("r")) && !r) { r = vn; continue; }
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
        const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, a));
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
    vial_unit_name(MD_STR_LIT("angstrom"), &length);
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

// A bond keyword's value (kind=, order=, origin=) as its expression, as filter_of has a selector's: order= takes ints and
// words alike, 1 to 4 (MD_BOND_ORDER_QUADRUPLE), 'aromatic' or 'delocalized'; kind= and origin= words of mdlib's names. 1 known, 0 not known,
// -1 refused (reported).
static int bond_filter_of(infer_t* x, uint8_t key, vial_node_id node, vial_filter_t* f) {
    const uint8_t want = key == VIAL_KEY_ORDER ? WANT_INT_OR_STR : VIAL_DT_STR;
    f->key = key;
    f->expr = NULL;
    vial_node_id bad = node;
    const int known = expr_of(x, node, want, &f->expr, &bad);
    static const char* orders = "1, 2, 3, 4, 'aromatic' or 'delocalized'";
    if (known < 0) {
        const T bt = x->out->types[(bad) - x->out->base];
        if (key == VIAL_KEY_ORDER) vx_diag_fix(x, "E0206", bad, NULL, "`order` takes %s, a list of them, or and/or/not over them; this is %s", orders, vx_tstr(&bt));
        else vx_diag_fix(x, "E0206", bad, NULL, "`%s` takes a string, a list of strings, or and/or/not over strings; this is %s", vial_select_key_word(key), vx_tstr(&bt));
        return -1;
    }
    for (size_t i = 0; known == 1 && i < md_array_size(f->expr); ++i) {
        const vial_match_t* m = &f->expr[i];
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), m->s);
        if (key == VIAL_KEY_ORDER && ((m->kind == VIAL_MATCH_STR && !vial_select_bond_order_known(m->s)) ||
                                      (m->kind == VIAL_MATCH_INT && (m->lo < 1 || m->lo > 4)))) {
            if (m->kind == VIAL_MATCH_STR) vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is no order of a bond; the orders are %s", (int)e.len, e.ptr, orders);
            else vx_diag_fix(x, "E0206", node, NULL, "%lld is no order of a bond; the orders are %s", (long long)m->lo, orders);
            return -1;
        }
        if (m->kind != VIAL_MATCH_STR) continue;
        if (key == VIAL_KEY_KIND && !vial_select_bond_kind_known(m->s)) {
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is no kind of bond; the kinds are %s", (int)e.len, e.ptr, vial_select_bond_kinds());
            return -1;
        }
        if (key == VIAL_KEY_ORIGIN && !vial_select_bond_origin_known(m->s)) {
            vx_diag_fix(x, "E0206", node, NULL, "'%.*s' is no origin of a bond; the origins are file, topology, user and inferred", (int)e.len, e.ptr);
            return -1;
        }
    }
    return known;
}

// bond(s=None, kind=, order=, origin=): per bond with both ends in s, its two particles as one-particle sets, set[n][2],
// of the kinds, orders and origins given (vial_select_bonds); without kind=, the strong bonds, covalent or coordinate.
// Known when compiling for a set known then; per frame for one that changes per frame.
static T call_bonds(infer_t* x, vial_node_id call) {
    vial_node_id sn = 0;
    vial_select_t q = {0};      // its filters, as a selector's
    bool failed = false, known = true;
    vial_node_id logic_at = 0;
    bool no_kind = false;       // kind= lets no kind through, said (W0601)
    const uint32_t nc = vx_num_children(x, call);
    uint32_t num_pos = 0;
    for (uint32_t i = 1; i < nc; ++i) {
        const vial_node_id arg = vx_child(x, call, i);
        if (!arg) continue;
        if (vx_kind(x, arg) != VIAL_N_KWARG) {
            if (++num_pos > 1) {
                vx_diag_fix(x, "E0203", arg, NULL, "`bond` takes one positional argument, the set");
                failed = true;
            } else {
                sn = arg;
            }
            continue;
        }
        const vial_node_id kn = vx_child(x, arg, 0), vn = vx_child(x, arg, 1);
        if (!kn || !vn) { failed = true; continue; }
        const md_str_t kw = vx_text(x, kn);
        if (md_str_eq(kw, MD_STR_LIT("s"))) { sn = vn; continue; }
        const uint8_t key = vial_select_key(kw);
        if (key != VIAL_KEY_KIND && key != VIAL_KEY_ORDER && key != VIAL_KEY_ORIGIN) {
            vx_diag_fix(x, "E0202", arg, NULL, "`bond` takes no keyword `%.*s`; it takes s, kind, order and origin", (int)kw.len, kw.ptr);
            failed = true;
            continue;
        }
        T vt;
        if (value_logic(x, vn, &vt)) { if (!logic_at) logic_at = vn; }
        else vt = vx_infer(x, vn);
        if (vial_type_is_error(&vt)) { failed = true; continue; }
        vial_filter_t f;
        const int r = bond_filter_of(x, key, vn, &f);
        if (r < 0) { failed = true; continue; }
        if (r == 0) known = false;
        if (r == 1 && key == VIAL_KEY_KIND && !vial_select_bond_kind_possible(&f) && !no_kind) {
            // Of no system: a bond is of one kind, so and between two kinds holds for none
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, vn));
            vx_diag_fix(x, "W0601", vn, NULL, "no bond is of the kind `%.*s`: a bond is of one kind, so and between two kinds holds for none", (int)e.len, e.ptr);
            no_kind = true;
        }
        md_array_push(q.filters, f, x->alloc);
    }
    if (failed) return err();
    T st = sn ? vx_infer(x, sn) : vial_type_scalar(VIAL_DT_NONE);
    if (vial_type_is_error(&st)) return err();
    if (st.dtype != VIAL_DT_SET && st.dtype != VIAL_DT_NONE) {
        vx_diag_fix(x, "E0313", sn, NULL, "`bond` takes a set; this is %s", vx_tstr(&st));
        return err();
    }
    if (st.dtype == VIAL_DT_SET && st.ndim) {
        vx_diag_fix(x, "E0313", sn, NULL, "`bond` takes one set; particle(...) merges an array of them");
        return err();
    }
    T t = vial_type_scalar(VIAL_DT_SET);
    t.flags |= VIAL_TYPE_ONE_PARTICLE;
    t.stage = st.stage;
    const vial_value_t sv = sn ? vx_value_of(x, sn) : (vial_value_t){0};
    const bool all = !sn || st.dtype == VIAL_DT_NONE;
    // Hydrogen bonds are found at each frame from the geometry (md_hbond), where the filters let one through, as of
    // kind='hydrogen' with no origin= or order= to keep them out: the bonds are then per frame
    if (vial_select_bonds_hydrogen(q.filters, md_array_size(q.filters))) {
        if (logic_at) {
            vx_diag_fix(x, "E0801", logic_at, NULL, "and/or/not in a keyword's value is resolved when compiling; hydrogen bonds, found per frame, "
                "do not take it yet: a list takes its place");
            return err();
        }
        if (x->sys) {
            md_hbond_query_t hq = {0};
            md_temp_scope_t temp = md_temp_begin_avoid(x->alloc);
            const bool none = md_hbond_query_init(&hq, &(md_hbond_desc_t){0}, x->sys->sys, md_temp_allocator(temp)) && (hq.flags & MD_HBOND_FLAG_NO_HYDROGENS);
            md_temp_end(temp);
            if (none) {
                char ex[VIAL_EXCERPT_MAX + 4];
                const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, call));
                vx_diag_fix(x, "W0601", call, NULL, "`%.*s` finds no hydrogen bonds: the system's N and O carry no hydrogens", (int)e.len, e.ptr);
            }
        }
        t.stage = VIAL_STAGE_FRAME;
        vial_type_push_dim(&t, (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
        vial_type_push_dim(&t, vial_extent_static(2));
        return t;
    }
    if (!x->sys || st.stage != VIAL_STAGE_STATIC || (!all && sv.kind != VIAL_VAL_SETS) || !known) {
        // and/or/not is resolved when compiling, as for a selector
        if (logic_at && x->sys && st.stage != VIAL_STAGE_STATIC) {
            vx_diag_fix(x, "E0801", logic_at, NULL, "and/or/not in a keyword's value is resolved when compiling; a selection that changes per frame does not take it yet");
            return err();
        }
        vial_type_push_dim(&t, (vial_extent_t){.kind = st.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
        vial_type_push_dim(&t, vial_extent_static(2));
        return t;
    }
    vial_sets_t* r = vial_select_bonds(x->sys, all ? NULL : sv.sets, q.filters, md_array_size(q.filters), x->alloc);
    if (r->count == 0 && md_array_size(q.filters) && !no_kind) {
        // Keywords that match no bond, as a selector's that match nothing; a set without bonds, an ion's, is none
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), vx_text(x, call));
        vx_diag_fix(x, "W0601", call, NULL, "`%.*s` matches nothing in its scope", (int)e.len, e.ptr);
    }
    vial_type_push_dim(&t, vial_extent_static(r->count / 2));
    vial_type_push_dim(&t, vial_extent_static(2));
    vx_set_value(x, call, sets_value(r));
    return t;
}

// formed(s=None, tolerance=0.1, breakable=['file', 'inferred']), broken(...): per frame, the bonds formed in it that the
// system does not have, or the system's broken in it, with an atom in s (md_bond_delta_compute), as bond() gives bonds,
// set[?][2]. A bond breaks beyond (1 + tolerance) times its cutoff and forms within (1 - tolerance) times it, and only
// the bonds of the origins breakable break.
static T call_bond_changes(infer_t* x, vial_node_id call, const char* fname) {
    static const vx_param_t p[] = {{"s", VX_KIND_WHOLE, false}, {"tolerance", VX_KIND_WHOLE, false}, {"breakable", VX_KIND_WHOLE, false}};
    vial_node_id a[3] = {0, 0, 0};
    if (!vx_bind(x, call, fname, p, 3, a)) return err();
    if (a[0]) {
        const T st = x->out->types[(a[0]) - x->out->base];
        if (vial_type_is_error(&st)) return err();
        if (st.dtype != VIAL_DT_SET && st.dtype != VIAL_DT_NONE) {
            vx_diag_fix(x, "E0313", a[0], NULL, "`%s` takes a set, the atoms whose bonds are looked at; this is %s", fname, vx_tstr(&st));
            return err();
        }
        if (st.dtype == VIAL_DT_SET && st.ndim) {
            vx_diag_fix(x, "E0313", a[0], NULL, "`%s` takes one set; particle(...) merges an array of them", fname);
            return err();
        }
    }
    if (a[1]) {
        const vial_value_t v = vx_value_of(x, a[1]);
        const double tol = v.kind == VIAL_VAL_FLOAT ? v.f : v.kind == VIAL_VAL_INT ? (double)v.i : -1.0;
        if (!(tol >= 0.0 && tol < 0.5)) {
            vx_diag_fix(x, "E0313", a[1], NULL, "tolerance= is a fraction of a bond's cutoff, at least 0 and less than 0.5, known when compiling");
            return err();
        }
    }
    if (a[2]) {
        // Origins, one or a list of them
        const T bt = x->out->types[(a[2]) - x->out->base];
        if (vial_type_is_error(&bt)) return err();
        vial_match_t* e = NULL;
        vial_node_id bad = a[2];
        const int known = bt.dtype == VIAL_DT_STR && bt.ndim <= 1 ? expr_of(x, a[2], VIAL_DT_STR, &e, &bad) : -1;
        if (known < 0) {
            vx_diag_fix(x, "E0313", a[2], NULL, "breakable= is an origin or a list of them, file, topology, user and inferred; this is %s", vx_tstr(&bt));
            return err();
        }
        for (size_t i = 0; known == 1 && i < md_array_size(e); ++i) {
            if (e[i].kind != VIAL_MATCH_STR || vial_select_bond_origin_known(e[i].s)) continue;
            char ex[VIAL_EXCERPT_MAX + 4];
            const md_str_t es = vial_excerpt(ex, sizeof(ex), e[i].s);
            vx_diag_fix(x, "E0206", a[2], NULL, "'%.*s' is no origin of a bond; the origins are file, topology, user and inferred", (int)es.len, es.ptr);
            return err();
        }
    }
    T t = vial_type_scalar(VIAL_DT_SET);
    t.flags |= VIAL_TYPE_ONE_PARTICLE;
    t.stage = VIAL_STAGE_FRAME;
    vial_type_push_dim(&t, (vial_extent_t){.kind = VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
    vial_type_push_dim(&t, vial_extent_static(2));
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
        if (c.kind != VIAL_VAL_STR || !(md_str_eq(c.s, MD_STR_LIT("mean")) || md_str_eq(c.s, MD_STR_LIT("max")) || md_str_eq(c.s, MD_STR_LIT("first")))) {
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

// ### PATTERNS ###

typedef struct match_acc_t {
    infer_t*        x;
    const md_match_query_t* q;
    bool            atoms;      // smiles_atoms: one particle per (tagged) pattern atom
    uint32_t        width;      // its columns
    uint32_t        cols[64];   // the query atoms of the columns, in class order or pattern order
    uint32_t*       idx;        // md_array: the particles of each match
    uint32_t*       offset;     // md_array
    md_hashmap32_t  seen;       // smiles_atoms with classes: one match per set of tagged atoms
    bool            overflow;
} match_acc_t;

static int cmp_u32_(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return x < y ? -1 : x > y;
}

static bool on_match(const int32_t* atom_idx, size_t width, uint32_t unit, void* user) {
    (void)unit;
    match_acc_t* m = user;
    if (md_array_size(m->offset) > 1000000) { m->overflow = true; return false; }
    if (m->atoms) {
        uint32_t row[64];
        for (uint32_t c = 0; c < m->width; ++c) {
            const int32_t a = m->cols[c] < width ? atom_idx[m->cols[c]] : -1;
            if (a < 0) return true;     // a pattern hydrogen with no hydrogen atom to map onto: no particle for the column
            row[c] = (uint32_t)a;
        }
        // Unique over the tagged atoms: the OXT and O of a C-terminus give one match
        uint32_t key[64];
        memcpy(key, row, m->width * sizeof(uint32_t));
        qsort(key, m->width, sizeof(uint32_t), cmp_u32_);
        const uint64_t h = md_hash64(key, m->width * sizeof(uint32_t), 0) & ~(1ull << 63);
        if (md_hashmap_get(&m->seen, h)) return true;
        md_hashmap_add(&m->seen, h, 1);
        for (uint32_t c = 0; c < m->width; ++c) {
            md_array_push(m->idx, row[c], m->x->alloc);
            md_array_push(m->offset, (uint32_t)md_array_size(m->idx), m->x->alloc);
        }
        return true;
    }
    const size_t mark = md_array_size(m->idx);
    for (size_t i = 0; i < width; ++i) if (atom_idx[i] >= 0) md_array_push(m->idx, (uint32_t)atom_idx[i], m->x->alloc);
    qsort(m->idx + mark, md_array_size(m->idx) - mark, sizeof(uint32_t), cmp_u32_);
    // The same particle twice in a match (a hydrogen mapped once per bond) is one particle of the set
    uint32_t w = (uint32_t)mark;
    for (size_t i = mark; i < md_array_size(m->idx); ++i) if (w == mark || m->idx[w - 1] != m->idx[i]) m->idx[w++] = m->idx[i];
    md_array_shrink(m->idx, w);
    md_array_push(m->offset, (uint32_t)md_array_size(m->idx), m->x->alloc);
    return true;
}

// smiles(pattern, scope=None, level='molecule', mode='unique') and smiles_atoms(...)
static T call_smiles(infer_t* x, vial_node_id call, bool atoms) {
    const char* fname = atoms ? "smiles_atoms" : "smiles";
    static const vx_param_t p[] = {{"pattern", VX_KIND_WHOLE, true}, {"scope", VX_KIND_WHOLE, false}, {"level", VX_KIND_WHOLE, false},
                                   {"mode", VX_KIND_WHOLE, false}};
    vial_node_id a[4];
    if (!vx_bind(x, call, fname, p, 4, a)) return err();
    for (int k = 0; k < 4; ++k) if (a[k] && vial_type_is_error(&x->out->types[(a[k]) - x->out->base])) return err();
    const vial_value_t pv = vx_value_of(x, a[0]);
    if (pv.kind != VIAL_VAL_STR) {
        vx_diag_fix(x, "E0313", a[0], NULL, "`%s` takes a pattern, a string known when compiling", fname);
        return err();
    }
    md_match_level_t level = MD_MATCH_LEVEL_MOLECULE;
    if (a[2]) {
        const vial_value_t v = vx_value_of(x, a[2]);
        if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("molecule"))) level = MD_MATCH_LEVEL_MOLECULE;
        else if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("residue"))) level = MD_MATCH_LEVEL_RESIDUE;
        else if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("part"))) level = MD_MATCH_LEVEL_PART;
        else { vx_diag_fix(x, "E0206", a[2], NULL, "`level` takes 'molecule', 'residue' or 'part'"); return err(); }
    }
    md_match_mode_t mode = MD_MATCH_MODE_UNIQUE;
    if (a[3]) {
        const vial_value_t v = vx_value_of(x, a[3]);
        if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("unique"))) mode = MD_MATCH_MODE_UNIQUE;
        else if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("all"))) mode = MD_MATCH_MODE_ALL;
        else if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("one_per_unit"))) mode = MD_MATCH_MODE_ONE_PER_UNIT;
        else if (v.kind == VIAL_VAL_STR && md_str_eq(v.s, MD_STR_LIT("disjoint"))) mode = MD_MATCH_MODE_DISJOINT;
        else { vx_diag_fix(x, "E0206", a[3], NULL, "`mode` takes 'unique', 'all', 'one_per_unit' or 'disjoint'"); return err(); }
    }
    T st = a[1] ? x->out->types[(a[1]) - x->out->base] : vial_type_scalar(VIAL_DT_NONE);
    if (st.dtype != VIAL_DT_SET && st.dtype != VIAL_DT_NONE) {
        vx_diag_fix(x, "E0206", a[1], NULL, "`scope` takes a set or system; this is %s", vx_tstr(&st));
        return err();
    }
    // The pattern, read when compiling
    md_match_query_t q = {0};
    md_smiles_error_t perr = {0};
    if (!md_match_query_init_smiles(&q, pv.s, x->alloc, &perr)) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), pv.s);
        vx_diag_fix(x, "E0206", a[0], NULL, "`%.*s` is not a SMILES pattern: %s at character %zu", (int)e.len, e.ptr, perr.message, perr.offset);
        return err();
    }
    // smiles_atoms: the tagged atoms in class order, or every atom but the pattern's hydrogens
    match_acc_t m = {.x = x, .q = &q, .atoms = atoms};
    m.seen.allocator = x->alloc;
    if (atoms) {
        bool tagged = false;
        for (size_t i = 0; i < q.num_atoms; ++i) tagged = tagged || q.atoms[i].tag;
        for (uint32_t tag = 1, done = 0; tagged && done < 64 && tag < UINT16_MAX; ++tag) {
            bool any = false;
            for (size_t i = 0; i < q.num_atoms && m.width < 64; ++i) if (q.atoms[i].tag == tag) { m.cols[m.width++] = (uint32_t)i; any = true; }
            if (any) ++done;
            if (m.width == 64) break;
            // Classes are sparse: stop once every tagged atom has its column
            uint32_t tagged_n = 0;
            for (size_t i = 0; i < q.num_atoms; ++i) tagged_n += q.atoms[i].tag != 0;
            if (m.width >= tagged_n) break;
        }
        for (size_t i = 0; !tagged && i < q.num_atoms && m.width < 64; ++i) {
            if (!((q.atoms[i].flags & MD_MATCH_ATOM_ELEMENT) && q.atoms[i].z == MD_Z_H)) m.cols[m.width++] = (uint32_t)i;
        }
    }
    T t = vial_type_scalar(VIAL_DT_SET);
    t.flags |= VIAL_TYPE_OVERLAP;
    const vial_value_t sv = a[1] ? vx_value_of(x, a[1]) : (vial_value_t){0};
    const bool scope_known = !a[1] || st.dtype == VIAL_DT_NONE || (sv.kind == VIAL_VAL_SETS && st.stage == VIAL_STAGE_STATIC);
    if (!x->sys || !scope_known) {
        vial_type_push_dim(&t, (vial_extent_t){.kind = st.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)call});
        if (atoms) { vial_type_push_dim(&t, vial_extent_static(m.width)); t.flags |= VIAL_TYPE_ONE_PARTICLE; }
        t.stage = st.stage;
        md_match_query_free(&q);
        return t;
    }
    md_bitset_t mask = {0};
    md_match_desc_t desc = {.query = &q, .level = level, .mode = mode, .max_steps = 50000000};
    if (a[1] && st.dtype == VIAL_DT_SET) {
        const vial_sets_t* s = sv.sets->single ? sv.sets : vial_sets_merge(sv.sets, x->alloc);
        mask = md_bitset_from_indices(vial_sets_items(s, 0), vial_sets_size(s, 0), x->alloc);
        desc.mask = &mask;
    }
    md_array_push(m.offset, 0, x->alloc);
    bool incomplete = false;
    const bool ok = md_match_for_each(NULL, &incomplete, &desc, x->sys->sys, on_match, &m, x->alloc);
    md_match_query_free(&q);
    if (!ok) {
        vx_diag_fix(x, "E0206", a[0], NULL, "`%s` cannot search for this pattern", fname);
        return err();
    }
    vial_sets_t* r = md_alloc(x->alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = (uint32_t)md_array_size(m.offset) - 1;
    r->offset = m.offset;
    r->idx = m.idx ? m.idx : md_alloc(x->alloc, sizeof(uint32_t));
    const uint32_t rows = atoms ? (m.width ? r->count / m.width : 0) : r->count;
    if (r->count == 0) {
        char ex[VIAL_EXCERPT_MAX + 4];
        const md_str_t e = vial_excerpt(ex, sizeof(ex), pv.s);
        vx_diag_fix(x, "W0601", call, NULL, "`%.*s` matches nothing in its scope", (int)e.len, e.ptr);
    }
    vial_type_push_dim(&t, vial_extent_static(rows));
    if (atoms) {
        vial_type_push_dim(&t, vial_extent_static(m.width));
        t.flags |= VIAL_TYPE_ONE_PARTICLE;
    } else {
        // Matches may share particles; say so only when they do
        const vial_sets_t* merged = vial_sets_merge(r, x->alloc);
        if (vial_sets_size(merged, 0) == r->offset[r->count]) t.flags &= ~VIAL_TYPE_OVERLAP;
    }
    vx_set_value(x, call, sets_value(r));
    return t;
}

// ### ATTRIBUTES ###

// The type of one attribute of the run: its numbers, its unit, its shape less the frame axis, and its components
static bool attribute_type(const md_attribute_t* at, T* t) {
    const md_attribute_format_t* f = &at->format;
    if (!md_attribute_type_is_numeric(f->type)) return false;
    *t = vial_type_scalar(f->type == MD_ATTRIBUTE_TYPE_F32 || f->type == MD_ATTRIBUTE_TYPE_F64 ? VIAL_DT_FLOAT : VIAL_DT_INT);
    t->unit = at->unit;
    const bool temporal = at->flags & MD_ATTRIBUTE_FLAG_TEMPORAL;
    t->stage = temporal ? VIAL_STAGE_FRAME : VIAL_STAGE_STATIC;
    for (uint32_t k = temporal ? 1 : 0; k < f->rank; ++k) if (!vial_type_push_dim(t, vial_extent_static(f->shape[k]))) return false;
    if (f->components > 1 && !vial_type_push_dim(t, vial_extent_static(f->components))) return false;
    return true;
}

// Whether a per-frame attribute has a value at every frame of the run: on the run's own axis, or on one of its own (an
// energy file's, sampled at its own rate) that every frame of the run lands on, as evaluation finds its row
// (md_attribute_axis_map). One that does not may be missing at a frame.
static bool every_frame_sampled(infer_t* x, const md_attribute_t* at) {
    if (!(at->flags & MD_ATTRIBUTE_FLAG_TEMPORAL)) return true;
    const md_attributes_t* attrs = &x->sys->sys->attributes;
    const md_str_t run = x->c->run;
    char p[512];
    snprintf(p, sizeof(p), "%.*s/time", (int)run.len, run.ptr);
    const md_attribute_t* frames = md_attributes_find(attrs, md_str_from_cstr(p));
    if (!frames || frames->format.rank < 1) return true;    // no frames to miss
    const size_t num = frames->format.shape[0];
    const md_attribute_t* axis = md_attributes_axis(attrs, at);
    if (!axis) return at->format.rank >= 1 && at->format.shape[0] >= num;     // by the frame's index
    if (axis->id == frames->id) return true;
    for (size_t i = 0; i < num; ++i) {
        size_t row;
        if (!md_attribute_axis_map(&row, attrs, frames, i, attrs, axis)) return false;
    }
    return true;
}

// attribute_type, and optional where a frame may have no value
static bool run_attribute_type(infer_t* x, const md_attribute_t* at, T* t) {
    if (!attribute_type(at, t)) return false;
    if (!every_frame_sampled(x, at)) t->flags |= VIAL_TYPE_MISSING;
    return true;
}

// A per-atom attribute of the system: under atom/, numeric, not per frame, and one value (of its components, at each
// variant) per atom, its atom axis last
static bool atom_attribute(const md_attribute_t* at, uint32_t num_atoms) {
    const md_attribute_format_t* f = &at->format;
    return md_attribute_type_is_numeric(f->type) && !(at->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) && f->rank >= 1 && f->rank <= 2 &&
           f->shape[f->rank - 1] == num_atoms;
}

// A per-atom attribute of the run: under <run>/atom/, numeric, per frame, and one value per atom at each frame (of its
// components, at each variant), its atom axis last
static bool run_atom_attribute(const md_attribute_t* at, uint32_t num_atoms) {
    const md_attribute_format_t* f = &at->format;
    return md_attribute_type_is_numeric(f->type) && (at->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) && f->rank >= 2 && f->rank <= 3 &&
           f->shape[f->rank - 1] == num_atoms;
}

// attr(name, of=s): one value per particle of s of the per-atom attribute <name>, as the tooltip of a particle names it:
// the run's per frame (<run>/atom/<name>) where it has one, else the system's (atom/<name>). Its variants and
// components are dimensions after the particles'. The system's is missing where it holds NAN.
static T attr_of(infer_t* x, vial_node_id at, md_str_t name, vial_node_id on) {
    const T s = x->out->types[(on) - x->out->base];
    if (vial_type_is_error(&s)) return err();
    char ex[VIAL_EXCERPT_MAX + 4];
    const md_str_t e = vial_excerpt(ex, sizeof(ex), name);
    if (s.dtype != VIAL_DT_SET) {
        vx_diag_fix(x, "E0313", on, NULL, "`attr` of= takes a set; this is %s", vx_tstr(&s));
        return err();
    }
    if (s.ndim) {
        char so[VIAL_EXCERPT_MAX + 4], fix[VIAL_EXCERPT_MAX + 16];
        const md_str_t o = vial_excerpt(so, sizeof(so), vx_text(x, on));
        snprintf(fix, sizeof(fix), "particle(%.*s)", (int)o.len, o.ptr);
        vx_diag_fix(x, "E0318", on, fix, "`%.*s` holds sets whose sizes may differ; their particles' values do not form one array",
            (int)o.len, o.ptr);
        return err();
    }
    if (vial_is_glob(name)) {
        vx_diag_fix(x, "E0206", at, NULL, "`attr` of= takes the name of one attribute; `%.*s` is a pattern", (int)e.len, e.ptr);
        return err();
    }
    if (name.len > 5 && strncmp(name.ptr, "atom/", 5) == 0) {
        char fix[VIAL_EXCERPT_MAX + 8];
        snprintf(fix, sizeof(fix), "'%.*s'", (int)(name.len - 5), name.ptr + 5);
        vx_diag_fix(x, "E0206", at, fix, "`attr` of= takes the name below atom/, as a particle's tooltip writes it");
        return err();
    }
    if (!x->sys) return err();      // nothing to read: typed when there is a system
    const md_attributes_t* attrs = &x->sys->sys->attributes;
    const uint32_t num_atoms = x->sys->num_atoms;
    const md_str_t run = x->c->run;
    char full[512];
    // The run's per-frame values of the particles, where it has them; else the system's
    const md_attribute_t* a = NULL;
    if (!md_str_empty(run)) {
        snprintf(full, sizeof(full), "%.*s/atom/%.*s", (int)run.len, run.ptr, (int)name.len, name.ptr);
        a = md_attributes_find(attrs, md_str_from_cstr(full));
        if (a && !run_atom_attribute(a, num_atoms)) a = NULL;
    }
    if (!a) {
        snprintf(full, sizeof(full), "atom/%.*s", (int)name.len, name.ptr);
        a = md_attributes_find(attrs, md_str_from_cstr(full));
        if (a && !atom_attribute(a, num_atoms)) a = NULL;
    }
    if (!a) {
        // The names there are, the run's and the system's, as the message's list
        char have[256] = "";
        size_t len = 0, num = 0;
        char group[512];
        snprintf(group, sizeof(group), "%.*s/atom", (int)run.len, run.ptr);
        for (int g = md_str_empty(run) ? 1 : 0; g < 2; ++g) {
            const md_str_t gs = g == 0 ? md_str_from_cstr(group) : MD_STR_LIT("atom");
            for (md_attribute_iter_t it = md_attributes_iter(attrs, gs); md_attributes_next(&it) && len < sizeof(have) - 8;) {
                if (g == 0 ? !run_atom_attribute(it.attr, num_atoms) : !atom_attribute(it.attr, num_atoms)) continue;
                const md_str_t rel = {it.attr->path.ptr + gs.len + 1, it.attr->path.len - gs.len - 1};
                const int w = snprintf(have + len, sizeof(have) - len, "%s%.*s", num ? ", " : "", (int)rel.len, rel.ptr);
                if (w < 0 || len + (size_t)w >= sizeof(have) - 8) { snprintf(have + len, sizeof(have) - len, ", ..."); len = sizeof(have); break; }
                len += (size_t)w;
                num += 1;
            }
        }
        if (num) vx_diag_fix(x, "E0206", at, NULL, "particles have no attribute `%.*s`; they have %s", (int)e.len, e.ptr, have);
        else vx_diag_fix(x, "E0206", at, NULL, "particles have no attribute `%.*s`; the data gives them none", (int)e.len, e.ptr);
        return err();
    }
    const md_attribute_format_t* f = &a->format;
    const bool temporal = a->flags & MD_ATTRIBUTE_FLAG_TEMPORAL;
    const bool real = f->type == MD_ATTRIBUTE_TYPE_F32 || f->type == MD_ATTRIBUTE_TYPE_F64;
    T t = vial_type_scalar(real ? VIAL_DT_FLOAT : VIAL_DT_INT);
    t.unit = a->unit;
    t.stage = temporal && s.stage < VIAL_STAGE_FRAME ? VIAL_STAGE_FRAME : s.stage;
    vial_type_push_dim(&t, vx_particle_extent(x, &s, on));
    if (f->rank == (temporal ? 3 : 2)) vial_type_push_dim(&t, vial_extent_static(f->shape[temporal ? 1 : 0]));
    if (f->components > 1) vial_type_push_dim(&t, vial_extent_static(f->components));
    if (temporal) {
        // Per frame: missing at a frame its axis has no sample at; a NaN in it is a value, as from arithmetic
        if (!every_frame_sampled(x, a)) t.flags |= VIAL_TYPE_MISSING;
    } else if (real) {
        // NAN marks an atom the system's attribute has no value for
        const size_t cnt = md_attribute_element_count(f);
        md_temp_scope_t temp = md_temp_begin_avoid(x->alloc);
        float* v = md_temp_alloc_array(temp, float, cnt);
        bool gaps = !v || md_attribute_extract_f32(v, cnt, attrs, a, md_attribute_slice_all(), md_unit_none()) != cnt;
        for (size_t i = 0; v && !gaps && i < cnt; ++i) gaps = isnan(v[i]);
        md_temp_end(temp);
        if (gaps) t.flags |= VIAL_TYPE_MISSING;
    }
    return t;
}

// attr(path): an attribute of the run, such as an energy term; a path with * gives one member per match. attr(name,
// of=s): a per-atom attribute of the system, one value per particle of s (attr_of)
static T call_attr(infer_t* x, vial_node_id call) {
    static const vx_param_t p[] = {{"path", VX_KIND_WHOLE, true}, {"of", VX_KIND_WHOLE, false}};
    vial_node_id args[2] = {0, 0};
    if (!vx_bind(x, call, "attr", p, 2, args)) return err();
    const vial_node_id a = args[0];
    const T pt = x->out->types[(a) - x->out->base];
    if (vial_type_is_error(&pt)) return err();
    const vial_value_t pv = vx_value_of(x, a);
    if (pv.kind != VIAL_VAL_STR) {
        vx_diag_fix(x, "E0313", a, NULL, "`attr` takes a path, a string known when compiling");
        return err();
    }
    if (args[1]) return attr_of(x, a, pv.s, args[1]);
    const md_str_t run = x->c->run;
    if (!x->sys || md_str_empty(run)) return err();    // no run to read: typed when there is one
    const md_attributes_t* attrs = &x->sys->sys->attributes;
    char full[512];
    snprintf(full, sizeof(full), "%.*s/%.*s", (int)run.len, run.ptr, (int)pv.s.len, pv.s.ptr);
    const bool glob = vial_is_glob(pv.s);
    char ex[VIAL_EXCERPT_MAX + 4];
    const md_str_t e = vial_excerpt(ex, sizeof(ex), pv.s);
    if (!glob) {
        const md_attribute_t* at = md_attributes_find(attrs, md_str_from_cstr(full));
        T t;
        if (!at || !run_attribute_type(x, at, &t)) {
            vx_diag_fix(x, "E0206", a, NULL, "the run has no numeric attribute `%.*s`", (int)e.len, e.ptr);
            return err();
        }
        return t;
    }
    // One member per match, which must agree in type, unit and shape
    char prefix[512];
    snprintf(prefix, sizeof(prefix), "%.*s/", (int)run.len, run.ptr);
    md_attribute_iter_t it = md_attributes_iter(attrs, run);
    uint32_t n = 0;
    T first = err();
    while (md_attributes_next(&it)) {
        const md_str_t path = it.attr->path;
        if (path.len <= strlen(prefix) || strncmp(path.ptr, prefix, strlen(prefix)) != 0) continue;
        const md_str_t rel = {path.ptr + strlen(prefix), path.len - strlen(prefix)};
        if (!vial_glob(pv.s, rel)) continue;
        T t;
        if (!run_attribute_type(x, it.attr, &t)) continue;
        if (n == 0) first = t;
        else first.flags |= t.flags & VIAL_TYPE_MISSING;
        if (n > 0 && (t.dtype != first.dtype || t.ndim != first.ndim || !md_unit_equal(t.unit, first.unit) || t.stage != first.stage)) {
            vx_diag_fix(x, "E0206", a, NULL, "the attributes `%.*s` matches differ in type, unit or shape", (int)e.len, e.ptr);
            return err();
        }
        n += 1;
    }
    if (n == 0) {
        vx_diag_fix(x, "E0206", a, NULL, "the run has no numeric attribute that `%.*s` matches", (int)e.len, e.ptr);
        return err();
    }
    T t = first;
    t.ndim = 0;
    vial_type_push_dim(&t, vial_extent_static(n));
    for (uint8_t k = 0; k < first.ndim; ++k) vial_type_push_dim(&t, first.dims[k]);
    return t;
}

bool vx_call_sets(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const md_str_t name = md_str_from_cstr(b->name);
    if ((b->flags & VIAL_BUILTIN_SELECTOR)) { *out = selector_call(x, call, name); return true; }
    if (md_str_eq(name, MD_STR_LIT("len")))    { *out = call_len(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("count")))  { *out = call_count(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("label"))) { *out = call_labels(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("within"))) { *out = call_within(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("bond")))  { *out = call_bonds(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("formed")) || md_str_eq(name, MD_STR_LIT("broken"))) { *out = call_bond_changes(x, call, b->name); return true; }
    if (md_str_eq(name, MD_STR_LIT("spread"))) { *out = call_spread(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("smiles")))       { *out = call_smiles(x, call, false); return true; }
    if (md_str_eq(name, MD_STR_LIT("attr")))         { *out = call_attr(x, call); return true; }
    if (md_str_eq(name, MD_STR_LIT("smiles_atoms"))) { *out = call_smiles(x, call, true); return true; }
    return false;
}

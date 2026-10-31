// The checker, step 3: literals, operators, broadcasting, units and indexing.
//
// Types are inferred bottom up and memoized by node. Every operation follows NumPy for dtypes and shapes, with the
// reference's additions: a bare number takes the unit it meets, `+`, `-` and comparisons convert compatible units and
// make a mismatch unknown (W0501), a quantity literal carries its unit, and two different dynamic extents never
// broadcast. Static values are folded where the operands are known, so `2 < n <= 500` is True when n is.

#include "infer_internal.h"
#include "ir.h"

#include <core/md_common.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef vial_type_t T;

// ### HELPERS ###

static inline vial_node_kind_t kind_of(const infer_t* x, vial_node_id id) { return (vial_node_kind_t)x->t->nodes[id].kind; }
static inline vial_node_id child(const infer_t* x, vial_node_id id, uint32_t i) { return vial_child(x->t, id, i); }
static inline uint32_t num_children(const infer_t* x, vial_node_id id) { return x->t->nodes[id].child_count; }
static inline uint16_t sub_of(const infer_t* x, vial_node_id id) { return x->t->nodes[id].sub; }

// In the prelude, a parameter's name stands for what the script wrote as its argument
static vial_node_id provenance(const infer_t* x, vial_node_id id) {
    if (!x->frame || !x->frame->prelude || !id || x->t->nodes[id].kind != VIAL_N_NAME) return 0;
    const uint32_t sym = x->n->ref[id];
    return sym ? vx_frame_prov(x, sym) : 0;
}

static inline str_t text(const infer_t* x, vial_node_id id) {
    const vial_node_id p = provenance(x, id);
    if (p) {
        const md_vial_span_t s = x->c->tree->nodes[p].span;
        return (str_t){x->c->src.ptr + s.beg, s.end - s.beg};
    }
    const md_vial_span_t s = x->t->nodes[id].span;
    return (str_t){x->src.ptr + s.beg, s.end - s.beg};
}

// A diagnostic, once. In the prelude it is reported in the script: at the argument a parameter stands for, or else at
// the call that led there, without a fix, since a fix would quote the prelude.
static void diag_fix(infer_t* x, const char code[6], vial_node_id at, const char* fix, const char* fmt, ...) {
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    md_vial_span_t span = x->t->nodes[at].span;
    if (x->frame && x->frame->prelude) {
        const vial_node_id p = provenance(x, at);
        if (p) span = x->c->tree->nodes[p].span;
        else { span = x->c->tree->nodes[x->frame->site].span; fix = NULL; }
    }
    if (vial_diag_at(x->diags, span.beg)) {
        const size_t n = md_array_size(x->diags->items);
        for (size_t i = n; i-- > 0;) {
            const md_vial_diagnostic_t* d = &x->diags->items[i];
            if (d->span.beg == span.beg && d->span.end == span.end && strncmp(d->code, code, 5) == 0) return;   // one mistake, one diagnostic
        }
    }
    vial_diag(x->diags, code, span, "%s", buf);
    if (fix) vial_diag_set_fix(x->diags, md_array_size(x->diags->items) - 1, str_from_cstr(fix));
}
#define diag(x, code, at, ...) diag_fix(x, code, at, NULL, __VA_ARGS__)

// A type as text for messages; four buffers, so one message can show several
static const char* tstr(const T* t) {
    static char bufs[4][128];
    static int k = 0;
    char* b = bufs[k++ & 3];
    vial_type_print(b, 128, t);
    return b;
}

static const char* ustr(const T* t) {
    static char bufs[4][64];
    static int k = 0;
    char* b = bufs[k++ & 3];
    if (t->unit_kind != VIAL_UNIT_KNOWN || md_unit_is_none(t->unit)) { snprintf(b, 64, "no unit"); return b; }
    vial_unit_print(b, 64, t->unit);
    return b;
}

void vx_diag_fix(infer_t* x, const char code[6], vial_node_id at, const char* fix, const char* fmt, ...) {
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    diag_fix(x, code, at, fix, "%s", buf);
}
const char* vx_tstr(const vial_type_t* t) { return tstr(t); }
str_t vx_text(const infer_t* x, vial_node_id id) { return text(x, id); }

uint32_t vx_desc(infer_t* x, uint8_t level, const vial_sets_t* sets) {
    if (level == VIAL_LEVEL_NONE) return 0;
    if (!md_array_size((*x->descs))) md_array_push((*x->descs), ((vial_desc_t){0}), x->alloc);
    // The same sets describe the same units: one description, so their dimensions are known to pair
    if (sets) {
        const size_t n = md_array_size(*x->descs);
        for (size_t i = n; i-- > 1 && i + 64 >= n;) {
            if ((*x->descs)[i].sets == sets && (*x->descs)[i].level == level) return (uint32_t)i;
        }
    }
    md_array_push((*x->descs), ((vial_desc_t){.level = level, .sets = sets}), x->alloc);
    return (uint32_t)md_array_size((*x->descs)) - 1;
}

static T scalar(vial_dtype_t d) { return vial_type_scalar(d); }
static T number(vial_dtype_t d) { T t = scalar(d); t.unit_kind = VIAL_UNIT_OPEN; return t; }
static T err(void) { return vial_type_error(); }

static bool is_number(const T* t) {
    return t->dtype == VIAL_DT_INT || t->dtype == VIAL_DT_FLOAT || t->dtype == VIAL_DT_FLOAT64 || t->dtype == VIAL_DT_BOOL;
}
static bool is_float(uint8_t d) { return d == VIAL_DT_FLOAT || d == VIAL_DT_FLOAT64; }

static const char* dtype_word(const T* t) {
    switch (t->dtype) {
    case VIAL_DT_BOOL: return "bool";   case VIAL_DT_INT: return "int";   case VIAL_DT_FLOAT: return "float";
    case VIAL_DT_FLOAT64: return "float64"; case VIAL_DT_STR: return "str"; case VIAL_DT_SET: return "set";
    case VIAL_DT_RANGE: return "range"; case VIAL_DT_DRAWING: return "drawing"; case VIAL_DT_NONE: return "None";
    case VIAL_DT_TUPLE: return "tuple"; case VIAL_DT_DEF: return "def";  default: return "?";
    }
}

// The later stage and pass of two values, and missing if either may be
// The stage and pass of a value computed from a and b. A value that exists per frame (frame stage, or a trajectory
// value that keeps the frame) combined with a trajectory value complete after pass p is computed per frame in pass
// p + 1: `ends - mean(ends, axis='frame')` is trajectory, per frame, pass 2.
void vx_join_stage(T* out, const T* a, const T* b) {
    const bool pa = a->stage == VIAL_STAGE_FRAME || (a->flags & VIAL_TYPE_PER_FRAME);
    const bool pb = b->stage == VIAL_STAGE_FRAME || (b->flags & VIAL_TYPE_PER_FRAME);
    const bool ta = a->stage == VIAL_STAGE_TRAJECTORY && !(a->flags & VIAL_TYPE_PER_FRAME);
    const bool tb = b->stage == VIAL_STAGE_TRAJECTORY && !(b->flags & VIAL_TYPE_PER_FRAME);
    const uint8_t qa = vx_pass(a), qb = vx_pass(b);
    out->stage = a->stage > b->stage ? a->stage : b->stage;
    out->pass  = qa > qb ? qa : qb;
    out->flags = (uint16_t)((out->flags & ~VIAL_TYPE_PER_FRAME) | ((a->flags | b->flags) & (VIAL_TYPE_MISSING | VIAL_TYPE_PER_FRAME)));
    if ((pa && tb) || (pb && ta)) {
        // A per-frame value that reads a finished whole-trajectory result of pass k is computed in pass k + 1
        const uint8_t after = (uint8_t)((ta ? qa : qb) + 1);
        out->stage = VIAL_STAGE_TRAJECTORY;
        out->flags |= VIAL_TYPE_PER_FRAME;
        if (after > out->pass) out->pass = after;
    }
}
static void join_stage(T* out, const T* a, const T* b) { vx_join_stage(out, a, b); }

// ### SHAPES ###

static bool extent_join(vial_extent_t a, vial_extent_t b, vial_extent_t* out) {
    if (a.kind == VIAL_EXT_STATIC && a.n == 1) { *out = b; return true; }
    if (b.kind == VIAL_EXT_STATIC && b.n == 1) { *out = a; return true; }
    if (a.kind == VIAL_EXT_UNKNOWN) { *out = b; return true; }
    if (b.kind == VIAL_EXT_UNKNOWN) { *out = a; return true; }
    if (a.kind != b.kind) return false;
    bool same = false;
    switch (a.kind) {
    case VIAL_EXT_STATIC:  same = a.n == b.n; break;
    case VIAL_EXT_DYNAMIC: same = a.id == b.id && a.offset == b.offset; break;
    case VIAL_EXT_RANGED:  same = a.n == b.n && a.m == b.m; break;
    default: break;
    }
    if (!same) return false;
    *out = a;
    if (!out->desc) out->desc = b.desc;
    return true;
}

// NumPy's broadcasting: shapes align from the right
static bool broadcast(const T* a, const T* b, T* out) {
    const uint8_t n = a->ndim > b->ndim ? a->ndim : b->ndim;
    vial_extent_t dims[VIAL_MAX_DIMS];
    for (uint8_t k = 0; k < n; ++k) {
        const int ia = (int)a->ndim - (int)n + k, ib = (int)b->ndim - (int)n + k;
        if (ia < 0) { dims[k] = b->dims[ib]; continue; }
        if (ib < 0) { dims[k] = a->dims[ia]; continue; }
        if (!extent_join(a->dims[ia], b->dims[ib], &dims[k])) return false;
    }
    out->ndim = n;
    memcpy(out->dims, dims, n * sizeof(dims[0]));
    return true;
}

// The fix for shapes that would broadcast with a new dimension on the right of the shorter: d[:, None]
static void fix_new_axis(infer_t* x, char* fix, size_t cap, vial_node_id short_node, const T* s, const T* l) {
    fix[0] = '\0';
    if (s->ndim == 0 || s->ndim >= l->ndim) return;
    for (uint8_t k = 0; k < s->ndim; ++k) {
        vial_extent_t e;
        if (!extent_join(s->dims[k], l->dims[k], &e)) return;
    }
    const str_t txt = text(x, short_node);
    int n = snprintf(fix, cap, "%.*s[", (int)txt.len, txt.ptr);
    for (uint8_t k = 0; k < l->ndim && n > 0 && (size_t)n < cap; ++k) {
        n += snprintf(fix + n, cap - n, "%s%s", k ? ", " : "", k < s->ndim ? ":" : "None");
    }
    if (n > 0 && (size_t)n < cap) snprintf(fix + n, cap - n, "]");
}

void vx_check_pairing(infer_t* x, vial_node_id at, const T* a, const T* b);

static bool broadcast_or_report(infer_t* x, vial_node_id at, vial_node_id an, const T* a, vial_node_id bn, const T* b, T* out) {
    if (broadcast(a, b, out)) { vx_check_pairing(x, at, a, b); return true; }
    char fix[160] = "";
    if (a->ndim < b->ndim) fix_new_axis(x, fix, sizeof(fix), an, a, b);
    else if (b->ndim < a->ndim) fix_new_axis(x, fix, sizeof(fix), bn, b, a);
    diag_fix(x, "E0301", at, fix[0] ? fix : NULL, "`%s` and `%s` do not broadcast", tstr(a), tstr(b));
    return false;
}

bool vx_broadcast(const T* a, const T* b, T* out) { return broadcast(a, b, out); }

// The unit at position k of what a description runs over: a particle for a set's particles, else the set's unit
static bool desc_unit(const vial_desc_t* d, uint32_t k, uint32_t* out) {
    const vial_sets_t* s = d->sets;
    if (!s) return false;
    if (s->single) { if (k >= vial_sets_size(s, 0)) return false; *out = vial_sets_items(s, 0)[k]; return true; }
    if (!s->unit || k >= s->count) return false;
    *out = s->unit[k];
    return true;
}

// W0605: dimensions of equal extent that run over different particles or units, paired by position. Slices of one
// array or set, and what common() paired, pair without it; so do two that describe the same units in the same order.
void vx_check_pairing(infer_t* x, vial_node_id at, const T* a, const T* b) {
    const uint8_t n = a->ndim < b->ndim ? a->ndim : b->ndim;
    for (uint8_t k = 1; k <= n; ++k) {
        const vial_extent_t ea = a->dims[a->ndim - k], eb = b->dims[b->ndim - k];
        if (ea.kind != VIAL_EXT_STATIC || eb.kind != VIAL_EXT_STATIC || ea.n != eb.n || ea.n < 2) continue;
        if (!ea.desc || !eb.desc || ea.desc == eb.desc) continue;
        if (ea.desc >= md_array_size(*x->descs) || eb.desc >= md_array_size(*x->descs)) continue;
        const vial_desc_t* da = &(*x->descs)[ea.desc];
        const vial_desc_t* db = &(*x->descs)[eb.desc];
        if (!da->sets || !db->sets || vial_sets_root(da->sets) == vial_sets_root(db->sets)) continue;
        for (uint32_t q = 0; q < ea.n; ++q) {
            uint32_t ua, ub;
            if (!desc_unit(da, q, &ua) || !desc_unit(db, q, &ub)) break;
            if (ua == ub && da->level == db->level) continue;
            char la[64] = "", lb[64] = "";
            if (x->sys) {
                vial_unit_label(la, sizeof(la), x->sys, da->sets->single ? VIAL_LEVEL_PARTICLE : da->level, ua);
                vial_unit_label(lb, sizeof(lb), x->sys, db->sets->single ? VIAL_LEVEL_PARTICLE : db->level, ub);
            }
            diag_fix(x, "W0605", at, NULL, "dimensions of %u that run over different units are paired by position: %u is %s against %s; "
                "common(a, b, by='id') pairs them by id", ea.n, q, la[0] ? la : "one", lb[0] ? lb : "another");
            return;
        }
    }
}
bool vx_extent_join(vial_extent_t a, vial_extent_t b, vial_extent_t* out) { return extent_join(a, b, out); }

// ### UNITS ###

// The unit of +, -, comparisons and the like: a bare number takes the other's; compatible units convert to the left's;
// others give W0501 once, and an unknown unit
static void unit_join(infer_t* x, vial_node_id at, const char* verb, const T* a, const T* b, T* out) {
    if (a->unit_kind == VIAL_UNIT_UNKNOWN || b->unit_kind == VIAL_UNIT_UNKNOWN) { out->unit_kind = VIAL_UNIT_UNKNOWN; return; }
    if (a->unit_kind == VIAL_UNIT_OPEN && b->unit_kind == VIAL_UNIT_OPEN) { out->unit_kind = VIAL_UNIT_OPEN; return; }
    if (a->unit_kind == VIAL_UNIT_OPEN) { out->unit_kind = VIAL_UNIT_KNOWN; out->unit = b->unit; return; }
    if (b->unit_kind == VIAL_UNIT_OPEN) { out->unit_kind = VIAL_UNIT_KNOWN; out->unit = a->unit; return; }
    if (md_unit_base_equal(a->unit, b->unit)) { out->unit_kind = VIAL_UNIT_KNOWN; out->unit = a->unit; return; }
    diag_fix(x, "W0501", at, "unitless(...)", "%s %s %s %s; the result's unit becomes unknown", verb, ustr(b),
        strcmp(verb, "compares") == 0 ? "with" : "to", ustr(a));
    out->unit_kind = VIAL_UNIT_UNKNOWN;
}

void vx_unit_join(infer_t* x, vial_node_id at, const char* verb, const T* a, const T* b, T* out) { unit_join(x, at, verb, a, b, out); }

// The unit of * and /: a bare number has none
static void unit_product(const T* a, const T* b, bool divide, T* out) {
    if (a->unit_kind == VIAL_UNIT_UNKNOWN || b->unit_kind == VIAL_UNIT_UNKNOWN) { out->unit_kind = VIAL_UNIT_UNKNOWN; return; }
    if (a->unit_kind == VIAL_UNIT_OPEN && b->unit_kind == VIAL_UNIT_OPEN) { out->unit_kind = VIAL_UNIT_OPEN; return; }
    const md_unit_t ua = a->unit_kind == VIAL_UNIT_KNOWN ? a->unit : md_unit_none();
    const md_unit_t ub = b->unit_kind == VIAL_UNIT_KNOWN ? b->unit : md_unit_none();
    out->unit_kind = VIAL_UNIT_KNOWN;
    out->unit = divide ? md_unit_div(ua, ub) : md_unit_mul(ua, ub);
}

// ### VALUES ###

static vial_value_t val_none(void) { return (vial_value_t){0}; }
static vial_value_t val_int(int64_t i) { return (vial_value_t){.kind = VIAL_VAL_INT, .i = i}; }
static vial_value_t val_float(double f) { return (vial_value_t){.kind = VIAL_VAL_FLOAT, .f = f}; }
static vial_value_t val_bool(bool b) { return (vial_value_t){.kind = VIAL_VAL_BOOL, .b = b}; }
static bool val_num(const vial_value_t* v, double* d) {
    switch (v->kind) {
    case VIAL_VAL_INT: *d = (double)v->i; return true;
    case VIAL_VAL_FLOAT: *d = v->f; return true;
    case VIAL_VAL_BOOL: *d = v->b ? 1.0 : 0.0; return true;
    default: return false;
    }
}

size_t vial_value_print(char* buf, size_t cap, const vial_value_t* v) {
    int n = 0;
    switch (v->kind) {
    case VIAL_VAL_BOOL:  n = snprintf(buf, cap, "%s", v->b ? "True" : "False"); break;
    case VIAL_VAL_INT:   n = snprintf(buf, cap, "%lld", (long long)v->i); break;
    case VIAL_VAL_FLOAT: n = snprintf(buf, cap, "%.10g", v->f); break;
    case VIAL_VAL_STR:   n = snprintf(buf, cap, "'%.*s'", (int)v->s.len, v->s.ptr); break;
    case VIAL_VAL_NULL:  n = snprintf(buf, cap, "None"); break;
    default:             n = snprintf(buf, cap, "?"); break;
    }
    return n > 0 ? (size_t)n : 0;
}

static int64_t parse_int(str_t s) {
    char buf[64];
    size_t n = 0;
    for (size_t i = 0; i < s.len && n + 1 < sizeof(buf); ++i) if (s.ptr[i] != '_') buf[n++] = s.ptr[i];
    buf[n] = '\0';
    int base = 10;
    const char* p = buf;
    if (n > 2 && buf[0] == '0' && (buf[1] == 'x' || buf[1] == 'X')) { base = 16; p += 2; }
    else if (n > 2 && buf[0] == '0' && (buf[1] == 'o' || buf[1] == 'O')) { base = 8; p += 2; }
    else if (n > 2 && buf[0] == '0' && (buf[1] == 'b' || buf[1] == 'B')) { base = 2; p += 2; }
    return strtoll(p, NULL, base);
}

static double parse_float(str_t s) {
    char buf[64];
    size_t n = 0;
    for (size_t i = 0; i < s.len && n + 1 < sizeof(buf); ++i) if (s.ptr[i] != '_') buf[n++] = s.ptr[i];
    buf[n] = '\0';
    return strtod(buf, NULL);
}

// The text of string literals, adjacent ones joined; escapes other than \\, \' and \" are kept as written
static str_t string_value(infer_t* x, vial_node_id id) {
    const str_t s = text(x, id);
    char* out = NULL;
    for (size_t i = 0; i < s.len;) {
        while (i < s.len && s.ptr[i] != '\'' && s.ptr[i] != '"') ++i;
        if (i >= s.len) break;
        const char q = s.ptr[i];
        const bool triple = i + 2 < s.len && s.ptr[i + 1] == q && s.ptr[i + 2] == q;
        i += triple ? 3 : 1;
        while (i < s.len) {
            if (s.ptr[i] == '\\' && i + 1 < s.len) { md_array_push(out, s.ptr[i + 1], x->alloc); i += 2; continue; }
            if (s.ptr[i] == q && (!triple || (i + 2 < s.len && s.ptr[i + 1] == q && s.ptr[i + 2] == q))) { i += triple ? 3 : 1; break; }
            md_array_push(out, s.ptr[i], x->alloc);
            ++i;
        }
    }
    return (str_t){out, md_array_size(out)};
}

// ### EXPRESSIONS ###

static T infer(infer_t* x, vial_node_id id);

static void set_value(infer_t* x, vial_node_id id, vial_value_t v) { x->out->values[(id) - x->out->base] = v; }
static vial_value_t value_of(const infer_t* x, vial_node_id id) { return x->out->values[(id) - x->out->base]; }
void vx_set_value(infer_t* x, vial_node_id id, vial_value_t v) { set_value(x, id, v); }
vial_value_t vx_value_of(const infer_t* x, vial_node_id id) { return value_of(x, id); }

static T infer_name(infer_t* x, vial_node_id id) {
    const vial_symbol_t* s = vial_names_symbol(x->n, id);
    if (!s) return err();
    const uint32_t si = x->n->ref[id];
    T ft;
    vial_value_t fv;
    if (x->frame && s->kind != VIAL_SYM_BUILTIN && vx_frame_get(x, si, &ft, &fv)) {
        set_value(x, id, fv);
        return ft;
    }
    switch (s->kind) {
    case VIAL_SYM_GLOBAL:
        if (x->t != x->c->tree) return err();
        set_value(x, id, x->top->sym_values[si]);
        return x->top->sym_types[si];
    case VIAL_SYM_DEF:
    case VIAL_SYM_PRELUDE:
        return scalar(VIAL_DT_DEF);
    case VIAL_SYM_BUILTIN: {
        T t;
        if (s->builtin && vx_name_sets(x, id, s->builtin, &t)) return t;
        return err();
    }
    default:
        return err();   // built-in values, parameters and locals arrive with later steps
    }
}

static T infer_quantity(infer_t* x, vial_node_id id) {
    const vial_node_id num = child(x, id, 0), unit = child(x, id, 1);
    T t = infer(x, num);
    if (vial_type_is_error(&t)) return t;
    md_unit_t u;
    const str_t ut = text(x, unit);
    bool ok;
    if (sub_of(x, unit) == 1) {
        // The text between the quotes, past a prefix
        size_t b = 0;
        while (b < ut.len && ut.ptr[b] != '\'' && ut.ptr[b] != '"') ++b;
        const size_t q = (b + 2 < ut.len && ut.ptr[b + 1] == ut.ptr[b] && ut.ptr[b + 2] == ut.ptr[b]) ? 3 : 1;
        ok = ut.len >= b + 2 * q && vial_unit_string(str_substr(ut, b + q, ut.len - b - 2 * q), &u);
    } else {
        ok = vial_unit_name(ut, &u);
    }
    if (!ok) return err();  // E0502, from the parser
    t.unit_kind = VIAL_UNIT_KNOWN;
    t.unit = u;
    set_value(x, id, value_of(x, num));
    return t;
}

static T infer_unary(infer_t* x, vial_node_id id) {
    const vial_node_id on = child(x, id, 0);
    const T a = infer(x, on);
    if (vial_type_is_error(&a)) return a;
    const uint16_t op = sub_of(x, id);
    if (!is_number(&a) || (op == VIAL_OP_TILDE && a.dtype != VIAL_DT_INT && a.dtype != VIAL_DT_BOOL)) {
        diag(x, "E0319", id, "`%s` takes %s; this is %s", vial_op_str((vial_op_t)op), op == VIAL_OP_TILDE ? "ints" : "numbers", dtype_word(&a));
        return err();
    }
    T t = a;
    if (t.dtype == VIAL_DT_BOOL) t.dtype = VIAL_DT_INT;
    t.flags &= ~VIAL_TYPE_PERIODIC;
    if (op == VIAL_OP_MINUS && (a.flags & VIAL_TYPE_PERIODIC)) {     // negation keeps a period, mirrored
        t.flags |= VIAL_TYPE_PERIODIC;
        t.period_beg = -a.period_end; t.period_end = -a.period_beg; t.seam = -a.seam;
    }
    const vial_value_t v = value_of(x, on);
    if (op == VIAL_OP_MINUS && v.kind == VIAL_VAL_INT) set_value(x, id, v.i > INT64_MIN ? val_int(-v.i) : val_none());
    else if (op == VIAL_OP_MINUS && v.kind == VIAL_VAL_FLOAT) set_value(x, id, val_float(-v.f));
    else if (op == VIAL_OP_PLUS) set_value(x, id, v);
    else if (op == VIAL_OP_TILDE && v.kind == VIAL_VAL_INT) set_value(x, id, val_int(~v.i));
    return t;
}

static uint8_t promote(const T* a, const T* b) {
    // A bare number takes the dtype of the array it meets when both are integers or both floats (NumPy 2)
    if (a->dtype == b->dtype) return a->dtype;
    uint8_t da = a->dtype == VIAL_DT_BOOL ? VIAL_DT_INT : a->dtype;
    uint8_t db = b->dtype == VIAL_DT_BOOL ? VIAL_DT_INT : b->dtype;
    if (da == db) return da;
    if (is_float(da) && is_float(db)) {
        // float64 meets float: 32 bits, unless the float side is a bare number
        if (da == VIAL_DT_FLOAT64) return b->unit_kind == VIAL_UNIT_OPEN && b->ndim == 0 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT;
        return a->unit_kind == VIAL_UNIT_OPEN && a->ndim == 0 ? VIAL_DT_FLOAT64 : VIAL_DT_FLOAT;
    }
    if (is_float(da)) return da;
    if (is_float(db)) return db;
    return VIAL_DT_INT;
}

uint8_t vx_promote(const T* a, const T* b) { return promote(a, b); }

// A value of b's type in a's unit, for + - % and comparisons; false when the units do not convert
static bool convert_value(const T* a, const T* b, vial_value_t* vb) {
    if (a->unit_kind != VIAL_UNIT_KNOWN || b->unit_kind != VIAL_UNIT_KNOWN) return true;   // a bare number takes the other's
    double f;
    if (!md_unit_conversion_factor(&f, b->unit, a->unit)) return false;
    if (f == 1.0) return true;
    double d;
    if (!val_num(vb, &d)) return true;
    *vb = val_float(d * f);
    return true;
}

// An int of a value: exact for ints and bools, and for a float only when it is whole and in range
static bool val_int_exact(const vial_value_t* v, int64_t* out) {
    switch (v->kind) {
    case VIAL_VAL_INT:  *out = v->i; return true;
    case VIAL_VAL_BOOL: *out = v->b; return true;
    case VIAL_VAL_FLOAT:
        if (!(v->f >= -9.0e15 && v->f <= 9.0e15) || v->f != floor(v->f)) return false;
        *out = (int64_t)v->f;
        return true;
    default: return false;
    }
}

// An int result is folded only when it fits in VIAL's 32 bits; beyond that it is R0102 when evaluated
static vial_value_t val_int32(int64_t i) { return i >= INT32_MIN && i <= INT32_MAX ? val_int(i) : val_none(); }

static vial_value_t fold(uint16_t op, const vial_value_t* a, const vial_value_t* b, uint8_t dtype) {
    double fa, fb;
    if (!val_num(a, &fa) || !val_num(b, &fb)) return val_none();
    if (dtype != VIAL_DT_INT) {
        double r;
        switch (op) {
        case VIAL_OP_PLUS:  r = fa + fb; break;
        case VIAL_OP_MINUS: r = fa - fb; break;
        case VIAL_OP_STAR:  r = fa * fb; break;
        case VIAL_OP_SLASH: if (fb == 0.0) return val_none(); r = fa / fb; break;
        case VIAL_OP_DSLASH: if (fb == 0.0) return val_none(); r = floor(fa / fb); break;
        case VIAL_OP_PERCENT: if (fb == 0.0) return val_none(); r = fa - fb * floor(fa / fb); break;
        case VIAL_OP_POW:   r = pow(fa, fb); break;
        default: return val_none();
        }
        return isfinite(r) ? val_float(r) : val_none();
    }
    // Ints, exactly; operands of 32 bits keep every result below in 64
    int64_t ia, ib;
    if (!val_int_exact(a, &ia) || !val_int_exact(b, &ib)) return val_none();
    if (ia < INT32_MIN || ia > INT32_MAX || ib < INT32_MIN || ib > INT32_MAX) return val_none();
    switch (op) {
    case VIAL_OP_PLUS:  return val_int32(ia + ib);
    case VIAL_OP_MINUS: return val_int32(ia - ib);
    case VIAL_OP_STAR:  return val_int32(ia * ib);
    case VIAL_OP_SLASH: return ib ? val_float((double)ia / (double)ib) : val_none();
    case VIAL_OP_DSLASH: return ib ? val_int32((ia - ((ia % ib) + ib) % ib) / ib) : val_none();
    case VIAL_OP_PERCENT: return ib ? val_int32(((ia % ib) + ib) % ib) : val_none();
    case VIAL_OP_POW: {
        if (ib < 0) return val_float(pow((double)ia, (double)ib));
        int64_t r = 1;
        for (int64_t k = 0; k < ib; ++k) {
            r *= ia;
            if (r < INT32_MIN || r > INT32_MAX) return val_none();
            if (r == 0 || r == 1) break;
        }
        return val_int32(r);
    }
    case VIAL_OP_AMP:    return val_int32(ia & ib);
    case VIAL_OP_PIPE:   return val_int32(ia | ib);
    case VIAL_OP_CARET:  return val_int32(ia ^ ib);
    case VIAL_OP_LSHIFT: return ib >= 0 && ib < 32 ? val_int32((int64_t)((uint64_t)ia << ib)) : val_none();
    case VIAL_OP_RSHIFT: return ib >= 0 && ib < 63 ? val_int32(ia >> ib) : val_none();
    default: return val_none();
    }
}

static T infer_matmul(infer_t* x, vial_node_id id, const T* a, const T* b) {
    // NumPy's matmul: a vector is promoted and its dimension removed again; leading dimensions broadcast
    if (a->ndim == 0 || b->ndim == 0) {
        diag(x, "E0319", id, "`@` takes vectors and matrices; this is %s and %s", tstr(a), tstr(b));
        return err();
    }
    T aa = *a, bb = *b;
    const bool va = a->ndim == 1, vb = b->ndim == 1;
    if (va) { aa.ndim = 2; aa.dims[1] = a->dims[0]; aa.dims[0] = vial_extent_static(1); }
    if (vb) { bb.ndim = 2; bb.dims[1] = vial_extent_static(1); bb.dims[0] = b->dims[0]; }
    vial_extent_t inner;
    if (!extent_join(aa.dims[aa.ndim - 1], bb.dims[bb.ndim - 2], &inner) ||
        (aa.dims[aa.ndim - 1].kind == VIAL_EXT_STATIC && aa.dims[aa.ndim - 1].n == 1 && !(bb.dims[bb.ndim - 2].kind == VIAL_EXT_STATIC && bb.dims[bb.ndim - 2].n == 1))) {
        diag(x, "E0301", id, "`%s` @ `%s`: the inner dimensions differ", tstr(a), tstr(b));
        return err();
    }
    T la = aa, lb = bb;
    la.ndim -= 2; lb.ndim -= 2;
    T t = scalar(promote(a, b));
    if (!broadcast(&la, &lb, &t)) {
        diag(x, "E0301", id, "`%s` and `%s` do not broadcast", tstr(a), tstr(b));
        return err();
    }
    if (!va) vial_type_push_dim(&t, aa.dims[aa.ndim - 2]);
    if (!vb) vial_type_push_dim(&t, bb.dims[bb.ndim - 1]);
    unit_product(a, b, false, &t);
    join_stage(&t, a, b);
    return t;
}

static T infer_binary(infer_t* x, vial_node_id id) {
    const vial_node_id an = child(x, id, 0), bn = child(x, id, 1);
    const T a = infer(x, an);
    const T b = infer(x, bn);
    if (vial_type_is_error(&a) || vial_type_is_error(&b)) return err();
    const uint16_t op = sub_of(x, id);

    // + on drawings combines them
    if (op == VIAL_OP_PLUS && a.dtype == VIAL_DT_DRAWING && b.dtype == VIAL_DT_DRAWING) {
        T t = scalar(VIAL_DT_DRAWING);
        if (!broadcast_or_report(x, id, an, &a, bn, &b, &t)) return err();
        join_stage(&t, &a, &b);
        return t;
    }
    const bool bits = op == VIAL_OP_AMP || op == VIAL_OP_PIPE || op == VIAL_OP_CARET || op == VIAL_OP_LSHIFT || op == VIAL_OP_RSHIFT;
    if (!is_number(&a) || !is_number(&b) || (bits && (is_float(a.dtype) || is_float(b.dtype)))) {
        if (a.dtype == VIAL_DT_SET && b.dtype == VIAL_DT_SET && (op == VIAL_OP_AMP || op == VIAL_OP_PIPE)) {
            diag_fix(x, "E0319", id, op == VIAL_OP_AMP ? "and" : "or", "sets combine with `and`, `or` and `not`, not `%s`", vial_op_str((vial_op_t)op));
        } else {
            diag(x, "E0319", id, "`%s` takes %s; this is %s and %s", vial_op_str((vial_op_t)op), bits ? "ints" : "numbers",
                dtype_word(&a), dtype_word(&b));
        }
        return err();
    }
    if (op == VIAL_OP_AT) return infer_matmul(x, id, &a, &b);

    T t = scalar(promote(&a, &b));
    if (t.dtype == VIAL_DT_BOOL) t.dtype = VIAL_DT_INT;    // arithmetic on bools counts
    if (!broadcast_or_report(x, id, an, &a, bn, &b, &t)) return err();
    join_stage(&t, &a, &b);
    switch (op) {
    case VIAL_OP_PLUS:    unit_join(x, id, "adds", &a, &b, &t); break;
    case VIAL_OP_MINUS:   unit_join(x, id, "subtracts", &a, &b, &t); break;
    case VIAL_OP_PERCENT: unit_join(x, id, "takes the remainder of", &a, &b, &t); break;
    case VIAL_OP_STAR:    unit_product(&a, &b, false, &t); break;
    case VIAL_OP_SLASH:   unit_product(&a, &b, true, &t); if (!is_float(t.dtype)) t.dtype = VIAL_DT_FLOAT; break;
    case VIAL_OP_DSLASH:  unit_product(&a, &b, true, &t); break;
    case VIAL_OP_POW: {
        const vial_value_t e = value_of(x, bn);
        if (b.unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(b.unit)) {
            diag_fix(x, "W0501", id, "unitless(...)", "a power of %s; an exponent has no unit", ustr(&b));
            t.unit_kind = VIAL_UNIT_UNKNOWN;
        } else if (a.unit_kind == VIAL_UNIT_KNOWN && !md_unit_is_none(a.unit)) {
            if (e.kind == VIAL_VAL_INT) { t.unit_kind = VIAL_UNIT_KNOWN; t.unit = md_unit_pow(a.unit, (int)e.i); }
            else if (e.kind == VIAL_VAL_FLOAT && e.f == floor(e.f)) { t.unit_kind = VIAL_UNIT_KNOWN; t.unit = md_unit_pow(a.unit, (int)e.f); }
            else {
                diag_fix(x, "W0501", id, "unitless(...)", "%s to a power not known when compiling; the result's unit becomes unknown", ustr(&a));
                t.unit_kind = VIAL_UNIT_UNKNOWN;
            }
        } else {
            t.unit_kind = a.unit_kind == VIAL_UNIT_UNKNOWN ? VIAL_UNIT_UNKNOWN : a.unit_kind;
        }
        if (e.kind == VIAL_VAL_INT && e.i < 0 && t.dtype == VIAL_DT_INT) t.dtype = VIAL_DT_FLOAT;
        break;
    }
    default:              t.unit_kind = VIAL_UNIT_OPEN; break;   // bits: integers without a unit
    }
    // Periods: moved or scaled by a constant, or the sum or difference of two of the same period
    const bool pa = a.flags & VIAL_TYPE_PERIODIC, pb = b.flags & VIAL_TYPE_PERIODIC;
    if ((op == VIAL_OP_PLUS || op == VIAL_OP_MINUS) && pa && pb && a.period_end - a.period_beg == b.period_end - b.period_beg) {
        t.flags |= VIAL_TYPE_PERIODIC;
        const float w = a.period_end - a.period_beg;
        t.period_beg = -w / 2; t.period_end = w / 2; t.seam = w / 2;    // a difference of angles is cut half a period from 0
    } else if ((op == VIAL_OP_PLUS || op == VIAL_OP_MINUS) && (pa != pb) && (pa ? b.stage : a.stage) == VIAL_STAGE_STATIC) {
        const T* p = pa ? &a : &b;
        double c;
        if (val_num(pa ? &x->out->values[(bn) - x->out->base] : &x->out->values[(an) - x->out->base], &c)) {
            const float s = (float)(op == VIAL_OP_MINUS && pa ? -c : c);
            t.flags |= VIAL_TYPE_PERIODIC;
            t.period_beg = p->period_beg + s; t.period_end = p->period_end + s; t.seam = p->seam + s;
        }
    } else if ((op == VIAL_OP_STAR || op == VIAL_OP_SLASH) && (pa != pb)) {
        double c;
        if (val_num(pa ? &x->out->values[(bn) - x->out->base] : &x->out->values[(an) - x->out->base], &c) && c != 0.0 && (op == VIAL_OP_STAR || pa)) {
            const T* p = pa ? &a : &b;
            const float k = (float)(op == VIAL_OP_SLASH ? 1.0 / c : c);
            t.flags |= VIAL_TYPE_PERIODIC;
            const float lo = p->period_beg * k, hi = p->period_end * k;
            t.period_beg = lo < hi ? lo : hi; t.period_end = lo < hi ? hi : lo; t.seam = p->seam * k;
        }
    }
    vial_value_t va = value_of(x, an), vb = value_of(x, bn);
    if (op == VIAL_OP_PLUS || op == VIAL_OP_MINUS || op == VIAL_OP_PERCENT) {
        if (t.unit_kind == VIAL_UNIT_UNKNOWN || !convert_value(&a, &b, &vb)) { set_value(x, id, val_none()); return t; }
    }
    set_value(x, id, fold(op, &va, &vb, t.dtype));
    return t;
}

static bool compare_values(uint16_t cmp, const vial_value_t* a, const vial_value_t* b, bool* out) {
    double fa, fb;
    if (a->kind == VIAL_VAL_STR && b->kind == VIAL_VAL_STR) {
        const bool eq = str_eq(a->s, b->s);
        if (cmp == VIAL_CMP_EQ) { *out = eq; return true; }
        if (cmp == VIAL_CMP_NE) { *out = !eq; return true; }
        return false;
    }
    if (!val_num(a, &fa) || !val_num(b, &fb)) return false;
    switch (cmp) {
    case VIAL_CMP_LT: *out = fa < fb; return true;
    case VIAL_CMP_GT: *out = fa > fb; return true;
    case VIAL_CMP_LE: *out = fa <= fb; return true;
    case VIAL_CMP_GE: *out = fa >= fb; return true;
    case VIAL_CMP_EQ: *out = fa == fb; return true;
    case VIAL_CMP_NE: *out = fa != fb; return true;
    default: return false;
    }
}

static T infer_compare(infer_t* x, vial_node_id id) {
    // [operand, CMPOP, operand, CMPOP, operand, ...]
    const uint32_t n = num_children(x, id);
    T result = scalar(VIAL_DT_BOOL);
    bool all_known = true, value = true;
    bool failed = false;
    for (uint32_t i = 0; i + 2 < n; i += 2) {
        const vial_node_id an = child(x, id, i), on = child(x, id, i + 1), bn = child(x, id, i + 2);
        const T a = infer(x, an);
        const T b = infer(x, bn);
        if (vial_type_is_error(&a) || vial_type_is_error(&b)) { failed = true; continue; }
        const uint16_t cmp = sub_of(x, on);
        T pair = scalar(VIAL_DT_BOOL);
        if (cmp == VIAL_CMP_IS || cmp == VIAL_CMP_IS_NOT) {
            // is compares static values, such as a parameter with None
            if (a.stage != VIAL_STAGE_STATIC || b.stage != VIAL_STAGE_STATIC) {
                diag_fix(x, "E0317", on, cmp == VIAL_CMP_IS ? "==" : "!=", "`is` compares values known when compiling, and this changes per frame");
                failed = true;
                continue;
            }
            join_stage(&pair, &a, &b);
        } else if ((cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) && b.dtype == VIAL_DT_TUPLE) {
            // One of a tuple's values, as assert op in ('closest', 'farthest')
            join_stage(&pair, &a, &b);
        } else if (cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) {
            if (b.dtype == VIAL_DT_RANGE) {
                if (!is_number(&a)) { diag(x, "E0319", on, "`in` a range takes numbers; this is %s", dtype_word(&a)); failed = true; continue; }
                pair.ndim = a.ndim; memcpy(pair.dims, a.dims, sizeof(a.dims));
                T u = pair; unit_join(x, on, "compares", &a, &b, &u);
            } else if (b.ndim >= 1) {
                pair.ndim = a.ndim; memcpy(pair.dims, a.dims, sizeof(a.dims));
            } else {
                diag(x, "E0319", on, "`in` takes a range or an array; this is %s", tstr(&b)); failed = true; continue;
            }
            join_stage(&pair, &a, &b);
        } else {
            const bool ordered = cmp != VIAL_CMP_EQ && cmp != VIAL_CMP_NE;
            const bool both_str = a.dtype == VIAL_DT_STR && b.dtype == VIAL_DT_STR;
            const bool with_none = !ordered && (a.dtype == VIAL_DT_NONE || b.dtype == VIAL_DT_NONE) && a.ndim == 0 && b.ndim == 0;
            if (with_none) {
                // None equals only None, as `weights == 'mass'` with weights=None
                join_stage(&pair, &a, &b);
                const vial_value_t va = value_of(x, an), vb = value_of(x, bn);
                if (all_known && va.kind != VIAL_VAL_NONE && vb.kind != VIAL_VAL_NONE) {
                    const bool eq = va.kind == vb.kind;
                    value = value && (cmp == VIAL_CMP_EQ ? eq : !eq);
                } else {
                    all_known = false;
                }
                T merged0 = result;
                if (broadcast(&result, &pair, &merged0)) { join_stage(&merged0, &result, &pair); result = merged0; }
                continue;
            }
            if (!(is_number(&a) && is_number(&b)) && !(both_str && !ordered) &&
                !(!ordered && a.dtype == b.dtype && (a.dtype == VIAL_DT_SET || a.dtype == VIAL_DT_BOOL))) {
                diag(x, "E0319", on, "`%s` takes %s; this is %s and %s", vial_cmp_str((vial_cmp_t)cmp), ordered ? "numbers" : "two values of one kind",
                    dtype_word(&a), dtype_word(&b));
                failed = true;
                continue;
            }
            if (!broadcast_or_report(x, on, an, &a, bn, &b, &pair)) { failed = true; continue; }
            join_stage(&pair, &a, &b);
            if (is_number(&a)) { T u = pair; unit_join(x, on, "compares", &a, &b, &u); }
        }
        T merged = result;
        if (!broadcast(&result, &pair, &merged)) {
            diag(x, "E0301", id, "`%s` and `%s` do not broadcast", tstr(&result), tstr(&pair));
            failed = true;
            continue;
        }
        join_stage(&merged, &result, &pair);
        result = merged;
        vial_value_t va = value_of(x, an), vb = value_of(x, bn);
        bool r;
        if (cmp == VIAL_CMP_IS || cmp == VIAL_CMP_IS_NOT) {
            // Known when either side is None and both are known, or both are static strings or bools
            const bool known = va.kind != VIAL_VAL_NONE && vb.kind != VIAL_VAL_NONE && (va.kind == VIAL_VAL_NULL || vb.kind == VIAL_VAL_NULL);
            if (all_known && known) {
                r = va.kind == vb.kind;
                value = value && (cmp == VIAL_CMP_IS ? r : !r);
            } else {
                all_known = false;
            }
            continue;
        }
        if ((cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) && b.dtype == VIAL_DT_TUPLE) {
            bool found = false, known = all_known && va.kind != VIAL_VAL_NONE && kind_of(x, bn) == VIAL_N_TUPLE;
            for (uint32_t k = 0; known && k < num_children(x, bn); ++k) {
                const vial_value_t e = value_of(x, child(x, bn, k));
                if (e.kind == VIAL_VAL_NONE) { known = false; break; }
                bool eq = false;
                if (e.kind == VIAL_VAL_NULL || va.kind == VIAL_VAL_NULL) eq = e.kind == va.kind;
                else if (!compare_values(VIAL_CMP_EQ, &va, &e, &eq)) eq = false;
                found = found || eq;
            }
            if (known) value = value && (cmp == VIAL_CMP_IN ? found : !found);
            else all_known = false;
            continue;
        }
        if (all_known && is_number(&a) && !convert_value(&a, &b, &vb)) all_known = false;
        if (all_known && compare_values(cmp, &va, &vb, &r)) value = value && r;
        else all_known = false;
    }
    if (failed) return err();
    result.dtype = VIAL_DT_BOOL;
    result.unit_kind = VIAL_UNIT_KNOWN;
    result.unit = md_unit_none();
    if (all_known && result.ndim == 0) set_value(x, id, val_bool(value));
    return result;
}

static T infer_boolop(infer_t* x, vial_node_id id) {
    const uint32_t n = num_children(x, id);
    T t = err();
    bool known = true, v = sub_of(x, id) == VIAL_KW_AND;
    vial_value_t sets = {0};
    for (uint32_t i = 0; i < n; ++i) {
        const vial_node_id cn = child(x, id, i);
        const T c = infer(x, cn);
        if (vial_type_is_error(&c)) return err();
        if (c.dtype != VIAL_DT_BOOL && c.dtype != VIAL_DT_SET) {
            diag(x, "E0319", cn, "`%s` takes bools or sets; this is %s", sub_of(x, id) == VIAL_KW_AND ? "and" : "or", dtype_word(&c));
            return err();
        }
        if (i == 0) { t = c; }
        else {
            if (c.dtype != t.dtype) {
                diag(x, "E0319", id, "`%s` takes two bools or two sets; this is %s and %s", sub_of(x, id) == VIAL_KW_AND ? "and" : "or",
                    dtype_word(&t), dtype_word(&c));
                return err();
            }
            T m = t;
            if (!broadcast_or_report(x, id, child(x, id, 0), &t, cn, &c, &m)) return err();
            join_stage(&m, &t, &c);
            t = m;
        }
        const vial_value_t cv = value_of(x, cn);
        if (cv.kind == VIAL_VAL_BOOL) v = sub_of(x, id) == VIAL_KW_AND ? (v && cv.b) : (v || cv.b);
        else known = false;
        if (c.dtype == VIAL_DT_SET && !(sub_of(x, id) == VIAL_KW_AND && kind_of(x, cn) == VIAL_N_NOT)) {
            sets = (i == 0 || sets.kind == VIAL_VAL_NULL) ? cv : vx_sets_boolop(x, sub_of(x, id) == VIAL_KW_AND, &sets, &cv);
        } else if (c.dtype == VIAL_DT_SET && i == 0) {
            sets = (vial_value_t){.kind = VIAL_VAL_NULL};     // only negated operands so far: subtracted below
        }
    }
    if (t.dtype == VIAL_DT_SET && sub_of(x, id) == VIAL_KW_AND) {
        // The negated operands, subtracted from what the others give, or from the scope when there are none
        for (uint32_t i = 0; i < n; ++i) {
            const vial_node_id cn = child(x, id, i);
            if (kind_of(x, cn) != VIAL_N_NOT) continue;
            const vial_value_t b = value_of(x, child(x, cn, 0));
            if (b.kind != VIAL_VAL_SETS || !b.sets->single || !x->sys) { sets = val_none(); break; }
            if (sets.kind == VIAL_VAL_NULL) {
                if (!x->all) x->all = vial_sets_all(x->sys->num_atoms, x->alloc);
                sets = (vial_value_t){.kind = VIAL_VAL_SETS, .sets = x->all};
            }
            if (sets.kind != VIAL_VAL_SETS || !sets.sets->single) { sets = val_none(); break; }
            sets = (vial_value_t){.kind = VIAL_VAL_SETS, .sets = vial_sets_not(sets.sets, b.sets, 0, x->alloc)};
        }
    }
    if (t.dtype == VIAL_DT_SET) {
        t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
        if (sets.kind == VIAL_VAL_SETS && t.stage == VIAL_STAGE_STATIC) {
            const T st = vx_sets_type(x, sets.sets, id);
            t.flags |= st.flags & (VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
            set_value(x, id, sets);
        }
    }
    if (known && t.dtype == VIAL_DT_BOOL && t.ndim == 0) set_value(x, id, val_bool(v));
    return t;
}

static T infer_not(infer_t* x, vial_node_id id) {
    const vial_node_id on = child(x, id, 0);
    const T a = infer(x, on);
    if (vial_type_is_error(&a)) return a;
    if (a.dtype != VIAL_DT_BOOL && a.dtype != VIAL_DT_SET) {
        diag(x, "E0319", id, "`not` takes bools or sets; this is %s", dtype_word(&a));
        return err();
    }
    const vial_value_t v = value_of(x, on);
    if (v.kind == VIAL_VAL_BOOL) set_value(x, id, val_bool(!v.b));
    T t = a;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    // Under an `and` of sets, `not b` is not made: the `and` subtracts b, so `protein and not backbone` costs the
    // size of the sets, not the system's
    const vial_node_id up = x->n->parent[id];
    const bool under_and = up && kind_of(x, up) == VIAL_N_BOOLOP && sub_of(x, up) == VIAL_KW_AND;
    if (a.dtype == VIAL_DT_SET && v.kind == VIAL_VAL_SETS && !under_and) {
        const vial_value_t nv = vx_sets_not(x, &v);
        if (nv.kind == VIAL_VAL_SETS) {
            const T st = vx_sets_type(x, nv.sets, id);
            t.flags |= st.flags & (VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
            set_value(x, id, nv);
        }
    }
    return t;
}

// Two values of one type: a list's items, a conditional's branches. Units convert to the first's.
static bool join_types(infer_t* x, vial_node_id at, const char* what, const T* a, const T* b, T* out) {
    if (a->dtype != b->dtype && !(is_number(a) && is_number(b))) {
        diag(x, "E0319", at, "%s have one type; these are %s and %s", what, tstr(a), tstr(b));
        return false;
    }
    T t = *a;
    if (is_number(a)) t.dtype = promote(a, b);
    if (a->ndim != b->ndim) {
        diag(x, "E0301", at, "%s have one shape; these are `%s` and `%s`", what, tstr(a), tstr(b));
        return false;
    }
    for (uint8_t k = 0; k < a->ndim; ++k) {
        if (!extent_join(a->dims[k], b->dims[k], &t.dims[k]) ||
            (a->dims[k].kind == VIAL_EXT_STATIC && b->dims[k].kind == VIAL_EXT_STATIC && a->dims[k].n != b->dims[k].n)) {
            diag(x, "E0301", at, "%s have one shape; these are `%s` and `%s`", what, tstr(a), tstr(b));
            return false;
        }
    }
    if (is_number(a)) unit_join(x, at, "joins", a, b, &t);
    join_stage(&t, a, b);
    *out = t;
    return true;
}

static T infer_list(infer_t* x, vial_node_id id) {
    const uint32_t n = num_children(x, id);
    if (n == 0) { T t = number(VIAL_DT_FLOAT); vial_type_push_dim(&t, vial_extent_static(0)); return t; }
    T item = infer(x, child(x, id, 0));
    if (vial_type_is_error(&item)) return err();
    for (uint32_t i = 1; i < n; ++i) {
        const T c = infer(x, child(x, id, i));
        if (vial_type_is_error(&c)) return err();
        T m;
        if (!join_types(x, child(x, id, i), "the items of a list", &item, &c, &m)) return err();
        item = m;
    }
    if (item.dtype == VIAL_DT_TUPLE || item.dtype == VIAL_DT_DEF || item.dtype == VIAL_DT_NONE || item.dtype == VIAL_DT_RANGE) {
        diag(x, "E0319", id, "a list holds numbers, strings, sets or drawings; this holds %s", dtype_word(&item));
        return err();
    }
    T t = item;
    t.ndim = 0;
    t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
    if (!vial_type_push_dim(&t, vial_extent_static(n))) goto too_many;
    for (uint8_t k = 0; k < item.ndim; ++k) if (!vial_type_push_dim(&t, item.dims[k])) goto too_many;
    if (item.dtype == VIAL_DT_SET && (item.flags & VIAL_TYPE_ONE_PARTICLE)) t.flags |= VIAL_TYPE_ONE_PARTICLE;
    return t;
too_many:
    diag(x, "E0308", id, "the result would have 7 dimensions; the limit is 6");
    return err();
}

static T infer_tuple(infer_t* x, vial_node_id id) {
    const uint32_t n = num_children(x, id);
    T* elems = md_alloc(x->alloc, (n ? n : 1) * sizeof(T));
    T t = scalar(VIAL_DT_TUPLE);
    bool failed = false;
    for (uint32_t i = 0; i < n; ++i) {
        elems[i] = infer(x, child(x, id, i));
        if (vial_type_is_error(&elems[i])) failed = true;
        join_stage(&t, &t, &elems[i]);
    }
    if (failed) return err();
    t.elems = elems;
    t.num_elems = n;
    return t;
}

static T infer_range(infer_t* x, vial_node_id id) {
    const vial_node_id an = child(x, id, 0), bn = child(x, id, 1);
    const T a = infer(x, an);
    const T b = infer(x, bn);
    if (vial_type_is_error(&a) || vial_type_is_error(&b)) return err();
    if (!is_number(&a) || !is_number(&b) || a.ndim || b.ndim) {
        diag(x, "E0319", id, "`..` takes two numbers; this is %s and %s", tstr(&a), tstr(&b));
        return err();
    }
    T t = scalar(VIAL_DT_RANGE);
    unit_join(x, id, "joins", &a, &b, &t);
    join_stage(&t, &a, &b);
    return t;
}

static T infer_ternary(infer_t* x, vial_node_id id) {
    const vial_node_id vn = child(x, id, 0), cn = child(x, id, 1), en = child(x, id, 2);
    const T c = infer(x, cn);
    if (vial_type_is_error(&c)) return err();
    if (c.dtype != VIAL_DT_BOOL || c.ndim) {
        diag_fix(x, "E0303", cn, "where(...)", "a conditional needs one bool; this is `%s`", tstr(&c));
        return err();
    }
    const vial_value_t cv = value_of(x, cn);
    if (cv.kind == VIAL_VAL_BOOL && c.stage == VIAL_STAGE_STATIC) {
        // A static condition: the branch it does not take is not type-checked
        const vial_node_id taken = cv.b ? vn : en;
        const T t = infer(x, taken);
        set_value(x, id, value_of(x, taken));
        return t;
    }
    // Lowering: the branches are regions of an if, so only the one taken is evaluated, as Python's conditional
    vial_ir_t* ir = x->ir;
    const vial_ssa_t cs = ir ? vx_ssa(x, cn) : 0;
    uint32_t r1 = 0, r2 = 0, saved = 0;
    if (ir && cs && c.ndim == 0) {
        saved = ir->insert;
        r1 = vial_ir_region(ir, VIAL_IR_NO_OP, ir->regions[saved].func, NULL, 0);
        ir->insert = r1;
    }
    const T a = infer(x, vn);
    if (r1) {
        const vial_ssa_t ya = vx_ssa(x, vn);
        vial_ir_yield(ir, &ya, ya ? 1 : 0, vx_span(x, vn));
        r2 = vial_ir_region(ir, VIAL_IR_NO_OP, ir->regions[saved].func, NULL, 0);
        ir->insert = r2;
    }
    const T b = infer(x, en);
    vial_ssa_t yb = 0;
    if (r1) {
        yb = vx_ssa(x, en);
        vial_ir_yield(ir, &yb, yb ? 1 : 0, vx_span(x, en));
        ir->insert = saved;
    }
    if (vial_type_is_error(&a) || vial_type_is_error(&b)) return err();
    T t;
    if (!join_types(x, id, "the branches of a conditional", &a, &b, &t)) return err();
    join_stage(&t, &t, &c);
    if (r1 && vx_ssa(x, vn) && yb && a.dtype == b.dtype) {
        const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_IF, .args = &cs, .num_args = 1, .types = &t, .num_results = 1,
                                                          .span = vx_span(x, id)});
        vial_ir_attach(ir, op, 0, r1);
        vial_ir_attach(ir, op, 1, r2);
        vx_set_ssa(x, id, ir->ops[op].result);
    }
    return t;
}

// ### INDEXING ###

static bool static_index(infer_t* x, vial_node_id id, int64_t* out) {
    const vial_value_t v = value_of(x, id);
    if (v.kind == VIAL_VAL_INT) { *out = v.i; return true; }
    return false;
}

static T infer_subscript(infer_t* x, vial_node_id id) {
    const vial_node_id on = child(x, id, 0), in = child(x, id, 1);
    const T a = infer(x, on);
    if (vial_type_is_error(&a)) { if (in) infer(x, in); return err(); }
    if (a.dtype == VIAL_DT_SET && a.ndim == 0) {
        // Positions of a set's particles: a set
        if (!in || kind_of(x, in) == VIAL_N_TUPLE) {
            if (in) infer(x, in);
            diag(x, "E0319", id, "a set takes one position, a slice or a list of positions");
            return err();
        }
        T t = a;
        t.flags &= ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
        if (kind_of(x, in) == VIAL_N_SLICE) {
            for (uint32_t k = 0; k < 3; ++k) {
                const vial_node_id b = child(x, in, k);
                if (!b) continue;
                const T bt = infer(x, b);
                if (vial_type_is_error(&bt)) return err();
                join_stage(&t, &t, &bt);
            }
        } else {
            const T ix = infer(x, in);
            if (vial_type_is_error(&ix)) return err();
            if (ix.dtype == VIAL_DT_INT && ix.ndim == 0) t.flags |= VIAL_TYPE_ONE_PARTICLE;
            else if (!(ix.dtype == VIAL_DT_INT && ix.ndim == 1)) {
                diag(x, "E0319", in, "a set's particles are taken by an int, a slice or a list of ints; this is %s", tstr(&ix));
                return err();
            }
            join_stage(&t, &t, &ix);
        }
        const vial_value_t av = value_of(x, on);
        if (av.kind == VIAL_VAL_SETS && t.stage == VIAL_STAGE_STATIC) {
            const vial_value_t v = vx_sets_subscript(x, &av, in);
            if (v.kind == VIAL_VAL_SETS) {
                const T st = vx_sets_type(x, v.sets, id);
                t.flags |= st.flags & (VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY);
                set_value(x, id, v);
            } else if (kind_of(x, in) != VIAL_N_SLICE && value_of(x, in).kind == VIAL_VAL_INT) {
                const int64_t k = value_of(x, in).i, n = vial_sets_size(av.sets, 0);
                if (k >= n || k < -n) {
                    diag(x, "E0304", in, "position %lld is outside a set of %lld particles", (long long)k, (long long)n);
                    return err();
                }
            }
        }
        return t;
    }
    if (a.dtype == VIAL_DT_TUPLE) {
        int64_t k;
        if (in && kind_of(x, in) != VIAL_N_SLICE) infer(x, in);
        if (in && static_index(x, in, &k) && k >= -(int64_t)a.num_elems && k < (int64_t)a.num_elems) {
            const uint32_t at = (uint32_t)(k < 0 ? k + a.num_elems : k);
            // shape(a)[k]: the extent, when it is static
            const vial_value_t sv = value_of(x, on);
            if (sv.kind == VIAL_VAL_SHAPE && sv.type && at < sv.type->ndim && sv.type->dims[at].kind == VIAL_EXT_STATIC) {
                set_value(x, id, val_int(sv.type->dims[at].n));
            }
            return a.elems[at];
        }
        diag(x, "E0319", id, "a tuple is indexed by a constant position");
        return err();
    }
    // The index items, a tuple of them or one
    vial_node_id items[16];
    uint32_t num = 0;
    if (in && kind_of(x, in) == VIAL_N_TUPLE) {
        for (uint32_t i = 0; i < num_children(x, in) && num < 16; ++i) items[num++] = child(x, in, i);
    } else if (in) {
        items[num++] = in;
    }
    // How many existing dimensions the items take: all but None and ...
    uint32_t takes = 0;
    for (uint32_t i = 0; i < num; ++i) {
        const vial_node_kind_t k = kind_of(x, items[i]);
        if (k == VIAL_N_ELLIPSIS) continue;
        if (!(k == VIAL_N_CONST && sub_of(x, items[i]) == VIAL_KW_NONE_)) takes += 1;
    }
    if (takes > a.ndim) {
        diag(x, "E0319", id, "%u positions for `%s`, which has %u dimensions", takes, tstr(&a), a.ndim);
        return err();
    }
    T t = a;
    t.ndim = 0;
    uint8_t d = 0;      // the next dimension of a to index
    for (uint32_t i = 0; i < num; ++i) {
        const vial_node_id it = items[i];
        const vial_node_kind_t k = kind_of(x, it);
        if (k == VIAL_N_ELLIPSIS) {
            const uint32_t skip = a.ndim - takes;
            for (uint32_t s = 0; s < skip; ++s) vial_type_push_dim(&t, a.dims[d++]);
            continue;
        }
        if (k == VIAL_N_CONST && sub_of(x, it) == VIAL_KW_NONE_) {
            if (!vial_type_push_dim(&t, vial_extent_static(1))) goto too_many;
            continue;
        }
        const vial_extent_t e = a.dims[d++];
        if (k == VIAL_N_SLICE) {
            int64_t lo = 0, hi = 0, step = 1;
            const vial_node_id ln = child(x, it, 0), hn = child(x, it, 1), sn = child(x, it, 2);
            bool klo = true, khi = true, kst = true;
            if (ln) { T lt = infer(x, ln); if (vial_type_is_error(&lt)) return err(); klo = static_index(x, ln, &lo); }
            if (hn) { T ht = infer(x, hn); if (vial_type_is_error(&ht)) return err(); khi = static_index(x, hn, &hi); }
            if (sn) { T st = infer(x, sn); if (vial_type_is_error(&st)) return err(); kst = static_index(x, sn, &step); }
            if (kst && step == 0) { diag(x, "E0319", it, "a slice's step is not 0"); return err(); }
            vial_extent_t r = e;
            r.desc = e.desc;
            if (e.kind == VIAL_EXT_STATIC && klo && khi && kst) {
                // Python's slice arithmetic
                const int64_t n = e.n;
                int64_t b, en;
                if (step > 0) {
                    b = !ln ? 0 : (lo < 0 ? (lo + n < 0 ? 0 : lo + n) : (lo > n ? n : lo));
                    en = !hn ? n : (hi < 0 ? (hi + n < 0 ? 0 : hi + n) : (hi > n ? n : hi));
                    r.n = en > b ? (uint32_t)((en - b + step - 1) / step) : 0;
                } else {
                    b = !ln ? n - 1 : (lo < 0 ? (lo + n < -1 ? -1 : lo + n) : (lo >= n ? n - 1 : lo));
                    en = !hn ? -1 : (hi < 0 ? (hi + n < -1 ? -1 : hi + n) : (hi >= n ? n - 1 : hi));
                    r.n = b > en ? (uint32_t)((b - en - step - 1) / (-step)) : 0;
                }
                r.kind = VIAL_EXT_STATIC;
            } else if (e.kind == VIAL_EXT_DYNAMIC && kst && step == 1 && klo && khi) {
                // A slice by fixed offsets keeps a related identity: c[1:] and c[:-1] pair
                const int64_t off_lo = ln ? lo : 0;
                const int64_t off_hi = hn ? (hi < 0 ? -hi : -1) : 0;
                if (off_lo >= 0 && off_hi >= 0 && (!hn || hi < 0)) { r.offset = e.offset + (int32_t)(off_lo + off_hi); }
                else { r.kind = VIAL_EXT_DYNAMIC; r.id = (uint32_t)id; r.offset = 0; }
            } else {
                r = (vial_extent_t){.kind = e.kind == VIAL_EXT_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)id, .desc = e.desc};
                if (e.kind == VIAL_EXT_UNKNOWN) r.kind = VIAL_EXT_UNKNOWN;
            }
            if (!vial_type_push_dim(&t, r)) goto too_many;
            continue;
        }
        const T ix = infer(x, it);
        if (vial_type_is_error(&ix)) return err();
        if ((ix.dtype == VIAL_DT_INT) && ix.ndim == 0) {
            int64_t k2;
            if (static_index(x, it, &k2) && e.kind == VIAL_EXT_STATIC && (k2 >= (int64_t)e.n || k2 < -(int64_t)e.n)) {
                char ext[32];
                snprintf(ext, sizeof(ext), "[%u]", e.n);
                diag(x, "E0304", it, "position %lld is outside `%s`", (long long)k2, ext);
                return err();
            }
            if (e.kind == VIAL_EXT_DYNAMIC && ix.stage == VIAL_STAGE_STATIC) {
                // A constant position may be past the end in some frame; a position computed per frame, such as the
                // index closest_pair gives, is one the array has
                t.flags |= VIAL_TYPE_MISSING;
                diag(x, "W0405", it, "a position in `[?]`, whose length changes per frame, is missing in the frames where it is outside (R0101)");
            }
            join_stage(&t, &t, &ix);
            continue;   // the dimension is removed
        }
        if (ix.dtype == VIAL_DT_INT && ix.ndim == 1) {
            // Positions listed, in that order
            vial_extent_t r = ix.dims[0];
            r.desc = e.desc;
            if (!vial_type_push_dim(&t, r)) goto too_many;
            join_stage(&t, &t, &ix);
            continue;
        }
        if (ix.dtype == VIAL_DT_BOOL && ix.ndim >= 1) {
            // A mask over the leading dimensions: a shorter dimension, per frame if the mask is
            vial_extent_t r = {.kind = ix.stage == VIAL_STAGE_STATIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC, .id = (uint32_t)it, .desc = e.desc};
            d += (uint8_t)(ix.ndim - 1);
            if (!vial_type_push_dim(&t, r)) goto too_many;
            join_stage(&t, &t, &ix);
            continue;
        }
        if (ix.dtype == VIAL_DT_SET) {
            // a[s]: the values of the particles in s, an intersection, along a's first dimension over particles
            if (ix.ndim || num != 1) {
                diag(x, "E0319", it, "a set as a position is one set, alone in the brackets; this is %s", tstr(&ix));
                return err();
            }
            int pd = -1;
            for (uint8_t q = 0; q < a.ndim && pd < 0; ++q) {
                const uint32_t dsc = a.dims[q].desc;
                if (dsc && dsc < md_array_size(*x->descs) && (*x->descs)[dsc].level == VIAL_LEVEL_PARTICLE) pd = q;
            }
            if (pd < 0) {
                diag_fix(x, "E0311", on, NULL, "a set picks values by particle, and `%s` has no dimension over particles", tstr(&a));
                return err();
            }
            T r = a;
            join_stage(&r, &a, &ix);
            const vial_desc_t* dd = &(*x->descs)[a.dims[pd].desc];
            const vial_value_t sv = value_of(x, it);
            if (dd->sets && dd->sets->single && sv.kind == VIAL_VAL_SETS && sv.sets->single && ix.stage == VIAL_STAGE_STATIC &&
                a.dims[pd].kind == VIAL_EXT_STATIC) {
                const vial_sets_t* both = vial_sets_and(dd->sets, 0, sv.sets, 0, x->alloc);
                r.dims[pd] = vial_extent_static(vial_sets_size(both, 0));
                r.dims[pd].desc = vx_desc(x, VIAL_LEVEL_PARTICLE, both);
            } else {
                r.dims[pd] = (vial_extent_t){.kind = ix.stage == VIAL_STAGE_STATIC && a.dims[pd].kind != VIAL_EXT_DYNAMIC ? VIAL_EXT_UNKNOWN : VIAL_EXT_DYNAMIC,
                                             .id = (uint32_t)id, .desc = vx_desc(x, VIAL_LEVEL_PARTICLE, NULL)};
            }
            return r;
        }
        diag(x, "E0319", it, "a position is an int, a slice, a list of ints or a bool mask; this is %s", tstr(&ix));
        return err();
    }
    while (d < a.ndim) if (!vial_type_push_dim(&t, a.dims[d++])) goto too_many;
    if (num == 1 && t.stage == VIAL_STAGE_STATIC) {
        // The value of one position, slice or list of an array of sets or of strings
        const vial_value_t av = value_of(x, on);
        if (av.kind == VIAL_VAL_SETS || av.kind == VIAL_VAL_STRS || av.kind == VIAL_VAL_LABELS) {
            const vial_value_t v = vx_sets_subscript(x, &av, items[0]);
            if (v.kind == VIAL_VAL_SETS) {
                const T st = vx_sets_type(x, v.sets, id);
                t.flags = (t.flags & ~(VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY)) | (st.flags & (VIAL_TYPE_ONE_PARTICLE | VIAL_TYPE_EMPTY));
                if (t.ndim >= 1 && st.ndim >= 1) t.dims[0] = st.dims[0];
                set_value(x, id, v);
            } else if (v.kind != VIAL_VAL_NONE) {
                set_value(x, id, v);
            }
        }
    }
    return t;
too_many:
    diag(x, "E0308", id, "the result would have 7 dimensions; the limit is 6");
    return err();
}

static T infer_attr(infer_t* x, vial_node_id id) {
    const vial_node_id on = child(x, id, 0), nn = child(x, id, 1);
    const T a = infer(x, on);
    if (vial_type_is_error(&a) || !nn) return err();
    const str_t name = text(x, nn);
    if (str_eq(name, STR_LIT("T"))) {
        if (a.dtype == VIAL_DT_TUPLE || a.dtype == VIAL_DT_DEF) { diag(x, "E0319", id, "`.T` takes an array; this is %s", dtype_word(&a)); return err(); }
        T t = a;
        for (uint8_t k = 0; k < a.ndim; ++k) t.dims[k] = a.dims[a.ndim - 1 - k];
        return t;
    }
    if (a.dtype == VIAL_DT_SET || str_eq(name, STR_LIT("com"))) {
        T t;
        if (vx_attr_geom(x, id, on, &a, name, &t)) return t;
        return err();
    }
    const bool comp = str_eq(name, STR_LIT("x")) || str_eq(name, STR_LIT("y")) || str_eq(name, STR_LIT("z"));
    if (comp || str_eq(name, STR_LIT("xyz"))) {
        const bool coords = is_float(a.dtype) && a.ndim >= 1 && a.dims[a.ndim - 1].kind == VIAL_EXT_STATIC && a.dims[a.ndim - 1].n == 3;
        if (!coords) {
            diag(x, "E0319", id, "`.%.*s` takes coordinates, a last dimension of 3, or a set; this is `%s`", (int)name.len, name.ptr, tstr(&a));
            return err();
        }
        T t = a;
        if (comp) t.ndim -= 1;
        return t;
    }
    return err();   // functions as methods: with the built-ins and defs
}

static T infer_uncached(infer_t* x, vial_node_id id) {
    switch (kind_of(x, id)) {
    case VIAL_N_INT: {
        // VIAL's ints are 32 bits: a literal beyond them is typed but not folded
        const int64_t i = parse_int(text(x, id));
        set_value(x, id, i >= INT32_MIN && i <= INT32_MAX ? val_int(i) : val_none());
        return number(VIAL_DT_INT);
    }
    case VIAL_N_FLOAT:    set_value(x, id, val_float(parse_float(text(x, id)))); return number(VIAL_DT_FLOAT);
    case VIAL_N_QUANTITY: return infer_quantity(x, id);
    case VIAL_N_STRING:   set_value(x, id, (vial_value_t){.kind = VIAL_VAL_STR, .s = string_value(x, id)}); return scalar(VIAL_DT_STR);
    case VIAL_N_CONST:
        switch (sub_of(x, id)) {
        case VIAL_KW_TRUE:  set_value(x, id, val_bool(true));  return scalar(VIAL_DT_BOOL);
        case VIAL_KW_FALSE: set_value(x, id, val_bool(false)); return scalar(VIAL_DT_BOOL);
        default:            set_value(x, id, (vial_value_t){.kind = VIAL_VAL_NULL}); return scalar(VIAL_DT_NONE);
        }
    case VIAL_N_NAME:     return infer_name(x, id);
    case VIAL_N_UNARY:    return infer_unary(x, id);
    case VIAL_N_NOT:      return infer_not(x, id);
    case VIAL_N_BINARY:   return infer_binary(x, id);
    case VIAL_N_BOOLOP:   return infer_boolop(x, id);
    case VIAL_N_COMPARE:  return infer_compare(x, id);
    case VIAL_N_RANGE:    return infer_range(x, id);
    case VIAL_N_TERNARY:  return infer_ternary(x, id);
    case VIAL_N_LIST:     return infer_list(x, id);
    case VIAL_N_TUPLE:    return infer_tuple(x, id);
    case VIAL_N_SUBSCRIPT: return infer_subscript(x, id);
    case VIAL_N_ATTR:     return infer_attr(x, id);
    case VIAL_N_LISTCOMP: return vx_listcomp(x, id);
    case VIAL_N_CALL: {
        // Built-ins typed so far; everything else still has its arguments checked
        const vial_node_id callee = child(x, id, 0);
        if (callee && kind_of(x, callee) == VIAL_N_NAME) {
            const vial_symbol_t* sym = vial_names_symbol(x->n, callee);
            T t;
            if (sym && sym->kind == VIAL_SYM_BUILTIN && sym->builtin && vx_call_sets(x, id, sym->builtin, &t)) return t;
            if (sym && sym->kind == VIAL_SYM_BUILTIN && sym->builtin && vx_call_geom(x, id, sym->builtin, &t)) return t;
            if (sym && sym->kind == VIAL_SYM_BUILTIN && sym->builtin && vx_call_math(x, id, sym->builtin, &t)) return t;
            if (sym && sym->kind == VIAL_SYM_BUILTIN && sym->builtin && vx_call_draw(x, id, sym->builtin, &t)) return t;
            if (sym && sym->kind == VIAL_SYM_BUILTIN && sym->builtin && vx_call_frames(x, id, sym->builtin, &t)) return t;
            if (sym && sym->kind == VIAL_SYM_PRELUDE && x->c->prelude) {
                vial_node_id d = sym->node;
                if (d && x->c->prelude->nodes[d].kind == VIAL_N_NAME) d = x->c->prelude_names.parent[d];
                if (d && x->c->prelude->nodes[d].kind == VIAL_N_DEF) {
                    return vx_call_def(x, id, x->c->prelude, &x->c->prelude_names, x->c->prelude_src, true, d);
                }
            }
            if (sym && sym->kind == VIAL_SYM_DEF && sym->node) {
                vial_node_id d = sym->node;
                if (kind_of(x, d) == VIAL_N_NAME) d = x->n->parent[d];
                if (d && kind_of(x, d) == VIAL_N_DEF) return vx_call_def(x, id, x->t, x->n, x->src, x->t == x->c->prelude, d);
            }
        }
    }
    // fallthrough
    default:
        // Calls, comprehensions and what arrives later: their parts are still checked
        for (uint32_t i = 0; i < num_children(x, id); ++i) {
            const vial_node_id c = child(x, id, i);
            if (!c) continue;
            const vial_node_kind_t k = kind_of(x, c);
            if (k == VIAL_N_KWARG) { if (child(x, c, 1)) infer(x, child(x, c, 1)); continue; }
            if (kind_of(x, id) == VIAL_N_CALL && i == 0) continue;
            if (k == VIAL_N_COMP_FOR || kind_of(x, id) == VIAL_N_LISTCOMP) continue;    // comprehensions: with defs
            if (k == VIAL_N_DEF || k == VIAL_N_UNSUPPORTED || k == VIAL_N_UNIT || k == VIAL_N_ERROR) continue;
            infer(x, c);
        }
        return err();
    }
}

T vx_infer(infer_t* x, vial_node_id id) { return infer(x, id); }

static T infer(infer_t* x, vial_node_id id) {
    if (!id) return err();
    if (x->out->typed[(id) - x->out->base]) return x->out->types[(id) - x->out->base];
    T t = infer_uncached(x, id);
    if (t.ndim > VIAL_MAX_DIMS) t = err();
    x->out->types[(id) - x->out->base] = t;
    x->out->typed[(id) - x->out->base] = 1;
    if (x->ir) vx_emit(x, id, &t);
    return t;
}

// ### STATEMENTS ###

static void assign_target(infer_t* x, vial_node_id target, const T* t, vial_node_id value) {
    const vial_node_kind_t k = kind_of(x, target);
    if (k == VIAL_N_NAME) {
        const vial_value_t v = value ? x->out->values[(value) - x->out->base] : val_none();
        x->out->types[(target) - x->out->base] = *t;
        x->out->typed[(target) - x->out->base] = 1;
        x->out->values[(target) - x->out->base] = v;
        const uint32_t s = x->n->ref[target];
        if (s && x->frame) {
            vx_frame_set(x, s, t, &v, provenance(x, value));
        } else if (s) {
            x->out->sym_types[s] = *t;
            x->out->sym_values[s] = v;
        }
        if (x->ir && s) {
            const vial_ssa_t sv = value ? vx_ssa(x, value) : 0;
            vx_frame_set_ssa(x, s, sv);
            vx_set_ssa(x, target, sv);
        }
        return;
    }
    if (k == VIAL_N_TUPLE || k == VIAL_N_LIST) {
        const uint32_t n = num_children(x, target);
        const vial_node_id vt = value && (kind_of(x, value) == VIAL_N_TUPLE || kind_of(x, value) == VIAL_N_LIST) &&
                                num_children(x, value) == n ? value : 0;
        for (uint32_t i = 0; i < n; ++i) {
            T e = err();
            if (!vial_type_is_error(t)) {
                if (t->dtype == VIAL_DT_TUPLE && t->num_elems == n) e = t->elems[i];
                else if (t->dtype != VIAL_DT_TUPLE && t->ndim >= 1 && t->dims[0].kind == VIAL_EXT_STATIC && t->dims[0].n == n) {
                    e = *t;
                    e.ndim -= 1;
                    memmove(e.dims, e.dims + 1, e.ndim * sizeof(e.dims[0]));
                } else if (i == 0) {
                    diag(x, "E0319", target, "%u names, and `%s` does not unpack into %u", n, tstr(t), n);
                }
            }
            if (x->ir && !vt && value && !vial_type_is_error(&e) && vx_ssa(x, value)) {
                // One element of the value: its i-th, of a tuple or along the first dimension
                const vial_node_id tn = child(x, target, i);
                const uint32_t sym = kind_of(x, tn) == VIAL_N_NAME ? x->n->ref[tn] : 0;
                assign_target(x, tn, &e, 0);
                if (sym) {
                    const vial_ssa_t sv = vx_unpack(x, target, vx_ssa(x, value), t, i, &e);
                    vx_frame_set_ssa(x, sym, sv);
                    vx_set_ssa(x, tn, sv);
                }
                continue;
            }
            assign_target(x, child(x, target, i), &e, vt ? child(x, vt, i) : 0);
        }
    }
}

static bool check_stmt(infer_t* x, vial_node_id s);

// The statements of a block in order; true when it surely returned, so what follows is never reached
static bool check_block(infer_t* x, vial_node_id b) {
    if (!b) return false;
    if (kind_of(x, b) != VIAL_N_BLOCK) return check_stmt(x, b);
    for (uint32_t i = 0; i < num_children(x, b); ++i) {
        if (check_stmt(x, child(x, b, i))) return true;
    }
    return false;
}

// An if: a static condition checks only the branch it takes. True when every branch that can run returned.
static bool check_if(infer_t* x, vial_node_id s);
static uint8_t ct_stage_of(const vial_ir_t* ir, vial_ssa_t v) { return ir->values[v].type.stage; }

// Whether a statement holds a return, outside the defs nested in it
static bool contains_return(const infer_t* x, vial_node_id n) {
    if (!n) return false;
    const vial_node_kind_t k = kind_of(x, n);
    if (k == VIAL_N_RETURN) return true;
    if (k == VIAL_N_DEF) return false;
    for (uint32_t i = 0; i < num_children(x, n); ++i) if (contains_return(x, child(x, n, i))) return true;
    return false;
}

static uint32_t* env_copy(const infer_t* x) {
    uint32_t* e = NULL;
    const size_t n = md_array_size(x->frame->sym_ssa);
    md_array_resize(e, n, x->alloc);
    if (n) memcpy(e, x->frame->sym_ssa, n * sizeof(uint32_t));
    return e;
}

static void env_set(infer_t* x, const uint32_t* e) {
    vx_frame_t* f = x->frame;
    const size_t n = md_array_size(f->sym_ssa), m = md_array_size(e);
    for (size_t i = 0; i < n; ++i) f->sym_ssa[i] = i < m ? e[i] : 0;
}

static uint32_t env_at(const uint32_t* e, size_t i) { return i < md_array_size(e) ? e[i] : 0; }

// An if in a def whose condition is known only per frame: an IF op whose regions are the branches. Both branches
// return, and the if's value is the def's; or neither holds a return, and the names they bind are merged. Anything
// else is not lowered yet: false, with nothing emitted that is reachable.
static bool lower_dynamic_if(infer_t* x, vial_node_id s, bool* returned) {
    vial_ir_t* ir = x->ir;
    const vial_node_id c = child(x, s, 0), then_n = child(x, s, 1), other = child(x, s, 2);
    const vial_ssa_t cv = vx_ssa(x, c);
    if (!cv) return false;
    const bool rt = contains_return(x, then_n), re = contains_return(x, other);
    const uint32_t saved = ir->insert, func = ir->regions[ir->insert].func;
    uint32_t* snap = env_copy(x);
    const uint32_t r1 = vial_ir_region(ir, VIAL_IR_NO_OP, func, NULL, 0);
    ir->insert = r1;
    const bool a = check_block(x, then_n);
    uint32_t* env1 = env_copy(x);
    env_set(x, snap);
    const uint32_t r2 = vial_ir_region(ir, VIAL_IR_NO_OP, func, NULL, 0);
    ir->insert = r2;
    const bool b = other ? (kind_of(x, other) == VIAL_N_IF ? check_if(x, other) : check_block(x, other)) : false;
    uint32_t* env2 = env_copy(x);
    ir->insert = saved;
    const md_vial_span_t sp = vx_span(x, s);
    if (a && b) {
        // Both return: each region ended in the yield of its return
        const T rt_t = x->frame->has_result ? x->frame->result : scalar(VIAL_DT_NONE);
        const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_IF, .args = &cv, .num_args = 1, .types = &rt_t, .num_results = 1, .span = sp});
        vial_ir_attach(ir, op, 0, r1);
        vial_ir_attach(ir, op, 1, r2);
        const vial_ssa_t r = ir->ops[op].result;
        vial_ir_yield(ir, &r, 1, sp);
        *returned = true;
        return true;
    }
    if (rt || re) {
        env_set(x, env2);
        vx_lower_fail(x, s, "an if per frame where one branch returns and the other does not");
        *returned = a && b;
        return true;
    }
    // Neither returns: the names either branch rebinds, from the branch taken
    vial_ssa_t y1[64], y2[64];
    T types[64];
    uint32_t slots[64], n = 0;
    const size_t ns = md_array_size(x->frame->sym_ssa);
    bool ok = true;
    for (size_t i = 0; i < ns && ok; ++i) {
        const uint32_t s0 = env_at(snap, i), s1 = env_at(env1, i), s2 = env_at(env2, i);
        if (s1 == s0 && s2 == s0) continue;
        if (!s1 || !s2 || n == 64) { ok = false; break; }
        const T* t1 = vial_ir_type(ir, s1);
        const T* t2 = vial_ir_type(ir, s2);
        if (t1->dtype != t2->dtype || t1->ndim != t2->ndim) { ok = false; break; }
        y1[n] = s1; y2[n] = s2; types[n] = *t1; slots[n] = (uint32_t)i;
        if (t2->stage > types[n].stage) types[n].stage = t2->stage;
        if (ct_stage_of(ir, cv) > types[n].stage) types[n].stage = ct_stage_of(ir, cv);
        n++;
    }
    if (!ok) {
        env_set(x, env2);
        vx_lower_fail(x, s, "an if per frame whose branches bind a name differently");
        *returned = false;
        return true;
    }
    ir->insert = r1;
    vial_ir_yield(ir, y1, n, sp);
    ir->insert = r2;
    vial_ir_yield(ir, y2, n, sp);
    ir->insert = saved;
    const uint32_t op = vial_ir_emit(ir, &(vial_ir_emit_t){.kind = VIAL_IR_IF, .args = &cv, .num_args = 1, .types = types, .num_results = n, .span = sp});
    vial_ir_attach(ir, op, 0, r1);
    vial_ir_attach(ir, op, 1, r2);
    env_set(x, snap);
    for (uint32_t k = 0; k < n; ++k) x->frame->sym_ssa[slots[k]] = ir->ops[op].result + k;
    *returned = false;
    return true;
}

static bool check_if(infer_t* x, vial_node_id s) {
    const vial_node_id c = child(x, s, 0);
    const T ct = infer(x, c);
    const vial_value_t v = value_of(x, c);
    const vial_node_id other = child(x, s, 2);
    if (!vial_type_is_error(&ct) && v.kind == VIAL_VAL_BOOL && ct.stage == VIAL_STAGE_STATIC) {
        if (v.b) return check_block(x, child(x, s, 1));
        return other ? (kind_of(x, other) == VIAL_N_IF ? check_if(x, other) : check_block(x, other)) : false;
    }
    if (!vial_type_is_error(&ct) && (ct.dtype != VIAL_DT_BOOL || ct.ndim)) {
        diag_fix(x, "E0303", c, "where(...)", "an if needs one bool; this is `%s`", tstr(&ct));
    } else if (!x->frame && !vial_type_is_error(&ct) && ct.stage != VIAL_STAGE_STATIC) {
        diag_fix(x, "E0403", c, "where(...)", "the condition changes per frame, and what a script publishes cannot; where() chooses per frame");
    }
    // Lowering: in a def, an if whose branch is chosen per frame is an if of the IR, its branches regions
    bool returned;
    if (x->ir && x->frame && !vial_type_is_error(&ct) && ct.dtype == VIAL_DT_BOOL && !ct.ndim && lower_dynamic_if(x, s, &returned)) return returned;
    // Otherwise its branches are checked without emitting, and the def it is in is not evaluated
    struct vial_ir_t* ir = x->ir;
    if (ir) {
        // An error was reported for the condition, E0303 or E0403, unless it is one bool in a def
        const bool reported = vial_type_is_error(&ct) || ct.dtype != VIAL_DT_BOOL || ct.ndim || !x->frame;
        if (reported) vx_emit_unreachable_mark(x);
        else vx_lower_fail(x, s, "an if whose condition is not known when compiling");
        x->ir = NULL;
    }
    const bool a = check_block(x, child(x, s, 1));
    const bool b = other ? (kind_of(x, other) == VIAL_N_IF ? check_if(x, other) : check_block(x, other)) : false;
    x->ir = ir;
    return a && b;
}

bool vx_check_block(infer_t* x, vial_node_id b) { return check_block(x, b); }

void vx_frame_return(infer_t* x, vial_node_id at, const T* t) {
    vx_frame_t* f = x->frame;
    if (vial_type_is_error(t)) { f->result = err(); f->has_result = true; return; }
    if (!f->has_result) { f->result = *t; f->has_result = true; f->value = value_of(x, at); return; }
    f->value = val_none();
    if (vial_type_is_error(&f->result)) return;
    T m;
    if (f->result.dtype == VIAL_DT_TUPLE || t->dtype == VIAL_DT_TUPLE) {
        if (f->result.dtype != t->dtype || f->result.num_elems != t->num_elems) {
            diag(x, "E0319", at, "the returns of a def have one type; these are %s and %s", tstr(&f->result), tstr(t));
            f->result = err();
        }
        return;
    }
    if (!join_types(x, at, "the returns of a def", &f->result, t, &m)) { f->result = err(); return; }
    f->result = m;
}

// What a top-level name publishes must keep its shape from frame to frame (E0401), selections apart, and should not
// take more than 64 MB over the run (W0404)
static void check_published(infer_t* x, vial_node_id target, const T* t) {
    if (x->frame || kind_of(x, target) != VIAL_N_NAME || vial_type_is_error(t)) return;
    const str_t name = text(x, target);
    if (name.len && name.ptr[0] == '_') return;
    const bool publishable = t->dtype == VIAL_DT_BOOL || t->dtype == VIAL_DT_INT || t->dtype == VIAL_DT_FLOAT ||
                             t->dtype == VIAL_DT_FLOAT64 || t->dtype == VIAL_DT_SET;
    if (!publishable || t->dtype == VIAL_DT_SET) return;
    bool dynamic = false;
    uint64_t n = 1;
    for (uint8_t k = 0; k < t->ndim; ++k) {
        if (t->dims[k].kind == VIAL_EXT_DYNAMIC) dynamic = true;
        if (t->dims[k].kind == VIAL_EXT_STATIC) n *= t->dims[k].n;
    }
    char fix[80];
    snprintf(fix, sizeof(fix), "_%.*s", (int)(name.len < 70 ? name.len : 70), name.ptr);
    if (dynamic) {
        diag_fix(x, "E0401", target, fix, "`%.*s` has a different length in each frame, so it cannot be published as a series; reduce it, or name it %s",
            (int)name.len, name.ptr, fix);
        return;
    }
    const bool per_frame = t->stage == VIAL_STAGE_FRAME || (t->flags & VIAL_TYPE_PER_FRAME);
    if (per_frame && x->c->num_frames && t->ndim < 3 && n * x->c->num_frames > 16u * 1024u * 1024u) {
        const double mb = (double)(n * x->c->num_frames) * 4.0 / (1024.0 * 1024.0);
        diag_fix(x, "W0404", target, fix, "`%.*s` takes %.0f MB over %llu frames; published values are kept for the whole run",
            (int)name.len, name.ptr, mb, (unsigned long long)x->c->num_frames);
    }
}

// W0405, once per top-level statement whose value may be missing, unless it reads a name that already may be, or
// the statement was already told
static void warn_missing(infer_t* x, vial_node_id s, vial_node_id value, const T* t) {
    if (x->frame || !(t->flags & VIAL_TYPE_MISSING) || !value) return;
    const md_vial_span_t span = x->t->nodes[s].span;
    for (size_t i = md_array_size(x->diags->items); i-- > 0;) {
        const md_vial_diagnostic_t* d = &x->diags->items[i];
        if (d->span.beg >= span.beg && d->span.end <= span.end && strncmp(d->code, "W0405", 5) == 0) return;
        if (d->span.end < span.beg) break;
    }
    vial_node_id* stack = NULL;
    md_array_push(stack, value, x->alloc);
    bool inherited = false;
    while (md_array_size(stack) && !inherited) {
        const vial_node_id n = *md_array_last(stack);
        md_array_pop(stack);
        if (kind_of(x, n) == VIAL_N_NAME) {
            const uint32_t sym = x->n->ref[n];
            const vial_symbol_t* sy = vial_names_symbol(x->n, n);
            if (sy && sy->kind == VIAL_SYM_GLOBAL && (x->top->sym_types[sym].flags & VIAL_TYPE_MISSING)) inherited = true;
        }
        for (uint32_t i = 0; i < num_children(x, n); ++i) if (child(x, n, i)) md_array_push(stack, child(x, n, i), x->alloc);
    }
    if (!inherited) diag(x, "W0405", value, "this may be missing in some frames; the series shows gaps there");
}

static bool check_stmt(infer_t* x, vial_node_id s) {
    switch (kind_of(x, s)) {
    case VIAL_N_ASSIGN: {
        const uint32_t n = num_children(x, s);
        if (!n) return false;
        const vial_node_id value = child(x, s, n - 1);
        const T t = infer(x, value);
        for (uint32_t i = 0; i + 1 < n; ++i) {
            assign_target(x, child(x, s, i), &t, value);
            check_published(x, child(x, s, i), &t);
        }
        warn_missing(x, s, value, &t);
        if (!x->frame && !vial_type_is_error(&t) && vx_pass(&t) > 3) {
            diag(x, "W0406", value, "this is pass %u; each pass reads every frame again", vx_pass(&t));
        }
        return false;
    }
    case VIAL_N_EXPR_STMT:
        if (child(x, s, 0)) {
            const T t = infer(x, child(x, s, 0));
            warn_missing(x, s, child(x, s, 0), &t);
        }
        return false;
    case VIAL_N_IF:
        return check_if(x, s);
    case VIAL_N_FOR:
        return vx_check_for(x, s);
    case VIAL_N_RETURN: {
        if (!x->frame) return true;
        const vial_node_id v = child(x, s, 0);
        const T t = v ? infer(x, v) : scalar(VIAL_DT_NONE);
        if (v) set_value(x, s, value_of(x, v));
        vx_frame_return(x, v ? v : s, &t);
        if (x->ir) vx_emit_return(x, v, s);
        return true;
    }
    case VIAL_N_ASSERT: {
        const vial_node_id c = child(x, s, 0), m = child(x, s, 1);
        const T ct = infer(x, c);
        if (m) infer(x, m);
        const vial_value_t v = value_of(x, c);
        if (!vial_type_is_error(&ct) && ct.stage == VIAL_STAGE_STATIC && v.kind == VIAL_VAL_BOOL && !v.b) {
            const vial_value_t mv = m ? value_of(x, m) : val_none();
            if (mv.kind == VIAL_VAL_STR) {
                char ex[256];
                const size_t len = mv.s.len < sizeof(ex) - 1 ? mv.s.len : sizeof(ex) - 1;
                memcpy(ex, mv.s.ptr, len);
                ex[len] = '\0';
                for (size_t i = 0; i < len; ++i) if (ex[i] == '\n' || ex[i] == '\r') ex[i] = ' ';
                diag(x, "E0314", c, "%s", ex);
            } else {
                diag(x, "E0314", c, "the assert fails for these arguments");
            }
            if (x->ir) vx_emit_unreachable(x);
            return true;    // nothing after it runs
        }
        return false;
    }
    default:
        return false;   // nested defs are checked at their calls; loops and the rest arrive with later steps
    }
}

void vial_check_types(vial_check_t* c, vial_types_t* out) {
    const size_t num = md_array_size(c->tree->nodes);
    const size_t nsym = md_array_size(c->names.symbols);
    memset(out, 0, sizeof(*out));
    md_array_resize(out->types, num, c->alloc);
    md_array_resize(out->values, num, c->alloc);
    md_array_resize(out->typed, num, c->alloc);
    md_array_resize(out->sym_types, nsym, c->alloc);
    md_array_resize(out->sym_values, nsym, c->alloc);
    memset(out->values, 0, num * sizeof(*out->values));
    memset(out->typed, 0, num);
    for (size_t i = 0; i < nsym; ++i) out->sym_types[i] = vial_type_error();
    memset(out->sym_values, 0, nsym * sizeof(*out->sym_values));
    for (size_t i = 0; i < num; ++i) out->types[i] = vial_type_error();

    vx_memo_t* memo = md_alloc(c->alloc, sizeof(vx_memo_t));
    memset(memo, 0, sizeof(*memo));
    memo->map.allocator = c->alloc;
    infer_t x = {.t = c->tree, .src = c->src, .diags = c->diags, .n = &c->names, .out = out, .top = out, .descs = &out->descs, .c = c, .memo = memo,
        .alloc = c->alloc};
    vial_sys_t* sys = md_alloc(c->alloc, sizeof(vial_sys_t));
    if (c->sys && vial_sys_init(sys, c->sys, c->alloc)) x.sys = sys;
    if (c->ir) {
        x.ir = c->ir;
        vx_sizes_ssa(&x, out, num, nsym);
    }
    const vial_node_id root = c->tree->root;
    for (uint32_t i = 0; i < num_children(&x, root); ++i) check_stmt(&x, child(&x, root, i));
    if (c->ir) vx_emit_finish(&x);
    vial_diags_sort(c->diags);
}

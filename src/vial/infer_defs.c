// The checker, step 7: defs. A def's result is inferred at each call: its arguments are bound to its parameters,
// split by their kinds into what one call takes and what it maps over, and its body is checked for those types, in a
// frame of its own. A static condition checks only the branch it takes, so `if between is None:` in distance checks
// one return. The result is memoized by the def and its arguments' types and static values.
//
// The prelude's defs are checked the same way. What they report is reported in the script: at the argument a
// parameter stands for, or at the call that led into the prelude.

#include "infer_internal.h"

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_hash.h>

#include <stdio.h>
#include <string.h>

typedef vial_type_t T;

static T err(void) { return vial_type_error(); }

static inline vial_node_kind_t kind_in(const vial_tree_t* t, vial_node_id id) { return (vial_node_kind_t)t->nodes[id].kind; }
static inline str_t text_in(const vial_tree_t* t, str_t src, vial_node_id id) {
    const md_vial_span_t s = t->nodes[id].span;
    return (str_t){src.ptr + s.beg, s.end - s.beg};
}

// ### FRAMES ###

static const uint32_t* slot_of(const vx_frame_t* f, uint32_t sym) {
    return f->syms.num_buckets ? md_hashmap_get((md_hashmap32_t*)&f->syms, (uint64_t)sym) : NULL;
}

bool vx_frame_get(const infer_t* x, uint32_t sym, T* t, vial_value_t* v) {
    for (const vx_frame_t* f = x->frame; f; f = f->parent) {
        if (f->t != x->t) continue;
        const uint32_t* s = slot_of(f, sym);
        if (s) { *t = f->sym_types[*s]; *v = f->sym_values[*s]; return true; }
    }
    return false;
}

vial_node_id vx_frame_prov(const infer_t* x, uint32_t sym) {
    for (const vx_frame_t* f = x->frame; f; f = f->parent) {
        if (f->t != x->t) continue;
        const uint32_t* s = slot_of(f, sym);
        if (s) return f->sym_prov[*s];
    }
    return 0;
}

void vx_frame_set(infer_t* x, uint32_t sym, const T* t, const vial_value_t* v, vial_node_id prov) {
    vx_frame_t* f = x->frame;
    const uint32_t* s = slot_of(f, sym);
    if (s) {
        f->sym_types[*s] = *t;
        f->sym_values[*s] = *v;
        f->sym_prov[*s] = prov;
        return;
    }
    const uint32_t slot = (uint32_t)md_array_size(f->sym_types);
    md_array_push(f->sym_types, *t, x->alloc);
    md_array_push(f->sym_values, *v, x->alloc);
    md_array_push(f->sym_prov, prov, x->alloc);
    md_hashmap_add(&f->syms, (uint64_t)sym, slot);
}

// ### MEMO ###

static uint64_t hash_type(uint64_t h, const T* t) {
    uint8_t buf[16 + VIAL_MAX_DIMS * 24];
    size_t n = 0;
    buf[n++] = t->dtype; buf[n++] = t->ndim; buf[n++] = t->stage; buf[n++] = t->pass; buf[n++] = t->unit_kind;
    memcpy(buf + n, &t->flags, 2); n += 2;
    for (uint8_t k = 0; k < t->ndim && k < VIAL_MAX_DIMS; ++k) {
        const vial_extent_t* e = &t->dims[k];
        buf[n++] = e->kind;
        memcpy(buf + n, &e->n, 4); n += 4;
        memcpy(buf + n, &e->m, 4); n += 4;
        memcpy(buf + n, &e->id, 4); n += 4;
        memcpy(buf + n, &e->offset, 4); n += 4;
        memcpy(buf + n, &e->desc, 4); n += 4;
    }
    h = md_hash64(buf, n, h);
    h = md_hash64(&t->unit.base, sizeof(t->unit.base), h);
    h = md_hash64(&t->unit.mult, sizeof(t->unit.mult), h);
    if (t->flags & VIAL_TYPE_PERIODIC) {
        const float p[3] = {t->period_beg, t->period_end, t->seam};
        h = md_hash64(p, sizeof(p), h);
    }
    h = md_hash64(&t->ident, sizeof(t->ident), h);
    for (uint32_t i = 0; i < t->num_elems && t->elems; ++i) h = hash_type(h, &t->elems[i]);
    return h;
}

static uint64_t hash_value(uint64_t h, const vial_value_t* v) {
    h = md_hash64(&v->kind, 1, h);
    switch (v->kind) {
    case VIAL_VAL_BOOL:  return md_hash64(&v->b, 1, h);
    case VIAL_VAL_INT:   return md_hash64(&v->i, 8, h);
    case VIAL_VAL_FLOAT: return md_hash64(&v->f, 8, h);
    case VIAL_VAL_STR:   return md_hash64_str(v->s, h);
    case VIAL_VAL_SETS:
    case VIAL_VAL_LABELS: return md_hash64(&v->sets, sizeof(v->sets), h);
    case VIAL_VAL_SHAPE: return v->type ? hash_type(h, v->type) : h;
    default: return h;
    }
}

// ### CALLS ###

static uint8_t kind_of_annotation(str_t a) {
    if (str_eq(a, STR_LIT("set")))    return VX_KIND_SET;
    if (str_eq(a, STR_LIT("group")))  return VX_KIND_GROUP;
    if (str_eq(a, STR_LIT("point")))  return VX_KIND_POINT;
    if (str_eq(a, STR_LIT("points"))) return VX_KIND_POINTS;
    if (str_eq(a, STR_LIT("float")))  return VX_KIND_FLOAT;
    if (str_eq(a, STR_LIT("int")))    return VX_KIND_INT;
    if (str_eq(a, STR_LIT("bool")))   return VX_KIND_BOOL;
    return VX_KIND_WHOLE;
}

// The smallest node id in a def's subtree: its nodes are made before it, so they lie in [min, def]
static vial_node_id subtree_min(const vial_tree_t* t, vial_node_id root, md_allocator_i* alloc) {
    vial_node_id lo = root;
    vial_node_id* stack = NULL;
    md_array_push(stack, root, alloc);
    while (md_array_size(stack)) {
        const vial_node_id n = *md_array_last(stack);
        md_array_pop(stack);
        if (n < lo) lo = n;
        for (uint32_t i = 0; i < t->nodes[n].child_count; ++i) {
            const vial_node_id c = vial_child(t, n, i);
            if (c) md_array_push(stack, c, alloc);
        }
    }
    md_array_free(stack, alloc);
    return lo;
}

// What a parameter's argument was written as in the script, for diagnostics in the prelude
static vial_node_id provenance_of(const infer_t* x, vial_node_id arg) {
    if (!arg) return 0;
    if (!x->frame || !x->frame->prelude) return x->t == x->c->tree ? arg : 0;
    if (x->t->nodes[arg].kind != VIAL_N_NAME) return 0;
    const uint32_t sym = x->n->ref[arg];
    return sym ? vx_frame_prov(x, sym) : 0;
}

static void infer_all_args(infer_t* x, vial_node_id call) {
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) { if (vx_child(x, a, 1)) vx_infer(x, vx_child(x, a, 1)); }
        else vx_infer(x, a);
    }
}

// The mapped dimensions, then the core's; a tuple gets them on each element
static bool with_loop(infer_t* x, vial_node_id call, const T* loop, const T* core, T* out) {
    if (core->dtype == VIAL_DT_TUPLE) {
        T* elems = md_alloc(x->alloc, (core->num_elems + 1) * sizeof(T));
        for (uint32_t i = 0; i < core->num_elems; ++i) if (!with_loop(x, call, loop, &core->elems[i], &elems[i])) return false;
        *out = *core;
        out->elems = elems;
        return true;
    }
    T t = *core;
    t.ndim = 0;
    for (uint8_t k = 0; k < loop->ndim; ++k) vial_type_push_dim(&t, loop->dims[k]);
    for (uint8_t k = 0; k < core->ndim; ++k) {
        if (!vial_type_push_dim(&t, core->dims[k])) {
            vx_diag_fix(x, "E0308", call, NULL, "the result would have 7 dimensions; the limit is 6");
            return false;
        }
    }
    if (loop->stage > t.stage) t.stage = loop->stage;
    *out = t;
    return true;
}

T vx_call_def(infer_t* x, vial_node_id call, const vial_tree_t* t, const vial_names_t* n, str_t src, bool prelude, vial_node_id def) {
    // A cycle of calls is E0802, from name resolution; checking it would not end
    const uint32_t depth = x->frame ? x->frame->depth + 1 : 1;
    bool cycle = depth > 32;
    for (const vx_frame_t* f = x->frame; f && !cycle; f = f->parent) cycle = f->def == def && f->t == t;
    if (cycle) { infer_all_args(x, call); return err(); }

    const vial_node_id name_node = vial_child(t, def, 0);
    const str_t name = name_node ? text_in(t, src, name_node) : STR_LIT("def");
    char fname[48];
    snprintf(fname, sizeof(fname), "%.*s", (int)(name.len < 40 ? name.len : 40), name.ptr);
    const vial_node_id params = vial_child(t, def, 1);
    const uint32_t np = params ? t->nodes[params].child_count : 0;
    if (np > 32) { infer_all_args(x, call); return err(); }
    vial_node_id pnodes[32], args[32];
    str_t pnames[32];
    for (uint32_t k = 0; k < np; ++k) {
        pnodes[k] = vial_child(t, params, k);
        const vial_node_id pn = pnodes[k] ? vial_child(t, pnodes[k], 0) : 0;
        pnames[k] = pn ? text_in(t, src, pn) : STR_LIT("");
        args[k] = 0;
    }

    // Bind: by position, then by keyword
    bool ok = true;
    uint32_t pos = 0;
    for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
        const vial_node_id a = vx_child(x, call, i);
        if (!a) continue;
        if (vx_kind(x, a) == VIAL_N_KWARG) {
            const vial_node_id kn = vx_child(x, a, 0), vn = vx_child(x, a, 1);
            if (vn) vx_infer(x, vn);
            const str_t kw = kn ? vx_text(x, kn) : STR_LIT("");
            uint32_t k = 0;
            while (k < np && !str_eq(kw, pnames[k])) ++k;
            if (k == np) {
                vx_diag_fix(x, "E0202", a, NULL, "`%s` takes no keyword `%.*s`", fname, (int)kw.len, kw.ptr);
                ok = false;
            } else if (args[k]) {
                vx_diag_fix(x, "E0204", a, NULL, "`%s` is given %.*s twice", fname, (int)kw.len, kw.ptr);
                ok = false;
            } else {
                args[k] = vn;
            }
            continue;
        }
        vx_infer(x, a);
        if (pos >= np) {
            vx_diag_fix(x, "E0204", a, NULL, "`%s` takes %u argument%s", fname, np, np == 1 ? "" : "s");
            ok = false;
            continue;
        }
        args[pos++] = a;
    }
    for (uint32_t k = 0; ok && k < np; ++k) {
        const vial_node_id dflt = pnodes[k] ? vial_child(t, pnodes[k], 2) : 0;
        if (!args[k] && !dflt) {
            vx_diag_fix(x, "E0204", call, NULL, "`%s` needs %.*s", fname, (int)pnames[k].len, pnames[k].ptr);
            ok = false;
        }
    }
    if (!ok) return err();

    // The frame: the def's nodes, its own symbols
    vx_frame_t f = {.parent = x->frame, .t = t, .def = def, .prelude = prelude, .depth = depth};
    f.site = prelude ? (x->frame && x->frame->prelude ? x->frame->site : call) : 0;
    f.syms.allocator = x->alloc;
    const vial_node_id lo = subtree_min(t, def, x->alloc);
    const uint32_t count = def - lo + 1;
    vial_types_t out = {.base = lo};
    md_array_resize(out.types, count, x->alloc);
    md_array_resize(out.values, count, x->alloc);
    md_array_resize(out.typed, count, x->alloc);
    for (uint32_t i = 0; i < count; ++i) out.types[i] = err();
    memset(out.values, 0, count * sizeof(*out.values));
    memset(out.typed, 0, count);
    infer_t y = *x;
    y.t = t; y.n = n; y.src = src; y.out = &out; y.frame = &f;

    // Each argument split by its parameter's kind; the mapped dimensions broadcast
    T loop = vial_type_scalar(VIAL_DT_INT);
    uint64_t key = md_hash64(&def, sizeof(def), (uint64_t)(uintptr_t)t);
    if (prelude) key = md_hash64(&f.site, sizeof(f.site), key);
    for (uint32_t k = 0; k < np; ++k) {
        const vial_node_id ann = pnodes[k] ? vial_child(t, pnodes[k], 1) : 0;
        const uint8_t kind = ann ? kind_of_annotation(text_in(t, src, ann)) : VX_KIND_WHOLE;
        T core, lp;
        vial_value_t v = {0};
        vial_node_id prov = 0;
        char pname[32];
        snprintf(pname, sizeof(pname), "%.*s", (int)(pnames[k].len < 24 ? pnames[k].len : 24), pnames[k].ptr);
        if (args[k]) {
            const T a = x->out->types[(args[k]) - x->out->base];
            if (vial_type_is_error(&a)) { ok = false; continue; }
            if (!vx_split(x, args[k], &a, kind, fname, pname, &core, &lp)) { ok = false; continue; }
            if (lp.ndim == 0) v = vx_value_of(x, args[k]);
            prov = prelude ? provenance_of(x, args[k]) : 0;
        } else {
            const vial_node_id dflt = vial_child(t, pnodes[k], 2);
            core = vx_infer(&y, dflt);
            lp = vial_type_scalar(VIAL_DT_INT);
            if (vial_type_is_error(&core)) { ok = false; continue; }
            v = vx_value_of(&y, dflt);
        }
        T m = loop;
        if (!vx_broadcast(&loop, &lp, &m)) {
            vx_diag_fix(x, "E0301", call, NULL, "the dimensions `%s` maps over do not broadcast: `%s` and `%s`", fname, vx_tstr(&loop), vx_tstr(&lp));
            return err();
        }
        loop.ndim = m.ndim;
        memcpy(loop.dims, m.dims, sizeof(m.dims));
        if (lp.stage > loop.stage) loop.stage = lp.stage;
        key = hash_type(key, &core);
        key = hash_value(key, &v);
        const vial_node_id pn = vial_child(t, pnodes[k], 0);
        const uint32_t sym = pn ? n->ref[pn] : 0;
        if (sym) vx_frame_set(&y, sym, &core, &v, prov);
    }
    if (!ok) return err();

    // Checked once per def and arguments
    key &= ~(1ull << 63);   // clear of the hashmap's reserved keys
    const uint32_t* hit = md_hashmap_get(&x->memo->map, key);
    T result;
    if (hit) {
        result = x->memo->results[*hit];
        if (loop.ndim == 0 && !prelude) vx_set_value(x, call, x->memo->values[*hit]);
    } else {
        const vial_node_id body = vial_child(t, def, 3);
        vx_check_block(&y, body);
        result = f.has_result ? f.result : vial_type_scalar(VIAL_DT_NONE);
        vial_value_t rv = {0};
        md_hashmap_add(&x->memo->map, key, (uint32_t)md_array_size(x->memo->results));
        md_array_push(x->memo->results, result, x->alloc);
        md_array_push(x->memo->values, rv, x->alloc);
    }
    if (vial_type_is_error(&result)) return result;
    T t2;
    if (!with_loop(x, call, &loop, &result, &t2)) return err();
    return t2;
}

// ### DRAWING, FOR NOW ###

bool vx_call_draw(infer_t* x, vial_node_id call, const vial_builtin_t* b, T* out) {
    const str_t name = str_from_cstr(b->name);
    static const char* prims[] = {"marker", "line", "arrow", "sphere", "arc", "disk", "label", "highlight"};
    if (str_eq(name, STR_LIT("shown"))) {
        // shown(x, d) is x, with d attached for hover. When x failed, d is not checked: one mistake, one diagnostic.
        vial_node_id first = 0;
        for (uint32_t i = 1; i < vx_num_children(x, call) && !first; ++i) {
            const vial_node_id a = vx_child(x, call, i);
            if (a && vx_kind(x, a) != VIAL_N_KWARG) first = a;
        }
        if (first && !vial_type_is_error((T[]){vx_infer(x, first)})) {
            for (uint32_t i = 1; i < vx_num_children(x, call); ++i) {
                const vial_node_id a = vx_child(x, call, i);
                if (!a || a == first) continue;
                vx_infer(x, vx_kind(x, a) == VIAL_N_KWARG ? vx_child(x, a, 1) : a);
            }
        }
        if (!first) { vx_diag_fix(x, "E0204", call, NULL, "`shown` takes a value and a drawing"); *out = err(); return true; }
        *out = x->out->types[(first) - x->out->base];
        if (!vial_type_is_error(out)) out->flags |= VIAL_TYPE_DRAWN;
        vx_set_value(x, call, vx_value_of(x, first));
        return true;
    }
    for (size_t i = 0; i < ARRAY_SIZE(prims); ++i) {
        if (!str_eq(name, str_from_cstr(prims[i]))) continue;
        // The leading points map: line(res_a.com, vlb.com) is one line per residue. The rest arrives with step 9.
        const uint32_t num_points = str_eq(name, STR_LIT("line")) ? 2 : str_eq(name, STR_LIT("highlight")) ? 0 : 1;
        T t = vial_type_scalar(VIAL_DT_DRAWING);
        T loop = vial_type_scalar(VIAL_DT_INT);
        uint32_t pos = 0;
        bool ok = true;
        for (uint32_t k = 1; k < vx_num_children(x, call); ++k) {
            const vial_node_id a = vx_child(x, call, k);
            if (!a) continue;
            const bool kw = vx_kind(x, a) == VIAL_N_KWARG;
            const vial_node_id v = kw ? vx_child(x, a, 1) : a;
            if (!v) continue;
            const T at = vx_infer(x, v);
            if (vial_type_is_error(&at)) { ok = false; continue; }
            T o = t;
            vx_join_stage(&t, &o, &at);
            if (kw || pos++ >= num_points) continue;
            T core, lp;
            if (!vx_split(x, v, &at, VX_KIND_POINT, prims[i], pos == 1 ? "a" : "b", &core, &lp)) { ok = false; continue; }
            T m = loop;
            if (!vx_broadcast(&loop, &lp, &m)) {
                vx_diag_fix(x, "E0301", call, NULL, "`%s` and `%s` do not broadcast", vx_tstr(&loop), vx_tstr(&lp));
                ok = false;
                continue;
            }
            loop.ndim = m.ndim;
            memcpy(loop.dims, m.dims, sizeof(m.dims));
        }
        if (!ok) { *out = err(); return true; }
        for (uint8_t k = 0; k < loop.ndim; ++k) vial_type_push_dim(&t, loop.dims[k]);
        *out = t;
        return true;
    }
    return false;
}

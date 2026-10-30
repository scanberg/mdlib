// The checker, step 2: names and scopes.
//
// Two walks over the tree. The first declares: it binds every top-level name in order (E0305 for a second binding,
// W0201 for one that hides a built-in), opens a scope for each def, top-level loop and comprehension, and binds
// parameters and locals. The second resolves every name a script reads. A top-level name is visible from the
// statement after the one that binds it, so the lines above a binding still see the built-in it hides; a def body
// sees every top-level name, since it runs when called. A name nothing binds is E0201, with a suggestion; a v1 name is
// E0701. Last, the calls between defs give the call graph, and a cycle in it is recursion, E0802.

#include "check.h"

#include <core/md_common.h>

#include <stdio.h>
#include <string.h>

typedef struct region_t {
    uint32_t if_beg;        // the first symbol bound inside this if
    uint32_t branch_beg;    // the first symbol bound inside the branch being declared
} region_t;

typedef struct resolver_t {
    const vial_tree_t*  t;
    str_t               src;
    vial_diags_t*       diags;
    vial_names_t*       n;
    md_allocator_i*     alloc;
    bool                is_prelude;
    const vial_names_t* prelude;        // the prelude's names, when resolving a script
    uint32_t            stmt;           // the top-level statement being resolved
    region_t*           regions;        // md_array: the top-level ifs being declared, innermost last
    md_hashset_t        shadow_warned;  // names W0201 has been given for
    uint32_t*           builtin_sym;    // md_array by built-in index: its symbol, 0 until used
    uint32_t            num_suggested;  // suggestions looked for, which is bounded
} resolver_t;

// ### HELPERS ###

static inline const vial_node_t* node(const resolver_t* r, vial_node_id id) { return &r->t->nodes[id]; }
static inline vial_node_kind_t kind_of(const resolver_t* r, vial_node_id id) { return (vial_node_kind_t)r->t->nodes[id].kind; }
static inline vial_node_id child(const resolver_t* r, vial_node_id id, uint32_t i) { return vial_child(r->t, id, i); }
static inline uint32_t num_children(const resolver_t* r, vial_node_id id) { return r->t->nodes[id].child_count; }

static inline str_t text(const resolver_t* r, vial_node_id id) {
    const md_vial_span_t s = r->t->nodes[id].span;
    return (str_t){r->src.ptr + s.beg, s.end - s.beg};
}

static inline uint32_t line_of(const resolver_t* r, vial_node_id id) {
    return vial_line_of(r->diags, r->t->nodes[id].span.beg);
}

static uint64_t name_key(str_t name) {
    uint64_t k = md_hash64_str(name, 0x5eed);
    return k >= MD_HASH_TOMBSTONE ? k - 2 : k;
}

static uint32_t scope_find(const resolver_t* r, uint32_t scope, str_t name) {
    const vial_scope_t* s = &r->n->scopes[scope];
    const uint32_t* v = md_hashmap_get(&s->names, name_key(name));
    if (v && str_eq(r->n->symbols[*v].name, name)) return *v;
    return 0;
}

static uint32_t new_symbol(resolver_t* r, uint32_t scope, str_t name, vial_sym_kind_t kind, vial_node_id at) {
    vial_symbol_t sym = {.name = name, .kind = (uint8_t)kind, .in_prelude = r->is_prelude, .scope = scope, .node = at,
                         .stmt = r->stmt};
    const uint32_t idx = (uint32_t)md_array_size(r->n->symbols);
    md_array_push(r->n->symbols, sym, r->alloc);
    return idx;
}

static uint32_t bind(resolver_t* r, uint32_t scope, str_t name, vial_sym_kind_t kind, vial_node_id at) {
    const uint32_t idx = new_symbol(r, scope, name, kind, at);
    md_hashmap_add(&r->n->scopes[scope].names, name_key(name), idx);
    if (at) r->n->ref[at] = idx;
    return idx;
}

static uint32_t new_scope(resolver_t* r, vial_scope_kind_t kind, uint32_t parent, vial_node_id at) {
    vial_scope_t s = {.kind = (uint8_t)kind, .parent = parent, .node = at};
    s.names.allocator = r->alloc;
    const uint32_t idx = (uint32_t)md_array_size(r->n->scopes);
    md_array_push(r->n->scopes, s, r->alloc);
    if (at) r->n->scope_of[at] = idx;
    return idx;
}

static void diag(resolver_t* r, const char code[6], vial_node_id at, const char* fmt, ...) {
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    vial_diag(r->diags, code, r->t->nodes[at].span, "%s", buf);
}

static void diag_fix(resolver_t* r, const char code[6], vial_node_id at, const char* fix, const char* fmt, ...) {
    char buf[400];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    vial_diag(r->diags, code, r->t->nodes[at].span, "%s", buf);
    if (fix) vial_diag_set_fix(r->diags, md_array_size(r->diags->items) - 1, str_from_cstr(fix));
}

static bool is_kind_name(str_t s) {
    static const char* kinds[] = {"set", "group", "point", "points", "float", "int", "bool"};
    for (size_t i = 0; i < ARRAY_SIZE(kinds); ++i) if (str_eq(s, str_from_cstr(kinds[i]))) return true;
    return false;
}

// A name the script cannot rebind without hiding something it may still mean
static bool is_predefined(const resolver_t* r, str_t name) {
    const vial_builtin_t* b = vial_builtin_find(name);
    if (b && (b->flags & (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_VALUE))) return true;
    if (r->prelude && scope_find(&(resolver_t){.n = (vial_names_t*)r->prelude}, 0, name)) return true;
    return false;
}

// ### DECLARING ###

static void declare_def(resolver_t* r, vial_node_id def, uint32_t parent);
static void declare_block(resolver_t* r, vial_node_id block, uint32_t scope);

// A second binding of a top-level name, unless the first is in another branch of an if being declared
static bool rebinds(const resolver_t* r, uint32_t existing) {
    for (size_t i = md_array_size(r->regions); i > 0; --i) {
        const region_t* g = &r->regions[i - 1];
        if (existing >= g->if_beg && existing < g->branch_beg) return false;
    }
    return true;
}

static void bind_top(resolver_t* r, vial_node_id name_node, vial_sym_kind_t kind) {
    const str_t name = text(r, name_node);
    if (str_eq(name, STR_LIT("__all__"))) {
        diag(r, "E0801", name_node, "`__all__` is reserved: a later version may list the published names with it");
        return;
    }
    const uint32_t existing = scope_find(r, 0, name);
    if (existing && (r->n->symbols[existing].kind == VIAL_SYM_GLOBAL || r->n->symbols[existing].kind == VIAL_SYM_DEF)) {
        if (rebinds(r, existing)) {
            diag(r, "E0305", name_node, "`%.*s` is already bound on line %u", (int)name.len, name.ptr,
                vial_line_of(r->diags, r->t->nodes[r->n->symbols[existing].node].span.beg));
            r->n->ref[name_node] = existing;
            return;
        }
        r->n->ref[name_node] = existing;     // the same name in another branch: one name
        return;
    }
    if (!r->is_prelude && is_predefined(r, name) && !md_hashset_get(&r->shadow_warned, name_key(name))) {
        md_hashset_add(&r->shadow_warned, name_key(name));
        diag_fix(r, "W0201", name_node, NULL, "`%.*s` now names a %s; the built-in is hidden from here on", (int)name.len,
            name.ptr, kind == VIAL_SYM_DEF ? "def" : "variable");
    }
    bind(r, 0, name, kind, name_node);
}

// The names an assignment target binds
static void bind_target(resolver_t* r, vial_node_id target, uint32_t scope, vial_sym_kind_t kind, bool top) {
    switch (kind_of(r, target)) {
    case VIAL_N_NAME:
        if (top) {
            bind_top(r, target, kind);
        } else {
            const str_t name = text(r, target);
            const uint32_t existing = scope_find(r, scope, name);
            if (existing) r->n->ref[target] = existing;     // a local bound again: one local
            else bind(r, scope, name, kind, target);
        }
        break;
    case VIAL_N_TUPLE:
    case VIAL_N_LIST:
        for (uint32_t i = 0; i < num_children(r, target); ++i) bind_target(r, child(r, target, i), scope, kind, top);
        break;
    default:
        break;  // a subscript or attribute target binds nothing; item assignment is refused by the parser
    }
}

static void declare_stmt(resolver_t* r, vial_node_id s, uint32_t scope) {
    const bool top = scope == 0;
    const vial_sym_kind_t local = r->n->scopes[scope].kind == VIAL_SCOPE_LOOP ? VIAL_SYM_LOCAL :
                                  top ? VIAL_SYM_GLOBAL : VIAL_SYM_LOCAL;
    switch (kind_of(r, s)) {
    case VIAL_N_ASSIGN: {
        const uint32_t n = num_children(r, s);
        for (uint32_t i = 0; i + 1 < n; ++i) bind_target(r, child(r, s, i), scope, local, top);
        break;
    }
    case VIAL_N_AUG_ASSIGN: {
        const vial_node_id target = child(r, s, 0);
        if (top && kind_of(r, target) == VIAL_N_NAME) {
            const str_t name = text(r, target);
            const uint32_t existing = scope_find(r, 0, name);
            if (existing) {
                diag(r, "E0305", s, "`%.*s` is already bound on line %u; at top level a name is bound once", (int)name.len,
                    name.ptr, vial_line_of(r->diags, r->t->nodes[r->n->symbols[existing].node].span.beg));
                r->n->ref[target] = existing;
            }
        } else if (!top) {
            bind_target(r, target, scope, local, false);
        }
        break;
    }
    case VIAL_N_DEF: {
        const vial_node_id name = child(r, s, 0);
        if (name && kind_of(r, name) == VIAL_N_NAME) {
            if (top) bind_top(r, name, VIAL_SYM_DEF);
            else bind_target(r, name, scope, VIAL_SYM_DEF, false);
        }
        declare_def(r, s, scope);
        break;
    }
    case VIAL_N_FOR: {
        if (top || r->n->scopes[scope].kind == VIAL_SCOPE_LOOP) {
            const uint32_t loop = new_scope(r, VIAL_SCOPE_LOOP, scope, s);
            bind_target(r, child(r, s, 0), loop, VIAL_SYM_LOOP, false);
            declare_block(r, child(r, s, 2), loop);
        } else {
            bind_target(r, child(r, s, 0), scope, VIAL_SYM_LOCAL, false);
            declare_block(r, child(r, s, 2), scope);
        }
        break;
    }
    case VIAL_N_WHILE:
        declare_block(r, child(r, s, 1), scope);
        break;
    case VIAL_N_IF: {
        // The branches of an if exclude each other: each may bind the same names
        region_t g = {.if_beg = (uint32_t)md_array_size(r->n->symbols)};
        g.branch_beg = g.if_beg;
        md_array_push(r->regions, g, r->alloc);
        for (vial_node_id at = s; at && kind_of(r, at) == VIAL_N_IF;) {
            md_array_last(r->regions)->branch_beg = (uint32_t)md_array_size(r->n->symbols);
            declare_block(r, child(r, at, 1), scope);
            const vial_node_id other = child(r, at, 2);
            if (!other) break;
            md_array_last(r->regions)->branch_beg = (uint32_t)md_array_size(r->n->symbols);
            if (kind_of(r, other) == VIAL_N_IF) { at = other; continue; }
            declare_block(r, other, scope);
            break;
        }
        md_array_pop(r->regions);
        break;
    }
    default:
        break;
    }
}

static void declare_block(resolver_t* r, vial_node_id block, uint32_t scope) {
    if (!block) return;
    if (kind_of(r, block) != VIAL_N_BLOCK && kind_of(r, block) != VIAL_N_MODULE) { declare_stmt(r, block, scope); return; }
    for (uint32_t i = 0; i < num_children(r, block); ++i) declare_stmt(r, child(r, block, i), scope);
}

static void declare_def(resolver_t* r, vial_node_id def, uint32_t parent) {
    const uint32_t scope = new_scope(r, VIAL_SCOPE_DEF, parent, def);
    const vial_node_id params = child(r, def, 1);
    for (uint32_t i = 0; params && i < num_children(r, params); ++i) {
        const vial_node_id p = child(r, params, i);
        if (kind_of(r, p) != VIAL_N_PARAM) continue;
        const vial_node_id name = child(r, p, 0);
        if (name && kind_of(r, name) == VIAL_N_NAME) bind(r, scope, text(r, name), VIAL_SYM_PARAM, name);
        const vial_node_id ann = child(r, p, 1);
        if (ann) {
            const str_t a = text(r, ann);
            if (kind_of(r, ann) != VIAL_N_NAME || !is_kind_name(a)) {
                char ex[VIAL_EXCERPT_MAX + 4];
                const str_t e = vial_excerpt(ex, sizeof(ex), a);
                diag_fix(r, "E0312", ann, "", "`%.*s` is not a parameter kind; the kinds are set, group, point, points, float, int and bool",
                    (int)e.len, e.ptr);
            }
        }
    }
    const vial_node_id ret = child(r, def, 2);
    if (ret) diag_fix(r, "E0312", ret, "", "a def's result is not annotated; it is inferred at each call");
    declare_block(r, child(r, def, 3), scope);
}

// ### RESOLVING ###

static uint32_t builtin_symbol(resolver_t* r, const vial_builtin_t* b) {
    const size_t i = (size_t)(b - vial_builtin_at(0));
    if (!r->builtin_sym[i]) {
        r->builtin_sym[i] = new_symbol(r, 0, str_from_cstr(b->name), VIAL_SYM_BUILTIN, 0);
        r->n->symbols[r->builtin_sym[i]].builtin = b;
    }
    return r->builtin_sym[i];
}

// A prelude def by name, as a symbol of this tree's table that points at it
static uint32_t prelude_symbol(resolver_t* r, str_t name) {
    if (!r->prelude) return 0;
    const vial_scope_t* s = &r->prelude->scopes[0];
    const uint32_t* v = md_hashmap_get(&s->names, name_key(name));
    if (!v || !str_eq(r->prelude->symbols[*v].name, name) || r->prelude->symbols[*v].kind != VIAL_SYM_DEF) return 0;
    const uint32_t idx = new_symbol(r, 0, name, VIAL_SYM_PRELUDE, r->prelude->symbols[*v].node);
    r->n->symbols[idx].in_prelude = true;
    return idx;
}

static bool in_def(const resolver_t* r, uint32_t scope) {
    for (uint32_t s = scope; s; s = r->n->scopes[s].parent) if (r->n->scopes[s].kind == VIAL_SCOPE_DEF) return true;
    return false;
}

// Looks a name up from a scope; *late is set when a top-level binding exists but comes later
static uint32_t lookup(resolver_t* r, uint32_t scope, str_t name, uint32_t* late) {
    *late = 0;
    const bool deferred = in_def(r, scope);
    for (uint32_t s = scope;; s = r->n->scopes[s].parent) {
        const uint32_t found = scope_find(r, s, name);
        if (found) {
            const vial_symbol_t* sym = &r->n->symbols[found];
            if (s != 0 || deferred || r->is_prelude || sym->stmt < r->stmt) return found;
            *late = found;
            break;
        }
        if (s == 0) break;
    }
    const uint32_t pre = prelude_symbol(r, name);
    if (pre) return pre;
    const vial_builtin_t* b = vial_builtin_find(name);
    if (b && (b->flags & (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_VALUE))) return builtin_symbol(r, b);
    return 0;
}

// Edit distance with transpositions, for suggestions, computed only near the diagonal: more than limit is limit + 1
static size_t edit_distance(str_t a, str_t b, size_t limit) {
    const size_t big = limit + 1;
    if (a.len > 32 || b.len > 32) return big;
    if (a.len > b.len + limit || b.len > a.len + limit) return big;
    if (a.len == 0 || b.len == 0) return a.len + b.len <= limit ? a.len + b.len : big;
    size_t d[34][34];
    // Only the band is read: each row's cells from one left of it to one right of it
    for (size_t j = 0; j <= b.len && j <= limit + 1; ++j) d[0][j] = j <= limit ? j : big;
    for (size_t i = 1; i <= a.len; ++i) {
        const size_t jlo = i > limit ? i - limit : 1;
        const size_t jhi = i + limit < b.len ? i + limit : b.len;
        d[i][jlo - 1] = (jlo == 1 && i <= limit) ? i : big;      // the cell left of the band
        if (i > 1 && jhi > i - 1 + limit) d[i - 1][jhi] = big;    // the cell above it, right of the row before's band
        size_t row_min = big;
        for (size_t j = jlo; j <= jhi; ++j) {
            const size_t cost = a.ptr[i - 1] == b.ptr[j - 1] ? 0 : 1;
            size_t v = d[i - 1][j - 1] + cost;
            if (d[i - 1][j] + 1 < v) v = d[i - 1][j] + 1;
            if (d[i][j - 1] + 1 < v) v = d[i][j - 1] + 1;
            if (i > 1 && j > 1 && a.ptr[i - 1] == b.ptr[j - 2] && a.ptr[i - 2] == b.ptr[j - 1] && d[i - 2][j - 2] + 1 < v) v = d[i - 2][j - 2] + 1;
            d[i][j] = v > big ? big : v;
            if (d[i][j] < row_min) row_min = d[i][j];
        }
        if (row_min > limit) return big;
    }
    return d[a.len][b.len];
}

typedef struct best_t { str_t name; size_t dist; } best_t;

static void consider(best_t* best, str_t word, str_t candidate) {
    if (str_empty(candidate) || str_eq(word, candidate)) return;
    if (word.len < 3) return;
    if (candidate.len + 2 < word.len || word.len + 2 < candidate.len) return;   // every short name is one edit from many others: a guess there is noise
    const size_t limit = word.len <= 5 ? 1 : 2;
    const size_t d = edit_distance(word, candidate, limit);
    if (d <= limit && d < best->dist) { best->dist = d; best->name = candidate; }
}

static str_t suggest(resolver_t* r, uint32_t scope, str_t word, bool functions_only) {
    const char* syn = vial_builtin_synonym(word);
    if (syn) return str_from_cstr(syn);
    // Each suggestion looks at every visible name; past the first twenty unknown names a script is mostly
    // mistakes, and more guesses would only make compiling slow
    if (++r->num_suggested > 20) return (str_t){0};
    best_t best = {.dist = SIZE_MAX};
    for (uint32_t s = scope;; s = r->n->scopes[s].parent) {
        const vial_scope_t* sc = &r->n->scopes[s];
        for (uint32_t i = 0; i < sc->names.num_buckets; ++i) {
            if (md_hashmap_skip_index(&sc->names, i)) continue;
            const vial_symbol_t* sym = &r->n->symbols[sc->names.values[i]];
            if (functions_only && sym->kind != VIAL_SYM_DEF) continue;
            consider(&best, word, sym->name);
        }
        if (s == 0) break;
    }
    if (r->prelude) {
        const vial_scope_t* sc = &r->prelude->scopes[0];
        for (uint32_t i = 0; i < sc->names.num_buckets; ++i) {
            if (md_hashmap_skip_index(&sc->names, i)) continue;
            const vial_symbol_t* sym = &r->prelude->symbols[sc->names.values[i]];
            if (sym->kind == VIAL_SYM_DEF && sym->name.len && sym->name.ptr[0] != '_') consider(&best, word, sym->name);
        }
    }
    for (size_t i = 0; i < vial_builtin_count(); ++i) {
        const vial_builtin_t* b = vial_builtin_at(i);
        if (b->flags == VIAL_BUILTIN_KIND) continue;
        if (functions_only && !(b->flags & (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_PROPERTY))) continue;
        consider(&best, word, str_from_cstr(b->name));
    }
    return best.name;
}

static void report_unknown(resolver_t* r, vial_node_id at, uint32_t scope, uint32_t late) {
    const str_t name = text(r, at);
    if (vial_diag_at(r->diags, r->t->nodes[at].span.beg)) return;   // the parser said more, as E0701 for residue(1:4)
    const char* v1_msg = NULL;
    const char* v1_fix = NULL;
    if (vial_builtin_v1(name, &v1_msg, &v1_fix)) {
        diag_fix(r, "E0701", at, v1_fix, "%s", v1_msg);
        return;
    }
    if (str_eq(name, STR_LIT("print"))) {
        diag(r, "E0803", at, "`print` is not part of VIAL; hover shows any value");
        return;
    }
    if (late) {
        diag(r, "E0201", at, "`%.*s` is used before it is bound, on line %u", (int)name.len, name.ptr,
            vial_line_of(r->diags, r->t->nodes[r->n->symbols[late].node].span.beg));
        return;
    }
    if (vial_unit_name(name, NULL)) {
        // 1.2 * nm: the fix writes the quantity literal
        const vial_node_id p = r->n->parent[at];
        char fix[96] = "";
        if (p && kind_of(r, p) == VIAL_N_BINARY && node(r, p)->sub == VIAL_OP_STAR && child(r, p, 1) == at) {
            const vial_node_id num = child(r, p, 0);
            if (kind_of(r, num) == VIAL_N_INT || kind_of(r, num) == VIAL_N_FLOAT) {
                const str_t ns = text(r, num);
                snprintf(fix, sizeof(fix), "%.*s %.*s", (int)ns.len, ns.ptr, (int)name.len, name.ptr);
            }
        }
        diag_fix(r, "E0201", at, fix[0] ? fix : NULL, "`%.*s` is a unit, not a name: write it after a number, as 10 %.*s",
            (int)name.len, name.ptr, (int)name.len, name.ptr);
        return;
    }
    const str_t s = suggest(r, scope, name, false);
    char fix[64] = "";
    if (!str_empty(s)) snprintf(fix, sizeof(fix), "%.*s", (int)s.len, s.ptr);
    diag_fix(r, "E0201", at, fix[0] ? fix : NULL, "`%.*s` is not defined", (int)name.len, name.ptr);
}

static void resolve(resolver_t* r, vial_node_id id, uint32_t scope);

static void resolve_children(resolver_t* r, vial_node_id id, uint32_t scope) {
    for (uint32_t i = 0; i < num_children(r, id); ++i) {
        const vial_node_id c = child(r, id, i);
        if (c) resolve(r, c, scope);
    }
}

// x.name: a built-in function or property, a prelude def or a def, called with x first
static void resolve_attr(resolver_t* r, vial_node_id id, uint32_t scope) {
    const vial_node_id obj = child(r, id, 0);
    const vial_node_id name_node = child(r, id, 1);
    if (obj) resolve(r, obj, scope);
    if (!name_node || kind_of(r, name_node) != VIAL_N_NAME) return;
    const str_t name = text(r, name_node);
    const vial_builtin_t* b = vial_builtin_find(name);
    if (b && (b->flags & (VIAL_BUILTIN_FUNC | VIAL_BUILTIN_PROPERTY))) { r->n->ref[name_node] = builtin_symbol(r, b); return; }
    uint32_t late;
    const uint32_t found = lookup(r, scope, name, &late);
    if (found && (r->n->symbols[found].kind == VIAL_SYM_DEF || r->n->symbols[found].kind == VIAL_SYM_PRELUDE)) {
        r->n->ref[name_node] = found;
        return;
    }
    if (vial_diag_at(r->diags, r->t->nodes[name_node].span.beg)) return;
    // .xy and the like: VIAL has no swizzles
    bool swizzle = name.len >= 2 && name.len <= 4;
    for (size_t i = 0; i < name.len && swizzle; ++i) swizzle = name.ptr[i] == 'x' || name.ptr[i] == 'y' || name.ptr[i] == 'z';
    if (swizzle) {
        const str_t o = text(r, obj);
        char fix[160];
        const int a = name.ptr[0] - 'x', b2 = name.ptr[1] - 'x';
        if (name.len == 2 && b2 == a + 1) snprintf(fix, sizeof(fix), "%.*s[..., %d:%d]", (int)o.len, o.ptr, a, a + 2);
        else {
            int n = snprintf(fix, sizeof(fix), "%.*s[..., [", (int)o.len, o.ptr);
            for (size_t i = 0; i < name.len && n > 0 && (size_t)n < sizeof(fix); ++i) n += snprintf(fix + n, sizeof(fix) - n, "%s%d", i ? ", " : "", name.ptr[i] - 'x');
            if (n > 0 && (size_t)n < sizeof(fix)) snprintf(fix + n, sizeof(fix) - n, "]]");
        }
        diag_fix(r, "E0201", name_node, fix, "`.%.*s` is not defined: VIAL has no swizzles; write %s", (int)name.len, name.ptr, fix);
        return;
    }
    const str_t s = suggest(r, scope, name, true);
    char fix[64] = "";
    if (!str_empty(s)) snprintf(fix, sizeof(fix), "%.*s", (int)s.len, s.ptr);
    diag_fix(r, "E0201", name_node, fix[0] ? fix : NULL, "`.%.*s` is not a function or property", (int)name.len, name.ptr);
}

static void resolve_comp(resolver_t* r, vial_node_id id, uint32_t scope) {
    const uint32_t comp = new_scope(r, VIAL_SCOPE_COMP, scope, id);
    const uint32_t n = num_children(r, id);
    for (uint32_t i = 1; i < n; ++i) {
        const vial_node_id f = child(r, id, i);
        if (!f || kind_of(r, f) != VIAL_N_COMP_FOR) continue;
        // The first iterable is read where the comprehension is; the rest, and the conditions, inside it
        const vial_node_id it = child(r, f, 1);
        if (it) resolve(r, it, i == 1 ? scope : comp);
        bind_target(r, child(r, f, 0), comp, VIAL_SYM_COMP, false);
        for (uint32_t k = 2; k < num_children(r, f); ++k) if (child(r, f, k)) resolve(r, child(r, f, k), comp);
    }
    if (child(r, id, 0)) resolve(r, child(r, id, 0), comp);
}

static void resolve(resolver_t* r, vial_node_id id, uint32_t scope) {
    switch (kind_of(r, id)) {
    case VIAL_N_NAME: {
        if (r->n->ref[id]) return;      // bound here
        uint32_t late;
        const uint32_t found = lookup(r, scope, text(r, id), &late);
        if (found) r->n->ref[id] = found;
        else report_unknown(r, id, scope, late);
        return;
    }
    case VIAL_N_ATTR:
        resolve_attr(r, id, scope);
        return;
    case VIAL_N_KWARG:
        if (child(r, id, 1)) resolve(r, child(r, id, 1), scope);
        return;
    case VIAL_N_LISTCOMP:
        resolve_comp(r, id, scope);
        return;
    case VIAL_N_UNSUPPORTED:
    case VIAL_N_UNIT:
    case VIAL_N_ERROR:
        return;     // reported already, or nothing to resolve
    case VIAL_N_DEF: {
        const uint32_t def_scope = r->n->scope_of[id];
        const vial_node_id params = child(r, id, 1);
        for (uint32_t i = 0; params && i < num_children(r, params); ++i) {
            const vial_node_id p = child(r, params, i);
            if (kind_of(r, p) == VIAL_N_PARAM && child(r, p, 2)) resolve(r, child(r, p, 2), scope);  // defaults
        }
        if (child(r, id, 3)) resolve(r, child(r, id, 3), def_scope);
        return;
    }
    case VIAL_N_FOR: {
        if (child(r, id, 1)) resolve(r, child(r, id, 1), scope);
        const uint32_t loop = r->n->scope_of[id] ? r->n->scope_of[id] : scope;
        if (child(r, id, 0)) resolve(r, child(r, id, 0), loop);
        if (child(r, id, 2)) resolve(r, child(r, id, 2), loop);
        return;
    }
    case VIAL_N_ASSIGN: {
        const uint32_t n = num_children(r, id);
        if (n) resolve(r, child(r, id, n - 1), scope);
        for (uint32_t i = 0; i + 1 < n; ++i) resolve(r, child(r, id, i), scope);
        return;
    }
    default:
        resolve_children(r, id, scope);
        return;
    }
}

// ### RECURSION ###

typedef struct callgraph_t {
    vial_node_id*   defs;       // md_array: the DEF nodes
    md_hashmap32_t  def_index;  // DEF node -> its index in defs
    uint8_t*        color;      // 0 unseen, 1 on the stack, 2 done
} callgraph_t;

// The calls in a def's own body, not in the defs nested in it
static void calls_in(resolver_t* r, vial_node_id id, vial_node_id** out) {
    for (uint32_t i = 0; i < num_children(r, id); ++i) {
        const vial_node_id c = child(r, id, i);
        if (!c || kind_of(r, c) == VIAL_N_DEF || kind_of(r, c) == VIAL_N_UNSUPPORTED) continue;
        if (kind_of(r, c) == VIAL_N_CALL) md_array_push(*out, c, r->alloc);
        calls_in(r, c, out);
    }
}

static uint32_t callee_symbol(resolver_t* r, vial_node_id call) {
    vial_node_id callee = child(r, call, 0);
    if (callee && kind_of(r, callee) == VIAL_N_ATTR) callee = child(r, callee, 1);
    if (!callee || kind_of(r, callee) != VIAL_N_NAME) return 0;
    const uint32_t s = r->n->ref[callee];
    return (s && r->n->symbols[s].kind == VIAL_SYM_DEF) ? s : 0;
}

// The def a call calls, as an index into g->defs, or UINT32_MAX
static uint32_t callee_def(resolver_t* r, const callgraph_t* g, vial_node_id call) {
    const uint32_t s = callee_symbol(r, call);
    if (!s) return UINT32_MAX;
    const vial_node_id target = r->n->parent[r->n->symbols[s].node];
    if (!target || r->t->nodes[target].kind != VIAL_N_DEF) return UINT32_MAX;
    const uint32_t* idx = md_hashmap_get(&g->def_index, (uint64_t)target);
    return idx ? *idx : UINT32_MAX;
}

static void report_cycle(resolver_t* r, const callgraph_t* g, uint32_t d, uint32_t e, vial_node_id call) {
    const str_t a = text(r, child(r, g->defs[d], 0));
    const str_t b = text(r, child(r, g->defs[e], 0));
    if (e == d) diag(r, "E0802", call, "`%.*s` calls itself; VIAL has no recursion", (int)a.len, a.ptr);
    else diag(r, "E0802", call, "`%.*s` calls `%.*s`, which calls `%.*s` again; VIAL has no recursion",
        (int)a.len, a.ptr, (int)b.len, b.ptr, (int)a.len, a.ptr);
}

// Depth first over the call graph, with an explicit stack: a script may chain thousands of defs
typedef struct frame_t { uint32_t def; uint32_t next; vial_node_id* calls; } frame_t;

static void check_recursion(resolver_t* r) {
    callgraph_t g = {.def_index = {.allocator = r->alloc}};
    for (size_t id = 1; id < md_array_size(r->t->nodes); ++id) {
        if (r->t->nodes[id].kind == VIAL_N_DEF && child(r, (vial_node_id)id, 0)) {
            md_hashmap_add(&g.def_index, (uint64_t)id, (uint32_t)md_array_size(g.defs));
            md_array_push(g.defs, (vial_node_id)id, r->alloc);
        }
    }
    const size_t n = md_array_size(g.defs);
    if (!n) return;
    g.color = md_alloc(r->alloc, n);
    memset(g.color, 0, n);
    frame_t* stack = NULL;
    for (size_t root = 0; root < n; ++root) {
        if (g.color[root]) continue;
        frame_t f = {.def = (uint32_t)root};
        calls_in(r, child(r, g.defs[root], 3), &f.calls);
        g.color[root] = 1;
        md_array_push(stack, f, r->alloc);
        while (md_array_size(stack)) {
            frame_t* top = md_array_last(stack);
            if (top->next >= md_array_size(top->calls)) {
                g.color[top->def] = 2;
                md_array_pop(stack);
                continue;
            }
            const vial_node_id call = top->calls[top->next++];
            const uint32_t e = callee_def(r, &g, call);
            if (e == UINT32_MAX) continue;
            if (g.color[e] == 1) {
                report_cycle(r, &g, top->def, e, call);
            } else if (g.color[e] == 0) {
                frame_t nf = {.def = e};
                calls_in(r, child(r, g.defs[e], 3), &nf.calls);
                g.color[e] = 1;
                md_array_push(stack, nf, r->alloc);
            }
        }
    }
}

// ### ENTRY ###

static void names_init(vial_names_t* n, const vial_tree_t* t, md_allocator_i* alloc) {
    memset(n, 0, sizeof(*n));
    const size_t num = md_array_size(t->nodes);
    md_array_resize(n->ref, num, alloc);
    md_array_resize(n->scope_of, num, alloc);
    md_array_resize(n->parent, num, alloc);
    memset(n->ref, 0, num * sizeof(*n->ref));
    memset(n->scope_of, 0, num * sizeof(*n->scope_of));
    memset(n->parent, 0, num * sizeof(*n->parent));
    for (size_t id = 1; id < num; ++id) {
        const vial_node_t* nd = &t->nodes[id];
        for (uint32_t i = 0; i < nd->child_count; ++i) {
            const vial_node_id c = t->children[nd->child_beg + i];
            if (c) n->parent[c] = (vial_node_id)id;
        }
    }
    vial_symbol_t none = {0};
    md_array_push(n->symbols, none, alloc);
}

static void run(resolver_t* r) {
    md_array_resize(r->builtin_sym, vial_builtin_count(), r->alloc);
    memset(r->builtin_sym, 0, vial_builtin_count() * sizeof(*r->builtin_sym));
    r->shadow_warned.allocator = r->alloc;

    const vial_node_id root = r->t->root;
    new_scope(r, VIAL_SCOPE_MODULE, 0, 0);
    r->n->scope_of[root] = 0;
    for (uint32_t i = 0; i < num_children(r, root); ++i) {
        r->stmt = i;
        declare_stmt(r, child(r, root, i), 0);
    }
    for (uint32_t i = 0; i < num_children(r, root); ++i) {
        r->stmt = i;
        resolve(r, child(r, root, i), 0);
    }
    check_recursion(r);
}

void vial_check_names(vial_check_t* c) {
    if (c->prelude) {
        names_init(&c->prelude_names, c->prelude, c->alloc);
        resolver_t pr = {.t = c->prelude, .src = c->prelude_src, .diags = c->prelude_diags, .n = &c->prelude_names,
                         .alloc = c->alloc, .is_prelude = true};
        run(&pr);
    }
    names_init(&c->names, c->tree, c->alloc);
    resolver_t r = {.t = c->tree, .src = c->src, .diags = c->diags, .n = &c->names, .alloc = c->alloc,
                    .prelude = c->prelude ? &c->prelude_names : NULL};
    run(&r);
    vial_diags_sort(c->diags);
}

const vial_symbol_t* vial_names_symbol(const vial_names_t* n, vial_node_id id) {
    if (id >= md_array_size(n->ref) || !n->ref[id]) return NULL;
    return &n->symbols[n->ref[id]];
}

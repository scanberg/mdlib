// A SMILES of atoms of a system, canonical: the same graph is written the same way whatever the order of the atoms.
//
// The graph is the atoms given and the bonds among them, with their elements, charges, bond orders and aromaticity
// as the system has them. The order is found as Weininger's CANON does: atoms are ranked by what they are, the ranks
// are refined by the ranks of their neighbours until they no longer split, and ties that remain (symmetric atoms) are
// broken one at a time, refining again after each. The string is then written depth first from the lowest rank, each
// atom's neighbours in rank order, its last one on the main chain and the others as branches.
//
// What is written is a pattern as much as a molecule: hydrogens are left implicit unless given as atoms, a bracket
// atom carries no hydrogen count, and a bond whose order is not known is '~', any bond, which md_smiles_parse reads
// as a query does.

#include <md_smiles.h>
#include <md_system.h>
#include <md_element.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_sort.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct sw_atom_t {
    int32_t  idx;           // in the system
    uint8_t  z;
    bool     aromatic;
    int8_t   charge;
    uint32_t cls;           // atom class, 0 for none
    uint32_t h;             // hydrogens not written as atoms: what a bracket atom states
    uint32_t num_nbr;
    uint32_t nbr_beg;       // into nbr and nbr_bond
    uint32_t rank;
} sw_atom_t;

typedef struct sw_t {
    sw_atom_t* atoms;
    uint32_t   n;
    uint32_t*  nbr;         // local atom indices
    uint8_t*   nbr_bond;    // md_smiles_bond_order_t
    // writing
    char*      buf;
    size_t     cap;
    size_t     len;
    bool*      visited;
    uint32_t*  parent;
    // ring closures: for each atom, the partners it opens or closes a ring bond with, and the digit
    uint32_t*  ring_beg;    // [n + 1] into ring_partner
    uint32_t*  ring_partner;
    uint8_t*   ring_bond;
    int*       ring_digit;  // per ring-closure entry, -1 until assigned
    bool       digits[100];
    uint32_t   depth;
    bool       too_deep;
    md_allocator_i* alloc;
} sw_t;

static void put(sw_t* w, const char* s) {
    for (; *s; ++s) {
        if (w->len + 1 < w->cap) w->buf[w->len] = *s;
        w->len += 1;
    }
    if (w->cap) w->buf[MD_MIN(w->len, w->cap - 1)] = '\0';
}

// ### RANKS ###

// What an atom is, before its neighbours are looked at, as a key in that order: its class, element, aromatic or not,
// charge and the number of its neighbours (at most 64)
static uint64_t invariant_key(const sw_atom_t* a) {
    return ((uint64_t)a->cls << 24) | ((uint64_t)a->z << 16) | ((uint64_t)a->aromatic << 15) | ((uint64_t)(uint8_t)(a->charge + 128) << 7) | a->num_nbr;
}

typedef struct sw_key_t {
    uint32_t  atom;
    uint32_t  rank;
    uint64_t* nbr;          // sorted (rank << 8 | bond) of the neighbours
    uint32_t  num;
} sw_key_t;

// A few at most: a list of neighbours
static void sort_small_u64(uint64_t* v, uint32_t n) {
    for (uint32_t i = 1; i < n; ++i) {
        const uint64_t x = v[i];
        uint32_t j = i;
        while (j > 0 && v[j - 1] > x) { v[j] = v[j - 1]; --j; }
        v[j] = x;
    }
}

static int cmp_key(const void* pa, const void* pb) {
    const sw_key_t* a = pa;
    const sw_key_t* b = pb;
    if (a->rank != b->rank) return a->rank < b->rank ? -1 : 1;
    if (a->num != b->num) return a->num < b->num ? -1 : 1;
    for (uint32_t i = 0; i < a->num; ++i) {
        if (a->nbr[i] != b->nbr[i]) return a->nbr[i] < b->nbr[i] ? -1 : 1;
    }
    return 0;
}

static int cmp_key_atom(const void* pa, const void* pb) {
    const sw_key_t* a = pa;
    const sw_key_t* b = pb;
    const int c = cmp_key(a, b);
    if (c) return c;
    return a->atom < b->atom ? -1 : a->atom > b->atom;
}

// The end of the class (atoms of one rank) which starts at place s of order
static inline uint32_t class_end(const sw_t* w, const uint32_t* order, uint32_t s) {
    uint32_t e = s + 1;
    while (e < w->n && w->atoms[order[e]].rank == w->atoms[order[s]].rank) ++e;
    return e;
}

// Refines the ranks by the neighbours' until the classes (atoms of one rank) no longer split; returns their number. A
// rank is the place of the first atom of its class in the order of (rank, the neighbours' ranks, atom): the atoms are
// kept in that order, and a class of one atom stays where it is, so a round sorts the classes of several only. order
// and fresh are of n atoms, count of n + 2.
static uint32_t refine(sw_t* w, sw_key_t* keys, uint32_t* order, uint32_t* fresh, uint32_t* count) {
    const uint32_t n = w->n;
    // The atoms by rank (a broken tie may have moved one up to n), and each rank the place of its class
    md_memset(count, 0, (n + 2) * sizeof(uint32_t));
    for (uint32_t i = 0; i < n; ++i) count[w->atoms[i].rank + 1] += 1;
    for (uint32_t r = 0; r <= n; ++r) count[r + 1] += count[r];
    for (uint32_t i = 0; i < n; ++i) order[count[w->atoms[i].rank]++] = i;
    for (uint32_t i = 0; i < n; ++i) fresh[i] = (i && w->atoms[order[i]].rank == w->atoms[order[i - 1]].rank) ? fresh[i - 1] : i;
    for (uint32_t i = 0; i < n; ++i) w->atoms[order[i]].rank = fresh[i];

    for (;;) {
        // Each class of several: the keys of its atoms from the ranks of this round, sorted, ranked anew by place
        bool split = false;
        for (uint32_t s = 0, e; s < n; s = e) {
            e = class_end(w, order, s);
            if (e - s < 2) continue;
            for (uint32_t i = s; i < e; ++i) {
                sw_key_t* k = &keys[i];
                const sw_atom_t* a = &w->atoms[order[i]];
                k->atom = order[i];
                k->rank = a->rank;
                k->num  = a->num_nbr;
                for (uint32_t j = 0; j < a->num_nbr; ++j) {
                    const uint32_t o = w->nbr[a->nbr_beg + j];
                    k->nbr[j] = ((uint64_t)w->atoms[o].rank << 8) | w->nbr_bond[a->nbr_beg + j];
                }
                sort_small_u64(k->nbr, k->num);
            }
            qsort(keys + s, e - s, sizeof(sw_key_t), cmp_key_atom);
            for (uint32_t i = s; i < e; ++i) {
                const bool same = i > s && cmp_key(&keys[i - 1], &keys[i]) == 0;
                fresh[i] = same ? fresh[i - 1] : i;
                split |= i > s && !same;
            }
        }
        if (!split) break;
        // The new ranks only now, as every key of the round is of the ranks before it
        for (uint32_t s = 0, e; s < n; s = e) {
            e = class_end(w, order, s);
            if (e - s < 2) continue;
            for (uint32_t i = s; i < e; ++i) order[i] = keys[i].atom;
            for (uint32_t i = s; i < e; ++i) w->atoms[order[i]].rank = fresh[i];
        }
    }
    uint32_t classes = 0;
    for (uint32_t i = 0; i < n; ++i) classes += !i || w->atoms[order[i]].rank != w->atoms[order[i - 1]].rank;
    return classes;
}

static void rank_atoms(sw_t* w) {
    // Initial ranks by what each atom is: the place of the first of its kind
    uint64_t* inv   = md_alloc(w->alloc, w->n * sizeof(uint64_t));
    uint32_t* order = md_alloc(w->alloc, w->n * sizeof(uint32_t));
    for (uint32_t i = 0; i < w->n; ++i) inv[i] = invariant_key(&w->atoms[i]);
    md_argsort_radix_u64(order, inv, w->n);
    for (uint32_t i = 0; i < w->n; ++i) {
        const bool same = i && inv[order[i - 1]] == inv[order[i]];
        w->atoms[order[i]].rank = same ? w->atoms[order[i - 1]].rank : i;
    }
    // Each key's list holds as many as any atom has: the keys are sorted, and so change atoms, every round
    uint32_t max_nbr = 1;
    for (uint32_t i = 0; i < w->n; ++i) max_nbr = MD_MAX(max_nbr, w->atoms[i].num_nbr);
    sw_key_t* keys = md_alloc(w->alloc, w->n * sizeof(sw_key_t));
    for (uint32_t i = 0; i < w->n; ++i) keys[i].nbr = md_alloc(w->alloc, max_nbr * sizeof(uint64_t));
    uint32_t* held  = md_alloc(w->alloc, w->n * sizeof(uint32_t));      // Atoms of each rank (ranks are below n)
    uint32_t* fresh = md_alloc(w->alloc, w->n * sizeof(uint32_t));
    uint32_t* count = md_alloc(w->alloc, (w->n + 2) * sizeof(uint32_t));
    uint32_t classes = refine(w, keys, order, fresh, count);
    while (classes < w->n) {
        // Break the lowest tie: the first atom of the lowest rank held by several moves ahead of the others
        md_memset(held, 0, w->n * sizeof(uint32_t));
        for (uint32_t i = 0; i < w->n; ++i) held[w->atoms[i].rank] += 1;
        uint32_t best_rank = 0;
        while (best_rank < w->n && held[best_rank] < 2) ++best_rank;
        if (best_rank == w->n) break;
        uint32_t pick = 0;
        while (w->atoms[pick].rank != best_rank) ++pick;
        // The tied atoms move up by one place, the one picked stays: it is now ranked before them
        for (uint32_t i = 0; i < w->n; ++i) if (i != pick && w->atoms[i].rank >= best_rank) w->atoms[i].rank += 1;
        classes = refine(w, keys, order, fresh, count);
    }
}

// ### WRITING ###

static bool organic(uint8_t z, bool aromatic) {
    switch (z) {
    case 5: case 6: case 7: case 8: case 15: case 16: return true;        // B C N O P S, also aromatic
    case 9: case 17: case 35: case 53: return !aromatic;                   // F Cl Br I
    default: return false;
    }
}

static void put_atom(sw_t* w, const sw_atom_t* a) {
    char sym[8];
    if (a->z == 0) snprintf(sym, sizeof(sym), "*");
    else {
        const md_str_t s = md_atomic_number_symbol((md_atomic_number_t)a->z);
        snprintf(sym, sizeof(sym), "%.*s", (int)MD_MIN(s.len, 3), s.ptr);
        if (a->aromatic && sym[0]) {
            sym[0] = (char)(sym[0] - 'A' + 'a');
        }
    }
    const bool bracket = a->cls || a->charge || a->z == 1 || (a->z && !organic(a->z, a->aromatic));
    if (!bracket) { put(w, sym); return; }
    // A bracket states the hydrogens and the charge, as a query reads it, unless it is there only for the class of an
    // atom of the organic subset, which leaves them open
    const bool states = a->charge || a->z == 1 || (a->z && !organic(a->z, a->aromatic));
    char b[48];
    int n = snprintf(b, sizeof(b), "[%s", sym);
    if (states && a->h) n += a->h == 1 ? snprintf(b + n, sizeof(b) - n, "H") : snprintf(b + n, sizeof(b) - n, "H%u", a->h);
    if (a->charge) {
        const int c = a->charge < 0 ? -a->charge : a->charge;
        n += c == 1 ? snprintf(b + n, sizeof(b) - n, "%c", a->charge < 0 ? '-' : '+')
                    : snprintf(b + n, sizeof(b) - n, "%c%d", a->charge < 0 ? '-' : '+', c);
    }
    if (a->cls) n += snprintf(b + n, sizeof(b) - n, ":%u", a->cls);
    snprintf(b + n, sizeof(b) - n, "]");
    put(w, b);
}

static void put_bond(sw_t* w, uint8_t bond, const sw_atom_t* a, const sw_atom_t* b) {
    switch (bond) {
    case MD_SMILES_BOND_SINGLE:    if (a->aromatic && b->aromatic) put(w, "-"); break;     // between rings, not in one
    case MD_SMILES_BOND_DOUBLE:    put(w, "="); break;
    case MD_SMILES_BOND_TRIPLE:    put(w, "#"); break;
    case MD_SMILES_BOND_QUADRUPLE: put(w, "$"); break;
    case MD_SMILES_BOND_AROMATIC:  if (!(a->aromatic && b->aromatic)) put(w, ":"); break;
    case MD_SMILES_BOND_ANY:       put(w, "~"); break;
    default: break;
    }
}

static uint8_t bond_between(const sw_t* w, uint32_t a, uint32_t b) {
    const sw_atom_t* x = &w->atoms[a];
    for (uint32_t j = 0; j < x->num_nbr; ++j) if (w->nbr[x->nbr_beg + j] == b) return w->nbr_bond[x->nbr_beg + j];
    return MD_SMILES_BOND_ANY;
}

// The neighbours of an atom in rank order
static void sorted_nbrs(const sw_t* w, uint32_t a, uint32_t* out) {
    const sw_atom_t* x = &w->atoms[a];
    for (uint32_t j = 0; j < x->num_nbr; ++j) out[j] = w->nbr[x->nbr_beg + j];
    for (uint32_t i = 1; i < x->num_nbr; ++i) {
        const uint32_t v = out[i];
        uint32_t j = i;
        while (j > 0 && w->atoms[out[j - 1]].rank > w->atoms[v].rank) { out[j] = out[j - 1]; --j; }
        out[j] = v;
    }
}

#define SW_MAX_DEPTH 2048

// First pass: the depth-first tree, and the bonds that close rings, in the order the second pass meets them
typedef struct sw_ring_t { uint32_t open, close; } sw_ring_t;

static void find_rings(sw_t* w, uint32_t a, uint32_t from, sw_ring_t** rings) {
    if (++w->depth > SW_MAX_DEPTH) { w->too_deep = true; --w->depth; return; }
    w->visited[a] = true;
    w->parent[a] = from;
    uint32_t nb[64];
    const uint32_t nn = MD_MIN(w->atoms[a].num_nbr, 64u);
    sorted_nbrs(w, a, nb);
    for (uint32_t j = 0; j < nn; ++j) {
        const uint32_t o = nb[j];
        if (o == from) continue;
        if (w->visited[o]) {
            // An earlier atom: a ring bond, recorded once, from the atom written first
            bool seen = false;
            for (size_t r = 0; r < md_array_size(*rings) && !seen; ++r) seen = ((*rings)[r].open == o && (*rings)[r].close == a) || ((*rings)[r].open == a && (*rings)[r].close == o);
            if (!seen) md_array_push(*rings, ((sw_ring_t){o, a}), w->alloc);
            continue;
        }
        find_rings(w, o, a, rings);
    }
    --w->depth;
}

static void write_from(sw_t* w, uint32_t a, uint32_t from, const sw_ring_t* rings, int* ring_digit, bool* done) {
    if (++w->depth > SW_MAX_DEPTH) { w->too_deep = true; --w->depth; return; }
    done[a] = true;
    put_atom(w, &w->atoms[a]);
    // Ring bonds at this atom: closing the ones opened before, then opening new ones, each in the order of the partner's
    // rank; a bond symbol goes with the opening digit
    const size_t nr = md_array_size(rings);
    for (int pass = 0; pass < 2; ++pass) {
        for (;;) {
            size_t best = SIZE_MAX;
            for (size_t r = 0; r < nr; ++r) {
                const bool closes = pass == 0 && rings[r].close == a && ring_digit[r] >= 0;
                const bool opens = pass == 1 && rings[r].open == a && ring_digit[r] == -1;
                if (!closes && !opens) continue;
                const uint32_t partner = closes ? rings[r].open : rings[r].close;
                const uint32_t pb = best == SIZE_MAX ? 0 : (pass == 0 ? rings[best].open : rings[best].close);
                if (best == SIZE_MAX || w->atoms[partner].rank < w->atoms[pb].rank) best = r;
            }
            if (best == SIZE_MAX) break;
            char d[8];
            if (pass == 0) {
                const int dg = ring_digit[best];
                snprintf(d, sizeof(d), dg < 10 ? "%d" : "%%%02d", dg);
                put(w, d);
                w->digits[dg] = false;
                ring_digit[best] = -2;      // closed
            } else {
                int dg = 1;
                while (dg < 100 && w->digits[dg]) ++dg;
                if (dg >= 100) { w->too_deep = true; break; }
                w->digits[dg] = true;
                ring_digit[best] = dg;
                put_bond(w, bond_between(w, rings[best].open, rings[best].close), &w->atoms[rings[best].open], &w->atoms[rings[best].close]);
                snprintf(d, sizeof(d), dg < 10 ? "%d" : "%%%02d", dg);
                put(w, d);
            }
        }
    }
    // The children of the tree: all but the last as branches
    uint32_t nb[64], kids[64];
    uint32_t nk = 0;
    const uint32_t nn = MD_MIN(w->atoms[a].num_nbr, 64u);
    sorted_nbrs(w, a, nb);
    for (uint32_t j = 0; j < nn; ++j) {
        const uint32_t o = nb[j];
        if (o != from && w->parent[o] == a && !done[o]) kids[nk++] = o;
    }
    for (uint32_t k = 0; k < nk; ++k) {
        const bool branch = k + 1 < nk;
        if (branch) put(w, "(");
        put_bond(w, bond_between(w, a, kids[k]), &w->atoms[a], &w->atoms[kids[k]]);
        write_from(w, kids[k], a, rings, ring_digit, done);
        if (branch) put(w, ")");
    }
    --w->depth;
}

// The place of atom o among the atoms given (lookup: atom << 32 | place, sorted), UINT32_MAX when it is not one
static uint32_t local_of(const uint64_t* lookup, uint32_t n, int32_t o) {
    uint32_t lo = 0, hi = n;
    while (lo < hi) {
        const uint32_t mid = (lo + hi) / 2;
        if ((uint32_t)(lookup[mid] >> 32) < (uint32_t)o) lo = mid + 1;
        else hi = mid;
    }
    return (lo < n && (uint32_t)(lookup[lo] >> 32) == (uint32_t)o) ? (uint32_t)lookup[lo] : UINT32_MAX;
}

size_t md_smiles_write(char* buf, size_t cap, const md_system_t* sys, const int32_t* atoms, size_t count, const uint32_t* classes) {
    if (buf && cap) buf[0] = '\0';
    if (!sys || !atoms || count == 0 || count > MD_SMILES_WRITE_MAX) return 0;
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    sw_t w = {.buf = buf, .cap = buf ? cap : 0, .alloc = alloc};
    w.n = (uint32_t)count;
    w.atoms = md_alloc(alloc, count * sizeof(sw_atom_t));
    memset(w.atoms, 0, count * sizeof(sw_atom_t));
    size_t result = 0;

    // The atoms given, in range and each once, sorted with their places to look them up
    uint64_t* lookup = md_alloc(alloc, count * sizeof(uint64_t));
    for (uint32_t i = 0; i < w.n; ++i) {
        const int32_t ai = atoms[i];
        if (ai < 0 || (size_t)ai >= sys->atom.count) goto done;
        lookup[i] = ((uint64_t)(uint32_t)ai << 32) | i;
    }
    md_sort_radix_u64(lookup, w.n);
    for (uint32_t i = 1; i < w.n; ++i) if ((lookup[i] >> 32) == (lookup[i - 1] >> 32)) goto done;

    // The atoms, and their bonds among themselves
    size_t num_conn = 0;
    for (uint32_t i = 0; i < w.n; ++i) {
        const int32_t ai = atoms[i];
        const md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)ai);
        num_conn += it.end_idx - it.i;
        sw_atom_t* a = &w.atoms[i];
        a->idx = ai;
        a->z = (uint8_t)md_atom_atomic_number(&sys->atom, (size_t)ai);
        a->charge = (int8_t)md_atom_formal_charge(&sys->atom, (size_t)ai);
        a->cls = classes ? classes[i] : 0;
    }
    uint32_t* nbr      = md_alloc(alloc, MD_MAX(1, num_conn) * sizeof(uint32_t));
    uint8_t*  nbr_bond = md_alloc(alloc, MD_MAX(1, num_conn));
    uint32_t  num_nbr  = 0;
    for (uint32_t i = 0; i < w.n; ++i) {
        sw_atom_t* a = &w.atoms[i];
        a->nbr_beg = num_nbr;
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)a->idx);
        while (md_bond_iter_has_next(&it)) {
            const int32_t o = md_bond_iter_atom_index(&it);
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            md_bond_iter_next(&it);
            const md_bond_kind_t kind = md_bond_kind(f);
            if (kind != MD_BOND_KIND_COVALENT && kind != MD_BOND_KIND_COORDINATE) continue;
            // Aromatic as the system has it, through any of its bonds: a ring atom given without its ring is still one
            if (f & MD_BOND_FLAG_AROMATIC) a->aromatic = true;
            const uint32_t local = local_of(lookup, w.n, o);
            if (local == UINT32_MAX) continue;
            uint8_t b;
            if (f & MD_BOND_FLAG_AROMATIC) b = MD_SMILES_BOND_AROMATIC;
            else switch (md_bond_order(f)) {
                case 1: b = MD_SMILES_BOND_SINGLE; break;
                case 2: b = MD_SMILES_BOND_DOUBLE; break;
                case 3: b = MD_SMILES_BOND_TRIPLE; break;
                case 4: b = MD_SMILES_BOND_QUADRUPLE; break;
                default: b = MD_SMILES_BOND_ANY; break;
            }
            nbr[num_nbr]      = local;
            nbr_bond[num_nbr] = b;
            num_nbr++;
        }
        a->num_nbr = num_nbr - a->nbr_beg;
        if (a->num_nbr > 64) goto done;
        // Its hydrogens, explicit and implicit as the system has them, less those written as atoms
        uint32_t bonded_h = 0, written_h = 0;
        md_bond_iter_t hit = md_bond_iter(&sys->bond, (size_t)a->idx);
        while (md_bond_iter_has_next(&hit)) {
            const int32_t o = md_bond_iter_atom_index(&hit);
            md_bond_iter_next(&hit);
            if (md_atom_atomic_number(&sys->atom, (size_t)o) != 1) continue;
            bonded_h += 1;
            written_h += local_of(lookup, w.n, o) != UINT32_MAX;
        }
        const uint32_t total_h = sys->atom.hydrogen_count ? MD_MAX((uint32_t)sys->atom.hydrogen_count[a->idx], bonded_h) : bonded_h;
        a->h = total_h > written_h ? total_h - written_h : 0;
    }
    w.nbr = nbr;
    w.nbr_bond = nbr_bond;

    rank_atoms(&w);

    // Depth first from the lowest rank, a component at a time
    w.visited = md_alloc(alloc, w.n * sizeof(bool));
    w.parent = md_alloc(alloc, w.n * sizeof(uint32_t));
    bool* done_ = md_alloc(alloc, w.n * sizeof(bool));
    memset(w.visited, 0, w.n * sizeof(bool));
    memset(done_, 0, w.n * sizeof(bool));
    uint32_t* by_rank = md_alloc(alloc, w.n * sizeof(uint32_t));
    for (uint32_t i = 0; i < w.n; ++i) by_rank[w.atoms[i].rank] = i;
    bool first = true;
    for (uint32_t r = 0; r < w.n; ++r) {
        const uint32_t start = by_rank[r];
        if (w.visited[start]) continue;
        md_array(sw_ring_t) rings = 0;
        find_rings(&w, start, UINT32_MAX, &rings);
        if (w.too_deep) goto done;
        int* ring_digit = md_alloc(alloc, (md_array_size(rings) + 1) * sizeof(int));
        for (size_t k = 0; k < md_array_size(rings); ++k) ring_digit[k] = -1;
        if (!first) put(&w, ".");
        first = false;
        write_from(&w, start, UINT32_MAX, rings, ring_digit, done_);
        if (w.too_deep) goto done;
    }
    result = w.len;

done:
    if (!result && buf && cap) buf[0] = '\0';
    md_temp_end(temp);
    return result;
}

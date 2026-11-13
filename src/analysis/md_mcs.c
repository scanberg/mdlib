// Maximum common substructure (md_mcs.h)
//
// Every scratch array is from one temp scope, opened in md_mcs and passed down; nothing here is freed on its own.

#include <md_mcs.h>
#include <md_system.h>
#include <md_match.h>

#include <core/md_allocator.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_geom.h>
#include <core/md_label.h>
#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_sort.h>
#include <core/md_vec_math.h>

#include <limits.h>
#include <math.h>
#include <stdlib.h>

#define MCS_DEFAULT_MAX_STEPS 10000000

// ### ATOMS ###

// What an atom is matched by: its element, and with MD_MCS_FLAG_BY_NAME its name as well (the 8 bytes of its label)
typedef struct mcs_label_t {
    uint64_t name;
    uint32_t z;
} mcs_label_t;

static inline uint64_t name_key(md_str_t name) {
    return md_label_key(md_label_make(name));
}

static inline mcs_label_t atom_label(const md_system_t* sys, int32_t a, bool by_name) {
    return (mcs_label_t){by_name ? name_key(md_atom_name(&sys->atom, a)) : 0, md_atom_atomic_number(&sys->atom, a)};
}

static inline bool heavy(const md_system_t* sys, int32_t a) {
    return md_atom_atomic_number(&sys->atom, a) != 1;
}

static inline bool in_range(const md_system_t* sys, int32_t a) {
    return a >= 0 && (size_t)a < sys->atom.count;
}

// The bonds of a, an upper bound of its neighbours within a set
static inline uint32_t bond_count(const md_bond_data_t* bond, uint32_t a) {
    const md_bond_iter_t it = md_bond_iter(bond, a);
    return it.end_idx - it.i;
}

// Sorted arr[n] each once; how many are left
static size_t unique_u32(uint32_t* arr, size_t n) {
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) if (!k || arr[k - 1] != arr[i]) arr[k++] = arr[i];
    return k;
}

// The index of x in sorted arr[n], -1 when it is not there
static int find_u32(const uint32_t* arr, int n, uint32_t x) {
    int lo = 0, hi = n;
    while (lo < hi) { const int mid = (lo + hi) / 2; if (arr[mid] < x) lo = mid + 1; else hi = mid; }
    return (lo < n && arr[lo] == x) ? lo : -1;
}

// ### GRAPHS ###

// A set of vertices: a bit for each, in words of 64
#define MCS_MAX_WORDS ((MD_MCS_GRAPH_MAX + 63) / 64)

static inline int words_of(int n) { return (n + 63) / 64; }

typedef struct mcs_graph_t {
    int             n;
    int             words;      // of a set of its vertices
    int32_t*        atom;       // [n] the atom of each vertex, by degree, highest first
    mcs_label_t*    label;      // [n]
    uint64_t*       adj;        // [n * words] the neighbours of each vertex, as a set
} mcs_graph_t;

// The heavy atoms of atoms[count], each once, and the bonds between them. The vertices are by their degree within the
// set, highest first and ties by atom: McSplit branches on the lowest vertex. False for more than MD_MCS_GRAPH_MAX.
static bool graph_init(mcs_graph_t* g, const md_system_t* sys, const int32_t* atoms, size_t count, bool by_name, md_temp_scope_t temp) {
    *g = (mcs_graph_t){0};
    uint32_t* h = md_temp_alloc_array(temp, uint32_t, count);
    size_t num = 0;
    for (size_t i = 0; i < count; ++i) if (in_range(sys, atoms[i]) && heavy(sys, atoms[i])) h[num++] = (uint32_t)atoms[i];
    md_sort_radix_u32(h, num);
    num = unique_u32(h, num);
    if (num > MD_MCS_GRAPH_MAX) return false;
    const int n = (int)num;
    if (!n) return true;

    // The neighbours of each vertex within the set, as rows nbr[off[i], off[i + 1])
    size_t cap = 0;
    for (int i = 0; i < n; ++i) cap += bond_count(&sys->bond, h[i]);
    uint32_t* off = md_temp_alloc_array(temp, uint32_t, n + 1);
    int*      nbr = md_temp_alloc_array(temp, int, cap);
    uint32_t k = 0, max_deg = 0;
    for (int i = 0; i < n; ++i) {
        off[i] = k;
        md_bond_iter_t it = md_bond_iter(&sys->bond, h[i]);
        while (md_bond_iter_has_next(&it)) {
            const int j = find_u32(h, n, md_bond_iter_atom_index(&it));
            md_bond_iter_next(&it);
            if (j < 0 || j == i) continue;
            bool listed = false;
            for (uint32_t m = off[i]; m < k; ++m) listed |= nbr[m] == j;
            if (!listed) nbr[k++] = j;
        }
        max_deg = MD_MAX(max_deg, k - off[i]);
    }
    off[n] = k;

    // Each vertex's place: a counting sort by degree, highest first, stable
    uint32_t* first = md_temp_alloc_zero_array(temp, uint32_t, max_deg + 1);
    for (int i = 0; i < n; ++i) first[off[i + 1] - off[i]]++;
    for (uint32_t d = max_deg + 1, sum = 0; d-- > 0;) { const uint32_t c = first[d]; first[d] = sum; sum += c; }
    int* pos = md_temp_alloc_array(temp, int, n);
    for (int i = 0; i < n; ++i) pos[i] = (int)first[off[i + 1] - off[i]]++;

    const int w = words_of(n);
    g->n     = n;
    g->words = w;
    g->atom  = md_temp_alloc_array(temp, int32_t, n);
    g->label = md_temp_alloc_array(temp, mcs_label_t, n);
    g->adj   = md_temp_alloc_zero_array(temp, uint64_t, (size_t)n * w);
    for (int i = 0; i < n; ++i) {
        const int p = pos[i];
        g->atom[p]  = (int32_t)h[i];
        g->label[p] = atom_label(sys, (int32_t)h[i], by_name);
        for (uint32_t m = off[i]; m < off[i + 1]; ++m) {
            const int q = pos[nbr[m]];
            g->adj[(size_t)p * w + q / 64] |= 1ull << (q % 64);
            g->adj[(size_t)q * w + p / 64] |= 1ull << (p % 64);
        }
    }
    return true;
}

// ### McSPLIT ###
// Classes (bidomains): the vertices of each graph that may still be matched to each other, alike in label and in
// their adjacency to every pair matched so far. A bound is the matched pairs and the least of each class's two sides.
//
// The two sides of a class are sets of vertices, a bit each: splitting a class by its adjacency to a pair is two ANDs
// a word, the lowest vertex of a side is its first bit, and what the search tries is in ascending order of vertex,
// as McSplit's. A class is a record of words, its first word then its two sides of W words each; W is a constant of each
// instance of the search (1, 2, 4, 8 or 16: up to 64, 128, 256, 512 and 1024 vertices).

typedef struct vpair_t { int v, w; } vpair_t;

typedef struct search_t {
    const mcs_graph_t*  g0;
    const mcs_graph_t*  g1;
    const uint64_t* adj0;       // [g0->n * W] the neighbours of each vertex of g0, as sets of W words
    const uint64_t* adj1;       // [g1->n * W]
    uint64_t*   cls;            // The classes of each level of the search, a stack of records of 1 + 2W words
    size_t      top, cap;       // In classes
    vpair_t*    cur;
    int         cur_n;
    vpair_t*    best;
    int         best_n;
    uint64_t    steps, max_steps;
    bool        gave_up;
    md_allocator_i* arena;      // What the stack grows in
} search_t;

// The first word of a class: how many on the left, on the right, and whether it is adjacent to a pair matched (so
// that what is matched stays connected through it). The sizes are kept for sides of more than one word; a side of one
// is counted when its size is wanted, a POPCNT, which is less than keeping it up to date.
#define KEEP_SIZES(W) ((W) > 1)

static inline uint64_t cls_head(uint32_t left, uint32_t right, bool adjacent) {
    return (uint64_t)left | ((uint64_t)right << 16) | ((uint64_t)adjacent << 32);
}
static inline bool cls_adj(uint64_t h) { return (h >> 32) != 0; }

static MD_FORCE_INLINE uint64_t* cls_at(const search_t* s, size_t c, const int W) {
    return s->cls + c * (size_t)(1 + 2 * W);
}

static MD_FORCE_INLINE int set_count(const uint64_t* set, const int W) {
    int n = 0;
    for (int j = 0; j < W; ++j) n += (int)md_popcnt64(set[j]);
    return n;
}

// The sizes of the two sides of class c
static MD_FORCE_INLINE int cls_left(const uint64_t* c, const int W) {
    return KEEP_SIZES(W) ? (int)(c[0] & 0xFFFF) : set_count(c + 1, W);
}
static MD_FORCE_INLINE int cls_right(const uint64_t* c, const int W) {
    return KEEP_SIZES(W) ? (int)((c[0] >> 16) & 0xFFFF) : set_count(c + 1 + W, W);
}

static MD_FORCE_INLINE int first_bit(const uint64_t* set, const int W) {
    for (int j = 0; j < W; ++j) if (set[j]) return j * 64 + (int)md_ctz64(set[j]);
    return -1;
}

// Room for need classes; the stack is addressed by index, so it may move
static void reserve(search_t* s, size_t need, int W) {
    if (need <= s->cap) return;
    const size_t cap = MD_MAX(need, 2 * s->cap);
    const size_t rec = (size_t)(1 + 2 * W) * sizeof(uint64_t);
    s->cls = md_realloc(s->arena, s->cls, s->cap * rec, cap * rec);
    s->cap = cap;
}

// The classes d[nd] split by adjacency to v and to w, onto the top of the stack; how many. Each class gives its part
// adjacent to neither, then its part adjacent to both, where each has vertices on either side.
static MD_FORCE_INLINE int filter(search_t* s, size_t d, int nd, int v, int w, const int W) {
    const uint64_t* av = s->adj0 + (size_t)v * W;
    const uint64_t* aw = s->adj1 + (size_t)w * W;
    const size_t out = s->top;
    int k = 0;
    for (int i = 0; i < nd; ++i) {
        const uint64_t* c = cls_at(s, d + i, W);
        const uint64_t* L = c + 1;
        const uint64_t* R = c + 1 + W;
        uint64_t* o = cls_at(s, out + k, W);
        // ln, rn: the sizes where they are kept, otherwise not zero where the side is not empty
        uint64_t ln = 0, rn = 0;
        for (int j = 0; j < W; ++j) { o[1 + j]     = L[j] & ~av[j]; ln = KEEP_SIZES(W) ? ln + md_popcnt64(o[1 + j])     : (ln | o[1 + j]); }
        for (int j = 0; j < W; ++j) { o[1 + W + j] = R[j] & ~aw[j]; rn = KEEP_SIZES(W) ? rn + md_popcnt64(o[1 + W + j]) : (rn | o[1 + W + j]); }
        if (ln && rn) {
            o[0] = KEEP_SIZES(W) ? cls_head((uint32_t)ln, (uint32_t)rn, cls_adj(c[0])) : c[0];
            o = cls_at(s, out + ++k, W);
        }
        // The rest of each side is adjacent; with the sizes kept, whether it is empty is known before it is made
        if (KEEP_SIZES(W)) {
            const int la = cls_left(c, W) - (int)ln, ra = cls_right(c, W) - (int)rn;
            if (!la || !ra) continue;
            for (int j = 0; j < W; ++j) o[1 + j]     = L[j] & av[j];
            for (int j = 0; j < W; ++j) o[1 + W + j] = R[j] & aw[j];
            o[0] = cls_head((uint32_t)la, (uint32_t)ra, true);
            k++;
        } else {
            uint64_t la = 0, ra = 0;
            for (int j = 0; j < W; ++j) { o[1 + j]     = L[j] & av[j]; la |= o[1 + j]; }
            for (int j = 0; j < W; ++j) { o[1 + W + j] = R[j] & aw[j]; ra |= o[1 + W + j]; }
            if (la && ra) { o[0] = cls_head(0, 0, true); k++; }
        }
    }
    return k;
}

// The class with the fewest on its larger side, ties to the one with the lowest vertex; once something is matched,
// only classes adjacent to it, so what is matched stays connected
static MD_FORCE_INLINE int select_class(const search_t* s, size_t d, int nd, const int W) {
    int min_size = INT_MAX, min_tie = INT_MAX, best = -1;
    for (int i = 0; i < nd; ++i) {
        const uint64_t* c = cls_at(s, d + i, W);
        if (s->cur_n > 0 && !cls_adj(c[0])) continue;
        const int len = MD_MAX(cls_left(c, W), cls_right(c, W));
        if (len > min_size) continue;
        const int tie = first_bit(c + 1, W);
        if (len < min_size || tie < min_tie) { min_size = len; min_tie = tie; best = i; }
    }
    return best;
}

typedef void (*solve_fn)(search_t* s, size_t d, int nd);

// The classes d[nd] of the stack. The lowest vertex v of a class is matched to each of the other side's in turn,
// lowest first; then the class is searched without v, which is the loop (in place of a call in the tail). recurse is
// the instance of this of the same W.
static MD_FORCE_INLINE void solve_w(search_t* s, size_t d, int nd, const int W, solve_fn recurse) {
    for (;;) {
        if (++s->steps > s->max_steps) { s->gave_up = true; return; }
        if (s->cur_n > s->best_n) { md_memcpy(s->best, s->cur, s->cur_n * sizeof(vpair_t)); s->best_n = s->cur_n; }
        int bound = 0;
        for (int i = 0; i < nd; ++i) {
            const uint64_t* c = cls_at(s, d + i, W);
            bound += MD_MIN(cls_left(c, W), cls_right(c, W));
        }
        if (s->cur_n + bound <= s->best_n) return;
        const int bi = select_class(s, d, nd, W);
        if (bi < 0) return;
        const size_t b = d + (size_t)bi;

        // v out of the class while it is matched, and each w in turn: the classes below split what is left of it
        uint64_t* c = cls_at(s, b, W);
        const int v = first_bit(c + 1, W);
        c[1 + v / 64] &= ~(1ull << (v % 64));
        if (KEEP_SIZES(W)) c[0] -= 1 + (1 << 16);
        uint64_t todo[MCS_MAX_WORDS];
        for (int j = 0; j < W; ++j) todo[j] = c[1 + W + j];
        for (int j = 0; j < W; ++j) {
            for (uint64_t t = todo[j]; t; t &= t - 1) {
                const uint64_t bit = t & -t;
                const int w = j * 64 + (int)md_ctz64(t);
                cls_at(s, b, W)[1 + W + j] &= ~bit;
                reserve(s, s->top + 2 * (size_t)nd, W);
                const int k = filter(s, d, nd, v, w, W);
                const size_t below = s->top;
                s->top += k;
                s->cur[s->cur_n++] = (vpair_t){v, w};
                recurse(s, below, k);
                s->cur_n--;
                s->top = below;
                cls_at(s, b, W)[1 + W + j] |= bit;
                if (s->gave_up) return;
            }
        }
        c = cls_at(s, b, W);
        if (KEEP_SIZES(W)) c[0] += 1 << 16;
        if (cls_left(c, W) == 0) {
            // The last class in its place
            const size_t last = d + (size_t)--nd;
            if (b != last) md_memcpy(c, cls_at(s, last, W), (size_t)(1 + 2 * W) * sizeof(uint64_t));
        }
    }
}

static void solve_1 (search_t* s, size_t d, int nd) { solve_w(s, d, nd,  1, solve_1); }
static void solve_2 (search_t* s, size_t d, int nd) { solve_w(s, d, nd,  2, solve_2); }
static void solve_4 (search_t* s, size_t d, int nd) { solve_w(s, d, nd,  4, solve_4); }
static void solve_8 (search_t* s, size_t d, int nd) { solve_w(s, d, nd,  8, solve_8); }
static void solve_16(search_t* s, size_t d, int nd) { solve_w(s, d, nd, 16, solve_16); }

typedef struct label_entry_t {
    uint64_t    name;
    uint32_t    z;
    uint32_t    id;         // the vertex, the top bit set for one of g1
} label_entry_t;

static int cmp_label_entry(const void* x, const void* y) {
    const label_entry_t* a = x;
    const label_entry_t* b = y;
    if (a->z != b->z) return a->z < b->z ? -1 : 1;
    if (a->name != b->name) return a->name < b->name ? -1 : 1;
    return (a->id > b->id) - (a->id < b->id);
}

// The neighbours of each vertex of g as sets of W words, W at least g's own
static const uint64_t* adj_words(const mcs_graph_t* g, int W, md_temp_scope_t temp) {
    if (W == g->words) return g->adj;
    uint64_t* adj = md_temp_alloc_zero_array(temp, uint64_t, (size_t)g->n * W);
    for (int v = 0; v < g->n; ++v) md_memcpy(adj + (size_t)v * W, g->adj + (size_t)v * g->words, g->words * sizeof(uint64_t));
    return adj;
}

// The largest connected common subgraph of g0 and g1: its pairs into out[min(g0->n, g1->n)] (vertices), how many;
// gave_up when the search ran out of steps
static int graph_mcs(vpair_t* out, const mcs_graph_t* g0, const mcs_graph_t* g1, uint64_t max_steps, bool* gave_up, md_temp_scope_t temp) {
    *gave_up = false;
    const int n0 = g0->n, n1 = g1->n, nmin = MD_MIN(n0, n1);
    if (!nmin) return 0;
    // The words of a set, a power of two
    int W = 1;
    while (W < MD_MAX(g0->words, g1->words)) W *= 2;
    MD_ASSERT(W <= MCS_MAX_WORDS);
    search_t s = {.g0 = g0, .g1 = g1, .max_steps = max_steps, .best = out, .arena = md_temp_allocator(temp)};
    s.adj0 = adj_words(g0, W, temp);
    s.adj1 = adj_words(g1, W, temp);
    s.cur  = md_temp_alloc_array(temp, vpair_t, nmin);

    // The first classes: the vertices of both by label
    const int ne = n0 + n1;
    label_entry_t* e = md_temp_alloc_array(temp, label_entry_t, ne);
    for (int v = 0; v < n0; ++v) e[v]      = (label_entry_t){g0->label[v].name, g0->label[v].z, (uint32_t)v};
    for (int w = 0; w < n1; ++w) e[n0 + w] = (label_entry_t){g1->label[w].name, g1->label[w].z, (uint32_t)w | 0x80000000u};
    qsort(e, ne, sizeof(label_entry_t), cmp_label_entry);

    // The stack last: it grows in place as the search goes deeper, the arena's last allocation
    reserve(&s, 4 * (size_t)nmin + 16, W);
    int nd = 0;
    for (int i = 0, j; i < ne; i = j) {
        uint64_t* c = cls_at(&s, nd, W);
        md_memset(c, 0, (size_t)(1 + 2 * W) * sizeof(uint64_t));
        uint32_t nl = 0, nr = 0;
        for (j = i; j < ne && e[j].z == e[i].z && e[j].name == e[i].name; ++j) {
            const uint32_t x = e[j].id & 0x7FFFFFFFu;
            if (e[j].id & 0x80000000u) { c[1 + W + x / 64] |= 1ull << (x % 64); nr++; }
            else                       { c[1 + x / 64]     |= 1ull << (x % 64); nl++; }
        }
        if (nl && nr) { c[0] = cls_head(nl, nr, false); nd++; }
    }
    s.top = nd;
    switch (W) {
    case 1:  solve_1 (&s, 0, nd); break;
    case 2:  solve_2 (&s, 0, nd); break;
    case 4:  solve_4 (&s, 0, nd); break;
    case 8:  solve_8 (&s, 0, nd); break;
    default: solve_16(&s, 0, nd); break;
    }
    *gave_up = s.gave_up;
    return s.best_n;
}

// ### SUPERPOSITION ###

// The positions of atoms idx[n] of xyz less their mean
static void centred(md_float4_t* out, const md_float3_t* xyz, const int32_t* idx, size_t n) {
    double c[3] = {0};
    for (size_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) c[k] += xyz[idx[i]].elem[k];
    for (int k = 0; k < 3; ++k) c[k] /= n;
    for (size_t i = 0; i < n; ++i) {
        const md_float3_t p = xyz[idx[i]];
        out[i] = md_float4_set((float)(p.x - c[0]), (float)(p.y - c[1]), (float)(p.z - c[2]), 1.0f);
    }
}

// RMSD of q superposed onto p, both centred, equal weights
static double rmsd_centred(const md_float4_t* p, const md_float4_t* q, size_t n) {
    const md_float3x3_t R = md_geom_optimal_rotation_rel(p, q, n);
    double ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const md_float3_t r = md_float3x3_mul_float3(R, md_float3_from_float4(q[i]));
        for (int k = 0; k < 3; ++k) { const double d = (double)r.elem[k] - p[i].elem[k]; ss += d * d; }
    }
    return sqrt(ss / n);
}

static double pairs_rmsd(const md_float3_t* xa, const int32_t* a, const md_float3_t* xb, const int32_t* b, size_t n, md_temp_scope_t temp) {
    if (n < 3) return NAN;
    md_float4_t* pa = md_temp_alloc_array(temp, md_float4_t, n);
    md_float4_t* pb = md_temp_alloc_array(temp, md_float4_t, n);
    centred(pa, xa, a, n);
    centred(pb, xb, b, n);
    return rmsd_centred(pa, pb, n);
}

// Among the mappings of a's common atoms onto b's heavy atoms that the bonds (and the labels) allow, the one of the
// least RMSD
typedef struct sym_t {
    const md_float3_t*  xb;
    const md_float4_t*  pa;     // [n] a's atoms, centred: the same for every mapping
    md_float4_t*        pb;     // [n] the mapping's atoms of b, centred
    int32_t*            best;   // [n]
    size_t              n;
    double              best_rmsd;
} sym_t;

// Called from within md_match_for_each, given no callback_alloc: it allocates nothing
static bool sym_cb(const int32_t* atom_idx, size_t width, uint32_t unit, void* user) {
    (void)unit;
    sym_t* s = user;
    if (width != s->n) return true;
    for (size_t i = 0; i < width; ++i) if (atom_idx[i] < 0) return true;
    centred(s->pb, s->xb, atom_idx, width);
    const double r = rmsd_centred(s->pa, s->pb, width);
    if (r < s->best_rmsd) { s->best_rmsd = r; md_memcpy(s->best, atom_idx, width * sizeof(int32_t)); }
    return true;
}

static void least_rmsd(const md_mcs_desc_t* desc, const int32_t* a, int32_t* b, size_t n, const mcs_graph_t* g1, bool by_name, md_temp_scope_t temp) {
    md_allocator_i* arena = md_temp_allocator(temp);
    // The query and the unit are of the scope: they go with it
    md_match_query_t q;
    if (!md_match_query_init_atoms(&q, a, n, by_name ? MD_MATCH_LABEL_NAME : MD_MATCH_LABEL_ELEMENT, desc->a.sys, arena)) return;
    const md_bitset_t unit = md_bitset_from_indices((const uint32_t*)g1->atom, (size_t)g1->n, arena);

    md_float4_t* pa = md_temp_alloc_array(temp, md_float4_t, n);
    sym_t s = {.xb = desc->b.xyz, .pa = pa, .n = n};
    s.pb   = md_temp_alloc_array(temp, md_float4_t, n);
    s.best = md_temp_alloc_array(temp, int32_t, n);
    centred(pa, desc->a.xyz, a, n);
    centred(s.pb, desc->b.xyz, b, n);
    s.best_rmsd = rmsd_centred(pa, s.pb, n);
    md_memcpy(s.best, b, n * sizeof(int32_t));

    const md_match_desc_t md = {.query = &q, .mode = MD_MATCH_MODE_ALL, .hydrogens = MD_MATCH_RESOLVE_NEVER, .bond_orders = MD_MATCH_RESOLVE_NEVER,
                                .charges = MD_MATCH_RESOLVE_NEVER, .units = &unit, .num_units = 1, .limit = 4096, .max_steps = 1000000};
    md_match_for_each(NULL, NULL, &md, desc->b.sys, sym_cb, &s, NULL);
    md_memcpy(b, s.best, n * sizeof(int32_t));
}

// ### BY RESIDUE ###

typedef struct res_list_t {
    uint32_t    n;          // Residues
    uint32_t*   res;        // [n] the residues of the set, in order
    uint64_t*   res_name;   // [n] their names, as label keys
    uint32_t*   beg;        // [n + 1] the atoms of each, atoms[beg[i], beg[i + 1])
    int32_t*    atoms;      // [beg[n]] the atoms of the set, ascending, each once
    uint64_t*   name;       // [beg[n]] the name of each, as a label key
    bool*       heavy;      // [beg[n]]
} res_list_t;

// The residues of a side, when it is two or more of them, all amino acids or nucleotides
static bool residues_of(res_list_t* out, const md_mcs_side_t* s, md_temp_scope_t temp) {
    *out = (res_list_t){0};
    const md_system_t* sys = s->sys;
    const uint32_t* offset = sys->residue.atom_offset;
    if (!sys->residue.count || !offset) return false;
    uint32_t* at = md_temp_alloc_array(temp, uint32_t, s->count);
    size_t n = 0;
    for (size_t i = 0; i < s->count; ++i) if (in_range(sys, s->atoms[i])) at[n++] = (uint32_t)s->atoms[i];
    md_sort_radix_u32(at, n);
    n = unique_u32(at, n);
    if (n < 2) return false;
    const int32_t r0 = md_system_residue_find_by_atom_idx(sys, at[0]);
    if (r0 < 0) return false;

    uint32_t* res = md_temp_alloc_array(temp, uint32_t, n);
    uint32_t* beg = md_temp_alloc_array(temp, uint32_t, n + 1);
    uint32_t nr = 0;
    for (uint32_t i = 0, r = (uint32_t)r0; i < n; ++i) {
        while (r + 1 < sys->residue.count && offset[r + 1] <= at[i]) ++r;
        if (at[i] < offset[r] || at[i] >= offset[r + 1]) return false;
        if (nr && res[nr - 1] == r) continue;
        const md_residue_kind_t k = md_system_residue_kind(sys, r);
        if (k != MD_RESIDUE_KIND_AMINO_ACID && k != MD_RESIDUE_KIND_NUCLEOTIDE) return false;
        res[nr] = r;
        beg[nr] = i;
        nr++;
    }
    beg[nr] = (uint32_t)n;
    if (nr < 2) return false;

    out->n = nr;
    out->res = res;
    out->beg = beg;
    out->atoms = (int32_t*)at;
    out->res_name = md_temp_alloc_array(temp, uint64_t, nr);
    out->name     = md_temp_alloc_array(temp, uint64_t, n);
    out->heavy    = md_temp_alloc_array(temp, bool, n);
    for (uint32_t i = 0; i < nr; ++i) out->res_name[i] = name_key(md_system_residue_name(sys, res[i]));
    for (size_t i = 0; i < n; ++i) {
        out->name[i]  = name_key(md_atom_name(&sys->atom, at[i]));
        out->heavy[i] = heavy(sys, (int32_t)at[i]);
    }
    return true;
}

// The first atom of atoms [beg, end) of a residue list of the name, -1 when none is
static int32_t atom_of_name(const res_list_t* rl, uint32_t beg, uint32_t end, uint64_t name) {
    for (uint32_t p = beg; p < end; ++p) if (rl->name[p] == name) return rl->atoms[p];
    return -1;
}

// Small ids of the residue names of both lists into ida[ra->n] and idb[rb->n], one where the names are one
static void name_ids(int32_t* ida, int32_t* idb, const res_list_t* ra, const res_list_t* rb, md_temp_scope_t temp) {
    const uint32_t n = ra->n + rb->n;
    uint64_t* key = md_temp_alloc_array(temp, uint64_t, n);
    md_memcpy(key, ra->res_name, ra->n * sizeof(uint64_t));
    md_memcpy(key + ra->n, rb->res_name, rb->n * sizeof(uint64_t));
    md_sort_radix_u64(key, n);
    uint32_t u = 0;
    for (uint32_t i = 0; i < n; ++i) if (!u || key[u - 1] != key[i]) key[u++] = key[i];
    for (uint32_t k = 0; k < n; ++k) {
        const uint64_t x = k < ra->n ? ra->res_name[k] : rb->res_name[k - ra->n];
        uint32_t lo = 0, hi = u;
        while (lo < hi) { const uint32_t mid = (lo + hi) / 2; if (key[mid] < x) lo = mid + 1; else hi = mid; }
        if (k < ra->n) ida[k] = (int32_t)lo;
        else idb[k - ra->n] = (int32_t)lo;
    }
}

// The residues aligned by name, in order (Needleman-Wunsch), and the atoms of each aligned pair matched: by name when
// the two are of one name, as graphs when not, once per pair of names. The pairs into pa and pb[ra->beg[ra->n]], how
// many.
static size_t by_residue(int32_t* pa, int32_t* pb, bool* incomplete, const md_mcs_desc_t* desc, const res_list_t* ra, const res_list_t* rb, md_temp_scope_t temp) {
    const md_system_t* sa = desc->a.sys;
    const md_system_t* sb = desc->b.sys;
    const uint32_t na = ra->n, nb = rb->n;
    int32_t* ida = md_temp_alloc_array(temp, int32_t, na);
    int32_t* idb = md_temp_alloc_array(temp, int32_t, (nb + 7) & ~7u);
    name_ids(ida, idb, ra, rb, temp);

    // The scores: a residue of the same name 4, of another -1, a gap -3. A score is the most of the one up and left
    // and the pair's, the one up and a gap, and the one left and a gap: the first two are of the row above, the last
    // a prefix maximum along the row (a gap less for each step). Rows are made eight at a time, padded to a multiple
    // of eight: what is past nb only ever reaches further past it.
    const int32_t gap = -3;
    const uint32_t nb8 = (nb + 7) & ~7u;
    const size_t stride = (size_t)nb8 + 1;
    int32_t* sc = md_temp_alloc_array(temp, int32_t, (na + 1) * stride);
    for (uint32_t j = 0; j <= nb8; ++j) sc[j] = (int32_t)j * gap;
    for (uint32_t j = nb; j < nb8; ++j) idb[j] = -1;
    // Shifts by 1, 2 and 4 lanes, and what is added to the lanes shifted: the gaps of the steps, far below any score for
    // lanes shifted in (md_mm256_set_epi32 takes the lanes highest first)
    const int32_t far = -(1 << 29);
    const md_256i sh1 = md_mm256_set_epi32(6, 5, 4, 3, 2, 1, 0, 0);
    const md_256i sh2 = md_mm256_set_epi32(5, 4, 3, 2, 1, 0, 0, 0);
    const md_256i sh4 = md_mm256_set_epi32(3, 2, 1, 0, 0, 0, 0, 0);
    const md_256i add1 = md_mm256_set_epi32(gap, gap, gap, gap, gap, gap, gap, far);
    const md_256i add2 = md_mm256_set_epi32(2 * gap, 2 * gap, 2 * gap, 2 * gap, 2 * gap, 2 * gap, far, far);
    const md_256i add4 = md_mm256_set_epi32(4 * gap, 4 * gap, 4 * gap, 4 * gap, far, far, far, far);
    const md_256i addc = md_mm256_set_epi32(8 * gap, 7 * gap, 6 * gap, 5 * gap, 4 * gap, 3 * gap, 2 * gap, gap);
    const md_256i last = md_mm256_set1_epi32(7);
    const md_256i vgap = md_mm256_set1_epi32(gap);
    const md_256i miss = md_mm256_set1_epi32(-1);
    const md_256i five = md_mm256_set1_epi32(5);
    for (uint32_t i = 1; i <= na; ++i) {
        const int32_t* up = sc + (i - 1) * stride;
        int32_t* row = sc + i * stride;
        row[0] = (int32_t)i * gap;
        const md_256i x = md_mm256_set1_epi32(ida[i - 1]);
        md_256i carry = md_mm256_set1_epi32(row[0]);
        for (uint32_t j = 1; j <= nb8; j += 8) {
            const md_256i same = md_mm256_cmpeq_epi32(md_mm256_loadu_si256((const md_256i*)(idb + j - 1)), x);
            const md_256i diag = md_mm256_add_epi32(md_mm256_loadu_si256((const md_256i*)(up + j - 1)), md_mm256_add_epi32(miss, md_mm256_and_si256(same, five)));
            const md_256i vert = md_mm256_add_epi32(md_mm256_loadu_si256((const md_256i*)(up + j)), vgap);
            md_256i h = md_mm256_max_epi32(diag, vert);
            h = md_mm256_max_epi32(h, md_mm256_add_epi32(md_mm256_permutevar8x32_epi32(h, sh1), add1));
            h = md_mm256_max_epi32(h, md_mm256_add_epi32(md_mm256_permutevar8x32_epi32(h, sh2), add2));
            h = md_mm256_max_epi32(h, md_mm256_add_epi32(md_mm256_permutevar8x32_epi32(h, sh4), add4));
            h = md_mm256_max_epi32(h, md_mm256_add_epi32(carry, addc));
            md_mm256_storeu_si256((md_256i*)(row + j), h);
            carry = md_mm256_permutevar8x32_epi32(h, last);
        }
    }
    // The aligned pairs, last first: back along the scores, by the first of the three a score is, in that order
    const uint32_t cap = MD_MIN(na, nb);
    uint32_t* ai = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* bi = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t np = 0;
    for (uint32_t i = na, j = nb; i > 0 || j > 0;) {
        const int32_t h = sc[i * stride + j];
        if (i && j && h == sc[(i - 1) * stride + j - 1] + (ida[i - 1] == idb[j - 1] ? 4 : -1)) { ai[np] = --i; bi[np] = --j; np++; }
        else if (i && (!j || h == sc[(i - 1) * stride + j] + gap)) i--;
        else j--;
    }

    // Each residue of a is in one pair at most, and a pair gives at most its heavy atoms of a: the cache of the pairs of
    // names and the pairs themselves are no more than the atoms of a
    typedef struct name_pair_t { uint64_t a, b; uint32_t beg, end; } name_pair_t;
    const size_t max_pairs = ra->beg[na];
    name_pair_t* cache = md_temp_alloc_array(temp, name_pair_t, np);
    uint64_t* cache_a = md_temp_alloc_array(temp, uint64_t, max_pairs);
    uint64_t* cache_b = md_temp_alloc_array(temp, uint64_t, max_pairs);
    uint32_t num_cache = 0, num_cached = 0;
    size_t count = 0;
    const bool by_name = desc->flags & MD_MCS_FLAG_BY_NAME;
    const uint64_t max_steps = desc->max_steps ? desc->max_steps : MCS_DEFAULT_MAX_STEPS;

    for (uint32_t k = np; k-- > 0;) {
        const uint32_t x = ai[k], y = bi[k];
        const uint32_t xb = ra->beg[x], xe = ra->beg[x + 1];
        const uint32_t yb = rb->beg[y], ye = rb->beg[y + 1];
        const uint64_t rna = ra->res_name[x], rnb = rb->res_name[y];
        if (rna == rnb) {
            for (uint32_t p = xb; p < xe; ++p) {
                if (!ra->heavy[p]) continue;
                for (uint32_t q = yb; q < ye; ++q) {
                    if (!rb->heavy[q] || rb->name[q] != ra->name[p]) continue;
                    pa[count] = ra->atoms[p];
                    pb[count] = rb->atoms[q];
                    count++;
                    break;
                }
            }
            continue;
        }
        // Of different names: the graph of the two, remembered as pairs of atom names
        const name_pair_t* hit = NULL;
        for (uint32_t c = 0; c < num_cache && !hit; ++c) if (cache[c].a == rna && cache[c].b == rnb) hit = &cache[c];
        if (!hit) {
            name_pair_t* e = &cache[num_cache++];
            *e = (name_pair_t){.a = rna, .b = rnb, .beg = num_cached};
            // The graphs are scratch of this pair alone
            const md_temp_mark_t mark = md_temp_mark(temp);
            mcs_graph_t g0, g1;
            if (graph_init(&g0, sa, ra->atoms + xb, xe - xb, by_name, temp) && graph_init(&g1, sb, rb->atoms + yb, ye - yb, by_name, temp)) {
                vpair_t* m = md_temp_alloc_array(temp, vpair_t, MD_MIN(g0.n, g1.n));
                bool gave_up = false;
                const int c = graph_mcs(m, &g0, &g1, max_steps, &gave_up, temp);
                *incomplete |= gave_up;
                for (int i = 0; i < c; ++i) {
                    cache_a[num_cached] = name_key(md_atom_name(&sa->atom, g0.atom[m[i].v]));
                    cache_b[num_cached] = name_key(md_atom_name(&sb->atom, g1.atom[m[i].w]));
                    num_cached++;
                }
            }
            md_temp_reset(temp, mark);
            e->end = num_cached;
            hit = e;
        }
        for (uint32_t c = hit->beg; c < hit->end; ++c) {
            const int32_t a = atom_of_name(ra, xb, xe, cache_a[c]);
            const int32_t b = atom_of_name(rb, yb, ye, cache_b[c]);
            if (a < 0 || b < 0) continue;
            pa[count] = a;
            pb[count] = b;
            count++;
        }
    }
    MD_ASSERT(count <= max_pairs);
    return count;
}

// ### THE CALL ###

static size_t count_heavy(const md_mcs_side_t* s) {
    size_t n = 0;
    for (size_t i = 0; i < s->count; ++i) if (in_range(s->sys, s->atoms[i]) && heavy(s->sys, s->atoms[i])) n++;
    return n;
}

bool md_mcs(md_mcs_result_t* out, const md_mcs_desc_t* desc, md_allocator_i* alloc) {
    if (!out) return false;
    md_memset(out, 0, sizeof(*out));
    out->rmsd = NAN;
    if (!desc || !alloc || !desc->a.sys || !desc->b.sys || (desc->a.count && !desc->a.atoms) || (desc->b.count && !desc->b.atoms)) {
        MD_LOG_ERROR("md_mcs: an invalid description");
        return false;
    }
    out->alloc = alloc;
    out->a_heavy = count_heavy(&desc->a);
    out->b_heavy = count_heavy(&desc->b);
    const bool by_name = desc->flags & MD_MCS_FLAG_BY_NAME;
    const bool xyz = desc->a.xyz && desc->b.xyz;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    int32_t* pa = NULL;
    int32_t* pb = NULL;
    size_t n = 0;
    bool ok = true;
    res_list_t ra, rb;
    if (!(desc->flags & MD_MCS_FLAG_GRAPH) && residues_of(&ra, &desc->a, temp) && residues_of(&rb, &desc->b, temp)) {
        pa = md_temp_alloc_array(temp, int32_t, ra.beg[ra.n]);
        pb = md_temp_alloc_array(temp, int32_t, ra.beg[ra.n]);
        n = by_residue(pa, pb, &out->incomplete, desc, &ra, &rb, temp);
        out->by_residue = true;
    } else {
        mcs_graph_t g0, g1;
        if (!graph_init(&g0, desc->a.sys, desc->a.atoms, desc->a.count, by_name, temp) || !graph_init(&g1, desc->b.sys, desc->b.atoms, desc->b.count, by_name, temp)) {
            MD_LOG_ERROR("md_mcs: more than %d heavy atoms is too many to search as one graph", MD_MCS_GRAPH_MAX);
            ok = false;
        } else {
            vpair_t* m = md_temp_alloc_array(temp, vpair_t, MD_MIN(g0.n, g1.n));
            const uint64_t max_steps = desc->max_steps ? desc->max_steps : MCS_DEFAULT_MAX_STEPS;
            n = (size_t)graph_mcs(m, &g0, &g1, max_steps, &out->incomplete, temp);
            pa = md_temp_alloc_array(temp, int32_t, n);
            pb = md_temp_alloc_array(temp, int32_t, n);
            for (size_t i = 0; i < n; ++i) { pa[i] = g0.atom[m[i].v]; pb[i] = g1.atom[m[i].w]; }
            if (n >= 3 && xyz) least_rmsd(desc, pa, pb, n, &g1, by_name, temp);
        }
    }
    if (ok && n) {
        // a's atoms ascending, each with its partner: the two as one key
        uint64_t* key = md_temp_alloc_array(temp, uint64_t, n);
        for (size_t i = 0; i < n; ++i) key[i] = ((uint64_t)(uint32_t)pa[i] << 32) | (uint32_t)pb[i];
        md_sort_radix_u64(key, n);
        out->count = n;
        out->a = md_alloc(alloc, n * sizeof(int32_t));
        out->b = md_alloc(alloc, n * sizeof(int32_t));
        for (size_t i = 0; i < n; ++i) { out->a[i] = (int32_t)(key[i] >> 32); out->b[i] = (int32_t)(uint32_t)key[i]; }
        if (xyz) out->rmsd = pairs_rmsd(desc->a.xyz, out->a, desc->b.xyz, out->b, n, temp);
    }
    md_temp_end(temp);
    return ok;
}

void md_mcs_result_free(md_mcs_result_t* r) {
    if (!r || !r->alloc) return;
    if (r->count) {
        md_free(r->alloc, r->a, r->count * sizeof(int32_t));
        md_free(r->alloc, r->b, r->count * sizeof(int32_t));
    }
    md_memset(r, 0, sizeof(*r));
}

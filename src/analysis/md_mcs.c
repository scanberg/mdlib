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

typedef struct mcs_graph_t {
    int             n;
    int32_t*        atom;       // [n] the atom of each vertex, by degree, highest first
    mcs_label_t*    label;      // [n]
    uint8_t*        adj;        // [n * n]
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

    g->n     = n;
    g->atom  = md_temp_alloc_array(temp, int32_t, n);
    g->label = md_temp_alloc_array(temp, mcs_label_t, n);
    g->adj   = md_temp_alloc_zero(temp, (size_t)n * n);
    for (int i = 0; i < n; ++i) {
        const int p = pos[i];
        g->atom[p]  = (int32_t)h[i];
        g->label[p] = atom_label(sys, (int32_t)h[i], by_name);
        for (uint32_t m = off[i]; m < off[i + 1]; ++m) {
            const int q = pos[nbr[m]];
            g->adj[(size_t)p * n + q] = g->adj[(size_t)q * n + p] = 1;
        }
    }
    return true;
}

// ### McSPLIT ###
// Bidomains: a class of vertices of each graph that may still be matched to each other, alike in label and in their
// adjacency to every pair matched so far. A bound is the matched pairs and the least of each class's two sides.

typedef struct bidomain_t { int l, r, left_len, right_len; bool adjacent; } bidomain_t;
typedef struct vpair_t { int v, w; } vpair_t;

typedef struct search_t {
    const mcs_graph_t*  g0;
    const mcs_graph_t*  g1;
    int*        left;       // [g0->n] the vertices of the classes, each class's a range
    int*        right;      // [g1->n]
    vpair_t*    cur;
    int         cur_n;
    vpair_t*    best;
    int         best_n;
    uint64_t    steps, max_steps;
    bool        gave_up;
    bidomain_t* pool;       // The classes of each level of the search, a stack
    size_t      pool_top, pool_cap;
} search_t;

static int calc_bound(const bidomain_t* d, int nd) {
    int b = 0;
    for (int i = 0; i < nd; ++i) b += MD_MIN(d[i].left_len, d[i].right_len);
    return b;
}

static int min_value(const int* arr, int len) {
    int m = INT_MAX;
    for (int i = 0; i < len; ++i) m = MD_MIN(m, arr[i]);
    return m;
}

// The class with the fewest on its larger side, ties to the one with the lowest vertex; once something is matched,
// only classes adjacent to it, so what is matched stays connected
static int select_bidomain(const search_t* s, const bidomain_t* d, int nd) {
    int min_size = INT_MAX, min_tie = INT_MAX, best = -1;
    for (int i = 0; i < nd; ++i) {
        if (s->cur_n > 0 && !d[i].adjacent) continue;
        const int len = MD_MAX(d[i].left_len, d[i].right_len);
        if (len > min_size) continue;
        const int tie = min_value(s->left + d[i].l, d[i].left_len);
        if (len < min_size || tie < min_tie) { min_size = len; min_tie = tie; best = i; }
    }
    return best;
}

// The vertices of arr[len] adjacent to the one of row first; how many
static int partition(int* arr, int len, const uint8_t* row) {
    int i = 0;
    for (int j = 0; j < len; ++j) {
        if (row[arr[j]]) { const int t = arr[i]; arr[i] = arr[j]; arr[j] = t; ++i; }
    }
    return i;
}

// The classes of d split by adjacency to v and to w, into out; how many
static int filter_domains(search_t* s, const bidomain_t* d, int nd, bidomain_t* out, int v, int w) {
    const uint8_t* rv = s->g0->adj + (size_t)v * s->g0->n;
    const uint8_t* rw = s->g1->adj + (size_t)w * s->g1->n;
    int k = 0;
    for (int i = 0; i < nd; ++i) {
        const bidomain_t* o = &d[i];
        const int ll = partition(s->left + o->l, o->left_len, rv);
        const int rl = partition(s->right + o->r, o->right_len, rw);
        const int ln = o->left_len - ll, rn = o->right_len - rl;
        if (ln && rn) out[k++] = (bidomain_t){o->l + ll, o->r + rl, ln, rn, o->adjacent};
        if (ll && rl) out[k++] = (bidomain_t){o->l, o->r, ll, rl, true};
    }
    return k;
}

// The index of the least of arr[len] above w, -1 when none is
static int next_smallest(const int* arr, int len, int w) {
    int idx = -1, smallest = INT_MAX;
    for (int i = 0; i < len; ++i) if (arr[i] > w && arr[i] < smallest) { smallest = arr[i]; idx = i; }
    return idx;
}

// The lowest vertex v of a class is matched to each of the other side's in turn, lowest first; then the class is
// searched without v, which is the loop (in place of a call in the tail)
static void solve(search_t* s, bidomain_t* d, int nd) {
    for (;;) {
        if (++s->steps > s->max_steps) { s->gave_up = true; return; }
        if (s->cur_n > s->best_n) { md_memcpy(s->best, s->cur, s->cur_n * sizeof(vpair_t)); s->best_n = s->cur_n; }
        if (s->cur_n + calc_bound(d, nd) <= s->best_n) return;
        const int bi = select_bidomain(s, d, nd);
        if (bi < 0) return;
        bidomain_t* bd = &d[bi];

        // v to the end of the left side, out of the class while it is matched
        int* L = s->left + bd->l;
        int vi = 0;
        for (int i = 1; i < bd->left_len; ++i) if (L[i] < L[vi]) vi = i;
        const int v = L[vi];
        L[vi] = L[bd->left_len - 1];
        L[bd->left_len - 1] = v;
        bd->left_len--;

        // Each w to the end of the right side in turn: the classes below are partitions of what is before it
        int* R = s->right + bd->r;
        const int rn = --bd->right_len;
        int w = -1;
        for (int i = 0; i <= rn; ++i) {
            const int idx = next_smallest(R, rn + 1, w);
            w = R[idx];
            R[idx] = R[rn];
            R[rn] = w;
            bidomain_t* below = s->pool + s->pool_top;
            const int k = filter_domains(s, d, nd, below, v, w);
            MD_ASSERT(s->pool_top + k <= s->pool_cap);
            s->pool_top += k;
            s->cur[s->cur_n++] = (vpair_t){v, w};
            solve(s, below, k);
            s->cur_n--;
            s->pool_top -= k;
            if (s->gave_up) return;
        }
        bd->right_len++;
        if (bd->left_len == 0) d[bi] = d[--nd];
    }
}

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

// The largest connected common subgraph of g0 and g1: its pairs into out[min(g0->n, g1->n)] (vertices), how many;
// gave_up when the search ran out of steps
static int graph_mcs(vpair_t* out, const mcs_graph_t* g0, const mcs_graph_t* g1, uint64_t max_steps, bool* gave_up, md_temp_scope_t temp) {
    *gave_up = false;
    const int n0 = g0->n, n1 = g1->n, nmin = MD_MIN(n0, n1);
    if (!nmin) return 0;
    search_t s = {.g0 = g0, .g1 = g1, .max_steps = max_steps, .best = out};
    s.left  = md_temp_alloc_array(temp, int, n0);
    s.right = md_temp_alloc_array(temp, int, n1);
    s.cur   = md_temp_alloc_array(temp, vpair_t, nmin);
    // The classes of a level are disjoint and none is empty on either side: below k pairs there are at most nmin - k,
    // and a path of the search holds at most nmin (nmin + 1) / 2
    s.pool_cap = (size_t)nmin * (nmin + 1) / 2;
    s.pool = md_temp_alloc_array(temp, bidomain_t, s.pool_cap);

    // The first classes: the vertices of both by label, each side's ascending
    const int ne = n0 + n1;
    label_entry_t* e = md_temp_alloc_array(temp, label_entry_t, ne);
    for (int v = 0; v < n0; ++v) e[v]      = (label_entry_t){g0->label[v].name, g0->label[v].z, (uint32_t)v};
    for (int w = 0; w < n1; ++w) e[n0 + w] = (label_entry_t){g1->label[w].name, g1->label[w].z, (uint32_t)w | 0x80000000u};
    qsort(e, ne, sizeof(label_entry_t), cmp_label_entry);
    int nd = 0, nl = 0, nr = 0;
    for (int i = 0, j; i < ne; i = j) {
        const int l0 = nl, r0 = nr;
        for (j = i; j < ne && e[j].z == e[i].z && e[j].name == e[i].name; ++j) {
            if (e[j].id & 0x80000000u) s.right[nr++] = (int)(e[j].id & 0x7FFFFFFFu);
            else s.left[nl++] = (int)e[j].id;
        }
        if (nl > l0 && nr > r0) s.pool[nd++] = (bidomain_t){l0, r0, nl - l0, nr - r0, false};
        else { nl = l0; nr = r0; }
    }
    s.pool_top = nd;
    solve(&s, s.pool, nd);
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

// Called from within md_match_for_each, whose temp scope is open above ours: it allocates nothing
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
    md_match_for_each(NULL, NULL, &md, desc->b.sys, sym_cb, &s);
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

// The residues aligned by name, in order (Needleman-Wunsch), and the atoms of each aligned pair matched: by name when
// the two are of one name, as graphs when not, once per pair of names. The pairs into pa and pb[ra->beg[ra->n]], how
// many.
static size_t by_residue(int32_t* pa, int32_t* pb, bool* incomplete, const md_mcs_desc_t* desc, const res_list_t* ra, const res_list_t* rb, md_temp_scope_t temp) {
    const md_system_t* sa = desc->a.sys;
    const md_system_t* sb = desc->b.sys;
    const uint32_t na = ra->n, nb = rb->n;
    const size_t stride = (size_t)nb + 1;
    int32_t* sc = md_temp_alloc_array(temp, int32_t, (na + 1) * stride);
    uint8_t* tb = md_temp_alloc_array(temp, uint8_t, (na + 1) * stride);
    const int32_t gap = -3;
    for (uint32_t i = 0; i <= na; ++i) { sc[i * stride] = (int32_t)i * gap; tb[i * stride] = 1; }
    for (uint32_t j = 0; j <= nb; ++j) { sc[j] = (int32_t)j * gap; tb[j] = 2; }
    for (uint32_t i = 1; i <= na; ++i) {
        const uint64_t ni = ra->res_name[i - 1];
        const int32_t* up = sc + (i - 1) * stride;
        int32_t* row = sc + i * stride;
        uint8_t* trow = tb + i * stride;
        for (uint32_t j = 1; j <= nb; ++j) {
            int32_t best = up[j - 1] + (ni == rb->res_name[j - 1] ? 4 : -1);
            uint8_t t = 0;
            if (up[j] + gap > best)      { best = up[j] + gap; t = 1; }
            if (row[j - 1] + gap > best) { best = row[j - 1] + gap; t = 2; }
            row[j] = best;
            trow[j] = t;
        }
    }
    // The aligned pairs, last first
    const uint32_t cap = MD_MIN(na, nb);
    uint32_t* ai = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* bi = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t np = 0;
    for (uint32_t i = na, j = nb; i > 0 || j > 0;) {
        const uint8_t t = (i > 0 && j > 0) ? tb[i * stride + j] : (i > 0 ? 1 : 2);
        if (t == 0) { ai[np] = --i; bi[np] = --j; np++; }
        else if (t == 1) i--;
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

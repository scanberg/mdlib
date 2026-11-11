// Maximum common substructure (md_mcs.h)

#include <md_mcs.h>
#include <md_system.h>
#include <md_match.h>
#include <md_util.h>
#include <md_util_geom.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ### GRAPHS ###

typedef struct mcs_graph_t {
    int         n;
    int32_t*    atom;       // [n] the atom of each vertex, by degree, highest first
    uint32_t*   label;      // [n]
    uint8_t*    adj;        // [n * n]
} mcs_graph_t;

static uint32_t atom_label(const md_system_t* sys, int32_t a, bool by_name) {
    const uint32_t z = md_atom_atomic_number(&sys->atom, a);
    if (!by_name) return z;
    const md_str_t name = md_atom_name(&sys->atom, a);
    return (z << 20) | ((uint32_t)md_hash64(name.ptr, name.len, 0) & 0xFFFFF);
}

static bool heavy(const md_system_t* sys, int32_t a) {
    return md_atom_atomic_number(&sys->atom, a) != 1;
}

// The heavy atoms of atoms[count] and their bonds, vertices by their degree within the set, highest first: McSplit
// branches on the lowest vertex
static bool graph_init(mcs_graph_t* g, const md_system_t* sys, const int32_t* atoms, size_t count, bool by_name, md_allocator_i* alloc) {
    memset(g, 0, sizeof(*g));
    int32_t* h = md_alloc(alloc, (count ? count : 1) * sizeof(int32_t));
    for (size_t i = 0; i < count; ++i) if (atoms[i] >= 0 && (size_t)atoms[i] < sys->atom.count && heavy(sys, atoms[i])) h[g->n++] = atoms[i];
    if (g->n > MD_MCS_GRAPH_MAX) return false;
    const int n = g->n;
    // Sorted by atom, for finding a neighbour's vertex
    for (int i = 1; i < n; ++i) { const int32_t x = h[i]; int j = i; while (j > 0 && h[j - 1] > x) { h[j] = h[j - 1]; --j; } h[j] = x; }
    int* deg = md_alloc(alloc, (n ? n : 1) * sizeof(int));
    memset(deg, 0, (n ? n : 1) * sizeof(int));
    uint8_t* adj0 = md_alloc(alloc, (size_t)(n ? n : 1) * (n ? n : 1));
    memset(adj0, 0, (size_t)(n ? n : 1) * (n ? n : 1));
    for (int i = 0; i < n; ++i) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, h[i]);
        while (md_bond_iter_has_next(&it)) {
            const int32_t o = (int32_t)md_bond_iter_atom_index(&it);
            md_bond_iter_next(&it);
            int lo = 0, hi = n;
            while (lo < hi) { const int mid = (lo + hi) / 2; if (h[mid] < o) lo = mid + 1; else hi = mid; }
            if (lo < n && h[lo] == o && lo != i && !adj0[(size_t)i * n + lo]) {
                adj0[(size_t)i * n + lo] = adj0[(size_t)lo * n + i] = 1;
                deg[i]++; deg[lo]++;
            }
        }
    }
    int* order = md_alloc(alloc, (n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; ++i) order[i] = i;
    for (int i = 1; i < n; ++i) { const int x = order[i]; int j = i; while (j > 0 && deg[order[j - 1]] < deg[x]) { order[j] = order[j - 1]; --j; } order[j] = x; }
    int* pos = md_alloc(alloc, (n ? n : 1) * sizeof(int));
    for (int i = 0; i < n; ++i) pos[order[i]] = i;
    g->atom = md_alloc(alloc, (n ? n : 1) * sizeof(int32_t));
    g->label = md_alloc(alloc, (n ? n : 1) * sizeof(uint32_t));
    g->adj = md_alloc(alloc, (size_t)(n ? n : 1) * (n ? n : 1));
    memset(g->adj, 0, (size_t)(n ? n : 1) * (n ? n : 1));
    for (int i = 0; i < n; ++i) {
        g->atom[i] = h[order[i]];
        g->label[i] = atom_label(sys, g->atom[i], by_name);
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) if (adj0[(size_t)i * n + j]) g->adj[(size_t)pos[i] * n + pos[j]] = 1;
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
    int*        left;
    int*        right;
    vpair_t*    cur;
    int         cur_n;
    vpair_t*    best;
    int         best_n;
    uint64_t    steps, max_steps;
    bool        gave_up;
    int         goal;
    bidomain_t* pool;
    size_t      pool_cap, pool_top, per_level;
} search_t;

static int calc_bound(const bidomain_t* d, int nd) {
    int b = 0;
    for (int i = 0; i < nd; ++i) b += d[i].left_len < d[i].right_len ? d[i].left_len : d[i].right_len;
    return b;
}

static int min_value(const int* arr, int start, int len) {
    int m = INT_MAX;
    for (int i = 0; i < len; ++i) if (arr[start + i] < m) m = arr[start + i];
    return m;
}

// The class with the fewest on its larger side, ties to the one with the lowest vertex; once something is matched,
// only classes adjacent to it, so what is matched stays connected
static int select_bidomain(const search_t* s, const bidomain_t* d, int nd) {
    int min_size = INT_MAX, min_tie = INT_MAX, best = -1;
    for (int i = 0; i < nd; ++i) {
        if (s->cur_n > 0 && !d[i].adjacent) continue;
        const int len = d[i].left_len > d[i].right_len ? d[i].left_len : d[i].right_len;
        const int tie = min_value(s->left, d[i].l, d[i].left_len);
        if (len < min_size || (len == min_size && tie < min_tie)) { min_size = len; min_tie = tie; best = i; }
    }
    return best;
}

static int partition(int* arr, int start, int len, const uint8_t* row) {
    int i = 0;
    for (int j = 0; j < len; ++j) {
        if (row[arr[start + j]]) { const int t = arr[start + i]; arr[start + i] = arr[start + j]; arr[start + j] = t; ++i; }
    }
    return i;
}

static int filter_domains(search_t* s, const bidomain_t* d, int nd, bidomain_t* out, int v, int w) {
    const uint8_t* rv = s->g0->adj + (size_t)v * s->g0->n;
    const uint8_t* rw = s->g1->adj + (size_t)w * s->g1->n;
    int k = 0;
    for (int i = 0; i < nd; ++i) {
        const bidomain_t* o = &d[i];
        const int ll = partition(s->left, o->l, o->left_len, rv);
        const int rl = partition(s->right, o->r, o->right_len, rw);
        const int ln = o->left_len - ll, rn = o->right_len - rl;
        if (ln && rn) out[k++] = (bidomain_t){o->l + ll, o->r + rl, ln, rn, o->adjacent};
        if (ll && rl) out[k++] = (bidomain_t){o->l, o->r, ll, rl, true};
    }
    return k;
}

static int next_smallest(const int* arr, int start, int len, int w) {
    int idx = -1, smallest = INT_MAX;
    for (int i = 0; i < len; ++i) if (arr[start + i] > w && arr[start + i] < smallest) { smallest = arr[start + i]; idx = i; }
    return idx;
}

static void solve(search_t* s, bidomain_t* d, int nd) {
    if (s->gave_up) return;
    if (++s->steps > s->max_steps) { s->gave_up = true; return; }
    if (s->cur_n > s->best_n) { memcpy(s->best, s->cur, s->cur_n * sizeof(vpair_t)); s->best_n = s->cur_n; }
    const int bound = s->cur_n + calc_bound(d, nd);
    if (bound <= s->best_n || bound < s->goal) return;
    const int bi = select_bidomain(s, d, nd);
    if (bi < 0) return;
    bidomain_t* bd = &d[bi];
    const int v = min_value(s->left, bd->l, bd->left_len);
    for (int i = 0; i < bd->left_len; ++i) {
        if (s->left[bd->l + i] != v) continue;
        s->left[bd->l + i] = s->left[bd->l + bd->left_len - 1];
        s->left[bd->l + bd->left_len - 1] = v;
        break;
    }
    bd->left_len--;
    int w = -1;
    bd->right_len--;
    for (int i = 0; i <= bd->right_len; ++i) {
        const int idx = next_smallest(s->right, bd->r, bd->right_len + 1, w);
        w = s->right[bd->r + idx];
        s->right[bd->r + idx] = s->right[bd->r + bd->right_len];
        s->right[bd->r + bd->right_len] = w;
        const size_t top = s->pool_top;
        MD_ASSERT(top + s->per_level <= s->pool_cap);
        bidomain_t* nd2 = s->pool + top;
        s->pool_top += s->per_level;
        const int k = filter_domains(s, d, nd, nd2, v, w);
        s->cur[s->cur_n++] = (vpair_t){v, w};
        solve(s, nd2, k);
        s->cur_n--;
        s->pool_top = top;
        if (s->gave_up) break;
    }
    bd->right_len++;
    if (bd->left_len == 0) { d[bi] = d[nd - 1]; nd--; }
    solve(s, d, nd);
}

static int cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

// The largest connected common subgraph of g0 and g1: its pairs into out (vertices), how many; gave_up when the search
// ran out of steps
static int graph_mcs(vpair_t* out, const mcs_graph_t* g0, const mcs_graph_t* g1, uint64_t max_steps, bool* gave_up, md_allocator_i* alloc) {
    *gave_up = false;
    if (!g0->n || !g1->n) return 0;
    search_t s = {.g0 = g0, .g1 = g1, .max_steps = max_steps};
    s.left = md_alloc(alloc, (g0->n + 1) * sizeof(int));
    s.right = md_alloc(alloc, (g1->n + 1) * sizeof(int));
    s.cur = md_alloc(alloc, (g0->n + 1) * sizeof(vpair_t));
    s.best = out;
    // At most as many classes as unmatched vertices of g0, a level for each pair matched
    s.per_level = (size_t)g0->n + 2;
    const int nmin = g0->n < g1->n ? g0->n : g1->n;
    s.pool_cap = s.per_level * ((size_t)nmin + 2);
    s.pool = md_alloc(alloc, s.pool_cap * sizeof(bidomain_t));
    // The first classes: by label, the labels of both
    uint32_t* labels = md_alloc(alloc, (g0->n + 1) * sizeof(uint32_t));
    memcpy(labels, g0->label, g0->n * sizeof(uint32_t));
    qsort(labels, g0->n, sizeof(uint32_t), cmp_u32);
    bidomain_t* d0 = s.pool;
    s.pool_top = s.per_level;
    int nd0 = 0, nl = 0, nr = 0;
    for (int i = 0; i < g0->n; ++i) {
        if (i && labels[i] == labels[i - 1]) continue;
        const uint32_t lab = labels[i];
        const int l0 = nl, r0 = nr;
        for (int v = 0; v < g0->n; ++v) if (g0->label[v] == lab) s.left[nl++] = v;
        for (int w = 0; w < g1->n; ++w) if (g1->label[w] == lab) s.right[nr++] = w;
        if (nr > r0) d0[nd0++] = (bidomain_t){l0, r0, nl - l0, nr - r0, false};
        else nl = l0;
    }
    solve(&s, d0, nd0);
    *gave_up = s.gave_up;
    return s.best_n;
}

// ### SUPERPOSITION ###

// RMSD of b's positions superposed onto a's, equal weights
static double pairs_rmsd(const md_f32x3_t* xa, const int32_t* a, const md_f32x3_t* xb, const int32_t* b, size_t n, md_allocator_i* alloc) {
    if (n < 3) return NAN;
    md_f32x4_t* pa = md_alloc(alloc, n * sizeof(md_f32x4_t));
    md_f32x4_t* pb = md_alloc(alloc, n * sizeof(md_f32x4_t));
    double ca[3] = {0}, cb[3] = {0};
    for (size_t i = 0; i < n; ++i) for (int k = 0; k < 3; ++k) { ca[k] += xa[a[i]].elem[k]; cb[k] += xb[b[i]].elem[k]; }
    for (int k = 0; k < 3; ++k) { ca[k] /= n; cb[k] /= n; }
    for (size_t i = 0; i < n; ++i) {
        pa[i] = md_f32x4_set((float)(xa[a[i]].x - ca[0]), (float)(xa[a[i]].y - ca[1]), (float)(xa[a[i]].z - ca[2]), 1.0f);
        pb[i] = md_f32x4_set((float)(xb[b[i]].x - cb[0]), (float)(xb[b[i]].y - cb[1]), (float)(xb[b[i]].z - cb[2]), 1.0f);
    }
    const md_f32x3x3_t R = md_util_optimal_rotation_rel_vec4(pa, pb, n);
    double ss = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const md_f32x3_t q = md_f32x3x3_mul_f32x3(R, md_f32x3_from_f32x4(pb[i]));
        for (int k = 0; k < 3; ++k) { const double d = (double)q.elem[k] - pa[i].elem[k]; ss += d * d; }
    }
    return sqrt(ss / n);
}

// Among the mappings of a's common atoms onto b's atoms that the bonds allow, the one of the least RMSD
typedef struct sym_t {
    const md_f32x3_t* xa;
    const md_f32x3_t* xb;
    const int32_t* a;
    int32_t* best;
    size_t n;
    double best_rmsd;
    md_allocator_i* temp;
} sym_t;

static bool sym_cb(const int32_t* atom_idx, size_t width, uint32_t unit, void* user) {
    (void)unit;
    sym_t* s = user;
    if (width != s->n) return true;
    for (size_t i = 0; i < width; ++i) if (atom_idx[i] < 0) return true;
    const double r = pairs_rmsd(s->xa, s->a, s->xb, atom_idx, width, s->temp);
    if (r < s->best_rmsd) { s->best_rmsd = r; memcpy(s->best, atom_idx, width * sizeof(int32_t)); }
    return true;
}

static void least_rmsd(const md_mcs_desc_t* desc, const int32_t* a, int32_t* b, size_t n, const mcs_graph_t* g1, md_allocator_i* temp) {
    md_match_query_t q;
    if (!md_match_query_init_atoms(&q, a, n, MD_MATCH_LABEL_ELEMENT, desc->a.sys, temp)) return;
    uint32_t* ua = md_alloc(temp, (g1->n ? g1->n : 1) * sizeof(uint32_t));
    for (int i = 0; i < g1->n; ++i) ua[i] = (uint32_t)g1->atom[i];
    qsort(ua, g1->n, sizeof(uint32_t), cmp_u32);
    md_bitset_t unit = md_bitset_from_indices(ua, g1->n, temp);
    sym_t s = {.xa = desc->a.xyz, .xb = desc->b.xyz, .a = a, .n = n, .best_rmsd = pairs_rmsd(desc->a.xyz, a, desc->b.xyz, b, n, temp), .temp = temp};
    s.best = md_alloc(temp, n * sizeof(int32_t));
    memcpy(s.best, b, n * sizeof(int32_t));
    const md_match_desc_t md = {.query = &q, .mode = MD_MATCH_MODE_ALL, .hydrogens = MD_MATCH_RESOLVE_NEVER, .bond_orders = MD_MATCH_RESOLVE_NEVER,
                                .charges = MD_MATCH_RESOLVE_NEVER, .units = &unit, .num_units = 1, .limit = 4096, .max_steps = 1000000};
    md_match_for_each(NULL, NULL, &md, desc->b.sys, sym_cb, &s);
    memcpy(b, s.best, n * sizeof(int32_t));
    md_match_query_free(&q);
}

// ### BY RESIDUE ###

typedef struct res_list_t {
    uint32_t*   res;        // the residues the atoms are in, in order
    uint32_t*   beg;        // [n + 1] their atoms of the set, atoms[beg[i], beg[i + 1]), ascending
    int32_t*    atoms;
    uint32_t    n;
} res_list_t;

static bool residues_of(res_list_t* out, const md_mcs_side_t* s, md_allocator_i* alloc) {
    memset(out, 0, sizeof(*out));
    const md_system_t* sys = s->sys;
    if (!sys->residue.count || !sys->residue.atom_offset) return false;
    int32_t* at = md_alloc(alloc, (s->count ? s->count : 1) * sizeof(int32_t));
    size_t n = 0;
    for (size_t i = 0; i < s->count; ++i) if (s->atoms[i] >= 0 && (size_t)s->atoms[i] < sys->atom.count) at[n++] = s->atoms[i];
    qsort(at, n, sizeof(int32_t), cmp_u32);
    out->atoms = at;
    md_array(uint32_t) res = NULL;
    md_array(uint32_t) beg = NULL;
    uint32_t r = 0;
    for (size_t i = 0; i < n; ++i) {
        while (r + 1 < sys->residue.count && sys->residue.atom_offset[r + 1] <= (uint32_t)at[i]) ++r;
        if (!(sys->residue.atom_offset[r] <= (uint32_t)at[i] && (uint32_t)at[i] < sys->residue.atom_offset[r + 1])) return false;
        if (!md_array_size(res) || *md_array_last(res) != r) {
            const md_residue_kind_t k = md_system_residue_kind(sys, r);
            if (k != MD_RESIDUE_KIND_AMINO_ACID && k != MD_RESIDUE_KIND_NUCLEOTIDE) return false;
            md_array_push(res, r, alloc);
            md_array_push(beg, (uint32_t)i, alloc);
        }
    }
    md_array_push(beg, (uint32_t)n, alloc);
    out->res = res;
    out->beg = beg;
    out->n = (uint32_t)md_array_size(res);
    return out->n >= 2;
}

typedef struct name_pair_t {
    uint64_t    key;        // the two residue names
    uint32_t    beg, end;   // into the cache's pairs of atom names
} name_pair_t;

static uint64_t name_key(md_str_t a, md_str_t b) {
    return md_hash64(b.ptr, b.len, md_hash64(a.ptr, a.len, 0x9e3779b9));
}

static bool by_residue(md_mcs_result_t* out, const md_mcs_desc_t* desc, const res_list_t* ra, const res_list_t* rb, md_allocator_i* temp,
                       md_array(int32_t)* pa, md_array(int32_t)* pb) {
    const md_system_t* sa = desc->a.sys;
    const md_system_t* sb = desc->b.sys;
    const uint32_t na = ra->n, nb = rb->n;
    // Needleman-Wunsch over the residue names
    int32_t* sc = md_alloc(temp, (size_t)(na + 1) * (nb + 1) * sizeof(int32_t));
    uint8_t* tb = md_alloc(temp, (size_t)(na + 1) * (nb + 1));
    const int gap = -3;
#define SC(i, j) sc[(size_t)(i) * (nb + 1) + (j)]
#define TB(i, j) tb[(size_t)(i) * (nb + 1) + (j)]
    for (uint32_t i = 0; i <= na; ++i) { SC(i, 0) = (int32_t)i * gap; TB(i, 0) = 1; }
    for (uint32_t j = 0; j <= nb; ++j) { SC(0, j) = (int32_t)j * gap; TB(0, j) = 2; }
    for (uint32_t i = 1; i <= na; ++i) {
        const md_str_t ni = md_system_residue_name(sa, ra->res[i - 1]);
        for (uint32_t j = 1; j <= nb; ++j) {
            const bool same = md_str_eq(ni, md_system_residue_name(sb, rb->res[j - 1]));
            int32_t best = SC(i - 1, j - 1) + (same ? 4 : -1);
            uint8_t t = 0;
            if (SC(i - 1, j) + gap > best) { best = SC(i - 1, j) + gap; t = 1; }
            if (SC(i, j - 1) + gap > best) { best = SC(i, j - 1) + gap; t = 2; }
            SC(i, j) = best;
            TB(i, j) = t;
        }
    }
    uint32_t* ai = md_alloc(temp, (na + nb + 1) * sizeof(uint32_t));
    uint32_t* bi = md_alloc(temp, (na + nb + 1) * sizeof(uint32_t));
    uint32_t np = 0;
    for (uint32_t i = na, j = nb; i > 0 || j > 0;) {
        const uint8_t t = (i > 0 && j > 0) ? TB(i, j) : (i > 0 ? 1 : 2);
        if (t == 0) { ai[np] = i - 1; bi[np] = j - 1; np++; i--; j--; }
        else if (t == 1) i--;
        else j--;
    }
#undef SC
#undef TB
    // Each aligned pair: by atom name when of one name, by graph when not, once per pair of names
    const bool by_name = desc->flags & MD_MCS_FLAG_BY_NAME;
    md_array(name_pair_t) cache = NULL;
    md_array(md_label_t) cache_a = NULL;
    md_array(md_label_t) cache_b = NULL;
    for (uint32_t k = np; k-- > 0;) {
        const uint32_t x = ai[k], y = bi[k];
        const int32_t* xa = ra->atoms + ra->beg[x];
        const uint32_t xn = ra->beg[x + 1] - ra->beg[x];
        const int32_t* yb = rb->atoms + rb->beg[y];
        const uint32_t yn = rb->beg[y + 1] - rb->beg[y];
        const md_str_t rna = md_system_residue_name(sa, ra->res[x]), rnb = md_system_residue_name(sb, rb->res[y]);
        if (md_str_eq(rna, rnb)) {
            for (uint32_t p = 0; p < xn; ++p) {
                if (!heavy(sa, xa[p])) continue;
                const md_str_t an = md_atom_name(&sa->atom, xa[p]);
                for (uint32_t q = 0; q < yn; ++q) {
                    if (!heavy(sb, yb[q]) || !md_str_eq(an, md_atom_name(&sb->atom, yb[q]))) continue;
                    md_array_push(*pa, xa[p], temp);
                    md_array_push(*pb, yb[q], temp);
                    break;
                }
            }
            continue;
        }
        // Of different names: the graph of the two, remembered as pairs of atom names
        const uint64_t key = name_key(rna, rnb);
        const name_pair_t* hit = NULL;
        for (size_t c = 0; c < md_array_size(cache) && !hit; ++c) if (cache[c].key == key) hit = &cache[c];
        if (!hit) {
            mcs_graph_t g0, g1;
            name_pair_t e = {.key = key, .beg = (uint32_t)md_array_size(cache_a)};
            if (graph_init(&g0, sa, xa, xn, by_name, temp) && graph_init(&g1, sb, yb, yn, by_name, temp)) {
                vpair_t* m = md_alloc(temp, (size_t)(MD_MIN(g0.n, g1.n) + 1) * sizeof(vpair_t));
                bool gave_up = false;
                const int c = graph_mcs(m, &g0, &g1, desc->max_steps ? desc->max_steps : 10000000, &gave_up, temp);
                if (gave_up) out->incomplete = true;
                for (int i = 0; i < c; ++i) {
                    md_array_push(cache_a, md_label_make(md_atom_name(&sa->atom, g0.atom[m[i].v])), temp);
                    md_array_push(cache_b, md_label_make(md_atom_name(&sb->atom, g1.atom[m[i].w])), temp);
                }
            }
            e.end = (uint32_t)md_array_size(cache_a);
            md_array_push(cache, e, temp);
            hit = md_array_last(cache);
        }
        for (uint32_t c = hit->beg; c < hit->end; ++c) {
            int32_t pa_atom = -1, pb_atom = -1;
            for (uint32_t p = 0; p < xn && pa_atom < 0; ++p) if (md_str_eq(md_atom_name(&sa->atom, xa[p]), MD_LABEL_STR(cache_a[c]))) pa_atom = xa[p];
            for (uint32_t q = 0; q < yn && pb_atom < 0; ++q) if (md_str_eq(md_atom_name(&sb->atom, yb[q]), MD_LABEL_STR(cache_b[c]))) pb_atom = yb[q];
            if (pa_atom < 0 || pb_atom < 0) continue;
            md_array_push(*pa, pa_atom, temp);
            md_array_push(*pb, pb_atom, temp);
        }
    }
    out->by_residue = true;
    return true;
}

// ### THE CALL ###

typedef struct pair_sort_t { int32_t a, b; } pair_sort_t;

static int cmp_pair(const void* x, const void* y) {
    const int32_t a = ((const pair_sort_t*)x)->a, b = ((const pair_sort_t*)y)->a;
    return (a > b) - (a < b);
}

static size_t count_heavy(const md_mcs_side_t* s) {
    size_t n = 0;
    for (size_t i = 0; i < s->count; ++i) if (s->atoms[i] >= 0 && (size_t)s->atoms[i] < s->sys->atom.count && heavy(s->sys, s->atoms[i])) n++;
    return n;
}

bool md_mcs(md_mcs_result_t* out, const md_mcs_desc_t* desc, md_allocator_i* alloc) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    out->rmsd = NAN;
    if (!desc || !alloc || !desc->a.sys || !desc->b.sys || (desc->a.count && !desc->a.atoms) || (desc->b.count && !desc->b.atoms)) {
        MD_LOG_ERROR("md_mcs: an invalid description");
        return false;
    }
    out->alloc = alloc;
    out->a_heavy = count_heavy(&desc->a);
    out->b_heavy = count_heavy(&desc->b);
    md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_array(int32_t) pa = NULL;
    md_array(int32_t) pb = NULL;
    bool ok = true;
    res_list_t ra, rb;
    const bool residues = !(desc->flags & MD_MCS_FLAG_GRAPH) && residues_of(&ra, &desc->a, temp) && residues_of(&rb, &desc->b, temp);
    if (residues) {
        ok = by_residue(out, desc, &ra, &rb, temp, &pa, &pb);
    } else {
        mcs_graph_t g0, g1;
        const bool by_name = desc->flags & MD_MCS_FLAG_BY_NAME;
        if (!graph_init(&g0, desc->a.sys, desc->a.atoms, desc->a.count, by_name, temp) || !graph_init(&g1, desc->b.sys, desc->b.atoms, desc->b.count, by_name, temp)) {
            MD_LOG_ERROR("md_mcs: more than %d heavy atoms is too many to search as one graph", MD_MCS_GRAPH_MAX);
            ok = false;
        } else {
            vpair_t* m = md_alloc(temp, (size_t)(MD_MIN(g0.n, g1.n) + 1) * sizeof(vpair_t));
            bool gave_up = false;
            const int c = graph_mcs(m, &g0, &g1, desc->max_steps ? desc->max_steps : 10000000, &gave_up, temp);
            out->incomplete = gave_up;
            for (int i = 0; i < c; ++i) {
                md_array_push(pa, g0.atom[m[i].v], temp);
                md_array_push(pb, g1.atom[m[i].w], temp);
            }
            if (c >= 3 && desc->a.xyz && desc->b.xyz) least_rmsd(desc, pa, pb, (size_t)c, &g1, temp);
        }
    }
    if (ok) {
        const size_t n = md_array_size(pa);
        pair_sort_t* ps = md_alloc(temp, (n ? n : 1) * sizeof(pair_sort_t));
        for (size_t i = 0; i < n; ++i) ps[i] = (pair_sort_t){pa[i], pb[i]};
        qsort(ps, n, sizeof(pair_sort_t), cmp_pair);
        out->count = n;
        out->a = md_alloc(alloc, (n ? n : 1) * sizeof(int32_t));
        out->b = md_alloc(alloc, (n ? n : 1) * sizeof(int32_t));
        for (size_t i = 0; i < n; ++i) { out->a[i] = ps[i].a; out->b[i] = ps[i].b; }
        if (desc->a.xyz && desc->b.xyz) out->rmsd = pairs_rmsd(desc->a.xyz, out->a, desc->b.xyz, out->b, n, temp);
    }
    md_arena_allocator_destroy(temp);
    return ok;
}

void md_mcs_result_free(md_mcs_result_t* r) {
    if (!r || !r->alloc) return;
    md_free(r->alloc, r->a, (r->count ? r->count : 1) * sizeof(int32_t));
    md_free(r->alloc, r->b, (r->count ? r->count : 1) * sizeof(int32_t));
    memset(r, 0, sizeof(*r));
}

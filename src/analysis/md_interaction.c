#include <md_interaction.h>
#include <md_hbond.h>
#include <md_system.h>
#include <md_element.h>
#include <core/md_geom.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_log.h>
#include <core/md_spatial_acc.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <stdlib.h>

md_interaction_params_t md_interaction_params_default(void) {
    return (md_interaction_params_t){
        .kinds             = MD_INTERACTION_IONIC | MD_INTERACTION_HALOGEN,
        .exclude_bonds     = 3,
        .ionic_max_dist    = 4.0f,
        .halogen_tolerance = 0.0f,
        .halogen_min_cxa   = 140.0f,
    };
}

static int compare_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return x < y ? -1 : x > y;
}

// The root of a union-find forest over the atoms, with path halving
static uint32_t uf_root(uint32_t* parent, uint32_t i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
}

static bool is_contact_element(md_atomic_number_t z) {
    return z == MD_Z_N || z == MD_Z_O || z == MD_Z_S;
}

// Row k of the exclusions: the atoms within max_bonds of atom, sorted
static void push_exclusions(md_interaction_query_t* q, const md_system_t* sys, uint32_t atom, uint32_t max_bonds, uint8_t* depth,
                            md_allocator_i* temp) {
    md_array(uint32_t) queue = 0;
    const size_t row = md_array_size(q->excl_atom);
    depth[atom] = 1;
    md_array_push(queue, atom, temp);
    for (size_t head = 0; head < md_array_size(queue); ++head) {
        const uint32_t cur = queue[head];
        if (depth[cur] > max_bonds) continue;
        md_bond_iter_t it = md_bond_iter(&sys->bond, cur);
        while (md_bond_iter_has_next(&it)) {
            const uint32_t next = (uint32_t)md_bond_iter_atom_index(&it);
            md_bond_iter_next(&it);
            if (next >= q->num_atoms || depth[next]) continue;
            depth[next] = depth[cur] + 1;
            md_array_push(queue, next, temp);
            md_array_push(q->excl_atom, next, q->alloc);
        }
    }
    for (size_t v = 0; v < md_array_size(queue); ++v) depth[queue[v]] = 0;
    const size_t n = md_array_size(q->excl_atom) - row;
    if (n > 1) qsort(q->excl_atom + row, n, sizeof(uint32_t), compare_u32);
    md_array_push(q->excl_off, (uint32_t)md_array_size(q->excl_atom), q->alloc);
}

static bool is_excluded(const md_interaction_query_t* q, size_t row, uint32_t atom) {
    if (!q->excl_off) return false;
    size_t lo = q->excl_off[row], hi = q->excl_off[row + 1];
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (q->excl_atom[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    return lo < q->excl_off[row + 1] && q->excl_atom[lo] == atom;
}

void md_interaction_query_free(md_interaction_query_t* q) {
    if (!q || !q->alloc) return;
    md_allocator_i* alloc = q->alloc;
    md_array_free(q->charged_atom, alloc);
    md_array_free(q->charged_group, alloc);
    md_array_free(q->charged_sign, alloc);
    md_array_free(q->halogen, alloc);
    md_array_free(q->halogen_c, alloc);
    md_array_free(q->halogen_vdw, alloc);
    md_array_free(q->acceptor, alloc);
    md_array_free(q->acceptor_vdw, alloc);
    md_array_free(q->excl_off, alloc);
    md_array_free(q->excl_atom, alloc);
    if (q->sel) md_free(alloc, q->sel, q->num_atoms);
    md_memset(q, 0, sizeof(*q));
}

bool md_interaction_query_init(md_interaction_query_t* q, const md_interaction_desc_t* desc, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(q);
    md_memset(q, 0, sizeof(*q));
    if (!sys || !alloc) {
        MD_LOG_ERROR("Interactions: missing system or allocator");
        return false;
    }
    q->params = desc && desc->params ? *desc->params : md_interaction_params_default();
    q->num_atoms = (uint32_t)sys->atom.count;
    q->alloc = alloc;
    const uint32_t N = q->num_atoms;
    if (desc && desc->set) {
        q->sel = md_alloc(alloc, N ? N : 1);
        for (uint32_t i = 0; i < N; ++i) q->sel[i] = md_bitset_test(*desc->set, i) ? 1 : 0;
    }
    if (N == 0) return true;

    md_allocator_i* temp = md_vm_arena_create(MD_GIGABYTES(1));
    const md_bond_data_t* bd = &sys->bond;

    if (q->params.kinds & MD_INTERACTION_IONIC) {
        // The groups: the atoms joined by delocalized bonds of the graph, a component each
        uint32_t* parent = md_alloc(temp, N * sizeof(uint32_t));
        for (uint32_t i = 0; i < N; ++i) parent[i] = i;
        for (size_t b = 0; b < bd->count; ++b) {
            const md_bond_flags_t f = bd->flags ? bd->flags[b] : 0;
            if (!md_bond_is_graph(f) || (f & MD_BOND_FLAG_RESONANCE_MASK) != MD_BOND_FLAG_DELOCALIZED) continue;
            const md_atom_pair_t pr = bd->pairs[b];
            if (pr.idx[0] < 0 || pr.idx[1] < 0 || (uint32_t)pr.idx[0] >= N || (uint32_t)pr.idx[1] >= N) continue;
            const uint32_t a = uf_root(parent, (uint32_t)pr.idx[0]), c = uf_root(parent, (uint32_t)pr.idx[1]);
            if (a != c) parent[a] = c;
        }
        int32_t* net = md_alloc(temp, N * sizeof(int32_t));
        uint32_t* size = md_alloc(temp, N * sizeof(uint32_t));
        md_memset(net, 0, N * sizeof(int32_t));
        md_memset(size, 0, N * sizeof(uint32_t));
        for (uint32_t i = 0; i < N; ++i) {
            const uint32_t r = uf_root(parent, i);
            net[r] += md_atom_formal_charge(&sys->atom, i);
            size[r] += 1;
        }
        // Each group a number, in the order of its first atom; its sign counted once
        uint32_t* group = md_alloc(temp, N * sizeof(uint32_t));
        for (uint32_t i = 0; i < N; ++i) group[i] = UINT32_MAX;
        uint32_t num_groups = 0;
        for (uint32_t i = 0; i < N; ++i) {
            const uint32_t r = uf_root(parent, i);
            if (!net[r]) continue;
            if (size[r] > 1 && !is_contact_element(md_atom_atomic_number(&sys->atom, i))) continue;
            if (group[r] == UINT32_MAX) {
                group[r] = num_groups++;
                if (net[r] > 0) q->num_cations += 1;
                else q->num_anions += 1;
            }
            md_array_push(q->charged_atom, i, alloc);
            md_array_push(q->charged_group, group[r], alloc);
            md_array_push(q->charged_sign, (int8_t)(net[r] > 0 ? 1 : -1), alloc);
        }
        q->num_charged = md_array_size(q->charged_atom);
    }

    if (q->params.kinds & MD_INTERACTION_HALOGEN) {
        // A Cl, Br or I of one bond, to a C
        for (uint32_t i = 0; i < N; ++i) {
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
            if (z != MD_Z_Cl && z != MD_Z_Br && z != MD_Z_I) continue;
            if (md_bond_conn_count(bd, i) != 1 || !bd->conn.offset) continue;
            md_bond_iter_t it = md_bond_iter(bd, i);
            const uint32_t c = (uint32_t)md_bond_iter_atom_index(&it);
            if (c >= N || md_atom_atomic_number(&sys->atom, c) != MD_Z_C) continue;
            md_array_push(q->halogen, i, alloc);
            md_array_push(q->halogen_c, c, alloc);
            md_array_push(q->halogen_vdw, md_atomic_number_vdw_radius(z), alloc);
        }
        q->num_halogens = md_array_size(q->halogen);
        if (q->num_halogens) {
            uint8_t* role = md_alloc(temp, N);
            md_memset(role, 0, N);
            md_hbond_perceive_roles(role, NULL, sys, NULL, MD_HBOND_ROLES_SULFUR);
            for (uint32_t i = 0; i < N; ++i) {
                if (!(role[i] & MD_HBOND_ROLE_ACCEPTOR)) continue;
                md_array_push(q->acceptor, i, alloc);
                md_array_push(q->acceptor_vdw, md_atomic_number_vdw_radius(md_atom_atomic_number(&sys->atom, i)), alloc);
            }
            q->num_acceptors = md_array_size(q->acceptor);
        }
    }

    // The exclusions, of the charged atoms then of the halogens
    if (q->params.exclude_bonds > 0 && bd->conn.offset && (q->num_charged || q->num_halogens)) {
        uint8_t* depth = md_alloc(temp, N);
        md_memset(depth, 0, N);
        const uint32_t max_bonds = MD_MIN(q->params.exclude_bonds, 254);
        md_array_push(q->excl_off, 0, alloc);
        for (size_t k = 0; k < q->num_charged; ++k) push_exclusions(q, sys, q->charged_atom[k], max_bonds, depth, temp);
        for (size_t k = 0; k < q->num_halogens; ++k) push_exclusions(q, sys, q->halogen[k], max_bonds, depth, temp);
    }
    md_vm_arena_destroy(temp);
    return true;
}

bool md_interaction_query_can_find(const md_interaction_query_t* q, md_interaction_kind_t kind) {
    if (!q || !(q->params.kinds & kind)) return false;
    if (kind == MD_INTERACTION_IONIC) return q->num_cations && q->num_anions;
    if (kind == MD_INTERACTION_HALOGEN) return q->num_halogens && q->num_acceptors;
    return false;
}

void md_interaction_set_free(md_interaction_set_t* set) {
    if (!set || !set->alloc) return;
    md_array_free(set->a, set->alloc);
    md_array_free(set->b, set->alloc);
    md_array_free(set->kind, set->alloc);
    md_array_free(set->dist, set->alloc);
    md_memset(set, 0, sizeof(*set));
}

typedef struct cand_t {
    uint32_t kind, ga, gb;      // the kind, the groups (the halogen and acceptor of a halogen bond)
    uint32_t a, b;
    float    d;
} cand_t;

typedef struct eval_t {
    const md_interaction_query_t* q;
    const md_float3_t* xyz;
    const md_unitcell_t* cell;
    const uint32_t* ext;        // the query's indices of the external points
    const uint32_t* in;         // and of the internal ones
    md_array(cand_t) cands;
    md_allocator_i* alloc;
} eval_t;

static md_float3_t diff(const eval_t* e, uint32_t from, uint32_t to) {
    md_float3_t d = md_float3_sub(e->xyz[to], e->xyz[from]);
    md_geom_min_image(d.elem, 0, 1, e->cell);
    return d;
}

static bool in_sel(const md_interaction_query_t* q, uint32_t a, uint32_t b) {
    return !q->sel || (q->sel[a] && q->sel[b]);
}

static void ionic_pairs(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t n, void* user) {
    eval_t* e = user;
    const md_interaction_query_t* q = e->q;
    for (size_t k = 0; k < n; ++k) {
        const uint32_t ci = e->ext[i_idx[k]], ai = e->in[j_idx[k]];
        const uint32_t a = q->charged_atom[ci], b = q->charged_atom[ai];
        if (!in_sel(q, a, b) || is_excluded(q, ci, b)) continue;
        const cand_t c = {MD_INTERACTION_IONIC, q->charged_group[ci], q->charged_group[ai], a, b, sqrtf(d2[k])};
        md_array_push(e->cands, c, e->alloc);
    }
}

static void halogen_pairs(const uint32_t* i_idx, const uint32_t* j_idx, const float* d2, size_t n, void* user) {
    eval_t* e = user;
    const md_interaction_query_t* q = e->q;
    const float cos_min = cosf((float)MD_DEG_TO_RAD(q->params.halogen_min_cxa));
    for (size_t k = 0; k < n; ++k) {
        const uint32_t hi = e->ext[i_idx[k]], ai = e->in[j_idx[k]];
        const uint32_t x = q->halogen[hi], acc = q->acceptor[ai], c = q->halogen_c[hi];
        if (acc == c || !in_sel(q, x, acc) || is_excluded(q, q->num_charged + hi, acc)) continue;
        const float d = sqrtf(d2[k]);
        if (d > q->halogen_vdw[hi] + q->acceptor_vdw[ai] + q->params.halogen_tolerance) continue;
        const md_float3_t xc = diff(e, x, c), xa = diff(e, x, acc);
        const float lc = md_float3_length(xc), la = md_float3_length(xa);
        if (lc <= 0.0f || la <= 0.0f) continue;
        // The angle C-X...A at X, between X->C and X->A: its cosine -1 when linear, and above the gate's when under it
        const float cos_cxa = md_float3_dot(xc, xa) / (lc * la);
        if (cos_cxa > cos_min) continue;
        const cand_t cd = {MD_INTERACTION_HALOGEN, hi, ai, x, acc, d};
        md_array_push(e->cands, cd, e->alloc);
    }
}

static int compare_cand(const void* pa, const void* pb) {
    const cand_t* a = pa;
    const cand_t* b = pb;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->ga != b->ga) return a->ga < b->ga ? -1 : 1;
    if (a->gb != b->gb) return a->gb < b->gb ? -1 : 1;
    if (a->d != b->d) return a->d < b->d ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : a->b > b->b;
}

static int compare_bond(const void* pa, const void* pb) {
    const cand_t* a = pa;
    const cand_t* b = pb;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return a->b < b->b ? -1 : a->b > b->b;
}

// The positions of the entries of a query's atoms that take part, gathered for the search, and their indices
static size_t gather(md_float3_t** out_xyz, uint32_t** out_idx, const md_float3_t* xyz, const uint32_t* atoms, size_t n,
                     const int8_t* sign, int8_t want, md_allocator_i* alloc) {
    md_float3_t* p = md_alloc(alloc, (n ? n : 1) * sizeof(md_float3_t));
    uint32_t* idx = md_alloc(alloc, (n ? n : 1) * sizeof(uint32_t));
    size_t m = 0;
    for (size_t k = 0; k < n; ++k) {
        if (sign && sign[k] != want) continue;
        p[m] = xyz[atoms[k]];
        idx[m] = (uint32_t)k;
        m += 1;
    }
    *out_xyz = p;
    *out_idx = idx;
    return m;
}

bool md_interaction_query_eval(md_interaction_set_t* out, const md_interaction_query_t* q, const md_system_state_t* state, md_allocator_i* alloc) {
    MD_ASSERT(out);
    md_memset(out, 0, sizeof(*out));
    if (!q || !q->alloc || !state || !alloc) {
        MD_LOG_ERROR("Interactions: missing query, state or allocator");
        return false;
    }
    if (state->num_atoms != q->num_atoms || (q->num_atoms && !state->xyz)) {
        MD_LOG_ERROR("Interactions: the state does not match the system the query was prepared for");
        return false;
    }
    out->alloc = alloc;
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    eval_t e = {.q = q, .xyz = state->xyz, .cell = &state->unitcell, .alloc = arena};

    if (md_interaction_query_can_find(q, MD_INTERACTION_IONIC)) {
        md_float3_t *cat_xyz, *an_xyz;
        uint32_t *cat_idx, *an_idx;
        const size_t nc = gather(&cat_xyz, &cat_idx, state->xyz, q->charged_atom, q->num_charged, q->charged_sign, 1, arena);
        const size_t na = gather(&an_xyz, &an_idx, state->xyz, q->charged_atom, q->num_charged, q->charged_sign, -1, arena);
        if (nc && na) {
            const double r = q->params.ionic_max_dist;
            md_spatial_acc_t sa = {.alloc = arena};
            md_spatial_acc_init(&sa, &(md_spatial_acc_desc_t){.xyz = (const float*)an_xyz, .stride = sizeof(md_float3_t), .count = na, .cutoff = r, .unitcell = &state->unitcell});
            e.ext = cat_idx;
            e.in = an_idx;
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&sa, (const float*)cat_xyz, sizeof(md_float3_t), NULL, nc, r, ionic_pairs, &e, 0);
        }
    }
    if (md_interaction_query_can_find(q, MD_INTERACTION_HALOGEN)) {
        md_float3_t *x_xyz, *a_xyz;
        uint32_t *x_idx, *a_idx;
        const size_t nx = gather(&x_xyz, &x_idx, state->xyz, q->halogen, q->num_halogens, NULL, 0, arena);
        const size_t na = gather(&a_xyz, &a_idx, state->xyz, q->acceptor, q->num_acceptors, NULL, 0, arena);
        float r = 0.0f;
        for (size_t k = 0; k < nx; ++k) r = MD_MAX(r, q->halogen_vdw[k]);
        float ra = 0.0f;
        for (size_t k = 0; k < na; ++k) ra = MD_MAX(ra, q->acceptor_vdw[k]);
        r += ra + MD_MAX(0.0f, q->params.halogen_tolerance);
        if (nx && na && r > 0.0f) {
            md_spatial_acc_t sa = {.alloc = arena};
            md_spatial_acc_init(&sa, &(md_spatial_acc_desc_t){.xyz = (const float*)a_xyz, .stride = sizeof(md_float3_t), .count = na, .cutoff = r, .unitcell = &state->unitcell});
            e.ext = x_idx;
            e.in = a_idx;
            md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&sa, (const float*)x_xyz, sizeof(md_float3_t), NULL, nx, r, halogen_pairs, &e, 0);
        }
    }

    // An ionic bond per pair of groups, of their closest atoms; every halogen bond. Then in the order of their keys
    const size_t nc = md_array_size(e.cands);
    if (nc > 1) qsort(e.cands, nc, sizeof(cand_t), compare_cand);
    size_t m = 0;
    for (size_t k = 0; k < nc; ++k) {
        const cand_t* c = &e.cands[k];
        if (m && c->kind == MD_INTERACTION_IONIC && e.cands[m - 1].kind == c->kind && e.cands[m - 1].ga == c->ga && e.cands[m - 1].gb == c->gb) continue;
        e.cands[m++] = *c;
    }
    if (m > 1) qsort(e.cands, m, sizeof(cand_t), compare_bond);
    for (size_t k = 0; k < m; ++k) {
        md_array_push(out->a, e.cands[k].a, alloc);
        md_array_push(out->b, e.cands[k].b, alloc);
        md_array_push(out->kind, (uint8_t)e.cands[k].kind, alloc);
        md_array_push(out->dist, e.cands[k].d, alloc);
    }
    out->count = m;
    md_vm_arena_destroy(arena);
    return true;
}

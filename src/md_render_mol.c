#include "md_render_mol.h"

#include <md_system.h>
#include <md_bond.h>
#include <md_system_state.h>
#include <md_unitcell.h>
#include <core/md_log.h>
#include <core/md_simd.h>

#include <float.h>
#include <string.h>

#include "md_render_mol_shaders.inl"

static uint32_t bond_class(md_bond_flags_t flags) {
    switch (md_bond_kind(flags)) {
    case MD_BOND_KIND_COORDINATE:  return MD_RENDER_MOL_BOND_COORDINATE;
    case MD_BOND_KIND_ELASTIC:     return MD_RENDER_MOL_BOND_ELASTIC;
    case MD_BOND_KIND_METALLIC:    return MD_RENDER_MOL_BOND_METALLIC;
    case MD_BOND_KIND_INTERACTION: return MD_RENDER_MOL_BOND_INTERACTION;
    case MD_BOND_KIND_HYDROGEN:    return MD_RENDER_MOL_BOND_HYDROGEN;
    default:                       return MD_RENDER_MOL_BOND_COVALENT;
    }
}

static uint64_t align256(uint64_t x) { return (x + 255u) & ~(uint64_t)255u; }

#define CLUSTER     64u     /* atoms or bonds per static cluster, at most       */
#define CELL_ATOMS  48u     /* grid capacity: one cell per this many small atoms (the
                               kernels aim at 64 per cell, so about 3/4 are used) */

/* ---- Layouts of the blocks ------------------------------------------------------- */

/* The topology block, as offsets */
typedef struct topo_layout_t {
    uint64_t state, atom_type, type_radius, atom_order, atom_point, resident, atom_cluster_offset, bond_pairs, bond_ids;
    uint64_t bond_cluster_offset[MD_RENDER_MOL_BOND_CLASS_COUNT];
    uint64_t small_tables, size;
} topo_layout_t;

/* The small tables, as offsets from small_tables */
typedef struct small_layout_t {
    uint64_t mol_atom_offset, mol_atoms, mol_bond_offset, mol_bond_local, mol_bond_ids, mixed_pos, mixed_atoms, size;
} small_layout_t;

/* Per small molecule: its atoms, and its covalent bonds as the two atoms'
   places in the molecule (k0 | k1 << 16). Then the mixed bonds: those that
   touch a small molecule without being within one (their place in the bond
   arrays, and their atoms), rewritten by every regroup. */
static small_layout_t small_layout(uint32_t n_mol, uint32_t n_atoms, uint32_t n_bonds, uint32_t n_mixed) {
    small_layout_t L;
    uint64_t o = 0;
    L.mol_atom_offset = o; o = align256(o + (uint64_t)(n_mol + 1) * 4u);
    L.mol_atoms       = o; o = align256(o + (uint64_t)n_atoms * 4u);
    L.mol_bond_offset = o; o = align256(o + (uint64_t)(n_mol + 1) * 4u);
    L.mol_bond_local  = o; o = align256(o + (uint64_t)n_bonds * 4u);
    L.mol_bond_ids    = o; o = align256(o + (uint64_t)n_bonds * 4u);
    L.mixed_pos       = o; o = align256(o + (uint64_t)n_mixed * 4u);
    L.mixed_atoms     = o; o = align256(o + (uint64_t)n_mixed * 8u);
    L.size = o;
    return L;
}

/* The per-frame block (scratch of the update) */
#define SCAN_BLOCK 1024u
typedef struct dyn_layout_t {
    uint64_t points, velocity, marker, params, mol_slot, cnt, sums, size, size_velocity;
    uint32_t scan_n, scan_blocks;
} dyn_layout_t;

/* The points (float4: position, radius) in cluster order: what is drawn.
   Then the scratch of the regroup. Last, only in a block allocated with
   them, the points' velocities. */
static dyn_layout_t dyn_layout(uint32_t n_atoms, uint32_t n_mol, uint32_t cells) {
    dyn_layout_t L;
    L.scan_n      = cells + 1;
    L.scan_blocks = (L.scan_n + SCAN_BLOCK - 1) / SCAN_BLOCK;
    uint64_t o = 0;
    L.points   = o; o = align256(o + (uint64_t)n_atoms * 16u);
    L.marker   = o; o = align256(o + 32u);
    L.params   = o; o = align256(o + 16u * 4u);
    L.mol_slot = o; o = align256(o + (uint64_t)n_mol * 16u);
    L.cnt      = o; o = align256(o + (uint64_t)L.scan_n * 16u);
    L.sums     = o; o = align256(o + (uint64_t)L.scan_blocks * 16u);
    L.size = o;
    L.velocity = o; o = align256(o + (uint64_t)n_atoms * 8u);
    L.size_velocity = o;
    return L;
}

/* The frames block: each slot's parameters (32 words), the slots' frames
   (x, y, z per point in 16 bits each, three words to two points), then the
   old points of the small atoms (a regroup's) */
typedef struct frames_layout_t {
    uint64_t params, frame, buffers, src_point, size;
} frames_layout_t;

static uint64_t frame_bytes(uint32_t n_atoms) { return ((3ull * n_atoms + 1u) / 2u) * 4u; }

static frames_layout_t frames_layout(uint32_t n_atoms, uint32_t n_small, uint32_t slots) {
    frames_layout_t L;
    uint64_t o = 0;
    L.params    = o; o = align256(o + (uint64_t)slots * 128u);
    L.frame     = align256(frame_bytes(n_atoms));
    L.buffers   = o; o += (uint64_t)slots * L.frame;
    L.src_point = o; o = align256(o + (uint64_t)n_small * 4u);
    L.size = o;
    return L;
}

/* ---- Plan ------------------------------------------------------------------------- */

/* Greedy packing of units (runs of atoms) into clusters of at most CLUSTER;
   a unit larger than that is split into chunks */
typedef struct packer_t {
    uint32_t* order;       /* atoms out                 */
    uint32_t* offset;      /* cluster offsets out       */
    uint32_t  count;       /* clusters closed           */
    uint32_t  fill;        /* atoms in the open cluster */
    uint32_t  pos;         /* atoms written             */
} packer_t;

static void pack_close(packer_t* p) {
    if (p->fill) { p->count++; p->offset[p->count] = p->pos; p->fill = 0; }
}

static void pack_unit(packer_t* p, const uint32_t* atoms, uint32_t n) {
    if (n > CLUSTER) {
        pack_close(p);
        for (uint32_t i = 0; i < n; i += CLUSTER) {
            uint32_t m = n - i < CLUSTER ? n - i : CLUSTER;
            memcpy(p->order + p->pos, atoms + i, m * 4u);
            p->pos += m; p->fill = m;
            pack_close(p);
        }
        return;
    }
    if (p->fill + n > CLUSTER) pack_close(p);
    memcpy(p->order + p->pos, atoms, n * 4u);
    p->pos += n; p->fill += n;
}

bool md_render_mol_plan(md_render_mol_plan_t* out, const md_system_t* sys, const md_bitset_t* resident, md_temp_scope_t temp) {
    if (!out || !sys) return false;
    memset(out, 0, sizeof(*out));
    md_render_mol_topo_t* T = &out->layout;
    const uint32_t n      = (uint32_t)sys->atom.count;
    const uint32_t nb_sys = (uint32_t)(sys->bond.count + sys->bond.noncovalent_count);
    T->topology_version = sys->topology_version;
    T->atom_count       = n;

    /* The resident atoms: all, or the members of `resident` below n */
    uint8_t* res = NULL;
    uint32_t n_res = n;
    if (resident) {
        res = md_temp_alloc_zero(temp, (size_t)n + 1);
        if (!res) return false;
        n_res = 0;
        md_bitset_iter_t it = md_bitset_iter_range(*resident, 0, n);
        while (md_bitset_iter_next(&it)) { res[it.idx] = 1; n_res++; }
        if (n_res == n) res = NULL;          /* all of them: no rank needed */
    }
    #define RES(a) (!res || res[a])
    T->point_count = n_res;

    /* The bonds kept: both ends resident (and atoms at all) */
    uint32_t* kept = md_temp_alloc(temp, ((size_t)nb_sys + 1) * 4u);
    if (!kept) return false;
    uint32_t nb = 0;
    for (uint32_t i = 0; i < nb_sys; ++i) {
        uint32_t a = (uint32_t)sys->bond.pairs[i].idx[0], b = (uint32_t)sys->bond.pairs[i].idx[1];
        if (a < n && b < n && RES(a) && RES(b)) kept[nb++] = i;
    }
    T->bond_count = nb;
    #define BA(j) ((uint32_t)sys->bond.pairs[kept[j]].idx[0])
    #define BB(j) ((uint32_t)sys->bond.pairs[kept[j]].idx[1])

    /* Molecules, and which are small */
    const md_molecule_data_t* M = &sys->molecule;
    const uint32_t n_mol = (uint32_t)M->count;
    int32_t*  mol_of = md_temp_alloc(temp, ((size_t)n + 1) * sizeof(int32_t));
    uint8_t*  small  = md_temp_alloc_zero(temp, (size_t)n_mol + 1);
    if (!mol_of || !small) return false;
    for (uint32_t i = 0; i < n; ++i) mol_of[i] = -1;
    /* small[m]: 1 for a small molecule (by its size in the system), 2 for one
       with resident atoms (those are in the tables) */
    uint32_t small_mols = 0, small_atoms = 0;
    for (uint32_t m = 0; m < n_mol; ++m) {
        uint32_t b = M->offset[m], e = M->offset[m + 1], r = 0;
        for (uint32_t k = b; k < e; ++k) {
            int32_t a = M->atom_idx[k];
            if (a >= 0 && (uint32_t)a < n) { mol_of[a] = (int32_t)m; r += RES(a); }
        }
        if (e > b && e - b <= MD_RENDER_MOL_SMALL) {
            small[m] = r ? 2 : 1;
            if (r) { small_mols++; small_atoms += r; }
        }
    }
    #define IS_SMALL(a) (mol_of[a] >= 0 && small[mol_of[a]])

    /* Static atoms: whole residues, in file order, packed up to CLUSTER, never
       across molecules; then whatever no residue holds, in file order */
    uint32_t* order  = md_temp_alloc(temp, ((size_t)n + 1) * 4u);
    uint32_t* coff   = md_temp_alloc_zero(temp, ((size_t)n + 2) * 4u);
    uint8_t*  placed = md_temp_alloc_zero(temp, (size_t)n + 1);
    uint32_t* tmp    = md_temp_alloc(temp, ((size_t)n + 1) * 4u);
    uint32_t* cl_of  = md_temp_alloc(temp, ((size_t)n + 1) * 4u);
    if (!order || !coff || !placed || !tmp || !cl_of) return false;
    packer_t P = { order, coff, 0, 0, 0 };
    const md_residue_data_t* R = &sys->residue;
    int32_t prev_mol = -2;
    for (size_t r = 0; r < R->count; ++r) {
        uint32_t b = R->atom_offset[r], e = R->atom_offset[r + 1], k = 0;
        if (e > n) e = n;
        for (uint32_t i = b; i < e; ++i) if (!placed[i] && !IS_SMALL(i) && RES(i)) { placed[i] = 1; tmp[k++] = i; }
        if (!k) continue;
        int32_t mol = mol_of[tmp[0]];
        if (mol != prev_mol) pack_close(&P);
        prev_mol = mol;
        pack_unit(&P, tmp, k);
    }
    pack_close(&P);
    for (uint32_t i = 0; i < n; ++i) if (!placed[i] && !IS_SMALL(i) && RES(i)) { placed[i] = 1; pack_unit(&P, &i, 1); }
    pack_close(&P);
    const uint32_t static_atoms = P.pos, static_clusters = P.count;
    for (uint32_t c = 0; c < static_clusters; ++c)
        for (uint32_t k = coff[c]; k < coff[c + 1]; ++k) cl_of[order[k]] = c;

    /* Small molecules: their atoms by molecule (slot order, the first is the
       one the grid places), and their covalent bonds */
    uint32_t* sm_index = md_temp_alloc(temp, ((size_t)n_mol + 1) * 4u);   /* molecule -> small index */
    uint32_t* sm_aoff  = md_temp_alloc(temp, ((size_t)small_mols + 1) * 4u);
    uint32_t* sm_atoms = md_temp_alloc(temp, ((size_t)small_atoms + 1) * 4u);
    uint32_t* sm_boff  = md_temp_alloc_zero(temp, ((size_t)small_mols + 2) * 4u);
    uint8_t*  local    = md_temp_alloc_zero(temp, (size_t)n + 1);     /* a small atom's place in its molecule */
    if (!sm_index || !sm_aoff || !sm_atoms || !sm_boff || !local) return false;
    {
        uint32_t k = 0, a = 0;
        for (uint32_t m = 0; m < n_mol; ++m) {
            if (small[m] != 2) continue;
            sm_index[m] = k;
            sm_aoff[k++] = a;
            uint32_t l = 0;                  /* its place among the molecule's resident atoms */
            for (uint32_t s = M->offset[m]; s < M->offset[m + 1]; ++s) {
                uint32_t at = (uint32_t)M->atom_idx[s];
                if (!RES(at)) continue;
                local[at] = (uint8_t)l++;
                sm_atoms[a++] = at;
            }
        }
        sm_aoff[k] = a;
    }

    /* Bonds by class; within a class the static ones by the cluster of their
       first static atom (stable), then those within small molecules
       (covalent only: other kinds between small molecules are rare and stay
       static, in chunks) */
    uint32_t* cls  = md_temp_alloc(temp, ((size_t)nb + 1) * 4u);
    uint32_t* key  = md_temp_alloc(temp, ((size_t)nb + 1) * 4u);
    if (!cls || !key) return false;
    const uint32_t NONE = static_clusters;            /* no static atom: sorted last */
    const uint32_t nkeys = static_clusters + 2;       /* + NONE, + small            */
    const uint32_t SMALL_KEY = static_clusters + 1;
    uint32_t small_bonds = 0;
    for (uint32_t i = 0; i < nb; ++i) {
        cls[i] = bond_class(sys->bond.flags ? sys->bond.flags[kept[i]] : 0);
        uint32_t a = BA(i), b = BB(i);
        if (cls[i] == MD_RENDER_MOL_BOND_COVALENT && IS_SMALL(a) && mol_of[a] == mol_of[b]) {
            key[i] = SMALL_KEY;
            sm_boff[sm_index[mol_of[a]] + 1]++;
            small_bonds++;
        } else {
            key[i] = !IS_SMALL(a) ? cl_of[a] : !IS_SMALL(b) ? cl_of[b] : NONE;
        }
    }
    for (uint32_t k = 0; k < small_mols; ++k) sm_boff[k + 1] += sm_boff[k];

    /* Counting sort by (class, key) */
    const size_t nbuckets = (size_t)MD_RENDER_MOL_BOND_CLASS_COUNT * nkeys;
    uint32_t* bucket = md_temp_alloc_zero(temp, (nbuckets + 1) * 4u);
    uint32_t* sorted = md_temp_alloc(temp, ((size_t)nb + 1) * 4u);
    uint32_t* sm_bfill = md_temp_alloc_zero(temp, ((size_t)small_mols + 1) * 4u);
    uint32_t* sm_blocal = md_temp_alloc(temp, ((size_t)small_bonds + 1) * 4u);
    uint32_t* sm_bids   = md_temp_alloc(temp, ((size_t)small_bonds + 1) * 4u);
    if (!bucket || !sorted || !sm_bfill || !sm_blocal || !sm_bids) return false;
    for (uint32_t i = 0; i < nb; ++i) bucket[(size_t)cls[i] * nkeys + key[i] + 1]++;
    for (size_t k = 0; k < nbuckets; ++k) bucket[k + 1] += bucket[k];
    for (uint32_t i = 0; i < nb; ++i) sorted[bucket[(size_t)cls[i] * nkeys + key[i]]++] = i;
    for (uint32_t i = 0; i < nb; ++i) {
        if (key[i] != SMALL_KEY) continue;
        uint32_t k = sm_index[mol_of[BA(i)]];
        uint32_t j = sm_boff[k] + sm_bfill[k]++;
        sm_blocal[j] = (uint32_t)local[BA(i)] | ((uint32_t)local[BB(i)] << 16);
        sm_bids[j] = kept[i];
    }

    /* The grid's capacity, and the clusters of every class */
    const uint32_t cells = small_atoms ? (small_atoms + CELL_ATOMS - 1) / CELL_ATOMS : 0;
    const uint32_t bond_cells = small_bonds ? cells : 0;
    uint32_t class_count[MD_RENDER_MOL_BOND_CLASS_COUNT] = {0};
    for (uint32_t i = 0; i < nb; ++i) class_count[cls[i]]++;
    {
        uint32_t off = 0;
        for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) { T->bond_class_offset[c] = off; off += class_count[c]; }
        T->bond_class_offset[MD_RENDER_MOL_BOND_CLASS_COUNT] = off;
    }
    /* Static bond clusters: a new one where the key changes or at CLUSTER */
    uint32_t* boff = md_temp_alloc(temp, ((size_t)nb + MD_RENDER_MOL_BOND_CLASS_COUNT * ((size_t)cells + 2) + 1) * 4u);
    if (!boff) return false;
    uint32_t boff_base[MD_RENDER_MOL_BOND_CLASS_COUNT];
    uint32_t nboff = 0;
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
        const uint32_t beg = T->bond_class_offset[c];
        const uint32_t n_static = class_count[c] - (c == MD_RENDER_MOL_BOND_COVALENT ? small_bonds : 0);
        boff_base[c] = nboff;
        uint32_t cnt = 0, run = 0, prev = UINT32_MAX;
        for (uint32_t j = 0; j < n_static; ++j) {
            uint32_t k = key[sorted[beg + j]];
            if (j == 0 || k != prev || run == CLUSTER) { boff[nboff + cnt++] = j; run = 0; }
            prev = k; run++;
        }
        if (c == MD_RENDER_MOL_BOND_COVALENT) T->bond_cluster_count[c] = cnt + bond_cells;
        else                                  T->bond_cluster_count[c] = cnt;
        nboff += T->bond_cluster_count[c] + 1;
        /* closing entry; the covalent class's cell entries follow below */
        if (c != MD_RENDER_MOL_BOND_COVALENT || !bond_cells) boff[nboff - 1] = n_static;
    }

    /* The small molecules' initial grouping (before any update): in molecule
       order, cell k holding molecules [k * mols / cells, (k + 1) * mols / cells) */
    const uint32_t c0 = boff_base[MD_RENDER_MOL_BOND_COVALENT] + (T->bond_cluster_count[0] - bond_cells);
    const uint32_t static_bonds0 = class_count[0] - small_bonds;
    {
        uint32_t m = 0, a = 0, b = 0;
        for (uint32_t k = 0; k < cells; ++k) {
            coff[static_clusters + k] = static_atoms + a;
            if (bond_cells) boff[c0 + k] = static_bonds0 + b;
            uint32_t end = (uint32_t)(((uint64_t)(k + 1) * small_mols) / cells);
            for (; m < end; ++m) {
                for (uint32_t s = sm_aoff[m]; s < sm_aoff[m + 1]; ++s) order[static_atoms + a++] = sm_atoms[s];
                b += sm_boff[m + 1] - sm_boff[m];
            }
        }
        if (cells) coff[static_clusters + cells] = static_atoms + a;
        if (bond_cells) boff[c0 + cells] = static_bonds0 + b;
    }
    /* Every atom's point (the inverse of the order), and the mixed bonds */
    uint32_t* inv = md_temp_alloc(temp, ((size_t)n + 1) * 4u);
    if (!inv) return false;
    for (uint32_t k = 0; k < n; ++k) inv[k] = UINT32_MAX;          /* not resident: no point */
    for (uint32_t k = 0; k < n_res; ++k) inv[order[k]] = k;
    uint32_t n_mixed = 0;
    for (uint32_t i = 0; i < nb; ++i) {
        if (key[i] != SMALL_KEY && (IS_SMALL(BA(i)) || IS_SMALL(BB(i)))) n_mixed++;
    }
    T->atom_cluster_count   = static_clusters + cells;
    T->mixed_bond_count     = n_mixed;
    T->small_molecule_count = small_mols;
    T->small_atom_count     = small_atoms;
    T->small_bond_count     = small_bonds;
    T->grid_cells           = cells;
    #undef IS_SMALL

    /* Layout of the block */
    topo_layout_t L;
    small_layout_t SL = small_layout(small_mols, small_atoms, small_bonds, n_mixed);
    uint64_t o = 0;
    L.state               = o; o = align256(o + 4u);            /* the regroups' count: word 0 */
    L.atom_type           = o; o = align256(o + ((uint64_t)n + 1) / 2 * 4u);
    L.type_radius         = o; o = align256(o + ((uint64_t)sys->atom.type.count + 1) * 4u);
    L.atom_order          = o; o = align256(o + (uint64_t)n_res * 4u);
    L.atom_point          = o; o = align256(o + (uint64_t)n * 4u);
    const uint32_t rwords = (n + 31) / 32;
    L.resident            = res ? o : 0; o = res ? align256(o + (uint64_t)rwords * 8u) : o;
    L.atom_cluster_offset = o; o = align256(o + (uint64_t)(T->atom_cluster_count + 1) * 4u);
    L.bond_pairs          = o; o = align256(o + (uint64_t)nb * 8u);
    L.bond_ids            = o; o = align256(o + (uint64_t)nb * 4u);
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) L.bond_cluster_offset[c] = o + (uint64_t)boff_base[c] * 4u;
    o = align256(o + (uint64_t)nboff * 4u);
    L.small_tables = o; o += small_mols ? SL.size : 0;
    L.size = o ? o : 256;

    uint8_t* img = md_temp_alloc_zero(temp, L.size);
    if (!img) return false;
    {
        const md_particle_type_data_t* PT = &sys->atom.type;
        uint16_t* type   = (uint16_t*)(img + L.atom_type);
        float*    radius = (float*)(img + L.type_radius);
        for (uint32_t i = 0; i < n; ++i) type[i] = sys->atom.type_idx ? sys->atom.type_idx[i] : 0;
        for (size_t t = 0; t < PT->count; ++t) radius[t] = PT->radius ? PT->radius[t] : 0.0f;
        T->type_count = (uint32_t)PT->count;
    }
    memcpy(img + L.atom_order, order, (size_t)n_res * 4u);
    memcpy(img + L.atom_point, inv, (size_t)n * 4u);
    if (res) {
        /* the resident atoms as bits, then each word's count of those before it */
        uint32_t* bits = (uint32_t*)(img + L.resident);
        uint32_t* rank = bits + rwords;
        for (uint32_t i = 0; i < n; ++i) if (res[i]) bits[i >> 5] |= 1u << (i & 31);
        uint32_t sum = 0;
        for (uint32_t w = 0; w < rwords; ++w) { rank[w] = sum; sum += (uint32_t)md_popcnt32(bits[w]); }
    }
    memcpy(img + L.atom_cluster_offset, coff, (size_t)(T->atom_cluster_count + 1) * 4u);
    uint32_t* mixed_pos   = small_mols ? (uint32_t*)(img + L.small_tables + SL.mixed_pos)   : NULL;
    uint32_t* mixed_atoms = small_mols ? (uint32_t*)(img + L.small_tables + SL.mixed_atoms) : NULL;
    {
        /* Bonds between points: an atom's point is inv[atom] */
        uint32_t* pairs = (uint32_t*)(img + L.bond_pairs);
        uint32_t* ids   = (uint32_t*)(img + L.bond_ids);
        uint32_t  nm    = 0;
        for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
            const uint32_t beg = T->bond_class_offset[c];
            const uint32_t n_static = class_count[c] - (c == MD_RENDER_MOL_BOND_COVALENT ? small_bonds : 0);
            for (uint32_t j = 0; j < n_static; ++j) {
                uint32_t i = sorted[beg + j];
                uint32_t a = BA(i), b = BB(i);
                pairs[2 * (beg + j) + 0] = inv[a];
                pairs[2 * (beg + j) + 1] = inv[b];
                ids[beg + j] = kept[i];
                if ((mol_of[a] >= 0 && small[mol_of[a]]) || (mol_of[b] >= 0 && small[mol_of[b]])) {
                    mixed_pos[nm] = beg + j;
                    mixed_atoms[2 * nm + 0] = a;
                    mixed_atoms[2 * nm + 1] = b;
                    nm++;
                }
            }
        }
        /* the small molecules' bonds, in the initial grouping: molecule order */
        for (uint32_t k = 0; k < small_bonds; ++k) {
            uint32_t i = sm_bids[k];         /* an index in sys->bond */
            pairs[2 * (static_bonds0 + k) + 0] = inv[sys->bond.pairs[i].idx[0]];
            pairs[2 * (static_bonds0 + k) + 1] = inv[sys->bond.pairs[i].idx[1]];
            ids[static_bonds0 + k] = i;
        }
    }
    memcpy(img + L.bond_cluster_offset[0], boff, (size_t)nboff * 4u);
    if (small_mols) {
        uint8_t* st = img + L.small_tables;
        memcpy(st + SL.mol_atom_offset, sm_aoff,  ((size_t)small_mols + 1) * 4u);
        memcpy(st + SL.mol_atoms,       sm_atoms, (size_t)small_atoms * 4u);
        memcpy(st + SL.mol_bond_offset, sm_boff,  ((size_t)small_mols + 1) * 4u);
        memcpy(st + SL.mol_bond_local,  sm_blocal, (size_t)small_bonds * 4u);
        memcpy(st + SL.mol_bond_ids,    sm_bids,   (size_t)small_bonds * 4u);
    }

    T->size                = L.size;
    T->dynamic_size           = dyn_layout(n_res, small_mols, cells).size;
    T->dynamic_size_velocity  = dyn_layout(n_res, small_mols, cells).size_velocity;
    T->resident            = L.resident;
    T->regroups            = L.state;
    T->atom_type           = L.atom_type;
    T->type_radius         = L.type_radius;
    T->atom_order          = L.atom_order;
    T->atom_point          = L.atom_point;
    T->atom_cluster_offset = L.atom_cluster_offset;
    T->bond_pairs          = L.bond_pairs;
    T->bond_ids            = L.bond_ids;
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) T->bond_cluster_offset[c] = L.bond_cluster_offset[c];
    T->small_tables        = small_mols ? L.small_tables : 0;
    #undef RES
    #undef BA
    #undef BB
    out->size = L.size;
    out->data = img;
    return true;
}

bool md_render_mol_write(md_gpu_stream_t s, md_render_mol_topo_t* out, const md_render_mol_plan_t* plan,
                         md_gpu_addr_t base, uint64_t content[MD_RENDER_CONTENT_COUNT]) {
    if (!s || !out || !plan || !base || !plan->data) return false;
    if (!md_gpu_upload(s, base, plan->data, (size_t)plan->size)) return false;
    *out = plan->layout;
    out->base                = base;
    out->atom_type           = base + plan->layout.atom_type;
    out->type_radius         = base + plan->layout.type_radius;
    out->atom_order          = base + plan->layout.atom_order;
    out->atom_point          = base + plan->layout.atom_point;
    out->resident            = plan->layout.resident ? base + plan->layout.resident : 0;
    out->regroups            = base + plan->layout.regroups;
    out->atom_cluster_offset = base + plan->layout.atom_cluster_offset;
    out->bond_pairs          = base + plan->layout.bond_pairs;
    out->bond_ids            = base + plan->layout.bond_ids;
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) out->bond_cluster_offset[c] = base + plan->layout.bond_cluster_offset[c];
    out->small_tables        = plan->layout.small_molecule_count ? base + plan->layout.small_tables : 0;
    if (content) content[MD_RENDER_CONTENT_GEOMETRY]++;
    return true;
}

/* ---- Kernels and the update ----------------------------------------------------- */

struct md_render_mol {
    md_allocator_i* alloc;
    md_gpu_kernel_t k_clear, k_bounds, k_dims, k_count, k_scan, k_write, k_offsets, k_mixed, k_gather, k_mark;
    md_gpu_kernel_t k_fr_gather, k_fr_reorder, k_fr_permute;
};

/* Mirrors Args of md_render_mol.slang */
typedef struct gpu_mol_args_t {
    md_gpu_addr_t position, velocity, radius_override, points, vel_out;
    md_gpu_addr_t params, mol_slot, cnt, sums, types, type_radius;
    md_gpu_addr_t mol_atom_offset, mol_atoms, mol_bond_offset, mol_bond_local, mol_bond_ids, mixed_pos, mixed_atoms;
    md_gpu_addr_t atom_order, atom_point, atom_cluster_offset, bond_pairs, bond_ids, bond_cluster_offset, marker;
    md_gpu_addr_t frame[4], fpar[4], qsrc, qdst, qpar, src_point, resident, regroups;
    uint32_t stride, n_mol, cells, n_atoms, static_atoms, static_atom_clusters, static_bonds, static_bond_clusters;
    uint32_t has_bonds, mode, n, atom_count, type_count, vel_stride, has_velocity, has_radius, mixed_count;
    uint32_t full, ident0, ident1, ident2, _pad;
    uint32_t frame_count, quantized, has_src, has_fvel, has_resident, resident_words;
    float    w[4], dw[4];
} gpu_mol_args_t;
_Static_assert(sizeof(gpu_mol_args_t) == 39 * 8 + 36 * 4, "gpu_mol_args_t mirrors Args");

md_render_mol_t md_render_mol_create(md_gpu_device_t device, md_allocator_i* alloc) {
    if (!device) return NULL;
    if (!alloc) alloc = md_get_heap_allocator();
    md_render_mol_t r = md_alloc(alloc, sizeof(struct md_render_mol));
    if (!r) return NULL;
    memset(r, 0, sizeof(*r));
    r->alloc = alloc;
    md_gpu_kernel_desc_t kd;
    kd = md_shader_md_render_mol_sol_clear_kernel();   r->k_clear   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_bounds_kernel();  r->k_bounds  = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_dims_kernel();    r->k_dims    = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_count_kernel();   r->k_count   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_scan_kernel();    r->k_scan    = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_write_kernel();   r->k_write   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_offsets_kernel(); r->k_offsets = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_mixed_kernel();   r->k_mixed   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_gather_kernel();  r->k_gather  = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_mark_kernel();    r->k_mark    = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_fr_gather_kernel();   r->k_fr_gather  = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_fr_reorder_kernel();  r->k_fr_reorder = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_fr_permute_kernel();  r->k_fr_permute = md_gpu_kernel_create(device, &kd);
    if (!r->k_clear || !r->k_bounds || !r->k_dims || !r->k_count || !r->k_scan || !r->k_write || !r->k_offsets ||
        !r->k_mixed || !r->k_gather || !r->k_mark || !r->k_fr_gather || !r->k_fr_reorder || !r->k_fr_permute) {
        MD_LOG_ERROR("md_render_mol: kernel creation failed: %s", md_gpu_last_error());
        md_render_mol_destroy(r);
        return NULL;
    }
    return r;
}

void md_render_mol_destroy(md_render_mol_t r) {
    if (!r) return;
    md_gpu_kernel_t* k[] = { &r->k_clear, &r->k_bounds, &r->k_dims, &r->k_count, &r->k_scan, &r->k_write, &r->k_offsets,
                             &r->k_mixed, &r->k_gather, &r->k_mark, &r->k_fr_gather, &r->k_fr_reorder, &r->k_fr_permute };
    for (size_t i = 0; i < sizeof(k) / sizeof(k[0]); ++i) if (*k[i]) md_gpu_kernel_destroy(*k[i]);
    md_free(r->alloc, r, sizeof(*r));
}

static md_gpu_grid_t groups_for(uint32_t threads, uint32_t per_group) {
    uint32_t g = (threads + per_group - 1) / per_group;
    if (g > 65535) g = 65535;
    if (g == 0) g = 1;
    return md_gpu_grid(g, 1, 1);
}

/* The update's launches, with a hook after each (the bench's per-pass profile) */
typedef void (*pass_hook_t)(void* user, const char* pass);

static bool launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t g, const gpu_mol_args_t* a, pass_hook_t hook, void* user, const char* name) {
    if (!md_gpu_launch(s, k, g, a, sizeof(*a))) return false;
    if (hook) hook(user, name);
    return true;
}

/* The arguments every kernel shares: the topology, the per-frame block */
static void base_args(gpu_mol_args_t* a, const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic) {
    const uint32_t cells = topo->grid_cells, n_mol = topo->small_molecule_count;
    const dyn_layout_t   D  = dyn_layout(topo->point_count, n_mol, cells);
    const small_layout_t SL = small_layout(n_mol, topo->small_atom_count, topo->small_bond_count, topo->mixed_bond_count);
    memset(a, 0, sizeof(*a));
    if (dynamic) {
        a->points         = dynamic + D.points;
        a->vel_out        = dynamic + D.velocity;
        a->params         = dynamic + D.params;
        a->mol_slot       = dynamic + D.mol_slot;
        a->cnt            = dynamic + D.cnt;
        a->sums           = dynamic + D.sums;
        a->marker         = dynamic + D.marker;
    }
    a->atom_count         = topo->point_count;      /* the kernels go over points */
    a->resident           = topo->resident;
    a->regroups           = topo->regroups;
    a->has_resident       = topo->resident ? 1u : 0u;
    a->resident_words     = (topo->atom_count + 31) / 32;
    a->type_count         = topo->type_count;
    a->types              = topo->atom_type;
    a->type_radius        = topo->type_radius;
    if (topo->small_tables) {
        a->mol_atom_offset = topo->small_tables + SL.mol_atom_offset;
        a->mol_atoms       = topo->small_tables + SL.mol_atoms;
        a->mol_bond_offset = topo->small_tables + SL.mol_bond_offset;
        a->mol_bond_local  = topo->small_tables + SL.mol_bond_local;
        a->mol_bond_ids    = topo->small_tables + SL.mol_bond_ids;
        a->mixed_pos       = topo->small_tables + SL.mixed_pos;
        a->mixed_atoms     = topo->small_tables + SL.mixed_atoms;
    }
    a->atom_order          = topo->atom_order;
    a->atom_point          = topo->atom_point;
    a->atom_cluster_offset = topo->atom_cluster_offset;
    a->bond_pairs          = topo->bond_pairs;      /* the covalent class is first: its bonds start at 0 */
    a->bond_ids            = topo->bond_ids;
    a->bond_cluster_offset = topo->bond_cluster_offset[MD_RENDER_MOL_BOND_COVALENT];
    a->n_mol               = n_mol;
    a->cells               = cells;
    a->n_atoms             = topo->small_atom_count;
    a->static_atoms        = topo->point_count - topo->small_atom_count;
    a->static_atom_clusters= topo->atom_cluster_count - cells;
    a->has_bonds           = topo->small_bond_count ? 1u : 0u;
    a->static_bonds        = (topo->bond_class_offset[MD_RENDER_MOL_BOND_COVALENT + 1] - topo->bond_class_offset[MD_RENDER_MOL_BOND_COVALENT]) - topo->small_bond_count;
    a->static_bond_clusters= a->has_bonds ? topo->bond_cluster_count[MD_RENDER_MOL_BOND_COVALENT] - cells : 0;
    a->mixed_count         = topo->mixed_bond_count;
    a->n                   = D.scan_n;
    a->ident0              = (uint32_t)topo->base;
    a->ident1              = (uint32_t)(topo->base >> 32);
    a->ident2              = (uint32_t)topo->topology_version;
}

/* The small molecules into grid cells, by the positions `a` points at */
static bool run_regroup(md_render_mol_t r, md_gpu_stream_t s, gpu_mol_args_t* a, const md_render_mol_topo_t* topo, pass_hook_t hook, void* user) {
    const uint32_t cells = topo->grid_cells, n_mol = topo->small_molecule_count;
    const dyn_layout_t D = dyn_layout(topo->point_count, n_mol, cells);
    bool ok = true;
    ok = ok && launch(s, r->k_clear,  groups_for(4 * (cells + 1) > 16 ? 4 * (cells + 1) : 16, 64), a, hook, user, "clear");
    ok = ok && launch(s, r->k_bounds, groups_for(n_mol, 64), a, hook, user, "bounds");
    ok = ok && launch(s, r->k_dims,   md_gpu_grid(1, 1, 1), a, hook, user, "dims");
    ok = ok && launch(s, r->k_count,  groups_for(n_mol, 64), a, hook, user, "count");
    a->mode = 0;
    ok = ok && launch(s, r->k_scan, md_gpu_grid(D.scan_blocks, 1, 1), a, hook, user, "scan");
    if (D.scan_blocks > 1) {
        a->mode = 1;
        ok = ok && launch(s, r->k_scan, md_gpu_grid(1, 1, 1), a, hook, user, "scan sums");
        a->mode = 2;
        ok = ok && launch(s, r->k_scan, md_gpu_grid(D.scan_blocks, 1, 1), a, hook, user, "scan add");
    }
    ok = ok && launch(s, r->k_write,   groups_for(n_mol, 64), a, hook, user, "write");
    ok = ok && launch(s, r->k_offsets, groups_for(cells + 1, 64), a, hook, user, "offsets");
    if (topo->mixed_bond_count) ok = ok && launch(s, r->k_mixed, groups_for(topo->mixed_bond_count, 64), a, hook, user, "mixed");
    return ok;
}

static bool frames_match(const md_render_mol_frames_t* fr, const md_render_mol_topo_t* topo) {
    return fr->base && fr->topology_base == topo->base && fr->topology_version == topo->topology_version && fr->point_count == topo->point_count;
}

static md_gpu_addr_t frame_buffer(const md_render_mol_frames_t* fr, uint32_t b) {
    const frames_layout_t L = frames_layout(fr->point_count, 0, fr->slots);
    return fr->base + L.buffers + (uint64_t)b * L.frame;
}

static md_gpu_addr_t frame_params(const md_render_mol_frames_t* fr, uint32_t slot) {
    return fr->base + frames_layout(fr->point_count, 0, fr->slots).params + (uint64_t)slot * 128u;
}

/* Velocities are written (and drawn) when given, or from frames with a
   derivative */
static bool has_velocity(const md_render_mol_atoms_t* atoms) {
    if (atoms->frames) {
        if (atoms->frame_count < 2) return false;
        for (uint32_t k = 0; k < atoms->frame_count && k < 4; ++k) if (atoms->frame_dw[k] != 0.0f) return true;
        return false;
    }
    return atoms->velocity != 0;
}

static bool update_impl(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo,
                        md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                        md_render_mol_update_flags_t flags, uint64_t content[MD_RENDER_CONTENT_COUNT],
                        pass_hook_t hook, void* user) {
    const md_render_mol_frames_t* fr = atoms ? atoms->frames : NULL;
    if (!r || !s || !topo || !atoms || !topo->base || !dynamic || (!fr && !atoms->position)) {
        MD_LOG_ERROR("md_render_mol_update: missing context, stream, topology, per-frame block or positions");
        return false;
    }
    const bool regroup = (flags & MD_RENDER_MOL_UPDATE_REGROUP) && topo->grid_cells;
    const bool radii   = (flags & MD_RENDER_MOL_UPDATE_RADII) != 0;
    if (fr) {
        if (regroup) {
            MD_LOG_ERROR("md_render_mol_update: with frames, the regroup belongs to md_render_mol_frame_upload");
            return false;
        }
        if (!frames_match(fr, topo)) {
            MD_LOG_ERROR("md_render_mol_update: the frames are not of this topology");
            return false;
        }
        const uint32_t n = atoms->frame_count;
        if (n != 1 && n != 2 && n != 4) {
            MD_LOG_ERROR("md_render_mol_update: frame_count must be 1, 2 or 4 (is %u)", n);
            return false;
        }
        for (uint32_t k = 0; k < n; ++k) {
            uint32_t sl = atoms->frame_slot[k];
            if (sl >= fr->slots || !(fr->valid & (1u << sl))) {
                MD_LOG_ERROR("md_render_mol_update: frame_slot[%u] = %u holds no frame", k, sl);
                return false;
            }
        }
    }

    gpu_mol_args_t a;
    base_args(&a, topo, dynamic);
    a.position            = atoms->position;
    a.stride              = atoms->stride ? atoms->stride : 12u;
    a.velocity            = atoms->velocity;
    a.vel_stride          = atoms->velocity_stride ? atoms->velocity_stride : 8u;
    a.has_velocity        = !fr && atoms->velocity ? 1u : 0u;
    a.radius_override     = atoms->radius;
    a.has_radius          = atoms->radius ? 1u : 0u;
    a.full                = (regroup || radii) ? 1u : 0u;

    bool ok = true;
    if (radii && atoms->radius_table && topo->type_count) ok = md_gpu_copy(s, topo->type_radius, atoms->radius_table, (size_t)topo->type_count * 4u);
    if (ok && regroup) ok = run_regroup(r, s, &a, topo, hook, user);
    if (topo->point_count) {
        if (fr) {
            a.frame_count = atoms->frame_count;
            for (uint32_t k = 0; k < 4; ++k) {
                uint32_t sl = atoms->frame_slot[k < atoms->frame_count ? k : 0];
                a.frame[k] = frame_buffer(fr, fr->buffer[sl]);
                a.fpar[k]  = frame_params(fr, sl);
                a.w[k]     = k < atoms->frame_count ? atoms->frame_w[k] : 0.0f;
                a.dw[k]    = k < atoms->frame_count ? atoms->frame_dw[k] : 0.0f;
            }
            a.has_fvel = has_velocity(atoms) ? 1u : 0u;
            ok = ok && launch(s, r->k_fr_gather, groups_for(topo->point_count, 64), &a, hook, user, "gather");
        } else {
            ok = ok && launch(s, r->k_gather, groups_for(topo->point_count, 64), &a, hook, user, "gather");
        }
        /* Every gather reads the radii from their source until a full one has
           marked the copy's w valid; after that only the flags make it full */
        ok = ok && launch(s, r->k_mark, md_gpu_grid(1, 1, 1), &a, hook, user, "mark");
    }
    if (ok && content) content[MD_RENDER_CONTENT_GEOMETRY]++;
    return ok;
}

bool md_render_mol_update(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo,
                          md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                          md_render_mol_update_flags_t flags, uint64_t content[MD_RENDER_CONTENT_COUNT]) {
    return update_impl(r, s, topo, dynamic, atoms, flags, content, NULL, NULL);
}

/* Not in the header: for md_render_bench, which declares it. The same work
   with hook(user, pass) after each launch. */
bool md_render_mol_update_hooked(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo,
                                 md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                                 md_render_mol_update_flags_t flags, void (*hook)(void* user, const char* pass), void* user) {
    return update_impl(r, s, topo, dynamic, atoms, flags, NULL, hook, user);
}

bool md_render_mol_upload_state(md_gpu_stream_t s, md_gpu_addr_t dst, const md_system_state_t* state,
                                uint64_t content[MD_RENDER_CONTENT_COUNT]) {
    if (!s || !dst || !state) return false;
    _Static_assert(sizeof(md_float3_t) == 12, "packed float3");
    if (state->num_atoms && !md_gpu_upload(s, dst, state->xyz, state->num_atoms * 12u)) return false;
    if (content) content[MD_RENDER_CONTENT_GEOMETRY]++;
    return true;
}

/* ---- Frames ---------------------------------------------------------------------- */

static inline uint16_t quant1(float v, float lo, float inv) {
    float t = (v - lo) * inv + 0.5f;
    if (!(t >= 0.0f)) t = 0.0f;                 /* NaN too */
    if (t > 65535.0f) t = 65535.0f;
    return (uint16_t)t;
}

/* The quantize goes over runs of consecutive members: a run of float3
   positions is four atoms to three vectors (each with its own pattern of
   axes: xyzx, yzxy, zxyz), the rest one by one */
typedef struct quant_ctx_t {
    const float* xyz;
    size_t       step;              /* floats from one position to the next */
    md_128       mn[3], mx[3];      /* the bounds, in the vector patterns    */
    float        lo[3], hi[3];      /* and one by one                        */
    md_128       L[3], I[3];        /* origin and 1 / scale, in the patterns */
    float        inv[3];
    uint16_t*    q;                 /* next value out                        */
} quant_ctx_t;

static inline void run_bounds(quant_ctx_t* c, uint32_t beg, uint32_t end) {
    uint32_t i = beg;
    if (c->step == 3 && end - beg >= 4) {
        md_128 mn0 = c->mn[0], mn1 = c->mn[1], mn2 = c->mn[2];     /* in registers: c may alias */
        md_128 mx0 = c->mx[0], mx1 = c->mx[1], mx2 = c->mx[2];
        for (; i + 4 <= end; i += 4) {
            const float* p = c->xyz + 3 * (size_t)i;
            md_128 v0 = md_mm_loadu_ps(p), v1 = md_mm_loadu_ps(p + 4), v2 = md_mm_loadu_ps(p + 8);
            mn0 = md_mm_min_ps(v0, mn0); mn1 = md_mm_min_ps(v1, mn1); mn2 = md_mm_min_ps(v2, mn2);   /* a NaN in the first operand is ignored */
            mx0 = md_mm_max_ps(v0, mx0); mx1 = md_mm_max_ps(v1, mx1); mx2 = md_mm_max_ps(v2, mx2);
        }
        c->mn[0] = mn0; c->mn[1] = mn1; c->mn[2] = mn2;
        c->mx[0] = mx0; c->mx[1] = mx1; c->mx[2] = mx2;
    }
    for (; i < end; ++i) {
        const float* p = c->xyz + c->step * i;
        for (int d = 0; d < 3; ++d) {                /* NaN compares false: skipped */
            if (p[d] < c->lo[d]) c->lo[d] = p[d];
            if (p[d] > c->hi[d]) c->hi[d] = p[d];
        }
    }
}

static inline void run_quant(quant_ctx_t* c, uint32_t beg, uint32_t end) {
    uint32_t i = beg;
    uint16_t* q = c->q;
    if (c->step == 3 && end - beg >= 4) {
        const md_128 half = md_mm_set1_ps(0.5f), zero = md_mm_set1_ps(0.0f), top = md_mm_set1_ps(65535.0f);
        const md_128 L0 = c->L[0], L1 = c->L[1], L2 = c->L[2], I0 = c->I[0], I1 = c->I[1], I2 = c->I[2];
        for (; i + 4 <= end; i += 4) {
            const float* p = c->xyz + 3 * (size_t)i;
            md_128 t0 = md_mm_add_ps(md_mm_mul_ps(md_mm_sub_ps(md_mm_loadu_ps(p),     L0), I0), half);
            md_128 t1 = md_mm_add_ps(md_mm_mul_ps(md_mm_sub_ps(md_mm_loadu_ps(p + 4), L1), I1), half);
            md_128 t2 = md_mm_add_ps(md_mm_mul_ps(md_mm_sub_ps(md_mm_loadu_ps(p + 8), L2), I2), half);
            t0 = md_mm_min_ps(md_mm_max_ps(t0, zero), top);              /* max first: a NaN becomes 0 */
            t1 = md_mm_min_ps(md_mm_max_ps(t1, zero), top);
            t2 = md_mm_min_ps(md_mm_max_ps(t2, zero), top);
            md_128i r0 = md_mm_cvttps_epi32(t0), r1 = md_mm_cvttps_epi32(t1), r2 = md_mm_cvttps_epi32(t2);
            md_mm_storeu_si128((md_128i*)q, simde_mm_packus_epi32(r0, r1));
            md_mm_storel_epi64((md_128i*)(q + 8), simde_mm_packus_epi32(r2, r2));
            q += 12;
        }
    }
    const float lo0 = c->lo[0], lo1 = c->lo[1], lo2 = c->lo[2], in0 = c->inv[0], in1 = c->inv[1], in2 = c->inv[2];
    for (; i < end; ++i) {
        const float* p = c->xyz + c->step * i;
        q[0] = quant1(p[0], lo0, in0); q[1] = quant1(p[1], lo1, in1); q[2] = quant1(p[2], lo2, in2);
        q += 3;
    }
    c->q = q;
}

/* fn over the maximal runs of members of s below count, in order (inlined
   with its fn: a run of three atoms is one call of nothing) */
static inline void each_run(md_bitset_t s, uint32_t count, quant_ctx_t* c, void (*fn)(quant_ctx_t*, uint32_t, uint32_t)) {
    const uint32_t last = s.end < count ? s.end : count;
    if (s.beg >= last) return;
    if (md_bitset_is_run(s)) { fn(c, s.beg, last); return; }
    uint32_t rb = 0, re = 0;                         /* the run being gathered */
    for (uint32_t base = s.beg & ~63u; base < last; base += 64) {
        uint64_t m = md_bitset_word(s, base);
        if (last - base < 64) m &= (1ull << (last - base)) - 1;
        while (m) {
            const uint32_t b = (uint32_t)md_ctz64(m);
            const uint64_t gap = ~(m >> b);          /* ones from the first non-member on */
            const uint32_t len = gap ? (uint32_t)md_ctz64(gap) : 64 - b;
            const uint32_t sb = base + b, se = sb + len;
            if (sb == re) re = se;
            else { if (re > rb) fn(c, rb, re); rb = sb; re = se; }
            m = b + len >= 64 ? 0 : m & (~0ull << (b + len));
        }
    }
    if (re > rb) fn(c, rb, re);
}

void md_render_mol_frame_quantize(md_render_mol_frame_t* out, uint16_t* q, const float* xyz, uint32_t stride,
                                  uint32_t count, const float cell[9], const md_bitset_t* resident) {
    if (!out) return;
    memset(out, 0, sizeof(*out));
    out->q = q;
    if (cell) memcpy(out->cell, cell, sizeof(out->cell));
    if (!q || (count && !xyz)) return;
    const md_bitset_t set = resident ? *resident : md_bitset_range(0, count);
    quant_ctx_t c;
    memset(&c, 0, sizeof(c));
    c.xyz  = xyz;
    c.step = (stride ? stride : 12u) / sizeof(float);
    for (int k = 0; k < 3; ++k) {
        c.mn[k] = md_mm_set1_ps(FLT_MAX); c.mx[k] = md_mm_set1_ps(-FLT_MAX);
        c.lo[k] = FLT_MAX; c.hi[k] = -FLT_MAX;
    }
    each_run(set, count, &c, run_bounds);
    float a[3][4], b[3][4];
    for (int k = 0; k < 3; ++k) { md_mm_storeu_ps(a[k], c.mn[k]); md_mm_storeu_ps(b[k], c.mx[k]); }
    for (int k = 0; k < 3; ++k) for (int l = 0; l < 4; ++l) {
        const int d = (4 * k + l) % 3;
        if (a[k][l] < c.lo[d]) c.lo[d] = a[k][l];
        if (b[k][l] > c.hi[d]) c.hi[d] = b[k][l];
    }
    for (int d = 0; d < 3; ++d) {
        if (!(c.lo[d] <= c.hi[d])) c.lo[d] = c.hi[d] = 0.0f;
        out->origin[d] = c.lo[d];
        out->scale[d]  = (c.hi[d] - c.lo[d]) / 65535.0f;
        c.inv[d] = out->scale[d] > 0.0f ? 65535.0f / (c.hi[d] - c.lo[d]) : 0.0f;
    }
    for (int k = 0; k < 3; ++k) {
        c.L[k] = md_mm_set_ps(c.lo[(4 * k + 3) % 3], c.lo[(4 * k + 2) % 3], c.lo[(4 * k + 1) % 3], c.lo[(4 * k) % 3]);
        c.I[k] = md_mm_set_ps(c.inv[(4 * k + 3) % 3], c.inv[(4 * k + 2) % 3], c.inv[(4 * k + 1) % 3], c.inv[(4 * k) % 3]);
    }
    c.q = q;
    each_run(set, count, &c, run_quant);
}

uint64_t md_render_mol_frame_size(const md_render_mol_topo_t* topo) {
    return topo ? frame_bytes(topo->point_count) : 0;
}

void md_render_mol_frame_cell(float out[9], const md_unitcell_t* cell) {
    if (!out) return;
    memset(out, 0, 9 * sizeof(float));
    if (!cell) return;
    float A[3][3];
    md_unitcell_A_extract_float(A, cell);
    const md_unitcell_flags_t pbc[3] = { MD_UNITCELL_PBC_X, MD_UNITCELL_PBC_Y, MD_UNITCELL_PBC_Z };
    for (int c = 0; c < 3; ++c) {
        if (!(cell->flags & pbc[c])) continue;
        for (int k = 0; k < 3; ++k) out[3 * c + k] = A[c][k];
    }
}

/* A frame's 32 words on the GPU: origin, scale, the cell's columns, the rows
   of its inverse (zero for an axis that is not periodic), periodic */
static void frame_param_words(uint32_t out[32], const md_render_mol_frame_t* f) {
    float o[32];
    memset(o, 0, sizeof(o));
    for (int d = 0; d < 3; ++d) { o[0 + d] = f->origin[d]; o[3 + d] = f->scale[d]; }
    memcpy(o + 6, f->cell, 9 * sizeof(float));
    bool col[3];
    double M[3][3];                 /* columns; a missing one replaced by its axis */
    for (int c = 0; c < 3; ++c) {
        col[c] = f->cell[3 * c] != 0.0f || f->cell[3 * c + 1] != 0.0f || f->cell[3 * c + 2] != 0.0f;
        for (int k = 0; k < 3; ++k) M[c][k] = col[c] ? f->cell[3 * c + k] : (k == c ? 1.0 : 0.0);
    }
    uint32_t periodic = 0;
    if (col[0] || col[1] || col[2]) {
        /* rows of the inverse: the cofactors over the determinant */
        double r[3][3];
        r[0][0] = M[1][1] * M[2][2] - M[1][2] * M[2][1]; r[0][1] = M[1][2] * M[2][0] - M[1][0] * M[2][2]; r[0][2] = M[1][0] * M[2][1] - M[1][1] * M[2][0];
        r[1][0] = M[2][1] * M[0][2] - M[2][2] * M[0][1]; r[1][1] = M[2][2] * M[0][0] - M[2][0] * M[0][2]; r[1][2] = M[2][0] * M[0][1] - M[2][1] * M[0][0];
        r[2][0] = M[0][1] * M[1][2] - M[0][2] * M[1][1]; r[2][1] = M[0][2] * M[1][0] - M[0][0] * M[1][2]; r[2][2] = M[0][0] * M[1][1] - M[0][1] * M[1][0];
        const double det = M[0][0] * r[0][0] + M[0][1] * r[0][1] + M[0][2] * r[0][2];
        double scale = 0.0;
        for (int c = 0; c < 3; ++c) for (int k = 0; k < 3; ++k) scale = M[c][k] * M[c][k] > scale ? M[c][k] * M[c][k] : scale;
        if (det * det > 1e-12 * scale * scale * scale) {
            for (int i = 0; i < 3; ++i) {
                if (!col[i]) continue;
                for (int k = 0; k < 3; ++k) o[15 + 3 * i + k] = (float)(r[i][k] / det);
            }
            periodic = 1;
        }
    }
    memcpy(out, o, sizeof(o));
    out[24] = periodic;
}

uint64_t md_render_mol_frames_size(const md_render_mol_topo_t* topo, uint32_t slots) {
    if (!topo || slots < 1 || slots > MD_RENDER_MOL_FRAME_SLOTS) return 0;
    return frames_layout(topo->point_count, topo->small_atom_count, slots).size;
}

bool md_render_mol_frames_init(md_render_mol_frames_t* out, const md_render_mol_topo_t* topo, md_gpu_addr_t base, uint32_t slots) {
    if (!out || !topo || !base || slots < 1 || slots > MD_RENDER_MOL_FRAME_SLOTS) {
        MD_LOG_ERROR("md_render_mol_frames_init: missing block or topology, or slots not in [1, %d]", MD_RENDER_MOL_FRAME_SLOTS);
        return false;
    }
    memset(out, 0, sizeof(*out));
    out->base             = base;
    out->size             = md_render_mol_frames_size(topo, slots);
    out->topology_base    = topo->base;
    out->topology_version = topo->topology_version;
    out->slots            = slots;
    out->point_count      = topo->point_count;
    for (uint32_t b = 0; b < slots; ++b) out->buffer[b] = b;
    return true;
}

static bool frame_upload_impl(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic,
                              md_render_mol_frames_t* fr, uint32_t slot, const md_render_mol_frame_t* f,
                              md_render_mol_update_flags_t flags, pass_hook_t hook, void* user) {
    if (!r || !s || !topo || !topo->base || !fr || !f || (topo->point_count && !f->q && !f->q_gpu)) {
        MD_LOG_ERROR("md_render_mol_frame_upload: missing context, stream, topology, frames or frame data");
        return false;
    }
    if (!frames_match(fr, topo) || slot >= fr->slots) {
        MD_LOG_ERROR("md_render_mol_frame_upload: frames not of this topology, or slot %u not in [0, %u)", slot, fr->slots);
        return false;
    }
    if (f->q_gpu & 3u) {
        MD_LOG_ERROR("md_render_mol_frame_upload: q_gpu must be 4-byte aligned");
        return false;
    }
    const bool regroup = (flags & MD_RENDER_MOL_UPDATE_REGROUP) && topo->grid_cells;
    if (regroup && !dynamic) {
        MD_LOG_ERROR("md_render_mol_frame_upload: a regroup needs the per-frame block (its scratch)");
        return false;
    }
    const uint32_t n = topo->point_count;
    const frames_layout_t L = frames_layout(n, topo->small_atom_count, fr->slots);
    uint32_t par[32];
    frame_param_words(par, f);
    bool ok = md_gpu_upload(s, frame_params(fr, slot), par, sizeof(par));
    if (!ok || !n) {
        if (ok) fr->valid |= 1u << slot;
        return ok;
    }

    /* The frame where the GPU reads it: the caller's, or a copy of the host's,
       freed in stream order once the reorder has read it */
    md_gpu_addr_t src = f->q_gpu, tmp = 0;
    if (!src) {
        const size_t bytes = (size_t)n * 6u, head = bytes & ~(size_t)3u;
        tmp = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, frame_bytes(n)).gpu;
        if (!tmp) return false;
        ok = md_gpu_upload(s, tmp, f->q, head);
        if (ok && head < bytes) {                  /* an odd count: the last value in a word of its own */
            uint32_t last = f->q[3 * (size_t)n - 1];
            ok = md_gpu_upload(s, tmp + head, &last, 4);
        }
        src = tmp;
    }

    gpu_mol_args_t a;
    base_args(&a, topo, dynamic);
    uint32_t free_buf = fr->buffer[slot];          /* its frame is replaced: free */
    if (ok && regroup) {
        const uint32_t others = fr->valid & ~(1u << slot);
        a.quantized = 1;
        a.qsrc      = src;
        a.qpar      = frame_params(fr, slot);
        a.has_src   = others ? 1u : 0u;
        a.src_point = fr->base + L.src_point;
        ok = run_regroup(r, s, &a, topo, hook, user);
        a.quantized = 0;
        for (uint32_t j = 0; ok && j < fr->slots; ++j) {
            if (!(others & (1u << j))) continue;
            a.qsrc = frame_buffer(fr, fr->buffer[j]);
            a.qdst = frame_buffer(fr, free_buf);
            ok = launch(s, r->k_fr_permute, groups_for((n + 1) / 2, 64), &a, hook, user, "permute");
            uint32_t t = fr->buffer[j]; fr->buffer[j] = free_buf; free_buf = t;
        }
    }
    a.qsrc = src;
    a.qdst = frame_buffer(fr, free_buf);
    ok = ok && launch(s, r->k_fr_reorder, groups_for((n + 1) / 2, 64), &a, hook, user, "reorder");
    fr->buffer[slot] = free_buf;
    if (tmp) md_gpu_free(s, tmp);
    if (ok) fr->valid |= 1u << slot;
    else    fr->valid = 0;                         /* the slots' contents are unknown */
    return ok;
}

bool md_render_mol_frame_upload(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic,
                                md_render_mol_frames_t* fr, uint32_t slot, const md_render_mol_frame_t* f,
                                md_render_mol_update_flags_t flags) {
    return frame_upload_impl(r, s, topo, dynamic, fr, slot, f, flags, NULL, NULL);
}

/* Not in the header: for md_render_bench, as md_render_mol_update_hooked */
bool md_render_mol_frame_upload_hooked(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic,
                                       md_render_mol_frames_t* fr, uint32_t slot, const md_render_mol_frame_t* f,
                                       md_render_mol_update_flags_t flags, void (*hook)(void* user, const char* pass), void* user) {
    return frame_upload_impl(r, s, topo, dynamic, fr, slot, f, flags, hook, user);
}

void md_render_mol_frame_weights(float w[4], float dw[4], uint32_t frame_count, float t, float dt, float tension) {
    if (!w || !dw) return;
    memset(w, 0, 4 * sizeof(float));
    memset(dw, 0, 4 * sizeof(float));
    t = t < 0.0f ? 0.0f : t > 1.0f ? 1.0f : t;
    const float id = dt > 0.0f ? 1.0f / dt : 0.0f;
    if (frame_count == 2) {
        w[0] = 1.0f - t; w[1] = t;
        dw[0] = -id; dw[1] = id;
    } else if (frame_count == 4) {
        /* Hermite: p = h00 p1 + h10 m1 + h01 p2 + h11 m2, m1 = s (p2 - p0), m2 = s (p3 - p1) */
        const float s = tension, t2 = t * t, t3 = t2 * t;
        const float h00 =  2.0f * t3 - 3.0f * t2 + 1.0f, h10 = t3 - 2.0f * t2 + t;
        const float h01 = -2.0f * t3 + 3.0f * t2,        h11 = t3 - t2;
        const float d00 =  6.0f * t2 - 6.0f * t,         d10 = 3.0f * t2 - 4.0f * t + 1.0f;
        const float d01 = -6.0f * t2 + 6.0f * t,         d11 = 3.0f * t2 - 2.0f * t;
        w[0] = -s * h10;  w[1] = h00 - s * h11;  w[2] = h01 + s * h10;  w[3] = s * h11;
        dw[0] = -s * d10 * id; dw[1] = (d00 - s * d11) * id; dw[2] = (d01 + s * d10) * id; dw[3] = s * d11 * id;
    } else {
        w[0] = 1.0f;
    }
}

/* ---- Emit ------------------------------------------------------------------------ */

/* The copy in the per-frame block: points in cluster order, each its atom's
   attributes through atom_order */
static md_render_points_t atom_points(const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms, const md_render_mol_rep_t* rep) {
    const dyn_layout_t D = dyn_layout(topo->point_count, topo->small_molecule_count, topo->grid_cells);
    md_render_points_t p;
    memset(&p, 0, sizeof(p));
    p.count           = topo->point_count;
    p.position        = dynamic + D.points;
    p.stride          = 16;
    p.radius          = dynamic + D.points + 12;
    p.radius_stride   = 16;
    p.velocity        = has_velocity(atoms) ? dynamic + D.velocity : 0;
    p.velocity_stride = 8;
    p.attr_count      = topo->atom_count;
    p.attr_index      = topo->atom_order;
    p.color           = rep->atom_color;
    p.palette_index   = rep->atom_palette_index;
    p.visible[0]      = rep->visible;
    p.visible[1]      = atoms->visible;
    p.state[0]        = atoms->selection;
    p.state[1]        = atoms->highlight;
    return p;
}

static void put(md_render_batch_t* out, uint32_t cap, uint32_t* n, const md_render_batch_t* b) {
    if (*n < cap) out[*n] = *b;
    (*n)++;
}

uint32_t md_render_mol_emit(const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                            const md_render_mol_rep_t* rep, md_render_batch_t* out, uint32_t cap) {
    if (!topo || !atoms || !rep || (cap && !out)) return 0;
    if (!dynamic) {
        MD_LOG_ERROR("md_render_mol: no per-frame block (the copy md_render_mol_update writes)");
        return 0;
    }
    if (rep->key == 0 || rep->key >= (1ull << 56)) {
        MD_LOG_ERROR("md_render_mol: representation key must be in [1, 2^56)");
        return 0;
    }
    md_render_batch_t base;
    memset(&base, 0, sizeof(base));
    base.points          = atom_points(topo, dynamic, atoms, rep);
    base.color           = rep->color;
    base.transform       = rep->transform;
    base.transform_count = rep->transform_count;

    bool  spheres = false, sticks = false;
    float sphere_scale = 0.0f, sphere_add = 0.0f, stick_radius = 0.0f;
    uint32_t classes = rep->bond_classes;
    const uint32_t default_classes = MD_RENDER_MOL_BONDS(MD_RENDER_MOL_BOND_COVALENT) | MD_RENDER_MOL_BONDS(MD_RENDER_MOL_BOND_COORDINATE);

    switch (rep->type) {
    case MD_RENDER_MOL_REP_SPACEFILL:
        spheres = true;
        sphere_scale = rep->params.spacefill.radius_scale;
        break;
    case MD_RENDER_MOL_REP_LICORICE:
        spheres = sticks = true;
        sphere_add = stick_radius = rep->params.licorice.radius;
        if (!classes) classes = default_classes;
        break;
    case MD_RENDER_MOL_REP_BALL_AND_STICK:
        spheres = sticks = true;
        sphere_scale = rep->params.ball_and_stick.ball_scale;
        stick_radius = rep->params.ball_and_stick.stick_radius;
        if (!classes) classes = default_classes;
        break;
    case MD_RENDER_MOL_REP_BONDS:
        sticks = true;
        stick_radius = rep->params.bonds.radius;
        if (!classes) classes = MD_RENDER_MOL_BONDS_ALL;
        break;
    default:
        MD_LOG_ERROR("md_render_mol: representation type %d is not implemented yet", (int)rep->type);
        return 0;
    }

    /* At licorice joints the atom's sphere and the stick's round end are the
       same surface; md_render gives such ties to the sphere, so a joint picks
       as its atom whatever the order here. */
    uint32_t n = 0;
    if (sticks) {
        for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) {
            if (!(classes & MD_RENDER_MOL_BONDS(c))) continue;
            uint32_t beg = topo->bond_class_offset[c], end = topo->bond_class_offset[c + 1];
            if (end <= beg) continue;
            md_render_batch_t b = base;
            b.key   = (rep->key << 8) | (MD_RENDER_MOL_PART_BONDS + c);
            b.kind  = MD_RENDER_CAPSULES;
            b.count = end - beg;
            b.index = topo->bond_pairs + (md_gpu_addr_t)beg * 8u;
            b.id    = topo->bond_ids   + (md_gpu_addr_t)beg * 4u;
            b.cluster_offset = topo->bond_cluster_offset[c];
            b.cluster_count  = topo->bond_cluster_count[c];
            b.style.capsules.radius = stick_radius;
            put(out, cap, &n, &b);
        }
    }
    if (spheres && topo->point_count) {
        md_render_batch_t b = base;
        b.key   = (rep->key << 8) | MD_RENDER_MOL_PART_ATOMS;
        b.kind  = MD_RENDER_SPHERES;
        b.count = topo->point_count;         /* the points are in cluster order: no index */
        b.cluster_offset = topo->atom_cluster_offset;
        b.cluster_count  = topo->atom_cluster_count;
        b.style.spheres.radius_scale = sphere_scale;
        b.style.spheres.radius_add   = sphere_add;
        put(out, cap, &n, &b);
    }
    return n;
}

bool md_render_mol_decode_pick(const md_render_pick_t* pick, md_render_mol_hit_t* hit) {
    if (!pick || !hit) return false;
    memset(hit, 0, sizeof(*hit));
    hit->bond = UINT32_MAX;
    if (pick->key == 0) return false;
    uint32_t part = (uint32_t)(pick->key & 0xFFu);
    hit->rep_key  = pick->key >> 8;
    hit->instance = pick->instance;
    memcpy(hit->position, pick->position, sizeof(hit->position));
    if (part == MD_RENDER_MOL_PART_ATOMS && pick->kind == MD_RENDER_SPHERES) {
        hit->kind = MD_RENDER_MOL_HIT_ATOM;
        hit->atom = pick->attr;
        return true;
    }
    if (part >= MD_RENDER_MOL_PART_BONDS && part < MD_RENDER_MOL_PART_BONDS + MD_RENDER_MOL_BOND_CLASS_COUNT && pick->kind == MD_RENDER_CAPSULES) {
        hit->kind = MD_RENDER_MOL_HIT_BOND;
        hit->bond = pick->id;
        hit->atom = pick->attr;
        return true;
    }
    return false;
}

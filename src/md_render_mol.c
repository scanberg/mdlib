#include "md_render_mol.h"

#include <md_system.h>
#include <md_bond.h>
#include <md_system_state.h>
#include <core/md_log.h>

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
    uint64_t atom_type, atom_radius, atom_order, atom_point, atom_cluster_offset, bond_pairs, bond_ids;
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
    uint64_t points, velocity, marker, params, mol_slot, cnt, sums, size;
    uint32_t scan_n, scan_blocks;
} dyn_layout_t;

/* The points (float4: position, radius) and their velocities in cluster
   order: what is drawn. Then the scratch of the regroup. */
static dyn_layout_t dyn_layout(uint32_t n_atoms, uint32_t n_mol, uint32_t cells) {
    dyn_layout_t L;
    L.scan_n      = cells + 1;
    L.scan_blocks = (L.scan_n + SCAN_BLOCK - 1) / SCAN_BLOCK;
    uint64_t o = 0;
    L.points   = o; o = align256(o + (uint64_t)n_atoms * 16u);
    L.velocity = o; o = align256(o + (uint64_t)n_atoms * 8u);
    L.marker   = o; o = align256(o + 16u);
    L.params   = o; o = align256(o + 16u * 4u);
    L.mol_slot = o; o = align256(o + (uint64_t)n_mol * 16u);
    L.cnt      = o; o = align256(o + (uint64_t)L.scan_n * 16u);
    L.sums     = o; o = align256(o + (uint64_t)L.scan_blocks * 16u);
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

bool md_render_mol_plan(md_render_mol_plan_t* out, const md_system_t* sys, md_temp_scope_t temp) {
    if (!out || !sys) return false;
    memset(out, 0, sizeof(*out));
    md_render_mol_topo_t* T = &out->layout;
    const uint32_t n  = (uint32_t)sys->atom.count;
    const uint32_t nb = (uint32_t)(sys->bond.count + sys->bond.noncovalent_count);
    T->topology_version = sys->topology_version;
    T->atom_count       = n;
    T->bond_count       = nb;

    /* Molecules, and which are small */
    const md_molecule_data_t* M = &sys->molecule;
    const uint32_t n_mol = (uint32_t)M->count;
    int32_t*  mol_of = md_temp_alloc(temp, ((size_t)n + 1) * sizeof(int32_t));
    uint8_t*  small  = md_temp_alloc_zero(temp, (size_t)n_mol + 1);
    if (!mol_of || !small) return false;
    for (uint32_t i = 0; i < n; ++i) mol_of[i] = -1;
    uint32_t small_mols = 0, small_atoms = 0;
    for (uint32_t m = 0; m < n_mol; ++m) {
        uint32_t b = M->offset[m], e = M->offset[m + 1];
        for (uint32_t k = b; k < e; ++k) {
            int32_t a = M->atom_idx[k];
            if (a >= 0 && (uint32_t)a < n) mol_of[a] = (int32_t)m;
        }
        if (e > b && e - b <= MD_RENDER_MOL_SMALL) { small[m] = 1; small_mols++; small_atoms += e - b; }
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
        for (uint32_t i = b; i < e; ++i) if (!placed[i] && !IS_SMALL(i)) { placed[i] = 1; tmp[k++] = i; }
        if (!k) continue;
        int32_t mol = mol_of[tmp[0]];
        if (mol != prev_mol) pack_close(&P);
        prev_mol = mol;
        pack_unit(&P, tmp, k);
    }
    pack_close(&P);
    for (uint32_t i = 0; i < n; ++i) if (!placed[i] && !IS_SMALL(i)) { placed[i] = 1; pack_unit(&P, &i, 1); }
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
            if (!small[m]) continue;
            sm_index[m] = k;
            sm_aoff[k++] = a;
            for (uint32_t s = M->offset[m]; s < M->offset[m + 1]; ++s) {
                local[M->atom_idx[s]] = (uint8_t)(s - M->offset[m]);
                sm_atoms[a++] = (uint32_t)M->atom_idx[s];
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
        cls[i] = bond_class(sys->bond.flags ? sys->bond.flags[i] : 0);
        uint32_t a = (uint32_t)sys->bond.pairs[i].idx[0], b = (uint32_t)sys->bond.pairs[i].idx[1];
        if (a >= n || b >= n) { key[i] = NONE; continue; }
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
        uint32_t k = sm_index[mol_of[sys->bond.pairs[i].idx[0]]];
        uint32_t j = sm_boff[k] + sm_bfill[k]++;
        sm_blocal[j] = (uint32_t)local[sys->bond.pairs[i].idx[0]] | ((uint32_t)local[sys->bond.pairs[i].idx[1]] << 16);
        sm_bids[j] = i;
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
    for (uint32_t k = 0; k < n; ++k) inv[order[k]] = k;
    uint32_t n_mixed = 0;
    for (uint32_t i = 0; i < nb; ++i) {
        uint32_t a = (uint32_t)sys->bond.pairs[i].idx[0], b = (uint32_t)sys->bond.pairs[i].idx[1];
        if (key[i] != SMALL_KEY && a < n && b < n && (IS_SMALL(a) || IS_SMALL(b))) n_mixed++;
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
    L.atom_type           = o; o = align256(o + ((uint64_t)n + 1) / 2 * 4u);
    L.atom_radius         = o; o = align256(o + (uint64_t)n * 4u);
    L.atom_order          = o; o = align256(o + (uint64_t)n * 4u);
    L.atom_point          = o; o = align256(o + (uint64_t)n * 4u);
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
        float*    radius = (float*)(img + L.atom_radius);
        for (uint32_t i = 0; i < n; ++i) {
            uint16_t t = sys->atom.type_idx ? sys->atom.type_idx[i] : 0;
            type[i]   = t;
            radius[i] = (PT->radius && t < PT->count) ? PT->radius[t] : 0.0f;
        }
        T->type_count = (uint32_t)PT->count;
    }
    memcpy(img + L.atom_order, order, (size_t)n * 4u);
    memcpy(img + L.atom_point, inv, (size_t)n * 4u);
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
                uint32_t a = (uint32_t)sys->bond.pairs[i].idx[0], b = (uint32_t)sys->bond.pairs[i].idx[1];
                bool ok = a < n && b < n;
                pairs[2 * (beg + j) + 0] = ok ? inv[a] : 0;
                pairs[2 * (beg + j) + 1] = ok ? inv[b] : 0;
                ids[beg + j] = i;
                if (ok && ((mol_of[a] >= 0 && small[mol_of[a]]) || (mol_of[b] >= 0 && small[mol_of[b]]))) {
                    mixed_pos[nm] = beg + j;
                    mixed_atoms[2 * nm + 0] = a;
                    mixed_atoms[2 * nm + 1] = b;
                    nm++;
                }
            }
        }
        /* the small molecules' bonds, in the initial grouping: molecule order */
        for (uint32_t k = 0; k < small_bonds; ++k) {
            uint32_t i = sm_bids[k];
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
    T->dynamic_size        = dyn_layout(n, small_mols, cells).size;
    T->atom_type           = L.atom_type;
    T->atom_radius         = L.atom_radius;
    T->atom_order          = L.atom_order;
    T->atom_point          = L.atom_point;
    T->atom_cluster_offset = L.atom_cluster_offset;
    T->bond_pairs          = L.bond_pairs;
    T->bond_ids            = L.bond_ids;
    for (uint32_t c = 0; c < MD_RENDER_MOL_BOND_CLASS_COUNT; ++c) T->bond_cluster_offset[c] = L.bond_cluster_offset[c];
    T->small_tables        = small_mols ? L.small_tables : 0;
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
    out->atom_radius         = base + plan->layout.atom_radius;
    out->atom_order          = base + plan->layout.atom_order;
    out->atom_point          = base + plan->layout.atom_point;
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
    md_gpu_kernel_t k_clear, k_bounds, k_dims, k_count, k_scan, k_write, k_offsets, k_radii, k_mixed, k_gather, k_mark;
};

/* Mirrors Args of md_render_mol.slang */
typedef struct gpu_mol_args_t {
    md_gpu_addr_t position, velocity, radius_override, points, vel_out;
    md_gpu_addr_t params, mol_slot, cnt, sums, types, radius_table, atom_radius;
    md_gpu_addr_t mol_atom_offset, mol_atoms, mol_bond_offset, mol_bond_local, mol_bond_ids, mixed_pos, mixed_atoms;
    md_gpu_addr_t atom_order, atom_point, atom_cluster_offset, bond_pairs, bond_ids, bond_cluster_offset, marker;
    uint32_t stride, n_mol, cells, n_atoms, static_atoms, static_atom_clusters, static_bonds, static_bond_clusters;
    uint32_t has_bonds, mode, n, atom_count, type_count, vel_stride, has_velocity, has_radius, mixed_count;
    uint32_t full, ident0, ident1, ident2, _pad;
} gpu_mol_args_t;
_Static_assert(sizeof(gpu_mol_args_t) == 26 * 8 + 22 * 4, "gpu_mol_args_t mirrors Args");

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
    kd = md_shader_md_render_mol_sol_radii_kernel();   r->k_radii   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_mixed_kernel();   r->k_mixed   = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_gather_kernel();  r->k_gather  = md_gpu_kernel_create(device, &kd);
    kd = md_shader_md_render_mol_sol_mark_kernel();    r->k_mark    = md_gpu_kernel_create(device, &kd);
    if (!r->k_clear || !r->k_bounds || !r->k_dims || !r->k_count || !r->k_scan || !r->k_write || !r->k_offsets || !r->k_radii ||
        !r->k_mixed || !r->k_gather || !r->k_mark) {
        MD_LOG_ERROR("md_render_mol: kernel creation failed: %s", md_gpu_last_error());
        md_render_mol_destroy(r);
        return NULL;
    }
    return r;
}

void md_render_mol_destroy(md_render_mol_t r) {
    if (!r) return;
    md_gpu_kernel_t* k[] = { &r->k_clear, &r->k_bounds, &r->k_dims, &r->k_count, &r->k_scan, &r->k_write, &r->k_offsets, &r->k_radii,
                             &r->k_mixed, &r->k_gather, &r->k_mark };
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

static bool update_impl(md_render_mol_t r, md_gpu_stream_t s, const md_render_mol_topo_t* topo,
                        md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                        md_render_mol_update_flags_t flags, uint64_t content[MD_RENDER_CONTENT_COUNT],
                        pass_hook_t hook, void* user) {
    if (!r || !s || !topo || !atoms || !topo->base || !dynamic || !atoms->position) {
        MD_LOG_ERROR("md_render_mol_update: missing context, stream, topology, per-frame block or positions");
        return false;
    }
    const bool regroup = (flags & MD_RENDER_MOL_UPDATE_REGROUP) && topo->grid_cells;
    const bool radii   = (flags & MD_RENDER_MOL_UPDATE_RADII) != 0;
    const uint32_t cells = topo->grid_cells, n_mol = topo->small_molecule_count;
    const dyn_layout_t   D  = dyn_layout(topo->atom_count, n_mol, cells);
    const small_layout_t SL = small_layout(n_mol, topo->small_atom_count, topo->small_bond_count, topo->mixed_bond_count);

    gpu_mol_args_t a;
    memset(&a, 0, sizeof(a));
    a.position            = atoms->position;
    a.stride              = atoms->stride ? atoms->stride : 12u;
    a.velocity            = atoms->velocity;
    a.vel_stride          = atoms->velocity_stride ? atoms->velocity_stride : 8u;
    a.has_velocity        = atoms->velocity ? 1u : 0u;
    a.radius_override     = atoms->radius;
    a.has_radius          = atoms->radius ? 1u : 0u;
    a.points              = dynamic + D.points;
    a.vel_out             = dynamic + D.velocity;
    a.atom_count          = topo->atom_count;
    a.type_count          = topo->type_count;
    a.types               = topo->atom_type;
    a.radius_table        = atoms->radius_table;
    a.atom_radius         = topo->atom_radius;
    a.params              = dynamic + D.params;
    a.mol_slot            = dynamic + D.mol_slot;
    a.cnt                 = dynamic + D.cnt;
    a.sums                = dynamic + D.sums;
    if (topo->small_tables) {
        a.mol_atom_offset = topo->small_tables + SL.mol_atom_offset;
        a.mol_atoms       = topo->small_tables + SL.mol_atoms;
        a.mol_bond_offset = topo->small_tables + SL.mol_bond_offset;
        a.mol_bond_local  = topo->small_tables + SL.mol_bond_local;
        a.mol_bond_ids    = topo->small_tables + SL.mol_bond_ids;
        a.mixed_pos       = topo->small_tables + SL.mixed_pos;
        a.mixed_atoms     = topo->small_tables + SL.mixed_atoms;
    }
    a.atom_order          = topo->atom_order;
    a.atom_point          = topo->atom_point;
    a.atom_cluster_offset = topo->atom_cluster_offset;
    a.bond_pairs          = topo->bond_pairs;      /* the covalent class is first: its bonds start at 0 */
    a.bond_ids            = topo->bond_ids;
    a.bond_cluster_offset = topo->bond_cluster_offset[MD_RENDER_MOL_BOND_COVALENT];
    a.n_mol               = n_mol;
    a.cells               = cells;
    a.n_atoms             = topo->small_atom_count;
    a.static_atoms        = topo->atom_count - topo->small_atom_count;
    a.static_atom_clusters= topo->atom_cluster_count - cells;
    a.has_bonds           = topo->small_bond_count ? 1u : 0u;
    a.static_bonds        = (topo->bond_class_offset[MD_RENDER_MOL_BOND_COVALENT + 1] - topo->bond_class_offset[MD_RENDER_MOL_BOND_COVALENT]) - topo->small_bond_count;
    a.static_bond_clusters= a.has_bonds ? topo->bond_cluster_count[MD_RENDER_MOL_BOND_COVALENT] - cells : 0;
    a.mixed_count         = topo->mixed_bond_count;
    a.n                   = D.scan_n;
    a.marker              = dynamic + D.marker;
    a.full                = (regroup || radii) ? 1u : 0u;
    a.ident0              = (uint32_t)topo->base;
    a.ident1              = (uint32_t)(topo->base >> 32);
    a.ident2              = (uint32_t)topo->topology_version;

    bool ok = true;
    if (radii && atoms->radius_table && topo->atom_count) ok = launch(s, r->k_radii, groups_for(topo->atom_count, 64), &a, hook, user, "radii");
    if (ok && regroup) {
        ok = ok && launch(s, r->k_clear,  groups_for(4 * (cells + 1) > 16 ? 4 * (cells + 1) : 16, 64), &a, hook, user, "clear");
        ok = ok && launch(s, r->k_bounds, groups_for(n_mol, 64), &a, hook, user, "bounds");
        ok = ok && launch(s, r->k_dims,   md_gpu_grid(1, 1, 1), &a, hook, user, "dims");
        ok = ok && launch(s, r->k_count,  groups_for(n_mol, 64), &a, hook, user, "count");
        a.mode = 0;
        ok = ok && launch(s, r->k_scan, md_gpu_grid(D.scan_blocks, 1, 1), &a, hook, user, "scan");
        if (D.scan_blocks > 1) {
            a.mode = 1;
            ok = ok && launch(s, r->k_scan, md_gpu_grid(1, 1, 1), &a, hook, user, "scan sums");
            a.mode = 2;
            ok = ok && launch(s, r->k_scan, md_gpu_grid(D.scan_blocks, 1, 1), &a, hook, user, "scan add");
        }
        ok = ok && launch(s, r->k_write,   groups_for(n_mol, 64), &a, hook, user, "write");
        ok = ok && launch(s, r->k_offsets, groups_for(cells + 1, 64), &a, hook, user, "offsets");
        if (topo->mixed_bond_count) ok = ok && launch(s, r->k_mixed, groups_for(topo->mixed_bond_count, 64), &a, hook, user, "mixed");
    }
    if (topo->atom_count) {
        ok = ok && launch(s, r->k_gather, groups_for(topo->atom_count, 64), &a, hook, user, "gather");
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

/* ---- Emit ------------------------------------------------------------------------ */

/* The copy in the per-frame block: points in cluster order, each its atom's
   attributes through atom_order */
static md_render_points_t atom_points(const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms, const md_render_mol_rep_t* rep) {
    const dyn_layout_t D = dyn_layout(topo->atom_count, topo->small_molecule_count, topo->grid_cells);
    md_render_points_t p;
    memset(&p, 0, sizeof(p));
    p.count           = topo->atom_count;
    p.position        = dynamic + D.points;
    p.stride          = 16;
    p.radius          = dynamic + D.points + 12;
    p.radius_stride   = 16;
    p.velocity        = atoms->velocity ? dynamic + D.velocity : 0;
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
    if (spheres && topo->atom_count) {
        md_render_batch_t b = base;
        b.key   = (rep->key << 8) | MD_RENDER_MOL_PART_ATOMS;
        b.kind  = MD_RENDER_SPHERES;
        b.count = topo->atom_count;          /* the points are in cluster order: no index */
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

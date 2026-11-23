// The system as the selectors see it, the static values of sets, and the selectors themselves.

#include "select.h"
#include <md_vial.h>

#include <md_system.h>
#include <md_bond.h>
#include <md_attributes.h>
#include <md_fragment.h>
#include <md_element.h>
#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NONE32 UINT32_MAX

// ### THE SYSTEM ###

static int cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return x < y ? -1 : x > y;
}

static uint8_t part_kind(const md_system_t* sys, uint32_t p) {
    const md_urange_t res = md_system_part_residue_range(sys, p);
    switch (md_system_part_kind(sys, p)) {
    case MD_PART_KIND_PEPTIDE:  return VIAL_PART_PEPTIDE;
    case MD_PART_KIND_DNA:      return VIAL_PART_DNA;
    case MD_PART_KIND_RNA:      return VIAL_PART_RNA;
    case MD_PART_KIND_NUCLEIC:  return VIAL_PART_NUCLEIC;
    case MD_PART_KIND_POLYMER:  return VIAL_PART_POLYMER;
    case MD_PART_KIND_BRANCHED: return VIAL_PART_GLYCAN;
    case MD_PART_KIND_WATER:    return VIAL_PART_WATER;
    default: break;
    }
    // Not a polymer: a lipid when its residues are, an ion when it is one, else a ligand
    bool lipid = res.end > res.beg;
    for (uint32_t r = res.beg; r < res.end; ++r) lipid = lipid && md_system_residue_kind(sys, r) == MD_RESIDUE_KIND_LIPID;
    if (lipid) return VIAL_PART_LIPID;
    if (res.end == res.beg + 1 && md_system_residue_kind(sys, res.beg) == MD_RESIDUE_KIND_ION) return VIAL_PART_ION;
    if (res.end == res.beg + 1 && md_system_residue_kind(sys, res.beg) == MD_RESIDUE_KIND_WATER) return VIAL_PART_WATER;
    return VIAL_PART_LIGAND;
}

static uint32_t atom_res(const vial_sys_t* s, uint32_t a);

int32_t vial_sys_atom_res(const vial_sys_t* s, uint32_t a) {
    const uint32_t r = a < s->num_atoms ? atom_res(s, a) : NONE32;
    return r != NONE32 && r < s->num_res ? (int32_t)r : -1;
}

static void need_res_segs(vial_sys_t* s) {
    if (s->have_res_segs) return;
    s->have_res_segs = true;
    s->res_seg = md_alloc(s->alloc, (s->num_res + 1) * sizeof(int32_t));
    for (uint32_t r = 0; r < s->num_res; ++r) s->res_seg[r] = -1;
    const md_protein_backbone_data_t* bb = &s->sys->protein_backbone;
    for (size_t i = 0; bb->segment.res_idx && i < bb->segment.count; ++i) {
        const int32_t r = bb->segment.res_idx[i];
        if (r >= 0 && (uint32_t)r < s->num_res) s->res_seg[r] = (int32_t)i;
    }
}

int32_t vial_sys_res_segment(vial_sys_t* s, uint32_t r) {
    need_res_segs(s);
    return r < s->num_res ? s->res_seg[r] : -1;
}

int32_t vial_sys_atom_segment(vial_sys_t* s, uint32_t a) {
    const int32_t r = vial_sys_atom_res(s, a);
    return r >= 0 ? vial_sys_res_segment(s, (uint32_t)r) : -1;
}

int32_t vial_sys_res_chain(const vial_sys_t* s, uint32_t r) {
    const md_system_t* sys = s->sys;
    if (!s->num_chains || !sys->chain.res_chain || r >= s->num_res) return -1;
    const int32_t c = sys->chain.res_chain[r];
    return c >= 0 && (uint32_t)c < s->num_chains ? c : -1;
}

md_str_t vial_sys_chain_id(const vial_sys_t* s, int32_t chain) {
    return chain >= 0 && (uint32_t)chain < s->num_chains ? md_system_chain_id(s->sys, (size_t)chain) : MD_STR_LIT("");
}

static uint32_t atom_res(const vial_sys_t* s, uint32_t a) {
    const md_residue_idx_t r = md_system_residue_find_by_atom_idx(s->sys, a);
    return r >= 0 ? (uint32_t)r : NONE32;
}

bool vial_sys_init(vial_sys_t* s, const md_system_t* sys, md_allocator_i* alloc) {
    memset(s, 0, sizeof(*s));
    if (!sys) return false;
    s->sys = sys;
    s->alloc = alloc;
    s->num_atoms  = (uint32_t)md_system_atom_count(sys);
    s->num_res    = sys->residue.atom_offset ? (uint32_t)md_system_residue_count(sys) : 0;
    s->num_parts  = sys->part.res_offset ? (uint32_t)md_system_part_count(sys) : 0;
    s->num_mols   = sys->molecule.offset ? (uint32_t)md_system_molecule_count(sys) : 0;
    s->num_chains = (uint32_t)md_system_chain_count(sys);
    s->num_types  = sys->atom.type_idx ? (uint32_t)sys->atom.type.count : 0;
    s->part_kind = md_alloc(alloc, s->num_parts + 1);
    s->part_chain = md_alloc(alloc, (s->num_parts + 1) * sizeof(int32_t));
    for (uint32_t p = 0; p < s->num_parts; ++p) {
        const md_urange_t res = md_system_part_residue_range(sys, p);
        s->part_kind[p] = part_kind(sys, p);
        s->part_chain[p] = res.beg < res.end ? vial_sys_res_chain(s, res.beg) : -1;
    }
    return true;
}

static void need_mols(vial_sys_t* s) {
    if (s->have_mols) return;
    s->have_mols = true;
    md_allocator_i* alloc = s->alloc;
    s->mol_offset = md_alloc(alloc, (s->num_mols + 1) * sizeof(uint32_t));
    s->mol_chain = md_alloc(alloc, (s->num_mols + 1) * sizeof(int32_t));
    const md_molecule_data_t* m = &s->sys->molecule;
    const uint32_t slots = s->num_mols ? m->offset[s->num_mols] : 0;
    s->mol_atoms = md_alloc(alloc, (slots + 1) * sizeof(uint32_t));
    uint32_t w = 0;
    for (uint32_t k = 0; k < s->num_mols; ++k) {
        s->mol_offset[k] = w;
        for (uint32_t q = m->offset[k]; q < m->offset[k + 1]; ++q) {
            const int32_t a = m->atom_idx[q];
            if (a >= 0 && (uint32_t)a < s->num_atoms) s->mol_atoms[w++] = (uint32_t)a;
        }
        qsort(s->mol_atoms + s->mol_offset[k], w - s->mol_offset[k], sizeof(uint32_t), cmp_u32);
        const uint32_t first = w > s->mol_offset[k] ? s->mol_atoms[s->mol_offset[k]] : NONE32;
        s->mol_chain[k] = first != NONE32 ? vial_sys_res_chain(s, atom_res(s, first)) : -1;
    }
    s->mol_offset[s->num_mols] = w;
}

static void need_tails(vial_sys_t* s) {
    if (s->have_tails) return;
    s->have_tails = true;
    const md_lipid_data_t* lip = &s->sys->lipid;
    s->num_tails = (uint32_t)lip->tail.count;
    s->tail_offset = md_alloc(s->alloc, (s->num_tails + 1) * sizeof(uint32_t));
    uint32_t total = 0;
    s->tail_offset[0] = 0;
    for (uint32_t t = 0; t < s->num_tails; ++t) s->tail_offset[t + 1] = total += (uint32_t)md_lipid_tail_atom_count(lip, t);
    s->tail_atoms = md_alloc(s->alloc, (total + 1) * sizeof(uint32_t));
    for (uint32_t t = 0; t < s->num_tails; ++t) {
        const md_atom_idx_t* a = md_lipid_tail_atoms(lip, t);
        uint32_t* dst = s->tail_atoms + s->tail_offset[t];
        const uint32_t n = s->tail_offset[t + 1] - s->tail_offset[t];
        for (uint32_t k = 0; k < n; ++k) dst[k] = (uint32_t)a[k];
        qsort(dst, n, sizeof(uint32_t), cmp_u32);
    }
}

// A head is its heavy atoms (md_lipid_data_t keeps those) and what else of its lipid is no tail and bonded to them, the
// hydrogens; in a lipid of one head, all of it that is no tail. As v1's head() takes it, hydrogens included.
static void need_heads(vial_sys_t* s) {
    if (s->have_heads) return;
    s->have_heads = true;
    const md_system_t* sys = s->sys;
    const md_lipid_data_t* lip = &sys->lipid;
    s->num_heads = (uint32_t)lip->head.count;
    int32_t* head_of = md_alloc(s->alloc, (s->num_atoms + 1) * sizeof(int32_t));
    for (uint32_t a = 0; a < s->num_atoms; ++a) head_of[a] = -1;
    for (uint32_t h = 0; h < s->num_heads; ++h) {
        const md_atom_idx_t* a = md_lipid_head_atoms(lip, h);
        const size_t n = md_lipid_head_atom_count(lip, h);
        for (size_t k = 0; k < n; ++k) if ((uint32_t)a[k] < s->num_atoms) head_of[a[k]] = (int32_t)h;
    }
    uint8_t* in_tail = md_alloc(s->alloc, s->num_atoms + 1);
    md_lipid_tail_atom_mask(in_tail, lip, s->num_atoms);
    for (size_t l = 0; l < lip->count; ++l) {
        const uint32_t h0 = lip->head_offset[l], h1 = lip->head_offset[l + 1];
        if (h0 == h1) continue;
        const int32_t* res = md_lipid_residues(lip, l);
        const size_t nr = md_lipid_residue_count(lip, l);
        for (size_t c = 0; c < nr; ++c) {
            const md_urange_t r = md_system_residue_atom_range(sys, (size_t)res[c]);
            for (uint32_t a = r.beg; a < r.end && a < s->num_atoms; ++a) {
                if (head_of[a] >= 0 || in_tail[a]) continue;
                if (h1 - h0 == 1) { head_of[a] = (int32_t)h0; continue; }
                for (md_bond_iter_t it = md_bond_iter(&sys->bond, a); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
                    const md_atom_idx_t o = md_bond_iter_atom_index(&it);
                    if (o >= 0 && (uint32_t)o < s->num_atoms && head_of[o] >= (int32_t)h0 && head_of[o] < (int32_t)h1) { head_of[a] = head_of[o]; break; }
                }
            }
        }
    }
    s->head_offset = md_alloc(s->alloc, (s->num_heads + 2) * sizeof(uint32_t));
    memset(s->head_offset, 0, (s->num_heads + 2) * sizeof(uint32_t));
    for (uint32_t a = 0; a < s->num_atoms; ++a) if (head_of[a] >= 0) s->head_offset[head_of[a] + 1] += 1;
    for (uint32_t h = 0; h < s->num_heads; ++h) s->head_offset[h + 1] += s->head_offset[h];
    s->head_atoms = md_alloc(s->alloc, (s->head_offset[s->num_heads] + 1) * sizeof(uint32_t));
    uint32_t* fill = md_alloc(s->alloc, (s->num_heads + 1) * sizeof(uint32_t));
    memcpy(fill, s->head_offset, (s->num_heads + 1) * sizeof(uint32_t));
    for (uint32_t a = 0; a < s->num_atoms; ++a) if (head_of[a] >= 0) s->head_atoms[fill[head_of[a]]++] = a;   // ascending
}

static void need_rings(vial_sys_t* s) {
    if (s->have_rings) return;
    s->have_rings = true;
    const md_ring_data_t* rd = &s->sys->ring;
    s->num_rings = (uint32_t)md_ring_data_count(rd);
    s->ring_offset = md_alloc(s->alloc, (s->num_rings + 1) * sizeof(uint32_t));
    uint32_t total = 0;
    s->ring_offset[0] = 0;
    for (uint32_t r = 0; r < s->num_rings; ++r) s->ring_offset[r + 1] = total += (uint32_t)md_ring_size(rd, r);
    s->ring_atoms = md_alloc(s->alloc, (total + 1) * sizeof(uint32_t));
    for (uint32_t r = 0; r < s->num_rings; ++r) {
        const int32_t* a = md_ring_atoms(rd, r);
        uint32_t* dst = s->ring_atoms + s->ring_offset[r];
        const uint32_t n = s->ring_offset[r + 1] - s->ring_offset[r];
        for (uint32_t k = 0; k < n; ++k) dst[k] = (uint32_t)a[k];
        qsort(dst, n, sizeof(uint32_t), cmp_u32);
    }
}

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return x < y ? -1 : x > y;
}

static void need_frags(vial_sys_t* s) {
    if (s->have_frags) return;
    s->have_frags = true;
    const md_system_t* sys = s->sys;
    const uint32_t num = (uint32_t)md_system_fragment_count(sys);
    // mdlib's fragments run by molecule and then as its template has them: put them in the order of their first
    // particle, each fragment's least, so fragment('phenyl')[0] is the first phenyl along the chain
    uint64_t* order = md_alloc(s->alloc, (num + 1) * sizeof(uint64_t));
    for (uint32_t f = 0; f < num; ++f) {
        uint32_t least = UINT32_MAX;
        for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, f); md_fragment_atom_iter_next(&it);) {
            if ((uint32_t)it.atom < least) least = (uint32_t)it.atom;
        }
        order[f] = (uint64_t)least << 32 | f;
    }
    qsort(order, num, sizeof(uint64_t), cmp_u64);
    s->num_frags = num;
    s->frag_idx = md_alloc(s->alloc, (num + 1) * sizeof(uint32_t));
    s->frag_offset = md_alloc(s->alloc, (num + 1) * sizeof(uint32_t));
    uint32_t total = 0;
    s->frag_offset[0] = 0;
    for (uint32_t u = 0; u < num; ++u) {
        s->frag_idx[u] = (uint32_t)order[u];
        s->frag_offset[u + 1] = total += (uint32_t)md_system_fragment_atom_count(sys, s->frag_idx[u]);
    }
    s->frag_atoms = md_alloc(s->alloc, (total + 1) * sizeof(uint32_t));
    for (uint32_t u = 0; u < num; ++u) {
        uint32_t* dst = s->frag_atoms + s->frag_offset[u];
        uint32_t n = 0;
        const uint32_t room = s->frag_offset[u + 1] - s->frag_offset[u];
        for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, s->frag_idx[u]); md_fragment_atom_iter_next(&it) && n < room;) {
            dst[n++] = (uint32_t)it.atom;
        }
        qsort(dst, n, sizeof(uint32_t), cmp_u32);
    }
}

static void need_groups(vial_sys_t* s) {
    if (s->have_groups) return;
    s->have_groups = true;
    const md_group_data_t* g = &s->sys->group;
    s->num_groups = (uint32_t)g->count;
    s->group_offset = md_alloc(s->alloc, (s->num_groups + 1) * sizeof(uint32_t));
    uint32_t total = 0;
    s->group_offset[0] = 0;
    for (uint32_t i = 0; i < s->num_groups; ++i) s->group_offset[i + 1] = total += (uint32_t)md_bitset_count_range(g->atoms[i], 0, s->num_atoms);
    s->group_atoms = md_alloc(s->alloc, (total + 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < s->num_groups; ++i) {
        md_bitset_iter_t it = md_bitset_iter_range(g->atoms[i], 0, s->num_atoms);
        md_bitset_iter_extract_indices((int32_t*)(s->group_atoms + s->group_offset[i]), s->group_offset[i + 1] - s->group_offset[i], it);
    }
}

static void need_chains(vial_sys_t* s) {
    if (s->have_chains) return;
    s->have_chains = true;
    md_allocator_i* alloc = s->alloc;
    const md_system_t* sys = s->sys;
    s->chain_offset = md_alloc(alloc, (s->num_chains + 2) * sizeof(uint32_t));
    memset(s->chain_offset, 0, (s->num_chains + 2) * sizeof(uint32_t));
    for (uint32_t r = 0; r < s->num_res; ++r) {
        const int32_t c = vial_sys_res_chain(s, r);
        if (c >= 0) s->chain_offset[c + 1] += (uint32_t)md_system_residue_atom_count(sys, r);
    }
    for (uint32_t c = 0; c < s->num_chains; ++c) s->chain_offset[c + 1] += s->chain_offset[c];
    s->chain_atoms = md_alloc(alloc, (s->chain_offset[s->num_chains] + 1) * sizeof(uint32_t));
    uint32_t* fill = md_alloc(alloc, (s->num_chains + 1) * sizeof(uint32_t));
    memcpy(fill, s->chain_offset, (s->num_chains + 1) * sizeof(uint32_t));
    for (uint32_t r = 0; r < s->num_res; ++r) {
        const int32_t c = vial_sys_res_chain(s, r);
        if (c < 0) continue;
        const md_urange_t a = md_system_residue_atom_range(sys, r);
        for (uint32_t i = a.beg; i < a.end; ++i) s->chain_atoms[fill[c]++] = i;     // residues ascend, so atoms do
    }
}

// ### SETS ###

static vial_sets_t* sets_new(uint32_t count, uint32_t total, md_allocator_i* alloc) {
    vial_sets_t* s = md_alloc(alloc, sizeof(vial_sets_t));
    memset(s, 0, sizeof(*s));
    s->count = count;
    s->offset = md_alloc(alloc, (count + 1) * sizeof(uint32_t));
    s->idx = md_alloc(alloc, (total + 1) * sizeof(uint32_t));
    s->offset[0] = 0;
    return s;
}

vial_sets_t* vial_sets_one(const uint32_t* idx, uint32_t n, md_allocator_i* alloc) {
    vial_sets_t* s = sets_new(1, n, alloc);
    s->single = true;
    if (n) memcpy(s->idx, idx, n * sizeof(uint32_t));
    s->offset[1] = n;
    return s;
}

vial_sets_t* vial_sets_particles(const int32_t* atoms, uint32_t n, md_allocator_i* alloc) {
    vial_sets_t* s = sets_new(n, n, alloc);
    s->level = VIAL_LEVEL_PARTICLE;
    s->unit = md_alloc(alloc, (n ? n : 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < n; ++i) {
        s->idx[i] = (uint32_t)atoms[i];
        s->unit[i] = (uint32_t)atoms[i];
        s->offset[i + 1] = i + 1;
    }
    return s;
}

vial_sets_t* vial_sets_all(uint32_t num_atoms, md_allocator_i* alloc) {
    vial_sets_t* s = sets_new(1, num_atoms, alloc);
    s->single = true;
    for (uint32_t i = 0; i < num_atoms; ++i) s->idx[i] = i;
    s->offset[1] = num_atoms;
    return s;
}

vial_sets_t* vial_sets_and(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc) {
    const uint32_t* x = vial_sets_items(a, ia); const uint32_t nx = vial_sets_size(a, ia);
    const uint32_t* y = vial_sets_items(b, ib); const uint32_t ny = vial_sets_size(b, ib);
    vial_sets_t* s = sets_new(1, nx < ny ? nx : ny, alloc);
    s->single = true;
    uint32_t i = 0, j = 0, n = 0;
    while (i < nx && j < ny) {
        if (x[i] < y[j]) ++i;
        else if (y[j] < x[i]) ++j;
        else { s->idx[n++] = x[i]; ++i; ++j; }
    }
    s->offset[1] = n;
    return s;
}

vial_sets_t* vial_sets_units_at(const vial_sets_t* a, const vial_sets_t* u, const int64_t* pos, const bool* skip, uint32_t m,
                                uint32_t* outside, md_allocator_i* alloc) {
    *outside = 0;
    // Each particle's units, as offsets into one list: units may share particles, as fused rings do
    uint32_t natoms = 0;
    for (uint32_t k = 0; k < u->offset[u->count]; ++k) natoms = u->idx[k] + 1 > natoms ? u->idx[k] + 1 : natoms;
    uint32_t* off = md_alloc(alloc, (natoms + 1) * sizeof(uint32_t));
    memset(off, 0, (natoms + 1) * sizeof(uint32_t));
    for (uint32_t k = 0; k < u->offset[u->count]; ++k) off[u->idx[k] + 1] += 1;
    for (uint32_t p = 0; p < natoms; ++p) off[p + 1] += off[p];
    const uint32_t total_refs = u->offset[u->count];
    uint32_t* of = md_alloc(alloc, (total_refs ? total_refs : 1) * sizeof(uint32_t));
    uint32_t* fill = md_alloc(alloc, (natoms ? natoms : 1) * sizeof(uint32_t));
    if (natoms) memcpy(fill, off, natoms * sizeof(uint32_t));
    for (uint32_t j = 0; j < u->count; ++j) {
        for (uint32_t k = u->offset[j]; k < u->offset[j + 1]; ++k) of[fill[u->idx[k]]++] = j;
    }
    uint32_t* stamp = md_alloc(alloc, (u->count ? u->count : 1) * sizeof(uint32_t));
    for (uint32_t j = 0; j < u->count; ++j) stamp[j] = UINT32_MAX;
    uint32_t* touched = md_alloc(alloc, (u->count ? u->count : 1) * sizeof(uint32_t));

    const uint32_t count = a->count * m;
    vial_sets_t* s = md_alloc(alloc, sizeof(vial_sets_t));
    memset(s, 0, sizeof(*s));
    s->count = count;
    s->offset = md_alloc(alloc, (count + 1) * sizeof(uint32_t));
    s->offset[0] = 0;
    uint32_t* units = md_alloc(alloc, (count ? count : 1) * sizeof(uint32_t));
    md_array(uint32_t) idx = NULL;
    bool whole = u->level != VIAL_LEVEL_NONE && u->unit != NULL;
    uint32_t c = 0;
    for (uint32_t i = 0; i < a->count; ++i) {
        const uint32_t* x = vial_sets_items(a, i);
        const uint32_t nx = vial_sets_size(a, i);
        uint32_t n = 0;
        for (uint32_t k = 0; k < nx; ++k) {
            if (x[k] >= natoms) continue;
            for (uint32_t r = off[x[k]]; r < off[x[k] + 1]; ++r) {
                const uint32_t j = of[r];
                if (stamp[j] != i) { stamp[j] = i; touched[n++] = j; }
            }
        }
        qsort(touched, n, sizeof(uint32_t), cmp_u32);
        bool short_set = false;
        for (uint32_t e = 0; e < m; ++e) {
            const int64_t q = pos[e] < 0 ? pos[e] + (int64_t)n : pos[e];
            const bool none = (skip && skip[e]) || q < 0 || q >= (int64_t)n;
            if (none) {
                if (!(skip && skip[e])) short_set = true;
                whole = false;
                units[c] = 0;
            } else {
                const uint32_t j = touched[q];
                const uint32_t* y = vial_sets_items(u, j);
                const uint32_t ny = vial_sets_size(u, j);
                uint32_t p = 0, t = 0, shared = 0;
                while (p < nx && t < ny) {
                    if (x[p] < y[t]) ++p;
                    else if (y[t] < x[p]) ++t;
                    else { md_array_push(idx, x[p], alloc); ++shared; ++p; ++t; }
                }
                if (shared != ny) whole = false;
                units[c] = u->unit ? u->unit[j] : 0;
            }
            s->offset[++c] = (uint32_t)md_array_size(idx);
        }
        *outside += short_set;
    }
    s->idx = idx ? idx : md_alloc(alloc, sizeof(uint32_t));
    if (whole && count) { s->level = u->level; s->unit = units; }
    return s;
}

vial_sets_t* vial_sets_or(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc) {
    const uint32_t* x = vial_sets_items(a, ia); const uint32_t nx = vial_sets_size(a, ia);
    const uint32_t* y = vial_sets_items(b, ib); const uint32_t ny = vial_sets_size(b, ib);
    vial_sets_t* s = sets_new(1, nx + ny, alloc);
    s->single = true;
    uint32_t i = 0, j = 0, n = 0;
    while (i < nx || j < ny) {
        if (j >= ny || (i < nx && x[i] < y[j])) s->idx[n++] = x[i++];
        else if (i >= nx || y[j] < x[i]) s->idx[n++] = y[j++];
        else { s->idx[n++] = x[i]; ++i; ++j; }
    }
    s->offset[1] = n;
    return s;
}

vial_sets_t* vial_sets_not(const vial_sets_t* scope, const vial_sets_t* a, uint32_t ia, md_allocator_i* alloc) {
    const uint32_t* x = vial_sets_items(scope, 0); const uint32_t nx = vial_sets_size(scope, 0);
    const uint32_t* y = vial_sets_items(a, ia); const uint32_t ny = vial_sets_size(a, ia);
    vial_sets_t* s = sets_new(1, nx, alloc);
    s->single = true;
    uint32_t i = 0, j = 0, n = 0;
    while (i < nx) {
        while (j < ny && y[j] < x[i]) ++j;
        if (j < ny && y[j] == x[i]) { ++i; continue; }
        s->idx[n++] = x[i++];
    }
    s->offset[1] = n;
    return s;
}

vial_sets_t* vial_sets_boolop_each(const vial_sets_t* arr, const vial_sets_t* one, bool is_and, md_allocator_i* alloc) {
    // Each element on its own, then gathered, as the runtime's set algebra gathers them: no level, as an element is no
    // longer a whole unit
    vial_sets_t** parts = md_alloc(alloc, (arr->count ? arr->count : 1) * sizeof(vial_sets_t*));
    uint32_t total = 0;
    for (uint32_t i = 0; i < arr->count; ++i) {
        parts[i] = is_and ? vial_sets_and(arr, i, one, 0, alloc) : vial_sets_or(arr, i, one, 0, alloc);
        total += parts[i]->offset[1];
    }
    vial_sets_t* s = sets_new(arr->count, total, alloc);
    uint32_t at = 0;
    for (uint32_t i = 0; i < arr->count; ++i) {
        md_memcpy(s->idx + at, parts[i]->idx, parts[i]->offset[1] * sizeof(uint32_t));
        at += parts[i]->offset[1];
        s->offset[i + 1] = at;
    }
    return s;
}

vial_sets_t* vial_sets_concat(const vial_sets_t* const* parts, uint32_t n, md_allocator_i* alloc) {
    uint32_t count = 0, total = 0;
    for (uint32_t i = 0; i < n; ++i) { count += parts[i]->count; total += parts[i]->offset[parts[i]->count]; }
    vial_sets_t* s = sets_new(count, total, alloc);
    uint32_t c = 0, at = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const vial_sets_t* p = parts[i];
        for (uint32_t j = 0; j < p->count; ++j) {
            const uint32_t sz = vial_sets_size(p, j);
            memcpy(s->idx + at, vial_sets_items(p, j), sz * sizeof(uint32_t));
            at += sz;
            s->offset[++c] = at;
        }
    }
    return s;
}

vial_sets_t* vial_sets_merge(const vial_sets_t* a, md_allocator_i* alloc) {
    const uint32_t total = a->offset[a->count];
    vial_sets_t* s = sets_new(1, total, alloc);
    s->single = true;
    if (total) memcpy(s->idx, a->idx, total * sizeof(uint32_t));
    qsort(s->idx, total, sizeof(uint32_t), cmp_u32);
    uint32_t n = 0;
    for (uint32_t i = 0; i < total; ++i) if (n == 0 || s->idx[n - 1] != s->idx[i]) s->idx[n++] = s->idx[i];
    s->offset[1] = n;
    return s;
}

vial_sets_t* vial_sets_pick(const vial_sets_t* a, const uint32_t* pos, uint32_t n, md_allocator_i* alloc) {
    uint32_t total = 0;
    for (uint32_t k = 0; k < n; ++k) total += vial_sets_size(a, pos[k]);
    vial_sets_t* s = sets_new(n, total, alloc);
    s->level = a->level;
    s->root = vial_sets_root(a);
    if (a->unit) s->unit = md_alloc(alloc, (n + 1) * sizeof(uint32_t));
    uint32_t w = 0;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t sz = vial_sets_size(a, pos[k]);
        memcpy(s->idx + w, vial_sets_items(a, pos[k]), sz * sizeof(uint32_t));
        w += sz;
        s->offset[k + 1] = w;
        if (s->unit) s->unit[k] = a->unit[pos[k]];
    }
    return s;
}

vial_sets_t* vial_sets_pick_particles(const vial_sets_t* a, uint32_t i, const uint32_t* pos, uint32_t n, md_allocator_i* alloc) {
    vial_sets_t* s = sets_new(1, n, alloc);
    s->single = true;
    if (a->single) s->root = vial_sets_root(a);
    const uint32_t* items = vial_sets_items(a, i);
    for (uint32_t k = 0; k < n; ++k) s->idx[k] = items[pos[k]];
    qsort(s->idx, n, sizeof(uint32_t), cmp_u32);
    uint32_t w = 0;
    for (uint32_t k = 0; k < n; ++k) if (w == 0 || s->idx[w - 1] != s->idx[k]) s->idx[w++] = s->idx[k];
    s->offset[1] = w;
    return s;
}

// ### LABELS ###

static size_t residue_label(char* buf, size_t cap, const vial_sys_t* s, uint32_t r) {
    const md_str_t name = md_system_residue_name(s->sys, r);
    const int32_t c = vial_sys_res_chain(s, r);
    if (c >= 0) {
        const md_str_t ch = md_system_chain_id(s->sys, (size_t)c);
        return (size_t)snprintf(buf, cap, "%.*s%d/%.*s", (int)name.len, name.ptr, (int)md_system_residue_seq_id(s->sys, r), (int)ch.len, ch.ptr);
    }
    return (size_t)snprintf(buf, cap, "%.*s%d", (int)name.len, name.ptr, (int)md_system_residue_seq_id(s->sys, r));
}

// The word of a part's kind (VIAL_PART_*), as part(kind=) takes it: the kind proper, not polymer for a peptide
const char* vial_part_kind_word(uint8_t kind) {
    static const char* words[] = {"ligand", "peptide", "dna", "rna", "nucleic", "polymer", "glycan", "lipid", "ion", "water"};
    return kind < MD_COUNTOF(words) ? words[kind] : "";
}

size_t vial_unit_label(char* buf, size_t cap, vial_sys_t* s, uint8_t level, uint32_t unit) {
    if (cap) buf[0] = '\0';
    const md_system_t* sys = s->sys;
    int n = 0;
    switch (level) {
    case VIAL_LEVEL_PARTICLE: {
        const bool has = sys->atom.type_idx && sys->atom.type_idx[unit] < sys->atom.type.count;
        const md_str_t name = has ? MD_LABEL_STR(sys->atom.type.name[sys->atom.type_idx[unit]]) : MD_STR_LIT("");
        n = snprintf(buf, cap, "%.*s/", (int)name.len, name.ptr);
        const uint32_t r = atom_res(s, unit);
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        break;
    }
    case VIAL_LEVEL_RESIDUE: return residue_label(buf, cap, s, unit);
    case VIAL_LEVEL_RING: {
        // Its number among the rings of the residue it starts in: ring2/TRP5
        need_rings(s);
        if (unit >= s->num_rings) break;
        const uint32_t first = s->ring_offset[unit] < s->ring_offset[unit + 1] ? s->ring_atoms[s->ring_offset[unit]] : NONE32;
        const uint32_t r = first != NONE32 ? atom_res(s, first) : NONE32;
        uint32_t k = 0;
        for (uint32_t q = unit; q > 0 && r != NONE32; --q) {
            const uint32_t f = s->ring_offset[q - 1] < s->ring_offset[q] ? s->ring_atoms[s->ring_offset[q - 1]] : NONE32;
            if (f == NONE32 || atom_res(s, f) != r) break;
            ++k;
        }
        n = snprintf(buf, cap, "ring%u/", k + 1);
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        break;
    }
    case VIAL_LEVEL_FRAGMENT: {
        // Its name and the residue its first particle is in: phenyl/TYR66, peptide/ALA5
        need_frags(s);
        if (unit >= s->num_frags) break;
        const md_str_t name = md_system_fragment_name(sys, s->frag_idx[unit]);
        n = snprintf(buf, cap, "%.*s/", (int)name.len, name.ptr);
        const uint32_t first = s->frag_offset[unit] < s->frag_offset[unit + 1] ? s->frag_atoms[s->frag_offset[unit]] : NONE32;
        const uint32_t r = first != NONE32 ? atom_res(s, first) : NONE32;
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        else if (n > 0 && (size_t)n <= cap) buf[--n] = '\0';    // in no residue: its name alone
        break;
    }
    case VIAL_LEVEL_GROUP: {
        // Its name; a group that shares its name with another, its roles too: rest/tc
        const md_group_data_t* g = &sys->group;
        if (unit >= g->count) break;
        const md_str_t name = g->name[unit];
        n = snprintf(buf, cap, "%.*s", (int)name.len, name.ptr);
        bool shared = false;
        for (size_t i = 0; i < g->count && !shared; ++i) shared = i != unit && md_str_eq(g->name[i], name);
        for (int r = 0; shared && r < MD_GROUP_ROLE_COUNT; ++r) {
            if (!(g->roles[unit] & (1u << r))) continue;
            if (n > 0 && (size_t)n < cap) n += snprintf(buf + n, cap - n, "/%s", md_group_role_name((md_group_role_t)r));
        }
        break;
    }
    case VIAL_LEVEL_HEAD: {
        // Its number within its lipid and the residue it starts in: head1/DPPC1
        need_heads(s);
        if (unit >= s->num_heads) break;
        const md_lipid_data_t* lip = &sys->lipid;
        uint32_t k = 0;
        for (size_t l = 0; l < lip->count; ++l) {
            if (unit >= lip->head_offset[l] && unit < lip->head_offset[l + 1]) { k = unit - lip->head_offset[l]; break; }
        }
        n = snprintf(buf, cap, "head%u/", k + 1);
        const uint32_t first = s->head_offset[unit] < s->head_offset[unit + 1] ? s->head_atoms[s->head_offset[unit]] : NONE32;
        const uint32_t r = first != NONE32 ? atom_res(s, first) : NONE32;
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        break;
    }
    case VIAL_LEVEL_TAIL: {
        // Its number within its lipid and the residue it starts in: tail2/DPPC1
        need_tails(s);
        if (unit >= s->num_tails) break;
        const md_lipid_data_t* lip = &sys->lipid;
        uint32_t k = 0;
        for (size_t l = 0; l < lip->count; ++l) {
            if (unit >= lip->tail_offset[l] && unit < lip->tail_offset[l + 1]) { k = unit - lip->tail_offset[l]; break; }
        }
        n = snprintf(buf, cap, "tail%u/", k + 1);
        const uint32_t first = s->tail_offset[unit] < s->tail_offset[unit + 1] ? s->tail_atoms[s->tail_offset[unit]] : NONE32;
        const uint32_t r = first != NONE32 ? atom_res(s, first) : NONE32;
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        break;
    }
    case VIAL_LEVEL_CHAIN: {
        const md_str_t id = md_system_chain_id(sys, unit);
        n = snprintf(buf, cap, "%.*s", (int)id.len, id.ptr);
        break;
    }
    case VIAL_LEVEL_PART:
    case VIAL_LEVEL_MOLECULE: {
        // One residue: that residue; else what it is, and its chain
        uint32_t r0 = NONE32, nres = 0;
        int32_t chain = -1;
        md_str_t what = MD_STR_LIT("");
        if (level == VIAL_LEVEL_PART) {
            const md_urange_t res = md_system_part_residue_range(sys, unit);
            r0 = res.beg; nres = res.end - res.beg;
            chain = s->part_chain[unit];
            what = md_str_from_cstr(vial_part_kind_word(s->part_kind[unit]));
        } else {
            need_mols(s);
            const uint32_t b = s->mol_offset[unit], e = s->mol_offset[unit + 1];
            if (b < e) {
                r0 = atom_res(s, s->mol_atoms[b]);
                nres = r0 != NONE32 && atom_res(s, s->mol_atoms[e - 1]) == r0 ? 1 : 2;
            }
            chain = s->mol_chain[unit];
            const md_molecule_type_idx_t t = md_system_molecule_type_idx(sys, unit);
            if (t >= 0) what = md_system_molecule_type_name(sys, (size_t)t);
        }
        if (nres == 1 && r0 != NONE32) return residue_label(buf, cap, s, r0);
        if (chain >= 0) {
            const md_str_t ch = md_system_chain_id(sys, (size_t)chain);
            n = snprintf(buf, cap, "%.*s/%.*s", (int)what.len, what.ptr, (int)ch.len, ch.ptr);
        } else {
            n = snprintf(buf, cap, "%.*s%u", (int)what.len, what.ptr, unit + 1);
        }
        break;
    }
    default: break;
    }
    return n > 0 ? ((size_t)n < cap ? (size_t)n : cap ? cap - 1 : 0) : 0;
}

// ### THE RUN'S FRAMES ###

double* vial_run_times(const md_system_t* sys, md_str_t run, md_allocator_i* alloc, size_t* num_frames) {
    *num_frames = 0;
    if (!sys || md_str_empty(run)) return NULL;
    char path[512];
    snprintf(path, sizeof(path), "%.*s/time", (int)run.len, run.ptr);
    const md_attribute_t* time = md_attributes_find(&sys->attributes, md_str_from_cstr(path));
    if (!time || time->format.rank < 1 || !time->format.shape[0]) return NULL;
    const size_t n = time->format.shape[0];
    double* t = md_alloc(alloc, n * sizeof(double));
    double f;
    const md_unit_t unit = !md_unit_is_none(time->unit) && md_unit_conversion_factor(&f, time->unit, md_unit_picosecond()) ? md_unit_picosecond() : md_unit_none();
    if (md_attribute_extract_f64(t, n, &sys->attributes, time, md_attribute_slice_all(), unit) != n) return NULL;
    *num_frames = n;
    return t;
}

int vial_run_frame_at(const md_system_t* sys, md_str_t run, double t, int64_t* frame, double* frame_time) {
    md_temp_scope_t temp = md_temp_begin();
    size_t n = 0;
    const double* times = vial_run_times(sys, run, md_temp_allocator(temp), &n);
    int r = 0;
    if (times) {
        // The axis does not decrease (md_attributes): the first frame not before t, and the one before it
        size_t lo = 0, hi = n;
        while (lo < hi) { const size_t mid = lo + (hi - lo) / 2; if (times[mid] < t) lo = mid + 1; else hi = mid; }
        size_t k = lo < n ? lo : n - 1;
        if (lo > 0 && (lo == n || t - times[lo - 1] <= times[lo] - t)) k = lo - 1;
        const double half = n > 1 ? 0.5 * (times[n - 1] - times[0]) / (double)(n - 1) : 0.0;
        const double slack = 1e-9 * (fabs(times[k]) > 1.0 ? fabs(times[k]) : 1.0);
        if (t < times[0] - half - slack || t > times[n - 1] + half + slack) {
            r = -1;
        } else {
            *frame = (int64_t)k;
            *frame_time = times[k];
            r = 1;
        }
    }
    md_temp_end(temp);
    return r;
}

// ### SELECTORS ###

static const struct { const char* word; uint8_t key; } keywords[] = {
    {"name", VIAL_KEY_NAME}, {"serial", VIAL_KEY_SERIAL}, {"element", VIAL_KEY_ELEMENT}, {"kind", VIAL_KEY_KIND},
    {"type", VIAL_KEY_TYPE}, {"resid", VIAL_KEY_RESID}, {"resname", VIAL_KEY_RESNAME}, {"chain", VIAL_KEY_CHAIN}, {"id", VIAL_KEY_ID},
    {"scope", VIAL_KEY_SCOPE}, {"size", VIAL_KEY_SIZE}, {"role", VIAL_KEY_ROLE},
    {"hybridization", VIAL_KEY_HYBRIDIZATION}, {"aromatic", VIAL_KEY_AROMATIC}, {"metallic", VIAL_KEY_METALLIC},
    {"hydrogens", VIAL_KEY_HYDROGENS}, {"formal_charge", VIAL_KEY_FORMAL_CHARGE}, {"order", VIAL_KEY_ORDER}, {"origin", VIAL_KEY_ORIGIN},
    {"gap", VIAL_KEY_GAP}, {"resolved", VIAL_KEY_RESOLVED}, {"ss", VIAL_KEY_SS}, {"pucker", VIAL_KEY_PUCKER},
};

uint8_t vial_select_key(md_str_t keyword) {
    for (size_t i = 0; i < MD_COUNTOF(keywords); ++i) if (md_str_eq(keyword, md_str_from_cstr(keywords[i].word))) return keywords[i].key;
    return VIAL_KEY_NONE;
}

md_str_t vial_match_encode(const vial_match_t* m, md_allocator_i* alloc) {
    switch (m->kind) {
    case VIAL_MATCH_STR:   return md_alloc_printf(alloc, "s%.*s", (int)m->s.len, m->s.ptr);
    case VIAL_MATCH_INT:   return md_alloc_printf(alloc, "i%lld", (long long)m->lo);
    case VIAL_MATCH_RANGE: return md_alloc_printf(alloc, "r%lld:%lld", (long long)m->lo, (long long)m->hi);
    case VIAL_MATCH_BOOL:  return md_alloc_printf(alloc, "b%d", m->lo ? 1 : 0);
    case VIAL_MATCH_AND:   return md_alloc_printf(alloc, "&%u", m->num);
    case VIAL_MATCH_OR:    return md_alloc_printf(alloc, "|%u", m->num);
    case VIAL_MATCH_NOT:   return md_alloc_printf(alloc, "!");
    default:               return (md_str_t){0};
    }
}

bool vial_match_decode(md_str_t s, vial_match_t* m) {
    *m = (vial_match_t){0};
    if (!s.len) return false;
    const md_str_t rest = {s.ptr + 1, s.len - 1};
    switch (s.ptr[0]) {
    case 's': m->kind = VIAL_MATCH_STR; m->s = rest; return true;
    case 'i': m->kind = VIAL_MATCH_INT; m->lo = md_parse_int(rest); return md_parse_is_int(rest);
    case 'b': m->kind = VIAL_MATCH_BOOL; m->lo = md_parse_int(rest); return md_parse_is_int(rest);
    case '&': m->kind = VIAL_MATCH_AND; m->num = (uint32_t)md_parse_int(rest); return md_parse_is_int(rest);
    case '|': m->kind = VIAL_MATCH_OR; m->num = (uint32_t)md_parse_int(rest); return md_parse_is_int(rest);
    case '!': m->kind = VIAL_MATCH_NOT; m->num = 1; return rest.len == 0;
    case 'r': {
        size_t c = 0;
        while (c < rest.len && rest.ptr[c] != ':') ++c;
        const md_str_t lo = {rest.ptr, c}, hi = {rest.ptr + MD_MIN(c + 1, rest.len), rest.len - MD_MIN(c + 1, rest.len)};
        m->kind = VIAL_MATCH_RANGE;
        m->lo = md_parse_int(lo);
        m->hi = md_parse_int(hi);
        return c < rest.len && md_parse_is_int(lo) && md_parse_is_int(hi);
    }
    default: return false;
    }
}

const char* vial_select_key_word(uint8_t key) {
    for (size_t i = 0; i < MD_COUNTOF(keywords); ++i) if (keywords[i].key == key) return keywords[i].word;
    return "";
}

#define K(k) (1u << (k))
static uint32_t takes(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return K(VIAL_KEY_NAME) | K(VIAL_KEY_SERIAL) | K(VIAL_KEY_ELEMENT) | K(VIAL_KEY_KIND) | K(VIAL_KEY_TYPE) |
                                     K(VIAL_KEY_RESID) | K(VIAL_KEY_RESNAME) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE) |
                                     K(VIAL_KEY_ROLE) | K(VIAL_KEY_HYBRIDIZATION) | K(VIAL_KEY_AROMATIC) | K(VIAL_KEY_METALLIC) |
                                     K(VIAL_KEY_HYDROGENS) | K(VIAL_KEY_FORMAL_CHARGE);
    case VIAL_LEVEL_RESIDUE:  return K(VIAL_KEY_ID) | K(VIAL_KEY_NAME) | K(VIAL_KEY_KIND) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_ROLE) |
                                     K(VIAL_KEY_GAP) | K(VIAL_KEY_RESOLVED) | K(VIAL_KEY_SS) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_PART:     return K(VIAL_KEY_NAME) | K(VIAL_KEY_KIND) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_MOLECULE: return K(VIAL_KEY_NAME) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_CHAIN:    return K(VIAL_KEY_ID) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_TAIL:     return K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_HEAD:     return K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_RING:     return K(VIAL_KEY_KIND) | K(VIAL_KEY_SIZE) | K(VIAL_KEY_PUCKER) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_GROUP:    return K(VIAL_KEY_NAME) | K(VIAL_KEY_ROLE) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_FRAGMENT: return K(VIAL_KEY_NAME) | K(VIAL_KEY_SCOPE);
    default: return 0;
    }
}
#undef K

bool vial_select_takes(uint8_t level, uint8_t key) { return key < 32 && (takes(level) >> key) & 1; }

const char* vial_select_keywords(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "name, serial, element, kind, type, resid, resname, chain, role, hybridization, aromatic, metallic, "
                                     "hydrogens, formal_charge and scope";
    case VIAL_LEVEL_RESIDUE:  return "id, name, kind, chain, role, gap, resolved, ss and scope";
    case VIAL_LEVEL_PART:     return "name, kind, chain and scope";
    case VIAL_LEVEL_MOLECULE: return "name, chain and scope";
    case VIAL_LEVEL_CHAIN:    return "id and scope";
    case VIAL_LEVEL_TAIL:     return "scope";
    case VIAL_LEVEL_HEAD:     return "scope";
    case VIAL_LEVEL_RING:     return "kind, size, pucker and scope";
    case VIAL_LEVEL_GROUP:    return "name, role and scope";
    case VIAL_LEVEL_FRAGMENT: return "name and scope";
    default: return "";
    }
}

const char* vial_select_name(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "particle";
    case VIAL_LEVEL_RESIDUE:  return "residue";
    case VIAL_LEVEL_PART:     return "part";
    case VIAL_LEVEL_MOLECULE: return "molecule";
    case VIAL_LEVEL_CHAIN:    return "chain";
    case VIAL_LEVEL_TAIL:     return "tail";
    case VIAL_LEVEL_HEAD:     return "head";
    case VIAL_LEVEL_RING:     return "ring";
    case VIAL_LEVEL_GROUP:    return "group";
    case VIAL_LEVEL_FRAGMENT: return "fragment";
    default: return "sets";
    }
}

static const char* residue_kinds[] = {"other", "amino_acid", "nucleotide", "water", "ion", "lipid"};   // md_residue_kind_t
static const char* particle_kinds[] = {"atom", "bead", "virtual"};                                    // md_particle_kind_t
static const char* ring_kinds[] = {"aromatic", "hetero", "fused", "fused_aromatic", "sugar"};
static const char* part_kinds[] = {"polymer", "peptide", "nucleic", "dna", "rna", "glycan", "ligand", "lipid", "ion", "water"};

static bool in_list(const char** list, size_t n, md_str_t s) {
    for (size_t i = 0; i < n; ++i) if (md_str_eq(s, md_str_from_cstr(list[i]))) return true;
    return false;
}

bool vial_select_kind_known(uint8_t level, md_str_t kind) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return in_list(particle_kinds, MD_COUNTOF(particle_kinds), kind);
    case VIAL_LEVEL_RESIDUE:  return in_list(residue_kinds, MD_COUNTOF(residue_kinds), kind);
    case VIAL_LEVEL_PART:     return in_list(part_kinds, MD_COUNTOF(part_kinds), kind);
    case VIAL_LEVEL_RING:     return in_list(ring_kinds, MD_COUNTOF(ring_kinds), kind);
    default: return false;
    }
}

const char* vial_select_kinds(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "atom, bead and virtual";
    case VIAL_LEVEL_RESIDUE:  return "amino_acid, nucleotide, water, ion, lipid and other";
    case VIAL_LEVEL_PART:     return "polymer, peptide, nucleic, dna, rna, glycan, ligand, lipid, ion and water";
    case VIAL_LEVEL_RING:     return "aromatic, hetero, fused, fused_aromatic and sugar";
    default: return "";
    }
}

// The class that starts at p.ptr[i], a '[', as Python's fnmatch has it: one character of the set, [abc], with ranges,
// [a-z], or of none of it, [!abc]; a ']' first is a member, as a '-' first or last is. Whether c is one of it in *hit,
// and its end, past the ']', in *end. False where no ']' closes it: the '[' is a character like any other.
static bool glob_class(md_str_t p, size_t i, char c, bool* hit, size_t* end) {
    size_t k = i + 1;
    const bool neg = k < p.len && p.ptr[k] == '!';
    if (neg) ++k;
    const size_t first = k;
    const unsigned char u = (unsigned char)c;
    bool in = false;
    while (k < p.len && (p.ptr[k] != ']' || k == first)) {
        if (k + 2 < p.len && p.ptr[k + 1] == '-' && p.ptr[k + 2] != ']') {
            in = in || ((unsigned char)p.ptr[k] <= u && u <= (unsigned char)p.ptr[k + 2]);
            k += 3;
        } else {
            in = in || p.ptr[k] == c;
            k += 1;
        }
    }
    if (k >= p.len) return false;
    *hit = in != neg;
    *end = k + 1;
    return true;
}

bool vial_glob(md_str_t p, md_str_t t) {
    // Iterative, with one backtrack point for the last *: linear in practice. A * is always a wildcard, also where the
    // text has a * there: [*] is a literal one
    size_t i = 0, j = 0, star = SIZE_MAX, mark = 0;
    while (j < t.len) {
        if (i < p.len && p.ptr[i] == '*') { star = i++; mark = j; continue; }
        if (i < p.len) {
            bool hit = false;
            size_t end = 0;
            if (p.ptr[i] == '[' && glob_class(p, i, t.ptr[j], &hit, &end)) {
                if (hit) { i = end; ++j; continue; }
            } else if (p.ptr[i] == '?' || p.ptr[i] == t.ptr[j]) {
                ++i; ++j;
                continue;
            }
        }
        if (star == SIZE_MAX) return false;
        i = star + 1;
        j = ++mark;
    }
    while (i < p.len && p.ptr[i] == '*') ++i;
    return i == p.len;
}

bool vial_is_glob(md_str_t p) {
    for (size_t i = 0; i < p.len; ++i) {
        bool hit;
        size_t end;
        if (p.ptr[i] == '*' || p.ptr[i] == '?' || (p.ptr[i] == '[' && glob_class(p, i, 0, &hit, &end))) return true;
    }
    return false;
}

// The node at *i and its operands, *i left past them: every operand is walked, also past the one that decides
static bool holds_at(const vial_match_t* e, uint32_t n, uint32_t* i, bool (*test)(const vial_match_t*, const void*), const void* user) {
    if (*i >= n) return false;
    const vial_match_t* m = &e[(*i)++];
    switch (m->kind) {
    case VIAL_MATCH_NOT: return !holds_at(e, n, i, test, user);
    case VIAL_MATCH_AND:
    case VIAL_MATCH_OR: {
        bool r = m->kind == VIAL_MATCH_AND;
        for (uint32_t k = 0; k < m->num; ++k) {
            const bool b = holds_at(e, n, i, test, user);
            r = m->kind == VIAL_MATCH_AND ? r && b : r || b;
        }
        return r;
    }
    default: return test(m, user);
    }
}

bool vial_filter_holds(const vial_filter_t* f, bool (*test)(const vial_match_t* m, const void* user), const void* user) {
    uint32_t i = 0;
    return holds_at(f->expr, (uint32_t)md_array_size(f->expr), &i, test, user);
}

typedef struct str_test_t { md_str_t s; bool glob; } str_test_t;

static bool test_str(const vial_match_t* m, const void* user) {
    const str_test_t* t = user;
    return m->kind == VIAL_MATCH_STR && (t->glob ? vial_glob(m->s, t->s) : md_str_eq(m->s, t->s));
}

static bool match_str(const vial_filter_t* f, md_str_t s, bool glob) {
    const str_test_t t = {s, glob};
    return vial_filter_holds(f, test_str, &t);
}

// Element symbols ignore case: MG, Mg and mg are magnesium, and no two symbols differ only in case
static bool test_element(const vial_match_t* m, const void* user) {
    return m->kind == VIAL_MATCH_STR && md_str_eq_ignore_case(m->s, *(const md_str_t*)user);
}

static bool match_element(const vial_filter_t* f, md_str_t s) {
    return vial_filter_holds(f, test_element, &s);
}

static bool test_int(const vial_match_t* m, const void* user) {
    const int64_t v = *(const int64_t*)user;
    return (m->kind == VIAL_MATCH_INT && m->lo == v) || (m->kind == VIAL_MATCH_RANGE && m->lo <= v && v <= m->hi);
}

static bool match_int(const vial_filter_t* f, int64_t v) {
    return vial_filter_holds(f, test_int, &v);
}

static bool test_part_kind(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const uint8_t k = *(const uint8_t*)user;
    const md_str_t w = m->s;
    if (md_str_eq(w, MD_STR_LIT("polymer"))) return k >= VIAL_PART_PEPTIDE && k <= VIAL_PART_POLYMER;
    if (md_str_eq(w, MD_STR_LIT("nucleic"))) return k == VIAL_PART_DNA || k == VIAL_PART_RNA || k == VIAL_PART_NUCLEIC;
    return md_str_eq(w, md_str_from_cstr(vial_part_kind_word(k)));
}

static bool match_part_kind(const vial_filter_t* f, uint8_t k) {
    return vial_filter_holds(f, test_part_kind, &k);
}

static bool chain_matches(const vial_sys_t* s, const vial_filter_t* f, int32_t c) {
    return c >= 0 && match_str(f, md_system_chain_id(s->sys, (size_t)c), false);
}

// The residue-level filters of a particle selector: resid, resname, chain. A residue that fails is skipped whole.
static bool residue_passes_for_particles(const vial_sys_t* s, const vial_select_t* q, uint32_t r) {
    const md_system_t* sys = s->sys;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_RESID:   if (r == NONE32 || !match_int(f, md_system_residue_seq_id(sys, r))) return false; break;
        case VIAL_KEY_RESNAME: if (r == NONE32 || !match_str(f, md_system_residue_name(sys, r), true)) return false; break;
        case VIAL_KEY_CHAIN:   if (r == NONE32 || !chain_matches(s, f, vial_sys_res_chain(s, r))) return false; break;
        default: break;
        }
    }
    return true;
}

// The particle-level filters, name, element, kind and type, are properties of a particle's type: decided once per
// type, so a particle costs one lookup. NULL when no filter looks at the type.
static uint8_t* type_mask(const vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    bool any = q->particle_kind != 0;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const uint8_t k = q->filters[i].key;
        any = any || k == VIAL_KEY_NAME || k == VIAL_KEY_ELEMENT || k == VIAL_KEY_KIND || k == VIAL_KEY_TYPE;
    }
    if (!any) return NULL;
    const md_system_t* sys = s->sys;
    uint8_t* ok = md_alloc(alloc, s->num_types + 1);
    for (uint32_t t = 0; t < s->num_types; ++t) {
        const uint32_t kind = (uint32_t)md_particle_type_flags_particle_kind(sys->atom.type.flags ? sys->atom.type.flags[t] : 0);
        bool pass = !q->particle_kind || kind + 1 == q->particle_kind;
        for (size_t i = 0; pass && i < md_array_size(q->filters); ++i) {
            const vial_filter_t* f = &q->filters[i];
            switch (f->key) {
            case VIAL_KEY_NAME: pass = match_str(f, MD_LABEL_STR(sys->atom.type.name[t]), true); break;
            case VIAL_KEY_ELEMENT: {
                // A bead or a virtual site has no element: it is none of element='s values, also under not
                const md_atomic_number_t z = sys->atom.type.z ? sys->atom.type.z[t] : 0;
                pass = z != 0 && match_element(f, md_atomic_number_symbol(z));
                break;
            }
            case VIAL_KEY_KIND: pass = kind < MD_COUNTOF(particle_kinds) && match_str(f, md_str_from_cstr(particle_kinds[kind]), false); break;
            case VIAL_KEY_TYPE: {
                // The force field's type, of an atom as of a bead; a particle the data gives no type is none of
                // type='s values, also under not
                const md_str_t ff = md_particle_type_ff_type(&sys->atom.type, t);
                pass = !md_str_empty(ff) && match_str(f, ff, true);
                break;
            }
            default: break;
            }
        }
        ok[t] = pass;
    }
    return ok;
}

static const char* hybridizations[] = {"sp", "sp2", "sp3"};   // md_hybridization_t from 1, as md_hybridization_name spells it

bool vial_select_hybridization_known(md_str_t hyb) {
    return in_list(hybridizations, MD_COUNTOF(hybridizations), hyb);
}

static bool test_bool(const vial_match_t* m, const void* user) {
    return m->kind == VIAL_MATCH_BOOL && (m->lo != 0) == *(const bool*)user;
}

static const char* particle_roles[] = {"backbone", "main", "side", "base", "sugar", "terminus", "n_terminus", "c_terminus",
                                      "5_terminus", "3_terminus", "qm"};

bool vial_select_particle_role_known(md_str_t role) {
    return in_list(particle_roles, MD_COUNTOF(particle_roles), role);
}

static const char* residue_roles[] = {"terminus", "n_terminus", "c_terminus", "5_terminus", "3_terminus"};
static const char* gaps[] = {"before", "after"};
// A residue's secondary structure, DSSP's states as mdlib names them (md_secondary_structure_t), and helix any of the three
static const char* sss[] = {"helix", "helix_alpha", "helix_310", "helix_pi", "sheet", "bridge", "turn", "bend", "coil"};
// A ring's conformation: the kinds of md_pucker_kind_t from planar, md_pucker_kind_name's as words ("half_chair")
static const char* puckers[] = {"planar", "chair", "boat", "skew", "half_chair", "envelope", "twist"};

bool vial_select_residue_role_known(md_str_t role) { return in_list(residue_roles, MD_COUNTOF(residue_roles), role); }
bool vial_select_gap_known(md_str_t gap) { return in_list(gaps, MD_COUNTOF(gaps), gap); }
bool vial_select_ss_known(md_str_t ss) { return in_list(sss, MD_COUNTOF(sss), ss); }
bool vial_select_pucker_known(md_str_t pucker) { return in_list(puckers, MD_COUNTOF(puckers), pucker); }

size_t vial_select_words(uint8_t level, uint8_t key, const char* const** words) {
    #define W(list) (*words = (list), MD_COUNTOF(list))
    switch (key) {
    case VIAL_KEY_KIND:
        switch (level) {
        case VIAL_LEVEL_PARTICLE: return W(particle_kinds);
        case VIAL_LEVEL_RESIDUE:  return W(residue_kinds);
        case VIAL_LEVEL_PART:     return W(part_kinds);
        case VIAL_LEVEL_RING:     return W(ring_kinds);
        default: break;
        }
        break;
    case VIAL_KEY_ROLE:
        if (level == VIAL_LEVEL_PARTICLE) return W(particle_roles);
        if (level == VIAL_LEVEL_RESIDUE) return W(residue_roles);
        break;
    case VIAL_KEY_GAP:           if (level == VIAL_LEVEL_RESIDUE) return W(gaps); break;
    case VIAL_KEY_SS:            if (level == VIAL_LEVEL_RESIDUE) return W(sss); break;
    case VIAL_KEY_PUCKER:        if (level == VIAL_LEVEL_RING) return W(puckers); break;
    case VIAL_KEY_HYBRIDIZATION: if (level == VIAL_LEVEL_PARTICLE) return W(hybridizations); break;
    default: break;
    }
    #undef W
    *words = NULL;
    return 0;
}

// A residue's roles are tags, as a particle's: the ends of its strand, a free amino or carboxyl group of an amino acid
// (n_terminus, c_terminus), a free 5' or 3' end of a nucleotide (MD_RESIDUE_FLAG_TERMINAL_BEG, _END)
typedef struct residue_role_test_t {
    md_residue_flags_t flags;
    bool               nucleotide;
} residue_role_test_t;

static bool test_residue_role(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const residue_role_test_t* r = user;
    const bool beg = r->flags & MD_RESIDUE_FLAG_TERMINAL_BEG, end = r->flags & MD_RESIDUE_FLAG_TERMINAL_END;
    const md_str_t w = m->s;
    if (md_str_eq(w, MD_STR_LIT("terminus")))   return beg || end;
    if (md_str_eq(w, MD_STR_LIT("n_terminus"))) return beg && !r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("c_terminus"))) return end && !r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("5_terminus"))) return beg && r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("3_terminus"))) return end && r->nucleotide;
    return false;
}

// A residue's gaps, tags too: a break of its chain before it, after it (MD_RESIDUE_FLAG_BREAK_BEFORE, _AFTER)
static bool test_gap(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const md_residue_flags_t f = *(const md_residue_flags_t*)user;
    if (md_str_eq(m->s, MD_STR_LIT("before"))) return f & MD_RESIDUE_FLAG_BREAK_BEFORE;
    if (md_str_eq(m->s, MD_STR_LIT("after")))  return f & MD_RESIDUE_FLAG_BREAK_AFTER;
    return false;
}

// A residue's secondary structure, one state with a group: helix is any of alpha, 3-10 and pi
static bool test_ss(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const md_secondary_structure_t ss = *(const md_secondary_structure_t*)user;
    static const char* words[] = {"", "coil", "turn", "bend", "helix_310", "helix_alpha", "helix_pi", "sheet", "bridge"};    // md_secondary_structure_t
    if (md_str_eq(m->s, MD_STR_LIT("helix"))) return ss == MD_SECONDARY_STRUCTURE_HELIX_310 || ss == MD_SECONDARY_STRUCTURE_HELIX_ALPHA || ss == MD_SECONDARY_STRUCTURE_HELIX_PI;
    return (size_t)ss < MD_COUNTOF(words) && ss != MD_SECONDARY_STRUCTURE_UNKNOWN && md_str_eq(m->s, md_str_from_cstr(words[ss]));
}

typedef struct role_test_t {
    md_atom_flags_t flags;
    bool            nucleotide;     // of a nucleotide, whose ends are 5' and 3' rather than N and C
    bool            qm;
} role_test_t;

// A particle's roles are tags, several to a particle: its part in an amino acid or a nucleotide of a strand
// (md_atom_flags_t: what the chain runs through, the backbone, and the side: an amino acid's side chain, a nucleotide's
// nucleoside, its sugar and base), the ends of the strand, and its part in the simulation as the data gives it, the QM
// region. main is the backbone and what of its residue hangs off it that is no side (MD_ATOM_FLAG_MAIN_CHAIN). A
// nucleotide's levels nest, nucleoside over base, so base and sugar are finer than side; its side overlaps the
// backbone (C3' and C4' are both). An amino acid's side starts at CB: HA, on CA, is main, though GROMACS's SideChain,
// every atom of the protein but its main chain, holds it. terminus is either end, of any strand. Every particle has a
// set of roles, empty for most of a solvent, so not takes the complement among all of them: role=not 'qm' is the MM
// part, waters included.
static bool test_particle_role(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const role_test_t* r = user;
    const md_atom_flags_t f = r->flags;
    const md_str_t w = m->s;
    if (md_str_eq(w, MD_STR_LIT("backbone")))   return f & MD_ATOM_FLAG_BACKBONE;
    if (md_str_eq(w, MD_STR_LIT("main")))       return f & MD_ATOM_FLAG_MAIN_CHAIN;
    if (md_str_eq(w, MD_STR_LIT("side")))       return f & (MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_NUCLEOSIDE);
    if (md_str_eq(w, MD_STR_LIT("base")))       return f & MD_ATOM_FLAG_NUCLEOBASE;
    if (md_str_eq(w, MD_STR_LIT("sugar")))      return (f & MD_ATOM_FLAG_NUCLEOSIDE) && !(f & MD_ATOM_FLAG_NUCLEOBASE);
    if (md_str_eq(w, MD_STR_LIT("terminus")))   return f & (MD_ATOM_FLAG_TERMINAL_BEG | MD_ATOM_FLAG_TERMINAL_END);
    if (md_str_eq(w, MD_STR_LIT("n_terminus"))) return (f & MD_ATOM_FLAG_TERMINAL_BEG) && !r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("c_terminus"))) return (f & MD_ATOM_FLAG_TERMINAL_END) && !r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("5_terminus"))) return (f & MD_ATOM_FLAG_TERMINAL_BEG) && r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("3_terminus"))) return (f & MD_ATOM_FLAG_TERMINAL_END) && r->nucleotide;
    if (md_str_eq(w, MD_STR_LIT("qm")))         return r->qm;
    return false;
}

// Whether a particle selector has filters of each particle, role and chemistry, and what they read made ready: the
// metallic phases' atoms
static bool particle_filters(vial_sys_t* s, const vial_select_t* q) {
    bool any = false;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const uint8_t k = q->filters[i].key;
        if (k != VIAL_KEY_ROLE && (k < VIAL_KEY_HYBRIDIZATION || k > VIAL_KEY_FORMAL_CHARGE)) continue;
        any = true;
        if (k == VIAL_KEY_METALLIC && !s->have_metallic) {
            s->have_metallic = true;
            s->metallic = md_alloc(s->alloc, s->num_atoms + 1);
            md_bond_metallic_atoms(s->metallic, &s->sys->bond, s->num_atoms);
        }
    }
    return any;
}

// The filters of a particle selector of each particle: its roles (test_particle_role), and its chemistry, what
// md_chem_perceive found (hybridization, aromatic, the hydrogens, the formal charge) and the metallic phases
// (md_bond_metallic_atoms). Only what has an element has chemistry, and a hybridization or a count of hydrogens only
// where one was perceived (a hydrogen or a metal has no hybridization): what has none is none of a chemistry filter's
// values, also under not, so hybridization=not 'sp3' is the sp and sp2 atoms and never a hydrogen, a bead or a
// virtual site.
static bool atom_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t a) {
    const md_system_t* sys = s->sys;
    const bool chem = md_atom_atomic_number(&sys->atom, a) != 0;
    const md_atom_flags_t flags = md_system_atom_flags(sys, a);
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_ROLE: {
            const uint32_t r = atom_res(s, a);
            const role_test_t t = {
                .flags = flags,
                .nucleotide = r != NONE32 && md_system_residue_kind(sys, r) == MD_RESIDUE_KIND_NUCLEOTIDE,
                .qm = md_atom_is_qm(&sys->atom, a),
            };
            if (!vial_filter_holds(f, test_particle_role, &t)) return false;
            break;
        }
        case VIAL_KEY_HYBRIDIZATION: {
            const md_hybridization_t h = md_atom_flags_hybridization(flags);
            if (!chem || h == MD_HYBRIDIZATION_UNKNOWN || !match_str(f, md_str_from_cstr(md_hybridization_name(h)), false)) return false;
            break;
        }
        case VIAL_KEY_AROMATIC: {
            const bool v = (flags & MD_ATOM_FLAG_AROMATIC) != 0;
            if (!chem || !vial_filter_holds(f, test_bool, &v)) return false;
            break;
        }
        case VIAL_KEY_METALLIC: {
            const bool v = s->metallic && s->metallic[a];
            if (!chem || !vial_filter_holds(f, test_bool, &v)) return false;
            break;
        }
        case VIAL_KEY_HYDROGENS: {
            const int h = md_atom_hydrogen_count(&sys->atom, a);
            if (!chem || h < 0 || !match_int(f, h)) return false;
            break;
        }
        case VIAL_KEY_FORMAL_CHARGE: if (!chem || !match_int(f, md_atom_formal_charge(&sys->atom, a))) return false; break;
        default: break;
        }
    }
    return true;
}

static bool residue_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t r) {
    const md_system_t* sys = s->sys;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_ID:    if (!match_int(f, md_system_residue_seq_id(sys, r))) return false; break;
        case VIAL_KEY_NAME:  if (!match_str(f, md_system_residue_name(sys, r), true)) return false; break;
        case VIAL_KEY_KIND: {
            const uint32_t k = (uint32_t)md_system_residue_kind(sys, r);
            if (k >= MD_COUNTOF(residue_kinds) || !match_str(f, md_str_from_cstr(residue_kinds[k]), false)) return false;
            break;
        }
        case VIAL_KEY_CHAIN: if (!chain_matches(s, f, vial_sys_res_chain(s, r))) return false; break;
        case VIAL_KEY_ROLE: {
            const residue_role_test_t t = {md_system_residue_flags(sys, r), md_system_residue_kind(sys, r) == MD_RESIDUE_KIND_NUCLEOTIDE};
            if (!vial_filter_holds(f, test_residue_role, &t)) return false;
            break;
        }
        case VIAL_KEY_GAP: {
            const md_residue_flags_t flags = md_system_residue_flags(sys, r);
            if (!vial_filter_holds(f, test_gap, &flags)) return false;
            break;
        }
        case VIAL_KEY_RESOLVED: {
            // Of the amino acids and nucleotides, which are resolved or not; any other residue is neither
            const md_residue_kind_t k = md_system_residue_kind(sys, r);
            const bool v = (md_system_residue_flags(sys, r) & MD_RESIDUE_FLAG_RESOLVED) != 0;
            if ((k != MD_RESIDUE_KIND_AMINO_ACID && k != MD_RESIDUE_KIND_NUCLEOTIDE) || !vial_filter_holds(f, test_bool, &v)) return false;
            break;
        }
        case VIAL_KEY_SS: {
            // Of a residue of the protein backbone, in the frame; any other has none, so not 'helix' never reaches it
            const int32_t seg = s->have_res_segs ? s->res_seg[r] : -1;
            const md_secondary_structure_t ss = q->ss && seg >= 0 ? (md_secondary_structure_t)q->ss[seg] : MD_SECONDARY_STRUCTURE_UNKNOWN;
            if (ss == MD_SECONDARY_STRUCTURE_UNKNOWN || !vial_filter_holds(f, test_ss, &ss)) return false;
            break;
        }
        default: break;
        }
    }
    return true;
}

// Its type's description, as the data gives it; an inferred type's is the kind of a strand or water, as part(kind=)
// spells it (md_part_kind_name), and the residue names of anything else
md_str_t vial_part_name(const md_system_t* sys, uint32_t p) {
    const md_part_type_idx_t t = md_system_part_type_idx(sys, p);
    const md_str_t d = t >= 0 ? md_system_part_type_description(sys, (size_t)t) : (md_str_t){0};
    return md_str_empty(d) ? md_str_from_cstr(vial_part_kind_word(part_kind(sys, p))) : d;
}

const char* md_vial_part_kind(const md_system_t* sys, size_t part_idx) {
    return sys && part_idx < sys->part.count ? vial_part_kind_word(part_kind(sys, (uint32_t)part_idx)) : "";
}

md_str_t md_vial_part_name(const md_system_t* sys, size_t part_idx) {
    return sys && part_idx < sys->part.count ? vial_part_name(sys, (uint32_t)part_idx) : (md_str_t){0};
}

static bool part_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t p) {
    if (q->particle_kind && !(s->part_kind[p] >= VIAL_PART_PEPTIDE && s->part_kind[p] <= VIAL_PART_POLYMER)) return false;    // polymer()
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_KIND:  if (!match_part_kind(f, s->part_kind[p])) return false; break;
        case VIAL_KEY_CHAIN: if (!chain_matches(s, f, s->part_chain[p])) return false; break;
        case VIAL_KEY_NAME:  if (!match_str(f, vial_part_name(s->sys, p), true)) return false; break;
        default: break;
        }
    }
    return true;
}

// A ring's kinds are tags, several to a ring (md_ring_flags_t): aromatic (every bond aromatic), hetero (an atom other than
// C), fused (a bond shared with another ring), fused_aromatic (an aromatic ring of a fused aromatic system, which shares a
// bond with another aromatic ring: naphthalene, indole, a purine, whose sextets are not assigned to rings), sugar (a
// furanose or pyranose). What the others combine to, an aliphatic ring or a heteroaromatic one, is written with
// and/or/not: not 'aromatic', 'hetero' and 'aromatic'.
static bool test_ring_kind(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const uint32_t flags = *(const uint32_t*)user;
    const md_str_t w = m->s;
    if (md_str_eq(w, MD_STR_LIT("aromatic")))       return flags & MD_RING_FLAG_AROMATIC;
    if (md_str_eq(w, MD_STR_LIT("hetero")))         return flags & MD_RING_FLAG_HETERO;
    if (md_str_eq(w, MD_STR_LIT("fused")))          return flags & MD_RING_FLAG_FUSED;
    if (md_str_eq(w, MD_STR_LIT("fused_aromatic"))) return (flags & MD_RING_FLAG_AROMATIC) && (flags & MD_RING_FLAG_FUSED_AROMATIC);
    if (md_str_eq(w, MD_STR_LIT("sugar")))          return flags & MD_RING_FLAG_SUGAR;
    return false;
}

// kind= its kinds (test_ring_kind); size= the number of atoms
// A ring's conformation, one kind (md_pucker_kind_t, from planar at 1)
static bool test_pucker(const vial_match_t* m, const void* user) {
    const uint8_t k = *(const uint8_t*)user;
    return m->kind == VIAL_MATCH_STR && k >= 1 && k <= MD_COUNTOF(puckers) && md_str_eq(m->s, md_str_from_cstr(puckers[k - 1]));
}

static bool ring_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t r) {
    const uint32_t flags = md_system_ring_flags(s->sys, r);
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        if (f->key == VIAL_KEY_SIZE) {
            if (!match_int(f, (int64_t)(s->ring_offset[r + 1] - s->ring_offset[r]))) return false;
        } else if (f->key == VIAL_KEY_KIND) {
            if (!vial_filter_holds(f, test_ring_kind, &flags)) return false;
        } else if (f->key == VIAL_KEY_PUCKER) {
            // Of a ring of 5 or 6 atoms, in the frame; another has none, so not 'chair' never reaches it
            const uint8_t k = q->pucker ? q->pucker[r] : 0;
            if (!k || !vial_filter_holds(f, test_pucker, &k)) return false;
        }
    }
    return true;
}

// A fragment by its name (md_system_fragment_name), the one keyword it takes besides scope
static bool fragment_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t f) {
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* fl = &q->filters[i];
        if (fl->key == VIAL_KEY_NAME && !match_str(fl, md_system_fragment_name(s->sys, s->frag_idx[f]), true)) return false;
    }
    return true;
}

// A group's roles are tags, several to a group (md_group_role_name)
static bool test_group_role(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const uint32_t roles = *(const uint32_t*)user;
    for (int r = 0; r < MD_GROUP_ROLE_COUNT; ++r) {
        if ((roles & (1u << r)) && md_str_eq(m->s, md_str_from_cstr(md_group_role_name((md_group_role_t)r)))) return true;
    }
    return false;
}

// name= its name, a glob; role= its roles (test_group_role). A group of an index file has no role, so no role= holds
// for it, also under not
static bool group_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t u) {
    const md_group_data_t* g = &s->sys->group;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        if (f->key == VIAL_KEY_NAME) {
            if (!match_str(f, g->name[u], true)) return false;
        } else if (f->key == VIAL_KEY_ROLE) {
            const uint32_t roles = g->roles[u];
            if (!roles || !vial_filter_holds(f, test_group_role, &roles)) return false;
        }
    }
    return true;
}

bool vial_select_role_known(md_str_t role) {
    for (int r = 0; r < MD_GROUP_ROLE_COUNT; ++r) if (md_str_eq(role, md_str_from_cstr(md_group_role_name((md_group_role_t)r)))) return true;
    return false;
}

// The order the groups are taken in: the system's, or a role's when role= names one, as GROMACS orders them (the order of
// energygrps is the order its energy file names the pairs in). n receives the number of units.
static uint32_t* group_order(vial_sys_t* s, const vial_select_t* q, uint32_t* n, md_allocator_i* alloc) {
    need_groups(s);
    const md_group_data_t* g = &s->sys->group;
    uint32_t* order = md_alloc(alloc, (s->num_groups + 1) * sizeof(uint32_t));
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        if (f->key != VIAL_KEY_ROLE || md_array_size(f->expr) != 1 || f->expr[0].kind != VIAL_MATCH_STR) continue;
        for (int r = 0; r < MD_GROUP_ROLE_COUNT; ++r) {
            if (!md_str_eq(f->expr[0].s, md_str_from_cstr(md_group_role_name((md_group_role_t)r)))) continue;
            size_t num;
            const uint32_t* role_order = md_system_group_role_order(g, (md_group_role_t)r, &num);
            uint32_t k = 0;
            for (size_t j = 0; j < num && k < s->num_groups; ++j) order[k++] = role_order[j];
            *n = k;
            return order;
        }
    }
    for (uint32_t i = 0; i < s->num_groups; ++i) order[i] = i;
    *n = s->num_groups;
    return order;
}

static bool molecule_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t m) {
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_NAME: {
            const md_molecule_type_idx_t t = md_system_molecule_type_idx(s->sys, m);
            if (t < 0 || !match_str(f, md_system_molecule_type_name(s->sys, (size_t)t), true)) return false;
            break;
        }
        case VIAL_KEY_CHAIN: if (!chain_matches(s, f, s->mol_chain[m])) return false; break;     // need_mols ran
        default: break;
        }
    }
    return true;
}

static bool chain_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t c) {
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        if (f->key == VIAL_KEY_ID && !match_str(f, md_system_chain_id(s->sys, c), false)) return false;
    }
    return true;
}

// A bitmap of particles, for scope and positional membership
static uint64_t* bitmap(const vial_sets_t* s, uint32_t num_atoms, md_allocator_i* alloc) {
    const size_t words = num_atoms / 64 + 1;
    uint64_t* b = md_alloc(alloc, words * sizeof(uint64_t));
    memset(b, 0, words * sizeof(uint64_t));
    const uint32_t total = s->offset[s->count];
    for (uint32_t i = 0; i < total; ++i) if (s->idx[i] < num_atoms) b[s->idx[i] / 64] |= 1ull << (s->idx[i] % 64);
    return b;
}
static inline bool bit(const uint64_t* b, uint32_t i) { return (b[i / 64] >> (i % 64)) & 1; }

// particle(s, ...) of one set: only s's particles are looked at, so selecting within a residue costs its size, not
// the system's. Residues are found by walking s, which is ascending.
static vial_sets_t* select_particles_of(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    const md_system_t* sys = s->sys;
    const uint32_t* items = vial_sets_items(q->within, 0);
    const uint32_t n = vial_sets_size(q->within, 0);
    const uint8_t* types = type_mask(s, q, alloc);
    const bool each = particle_filters(s, q);
    const uint64_t* scope = q->scope ? bitmap(q->scope, s->num_atoms, alloc) : NULL;
    vial_sets_t* out = sets_new(1, n, alloc);
    out->single = true;
    uint32_t w = 0, r = NONE32, r_end = 0;
    bool r_pass = false;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t a = items[k];
        if (a >= s->num_atoms || (scope && !bit(scope, a))) continue;
        if (r == NONE32 || a >= r_end) {
            r = atom_res(s, a);
            r_end = r != NONE32 ? md_system_residue_atom_range(sys, r).end : a + 1;
            r_pass = residue_passes_for_particles(s, q, s->num_res ? r : NONE32);
        }
        if (!r_pass) continue;
        if (types) {
            const md_particle_type_idx_t t = sys->atom.type_idx[a];
            if ((uint32_t)t >= s->num_types || !types[t]) continue;
        }
        if (each && !atom_passes(s, q, a)) continue;
        out->idx[w++] = a;
    }
    out->offset[1] = w;
    return out;
}

vial_sets_t* vial_select(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    if (q->ss) need_res_segs(s);
    if (q->level == VIAL_LEVEL_PARTICLE && q->within && q->within->single) return select_particles_of(s, q, alloc);
    const uint64_t* scope = q->scope ? bitmap(q->scope, s->num_atoms, alloc) : NULL;
    const uint64_t* within = q->within ? bitmap(q->within, s->num_atoms, alloc) : NULL;
    #define IN_SCOPE(a) (!scope || bit(scope, (a)))

    if (q->level == VIAL_LEVEL_PARTICLE) {
        // Residue by residue: a residue that fails resid, resname or chain is skipped whole, and a particle costs a
        // lookup of its type's mask
        const md_system_t* sys = s->sys;
        const uint8_t* types = type_mask(s, q, alloc);
        const bool each = particle_filters(s, q);
        const uint32_t num_units = s->num_res ? s->num_res : 1;
        // The residues that pass, and room for their particles: a selection costs its size, not the system's
        uint32_t room = 0;
        uint8_t* pass = md_alloc(alloc, num_units);
        for (uint32_t r = 0; r < num_units; ++r) {
            pass[r] = residue_passes_for_particles(s, q, s->num_res ? r : NONE32);
            if (pass[r]) room += s->num_res ? (uint32_t)md_system_residue_atom_count(sys, r) : s->num_atoms;
        }
        vial_sets_t* out = sets_new(1, room, alloc);
        out->single = true;
        uint32_t n = 0;
        for (uint32_t r = 0; r < num_units; ++r) {
            uint32_t beg = 0, end = s->num_atoms;
            if (s->num_res) {
                const md_urange_t range = md_system_residue_atom_range(sys, r);
                beg = range.beg; end = range.end < s->num_atoms ? range.end : s->num_atoms;
            }
            if (!pass[r]) continue;
            for (uint32_t a = beg; a < end && n < room; ++a) {
                if (!IN_SCOPE(a) || (within && !bit(within, a))) continue;
                if (types) {
                    const md_particle_type_idx_t t = sys->atom.type_idx[a];
                    if ((uint32_t)t >= s->num_types || !types[t]) continue;
                }
                if (each && !atom_passes(s, q, a)) continue;
                out->idx[n++] = a;
            }
        }
        out->offset[1] = n;
        return out;
    }

    // The units, each with its particles inside the scope
    uint32_t num_units = 0;
    switch (q->level) {
    case VIAL_LEVEL_RESIDUE:  num_units = s->num_res; break;
    case VIAL_LEVEL_PART:     num_units = s->num_parts; break;
    case VIAL_LEVEL_MOLECULE: num_units = s->num_mols; need_mols(s); break;
    case VIAL_LEVEL_CHAIN:    num_units = s->num_chains; need_chains(s); break;
    case VIAL_LEVEL_TAIL:     need_tails(s); num_units = s->num_tails; break;
    case VIAL_LEVEL_HEAD:     need_heads(s); num_units = s->num_heads; break;
    case VIAL_LEVEL_RING:     need_rings(s); num_units = s->num_rings; break;
    case VIAL_LEVEL_FRAGMENT: need_frags(s); num_units = s->num_frags; break;
    case VIAL_LEVEL_GROUP:    break;
    default: break;
    }
    const uint32_t* order = q->level == VIAL_LEVEL_GROUP ? group_order(s, q, &num_units, alloc) : NULL;
    // The units that pass, and room for their particles: the units of a level never share a particle
    uint32_t* passing = md_alloc(alloc, (num_units + 1) * sizeof(uint32_t));
    uint32_t num_pass = 0, room = 0;
    for (uint32_t ui = 0; ui < num_units; ++ui) {
        const uint32_t u = order ? order[ui] : ui;
        bool pass = false;
        uint32_t size = 0;
        switch (q->level) {
        case VIAL_LEVEL_RESIDUE:  pass = residue_passes(s, q, u);  size = (uint32_t)md_system_residue_atom_count(s->sys, u); break;
        case VIAL_LEVEL_PART:     pass = part_passes(s, q, u);     { const md_urange_t r = md_system_part_atom_range(s->sys, u); size = r.end - r.beg; } break;
        case VIAL_LEVEL_MOLECULE: pass = molecule_passes(s, q, u); size = s->mol_offset[u + 1] - s->mol_offset[u]; break;
        case VIAL_LEVEL_CHAIN:    pass = chain_passes(s, q, u);    size = s->chain_offset[u + 1] - s->chain_offset[u]; break;
        case VIAL_LEVEL_TAIL:     pass = true;                     size = s->tail_offset[u + 1] - s->tail_offset[u]; break;
        case VIAL_LEVEL_HEAD:     pass = true;                     size = s->head_offset[u + 1] - s->head_offset[u]; break;
        case VIAL_LEVEL_RING:     pass = ring_passes(s, q, u);     size = s->ring_offset[u + 1] - s->ring_offset[u]; break;
        case VIAL_LEVEL_GROUP:    pass = group_passes(s, q, u);    size = s->group_offset[u + 1] - s->group_offset[u]; break;
        case VIAL_LEVEL_FRAGMENT: pass = fragment_passes(s, q, u); size = s->frag_offset[u + 1] - s->frag_offset[u]; break;
        default: break;
        }
        if (!pass) continue;
        passing[num_pass++] = u;
        room += size;
    }
    vial_sets_t* out = sets_new(0, room, alloc);
    out->offset = md_alloc(alloc, (num_pass + 1) * sizeof(uint32_t));
    out->unit = md_alloc(alloc, (num_pass + 1) * sizeof(uint32_t));
    out->offset[0] = 0;
    out->level = q->level;
    uint32_t n = 0, w = 0;
    for (uint32_t pi = 0; pi < num_pass; ++pi) {
        const uint32_t u = passing[pi];
        // Its particles: a range, or a list for molecules and chains
        const uint32_t* list = NULL;
        uint32_t beg = 0, end = 0;
        switch (q->level) {
        case VIAL_LEVEL_RESIDUE: { const md_urange_t r = md_system_residue_atom_range(s->sys, u); beg = r.beg; end = r.end; break; }
        case VIAL_LEVEL_PART:    { const md_urange_t r = md_system_part_atom_range(s->sys, u); beg = r.beg; end = r.end; break; }
        case VIAL_LEVEL_MOLECULE: list = s->mol_atoms; beg = s->mol_offset[u]; end = s->mol_offset[u + 1]; break;
        case VIAL_LEVEL_CHAIN:    list = s->chain_atoms; beg = s->chain_offset[u]; end = s->chain_offset[u + 1]; break;
        case VIAL_LEVEL_TAIL:     list = s->tail_atoms; beg = s->tail_offset[u]; end = s->tail_offset[u + 1]; break;
        case VIAL_LEVEL_HEAD:     list = s->head_atoms; beg = s->head_offset[u]; end = s->head_offset[u + 1]; break;
        case VIAL_LEVEL_RING:     list = s->ring_atoms; beg = s->ring_offset[u]; end = s->ring_offset[u + 1]; break;
        case VIAL_LEVEL_GROUP:    list = s->group_atoms; beg = s->group_offset[u]; end = s->group_offset[u + 1]; break;
        case VIAL_LEVEL_FRAGMENT: list = s->frag_atoms; beg = s->frag_offset[u]; end = s->frag_offset[u + 1]; break;
        default: break;
        }
        const uint32_t mark = w;
        bool hit = !within;
        for (uint32_t k = beg; k < end && w < room; ++k) {
            const uint32_t a = list ? list[k] : k;
            if (a >= s->num_atoms || !IN_SCOPE(a)) continue;
            if (within && bit(within, a)) hit = true;
            out->idx[w++] = a;
        }
        if (w == mark || !hit) { w = mark; continue; }    // outside the scope, or not containing
        out->unit[n] = u;
        out->offset[++n] = w;
    }
    #undef IN_SCOPE
    out->count = n;
    return out;
}

// ### BONDS ###

bool vial_select_bond_kind_known(md_str_t kind) {
    size_t n = 0;
    const md_bond_kind_name_t* k = md_bond_kind_names(&n);
    for (size_t i = 0; i < n; ++i) if (md_str_eq(kind, md_str_from_cstr(k[i].name))) return true;
    return false;
}

const char* vial_select_bond_kinds(void) { return "covalent, coordinate, elastic, metallic, interaction, contact, hydrogen, ionic and halogen"; }

static const char* bond_order_words[] = {"aromatic", "delocalized"};    // the resonance, MD_BOND_FLAG_AROMATIC, _DELOCALIZED

bool vial_select_bond_order_known(md_str_t word) { return in_list(bond_order_words, MD_COUNTOF(bond_order_words), word); }

bool vial_select_bond_origin_known(md_str_t origin) {
    for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED; ++o) if (md_str_eq(origin, md_str_from_cstr(md_bond_origin_name((md_bond_origin_t)o)))) return true;
    return false;
}

// A bond's kind, by md_bond_kind_names: 'covalent' is the bond graph, covalent and coordinate bonds, and 'interaction'
// the four interactions
static bool test_bond_kind(const vial_match_t* m, const void* user) {
    if (m->kind != VIAL_MATCH_STR) return false;
    const md_bond_flags_t f = *(const md_bond_flags_t*)user;
    size_t n = 0;
    const md_bond_kind_name_t* k = md_bond_kind_names(&n);
    for (size_t i = 0; i < n; ++i) if (md_str_eq(m->s, md_str_from_cstr(k[i].name))) return (f & k[i].mask) == k[i].value;
    return false;
}

bool vial_select_bond_kind_possible(const vial_filter_t* f) {
    for (unsigned int v = 0; v <= (unsigned int)MD_BOND_FLAG_KIND_MASK; v += 0x20u) {
        const md_bond_flags_t kind = (md_bond_flags_t)v;
        if (vial_filter_holds(f, test_bond_kind, &kind)) return true;
    }
    return false;
}

// A bond's order: a word for one of resonance, as a file writes an aromatic bond's order (SDF's 4, MOL2's ar) and
// RDKit's types have it, else its order, 1 to 4. The order mdlib holds beneath the resonance is one Kekule structure,
// and no order of the bond.
static bool test_bond_order(const vial_match_t* m, const void* user) {
    const md_bond_flags_t f = *(const md_bond_flags_t*)user;
    const md_bond_flags_t res = f & MD_BOND_FLAG_RESONANCE_MASK;
    if (m->kind == VIAL_MATCH_STR) {
        return (res == MD_BOND_FLAG_AROMATIC && md_str_eq(m->s, MD_STR_LIT("aromatic"))) ||
               (res == MD_BOND_FLAG_DELOCALIZED && md_str_eq(m->s, MD_STR_LIT("delocalized")));
    }
    if (res) return false;
    const int64_t o = md_bond_order(f);
    return (m->kind == VIAL_MATCH_INT && m->lo == o) || (m->kind == VIAL_MATCH_RANGE && m->lo <= o && o <= m->hi);
}

static bool test_bond_origin(const vial_match_t* m, const void* user) {
    const md_bond_flags_t f = *(const md_bond_flags_t*)user;
    return m->kind == VIAL_MATCH_STR && md_str_eq(m->s, md_str_from_cstr(md_bond_origin_name(md_bond_origin(f))));
}

// Whether a bond passes the filters. not takes the complement among the bonds with the property: every bond has a kind
// and an origin, and a bond of unknown order, as an interaction, has none, so order=not 1 is never one of those.
static bool bond_passes(md_bond_flags_t f, const vial_filter_t* filters, size_t num_filters) {
    bool kind_given = false;
    for (size_t i = 0; i < num_filters; ++i) {
        const vial_filter_t* fl = &filters[i];
        switch (fl->key) {
        case VIAL_KEY_KIND:   kind_given = true; if (!vial_filter_holds(fl, test_bond_kind, &f)) return false; break;
        case VIAL_KEY_ORDER:  if ((!(f & MD_BOND_FLAG_RESONANCE_MASK) && md_bond_order(f) == 0) || !vial_filter_holds(fl, test_bond_order, &f)) return false; break;
        case VIAL_KEY_ORIGIN: if (!vial_filter_holds(fl, test_bond_origin, &f)) return false; break;
        default: break;
        }
    }
    // The strong bonds by default: the bond graph, which the molecules, rings and chemistry are made of
    return kind_given || md_bond_is_graph(f);
}

uint32_t vial_select_bonds_found(const vial_filter_t* filters, size_t num_filters) {
    static const struct { md_bond_kind_t kind; uint32_t bit; } found[] = {
        {MD_BOND_KIND_HYDROGEN, VIAL_FOUND_HYDROGEN}, {MD_BOND_KIND_IONIC, VIAL_FOUND_IONIC}, {MD_BOND_KIND_HALOGEN, VIAL_FOUND_HALOGEN},
    };
    uint32_t bits = 0;
    for (size_t i = 0; i < MD_COUNTOF(found); ++i) {
        const md_bond_flags_t f = md_bond_flags_set_origin(md_bond_flags_set_kind(MD_BOND_FLAG_NONE, found[i].kind), MD_BOND_ORIGIN_INFERRED);
        if (bond_passes(f, filters, num_filters)) bits |= found[i].bit;
    }
    return bits;
}

vial_sets_t* vial_select_bonds(const vial_sys_t* s, const vial_sets_t* within, const vial_filter_t* filters, size_t num_filters,
                               md_allocator_i* alloc) {
    const md_system_t* sys = s->sys;
    const md_bond_data_t* bd = &sys->bond;
    const size_t total = bd->count + bd->interaction_count;
    // Both ends in the set: by a mark per particle
    uint8_t* in = NULL;
    if (within) {
        in = md_alloc(alloc, s->num_atoms + 1);
        memset(in, 0, s->num_atoms + 1);
        for (uint32_t u = 0; u < within->count; ++u) {
            const uint32_t* it = vial_sets_items(within, u);
            for (uint32_t k = 0, n = vial_sets_size(within, u); k < n; ++k) if (it[k] < s->num_atoms) in[it[k]] = 1;
        }
    }
    md_array(uint32_t) idx = 0;
    for (size_t i = 0; i < total; ++i) {
        const md_atom_pair_t pr = bd->pairs[i];
        if (pr.idx[0] < 0 || pr.idx[1] < 0 || (uint32_t)pr.idx[0] >= s->num_atoms || (uint32_t)pr.idx[1] >= s->num_atoms) continue;
        if (in && (!in[pr.idx[0]] || !in[pr.idx[1]])) continue;
        if (!bond_passes(bd->flags ? bd->flags[i] : 0, filters, num_filters)) continue;
        md_array_push(idx, (uint32_t)pr.idx[0], alloc);
        md_array_push(idx, (uint32_t)pr.idx[1], alloc);
    }
    const uint32_t n2 = (uint32_t)md_array_size(idx);
    vial_sets_t* r = md_alloc(alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = n2;
    r->idx = idx ? idx : md_alloc(alloc, sizeof(uint32_t));
    r->offset = md_alloc(alloc, (n2 + 1) * sizeof(uint32_t));
    for (uint32_t k = 0; k <= n2; ++k) r->offset[k] = k;
    return r;
}

vial_sets_t* vial_select_stored_bonds(const vial_sys_t* s, md_allocator_i* alloc) {
    const md_bond_data_t* bd = &s->sys->bond;
    const uint32_t n2 = (uint32_t)(2 * (bd->count + bd->interaction_count));
    vial_sets_t* r = md_alloc(alloc, sizeof(vial_sets_t));
    memset(r, 0, sizeof(*r));
    r->count = n2;
    r->idx = md_alloc(alloc, (n2 + 1) * sizeof(uint32_t));
    r->offset = md_alloc(alloc, (n2 + 1) * sizeof(uint32_t));
    for (uint32_t i = 0; i < n2 / 2; ++i) {
        r->idx[2 * i]     = (uint32_t)bd->pairs[i].idx[0];
        r->idx[2 * i + 1] = (uint32_t)bd->pairs[i].idx[1];
    }
    for (uint32_t k = 0; k <= n2; ++k) r->offset[k] = k;
    return r;
}

vial_sets_t* vial_select_predefined(vial_sys_t* s, md_str_t name, md_allocator_i* alloc) {
    if (md_str_eq(name, MD_STR_LIT("system"))) return vial_sets_all(s->num_atoms, alloc);
    const md_system_t* sys = s->sys;
    // By the flags of the particles: what the structure's perception marked
    md_atom_flags_t flag = 0;
    if (md_str_eq(name, MD_STR_LIT("backbone")))        flag = MD_ATOM_FLAG_BACKBONE;
    else if (md_str_eq(name, MD_STR_LIT("sidechain")))  flag = MD_ATOM_FLAG_SIDE_CHAIN;
    else if (md_str_eq(name, MD_STR_LIT("nucleobase"))) flag = MD_ATOM_FLAG_NUCLEOBASE;
    else if (md_str_eq(name, MD_STR_LIT("nucleoside"))) flag = MD_ATOM_FLAG_NUCLEOSIDE;
    int part = -1;
    if (md_str_eq(name, MD_STR_LIT("protein"))) part = 0;
    else if (md_str_eq(name, MD_STR_LIT("nucleic"))) part = 1;
    else if (md_str_eq(name, MD_STR_LIT("water"))) part = 2;
    else if (md_str_eq(name, MD_STR_LIT("ion"))) part = 3;
    else if (md_str_eq(name, MD_STR_LIT("lipid"))) part = 4;
    if (!flag && part < 0) return NULL;
    vial_sets_t* out = sets_new(1, s->num_atoms, alloc);
    out->single = true;
    uint32_t n = 0;
    if (flag) {
        for (uint32_t a = 0; a < s->num_atoms; ++a) if (md_system_atom_flags(sys, a) & flag) out->idx[n++] = a;
    } else {
        // The particles of the parts of that kind: each part is one range, and parts ascend
        for (uint32_t p = 0; p < s->num_parts; ++p) {
            const uint8_t k = s->part_kind[p];
            bool in = false;
            switch (part) {
            case 0: in = k == VIAL_PART_PEPTIDE; break;
            case 1: in = k == VIAL_PART_DNA || k == VIAL_PART_RNA || k == VIAL_PART_NUCLEIC; break;
            case 2: in = k == VIAL_PART_WATER; break;
            case 3: in = k == VIAL_PART_ION; break;
            case 4: in = k == VIAL_PART_LIPID; break;
            default: break;
            }
            if (!in) continue;
            const md_urange_t r = md_system_part_atom_range(sys, p);
            for (uint32_t a = r.beg; a < r.end && a < s->num_atoms; ++a) out->idx[n++] = a;
        }
    }
    out->offset[1] = n;
    return out;
}

uint32_t vial_select_resid_chains(const vial_sys_t* s, const vial_filter_t* f, const vial_sets_t* scope, int32_t out[2]) {
    md_allocator_i* temp = md_get_heap_allocator();
    uint64_t* in = scope ? bitmap(scope, s->num_atoms, temp) : NULL;
    uint32_t n = 0;
    for (uint32_t r = 0; r < s->num_res && n < 2; ++r) {
        if (!match_int(f, md_system_residue_seq_id(s->sys, r))) continue;
        if (in) {
            // Inside the scope when any of its particles is
            const md_urange_t a = md_system_residue_atom_range(s->sys, r);
            bool hit = false;
            for (uint32_t k = a.beg; k < a.end && !hit; ++k) hit = bit(in, k);
            if (!hit) continue;
        }
        const int32_t c = vial_sys_res_chain(s, r);
        if (c < 0) continue;
        if (n == 0 || out[0] != c) out[n++] = c;
    }
    if (in) md_free(temp, in, (s->num_atoms / 64 + 1) * sizeof(uint64_t));
    return n;
}

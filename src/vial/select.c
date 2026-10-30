// The system as the selectors see it, the static values of sets, and the selectors themselves.

#include "select.h"

#include <md_system.h>
#include <md_util.h>
#include <core/md_common.h>
#include <core/md_array.h>

#include <stdio.h>
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

int32_t vial_sys_res_chain(const vial_sys_t* s, uint32_t r) {
    const md_system_t* sys = s->sys;
    if (!s->num_chains || !sys->chain.res_chain || r >= s->num_res) return -1;
    const int32_t c = sys->chain.res_chain[r];
    return c >= 0 && (uint32_t)c < s->num_chains ? c : -1;
}

str_t vial_sys_chain_id(const vial_sys_t* s, int32_t chain) {
    return chain >= 0 && (uint32_t)chain < s->num_chains ? md_system_chain_id(s->sys, (size_t)chain) : STR_LIT("");
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
    const str_t name = md_system_residue_name(s->sys, r);
    const int32_t c = vial_sys_res_chain(s, r);
    if (c >= 0) {
        const str_t ch = md_system_chain_id(s->sys, (size_t)c);
        return (size_t)snprintf(buf, cap, "%.*s%d/%.*s", (int)name.len, name.ptr, (int)md_system_residue_seq_id(s->sys, r), (int)ch.len, ch.ptr);
    }
    return (size_t)snprintf(buf, cap, "%.*s%d", (int)name.len, name.ptr, (int)md_system_residue_seq_id(s->sys, r));
}

static const char* part_kind_name(uint8_t k) {
    static const char* names[] = {"ligand", "peptide", "dna", "rna", "nucleic", "polymer", "glycan", "lipid", "ion", "water"};
    return k < ARRAY_SIZE(names) ? names[k] : "part";
}

size_t vial_unit_label(char* buf, size_t cap, vial_sys_t* s, uint8_t level, uint32_t unit) {
    if (cap) buf[0] = '\0';
    const md_system_t* sys = s->sys;
    int n = 0;
    switch (level) {
    case VIAL_LEVEL_PARTICLE: {
        const md_particle_type_idx_t t = sys->atom.type_idx ? sys->atom.type_idx[unit] : -1;
        const str_t name = t >= 0 ? LBL_TO_STR(sys->atom.type.name[t]) : STR_LIT("");
        n = snprintf(buf, cap, "%.*s/", (int)name.len, name.ptr);
        const uint32_t r = atom_res(s, unit);
        if (n > 0 && (size_t)n < cap && r != NONE32) n += (int)residue_label(buf + n, cap - n, s, r);
        break;
    }
    case VIAL_LEVEL_RESIDUE: return residue_label(buf, cap, s, unit);
    case VIAL_LEVEL_CHAIN: {
        const str_t id = md_system_chain_id(sys, unit);
        n = snprintf(buf, cap, "%.*s", (int)id.len, id.ptr);
        break;
    }
    case VIAL_LEVEL_PART:
    case VIAL_LEVEL_MOLECULE: {
        // One residue: that residue; else what it is, and its chain
        uint32_t r0 = NONE32, nres = 0;
        int32_t chain = -1;
        str_t what = STR_LIT("");
        if (level == VIAL_LEVEL_PART) {
            const md_urange_t res = md_system_part_residue_range(sys, unit);
            r0 = res.beg; nres = res.end - res.beg;
            chain = s->part_chain[unit];
            what = str_from_cstr(part_kind_name(s->part_kind[unit]));
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
            const str_t ch = md_system_chain_id(sys, (size_t)chain);
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

// ### SELECTORS ###

static const struct { const char* word; uint8_t key; } keywords[] = {
    {"name", VIAL_KEY_NAME}, {"serial", VIAL_KEY_SERIAL}, {"element", VIAL_KEY_ELEMENT}, {"kind", VIAL_KEY_KIND},
    {"resid", VIAL_KEY_RESID}, {"resname", VIAL_KEY_RESNAME}, {"chain", VIAL_KEY_CHAIN}, {"id", VIAL_KEY_ID},
    {"scope", VIAL_KEY_SCOPE},
};

uint8_t vial_select_key(str_t keyword) {
    for (size_t i = 0; i < ARRAY_SIZE(keywords); ++i) if (str_eq(keyword, str_from_cstr(keywords[i].word))) return keywords[i].key;
    return VIAL_KEY_NONE;
}

#define K(k) (1u << (k))
static uint32_t takes(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return K(VIAL_KEY_NAME) | K(VIAL_KEY_SERIAL) | K(VIAL_KEY_ELEMENT) | K(VIAL_KEY_KIND) |
                                     K(VIAL_KEY_RESID) | K(VIAL_KEY_RESNAME) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_RESIDUE:  return K(VIAL_KEY_ID) | K(VIAL_KEY_NAME) | K(VIAL_KEY_KIND) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_PART:     return K(VIAL_KEY_KIND) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_MOLECULE: return K(VIAL_KEY_NAME) | K(VIAL_KEY_CHAIN) | K(VIAL_KEY_SCOPE);
    case VIAL_LEVEL_CHAIN:    return K(VIAL_KEY_ID) | K(VIAL_KEY_SCOPE);
    default: return 0;
    }
}
#undef K

bool vial_select_takes(uint8_t level, uint8_t key) { return key < 32 && (takes(level) >> key) & 1; }

const char* vial_select_keywords(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "name, serial, element, kind, resid, resname, chain and scope";
    case VIAL_LEVEL_RESIDUE:  return "id, name, kind, chain and scope";
    case VIAL_LEVEL_PART:     return "kind, chain and scope";
    case VIAL_LEVEL_MOLECULE: return "name, chain and scope";
    case VIAL_LEVEL_CHAIN:    return "id and scope";
    default: return "";
    }
}

const char* vial_select_name(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "particles";
    case VIAL_LEVEL_RESIDUE:  return "residues";
    case VIAL_LEVEL_PART:     return "parts";
    case VIAL_LEVEL_MOLECULE: return "molecules";
    case VIAL_LEVEL_CHAIN:    return "chains";
    default: return "sets";
    }
}

static const char* residue_kinds[] = {"other", "amino_acid", "nucleotide", "water", "ion", "lipid"};   // md_residue_kind_t
static const char* particle_kinds[] = {"atom", "bead", "virtual"};                                    // md_particle_kind_t
static const char* part_kinds[] = {"polymer", "peptide", "nucleic", "dna", "rna", "glycan", "ligand", "lipid", "ion", "water"};

static bool in_list(const char** list, size_t n, str_t s) {
    for (size_t i = 0; i < n; ++i) if (str_eq(s, str_from_cstr(list[i]))) return true;
    return false;
}

bool vial_select_kind_known(uint8_t level, str_t kind) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return in_list(particle_kinds, ARRAY_SIZE(particle_kinds), kind);
    case VIAL_LEVEL_RESIDUE:  return in_list(residue_kinds, ARRAY_SIZE(residue_kinds), kind);
    case VIAL_LEVEL_PART:     return in_list(part_kinds, ARRAY_SIZE(part_kinds), kind);
    default: return false;
    }
}

const char* vial_select_kinds(uint8_t level) {
    switch (level) {
    case VIAL_LEVEL_PARTICLE: return "atom, bead and virtual";
    case VIAL_LEVEL_RESIDUE:  return "amino_acid, nucleotide, water, ion, lipid and other";
    case VIAL_LEVEL_PART:     return "polymer, peptide, nucleic, dna, rna, glycan, ligand, lipid, ion and water";
    default: return "";
    }
}

bool vial_glob(str_t p, str_t t) {
    // Iterative, with one backtrack point for the last *: linear in practice
    size_t i = 0, j = 0, star = SIZE_MAX, mark = 0;
    while (j < t.len) {
        if (i < p.len && (p.ptr[i] == '?' || p.ptr[i] == t.ptr[j])) { ++i; ++j; }
        else if (i < p.len && p.ptr[i] == '*') { star = i++; mark = j; }
        else if (star != SIZE_MAX) { i = star + 1; j = ++mark; }
        else return false;
    }
    while (i < p.len && p.ptr[i] == '*') ++i;
    return i == p.len;
}

static bool match_str(const vial_filter_t* f, str_t s, bool glob) {
    for (size_t k = 0; k < md_array_size(f->any); ++k) {
        const vial_match_t* m = &f->any[k];
        if (m->kind != 0) continue;
        if (glob ? vial_glob(m->s, s) : str_eq(m->s, s)) return true;
    }
    return false;
}

// Element symbols ignore case: MG, Mg and mg are magnesium, and no two symbols differ only in case
static bool match_element(const vial_filter_t* f, str_t s) {
    for (size_t k = 0; k < md_array_size(f->any); ++k) {
        const vial_match_t* m = &f->any[k];
        if (m->kind == 0 && str_eq_ignore_case(m->s, s)) return true;
    }
    return false;
}

static bool match_int(const vial_filter_t* f, int64_t v) {
    for (size_t k = 0; k < md_array_size(f->any); ++k) {
        const vial_match_t* m = &f->any[k];
        if (m->kind == 1 && m->lo == v) return true;
        if (m->kind == 2 && m->lo <= v && v <= m->hi) return true;
    }
    return false;
}

static bool match_part_kind(const vial_filter_t* f, uint8_t k) {
    for (size_t q = 0; q < md_array_size(f->any); ++q) {
        const vial_match_t* m = &f->any[q];
        if (m->kind != 0) continue;
        const str_t w = m->s;
        if (str_eq(w, STR_LIT("polymer")) && k >= VIAL_PART_PEPTIDE && k <= VIAL_PART_POLYMER) return true;
        if (str_eq(w, STR_LIT("nucleic")) && (k == VIAL_PART_DNA || k == VIAL_PART_RNA || k == VIAL_PART_NUCLEIC)) return true;
        if (str_eq(w, str_from_cstr(part_kind_name(k))) && !str_eq(w, STR_LIT("nucleic"))) return true;
    }
    return false;
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

// The particle-level filters, name, element and kind, are properties of a particle's type: decided once per type,
// so a particle costs one lookup. NULL when no filter looks at the type.
static uint8_t* type_mask(const vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    bool any = q->particle_kind != 0;
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const uint8_t k = q->filters[i].key;
        any = any || k == VIAL_KEY_NAME || k == VIAL_KEY_ELEMENT || k == VIAL_KEY_KIND;
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
            case VIAL_KEY_NAME: pass = match_str(f, LBL_TO_STR(sys->atom.type.name[t]), true); break;
            case VIAL_KEY_ELEMENT: {
                const str_t sym = sys->atom.type.z ? md_util_element_symbol(sys->atom.type.z[t]) : STR_LIT("");
                pass = match_element(f, sym);
                break;
            }
            case VIAL_KEY_KIND: pass = kind < ARRAY_SIZE(particle_kinds) && match_str(f, str_from_cstr(particle_kinds[kind]), false); break;
            default: break;
            }
        }
        ok[t] = pass;
    }
    return ok;
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
            if (k >= ARRAY_SIZE(residue_kinds) || !match_str(f, str_from_cstr(residue_kinds[k]), false)) return false;
            break;
        }
        case VIAL_KEY_CHAIN: if (!chain_matches(s, f, vial_sys_res_chain(s, r))) return false; break;
        default: break;
        }
    }
    return true;
}

static bool part_passes(const vial_sys_t* s, const vial_select_t* q, uint32_t p) {
    for (size_t i = 0; i < md_array_size(q->filters); ++i) {
        const vial_filter_t* f = &q->filters[i];
        switch (f->key) {
        case VIAL_KEY_KIND:  if (!match_part_kind(f, s->part_kind[p])) return false; break;
        case VIAL_KEY_CHAIN: if (!chain_matches(s, f, s->part_chain[p])) return false; break;
        default: break;
        }
    }
    return true;
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

// particles(s, ...) of one set: only s's particles are looked at, so selecting within a residue costs its size, not
// the system's. Residues are found by walking s, which is ascending.
static vial_sets_t* select_particles_of(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    const md_system_t* sys = s->sys;
    const uint32_t* items = vial_sets_items(q->within, 0);
    const uint32_t n = vial_sets_size(q->within, 0);
    const uint8_t* types = type_mask(s, q, alloc);
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
            if (t < 0 || (uint32_t)t >= s->num_types || !types[t]) continue;
        }
        out->idx[w++] = a;
    }
    out->offset[1] = w;
    return out;
}

vial_sets_t* vial_select(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc) {
    if (q->level == VIAL_LEVEL_PARTICLE && q->within && q->within->single) return select_particles_of(s, q, alloc);
    const uint64_t* scope = q->scope ? bitmap(q->scope, s->num_atoms, alloc) : NULL;
    const uint64_t* within = q->within ? bitmap(q->within, s->num_atoms, alloc) : NULL;
    #define IN_SCOPE(a) (!scope || bit(scope, (a)))

    if (q->level == VIAL_LEVEL_PARTICLE) {
        // Residue by residue: a residue that fails resid, resname or chain is skipped whole, and a particle costs a
        // lookup of its type's mask
        const md_system_t* sys = s->sys;
        const uint8_t* types = type_mask(s, q, alloc);
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
                    if (t < 0 || (uint32_t)t >= s->num_types || !types[t]) continue;
                }
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
    default: break;
    }
    // The units that pass, and room for their particles: the units of a level never share a particle
    uint32_t* passing = md_alloc(alloc, (num_units + 1) * sizeof(uint32_t));
    uint32_t num_pass = 0, room = 0;
    for (uint32_t u = 0; u < num_units; ++u) {
        bool pass = false;
        uint32_t size = 0;
        switch (q->level) {
        case VIAL_LEVEL_RESIDUE:  pass = residue_passes(s, q, u);  size = (uint32_t)md_system_residue_atom_count(s->sys, u); break;
        case VIAL_LEVEL_PART:     pass = part_passes(s, q, u);     { const md_urange_t r = md_system_part_atom_range(s->sys, u); size = r.end - r.beg; } break;
        case VIAL_LEVEL_MOLECULE: pass = molecule_passes(s, q, u); size = s->mol_offset[u + 1] - s->mol_offset[u]; break;
        case VIAL_LEVEL_CHAIN:    pass = chain_passes(s, q, u);    size = s->chain_offset[u + 1] - s->chain_offset[u]; break;
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

vial_sets_t* vial_select_predefined(vial_sys_t* s, str_t name, md_allocator_i* alloc) {
    if (str_eq(name, STR_LIT("system"))) return vial_sets_all(s->num_atoms, alloc);
    const md_system_t* sys = s->sys;
    const bool backbone = str_eq(name, STR_LIT("backbone")), sidechain = str_eq(name, STR_LIT("sidechain"));
    int part = -1;
    if (str_eq(name, STR_LIT("protein"))) part = 0;
    else if (str_eq(name, STR_LIT("nucleic"))) part = 1;
    else if (str_eq(name, STR_LIT("water"))) part = 2;
    else if (str_eq(name, STR_LIT("ion"))) part = 3;
    else if (str_eq(name, STR_LIT("lipid"))) part = 4;
    if (!backbone && !sidechain && part < 0) return NULL;
    vial_sets_t* out = sets_new(1, s->num_atoms, alloc);
    out->single = true;
    uint32_t n = 0;
    if (backbone || sidechain) {
        const md_atom_flags_t want = backbone ? MD_ATOM_FLAG_BACKBONE : MD_ATOM_FLAG_SIDE_CHAIN;
        for (uint32_t a = 0; a < s->num_atoms; ++a) if (md_system_atom_flags(sys, a) & want) out->idx[n++] = a;
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

#include <md_chem.h>

#include <md_system.h>
#include <md_util.h>
#include <md_attributes.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ### ELEMENTS ###

// Valence of the neutral atom in its lowest common state, 0 for elements outside of what perception handles
// (metals, noble gases, unknown): those take no part.
static inline int base_valence(md_atomic_number_t z) {
    switch (z) {
    case MD_Z_H:  return 1;
    case MD_Z_B:  return 3;
    case MD_Z_C:  return 4;
    case MD_Z_N:  return 3;
    case MD_Z_O:  return 2;
    case MD_Z_F:  return 1;
    case MD_Z_Si: return 4;
    case MD_Z_P:  return 3;
    case MD_Z_S:  return 2;
    case MD_Z_Cl: return 1;
    case 33:      return 3; // As
    case 34:      return 2; // Se
    case MD_Z_Br: return 1;
    case 52:      return 2; // Te
    case MD_Z_I:  return 1;
    default:      return 0;
    }
}

static inline bool is_chalcogen(md_atomic_number_t z) { return z == MD_Z_O || z == MD_Z_S || z == 34 || z == 52; }
static inline bool is_heavy_chalcogen(md_atomic_number_t z) { return z == MD_Z_S || z == 34 || z == 52; }
static inline bool is_pnictogen(md_atomic_number_t z) { return z == MD_Z_N || z == MD_Z_P || z == 33; }
static inline bool is_halogen_z(md_atomic_number_t z) { return z == MD_Z_F || z == MD_Z_Cl || z == MD_Z_Br || z == MD_Z_I; }

// Charge of a lone ion of the element, as found in biomolecular structures and simulations
static int monatomic_ion_charge(md_atomic_number_t z) {
    switch (z) {
    case 3: case 11: case 19: case 37: case 55: case 47:    return 1;   // Li Na K Rb Cs Ag
    case 4: case 12: case 20: case 38: case 56:             return 2;   // Be Mg Ca Sr Ba
    case 25: case 26: case 27: case 28: case 29: case 30:   return 2;   // Mn Fe Co Ni Cu Zn
    case 48: case 80: case 82:                              return 2;   // Cd Hg Pb
    case 13:                                                return 3;   // Al
    case MD_Z_F: case MD_Z_Cl: case MD_Z_Br: case MD_Z_I:   return -1;
    default:                                                return 0;
    }
}

enum {
    HMODE_ALL   = 0,    // Hydrogens on carbons present: all hydrogens explicit
    HMODE_POLAR = 1,    // Hydrogens on N, O, S only: carbons implicit (united atom)
    HMODE_NONE  = 2,    // No hydrogens
};

typedef struct chem_t {
    md_system_t* sys;
    const vec3_t* xyz;              // NULL without coordinates
    const md_unitcell_t* cell;
    size_t N;
    uint32_t flags;
    md_allocator_i* alloc;

    // The covalent graph between atoms that take part (no metals, no virtual sites, no coordination bonds), CSR
    uint32_t* off;
    uint32_t* nbr;
    uint32_t* nbr_bond;

    uint8_t* z;
    uint8_t* part;                  // Takes part
    uint8_t* mode;
    uint8_t* heavy;                 // Heavy neighbours
    uint8_t* h_exp;
    uint8_t* h_imp;
    int8_t*  charge;
    uint8_t* charge_fixed;
    uint8_t* metal_bound;           // Coordinates a metal: without hydrogens the metal takes the place of a proton
    uint8_t* hyb;                   // From the geometry: 1 sp, 2 sp2, 3 sp3, 0 unknown
    uint8_t* in_ring;
    uint8_t* ring_bond;             // [bond] The bond lies in a ring
    uint8_t* cap;                   // Pi bonds the atom wants, beyond the fixed ones
    uint8_t* left;                  // Pi bonds the atom still wants after matching

    size_t   num_bonds;
    uint8_t* bond_part;
    uint8_t* pi_fixed;              // Extra order given by the file, or settled before the matching (triple bonds)
    uint8_t* given;                 // The order is the file's: the matching leaves the bond alone
    uint8_t* pi_match;              // Extra order from the matching
    float*   ratio;                 // Length relative to the sum of the covalent radii, 1 without coordinates
    uint8_t* aromatic;
    uint8_t* delocalized;

    uint32_t* mark;                 // Visited, when equal to mark_gen
    uint32_t  mark_gen;
} chem_t;

// Starts a new set of marks, without clearing them all
static inline void mark_begin(chem_t* c) {
    c->mark_gen += 1;
    if (c->mark_gen == 0) {
        MEMSET(c->mark, 0, sizeof(uint32_t) * c->N);
        c->mark_gen = 1;
    }
}

// Built with fast math: NAN is told by its bits, not by v != v
static inline bool float_bits_nan(float v) {
    uint32_t u;
    MEMCPY(&u, &v, sizeof(u));
    return (u & 0x7F800000u) == 0x7F800000u && (u & 0x007FFFFFu) != 0;
}

static inline int sigma(const chem_t* c, uint32_t i) {
    return (int)c->heavy[i] + (int)c->h_exp[i] + (int)c->h_imp[i];
}

static inline int pi_of_bond(const chem_t* c, uint32_t b) {
    return (int)c->pi_fixed[b] + (int)c->pi_match[b];
}

static inline int pi_sum(const chem_t* c, uint32_t i) {
    int s = 0;
    for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) s += pi_of_bond(c, c->nbr_bond[k]);
    return s;
}

static inline vec3_t bond_vec(const chem_t* c, uint32_t from, uint32_t to) {
    vec3_t d = vec3_sub(c->xyz[to], c->xyz[from]);
    // A bond is far shorter than 3 Å: only one split by the periodic boundary pays for the minimum image
    if (vec3_dot(d, d) > 9.0f) md_util_min_image_vec3(&d, 1, c->cell);
    return d;
}

// ### GEOMETRY ###

// Hybridization from the geometry around a heavy atom: the angles between its bonds, and for a terminal atom the
// length of its bond. 1 sp, 2 sp2, 3 sp3, 0 unknown.
static uint8_t geometric_hybridization(const chem_t* c, uint32_t i) {
    if (!c->xyz) return 0;
    vec3_t v[4];
    uint32_t nb[4];
    int n = 0;
    for (uint32_t k = c->off[i]; k < c->off[i + 1] && n < 4; ++k) {
        nb[n] = c->nbr[k];
        v[n] = bond_vec(c, i, c->nbr[k]);
        n += 1;
    }
    const int total = (int)(c->off[i + 1] - c->off[i]) + (int)c->h_imp[i];
    if (total >= 4) return 3;
    if (n == 3) {
        const float sum = (float)RAD_TO_DEG(vec3_angle(v[0], v[1]) + vec3_angle(v[1], v[2]) + vec3_angle(v[0], v[2]));
        return sum >= 350.0f ? 2 : 3;
    }
    if (n == 2) {
        const float a = (float)RAD_TO_DEG(vec3_angle(v[0], v[1]));
        if (a > 165.0f) return 1;
        // sp2 has a short bond (a double or aromatic one), whatever the angle: about 108 in a five membered ring.
        // Without one only a wide angle says sp2: a CH2 in a ring is often opened to 117.
        float rmin = FLT_MAX;
        for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
            if (c->z[c->nbr[k]] != MD_Z_H) rmin = MIN(rmin, c->ratio[c->nbr_bond[k]]);
        }
        if (rmin < 0.83f) return 1;     // A triple bond, also bent (in a strained ring)
        if (rmin < 0.94f) return 2;
        return a > 128.0f ? 2 : 3;
    }
    if (n == 1) {
        if (c->z[nb[0]] == MD_Z_H) return 3;
        const float r = c->ratio[c->nbr_bond[c->off[i]]];
        if (r < 0.83f) return 1;
        if (r < 0.935f) return 2;
        return 3;
    }
    return 3;
}

static float planar_angle_sum(const chem_t* c, uint32_t i) {
    if (!c->xyz || c->off[i + 1] - c->off[i] != 3) return 0.0f;
    vec3_t v[3];
    for (int k = 0; k < 3; ++k) v[k] = bond_vec(c, i, c->nbr[c->off[i] + k]);
    return (float)RAD_TO_DEG(vec3_angle(v[0], v[1]) + vec3_angle(v[1], v[2]) + vec3_angle(v[0], v[2]));
}

// ### SETUP ###

static void build_graph(chem_t* c) {
    md_system_t* sys = c->sys;
    const size_t N = c->N;
    c->z     = md_alloc(c->alloc, N);
    c->part  = md_alloc(c->alloc, N);
    for (size_t i = 0; i < N; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
        const md_flags_t f = md_atom_flags(&sys->atom, i) | md_atom_type_flags(&sys->atom.type, sys->atom.type_idx[i]);
        c->z[i] = (uint8_t)z;
        c->part[i] = base_valence(z) > 0 && !(f & (MD_FLAG_VIRTUAL_SITE | MD_FLAG_COARSE_GRAINED));
    }

    c->num_bonds = sys->bond.count;
    c->metal_bound = md_alloc(c->alloc, MAX(N, 1));
    MEMSET(c->metal_bound, 0, MAX(N, 1));
    for (size_t b = 0; b < c->num_bonds; ++b) {
        const md_atom_pair_t p = sys->bond.pairs[b];
        if (p.idx[0] < 0 || p.idx[1] < 0 || (size_t)p.idx[0] >= N || (size_t)p.idx[1] >= N) continue;
        const bool m0 = !c->part[p.idx[0]] && c->z[p.idx[0]] > 2;
        const bool m1 = !c->part[p.idx[1]] && c->z[p.idx[1]] > 2;
        if (m0 && c->part[p.idx[1]]) c->metal_bound[p.idx[1]] = 1;
        if (m1 && c->part[p.idx[0]]) c->metal_bound[p.idx[0]] = 1;
    }
    c->bond_part = md_alloc(c->alloc, MAX(c->num_bonds, 1));
    uint32_t* deg = md_alloc(c->alloc, sizeof(uint32_t) * (N + 1));
    MEMSET(deg, 0, sizeof(uint32_t) * (N + 1));
    for (size_t b = 0; b < c->num_bonds; ++b) {
        const md_atom_pair_t p = sys->bond.pairs[b];
        const md_bond_flags_t f = sys->bond.flags[b];
        const bool ok = p.idx[0] >= 0 && p.idx[1] >= 0 && (size_t)p.idx[0] < N && (size_t)p.idx[1] < N && p.idx[0] != p.idx[1] &&
                        c->part[p.idx[0]] && c->part[p.idx[1]] && !(f & MD_BOND_FLAG_COORDINATE);
        c->bond_part[b] = ok;
        if (ok) {
            deg[p.idx[0]] += 1;
            deg[p.idx[1]] += 1;
        }
    }
    c->off = md_alloc(c->alloc, sizeof(uint32_t) * (N + 1));
    c->off[0] = 0;
    for (size_t i = 0; i < N; ++i) c->off[i + 1] = c->off[i] + deg[i];
    const uint32_t total = c->off[N];
    c->nbr      = md_alloc(c->alloc, sizeof(uint32_t) * MAX(total, 1));
    c->nbr_bond = md_alloc(c->alloc, sizeof(uint32_t) * MAX(total, 1));
    MEMSET(deg, 0, sizeof(uint32_t) * (N + 1));
    for (size_t b = 0; b < c->num_bonds; ++b) {
        if (!c->bond_part[b]) continue;
        const uint32_t x = (uint32_t)sys->bond.pairs[b].idx[0];
        const uint32_t y = (uint32_t)sys->bond.pairs[b].idx[1];
        c->nbr[c->off[x] + deg[x]] = y; c->nbr_bond[c->off[x] + deg[x]] = (uint32_t)b; deg[x] += 1;
        c->nbr[c->off[y] + deg[y]] = x; c->nbr_bond[c->off[y] + deg[y]] = (uint32_t)b; deg[y] += 1;
    }

    c->ratio = md_alloc(c->alloc, sizeof(float) * MAX(c->num_bonds, 1));
    for (size_t b = 0; b < c->num_bonds; ++b) {
        c->ratio[b] = 1.0f;
        if (!c->bond_part[b] || !c->xyz) continue;
        // Bonds to hydrogens are never double, and their lengths say nothing the angles do not
        if (c->z[sys->bond.pairs[b].idx[0]] == MD_Z_H || c->z[sys->bond.pairs[b].idx[1]] == MD_Z_H) continue;
        const uint32_t x = (uint32_t)sys->bond.pairs[b].idx[0];
        const uint32_t y = (uint32_t)sys->bond.pairs[b].idx[1];
        const float sum = md_util_element_covalent_radius(c->z[x]) + md_util_element_covalent_radius(c->z[y]);
        if (sum > 0) c->ratio[b] = vec3_length(bond_vec(c, x, y)) / sum;
    }
}

// Hydrogen mode of each atom: by residue where there are residues, by molecule otherwise
static uint32_t uf_find(uint32_t* p, uint32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

static void assign_modes(chem_t* c) {
    md_system_t* sys = c->sys;
    const size_t N = c->N;
    c->mode = md_alloc(c->alloc, N);
    uint32_t* group = md_alloc(c->alloc, sizeof(uint32_t) * MAX(N, 1));

    if (sys->component.count && sys->component.atom_offset) {
        for (size_t i = 0; i < N; ++i) group[i] = UINT32_MAX;
        for (size_t ci = 0; ci < sys->component.count; ++ci) {
            const md_urange_t r = md_component_atom_range(&sys->component, ci);
            for (uint32_t i = r.beg; i < r.end && i < N; ++i) group[i] = (uint32_t)ci;
        }
        // Atoms outside of every residue form their own groups by molecule
        for (size_t i = 0; i < N; ++i) if (group[i] == UINT32_MAX) group[i] = (uint32_t)(sys->component.count + i);
    } else {
        for (size_t i = 0; i < N; ++i) group[i] = (uint32_t)i;
        for (size_t i = 0; i < N; ++i) {
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                const uint32_t a = uf_find(group, (uint32_t)i);
                const uint32_t b = uf_find(group, c->nbr[k]);
                if (a != b) group[MAX(a, b)] = MIN(a, b);
            }
        }
        for (size_t i = 0; i < N; ++i) group[i] = uf_find(group, (uint32_t)i);
    }

    // Per group: any hydrogen, any hydrogen on a carbon
    const size_t num_groups = sys->component.count + N;
    uint8_t* has = md_alloc(c->alloc, MAX(num_groups, 1));
    MEMSET(has, 0, num_groups);
    for (size_t i = 0; i < N; ++i) {
        if (c->z[i] != MD_Z_H) continue;
        has[group[i]] |= 1;
        for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
            if (c->z[c->nbr[k]] == MD_Z_C) has[group[i]] |= 2;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        const uint8_t h = has[group[i]];
        c->mode[i] = (h & 2) ? HMODE_ALL : (h & 1) ? HMODE_POLAR : HMODE_NONE;
    }
}

static int find_bond(const chem_t* c, uint32_t a, uint32_t b) {
    for (uint32_t k = c->off[a]; k < c->off[a + 1]; ++k) if (c->nbr[k] == b) return (int)c->nbr_bond[k];
    return -1;
}

static void mark_rings(chem_t* c) {
    c->in_ring = md_alloc(c->alloc, c->N);
    MEMSET(c->in_ring, 0, c->N);
    c->ring_bond = md_alloc(c->alloc, MAX(c->num_bonds, 1));
    MEMSET(c->ring_bond, 0, MAX(c->num_bonds, 1));
    const md_index_data_t* rings = &c->sys->ring;
    const size_t num = md_index_data_num_ranges(rings);
    for (size_t r = 0; r < num; ++r) {
        const int32_t* ring = md_index_range_beg(rings, r);
        const size_t n = md_index_range_size(rings, r);
        bool valid = true;
        for (size_t k = 0; k < n; ++k) valid &= ring[k] >= 0 && (size_t)ring[k] < c->N;
        if (!valid) continue;
        for (size_t k = 0; k < n; ++k) {
            c->in_ring[ring[k]] = 1;
            const int b = find_bond(c, (uint32_t)ring[k], (uint32_t)ring[(k + 1) % n]);
            if (b >= 0) c->ring_bond[b] = 1;
        }
    }
}

// Valence the atom has to fill with its sigma bonds and pi bonds, given its sigma bonds and charge. Sets the charge
// of atoms whose sigma bonds alone call for one (ammonium N+, oxonium O+, borate B-).
static int target_valence(chem_t* c, uint32_t i) {
    const md_atomic_number_t z = c->z[i];
    const int s = sigma(c, i);
    const int base = base_valence(z);
    // Clusters (boranes, carboranes) have more sigma bonds than any classical valence: nothing to place
    if ((z == MD_Z_B || z == MD_Z_C) && s >= 5) return s;
    if (c->charge_fixed[i] || c->charge[i]) {
        const int q = c->charge[i];
        if (z == MD_Z_C || z == MD_Z_Si) return 4 - abs(q);
        if (z == MD_Z_B) return 3 - q;
        if (is_pnictogen(z) || is_chalcogen(z)) {
            // Hypervalent with a charge: keep the expanded valence
            if (z != MD_Z_N && z != MD_Z_O && s > base + q) return s + ((s - base - q) & 1);
            return base + q;
        }
        return base - q;
    }
    switch (z) {
    case MD_Z_N:
        if (s >= 4) { c->charge[i] = 1; return 4; }
        return 3;
    case MD_Z_O:
        if (s >= 3) { c->charge[i] = 1; return 3; }
        return 2;
    case MD_Z_B:
        if (s >= 4) { c->charge[i] = -1; return 4; }
        return 3;
    case MD_Z_P: case 33:
        if (s >= 4) return 5;   // Phosphate, phosphine oxide: P=O
        return 3;
    case MD_Z_S: case 34: case 52:
        if (s == 3) return 4;   // Sulfoxide
        if (s >= 4) return 6;   // Sulfone, sulfonamide, sulfate
        return 2;
    case MD_Z_Cl: case MD_Z_Br: case MD_Z_I:
        if (s >= 2) return 7;   // Perchlorate and kin
        return 1;
    default:
        return base;
    }
}

// Implicit hydrogens from the geometry and the valence. Carbons in residues without hydrogens on carbon, and every
// heavy atom in residues without any hydrogen. N, O and S start with those the geometry demands (an sp3 amine, an
// alcohol); whether the rest are protonated is settled after the matching.
static void implicit_hydrogens(chem_t* c) {
    for (size_t i = 0; i < c->N; ++i) {
        c->h_imp[i] = 0;
        if (!c->part[i] || c->z[i] == MD_Z_H || is_halogen_z(c->z[i])) continue;
        const int mode = c->mode[i];
        if (mode == HMODE_ALL) continue;
        const md_atomic_number_t z = c->z[i];
        if (mode == HMODE_POLAR && z != MD_Z_C && z != MD_Z_B && z != MD_Z_Si) continue;

        const int n = (int)c->heavy[i] + (int)c->h_exp[i];
        const int hyb = c->hyb[i] ? c->hyb[i] : 3;
        const int pi = 3 - hyb;
        if (c->charge_fixed[i] && mode == HMODE_NONE) {
            // A charge the file gives sets the valence; the geometry says how much of it is pi bonds (a flat N+
            // has one: the protonated N of a cationic histidine, a pyridinium)
            const int q = c->charge[i];
            int v = base_valence(z);
            if (z == MD_Z_C || z == MD_Z_Si) v = 4 - abs(q);
            else if (z == MD_Z_B) v = 3 - q;
            else if (is_pnictogen(z) || is_chalcogen(z)) v = base_valence(z) + q;
            c->h_imp[i] = (uint8_t)MAX(0, v - n - pi);
            continue;
        }
        if (z == MD_Z_C || z == MD_Z_Si) {
            c->h_imp[i] = (uint8_t)MAX(0, 4 - n - pi);
        } else if (z == MD_Z_B) {
            c->h_imp[i] = (uint8_t)MAX(0, 3 - n);
        } else if (z == MD_Z_N) {
            // An sp3 amine. A flat terminal N has at least one (=NH, or -NH2 conjugated: the matching decides), a
            // linear one none (nitrile); anything else flat waits for the matching.
            if (hyb == 3 && c->hyb[i]) c->h_imp[i] = (uint8_t)MAX(0, 3 - n);
            else if (n == 1 && c->hyb[i] == 2) c->h_imp[i] = 1;
            if (n == 0) c->h_imp[i] = 3;
        } else if (z == MD_Z_O) {
            if (n == 0) c->h_imp[i] = 2;                    // Water
            else if (n == 1 && c->hyb[i]) {
                // A terminal oxygen on an atom which cannot take a double bond (an sp3 carbon), or bonded to a
                // carbon by a single bond (a phenol, enol, the OH of an acid), is a hydroxyl: the matching would
                // otherwise rather make quinones. The oxygens of P, S, As and Cl are left to the matching (P=O).
                uint32_t nb = UINT32_MAX;
                float r = 0;
                for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                    if (c->z[c->nbr[k]] != MD_Z_H) { nb = c->nbr[k]; r = c->ratio[c->nbr_bond[k]]; }
                }
                const md_atomic_number_t zn = nb != UINT32_MAX ? c->z[nb] : 0;
                if (zn == MD_Z_C || zn == MD_Z_Si || zn == MD_Z_B) {
                    const int nn = (int)c->heavy[nb] + (int)c->h_exp[nb];
                    // The C-O of a phenol is 1.36 A (0.96 of the radii), of an acid 1.31 (0.92), a C=O 1.23 (0.87).
                    // Acids (another terminal O or S on the carbon) are left to the matching and to the pH.
                    int terminal = 0;
                    for (uint32_t k = c->off[nb]; k < c->off[nb + 1]; ++k) terminal += is_chalcogen(c->z[c->nbr[k]]) && c->heavy[c->nbr[k]] == 1;
                    if ((hyb == 3 && nn >= 4) || (r >= 0.945f && terminal == 1)) c->h_imp[i] = 1;
                }
            }
        } else if (is_heavy_chalcogen(z)) {
            if (n == 0) c->h_imp[i] = 2;
        }
    }
}

// N-oxides and nitro groups in their charge separated form: an N with three sigma bonds and a terminal O without
// hydrogen is N+ and the O, of the longest such bond, O-. The N then takes its pi bond where it is conjugated
// (pyridine N-oxide, nitrone, the other O of a nitro group). Without hydrogens the O must be short enough not to be
// an N-OH (hydroxylamine, hydroxamic acid).
static void charge_n_oxides(chem_t* c) {
    for (size_t i = 0; i < c->N; ++i) {
        if (!c->part[i] || c->z[i] != MD_Z_N || sigma(c, (uint32_t)i) != 3 || c->charge[i] || c->charge_fixed[i]) continue;
        int best = -1;
        float best_ratio = 0;
        for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
            const uint32_t o = c->nbr[k];
            if (c->z[o] != MD_Z_O || c->heavy[o] != 1 || c->h_exp[o] + c->h_imp[o] != 0 || c->charge[o] || c->charge_fixed[o]) continue;
            const float r = c->ratio[c->nbr_bond[k]];
            if (c->mode[o] == HMODE_NONE && c->xyz && r > 0.98f) continue;
            if (best < 0 || r > best_ratio) {
                best = (int)o;
                best_ratio = r;
            }
        }
        if (best >= 0) {
            c->charge[i] = 1;
            c->charge[best] = -1;
        }
    }
}

// ### MATCHING ###
// Maximum matching in a general graph (Edmonds), over vertices which stand for the pi bonds an atom wants: an atom
// wanting two has two vertices. The graphs are the conjugated parts of molecules, small.

typedef struct matcher_t {
    int n;
    int* adj_off;
    int* adj;
    int* match;
    int* p;
    int* base;
    int* q;
    uint8_t* used;
    uint8_t* blossom;
    uint8_t* lca_used;
} matcher_t;

static int mt_lca(matcher_t* m, int a, int b) {
    MEMSET(m->lca_used, 0, (size_t)m->n);
    for (;;) {
        a = m->base[a];
        m->lca_used[a] = 1;
        if (m->match[a] == -1) break;
        a = m->p[m->match[a]];
    }
    for (;;) {
        b = m->base[b];
        if (m->lca_used[b]) return b;
        b = m->p[m->match[b]];
    }
}

static void mt_mark_path(matcher_t* m, int v, int b, int children) {
    while (m->base[v] != b) {
        m->blossom[m->base[v]] = m->blossom[m->base[m->match[v]]] = 1;
        m->p[v] = children;
        children = m->match[v];
        v = m->p[m->match[v]];
    }
}

static int mt_find_path(matcher_t* m, int root) {
    const int n = m->n;
    MEMSET(m->used, 0, (size_t)n);
    for (int i = 0; i < n; ++i) { m->p[i] = -1; m->base[i] = i; }
    m->used[root] = 1;
    int qh = 0, qt = 0;
    m->q[qt++] = root;
    while (qh < qt) {
        const int v = m->q[qh++];
        for (int k = m->adj_off[v]; k < m->adj_off[v + 1]; ++k) {
            int to = m->adj[k];
            if (m->base[v] == m->base[to] || m->match[v] == to) continue;
            if (to == root || (m->match[to] != -1 && m->p[m->match[to]] != -1)) {
                const int curbase = mt_lca(m, v, to);
                MEMSET(m->blossom, 0, (size_t)n);
                mt_mark_path(m, v, curbase, to);
                mt_mark_path(m, to, curbase, v);
                for (int i = 0; i < n; ++i) {
                    if (m->blossom[m->base[i]]) {
                        m->base[i] = curbase;
                        if (!m->used[i]) {
                            m->used[i] = 1;
                            m->q[qt++] = i;
                        }
                    }
                }
            } else if (m->p[to] == -1) {
                m->p[to] = v;
                if (m->match[to] == -1) return to;
                to = m->match[to];
                m->used[to] = 1;
                m->q[qt++] = to;
            }
        }
    }
    return -1;
}

typedef struct vedge_t {
    float key;
    int a, b;
} vedge_t;

static int compare_vedge(const void* x, const void* y) {
    const vedge_t* a = (const vedge_t*)x;
    const vedge_t* b = (const vedge_t*)y;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    return (a->b > b->b) - (a->b < b->b);
}

// Collects the atoms that want pi bonds and are connected to 'seed' through bonds between such atoms
static size_t collect_component(chem_t* c, uint32_t seed, uint32_t* out) {
    size_t n = 0;
    if (!c->cap[seed] || c->mark[seed] == c->mark_gen) return 0;
    out[n++] = seed;
    c->mark[seed] = c->mark_gen;
    for (size_t h = 0; h < n; ++h) {
        const uint32_t a = out[h];
        for (uint32_t k = c->off[a]; k < c->off[a + 1]; ++k) {
            const uint32_t b = c->nbr[k];
            if (c->cap[b] && c->mark[b] != c->mark_gen) {
                c->mark[b] = c->mark_gen;
                out[n++] = b;
            }
        }
    }
    return n;
}

enum { RELOCATE_NONE = 0, RELOCATE_GEOMETRY = 1, RELOCATE_SCORE = 2 };

// How well an atom resolves a pi bond it cannot have: a terminal O or S as an anion (or protonated), an N charged or
// protonated, a P or hypervalent S as a cation. A carbon left wanting is a radical or carbocation.
static int free_score(const chem_t* c, uint32_t a) {
    const md_atomic_number_t z = c->z[a];
    // Without hydrogens what is left is protonated, and which atom is a tautomer. Carbons are avoided, and an N is
    // protonated before an O or S: amide, lactam, thioamide and thiourea over imidic acid and isothiourea; a
    // terminal N most of all: amino over imino (aminopyridine, adenine, cytosine). Otherwise the geometry decides.
    if (c->mode[a] == HMODE_NONE) {
        if (z == MD_Z_C || z == MD_Z_Si || z == MD_Z_B) return 0;
        if (z == MD_Z_N) return c->heavy[a] <= 1 ? 3 : 2;
        return 1;
    }
    if (is_chalcogen(z)) return c->heavy[a] <= 1 ? 3 : 1;
    if (z == MD_Z_N) return 2;
    if (z == MD_Z_P || z == 33) return 1;
    return 0;
}

// Matches the atoms of one component; returns the number of pi bonds left wanting
static int match_component(chem_t* c, const uint32_t* atoms, size_t num, int32_t* vbase, int relocate) {
    // Two atoms wanting one pi bond each (a carbonyl, most of them): nothing to choose
    if (num == 2 && c->cap[atoms[0]] == 1 && c->cap[atoms[1]] == 1) {
        const uint32_t a = atoms[0], b = atoms[1];
        int bond = -1;
        for (uint32_t k = c->off[a]; k < c->off[a + 1]; ++k) {
            c->pi_match[c->nbr_bond[k]] = 0;
            if (c->nbr[k] == b) bond = (int)c->nbr_bond[k];
        }
        for (uint32_t k = c->off[b]; k < c->off[b + 1]; ++k) c->pi_match[c->nbr_bond[k]] = 0;
        const bool too_long = bond >= 0 && ((c->xyz && c->ratio[bond] >= 0.95f && (c->mode[a] == HMODE_NONE || c->mode[b] == HMODE_NONE)) || c->given[bond]);
        if (bond >= 0 && !too_long) {
            c->pi_match[bond] = 1;
            c->left[a] = c->left[b] = 0;
            return 0;
        }
        c->left[a] = c->left[b] = 1;
        return 2;
    }

    md_allocator_i* arena = c->alloc;
    const size_t arena_pos = md_vm_arena_get_pos(arena);

    // Vertices
    int nv = 0;
    for (size_t k = 0; k < num; ++k) {
        vbase[atoms[k]] = nv;
        nv += c->cap[atoms[k]];
        c->left[atoms[k]] = c->cap[atoms[k]];
    }
    // Reset the matching of the component's bonds
    for (size_t k = 0; k < num; ++k) {
        const uint32_t a = atoms[k];
        for (uint32_t j = c->off[a]; j < c->off[a + 1]; ++j) c->pi_match[c->nbr_bond[j]] = 0;
    }

    md_array(vedge_t) edges = 0;
    for (size_t k = 0; k < num; ++k) {
        const uint32_t a = atoms[k];
        for (uint32_t j = c->off[a]; j < c->off[a + 1]; ++j) {
            const uint32_t b = c->nbr[j];
            if (b <= a || !c->cap[b] || vbase[b] < 0 || c->given[c->nbr_bond[j]]) continue;
            const float key = c->ratio[c->nbr_bond[j]];
            // Without hydrogens the lengths are the evidence: a bond as long as a single one is not made double (the
            // N-N of a hydrazide, the link of a biaryl)
            if (c->xyz && key >= 0.95f && (c->mode[a] == HMODE_NONE || c->mode[b] == HMODE_NONE)) continue;
            for (int x = 0; x < c->cap[a]; ++x) {
                for (int y = 0; y < c->cap[b]; ++y) {
                    vedge_t e = { key, vbase[a] + x, vbase[b] + y };
                    md_array_push(edges, e, arena);
                }
            }
        }
    }
    const size_t ne = md_array_size(edges);

    matcher_t m = { .n = nv };
    m.adj_off = md_alloc(arena, sizeof(int) * (nv + 1));
    m.adj     = md_alloc(arena, sizeof(int) * MAX(2 * ne, 1));
    float* adj_key = md_alloc(arena, sizeof(float) * MAX(2 * ne, 1));
    m.match   = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    m.p       = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    m.base    = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    m.q       = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    m.used    = md_alloc(arena, MAX(nv, 1));
    m.blossom = md_alloc(arena, MAX(nv, 1));
    m.lca_used = md_alloc(arena, MAX(nv, 1));

    if (ne) qsort(edges, ne, sizeof(vedge_t), compare_vedge);
    int* deg = md_alloc(arena, sizeof(int) * (nv + 1));
    MEMSET(deg, 0, sizeof(int) * (nv + 1));
    for (size_t e = 0; e < ne; ++e) { deg[edges[e].a] += 1; deg[edges[e].b] += 1; }
    m.adj_off[0] = 0;
    for (int v = 0; v < nv; ++v) m.adj_off[v + 1] = m.adj_off[v] + deg[v];
    MEMSET(deg, 0, sizeof(int) * (nv + 1));
    for (size_t e = 0; e < ne; ++e) {
        const int a = edges[e].a, b = edges[e].b;
        adj_key[m.adj_off[a] + deg[a]] = edges[e].key;
        adj_key[m.adj_off[b] + deg[b]] = edges[e].key;
        m.adj[m.adj_off[a] + deg[a]++] = b;
        m.adj[m.adj_off[b] + deg[b]++] = a;
    }

    // Greedy start: the shortest bonds (relative to single) first
    for (int v = 0; v < nv; ++v) m.match[v] = -1;
    for (size_t e = 0; e < ne; ++e) {
        const int a = edges[e].a, b = edges[e].b;
        if (m.match[a] == -1 && m.match[b] == -1) {
            m.match[a] = b;
            m.match[b] = a;
        }
    }
    // Vertex back to atom
    uint32_t* vatom = md_alloc(arena, sizeof(uint32_t) * MAX(nv, 1));
    for (size_t k = 0; k < num; ++k) {
        for (int x = 0; x < c->cap[atoms[k]]; ++x) vatom[vbase[atoms[k]] + x] = atoms[k];
    }

    // Augment to a maximum matching. A vertex once matched stays matched, so the roots are taken carbons first and
    // the atoms which resolve a missing pi bond best (free_score) last: those are the ones left over.
    for (int level = 0; level <= 3; ++level) {
        for (int v = 0; v < nv; ++v) {
            if (m.match[v] != -1 || free_score(c, vatom[v]) != level) continue;
            int u = mt_find_path(&m, v);
            while (u != -1) {
                const int pv = m.p[u];
                const int ppv = m.match[pv];
                m.match[u] = pv;
                m.match[pv] = u;
                u = ppv;
            }
        }
    }


    // A maximum matching may leave any of several vertices unmatched, and the augmentation does not look at the
    // bond lengths. Unmatched vertices are moved along even alternating paths u - x1 = y1 - x2 = y2 ... yk, which
    // match u and free yk: to an atom that resolves the missing pi bond better (RELOCATE_SCORE: a terminal O becomes
    // an anion, an N is protonated or charged; a carbon is avoided), or, as well, where the double bonds come out
    // shorter (RELOCATE_GEOMETRY: the C=O of an amide rather than a quinoid C=N).
    int*   via   = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    int*   from  = md_alloc(arena, sizeof(int) * MAX(nv, 1));
    float* delta = md_alloc(arena, sizeof(float) * MAX(nv, 1));
    for (int u = 0; u < nv && relocate; ++u) {
        if (m.match[u] != -1) continue;
        const int su = free_score(c, vatom[u]);
        MEMSET(m.used, 0, (size_t)nv);
        int qh = 0, qt = 0;
        m.q[qt++] = u;
        m.used[u] = 1;
        delta[u] = 0;
        int best = -1, best_gain = 0;
        float best_delta = -0.02f;      // Shorter by more than noise
        while (qh < qt) {
            const int y = m.q[qh++];
            for (int k = m.adj_off[y]; k < m.adj_off[y + 1]; ++k) {
                const int x = m.adj[k];
                if (m.used[x] || m.match[x] == -1 || m.match[y] == x) continue;
                const int y2 = m.match[x];
                if (m.used[y2]) continue;
                // Length of the matched edge x = y2
                float matched = 1.0f;
                for (int j = m.adj_off[x]; j < m.adj_off[x + 1]; ++j) if (m.adj[j] == y2) { matched = adj_key[j]; break; }
                m.used[x] = m.used[y2] = 1;
                from[x] = y;
                via[y2] = x;
                delta[y2] = delta[y] + adj_key[k] - matched;
                m.q[qt++] = y2;
                const int gain = free_score(c, vatom[y2]) - su;
                // Before the fixes a carbon keeps what it lacks (a cation nearby may complete it), except without
                // hydrogens, where the only fix is the protonation of a heteroatom anyway
                if (gain < 0 || (relocate == RELOCATE_GEOMETRY && gain > 0 && c->mode[vatom[u]] != HMODE_NONE)) continue;
                if (gain > best_gain || (gain == best_gain && delta[y2] < best_delta)) {
                    best = y2;
                    best_gain = gain;
                    best_delta = delta[y2];
                }
            }
        }
        if (best == -1) continue;
        int cur = best;
        m.match[cur] = -1;
        for (;;) {
            const int x = via[cur];
            const int y = from[x];
            m.match[x] = y;
            m.match[y] = x;
            if (y == u) break;
            cur = y;
        }
    }
    int left = 0;
    for (int v = 0; v < nv; ++v) {
        const int u = m.match[v];
        if (u == -1) { left += 1; continue; }
        if (u < v) continue;
        const uint32_t a = vatom[v], b = vatom[u];
        for (uint32_t j = c->off[a]; j < c->off[a + 1]; ++j) {
            if (c->nbr[j] == b) { c->pi_match[c->nbr_bond[j]] += 1; break; }
        }
        c->left[a] -= 1;
        c->left[b] -= 1;
    }
    for (size_t k = 0; k < num; ++k) vbase[atoms[k]] = -1;
    md_vm_arena_set_pos_back(arena, arena_pos);
    return left;
}

static int match_all(chem_t* c, uint32_t* buf, int32_t* vbase, int relocate) {
    mark_begin(c);
    int left = 0;
    for (size_t i = 0; i < c->N; ++i) {
        if (!c->cap[i] || c->mark[i] == c->mark_gen) continue;
        const size_t n = collect_component(c, (uint32_t)i, buf);
        if (relocate == RELOCATE_SCORE) {
            // Only what is left wanting is moved: components matched completely stay as they are
            int want = 0;
            for (size_t k = 0; k < n; ++k) want += c->left[buf[k]];
            if (!want) continue;
        }
        left += match_component(c, buf, n, vbase, relocate);
    }
    return left;
}

// An oxygen (or S) on the end of a group with another such oxygen: carboxylate, phosphate, sulfonate, nitro
static bool in_anionic_group(const chem_t* c, uint32_t o) {
    if (c->heavy[o] != 1 || c->h_exp[o] + c->h_imp[o] != 0) return false;
    uint32_t center = UINT32_MAX;
    for (uint32_t k = c->off[o]; k < c->off[o + 1]; ++k) if (c->z[c->nbr[k]] != MD_Z_H) center = c->nbr[k];
    if (center == UINT32_MAX) return false;
    const md_atomic_number_t zc = c->z[center];
    if (!(zc == MD_Z_C || zc == MD_Z_N || zc == MD_Z_P || zc == MD_Z_S || zc == 33 || zc == 34 || zc == MD_Z_Cl)) return false;
    int terminal = 0;
    for (uint32_t k = c->off[center]; k < c->off[center + 1]; ++k) {
        const uint32_t x = c->nbr[k];
        if (is_chalcogen(c->z[x]) && c->heavy[x] == 1 && c->h_exp[x] + c->h_imp[x] == 0) terminal += 1;
    }
    return terminal >= 2;
}

// The components a change at atom x touches: those of x, of its neighbours and of 'seed'. Returns the pi bonds they
// leave wanting, after matching them again if 'rematch'.
static int region_seed(chem_t* c, uint32_t s, bool rematch, uint32_t* buf, int32_t* vbase) {
    if (!c->cap[s] || c->mark[s] == c->mark_gen) return 0;
    const size_t n = collect_component(c, s, buf);
    if (rematch) return match_component(c, buf, n, vbase, RELOCATE_GEOMETRY);
    int left = 0;
    for (size_t k = 0; k < n; ++k) left += c->left[buf[k]];
    return left;
}

static int region(chem_t* c, uint32_t seed, uint32_t x, bool rematch, uint32_t* buf, int32_t* vbase) {
    if (rematch) {
        for (uint32_t k = c->off[x]; k < c->off[x + 1]; ++k) c->pi_match[c->nbr_bond[k]] = 0;
        if (!c->cap[x]) c->left[x] = 0;
    }
    mark_begin(c);
    int left = region_seed(c, seed, rematch, buf, vbase) + region_seed(c, x, rematch, buf, vbase);
    for (uint32_t k = c->off[x]; k < c->off[x + 1]; ++k) left += region_seed(c, c->nbr[k], rematch, buf, vbase);
    return left;
}

static bool expandable(const chem_t* c, uint32_t x) {
    const md_atomic_number_t z = c->z[x];
    const int s = sigma(c, x);
    return (z == MD_Z_N && s == 3) || ((z == MD_Z_O || is_heavy_chalcogen(z)) && s == 2 && c->in_ring[x]);
}

#define PROTONATE_BIT 0x80000000u
#define MAX_FIX_CANDIDATES 8

// An N next to a carbonyl (amide, lactam, urea) has its lone pair taken: a poor place for a positive charge
static bool next_to_carbonyl(const chem_t* c, uint32_t x) {
    for (uint32_t k = c->off[x]; k < c->off[x + 1]; ++k) {
        const uint32_t y = c->nbr[k];
        if (c->z[y] != MD_Z_C) continue;
        for (uint32_t j = c->off[y]; j < c->off[y + 1]; ++j) {
            const uint32_t o = c->nbr[j];
            if (is_chalcogen(c->z[o]) && c->heavy[o] == 1 && pi_of_bond(c, c->nbr_bond[j]) > 0) return true;
        }
    }
    return false;
}

// Preference among fixes that resolve as much: a cation on N before O or S (pyrylium, thiophenium), an N whose lone
// pair is free before an amide N, and for a leftover in a ring a cation in the ring (an aromatic N-alkyl pyridinium
// or adeninium, not an exocyclic iminium).
static int fix_penalty(const chem_t* c, uint32_t cand, uint32_t leftover) {
    if (cand & PROTONATE_BIT) return 0;
    const uint32_t x = cand;
    int p = 0;
    if (c->z[x] != MD_Z_N) p += 4;
    if (next_to_carbonyl(c, x)) p += 2;
    if (c->in_ring[leftover]) {
        // Joined to the conjugated atoms through a ring bond
        bool ring = false;
        for (uint32_t k = c->off[x]; k < c->off[x + 1]; ++k) ring |= c->cap[c->nbr[k]] && c->ring_bond[c->nbr_bond[k]];
        if (!ring) p += 1;
    }
    return p;
}

static void apply_fix(chem_t* c, uint32_t cand, bool undo) {
    const uint32_t x = cand & ~PROTONATE_BIT;
    if (cand & PROTONATE_BIT) {
        c->cap[x] = undo ? 1 : 0;
        c->h_imp[x] = (uint8_t)(c->h_imp[x] + (undo ? -1 : 1));
    } else {
        c->cap[x] = (uint8_t)(c->cap[x] + (undo ? -1 : 1));
        c->charge[x] = undo ? 0 : 1;
    }
}

// A linear N with two sigma bonds which can take a second pi bond as a cation: the middle N of an azide or diazo
static bool linear_cation(const chem_t* c, uint32_t a) {
    return c->z[a] == MD_Z_N && c->cap[a] == 1 && !c->charge_fixed[a] && !c->charge[a] && sigma(c, a) == 2 && (c->hyb[a] == 1 || !c->xyz);
}

// What a component's leftovers can be resolved with by changing one atom: an N (or ring O, S) with three sigma bonds
// that becomes cationic and takes a pi bond (nitro, N-oxide, guanidinium, pyridinium, protonated imidazole), a
// linear N that takes a second one (azide, diazo), or,
// without hydrogens, an N that is protonated and gives one up (the N-H of pyrrole, imidazole, indole, amide). Every
// candidate is tried; the one leaving the least wanting is kept, an amide N last, then the nearest. A terminal O or S
// only takes the atom it is bonded to (nitro, N-oxide): a cation further away would turn a phenolate into a
// quinoid with two more charges.
static bool try_fix(chem_t* c, uint32_t leftover, uint32_t* buf, uint32_t* buf2, int32_t* vbase) {
    mark_begin(c);
    const size_t n = collect_component(c, leftover, buf);
    if (!n) return false;

    uint32_t cand[MAX_FIX_CANDIDATES];
    size_t nc = 0;
    const bool none = c->mode[leftover] == HMODE_NONE;
    const bool chalcogen = is_chalcogen(c->z[leftover]) && c->heavy[leftover] == 1;
    for (size_t h = 0; h < n && nc < MAX_FIX_CANDIDATES; ++h) {
        const uint32_t a = buf[h];
        if (chalcogen && a != leftover) break;
        if (none && c->z[a] == MD_Z_N && c->cap[a] == 1 && !c->charge_fixed[a] && !c->charge[a] && !c->metal_bound[a] && sigma(c, a) == 2 && c->hyb[a] != 1) {
            cand[nc++] = a | PROTONATE_BIT;
        }
        if (linear_cation(c, a) && nc < MAX_FIX_CANDIDATES) cand[nc++] = a;
        // Expandable neighbours of the atom, the shortest bond first
        uint32_t loc[8];
        float key[8];
        int nl = 0;
        for (uint32_t k = c->off[a]; k < c->off[a + 1] && nl < 8; ++k) {
            const uint32_t x = c->nbr[k];
            if (c->cap[x] || c->charge_fixed[x] || c->charge[x] || !expandable(c, x)) continue;
            bool dup = false;
            for (size_t j = 0; j < nc; ++j) dup |= (cand[j] & ~PROTONATE_BIT) == x;
            if (dup) continue;
            const float r = c->ratio[c->nbr_bond[k]];
            int j = nl++;
            while (j > 0 && key[j - 1] > r) { loc[j] = loc[j - 1]; key[j] = key[j - 1]; --j; }
            loc[j] = x;
            key[j] = r;
        }
        for (int j = 0; j < nl && nc < MAX_FIX_CANDIDATES; ++j) cand[nc++] = loc[j];
    }

    int best = -1, best_after = INT32_MAX;
    float best_penalty = FLT_MAX;
    for (size_t k = 0; k < nc; ++k) {
        const uint32_t x = cand[k] & ~PROTONATE_BIT;
        const int before = region(c, leftover, x, false, buf2, vbase);
        float penalty = (float)fix_penalty(c, cand[k], leftover);
        if (cand[k] & PROTONATE_BIT) {
            // The N whose pi bond is the longest (an amide or anilide N rather than that of a pyridine)
            for (uint32_t j = c->off[x]; j < c->off[x + 1]; ++j) if (c->pi_match[c->nbr_bond[j]]) penalty += 1.0f - c->ratio[c->nbr_bond[j]];
        }
        apply_fix(c, cand[k], false);
        const int after = region(c, leftover, x, true, buf2, vbase);
        if (after < before && c->left[x] == 0 && (after < best_after || (after == best_after && penalty < best_penalty))) {
            best = (int)k;
            best_after = after;
            best_penalty = penalty;
        }
        apply_fix(c, cand[k], true);
        region(c, leftover, x, true, buf2, vbase);
    }
    if (best < 0) return false;
    apply_fix(c, cand[best], false);
    region(c, leftover, cand[best] & ~PROTONATE_BIT, true, buf2, vbase);
    return true;
}

// ### RINGS ###

// An exocyclic double bond from a ring atom to x takes the ring atom's electron away when x is the more
// electronegative: O, N, S or Se from a carbon (C=O, C=N, C=S, C=Se), O from a nitrogen
static bool exocyclic_acceptor(md_atomic_number_t ring, md_atomic_number_t x) {
    if (ring == MD_Z_C || ring == MD_Z_B) return x == MD_Z_O || x == MD_Z_N || is_heavy_chalcogen(x);
    if (ring == MD_Z_N || ring == MD_Z_P || ring == 33) return x == MD_Z_O;
    return false;
}

// Pi electrons an atom gives to a ring (or fused ring system) whose atoms are marked in 'member', -1 if it cannot
// take part. This is the aromaticity model of RDKit, which most SMILES are written in, so that a lowercase atom of a
// SMILES pattern finds the atom here:
//   - a double bond within the ring gives one electron, as does one leaving the ring into a fused ring (the Kekule
//     structure of naphthalene or anthracene)
//   - a double bond leaving the ring to a more electronegative atom takes the electron away (the C=O of a pyridone,
//     uracil, guanine or flavin): none, but the ring may still be aromatic. Any other exocyclic double bond (C=C)
//     breaks the ring.
//   - a lone pair gives two (pyrrole N-H, furan O, thiophene S), an empty orbital none (B, a carbocation)
// The PDB chemical component dictionary differs in the second point: there a ring with an exocyclic C=O is not aromatic.
static int ring_electrons(const chem_t* c, uint32_t a, const uint8_t* member) {
    const md_atomic_number_t z = c->z[a];
    if (!(z == MD_Z_C || z == MD_Z_N || z == MD_Z_O || is_heavy_chalcogen(z) || z == MD_Z_B || z == MD_Z_P || z == 33)) return -1;
    int in_double = 0, exo_ring_double = 0, exo_electronegative = 0;
    for (uint32_t k = c->off[a]; k < c->off[a + 1]; ++k) {
        const uint32_t b = c->nbr_bond[k];
        const int pi = pi_of_bond(c, b);
        if (pi <= 0) continue;
        if (pi >= 2) return -1;
        const uint32_t x = c->nbr[k];
        if (member[x]) in_double += 1;
        else if (c->ring_bond[b]) exo_ring_double += 1;
        else if (exocyclic_acceptor(z, c->z[x])) exo_electronegative += 1;
        else return -1;
    }
    if (in_double) return (in_double == 1 && !exo_electronegative) ? 1 : -1;
    if (exo_ring_double) return exo_electronegative ? -1 : 1;
    if (exo_electronegative) return exo_electronegative == 1 ? 0 : -1;
    const int q = c->charge[a];
    const int s = sigma(c, a);
    if (z == MD_Z_C) return q < 0 ? 2 : (q > 0 ? 0 : -1);
    if (z == MD_Z_B) return s == 3 ? 0 : -1;
    if (z == MD_Z_N || z == MD_Z_P || z == 33) return ((s == 3 && q == 0) || (s == 2 && q == -1)) ? 2 : -1;
    if (z == MD_Z_O || is_heavy_chalcogen(z)) return (s == 2 && q == 0) ? 2 : -1;
    return -1;
}

static bool huckel(const chem_t* c, const int32_t* atoms, size_t n, uint8_t* member) {
    for (size_t k = 0; k < n; ++k) member[atoms[k]] = 1;
    int e = 0;
    bool ok = true;
    for (size_t k = 0; k < n && ok; ++k) {
        const int x = ring_electrons(c, (uint32_t)atoms[k], member);
        if (x < 0) ok = false; else e += x;
    }
    for (size_t k = 0; k < n; ++k) member[atoms[k]] = 0;
    return ok && e % 4 == 2;
}

static void mark_aromatic_ring(chem_t* c, const int32_t* ring, size_t n, uint8_t* atom_aromatic) {
    for (size_t k = 0; k < n; ++k) {
        const uint32_t a = (uint32_t)ring[k];
        const uint32_t b = (uint32_t)ring[(k + 1) % n];
        atom_aromatic[a] = 1;
        const int bi = find_bond(c, a, b);
        if (bi >= 0) c->aromatic[bi] = 1;
    }
}

static void perceive_aromaticity(chem_t* c, uint8_t* atom_aromatic) {
    const md_index_data_t* rings = &c->sys->ring;
    const size_t num = md_index_data_num_ranges(rings);
    if (!num) return;
    uint8_t* member = md_alloc(c->alloc, c->N);
    MEMSET(member, 0, c->N);
    uint8_t* aromatic_ring = md_alloc(c->alloc, num);
    MEMSET(aromatic_ring, 0, num);

    for (size_t r = 0; r < num; ++r) {
        const int32_t* ring = md_index_range_beg(rings, r);
        const size_t n = md_index_range_size(rings, r);
        if (n < 4) continue;
        if (huckel(c, ring, n, member)) aromatic_ring[r] = 1;
    }
    // Fused pairs: rings that are not aromatic alone but are as a system with a neighbour sharing a bond. The
    // neighbours are found through the rings of each atom.
    uint32_t* ring_off = md_alloc(c->alloc, sizeof(uint32_t) * (c->N + 1));
    MEMSET(ring_off, 0, sizeof(uint32_t) * (c->N + 1));
    for (size_t r = 0; r < num; ++r) {
        const int32_t* ring = md_index_range_beg(rings, r);
        const size_t n = md_index_range_size(rings, r);
        for (size_t k = 0; k < n; ++k) if (ring[k] >= 0 && (size_t)ring[k] < c->N) ring_off[ring[k] + 1] += 1;
    }
    for (size_t i = 0; i < c->N; ++i) ring_off[i + 1] += ring_off[i];
    uint32_t* ring_idx = md_alloc(c->alloc, sizeof(uint32_t) * MAX(ring_off[c->N], 1));
    uint32_t* fill = md_alloc(c->alloc, sizeof(uint32_t) * MAX(c->N, 1));
    MEMCPY(fill, ring_off, sizeof(uint32_t) * c->N);
    for (size_t r = 0; r < num; ++r) {
        const int32_t* ring = md_index_range_beg(rings, r);
        const size_t n = md_index_range_size(rings, r);
        for (size_t k = 0; k < n; ++k) if (ring[k] >= 0 && (size_t)ring[k] < c->N) ring_idx[fill[ring[k]]++] = (uint32_t)r;
    }

    int32_t buf[64];
    for (size_t r = 0; r < num; ++r) {
        if (aromatic_ring[r]) continue;
        const int32_t* ra = md_index_range_beg(rings, r);
        const size_t na = md_index_range_size(rings, r);
        if (na < 4) continue;
        for (size_t i = 0; i < na && !aromatic_ring[r]; ++i) {
            if (ra[i] < 0 || (size_t)ra[i] >= c->N) continue;
            for (uint32_t t = ring_off[ra[i]]; t < ring_off[ra[i] + 1]; ++t) {
                const size_t s = ring_idx[t];
                // Each pair once from the ring of lower index; skip pairs met through an earlier shared atom
                if (s == r) continue;
                const int32_t* rb = md_index_range_beg(rings, s);
                const size_t nb = md_index_range_size(rings, s);
                if (nb < 4 || na + nb > ARRAY_SIZE(buf)) continue;
                bool earlier = false;
                for (size_t j = 0; j < i && !earlier; ++j) for (size_t l = 0; l < nb; ++l) earlier |= ra[j] == rb[l];
                if (earlier) continue;
                size_t shared = 0;
                for (size_t j = 0; j < na; ++j) for (size_t l = 0; l < nb; ++l) shared += ra[j] == rb[l];
                if (shared < 2) continue;
                size_t n = 0;
                for (size_t j = 0; j < na; ++j) buf[n++] = ra[j];
                for (size_t l = 0; l < nb; ++l) {
                    bool dup = false;
                    for (size_t j = 0; j < na; ++j) dup |= rb[l] == ra[j];
                    if (!dup) buf[n++] = rb[l];
                }
                if (huckel(c, buf, n, member)) {
                    aromatic_ring[r] = 1;
                    aromatic_ring[s] = 1;
                    break;
                }
            }
        }
    }
    for (size_t r = 0; r < num; ++r) {
        if (aromatic_ring[r]) mark_aromatic_ring(c, md_index_range_beg(rings, r), md_index_range_size(rings, r), atom_aromatic);
    }
}

// ### DELOCALIZED GROUPS ###

static void perceive_delocalized(chem_t* c, const uint8_t* atom_aromatic) {
    for (size_t i = 0; i < c->N; ++i) {
        if (!c->part[i] || c->z[i] == MD_Z_H) continue;
        const md_atomic_number_t z = c->z[i];
        // Anionic: two or more terminal O/S on one centre, one of them double bonded, one charged
        int terminal = 0, dbl = 0, neg = 0;
        for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
            const uint32_t x = c->nbr[k];
            if (!is_chalcogen(c->z[x]) || c->heavy[x] != 1 || c->h_exp[x] + c->h_imp[x] != 0) continue;
            terminal += 1;
            if (pi_of_bond(c, c->nbr_bond[k]) == 1) dbl += 1;
            if (c->charge[x] < 0) neg += 1;
        }
        if (terminal >= 2 && dbl >= 1 && neg >= 1) {
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                const uint32_t x = c->nbr[k];
                if (is_chalcogen(c->z[x]) && c->heavy[x] == 1 && c->h_exp[x] + c->h_imp[x] == 0) c->delocalized[c->nbr_bond[k]] = 1;
            }
            continue;
        }
        // 1,3-dipoles: an sp N+ between two atoms, one of them negative (azide, diazo)
        if (z == MD_Z_N && c->charge[i] > 0 && c->heavy[i] == 2 && sigma(c, (uint32_t)i) == 2 && pi_sum(c, (uint32_t)i) == 2) {
            bool neg = false;
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) neg |= c->charge[c->nbr[k]] < 0;
            if (neg) {
                for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) c->delocalized[c->nbr_bond[k]] = 1;
                continue;
            }
        }
        // Cationic: a carbon with two or three flat N, one of them N+ double bonded (amidinium, guanidinium)
        if (z == MD_Z_C && !atom_aromatic[i]) {
            int n_flat = 0, n_plus = 0;
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                const uint32_t x = c->nbr[k];
                if (c->z[x] != MD_Z_N || atom_aromatic[x]) continue;
                const int pi = pi_of_bond(c, c->nbr_bond[k]);
                if (pi == 1 && c->charge[x] > 0) n_plus += 1;
                else if (pi == 0 && sigma(c, x) == 3 && pi_sum(c, x) == 0) n_flat += 1;
            }
            if (n_plus == 1 && n_flat >= 1) {
                for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                    const uint32_t x = c->nbr[k];
                    if (c->z[x] != MD_Z_N || atom_aromatic[x]) continue;
                    const int pi = pi_of_bond(c, c->nbr_bond[k]);
                    if ((pi == 1 && c->charge[x] > 0) || (pi == 0 && sigma(c, x) == 3 && pi_sum(c, x) == 0)) {
                        c->delocalized[c->nbr_bond[k]] = 1;
                    }
                }
            }
        }
    }
}

// ### PROTONATION AT PH 7 (no hydrogens) ###

static bool adjacent_to_pi(const chem_t* c, uint32_t i, const uint8_t* atom_aromatic) {
    for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
        const uint32_t x = c->nbr[k];
        if (atom_aromatic[x] || pi_sum(c, x) > 0) return true;
    }
    return false;
}

static void protonate_ph7(chem_t* c, const uint8_t* atom_aromatic) {
    for (size_t i = 0; i < c->N; ++i) {
        if (c->mode[i] != HMODE_NONE || c->z[i] != MD_Z_N || c->charge_fixed[i] || c->charge[i] || c->metal_bound[i]) continue;
        const int s = sigma(c, (uint32_t)i);
        const int pi = pi_sum(c, (uint32_t)i);
        if (atom_aromatic[i]) continue;
        // Aliphatic amine
        if (pi == 0 && s == 3 && !adjacent_to_pi(c, (uint32_t)i, atom_aromatic)) {
            bool all_carbon = true;
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) all_carbon &= c->z[c->nbr[k]] == MD_Z_C || c->z[c->nbr[k]] == MD_Z_H;
            if (all_carbon) {
                c->h_imp[i] += 1;
                c->charge[i] = 1;
            }
            continue;
        }
        // Imine N of an acyclic amidine or guanidine
        if (pi == 1 && s == 2 && !c->in_ring[i]) {
            for (uint32_t k = c->off[i]; k < c->off[i + 1]; ++k) {
                const uint32_t x = c->nbr[k];
                if (c->z[x] != MD_Z_C || c->in_ring[x] || pi_of_bond(c, c->nbr_bond[k]) != 1) continue;
                int amino = 0;
                for (uint32_t j = c->off[x]; j < c->off[x + 1]; ++j) {
                    const uint32_t y = c->nbr[j];
                    if (y != i && c->z[y] == MD_Z_N && pi_sum(c, y) == 0 && !atom_aromatic[y]) amino += 1;
                }
                if (amino >= 1) {
                    c->h_imp[i] += 1;
                    c->charge[i] = 1;
                }
                break;
            }
        }
    }
}

// ### PERCEPTION ###

bool md_chem_perceive(md_system_t* sys, const md_system_state_t* reference, uint32_t flags) {
    if (!sys || !sys->alloc) {
        MD_LOG_ERROR("Chemistry perception: missing system");
        return false;
    }
    const size_t N = sys->atom.count;
    if (N == 0) return true;
    if (N > INT32_MAX / 4) {
        MD_LOG_ERROR("Chemistry perception: too many atoms");
        return false;
    }

    if (sys->bond.count && md_index_data_num_ranges(&sys->ring) == 0) {
        if (!sys->bond.conn.offset) md_bond_build_connectivity(&sys->bond, N, sys->alloc);
        md_util_system_infer_rings(sys);
    }

    if (!sys->atom.flags) {
        md_array_resize(sys->atom.flags, N, sys->alloc);
        MEMSET(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* arena = md_temp_allocator(temp);
    chem_t c = {
        .sys   = sys,
        .xyz   = (md_system_state_has_coords(reference) && reference->num_atoms == N) ? reference->xyz : NULL,
        .cell  = reference ? &reference->unitcell : NULL,
        .N     = N,
        .flags = flags,
        .alloc = arena,
    };
    build_graph(&c);
    assign_modes(&c);
    mark_rings(&c);

    // The formal charges a file gives: the atom/formal_charge column its loader publishes, NAN where it gives none.
    // The charges in md_atom_data_t are this function's output and are not read back, so running it again starts
    // from the same place.
    const md_attribute_t* q_attr = md_attributes_find(&sys->attributes, STR_LIT("atom/formal_charge"));
    const float* q_given = (const float*)md_attribute_view(q_attr, MD_ATTRIBUTE_TYPE_F32, 1, 1);
    if (q_given && md_attribute_slice_count(q_attr, md_attribute_slice_all()) != N) q_given = NULL;

    c.heavy        = md_alloc(arena, N);
    c.h_exp        = md_alloc(arena, N);
    c.h_imp        = md_alloc(arena, N);
    c.charge       = md_alloc(arena, N);
    c.charge_fixed = md_alloc(arena, N);
    c.hyb          = md_alloc(arena, N);
    c.cap          = md_alloc(arena, N);
    c.left         = md_alloc(arena, N);
    MEMSET(c.h_imp, 0, N);
    MEMSET(c.cap, 0, N);
    MEMSET(c.left, 0, N);
    for (size_t i = 0; i < N; ++i) {
        uint8_t nh = 0, nx = 0;
        for (uint32_t k = c.off[i]; k < c.off[i + 1]; ++k) {
            if (c.z[c.nbr[k]] == MD_Z_H) nh += 1; else nx += 1;
        }
        c.heavy[i] = nx;
        c.h_exp[i] = nh;
        c.charge[i] = 0;
        c.charge_fixed[i] = 0;
        if (q_given && !float_bits_nan(q_given[i])) {
            const long q = lroundf(q_given[i]);
            c.charge[i] = (int8_t)CLAMP(q, -8, 8);
            c.charge_fixed[i] = q != 0;
        }
    }

    const size_t nb = c.num_bonds;
    c.pi_fixed    = md_alloc(arena, MAX(nb, 1));
    c.pi_match    = md_alloc(arena, MAX(nb, 1));
    c.given       = md_alloc(arena, MAX(nb, 1));
    MEMSET(c.given, 0, MAX(nb, 1));
    c.aromatic    = md_alloc(arena, MAX(nb, 1));
    c.delocalized = md_alloc(arena, MAX(nb, 1));
    MEMSET(c.pi_fixed, 0, MAX(nb, 1));
    MEMSET(c.pi_match, 0, MAX(nb, 1));
    MEMSET(c.aromatic, 0, MAX(nb, 1));
    MEMSET(c.delocalized, 0, MAX(nb, 1));

    // Orders the file gives are kept
    for (size_t b = 0; b < nb; ++b) {
        if (!c.bond_part[b]) continue;
        const md_bond_flags_t f = sys->bond.flags[b];
        const int order = md_bond_order(f);
        if (order > 0 && !(f & MD_BOND_FLAG_ORDER_PERCEIVED)) {
            c.pi_fixed[b] = (uint8_t)(order - 1);
            c.given[b] = 1;
            if (f & MD_BOND_FLAG_AROMATIC)    c.aromatic[b] = 1;
            if (f & MD_BOND_FLAG_DELOCALIZED) c.delocalized[b] = 1;
        }
    }

    // Geometry, then implicit hydrogens (the geometry of an atom counts its implicit hydrogens as neighbours, so twice)
    for (size_t i = 0; i < N; ++i) {
        c.hyb[i] = 0;
        if (!c.part[i] || c.z[i] == MD_Z_H) continue;
        // With all hydrogens explicit an atom saturated by its sigma bonds has no use for its geometry (water, CH3)
        if (c.mode[i] == HMODE_ALL && sigma(&c, (uint32_t)i) >= base_valence(c.z[i]) && !c.charge[i]) { c.hyb[i] = 3; continue; }
        c.hyb[i] = geometric_hybridization(&c, (uint32_t)i);
    }
    implicit_hydrogens(&c);
    charge_n_oxides(&c);

    // Pi bonds wanted
    for (size_t i = 0; i < N; ++i) {
        if (!c.part[i] || c.z[i] == MD_Z_H) continue;
        const int v = target_valence(&c, (uint32_t)i);
        int fixed = 0;
        for (uint32_t k = c.off[i]; k < c.off[i + 1]; ++k) fixed += c.pi_fixed[c.nbr_bond[k]];
        c.cap[i] = (uint8_t)CLAMP(v - sigma(&c, (uint32_t)i) - fixed, 0, 3);
    }

    // Triple bonds first: two atoms that each want two, joined by a short bond
    for (size_t b = 0; b < nb; ++b) {
        if (!c.bond_part[b] || c.given[b]) continue;
        const uint32_t x = (uint32_t)sys->bond.pairs[b].idx[0];
        const uint32_t y = (uint32_t)sys->bond.pairs[b].idx[1];
        if (c.cap[x] >= 2 && c.cap[y] >= 2 && (!c.xyz || c.ratio[b] < 0.86f || (c.hyb[x] == 1 && c.hyb[y] == 1))) {
            c.pi_fixed[b] += 2;
            c.cap[x] -= 2;
            c.cap[y] -= 2;
        }
    }

    uint32_t* buf  = md_alloc(arena, sizeof(uint32_t) * N);
    uint32_t* buf2 = md_alloc(arena, sizeof(uint32_t) * N);
    c.mark = md_alloc(arena, sizeof(uint32_t) * N);
    MEMSET(c.mark, 0, sizeof(uint32_t) * N);
    c.mark_gen = 0;
    int32_t*  vbase = md_alloc(arena, sizeof(int32_t) * N);
    for (size_t i = 0; i < N; ++i) vbase[i] = -1;

    match_all(&c, buf, vbase, RELOCATE_GEOMETRY);

    // Leftovers: fix by a cation or a protonation where that completes the matching, else by the atom itself
    for (int pass = 0; pass < 2; ++pass) {
        // What pass 0 leaves is moved off carbons, to the atoms that resolve it themselves
        if (pass == 1) match_all(&c, buf, vbase, RELOCATE_SCORE);
        for (size_t i = 0; i < N; ++i) {
            if (!c.left[i]) continue;
            const md_atomic_number_t z = c.z[i];
            const bool none = c.mode[i] == HMODE_NONE;
            if (pass == 0) {
                // Charge separation or protonation elsewhere is only tried for what the atom cannot resolve
                // itself: a carbon, or (with hydrogens) an N; a terminal O or S is an anion or, without
                // hydrogens, protonated, unless an N+ (nitro) completes it.
                const bool terminal_chalcogen = is_chalcogen(z) && c.heavy[i] == 1;
                if (terminal_chalcogen && none && !in_anionic_group(&c, (uint32_t)i)) continue;
                try_fix(&c, (uint32_t)i, buf, buf2, vbase);
                continue;
            }
            // Pass 1: the atom resolves what is left itself
            while (c.left[i]) {
                // Without hydrogens a metal bound atom is the anion (the thiolate of a zinc finger cysteine, the
                // N of a porphyrin, an imidazolate bridging two metals)
                const bool protonate = none && !c.metal_bound[i];
                if (is_chalcogen(z) && c.heavy[i] <= 1) {
                    if (protonate && !((c.flags & MD_CHEM_FLAG_PROTONATE_PH7) && in_anionic_group(&c, (uint32_t)i))) {
                        c.h_imp[i] += 1;
                    } else {
                        c.charge[i] -= 1;
                    }
                } else if (z == MD_Z_N) {
                    if (protonate) c.h_imp[i] += 1;
                    else c.charge[i] -= 1;
                } else if ((z == MD_Z_P || is_heavy_chalcogen(z) || z == 33) && c.left[i] >= 2) {
                    // An expanded valence with no partners for its pi bonds: the lower one (S(IV), not S(VI)2+)
                    c.left[i] -= 1;
                    c.cap[i] -= 1;
                } else if (z == MD_Z_P || is_heavy_chalcogen(z) || z == 33) {
                    c.charge[i] += 1;       // Phosphonium, sulfonium
                } else if (z == MD_Z_C || z == MD_Z_Si) {
                    if (none) c.h_imp[i] += 1;
                    else break;             // A carbocation or radical: left as it is
                } else {
                    break;
                }
                c.left[i] -= 1;
                c.cap[i] -= 1;
            }
        }
    }

    uint8_t* atom_aromatic = md_alloc(arena, N);
    MEMSET(atom_aromatic, 0, N);
    perceive_aromaticity(&c, atom_aromatic);
    if (flags & MD_CHEM_FLAG_PROTONATE_PH7) protonate_ph7(&c, atom_aromatic);
    perceive_delocalized(&c, atom_aromatic);

    // Monatomic ions
    for (size_t ci = 0; ci < sys->component.count && sys->component.atom_offset; ++ci) {
        const md_urange_t r = md_component_atom_range(&sys->component, ci);
        if (r.end - r.beg != 1 || r.beg >= N) continue;
        const uint32_t a = r.beg;
        if (c.charge_fixed[a] || (sys->bond.conn.offset && md_bond_conn_count(&sys->bond, a))) continue;
        const int q = monatomic_ion_charge(c.z[a] ? c.z[a] : md_atom_atomic_number(&sys->atom, a));
        if (q) c.charge[a] = (int8_t)q;
    }

    // ### WRITE ###
    for (size_t b = 0; b < nb; ++b) {
        if (!c.bond_part[b]) continue;
        md_bond_flags_t f = sys->bond.flags[b];
        const bool given = md_bond_order(f) > 0 && !(f & MD_BOND_FLAG_ORDER_PERCEIVED);
        if (given) continue;
        f = md_bond_flags_set_order(f, 1 + pi_of_bond(&c, (uint32_t)b));
        f &= ~(MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED);
        if (c.aromatic[b])    f |= MD_BOND_FLAG_AROMATIC;
        if (c.delocalized[b]) f |= MD_BOND_FLAG_DELOCALIZED;
        f |= MD_BOND_FLAG_ORDER_PERCEIVED;
        sys->bond.flags[b] = f;
    }

    md_array_resize(sys->atom.formal_charge, N, sys->alloc);
    md_array_resize(sys->atom.hydrogen_count, N, sys->alloc);
    for (size_t i = 0; i < N; ++i) {
        sys->atom.formal_charge[i] = c.charge[i];
        sys->atom.hydrogen_count[i] = (c.part[i] && c.z[i] != MD_Z_H) ? (uint8_t)(c.h_exp[i] + c.h_imp[i]) : 0;
    }

    // Hybridization
    for (size_t i = 0; i < N; ++i) {
        md_flags_t f = sys->atom.flags[i] & ~(MD_FLAG_SP | MD_FLAG_SP2 | MD_FLAG_SP3 | MD_FLAG_AROMATIC);
        const md_atomic_number_t z = c.z[i];
        if (c.part[i] && z != MD_Z_H && !is_halogen_z(z)) {
            const int pi = pi_sum(&c, (uint32_t)i);
            const int s = sigma(&c, (uint32_t)i);
            if (atom_aromatic[i]) {
                f |= MD_FLAG_SP2 | MD_FLAG_AROMATIC;
            } else if ((z == MD_Z_P || is_heavy_chalcogen(z) || z == 33) && s >= 3) {
                f |= MD_FLAG_SP3;                       // Tetrahedral, whatever the formal double bonds
            } else if (pi >= 2) {
                f |= MD_FLAG_SP;
            } else if (pi == 1 || z == MD_Z_B) {
                f |= MD_FLAG_SP2;
            } else if (z == MD_Z_N && adjacent_to_pi(&c, (uint32_t)i, atom_aromatic) && (s < 3 || !c.xyz || planar_angle_sum(&c, (uint32_t)i) >= 345.0f || c.off[i + 1] - c.off[i] < 3)) {
                f |= MD_FLAG_SP2;                       // Amide, aniline, sulfonamide N: conjugated and flat
            } else {
                f |= MD_FLAG_SP3;
            }
        }
        sys->atom.flags[i] = f;
    }

    md_temp_end(temp);
    return true;
}

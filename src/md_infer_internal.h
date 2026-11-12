#pragma once

// What the inference of a system (md_system_infer and its steps, in md_infer.c, md_bond.c, md_ring.c, md_atom.c,
// md_hierarchy.c, md_backbone.c, md_lipid.c and md_martini.c) shares between its files. Not a public header.

#include <md_infer.h>
#include <md_util_geom.h>
#include <md_system.h>
#include "md_system_fields.h"
#include <md_chem.h>
#include <md_dssp.h>
#include <md_element.h>
#include <core/md_coord_stream.h>
#include <core/md_format.h>
#include <core/md_compiler.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_simd.h>
#include <core/md_spatial_acc.h>
#include <core/md_executor.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>
#include <core/md_fifo.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

// Tables
// The zero slot is 'Xx', not 'Unknown': winioctl.h declares an enumerator by the latter name at
// file scope, so on any build whose include chain reaches windows.h the two collide. 'Xx' is what
// element_symbols[0] holds anyway.
enum {
    Xx = 0,
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V,
    Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru,
    Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb,
    Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh,
    Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og, Num_Elements
};


static inline md_array(uint64_t) make_bitfield(size_t num_bits, md_allocator_i* alloc) {
    size_t num_elem = MD_DIV_UP(num_bits, 64);
    if (num_elem == 0) {
        return NULL;
    }
    md_array(uint64_t) bits = md_array_create(uint64_t, num_elem, alloc);
    md_memset(bits, 0, md_array_bytes(bits));
    return bits;
}

static inline void bitfield_set_all(uint64_t* bits, size_t num_bits) {
    // No bits is a NULL bitfield (make_bitfield), which memset must not be handed even for zero bytes
    if (num_bits == 0) return;
    size_t num_bytes = MD_DIV_UP(num_bits, 8);
    md_memset(bits, 0xFF, num_bytes);
}

static inline bool bitfield_test_bit(const uint64_t* bits, int64_t idx) {
    return (bits[idx >> 6] & (1ULL << (idx & 63)));
}

static inline void bitfield_set_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] |= (1ULL << (idx & 63));
}

static inline void bitfield_clear_bit(uint64_t* bits, int64_t idx) {
    bits[idx >> 6] &= ~(1ULL << (idx & 63));
}

static inline size_t bitfield_popcount(const uint64_t* bits, size_t num_bits) {
    if (bits == NULL || num_bits == 0) {
        return 0;
    }
    const uint64_t* it = bits;
    const uint64_t* end = bits + MD_DIV_UP(num_bits, 64) - 1;
    size_t count = 0;
    while (it != end) {
        count += (size_t)md_popcnt64(*it);
        it++;
    }
    const uint64_t mask = (1ULL << (num_bits & 63)) - 1;
    count += (size_t)md_popcnt64(*end & mask);
    return count;
}

static inline float angle(md_f32x3_t a, md_f32x3_t b, md_f32x3_t c) {
    md_f32x3_t v0 = md_f32x3_normalize(md_f32x3_sub(a, b));
    md_f32x3_t v1 = md_f32x3_normalize(md_f32x3_sub(c, b));
    return (float)MD_RAD_TO_DEG(acosf(md_f32x3_dot(v0, v1)));
}


typedef struct {
    int atom_i;
    int atom_j;
    float dist;
} bond_pair_t;

// A water molecule whose atoms beyond its O and two H are virtual sites (4 and 5 site models)
// One oxygen, two hydrogens and possibly virtual sites (4 and 5 site models): a water molecule, whatever it is called
static inline bool atoms_are_water(const md_system_t* sys, md_urange_t range) {
    uint32_t num_o = 0, num_h = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
        if (z == MD_Z_O) num_o += 1;
        else if (z == MD_Z_H) num_h += 1;
        else if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_VIRTUAL_SITE) return false;
    }
    return num_o == 1 && num_h == 2;
}

// ### MARTINI ###

// The versions of Martini, which give beads of the same names different masses and sizes, and some residues different
// beads. Which one a system is from its residue and bead names tell (mdi_martini_system).
enum {
    MARTINI_2   = 1,    // martini22, elnedyn21, elnedyn22
    MARTINI_2P  = 2,    // Polarizable Martini 2: martini22p, elnedyn22p
    MARTINI_3   = 4,    // martini3001
    MARTINI_ANY = MARTINI_2 | MARTINI_2P | MARTINI_3,
};

// The radius of a regular Martini bead (sigma 4.7 A), as md_tpr gives it: half the distance of the minimum of the
// Lennard-Jones potential, 2^(1/6) sigma / 2
#define MARTINI_R   2.638f

typedef struct {
    const char*         name;
    float               mass;
    float               radius;     // From the sigma of the bead type with itself, as md_tpr gives it
    md_particle_kind_t  particle;   // A massless bead is a virtual site; one without Lennard-Jones a charge site
} martini_db_bead_t;

typedef struct {
    const char*         name;
    uint32_t            versions;
    md_residue_kind_t kind;
    uint16_t            bead_offset;    // Into martini_db_residue_beads
    uint8_t             num_beads;
    uint16_t            bond_offset;    // Into martini_db_residue_bonds, pairs of beads of the residue
    uint8_t             num_bonds;
} martini_db_residue_t;

typedef struct {
    const char* name;       // A residue of a polymer (polyply), and its bead bonded to
    const char* bead;
    const char* next_name;  // the bead of the residue after it
    const char* next_bead;
    uint32_t    versions;
} martini_db_link_t;

typedef struct {
    uint32_t version;
    bool     martini;   // A Martini system: residues of several beads which are Martini's, Martini water or polymers
} martini_system_t;

static inline md_str_t label_of_array(const void* data, size_t i) {
    return ((const md_str_t*)data)[i];
}

static inline md_str_t label_of_atom(const void* data, size_t i) {
    return md_atom_name((const md_atom_data_t*)data, i);
}

static inline void set_residue_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind) {
    md_system_set_residue_kind(sys, res_idx, kind);
}

static inline int32_t uf_find(int32_t* parent, int32_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

// Returns true if a and b were in different sets
static inline bool uf_union(int32_t* parent, int32_t a, int32_t b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a == b) return false;
    if (a < b) parent[b] = a; else parent[a] = b;
    return true;
}

static inline md_f32x3_t hcl_to_rgb(float h, float c, float l) {
    float r = 0, g = 0, b = 0;

    if (c == 0) {
        r = g = b = l;
    } else {
        float h_prime = h * 6.0f;
        float x = c * (1.0f - fabsf(fmodf(h_prime, 2.0f) - 1.0f));
        float m = l - c * 0.5f;

        if (0.0f <= h_prime && h_prime < 1.0f) {
            r = c; g = x; b = 0;
        } else if (1.0f <= h_prime && h_prime < 2.0f) {
            r = x; g = c; b = 0;
        } else if (2.0f <= h_prime && h_prime < 3.0f) {
            r = 0; g = c; b = x;
        } else if (3.0f <= h_prime && h_prime < 4.0f) {
            r = 0; g = x; b = c;
        } else if (4.0f <= h_prime && h_prime < 5.0f) {
            r = x; g = 0; b = c;
        } else if (5.0f <= h_prime && h_prime < 6.0f) {
            r = c; g = 0; b = x;
        }

        r += m;
        g += m;
        b += m;
    }

    return (md_f32x3_t){r, g, b};
}

// ### SHARED BETWEEN THE FILES OF INFERENCE ###

// md_atom.c
md_particle_type_idx_t mdi_find_or_add_particle_type(md_particle_type_data_t* type, md_str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_particle_type_flags_t flags, md_allocator_i* alloc);

// md_backbone.c
md_atom_idx_t mdi_bonded_atom_with(const md_system_t* sys, md_atom_idx_t a, md_atomic_number_t z, md_atom_flags_t with, md_atom_flags_t without, md_atom_idx_t x, md_atom_idx_t y);
bool mdi_residue_amino_acid_atoms(md_amino_acid_atoms_t* out, const md_system_t* sys, md_urange_t range);
bool mdi_residue_nucleotide_atoms(md_nucleic_acid_atoms_t* out, const md_system_t* sys, md_urange_t range);
bool mdi_polymer_chain_ends_between(const md_system_t* sys, size_t a, size_t b);

// md_martini.c
const martini_db_bead_t* mdi_martini_db_bead(const martini_db_residue_t* res, size_t k);
const martini_db_residue_t* mdi_martini_db_match_molecule(uint32_t version, md_str_t res_name, const md_str_t labels[], size_t num_labels, uint8_t map[]);
const martini_db_residue_t* mdi_martini_db_match_beads(uint32_t version, const md_str_t labels[], const float masses[], size_t num_labels, uint8_t map[]);
const martini_db_residue_t* mdi_martini_db_unknown_variant(md_str_t res_name, const md_str_t labels[], size_t num_labels);
bool mdi_martini_db_linked(md_str_t name);
const martini_db_bead_t* mdi_martini_db_biopolymer_bead(uint32_t version, md_str_t res_name, const md_str_t labels[], size_t num_labels, size_t idx);
martini_system_t mdi_martini_system_ex(const md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* data, const md_executor_t* ex);
martini_system_t mdi_martini_system(const md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* data);
void mdi_infer_coarse_grained_bonds(md_array(bond_pair_t)* candidates, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc);

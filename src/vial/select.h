#pragma once

// The system as the selectors see it, and the static values of sets.
//
// A set is a list of particle indices, ascending. An array of sets is one list per set, back to back, with offsets
// (CSR), and what each set is: a unit of a level of the system (residue 12, chain A), or nothing in particular. Lists
// rather than bitsets, since an array of 880 residues over 7,000 particles would be 880 bitsets of 7,000 bits; the
// lists of an array of units hold each particle once.

#include "vial_internal.h"

struct md_system_t;

typedef enum vial_level_t {
    VIAL_LEVEL_NONE = 0,        // sets that are no unit of the system
    VIAL_LEVEL_PARTICLE,
    VIAL_LEVEL_RESIDUE,
    VIAL_LEVEL_PART,
    VIAL_LEVEL_MOLECULE,
    VIAL_LEVEL_CHAIN,
} vial_level_t;

// The kinds of a part, as the reference names them. mdlib's non-polymer parts are told apart here: a lipid when its
// residues are lipids, an ion when it is one ion, else a ligand.
typedef enum vial_part_kind_t {
    VIAL_PART_LIGAND = 0,
    VIAL_PART_PEPTIDE,
    VIAL_PART_DNA,
    VIAL_PART_RNA,
    VIAL_PART_NUCLEIC,          // a hybrid, or one not told apart
    VIAL_PART_POLYMER,          // another polymer
    VIAL_PART_GLYCAN,
    VIAL_PART_LIPID,
    VIAL_PART_ION,
    VIAL_PART_WATER,
} vial_part_kind_t;

// What is derived from the system is built when first needed: a script that selects no molecules never sorts their
// particles. Built on the compile's arena, so nothing outlives the program.
typedef struct vial_sys_t {
    const struct md_system_t* sys;
    md_allocator_i* alloc;
    uint32_t    num_atoms, num_res, num_parts, num_mols, num_chains, num_types;
    uint8_t*    part_kind;      // [part] vial_part_kind_t
    int32_t*    part_chain;     // [part] the chain of its first residue, -1 for none
    // Molecules, once needed: their particles need not be consecutive
    bool        have_mols;
    int32_t*    mol_chain;      // [molecule] the chain of its first particle, -1 for none
    uint32_t*   mol_offset;     // [molecule + 1] the particles of each molecule, ascending, in mol_atoms
    uint32_t*   mol_atoms;
    // Chains, once needed: the particles filed under each
    bool        have_chains;
    uint32_t*   chain_offset;   // [chain + 1]
    uint32_t*   chain_atoms;
} vial_sys_t;

str_t   vial_sys_chain_id(const vial_sys_t* s, int32_t chain);   // empty for -1
int32_t vial_sys_res_chain(const vial_sys_t* s, uint32_t res);   // -1 for none

// Builds the view of a system; false without one
bool vial_sys_init(vial_sys_t* s, const struct md_system_t* sys, md_allocator_i* alloc);

// One set, or an array of sets
typedef struct vial_sets_t {
    uint32_t    count;          // sets
    bool        single;         // one set, not an array of one
    uint8_t     level;          // vial_level_t: what each set is
    uint32_t*   unit;           // [count] the unit each set is, when level is not NONE
    const struct vial_sets_t* root; // what it was sliced or picked from, transitively, or what common() paired it
                                    // with; NULL for an original. Dimensions over sets of one root pair without W0605.
    uint32_t*   offset;         // [count + 1]
    uint32_t*   idx;            // [offset[count]] particle indices, ascending within a set
} vial_sets_t;

static inline uint32_t vial_sets_size(const vial_sets_t* s, uint32_t i) { return s->offset[i + 1] - s->offset[i]; }
static inline const uint32_t* vial_sets_items(const vial_sets_t* s, uint32_t i) { return s->idx + s->offset[i]; }

vial_sets_t* vial_sets_one(const uint32_t* idx, uint32_t n, md_allocator_i* alloc);    // copies idx
vial_sets_t* vial_sets_all(uint32_t num_atoms, md_allocator_i* alloc);                  // every particle

static inline const struct vial_sets_t* vial_sets_root(const struct vial_sets_t* s) { return s && s->root ? s->root : s; }

// Set algebra on single sets: and, or, and the scope without a
vial_sets_t* vial_sets_and(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc);
vial_sets_t* vial_sets_or(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc);
vial_sets_t* vial_sets_not(const vial_sets_t* scope, const vial_sets_t* a, uint32_t ia, md_allocator_i* alloc);

// Every particle of every set, merged
vial_sets_t* vial_sets_merge(const vial_sets_t* a, md_allocator_i* alloc);

// The sets at positions pos[0..n), in that order; level and units kept
vial_sets_t* vial_sets_pick(const vial_sets_t* a, const uint32_t* pos, uint32_t n, md_allocator_i* alloc);

// The particles of one set at positions pos[0..n): a set, the positions in ascending order
vial_sets_t* vial_sets_pick_particles(const vial_sets_t* a, uint32_t i, const uint32_t* pos, uint32_t n, md_allocator_i* alloc);

// A unit's label as hover and exports write it: MET1/A for a residue, CA/MET1/A for a particle, A for a chain; a part or
// a molecule of one residue as that residue, else its kind or name and chain, peptide/A
size_t vial_unit_label(char* buf, size_t cap, vial_sys_t* s, uint8_t level, uint32_t unit);

// ### SELECTORS ###

typedef enum vial_key_t {
    VIAL_KEY_NONE = 0,
    VIAL_KEY_NAME,
    VIAL_KEY_SERIAL,
    VIAL_KEY_ELEMENT,
    VIAL_KEY_KIND,
    VIAL_KEY_RESID,
    VIAL_KEY_RESNAME,
    VIAL_KEY_CHAIN,
    VIAL_KEY_ID,
    VIAL_KEY_SCOPE,
    VIAL_KEY_COUNT
} vial_key_t;

// One value a keyword matches: a string (a glob where the keyword takes one), an int, or an inclusive range
typedef struct vial_match_t {
    uint8_t     kind;           // 0 string, 1 int, 2 range
    str_t       s;
    int64_t     lo, hi;
} vial_match_t;

typedef struct vial_filter_t {
    uint8_t         key;        // vial_key_t
    vial_match_t*   any;        // md_array: a list means any of its values
} vial_filter_t;

typedef struct vial_select_t {
    uint8_t                 level;          // the selector's: particles is VIAL_LEVEL_PARTICLE
    uint8_t                 particle_kind;  // atoms(), beads(): 1 + md_particle_kind_t; 0 for any
    vial_filter_t*          filters;        // md_array
    const vial_sets_t*      within;         // a positional set or array of sets, or NULL
    const vial_sets_t*      scope;          // NULL for the system
} vial_select_t;

// The keywords of a selector, in the order the reference lists them; false for a keyword it does not take
bool        vial_select_takes(uint8_t level, uint8_t key);
const char* vial_select_keywords(uint8_t level);    // "id" for chains, for E0202's message
uint8_t     vial_select_key(str_t keyword);         // VIAL_KEY_NONE when it is none of them
const char* vial_select_name(uint8_t level);        // "residues"

// Whether a kind= value is known for the level; the known values for E0206's message
bool        vial_select_kind_known(uint8_t level, str_t kind);
const char* vial_select_kinds(uint8_t level);

// Evaluates a selector: one set for particles, one set per unit otherwise, in the system's order
vial_sets_t* vial_select(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc);

// The predefined sets: protein, nucleic, water, ion, lipid, backbone, sidechain, system; NULL for another name
vial_sets_t* vial_select_predefined(vial_sys_t* s, str_t name, md_allocator_i* alloc);

// For W0602: the chains among the residues a resid filter matches within the scope, up to two
uint32_t vial_select_resid_chains(const vial_sys_t* s, const vial_filter_t* f, const vial_sets_t* scope, int32_t out[2]);

// Glob matching with * and ?, as the name keywords take it
bool vial_glob(str_t pattern, str_t text);

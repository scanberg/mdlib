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
    VIAL_LEVEL_TAIL,            // a tail of a lipid (md_lipid_data_t), with its branches and hydrogens
    VIAL_LEVEL_HEAD,            // a head of a lipid: a connected part of what is not a tail, with its hydrogens
    VIAL_LEVEL_RING,            // a ring of the bonds (md_ring_data_t): 3 to 7 atoms, the smallest rings
    VIAL_LEVEL_GROUP,           // a group the data gives (md_group_data_t): any particles, overlapping other groups
    VIAL_LEVEL_FRAGMENT,        // a fragment of a molecule (md_fragment.h): a ring system, a functional group, a piece of
                                // the backbone, a tail; one to an atom, and not within one residue
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
    // Lipid tails, once needed: the particles of each, ascending
    bool        have_tails;
    uint32_t    num_tails;
    uint32_t*   tail_offset;    // [tail + 1]
    uint32_t*   tail_atoms;
    // Lipid heads, once needed: the particles of each, ascending
    bool        have_heads;
    uint32_t    num_heads;
    uint32_t*   head_offset;    // [head + 1]
    uint32_t*   head_atoms;
    // Rings, once needed: the particles of each, ascending (rings that are fused share some)
    bool        have_rings;
    uint32_t    num_rings;
    uint32_t*   ring_offset;    // [ring + 1]
    uint32_t*   ring_atoms;
    // Groups, once needed: the particles of each, ascending
    bool        have_groups;
    uint32_t    num_groups;
    uint32_t*   group_offset;   // [group + 1]
    uint32_t*   group_atoms;
    // Fragments, once needed, in the system's order (by their first particle) as the other units are: the particles of
    // each, ascending, and which of mdlib's each is
    bool        have_frags;
    uint32_t    num_frags;
    uint32_t*   frag_offset;    // [fragment + 1]
    uint32_t*   frag_atoms;
    uint32_t*   frag_idx;       // [fragment] md_fragment_idx_t
    // The atoms of the metallic phases, once needed: 1 for such an atom (md_bond_metallic_atoms)
    bool        have_metallic;
    uint8_t*    metallic;       // [atom]
} vial_sys_t;

md_str_t   vial_sys_chain_id(const vial_sys_t* s, int32_t chain);   // empty for -1
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
// An array of one particle each, in the order given: particles as units, so they pair in that order
vial_sets_t* vial_sets_particles(const int32_t* atoms, uint32_t n, md_allocator_i* alloc);

static inline const struct vial_sets_t* vial_sets_root(const struct vial_sets_t* s) { return s && s->root ? s->root : s; }

// Set algebra on single sets: and, or, and the scope without a
vial_sets_t* vial_sets_and(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc);
vial_sets_t* vial_sets_or(const vial_sets_t* a, uint32_t ia, const vial_sets_t* b, uint32_t ib, md_allocator_i* alloc);
vial_sets_t* vial_sets_not(const vial_sets_t* scope, const vial_sets_t* a, uint32_t ia, md_allocator_i* alloc);

// and, or of each set of an array with one set: an array of as many
vial_sets_t* vial_sets_boolop_each(const vial_sets_t* arr, const vial_sets_t* one, bool is_and, md_allocator_i* alloc);

// Every particle of every set, merged
vial_sets_t* vial_sets_merge(const vial_sets_t* a, md_allocator_i* alloc);

// The sets at positions pos[0..n), in that order; level and units kept
vial_sets_t* vial_sets_pick(const vial_sets_t* a, const uint32_t* pos, uint32_t n, md_allocator_i* alloc);

// The particles of one set at positions pos[0..n): a set, the positions in ascending order
vial_sets_t* vial_sets_pick_particles(const vial_sets_t* a, uint32_t i, const uint32_t* pos, uint32_t n, md_allocator_i* alloc);

// s.residue[i]: within each set of a, the units of u it touches, in u's order, and of those the ones at positions
// pos[0..m), counted from the end when negative, each with only the particles it shares with the set; skip[j] (may be
// NULL) leaves entry j empty. a->count * m sets, set i's entries together. *outside counts the sets that lack a unit at
// some position, whose entries there are empty. The level and units of u where every entry is a whole unit of it.
// The sets of several arrays one after the other, as one array (single sets count as one each)
vial_sets_t* vial_sets_concat(const vial_sets_t* const* parts, uint32_t n, md_allocator_i* alloc);
vial_sets_t* vial_sets_units_at(const vial_sets_t* a, const vial_sets_t* u, const int64_t* pos, const bool* skip, uint32_t m,
                                uint32_t* outside, md_allocator_i* alloc);

// A unit's label as hover and exports write it: MET1/A for a residue, CA/MET1/A for a particle, A for a chain; a part or
// a molecule of one residue as that residue, else its kind or name and chain, peptide/A
size_t vial_unit_label(char* buf, size_t cap, vial_sys_t* s, uint8_t level, uint32_t unit);

// ### THE RUN'S FRAMES ###

// The run's time axis, <run>/time: the time of each frame in ps where it has a unit of time, else its own numbers;
// NULL, and 0 frames, where the run has none
double* vial_run_times(const struct md_system_t* sys, md_str_t run, md_allocator_i* alloc, size_t* num_frames);
// The frame of the run closest to time t (ps, as vial_run_times gives them), and that frame's time: 1 found; 0 where the
// run has no time axis; -1 where t lies outside the run, before its first frame or past its last by more than half
// the spacing of its frames (a run of one frame has only its own time). Equally close frames give the earlier.
int     vial_run_frame_at(const struct md_system_t* sys, md_str_t run, double t, int64_t* frame, double* frame_time);

// ### SELECTORS ###

typedef enum vial_key_t {
    VIAL_KEY_NONE = 0,
    VIAL_KEY_NAME,
    VIAL_KEY_SERIAL,
    VIAL_KEY_ELEMENT,
    VIAL_KEY_KIND,
    VIAL_KEY_TYPE,              // particles: the force field's type, where the data gives one (md_particle_type_ff_type)
    VIAL_KEY_RESID,
    VIAL_KEY_RESNAME,
    VIAL_KEY_CHAIN,
    VIAL_KEY_ID,
    VIAL_KEY_SCOPE,
    VIAL_KEY_SIZE,              // rings: the number of atoms
    VIAL_KEY_ROLE,              // groups: what the simulation used them for (md_group_role_name); particles: their roles
    // A particle's chemistry, of each atom (md_chem_perceive, md_bond_metallic_atoms)
    VIAL_KEY_HYBRIDIZATION,     // sp, sp2 or sp3, where one was perceived
    VIAL_KEY_AROMATIC,          // in an aromatic ring
    VIAL_KEY_METALLIC,          // of a metallic phase: an atom of a metallic bond
    VIAL_KEY_HYDROGENS,         // the hydrogens attached, explicit and implicit
    VIAL_KEY_FORMAL_CHARGE,
    VIAL_KEY_COUNT
} vial_key_t;

typedef enum vial_match_kind_t {
    VIAL_MATCH_STR = 0,         // a string, a glob where the keyword takes one
    VIAL_MATCH_INT,
    VIAL_MATCH_RANGE,           // inclusive
    VIAL_MATCH_BOOL,            // lo: 0 or 1
    VIAL_MATCH_AND,             // of the num operands after it
    VIAL_MATCH_OR,              // of the num operands after it; a list is one
    VIAL_MATCH_NOT,             // of the one operand after it
} vial_match_kind_t;

// One node of what a keyword matches: a value, or and/or/not over the nodes after it
typedef struct vial_match_t {
    uint8_t     kind;           // vial_match_kind_t
    uint32_t    num;            // and, or: the number of operands
    md_str_t    s;
    int64_t     lo, hi;
} vial_match_t;

// A keyword and what it matches, in prefix order: name='CA' is one value; name=['CA', 'CB'] and name='CA' or 'CB' are an or
// of two; name='H*' and not 'HA' an and of a value and a not. Every node stands for the values it accepts: or is their
// union, and their intersection, and not the complement among the values of what has the keyword's property at all, so
// not's complement never reaches what has none (chain=not 'A' is the other chains, not also what is in no chain).
typedef struct vial_filter_t {
    uint8_t         key;        // vial_key_t
    vial_match_t*   expr;       // md_array, prefix order
} vial_filter_t;

// Whether a value satisfies a filter's expression: test says whether one value node holds for it
bool vial_filter_holds(const vial_filter_t* f, bool (*test)(const vial_match_t* m, const void* user), const void* user);

typedef struct vial_select_t {
    uint8_t                 level;          // the selector's: particles is VIAL_LEVEL_PARTICLE
    uint8_t                 particle_kind;  // atom(), bead(): 1 + md_particle_kind_t; 0 for any. polymer(): 1, the
                                            // parts of a polymer kind
    vial_filter_t*          filters;        // md_array
    const vial_sets_t*      within;         // a positional set or array of sets, or NULL
    const vial_sets_t*      scope;          // NULL for the system
} vial_select_t;

// The keywords of a selector, in the order the reference lists them; false for a keyword it does not take
bool        vial_select_takes(uint8_t level, uint8_t key);
const char* vial_select_keywords(uint8_t level);    // "id" for chains, for E0202's message
uint8_t     vial_select_key(md_str_t keyword);         // VIAL_KEY_NONE when it is none of them
const char* vial_select_name(uint8_t level);        // "residues"

// Whether a kind= value is known for the level; the known values for E0206's message
bool        vial_select_kind_known(uint8_t level, md_str_t kind);
const char* vial_select_kinds(uint8_t level);

// Whether a role= value is the name of a role of groups (md_group_role_name)
bool        vial_select_role_known(md_str_t role);
// Whether a value is a role of particles: backbone, main, side, base, sugar, terminus, n_terminus, c_terminus,
// 5_terminus, 3_terminus or qm
bool        vial_select_particle_role_known(md_str_t role);
// Whether a value is a hybridization a particle selector takes: sp, sp2 or sp3
bool        vial_select_hybridization_known(md_str_t hyb);

// Evaluates a selector: one set for particles, one set per unit otherwise, in the system's order
vial_sets_t* vial_select(vial_sys_t* s, const vial_select_t* q, md_allocator_i* alloc);

// The predefined sets: protein, nucleic, water, ion, lipid, backbone, sidechain, system; NULL for another name
vial_sets_t* vial_select_predefined(vial_sys_t* s, md_str_t name, md_allocator_i* alloc);

// For W0602: the chains among the residues a resid filter matches within the scope, up to two
uint32_t vial_select_resid_chains(const vial_sys_t* s, const vial_filter_t* f, const vial_sets_t* scope, int32_t out[2]);

// Glob matching as the name keywords take it, Python's fnmatch: * any run of characters, ? any one, [abc] one of a set,
// with ranges [a-z], and [!abc] one of none of it; [*], [?] and [[] are those characters themselves. Case matters.
bool vial_glob(md_str_t pattern, md_str_t text);
// Whether a pattern has any of that, so it may match more than one text
bool vial_is_glob(md_str_t pattern);

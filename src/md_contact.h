#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitfield_t;
struct md_system_t;
struct md_system_state_t;

// ### CONTACTS ###
// Two layers.
//
// PARTICLE PAIRS (md_contact_pairs_*) finds every pair of particles within a radius of each other, minus the
// pairs the topology excludes. It knows nothing about groups or what a contact means. The pairs are streamed to
// a callback in batches and never collected: in a large system (10^6 particles with a few tens of neighbours
// each) they would not fit, and every consumer reduces them to something much smaller anyway.
//
// GROUP CONTACTS (md_contact_query_*) is one such consumer. It applies a contact criterion to the pairs and
// reduces them to the pairs of groups in contact.
//
// CONTACT REGIONS (md_contact_regions) builds on group contacts: with the groups as units of larger bodies, it
// joins the unit pairs in contact into connected patches of contact between the bodies.
//
// Both split their work by what it depends on: *_init does what depends on the topology only (membership,
// exclusions), the evaluation what depends on the coordinates. Many frames of a trajectory are prepared once
// and evaluated per frame. Distances honour the unit cell of the evaluated state (minimum image).

typedef enum md_contact_flags_t {
    MD_CONTACT_FLAG_NONE = 0,
    MD_CONTACT_FLAG_SELF = 1,           // Within one set: every unordered pair once
} md_contact_flags_t;

// ### PARTICLE PAIRS ###

typedef struct md_contact_pairs_desc_t {
    const struct md_bitfield_t* set_a;

    // Optional. NULL for the pairs within set_a.
    const struct md_bitfield_t* set_b;

    double radius;                      // Ångström. Pairs closer than this are reported.

    // Exclude pairs joined by a path of at most this many bonds, 0 excludes nothing. This is what a force field
    // calls the exclusions of a molecule (nrexcl in GROMACS): 3 excludes 1-2, 1-3 and 1-4 pairs.
    uint32_t exclude_bonds;

    // Optional, one per particle of the system: pairs of particles with the same label are skipped. With the
    // molecule (or structure) as label, only pairs between molecules remain. In a system of packed fibrils most
    // neighbours of a particle are in its own fibril, so this removes most of the stream at its source.
    const uint32_t* particle_label;
} md_contact_pairs_desc_t;

// Prepared particle pair query. Treat as opaque.
typedef struct md_contact_pairs_t {
    uint32_t num_atoms;
    uint32_t flags;
    float    radius;

    int32_t* a_atoms;       // The particles of set_a (and set_b), ascending
    size_t   num_a_atoms;
    int32_t* b_atoms;
    size_t   num_b_atoms;

    // Excluded partners per particle, in compressed rows over all particles, each row sorted. NULL if none.
    uint32_t* excl_off;
    uint32_t* excl_atom;

    uint32_t* label;        // [num_atoms] or NULL

    struct md_allocator_i* alloc;
} md_contact_pairs_t;

// Receives a batch of pairs: particle a[k] and b[k] at distance r[k], for k < count.
//   - Within one set: every unordered pair once, as a < b.
//   - Between two sets: a from set_a and b from set_b. A particle in both sets is never paired with itself, but
//     two particles which are both in both sets are reported in both orders, as each order is a pair of the two sets.
// Batches come in no particular order and have no particular size; the arrays are only valid during the call.
// The stream keeps its own scratch memory apart, so the callback may allocate from any allocator, temp arenas included.
typedef void (*md_contact_pair_callback_t)(const uint32_t* a, const uint32_t* b, const float* r, size_t count, void* user_param);

// ### GROUP CONTACTS ###
// A group is any set of particles (a residue, a chain, a slice, a fibril, a ligand, ...), given as a bitfield;
// groups may overlap. Two groups are in contact when some pair of their particles satisfies the criterion. The
// result is the sparse set of group pairs in contact, with the number of particle pairs behind each and their
// shortest distance.
//
//   - Between two sets of groups A and B: ordered pairs (i, j), group i of A and group j of B.
//   - Within one set of groups A (B omitted): unordered pairs of distinct groups, reported once as i < j.
//     Particle pairs within a single group never count.

typedef enum md_contact_criterion_t {
    MD_CONTACT_CRITERION_DISTANCE  = 0, // r < cutoff
    MD_CONTACT_CRITERION_RADII     = 1, // r < r_a + r_b + cutoff, the cutoff being a tolerance (which may be negative).
                                        // Suits atomistic systems, where a particle has a size of its own.
    MD_CONTACT_CRITERION_TYPE_PAIR = 2, // r < type_cutoff[type_a * num_types + type_b]. For particles whose contact
                                        // distance depends on the pair, not additively on each, as coarse grained beads
                                        // with per pair Lennard-Jones parameters: r_c = lambda 2^(1/6) sigma_ab.
} md_contact_criterion_t;

typedef struct md_contact_desc_t {
    const struct md_bitfield_t* group_a;
    size_t num_a;

    // Optional. NULL for contacts within group_a. Given with num_b == 0 it is an empty set, which has no contacts.
    const struct md_bitfield_t* group_b;
    size_t num_b;

    md_contact_criterion_t criterion;
    double cutoff;                      // Ångström

    // MD_CONTACT_CRITERION_RADII only, optional: one radius per particle of the system. NULL uses the radii of the atom types.
    const float* atom_radius;

    // MD_CONTACT_CRITERION_TYPE_PAIR only: a type per particle of the system, and a symmetric num_types x num_types
    // table of contact distances (Ångström) indexed by type pairs. The types are the caller's: typically the non-bonded
    // types of the force field, which need not map one to one onto the atom types of the system.
    const uint32_t* particle_type;
    uint32_t num_types;
    const float* type_cutoff;

    // See md_contact_pairs_desc_t
    uint32_t exclude_bonds;
    const uint32_t* particle_label;

    // Contacts within one set only: exclude pairs of groups of the same parent whose indices differ by less than
    // this; 0 and 1 exclude nothing. For groups in sequence order along their parent (residues along a chain,
    // slices along a fibril), 3 ignores contacts between neighbours.
    uint32_t min_separation;

    // Optional, one per group of group_a: the parent (chain, fibril, ...) of each group. Groups of different parents
    // are never 'neighbours', whatever their indices. NULL treats all groups as one sequence.
    const uint32_t* group_parent;
} md_contact_desc_t;

// A set of group pairs in contact, sorted by (i, j).
typedef struct md_contact_set_t {
    uint32_t num_a;
    uint32_t num_b;         // Equal to num_a for contacts within one set
    uint32_t flags;         // md_contact_flags_t

    size_t    count;        // Number of group pairs in contact
    uint32_t* i;            // [count] Index of the group in A
    uint32_t* j;            // [count] Index of the group in B, or of the other group in A (then i < j)
    uint32_t* atom_pairs;   // [count] Number of (particle of i, particle of j) pairs which satisfy the criterion.
                            //         For disjoint groups, the number of particle pairs. With overlapping groups a
                            //         particle pair counts once for every way it joins i and j.
    float*    d_min;        // [count] Shortest distance among those pairs

    struct md_allocator_i* alloc;
} md_contact_set_t;

// Prepared group contact query. Treat as opaque.
typedef struct md_contact_query_t {
    md_contact_pairs_t pairs;

    uint32_t num_atoms;
    uint32_t num_a;
    uint32_t num_b;
    uint32_t flags;

    md_contact_criterion_t criterion;
    float    cutoff;
    uint32_t min_separation;
    uint32_t* parent;       // [num_a] or NULL

    // Particle to groups, in compressed rows over all particles of the system: the groups of particle k are
    // grp[off[k]] .. grp[off[k+1]-1]. The B set is NULL for contacts within one set.
    uint32_t* a_off;
    uint32_t* a_grp;
    uint32_t* b_off;
    uint32_t* b_grp;

    float* radius;          // RADII only: per particle

    uint32_t* type;         // TYPE_PAIR only: per particle
    uint32_t  num_types;
    float*    type_cutoff;  // TYPE_PAIR only: num_types^2

    struct md_allocator_i* alloc;
} md_contact_query_t;

// ### CONTACT REGIONS ###
// A set of group contacts says which groups touch. A region says where two larger bodies touch: one connected
// patch of contact between them. The groups are the units of the bodies: residues of chains, monomers of polymers,
// slices of fibrils, lipids of a leaflet, ... A body is given per unit.
//
// Two unit pairs in contact belong to the same region when they join the same two bodies and their units are
// neighbours on both sides: within 'reach' steps of each other in the unit adjacency of their body. For units in
// sequence along two chains, a region is a connected patch of the contact map of the pair: a parallel pairing is a
// diagonal streak, an antiparallel one an anti-diagonal, a crossing a compact blob, and two chains touching twice
// have two regions. Within one body (a chain folding onto itself) regions are patches of its own contact map.
//
// The adjacency says which units neighbour each other. md_contact_unit_adjacency_from_bonds derives it from the
// bonds between units, which depends on no numbering and holds for branched and cyclic bodies alike. Units without
// neighbours form a region per unit pair.

typedef struct md_contact_units_t {
    size_t count;               // Number of units: the groups of one side of the contact set
    const uint32_t* body;       // [count] The body of each unit

    // Optional unit adjacency in compressed rows: the neighbours of unit u are adj[adj_off[u]] .. adj[adj_off[u+1]-1].
    // Neighbours in another body are ignored.
    const uint32_t* adj_off;    // [count + 1]
    const uint32_t* adj;
} md_contact_units_t;

typedef struct md_contact_region_set_t {
    size_t    count;            // Number of regions, ordered by their first unit pair in the contact set
    uint32_t* body_i;           // [count] Body of the A side (within one set: the lower of the two)
    uint32_t* body_j;           // [count] Body of the B side (within one set: the higher, or the same body)
    uint32_t* unit_pairs;       // [count] Number of unit pairs in contact
    uint32_t* units_i;          // [count] Number of distinct units on the i side...
    uint32_t* units_j;          // [count] ... and on the j side. Within one body: the units on either side of the pairs.
    uint32_t* atom_pairs;       // [count] Sum of the particle pairs of the unit pairs
    float*    d_min;            // [count] Shortest particle distance

    size_t    num_pairs;        // The count of the contact set...
    uint32_t* region;           // [num_pairs] ... and the region of each of its unit pairs

    struct md_allocator_i* alloc;
} md_contact_region_set_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### PARTICLE PAIRS ###
bool md_contact_pairs_init(md_contact_pairs_t* pairs, const md_contact_pairs_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_contact_pairs_free(md_contact_pairs_t* pairs);

// Streams the pairs for the coordinates and unit cell of state to callback. Returns false if the state does not
// match the system the pairs were prepared for.
bool md_contact_pairs_for_each(const md_contact_pairs_t* pairs, const struct md_system_state_t* state, md_contact_pair_callback_t callback, void* user_param);

// ### GROUP CONTACTS ###
bool md_contact_query_init(md_contact_query_t* query, const md_contact_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_contact_query_free(md_contact_query_t* query);

// Evaluates the query for the coordinates and unit cell of state. The arrays of out are allocated from alloc.
bool md_contact_query_eval(md_contact_set_t* out, const md_contact_query_t* query, const struct md_system_state_t* state, struct md_allocator_i* alloc);

// Prepare, evaluate and release in one go, for a single evaluation
bool md_contact_compute(md_contact_set_t* out, const md_contact_desc_t* desc, const struct md_system_t* sys, const struct md_system_state_t* state, struct md_allocator_i* alloc);

void md_contact_set_free(md_contact_set_t* set);

// ### CONTACT REGIONS ###

// Unit adjacency from the bonds of the system: units are neighbours when a bond joins a particle of one with a
// particle of the other. A particle in several units counts for the first. Writes compressed rows (see
// md_contact_units_t), each row sorted and unique, allocated from alloc: out_off [num_units + 1], out_adj [out_off[num_units]].
bool md_contact_unit_adjacency_from_bonds(uint32_t** out_off, uint32_t** out_adj, const struct md_bitfield_t* units, size_t num_units, const struct md_system_t* sys, struct md_allocator_i* alloc);

// The regions of a contact set between units. units_b is NULL for a contact set within one set of units.
// reach >= 1 is the number of adjacency steps within which units count as neighbours: 1 joins contacts of adjacent
// units only, 2 bridges a single unit out of contact, and so on.
bool md_contact_regions(md_contact_region_set_t* out, const md_contact_set_t* contacts, const md_contact_units_t* units_a, const md_contact_units_t* units_b, uint32_t reach, struct md_allocator_i* alloc);

void md_contact_region_set_free(md_contact_region_set_t* set);

#ifdef __cplusplus
}
#endif

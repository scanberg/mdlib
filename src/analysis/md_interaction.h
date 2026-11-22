#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitset_t;
struct md_system_t;
struct md_system_state_t;

// ### IONIC AND HALOGEN BONDS ###
//
// The interactions besides hydrogen bonds (md_hbond) that are found from the geometry of a frame, of the kinds of
// md_bond_kind_t: ionic bonds (salt bridges) and halogen bonds. As md_hbond and md_contact, the work is split by what it
// depends on: md_interaction_query_init the topology (charged groups, halogens, acceptors, exclusions), and
// md_interaction_query_eval the coordinates of one state. Distances and angles use the minimum image of its unit cell.
//
// IONIC. Between a cationic and an anionic group, by the formal charges the chemistry perceived (md_chem_perceive, at
// pH 7 for a residue without hydrogens) and the ions' charges. A group is the atoms its charge is spread over: the atoms
// joined by delocalized bonds (MD_BOND_FLAG_DELOCALIZED: a carboxylate, a guanidinium, a phosphate), or a charged atom
// alone (an ammonium N, an ion). Its net charge is the sum of its atoms', and a group whose charges cancel, a nitro, is
// none. Its charged atoms are its N, O and S, or the atom of a group of one. Two groups of opposite charges are in an
// ionic bond when two of their charged atoms are within ionic_max_dist, 4.0 Å as Barlow and Thornton (J. Mol. Biol. 168,
// 867 (1983)) and Kumar and Nussinov (ChemBioChem 3, 604 (2002)) take a salt bridge, and the bond is between the
// closest two: one bond per pair of groups, the cation's atom first. A charge on one atom of a ring, a histidinium's, is
// that atom's alone, as the ring's bonds are aromatic rather than delocalized.
//
// HALOGEN. C-X...A: a Cl, Br or I bonded to one C, X, towards an acceptor A as md_hbond perceives the acceptors (O and N
// with a free lone pair, divalent S). X...A within the sum of their van der Waals radii (Bondi's: 3.27 Å for Cl...O) plus
// halogen_tolerance, and the angle C-X...A at least halogen_min_cxa: 140 degrees, the edge of the sigma hole whose bonds
// centre at 165 (Auffinger, Hays, Westhof and Ho, PNAS 101, 16789 (2004)). X first.
//
// EXCLUSIONS. Two atoms joined by a path of at most exclude_bonds bonds are in no bond of either kind, as the ends of a
// zwitterion of one residue. SELECTION. With a set, the bonds with both atoms in it.

typedef enum md_interaction_kind_t {
    MD_INTERACTION_IONIC   = 1u << 0,
    MD_INTERACTION_HALOGEN = 1u << 1,
} md_interaction_kind_t;

typedef struct md_interaction_params_t {
    uint32_t kinds;             // md_interaction_kind_t, the kinds to find
    uint32_t exclude_bonds;     // 3
    float    ionic_max_dist;    // Å, between charged atoms of the two groups: 4.0
    float    halogen_tolerance; // Å, added to the sum of the van der Waals radii of X and A: 0
    float    halogen_min_cxa;   // degrees, the angle C-X...A at X (180 is linear): 140
} md_interaction_params_t;

typedef struct md_interaction_desc_t {
    const md_interaction_params_t* params;  // NULL for the defaults, every kind
    const struct md_bitset_t* set;          // NULL for every bond of the system
} md_interaction_desc_t;

// Prepared query. Treat as opaque.
typedef struct md_interaction_query_t {
    md_interaction_params_t params;
    uint32_t  num_atoms;

    // The charged atoms of the charged groups
    size_t    num_charged;
    uint32_t* charged_atom;     // [num_charged]
    uint32_t* charged_group;    // [num_charged] the group, a number of its own
    int8_t*   charged_sign;     // [num_charged] +1 of a cation, -1 of an anion
    uint32_t  num_cations, num_anions;  // groups

    // The halogens bonded to one C, and the acceptors
    size_t    num_halogens;
    uint32_t* halogen;          // [num_halogens] X
    uint32_t* halogen_c;        // [num_halogens] its C
    float*    halogen_vdw;      // [num_halogens] Å, its van der Waals radius
    size_t    num_acceptors;
    uint32_t* acceptor;         // [num_acceptors]
    float*    acceptor_vdw;     // [num_acceptors] Å

    // The atoms within exclude_bonds of each charged atom, then of each halogen, as sorted rows:
    // [num_charged + num_halogens + 1] offsets
    uint32_t* excl_off;
    uint32_t* excl_atom;

    uint8_t*  sel;              // [num_atoms] in the set, NULL without one
    struct md_allocator_i* alloc;
} md_interaction_query_t;

// The bonds of a state, sorted by (kind, a, b): the same bond has the same key in every frame
typedef struct md_interaction_set_t {
    size_t    count;
    uint32_t* a;                // [count] the cation's atom, or X
    uint32_t* b;                // [count] the anion's atom, or A
    uint8_t*  kind;             // [count] md_interaction_kind_t, one bit
    float*    dist;             // [count] Å
    struct md_allocator_i* alloc;
} md_interaction_set_t;

#ifdef __cplusplus
extern "C" {
#endif

md_interaction_params_t md_interaction_params_default(void);

bool md_interaction_query_init(md_interaction_query_t* query, const md_interaction_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_interaction_query_free(md_interaction_query_t* query);

// Whether the query finds any bond of the kind at all: groups of both signs, or halogens and acceptors
bool md_interaction_query_can_find(const md_interaction_query_t* query, md_interaction_kind_t kind);

// Evaluates the query for the coordinates and unit cell of state. The arrays of out are allocated from alloc.
bool md_interaction_query_eval(md_interaction_set_t* out, const md_interaction_query_t* query, const struct md_system_state_t* state, struct md_allocator_i* alloc);

void md_interaction_set_free(md_interaction_set_t* set);

#ifdef __cplusplus
}
#endif

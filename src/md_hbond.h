#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitfield_t;
struct md_system_t;
struct md_system_state_t;

// ### HYDROGEN BONDS ###
//
// A hydrogen bond D-H...A joins a donor D (a heavy atom carrying the hydrogen H) to an acceptor A. Which bonds are
// reported is governed by md_hbond_params_t, applied as a pipeline of independent stages in a fixed order:
//
//   1. ROLES        Which atoms donate and accept, and how many bonds an acceptor can take (its free lone pairs).
//                   By default perceived from the chemistry: an N accepts only if its lone pair is free, so amide,
//                   guanidinium, pyrrole type and ammonium N do not. MD_HBOND_ROLES_ALL_N_O makes every N and O an
//                   acceptor, as most analysis tools do.
//   2. GEOMETRY     Independent gates on distances and angles; a gate at 0 is off.
//   3. STRENGTH     Every bond that passes the gates gets a strength in [0, 1], the product of a distance and an
//                   angle term. Each is 1 at the ideal (H...A 1.9 Å or D...A 2.8 Å, and linear D-H...A) and falls
//                   smoothly to 0 at its gate, so a bond fades out before it crosses a gate rather than vanishing
//                   abruptly. Bonds below min_strength are dropped. The strength measures how clearly a bond meets
//                   the criterion in use, not an energy: the same bond is weaker under tighter gates.
//   4. COMPETITION  Optional limits on the number of bonds per hydrogen and per acceptor. Bonds are taken in order of
//                   decreasing strength and a bond is kept while both its hydrogen and its acceptor have room.
//                   Without limits (capacity 0 and UNLIMITED) every bond that passes the gates is reported.
//
// The strength threshold comes before the competition so that weak bonds never take up room.
//
// Like md_contact, the work is split by what it depends on: md_hbond_query_init depends on the topology only (roles,
// exclusions, selections) and md_hbond_query_eval on the coordinates of one state. Prepare once, evaluate per frame.
// Evaluation does not modify the query, so frames may be evaluated concurrently.
//
// SELECTIONS. The reported bonds can be restricted to a set (bonds within it) or to two sets (bonds between them).
// The result does not depend on the selection: competing donors and acceptors outside the selection are taken into
// account (all candidates within three times the search radius of the selected ones, which resolves the competition
// exactly up to two steps away), so a bond is reported for a selection exactly when it is reported for the whole system.
//
// PERIODIC BOUNDARIES. All distances and angles use the minimum image of the unit cell of the evaluated state.
//
// HYDROGENS. Bonds are only found for donors with explicit hydrogens. A system without any hydrogen on its N and O
// atoms (most crystal structures) has no donors; the query reports this with MD_HBOND_FLAG_NO_HYDROGENS.

typedef enum md_hbond_role_t {
    MD_HBOND_ROLE_NONE     = 0,
    MD_HBOND_ROLE_DONOR    = 1,
    MD_HBOND_ROLE_ACCEPTOR = 2,
} md_hbond_role_t;

typedef enum md_hbond_role_flags_t {
    MD_HBOND_ROLES_DEFAULT      = 0,        // N and O, acceptors perceived from the chemistry
    MD_HBOND_ROLES_ALL_N_O      = 1u << 0,  // Every N and O is an acceptor, no perception (gmx hbond, MDTraj, MDAnalysis)
    MD_HBOND_ROLES_SULFUR       = 1u << 1,  // S-H as donor, divalent S as acceptor
    MD_HBOND_ROLES_FLUORINE     = 1u << 2,  // Covalently bound F as acceptor
    MD_HBOND_ROLES_HALIDE_IONS  = 1u << 3,  // F-, Cl-, Br-, I- as acceptors without a capacity limit
} md_hbond_role_flags_t;

typedef enum md_hbond_capacity_mode_t {
    MD_HBOND_CAPACITY_LONE_PAIRS = 0,       // Per acceptor: its free lone pairs (O 2, N 1, S 2, F 3, halide ions unlimited)
    MD_HBOND_CAPACITY_FIXED      = 1,       // The same for every acceptor: acc_capacity_fixed
    MD_HBOND_CAPACITY_UNLIMITED  = 2,
} md_hbond_capacity_mode_t;

// Capacity value for an acceptor without a limit (md_hbond_perceive_roles)
#define MD_HBOND_CAPACITY_NO_LIMIT 0xFF

typedef struct md_hbond_params_t {
    // 1. Roles, md_hbond_role_flags_t
    uint32_t roles;

    // Donor and acceptor joined by a path of at most this many bonds never form a hydrogen bond (3 excludes 1-2, 1-3
    // and 1-4 pairs, such as the N-H...O=C of a single residue). A donor never bonds to itself whatever the value.
    uint32_t exclude_bonds;

    // 2. Geometry gates, 0 is off. At least one of the distance gates must be on.
    float max_ha;           // Å, distance H...A
    float max_da;           // Å, distance D...A
    float min_dha;          // degrees, angle D-H...A at the hydrogen (180 is linear)
    float max_hda;          // degrees, angle H-D...A at the donor (0 is linear), the angle gmx hbond uses
    float min_xah;          // degrees, angle X-A...H at the acceptor for every atom X bonded to A. Rejects hydrogens
                            // approaching an acceptor from behind, away from its lone pairs.

    // 3. Strength
    float min_strength;     // [0, 1], 0 is off

    // 4. Competition
    uint32_t h_capacity;            // Bonds per hydrogen, 0 is unlimited
    float    bifurcation_tol;       // A hydrogen's bonds after its first are kept only with a strength of at least this
                                    // fraction of its strongest, 0 is off. Matters only for h_capacity != 1.
    md_hbond_capacity_mode_t acc_capacity_mode;
    uint32_t acc_capacity_fixed;    // MD_HBOND_CAPACITY_FIXED only, 0 is unlimited
} md_hbond_params_t;

typedef enum md_hbond_preset_t {
    MD_HBOND_PRESET_REALISTIC = 0,  // Perceived roles, Baker-Hubbard geometry, acceptor angle, one bond per hydrogen,
                                    // lone pair capacity: H...A <= 2.5 Å, D-H...A >= 120°, X-A...H >= 90°
    MD_HBOND_PRESET_MDTRAJ,         // Baker-Hubbard as in MDTraj: H...A <= 2.5 Å, D-H...A >= 120°
    MD_HBOND_PRESET_MDANALYSIS,     // MDAnalysis HydrogenBondAnalysis defaults: D...A <= 3.0 Å, D-H...A >= 150°
    MD_HBOND_PRESET_GROMACS,        // gmx hbond defaults: D...A <= 3.5 Å, H-D...A <= 30°
    MD_HBOND_PRESET_VMD,            // VMD defaults: D...A <= 3.0 Å, D-H...A within 20° of linear
    MD_HBOND_PRESET_COUNT,
} md_hbond_preset_t;

// Flags of a query and of the sets it produces
typedef enum md_hbond_flags_t {
    MD_HBOND_FLAG_NONE          = 0,
    MD_HBOND_FLAG_SELECTION     = 1u << 0,  // Restricted to set_a (and set_b)
    MD_HBOND_FLAG_BETWEEN       = 1u << 1,  // Between set_a and set_b, rather than within set_a
    MD_HBOND_FLAG_NO_HYDROGENS  = 1u << 2,  // The system has N or O atoms but none of them carries a hydrogen
} md_hbond_flags_t;

typedef struct md_hbond_desc_t {
    // NULL for MD_HBOND_PRESET_REALISTIC
    const md_hbond_params_t* params;

    // Optional. NULL for every bond of the system. Given alone: the bonds with donor and acceptor both in set_a.
    // Given with set_b: the bonds with the donor in one set and the acceptor in the other. A donor is in a set when
    // its heavy atom or its hydrogen is.
    const struct md_bitfield_t* set_a;
    const struct md_bitfield_t* set_b;

    // Optional overrides of the perceived roles. donors: the heavy atoms which donate through the hydrogens bonded to
    // them. acceptors: the atoms which accept, with their perceived capacity, or that of their element if they were
    // not perceived as acceptors.
    const struct md_bitfield_t* donors;
    const struct md_bitfield_t* acceptors;

    // Optional coordinates for perceiving the roles (whether an N is planar, and thereby has no free lone pair).
    // Without them the atom flags of md_util_system_infer are used where present, and the bond graph otherwise.
    const struct md_system_state_t* reference;
} md_hbond_desc_t;

// Prepared query. Treat as opaque.
typedef struct md_hbond_query_t {
    md_hbond_params_t params;
    uint32_t num_atoms;
    uint32_t flags;             // md_hbond_flags_t

    // Donors, one per D-H pair
    size_t    num_donors;
    uint32_t* donor_d;          // [num_donors]
    uint32_t* donor_h;          // [num_donors]
    uint32_t* excl_off;         // [num_donors + 1] Acceptors within exclude_bonds of the donor, sorted rows, or NULL
    uint32_t* excl_atom;

    // Acceptors
    size_t    num_acceptors;
    uint32_t* acceptor;         // [num_acceptors]
    uint8_t*  acceptor_cap;     // [num_acceptors] MD_HBOND_CAPACITY_NO_LIMIT for none
    uint32_t* acc_nbr_off;      // [num_acceptors + 1] Atoms bonded to each acceptor
    uint32_t* acc_nbr;

    uint8_t*  sel;              // [num_atoms] bit 0: in set_a, bit 1: in set_b. NULL without a selection.

    struct md_allocator_i* alloc;
} md_hbond_query_t;

// A set of hydrogen bonds, sorted by (hydrogen, acceptor): the same bond has the same key in every frame.
typedef struct md_hbond_set_t {
    size_t    count;
    uint32_t* donor;            // [count] Heavy atom of the donor
    uint32_t* hydrogen;         // [count]
    uint32_t* acceptor;         // [count]
    float*    strength;         // [count] [0, 1]
    float*    dist_da;          // [count] Å
    float*    dist_ha;          // [count] Å
    float*    angle_dha;        // [count] degrees, D-H...A

    md_hbond_params_t params;   // The parameters which produced the set
    uint32_t flags;             // md_hbond_flags_t of the query

    struct md_allocator_i* alloc;
} md_hbond_set_t;

#ifdef __cplusplus
extern "C" {
#endif

md_hbond_params_t md_hbond_params_preset(md_hbond_preset_t preset);
const char*       md_hbond_preset_name(md_hbond_preset_t preset);

// Roles of the atoms of a system: out_role [num_atoms] md_hbond_role_t bits, out_capacity [num_atoms] (optional) the
// number of bonds an acceptor can take, MD_HBOND_CAPACITY_NO_LIMIT for none, 0 for atoms which do not accept.
// reference is optional, see md_hbond_desc_t.
bool md_hbond_perceive_roles(uint8_t* out_role, uint8_t* out_capacity, const struct md_system_t* sys, const struct md_system_state_t* reference, uint32_t role_flags);

// Sets MD_FLAG_HBOND_DONOR and MD_FLAG_HBOND_ACCEPTOR on the atoms of the system from the default roles
// (MD_HBOND_ROLES_DEFAULT | MD_HBOND_ROLES_HALIDE_IONS). Called by md_util_system_infer for MD_UTIL_INFER_HBOND_BIT.
void md_hbond_infer_atom_flags(struct md_system_t* sys, const struct md_system_state_t* reference);

bool md_hbond_query_init(md_hbond_query_t* query, const md_hbond_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_hbond_query_free(md_hbond_query_t* query);

// Evaluates the query for the coordinates and unit cell of state. The arrays of out are allocated from alloc.
bool md_hbond_query_eval(md_hbond_set_t* out, const md_hbond_query_t* query, const struct md_system_state_t* state, struct md_allocator_i* alloc);

// Prepare, evaluate and release in one go, for a single evaluation
bool md_hbond_compute(md_hbond_set_t* out, const md_hbond_desc_t* desc, const struct md_system_t* sys, const struct md_system_state_t* state, struct md_allocator_i* alloc);

void md_hbond_set_free(md_hbond_set_t* set);

#ifdef __cplusplus
}
#endif

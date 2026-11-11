#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <core/md_array.h>
#include <md_atom.h>
#include <core/md_bitset.h>

struct md_system_t;
struct md_system_state_t;
struct md_executor_t;

// Bonds: the bond graph of the atoms and the non-covalent bonds (md_bond_data_t), and iterating them.

// Bond flags: what a bond is (the CHEMISTRY in the low byte: order, resonance and kind) and where it and its chemistry
// came from (the PROVENANCE above it). Each field is one value or one bit, so every combination describes a case (bar
// the unused values of a field). Only the low byte is part of a bond's identity when comparing structures.
//
//   bit   11         10               9-8     7-5    4-3        2-0
//         KIND_GIVEN ORDER_PERCEIVED  ORIGIN  KIND   RESONANCE  ORDER
//
// CHEMISTRY
//   ORDER      A value in bits 0-2: read it with md_bond_order and write it with md_bond_flags_set_order. 0 is unknown
//              (not perceived, and not given by the file), which is not the same as single. Always 0 for a bond outside
//              of the bond graph.
//   RESONANCE  A value in bits 3-4, at most one of:
//     AROMATIC     In an aromatic ring.
//     DELOCALIZED  Resonance outside of rings: carboxylate, guanidinium, nitro, phosphate, ...
//              Either comes on top of an order which holds one Kekule structure: a benzene bond is AROMATIC with order 1
//              or 2, the two C-O of a carboxylate are DELOCALIZED with orders 2 and 1. Neither is localized. Always
//              localized for a bond outside of the bond graph.
//   KIND       A value in bits 5-7 (md_bond_kind_t): read it with md_bond_kind and write it with md_bond_flags_set_kind.
//              Bit 7 set: outside of the bond graph (md_bond_is_graph). Bits 7 and 6 set: non-covalent, of the second
//              range of md_bond_data_t (md_bond_is_noncovalent).
//     COVALENT     0x00  A bond of the molecule, also between two metals of a molecule (the Mn-Mn of Mn2(CO)10, the
//                        Mo-Mo of a paddlewheel): whether a bond involves a metal is a question for the elements of its
//                        atoms. Two metal atoms of a metallic phase have a METALLIC bond instead.
//     COORDINATE   0x20  A metal-ligand bond, which counts towards the valence of neither atom: every bond between a
//                        metal and a non-metal unless its kind was given (md_system_infer_coordination), and any
//                        other bond a file says is one. Whether the ligand is then an anion (the thiolate of a zinc
//                        finger cysteine) is for the perception of its atom to say (md_chem_perceive), not the bond.
//     ELASTIC      0x80  Not a bond of the molecule but a spring of its force field which holds two particles at a
//                        distance: the elastic network of a coarse grained protein or nucleic acid (the rubber bands of
//                        martinize, ElNeDyn), the restraints of its secondary structure, a restraint between two atoms.
//                        GROMACS writes them as bonds which make no connection (bond types 6, 9 and 10, constraints of
//                        type 2), or as plain bonds (md_system_infer_elastic_bonds). They are kept, as
//                        information, in the first range of md_bond_data_t, and left out of everything which follows
//                        the bonds of a molecule: the connectivity (md_bond_conn_data_t, so md_bond_iter and
//                        md_bond_find), the molecules, rings, chemistry and matching. The representations do not draw
//                        them unless asked to.
//     METALLIC     0xA0  Between two atoms of a metallic phase (md_bond_metallic_atoms: a metal in bulk, at a surface,
//                        in a particle), neighbours within its first coordination shell. The bonding of a metal is not
//                        that of a molecule, so like the elastic bonds they are kept in the first range of
//                        md_bond_data_t and left out of everything which follows the bonds of a molecule (the
//                        connectivity, molecules, rings, chemistry, matching, the bonds formed and broken over a run).
//                        Licorice and ball and stick do not draw them; a Bonds representation does (bond(kind=
//                        "metallic")). Inferred only (METALLIC PHASES in md_bond.c).
//     INTERACTION  0xC0  A non-covalent bond of no particular kind: one a user marks between two atoms, a salt bridge,
//                        the contact of a halogen with its acceptor.
//     HYDROGEN     0xE0  A hydrogen bond, D-H...A. Stored with the hydrogen (or the donor without one) as the first
//                        atom.
//              COVALENT and COORDINATE bonds make up the bond graph, from which the connectivity, the rings, the
//              molecules and the perceived chemistry are derived. INTERACTION and HYDROGEN are the non-covalent bonds,
//              which md_bond_data_t keeps in a range of their own. 0x40 and 0x60 are unused.
//
// PROVENANCE
//   ORIGIN           A value in bits 8-9 (md_bond_origin_t): read it with md_bond_origin and write it with
//                    md_bond_flags_set_origin.
//   ORDER_PERCEIVED  The order and the resonance were perceived (md_chem_perceive), not given. A known order without
//                    it was given, by the file or by hand, and stays as it is when the chemistry is perceived again;
//                    one with it is perceived anew. A bond without an order is perceived whether or not it is set, so
//                    a loader which gives orders gets them kept without saying so.
//   KIND_GIVEN       The kind was given, by the file or by hand, and stays as it is when bonds are perceived again: a
//                    bond between a metal and a non-metal a user says is covalent stays covalent. The other way around
//                    from ORDER_PERCEIVED, as a kind has no 'unknown' to tell a given one by: a bond without it has the
//                    kind perception gives it.
typedef enum {
    MD_BOND_FLAG_NONE            = 0,
    MD_BOND_FLAG_ORDER_MASK      = 0x7,     // See md_bond_order
    MD_BOND_FLAG_AROMATIC        = 0x8,     // Resonance: in an aromatic ring
    MD_BOND_FLAG_DELOCALIZED     = 0x10,    // Resonance: outside of rings
    MD_BOND_FLAG_RESONANCE_MASK  = 0x18,
    MD_BOND_FLAG_KIND_MASK       = 0xE0,    // See md_bond_kind
    MD_BOND_FLAG_ORIGIN_MASK     = 0x300,   // See md_bond_origin
    MD_BOND_FLAG_ORDER_PERCEIVED = 0x400,   // The order and the resonance were perceived, not given
    MD_BOND_FLAG_KIND_GIVEN      = 0x800,   // The kind was given, not perceived
} md_bond_flags_t;

MD_ENUM_FLAGS(md_bond_flags_t)

enum {
    MD_BOND_ORDER_UNKNOWN   = 0,
    MD_BOND_ORDER_SINGLE    = 1,
    MD_BOND_ORDER_DOUBLE    = 2,
    MD_BOND_ORDER_TRIPLE    = 3,
    MD_BOND_ORDER_QUADRUPLE = 4,
};

// Where a bond came from. FILE is 0, so a bond added without saying where it came from counts as given.
typedef enum md_bond_origin_t {
    MD_BOND_ORIGIN_FILE     = 0,    // Given by the file, for some of its atoms (CONECT records, a residue's bonds)
    MD_BOND_ORIGIN_TOPOLOGY = 1,    // A force field topology: all of the bonds of the atoms it covers
    MD_BOND_ORIGIN_USER     = 2,    // Added by hand
    MD_BOND_ORIGIN_INFERRED = 3,    // Inferred from the geometry, and replaced when inferred again
} md_bond_origin_t;

#define MD_BOND_ORDER_SHIFT  0
#define MD_BOND_ORIGIN_SHIFT 8

static inline int md_bond_order(md_bond_flags_t flags) {
    return (int)((flags & MD_BOND_FLAG_ORDER_MASK) >> MD_BOND_ORDER_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_order(md_bond_flags_t flags, int order) {
    const unsigned int value = (order < 0 ? 0u : (order > 4 ? 4u : (unsigned int)order)) << MD_BOND_ORDER_SHIFT;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORDER_MASK) | value);
}

static inline md_bond_origin_t md_bond_origin(md_bond_flags_t flags) {
    return (md_bond_origin_t)(((unsigned int)flags & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK) >> MD_BOND_ORIGIN_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_origin(md_bond_flags_t flags, md_bond_origin_t origin) {
    const unsigned int value = ((unsigned int)origin << MD_BOND_ORIGIN_SHIFT) & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORIGIN_MASK) | value);
}

// What a bond is, the KIND of md_bond_flags_t, in place (not shifted down): see there.
typedef enum md_bond_kind_t {
    MD_BOND_KIND_COVALENT    = 0x00,
    MD_BOND_KIND_COORDINATE  = 0x20,
    MD_BOND_KIND_ELASTIC     = 0x80,
    MD_BOND_KIND_METALLIC    = 0xA0,
    MD_BOND_KIND_INTERACTION = 0xC0,
    MD_BOND_KIND_HYDROGEN    = 0xE0,
} md_bond_kind_t;

static inline md_bond_kind_t md_bond_kind(md_bond_flags_t flags) {
    return (md_bond_kind_t)((unsigned int)flags & (unsigned int)MD_BOND_FLAG_KIND_MASK);
}

static inline md_bond_flags_t md_bond_flags_set_kind(md_bond_flags_t flags, md_bond_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_BOND_FLAG_KIND_MASK;
    return (md_bond_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

// A bond of the bond graph: covalent or coordinate. Elastic, metallic and non-covalent bonds are not.
static inline bool md_bond_is_graph(md_bond_flags_t flags) {
    return ((unsigned int)flags & 0x80u) == 0;
}

// A non-covalent bond, of the second range of md_bond_data_t: an interaction or a hydrogen bond.
static inline bool md_bond_is_noncovalent(md_bond_flags_t flags) {
    return ((unsigned int)flags & 0xC0u) == 0xC0u;
}

#ifdef __cplusplus
extern "C" {
#endif

const char* md_bond_origin_name(md_bond_origin_t origin);   // "" for a value out of range

// The kinds of bonds by name, as the script selects them and a workspace stores them: a bond is of a kind when
// (flags & mask) == value. Each kind of md_bond_kind_t has a name ('covalent', 'coordinate', 'elastic', 'interaction',
// 'hydrogen'), and 'noncovalent' is that of the two non-covalent kinds together.
typedef struct md_bond_kind_name_t {
    const char*     name;
    md_bond_flags_t mask;
    md_bond_flags_t value;
} md_bond_kind_name_t;

// The table, in the order above. count receives its length.
const md_bond_kind_name_t* md_bond_kind_names(size_t* count);

#ifdef __cplusplus
}
#endif

// Atom centric representation of bonds: the graph of the molecules, every bond but the elastic ones
// (MD_BOND_KIND_ELASTIC), which are springs of the force field and connect nothing. Whatever follows bonds through it
// (md_bond_iter, md_bond_find, the molecules, rings and chemistry) leaves them out.
typedef struct md_bond_conn_data_t {
    size_t count;
    md_atom_idx_t* atom_idx; // Indices to the 'other' atoms
    md_bond_idx_t* bond_idx; // Indices to the bonds
    // The offsets into the atom_idx and bond_idx for each atom.
    // Consequently offset_count should be atom count + 1
    size_t offset_count;
    uint32_t* offset;
} md_bond_conn_data_t;

// Bond centric representation, in two ranges of the same arrays:
//   [0, count)                              the bond graph, covalent and coordinate bonds, and the bonds outside of it
//                                           which the connectivity leaves out: the elastic bonds (MD_BOND_KIND_ELASTIC),
//                                           springs of the force field, and the metallic bonds (MD_BOND_KIND_METALLIC)
//   [count, count + noncovalent_count)      the persistent non-covalent bonds (md_bond_is_noncovalent): hydrogen
//                                           bonds and interactions a file gives or a user marks
// Whatever walks the bonds as the graph (the connectivity, rings, molecules, perception, matching, bonds formed and
// broken over a run) reads the first range only, by count, and leaves the elastic and metallic bonds out. The arrays
// hold both ranges, so they are uploaded as one. A bond index names a bond in either range: md_bond_insert,
// md_bond_remove, md_bond_pair and md_system_bond_flags take both, and keep every bond in the range of its kind.
// md_bond_find finds a bond of the graph, md_bond_find_noncovalent one of the second range. Bonds computed for a frame
// (bonds formed, hydrogen bonds found) are never here: they are a bond set of their own (md_bond_delta_t, the script).
typedef struct md_bond_data_t {
    size_t count;                   // Bonds of the first range: the graph, the elastic and the metallic bonds
    size_t noncovalent_count;       // Non-covalent bonds, the second range
    md_atom_pair_t*  pairs;         // [count + noncovalent_count]
    md_bond_flags_t* flags;         // [count + noncovalent_count]
    md_bond_conn_data_t  conn;      // Connectivity of the graph, without the elastic bonds
} md_bond_data_t;

typedef struct md_bond_iter_t {
    const md_bond_data_t* data;
	uint32_t i;
    uint32_t end_idx;
} md_bond_iter_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline md_bond_iter_t md_bond_iter(const md_bond_data_t* bond_data, size_t atom_idx) {
    md_bond_iter_t it = {0};
    if (bond_data && bond_data->conn.offset && atom_idx < bond_data->conn.offset_count) {
        it.data = bond_data;
        it.i = bond_data->conn.offset[atom_idx];
		it.end_idx = bond_data->conn.offset[atom_idx + 1];
    }
    return it;
}

// This is not something which should be done frequently. The elastic and metallic bonds (md_bond_is_graph) are left out.
void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc);

// Frees the arrays of the bond data, allocated with alloc, and empties it
void md_bond_data_free(md_bond_data_t* bond, md_allocator_i* alloc);

static inline size_t md_bond_conn_count(const md_bond_data_t* bond_data, size_t atom_idx) {
    MD_ASSERT(bond_data);
    return bond_data->conn.offset[atom_idx + 1] - bond_data->conn.offset[atom_idx];
}

static inline md_atom_idx_t md_bond_conn_atom_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    MD_ASSERT(bond_data);
    return bond_data->conn.atom_idx[atom_conn_idx + idx];
}

static inline md_bond_idx_t md_bond_conn_bond_idx(const md_bond_data_t* bond_data, uint32_t atom_conn_idx, uint32_t idx) {
    MD_ASSERT(bond_data);
    return bond_data->conn.bond_idx[atom_conn_idx + idx];
}

static inline bool md_bond_iter_has_next(const md_bond_iter_t* it) {
    MD_ASSERT(it);
    return it->i < it->end_idx;
}

static inline void md_bond_iter_next(md_bond_iter_t* it) {
    MD_ASSERT(it);
    it->i += 1;
}

static inline md_atom_idx_t md_bond_iter_atom_index(const md_bond_iter_t* it) {
    MD_ASSERT(it);
	return it->data->conn.atom_idx[it->i];
}

static inline md_atom_idx_t md_bond_iter_bond_index(const md_bond_iter_t* it) {
    MD_ASSERT(it);
    return it->data->conn.bond_idx[it->i];
}

static inline uint32_t md_bond_iter_bond_flags(const md_bond_iter_t* it) {
    MD_ASSERT(it);
    return it->data->flags[it->data->conn.bond_idx[it->i]];
}

// The bond of the graph between two atoms (covalent or coordinate), or -1: whether the two are bonded. Found through
// the connectivity, so neither an elastic bond nor a bond inserted since it was last built is found.
static inline md_bond_idx_t md_bond_find(const md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    MD_ASSERT(bond_data);
	md_bond_iter_t it = md_bond_iter(bond_data, atom_idx_a);
    while (md_bond_iter_has_next(&it)) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        if (other_atom_idx == atom_idx_b) {
            return md_bond_iter_bond_index(&it);
        }
        md_bond_iter_next(&it);
	}
	return -1;
}

// The non-covalent bond between two atoms (an index in the second range of md_bond_data_t), or -1. An atom pair has
// one bond at most, of the graph (md_bond_find) or non-covalent.
static inline md_bond_idx_t md_bond_find_noncovalent(const md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b) {
    MD_ASSERT(bond_data);
    // They are few and have no connectivity
    const size_t end = bond_data->count + bond_data->noncovalent_count;
    for (size_t i = bond_data->count; i < end; ++i) {
        const md_atom_pair_t p = bond_data->pairs[i];
        if ((p.idx[0] == atom_idx_a && p.idx[1] == atom_idx_b) || (p.idx[0] == atom_idx_b && p.idx[1] == atom_idx_a)) {
            return (md_bond_idx_t)i;
        }
    }
	return -1;
}

// Adds a bond between two atoms unless they have one, into the range of its kind: the bond graph, or the non-covalent
// bonds (md_bond_is_noncovalent). A bond of the graph inserted where there are non-covalent bonds moves the first of
// them to the end, which changes its index.
static inline void md_bond_insert(md_bond_data_t* bond_data, md_atom_idx_t atom_idx_a, md_atom_idx_t atom_idx_b, md_bond_flags_t flags, md_allocator_i* alloc) {
    MD_ASSERT(bond_data);
    MD_ASSERT(alloc);
    // A bond outside of the graph has neither an order nor resonance
    MD_ASSERT(md_bond_is_graph(flags) || !(flags & (MD_BOND_FLAG_ORDER_MASK | MD_BOND_FLAG_RESONANCE_MASK)));

    // An atom pair has one bond at most, but for an elastic bond: it is not found, and does not stand in for a bond of the
    // molecule
    if (md_bond_find(bond_data, atom_idx_a, atom_idx_b) != -1 || md_bond_find_noncovalent(bond_data, atom_idx_a, atom_idx_b) != -1) {
        return;
    }

    // Add bond
    md_atom_pair_t pair = {atom_idx_a, atom_idx_b};
    md_array_push(bond_data->pairs, pair, alloc);
    md_array_push(bond_data->flags, flags, alloc);
    if (md_bond_is_noncovalent(flags)) {
        bond_data->noncovalent_count++;
        return;
    }
    if (bond_data->noncovalent_count) {
        // Into the place of the first non-covalent bond, which goes to the end
        const size_t first = bond_data->count;
        const size_t last  = bond_data->count + bond_data->noncovalent_count;
        bond_data->pairs[last] = bond_data->pairs[first];
        bond_data->flags[last] = bond_data->flags[first];
        bond_data->pairs[first] = pair;
        bond_data->flags[first] = flags;
    }
    bond_data->count++;
}

// Removes a bond of either range. The last bond of its range takes its place, so that one changes its index (and so
// does the last non-covalent bond, when a bond of the graph is removed). A bond of the graph removed requires a rebuild
// of the connectivity to be valid again.
static inline void md_bond_remove(md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    MD_ASSERT(bond_data);
    const size_t total = bond_data->count + bond_data->noncovalent_count;
    MD_ASSERT(bond_idx >= 0 && (size_t)bond_idx < total);

    if ((size_t)bond_idx < bond_data->count) {
        // The last bond of the graph takes its place, and the last non-covalent bond that of the last of the graph
        const size_t last = bond_data->count - 1;
        bond_data->pairs[bond_idx] = bond_data->pairs[last];
        bond_data->flags[bond_idx] = bond_data->flags[last];
        if (bond_data->noncovalent_count) {
            bond_data->pairs[last] = bond_data->pairs[total - 1];
            bond_data->flags[last] = bond_data->flags[total - 1];
        }
        bond_data->count--;
    } else {
        bond_data->pairs[bond_idx] = bond_data->pairs[total - 1];
        bond_data->flags[bond_idx] = bond_data->flags[total - 1];
        bond_data->noncovalent_count--;
    }
    md_array_shrink(bond_data->pairs, total - 1);
    md_array_shrink(bond_data->flags, total - 1);
}

static inline md_atom_pair_t md_bond_pair(const md_bond_data_t* bond_data, md_bond_idx_t bond_idx) {
    MD_ASSERT(bond_data);
    if (bond_idx >= 0 && (size_t)bond_idx < bond_data->count + bond_data->noncovalent_count) {
        return bond_data->pairs[bond_idx];
    }
    md_atom_pair_t invalid_pair = { -1, -1 };
    return invalid_pair;
}

static inline void md_bond_conn_clear(md_bond_conn_data_t* conn_data) {
    MD_ASSERT(conn_data);
    conn_data->count = 0;
    md_array_shrink(conn_data->atom_idx, 0);
    md_array_shrink(conn_data->bond_idx, 0);

    conn_data->offset_count = 0;
    md_array_shrink(conn_data->offset, 0);
}

static inline void md_bond_data_clear(md_bond_data_t* bond_data) {
    MD_ASSERT(bond_data);

    bond_data->count = 0;
    bond_data->noncovalent_count = 0;
    md_array_shrink(bond_data->pairs, 0);
    md_array_shrink(bond_data->flags, 0);
    
    bond_data->conn.count = 0;
    md_array_shrink(bond_data->conn.atom_idx, 0);
    md_array_shrink(bond_data->conn.bond_idx, 0);

    md_bond_conn_clear(&bond_data->conn);
}

static inline bool md_atom_is_connected_to_atomic_numbers(const md_atom_data_t* atom_data, const md_bond_data_t* bond_data, size_t atom_idx, const md_atomic_number_t z_list[], size_t z_count) {
    MD_ASSERT(bond_data);
    MD_ASSERT(atom_data);
    MD_ASSERT(atom_idx < atom_data->count);
    bool found = false;
    md_bond_iter_t it = md_bond_iter(bond_data, atom_idx);
    while (md_bond_iter_has_next(&it) && !found) {
        md_atom_idx_t other_atom_idx = md_bond_iter_atom_index(&it);
        md_atomic_number_t other_z = md_atom_atomic_number(atom_data, other_atom_idx);
        for (size_t i = 0; i < z_count; ++i) {
            if (other_z == z_list[i]) {
                found = true;
                break;
            }
        }
        md_bond_iter_next(&it);
    }
    return found;
}

// ### INFERENCE ###
// What is written is out, allocated with alloc; the system is read and never written. out is not the system's own
// (&sys->bond): md_util_system_infer gives the system the result, or the caller does (md_bond_data_free the old, assign).

// The bonds of the atoms of the system from the coordinates of state alone, every one of origin
// MD_BOND_ORIGIN_INFERRED, with their connectivity: the bonds the system has are not looked at (md_bond_merge puts the
// two together). Atoms are bonded by their covalent radii (derived from the element). The atoms of a metallic phase (a
// metal in bulk, at a surface or in a particle) are bonded to the metal atoms of their first coordination shell by
// metallic bonds and to nothing else, and a halogen bonded to a non-metal coordinates no metal (METALLIC PHASES in
// md_bond.c). Coarse grained systems (md_system_is_coarse_grained) are bonded as their Martini topologies are: a residue
// whose beads are those of a residue of the Martini force fields has the bonds of its topology, any other amino acid or
// nucleotide a minimum spanning tree of its beads with its rings. The backbone bead of each residue is bonded to that of
// the next unless the numbering starts over, and cysteines are bridged (infer_coarse_grained_bonds in md_martini.c). The
// pairs within reach are found on the threads of the executor (NULL: the calling thread), and the bonds come out the
// same whatever the threads. False, and out empty, when state has no coordinates.
bool md_bond_infer(md_bond_data_t* out, const struct md_system_t* sys, const struct md_system_state_t* state, md_allocator_i* alloc, const struct md_executor_t* executor);

// The bonds given (those a system has) with the inferred ones (md_bond_infer) in place of those inferred before, with
// their connectivity over atom_count atoms. Every bond of given not of origin MD_BOND_ORIGIN_INFERRED is kept, in its
// order and after the inferred ones, and so are the non-covalent bonds. An inferred bond is left out where it would
// duplicate a pair of given (a non-covalent bond included: a pair has one bond, and that one was given), and between two
// atoms a topology covers (bonds of origin MD_BOND_ORIGIN_TOPOLOGY, which is authoritative where it has bonds). An
// elastic bond does neither: it is no bond of the molecule. out is neither given nor inferred.
void md_bond_merge(md_bond_data_t* out, const md_bond_data_t* given, const md_bond_data_t* inferred, size_t atom_count, md_allocator_i* alloc);

// The atoms of the metallic phases: those of the metallic bonds (MD_BOND_KIND_METALLIC, METALLIC PHASES in md_bond.c).
// out[i] is 1 for such an atom and 0 for any other, of atom_count atoms. Returns how many there are.
size_t md_bond_metallic_atoms(uint8_t* out, const md_bond_data_t* bond, size_t atom_count);

// Marks every bond between a metal and a non-metal MD_BOND_KIND_COORDINATE, whatever its origin, unless its kind was
// given (MD_BOND_FLAG_KIND_GIVEN), and returns how many it marked. md_util_system_infer does it for each system; a bond
// added afterwards (by hand) wants it too.
size_t md_system_infer_coordination(struct md_system_t* sys);

// Marks the bonds of a coarse grained topology which hold a protein or nucleic acid in shape without being bonds of it
// MD_BOND_KIND_ELASTIC, and returns how many it marked: an elastic network (ElNeDyn, the rubber bands of martinize) and
// the restraints of secondary structure, when they are written as plain bonds (bond type 1), as they are as often as
// not. Between two residues of a chain (amino acids or nucleotides, beads) a Martini topology has the link from the last
// backbone bead of one to the first of the next (BB to BB, BB3 to BB1, and BB2 to BB1 as well in Martini 3 RNA), a
// disulfide (the SC1 of two cysteines), and the bond which closes a cyclic chain; every other bond between beads of
// two residues of a chain is elastic. Bonds within a residue, and to anything other than amino acids and nucleotides,
// are left as they are. Needs the kinds of the residues and the backbone beads (md_util_system_infer_residue_flags),
// and rebuilds the connectivity when it marks any. The topology readers (md_tpr, md_itp) call it.
size_t md_system_infer_elastic_bonds(struct md_system_t* sys);

// ### BONDS FORMED AND BROKEN ###
// The bonds of a state against the bond graph of the system, which is the reference topology: the bonds of the graph
// broken in the state, and the bonds formed in it which the system does not have. What the bonds of the state are is
// the graph without the broken bonds and with the formed ones; the system is not changed. Only the bonds of a scope are
// looked at, those with an atom in it (so a bond to an atom just outside of it forms too): outside of it the bonds of
// the state are those of the graph.
//
// The rules are those of md_bond_infer, with a band around their cutoffs so that a bond does not flicker
// from frame to frame as it vibrates near one: a bond of the graph breaks only beyond (1 + tolerance) times its cutoff,
// and a pair forms a bond only within (1 - tolerance) times it. Measured against the graph rather than an earlier state,
// so the bonds of a state are those of that state alone, whatever was evaluated before it.
//   - Only the bonds of the breakable origins break. A bond a user added stays, and by default so does one of a force
//     field topology, which holds throughout a classical simulation.
//   - A formed bond is limited by what its atoms can take (md_bond_infer prunes the same way). A bond of
//     the graph counts against that while it is within its cutoff, but not once it is stretched beyond it on its way to
//     breaking: a bond forming as another breaks (a proton moving, the transition state of a substitution) shows both
//     for the frames it takes, while two atoms which clash in the reference do not bond in every frame.
//   - A pair with a non-covalent bond in the system forms no bond: what was given of it stands.
//   - The metallic bonds neither break nor form: two atoms of a metallic phase (atoms of metallic bonds, md_bond_metallic_atoms) are
//     left out.
//     What touches the phase has its bonds to it here only, formed in the frames where it is close enough, as it
//     moves across a surface.
//   - A halogen bonded to a non-metal in the state, by a bond of the graph which is not broken or by one formed,
//     forms no coordination bond.
//   - A coarse grained system has no bonds formed or broken: its bonds come from templates, not distances.
typedef struct md_bond_delta_t {
    md_bitset_t      broken;        // Bonds of the graph broken in the state: indices into sys->bond, in [0, count)
    size_t           num_formed;
    md_atom_pair_t*  formed_pair;   // [num_formed] Pairs bonded in the state which have no bond in the system, each the
                                    // lower atom first, sorted
    md_bond_flags_t* formed_flags;  // [num_formed] Covalent or coordinate, of order unknown and origin inferred
} md_bond_delta_t;

#define MD_BOND_DELTA_TOLERANCE 0.1f
#define MD_BOND_DELTA_BREAKABLE ((1u << MD_BOND_ORIGIN_FILE) | (1u << MD_BOND_ORIGIN_INFERRED))

typedef struct md_bond_delta_desc_t {
    const md_bitset_t* scope;       // The atoms whose bonds are looked at, NULL for all of them
    float    tolerance;             // Half the width of the band around the cutoffs, as a fraction of them: [0, 0.5).
                                    // MD_BOND_DELTA_TOLERANCE is what to use without a reason for another
    uint32_t breakable;             // The origins whose bonds may break, a mask of 1 << md_bond_origin_t:
                                    // MD_BOND_DELTA_BREAKABLE (a file's and those inferred) without a reason for another
} md_bond_delta_desc_t;

// The delta of a state, written to out and allocated with alloc. False, with an empty delta, when the state has no
// coordinates (or not those of the system) or the system is coarse grained.
bool md_bond_delta_compute(md_bond_delta_t* out, const struct md_system_t* sys, const struct md_system_state_t* state, const md_bond_delta_desc_t* desc, md_allocator_i* alloc);
void md_bond_delta_free(md_bond_delta_t* delta, md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

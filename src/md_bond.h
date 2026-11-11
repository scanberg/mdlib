#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>
#include <md_atom.h>

// Bonds: the bond graph of the atoms and the non-covalent bonds (md_bond_data_t), and iterating them.

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

#ifdef __cplusplus
}
#endif

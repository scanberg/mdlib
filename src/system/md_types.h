#pragma once

#include <stdint.h>

// The plain types of a system: its indices, ranges and atom pairs. What a system stores is declared with what
// stores it: md_atom.h, md_hierarchy.h, md_bond.h, md_backbone.h, md_lipid.h, md_ring.h, the elements in md_element.h
// and the standard monomers in md_monomer.h.

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_residue_idx_t;
typedef int32_t     md_instance_idx_t;
typedef int32_t     md_chain_idx_t;
typedef int32_t     md_part_idx_t;
typedef int32_t     md_part_type_idx_t;
typedef int32_t     md_molecule_idx_t;
typedef int32_t     md_molecule_type_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_sequence_id_t;
typedef int32_t     md_bond_idx_t;
typedef uint16_t    md_particle_type_idx_t;
typedef int32_t     md_residue_type_idx_t;
typedef uint8_t     md_atomic_number_t;      // Atomic number (1-118), 0 is unknown

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_irange_t {
    int32_t beg;
    int32_t end;
} md_irange_t;

typedef struct md_urange_t {
    uint32_t beg;
    uint32_t end;
} md_urange_t;

typedef struct md_atom_pair_t {
    md_atom_idx_t idx[2];
} md_atom_pair_t;

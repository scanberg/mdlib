#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef struct md_molecule_o md_molecule_o; // Opaque data blob

typedef int32_t                     md_atom_idx_t;
typedef int32_t                     md_residue_idx_t;
typedef int32_t                     md_backbone_idx_t;
typedef int32_t                     md_residue_id_t;
typedef int32_t                     md_chain_idx_t;
typedef int32_t                     md_molecule_idx_t;
typedef uint32_t                    md_secondary_structure_t;
typedef uint8_t                     md_flag_t;
typedef uint8_t                     md_element_t;
typedef uint8_t                     md_ramachandran_type_t;

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN    = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL    = 1,
    MD_RAMACHANDRAN_TYPE_GLYCINE    = 2,
    MD_RAMACHANDRAN_TYPE_PROLINE    = 3,
    MD_RAMACHANDRAN_TYPE_PREPROL    = 4,
};

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_range_t {
    int32_t beg;
    int32_t end;
} md_range_t;

// Single bond between two entities represented by the indices
typedef struct md_bond_t {
    md_atom_idx_t idx[2];
} md_bond_t;

typedef struct md_backbone_atoms_t {
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
    md_atom_idx_t o;
} md_backbone_atoms_t;

typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

typedef struct md_molecule_t {
    md_molecule_o* inst;

    struct {
        int64_t             count;
        float*              x;
        float*              y;
        float*              z;
        float*              radius;
        float*              mass;
        md_element_t*       element;
        const char**        name;
        md_flag_t*          flags;                          // Auxillary bit buffer for flagging individual atoms
        md_residue_idx_t*   residue_idx;
        md_chain_idx_t*     chain_idx;
    } atom;

    struct {
        int64_t             count;
        const char**        name;
        md_residue_id_t*    id;
        md_range_t*         atom_range;
        md_range_t*         internal_covalent_bond_range;   // Range of covalent bonds within the resuidue
        md_range_t*         complete_covalent_bond_range;   // Range of covalent bonds that in anyway is part of the residue
    } residue;

    struct {
        int64_t       count;
        const char**  id;
        md_range_t*   residue_range;
        md_range_t*   atom_range;
        md_range_t*   backbone_range;
    } chain;

    struct {
        int64_t                   count;
        md_backbone_atoms_t*      atoms;
        md_backbone_angles_t*     angle;
        md_secondary_structure_t* secondary_structure;
        md_ramachandran_type_t*   ramachandran_type;
        md_residue_idx_t*         residue_idx;            // Index to the residue which contains the backbone
    } backbone;

    struct {
        int64_t  count;
        md_bond_t* bond;
    } covalent_bond;
} md_molecule_t;

#if 0
// This is not really well defined yet since it has no real use case so far.
struct md_macro_molecule {
    int64_t num_molecules;
    md_molecule* molecules;

    struct {
        int64_t             count;
        md_molecule_idx_t*    idx;
        const char**        label;
        float               (*transform)[4][4];
    } instance;
};
#endif

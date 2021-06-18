#ifndef _MD_MOLECULE_H
#define _MD_MOLECULE_H

#include <stdint.h>
#include <stdbool.h>

struct md_allocator;

struct md_molecule_o; // Opaque data blob

typedef int32_t                     md_atom_idx;
typedef int32_t                     md_residue_idx;
typedef int32_t                     md_backbone_idx;
typedef int32_t                     md_residue_id;
typedef int32_t                     md_chain_idx;
typedef int32_t                     md_molecule_idx;
typedef uint32_t                    md_secondary_structure;
typedef uint8_t                     md_flag;
typedef uint8_t                     md_element;
typedef struct md_range             md_range;
typedef struct md_bond              md_bond;
typedef struct md_bond_entry        md_bond_entry;
typedef struct md_backbone_atoms    md_backbone_atoms;
typedef struct md_backbone_angles   md_backbone_angles;
typedef struct md_molecule          md_molecule;
typedef struct md_macro_molecule    md_macro_molecule;
typedef struct md_molecule_o        md_molecule_o;

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
struct md_range {
    int32_t beg;
    int32_t end;
};

// Single bond between two entities represented by the indices
struct md_bond {
    md_atom_idx idx[2];
};

struct md_backbone_atoms {
    md_atom_idx n;
    md_atom_idx ca;
    md_atom_idx c;
    md_atom_idx o;
};

struct md_backbone_angles {
    float phi;
    float psi;
};

typedef struct md_molecule {
    md_molecule_o* inst;

    struct {
        int64_t         count;
        float*          x;
        float*          y;
        float*          z;
        float*          radius;
        float*          mass;
        md_element*     element;
        const char**    name;
        md_flag*        flags;       // Auxillary bit buffer for flagging individual atoms
        md_residue_idx* residue_idx;
        md_chain_idx*   chain_idx;
    } atom;

    struct {
        int64_t                 count;
        const char**            name;
        md_residue_id*          id;
        md_range*               atom_range;
        md_range*               internal_covalent_bond_range;   // Range of covalent bonds within the resuidue
        md_range*               complete_covalent_bond_range;   // Range of covalent bonds that in anyway is part of the residue
    } residue;

    struct {
        int64_t         count;
        const char**    id;
        md_range*       residue_range;
        md_range*       atom_range;
        md_range*       backbone_range;
    } chain;

    struct {
        int64_t                 count;
        md_backbone_atoms*      atoms;
        md_backbone_angles*     angle;
        md_secondary_structure* secondary_structure;
        md_residue_idx*         residue_idx;            // Index to the residue which contains the backbone
    } backbone;

    struct {
        int64_t  count;
        md_bond* bond;
    } covalent_bond;
} md_molecule_t;

#if 0
// This is not really well defined yet since it has no real use case.
struct md_macro_molecule {
    int64_t num_molecules;
    md_molecule* molecules;

    struct {
        int64_t             count;
        md_molecule_idx*    idx;
        const char**        label;
        float               (*transform)[4][4];
    } instance;
};
#endif

#endif

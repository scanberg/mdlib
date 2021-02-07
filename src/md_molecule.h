#ifndef _MD_MOLECULE_H
#define _MD_MOLECULE_H

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;

typedef uint32_t                    md_size;
typedef uint32_t                    md_atom_idx;
typedef uint32_t                    md_residue_idx;
typedef uint32_t                    md_residue_id;
typedef uint32_t                    md_chain_idx;
typedef uint32_t                    md_secondary_structure;
typedef uint8_t                     md_flag;
typedef uint8_t                     md_element;
typedef struct md_label           md_label;
typedef struct md_range           md_range;
typedef struct md_bond            md_bond;
typedef struct md_bond_entry      md_bond_entry;
typedef struct md_backbone_atoms  md_backbone_atoms;
typedef struct md_molecule        md_molecule;
typedef struct md_model           md_model;

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

enum {
    MD_ATOM_FIELD_X           = (1 << 0),
    MD_ATOM_FIELD_Y           = (1 << 1),
    MD_ATOM_FIELD_Z           = (1 << 2),
    MD_ATOM_FIELD_RADIUS      = (1 << 3),
    MD_ATOM_FIELD_MASS        = (1 << 4),
    MD_ATOM_FIELD_ELEMENT     = (1 << 5),
    MD_ATOM_FIELD_NAME        = (1 << 6),
    MD_ATOM_FIELD_BFACTOR     = (1 << 7),
    MD_ATOM_FIELD_OCCUPANCY   = (1 << 8),
    MD_ATOM_FIELD_FLAGS       = (1 << 9),
    MD_ATOM_FIELD_RESIDUE_IDX = (1 << 10),
    MD_ATOM_FIELD_CHAIN_IDX   = (1 << 11),

    MD_ATOM_FIELD_POS = MD_ATOM_FIELD_X | MD_ATOM_FIELD_Y | MD_ATOM_FIELD_Z,
    MD_ATOM_FIELD_ALL = ~0,
};

enum {
    MD_RESIDUE_FIELD_NAME                 = (1 << 0),
    MD_RESIDUE_FIELD_ID                   = (1 << 1),
    MD_RESIDUE_FIELD_ATOM_RANGE           = (1 << 2),
    MD_RESIDUE_FIELD_BACKBONE_ATOMS       = (1 << 3),
    MD_RESIDUE_FIELD_SECONDARY_STRUCTURE  = (1 << 4),

    MD_RESIDUE_FIELD_ALL = ~0,
};

enum {
    MD_CHAIN_FIELD_ID            = (1 << 0),
    MD_CHAIN_FIELD_RESIDUE_RANGE = (1 << 1),
    MD_CHAIN_FIELD_ATOM_RANGE    = (1 << 2),

    MD_CHAIN_FIELD_ALL = ~0,
};

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
struct md_range {
    uint32_t beg;
    uint32_t end;
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

struct md_molecule {
    struct {
        md_size         count;
        float*            x;
        float*            y;
        float*            z;
        float*            radius;
        float*            mass;
        md_element*     element;
        const char**      name;
        float*            bfactor;
        float*            occupancy;
        md_flag*        flags;       // Auxillary bit buffer for flagging individual atoms, only 8-bit
        md_residue_idx* residue_idx;
        md_chain_idx*   chain_idx;
    } atom;

    struct {
        md_size                 count;
        const char**              name;
        md_residue_id*          id;
        md_range*               atom_range;
        md_backbone_atoms*      backbone_atoms;         // Only amino acids and nucleic acids
        md_secondary_structure* secondary_structure;    // Only amino acids
    } residue;

    struct {
        md_size    count;
        const char** id;
        md_range*  residue_range;
        md_range*  atom_range;
    } chain;

    struct {
        md_size  count;
        md_bond* bond;
    } bond;
};

struct md_model {
    md_molecule molecule;
    struct {
        uint32_t count;
        float (*transform)[4][4];
    } instance;
};

#endif

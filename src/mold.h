#ifndef _MOLD_H_
#define _MOLD_H_

#include <stdint.h>
#include <stdbool.h>

typedef void (mold_error_callback)(const char* str, int len);	// Function callback for receiving errors from mold

typedef int 						mold_error;
typedef uint32_t                    mold_size;
typedef uint32_t                    mold_atom_idx;
typedef uint32_t                    mold_residue_idx;
typedef uint32_t                    mold_residue_id;
typedef uint32_t                    mold_chain_idx;
typedef uint32_t                    mold_secondary_structure;
typedef uint8_t                     mold_flag;
typedef uint8_t                     mold_element;
typedef struct mold_range           mold_range;
typedef struct mold_bond            mold_bond;
typedef struct mold_bond_entry      mold_bond_entry;
typedef struct mold_backbone_atoms  mold_backbone_atoms;
typedef struct mold_molecule        mold_molecule;

// We are sneaky, we encode the secondary structure as a uint8x4 unorm where the the components encode the fraction of each secondary structure type
enum {
    MOLD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MOLD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MOLD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MOLD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

// Range of indices where end is exclusive (e.g. range(0,4) -> [0,1,2,3])
struct mold_range {
    uint32_t beg;
    uint32_t end;
};

// Single bond between two entities represented by the indices
struct mold_bond {
    mold_atom_idx idx[2];
};

// Used as a look up to find which bonds 
struct mold_bond_entry {
    uint32_t offset : 24;
    uint32_t length : 8;
};

struct mold_backbone_atoms {
    mold_atom_idx n;
    mold_atom_idx ca;
    mold_atom_idx c;
    mold_atom_idx o;
};

struct mold_molecule {
    struct {
        mold_size           count;
        float*              x;
        float*              y;
        float*              z;
        float*              radius;
        float*              mass;
        mold_element*       element;
        const char**        name;
        float*              bfactor;
        float*              occupancy;
        mold_residue_idx*   residue_idx;
        mold_chain_idx*     chain_idx;
        mold_bond_entry*    bond;
        mold_flag*          flags;       // Auxillary bit buffer for flagging individual atoms
    } atom;

    struct {
        mold_size                 count;
        const char**              name;
        mold_residue_id*          id;
        mold_range*               atom_range;
        mold_backbone_atoms*      backbone_atoms;         // Only amino acids and nucleic acids
        mold_secondary_structure* secondary_structure;    // Only amino acids
    } residue;

    struct {
        mold_size    count;
        const char** id;
        mold_range*  residue_range; 
        mold_range*  atom_range;
    } chain;

    struct {
        mold_size  count;
        mold_bond* bond;
    } bond;
};

#endif
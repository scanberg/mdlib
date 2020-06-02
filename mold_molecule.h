#ifndef _MOLD_MOLECULE_H_
#define _MOLD_MOLECULE_H_

#include <stdint.h>

typedef struct mold_molecule        mold_molecule;
typedef struct mold_range           mold_range;
typedef struct mold_bond            mold_bond;
typedef struct mold_bond_entry      mold_bond_entry;
typedef struct mold_backbone_atoms  mold_backbone_atoms;
typedef uint32_t                    mold_secondary_structure; 

// We are snieky, we encode the secondary structure as a uint8x4 probability where the the components encode the fraction of each secondary structure type
enum {
    MOLD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MOLD_SECONDARY_STRUCTURE_COIL    = 0x000000FF,
    MOLD_SECONDARY_STRUCTURE_HELIX   = 0x0000FF00,
    MOLD_SECONDARY_STRUCTURE_SHEET   = 0x00FF0000
};

// Range of atom indices where end is exclusive (e.g. range(0,4) -> [0,1,2,3])
struct mold_range {
    uint32_t beg;
    uint32_t end;
};

// Single bond between two entities represented by the indices
struct mold_bond {
    uint32_t idx[2];
};

// Used as a look up to find which bonds 
struct mold_bond_entry {
    uint32_t offset : 24;
    uint32_t length : 8;
};

struct mold_backbone_atoms {
    uint32_t n;
    uint32_t ca;
    uint32_t c;
    uint32_t o;
};


struct mold_molecule {
    struct {
        uint32_t count;
        float*           x;
        float*           y;
        float*           z;
        float*           radius;
        float*           mass;
        uint8_t*         element;
        char**           name;
        float*           bfactor;
        float*           occupancy;
        uint32_t*        residue_idx;
        uint32_t*        chain_idx;
        mold_bond_entry* bond;
        uint8_t*         flags;       // Auxillary bit buffer for flagging individual atoms
    } atom;

    struct {
        uint32_t count;
        char**                    name;
        uint32_t*                 id;
        mold_range*               atom_range;
        mold_backbone_atoms*      backbone_atoms;         // Only amino acids and nucleic acids
        mold_secondary_structure* secondary_structure;    // Only amino acids
    } residue;

    struct {
        uint32_t count;
        char**      id;
        mold_range* residue_range;
        mold_range* atom_range;
    } chain;

    struct {
        uint32_t count;
        mold_bond* bond;
    } bond;
};

#endif
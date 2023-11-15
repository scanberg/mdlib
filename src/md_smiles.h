#pragma once

#include <stdint.h>

// This is a non feature complete parser for SMILES strings
// It does not resolve the implicit hydrogens, nor does it resolve the aromaticity

enum {
    MD_SMILES_FLAG_CHIRAL_1 = 1,    // @
    MD_SMILES_FLAG_CHIRAL_2 = 3,    // @@
    MD_SMILES_FLAG_AROMATIC = 4,
};

typedef enum {
    MD_SMILES_NODE_UNKNOWN,
    MD_SMILES_NODE_ATOM,
    MD_SMILES_NODE_BOND,
    MD_SMILES_NODE_BRANCH_OPEN,
    MD_SMILES_NODE_BRANCH_CLOSE,
    MD_SMILES_NODE_RING_CLOSURE,
} md_smiles_type_t;

typedef struct md_smiles_node_t {
    md_smiles_type_t type;
    union {
        struct {
            char symbol[3];
            uint8_t flags;
            uint8_t isotope;
            uint8_t h_count;    // Hydrogen count
            int8_t  charge;     
            uint8_t atom_class; // Atom Class
        } atom;
        struct {
            uint32_t index;    // Index for ring closure
        } ring;
        struct {
            char symbol;
            uint8_t order;
        } bond;
    };
} md_smiles_node_t;

#ifdef __cplusplus
extern "C" {
#endif

// Returns the number of nodes parsed
// The nodes array must be large enough to hold the parsed nodes
// The upper bound of dst nodes can be estimated as the length of the string
int64_t md_smiles_parse(md_smiles_node_t* out_nodes, int64_t in_node_cap, const char* in_str_ptr, int64_t in_str_len);

#ifdef __cplusplus
}
#endif
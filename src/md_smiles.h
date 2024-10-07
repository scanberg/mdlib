#pragma once

#include <stdint.h>
#include <stddef.h>

// This is a non feature complete parser for SMILES strings
// It does not resolve the implicit hydrogens, nor does it resolve the aromaticity

enum {
    MD_SMILES_FLAG_NONE     = 0,
    MD_SMILES_FLAG_CHIRAL_1 = 1,    // @
    MD_SMILES_FLAG_CHIRAL_2 = 2,    // @@
    MD_SMILES_FLAG_AROMATIC = 4,
    MD_SMILES_FLAG_UP       = 8,    /* / */
    MD_SMILES_FLAG_DOWN     = 16,   /* \ */
};

enum {
    MD_SMILES_NODE_UNKNOWN,
    MD_SMILES_NODE_ATOM,
    MD_SMILES_NODE_BOND,
    MD_SMILES_NODE_BRANCH_OPEN,
    MD_SMILES_NODE_BRANCH_CLOSE,
    MD_SMILES_NODE_BRIDGE,
};

typedef int md_smiles_type_t;

typedef struct md_smiles_node_t {
    md_smiles_type_t type;
    union {
        struct {
            uint8_t  element;		    // Atom element, '*' is mapped to 0 and signifies any/unknown atom
            uint8_t  isotope;
            uint8_t  hydrogen_count;    // Explicit hydrogen count
            int8_t   charge;     
            uint16_t atom_class;
            uint8_t  flags;
        } atom;

        struct {
            uint32_t index;
        } bridge;

        struct {
            char     symbol;
            uint8_t  order;
            uint8_t  flags;
        } bond;
    };
} md_smiles_node_t;

#ifdef __cplusplus
extern "C" {
#endif

// Returns the number of nodes parsed
// The nodes array must be large enough to hold the parsed nodes
// The upper bound of dst nodes can be estimated as the length of the input string
size_t md_smiles_parse(md_smiles_node_t* out_node_arr, size_t node_cap, const char* str_ptr, size_t str_len);

#ifdef __cplusplus
}
#endif

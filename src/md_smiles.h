#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_types.h>

struct md_allocator_i;

// ### SMILES ###
// Parses a SMILES string (OpenSMILES) into a molecular graph: the atoms and bonds as written, nothing perceived.
// It is the reader behind md_match_query_init_smiles, and has no knowledge of matching.
//
// SUPPORTED
//   - Organic subset atoms:  B C N O P S F Cl Br I, aromatic b c n o p s, and '*'
//   - Bracket atoms:         [isotope symbol chirality hcount charge :class], every element, aromatic se as te too.
//                            Charge as +, -, ++, --, +n, -n. Chirality as @ and @@; the extended classes (@TH1, @SP1,
//                            @TB1, @OH1, ...) are accepted and kept as written but not interpreted.
//   - Bonds:                 - = # $ : / \, and none (implicit)
//   - Branches:              ( ), nested to any depth
//   - Ring closures:         0-9 and %nn, with a bond symbol on either or both sides (C=1CCCCC1, C1CCCCC=1). Both
//                            sides giving different symbols is an error.
//   - Disconnection:         '.', the parts are numbered in 'num_components'
//
// NOT DONE
//   - Implicit hydrogens of organic subset atoms are not computed: h_count is meaningful for bracket atoms only.
//     Matching reads an organic subset atom as having an open hydrogen count, so nothing here needs it. When something
//     does (building a molecule from SMILES), it belongs in a separate pass over this graph, with the valence model
//     that pass needs.
//   - No aromaticity perception, no kekulization. Lowercase atoms are flagged aromatic, ':' bonds are aromatic, and an
//     implicit bond between two aromatic atoms is reported as IMPLICIT: it is the consumer's call whether that means
//     aromatic (OpenSMILES says it does).
//   - Chirality is recorded, not resolved. The neighbour order it refers to is the order in which the atom's bonds
//     appear in 'bonds', which follows the string (ring closure bonds at the position of their digit).
//
// ERRORS
//   The parse stops at the first syntax error and reports where it is (offset into the string) and what was expected.
//   A failed parse leaves no partial graph behind: out is zeroed. Errors include unbalanced parentheses, a ring
//   closure left open, a ring closure bonding an atom to itself or duplicating an existing bond, an unknown element,
//   a bond symbol with no atom to follow, and trailing characters.
//   Leading and trailing whitespace is ignored. Anything after the first whitespace inside the string (the name
//   field of a SMILES file line) is an error, not ignored: a pattern given in a script should not silently lose half
//   of itself.

typedef enum md_smiles_atom_flags_t {
    MD_SMILES_ATOM_AROMATIC   = 0x01,   // Written lowercase
    MD_SMILES_ATOM_BRACKET    = 0x02,   // Written in brackets: h_count and charge are as given, absent meaning 0
    MD_SMILES_ATOM_CHIRAL_CCW = 0x04,   // @
    MD_SMILES_ATOM_CHIRAL_CW  = 0x08,   // @@
    MD_SMILES_ATOM_CHIRAL_EXT = 0x10,   // An extended chirality class (@TH, @AL, @SP, @TB, @OH), number in 'chiral_class'
} md_smiles_atom_flags_t;

typedef struct md_smiles_atom_t {
    md_atomic_number_t z;               // 0 for '*'
    uint8_t  flags;                     // md_smiles_atom_flags_t
    uint8_t  h_count;                   // Bracket atoms only
    int8_t   charge;                    // Bracket atoms only
    uint16_t isotope;                   // 0 when not given
    uint16_t atom_class;                // 0 when not given
    uint8_t  chiral_class;              // MD_SMILES_ATOM_CHIRAL_EXT only
    uint32_t offset;                    // Character offset of the atom in the string
} md_smiles_atom_t;

typedef enum md_smiles_bond_order_t {
    MD_SMILES_BOND_IMPLICIT  = 0,       // No symbol: single, or aromatic between two aromatic atoms
    MD_SMILES_BOND_SINGLE    = 1,       // '-', '/', '\'
    MD_SMILES_BOND_DOUBLE    = 2,       // '='
    MD_SMILES_BOND_TRIPLE    = 3,       // '#'
    MD_SMILES_BOND_QUADRUPLE = 4,       // '$'
    MD_SMILES_BOND_AROMATIC  = 5,       // ':'
} md_smiles_bond_order_t;

typedef enum md_smiles_bond_flags_t {
    MD_SMILES_BOND_UP    = 0x01,        // '/'
    MD_SMILES_BOND_DOWN  = 0x02,        // '\'
    MD_SMILES_BOND_RING  = 0x04,        // Formed by a ring closure
} md_smiles_bond_flags_t;

typedef struct md_smiles_bond_t {
    uint32_t a;                         // The atom written first
    uint32_t b;
    uint8_t  order;                     // md_smiles_bond_order_t
    uint8_t  flags;                     // md_smiles_bond_flags_t
} md_smiles_bond_t;

typedef struct md_smiles_t {
    size_t num_atoms;
    md_smiles_atom_t* atoms;            // In the order written

    size_t num_bonds;
    md_smiles_bond_t* bonds;            // In the order written

    size_t num_components;              // Parts separated by '.'

    struct md_allocator_i* alloc;
} md_smiles_t;

typedef struct md_smiles_error_t {
    size_t offset;                      // Character offset into the string given
    char   message[96];
} md_smiles_error_t;

#ifdef __cplusplus
extern "C" {
#endif

// Returns false on a syntax error, described in err (optional). The graph is allocated from alloc.
bool md_smiles_parse(md_smiles_t* out, str_t str, struct md_allocator_i* alloc, md_smiles_error_t* err);
void md_smiles_free(md_smiles_t* smiles);

#ifdef __cplusplus
}
#endif

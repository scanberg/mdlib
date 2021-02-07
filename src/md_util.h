#ifndef _MD_UTIL_H_
#define _MD_UTIL_H_

#include "md_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_util_backbone_args md_util_backbone_args;
typedef struct md_util_secondary_structure_args md_util_secondary_structure_args;

struct md_util_backbone_args {
    struct {
        uint32_t count;
        const char** name;
    } atom;

    struct {
        uint32_t count;
        const md_range* atom_range;
    } residue;
};

struct md_util_secondary_structure_args {
    struct {
        const float* x;
        const float* y;
        const float* z;
    } atom;

    struct {
        uint32_t count;
        const md_backbone_atoms* backbone_atoms;
    } residue;

    struct {
        uint32_t count;
        const md_range* residue_range;
    } chain;
};

inline bool md_util_protein_backbone_atoms_valid(md_backbone_atoms prot) {
    return (prot.ca != prot.c) && (prot.ca != prot.o) && (prot.c != prot.o);
}

bool md_util_extract_backbone_atoms(md_backbone_atoms* backbone_atoms, const char** residue_atom_name, md_range residue_atom_range);
void md_util_extract_backbone(md_backbone_atoms* backbone_atoms, const md_util_backbone_args* args);

// Computes secondary structures for residues with the supplied backbone_atoms.
// Does not allocate any data, it assumes that secondary_structures has the same length as args->residue.count
void md_util_extract_secondary_structure(md_secondary_structure* secondary_structures, const md_util_secondary_structure_args* args);

#ifdef __cplusplus
}
#endif

#endif
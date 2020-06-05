#ifndef _MOLD_UTIL_H_
#define _MOLD_UTIL_H_

#include <stdint.h>
#include <stdbool.h>
#include "mold_molecule.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mold_util_backbone_args mold_util_backbone_args;
typedef struct mold_util_secondary_structure_args mold_util_secondary_structure_args;

struct mold_util_backbone_args {
    // OUTPUT
    mold_backbone_atoms* backbone_atoms;

    // INPUT
    struct {
        uint32_t count;
        const char** name;
    } atom;

    struct {
        uint32_t count;
        const mold_range* atom_range;
    } residue;
};

struct mold_util_secondary_structure_args {
    // OUTPUT
    mold_secondary_structure* secondary_structure;

    // INPUT
    struct {
        const float* x;
        const float* y;
        const float* z;
    } atom;

    struct {
        uint32_t count;
        const mold_backbone_atoms* backbone_atoms;
    } residue;

    struct {
        uint32_t count;
        const mold_range* residue_range;
    } chain;
};

inline bool mold_util_protein_backbone_atoms_valid(mold_backbone_atoms prot) {
    return (prot.ca != prot.c  && prot.ca != prot.o) && (prot.c != prot.o);
}

bool mold_util_extract_backbone_atoms(mold_backbone_atoms* out_backbone_atoms, const char** residue_atom_name, mold_range residue_atom_range);
void mold_util_extract_backbone(const mold_util_backbone_args* args);
void mold_util_extract_secondary_structure(const mold_util_secondary_structure_args* args);

#ifdef __cplusplus
}
#endif

#endif
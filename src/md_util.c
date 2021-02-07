#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable : 4204)            // non-constant aggregate initialization (fully supported in proper c compilers)
#endif

#include "md_util.h"
#include "core/common.h"

#include <math.h>
#include <string.h>

#define EXTRACT_POS(i) {args->atom.x[i], args->atom.y[i], args->atom.z[i]}

#ifdef __cplusplus
extern "C" {
#endif

internal inline float distance(const float v[3], const float u[3]) {
    const float d[3] = {v[0] - u[0], v[1] - u[1], v[2] - u[2]};
    return sqrtf(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

internal inline bool cmp1(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == '\0';
}

internal inline bool cmp2(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == '\0';
}

bool md_util_extract_backbone_atoms(md_backbone_atoms* out_backbone_atoms, const char** residue_atom_name, md_range residue_range) {
    uint32_t bits = 0;
    md_backbone_atoms bb = {0};
    for (uint32_t i = residue_range.beg; i < residue_range.end; ++i) {
        if (!(bits & 1) && cmp1(residue_atom_name[i], "N"))  { bb.n  = i; bits |= 1;  continue; }
        if (!(bits & 2) && cmp2(residue_atom_name[i], "CA")) { bb.ca = i; bits |= 2;  continue; }
        if (!(bits & 4) && cmp1(residue_atom_name[i], "C"))  { bb.c  = i; bits |= 4;  continue; }
        if (!(bits & 8) && cmp1(residue_atom_name[i], "O"))  { bb.o  = i; bits |= 8;  continue; }
        if (!(bits & 8) && i == (residue_range.end - 1) && residue_atom_name[i][0] == 'O') {
            bb.o = i; bits |= 8; continue;
        }
    }

    // If we have CA, C and O, we have enough for computing the backbone
    if (bits & (2|4|8)) {
        if (out_backbone_atoms) *out_backbone_atoms = bb;
        return true;
    }
    return false;
}

void md_util_extract_backbone(md_backbone_atoms* backbone_atoms, const md_util_backbone_args* args) {
    ASSERT(backbone_atoms);
    memset(backbone_atoms, 0, args->residue.count * sizeof(md_backbone_atoms));
    for (uint32_t i = 0; i < args->residue.count; ++i) {
        md_util_extract_backbone_atoms(&backbone_atoms[i], args->atom.name, args->residue.atom_range[i]);
    }
}

static inline bool zhang_skolnick_ss(const md_util_secondary_structure_args* args, md_range res_range, int i, const float distances[3], float delta) {
    for (int j = MAX((int)res_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= (int)res_range.end) continue;
            const int ca_j = args->residue.backbone_atoms[j].ca;
            const int ca_k = args->residue.backbone_atoms[j + k].ca;
            const float pos_j[3] = {args->atom.x[ca_j], args->atom.y[ca_j], args->atom.z[ca_j]};
            const float pos_k[3] = {args->atom.x[ca_k], args->atom.y[ca_k], args->atom.z[ca_k]};
            const float d = distance(pos_j, pos_k);
            if (ABS(d - distances[k - 2]) > delta) {
                return false;
            }
        }
    }
    return true;
}

static inline bool is_helical(const md_util_secondary_structure_args* args, md_range res_range, int i) {
    const float distances[] = { 5.45, 5.18, 6.37 };
    const float delta = 2.1;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

static inline bool is_sheet(const md_util_secondary_structure_args* args, md_range res_range, int i) {
    const float distances[] = { 6.1, 10.4, 13.0 };
    const float delta = 1.42;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

// TM-align: a protein structure alignment algorithm based on the Tm-score
// doi:10.1093/nar/gki524
void md_util_extract_secondary_structure(md_secondary_structure* secondary_structure, const md_util_secondary_structure_args* args) {
    ASSERT(args);
    ASSERT(args->atom.x);
    ASSERT(args->atom.y);
    ASSERT(args->atom.z);
    ASSERT(args->residue.backbone_atoms);
    ASSERT(args->chain.residue_range);

    for (uint32_t chain_idx = 0; chain_idx < args->chain.count; ++chain_idx) {
        const md_range range = args->chain.residue_range[chain_idx];
        if (range.end - range.beg < 4) continue;
        
        // Classify residues
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_secondary_structure ss = MD_SECONDARY_STRUCTURE_COIL;
            if (is_sheet(args, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_SHEET;
            }
            else if (is_helical(args, range, i)) {
                ss = MD_SECONDARY_STRUCTURE_HELIX;
            }

            secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        md_secondary_structure* ss = secondary_structure;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] != MD_SECONDARY_STRUCTURE_COIL && ss[i] != ss[i-1] && ss[i-1] == ss[i+1]) ss[i] = ss[i-1];
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MD_SECONDARY_STRUCTURE_COIL;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i] != ss[i-1] && ss[i] != ss[i+1]) ss[i] = MD_SECONDARY_STRUCTURE_COIL;
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MD_SECONDARY_STRUCTURE_COIL;
    }
}

#ifdef __cplusplus
}
#endif
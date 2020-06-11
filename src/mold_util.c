#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#pragma warning( disable : 4204)            // non-constant aggregate initialization (fully supported in proper c compilers)
#endif

#include "mold_util.h"
#include <math.h>
#include <string.h>

#ifndef ASSERT
#include <assert.h>
#define ASSERT assert
#endif
#define DEG_TO_RAD (3.14159265f / 180.0f)
#define EXTRACT_POS(i) {args->atom.x[i], args->atom.y[i], args->atom.z[i]}
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define ABS(x) ((x) > 0 ? (x) : -(x))

#ifdef __cplusplus
extern "C" {
#endif

typedef struct vec3 vec3;
struct vec3 {
    float x, y, z;
};

static inline float dot(vec3 v, vec3 u) {
    return (v.x * u.x + v.y * u.y + v.z * u.z);
}

static inline float length(vec3 v) {
    return sqrtf(dot(v,v));
}

static inline vec3 normalize(vec3 v) {
    const float l = length(v);
    vec3 u = {v.x/l, v.y/l, v.z/l};
    return u;
}

static inline float distance(vec3 v, vec3 u) {
    const vec3 d = {v.x - u.x, v.y - u.y, v.z - u.z};
    return length(d);
}

static inline bool cmp1(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == '\0';
}

static inline bool cmp2(const char* str, const char* ref) {
    return str[0] == ref[0] && str[1] == ref[1] && str[2] == '\0';
}

bool mold_util_extract_backbone_atoms(mold_backbone_atoms* out_backbone_atoms, const char** residue_atom_name, mold_range residue_range) {
    uint32_t bits = 0;
    mold_backbone_atoms bb = {0};
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

void mold_util_extract_backbone(const mold_util_backbone_args* args) {
    ASSERT(args->backbone_atoms);
    memset(args->backbone_atoms, 0, args->residue.count * sizeof(mold_backbone_atoms));
    for (uint32_t i = 0; i < args->residue.count; ++i) {
        mold_util_extract_backbone_atoms(&args->backbone_atoms[i], args->atom.name, args->residue.atom_range[i]);
    }
}

static inline bool zhang_skolnick_ss(const mold_util_secondary_structure_args* args, mold_range res_range, int i, float distances[3], float delta) {
    const int res_count = res_range.end - res_range.beg;
    for (int j = MAX((int)res_range.beg, i - 2); j <= i; ++j) {
        for (int k = 2; k < 5; ++k) {
            if (j + k >= res_range.end) continue;
            const vec3 ca_j = EXTRACT_POS(args->residue.backbone_atoms[j].ca);
            const vec3 ca_k = EXTRACT_POS(args->residue.backbone_atoms[j + k].ca);
            const float d = distance(ca_j, ca_k);
            if (ABS(d - distances[k - 2]) > delta) {
                return false;
            }
        }
    }
    return true;
}

static inline bool is_helical(const mold_util_secondary_structure_args* args, mold_range res_range, int i) {
    const float distances[] = { 5.45, 5.18, 6.37 };
    const float delta = 2.1;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

static inline bool is_sheet(const mold_util_secondary_structure_args* args, mold_range res_range, int i) {
    const float distances[] = { 6.1, 10.4, 13.0 };
    const float delta = 1.42;
    return zhang_skolnick_ss(args, res_range, i, distances, delta);
}

// TM-align: a protein structure alignment algorithm based on the Tm-score
// doi:10.1093/nar/gki524
void mold_util_extract_secondary_structure(const mold_util_secondary_structure_args* args) {
    ASSERT(args);
    ASSERT(args->secondary_structure);
    ASSERT(args->atom.x);
    ASSERT(args->atom.y);
    ASSERT(args->atom.z);
    ASSERT(args->residue.backbone_atoms);
    ASSERT(args->chain.residue_range);

    for (uint32_t chain_idx = 0; chain_idx < args->chain.count; ++chain_idx) {
        const mold_range range = args->chain.residue_range[chain_idx];
        if (range.end - range.beg < 4) continue;
        
        // Classify residues
        for (uint32_t i = range.beg; i < range.end; ++i) {
            mold_secondary_structure ss = MOLD_SECONDARY_STRUCTURE_COIL;
            if (is_sheet(args, range, i)) {
                ss = MOLD_SECONDARY_STRUCTURE_SHEET;
            }
            else if (is_helical(args, range, i)) {
                ss = MOLD_SECONDARY_STRUCTURE_HELIX;
            }

            args->secondary_structure[i] = ss;
        }

        // Set squished isolated structures to the surrounding (only for sheets and helices)
        mold_secondary_structure* ss = args->secondary_structure;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i-1] != MOLD_SECONDARY_STRUCTURE_COIL && ss[i] != ss[i-1] && ss[i-1] == ss[i+1]) ss[i] = ss[i-1];
        }

        // Set remaining isolated structures to coil
        if (ss[range.beg] != ss[range.beg + 1]) ss[range.beg] = MOLD_SECONDARY_STRUCTURE_COIL;
        for (uint32_t i = range.beg + 1; i < range.end - 1; ++i) {
            if (ss[i] != ss[i-1] && ss[i] != ss[i+1]) ss[i] = MOLD_SECONDARY_STRUCTURE_COIL;
        }
        if (ss[range.end - 1] != ss[range.end - 2]) ss[range.end - 1] = MOLD_SECONDARY_STRUCTURE_COIL;
    }
}

#ifdef __cplusplus
}
#endif
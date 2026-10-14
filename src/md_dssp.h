#pragma once

#include <stdbool.h>
#include <stddef.h>

#include <md_types.h>

struct vec3_t;
struct md_unitcell_t;
struct md_protein_backbone_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### SECONDARY STRUCTURE ###
//
// Secondary structure (DSSP like) from the coordinates of a frame. Atomistic ranges get DSSP (Kabsch & Sander) from the
// hydrogen bonds of their backbone. Traces of coarse grained beads (md_protein_backbone_range_is_trace) get the same
// patterns from bonds between their peptide groups, which the beads give (CG-DSSP, see md_dssp.c). It is fit for Martini 3
// BB beads (the center of N, CA, C and O), and to tell helix, sheet and coil apart rather than every class of DSSP.
// Bonds and bends are taken between nearest periodic images (cell), so a protein split by the boundary gets what it has whole.
// Computed into the caller's array of one entry per segment of the backbone (capacity at least backbone->segment.count).
bool md_dssp_compute(md_secondary_structure_t secondary_structures[], size_t capacity, const struct vec3_t* xyz, const struct md_unitcell_t* cell, const struct md_protein_backbone_data_t* backbone);

#ifdef __cplusplus
}
#endif

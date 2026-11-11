#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <md_types.h>
#include <core/md_common.h>

struct md_system_t;
struct md_allocator_i;

// The protein and nucleic acid backbones of a system.

// Secondary structure types (DSSP-like)
typedef enum {
    MD_SECONDARY_STRUCTURE_UNKNOWN = 0,
    MD_SECONDARY_STRUCTURE_COIL,
    MD_SECONDARY_STRUCTURE_TURN,
    MD_SECONDARY_STRUCTURE_BEND,
    MD_SECONDARY_STRUCTURE_HELIX_310,
    MD_SECONDARY_STRUCTURE_HELIX_ALPHA,
    MD_SECONDARY_STRUCTURE_HELIX_PI,
    MD_SECONDARY_STRUCTURE_BETA_SHEET,
    MD_SECONDARY_STRUCTURE_BETA_BRIDGE,
} md_secondary_structure_t;

typedef uint8_t md_ramachandran_type_t;

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL,
    MD_RAMACHANDRAN_TYPE_GLYCINE,
    MD_RAMACHANDRAN_TYPE_PROLINE,
    MD_RAMACHANDRAN_TYPE_PREPROL,
};

// The atoms of an amino acid of a protein backbone, -1 where it has none. A coarse grained bead which stands for the
// residue's backbone (the BB bead of Martini) is held in ca, with every other index -1 (md_protein_backbone_data_t).
typedef struct md_amino_acid_atoms_t {
    // Backbone
    md_atom_idx_t n;
    md_atom_idx_t ca;
    md_atom_idx_t c;
	// Side chain and appendages
    md_atom_idx_t o;
    md_atom_idx_t cb;
    md_atom_idx_t hn;
} md_amino_acid_atoms_t;

// Backbone angles
// φ (phi) is the angle in the chain C' − N − Cα − C'
// ψ (psi) is the angle in the chain N − Cα − C' − N
// https://en.wikipedia.org/wiki/Dihedral_angle#/media/File:Protein_backbone_PhiPsiOmega_drawing.svg)
typedef struct md_backbone_angles_t {
    float phi;
    float psi;
} md_backbone_angles_t;

typedef struct md_nucleic_acid_atoms_t {
    // Backbone
    md_atom_idx_t p;
    md_atom_idx_t o5;
    md_atom_idx_t c5;
    md_atom_idx_t c4;
    md_atom_idx_t c3;
    md_atom_idx_t o3;
    // Nucleoside
    md_atom_idx_t c1;
    md_atom_idx_t c2;
    md_atom_idx_t o4;
} md_nucleic_acid_atoms_t;

#ifdef __cplusplus
extern "C" {
#endif

const char* md_secondary_structure_name(md_secondary_structure_t ss);   // "" for a value out of range

#ifdef __cplusplus
}
#endif

// The backbones hold only what is static with regard to the topology: which atoms of which residues form them, in
// consecutive ranges along a chain. What depends on a frame (the backbone angles, the secondary structure) or serves a
// single consumer (the ramachandran classification) is computed from them by whoever needs it, into arrays of its own
// with one entry per segment (md_util_backbone_angles_compute, md_dssp_compute,
// md_util_backbone_ramachandran_classify).
//
// A system without residues (an XYZ file) has its backbones found from the bonds alone (BACKBONES FROM THE BONDS in
// md_util.c): its segments have no residue (res_idx -1) and its ranges no part (part_idx -1).
//
// A protein range is ATOMISTIC, each segment the N, CA, C and O of its residue (and HN, CB where they are), or a TRACE:
// one bead per residue of a coarse grained model (the BB bead of Martini), held in ca, with n, c, o, hn and cb -1. A
// range is one or the other throughout (md_protein_backbone_range_is_trace). The secondary structure of a trace is
// inferred from its beads alone, and a trace has no backbone angles.
typedef struct md_protein_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the segments
        md_part_idx_t* part_idx; // The part (a strand of a polymer) the range is of, -1 for none
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_amino_acid_atoms_t* atoms;
        md_residue_idx_t* res_idx;                     // Index to the residue which contains the backbone, -1 for none
    } segment;
} md_protein_backbone_data_t;

// A nucleic range is ATOMISTIC, each segment the backbone and sugar atoms of its nucleotide, or COARSE GRAINED: the
// backbone beads of a Martini nucleotide in the fields of the atoms they hold, the phosphate bead (BB1) in p, the sugar
// beads with the C4' (BB2) and the C1' (BB3) in c4 and c1, the other fields -1 (p too at a 5' end). A range is one or
// the other throughout (md_nucleic_backbone_range_is_coarse_grained).
typedef struct md_nucleic_backbone_data_t {
    // This holds the consecutive ranges which form the backbones
    struct {
        size_t count;
        uint32_t* offset; // Offsets into the backbone fields stored bellow
        md_part_idx_t* part_idx; // The part (a strand of a polymer) the range is of, -1 for none
    } range;

    // These fields share the same length 'count'
    struct {
        size_t count;
        md_nucleic_acid_atoms_t* atoms;
        md_residue_idx_t* res_idx;                     // Index to the residue which contains the backbone segment, -1 for none
    } segment;
} md_nucleic_backbone_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// Whether a protein backbone range is a trace of coarse grained beads rather than atoms (see md_protein_backbone_data_t).
// A backbone without atoms, only the residues of its segments, has no traces.
static inline bool md_protein_backbone_range_is_trace(const md_protein_backbone_data_t* backbone, size_t range_idx) {
    MD_ASSERT(backbone);
    MD_ASSERT(range_idx < backbone->range.count);
    return backbone->segment.atoms && backbone->segment.atoms[backbone->range.offset[range_idx]].n < 0;
}

// Whether a nucleic backbone range is of coarse grained beads rather than atoms (see md_nucleic_backbone_data_t)
static inline bool md_nucleic_backbone_range_is_coarse_grained(const md_nucleic_backbone_data_t* backbone, size_t range_idx) {
    MD_ASSERT(backbone);
    MD_ASSERT(range_idx < backbone->range.count);
    return backbone->segment.atoms && backbone->segment.atoms[backbone->range.offset[range_idx]].o5 < 0;
}

// ### INFERENCE ###

// The backbones of the system (md_protein_backbone_data_t, md_nucleic_backbone_data_t) and the breaks of its chains:
// ranges of at least 3 residues along each chain of a polymer, ending where the chain breaks (CHAIN BREAKS in
// md_util.h).
//   ATOMS  From the covalent bonds alone, whatever the atoms and residues are called (BACKBONES FROM THE BONDS in
//          md_backbone.c): an atomistic system without bonds has no backbones. A segment is of the residue of its CA (C4')
//          where the system has residues, a range of the part of its residues.
//   BEADS  The traces of a coarse grained model by the backbone beads its residues were given (Martini's BB, BB1 to
//          BB3), along each strand of a polymer (the peptide and nucleic acid parts). Linked by their bonds, else by
//          consecutive sequence numbers.
// The ranges come in the order of their first atom. Requires the bonds and the roles of the atoms
// (md_util_system_infer_residue_flags), and the parts for the beads. Into out (allocated with alloc); the system is read
// and never written. md_util_system_infer (MD_UTIL_INFER_BACKBONE_BIT, after the parts) gives the system the backbones
// and marks the breaks (md_system_set_chain_breaks).
typedef struct md_backbone_data_t {
    md_protein_backbone_data_t protein;
    md_nucleic_backbone_data_t nucleic;
    size_t num_breaks;
    md_residue_idx_t* breaks;       // [num_breaks] The residues after which a chain breaks, ascending
} md_backbone_data_t;

bool md_backbone_infer(md_backbone_data_t* out, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_backbone_data_free(md_backbone_data_t* bb, struct md_allocator_i* alloc);

// Marks the breaks of the chains of the system (CHAIN BREAKS in md_util.h): residue breaks[k] MD_RESIDUE_FLAG_BREAK_AFTER
// and the one after it MD_RESIDUE_FLAG_BREAK_BEFORE, the flags of every other residue cleared. A break is no end of a
// chain: the residues on either side of one lose the terminal flags their classification gave them
// (MD_RESIDUE_FLAG_TERMINAL_BEG, _END, and those of their atoms).
void md_system_set_chain_breaks(struct md_system_t* sys, const md_residue_idx_t breaks[], size_t num_breaks);

#ifdef __cplusplus
}
#endif

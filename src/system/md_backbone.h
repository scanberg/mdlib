#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <md_types.h>
#include <core/md_unitcell.h>

struct md_float3_t;

struct md_system_t;
struct md_system_state_t;
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
// with one entry per segment (md_backbone_angles_compute, md_dssp_compute,
// md_backbone_ramachandran_classify).
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

// ### POLYMER ATOMS, CHAIN BREAKS, THE BACKBONE OF A STATE ###

// Gives the atoms of the amino acids and nucleotides of an atomistic system their roles (MD_ATOM_FLAG_BACKBONE,
// _SIDE_CHAIN, _NUCLEOSIDE, _NUCLEOBASE, _TERMINAL_BEG, _END), found from the covalent bonds alone, whatever the atoms
// are called (POLYMER ATOMS in md_backbone.c), also in a system without residues. Then classifies the residues
// (md_residue_kind_t: amino acid, nucleotide, water, ion, lipid). A residue which holds the backbone of one amino
// acid or nucleotide of a chain (a standard residue by name, or linked to another residue by its backbone) is one,
// MD_RESIDUE_FLAG_RESOLVED, with its terminal groups marked: the link keeps out the ligands which hold the atoms of
// one (GTP, ATP, NAD) while modified residues in a chain (MSE, PSU) come in. A residue which is only named like an
// amino acid or nucleotide gets the kind alone. A residue of coarse grained beads which is a molecule of the Martini
// force fields by its name and all its bead names (a lipid, water, an ion, a sugar, ...) is of the kind of that
// molecule; residues a predefined coarse grained type already classified as water, ion or lipid are left alone.
// Requires the bonds. state (optional) gives the coordinates which tell the carbonyl of a C terminus without its OXT
// from a hydroxyl.
bool md_system_infer_residue_flags(struct md_system_t* sys, const struct md_system_state_t* state);

// CHAIN BREAKS. A polymer chain is broken between two consecutive residues of it (amino acids, nucleotides) whose
// backbones are not linked: no peptide bond from the C of the one to the N of the next, no phosphodiester bond from the
// O3' of the one to the P of the next, no bond between the backbone beads of a coarse grained chain. The residues a
// crystal structure could not resolve leave such breaks, as does a model built in pieces.
// Where the chemistry says so unambiguously the chain ends there instead: a free carboxyl group (O and OXT) before it,
// a free amino group (two or three hydrogens on the N) after it, a free amino acid after it (a lone residue with a free
// carboxyl group of its own, which a PDB file can number on after the chain, as 1LAF does its bound arginine), a
// nucleotide without its P (or phosphate bead) after it (a 5' end), or a free hydroxyl (with its hydrogen) on either
// side. A crystal structure rarely has the hydrogens and often leaves out the OXT, so the ends of the pieces of a
// broken chain look like termini to the classification of the residues: the residues on either side of a break lose
// the terminal flags it gave them (MD_RESIDUE_FLAG_TERMINAL_BEG, _END, and those of their atoms).
// The bonds decide, not the sequence numbers: an insertion code (52, 52A) or numbering that skips (antibodies, the
// chymotrypsin numbering of proteases) does not break a chain whose backbone is bonded, and numbering that runs on
// does not join one whose backbone is not. Atoms without bonds (a structure file read without inferring them) have no
// backbones to break; beads without bonds are linked by consecutive sequence numbers. md_system_infer (MD_SYSTEM_INFER_BACKBONE_BIT) marks the residues on either side of each break,
// MD_RESIDUE_FLAG_BREAK_AFTER the one before it and MD_RESIDUE_FLAG_BREAK_BEFORE the one after it, and ends the
// backbone ranges there (md_protein_backbone_data_t, md_nucleic_backbone_data_t).

// The atoms a link over the break after the residue res_idx (MD_RESIDUE_FLAG_BREAK_AFTER) would join: the C and
// the N of a peptide bond, the O3' and the P of a phosphodiester bond, the backbone beads of coarse grained residues.
// False, and out_atoms -1, when there is no break after it.
bool md_system_chain_break_atoms(md_atom_idx_t out_atoms[2], const struct md_system_t* sys, size_t res_idx);

// The atoms segment seg_idx of the protein backbone stands for: those of its residue; for a segment of none (a system
// without residues, a residue which is a whole chain) those of its amino acid by the roles of the atoms: from its CA,
// everything bonded to it short of the backbone of another segment, a disulfide bond and a coordination bond, which is
// its N, C and O, its side chain, and the hydrogens and terminal groups on these. Writes at most cap of them to out and
// returns how many there are (at most 512 for a segment of no residue).
size_t md_system_protein_segment_atoms(md_atom_idx_t out[], size_t cap, const struct md_system_t* sys, size_t seg_idx);

// The per segment quantities of a protein backbone (see md_protein_backbone_data_t), computed into the caller's arrays
// of one entry per segment (capacity at least backbone->segment.count). They do not allocate. The secondary structure
// has its own header, md_dssp.h.

// Backbone angles (phi, psi) from the coordinates of a frame. A trace has none: its segments are left 0.
bool md_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const struct md_float3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone);

// Ramachandran type (General / Glycine / Proline / Preproline, as MolProbity assigns them) from the residue names of
// sys->protein_backbone, or for segments without a residue (res_idx -1) from their bonds. Writes min(capacity, segment count) types; the rest of the capacity is left UNKNOWN, as are the
// segments of traces, which have no phi and psi.
bool md_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const struct md_system_t* sys);

// THE BACKBONE OF A STATE. The angles and the secondary structure depend on the coordinates, so they belong to the
// state they were computed from and live in its attributes (md_system_state_t.attributes), one value per segment of
// sys->protein_backbone:
//   MD_BACKBONE_ANGLE_PATH                 F32 x 2 (phi, psi), radians
//   MD_BACKBONE_SECONDARY_STRUCTURE_PATH   I32, md_secondary_structure_t
// The same paths below a run (in the system's attributes) hold them for every frame of the run, with the frame as
// the first axis, and md_system_extract_frame slices a frame of them into a state when asked for them.
// A state carries them when its producer computed them, like the velocities of a TRR frame: a consumer handles their
// absence.
#define MD_BACKBONE_ANGLE_PATH                  "backbone/angle"
#define MD_BACKBONE_SECONDARY_STRUCTURE_PATH    "backbone/secondary_structure"

// The backbone of a state, NULL when it carries none or what it carries does not fit the backbone of sys
const md_backbone_angles_t*     md_system_state_backbone_angles(const struct md_system_state_t* state, const struct md_system_t* sys);
const md_secondary_structure_t* md_system_state_secondary_structure(const struct md_system_state_t* state, const struct md_system_t* sys);

// Storage in the state to write the backbone into, created when missing (or of another size) and reused otherwise,
// so that a producer which writes every frame does not allocate every frame. Each call marks the attribute as
// changed (md_attributes_touch): call it when about to write. NULL for a state without an allocator (a view) or a
// system without a protein backbone.
md_backbone_angles_t*     md_system_state_backbone_angles_write(struct md_system_state_t* state, const struct md_system_t* sys);
md_secondary_structure_t* md_system_state_secondary_structure_write(struct md_system_state_t* state, const struct md_system_t* sys);

// Computes the backbone angles and secondary structure from the state's own coordinates and stores them in it.
// False when the state has no coordinates, or the system no protein backbone.
bool md_system_state_backbone_compute(struct md_system_state_t* state, const struct md_system_t* sys);

// ### INFERENCE ###

// The backbones of the system (md_protein_backbone_data_t, md_nucleic_backbone_data_t) and the breaks of its chains:
// ranges of at least 3 residues along each chain of a polymer, ending where the chain breaks (CHAIN BREAKS above).
//   ATOMS  From the covalent bonds alone, whatever the atoms and residues are called (BACKBONES FROM THE BONDS in
//          md_backbone.c): an atomistic system without bonds has no backbones. A segment is of the residue of its CA (C4')
//          where the system has residues, a range of the part of its residues.
//   BEADS  The traces of a coarse grained model by the backbone beads its residues were given (Martini's BB, BB1 to
//          BB3), along each strand of a polymer (the peptide and nucleic acid parts). Linked by their bonds, else by
//          consecutive sequence numbers.
// The ranges come in the order of their first atom. Requires the bonds and the roles of the atoms
// (md_system_infer_residue_flags), and the parts for the beads. Into out (allocated with alloc); the system is read
// and never written. md_system_infer (MD_SYSTEM_INFER_BACKBONE_BIT, after the parts) gives the system the backbones
// and marks the breaks (md_system_set_chain_breaks).
typedef struct md_backbone_data_t {
    md_protein_backbone_data_t protein;
    md_nucleic_backbone_data_t nucleic;
    size_t num_breaks;
    md_residue_idx_t* breaks;       // [num_breaks] The residues after which a chain breaks, ascending
} md_backbone_data_t;

bool md_backbone_infer(md_backbone_data_t* out, const struct md_system_t* sys, struct md_allocator_i* alloc);
void md_backbone_data_free(md_backbone_data_t* bb, struct md_allocator_i* alloc);

// Marks the breaks of the chains of the system (CHAIN BREAKS above): residue breaks[k] MD_RESIDUE_FLAG_BREAK_AFTER
// and the one after it MD_RESIDUE_FLAG_BREAK_BEFORE, the flags of every other residue cleared. A break is no end of a
// chain: the residues on either side of one lose the terminal flags their classification gave them
// (MD_RESIDUE_FLAG_TERMINAL_BEG, _END, and those of their atoms).
void md_system_set_chain_breaks(struct md_system_t* sys, const md_residue_idx_t breaks[], size_t num_breaks);

#ifdef __cplusplus
}
#endif

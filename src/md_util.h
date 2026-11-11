#pragma once

#include <md_system.h>
#include <md_element.h>

#include <core/md_str.h>
#include <core/md_vec.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>

struct md_allocator_i;
struct md_executor_t;
struct md_bitset_t;
struct md_bitmask_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MD_UTIL_INFER_NONE                = 0,
    MD_UTIL_INFER_COLOR_BIT           = 1u << 1,
    MD_UTIL_INFER_BOND_BIT            = 1u << 2,
    MD_UTIL_INFER_PART_BIT            = 1u << 3,     // The chains and the parts with their types (md_chain_infer, _parts)
    MD_UTIL_INFER_BACKBONE_BIT        = 1u << 4,
    MD_UTIL_INFER_MOLECULE_BIT        = 1u << 5,
    MD_UTIL_INFER_RESIDUE_BIT         = 1u << 6,     // The kinds of the residues and the roles of their atoms (md_util_system_infer_residue_flags), after the bonds which tell them
	//MD_UTIL_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_UTIL_INFER_UNWRAP_MOLECULE_BIT = 1u << 7,
    MD_UTIL_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h
    MD_UTIL_INFER_LIPID_BIT           = 1u << 9,     // The heads and tails of the lipids (md_lipid_infer), and the kinds of their residues
    MD_UTIL_INFER_FRAGMENT_BIT        = 1u << 10,    // Ring systems, functional groups, backbones (md_util_system_infer_fragments, md_fragment.h), last

    MD_UTIL_INFER_ALL                 = -1,
} md_infer_flags_t;

MD_ENUM_FLAGS(md_infer_flags_t)

bool md_util_resname_dna(md_str_t str);
bool md_util_resname_rna(md_str_t str);
bool md_util_resname_acidic(md_str_t str);
bool md_util_resname_basic(md_str_t str);
bool md_util_resname_neutral(md_str_t str);
bool md_util_resname_water(md_str_t str);
bool md_util_resname_ion(md_str_t str);         // The common monatomic ions (NA, CL, ...)
bool md_util_resname_lipid(md_str_t str);      // The lipids of the Martini lipidomes (their names are mostly the atomistic ones)
bool md_util_resname_hydrophobic(md_str_t str);
bool md_util_resname_amino_acid(md_str_t str);
bool md_util_resname_nucleotide(md_str_t str);

void md_util_system_extract_xyzw_from_mask(md_f32x4_t* out_xyzw, const struct md_bitset_t* mask, const md_system_t* sys, const md_system_state_t* state);

// Gives the atoms of the amino acids and nucleotides of an atomistic system their roles (MD_ATOM_FLAG_BACKBONE,
// _SIDE_CHAIN, _NUCLEOSIDE, _NUCLEOBASE, _TERMINAL_BEG, _END), found from the covalent bonds alone, whatever the atoms
// are called (POLYMER ATOMS in md_util.c), also in a system without residues. Then classifies the residues
// (md_residue_kind_t: amino acid, nucleotide, water, ion, lipid). A residue which holds the backbone of one amino
// acid or nucleotide of a chain (a standard residue by name, or linked to another residue by its backbone) is one,
// MD_RESIDUE_FLAG_RESOLVED, with its terminal groups marked: the link keeps out the ligands which hold the atoms of
// one (GTP, ATP, NAD) while modified residues in a chain (MSE, PSU) come in. A residue which is only named like an
// amino acid or nucleotide gets the kind alone. A residue of coarse grained beads which is a molecule of the Martini
// force fields by its name and all its bead names (a lipid, water, an ion, a sugar, ...) is of the kind of that
// molecule; residues a predefined coarse grained type already classified as water, ion or lipid are left alone.
// Requires the bonds. state (optional) gives the coordinates which tell the carbonyl of a C terminus without its OXT
// from a hydroxyl.
bool md_util_system_infer_residue_flags(md_system_t* sys, const md_system_state_t* state);

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
// backbones to break; beads without bonds are linked by consecutive sequence numbers. md_util_system_infer (MD_UTIL_INFER_BACKBONE_BIT) marks the residues on either side of each break,
// MD_RESIDUE_FLAG_BREAK_AFTER the one before it and MD_RESIDUE_FLAG_BREAK_BEFORE the one after it, and ends the
// backbone ranges there (md_protein_backbone_data_t, md_nucleic_backbone_data_t).

// The atoms a link over the break after the residue res_idx (MD_RESIDUE_FLAG_BREAK_AFTER) would join: the C and
// the N of a peptide bond, the O3' and the P of a phosphodiester bond, the backbone beads of coarse grained residues.
// False, and out_atoms -1, when there is no break after it.
bool md_util_chain_break_atoms(md_atom_idx_t out_atoms[2], const md_system_t* sys, size_t res_idx);

// The atoms segment seg_idx of the protein backbone stands for: those of its residue; for a segment of none (a system
// without residues, a residue which is a whole chain) those of its amino acid by the roles of the atoms: from its CA,
// everything bonded to it short of the backbone of another segment, a disulfide bond and a coordination bond, which is
// its N, C and O, its side chain, and the hydrogens and terminal groups on these. Writes at most cap of them to out and
// returns how many there are (at most 512 for a segment of no residue).
size_t md_util_protein_segment_atoms(md_atom_idx_t out[], size_t cap, const md_system_t* sys, size_t seg_idx);

// The per segment quantities of a protein backbone (see md_protein_backbone_data_t), computed into the caller's arrays
// of one entry per segment (capacity at least backbone->segment.count). They do not allocate. The secondary structure
// has its own header, md_dssp.h.

// Backbone angles (phi, psi) from the coordinates of a frame. A trace has none: its segments are left 0.
bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const md_f32x3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone);

// Ramachandran type (General / Glycine / Proline / Preproline, as MolProbity assigns them) from the residue names of
// sys->protein_backbone, or for segments without a residue (res_idx -1) from their bonds. Writes min(capacity, segment count) types; the rest of the capacity is left UNKNOWN, as are the
// segments of traces, which have no phi and psi.
bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const struct md_system_t* sys);

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
const md_backbone_angles_t*     md_util_state_backbone_angles(const md_system_state_t* state, const struct md_system_t* sys);
const md_secondary_structure_t* md_util_state_secondary_structure(const md_system_state_t* state, const struct md_system_t* sys);

// Storage in the state to write the backbone into, created when missing (or of another size) and reused otherwise,
// so that a producer which writes every frame does not allocate every frame. Each call marks the attribute as
// changed (md_attributes_touch): call it when about to write. NULL for a state without an allocator (a view) or a
// system without a protein backbone.
md_backbone_angles_t*     md_util_state_backbone_angles_write(md_system_state_t* state, const struct md_system_t* sys);
md_secondary_structure_t* md_util_state_secondary_structure_write(md_system_state_t* state, const struct md_system_t* sys);

// Computes the backbone angles and secondary structure from the state's own coordinates and stores them in it.
// False when the state has no coordinates, or the system no protein backbone.
bool md_util_state_backbone_compute(md_system_state_t* state, const struct md_system_t* sys);

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_bonds(struct md_bitmask_t* mask, const struct md_system_t* sys, size_t extent, const struct md_bitset_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_radius(struct md_bitmask_t* mask, const md_system_state_t* state, double radius, const struct md_bitset_t* viable_mask);

// The type of an atom of the name and element, added when the system has none, with the mass, size and color of the
// element. map, which the caller keeps for a load (zeroed, with an allocator), finds a type it has given before in one
// lookup rather than a search of the types.
md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, md_str_t name, md_atomic_number_t z);

// Attempts to generate missing data such as covalent bonds, chains, secondary structures, backbone angles etc.
// Infers the derivable parts of a system (covalent bonds, rings, molecules, backbones, hydrogen
// bond roles) from the supplied state, and records that state as sys->reference. A system without residues has its
// backbones found from the covalent bonds alone (md_protein_backbone_data_t).
//
// The state is an explicit parameter rather than read off the system so that the recorded
// reference is by construction the input which produced the topology. Requires coordinates for
// anything geometric; a state without them (a topology only format such as PSF) still yields the
// purely graph derived parts, provided the bonds were supplied by the format.
//
// The phases run one after the other, each on the threads of the executor where it spreads its work (NULL: the calling
// thread); the result is the same whatever the threads.
bool md_util_system_infer(struct md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags, const struct md_executor_t* executor);

// INFERENCE WHILE THE SYSTEM IS IN USE. Inference reads the CORE of a system (the atoms, their types, the residues, the
// instances, the assemblies) and writes the rest, the DERIVED (md_system_fields.h). So a system can be shown and used
// while it is inferred on another thread, without a copy of it:
//
//   md_util_system_infer_prepare     The work: of the system, what inference writes, copied into alloc; the coordinates
//                                    in state, which become its reference; the per atom arrays of the core shared. On
//                                    the thread which owns the system.
//   md_util_system_infer_run         Inference on the work: on any thread, while the system is read (and its core
//                                    is not written: what changes the topology waits for the work, or discards it).
//                                    All but MD_UTIL_INFER_COLOR_BIT, which is of the types of the atoms: the system
//                                    is given its colors before, md_util_system_infer with that bit alone.
//   md_util_system_infer_publish     The derived arrays of the work become the system's, which keeps its core, its
//                                    attributes and its force field. On the thread which owns the system, nothing
//                                    else reading it. False, and nothing done, when the topology of the system changed
//                                    after the work was prepared: the work is of a system which is no more.
//
// The arrays of the work are in alloc, which the caller keeps with the system once published (they become the system's)
// and frees otherwise (the work is never md_system_free'd: it shares arrays of the system). The types of the atoms are
// copied too (a few per system), so that the system's may be changed meanwhile (their colors, their radii): they stay
// the system's when published.
typedef struct md_util_infer_work_t {
    struct md_system_t sys;     // Inferred on: md_util_system_infer_run
    uint64_t base_version;      // The topology version of the system it was prepared of
} md_util_infer_work_t;

bool md_util_system_infer_prepare(md_util_infer_work_t* work, const struct md_system_t* sys, const md_system_state_t* state, struct md_allocator_i* alloc);
bool md_util_system_infer_run(md_util_infer_work_t* work, md_infer_flags_t flags, const struct md_executor_t* executor);
bool md_util_system_infer_publish(struct md_system_t* sys, md_util_infer_work_t* work);

#ifdef __cplusplus
}
#endif

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
    MD_UTIL_INFER_PART_BIT            = 1u << 3,     // The chains and the parts with their types (md_util_system_infer_chains, _parts)
    MD_UTIL_INFER_BACKBONE_BIT        = 1u << 4,
    MD_UTIL_INFER_MOLECULE_BIT        = 1u << 5,
    MD_UTIL_INFER_RESIDUE_BIT         = 1u << 6,     // The kinds of the residues and the roles of their atoms (md_util_system_infer_residue_flags), after the bonds which tell them
	//MD_UTIL_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_UTIL_INFER_UNWRAP_MOLECULE_BIT = 1u << 7,
    MD_UTIL_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h
    MD_UTIL_INFER_LIPID_BIT           = 1u << 9,     // The heads and tails of the lipids (md_util_system_infer_lipids)
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

// The lipids of the system (sys->lipid, see md_lipid_data_t): every lipid molecule, its heads and tails, and the atoms
// of the tails (MD_ATOM_FLAG_LIPID_TAIL). A molecule (a run of consecutive residues bonded to each other, none of them an amino
// acid, a nucleotide, water or an ion) whose residues are lipids by name gets the tails found in it. One whose
// residues are not is a lipid by its structure when it has a polar head and a hydrocarbon chain of at least 8 carbons
// or the skeleton of a sterol (atomistic molecules only), and its residues become MD_RESIDUE_KIND_LIPID: lipids of
// names nothing knows (the 3 letter names of PDB files, the codes of the PDB, lipids split into several residues).
// Tails are found from the bonds alone for atoms, and from the names of the beads (the conventions of the Martini
// lipidomes) for coarse grained lipids. Requires the classified residues, the bonds for the heads and tails and for
// the molecules of several residues (without them a lipid is its residue), and for sterols the rings.
bool md_util_system_infer_lipids(md_system_t* sys);

// Instances from the chain ids of the residues, for a format which gives a chain id per residue and no units of its
// own (a PDB file): one instance per run of consecutive residues of the same id, its auth_id the id, and a new one
// where res_unit_begins says the data ends one before the residue (a TER record of a PDB file), whatever the ids: two
// chains of one id in a row are two (NULL: by the ids alone). A run without an id is an instance of no id (the
// instances are back to back). Replaces the instances of the system. The instances are what the data gives
// (md_instance_data_t): where it gives no chain ids at all the system has none.
bool md_util_system_instances_from_chain_ids(md_system_t* sys, const md_str_t res_chain_id[], const bool res_unit_begins[]);

// The chains of the system (sys->chain, see md_chain_data_t) from the author ids of its instances: one per id. A system
// without ids has no chains. Done by md_util_system_infer (MD_UTIL_INFER_PART_BIT) and by the loaders which give the
// ids; whatever changes the instances after calls it again.
bool md_util_system_infer_chains(md_system_t* sys);

// The parts of the system (sys->part, see md_part_data_t) and their types from its instances, its bonds and the kinds
// of its residues. A part never reaches over an instance; a system without instances is one span to divide. A strand
// goes on over a gap (two consecutive residues of it not bonded) within one instance the data names (an id or an author
// id: an asym, a chain of a PDB file; not a molecule of a topology) while the sequence ids increase, and elsewhere in
// a system without bonds while they are consecutive, unless the chain ends there (CHAIN BREAKS). A part which is
// all of an instance with a given part type is of that type; the others are of inferred types
// (MD_PART_TYPE_FLAG_INFERRED), one per kind and sequence of residue names, which replace those of the last inference.
// A given type of unknown kind is classified from its first part. Done by md_util_system_infer (MD_UTIL_INFER_PART_BIT,
// after the bonds); whatever changes the instances, the bonds or the kinds calls it after.
bool md_util_system_infer_parts(md_system_t* sys);

// The backbones of the system (sys->protein_backbone, sys->nucleic_backbone, see md_protein_backbone_data_t and
// md_nucleic_backbone_data_t), found anew: ranges of at least 3 residues along each chain of a polymer, ending where the
// chain breaks (CHAIN BREAKS), whose residues it marks MD_RESIDUE_FLAG_BREAK_AFTER and _BEFORE (cleared first).
//   ATOMS  From the covalent bonds alone, whatever the atoms and residues are called (BACKBONES FROM THE BONDS in
//          md_util.c): an atomistic system without bonds has no backbones. A segment is of the residue of its CA (C4')
//          where the system has residues, a range of the part of its residues.
//   BEADS  The traces of a coarse grained model by the backbone beads its residues were given (Martini's BB, BB1 to
//          BB3), along each strand of a polymer (the peptide and nucleic acid parts). Linked by their bonds, else by
//          consecutive sequence numbers (state unused).
// The ranges come in the order of their first atom. Done by md_util_system_infer (MD_UTIL_INFER_BACKBONE_BIT, after the
// parts); requires the bonds, and the parts for the beads.
bool md_util_system_infer_backbones(md_system_t* sys, const md_system_state_t* state);

// The kind of a residue is that of its type (md_residue_type_data_t): one per name and kind of residue, with the
// standard monomer it is a form of (md_monomer_from_residue_name) and a default color: that of the monomer, that of the
// element of a monatomic ion, else one made from its name. md_util_system_set_residue_kind makes the residue of the
// type of its name and the kind, added where the system has none; the types the system has already are kept (with
// the colors they were given), also one no residue is of any more.
void md_util_system_set_residue_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind);

// Gives every residue of no type yet (one never classified) the type of its name and kind OTHER. md_util_system_infer
// calls it last.
bool md_util_system_infer_residue_types(md_system_t* sys);

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

// Bonds of origin MD_BOND_ORIGIN_INFERRED, replacing those of an earlier inference and keeping every other: an inferred
// bond never duplicates one of another origin. Atoms are bonded by their covalent radii (derived from the element).
// The atoms of a metallic phase (a metal in bulk, at a surface or in a particle) are bonded to the metal atoms of their
// first coordination shell by metallic bonds and to nothing else, and a halogen bonded to a non-metal coordinates no
// metal (METALLIC PHASES in md_util.c). Coarse grained systems (md_system_is_coarse_grained) are bonded as their
// Martini topologies are: a residue whose beads are those of a residue of the Martini force fields has the bonds of its
// topology, any other amino acid or nucleotide a minimum spanning tree of its beads with its rings. The backbone bead of
// each residue is bonded to that of the next unless the numbering starts over, and cysteines are bridged (see
// infer_coarse_grained_bonds in md_util.c). The pairs within reach are found on the threads of the executor (NULL: the
// calling thread), and the bonds come out the same whatever the threads.
void md_util_infer_covalent_bonds(md_bond_data_t* out_bond, const md_system_state_t* state, const md_system_t* sys, struct md_allocator_i* alloc, const struct md_executor_t* executor);

// Infers the covalent bonds of the system (md_util_infer_covalent_bonds) anew. The bonds a user added by hand and those
// the file gives (MD_BOND_ORIGIN_FILE: the CONECT records of a PDB file) are kept, and so are the non-covalent bonds:
// all other bonds of the graph are replaced, as are the metallic bonds. The atoms of the metallic bonds are marked
// MD_ATOM_FLAG_METALLIC, and only they.
void md_util_system_infer_covalent_bonds(md_system_t* sys, const md_system_state_t* state, const struct md_executor_t* executor);

// Marks every bond between a metal and a non-metal MD_BOND_KIND_COORDINATE, whatever its origin, unless its kind was
// given (MD_BOND_FLAG_KIND_GIVEN), and returns how many it marked. md_util_system_infer does it for each system; a bond
// added afterwards (by hand) wants it too.
size_t md_util_system_infer_coordination(md_system_t* sys);

// Marks the bonds of a coarse grained topology which hold a protein or nucleic acid in shape without being bonds of it
// MD_BOND_KIND_ELASTIC, and returns how many it marked: an elastic network (ElNeDyn, the rubber bands of martinize) and
// the restraints of secondary structure, when they are written as plain bonds (bond type 1), as they are as often as
// not. Between two residues of a chain (amino acids or nucleotides, beads) a Martini topology has the link from the last
// backbone bead of one to the first of the next (BB to BB, BB3 to BB1, and BB2 to BB1 as well in Martini 3 RNA), a
// disulfide (the SC1 of two cysteines), and the bond which closes a cyclic chain; every other bond between beads of
// two residues of a chain is elastic. Bonds within a residue, and to anything other than amino acids and nucleotides,
// are left as they are. Needs the kinds of the residues and the backbone beads (md_util_system_infer_residue_flags),
// and rebuilds the connectivity when it marks any. The topology readers (md_tpr, md_itp) call it.
size_t md_util_system_infer_elastic_bonds(md_system_t* sys);

// ### BONDS FORMED AND BROKEN ###
// The bonds of a state against the bond graph of the system, which is the reference topology: the bonds of the graph
// broken in the state, and the bonds formed in it which the system does not have. What the bonds of the state are is
// the graph without the broken bonds and with the formed ones; the system is not changed. Only the bonds of a scope are
// looked at, those with an atom in it (so a bond to an atom just outside of it forms too): outside of it the bonds of
// the state are those of the graph.
//
// The rules are those of md_util_infer_covalent_bonds, with a band around their cutoffs so that a bond does not flicker
// from frame to frame as it vibrates near one: a bond of the graph breaks only beyond (1 + tolerance) times its cutoff,
// and a pair forms a bond only within (1 - tolerance) times it. Measured against the graph rather than an earlier state,
// so the bonds of a state are those of that state alone, whatever was evaluated before it.
//   - Only the bonds of the breakable origins break. A bond a user added stays, and by default so does one of a force
//     field topology, which holds throughout a classical simulation.
//   - A formed bond is limited by what its atoms can take (md_util_infer_covalent_bonds prunes the same way). A bond of
//     the graph counts against that while it is within its cutoff, but not once it is stretched beyond it on its way to
//     breaking: a bond forming as another breaks (a proton moving, the transition state of a substitution) shows both
//     for the frames it takes, while two atoms which clash in the reference do not bond in every frame.
//   - A pair with a non-covalent bond in the system forms no bond: what was given of it stands.
//   - The metallic bonds neither break nor form: two atoms of a metallic phase (MD_ATOM_FLAG_METALLIC) are left out.
//     What touches the phase has its bonds to it here only, formed in the frames where it is close enough, as it
//     moves across a surface.
//   - A halogen bonded to a non-metal in the state, by a bond of the graph which is not broken or by one formed,
//     forms no coordination bond.
//   - A coarse grained system has no bonds formed or broken: its bonds come from templates, not distances.
typedef struct md_bond_delta_t {
    md_bitset_t      broken;        // Bonds of the graph broken in the state: indices into sys->bond, in [0, count)
    size_t           num_formed;
    md_atom_pair_t*  formed_pair;   // [num_formed] Pairs bonded in the state which have no bond in the system, each the
                                    // lower atom first, sorted
    md_bond_flags_t* formed_flags;  // [num_formed] Covalent or coordinate, of order unknown and origin inferred
} md_bond_delta_t;

#define MD_BOND_DELTA_TOLERANCE 0.1f
#define MD_BOND_DELTA_BREAKABLE ((1u << MD_BOND_ORIGIN_FILE) | (1u << MD_BOND_ORIGIN_INFERRED))

typedef struct md_bond_delta_desc_t {
    const md_bitset_t* scope;       // The atoms whose bonds are looked at, NULL for all of them
    float    tolerance;             // Half the width of the band around the cutoffs, as a fraction of them: [0, 0.5).
                                    // MD_BOND_DELTA_TOLERANCE is what to use without a reason for another
    uint32_t breakable;             // The origins whose bonds may break, a mask of 1 << md_bond_origin_t:
                                    // MD_BOND_DELTA_BREAKABLE (a file's and those inferred) without a reason for another
} md_bond_delta_desc_t;

// The delta of a state, written to out and allocated with alloc. False, with an empty delta, when the state has no
// coordinates (or not those of the system) or the system is coarse grained.
bool md_util_bond_delta(md_bond_delta_t* out, const md_system_t* sys, const md_system_state_t* state, const md_bond_delta_desc_t* desc, struct md_allocator_i* alloc);
void md_bond_delta_free(md_bond_delta_t* delta, struct md_allocator_i* alloc);

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_bonds(struct md_bitmask_t* mask, const struct md_system_t* sys, size_t extent, const struct md_bitset_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_radius(struct md_bitmask_t* mask, const md_system_state_t* state, double radius, const struct md_bitset_t* viable_mask);

// Infer rings formed by covalent bonds: the rings of 3 to 7 atoms which are not sums of smaller rings (sys->ring: the
// smallest set of smallest rings, made symmetric), and their flags (md_ring_t::flags)
// (md_util_system_infer_ring_flags). Only atoms
// (MD_PARTICLE_ATOM) make rings: the cycles of coarse grained beads and virtual sites are not searched.
bool md_util_system_infer_rings(md_system_t* sys);

// The rings of a bond graph as md_util_system_infer_rings finds them, of every atom: bond, with its connectivity
// (md_bond_build_connectivity), over num_atoms atoms. Pushed to out_rings, cleared first and allocated with alloc.
void md_util_infer_rings(md_ring_data_t* out_rings, const md_bond_data_t* bond, size_t num_atoms, struct md_allocator_i* alloc);

// The flags of the rings (md_ring_t::flags, md_ring_flags_t) from the rings, the bonds and their flags, and the atoms of
// sugar rings in canonical order (MD_RING_FLAG_SUGAR). Call it again when the bond orders or aromaticity change
// (md_chem_perceive does).
void md_util_system_infer_ring_flags(md_system_t* sys);

// The molecules of the system (md_molecule_data_t), see there for what holds one together.
// For coarse grained systems (any particle a bead, md_system_is_coarse_grained), beads the bonds leave disconnected are
// joined through the residue hierarchy instead: each bead to its residue's backbone (or first) bead, and those
// anchors along consecutive polymer residues. No coordinates are consulted and sys->bond is not modified.
// The molecules get their types too (md_molecule_type_data_t): the given type of the instance a molecule is all of (a
// molecule of a topology), else an inferred one, which replace those of the last inference. Requires the parts.
bool md_util_system_infer_molecules(md_system_t* sys);

// The type of an atom of the name and element, added when the system has none, with the mass, size and color of the
// element. map, which the caller keeps for a load (zeroed, with an allocator), finds a type it has given before in one
// lookup rather than a search of the types.
md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, md_str_t name, md_atomic_number_t z);

// Identify particle types within the system from the atom labels and the names of their residues. Coarse grained beads
// come from predefined tables. Martini beads of amino acids (BB, SC#) and nucleotides (BB#, SC#) are recognized in any
// residue named as one, and get the masses and sizes of the version of Martini (2, polarizable 2 or 3) the names of
// the system's residues and beads tell (Martini 3 when they tell nothing); a massless bead is a virtual site. Atoms of
// one label with a different element, mass, size or particle kind are of different types.
void md_util_system_infer_particle_types(md_system_t* sys, const md_str_t atom_labels[], const struct md_executor_t* ex);
// The same, with the labels as labels (8 characters): what a loader keeps per atom at half the size of a string.
void md_util_system_infer_particle_types_from_labels(md_system_t* sys, const md_label_t atom_labels[], const struct md_executor_t* ex);

// Applies the predefined particle type tables (coarse grained beads) to atoms that already have a type,
// for loaders that know more about their atoms than a name does (a tpr carries the mass and the LJ
// parameters of every bead). Only types without an element (z == 0) are considered. The table gives
// the atoms their role (backbone, side chain), the residues their kind (amino acid, water) where
// they have none, and the types their particle kind where the loader did not say. Mass and radius
// are left alone: they are the loader's.
void md_util_system_augment_particle_types(md_system_t* sys);

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

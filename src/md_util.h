#pragma once

#include <md_system.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>

struct md_allocator_i;
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
	//MD_UTIL_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_UTIL_INFER_UNWRAP_MOLECULE_BIT = 1u << 7,
    MD_UTIL_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h
    MD_UTIL_INFER_LIPID_BIT           = 1u << 9,     // The heads and tails of the lipids (md_util_system_infer_lipids)

    MD_UTIL_INFER_ALL                 = -1,
} md_infer_flags_t;

ENUM_FLAGS(md_infer_flags_t)

// Access to the static arrays (preserved for direct access)
const str_t* md_util_element_symbols(void);
const str_t* md_util_element_names(void);
const float* md_util_element_vdw_radii(void);

// Element functions (now calling new atomic number API internally)
md_element_t md_util_element_lookup(str_t element_str, bool ignore_case);
str_t md_util_element_symbol(md_element_t element);
str_t md_util_element_name(md_element_t element);
float md_util_element_vdw_radius(md_element_t element);
float md_util_element_covalent_radius(md_element_t element);
float md_util_element_atomic_mass(md_element_t element);
int   md_util_element_max_valence(md_element_t element);
uint32_t md_util_element_cpk_color(md_element_t element);

bool md_util_resname_dna(str_t str);
bool md_util_resname_rna(str_t str);
bool md_util_resname_acidic(str_t str);
bool md_util_resname_basic(str_t str);
bool md_util_resname_neutral(str_t str);
bool md_util_resname_water(str_t str);
bool md_util_resname_lipid(str_t str);      // The lipids of the Martini lipidomes (their names are mostly the atomistic ones)
bool md_util_resname_hydrophobic(str_t str);
bool md_util_resname_amino_acid(str_t str);
bool md_util_resname_nucleotide(str_t str);

void md_util_system_extract_xyzw_from_mask(vec4_t* out_xyzw, const struct md_bitset_t* mask, const md_system_t* sys, const md_system_state_t* state);

// Classifies the residues (md_residue_kind_t: amino acid, nucleotide, water, ion, lipid) from their names, atoms
// and bonds. A residue of coarse grained beads which is a molecule of the Martini force fields by its name and all
// its bead names (a lipid, water, an ion, a sugar, ...) is of the kind of that molecule. An amino acid or nucleotide whose backbone atoms are found by name and verified by their bonds, and which
// is either a standard residue by name or linked into a chain (a peptide or phosphodiester bond to another
// residue), is MD_RESIDUE_FLAG_RESOLVED: its atoms are given their roles (MD_ATOM_FLAG_BACKBONE, _SIDE_CHAIN,
// ...) and its terminal groups are marked. The link keeps out the ligands which name their atoms like a residue
// (GTP, ATP, NAD) while modified residues in a chain (MSE, PSU) come in. A residue which is only named like an
// amino acid or nucleotide gets the kind alone. Residues a predefined coarse grained type already classified as
// water, ion or lipid are left alone. Requires the bonds.
bool md_util_system_infer_residue_flags(md_system_t* sys);

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
// own (a PDB file): one instance per run of consecutive residues of the same id, its auth_id the id. A run without an
// id is an instance of no id (the instances are back to back). Replaces the instances of the system. The instances are
// what the data gives (md_instance_data_t): where it gives no chain ids at all the system has none.
bool md_util_system_instances_from_chain_ids(md_system_t* sys, const str_t res_chain_id[]);

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

// The residue types of the system (sys->residue.type, see md_residue_type_data_t): one per name and kind of
// residue, with the standard monomer it is a form of (md_monomer_from_residue_name) and a default color: that of the
// monomer, that of the element of a monatomic ion, else one made from its name. The types the system has already are
// kept (with the colors they were given), the residues which match none get new ones. Requires the classified
// residues, lipids included; md_util_system_infer calls it last.
bool md_util_system_infer_residue_types(md_system_t* sys);

size_t md_util_element_from_mass(md_element_t out_element[], const float in_mass[], size_t count);

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
// does not join one whose backbone is not. Atoms without bonds (a structure file read without inferring them) are
// linked by their distance, the shortest way round a periodic cell, and without coordinates by consecutive sequence
// numbers. md_util_system_infer (MD_UTIL_INFER_BACKBONE_BIT) marks the residues on either side of each break,
// MD_RESIDUE_FLAG_BREAK_AFTER the one before it and MD_RESIDUE_FLAG_BREAK_BEFORE the one after it, and ends the
// backbone ranges there (md_protein_backbone_data_t, md_nucleic_backbone_data_t).

// The atoms a link over the break after the residue res_idx (MD_RESIDUE_FLAG_BREAK_AFTER) would join: the C and
// the N of a peptide bond, the O3' and the P of a phosphodiester bond, the backbone beads of coarse grained residues.
// False, and out_atoms -1, when there is no break after it.
bool md_util_chain_break_atoms(md_atom_idx_t out_atoms[2], const md_system_t* sys, size_t res_idx);

// The per segment quantities of a protein backbone (see md_protein_backbone_data_t), computed into the caller's arrays
// of one entry per segment (capacity at least backbone->segment.count). They do not allocate. The secondary structure
// has its own header, md_dssp.h.

// Backbone angles (phi, psi) from the coordinates of a frame. A trace has none: its segments are left 0.
bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone);

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
// infer_coarse_grained_bonds in md_util.c).
void md_util_infer_covalent_bonds(md_bond_data_t* out_bond, const md_system_state_t* state, const md_system_t* sys, struct md_allocator_i* alloc);

// Infers the covalent bonds of the system (md_util_infer_covalent_bonds) anew. The bonds a user added by hand and those
// the file gives (MD_BOND_ORIGIN_FILE: the CONECT records of a PDB file) are kept, and so are the non-covalent bonds:
// all other bonds of the graph are replaced, as are the metallic bonds. The atoms of the metallic bonds are marked
// MD_ATOM_FLAG_METALLIC, and only they.
void md_util_system_infer_covalent_bonds(md_system_t* sys, const md_system_state_t* state);

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
md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, str_t name, md_atomic_number_t z);

// Identify particle types within the system from the atom labels and the names of their residues. Coarse grained beads
// come from predefined tables. Martini beads of amino acids (BB, SC#) and nucleotides (BB#, SC#) are recognized in any
// residue named as one, and get the masses and sizes of the version of Martini (2, polarizable 2 or 3) the names of
// the system's residues and beads tell (Martini 3 when they tell nothing); a massless bead is a virtual site. Atoms of
// one label with a different element, mass, size or particle kind are of different types.
void md_util_system_infer_particle_types(md_system_t* sys, const str_t atom_labels[]);

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
bool md_util_system_infer(struct md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags);


// Computes an array of distances between two sets of coordinates in a periodic domain (cell)
// out_dist:  Output array of distances, must have length of (num_a * num_b)
// coord_a:   Array of coordinates (a)
// num_a:     Length of coord_a
// coord_b:   Array of coordinates (b)
// num_b:     Length of coord_b
// cell:      Periodic boundary cell
void md_util_distance_array(float* out_dist_arr, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

// The minimum (maximum) distance between each group of points of a and the points of b, under the minimum image
// convention of cell (NULL, or a cell without periodic axes: plain euclidean). The minimum image is exact in triclinic
// cells too, also where distances reach past half the cell (see MINIMUM AND MAXIMUM DISTANCE BETWEEN SETS in md_util.c).
// Group g is coord_a[a_offsets[g] .. a_offsets[g + 1]), a_offsets holds num_groups + 1 entries.
// out_dist:   The distance per group, 0 for an empty group or an empty b
// out_idx_a:  (optional) Per group, the index into coord_a of the nearest (farthest) pair, -1 where there is none
// out_idx_b:  (optional) Per group, the index into coord_b of the pair, -1 where there is none
// Each group starts from the pair of the one before, so groups which follow each other in space are cheaper.
void md_util_min_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);
void md_util_max_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

// The same for a single group. The minimum is FLT_MAX and the maximum 0 when a or b is empty, and the indices are then
// left unwritten (for the maximum also when it is 0).
float md_util_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);
float md_util_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

// Every pair whose distance is within tol of the minimum (maximum), in one pass with the same culling: for a caller
// which computes the distances another way, where they may differ from these by some ulps, and wants the extreme by
// its own. Up to cap pairs are written, as indices into coord_a and coord_b, in no particular order. Returns how many
// there are, or more than cap when they do not all fit (the ones written are then only some of them); 0 when a or b
// is empty.
size_t md_util_min_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);
size_t md_util_max_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

void md_util_min_image_vec3(vec3_t* in_out_dx, size_t count, const md_unitcell_t* cell);
void md_util_min_image_vec4(vec4_t* in_out_dx, size_t count, const md_unitcell_t* cell);

// Applies periodic boundary conditions to coordinates
bool md_util_pbc(vec3_t* in_out_xyz, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
bool md_util_pbc_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell);

// Applies periodic boundary conditions to all coordinates in a systems state (convenience function)
bool md_util_system_pbc(md_system_state_t* state);

// Unwraps a single molecule by walking the parent hierarchy md_util_system_infer_molecules
// already computed. Linear, allocation free, and does not consult bond connectivity.
void md_util_unwrap_molecule(md_system_state_t* state, const md_molecule_t* mol);

// Unwraps all molecules in a system
void md_util_unwrap_system(md_system_state_t* state, const md_system_t* sys);

// Batch deperiodize a set of coordinates (vec4) with respect to a given reference
bool md_util_deperiodize_vec4(vec4_t* xyzw, size_t count, vec3_t ref_xyz, const md_unitcell_t* cell);

// Converts absolute coordinates (in unitcell) to relative coordinates (with respect to given reference point)
void md_util_convert_to_relative_coordinates_vec4(vec4_t* in_out_rel_xyzw, vec3_t ref_xyz, size_t count, const md_unitcell_t* cell);

// PERIODIC IMAGE SELECTION
//
// These let the estimator choose each point's periodic image, using its own objective as the
// criterion. Prefer them over unwrapping followed by estimation: unwrap has to commit to an image
// assignment before it knows what is being measured, and it commits using topology, which is
// unrelated to the measurement. These do not need topology at all, which is also why they work for
// sparse selections and for points which are not atoms.
//
// @NOTE: md_util_unwrap_* remains the right tool for making a whole connected structure whole - for
// rendering or export - where there is no estimator and no reference to align against, and where the
// structure may be more extended than half a cell. That is the only case it is for.

// Places a set of points in mutually consistent images and reports its centre.
//
// Seeded from the circular mean (which requires no image assignment, so it cannot be skewed by the
// images the input arrived in) and then alternated to a fixed point: place every point in the image
// nearest the centre, recompute the centre, repeat. Settles in one or two passes.
//
// The set is left in the image it arrived in: the points are made mutually consistent and the
// cluster as a whole is NOT moved into the reference cell, so the reported centre can be combined
// with coordinates that were not passed in. Call md_util_pbc afterwards if the cell is wanted.
// A set straddling a boundary has no image to preserve and may come back on either side.
//
// Correct while the set spans less than half a cell. Beyond that no per point image choice can
// represent it and the structure has to be unwrapped topologically instead.
bool md_util_deperiodize_self_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell, vec3_t* out_com);

// Optimal rigid rotation between a reference set and a target set under periodic boundaries.
//
// Solves for the rotation, the target centre AND each target point's periodic image together:
//     min over R, c, n_k  of  sum_k w_k | R (q_k + A n_k - c) - (p_k - ref_com) |^2
// Given R and c the best n_k is the image nearest the predicted position, per point and in closed
// form; given the n_k the best R and c is ordinary Kabsch. Alternating the two decreases the
// objective monotonically and terminates when no point changes image.
//
// out_rot maps the TARGET frame onto the REFERENCE frame: R * (q - out_com) ~= p - ref_com.
// out_trg_xyzw optionally receives the placed target points (weights preserved); pass NULL if only
// the transform is wanted, and it may alias trg_xyzw.
//
// out_com and the placed points are expressed in the image the TARGET arrived in, not folded into
// the reference cell - see md_util_deperiodize_self_vec4. ref_com is used only as the reference
// set's own centre, so its image does not matter.
//
// Returns the largest per point residual, which is the margin to an ambiguous image choice. Small
// against half a cell means the assignment is nowhere near flipping; approaching it means this
// frame's alignment should not be trusted.
float md_util_optimal_rotation_pbc_vec4_iter(mat3_t* out_rot, vec3_t* out_com, vec4_t* out_trg_xyzw,
                                             const vec4_t* ref_xyzw, vec3_t ref_com,
                                             const vec4_t* trg_xyzw, size_t count,
                                             const md_unitcell_t* cell, int max_iter, float tol);

mat3_t md_util_optimal_rotation_rel_vec4(const vec4_t* ref_rel_xyzw, const vec4_t* trg_rel_xyzw, size_t count);

// Computes the minimum axis aligned bounding box for a set of points with a given radius
// Indices are optional and are used to select a subset of points, the count dictates the number of elements to process
void md_util_aabb_compute     (float out_ext_min[3], float out_ext_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count);
void md_util_aabb_compute_vec4(float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count);

// Computes an object oriented bounding box based on the PCA of the provided points (with optional radius)
void md_util_oobb_compute     (float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
void md_util_oobb_compute_vec4(float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);

// Computes the center of mass for a set of points with a given weight
// xyz:         Packed coordinates
// w:           Array of weights (optional): set as NULL to use equal weights
// indices:     Array of indices (optional): indices into the arrays (xyz,w)
// count:       Length of all arrays
// unit_cell:   The unit_cell of the system [Optional]
vec3_t md_util_com_compute(const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
vec3_t md_util_com_compute_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);

// Computes the similarity between two sets of points with given weights.
// One of the sets is rotated and translated to match the other set in an optimal fashion before the similarity is computed.
// The rmsd is the root mean squared deviation between the two sets of aligned vectors.
// xyz:     Packed coordinate arrays [2]
// com:     Center of mass [2] (xyz0), (xyz1)
// w:       Array of weights (optional): set as NULL to use equal weights
// count:   Length of all arrays (xyz0, xyz1, w)
double md_util_rmsd_compute(const vec3_t* const in_xyz[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);
double md_util_rmsd_compute_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);

// Computes linear shape descriptor weights (linear, planar, isotropic) from a covariance matrix
vec3_t md_util_shape_weights(const mat3_t* covariance_matrix);

// Perform linear interpolation of supplied coordinates
// out_xyz:     Destination, packed
// in_xyz:      Sources [2], packed
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t: interpolation factor (0..1)
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
bool md_util_interpolate_linear(vec3_t* out_xyz, const vec3_t* const in_xyz[2], size_t count, const md_unitcell_t* cell, float t);

// Perform cubic interpolation of supplied coordinates
// out_xyz:     Destination, packed
// in_xyz:      Sources [4], packed
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t:           Interpolation factor (0..1)
// s:           Scaling factor (0..1), 0 is jerky, 0.5 is catmul rom, 1.0 is silky smooth
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
bool md_util_interpolate_cubic_spline(vec3_t* out_xyz, const vec3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, float t, float s);

// The weights of four consecutive frames at t (0..1) for a cubic scheme over them, on the interval between the
// two middle ones: out = w[0] in[0] + w[1] in[1] + w[2] in[2] + w[3] in[3]. They sum to one.
// Cardinal spline with tangent scale s (0.5 is Catmull-Rom): C1, and passes through in[1] at t = 0 and in[2] at t = 1.
void md_util_cubic_spline_weights(float out_w[4], float t, float s);
// Uniform cubic B-spline with the four frames as its control points: C2, and damps frame to frame jitter, but it
// approximates rather than interpolates: at t = 0 it is (in[0] + 4 in[1] + in[2]) / 6, not in[1].
void md_util_cubic_bspline_weights(float out_w[4], float t);

// Weighted sum of four consecutive frames of coordinates, out = w[0] in[0] + w[1] in[1] + w[2] in[2] + w[3] in[3]:
// any cubic scheme over them at a given t (see the weights above). In a periodic cell in[0] and in[2] are taken at
// the minimum image of in[1], and in[3] at that of in[2].
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
bool md_util_interpolate_cubic_weighted(vec3_t* out_xyz, const vec3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, const float w[4]);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial(uint32_t* source_indices, const vec3_t* xyz, size_t count);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial_xyz(uint32_t* source_indices, const float* xyz, size_t stride_in_bytes, size_t count);

// Radix sorts of unsigned integer keys, in place (sort) or as the order of the keys (argsort): out_indices[0] is the index
// of the smallest key. The argsorts are stable: indices of equal keys come in increasing order. They take temporary memory
// from the temp arena of the calling thread, for the keys and indices between passes.
void md_util_sort_radix_u32(uint32_t* data, size_t count);
void md_util_sort_radix_u64(uint64_t* data, size_t count);
void md_util_argsort_radix_u32(uint32_t* out_indices, const uint32_t* keys, size_t count);
void md_util_argsort_radix_u64(uint32_t* out_indices, const uint64_t* keys, size_t count);

#ifdef __cplusplus
}
#endif

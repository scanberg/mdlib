#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

#include <stdint.h>

#if DEBUG
#include <core/md_log.h>
#endif

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

enum {
    MD_RAMACHANDRAN_TYPE_UNKNOWN = 0,
    MD_RAMACHANDRAN_TYPE_GENERAL,
    MD_RAMACHANDRAN_TYPE_GLYCINE,
    MD_RAMACHANDRAN_TYPE_PROLINE,
    MD_RAMACHANDRAN_TYPE_PREPROL,
};

// CLASSIFICATION
//
// What a part of a system is, is stored at the level of the hierarchy it describes, and only there:
//
//   atom type  md_atom_type_flags_t  what the particle is: an atom, a coarse grained bead or a virtual site
//   atom       md_atom_flags_t       the atom's role in its component (backbone, side chain, ...) and its chemistry
//   component  md_component_flags_t  what the component is (amino acid, nucleotide, water, ion), and where in its chain
//   entity     md_entity_flags_t     what the molecule is (peptide, DNA, water, ...)
//
// Components of the same name and kind share a component type (md_component_type_data_t), as atoms share an atom
// type. What an atom type is a form of is its element, what a component type is a form of is its standard monomer
// (md_monomer_t): HIE, HSD and MSE are forms of HIS, HIS and MET. Both are where the defaults of the type come from.
//
// Instances have no flags of their own: an instance is one occurrence of its entity (md_system_instance_entity_kind).
// Nothing is copied from one level to another. A property of an atom's type or of its component is read through the
// type or the component (md_atom_particle_kind, md_system_atom_component_kind), so there is one copy and it cannot go
// stale when the type or the component is reclassified.
//
// Options which exclude each other are a VALUE, packed into a few bits and read and written through the accessors
// below (md_component_flags_kind, md_component_flags_set_kind, ...). Independent properties are single bits. No
// meaning is carried by a combination of bits.

// ### ATOM TYPE ###

// What a particle is. A property of its type: every particle of a type is the same kind of particle.
typedef enum md_particle_kind_t {
    MD_PARTICLE_ATOM         = 0,   // An atom. Its element is the type's atomic number, 0 when it is not known.
    MD_PARTICLE_BEAD         = 1,   // A coarse grained bead which stands for a group of atoms. A system with beads is
                                    // coarse grained (md_system_is_coarse_grained): bonds and hydrogen bonds are not
                                    // inferred from the geometry, as atomic radii do not apply.
    MD_PARTICLE_VIRTUAL_SITE = 2,   // A massless interaction site without an element: the M site of 4 site water models
                                    // (TIP4P, OPC), the lone pairs of 5 site models (TIP5P). It belongs to its molecule
                                    // but forms no bonds.
} md_particle_kind_t;

typedef enum {
    MD_ATOM_TYPE_FLAG_NONE                  = 0,
    MD_ATOM_TYPE_FLAG_PARTICLE_KIND_MASK    = 0x3,      // md_particle_kind_t, see md_atom_type_flags_particle_kind
} md_atom_type_flags_t;

ENUM_FLAGS(md_atom_type_flags_t)

// ### ATOM ###

// Hybridization of a heavy atom, see md_chem.h. SP2 includes the N and O whose lone pair takes part in a pi system
// (amide and aniline N, pyrrole N, furan O).
typedef enum md_hybridization_t {
    MD_HYBRIDIZATION_UNKNOWN = 0,   // Not perceived, or not applicable (hydrogens, metals)
    MD_HYBRIDIZATION_SP      = 1,
    MD_HYBRIDIZATION_SP2     = 2,
    MD_HYBRIDIZATION_SP3     = 3,
} md_hybridization_t;

typedef enum {
    MD_ATOM_FLAG_NONE                   = 0,

    // The atom's role in its component, set when the component is classified (md_util_system_infer_comp_flags) and
    // for the beads of the predefined coarse grained types. Only the components whose backbone was verified
    // (MD_COMPONENT_FLAG_RESOLVED) have their atoms given a role.
    MD_ATOM_FLAG_BACKBONE               = 0x1,      // N, CA, C of an amino acid; P, O5', C5', C4', C3', O3' of a nucleotide
    MD_ATOM_FLAG_SIDE_CHAIN             = 0x2,      // Amino acid side chain, from CB outwards
    MD_ATOM_FLAG_NUCLEOSIDE             = 0x4,      // Sugar and base of a nucleotide
    MD_ATOM_FLAG_NUCLEOBASE             = 0x8,      // Base of a nucleotide
    MD_ATOM_FLAG_TERMINAL_BEG           = 0x10,     // In the terminal group which begins a chain (N-terminus, 5' end)
    MD_ATOM_FLAG_TERMINAL_END           = 0x20,     // In the terminal group which ends a chain (C-terminus, 3' end)
    // In a tail of a lipid (md_util_system_infer_lipids): an acyl, alkyl or sphingoid chain, or the ring skeleton of a
    // sterol, with its branches and hydrogens. The other atoms of a lipid are its head(s).
    MD_ATOM_FLAG_LIPID_TAIL             = 0x40,

    // Chemistry, see md_chem.h (md_chem_perceive)
    MD_ATOM_FLAG_HYBRIDIZATION_MASK     = 0x300,    // md_hybridization_t, see md_atom_flags_hybridization
    MD_ATOM_FLAG_AROMATIC               = 0x400,    // In an aromatic ring

    // Treated quantum mechanically: an atom of the QM calculation the system was loaded from, set by its loader.
    // A system from a QM calculation alone is QM throughout; one with an embedding or a QM/MM partition has its
    // QM region flagged, and the atoms without the flag are its environment. Nothing infers it: a system no
    // loader said it of has no QM region. Selected in scripts by 'qm' and 'environment'.
    MD_ATOM_FLAG_QM                     = 0x800,

    // Hydrogen bond roles are not flags: they depend on the hydrogen bond model and its options, and a donor is a
    // D-H pair rather than an atom. See md_hbond_perceive_roles and md_hbond_query_t.
} md_atom_flags_t;

ENUM_FLAGS(md_atom_flags_t)

// ### COMPONENT ###

// What a component (residue) is, see md_util_system_infer_comp_flags.
typedef enum md_component_kind_t {
    MD_COMPONENT_KIND_OTHER         = 0,    // Unclassified, or none of the below: a ligand, a sugar, ...
    MD_COMPONENT_KIND_AMINO_ACID    = 1,
    MD_COMPONENT_KIND_NUCLEOTIDE    = 2,
    MD_COMPONENT_KIND_WATER         = 3,
    MD_COMPONENT_KIND_ION           = 4,    // A monatomic ion
    MD_COMPONENT_KIND_LIPID         = 5,    // A lipid, sterols included: by name (md_util_resname_lipid), or by its structure
                                            // (a polar head with hydrocarbon tails, see md_util_system_infer_lipids)
} md_component_kind_t;

typedef enum {
    MD_COMPONENT_FLAG_NONE              = 0,
    MD_COMPONENT_FLAG_KIND_MASK         = 0x7,      // md_component_kind_t, see md_component_flags_kind

    // The backbone of the amino acid or nucleotide was found by atom names and verified by its bonds, and its atoms
    // were given their roles (MD_ATOM_FLAG_BACKBONE, ...). Without it the component is one by name only: an unusual
    // naming scheme, missing atoms, or a coarse grained model.
    MD_COMPONENT_FLAG_RESOLVED          = 0x10,
    MD_COMPONENT_FLAG_TERMINAL_BEG      = 0x20,     // Begins its chain: a free amino group (N-terminus), a free 5' end
    MD_COMPONENT_FLAG_TERMINAL_END      = 0x40,     // Ends its chain: a free carboxyl group (C-terminus), a free 3' end

    // The chain is broken next to it: its backbone is not linked to that of the residue before (after) it in its chain,
    // which carries BREAK_AFTER (BREAK_BEFORE) in turn. Residues missing in a crystal structure leave such a break. See
    // CHAIN BREAKS in md_util.h.
    MD_COMPONENT_FLAG_BREAK_BEFORE      = 0x80,
    MD_COMPONENT_FLAG_BREAK_AFTER       = 0x100,
} md_component_flags_t;

ENUM_FLAGS(md_component_flags_t)

// ### MONOMER ###

// The standard monomers of proteins and nucleic acids, the predefined table component types are matched against
// (md_monomer_from_comp_name), as atom types are matched against the elements. A monomer has its id in the Chemical
// Component Dictionary of the PDB, a one letter code, a kind and a default color (md_monomer_*). The values are not
// stable across versions: store a monomer by its id (md_monomer_id, md_monomer_from_id), not by its value.
typedef enum md_monomer_t {
    MD_MONOMER_NONE = 0,    // Not a form of a standard monomer: a ligand, water, an ion, a lipid, ...

    // Amino acids
    MD_MONOMER_ALA, MD_MONOMER_ARG, MD_MONOMER_ASN, MD_MONOMER_ASP, MD_MONOMER_CYS,
    MD_MONOMER_GLN, MD_MONOMER_GLU, MD_MONOMER_GLY, MD_MONOMER_HIS, MD_MONOMER_ILE,
    MD_MONOMER_LEU, MD_MONOMER_LYS, MD_MONOMER_MET, MD_MONOMER_PHE, MD_MONOMER_PRO,
    MD_MONOMER_SER, MD_MONOMER_THR, MD_MONOMER_TRP, MD_MONOMER_TYR, MD_MONOMER_VAL,
    MD_MONOMER_SEC,         // Selenocysteine (U)
    MD_MONOMER_PYL,         // Pyrrolysine (O)
    MD_MONOMER_UNK,         // An amino acid which is none of the above (a modified residue of no known parent), or a
                            // cap of a chain (ACE, NME)

    // Deoxyribonucleotides
    MD_MONOMER_DA, MD_MONOMER_DC, MD_MONOMER_DG, MD_MONOMER_DT, MD_MONOMER_DU, MD_MONOMER_DI,
    // Ribonucleotides
    MD_MONOMER_A, MD_MONOMER_C, MD_MONOMER_G, MD_MONOMER_U, MD_MONOMER_I,
    MD_MONOMER_N,           // A nucleotide which is none of the above

    MD_MONOMER_COUNT
} md_monomer_t;

// ### ENTITY ###

// What a molecule is. An entity is a molecule type (or a polymer sequence), each of its instances one molecule or one
// polymer chain of that type.
typedef enum md_entity_kind_t {
    MD_ENTITY_KIND_UNKNOWN      = 0,    // Not classified yet: a topology names its molecule types but does not say what
                                        // they are. md_util_system_infer classifies these from their components.
    MD_ENTITY_KIND_NON_POLYMER  = 1,    // A small molecule: a ligand, an ion, a lipid, a single residue
    MD_ENTITY_KIND_WATER        = 2,
    MD_ENTITY_KIND_BRANCHED     = 3,    // A branched oligosaccharide (mmCIF 'branched')
    // Polymers from here on, see md_entity_kind_is_polymer
    MD_ENTITY_KIND_PEPTIDE      = 4,    // A polypeptide chain (protein)
    MD_ENTITY_KIND_DNA          = 5,
    MD_ENTITY_KIND_RNA          = 6,
    MD_ENTITY_KIND_NUCLEIC      = 7,    // A nucleic acid which is both (a DNA/RNA hybrid), or which could not be told apart
    MD_ENTITY_KIND_POLYMER      = 8,    // Any other polymer
} md_entity_kind_t;

typedef enum {
    MD_ENTITY_FLAG_NONE         = 0,
    MD_ENTITY_FLAG_KIND_MASK    = 0xF,      // md_entity_kind_t, see md_entity_flags_kind

    // Not given by the source: the entity was inferred together with its instances (md_util_system_infer_entity_and_instance)
    MD_ENTITY_FLAG_INFERRED     = 0x10,
} md_entity_flags_t;

ENUM_FLAGS(md_entity_flags_t)

// What a tail of a lipid is (md_lipid_data_t)
typedef enum md_lipid_tail_kind_t {
    MD_LIPID_TAIL_CHAIN     = 0,    // A hydrocarbon chain: acyl, alkyl, alkenyl or sphingoid, branched or not (phytanyl).
                                    // It hangs off a head, or spans two (the chains of a bolalipid).
    MD_LIPID_TAIL_SKELETON  = 1,    // The fused ring system of a sterol, with its side chain
} md_lipid_tail_kind_t;

static inline bool md_entity_kind_is_polymer(md_entity_kind_t kind) {
    return kind >= MD_ENTITY_KIND_PEPTIDE;
}

static inline bool md_entity_kind_is_nucleic_acid(md_entity_kind_t kind) {
    return kind == MD_ENTITY_KIND_DNA || kind == MD_ENTITY_KIND_RNA || kind == MD_ENTITY_KIND_NUCLEIC;
}

// ### ACCESSORS for the packed values ###

static inline md_particle_kind_t md_atom_type_flags_particle_kind(md_atom_type_flags_t flags) {
    return (md_particle_kind_t)((unsigned int)flags & (unsigned int)MD_ATOM_TYPE_FLAG_PARTICLE_KIND_MASK);
}

static inline md_atom_type_flags_t md_atom_type_flags_set_particle_kind(md_atom_type_flags_t flags, md_particle_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_ATOM_TYPE_FLAG_PARTICLE_KIND_MASK;
    return (md_atom_type_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

#define MD_ATOM_HYBRIDIZATION_SHIFT 8

static inline md_hybridization_t md_atom_flags_hybridization(md_atom_flags_t flags) {
    return (md_hybridization_t)(((unsigned int)flags & (unsigned int)MD_ATOM_FLAG_HYBRIDIZATION_MASK) >> MD_ATOM_HYBRIDIZATION_SHIFT);
}

static inline md_atom_flags_t md_atom_flags_set_hybridization(md_atom_flags_t flags, md_hybridization_t hyb) {
    const unsigned int mask = (unsigned int)MD_ATOM_FLAG_HYBRIDIZATION_MASK;
    return (md_atom_flags_t)(((unsigned int)flags & ~mask) | (((unsigned int)hyb << MD_ATOM_HYBRIDIZATION_SHIFT) & mask));
}

static inline md_component_kind_t md_component_flags_kind(md_component_flags_t flags) {
    return (md_component_kind_t)((unsigned int)flags & (unsigned int)MD_COMPONENT_FLAG_KIND_MASK);
}

static inline md_component_flags_t md_component_flags_set_kind(md_component_flags_t flags, md_component_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_COMPONENT_FLAG_KIND_MASK;
    return (md_component_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

static inline md_entity_kind_t md_entity_flags_kind(md_entity_flags_t flags) {
    return (md_entity_kind_t)((unsigned int)flags & (unsigned int)MD_ENTITY_FLAG_KIND_MASK);
}

static inline md_entity_flags_t md_entity_flags_set_kind(md_entity_flags_t flags, md_entity_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_ENTITY_FLAG_KIND_MASK;
    return (md_entity_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}


// Bond flags, on three axes: the low byte describes the bond chemically, the bits above it where the bond and its
// chemistry came from, and the character bits what a non-covalent bond is. Only the low byte is part of a bond's
// identity when comparing structures (non-covalent bonds take no part in that: they are not in the bond graph). Each
// field is either one bit or one value, so every combination describes a case (two of them aside: AROMATIC with
// DELOCALIZED, and two kinds at once, which never occur together).
//
// CHEMISTRY
//   ORDER        A value in bits 0-2, not a set of bits: read it with md_bond_order and write it with
//                md_bond_flags_set_order. 0 is unknown (not perceived, and not given by the file), which is not the
//                same as single. Always 0 for a non-covalent bond.
//   AROMATIC     In an aromatic ring.
//   DELOCALIZED  Resonance outside of rings: carboxylate, guanidinium, nitro, phosphate, ... AROMATIC and DELOCALIZED
//                come on top of an order which holds one Kekule structure: a benzene bond is AROMATIC with order 1 or
//                2, the two C-O of a carboxylate are DELOCALIZED with orders 2 and 1.
//   KIND         COORDINATE, ELASTIC, NONCOVALENT or none of them (covalent), at most one:
//     covalent     None of the bits. Metallic between two metals: whether a bond involves a metal is a question for
//                  the elements of its atoms.
//     COORDINATE   A metal-ligand bond, which counts towards the valence of neither atom: every bond between a metal
//                  and a non-metal unless its kind was given (md_util_system_infer_coordination), and any other bond a
//                  file says is one. Whether the ligand is then an anion (the thiolate of a zinc finger cysteine) is
//                  for the perception of its atom to say (md_chem_perceive), not for the bond.
//     ELASTIC      Not a bond of the molecule but a spring of its force field which holds two particles at a distance:
//                  the elastic network of a coarse grained protein or nucleic acid (the rubber bands of martinize,
//                  ElNeDyn), the restraints of its secondary structure, a restraint between two atoms. GROMACS writes
//                  them as bonds which make no connection (bond types 6, 9 and 10, constraints of type 2), or as plain
//                  bonds (md_util_system_infer_elastic_bonds). They are kept, as information, in the first range of
//                  md_bond_data_t, and left out of everything which follows the bonds of a molecule: the connectivity
//                  (md_bond_conn_data_t, so md_bond_iter and md_bond_find), the structures, rings, chemistry and
//                  matching. Nothing else of the CHEMISTRY applies to them, and the representations do not draw them
//                  unless asked to.
//     NONCOVALENT  A bond which holds nothing together: a hydrogen bond, a salt bridge, an interaction a user marks.
//                  Never part of the bond graph (md_bond_data_t keeps these bonds in a range of their own), so never
//                  of an order, aromatic or delocalized. What it is, is said by the character bits.
//   COVALENT and COORDINATE bonds make up the bond graph, from which the connectivity, the rings, the structures
//   and the perceived chemistry are derived.
//
// CHARACTER (non-covalent bonds only, any combination of them; none is an interaction of no particular character)
//   HYDROGEN     A hydrogen bond, D-H...A. Stored with the hydrogen (or the donor without one) as the first atom.
//   IONIC        Between ions of opposite charge: a salt bridge. A charge-assisted hydrogen bond is both.
//   SIGMA_HOLE   From the sigma-hole of a halogen, chalcogen or pnictogen (the first atom) to an acceptor: a halogen
//                bond is one whose first atom is a halogen.
//
// PROVENANCE
//   ORIGIN           A value in bits 8-9 (md_bond_origin_t): read it with md_bond_origin and write it with
//                    md_bond_flags_set_origin.
//   ORDER_PERCEIVED  The order, AROMATIC and DELOCALIZED were perceived (md_chem_perceive), not given. A known order
//                    without it was given, by the file or by hand, and stays as it is when the chemistry is perceived
//                    again; one with it is perceived anew.
//   KIND_GIVEN       The kind (and the character) was given, by the file or by hand, and stays as it is when bonds are
//                    perceived again: a bond between a metal and a non-metal a user says is covalent stays covalent. The
//                    other way around from ORDER_PERCEIVED, as a kind has no 'unknown' to tell a given one by: a bond
//                    without it has the kind perception gives it.
typedef enum {
    MD_BOND_FLAG_NONE            = 0,
    MD_BOND_FLAG_ORDER_MASK      = 0x7,     // See md_bond_order
    MD_BOND_FLAG_AROMATIC        = 0x8,     // In an aromatic ring
    MD_BOND_FLAG_DELOCALIZED     = 0x10,    // Resonance outside of rings
    MD_BOND_FLAG_COORDINATE      = 0x20,    // Metal-ligand
    MD_BOND_FLAG_ELASTIC         = 0x40,    // A spring of the force field, not a bond of the molecule
    MD_BOND_FLAG_NONCOVALENT     = 0x80,    // Not in the bond graph, see the character bits
    MD_BOND_FLAG_KIND_MASK       = 0xE0,    // COORDINATE | ELASTIC | NONCOVALENT, none of them is covalent

    MD_BOND_FLAG_ORIGIN_MASK     = 0x300,   // See md_bond_origin
    MD_BOND_FLAG_ORDER_PERCEIVED = 0x400,   // The order (and AROMATIC, DELOCALIZED) was perceived, not given
    MD_BOND_FLAG_KIND_GIVEN      = 0x800,   // The kind (and the character) was given, not perceived

    MD_BOND_FLAG_HYDROGEN        = 0x1000,  // Non-covalent: a hydrogen bond
    MD_BOND_FLAG_IONIC           = 0x2000,  // Non-covalent: between ions of opposite charge
    MD_BOND_FLAG_SIGMA_HOLE      = 0x4000,  // Non-covalent: from the sigma-hole of its first atom
    MD_BOND_FLAG_CHARACTER_MASK  = 0x7000,
} md_bond_flags_t;

ENUM_FLAGS(md_bond_flags_t)

enum {
    MD_BOND_ORDER_UNKNOWN   = 0,
    MD_BOND_ORDER_SINGLE    = 1,
    MD_BOND_ORDER_DOUBLE    = 2,
    MD_BOND_ORDER_TRIPLE    = 3,
    MD_BOND_ORDER_QUADRUPLE = 4,
};

// The sizes of the rings of md_system_t::ring (md_util_system_infer_rings). Larger cycles are loops, not rings.
#define MD_RING_MIN_SIZE 3
#define MD_RING_MAX_SIZE 7

// Ring flags (md_system_t::ring_flags), what a ring of the system is. md_util_system_infer_ring_flags derives them
// from the rings, the bonds and their flags; md_util_system_infer_rings and md_chem_perceive call it.
//   AROMATIC  Every bond of the ring is aromatic (md_chem_perceive, or the bonds a file gives): an aromatic ring, alone
//             or fused with a neighbor.
//   HETERO    The ring has an atom other than carbon.
//   FUSED     The ring shares a bond with another ring.
//   SUGAR     The ring of a monosaccharide, a furanose or a pyranose: 5 or 6 atoms, one of them O and the rest C, no
//             bond of a known order above single or aromatic in or out of the ring, and at least two of the C carrying
//             an O or N outside of the ring (a hydroxyl, a glycosidic link, the N of a nucleoside), one of them next to
//             the ring O (the anomeric C). Deoxy sugars (the 2-deoxyribose of DNA) count; tetrahydrofuran does not.
//             The atoms of a sugar ring are in canonical order in md_system_t::ring: the ring O, the anomeric C, and on
//             around the ring away from the O (O5 C1 C2 C3 C4 C5 of a pyranose, O4' C1' C2' C3' C4' of a furanose).
//   FUSED_AROMATIC
//             The ring shares a bond with an aromatic ring. With AROMATIC: a ring of a fused aromatic system, whose
//             sextets are not assigned to particular rings (naphthalene, indole, a purine, azulene). An aromatic ring
//             without it holds a sextet of its own, alone or fused to rings which are not aromatic (benzene, the
//             benzene ring of tetralin or of estradiol, those of biphenylene).
typedef enum md_ring_flags_t {
    MD_RING_FLAG_NONE           = 0,
    MD_RING_FLAG_AROMATIC       = 1u << 0,
    MD_RING_FLAG_HETERO         = 1u << 1,
    MD_RING_FLAG_FUSED          = 1u << 2,
    MD_RING_FLAG_SUGAR          = 1u << 3,
    MD_RING_FLAG_FUSED_AROMATIC = 1u << 4,
} md_ring_flags_t;

ENUM_FLAGS(md_ring_flags_t)

// Where a bond came from. FILE is 0, so a bond added without saying where it came from counts as given.
typedef enum md_bond_origin_t {
    MD_BOND_ORIGIN_FILE     = 0,    // Given by the file, for some of its atoms (CONECT records, a component's bonds)
    MD_BOND_ORIGIN_TOPOLOGY = 1,    // A force field topology: all of the bonds of the atoms it covers
    MD_BOND_ORIGIN_USER     = 2,    // Added by hand
    MD_BOND_ORIGIN_INFERRED = 3,    // Inferred from the geometry, and replaced when inferred again
} md_bond_origin_t;

#define MD_BOND_ORDER_SHIFT  0
#define MD_BOND_ORIGIN_SHIFT 8

static inline int md_bond_order(md_bond_flags_t flags) {
    return (int)((flags & MD_BOND_FLAG_ORDER_MASK) >> MD_BOND_ORDER_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_order(md_bond_flags_t flags, int order) {
    const unsigned int value = (order < 0 ? 0u : (order > 4 ? 4u : (unsigned int)order)) << MD_BOND_ORDER_SHIFT;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORDER_MASK) | value);
}

static inline md_bond_origin_t md_bond_origin(md_bond_flags_t flags) {
    return (md_bond_origin_t)(((unsigned int)flags & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK) >> MD_BOND_ORIGIN_SHIFT);
}

static inline md_bond_flags_t md_bond_flags_set_origin(md_bond_flags_t flags, md_bond_origin_t origin) {
    const unsigned int value = ((unsigned int)origin << MD_BOND_ORIGIN_SHIFT) & (unsigned int)MD_BOND_FLAG_ORIGIN_MASK;
    return (md_bond_flags_t)(((unsigned int)flags & ~(unsigned int)MD_BOND_FLAG_ORIGIN_MASK) | value);
}

// Atomic number constants for all elements (Z values)
enum {
    MD_Z_X  = 0,   // Unknown
    MD_Z_H  = 1,   // Hydrogen
    MD_Z_He = 2,   // Helium
    MD_Z_Li = 3,   // Lithium
    MD_Z_Be = 4,   // Beryllium
    MD_Z_B  = 5,   // Boron
    MD_Z_C  = 6,   // Carbon
    MD_Z_N  = 7,   // Nitrogen
    MD_Z_O  = 8,   // Oxygen
    MD_Z_F  = 9,   // Fluorine
    MD_Z_Ne = 10,  // Neon
    MD_Z_Na = 11,  // Sodium
    MD_Z_Mg = 12,  // Magnesium
    MD_Z_Al = 13,  // Aluminium
    MD_Z_Si = 14,  // Silicon
    MD_Z_P  = 15,  // Phosphorus
    MD_Z_S  = 16,  // Sulfur
    MD_Z_Cl = 17,  // Chlorine
    MD_Z_Ar = 18,  // Argon
    MD_Z_K  = 19,  // Potassium
    MD_Z_Ca = 20,  // Calcium
    MD_Z_Sc = 21,  // Scandium
    MD_Z_Ti = 22,  // Titanium
    MD_Z_V  = 23,  // Vanadium
    MD_Z_Cr = 24,  // Chromium
    MD_Z_Mn = 25,  // Manganese
    MD_Z_Fe = 26,  // Iron
    MD_Z_Co = 27,  // Cobalt
    MD_Z_Ni = 28,  // Nickel
    MD_Z_Cu = 29,  // Copper
    MD_Z_Zn = 30,  // Zinc
    MD_Z_Ga = 31,  // Gallium
    MD_Z_Ge = 32,  // Germanium
    MD_Z_As = 33,  // Arsenic
    MD_Z_Se = 34,  // Selenium
    MD_Z_Br = 35,  // Bromine
    MD_Z_Kr = 36,  // Krypton
    MD_Z_Rb = 37,  // Rubidium
    MD_Z_Sr = 38,  // Strontium
    MD_Z_Y  = 39,  // Yttrium
    MD_Z_Zr = 40,  // Zirconium
    MD_Z_Nb = 41,  // Niobium
    MD_Z_Mo = 42,  // Molybdenum
    MD_Z_Tc = 43,  // Technetium
    MD_Z_Ru = 44,  // Ruthenium
    MD_Z_Rh = 45,  // Rhodium
    MD_Z_Pd = 46,  // Palladium
    MD_Z_Ag = 47,  // Silver
    MD_Z_Cd = 48,  // Cadmium
    MD_Z_In = 49,  // Indium
    MD_Z_Sn = 50,  // Tin
    MD_Z_Sb = 51,  // Antimony
    MD_Z_Te = 52,  // Tellurium
    MD_Z_I  = 53,  // Iodine
    MD_Z_Xe = 54,  // Xenon
    MD_Z_Cs = 55,  // Caesium
    MD_Z_Ba = 56,  // Barium
    MD_Z_La = 57,  // Lanthanum
    MD_Z_Ce = 58,  // Cerium
    MD_Z_Pr = 59,  // Praseodymium
    MD_Z_Nd = 60,  // Neodymium
    MD_Z_Pm = 61,  // Promethium
    MD_Z_Sm = 62,  // Samarium
    MD_Z_Eu = 63,  // Europium
    MD_Z_Gd = 64,  // Gadolinium
    MD_Z_Tb = 65,  // Terbium
    MD_Z_Dy = 66,  // Dysprosium
    MD_Z_Ho = 67,  // Holmium
    MD_Z_Er = 68,  // Erbium
    MD_Z_Tm = 69,  // Thulium
    MD_Z_Yb = 70,  // Ytterbium
    MD_Z_Lu = 71,  // Lutetium
    MD_Z_Hf = 72,  // Hafnium
    MD_Z_Ta = 73,  // Tantalum
    MD_Z_W  = 74,  // Tungsten
    MD_Z_Re = 75,  // Rhenium
    MD_Z_Os = 76,  // Osmium
    MD_Z_Ir = 77,  // Iridium
    MD_Z_Pt = 78,  // Platinum
    MD_Z_Au = 79,  // Gold
    MD_Z_Hg = 80,  // Mercury
    MD_Z_Tl = 81,  // Thallium
    MD_Z_Pb = 82,  // Lead
    MD_Z_Bi = 83,  // Bismuth
    MD_Z_Po = 84,  // Polonium
    MD_Z_At = 85,  // Astatine
    MD_Z_Rn = 86,  // Radon
    MD_Z_Fr = 87,  // Francium
    MD_Z_Ra = 88,  // Radium
    MD_Z_Ac = 89,  // Actinium
    MD_Z_Th = 90,  // Thorium
    MD_Z_Pa = 91,  // Protactinium
    MD_Z_U  = 92,  // Uranium
    MD_Z_Np = 93,  // Neptunium
    MD_Z_Pu = 94,  // Plutonium
    MD_Z_Am = 95,  // Americium
    MD_Z_Cm = 96,  // Curium
    MD_Z_Bk = 97,  // Berkelium
    MD_Z_Cf = 98,  // Californium
    MD_Z_Es = 99,  // Einsteinium
    MD_Z_Fm = 100, // Fermium
    MD_Z_Md = 101, // Mendelevium
    MD_Z_No = 102, // Nobelium
    MD_Z_Lr = 103, // Lawrencium
    MD_Z_Rf = 104, // Rutherfordium
    MD_Z_Db = 105, // Dubnium
    MD_Z_Sg = 106, // Seaborgium
    MD_Z_Bh = 107, // Bohrium
    MD_Z_Hs = 108, // Hassium
    MD_Z_Mt = 109, // Meitnerium
    MD_Z_Ds = 110, // Darmstadtium
    MD_Z_Rg = 111, // Roentgenium
    MD_Z_Cn = 112, // Copernicium
    MD_Z_Nh = 113, // Nihonium
    MD_Z_Fl = 114, // Flerovium
    MD_Z_Mc = 115, // Moscovium
    MD_Z_Lv = 116, // Livermorium
    MD_Z_Ts = 117, // Tennessine
    MD_Z_Og = 118, // Oganesson
    MD_Z_Count
};

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_component_idx_t;
typedef int32_t     md_instance_idx_t;
typedef int32_t     md_entity_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_sequence_id_t;
typedef int32_t     md_bond_idx_t;
typedef uint16_t    md_atom_type_idx_t;
typedef int32_t     md_component_type_idx_t;
typedef uint8_t     md_atomic_number_t;
typedef uint8_t     md_ramachandran_type_t;
typedef uint8_t     md_order_t;

// For backwards compatibility
typedef uint8_t md_element_t;  // Atomic number (1-118), 0 is unknown

// Open ended range of indices (e.g. range(0,4) -> [0,1,2,3])
typedef struct md_irange_t {
    int32_t beg;
    int32_t end;
} md_irange_t;

typedef struct md_urange_t {
    uint32_t beg;
    uint32_t end;
} md_urange_t;

typedef struct md_atom_pair_t {
    md_atom_idx_t idx[2];
} md_atom_pair_t;

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

// Miniature string buffer with explicit length
// It can store up to 6 characters + null terminator which makes it compatible as a C-string
// This is to store the atom symbol and other short strings common in molecular data
// The motiviation is that it saves an indirection
typedef struct md_label_t {
    char    buf[7];
    uint8_t len;

#ifdef __cplusplus
    constexpr operator str_t() { return {buf, len}; }
    constexpr operator const char*() { return buf; }
#endif
} md_label_t;

// Container structure for ranges of indices.
// It stores a set of indices for every stored entity. (Think of it as an Array of Arrays of integers)
// The indices are packed into a single array and we store explicit offsets to represent ranges within this index array.
// It is used to represent connected structures, rings and atom connectivity etc.

typedef struct md_index_data_t {
    struct md_allocator_i* alloc;
    md_array(uint32_t) offsets;
    md_array(int32_t)  indices;
} md_index_data_t;


// OPERATIONS ON THE TYPES

#ifdef __cplusplus
extern "C" {
#endif

// Names of the kinds for printing ("atom", "amino acid", "peptide", ...), "" for a value out of range
const char* md_particle_kind_name(md_particle_kind_t kind);
const char* md_component_kind_name(md_component_kind_t kind);
const char* md_entity_kind_name(md_entity_kind_t kind);
const char* md_secondary_structure_name(md_secondary_structure_t ss);
const char* md_hybridization_name(md_hybridization_t hyb);
const char* md_bond_origin_name(md_bond_origin_t origin);

// The kinds of bonds by name, as the script selects them and a workspace stores them: a bond is of a kind when
// (flags & mask) == value. 'covalent', 'coordinate', 'elastic' and 'noncovalent' are the kinds proper, 'hydrogen',
// 'ionic' and 'sigma_hole' the characters of a non-covalent bond, so a charge-assisted hydrogen bond is of
// 'noncovalent', 'hydrogen' and 'ionic' at once.
typedef struct md_bond_kind_name_t {
    const char*     name;
    md_bond_flags_t mask;
    md_bond_flags_t value;
} md_bond_kind_name_t;

// The table, in the order above. count receives its length.
const md_bond_kind_name_t* md_bond_kind_names(size_t* count);

// Element property functions
str_t md_atomic_number_name(md_atomic_number_t z);
str_t md_atomic_number_symbol(md_atomic_number_t z);
float md_atomic_number_mass(md_atomic_number_t z);
float md_atomic_number_vdw_radius(md_atomic_number_t z);
float md_atomic_number_covalent_radius(md_atomic_number_t z);
int   md_atomic_number_max_valence(md_atomic_number_t z);
uint32_t md_atomic_number_cpk_color(md_atomic_number_t z);

// Element symbol and name lookup functions
md_atomic_number_t md_atomic_number_from_symbol(str_t sym, bool ignore_case);

// Infer atomic number from other properties, this is a heuristic and may fail. 0 is returned on failure
// Per-atom inference from labels (atom name + residue)
// res_name is the name of the component of the atom [optional]
// res_size is the size of the component of the atom [optional]
md_atomic_number_t md_atomic_number_infer_from_label(str_t atom_name, str_t res_name, size_t res_size);
md_atomic_number_t md_atomic_number_infer_from_mass(float mass);

// Batch form wired to molecule structure
// Returns the number of successfully inferred atomic numbers
// All of the supplied arrays must have at least 'count' elements
size_t md_atomic_number_infer_from_mass_batch(md_atomic_number_t out_z[], const float masses[], size_t count);

// Standard monomers (md_monomer_t), "" / 0 for MD_MONOMER_NONE or a value out of range
str_t    md_monomer_id(md_monomer_t monomer);          // Id in the Chemical Component Dictionary: "ALA", "DA", "A"
str_t    md_monomer_name(md_monomer_t monomer);        // "Alanine", "Deoxyadenosine", "Adenosine"
char     md_monomer_letter(md_monomer_t monomer);      // One letter code: 'A' for ALA, DA and A
md_component_kind_t md_monomer_kind(md_monomer_t monomer);  // MD_COMPONENT_KIND_AMINO_ACID or _NUCLEOTIDE
uint32_t md_monomer_color(md_monomer_t monomer);       // Default color (RGBA, R in the lowest byte, as the CPK colors)

// The monomer of the given id ("HIS"), exactly and ignoring case. MD_MONOMER_NONE if there is none.
md_monomer_t md_monomer_from_id(str_t id);

// The standard monomer a component of the given name and kind is a form of. Only amino acids and nucleotides are forms
// of one: the variants of the force fields (HID, HIE, HSP, HISH, CYX, ASH, LYN, the Martini ARG0, ...), the termini of
// Amber (NALA, CALA, DA5, RA3), the names of CHARMM (ADE, CYT, ...) and the common modified residues (MSE, SEP, PSU,
// 5MC, ...). An amino acid or nucleotide which is not known is MD_MONOMER_UNK or MD_MONOMER_N, any other kind
// MD_MONOMER_NONE, except the caps of a peptide chain (ACE, NME, ...) which are MD_MONOMER_UNK.
md_monomer_t md_monomer_from_comp_name(str_t comp_name, md_component_kind_t kind);

#ifdef __cplusplus
}
#endif

// macro concatenate trick to assert that the input is a valid compile time C-string
#define MAKE_LABEL(cstr) {cstr"", sizeof(cstr)-1}

#ifdef __cplusplus
#define LBL_TO_STR(lbl) {lbl.buf, lbl.len}
inline bool operator==(const md_label_t& a, const md_label_t& b) {
    return MEMCMP(&a, &b, sizeof(md_label_t)) == 0;
}
inline bool operator!=(const md_label_t& a, const md_label_t& b) {
    return MEMCMP(&a, &b, sizeof(md_label_t)) != 0;
}
#else
#define LBL_TO_STR(lbl) (str_t){lbl.buf, lbl.len}
#endif // __cplusplus

static inline bool label_empty(md_label_t lbl) {
    uint64_t ref = 0;
    return MEMCMP(&lbl, &ref, sizeof(md_label_t)) == 0;
}

static inline md_label_t make_label(str_t str) {
    md_label_t lbl = {0};
    if (str.ptr) {
        const size_t len = MIN(str.len, sizeof(lbl.buf) - 1);
        MEMCPY(lbl.buf, str.ptr, len);
        lbl.len = (uint8_t)len;
    }
    return lbl;
}

// Access to substructure data
static inline void md_index_data_free (md_index_data_t* data) {
    ASSERT(data);
    if (data->offsets) {
        ASSERT(data->alloc);
        md_array_free(data->offsets,  data->alloc);
    }
    if (data->indices) {
        ASSERT(data->alloc);
        md_array_free(data->indices,  data->alloc);
    }
#if DEBUG
    MEMSET(data, 0, sizeof(md_index_data_t));
#endif
}

static inline size_t md_index_data_push_arr (md_index_data_t* data, const int32_t* index_data, size_t index_count) {
    ASSERT(data);
    ASSERT(data->alloc);
    ASSERT(index_count > 0);

    if (md_array_size(data->offsets) == 0) {
        md_array_push(data->offsets, 0, data->alloc);
    }

    size_t offset = *md_array_last(data->offsets);
    if (index_count > 0) {
        if (index_data) {
            md_array_push_array(data->indices, index_data, index_count, data->alloc);
        } else {
            md_array_grow(data->indices, md_array_size(data->indices) + index_count, data->alloc);
        }
        offset = md_array_size(data->indices);
        md_array_push(data->offsets, (uint32_t)offset, data->alloc);
    }

    return offset;
}

static inline void md_index_data_clear(md_index_data_t* data) {
    ASSERT(data);
    md_array_shrink(data->offsets, 0);
    md_array_shrink(data->indices, 0);
}

// The offsets of a cleared index data are empty rather than NULL: it has no ranges, as one never pushed to
static inline size_t md_index_data_num_ranges(const md_index_data_t* data) {
    ASSERT(data);
    const size_t num_offsets = md_array_size(data->offsets);
    return num_offsets > 0 ? num_offsets - 1 : 0;
}

// Access to individual substructures
static inline const int32_t* md_index_range_beg(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    const int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx + 1 < md_array_size(data->offsets)) {
        ptr = data->indices + data->offsets[range_idx];
    }
    return ptr;
}

static inline const int32_t* md_index_range_end(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    const int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx + 1 < md_array_size(data->offsets)) {
        ptr = data->indices + data->offsets[range_idx + 1];
    }
    return ptr;
}

// This is a mutable version of md_index_range_beg
static inline int32_t* md_index_range_ptr(md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    int32_t* ptr = NULL;
    if (data->indices && data->offsets && range_idx + 1 < md_array_size(data->offsets)) {
        ptr = data->indices + data->offsets[range_idx];
    }
    return ptr;
}

static inline size_t md_index_range_size(const md_index_data_t* data, size_t range_idx) {
    ASSERT(data);
    size_t size = 0;
    if (data->offsets && range_idx + 1 < md_array_size(data->offsets)) {
        size = data->offsets[range_idx+1] - data->offsets[range_idx];
    }
    return size;
}

static inline void md_index_data_merge(md_index_data_t* dest, const md_index_data_t* src) {
    ASSERT(dest);
    ASSERT(src);
    ASSERT(dest->alloc == src->alloc);
    size_t dest_num_indices = md_array_size(dest->indices);
    size_t src_num_ranges = md_index_data_num_ranges(src);
    for (size_t i = 0; i < src_num_ranges; ++i) {
        uint32_t beg_offset = src->offsets[i];
        uint32_t end_offset = src->offsets[i + 1];
        size_t range_size = end_offset - beg_offset;
        // Push adjusted offsets
        uint32_t new_offset = (uint32_t)(dest_num_indices + range_size);
        md_array_push(dest->offsets, new_offset, dest->alloc);
        // Push indices
        for (uint32_t j = beg_offset; j < end_offset; ++j) {
            int32_t idx = src->indices[j];
            md_array_push(dest->indices, idx, dest->alloc);
        }
        dest_num_indices += range_size;
    }
}

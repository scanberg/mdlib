#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_common.h>
#include <core/md_str.h>

// The index types of a system and the classification of its residues and parts. The rest of what a system stores is
// declared with what stores it: md_atom.h, md_hierarchy.h, md_bond.h, md_backbone.h, md_lipid.h, md_ring.h, and the
// elements in md_element.h.

typedef int32_t     md_atom_idx_t;
typedef int32_t     md_residue_idx_t;
typedef int32_t     md_instance_idx_t;
typedef int32_t     md_chain_idx_t;
typedef int32_t     md_part_idx_t;
typedef int32_t     md_part_type_idx_t;
typedef int32_t     md_molecule_idx_t;
typedef int32_t     md_molecule_type_idx_t;
typedef int32_t     md_backbone_idx_t;
typedef int32_t     md_sequence_id_t;
typedef int32_t     md_bond_idx_t;
typedef uint16_t    md_particle_type_idx_t;
typedef int32_t     md_residue_type_idx_t;
typedef uint8_t     md_atomic_number_t;      // Atomic number (1-118), 0 is unknown

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

// CLASSIFICATION
//
// What a part of a system is, is stored at the level of the hierarchy it describes, and only there:
//
//   particle type md_particle_type_flags_t what the particle is: an atom, a coarse grained bead or a virtual site
//   atom       md_atom_flags_t       the atom's role in its residue (backbone, side chain, ...) and its chemistry
//   residue    md_residue_flags_t    what the residue is (amino acid, nucleotide, water, ion), and where in its chain
//   part type  md_part_type_flags_t  what the part is (peptide, DNA, water, ...)
//
// Residues of the same name and kind share a residue type (md_residue_type_data_t), as atoms share a particle
// type. What a particle type is a form of is its element, what a residue type is a form of is its standard monomer
// (md_monomer_t): HIE, HSD and MSE are forms of HIS, HIS and MET. Both are where the defaults of the type come from.
//
// A part has no flags of its own: what it is is the kind of its type (md_system_part_kind).
// Nothing is copied from one level to another. A property of an atom's particle type or of its residue is read through the
// type or the residue (md_atom_particle_kind, md_system_atom_residue_kind), so there is one copy and it cannot go
// stale when the type or the residue is reclassified.
//
// Options which exclude each other are a VALUE, packed into a few bits and read and written through the accessors
// below (md_residue_flags_kind, md_residue_flags_set_kind, ...). Independent properties are single bits. No
// meaning is carried by a combination of bits.

// ### RESIDUE ###

// What a residue is, see md_util_system_infer_residue_flags.
typedef enum md_residue_kind_t {
    MD_RESIDUE_KIND_OTHER           = 0,    // Unclassified, or none of the below: a ligand, a sugar, ...
    MD_RESIDUE_KIND_AMINO_ACID      = 1,
    MD_RESIDUE_KIND_NUCLEOTIDE      = 2,
    MD_RESIDUE_KIND_WATER           = 3,
    MD_RESIDUE_KIND_ION             = 4,    // A monatomic ion
    MD_RESIDUE_KIND_LIPID           = 5,    // A lipid, sterols included: by name (md_util_resname_lipid), or by its structure
                                            // (a polar head with hydrocarbon tails, see md_util_system_infer_lipids)
} md_residue_kind_t;

typedef enum {
    MD_RESIDUE_FLAG_NONE                = 0,
    MD_RESIDUE_FLAG_KIND_MASK           = 0x7,      // md_residue_kind_t, see md_residue_flags_kind

    // The backbone of the amino acid or nucleotide was found by atom names and verified by its bonds, and its atoms
    // were given their roles (MD_ATOM_FLAG_BACKBONE, ...). Without it the residue is one by name only: an unusual
    // naming scheme, missing atoms, or a coarse grained model.
    MD_RESIDUE_FLAG_RESOLVED            = 0x10,
    MD_RESIDUE_FLAG_TERMINAL_BEG        = 0x20,     // Begins its chain: a free amino group (N-terminus), a free 5' end
    MD_RESIDUE_FLAG_TERMINAL_END        = 0x40,     // Ends its chain: a free carboxyl group (C-terminus), a free 3' end

    // The chain is broken next to it: its backbone is not linked to that of the residue before (after) it in its chain,
    // which carries BREAK_AFTER (BREAK_BEFORE) in turn. Residues missing in a crystal structure leave such a break. See
    // CHAIN BREAKS in md_util.h.
    MD_RESIDUE_FLAG_BREAK_BEFORE        = 0x80,
    MD_RESIDUE_FLAG_BREAK_AFTER         = 0x100,
} md_residue_flags_t;

MD_ENUM_FLAGS(md_residue_flags_t)

// ### MONOMER ###

// The standard monomers of proteins and nucleic acids, the predefined table residue types are matched against
// (md_monomer_from_residue_name), as particle types are matched against the elements. A monomer has its id in the Chemical
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

// ### PART TYPE, MOLECULE TYPE ###

// What a part is (md_part_data_t): the kind of its type. A part is one polymer strand, one glycan tree, one small
// molecule (a ligand, an ion, a lipid) or one water.
typedef enum md_part_kind_t {
    MD_PART_KIND_UNKNOWN      = 0,    // Not classified yet: a type the data names without saying what it is.
                                      // md_util_system_infer_parts classifies these from their residues.
    MD_PART_KIND_NON_POLYMER  = 1,    // A small molecule: a ligand, an ion, a lipid, a single residue
    MD_PART_KIND_WATER        = 2,
    MD_PART_KIND_BRANCHED     = 3,    // A branched oligosaccharide (mmCIF 'branched')
    // Polymers from here on, see md_part_kind_is_polymer
    MD_PART_KIND_PEPTIDE      = 4,    // A polypeptide chain (protein)
    MD_PART_KIND_DNA          = 5,
    MD_PART_KIND_RNA          = 6,
    MD_PART_KIND_NUCLEIC      = 7,    // A nucleic acid which is both (a DNA/RNA hybrid), or which could not be told apart
    MD_PART_KIND_POLYMER      = 8,    // Any other polymer
} md_part_kind_t;

typedef enum {
    MD_PART_TYPE_FLAG_NONE         = 0,
    MD_PART_TYPE_FLAG_KIND_MASK    = 0xF,      // md_part_kind_t, see md_part_type_flags_kind

    // Not given by the data (an mmCIF entity is): the type was inferred with the parts (md_util_system_infer_parts),
    // and is inferred anew with them
    MD_PART_TYPE_FLAG_INFERRED     = 0x10,
} md_part_type_flags_t;

MD_ENUM_FLAGS(md_part_type_flags_t)

typedef enum {
    MD_MOLECULE_TYPE_FLAG_NONE     = 0,

    // Not given by the data (a molecule type of a topology is): the type was inferred with the molecules
    // (md_util_system_infer_molecules), and is inferred anew with them
    MD_MOLECULE_TYPE_FLAG_INFERRED = 0x1,
} md_molecule_type_flags_t;

MD_ENUM_FLAGS(md_molecule_type_flags_t)

static inline bool md_part_kind_is_polymer(md_part_kind_t kind) {
    return kind >= MD_PART_KIND_PEPTIDE;
}

static inline bool md_part_kind_is_nucleic_acid(md_part_kind_t kind) {
    return kind == MD_PART_KIND_DNA || kind == MD_PART_KIND_RNA || kind == MD_PART_KIND_NUCLEIC;
}

static inline md_residue_kind_t md_residue_flags_kind(md_residue_flags_t flags) {
    return (md_residue_kind_t)((unsigned int)flags & (unsigned int)MD_RESIDUE_FLAG_KIND_MASK);
}

static inline md_residue_flags_t md_residue_flags_set_kind(md_residue_flags_t flags, md_residue_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_RESIDUE_FLAG_KIND_MASK;
    return (md_residue_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

static inline md_part_kind_t md_part_type_flags_kind(md_part_type_flags_t flags) {
    return (md_part_kind_t)((unsigned int)flags & (unsigned int)MD_PART_TYPE_FLAG_KIND_MASK);
}

static inline md_part_type_flags_t md_part_type_flags_set_kind(md_part_type_flags_t flags, md_part_kind_t kind) {
    const unsigned int mask = (unsigned int)MD_PART_TYPE_FLAG_KIND_MASK;
    return (md_part_type_flags_t)(((unsigned int)flags & ~mask) | ((unsigned int)kind & mask));
}

#ifdef __cplusplus
extern "C" {
#endif

const char* md_residue_kind_name(md_residue_kind_t kind);   // "" for a value out of range
const char* md_part_kind_name(md_part_kind_t kind);

// Standard monomers (md_monomer_t), "" / 0 for MD_MONOMER_NONE or a value out of range
md_str_t    md_monomer_id(md_monomer_t monomer);          // Id in the Chemical Component Dictionary: "ALA", "DA", "A"
md_str_t    md_monomer_name(md_monomer_t monomer);        // "Alanine", "Deoxyadenosine", "Adenosine"
char     md_monomer_letter(md_monomer_t monomer);      // One letter code: 'A' for ALA, DA and A
md_residue_kind_t md_monomer_kind(md_monomer_t monomer);    // MD_RESIDUE_KIND_AMINO_ACID or _NUCLEOTIDE
uint32_t md_monomer_color(md_monomer_t monomer);       // Default color (RGBA, R in the lowest byte, as the CPK colors)

// The monomer of the given id ("HIS"), exactly and ignoring case. MD_MONOMER_NONE if there is none.
md_monomer_t md_monomer_from_id(md_str_t id);

// The standard monomer a residue of the given name and kind is a form of. Only amino acids and nucleotides are forms
// of one: the variants of the force fields (HID, HIE, HSP, HISH, CYX, ASH, LYN, the Martini ARG0, ...), the termini of
// Amber (NALA, CALA, DA5, RA3), the names of CHARMM (ADE, CYT, ...) and the common modified residues (MSE, SEP, PSU,
// 5MC, ...). An amino acid or nucleotide which is not known is MD_MONOMER_UNK or MD_MONOMER_N, any other kind
// MD_MONOMER_NONE, except the caps of a peptide chain (ACE, NME, ...) which are MD_MONOMER_UNK.
md_monomer_t md_monomer_from_residue_name(md_str_t res_name, md_residue_kind_t kind);

#ifdef __cplusplus
}
#endif

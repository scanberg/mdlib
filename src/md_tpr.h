#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_types.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
typedef struct md_system_state_t md_system_state_t;

// GROMACS portable run input (.tpr)
//
// A tpr is what grompp writes and mdrun reads: the complete topology together with the starting
// coordinates, so unlike a .gro it describes the system rather than just where the atoms are.
// What is read here is the part a viewer can use:
//   - the molecule types (atoms with name, type, mass, charge, element, residue) and their bonds
//   - the molecule blocks that lay the molecule types out into the system
//   - the box, the coordinates and (when present) the velocities
// Force field parameters, exclusions, the simulation parameters and so on are skipped.
//
// Every file version GROMACS itself can still read is supported: tpx version 58 (GROMACS 4.0) and
// later, single and double precision, both the old layout (everything XDR encoded) and the one
// used since GROMACS 2020 (a big endian, unpadded body behind an XDR header). Newer files are read
// as long as their topology generation is one this parser knows, which is the same promise GROMACS
// makes for reading the topology of a file written by a newer version.
//
// Units are the file's own: nm, ps, e and u.

// GROMACS particle types
typedef enum md_tpr_ptype_t {
    MD_TPR_PTYPE_ATOM    = 0,
    MD_TPR_PTYPE_NUCLEUS = 1,
    MD_TPR_PTYPE_SHELL   = 2,
    MD_TPR_PTYPE_BOND    = 3,
    MD_TPR_PTYPE_VSITE   = 4,
} md_tpr_ptype_t;

// GROMACS periodic boundary types, as stored in the file
typedef enum md_tpr_pbc_t {
    MD_TPR_PBC_XYZ   = 0,
    MD_TPR_PBC_NO    = 1,
    MD_TPR_PBC_XY    = 2,
    MD_TPR_PBC_SCREW = 3,
    MD_TPR_PBC_UNSET = 4,
} md_tpr_pbc_t;

typedef struct md_tpr_atom_t {
    str_t    name;
    str_t    type;           // Force field atom type (e.g. 'opls_135')
    float    mass;
    float    charge;
    int32_t  residue;        // Index into the molecule type's residues
    int32_t  atomic_number;  // -1 when the force field does not give one (common for coarse grained)
    uint16_t type_idx;       // Index of the non-bonded type
    uint8_t  ptype;          // md_tpr_ptype_t
} md_tpr_atom_t;

typedef struct md_tpr_residue_t {
    str_t   name;
    int32_t nr;              // Residue number as given in the topology
    char    ic;              // Insertion code, ' ' when there is none
} md_tpr_residue_t;

typedef struct md_tpr_moltype_t {
    str_t name;
    size_t num_atoms;
    md_tpr_atom_t* atoms;
    size_t num_residues;
    md_tpr_residue_t* residues;
    // Connectivity in local, zero based atom indices, idx[0] < idx[1], sorted and unique: the pairs
    // GROMACS itself considers chemically bonded (bonds of every type that connects, constraints and
    // SETTLE), plus every virtual site to its first constructing atom, so a virtual site belongs to
    // the structure of the molecule it is part of.
    size_t num_bonds;
    md_atom_pair_t* bonds;
} md_tpr_moltype_t;

// Lennard-Jones parameters of a non-bonded type with itself, in kJ/mol nm^6 and kJ/mol nm^12.
// Both zero when the type has none (a charge-only site) or the force field is not Lennard-Jones.
typedef struct md_tpr_lj_t {
    float c6;
    float c12;
} md_tpr_lj_t;

// nmol consecutive molecules of type moltype
typedef struct md_tpr_molblock_t {
    int32_t moltype;
    int32_t nmol;
} md_tpr_molblock_t;

typedef struct md_tpr_data_t {
    int32_t file_version;     // tpx version
    int32_t file_generation;  // tpx topology generation
    bool    double_precision;
    str_t   version_string;   // e.g. 'VERSION 2023.3'
    str_t   name;             // The [ system ] name

    // Indexed by md_tpr_atom_t::type_idx
    size_t num_nb_types;
    md_tpr_lj_t* lj;

    size_t num_moltypes;
    md_tpr_moltype_t* moltypes;
    size_t num_molblocks;
    md_tpr_molblock_t* molblocks;

    // Bonds between atoms of different molecules ([ intermolecular_interactions ]), in global atom indices
    size_t num_intermolecular_bonds;
    md_atom_pair_t* intermolecular_bonds;

    size_t num_atoms;
    bool   has_box;
    float  box[3][3];         // Box vectors as rows
    int32_t pbc;              // md_tpr_pbc_t, MD_TPR_PBC_UNSET when the file does not say

    float* x;                 // num_atoms * 3, NULL if the file has no coordinates
    float* v;                 // num_atoms * 3, NULL if the file has no velocities

    // Backing storage of every string above, one allocation
    char*  str_data;
    size_t str_size;
} md_tpr_data_t;

bool md_tpr_data_parse_buffer(md_tpr_data_t* data, const void* buffer, size_t size, struct md_allocator_i* alloc);
bool md_tpr_data_parse_file(md_tpr_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_tpr_data_free(md_tpr_data_t* data, struct md_allocator_i* alloc);

// The van der Waals radius, in Ångström, implied by Lennard-Jones parameters: half the distance
// of the potential minimum, 2^(1/6) sigma / 2. Zero when there is no attraction to define it.
float md_tpr_lj_vdw_radius(md_tpr_lj_t lj);

// Builds a system from the tpr: atoms, residues, bonds (flagged MD_BOND_FLAG_TOPOLOGY), coordinates
// and box. Residues are numbered the way gmx numbers them when it writes the system out (single
// residue molecules such as water and ions are renumbered consecutively), so the numbers agree with
// a .gro written from the same tpr.
//
// Atom types are the particle types of the topology: atoms with the same name, element, force field
// type, mass and particle type share one. The type's mass is therefore exact for every atom of it,
// and the force field type is kept on it (md_atom_type_ff_type) to tell same named types apart.
// Atoms with an atomic number get their element's radius. Atoms without one are coarse grained
// beads (flagged MD_FLAG_COARSE_GRAINED) with the radius md_tpr_lj_vdw_radius gives for their
// non-bonded type; the predefined bead tables then add what they know about them (see
// md_util_system_augment_atom_types). A virtual site without Lennard-Jones parameters, such as the
// M site of TIP4P, has no element and is not a bead.
//
// Published attributes: 'atom/charge' (e) and, when present, 'atom/velocity' (nm/ps).
// The bonds are the topology's and are complete, so a caller should not infer covalent bonds on
// top of them.
bool md_tpr_system_init_from_data(struct md_system_t* sys, md_system_state_t* state, const md_tpr_data_t* data);
bool md_tpr_system_init_from_file(struct md_system_t* sys, md_system_state_t* state, str_t filename);

#ifdef __cplusplus
}
#endif

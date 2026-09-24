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
//   - the molecule types (atoms with name, type, mass, charge, element, residue), their bonds and
//     their non-bonded exclusions
//   - the molecule blocks that lay the molecule types out into the system
//   - the Lennard-Jones parameters of every pair of non-bonded types
//   - how mdrun computes the non-bonded interactions: potentials, modifiers and cut-offs
//   - the box, the coordinates and (when present) the velocities
// Other force field parameters (bonded ones, 1-4 pairs) and the rest of the simulation parameters
// are skipped.
//
// Every file version GROMACS itself can still read is supported: tpx version 58 (GROMACS 4.0) and
// later, single and double precision, both the old layout (everything XDR encoded) and the one
// used since GROMACS 2020 (a big endian, unpadded body behind an XDR header). Newer files are read
// as long as their topology generation is one this parser knows, which is the same promise GROMACS
// makes for reading the topology of a file written by a newer version.
//
// Units are the file's own: nm, ps, e, u and kJ/mol.

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

// The enumerations of the simulation parameters, as stored in the file (GROMACS' md_enums.h). The
// values have been stable since GROMACS 4.6; the reader converts the older files which differ.
typedef enum md_tpr_cutoff_scheme_t {
    MD_TPR_CUTOFF_SCHEME_VERLET = 0,
    MD_TPR_CUTOFF_SCHEME_GROUP  = 1,
} md_tpr_cutoff_scheme_t;

typedef enum md_tpr_vdw_type_t {
    MD_TPR_VDW_CUT         = 0,
    MD_TPR_VDW_SWITCH      = 1,     // Group scheme only: grompp turns it into CUT + POT_SWITCH for Verlet
    MD_TPR_VDW_SHIFT       = 2,     // Group scheme only: grompp turns it into CUT + FORCE_SWITCH for Verlet
    MD_TPR_VDW_USER        = 3,     // Tabulated: the tables are files given to mdrun, not part of the tpr
    MD_TPR_VDW_ENCAD_SHIFT = 4,     // Obsolete
    MD_TPR_VDW_PME         = 5,     // LJ-PME
} md_tpr_vdw_type_t;

typedef enum md_tpr_coulomb_type_t {
    MD_TPR_COULOMB_CUT             = 0,
    MD_TPR_COULOMB_RF              = 1,     // Reaction field
    MD_TPR_COULOMB_GRF             = 2,     // Obsolete
    MD_TPR_COULOMB_PME             = 3,
    MD_TPR_COULOMB_EWALD           = 4,
    MD_TPR_COULOMB_P3M_AD          = 5,
    MD_TPR_COULOMB_POISSON         = 6,
    MD_TPR_COULOMB_SWITCH          = 7,
    MD_TPR_COULOMB_SHIFT           = 8,
    MD_TPR_COULOMB_USER            = 9,     // Tabulated: the tables are files given to mdrun, not part of the tpr
    MD_TPR_COULOMB_GB              = 10,    // Obsolete
    MD_TPR_COULOMB_RF_NEC          = 11,    // Obsolete
    MD_TPR_COULOMB_ENCAD_SHIFT     = 12,    // Obsolete
    MD_TPR_COULOMB_PME_USER        = 13,
    MD_TPR_COULOMB_PME_SWITCH      = 14,
    MD_TPR_COULOMB_PME_USER_SWITCH = 15,
    MD_TPR_COULOMB_RF_ZERO         = 16,
    MD_TPR_COULOMB_FMM             = 17,
} md_tpr_coulomb_type_t;

// How a potential is brought to zero at its cut-off
typedef enum md_tpr_modifier_t {
    MD_TPR_MODIFIER_POT_SHIFT_VERLET_UNSUPPORTED = 0,
    MD_TPR_MODIFIER_POT_SHIFT    = 1,   // V(r) - V(rc)
    MD_TPR_MODIFIER_NONE         = 2,   // Plain cut-off
    MD_TPR_MODIFIER_POT_SWITCH   = 3,   // The potential switched to zero from the switch distance
    MD_TPR_MODIFIER_EXACT_CUTOFF = 4,
    MD_TPR_MODIFIER_FORCE_SWITCH = 5,   // The force switched to zero from the switch distance
} md_tpr_modifier_t;

typedef enum md_tpr_disp_corr_t {
    MD_TPR_DISP_CORR_NO            = 0,
    MD_TPR_DISP_CORR_ENER_PRES     = 1,
    MD_TPR_DISP_CORR_ENER          = 2,
    MD_TPR_DISP_CORR_ALL_ENER_PRES = 3,
    MD_TPR_DISP_CORR_ALL_ENER      = 4,
} md_tpr_disp_corr_t;

// How mdrun computes the non-bonded interactions of the system, which is what turns the parameters of a
// pair of types into a potential. Distances in nm. With the Verlet scheme, grompp has already replaced the
// group scheme's vdw switch and shift by a cut-off with a modifier.
typedef struct md_tpr_nonbonded_t {
    // False when the file has no simulation parameters, or they are of a version newer than this reader
    // knows how to read: then nothing below is set.
    bool    valid;

    int32_t cutoff_scheme;      // md_tpr_cutoff_scheme_t
    float   rlist;              // Pair list cut-off, including the buffer

    int32_t vdw_type;           // md_tpr_vdw_type_t
    int32_t vdw_modifier;       // md_tpr_modifier_t
    float   rvdw_switch;        // Where a switch starts
    float   rvdw;               // Cut-off

    int32_t coulomb_type;       // md_tpr_coulomb_type_t
    int32_t coulomb_modifier;   // md_tpr_modifier_t
    float   rcoulomb_switch;
    float   rcoulomb;
    float   epsilon_r;          // Relative dielectric constant (0 means infinity)
    float   epsilon_rf;         // Reaction field dielectric constant (0 means infinity)

    int32_t disp_corr;          // md_tpr_disp_corr_t. A mean field correction of the energy beyond rvdw, not per pair.

    float   ewald_rtol;         // Relative strength of the Coulomb potential at rcoulomb (PME, Ewald): sets the splitting
    float   ewald_rtol_lj;      // The same for LJ-PME at rvdw
    int32_t ljpme_comb_rule;    // LJ-PME grid combination rule: 0 geometric, 1 Lorentz-Berthelot
} md_tpr_nonbonded_t;

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
    // Non-bonded exclusions in local, zero based atom indices, as grompp wrote them: the atoms within nrexcl
    // bonds and the explicit [ exclusions ]. The atoms excluded from atom i are excl[excl_offset[i]] ..
    // excl[excl_offset[i + 1] - 1], ascending, without i itself. Symmetric. NULL when the molecule type has none.
    uint32_t* excl_offset;      // [num_atoms + 1]
    uint32_t* excl;
} md_tpr_moltype_t;

// Lennard-Jones parameters of a pair of non-bonded types, in kJ/mol nm^6 and kJ/mol nm^12.
// Both zero when the pair has none (a charge-only site) or the force field is not Lennard-Jones.
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

    // Lennard-Jones parameters of every pair of non-bonded types (md_tpr_atom_t::type_idx): a symmetric,
    // row major num_nb_types x num_nb_types table, lj[a * num_nb_types + b]. grompp resolves combination
    // rules and [ nonbond_params ] into this table, so an entry off the diagonal is the force field's own
    // value for the pair, which in general is not a combination of the two diagonal entries.
    // All zero when the non-bonded interactions are not Lennard-Jones (nb_is_lj false, e.g. Buckingham).
    size_t num_nb_types;
    md_tpr_lj_t* lj;
    bool nb_is_lj;
    double repulsion_power;   // The exponent of the repulsion: 12, unless the force field says otherwise
    float  fudge_qq;          // Scaling of the Coulomb interaction of 1-4 pairs

    md_tpr_nonbonded_t nonbonded;

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

// The Lennard-Jones parameters of a pair of non-bonded types, zero for types beyond the table
static inline md_tpr_lj_t md_tpr_lj_pair(const md_tpr_data_t* data, size_t type_a, size_t type_b) {
    if (!data || !data->lj || type_a >= data->num_nb_types || type_b >= data->num_nb_types) {
        md_tpr_lj_t zero = {0};
        return zero;
    }
    return data->lj[type_a * data->num_nb_types + type_b];
}

// Whether the non-bonded interaction between two atoms (global, zero based indices) is excluded. Exclusions are
// only ever between atoms of the same molecule, and an atom is excluded from itself.
bool md_tpr_atoms_excluded(const md_tpr_data_t* data, size_t atom_a, size_t atom_b);

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

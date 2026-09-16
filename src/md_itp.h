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

// GROMACS topology (.top / .itp)
//
// A topology does not describe coordinates, so it never creates a system. It SUPPLEMENTS one that a
// structure file (gro, pdb, ...) already loaded: its molecule definitions are matched onto the
// system's atoms by atom name, and what the structure file could not carry is filled in from them -
// bonds first of all, which is what gives a coarse grained system its structures (and so unwrap).
//
// Parsed: [ atomtypes ], [ moleculetype ], [ atoms ], [ bonds ], [ constraints ], [ settles ],
// [ virtual_sites1..4 ] / [ virtual_sitesn ] (and the old [ dummies* ] names), [ molecules ], [ system ].
// Everything else (parameters, angles, dihedrals, pairs, exclusions, restraints) is skipped.
// Preprocessor: #include "file" (relative to the including file, then $GMXLIB), #define, #undef,
// #ifdef, #ifndef, #else, #endif. Macros are only tested for being defined, never substituted - no
// field read here is ever a macro.

typedef struct md_itp_atomtype_t {
    str_t   name;
    int32_t atomic_number;  // -1 when the line has no atomic number column
    float   mass;
    float   charge;
} md_itp_atomtype_t;

typedef struct md_itp_atom_t {
    str_t   type;
    str_t   name;
    str_t   res_name;
    int32_t res_nr;
    float   charge;         // Only meaningful when has_charge
    float   mass;           // Only meaningful when has_mass
    bool    has_charge;
    bool    has_mass;
} md_itp_atom_t;

typedef struct md_itp_moleculetype_t {
    str_t   name;
    int32_t nrexcl;
    md_itp_atom_t*  atoms;    // md_array
    // Connectivity in local, zero based atom indices, i < j, sorted and unique. The union of [ bonds ]
    // (every function type), [ constraints ], [ settles ] and virtual sites to their first constructing
    // atom: everything that says two atoms belong together, which is what structures need.
    md_atom_pair_t* bonds;    // md_array
} md_itp_moleculetype_t;

typedef struct md_itp_molecules_t {
    str_t   name;
    int32_t count;
} md_itp_molecules_t;

typedef struct md_itp_data_t {
    str_t system_name;
    md_itp_atomtype_t*     atomtypes;      // md_array
    md_itp_moleculetype_t* moleculetypes;  // md_array
    md_itp_molecules_t*    molecules;      // md_array, empty for a plain .itp
} md_itp_data_t;

// One molecule of the topology placed onto the system: moleculetype 'type' covers the atoms
// [atom_offset, atom_offset + number of atoms in that moleculetype).
typedef struct md_itp_instance_t {
    uint32_t type;
    uint32_t atom_offset;
} md_itp_instance_t;

// Strings and arrays in 'data' are allocated from alloc. base_folder resolves #include in a string
// (may be empty, in which case only absolute includes and $GMXLIB are tried).
bool md_itp_data_parse_str (md_itp_data_t* data, str_t str, str_t base_folder, struct md_allocator_i* alloc);
bool md_itp_data_parse_file(md_itp_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_itp_data_free(md_itp_data_t* data, struct md_allocator_i* alloc);

// Places the topology's molecules onto the system's atoms.
// With a [ molecules ] section the molecules are laid down in that order from the first atom, and every
// one must match. Without one (a plain .itp), or when that fails, the system is scanned instead and each
// moleculetype is placed wherever its atom name sequence occurs, so an .itp describing only part of the
// system (just the cellulose, not the solvent) still applies to that part.
// Atom names match when equal, or when the system name is at least 4 characters and the topology name
// starts with it (structure formats truncate names). Returns the number of instances.
size_t md_itp_match_system(md_itp_instance_t** out_instances, const md_itp_data_t* data, const struct md_system_t* sys, struct md_allocator_i* alloc);

// Applies the topology to the system:
// - bonds: inside matched molecules the topology is authoritative, so inferred bonds between two matched
//   atoms are replaced by the topology's, flagged MD_BOND_FLAG_TOPOLOGY. Bonds touching unmatched atoms
//   and user defined bonds are kept.
// - 'atom/charge' (e) and 'atom/mass' (Da) are published for the matched atoms, taken from [ atoms ] and
//   falling back to [ atomtypes ]; unmatched atoms are left as gaps.
// - atom types with no mass (zero) get the topology's mass when all matched atoms of that type agree.
// - structures and rings are re-inferred from the new bonds.
// Returns false if nothing in the topology matched the system.
bool md_itp_system_supplement(struct md_system_t* sys, const md_itp_data_t* data);
bool md_itp_system_supplement_from_file(struct md_system_t* sys, str_t filename);

#ifdef __cplusplus
}
#endif

#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_io.h>

struct md_atom_pair_t;

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
struct md_system_state_t;

// GROMACS topology (.top / .itp)
//
// A topology is both a system of its own and data of a system's atoms (md_io.h, LOAD). Onto no atoms it builds the
// system of its molecules, without coordinates. Onto a system with atoms its molecule definitions are matched onto the
// atoms by atom name, and what the structure file could not carry is filled in from them - bonds first of all, which
// is what gives a coarse grained system its structures (and so unwrap). It may describe only part of the system.
//
// Parsed: [ atomtypes ], [ moleculetype ], [ atoms ], [ bonds ], [ constraints ], [ settles ],
// [ virtual_sites1..4 ] / [ virtual_sitesn ] (and the old [ dummies* ] names), [ molecules ], [ system ].
// Everything else (parameters, angles, dihedrals, pairs, exclusions, restraints) is skipped.
// Preprocessor: #include "file" (relative to the including file, then $GMXLIB), #define, #undef,
// #ifdef, #ifndef, #else, #endif. Macros are only tested for being defined, never substituted - no
// field read here is ever a macro.

typedef struct md_itp_atomtype_t {
    md_str_t   name;
    int32_t atomic_number;  // -1 when the line has no atomic number column
    float   mass;
    float   charge;
} md_itp_atomtype_t;

typedef struct md_itp_atom_t {
    md_str_t   type;
    md_str_t   name;
    md_str_t   res_name;
    int32_t res_nr;
    float   charge;         // Only meaningful when has_charge
    float   mass;           // Only meaningful when has_mass
    bool    has_charge;
    bool    has_mass;
} md_itp_atom_t;

typedef struct md_itp_moleculetype_t {
    md_str_t   name;
    int32_t nrexcl;
    md_itp_atom_t*  atoms;    // md_array
    // Connectivity in local, zero based atom indices, i < j, sorted and unique. The union of [ bonds ]
    // and [ constraints ] of the function types which connect, [ settles ] and virtual sites to their first
    // constructing atom: everything that says two atoms belong together, which is what structures need.
    struct md_atom_pair_t* bonds;    // md_array
    // Elastic bonds, the same way: the [ bonds ] of function types 6 (harmonic potential), 9 (tabulated, no
    // exclusions) and 10 (restraint potential) and the [ constraints ] of type 2 (no connection), which make no
    // connection: the elastic network of a coarse grained protein or nucleic acid, a restraint. A pair which is also
    // in bonds is not repeated here.
    struct md_atom_pair_t* elastic_bonds;    // md_array
} md_itp_moleculetype_t;

typedef struct md_itp_molecules_t {
    md_str_t   name;
    int32_t count;
} md_itp_molecules_t;

typedef struct md_itp_data_t {
    md_str_t system_name;
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
bool md_itp_data_parse_str (md_itp_data_t* data, md_str_t str, md_str_t base_folder, struct md_allocator_i* alloc);
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
// - bonds: the atoms of the matched molecules have the topology's bonds, of origin MD_BOND_ORIGIN_TOPOLOGY, and no
//   others (md_atom_data_t::bonds_complete): the bonds given between two of them are replaced, and none is inferred
//   between two of them. Bonds given which touch unmatched atoms and user defined bonds are kept. The elastic bonds of
//   the topology are bonds with MD_BOND_KIND_ELASTIC, and so are those of an elastic network written as plain bonds
//   (md_system_infer_elastic_bonds, by md_system_infer).
// - 'atom/charge' (e) and 'atom/mass' (Da) are published for the matched atoms, taken from [ atoms ] and
//   falling back to [ atomtypes ]; unmatched atoms are left as gaps.
// - particle types with no mass (zero) get the topology's mass when all matched atoms of that type agree.
// - molecules: where the structure gave no units of its own (a GRO file), the matched molecules become its instances,
//   each of its molecule type (md_instance_data_t, md_molecule_type_data_t), and the residues between them units of
//   no type. The molecules of the system are then those of the topology, as GROMACS has them.
// Nothing is inferred (md_io.h): the elastic network written as plain bonds is told by md_system_infer.
// Returns false, changing nothing, if nothing in the topology matched the system.
bool md_itp_system_supplement(struct md_system_t* sys, const md_itp_data_t* data);

// Onto no atoms: the system of the topology, its molecules those of [ molecules ], or one of each moleculetype for a
// plain .itp; the state has no coordinates. MD_LOAD_FAILED when [ molecules ] names a molecule it does not define (an
// #include not found). Onto atoms: md_itp_system_supplement, MD_LOAD_MISMATCH when nothing matched.
md_load_result_t md_itp_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

// Whether any of the topology's molecules match the system's atoms, changing nothing (md_io.h, CHECK)
md_load_result_t md_itp_system_check_file(const struct md_system_t* sys, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif

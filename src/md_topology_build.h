#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <md_types.h>
#include <md_unitcell.h>
#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// TOPOLOGY ONLY FORMATS
//
// What an atomistic topology without coordinates says (an Amber prmtop, a CHARMM/NAMD PSF), in one shape, and the
// two things done with it: building a system from it, and applying it to a system loaded from a structure file of
// the same atoms. The readers fill the description; everything after is the same for every format.

typedef struct md_topology_desc_t {
    const char* format;             // for the log, "Amber topology"

    size_t num_atoms;
    const md_str_t*   atom_name;       // required
    const md_str_t*   atom_type;       // force field types, NULL for none
    const double*  charge;          // e, NULL for none
    const double*  mass;            // Da, NULL for none
    const int64_t* atomic_number;   // NULL for none: inferred from names and masses. <= 0: unknown

    size_t num_residues;
    const md_str_t*    res_name;       // required
    const uint32_t* res_offset;     // num_residues + 1, required
    const int64_t*  res_number;     // NULL: numbered from 1
    const md_str_t*    res_chain;      // the chain or segment of each residue, NULL for none

    size_t num_bonds;
    const md_atom_pair_t* bonds;

    md_unitcell_t cell;             // the topology's own box, if any
    md_str_t         title;
} md_topology_desc_t;

// The system from the topology. The state receives no coordinates (md_system_state_has_coords is false), only the
// cell. Massless particles without an element are virtual sites, and bonds to them are left out (a site is tied to
// its atom by its construction). Published: 'atom/charge' (e) and 'atom/mass' (Da) where the topology has them.
bool md_topology_system_init(struct md_system_t* sys, struct md_system_state_t* state, const md_topology_desc_t* desc);

// Applies the topology to a loaded system of the same atoms in the same order: the atom count and every atom name
// (compared as structure files truncate them) must agree, and the residue boundaries where the residue counts do.
// The topology's bonds replace all but those the user made, 'atom/charge' and 'atom/mass' are published, and parts,
// molecules and backbones are inferred again. Nothing changes when they do not agree.
bool md_topology_system_supplement(struct md_system_t* sys, const md_topology_desc_t* desc, md_str_t filename);

#ifdef __cplusplus
}
#endif

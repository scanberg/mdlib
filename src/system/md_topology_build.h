#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_unitcell.h>
#include <core/md_str.h>

#include <md_io.h>

struct md_atom_pair_t;

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// TOPOLOGY ONLY FORMATS
//
// What an atomistic topology without coordinates says (an Amber prmtop, a CHARMM/NAMD PSF), in one shape, and its
// load (md_io.h, LOAD). The readers fill the description; everything after is the same for every format.

typedef struct md_topology_desc_t {
    const char* format;             // for the log, "Amber topology"

    size_t num_atoms;
    const md_str_t*   atom_name;    // required
    const md_str_t*   atom_type;    // force field types, NULL for none
    const double*  charge;          // e, NULL for none
    const double*  mass;            // Da, NULL for none
    const int64_t* atomic_number;   // > 0 the element, 0 none, < 0 inferred from name and mass. NULL: all inferred

    size_t num_residues;
    const md_str_t*    res_name;    // required
    const uint32_t* res_offset;     // num_residues + 1, required
    const int64_t*  res_number;     // NULL: numbered from 1
    const md_str_t*    res_chain;   // the chain or segment of each residue, NULL for none

    size_t num_bonds;
    const struct md_atom_pair_t* bonds;

    md_unitcell_t cell;             // the topology's own box, if any
    md_str_t         title;
} md_topology_desc_t;

// Builds the system of the topology onto a system without atoms (md_io.h, LOAD): its particle types, residues, bonds
// (those of every atom, so none is inferred: md_atom_data_t::bonds_complete), and chains where the topology has them.
// Massless particles without an element are virtual sites, and bonds to them are left out (a site is tied to its atom
// by its construction). The state is made for no atoms (md_system_state_has_coords is false): the coordinates come
// from another file. It holds the topology's cell. Published: 'atom/charge' (e) and 'atom/mass' (Da) where the
// topology has them. Not inferred.
//
// Onto a system with atoms it is MD_LOAD_FAILED, and nothing changes: the topology is a system of its own, which joins
// one of its atoms through md_system_merge.
md_load_result_t md_topology_system_load(struct md_system_t* sys, struct md_system_state_t* state, const md_topology_desc_t* desc, md_str_t filename);

#ifdef __cplusplus
}
#endif

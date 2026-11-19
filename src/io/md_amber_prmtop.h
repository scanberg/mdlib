#pragma once

#include <stdbool.h>
#include <stddef.h>

#include <core/md_str.h>
#include <md_io.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// AMBER TOPOLOGY (.prmtop, .parm7, .top)
//
// The %FLAG/%FORMAT topology of Amber 7 and later, as tleap and ParmEd write it, CHAMBER topologies (a CHARMM force
// field converted by chamber or ParmEd) included: what a system is made of is in the same blocks either way, and the
// force field terms which differ (CHARMM_UREY_BRADLEY, CHARMM_CMAP_*, ...) are not read.
//
// Read: ATOM_NAME, CHARGE (in units of e * 18.2223), MASS, ATOMIC_NUMBER, AMBER_ATOM_TYPE, RESIDUE_LABEL,
// RESIDUE_POINTER, BONDS_INC_HYDROGEN, BONDS_WITHOUT_HYDROGEN, BOX_DIMENSIONS and, as ParmEd writes them,
// RESIDUE_NUMBER and RESIDUE_CHAINID. Elements come from ATOMIC_NUMBER, and from names and masses in a topology
// without it. Massless particles without an element (TIP4P's M, the lone pairs of a force field) are virtual sites.
//
// A topology holds no coordinates. These come from another file - a restart (md_amber_rst7.h, or
// md_amber_netcdf.h), a structure file, or the first frame of a trajectory - in whatever order they are loaded.

// True when the file starts as a topology of this format does (%VERSION or %FLAG).
bool md_amber_prmtop_file_is_topology(md_str_t filename);

// The atom count of the topology (NATOM of POINTERS), 0 when the file is not one.
size_t md_amber_prmtop_file_atom_count(md_str_t filename);

// Loads the topology into the system (md_io.h, LOAD; md_topology_system_load): onto a system without atoms, its atoms
// without coordinates, and the cell of BOX_DIMENSIONS when the topology is periodic; onto a system of its atoms in its
// order (a PDB file of the same simulation), its bonds. Published: 'atom/charge' (e), 'atom/mass' (Da).
md_load_result_t md_amber_prmtop_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif

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

// PSF: the protein structure file of CHARMM, NAMD and X-PLOR (CHARMM's numbered atom types or X-PLOR's named ones,
// the standard or the EXT widths, VMD's output alike)
//
//     PSF EXT CMAP
//
//            2 !NTITLE
//      REMARKS ...
//
//         3341 !NATOM
//            1 U        1        MET      N        NH3     -0.300000       14.0070           0
//          ... id segment residue resname name type charge mass imove
//
//         3365 !NBOND: bonds
//            1       5       1       2 ...         pairs of atom ids, four pairs a line
//
// Read: the atoms (segments become the chains), their charges and masses, and the bonds. Elements are inferred from
// names and masses: a PSF has none. Massless particles (lone pairs) are virtual sites. The angles, dihedrals and the
// rest are not read.
//
// A PSF holds no coordinates. These come from another file - a PDB file of the same atoms, a NAMD binary coordinate
// file (md_namdbin.h), or the first frame of a DCD trajectory - in whatever order they are loaded.

bool md_psf_file_is_psf(md_str_t filename);

// The atom count of the file (!NATOM), 0 when the file is not a PSF
size_t md_psf_file_atom_count(md_str_t filename);

// Builds the system of the PSF (md_io.h, LOAD; md_topology_system_load): its atoms without coordinates, and its
// segments as chains. Published: 'atom/charge' (e), 'atom/mass' (Da). It joins a system of its atoms in its order (a PDB
// file of the same atoms) through md_system_merge.
md_load_result_t md_psf_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif

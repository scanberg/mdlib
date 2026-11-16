#pragma once

#include <stdbool.h>
#include <stddef.h>

#include <core/md_str.h>

struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// NAMD BINARY COORDINATES (.coor, binaryoutput on, the default): the atom count as a 32 bit integer, then x y z of
// every atom as doubles, in Angstrom, in the byte order of the machine that wrote it. Nothing else: no names, no
// cell. The atoms are those of the PSF (md_psf.h).
//
// NAMD's .vel restarts have the same layout; they are not coordinates, but nothing in the file says so.

// True when the size of the file is that of the atom count it starts with, in either byte order.
bool md_namdbin_file_is_coor(md_str_t filename);

size_t md_namdbin_file_atom_count(md_str_t filename);

// The coordinates into state; num_atoms 0 accepts any count. state->alloc must be set.
bool md_namdbin_state_init_from_file(struct md_system_state_t* state, md_str_t filename, size_t num_atoms);

#ifdef __cplusplus
}
#endif

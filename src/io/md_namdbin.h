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

// NAMD BINARY COORDINATES (.coor, binaryoutput on, the default): the atom count as a 32 bit integer, then x y z of
// every atom as doubles, in Angstrom, in the byte order of the machine that wrote it. Nothing else: no names, no
// cell. The atoms are those of the PSF (md_psf.h).
//
// NAMD's .vel restarts have the same layout; they are not coordinates, but nothing in the file says so.

// True when the size of the file is that of the atom count it starts with, in either byte order.
bool md_namdbin_file_is_coor(md_str_t filename);

size_t md_namdbin_file_atom_count(md_str_t filename);

// Builds the system of the file (md_io.h, LOAD): its atoms known by their count alone, their coordinates in state. It
// joins a topology of its atoms through md_system_merge.
md_load_result_t md_namdbin_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif

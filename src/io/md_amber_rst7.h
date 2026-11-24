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

// AMBER ASCII RESTART (.rst7, .inpcrd, .restrt, .rst)
//
//     title
//     NATOM [TIME]
//     coordinates      6F12.7, two atoms a line
//     [velocities]     6F12.7, as many lines
//     [box]            6F12.7: a b c alpha beta gamma (old files: a b c)
//
// Coordinates only: the atoms are the topology's (md_amber_prmtop.h), which says what they are. What follows the
// coordinates is told by its number of lines, which is unambiguous where the number of values is not. A restart
// in the NetCDF format (Conventions "AMBERRESTART") is read by md_amber_netcdf.h.

// True when the file is laid out as above: an atom count on the second line, coordinates in fields of 12 after.
bool md_amber_rst7_file_is_restart(md_str_t filename);

// NATOM of the second line, 0 when the file is not a restart.
size_t md_amber_rst7_file_atom_count(md_str_t filename);

// Builds the system of the restart (md_io.h, LOAD): its atoms known by their count alone, their coordinates (Angstrom)
// and cell in state. It joins a topology of its atoms through md_system_merge.
md_load_result_t md_amber_rst7_system_load_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t filename, const struct md_load_desc_t* desc);

#ifdef __cplusplus
}
#endif

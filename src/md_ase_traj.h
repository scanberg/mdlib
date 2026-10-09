#pragma once

#include <core/md_str.h>
#include <stdbool.h>
#include <stdint.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// ASE trajectories (.traj), the ULM files ase.io.Trajectory writes.
//
// What is read of one is what a system and a run can hold: a fixed set of atoms, and per frame their
// positions and a cell. The cell must be fully periodic with a along x and b in the xy plane, or not
// periodic at all, in which case it is dropped (it is then only the box ASE keeps around a
// molecule). Files whose atoms change, whose cell is rotated or periodic along only some axes, or
// that were written on a big endian machine are refused.

// SYSTEM
// The atoms, typed by their atomic numbers, with the positions and cell of the first frame, which is
// all that is read. The system is left as it was when the file cannot be read.
bool md_ase_traj_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);

// RUN
// Publishes the file as a run in the system's attribute table (see RUNS in md_system.h):
//
//     <run>/time              F64 {F}        ps, from info["time_ps"] when every frame has it; frame
//                                            ordinals without a unit otherwise
//     <run>/unitcell          F32 {F,3,3}    Angstrom, row i box vector i, zero for a frame without
//                                            periodicity
//     <run>/atom/position     F32 {F,N} c3   Angstrom, VIRTUAL: read from the file when a frame is
//                                            asked for, as float32 or float64 as it was written
//     <run>/source/path       STR            the file
//     <run>/source/offset     I64 {F}        byte offset of each frame's positions in it
//     <run>/source/size       I64 {F}        and their byte size, which tells float32 from float64
//
// False for a file of a single frame, which is a structure and not a trajectory. The header of every
// frame is read and checked; there is no index cache, since the file has a table of its frames, so
// flags has nothing to change. The provider reads through the io it is handed. Its
// user_data is sys, borrowed: the run goes before the system does. The atom count must match the
// system's.
bool md_ase_traj_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

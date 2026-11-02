#pragma once

#include <stdbool.h>
#include <stdint.h>

#include <core/md_str.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// AMBER NETCDF
//
// Trajectories (Conventions "AMBER") and restarts ("AMBERRESTART") as Amber, cpptraj, MDAnalysis,
// MDTraj and OpenMM write them: the classic NetCDF encodings, CDF-1 (32 bit offsets) and CDF-2 (64 bit
// offsets) which the convention requires, and CDF-5 for good measure. Read here directly: the classic
// format is a header followed by plain big endian arrays, so a frame is one read at a known offset, from
// any thread, without NetCDF-C.
//
// NetCDF-4 (HDF5 underneath, cpptraj's compressed "nc4") is outside the convention and recognized only
// to be refused with a clear message.
//
// See https://ambermd.org/netcdf/nctraj.xhtml

typedef enum md_amber_netcdf_kind_t {
    MD_AMBER_NETCDF_NONE = 0,       // not NetCDF, or NetCDF without the Amber convention
    MD_AMBER_NETCDF_TRAJECTORY,     // Conventions "AMBER"
    MD_AMBER_NETCDF_RESTART,        // Conventions "AMBERRESTART"
    MD_AMBER_NETCDF_HDF5,           // NetCDF-4 / HDF5: not supported
} md_amber_netcdf_kind_t;

// Tells from the header alone what the file is. Cheap: reads the first block of the file.
md_amber_netcdf_kind_t md_amber_netcdf_file_kind(str_t filename);

// The atom count of an Amber NetCDF trajectory or restart, 0 when it is neither.
size_t md_amber_netcdf_file_atom_count(str_t filename);

// Coordinates (Angstrom) and cell of a restart, or of the first frame of a trajectory, into state.
// num_atoms is what the system expects: the file must hold that many atoms (0 accepts any count).
// state->alloc must be set.
bool md_amber_netcdf_state_init_from_file(struct md_system_state_t* state, str_t filename, size_t num_atoms);

// RUN
// Publishes a trajectory as a run in the system's attribute table (see md_run_publish in md_system.h):
//
//     <run>/time, unitcell, source/...   as every run; time in ps, ordinals when the file has no time
//                                        (or one that does not increase), the cell zero without a box
//     <run>/atom/position   F32 {F,N} c3   Angstrom               VIRTUAL, read from the file when a
//     <run>/atom/velocity   F32 {F,N} c3   Angstrom/ps            frame is asked for, scale_factor
//     <run>/atom/force      F32 {F,N} c3   kcal/(mol Angstrom)    applied
//     <run>/remd/temp0      F64 {F}        K        replica exchange: the frame's target temperature
//     <run>/remd/repidx     I32 {F}                 ... the replica that wrote it
//     <run>/remd/crdidx     I32 {F}                 ... the coordinate index
//
// A file still being written, or cut short, publishes the frames that are complete. Providers read
// through the io they are handed and keep no state of their own; their user_data is sys, borrowed.
bool md_amber_netcdf_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

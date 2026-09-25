#pragma once

#include <core/md_str.h>

struct md_system_t;
struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// RUN
// Publishes the file as a run in the system's attribute table (see md_run_publish in md_system.h).
// The run's frames are the frames that carry coordinates:
//
//     <run>/time, step, unitcell, source/...   as every run; the cell in Angstrom, zero without a box
//     <run>/atom/position   F32 {F,N} c3   Angstrom            VIRTUAL, read from the file when a
//     <run>/atom/velocity   F32 {F,N} c3   Angstrom/ps         frame is asked for
//     <run>/atom/force      F32 {F,N} c3   kJ/(mol nm)
//
// Velocities and forces sit beside the positions when the file writes them in exactly the frames it
// writes coordinates. Written at another interval (nstvout, nstfout against nstxout) a section gets
// a group of its own with its own frame axis, so every value keeps the time it was written at:
//
//     <run>/trr/velocity/time, step, source/...   the frames that carry velocities
//     <run>/trr/velocity/atom/velocity            {F_v,N} c3
//     <run>/trr/force/...                         likewise
//
// Everything but the per atom data comes from the index kept beside the file ('<file>.cache',
// written unless flags carries MD_RUN_FLAG_DISABLE_CACHE_WRITE). The providers read only
// the section asked for, through the io they are handed: inside an extraction context the file
// stays open across frames. Their user_data is sys, borrowed: the run goes before the system does.
// The file's atom count must match the system's.
bool md_trr_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

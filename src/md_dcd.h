#pragma once

#include <core/md_str.h>

struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

// RUN
// Publishes the file as a run in the system's attribute table (see RUNS in md_system.h):
//
//     <run>/time              F64 {F}        ps; frame ordinals without a unit when the file has no
//                                            timestep
//     <run>/step              I64 {F}        istart + frame * nsavc
//     <run>/unitcell          F32 {F,3,3}    Angstrom, row i box vector i. VIRTUAL, read from each
//                                            frame, when the file carries a cell per frame; otherwise
//                                            the system's reference cell at every frame
//     <run>/atom/position     F32 {F,N} c3   Angstrom, VIRTUAL: read from the file when a frame is
//                                            asked for, the three planes interleaved
//     <run>/source/path       STR            the file
//     <run>/source/offset     I64 {F}        byte offset of each frame in it
//     <run>/source/size       I64 {F}        byte size of each frame
//     <run>/source/layout     I32 {3}        CHARMM flags, byte order reversed, number of fixed atoms
//     <run>/source/translation F32 c3        Angstrom, added to every coordinate: the cell's centre
//                                            when the first frame sits around it rather than around
//                                            the origin, zero otherwise
//     <run>/source/free_atoms I32 {N-fixed}  with fixed atoms only: the atoms frames after the first
//     <run>/source/first_frame F32 {N} c3    carry, and the first frame, which the others keep
//
// Nothing is scanned: DCD frames have a fixed size. The providers read through the io they are
// handed, so inside an extraction context the file stays open across frames. Their user_data is
// sys, borrowed: the run goes before the system does. The atom count must match the system's.
bool md_dcd_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

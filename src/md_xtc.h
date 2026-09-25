#pragma once

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_os.h>

struct md_allocator_i;
struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_xtc_header_t {
	int32_t natoms;
	int32_t step;
	float time;
	float box[3][3];
} md_xtc_header_t;

size_t md_xtc_read_frame_offsets_and_times(md_file_t xdr_file, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, struct md_allocator_i* alloc);

// This is an internal procedure exposed to enable testing and profiling.
// Returns the number of atoms decoded, or zero if the frame could not be decoded. The frame should be decoded into `out_header` and `out_xyz`, which should have capacity for at least three floats per atom.
// Note that the data is only decoded and length units are typically nm.
bool md_xtc_decode_frame_data(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_xyz, size_t num_atoms);
bool md_xtc_decode_frame_data_soa(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, size_t num_atoms);

// RUN
// Publishes the file as a run in the system's attribute table (see RUNS in md_system.h):
//
//     <run>/time            F64 {F}        ps
//     <run>/step            I64 {F}        the simulation step each frame was written at
//     <run>/unitcell        F32 {F,3,3}    Angstrom, the box of each frame, row i box vector i; zero
//                                          when the frame has no box
//     <run>/atom/position   F32 {F,N} c3   Angstrom, VIRTUAL: decoded from the file when a frame is
//                                          asked for, so nothing the size of the trajectory is held
//     <run>/source/path     STR            the file
//     <run>/source/offset   I64 {F}        byte offset of each frame in it
//     <run>/source/size     I64 {F}        byte size of each frame
//
// Everything but the coordinates comes from the index the trajectory keeps beside itself
// ('<file>.cache', written unless flags carries MD_RUN_FLAG_DISABLE_CACHE_WRITE), so this
// costs a scan of the headers at most and a cache read usually.
//
// The position provider reads the file through the source attributes, which is all the state it
// has, and through the io it is handed: inside an extraction context (md_system_extract_begin) the
// file stays open across frames, outside one it is opened per frame. Either way it is safe to call
// from any number of threads at once, each with its own context.
// Its user_data is sys, borrowed: the run has to be removed before the system goes, as it is when
// the prefix is removed with the trajectory. The file's atom count must match the system's.
bool md_xtc_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

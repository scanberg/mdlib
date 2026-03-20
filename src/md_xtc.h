#pragma once

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

// Returns the number of atoms decoded, or zero if the frame could not be decoded. The frame should be decoded into `out_header` and `out_xyz`, which should have capacity for at least three floats per atom.
// Note that the data is only decoded and length units are typically nm.
bool md_xtc_decode_frame_data(const uint8_t* frame_ptr, size_t frame_bytes, md_xtc_header_t* out_header, float* out_xyz, size_t num_atoms);

bool md_xtc_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags);

#ifdef __cplusplus
}
#endif

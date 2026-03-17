#pragma once

#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_os.h>

struct md_allocator_i;
struct md_trajectory_i;

#ifdef __cplusplus
extern "C" {
#endif

size_t md_xtc_read_frame_offsets_and_times(md_file_t  xdr_file, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, struct md_allocator_i* alloc);
bool   md_xtc_read_frame_header(md_file_t  xdr_file, int* natoms, int* step, float* time, float box[3][3]);
size_t md_xtc_read_frame_coords(md_file_t  xdr_file, float* out_xyz, size_t xyz_coords);

struct md_trajectory_i* md_xtc_trajectory_create(str_t filename, struct md_allocator_i* alloc, uint32_t flags);
void md_xtc_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_loader_i* md_xtc_trajectory_loader(void);

#ifdef __cplusplus
}
#endif

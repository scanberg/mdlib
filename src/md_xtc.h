#pragma once

#include <core/md_str.h>
#include <core/md_array.h>

struct md_allocator_i;
struct md_trajectory_i;
struct md_file_o;

#ifdef __cplusplus
extern "C" {
#endif

bool md_xtc_read_frame_offsets_and_times(struct md_file_o* xdr_file, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, struct md_allocator_i* alloc);
bool md_xtc_read_frame_header(struct md_file_o* xdr_file, int* natoms, int* step, float* time, float box[3][3]);
size_t md_xtc_read_frame_coords(struct md_file_o* xdr_file, float* out_xyz, size_t xyz_coords);

struct md_trajectory_i* md_xtc_trajectory_create(str_t filename, struct md_allocator_i* alloc, uint32_t flags);
void md_xtc_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_loader_i* md_xtc_trajectory_loader(void);

#ifdef __cplusplus
}
#endif

#pragma once

#include <core/md_str.h>

struct md_system_t;
#include <md_trajectory.h>

#ifdef __cplusplus
extern "C" {
#endif

bool md_dcd_attach_from_file(struct md_system_t* sys, str_t filename);

struct md_trajectory_i* md_dcd_trajectory_create(str_t filename, struct md_allocator_i* ext_alloc, uint32_t flags);
void md_dcd_trajectory_free(struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif

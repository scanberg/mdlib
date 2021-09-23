#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

struct md_allocator_i;
struct md_trajectory_i;

#ifdef __cplusplus
extern "C" {
#endif

bool md_xtc_trajectory_open(struct md_trajectory_i* traj, str_t filename, struct md_allocator_i* alloc);
bool md_xtc_trajectory_close(struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif

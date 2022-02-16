#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

struct md_allocator_i;
struct md_trajectory_i;

#ifdef __cplusplus
extern "C" {
#endif

struct md_trajectory_i* md_xtc_trajectory_create(str_t filename, struct md_allocator_i* alloc);
void md_xtc_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_api* md_xtc_trajectory_api(void);

#ifdef __cplusplus
}
#endif

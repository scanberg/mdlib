#ifndef _MD_XTC_H_
#define _MD_XTC_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_trajectory_i;

bool md_xtc_trajectory_open(struct md_trajectory_i* traj, str_t filename, struct md_allocator_i* alloc);
bool md_xtc_trajectory_close(struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif

#endif
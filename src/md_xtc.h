#ifndef _MD_XTC_H_
#define _MD_XTC_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_trajectory_i;

struct md_trajectory_i* md_xtc_trajectory_open(const char* filename_str, int64_t filename_len, struct md_allocator_i* alloc);
void md_xtc_trajectory_close(struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif

#endif
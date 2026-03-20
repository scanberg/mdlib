#pragma once

#include <core/md_str.h>

struct md_system_t;
struct md_trajectory_i;
struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

bool md_trr_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags);

struct md_trajectory_i* md_trr_trajectory_create(str_t filename, struct md_allocator_i* alloc, uint32_t flags);

#ifdef __cplusplus
}
#endif

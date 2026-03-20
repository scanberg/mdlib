#pragma once

#include <core/md_str.h>

struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

bool md_trr_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags);

#ifdef __cplusplus
}
#endif

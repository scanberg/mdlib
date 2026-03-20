#pragma once

#include <core/md_str.h>
#include <stdbool.h>

struct md_system_t;

// Utils for reading PDBX/mmCIF files
// https://mmcif.wwpdb.org/
// This is only meant to extract subset of data from the file, not the entire file and is thus
// not a complete implementation of the format.

#ifdef __cplusplus
extern "C" {
#endif

bool md_mmcif_system_init_from_file(struct md_system_t* sys, str_t filename);
bool md_mmcif_system_init_from_str (struct md_system_t* sys, str_t str);

#ifdef __cplusplus
}
#endif

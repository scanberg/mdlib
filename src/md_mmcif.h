#pragma once

#include <core/md_str.h>
#include <stdbool.h>

struct md_system_t;
struct md_executor_t;
// Coordinates are handed back through a state; only passed by pointer here.
typedef struct md_system_state_t md_system_state_t;


// Utils for reading PDBX/mmCIF files
// https://mmcif.wwpdb.org/
// This is only meant to extract subset of data from the file, not the entire file and is thus
// not a complete implementation of the format.

#ifdef __cplusplus
extern "C" {
#endif

// The rows of _atom_site are read over the threads of ex (md_executor.h), NULL for the calling thread alone. The system is
// the same however many there are.
bool md_mmcif_system_init_from_file(struct md_system_t* sys, md_system_state_t* state, str_t filename, const struct md_executor_t* ex);
bool md_mmcif_system_init_from_str (struct md_system_t* sys, md_system_state_t* state, str_t str, const struct md_executor_t* ex);

#ifdef __cplusplus
}
#endif

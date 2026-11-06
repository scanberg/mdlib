#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_system_t;
// Coordinates are handed back through a state; only passed by pointer here.
typedef struct md_system_state_t md_system_state_t;

struct md_system_loader_i;
struct md_executor_t;

// A GRO file is read into a system, and nothing else: there is no representation of the file of its own. The system: the
// atoms (of the first frame), where they are, their velocities where the file has them (atom/velocity, in nm/ps), their
// types, the residues and the cell. No bonds, nor what the residues are by them (amino acids, nucleotides, waters by
// their atoms): md_util_system_infer derives those, when the caller wants them. The executor (NULL: the calling thread)
// reads a file and its atom lines in parallel.
bool md_gro_system_init_from_file(struct md_system_t* sys, md_system_state_t* state, str_t filename, const struct md_executor_t* executor);
bool md_gro_system_init_from_str (struct md_system_t* sys, md_system_state_t* state, str_t str, const struct md_executor_t* executor);

#ifdef __cplusplus
}
#endif


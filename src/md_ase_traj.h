#pragma once

#include <core/md_str.h>
#include <stdbool.h>
#include <stdint.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

// Loads a modern ASE ULM trajectory with a fixed set of atoms.
bool md_ase_traj_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, str_t filename);
bool md_ase_traj_system_publish_run(struct md_system_t* sys, str_t filename, str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

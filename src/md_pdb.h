#pragma once

#include <core/md_str.h>
#include <core/md_file.h>

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_system_t;
// Coordinates are handed back through a state; only passed by pointer here.
typedef struct md_system_state_t md_system_state_t;

struct md_allocator_i;
struct md_executor_t;

enum {
    MD_PDB_OPTION_NONE                      = 0,
    MD_PDB_OPTION_CONCAT_MODELS             = 1,
    MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE  = 2, // Unused: a run is published with md_pdb_system_publish_run, which takes its own flags
};

typedef uint32_t md_pdb_options_t;

// A PDB file is read into a system, and nothing else: there is no representation of the file of its own. What the
// system has of it is the atoms (ATOM, HETATM) of the first model, the chains (and where TER records end them), the
// cell of the first CRYST1 record, the bonds of the CONECT records and the biomolecules of REMARK 350. The other
// records are passed over: this is a loader of systems, not a reader of every record of the format.

// SYSTEM
// The system: the atoms (of the first model, or of all with MD_PDB_OPTION_CONCAT_MODELS), where they are, their types,
// the residues, the units of the file (its chains, by the character of column 22 as the file has it: a digit is an id as
// a letter is, a blank none), the bonds of the CONECT records, the cell and the assemblies. No bonds beyond the file's,
// nor what the residues are by them: md_system_infer derives those, when the caller wants them. The executor
// (NULL: the calling thread) reads a file and its lines in parallel; of a file of several models, the first model and
// the CONECT records after the last are read, not the models between.
bool md_pdb_system_init_from_file(struct md_system_t* sys, md_system_state_t* state, md_str_t filename, md_pdb_options_t options, const struct md_executor_t* executor);
bool md_pdb_system_init_from_str (struct md_system_t* sys, md_system_state_t* state, md_str_t str,      md_pdb_options_t options, const struct md_executor_t* executor);

// EXPORT
bool md_pdb_system_write_state_to_file(md_file_t file, const struct md_system_t* sys, const md_system_state_t* state, const int32_t* atom_indices, size_t num_atoms, int model_num);

// RUN
// A file of several models of the same atoms, published as a run (see md_run_publish in md_system.h): one frame per
// model, time as model ordinals without a unit, the file's CRYST1 cell at every frame (zero without one) and
// atom/position read from the model's records when a frame is asked for: the atoms the system has of it (its first
// alternate locations). The index of where each model starts is found by a scan of the record names, and kept beside
// the file ('<file>.cache', written unless flags carries MD_RUN_FLAG_DISABLE_CACHE_WRITE). False for a file of one
// model.
bool md_pdb_system_publish_run(struct md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags);

#ifdef __cplusplus
}
#endif

#pragma once

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;
// Coordinates are handed back through a state; only passed by pointer here.
typedef struct md_system_state_t md_system_state_t;

struct md_system_loader_i;

typedef struct md_gro_atom_t {
    int32_t res_id;
    float x;
    float y;
    float z;
    // GROMACS writes three optional velocity columns after the coordinates. NAN when the line has
    // none, which is the common case for a plain structure file - see md_attributes_publish_atom_column
    // on why absence is a value here rather than a zero.
    float vx;
    float vy;
    float vz;
    char res_name[8];
    char atom_name[8];
} md_gro_atom_t;

typedef struct md_gro_data_t {
    char title[256];
    float box[3][3];

    // Field data
    size_t num_atoms;
    md_gro_atom_t* atom_data;
} md_gro_data_t;

// Parse a text-blob as GRO
bool md_gro_data_parse_str(md_gro_data_t* data, str_t string, struct md_allocator_i* alloc);
bool md_gro_data_parse_file(md_gro_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_gro_data_free(md_gro_data_t* data, struct md_allocator_i* alloc);

// Molecule
bool md_gro_system_init_from_data(struct md_system_t* sys, md_system_state_t* state, const md_gro_data_t* gro_data);
bool md_gro_system_init_from_file(struct md_system_t* sys, md_system_state_t* state, str_t filename);
bool md_gro_system_init_from_str (struct md_system_t* sys, md_system_state_t* state, str_t str);

#ifdef __cplusplus
}
#endif


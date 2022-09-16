#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_molecule_t;
struct md_molecule_api;

typedef struct md_gro_atom_t {
    int32_t res_id;
    char res_name[8];
    char atom_name[8];
    float x;
    float y;
    float z;
} md_gro_atom_t;

typedef struct md_gro_data_t {
    char title[256];
    float cell_ext[3];

    // Field data
    int64_t num_atoms;
    md_gro_atom_t* atom_data;
} md_gro_data_t;

// Parse a text-blob as GRO
bool md_gro_data_parse_str(md_gro_data_t* data, str_t string, struct md_allocator_i* alloc);
bool md_gro_data_parse_file(md_gro_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_gro_data_free(md_gro_data_t* data, struct md_allocator_i* alloc);

// Molecule
bool md_gro_molecule_init(struct md_molecule_t* mol, const md_gro_data_t* gro_data, struct md_allocator_i* alloc);

struct md_molecule_api* md_gro_molecule_api();

#ifdef __cplusplus
}
#endif


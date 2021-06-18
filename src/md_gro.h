#ifndef _MD_GRO_H_
#define _MD_GRO_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator;
struct md_molecule;

typedef struct md_gro_atom {
    int32_t res_id;
    char res_name[8];
    char atom_name[8];
    float x;
    float y;
    float z;
} md_gro_atom_t;

typedef struct md_gro_data {
    char title[256];
    float unit_cell[3];

    // Field data
    int64_t num_atoms;
    md_gro_atom_t* atom_data;
} md_gro_data_t;

// Parse a text-blob as GRO
bool md_gro_data_parse_str(str_t string, struct md_gro_data* data, struct md_allocator* alloc);
bool md_gro_data_parse_file(str_t filename, struct md_gro_data* data, struct md_allocator* alloc);

void md_gro_data_free(struct md_gro_data* data, struct md_allocator* alloc);

// Molecule
bool md_gro_molecule_init(struct md_molecule* mol, const struct md_gro_data* gro_data, struct md_allocator* alloc);
bool md_gro_molecule_free(struct md_molecule* mol, struct md_allocator* alloc);

#ifdef __cplusplus
}
#endif

#endif
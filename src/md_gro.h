#ifndef _MD_GRO_H_
#define _MD_GRO_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <md_molecule.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

typedef struct md_gro_atom_t md_gro_atom_t;

struct md_gro_atom_t {
    int32_t res_id;
    char res_name[8];
    char atom_name[8];
    float x;
    float y;
    float z;
};

struct md_gro_data_t {
    char title[256];
    float unit_cell[3];

    // Field data
    int64_t num_atoms;
    md_gro_atom_t* atom_data;
};

typedef struct md_gro_label_t {
    char str[8];
} md_gro_label_t;

typedef struct md_gro_molecule_t {
    struct md_molecule mol;

    struct {
        md_gro_label_t* labels;
        struct md_allocator_i* arena;
    } storage;
} md_gro_molecule_t;

// Parse a text-blob as GRO
bool md_gro_data_parse_str(str_t string, struct md_gro_data_t* data, struct md_allocator_i* alloc);
bool md_gro_data_parse_file(str_t filename, struct md_gro_data_t* data, struct md_allocator_i* alloc);

void md_gro_data_free(struct md_gro_data_t* data, struct md_allocator_i* alloc);

// Molecule
bool md_gro_molecule_init(struct md_gro_molecule_t* mol, const struct md_gro_data_t* gro_data, struct md_allocator_i* alloc);
bool md_gro_molecule_free(struct md_gro_molecule_t* mol, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

#endif
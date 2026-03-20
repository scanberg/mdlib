#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_system_t;

enum {
    MD_XYZ_OPTION_NONE = 0,
    MD_XYZ_OPTION_DISABLE_CACHE_WRITE = 1, // Only applies if the XYZ stems from file on disk and contains a trajectory.
};

typedef uint32_t md_xyz_options_t;

typedef struct md_xyz_coordinate_t {
	int  atom_index;
	char element_symbol[4];
	int  atomic_number;
	float x;
	float y;
	float z;
	int atom_type;			// Tinker only
	int connectivity[9];	// Tinker only
} md_xyz_coordinate_t;

typedef struct md_xyz_model_t {
	size_t byte_offset;
	uint32_t beg_coord_index;
	uint32_t end_coord_index;
	float cell[3][3];
	str_t comment;
} md_xyz_model_t;

typedef struct md_xyz_data_t {
	size_t num_coordinates;
	md_xyz_coordinate_t* coordinates;

	size_t num_models;
	md_xyz_model_t* models;
} md_xyz_data_t;

// RAW FUNCTIONS
// Parse a text-blob as XYZ
bool md_xyz_data_parse_str(md_xyz_data_t* data, str_t str, struct md_allocator_i* alloc);
bool md_xyz_data_parse_file(md_xyz_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc);

// SYSTEM
bool md_xyz_system_init_from_data(struct md_system_t* sys, const md_xyz_data_t* data, md_xyz_options_t options);
bool md_xyz_system_init_from_file(struct md_system_t* sys, str_t filename, md_xyz_options_t options);
bool md_xyz_system_init_from_str (struct md_system_t* sys, str_t str, md_xyz_options_t options);

#ifdef __cplusplus
}
#endif

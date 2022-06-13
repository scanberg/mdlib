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
struct md_trajectory_i;
struct md_trajectory_api;

typedef struct md_xyz_coordinate_t {
	int atom_index;
	char element_symbol[4];
	int  atomic_number;
	float x;
	float y;
	float z;
	int atom_type;			// Tinker only
	int connectivity[9];	// Tinker only
} md_xyz_coordinate_t;

typedef struct md_xyz_model_t {
	int64_t byte_offset;
	int32_t beg_coord_index;
	int32_t end_coord_index;
	float cell_extent[3];	// Tinker (ARC) only
	float cell_angle[3];	// Tinker (ARC) only
	char comment[88];
} md_xyz_model_t;

typedef struct md_xyz_data_t {
	int64_t num_coordinates;
	md_xyz_coordinate_t* coordinates;
	int64_t num_models;
	md_xyz_model_t* models;
} md_xyz_data_t;

// RAW FUNCTIONS
// Parse a text-blob as XYZ
bool md_xyz_data_parse_str(md_xyz_data_t* data, str_t str, struct md_allocator_i* alloc);
bool md_xyz_data_parse_file(md_xyz_data_t* data, str_t filename, struct md_allocator_i* alloc);
void md_xyz_data_free(md_xyz_data_t* data, struct md_allocator_i* alloc);

// MOLECULE
bool md_xyz_molecule_init(struct md_molecule_t* mol, const md_xyz_data_t* data, struct md_allocator_i* alloc);
bool md_xyz_molecule_free(struct md_molecule_t* mol, struct md_allocator_i* alloc);

struct md_molecule_api* md_xyz_molecule_api();

// TRAJECTORY
struct md_trajectory_i* md_xyz_trajectory_create(str_t filename, struct md_allocator_i* alloc);
void md_xyz_trajectory_free(struct md_trajectory_i* traj);

struct md_trajectory_api* md_xyz_trajectory_api();

#ifdef __cplusplus
}
#endif

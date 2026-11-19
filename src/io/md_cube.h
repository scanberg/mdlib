#pragma once

// Utils for (Gaussian) Cube files
// https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

#include <core/md_str.h>

struct md_allocator_i;
struct md_writer_i;
struct md_system_t;
struct md_system_state_t;
struct md_load_desc_t;

typedef float md_cube_v3[3];

typedef struct md_cube_t {
    // Gaussian files supports two lines for comments. with maximum length of 80 characters per line
    md_str_t title;
	md_str_t comment;
	
	md_cube_v3 origin;
	md_cube_v3 xaxis;
	md_cube_v3 yaxis;
	md_cube_v3 zaxis;

	struct {
        int num_x;	// Number of points along x-axis
        int num_y;	// Number of points along y-axis
        int num_z;	// Number of points along z-axis
		int num_m;	// Number of data values per voxel
        int*   id;  // Data ID [Optional] If defined, it should have num_m elements
		float* val; // Data is stored as data[x][y][z][m]
	} data;

	struct {
        int			count;	// Number of atoms
		int*		number;	// Atomic number
		float*		charge; // Nuclear charge
        md_cube_v3* coord;  // Coordinates
	} atom;
} md_cube_t;

#ifdef __cplusplus
extern "C" {
#endif

// A cube file, from the text of one or from a file. Writing a cube is the application's (viamd's density volume export)
bool md_cube_deserialize(md_cube_t* cube, md_str_t data, struct md_allocator_i* alloc);
bool md_cube_file_load(md_cube_t* cube, md_str_t path_to_file, struct md_allocator_i* alloc);

// A cube as a system of its atoms (system_init), or added to a loaded system whose atoms, or whose QM atoms
// (qm/atom/system_index), its atoms are in order (supplement; md_cube_system_is_file_supplemental tells, from the header
// alone). Either publishes the grid under cube/<stem of the file> (see md_cube.c, cube_publish): its values, where its
// voxels are in Å as a script volume's (origin, basis), and its atoms in its frame, by which it is placed.
bool md_cube_system_init_from_file(struct md_system_t* sys, struct md_system_state_t* state, md_str_t path, const struct md_load_desc_t* desc);
bool md_cube_system_supplement_from_file(struct md_system_t* sys, md_str_t path, const struct md_load_desc_t* desc);
bool md_cube_system_is_file_supplemental(const struct md_system_t* sys, md_str_t path);

// Only free using this function if the cube object was allocated through cube_file_load or cube_file_deserialize
void md_cube_free(md_cube_t* cube, struct md_allocator_i* alloc);


#ifdef __cplusplus
}
#endif

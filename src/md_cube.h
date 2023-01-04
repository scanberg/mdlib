#pragma once

// Utils for (Gaussian) Cube files
// https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

#include <core/md_str.h>

struct md_allocator_i;

typedef float md_cube_v3[3];

typedef struct md_cube_t {
    // Gaussian files supports two lines for comments. with maximum length of 80 characters per line
    str_t title;
	str_t comment;
	
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

str_t md_cube_serialize(const md_cube_t* cube, struct md_allocator_i* alloc);
bool  md_cube_deserialize(md_cube_t* cube, str_t data, struct md_allocator_i* alloc);
	
bool md_cube_file_store(const md_cube_t* cube, str_t path_to_file);
bool md_cube_file_load(md_cube_t* cube, str_t path_to_file, struct md_allocator_i* alloc);

bool md_cube_valid(const md_cube_t* cube);

// Only free using this function if the cube object was allocated through cube_file_load or cube_file_deserialize
void md_cube_free(md_cube_t* cube, struct md_allocator_i* alloc);


#ifdef __cplusplus
}
#endif

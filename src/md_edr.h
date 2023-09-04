#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_unit.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_edr_energy_t {
	str_t name;
	str_t unit_str;
    md_unit_t unit;
    float* values;  // Has length num_frames
} md_edr_energy_t;

typedef struct md_edr_energies_t {
	int64_t			 num_frames;
	double*			 frame_time;	// timestamp of frames

	int64_t			 num_energies;
	md_edr_energy_t* energy;
	
	struct md_allocator_i* alloc;
} md_edr_energies_t;

bool md_edr_energies_parse_file(md_edr_energies_t* energies, str_t filename, struct md_allocator_i* alloc);
void md_edr_energies_free(md_edr_energies_t* energies);

#if 0
// Represents the string block inside the edr file which gives us the names and units (encoded as text)
// of all energy fields present in the file.
typedef struct md_edr_str_block_t {
	int64_t count;
	str_t* names;
	str_t* units;
} md_edr_str_block_t;

bool md_edr_str_block_read_file(md_edr_str_block_t* str_block, str_t filename, struct md_allocator_i* alloc);
void md_edr_str_block_free(md_edr_str_block_t* str_block, struct md_allocator_i* alloc);
#endif


#ifdef __cplusplus
}
#endif

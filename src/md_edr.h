#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_unit.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_edr_energies_t {
	int64_t		count;
	str_t*		names;
	md_unit_t*	units;
	int64_t		num_values;
	double**	dvalues;
	float**		fvalues;
	struct md_allocator_i* alloc;
} md_edr_energies_t;

bool md_edr_energies_read_file(md_edr_energies_t* energies, str_t filename, struct md_allocator_i* alloc);
void md_edr_energies_free(md_edr_energies_t* energies);

// Represents the string block inside the edr file which gives us the names and units (encoded as text)
// of all energy fields present in the file.
typedef struct md_edr_str_block_t {
	int64_t count;
	str_t* names;
	str_t* units;
} md_edr_str_block_t;

bool md_edr_str_block_read_file(md_edr_str_block_t* str_block, str_t filename, struct md_allocator_i* alloc);
void md_edr_str_block_free(md_edr_str_block_t* str_block, struct md_allocator_i* alloc);


#ifdef __cplusplus
}
#endif

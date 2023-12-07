#include "ubench.h"

#include <md_lammps.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UBENCH_EX(lammps, load) {
	md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
	str_t path = STR(MD_BENCHMARK_DATA_DIR "/Water_Ethane_Cubic_Init.data");

	md_file_o* file = md_file_open(path, MD_FILE_READ);
	if (!file) {
		MD_LOG_ERROR("Could not open file '%.*s'", path.len, path.ptr);
		return;
	}
	UBENCH_SET_BYTES(md_file_size(file));
	md_file_close(file);

	UBENCH_DO_BENCHMARK() {
		md_arena_allocator_reset(alloc);
		md_lammps_data_t lammps = {0};
		md_lammps_data_parse_file(&lammps, path, MD_LAMMPS_ATOM_FORMAT_FULL, alloc);
	}

	md_arena_allocator_destroy(alloc);
}

#include <md_cube.h>

#include "core/md_allocator.h"
#include "core/md_parse.h"
#include "core/md_file.h"
#include "core/md_log.h"
#include "core/md_writer.h"
#include "core/md_vec_math.h"
#include "core/md_array.h"

#include <string.h>

static bool md_cube_valid(const md_cube_t* cube) {
	if (!cube) return false;

	md_float3_t x = md_float3_set(cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]);
	md_float3_t y = md_float3_set(cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]);
	md_float3_t z = md_float3_set(cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]);

	if (fabsf(md_float3_dot(x, md_float3_cross(y, z))) < 1.e-3f) {
		MD_LOG_ERROR("CUBE: The coordinate axes do not span a volume");
		return false;
	}

	if (cube->data.num_x <= 0 || cube->data.num_y <= 0 || cube->data.num_z <= 0 || cube->data.num_m <= 0) {
		MD_LOG_ERROR("CUBE: One or more data dimensions are zero");
		return false;
	}

	if (cube->atom.count < 0) {
		MD_LOG_ERROR("CUBE: Geom count is negative");
		return false;
	}

	if (cube->atom.count == 0) {
		// A grid without atoms is still a grid
		return true;
	}

	if (!cube->atom.coord) {
		MD_LOG_ERROR("CUBE: coord was NULL");
		return false;
	}

	if (!cube->atom.number) {
		MD_LOG_ERROR("CUBE: number was NULL");
		return false;
	}

	if (!cube->atom.charge) {
		MD_LOG_ERROR("CUBE: charge was NULL");
		return false;
	}

	return true;
}

md_str_t tidy_comment(md_str_t comment) {
    comment = md_str_substr(comment, 0, md_str_find_char(comment, '\n'));
    return comment;
}

static inline bool str_extract_v3(md_cube_v3 v, md_str_t* str) {
	return md_parse_next_f32(&v[0], str) && md_parse_next_f32(&v[1], str) && md_parse_next_f32(&v[2], str);
}

bool md_cube_deserialize(md_cube_t* cube, md_str_t str, struct md_allocator_i* alloc) {
	MD_ASSERT(cube);
	MD_ASSERT(alloc);

	md_str_t title, comment;
    if (!(md_parse_next_line(&title,	 &str) &&
		  md_parse_next_line(&comment, &str)))
	{
		MD_LOG_ERROR("CUBE: Could not extract comment section");
		return false;
	}
    cube->title   = md_alloc_str(alloc, title);
    cube->comment = md_alloc_str(alloc, comment);

	md_str_t line;
	
	// natoms, origin, [nval]
	if (!md_parse_next_line(&line, &str)) {
		goto line_error;
	}
    if (!(md_parse_next_i32(&cube->atom.count, &line) && str_extract_v3(cube->origin, &line))) {
		goto parse_error;
	}
	if (!md_parse_next_i32(&cube->data.num_m, &line)) {
		cube->data.num_m = 1;
	}

	// num_x, xaxis
	if (!md_parse_next_line(&line, &str)) {
		goto line_error;
	}
    if (!(md_parse_next_i32(&cube->data.num_x, &line) && str_extract_v3(cube->xaxis, &line))) {
		goto parse_error;
	}
	
	// num_y, yaxis
	if (!md_parse_next_line(&line, &str)) {
		goto line_error;
	}
	if (!(md_parse_next_i32(&cube->data.num_y, &line) && str_extract_v3(cube->yaxis, &line))) {
		goto parse_error;
	}

	// num_z, zaxis
	if (!md_parse_next_line(&line, &str)) {
		goto line_error;
	}
	if (!(md_parse_next_i32(&cube->data.num_z, &line) && str_extract_v3(cube->zaxis, &line))) {
		goto parse_error;
	}

	// A negative number of points (by convention that of the first axis) means that the lengths of the header
	// - the origin, the axes and the atoms - are given in Ångström rather than Bohr. They are kept in Bohr here.
	const bool angstrom = cube->data.num_x < 0;
	cube->data.num_x = abs(cube->data.num_x);
	cube->data.num_y = abs(cube->data.num_y);
	cube->data.num_z = abs(cube->data.num_z);
	const float length_scale = angstrom ? (float)(1.0 / 0.52917721092) : 1.0f;
	for (int i = 0; i < 3; ++i) {
		cube->origin[i] *= length_scale;
		cube->xaxis[i]  *= length_scale;
		cube->yaxis[i]  *= length_scale;
		cube->zaxis[i]  *= length_scale;
	}

	// geom
	const int count = abs(cube->atom.count);
	for (int i = 0; i < count; ++i) {
        int number;
		float charge;
        md_cube_v3 coord;
		
        if (!md_parse_next_line(&line, &str)) {
			goto line_error;
		}
        if (!(md_parse_next_i32(&number, &line) && md_parse_next_f32(&charge, &line) && str_extract_v3(coord, &line))) {
			goto parse_error;
		}
		coord[0] *= length_scale;
		coord[1] *= length_scale;
		coord[2] *= length_scale;
		
		md_array_push(cube->atom.number, number, alloc);
		md_array_push(cube->atom.charge, charge, alloc);
        md_array_push_array(cube->atom.coord, coord, 1, alloc);
	}

	if (cube->atom.count < 0) {
        cube->atom.count = -cube->atom.count;
        // data id
		if (!md_parse_next_i32(&cube->data.num_m, &str)) {
			goto parse_error;
		}
		for (int i = 0; i < cube->data.num_m; ++i) {
			int id;
			if (!md_parse_next_i32(&id, &str)) {
				goto parse_error;
			}
			md_array_push(cube->data.id, id, alloc);
		}
	}

	if (cube->data.num_x <= 0 || cube->data.num_y <= 0 || cube->data.num_z <= 0 || cube->data.num_m <= 0) {
		MD_LOG_ERROR("CUBE: One or more data dimensions are zero");
		return false;
	}

	// In 64 bits: a 1300^3 grid alone does not fit in an int
	const int64_t num_data = (int64_t)cube->data.num_x * cube->data.num_y * cube->data.num_z * cube->data.num_m;
	if ((uint64_t)num_data > str.len / 2 + 1) {
		// Every value takes at least a digit and a separator: the file cannot hold them
		MD_LOG_ERROR("CUBE: The file is too short for the %lld values its grid declares", (long long)num_data);
		return false;
	}
	md_array_ensure(cube->data.val, (size_t)num_data, alloc);
    for (int64_t i = 0; i < num_data; ++i) {
        float val;

        if (!md_parse_next_f32(&val, &str)) {
			goto parse_error;
		}
        md_array_push(cube->data.val, val, alloc);
    }

	return md_cube_valid(cube);

line_error:
	MD_LOG_ERROR("CUBE: Could not extract line");
	return false;
parse_error:
	MD_LOG_ERROR("CUBE: Could not parse expected value");
	return false;
}

bool md_cube_file_load(md_cube_t* cube, md_str_t path, md_allocator_i* alloc) {
	if (!cube) {
		MD_LOG_ERROR("CUBE: cube-object was NULL");
		return false;
	}

	bool success = false;
	md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
	md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    md_str_t str = md_load_textfile(path, temp_alloc);
	if (!md_str_empty(str)) {
        success = md_cube_deserialize(cube, str, alloc);
	}
	md_temp_end(temp_scope);
	
	return success;
}

void md_cube_free(md_cube_t* cube, md_allocator_i* alloc) {
	MD_ASSERT(cube);
	MD_ASSERT(alloc);
	if (cube->data.val)		md_array_free(cube->data.val,	 alloc);
	if (cube->data.id)		md_array_free(cube->data.id,	 alloc);
	if (cube->atom.coord)	md_array_free(cube->atom.coord,	 alloc); 
	if (cube->atom.number)	md_array_free(cube->atom.number, alloc); 
	if (cube->atom.charge)	md_array_free(cube->atom.charge, alloc); 
	if (!md_str_empty(cube->title)) md_free_str(alloc, cube->title);
	if (!md_str_empty(cube->comment)) md_free_str(alloc, cube->comment);
}

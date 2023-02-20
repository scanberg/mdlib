#include "md_cube.h"
#include "core/md_allocator.h"
#include "core/md_os.h"
#include "core/md_log.h"
#include "core/md_str_builder.h"
#include "core/md_vec_math.h"
#include "core/md_array.h"

#include <string.h>

bool md_cube_valid(const md_cube_t* cube) {
	if (!cube) return false;

	vec3_t x = vec3_set(cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]);
	vec3_t y = vec3_set(cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]);
	vec3_t z = vec3_set(cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]);

	if (fabsf(vec3_dot(x, vec3_cross(y, z))) < 1.e-3f) {
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

static bool open_file(md_file_o* file, str_t path) {
	ASSERT(file);

	if (str_empty(path)) {
		MD_LOG_ERROR("CUBE: Path was empty");
		return false;
	}

	file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("CUBE: Could not open file: '%.*s'", path.len, path.ptr);
		return false;
	}

	return true;
}

str_t tidy_comment(str_t comment) {
    int64_t new_line = str_find_char(comment, '\n');
	if (new_line != -1) {
		comment.len = new_line;
	}
    return comment;
}

str_t md_cube_serialize(const md_cube_t* cube, struct md_allocator_i* alloc) {
    ASSERT(alloc);
	str_t str = {0};

	if (!md_cube_valid(cube)) {
		MD_LOG_ERROR("CUBE: cube-object was invalid");
		goto done;
	}

	md_strb_t sb = {0};
	md_strb_init(&sb, default_allocator);
	
	md_strb_push_str(&sb, tidy_comment(cube->title));   md_strb_push_char(&sb, '\n');
	md_strb_push_str(&sb, tidy_comment(cube->comment)); md_strb_push_char(&sb, '\n');

	int atom_count = cube->data.id ? -cube->atom.count : cube->atom.count;
	int num_val = cube->data.id ? 1 : cube->data.num_m;
	md_strb_fmt(&sb, "%5d%12.6f%12.6f%12.6f%5d\n", atom_count, cube->origin[0], cube->origin[1], cube->origin[2], num_val);
    md_strb_fmt(&sb, "%5d%12.6f%12.6f%12.6f\n", cube->data.num_x, cube->xaxis[0], cube->xaxis[1], cube->xaxis[2]);
	md_strb_fmt(&sb, "%5d%12.6f%12.6f%12.6f\n", cube->data.num_y, cube->yaxis[0], cube->yaxis[1], cube->yaxis[2]);
	md_strb_fmt(&sb, "%5d%12.6f%12.6f%12.6f\n", cube->data.num_z, cube->zaxis[0], cube->zaxis[1], cube->zaxis[2]);
	
	for (int i = 0; i < cube->atom.count; ++i) {
		md_strb_fmt(&sb, "%5d%12.6f%12.6f%12.6f%12.6f\n", cube->atom.number[i], cube->atom.charge[i], cube->atom.coord[i][0], cube->atom.coord[i][1], cube->atom.coord[i][2]);
	}

	if (cube->data.id) {
		md_strb_fmt(&sb, "%5d", cube->data.num_m);
		for (int i = 0; i < cube->data.num_m; ++i) {
			md_strb_fmt(&sb, "%5d", cube->data.id[i]);
			if (i % 23 == 22) {
				md_strb_push_char(&sb, '\n');
			}
		}
		md_strb_push_char(&sb, '\n');
	}
	
	const int64_t data_count = cube->data.num_x * cube->data.num_y * cube->data.num_z;
	for (int64_t i = 0; i < data_count; ++i) {
        md_strb_fmt(&sb, " %12.5E", cube->data.val[i]);
        if (i % 6 == 5) {
			md_strb_push_char(&sb, '\n');
		}
	}
	
	str = str_copy(md_strb_to_str(&sb), alloc);
	md_strb_free(&sb);
done:
	return str;
}

static inline bool str_extract_v3(md_cube_v3 v, str_t* str) {
	return str_extract_f32(&v[0], str) && str_extract_f32(&v[1], str) && str_extract_f32(&v[2], str);
}

bool md_cube_deserialize(md_cube_t* cube, str_t str, struct md_allocator_i* alloc) {
	ASSERT(cube);
	ASSERT(alloc);

	str_t title, comment;
    if (!(str_extract_line(&title,	 &str) &&
		  str_extract_line(&comment, &str)))
	{
		MD_LOG_ERROR("CUBE: Could not extract comment section");
		return false;
	}
    cube->title   = str_copy(title, alloc);
    cube->comment = str_copy(comment, alloc);

	str_t line;
	
	// natoms, origin, [nval]
	if (!str_extract_line(&line, &str)) {
		goto line_error;
	}
    if (!(str_extract_i32(&cube->atom.count, &line) && str_extract_v3(cube->origin, &line))) {
		goto parse_error;
	}
	if (!str_extract_i32(&cube->data.num_m, &line)) {
		cube->data.num_m = 1;
	}

	// num_x, xaxis
	if (!str_extract_line(&line, &str)) {
		goto line_error;
	}
    if (!(str_extract_i32(&cube->data.num_x, &line) && str_extract_v3(cube->xaxis, &line))) {
		goto parse_error;
	}
	
	// num_y, yaxis
	if (!str_extract_line(&line, &str)) {
		goto line_error;
	}
	if (!(str_extract_i32(&cube->data.num_y, &line) && str_extract_v3(cube->yaxis, &line))) {
		goto parse_error;
	}

	// num_z, zaxis
	if (!str_extract_line(&line, &str)) {
		goto line_error;
	}
	if (!(str_extract_i32(&cube->data.num_z, &line) && str_extract_v3(cube->zaxis, &line))) {
		goto parse_error;
	}

	// geom
	const int count = abs(cube->atom.count);
	for (int i = 0; i < count; ++i) {
        int number;
		float charge;
        md_cube_v3 coord;
		
        if (!str_extract_line(&line, &str)) {
			goto line_error;
		}
        if (!(str_extract_i32(&number, &line) && str_extract_f32(&charge, &line) && str_extract_v3(coord, &line))) {
			goto parse_error;
		}
		
		md_array_push(cube->atom.number, number, alloc);
		md_array_push(cube->atom.charge, charge, alloc);
        md_array_push_array(cube->atom.coord, coord, 1, alloc);
	}

	if (cube->atom.count < 0) {
        cube->atom.count = -cube->atom.count;
        // data id
		if (!str_extract_i32(&cube->data.num_m, &str)) {
			goto parse_error;
		}
		for (int i = 0; i < cube->data.num_m; ++i) {
			int id;
			if (!str_extract_i32(&id, &str)) {
				goto parse_error;
			}
			md_array_push(cube->data.id, id, alloc);
		}
	}

	int64_t num_data = cube->data.num_x * cube->data.num_y * cube->data.num_z * cube->data.num_m;
    for (int64_t i = 0; i < num_data; ++i) {
        float val;
        if (!str_extract_f32(&val, &str)) {
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

bool md_cube_file_load(md_cube_t* cube, str_t path, md_allocator_i* alloc) {
	if (!cube) {
		MD_LOG_ERROR("CUBE: cube-object was NULL");
		return false;
	}

	bool success = false;
    str_t str = load_textfile(path, default_allocator);
	if (!str_empty(str)) {
        success = md_cube_deserialize(cube, str, alloc);
		str_free(str, default_allocator);
	}
	
	return success;
}

bool md_cube_file_store(const md_cube_t* cube, str_t path) {
	if (!md_cube_valid(cube)) {
		MD_LOG_ERROR("CUBE: cube-object was invalid");
		return false;
	}

	md_file_o* file = NULL;
	bool success = false;
	if (open_file(file, path)) {
		str_t str = md_cube_serialize(cube, default_allocator);
        if (!str_empty(str)) {
			success = md_file_write(file, str_ptr(str), str_len(str)) == str_len(str);
			str_free(str, default_allocator);
        }
		md_file_close(file);
	}
	
	return success;
}

void md_cube_free(md_cube_t* cube, md_allocator_i* alloc) {
	ASSERT(cube);
	ASSERT(alloc);
	if (cube->data.val)		md_array_free(cube->data.val,	 alloc);
	if (cube->atom.coord)	md_array_free(cube->atom.coord,	 alloc); 
	if (cube->atom.number)	md_array_free(cube->atom.number, alloc); 
	if (cube->atom.charge)	md_array_free(cube->atom.charge, alloc); 
	if (!str_empty(cube->title)) str_free(cube->title, alloc);
	if (!str_empty(cube->comment)) str_free(cube->comment, alloc);
}
#include <md_cube.h>
#include <md_io.h>
#include <md_system.h>
#include <md_qm.h>
#include <md_element.h>

#include "core/md_allocator.h"
#include "core/md_arena_allocator.h"
#include "core/md_path.h"
#include "core/md_unit.h"
#include "core/md_format.h"
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

// The header: the comments, the grid's axes, and the atoms, in bohr. Leaves str at what follows the atoms (the ids of an
// MO cube, or the values).
static bool cube_parse_header(md_cube_t* cube, md_str_t* str, struct md_allocator_i* alloc) {
	md_str_t title, comment;
    if (!(md_parse_next_line(&title,	 str) &&
		  md_parse_next_line(&comment, str)))
	{
		MD_LOG_ERROR("CUBE: Could not extract comment section");
		return false;
	}
    cube->title   = md_alloc_str(alloc, title);
    cube->comment = md_alloc_str(alloc, comment);

	md_str_t line;
	
	// natoms, origin, [nval]
	if (!md_parse_next_line(&line, str)) {
		goto line_error;
	}
    if (!(md_parse_next_i32(&cube->atom.count, &line) && str_extract_v3(cube->origin, &line))) {
		goto parse_error;
	}
	if (!md_parse_next_i32(&cube->data.num_m, &line)) {
		cube->data.num_m = 1;
	}

	// num_x, xaxis
	if (!md_parse_next_line(&line, str)) {
		goto line_error;
	}
    if (!(md_parse_next_i32(&cube->data.num_x, &line) && str_extract_v3(cube->xaxis, &line))) {
		goto parse_error;
	}
	
	// num_y, yaxis
	if (!md_parse_next_line(&line, str)) {
		goto line_error;
	}
	if (!(md_parse_next_i32(&cube->data.num_y, &line) && str_extract_v3(cube->yaxis, &line))) {
		goto parse_error;
	}

	// num_z, zaxis
	if (!md_parse_next_line(&line, str)) {
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
		
        if (!md_parse_next_line(&line, str)) {
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
	return true;

line_error:
	MD_LOG_ERROR("CUBE: Could not extract line");
	return false;
parse_error:
	MD_LOG_ERROR("CUBE: Could not parse expected value");
	return false;
}

bool md_cube_deserialize(md_cube_t* cube, md_str_t str, struct md_allocator_i* alloc) {
	MD_ASSERT(cube);
	MD_ASSERT(alloc);

	if (!cube_parse_header(cube, &str, alloc)) {
		return false;
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

// ### THE CUBE AS A SYSTEM ###

#define CUBE_BOHR_IN_ANGSTROM 0.52917721092

// The cube's name in attribute paths: the stem of its file, of letters, digits, '_' and '-'
static md_str_t cube_name(char* buf, size_t cap, md_str_t path) {
	const md_str_t stem = md_path_stem(path);
	const size_t n = MD_MIN(stem.len, cap - 1);
	for (size_t i = 0; i < n; ++i) {
		const char c = stem.ptr[i];
		const bool keep = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
		buf[i] = keep ? c : '_';
	}
	buf[n] = '\0';
	return (md_str_t){ buf, n };
}

// The header alone, without reading the values: the lines before the atoms, then what the atoms can take of the rest
static bool cube_load_header(md_cube_t* cube, md_str_t path, md_allocator_i* alloc) {
	md_file_t file = {0};
	if (!md_file_open(&file, path, MD_FILE_READ)) {
		return false;
	}
	const size_t size = (size_t)md_file_size(file);
	char head[1024];
	const size_t head_len = md_file_read(file, head, MD_MIN(sizeof(head), size));

	// The atom count, on line 3
	md_str_t str = { head, head_len };
	md_str_t line = {0};
	int32_t count = 0;
	bool ok = md_parse_next_line(&line, &str) && md_parse_next_line(&line, &str) && md_parse_next_line(&line, &str) &&
			  md_parse_next_i32(&count, &line);
	if (ok) {
		// An atom's line is some 60 characters
		const size_t cap = MD_MIN(size, sizeof(head) + (size_t)abs(count) * 256);
		char* buf = (char*)md_alloc(alloc, cap);
		md_memcpy(buf, head, head_len);
		const size_t len = head_len + (cap > head_len ? md_file_read(file, buf + head_len, cap - head_len) : 0);
		str = (md_str_t){ buf, len };
		ok = cube_parse_header(cube, &str, alloc);
		if (ok && cube->atom.count < 0) cube->atom.count = -cube->atom.count;
	}
	md_file_close(&file);
	return ok;
}

// Which atoms of the system the cube's are, in its order: all of them, or its QM atoms (qm/atom/system_index, as a
// supplemental VeloxChem file publishes them). False when the elements do not agree with either.
static bool cube_match(uint32_t* out_index, const md_system_t* sys, const md_cube_t* cube, md_allocator_i* temp) {
	const size_t n = (size_t)cube->atom.count;
	if (n == 0) return false;

	if (n == sys->atom.count) {
		bool match = true;
		for (size_t i = 0; i < n && match; ++i) {
			match = md_atom_atomic_number(&sys->atom, i) == cube->atom.number[i];
		}
		if (match) {
			for (size_t i = 0; i < n; ++i) out_index[i] = (uint32_t)i;
			return true;
		}
	}

	const md_attribute_t* attr = md_attributes_find(&sys->attributes, MD_STR_LIT("qm/atom/system_index"));
	if (attr && md_attribute_element_count(&attr->format) == n) {
		double* idx = (double*)md_alloc(temp, n * sizeof(double));
		if (md_attribute_extract_f64(idx, n, &sys->attributes, attr, md_attribute_slice_all(), attr->unit) == n) {
			bool match = true;
			for (size_t i = 0; i < n && match; ++i) {
				const size_t j = (size_t)idx[i];
				match = j < sys->atom.count && md_atom_atomic_number(&sys->atom, j) == cube->atom.number[i];
				out_index[i] = (uint32_t)j;
			}
			if (match) return true;
		}
	}
	return false;
}

// cube/<name>               the values: {x, y, z}, or {m, x, y, z} for m values a point (an MO cube), x outermost
// cube/<name>/origin        the corner of the first voxel, Å (a cube's own origin is that voxel's centre)
// cube/<name>/basis         one step along x, y and z, Å, rows
// cube/<name>/id            the ids of the m values (MO numbers), where the cube gives them
// cube/<name>/title, comment
// cube/<name>/atom/{atomic_number, coordinate, system_index}
//                           the cube's atoms in its frame (Å), and which of the system's they are
static bool cube_publish(md_system_t* sys, const md_cube_t* cube, md_str_t name, const uint32_t* system_index, md_allocator_i* temp) {
	char root_buf[256];
	const md_str_t root = md_qm_attribute_path(root_buf, sizeof(root_buf), MD_STR_LIT("cube"), name);
	if (md_str_empty(root)) return false;
	md_attributes_remove_prefix(&sys->attributes, root);

	const uint32_t nx = (uint32_t)cube->data.num_x, ny = (uint32_t)cube->data.num_y, nz = (uint32_t)cube->data.num_z;
	const uint32_t m = (uint32_t)cube->data.num_m;
	const size_t num_vox = (size_t)nx * ny * nz;

	// [x][y][z][m] to [m][x][y][z]
	const float* src = cube->data.val;
	float* values = (float*)src;
	if (m > 1) {
		values = (float*)md_alloc(temp, num_vox * m * sizeof(float));
		for (size_t v = 0; v < num_vox; ++v) {
			for (uint32_t k = 0; k < m; ++k) values[k * num_vox + v] = src[v * m + k];
		}
	}
	md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1 };
	if (m == 1) {
		format.rank = 3;
		format.shape[0] = nx; format.shape[1] = ny; format.shape[2] = nz;
	} else {
		format.rank = 4;
		format.shape[0] = m; format.shape[1] = nx; format.shape[2] = ny; format.shape[3] = nz;
	}
	const md_str_t title = md_str_trim(cube->title);
	if (md_qm_publish(sys, root, md_str_empty(title) ? name : title, md_unit_none(), format, values, num_vox * m * sizeof(float)) == MD_ATTRIBUTE_INVALID) {
		return false;
	}

	char path_buf[320];
	// root fits in root_buf, so root and a leaf fit here
#define CUBE_PATH(leaf) md_str_from_cstrn(path_buf, (size_t)md_snprintf(path_buf, sizeof(path_buf), MD_STR_FMT "/" leaf, MD_STR_ARG(root)))
	const float* ax[3] = { cube->xaxis, cube->yaxis, cube->zaxis };
	double basis[9];
	md_double3_t origin;
	double* o = &origin.x;
	for (int i = 0; i < 3; ++i) {
		o[i] = cube->origin[i];
		for (int a = 0; a < 3; ++a) {
			o[i] -= 0.5 * ax[a][i];
			basis[3 * a + i] = ax[a][i] * CUBE_BOHR_IN_ANGSTROM;
		}
		o[i] *= CUBE_BOHR_IN_ANGSTROM;
	}
	md_qm_publish_origin(sys, CUBE_PATH("origin"), origin);
	md_qm_publish_matrix(sys, CUBE_PATH("basis"), MD_STR_LIT("Basis"), md_unit_angstrom(), basis, 3, 3);
	if (cube->data.id) {
		const md_attribute_format_t f = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { m } };
		md_qm_publish(sys, CUBE_PATH("id"), MD_STR_LIT("Id"), md_unit_none(), f, cube->data.id, m * sizeof(int));
	}
	if (!md_str_empty(title)) md_qm_publish_str(sys, CUBE_PATH("title"), MD_STR_LIT("Title"), title);
	const md_str_t comment = md_str_trim(cube->comment);
	if (!md_str_empty(comment)) md_qm_publish_str(sys, CUBE_PATH("comment"), MD_STR_LIT("Comment"), comment);

	const size_t n = (size_t)cube->atom.count;
	if (n > 0) {
		md_double3_t* xyz = (md_double3_t*)md_alloc(temp, n * sizeof(md_double3_t));
		for (size_t i = 0; i < n; ++i) {
			xyz[i] = (md_double3_t){ cube->atom.coord[i][0] * CUBE_BOHR_IN_ANGSTROM, cube->atom.coord[i][1] * CUBE_BOHR_IN_ANGSTROM, cube->atom.coord[i][2] * CUBE_BOHR_IN_ANGSTROM };
		}
		const md_attribute_format_t fi = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
		md_qm_publish(sys, CUBE_PATH("atom/atomic_number"), MD_STR_LIT("Atomic Number"), md_unit_none(), fi, cube->atom.number, n * sizeof(int));
		md_qm_publish_vec3_series(sys, CUBE_PATH("atom/coordinate"), MD_STR_LIT("Coordinate"), md_unit_angstrom(), xyz, n);
		if (system_index) {
			const md_attribute_format_t fu = { .type = MD_ATTRIBUTE_TYPE_U32, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
			md_qm_publish(sys, CUBE_PATH("atom/system_index"), MD_STR_LIT("System Atom Index"), md_unit_none(), fu, system_index, n * sizeof(uint32_t));
		}
	}
#undef CUBE_PATH
	return true;
}

bool md_cube_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t path, const md_load_desc_t* desc) {
	(void)desc;
	MD_ASSERT(sys);
	if (!sys->alloc || !state || !state->alloc) {
		MD_LOG_ERROR("CUBE: System or state allocator not set");
		return false;
	}
	md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
	md_cube_t cube = {0};
	bool ok = md_cube_file_load(&cube, path, temp);
	const size_t n = ok ? (size_t)cube.atom.count : 0;
	if (ok && n == 0) {
		MD_LOG_ERROR("CUBE: '" MD_STR_FMT "' has no atoms to make a system of", MD_STR_ARG(path));
		ok = false;
	}
	if (ok) {
		md_system_reset(sys);
		md_system_state_init(state, n);
		md_array_resize(sys->atom.type_idx, MD_ROUND_UP(n, 16), sys->alloc);
		md_memset(sys->atom.type_idx, 0, md_array_bytes(sys->atom.type_idx));
		md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
		uint32_t* index = (uint32_t*)md_alloc(temp, n * sizeof(uint32_t));
		for (size_t i = 0; i < n; ++i) {
			const float s = (float)CUBE_BOHR_IN_ANGSTROM;
			state->xyz[i] = md_float3_set(cube.atom.coord[i][0] * s, cube.atom.coord[i][1] * s, cube.atom.coord[i][2] * s);
			const md_atomic_number_t z = (md_atomic_number_t)cube.atom.number[i];
			sys->atom.type_idx[i] = md_particle_type_find_or_add(&sys->atom.type, md_atomic_number_symbol(z), z,
				md_atomic_number_mass(z), md_atomic_number_vdw_radius(z), md_atomic_number_cpk_color(z), 0, sys->alloc);
			index[i] = (uint32_t)i;
		}
		sys->atom.count  = n;
		state->num_atoms = n;
		char name_buf[128];
		ok = cube_publish(sys, &cube, cube_name(name_buf, sizeof(name_buf), path), index, temp);
	}
	md_arena_allocator_destroy(temp);
	return ok;
}

bool md_cube_system_supplement_from_file(md_system_t* sys, md_str_t path, const md_load_desc_t* desc) {
	(void)desc;
	MD_ASSERT(sys);
	md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
	md_cube_t cube = {0};
	bool ok = md_cube_file_load(&cube, path, temp);
	uint32_t* index = NULL;
	if (ok) {
		index = (uint32_t*)md_alloc(temp, MD_MAX(1, (size_t)cube.atom.count) * sizeof(uint32_t));
		if (!cube_match(index, sys, &cube, temp)) {
			MD_LOG_ERROR("CUBE: The atoms of '" MD_STR_FMT "' are not those of the system, nor its QM atoms", MD_STR_ARG(path));
			ok = false;
		}
	}
	if (ok) {
		char name_buf[128];
		ok = cube_publish(sys, &cube, cube_name(name_buf, sizeof(name_buf), path), index, temp);
	}
	md_arena_allocator_destroy(temp);
	return ok;
}

bool md_cube_system_is_file_supplemental(const md_system_t* sys, md_str_t path) {
	if (!sys || sys->atom.count == 0) return false;
	md_allocator_i* temp = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
	md_cube_t cube = {0};
	bool ok = cube_load_header(&cube, path, temp);
	if (ok) {
		uint32_t* index = (uint32_t*)md_alloc(temp, MD_MAX(1, (size_t)cube.atom.count) * sizeof(uint32_t));
		ok = cube_match(index, sys, &cube, temp);
	}
	md_arena_allocator_destroy(temp);
	return ok;
}

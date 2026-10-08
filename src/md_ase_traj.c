#include <md_ase_traj.h>

#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_json.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <math.h>

// ASE trajectories are ULM files (ase/io/ulm.py, ase/io/trajectory.py). As far as they are read here:
//
//     0   "- of Ulm"              magic
//     8   "ASE-Trajectory  "      tag, padded with spaces to 16 bytes
//     24  int64                   ULM version, 3
//     32  int64                   number of items, one per frame
//     40  int64                   offset of the item table: an int64 offset per item
//
// An item is an int64 byte count followed by that many bytes of JSON, an object. Arrays are kept out
// of the JSON: a key ending in '.' holds {"ndarray": [shape, dtype, offset]}, the data at that
// absolute offset. All of it is little endian unless an item says "_little_endian": false, which a
// big endian machine writes and which is refused here.
//
// Every frame has "positions." and "cell", whose rows are the box vectors. The first frame also has
// the header - "version", "pbc", "numbers." - which ASE writes again only in a frame where something
// in it changed; a frame without one has the first frame's, as ASE reads it.
//
// Arrays are read as the host's numbers: like the rest of mdlib this assumes a little endian host.

#define ASE_HEADER_SIZE   48
#define ASE_MAX_ITEM_SIZE MEGABYTES(64)
#define ASE_CELL_EPS      1.0e-6

typedef struct ase_index_t {
    size_t              num_frames;
    size_t              num_atoms;
    md_atomic_number_t* numbers;    // num_atoms, from the first frame
    int64_t*            offset;     // num_frames: where the positions of each frame start
    int64_t*            size;       // num_frames: their byte size, num_atoms * 3 * (4 or 8)
    float*              cell;       // num_frames * 9: row i box vector i, zero without periodicity
    double*             time;       // num_frames
    md_unit_t           time_unit;
} ase_index_t;

// {"ndarray": [shape, dtype, offset]}
typedef struct ase_array_t {
    int64_t offset;
    size_t  dim[2];
    size_t  width;      // bytes per element
    bool    is_float;
} ase_array_t;

// The array an item holds under key, of the given rank. False when it is missing, of a dtype not
// read here, or does not fit in the file.
static bool ase_array(ase_array_t* out, md_json_val_t item, str_t key, size_t rank, int64_t file_size) {
    static const struct { str_t name; size_t width; bool is_float; } dtypes[] = {
        { STR_INIT("float64"), 8, true  },
        { STR_INIT("float32"), 4, true  },
        { STR_INIT("int64"),   8, false },
        { STR_INIT("int32"),   4, false },   // numpy's default integer on Windows before numpy 2
    };

    const md_json_val_t nd    = md_json_get(md_json_get(item, key), STR_LIT("ndarray"));
    const md_json_val_t shape = md_json_at(nd, 0);
    const md_json_val_t dtype = md_json_at(nd, 1);
    if (md_json_type(nd) != MD_JSON_TYPE_ARRAY || md_json_count(nd) != 3 ||
        md_json_type(shape) != MD_JSON_TYPE_ARRAY || md_json_count(shape) != rank ||
        !md_json_i64(&out->offset, md_json_at(nd, 2))) {
        return false;
    }

    out->width = 0;
    for (size_t i = 0; i < ARRAY_SIZE(dtypes); ++i) {
        if (md_json_string_eq(dtype, dtypes[i].name)) {
            out->width    = dtypes[i].width;
            out->is_float = dtypes[i].is_float;
            break;
        }
    }
    if (out->width == 0) return false;

    // Each dimension is held to what the file has room for, so the byte count cannot overflow
    int64_t bytes = (int64_t)out->width;
    size_t i = 0;
    for (md_json_val_t d = md_json_first(shape); md_json_valid(d); d = md_json_next(d), ++i) {
        int64_t n = 0;
        if (!md_json_i64(&n, d) || n <= 0 || n > file_size / bytes) return false;
        out->dim[i] = (size_t)n;
        bytes *= n;
    }
    return 0 <= out->offset && out->offset <= file_size - bytes;
}

// The JSON object of the item at offset, parsed in temp along with its text. NONE when it is not one.
static md_json_val_t ase_item(md_file_t file, int64_t file_size, int64_t offset, md_temp_scope_t temp) {
    const md_json_val_t none = {0};
    int64_t len = 0;
    if (offset < ASE_HEADER_SIZE || offset > file_size - 8 || md_file_read_at(file, offset, &len, sizeof(len)) != sizeof(len) ||
        len <= 0 || len > ASE_MAX_ITEM_SIZE || len > file_size - offset - 8) {
        return none;
    }
    char* text = md_temp_alloc(temp, (size_t)len);
    if (!text || md_file_read_at(file, offset + 8, text, (size_t)len) != (size_t)len) {
        return none;
    }
    const md_json_val_t root = md_json_root(md_json_parse((str_t){ text, (size_t)len }, md_temp_allocator(temp), NULL));
    return md_json_type(root) == MD_JSON_TYPE_OBJECT ? root : none;
}

// "numbers." of an item, num_atoms atomic numbers, into out. Scratch from temp.
static bool ase_numbers(md_atomic_number_t* out, md_json_val_t item, size_t num_atoms, md_file_t file, int64_t file_size, md_temp_scope_t temp) {
    ase_array_t arr;
    if (!ase_array(&arr, item, STR_LIT("numbers."), 1, file_size) || arr.is_float || arr.dim[0] != num_atoms) {
        return false;
    }
    const size_t bytes = num_atoms * arr.width;
    void* raw = md_temp_alloc(temp, bytes);
    if (!raw || md_file_read_at(file, arr.offset, raw, bytes) != bytes) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        const int64_t z = (arr.width == 8) ? ((const int64_t*)raw)[i] : ((const int32_t*)raw)[i];
        if (z < 0 || z >= MD_Z_Count) return false;
        out[i] = (md_atomic_number_t)z;
    }
    return true;
}

static bool ase_pbc(bool out[3], md_json_val_t pbc) {
    if (md_json_type(pbc) != MD_JSON_TYPE_ARRAY || md_json_count(pbc) != 3) return false;
    for (size_t i = 0; i < 3; ++i) {
        if (!md_json_bool(&out[i], md_json_at(pbc, i))) return false;
    }
    return true;
}

static bool ase_mat3(double out[3][3], md_json_val_t mat) {
    if (md_json_type(mat) != MD_JSON_TYPE_ARRAY || md_json_count(mat) != 3) return false;
    for (size_t i = 0; i < 3; ++i) {
        const md_json_val_t row = md_json_at(mat, i);
        if (md_json_type(row) != MD_JSON_TYPE_ARRAY || md_json_count(row) != 3 || md_json_extract_f64(out[i], 3, row) != 3) {
            return false;
        }
    }
    return true;
}

// The cell of a frame as a run holds it. md_unitcell_t has a along x and b in the xy plane, and is
// periodic along every axis its diagonal is non zero along. A fully periodic cell in that orientation
// is kept. Without periodicity the cell is dropped: the box ASE keeps around a molecule would
// otherwise be taken to be a periodic one. Anything else - a rotated cell, or one periodic along
// only some axes - is refused rather than shown wrong.
static bool ase_cell(float out[9], const double A[3][3], const bool pbc[3]) {
    MEMSET(out, 0, 9 * sizeof(float));
    if (!pbc[0] && !pbc[1] && !pbc[2]) {
        return true;
    }
    if (!pbc[0] || !pbc[1] || !pbc[2] ||
        fabs(A[0][1]) > ASE_CELL_EPS || fabs(A[0][2]) > ASE_CELL_EPS || fabs(A[1][2]) > ASE_CELL_EPS ||
        !(A[0][0] > 0.0 && A[1][1] > 0.0 && A[2][2] > 0.0)) {
        return false;
    }
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            out[i * 3 + j] = (float)A[i][j];
        }
    }
    return true;
}

// Frame i of the index from its item: the first frame sets the atoms and the header, the others
// must keep the atoms. What the index keeps goes in out, scratch in temp. NULL when the frame is fine,
// otherwise what is wrong with it.
static const char* ase_frame_read(ase_index_t* idx, size_t i, md_json_val_t item, bool pbc0[3], bool* has_time,
                                  md_file_t file, int64_t file_size, md_temp_scope_t out, md_temp_scope_t temp) {
    if (!md_json_valid(item)) {
        return "not a JSON object";
    }
    bool little_endian = true;
    md_json_bool(&little_endian, md_json_get(item, STR_LIT("_little_endian")));
    if (!little_endian) {
        return "written on a big endian machine, which is not supported";
    }

    ase_array_t pos;
    if (!ase_array(&pos, item, STR_LIT("positions."), 2, file_size) || !pos.is_float || pos.dim[1] != 3) {
        return "no positions";
    }

    if (i == 0) {
        int64_t version = 0;
        if (!md_json_i64(&version, md_json_get(item, STR_LIT("version"))) || version != 1) {
            return "not version 1 of the trajectory format";
        }
        idx->num_atoms = pos.dim[0];
        idx->numbers   = md_temp_alloc_array(out, md_atomic_number_t, idx->num_atoms);
        if (!idx->numbers || !ase_numbers(idx->numbers, item, idx->num_atoms, file, file_size, temp)) {
            return "no atomic numbers";
        }
    } else if (pos.dim[0] != idx->num_atoms) {
        return "a different number of atoms than the first frame, which a run cannot hold";
    } else if (md_json_valid(md_json_get(item, STR_LIT("numbers.")))) {
        md_atomic_number_t* numbers = md_temp_alloc_array(temp, md_atomic_number_t, idx->num_atoms);
        if (!numbers || !ase_numbers(numbers, item, idx->num_atoms, file, file_size, temp)) {
            return "malformed atomic numbers";
        }
        if (MEMCMP(numbers, idx->numbers, idx->num_atoms * sizeof(md_atomic_number_t)) != 0) {
            return "different atoms than the first frame, which a run cannot hold";
        }
    }

    const md_json_val_t pbc_val = md_json_get(item, STR_LIT("pbc"));
    bool pbc[3] = { pbc0[0], pbc0[1], pbc0[2] };
    if ((i == 0 || md_json_valid(pbc_val)) && !ase_pbc(pbc, pbc_val)) {
        return "no pbc";
    }
    if (i == 0) {
        MEMCPY(pbc0, pbc, sizeof(pbc));
    }

    double A[3][3];
    if (!ase_mat3(A, md_json_get(item, STR_LIT("cell")))) {
        return "no cell";
    }
    if (!ase_cell(idx->cell + i * 9, A, pbc)) {
        return "a cell that is rotated or periodic along only some axes, which mdlib cannot represent";
    }

    // Not an ASE convention: some workflows keep the simulation time in info
    *has_time = *has_time && md_json_f64(&idx->time[i], md_json_get(md_json_get(item, STR_LIT("info")), STR_LIT("time_ps")));

    idx->offset[i] = pos.offset;
    idx->size[i]   = (int64_t)(idx->num_atoms * 3 * pos.width);
    return NULL;
}

// The frames of an open file, as many as max_frames. The index goes in out, scratch in temp, which
// must be in another arena.
static bool ase_index_read_file(ase_index_t* idx, md_file_t file, str_t path, size_t max_frames, md_temp_scope_t out, md_temp_scope_t temp) {
    const int64_t file_size = (int64_t)md_file_size(file);
    uint8_t header[ASE_HEADER_SIZE];
    if (file_size < ASE_HEADER_SIZE || md_file_read_at(file, 0, header, sizeof(header)) != sizeof(header) ||
        MEMCMP(header, "- of Ulm", 8) != 0 || MEMCMP(header + 8, "ASE-Trajectory  ", 16) != 0) {
        MD_LOG_ERROR("ASE: '" STR_FMT "' is not an ASE trajectory", STR_ARG(path));
        return false;
    }

    int64_t version, num_items, table_offset;
    MEMCPY(&version,      header + 24, sizeof(int64_t));
    MEMCPY(&num_items,    header + 32, sizeof(int64_t));
    MEMCPY(&table_offset, header + 40, sizeof(int64_t));
    if (version != 3) {
        MD_LOG_ERROR("ASE: '" STR_FMT "' is of ULM version %lld, only version 3 is supported", STR_ARG(path), (long long)version);
        return false;
    }
    if (num_items <= 0 || table_offset < ASE_HEADER_SIZE || table_offset > file_size || num_items > (file_size - table_offset) / 8) {
        MD_LOG_ERROR("ASE: '" STR_FMT "' has no frames, or a damaged table of them", STR_ARG(path));
        return false;
    }

    const size_t F = MIN((size_t)num_items, max_frames);
    int64_t* items  = md_temp_alloc_array(temp, int64_t, F);
    idx->num_frames = F;
    idx->offset     = md_temp_alloc_array(out, int64_t, F);
    idx->size       = md_temp_alloc_array(out, int64_t, F);
    idx->cell       = md_temp_alloc_array(out, float, F * 9);
    idx->time       = md_temp_alloc_array(out, double, F);
    if (!items || !idx->offset || !idx->size || !idx->cell || !idx->time ||
        md_file_read_at(file, table_offset, items, F * sizeof(int64_t)) != F * sizeof(int64_t)) {
        MD_LOG_ERROR("ASE: could not read the table of frames of '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    bool pbc0[3] = {0};
    bool has_time = true;
    for (size_t i = 0; i < F; ++i) {
        // Each frame's text and document are let go of before the next one is read
        md_temp_scope_t frame_temp = md_temp_begin_in(temp.arena);
        const md_json_val_t item = ase_item(file, file_size, items[i], frame_temp);
        const char* error = ase_frame_read(idx, i, item, pbc0, &has_time, file, file_size, out, frame_temp);
        md_temp_end(frame_temp);
        if (error) {
            MD_LOG_ERROR("ASE: frame %zu of '" STR_FMT "': %s", i, STR_ARG(path), error);
            return false;
        }
    }

    idx->time_unit = md_unit_picosecond();
    if (!has_time) {
        // Nothing says when the frames are: ordinals
        for (size_t i = 0; i < F; ++i) {
            idx->time[i] = (double)i;
        }
        idx->time_unit = md_unit_none();
    }
    return true;
}

// The first max_frames frames of the file at path. The index goes in out, the caller's scope, and is
// gone when it ends; what is only needed while reading comes from the other temp arena.
static bool ase_index_read(ase_index_t* idx, str_t path, size_t max_frames, md_temp_scope_t out) {
    MEMSET(idx, 0, sizeof(ase_index_t));
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("ASE: could not open '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(out.arena);
    const bool result = ase_index_read_file(idx, file, path, max_frames, out, temp);
    md_temp_end(temp);
    md_file_close(&file);
    return result;
}

// Bytes per coordinate of a frame of num_atoms positions with byte size: 4 or 8, 0 for anything else
static size_t ase_frame_width(int64_t size, size_t num_atoms) {
    if (size == (int64_t)(num_atoms * 3 * sizeof(float)))  return sizeof(float);
    if (size == (int64_t)(num_atoms * 3 * sizeof(double))) return sizeof(double);
    return 0;
}

// count positions stored with width bytes per coordinate, at offset in the file at path, into dst as
// count * 3 floats: float32 straight into place, float64 through a buffer. Through io when given.
static bool ase_positions_read(float* dst, size_t count, size_t width, md_attribute_io_t* io, str_t path, int64_t offset) {
    const size_t n = count * 3;
    if (width == sizeof(float)) {
        return md_attribute_io_read_at(io, path, offset, dst, n * sizeof(float)) == n * sizeof(float);
    }
    if (width != sizeof(double)) {
        return false;
    }
    md_temp_scope_t temp = md_temp_begin();
    double* raw = md_temp_alloc_array(temp, double, n);
    const bool ok = raw && md_attribute_io_read_at(io, path, offset, raw, n * sizeof(double)) == n * sizeof(double);
    for (size_t i = 0; ok && i < n; ++i) {
        dst[i] = (float)raw[i];
    }
    md_temp_end(temp);
    return ok;
}

// <run>/atom/position: one read of the frame, or of the one atom asked for
static size_t ase_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

    const size_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (frame >= src.num_frames || cap != count * 3) return 0;

    const size_t width = ase_frame_width(src.size[frame], N);
    if (!ase_positions_read((float*)dst, count, width, io, src.path, src.offset[frame] + (int64_t)(first * 3 * width))) {
        MD_LOG_ERROR("ASE: failed to read frame %zu of '" STR_FMT "'", frame, STR_ARG(src.path));
        return 0;
    }
    return cap;
}

bool md_ase_traj_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);
    ASSERT(state);
    if (!sys->alloc || !state->alloc) {
        MD_LOG_ERROR("ASE: system or state allocator not set");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);

    // The first frame only, a run reads the others. All of it is read before the system is touched.
    ase_index_t idx;
    float* xyz = NULL;
    bool result = ase_index_read(&idx, filename, 1, temp);
    if (result) {
        xyz = md_temp_alloc_array(temp, float, idx.num_atoms * 3);
        result = xyz && ase_positions_read(xyz, idx.num_atoms, ase_frame_width(idx.size[0], idx.num_atoms), NULL, filename, idx.offset[0]);
    }
    if (result) {
        md_system_reset(sys);
        result = md_system_state_init(state, idx.num_atoms);
    }
    if (result) {
        const size_t N = idx.num_atoms;
        MEMCPY(state->xyz, xyz, N * sizeof(vec3_t));
        state->unitcell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(idx.cell));

        md_array_ensure(sys->atom.type_idx, N, sys->alloc);
        md_array_ensure(sys->atom.flags,    N, sys->alloc);
        md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
        for (size_t i = 0; i < N; ++i) {
            const md_atomic_number_t z = idx.numbers[i];
            const md_atom_type_idx_t type = md_atom_type_find_or_add(&sys->atom.type, md_atomic_number_symbol(z), z,
                md_atomic_number_mass(z), md_atomic_number_vdw_radius(z), md_atomic_number_cpk_color(z), 0, sys->alloc);
            md_array_push(sys->atom.type_idx, type, sys->alloc);
            md_array_push(sys->atom.flags, 0, sys->alloc);
        }
        sys->atom.count = N;
    }

    md_temp_end(temp);
    return result;
}

bool md_ase_traj_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    (void)flags;    // no index cache to write: the file has a table of its frames
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    if (path_len == 0) {
        MD_LOG_ERROR("ASE: could not resolve the path '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    const str_t path = { path_buf, path_len };

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);

    ase_index_t idx;
    bool result = ase_index_read(&idx, path, SIZE_MAX, temp);
    if (result && idx.num_frames < 2) {
        // A structure, as a PDB of one model is: md_ase_traj_system_init_from_file has all of it
        MD_LOG_INFO("ASE: '" STR_FMT "' has a single frame and is not read as a trajectory", STR_ARG(path));
        result = false;
    }
    if (result) {
        const md_attribute_virtual_t virt = { .provider = ase_position_provider, .user_data = sys };
        const md_run_desc_t desc = {
            .num_frames    = idx.num_frames,
            .num_atoms     = idx.num_atoms,
            .time          = idx.time,
            .time_unit     = idx.time_unit,
            .unitcell      = idx.cell,
            .source_path   = path,
            .source_offset = idx.offset,
            .source_size   = idx.size,
            .position_virt = &virt,
        };
        result = md_run_publish(sys, run, &desc);
    }

    md_temp_end(temp);
    return result;
}

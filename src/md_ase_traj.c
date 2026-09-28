#include <md_ase_traj.h>

#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>

#define JSMN_STATIC
#define JSMN_STRICT
#define JSMN_PARENT_LINKS
#include "../ext/jsmn/jsmn.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_JSON_BYTES = 16 * 1024 * 1024, MAX_JSON_TOKENS = 65536 };

typedef struct {
    uint64_t offset;
    uint8_t width;
    bool little_endian;
    double cell[3][3];
    uint8_t pbc;
} ase_frame_t;

typedef struct {
    md_allocator_i* arena;
    size_t num_frames;
    size_t num_atoms;
    int32_t* numbers;
    ase_frame_t* frames;
    double* times;
    md_unit_t time_unit;
} ase_traj_t;

typedef struct {
    char* text;
    jsmntok_t* tok;
    int count;
} json_doc_t;

typedef struct {
    uint64_t offset;
    size_t count;
    uint8_t width;
    bool is_float;
    bool is_signed;
} array_desc_t;

static uint64_t read_integer(const uint8_t* p, uint8_t width, bool little_endian) {
    uint64_t value = 0;
    for (uint8_t i = 0; i < width; ++i) {
        value = (value << 8) | p[little_endian ? width - i - 1 : i];
    }
    return value;
}

static bool read_at(md_file_t file, uint64_t file_size, uint64_t offset, void* dst, size_t size) {
    if (offset > file_size || size > file_size - offset || offset > INT64_MAX) return false;
    return md_file_seek(file, (md_file_offset_t)offset, MD_FILE_BEG) &&
           md_file_read(file, dst, size) == size;
}

static bool token_eq(const json_doc_t* doc, int idx, const char* text) {
    if (idx < 0 || idx >= doc->count) return false;
    const jsmntok_t* t = &doc->tok[idx];
    size_t size = strlen(text);
    return t->type == JSMN_STRING && t->end - t->start == (int)size &&
           memcmp(doc->text + t->start, text, size) == 0;
}

static int object_value(const json_doc_t* doc, int object, const char* key) {
    if (object < 0 || object >= doc->count || doc->tok[object].type != JSMN_OBJECT) return -1;
    for (int i = object + 1; i + 1 < doc->count && doc->tok[i].start < doc->tok[object].end; ++i) {
        if (doc->tok[i].parent == object && token_eq(doc, i, key)) return i + 1;
    }
    return -1;
}

static int array_value(const json_doc_t* doc, int array, int index) {
    if (array < 0 || array >= doc->count || doc->tok[array].type != JSMN_ARRAY || index < 0) return -1;
    int found = 0;
    for (int i = array + 1; i < doc->count && doc->tok[i].start < doc->tok[array].end; ++i) {
        if (doc->tok[i].parent == array && found++ == index) return i;
    }
    return -1;
}

static bool token_uint(const json_doc_t* doc, int idx, uint64_t* out) {
    if (idx < 0 || idx >= doc->count || doc->tok[idx].type != JSMN_PRIMITIVE) return false;
    const jsmntok_t* t = &doc->tok[idx];
    if (t->start == t->end) return false;
    uint64_t value = 0;
    for (int i = t->start; i < t->end; ++i) {
        char c = doc->text[i];
        if (c < '0' || c > '9' || value > (UINT64_MAX - (uint64_t)(c - '0')) / 10) return false;
        value = value * 10 + (uint64_t)(c - '0');
    }
    *out = value;
    return true;
}

static bool token_double(const json_doc_t* doc, int idx, double* out) {
    if (idx < 0 || idx >= doc->count || doc->tok[idx].type != JSMN_PRIMITIVE) return false;
    const jsmntok_t* t = &doc->tok[idx];
    int size = t->end - t->start;
    if (size < 1 || size > 63) return false;
    char text[64];
    memcpy(text, doc->text + t->start, (size_t)size);
    text[size] = '\0';
    char* end = NULL;
    double value = strtod(text, &end);
    if (end != text + size || !isfinite(value)) return false;
    *out = value;
    return true;
}

static bool token_bool(const json_doc_t* doc, int idx, bool* out) {
    if (idx < 0 || idx >= doc->count || doc->tok[idx].type != JSMN_PRIMITIVE) return false;
    const jsmntok_t* t = &doc->tok[idx];
    if (t->end - t->start == 4 && memcmp(doc->text + t->start, "true", 4) == 0) {
        *out = true;
        return true;
    }
    if (t->end - t->start == 5 && memcmp(doc->text + t->start, "false", 5) == 0) {
        *out = false;
        return true;
    }
    return false;
}

static void json_free(json_doc_t* doc, size_t text_size, unsigned token_cap) {
    md_allocator_i* heap = md_get_heap_allocator();
    if (doc->text) md_free(heap, doc->text, text_size);
    if (doc->tok) md_free(heap, doc->tok, (size_t)token_cap * sizeof(jsmntok_t));
    *doc = (json_doc_t){0};
}

static bool json_read(md_file_t file, uint64_t file_size, uint64_t offset,
                      json_doc_t* doc, size_t* text_size, unsigned* token_cap) {
    uint8_t size_bytes[8];
    if (!read_at(file, file_size, offset, size_bytes, sizeof(size_bytes))) return false;
    uint64_t json_size = read_integer(size_bytes, 8, true);
    if (json_size == 0 || json_size > MAX_JSON_BYTES || offset + 8 > file_size ||
        json_size > file_size - offset - 8) return false;
    md_allocator_i* heap = md_get_heap_allocator();
    *text_size = (size_t)json_size + 1;
    doc->text = md_alloc(heap, *text_size);
    if (!doc->text || !read_at(file, file_size, offset + 8, doc->text, (size_t)json_size)) return false;
    doc->text[json_size] = '\0';
    for (*token_cap = 256; *token_cap <= MAX_JSON_TOKENS; *token_cap *= 2) {
        doc->tok = md_alloc(heap, (size_t)*token_cap * sizeof(jsmntok_t));
        if (!doc->tok) return false;
        jsmn_parser parser;
        jsmn_init(&parser);
        doc->count = jsmn_parse(&parser, doc->text, (size_t)json_size, doc->tok, *token_cap);
        if (doc->count >= 1) return doc->tok[0].type == JSMN_OBJECT;
        if (doc->count != JSMN_ERROR_NOMEM) return false;
        md_free(heap, doc->tok, (size_t)*token_cap * sizeof(jsmntok_t));
        doc->tok = NULL;
    }
    return false;
}

static bool parse_array(const json_doc_t* doc, int object, const char* key,
                        uint8_t ndim, uint64_t file_size, array_desc_t* out,
                        size_t* dim0, size_t* dim1) {
    int desc = object_value(doc, object, key);
    int array = object_value(doc, desc, "ndarray");
    int shape = array_value(doc, array, 0);
    int dtype = array_value(doc, array, 1);
    int offset = array_value(doc, array, 2);
    uint64_t n0, n1 = 1, at;
    if (array < 0 || doc->tok[array].size != 3 || shape < 0 ||
        doc->tok[shape].type != JSMN_ARRAY || doc->tok[shape].size != ndim ||
        !token_uint(doc, array_value(doc, shape, 0), &n0) ||
        (ndim == 2 && !token_uint(doc, array_value(doc, shape, 1), &n1)) ||
        !token_uint(doc, offset, &at) || n0 == 0 || n0 > SIZE_MAX || n1 > SIZE_MAX ||
        n0 > SIZE_MAX / n1) return false;
    *dim0 = (size_t)n0;
    *dim1 = (size_t)n1;
    out->count = *dim0 * *dim1;
    out->offset = at;
    out->is_float = false;
    out->is_signed = false;
    if (token_eq(doc, dtype, "float64")) { out->width = 8; out->is_float = true; }
    else if (token_eq(doc, dtype, "float32")) { out->width = 4; out->is_float = true; }
    else if (token_eq(doc, dtype, "int64")) { out->width = 8; out->is_signed = true; }
    else if (token_eq(doc, dtype, "int32")) { out->width = 4; out->is_signed = true; }
    else if (token_eq(doc, dtype, "uint8")) { out->width = 1; }
    else return false;
    return out->count <= SIZE_MAX / out->width && at <= (uint64_t)INT64_MAX &&
           out->count * out->width <= (uint64_t)INT64_MAX - at &&
           at <= file_size && out->count * out->width <= file_size - at;
}

static bool parse_cell(const json_doc_t* doc, int root, double cell[3][3]) {
    int matrix = object_value(doc, root, "cell");
    if (matrix < 0 || doc->tok[matrix].type != JSMN_ARRAY || doc->tok[matrix].size != 3) return false;
    for (int i = 0; i < 3; ++i) {
        int row = array_value(doc, matrix, i);
        if (row < 0 || doc->tok[row].type != JSMN_ARRAY || doc->tok[row].size != 3) return false;
        for (int j = 0; j < 3; ++j) {
            if (!token_double(doc, array_value(doc, row, j), &cell[i][j])) return false;
        }
    }
    return true;
}

static bool parse_pbc(const json_doc_t* doc, int root, uint8_t* mask) {
    int array = object_value(doc, root, "pbc");
    if (array < 0) return true;
    if (doc->tok[array].type != JSMN_ARRAY || doc->tok[array].size != 3) return false;
    *mask = 0;
    for (int i = 0; i < 3; ++i) {
        bool periodic;
        if (!token_bool(doc, array_value(doc, array, i), &periodic)) return false;
        if (periodic) *mask |= (uint8_t)(1u << i);
    }
    return true;
}

static bool read_numbers(md_file_t file, uint64_t file_size, const array_desc_t* desc,
                         bool little_endian, int32_t* numbers, const int32_t* expected) {
    size_t bytes = desc->count * desc->width;
    uint8_t* raw = md_alloc(md_get_heap_allocator(), bytes);
    if (!raw) return false;
    bool ok = read_at(file, file_size, desc->offset, raw, bytes);
    for (size_t i = 0; ok && i < desc->count; ++i) {
        uint64_t value = read_integer(raw + i * desc->width, desc->width, little_endian);
        if (value > 118 || (expected && expected[i] != (int32_t)value)) ok = false;
        if (numbers) numbers[i] = (int32_t)value;
    }
    md_free(md_get_heap_allocator(), raw, bytes);
    return ok;
}

static bool parse_frame(md_file_t file, uint64_t file_size, uint64_t json_offset,
                        ase_traj_t* traj, size_t index, bool* has_time) {
    json_doc_t doc = {0};
    size_t text_size = 0;
    unsigned token_cap = 0;
    bool ok = false;
    if (!json_read(file, file_size, json_offset, &doc, &text_size, &token_cap)) goto done;

    ase_frame_t* frame = &traj->frames[index];
    if (index) *frame = traj->frames[index - 1];
    else {
        uint64_t version;
        if (!token_uint(&doc, object_value(&doc, 0, "version"), &version) || version != 1 ||
            object_value(&doc, 0, "pbc") < 0) goto done;
    }
    frame->little_endian = true;
    int endian = object_value(&doc, 0, "_little_endian");
    if (endian >= 0 && !token_bool(&doc, endian, &frame->little_endian)) goto done;
    if (!frame->little_endian) {
        MD_LOG_ERROR("ASE trajectory: big-endian arrays are not supported");
        goto done;
    }

    array_desc_t pos = {0};
    size_t num_atoms, width;
    if (!parse_array(&doc, 0, "positions.", 2, file_size, &pos, &num_atoms, &width) ||
        width != 3 || !pos.is_float || num_atoms > SIZE_MAX / 24 ||
        (index && num_atoms != traj->num_atoms) ||
        !parse_cell(&doc, 0, frame->cell) || !parse_pbc(&doc, 0, &frame->pbc)) goto done;
    frame->offset = pos.offset;
    frame->width = pos.width;
    // mdlib stores a triangular cell and derives periodicity from its diagonal.
    if (fabs(frame->cell[0][1]) > 1e-6 || fabs(frame->cell[0][2]) > 1e-6 ||
        fabs(frame->cell[1][2]) > 1e-6 ||
        (frame->pbc != 7 && (frame->pbc != 0 ||
            fabs(frame->cell[0][0]) > 1e-6 || fabs(frame->cell[1][0]) > 1e-6 ||
            fabs(frame->cell[1][1]) > 1e-6 || fabs(frame->cell[2][0]) > 1e-6 ||
            fabs(frame->cell[2][1]) > 1e-6 || fabs(frame->cell[2][2]) > 1e-6)) ||
        (frame->pbc == 7 && (frame->cell[0][0] <= 0 || frame->cell[1][1] <= 0 || frame->cell[2][2] <= 0))) {
        MD_LOG_ERROR("ASE trajectory: rotated cells or partial periodicity are not supported");
        goto done;
    }
    if (!index) traj->num_atoms = num_atoms;

    int number_field = object_value(&doc, 0, "numbers.");
    if (index == 0 && number_field < 0) goto done;
    if (number_field >= 0) {
        array_desc_t nums = {0};
        size_t count, unused;
        if (!parse_array(&doc, 0, "numbers.", 1, file_size, &nums, &count, &unused) ||
            nums.is_float || count != traj->num_atoms) goto done;
        if (index == 0) {
            traj->numbers = md_alloc(traj->arena, count * sizeof(int32_t));
            if (!traj->numbers) goto done;
        }
        if (!read_numbers(file, file_size, &nums, frame->little_endian,
                          index == 0 ? traj->numbers : NULL,
                          index == 0 ? NULL : traj->numbers)) goto done;
    }

    int info = object_value(&doc, 0, "info");
    int time = object_value(&doc, info, "time_ps");
    *has_time = time >= 0 && token_double(&doc, time, &traj->times[index]);
    if (time >= 0 && !*has_time) goto done;
    ok = true;
done:
    json_free(&doc, text_size, token_cap);
    return ok;
}

static double read_coordinate(const uint8_t* raw, uint8_t width, bool little_endian) {
    uint64_t bits = read_integer(raw, width, little_endian);
    if (width == 8) {
        double value;
        memcpy(&value, &bits, sizeof(value));
        return value;
    }
    uint32_t bits32 = (uint32_t)bits;
    float value;
    memcpy(&value, &bits32, sizeof(value));
    return value;
}

static ase_traj_t* ase_index_load(str_t filename, bool first_only) {
    md_allocator_i* backing = md_get_heap_allocator();
    md_allocator_i* arena = md_arena_allocator_create(backing, MEGABYTES(1));
    if (!arena) return NULL;
    ase_traj_t* traj = md_alloc(arena, sizeof(*traj));
    if (!traj) { md_arena_allocator_destroy(arena); return NULL; }
    memset(traj, 0, sizeof(*traj));
    traj->arena = arena;
    md_file_t file = {0};
    bool opened = md_file_open(&file, filename, MD_FILE_READ);
    if (!opened) goto fail;
    int64_t file_size_signed = md_file_size(file);
    if (file_size_signed < 56) goto fail;
    uint64_t file_size = (uint64_t)file_size_signed;
    uint8_t header[48];
    if (!read_at(file, file_size, 0, header, sizeof(header)) ||
        memcmp(header, "- of Ulm", 8) != 0 ||
        memcmp(header + 8, "ASE-Trajectory  ", 16) != 0 ||
        read_integer(header + 24, 8, true) != 3) goto fail;
    uint64_t count = read_integer(header + 32, 8, true);
    uint64_t table = read_integer(header + 40, 8, true);
    if (count < 1 || count > SIZE_MAX / sizeof(ase_frame_t) ||
        table > file_size || count > (file_size - table) / 8) goto fail;
    traj->num_frames = first_only ? 1 : (size_t)count;
    traj->frames = md_alloc(arena, traj->num_frames * sizeof(ase_frame_t));
    traj->times = md_alloc(arena, traj->num_frames * sizeof(double));
    if (!traj->frames || !traj->times) goto fail;
    memset(traj->frames, 0, traj->num_frames * sizeof(ase_frame_t));
    bool all_times = true;
    for (size_t i = 0; i < traj->num_frames; ++i) {
        uint8_t offset_bytes[8];
        if (!read_at(file, file_size, table + i * 8, offset_bytes, 8)) goto fail;
        uint64_t offset = read_integer(offset_bytes, 8, true);
        bool has_time = false;
        if (!parse_frame(file, file_size, offset, traj, i, &has_time)) {
            MD_LOG_ERROR("ASE trajectory: invalid or unsupported frame %zu", i);
            goto fail;
        }
        all_times &= has_time;
    }
    if (!all_times) {
        for (size_t i = 0; i < traj->num_frames; ++i) traj->times[i] = (double)i;
        traj->time_unit = md_unit_none();
    } else {
        traj->time_unit = md_unit_picosecond();
    }
    md_file_close(&file);
    return traj;
fail:
    if (opened) md_file_close(&file);
    MD_LOG_ERROR("ASE trajectory: cannot load '" STR_FMT "' as modern fixed-atom ULM", STR_ARG(filename));
    md_arena_allocator_destroy(arena);
    return NULL;
}

static bool decode_positions(float* dst, const uint8_t* raw, size_t count, uint8_t width) {
    for (size_t i = 0; i < count * 3; ++i) {
        float value = (float)read_coordinate(raw + i * width, width, true);
        if (!isfinite(value)) return false;
        dst[i] = value;
    }
    return true;
}

static size_t ase_position_provider(void* dst, size_t cap, const md_attribute_t* attr,
                                    const md_attribute_slice_t* slice, void* user_data,
                                    md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    if (!slice || slice->num_idx < 1 || slice->num_idx > 2) return 0;
    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;
    const size_t frame = slice->idx[0];
    const size_t n = attr->format.shape[1];
    if (frame >= src.num_frames || n == 0 || n > SIZE_MAX / 3) return 0;
    const size_t first = slice->num_idx == 2 ? slice->idx[1] : 0;
    const size_t count = slice->num_idx == 2 ? 1 : n;
    if (first >= n || cap != count * 3 || src.size[frame] <= 0 ||
        src.size[frame] % (int64_t)(n * 3) != 0) return 0;
    const size_t width = (size_t)(src.size[frame] / (int64_t)(n * 3));
    if (width != 4 && width != 8) return 0;
    const size_t bytes = count * 3 * width;
    md_temp_scope_t temp = md_temp_begin();
    uint8_t* raw = md_temp_alloc(temp, bytes);
    size_t written = 0;
    if (raw && md_attribute_io_read_at(io, src.path,
                                        src.offset[frame] + (int64_t)(first * 3 * width),
                                        raw, bytes) == bytes &&
        decode_positions((float*)dst, raw, count, (uint8_t)width)) {
        written = cap;
    } else {
        MD_LOG_ERROR("ASE trajectory: failed to read frame %zu", frame);
    }
    md_temp_end(temp);
    return written;
}

bool md_ase_traj_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    if (!sys || !sys->alloc || !state || !state->alloc) return false;
    // ponytail: inspect one frame here; publish_run validates the rest without scanning a large file twice.
    ase_traj_t* traj = ase_index_load(filename, true);
    if (!traj) return false;
    const size_t n = traj->num_atoms;
    const ase_frame_t* frame = &traj->frames[0];
    const size_t bytes = n * 3 * frame->width;
    uint8_t* raw = md_alloc(md_get_heap_allocator(), bytes);
    md_file_t file = {0};
    bool ok = raw && md_file_open(&file, filename, MD_FILE_READ);
    if (ok) ok = read_at(file, (uint64_t)md_file_size(file), frame->offset, raw, bytes);
    if (md_file_valid(file)) md_file_close(&file);
    if (ok) {
        md_system_reset(sys);
        ok = md_system_state_init(state, n);
        if (ok) ok = decode_positions((float*)state->xyz, raw, n, frame->width);
    }
    if (raw) md_free(md_get_heap_allocator(), raw, bytes);
    if (ok) {
        md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0, 0, 0, 0, sys->alloc);
        for (size_t i = 0; i < n; ++i) {
            const md_atomic_number_t z = (md_atomic_number_t)traj->numbers[i];
            const md_atom_type_idx_t type = md_atom_type_find_or_add(
                &sys->atom.type, md_atomic_number_symbol(z), z,
                md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
                md_atomic_number_cpk_color(z), 0, sys->alloc);
            md_array_push(sys->atom.type_idx, type, sys->alloc);
            md_array_push(sys->atom.flags, 0, sys->alloc);
        }
        sys->atom.count = n;
        state->unitcell = md_unitcell_from_matrix_double(frame->cell);
    }
    md_arena_allocator_destroy(traj->arena);
    return ok;
}

bool md_ase_traj_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    (void)flags;
    if (!sys || !sys->alloc) return false;
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    if (!path_len) return false;
    const str_t path = {path_buf, path_len};
    ase_traj_t* traj = ase_index_load(path, false);
    if (!traj) return false;
    const size_t f = traj->num_frames;
    int64_t* offsets = md_alloc(traj->arena, f * sizeof(int64_t));
    int64_t* sizes = md_alloc(traj->arena, f * sizeof(int64_t));
    float* cells = md_alloc(traj->arena, f * 9 * sizeof(float));
    bool ok = offsets && sizes && cells && (!sys->atom.count || sys->atom.count == traj->num_atoms);
    if (ok) {
        for (size_t i = 0; i < f; ++i) {
            const ase_frame_t* frame = &traj->frames[i];
            offsets[i] = (int64_t)frame->offset;
            sizes[i] = (int64_t)(traj->num_atoms * 3 * frame->width);
            for (size_t row = 0; row < 3; ++row) {
                for (size_t col = 0; col < 3; ++col) {
                    cells[i * 9 + row * 3 + col] = (float)frame->cell[row][col];
                }
            }
        }
        const md_attribute_virtual_t virt = {.provider = ase_position_provider, .user_data = sys};
        const md_run_desc_t desc = {
            .num_frames = f, .num_atoms = traj->num_atoms,
            .time = traj->times, .time_unit = traj->time_unit,
            .unitcell = cells, .source_path = path,
            .source_offset = offsets, .source_size = sizes,
            .position_virt = &virt,
        };
        ok = md_run_publish(sys, run, &desc);
    }
    md_arena_allocator_destroy(traj->arena);
    return ok;
}

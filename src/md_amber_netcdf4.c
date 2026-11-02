#include <md_amber_netcdf4.h>

#include <md_hdf5.h>
#include <md_system.h>
#include <md_unitcell.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <hdf5.h>
#include <zlib.h>

#include <math.h>
#include <string.h>

// ## NETCDF-4
//
// A NetCDF-4 file is an HDF5 file: each variable a dataset of the root group, each attribute an HDF5 attribute, each
// dimension a dimension scale (which the variables' own extents make unnecessary here). Writers chunk the per atom
// variables one frame to a chunk (netCDF4-python's default, and cpptraj's), deflated and shuffled when compressed;
// cpptraj's lossy "icompress" stores them as integers instead, as compressedpos, compressedvel and compressedfrc, the
// values multiplied by their icompressfac.
//
// HDF5 is not thread safe, and every call to it is made under the one mdlib HDF5 lock (md_hdf5.h). A frame is
// therefore not read through HDF5: where each frame's chunk lies in the file and how it is encoded is learned once,
// when the run is published, and a frame is then one read of its chunk through the io a provider is handed, inflated
// and unshuffled here, on whichever thread asks, without the lock. A layout this does not decode itself (a filter other
// than deflate and shuffle, chunks that split a frame) is read through HDF5 instead, under the lock.

enum {
    NC4_F32 = 1,
    NC4_F64,
    NC4_I32,
    NC4_I16,
};

enum {
    NC4_MODE_CONTIGUOUS = 1,    // frame f at base + f * frame_bytes
    NC4_MODE_CHUNKED,           // chunk c at its own offset, frames_per_chunk frames to a chunk
    NC4_MODE_HDF5,              // read through HDF5
};

#define NC4_PATH_MAX 4096
#define NC4_NAME_MAX 32

// How a per atom variable is laid out and encoded: what a provider needs, beside its frame table
typedef struct nc4_layout_t {
    uint32_t mode;
    uint32_t type;
    uint32_t big_endian;
    uint32_t shuffle;
    uint32_t deflate;
    uint32_t frames_per_chunk;
    uint32_t num_atoms;
    uint32_t num_chunks;
    uint64_t frame_bytes;       // stored, of one frame
    uint64_t chunk_bytes;       // stored and inflated, of one chunk
    uint64_t base;              // contiguous: where frame 0 is
    double   scale;             // what a stored value is multiplied by: scale_factor / icompressfac
} nc4_layout_t;

// A provider's own state, freed with its attribute: the layout, the dataset, the file, then the chunk table
typedef struct nc4_source_t {
    nc4_layout_t layout;
    char dataset[NC4_NAME_MAX];
    char path[NC4_PATH_MAX];
    int64_t table[];            // num_chunks offsets, then num_chunks sizes
} nc4_source_t;

static size_t nc4_type_size(uint32_t type) {
    switch (type) {
    case NC4_F32: case NC4_I32: return 4;
    case NC4_F64: return 8;
    case NC4_I16: return 2;
    default: return 0;
    }
}

// ## HDF5, with the lock held

static size_t h5_attr_text(hid_t obj, const char* name, char* buf, size_t cap) {
    buf[0] = '\0';
    if (H5Aexists(obj, name) <= 0) return 0;
    const hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
    if (a < 0) return 0;
    const hid_t t = H5Aget_type(a);
    size_t len = 0;
    if (H5Tget_class(t) == H5T_STRING) {
        if (H5Tis_variable_str(t) > 0) {
            const hid_t mt = H5Tcopy(H5T_C_S1);
            H5Tset_size(mt, H5T_VARIABLE);
            char* s = NULL;
            if (H5Aread(a, mt, &s) >= 0 && s) {
                len = MIN(strlen(s), cap - 1);
                MEMCPY(buf, s, len);
                H5free_memory(s);
            }
            H5Tclose(mt);
        } else {
            const hid_t space = H5Aget_space(a);
            const hssize_t n = H5Sget_simple_extent_npoints(space);
            H5Sclose(space);
            const size_t size = H5Tget_size(t) * (size_t)MAX(n, 1);
            if (size < cap && H5Aread(a, t, buf) >= 0) len = size;
        }
    }
    H5Tclose(t);
    H5Aclose(a);
    while (len && (buf[len - 1] == '\0' || buf[len - 1] == ' ')) len -= 1;
    buf[len] = '\0';
    return len;
}

static bool h5_attr_number(hid_t obj, const char* name, double* out) {
    if (H5Aexists(obj, name) <= 0) return false;
    const hid_t a = H5Aopen(obj, name, H5P_DEFAULT);
    if (a < 0) return false;
    const hid_t t = H5Aget_type(a);
    const hid_t space = H5Aget_space(a);
    const hssize_t n = H5Sget_simple_extent_npoints(space);
    bool ok = false;
    const H5T_class_t cls = H5Tget_class(t);
    if ((cls == H5T_FLOAT || cls == H5T_INTEGER) && n >= 1 && n <= 16) {
        double v[16];
        ok = H5Aread(a, H5T_NATIVE_DOUBLE, v) >= 0 && isfinite(v[0]);
        if (ok) *out = v[0];
    }
    H5Sclose(space);
    H5Tclose(t);
    H5Aclose(a);
    return ok;
}

static bool h5_has(hid_t file, const char* name) {
    return H5Lexists(file, name, H5P_DEFAULT) > 0;
}

// The extent of a dataset; its rank, 0 when it is not there
static int h5_dims(hid_t file, const char* name, hsize_t dims[3]) {
    if (!h5_has(file, name)) return 0;
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return 0;
    const hid_t space = H5Dget_space(d);
    const int rank = H5Sget_simple_extent_ndims(space);
    int result = 0;
    if (rank >= 0 && rank <= 3) {
        H5Sget_simple_extent_dims(space, dims, NULL);
        result = rank == 0 ? -1 : rank;     // a scalar is there, of rank 0
    }
    H5Sclose(space);
    H5Dclose(d);
    return result;
}

// count doubles from the start of a dataset, of any stored type: HDF5 converts and decodes
static bool h5_read_doubles(hid_t file, const char* name, double* out, size_t count) {
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return false;
    const hid_t space = H5Dget_space(d);
    bool ok = (size_t)H5Sget_simple_extent_npoints(space) >= count;
    if (ok) {
        const int rank = H5Sget_simple_extent_ndims(space);
        hsize_t dims[3] = {0}, start[3] = {0}, cnt[3] = {0};
        H5Sget_simple_extent_dims(space, dims, NULL);
        if (rank > 0) {
            // A leading run of the dataset: whole rows of its last dimensions, as many as count fills
            size_t row = 1;
            for (int i = 1; i < rank; ++i) { row *= dims[i]; cnt[i] = dims[i]; }
            cnt[0] = count / row;
            ok = cnt[0] * row == count && H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, cnt, NULL) >= 0;
        }
        const hsize_t n = count;
        const hid_t mem = H5Screate_simple(1, &n, NULL);
        ok = ok && H5Dread(d, H5T_NATIVE_DOUBLE, mem, rank > 0 ? space : H5S_ALL, H5P_DEFAULT, out) >= 0;
        H5Sclose(mem);
    }
    H5Sclose(space);
    H5Dclose(d);
    return ok;
}

// One frame's rows [first, first + count) of a per atom variable, as floats, through HDF5. The lock is held.
static bool h5_read_frame(hid_t file, const char* name, int rank, size_t frame, size_t first, size_t count, float* out) {
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return false;
    const hid_t space = H5Dget_space(d);
    const hsize_t start3[3] = { frame, first, 0 }, cnt3[3] = { 1, count, 3 };
    const hsize_t start2[2] = { first, 0 }, cnt2[2] = { count, 3 };
    bool ok = H5Sselect_hyperslab(space, H5S_SELECT_SET, rank == 3 ? start3 : start2, NULL, rank == 3 ? cnt3 : cnt2, NULL) >= 0;
    const hsize_t n = count * 3;
    const hid_t mem = H5Screate_simple(1, &n, NULL);
    ok = ok && H5Dread(d, H5T_NATIVE_FLOAT, mem, space, H5P_DEFAULT, out) >= 0;
    H5Sclose(mem);
    H5Sclose(space);
    H5Dclose(d);
    return ok;
}

// ## The Amber convention

static bool nc4_has_convention(const char* conventions, const char* name) {
    const str_t c = str_from_cstr(conventions);
    size_t i = 0;
    while (i < c.len) {
        while (i < c.len && (c.ptr[i] == ',' || c.ptr[i] == ' ')) ++i;
        size_t j = i;
        while (j < c.len && c.ptr[j] != ',' && c.ptr[j] != ' ') ++j;
        if (j > i && str_eq_cstr((str_t){c.ptr + i, j - i}, name)) return true;
        i = j;
    }
    return false;
}

static md_amber_netcdf_kind_t nc4_kind(hid_t file) {
    char conventions[128];
    h5_attr_text(file, "Conventions", conventions, sizeof(conventions));
    if (nc4_has_convention(conventions, "AMBERRESTART")) return MD_AMBER_NETCDF_RESTART;
    if (nc4_has_convention(conventions, "AMBER"))        return MD_AMBER_NETCDF_TRAJECTORY;
    return MD_AMBER_NETCDF_HDF5;
}

// The per atom variables: their plain names, and cpptraj's integer compressed ones
static const char* const nc4_section_names[3][2] = {
    { "coordinates", "compressedpos" },
    { "velocities",  "compressedvel" },
    { "forces",      "compressedfrc" },
};

// Which name of the section the file has, NULL for none
static const char* nc4_section_name(hid_t file, int k) {
    for (int i = 0; i < 2; ++i) {
        if (h5_has(file, nc4_section_names[k][i])) return nc4_section_names[k][i];
    }
    return NULL;
}

// What a stored value of the variable is multiplied by
static double nc4_scale(hid_t file, const char* name) {
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return 1.0;
    double scale = 1.0, v;
    if (h5_attr_number(d, "scale_factor", &v) && v != 0.0) scale *= v;
    if (h5_attr_number(d, "icompressfac", &v) && v != 0.0) scale /= v;
    H5Dclose(d);
    return scale;
}

// The layout of a per atom variable of rank 3 {F, N, 3}, and its chunk table into table (from alloc). The lock is held.
static bool nc4_layout(nc4_layout_t* out, int64_t** table, hid_t file, const char* name, md_allocator_i* alloc) {
    MEMSET(out, 0, sizeof(*out));
    *table = NULL;
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return false;
    bool ok = false;
    const hid_t space = H5Dget_space(d);
    const hid_t type  = H5Dget_type(d);
    const hid_t dcpl  = H5Dget_create_plist(d);
    hsize_t dims[3] = {0};
    if (H5Sget_simple_extent_ndims(space) != 3) goto done;
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[2] != 3 || dims[1] == 0 || dims[1] > UINT32_MAX) goto done;

    const H5T_class_t cls = H5Tget_class(type);
    const size_t size = H5Tget_size(type);
    if      (cls == H5T_FLOAT   && size == 4) out->type = NC4_F32;
    else if (cls == H5T_FLOAT   && size == 8) out->type = NC4_F64;
    else if (cls == H5T_INTEGER && size == 4 && H5Tget_sign(type) == H5T_SGN_2) out->type = NC4_I32;
    else if (cls == H5T_INTEGER && size == 2 && H5Tget_sign(type) == H5T_SGN_2) out->type = NC4_I16;
    out->big_endian = H5Tget_order(type) == H5T_ORDER_BE;
    out->num_atoms = (uint32_t)dims[1];
    out->frame_bytes = dims[1] * 3 * (out->type ? nc4_type_size(out->type) : 4);
    out->scale = 1.0;
    out->mode = NC4_MODE_HDF5;  // until shown to be one decoded here
    ok = true;
    if (!out->type) goto done;  // HDF5 converts what this does not

    hsize_t userblock = 0;
    {
        const hid_t fcpl = H5Fget_create_plist(file);
        H5Pget_userblock(fcpl, &userblock);
        H5Pclose(fcpl);
    }

    const H5D_layout_t layout = H5Pget_layout(dcpl);
    if (layout == H5D_CONTIGUOUS) {
        const haddr_t addr = H5Dget_offset(d);
        if (addr != HADDR_UNDEF && H5Pget_nfilters(dcpl) == 0) {
            out->mode = NC4_MODE_CONTIGUOUS;
            out->base = userblock + addr;
        }
    } else if (layout == H5D_CHUNKED) {
        hsize_t chunk[3] = {0};
        if (H5Pget_chunk(dcpl, 3, chunk) != 3 || chunk[1] != dims[1] || chunk[2] != 3 || chunk[0] == 0) goto done;
        // Only the filters decoded here, in the order they are applied: shuffle, then deflate
        const int num_filters = H5Pget_nfilters(dcpl);
        bool known = true;
        for (int i = 0; i < num_filters && known; ++i) {
            unsigned int flags = 0, values[8];
            size_t num_values = ARRAY_SIZE(values);
            char filter_name[64];
            const H5Z_filter_t f = H5Pget_filter2(dcpl, (unsigned)i, &flags, &num_values, values, sizeof(filter_name), filter_name, NULL);
            if (f == H5Z_FILTER_SHUFFLE && !out->deflate) out->shuffle = 1;
            else if (f == H5Z_FILTER_DEFLATE) out->deflate = 1;
            else known = false;
        }
        if (!known) goto done;

        const size_t F = (size_t)dims[0];
        const size_t fpc = (size_t)chunk[0];
        const size_t C = (F + fpc - 1) / fpc;
        int64_t* t = md_alloc(alloc, MAX(C, 1) * 2 * sizeof(int64_t));
        for (size_t c = 0; c < C; ++c) {
            const hsize_t offset[3] = { c * fpc, 0, 0 };
            unsigned filter_mask = 0;
            haddr_t addr = HADDR_UNDEF;
            hsize_t bytes = 0;
            if (H5Dget_chunk_info_by_coord(d, offset, &filter_mask, &addr, &bytes) < 0 || addr == HADDR_UNDEF || filter_mask != 0) {
                goto done;  // a chunk not written, or one written without its filters: HDF5 knows what to do with it
            }
            t[c] = (int64_t)(userblock + addr);
            t[C + c] = (int64_t)bytes;
        }
        out->mode = NC4_MODE_CHUNKED;
        out->frames_per_chunk = (uint32_t)fpc;
        out->chunk_bytes = fpc * out->frame_bytes;
        out->num_chunks = (uint32_t)C;
        *table = t;
    }

done:
    H5Pclose(dcpl);
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(d);
    return ok;
}

// ## Decoding a frame

static void nc4_unshuffle(uint8_t* dst, const uint8_t* src, size_t bytes, size_t elem) {
    const size_t n = bytes / elem;
    for (size_t b = 0; b < elem; ++b) {
        const uint8_t* s = src + b * n;
        for (size_t i = 0; i < n; ++i) dst[i * elem + b] = s[i];
    }
    // What does not fill an element is left as it is
    MEMCPY(dst + n * elem, src + n * elem, bytes - n * elem);
}

static void nc4_convert(float* out, const uint8_t* raw, size_t count, uint32_t type, bool big_endian, double scale) {
    for (size_t i = 0; i < count; ++i) {
        double v = 0.0;
        switch (type) {
        case NC4_F32: { uint32_t u; MEMCPY(&u, raw + i * 4, 4); if (big_endian) u = BSWAP32(u); float f; MEMCPY(&f, &u, 4); v = f; break; }
        case NC4_F64: { uint64_t u; MEMCPY(&u, raw + i * 8, 8); if (big_endian) u = BSWAP64(u); MEMCPY(&v, &u, 8); break; }
        case NC4_I32: { uint32_t u; MEMCPY(&u, raw + i * 4, 4); if (big_endian) u = BSWAP32(u); int32_t s; MEMCPY(&s, &u, 4); v = s; break; }
        case NC4_I16: { uint16_t u; MEMCPY(&u, raw + i * 2, 2); if (big_endian) u = BSWAP16(u); int16_t s; MEMCPY(&s, &u, 2); v = s; break; }
        default: break;
        }
        out[i] = (float)(v * scale);
    }
}

static size_t nc4_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const nc4_source_t* src = (const nc4_source_t*)user_data;
    if (!src || !slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;
    const nc4_layout_t* l = &src->layout;
    const size_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    if (frame >= attr->format.shape[0] || N != l->num_atoms) return 0;
    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) return 0;
    const str_t path = str_from_cstr(src->path);

    if (l->mode == NC4_MODE_HDF5) {
        const md_hdf5_lock_t lock = md_hdf5_lock();
        bool ok = false;
        const hid_t file = H5Fopen(src->path, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file >= 0) {
            ok = h5_read_frame(file, src->dataset, 3, frame, first, count, (float*)dst);
            H5Fclose(file);
        }
        md_hdf5_unlock(lock);
        if (!ok) {
            MD_LOG_ERROR("Amber NetCDF-4: failed to read frame %zu of '%s' from '" STR_FMT "'", frame, src->dataset, STR_ARG(path));
            return 0;
        }
        if (l->scale != 1.0) {
            float* out = (float*)dst;
            for (size_t i = 0; i < cap; ++i) out[i] = (float)(out[i] * l->scale);
        }
        return cap;
    }

    const size_t elem = nc4_type_size(l->type);
    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;
    const uint8_t* values = NULL;
    if (l->mode == NC4_MODE_CONTIGUOUS) {
        const size_t bytes = count * 3 * elem;
        uint8_t* raw = md_temp_alloc(temp, bytes);
        if (raw && md_attribute_io_read_at(io, path, (int64_t)(l->base + frame * l->frame_bytes + first * 3 * elem), raw, bytes) == bytes) {
            values = raw;
        }
    } else {
        const size_t c = frame / l->frames_per_chunk;
        const size_t k = frame % l->frames_per_chunk;
        const int64_t offset = src->table[c];
        const size_t stored = (size_t)src->table[l->num_chunks + c];
        uint8_t* raw = md_temp_alloc(temp, stored);
        if (raw && md_attribute_io_read_at(io, path, offset, raw, stored) == stored) {
            uint8_t* chunk = raw;
            bool ok = true;
            if (l->deflate) {
                chunk = md_temp_alloc(temp, l->chunk_bytes);
                uLongf len = (uLongf)l->chunk_bytes;
                ok = chunk && uncompress(chunk, &len, raw, (uLong)stored) == Z_OK && len == l->chunk_bytes;
            } else {
                ok = stored == l->chunk_bytes;
            }
            if (ok && l->shuffle) {
                uint8_t* plain = md_temp_alloc(temp, l->chunk_bytes);
                nc4_unshuffle(plain, chunk, l->chunk_bytes, elem);
                chunk = plain;
            }
            if (ok) values = chunk + k * l->frame_bytes + first * 3 * elem;
        }
    }
    if (values) {
        nc4_convert((float*)dst, values, count * 3, l->type, l->big_endian, l->scale);
        written = cap;
    } else {
        MD_LOG_ERROR("Amber NetCDF-4: failed to read frame %zu of '%s' from '" STR_FMT "'", frame, src->dataset, STR_ARG(path));
    }
    md_temp_end(temp);
    return written;
}

// ## Entry points

static bool nc4_is_angstrom(hid_t file, const char* name) {
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return false;
    char units[64];
    const size_t len = h5_attr_text(d, "units", units, sizeof(units));
    H5Dclose(d);
    return len == 0 || str_eq_cstr_ignore_case(str_from_cstr(units), "angstrom") || str_eq_cstr_ignore_case(str_from_cstr(units), "angstroms");
}

md_amber_netcdf_kind_t md_amber_nc4_file_kind(str_t filename, size_t* out_num_atoms) {
    char path[NC4_PATH_MAX];
    str_copy_to_char_buf(path, sizeof(path), filename);
    if (out_num_atoms) *out_num_atoms = 0;
    const md_hdf5_lock_t lock = md_hdf5_lock();
    md_amber_netcdf_kind_t kind = MD_AMBER_NETCDF_NONE;
    const hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file >= 0) {
        kind = nc4_kind(file);
        const char* name = nc4_section_name(file, 0);
        hsize_t dims[3] = {0};
        const int rank = name ? h5_dims(file, name, dims) : 0;
        if (out_num_atoms && (kind == MD_AMBER_NETCDF_TRAJECTORY || kind == MD_AMBER_NETCDF_RESTART)) {
            *out_num_atoms = rank == 3 ? (size_t)dims[1] : rank == 2 ? (size_t)dims[0] : 0;
        }
        H5Fclose(file);
    }
    md_hdf5_unlock(lock);
    return kind;
}

// The cell of lengths (and angles, NULL for 90): zero lengths are no box
static md_unitcell_t nc4_cell(const double lengths[3], const double* angles) {
    if (!(lengths[0] > 0.0 && lengths[1] > 0.0 && lengths[2] > 0.0)) return md_unitcell_none();
    double a[3] = { 90.0, 90.0, 90.0 };
    if (angles && angles[0] > 0.0 && angles[0] < 180.0 && angles[1] > 0.0 && angles[1] < 180.0 && angles[2] > 0.0 && angles[2] < 180.0) {
        a[0] = angles[0]; a[1] = angles[1]; a[2] = angles[2];
    }
    return md_unitcell_from_extent_and_angles(lengths[0], lengths[1], lengths[2], a[0], a[1], a[2]);
}

bool md_amber_nc4_state_init_from_file(md_system_state_t* state, str_t filename, size_t num_atoms) {
    ASSERT(state && state->alloc);
    char path[NC4_PATH_MAX];
    str_copy_to_char_buf(path, sizeof(path), filename);

    md_temp_scope_t temp = md_temp_begin_avoid(state->alloc);
    bool ok = false;
    float* xyz = NULL;
    size_t N = 0;
    md_unitcell_t cell = md_unitcell_none();

    const md_hdf5_lock_t lock = md_hdf5_lock();
    const hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file >= 0) {
        const md_amber_netcdf_kind_t kind = nc4_kind(file);
        const char* name = nc4_section_name(file, 0);
        hsize_t dims[3] = {0};
        const int rank = name ? h5_dims(file, name, dims) : 0;
        N = rank == 3 ? (size_t)dims[1] : rank == 2 ? (size_t)dims[0] : 0;
        if (kind != MD_AMBER_NETCDF_TRAJECTORY && kind != MD_AMBER_NETCDF_RESTART) {
            MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' does not follow the Amber convention", STR_ARG(filename));
        } else if (!N || dims[rank - 1] != 3 || (rank == 3 && dims[0] == 0)) {
            MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' has no coordinates", STR_ARG(filename));
        } else if (num_atoms && N != num_atoms) {
            MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' holds %zu atoms, the system %zu", STR_ARG(filename), N, num_atoms);
        } else if (!nc4_is_angstrom(file, name)) {
            MD_LOG_ERROR("Amber NetCDF-4: the coordinates of '" STR_FMT "' are not in angstrom", STR_ARG(filename));
        } else {
            xyz = md_temp_alloc_array(temp, float, N * 3);
            ok = h5_read_frame(file, name, rank, 0, 0, N, xyz);
            const double scale = nc4_scale(file, name);
            for (size_t i = 0; ok && scale != 1.0 && i < N * 3; ++i) xyz[i] = (float)(xyz[i] * scale);
            if (ok && h5_has(file, "cell_lengths")) {
                double lengths[3], angles[3];
                const bool has_angles = h5_has(file, "cell_angles") && h5_read_doubles(file, "cell_angles", angles, 3);
                if (h5_read_doubles(file, "cell_lengths", lengths, 3)) cell = nc4_cell(lengths, has_angles ? angles : NULL);
            }
            if (!ok) MD_LOG_ERROR("Amber NetCDF-4: failed to read the coordinates of '" STR_FMT "'", STR_ARG(filename));
        }
        H5Fclose(file);
    } else {
        MD_LOG_ERROR("Amber NetCDF-4: cannot open '" STR_FMT "'", STR_ARG(filename));
    }
    md_hdf5_unlock(lock);

    if (ok) ok = md_system_state_init(state, N);
    if (ok) {
        MEMCPY(state->xyz, xyz, N * 3 * sizeof(float));
        state->num_atoms = N;
        state->unitcell = cell;
    }
    md_temp_end(temp);
    return ok;
}

static md_unit_t nc4_time_unit(const char* units, bool* known) {
    const str_t u = str_from_cstr(units);
    *known = true;
    if (str_empty(u) || str_eq_cstr_ignore_case(u, "picosecond") || str_eq_cstr_ignore_case(u, "picoseconds") || str_eq_cstr_ignore_case(u, "ps")) return md_unit_picosecond();
    if (str_eq_cstr_ignore_case(u, "nanosecond") || str_eq_cstr_ignore_case(u, "nanoseconds") || str_eq_cstr_ignore_case(u, "ns")) return md_unit_nanosecond();
    if (str_eq_cstr_ignore_case(u, "femtosecond") || str_eq_cstr_ignore_case(u, "femtoseconds") || str_eq_cstr_ignore_case(u, "fs")) return md_unit_femtosecond();
    *known = false;
    return md_unit_none();
}

bool md_amber_nc4_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    (void)flags;
    char path[NC4_PATH_MAX];
    const size_t path_len = md_path_write_canonical(path, sizeof(path), filename);
    if (!path_len) {
        MD_LOG_ERROR("Amber NetCDF-4: cannot resolve '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    const str_t spath = {path, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;

    // ## Everything the run needs from HDF5, at once and under the lock
    nc4_layout_t layout[3] = {0};
    int64_t* table[3] = {0};
    const char* names[3] = {0};
    size_t F = 0, N = 0;
    double* times = NULL;
    double* lengths = NULL;
    double* angles = NULL;
    double* temp0 = NULL;
    double* repidx = NULL;
    double* crdidx = NULL;
    char time_units[64] = "";
    bool ok = false;
    {
        const md_hdf5_lock_t lock = md_hdf5_lock();
        const hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0) {
            MD_LOG_ERROR("Amber NetCDF-4: cannot open '" STR_FMT "'", STR_ARG(spath));
        } else {
            const md_amber_netcdf_kind_t kind = nc4_kind(file);
            hsize_t dims[3] = {0};
            names[0] = nc4_section_name(file, 0);
            const int rank = names[0] ? h5_dims(file, names[0], dims) : 0;
            if (kind == MD_AMBER_NETCDF_RESTART) {
                MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' is a restart, a single set of coordinates rather than a trajectory", STR_ARG(spath));
            } else if (kind != MD_AMBER_NETCDF_TRAJECTORY) {
                MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' does not follow the Amber convention (no Conventions \"AMBER\")", STR_ARG(spath));
            } else if (rank != 3 || dims[0] == 0 || dims[1] == 0 || dims[2] != 3 || dims[0] > UINT32_MAX || dims[1] > UINT32_MAX) {
                MD_LOG_ERROR("Amber NetCDF-4: '" STR_FMT "' has no coordinates of {frame, atom, spatial}", STR_ARG(spath));
            } else if (!nc4_is_angstrom(file, names[0])) {
                MD_LOG_ERROR("Amber NetCDF-4: the coordinates of '" STR_FMT "' are not in angstrom", STR_ARG(spath));
            } else {
                F = (size_t)dims[0];
                N = (size_t)dims[1];
                ok = true;
                for (int k = 0; ok && k < 3; ++k) {
                    if (k > 0) {
                        names[k] = nc4_section_name(file, k);
                        hsize_t vd[3] = {0};
                        if (names[k] && (h5_dims(file, names[k], vd) != 3 || vd[0] != F || vd[1] != N || vd[2] != 3)) {
                            MD_LOG_INFO("Amber NetCDF-4: '%s' of '" STR_FMT "' is not of the frames and atoms of the coordinates, and is left out", names[k], STR_ARG(spath));
                            names[k] = NULL;
                        }
                    }
                    if (!names[k]) continue;
                    ok = nc4_layout(&layout[k], &table[k], file, names[k], arena);
                    layout[k].scale = nc4_scale(file, names[k]);
                }
                // The values per frame, through HDF5 whatever their encoding
                hsize_t d1[3];
                if (ok && h5_dims(file, "time", d1) == 1 && d1[0] == F) {
                    times = md_alloc(arena, F * sizeof(double));
                    if (!h5_read_doubles(file, "time", times, F)) times = NULL;
                    const hid_t t = H5Dopen2(file, "time", H5P_DEFAULT);
                    if (t >= 0) { h5_attr_text(t, "units", time_units, sizeof(time_units)); H5Dclose(t); }
                }
                if (ok && h5_dims(file, "cell_lengths", d1) == 2 && d1[0] == F && d1[1] == 3) {
                    lengths = md_alloc(arena, F * 3 * sizeof(double));
                    if (!h5_read_doubles(file, "cell_lengths", lengths, F * 3)) lengths = NULL;
                    if (lengths && h5_dims(file, "cell_angles", d1) == 2 && d1[0] == F && d1[1] == 3) {
                        angles = md_alloc(arena, F * 3 * sizeof(double));
                        if (!h5_read_doubles(file, "cell_angles", angles, F * 3)) angles = NULL;
                    }
                }
                double** series[3] = { &temp0, &repidx, &crdidx };
                const char* series_names[3] = { "temp0", "remd_repidx", "remd_crdidx" };
                for (int s = 0; ok && s < 3; ++s) {
                    if (h5_dims(file, series_names[s], d1) == 1 && d1[0] == F) {
                        *series[s] = md_alloc(arena, F * sizeof(double));
                        if (!h5_read_doubles(file, series_names[s], *series[s], F)) *series[s] = NULL;
                    }
                }
            }
            H5Fclose(file);
        }
        md_hdf5_unlock(lock);
    }
    if (!ok) goto done;

    // ## Time
    bool use_time = times != NULL;
    md_unit_t time_unit = md_unit_none();
    if (use_time) {
        bool known = false;
        time_unit = nc4_time_unit(time_units, &known);
        use_time = known;
        for (size_t f = 0; use_time && f < F; ++f) {
            if (!isfinite(times[f]) || (f > 0 && times[f] <= times[f - 1])) use_time = false;
        }
        if (!use_time) MD_LOG_INFO("Amber NetCDF-4: the time of '" STR_FMT "' does not increase from frame to frame, or is in units not known here; the frames are numbered instead", STR_ARG(spath));
    }
    if (!use_time) {
        times = md_alloc(arena, F * sizeof(double));
        for (size_t f = 0; f < F; ++f) times[f] = (double)f;
        time_unit = md_unit_none();
    }

    // ## Cells, zero without a box
    float* boxes = md_alloc(arena, F * 9 * sizeof(float));
    MEMSET(boxes, 0, F * 9 * sizeof(float));
    for (size_t f = 0; lengths && f < F; ++f) {
        const md_unitcell_t cell = nc4_cell(lengths + f * 3, angles ? angles + f * 3 : NULL);
        md_unitcell_A_extract_float((float(*)[3])(boxes + f * 9), &cell);
    }

    // ## The providers' state, one per variable, owned by its attribute
    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;
    md_attribute_virtual_t virt[3] = {0};
    for (int k = 0; k < 3; ++k) {
        if (!names[k]) continue;
        const size_t size = sizeof(nc4_source_t) + 2 * (size_t)layout[k].num_chunks * sizeof(int64_t);
        nc4_source_t* src = md_attributes_alloc_user_data(&sys->attributes, size);
        if (!src) goto done;
        MEMSET(src, 0, size);
        src->layout = layout[k];
        str_copy_to_char_buf(src->dataset, sizeof(src->dataset), str_from_cstr(names[k]));
        str_copy_to_char_buf(src->path, sizeof(src->path), spath);
        if (layout[k].num_chunks) MEMCPY(src->table, table[k], 2 * (size_t)layout[k].num_chunks * sizeof(int64_t));
        virt[k] = (md_attribute_virtual_t){ .provider = nc4_provider, .user_data = src, .user_data_size = size };
    }

    // The frame table, for what reads the run's file as a whole: where each frame's coordinates are
    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    for (size_t f = 0; f < F; ++f) {
        const nc4_layout_t* l = &layout[0];
        if (l->mode == NC4_MODE_CHUNKED) {
            const size_t c = f / l->frames_per_chunk;
            offsets[f] = table[0][c];
            sizes[f]   = table[0][l->num_chunks + c];
        } else if (l->mode == NC4_MODE_CONTIGUOUS) {
            offsets[f] = (int64_t)(l->base + f * l->frame_bytes);
            sizes[f]   = (int64_t)l->frame_bytes;
        } else {
            offsets[f] = -1;
            sizes[f]   = 0;
        }
    }

    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = times,
        .time_unit     = time_unit,
        .unitcell      = boxes,
        .source_path   = spath,
        .source_offset = offsets,
        .source_size   = sizes,
        .position_virt = &virt[0],
    };
    if (!md_run_publish(sys, run, &desc)) {
        for (int k = 0; k < 3; ++k) if (virt[k].user_data) md_free(sys->attributes.alloc, virt[k].user_data, virt[k].user_data_size);
        goto done;
    }

    md_attributes_t* attributes = &sys->attributes;
    char buf[512];
    const md_unit_t kcal_per_mol = md_unit_div(md_unit_scl(md_unit_calorie(), 1.0e3), md_unit_mole());
    const struct { int k; const char* leaf; md_unit_t unit; str_t label; } per_atom[] = {
        { 1, "atom/velocity", md_unit_div(md_unit_angstrom(), md_unit_picosecond()), STR_INIT("Velocity") },
        { 2, "atom/force",    md_unit_div(kcal_per_mol, md_unit_angstrom()),         STR_INIT("Force") },
    };
    bool published = true;
    for (size_t i = 0; published && i < ARRAY_SIZE(per_atom); ++i) {
        const int k = per_atom[i].k;
        if (!names[k]) continue;
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, str_from_cstr(per_atom[i].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)F, (uint32_t)N } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = per_atom[i].unit, .label = per_atom[i].label,
            .virt = &virt[k]}) != MD_ATTRIBUTE_INVALID;
    }
    if (published && temp0) {
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, STR_LIT("remd/temp0")),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_kelvin(), .label = STR_INIT("Target temperature"),
            .data = temp0, .byte_size = F * sizeof(double)}) != MD_ATTRIBUTE_INVALID;
    }
    const struct { const double* v; const char* leaf; str_t label; } indices[] = {
        { repidx, "remd/repidx", STR_INIT("Replica index") },
        { crdidx, "remd/crdidx", STR_INIT("Coordinate index") },
    };
    for (size_t i = 0; published && i < ARRAY_SIZE(indices); ++i) {
        if (!indices[i].v) continue;
        int32_t* ints = md_alloc(arena, F * sizeof(int32_t));
        for (size_t f = 0; f < F; ++f) ints[f] = (int32_t)indices[i].v[f];
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, str_from_cstr(indices[i].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .label = indices[i].label,
            .data = ints, .byte_size = F * sizeof(int32_t)}) != MD_ATTRIBUTE_INVALID;
    }
    if (!published) {
        MD_LOG_ERROR("Amber NetCDF-4: failed to publish '" STR_FMT "' as '" STR_FMT "'", STR_ARG(spath), STR_ARG(run));
        md_attributes_remove_prefix(attributes, run);
        goto done;
    }
    if (layout[0].mode == NC4_MODE_HDF5) {
        MD_LOG_INFO("Amber NetCDF-4: the coordinates of '" STR_FMT "' are encoded in a way read through HDF5, one frame at a time", STR_ARG(spath));
    }
    result = true;

done:
    md_arena_allocator_destroy(arena);
    return result;
}

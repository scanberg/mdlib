#include <md_h5md.h>

#include <md_system.h>
#include <md_tpr.h>
#include <md_types.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>
#include <core/md_vec_math.h>

#include <hdf5.h>

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define H5MD_PATH_MAX 512

// ### LOCK ###
// HDF5 is not thread safe unless it was built to be, and a provider may be entered from any number
// of threads. Every HDF5 call this reader makes is therefore made holding this lock. Initialised on
// first use, as md_allocator.c does its thread key: the first call into this reader must not race
// another first call, which a loader never does.
static md_mutex_t h5md_mutex;
static bool       h5md_mutex_ready = false;

typedef struct h5md_lock_t {
    H5E_auto2_t func;
    void*       client_data;
} h5md_lock_t;

// Takes the lock and silences HDF5's own error printing for as long as it is held: it prints its
// stack to stderr on every failed probe, and probing for optional objects is most of what H5MD
// reading is.
static h5md_lock_t h5md_lock(void) {
    if (!h5md_mutex_ready) {
        md_mutex_init(&h5md_mutex);
        h5md_mutex_ready = true;
    }
    md_mutex_lock(&h5md_mutex);
    h5md_lock_t lock = {0};
    H5Eget_auto2(H5E_DEFAULT, &lock.func, &lock.client_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    return lock;
}

static void h5md_unlock(h5md_lock_t lock) {
    H5Eset_auto2(H5E_DEFAULT, lock.func, lock.client_data);
    md_mutex_unlock(&h5md_mutex);
}

// ### HDF5 HELPERS ###

// What is at path: H5I_GROUP, H5I_DATASET, or H5I_BADID for nothing. H5Lexists fails rather than
// answering false for a path through a group that is not there, so the object is opened instead.
static H5I_type_t h5_kind(hid_t loc, const char* path) {
    hid_t obj = H5Oopen(loc, path, H5P_DEFAULT);
    if (obj < 0) {
        return H5I_BADID;
    }
    const H5I_type_t kind = H5Iget_type(obj);
    H5Oclose(obj);
    return kind;
}

// The rank of a dataset and the extent of each axis, -1 when it is not a simple dataspace.
static int h5_shape(hsize_t dims[H5S_MAX_RANK], hid_t dset) {
    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        return -1;
    }
    int rank = -1;
    const H5S_class_t cls = H5Sget_simple_extent_type(space);
    if (cls == H5S_SCALAR) {
        rank = 0;
    } else if (cls == H5S_SIMPLE) {
        rank = H5Sget_simple_extent_ndims(space);
        if (rank < 0 || rank > H5S_MAX_RANK || H5Sget_simple_extent_dims(space, dims, NULL) != rank) {
            rank = -1;
        }
    }
    H5Sclose(space);
    return rank;
}

static size_t h5_product(const hsize_t dims[], int beg, int end) {
    size_t n = 1;
    for (int i = beg; i < end; ++i) {
        n *= (size_t)dims[i];
    }
    return n;
}

// Row 'row' of the outermost axis of a dataset, all of it when row is negative, as mem_type.
// count is the number of elements that makes: what dst has room for, and a mismatch is a failure.
static bool h5_read_rows(void* dst, hid_t dset, hid_t mem_type, int64_t row, size_t count) {
    hsize_t dims[H5S_MAX_RANK];
    const int rank = h5_shape(dims, dset);
    if (rank < 0) {
        return false;
    }
    if (row < 0) {
        if (h5_product(dims, 0, rank) != count) {
            return false;
        }
        return count == 0 || H5Dread(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst) >= 0;
    }
    if (rank == 0 || (hsize_t)row >= dims[0] || h5_product(dims, 1, rank) != count) {
        return false;
    }
    hsize_t start[H5S_MAX_RANK] = {0};
    hsize_t extent[H5S_MAX_RANK];
    start[0] = (hsize_t)row;
    extent[0] = 1;
    for (int i = 1; i < rank; ++i) {
        extent[i] = dims[i];
    }
    const hsize_t mem_dims = (hsize_t)count;
    hid_t file_space = H5Dget_space(dset);
    hid_t mem_space  = H5Screate_simple(1, &mem_dims, NULL);
    const bool ok = file_space >= 0 && mem_space >= 0 &&
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, extent, NULL) >= 0 &&
        H5Dread(dset, mem_type, mem_space, file_space, H5P_DEFAULT, dst) >= 0;
    if (mem_space >= 0) H5Sclose(mem_space);
    if (file_space >= 0) H5Sclose(file_space);
    return ok;
}

static bool h5_read_rows_at(void* dst, hid_t loc, const char* path, hid_t mem_type, int64_t row, size_t count) {
    hid_t dset = H5Dopen(loc, path, H5P_DEFAULT);
    if (dset < 0) {
        return false;
    }
    const bool ok = h5_read_rows(dst, dset, mem_type, row, count);
    H5Dclose(dset);
    return ok;
}

// A numeric attribute, or its first element, as mem_type
static bool h5_attr_read(void* dst, hid_t loc, const char* name, hid_t mem_type) {
    if (H5Aexists(loc, name) <= 0) {
        return false;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        return false;
    }
    bool ok = false;
    hid_t space = H5Aget_space(attr);
    if (space >= 0) {
        const hssize_t n = H5Sget_simple_extent_npoints(space);
        if (n == 1) {
            ok = H5Aread(attr, mem_type, dst) >= 0;
        } else if (n > 1) {
            md_temp_scope_t temp = md_temp_begin();
            void* buf = md_temp_alloc(temp, (size_t)n * H5Tget_size(mem_type));
            if (buf && H5Aread(attr, mem_type, buf) >= 0) {
                MEMCPY(dst, buf, H5Tget_size(mem_type));
                ok = true;
            }
            md_temp_end(temp);
        }
        H5Sclose(space);
    }
    H5Aclose(attr);
    return ok;
}

// Every element of an integer attribute, at most cap of them. Returns how many it has.
static size_t h5_attr_i32s(int32_t* dst, size_t cap, hid_t loc, const char* name) {
    if (H5Aexists(loc, name) <= 0) {
        return 0;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        return 0;
    }
    size_t count = 0;
    hid_t space = H5Aget_space(attr);
    const hssize_t n = space >= 0 ? H5Sget_simple_extent_npoints(space) : -1;
    if (n > 0) {
        md_temp_scope_t temp = md_temp_begin();
        int32_t* buf = md_temp_alloc_array(temp, int32_t, (size_t)n);
        if (buf && H5Aread(attr, H5T_NATIVE_INT32, buf) >= 0) {
            count = (size_t)n;
            if (cap) MEMCPY(dst, buf, MIN(count, cap) * sizeof(int32_t));
        }
        md_temp_end(temp);
    }
    if (space >= 0) H5Sclose(space);
    H5Aclose(attr);
    return count;
}

// Every string of an attribute (is_attr) or a dataset, copied into alloc, at most cap of them.
// Fixed and variable length alike; fixed length strings lose their padding. Returns how many
// there are, 0 for anything that is not text.
static size_t h5_read_strs(str_t* dst, size_t cap, hid_t obj, bool is_attr, md_allocator_i* alloc) {
    hid_t type  = is_attr ? H5Aget_type(obj)  : H5Dget_type(obj);
    hid_t space = is_attr ? H5Aget_space(obj) : H5Dget_space(obj);
    size_t count = 0;
    // The strings are copied into alloc while this scope is open, so it must not be alloc's arena
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    if (type >= 0 && space >= 0 && H5Tget_class(type) == H5T_STRING) {
        const hssize_t n = H5Sget_simple_extent_npoints(space);
        if (n > 0) {
            if (H5Tis_variable_str(type) > 0) {
                hid_t mem = H5Tcopy(H5T_C_S1);
                H5Tset_size(mem, H5T_VARIABLE);
                H5Tset_cset(mem, H5Tget_cset(type));
                char** buf = md_temp_alloc_array(temp, char*, (size_t)n);
                if (buf) MEMSET(buf, 0, (size_t)n * sizeof(char*));
                const herr_t err = !buf ? -1 : is_attr ? H5Aread(obj, mem, buf) : H5Dread(obj, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
                if (err >= 0) {
                    count = (size_t)n;
                    for (size_t i = 0; i < (size_t)n; ++i) {
                        if (i < cap) dst[i] = buf[i] ? str_copy(str_from_cstr(buf[i]), alloc) : (str_t){0};
                        if (buf[i]) H5free_memory(buf[i]);
                    }
                }
                H5Tclose(mem);
            } else {
                const size_t size = H5Tget_size(type);
                char* buf = md_temp_alloc(temp, (size_t)n * size);
                const herr_t err = !buf ? -1 : is_attr ? H5Aread(obj, type, buf) : H5Dread(obj, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
                if (err >= 0) {
                    count = (size_t)n;
                    for (size_t i = 0; i < MIN((size_t)n, cap); ++i) {
                        const char* s = buf + i * size;
                        size_t len = 0;
                        while (len < size && s[len] != '\0') ++len;
                        while (len > 0 && s[len - 1] == ' ') --len;
                        dst[i] = str_copy((str_t){s, len}, alloc);
                    }
                }
            }
        }
    }
    md_temp_end(temp);
    if (space >= 0) H5Sclose(space);
    if (type >= 0) H5Tclose(type);
    return count;
}

static size_t h5_attr_strs(str_t* dst, size_t cap, hid_t loc, const char* name, md_allocator_i* alloc) {
    if (H5Aexists(loc, name) <= 0) {
        return 0;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr < 0) {
        return 0;
    }
    const size_t count = h5_read_strs(dst, cap, attr, true, alloc);
    H5Aclose(attr);
    return count;
}

static str_t h5_attr_str(hid_t loc, const char* name, md_allocator_i* alloc) {
    str_t s = {0};
    h5_attr_strs(&s, 1, loc, name, alloc);
    return s;
}

// The names of the links in a group, in name order, copied into alloc
static size_t h5_children(md_array(str_t)* out, hid_t loc, const char* path, md_allocator_i* alloc) {
    hid_t group = H5Gopen(loc, path, H5P_DEFAULT);
    if (group < 0) {
        return 0;
    }
    H5G_info_t info;
    if (H5Gget_info(group, &info) >= 0) {
        for (hsize_t i = 0; i < info.nlinks; ++i) {
            char name[H5MD_PATH_MAX];
            const ssize_t len = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, name, sizeof(name), H5P_DEFAULT);
            if (len > 0 && (size_t)len < sizeof(name)) {
                md_array_push(*out, str_copy((str_t){name, (size_t)len}, alloc), alloc);
            }
        }
    }
    H5Gclose(group);
    return md_array_size(*out);
}

static bool h5_path(char* buf, size_t cap, const char* parent, str_t child) {
    const int len = snprintf(buf, cap, "%s/" STR_FMT, parent, STR_ARG(child));
    return len > 0 && (size_t)len < cap;
}

// ### UNITS ###

// A unit as the H5MD units module writes it: factors separated by spaces, each a number or a
// symbol, raised to a power written as a signed integer straight after it - "nm+3", "kJ mol-1
// nm-1", "10+3 m". The symbols are SI and parsed by md_unit_parse, so "A" is an ampere.
static bool h5md_unit_parse(md_unit_t* out, str_t str) {
    md_unit_t unit = md_unit_none();
    str = str_trim(str);
    while (!str_empty(str)) {
        size_t end = 0;
        while (end < str.len && str.ptr[end] != ' ') ++end;
        str_t tok = str_substr(str, 0, end);
        str = str_trim(str_substr(str, end, SIZE_MAX));

        int power = 1;
        size_t k = tok.len;
        while (k > 0 && is_digit(tok.ptr[k - 1])) --k;
        if (k > 1 && k < tok.len && (tok.ptr[k - 1] == '+' || tok.ptr[k - 1] == '-')) {
            power = 0;
            for (size_t i = k; i < tok.len; ++i) power = power * 10 + (tok.ptr[i] - '0');
            if (tok.ptr[k - 1] == '-') power = -power;
            tok = str_substr(tok, 0, k - 1);
        }

        md_unit_t factor;
        if (is_digit(tok.ptr[0]) || tok.ptr[0] == '.') {
            char num[64];
            str_copy_to_char_buf(num, sizeof(num), tok);
            char* num_end = NULL;
            const double value = strtod(num, &num_end);
            if (!num_end || *num_end != '\0' || value == 0.0) {
                return false;
            }
            factor = md_unit_scl(md_unit_none(), pow(value, power));
        } else {
            if (!md_unit_parse(&factor, tok)) {
                return false;
            }
            factor = md_unit_pow(factor, power);
        }
        unit = md_unit_mul(unit, factor);
    }
    *out = unit;
    return true;
}

// The 'unit' attribute of a dataset. False, and none, when it has none or it cannot be read.
static bool h5md_unit_of(md_unit_t* out, hid_t loc, const char* path, md_allocator_i* alloc) {
    *out = md_unit_none();
    hid_t obj = H5Oopen(loc, path, H5P_DEFAULT);
    if (obj < 0) {
        return false;
    }
    const str_t text = h5_attr_str(obj, "unit", alloc);
    H5Oclose(obj);
    if (str_empty(text)) {
        return false;
    }
    if (!h5md_unit_parse(out, text)) {
        MD_LOG_INFO("H5MD: the unit '" STR_FMT "' of '%s' is not one this reader knows; it is read without a unit", STR_ARG(text), path);
        *out = md_unit_none();
        return false;
    }
    return true;
}

// ### ELEMENTS ###

// An H5MD element: a dataset when it is time independent, a group of value, step and optionally
// time when it is not.
typedef struct h5md_element_t {
    char      value[H5MD_PATH_MAX]; // the dataset holding the values: the element itself, or <element>/value
    bool      temporal;
    bool      has_unit;
    int       rank;                 // of the value dataset, the frame axis included
    hsize_t   dims[H5S_MAX_RANK];
    md_unit_t unit;

    // Time dependent only
    size_t    num_frames;
    int64_t*  step;                 // num_frames
    double*   time;                 // num_frames, NULL when the element has no time
    md_unit_t time_unit;
} h5md_element_t;

// One of the two axes of a time dependent element, as H5MD stores them either way: a dataset with a
// value per frame, or a scalar increment whose first value is its 'offset' attribute (zero when it
// has none). mem_type is H5T_NATIVE_INT64 or H5T_NATIVE_DOUBLE, and dst has num_frames of them.
static bool h5md_read_axis(void* dst, hid_t mem_type, size_t num_frames, hid_t group, const char* name) {
    hid_t dset = H5Dopen(group, name, H5P_DEFAULT);
    if (dset < 0) {
        return false;
    }
    bool ok = false;
    hsize_t dims[H5S_MAX_RANK];
    const int rank = h5_shape(dims, dset);
    const bool is_int = H5Tequal(mem_type, H5T_NATIVE_INT64) > 0;
    if (rank == 1) {
        ok = dims[0] == num_frames && h5_read_rows(dst, dset, mem_type, -1, num_frames);
    } else if (rank == 0) {
        if (is_int) {
            int64_t inc = 0, offset = 0;
            ok = h5_read_rows(&inc, dset, mem_type, -1, 1);
            h5_attr_read(&offset, dset, "offset", mem_type);
            for (size_t i = 0; ok && i < num_frames; ++i) ((int64_t*)dst)[i] = offset + (int64_t)i * inc;
        } else {
            double inc = 0, offset = 0;
            ok = h5_read_rows(&inc, dset, mem_type, -1, 1);
            h5_attr_read(&offset, dset, "offset", mem_type);
            for (size_t i = 0; ok && i < num_frames; ++i) ((double*)dst)[i] = offset + (double)i * inc;
        }
    }
    H5Dclose(dset);
    return ok;
}

// Opens the element at path. False when there is none; logged when there is something there that
// is not an element.
static bool h5md_element_open(h5md_element_t* e, hid_t file, const char* path, md_allocator_i* alloc) {
    MEMSET(e, 0, sizeof(*e));
    const H5I_type_t kind = h5_kind(file, path);
    if (kind == H5I_DATASET) {
        snprintf(e->value, sizeof(e->value), "%s", path);
    } else if (kind == H5I_GROUP) {
        e->temporal = true;
        const int len = snprintf(e->value, sizeof(e->value), "%s/value", path);
        if (len <= 0 || (size_t)len >= sizeof(e->value) || h5_kind(file, e->value) != H5I_DATASET) {
            return false;   // a group, but not an element: a subgroup of observables, say
        }
    } else {
        return false;
    }

    hid_t dset = H5Dopen(file, e->value, H5P_DEFAULT);
    if (dset < 0) {
        return false;
    }
    e->rank = h5_shape(e->dims, dset);
    H5Dclose(dset);
    if (e->rank < 0 || (e->temporal && e->rank == 0)) {
        MD_LOG_ERROR("H5MD: '%s' does not have the shape of an element", path);
        return false;
    }
    e->has_unit = h5md_unit_of(&e->unit, file, e->value, alloc);

    if (e->temporal) {
        e->num_frames = (size_t)e->dims[0];
        hid_t group = H5Gopen(file, path, H5P_DEFAULT);
        if (group < 0) {
            return false;
        }
        bool ok = true;
        if (e->num_frames > 0) {
            e->step = md_alloc(alloc, e->num_frames * sizeof(int64_t));
            ok = h5md_read_axis(e->step, H5T_NATIVE_INT64, e->num_frames, group, "step");
            if (!ok) {
                MD_LOG_ERROR("H5MD: '%s' has no step for each of its %zu frames", path, e->num_frames);
            } else if (h5_kind(group, "time") == H5I_DATASET) {
                e->time = md_alloc(alloc, e->num_frames * sizeof(double));
                ok = h5md_read_axis(e->time, H5T_NATIVE_DOUBLE, e->num_frames, group, "time");
                if (!ok) {
                    MD_LOG_ERROR("H5MD: '%s' has a time that does not match its %zu frames", path, e->num_frames);
                } else {
                    char time_path[H5MD_PATH_MAX];
                    snprintf(time_path, sizeof(time_path), "%s/time", path);
                    h5md_unit_of(&e->time_unit, file, time_path, alloc);
                }
            }
        }
        H5Gclose(group);
        return ok;
    }
    return true;
}

// The factor taking an element's values to Angstrom. An element without a unit is taken to be in
// Angstrom already, which is what a file written without the units module most likely means.
static bool h5md_length_scale(double* out, const h5md_element_t* e, const char* what) {
    *out = 1.0;
    if (e->has_unit && !md_unit_conversion_factor(out, e->unit, md_unit_angstrom())) {
        MD_LOG_ERROR("H5MD: the %s is not in a unit of length", what);
        return false;
    }
    return true;
}

// ### PARTICLES ###

// What is read of the particle group
typedef struct h5md_particles_t {
    char   path[H5MD_PATH_MAX];     // "/particles/<name>"
    size_t num_particles;
    size_t num_groups;              // particle groups in the file

    h5md_element_t position;        // num_particles == 0 when there is no position
    bool           has_position;

    bool           periodic[3];
    h5md_element_t edges;
    bool           has_edges;
} h5md_particles_t;

// The particles in a group: the extent of the particle axis of position, or of whichever standard
// element it has instead. Only the shape is looked at.
static size_t h5md_group_size(hid_t file, const char* group) {
    static const char* names[] = { "position", "species", "mass", "charge", "id", "velocity", "force" };
    for (size_t i = 0; i < ARRAY_SIZE(names); ++i) {
        char path[H5MD_PATH_MAX];
        snprintf(path, sizeof(path), "%s/%s", group, names[i]);
        const H5I_type_t kind = h5_kind(file, path);
        if (kind == H5I_GROUP) {
            snprintf(path, sizeof(path), "%s/%s/value", group, names[i]);
        } else if (kind != H5I_DATASET) {
            continue;
        }
        hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
        if (dset < 0) continue;
        hsize_t dims[H5S_MAX_RANK];
        const int rank = h5_shape(dims, dset);
        H5Dclose(dset);
        const int axis = kind == H5I_GROUP ? 1 : 0;
        if (rank > axis) {
            return (size_t)dims[axis];
        }
    }
    return 0;
}

// Opens the file's H5MD structure: the metadata group, which is what makes an HDF5 file an H5MD
// one, and the particle group read - the one with the most particles, the first by name of equals.
static bool h5md_open(h5md_particles_t* p, hid_t file, md_allocator_i* alloc) {
    MEMSET(p, 0, sizeof(*p));

    int32_t version[2] = {0};
    hid_t meta = H5Gopen(file, "h5md", H5P_DEFAULT);
    const size_t num_version = meta >= 0 ? h5_attr_i32s(version, 2, meta, "version") : 0;
    if (meta >= 0) H5Gclose(meta);
    if (num_version != 2) {
        MD_LOG_ERROR("H5MD: the file has no /h5md group with a version, which every H5MD file has at its root");
        return false;
    }
    if (version[0] != 1) {
        MD_LOG_ERROR("H5MD: the file is H5MD %d.%d, and only 1.x is read", version[0], version[1]);
        return false;
    }

    md_array(str_t) groups = 0;
    h5_children(&groups, file, "particles", alloc);
    for (size_t i = 0; i < md_array_size(groups); ++i) {
        char path[H5MD_PATH_MAX];
        if (!h5_path(path, sizeof(path), "/particles", groups[i]) || h5_kind(file, path) != H5I_GROUP) {
            continue;
        }
        p->num_groups += 1;
        const size_t n = h5md_group_size(file, path);
        if (n > p->num_particles) {
            p->num_particles = n;
            snprintf(p->path, sizeof(p->path), "%s", path);
        }
    }
    if (p->num_particles == 0) {
        MD_LOG_ERROR("H5MD: the file has no particles");
        return false;
    }

    char path[H5MD_PATH_MAX];
    snprintf(path, sizeof(path), "%s/id", p->path);
    h5md_element_t id;
    if (h5md_element_open(&id, file, path, alloc) && id.temporal) {
        MD_LOG_ERROR("H5MD: the particles of '%s' change between frames (a time dependent id), which is not read", p->path);
        return false;
    }

    snprintf(path, sizeof(path), "%s/position", p->path);
    p->has_position = h5md_element_open(&p->position, file, path, alloc);
    if (p->has_position) {
        const int axis = p->position.temporal ? 1 : 0;
        if (p->position.rank != axis + 2 || p->position.dims[axis] != p->num_particles || p->position.dims[axis + 1] != 3) {
            MD_LOG_ERROR("H5MD: '%s' is not a position per particle in three dimensions", path);
            return false;
        }
        if (p->position.temporal && p->position.num_frames == 0) {
            MD_LOG_ERROR("H5MD: '%s' has no frames", path);
            return false;
        }
    }

    // The box: a periodic dimension has a box vector, a dimension without a boundary does not
    snprintf(path, sizeof(path), "%s/box", p->path);
    hid_t box = H5Gopen(file, path, H5P_DEFAULT);
    if (box >= 0) {
        int32_t dimension = 3;
        h5_attr_read(&dimension, box, "dimension", H5T_NATIVE_INT32);
        str_t boundary[3] = {0};
        const size_t num_boundary = h5_attr_strs(boundary, 3, box, "boundary", alloc);
        H5Gclose(box);
        if (dimension != 3 || num_boundary != 3) {
            MD_LOG_ERROR("H5MD: the box of '%s' is not three dimensional", p->path);
            return false;
        }
        for (int k = 0; k < 3; ++k) {
            p->periodic[k] = str_eq_cstr(boundary[k], "periodic");
        }
        snprintf(path, sizeof(path), "%s/box/edges", p->path);
        p->has_edges = h5md_element_open(&p->edges, file, path, alloc);
        if (p->has_edges) {
            const int tensor = p->edges.rank - (p->edges.temporal ? 1 : 0);
            const hsize_t* d = p->edges.dims + (p->edges.temporal ? 1 : 0);
            if (!((tensor == 1 && d[0] == 3) || (tensor == 2 && d[0] == 3 && d[1] == 3))) {
                MD_LOG_ERROR("H5MD: the box edges of '%s' are neither a vector nor a 3x3 matrix", p->path);
                return false;
            }
        }
    }
    return true;
}

// The cell at every position frame (one cell for a position without frames), box vectors as rows,
// Angstrom: from the edges row at the same step, or the time independent edges throughout. A frame
// with no row at its step has no cell, and so does a dimension that is not periodic - its row is
// zero, which is what makes it not periodic in the cell.
static bool h5md_read_cells(float* out, hid_t file, const h5md_particles_t* p, md_allocator_i* alloc) {
    const size_t F = p->position.temporal ? p->position.num_frames : 1;
    MEMSET(out, 0, F * 9 * sizeof(float));
    if (!p->has_edges || !(p->periodic[0] || p->periodic[1] || p->periodic[2])) {
        return true;
    }
    double scale;
    if (!h5md_length_scale(&scale, &p->edges, "box")) {
        return false;
    }

    const bool   matrix = p->edges.rank - (p->edges.temporal ? 1 : 0) == 2;
    const size_t width  = matrix ? 9 : 3;
    const size_t rows   = p->edges.temporal ? p->edges.num_frames : 1;
    double* edges = md_alloc(alloc, MAX(rows, 1) * width * sizeof(double));
    if (rows > 0 && !h5_read_rows_at(edges, file, p->edges.value, H5T_NATIVE_DOUBLE, -1, rows * width)) {
        MD_LOG_ERROR("H5MD: failed to read the box of '%s'", p->path);
        return false;
    }

    size_t missing = 0;
    size_t row = 0;
    for (size_t f = 0; f < F; ++f) {
        if (p->edges.temporal) {
            // Without frames of its own the position is matched to the first row
            const int64_t step = p->position.temporal ? p->position.step[f] : (rows ? p->edges.step[0] : 0);
            // Both are in step order, so the search only ever moves forward
            while (row < rows && p->edges.step[row] < step) ++row;
            if (row == rows || p->edges.step[row] != step) {
                missing += 1;
                continue;
            }
        }
        const double* m = edges + row * width;
        float* cell = out + f * 9;
        for (int i = 0; i < 3; ++i) {
            if (!p->periodic[i]) continue;
            for (int j = 0; j < 3; ++j) {
                cell[i * 3 + j] = (float)((matrix ? m[i * 3 + j] : (i == j ? m[i] : 0.0)) * scale);
            }
        }
    }
    if (missing) {
        MD_LOG_INFO("H5MD: %zu frames of '%s' have no box at their step", missing, p->path);
    }
    return true;
}

// ### SYSTEM ###

// The index of each particle id in the group, for tuples that name particles by id. NULL, and
// tuples naming particles by index, when the group has no id.
static bool h5md_read_ids(int64_t** out, hid_t file, const h5md_particles_t* p, md_allocator_i* alloc) {
    *out = NULL;
    char path[H5MD_PATH_MAX];
    snprintf(path, sizeof(path), "%s/id", p->path);
    if (h5_kind(file, path) != H5I_DATASET) {
        return true;
    }
    int64_t* ids = md_alloc(alloc, p->num_particles * sizeof(int64_t));
    if (!h5_read_rows_at(ids, file, path, H5T_NATIVE_INT64, -1, p->num_particles)) {
        MD_LOG_ERROR("H5MD: '%s' is not one integer per particle", path);
        return false;
    }
    *out = ids;
    return true;
}

// Whether a tuple list belongs to the particle group read: it names the group in its particles_group
// attribute, or it has none and the group is the only one.
static bool h5md_tuples_belong(hid_t dset, const h5md_particles_t* p) {
    if (H5Aexists(dset, "particles_group") <= 0) {
        return p->num_groups == 1;
    }
    bool result = false;
    hid_t group = -1;
    hid_t attr = H5Aopen(dset, "particles_group", H5P_DEFAULT);
#if H5_VERSION_GE(1, 12, 0)
    // The references API of 1.12, which reads the object references H5MD writes as well
    H5R_ref_t ref;
    if (attr >= 0 && H5Aread(attr, H5T_STD_REF, &ref) >= 0) {
        group = H5Ropen_object(&ref, H5P_DEFAULT, H5P_DEFAULT);
        H5Rdestroy(&ref);
    }
#else
    hobj_ref_t ref;
    if (attr >= 0 && H5Aread(attr, H5T_STD_REF_OBJ, &ref) >= 0) {
        group = H5Rdereference2(dset, H5P_DEFAULT, H5R_OBJECT, &ref);
    }
#endif
    if (group >= 0) {
        char name[H5MD_PATH_MAX];
        const ssize_t len = H5Iget_name(group, name, sizeof(name));
        result = len > 0 && strcmp(name, p->path) == 0;
        H5Oclose(group);
    }
    if (attr >= 0) H5Aclose(attr);
    return result;
}

// The pairs in /connectivity that belong to the group read, as particle indices, appended to out.
// Pairs holding the dataset's fill value, or an index or id outside the group, are left out.
static void h5md_read_bonds(md_array(md_atom_pair_t)* out, hid_t file, const h5md_particles_t* p, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* scratch = md_temp_allocator(temp);

    int64_t* ids = NULL;
    md_hashmap32_t id_map = { .allocator = scratch };
    if (h5md_read_ids(&ids, file, p, scratch) && ids) {
        for (size_t i = 0; i < p->num_particles; ++i) {
            if ((uint64_t)ids[i] < MD_HASH_TOMBSTONE) md_hashmap_add(&id_map, (uint64_t)ids[i], (uint32_t)i);
        }
    }

    md_array(str_t) names = 0;
    h5_children(&names, file, "connectivity", scratch);
    for (size_t n = 0; n < md_array_size(names); ++n) {
        char path[H5MD_PATH_MAX];
        if (!h5_path(path, sizeof(path), "/connectivity", names[n])) continue;
        const H5I_type_t kind = h5_kind(file, path);
        if (kind == H5I_GROUP) {
            MD_LOG_INFO("H5MD: '%s' is connectivity that changes over time, which is not read", path);
            continue;
        }
        if (kind != H5I_DATASET) continue;

        hid_t dset = H5Dopen(file, path, H5P_DEFAULT);
        if (dset < 0) continue;
        hsize_t dims[H5S_MAX_RANK];
        hid_t type = H5Dget_type(dset);
        const bool pairs = h5_shape(dims, dset) == 2 && dims[1] == 2 && type >= 0 && H5Tget_class(type) == H5T_INTEGER;
        if (type >= 0) H5Tclose(type);
        if (!pairs || !h5md_tuples_belong(dset, p)) {
            H5Dclose(dset);
            continue;
        }

        bool has_fill = false;
        int64_t fill = 0;
        hid_t dcpl = H5Dget_create_plist(dset);
        H5D_fill_value_t fill_status;
        if (dcpl >= 0 && H5Pfill_value_defined(dcpl, &fill_status) >= 0 && fill_status == H5D_FILL_VALUE_USER_DEFINED) {
            has_fill = H5Pget_fill_value(dcpl, H5T_NATIVE_INT64, &fill) >= 0;
        }
        if (dcpl >= 0) H5Pclose(dcpl);

        const size_t count = (size_t)dims[0];
        int64_t* idx = md_alloc(scratch, MAX(count, 1) * 2 * sizeof(int64_t));
        if (count && h5_read_rows(idx, dset, H5T_NATIVE_INT64, -1, count * 2)) {
            size_t skipped = 0;
            for (size_t i = 0; i < count; ++i) {
                int64_t a = idx[i * 2 + 0];
                int64_t b = idx[i * 2 + 1];
                if (has_fill && (a == fill || b == fill)) continue;
                if (ids) {
                    const uint32_t* ia = (uint64_t)a < MD_HASH_TOMBSTONE ? md_hashmap_get(&id_map, (uint64_t)a) : NULL;
                    const uint32_t* ib = (uint64_t)b < MD_HASH_TOMBSTONE ? md_hashmap_get(&id_map, (uint64_t)b) : NULL;
                    a = ia ? (int64_t)*ia : -1;
                    b = ib ? (int64_t)*ib : -1;
                }
                if (a < 0 || b < 0 || (size_t)a >= p->num_particles || (size_t)b >= p->num_particles || a == b) {
                    skipped += 1;
                    continue;
                }
                const md_atom_pair_t pair = { { (md_atom_idx_t)MIN(a, b), (md_atom_idx_t)MAX(a, b) } };
                md_array_push(*out, pair, alloc);
            }
            if (skipped) {
                MD_LOG_INFO("H5MD: %zu pairs of '%s' name no particle of '%s' and were left out", skipped, path, p->path);
            }
        }
        H5Dclose(dset);
    }
    md_temp_end(temp);
}

// Per particle values of a standard element at the first frame (or of a time independent one), as
// floats. False when the group has no such element or it is not 'components' values per particle.
static bool h5md_read_particle_values(float* out, h5md_element_t* e, hid_t file, const h5md_particles_t* p, const char* name, size_t components, md_allocator_i* alloc) {
    char path[H5MD_PATH_MAX];
    snprintf(path, sizeof(path), "%s/%s", p->path, name);
    if (!h5md_element_open(e, file, path, alloc)) {
        return false;
    }
    const int axis = e->temporal ? 1 : 0;
    const size_t count = h5_product(e->dims, axis, e->rank);
    if (e->rank <= axis || e->dims[axis] != p->num_particles || count != p->num_particles * components || (e->temporal && e->num_frames == 0)) {
        MD_LOG_INFO("H5MD: '%s' is not %zu value%s per particle and is not read", path, components, components == 1 ? "" : "s");
        return false;
    }
    return h5_read_rows_at(out, file, e->value, H5T_NATIVE_FLOAT, e->temporal ? 0 : -1, count);
}

// The GROMACS topology module: molecule types, and the blocks laying them out, as the md_tpr_data_t
// md_tpr_system_init_from_data builds a system from - which is the point, since it is the same
// topology and has to become the same system. False when the file has no module, or one that does
// not describe the particle group read.
static bool h5md_read_gromacs_topology(md_tpr_data_t* data, hid_t file, const h5md_particles_t* p, md_allocator_i* alloc) {
    hid_t top = H5Gopen(file, "h5md/modules/gromacs_topology", H5P_DEFAULT);
    if (top < 0) {
        return false;
    }

    bool result = false;
    int32_t version[2] = {0};
    const size_t num_version = h5_attr_i32s(version, 2, top, "version");
    const size_t num_blocks  = h5_attr_strs(NULL, 0, top, "molecule_block_names", alloc);
    str_t*   block_names  = md_alloc(alloc, MAX(num_blocks, 1) * sizeof(str_t));
    int32_t* block_counts = md_alloc(alloc, MAX(num_blocks, 1) * sizeof(int32_t));
    h5_attr_strs(block_names, num_blocks, top, "molecule_block_names", alloc);
    const size_t num_counts  = h5_attr_i32s(block_counts, num_blocks, top, "molecule_block_counts");

    if (num_version != 2 || version[0] != 0) {
        MD_LOG_INFO("H5MD: the GROMACS topology module is version %d.%d, which this reader does not know; the core elements are read instead", version[0], version[1]);
        goto done;
    }
    if (num_blocks == 0 || num_blocks != num_counts) {
        MD_LOG_INFO("H5MD: the GROMACS topology module has no molecule blocks it can be read by; the core elements are read instead");
        goto done;
    }

    MEMSET(data, 0, sizeof(*data));
    data->name = h5_attr_str(top, "system_name", alloc);
    data->molblocks  = md_alloc(alloc, num_blocks * sizeof(md_tpr_molblock_t));
    data->moltypes   = md_alloc(alloc, num_blocks * sizeof(md_tpr_moltype_t));   // at most one per block
    data->num_molblocks = num_blocks;
    data->nonbonded.valid = false;
    data->repulsion_power = 12.0;
    data->pbc = MD_TPR_PBC_UNSET;

    size_t num_atoms = 0;
    for (size_t b = 0; b < num_blocks; ++b) {
        // A molecule type is its name; blocks of the same type share it
        size_t t = 0;
        while (t < data->num_moltypes && !str_eq(data->moltypes[t].name, block_names[b])) ++t;
        if (t == data->num_moltypes) {
            md_tpr_moltype_t* mt = &data->moltypes[t];
            MEMSET(mt, 0, sizeof(*mt));
            mt->name = block_names[b];

            char group_path[H5MD_PATH_MAX];
            hid_t group = -1;
            if (str_copy_to_char_buf(group_path, sizeof(group_path), block_names[b]) == block_names[b].len) {
                group = H5Gopen(top, group_path, H5P_DEFAULT);
            }
            int64_t count = 0;
            if (group < 0 || !h5_attr_read(&count, group, "particle_count", H5T_NATIVE_INT64) || count <= 0) {
                MD_LOG_INFO("H5MD: the GROMACS topology module has no molecule type '" STR_FMT "'", STR_ARG(block_names[b]));
                if (group >= 0) H5Gclose(group);
                goto done;
            }
            const size_t n = (size_t)count;

            int32_t* species   = md_alloc(alloc, n * sizeof(int32_t));
            float*   mass      = md_alloc(alloc, n * sizeof(float));
            float*   charge    = md_alloc(alloc, n * sizeof(float));
            int32_t* name_idx  = md_alloc(alloc, n * sizeof(int32_t));
            int32_t* res_id    = md_alloc(alloc, n * sizeof(int32_t));
            int32_t* res_idx   = md_alloc(alloc, n * sizeof(int32_t));
            md_array(str_t) names = 0;
            md_array(str_t) res_names = 0;

            bool ok = h5_read_rows_at(species,  group, "species",       H5T_NATIVE_INT32, -1, n)
                   && h5_read_rows_at(mass,     group, "mass",          H5T_NATIVE_FLOAT, -1, n)
                   && h5_read_rows_at(charge,   group, "charge",        H5T_NATIVE_FLOAT, -1, n)
                   && h5_read_rows_at(name_idx, group, "particle_name", H5T_NATIVE_INT32, -1, n)
                   && h5_read_rows_at(res_id,   group, "residue_id",    H5T_NATIVE_INT32, -1, n)
                   && h5_read_rows_at(res_idx,  group, "residue_name",  H5T_NATIVE_INT32, -1, n);
            for (int k = 0; ok && k < 2; ++k) {
                md_array(str_t)* table = k == 0 ? &names : &res_names;
                hid_t dset = H5Dopen(group, k == 0 ? "particle_name_table" : "residue_name_table", H5P_DEFAULT);
                const size_t num = dset >= 0 ? h5_read_strs(NULL, 0, dset, false, alloc) : 0;
                if (num > 0) {
                    md_array_resize(*table, num, alloc);
                    h5_read_strs(*table, num, dset, false, alloc);
                }
                if (dset >= 0) H5Dclose(dset);
                ok = num > 0;
            }
            H5Gclose(group);
            if (!ok) {
                MD_LOG_INFO("H5MD: the molecule type '" STR_FMT "' of the GROMACS topology module is incomplete", STR_ARG(block_names[b]));
                goto done;
            }

            mt->num_atoms = n;
            mt->atoms     = md_alloc(alloc, n * sizeof(md_tpr_atom_t));
            mt->residues  = md_alloc(alloc, n * sizeof(md_tpr_residue_t));
            for (size_t i = 0; i < n; ++i) {
                if (name_idx[i] < 0 || (size_t)name_idx[i] >= md_array_size(names) || res_idx[i] < 0 || (size_t)res_idx[i] >= md_array_size(res_names)) {
                    MD_LOG_INFO("H5MD: the molecule type '" STR_FMT "' of the GROMACS topology module names a particle or residue it has no name for", STR_ARG(block_names[b]));
                    goto done;
                }
                // A residue is a run of particles with one residue id and name
                if (i == 0 || res_id[i] != res_id[i - 1] || res_idx[i] != res_idx[i - 1]) {
                    md_tpr_residue_t* res = &mt->residues[mt->num_residues++];
                    res->name = res_names[res_idx[i]];
                    // GROMACS 2026 writes one more than the residue's number in the topology (the
                    // module's version 0.1 writes resinfo.nr + 1), so this is the topology's number
                    res->nr = version[1] == 1 ? res_id[i] - 1 : res_id[i];
                    res->ic = ' ';
                }
                md_tpr_atom_t* atom = &mt->atoms[i];
                MEMSET(atom, 0, sizeof(*atom));
                atom->name          = names[name_idx[i]];
                atom->type          = STR_LIT("");
                atom->mass          = mass[i];
                atom->charge        = charge[i];
                atom->residue       = (int32_t)(mt->num_residues - 1);
                atom->atomic_number = species[i] > 0 ? species[i] : -1;
                // The module has no particle types. A particle without mass is a virtual site: it
                // is how GROMACS builds one, and it keeps its element from being guessed
                atom->ptype         = mass[i] == 0.0f ? MD_TPR_PTYPE_VSITE : MD_TPR_PTYPE_ATOM;
            }
            data->num_moltypes += 1;
        }
        data->molblocks[b].moltype = (int32_t)t;
        data->molblocks[b].nmol    = block_counts[b];
        num_atoms += (size_t)MAX(block_counts[b], 0) * data->moltypes[t].num_atoms;
    }

    if (num_atoms != p->num_particles) {
        MD_LOG_INFO("H5MD: the GROMACS topology module describes %zu atoms and '%s' holds %zu; the core elements are read instead", num_atoms, p->path, p->num_particles);
        goto done;
    }
    data->num_atoms = num_atoms;
    result = true;

done:
    H5Gclose(top);
    return result;
}

// A system from the core elements alone: species names the atom types, mass gives them their mass.
static bool h5md_system_from_core(md_system_t* sys, md_system_state_t* state, hid_t file, const h5md_particles_t* p, md_allocator_i* scratch) {
    const size_t N = p->num_particles;
    md_allocator_i* alloc = sys->alloc;
    md_system_reset(sys);
    md_system_state_init(state, N);

    h5md_element_t e;
    float* mass = md_alloc(scratch, N * sizeof(float));
    const bool has_mass = h5md_read_particle_values(mass, &e, file, p, "mass", 1, scratch);

    // Species: the name of an enumeration member, or the number of an integer species
    str_t* names = md_alloc(scratch, N * sizeof(str_t));
    bool from_label = false;
    for (size_t i = 0; i < N; ++i) names[i] = STR_LIT("X");
    char path[H5MD_PATH_MAX];
    snprintf(path, sizeof(path), "%s/species", p->path);
    if (h5md_element_open(&e, file, path, scratch)) {
        const int axis = e.temporal ? 1 : 0;
        hid_t dset = H5Dopen(file, e.value, H5P_DEFAULT);
        hid_t type = dset >= 0 ? H5Dget_type(dset) : -1;
        const H5T_class_t cls = type >= 0 ? H5Tget_class(type) : H5T_NO_CLASS;
        if (e.rank != axis + 1 || e.dims[axis] != N || (e.temporal && e.num_frames == 0) || (cls != H5T_ENUM && cls != H5T_INTEGER)) {
            MD_LOG_INFO("H5MD: '%s' is not one species per particle and is not read", path);
        } else if (cls == H5T_ENUM) {
            hid_t mem = H5Tget_native_type(type, H5T_DIR_ASCEND);
            const size_t size = H5Tget_size(mem);
            uint8_t* raw = md_alloc(scratch, N * size);
            if (h5_read_rows(raw, dset, mem, e.temporal ? 0 : -1, N)) {
                from_label = true;
                for (size_t i = 0; i < N; ++i) {
                    if (i > 0 && MEMCMP(raw + i * size, raw + (i - 1) * size, size) == 0) {
                        names[i] = names[i - 1];
                        continue;
                    }
                    char name[64] = "";
                    if (H5Tenum_nameof(mem, raw + i * size, name, sizeof(name)) < 0) name[0] = '\0';
                    names[i] = str_copy(str_from_cstr(name), scratch);
                }
            }
            H5Tclose(mem);
        } else {
            int64_t* species = md_alloc(scratch, N * sizeof(int64_t));
            if (h5_read_rows(species, dset, H5T_NATIVE_INT64, e.temporal ? 0 : -1, N)) {
                for (size_t i = 0; i < N; ++i) {
                    names[i] = (i > 0 && species[i] == species[i - 1]) ? names[i - 1] : str_printf(scratch, "%" PRId64, species[i]);
                }
            }
        }
        if (type >= 0) H5Tclose(type);
        if (dset >= 0) H5Dclose(dset);
    }

    const size_t capacity = ROUND_UP(N, 16);
    md_array_resize(sys->atom.type_idx, capacity, alloc);
    md_array_resize(sys->atom.flags, capacity, alloc);
    MEMSET(sys->atom.type_idx, 0, capacity * sizeof(md_atom_type_idx_t));
    MEMSET(sys->atom.flags, 0, capacity * sizeof(md_flags_t));
    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, alloc);

    // Particles of one species are identical (the specification says so), so a type is a species
    // and the element is decided once per species
    md_hashmap32_t type_map = { .allocator = scratch };
    for (size_t i = 0; i < N; ++i) {
        const float m = has_mass ? mass[i] : 0.0f;
        const uint64_t key = md_hash64_str(names[i], md_hash64(&m, sizeof(m), 0)) & ~(3ULL << 62);
        const uint32_t* cached = md_hashmap_get(&type_map, key);
        md_atom_type_idx_t type;
        if (cached) {
            type = (md_atom_type_idx_t)*cached;
        } else {
            md_atomic_number_t z = from_label ? md_atomic_number_infer_from_label(names[i], STR_LIT(""), 0) : 0;
            if (!z && m > 0.0f) {
                z = md_atomic_number_infer_from_mass(m);
            }
            const float type_mass = has_mass ? m : md_atomic_number_mass(z);
            type = md_atom_type_add(&sys->atom.type, names[i], STR_LIT(""), z, type_mass, md_atomic_number_vdw_radius(z), md_atomic_number_cpk_color(z), 0, alloc);
            md_hashmap_add(&type_map, key, (uint32_t)type);
        }
        sys->atom.type_idx[i] = type;
    }
    sys->atom.count = N;

    md_array(md_atom_pair_t) bonds = 0;
    h5md_read_bonds(&bonds, file, p, scratch);
    for (size_t i = 0; i < md_array_size(bonds); ++i) {
        md_array_push(sys->bond.pairs, bonds[i], alloc);
        md_array_push(sys->bond.flags, MD_BOND_FLAG_COVALENT | MD_BOND_FLAG_TOPOLOGY, alloc);
    }
    sys->bond.count = md_array_size(bonds);
    if (sys->bond.count) {
        md_bond_build_connectivity(&sys->bond, N, alloc);
    }

    float* charge = md_alloc(scratch, N * sizeof(float));
    if (h5md_read_particle_values(charge, &e, file, p, "charge", 1, scratch)) {
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), e.has_unit ? e.unit : md_unit_elementary_charge(), 1, charge, N);
    }
    return true;
}

static void h5md_publish_str(md_attributes_t* attributes, const char* path, str_t value) {
    if (str_empty(value)) return;
    md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = str_from_cstr(path), .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 } },
        .unit = md_unit_none(), .data = &value, .byte_size = sizeof(str_t)});
}

bool md_h5md_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);
    ASSERT(state);
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    char path[4096];
    str_copy_to_char_buf(path, sizeof(path), filename);

    h5md_lock_t lock = h5md_lock();
    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        MD_LOG_ERROR("H5MD: could not open '" STR_FMT "' as an HDF5 file", STR_ARG(filename));
        h5md_unlock(lock);
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* scratch = md_temp_allocator(temp);
    bool result = false;

    h5md_particles_t p;
    if (!h5md_open(&p, file, scratch)) {
        goto done;
    }
    const size_t N = p.num_particles;

    // The first frame, in Angstrom
    float* xyz  = NULL;
    float  cell[3][3] = {0};
    if (p.has_position) {
        double scale;
        xyz = md_alloc(scratch, N * 3 * sizeof(float));
        if (!h5md_length_scale(&scale, &p.position, "position") ||
            !h5_read_rows_at(xyz, file, p.position.value, H5T_NATIVE_FLOAT, p.position.temporal ? 0 : -1, N * 3)) {
            MD_LOG_ERROR("H5MD: failed to read the positions of '%s'", p.path);
            goto done;
        }
        for (size_t i = 0; i < N * 3; ++i) xyz[i] = (float)(xyz[i] * scale);

        float* cells = md_alloc(scratch, (p.position.temporal ? p.position.num_frames : 1) * 9 * sizeof(float));
        if (!h5md_read_cells(cells, file, &p, scratch)) {
            goto done;
        }
        MEMCPY(cell, cells, sizeof(cell));
    } else {
        MD_LOG_INFO("H5MD: '%s' has no positions", p.path);
    }

    md_tpr_data_t top;
    if (h5md_read_gromacs_topology(&top, file, &p, scratch)) {
        md_array(md_atom_pair_t) bonds = 0;
        h5md_read_bonds(&bonds, file, &p, scratch);
        top.num_intermolecular_bonds = md_array_size(bonds);
        top.intermolecular_bonds     = bonds;
        if (xyz) {
            // md_tpr_system_init_from_data takes nm; the coordinates are set exactly below
            top.x = md_alloc(scratch, N * 3 * sizeof(float));
            for (size_t i = 0; i < N * 3; ++i) top.x[i] = xyz[i] * 0.1f;
        }
        if (!md_tpr_system_init_from_data(sys, state, &top)) {
            goto done;
        }
    } else if (!h5md_system_from_core(sys, state, file, &p, scratch)) {
        goto done;
    }

    if (xyz) {
        for (size_t i = 0; i < N; ++i) {
            state->xyz[i] = vec3_set(xyz[i * 3 + 0], xyz[i * 3 + 1], xyz[i * 3 + 2]);
        }
        state->unitcell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(cell));
    }
    // A file without positions is a topology: no coordinates, which is what zero atoms says
    state->num_atoms = xyz ? N : 0;

    // Per particle vectors written without frames belong to the structure
    static const char* columns[] = { "velocity", "force", "image" };
    for (size_t k = 0; k < ARRAY_SIZE(columns); ++k) {
        char probe[H5MD_PATH_MAX];
        snprintf(probe, sizeof(probe), "%s/%s", p.path, columns[k]);
        h5md_element_t e;
        if (h5_kind(file, probe) != H5I_DATASET) continue;
        float* values = md_alloc(scratch, N * 3 * sizeof(float));
        if (h5md_read_particle_values(values, &e, file, &p, columns[k], 3, scratch)) {
            char attr_path[64];
            snprintf(attr_path, sizeof(attr_path), "atom/%s", columns[k]);
            md_attributes_publish_atom_column(&sys->attributes, str_from_cstr(attr_path), e.unit, 3, values, N);
        }
    }

    // Where the file came from
    {
        md_attributes_t* attributes = &sys->attributes;
        int32_t version[2] = {0};
        hid_t meta = H5Gopen(file, "h5md", H5P_DEFAULT);
        h5_attr_i32s(version, 2, meta, "version");
        md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = STR_LIT("h5md/version"), .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { 2 } },
            .unit = md_unit_none(), .data = version, .byte_size = sizeof(version)});
        hid_t author  = H5Gopen(meta, "author", H5P_DEFAULT);
        hid_t creator = H5Gopen(meta, "creator", H5P_DEFAULT);
        if (author >= 0) {
            h5md_publish_str(attributes, "h5md/author/name",  h5_attr_str(author, "name", scratch));
            h5md_publish_str(attributes, "h5md/author/email", h5_attr_str(author, "email", scratch));
            H5Gclose(author);
        }
        if (creator >= 0) {
            h5md_publish_str(attributes, "h5md/creator/name",    h5_attr_str(creator, "name", scratch));
            h5md_publish_str(attributes, "h5md/creator/version", h5_attr_str(creator, "version", scratch));
            H5Gclose(creator);
        }
        H5Gclose(meta);
    }
    result = true;

done:
    md_temp_end(temp);
    H5Fclose(file);
    h5md_unlock(lock);
    return result;
}

// ### RUN ###

// How the bytes of a frame read when they can be read without HDF5
typedef enum h5md_raw_t {
    H5MD_RAW_NONE = 0,  // not one contiguous run of IEEE floats: read through HDF5
    H5MD_RAW_F32LE,
    H5MD_RAW_F32BE,
    H5MD_RAW_F64LE,
    H5MD_RAW_F64BE,
} h5md_raw_t;

// A provider's own state, one per attribute and released with it
typedef struct h5md_source_t {
    const md_system_t* sys;         // borrowed: the run goes before the system does
    md_attribute_id_t  file_id;     // <run>/source/path
    md_attribute_id_t  offset_id;   // where each frame is in the file
    uint32_t           raw;         // h5md_raw_t
    uint32_t           num_atoms;
    uint32_t           components;
    float              scale;       // the file's unit to the published one
    char               dataset[H5MD_PATH_MAX];
} h5md_source_t;

#if H5_VERSION_GE(1, 14, 0)
typedef struct h5md_chunk_walk_t {
    int64_t* out;
    size_t   num_frames;
    size_t   frames_per_chunk;
    size_t   frame_bytes;
} h5md_chunk_walk_t;

static int h5md_chunk_visit(const hsize_t* offset, unsigned filter_mask, haddr_t addr, hsize_t size, void* data) {
    (void)filter_mask;
    (void)size;
    const h5md_chunk_walk_t* walk = (const h5md_chunk_walk_t*)data;
    const size_t first = (size_t)offset[0];
    for (size_t f = first; f < MIN(walk->num_frames, first + walk->frames_per_chunk); ++f) {
        walk->out[f] = (int64_t)(addr + (f - first) * walk->frame_bytes);
    }
    return H5_ITER_CONT;
}
#endif

// Where each frame of a [F, ...] dataset starts in the file, when every frame is one contiguous
// run of bytes of IEEE floats: contiguous storage, or chunks spanning whole frames with no filter.
// H5MD_RAW_NONE otherwise, and then out is -1 throughout.
static h5md_raw_t h5md_frame_offsets(int64_t* out, hid_t dset, size_t num_frames) {
    for (size_t f = 0; f < num_frames; ++f) out[f] = -1;

    hsize_t dims[H5S_MAX_RANK];
    const int rank = h5_shape(dims, dset);
    hid_t type = H5Dget_type(dset);
    h5md_raw_t raw = H5MD_RAW_NONE;
    if (type >= 0) {
        if      (H5Tequal(type, H5T_IEEE_F32LE) > 0) raw = H5MD_RAW_F32LE;
        else if (H5Tequal(type, H5T_IEEE_F32BE) > 0) raw = H5MD_RAW_F32BE;
        else if (H5Tequal(type, H5T_IEEE_F64LE) > 0) raw = H5MD_RAW_F64LE;
        else if (H5Tequal(type, H5T_IEEE_F64BE) > 0) raw = H5MD_RAW_F64BE;
        H5Tclose(type);
    }
    if (raw == H5MD_RAW_NONE || rank < 1 || dims[0] != num_frames) {
        return H5MD_RAW_NONE;
    }
    const size_t elem = (raw == H5MD_RAW_F32LE || raw == H5MD_RAW_F32BE) ? 4 : 8;
    const size_t frame_bytes = h5_product(dims, 1, rank) * elem;

    hid_t dcpl = H5Dget_create_plist(dset);
    if (dcpl < 0) {
        return H5MD_RAW_NONE;
    }
    bool ok = false;
    const H5D_layout_t layout = H5Pget_layout(dcpl);
    if (layout == H5D_CONTIGUOUS) {
        const haddr_t addr = H5Dget_offset(dset);
        if (addr != HADDR_UNDEF) {
            for (size_t f = 0; f < num_frames; ++f) out[f] = (int64_t)(addr + f * frame_bytes);
            ok = true;
        }
    }
#if H5_VERSION_GE(1, 10, 5)
    else if (layout == H5D_CHUNKED && H5Pget_nfilters(dcpl) == 0) {
        // A chunk holds chunk[0] whole frames, one after the other, when it spans the other axes
        hsize_t chunk[H5S_MAX_RANK];
        ok = H5Pget_chunk(dcpl, rank, chunk) == rank && chunk[0] > 0;
        for (int i = 1; ok && i < rank; ++i) ok = chunk[i] == dims[i];
#if H5_VERSION_GE(1, 14, 0)
        // One walk over the chunk index. Asking for each chunk by its coordinate costs a lookup of
        // a third of a millisecond, which for a long trajectory is most of the time publishing takes
        h5md_chunk_walk_t walk = { out, num_frames, (size_t)chunk[0], frame_bytes };
        ok = ok && H5Dchunk_iter(dset, H5P_DEFAULT, h5md_chunk_visit, &walk) >= 0;
#else
        for (size_t f = 0; ok && f < num_frames; f += (size_t)chunk[0]) {
            hsize_t coord[H5S_MAX_RANK] = {0};
            coord[0] = (hsize_t)f;
            unsigned mask = 0;
            haddr_t  addr = HADDR_UNDEF;
            hsize_t  size = 0;
            ok = H5Dget_chunk_info_by_coord(dset, coord, &mask, &addr, &size) >= 0 && addr != HADDR_UNDEF;
            for (size_t i = f; ok && i < MIN(num_frames, f + (size_t)chunk[0]); ++i) {
                out[i] = (int64_t)(addr + (i - f) * frame_bytes);
            }
        }
#endif
        // A chunk never written has no place in the file
        for (size_t f = 0; ok && f < num_frames; ++f) ok = out[f] >= 0;
    }
#endif
    H5Pclose(dcpl);
    if (!ok) {
        for (size_t f = 0; f < num_frames; ++f) out[f] = -1;
        return H5MD_RAW_NONE;
    }
    return raw;
}

static uint32_t h5md_load_u32_be(const uint8_t* p) {
    return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | (uint32_t)p[3];
}

static uint32_t h5md_load_u32_le(const uint8_t* p) {
    return ((uint32_t)p[3] << 24) | ((uint32_t)p[2] << 16) | ((uint32_t)p[1] << 8) | (uint32_t)p[0];
}

// Values of the frame as floats, from the bytes of the file
static void h5md_decode(float* out, const uint8_t* raw, size_t count, h5md_raw_t kind, float scale) {
    for (size_t i = 0; i < count; ++i) {
        switch (kind) {
        case H5MD_RAW_F32LE:
        case H5MD_RAW_F32BE: {
            const uint32_t bits = kind == H5MD_RAW_F32LE ? h5md_load_u32_le(raw + i * 4) : h5md_load_u32_be(raw + i * 4);
            float v;
            MEMCPY(&v, &bits, 4);
            out[i] = v * scale;
            break;
        }
        case H5MD_RAW_F64LE:
        case H5MD_RAW_F64BE: {
            const uint8_t* p = raw + i * 8;
            const uint64_t lo = kind == H5MD_RAW_F64LE ? h5md_load_u32_le(p)     : h5md_load_u32_be(p + 4);
            const uint64_t hi = kind == H5MD_RAW_F64LE ? h5md_load_u32_le(p + 4) : h5md_load_u32_be(p);
            const uint64_t bits = (hi << 32) | lo;
            double v;
            MEMCPY(&v, &bits, 8);
            out[i] = (float)(v * scale);
            break;
        }
        default:
            out[i] = 0.0f;
            break;
        }
    }
}

// Values of a frame through HDF5, for storage that is not plain bytes: particles first to first +
// count of frame 'frame', as floats.
static bool h5md_read_hdf5(float* out, str_t file_path, const char* dataset, uint32_t frame, size_t first, size_t count, size_t components) {
    char path[4096];
    str_copy_to_char_buf(path, sizeof(path), file_path);

    h5md_lock_t lock = h5md_lock();
    bool ok = false;
    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dset = file >= 0 ? H5Dopen(file, dataset, H5P_DEFAULT) : -1;
    if (dset >= 0) {
        hsize_t dims[H5S_MAX_RANK];
        const int rank = h5_shape(dims, dset);
        if (rank == 2 || rank == 3) {
            const hsize_t start[3]  = { frame, first, 0 };
            const hsize_t extent[3] = { 1, count, rank == 3 ? dims[2] : 1 };
            const hsize_t mem_dims  = (hsize_t)(count * components);
            hid_t file_space = H5Dget_space(dset);
            hid_t mem_space  = H5Screate_simple(1, &mem_dims, NULL);
            ok = h5_product(extent, 0, rank) == count * components &&
                 H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, extent, NULL) >= 0 &&
                 H5Dread(dset, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT, out) >= 0;
            H5Sclose(mem_space);
            H5Sclose(file_space);
        }
        H5Dclose(dset);
    }
    if (file >= 0) H5Fclose(file);
    h5md_unlock(lock);
    return ok;
}

// Every per particle quantity of a run: <run>/atom/position and the elements beside it
static size_t h5md_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    (void)attr;
    const h5md_source_t* src = (const h5md_source_t*)user_data;
    // A whole temporal virtual attribute is refused before it gets here, so the frame is fixed
    if (!src || !slice || slice->num_idx == 0 || slice->num_idx > 2) {
        return 0;
    }
    const md_attributes_t* attributes = &src->sys->attributes;
    const md_attribute_t* file   = md_attributes_get(attributes, src->file_id);
    const md_attribute_t* offset = md_attributes_get(attributes, src->offset_id);
    const int64_t* offsets = (const int64_t*)md_attribute_view(offset, MD_ATTRIBUTE_TYPE_I64, 1, 1);
    if (!file || !offsets) {
        MD_LOG_ERROR("H5MD: the run has lost its source attributes");
        return 0;
    }

    const uint32_t frame = slice->idx[0];
    const size_t K = src->components;
    size_t first = 0;
    size_t count = src->num_atoms;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= src->num_atoms) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (frame >= offset->format.shape[0] || cap != count * K) {
        return 0;
    }

    const str_t path = md_attribute_str(attributes, file, 0);
    const int64_t at = offsets[frame];
    float* out = (float*)dst;

    if (src->raw != H5MD_RAW_NONE && at >= 0) {
        const size_t elem  = (src->raw == H5MD_RAW_F32LE || src->raw == H5MD_RAW_F32BE) ? 4 : 8;
        const size_t bytes = count * K * elem;
        md_temp_scope_t temp = md_temp_begin();
        uint8_t* raw = md_temp_alloc(temp, bytes);
        const bool ok = raw && md_attribute_io_read_at(io, path, at + (int64_t)(first * K * elem), raw, bytes) == bytes;
        if (ok) {
            h5md_decode(out, raw, count * K, (h5md_raw_t)src->raw, src->scale);
        }
        md_temp_end(temp);
        if (!ok) {
            MD_LOG_ERROR("H5MD: failed to read frame %u of '%s' from '" STR_FMT "'", frame, src->dataset, STR_ARG(path));
            return 0;
        }
        return cap;
    }

    if (!h5md_read_hdf5(out, path, src->dataset, frame, first, count, K)) {
        MD_LOG_ERROR("H5MD: failed to read frame %u of '%s' from '" STR_FMT "'", frame, src->dataset, STR_ARG(path));
        return 0;
    }
    if (src->scale != 1.0f) {
        for (size_t i = 0; i < count * K; ++i) out[i] *= src->scale;
    }
    return cap;
}

// Publishes a virtual attribute read by h5md_provider, with src as its own state. Freed again when
// the attribute cannot be published, since only a published attribute releases its user_data.
static bool h5md_publish_virtual(md_attributes_t* attributes, md_attribute_desc_t desc, const h5md_source_t* src) {
    h5md_source_t* state = md_attributes_alloc_user_data(attributes, sizeof(h5md_source_t));
    if (!state) {
        return false;
    }
    *state = *src;
    const md_attribute_virtual_t virt = { .provider = h5md_provider, .user_data = state, .user_data_size = sizeof(h5md_source_t) };
    desc.virt = &virt;
    if (md_attributes_replace(attributes, &desc) == MD_ATTRIBUTE_INVALID) {
        md_free(attributes->alloc, state, sizeof(h5md_source_t));
        return false;
    }
    return true;
}

// The time of each step, for an element without a time of its own: along the line the run's steps
// and times lie on, which a constant time step puts them on. False when they do not lie on one.
static bool h5md_time_from_step(double* out, const int64_t* step, size_t n, const int64_t* run_step, const double* run_time, size_t num_frames) {
    if (num_frames < 2 || run_step[num_frames - 1] == run_step[0]) {
        return false;
    }
    const double dt = (run_time[num_frames - 1] - run_time[0]) / (double)(run_step[num_frames - 1] - run_step[0]);
    for (size_t f = 0; f < num_frames; ++f) {
        const double t = run_time[0] + dt * (double)(run_step[f] - run_step[0]);
        if (fabs(t - run_time[f]) > 1.0e-6 * MAX(fabs(run_time[f]), 1.0)) {
            return false;
        }
    }
    for (size_t i = 0; i < n; ++i) {
        out[i] = run_time[0] + dt * (double)(step[i] - run_step[0]);
    }
    return true;
}

// The frame axis of an element with frames of its own: '<group>/time' and '<group>/step'. Its time
// when it has one, derived from its step through the run's otherwise; and the step itself, without
// a unit, when neither can be had - which only ever lines up with a run that has no time either.
static bool h5md_publish_axis(md_attributes_t* attributes, str_t group, const h5md_element_t* e, const int64_t* run_step, const double* run_time, md_unit_t run_time_unit, size_t num_frames, md_allocator_i* scratch) {
    const size_t n = e->num_frames;
    double* time = e->time;
    md_unit_t unit = e->time_unit;
    if (!time) {
        time = md_alloc(scratch, n * sizeof(double));
        unit = run_time_unit;
        if (md_unit_is_none(run_time_unit) || !h5md_time_from_step(time, e->step, n, run_step, run_time, num_frames)) {
            for (size_t i = 0; i < n; ++i) time[i] = (double)e->step[i];
            unit = md_unit_none();
        }
    }
    char buf[512];
    const md_attribute_format_t series_f64 = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
    const md_attribute_format_t series_i64 = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
    return md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, STR_LIT("time")), .format = series_f64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = unit, .label = STR_INIT("Time"), .data = time, .byte_size = n * sizeof(double)}) != MD_ATTRIBUTE_INVALID
        && md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, STR_LIT("step")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_none(), .label = STR_INIT("Step"), .data = e->step, .byte_size = n * sizeof(int64_t)}) != MD_ATTRIBUTE_INVALID;
}

// Every time dependent per particle element but the position (velocity, force, image, ...): beside
// the positions when sampled at their steps, in '<run>/h5md/particles/<name>' with its own frame
// axis otherwise. Elements that are not a value per particle, or not numbers, are passed over.
static bool h5md_publish_particle_elements(md_system_t* sys, hid_t file, const h5md_particles_t* p, str_t run, const h5md_source_t* position, const double* run_time, md_unit_t run_time_unit, md_allocator_i* scratch) {
    md_attributes_t* attributes = &sys->attributes;
    const size_t F = p->position.num_frames;
    const size_t N = p->num_particles;

    md_array(str_t) names = 0;
    h5_children(&names, file, p->path, scratch);
    for (size_t n = 0; n < md_array_size(names); ++n) {
        const str_t name = names[n];
        if (str_eq_cstr(name, "position") || str_eq_cstr(name, "box") || str_eq_cstr(name, "id")) continue;

        char path[H5MD_PATH_MAX];
        h5md_element_t e;
        if (!h5_path(path, sizeof(path), p->path, name) || !h5md_element_open(&e, file, path, scratch) || !e.temporal || e.num_frames == 0) continue;
        if (!((e.rank == 2 || e.rank == 3) && e.dims[1] == N)) continue;
        hid_t dset = H5Dopen(file, e.value, H5P_DEFAULT);
        if (dset < 0) continue;
        hid_t type = H5Dget_type(dset);
        const H5T_class_t cls = type >= 0 ? H5Tget_class(type) : H5T_NO_CLASS;
        if (type >= 0) H5Tclose(type);
        if (cls != H5T_FLOAT && cls != H5T_INTEGER) {
            H5Dclose(dset);
            continue;
        }

        const size_t K = e.rank == 3 ? (size_t)e.dims[2] : 1;
        int64_t* offsets = md_alloc(scratch, e.num_frames * sizeof(int64_t));
        const h5md_raw_t raw = h5md_frame_offsets(offsets, dset, e.num_frames);
        H5Dclose(dset);

        const bool aligned = e.num_frames == F && MEMCMP(e.step, p->position.step, F * sizeof(int64_t)) == 0;
        char group_buf[512];
        const int group_len = snprintf(group_buf, sizeof(group_buf), STR_FMT "/h5md/particles/" STR_FMT, STR_ARG(run), STR_ARG(name));
        char leaf_buf[256];
        const int leaf_len = snprintf(leaf_buf, sizeof(leaf_buf), "atom/" STR_FMT, STR_ARG(name));
        if (group_len <= 0 || group_len >= (int)sizeof(group_buf) || leaf_len <= 0 || leaf_len >= (int)sizeof(leaf_buf)) continue;
        const str_t group = { group_buf, (size_t)group_len };
        const str_t leaf  = { leaf_buf, (size_t)leaf_len };

        char buf[512];
        bool ok = aligned || h5md_publish_axis(attributes, group, &e, p->position.step, run_time, run_time_unit, F, scratch);
        const str_t offset_path = md_run_path(buf, sizeof(buf), group, STR_LIT("source/offset"));
        ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = offset_path, .format = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)e.num_frames } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .data = offsets, .byte_size = e.num_frames * sizeof(int64_t)}) != MD_ATTRIBUTE_INVALID;

        h5md_source_t src = *position;
        src.offset_id  = md_attributes_id_from_path(offset_path);
        src.raw        = raw;
        src.components = (uint32_t)K;
        src.scale      = 1.0f;
        snprintf(src.dataset, sizeof(src.dataset), "%s", e.value);
        ok = ok && h5md_publish_virtual(attributes, (md_attribute_desc_t){
            .path = aligned ? md_run_path(buf, sizeof(buf), run, leaf) : md_run_path(buf, sizeof(buf), group, leaf),
            .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = (uint32_t)K, .rank = 2, .shape = { (uint32_t)e.num_frames, (uint32_t)N } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = e.unit}, &src);
        if (!ok) {
            MD_LOG_ERROR("H5MD: failed to publish '%s'", path);
            return false;
        }
    }
    return true;
}

// The observables below 'dir' in the file, each to '<run>/h5md/observables/<path>/value' as the file
// has it, with its own frame axis beside it when it is time dependent. Anything that is not numbers
// is passed over.
static bool h5md_publish_observables(md_system_t* sys, hid_t file, const char* dir, str_t run, const int64_t* run_step, const double* run_time, md_unit_t run_time_unit, size_t num_frames, md_allocator_i* scratch) {
    md_attributes_t* attributes = &sys->attributes;
    md_array(str_t) names = 0;
    h5_children(&names, file, dir, scratch);
    for (size_t n = 0; n < md_array_size(names); ++n) {
        char path[H5MD_PATH_MAX];
        if (!h5_path(path, sizeof(path), dir, names[n])) continue;

        h5md_element_t e;
        if (!h5md_element_open(&e, file, path, scratch)) {
            // Not an element: a group of them
            if (h5_kind(file, path) == H5I_GROUP && !h5md_publish_observables(sys, file, path, run, run_step, run_time, run_time_unit, num_frames, scratch)) {
                return false;
            }
            continue;
        }
        if (e.rank > MD_ATTRIBUTE_MAX_RANK || (e.temporal && e.num_frames == 0)) continue;
        bool numeric = false;
        hid_t dset = H5Dopen(file, e.value, H5P_DEFAULT);
        if (dset >= 0) {
            hid_t type = H5Dget_type(dset);
            numeric = type >= 0 && (H5Tget_class(type) == H5T_FLOAT || H5Tget_class(type) == H5T_INTEGER);
            if (type >= 0) H5Tclose(type);
            H5Dclose(dset);
        }
        if (!numeric) continue;

        md_attribute_format_t format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = (uint32_t)e.rank };
        size_t count = 1;
        for (int i = 0; i < e.rank; ++i) {
            format.shape[i] = (uint32_t)e.dims[i];
            count *= (size_t)e.dims[i];
        }
        if (count == 0) continue;
        double* values = md_alloc(scratch, count * sizeof(double));
        if (!h5_read_rows_at(values, file, e.value, H5T_NATIVE_DOUBLE, -1, count)) {
            MD_LOG_INFO("H5MD: failed to read '%s'", path);
            continue;
        }

        // "/observables/a/b" -> "<run>/h5md/observables/a/b"
        char group_buf[512];
        const int group_len = snprintf(group_buf, sizeof(group_buf), STR_FMT "/h5md%s", STR_ARG(run), path);
        if (group_len <= 0 || group_len >= (int)sizeof(group_buf)) continue;
        const str_t group = { group_buf, (size_t)group_len };

        char buf[512];
        bool ok = !e.temporal || h5md_publish_axis(attributes, group, &e, run_step, run_time, run_time_unit, num_frames, scratch);
        ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, STR_LIT("value")), .format = format,
            .flags = e.temporal ? MD_ATTRIBUTE_FLAG_TEMPORAL : MD_ATTRIBUTE_FLAG_NONE, .unit = e.unit,
            .data = values, .byte_size = count * sizeof(double)}) != MD_ATTRIBUTE_INVALID;
        if (!ok) {
            MD_LOG_ERROR("H5MD: failed to publish '%s'", path);
            return false;
        }
    }
    return true;
}

bool md_h5md_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    (void)flags;
    if (str_empty(run)) {
        MD_LOG_ERROR("H5MD: no run to publish into");
        return false;
    }
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    if (path_len == 0) {
        MD_LOG_ERROR("H5MD: could not resolve '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    const str_t path = { path_buf, path_len };
    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) {
        attributes->alloc = sys->alloc;
    }

    h5md_lock_t lock = h5md_lock();
    hid_t file = H5Fopen(path_buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        MD_LOG_ERROR("H5MD: could not open '" STR_FMT "' as an HDF5 file", STR_ARG(path));
        h5md_unlock(lock);
        return false;
    }

    md_allocator_i* scratch = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;
    bool published = false;

    h5md_particles_t p;
    if (!h5md_open(&p, file, scratch)) {
        goto done;
    }
    if (!p.has_position || !p.position.temporal) {
        MD_LOG_ERROR("H5MD: '%s' has no position frames to make a run of", p.path);
        goto done;
    }
    const size_t F = p.position.num_frames;
    const size_t N = p.num_particles;
    if (sys->atom.count != 0 && sys->atom.count != N) {
        MD_LOG_ERROR("H5MD: '%s' holds %zu particles and the system %zu atoms", p.path, N, sys->atom.count);
        goto done;
    }

    // The run's time: the position's, or its step without a unit
    double* time = p.position.time;
    md_unit_t time_unit = p.position.time_unit;
    if (!time) {
        time = md_alloc(scratch, F * sizeof(double));
        for (size_t f = 0; f < F; ++f) time[f] = (double)p.position.step[f];
        time_unit = md_unit_none();
    }

    float* cells = md_alloc(scratch, F * 9 * sizeof(float));
    double scale;
    if (!h5md_read_cells(cells, file, &p, scratch) || !h5md_length_scale(&scale, &p.position, "position")) {
        goto done;
    }

    int64_t* offsets = md_alloc(scratch, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(scratch, F * sizeof(int64_t));
    hid_t dset = H5Dopen(file, p.position.value, H5P_DEFAULT);
    const h5md_raw_t raw = dset >= 0 ? h5md_frame_offsets(offsets, dset, F) : H5MD_RAW_NONE;
    if (dset >= 0) H5Dclose(dset);
    const size_t elem = (raw == H5MD_RAW_F64LE || raw == H5MD_RAW_F64BE) ? 8 : 4;
    for (size_t f = 0; f < F; ++f) sizes[f] = raw != H5MD_RAW_NONE ? (int64_t)(N * 3 * elem) : 0;

    char buf[512];
    h5md_source_t src = {
        .sys        = sys,
        .file_id    = md_attributes_id_from_path(md_run_path(buf, sizeof(buf), run, STR_LIT("source/path"))),
        .offset_id  = md_attributes_id_from_path(md_run_path(buf, sizeof(buf), run, STR_LIT("source/offset"))),
        .raw        = raw,
        .num_atoms  = (uint32_t)N,
        .components = 3,
        .scale      = (float)scale,
    };
    snprintf(src.dataset, sizeof(src.dataset), "%s", p.position.value);

    h5md_source_t* state = md_attributes_alloc_user_data(attributes, sizeof(h5md_source_t));
    if (!state) {
        goto done;
    }
    *state = src;
    const md_attribute_virtual_t virt = { .provider = h5md_provider, .user_data = state, .user_data_size = sizeof(h5md_source_t) };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = time,
        .time_unit     = time_unit,
        .step          = p.position.step,
        .unitcell      = (p.has_edges && (p.periodic[0] || p.periodic[1] || p.periodic[2])) ? cells : NULL,
        .source_path   = path,
        .source_offset = offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    if (!md_run_publish(sys, run, &desc)) {
        // The position is the last thing published, so its state was never taken
        md_free(attributes->alloc, state, sizeof(h5md_source_t));
        goto done;
    }
    published = true;

    if (!h5md_publish_particle_elements(sys, file, &p, run, &src, time, time_unit, scratch) ||
        !h5md_publish_observables(sys, file, "/observables", run, p.position.step, time, time_unit, F, scratch)) {
        goto done;
    }
    result = true;

done:
    if (published && !result) {
        md_attributes_remove_prefix(attributes, run);
    }
    md_arena_allocator_destroy(scratch);
    H5Fclose(file);
    h5md_unlock(lock);
    return result;
}

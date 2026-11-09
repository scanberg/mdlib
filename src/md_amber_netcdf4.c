#include <md_amber_netcdf4.h>

#include <md_hdf5.h>
#include <md_hdf5_direct.h>
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
// HDF5 is not concurrent, and every call to it is made under the one mdlib HDF5 lock (md_hdf5.h). A frame is therefore
// not read through HDF5 but through md_hdf5_direct.h: how each per atom variable is stored is learned once, when the
// run is published, and a frame is then read and decoded on whichever thread asks, without the lock.

#define NC4_PATH_MAX 4096
#define NC4_NAME_MAX 32

// A provider's own state, freed with its attribute: this, then the storage of the variable (md_hdf5_direct_t)
typedef struct nc4_source_t {
    char   dataset[NC4_NAME_MAX];
    char   path[NC4_PATH_MAX];
    double scale;               // what a stored value is multiplied by: scale_factor / icompressfac
} nc4_source_t;

MD_STATIC_ASSERT(sizeof(nc4_source_t) % 8 == 0, "the storage after it must be aligned");

static const md_hdf5_direct_t* nc4_source_direct(const nc4_source_t* src) {
    return (const md_hdf5_direct_t*)(src + 1);
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
                len = MD_MIN(strlen(s), cap - 1);
                MD_MEMCPY(buf, s, len);
                H5free_memory(s);
            }
            H5Tclose(mt);
        } else {
            const hid_t space = H5Aget_space(a);
            const hssize_t n = H5Sget_simple_extent_npoints(space);
            H5Sclose(space);
            const size_t size = H5Tget_size(t) * (size_t)MD_MAX(n, 1);
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
    const md_str_t c = md_str_from_cstr(conventions);
    size_t i = 0;
    while (i < c.len) {
        while (i < c.len && (c.ptr[i] == ',' || c.ptr[i] == ' ')) ++i;
        size_t j = i;
        while (j < c.len && c.ptr[j] != ',' && c.ptr[j] != ' ') ++j;
        if (j > i && md_str_eq_cstr((md_str_t){c.ptr + i, j - i}, name)) return true;
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

// ## A frame

static size_t nc4_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const nc4_source_t* src = (const nc4_source_t*)user_data;
    if (!src || !slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;
    const size_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    if (frame >= attr->format.shape[0]) return 0;
    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) return 0;
    const md_str_t path = md_str_from_cstr(src->path);
    if (!md_hdf5_direct_read_f32((float*)dst, nc4_source_direct(src), path, src->dataset, frame, first * 3, count * 3, src->scale, io)) {
        MD_LOG_ERROR("Amber NetCDF-4: failed to read frame %zu of '%s' from '" MD_STR_FMT "'", frame, src->dataset, MD_STR_ARG(path));
        return 0;
    }
    return cap;
}

// ## Entry points

static bool nc4_is_angstrom(hid_t file, const char* name) {
    const hid_t d = H5Dopen2(file, name, H5P_DEFAULT);
    if (d < 0) return false;
    char units[64];
    const size_t len = h5_attr_text(d, "units", units, sizeof(units));
    H5Dclose(d);
    return len == 0 || md_str_eq_cstr_ignore_case(md_str_from_cstr(units), "angstrom") || md_str_eq_cstr_ignore_case(md_str_from_cstr(units), "angstroms");
}

md_amber_netcdf_kind_t md_amber_nc4_file_kind(md_str_t filename, size_t* out_num_atoms) {
    char path[NC4_PATH_MAX];
    md_str_copy_to_char_buf(path, sizeof(path), filename);
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

bool md_amber_nc4_state_init_from_file(md_system_state_t* state, md_str_t filename, size_t num_atoms) {
    MD_ASSERT(state && state->alloc);
    char path[NC4_PATH_MAX];
    md_str_copy_to_char_buf(path, sizeof(path), filename);

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
            MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' does not follow the Amber convention", MD_STR_ARG(filename));
        } else if (!N || dims[rank - 1] != 3 || (rank == 3 && dims[0] == 0)) {
            MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' has no coordinates", MD_STR_ARG(filename));
        } else if (num_atoms && N != num_atoms) {
            MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' holds %zu atoms, the system %zu", MD_STR_ARG(filename), N, num_atoms);
        } else if (!nc4_is_angstrom(file, name)) {
            MD_LOG_ERROR("Amber NetCDF-4: the coordinates of '" MD_STR_FMT "' are not in angstrom", MD_STR_ARG(filename));
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
            if (!ok) MD_LOG_ERROR("Amber NetCDF-4: failed to read the coordinates of '" MD_STR_FMT "'", MD_STR_ARG(filename));
        }
        H5Fclose(file);
    } else {
        MD_LOG_ERROR("Amber NetCDF-4: cannot open '" MD_STR_FMT "'", MD_STR_ARG(filename));
    }
    md_hdf5_unlock(lock);

    if (ok) ok = md_system_state_init(state, N);
    if (ok) {
        MD_MEMCPY(state->xyz, xyz, N * 3 * sizeof(float));
        state->num_atoms = N;
        state->unitcell = cell;
    }
    md_temp_end(temp);
    return ok;
}

static md_unit_t nc4_time_unit(const char* units, bool* known) {
    const md_str_t u = md_str_from_cstr(units);
    *known = true;
    if (md_str_empty(u) || md_str_eq_cstr_ignore_case(u, "picosecond") || md_str_eq_cstr_ignore_case(u, "picoseconds") || md_str_eq_cstr_ignore_case(u, "ps")) return md_unit_picosecond();
    if (md_str_eq_cstr_ignore_case(u, "nanosecond") || md_str_eq_cstr_ignore_case(u, "nanoseconds") || md_str_eq_cstr_ignore_case(u, "ns")) return md_unit_nanosecond();
    if (md_str_eq_cstr_ignore_case(u, "femtosecond") || md_str_eq_cstr_ignore_case(u, "femtoseconds") || md_str_eq_cstr_ignore_case(u, "fs")) return md_unit_femtosecond();
    *known = false;
    return md_unit_none();
}

bool md_amber_nc4_system_publish_run(md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags) {
    MD_ASSERT(sys);
    (void)flags;
    char path[NC4_PATH_MAX];
    const size_t path_len = md_path_write_canonical(path, sizeof(path), filename);
    if (!path_len) {
        MD_LOG_ERROR("Amber NetCDF-4: cannot resolve '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }
    const md_str_t spath = {path, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    bool result = false;

    // ## Everything the run needs from HDF5, at once and under the lock
    md_hdf5_direct_t* direct[3] = {0};
    double scale[3] = { 1.0, 1.0, 1.0 };
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
            MD_LOG_ERROR("Amber NetCDF-4: cannot open '" MD_STR_FMT "'", MD_STR_ARG(spath));
        } else {
            const md_amber_netcdf_kind_t kind = nc4_kind(file);
            hsize_t dims[3] = {0};
            names[0] = nc4_section_name(file, 0);
            const int rank = names[0] ? h5_dims(file, names[0], dims) : 0;
            if (kind == MD_AMBER_NETCDF_RESTART) {
                MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' is a restart, a single set of coordinates rather than a trajectory", MD_STR_ARG(spath));
            } else if (kind != MD_AMBER_NETCDF_TRAJECTORY) {
                MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' does not follow the Amber convention (no Conventions \"AMBER\")", MD_STR_ARG(spath));
            } else if (rank != 3 || dims[0] == 0 || dims[1] == 0 || dims[2] != 3 || dims[0] > UINT32_MAX || dims[1] > UINT32_MAX) {
                MD_LOG_ERROR("Amber NetCDF-4: '" MD_STR_FMT "' has no coordinates of {frame, atom, spatial}", MD_STR_ARG(spath));
            } else if (!nc4_is_angstrom(file, names[0])) {
                MD_LOG_ERROR("Amber NetCDF-4: the coordinates of '" MD_STR_FMT "' are not in angstrom", MD_STR_ARG(spath));
            } else {
                F = (size_t)dims[0];
                N = (size_t)dims[1];
                ok = true;
                for (int k = 0; ok && k < 3; ++k) {
                    if (k > 0) {
                        names[k] = nc4_section_name(file, k);
                        hsize_t vd[3] = {0};
                        if (names[k] && (h5_dims(file, names[k], vd) != 3 || vd[0] != F || vd[1] != N || vd[2] != 3)) {
                            MD_LOG_INFO("Amber NetCDF-4: '%s' of '" MD_STR_FMT "' is not of the frames and atoms of the coordinates, and is left out", names[k], MD_STR_ARG(spath));
                            names[k] = NULL;
                        }
                    }
                    if (!names[k]) continue;
                    const hid_t d = H5Dopen2(file, names[k], H5P_DEFAULT);
                    direct[k] = d >= 0 ? md_hdf5_direct_describe(file, d, arena) : NULL;
                    if (d >= 0) H5Dclose(d);
                    ok = direct[k] != NULL;
                    scale[k] = nc4_scale(file, names[k]);
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
        if (!use_time) MD_LOG_INFO("Amber NetCDF-4: the time of '" MD_STR_FMT "' does not increase from frame to frame, or is in units not known here; the frames are numbered instead", MD_STR_ARG(spath));
    }
    if (!use_time) {
        times = md_alloc(arena, F * sizeof(double));
        for (size_t f = 0; f < F; ++f) times[f] = (double)f;
        time_unit = md_unit_none();
    }

    // ## Cells, zero without a box
    float* boxes = md_alloc(arena, F * 9 * sizeof(float));
    MD_MEMSET(boxes, 0, F * 9 * sizeof(float));
    for (size_t f = 0; lengths && f < F; ++f) {
        const md_unitcell_t cell = nc4_cell(lengths + f * 3, angles ? angles + f * 3 : NULL);
        md_unitcell_A_extract_float((float(*)[3])(boxes + f * 9), &cell);
    }

    // ## The providers' state, one per variable, owned by its attribute
    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;
    md_attribute_virtual_t virt[3] = {0};
    for (int k = 0; k < 3; ++k) {
        if (!names[k]) continue;
        const size_t size = sizeof(nc4_source_t) + (size_t)direct[k]->size;
        nc4_source_t* src = md_attributes_alloc_user_data(&sys->attributes, size);
        if (!src) goto done;
        MD_MEMSET(src, 0, size);
        src->scale = scale[k];
        md_str_copy_to_char_buf(src->dataset, sizeof(src->dataset), md_str_from_cstr(names[k]));
        md_str_copy_to_char_buf(src->path, sizeof(src->path), spath);
        MD_MEMCPY(src + 1, direct[k], (size_t)direct[k]->size);
        virt[k] = (md_attribute_virtual_t){ .provider = nc4_provider, .user_data = src, .user_data_size = size };
    }

    // The frame table, for what reads the run's file as a whole: where each frame's coordinates are
    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    md_hdf5_direct_frame_spans(offsets, sizes, direct[0], F);

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
    const struct { int k; const char* leaf; md_unit_t unit; md_str_t label; } per_atom[] = {
        { 1, "atom/velocity", md_unit_div(md_unit_angstrom(), md_unit_picosecond()), MD_STR_INIT("Velocity") },
        { 2, "atom/force",    md_unit_div(kcal_per_mol, md_unit_angstrom()),         MD_STR_INIT("Force") },
    };
    bool published = true;
    for (size_t i = 0; published && i < MD_COUNTOF(per_atom); ++i) {
        const int k = per_atom[i].k;
        if (!names[k]) continue;
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, md_str_from_cstr(per_atom[i].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)F, (uint32_t)N } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = per_atom[i].unit, .label = per_atom[i].label,
            .virt = &virt[k]}) != MD_ATTRIBUTE_INVALID;
    }
    if (published && temp0) {
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("remd/temp0")),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_kelvin(), .label = MD_STR_INIT("Target temperature"),
            .data = temp0, .byte_size = F * sizeof(double)}) != MD_ATTRIBUTE_INVALID;
    }
    const struct { const double* v; const char* leaf; md_str_t label; } indices[] = {
        { repidx, "remd/repidx", MD_STR_INIT("Replica index") },
        { crdidx, "remd/crdidx", MD_STR_INIT("Coordinate index") },
    };
    for (size_t i = 0; published && i < MD_COUNTOF(indices); ++i) {
        if (!indices[i].v) continue;
        int32_t* ints = md_alloc(arena, F * sizeof(int32_t));
        for (size_t f = 0; f < F; ++f) ints[f] = (int32_t)indices[i].v[f];
        published = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, md_str_from_cstr(indices[i].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .label = indices[i].label,
            .data = ints, .byte_size = F * sizeof(int32_t)}) != MD_ATTRIBUTE_INVALID;
    }
    if (!published) {
        MD_LOG_ERROR("Amber NetCDF-4: failed to publish '" MD_STR_FMT "' as '" MD_STR_FMT "'", MD_STR_ARG(spath), MD_STR_ARG(run));
        md_attributes_remove_prefix(attributes, run);
        goto done;
    }
    if (!md_hdf5_direct_is_direct(direct[0])) {
        MD_LOG_INFO("Amber NetCDF-4: the coordinates of '" MD_STR_FMT "' are encoded in a way read through HDF5, one frame at a time", MD_STR_ARG(spath));
    }
    result = true;

done:
    md_arena_allocator_destroy(arena);
    return result;
}

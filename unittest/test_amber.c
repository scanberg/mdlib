#include "utest.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <md_amber_netcdf.h>
#include <md_amber_prmtop.h>
#include <md_amber_rst7.h>
#include <md_amber_energy.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

// ## A classic NetCDF writer, enough of one for the Amber convention: the files are written here rather than kept as
// samples, so that each layout the format allows (record or fixed frames, CDF-1 or CDF-2, a scale factor, a file cut
// short) is exercised on its own.

enum { NC_INT_ = 4, NC_FLOAT_ = 5, NC_DOUBLE_ = 6, NC_CHAR_ = 2 };

typedef struct ncw_att_t { const char* name; const char* text; double number; int type; } ncw_att_t;
typedef struct ncw_dim_t { const char* name; uint32_t len; } ncw_dim_t;   // len 0: the record dimension
typedef struct ncw_var_t {
    const char* name;
    int         type;
    int         rank;
    int         dim[3];
    ncw_att_t   att[2];
    int         num_att;
} ncw_var_t;

typedef struct ncw_t {
    int        version;             // 1 or 2
    uint32_t   numrecs;             // written to the header; 0xFFFFFFFF streaming
    uint32_t   records;             // written as data
    ncw_dim_t  dim[8];   int num_dim;
    ncw_att_t  gatt[4];  int num_gatt;
    ncw_var_t  var[10];  int num_var;
    double   (*value)(const char* var, uint32_t record, size_t index);
} ncw_t;

typedef struct buf_t { uint8_t* data; size_t size, cap; } buf_t;

static void put(buf_t* b, const void* p, size_t n) {
    if (b->size + n > b->cap) {
        b->cap = (b->size + n) * 2 + 1024;
        b->data = realloc(b->data, b->cap);
    }
    MEMCPY(b->data + b->size, p, n);
    b->size += n;
}
static void put_u32(buf_t* b, uint32_t v) { v = BSWAP32(v); put(b, &v, 4); }
static void put_u64(buf_t* b, uint64_t v) { v = BSWAP64(v); put(b, &v, 8); }
static void put_name(buf_t* b, const char* s) {
    const size_t n = strlen(s);
    put_u32(b, (uint32_t)n);
    put(b, s, n);
    const uint8_t zero[4] = {0};
    put(b, zero, (4 - n % 4) % 4);
}
static void put_value(buf_t* b, int type, double v) {
    if (type == NC_FLOAT_) { float f = (float)v; uint32_t u; MEMCPY(&u, &f, 4); put_u32(b, u); }
    else if (type == NC_DOUBLE_) { uint64_t u; MEMCPY(&u, &v, 8); put_u64(b, u); }
    else { put_u32(b, (uint32_t)(int32_t)v); }
}
static size_t type_size(int type) { return type == NC_DOUBLE_ ? 8 : type == NC_CHAR_ ? 1 : 4; }

static void put_atts(buf_t* b, const ncw_att_t* att, int n) {
    if (n == 0) { put_u32(b, 0); put_u32(b, 0); return; }
    put_u32(b, 0x0C);
    put_u32(b, (uint32_t)n);
    for (int i = 0; i < n; ++i) {
        put_name(b, att[i].name);
        if (att[i].text) {
            put_u32(b, NC_CHAR_);
            put_name(b, att[i].text);
        } else {
            put_u32(b, (uint32_t)att[i].type);
            put_u32(b, 1);
            put_value(b, att[i].type, att[i].number);
        }
    }
}

static size_t var_slab(const ncw_t* w, const ncw_var_t* v, bool* record) {
    size_t n = type_size(v->type);
    *record = false;
    for (int d = 0; d < v->rank; ++d) {
        if (w->dim[v->dim[d]].len == 0) { *record = true; continue; }
        n *= w->dim[v->dim[d]].len;
    }
    return n;
}

// Writes the file; cut, when not 0, is how many bytes to leave off its end
static bool ncw_write(const ncw_t* w, str_t path, size_t cut) {
    // The header twice: once to know its size (and so where the data begins), once for real
    uint64_t begin[10] = {0};
    buf_t h = {0};
    size_t recsize = 0, num_rec = 0;
    for (int pass = 0; pass < 2; ++pass) {
        h.size = 0;
        put(&h, "CDF", 3);
        const uint8_t ver = (uint8_t)w->version;
        put(&h, &ver, 1);
        put_u32(&h, w->numrecs);
        put_u32(&h, 0x0A);
        put_u32(&h, (uint32_t)w->num_dim);
        for (int i = 0; i < w->num_dim; ++i) { put_name(&h, w->dim[i].name); put_u32(&h, w->dim[i].len); }
        put_atts(&h, w->gatt, w->num_gatt);
        put_u32(&h, 0x0B);
        put_u32(&h, (uint32_t)w->num_var);
        for (int i = 0; i < w->num_var; ++i) {
            const ncw_var_t* v = &w->var[i];
            put_name(&h, v->name);
            put_u32(&h, (uint32_t)v->rank);
            for (int d = 0; d < v->rank; ++d) put_u32(&h, (uint32_t)v->dim[d]);
            put_atts(&h, v->att, v->num_att);
            put_u32(&h, (uint32_t)v->type);
            bool rec;
            put_u32(&h, (uint32_t)((var_slab(w, v, &rec) + 3) & ~3u));
            if (w->version == 1) put_u32(&h, (uint32_t)begin[i]); else put_u64(&h, begin[i]);
        }
        // Where everything goes: the fixed variables after the header, then the records
        uint64_t at = h.size;
        recsize = 0; num_rec = 0;
        for (int i = 0; i < w->num_var; ++i) {
            bool rec;
            const size_t slab = var_slab(w, &w->var[i], &rec);
            if (!rec) { begin[i] = at; at += (slab + 3) & ~3u; }
        }
        for (int i = 0; i < w->num_var; ++i) {
            bool rec;
            const size_t slab = var_slab(w, &w->var[i], &rec);
            if (rec) { begin[i] = at + recsize; recsize += (slab + 3) & ~3u; num_rec += 1; }
        }
    }
    if (num_rec == 1) {
        for (int i = 0; i < w->num_var; ++i) { bool rec; const size_t s = var_slab(w, &w->var[i], &rec); if (rec) recsize = s; }
    }

    buf_t b = {0};
    put(&b, h.data, h.size);
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < w->num_var; ++i) {
            const ncw_var_t* v = &w->var[i];
            bool rec;
            const size_t slab = var_slab(w, v, &rec);
            if (rec != (pass == 1)) continue;
            if (!rec) {
                const size_t n = slab / type_size(v->type);
                for (size_t k = 0; k < n; ++k) {
                    if (v->type == NC_CHAR_) { const char c = "xyzabc"[k % 6]; put(&b, &c, 1); }
                    else put_value(&b, v->type, w->value(v->name, 0, k));
                }
                const uint8_t zero[4] = {0};
                put(&b, zero, (4 - slab % 4) % 4);
            }
        }
    }
    for (uint32_t r = 0; r < w->records; ++r) {
        for (int i = 0; i < w->num_var; ++i) {
            const ncw_var_t* v = &w->var[i];
            bool rec;
            const size_t slab = var_slab(w, v, &rec);
            if (!rec) continue;
            const size_t n = slab / type_size(v->type);
            for (size_t k = 0; k < n; ++k) put_value(&b, v->type, w->value(v->name, r, k));
            if (num_rec > 1) { const uint8_t zero[4] = {0}; put(&b, zero, (4 - slab % 4) % 4); }
        }
    }
    md_file_t f = {0};
    bool ok = md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE);
    if (ok) {
        md_file_write(f, b.data, b.size - cut);
        md_file_close(&f);
    }
    free(h.data);
    free(b.data);
    return ok;
}

#define NATOMS 5

static double coord(uint32_t f, size_t i, int d) { return 10.0 + (double)i + 0.5 * d + 0.25 * f; }

static double traj_value(const char* var, uint32_t r, size_t k) {
    if (!strcmp(var, "coordinates")) return coord(r, k / 3, (int)(k % 3));
    if (!strcmp(var, "velocities"))  return (double)(k + 1) / 20.455;     // stored in AKMA, scaled on reading
    if (!strcmp(var, "forces"))      return -(double)k;
    if (!strcmp(var, "time"))        return 2.0 * (r + 1);
    if (!strcmp(var, "cell_lengths")) return 30.0 + r + k;
    if (!strcmp(var, "cell_angles")) return 90.0;
    if (!strcmp(var, "temp0"))       return 300.0 + 10.0 * r;
    return 0.0;
}

// The fixed frame variant: a variable along a fixed frame dimension is one array; k counts over all of it
static double fixed_value(const char* var, uint32_t r, size_t k) {
    (void)r;
    if (!strcmp(var, "coordinates")) return coord((uint32_t)(k / (NATOMS * 3)), (k / 3) % NATOMS, (int)(k % 3));
    return 0.0;
}

static ncw_t amber_trajectory(int version, uint32_t frames) {
    ncw_t w = { .version = version, .numrecs = frames, .records = frames, .value = traj_value };
    w.dim[w.num_dim++] = (ncw_dim_t){ "frame", 0 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "spatial", 3 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "atom", NATOMS };
    w.dim[w.num_dim++] = (ncw_dim_t){ "cell_spatial", 3 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "cell_angular", 3 };
    w.gatt[w.num_gatt++] = (ncw_att_t){ "Conventions", "AMBER" };
    w.gatt[w.num_gatt++] = (ncw_att_t){ "ConventionVersion", "1.0" };
    w.var[w.num_var++] = (ncw_var_t){ "spatial", NC_CHAR_, 1, {1}, {{0}}, 0 };
    w.var[w.num_var++] = (ncw_var_t){ "time", NC_FLOAT_, 1, {0}, {{ "units", "picosecond" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "coordinates", NC_FLOAT_, 3, {0, 2, 1}, {{ "units", "angstrom" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "cell_lengths", NC_DOUBLE_, 2, {0, 3}, {{ "units", "angstrom" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "cell_angles", NC_DOUBLE_, 2, {0, 4}, {{ "units", "degree" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "velocities", NC_FLOAT_, 3, {0, 2, 1}, {{ "units", "angstrom/picosecond" }, { "scale_factor", NULL, 20.455, NC_DOUBLE_ }}, 2 };
    w.var[w.num_var++] = (ncw_var_t){ "forces", NC_FLOAT_, 3, {0, 2, 1}, {{ "units", "kilocalorie/mole/angstrom" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "temp0", NC_DOUBLE_, 1, {0}, {{ "units", "kelvin" }}, 1 };
    return w;
}

static void check_positions(int* utest_result, md_system_t* sys, str_t run, uint32_t frames, md_allocator_i* arena) {
    const str_t paths[] = { STR_LIT("atom/position") };
    md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 1, arena);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t s = { .alloc = arena };
    md_system_state_init(&s, NATOMS);
    for (uint32_t f = 0; f < frames; ++f) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &s));
        for (size_t i = 0; i < NATOMS; ++i) {
            EXPECT_EQ((float)coord(f, i, 0), s.xyz[i].x);
            EXPECT_EQ((float)coord(f, i, 1), s.xyz[i].y);
            EXPECT_EQ((float)coord(f, i, 2), s.xyz[i].z);
        }
    }
    md_system_extract_end(ex);
}

UTEST(amber_netcdf, trajectory_with_records) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber.nc");
    for (int version = 1; version <= 2; ++version) {
        ncw_t w = amber_trajectory(version, 4);
        ASSERT_TRUE(ncw_write(&w, path, 0));
        EXPECT_EQ(MD_AMBER_NETCDF_TRAJECTORY, md_amber_netcdf_file_kind(path));
        EXPECT_EQ((size_t)NATOMS, md_amber_netcdf_file_atom_count(path));

        md_system_t sys = { .alloc = arena };
        ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, path, STR_LIT("run/a"), 0));
        const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/a/time"));
        ASSERT_TRUE(time != NULL);
        ASSERT_EQ(4u, time->format.shape[0]);
        EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
        EXPECT_EQ(2.0, ((const double*)time->data)[0]);
        EXPECT_EQ(8.0, ((const double*)time->data)[3]);
        check_positions(utest_result, &sys, STR_LIT("run/a"), 4, arena);

        // The cell of frame 2: lengths 32, 33, 34
        const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/a/unitcell"));
        ASSERT_TRUE(cell != NULL);
        const float* box = (const float*)cell->data + 2 * 9;
        EXPECT_NEAR(32.0f, box[0], 1e-4f);
        EXPECT_NEAR(33.0f, box[4], 1e-4f);
        EXPECT_NEAR(34.0f, box[8], 1e-4f);

        // Velocities with their scale factor, in Angstrom/ps; forces as they are
        const md_attribute_t* vel = md_attributes_find(&sys.attributes, STR_LIT("run/a/atom/velocity"));
        const md_attribute_t* frc = md_attributes_find(&sys.attributes, STR_LIT("run/a/atom/force"));
        ASSERT_TRUE(vel && frc);
        float v[3];
        ASSERT_EQ(3u, md_attribute_extract_f32(v, 3, vel, md_attribute_slice_2(1, 2), md_unit_none()));
        EXPECT_NEAR(7.0f, v[0], 1e-4f);     // k = 6: (6 + 1) / 20.455 * 20.455
        EXPECT_NEAR(9.0f, v[2], 1e-4f);
        ASSERT_EQ(3u, md_attribute_extract_f32(v, 3, frc, md_attribute_slice_2(3, 1), md_unit_none()));
        EXPECT_EQ(-3.0f, v[0]);

        const md_attribute_t* temp0 = md_attributes_find(&sys.attributes, STR_LIT("run/a/remd/temp0"));
        ASSERT_TRUE(temp0 != NULL);
        EXPECT_EQ(330.0, ((const double*)temp0->data)[3]);

        // The first frame as a state, as a topology takes its coordinates
        md_system_state_t s = { .alloc = arena };
        ASSERT_TRUE(md_amber_netcdf_state_init_from_file(&s, path, NATOMS));
        EXPECT_EQ((float)coord(0, 4, 2), s.xyz[4].z);
        EXPECT_NEAR(30.0f, s.unitcell.x, 1e-4f);
        EXPECT_FALSE(md_amber_netcdf_state_init_from_file(&s, path, NATOMS + 1));
        md_system_free(&sys);
    }
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// A file still being written: the last record is incomplete, and with the header's record count unknown
UTEST(amber_netcdf, trajectory_cut_short) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber_cut.nc");
    for (int streaming = 0; streaming < 2; ++streaming) {
        ncw_t w = amber_trajectory(2, 5);
        if (streaming) w.numrecs = 0xFFFFFFFFu;
        ASSERT_TRUE(ncw_write(&w, path, 20));
        md_system_t sys = { .alloc = arena };
        ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, path, STR_LIT("run/a"), 0));
        const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/a/time"));
        ASSERT_TRUE(time != NULL);
        EXPECT_EQ(4u, time->format.shape[0]);
        check_positions(utest_result, &sys, STR_LIT("run/a"), 4, arena);
        md_system_free(&sys);
    }
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// Frames along a fixed dimension rather than records, no time, no cell: CDF-1 as some converters write it
UTEST(amber_netcdf, trajectory_fixed_frames) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber_fixed.nc");
    ncw_t w = { .version = 1, .numrecs = 0, .records = 0, .value = fixed_value };
    w.dim[w.num_dim++] = (ncw_dim_t){ "frame", 3 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "atom", NATOMS };
    w.dim[w.num_dim++] = (ncw_dim_t){ "spatial", 3 };
    w.gatt[w.num_gatt++] = (ncw_att_t){ "Conventions", "AMBER" };
    w.var[w.num_var++] = (ncw_var_t){ "coordinates", NC_FLOAT_, 3, {0, 1, 2}, {{ "units", "angstrom" }}, 1 };
    ASSERT_TRUE(ncw_write(&w, path, 0));

    md_system_t sys = { .alloc = arena };
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, path, STR_LIT("run/a"), 0));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/a/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_EQ(3u, time->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_none()));
    EXPECT_EQ(2.0, ((const double*)time->data)[2]);
    // No box: a cell of zeros, which every run has
    const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/a/unitcell"));
    ASSERT_TRUE(cell != NULL);
    EXPECT_EQ(0.0f, ((const float*)cell->data)[0]);
    check_positions(utest_result, &sys, STR_LIT("run/a"), 3, arena);
    md_system_free(&sys);
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

static double restart_value(const char* var, uint32_t r, size_t k) {
    (void)r;
    if (!strcmp(var, "coordinates")) return coord(7, k / 3, (int)(k % 3));
    if (!strcmp(var, "cell_lengths")) return 40.0;
    if (!strcmp(var, "cell_angles")) return 109.4712206;
    if (!strcmp(var, "time")) return 123.0;
    return 0.0;
}

UTEST(amber_netcdf, restart_and_what_is_not_amber) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber.ncrst");
    ncw_t w = { .version = 2, .value = restart_value };
    w.dim[w.num_dim++] = (ncw_dim_t){ "spatial", 3 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "atom", NATOMS };
    w.dim[w.num_dim++] = (ncw_dim_t){ "cell_spatial", 3 };
    w.dim[w.num_dim++] = (ncw_dim_t){ "cell_angular", 3 };
    w.gatt[w.num_gatt++] = (ncw_att_t){ "Conventions", "AMBERRESTART" };
    w.var[w.num_var++] = (ncw_var_t){ "time", NC_DOUBLE_, 0, {0}, {{ "units", "picosecond" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "coordinates", NC_DOUBLE_, 2, {1, 0}, {{ "units", "angstrom" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "cell_lengths", NC_DOUBLE_, 1, {2}, {{ "units", "angstrom" }}, 1 };
    w.var[w.num_var++] = (ncw_var_t){ "cell_angles", NC_DOUBLE_, 1, {3}, {{ "units", "degree" }}, 1 };
    ASSERT_TRUE(ncw_write(&w, path, 0));
    EXPECT_EQ(MD_AMBER_NETCDF_RESTART, md_amber_netcdf_file_kind(path));

    md_system_state_t s = { .alloc = arena };
    ASSERT_TRUE(md_amber_netcdf_state_init_from_file(&s, path, 0));
    EXPECT_EQ((size_t)NATOMS, s.num_atoms);
    EXPECT_EQ((float)coord(7, 3, 1), s.xyz[3].y);
    EXPECT_TRUE(md_system_state_has_unitcell(&s));
    // A restart is not a run
    md_system_t sys = { .alloc = arena };
    EXPECT_FALSE(md_amber_netcdf_system_publish_run(&sys, path, STR_LIT("run/a"), 0));

    // NetCDF without the convention, NetCDF-4 and something else altogether
    w.gatt[0].text = "CF-1.6";
    ASSERT_TRUE(ncw_write(&w, path, 0));
    EXPECT_EQ(MD_AMBER_NETCDF_NONE, md_amber_netcdf_file_kind(path));
    EXPECT_FALSE(md_amber_netcdf_state_init_from_file(&s, path, 0));

    md_file_t f = {0};
    ASSERT_TRUE(md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE));
    md_file_write(f, "\211HDF\r\n\032\n and then some", 24);
    md_file_close(&f);
    EXPECT_EQ(MD_AMBER_NETCDF_HDF5, md_amber_netcdf_file_kind(path));
    EXPECT_FALSE(md_amber_netcdf_system_publish_run(&sys, path, STR_LIT("run/a"), 0));

    md_system_free(&sys);
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// ## Topology

static bool write_text(str_t path, const char* text) {
    md_file_t f = {0};
    if (!md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    md_file_write(f, text, strlen(text));
    md_file_close(&f);
    return true;
}

// Methanol and a TIP4P-Ew water with its extra point, periodic. Charges in Amber units (e * 18.2223).
static const char* prmtop_text =
    "%VERSION  VERSION_STAMP = V0001.000  DATE = 11/02/26  12:00:00\n"
    "%FLAG TITLE\n"
    "%FORMAT(20a4)\n"
    "methanol in water\n"
    "%FLAG POINTERS\n"
    "%FORMAT(10I8)\n"
    "      10       4       5       2      0       0       0       0       0       0\n"
    "       0       2       2       0       0       2       0       0       4       0\n"
    "       0       0       0       0       0       0       0       1       6       0\n"
    "       1\n"
    "%FLAG ATOM_NAME\n"
    "%FORMAT(20a4)\n"
    "C1  H11 H12 H13 O1  HO1 O   H1  H2  EPW \n"
    "%FLAG CHARGE\n"
    "%FORMAT(5E16.8)\n"
    "  2.18223000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00 -1.18223000E+01\n"
    "  9.64000000E+00  0.00000000E+00  9.64000000E+00  9.64000000E+00 -1.92800000E+01\n"
    "%FLAG ATOMIC_NUMBER\n"
    "%FORMAT(10I8)\n"
    "       6       1       1       1       8       1       8       1       1       0\n"
    "%FLAG MASS\n"
    "%FORMAT(5E16.8)\n"
    "  1.20100000E+01  1.00800000E+00  1.00800000E+00  1.00800000E+00  1.60000000E+01\n"
    "  1.00800000E+00  1.60000000E+01  1.00800000E+00  1.00800000E+00  0.00000000E+00\n"
    "%FLAG RESIDUE_LABEL\n"
    "%FORMAT(20a4)\n"
    "MOH WAT \n"
    "%FLAG RESIDUE_POINTER\n"
    "%FORMAT(10I8)\n"
    "       1       7\n"
    "%FLAG AMBER_ATOM_TYPE\n"
    "%FORMAT(20a4)\n"
    "CT  HC  HC  HC  OH  HO  OW  HW  HW  EP  \n"
    "%FLAG BONDS_INC_HYDROGEN\n"
    "%FORMAT(10I8)\n"
    "       0       3       1       0       6       1       0       9       1      12\n"
    "      15       2      18      21       3      18      24       3      18      27\n"
    "       4\n"
    "%FLAG BONDS_WITHOUT_HYDROGEN\n"
    "%FORMAT(10I8)\n"
    "       0      12       2      21      24       3\n"
    "%FLAG RESIDUE_CHAINID\n"
    "%FORMAT(20a4)\n"
    "A   B   \n"
    "%FLAG BOX_DIMENSIONS\n"
    "%FORMAT(5E16.8)\n"
    "  9.00000000E+01  2.00000000E+01  2.10000000E+01  2.20000000E+01\n";

UTEST(amber_prmtop, builds_a_system_without_coordinates) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber.prmtop");
    ASSERT_TRUE(write_text(path, prmtop_text));
    EXPECT_TRUE(md_amber_prmtop_file_is_topology(path));
    EXPECT_EQ(10u, md_amber_prmtop_file_atom_count(path));

    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_amber_prmtop_system_init_from_file(&sys, &state, path));
    EXPECT_EQ(10u, sys.atom.count);
    EXPECT_EQ(2u, sys.residue.count);
    EXPECT_FALSE(md_system_state_has_coords(&state));
    EXPECT_TRUE(md_system_state_has_unitcell(&state));
    EXPECT_NEAR(20.0f, state.unitcell.x, 1e-4f);

    // Elements, and the extra point as a virtual site with its bond left out: 8 of the 9 bonds
    EXPECT_EQ(6, sys.atom.type.z[sys.atom.type_idx[0]]);
    EXPECT_EQ(8, sys.atom.type.z[sys.atom.type_idx[6]]);
    EXPECT_EQ(MD_PARTICLE_VIRTUAL_SITE, md_atom_particle_kind(&sys.atom, 9));
    EXPECT_EQ(8u, sys.bond.count);
    EXPECT_TRUE(str_eq_cstr(md_label16_str(&sys.atom.type.ff_type[sys.atom.type_idx[0]]), "CT"));

    const md_attribute_t* charge = md_attributes_find(&sys.attributes, STR_LIT("atom/charge"));
    ASSERT_TRUE(charge != NULL);
    EXPECT_NEAR(0.11975f, ((const float*)charge->data)[0], 1e-4f);
    EXPECT_NEAR(-1.05804f, ((const float*)charge->data)[9], 1e-4f);
    EXPECT_TRUE(str_eq_cstr(md_system_residue_name(&sys, 1), "WAT"));

    // CHAMBER: the same blocks, and CHARMM terms that are not read
    char chamber[8192];
    snprintf(chamber, sizeof(chamber), "%s%%FLAG CHARMM_UREY_BRADLEY_COUNT\n%%FORMAT(2I8)\n       0       0\n", prmtop_text);
    ASSERT_TRUE(write_text(path, chamber));
    md_system_t sys2 = { .alloc = arena };
    ASSERT_TRUE(md_amber_prmtop_system_init_from_file(&sys2, &state, path));
    EXPECT_EQ(10u, sys2.atom.count);

    md_system_free(&sys);
    md_system_free(&sys2);
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

UTEST(amber_prmtop, supplements_a_system_of_the_same_atoms) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber_supp.prmtop");
    ASSERT_TRUE(write_text(path, prmtop_text));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_amber_prmtop_system_init_from_file(&sys, &state, path));
    md_bond_data_clear(&sys.bond);
    md_bond_build_connectivity(&sys.bond, sys.atom.count, arena);
    EXPECT_EQ(0u, sys.bond.count);
    ASSERT_TRUE(md_amber_prmtop_system_supplement_from_file(&sys, path));
    EXPECT_EQ(8u, sys.bond.count);

    // Other atoms: refused, nothing changed
    char other[8192];
    snprintf(other, sizeof(other), "%s", prmtop_text);
    char* name = strstr(other, "HO1 ");
    MEMCPY(name, "XX1 ", 4);
    ASSERT_TRUE(write_text(path, other));
    EXPECT_FALSE(md_amber_prmtop_system_supplement_from_file(&sys, path));
    EXPECT_EQ(8u, sys.bond.count);

    md_system_free(&sys);
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// ## ASCII restart

UTEST(amber_rst7, coordinates_velocities_and_box) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t path = STR_LIT("md_unittest_amber.rst7");
    // 3 atoms: two lines of coordinates, two of velocities, a box
    const char* text =
        "default_name\n"
        "     3  0.1000000E+02\n"
        "   1.0000000   2.0000000   3.0000000   4.0000000   5.0000000   6.0000000\n"
        "  -7.0000000-100.0000000   9.0000000\n"
        "   0.1000000   0.2000000   0.3000000   0.4000000   0.5000000   0.6000000\n"
        "   0.7000000   0.8000000   0.9000000\n"
        "  30.0000000  31.0000000  32.0000000  90.0000000  90.0000000  90.0000000\n";
    ASSERT_TRUE(write_text(path, text));
    EXPECT_TRUE(md_amber_rst7_file_is_restart(path));
    EXPECT_EQ(3u, md_amber_rst7_file_atom_count(path));
    md_system_state_t s = { .alloc = arena };
    ASSERT_TRUE(md_amber_rst7_state_init_from_file(&s, path, 3));
    EXPECT_EQ(3u, s.num_atoms);
    EXPECT_EQ(-100.0f, s.xyz[2].y);   // fields of 12 that run together
    EXPECT_NEAR(31.0f, s.unitcell.y, 1e-4f);
    EXPECT_FALSE(md_amber_rst7_state_init_from_file(&s, path, 4));

    // Coordinates alone: no box
    const char* bare =
        "bare\n"
        "     2\n"
        "   1.0000000   2.0000000   3.0000000   4.0000000   5.0000000   6.0000000\n";
    ASSERT_TRUE(write_text(path, bare));
    ASSERT_TRUE(md_amber_rst7_state_init_from_file(&s, path, 0));
    EXPECT_FALSE(md_system_state_has_unitcell(&s));
    EXPECT_EQ(6.0f, s.xyz[1].z);

    // Not a restart
    ASSERT_TRUE(write_text(path, "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"));
    EXPECT_FALSE(md_amber_rst7_file_is_restart(path));
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

// ## Energies

static const char* mdout_text =
    "          -------------------------------------------------------\n"
    "          Amber 24 PMEMD                              2024\n"
    "          -------------------------------------------------------\n"
    "| PMEMD implementation of SANDER, Release 24\n"
    "File Assignments:\n"
    "|   MDIN: mdin\n"
    "|  MDOUT: mdout\n"
    "--------------------------------------------------------------------------------\n"
    "   4.  RESULTS\n"
    "--------------------------------------------------------------------------------\n"
    "\n"
    " NSTEP =     1000   TIME(PS) =       2.000  TEMP(K) =   300.12  PRESS =     0.0\n"
    " Etot   =     -1000.0000  EKtot   =       100.0000  EPtot      =     -1100.0000\n"
    " BOND   =         2.4475  ANGLE   =         2.8681  DIHED      =         3.3755\n"
    " 1-4 NB =         0.5154  1-4 EEL =       -19.8735  VDWAALS    =        -0.2550\n"
    " EELEC  =       -11.2733  EHBOND  =         0.0000  RESTRAINT  =         0.0000\n"
    " ------------------------------------------------------------------------------\n"
    "\n"
    "   KE Trans =     0.0000   KE Rot =     0.0000   C.O.M. Vel =    0.000000\n"
    "\n"
    " NSTEP =     2000   TIME(PS) =       4.000  TEMP(K) =   301.12  PRESS =     0.0\n"
    " Etot   =     -1001.0000  EKtot   =       101.0000  EPtot      =     -1102.0000\n"
    " BOND   =         2.4475  ANGLE   =         2.8681  DIHED      =         3.3755\n"
    " 1-4 NB =         0.5154  1-4 EEL =       -19.8735  VDWAALS    =        -0.2550\n"
    " EELEC  =       -11.2733  EHBOND  =         0.0000  RESTRAINT  =         0.0000\n"
    " ------------------------------------------------------------------------------\n"
    "\n"
    "      A V E R A G E S   O V E R       2 S T E P S\n"
    "\n"
    " NSTEP =     2000   TIME(PS) =       4.000  TEMP(K) =   300.62  PRESS =     0.0\n"
    " Etot   =     -1000.5000  EKtot   =       100.5000  EPtot      =     -1101.0000\n"
    " ------------------------------------------------------------------------------\n"
    "\n"
    "   5.  TIMINGS\n";

static double two_frames(const char* var, uint32_t r, size_t k) {
    if (!strcmp(var, "time")) return 2.0 * (r + 1);
    return traj_value(var, r, k);
}

UTEST(amber_energy, mdout_and_mden_along_the_run) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t nc = STR_LIT("md_unittest_amber_energy.nc");
    const str_t out = STR_LIT("md_unittest_amber.mdout");
    ncw_t w = amber_trajectory(2, 2);
    w.value = two_frames;
    ASSERT_TRUE(ncw_write(&w, nc, 0));
    ASSERT_TRUE(write_text(out, mdout_text));
    EXPECT_TRUE(md_amber_energy_file_is_mdout(out));
    EXPECT_FALSE(md_amber_energy_file_is_mden(out));

    md_system_t sys = { .alloc = arena };
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, nc, STR_LIT("run/a"), 0));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, out, STR_LIT("run/a")));
    const md_attribute_t* ep = md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/eptot"));
    const md_attribute_t* nb = md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/1_4_nb"));
    const md_attribute_t* tk = md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/temp_k"));
    ASSERT_TRUE(ep && nb && tk);
    EXPECT_EQ(2u, ep->format.shape[0]);     // the averages are not a record
    EXPECT_EQ(-1102.0, ((const double*)ep->data)[1]);
    EXPECT_TRUE(md_unit_equal(tk->unit, md_unit_kelvin()));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/ke_trans")) == NULL);

    // mden, of the same records: under the same names
    const char* mden =
        "L0  Nsteps           time(ps)         Etot             EKinetic\n"
        "L1  Temp             T_solute         T_solv           Pres_scal_solu\n"
        "L2  E_pot            E_vdw            E_el\n"
        "L0     1000  0.2000000000E+01  -.1000000000E+04  0.1000000000E+03\n"
        "L1 0.3001200000E+03  0.3000000000E+03  0.0000000000E+00  0.1000000000E+01\n"
        "L2 -.1100000000E+04  -.2550000000E+00  -.1127330000E+02\n"
        "L0     2000  0.4000000000E+01  -.1001000000E+04  0.1010000000E+03\n"
        "L1 0.3011200000E+03  0.3000000000E+03  0.0000000000E+00  0.1000000000E+01\n"
        "L2 -.1102000000E+04  -.2550000000E+00  -.1127330000E+02\n";
    ASSERT_TRUE(write_text(out, mden));
    EXPECT_TRUE(md_amber_energy_file_is_mden(out));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, out, STR_LIT("run/a")));
    ep = md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/eptot"));
    ASSERT_TRUE(ep != NULL);
    EXPECT_EQ(-1102.0, ((const double*)ep->data)[1]);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/a/amber_energy/vdwaals")) != NULL);

    // A run the file does not cover: refused, nothing changed
    md_system_t other = { .alloc = arena };
    ncw_t w3 = amber_trajectory(2, 3);
    ASSERT_TRUE(ncw_write(&w3, nc, 0));
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&other, nc, STR_LIT("run/b"), 0));
    EXPECT_FALSE(md_amber_energy_system_supplement_from_file(&other, out, STR_LIT("run/b")));
    EXPECT_TRUE(md_attributes_find(&other.attributes, STR_LIT("run/b/amber_energy/time")) == NULL);

    md_system_free(&sys);
    md_system_free(&other);
    remove(nc.ptr);
    remove(out.ptr);
    md_arena_allocator_destroy(arena);
}

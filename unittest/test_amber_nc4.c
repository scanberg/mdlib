#include "utest.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <hdf5.h>

#include <md_amber_netcdf.h>
#include <md_system.h>
#include <md_run.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_unit.h>

// Amber NetCDF in the NetCDF-4 encoding: HDF5 files laid out as the NetCDF library lays them out, written here with
// the HDF5 API. Each storage the reader decodes itself (contiguous, chunks with shuffle and deflate, several frames a
// chunk, big endian, cpptraj's integer compression) and two it leaves to HDF5 (a chunk which splits a frame, a filter
// it does not know), each held against the same values.

#define NC4_ATOMS  37
#define NC4_FRAMES 5

static double nc4_coord(size_t f, size_t i, int d) { return 10.0 + 0.5 * (double)i + 0.25 * d + 1.5 * (double)f; }
static double nc4_vel(size_t f, size_t i, int d)   { return 0.01 * (double)(i + 1) - 0.02 * d + 0.001 * (double)f; }

typedef struct nc4_case_t {
    const char* name;
    bool     contiguous;
    bool     big_endian;
    bool     shuffle;
    int      deflate;           // level, 0 for none
    bool     fletcher;          // a filter the reader does not decode itself
    hsize_t  chunk_frames;
    hsize_t  chunk_atoms;       // NC4_ATOMS for whole frames
    bool     icompress;         // compressedpos as int32 with icompressfac
} nc4_case_t;

static void nc4_text_attr(hid_t obj, const char* name, const char* text) {
    const hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, strlen(text));
    H5Tset_strpad(type, H5T_STR_NULLTERM);
    const hid_t space = H5Screate(H5S_SCALAR);
    const hid_t a = H5Acreate2(obj, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, type, text);
    H5Aclose(a);
    H5Sclose(space);
    H5Tclose(type);
}

static void nc4_number_attr(hid_t obj, const char* name, hid_t type, double value) {
    const hid_t space = H5Screate(H5S_SCALAR);
    const hid_t a = H5Acreate2(obj, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(a);
    H5Sclose(space);
}

static hid_t nc4_dataset(hid_t file, const char* name, hid_t file_type, int rank, const hsize_t* dims, const nc4_case_t* c, const void* data, hid_t mem_type) {
    hsize_t max[3];
    for (int i = 0; i < rank; ++i) max[i] = dims[i];
    if (!c->contiguous && rank == 3) max[0] = H5S_UNLIMITED;
    const hid_t space = H5Screate_simple(rank, dims, max);
    const hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    if (rank == 3 && !c->contiguous) {
        const hsize_t chunk[3] = { c->chunk_frames, c->chunk_atoms, 3 };
        H5Pset_chunk(dcpl, 3, chunk);
        if (c->fletcher) H5Pset_fletcher32(dcpl);
        if (c->shuffle)  H5Pset_shuffle(dcpl);
        if (c->deflate)  H5Pset_deflate(dcpl, (unsigned)c->deflate);
    }
    const hid_t d = H5Dcreate2(file, name, file_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(d, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Pclose(dcpl);
    H5Sclose(space);
    return d;
}

static bool nc4_write(md_str_t path, const nc4_case_t* c, bool restart) {
    const hid_t file = H5Fcreate(path.ptr, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return false;
    nc4_text_attr(file, "Conventions", restart ? "AMBERRESTART" : "AMBER");
    nc4_text_attr(file, "ConventionVersion", "1.0");
    nc4_text_attr(file, "program", "md_unittest");

    const size_t F = restart ? 1 : NC4_FRAMES;
    static float  pos[NC4_FRAMES][NC4_ATOMS][3];
    static float  vel[NC4_FRAMES][NC4_ATOMS][3];
    static int    ipos[NC4_FRAMES][NC4_ATOMS][3];
    for (size_t f = 0; f < F; ++f) for (size_t i = 0; i < NC4_ATOMS; ++i) for (int d = 0; d < 3; ++d) {
        pos[f][i][d]  = (float)nc4_coord(f, i, d);
        vel[f][i][d]  = (float)(nc4_vel(f, i, d) / 20.455);
        ipos[f][i][d] = (int)(nc4_coord(f, i, d) * 1000.0 + 0.5);
    }
    const hid_t ftype = c->big_endian ? H5T_IEEE_F32BE : H5T_IEEE_F32LE;

    if (restart) {
        const hsize_t dims[2] = { NC4_ATOMS, 3 };
        hid_t d = nc4_dataset(file, "coordinates", H5T_IEEE_F64LE, 2, dims, c, pos, H5T_NATIVE_FLOAT);
        nc4_text_attr(d, "units", "angstrom");
        H5Dclose(d);
        const hsize_t three = 3;
        const double lengths[3] = { 40.0, 41.0, 42.0 }, angles[3] = { 90.0, 90.0, 90.0 };
        d = nc4_dataset(file, "cell_lengths", H5T_IEEE_F64LE, 1, &three, c, lengths, H5T_NATIVE_DOUBLE); H5Dclose(d);
        d = nc4_dataset(file, "cell_angles",  H5T_IEEE_F64LE, 1, &three, c, angles,  H5T_NATIVE_DOUBLE); H5Dclose(d);
        H5Fclose(file);
        return true;
    }

    const hsize_t dims[3] = { F, NC4_ATOMS, 3 };
    hid_t d;
    if (c->icompress) {
        d = nc4_dataset(file, "compressedpos", H5T_STD_I32LE, 3, dims, c, ipos, H5T_NATIVE_INT);
        nc4_number_attr(d, "icompressfac", H5T_IEEE_F64LE, 1000.0);
    } else {
        d = nc4_dataset(file, "coordinates", ftype, 3, dims, c, pos, H5T_NATIVE_FLOAT);
    }
    nc4_text_attr(d, "units", "angstrom");
    H5Dclose(d);
    d = nc4_dataset(file, "velocities", ftype, 3, dims, c, vel, H5T_NATIVE_FLOAT);
    nc4_number_attr(d, "scale_factor", H5T_IEEE_F64LE, 20.455);
    H5Dclose(d);

    double time[NC4_FRAMES], lengths[NC4_FRAMES][3], angles[NC4_FRAMES][3], temp0[NC4_FRAMES];
    for (size_t f = 0; f < F; ++f) {
        time[f] = 2.0 * (double)(f + 1);
        temp0[f] = 300.0 + 10.0 * (double)f;
        for (int k = 0; k < 3; ++k) { lengths[f][k] = 30.0 + (double)f + k; angles[f][k] = 90.0; }
    }
    const nc4_case_t plain = { 0 };
    const hsize_t d1 = F, d2[2] = { F, 3 };
    d = nc4_dataset(file, "time", H5T_IEEE_F32LE, 1, &d1, &plain, time, H5T_NATIVE_DOUBLE);
    nc4_text_attr(d, "units", "picosecond");
    H5Dclose(d);
    d = nc4_dataset(file, "cell_lengths", H5T_IEEE_F64LE, 2, d2, &plain, lengths, H5T_NATIVE_DOUBLE); H5Dclose(d);
    d = nc4_dataset(file, "cell_angles",  H5T_IEEE_F64LE, 2, d2, &plain, angles,  H5T_NATIVE_DOUBLE); H5Dclose(d);
    d = nc4_dataset(file, "temp0",        H5T_IEEE_F64LE, 1, &d1, &plain, temp0,  H5T_NATIVE_DOUBLE); H5Dclose(d);
    H5Fclose(file);
    return true;
}

UTEST(amber_netcdf4, trajectories_in_every_storage) {
    const nc4_case_t cases[] = {
        { "deflate+shuffle",   .shuffle = true, .deflate = 4, .chunk_frames = 1, .chunk_atoms = NC4_ATOMS },
        { "chunked raw",       .chunk_frames = 1, .chunk_atoms = NC4_ATOMS },
        { "contiguous",        .contiguous = true },
        { "big endian",        .big_endian = true, .deflate = 1, .chunk_frames = 1, .chunk_atoms = NC4_ATOMS },
        { "two frames a chunk",.shuffle = true, .deflate = 2, .chunk_frames = 2, .chunk_atoms = NC4_ATOMS },
        { "icompress",         .icompress = true, .shuffle = true, .deflate = 4, .chunk_frames = 1, .chunk_atoms = NC4_ATOMS },
        { "split frame",       .deflate = 4, .chunk_frames = 1, .chunk_atoms = 16 },
        { "fletcher32",        .fletcher = true, .deflate = 4, .chunk_frames = 1, .chunk_atoms = NC4_ATOMS },
    };
    const md_str_t path = MD_STR_LIT("md_unittest_amber4.nc");
    for (size_t ci = 0; ci < MD_COUNTOF(cases); ++ci) {
        const nc4_case_t* c = &cases[ci];
        md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
        ASSERT_TRUE(nc4_write(path, c, false));
        EXPECT_EQ(MD_AMBER_NETCDF_TRAJECTORY, md_amber_netcdf_file_kind(path));
        EXPECT_EQ((size_t)NC4_ATOMS, md_amber_netcdf_file_atom_count(path));

        md_system_t sys = { .alloc = arena };
        ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, path, MD_STR_LIT("run/a"), 0));
        const md_attribute_t* time = md_attributes_find(&sys.attributes, MD_STR_LIT("run/a/time"));
        ASSERT_TRUE(time != NULL);
        ASSERT_EQ((size_t)NC4_FRAMES, time->format.shape[0]);
        EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
        EXPECT_EQ(4.0, ((const double*)time->data)[1]);

        const md_attribute_t* cell = md_attributes_find(&sys.attributes, MD_STR_LIT("run/a/unitcell"));
        ASSERT_TRUE(cell != NULL);
        EXPECT_NEAR(33.0f, ((const float*)cell->data)[3 * 9 + 0], 1e-4f);
        EXPECT_NEAR(35.0f, ((const float*)cell->data)[3 * 9 + 8], 1e-4f);
        const md_attribute_t* temp0 = md_attributes_find(&sys.attributes, MD_STR_LIT("run/a/remd/temp0"));
        ASSERT_TRUE(temp0 != NULL);
        EXPECT_EQ(320.0, ((const double*)temp0->data)[2]);

        // Every frame, read back in reverse so that no chunk is met in the order it was written
        const md_str_t paths[] = { MD_STR_LIT("atom/position") };
        md_system_extract_t* ex = md_system_extract_begin(&sys, MD_STR_LIT("run/a"), paths, 1, arena);
        ASSERT_TRUE(ex != NULL);
        md_system_state_t s = { .alloc = arena };
        md_system_state_init(&s, NC4_ATOMS);
        int wrong = 0;
        for (size_t f = NC4_FRAMES; f-- > 0;) {
            ASSERT_TRUE(md_system_extract_frame(ex, (uint32_t)f, &s));
            for (size_t i = 0; i < NC4_ATOMS; ++i) {
                const float tol = c->icompress ? 1e-3f : 0.0f;
                wrong += fabsf(s.xyz[i].x - (float)nc4_coord(f, i, 0)) > tol;
                wrong += fabsf(s.xyz[i].y - (float)nc4_coord(f, i, 1)) > tol;
                wrong += fabsf(s.xyz[i].z - (float)nc4_coord(f, i, 2)) > tol;
            }
        }
        md_system_extract_end(ex);
        EXPECT_EQ(0, wrong);
        if (wrong) printf("%s: %d values wrong\n", c->name, wrong);

        // Velocities with their scale factor; one atom of a frame alone
        const md_attribute_t* vel = md_attributes_find(&sys.attributes, MD_STR_LIT("run/a/atom/velocity"));
        ASSERT_TRUE(vel != NULL);
        float v[3];
        ASSERT_EQ(3u, md_attribute_extract_f32(v, 3, vel, md_attribute_slice_2(3, 20), md_unit_none()));
        EXPECT_NEAR((float)nc4_vel(3, 20, 0), v[0], 1e-5f);
        EXPECT_NEAR((float)nc4_vel(3, 20, 2), v[2], 1e-5f);

        // The first frame as a state
        md_system_state_t first = { .alloc = arena };
        ASSERT_TRUE(md_amber_netcdf_state_init_from_file(&first, path, NC4_ATOMS));
        EXPECT_EQ((float)nc4_coord(0, 36, 1), first.xyz[36].y);
        EXPECT_FALSE(md_amber_netcdf_state_init_from_file(&first, path, NC4_ATOMS + 1));

        md_system_free(&sys);
        md_arena_allocator_destroy(arena);
    }
    remove(path.ptr);
}

UTEST(amber_netcdf4, restart_and_what_is_not_amber) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_str_t path = MD_STR_LIT("md_unittest_amber4.ncrst");
    const nc4_case_t c = { "restart", .contiguous = true };
    ASSERT_TRUE(nc4_write(path, &c, true));
    EXPECT_EQ(MD_AMBER_NETCDF_RESTART, md_amber_netcdf_file_kind(path));
    EXPECT_EQ((size_t)NC4_ATOMS, md_amber_netcdf_file_atom_count(path));
    md_system_state_t s = { .alloc = arena };
    ASSERT_TRUE(md_amber_netcdf_state_init_from_file(&s, path, NC4_ATOMS));
    EXPECT_EQ((float)nc4_coord(0, 5, 2), s.xyz[5].z);
    EXPECT_NEAR(41.0f, s.unitcell.y, 1e-4f);

    // An HDF5 file without the convention is not Amber
    const hid_t file = H5Fcreate(path.ptr, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ASSERT_TRUE(file >= 0);
    nc4_text_attr(file, "Conventions", "CF-1.8");
    H5Fclose(file);
    EXPECT_EQ(MD_AMBER_NETCDF_NONE, md_amber_netcdf_file_kind(path));
    md_system_t sys = { .alloc = arena };
    EXPECT_FALSE(md_amber_netcdf_system_publish_run(&sys, path, MD_STR_LIT("run/a"), 0));
    remove(path.ptr);
    md_arena_allocator_destroy(arena);
}

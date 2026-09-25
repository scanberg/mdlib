#include "utest.h"

#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_xyz.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_lammps.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include <stdio.h>
#include <stddef.h>

#include "run_check.h"

// A '<file>.cache' is only as good as its match with the file it indexes. It carries the file's size
// and its modification time from the OS, as they were before the scan that made it, and is accepted
// only while both still match.

#define TEST_MAGIC   0x7e57ca4eull
#define TEST_VERSION 3

static bool copy_file(str_t dst, str_t src, md_allocator_i* alloc) {
    md_file_t in = {0}, out = {0};
    if (!md_file_open(&in, src, MD_FILE_READ)) return false;
    if (!md_file_open(&out, dst, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        md_file_close(&in);
        return false;
    }
    const size_t size = md_file_size(in);
    void* bytes = md_alloc(alloc, size);
    const bool ok = md_file_read(in, bytes, size) == size && md_file_write(out, bytes, size) == size;
    md_free(alloc, bytes, size);
    md_file_close(&in);
    md_file_close(&out);
    return ok;
}

static bool write_text(str_t path, const char* text) {
    md_file_t f = {0};
    if (!md_file_open(&f, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    const size_t len = strlen(text);
    const bool ok = md_file_write(f, text, len) == len;
    md_file_close(&f);
    return ok;
}

static bool read_header(md_run_cache_header_t* out, str_t cache_path) {
    md_file_t f = {0};
    if (!md_file_open(&f, cache_path, MD_FILE_READ)) return false;
    const bool ok = md_file_read(f, out, sizeof(*out)) == sizeof(*out);
    md_file_close(&f);
    return ok;
}

// Moves the stamp the cache was made with, as if the file had been changed since
static bool shift_stamp(str_t cache_path) {
    md_file_t f = {0};
    if (!md_file_open(&f, cache_path, MD_FILE_READ | MD_FILE_WRITE)) return false;
    md_run_cache_header_t h;
    bool ok = md_file_read(f, &h, sizeof(h)) == sizeof(h);
    h.source_modified += 1;
    ok = ok && md_file_write_at(f, offsetof(md_run_cache_header_t, source_modified), &h.source_modified, sizeof(h.source_modified)) == sizeof(h.source_modified);
    md_file_close(&f);
    return ok;
}

// A cache made with a stamp: the header, then one payload word
static bool make_cache(str_t src, const md_file_info_t* stamp, uint64_t magic, uint64_t version) {
    md_file_t f = {0};
    if (!md_run_cache_create(&f, src, stamp, magic, version, 7, 2)) return false;
    const uint64_t payload = 0x1234;
    const bool ok = md_file_write(f, &payload, sizeof(payload)) == sizeof(payload);
    md_file_close(&f);
    return ok;
}

static bool cache_opens(str_t src, uint64_t magic, uint64_t version) {
    md_file_t f = {0};
    md_run_cache_header_t h = {0};
    if (!md_run_cache_open(&f, &h, src, magic, version)) return false;
    md_file_close(&f);
    return true;
}

UTEST(run_cache, accepts_the_file_it_was_made_from) {
    const str_t src   = STR_INIT("md_unittest_run_cache.dat");
    const str_t cache = STR_INIT("md_unittest_run_cache.dat.cache");
    ASSERT_TRUE(write_text(src, "some frames"));
    md_file_info_t info = {0};
    ASSERT_TRUE(md_file_info_extract_from_path(src, &info));
    ASSERT_TRUE(make_cache(src, &info, TEST_MAGIC, TEST_VERSION));

    md_file_t f = {0};
    md_run_cache_header_t h = {0};
    ASSERT_TRUE(md_run_cache_open(&f, &h, src, TEST_MAGIC, TEST_VERSION));
    EXPECT_EQ((uint64_t)TEST_MAGIC, h.magic);
    EXPECT_EQ((uint64_t)TEST_VERSION, h.version);
    EXPECT_EQ((uint64_t)info.size, h.source_size);
    EXPECT_EQ(info.modified_time, h.source_modified);
    EXPECT_EQ(7u, h.num_atoms);
    EXPECT_EQ(2u, h.num_frames);
    // Positioned at what follows the header
    uint64_t payload = 0;
    EXPECT_EQ(sizeof(payload), md_file_read(f, &payload, sizeof(payload)));
    EXPECT_EQ(0x1234u, payload);
    md_file_close(&f);

    // Another format's, or another version's
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC + 1, TEST_VERSION));
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC, TEST_VERSION + 1));

    remove(cache.ptr);
    remove(src.ptr);
}

// The stamp is what the file was when the scan began, so a file written during the scan leaves a
// cache that is already stale.
UTEST(run_cache, rejects_another_size_or_time) {
    const str_t src   = STR_INIT("md_unittest_run_cache_stale.dat");
    const str_t cache = STR_INIT("md_unittest_run_cache_stale.dat.cache");
    ASSERT_TRUE(write_text(src, "some frames"));
    md_file_info_t info = {0};
    ASSERT_TRUE(md_file_info_extract_from_path(src, &info));

    md_file_info_t other = info;
    other.size += 1;
    ASSERT_TRUE(make_cache(src, &other, TEST_MAGIC, TEST_VERSION));
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC, TEST_VERSION));

    other = info;
    other.modified_time -= 1;
    ASSERT_TRUE(make_cache(src, &other, TEST_MAGIC, TEST_VERSION));
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC, TEST_VERSION));

    // The same size, written again: only the time tells
    ASSERT_TRUE(make_cache(src, &info, TEST_MAGIC, TEST_VERSION));
    EXPECT_TRUE(cache_opens(src, TEST_MAGIC, TEST_VERSION));
    ASSERT_TRUE(shift_stamp(cache));
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC, TEST_VERSION));

    // No file, no cache
    remove(src.ptr);
    EXPECT_FALSE(cache_opens(src, TEST_MAGIC, TEST_VERSION));
    remove(cache.ptr);
}

typedef bool (*publish_fn)(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena);

static bool publish_xtc(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena) {
    md_system_state_t st = {.alloc = arena};
    return md_gro_system_init_from_file(sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro")) && md_xtc_system_publish_run(sys, path, run, 0);
}
static bool publish_trr(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena) {
    md_system_state_t st = {.alloc = arena};
    return md_gro_system_init_from_file(sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro")) && md_trr_system_publish_run(sys, path, run, 0);
}
static bool publish_pdb(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena) {
    md_system_state_t st = {.alloc = arena};
    return md_pdb_system_init_from_file(sys, &st, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) && md_pdb_system_publish_run(sys, path, run, 0);
}
static bool publish_xyz(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena) {
    md_system_state_t st = {.alloc = arena};
    return md_xyz_system_init_from_file(sys, &st, path, MD_XYZ_OPTION_DISABLE_CACHE_WRITE) && md_xyz_system_publish_run(sys, path, run, 0);
}
static bool publish_lammps(md_system_t* sys, str_t path, str_t run, md_allocator_i* arena) {
    (void)arena;
    return md_lammps_system_publish_run(sys, path, run, 0);
}

// Each format end to end: the first publish writes the cache stamped with the file, the second
// reads it and leaves it be, and once the stamp no longer matches the file the next publish scans
// again and writes a cache that does.
static void check_format(int* utest_result, str_t data, str_t copy, publish_fn publish, size_t expected_frames) {
    const str_t run = STR_INIT("run/c");
    char cache_buf[512];
    const str_t cache = {cache_buf, (size_t)snprintf(cache_buf, sizeof(cache_buf), STR_FMT ".cache", STR_ARG(copy))};

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    ASSERT_TRUE(copy_file(copy, data, arena));
    remove(cache.ptr);
    md_file_info_t info = {0};
    ASSERT_TRUE(md_file_info_extract_from_path(copy, &info));

    for (int pass = 0; pass < 3; ++pass) {
        md_file_info_t before = {0};
        const bool had_cache = pass == 1 && md_file_info_extract_from_path(cache, &before);

        md_system_t sys = {.alloc = arena};
        sys.attributes.alloc = arena;
        EXPECT_TRUE(publish(&sys, copy, run, arena));
        EXPECT_EQ(expected_frames, run_num_frames(&sys, run));
        md_system_free(&sys);

        md_run_cache_header_t h = {0};
        ASSERT_TRUE(read_header(&h, cache));
        EXPECT_EQ((uint64_t)info.size, h.source_size);
        EXPECT_EQ(info.modified_time, h.source_modified);
        EXPECT_EQ((uint64_t)expected_frames, h.num_frames);

        if (pass == 1) {
            // A current cache is read, not written again
            md_file_info_t after = {0};
            ASSERT_TRUE(had_cache);
            ASSERT_TRUE(md_file_info_extract_from_path(cache, &after));
            EXPECT_EQ(before.modified_time, after.modified_time);
            ASSERT_TRUE(shift_stamp(cache));
        }
    }
    remove(cache.ptr);
    remove(copy.ptr);
    md_vm_arena_destroy(arena);
}

UTEST(run_cache, xtc) {
    check_format(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc"), STR_LIT("md_unittest_run_cache.xtc"), publish_xtc, 501);
}
UTEST(run_cache, trr) {
    check_format(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), STR_LIT("md_unittest_run_cache.trr"), publish_trr, 101);
}
UTEST(run_cache, pdb) {
    check_format(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), STR_LIT("md_unittest_run_cache.pdb"), publish_pdb, 38);
}
UTEST(run_cache, xyz) {
    check_format(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/traj-30-P_10.xyz"), STR_LIT("md_unittest_run_cache.xyz"), publish_xyz, 10);
}
UTEST(run_cache, lammps) {
    check_format(utest_result, STR_LIT(MD_UNITTEST_DATA_DIR "/cubic_standardASCII.lammpstrj"), STR_LIT("md_unittest_run_cache.lammpstrj"), publish_lammps, 10);
}

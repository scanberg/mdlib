#include "utest.h"

#include <stdio.h>

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_platform.h>

UTEST(os, path_canonical) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    {
        char buf[1024];
        md_path_write_cwd(buf, sizeof(buf));
        md_str_t path = md_str_printf(temp_arena, "%s/../../", buf);
        md_str_t result = md_path_make_canonical(path, temp_arena);
        printf("result: '%.*s'\n", (int)result.len, result.ptr);
    }
    {
        /*
        // FAILS ON UNIX
        // Probably because the file does not exist
        md_str_t path = MD_STR_LIT("cool/fool/../bool/../file.txt");
        md_str_t result = md_path_make_canonical(path, temp_arena);
        md_str_t ref = MD_STR_LIT("cool/file.txt");
        EXPECT_TRUE(str_equal(result, ref));
        */
    }
    {
        md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t result = md_path_make_canonical(path, temp_arena);
        printf("canonical: '%.*s'\n", (int)result.len, result.ptr);
    }
    {
        // '.' and '..' are resolved and the path does not have to exist
        md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/./../dir/subdir/does_not_exist.txt");
        md_str_t result = md_path_make_canonical(path, temp_arena);
        EXPECT_TRUE(md_str_ends_with(result, MD_STR_LIT("/dir/subdir/does_not_exist.txt")));
    }
    {
        // A directory is reported with a trailing separator
        md_str_t result = md_path_make_canonical(MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir"), temp_arena);
        EXPECT_TRUE(md_str_ends_with(result, MD_STR_LIT("/dir/subdir/")));
    }
    md_temp_end(temp);
}

UTEST(os, path_relative) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("../../40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./dir/subdir/", result.ptr);
    }
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/" );
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t to = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./file.txt", result.ptr);
    }
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        md_str_t to = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./", result.ptr);
    }
    {
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        md_str_t to = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.dat");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./file.dat", result.ptr);
    }
    {
        // Whole components are compared: 'sub' is not a prefix folder of 'subdir'
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/sub/file.txt");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("../sub/file.txt", result.ptr);
    }
    {
        // Neither path has to exist: this is the case when saving a new workspace
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/does_not_exist.vwsp");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        md_str_t result = md_path_make_relative(from, to, temp_arena);
        EXPECT_STREQ("./40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        // folder(from) + relative(from, to) has to resolve back to to
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        md_str_t rel  = md_path_make_relative(from, to, temp_arena);
        md_str_t folder = {0};
        EXPECT_TRUE(md_extract_folder_path(&folder, from));
        md_str_t joined = md_str_printf(temp_arena, MD_STR_FMT MD_STR_FMT, MD_STR_ARG(folder), MD_STR_ARG(rel));
        EXPECT_TRUE(md_str_eq(md_path_make_canonical(joined, temp_arena), md_path_make_canonical(to, temp_arena)));
    }
    {
        // An empty path yields an empty result, not whatever happened to be on the stack
        md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        EXPECT_TRUE(md_str_empty(md_path_make_relative(path, (md_str_t){0}, temp_arena)));
        EXPECT_TRUE(md_str_empty(md_path_make_relative((md_str_t){0}, path, temp_arena)));
    }
    {
        // A buffer which cannot hold the result fails cleanly
        char buf[8] = {'x'};
        md_str_t from = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        md_str_t to   = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        EXPECT_EQ(0u, md_path_write_relative(buf, sizeof(buf), from, to));
        EXPECT_EQ('\0', buf[0]);
    }
#if MD_PLATFORM_WINDOWS
    {
        // Separate volumes and UNC shares have no relative path between them
        EXPECT_TRUE(md_str_empty(md_path_make_relative(MD_STR_LIT("C:/some/dir/file.txt"), MD_STR_LIT("D:/other/dir/file.txt"), temp_arena)));
        EXPECT_TRUE(md_str_empty(md_path_make_relative(MD_STR_LIT("//server_a/share/file.txt"), MD_STR_LIT("//server_b/share/file.txt"), temp_arena)));
    }
#endif
    md_temp_end(temp);
}

UTEST(os, sys_info) {
	md_os_sys_info_t info = { 0 };
	EXPECT_TRUE(md_os_sys_info_query(&info));
    printf("total physical ram %i GB\n", (int)(info.physical_ram_bytes / MD_GIGABYTES(1)));
	printf("physical cores: %d\n", info.num_physical_cores);
	printf("virtual cores: %d\n", info.num_virtual_cores);
    EXPECT_GT(info.physical_ram_bytes, 0);
	EXPECT_GT(info.num_physical_cores, 0);
    EXPECT_GT(info.num_virtual_cores,  0);
}

UTEST(os, mem) {
    const uint64_t reserve_size = MD_GIGABYTES(8);
    void* ptr = md_vm_reserve(reserve_size);
    EXPECT_NE(ptr, NULL);

    md_vm_commit(ptr, MD_MEGABYTES(1));
    md_vm_commit((char*)ptr + MD_MEGABYTES(1), MD_MEGABYTES(2));

    md_vm_decommit((char*)ptr + MD_MEGABYTES(1), MD_MEGABYTES(2));
    md_vm_decommit(ptr, MD_MEGABYTES(1));

    const uint64_t commit_size = MD_GIGABYTES(1);
    md_vm_commit(ptr, commit_size);
    uint64_t* c = ptr;
    for (uint64_t i = 0; i < commit_size / sizeof(uint64_t); ++i) {
        c[i] = i;
    }

    md_vm_release(ptr, reserve_size);
}

static void other_func(int* value) {
    *value = 10;
}

static void function(void* user_data) {
    int* value = (int*)user_data;
    other_func(value);
    *value = 5;
}

UTEST(os, thread) {
    int number_of_the_beast = 666;
    md_thread_t* thread = md_thread_create(function, &number_of_the_beast);
    EXPECT_TRUE(md_thread_join(thread));
    EXPECT_EQ(5, number_of_the_beast);
}

UTEST(os, mutex) {
    md_mutex_t mutex;

    EXPECT_TRUE(md_mutex_init(&mutex));

    md_mutex_lock(&mutex);
    md_mutex_unlock(&mutex);

    EXPECT_TRUE(md_mutex_try_lock(&mutex));
    md_mutex_unlock(&mutex);

    EXPECT_TRUE(md_mutex_destroy(&mutex));
}

UTEST(os, semaphore) {
    md_semaphore_t sema = {0};

    EXPECT_TRUE(md_semaphore_init(&sema, 8));

    md_semaphore_aquire(&sema);
    md_semaphore_aquire(&sema);
    md_semaphore_aquire(&sema);
    md_semaphore_aquire(&sema);

    EXPECT_TRUE(md_semaphore_try_aquire(&sema));
    EXPECT_TRUE(md_semaphore_try_aquire(&sema));
    EXPECT_TRUE(md_semaphore_try_aquire(&sema));
    EXPECT_TRUE(md_semaphore_try_aquire(&sema));

    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));

    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));
    EXPECT_TRUE(md_semaphore_release(&sema));

    EXPECT_TRUE(md_semaphore_destroy(&sema));
}

// The three converters are one scale expressed three ways, so they are checked against each
// other and against a known wall-clock duration. The ratios alone would not catch a factor
// applied to all three, and the duration alone would not catch one of them drifting from the
// others - which is exactly how md_tick_to_nanoseconds returned microseconds on windows.
UTEST(os, time_converters_agree) {
    const md_tick_t t0 = md_tick_now();
    md_thread_sleep(50);
    const md_tick_t dt = md_tick_now() - t0;
    ASSERT_GT(dt, 0);

    const double ns = md_tick_to_nanoseconds(dt);
    const double ms = md_tick_to_milliseconds(dt);
    const double s  = md_tick_to_seconds(dt);

    // Same quantity, different units. 1e-9 relative slack is far tighter than a unit slip.
    EXPECT_NEAR(ns, ms * 1.0e6, ns * 1.0e-9);
    EXPECT_NEAR(ms, s  * 1.0e3, ms * 1.0e-9);

    // And the scale itself is right. The band is wide because a sleep only guarantees a lower
    // bound and a loaded machine can overshoot badly; it is still a thousand times tighter
    // than the error it exists to catch.
    EXPECT_GE(ms, 25.0);
    EXPECT_LT(ms, 5000.0);
}

UTEST(os, path_validity) {
    const md_str_t file_path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    const md_str_t dir_path  = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
    const md_str_t bad_path  = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/missing.txt");

    EXPECT_TRUE(md_path_is_valid(file_path));
    EXPECT_TRUE(md_path_is_valid(dir_path));
    EXPECT_FALSE(md_path_is_valid(bad_path));

    EXPECT_FALSE(md_path_is_directory(file_path));
    EXPECT_TRUE(md_path_is_directory(dir_path));
    EXPECT_FALSE(md_path_is_directory(bad_path));

    // An empty path is never valid and never a directory (and must not read past the str)
    EXPECT_FALSE(md_path_is_valid((md_str_t){0}));
    EXPECT_FALSE(md_path_is_directory((md_str_t){0}));

    // Absolute paths are recognized the same way on every platform
    EXPECT_TRUE(md_path_is_absolute(MD_STR_LIT("/home/user/file.txt")));
    EXPECT_TRUE(md_path_is_absolute(MD_STR_LIT("C:/Users/user/file.txt")));
    EXPECT_TRUE(md_path_is_absolute(MD_STR_LIT("C:\\Users\\user\\file.txt")));
    EXPECT_TRUE(md_path_is_absolute(MD_STR_LIT("//server/share/file.txt")));
    EXPECT_FALSE(md_path_is_absolute(MD_STR_LIT("./file.txt")));
    EXPECT_FALSE(md_path_is_absolute(MD_STR_LIT("../dir/file.txt")));
    EXPECT_FALSE(md_path_is_absolute(MD_STR_LIT("dir/file.txt")));
    EXPECT_FALSE(md_path_is_absolute((md_str_t){0}));
}

UTEST(os, file_info_extract) {
    const md_str_t file_path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    const md_str_t dir_path  = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, file_path, MD_FILE_READ));

    md_file_info_t info = {0};
    ASSERT_TRUE(md_file_info_extract(file, &info));
    EXPECT_EQ(info.size, 5);
    EXPECT_TRUE((info.flags & MD_FILE_INFO_IS_REGULAR) != 0);
    EXPECT_FALSE((info.flags & MD_FILE_INFO_IS_DIRECTORY) != 0);

    md_file_info_t path_info = {0};
    ASSERT_TRUE(md_file_info_extract_from_path(file_path, &path_info));
    EXPECT_EQ(path_info.size, 5);
    EXPECT_TRUE((path_info.flags & MD_FILE_INFO_IS_REGULAR) != 0);

    md_file_info_t dir_info = {0};
    ASSERT_TRUE(md_file_info_extract_from_path(dir_path, &dir_info));
    EXPECT_TRUE((dir_info.flags & MD_FILE_INFO_IS_DIRECTORY) != 0);

    EXPECT_TRUE(md_file_close(&file));
    EXPECT_FALSE(md_file_valid(file));
}

UTEST(os, file_read_seek_tell_size) {
    const md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    EXPECT_EQ(md_file_size(file), 5);
    EXPECT_EQ(md_file_tell(file), 0);

    char buf[6] = {0};
    EXPECT_EQ(md_file_read(file, buf, 2), 2);
    EXPECT_EQ(md_file_tell(file), 2);
    EXPECT_STREQ(buf, "he");

    EXPECT_TRUE(md_file_seek(file, 0, MD_FILE_BEG));
    md_memset(buf, 0, sizeof(buf));
    EXPECT_EQ(md_file_read(file, buf, 5), 5);
    EXPECT_STREQ(buf, "hello");
    EXPECT_EQ(md_file_tell(file), 5);

    EXPECT_TRUE(md_file_seek(file, -2, MD_FILE_END));
    EXPECT_EQ(md_file_tell(file), 3);
    md_memset(buf, 0, sizeof(buf));
    EXPECT_EQ(md_file_read(file, buf, 2), 2);
    EXPECT_STREQ(buf, "lo");

    EXPECT_TRUE(md_file_close(&file));
    EXPECT_FALSE(md_file_valid(file));
}

UTEST(os, env_set_get_remove) {
    const char* name = "MD_OS_TEST_ENV_VARIABLE";
    ASSERT_TRUE(md_os_env_set(name, NULL));     // removing what is not there is no failure
    EXPECT_TRUE(md_os_env_get(name) == NULL);
    ASSERT_TRUE(md_os_env_set(name, "a value"));
    ASSERT_TRUE(md_os_env_get(name) != NULL);
    EXPECT_STREQ(md_os_env_get(name), "a value");
    ASSERT_TRUE(md_os_env_set(name, "another"));
    EXPECT_STREQ(md_os_env_get(name), "another");
    ASSERT_TRUE(md_os_env_set(name, NULL));
    EXPECT_TRUE(md_os_env_get(name) == NULL);
    EXPECT_FALSE(md_os_env_set("A=B", "c"));    // not a name
    EXPECT_FALSE(md_os_env_set("", "c"));
}

// On every platform: on Windows the file pointer which a positional read moves is put back
UTEST(os, file_read_at_preserves_offset) {
    const md_str_t path = MD_STR_INIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    EXPECT_TRUE(md_file_seek(file, 2, MD_FILE_BEG));
    EXPECT_EQ(md_file_tell(file), 2);

    char whole[6] = {0};
    EXPECT_EQ(md_file_read_at(file, 0, whole, 5), 5);
    EXPECT_STREQ(whole, "hello");
    EXPECT_EQ(md_file_tell(file), 2);

    char tail[4] = {0};
    EXPECT_EQ(md_file_read(file, tail, 3), 3);
    EXPECT_STREQ(tail, "llo");

    EXPECT_TRUE(md_file_close(&file));
}

// A file of n bytes, byte i of which is a function of i: what a read at any offset must give
static uint8_t file_pattern(size_t i) {
    return (uint8_t)((i * 2654435761u) >> 13);
}

static bool write_pattern_file(const char* path, size_t n) {
    md_file_t file = {0};
    if (!md_file_open(&file, md_str_from_cstr(path), MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;
    uint8_t buf[4096];
    bool ok = true;
    for (size_t at = 0; at < n && ok; at += sizeof(buf)) {
        const size_t len = MD_MIN(sizeof(buf), n - at);
        for (size_t i = 0; i < len; ++i) buf[i] = file_pattern(at + i);
        ok = md_file_write(file, buf, len) == len;
    }
    md_file_close(&file);
    return ok;
}

typedef struct concurrent_reads_t {
    md_file_t file;
    size_t    size;
    uint32_t  seed;
    size_t    wrong;    // bytes which differ from the pattern, or reads short of what was asked
} concurrent_reads_t;

static void concurrent_reads(void* arg) {
    concurrent_reads_t* r = (concurrent_reads_t*)arg;
    uint8_t buf[70000];
    uint32_t x = r->seed;
    for (int k = 0; k < 400; ++k) {
        x = x * 1664525u + 1013904223u;
        const size_t len = 1 + (x >> 8) % sizeof(buf);
        x = x * 1664525u + 1013904223u;
        const size_t at = (size_t)x % (r->size - len);
        if (md_file_read_at(r->file, (md_file_offset_t)at, buf, len) != len) {
            r->wrong += 1;
            continue;
        }
        for (size_t i = 0; i < len; ++i) r->wrong += buf[i] != file_pattern(at + i);
    }
}

// MD_FILE_CONCURRENT: one handle read at offsets by several threads at once, each read what is at its offset; reads past
// the end come back short; and the calls which take a position fail
UTEST(os, file_concurrent_reads) {
    const char* path = "md_unittest_concurrent.bin";
    const size_t size = MD_MEGABYTES(4) + 123;
    ASSERT_TRUE(write_pattern_file(path, size));

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_READ | MD_FILE_CONCURRENT | MD_FILE_SEQUENTIAL));
    EXPECT_EQ((md_file_offset_t)size, md_file_size(file));

    concurrent_reads_t reads[8];
    md_thread_t* threads[8];
    for (int t = 0; t < 8; ++t) {
        reads[t] = (concurrent_reads_t){ .file = file, .size = size, .seed = 17u + (uint32_t)t };
        threads[t] = md_thread_create(concurrent_reads, &reads[t]);
        ASSERT_TRUE(threads[t] != NULL);
    }
    size_t wrong = 0;
    for (int t = 0; t < 8; ++t) {
        md_thread_join(threads[t]);
        wrong += reads[t].wrong;
    }
    EXPECT_EQ((size_t)0, wrong);

    // The end
    uint8_t tail[64];
    EXPECT_EQ((size_t)23, md_file_read_at(file, (md_file_offset_t)(size - 23), tail, sizeof(tail)));
    EXPECT_EQ(file_pattern(size - 1), tail[22]);
    EXPECT_EQ((size_t)0, md_file_read_at(file, (md_file_offset_t)size + 10, tail, sizeof(tail)));

    // No position
    EXPECT_EQ((size_t)0, md_file_read(file, tail, 8));
    EXPECT_FALSE(md_file_seek(file, 0, MD_FILE_BEG));
    EXPECT_EQ((md_file_offset_t)-1, md_file_tell(file));
    md_file_close(&file);

    // Writes at offsets on a concurrent handle, read back through an ordinary one
    ASSERT_TRUE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_READ | MD_FILE_WRITE | MD_FILE_CONCURRENT));
    const char text[] = "concurrent";
    EXPECT_EQ(sizeof(text), md_file_write_at(file, 1000, text, sizeof(text)));
    EXPECT_EQ((size_t)0, md_file_write(file, text, sizeof(text)));
    md_file_close(&file);
    ASSERT_TRUE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_READ));
    char back[sizeof(text)] = {0};
    EXPECT_EQ(sizeof(text), md_file_read_at(file, 1000, back, sizeof(back)));
    EXPECT_STREQ(text, back);
    md_file_close(&file);

    // What does not go together
    EXPECT_FALSE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_APPEND | MD_FILE_CONCURRENT));
    EXPECT_FALSE(md_file_open(&file, md_str_from_cstr(path), MD_FILE_READ | MD_FILE_SEQUENTIAL | MD_FILE_RANDOM));
    remove(path);
}

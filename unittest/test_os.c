#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

UTEST(os, path_canonical) {
    {
        char buf[1024];
        md_path_write_cwd(buf, sizeof(buf));
        str_t path = str_printf(md_get_temp_allocator(), "%s/../../", buf);
        str_t result = md_path_make_canonical(path, md_get_temp_allocator());
        printf("result: '%.*s'\n", (int)result.len, result.ptr);
    }
    {
        /*
        // FAILS ON UNIX
        // Probably because the file does not exist
        str_t path = STR_LIT("cool/fool/../bool/../file.txt");
        str_t result = md_path_make_canonical(path, md_get_temp_allocator());
        str_t ref = STR_LIT("cool/file.txt");
        EXPECT_TRUE(str_equal(result, ref));
        */
    }
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_canonical(path, md_get_temp_allocator());
        printf("canonical: '%.*s'\n", (int)result.len, result.ptr);
    }
}

UTEST(os, path_relative) {
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("../../40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("./dir/subdir/", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/" );
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("./40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("./file.txt", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("./", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.dat");
        str_t result = md_path_make_relative(from, to, md_get_temp_allocator());
        EXPECT_STREQ("./file.dat", result.ptr);
    }
}

UTEST(os, sys_info) {
	md_os_sys_info_t info = { 0 };
	EXPECT_TRUE(md_os_sys_info_query(&info));
    printf("total physical ram %i GB\n", (int)(info.physical_ram_bytes / GIGABYTES(1)));
	printf("physical cores: %d\n", info.num_physical_cores);
	printf("virtual cores: %d\n", info.num_virtual_cores);
    EXPECT_GT(info.physical_ram_bytes, 0);
	EXPECT_GT(info.num_physical_cores, 0);
    EXPECT_GT(info.num_virtual_cores,  0);
}

UTEST(os, mem) {
    const uint64_t reserve_size = GIGABYTES(8);
    void* ptr = md_vm_reserve(reserve_size);
    EXPECT_NE(ptr, NULL);

    md_vm_commit(ptr, MEGABYTES(1));
    md_vm_commit((char*)ptr + MEGABYTES(1), MEGABYTES(2));

    md_vm_decommit((char*)ptr + MEGABYTES(1), MEGABYTES(2));
    md_vm_decommit(ptr, MEGABYTES(1));

    const uint64_t commit_size = GIGABYTES(1);
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

UTEST(os, path_validity) {
    const str_t file_path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    const str_t dir_path  = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
    const str_t bad_path  = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/missing.txt");

    EXPECT_TRUE(md_path_is_valid(file_path));
    EXPECT_TRUE(md_path_is_valid(dir_path));
    EXPECT_FALSE(md_path_is_valid(bad_path));

    EXPECT_FALSE(md_path_is_directory(file_path));
    EXPECT_TRUE(md_path_is_directory(dir_path));
    EXPECT_FALSE(md_path_is_directory(bad_path));
}

UTEST(os, file_info_extract) {
    const str_t file_path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    const str_t dir_path  = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");

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
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, path, MD_FILE_READ));

    EXPECT_EQ(md_file_size(file), 5);
    EXPECT_EQ(md_file_tell(file), 0);

    char buf[6] = {0};
    EXPECT_EQ(md_file_read(file, buf, 2), 2);
    EXPECT_EQ(md_file_tell(file), 2);
    EXPECT_STREQ(buf, "he");

    EXPECT_TRUE(md_file_seek(file, 0, MD_FILE_BEG));
    MEMSET(buf, 0, sizeof(buf));
    EXPECT_EQ(md_file_read(file, buf, 5), 5);
    EXPECT_STREQ(buf, "hello");
    EXPECT_EQ(md_file_tell(file), 5);

    EXPECT_TRUE(md_file_seek(file, -2, MD_FILE_END));
    EXPECT_EQ(md_file_tell(file), 3);
    MEMSET(buf, 0, sizeof(buf));
    EXPECT_EQ(md_file_read(file, buf, 2), 2);
    EXPECT_STREQ(buf, "lo");

    EXPECT_TRUE(md_file_close(&file));
    EXPECT_FALSE(md_file_valid(file));
}

#if 0
// This cannot be guaranteed on windows because the underlying implementation uses the synchronous path.
UTEST(os, file_read_at_preserves_offset) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
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
#endif

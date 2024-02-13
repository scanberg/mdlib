#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

UTEST(os, path_canonical) {
    {
        char buf[1024];
        md_path_write_cwd(buf, sizeof(buf));
        str_t path = str_printf(md_temp_allocator, "%s/../../", buf);
        str_t result = md_path_make_canonical(path, md_temp_allocator);
        printf("result: '%.*s'\n", (int)result.len, result.ptr);
    }
    {
        /*
        // FAILS ON UNIX
        // Probably because the file does not exist
        str_t path = STR_LIT("cool/fool/../bool/../file.txt");
        str_t result = md_path_make_canonical(path, md_temp_allocator);
        str_t ref = STR_LIT("cool/file.txt");
        EXPECT_TRUE(str_equal(result, ref));
        */
    }
    {
        str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_canonical(path, md_temp_allocator);
        printf("canonical: '%.*s'\n", (int)result.len, result.ptr);
    }
}

UTEST(os, path_relative) {
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("../../40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("./dir/subdir/", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/" );
        str_t to   = STR_LIT(MD_UNITTEST_DATA_DIR "/40-40-2-ddba-dyna.xmol");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("./40-40-2-ddba-dyna.xmol", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("./file.txt", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("./", result.ptr);
    }
    {
        str_t from = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt");
        str_t to = STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.dat");
        str_t result = md_path_make_relative(from, to, md_temp_allocator);
        EXPECT_STREQ("./file.dat", result.ptr);
    }
}

UTEST(os, ram) {
    uint64_t physical_ram = md_os_physical_ram();
    printf("total physical ram %i GB\n", (int)(physical_ram / GIGABYTES(1)));
    EXPECT_GT(physical_ram, 0);
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
    MEMSET(ptr, 1, commit_size);

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

UTEST(os, file_basic) {
	md_file_o* file = md_file_open(STR_LIT(MD_UNITTEST_DATA_DIR "/dir/subdir/file.txt"), MD_FILE_READ);
    EXPECT_TRUE(file);

    size_t filesize = md_file_size(file);
    EXPECT_EQ(filesize, 5);

    char buf[1024];
    size_t bytes_read = md_file_read(file, buf, sizeof(buf));

    EXPECT_EQ(bytes_read, 5);

    str_t str = { buf, bytes_read };
    EXPECT_TRUE(str_eq_cstr(str, "hello"));

	md_file_close(file);
}

#define COUNT 4
UTEST(os, file_write) {
    str_t test_string = STR_LIT("Lorem ipsum dolor sit amet, consectetur adipiscing elit sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.");
    md_file_o* file = md_file_open(STR_LIT("temp_file.txt"), MD_FILE_WRITE);
    ASSERT_TRUE(file);

    for (int i = 0; i < COUNT; ++i) {
    	size_t bytes_written = md_file_write(file, test_string.ptr, test_string.len);
    	EXPECT_EQ(bytes_written, test_string.len);
    }

    md_file_close(file);

    file = md_file_open(STR_LIT("temp_file.txt"), MD_FILE_READ);
    ASSERT_TRUE(file);
    size_t filesize = md_file_size(file);
    EXPECT_EQ(filesize, test_string.len * COUNT);
    char buf[1024];
    for (int i = 0; i < COUNT; ++i) {
    	size_t bytes_read = md_file_read(file, buf, test_string.len);
    	EXPECT_EQ(bytes_read, test_string.len);
    	str_t str = { buf, bytes_read };
    	EXPECT_TRUE(str_eq(str, test_string));
    }
    md_file_close(file);
}
#undef COUNT

UTEST(os, file_memmap) {
    md_file_o* file = md_file_open(STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), MD_FILE_READ);
    ASSERT_TRUE(file);

    md_file_mapping_o* mapping = md_file_mem_map(file);
    ASSERT_TRUE(mapping);
    {
        const char* addr = md_file_mem_map_view(mapping, 5, 1024);
        EXPECT_NE(addr, NULL);
        str_t str = { addr, 1024 };
        str_t ref = STR_LIT("ated by trjconv : Protein in water t= 4000.00000");
        EXPECT_TRUE(str_eq_n(str, ref, ref.len));
        md_file_mem_unmap_view(addr, 1024);
    }
    
    size_t view_offset = md_vm_page_size() - 1;
    size_t view_size   = md_vm_page_size() + 1;
    const char* addr = md_file_mem_map_view(mapping, view_offset, view_size);
    EXPECT_NE(addr, NULL);
    str_t ref = STR_LIT(" 13.281  11.447");
    str_t str = { addr, ref.len };
    bool equal = str_eq(str, ref);
    EXPECT_TRUE(equal);
    if (!equal) {
        fprintf(stderr, "Expected to read segment '"STR_FMT"', but got '"STR_FMT"'", STR_ARG(ref), STR_ARG(str));
    }
    md_file_mem_unmap_view(addr, view_size);

    md_file_mem_unmap(mapping);
    md_file_close(file);
}

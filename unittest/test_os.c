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

#include <core/md_str_builder.h>

UTEST(os, file_read_lines) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro");
    str_t ref = load_textfile(path, md_heap_allocator);

    
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    
    md_strb_t sb = md_strb_create(md_heap_allocator);
    
    const int64_t cap = KILOBYTES(1);
    char* buf = md_alloc(md_heap_allocator, cap);
    
    int64_t bytes_read;
    while (bytes_read = md_file_read_lines(file, buf, cap)) {
        md_strb_push_str(&sb, (str_t){buf, bytes_read});
    }

    str_t str = md_strb_to_str(sb);
    EXPECT_TRUE(str_eq(str, ref));

    md_free(md_heap_allocator, buf, cap);
}
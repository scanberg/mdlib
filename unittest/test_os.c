#include "utest.h"

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <core/md_os.h>

UTEST(os, path) {
    {
        char buf[4096];
        str_t cwd = md_path_cwd();
        snprintf(buf, sizeof(buf), "%.*s/../../", (int)cwd.len, cwd.ptr);
        str_t path = {buf, strnlen(buf, sizeof(buf))};
        str_t result = md_path_make_canonical(path, default_temp_allocator);
        printf("result: '%.*s'\n", (int)result.len, result.ptr);
    }
    {
        /*
        // FAILS ON UNIX
        str_t path = STR("cool/fool/../bool/../file.txt");
        str_t result = md_path_make_canonical(path, default_temp_allocator);
        str_t ref = STR("cool/file.txt");
        EXPECT_TRUE(str_equal(result, ref));
        */
    }
}

UTEST(os, ram) {
    uint64_t physical_ram = md_physical_ram();
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
    str_t path = STR(MD_UNITTEST_DATA_DIR "/pftaa.gro");
    str_t ref = load_textfile(path, default_allocator);

    
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    
    md_strb_t sb = md_strb_create(default_allocator);
    
    const int64_t cap = KILOBYTES(1);
    char* buf = md_alloc(default_allocator, cap);
    
    int64_t bytes_read;
    while (bytes_read = md_file_read_lines(file, buf, cap)) {
        md_strb_push_str(&sb, (str_t){buf, bytes_read});
    }

    str_t str = md_strb_to_str(&sb);
    EXPECT_TRUE(str_equal(str, ref));

    md_free(default_allocator, buf, cap);
}
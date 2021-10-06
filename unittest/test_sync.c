#include "utest.h"

#include "core/md_sync.h"

#include <stdio.h>

void other_func(int* value) {
    printf("Printing yo!");
    *value = 10;
}

void function(void* user_data) {
    int* value = (int*)user_data;
    other_func(value);
    *value = 5;
}

UTEST(sync, thread) {
    int value = 666;
    md_thread_t* thread = md_thread_create(function, &value);
    EXPECT_TRUE(md_thread_join(thread));
    EXPECT_EQ(5, value);
}

// @NOTE: Pool is an outlier here, since it is meant for allocations of a fixed size, thus cannot be tested with the common allocator test

UTEST(sync, mutex) {
    md_mutex_t mutex;

    EXPECT_TRUE(md_mutex_init(&mutex));

    md_mutex_lock(&mutex);
    md_mutex_unlock(&mutex);

    EXPECT_TRUE(md_mutex_try_lock(&mutex));
    md_mutex_unlock(&mutex);

    EXPECT_TRUE(md_mutex_destroy(&mutex));
}

UTEST(sync, semaphore) {
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

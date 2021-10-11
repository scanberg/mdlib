#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef struct md_semaphore_t md_semaphore_t;
typedef struct md_mutex_t md_mutex_t;
typedef struct md_thread_t md_thread_t;
typedef uint64_t md_thread_id_t;
typedef void (md_thread_func)(void *user_data);

#ifdef __cplusplus
extern "C" {
#endif

struct md_mutex_t {
    union {
        void* _align;
        char _data[64];
    };
};

struct md_semaphore_t {
    void* _data[4];
};

// Thread
md_thread_t* md_thread_create(md_thread_func func, void* user_data);

void md_thread_detach(md_thread_t* thread);
bool md_thread_join(md_thread_t* thread);

md_thread_id_t md_thread_get_id(md_thread_t* thread);
md_thread_id_t md_thread_id(void);

// Mutex
md_mutex_t md_mutex_create(void);
bool md_mutex_init(md_mutex_t* mutex);
bool md_mutex_destroy(md_mutex_t* mutex);

bool md_mutex_lock(md_mutex_t* mutex);

bool md_mutex_try_lock(md_mutex_t* mutex);
bool md_mutex_unlock(md_mutex_t* mutex);

// Semaphore
md_semaphore_t md_semaphore_create(int32_t initial_count);
bool md_semaphore_init(md_semaphore_t* semaphore, int32_t initial_count);
bool md_semaphore_destroy(md_semaphore_t* semaphore);

bool md_semaphore_aquire(md_semaphore_t* semaphore);
bool md_semaphore_try_aquire(md_semaphore_t* semaphore);
bool md_semaphore_try_aquire_n(md_semaphore_t* semaphore, int32_t count);

bool md_semaphore_query_count(md_semaphore_t* semaphore, int32_t* count);

bool md_semaphore_release(md_semaphore_t* semaphore);
bool md_semaphore_release_n(md_semaphore_t* semaphore, int32_t count);

#ifdef __cplusplus
}
#endif

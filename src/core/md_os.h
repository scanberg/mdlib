#pragma once

#include <stdint.h>
#include "md_str.h"

// Collection of OS specific things

struct  md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// ### PATH ###
// Gets the current working directory
// Returns a non owning reference to a static string buffer
str_t       md_os_path_cwd(void);

str_t       md_os_path_make_canonical(str_t path, struct md_allocator_i* alloc);
str_t       md_os_path_make_relative(str_t path_from, str_t path_to, struct md_allocator_i* alloc);

// Checks if a path to a file or directory is valid: i.e. it points to an actual file / dir on disk.
bool		md_os_path_is_valid(str_t path);
bool		md_os_path_is_directory(str_t path);

// ### TIME ###
// This represents a os specific time stamp with the highest precision available (usually nanoseconds)
// It is possible to subtract or perform simple arithmetic directly on the timestamp then
// Extract seconds or milliseconds from it.
typedef int64_t md_timestamp_t;

md_timestamp_t md_os_time_current(void);

double      md_os_time_as_nanoseconds(md_timestamp_t t);
double      md_os_time_as_milliseconds(md_timestamp_t t);
double      md_os_time_as_seconds(md_timestamp_t t);

// ### MEMORY ###
// This exposes system calls to the OS virtual memory allocator.
// If you are using this, the assumption is that you know what you are doing
// pointer and size should be aligned on page_size, otherwise it will most certainly crash on some OS

// returns total physical ram available on the machine
uint64_t    md_physical_ram();

// Virtual Memory
uint64_t	md_vm_page_size(void);
void*		md_vm_reserve(uint64_t size);
void		md_vm_release(void* ptr, uint64_t size);
void		md_vm_commit(void* ptr, uint64_t size);
void		md_vm_decommit(void* ptr, uint64_t size);

// ### THREAD ###
typedef void (*md_thread_exit) (void* data);
typedef void (*md_thread_entry)(void *data);
typedef struct md_thread_t md_thread_t;
typedef uint64_t md_thread_id_t;

md_thread_t*	md_thread_create(md_thread_entry func, void* data);
void			md_thread_detach(md_thread_t* thread);
bool			md_thread_join(md_thread_t* thread);

// Register a callback 
bool			md_thread_on_exit(md_thread_exit callback);

// Get ID from supplied thread object
md_thread_id_t	md_thread_get_id(md_thread_t* thread);

// Get ID of current thread
md_thread_id_t	md_thread_id(void);

// Put the current thread to sleep for a supplied number of milliseconds
void			md_thread_sleep(uint64_t milliseconds);

// ### Mutex ###
typedef struct md_mutex_t {
    union {
        void* _align;
        char _data[64];
    };
} md_mutex_t;

md_mutex_t md_mutex_create(void);
bool md_mutex_init(md_mutex_t* mutex);
bool md_mutex_destroy(md_mutex_t* mutex);

bool md_mutex_lock(md_mutex_t* mutex);

bool md_mutex_try_lock(md_mutex_t* mutex);
bool md_mutex_unlock(md_mutex_t* mutex);

// ### Semaphore ###
typedef struct md_semaphore_t {
    void* _data[4];
} md_semaphore_t;

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

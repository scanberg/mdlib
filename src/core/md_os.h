#pragma once

#include <core/md_platform.h>

#if MD_PLATFORM_UNIX
#include <pthread.h>
#endif

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

// Collection of OS specific things

struct  md_allocator_i;

// ### TICKS ###
// A reading of the platform's monotonic high-resolution counter, in ticks whose rate is
// machine-defined: nanoseconds on unix, QueryPerformanceCounter ticks on windows. Deliberately
// not called a timestamp - it is not anchored to any epoch and cannot be rendered as a date,
// so an individual reading carries no meaning on its own. Subtract two of them and convert the
// difference with md_tick_to_*; those are the only calls that turn ticks into time.
//
// Plain arithmetic on ticks is intended. Do not assume the rate: a delay in milliseconds is not
// a number of ticks and adding one to the other is a bug, however plausible it reads.
typedef int64_t md_tick_t;

// ### THREAD ###
typedef void (*md_thread_exit) (void *data);
typedef void (*md_thread_entry)(void *data);
typedef struct md_thread_t md_thread_t;
typedef uint64_t md_thread_id_t;

#if MD_PLATFORM_WINDOWS
#define md_thread_key_t uint32_t
#elif MD_PLATFORM_UNIX
#define md_thread_key_t pthread_key_t
#endif

// ### MUTEX ###
typedef struct md_mutex_t {
    union {
        void* _align;
        char _data[64];
    };
} md_mutex_t;

// ### Semaphore ###
typedef struct md_semaphore_t {
    void* _data[4];
} md_semaphore_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### OS INFO ###
typedef struct md_os_sys_info_t {
    int num_physical_cores;
    int num_virtual_cores;
    size_t physical_ram_bytes;
    // Add more fields as needed
} md_os_sys_info_t;

// Fills the md_os_info_t struct with system information. Returns true on success.
bool md_os_sys_info_query(md_os_sys_info_t* info);

// ### ENVIRONMENT ###
// The value of a variable of the process's environment, NULL when it is not set. The pointer is the environment's:
// valid until the variable is set again.
const char* md_os_env_get(const char* name);

// Sets a variable of the process's environment to value, or removes it (value NULL), for this process and what it
// starts from now on. True on success. Not thread safe with respect to other changes of or reads from the environment:
// for setup, tests and tools, not for state that threads share.
bool md_os_env_set(const char* name, const char* value);

// ### PATH ###
// Paths, those that ask the file system included, are in md_path.h.

// ### TICKS ###
md_tick_t md_tick_now(void);

// Each takes a difference between two readings, not a single reading.
double  md_tick_to_nanoseconds(md_tick_t ticks);
double  md_tick_to_milliseconds(md_tick_t ticks);
double  md_tick_to_seconds(md_tick_t ticks);

// ### MEMORY ###
// This exposes system calls to the OS virtual memory allocator.
// If you are using this, the assumption is that you know what you are doing
// pointer and size should be aligned on page_size, otherwise it will most certainly crash on some OS

// ### VIRTUAL MEMORY ###
size_t  md_vm_page_size(void);
void*	md_vm_reserve(size_t size);
void	md_vm_release (void* ptr, size_t size);
void	md_vm_commit  (void* ptr, size_t size);
void	md_vm_decommit(void* ptr, size_t size);

// ### THREAD ###
md_thread_t*    md_thread_create(md_thread_entry func, void* data);
void			md_thread_detach(md_thread_t* thread);
bool			md_thread_join(md_thread_t* thread);

// Register callbacks.
bool			md_thread_key_create(md_thread_key_t* key, md_thread_exit callback);
void            md_thread_key_delete(md_thread_key_t key);
void            md_thread_key_set_value(md_thread_key_t key, void* value);

// Get ID from supplied thread object.
md_thread_id_t	md_thread_get_id(md_thread_t* thread);

// Get ID of current thread.
md_thread_id_t	md_thread_id(void);

// Put the current thread to sleep for a supplied number of milliseconds.
void			md_thread_sleep(size_t milliseconds);

// ### MUTEX ###
md_mutex_t md_mutex_create(void);
bool md_mutex_init(md_mutex_t* mutex);
bool md_mutex_destroy(md_mutex_t* mutex);

bool md_mutex_lock(md_mutex_t* mutex);
bool md_mutex_try_lock(md_mutex_t* mutex);
bool md_mutex_unlock(md_mutex_t* mutex);

// ### FILE ###
// Files are in md_file.h.

// ### Semaphore ###
// @NOTE: The size type does not reflect the actual size of the semaphore which is probably (int).
// But it is used to make the API more consistent.
md_semaphore_t md_semaphore_create(size_t initial_count);
bool md_semaphore_init(md_semaphore_t* semaphore, size_t initial_count);
bool md_semaphore_destroy(md_semaphore_t* semaphore);

bool md_semaphore_aquire(md_semaphore_t* semaphore);
bool md_semaphore_try_aquire(md_semaphore_t* semaphore);
bool md_semaphore_try_aquire_n(md_semaphore_t* semaphore, size_t count);

bool md_semaphore_query_count(md_semaphore_t* semaphore, size_t* count);

bool md_semaphore_release(md_semaphore_t* semaphore);
bool md_semaphore_release_n(md_semaphore_t* semaphore, size_t count);

#ifdef __cplusplus
}
#endif

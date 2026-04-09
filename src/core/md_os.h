#pragma once

#include <core/md_platform.h>
#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

// Collection of OS specific things

struct  md_allocator_i;

// ### TIME ###
// This represents a os specific time stamp with the highest precision available (usually nanoseconds)
// It is possible to subtract or perform simple arithmetic directly on the timestamp then
// Extract seconds or milliseconds from it.
typedef int64_t md_timestamp_t;

// ### THREAD ###
typedef void (*md_thread_exit) (void *data);
typedef void (*md_thread_entry)(void *data);
typedef struct md_thread_t md_thread_t;
typedef uint64_t md_thread_id_t;

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

// ### FILE ###
typedef int64_t md_file_offset_t;
typedef int64_t md_file_time_t;

typedef enum md_file_info_flags_t {
    MD_FILE_INFO_IS_DIRECTORY = 1u << 0,
    MD_FILE_INFO_IS_REGULAR   = 1u << 1,
    MD_FILE_INFO_IS_READONLY  = 1u << 2,
} md_file_info_flags_t;
ENUM_FLAGS(md_file_info_flags_t)

typedef struct md_file_info_t {
    md_file_offset_t size;
    md_file_time_t modified_time;
    md_file_time_t accessed_time;
    uint32_t flags;
} md_file_info_t;

typedef struct md_file_t {
#if MD_PLATFORM_WINDOWS
    void* handle;
#elif MD_PLATFORM_UNIX
    int fd;
#else
    #error "Unsupported platform"
#endif
    uint32_t flags;
} md_file_t;

typedef enum md_file_flags_t {
    MD_FILE_READ     = 1u << 0,
    MD_FILE_WRITE    = 1u << 1,
    MD_FILE_APPEND   = 1u << 2,
    MD_FILE_CREATE   = 1u << 3,
    MD_FILE_TRUNCATE = 1u << 4,
} md_file_flags_t;
ENUM_FLAGS(md_file_flags_t)

typedef enum md_file_seek_origin_t {
    MD_FILE_BEG = 0,
    MD_FILE_CUR = 1,
    MD_FILE_END = 2,
} md_file_seek_origin_t;

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

// ### PATH ###
// Gets the current working directory
// Writes into a buffer (up to buf_cap) and returns the number of characters written
size_t md_path_write_cwd(char* buf, size_t buf_cap);
bool  md_path_set_cwd(str_t path);

// Gets the path to the executable
// Writes into a buffer (up to buf_cap) and returns the number of characters written
size_t md_path_write_exe(char* buf, size_t buf_cap);

// Gets the path to the user directory
// Writes into a buffer (up to buf_cap) and returns the number of characters written
// On windows this would correspond to the user's home directory e.g. C:\Users\username
// On linux this would correspond to the user's home directory e.g. /home/username
// On mac this would correspond to the user's home directory e.g. /Users/username
size_t md_path_write_user_dir(char* buf, size_t buf_cap);

// Writes a canonical path to buf.
size_t md_path_write_canonical(char* buf, size_t buf_cap, str_t path);
str_t md_path_make_canonical(str_t path, struct md_allocator_i* alloc);

// Writes a relative path from path_from to path_to to buf.
size_t md_path_write_relative(char* buf, size_t buf_cap, str_t path_from, str_t path_to);
str_t md_path_make_relative(str_t path_from, str_t path_to, struct md_allocator_i* alloc);

// Checks if a path to a file or directory is valid: i.e. it points to an actual file / dir on disk.
bool md_path_is_valid(str_t path);
bool md_path_is_directory(str_t path);

// ### TIME ###
md_timestamp_t md_time_now(void);

double  md_time_as_nanoseconds(md_timestamp_t t);
double  md_time_as_milliseconds(md_timestamp_t t);
double  md_time_as_seconds(md_timestamp_t t);

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
bool			md_thread_on_exit(md_thread_exit callback);

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
// Open semantics:
// - MD_FILE_APPEND implies write access and forces sequential writes to append.
// - MD_FILE_CREATE creates the file if it does not exist.
// - MD_FILE_TRUNCATE truncates an existing file and requires write access.
// - MD_FILE_APPEND and MD_FILE_TRUNCATE are mutually exclusive.
// - A zero-initialized md_file_t represents an invalid or unopened file.
// returns true on success and false on failure, in which case out_file is not modified.
bool md_file_open(md_file_t* out_file, str_t filename, md_file_flags_t flags);

// Checks if the file handle is valid
bool md_file_valid(md_file_t file);

// Closes the file and invalidates the handle on success.
bool md_file_close(md_file_t* file);

// Queries and stream-position operations.
bool md_file_info_extract(md_file_t file, md_file_info_t* out_info);
bool md_file_info_extract_from_path(str_t filename, md_file_info_t* out_info);

md_file_offset_t md_file_tell(md_file_t file);
md_file_offset_t md_file_size(md_file_t file);
bool md_file_seek(md_file_t file, md_file_offset_t offset, md_file_seek_origin_t origin);

// Sequential I/O. Returns the number of bytes transferred and may return less than requested.
size_t md_file_read(md_file_t file, void* ptr, size_t num_bytes);
size_t md_file_write(md_file_t file, const void* ptr, size_t num_bytes);

// Positional I/O. Does not modify the current file position.
// md_file_write_at rejects handles opened with MD_FILE_APPEND because append-only and positional
// write semantics conflict across platforms.
size_t md_file_read_at(md_file_t file, md_file_offset_t offset, void* ptr, size_t num_bytes);
size_t md_file_write_at(md_file_t file, md_file_offset_t offset, const void* ptr, size_t num_bytes);

size_t md_file_printf(md_file_t file, const char* format, ...);

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

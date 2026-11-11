#pragma once

#include <core/md_platform.h>
#include <core/md_str.h>

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

// Files, read and written through the platform's own API. Implemented in md_os.c.

struct md_allocator_i;

typedef int64_t md_file_offset_t;
typedef int64_t md_file_time_t;

typedef enum md_file_info_flags_t {
    MD_FILE_INFO_IS_DIRECTORY = 1u << 0,
    MD_FILE_INFO_IS_REGULAR   = 1u << 1,
    MD_FILE_INFO_IS_READONLY  = 1u << 2,
} md_file_info_flags_t;
MD_ENUM_FLAGS(md_file_info_flags_t)

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
    MD_FILE_READ       = 1u << 0,
    MD_FILE_WRITE      = 1u << 1,
    MD_FILE_APPEND     = 1u << 2,
    MD_FILE_CREATE     = 1u << 3,
    MD_FILE_TRUNCATE   = 1u << 4,

    // Positional only: the handle has no position of its own. md_file_read_at and md_file_write_at may be called on it
    // from several threads at once, and run concurrently; md_file_read, md_file_write, md_file_seek and md_file_tell
    // fail on it. Not with MD_FILE_APPEND.
    MD_FILE_CONCURRENT = 1u << 5,

    // Hints of how the file will be read, which change no result: front to back (the OS reads further ahead), or in
    // jumps (it reads no further than asked). At most one of them.
    MD_FILE_SEQUENTIAL = 1u << 6,
    MD_FILE_RANDOM     = 1u << 7,
} md_file_flags_t;
MD_ENUM_FLAGS(md_file_flags_t)

typedef enum md_file_seek_origin_t {
    MD_FILE_BEG = 0,
    MD_FILE_CUR = 1,
    MD_FILE_END = 2,
} md_file_seek_origin_t;

#ifdef __cplusplus
extern "C" {
#endif

// Open semantics:
// - MD_FILE_APPEND implies write access and forces sequential writes to append.
// - MD_FILE_CREATE creates the file if it does not exist.
// - MD_FILE_TRUNCATE truncates an existing file and requires write access.
// - MD_FILE_APPEND and MD_FILE_TRUNCATE are mutually exclusive.
// - MD_FILE_CONCURRENT makes a handle positional only, shared by threads (see md_file_flags_t); not with MD_FILE_APPEND.
// - MD_FILE_SEQUENTIAL and MD_FILE_RANDOM are mutually exclusive hints.
// - A zero-initialized md_file_t represents an invalid or unopened file.
// returns true on success and false on failure, in which case out_file is not modified.
bool md_file_open(md_file_t* out_file, md_str_t filename, md_file_flags_t flags);

// Checks if the file handle is valid
bool md_file_valid(md_file_t file);

// Closes the file and invalidates the handle on success.
bool md_file_close(md_file_t* file);

// Queries and stream-position operations.
bool md_file_info_extract(md_file_t file, md_file_info_t* out_info);
bool md_file_info_extract_from_path(md_str_t filename, md_file_info_t* out_info);

// The position: not of a handle opened with MD_FILE_CONCURRENT, which has none (-1 and false)
md_file_offset_t md_file_tell(md_file_t file);
md_file_offset_t md_file_size(md_file_t file);
bool md_file_seek(md_file_t file, md_file_offset_t offset, md_file_seek_origin_t origin);

// Sequential I/O. Returns the number of bytes transferred and may return less than requested. Not on a handle opened
// with MD_FILE_CONCURRENT (0).
size_t md_file_read(md_file_t file, void* ptr, size_t num_bytes);
size_t md_file_write(md_file_t file, const void* ptr, size_t num_bytes);

// Positional I/O. Does not modify the current file position, on any platform.
// On a handle opened with MD_FILE_CONCURRENT they may be called from several threads at once and run concurrently (on
// Windows the handle is overlapped; elsewhere pread and pwrite are). On any other handle, one call at a time.
// md_file_write_at rejects handles opened with MD_FILE_APPEND because append-only and positional
// write semantics conflict across platforms.
size_t md_file_read_at(md_file_t file, md_file_offset_t offset, void* ptr, size_t num_bytes);
size_t md_file_write_at(md_file_t file, md_file_offset_t offset, const void* ptr, size_t num_bytes);

size_t md_file_printf(md_file_t file, const char* format, ...);

// Reads a whole file into memory allocated from alloc, zero terminated (not counted in len).
// Returns an empty string on failure.
md_str_t md_load_textfile(md_str_t filename, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

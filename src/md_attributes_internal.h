#pragma once

// Shared between the attribute table (md_attributes.c) and the extraction context (md_system.c).
// Not part of the public API.

#include <md_attributes.h>
#include <core/md_os.h>

#define MD_ATTRIBUTE_IO_MAX_FILES 16
#define MD_ATTRIBUTE_IO_MAX_BLOCKS 4

// A small cache of open files, keyed by path. A handful is enough: an extraction context reads one
// run, and a run is one file per source. The limit is what keeps a long ensemble, one context per
// thread, from walking into the process limit on open descriptors.
struct md_attribute_io_t {
    struct {
        uint64_t  hash;       // of the path; 0 marks an empty slot
        md_str_t     path;       // owned by alloc
        md_file_t file;
        uint64_t  last_use;
    } slot[MD_ATTRIBUTE_IO_MAX_FILES];
    // Decoded blocks of the files (md_attribute_io_block), from the heap: a block of a chunked HDF5 dataset can be
    // far larger than what the context's arena is meant to hold
    struct {
        uint64_t  key;        // 0 marks an empty slot
        void*     data;
        size_t    bytes;      // valid
        size_t    capacity;   // allocated
        uint64_t  last_use;
    } block[MD_ATTRIBUTE_IO_MAX_BLOCKS];
    uint64_t               tick;
    struct md_allocator_i* alloc;
    const struct md_executor_t* executor;   // what providers spread heavy work over (md_attribute_io_executor)
};

void md_attribute_io_close_all(md_attribute_io_t* io);

// The contiguous window a slice selects, in elements. Logs and returns false when it does not apply.
bool md_attribute_slice_window(size_t* out_first, size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice);

// The slice in its STORED type and unit, unconverted: resident storage is copied, a provider asked.
bool md_attribute_read_stored(void* dst, const md_attribute_t* attr, md_attribute_slice_t slice, md_attribute_io_t* io);

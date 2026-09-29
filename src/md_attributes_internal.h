#pragma once

// Shared between the attribute table (md_attributes.c) and the extraction context (md_system.c).
// Not part of the public API.

#include <md_attributes.h>
#include <core/md_os.h>

#define MD_ATTRIBUTE_IO_MAX_FILES 16

// A small cache of open files, keyed by path. A handful is enough: an extraction context reads one
// run, and a run is one file per source. The limit is what keeps a long ensemble, one context per
// thread, from walking into the process limit on open descriptors.
struct md_attribute_io_t {
    struct {
        uint64_t  hash;       // of the path; 0 marks an empty slot
        str_t     path;       // owned by alloc
        md_file_t file;
        uint64_t  last_use;
    } slot[MD_ATTRIBUTE_IO_MAX_FILES];
    uint64_t               tick;
    struct md_allocator_i* alloc;
};

void md_attribute_io_close_all(md_attribute_io_t* io);

// The contiguous window a slice selects, in elements. Logs and returns false when it does not apply.
bool md_attribute_slice_window(size_t* out_first, size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice);

// md_attribute_extract_f32 with an io to hand to a provider.
size_t md_attribute_extract_io_f32(float dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, md_attribute_io_t* io);

// The slice in its STORED type and unit, unconverted: resident storage is copied, a provider asked.
bool md_attribute_read_stored(void* dst, const md_attribute_t* attr, md_attribute_slice_t slice, md_attribute_io_t* io);

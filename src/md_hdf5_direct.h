#pragma once

// Reading the frames of an HDF5 dataset without HDF5.
//
// The HDF5 library is not concurrent: every call is made under the one mdlib lock (md_hdf5.h), and even a build of it
// that is thread safe takes a global lock of its own around each call. A trajectory read frame by frame through it is
// read one frame at a time, however many threads ask.
//
// What a frame costs is not HDF5's bookkeeping, though, but reading its bytes and decoding them, and neither needs
// HDF5. So the storage of a dataset is described once, under the lock: its element type, its layout, the filters its
// chunks pass through, and where each chunk is in the file. A frame is then read from that description alone, on
// whichever thread asks: the chunks it touches are read through the provider's io, the filters undone here, and the
// elements converted to float. A chunk holding several frames is decoded once and kept by the io for the frames after.
//
// Decoded here: contiguous storage and chunked storage of any chunk shape, IEEE floats and integers of either byte
// order, and the filters shuffle, deflate and fletcher32 (whose checksum is verified, as HDF5 does). A dataset stored any other
// way - another filter, a compact or virtual layout, external storage, an element type that is not one of these - is
// read through HDF5 under the lock instead, as is a chunk that was never written.
//
// Requires MD_HDF5.

#include <hdf5.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>

struct md_allocator_i;
struct md_attribute_io_t;

#ifdef __cplusplus
extern "C" {
#endif

#define MD_HDF5_DIRECT_MAX_RANK    4
#define MD_HDF5_DIRECT_MAX_FILTERS 6

typedef enum md_hdf5_direct_mode_t {
    MD_HDF5_DIRECT_HDF5 = 0,        // read through HDF5, under the lock
    MD_HDF5_DIRECT_CONTIGUOUS,
    MD_HDF5_DIRECT_CHUNKED,
} md_hdf5_direct_mode_t;

typedef struct md_hdf5_direct_chunk_t {
    int64_t  offset;                // in the file; -1 for a chunk never written
    uint32_t size;                  // stored bytes
    uint32_t filter_mask;           // bit i set: filter i was skipped for this chunk
} md_hdf5_direct_chunk_t;

// How a dataset is stored, all in one block that can be copied as it is (into a provider's own state): this header,
// then the chunk table, one entry per chunk of the grid in row major order.
typedef struct md_hdf5_direct_t {
    uint64_t size;                  // bytes of the whole block, header and table
    uint32_t mode;                  // md_hdf5_direct_mode_t
    uint32_t type;                  // the element type, internal
    uint32_t elem_size;             // bytes
    uint32_t big_endian;
    uint32_t rank;
    uint32_t num_filters;
    uint64_t dims[MD_HDF5_DIRECT_MAX_RANK];
    uint64_t chunk[MD_HDF5_DIRECT_MAX_RANK];       // CHUNKED: the extent of a chunk
    uint64_t grid[MD_HDF5_DIRECT_MAX_RANK];        // CHUNKED: chunks along each axis
    uint32_t filter[MD_HDF5_DIRECT_MAX_FILTERS];   // H5Z filter ids, in the order they are applied when writing
    uint64_t base;                  // CONTIGUOUS: where element 0 is
    uint64_t num_chunks;            // CHUNKED: entries of the table
    md_hdf5_direct_chunk_t table[];
} md_hdf5_direct_t;

// The storage of dataset dset of file, allocated from alloc (and freed by its size). The lock is held. Never NULL for
// a dataset of a rank from 1 to MD_HDF5_DIRECT_MAX_RANK: what is not decoded here is described as MD_HDF5_DIRECT_HDF5.
// NULL for anything else.
md_hdf5_direct_t* md_hdf5_direct_describe(hid_t file, hid_t dset, struct md_allocator_i* alloc);

// Elements [first, first + count) of frame 'frame' - index 'frame' along the first axis, the axes after it flattened in
// row major order - as floats, each multiplied by scale. Without the lock unless the storage is read through HDF5, from
// the file at path and its dataset 'dataset'. io may be NULL; without it nothing is kept between calls.
bool md_hdf5_direct_read_f32(float* out, const md_hdf5_direct_t* d, str_t path, const char* dataset, size_t frame, size_t first, size_t count, double scale, struct md_attribute_io_t* io);

// Where the stored bytes of each frame are, for a run's frame table (md_run_desc_t): the frame's own bytes when they are
// stored as they are (contiguous, or in a chunk of whole frames without a filter), the chunk holding the frame when it
// is a filtered chunk of whole frames; -1 and 0 otherwise.
void md_hdf5_direct_frame_spans(int64_t* offsets, int64_t* sizes, const md_hdf5_direct_t* d, size_t num_frames);

// The Fletcher-32 checksum of HDF5's fletcher32 filter (H5_checksum_fletcher32): the same value, 16 words at a time
uint32_t md_hdf5_fletcher32(const void* data, size_t bytes);

// Whether the frames of the dataset are read without HDF5
static inline bool md_hdf5_direct_is_direct(const md_hdf5_direct_t* d) {
    return d && d->mode != MD_HDF5_DIRECT_HDF5;
}

#ifdef __cplusplus
}
#endif

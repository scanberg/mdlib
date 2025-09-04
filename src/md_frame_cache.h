#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_trajectory.h>

struct md_allocator_i;
struct md_frame_cache_lock_t;

// This is a simple N-way associative cache set for loading trajectory frames.
// N-way associative meaning that each frame index will end up in 1 of N possible locations.
// The reasoning behind this is to avoid having frame cache evictions from when accessing frames in parallel.

typedef struct md_frame_data_t {
    md_trajectory_frame_header_t header;
    float* x;
    float* y;
    float* z;
} md_frame_data_t;

struct md_slot_header_t;

typedef struct md_frame_cache_t {
    struct {
        size_t count;
        struct md_semaphore_t*   lock;
        struct md_slot_header_t* header;
        struct md_frame_data_t*  data;
    } slot;

    struct md_trajectory_i* traj;
    struct md_allocator_i* alloc;
    char* buf;

    uint64_t magic;
} md_frame_cache_t;

#ifdef __cplusplus
extern "C" {
#endif


// Initialize a frame cache with supplied trajectory, allocator and number of cache frames it should at least contain.
// The actual number of frames which is allocated and used may differ from the supplied number depending on the associativity of the cache.
bool md_frame_cache_init(md_frame_cache_t* cache, struct md_trajectory_i* traj, struct md_allocator_i* alloc, size_t num_cache_frames);
void md_frame_cache_free(md_frame_cache_t* cache);

// Clears the cache by whiping the headers
void md_frame_cache_clear(md_frame_cache_t* cache);

size_t md_frame_cache_num_frames(const md_frame_cache_t* cache);

// This will load the data from a frame into the memory of the supplied pointers.
// The pointers are each optional and if the value of NULL is passed for that pointer, then no data will be written to that adress.
// If the frame is not already in the cache, it will be loaded into the cache and if the postprocess function is supplied, that will be called.
bool md_frame_cache_load_frame_data(md_frame_cache_t* cache, int64_t frame_idx, float* x, float* y, float* z, md_unitcell_t* cell, double* timestamp);

// ### DANGER ZONE ###
// // These are operations which should be handled with care since they expose explicit locks for frames
// This is essentially the same operation as load_frame_data, but without extracting any data
// This also blocks the current thread from executing code. So this is more intended to be used in a parallel
// thread to load 'upcomming' frame data into the cache while operating on other frames.

// Each frame index has an individual lock associated with it in order to synchronize access.
// It is crucial to release the lock when you are done operating on the frame_data.
void md_frame_cache_frame_lock_release(struct md_frame_cache_lock_t* lock);

// Try to find a frame within the cache, if the frame is missing, it reserves space for the frame.
// Returns True if the frame already exists
bool md_frame_cache_find_or_reserve(md_frame_cache_t* cache, int64_t frame_idx, md_frame_data_t** frame_data, struct md_frame_cache_lock_t** frame_lock);

// Try to find a frame within the cache
// Returns true if the frame exists, and its data is provided in frame data together with an active frame lock that has to be released by the callee.
// If the frame does not exist, it simply returns false and no action is required by the user
bool md_frame_cache_find(md_frame_cache_t* cache, int64_t frame_idx, md_frame_data_t** frame_data, struct md_frame_cache_lock_t** frame_lock);

#ifdef __cplusplus
}
#endif

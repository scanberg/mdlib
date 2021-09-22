#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_trajectory_i;
struct md_frame_cache_lock_t;

// This is a simple N-way associative cache set for loading trajectory frames.
// N-way associative meaning that each frame index will end up in 1 of N possible locations.
// The reasoning behind this is to avoid having frame cache evictions from when accessing frames in parallel.

typedef struct md_frame_data_t {
    int frame_index;
    int num_atoms;
    float* x;
    float* y;
    float* z;
    float box[3][3];
    double timestamp;
} md_frame_data_t;

typedef void (md_frame_cache_load_frame_fn)(struct md_trajectory_i* traj, int64_t frame_idx, md_frame_data_t* frame_data, struct md_frame_cache_lock_t* frame_lock, void* user_data);

typedef void (md_frame_cache_postprocess_frame_fn)(md_frame_data_t* frame_data, void* user_data);

struct md_slot_header_t;

typedef struct md_frame_cache_t {
    struct {
        int64_t count;
        struct md_semaphore_t*   lock;
        struct md_slot_header_t* header;
        struct md_frame_data_t*  data;
    } slot;

    struct md_trajectory_i* traj;
    struct md_allocator_i* alloc;
    void*   mem_ptr;
    int64_t mem_bytes;

    uint64_t magic;
} md_frame_cache_t;

void md_frame_cache_init(md_frame_cache_t* cache, struct md_trajectory_i* traj, struct md_allocator_i* alloc, int64_t cache_byte_size);
void md_frame_cache_free(md_frame_cache_t* cache);

// This will load the data from a frame into the memory of the supplied pointers.
// The pointers are each optional and if the value of NULL is passed for that pointer, then no data will be written to that adress.
// If the frame is not already in the cache, it will be loaded into the cache and if the postprocess function is supplied, that will be called.
bool md_frame_cache_load_frame_data(md_frame_cache_t* cache, int64_t frame_idx, float* x, float* y, float* z, float box[3][3], double* timestamp, md_frame_cache_postprocess_frame_fn* postprocess_fn, void* user_data);

// ### DANGER ZONE ###
// // These are operations which should be handled with care since they expose explicit locks for frames
// This is essentially the same operation as load_frame_data, but without extracting any data
// This also blocks the current thread from executing code. So this is more intended to be used in a parallel
// thread to load 'upcomming' frame data into the cache while operating on other frames.
// load_fn is an optional callback function and if set, it will be called when it is time to load the frame from the trajectory and store it into
// the supplied frame_data pointer.

// Each frame index has an individual lock associated with it in order to synchronize access.
// It is crucial to release the lock when you are done operating on the frame_data.
void md_frame_cache_release_frame_lock(struct md_frame_cache_lock_t* lock);

// Reserves a frame in the cache
// Returns true if the data was not already in cache, the user can perform operations on the data, before releasing the frame lock.
bool md_frame_cache_reserve_frame(md_frame_cache_t* cache, int64_t frame_idx, md_frame_data_t** frame_data, struct md_frame_cache_lock_t** frame_lock);

//bool md_frame_cache_fetch_frame_range(md_frame_cache_t* cache, int64_t frame_beg_idx, int64_t frame_end_idx, md_frame_cache_load_frame_fn* load_fn, void* user_data);

#ifdef __cplusplus
}
#endif
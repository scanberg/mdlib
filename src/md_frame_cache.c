#include "md_frame_cache.h"

#include "core/md_allocator.h"
#include "core/md_simd.h"
#include "core/md_log.h"
#include "core/md_sync.h"
#include "md_trajectory.h"
#include "md_util.h"

#include <string.h>

#define CACHE_MAGIC 0x7bacacababca8181
#define CACHE_ASSOCIATIVITY 2
#define CACHE_MEM_ALIGNMENT 16

typedef struct md_slot_header_t {
    uint32_t frame_index;
    uint32_t access_count;
} md_slot_header_t;

void md_frame_cache_init(md_frame_cache_t* cache, md_trajectory_i* traj, md_allocator_i* alloc, int64_t num_cached_frames) {
    ASSERT(cache);
    ASSERT(traj);
    ASSERT(num_cached_frames >= 0);

    if (num_cached_frames == 0) {
        num_cached_frames = traj->num_frames;
    }

    // We want to ensure that there are enough padding for each frame to avoid overlap if one wants to do full-width simd stores.
    const int64_t num_atoms = ROUND_UP(traj->num_atoms, md_simd_width); 
    const int64_t bytes_per_frame = sizeof(md_semaphore_t) + sizeof(md_slot_header_t) + sizeof(md_frame_data_t) + num_atoms * sizeof(float) * 3;
    const int64_t num_slots = ROUND_UP(num_cached_frames, CACHE_ASSOCIATIVITY); // This needs to be divisible by N for N-way associativity.

    cache->traj = traj;
    cache->mem_bytes = num_slots * bytes_per_frame;
    cache->mem_ptr = md_alloc(alloc, num_slots * bytes_per_frame + CACHE_MEM_ALIGNMENT);
    ASSERT(cache->mem_ptr);
    cache->slot.count  = num_slots;
    cache->slot.lock   = (md_semaphore_t*)NEXT_ALIGNED_ADRESS(cache->mem_ptr, CACHE_MEM_ALIGNMENT);
    cache->slot.header = (md_slot_header_t*)(cache->slot.lock + num_slots);
    cache->slot.data   = (md_frame_data_t*)(cache->slot.header + num_slots);

    cache->magic = CACHE_MAGIC;

    float* coord_data = (float*)(cache->slot.data + num_slots);
    for (int64_t i = 0; i < num_slots; ++i) {
        md_semaphore_init(&cache->slot.lock[i], 1);

        cache->slot.header[i].frame_index = 0xFFFFFFFF;
        cache->slot.header[i].access_count = 0;

        cache->slot.data[i].num_atoms = num_atoms;
        cache->slot.data[i].x = coord_data + (num_slots * num_atoms) * 0 + num_atoms * i;
        cache->slot.data[i].y = coord_data + (num_slots * num_atoms) * 1 + num_atoms * i;
        cache->slot.data[i].z = coord_data + (num_slots * num_atoms) * 2 + num_atoms * i;
    }
}

void md_frame_cache_free(md_frame_cache_t* cache) {
    ASSERT(cache);
    ASSERT(cache->magic == CACHE_MAGIC);
    if (cache->mem_ptr) {
        for (int64_t i = 0; i < cache->slot.count; ++i) {
            md_semaphore_destroy(&cache->slot.lock[i]);
        }

        ASSERT(cache->alloc);
        md_free(cache->alloc, cache->mem_ptr, cache->mem_bytes);
    }
    memset(cache, 0, sizeof(md_frame_cache_t));
}

#include <stdio.h>

static inline void md_frame_cache_aquire_frame_lock(struct md_frame_cache_lock_t* lock) {
    ASSERT(lock);
    bool success = md_semaphore_aquire((md_semaphore_t*)lock);
#if 0
    char buf[512] = {0};
    md_thread_id_t id = md_thread_id();
    snprintf(buf, 512, "thread-id: %llu", id);
    md_printf(MD_LOG_TYPE_DEBUG, "AQUIRE LOCK \t%llu: %s", (uint64_t)lock, buf);
#endif
    ASSERT(success);
}

void md_frame_cache_release_frame_lock(struct md_frame_cache_lock_t* lock) {
    ASSERT(lock);
    bool success = md_semaphore_release((md_semaphore_t*)lock);
#if 0
    md_printf(MD_LOG_TYPE_DEBUG, "RELEASE\t%llu\n", (uint64_t)lock);
#endif
    ASSERT(success);
}

static bool find_frame_or_reserve_slot(md_frame_cache_t* cache, int64_t frame_idx, md_frame_data_t** slot_data, struct md_frame_cache_lock_t** slot_lock) {
    ASSERT(cache);
    ASSERT(cache->magic == CACHE_MAGIC);
    ASSERT(cache->traj);
    ASSERT(frame_idx >= 0);

    const int64_t start_slot = ((frame_idx % cache->slot.count) / CACHE_ASSOCIATIVITY) * CACHE_ASSOCIATIVITY;

    int64_t slot_idx = -1;
    int64_t max_count_idx = -1;
    int64_t max_count = -1;

    // If the frame is already in cache -> return data
    for (int64_t i = start_slot; i < start_slot + CACHE_ASSOCIATIVITY; ++i) {
        cache->slot.header[i].access_count += 1;
        if (cache->slot.header[i].access_count > max_count) {
            max_count = cache->slot.header[i].access_count;
            max_count_idx = i;
        }
        if (cache->slot.header[i].frame_index == (uint32_t)frame_idx) {
            slot_idx = i;
        }
    }

    bool found = slot_idx != -1;

    if (slot_idx == -1) {
        slot_idx = max_count_idx;
    }

    ASSERT(slot_idx != -1);

    md_frame_cache_aquire_frame_lock((struct md_frame_cache_lock_t*)&cache->slot.lock[slot_idx]);

    cache->slot.header[slot_idx].access_count = 0;
    cache->slot.header[slot_idx].frame_index = (uint32_t)frame_idx;

    *slot_lock = (struct md_frame_cache_lock_t*)&cache->slot.lock[slot_idx];
    if (slot_data) *slot_data = &cache->slot.data[slot_idx];

    return found;
}

static inline bool load_frame_data(md_trajectory_i* traj, int64_t frame_idx, md_frame_data_t* frame_data) {
    md_trajectory_frame_header_t header = {0};
    if (md_trajectory_load_frame(traj, frame_idx, &header, frame_data->x, frame_data->y, frame_data->z, frame_data->num_atoms)) {
        memcpy(&frame_data->box, header.box, sizeof(frame_data->box));
        frame_data->timestamp = header.timestamp;
        frame_data->num_atoms = header.num_atoms;
        return true;
    }
    return false;
}

bool md_frame_cache_load_frame_data(md_frame_cache_t* cache, int64_t frame_idx, float* x, float* y, float* z, float box[3][3], double* timestamp) {
    ASSERT(cache);
    ASSERT(cache->magic == CACHE_MAGIC);
    ASSERT(cache->traj);
    ASSERT(0 <= frame_idx && frame_idx < cache->traj->num_frames);

    md_frame_data_t* data = NULL;
    struct md_frame_cache_lock_t* lock = NULL;
    if (!find_frame_or_reserve_slot(cache, frame_idx, &data, &lock)) {
        load_frame_data(cache->traj, frame_idx, data);
    }

    bool result = false;
    if (x || y || z || box || timestamp) {
        ASSERT(data);
        ASSERT(lock);

        if (x)   memcpy(x,   data->x,   data->num_atoms * sizeof(float));
        if (y)   memcpy(y,   data->y,   data->num_atoms * sizeof(float));
        if (z)   memcpy(z,   data->z,   data->num_atoms * sizeof(float));
        if (box) memcpy(box, data->box, sizeof(data->box));
        if (timestamp) *timestamp = data->timestamp;

        result = true;
    }

    md_frame_cache_release_frame_lock(lock);
    return result;
}

bool md_frame_cache_fetch_frame(md_frame_cache_t* cache, int64_t frame_idx, md_frame_cache_load_frame_fn* load_fn, void* user_data) {
    ASSERT(cache);
    ASSERT(cache->magic == CACHE_MAGIC);
    ASSERT(cache->traj);
    ASSERT(0 <= frame_idx && frame_idx < cache->traj->num_frames);

    md_frame_data_t* data = NULL;
    struct md_frame_cache_lock_t* lock = NULL;
    bool result = true;
    if (find_frame_or_reserve_slot(cache, frame_idx, &data, &lock)) {
        // Nothing to do here!
        md_frame_cache_release_frame_lock(lock);
    }
    else {
        if (load_fn) {
            load_fn(cache->traj, frame_idx, data, (struct md_frame_cache_lock_t*)lock, user_data);
        }
        else {
            load_frame_data(cache->traj, frame_idx, data);
            md_frame_cache_release_frame_lock(lock);
        }
    }
    return result;
}

bool md_frame_cache_fetch_frame_range(md_frame_cache_t* cache, int64_t frame_beg_idx, int64_t frame_end_idx, md_frame_cache_load_frame_fn* load_fn, void* user_data) {
    ASSERT(cache);
    ASSERT(cache->magic == CACHE_MAGIC);
    ASSERT(cache->traj);
    ASSERT(0 <= frame_beg_idx && frame_end_idx <= cache->traj->num_frames);

    md_frame_data_t* data = NULL;
    struct md_frame_cache_lock_t* lock = NULL;
    bool result = true;
    for (int64_t i = frame_beg_idx; i < frame_end_idx; ++i) {
        if (find_frame_or_reserve_slot(cache, i, &data, &lock)) {
            // Nothing to do here!
            md_frame_cache_release_frame_lock(lock);
        }
        else {
            if (load_fn) {
                load_fn(cache->traj, i, data, (struct md_frame_cache_lock_t*)lock, user_data);
            }
            else {
                load_frame_data(cache->traj, i, data);
                md_frame_cache_release_frame_lock(lock);
            }
        }
    }
    return result;
}
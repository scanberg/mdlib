#include "utest.h"
#include <string.h>

#include <md_trajectory.h>
#include <md_molecule.h>
#include <md_frame_cache.h>
#include <md_xtc.h>
#include <md_gro.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#define NUM_THREADS 8

typedef struct thread_data_t {
    md_frame_cache_t* cache;
    float* ref_coords;
    md_unit_cell_t ref_cell;
    int corrupt_count;
    int thread_rank;
} thread_data_t;

void thread_func(void* user_data) {
    ASSERT(user_data);
    thread_data_t* data = (thread_data_t*)user_data;

    const int64_t num_atoms = md_trajectory_num_atoms(data->cache->traj);
    const int64_t num_frames = md_trajectory_num_frames(data->cache->traj);
    const int64_t frame_stride = num_atoms * 3;
    
    const int64_t mem_size = num_atoms * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_temp_allocator(), mem_size);
    float *x = (float*)mem_ptr + num_atoms * 0;
    float *y = (float*)mem_ptr + num_atoms * 1;
    float *z = (float*)mem_ptr + num_atoms * 2;

    md_unit_cell_t cell;

    // Load frame by frame and memcmp to reference
    // Offset the frame_index 'randomly' so we really stress the cache
    const int64_t offset = (int64_t)(md_thread_id() % num_frames);
    for (int64_t i = 0; i < (int64_t)num_frames; ++i) {
        int64_t frame_idx = (offset + i) % num_frames;
        md_frame_cache_load_frame_data(data->cache, frame_idx, x, y, z, &cell, NULL);

        const float* ref_x = data->ref_coords + frame_stride * frame_idx + num_atoms * 0;
        const float* ref_y = data->ref_coords + frame_stride * frame_idx + num_atoms * 1;
        const float* ref_z = data->ref_coords + frame_stride * frame_idx + num_atoms * 2;

        bool corrupt = false;
        if (memcmp(ref_x, x, num_atoms * sizeof(float)) != 0) corrupt = true;
        if (memcmp(ref_y, y, num_atoms * sizeof(float)) != 0) corrupt = true;
        if (memcmp(ref_z, z, num_atoms * sizeof(float)) != 0) corrupt = true;
        if (memcmp(&data->ref_cell, &cell, sizeof(cell)) != 0) corrupt = true;

        if (corrupt) {
            data->corrupt_count += 1;
        }
    }

    md_free(md_get_temp_allocator(), mem_ptr, mem_size);
}

UTEST(frame_cache, parallel_workload) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));

    md_molecule_t mol = {0};
    md_gro_data_t gro = {0};
    md_frame_cache_t cache = {0};
    md_trajectory_i* traj = md_xtc_trajectory_create(STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc"), alloc, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);

    ASSERT_TRUE(md_gro_data_parse_file(&gro, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"), alloc));
    ASSERT_TRUE(md_gro_molecule_init(&mol, &gro, alloc));
    ASSERT_TRUE(traj);

    const int64_t num_atoms = md_trajectory_num_atoms(traj);
    const int64_t num_frames = md_trajectory_num_frames(traj);
    // We cache much less of the frames so we get some eviction and contention going
    ASSERT_TRUE(md_frame_cache_init(&cache, traj, alloc, num_frames / 8));

    const int64_t frame_stride = num_atoms * 3;
    float* ref_coords = md_alloc(alloc, num_frames * num_atoms * 3 * sizeof(float));
    md_unit_cell_t* ref_cells = md_alloc(alloc, num_frames * sizeof(md_unit_cell_t));

    for (int64_t i = 0; i < num_frames; ++i) {
        float* ref_x = ref_coords + frame_stride * i + num_atoms * 0;
        float* ref_y = ref_coords + frame_stride * i + num_atoms * 1;
        float* ref_z = ref_coords + frame_stride * i + num_atoms * 2;
        EXPECT_TRUE(md_frame_cache_load_frame_data(&cache, i, ref_x, ref_y, ref_z, &ref_cells[i], NULL));
    }

    thread_data_t thread_data[NUM_THREADS] = {0};
    md_thread_t* threads[NUM_THREADS] = {0};

    for (int pass = 0; pass < 10; ++pass) {
        for (int i = 0; i < NUM_THREADS; ++i) {
            thread_data[i].cache = &cache;
            thread_data[i].ref_coords = ref_coords;
            thread_data[i].ref_cell   = ref_cells[i];
            thread_data[i].corrupt_count = 0;
            threads[i] = md_thread_create(thread_func, &thread_data[i]);
        }

        int corrupt_count = 0;
        for (int i = 0; i < NUM_THREADS; ++i) {
            md_thread_join(threads[i]);
            corrupt_count += thread_data[i].corrupt_count;
        }

        EXPECT_EQ(0, corrupt_count);
    }

    md_xtc_trajectory_free(traj);
    md_arena_allocator_destroy(alloc);
}

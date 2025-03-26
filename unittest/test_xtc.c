#include "utest.h"
#include <string.h>

#include <md_xtc.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

UTEST(xtc, trajectory_i) {
    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_trajectory_i* traj = md_xtc_trajectory_create(path, md_get_heap_allocator(), MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
    ASSERT_TRUE(traj);

    const int64_t num_atoms = md_trajectory_num_atoms(traj);
    const int64_t num_frames = md_trajectory_num_frames(traj);

    EXPECT_EQ(num_atoms, 1336);
    EXPECT_EQ(num_frames, 501);

    const int64_t mem_size = num_atoms * 3 * sizeof(float);
    void* mem_ptr = md_alloc(md_get_heap_allocator(), mem_size);
    float *x = (float*)mem_ptr;
    float *y = (float*)mem_ptr + num_atoms * 1;
    float *z = (float*)mem_ptr + num_atoms * 2;

    md_trajectory_frame_header_t header;

    for (int64_t i = 0; i < num_frames; ++i) {
        EXPECT_TRUE(md_trajectory_load_frame(traj, i, &header, x, y, z));
    }

    md_free(md_get_heap_allocator(), mem_ptr, mem_size);
    md_xtc_trajectory_free(traj);
}

UTEST(xtc, read_frame_data) {
    md_allocator_i* alloc = md_get_heap_allocator();

    const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);
    ASSERT_TRUE(file);

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, alloc);
    size_t num_frames = frame_offsets ? md_array_size(frame_offsets) - 1 : 0;

    EXPECT_EQ(num_frames, 501);

    int natoms, step;
    float time, box[3][3];

    md_array(vec3_t) coords = 0;

    for (size_t i = 0; i < num_frames; ++i) {
        md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
        EXPECT_TRUE(md_xtc_read_frame_header(file, &natoms, &step, &time, box));

        if (natoms > 0) {
            if (!coords) {
                md_array_resize(coords, natoms, alloc);
            }
            EXPECT_TRUE(md_xtc_read_frame_coords(file, (float*)coords, md_array_size(coords)));
        }
    }

    md_array_free(frame_offsets, alloc);
    md_array_free(frame_times, alloc);
    md_array_free(coords, alloc);
}


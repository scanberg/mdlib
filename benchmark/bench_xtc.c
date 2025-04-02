#include "ubench.h"

#include <md_xtc.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_arena_allocator.h>

UBENCH_EX(xtc, catalyst) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_LIT(MD_BENCHMARK_DATA_DIR "/catalyst.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    if (!md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena)) {
        MD_LOG_ERROR("Bad");
        return;
    }
    
    static const size_t num_frames = 501;
    static const size_t num_atoms  = 1336;

    int natoms, step;
    float time, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);
    
    UBENCH_DO_BENCHMARK() {
        for (size_t i = 0; i < num_frames; ++i) {
            md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
            md_xtc_read_frame_header(file, &natoms, &step, &time, box);
            md_xtc_read_frame_coords(file, coords, num_atoms);
        }
    }

    md_vm_arena_destroy(arena);
}

#if 1
UBENCH_EX(xtc, amyloid) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_LIT("E:/data/md/amyloid-6T/prod-centered.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    if (!md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena)) {
        MD_LOG_ERROR("Bad");
        return;
    }

    static const size_t num_frames = 5701;
    static const size_t num_atoms  = 161271;

    int natoms, step;
    float time, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);

    UBENCH_DO_BENCHMARK() {
        for (size_t i = 0; i < num_frames; ++i) {
            md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
            md_xtc_read_frame_header(file, &natoms, &step, &time, box);
            md_xtc_read_frame_coords(file, coords, num_atoms);
        }
    }

    md_vm_arena_destroy(arena);
}

UBENCH_EX(xtc, H1N1) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    const str_t path = STR_LIT("E:/data/md/H1N1/H1N1-Mich2015-TRAJECTORY-not_water_not_ions-sk100.xtc");
    md_file_o* file = md_file_open(path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));

    md_array(int64_t) frame_offsets = 0;
    md_array(double)  frame_times = 0;

    if (!md_xtc_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, arena)) {
        MD_LOG_ERROR("Bad");
        return;
    }

    static const size_t num_frames = 71;
    static const size_t num_atoms  = 14009213;

    int natoms, step;
    float time, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);

    UBENCH_DO_BENCHMARK() {
        for (size_t i = 0; i < num_frames; ++i) {
            md_file_seek(file, frame_offsets[i], MD_FILE_BEG);
            md_xtc_read_frame_header(file, &natoms, &step, &time, box);
            md_xtc_read_frame_coords(file, coords, num_atoms);
        }
    }

    md_vm_arena_destroy(arena);
}
#endif

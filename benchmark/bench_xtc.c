#include "ubench.h"

#include <md_xtc.h>
#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_arena_allocator.h>

#include <xdrfile_xtc.h>

#define STR(x) {x"", sizeof(x"")-1};

#define FULL_TEST 0

static const str_t cat_path = STR(MD_BENCHMARK_DATA_DIR "/catalyst.xtc");
static const str_t amy_path = STR(MD_BENCHMARK_DATA_DIR "/amyloid-pftaa.xtc");
static const str_t asp_path = STR(MD_BENCHMARK_DATA_DIR "/aspirin-phospholipase.xtc");
static const str_t ion_path = STR(MD_BENCHMARK_DATA_DIR "/ef.xtc");

UBENCH_EX(xtc, xdr_catalyst) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(cat_path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    XDRFILE* xdr = xdrfile_open(cat_path.ptr, "rb");
    
    static const size_t num_frames = 501;
    static const size_t num_atoms  = 1336;

    int natoms, step;
    float time, prec, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);
    
    UBENCH_DO_BENCHMARK() {
        xdr_seek(xdr, 0, SEEK_SET);
        for (size_t i = 0; i < num_frames; ++i) {
            /* Read one frame of an open xtc file */
            int res = read_xtc(xdr, &natoms, &step, &time, box, (rvec*)coords, &prec);
        }
    }

    xdrfile_close(xdr);

    md_vm_arena_destroy(arena);
}

#if FULL_TEST
UBENCH_EX(xtc, xdr_amyloid) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(amy_path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    XDRFILE* xdr = xdrfile_open(amy_path.ptr, "rb");
    
    static const size_t num_frames = 2345;
    static const size_t num_atoms  = 161742;

    int natoms, step;
    float time, prec, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);
    
    UBENCH_DO_BENCHMARK() {
        xdr_seek(xdr, 0, SEEK_SET);
        for (size_t i = 0; i < num_frames; ++i) {
            /* Read one frame of an open xtc file */
            int res = read_xtc(xdr, &natoms, &step, &time, box, (rvec*)coords, &prec);
        }
    }

    xdrfile_close(xdr);

    md_vm_arena_destroy(arena);
}

UBENCH_EX(xtc, xdr_aspirin) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(asp_path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    XDRFILE* xdr = xdrfile_open(asp_path.ptr, "rb");

    static const size_t num_frames = 601;
    static const size_t num_atoms  = 5015;

    int natoms, step;
    float time, prec, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);
    
    UBENCH_DO_BENCHMARK() {
        xdr_seek(xdr, 0, SEEK_SET);
        for (size_t i = 0; i < num_frames; ++i) {
            /* Read one frame of an open xtc file */
            int res = read_xtc(xdr, &natoms, &step, &time, box, (rvec*)coords, &prec);
        }
    }

    xdrfile_close(xdr);

    md_vm_arena_destroy(arena);
}

UBENCH_EX(xtc, xdr_ion_channel) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(ion_path, MD_FILE_READ | MD_FILE_BINARY);

    if (!file) {
        MD_LOG_ERROR("Bad");
        return;
    }

    UBENCH_SET_BYTES(md_file_size(file));
    md_file_close(file);

    XDRFILE* xdr = xdrfile_open(ion_path.ptr, "rb");

    static const size_t num_frames = 1989;
    static const size_t num_atoms  = 222387;

    int natoms, step;
    float time, prec, box[3][3];

    float* coords = md_vm_arena_push(arena, sizeof(float) * num_atoms * 3);
    
    UBENCH_DO_BENCHMARK() {
        xdr_seek(xdr, 0, SEEK_SET);
        for (size_t i = 0; i < num_frames; ++i) {
            /* Read one frame of an open xtc file */
            int res = read_xtc(xdr, &natoms, &step, &time, box, (rvec*)coords, &prec);
        }
    }

    xdrfile_close(xdr);

    md_vm_arena_destroy(arena);
}
#endif

UBENCH_EX(xtc, xtc_catalyst) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(cat_path, MD_FILE_READ | MD_FILE_BINARY);

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

#if FULL_TEST
UBENCH_EX(xtc, xtc_amyloid) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(amy_path, MD_FILE_READ | MD_FILE_BINARY);

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

    static const size_t num_frames = 2345;
    static const size_t num_atoms  = 161742;

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

UBENCH_EX(xtc, xtc_aspirin) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(asp_path, MD_FILE_READ | MD_FILE_BINARY);

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

    static const size_t num_frames = 601;
    static const size_t num_atoms  = 50515;

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

UBENCH_EX(xtc, xtc_ion_channel) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(1));
    md_file_o* file = md_file_open(ion_path, MD_FILE_READ | MD_FILE_BINARY);

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

    static const size_t num_frames = 1989;
    static const size_t num_atoms  = 222387;

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

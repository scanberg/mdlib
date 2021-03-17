#ifndef _MD_XTC_H_
#define _MD_XTC_H_

#include <stdint.h>
#include <stdbool.h>

struct md_trajectory_i;

#ifdef __cplusplus
extern "C" {
#endif

struct md_xtc_frame {
    void* ptr;
    uint64_t size;
};

bool md_xtc_decompress_frame(float* x, float* y, float* z, float* box[3][3], double* time, void* data_ptr, uint64_t data_len);

bool md_xtc_load_frame(struct md_trajectory_i* traj, void* ptr, uint64_t* size);

struct md_trajectory_i* md_xtc_create_trajectory_from_file(const char* filename_str, uint64_t filename_len);
void md_xtc_close_trajectory(struct md_trajectory_i* traj);

#ifdef __cplusplus
}
#endif

#endif
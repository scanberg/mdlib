#ifndef _MD_GRO_H_
#define _MD_GRO_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

struct md_parse_result {
    int result;
    const char* title;
    float unit_cell[3];

    // Field data
    uint32_t count;
    uint32_t* res_idx;
    char (*res_name[8]);
    char (*name[8]);
    uint32_t* idx;
    float* x;
    float* y;
    float* z;
    float* vx;
    float* vy;
    float* vz;
};

// Parse a text-blob as GRO
struct md_parse_result md_parse_gro(const char* str, uint32_t len, struct md_allocator_i*);

#ifdef __cplusplus
}
#endif

#endif
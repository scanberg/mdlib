#ifndef _MD_PDB_H_
#define _MD_PDB_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

struct md_pdb_result {
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

// Parse a text-blob as PDB
struct md_pdb_result md_parse_pdb(const char* str, uint32_t len, struct md_allocator_i*);

#ifdef __cplusplus
}
#endif

#endif
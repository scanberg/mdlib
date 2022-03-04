#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>

struct md_allocator_i;
struct md_bitfield_t;
struct md_molecule_t;
struct md_script_ir_t;

typedef struct md_filter_result_t {
    int64_t num_bitfields;
    struct md_bitfield_t* bitfields;
    char error_buf[256];
    bool is_dynamic;
} md_filter_result_t;

#ifdef __cplusplus
    extern "C" {
#endif

bool md_filter_evaluate(md_filter_result_t* result, str_t expression, const struct md_molecule_t* mol, const struct md_script_ir_t* ir, struct md_allocator_i* alloc);
bool md_filter_free(md_filter_result_t* result, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

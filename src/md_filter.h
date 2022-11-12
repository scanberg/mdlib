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

// This is more advanced and ugly
bool md_filter_evaluate(md_filter_result_t* result, str_t expression, const struct md_molecule_t* mol, const struct md_script_ir_t* ctx_ir, struct md_allocator_i* alloc);
bool md_filter_free(md_filter_result_t* result, struct md_allocator_i* alloc);

// This is simpler, it evaluates an expression with a molecule.
// ctx_ir is an optional evaluation context (which contains additional symbols and stuff it shold be aware of from some other script compilation)
// err_buf and err_cap are optional and are pointer + length to a string buffer where it can write any evaluation errors of the expression
bool md_filter(struct md_bitfield_t* dst_bf, str_t expression, const struct md_molecule_t* mol, const struct md_script_ir_t* ctx_ir, char* err_buf, int err_cap);

#ifdef __cplusplus
}
#endif

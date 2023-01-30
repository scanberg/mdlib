#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>
#include <core/md_array.h>

struct md_allocator_i;
struct md_bitfield_t;
struct md_molecule_t;
struct md_script_ir_t;

#ifdef __cplusplus 
    extern "C" {
#endif

// This is more advanced and ugly
bool md_filter_evaluate(md_array(struct md_bitfield_t)* bitfields, str_t expression, const struct md_molecule_t* mol, const struct md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, int err_cap, struct md_allocator_i* alloc);

// Simpler version, it evaluates a filter expression for a given molecule and concatenates the result down into a single expression.
// ctx_ir is an optional evaluation context (which contains additional symbols and stuff it shold be aware of from some other script compilation)
// err_buf and err_cap are optional and are pointer + length to a string buffer where it can write any evaluation errors of the expression
bool md_filter(struct md_bitfield_t* bf, str_t expression, const struct md_molecule_t* mol, const struct md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, int err_cap);

#ifdef __cplusplus
}
#endif

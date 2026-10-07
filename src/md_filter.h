#pragma once

#include <core/md_array.h>
#include <core/md_str.h>

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitset_t;
struct md_system_t;
struct md_system_state_t;
struct md_script_ir_t;

#ifdef __cplusplus 
    extern "C" {
#endif

// This is more advanced and ugly
// The sets are pushed to *sets, both array and payloads allocated with alloc.
bool md_filter_evaluate(md_array(struct md_bitset_t)* sets, str_t expression, const struct md_system_t* sys, const struct md_system_state_t* state, const struct md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, size_t err_cap, struct md_allocator_i* alloc);

// Simpler version, it evaluates a filter expression for a given molecule and concatenates the result down into a single expression.
// ctx_ir is an optional evaluation context (which contains additional symbols and stuff it shold be aware of from some other script compilation)
// err_buf and err_cap are optional and are pointer + length to a string buffer where it can write any evaluation errors of the expression
// The result replaces *set: its previous payload is released with alloc, so it must have been produced with
// alloc (or be empty), and the new one is allocated with it.
bool md_filter(struct md_bitset_t* set, str_t expression, const struct md_system_t* sys, const struct md_system_state_t* state, const struct md_script_ir_t* ctx_ir, bool* is_dynamic, char* err_buf, size_t err_cap, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

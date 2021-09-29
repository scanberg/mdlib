#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;
struct md_exp_bitfield_t;
struct md_molecule_t;
struct md_script_ir_t;

typedef struct md_filter_stored_selection_t {
    str_t ident;
    const struct md_exp_bitfield_t*  bitfield;
} md_filter_stored_selection_t;

typedef struct md_filter_context_t {
    const struct md_script_ir_t* ir;
    const struct md_molecule_t* mol;
    struct md_allocator_i* alloc;
    struct {
        int64_t count;
        const md_filter_stored_selection_t* ptr;
    } selection;
} md_filter_context_t;

bool md_filter_evaluate(str_t expression, struct md_exp_bitfield_t* target, md_filter_context_t ctx);

#ifdef __cplusplus
}
#endif

#ifndef _MD_FILTER_H_
#define _MD_FILTER_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_bitfield.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule;

typedef struct md_filter_stored_selection {
    str_t           ident;
    md_bitfield_t   bitfield;
} md_filter_stored_selection_t;

typedef struct md_filter_context {
    const struct md_molecule* mol;
    struct {
        int64_t count;
        const md_filter_stored_selection_t* ptr;
    } selection;
} md_filter_context_t;

bool md_filter_evaluate(str_t expression, md_bitfield_t target, md_filter_context_t ctx);

#ifdef __cplusplus
}
#endif

#endif
#ifndef _MD_FILTER_H_
#define _MD_FILTER_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_log_i;
struct md_molecule;

typedef struct md_filter_context md_filter_context;
typedef struct md_filter_stored_selection md_filter_stored_selection;

struct md_bitfield {
    uint64_t* bits;
    uint64_t bit_count;
};

struct md_filter_stored_selection {
    const char*     ident;
    const uint64_t* bits;
};

struct md_filter_context {
    const struct md_molecule*         mol;
    uint64_t                            bit_count;
    uint64_t*                           bits;
    uint32_t                            sel_count;
    const md_filter_stored_selection* sel;
    const struct md_log_i*            err_log;
};

bool md_filter_apply(const char* expression, const md_filter_context* context);

#ifdef __cplusplus
}
#endif

#endif
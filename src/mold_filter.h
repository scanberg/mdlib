#ifndef _MOLD_FILTER_H_
#define _MOLD_FILTER_H_

#include "mold.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mold_filter_context mold_filter_context;
typedef struct mold_filter_stored_selection mold_filter_stored_selection;

struct mold_filter_stored_selection {
    const char*     ident;
    const uint64_t* bits;
};

struct mold_filter_context {
    const mold_molecule*                mol;
    uint64_t                            bit_count;
    uint64_t*                           bits;
    uint32_t                            sel_count;
    const mold_filter_stored_selection* sel;
    mold_error_callback*                err_cb;
};

bool mold_filter_parse(const char* expression, const mold_filter_context* context);

#ifdef __cplusplus
}
#endif

#endif
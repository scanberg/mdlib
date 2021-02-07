#ifndef _MD_CONTEXT_H_
#define _MD_CONTEXT_H_

#include <stdint.h>

struct md_allocator_i;
struct md_log_i;

struct md_context {
    struct md_allocator_i* alloc;
    struct md_allocator_i* temp_alloc;
    struct md_log_i* log;
};

#endif
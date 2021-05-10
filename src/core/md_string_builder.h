#ifndef __MD_STRING_BUILDER_H__
#define __MD_STRING_BUILDER_H__

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

typedef struct md_string_builder_t {
    struct md_string_builder_page_t {
        char buf[4096];
    };

    struct md_string_builder_page_t base_page;
    md_string_builder_t* pages;

    struct md_allocator_i* alloc;
} md_string_builder_t;

bool md_string_builder_init(md_string_builder_t* sb, struct md_allocator_i* alloc);
void md_string_builder_free(md_string_builder_t* sb);

void md_string_builder_append(md_string_builder_t* sb, const char* format, ...);
void md_string_builder_reset(md_string_builder_t* sb);

#ifdef __cplusplus
}
#endif

#endif

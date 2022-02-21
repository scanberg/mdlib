#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <core/md_str.h>

struct md_allocator_i;

typedef struct md_string_builder_t {
    char* buf;
    struct md_allocator_i* alloc;
} md_string_builder_t;

#ifdef __cplusplus
extern "C" {
#endif

void md_string_builder_init(md_string_builder_t* sb, struct md_allocator_i* alloc);
void md_string_builder_free(md_string_builder_t* sb);

void md_string_builder_print(md_string_builder_t* sb, const char* string);
void md_string_builder_printf(md_string_builder_t* sb, const char* format, ...);
void md_string_builder_append_str(md_string_builder_t* sb, str_t str);
void md_string_builder_reset(md_string_builder_t* sb);

str_t md_string_builder_to_string(md_string_builder_t* sb);

#ifdef __cplusplus
}
#endif

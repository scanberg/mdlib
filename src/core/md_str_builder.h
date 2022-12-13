#pragma once

#include <core/md_str.h>

struct md_allocator_i;

typedef struct md_str_builder_t {
    char* buf;
    struct md_allocator_i* alloc;
} md_str_builder_t;

#ifdef __cplusplus
extern "C" {
#endif

void md_str_builder_init(md_str_builder_t* sb, struct md_allocator_i* alloc);
void md_str_builder_free(md_str_builder_t* sb);

void md_str_builder_printf(md_str_builder_t* sb, const char* format, ...);
void md_str_builder_append_char(md_str_builder_t* sb, char c);
void md_str_builder_append_cstr(md_str_builder_t* sb, const char* cstr);
void md_str_builder_append_cstr_len(md_str_builder_t* sb, const char* cstr, int64_t len);
void md_str_builder_append_str(md_str_builder_t* sb, str_t str);
void md_str_builder_reset(md_str_builder_t* sb);

void md_str_builder_pop(md_str_builder_t* sb, int64_t n);

const char* md_str_builder_cstr(const md_str_builder_t* sb);
int64_t md_str_builder_len(const md_str_builder_t* sb);
str_t md_str_builder_to_str(const md_str_builder_t* sb);

#ifdef __cplusplus
}
#endif

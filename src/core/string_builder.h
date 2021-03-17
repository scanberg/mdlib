#ifndef __MD_STRING_BUILDER_H__
#define __MD_STRING_BUILDER_H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

typedef struct md_string_builder_o md_string_builder_o;
typedef struct md_string_builder_i {
    md_string_builder_o* inst;

    void (*append)(const char* str, uint64_t len);
    void (*print)(const char* format, ...);

    // Resets the internal buffers
    void (*reset)();
    
} md_string_builder_i;

md_string_builder_i* md_create_string_builder(struct md_allocator_i* alloc);
void md_destroy_string_builder(md_string_builder_i* sb);

#ifdef __cplusplus
}
#endif

#endif

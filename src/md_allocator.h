#ifndef _MD_ALLOCATOR_H_
#define _MD_ALLOCATOR_H_

#include <stdint.h>

#define md_alloc(a, sz) a->realloc(a, NULL, 0, sz, __FILE__, __LINE__)
#define md_realloc(a, ptr, old_sz, new_sz) a->realloc(a, ptr, old_sz, new_sz, __FILE__, __LINE__)
#define md_alloc_at(a, sz, file, line) a->realloc(a, NULL, 0, sz, file, line)
#define md_free(a, p, sz) a->realloc(a, p, sz, 0, __FILE__, __LINE__)

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_o;

struct md_allocator_i {
    struct md_allocator_o* inst; // Opaque data associated with allocator
    void *(*realloc)(struct md_allocator_i *a, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line);
};

extern struct md_allocator_i* default_allocator;
extern struct md_allocator_i* default_temp_allocator;

#ifdef __cplusplus
}
#endif

#endif
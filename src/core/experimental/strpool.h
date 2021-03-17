#ifndef _MD_STRPOOL_H_
#define _MD_STRPOOL_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

// This is a string pool implementation for interning of small strings
// which are assumed to fit into a single page.
// It is meant to be used to avoid duplication of string entries.
// Items can be inserted but never individually removed.

typedef struct md_strpool md_strpool;
typedef struct md_strpool_entry md_strpool_entry;

struct md_strpool_page;
struct md_strpool_entry {
    uint32_t offset : 22;   //  4M
    uint32_t length : 10;   //  1K
};

struct md_strpool {
    struct {
        struct md_strpool_page* ptr;
        uint32_t                  len;
    } page;

    struct {
        md_strpool_entry* ptr;
        uint32_t            len;
        uint32_t            cap;
    } entry;

    struct md_allocator_i* alloc;
};

void md_strpool_init(md_strpool* pool, struct md_allocator_i* alloc);
void md_strpool_free(md_strpool* pool);
void md_strpool_clear(md_strpool* pool);

md_strpool_entry md_strpool_insert(md_strpool* pool, const char* str, int len);
const char* md_strpool_cstr(md_strpool_entry entry);

#ifdef __cplusplus
}
#endif

#endif
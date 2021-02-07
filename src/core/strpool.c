#define _CRT_SECURE_NO_WARNINGS

#include "strpool.h"
#include <core/common.h>

#include <string.h> // memset
#include <stdlib.h> // malloc
#include <stdbool.h>

#include <md_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MD_STRPOOL_PAGE_SIZE 1024
#define MD_STRPOOL_ENTRY_CAP_INIT 16
#define MD_STRPOOL_ENTRY_CAP_GROW 16      // For the purpose of short strings (usually even less than 4 characters, we will not have alot of new unique entries. Therefore using linear growth makes sense

typedef struct md_strpool_page md_strpool_page;

struct md_strpool_page {
    char data[MD_STRPOOL_PAGE_SIZE];
};

internal inline bool realloc_pages(md_strpool* pool, uint32_t new_len) {
    ASSERT(pool);
    struct md_allocator_i* alloc = pool->alloc;
    if (!alloc) alloc = default_allocator;

    pool->page.ptr = md_realloc(alloc, pool->page.ptr, sizeof(md_strpool_page) * pool->page.len, sizeof(md_strpool_page) * new_len);
    pool->page.len = new_len;

    ASSERT(pool->page.ptr);
    return true;
}

internal inline bool realloc_entries(md_strpool* pool, uint32_t new_cap) {
    ASSERT(pool);
    struct md_allocator_i* alloc = pool->alloc;
    if (!alloc) alloc = default_allocator;

    pool->entry.ptr = md_realloc(alloc, pool->entry.ptr, sizeof(md_strpool_entry) * pool->entry.cap, sizeof(md_strpool_entry) * new_cap);
    pool->entry.cap = new_cap;

    ASSERT(pool->entry.ptr);
    return true;
}

internal inline const char* lookup_cstr(md_strpool* pool, md_strpool_entry entry) {
    unsigned int page_idx = entry.offset / MD_STRPOOL_PAGE_SIZE;
    unsigned int page_off = entry.offset & (MD_STRPOOL_PAGE_SIZE-1); // fast modulo
    return pool->page.ptr[page_idx].data + page_off;
}

internal inline md_strpool_entry* find_entry(md_strpool* pool, const char* str, uint32_t len) {
    ASSERT(pool);
    for (uint32_t i = 0; i < pool->entry.len; ++i) {
        md_strpool_entry* entry = &pool->entry.ptr[i];
        const char* entry_str = lookup_cstr(pool, *entry);
        if (entry->length == len && strncmp(str, entry_str, len) == 0) {
            return entry;
        }
    }
    return NULL;
}


void md_strpool_init(md_strpool* pool, struct md_allocator_i* alloc) {
    ASSERT(pool);
    memset(pool, 0, sizeof(*pool));
    md_strpool_clear(pool);
}

void md_strpool_free(md_strpool* pool) {
    ASSERT(pool);
    struct md_allocator_i* alloc = pool->alloc;
    if (!alloc) alloc = default_allocator;
    md_free(alloc, pool->page.ptr, pool->page.len * sizeof(md_strpool_page));
    md_free(alloc, pool->entry.ptr, pool->entry.cap * sizeof(md_strpool_entry));
}

void md_strpool_clear(md_strpool* pool) {
    ASSERT(pool);
    realloc_pages(pool, 1);
    realloc_entries(pool, MD_STRPOOL_ENTRY_CAP_INIT);
    pool->entry.len = 0;
}

md_strpool_entry md_strpool_insert(md_strpool* pool, const char* str, int len) {
    ASSERT(pool);
    ASSERT(len < MD_STRPOOL_PAGE_SIZE);

    md_strpool_entry* entry = find_entry(pool, str, len);
    if (entry) {
        return *entry;
    }

    unsigned int offset = pool->entry.len > 0 ? (pool->entry.ptr[pool->entry.len - 1].offset + pool->entry.ptr[pool->entry.len - 1].length + 1) : 0;
    unsigned int page_idx = offset / MD_STRPOOL_PAGE_SIZE;
    unsigned int page_off = offset & (MD_STRPOOL_PAGE_SIZE-1);
    unsigned int page_cap = (MD_STRPOOL_PAGE_SIZE - page_off);
    if ((len + 1) > page_cap) {
        // Need to allocate new page
        realloc_pages(pool, pool->page.len + 1);
        page_idx += 1;
        page_off = 0;
        page_cap = MD_STRPOOL_PAGE_SIZE;
        offset = page_idx * MD_STRPOOL_PAGE_SIZE;
    }

    strncpy(pool->page.ptr[page_idx].data + page_off, str, len);
    pool->page.ptr[page_idx].data[page_off + len] = '\0';

    if (pool->entry.len == pool->entry.cap) {
        realloc_entries(pool, pool->entry.cap + MD_STRPOOL_ENTRY_CAP_GROW);
    }
    pool->entry.ptr[pool->entry.len].offset = offset;
    pool->entry.ptr[pool->entry.len].length = len;
    pool->entry.len += 1;
}

const char* md_strpool_cstr(md_strpool* pool, md_strpool_entry entry) {
    ASSERT(pool);
    if (entry.offset > pool->entry.len) return NULL;
    unsigned int page_idx = entry.offset / MD_STRPOOL_PAGE_SIZE;
    unsigned int page_off = entry.offset & (MD_STRPOOL_PAGE_SIZE-1); // fast modulo for power of 2
    return pool->page.ptr[page_idx].data + page_off;
}

#ifdef __cplusplus
}
#endif
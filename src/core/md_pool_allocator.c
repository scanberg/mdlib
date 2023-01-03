#include "md_pool_allocator.h"

#include "md_allocator.h"
#include "md_common.h"
#include "md_array.h"
#include "md_intrinsics.h"

#include <stdbool.h>

#define MAGIC_NUMBER 0xf00b182c847bc69a
#define PAGE_SLOT_CAPACITY 64

typedef struct page_t {
    uint64_t free_slots; // This is a bitfield 64 slots, resort to intrinsics to find the next available bit
    void* mem;
} page_t;

// TODO:
// We need a free list for pages, where we insert pages that have one or more available slots.
typedef struct pool_t {
    uint64_t magic_number;
    uint64_t slot_size;
    page_t* pages;
    md_allocator_i* alloc;
} pool_t;

static inline uint64_t pool_page_size(const pool_t* pool) {
    return pool->slot_size * PAGE_SLOT_CAPACITY;
}

static inline void* page_slot_ptr(const page_t* page, uint64_t idx, uint64_t slot_size) {
    return (char*)page->mem + idx * slot_size;
}

static inline page_t* pool_new_page(pool_t* pool) {
    page_t page = {
        .free_slots = ~0ULL,
        .mem = md_alloc(pool->alloc, pool_page_size(pool))
    };
    ASSERT(page.mem);
    return md_array_push(pool->pages, page, pool->alloc);
}

static void* pool_new_slot(pool_t* pool, uint64_t size) {
    (void)size;
    ASSERT(size <= pool->slot_size);
    
    page_t* page = NULL;
    for (int64_t i = 0; i < md_array_size(pool->pages); ++i) {
        if (pool->pages[i].free_slots) {
            page = &pool->pages[i];
            break;
        }
    }
    if (!page) page = pool_new_page(pool);
    ASSERT(page);
    
    const uint64_t bit = bsf64(page->free_slots);
    ASSERT(bit);
    const uint64_t idx = bit - 1;

    page->free_slots &= ~(1ULL << idx); // Clear bit
    return page_slot_ptr(page, idx, pool->slot_size);
}

static void pool_free_slot(pool_t* pool, void* mem) {
    const uint64_t page_size = pool_page_size(pool);
    for (int64_t i = 0; i < md_array_size(pool->pages); ++i) {
        const char* base = (const char*)pool->pages[i].mem;
        if (base <= (char*)mem && (char*)mem < (base + page_size)) {
            uint64_t idx = ((uint64_t)mem - (uint64_t)pool->pages[i].mem) / pool->slot_size;
            pool->pages[i].free_slots |= (1ULL << idx);
            return;
        }
    }
    ASSERT(false); // failed to find page
}

static void pool_init(pool_t* pool, md_allocator_i* alloc, int64_t slot_size) {
    pool->magic_number = MAGIC_NUMBER;
    pool->slot_size = (uint64_t)slot_size;
    pool->pages = NULL;
    pool->alloc = alloc;
    pool_new_page(pool);
}

static void pool_free(pool_t* pool) {
    const uint64_t page_size = pool_page_size(pool);
    for (int64_t i = 0; i < md_array_size(pool->pages); ++i) {
        md_free(pool->alloc, pool->pages[i].mem, page_size);
    }
    md_array_free(pool->pages, pool->alloc);
}

static void* pool_realloc(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    (void)old_size;

    pool_t* pool = (pool_t*)inst;
    ASSERT(pool && pool->magic_number == MAGIC_NUMBER);
    if (new_size == 0) {
        pool_free_slot(pool, ptr);
        return NULL;
    }

    ASSERT(!ptr); // WE DO NOT DEAL PROPERLY WITH REALLOCATIONS, LET IT BE KNOWN!

    return pool_new_slot(pool, new_size);
}

struct md_allocator_i* md_pool_allocator_create(struct md_allocator_i* backing, uint64_t slot_size) {
    ASSERT(backing);
    ASSERT(slot_size > 0);
    uint64_t mem_size = sizeof(md_allocator_i) + sizeof(pool_t);
    void* mem = md_alloc(backing, mem_size);
    MEMSET(mem, 0, mem_size);

    md_allocator_i* pool_alloc = mem;
    md_allocator_o* inst = (md_allocator_o*)((char*)mem + sizeof(md_allocator_i));

    pool_t* pool = (pool_t*)inst;
    pool_init(pool, backing, slot_size);

    pool_alloc->inst = inst;
    pool_alloc->realloc = pool_realloc;

    return pool_alloc;
}

void md_pool_allocator_destroy(struct md_allocator_i* a) {
    pool_t* pool = (pool_t*)a->inst;
    ASSERT(pool->magic_number == MAGIC_NUMBER); // Make sure this allocator is a pool allocator
    pool_free(pool);
    md_free(pool->alloc, a, sizeof(md_allocator_i) + sizeof(pool_t));
}

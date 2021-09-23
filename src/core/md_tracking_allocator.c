#include "md_tracking_allocator.h"

#include <stdbool.h>

#include "md_allocator.h"
#include "md_array.inl"
#include "md_common.h"
#include "md_log.h"

#define MAGIC_NUMBER 0xdc1728367bca6273

typedef struct allocation {
    void* ptr;
    uint64_t size;
    const char* file;
    uint32_t line;
} allocation_t;

typedef struct tracking {
    struct md_allocator_i* backing;
    allocation_t* allocations;
    uint64_t magic;
} tracking_t;

allocation_t* find_allocation(tracking_t* tracking, void* ptr) {
    for (int64_t i = 0; i < md_array_size(tracking->allocations); ++i) {
        if (tracking->allocations[i].ptr == ptr) return &tracking->allocations[i];
    }
    return NULL;
}

allocation_t* new_allocation(tracking_t* tracking) {
    allocation_t item = {};
    return md_array_push(tracking->allocations, item, default_allocator);
}

void delete_allocation(tracking_t* tracking, allocation_t* alloc) {
    ASSERT(find_allocation(tracking, alloc->ptr));
    int64_t idx = alloc - tracking->allocations;
    *alloc = *md_array_last(tracking->allocations);
    tracking->allocations[idx] = *alloc;
    md_array_pop(tracking->allocations);
}

static void* tracking_realloc(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    tracking_t* tracking = (tracking_t*)inst;
    ASSERT(tracking && tracking->magic == MAGIC_NUMBER);

    if (ptr && old_size) {
        // REALLOC OR FREE
        // MAKE SURE OLD PTR AND SIZE ARE REPRESENTED!

        allocation_t* alloc = find_allocation(tracking, ptr);
        ASSERT(alloc);
        ASSERT(alloc->size == old_size);

        if (new_size) {
            // REALLOC
            alloc->ptr = tracking->backing->realloc(tracking->backing->inst, ptr, old_size, new_size, file, line);
            alloc->size = new_size;
            alloc->file = file;
            alloc->line = line;

            return alloc->ptr;
        }
        else {
            // FREE
            delete_allocation(tracking, alloc);
            return NULL;
        }
    } else {
        // MALLOC
        allocation_t* alloc = new_allocation(tracking);
        alloc->ptr = tracking->backing->realloc(tracking->backing->inst, ptr, old_size, new_size, file, line);
        alloc->size = new_size;
        alloc->file = file;
        alloc->line = line;

        return alloc->ptr;
    }
}

struct md_allocator_i* md_tracking_allocator_create(struct md_allocator_i* backing) {
    ASSERT(backing);
    tracking_t* inst = (tracking_t*)md_alloc(backing, sizeof(tracking_t) + sizeof(md_allocator_i));
    inst->backing = backing;
    inst->allocations = 0;
    inst->magic = MAGIC_NUMBER;

    md_allocator_i* alloc = (md_allocator_i*)((char*)inst + sizeof(tracking_t));
    alloc->inst = (md_allocator_o*)inst;
    alloc->realloc = tracking_realloc;

    return alloc;
}

void md_tracking_allocator_destroy(struct md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(alloc->inst);
    tracking_t* tracking = (tracking_t*)alloc->inst;
    ASSERT(tracking->magic == MAGIC_NUMBER);

    int64_t size = md_array_size(tracking->allocations);
    for (int64_t i = 0; i < size; ++i) {
        md_printf(MD_LOG_TYPE_DEBUG, "Allocation never freed, in file '%s', at line '%i'.", tracking->allocations[i].file, tracking->allocations[i].line);
    }

    md_free(tracking->backing, tracking, sizeof(tracking_t) + sizeof(md_allocator_i));
}


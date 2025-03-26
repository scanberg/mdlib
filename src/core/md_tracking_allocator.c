#include <core/md_tracking_allocator.h>

#include <core/md_os.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <stdbool.h>

// We need to protect this via a mutex lock since we don't know in which context this allocator is used.
// Better safe than sorry.

#define MAGIC_NUMBER 0xdc1728367bca6273

typedef struct {
    void* ptr;
    size_t size;
    const char* file;
    size_t line;
} allocation_t;

typedef struct {
    struct md_allocator_i* backing;
    allocation_t* allocations;
    md_mutex_t mutex;
    uint64_t magic;
} tracking_t;

allocation_t* find_allocation(tracking_t* tracking, void* ptr) {
    for (size_t i = 0; i < md_array_size(tracking->allocations); ++i) {
        if (tracking->allocations[i].ptr == ptr) return &tracking->allocations[i];
    }
    return NULL;
}

allocation_t* new_allocation(tracking_t* tracking) {
    allocation_t item = {0};
    md_array_push(tracking->allocations, item, md_get_heap_allocator());
    return md_array_last(tracking->allocations);
}

void register_deallocation(tracking_t* tracking, allocation_t* alloc) {
    ASSERT(find_allocation(tracking, alloc->ptr));
    size_t idx = alloc - tracking->allocations;
    md_array_swap_back_and_pop(tracking->allocations, idx);
    //*alloc = *md_array_last(tracking->allocations);
    //tracking->allocations[idx] = *alloc;
    //md_array_pop(tracking->allocations);
}

static void* tracking_realloc(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file;
    (void)line;
    tracking_t* tracking = (tracking_t*)inst;
    ASSERT(tracking && tracking->magic == MAGIC_NUMBER);

    void* result = 0;

    md_mutex_lock(&tracking->mutex);
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
            result = alloc->ptr;
        }
        else {
            // FREE
            register_deallocation(tracking, alloc);
            result = NULL;
        }
    } else if (new_size) {
        // MALLOC
        allocation_t* alloc = new_allocation(tracking);
        alloc->ptr = tracking->backing->realloc(tracking->backing->inst, ptr, old_size, new_size, file, line);
        alloc->size = new_size;
        alloc->file = file;
        alloc->line = line;
        result = alloc->ptr;
    }
    md_mutex_unlock(&tracking->mutex);
    return result;
}

struct md_allocator_i* md_tracking_allocator_create(struct md_allocator_i* backing) {
    ASSERT(backing);
    tracking_t* inst = (tracking_t*)md_alloc(backing, sizeof(tracking_t) + sizeof(md_allocator_i));
    inst->backing = backing;
    inst->allocations = 0;
    inst->mutex = md_mutex_create();
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

    md_mutex_lock(&tracking->mutex);
    size_t size = md_array_size(tracking->allocations);
    for (size_t i = 0; i < size; ++i) {
        MD_LOG_DEBUG("Allocation never freed, in file '%s', at line '%i'.", tracking->allocations[i].file, tracking->allocations[i].line);
    }
    md_mutex_unlock(&tracking->mutex);
    md_mutex_destroy(&tracking->mutex);

    md_free(tracking->backing, tracking, sizeof(tracking_t) + sizeof(md_allocator_i));
}

void md_tracking_allocator_print(struct md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(alloc->inst);
    tracking_t* tracking = (tracking_t*)alloc->inst;
    ASSERT(tracking->magic == MAGIC_NUMBER);

    md_mutex_lock(&tracking->mutex);

    MD_LOG_DEBUG("### Beg of Allocated Data ###");
    size_t size = md_array_size(tracking->allocations);
    for (size_t i = 0; i < size; ++i) {
        MD_LOG_DEBUG("'%s':%i", tracking->allocations[i].file, tracking->allocations[i].line);
    }
    MD_LOG_DEBUG("### End of Allocated Data ###");

    md_mutex_unlock(&tracking->mutex);
}


#include <core/md_allocator.h>

#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <stdbool.h>
#include <stdlib.h>

#ifndef MD_TEMP_ARENA_RESERVATION_SIZE
#ifdef MD_TEMP_ALLOC_SIZE
#define MD_TEMP_ARENA_RESERVATION_SIZE MD_TEMP_ALLOC_SIZE
#else
#define MD_TEMP_ARENA_RESERVATION_SIZE GIGABYTES(8)
#endif
#endif

#ifndef MD_TEMP_ALLOCATOR_COUNT
#define MD_TEMP_ALLOCATOR_COUNT 2
#endif

STATIC_ASSERT(MD_TEMP_ALLOCATOR_COUNT > 0, "At least one temporary allocator is required");

size_t md_temp_arena_reservation_size(void) {
    return MD_TEMP_ARENA_RESERVATION_SIZE;
}

THREAD_LOCAL md_allocator_i* _temp_allocators[MD_TEMP_ALLOCATOR_COUNT];
THREAD_LOCAL int _temp_exit_callback_registered;
THREAD_LOCAL md_allocator_i* _active_temp_allocator;
THREAD_LOCAL size_t _active_temp_depth;

static void release_temp_allocators(void* data) {
    (void)data;
    for (size_t i = 0; i < ARRAY_SIZE(_temp_allocators); ++i) {
        if (_temp_allocators[i]) {
            md_vm_arena_destroy(_temp_allocators[i]);
            _temp_allocators[i] = 0;
        }
    }
    _temp_exit_callback_registered = 0;
}

static void ensure_temp_exit_callback(void) {
    if (!_temp_exit_callback_registered) {
        bool registered = md_thread_on_exit(release_temp_allocators);
        ASSERT(registered);
        (void)registered;
        _temp_exit_callback_registered = 1;
    }
}

static void* realloc_internal(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)inst;
    (void)old_size;
    (void)file;
    (void)line;
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    return realloc(ptr, (size_t)new_size);
}

static struct md_allocator_i _heap_allocator = {
    NULL,
    realloc_internal,
};

static md_allocator_i* create_temp_allocator(void) {
    md_allocator_i* allocator = md_vm_arena_create(MD_TEMP_ARENA_RESERVATION_SIZE);
    ASSERT(allocator);
    return allocator;
}

static md_allocator_i* thread_temp_allocator_at(size_t index) {
    ASSERT(index < ARRAY_SIZE(_temp_allocators));
    ensure_temp_exit_callback();
    if (!_temp_allocators[index]) {
        _temp_allocators[index] = create_temp_allocator();
    }
    return _temp_allocators[index];
}

static bool allocator_matches(md_allocator_i* allocator, md_allocator_i* conflict) {
    return allocator && conflict && (allocator == conflict || allocator->inst == conflict->inst);
}

static bool allocator_is_conflicted(md_allocator_i* allocator, md_allocator_i* const* conflicts, size_t conflict_count) {
    if (!conflicts) {
        return false;
    }
    for (size_t i = 0; i < conflict_count; ++i) {
        if (allocator_matches(allocator, conflicts[i])) {
            return true;
        }
    }
    return false;
}

// Get general allocator interface to heap (malloc)
md_allocator_i* md_get_heap_allocator(void) {
    return &_heap_allocator;
}

md_allocator_i* md_get_temp_arena(void) {
    return thread_temp_allocator_at(0);
}

md_allocator_i* md_get_temp_arena_avoid(md_allocator_i* const* conflicts, size_t conflict_count) {
    for (size_t i = 0; i < ARRAY_SIZE(_temp_allocators); ++i) {
        md_allocator_i* allocator = thread_temp_allocator_at(i);
        if (!allocator_is_conflicted(allocator, conflicts, conflict_count)) {
            return allocator;
        }
    }

    MD_LOG_ERROR("Could not find a non-conflicting temporary allocator");
    ASSERT(false);
    return NULL;
}

static md_temp_t temp_begin(md_allocator_i* arena) {
    ASSERT(arena);
    md_temp_t temp = {
        .arena = arena,
        .pos = md_vm_arena_get_pos(arena),
        .prev_arena = _active_temp_allocator,
        .depth = _active_temp_depth + 1,
    };
    _active_temp_allocator = arena;
    _active_temp_depth = temp.depth;
    return temp;
}

md_temp_t md_temp_begin(void) {
    return temp_begin(md_get_temp_arena());
}

md_temp_t md_temp_begin_avoid(md_allocator_i* const* conflicts, size_t conflict_count) {
    return temp_begin(md_get_temp_arena_avoid(conflicts, conflict_count));
}

md_temp_t md_temp_begin_arena(md_allocator_i* arena) {
    return temp_begin(arena);
}

void md_temp_end(md_temp_t temp) {
    ASSERT(temp.arena);
    if (_active_temp_allocator != temp.arena || _active_temp_depth != temp.depth) {
        MD_LOG_ERROR("Temporary arena scopes must be ended in LIFO order");
        ASSERT(false);
        return;
    }

    md_vm_arena_set_pos_back(temp.arena, temp.pos);
    _active_temp_allocator = temp.prev_arena;
    _active_temp_depth = temp.depth - 1;
}

md_allocator_i* md_temp_allocator(md_temp_t temp) {
    ASSERT(temp.arena);
    return temp.arena;
}

static md_allocator_i* current_temp_allocator(void) {
    if (!_active_temp_allocator) {
        MD_LOG_ERROR("No active temporary arena scope");
        ASSERT(false);
        return NULL;
    }
    return _active_temp_allocator;
}

void* md_temp_push(size_t size) {
    return md_vm_arena_push(current_temp_allocator(), size);
}

void* md_temp_push_zero(size_t size) {
    return md_vm_arena_push_zero(current_temp_allocator(), size);
}

void* md_temp_push_aligned(size_t size, size_t alignment) {
    return md_vm_arena_push_aligned(current_temp_allocator(), size, alignment);
}

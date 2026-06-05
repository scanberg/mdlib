#include <core/md_allocator.h>

#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

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

static bool g_arena_key_initialized = false;
static md_thread_key_t g_arena_key;

typedef struct thread_local_data_t{
    md_allocator_i* allocator[MD_TEMP_ALLOCATOR_COUNT];
#if DEBUG
    md_temp_scope_t* scope_stack;
    size_t depth;
    size_t capacity;
#endif
} thread_local_data_t;

static THREAD_LOCAL thread_local_data_t* _thread_local_data = NULL;

static thread_local_data_t* thread_local_data_get(void);

static void release_thread_local_data(void* data) {
    thread_local_data_t* thread_data = (thread_local_data_t*)data;
    if (!thread_data) {
        return;
    }

    for (size_t i = 0; i < MD_TEMP_ALLOCATOR_COUNT; ++i) {
        if (thread_data->allocator[i]) {
            md_vm_arena_destroy(thread_data->allocator[i]);
            thread_data->allocator[i] = NULL;
        }
    }
#if DEBUG
    free(thread_data->scope_stack);
#endif

    if (_thread_local_data == thread_data) {
        _thread_local_data = NULL;
    }

    free(thread_data);
}

static bool ensure_arena_key_initialized(void) {
    if (!g_arena_key_initialized) {
        if (!md_thread_key_create(&g_arena_key, release_thread_local_data)) {
            MD_LOG_ERROR("Failed to create temporary arena thread local storage key");
            ASSERT(false);
            return false;
        }
        g_arena_key_initialized = true;
    }
    return true;
}

void md_temp_arena_system_init(void) {
    ensure_arena_key_initialized();
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

static inline bool allocator_matches(md_allocator_i* allocator, md_allocator_i* conflict) {
    return allocator && conflict && (allocator == conflict || allocator->inst == conflict->inst);
}

static struct md_allocator_i _heap_allocator = {
    NULL,
    realloc_internal,
};

static thread_local_data_t* thread_local_data_get(void) {
    if (_thread_local_data) {
        return _thread_local_data;
    }

    if (!ensure_arena_key_initialized()) {
        return NULL;
    }

    thread_local_data_t* thread_data = (thread_local_data_t*)calloc(1, sizeof(thread_local_data_t));
    ASSERT(thread_data);
    if (!thread_data) {
        return NULL;
    }

    md_thread_key_set_value(g_arena_key, thread_data);
    _thread_local_data = thread_data;
    return thread_data;
}

static bool ensure_scope_stack_capacity(thread_local_data_t* thread_data, size_t capacity) {
    ASSERT(thread_data);
#if DEBUG
    if (thread_data->capacity >= capacity) {
        return true;
    }

    size_t new_capacity = thread_data->capacity ? thread_data->capacity : 8;
    while (new_capacity < capacity) {
        new_capacity *= 2;
    }

    md_temp_scope_t* new_stack = (md_temp_scope_t*)realloc(thread_data->scope_stack, sizeof(md_temp_scope_t) * new_capacity);
    if (!new_stack) {
        MD_LOG_ERROR("Failed to allocate temporary scope stack");
        ASSERT(false);
        return false;
    }

    thread_data->scope_stack = new_stack;
    thread_data->capacity = new_capacity;
    return true;
#else
    (void)capacity;
    return true;
#endif
}

static md_allocator_i* create_temp_allocator(void) {
    md_allocator_i* allocator = md_vm_arena_create(MD_TEMP_ARENA_RESERVATION_SIZE);
    ASSERT(allocator);
    return allocator;
}

static md_allocator_i* thread_temp_allocator_at(size_t index) {
    ASSERT(index < MD_TEMP_ALLOCATOR_COUNT);

    thread_local_data_t* thread_data = thread_local_data_get();
    ASSERT(thread_data);
    if (!thread_data) {
        return NULL;
    }

    if (!thread_data->allocator[index]) {
        thread_data->allocator[index] = create_temp_allocator();
    }

    return thread_data->allocator[index];
}

static md_allocator_i* thread_temp_allocator_avoid(const md_allocator_i* conflict) {
    if (!conflict) {
        return thread_temp_allocator_at(0);
    }

    for (size_t i = 0; i < MD_TEMP_ALLOCATOR_COUNT; ++i) {
        md_allocator_i* allocator = thread_temp_allocator_at(i);
        if (!allocator_matches(allocator, (md_allocator_i*)conflict)) {
            return allocator;
        }
    }

    MD_LOG_ERROR("Could not find a non-conflicting temporary allocator");
    ASSERT(false);
    return NULL;
}

// Get general allocator interface to heap (malloc)
md_allocator_i* md_get_heap_allocator(void) {
    return &_heap_allocator;
}

static md_temp_scope_t temp_begin(md_allocator_i* arena) {
    ASSERT(arena);

#if DEBUG
    thread_local_data_t* thread_data = thread_local_data_get();
    ASSERT(thread_data);
    if (!thread_data || !ensure_scope_stack_capacity(thread_data, thread_data->depth + 1)) {
        return (md_temp_scope_t){0};
    }
#endif

    md_temp_scope_t temp = {
        .arena = arena,
        .pos = md_vm_arena_get_pos(arena),
    };

#if DEBUG
    thread_data->scope_stack[thread_data->depth++] = temp;
#endif
    return temp;
}

md_temp_scope_t md_temp_begin(void) {
    return temp_begin(thread_temp_allocator_at(0));
}

md_temp_scope_t md_temp_begin_avoid(const md_allocator_i* avoid) {
    return temp_begin(thread_temp_allocator_avoid(avoid));
}

md_temp_scope_t md_temp_begin_in(md_allocator_i* arena) {
    return temp_begin(arena);
}

void md_temp_end(md_temp_scope_t temp) {
    ASSERT(temp.arena);

#if DEBUG
    thread_local_data_t* thread_data = thread_local_data_get();
    ASSERT(thread_data);
    if (!thread_data || thread_data->depth == 0) {
        MD_LOG_ERROR("No active temporary arena scope");
        ASSERT(false);
        return;
    }

    md_temp_scope_t active_scope = thread_data->scope_stack[thread_data->depth - 1];
    if (active_scope.arena != temp.arena || active_scope.pos != temp.pos) {
        MD_LOG_ERROR("Temporary arena scopes must be ended in LIFO order");
        ASSERT(false);
        return;
    }

    thread_data->depth -= 1;
#endif
    md_vm_arena_set_pos_back(temp.arena, temp.pos);
}

md_allocator_i* md_temp_allocator(md_temp_scope_t temp) {
    ASSERT(temp.arena);
    return temp.arena;
}

void* md_temp_alloc(md_temp_scope_t temp, size_t size) {
    ASSERT(temp.arena);
    return md_vm_arena_push(temp.arena, size);
}

void* md_temp_alloc_zero(md_temp_scope_t temp, size_t size) {
    ASSERT(temp.arena);
    return md_vm_arena_push_zero(temp.arena, size);
}

void* md_temp_alloc_aligned(md_temp_scope_t temp, size_t size, size_t alignment) {
    ASSERT(temp.arena);
    return md_vm_arena_push_aligned(temp.arena, size, alignment);
}

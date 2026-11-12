#include <core/md_allocator.h>
#include <core/md_format.h>

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
#define MD_TEMP_ARENA_RESERVATION_SIZE MD_GIGABYTES(8)
#endif
#endif

#ifndef MD_TEMP_ALLOCATOR_COUNT
#define MD_TEMP_ALLOCATOR_COUNT 2
#endif

MD_STATIC_ASSERT(MD_TEMP_ALLOCATOR_COUNT > 0, "At least one temporary allocator is required");

size_t md_temp_arena_reservation_size(void) {
    return MD_TEMP_ARENA_RESERVATION_SIZE;
}

static bool g_arena_key_initialized = false;
static md_thread_key_t g_arena_key;

typedef struct thread_local_data_t{
    md_allocator_i* allocator[MD_TEMP_ALLOCATOR_COUNT];
#if MD_DEBUG
    md_temp_scope_t* scope_stack;
    size_t depth;
    size_t capacity;
#endif
} thread_local_data_t;

static MD_THREAD_LOCAL thread_local_data_t* _thread_local_data = NULL;

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
#if MD_DEBUG
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
            MD_ASSERT(false);
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
    MD_ASSERT(thread_data);
    if (!thread_data) {
        return NULL;
    }

    md_thread_key_set_value(g_arena_key, thread_data);
    _thread_local_data = thread_data;
    return thread_data;
}

static bool ensure_scope_stack_capacity(thread_local_data_t* thread_data, size_t capacity) {
    MD_ASSERT(thread_data);
#if MD_DEBUG
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
        MD_ASSERT(false);
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
    MD_ASSERT(allocator);
    return allocator;
}

static md_allocator_i* thread_temp_allocator_at(size_t index) {
    MD_ASSERT(index < MD_TEMP_ALLOCATOR_COUNT);

    thread_local_data_t* thread_data = thread_local_data_get();
    MD_ASSERT(thread_data);
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
    MD_ASSERT(false);
    return NULL;
}

// Get general allocator interface to heap (malloc)
md_allocator_i* md_get_heap_allocator(void) {
    return &_heap_allocator;
}

static md_temp_scope_t temp_begin(md_allocator_i* arena) {
    MD_ASSERT(arena);

#if MD_DEBUG
    thread_local_data_t* thread_data = thread_local_data_get();
    MD_ASSERT(thread_data);
    if (!thread_data || !ensure_scope_stack_capacity(thread_data, thread_data->depth + 1)) {
        return (md_temp_scope_t){0};
    }
#endif

    md_temp_scope_t temp = {
        .arena = arena,
        .pos = md_vm_arena_get_pos(arena),
    };

#if MD_DEBUG
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
    MD_ASSERT(temp.arena);

#if MD_DEBUG
    thread_local_data_t* thread_data = thread_local_data_get();
    MD_ASSERT(thread_data);
    if (!thread_data || thread_data->depth == 0) {
        MD_LOG_ERROR("No active temporary arena scope");
        MD_ASSERT(false);
        return;
    }

    md_temp_scope_t active_scope = thread_data->scope_stack[thread_data->depth - 1];
    if (active_scope.arena != temp.arena || active_scope.pos != temp.pos) {
        MD_LOG_ERROR("Temporary arena scopes must be ended in LIFO order");
        MD_ASSERT(false);
        return;
    }

    thread_data->depth -= 1;
#endif
    md_vm_arena_set_pos_back(temp.arena, temp.pos);
}

md_temp_mark_t md_temp_mark(md_temp_scope_t temp) {
    MD_ASSERT(temp.arena);
    md_temp_mark_t mark = { md_vm_arena_get_pos(temp.arena) };
    return mark;
}

void md_temp_reset(md_temp_scope_t temp, md_temp_mark_t mark) {
    MD_ASSERT(temp.arena);
    if (mark.pos < temp.pos) {
        MD_LOG_ERROR("A temporary scope is reset to a mark before its own beginning");
        MD_ASSERT(false);
        return;
    }
#if MD_DEBUG
    thread_local_data_t* thread_data = thread_local_data_get();
    if (thread_data && thread_data->depth > 0) {
        const md_temp_scope_t active = thread_data->scope_stack[thread_data->depth - 1];
        if (active.arena != temp.arena || active.pos != temp.pos) {
            MD_LOG_ERROR("Only the innermost temporary scope is reset to a mark: the scopes begun after it end first");
            MD_ASSERT(false);
            return;
        }
    }
#endif
    if (mark.pos <= md_vm_arena_get_pos(temp.arena)) {
        md_vm_arena_set_pos_back(temp.arena, mark.pos);
    }
}

void md_temp_trim(void) {
    thread_local_data_t* thread_data = _thread_local_data;
    if (!thread_data) return;
    for (size_t i = 0; i < MD_TEMP_ALLOCATOR_COUNT; ++i) {
        if (thread_data->allocator[i]) md_vm_arena_trim(thread_data->allocator[i]);
    }
}

md_allocator_i* md_temp_allocator(md_temp_scope_t temp) {
    MD_ASSERT(temp.arena);
    return temp.arena;
}

void* md_temp_alloc(md_temp_scope_t temp, size_t size) {
    MD_ASSERT(temp.arena);
    return md_vm_arena_push(temp.arena, size);
}

void* md_temp_alloc_zero(md_temp_scope_t temp, size_t size) {
    MD_ASSERT(temp.arena);
    return md_vm_arena_push_zero(temp.arena, size);
}

void* md_temp_alloc_aligned(md_temp_scope_t temp, size_t size, size_t alignment) {
    MD_ASSERT(temp.arena);
    return md_vm_arena_push_aligned(temp.arena, size, alignment);
}

// ### STRINGS ###

md_str_t md_alloc_str(md_allocator_i* alloc, md_str_t str) {
    MD_ASSERT(alloc);
    md_str_t result = {0};
    if (str.ptr && str.len > 0) {
        char* data = md_alloc(alloc, str.len + 1);
        md_memcpy(data, str.ptr, str.len);
        data[str.len] = '\0';
        result.ptr = data;
        result.len = str.len;
    }
    return result;
}

static md_str_t md_alloc_vprintf(md_allocator_i* alloc, const char* format, va_list args) {
    MD_ASSERT(alloc);
    MD_ASSERT(format);
    va_list measure;
    va_copy(measure, args);
    const int len = md_vsnprintf(NULL, 0, format, measure);
    va_end(measure);

    if (len <= 0) {
        return (md_str_t){0};
    }

    char* buf = md_alloc(alloc, (size_t)len + 1);
    md_vsnprintf(buf, (size_t)len + 1, format, args);
    return (md_str_t){buf, (size_t)len};
}

md_str_t md_alloc_printf(md_allocator_i* alloc, const char* format, ...) {
    va_list args;
    va_start(args, format);
    const md_str_t result = md_alloc_vprintf(alloc, format, args);
    va_end(args);
    return result;
}

void md_free_str(md_allocator_i* alloc, md_str_t str) {
    MD_ASSERT(alloc);
    if (str.ptr) {
        md_free(alloc, (void*)str.ptr, str.len + 1);
    }
}

#include <core/md_arena_allocator.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_os.h>

#define ARENA_MAGIC 0xfdc1728d827856cb
#define DEFAULT_ALIGNMENT (sizeof(void*)*2) // Should be 16 when compiled for x64

#define VM_MAGIC 0x87b716a78ccb2813
#define VM_COMMIT_SIZE MEGABYTES(1)
#define VM_DECOMMIT_THRESHOLD MEGABYTES(1)

#define VALIDATE_AND_EXTRACT_ARENA(alloc) \
    ASSERT(alloc); \
    ASSERT(alloc->inst); \
    arena_t* arena = (arena_t*)alloc->inst; \
    ASSERT(arena->magic == ARENA_MAGIC);

typedef struct page_t {
    void* mem;
    size_t size;
    size_t curr; // current offset
    struct page_t* next;
} page_t;

typedef struct arena_t {
    struct md_allocator_i* backing;
    page_t *base_page;
    page_t *curr_page;
    size_t default_page_size;
    uint64_t magic;
} arena_t;

static inline void arena_reset(arena_t* arena) {
    page_t* p = arena->base_page;
    while (p) {
        page_t* next = p->next;
        md_free(arena->backing, p, sizeof(page_t) + p->size); // @NOTE: we free the page + its memory
        p = next;
    }
    arena->base_page = NULL;
    arena->curr_page = NULL;
}

static inline page_t* arena_new_page(arena_t* arena, size_t size) {
    // @NOTE: Allocate page + memory
    page_t* page = md_alloc(arena->backing, sizeof(page_t) + size);
    page->mem = (char*)page + sizeof(page_t);
    page->size = size;
    page->curr = 0;
    page->next = NULL;
    
    if (arena->curr_page) {
        arena->curr_page->next = page;
    }
    arena->curr_page = page;
    
    if (!arena->base_page) {
        arena->base_page = page;
    }

    return page;
}

static inline void* arena_push(arena_t* arena, size_t size, size_t alignment) {
    page_t* p = arena->curr_page;
    if (!p || p->curr + size + alignment > p->size) {
        // We add the alignment here to be conservative in that the new page will have room for the size
        p = arena_new_page(arena, MAX(size + alignment, arena->default_page_size));
    }

    const size_t addr = (size_t)p->mem + p->curr;
    const size_t mod = addr & (alignment - 1);      // modulo: addr % alignment, but fast for power of 2
    const size_t aligned_curr = mod ? p->curr + alignment - mod : p->curr;

    p->curr = aligned_curr + size;
    return (char*)p->mem + aligned_curr;
}

static inline void arena_pop(arena_t* arena, size_t size) {
    // We cheat a bit, we only shrink within the current page and not beyond that.
    page_t* p = arena->curr_page;
    p->curr = size > p->curr ? 0 : p->curr - size;
}

static void* arena_realloc(struct md_allocator_o *inst, void *ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file;
    (void)line;
    arena_t* arena = (arena_t*)inst;
    ASSERT(arena && arena->magic == ARENA_MAGIC);

    // Realloc or Free
    if (ptr && old_size) {
        if (new_size == 0) {
            // Free
            return NULL;
        }
        if ((char*)ptr + old_size == (char*)arena->curr_page->mem + arena->curr_page->curr) {
            // This is the last allocation that occured, then we can shrink or grow that sucker
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            int64_t new_curr = (int64_t)arena->curr_page->curr + diff;
            ASSERT(new_curr >= 0);
            if (new_curr < (int64_t)arena->curr_page->size) {
                arena->curr_page->curr = new_curr;
                return ptr;
            }
        }
        size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
        // ptr is not the last allocation or the new size did not fit into the existing page.
        void* new_ptr = arena_push(arena, new_size, alignment);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }
    // Alloc
    size_t alignment = new_size <= 2 ? new_size : DEFAULT_ALIGNMENT;
    return arena_push(arena, new_size, alignment);
}

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, size_t page_size) {
    ASSERT(backing);
    arena_t* arena = (arena_t*)md_alloc(backing, sizeof(arena_t) + sizeof(md_allocator_i));
    arena->backing = backing;
    arena->base_page = NULL;
    arena->curr_page = NULL;
    arena->default_page_size = page_size;
    arena->magic = ARENA_MAGIC;

    md_allocator_i* arena_alloc = (md_allocator_i*)((char*)arena + sizeof(arena_t));
    arena_alloc->inst = (md_allocator_o*)arena;
    arena_alloc->realloc = arena_realloc;

    return arena_alloc;
}

void md_arena_allocator_reset(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    arena_reset(arena);
}

md_allocator_i* md_arena_allocator_backing(md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    return arena->backing;
}

void md_arena_allocator_destroy(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    arena_reset(arena);
    md_free(arena->backing, arena, sizeof(arena_t) + sizeof(md_allocator_i));
}

void* md_arena_allocator_push(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    size_t alignment = size <= 2 ? size : DEFAULT_ALIGNMENT;
    return arena_push(arena, size, alignment);
}

void* md_arena_allocator_push_zero(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    size_t alignment = size <= 2 ? size : DEFAULT_ALIGNMENT;
    void* mem = arena_push(arena, size, alignment);
    MEMSET(mem, 0, size);
    return mem;
}

void* md_arena_allocator_push_aligned(struct md_allocator_i* alloc, size_t size, size_t alignment) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    ASSERT(IS_POW2(alignment));
    return arena_push(arena, size, alignment);
}

void md_arena_allocator_pop(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_ARENA(alloc);
    arena_pop(arena, size);
}


// VM
#define VALIDATE_AND_EXTRACT_VM_ARENA(alloc) \
    ASSERT(alloc); \
    ASSERT(alloc->inst); \
    vm_arena_t* arena = (vm_arena_t*)alloc->inst; \
    ASSERT(arena->magic == VM_MAGIC);
#define MIN_VM_POS sizeof(vm_arena_t) + sizeof(md_allocator_i)

typedef struct vm_arena_t {
    void* base;   // Base adress of reserved memory space
    size_t size;  // Reservation size
    size_t pos;
    size_t commit_pos;
    size_t align;
    uint64_t magic;
} vm_arena_t;

static inline void* vm_push(vm_arena_t* arena, size_t size, size_t alignment) {
    ASSERT(IS_POW2(alignment));

    void* mem = 0;
    alignment = MAX(arena->align, alignment);

    size_t pos = arena->pos;
    size_t pos_address = (size_t)arena->base + pos;
    size_t alignment_size = ALIGN_TO(pos_address, alignment) - pos_address;

    if (pos + alignment_size + size <= arena->size) {
        mem = (char*)arena->base + pos + alignment_size;
        size_t new_pos = pos + alignment_size + size;
        arena->pos = new_pos;

        if (new_pos > arena->commit_pos) {
            size_t commit_size = ALIGN_TO(new_pos - arena->commit_pos, VM_COMMIT_SIZE);
            md_vm_commit((char*)arena->base + arena->commit_pos, commit_size);
            arena->commit_pos += commit_size;
        }
    }

    return mem;
}

static inline void vm_set_pos(vm_arena_t* arena, size_t pos) {
    ASSERT(pos <= arena->pos);
    arena->pos = pos;
#if 0
    size_t decommit_pos = ALIGN_TO(pos, VM_COMMIT_SIZE);
    size_t over_commited = arena->commit_pos - decommit_pos;
    if (decommit_pos > 0 && over_commited >= VM_DECOMMIT_THRESHOLD) {
        md_vm_decommit((char*)arena->base + decommit_pos, over_commited);
        arena->commit_pos -= over_commited;
    }
#endif
}

static inline void vm_pop(vm_arena_t* arena, size_t size) {
    // Ensure that the new position is not within the meta information of the allocator
    // That would be very bad...
    int64_t new_pos = arena->pos - size;
    ASSERT(new_pos >= MIN_VM_POS);
    vm_set_pos(arena, (size_t)new_pos);
}

static void* vm_arena_realloc(struct md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file;
    (void)line;
    vm_arena_t* arena = (vm_arena_t*)inst;
    ASSERT(arena);
    ASSERT(arena->magic == VM_MAGIC);

    if (new_size == 0) {
        if ((char*)arena->base + arena->pos == (char*)ptr + old_size) {
            vm_pop(arena, old_size);
        }
        return NULL;
    }

    if (old_size) {
        ASSERT(ptr);
        if ((char*)arena->base + arena->pos == (char*)ptr + old_size) {
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            const int64_t new_pos = arena->pos + diff;
            ASSERT(0 <= new_pos && new_pos < (int64_t)arena->size);
            if (diff < 0) {
                vm_set_pos(arena, new_pos);
            }
            else {
                vm_push(arena, (size_t)diff, arena->align);
            }
            return ptr;
        }
        void* new_ptr = vm_push(arena, new_size, arena->align);
        MEMCPY(new_ptr, ptr, old_size);
        return new_ptr;
    }

    return vm_push(arena, new_size, arena->align);
}

struct md_allocator_i* md_vm_arena_create(size_t reservation_size) {
    reservation_size = ALIGN_TO(reservation_size, MEGABYTES(1));
    void* base = md_vm_reserve(reservation_size);
    md_vm_commit(base, VM_COMMIT_SIZE);
    vm_arena_t* arena = base;
    arena->base = base;
    arena->size = reservation_size;
    arena->pos = MIN_VM_POS;
    arena->commit_pos = VM_COMMIT_SIZE;
    arena->align = DEFAULT_ALIGNMENT;
    arena->magic = VM_MAGIC;
    md_allocator_i* alloc = (md_allocator_i*)((char*)base + sizeof(vm_arena_t));
    alloc->inst = (md_allocator_o*)arena;
    alloc->realloc = vm_arena_realloc;

    return alloc;
}

void md_vm_arena_destroy(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    md_vm_release(arena->base, arena->size);
}

void* md_vm_arena_push_aligned(struct md_allocator_i* alloc, size_t size, size_t align) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    ASSERT(IS_POW2(align));

    void* mem = 0;
    align = MAX(arena->align, align);

    size_t pos = arena->pos;
    size_t pos_address = (size_t)arena->base + pos;
    size_t alignment_size = ALIGN_TO(pos_address, align) - pos_address;

    if (pos + alignment_size + size <= arena->size) {
        mem = (char*)arena->base + pos + alignment_size;
        size_t new_pos = pos + alignment_size + size;
        arena->pos = new_pos;

        if (new_pos > arena->commit_pos) {
            size_t commit_size = ALIGN_TO(new_pos - arena->commit_pos, VM_COMMIT_SIZE);
            md_vm_commit((char*)arena->base + arena->commit_pos, commit_size);
            arena->commit_pos += commit_size;
        }
    }

    return mem;
}

void* md_vm_arena_push(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    return vm_push(arena, size, arena->align);
}

void* md_vm_arena_push_zero(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    void* mem = vm_push(arena, size, arena->align);
    MEMSET(mem, 0, size);
    return mem;
}

void md_vm_arena_pop(struct md_allocator_i* alloc, size_t size) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    vm_pop(arena, size);
}

void md_vm_arena_reset(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    vm_set_pos(arena, sizeof(vm_arena_t) + sizeof(md_allocator_i));
}

void md_vm_arena_set_pos_back(struct md_allocator_i* alloc, size_t pos) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    vm_set_pos(arena, pos);
}

size_t md_vm_arena_get_pos(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    return arena->pos;
}

md_vm_arena_temp_t md_vm_arena_temp_begin(struct md_allocator_i* alloc) {
    VALIDATE_AND_EXTRACT_VM_ARENA(alloc)
    return (md_vm_arena_temp_t) {alloc, arena->pos};
}

void md_vm_arena_temp_end(md_vm_arena_temp_t temp) {
    md_vm_arena_set_pos_back(temp.vm_arena, temp.pos);
}

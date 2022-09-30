#include "md_arena_allocator.h"

#include "md_allocator.h"
#include "md_array.inl"
#include "md_os.h"

#define MAGIC_NUMBER 0xfdc1728d827856cb
#define DEFAULT_ALIGNMENT (sizeof(void*)*2) // Should be 16 when compiled for x64

#define VM_MAGIC 0x87b716a78ccb2813
#define VM_COMMIT_SIZE KILOBYTES(64)
#define VM_DECOMMIT_THRESHOLD KILOBYTES(64)

typedef struct page_t {
    void* mem;
    uint64_t size;
    uint64_t curr; // current offset
    struct page_t* next;
} page_t;

typedef struct arena_t {
    struct md_allocator_i* alloc;
    page_t *base_page;
    page_t *curr_page;
    uint64_t default_page_size;
    uint64_t magic;
} arena_t;

static inline void arena_reset(arena_t* arena) {
    page_t* p = arena->base_page;
    while (p) {
        page_t* next = p->next;
        md_free(arena->alloc, p, sizeof(page_t) + p->size); // @NOTE: we free the page + its memory
        p = next;
    }
    arena->base_page = NULL;
    arena->curr_page = NULL;
}

static inline page_t* arena_new_page(arena_t* arena, uint64_t size) {
    // @NOTE: Allocate page + memory
    page_t* page = md_alloc(arena->alloc, sizeof(page_t) + size);
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

static inline void* arena_alloc(arena_t* arena, uint64_t size) {
    // We want to make sure that we maintain some type of alignment for allocations
    // By default on x64, the alignment should be 16 bytes when using malloc or stack allocations (alloca)
    // For single byte and double byte allocations, we could skip this default alignment, but for 4 byte and above
    // we should enforce this as x64 calling convention uses xmm registers for passing floats, which could result in a
    // performance penalty if the data is not 16 byte aligned.

    const uint64_t alignment = (size <= 2) ? size : DEFAULT_ALIGNMENT;

    page_t* p = arena->curr_page;
    if (!p || p->curr + size + alignment > p->size) {
        // We add the alignment here to be conservative in that the new page will have room for the size
        p = arena_new_page(arena, MAX(size + alignment, arena->default_page_size));
    }

    const uint64_t addr = (uint64_t)p->mem + p->curr;
    const uint64_t mod = addr & (alignment - 1);      // modulo: addr % alignment, but fast for power of 2
    const uint64_t aligned_curr = mod ? p->curr + alignment - mod : p->curr;

    p->curr = aligned_curr + size;
    return (char*)p->mem + aligned_curr;
}

static void* arena_realloc(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    arena_t* arena = (arena_t*)inst;
    ASSERT(arena && arena->magic == MAGIC_NUMBER);

    // Free
    if (new_size == 0) {
        // We cannot assert here since new_size may be NULL in case of a new allocation of size 0 as well.
        //ASSERT(ptr);
        //ASSERT(old_size);
        return NULL;
    }
    // Realloc
    if (ptr && old_size) {
        if ((char*)ptr + old_size == (char*)arena->curr_page->mem + arena->curr_page->curr) {
            // This is the last allocation that occured, then we can shrink or grow that sucker
            const int64_t diff = (int64_t)new_size - (int64_t)old_size;
            int64_t new_curr = (int64_t)arena->curr_page->curr + diff;
            ASSERT(new_curr > 0);
            if (new_curr < (int64_t)arena->curr_page->size) {
                arena->curr_page->curr = new_curr;
                return ptr;
            }
        }
        // ptr is not the last allocation or the new size did not fit into the existing page.
        void* new_ptr = arena_alloc(arena, new_size);
        memcpy(new_ptr, ptr, old_size);
        return new_ptr;
    }
    // Alloc
    return arena_alloc(arena, new_size);
}

struct md_allocator_i* md_arena_allocator_create(struct md_allocator_i* backing, uint64_t page_size) {
    ASSERT(backing);
    arena_t* arena = (arena_t*)md_alloc(backing, sizeof(arena_t) + sizeof(md_allocator_i));
    arena->alloc = backing;
    arena->base_page = NULL;
    arena->curr_page = NULL;
    arena->default_page_size = page_size;
    arena->magic = MAGIC_NUMBER;

    md_allocator_i* arena_alloc = (md_allocator_i*)((char*)arena + sizeof(arena_t));
    arena_alloc->inst = (md_allocator_o*)arena;
    arena_alloc->realloc = arena_realloc;

    return arena_alloc;
}

void md_arena_allocator_reset(struct md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(alloc->inst);
    arena_t* arena = (arena_t*)alloc->inst;
    ASSERT(arena->magic == MAGIC_NUMBER);
    arena_reset(arena);
}

void md_arena_allocator_destroy(struct md_allocator_i* alloc) {
    ASSERT(alloc);
    ASSERT(alloc->inst);
    arena_t* arena = (arena_t*)alloc->inst;
    ASSERT(arena->magic == MAGIC_NUMBER);
    arena_reset(arena);
    md_free(alloc, arena, sizeof(arena_t) + sizeof(md_allocator_i));
}




// VM

void md_vm_arena_init(md_vm_arena_t* arena, uint64_t reservation_size) {
    reservation_size = ALIGN_TO(reservation_size, GIGABYTES(1));
    arena->base = md_os_reserve(reservation_size);
    arena->size = reservation_size;
    arena->commit_pos = VM_COMMIT_SIZE;
    arena->pos = 0;
    arena->align = DEFAULT_ALIGNMENT;
    arena->magic = VM_MAGIC;
    md_os_commit(arena->base, VM_COMMIT_SIZE);
}

void md_vm_arena_free(md_vm_arena_t* arena) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    md_os_release(arena->base);
    memset(arena, 0, sizeof(md_vm_arena_t));
}

void* md_vm_arena_push_aligned(md_vm_arena_t* arena, uint64_t size, uint64_t align) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    ASSERT(IS_POW2(align));

    void* mem = 0;
    align = MAX(arena->align, align);

    uint64_t pos = arena->pos;
    uint64_t pos_address = (uint64_t)arena->base + pos;
    uint64_t alignment_size = ALIGN_TO(pos_address, align) - pos_address;

    if (pos + alignment_size + size <= arena->size) {
        mem = (char*)arena->base + pos + alignment_size;
        uint64_t new_pos = pos + alignment_size + size;
        arena->pos = new_pos;

        if (new_pos > arena->commit_pos) {
            uint64_t commit_size = ALIGN_TO(new_pos - arena->commit_pos, VM_COMMIT_SIZE);
            md_os_commit((char*)arena->base + arena->commit_pos, commit_size);
            arena->commit_pos += commit_size;
        }
    }

    return mem;
}

void* md_vm_arena_push(md_vm_arena_t* arena, uint64_t size) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    return md_vm_arena_push_aligned(arena, size, arena->align);
}

void* md_vm_arena_push_zero(md_vm_arena_t* arena, uint64_t size) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    void* mem = md_vm_arena_push(arena, size);
    memset(mem, 0, size);
    return mem;
}

void md_vm_arena_pop(md_vm_arena_t* arena, uint64_t size) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    md_vm_arena_set_pos(arena, arena->pos - size);
}

void md_vm_arena_reset(md_vm_arena_t* arena) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    md_vm_arena_pop(arena, arena->pos);
}

void md_vm_arena_set_pos(md_vm_arena_t* arena, uint64_t pos) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    ASSERT(pos < arena->pos);
    arena->pos = pos;
    uint64_t decommit_pos = ALIGN_TO(pos, VM_COMMIT_SIZE);
    uint64_t over_commited = arena->commit_pos - decommit_pos;
    if (decommit_pos > 0 && over_commited >= VM_DECOMMIT_THRESHOLD) {
        md_os_decommit((char*)arena->base, over_commited);
        arena->commit_pos -= over_commited;
    }
}

uint64_t md_vm_arena_get_pos(md_vm_arena_t* arena) {
    return arena->pos;
}

md_vm_arena_temp_t md_vm_arena_temp_begin(md_vm_arena_t* arena) {
    ASSERT(arena && arena->magic == VM_MAGIC);
    return (md_vm_arena_temp_t) {arena, arena->pos};
}

void md_vm_arena_temp_end(md_vm_arena_temp_t temp) {
    md_vm_arena_set_pos(temp.arena, temp.pos);
}


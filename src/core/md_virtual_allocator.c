#include "md_virtual_allocator.h"

#include <core/md_os.h>
#include <string.h>

#define MAGIC 0xfdc1274a

static void* vm_realloc(struct md_allocator_o *inst, void *ptr, uint64_t old_size, uint64_t new_size, const char* file, uint32_t line) {
    (void)file;
    (void)line;
    md_vm_allocator_t* vm = (md_vm_allocator_t*)inst;
    ASSERT(vm && vm->magic == MAGIC);

    uint64_t old_page_aligned_size = ROUND_UP(old_size, vm->page_size);
    uint64_t new_page_aligned_size = ROUND_UP(new_size, vm->page_size);

    // FREE
    if (new_size == 0) {
        // Only decommit full pages
        md_os_decommit(ptr, old_page_aligned_size);

        // If the allocation was the last one commited, we move the position back
        if ((char*)vm->base + vm->commit_pos == (char*)ptr + old_page_aligned_size) {
            vm->commit_pos -= old_page_aligned_size;
        }
    }

    // REALLOC
    if (ptr && old_size) {
        if (old_page_aligned_size == new_page_aligned_size) {
            // The reallocation fit into the already commited page space
            return ptr;
        }

        if (new_page_aligned_size > old_page_aligned_size) {
            // GROW
            bool last = ((char*)vm->base + vm->commit_pos == (char*)ptr + old_page_aligned_size);
            void* new_ptr = last ? ptr : (char*)vm->base + vm->commit_pos;
            uint64_t commit_size = last ? new_page_aligned_size - old_page_aligned_size : new_page_aligned_size;
            md_os_commit((char*)vm->base + vm->commit_pos, commit_size);
            vm->commit_pos += commit_size;
            if (!last) {
                memcpy(new_ptr, ptr, old_size);
            }
            return new_ptr;
        } else {
            uint64_t decommit_size = old_page_aligned_size - new_page_aligned_size;
            md_os_decommit((char*)vm->base + new_page_aligned_size, decommit_size);
            // If the allocation was the last one commited, we move the position back
            if ((char*)vm->base + vm->commit_pos == (char*)ptr + old_page_aligned_size) {
                vm->commit_pos -= decommit_size;
            }
        }
    }

    // ALLOC
    void* new_ptr = (char*)vm->base + vm->commit_pos;
    uint64_t commit_size = ROUND_UP(new_size, vm->page_size);
    md_os_commit((char*)vm->base + vm->commit_pos, commit_size);
    vm->commit_pos += commit_size;
    return new_ptr;
}

void md_vm_allocator_init(md_vm_allocator_t* vm, uint64_t reservation_size) {
    void* ptr = md_os_reserve(reservation_size);
    uint32_t page_size = (uint32_t)md_os_page_size();
    md_os_commit(ptr, page_size);

    vm->base = ptr;
    vm->commit_pos = page_size;
    vm->size = reservation_size;
    vm->page_size = (uint32_t)md_os_page_size();
    vm->magic = MAGIC;
}

void md_vm_allocator_free(md_vm_allocator_t* vm) {
    ASSERT(vm && vm->magic == MAGIC);
    md_os_release(vm->base);
    vm->base = 0;
    vm->size = 0;
    vm->commit_pos = 0;
}

void* md_vm_allocator_base(md_vm_allocator_t* vm) {
    ASSERT(vm && vm->magic == MAGIC);
    return vm->base;
}

uint64_t md_vm_allocator_size(md_vm_allocator_t* vm) {
    ASSERT(vm && vm->magic == MAGIC);
    return vm->size;
}

md_allocator_i md_vm_allocator_create_interface(md_vm_allocator_t* vm) {
    md_allocator_i alloc = {
        .inst = (md_allocator_o*)vm,
        .realloc = vm_realloc,
    };
    return alloc;
}

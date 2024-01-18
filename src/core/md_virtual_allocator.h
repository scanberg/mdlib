#pragma once

#include <stdint.h>
#include <core/md_allocator.h>

typedef struct md_vm_allocator_t md_vm_allocator_t;
struct md_vm_allocator_t {
    void* base;             // Reservation base address
    size_t size;          // Reservation size
    size_t commit_pos;    // Commited pos in bytes (always aligned to page size)
    size_t page_size;
    uint64_t magic;
};

#ifdef __cplusplus
extern "C" {
#endif

// This is the virtual memory allocator which is exposed through Windows (VirtualAlloc) and Unix (mmap / mprotect)
// In principle, you reserve a range of the the user memory adress space (8 TB on Windows?) and commit to chunks within that space that you can use.
// When commiting a chunk of the memory space, only then is it backed by physical memory.
// Each allocation is rounded up to the page size of the OS, so this is not for small allocations
// Note that each allocation results in a systems call and clearing of pages, so this is essentially a building block for other allocators
// The reservation size is going to be rounded up to the nearest GIGABYTE (It makes little sense to allocate small portions within such vast address space)

void md_vm_allocator_init(md_vm_allocator_t* vm, size_t reservation_size);
void md_vm_allocator_free(md_vm_allocator_t* vm);

// Returns the base adress of the reserved memory space
void* md_vm_allocator_base(md_vm_allocator_t* vm);

// Returns the size of the reserved memory space
size_t md_vm_allocator_size(md_vm_allocator_t* vm);

// Create generic interface for vm allocator
md_allocator_i md_vm_allocator_create_interface(md_vm_allocator_t* vm);

#ifdef __cplusplus
}
#endif

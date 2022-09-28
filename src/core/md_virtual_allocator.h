#pragma once

#include "md_allocator.h"

#define MD_VM_ALLOCATOR_RESERVATION_SIZE (4 * 1024ULL * 1024ULL)

typedef struct md_vm_allocator_t md_vm_allocator_t;
struct md_vm_allocator_t {
    void* base;             // Reservation base address
    uint64_t size;          // Reservation size
    uint64_t commit_pos;    // Commited pos in bytes (always aligned to page size)
    uint32_t page_size;
    uint32_t magic;
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

void md_vm_allocator_init(md_vm_allocator_t* vm, uint64_t reservation_size);
void md_vm_allocator_free(md_vm_allocator_t* vm);

// Returns the base adress of the reserved memory space
void* md_vm_allocator_base(md_vm_allocator_t* vm);

// Returns the size of the reserved memory space
uint64_t md_vm_allocator_size(md_vm_allocator_t* vm);

// Create generic interface for vm allocator
md_allocator_i md_vm_allocator_create_interface(md_vm_allocator_t* vm);

#ifdef __cplusplus
}
#endif

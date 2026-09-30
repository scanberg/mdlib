/*
md_gpu_tlsf.h -- internal to the md_gpu backends.

A two-level segregated-fit (TLSF) sub-allocator over offset ranges. It knows
nothing about GPUs: a "region" is an opaque pointer the caller associates with
a range [0, size), in practice one device buffer. Bookkeeping lives entirely
in host memory, because the memory being managed may not be host-visible.

Allocation and free are O(1): a size maps to one of 64 x 32 free lists, and a
pair of bitmaps finds the first non-empty list that is guaranteed to fit. With
32 second-level lists the rounding a request suffers when searching is under
1/32 of its size. Freed nodes merge with free physical neighbours at once, so
a region whose allocations have all been freed becomes one free node again,
which is how the caller learns it may give the region back.

All sizes and offsets are multiples of `granularity` (a power of two), which
is therefore also the alignment of every allocation.

Not thread-safe; the backends call it under their device lock.
*/

#ifndef MD_GPU_TLSF_H
#define MD_GPU_TLSF_H

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;

#define MD_TLSF_SL_LOG2  5
#define MD_TLSF_SL_COUNT (1u << MD_TLSF_SL_LOG2)
#define MD_TLSF_FL_COUNT 64

typedef struct md_tlsf_node_t {
    uint64_t offset;                    /* within its region                    */
    uint64_t size;
    void*    region;
    struct md_tlsf_node_t* prev_phys;   /* neighbours in the region, or NULL    */
    struct md_tlsf_node_t* next_phys;
    struct md_tlsf_node_t* prev_free;   /* free-list links while free           */
    struct md_tlsf_node_t* next_free;
    bool     is_free;
} md_tlsf_node_t;

typedef struct md_tlsf_slab_t md_tlsf_slab_t;

typedef struct md_tlsf_t {
    struct md_allocator_i* alloc;
    uint64_t        granularity;
    uint64_t        fl_bitmap;
    uint32_t        sl_bitmap[MD_TLSF_FL_COUNT];
    md_tlsf_node_t* heads[MD_TLSF_FL_COUNT][MD_TLSF_SL_COUNT];
    md_tlsf_node_t* spare;              /* recycled node structs               */
    md_tlsf_slab_t* slabs;              /* node storage, freed on destroy      */
    uint64_t        free_bytes;
} md_tlsf_t;

void md_tlsf_init(md_tlsf_t* t, struct md_allocator_i* alloc, uint64_t granularity);

/* Frees the bookkeeping. Regions are the caller's to release. */
void md_tlsf_destroy(md_tlsf_t* t);

/* Adds [0, size) of `region` as one free node and returns it. `size` is
   rounded down to the granularity. NULL on host out-of-memory. */
md_tlsf_node_t* md_tlsf_add_region(md_tlsf_t* t, void* region, uint64_t size);

/* A node of at least `size` bytes (rounded up to the granularity), or NULL
   when nothing fits. The remainder of a larger free node is split off and
   stays free. */
md_tlsf_node_t* md_tlsf_alloc(md_tlsf_t* t, uint64_t size);

/* Frees `node`, merging it with free neighbours. Returns the resulting free
   node, which may be a neighbour that absorbed it; `node` must not be used
   afterwards. */
md_tlsf_node_t* md_tlsf_free(md_tlsf_t* t, md_tlsf_node_t* node);

/* True for a free node covering its whole region: the region is empty. */
static inline bool md_tlsf_region_empty(const md_tlsf_node_t* free_node) {
    return free_node && free_node->is_free && !free_node->prev_phys && !free_node->next_phys;
}

/* Removes an empty region (given its single free node) from the allocator. */
void md_tlsf_remove_region(md_tlsf_t* t, md_tlsf_node_t* free_node);

#endif /* MD_GPU_TLSF_H */

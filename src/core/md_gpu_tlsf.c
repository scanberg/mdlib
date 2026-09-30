#include "md_gpu_tlsf.h"

#include <core/md_allocator.h>
#include <core/md_intrinsics.h>

#include <string.h>

#define MD_TLSF_SLAB_NODES 255

struct md_tlsf_slab_t {
    md_tlsf_slab_t* next;
    md_tlsf_node_t  nodes[MD_TLSF_SLAB_NODES];
};

static inline uint32_t md_tlsf_msb(uint64_t v) {
    return 63u - (uint32_t)clz64(v);
}

/* The (fl, sl) list a free node of `size` lives in. */
static void md_tlsf_mapping(uint64_t size, uint32_t* fl, uint32_t* sl) {
    if (size < MD_TLSF_SL_COUNT) {
        *fl = 0;
        *sl = (uint32_t)size;
    } else {
        const uint32_t m = md_tlsf_msb(size);
        *fl = m;
        *sl = (uint32_t)(size >> (m - MD_TLSF_SL_LOG2)) - MD_TLSF_SL_COUNT;
    }
}

static md_tlsf_node_t* md_tlsf_node_new(md_tlsf_t* t) {
    if (!t->spare) {
        md_tlsf_slab_t* slab = (md_tlsf_slab_t*)md_alloc(t->alloc, sizeof(md_tlsf_slab_t));
        if (!slab) return NULL;
        slab->next = t->slabs;
        t->slabs   = slab;
        for (int i = MD_TLSF_SLAB_NODES - 1; i >= 0; --i) {
            slab->nodes[i].next_free = t->spare;
            t->spare = &slab->nodes[i];
        }
    }
    md_tlsf_node_t* n = t->spare;
    t->spare = n->next_free;
    memset(n, 0, sizeof(*n));
    return n;
}

static void md_tlsf_node_release(md_tlsf_t* t, md_tlsf_node_t* n) {
    n->next_free = t->spare;
    t->spare     = n;
}

static void md_tlsf_insert_free(md_tlsf_t* t, md_tlsf_node_t* n) {
    uint32_t fl, sl;
    md_tlsf_mapping(n->size, &fl, &sl);
    n->is_free   = true;
    n->prev_free = NULL;
    n->next_free = t->heads[fl][sl];
    if (n->next_free) n->next_free->prev_free = n;
    t->heads[fl][sl] = n;
    t->fl_bitmap    |= 1ull << fl;
    t->sl_bitmap[fl] |= 1u << sl;
    t->free_bytes   += n->size;
}

static void md_tlsf_remove_free(md_tlsf_t* t, md_tlsf_node_t* n) {
    uint32_t fl, sl;
    md_tlsf_mapping(n->size, &fl, &sl);
    if (n->prev_free) n->prev_free->next_free = n->next_free;
    else              t->heads[fl][sl]        = n->next_free;
    if (n->next_free) n->next_free->prev_free = n->prev_free;
    if (!t->heads[fl][sl]) {
        t->sl_bitmap[fl] &= ~(1u << sl);
        if (!t->sl_bitmap[fl]) t->fl_bitmap &= ~(1ull << fl);
    }
    n->is_free   = false;
    n->prev_free = NULL;
    n->next_free = NULL;
    t->free_bytes -= n->size;
}

void md_tlsf_init(md_tlsf_t* t, struct md_allocator_i* alloc, uint64_t granularity) {
    memset(t, 0, sizeof(*t));
    t->alloc       = alloc;
    t->granularity = granularity ? granularity : 1;
}

void md_tlsf_destroy(md_tlsf_t* t) {
    md_tlsf_slab_t* s = t->slabs;
    while (s) {
        md_tlsf_slab_t* next = s->next;
        md_free(t->alloc, s, sizeof(md_tlsf_slab_t));
        s = next;
    }
    struct md_allocator_i* alloc = t->alloc;
    uint64_t gran = t->granularity;
    md_tlsf_init(t, alloc, gran);
}

md_tlsf_node_t* md_tlsf_add_region(md_tlsf_t* t, void* region, uint64_t size) {
    size &= ~(t->granularity - 1);
    if (size == 0) return NULL;
    md_tlsf_node_t* n = md_tlsf_node_new(t);
    if (!n) return NULL;
    n->offset = 0;
    n->size   = size;
    n->region = region;
    md_tlsf_insert_free(t, n);
    return n;
}

/* First node in the class that is certain to fit `size` (good fit), or, when
   every such class is empty, a node of exactly the right class that happens
   to be big enough. The second step matters for a region sized exactly for
   one request, which the rounded-up search would otherwise skip. */
static md_tlsf_node_t* md_tlsf_find(md_tlsf_t* t, uint64_t size) {
    uint64_t search = size;
    if (search >= MD_TLSF_SL_COUNT) {
        search += (1ull << (md_tlsf_msb(search) - MD_TLSF_SL_LOG2)) - 1;
    }
    uint32_t fl, sl;
    md_tlsf_mapping(search, &fl, &sl);
    if (fl < MD_TLSF_FL_COUNT) {
        uint32_t sl_map = t->sl_bitmap[fl] & (~0u << sl);
        if (!sl_map) {
            const uint64_t fl_map = (fl + 1 < MD_TLSF_FL_COUNT) ? (t->fl_bitmap & (~0ull << (fl + 1))) : 0;
            if (fl_map) {
                fl     = (uint32_t)ctz64(fl_map);
                sl_map = t->sl_bitmap[fl];
            }
        }
        if (sl_map) return t->heads[fl][(uint32_t)ctz32(sl_map)];
    }
    md_tlsf_mapping(size, &fl, &sl);
    for (md_tlsf_node_t* n = t->heads[fl][sl]; n; n = n->next_free) {
        if (n->size >= size) return n;
    }
    return NULL;
}

md_tlsf_node_t* md_tlsf_alloc(md_tlsf_t* t, uint64_t size) {
    if (size == 0) size = 1;
    const uint64_t g = t->granularity;
    if (size > UINT64_MAX - g) return NULL;
    size = (size + g - 1) & ~(g - 1);

    md_tlsf_node_t* n = md_tlsf_find(t, size);
    if (!n) return NULL;
    md_tlsf_remove_free(t, n);

    if (n->size - size >= g) {
        md_tlsf_node_t* rest = md_tlsf_node_new(t);
        if (rest) {
            rest->offset    = n->offset + size;
            rest->size      = n->size - size;
            rest->region    = n->region;
            rest->prev_phys = n;
            rest->next_phys = n->next_phys;
            if (rest->next_phys) rest->next_phys->prev_phys = rest;
            n->next_phys = rest;
            n->size      = size;
            md_tlsf_insert_free(t, rest);
        }
        /* Without a node for the remainder the whole node is handed out --
           wasteful but correct. */
    }
    return n;
}

md_tlsf_node_t* md_tlsf_free(md_tlsf_t* t, md_tlsf_node_t* n) {
    md_tlsf_node_t* prev = n->prev_phys;
    if (prev && prev->is_free) {
        md_tlsf_remove_free(t, prev);
        prev->size     += n->size;
        prev->next_phys = n->next_phys;
        if (prev->next_phys) prev->next_phys->prev_phys = prev;
        md_tlsf_node_release(t, n);
        n = prev;
    }
    md_tlsf_node_t* next = n->next_phys;
    if (next && next->is_free) {
        md_tlsf_remove_free(t, next);
        n->size     += next->size;
        n->next_phys = next->next_phys;
        if (n->next_phys) n->next_phys->prev_phys = n;
        md_tlsf_node_release(t, next);
    }
    md_tlsf_insert_free(t, n);
    return n;
}

void md_tlsf_remove_region(md_tlsf_t* t, md_tlsf_node_t* n) {
    md_tlsf_remove_free(t, n);
    md_tlsf_node_release(t, n);
}

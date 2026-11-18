/*
md_gpu_heap.c -- the memory behind md_gpu_malloc, md_gpu_free and the temp
scopes, shared by the md_gpu backends. See md_gpu_heap.h.

md_gpu_malloc draws from one heap per memory kind: a TLSF sub-allocator over
chunks the backend makes. Chunks start at MD_GPU_HEAP_CHUNK_MIN and double up
to MD_GPU_HEAP_CHUNK_MAX; a request bigger than that gets a chunk of its own.
So a thousand small buffers cost a handful of driver allocations, not a
thousand -- a driver's allocation limit (Vulkan's maxMemoryAllocationCount,
commonly 4096) is never in reach.

md_gpu_free is stream-ordered. At a point the GPU has already passed, the node
goes straight back to the heap. Otherwise it waits on pending_frees until the
stream gets there (checked at malloc and poll). While it waits, the freeing
stream itself may reuse it for DEVICE memory: stream order puts the new work
after the old, and the backend's reuse_early covers EXPLICIT mode.
Host-visible memory is never reused early, because the CPU would write it at
once.

A chunk whose last allocation is released counts as empty. Empty chunks are
kept, up to cache_limit bytes per kind, and beyond that returned to the
backend. That is safe at once, because a node is released only after the GPU
has passed its free.

Debug memory gives every allocation a chunk of its own, exactly its size, and
never reuses it; the backend makes such chunks however lets its validation
see past their ends.
*/

#include "md_gpu_heap.h"

#include <core/md_allocator.h>
#include <core/md_format.h>
#include <core/md_os.h>

#include <stdarg.h>
#include <string.h>

/* A registry entry: one live md_gpu_malloc allocation, or one whole temp
   chunk. What md_gpu_heap_resolve checks addresses against. */
typedef struct md_gpu_heap_range_t {
    uint64_t          address;
    uint64_t          size;          /* requested bytes, or the chunk's size    */
    md_gpu_chunk_t*   chunk;
    uint64_t          chunk_offset;  /* of `address` within the chunk           */
    md_tlsf_node_t*   node;          /* heap allocation; NULL for a temp chunk  */
    md_gpu_chunk_t*   own;           /* debug memory: the allocation's own chunk */
    md_gpu_mem_kind_t kind;
} md_gpu_heap_range_t;

/* An allocation freed at a point its stream has not yet reached. */
typedef struct md_gpu_heap_pending_free_t {
    md_tlsf_node_t*   node;          /* a heap allocation, or               */
    md_gpu_chunk_t*   own;           /* a debug memory allocation's chunk   */
    md_gpu_mem_kind_t kind;
    md_gpu_stream_t   stream;
    uint64_t          value;
} md_gpu_heap_pending_free_t;

static inline uint64_t md_gpu_heap_align_up(uint64_t v, uint64_t a) {
    return (v + a - 1) & ~(a - 1);
}

static bool md_gpu_heap_fail(const md_gpu_heap_t* h, const char* fmt, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    md_vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    h->backend.fail(buf);
    return false;
}

/* ---- Vectors ------------------------------------------------------------------- */

static void md_gpu_heap_vec_init(md_gpu_heap_vec_t* v, size_t stride) {
    v->data = NULL; v->count = 0; v->capacity = 0; v->stride = stride;
}

static bool md_gpu_heap_vec_reserve(md_gpu_heap_vec_t* v, struct md_allocator_i* alloc, size_t n) {
    if (n <= v->capacity) return true;
    size_t cap = v->capacity ? v->capacity * 2 : 16;
    while (cap < n) cap *= 2;
    void* mem = md_alloc(alloc, cap * v->stride);
    if (!mem) return false;
    if (v->data) {
        memcpy(mem, v->data, v->count * v->stride);
        md_free(alloc, v->data, v->capacity * v->stride);
    }
    v->data = mem;
    v->capacity = cap;
    return true;
}

static void* md_gpu_heap_vec_push(md_gpu_heap_vec_t* v, struct md_allocator_i* alloc) {
    if (!md_gpu_heap_vec_reserve(v, alloc, v->count + 1)) return NULL;
    void* slot = (char*)v->data + v->count * v->stride;
    memset(slot, 0, v->stride);
    v->count++;
    return slot;
}

/* Remove element i, preserving order. */
static void md_gpu_heap_vec_remove(md_gpu_heap_vec_t* v, size_t i) {
    char* base = (char*)v->data;
    memmove(base + i * v->stride, base + (i + 1) * v->stride, (v->count - i - 1) * v->stride);
    v->count--;
}

/* Remove the first element equal to the pointer `p` (vectors of pointers). */
static void md_gpu_heap_vec_remove_ptr(md_gpu_heap_vec_t* v, const void* p) {
    void** arr = (void**)v->data;
    for (size_t i = 0; i < v->count; ++i) {
        if (arr[i] == p) { md_gpu_heap_vec_remove(v, i); return; }
    }
}

static void md_gpu_heap_vec_free(md_gpu_heap_vec_t* v, struct md_allocator_i* alloc) {
    if (v->data) md_free(alloc, v->data, v->capacity * v->stride);
    v->data = NULL; v->count = 0; v->capacity = 0;
}

#define MD_GPU_HEAP_VEC_AT(v, type, i) (((type*)(v).data)[i])

/* ---- Registry ------------------------------------------------------------------ */

/* Binary search for the range containing `address`. */
static md_gpu_heap_range_t* md_gpu_heap_registry_find_locked(const md_gpu_heap_t* h, uint64_t address) {
    size_t lo = 0, hi = h->registry.count;
    md_gpu_heap_range_t** arr = (md_gpu_heap_range_t**)h->registry.data;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        const md_gpu_heap_range_t* r = arr[mid];
        if (address < r->address)                hi = mid;
        else if (address >= r->address + r->size) lo = mid + 1;
        else                                      return arr[mid];
    }
    return NULL;
}

/* Where `address` goes, or is, in the registry. */
static size_t md_gpu_heap_registry_lower_bound(const md_gpu_heap_t* h, uint64_t address) {
    size_t lo = 0, hi = h->registry.count;
    md_gpu_heap_range_t** arr = (md_gpu_heap_range_t**)h->registry.data;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (arr[mid]->address < address) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static bool md_gpu_heap_registry_insert_locked(md_gpu_heap_t* h, md_gpu_heap_range_t* r) {
    if (!md_gpu_heap_vec_reserve(&h->registry, h->alloc, h->registry.count + 1)) return false;
    md_gpu_heap_range_t** arr = (md_gpu_heap_range_t**)h->registry.data;
    const size_t at = md_gpu_heap_registry_lower_bound(h, r->address);
    memmove(arr + at + 1, arr + at, (h->registry.count - at) * sizeof(*arr));
    arr[at] = r;
    h->registry.count++;
    return true;
}

static void md_gpu_heap_registry_remove_locked(md_gpu_heap_t* h, md_gpu_heap_range_t* r) {
    const size_t at = md_gpu_heap_registry_lower_bound(h, r->address);
    if (at < h->registry.count && MD_GPU_HEAP_VEC_AT(h->registry, md_gpu_heap_range_t*, at) == r) {
        md_gpu_heap_vec_remove(&h->registry, at);
    }
}

/* The registry is mutated by malloc on other threads, so the lookup takes the
   lock; what it found is copied out before the lock is let go. */
bool md_gpu_heap_resolve(md_gpu_heap_t* h, md_gpu_addr_t addr, uint64_t size, md_gpu_heap_span_t* out, const char* what) {
    md_mutex_lock(h->mutex);
    const md_gpu_heap_range_t* r = md_gpu_heap_registry_find_locked(h, addr);
    md_gpu_heap_range_t copy;
    if (r) copy = *r;
    md_mutex_unlock(h->mutex);
    if (!r) return md_gpu_heap_fail(h, "%s: 0x%llx is not a live md_gpu allocation", what, (unsigned long long)addr);
    const uint64_t off = addr - copy.address;
    if (size > copy.size || off > copy.size - size) {
        return md_gpu_heap_fail(h, "%s: range [0x%llx, +%llu) overruns its %llu-byte allocation",
                                what, (unsigned long long)addr, (unsigned long long)size, (unsigned long long)copy.size);
    }
    if (out) {
        out->chunk  = copy.chunk;
        out->offset = copy.chunk_offset + off;
        out->host   = copy.chunk->host ? copy.chunk->host + copy.chunk_offset + off : NULL;
    }
    return true;
}

/* ---- Heaps ---------------------------------------------------------------------- */

void md_gpu_heap_init(md_gpu_heap_t* h, struct md_allocator_i* alloc, struct md_mutex_t* mutex,
                      const md_gpu_heap_backend_t* backend, uint64_t cache_limit, bool debug) {
    memset(h, 0, sizeof(*h));
    h->alloc       = alloc;
    h->mutex       = mutex;
    h->backend     = *backend;
    h->cache_limit = cache_limit ? cache_limit : MD_GPU_HEAP_CACHE_DEFAULT;
    h->debug       = debug;
    md_gpu_heap_vec_init(&h->registry, sizeof(md_gpu_heap_range_t*));
    md_gpu_heap_vec_init(&h->pending_frees, sizeof(md_gpu_heap_pending_free_t));
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_tlsf_init(&h->kinds[k].tlsf, alloc, MD_GPU_HEAP_ALIGN);
        md_gpu_heap_vec_init(&h->kinds[k].chunks, sizeof(md_gpu_chunk_t*));
    }
}

void md_gpu_heap_destroy(md_gpu_heap_t* h) {
    md_mutex_lock(h->mutex);
    for (size_t i = 0; i < h->registry.count; ++i) {
        md_free(h->alloc, MD_GPU_HEAP_VEC_AT(h->registry, md_gpu_heap_range_t*, i), sizeof(md_gpu_heap_range_t));
    }
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_gpu_heap_kind_t* hk = &h->kinds[k];
        for (size_t i = 0; i < hk->chunks.count; ++i) {
            h->backend.chunk_destroy(h->backend.user, MD_GPU_HEAP_VEC_AT(hk->chunks, md_gpu_chunk_t*, i));
        }
        md_gpu_heap_vec_free(&hk->chunks, h->alloc);
        md_tlsf_destroy(&hk->tlsf);
    }
    md_gpu_heap_vec_free(&h->registry, h->alloc);
    md_gpu_heap_vec_free(&h->pending_frees, h->alloc);
    md_mutex_unlock(h->mutex);
}

static uint64_t md_gpu_heap_chunk_size(const md_gpu_heap_kind_t* hk, uint64_t need) {
    uint64_t size = MD_GPU_HEAP_CHUNK_MIN;
    for (size_t i = 0; i < hk->chunks.count && size < MD_GPU_HEAP_CHUNK_MAX; ++i) size *= 2;
    if (need > size) size = md_gpu_heap_align_up(need, MD_GPU_HEAP_LARGE_ALIGN);
    return size;
}

/* Give a heap node back to its heap: the GPU is done with it. */
static void md_gpu_heap_release_node_locked(md_gpu_heap_t* h, md_gpu_mem_kind_t kind, md_tlsf_node_t* node) {
    md_gpu_heap_kind_t* hk = &h->kinds[kind];
    md_tlsf_node_t* f = md_tlsf_free(&hk->tlsf, node);
    if (!md_tlsf_region_empty(f)) return;
    md_gpu_chunk_t* c = (md_gpu_chunk_t*)f->region;
    c->empty      = true;
    c->empty_node = f;
    hk->empty_bytes += c->size;
    if (hk->empty_bytes > h->cache_limit) {
        md_tlsf_remove_region(&hk->tlsf, f);
        md_gpu_heap_vec_remove_ptr(&hk->chunks, c);
        hk->empty_bytes -= c->size;
        hk->reserved    -= c->size;
        h->backend.chunk_destroy(h->backend.user, c);
    }
}

/* Give a debug memory allocation's chunk back: the GPU is done with it. */
static void md_gpu_heap_own_release_locked(md_gpu_heap_t* h, md_gpu_mem_kind_t kind, md_gpu_chunk_t* c) {
    md_gpu_heap_kind_t* hk = &h->kinds[kind];
    md_gpu_heap_vec_remove_ptr(&hk->chunks, c);
    hk->reserved -= c->size;
    h->backend.chunk_destroy(h->backend.user, c);
}

static void md_gpu_heap_release_free_locked(md_gpu_heap_t* h, const md_gpu_heap_pending_free_t* pf) {
    if (pf->own) md_gpu_heap_own_release_locked(h, pf->kind, pf->own);
    else         md_gpu_heap_release_node_locked(h, pf->kind, pf->node);
}

void md_gpu_heap_poll_locked(md_gpu_heap_t* h) {
    for (size_t i = 0; i < h->pending_frees.count;) {
        const md_gpu_heap_pending_free_t pf = MD_GPU_HEAP_VEC_AT(h->pending_frees, md_gpu_heap_pending_free_t, i);
        if (h->backend.stream_completed(pf.stream) >= pf.value) {
            md_gpu_heap_vec_remove(&h->pending_frees, i);
            md_gpu_heap_release_free_locked(h, &pf);
        } else {
            ++i;
        }
    }
}

void md_gpu_heap_forget_stream_locked(md_gpu_heap_t* h, md_gpu_stream_t s) {
    for (size_t i = 0; i < h->pending_frees.count;) {
        const md_gpu_heap_pending_free_t pf = MD_GPU_HEAP_VEC_AT(h->pending_frees, md_gpu_heap_pending_free_t, i);
        if (pf.stream == s) {
            md_gpu_heap_vec_remove(&h->pending_frees, i);
            md_gpu_heap_release_free_locked(h, &pf);
        } else {
            ++i;
        }
    }
}

/* A node for `need` bytes: the stream's own pending frees first (DEVICE only,
   and only when little of the node would be wasted), then the heap, then a new
   chunk. */
static md_tlsf_node_t* md_gpu_heap_alloc_locked(md_gpu_heap_t* h, md_gpu_stream_t s,
                                                md_gpu_mem_kind_t kind, uint64_t need) {
    md_gpu_heap_kind_t* hk = &h->kinds[kind];
    md_tlsf_node_t* node = NULL;

    if (kind == MD_GPU_MEM_DEVICE) {
        for (size_t i = 0; i < h->pending_frees.count; ++i) {
            const md_gpu_heap_pending_free_t* pf = &MD_GPU_HEAP_VEC_AT(h->pending_frees, md_gpu_heap_pending_free_t, i);
            if (pf->stream != s || pf->kind != kind || pf->own) continue;
            if (pf->node->size < need || pf->node->size > need + need / 4 + MD_GPU_HEAP_ALIGN) continue;
            if (!h->backend.reuse_early(s)) break;
            node = pf->node;
            md_gpu_heap_vec_remove(&h->pending_frees, i);
            break;
        }
    }
    if (!node) node = md_tlsf_alloc(&hk->tlsf, need);
    if (!node) {
        const uint64_t size = md_gpu_heap_chunk_size(hk, need);
        md_gpu_chunk_t* c = h->backend.chunk_create(h->backend.user, kind, size);
        if (!c) return NULL;
        md_gpu_chunk_t** slot = (md_gpu_chunk_t**)md_gpu_heap_vec_push(&hk->chunks, h->alloc);
        if (!slot || !md_tlsf_add_region(&hk->tlsf, c, size)) {
            if (slot) hk->chunks.count--;
            h->backend.chunk_destroy(h->backend.user, c);
            md_gpu_heap_fail(h, "out of memory");
            return NULL;
        }
        *slot = c;
        hk->reserved += size;
        node = md_tlsf_alloc(&hk->tlsf, need);
        if (!node) { md_gpu_heap_fail(h, "md_gpu_malloc: internal error, a fresh chunk did not fit"); return NULL; }
    }
    md_gpu_chunk_t* c = (md_gpu_chunk_t*)node->region;
    if (c->empty) {
        c->empty      = false;
        c->empty_node = NULL;
        hk->empty_bytes -= c->size;
    }
    return node;
}

md_gpu_mem_t md_gpu_heap_malloc(md_gpu_heap_t* h, md_gpu_stream_t s, md_gpu_mem_kind_t kind, size_t size) {
    md_gpu_mem_t out = {0, NULL};
    if (!s) { md_gpu_heap_fail(h, "md_gpu_malloc: null stream"); return out; }
    if ((unsigned)kind >= (unsigned)MD_GPU_MEM_KIND_COUNT) { md_gpu_heap_fail(h, "md_gpu_malloc: invalid memory kind %d", (int)kind); return out; }
    if (size == 0) return out;
    const uint64_t need = md_gpu_heap_align_up(size, MD_GPU_HEAP_ALIGN);

    md_gpu_heap_range_t* r = (md_gpu_heap_range_t*)md_alloc(h->alloc, sizeof(md_gpu_heap_range_t));
    if (!r) { md_gpu_heap_fail(h, "out of memory"); return out; }

    md_mutex_lock(h->mutex);
    md_gpu_heap_poll_locked(h);
    md_gpu_heap_kind_t* hk = &h->kinds[kind];
    md_tlsf_node_t* node = NULL;
    md_gpu_chunk_t* c = NULL;
    if (h->debug) {
        /* Memory of its own, exactly `size`, never reused. */
        c = h->backend.chunk_create(h->backend.user, kind, size);
        md_gpu_chunk_t** slot = c ? (md_gpu_chunk_t**)md_gpu_heap_vec_push(&hk->chunks, h->alloc) : NULL;
        if (slot) { *slot = c; hk->reserved += size; }
        else if (c) { h->backend.chunk_destroy(h->backend.user, c); c = NULL; md_gpu_heap_fail(h, "out of memory"); }
    } else {
        node = md_gpu_heap_alloc_locked(h, s, kind, need);
        c = node ? (md_gpu_chunk_t*)node->region : NULL;
    }
    if (!c) {
        md_mutex_unlock(h->mutex);
        md_free(h->alloc, r, sizeof(*r));
        return out;
    }
    const uint64_t offset = node ? node->offset : 0;
    r->address      = c->address + offset;
    r->size         = size;
    r->chunk        = c;
    r->chunk_offset = offset;
    r->node         = node;
    r->own          = node ? NULL : c;
    r->kind         = kind;
    if (!md_gpu_heap_registry_insert_locked(h, r)) {
        if (node) md_gpu_heap_release_node_locked(h, kind, node);
        else      md_gpu_heap_own_release_locked(h, kind, c);
        md_mutex_unlock(h->mutex);
        md_free(h->alloc, r, sizeof(*r));
        md_gpu_heap_fail(h, "out of memory");
        return out;
    }
    hk->in_use += node ? node->size : size;
    if (hk->in_use > hk->peak_in_use) hk->peak_in_use = hk->in_use;
    hk->allocations++;
    md_mutex_unlock(h->mutex);

    out.gpu = r->address;
    out.cpu = c->host ? c->host + offset : NULL;
    return out;
}

void md_gpu_heap_free(md_gpu_heap_t* h, md_gpu_stream_t s, md_gpu_addr_t addr) {
    if (!addr) return;
    if (!s) { md_gpu_heap_fail(h, "md_gpu_free: a stream is required"); return; }

    md_mutex_lock(h->mutex);
    md_gpu_heap_range_t* r = md_gpu_heap_registry_find_locked(h, addr);
    if (!r || !(r->node || r->own) || r->address != addr) {
        md_mutex_unlock(h->mutex);
        md_gpu_heap_fail(h, "md_gpu_free: 0x%llx is not the start of a live md_gpu_malloc allocation", (unsigned long long)addr);
        return;
    }
    md_gpu_heap_registry_remove_locked(h, r);
    md_gpu_heap_kind_t* hk = &h->kinds[r->kind];
    const uint64_t bytes = r->node ? r->node->size : r->size;
    hk->in_use -= bytes;
    hk->allocations--;

    const md_gpu_heap_pending_free_t free_ = {r->node, r->own, r->kind, s, h->backend.stream_position(s)};
    bool queued = false;
    if (free_.value != 0 && h->backend.stream_completed(s) < free_.value) {
        md_gpu_heap_pending_free_t* pf = (md_gpu_heap_pending_free_t*)md_gpu_heap_vec_push(&h->pending_frees, h->alloc);
        if (pf) {
            *pf = free_;
        } else {
            /* Cannot remember it: leak the memory rather than free it under the GPU. */
            md_gpu_heap_fail(h, "out of memory recording a free; %llu bytes leaked", (unsigned long long)bytes);
        }
        queued = true;
    }
    if (!queued) md_gpu_heap_release_free_locked(h, &free_);
    md_mutex_unlock(h->mutex);
    md_free(h->alloc, r, sizeof(*r));
}

void md_gpu_heap_stats_locked(const md_gpu_heap_t* h, md_gpu_mem_kind_t kind, md_gpu_memory_stats_t* out) {
    const md_gpu_heap_kind_t* hk = &h->kinds[kind];
    out->bytes_in_use      = hk->in_use;
    out->bytes_peak_in_use = hk->peak_in_use;
    out->bytes_reserved    = hk->reserved;
    out->bytes_temp        = hk->temp_bytes;
    out->allocations       = hk->allocations;
    out->chunks            = (uint32_t)hk->chunks.count;
}

/* ---- Temp arenas -----------------------------------------------------------------

   Per stream and per kind, a stack of chunks being filled by the open scopes.
   A chunk belongs to the scope depth that acquired it. Inner scopes bump into
   an outer scope's chunk as well; when an inner scope ends, only chunks it
   acquired itself are retired, and what it used of an outer chunk stays used
   until the outer scope ends. That costs at most one chunk's tail per level
   and needs no per-allocation bookkeeping. A retired chunk is stamped with
   the stream's position and reused only once the GPU has passed it -- the CPU
   may write it at once, so stream order alone is not enough. */

static int md_gpu_heap_temp_index(md_gpu_mem_kind_t kind) {
    return kind == MD_GPU_MEM_DEVICE ? 0 : (kind == MD_GPU_MEM_HOST_WRITE ? 1 : -1);
}

static void md_gpu_heap_temp_chunk_destroy_locked(md_gpu_heap_t* h, md_gpu_chunk_t* c) {
    if (c->range) {
        md_gpu_heap_registry_remove_locked(h, c->range);
        md_free(h->alloc, c->range, sizeof(md_gpu_heap_range_t));
    }
    h->kinds[c->kind].temp_bytes -= c->size;
    h->backend.chunk_destroy(h->backend.user, c);
}

static md_gpu_chunk_t* md_gpu_heap_temp_acquire(md_gpu_heap_t* h, md_gpu_heap_temps_t* t, md_gpu_stream_t s,
                                                md_gpu_mem_kind_t kind, uint64_t need) {
    md_gpu_heap_vec_t* active = &t->arena[md_gpu_heap_temp_index(kind)].active;
    md_gpu_heap_vec_t* spare  = &t->arena[md_gpu_heap_temp_index(kind)].spare;
    const uint64_t done = h->backend.stream_completed(s);

    /* 1. A retired chunk the GPU is done with: start it over.
       2. A retired chunk with room left after its cursor: carry on filling it.
          The bytes before the cursor may still be in use, but the tail never
          was; the chunk is restamped when this scope ends, and only a later,
          complete retirement ever rewinds it. Small scopes (a frame's worth of
          constants) thus share a chunk instead of taking one each.
       3. A new chunk.
       Debug memory takes a new chunk, exactly `need`, every time; the
       retired ones the GPU is done with go back to the backend. */
    md_gpu_chunk_t* c = NULL;
    if (h->debug) {
        md_mutex_lock(h->mutex);
        for (size_t i = 0; i < spare->count;) {
            md_gpu_chunk_t* sp = MD_GPU_HEAP_VEC_AT(*spare, md_gpu_chunk_t*, i);
            if (sp->retire_value <= done) { md_gpu_heap_vec_remove(spare, i); md_gpu_heap_temp_chunk_destroy_locked(h, sp); }
            else ++i;
        }
        md_mutex_unlock(h->mutex);
    } else {
        for (size_t i = 0; i < spare->count && !c; ++i) {
            md_gpu_chunk_t* sp = MD_GPU_HEAP_VEC_AT(*spare, md_gpu_chunk_t*, i);
            if (sp->retire_value <= done && sp->size >= need) {
                c = sp;
                c->cursor = 0;
                md_gpu_heap_vec_remove(spare, i);
            }
        }
        for (size_t i = spare->count; i-- > 0 && !c;) {
            md_gpu_chunk_t* sp = MD_GPU_HEAP_VEC_AT(*spare, md_gpu_chunk_t*, i);
            if (sp->cursor + need <= sp->size) {
                c = sp;
                md_gpu_heap_vec_remove(spare, i);
            }
        }
    }
    if (!c) {
        const uint64_t size = h->debug ? need
                            : need > MD_GPU_HEAP_TEMP_CHUNK_MIN ? md_gpu_heap_align_up(need, MD_GPU_HEAP_LARGE_ALIGN)
                                                                : MD_GPU_HEAP_TEMP_CHUNK_MIN;
        md_gpu_heap_range_t* r = (md_gpu_heap_range_t*)md_alloc(h->alloc, sizeof(md_gpu_heap_range_t));
        if (!r) { md_gpu_heap_fail(h, "out of memory"); return NULL; }
        md_mutex_lock(h->mutex);
        c = h->backend.chunk_create(h->backend.user, kind, size);
        bool ok = c != NULL;
        if (ok) {
            memset(r, 0, sizeof(*r));
            r->address = c->address;
            r->size    = c->size;
            r->chunk   = c;
            r->kind    = kind;
            ok = md_gpu_heap_registry_insert_locked(h, r);
            if (ok) { c->range = r; h->kinds[kind].temp_bytes += c->size; }
            else    { h->backend.chunk_destroy(h->backend.user, c); md_gpu_heap_fail(h, "out of memory"); }
        }
        md_mutex_unlock(h->mutex);
        if (!ok) { md_free(h->alloc, r, sizeof(*r)); return NULL; }
    }
    md_gpu_chunk_t** slot = (md_gpu_chunk_t**)md_gpu_heap_vec_push(active, h->alloc);
    if (!slot) {
        md_mutex_lock(h->mutex);
        md_gpu_heap_temp_chunk_destroy_locked(h, c);
        md_mutex_unlock(h->mutex);
        md_gpu_heap_fail(h, "out of memory");
        return NULL;
    }
    *slot = c;
    c->retire_value = 0;
    c->owner_depth  = t->depth;
    return c;
}

void md_gpu_heap_temps_init(md_gpu_heap_temps_t* t) {
    for (uint32_t k = 0; k < MD_GPU_HEAP_TEMP_KINDS; ++k) {
        md_gpu_heap_vec_init(&t->arena[k].active, sizeof(md_gpu_chunk_t*));
        md_gpu_heap_vec_init(&t->arena[k].spare,  sizeof(md_gpu_chunk_t*));
    }
    t->depth = 0;
}

md_gpu_temp_t md_gpu_heap_temp_begin(md_gpu_heap_temps_t* t, md_gpu_stream_t s) {
    md_gpu_temp_t scope = {s, ++t->depth};
    return scope;
}

md_gpu_mem_t md_gpu_heap_temp_alloc(md_gpu_heap_t* h, md_gpu_heap_temps_t* t, md_gpu_temp_t scope,
                                    md_gpu_mem_kind_t kind, size_t size, const char* label) {
    md_gpu_mem_t out = {0, NULL};
    /* Chunks belong to the innermost scope, so memory for an outer one would be
       reclaimed by the inner one's end. */
    if (scope.depth == 0 || scope.depth != t->depth) {
        md_gpu_heap_fail(h, "md_gpu_temp_alloc: stream '%s': the scope (depth %u) is not the innermost open one (%u)",
                         label, scope.depth, t->depth);
        return out;
    }
    const int ti = md_gpu_heap_temp_index(kind);
    if (ti < 0) { md_gpu_heap_fail(h, "md_gpu_temp_alloc: only MD_GPU_MEM_DEVICE and MD_GPU_MEM_HOST_WRITE memory is temporary"); return out; }
    if (size == 0) return out;
    /* Debug memory: every allocation takes a chunk of exactly its size. */
    const uint64_t need = h->debug ? size : md_gpu_heap_align_up(size, MD_GPU_HEAP_ALIGN);

    md_gpu_heap_vec_t* active = &t->arena[ti].active;
    md_gpu_chunk_t* c = active->count ? MD_GPU_HEAP_VEC_AT(*active, md_gpu_chunk_t*, active->count - 1) : NULL;
    if (!c || c->cursor + need > c->size || h->debug) {
        c = md_gpu_heap_temp_acquire(h, t, scope.stream, kind, need);
        if (!c) return out;
    }
    out.gpu = c->address + c->cursor;
    out.cpu = c->host ? c->host + c->cursor : NULL;
    c->cursor += need;
    return out;
}

void md_gpu_heap_temp_end(md_gpu_heap_t* h, md_gpu_heap_temps_t* t, md_gpu_temp_t scope) {
    if (scope.depth == 0 || scope.depth != t->depth) {
        md_gpu_heap_fail(h, "md_gpu_temp_end: scopes must end in reverse order (ending depth %u, innermost open is %u)",
                         scope.depth, t->depth);
        return;
    }
    const uint64_t at = h->backend.stream_position(scope.stream);
    for (uint32_t k = 0; k < MD_GPU_HEAP_TEMP_KINDS; ++k) {
        md_gpu_heap_vec_t* active = &t->arena[k].active;
        while (active->count > 0) {
            md_gpu_chunk_t* c = MD_GPU_HEAP_VEC_AT(*active, md_gpu_chunk_t*, active->count - 1);
            if (c->owner_depth < scope.depth) break;
            active->count--;
            c->retire_value = at;
            md_gpu_chunk_t** slot = (md_gpu_chunk_t**)md_gpu_heap_vec_push(&t->arena[k].spare, h->alloc);
            if (slot) {
                *slot = c;
            } else {
                /* Cannot keep it for reuse, and cannot free it under the GPU: leak. */
                md_gpu_heap_fail(h, "out of memory retiring a temp chunk; %llu bytes leaked", (unsigned long long)c->size);
            }
        }
    }
    t->depth--;
}

void md_gpu_heap_temps_free(md_gpu_heap_t* h, md_gpu_heap_temps_t* t) {
    md_mutex_lock(h->mutex);
    for (uint32_t k = 0; k < MD_GPU_HEAP_TEMP_KINDS; ++k) {
        md_gpu_heap_vec_t* active = &t->arena[k].active;
        md_gpu_heap_vec_t* spare  = &t->arena[k].spare;
        for (size_t i = 0; i < active->count; ++i) md_gpu_heap_temp_chunk_destroy_locked(h, MD_GPU_HEAP_VEC_AT(*active, md_gpu_chunk_t*, i));
        for (size_t i = 0; i < spare->count;  ++i) md_gpu_heap_temp_chunk_destroy_locked(h, MD_GPU_HEAP_VEC_AT(*spare,  md_gpu_chunk_t*, i));
        md_gpu_heap_vec_free(active, h->alloc);
        md_gpu_heap_vec_free(spare,  h->alloc);
    }
    md_mutex_unlock(h->mutex);
    t->depth = 0;
}

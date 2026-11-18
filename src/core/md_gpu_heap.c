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
#include <core/md_common.h>
#include <core/md_format.h>
#include <core/md_os.h>

#include <stdarg.h>

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

static bool md_gpu_heap_fail(const md_gpu_heap_t* h, const char* fmt, ...) {
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    md_vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    h->backend.fail(buf);
    return false;
}

/* ---- Arrays --------------------------------------------------------------------
   md_array.h has no ordered removal: the registry stays sorted and the temp
   arenas are stacks. */

#define MD_GPU_HEAP_REMOVE_AT(arr, i) \
    (md_memmove((arr).data + (i), (arr).data + (i) + 1, ((arr).size - (i) - 1) * sizeof(*(arr).data)), (arr).size--)

static void md_gpu_heap_remove_chunk(md_gpu_chunk_array_t* arr, const md_gpu_chunk_t* c) {
    for (size_t i = 0; i < arr->size; ++i) {
        if (arr->data[i] == c) { MD_GPU_HEAP_REMOVE_AT(*arr, i); return; }
    }
}

/* ---- Registry ------------------------------------------------------------------ */

/* Binary search for the range containing `address`. */
static md_gpu_heap_range_t* md_gpu_heap_registry_find_locked(const md_gpu_heap_t* h, uint64_t address) {
    size_t lo = 0, hi = h->registry.size;
    md_gpu_heap_range_t* const* arr = h->registry.data;
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
    size_t lo = 0, hi = h->registry.size;
    md_gpu_heap_range_t* const* arr = h->registry.data;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (arr[mid]->address < address) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static bool md_gpu_heap_registry_insert_locked(md_gpu_heap_t* h, md_gpu_heap_range_t* r) {
    const size_t at = md_gpu_heap_registry_lower_bound(h, r->address);
    if (!MD_ARRAY_TRY_EXTEND(h->registry, 1)) return false;
    md_memmove(h->registry.data + at + 1, h->registry.data + at, (h->registry.size - 1 - at) * sizeof(*h->registry.data));
    h->registry.data[at] = r;
    return true;
}

static void md_gpu_heap_registry_remove_locked(md_gpu_heap_t* h, md_gpu_heap_range_t* r) {
    const size_t at = md_gpu_heap_registry_lower_bound(h, r->address);
    if (at < h->registry.size && h->registry.data[at] == r) MD_GPU_HEAP_REMOVE_AT(h->registry, at);
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
    md_memset(h, 0, sizeof(*h));
    h->alloc       = alloc;
    h->mutex       = mutex;
    h->backend     = *backend;
    h->cache_limit = cache_limit ? cache_limit : MD_GPU_HEAP_CACHE_DEFAULT;
    h->debug       = debug;
    MD_ARRAY_INIT(h->registry, alloc);
    MD_ARRAY_INIT(h->pending_frees, alloc);
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_tlsf_init(&h->kinds[k].tlsf, alloc, MD_GPU_HEAP_ALIGN);
        MD_ARRAY_INIT(h->kinds[k].chunks, alloc);
    }
}

void md_gpu_heap_destroy(md_gpu_heap_t* h) {
    md_mutex_lock(h->mutex);
    for (size_t i = 0; i < h->registry.size; ++i) md_free(h->alloc, h->registry.data[i], sizeof(md_gpu_heap_range_t));
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_gpu_heap_kind_t* hk = &h->kinds[k];
        for (size_t i = 0; i < hk->chunks.size; ++i) h->backend.chunk_destroy(h->backend.user, hk->chunks.data[i]);
        (void)MD_ARRAY_TRY_FREE(hk->chunks);
        md_tlsf_destroy(&hk->tlsf);
    }
    (void)MD_ARRAY_TRY_FREE(h->registry);
    (void)MD_ARRAY_TRY_FREE(h->pending_frees);
    md_mutex_unlock(h->mutex);
}

static uint64_t md_gpu_heap_chunk_size(const md_gpu_heap_kind_t* hk, uint64_t need) {
    uint64_t size = MD_GPU_HEAP_CHUNK_MIN;
    for (size_t i = 0; i < hk->chunks.size && size < MD_GPU_HEAP_CHUNK_MAX; ++i) size *= 2;
    if (need > size) size = MD_ALIGN_TO(need, MD_GPU_HEAP_LARGE_ALIGN);
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
        md_gpu_heap_remove_chunk(&hk->chunks, c);
        hk->empty_bytes -= c->size;
        hk->reserved    -= c->size;
        h->backend.chunk_destroy(h->backend.user, c);
    }
}

/* Give a debug memory allocation's chunk back: the GPU is done with it. */
static void md_gpu_heap_own_release_locked(md_gpu_heap_t* h, md_gpu_mem_kind_t kind, md_gpu_chunk_t* c) {
    md_gpu_heap_kind_t* hk = &h->kinds[kind];
    md_gpu_heap_remove_chunk(&hk->chunks, c);
    hk->reserved -= c->size;
    h->backend.chunk_destroy(h->backend.user, c);
}

static void md_gpu_heap_release_free_locked(md_gpu_heap_t* h, const md_gpu_heap_pending_free_t* pf) {
    if (pf->own) md_gpu_heap_own_release_locked(h, pf->kind, pf->own);
    else         md_gpu_heap_release_node_locked(h, pf->kind, pf->node);
}

void md_gpu_heap_poll_locked(md_gpu_heap_t* h) {
    for (size_t i = 0; i < h->pending_frees.size;) {
        const md_gpu_heap_pending_free_t pf = h->pending_frees.data[i];
        if (h->backend.stream_completed(pf.stream) >= pf.value) {
            MD_GPU_HEAP_REMOVE_AT(h->pending_frees, i);
            md_gpu_heap_release_free_locked(h, &pf);
        } else {
            ++i;
        }
    }
}

void md_gpu_heap_forget_stream_locked(md_gpu_heap_t* h, md_gpu_stream_t s) {
    for (size_t i = 0; i < h->pending_frees.size;) {
        const md_gpu_heap_pending_free_t pf = h->pending_frees.data[i];
        if (pf.stream == s) {
            MD_GPU_HEAP_REMOVE_AT(h->pending_frees, i);
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
        for (size_t i = 0; i < h->pending_frees.size; ++i) {
            const md_gpu_heap_pending_free_t* pf = &h->pending_frees.data[i];
            if (pf->stream != s || pf->kind != kind || pf->own) continue;
            if (pf->node->size < need || pf->node->size > need + need / 4 + MD_GPU_HEAP_ALIGN) continue;
            if (!h->backend.reuse_early(s)) break;
            node = pf->node;
            MD_GPU_HEAP_REMOVE_AT(h->pending_frees, i);
            break;
        }
    }
    if (!node) node = md_tlsf_alloc(&hk->tlsf, need);
    if (!node) {
        const uint64_t size = md_gpu_heap_chunk_size(hk, need);
        md_gpu_chunk_t* c = h->backend.chunk_create(h->backend.user, kind, size);
        if (!c) return NULL;
        const bool listed = MD_ARRAY_TRY_PUSH(hk->chunks, c);
        if (!listed || !md_tlsf_add_region(&hk->tlsf, c, size)) {
            if (listed) MD_ARRAY_POP(hk->chunks);
            h->backend.chunk_destroy(h->backend.user, c);
            md_gpu_heap_fail(h, "out of memory");
            return NULL;
        }
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
    const uint64_t need = MD_ALIGN_TO(size, MD_GPU_HEAP_ALIGN);

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
        if (c && MD_ARRAY_TRY_PUSH(hk->chunks, c)) hk->reserved += size;
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
        if (!MD_ARRAY_TRY_PUSH(h->pending_frees, free_)) {
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
    out->chunks            = (uint32_t)hk->chunks.size;
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
    md_gpu_chunk_array_t* active = &t->arena[md_gpu_heap_temp_index(kind)].active;
    md_gpu_chunk_array_t* spare  = &t->arena[md_gpu_heap_temp_index(kind)].spare;
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
        for (size_t i = 0; i < spare->size;) {
            md_gpu_chunk_t* sp = spare->data[i];
            if (sp->retire_value <= done) { MD_GPU_HEAP_REMOVE_AT(*spare, i); md_gpu_heap_temp_chunk_destroy_locked(h, sp); }
            else ++i;
        }
        md_mutex_unlock(h->mutex);
    } else {
        for (size_t i = 0; i < spare->size && !c; ++i) {
            md_gpu_chunk_t* sp = spare->data[i];
            if (sp->retire_value <= done && sp->size >= need) {
                c = sp;
                c->cursor = 0;
                MD_GPU_HEAP_REMOVE_AT(*spare, i);
            }
        }
        for (size_t i = spare->size; i-- > 0 && !c;) {
            md_gpu_chunk_t* sp = spare->data[i];
            if (sp->cursor + need <= sp->size) {
                c = sp;
                MD_GPU_HEAP_REMOVE_AT(*spare, i);
            }
        }
    }
    if (!c) {
        const uint64_t size = h->debug ? need
                            : need > MD_GPU_HEAP_TEMP_CHUNK_MIN ? MD_ALIGN_TO(need, MD_GPU_HEAP_LARGE_ALIGN)
                                                                : MD_GPU_HEAP_TEMP_CHUNK_MIN;
        md_gpu_heap_range_t* r = (md_gpu_heap_range_t*)md_alloc(h->alloc, sizeof(md_gpu_heap_range_t));
        if (!r) { md_gpu_heap_fail(h, "out of memory"); return NULL; }
        md_mutex_lock(h->mutex);
        c = h->backend.chunk_create(h->backend.user, kind, size);
        bool ok = c != NULL;
        if (ok) {
            md_memset(r, 0, sizeof(*r));
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
    if (!MD_ARRAY_TRY_PUSH(*active, c)) {
        md_mutex_lock(h->mutex);
        md_gpu_heap_temp_chunk_destroy_locked(h, c);
        md_mutex_unlock(h->mutex);
        md_gpu_heap_fail(h, "out of memory");
        return NULL;
    }
    c->retire_value = 0;
    c->owner_depth  = t->depth;
    return c;
}

void md_gpu_heap_temps_init(const md_gpu_heap_t* h, md_gpu_heap_temps_t* t) {
    for (uint32_t k = 0; k < MD_GPU_HEAP_TEMP_KINDS; ++k) {
        MD_ARRAY_INIT(t->arena[k].active, h->alloc);
        MD_ARRAY_INIT(t->arena[k].spare,  h->alloc);
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
    const uint64_t need = h->debug ? size : MD_ALIGN_TO(size, MD_GPU_HEAP_ALIGN);

    const md_gpu_chunk_array_t* active = &t->arena[ti].active;
    md_gpu_chunk_t* c = active->size ? MD_ARRAY_BACK(*active) : NULL;
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
        md_gpu_chunk_array_t* active = &t->arena[k].active;
        while (active->size > 0) {
            md_gpu_chunk_t* c = MD_ARRAY_BACK(*active);
            if (c->owner_depth < scope.depth) break;
            MD_ARRAY_POP(*active);
            c->retire_value = at;
            if (!MD_ARRAY_TRY_PUSH(t->arena[k].spare, c)) {
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
        md_gpu_chunk_array_t* active = &t->arena[k].active;
        md_gpu_chunk_array_t* spare  = &t->arena[k].spare;
        for (size_t i = 0; i < active->size; ++i) md_gpu_heap_temp_chunk_destroy_locked(h, active->data[i]);
        for (size_t i = 0; i < spare->size;  ++i) md_gpu_heap_temp_chunk_destroy_locked(h, spare->data[i]);
        (void)MD_ARRAY_TRY_FREE(*active);
        (void)MD_ARRAY_TRY_FREE(*spare);
    }
    md_mutex_unlock(h->mutex);
    t->depth = 0;
}

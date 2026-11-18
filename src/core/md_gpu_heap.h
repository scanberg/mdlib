/*
md_gpu_heap.h -- internal to the md_gpu backends.

The memory behind md_gpu_malloc, md_gpu_free and the temp scopes, the same on
every backend:

  - one heap per memory kind: a TLSF sub-allocator (md_gpu_tlsf.h) over large
    chunks, so that a thousand small allocations cost a handful of driver
    allocations;
  - stream-ordered frees, released once the GPU has passed them, and reused
    early by the freeing stream itself where the backend can order that;
  - per-stream temp arenas, bump-allocated within nested scopes;
  - the registry that maps a device address back to its allocation;
  - debug memory, where every allocation is a chunk of its own.

What a chunk *is* -- a VkBuffer with its VkDeviceMemory, an MTLBuffer over a
placement heap -- is the backend's. It makes and destroys chunks through
md_gpu_heap_backend_t, and reaches its own data by putting md_gpu_chunk_t
first in its chunk struct and casting.

Locking: everything here runs under the backend's device lock, handed over at
md_gpu_heap_init. Functions named _locked expect it held; the others take it.
The backend's callbacks are called with it held.
*/

#ifndef MD_GPU_HEAP_H
#define MD_GPU_HEAP_H

#include "md_gpu.h"
#include "md_gpu_tlsf.h"

#include <core/md_array.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct md_mutex_t;

#define MD_GPU_HEAP_ALIGN        256u              /* every md_gpu_malloc / temp allocation */
#define MD_GPU_HEAP_CHUNK_MIN    (4ull << 20)       /* first heap chunk of a kind            */
#define MD_GPU_HEAP_CHUNK_MAX    (64ull << 20)      /* chunks grow by doubling up to this    */
#define MD_GPU_HEAP_LARGE_ALIGN  (64ull << 10)      /* granularity of oversized chunks       */
#define MD_GPU_HEAP_CACHE_DEFAULT (256ull << 20)    /* empty chunks kept per kind            */
#define MD_GPU_HEAP_TEMP_CHUNK_MIN (4ull << 20)
#define MD_GPU_HEAP_TEMP_KINDS   2u                 /* DEVICE, HOST_WRITE                    */

/* The part of a chunk md_gpu_heap keeps. First in the backend's chunk struct. */
typedef struct md_gpu_chunk_t {
    /* Filled in by the backend's chunk_create. */
    uint64_t          address;
    uint8_t*          host;          /* mapped pointer, or NULL                 */
    uint64_t          size;
    md_gpu_mem_kind_t kind;
    /* Heap chunks. */
    bool              empty;         /* no live allocation in it                */
    md_tlsf_node_t*   empty_node;    /* its single free node while empty        */
    /* Temp chunks. */
    uint64_t          cursor;
    uint32_t          owner_depth;   /* scope depth that acquired it            */
    uint64_t          retire_value;  /* reusable once the stream reaches this   */
    struct md_gpu_heap_range_t* range;  /* registry entry                       */
} md_gpu_chunk_t;

/* Arrays (md_array.h, grown with the TRY_ operations only: running out of
   host memory is reported, not fatal). */
MD_ARRAY_RAW_DEFINE_TYPE(md_gpu_chunk_array_t, md_gpu_chunk_t*);
MD_ARRAY_RAW_DEFINE_TYPE(md_gpu_heap_range_array_t, struct md_gpu_heap_range_t*);
MD_ARRAY_RAW_DEFINE_TYPE(md_gpu_heap_pending_free_array_t, struct md_gpu_heap_pending_free_t);

typedef struct md_gpu_heap_backend_t {
    void* user;   /* passed back to chunk_create and chunk_destroy */
    /* A chunk of `size` bytes of `kind`: a zeroed struct of the backend's with
       md_gpu_chunk_t first and its address, host, size and kind filled in.
       NULL on failure, reported (md_gpu_fail). */
    md_gpu_chunk_t* (*chunk_create)(void* user, md_gpu_mem_kind_t kind, uint64_t size);
    void            (*chunk_destroy)(void* user, md_gpu_chunk_t* chunk);
    /* The point a free or a temp scope's end issued now waits for, and how far
       the GPU has got, on the same timeline. */
    uint64_t        (*stream_position)(md_gpu_stream_t s);
    uint64_t        (*stream_completed)(md_gpu_stream_t s);
    /* `s` is about to reuse DEVICE memory it freed before the GPU passed the
       free. Order what follows after what came before (a barrier in EXPLICIT
       mode) and return true, or return false when that cannot be done now. */
    bool            (*reuse_early)(md_gpu_stream_t s);
} md_gpu_heap_backend_t;

/* The device-wide heap for one memory kind. */
typedef struct md_gpu_heap_kind_t {
    md_tlsf_t            tlsf;
    md_gpu_chunk_array_t chunks;
    uint64_t          reserved;
    uint64_t          empty_bytes;   /* bytes of chunks with nothing in them    */
    uint64_t          in_use;
    uint64_t          peak_in_use;
    uint32_t          allocations;
    uint64_t          temp_bytes;    /* temp arenas of every stream, this kind  */
} md_gpu_heap_kind_t;

typedef struct md_gpu_heap_t {
    struct md_allocator_i* alloc;
    struct md_mutex_t*     mutex;
    md_gpu_heap_backend_t  backend;
    md_gpu_heap_kind_t     kinds[MD_GPU_MEM_KIND_COUNT];
    md_gpu_heap_pending_free_array_t pending_frees;
    md_gpu_heap_range_array_t        registry;   /* sorted by address */
    uint64_t               cache_limit;
    bool                   debug;           /* debug memory: a chunk per allocation */
} md_gpu_heap_t;

/* A stream's temp arenas. Per kind, `active` is a stack of the chunks the
   open scopes are filling, in the order they were acquired; `spare` holds
   retired ones. */
typedef struct md_gpu_heap_temps_t {
    struct {
        md_gpu_chunk_array_t active;
        md_gpu_chunk_array_t spare;
    } arena[MD_GPU_HEAP_TEMP_KINDS];
    uint32_t depth;                  /* open scopes */
} md_gpu_heap_temps_t;

/* Where a resolved address lives. */
typedef struct md_gpu_heap_span_t {
    md_gpu_chunk_t* chunk;
    uint64_t        offset;          /* of the address within the chunk */
    uint8_t*        host;            /* host pointer of the address, or NULL */
} md_gpu_heap_span_t;

/* `cache_limit` 0 selects MD_GPU_HEAP_CACHE_DEFAULT. */
void md_gpu_heap_init(md_gpu_heap_t* h, struct md_allocator_i* alloc, struct md_mutex_t* mutex,
                      const md_gpu_heap_backend_t* backend, uint64_t cache_limit, bool debug);
/* Every chunk back to the backend. Every stream is idle and its temps freed. */
void md_gpu_heap_destroy(md_gpu_heap_t* h);

/* md_gpu_malloc and md_gpu_free. */
md_gpu_mem_t md_gpu_heap_malloc(md_gpu_heap_t* h, md_gpu_stream_t s, md_gpu_mem_kind_t kind, size_t size);
void         md_gpu_heap_free(md_gpu_heap_t* h, md_gpu_stream_t s, md_gpu_addr_t addr);

/* [addr, addr + size) within one live allocation, or a failure naming `what`. */
bool md_gpu_heap_resolve(md_gpu_heap_t* h, md_gpu_addr_t addr, uint64_t size, md_gpu_heap_span_t* out, const char* what);

/* Release the frees the GPU has passed. */
void md_gpu_heap_poll_locked(md_gpu_heap_t* h);
/* `s` is going away, synchronised: release what it freed. */
void md_gpu_heap_forget_stream_locked(md_gpu_heap_t* h, md_gpu_stream_t s);
/* All but bytes_textures, which is the backend's. */
void md_gpu_heap_stats_locked(const md_gpu_heap_t* h, md_gpu_mem_kind_t kind, md_gpu_memory_stats_t* out);

/* The temp scopes of a stream. `label` names the stream in errors. */
void          md_gpu_heap_temps_init(const md_gpu_heap_t* h, md_gpu_heap_temps_t* t);
md_gpu_temp_t md_gpu_heap_temp_begin(md_gpu_heap_temps_t* t, md_gpu_stream_t s);
md_gpu_mem_t  md_gpu_heap_temp_alloc(md_gpu_heap_t* h, md_gpu_heap_temps_t* t, md_gpu_temp_t scope,
                                     md_gpu_mem_kind_t kind, size_t size, const char* label);
void          md_gpu_heap_temp_end(md_gpu_heap_t* h, md_gpu_heap_temps_t* t, md_gpu_temp_t scope);
/* Every chunk of the stream's temp arenas. The stream is idle. */
void          md_gpu_heap_temps_free(md_gpu_heap_t* h, md_gpu_heap_temps_t* t);

#endif

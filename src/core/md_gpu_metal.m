/*
md_gpu_metal.m — Metal backend for md_gpu.h

Mirrors md_gpu_vulkan.c section for section. Where the two backends differ,
the difference is called out in a comment.

Three structural differences from Vulkan:

  1. Program order is a memoryBarrierWithScope: inside a compute encoder, and
     an encoder boundary when switching between compute and blit work. Metal
     orders encoders within a command buffer, so the boundary is free.

  2. There are no image layouts, so the whole layout question disappears
     rather than being answered with VK_IMAGE_LAYOUT_GENERAL.

  3. MTLCommandBuffer is single-use, so a graph cannot be a pre-recorded
     command buffer the way it is on Vulkan. A graph here is a recorded list
     of operations that is re-encoded on each launch. Encoding is cheap on
     Metal by design, and the CPU work a graph actually eliminates -- building
     argument blocks, hashing, allocating -- is eliminated either way, because
     argument blocks live in graph-owned buffers and are never rebuilt.

NOTE ON RESIDENCY: with no per-dispatch resource declarations there is nothing
to derive residency from, so every live allocation goes into a device-wide
MTLResidencySet attached to both queues (macOS 15 / iOS 18 and later). On
older systems the backend falls back to calling useResources: once per encoder
with the full live set, which is correct but O(live allocations) per encoder.
*/

#include "md_gpu.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
   1. Configuration, error handling, utilities
   ========================================================================= */

/* Metal needs no descriptor heap at all: Slang lowers DescriptorHandle<T> to
   the texture or sampler object sitting directly in the argument struct (a
   Metal 3 argument buffer entry), so the handle a shader receives is simply
   the resource's MTLResourceID. There is no table, no index and no binding.

   Requires argument buffers tier 2, i.e. Metal 3 -- true of all Apple silicon. */
#define MD_MTL_ARENA_PAGE_SIZE (256u * 1024u)
#define MD_MTL_GRAPH_ARG_PAGE   (64u * 1024u)
#define MD_MTL_ARG_ALIGN          64u
#define MD_MTL_MAX_PENDING_WAITS  16u
#define MD_MTL_ERROR_BUF         512u

/* Slang emits the push-constant block as buffer argument 0 on Metal. */
#define MD_MTL_ARG_BUFFER_INDEX    0

static __thread char md_mtl_error_buf[MD_MTL_ERROR_BUF];
static __thread bool md_mtl_has_error;

static bool md_mtl_fail(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(md_mtl_error_buf, sizeof(md_mtl_error_buf), fmt, ap);
    va_end(ap);
    md_mtl_has_error = true;
    MD_LOG_ERROR("md_gpu: %s", md_mtl_error_buf);
    return false;
}

const char* md_gpu_last_error(void) {
    return md_mtl_has_error ? md_mtl_error_buf : NULL;
}

static inline uint64_t md_mtl_align_up(uint64_t v, uint64_t a) { return (v + a - 1) & ~(a - 1); }

static inline uint64_t md_mtl_next_pow2(uint64_t v) {
    if (v < 256) return 256;
    v--;
    v |= v >> 1;  v |= v >> 2;  v |= v >> 4;
    v |= v >> 8;  v |= v >> 16; v |= v >> 32;
    return v + 1;
}

typedef struct md_mtl_vec_t {
    void*  data;
    size_t count;
    size_t capacity;
    size_t stride;
} md_mtl_vec_t;

static void md_mtl_vec_init(md_mtl_vec_t* v, size_t stride) {
    v->data = NULL; v->count = 0; v->capacity = 0; v->stride = stride;
}

static bool md_mtl_vec_reserve(md_mtl_vec_t* v, struct md_allocator_i* alloc, size_t n) {
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

static void* md_mtl_vec_push(md_mtl_vec_t* v, struct md_allocator_i* alloc) {
    if (!md_mtl_vec_reserve(v, alloc, v->count + 1)) return NULL;
    void* slot = (char*)v->data + v->count * v->stride;
    memset(slot, 0, v->stride);
    v->count++;
    return slot;
}

static void md_mtl_vec_free(md_mtl_vec_t* v, struct md_allocator_i* alloc) {
    if (v->data) md_free(alloc, v->data, v->capacity * v->stride);
    v->data = NULL; v->count = 0; v->capacity = 0;
}

#define MD_MTL_VEC_AT(v, type, i) (((type*)(v).data)[i])

/* =========================================================================
   2. Types
   ========================================================================= */

typedef struct md_mtl_block_t {
    __unsafe_unretained id<MTLBuffer> buffer;
    uint64_t            address;
    void*               host;
    uint64_t            capacity;
    uint64_t            size;
    md_gpu_mem_flags_t flags;
    md_gpu_pool_t      pool;
    bool                in_use;
    md_gpu_stream_t    free_stream;
    uint64_t            free_value;
} md_mtl_block_t;

typedef struct md_gpu_pool {
    md_gpu_device_t device;
    md_gpu_mem_flags_t flags;
    uint64_t         block_size;
    uint64_t         release_threshold;
    md_mtl_vec_t     blocks;          /* md_mtl_block_t* */
    uint64_t         in_use_bytes;
    uint64_t         reserved_bytes;
    uint64_t         peak_in_use_bytes;
    uint64_t         alloc_count;
    uint64_t         reuse_count;
    char             label[64];
} md_gpu_pool;

typedef struct md_mtl_page_t {
    __unsafe_unretained id<MTLBuffer> buffer;
    uint64_t address;
    uint8_t* host;
    uint64_t capacity;
    uint64_t cursor;
    uint64_t retire_value;
} md_mtl_page_t;

typedef struct md_mtl_arena_t {
    md_mtl_vec_t pages;    /* md_mtl_page_t* */
    size_t       current;
} md_mtl_arena_t;

/* Recorded operations, used to replay a graph. */
typedef enum md_mtl_op_kind_t {
    MD_MTL_OP_DISPATCH,
    MD_MTL_OP_DISPATCH_INDIRECT,
    MD_MTL_OP_COPY_BUFFER,
    MD_MTL_OP_FILL_BUFFER,
    MD_MTL_OP_COPY_BUFFER_TO_TEX,
    MD_MTL_OP_COPY_TEX_TO_BUFFER,
} md_mtl_op_kind_t;

typedef struct md_mtl_op_t {
    md_mtl_op_kind_t kind;

    /* dispatch */
    md_gpu_kernel_t kernel;
    uint64_t         arg_addr;
    __unsafe_unretained id<MTLBuffer> arg_buffer;
    uint64_t         arg_offset;
    md_gpu_grid_t   grid;
    __unsafe_unretained id<MTLBuffer> indirect_buffer;
    uint64_t         indirect_offset;

    /* copies and fills */
    __unsafe_unretained id<MTLBuffer>  src_buffer;
    __unsafe_unretained id<MTLBuffer>  dst_buffer;
    __unsafe_unretained id<MTLTexture> texture;
    uint64_t src_offset, dst_offset, size;
    uint8_t  fill_value;
    MTLOrigin tex_origin;
    MTLSize   tex_size;
    uint64_t  bytes_per_row, bytes_per_image;
} md_mtl_op_t;

typedef struct md_mtl_graph_launch_t {
    size_t   page_index;
    uint64_t offset;
    size_t   size;
} md_mtl_graph_launch_t;

typedef struct md_gpu_graph {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    md_mtl_vec_t          arg_pages;   /* md_mtl_page_t* — owned, persistent */
    md_mtl_vec_t          launches;    /* md_mtl_graph_launch_t */
    md_mtl_vec_t          ops;         /* md_mtl_op_t */
    char                  label[64];
} md_gpu_graph;

typedef struct md_gpu_stream {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    __unsafe_unretained id<MTLCommandQueue>        queue;
    __unsafe_unretained id<MTLSharedEvent>         timeline;
    __unsafe_unretained id<MTLCommandBuffer>       cmd;
    __unsafe_unretained id<MTLComputeCommandEncoder> compute_enc;
    __unsafe_unretained id<MTLBlitCommandEncoder>    blit_enc;

    uint64_t next_value;
    uint64_t submitted_value;
    bool     has_work;
    bool     needs_barrier;
    bool     auto_order;

    __unsafe_unretained id<MTLSharedEvent> wait_events[MD_MTL_MAX_PENDING_WAITS];
    uint64_t wait_values[MD_MTL_MAX_PENDING_WAITS];
    uint32_t wait_count;

    md_mtl_arena_t arena;

    bool             upload_open;
    bool             upload_direct;
    md_gpu_ptr_t    upload_dst;
    uint64_t         upload_src_addr;
    size_t           upload_size;

    md_gpu_graph_t  capture;
    bool             owns_queue;
    bool             is_default;
    char             label[64];
} md_gpu_stream;

typedef struct md_mtl_tex_t {
    __unsafe_unretained id<MTLTexture> texture;
    md_gpu_tex_desc_t desc;
    uint64_t handle;   /* gpuResourceID -- the value shaders receive */
} md_mtl_tex_t;

typedef struct md_mtl_sampler_t {
    __unsafe_unretained id<MTLSamplerState> sampler;
    uint64_t handle;
} md_mtl_sampler_t;

typedef struct md_gpu_kernel {
    md_gpu_device_t device;
    __unsafe_unretained id<MTLComputePipelineState> pso;
    __unsafe_unretained id<MTLFunction>             fn;
    uint32_t group_size[3];
    char     label[64];
} md_gpu_kernel;

typedef struct md_mtl_hostfn_t {
    md_gpu_sync_t sync;
    md_gpu_host_fn fn;
    void* user;
    bool  internal;
} md_mtl_hostfn_t;

typedef struct md_mtl_retire_t {
    __unsafe_unretained id<MTLTexture> texture;
    md_gpu_sync_t waits[8];
    uint32_t       wait_count;
} md_mtl_retire_t;

typedef struct md_gpu_device {
    struct md_allocator_i* alloc;
    __unsafe_unretained id<MTLDevice>       device;
    /* Kept only so the first default stream of each kind has a queue to adopt;
       every stream creates its own. MTLCommandQueue is cheap. */
    __unsafe_unretained id<MTLCommandQueue> compute_queue;
    __unsafe_unretained id<MTLCommandQueue> transfer_queue;
    __unsafe_unretained id                  residency_set;   /* id<MTLResidencySet> */
    bool     has_residency_set;

    md_mutex_t queue_mutex;
    md_mutex_t device_mutex;

    /* Live texture / sampler records, looked up by resource id. Both stay
       small in practice, so a linear scan is cheaper than a hash. */
    md_mtl_vec_t textures;   /* md_mtl_tex_t*     */
    md_mtl_vec_t samplers;   /* md_mtl_sampler_t* */

    md_mtl_vec_t registry;   /* md_mtl_block_t*, sorted by address */
    md_mtl_vec_t pools;      /* md_gpu_pool_t  */
    md_mtl_vec_t streams;    /* md_gpu_stream_t */
    md_mtl_vec_t hostfns;    /* md_mtl_hostfn_t */
    md_mtl_vec_t retires;    /* md_mtl_retire_t */
    md_mtl_vec_t live_bufs;  /* id<MTLBuffer>, for the no-residency-set path */

    md_gpu_kernel_t make_grid_kernel;
    bool             is_discrete;
} md_gpu_device;

/* On Metal the handle is the resource id itself. Zero is the null handle. */
static inline bool md_mtl_tex_valid(uint64_t h) { return h != 0; }

static md_gpu_device_t md_mtl_devices[8];
static uint32_t         md_mtl_device_count;

/* Objective-C objects are held in CFRetain'd slots so that plain C structs can
   own them without ARC fighting us. */
#define MD_MTL_RETAIN(obj)  ((void)CFBridgingRetain(obj))
#define MD_MTL_RELEASE(obj) do { if (obj) CFRelease((__bridge CFTypeRef)(obj)); } while (0)

static uint32_t md_mtl_poll_internal(md_gpu_device_t dev, bool fire_user);
static bool     md_mtl_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host, id<MTLBuffer>* out_buf, uint64_t* out_off);
static bool     md_mtl_stream_submit(md_gpu_stream_t s);
static uint64_t md_mtl_stream_completed(md_gpu_stream_t s);

/* =========================================================================
   3. Allocation registry
   ========================================================================= */

static md_mtl_block_t* md_mtl_registry_find(md_gpu_device_t dev, uint64_t address) {
    size_t lo = 0, hi = dev->registry.count;
    md_mtl_block_t** arr = (md_mtl_block_t**)dev->registry.data;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        md_mtl_block_t* b = arr[mid];
        if (address < b->address)                        hi = mid;
        else if (address >= b->address + b->capacity)    lo = mid + 1;
        else                                             return b;
    }
    return NULL;
}

static bool md_mtl_registry_insert(md_gpu_device_t dev, md_mtl_block_t* blk) {
    if (!md_mtl_vec_reserve(&dev->registry, dev->alloc, dev->registry.count + 1)) return false;
    md_mtl_block_t** arr = (md_mtl_block_t**)dev->registry.data;
    size_t i = dev->registry.count;
    while (i > 0 && arr[i - 1]->address > blk->address) { arr[i] = arr[i - 1]; i--; }
    arr[i] = blk;
    dev->registry.count++;
    return true;
}

static void md_mtl_registry_remove(md_gpu_device_t dev, md_mtl_block_t* blk) {
    md_mtl_block_t** arr = (md_mtl_block_t**)dev->registry.data;
    for (size_t i = 0; i < dev->registry.count; ++i) {
        if (arr[i] == blk) {
            memmove(&arr[i], &arr[i + 1], (dev->registry.count - i - 1) * sizeof(*arr));
            dev->registry.count--;
            return;
        }
    }
}

static md_mtl_block_t* md_mtl_find_anywhere(uint64_t addr) {
    for (uint32_t i = 0; i < md_mtl_device_count; ++i) {
        md_gpu_device_t d = md_mtl_devices[i];
        if (!d) continue;
        md_mtl_block_t* b = md_mtl_registry_find(d, addr);
        if (b) return b;
    }
    return NULL;
}

/* =========================================================================
   4. Buffers, residency and transient arenas
   ========================================================================= */

static void md_mtl_make_resident(md_gpu_device_t dev, id<MTLAllocation> res) {
    if (dev->has_residency_set) {
        id set = dev->residency_set;
        [set performSelector:@selector(addAllocation:) withObject:res];
        [set performSelector:@selector(commit)];
    }
}

static void md_mtl_end_residency(md_gpu_device_t dev, id<MTLAllocation> res) {
    if (dev->has_residency_set) {
        id set = dev->residency_set;
        [set performSelector:@selector(removeAllocation:) withObject:res];
        [set performSelector:@selector(commit)];
    }
}

static bool md_mtl_create_raw_buffer(md_gpu_device_t dev, uint64_t size, md_gpu_mem_flags_t flags,
                                     id<MTLBuffer>* out_buf, uint64_t* out_addr, void** out_host) {
    MTLResourceOptions opts;
    if (flags & (MD_GPU_MEM_HOST_WRITE | MD_GPU_MEM_HOST_READ)) {
        opts = MTLResourceStorageModeShared;
    } else {
        opts = MTLResourceStorageModePrivate;
    }
    id<MTLBuffer> buf = [dev->device newBufferWithLength:(NSUInteger)size options:opts];
    if (!buf) return md_mtl_fail("newBufferWithLength failed for %llu bytes", (unsigned long long)size);

    MD_MTL_RETAIN(buf);
    md_mtl_make_resident(dev, buf);

    *out_buf  = buf;
    *out_addr = (uint64_t)[buf gpuAddress];
    *out_host = (opts == MTLResourceStorageModeShared) ? [buf contents] : NULL;
    return true;
}

static void md_mtl_destroy_raw_buffer(md_gpu_device_t dev, id<MTLBuffer> buf) {
    if (!buf) return;
    md_mtl_end_residency(dev, buf);
    MD_MTL_RELEASE(buf);
}

static md_mtl_page_t* md_mtl_page_create(md_gpu_device_t dev, uint64_t size) {
    md_mtl_page_t* p = (md_mtl_page_t*)md_alloc(dev->alloc, sizeof(md_mtl_page_t));
    if (!p) return NULL;
    memset(p, 0, sizeof(*p));
    void* host = NULL;
    id<MTLBuffer> buf = nil;
    if (!md_mtl_create_raw_buffer(dev, size, MD_GPU_MEM_HOST_WRITE, &buf, &p->address, &host)) {
        md_free(dev->alloc, p, sizeof(md_mtl_page_t));
        return NULL;
    }
    p->buffer   = buf;
    p->host     = (uint8_t*)host;
    p->capacity = size;
    return p;
}

static void md_mtl_page_destroy(md_gpu_device_t dev, md_mtl_page_t* p) {
    md_mtl_destroy_raw_buffer(dev, p->buffer);
    md_free(dev->alloc, p, sizeof(md_mtl_page_t));
}

static bool md_mtl_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host,
                               id<MTLBuffer>* out_buf, uint64_t* out_off) {
    md_gpu_device_t dev = s->device;
    md_mtl_arena_t*  a   = &s->arena;
    uint64_t need = md_mtl_align_up(size, MD_MTL_ARG_ALIGN);

    md_mtl_page_t* chosen = NULL;
    if (a->pages.count > 0) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, a->current);
        if (p->cursor + need <= p->capacity) chosen = p;
    }
    if (!chosen) {
        uint64_t done = md_mtl_stream_completed(s);
        for (size_t i = 0; i < a->pages.count && !chosen; ++i) {
            md_mtl_page_t* p = MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, i);
            if (p->retire_value != 0 && p->retire_value <= done) { p->cursor = 0; p->retire_value = 0; }
            if (p->retire_value == 0 && p->cursor + need <= p->capacity) { a->current = i; chosen = p; }
        }
    }
    if (!chosen) {
        uint64_t page_size = MD_MTL_ARENA_PAGE_SIZE;
        while (page_size < need) page_size *= 2;
        chosen = md_mtl_page_create(dev, page_size);
        if (!chosen) return md_mtl_fail("failed to allocate a transient page");
        md_mtl_page_t** slot = (md_mtl_page_t**)md_mtl_vec_push(&a->pages, dev->alloc);
        if (!slot) { md_mtl_page_destroy(dev, chosen); return false; }
        *slot = chosen;
        a->current = a->pages.count - 1;
    }

    *out_addr = chosen->address + chosen->cursor;
    *out_host = chosen->host + chosen->cursor;
    if (out_buf) *out_buf = chosen->buffer;
    if (out_off) *out_off = chosen->cursor;
    chosen->cursor += need;
    return true;
}

/* Stamp every page holding data with the value that releases it. Must
   overwrite an existing stamp -- see the note in md_gpu_vulkan.c: a page
   appended to across several submissions would otherwise be recycled when the
   first of them completed, while a later one was still reading it. */
static void md_mtl_arena_retire(md_gpu_stream_t s, uint64_t value) {
    md_mtl_arena_t* a = &s->arena;
    for (size_t i = 0; i < a->pages.count; ++i) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, i);
        if (p->cursor > 0) p->retire_value = value;
    }
}

static void md_mtl_arena_free(md_gpu_device_t dev, md_mtl_arena_t* a) {
    for (size_t i = 0; i < a->pages.count; ++i) {
        md_mtl_page_destroy(dev, MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, i));
    }
    md_mtl_vec_free(&a->pages, dev->alloc);
}

/* =========================================================================
   5. Encoders and program order
   ========================================================================= */

static bool md_mtl_stream_ensure_cmd(md_gpu_stream_t s) {
    if (s->cmd) return true;
    id<MTLCommandBuffer> cb = [s->queue commandBuffer];
    if (!cb) return md_mtl_fail("commandBuffer failed on stream '%s'", s->label);
    MD_MTL_RETAIN(cb);
    cb.label = [NSString stringWithUTF8String:s->label];
    s->cmd      = cb;
    s->has_work = false;
    /* Unlike Vulkan, Metal executes command buffers on a queue in the order
       they were committed, and each stream owns its queue, so work in a new
       command buffer is already ordered after the previous one. No barrier is
       needed at the boundary -- only within an encoder. */
    s->needs_barrier = false;

    /* Pending cross-stream waits are encoded at the head of the buffer. */
    for (uint32_t i = 0; i < s->wait_count; ++i) {
        [cb encodeWaitForEvent:s->wait_events[i] value:s->wait_values[i]];
    }
    s->wait_count = 0;
    return true;
}

static void md_mtl_end_encoders(md_gpu_stream_t s) {
    if (s->compute_enc) { [s->compute_enc endEncoding]; MD_MTL_RELEASE(s->compute_enc); s->compute_enc = nil; }
    if (s->blit_enc)    { [s->blit_enc    endEncoding]; MD_MTL_RELEASE(s->blit_enc);    s->blit_enc    = nil; }
}

/* Bind the bindless texture table on every compute encoder we open. */
static id<MTLComputeCommandEncoder> md_mtl_compute_encoder(md_gpu_stream_t s) {
    if (!md_mtl_stream_ensure_cmd(s)) return nil;
    if (s->compute_enc) return s->compute_enc;
    if (s->blit_enc) { [s->blit_enc endEncoding]; MD_MTL_RELEASE(s->blit_enc); s->blit_enc = nil; }

    id<MTLComputeCommandEncoder> enc = [s->cmd computeCommandEncoder];
    if (!enc) { md_mtl_fail("computeCommandEncoder failed"); return nil; }
    MD_MTL_RETAIN(enc);
    md_gpu_device_t dev = s->device;
    if (!dev->has_residency_set && dev->live_bufs.count > 0) {
        /* Fallback residency: declare everything once per encoder. */
        [enc useResources:(__unsafe_unretained id<MTLResource>*)dev->live_bufs.data
                    count:dev->live_bufs.count
                    usage:MTLResourceUsageRead | MTLResourceUsageWrite];
    }
    s->compute_enc = enc;
    return enc;
}

static id<MTLBlitCommandEncoder> md_mtl_blit_encoder(md_gpu_stream_t s) {
    if (!md_mtl_stream_ensure_cmd(s)) return nil;
    if (s->blit_enc) return s->blit_enc;
    if (s->compute_enc) { [s->compute_enc endEncoding]; MD_MTL_RELEASE(s->compute_enc); s->compute_enc = nil; }

    id<MTLBlitCommandEncoder> enc = [s->cmd blitCommandEncoder];
    if (!enc) { md_mtl_fail("blitCommandEncoder failed"); return nil; }
    MD_MTL_RETAIN(enc);
    s->blit_enc = enc;
    return enc;
}

/* Program order. Within a compute encoder this is an explicit memory barrier;
   across encoders Metal already orders the work, so nothing is needed. */
static void md_mtl_barrier_if_needed(md_gpu_stream_t s) {
    if (!s->needs_barrier || !s->auto_order) return;
    if (s->compute_enc) {
        [s->compute_enc memoryBarrierWithScope:MTLBarrierScopeBuffers | MTLBarrierScopeTextures];
    }
    s->needs_barrier = false;
}

void md_gpu_stream_barrier(md_gpu_stream_t s) {
    if (!s) return;
    bool saved = s->auto_order;
    s->auto_order    = true;
    s->needs_barrier = true;
    md_mtl_barrier_if_needed(s);
    s->auto_order = saved;
}

void md_gpu_stream_set_auto_order(md_gpu_stream_t s, bool enabled) { if (s) s->auto_order = enabled; }

static void md_mtl_did_op(md_gpu_stream_t s) {
    s->needs_barrier = true;
    if (!s->capture) s->has_work = true;
}

/* =========================================================================
   6. Streams
   ========================================================================= */

static uint64_t md_mtl_stream_completed(md_gpu_stream_t s) {
    return (uint64_t)s->timeline.signaledValue;
}

static md_gpu_stream_t md_mtl_stream_create_internal(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label, bool is_default) {
    md_gpu_stream_t s = (md_gpu_stream_t)md_alloc(dev->alloc, sizeof(md_gpu_stream));
    if (!s) { md_mtl_fail("out of memory"); return NULL; }
    memset(s, 0, sizeof(*s));
    s->device     = dev;
    s->kind       = kind;
    s->auto_order = true;
    s->next_value = 1;
    s->is_default = is_default;
    snprintf(s->label, sizeof(s->label), "%s", label ? label : "stream");

    /* Its own hardware queue, so two streams can genuinely run at once. */
    id<MTLCommandQueue> q = [dev->device newCommandQueue];
    if (!q) { md_free(dev->alloc, s, sizeof(*s)); md_mtl_fail("newCommandQueue failed"); return NULL; }
    MD_MTL_RETAIN(q);
    q.label = [NSString stringWithUTF8String:s->label];
    s->queue     = q;
    s->owns_queue = true;
    if (dev->has_residency_set) {
        [q performSelector:@selector(addResidencySet:) withObject:dev->residency_set];
    }

    md_mtl_vec_init(&s->arena.pages, sizeof(md_mtl_page_t*));

    id<MTLSharedEvent> ev = [dev->device newSharedEvent];
    if (!ev) { md_free(dev->alloc, s, sizeof(*s)); md_mtl_fail("newSharedEvent failed"); return NULL; }
    MD_MTL_RETAIN(ev);
    s->timeline = ev;

    md_mutex_lock(&dev->device_mutex);
    md_gpu_stream_t* slot = (md_gpu_stream_t*)md_mtl_vec_push(&dev->streams, dev->alloc);
    if (slot) *slot = s;
    md_mutex_unlock(&dev->device_mutex);
    return s;
}

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label) {
    if (!dev) return NULL;
    return md_mtl_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
    for (size_t i = 0; i < dev->streams.count; ++i) {
        if (arr[i]->is_default && arr[i]->kind == kind) return arr[i];
    }
    return NULL;
}

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t s) { return s ? s->device : NULL; }

static bool md_mtl_stream_submit(md_gpu_stream_t s) {
    if (!s->cmd || !s->has_work) return true;
    md_mtl_end_encoders(s);

    uint64_t signal_value = s->next_value;
    [s->cmd encodeSignalEvent:s->timeline value:signal_value];

    md_mutex_lock(&s->device->queue_mutex);
    [s->cmd commit];
    md_mutex_unlock(&s->device->queue_mutex);

    MD_MTL_RELEASE(s->cmd);
    s->cmd = nil;

    md_mtl_arena_retire(s, signal_value);
    s->submitted_value = signal_value;
    s->next_value      = signal_value + 1;
    s->has_work        = false;
    s->needs_barrier   = false;
    return true;
}

md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t s) {
    md_gpu_sync_t out = md_gpu_sync_none();
    if (!s) return out;
    if (s->capture) { md_mtl_fail("md_gpu_stream_record is not allowed during capture"); return out; }
    md_mtl_stream_submit(s);
    if (s->submitted_value == 0) return out;
    out.stream = s;
    out.value  = s->submitted_value;
    return out;
}

void md_gpu_stream_wait(md_gpu_stream_t s, md_gpu_sync_t sync) {
    if (!s || !md_gpu_sync_is_valid(sync)) return;
    if (sync.stream == s) return;
    if (md_gpu_sync_is_complete(sync)) return;

    if (s->has_work) md_mtl_stream_submit(s);

    for (uint32_t i = 0; i < s->wait_count; ++i) {
        if (s->wait_events[i] == sync.stream->timeline) {
            if (sync.value > s->wait_values[i]) s->wait_values[i] = sync.value;
            return;
        }
    }
    if (s->wait_count >= MD_MTL_MAX_PENDING_WAITS) {
        md_mtl_fail("more than %u pending waits on stream '%s'", MD_MTL_MAX_PENDING_WAITS, s->label);
        return;
    }
    s->wait_events[s->wait_count] = sync.stream->timeline;
    s->wait_values[s->wait_count] = sync.value;
    s->wait_count++;
}

void md_gpu_stream_flush(md_gpu_stream_t s) {
    if (!s || s->capture) return;
    md_mtl_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    if (s->capture) { md_mtl_fail("md_gpu_stream_sync is not allowed during capture"); return; }
    md_mtl_stream_submit(s);
    while (s->submitted_value > 0 && md_mtl_stream_completed(s) < s->submitted_value) {
        /* MTLSharedEvent has no blocking wait in C; spin with a short yield.
           Submissions are already in flight, so this is bounded. */
        md_thread_sleep(0);
    }
    md_mtl_poll_internal(s->device, false);
}

bool md_gpu_sync_is_complete(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    return md_mtl_stream_completed(sync.stream) >= sync.value;
}

void md_gpu_sync_wait(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return;
    while (md_mtl_stream_completed(sync.stream) < sync.value) md_thread_sleep(0);
}

static void md_mtl_stream_destroy_internal(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    md_gpu_stream_sync(s);
    md_mtl_end_encoders(s);
    if (s->cmd) { MD_MTL_RELEASE(s->cmd); s->cmd = nil; }
    md_mtl_arena_free(dev, &s->arena);
    if (s->owns_queue) MD_MTL_RELEASE(s->queue);
    MD_MTL_RELEASE(s->timeline);
    md_free(dev->alloc, s, sizeof(*s));
}

void md_gpu_stream_destroy(md_gpu_stream_t s) {
    if (!s || s->is_default) return;
    md_gpu_device_t dev = s->device;
    md_mutex_lock(&dev->device_mutex);
    md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
    for (size_t i = 0; i < dev->streams.count; ++i) {
        if (arr[i] == s) {
            memmove(&arr[i], &arr[i + 1], (dev->streams.count - i - 1) * sizeof(*arr));
            dev->streams.count--;
            break;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    md_mtl_stream_destroy_internal(s);
}

/* =========================================================================
   7. Memory
   ========================================================================= */

md_gpu_pool_t md_gpu_pool_create(md_gpu_device_t dev, const md_gpu_pool_desc_t* desc) {
    if (!dev || !desc) return NULL;
    const uint64_t release_threshold = desc->release_threshold;
    const char* label = desc->label;
    md_gpu_pool_t p = (md_gpu_pool_t)md_alloc(dev->alloc, sizeof(md_gpu_pool));
    if (!p) { md_mtl_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device            = dev;
    p->release_threshold = release_threshold;
    p->flags             = desc->flags;
    p->block_size        = desc->block_size;
    snprintf(p->label, sizeof(p->label), "%s", label ? label : "pool");
    md_mtl_vec_init(&p->blocks, sizeof(md_mtl_block_t*));

    md_mutex_lock(&dev->device_mutex);
    md_gpu_pool_t* slot = (md_gpu_pool_t*)md_mtl_vec_push(&dev->pools, dev->alloc);
    if (slot) *slot = p;
    md_mutex_unlock(&dev->device_mutex);
    return p;
}

void md_gpu_pool_stats(md_gpu_pool_t p, md_gpu_pool_stats_t* out) {
    if (!p || !out) return;
    memset(out, 0, sizeof(*out));
    md_mutex_lock(&p->device->device_mutex);
    out->bytes_in_use      = p->in_use_bytes;
    out->bytes_reserved    = p->reserved_bytes;
    out->bytes_cached      = p->reserved_bytes - p->in_use_bytes;
    out->bytes_peak_in_use = p->peak_in_use_bytes;
    out->alloc_count       = p->alloc_count;
    out->reuse_count       = p->reuse_count;
    md_mtl_block_t** arr = (md_mtl_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (arr[i]->in_use) out->blocks_in_use++;
        else                out->blocks_cached++;
    }
    md_mutex_unlock(&p->device->device_mutex);
}

/* Mark a block free at the current point in `stream`. Shared by
   md_gpu_free and md_gpu_pool_reset. Caller holds device_mutex. */
static void md_mtl_block_release(md_mtl_block_t* b, md_gpu_stream_t stream) {
    b->in_use      = false;
    b->free_stream = stream;
    b->free_value  = stream->has_work ? stream->next_value : stream->submitted_value;
    if (b->free_value == 0) b->free_stream = NULL;
    b->pool->in_use_bytes -= b->capacity;
}

void md_gpu_pool_reset(md_gpu_pool_t p, md_gpu_stream_t stream) {
    if (!p || !stream) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_block_t** arr = (md_mtl_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (arr[i]->in_use) md_mtl_block_release(arr[i], stream);
    }
    md_mutex_unlock(&dev->device_mutex);
}

static void md_mtl_block_destroy(md_gpu_pool_t p, md_mtl_block_t* b) {
    md_gpu_device_t dev = p->device;
    md_mtl_registry_remove(dev, b);
    /* Drop it from the fallback residency list too. */
    id<MTLBuffer>* arr = (id<MTLBuffer>*)dev->live_bufs.data;
    for (size_t i = 0; i < dev->live_bufs.count; ++i) {
        if (arr[i] == b->buffer) {
            memmove(&arr[i], &arr[i + 1], (dev->live_bufs.count - i - 1) * sizeof(*arr));
            dev->live_bufs.count--;
            break;
        }
    }
    md_mtl_destroy_raw_buffer(dev, b->buffer);
    p->reserved_bytes -= b->capacity;
    md_free(dev->alloc, b, sizeof(*b));
}

static void md_mtl_pool_trim_locked(md_gpu_pool_t p, uint64_t keep_bytes) {
    md_mtl_block_t** arr = (md_mtl_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count;) {
        md_mtl_block_t* b = arr[i];
        bool idle = !b->in_use &&
            (b->free_stream == NULL || md_mtl_stream_completed(b->free_stream) >= b->free_value);
        if (idle && p->reserved_bytes - p->in_use_bytes > keep_bytes) {
            md_mtl_block_destroy(p, b);
            memmove(&arr[i], &arr[i + 1], (p->blocks.count - i - 1) * sizeof(*arr));
            p->blocks.count--;
            continue;
        }
        ++i;
    }
}

void md_gpu_pool_trim(md_gpu_pool_t p, uint64_t keep_bytes) {
    if (!p) return;
    md_mutex_lock(&p->device->device_mutex);
    md_mtl_pool_trim_locked(p, keep_bytes);
    md_mutex_unlock(&p->device->device_mutex);
}

void md_gpu_pool_destroy(md_gpu_pool_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_block_t** arr = (md_mtl_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) md_mtl_block_destroy(p, arr[i]);
    md_mtl_vec_free(&p->blocks, dev->alloc);

    md_gpu_pool_t* parr = (md_gpu_pool_t*)dev->pools.data;
    for (size_t i = 0; i < dev->pools.count; ++i) {
        if (parr[i] == p) {
            memmove(&parr[i], &parr[i + 1], (dev->pools.count - i - 1) * sizeof(*parr));
            dev->pools.count--;
            break;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, p, sizeof(*p));
}

md_gpu_ptr_t md_gpu_malloc(md_gpu_pool_t p, size_t size, md_gpu_stream_t stream) {
    if (!p || size == 0) return NULL;
    const md_gpu_mem_flags_t flags = p->flags;
    md_gpu_device_t dev = p->device;
    uint64_t capacity = md_mtl_next_pow2(size);

    md_mutex_lock(&dev->device_mutex);

    md_mtl_block_t** arr = (md_mtl_block_t**)p->blocks.data;
    md_mtl_block_t* best = NULL;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_mtl_block_t* b = arr[i];
        if (b->in_use || b->flags != flags || b->capacity < size) continue;
        bool safe = (b->free_stream == NULL) || (b->free_stream == stream)
                 || (md_mtl_stream_completed(b->free_stream) >= b->free_value);
        if (!safe) continue;
        if (!best || b->capacity < best->capacity) best = b;
    }
    if (best) {
        best->in_use      = true;
        best->size        = size;
        best->free_stream = NULL;
        best->free_value  = 0;
        p->in_use_bytes  += best->capacity;
        if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
        p->alloc_count++;
        p->reuse_count++;
        md_mutex_unlock(&dev->device_mutex);
        return (md_gpu_ptr_t)(uintptr_t)best->address;
    }

    md_mtl_block_t* b = (md_mtl_block_t*)md_alloc(dev->alloc, sizeof(md_mtl_block_t));
    if (!b) { md_mutex_unlock(&dev->device_mutex); md_mtl_fail("out of memory"); return NULL; }
    memset(b, 0, sizeof(*b));
    id<MTLBuffer> buf = nil;
    if (!md_mtl_create_raw_buffer(dev, capacity, flags, &buf, &b->address, &b->host)) {
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        return NULL;
    }
    b->buffer   = buf;
    b->capacity = capacity;
    b->size     = size;
    b->flags    = flags;
    b->pool     = p;
    b->in_use   = true;

    md_mtl_block_t** slot = (md_mtl_block_t**)md_mtl_vec_push(&p->blocks, dev->alloc);
    if (!slot || !md_mtl_registry_insert(dev, b)) {
        md_mtl_destroy_raw_buffer(dev, buf);
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        return NULL;
    }
    *slot = b;
    if (!dev->has_residency_set) {
        id<MTLBuffer>* lslot = (id<MTLBuffer>*)md_mtl_vec_push(&dev->live_bufs, dev->alloc);
        if (lslot) *lslot = buf;
    }
    p->reserved_bytes += capacity;
    p->in_use_bytes   += capacity;
    if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
    p->alloc_count++;

    md_mutex_unlock(&dev->device_mutex);
    return (md_gpu_ptr_t)(uintptr_t)b->address;
}

void md_gpu_free(md_gpu_ptr_t ptr, md_gpu_stream_t stream) {
    if (!ptr || !stream) return;
    md_gpu_device_t dev = stream->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_block_t* b = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)ptr);
    if (!b || !b->in_use) { md_mutex_unlock(&dev->device_mutex); return; }
    md_mtl_block_release(b, stream);
    if (b->pool->release_threshold == 0) md_mtl_pool_trim_locked(b->pool, 0);
    md_mutex_unlock(&dev->device_mutex);
}

void* md_gpu_host_ptr(md_gpu_ptr_t ptr) {
    if (!ptr) return NULL;
    uint64_t addr = (uint64_t)(uintptr_t)ptr;
    md_mtl_block_t* b = md_mtl_find_anywhere(addr);
    if (!b || !b->host) return NULL;
    return (uint8_t*)b->host + (addr - b->address);
}

size_t md_gpu_ptr_size(md_gpu_ptr_t ptr) {
    if (!ptr) return 0;
    uint64_t addr = (uint64_t)(uintptr_t)ptr;
    md_mtl_block_t* b = md_mtl_find_anywhere(addr);
    if (!b || !b->in_use) return 0;
    uint64_t off = addr - b->address;
    return (size_t)(b->size > off ? b->size - off : 0);
}

md_gpu_ptr_t md_gpu_ptr_base(md_gpu_ptr_t ptr) {
    if (!ptr) return NULL;
    md_mtl_block_t* b = md_mtl_find_anywhere((uint64_t)(uintptr_t)ptr);
    return b ? (md_gpu_ptr_t)(uintptr_t)b->address : NULL;
}

/* Record or execute a buffer copy. */
static bool md_mtl_emit_copy(md_gpu_stream_t s, id<MTLBuffer> src, uint64_t src_off,
                             id<MTLBuffer> dst, uint64_t dst_off, uint64_t size) {
    if (s->capture) {
        md_mtl_op_t* op = (md_mtl_op_t*)md_mtl_vec_push(&s->capture->ops, s->device->alloc);
        if (!op) return false;
        op->kind = MD_MTL_OP_COPY_BUFFER;
        op->src_buffer = src; op->src_offset = src_off;
        op->dst_buffer = dst; op->dst_offset = dst_off;
        op->size = size;
        md_mtl_did_op(s);
        return true;
    }
    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    [enc copyFromBuffer:src sourceOffset:src_off toBuffer:dst destinationOffset:dst_off size:size];
    md_mtl_did_op(s);
    return true;
}

typedef struct md_mtl_d2h_t {
    md_gpu_device_t dev;
    void*            dst;
    const void*      staging_host;
    size_t           size;
} md_mtl_d2h_t;

static void md_mtl_d2h_finish(void* user) {
    md_mtl_d2h_t* c = (md_mtl_d2h_t*)user;
    memcpy(c->dst, c->staging_host, c->size);
    md_free(c->dev->alloc, c, sizeof(*c));
}

bool md_gpu_memcpy_async(void* dst, const void* src, size_t size, md_gpu_stream_t s) {
    if (!s) return md_mtl_fail("md_gpu_memcpy_async requires a stream");
    if (size == 0) return true;
    md_gpu_device_t dev = s->device;

    md_mtl_block_t* dblk = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)dst);
    md_mtl_block_t* sblk = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)src);

    if (dblk && sblk) {
        return md_mtl_emit_copy(s, sblk->buffer, (uint64_t)(uintptr_t)src - sblk->address,
                                   dblk->buffer, (uint64_t)(uintptr_t)dst - dblk->address, size);
    }

    if (dblk && !sblk) {
        uint64_t dst_off = (uint64_t)(uintptr_t)dst - dblk->address;
        if (dblk->host && !s->capture && !s->has_work &&
            md_mtl_stream_completed(s) >= s->submitted_value) {
            memcpy((uint8_t*)dblk->host + dst_off, src, size);
            return true;
        }
        uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
        if (!md_mtl_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
        memcpy(host, src, size);
        return md_mtl_emit_copy(s, buf, off, dblk->buffer, dst_off, size);
    }

    if (!dblk && sblk) {
        uint64_t src_off = (uint64_t)(uintptr_t)src - sblk->address;
        if (sblk->host && !s->has_work && md_mtl_stream_completed(s) >= s->submitted_value) {
            memcpy(dst, (uint8_t*)sblk->host + src_off, size);
            return true;
        }
        uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
        if (!md_mtl_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
        if (!md_mtl_emit_copy(s, sblk->buffer, src_off, buf, off, size)) return false;

        md_mtl_d2h_t* c = (md_mtl_d2h_t*)md_alloc(dev->alloc, sizeof(md_mtl_d2h_t));
        if (!c) return md_mtl_fail("out of memory");
        c->dev = dev; c->dst = dst; c->staging_host = host; c->size = size;

        md_gpu_sync_t sync = md_gpu_stream_record(s);
        md_mutex_lock(&dev->device_mutex);
        md_mtl_hostfn_t* h = (md_mtl_hostfn_t*)md_mtl_vec_push(&dev->hostfns, dev->alloc);
        if (h) { h->sync = sync; h->fn = md_mtl_d2h_finish; h->user = c; h->internal = true; }
        md_mutex_unlock(&dev->device_mutex);
        return h != NULL;
    }

    memcpy(dst, src, size);
    return true;
}

bool md_gpu_memset_async(md_gpu_ptr_t dst, uint8_t value, size_t size, md_gpu_stream_t s) {
    if (!s) return md_mtl_fail("md_gpu_memset_async requires a stream");
    if (!dst || size == 0) return true;
    md_mtl_block_t* b = md_mtl_registry_find(s->device, (uint64_t)(uintptr_t)dst);
    if (!b) return md_mtl_fail("md_gpu_memset_async: not a device pointer");
    uint64_t off = (uint64_t)(uintptr_t)dst - b->address;

    /* Metal's fillBuffer works at byte granularity, so there is no aligned /
       unaligned split to do here (unlike vkCmdFillBuffer). */
    if (s->capture) {
        md_mtl_op_t* op = (md_mtl_op_t*)md_mtl_vec_push(&s->capture->ops, s->device->alloc);
        if (!op) return false;
        op->kind = MD_MTL_OP_FILL_BUFFER;
        op->dst_buffer = b->buffer;
        op->dst_offset = off;
        op->size = size;
        op->fill_value = value;
        md_mtl_did_op(s);
        return true;
    }
    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    [enc fillBuffer:b->buffer range:NSMakeRange((NSUInteger)off, (NSUInteger)size) value:value];
    md_mtl_did_op(s);
    return true;
}

void* md_gpu_upload_begin(md_gpu_stream_t s, md_gpu_ptr_t dst, size_t size) {
    if (!s || !dst || size == 0) return NULL;
    if (s->upload_open) { md_mtl_fail("an upload is already open on stream '%s'", s->label); return NULL; }
    md_mtl_block_t* b = md_mtl_registry_find(s->device, (uint64_t)(uintptr_t)dst);
    if (!b) { md_mtl_fail("md_gpu_upload_begin: not a device pointer"); return NULL; }
    uint64_t dst_off = (uint64_t)(uintptr_t)dst - b->address;

    if (b->host && !s->capture && !s->has_work && md_mtl_stream_completed(s) >= s->submitted_value) {
        s->upload_open = true; s->upload_direct = true;
        s->upload_dst = dst; s->upload_size = size;
        return (uint8_t*)b->host + dst_off;
    }
    uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
    if (!md_mtl_arena_alloc(s, size, &addr, &host, &buf, &off)) return NULL;
    s->upload_open     = true;
    s->upload_direct   = false;
    s->upload_dst      = dst;
    s->upload_src_addr = addr;
    s->upload_size     = size;
    return host;
}

bool md_gpu_upload_end(md_gpu_stream_t s) {
    if (!s || !s->upload_open) return md_mtl_fail("no upload is open");
    s->upload_open = false;
    if (s->upload_direct) return true;

    md_gpu_device_t dev = s->device;
    md_mtl_block_t* b = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)s->upload_dst);
    if (!b) return md_mtl_fail("upload destination is no longer live");
    uint64_t dst_off = (uint64_t)(uintptr_t)s->upload_dst - b->address;

    md_mtl_page_t* page = NULL; uint64_t page_off = 0;
    for (size_t i = 0; i < s->arena.pages.count; ++i) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(s->arena.pages, md_mtl_page_t*, i);
        if (s->upload_src_addr >= p->address && s->upload_src_addr < p->address + p->capacity) {
            page = p; page_off = s->upload_src_addr - p->address; break;
        }
    }
    if (!page) return md_mtl_fail("staging address not found in arena");
    return md_mtl_emit_copy(s, page->buffer, page_off, b->buffer, dst_off, s->upload_size);
}

/* =========================================================================
   8. Textures and samplers
   ========================================================================= */

typedef struct { MTLPixelFormat fmt; uint32_t bytes; } md_mtl_fmt_t;

static md_mtl_fmt_t md_mtl_format_info(md_gpu_format_t f) {
    md_mtl_fmt_t i = {MTLPixelFormatInvalid, 0};
    switch (f) {
    case MD_GPU_FORMAT_R32_FLOAT:    i.fmt = MTLPixelFormatR32Float;    i.bytes = 4;  break;
    case MD_GPU_FORMAT_R32_UINT:     i.fmt = MTLPixelFormatR32Uint;     i.bytes = 4;  break;
    case MD_GPU_FORMAT_RGBA32_FLOAT: i.fmt = MTLPixelFormatRGBA32Float; i.bytes = 16; break;
    case MD_GPU_FORMAT_RGBA8_UNORM:  i.fmt = MTLPixelFormatRGBA8Unorm;  i.bytes = 4;  break;
    default: break;
    }
    return i;
}

/* Resolve a handle back to its record. Both tables stay small. */
static md_gpu_device_t md_mtl_tex_device(md_gpu_tex_t tex, md_mtl_tex_t** out) {
    if (!md_mtl_tex_valid(tex)) return NULL;
    for (uint32_t i = 0; i < md_mtl_device_count; ++i) {
        md_gpu_device_t d = md_mtl_devices[i];
        if (!d) continue;
        md_mtl_tex_t** arr = (md_mtl_tex_t**)d->textures.data;
        for (size_t j = 0; j < d->textures.count; ++j) {
            if (arr[j]->handle == tex) { if (out) *out = arr[j]; return d; }
        }
    }
    return NULL;
}

md_gpu_tex_t md_gpu_tex_create(md_gpu_device_t dev, const md_gpu_tex_desc_t* desc) {
    if (!dev || !desc) return 0;
    md_mtl_fmt_t fi = md_mtl_format_info(desc->format);
    if (fi.fmt == MTLPixelFormatInvalid) { md_mtl_fail("unsupported texture format %d", (int)desc->format); return 0; }
    if (!(desc->flags & (MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED))) {
        md_mtl_fail("texture needs at least one of MD_GPU_TEX_STORAGE / MD_GPU_TEX_SAMPLED");
        return 0;
    }

    bool is_3d = desc->depth > 1;

    MTLTextureDescriptor* td = [[MTLTextureDescriptor alloc] init];
    td.textureType      = is_3d ? MTLTextureType3D : MTLTextureType2D;
    td.pixelFormat      = fi.fmt;
    td.width            = desc->width  ? desc->width  : 1;
    td.height           = desc->height ? desc->height : 1;
    td.depth            = is_3d ? desc->depth : 1;
    td.mipmapLevelCount = 1;
    td.storageMode      = MTLStorageModePrivate;
    td.usage            = 0;
    if (desc->flags & MD_GPU_TEX_STORAGE) td.usage |= MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;
    if (desc->flags & MD_GPU_TEX_SAMPLED) td.usage |= MTLTextureUsageShaderRead;

    id<MTLTexture> tex = [dev->device newTextureWithDescriptor:td];
    if (!tex) { md_mtl_fail("newTextureWithDescriptor failed"); return 0; }
    MD_MTL_RETAIN(tex);
    md_mtl_make_resident(dev, tex);

    md_mtl_tex_t* rec = (md_mtl_tex_t*)md_alloc(dev->alloc, sizeof(md_mtl_tex_t));
    if (!rec) { MD_MTL_RELEASE(tex); md_mtl_fail("out of memory"); return 0; }
    memset(rec, 0, sizeof(*rec));
    rec->texture = tex;
    rec->desc    = *desc;
    rec->handle  = (uint64_t)tex.gpuResourceID._impl;

    md_mutex_lock(&dev->device_mutex);
    md_mtl_tex_t** slot = (md_mtl_tex_t**)md_mtl_vec_push(&dev->textures, dev->alloc);
    if (slot) *slot = rec;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { MD_MTL_RELEASE(tex); md_free(dev->alloc, rec, sizeof(*rec)); return 0; }

    return rec->handle;
}

/* On Metal one texture serves both access::read_write and access::sample, so
   there is no second handle -- unlike Vulkan, where a mutable descriptor slot
   holds one descriptor type at a time. */
md_gpu_tex_t md_gpu_tex_sampled(md_gpu_tex_t tex) {
    return tex;
}

bool md_gpu_tex_desc(md_gpu_tex_t tex, md_gpu_tex_desc_t* out) {
    md_mtl_tex_t* t = NULL;
    if (!out || !md_mtl_tex_device(tex, &t) || !t) return false;
    *out = t->desc;
    return true;
}

void md_gpu_tex_destroy(md_gpu_tex_t tex, md_gpu_stream_t stream) {
    md_mtl_tex_t* t = NULL;
    md_gpu_device_t dev = md_mtl_tex_device(tex, &t);
    if (!dev || !t) return;
    (void)stream;

    md_mutex_lock(&dev->device_mutex);
    md_mtl_retire_t* r = (md_mtl_retire_t*)md_mtl_vec_push(&dev->retires, dev->alloc);
    if (r) {
        r->texture    = t->texture;
        r->wait_count = 0;
        md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
        for (size_t i = 0; i < dev->streams.count && r->wait_count < 8; ++i) {
            md_gpu_stream_t s = arr[i];
            uint64_t v = s->has_work ? s->next_value : s->submitted_value;
            if (v > 0) {
                r->waits[r->wait_count].stream = s;
                r->waits[r->wait_count].value  = v;
                r->wait_count++;
            }
        }
    }
    md_mtl_tex_t** arr = (md_mtl_tex_t**)dev->textures.data;
    for (size_t i = 0; i < dev->textures.count; ++i) {
        if (arr[i] == t) {
            memmove(&arr[i], &arr[i + 1], (dev->textures.count - i - 1) * sizeof(*arr));
            dev->textures.count--;
            break;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, t, sizeof(*t));
}

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t dev, const md_gpu_sampler_desc_t* desc) {
    if (!dev) return 0;
    md_gpu_sampler_desc_t d = {0};
    if (desc) d = *desc;

    static const MTLSamplerAddressMode modes[] = {
        MTLSamplerAddressModeClampToEdge,
        MTLSamplerAddressModeRepeat,
        MTLSamplerAddressModeMirrorRepeat,
    };

    MTLSamplerDescriptor* sd = [[MTLSamplerDescriptor alloc] init];
    sd.minFilter    = d.min_filter == MD_GPU_FILTER_LINEAR ? MTLSamplerMinMagFilterLinear : MTLSamplerMinMagFilterNearest;
    sd.magFilter    = d.mag_filter == MD_GPU_FILTER_LINEAR ? MTLSamplerMinMagFilterLinear : MTLSamplerMinMagFilterNearest;
    sd.sAddressMode = modes[d.address_u % 3];
    sd.tAddressMode = modes[d.address_v % 3];
    sd.rAddressMode = modes[d.address_w % 3];
    /* Required for gpuResourceID to be usable from an argument buffer. */
    sd.supportArgumentBuffers = YES;

    id<MTLSamplerState> smp = [dev->device newSamplerStateWithDescriptor:sd];
    if (!smp) { md_mtl_fail("newSamplerStateWithDescriptor failed"); return 0; }
    MD_MTL_RETAIN(smp);

    md_mtl_sampler_t* rec = (md_mtl_sampler_t*)md_alloc(dev->alloc, sizeof(md_mtl_sampler_t));
    if (!rec) { MD_MTL_RELEASE(smp); md_mtl_fail("out of memory"); return 0; }
    memset(rec, 0, sizeof(*rec));
    rec->sampler = smp;
    rec->handle  = (uint64_t)smp.gpuResourceID._impl;

    md_mutex_lock(&dev->device_mutex);
    md_mtl_sampler_t** slot = (md_mtl_sampler_t**)md_mtl_vec_push(&dev->samplers, dev->alloc);
    if (slot) *slot = rec;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { MD_MTL_RELEASE(smp); md_free(dev->alloc, rec, sizeof(*rec)); return 0; }

    return rec->handle;
}

void md_gpu_sampler_destroy(md_gpu_sampler_t sampler) {
    if (!md_mtl_tex_valid(sampler)) return;
    for (uint32_t i = 0; i < md_mtl_device_count; ++i) {
        md_gpu_device_t d = md_mtl_devices[i];
        if (!d) continue;
        md_mutex_lock(&d->device_mutex);
        md_mtl_sampler_t** arr = (md_mtl_sampler_t**)d->samplers.data;
        for (size_t j = 0; j < d->samplers.count; ++j) {
            if (arr[j]->handle != sampler) continue;
            md_mtl_sampler_t* rec = arr[j];
            memmove(&arr[j], &arr[j + 1], (d->samplers.count - j - 1) * sizeof(*arr));
            d->samplers.count--;
            md_mutex_unlock(&d->device_mutex);
            MD_MTL_RELEASE(rec->sampler);
            md_free(d->alloc, rec, sizeof(*rec));
            return;
        }
        md_mutex_unlock(&d->device_mutex);
    }
}

static void md_mtl_resolve_region(const md_gpu_tex_desc_t* d, const md_gpu_tex_region_t* r,
                                  MTLOrigin* origin, MTLSize* size) {
    uint32_t w = d->width ? d->width : 1;
    uint32_t h = d->height ? d->height : 1;
    uint32_t z = d->depth ? d->depth : 1;
    if (r) {
        *origin = MTLOriginMake(r->offset[0], r->offset[1], r->offset[2]);
        *size   = MTLSizeMake(r->extent[0] ? r->extent[0] : w - r->offset[0],
                              r->extent[1] ? r->extent[1] : h - r->offset[1],
                              r->extent[2] ? r->extent[2] : z - r->offset[2]);
    } else {
        *origin = MTLOriginMake(0, 0, 0);
        *size   = MTLSizeMake(w, h, z);
    }
}

bool md_gpu_memcpy_to_tex_async(md_gpu_tex_t tex, const md_gpu_tex_region_t* region,
                                 const void* src, size_t size, md_gpu_stream_t s) {
    if (!s) return md_mtl_fail("md_gpu_memcpy_to_tex_async requires a stream");
    md_mtl_tex_t* t = NULL;
    md_gpu_device_t dev = md_mtl_tex_device(tex, &t);
    if (!dev || !t) return md_mtl_fail("invalid texture handle");

    MTLOrigin origin; MTLSize extent;
    md_mtl_resolve_region(&t->desc, region, &origin, &extent);
    md_mtl_fmt_t fi = md_mtl_format_info(t->desc.format);
    uint64_t bpr = (uint64_t)extent.width * fi.bytes;
    uint64_t bpi = bpr * extent.height;

    id<MTLBuffer> src_buf = nil; uint64_t src_off = 0;
    md_mtl_block_t* sblk = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)src);
    if (sblk) {
        src_buf = sblk->buffer;
        src_off = (uint64_t)(uintptr_t)src - sblk->address;
    } else {
        uint64_t addr; void* host;
        if (!md_mtl_arena_alloc(s, size, &addr, &host, &src_buf, &src_off)) return false;
        memcpy(host, src, size);
    }

    if (s->capture) {
        md_mtl_op_t* op = (md_mtl_op_t*)md_mtl_vec_push(&s->capture->ops, dev->alloc);
        if (!op) return false;
        op->kind = MD_MTL_OP_COPY_BUFFER_TO_TEX;
        op->src_buffer = src_buf; op->src_offset = src_off;
        op->texture = t->texture;
        op->tex_origin = origin; op->tex_size = extent;
        op->bytes_per_row = bpr; op->bytes_per_image = bpi;
        md_mtl_did_op(s);
        return true;
    }
    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    [enc copyFromBuffer:src_buf sourceOffset:src_off
      sourceBytesPerRow:bpr sourceBytesPerImage:bpi sourceSize:extent
              toTexture:t->texture destinationSlice:0 destinationLevel:0 destinationOrigin:origin];
    md_mtl_did_op(s);
    return true;
}

bool md_gpu_memcpy_from_tex_async(void* dst, md_gpu_tex_t tex, const md_gpu_tex_region_t* region,
                                   size_t size, md_gpu_stream_t s) {
    if (!s) return md_mtl_fail("md_gpu_memcpy_from_tex_async requires a stream");
    md_mtl_tex_t* t = NULL;
    md_gpu_device_t dev = md_mtl_tex_device(tex, &t);
    if (!dev || !t) return md_mtl_fail("invalid texture handle");

    MTLOrigin origin; MTLSize extent;
    md_mtl_resolve_region(&t->desc, region, &origin, &extent);
    md_mtl_fmt_t fi = md_mtl_format_info(t->desc.format);
    uint64_t bpr = (uint64_t)extent.width * fi.bytes;
    uint64_t bpi = bpr * extent.height;

    md_mtl_block_t* dblk = md_mtl_registry_find(dev, (uint64_t)(uintptr_t)dst);
    id<MTLBuffer> dst_buf = nil; uint64_t dst_off = 0; void* staging_host = NULL;
    if (dblk) {
        dst_buf = dblk->buffer;
        dst_off = (uint64_t)(uintptr_t)dst - dblk->address;
    } else {
        uint64_t addr;
        if (!md_mtl_arena_alloc(s, size, &addr, &staging_host, &dst_buf, &dst_off)) return false;
    }

    if (s->capture) {
        md_mtl_op_t* op = (md_mtl_op_t*)md_mtl_vec_push(&s->capture->ops, dev->alloc);
        if (!op) return false;
        op->kind = MD_MTL_OP_COPY_TEX_TO_BUFFER;
        op->texture = t->texture;
        op->dst_buffer = dst_buf; op->dst_offset = dst_off;
        op->tex_origin = origin; op->tex_size = extent;
        op->bytes_per_row = bpr; op->bytes_per_image = bpi;
        md_mtl_did_op(s);
        return true;
    }

    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    [enc copyFromTexture:t->texture sourceSlice:0 sourceLevel:0
            sourceOrigin:origin sourceSize:extent
                toBuffer:dst_buf destinationOffset:dst_off
       destinationBytesPerRow:bpr destinationBytesPerImage:bpi];
    md_mtl_did_op(s);

    if (!dblk) {
        md_mtl_d2h_t* c = (md_mtl_d2h_t*)md_alloc(dev->alloc, sizeof(md_mtl_d2h_t));
        if (!c) return md_mtl_fail("out of memory");
        c->dev = dev; c->dst = dst; c->staging_host = staging_host; c->size = size;
        md_gpu_sync_t sync = md_gpu_stream_record(s);
        md_mutex_lock(&dev->device_mutex);
        md_mtl_hostfn_t* h = (md_mtl_hostfn_t*)md_mtl_vec_push(&dev->hostfns, dev->alloc);
        if (h) { h->sync = sync; h->fn = md_mtl_d2h_finish; h->user = c; h->internal = true; }
        md_mutex_unlock(&dev->device_mutex);
        return h != NULL;
    }
    return true;
}

/* =========================================================================
   9. Kernels and launches
   ========================================================================= */

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t dev, const md_gpu_kernel_desc_t* desc) {
    if (!dev || !desc || !desc->code || desc->code_size == 0) {
        md_mtl_fail("md_gpu_kernel_create: missing code");
        return NULL;
    }
    NSError* err = nil;
    dispatch_data_t data = dispatch_data_create(desc->code, desc->code_size,
                                                dispatch_get_main_queue(),
                                                DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> lib = [dev->device newLibraryWithData:data error:&err];
    if (!lib) {
        md_mtl_fail("newLibraryWithData failed: %s", err ? [[err localizedDescription] UTF8String] : "?");
        return NULL;
    }

    /* Metal reserves 'main', so Slang renames entry points to '<name>_0'. Try
       the requested name, then the Slang-mangled form, then the sole function. */
    const char* want = desc->entry_point ? desc->entry_point : "main";
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:want]];
    if (!fn) {
        NSString* mangled = [NSString stringWithFormat:@"%s_0", want];
        fn = [lib newFunctionWithName:mangled];
    }
    if (!fn && lib.functionNames.count == 1) {
        fn = [lib newFunctionWithName:lib.functionNames[0]];
    }
    if (!fn) { md_mtl_fail("entry point '%s' not found in metallib", want); return NULL; }

    id<MTLComputePipelineState> pso = [dev->device newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) {
        md_mtl_fail("newComputePipelineStateWithFunction failed: %s", err ? [[err localizedDescription] UTF8String] : "?");
        return NULL;
    }
    MD_MTL_RETAIN(fn);
    MD_MTL_RETAIN(pso);

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) { md_mtl_fail("out of memory"); return NULL; }
    memset(k, 0, sizeof(*k));
    k->device = dev;
    k->fn     = fn;
    k->pso    = pso;
    snprintf(k->label, sizeof(k->label), "%s", desc->label ? desc->label : "kernel");

    /* Metal cannot report the threadgroup size a kernel declares, so the
       caller must supply it. This is the one place md_gpu_kernel_desc_t's
       group_size field is mandatory rather than advisory. */
    k->group_size[0] = desc->group_size[0] ? desc->group_size[0] : 1;
    k->group_size[1] = desc->group_size[1] ? desc->group_size[1] : 1;
    k->group_size[2] = desc->group_size[2] ? desc->group_size[2] : 1;
    return k;
}

void md_gpu_kernel_destroy(md_gpu_kernel_t k) {
    if (!k) return;
    MD_MTL_RELEASE(k->pso);
    MD_MTL_RELEASE(k->fn);
    md_free(k->device->alloc, k, sizeof(*k));
}

bool md_gpu_kernel_info(md_gpu_kernel_t k, md_gpu_kernel_info_t* info) {
    if (!k || !info) return false;
    info->max_threads_per_group    = (uint32_t)k->pso.maxTotalThreadsPerThreadgroup;
    info->preferred_group_multiple = (uint32_t)k->pso.threadExecutionWidth;
    info->group_size[0] = k->group_size[0];
    info->group_size[1] = k->group_size[1];
    info->group_size[2] = k->group_size[2];
    return true;
}

static bool md_mtl_place_args(md_gpu_stream_t s, const void* args, size_t size,
                              uint64_t* out_addr, id<MTLBuffer>* out_buf, uint64_t* out_off) {
    if (size == 0) { *out_addr = 0; *out_buf = nil; *out_off = 0; return true; }

    if (s->capture) {
        md_gpu_graph_t g = s->capture;
        md_gpu_device_t dev = s->device;
        uint64_t need = md_mtl_align_up(size, MD_MTL_ARG_ALIGN);

        md_mtl_page_t* page = NULL;
        size_t page_index = 0;
        for (size_t i = 0; i < g->arg_pages.count; ++i) {
            md_mtl_page_t* p = ((md_mtl_page_t**)g->arg_pages.data)[i];
            if (p->cursor + need <= p->capacity) { page = p; page_index = i; break; }
        }
        if (!page) {
            uint64_t page_size = MD_MTL_GRAPH_ARG_PAGE;
            while (page_size < need) page_size *= 2;
            page = md_mtl_page_create(dev, page_size);
            if (!page) return md_mtl_fail("failed to allocate graph argument page");
            md_mtl_page_t** slot = (md_mtl_page_t**)md_mtl_vec_push(&g->arg_pages, dev->alloc);
            if (!slot) { md_mtl_page_destroy(dev, page); return false; }
            *slot = page;
            page_index = g->arg_pages.count - 1;
        }
        uint64_t offset = page->cursor;
        memcpy(page->host + offset, args, size);
        page->cursor += need;

        md_mtl_graph_launch_t* l = (md_mtl_graph_launch_t*)md_mtl_vec_push(&g->launches, dev->alloc);
        if (!l) return false;
        l->page_index = page_index;
        l->offset     = offset;
        l->size       = size;

        *out_addr = page->address + offset;
        *out_buf  = page->buffer;
        *out_off  = offset;
        return true;
    }

    void* host;
    return md_mtl_arena_alloc(s, size, out_addr, &host, out_buf, out_off)
        && (memcpy(host, args, size), true);
}

static void md_mtl_encode_dispatch(md_gpu_stream_t s, id<MTLComputeCommandEncoder> enc,
                                   const md_mtl_op_t* op) {
    [enc setComputePipelineState:op->kernel->pso];
    if (op->arg_buffer) {
        [enc setBuffer:op->arg_buffer offset:op->arg_offset atIndex:MD_MTL_ARG_BUFFER_INDEX];
    }
    MTLSize tg = MTLSizeMake(op->kernel->group_size[0], op->kernel->group_size[1], op->kernel->group_size[2]);
    if (op->kind == MD_MTL_OP_DISPATCH_INDIRECT) {
        [enc dispatchThreadgroupsWithIndirectBuffer:op->indirect_buffer
                               indirectBufferOffset:op->indirect_offset
                              threadsPerThreadgroup:tg];
    } else {
        MTLSize grid = MTLSizeMake(op->grid.x, op->grid.y, op->grid.z);
        [enc dispatchThreadgroups:grid threadsPerThreadgroup:tg];
    }
    (void)s;
}

static bool md_mtl_launch_common(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid,
                                 const void* args, size_t args_size,
                                 id<MTLBuffer> indirect_buf, uint64_t indirect_off, bool is_indirect) {
    if (k->device != s->device) return md_mtl_fail("kernel belongs to a different device");

    uint64_t arg_addr = 0, arg_off = 0;
    id<MTLBuffer> arg_buf = nil;
    if (!md_mtl_place_args(s, args, args_size, &arg_addr, &arg_buf, &arg_off)) return false;

    md_mtl_op_t op = {0};
    op.kind            = is_indirect ? MD_MTL_OP_DISPATCH_INDIRECT : MD_MTL_OP_DISPATCH;
    op.kernel          = k;
    op.arg_addr        = arg_addr;
    op.arg_buffer      = arg_buf;
    op.arg_offset      = arg_off;
    op.grid            = grid;
    op.indirect_buffer = indirect_buf;
    op.indirect_offset = indirect_off;

    if (s->capture) {
        md_mtl_op_t* slot = (md_mtl_op_t*)md_mtl_vec_push(&s->capture->ops, s->device->alloc);
        if (!slot) return false;
        *slot = op;
        md_mtl_did_op(s);
        return true;
    }

    id<MTLComputeCommandEncoder> enc = md_mtl_compute_encoder(s);
    if (!enc) return false;
    md_mtl_barrier_if_needed(s);
    md_mtl_encode_dispatch(s, enc, &op);
    md_mtl_did_op(s);
    return true;
}

bool md_gpu_launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid,
                    const void* args, size_t args_size) {
    if (!s || !k) return md_mtl_fail("md_gpu_launch: null stream or kernel");
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;
    return md_mtl_launch_common(s, k, grid, args, args_size, nil, 0, false);
}

bool md_gpu_launch_indirect(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_ptr_t grid_ptr,
                             const void* args, size_t args_size) {
    if (!s || !k || !grid_ptr) return md_mtl_fail("md_gpu_launch_indirect: null argument");
    md_mtl_block_t* b = md_mtl_registry_find(s->device, (uint64_t)(uintptr_t)grid_ptr);
    if (!b) return md_mtl_fail("md_gpu_launch_indirect: grid_ptr is not a device pointer");
    uint64_t off = (uint64_t)(uintptr_t)grid_ptr - b->address;
    return md_mtl_launch_common(s, k, md_gpu_grid(1, 1, 1), args, args_size, b->buffer, off, true);
}

typedef struct md_mtl_make_grid_args_t {
    uint64_t count;
    uint64_t out_grid;
    uint32_t local[3];
    uint32_t _pad;
} md_mtl_make_grid_args_t;

bool md_gpu_make_grid(md_gpu_stream_t s, md_gpu_ptr_t out_grid, const md_gpu_ptr_t count, const uint32_t local[3]) {
    if (!s || !out_grid || !count) return md_mtl_fail("md_gpu_make_grid: null argument");
    md_mtl_make_grid_args_t a = {0};
    a.count    = (uint64_t)(uintptr_t)count;
    a.out_grid = (uint64_t)(uintptr_t)out_grid;
    a.local[0] = local && local[0] ? local[0] : 1;
    a.local[1] = local && local[1] ? local[1] : 1;
    a.local[2] = local && local[2] ? local[2] : 1;
    return md_gpu_launch(s, s->device->make_grid_kernel, md_gpu_grid(1, 1, 1), &a, sizeof(a));
}

/* =========================================================================
   10. Graphs
   ========================================================================= */

bool md_gpu_capture_begin(md_gpu_stream_t s, const char* label) {
    if (!s) return md_mtl_fail("md_gpu_capture_begin: null stream");
    if (s->capture) return md_mtl_fail("stream '%s' is already capturing", s->label);
    md_gpu_device_t dev = s->device;
    md_mtl_stream_submit(s);

    md_gpu_graph_t g = (md_gpu_graph_t)md_alloc(dev->alloc, sizeof(md_gpu_graph));
    if (!g) return md_mtl_fail("out of memory");
    memset(g, 0, sizeof(*g));
    g->device = dev;
    g->kind   = s->kind;
    snprintf(g->label, sizeof(g->label), "%s", label ? label : "graph");
    md_mtl_vec_init(&g->arg_pages, sizeof(md_mtl_page_t*));
    md_mtl_vec_init(&g->launches,  sizeof(md_mtl_graph_launch_t));
    md_mtl_vec_init(&g->ops,       sizeof(md_mtl_op_t));

    s->capture       = g;
    s->needs_barrier = false;
    return true;
}

bool md_gpu_is_capturing(md_gpu_stream_t s) { return s && s->capture != NULL; }

uint32_t md_gpu_capture_next_index(md_gpu_stream_t s) {
    if (!s || !s->capture) return 0;
    return (uint32_t)s->capture->launches.count;
}

md_gpu_graph_t md_gpu_capture_end(md_gpu_stream_t s) {
    if (!s || !s->capture) { md_mtl_fail("stream is not capturing"); return NULL; }
    md_gpu_graph_t g = s->capture;
    s->capture = NULL;
    s->needs_barrier = false;
    return g;
}

uint32_t md_gpu_graph_launch_count(md_gpu_graph_t g) { return g ? (uint32_t)g->launches.count : 0; }

void* md_gpu_graph_args(md_gpu_graph_t g, uint32_t index) {
    if (!g || index >= g->launches.count) return NULL;
    md_mtl_graph_launch_t* l = &((md_mtl_graph_launch_t*)g->launches.data)[index];
    md_mtl_page_t* p = ((md_mtl_page_t**)g->arg_pages.data)[l->page_index];
    return p->host + l->offset;
}

/* Replay: re-encode the recorded operations. Argument blocks are not touched,
   so patching them through md_gpu_graph_args is picked up automatically. */
bool md_gpu_graph_launch(md_gpu_graph_t g, md_gpu_stream_t s) {
    if (!g || !s) return md_mtl_fail("md_gpu_graph_launch: null argument");
    if (s->capture) return md_mtl_fail("cannot launch a graph while capturing");
    if (g->kind != s->kind) return md_mtl_fail("graph was captured on a different stream kind");

    const md_mtl_op_t* ops = (const md_mtl_op_t*)g->ops.data;
    for (size_t i = 0; i < g->ops.count; ++i) {
        const md_mtl_op_t* op = &ops[i];
        switch (op->kind) {
        case MD_MTL_OP_DISPATCH:
        case MD_MTL_OP_DISPATCH_INDIRECT: {
            id<MTLComputeCommandEncoder> enc = md_mtl_compute_encoder(s);
            if (!enc) return false;
            md_mtl_barrier_if_needed(s);
            md_mtl_encode_dispatch(s, enc, op);
            break;
        }
        case MD_MTL_OP_COPY_BUFFER: {
            id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
            if (!enc) return false;
            [enc copyFromBuffer:op->src_buffer sourceOffset:op->src_offset
                       toBuffer:op->dst_buffer destinationOffset:op->dst_offset size:op->size];
            break;
        }
        case MD_MTL_OP_FILL_BUFFER: {
            id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
            if (!enc) return false;
            [enc fillBuffer:op->dst_buffer
                      range:NSMakeRange((NSUInteger)op->dst_offset, (NSUInteger)op->size)
                      value:op->fill_value];
            break;
        }
        case MD_MTL_OP_COPY_BUFFER_TO_TEX: {
            id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
            if (!enc) return false;
            [enc copyFromBuffer:op->src_buffer sourceOffset:op->src_offset
              sourceBytesPerRow:op->bytes_per_row sourceBytesPerImage:op->bytes_per_image
                     sourceSize:op->tex_size
                      toTexture:op->texture destinationSlice:0 destinationLevel:0
              destinationOrigin:op->tex_origin];
            break;
        }
        case MD_MTL_OP_COPY_TEX_TO_BUFFER: {
            id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
            if (!enc) return false;
            [enc copyFromTexture:op->texture sourceSlice:0 sourceLevel:0
                    sourceOrigin:op->tex_origin sourceSize:op->tex_size
                        toBuffer:op->dst_buffer destinationOffset:op->dst_offset
          destinationBytesPerRow:op->bytes_per_row destinationBytesPerImage:op->bytes_per_image];
            break;
        }
        }
        md_mtl_did_op(s);
    }
    return true;
}

void md_gpu_graph_destroy(md_gpu_graph_t g) {
    if (!g) return;
    md_gpu_device_t dev = g->device;
    for (size_t i = 0; i < g->arg_pages.count; ++i) {
        md_mtl_page_destroy(dev, ((md_mtl_page_t**)g->arg_pages.data)[i]);
    }
    md_mtl_vec_free(&g->arg_pages, dev->alloc);
    md_mtl_vec_free(&g->launches,  dev->alloc);
    md_mtl_vec_free(&g->ops,       dev->alloc);
    md_free(dev->alloc, g, sizeof(*g));
}

/* =========================================================================
   11. Host callbacks and polling
   ========================================================================= */

bool md_gpu_launch_host_fn(md_gpu_stream_t s, md_gpu_host_fn fn, void* user) {
    if (!s || !fn) return md_mtl_fail("md_gpu_launch_host_fn: null argument");
    if (s->capture) return md_mtl_fail("host functions cannot be captured into a graph");
    md_gpu_sync_t sync = md_gpu_stream_record(s);
    return md_gpu_sync_on_complete(sync, fn, user);
}

bool md_gpu_sync_on_complete(md_gpu_sync_t sync, md_gpu_host_fn fn, void* user) {
    if (!fn) return false;
    md_gpu_device_t dev = NULL;
    if (md_gpu_sync_is_valid(sync)) dev = sync.stream->device;
    else if (md_mtl_device_count > 0) dev = md_mtl_devices[0];
    if (!dev) return md_mtl_fail("md_gpu_sync_on_complete: no device");

    md_mutex_lock(&dev->device_mutex);
    md_mtl_hostfn_t* h = (md_mtl_hostfn_t*)md_mtl_vec_push(&dev->hostfns, dev->alloc);
    if (h) { h->sync = sync; h->fn = fn; h->user = user; h->internal = false; }
    md_mutex_unlock(&dev->device_mutex);
    return h != NULL;
}

/* Timeline values sampled once per poll pass. Callbacks must fire in
   registration order; re-reading signaledValue per entry breaks that, because
   two callbacks sharing a sync value can resolve differently if the timeline
   advances between the reads, letting a later one overtake an earlier one.
   See the fuller note in md_gpu_vulkan.c. */
typedef struct md_mtl_timeline_snapshot_t {
    md_gpu_stream_t stream[16];
    uint64_t         completed[16];
    uint32_t         count;
} md_mtl_timeline_snapshot_t;

static bool md_mtl_snapshot_complete(md_mtl_timeline_snapshot_t* snap, md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    for (uint32_t i = 0; i < snap->count; ++i) {
        if (snap->stream[i] == sync.stream) return snap->completed[i] >= sync.value;
    }
    uint64_t completed = md_mtl_stream_completed(sync.stream);
    if (snap->count < (uint32_t)(sizeof(snap->stream) / sizeof(snap->stream[0]))) {
        snap->stream[snap->count]    = sync.stream;
        snap->completed[snap->count] = completed;
        snap->count++;
    }
    return completed >= sync.value;
}

static uint32_t md_mtl_poll_internal(md_gpu_device_t dev, bool fire_user) {
    uint32_t fired = 0;
    md_mtl_timeline_snapshot_t snap = {0};
    for (;;) {
        md_mtl_hostfn_t ready;
        bool have = false;
        md_mutex_lock(&dev->device_mutex);
        md_mtl_hostfn_t* arr = (md_mtl_hostfn_t*)dev->hostfns.data;
        for (size_t i = 0; i < dev->hostfns.count; ++i) {
            if (!fire_user && !arr[i].internal) continue;
            if (!md_mtl_snapshot_complete(&snap, arr[i].sync)) continue;
            ready = arr[i];
            memmove(&arr[i], &arr[i + 1], (dev->hostfns.count - i - 1) * sizeof(*arr));
            dev->hostfns.count--;
            have = true;
            break;
        }
        md_mutex_unlock(&dev->device_mutex);
        if (!have) break;
        ready.fn(ready.user);
        fired++;
    }

    md_mutex_lock(&dev->device_mutex);
    md_mtl_retire_t* rr = (md_mtl_retire_t*)dev->retires.data;
    for (size_t i = 0; i < dev->retires.count;) {
        md_mtl_retire_t* r = &rr[i];
        bool done = true;
        for (uint32_t w = 0; w < r->wait_count && done; ++w) {
            if (!md_gpu_sync_is_complete(r->waits[w])) done = false;
        }
        if (!done) { ++i; continue; }
        if (r->texture) {
            md_mtl_end_residency(dev, r->texture);
            MD_MTL_RELEASE(r->texture);
        }
        memmove(&rr[i], &rr[i + 1], (dev->retires.count - i - 1) * sizeof(*rr));
        dev->retires.count--;
    }
    md_mutex_unlock(&dev->device_mutex);
    return fired;
}

uint32_t md_gpu_device_poll(md_gpu_device_t dev) {
    if (!dev) return 0;
    return md_mtl_poll_internal(dev, true);
}

/* =========================================================================
   12. Device
   ========================================================================= */

#include "md_gpu_builtin_msl.inl"

static bool md_mtl_create_builtin_kernels(md_gpu_device_t dev) {
    NSError* err = nil;
    NSString* src = [NSString stringWithUTF8String:md_gpu_make_grid_msl];
    id<MTLLibrary> lib = [dev->device newLibraryWithSource:src options:nil error:&err];
    if (!lib) return md_mtl_fail("built-in make_grid failed to compile: %s",
                                 err ? [[err localizedDescription] UTF8String] : "?");
    id<MTLFunction> fn = [lib newFunctionWithName:@"md_gpu_make_grid"];
    if (!fn) return md_mtl_fail("built-in make_grid entry point missing");
    id<MTLComputePipelineState> pso = [dev->device newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) return md_mtl_fail("built-in make_grid pipeline failed: %s",
                                 err ? [[err localizedDescription] UTF8String] : "?");
    MD_MTL_RETAIN(fn);
    MD_MTL_RETAIN(pso);

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) return md_mtl_fail("out of memory");
    memset(k, 0, sizeof(*k));
    k->device = dev;
    k->fn     = fn;
    k->pso    = pso;
    k->group_size[0] = k->group_size[1] = k->group_size[2] = 1;
    snprintf(k->label, sizeof(k->label), "md_gpu make_grid");
    dev->make_grid_kernel = k;
    return true;
}

md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc) {
    md_mtl_has_error = false;
    struct md_allocator_i* alloc = (desc && desc->alloc) ? desc->alloc : md_get_heap_allocator();

    id<MTLDevice> mtl = MTLCreateSystemDefaultDevice();
    if (!mtl) { md_mtl_fail("MTLCreateSystemDefaultDevice returned nil"); return NULL; }

    md_gpu_device_t dev = (md_gpu_device_t)md_alloc(alloc, sizeof(md_gpu_device));
    if (!dev) { md_mtl_fail("out of memory"); return NULL; }
    memset(dev, 0, sizeof(*dev));
    dev->alloc  = alloc;
    dev->device = mtl;
    MD_MTL_RETAIN(mtl);
    dev->is_discrete = ![mtl hasUnifiedMemory];

    md_mtl_vec_init(&dev->registry,  sizeof(md_mtl_block_t*));
    md_mtl_vec_init(&dev->pools,     sizeof(md_gpu_pool_t));
    md_mtl_vec_init(&dev->streams,   sizeof(md_gpu_stream_t));
    md_mtl_vec_init(&dev->hostfns,   sizeof(md_mtl_hostfn_t));
    md_mtl_vec_init(&dev->retires,   sizeof(md_mtl_retire_t));
    md_mtl_vec_init(&dev->live_bufs, sizeof(id));
    md_mtl_vec_init(&dev->textures,  sizeof(md_mtl_tex_t*));
    md_mtl_vec_init(&dev->samplers,  sizeof(md_mtl_sampler_t*));

    md_mutex_init(&dev->queue_mutex);
    md_mutex_init(&dev->device_mutex);

    dev->compute_queue  = [mtl newCommandQueue];
    dev->transfer_queue = [mtl newCommandQueue];
    MD_MTL_RETAIN(dev->compute_queue);
    MD_MTL_RETAIN(dev->transfer_queue);

    /* Device-wide residency, where available. */
    Class rsd = NSClassFromString(@"MTLResidencySetDescriptor");
    if (rsd && [mtl respondsToSelector:@selector(newResidencySetWithDescriptor:error:)]) {
        id descriptor = [[rsd alloc] init];
        NSError* err = nil;
        id set = [mtl performSelector:@selector(newResidencySetWithDescriptor:error:)
                           withObject:descriptor withObject:(id)&err];
        if (set) {
            MD_MTL_RETAIN(set);
            dev->residency_set     = set;
            dev->has_residency_set = true;
            [dev->compute_queue  performSelector:@selector(addResidencySet:) withObject:set];
            [dev->transfer_queue performSelector:@selector(addResidencySet:) withObject:set];
        }
    }
    if (!dev->has_residency_set) {
        MD_LOG_DEBUG("md_gpu: MTLResidencySet unavailable, falling back to per-encoder useResources");
    }

    if (md_mtl_device_count >= (uint32_t)(sizeof(md_mtl_devices) / sizeof(md_mtl_devices[0]))) {
        md_mtl_fail("too many simultaneous md_gpu devices");
        goto fail;
    }
    md_mtl_devices[md_mtl_device_count++] = dev;

    if (!md_mtl_stream_create_internal(dev, MD_GPU_STREAM_COMPUTE,  "default compute",  true)) goto fail;
    if (!md_mtl_stream_create_internal(dev, MD_GPU_STREAM_TRANSFER, "default transfer", true)) goto fail;
    if (!md_mtl_create_builtin_kernels(dev)) goto fail;

    MD_LOG_DEBUG("md_gpu: device '%s'", [[mtl name] UTF8String]);
    return dev;

fail:
    md_gpu_device_destroy(dev);
    return NULL;
}

bool md_gpu_device_info(md_gpu_device_t dev, md_gpu_device_info_t* info) {
    if (!dev || !info) return false;
    memset(info, 0, sizeof(*info));
    info->is_discrete              = dev->is_discrete;
    info->max_threads_per_group    = (uint32_t)dev->device.maxThreadsPerThreadgroup.width;
    info->preferred_group_multiple = 32;
    info->timestamp_period_ns_num  = 1;
    info->timestamp_period_ns_den  = 1;
    snprintf(info->name, sizeof(info->name), "%s", [[dev->device name] UTF8String]);
    return true;
}

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;

    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = ((md_gpu_stream_t*)dev->streams.data)[i];
        s->capture = NULL;
        md_mtl_stream_submit(s);
    }
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = ((md_gpu_stream_t*)dev->streams.data)[i];
        while (s->submitted_value > 0 && md_mtl_stream_completed(s) < s->submitted_value) md_thread_sleep(0);
    }
    md_mtl_poll_internal(dev, true);

    if (dev->make_grid_kernel) md_gpu_kernel_destroy(dev->make_grid_kernel);

    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_mtl_stream_destroy_internal(((md_gpu_stream_t*)dev->streams.data)[i]);
    }
    dev->streams.count = 0;
    md_mtl_vec_free(&dev->streams, alloc);

    while (dev->pools.count > 0) md_gpu_pool_destroy(((md_gpu_pool_t*)dev->pools.data)[0]);
    md_mtl_vec_free(&dev->pools, alloc);

    for (size_t i = 0; i < dev->textures.count; ++i) {
        md_mtl_tex_t* t = ((md_mtl_tex_t**)dev->textures.data)[i];
        if (t->texture) MD_MTL_RELEASE(t->texture);
        md_free(alloc, t, sizeof(*t));
    }
    md_mtl_vec_free(&dev->textures, alloc);
    for (size_t i = 0; i < dev->samplers.count; ++i) {
        md_mtl_sampler_t* sm = ((md_mtl_sampler_t**)dev->samplers.data)[i];
        if (sm->sampler) MD_MTL_RELEASE(sm->sampler);
        md_free(alloc, sm, sizeof(*sm));
    }
    md_mtl_vec_free(&dev->samplers, alloc);

    md_mtl_vec_free(&dev->hostfns,   alloc);
    md_mtl_vec_free(&dev->retires,   alloc);
    md_mtl_vec_free(&dev->registry,  alloc);
    md_mtl_vec_free(&dev->live_bufs, alloc);

    if (dev->residency_set) MD_MTL_RELEASE(dev->residency_set);
    if (dev->compute_queue)  MD_MTL_RELEASE(dev->compute_queue);
    if (dev->transfer_queue) MD_MTL_RELEASE(dev->transfer_queue);

    md_mutex_destroy(&dev->queue_mutex);
    md_mutex_destroy(&dev->device_mutex);

    if (dev->device) MD_MTL_RELEASE(dev->device);

    for (uint32_t i = 0; i < md_mtl_device_count; ++i) {
        if (md_mtl_devices[i] == dev) {
            md_mtl_devices[i] = md_mtl_devices[md_mtl_device_count - 1];
            md_mtl_device_count--;
            break;
        }
    }
    md_free(alloc, dev, sizeof(*dev));
}

md_gpu_mem_flags_t md_gpu_pool_flags(md_gpu_pool_t p) { return p ? p->flags : 0; }

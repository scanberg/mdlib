/*
md_gpu_metal.m -- Metal backend for md_gpu.h

Mirrors md_gpu_vulkan.c section for section. Where the two backends differ,
the difference is called out in a comment.

Structural differences from Vulkan:

  1. Program order inside a command buffer costs nothing to express. A serial
     compute encoder orders its dispatches, and an MTLFence orders work across
     encoder boundaries (see md_mtl_close_encoder), so this backend emits no
     per-operation barriers.

     Across command buffers it is not free, and for a reason worth stating up
     front: see NOTE ON RESIDENCY below. Because nothing is ever declared to a
     dispatch, Metal has no hazard tracking to drive its implicit cross-submission
     ordering, so each new command buffer explicitly waits on the stream's own
     previous signal value. See md_mtl_stream_ensure_cmd.

  2. EXPLICIT ordering and md_gpu_barrier are accepted and do nothing on this
     (Metal 3) path: serial encoders already order everything, so explicit code
     is correct here, merely not faster. A Metal 4 path with real stage barriers
     is a separate piece of work.

  3. There are no image layouts and no descriptor heap. A shader handle is the
     resource's gpuResourceID, sitting directly in the argument struct.

  4. Kernels need no deferred destruction: a command buffer retains the
     pipeline states set on its encoders. Buffers and textures do, because they
     are reached by address / resource id and nothing retains them.

NOTE ON RESIDENCY: with no per-dispatch resource declarations there is nothing
to derive residency from, so every live allocation goes into a device-wide
MTLResidencySet attached to every queue (macOS 15 / iOS 18 and later). On
older systems the backend falls back to useResources: with the full live set
on each compute encoder, repeated before a dispatch whenever something became
resident since -- correct, but O(live allocations) each time.

NOTE ON BLIT ALIGNMENT: on macOS, buffer-to-buffer blits and fillBuffer want
offsets and sizes in multiples of 4. md_gpu promises byte granularity, so
unaligned copies and the ragged ends of fills go through a small built-in
compute kernel (md_gpu_byte_op). It runs on transfer streams too.

NOTE ON ROOT ALIGNMENT: the argument root cell is bound with setBuffer:offset:
in the constant address space, which on macOS needs 256-byte offsets; root
cells are allocated at that alignment, argument blocks at 64.

The second, less obvious consequence: a residency set makes resources resident
and does explicitly nothing for hazards. Combined with reaching every buffer by
raw gpuAddress, it means Metal never learns which resources a dispatch touches
and so cannot insert the implicit memory barriers it normally would. Anywhere
ordering is not already structural -- between encoders and between command
buffers -- this backend states the dependency itself.

NOTE ON OWNERSHIP: this file must be correct whether or not it is compiled with
ARC, so Objective-C objects held in C structs go through three macros:
MD_MTL_OWN for a +1 result (new/alloc/copy/Create), MD_MTL_RETAIN for an
autoreleased one, MD_MTL_RELEASE to let go. Temporaries created +1 are dropped
with MD_MTL_DROP_NEW. Entry points that create autoreleased objects wrap their
work in @autoreleasepool, since they may be called from threads without one.

Nothing here blocks the calling thread except md_gpu_stream_sync,
md_gpu_sync_wait, md_gpu_stream_destroy (on its own stream) and device
destruction.
*/

#include "md_gpu.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <objc/message.h>

/* This file targets macOS 13+ at runtime but has, historically, only compiled
   against the macOS 15 SDK. Everything genuinely macOS-15-only is reached
   dynamically (NSClassFromString / respondsToSelector: / objc_msgSend), so
   the only thing an older SDK actually lacks is this spelling of the enum --
   renamed from MTLPipelineOptionArgumentInfo in Xcode 16, same value, 1 << 0.
   Guard on the SDK version (MAX_ALLOWED), never on the deployment target. */
#if !defined(__MAC_OS_X_VERSION_MAX_ALLOWED) || __MAC_OS_X_VERSION_MAX_ALLOWED < 150000
#define MTLPipelineOptionBindingInfo MTLPipelineOptionArgumentInfo
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =========================================================================
   1. Configuration, error handling, utilities
   ========================================================================= */

#define MD_MTL_ARENA_PAGE_SIZE (256u * 1024u)
#define MD_MTL_ARG_ALIGN          64u
/* setBuffer:offset: for the constant address space wants 256-byte offsets on
   macOS (Mac2 family GPUs); the root cell is bound that way, so it gets them. */
#define MD_MTL_ROOT_ALIGN        256u
/* Bytes each thread of the built-in byte kernel handles. */
#define MD_MTL_BYTE_OP_SPAN       16u
#define MD_MTL_MAX_SAMPLERS      256u
#define MD_MTL_ERROR_BUF         512u

/* The root buffer index is NOT fixed, and it is NOT always 0.

   Slang assigns Metal buffer indices per *file*, in declaration order of the
   push-constant globals -- not per entry point. A file with one kernel always
   gets 0; unittest/shaders/gpu_test.slang declares a root per entry point, so
   its kernels land on 0, 1, 2, ...

   So the index is a property of the compiled kernel and is read out of the
   pipeline's binding reflection at create time. This constant is only the
   fallback for when reflection is unavailable. */
#define MD_MTL_ARG_BUFFER_INDEX    0

#if __has_feature(objc_arc)
#  define MD_MTL_OWN(obj)      ((void)CFBridgingRetain(obj))
#  define MD_MTL_RETAIN(obj)   ((void)CFBridgingRetain(obj))
#  define MD_MTL_DROP_NEW(obj) ((void)(obj))
#else
#  define MD_MTL_OWN(obj)      ((void)(obj))
#  define MD_MTL_RETAIN(obj)   ((void)CFRetain((__bridge CFTypeRef)(obj)))
#  define MD_MTL_DROP_NEW(obj) do { if (obj) CFRelease((__bridge CFTypeRef)(obj)); } while (0)
#endif
#define MD_MTL_RELEASE(obj) do { if (obj) CFRelease((__bridge CFTypeRef)(obj)); } while (0)

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

static void md_mtl_vec_remove(md_mtl_vec_t* v, size_t i) {
    char* base = (char*)v->data;
    memmove(base + i * v->stride, base + (i + 1) * v->stride, (v->count - i - 1) * v->stride);
    v->count--;
}

/* For vectors of pointers: remove the first element equal to `p`. */
static void md_mtl_vec_remove_ptr(md_mtl_vec_t* v, const void* p) {
    void** arr = (void**)v->data;
    for (size_t i = 0; i < v->count; ++i) {
        if (arr[i] == p) { md_mtl_vec_remove(v, i); return; }
    }
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
    uint64_t          address;
    void*             host;
    uint64_t          capacity;
    uint64_t          size;
    md_gpu_mem_kind_t kind;
    md_gpu_pool_t     pool;
    bool              in_use;
    md_gpu_stream_t   free_stream;
    uint64_t          free_value;
} md_mtl_block_t;

typedef struct md_gpu_pool {
    md_gpu_device_t   device;
    md_gpu_mem_kind_t kind;
    uint64_t          cache_limit;   /* 0 = unlimited */
    md_mtl_vec_t      blocks;        /* md_mtl_block_t*  */
    md_mtl_vec_t      textures;      /* md_gpu_texture_t */
    uint64_t          in_use_bytes;
    uint64_t          reserved_bytes;
    uint64_t          peak_in_use_bytes;
    uint64_t          alloc_count;
    uint64_t          reuse_count;
    char              label[64];
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

typedef struct md_mtl_wait_t {
    md_gpu_stream_t stream;
    uint64_t        value;
} md_mtl_wait_t;

typedef struct md_gpu_stream {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    __unsafe_unretained id<MTLCommandQueue>          queue;
    __unsafe_unretained id<MTLSharedEvent>           timeline;
    __unsafe_unretained id<MTLCommandBuffer>         cmd;
    __unsafe_unretained id<MTLComputeCommandEncoder> compute_enc;
    __unsafe_unretained id<MTLBlitCommandEncoder>    blit_enc;
    /* Orders work across an encoder boundary. See md_mtl_close_encoder. */
    __unsafe_unretained id<MTLFence>                 fence;
    bool                                             fence_valid;
    uint64_t                                         res_gen;       /* residency generation declared to compute_enc */

    uint64_t          next_value;
    uint64_t          submitted_value;
    bool              has_work;
    md_gpu_ordering_t ordering;       /* recorded; see note 2 at the top */

    md_mtl_vec_t      waits;          /* md_mtl_wait_t, encoded at the next command buffer */

    md_mtl_arena_t    arena;

    bool              upload_open;
    bool              upload_direct;
    md_gpu_addr_t     upload_dst;
    uint64_t          upload_src_addr;
    size_t            upload_size;

    bool              is_default;
    char              label[64];
} md_gpu_stream;

typedef struct md_mtl_fmt_t {
    MTLPixelFormat fmt;
    uint32_t       bytes;       /* texel size in buffer copies (depth plane only for D32S8) */
    bool           depth;
    bool           stencil;
    bool           srgb;
    const char*    name;
} md_mtl_fmt_t;

typedef struct md_gpu_texture {
    md_gpu_device_t       device;
    md_gpu_pool_t         pool;
    __unsafe_unretained id<MTLTexture> texture;
    /* One view per mip when the texture has several and STORAGE usage: a
       read_write texture in a shader addresses a single level. NULL when the
       texture itself serves (one mip). */
    __unsafe_unretained id<MTLTexture>* mip_views;
    uint64_t*             storage_handles;   /* per mip, or NULL */
    uint64_t              sampled_handle;    /* 0 if not SAMPLED */
    uint64_t              bytes;
    md_mtl_fmt_t          fi;
    md_gpu_texture_desc_t desc;              /* normalised */
    char                  label[64];
} md_gpu_texture;

typedef struct md_mtl_sampler_entry_t {
    md_gpu_sampler_desc_t desc;
    __unsafe_unretained id<MTLSamplerState> sampler;
    uint64_t              handle;
} md_mtl_sampler_entry_t;

typedef struct md_gpu_kernel {
    md_gpu_device_t device;
    __unsafe_unretained id<MTLComputePipelineState> pso;
    uint32_t group_size[3];
    uint32_t args_size;
    uint32_t arg_buffer_index;   /* from binding reflection; see the note above */
    char     label[64];
} md_gpu_kernel;

typedef struct md_mtl_hostfn_t {
    md_gpu_sync_t  sync;
    md_gpu_host_fn fn;
    void*          user;
} md_mtl_hostfn_t;

typedef enum md_mtl_retire_kind_t {
    MD_MTL_RETIRE_BLOCK,
    MD_MTL_RETIRE_TEXTURE,
} md_mtl_retire_kind_t;

typedef struct md_mtl_retire_t {
    md_mtl_retire_kind_t kind;
    void*                object;
    md_mtl_wait_t*       waits;
    uint32_t             wait_count;
    uint32_t             wait_capacity;
} md_mtl_retire_t;

typedef struct md_gpu_device {
    struct md_allocator_i* alloc;
    __unsafe_unretained id<MTLDevice> device;
    __unsafe_unretained id            residency_set;   /* id<MTLResidencySet> */
    bool     has_residency_set;

    md_mutex_t queue_mutex;
    md_mutex_t device_mutex;

    md_mtl_vec_t registry;   /* md_mtl_block_t*, sorted by address */
    md_mtl_vec_t pools;      /* md_gpu_pool_t   */
    md_mtl_vec_t kernels;    /* md_gpu_kernel_t */
    md_mtl_vec_t streams;    /* md_gpu_stream_t */
    md_mtl_vec_t hostfns;    /* md_mtl_hostfn_t */
    md_mtl_vec_t retires;    /* md_mtl_retire_t */
    md_mtl_vec_t live_res;   /* id<MTLResource>, for the no-residency-set path */
    uint64_t     res_gen;    /* bumped whenever live_res gains a resource */

    md_gpu_stream_t default_compute;
    md_gpu_stream_t default_transfer;

    md_mtl_sampler_entry_t samplers[MD_MTL_MAX_SAMPLERS];
    uint32_t               sampler_count;

    md_gpu_kernel_t make_grid_kernel;
    md_gpu_kernel_t byte_op_kernel;   /* unaligned copies and fills */
    bool            is_discrete;
} md_gpu_device;

static bool     md_mtl_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t align, uint64_t* out_addr, void** out_host, id<MTLBuffer>* out_buf, uint64_t* out_off);
static bool     md_mtl_stream_submit(md_gpu_stream_t s);
static bool     md_mtl_byte_op(md_gpu_stream_t s, uint64_t dst, uint64_t src, uint64_t size, uint8_t value);
static uint64_t md_mtl_stream_completed(md_gpu_stream_t s);
static void     md_mtl_block_free(md_gpu_device_t dev, md_mtl_block_t* b);
static void     md_mtl_texture_free(md_gpu_device_t dev, md_gpu_texture_t t);

/* =========================================================================
   3. Allocation registry
   ========================================================================= */

/* Caller holds device_mutex. */
static md_mtl_block_t* md_mtl_registry_find_locked(md_gpu_device_t dev, uint64_t address) {
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

static bool md_mtl_registry_insert_locked(md_gpu_device_t dev, md_mtl_block_t* blk) {
    if (!md_mtl_vec_reserve(&dev->registry, dev->alloc, dev->registry.count + 1)) return false;
    md_mtl_block_t** arr = (md_mtl_block_t**)dev->registry.data;
    size_t i = dev->registry.count;
    while (i > 0 && arr[i - 1]->address > blk->address) { arr[i] = arr[i - 1]; i--; }
    arr[i] = blk;
    dev->registry.count++;
    return true;
}

static void md_mtl_registry_remove_locked(md_gpu_device_t dev, md_mtl_block_t* blk) {
    md_mtl_vec_remove_ptr(&dev->registry, blk);
}

/* Resolve [addr, addr + size) to a live allocation, under the device lock --
   malloc on another thread mutates the registry. */
static md_mtl_block_t* md_mtl_resolve(md_gpu_device_t dev, md_gpu_addr_t addr, uint64_t size,
                                      uint64_t* out_offset, const char* what) {
    md_mutex_lock(&dev->device_mutex);
    md_mtl_block_t* b = md_mtl_registry_find_locked(dev, addr);
    md_mutex_unlock(&dev->device_mutex);
    if (!b || !b->in_use) {
        md_mtl_fail("%s: 0x%llx is not a live md_gpu allocation", what, (unsigned long long)addr);
        return NULL;
    }
    uint64_t off = addr - b->address;
    if (size > b->size || off > b->size - size) {
        md_mtl_fail("%s: range [0x%llx, +%llu) overruns its %llu-byte allocation",
                    what, (unsigned long long)addr, (unsigned long long)size, (unsigned long long)b->size);
        return NULL;
    }
    if (out_offset) *out_offset = off;
    return b;
}

/* =========================================================================
   4. Residency, buffers and transient arenas
   ========================================================================= */

/* MTLResidencySet is reached through typed objc_msgSend casts rather than
   performSelector: -- the latter is wrong for void returns and, under ARC,
   for a 'new' method's +1 result. */
typedef void (*md_mtl_msg_v_t)(id, SEL);
typedef void (*md_mtl_msg_vo_t)(id, SEL, id);

static void md_mtl_msg0(id obj, SEL sel)        { ((md_mtl_msg_v_t)objc_msgSend)(obj, sel); }
static void md_mtl_msg1(id obj, SEL sel, id arg) { ((md_mtl_msg_vo_t)objc_msgSend)(obj, sel, arg); }

/* Caller holds device_mutex (it guards live_res). */
static void md_mtl_make_resident_locked(md_gpu_device_t dev, id<MTLResource> res) {
    if (dev->has_residency_set) {
        md_mtl_msg1(dev->residency_set, @selector(addAllocation:), res);
        md_mtl_msg0(dev->residency_set, @selector(commit));
    } else {
        __unsafe_unretained id<MTLResource>* slot = (__unsafe_unretained id<MTLResource>*)md_mtl_vec_push(&dev->live_res, dev->alloc);
        if (slot) *slot = res;
        dev->res_gen++;
    }
}

/* Caller holds device_mutex. */
static void md_mtl_end_residency_locked(md_gpu_device_t dev, id<MTLResource> res) {
    if (dev->has_residency_set) {
        md_mtl_msg1(dev->residency_set, @selector(removeAllocation:), res);
        md_mtl_msg0(dev->residency_set, @selector(commit));
    } else {
        md_mtl_vec_remove_ptr(&dev->live_res, (__bridge const void*)res);
    }
}

static bool md_mtl_create_raw_buffer(md_gpu_device_t dev, uint64_t size, md_gpu_mem_kind_t kind,
                                     id<MTLBuffer>* out_buf, uint64_t* out_addr, void** out_host) {
    const MTLResourceOptions opts = (kind != MD_GPU_MEM_DEVICE)
        ? MTLResourceStorageModeShared
        : MTLResourceStorageModePrivate;
    id<MTLBuffer> buf = [dev->device newBufferWithLength:(NSUInteger)size options:opts];
    if (!buf) return md_mtl_fail("newBufferWithLength failed for %llu bytes", (unsigned long long)size);
    MD_MTL_OWN(buf);

    *out_buf  = buf;
    *out_addr = (uint64_t)[buf gpuAddress];
    *out_host = (kind != MD_GPU_MEM_DEVICE) ? [buf contents] : NULL;
    return true;
}

/* Residency is the caller's: made resident by whoever creates the buffer,
   ended here. Caller holds device_mutex. */
static void md_mtl_destroy_raw_buffer_locked(md_gpu_device_t dev, id<MTLBuffer> buf) {
    if (!buf) return;
    md_mtl_end_residency_locked(dev, buf);
    MD_MTL_RELEASE(buf);
}

static md_mtl_page_t* md_mtl_page_create(md_gpu_device_t dev, uint64_t size) {
    md_mtl_page_t* p = (md_mtl_page_t*)md_alloc(dev->alloc, sizeof(md_mtl_page_t));
    if (!p) return NULL;
    memset(p, 0, sizeof(*p));
    void* host = NULL;
    id<MTLBuffer> buf = nil;
    /* Argument structs live in these pages and the shader reaches them by raw
       device address, so a page is made resident exactly like a pool block. */
    if (!md_mtl_create_raw_buffer(dev, size, MD_GPU_MEM_HOST_WRITE, &buf, &p->address, &host)) {
        md_free(dev->alloc, p, sizeof(md_mtl_page_t));
        return NULL;
    }
    p->buffer   = buf;
    p->host     = (uint8_t*)host;
    p->capacity = size;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_make_resident_locked(dev, buf);
    md_mutex_unlock(&dev->device_mutex);
    return p;
}

static void md_mtl_page_destroy(md_gpu_device_t dev, md_mtl_page_t* p) {
    md_mutex_lock(&dev->device_mutex);
    md_mtl_destroy_raw_buffer_locked(dev, p->buffer);
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, p, sizeof(md_mtl_page_t));
}

/* A page taken for new data is unstamped: it now holds work that has not been
   submitted, and md_mtl_arena_retire restamps it at the next submit. Leaving
   the old (possibly completed) stamp would let a later allocation recycle the
   page -- cursor back to 0 -- over data the open command buffer still needs. */
static bool md_mtl_page_fits(const md_mtl_page_t* p, uint64_t need, uint64_t align) {
    return md_mtl_align_up(p->cursor, align) + need <= p->capacity;
}

static void md_mtl_page_take(md_mtl_page_t* p, uint64_t need, uint64_t align, uint64_t* out_addr, void** out_host,
                             id<MTLBuffer>* out_buf, uint64_t* out_off) {
    const uint64_t at = md_mtl_align_up(p->cursor, align);
    *out_addr = p->address + at;
    if (out_host) *out_host = p->host + at;
    if (out_buf)  *out_buf  = p->buffer;
    if (out_off)  *out_off  = at;
    p->cursor       = at + need;
    p->retire_value = 0;
}

static bool md_mtl_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t align, uint64_t* out_addr, void** out_host,
                               id<MTLBuffer>* out_buf, uint64_t* out_off) {
    md_gpu_device_t dev = s->device;
    md_mtl_arena_t*  a   = &s->arena;
    if (align < MD_MTL_ARG_ALIGN) align = MD_MTL_ARG_ALIGN;
    const uint64_t need = md_mtl_align_up(size, MD_MTL_ARG_ALIGN);

    md_mtl_page_t* chosen = NULL;
    if (a->pages.count > 0) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, a->current);
        if (md_mtl_page_fits(p, need, align)) chosen = p;
    }
    if (!chosen) {
        uint64_t done = md_mtl_stream_completed(s);
        for (size_t i = 0; i < a->pages.count && !chosen; ++i) {
            md_mtl_page_t* p = MD_MTL_VEC_AT(a->pages, md_mtl_page_t*, i);
            if (p->retire_value != 0 && p->retire_value <= done) { p->cursor = 0; p->retire_value = 0; }
            if (p->retire_value == 0 && md_mtl_page_fits(p, need, align)) { a->current = i; chosen = p; }
        }
    }
    if (!chosen) {
        uint64_t page_size = MD_MTL_ARENA_PAGE_SIZE;
        while (page_size < need) page_size *= 2;
        chosen = md_mtl_page_create(dev, page_size);
        if (!chosen) return md_mtl_fail("failed to allocate a transient page");
        md_mtl_page_t** slot = (md_mtl_page_t**)md_mtl_vec_push(&a->pages, dev->alloc);
        if (!slot) { md_mtl_page_destroy(dev, chosen); return md_mtl_fail("out of memory"); }
        *slot = chosen;
        a->current = a->pages.count - 1;
    }
    md_mtl_page_take(chosen, need, align, out_addr, out_host, out_buf, out_off);
    return true;
}

/* Stamp every page holding data with the value that releases it. Must
   overwrite an existing stamp -- see the note in md_gpu_vulkan.c. */
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
   5. Deferred destruction
   ========================================================================= */

static uint64_t md_mtl_stream_position(md_gpu_stream_t s) {
    return s->has_work ? s->next_value : s->submitted_value;
}

/* Queue `object` for release once every stream has passed its current
   position. Caller holds device_mutex. Never blocks, except when the wait
   list itself cannot be allocated. */
static void md_mtl_retire_locked(md_gpu_device_t dev, md_mtl_retire_kind_t kind, void* object) {
    md_mtl_retire_t r;
    memset(&r, 0, sizeof(r));
    r.kind   = kind;
    r.object = object;

    const uint32_t n = (uint32_t)dev->streams.count;
    if (n > 0) r.waits = (md_mtl_wait_t*)md_alloc(dev->alloc, n * sizeof(md_mtl_wait_t));
    if (r.waits) {
        r.wait_capacity = n;
        for (uint32_t i = 0; i < n; ++i) {
            md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
            uint64_t v = md_mtl_stream_position(s);
            if (v > 0 && md_mtl_stream_completed(s) < v) {
                r.waits[r.wait_count].stream = s;
                r.waits[r.wait_count].value  = v;
                r.wait_count++;
            }
        }
    }

    md_mtl_retire_t* slot = NULL;
    if (n == 0 || r.waits) slot = (md_mtl_retire_t*)md_mtl_vec_push(&dev->retires, dev->alloc);
    if (!slot) {
        md_mtl_fail("out of memory recording a deferred destruction; waiting for the device");
        if (r.waits) md_free(dev->alloc, r.waits, n * sizeof(md_mtl_wait_t));
        for (uint32_t i = 0; i < n; ++i) {
            md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
            while (s->submitted_value > 0 && md_mtl_stream_completed(s) < s->submitted_value) md_thread_sleep(0);
        }
        if (kind == MD_MTL_RETIRE_BLOCK) md_mtl_block_free(dev, (md_mtl_block_t*)object);
        else                             md_mtl_texture_free(dev, (md_gpu_texture_t)object);
        return;
    }
    *slot = r;
}

/* Caller holds device_mutex. */
static void md_mtl_process_retires_locked(md_gpu_device_t dev, bool force) {
    for (size_t i = 0; i < dev->retires.count;) {
        md_mtl_retire_t* r = &MD_MTL_VEC_AT(dev->retires, md_mtl_retire_t, i);
        bool done = true;
        for (uint32_t w = 0; w < r->wait_count && done && !force; ++w) {
            if (md_mtl_stream_completed(r->waits[w].stream) < r->waits[w].value) done = false;
        }
        if (!done) { ++i; continue; }
        md_mtl_retire_t e = *r;
        md_mtl_vec_remove(&dev->retires, i);
        if (e.kind == MD_MTL_RETIRE_BLOCK) md_mtl_block_free(dev, (md_mtl_block_t*)e.object);
        else                               md_mtl_texture_free(dev, (md_gpu_texture_t)e.object);
        if (e.waits) md_free(dev->alloc, e.waits, e.wait_capacity * sizeof(md_mtl_wait_t));
    }
}

/* A stream is going away after being synchronised: drop every reference to
   it. Caller holds device_mutex. */
static void md_mtl_forget_stream_locked(md_gpu_device_t dev, md_gpu_stream_t s) {
    for (size_t i = 0; i < dev->pools.count; ++i) {
        md_gpu_pool_t pool = MD_MTL_VEC_AT(dev->pools, md_gpu_pool_t, i);
        for (size_t j = 0; j < pool->blocks.count; ++j) {
            md_mtl_block_t* b = MD_MTL_VEC_AT(pool->blocks, md_mtl_block_t*, j);
            if (b->free_stream == s) { b->free_stream = NULL; b->free_value = 0; }
        }
    }
    for (size_t i = 0; i < dev->hostfns.count; ++i) {
        md_mtl_hostfn_t* h = &MD_MTL_VEC_AT(dev->hostfns, md_mtl_hostfn_t, i);
        if (h->sync.stream == s) h->sync = md_gpu_sync_none();
    }
    for (size_t i = 0; i < dev->retires.count; ++i) {
        md_mtl_retire_t* r = &MD_MTL_VEC_AT(dev->retires, md_mtl_retire_t, i);
        for (uint32_t w = 0; w < r->wait_count;) {
            if (r->waits[w].stream == s) r->waits[w] = r->waits[--r->wait_count];
            else ++w;
        }
    }
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t o = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
        for (size_t w = 0; w < o->waits.count;) {
            if (MD_MTL_VEC_AT(o->waits, md_mtl_wait_t, w).stream == s) md_mtl_vec_remove(&o->waits, w);
            else ++w;
        }
    }
}

/* =========================================================================
   6. Encoders and program order
   ========================================================================= */

static bool md_mtl_stream_ensure_cmd(md_gpu_stream_t s) {
    if (s->cmd) return true;
    @autoreleasepool {
        id<MTLCommandBuffer> cb = [s->queue commandBuffer];
        if (!cb) return md_mtl_fail("commandBuffer failed on stream '%s'", s->label);
        MD_MTL_RETAIN(cb);
        cb.label = [NSString stringWithUTF8String:s->label];
        s->cmd         = cb;
        s->has_work    = false;
        /* An MTLFence is scoped to a command buffer: waiting on one that has
           not been updated in this buffer is undefined, so re-arm the flag. */
        s->fence_valid = false;

        /* Chain onto this stream's own timeline. Metal starts command buffers
           on a queue in commit order, but the memory ordering between them is
           driven by hazard tracking, which this backend gives Metal nothing to
           do (see NOTE ON RESIDENCY). One event wait per command buffer, only
           at a submission boundary. */
        if (s->submitted_value > 0) {
            [cb encodeWaitForEvent:s->timeline value:s->submitted_value];
        }
        /* Pending cross-stream waits go at the head of the buffer. */
        for (size_t i = 0; i < s->waits.count; ++i) {
            md_mtl_wait_t w = MD_MTL_VEC_AT(s->waits, md_mtl_wait_t, i);
            [cb encodeWaitForEvent:w.stream->timeline value:w.value];
        }
        s->waits.count = 0;
    }
    return true;
}

/* Close whichever encoder is open, signalling the stream's fence on the way out
   so the next encoder can wait on it.

   Encoders inside one command buffer are NOT guaranteed to run one after
   another: Metal decides whether they may overlap from hazard tracking, which
   this backend does not give it. A dispatch writing a texture through a
   resource id followed by a blit reading it is, as far as Metal can tell, two
   independent pieces of work. MTLFence orders untracked access across encoders
   and is cheaper than splitting command buffers. `signal_fence` is false only
   when the command buffer itself is ending. */
static void md_mtl_close_encoder(md_gpu_stream_t s, bool signal_fence) {
    if (!s->compute_enc && !s->blit_enc) return;
    if (signal_fence && s->fence) {
        if (s->compute_enc) [s->compute_enc updateFence:s->fence];
        else                [s->blit_enc    updateFence:s->fence];
        s->fence_valid = true;
    }
    if (s->compute_enc) { [s->compute_enc endEncoding]; MD_MTL_RELEASE(s->compute_enc); s->compute_enc = nil; }
    if (s->blit_enc)    { [s->blit_enc    endEncoding]; MD_MTL_RELEASE(s->blit_enc);    s->blit_enc    = nil; }
}

/* An operation may not be recorded while an upload is open: the staged copy
   the upload ends with must land where the upload began. */
static bool md_mtl_check_no_upload(md_gpu_stream_t s, const char* what) {
    if (s->upload_open) return md_mtl_fail("%s: stream '%s' has an open upload; call md_gpu_upload_end first", what, s->label);
    return true;
}

static id<MTLComputeCommandEncoder> md_mtl_compute_encoder(md_gpu_stream_t s) {
    if (!md_mtl_stream_ensure_cmd(s)) return nil;
    if (s->compute_enc) return s->compute_enc;
    md_mtl_close_encoder(s, true);

    /* MTLDispatchTypeSerial -- the default, and deliberately so: it is exactly
       md_gpu's IMPLICIT ordering, and Metal gets it right by construction.
       memoryBarrierWithScope: is only defined on concurrent encoders; calling
       it on a serial one corrupts execution order, which is how this backend
       once landed 6 of 40 dependent dispatches. */
    @autoreleasepool {
        id<MTLComputeCommandEncoder> enc = [s->cmd computeCommandEncoder];
        if (!enc) { md_mtl_fail("computeCommandEncoder failed"); return nil; }
        MD_MTL_RETAIN(enc);
        if (s->fence_valid) [enc waitForFence:s->fence];
        s->compute_enc = enc;
        s->res_gen     = 0;
    }
    return s->compute_enc;
}

/* Fallback residency (no MTLResidencySet): the open compute encoder must have
   declared every live resource before a dispatch that may touch it. Declared
   in full when the encoder opens and again whenever something became resident
   since -- a buffer malloc'd, a texture created, an arena page added between
   two launches into the same encoder. Call right before each dispatch. */
static void md_mtl_declare_residency(md_gpu_stream_t s, id<MTLComputeCommandEncoder> enc) {
    md_gpu_device_t dev = s->device;
    if (dev->has_residency_set) return;
    md_mutex_lock(&dev->device_mutex);
    if (s->res_gen != dev->res_gen + 1) {
        if (dev->live_res.count > 0) {
            [enc useResources:(__unsafe_unretained id<MTLResource>*)dev->live_res.data
                        count:dev->live_res.count
                        usage:MTLResourceUsageRead | MTLResourceUsageWrite];
        }
        s->res_gen = dev->res_gen + 1;   /* +1 so that 0 always means "not yet declared" */
    }
    md_mutex_unlock(&dev->device_mutex);
}

static id<MTLBlitCommandEncoder> md_mtl_blit_encoder(md_gpu_stream_t s) {
    if (!md_mtl_stream_ensure_cmd(s)) return nil;
    if (s->blit_enc) return s->blit_enc;
    md_mtl_close_encoder(s, true);
    @autoreleasepool {
        id<MTLBlitCommandEncoder> enc = [s->cmd blitCommandEncoder];
        if (!enc) { md_mtl_fail("blitCommandEncoder failed"); return nil; }
        MD_MTL_RETAIN(enc);
        if (s->fence_valid) [enc waitForFence:s->fence];
        s->blit_enc = enc;
    }
    return s->blit_enc;
}

static void md_mtl_did_op(md_gpu_stream_t s) {
    s->has_work = true;
}

void md_gpu_stream_set_ordering(md_gpu_stream_t s, md_gpu_ordering_t ordering) {
    if (s) s->ordering = ordering;
}

md_gpu_ordering_t md_gpu_stream_ordering(md_gpu_stream_t s) {
    return s ? s->ordering : MD_GPU_ORDER_IMPLICIT;
}

/* Metal 3: serial encoders and the fence already order everything, so a stage
   barrier has nothing to add. Code written for EXPLICIT stays correct here. */
void md_gpu_barrier(md_gpu_stream_t s, md_gpu_stage_flags_t producers, md_gpu_stage_flags_t consumers) {
    (void)producers; (void)consumers;
    if (s) md_mtl_check_no_upload(s, "md_gpu_barrier");
}

/* =========================================================================
   7. Streams
   ========================================================================= */

static uint64_t md_mtl_stream_completed(md_gpu_stream_t s) {
    return (uint64_t)s->timeline.signaledValue;
}

/* Block until `ev` reaches `value`. MTLSharedEvent has a real blocking wait
   (macOS 12+); spinning is only the fallback for a runtime without it. */
static void md_mtl_event_wait(id<MTLSharedEvent> ev, uint64_t value) {
    if ((uint64_t)ev.signaledValue >= value) return;
    if ([ev respondsToSelector:@selector(waitUntilSignaledValue:timeoutMS:)]) {
        while (![ev waitUntilSignaledValue:value timeoutMS:1000]) {}
    } else {
        while ((uint64_t)ev.signaledValue < value) md_thread_sleep(0);
    }
}

static md_gpu_stream_t md_mtl_stream_create_internal(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label, bool is_default) {
    md_gpu_stream_t s = (md_gpu_stream_t)md_alloc(dev->alloc, sizeof(md_gpu_stream));
    if (!s) { md_mtl_fail("out of memory"); return NULL; }
    memset(s, 0, sizeof(*s));
    s->device     = dev;
    s->kind       = kind;
    s->ordering   = MD_GPU_ORDER_IMPLICIT;
    s->next_value = 1;
    s->is_default = is_default;
    snprintf(s->label, sizeof(s->label), "%s", label ? label : "stream");
    md_mtl_vec_init(&s->arena.pages, sizeof(md_mtl_page_t*));
    md_mtl_vec_init(&s->waits, sizeof(md_mtl_wait_t));

    @autoreleasepool {
        /* Its own hardware queue, so two streams can genuinely run at once. */
        id<MTLCommandQueue> q = [dev->device newCommandQueue];
        if (!q) { md_free(dev->alloc, s, sizeof(*s)); md_mtl_fail("newCommandQueue failed"); return NULL; }
        MD_MTL_OWN(q);
        q.label = [NSString stringWithUTF8String:s->label];
        s->queue = q;
        if (dev->has_residency_set) {
            md_mtl_msg1(q, @selector(addResidencySet:), dev->residency_set);
        }

        id<MTLSharedEvent> ev = [dev->device newSharedEvent];
        if (!ev) {
            MD_MTL_RELEASE(q);
            md_free(dev->alloc, s, sizeof(*s));
            md_mtl_fail("newSharedEvent failed");
            return NULL;
        }
        MD_MTL_OWN(ev);
        s->timeline = ev;

        id<MTLFence> fence = [dev->device newFence];
        if (!fence) {
            MD_MTL_RELEASE(ev);
            MD_MTL_RELEASE(q);
            md_free(dev->alloc, s, sizeof(*s));
            md_mtl_fail("newFence failed");
            return NULL;
        }
        MD_MTL_OWN(fence);
        s->fence = fence;
    }

    md_mutex_lock(&dev->device_mutex);
    md_gpu_stream_t* slot = (md_gpu_stream_t*)md_mtl_vec_push(&dev->streams, dev->alloc);
    if (slot) *slot = s;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) {
        MD_MTL_RELEASE(s->fence);
        MD_MTL_RELEASE(s->timeline);
        MD_MTL_RELEASE(s->queue);
        md_free(dev->alloc, s, sizeof(*s));
        md_mtl_fail("out of memory");
        return NULL;
    }
    return s;
}

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label) {
    if (!dev) { md_mtl_fail("md_gpu_stream_create: null device"); return NULL; }
    return md_mtl_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    return kind == MD_GPU_STREAM_TRANSFER ? dev->default_transfer : dev->default_compute;
}

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t s) { return s ? s->device : NULL; }

static bool md_mtl_stream_submit(md_gpu_stream_t s) {
    if (!s->cmd || !s->has_work) return true;
    /* The command buffer is ending; cross-submission ordering is the timeline
       wait in ensure_cmd, so no fence is needed here. */
    md_mtl_close_encoder(s, false);

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
    return true;
}

md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t s) {
    md_gpu_sync_t out = md_gpu_sync_none();
    if (!s) return out;
    if (!md_mtl_check_no_upload(s, "md_gpu_stream_record")) return out;
    md_mtl_stream_submit(s);
    if (s->submitted_value == 0) return out;
    out.stream = s;
    out.value  = s->submitted_value;
    return out;
}

void md_gpu_stream_wait(md_gpu_stream_t s, md_gpu_sync_t sync) {
    if (!s || !md_gpu_sync_is_valid(sync)) return;
    if (sync.stream == s) return;
    if (sync.stream->device != s->device) { md_mtl_fail("md_gpu_stream_wait: sync from another device"); return; }
    if (md_gpu_sync_is_complete(sync)) return;

    /* Work already issued must not be retroactively delayed: close it first. */
    if (s->has_work) md_mtl_stream_submit(s);
    /* A command buffer that is open but empty (an op failed after opening it)
       already carries the waits ensure_cmd encoded; dropping it would lose
       them. Encode this wait onto it too -- legal between encoders. */
    if (s->cmd && !s->has_work) {
        md_mtl_close_encoder(s, false);
        [s->cmd encodeWaitForEvent:sync.stream->timeline value:sync.value];
        return;
    }

    for (size_t i = 0; i < s->waits.count; ++i) {
        md_mtl_wait_t* w = &MD_MTL_VEC_AT(s->waits, md_mtl_wait_t, i);
        if (w->stream == sync.stream) {
            if (sync.value > w->value) w->value = sync.value;
            return;
        }
    }
    md_mtl_wait_t* w = (md_mtl_wait_t*)md_mtl_vec_push(&s->waits, s->device->alloc);
    if (!w) { md_mtl_fail("out of memory queueing a stream wait"); return; }
    w->stream = sync.stream;
    w->value  = sync.value;
}

void md_gpu_stream_flush(md_gpu_stream_t s) {
    if (s) md_mtl_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    md_mtl_stream_submit(s);
    if (s->submitted_value > 0) md_mtl_event_wait(s->timeline, s->submitted_value);
}

bool md_gpu_sync_is_complete(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    return md_mtl_stream_completed(sync.stream) >= sync.value;
}

void md_gpu_sync_wait(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return;
    md_mtl_event_wait(sync.stream->timeline, sync.value);
}

/* Free a stream's own objects. The stream must be idle. */
static void md_mtl_stream_free(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    md_mtl_close_encoder(s, false);
    if (s->cmd) { MD_MTL_RELEASE(s->cmd); s->cmd = nil; }
    md_mtl_arena_free(dev, &s->arena);
    md_mtl_vec_free(&s->waits, dev->alloc);
    MD_MTL_RELEASE(s->queue);
    MD_MTL_RELEASE(s->fence);
    MD_MTL_RELEASE(s->timeline);
    md_free(dev->alloc, s, sizeof(*s));
}

void md_gpu_stream_destroy(md_gpu_stream_t s) {
    if (!s || s->is_default) return;
    md_gpu_device_t dev = s->device;

    /* Only this stream's own work. Anything waiting on it is thereby
       satisfied, which is what makes forgetting it below safe. */
    s->upload_open = false;
    md_gpu_stream_sync(s);

    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->streams, s);
    md_mtl_forget_stream_locked(dev, s);
    md_mutex_unlock(&dev->device_mutex);

    md_mtl_stream_free(s);
}

/* =========================================================================
   8. Memory
   ========================================================================= */

md_gpu_pool_t md_gpu_pool_create(md_gpu_device_t dev, const md_gpu_pool_desc_t* desc) {
    if (!dev || !desc) { md_mtl_fail("md_gpu_pool_create: null argument"); return NULL; }
    if ((unsigned)desc->kind > (unsigned)MD_GPU_MEM_HOST_READ) { md_mtl_fail("md_gpu_pool_create: invalid memory kind %d", (int)desc->kind); return NULL; }
    md_gpu_pool_t p = (md_gpu_pool_t)md_alloc(dev->alloc, sizeof(md_gpu_pool));
    if (!p) { md_mtl_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device      = dev;
    p->kind        = desc->kind;
    p->cache_limit = desc->cache_limit;
    snprintf(p->label, sizeof(p->label), "%s", desc->label ? desc->label : "pool");
    md_mtl_vec_init(&p->blocks,   sizeof(md_mtl_block_t*));
    md_mtl_vec_init(&p->textures, sizeof(md_gpu_texture_t));

    md_mutex_lock(&dev->device_mutex);
    md_gpu_pool_t* slot = (md_gpu_pool_t*)md_mtl_vec_push(&dev->pools, dev->alloc);
    if (slot) *slot = p;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { md_free(dev->alloc, p, sizeof(*p)); md_mtl_fail("out of memory"); return NULL; }
    return p;
}

md_gpu_mem_kind_t md_gpu_pool_kind(md_gpu_pool_t p) { return p ? p->kind : MD_GPU_MEM_DEVICE; }

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
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (MD_MTL_VEC_AT(p->blocks, md_mtl_block_t*, i)->in_use) out->blocks_in_use++;
        else                                                     out->blocks_cached++;
    }
    md_mutex_unlock(&p->device->device_mutex);
}

/* Caller holds device_mutex. */
static void md_mtl_block_release(md_mtl_block_t* b, md_gpu_stream_t stream) {
    b->in_use      = false;
    b->free_stream = stream;
    b->free_value  = md_mtl_stream_position(stream);
    if (b->free_value == 0 || md_mtl_stream_completed(stream) >= b->free_value) {
        b->free_stream = NULL;
        b->free_value  = 0;
    }
    b->pool->in_use_bytes -= b->capacity;
}

static bool md_mtl_block_idle(md_mtl_block_t* b) {
    return !b->in_use && (b->free_stream == NULL || md_mtl_stream_completed(b->free_stream) >= b->free_value);
}

/* Out of the registry and its pool already. Caller holds device_mutex. */
static void md_mtl_block_free(md_gpu_device_t dev, md_mtl_block_t* b) {
    md_mtl_destroy_raw_buffer_locked(dev, b->buffer);
    md_free(dev->alloc, b, sizeof(*b));
}

/* Caller holds device_mutex. */
static void md_mtl_pool_trim_locked(md_gpu_pool_t p, uint64_t keep_bytes) {
    md_gpu_device_t dev = p->device;
    for (size_t i = 0; i < p->blocks.count && p->reserved_bytes - p->in_use_bytes > keep_bytes;) {
        md_mtl_block_t* b = MD_MTL_VEC_AT(p->blocks, md_mtl_block_t*, i);
        if (md_mtl_block_idle(b)) {
            md_mtl_registry_remove_locked(dev, b);
            p->reserved_bytes -= b->capacity;
            md_mtl_vec_remove(&p->blocks, i);
            md_mtl_block_free(dev, b);
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

/* Caller holds device_mutex and has removed `t` from its pool's list. */
static void md_mtl_texture_retire_locked(md_gpu_device_t dev, md_gpu_texture_t t) {
    if (t->pool) {
        t->pool->in_use_bytes   -= t->bytes;
        t->pool->reserved_bytes -= t->bytes;
        t->pool = NULL;
    }
    md_mtl_retire_locked(dev, MD_MTL_RETIRE_TEXTURE, t);
}

void md_gpu_pool_reset(md_gpu_stream_t stream, md_gpu_pool_t p) {
    if (!p || !stream) { md_mtl_fail("md_gpu_pool_reset: null argument"); return; }
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_mtl_block_t* b = MD_MTL_VEC_AT(p->blocks, md_mtl_block_t*, i);
        if (b->in_use) md_mtl_block_release(b, stream);
    }
    while (p->textures.count > 0) {
        md_gpu_texture_t t = MD_MTL_VEC_AT(p->textures, md_gpu_texture_t, p->textures.count - 1);
        p->textures.count--;
        md_mtl_texture_retire_locked(dev, t);
    }
    md_mutex_unlock(&dev->device_mutex);
}

void md_gpu_pool_destroy(md_gpu_pool_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    /* Nothing is freed here: every block and texture goes onto the retire list
       and is released by md_gpu_device_poll once every stream has passed the
       work issued before this call. */
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_mtl_block_t* b = MD_MTL_VEC_AT(p->blocks, md_mtl_block_t*, i);
        md_mtl_registry_remove_locked(dev, b);
        b->pool = NULL;
        md_mtl_retire_locked(dev, MD_MTL_RETIRE_BLOCK, b);
    }
    md_mtl_vec_free(&p->blocks, dev->alloc);
    while (p->textures.count > 0) {
        md_gpu_texture_t t = MD_MTL_VEC_AT(p->textures, md_gpu_texture_t, p->textures.count - 1);
        p->textures.count--;
        md_mtl_texture_retire_locked(dev, t);
    }
    md_mtl_vec_free(&p->textures, dev->alloc);
    md_mtl_vec_remove_ptr(&dev->pools, p);
    md_mtl_process_retires_locked(dev, false);
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, p, sizeof(*p));
}

md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t stream, md_gpu_pool_t p, size_t size) {
    md_gpu_mem_t out = {0, NULL};
    if (!stream || !p) { md_mtl_fail("md_gpu_malloc: null stream or pool"); return out; }
    if (size == 0) return out;
    md_gpu_device_t dev = p->device;
    if (stream->device != dev) { md_mtl_fail("md_gpu_malloc: stream and pool belong to different devices"); return out; }
    uint64_t capacity = md_mtl_next_pow2(size);

    md_mutex_lock(&dev->device_mutex);

    md_mtl_block_t* best = NULL;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_mtl_block_t* b = MD_MTL_VEC_AT(p->blocks, md_mtl_block_t*, i);
        if (b->in_use || b->capacity < size) continue;
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
        out.gpu = best->address;
        out.cpu = best->host;
        return out;
    }

    md_mtl_block_t* b = (md_mtl_block_t*)md_alloc(dev->alloc, sizeof(md_mtl_block_t));
    if (!b) { md_mutex_unlock(&dev->device_mutex); md_mtl_fail("out of memory"); return out; }
    memset(b, 0, sizeof(*b));
    id<MTLBuffer> buf = nil;
    if (!md_mtl_create_raw_buffer(dev, capacity, p->kind, &buf, &b->address, &b->host)) {
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        return out;
    }
    md_mtl_make_resident_locked(dev, buf);
    b->buffer   = buf;
    b->capacity = capacity;
    b->size     = size;
    b->kind     = p->kind;
    b->pool     = p;
    b->in_use   = true;

    md_mtl_block_t** slot = (md_mtl_block_t**)md_mtl_vec_push(&p->blocks, dev->alloc);
    if (!slot || !md_mtl_registry_insert_locked(dev, b)) {
        if (slot) p->blocks.count--;
        md_mtl_block_free(dev, b);
        md_mutex_unlock(&dev->device_mutex);
        md_mtl_fail("out of memory");
        return out;
    }
    *slot = b;
    p->reserved_bytes += capacity;
    p->in_use_bytes   += capacity;
    if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
    p->alloc_count++;

    md_mutex_unlock(&dev->device_mutex);
    out.gpu = b->address;
    out.cpu = b->host;
    return out;
}

void md_gpu_free(md_gpu_stream_t stream, md_gpu_addr_t addr) {
    if (!addr) return;
    if (!stream) { md_mtl_fail("md_gpu_free: a stream is required"); return; }
    md_gpu_device_t dev = stream->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_block_t* b = md_mtl_registry_find_locked(dev, addr);
    if (!b || !b->in_use || b->address != addr) {
        md_mutex_unlock(&dev->device_mutex);
        md_mtl_fail("md_gpu_free: 0x%llx is not the start of a live allocation", (unsigned long long)addr);
        return;
    }
    md_mtl_block_release(b, stream);
    md_gpu_pool_t p = b->pool;
    if (p->cache_limit != 0) md_mtl_pool_trim_locked(p, p->cache_limit);
    md_mutex_unlock(&dev->device_mutex);
}

/* ---- Copies ------------------------------------------------------------------ */

/* On macOS a blit buffer copy needs offsets and size in multiples of 4 (the
   documented rule; Apple GPUs tolerate less, Intel/AMD do not). Anything else
   goes through the built-in byte kernel, by GPU address. */
static bool md_mtl_record_buffer_copy(md_gpu_stream_t s, id<MTLBuffer> src, uint64_t src_off, uint64_t src_addr,
                                      id<MTLBuffer> dst, uint64_t dst_off, uint64_t dst_addr, uint64_t size) {
    if (!md_mtl_check_no_upload(s, "copy")) return false;
    if ((src_off | dst_off | size) & 3u) return md_mtl_byte_op(s, dst_addr, src_addr, size, 0);
    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    [enc copyFromBuffer:src sourceOffset:src_off toBuffer:dst destinationOffset:dst_off size:size];
    md_mtl_did_op(s);
    return true;
}

bool md_gpu_copy(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_addr_t src, size_t size) {
    if (!s) return md_mtl_fail("md_gpu_copy: null stream");
    if (size == 0) return true;
    uint64_t doff, soff;
    md_mtl_block_t* d  = md_mtl_resolve(s->device, dst, size, &doff, "md_gpu_copy (dst)");
    if (!d) return false;
    md_mtl_block_t* sb = md_mtl_resolve(s->device, src, size, &soff, "md_gpu_copy (src)");
    if (!sb) return false;
    return md_mtl_record_buffer_copy(s, sb->buffer, soff, src, d->buffer, doff, dst, size);
}

/* Idle means a host write now cannot race anything this stream is ordered
   after: nothing recorded, nothing in flight, and no cross-stream wait still
   pending -- whether queued, or already encoded into an open command buffer
   (the reason an open buffer disqualifies even when it holds no work). */
static bool md_mtl_stream_idle(md_gpu_stream_t s) {
    if (s->has_work || s->cmd) return false;
    if (md_mtl_stream_completed(s) < s->submitted_value) return false;
    for (size_t i = 0; i < s->waits.count; ++i) {
        md_mtl_wait_t w = MD_MTL_VEC_AT(s->waits, md_mtl_wait_t, i);
        if (md_mtl_stream_completed(w.stream) < w.value) return false;
    }
    return true;
}

bool md_gpu_upload(md_gpu_stream_t s, md_gpu_addr_t dst, const void* src, size_t size) {
    if (!s) return md_mtl_fail("md_gpu_upload: null stream");
    if (size == 0) return true;
    if (!src) return md_mtl_fail("md_gpu_upload: null source");
    if (!md_mtl_check_no_upload(s, "md_gpu_upload")) return false;
    uint64_t doff;
    md_mtl_block_t* d = md_mtl_resolve(s->device, dst, size, &doff, "md_gpu_upload");
    if (!d) return false;
    if (d->host && md_mtl_stream_idle(s)) {
        memcpy((uint8_t*)d->host + doff, src, size);
        return true;
    }
    uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
    if (!md_mtl_arena_alloc(s, size, MD_MTL_ARG_ALIGN, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    return md_mtl_record_buffer_copy(s, buf, off, addr, d->buffer, doff, dst, size);
}

bool md_gpu_memset(md_gpu_stream_t s, md_gpu_addr_t dst, uint8_t value, size_t size) {
    if (!s) return md_mtl_fail("md_gpu_memset: null stream");
    if (size == 0) return true;
    if (!md_mtl_check_no_upload(s, "md_gpu_memset")) return false;
    uint64_t off;
    md_mtl_block_t* b = md_mtl_resolve(s->device, dst, size, &off, "md_gpu_memset");
    if (!b) return false;
    /* fillBuffer's range must be 4-byte aligned on macOS, as for copies: fill
       the aligned middle, and the (at most 3-byte) head and tail by kernel. */
    const uint64_t begin = off, end = off + size;
    const uint64_t abeg  = md_mtl_align_up(begin, 4);
    const uint64_t aend  = end & ~3ull;
    if (aend > abeg) {
        id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
        if (!enc) return false;
        [enc fillBuffer:b->buffer range:NSMakeRange((NSUInteger)abeg, (NSUInteger)(aend - abeg)) value:value];
        md_mtl_did_op(s);
        if (abeg > begin && !md_mtl_byte_op(s, dst, 0, abeg - begin, value)) return false;
        if (end  > aend  && !md_mtl_byte_op(s, dst + (aend - begin), 0, end - aend, value)) return false;
        return true;
    }
    return md_mtl_byte_op(s, dst, 0, size, value);
}

void* md_gpu_upload_begin(md_gpu_stream_t s, md_gpu_addr_t dst, size_t size) {
    if (!s || !dst || size == 0) { md_mtl_fail("md_gpu_upload_begin: null argument"); return NULL; }
    if (s->upload_open) { md_mtl_fail("an upload is already open on stream '%s'", s->label); return NULL; }
    uint64_t doff;
    md_mtl_block_t* b = md_mtl_resolve(s->device, dst, size, &doff, "md_gpu_upload_begin");
    if (!b) return NULL;
    if (b->host && md_mtl_stream_idle(s)) {
        s->upload_open = true; s->upload_direct = true;
        s->upload_dst = dst; s->upload_size = size;
        return (uint8_t*)b->host + doff;
    }
    uint64_t addr; void* host;
    if (!md_mtl_arena_alloc(s, size, MD_MTL_ARG_ALIGN, &addr, &host, NULL, NULL)) return NULL;
    s->upload_open     = true;
    s->upload_direct   = false;
    s->upload_dst      = dst;
    s->upload_src_addr = addr;
    s->upload_size     = size;
    return host;
}

bool md_gpu_upload_end(md_gpu_stream_t s) {
    if (!s || !s->upload_open) return md_mtl_fail("md_gpu_upload_end: no upload is open");
    s->upload_open = false;
    if (s->upload_direct) return true;

    uint64_t doff;
    md_mtl_block_t* d = md_mtl_resolve(s->device, s->upload_dst, s->upload_size, &doff, "md_gpu_upload_end");
    if (!d) return false;
    for (size_t i = 0; i < s->arena.pages.count; ++i) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(s->arena.pages, md_mtl_page_t*, i);
        if (s->upload_src_addr >= p->address && s->upload_src_addr < p->address + p->capacity) {
            return md_mtl_record_buffer_copy(s, p->buffer, s->upload_src_addr - p->address, s->upload_src_addr,
                                             d->buffer, doff, s->upload_dst, s->upload_size);
        }
    }
    return md_mtl_fail("md_gpu_upload_end: staging page not found");
}

/* =========================================================================
   9. Textures and samplers
   ========================================================================= */

static md_mtl_fmt_t md_mtl_format_info(md_gpu_format_t f) {
    md_mtl_fmt_t i = {MTLPixelFormatInvalid, 0, false, false, false, "invalid"};
#define MD_MTL_FMT(e, mf, b, dep, sten, srgb) case e: i = (md_mtl_fmt_t){mf, b, dep, sten, srgb, #e}; break
    switch (f) {
    MD_MTL_FMT(MD_GPU_FORMAT_R8_UNORM,          MTLPixelFormatR8Unorm,               1, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RG8_UNORM,         MTLPixelFormatRG8Unorm,              2, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGBA8_UNORM,       MTLPixelFormatRGBA8Unorm,            4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGBA8_SRGB,        MTLPixelFormatRGBA8Unorm_sRGB,       4, false, false, true);
    MD_MTL_FMT(MD_GPU_FORMAT_BGRA8_UNORM,       MTLPixelFormatBGRA8Unorm,            4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_BGRA8_SRGB,        MTLPixelFormatBGRA8Unorm_sRGB,       4, false, false, true);
    MD_MTL_FMT(MD_GPU_FORMAT_R16_FLOAT,         MTLPixelFormatR16Float,              2, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RG16_FLOAT,        MTLPixelFormatRG16Float,             4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGBA16_FLOAT,      MTLPixelFormatRGBA16Float,           8, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_R32_FLOAT,         MTLPixelFormatR32Float,              4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RG32_FLOAT,        MTLPixelFormatRG32Float,             8, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGBA32_FLOAT,      MTLPixelFormatRGBA32Float,          16, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_R32_UINT,          MTLPixelFormatR32Uint,               4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RG32_UINT,         MTLPixelFormatRG32Uint,              8, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGBA32_UINT,       MTLPixelFormatRGBA32Uint,           16, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RG11B10_FLOAT,     MTLPixelFormatRG11B10Float,          4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_RGB10A2_UNORM,     MTLPixelFormatRGB10A2Unorm,          4, false, false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_D32_FLOAT,         MTLPixelFormatDepth32Float,          4, true,  false, false);
    MD_MTL_FMT(MD_GPU_FORMAT_D32_FLOAT_S8_UINT, MTLPixelFormatDepth32Float_Stencil8, 4, true,  true,  false);
    default: break;
    }
#undef MD_MTL_FMT
    return i;
}

uint32_t md_gpu_format_texel_size(md_gpu_format_t format) {
    return md_mtl_format_info(format).bytes;
}

static MTLTextureType md_mtl_texture_type(md_gpu_tex_type_t t) {
    switch (t) {
    case MD_GPU_TEX_2D_ARRAY: return MTLTextureType2DArray;
    case MD_GPU_TEX_3D:       return MTLTextureType3D;
    default:                  return MTLTextureType2D;
    }
}

/* Release a texture and everything it holds. Caller holds device_mutex. */
static void md_mtl_texture_free(md_gpu_device_t dev, md_gpu_texture_t t) {
    const uint32_t mips = t->desc.mip_levels;
    if (t->mip_views) {
        for (uint32_t m = 0; m < mips; ++m) MD_MTL_RELEASE(t->mip_views[m]);
        md_free(dev->alloc, (void*)t->mip_views, mips * sizeof(id));
    }
    if (t->storage_handles) md_free(dev->alloc, t->storage_handles, mips * sizeof(uint64_t));
    if (t->texture) {
        md_mtl_end_residency_locked(dev, t->texture);
        MD_MTL_RELEASE(t->texture);
    }
    md_free(dev->alloc, t, sizeof(*t));
}

md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t s, md_gpu_pool_t pool, const md_gpu_texture_desc_t* desc) {
    if (!s || !pool || !desc) { md_mtl_fail("md_gpu_texture_create: null argument"); return NULL; }
    md_gpu_device_t dev = s->device;
    if (pool->device != dev)             { md_mtl_fail("md_gpu_texture_create: stream and pool belong to different devices"); return NULL; }
    if (pool->kind != MD_GPU_MEM_DEVICE) { md_mtl_fail("md_gpu_texture_create: pool '%s' is not an MD_GPU_MEM_DEVICE pool", pool->label); return NULL; }
    if (!md_mtl_check_no_upload(s, "md_gpu_texture_create")) return NULL;

    md_gpu_texture_desc_t d = *desc;
    const char* label = d.label ? d.label : "texture";
    md_mtl_fmt_t fi = md_mtl_format_info(d.format);
    if (fi.fmt == MTLPixelFormatInvalid) { md_mtl_fail("texture '%s': invalid format %d", label, (int)d.format); return NULL; }
    if (d.type != MD_GPU_TEX_2D && d.type != MD_GPU_TEX_2D_ARRAY && d.type != MD_GPU_TEX_3D) {
        md_mtl_fail("texture '%s': invalid type %d (a zero-initialised desc has no type)", label, (int)d.type);
        return NULL;
    }
    if (!(d.usage & (MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED | MD_GPU_TEX_RENDER_TARGET))) {
        md_mtl_fail("texture '%s': usage is empty", label);
        return NULL;
    }
    if (d.width == 0 || d.height == 0) { md_mtl_fail("texture '%s': zero width or height", label); return NULL; }
    if (d.type == MD_GPU_TEX_2D && d.depth_or_layers > 1) {
        md_mtl_fail("texture '%s': a 2D texture has depth_or_layers %u; use MD_GPU_TEX_3D or MD_GPU_TEX_2D_ARRAY",
                    label, d.depth_or_layers);
        return NULL;
    }
    /* Metal has no format-capability query; reject what it never supports,
       by name, and leave the rest to newTextureWithDescriptor:. */
    if ((d.usage & MD_GPU_TEX_STORAGE) && (fi.depth || fi.srgb)) {
        md_mtl_fail("texture '%s': the device does not support %s with usage STORAGE", label, fi.name);
        return NULL;
    }
    if (fi.depth && d.type == MD_GPU_TEX_3D) {
        md_mtl_fail("texture '%s': the device does not support 3D %s", label, fi.name);
        return NULL;
    }
    if (d.depth_or_layers == 0) d.depth_or_layers = 1;
    if (d.mip_levels == 0) d.mip_levels = 1;

    md_gpu_texture_t t = (md_gpu_texture_t)md_alloc(dev->alloc, sizeof(md_gpu_texture));
    if (!t) { md_mtl_fail("out of memory"); return NULL; }
    memset(t, 0, sizeof(*t));
    t->device = dev;
    t->fi     = fi;
    t->desc   = d;
    snprintf(t->label, sizeof(t->label), "%s", label);
    t->desc.label = t->label;

    bool ok = true;
    @autoreleasepool {
        MTLTextureDescriptor* td = [[MTLTextureDescriptor alloc] init];
        td.textureType      = md_mtl_texture_type(d.type);
        td.pixelFormat      = fi.fmt;
        td.width            = d.width;
        td.height           = d.height;
        td.depth            = d.type == MD_GPU_TEX_3D       ? d.depth_or_layers : 1;
        td.arrayLength      = d.type == MD_GPU_TEX_2D_ARRAY ? d.depth_or_layers : 1;
        td.mipmapLevelCount = d.mip_levels;
        td.storageMode      = MTLStorageModePrivate;
        MTLTextureUsage usage = 0;
        if (d.usage & MD_GPU_TEX_STORAGE)       usage |= MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;
        if (d.usage & MD_GPU_TEX_SAMPLED)       usage |= MTLTextureUsageShaderRead;
        if (d.usage & MD_GPU_TEX_RENDER_TARGET) usage |= MTLTextureUsageRenderTarget;
        td.usage = usage;

        id<MTLTexture> tex = [dev->device newTextureWithDescriptor:td];
        MD_MTL_DROP_NEW(td);
        if (!tex) {
            ok = md_mtl_fail("texture '%s': newTextureWithDescriptor failed for %s", label, fi.name);
        } else {
            MD_MTL_OWN(tex);
            tex.label = [NSString stringWithUTF8String:t->label];
            t->texture = tex;
            t->bytes   = (uint64_t)[tex allocatedSize];
            if (d.usage & MD_GPU_TEX_SAMPLED) t->sampled_handle = (uint64_t)tex.gpuResourceID._impl;
        }

        if (ok && (d.usage & MD_GPU_TEX_STORAGE)) {
            t->storage_handles = (uint64_t*)md_alloc(dev->alloc, d.mip_levels * sizeof(uint64_t));
            if (!t->storage_handles) ok = md_mtl_fail("out of memory");
            else if (d.mip_levels == 1) {
                t->storage_handles[0] = (uint64_t)tex.gpuResourceID._impl;
            } else {
                /* A read_write texture addresses one level, so each mip gets a
                   view of its own. Same format, so no PixelFormatView usage. */
                t->mip_views = (__unsafe_unretained id<MTLTexture>*)md_alloc(dev->alloc, d.mip_levels * sizeof(id));
                if (!t->mip_views) ok = md_mtl_fail("out of memory");
                else memset((void*)t->mip_views, 0, d.mip_levels * sizeof(id));
                const NSUInteger slices = d.type == MD_GPU_TEX_2D_ARRAY ? d.depth_or_layers : 1;
                for (uint32_t m = 0; ok && m < d.mip_levels; ++m) {
                    id<MTLTexture> v = [tex newTextureViewWithPixelFormat:fi.fmt
                                                              textureType:md_mtl_texture_type(d.type)
                                                                   levels:NSMakeRange(m, 1)
                                                                   slices:NSMakeRange(0, slices)];
                    if (!v) { ok = md_mtl_fail("texture '%s': view of mip %u failed", label, m); break; }
                    MD_MTL_OWN(v);
                    t->mip_views[m]       = v;
                    t->storage_handles[m] = (uint64_t)v.gpuResourceID._impl;
                }
            }
        }
    }

    md_mutex_lock(&dev->device_mutex);
    if (ok) {
        md_gpu_texture_t* slot = (md_gpu_texture_t*)md_mtl_vec_push(&pool->textures, dev->alloc);
        if (!slot) ok = md_mtl_fail("out of memory");
        else {
            *slot   = t;
            t->pool = pool;
            md_mtl_make_resident_locked(dev, t->texture);
            pool->in_use_bytes   += t->bytes;
            pool->reserved_bytes += t->bytes;
            if (pool->in_use_bytes > pool->peak_in_use_bytes) pool->peak_in_use_bytes = pool->in_use_bytes;
        }
    }
    if (!ok) {
        /* Never made resident; drop the texture without ending residency. */
        id<MTLTexture> tex = t->texture;
        t->texture = nil;
        md_mtl_texture_free(dev, t);
        MD_MTL_RELEASE(tex);
        t = NULL;
    }
    md_mutex_unlock(&dev->device_mutex);
    (void)s;   /* Metal has no layout transition to order; creation is complete now. */
    return t;
}

void md_gpu_texture_destroy(md_gpu_texture_t t) {
    if (!t) return;
    md_gpu_device_t dev = t->device;
    md_mutex_lock(&dev->device_mutex);
    if (t->pool) md_mtl_vec_remove_ptr(&t->pool->textures, t);
    md_mtl_texture_retire_locked(dev, t);
    md_mutex_unlock(&dev->device_mutex);
}

const md_gpu_texture_desc_t* md_gpu_texture_desc(md_gpu_texture_t t) {
    return t ? &t->desc : NULL;
}

md_gpu_storage_tex_t md_gpu_texture_storage(md_gpu_texture_t t, uint32_t mip) {
    md_gpu_storage_tex_t h = {0};
    if (t && t->storage_handles && mip < t->desc.mip_levels) h.handle = t->storage_handles[mip];
    return h;
}

md_gpu_sampled_tex_t md_gpu_texture_sampled(md_gpu_texture_t t) {
    md_gpu_sampled_tex_t h = {0};
    if (t) h.handle = t->sampled_handle;
    return h;
}

/* ---- Samplers ------------------------------------------------------------------ */

static bool md_mtl_sampler_desc_eq(const md_gpu_sampler_desc_t* a, const md_gpu_sampler_desc_t* b) {
    return a->min_filter == b->min_filter && a->mag_filter == b->mag_filter && a->mip_filter == b->mip_filter &&
           a->address_u == b->address_u && a->address_v == b->address_v && a->address_w == b->address_w;
}

md_gpu_sampler_t md_gpu_sampler(md_gpu_device_t dev, const md_gpu_sampler_desc_t* desc) {
    md_gpu_sampler_t h = {0};
    if (!dev) { md_mtl_fail("md_gpu_sampler: null device"); return h; }
    md_gpu_sampler_desc_t d;
    memset(&d, 0, sizeof(d));
    if (desc) d = *desc;

    md_mutex_lock(&dev->device_mutex);
    for (uint32_t i = 0; i < dev->sampler_count; ++i) {
        if (md_mtl_sampler_desc_eq(&dev->samplers[i].desc, &d)) {
            h.handle = dev->samplers[i].handle;
            md_mutex_unlock(&dev->device_mutex);
            return h;
        }
    }
    if (dev->sampler_count >= MD_MTL_MAX_SAMPLERS) {
        md_mutex_unlock(&dev->device_mutex);
        md_mtl_fail("md_gpu_sampler: more than %u distinct samplers", MD_MTL_MAX_SAMPLERS);
        return h;
    }

    static const MTLSamplerAddressMode modes[] = {
        MTLSamplerAddressModeClampToEdge,
        MTLSamplerAddressModeRepeat,
        MTLSamplerAddressModeMirrorRepeat,
    };
    @autoreleasepool {
        MTLSamplerDescriptor* sd = [[MTLSamplerDescriptor alloc] init];
        sd.minFilter    = d.min_filter == MD_GPU_FILTER_LINEAR ? MTLSamplerMinMagFilterLinear : MTLSamplerMinMagFilterNearest;
        sd.magFilter    = d.mag_filter == MD_GPU_FILTER_LINEAR ? MTLSamplerMinMagFilterLinear : MTLSamplerMinMagFilterNearest;
        sd.mipFilter    = d.mip_filter == MD_GPU_FILTER_LINEAR ? MTLSamplerMipFilterLinear    : MTLSamplerMipFilterNearest;
        sd.sAddressMode = modes[(unsigned)d.address_u % 3u];
        sd.tAddressMode = modes[(unsigned)d.address_v % 3u];
        sd.rAddressMode = modes[(unsigned)d.address_w % 3u];
        /* Required for gpuResourceID to be usable from an argument buffer. */
        sd.supportArgumentBuffers = YES;

        id<MTLSamplerState> smp = [dev->device newSamplerStateWithDescriptor:sd];
        MD_MTL_DROP_NEW(sd);
        if (smp) {
            MD_MTL_OWN(smp);
            md_mtl_sampler_entry_t* e = &dev->samplers[dev->sampler_count++];
            e->desc    = d;
            e->sampler = smp;
            e->handle  = (uint64_t)smp.gpuResourceID._impl;
            h.handle   = e->handle;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    if (!h.handle) md_mtl_fail("md_gpu_sampler: newSamplerStateWithDescriptor failed");
    return h;
}

/* ---- Texture copies ------------------------------------------------------------ */

typedef struct md_mtl_copy_region_t {
    MTLOrigin origin;          /* z unused for 2D arrays */
    MTLSize   size;            /* depth 1 for 2D arrays  */
    uint32_t  mip;
    uint32_t  first_layer;     /* 2D arrays only */
    uint32_t  layer_count;     /* 1 otherwise    */
    uint64_t  bytes_per_row;
    uint64_t  bytes_per_image; /* 0 for 2D / 2D-array slices, as Metal requires */
    uint64_t  layer_stride;    /* bytes between consecutive layers in the buffer */
    uint64_t  bytes;
} md_mtl_copy_region_t;

static bool md_mtl_resolve_region(md_gpu_texture_t t, const md_gpu_tex_region_t* r, md_mtl_copy_region_t* out, const char* what) {
    md_gpu_tex_region_t z;
    memset(&z, 0, sizeof(z));
    if (!r) r = &z;
    const md_gpu_texture_desc_t* d = &t->desc;
    if (r->mip >= d->mip_levels) return md_mtl_fail("%s: mip %u out of range (texture '%s' has %u)", what, r->mip, t->label, d->mip_levels);

    uint32_t dim[3];
    dim[0] = d->width  >> r->mip; if (!dim[0]) dim[0] = 1;
    dim[1] = d->height >> r->mip; if (!dim[1]) dim[1] = 1;
    if (d->type == MD_GPU_TEX_3D) { dim[2] = d->depth_or_layers >> r->mip; if (!dim[2]) dim[2] = 1; }
    else                          { dim[2] = d->depth_or_layers; }

    uint32_t ext[3];
    for (int i = 0; i < 3; ++i) {
        if (r->offset[i] >= dim[i]) {
            return md_mtl_fail("%s: offset[%d] = %u is outside texture '%s' (extent %u at mip %u)",
                               what, i, r->offset[i], t->label, dim[i], r->mip);
        }
        ext[i] = r->extent[i] ? r->extent[i] : dim[i] - r->offset[i];
        if (ext[i] > dim[i] - r->offset[i]) {
            return md_mtl_fail("%s: region [%u, +%u) on axis %d overruns texture '%s' (extent %u at mip %u)",
                               what, r->offset[i], ext[i], i, t->label, dim[i], r->mip);
        }
    }

    memset(out, 0, sizeof(*out));
    out->mip           = r->mip;
    out->bytes_per_row = (uint64_t)ext[0] * t->fi.bytes;
    if (d->type == MD_GPU_TEX_3D) {
        out->origin          = MTLOriginMake(r->offset[0], r->offset[1], r->offset[2]);
        out->size            = MTLSizeMake(ext[0], ext[1], ext[2]);
        out->bytes_per_image = out->bytes_per_row * ext[1];
        out->first_layer     = 0;
        out->layer_count     = 1;
    } else {
        out->origin          = MTLOriginMake(r->offset[0], r->offset[1], 0);
        out->size            = MTLSizeMake(ext[0], ext[1], 1);
        out->bytes_per_image = 0;
        out->first_layer     = d->type == MD_GPU_TEX_2D_ARRAY ? r->offset[2] : 0;
        out->layer_count     = d->type == MD_GPU_TEX_2D_ARRAY ? ext[2] : 1;
    }
    out->layer_stride = out->bytes_per_row * ext[1];
    out->bytes        = (uint64_t)ext[0] * ext[1] * ext[2] * t->fi.bytes;
    return true;
}

size_t md_gpu_texture_region_size(md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!t) return 0;
    md_mtl_copy_region_t cr;
    if (!md_mtl_resolve_region(t, region, &cr, "md_gpu_texture_region_size")) return 0;
    return (size_t)cr.bytes;
}

static bool md_mtl_check_buffer_offset(md_gpu_texture_t t, uint64_t off, const char* what) {
    const uint64_t a = t->fi.depth ? 4 : t->fi.bytes;
    if (off % a != 0) return md_mtl_fail("%s: buffer address must be %llu-byte aligned for %s", what, (unsigned long long)a, t->fi.name);
    return true;
}

static bool md_mtl_record_texture_copy(md_gpu_stream_t s, md_gpu_texture_t t, id<MTLBuffer> buffer, uint64_t offset,
                                       const md_mtl_copy_region_t* cr, bool to_texture) {
    if (!md_mtl_check_no_upload(s, "texture copy")) return false;
    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    /* Depth-stencil textures move their depth plane only. */
    const MTLBlitOption opt = t->fi.stencil ? MTLBlitOptionDepthFromDepthStencil : MTLBlitOptionNone;
    for (uint32_t l = 0; l < cr->layer_count; ++l) {
        const uint64_t boff = offset + (uint64_t)l * cr->layer_stride;
        const NSUInteger slice = cr->first_layer + l;
        if (to_texture) {
            [enc copyFromBuffer:buffer sourceOffset:boff
              sourceBytesPerRow:cr->bytes_per_row sourceBytesPerImage:cr->bytes_per_image
                     sourceSize:cr->size
                      toTexture:t->texture destinationSlice:slice destinationLevel:cr->mip
              destinationOrigin:cr->origin options:opt];
        } else {
            [enc copyFromTexture:t->texture sourceSlice:slice sourceLevel:cr->mip
                    sourceOrigin:cr->origin sourceSize:cr->size
                        toBuffer:buffer destinationOffset:boff
          destinationBytesPerRow:cr->bytes_per_row destinationBytesPerImage:cr->bytes_per_image
                         options:opt];
        }
    }
    md_mtl_did_op(s);
    return true;
}

bool md_gpu_copy_to_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, md_gpu_addr_t src) {
    if (!s || !t) return md_mtl_fail("md_gpu_copy_to_texture: null argument");
    md_mtl_copy_region_t cr;
    if (!md_mtl_resolve_region(t, region, &cr, "md_gpu_copy_to_texture")) return false;
    uint64_t off;
    md_mtl_block_t* b = md_mtl_resolve(s->device, src, cr.bytes, &off, "md_gpu_copy_to_texture");
    if (!b) return false;
    if (!md_mtl_check_buffer_offset(t, off, "md_gpu_copy_to_texture")) return false;
    return md_mtl_record_texture_copy(s, t, b->buffer, off, &cr, true);
}

bool md_gpu_copy_from_texture(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!s || !t) return md_mtl_fail("md_gpu_copy_from_texture: null argument");
    md_mtl_copy_region_t cr;
    if (!md_mtl_resolve_region(t, region, &cr, "md_gpu_copy_from_texture")) return false;
    uint64_t off;
    md_mtl_block_t* b = md_mtl_resolve(s->device, dst, cr.bytes, &off, "md_gpu_copy_from_texture");
    if (!b) return false;
    if (!md_mtl_check_buffer_offset(t, off, "md_gpu_copy_from_texture")) return false;
    return md_mtl_record_texture_copy(s, t, b->buffer, off, &cr, false);
}

bool md_gpu_upload_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, const void* src, size_t size) {
    if (!s || !t || !src) return md_mtl_fail("md_gpu_upload_texture: null argument");
    if (!md_mtl_check_no_upload(s, "md_gpu_upload_texture")) return false;
    md_mtl_copy_region_t cr;
    if (!md_mtl_resolve_region(t, region, &cr, "md_gpu_upload_texture")) return false;
    if ((uint64_t)size != cr.bytes) {
        return md_mtl_fail("md_gpu_upload_texture: region of texture '%s' is %llu bytes but %zu were given",
                           t->label, (unsigned long long)cr.bytes, size);
    }
    uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
    if (!md_mtl_arena_alloc(s, size, MD_MTL_ARG_ALIGN, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    return md_mtl_record_texture_copy(s, t, buf, off, &cr, true);
}

/* =========================================================================
   10. Kernels and launches
   ========================================================================= */

/* --- Obtaining an MTLLibrary ------------------------------------------------

   A kernel blob is whatever compile_gpu_shaders() embedded for it, and which of
   the two it is depends on the machine the build ran on:

     Apple's offline compiler present -> .metallib bytes -> newLibraryWithData:
     absent (no Xcode Metal toolchain) -> MSL source text -> newLibraryWithSource:

   The blob says which it is: every metallib starts with 'MTLB' and MSL source
   cannot. There is deliberately no fallback from one path to the other; both
   failures are real errors. Returns a +1 library the caller drops. */

static bool md_mtl_blob_is_metallib(const void* code, size_t size) {
    return size >= 4 && memcmp(code, "MTLB", 4) == 0;
}

static id<MTLLibrary> md_mtl_library_from_source(md_gpu_device_t dev, const void* code, size_t size, const char* who) {
    /* The embedded MSL has no terminator; copy it into one so the autoreleased
       +stringWithUTF8String: can be used, which also rejects non-UTF-8. */
    char* buf = (char*)md_alloc(dev->alloc, size + 1);
    if (!buf) { md_mtl_fail("out of memory preparing Metal source for '%s'", who); return nil; }
    memcpy(buf, code, size);
    buf[size] = '\0';
    NSString* src = [NSString stringWithUTF8String:buf];
    md_free(dev->alloc, buf, size + 1);
    if (!src) {
        md_mtl_fail("shader blob for '%s' is neither a metallib nor valid UTF-8 Metal source", who);
        return nil;
    }

    NSError* err = nil;
    id<MTLLibrary> lib = [dev->device newLibraryWithSource:src options:nil error:&err];
    if (!lib) {
        MD_LOG_ERROR("md_gpu: failed to compile Metal shaders at runtime.\n"
                     "Shader: %s\n\nMetal compiler error:\n%s",
                     who, err ? [[err localizedDescription] UTF8String] : "(no diagnostic returned)");
        md_mtl_fail("runtime Metal compilation of '%s' failed (compiler error above)", who);
        return nil;
    }
    MD_LOG_DEBUG("md_gpu: compiled '%s' from embedded MSL at runtime", who);
    return lib;
}

static id<MTLLibrary> md_mtl_library_from_blob(md_gpu_device_t dev, const void* code, size_t size, const char* label) {
    const char* who = label ? label : "kernel";
    if (!md_mtl_blob_is_metallib(code, size)) {
        return md_mtl_library_from_source(dev, code, size, who);
    }
    NSError* err = nil;
    dispatch_data_t data = dispatch_data_create(code, size, dispatch_get_main_queue(),
                                                DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> lib = [dev->device newLibraryWithData:data error:&err];
#if !__has_feature(objc_arc)
    dispatch_release(data);
#endif
    if (!lib) {
        md_mtl_fail("newLibraryWithData failed for '%s': %s", who,
                    err ? [[err localizedDescription] UTF8String] : "?");
        return nil;
    }
    return lib;
}

/* A kernel takes exactly one buffer -- the root -- so the first buffer binding
   reflection reports is the one wanted. */
static uint32_t md_mtl_reflect_arg_buffer_index(MTLComputePipelineReflection* refl, const char* label) {
    if (refl) {
        if (@available(macOS 13.0, iOS 16.0, *)) {
            for (id<MTLBinding> b in refl.bindings) {
                if (b.type == MTLBindingTypeBuffer) return (uint32_t)b.index;
            }
            return MD_MTL_ARG_BUFFER_INDEX;
        }
    }
    MD_LOG_DEBUG("md_gpu: no binding reflection for kernel '%s', assuming buffer(%d)",
                 label ? label : "kernel", MD_MTL_ARG_BUFFER_INDEX);
    return MD_MTL_ARG_BUFFER_INDEX;
}

/* Build a pipeline from a library entry point. Caller registers the kernel. */
static md_gpu_kernel_t md_mtl_kernel_from_library(md_gpu_device_t dev, id<MTLLibrary> lib, const char* entry,
                                                  const char* label, const uint32_t group_size[3], uint32_t args_size) {
    /* Metal reserves 'main', so Slang renames entry points to '<name>_0'. Try
       the requested name, then the Slang-mangled form, then the sole function. */
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:entry]];
    if (!fn) fn = [lib newFunctionWithName:[NSString stringWithFormat:@"%s_0", entry]];
    if (!fn && lib.functionNames.count == 1) fn = [lib newFunctionWithName:lib.functionNames[0]];
    if (!fn) { md_mtl_fail("entry point '%s' not found in shader library for '%s'", entry, label); return NULL; }

    NSError* err = nil;
    MTLComputePipelineReflection* refl = nil;
    id<MTLComputePipelineState> pso = [dev->device newComputePipelineStateWithFunction:fn
                                                                              options:MTLPipelineOptionBindingInfo
                                                                           reflection:&refl
                                                                                error:&err];
    MD_MTL_DROP_NEW(fn);   /* the pipeline keeps what it needs */
    if (!pso) {
        md_mtl_fail("kernel '%s': newComputePipelineStateWithFunction failed: %s", label,
                    err ? [[err localizedDescription] UTF8String] : "?");
        return NULL;
    }
    MD_MTL_OWN(pso);

    const uint64_t threads = (uint64_t)group_size[0] * group_size[1] * group_size[2];
    if (threads > (uint64_t)pso.maxTotalThreadsPerThreadgroup) {
        md_mtl_fail("kernel '%s': %llu threads per group exceeds this pipeline's limit of %llu",
                    label, (unsigned long long)threads, (unsigned long long)pso.maxTotalThreadsPerThreadgroup);
        MD_MTL_RELEASE(pso);
        return NULL;
    }

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) { MD_MTL_RELEASE(pso); md_mtl_fail("out of memory"); return NULL; }
    memset(k, 0, sizeof(*k));
    k->device           = dev;
    k->pso              = pso;
    k->group_size[0]    = group_size[0];
    k->group_size[1]    = group_size[1];
    k->group_size[2]    = group_size[2];
    k->args_size        = args_size;
    k->arg_buffer_index = md_mtl_reflect_arg_buffer_index(refl, label);
    snprintf(k->label, sizeof(k->label), "%s", label);
    return k;
}

static void md_mtl_kernel_free(md_gpu_device_t dev, md_gpu_kernel_t k) {
    MD_MTL_RELEASE(k->pso);
    md_free(dev->alloc, k, sizeof(*k));
}

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t dev, const md_gpu_kernel_desc_t* desc) {
    if (!dev || !desc || !desc->code || desc->code_size == 0) {
        md_mtl_fail("md_gpu_kernel_create: missing code");
        return NULL;
    }
    const char* label = desc->label ? desc->label : "kernel";
    /* Metal cannot recover the threadgroup size from a library, and a silent
       {1,1,1} dispatches a fraction of the threads. So it is required. */
    if (desc->group_size[0] == 0 || desc->group_size[1] == 0 || desc->group_size[2] == 0) {
        md_mtl_fail("kernel '%s': group_size is required ({%u,%u,%u} given); use the generated kernel descriptor",
                    label, desc->group_size[0], desc->group_size[1], desc->group_size[2]);
        return NULL;
    }

    md_gpu_kernel_t k = NULL;
    @autoreleasepool {
        id<MTLLibrary> lib = md_mtl_library_from_blob(dev, desc->code, desc->code_size, label);
        if (lib) {
            k = md_mtl_kernel_from_library(dev, lib, desc->entry_point ? desc->entry_point : "main",
                                           label, desc->group_size, desc->args_size);
            MD_MTL_DROP_NEW(lib);
        }
    }
    if (!k) return NULL;

    md_mutex_lock(&dev->device_mutex);
    md_gpu_kernel_t* slot = (md_gpu_kernel_t*)md_mtl_vec_push(&dev->kernels, dev->alloc);
    if (slot) *slot = k;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { md_mtl_kernel_free(dev, k); md_mtl_fail("out of memory"); return NULL; }
    return k;
}

/* Immediate: a command buffer retains the pipeline states set on its encoders,
   so work in flight keeps what it uses alive. */
void md_gpu_kernel_destroy(md_gpu_kernel_t k) {
    if (!k) return;
    md_gpu_device_t dev = k->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->kernels, k);
    md_mutex_unlock(&dev->device_mutex);
    md_mtl_kernel_free(dev, k);
}

bool md_gpu_kernel_info(md_gpu_kernel_t k, md_gpu_kernel_info_t* info) {
    if (!k || !info) return false;
    memset(info, 0, sizeof(*info));
    info->group_size[0]            = k->group_size[0];
    info->group_size[1]            = k->group_size[1];
    info->group_size[2]            = k->group_size[2];
    info->args_size                = k->args_size;
    info->max_threads_per_group    = (uint32_t)k->pso.maxTotalThreadsPerThreadgroup;
    info->preferred_group_multiple = (uint32_t)k->pso.threadExecutionWidth;
    return true;
}

md_gpu_grid_t md_gpu_grid_for(md_gpu_kernel_t k, uint32_t nx, uint32_t ny, uint32_t nz) {
    md_gpu_grid_t g = {0, 0, 0};
    if (!k) return g;
    g.x = (uint32_t)(((uint64_t)nx + k->group_size[0] - 1) / k->group_size[0]);
    g.y = (uint32_t)(((uint64_t)ny + k->group_size[1] - 1) / k->group_size[1]);
    g.z = (uint32_t)(((uint64_t)nz + k->group_size[2] - 1) / k->group_size[2]);
    return g;
}

/* What the root buffer holds.

   Slang lowers `struct Root { Args* args; }; ConstantBuffer<Root> root;` to a
   kernel parameter `Root_0 constant* root_0 [[buffer(N)]]` whose one member is
   a device pointer: the bound buffer holds an 8-byte *address*, and the
   argument struct lives wherever that points. So the argument block goes into
   the arena, its address into an 8-byte root cell of its own, and the cell is
   what gets bound. (setBytes: for the 8 bytes proved unreliable when called
   repeatedly at one index in a single encoder.)

   Every dispatch rebinds its pipeline and root cell unconditionally; caching
   either saves nothing measurable and fails silently when wrong. */
static bool md_mtl_launch_common(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid,
                                 const void* args, size_t args_size,
                                 id<MTLBuffer> indirect_buf, uint64_t indirect_off, bool is_indirect,
                                 bool internal) {
    if (k->device != s->device) return md_mtl_fail("kernel '%s' belongs to a different device", k->label);
    if (s->kind == MD_GPU_STREAM_TRANSFER && !internal) {
        return md_mtl_fail("kernel '%s' launched into '%s', a transfer stream; kernels run on compute streams", k->label, s->label);
    }
    if (k->args_size != 0 && args_size != k->args_size) {
        return md_mtl_fail("kernel '%s' expects a %u-byte argument struct but %zu bytes were passed", k->label, k->args_size, args_size);
    }
    if (args_size > 0 && !args) return md_mtl_fail("kernel '%s': null args with non-zero size", k->label);
    if (!md_mtl_check_no_upload(s, "launch")) return false;

    id<MTLBuffer> root_buf = nil;
    uint64_t root_off = 0;
    if (args_size > 0) {
        uint64_t arg_addr, root_addr;
        void* arg_host; void* root_host;
        if (!md_mtl_arena_alloc(s, args_size, MD_MTL_ARG_ALIGN, &arg_addr, &arg_host, NULL, NULL)) return false;
        memcpy(arg_host, args, args_size);
        if (!md_mtl_arena_alloc(s, sizeof(uint64_t), MD_MTL_ROOT_ALIGN, &root_addr, &root_host, &root_buf, &root_off)) return false;
        memcpy(root_host, &arg_addr, sizeof(arg_addr));
    }

    id<MTLComputeCommandEncoder> enc = md_mtl_compute_encoder(s);
    if (!enc) return false;
    md_mtl_declare_residency(s, enc);
    [enc setComputePipelineState:k->pso];
    if (root_buf) [enc setBuffer:root_buf offset:root_off atIndex:k->arg_buffer_index];
    MTLSize tg = MTLSizeMake(k->group_size[0], k->group_size[1], k->group_size[2]);
    if (is_indirect) {
        [enc dispatchThreadgroupsWithIndirectBuffer:indirect_buf indirectBufferOffset:indirect_off threadsPerThreadgroup:tg];
    } else {
        [enc dispatchThreadgroups:MTLSizeMake(grid.x, grid.y, grid.z) threadsPerThreadgroup:tg];
    }
    md_mtl_did_op(s);
    return true;
}

bool md_gpu_launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid, const void* args, size_t args_size) {
    if (!s || !k) return md_mtl_fail("md_gpu_launch: null stream or kernel");
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;
    return md_mtl_launch_common(s, k, grid, args, args_size, nil, 0, false, false);
}

bool md_gpu_launch_indirect(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_addr_t grid, const void* args, size_t args_size) {
    if (!s || !k || !grid) return md_mtl_fail("md_gpu_launch_indirect: null argument");
    uint64_t off;
    md_mtl_block_t* b = md_mtl_resolve(s->device, grid, 3 * sizeof(uint32_t), &off, "md_gpu_launch_indirect");
    if (!b) return false;
    if (off % 4 != 0) return md_mtl_fail("md_gpu_launch_indirect: grid address must be 4-byte aligned");
    return md_mtl_launch_common(s, k, md_gpu_grid(1, 1, 1), args, args_size, b->buffer, off, true, false);
}

/* Mirrors MdMakeGridArgs in md_gpu_builtin_msl.inl. */
typedef struct md_mtl_make_grid_args_t {
    md_gpu_addr_t count;
    md_gpu_addr_t out_grid;
    md_gpu_uint4  local;     /* xyz = threads per group, w unused */
} md_mtl_make_grid_args_t;

bool md_gpu_make_grid(md_gpu_stream_t s, md_gpu_addr_t out_grid, md_gpu_addr_t count, md_gpu_kernel_t k) {
    if (!s || !out_grid || !count || !k) return md_mtl_fail("md_gpu_make_grid: null argument");
    if (!md_mtl_resolve(s->device, out_grid, 3 * sizeof(uint32_t), NULL, "md_gpu_make_grid (out_grid)")) return false;
    if (!md_mtl_resolve(s->device, count, sizeof(uint32_t), NULL, "md_gpu_make_grid (count)")) return false;
    md_mtl_make_grid_args_t a;
    memset(&a, 0, sizeof(a));
    a.count    = count;
    a.out_grid = out_grid;
    a.local.x  = k->group_size[0];
    a.local.y  = k->group_size[1];
    a.local.z  = k->group_size[2];
    return md_gpu_launch(s, s->device->make_grid_kernel, md_gpu_grid(1, 1, 1), &a, sizeof(a));
}

/* Mirrors MdByteOpArgs in md_gpu_builtin_msl.inl. src == 0 means fill. */
typedef struct md_mtl_byte_op_args_t {
    md_gpu_addr_t dst;
    md_gpu_addr_t src;
    uint32_t      size;
    uint32_t      value;
} md_mtl_byte_op_args_t;

#define MD_MTL_BYTE_OP_GROUP 64u
#define MD_MTL_BYTE_OP_CHUNK (1ull << 30)

/* Byte-granular copy or fill through a compute dispatch, for what blits may
   not do on macOS. Internal: allowed on transfer streams too, since a Metal
   queue runs compute encoders whatever md_gpu calls the stream. */
static bool md_mtl_byte_op(md_gpu_stream_t s, uint64_t dst, uint64_t src, uint64_t size, uint8_t value) {
    md_gpu_kernel_t k = s->device->byte_op_kernel;
    if (!k) return md_mtl_fail("byte copy kernel unavailable");
    for (uint64_t done = 0; done < size;) {
        const uint64_t n = (size - done) < MD_MTL_BYTE_OP_CHUNK ? (size - done) : MD_MTL_BYTE_OP_CHUNK;
        md_mtl_byte_op_args_t a;
        a.dst   = dst + done;
        a.src   = src ? src + done : 0;
        a.size  = (uint32_t)n;
        a.value = value;
        const uint64_t threads = (n + MD_MTL_BYTE_OP_SPAN - 1) / MD_MTL_BYTE_OP_SPAN;
        const uint32_t groups  = (uint32_t)((threads + MD_MTL_BYTE_OP_GROUP - 1) / MD_MTL_BYTE_OP_GROUP);
        if (!md_mtl_launch_common(s, k, md_gpu_grid(groups, 1, 1), &a, sizeof(a), nil, 0, false, true)) return false;
        done += n;
    }
    return true;
}

/* =========================================================================
   11. Host callbacks and polling
   ========================================================================= */

bool md_gpu_sync_on_complete(md_gpu_device_t dev, md_gpu_sync_t sync, md_gpu_host_fn fn, void* user) {
    if (!dev || !fn) return md_mtl_fail("md_gpu_sync_on_complete: null argument");
    if (md_gpu_sync_is_valid(sync) && sync.stream->device != dev) return md_mtl_fail("md_gpu_sync_on_complete: sync from another device");
    md_mutex_lock(&dev->device_mutex);
    md_mtl_hostfn_t* h = (md_mtl_hostfn_t*)md_mtl_vec_push(&dev->hostfns, dev->alloc);
    if (h) { h->sync = sync; h->fn = fn; h->user = user; }
    md_mutex_unlock(&dev->device_mutex);
    return h ? true : md_mtl_fail("out of memory");
}

bool md_gpu_launch_host_fn(md_gpu_stream_t s, md_gpu_host_fn fn, void* user) {
    if (!s || !fn) return md_mtl_fail("md_gpu_launch_host_fn: null argument");
    md_gpu_sync_t sync = md_gpu_stream_record(s);
    return md_gpu_sync_on_complete(s->device, sync, fn, user);
}

/* Timeline values sampled once per poll pass, so that callbacks fire in
   registration order. See the fuller note in md_gpu_vulkan.c. */
typedef struct md_mtl_snapshot_t {
    md_gpu_stream_t stream;
    uint64_t        completed;
} md_mtl_snapshot_t;

static bool md_mtl_snapshot_complete(md_gpu_device_t dev, md_mtl_vec_t* snap, md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    for (size_t i = 0; i < snap->count; ++i) {
        md_mtl_snapshot_t* e = &MD_MTL_VEC_AT(*snap, md_mtl_snapshot_t, i);
        if (e->stream == sync.stream) return e->completed >= sync.value;
    }
    uint64_t completed = md_mtl_stream_completed(sync.stream);
    md_mtl_snapshot_t* e = (md_mtl_snapshot_t*)md_mtl_vec_push(snap, dev->alloc);
    if (e) { e->stream = sync.stream; e->completed = completed; }
    return completed >= sync.value;
}

uint32_t md_gpu_device_poll(md_gpu_device_t dev) {
    if (!dev) return 0;
    uint32_t fired = 0;
    md_mtl_vec_t snap;
    md_mtl_vec_init(&snap, sizeof(md_mtl_snapshot_t));
    for (;;) {
        md_mtl_hostfn_t ready;
        bool have = false;
        md_mutex_lock(&dev->device_mutex);
        for (size_t i = 0; i < dev->hostfns.count; ++i) {
            md_mtl_hostfn_t* h = &MD_MTL_VEC_AT(dev->hostfns, md_mtl_hostfn_t, i);
            if (!md_mtl_snapshot_complete(dev, &snap, h->sync)) continue;
            ready = *h;
            md_mtl_vec_remove(&dev->hostfns, i);
            have = true;
            break;
        }
        md_mutex_unlock(&dev->device_mutex);
        if (!have) break;
        ready.fn(ready.user);
        fired++;
    }
    md_mtl_vec_free(&snap, dev->alloc);

    md_mutex_lock(&dev->device_mutex);
    md_mtl_process_retires_locked(dev, false);
    for (size_t i = 0; i < dev->pools.count; ++i) {
        md_gpu_pool_t p = MD_MTL_VEC_AT(dev->pools, md_gpu_pool_t, i);
        if (p->cache_limit != 0) md_mtl_pool_trim_locked(p, p->cache_limit);
    }
    md_mutex_unlock(&dev->device_mutex);
    return fired;
}

/* =========================================================================
   12. Device
   ========================================================================= */

#include "md_gpu_builtin_msl.inl"

static bool md_mtl_create_builtin_kernels(md_gpu_device_t dev) {
    @autoreleasepool {
        /* Hand-written MSL, never through slangc, so always compiled from source. */
        id<MTLLibrary> lib = md_mtl_library_from_blob(dev, md_gpu_make_grid_msl, strlen(md_gpu_make_grid_msl), "md_gpu make_grid");
        if (!lib) return false;
        const uint32_t gs[3] = {1, 1, 1};
        dev->make_grid_kernel = md_mtl_kernel_from_library(dev, lib, "md_gpu_make_grid", "md_gpu make_grid",
                                                           gs, (uint32_t)sizeof(md_mtl_make_grid_args_t));
        MD_MTL_DROP_NEW(lib);

        id<MTLLibrary> blib = md_mtl_library_from_blob(dev, md_gpu_byte_op_msl, strlen(md_gpu_byte_op_msl), "md_gpu byte_op");
        if (!blib) return false;
        const uint32_t bgs[3] = {MD_MTL_BYTE_OP_GROUP, 1, 1};
        dev->byte_op_kernel = md_mtl_kernel_from_library(dev, blib, "md_gpu_byte_op", "md_gpu byte_op",
                                                         bgs, (uint32_t)sizeof(md_mtl_byte_op_args_t));
        MD_MTL_DROP_NEW(blib);
    }
    return dev->make_grid_kernel != NULL && dev->byte_op_kernel != NULL;
}

md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc) {
    md_mtl_has_error = false;
    struct md_allocator_i* alloc = (desc && desc->alloc) ? desc->alloc : md_get_heap_allocator();

    /* Metal API validation cannot be switched on from here: Metal reads
       METAL_DEVICE_WRAPPER_TYPE once, before the first Metal call. So say so,
       rather than letting the flag look like it did something. */
    if (desc && desc->enable_validation) {
        const char* wrapper = getenv("METAL_DEVICE_WRAPPER_TYPE");
        if (!wrapper || wrapper[0] == '0') {
            MD_LOG_INFO("md_gpu: validation requested but Metal API validation is a launch-time "
                        "setting. Re-run with METAL_DEVICE_WRAPPER_TYPE=1 (and "
                        "MTL_DEBUG_LAYER=1 for the shader validation layer).");
        }
    }

    md_gpu_device_t dev = (md_gpu_device_t)md_alloc(alloc, sizeof(md_gpu_device));
    if (!dev) { md_mtl_fail("out of memory"); return NULL; }
    memset(dev, 0, sizeof(*dev));
    dev->alloc = alloc;
    md_mtl_vec_init(&dev->registry, sizeof(md_mtl_block_t*));
    md_mtl_vec_init(&dev->pools,    sizeof(md_gpu_pool_t));
    md_mtl_vec_init(&dev->kernels,  sizeof(md_gpu_kernel_t));
    md_mtl_vec_init(&dev->streams,  sizeof(md_gpu_stream_t));
    md_mtl_vec_init(&dev->hostfns,  sizeof(md_mtl_hostfn_t));
    md_mtl_vec_init(&dev->retires,  sizeof(md_mtl_retire_t));
    md_mtl_vec_init(&dev->live_res, sizeof(id));
    md_mutex_init(&dev->queue_mutex);
    md_mutex_init(&dev->device_mutex);

    @autoreleasepool {
        id<MTLDevice> mtl = MTLCreateSystemDefaultDevice();
        if (!mtl) {
            md_mtl_fail("MTLCreateSystemDefaultDevice returned nil");
        } else {
            MD_MTL_OWN(mtl);
            dev->device      = mtl;
            dev->is_discrete = ![mtl hasUnifiedMemory];

            /* Device-wide residency, where available. Reached dynamically so
               the file still builds against SDKs that predate it. */
            Class rsd = NSClassFromString(@"MTLResidencySetDescriptor");
            if (rsd && [mtl respondsToSelector:@selector(newResidencySetWithDescriptor:error:)]) {
                id descriptor = [[rsd alloc] init];
                /* Returned as a raw +1 pointer so ARC never sees it: ownership
                   is held manually and dropped with MD_MTL_RELEASE on destroy. */
                typedef void* (*md_mtl_new_rs_t)(id, SEL, id, void*);
                void* set = ((md_mtl_new_rs_t)objc_msgSend)(mtl, @selector(newResidencySetWithDescriptor:error:),
                                                            descriptor, NULL);
                MD_MTL_DROP_NEW(descriptor);
                if (set) {
                    dev->residency_set     = (__bridge id)set;
                    dev->has_residency_set = true;
                }
            }
            if (!dev->has_residency_set) {
                MD_LOG_DEBUG("md_gpu: MTLResidencySet unavailable, falling back to per-encoder useResources");
            }
        }
    }
    if (!dev->device) goto fail;

    dev->default_compute  = md_mtl_stream_create_internal(dev, MD_GPU_STREAM_COMPUTE,  "default compute",  true);
    dev->default_transfer = md_mtl_stream_create_internal(dev, MD_GPU_STREAM_TRANSFER, "default transfer", true);
    if (!dev->default_compute || !dev->default_transfer) goto fail;
    if (!md_mtl_create_builtin_kernels(dev)) goto fail;

    MD_LOG_DEBUG("md_gpu: device '%s'", [[dev->device name] UTF8String]);
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
    info->preferred_group_multiple = dev->make_grid_kernel ? (uint32_t)dev->make_grid_kernel->pso.threadExecutionWidth : 32;
    snprintf(info->name, sizeof(info->name), "%s", [[dev->device name] UTF8String]);
    return true;
}

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;

    /* Flush and idle every stream, then fire what is pending. */
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
        s->upload_open = false;
        md_mtl_stream_submit(s);
    }
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
        if (s->submitted_value > 0) md_mtl_event_wait(s->timeline, s->submitted_value);
    }
    md_gpu_device_poll(dev);

    /* Everything the caller did not destroy, the device does. */
    while (dev->pools.count > 0)   md_gpu_pool_destroy(MD_MTL_VEC_AT(dev->pools, md_gpu_pool_t, 0));
    while (dev->kernels.count > 0) md_gpu_kernel_destroy(MD_MTL_VEC_AT(dev->kernels, md_gpu_kernel_t, 0));
    if (dev->make_grid_kernel) md_mtl_kernel_free(dev, dev->make_grid_kernel);
    if (dev->byte_op_kernel)   md_mtl_kernel_free(dev, dev->byte_op_kernel);

    md_mutex_lock(&dev->device_mutex);
    md_mtl_process_retires_locked(dev, true);
    md_mutex_unlock(&dev->device_mutex);

    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_mtl_stream_free(MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i));
    }
    md_mtl_vec_free(&dev->streams, alloc);
    md_mtl_vec_free(&dev->pools,   alloc);
    md_mtl_vec_free(&dev->kernels, alloc);

    for (uint32_t i = 0; i < dev->sampler_count; ++i) MD_MTL_RELEASE(dev->samplers[i].sampler);
    md_mtl_vec_free(&dev->hostfns,  alloc);
    md_mtl_vec_free(&dev->retires,  alloc);
    md_mtl_vec_free(&dev->registry, alloc);
    md_mtl_vec_free(&dev->live_res, alloc);
    if (dev->residency_set) MD_MTL_RELEASE(dev->residency_set);
    md_mutex_destroy(&dev->queue_mutex);
    md_mutex_destroy(&dev->device_mutex);
    if (dev->device) MD_MTL_RELEASE(dev->device);
    md_free(alloc, dev, sizeof(*dev));
}

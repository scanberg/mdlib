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

     A render pass is one render command encoder, fenced like the others
     (updateFence:afterStages: / waitForFence:beforeStages:), so it too is
     ordered after everything before it and before everything after it.

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
#include "md_gpu_tlsf.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>
#include <TargetConditionals.h>

#include "md_gpu_select.inl"
#include <TargetConditionals.h>
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
#define MD_MTL_HEAP_ALIGN        256u
#define MD_MTL_HEAP_CHUNK_MIN    (4ull << 20)
#define MD_MTL_HEAP_CHUNK_MAX    (64ull << 20)
#define MD_MTL_HEAP_LARGE_ALIGN  (64ull << 10)
#define MD_MTL_HEAP_CACHE_DEFAULT (256ull << 20)
#define MD_MTL_TEMP_CHUNK_MIN    (4ull << 20)
#define MD_MTL_TEMP_KINDS           2u
/* setBuffer:offset: for the constant address space wants 256-byte offsets on
   macOS (Mac2 family GPUs); the root cell is bound that way, so it gets them. */
#define MD_MTL_ROOT_ALIGN        256u
/* Bytes each thread of the built-in byte kernel handles. */
#define MD_MTL_BYTE_OP_SPAN       16u
#define MD_MTL_MAX_SAMPLERS      256u
#define MD_MTL_ERROR_BUF         2560u

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

/* One MTLBuffer: a region of a heap, or a chunk of a stream's temp arena. */
typedef struct md_mtl_chunk_t {
    __unsafe_unretained id<MTLBuffer> buffer;
    uint64_t          address;
    uint8_t*          host;
    uint64_t          size;
    md_gpu_mem_kind_t kind;
    bool              empty;
    md_tlsf_node_t*   empty_node;
    uint64_t          cursor;
    uint32_t          owner_depth;
    uint64_t          retire_value;
    struct md_mtl_range_t* range;
} md_mtl_chunk_t;

/* Registry entry: a live md_gpu_malloc allocation, or a whole temp chunk. */
typedef struct md_mtl_range_t {
    uint64_t          address;
    uint64_t          size;
    md_mtl_chunk_t*   chunk;
    uint64_t          chunk_offset;
    md_tlsf_node_t*   node;          /* NULL for a temp chunk */
    md_gpu_mem_kind_t kind;
} md_mtl_range_t;

typedef struct md_mtl_span_t {
    __unsafe_unretained id<MTLBuffer> buffer;
    uint64_t offset;
    uint8_t* host;
} md_mtl_span_t;

typedef struct md_mtl_heap_t {
    md_tlsf_t    tlsf;
    md_mtl_vec_t chunks;             /* md_mtl_chunk_t* */
    uint64_t     reserved;
    uint64_t     empty_bytes;
    uint64_t     in_use;
    uint64_t     peak_in_use;
    uint32_t     allocations;
    uint64_t     temp_bytes;
} md_mtl_heap_t;

typedef struct md_mtl_pending_free_t {
    md_tlsf_node_t*   node;
    md_gpu_mem_kind_t kind;
    md_gpu_stream_t   stream;
    uint64_t          value;
} md_mtl_pending_free_t;

typedef struct md_mtl_temp_arena_t {
    md_mtl_vec_t active;             /* md_mtl_chunk_t* */
    md_mtl_vec_t spare;              /* md_mtl_chunk_t* */
} md_mtl_temp_arena_t;

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
    __unsafe_unretained id<MTLRenderCommandEncoder>  render_enc;    /* the open pass, if any */
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

    md_mtl_temp_arena_t temp[MD_MTL_TEMP_KINDS];
    uint32_t          temp_depth;

    bool              upload_open;
    bool              upload_direct;
    md_gpu_addr_t     upload_dst;
    uint64_t          upload_src_addr;
    size_t            upload_size;

    /* The open render pass, if any. */
    bool              in_pass;
    uint32_t          pass_color_count;
    md_gpu_format_t   pass_color[MD_GPU_MAX_COLOR_TARGETS];
    md_gpu_format_t   pass_depth;
    uint32_t          pass_width, pass_height;
    struct md_gpu_pipeline* bound_pipeline;
    md_gpu_draw_state_t draw_state;   /* as set on render_enc */

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
    bool                  resident;          /* in the residency set / live_res */
    bool                  external;          /* a drawable's texture; not the caller's */
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

typedef struct md_gpu_pipeline {
    md_gpu_device_t   device;
    __unsafe_unretained id<MTLRenderPipelineState> pso;
    MTLPrimitiveType  primitive;
    uint32_t          vs_root;          /* buffer index of the root, or UINT32_MAX */
    uint32_t          fs_root;
    uint32_t          color_count;
    md_gpu_format_t   color[MD_GPU_MAX_COLOR_TARGETS];
    md_gpu_format_t   depth;
    uint32_t          args_size;
    char              label[64];
} md_gpu_pipeline;

typedef struct md_gpu_surface {
    md_gpu_device_t        device;
    __unsafe_unretained CAMetalLayer* layer;
    md_gpu_surface_desc_t  desc;
    char                   label[64];
    uint32_t               want_width, want_height;
    uint32_t               width, height;       /* drawableSize last set */
    uint64_t               frame;
    /* Between acquire and present. */
    __unsafe_unretained id<CAMetalDrawable> drawable;
    md_gpu_texture_t       texture;
    md_gpu_stream_t        acquired_stream;
} md_gpu_surface;

typedef struct md_mtl_hostfn_t {
    md_gpu_sync_t  sync;
    md_gpu_host_fn fn;
    void*          user;
} md_mtl_hostfn_t;

typedef enum md_mtl_retire_kind_t {
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

    md_mtl_vec_t registry;   /* md_mtl_range_t*, sorted by address */
    md_mtl_heap_t heaps[MD_GPU_MEM_KIND_COUNT];
    md_mtl_vec_t pending_frees;   /* md_mtl_pending_free_t */
    uint64_t     heap_cache_limit;
    md_mtl_vec_t textures;   /* md_gpu_texture_t, live */
    uint64_t     texture_bytes;
    md_mtl_vec_t kernels;    /* md_gpu_kernel_t */
    md_mtl_vec_t pipelines;  /* md_gpu_pipeline_t */
    md_mtl_vec_t surfaces;   /* md_gpu_surface_t */
    md_mtl_vec_t streams;    /* md_gpu_stream_t */
    md_mtl_vec_t hostfns;    /* md_mtl_hostfn_t */
    md_mtl_vec_t retires;    /* md_mtl_retire_t */
    md_mtl_vec_t live_res;   /* id<MTLResource>, for the no-residency-set path */
    uint64_t     res_gen;    /* bumped whenever live_res gains a resource */

    md_gpu_stream_t default_compute;
    md_gpu_stream_t default_transfer;
    md_gpu_stream_t default_graphics;

    /* Depth state is an object in Metal; one per (compare, write), made on
       first use. */
    __unsafe_unretained id<MTLDepthStencilState> depth_states[8][2];
    bool            can_present;

    md_mtl_sampler_entry_t samplers[MD_MTL_MAX_SAMPLERS];
    uint32_t               sampler_count;

    md_gpu_kernel_t make_grid_kernel;
    md_gpu_kernel_t byte_op_kernel;   /* unaligned copies and fills */
    bool            is_discrete;
    uint32_t        adapter_index;
} md_gpu_device;

static bool     md_mtl_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t align, uint64_t* out_addr, void** out_host, id<MTLBuffer>* out_buf, uint64_t* out_off);
static bool     md_mtl_stream_submit(md_gpu_stream_t s);
static bool     md_mtl_byte_op(md_gpu_stream_t s, uint64_t dst, uint64_t src, uint64_t size, uint8_t value);
static uint64_t md_mtl_stream_completed(md_gpu_stream_t s);
static void     md_mtl_heap_release_node_locked(md_gpu_device_t dev, md_gpu_mem_kind_t kind, md_tlsf_node_t* node);
static void     md_mtl_temp_free_all(md_gpu_stream_t s);
static void     md_mtl_texture_free(md_gpu_device_t dev, md_gpu_texture_t t);
static void     md_mtl_abandon_pass(md_gpu_stream_t s);

/* =========================================================================
   3. Allocation registry
   ========================================================================= */

/* Caller holds device_mutex. */
static md_mtl_range_t* md_mtl_registry_find_locked(md_gpu_device_t dev, uint64_t address) {
    size_t lo = 0, hi = dev->registry.count;
    md_mtl_range_t** arr = (md_mtl_range_t**)dev->registry.data;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        md_mtl_range_t* r = arr[mid];
        if (address < r->address)                 hi = mid;
        else if (address >= r->address + r->size) lo = mid + 1;
        else                                      return r;
    }
    return NULL;
}

static size_t md_mtl_registry_lower_bound(md_gpu_device_t dev, uint64_t address) {
    size_t lo = 0, hi = dev->registry.count;
    md_mtl_range_t** arr = (md_mtl_range_t**)dev->registry.data;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (arr[mid]->address < address) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static bool md_mtl_registry_insert_locked(md_gpu_device_t dev, md_mtl_range_t* r) {
    if (!md_mtl_vec_reserve(&dev->registry, dev->alloc, dev->registry.count + 1)) return false;
    md_mtl_range_t** arr = (md_mtl_range_t**)dev->registry.data;
    const size_t at = md_mtl_registry_lower_bound(dev, r->address);
    memmove(arr + at + 1, arr + at, (dev->registry.count - at) * sizeof(*arr));
    arr[at] = r;
    dev->registry.count++;
    return true;
}

static void md_mtl_registry_remove_locked(md_gpu_device_t dev, md_mtl_range_t* r) {
    const size_t at = md_mtl_registry_lower_bound(dev, r->address);
    if (at < dev->registry.count && MD_MTL_VEC_AT(dev->registry, md_mtl_range_t*, at) == r) md_mtl_vec_remove(&dev->registry, at);
}

/* Resolve [addr, addr + size), under the device lock -- malloc on another
   thread mutates the registry. */
static bool md_mtl_resolve(md_gpu_device_t dev, md_gpu_addr_t addr, uint64_t size,
                           md_mtl_span_t* out, const char* what) {
    md_mutex_lock(&dev->device_mutex);
    md_mtl_range_t* r = md_mtl_registry_find_locked(dev, addr);
    md_mtl_range_t copy;
    if (r) copy = *r;
    md_mutex_unlock(&dev->device_mutex);
    if (!r) return md_mtl_fail("%s: 0x%llx is not a live md_gpu allocation", what, (unsigned long long)addr);
    const uint64_t off = addr - copy.address;
    if (size > copy.size || off > copy.size - size) {
        return md_mtl_fail("%s: range [0x%llx, +%llu) overruns its %llu-byte allocation",
                           what, (unsigned long long)addr, (unsigned long long)size, (unsigned long long)copy.size);
    }
    if (out) {
        out->buffer = copy.chunk->buffer;
        out->offset = copy.chunk_offset + off;
        out->host   = copy.chunk->host ? copy.chunk->host + copy.chunk_offset + off : NULL;
    }
    return true;
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
       device address, so a page is made resident exactly like a heap chunk. */
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
        (void)kind;
        md_mtl_texture_free(dev, (md_gpu_texture_t)object);
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
        md_mtl_texture_free(dev, (md_gpu_texture_t)e.object);
        if (e.waits) md_free(dev->alloc, e.waits, e.wait_capacity * sizeof(md_mtl_wait_t));
    }
}

/* A stream is going away after being synchronised: drop every reference to
   it. Caller holds device_mutex. */
static void md_mtl_forget_stream_locked(md_gpu_device_t dev, md_gpu_stream_t s) {
    for (size_t i = 0; i < dev->pending_frees.count;) {
        md_mtl_pending_free_t pf = MD_MTL_VEC_AT(dev->pending_frees, md_mtl_pending_free_t, i);
        if (pf.stream == s) {
            md_mtl_vec_remove(&dev->pending_frees, i);
            md_mtl_heap_release_node_locked(dev, pf.kind, pf.node);
        } else {
            ++i;
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
    if (!s->compute_enc && !s->blit_enc && !s->render_enc) return;
    if (signal_fence && s->fence) {
        if (s->compute_enc)     [s->compute_enc updateFence:s->fence];
        else if (s->blit_enc)   [s->blit_enc    updateFence:s->fence];
        else                    [s->render_enc  updateFence:s->fence afterStages:MTLRenderStageVertex | MTLRenderStageFragment];
        s->fence_valid = true;
    }
    if (s->compute_enc) { [s->compute_enc endEncoding]; MD_MTL_RELEASE(s->compute_enc); s->compute_enc = nil; }
    if (s->blit_enc)    { [s->blit_enc    endEncoding]; MD_MTL_RELEASE(s->blit_enc);    s->blit_enc    = nil; }
    if (s->render_enc)  { [s->render_enc  endEncoding]; MD_MTL_RELEASE(s->render_enc);  s->render_enc  = nil; }
}

/* An operation may not be recorded while an upload is open: the staged copy
   the upload ends with must land where the upload began. */
static bool md_mtl_check_no_upload(md_gpu_stream_t s, const char* what) {
    if (s->upload_open) return md_mtl_fail("%s: stream '%s' has an open upload; call md_gpu_upload_end first", what, s->label);
    return true;
}

/* Inside a render pass only draws and dynamic state may be recorded. */
static bool md_mtl_not_in_pass(md_gpu_stream_t s, const char* what) {
    if (!s->in_pass) return true;
    return md_mtl_fail("%s: stream '%s' is inside a render pass; call md_gpu_render_end first", what, s->label);
}

static id<MTLComputeCommandEncoder> md_mtl_compute_encoder(md_gpu_stream_t s) {
    if (!md_mtl_not_in_pass(s, "operation")) return nil;
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
    if (!md_mtl_not_in_pass(s, "operation")) return nil;
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
    if (s && md_mtl_check_no_upload(s, "md_gpu_barrier")) md_mtl_not_in_pass(s, "md_gpu_barrier");
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
    for (uint32_t k = 0; k < MD_MTL_TEMP_KINDS; ++k) {
        md_mtl_vec_init(&s->temp[k].active, sizeof(md_mtl_chunk_t*));
        md_mtl_vec_init(&s->temp[k].spare,  sizeof(md_mtl_chunk_t*));
    }
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
    if (kind != MD_GPU_STREAM_COMPUTE && kind != MD_GPU_STREAM_TRANSFER && kind != MD_GPU_STREAM_GRAPHICS) {
        md_mtl_fail("md_gpu_stream_create: invalid stream kind %d", (int)kind);
        return NULL;
    }
    /* Metal queues are universal: the kind only decides what is validated. */
    return md_mtl_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    switch (kind) {
    case MD_GPU_STREAM_TRANSFER: return dev->default_transfer;
    case MD_GPU_STREAM_GRAPHICS: return dev->default_graphics;
    default:                     return dev->default_compute;
    }
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
    if (!md_mtl_not_in_pass(s, "md_gpu_stream_record")) return out;
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
    if (!md_mtl_not_in_pass(s, "md_gpu_stream_wait")) return;
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
    if (s && md_mtl_not_in_pass(s, "md_gpu_stream_flush")) md_mtl_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    if (!md_mtl_not_in_pass(s, "md_gpu_stream_sync")) return;
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
    md_mtl_temp_free_all(s);
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
    md_mtl_abandon_pass(s);
    md_gpu_stream_sync(s);

    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->streams, s);
    md_mtl_forget_stream_locked(dev, s);
    md_mutex_unlock(&dev->device_mutex);

    md_mtl_stream_free(s);
}

/* =========================================================================
   8. Memory
   =========================================================================

   The same design as md_gpu_vulkan.c section 9: one TLSF heap per memory
   kind over large MTLBuffer chunks, stream-ordered frees parked on
   pending_frees until their stream passes them, DEVICE memory reusable at
   once by the freeing stream, and empty chunks kept up to heap_cache_limit.
   Chunks are made resident like any other buffer. Metal 3 orders everything
   within a stream, so same-stream reuse needs no extra barrier here. */

static uint64_t md_mtl_heap_chunk_size(const md_mtl_heap_t* h, uint64_t need) {
    uint64_t size = MD_MTL_HEAP_CHUNK_MIN;
    for (size_t i = 0; i < h->chunks.count && size < MD_MTL_HEAP_CHUNK_MAX; ++i) size *= 2;
    if (need > size) size = md_mtl_align_up(need, MD_MTL_HEAP_LARGE_ALIGN);
    return size;
}

/* Caller holds device_mutex (residency). */
static md_mtl_chunk_t* md_mtl_chunk_create_locked(md_gpu_device_t dev, md_gpu_mem_kind_t kind, uint64_t size) {
    md_mtl_chunk_t* c = (md_mtl_chunk_t*)md_alloc(dev->alloc, sizeof(md_mtl_chunk_t));
    if (!c) { md_mtl_fail("out of memory"); return NULL; }
    memset(c, 0, sizeof(*c));
    void* host = NULL;
    id<MTLBuffer> buf = nil;
    @autoreleasepool {
        if (!md_mtl_create_raw_buffer(dev, size, kind, &buf, &c->address, &host)) {
            md_free(dev->alloc, c, sizeof(*c));
            return NULL;
        }
    }
    c->buffer = buf;
    c->host   = (uint8_t*)host;
    c->size   = size;
    c->kind   = kind;
    md_mtl_make_resident_locked(dev, buf);
    return c;
}

/* Caller holds device_mutex. */
static void md_mtl_chunk_destroy_locked(md_gpu_device_t dev, md_mtl_chunk_t* c) {
    md_mtl_destroy_raw_buffer_locked(dev, c->buffer);
    md_free(dev->alloc, c, sizeof(*c));
}

/* Caller holds device_mutex. */
static void md_mtl_heap_release_node_locked(md_gpu_device_t dev, md_gpu_mem_kind_t kind, md_tlsf_node_t* node) {
    md_mtl_heap_t* h = &dev->heaps[kind];
    md_tlsf_node_t* f = md_tlsf_free(&h->tlsf, node);
    if (!md_tlsf_region_empty(f)) return;
    md_mtl_chunk_t* c = (md_mtl_chunk_t*)f->region;
    c->empty      = true;
    c->empty_node = f;
    h->empty_bytes += c->size;
    if (h->empty_bytes > dev->heap_cache_limit) {
        md_tlsf_remove_region(&h->tlsf, f);
        md_mtl_vec_remove_ptr(&h->chunks, c);
        h->empty_bytes -= c->size;
        h->reserved    -= c->size;
        md_mtl_chunk_destroy_locked(dev, c);
    }
}

/* Caller holds device_mutex. */
static void md_mtl_process_pending_frees_locked(md_gpu_device_t dev) {
    for (size_t i = 0; i < dev->pending_frees.count;) {
        md_mtl_pending_free_t pf = MD_MTL_VEC_AT(dev->pending_frees, md_mtl_pending_free_t, i);
        if (md_mtl_stream_completed(pf.stream) >= pf.value) {
            md_mtl_vec_remove(&dev->pending_frees, i);
            md_mtl_heap_release_node_locked(dev, pf.kind, pf.node);
        } else {
            ++i;
        }
    }
}

/* Caller holds device_mutex. */
static md_tlsf_node_t* md_mtl_heap_alloc_locked(md_gpu_device_t dev, md_gpu_stream_t s,
                                                md_gpu_mem_kind_t kind, uint64_t need) {
    md_mtl_heap_t* h = &dev->heaps[kind];
    md_tlsf_node_t* node = NULL;
    if (kind == MD_GPU_MEM_DEVICE) {
        for (size_t i = 0; i < dev->pending_frees.count; ++i) {
            md_mtl_pending_free_t* pf = &MD_MTL_VEC_AT(dev->pending_frees, md_mtl_pending_free_t, i);
            if (pf->stream != s || pf->kind != kind) continue;
            if (pf->node->size < need || pf->node->size > need + need / 4 + MD_MTL_HEAP_ALIGN) continue;
            node = pf->node;
            md_mtl_vec_remove(&dev->pending_frees, i);
            break;
        }
    }
    if (!node) node = md_tlsf_alloc(&h->tlsf, need);
    if (!node) {
        const uint64_t size = md_mtl_heap_chunk_size(h, need);
        md_mtl_chunk_t* c = md_mtl_chunk_create_locked(dev, kind, size);
        if (!c) return NULL;
        md_mtl_chunk_t** slot = (md_mtl_chunk_t**)md_mtl_vec_push(&h->chunks, dev->alloc);
        if (!slot || !md_tlsf_add_region(&h->tlsf, c, size)) {
            if (slot) h->chunks.count--;
            md_mtl_chunk_destroy_locked(dev, c);
            md_mtl_fail("out of memory");
            return NULL;
        }
        *slot = c;
        h->reserved += size;
        node = md_tlsf_alloc(&h->tlsf, need);
        if (!node) { md_mtl_fail("md_gpu_malloc: internal error, a fresh chunk did not fit"); return NULL; }
    }
    md_mtl_chunk_t* c = (md_mtl_chunk_t*)node->region;
    if (c->empty) {
        c->empty      = false;
        c->empty_node = NULL;
        h->empty_bytes -= c->size;
    }
    return node;
}

md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t stream, md_gpu_mem_kind_t kind, size_t size) {
    md_gpu_mem_t out = {0, NULL};
    if (!stream) { md_mtl_fail("md_gpu_malloc: null stream"); return out; }
    if ((unsigned)kind >= (unsigned)MD_GPU_MEM_KIND_COUNT) { md_mtl_fail("md_gpu_malloc: invalid memory kind %d", (int)kind); return out; }
    if (size == 0) return out;
    md_gpu_device_t dev = stream->device;
    const uint64_t need = md_mtl_align_up(size, MD_MTL_HEAP_ALIGN);

    md_mtl_range_t* r = (md_mtl_range_t*)md_alloc(dev->alloc, sizeof(md_mtl_range_t));
    if (!r) { md_mtl_fail("out of memory"); return out; }

    md_mutex_lock(&dev->device_mutex);
    md_mtl_process_pending_frees_locked(dev);
    md_tlsf_node_t* node = md_mtl_heap_alloc_locked(dev, stream, kind, need);
    if (!node) {
        md_mutex_unlock(&dev->device_mutex);
        md_free(dev->alloc, r, sizeof(*r));
        return out;
    }
    md_mtl_chunk_t* c = (md_mtl_chunk_t*)node->region;
    r->address      = c->address + node->offset;
    r->size         = size;
    r->chunk        = c;
    r->chunk_offset = node->offset;
    r->node         = node;
    r->kind         = kind;
    if (!md_mtl_registry_insert_locked(dev, r)) {
        md_mtl_heap_release_node_locked(dev, kind, node);
        md_mutex_unlock(&dev->device_mutex);
        md_free(dev->alloc, r, sizeof(*r));
        md_mtl_fail("out of memory");
        return out;
    }
    md_mtl_heap_t* h = &dev->heaps[kind];
    h->in_use += node->size;
    if (h->in_use > h->peak_in_use) h->peak_in_use = h->in_use;
    h->allocations++;
    md_mutex_unlock(&dev->device_mutex);

    out.gpu = r->address;
    out.cpu = c->host ? c->host + node->offset : NULL;
    return out;
}

void md_gpu_free(md_gpu_stream_t stream, md_gpu_addr_t addr) {
    if (!addr) return;
    if (!stream) { md_mtl_fail("md_gpu_free: a stream is required"); return; }
    md_gpu_device_t dev = stream->device;

    md_mutex_lock(&dev->device_mutex);
    md_mtl_range_t* r = md_mtl_registry_find_locked(dev, addr);
    if (!r || !r->node || r->address != addr) {
        md_mutex_unlock(&dev->device_mutex);
        md_mtl_fail("md_gpu_free: 0x%llx is not the start of a live md_gpu_malloc allocation", (unsigned long long)addr);
        return;
    }
    md_mtl_registry_remove_locked(dev, r);
    md_mtl_heap_t* h = &dev->heaps[r->kind];
    h->in_use -= r->node->size;
    h->allocations--;

    const uint64_t at = md_mtl_stream_position(stream);
    bool queued = false;
    if (at != 0 && md_mtl_stream_completed(stream) < at) {
        md_mtl_pending_free_t* pf = (md_mtl_pending_free_t*)md_mtl_vec_push(&dev->pending_frees, dev->alloc);
        if (pf) {
            pf->node   = r->node;
            pf->kind   = r->kind;
            pf->stream = stream;
            pf->value  = at;
        } else {
            md_mtl_fail("out of memory recording a free; %llu bytes leaked", (unsigned long long)r->node->size);
        }
        queued = true;
    }
    if (!queued) md_mtl_heap_release_node_locked(dev, r->kind, r->node);
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, r, sizeof(*r));
}

/* ---- Temp arenas -------------------------------------------------------------
   As in md_gpu_vulkan.c: per stream and kind, a stack of chunks owned by the
   scope depth that acquired them; retired chunks are restarted once the GPU
   has passed them, or continued past their cursor before that. */

static int md_mtl_temp_index(md_gpu_mem_kind_t kind) {
    return kind == MD_GPU_MEM_DEVICE ? 0 : (kind == MD_GPU_MEM_HOST_WRITE ? 1 : -1);
}

/* Caller holds device_mutex. */
static void md_mtl_temp_chunk_destroy_locked(md_gpu_device_t dev, md_mtl_chunk_t* c) {
    if (c->range) {
        md_mtl_registry_remove_locked(dev, c->range);
        md_free(dev->alloc, c->range, sizeof(md_mtl_range_t));
    }
    dev->heaps[c->kind].temp_bytes -= c->size;
    md_mtl_chunk_destroy_locked(dev, c);
}

static md_mtl_chunk_t* md_mtl_temp_acquire(md_gpu_stream_t s, md_gpu_mem_kind_t kind, uint64_t need) {
    md_gpu_device_t dev = s->device;
    md_mtl_temp_arena_t* a = &s->temp[md_mtl_temp_index(kind)];
    const uint64_t done = md_mtl_stream_completed(s);

    md_mtl_chunk_t* c = NULL;
    for (size_t i = 0; i < a->spare.count && !c; ++i) {
        md_mtl_chunk_t* sp = MD_MTL_VEC_AT(a->spare, md_mtl_chunk_t*, i);
        if (sp->retire_value <= done && sp->size >= need) {
            c = sp;
            c->cursor = 0;
            md_mtl_vec_remove(&a->spare, i);
        }
    }
    for (size_t i = a->spare.count; i-- > 0 && !c;) {
        md_mtl_chunk_t* sp = MD_MTL_VEC_AT(a->spare, md_mtl_chunk_t*, i);
        if (sp->cursor + need <= sp->size) {
            c = sp;
            md_mtl_vec_remove(&a->spare, i);
        }
    }
    if (!c) {
        const uint64_t size = need > MD_MTL_TEMP_CHUNK_MIN ? md_mtl_align_up(need, MD_MTL_HEAP_LARGE_ALIGN) : MD_MTL_TEMP_CHUNK_MIN;
        md_mtl_range_t* r = (md_mtl_range_t*)md_alloc(dev->alloc, sizeof(md_mtl_range_t));
        if (!r) { md_mtl_fail("out of memory"); return NULL; }
        md_mutex_lock(&dev->device_mutex);
        c = md_mtl_chunk_create_locked(dev, kind, size);
        bool ok = c != NULL;
        if (ok) {
            r->address      = c->address;
            r->size         = c->size;
            r->chunk        = c;
            r->chunk_offset = 0;
            r->node         = NULL;
            r->kind         = kind;
            ok = md_mtl_registry_insert_locked(dev, r);
            if (ok) { c->range = r; dev->heaps[kind].temp_bytes += c->size; }
            else    { md_mtl_chunk_destroy_locked(dev, c); md_mtl_fail("out of memory"); }
        }
        md_mutex_unlock(&dev->device_mutex);
        if (!ok) { md_free(dev->alloc, r, sizeof(*r)); return NULL; }
    }
    md_mtl_chunk_t** slot = (md_mtl_chunk_t**)md_mtl_vec_push(&a->active, dev->alloc);
    if (!slot) {
        md_mutex_lock(&dev->device_mutex);
        md_mtl_temp_chunk_destroy_locked(dev, c);
        md_mutex_unlock(&dev->device_mutex);
        md_mtl_fail("out of memory");
        return NULL;
    }
    *slot = c;
    c->retire_value = 0;
    c->owner_depth  = s->temp_depth;
    return c;
}

md_gpu_temp_t md_gpu_temp_begin(md_gpu_stream_t stream) {
    md_gpu_temp_t t = {NULL, 0};
    if (!stream) { md_mtl_fail("md_gpu_temp_begin: null stream"); return t; }
    t.stream = stream;
    t.depth  = ++stream->temp_depth;
    return t;
}

md_gpu_mem_t md_gpu_temp_alloc(md_gpu_stream_t s, md_gpu_mem_kind_t kind, size_t size) {
    md_gpu_mem_t out = {0, NULL};
    if (!s) { md_mtl_fail("md_gpu_temp_alloc: null stream"); return out; }
    if (s->temp_depth == 0) { md_mtl_fail("md_gpu_temp_alloc: stream '%s' has no open temp scope", s->label); return out; }
    const int ti = md_mtl_temp_index(kind);
    if (ti < 0) { md_mtl_fail("md_gpu_temp_alloc: only MD_GPU_MEM_DEVICE and MD_GPU_MEM_HOST_WRITE memory is temporary"); return out; }
    if (size == 0) return out;
    const uint64_t need = md_mtl_align_up(size, MD_MTL_HEAP_ALIGN);

    md_mtl_temp_arena_t* a = &s->temp[ti];
    md_mtl_chunk_t* c = a->active.count ? MD_MTL_VEC_AT(a->active, md_mtl_chunk_t*, a->active.count - 1) : NULL;
    if (!c || c->cursor + need > c->size) {
        c = md_mtl_temp_acquire(s, kind, need);
        if (!c) return out;
    }
    out.gpu = c->address + c->cursor;
    out.cpu = c->host ? c->host + c->cursor : NULL;
    c->cursor += need;
    return out;
}

void md_gpu_temp_end(md_gpu_stream_t s, md_gpu_temp_t scope) {
    if (!s) { md_mtl_fail("md_gpu_temp_end: null stream"); return; }
    if (scope.stream != s) { md_mtl_fail("md_gpu_temp_end: the scope was begun on another stream"); return; }
    if (scope.depth == 0 || scope.depth != s->temp_depth) {
        md_mtl_fail("md_gpu_temp_end: scopes must end in reverse order (ending depth %u, innermost open is %u)",
                    scope.depth, s->temp_depth);
        return;
    }
    const uint64_t at = md_mtl_stream_position(s);
    for (uint32_t k = 0; k < MD_MTL_TEMP_KINDS; ++k) {
        md_mtl_temp_arena_t* a = &s->temp[k];
        while (a->active.count > 0) {
            md_mtl_chunk_t* c = MD_MTL_VEC_AT(a->active, md_mtl_chunk_t*, a->active.count - 1);
            if (c->owner_depth < scope.depth) break;
            a->active.count--;
            c->retire_value = at;
            md_mtl_chunk_t** slot = (md_mtl_chunk_t**)md_mtl_vec_push(&a->spare, s->device->alloc);
            if (slot) *slot = c;
            else md_mtl_fail("out of memory retiring a temp chunk; %llu bytes leaked", (unsigned long long)c->size);
        }
    }
    s->temp_depth--;
}

/* Every chunk of a stream's temp arenas. The stream must be idle. */
static void md_mtl_temp_free_all(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    md_mutex_lock(&dev->device_mutex);
    for (uint32_t k = 0; k < MD_MTL_TEMP_KINDS; ++k) {
        md_mtl_temp_arena_t* a = &s->temp[k];
        for (size_t i = 0; i < a->active.count; ++i) md_mtl_temp_chunk_destroy_locked(dev, MD_MTL_VEC_AT(a->active, md_mtl_chunk_t*, i));
        for (size_t i = 0; i < a->spare.count;  ++i) md_mtl_temp_chunk_destroy_locked(dev, MD_MTL_VEC_AT(a->spare,  md_mtl_chunk_t*, i));
        md_mtl_vec_free(&a->active, dev->alloc);
        md_mtl_vec_free(&a->spare,  dev->alloc);
    }
    md_mutex_unlock(&dev->device_mutex);
    s->temp_depth = 0;
}

bool md_gpu_memory_stats(md_gpu_device_t dev, md_gpu_mem_kind_t kind, md_gpu_memory_stats_t* out) {
    if (!dev || !out || (unsigned)kind >= (unsigned)MD_GPU_MEM_KIND_COUNT) return md_mtl_fail("md_gpu_memory_stats: invalid argument");
    memset(out, 0, sizeof(*out));
    md_mutex_lock(&dev->device_mutex);
    const md_mtl_heap_t* h = &dev->heaps[kind];
    out->bytes_in_use      = h->in_use;
    out->bytes_peak_in_use = h->peak_in_use;
    out->bytes_reserved    = h->reserved;
    out->bytes_temp        = h->temp_bytes;
    out->bytes_textures    = kind == MD_GPU_MEM_DEVICE ? dev->texture_bytes : 0;
    out->allocations       = h->allocations;
    out->chunks            = (uint32_t)h->chunks.count;
    md_mutex_unlock(&dev->device_mutex);
    return true;
}

/* At device destruction. Caller holds device_mutex; streams are idle and gone. */
static void md_mtl_heaps_free_locked(md_gpu_device_t dev) {
    for (size_t i = 0; i < dev->registry.count; ++i) {
        md_free(dev->alloc, MD_MTL_VEC_AT(dev->registry, md_mtl_range_t*, i), sizeof(md_mtl_range_t));
    }
    dev->registry.count = 0;
    dev->pending_frees.count = 0;
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_mtl_heap_t* h = &dev->heaps[k];
        for (size_t i = 0; i < h->chunks.count; ++i) md_mtl_chunk_destroy_locked(dev, MD_MTL_VEC_AT(h->chunks, md_mtl_chunk_t*, i));
        md_mtl_vec_free(&h->chunks, dev->alloc);
        md_tlsf_destroy(&h->tlsf);
    }
    md_mtl_vec_free(&dev->pending_frees, dev->alloc);
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
    md_mtl_span_t d, sp;
    if (!md_mtl_resolve(s->device, dst, size, &d,  "md_gpu_copy (dst)")) return false;
    if (!md_mtl_resolve(s->device, src, size, &sp, "md_gpu_copy (src)")) return false;
    return md_mtl_record_buffer_copy(s, sp.buffer, sp.offset, src, d.buffer, d.offset, dst, size);
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
    if (!md_mtl_not_in_pass(s, "md_gpu_upload")) return false;
    md_mtl_span_t d;
    if (!md_mtl_resolve(s->device, dst, size, &d, "md_gpu_upload")) return false;
    if (d.host && md_mtl_stream_idle(s)) {
        memcpy(d.host, src, size);
        return true;
    }
    uint64_t addr, off; void* host; id<MTLBuffer> buf = nil;
    if (!md_mtl_arena_alloc(s, size, MD_MTL_ARG_ALIGN, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    return md_mtl_record_buffer_copy(s, buf, off, addr, d.buffer, d.offset, dst, size);
}

bool md_gpu_memset(md_gpu_stream_t s, md_gpu_addr_t dst, uint8_t value, size_t size) {
    if (!s) return md_mtl_fail("md_gpu_memset: null stream");
    if (size == 0) return true;
    if (!md_mtl_check_no_upload(s, "md_gpu_memset")) return false;
    if (!md_mtl_not_in_pass(s, "md_gpu_memset")) return false;
    md_mtl_span_t b;
    if (!md_mtl_resolve(s->device, dst, size, &b, "md_gpu_memset")) return false;
    const uint64_t off = b.offset;
    /* fillBuffer's range must be 4-byte aligned on macOS, as for copies: fill
       the aligned middle, and the (at most 3-byte) head and tail by kernel. */
    const uint64_t begin = off, end = off + size;
    const uint64_t abeg  = md_mtl_align_up(begin, 4);
    const uint64_t aend  = end & ~3ull;
    if (aend > abeg) {
        id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
        if (!enc) return false;
        [enc fillBuffer:b.buffer range:NSMakeRange((NSUInteger)abeg, (NSUInteger)(aend - abeg)) value:value];
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
    if (!md_mtl_not_in_pass(s, "md_gpu_upload_begin")) return NULL;
    md_mtl_span_t b;
    if (!md_mtl_resolve(s->device, dst, size, &b, "md_gpu_upload_begin")) return NULL;
    if (b.host && md_mtl_stream_idle(s)) {
        s->upload_open = true; s->upload_direct = true;
        s->upload_dst = dst; s->upload_size = size;
        return b.host;
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

    md_mtl_span_t d;
    if (!md_mtl_resolve(s->device, s->upload_dst, s->upload_size, &d, "md_gpu_upload_end")) return false;
    for (size_t i = 0; i < s->arena.pages.count; ++i) {
        md_mtl_page_t* p = MD_MTL_VEC_AT(s->arena.pages, md_mtl_page_t*, i);
        if (s->upload_src_addr >= p->address && s->upload_src_addr < p->address + p->capacity) {
            return md_mtl_record_buffer_copy(s, p->buffer, s->upload_src_addr - p->address, s->upload_src_addr,
                                             d.buffer, d.offset, s->upload_dst, s->upload_size);
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
        if (t->resident) md_mtl_end_residency_locked(dev, t->texture);
        MD_MTL_RELEASE(t->texture);
    }
    md_free(dev->alloc, t, sizeof(*t));
}

md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t s, const md_gpu_texture_desc_t* desc) {
    if (!s || !desc) { md_mtl_fail("md_gpu_texture_create: null argument"); return NULL; }
    md_gpu_device_t dev = s->device;
    if (!md_mtl_check_no_upload(s, "md_gpu_texture_create")) return NULL;
    if (!md_mtl_not_in_pass(s, "md_gpu_texture_create")) return NULL;

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
    if ((d.usage & MD_GPU_TEX_RENDER_TARGET) && d.type == MD_GPU_TEX_3D) {
        md_mtl_fail("texture '%s': RENDER_TARGET needs a 2D or 2D_ARRAY texture", label);
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
        md_gpu_texture_t* slot = (md_gpu_texture_t*)md_mtl_vec_push(&dev->textures, dev->alloc);
        if (!slot) ok = md_mtl_fail("out of memory");
        else {
            *slot = t;
            md_mtl_make_resident_locked(dev, t->texture);
            t->resident = true;
            dev->texture_bytes += t->bytes;
        }
    }
    if (!ok) {
        /* Never made resident, so texture_free ends no residency. */
        md_mtl_texture_free(dev, t);
        t = NULL;
    }
    md_mutex_unlock(&dev->device_mutex);
    (void)s;   /* Metal has no layout transition to order; creation is complete now. */
    return t;
}

void md_gpu_texture_destroy(md_gpu_texture_t t) {
    if (!t) return;
    if (t->external) { md_mtl_fail("md_gpu_texture_destroy: '%s' belongs to a surface and is not destroyed by the caller", t->label); return; }
    md_gpu_device_t dev = t->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->textures, t);
    dev->texture_bytes -= t->bytes;
    md_mtl_retire_locked(dev, MD_MTL_RETIRE_TEXTURE, t);
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
    md_mtl_span_t b;
    if (!md_mtl_resolve(s->device, src, cr.bytes, &b, "md_gpu_copy_to_texture")) return false;
    if (!md_mtl_check_buffer_offset(t, b.offset, "md_gpu_copy_to_texture")) return false;
    return md_mtl_record_texture_copy(s, t, b.buffer, b.offset, &cr, true);
}

bool md_gpu_copy_from_texture(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!s || !t) return md_mtl_fail("md_gpu_copy_from_texture: null argument");
    md_mtl_copy_region_t cr;
    if (!md_mtl_resolve_region(t, region, &cr, "md_gpu_copy_from_texture")) return false;
    md_mtl_span_t b;
    if (!md_mtl_resolve(s->device, dst, cr.bytes, &b, "md_gpu_copy_from_texture")) return false;
    if (!md_mtl_check_buffer_offset(t, b.offset, "md_gpu_copy_from_texture")) return false;
    return md_mtl_record_texture_copy(s, t, b.buffer, b.offset, &cr, false);
}

bool md_gpu_upload_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, const void* src, size_t size) {
    if (!s || !t || !src) return md_mtl_fail("md_gpu_upload_texture: null argument");
    if (!md_mtl_check_no_upload(s, "md_gpu_upload_texture")) return false;
    if (!md_mtl_not_in_pass(s, "md_gpu_upload_texture")) return false;
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
    md_mtl_span_t b;
    if (!md_mtl_resolve(s->device, grid, 3 * sizeof(uint32_t), &b, "md_gpu_launch_indirect")) return false;
    if (b.offset % 4 != 0) return md_mtl_fail("md_gpu_launch_indirect: grid address must be 4-byte aligned");
    return md_mtl_launch_common(s, k, md_gpu_grid(1, 1, 1), args, args_size, b.buffer, b.offset, true, false);
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
   11. Rendering
   =========================================================================

   A pass is one MTLRenderCommandEncoder, opened by md_gpu_render_begin and
   fenced like every other encoder. Pipelines bake shaders, formats and
   blend; depth test and write live in MTLDepthStencilState objects, cached
   on the device per (compare, write); the rest is encoder state. Metal's
   clip space is already +Y up with a top-left framebuffer origin, so the
   viewport is passed through unflipped. */

static bool md_mtl_format_is_uint(md_gpu_format_t f) {
    return f == MD_GPU_FORMAT_R32_UINT || f == MD_GPU_FORMAT_RG32_UINT || f == MD_GPU_FORMAT_RGBA32_UINT;
}

static MTLBlendFactor md_mtl_blend_factor(md_gpu_blend_factor_t f, bool* ok) {
    switch (f) {
    case MD_GPU_BLEND_ZERO:                return MTLBlendFactorZero;
    case MD_GPU_BLEND_ONE:                 return MTLBlendFactorOne;
    case MD_GPU_BLEND_SRC_COLOR:           return MTLBlendFactorSourceColor;
    case MD_GPU_BLEND_ONE_MINUS_SRC_COLOR: return MTLBlendFactorOneMinusSourceColor;
    case MD_GPU_BLEND_SRC_ALPHA:           return MTLBlendFactorSourceAlpha;
    case MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA: return MTLBlendFactorOneMinusSourceAlpha;
    case MD_GPU_BLEND_DST_COLOR:           return MTLBlendFactorDestinationColor;
    case MD_GPU_BLEND_ONE_MINUS_DST_COLOR: return MTLBlendFactorOneMinusDestinationColor;
    case MD_GPU_BLEND_DST_ALPHA:           return MTLBlendFactorDestinationAlpha;
    case MD_GPU_BLEND_ONE_MINUS_DST_ALPHA: return MTLBlendFactorOneMinusDestinationAlpha;
    case MD_GPU_BLEND_CONSTANT:            return MTLBlendFactorBlendColor;
    case MD_GPU_BLEND_ONE_MINUS_CONSTANT:  return MTLBlendFactorOneMinusBlendColor;
    case MD_GPU_BLEND_SRC_ALPHA_SATURATE:  return MTLBlendFactorSourceAlphaSaturated;
    default: *ok = false;                  return MTLBlendFactorZero;
    }
}

static MTLBlendOperation md_mtl_blend_op(md_gpu_blend_op_t op, bool* ok) {
    switch (op) {
    case MD_GPU_BLEND_OP_ADD:              return MTLBlendOperationAdd;
    case MD_GPU_BLEND_OP_SUBTRACT:         return MTLBlendOperationSubtract;
    case MD_GPU_BLEND_OP_REVERSE_SUBTRACT: return MTLBlendOperationReverseSubtract;
    case MD_GPU_BLEND_OP_MIN:              return MTLBlendOperationMin;
    case MD_GPU_BLEND_OP_MAX:              return MTLBlendOperationMax;
    default: *ok = false;                  return MTLBlendOperationAdd;
    }
}

/* The entry point of a vertex or fragment library: the name asked for, its
   Slang-mangled form, or the library's only function. +1 or nil. */
static id<MTLFunction> md_mtl_find_function(id<MTLLibrary> lib, const char* entry) {
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:entry]];
    if (!fn) fn = [lib newFunctionWithName:[NSString stringWithFormat:@"%s_0", entry]];
    if (!fn && lib.functionNames.count == 1) fn = [lib newFunctionWithName:lib.functionNames[0]];
    return fn;
}

/* The first buffer a stage binds is its root; with MD_SHADER_ARGS it is the
   only one. UINT32_MAX when the stage binds none (no argument struct). */
static uint32_t md_mtl_reflect_root(NSArray<id<MTLBinding>>* bindings) {
    for (id<MTLBinding> b in bindings) {
        if (b.type == MTLBindingTypeBuffer) return (uint32_t)b.index;
    }
    return UINT32_MAX;
}

static void md_mtl_pipeline_free(md_gpu_pipeline_t p) {
    MD_MTL_RELEASE(p->pso);
    md_free(p->device->alloc, p, sizeof(*p));
}

md_gpu_pipeline_t md_gpu_pipeline_create(md_gpu_device_t dev, const md_gpu_pipeline_desc_t* desc) {
    if (!dev || !desc) { md_mtl_fail("md_gpu_pipeline_create: null argument"); return NULL; }
    const char* label = desc->label ? desc->label : "pipeline";
    if (!desc->vertex.code || desc->vertex.code_size == 0) { md_mtl_fail("pipeline '%s': no vertex shader", label); return NULL; }
    const bool has_fs = desc->fragment.code != NULL && desc->fragment.code_size > 0;
    if (desc->color_count > MD_GPU_MAX_COLOR_TARGETS) {
        md_mtl_fail("pipeline '%s': %u colour targets, at most %u", label, desc->color_count, MD_GPU_MAX_COLOR_TARGETS);
        return NULL;
    }
    if (!has_fs && desc->color_count > 0) { md_mtl_fail("pipeline '%s': colour targets without a fragment shader", label); return NULL; }
    MTLPrimitiveType prim;
    MTLPrimitiveTopologyClass topo_class;
    switch (desc->topology) {
    case MD_GPU_TOPOLOGY_TRIANGLES:      prim = MTLPrimitiveTypeTriangle;      topo_class = MTLPrimitiveTopologyClassTriangle; break;
    case MD_GPU_TOPOLOGY_TRIANGLE_STRIP: prim = MTLPrimitiveTypeTriangleStrip; topo_class = MTLPrimitiveTopologyClassTriangle; break;
    case MD_GPU_TOPOLOGY_LINES:          prim = MTLPrimitiveTypeLine;          topo_class = MTLPrimitiveTopologyClassLine;     break;
    case MD_GPU_TOPOLOGY_LINE_STRIP:     prim = MTLPrimitiveTypeLineStrip;     topo_class = MTLPrimitiveTopologyClassLine;     break;
    case MD_GPU_TOPOLOGY_POINTS:         prim = MTLPrimitiveTypePoint;         topo_class = MTLPrimitiveTopologyClassPoint;    break;
    default: md_mtl_fail("pipeline '%s': invalid topology %d", label, (int)desc->topology); return NULL;
    }
    const uint32_t vs_args = desc->vertex.args_size;
    const uint32_t fs_args = has_fs ? desc->fragment.args_size : 0;
    if (vs_args && fs_args && vs_args != fs_args) {
        md_mtl_fail("pipeline '%s': the vertex shader takes a %u-byte argument struct and the fragment shader %u bytes; "
                    "both stages read the same one", label, vs_args, fs_args);
        return NULL;
    }
    for (uint32_t i = 0; i < desc->color_count; ++i) {
        const md_mtl_fmt_t fi = md_mtl_format_info(desc->color[i].format);
        if (fi.fmt == MTLPixelFormatInvalid || fi.depth) {
            md_mtl_fail("pipeline '%s': colour target %u has %s, which is not a colour format", label, i,
                        fi.fmt == MTLPixelFormatInvalid ? "no format" : fi.name);
            return NULL;
        }
        if (desc->color[i].blend.enable && md_mtl_format_is_uint(desc->color[i].format)) {
            md_mtl_fail("pipeline '%s': colour target %u is %s, an integer format, and cannot blend", label, i, fi.name);
            return NULL;
        }
    }
    md_mtl_fmt_t dfi = md_mtl_format_info(desc->depth_format);
    if (desc->depth_format != MD_GPU_FORMAT_INVALID && !dfi.depth) {
        md_mtl_fail("pipeline '%s': depth_format %s is not a depth format", label,
                    dfi.fmt == MTLPixelFormatInvalid ? "(invalid)" : dfi.name);
        return NULL;
    }

    md_gpu_pipeline_t p = (md_gpu_pipeline_t)md_alloc(dev->alloc, sizeof(md_gpu_pipeline));
    if (!p) { md_mtl_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device      = dev;
    p->primitive   = prim;
    p->color_count = desc->color_count;
    for (uint32_t i = 0; i < desc->color_count; ++i) p->color[i] = desc->color[i].format;
    p->depth       = desc->depth_format;
    p->args_size   = vs_args ? vs_args : fs_args;
    p->vs_root     = UINT32_MAX;
    p->fs_root     = UINT32_MAX;
    snprintf(p->label, sizeof(p->label), "%s", label);

    bool ok = true;
    @autoreleasepool {
        id<MTLLibrary>  vlib = md_mtl_library_from_blob(dev, desc->vertex.code, desc->vertex.code_size, label);
        id<MTLLibrary>  flib = has_fs && vlib ? md_mtl_library_from_blob(dev, desc->fragment.code, desc->fragment.code_size, label) : nil;
        id<MTLFunction> vfn  = vlib ? md_mtl_find_function(vlib, desc->vertex.entry_point ? desc->vertex.entry_point : "main") : nil;
        id<MTLFunction> ffn  = flib ? md_mtl_find_function(flib, desc->fragment.entry_point ? desc->fragment.entry_point : "main") : nil;
        if (!vlib || (has_fs && !flib)) {
            ok = false;   /* reported by md_mtl_library_from_blob */
        } else if (!vfn || (has_fs && !ffn)) {
            ok = md_mtl_fail("pipeline '%s': %s entry point not found", label, vfn ? "fragment" : "vertex");
        }

        if (ok) {
            MTLRenderPipelineDescriptor* rpd = [[MTLRenderPipelineDescriptor alloc] init];
            rpd.label                  = [NSString stringWithUTF8String:p->label];
            rpd.vertexFunction         = vfn;
            rpd.fragmentFunction       = ffn;
            rpd.inputPrimitiveTopology = topo_class;
            rpd.rasterSampleCount      = 1;
            for (uint32_t i = 0; i < desc->color_count && ok; ++i) {
                const md_gpu_color_target_t* ct = &desc->color[i];
                MTLRenderPipelineColorAttachmentDescriptor* ca = rpd.colorAttachments[i];
                ca.pixelFormat = md_mtl_format_info(ct->format).fmt;
                MTLColorWriteMask mask = MTLColorWriteMaskNone;
                if (!(ct->write_disable & MD_GPU_COLOR_R)) mask |= MTLColorWriteMaskRed;
                if (!(ct->write_disable & MD_GPU_COLOR_G)) mask |= MTLColorWriteMaskGreen;
                if (!(ct->write_disable & MD_GPU_COLOR_B)) mask |= MTLColorWriteMaskBlue;
                if (!(ct->write_disable & MD_GPU_COLOR_A)) mask |= MTLColorWriteMaskAlpha;
                ca.writeMask = mask;
                if (ct->blend.enable) {
                    bool valid = true;
                    ca.blendingEnabled             = YES;
                    ca.sourceRGBBlendFactor        = md_mtl_blend_factor(ct->blend.src_color, &valid);
                    ca.destinationRGBBlendFactor   = md_mtl_blend_factor(ct->blend.dst_color, &valid);
                    ca.rgbBlendOperation           = md_mtl_blend_op(ct->blend.color_op, &valid);
                    ca.sourceAlphaBlendFactor      = md_mtl_blend_factor(ct->blend.src_alpha, &valid);
                    ca.destinationAlphaBlendFactor = md_mtl_blend_factor(ct->blend.dst_alpha, &valid);
                    ca.alphaBlendOperation         = md_mtl_blend_op(ct->blend.alpha_op, &valid);
                    if (!valid) ok = md_mtl_fail("pipeline '%s': invalid blend factor or op on colour target %u", label, i);
                }
            }
            if (desc->depth_format != MD_GPU_FORMAT_INVALID) rpd.depthAttachmentPixelFormat = dfi.fmt;

            if (ok) {
                NSError* err = nil;
                MTLRenderPipelineReflection* refl = nil;
                id<MTLRenderPipelineState> pso = [dev->device newRenderPipelineStateWithDescriptor:rpd
                                                                                           options:MTLPipelineOptionBindingInfo
                                                                                        reflection:&refl
                                                                                             error:&err];
                if (!pso) {
                    ok = md_mtl_fail("pipeline '%s': newRenderPipelineStateWithDescriptor failed: %s", label,
                                     err ? [[err localizedDescription] UTF8String] : "?");
                } else {
                    MD_MTL_OWN(pso);
                    p->pso = pso;
                    if (refl) {
                        if (@available(macOS 13.0, iOS 16.0, *)) {
                            p->vs_root = md_mtl_reflect_root(refl.vertexBindings);
                            p->fs_root = md_mtl_reflect_root(refl.fragmentBindings);
                        }
                    }
                    /* Without reflection, assume what Slang emits for a file
                       with one MD_SHADER_ARGS root. */
                    if (!refl) {
                        p->vs_root = vs_args ? MD_MTL_ARG_BUFFER_INDEX : UINT32_MAX;
                        p->fs_root = fs_args ? MD_MTL_ARG_BUFFER_INDEX : UINT32_MAX;
                    }
                }
            }
            MD_MTL_DROP_NEW(rpd);
        }
        MD_MTL_DROP_NEW(vfn);
        MD_MTL_DROP_NEW(ffn);
        MD_MTL_DROP_NEW(vlib);
        MD_MTL_DROP_NEW(flib);
    }
    if (!ok) { md_mtl_pipeline_free(p); return NULL; }

    md_mutex_lock(&dev->device_mutex);
    md_gpu_pipeline_t* slot = (md_gpu_pipeline_t*)md_mtl_vec_push(&dev->pipelines, dev->alloc);
    if (slot) *slot = p;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { md_mtl_pipeline_free(p); md_mtl_fail("out of memory"); return NULL; }
    return p;
}

/* Immediate, like kernels: command buffers retain the pipeline states set on
   their encoders. */
void md_gpu_pipeline_destroy(md_gpu_pipeline_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->pipelines, p);
    md_mutex_unlock(&dev->device_mutex);
    md_mtl_pipeline_free(p);
}

/* ---- Render passes ---------------------------------------------------------- */

static bool md_mtl_check_attachment(md_gpu_stream_t s, md_gpu_texture_t t, uint32_t mip, uint32_t layer,
                                    bool depth, uint32_t index, uint32_t* w, uint32_t* h) {
    char what[32];
    if (depth) snprintf(what, sizeof(what), "depth attachment");
    else       snprintf(what, sizeof(what), "colour attachment %u", index);
    if (!t) return md_mtl_fail("md_gpu_render_begin: %s has no texture", what);
    if (t->device != s->device) return md_mtl_fail("md_gpu_render_begin: %s '%s' belongs to another device", what, t->label);
    if (!(t->desc.usage & MD_GPU_TEX_RENDER_TARGET)) {
        return md_mtl_fail("md_gpu_render_begin: %s '%s' lacks MD_GPU_TEX_RENDER_TARGET usage", what, t->label);
    }
    if (t->fi.depth != depth) {
        return md_mtl_fail("md_gpu_render_begin: %s '%s' is %s; %s", what, t->label, t->fi.name,
                           depth ? "a depth attachment needs a depth format" : "depth formats go in md_gpu_render_desc_t.depth");
    }
    if (mip >= t->desc.mip_levels) {
        return md_mtl_fail("md_gpu_render_begin: %s '%s': mip %u out of range (%u levels)", what, t->label, mip, t->desc.mip_levels);
    }
    const uint32_t layers = t->desc.type == MD_GPU_TEX_2D_ARRAY ? t->desc.depth_or_layers : 1u;
    if (layer >= layers) {
        return md_mtl_fail("md_gpu_render_begin: %s '%s': layer %u out of range (%u layers)", what, t->label, layer, layers);
    }
    uint32_t mw = t->desc.width  >> mip; if (!mw) mw = 1;
    uint32_t mh = t->desc.height >> mip; if (!mh) mh = 1;
    if (*w == 0) { *w = mw; *h = mh; }
    else if (*w != mw || *h != mh) {
        return md_mtl_fail("md_gpu_render_begin: %s '%s' is %ux%u at mip %u but the pass is %ux%u; every attachment must be the same size",
                           what, t->label, mw, mh, mip, *w, *h);
    }
    return true;
}

static MTLLoadAction md_mtl_load_action(md_gpu_load_t l) {
    switch (l) {
    case MD_GPU_LOAD_CLEAR:     return MTLLoadActionClear;
    case MD_GPU_LOAD_DONT_CARE: return MTLLoadActionDontCare;
    default:                    return MTLLoadActionLoad;
    }
}

static MTLCompareFunction md_mtl_compare(md_gpu_compare_t c) {
    switch (c) {
    case MD_GPU_COMPARE_NEVER:         return MTLCompareFunctionNever;
    case MD_GPU_COMPARE_LESS:          return MTLCompareFunctionLess;
    case MD_GPU_COMPARE_LESS_EQUAL:    return MTLCompareFunctionLessEqual;
    case MD_GPU_COMPARE_EQUAL:         return MTLCompareFunctionEqual;
    case MD_GPU_COMPARE_NOT_EQUAL:     return MTLCompareFunctionNotEqual;
    case MD_GPU_COMPARE_GREATER_EQUAL: return MTLCompareFunctionGreaterEqual;
    case MD_GPU_COMPARE_GREATER:       return MTLCompareFunctionGreater;
    default:                           return MTLCompareFunctionAlways;
    }
}

/* The cached depth state for (compare, write), or nil on failure. */
static id<MTLDepthStencilState> md_mtl_depth_state(md_gpu_device_t dev, md_gpu_compare_t compare, bool write) {
    const unsigned c = (unsigned)compare < 8u ? (unsigned)compare : 0u;
    const unsigned w = write ? 1u : 0u;
    md_mutex_lock(&dev->device_mutex);
    id<MTLDepthStencilState> ds = dev->depth_states[c][w];
    if (!ds) {
        @autoreleasepool {
            MTLDepthStencilDescriptor* d = [[MTLDepthStencilDescriptor alloc] init];
            d.depthCompareFunction = md_mtl_compare(compare);
            d.depthWriteEnabled    = write ? YES : NO;
            ds = [dev->device newDepthStencilStateWithDescriptor:d];
            MD_MTL_DROP_NEW(d);
            if (ds) {
                MD_MTL_OWN(ds);
                dev->depth_states[c][w] = ds;
            }
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    if (!ds) md_mtl_fail("newDepthStencilStateWithDescriptor failed");
    return ds;
}

static void md_mtl_apply_draw_state(md_gpu_stream_t s, const md_gpu_draw_state_t* n, bool all) {
    id<MTLRenderCommandEncoder> enc = s->render_enc;
    const md_gpu_draw_state_t* o = &s->draw_state;
    if (all || n->depth_compare != o->depth_compare || n->depth_write != o->depth_write) {
        id<MTLDepthStencilState> ds = md_mtl_depth_state(s->device, n->depth_compare, n->depth_write);
        if (ds) [enc setDepthStencilState:ds];
    }
    if (all || n->cull != o->cull) {
        [enc setCullMode:n->cull == MD_GPU_CULL_BACK  ? MTLCullModeBack :
                         n->cull == MD_GPU_CULL_FRONT ? MTLCullModeFront : MTLCullModeNone];
    }
    if (all || n->front_clockwise != o->front_clockwise) {
        [enc setFrontFacingWinding:n->front_clockwise ? MTLWindingClockwise : MTLWindingCounterClockwise];
    }
    if (all || n->depth_bias != o->depth_bias || n->depth_bias_slope != o->depth_bias_slope || n->depth_bias_clamp != o->depth_bias_clamp) {
        [enc setDepthBias:n->depth_bias slopeScale:n->depth_bias_slope clamp:n->depth_bias_clamp];
    }
    if (all || memcmp(n->blend_constant, o->blend_constant, sizeof(n->blend_constant)) != 0) {
        [enc setBlendColorRed:n->blend_constant[0] green:n->blend_constant[1] blue:n->blend_constant[2] alpha:n->blend_constant[3]];
    }
    s->draw_state = *n;
}

static void md_mtl_apply_viewport(md_gpu_stream_t s, const md_gpu_viewport_t* v) {
    MTLViewport vp;
    vp.originX = v->x;
    vp.originY = v->y;
    vp.width   = v->width;
    vp.height  = v->height;
    vp.znear   = v->min_depth;
    vp.zfar    = (v->min_depth == 0.0f && v->max_depth == 0.0f) ? 1.0 : v->max_depth;
    [s->render_enc setViewport:vp];
}

static void md_mtl_apply_scissor(md_gpu_stream_t s, const md_gpu_rect_t* r) {
    /* Metal requires the rectangle inside the render target. */
    MTLScissorRect sc;
    sc.x      = r->x < s->pass_width  ? r->x : s->pass_width;
    sc.y      = r->y < s->pass_height ? r->y : s->pass_height;
    sc.width  = r->width  < s->pass_width  - sc.x ? r->width  : s->pass_width  - sc.x;
    sc.height = r->height < s->pass_height - sc.y ? r->height : s->pass_height - sc.y;
    [s->render_enc setScissorRect:sc];
}

/* Fallback residency for a render encoder, as md_mtl_declare_residency does
   for compute. */
static void md_mtl_declare_residency_render(md_gpu_stream_t s, id<MTLRenderCommandEncoder> enc) {
    md_gpu_device_t dev = s->device;
    if (dev->has_residency_set) return;
    md_mutex_lock(&dev->device_mutex);
    if (s->res_gen != dev->res_gen + 1) {
        if (dev->live_res.count > 0) {
            [enc useResources:(__unsafe_unretained id<MTLResource>*)dev->live_res.data
                        count:dev->live_res.count
                        usage:MTLResourceUsageRead | MTLResourceUsageWrite
                       stages:MTLRenderStageVertex | MTLRenderStageFragment];
        }
        s->res_gen = dev->res_gen + 1;
    }
    md_mutex_unlock(&dev->device_mutex);
}

bool md_gpu_render_begin(md_gpu_stream_t s, const md_gpu_render_desc_t* desc) {
    if (!s || !desc) return md_mtl_fail("md_gpu_render_begin: null argument");
    if (s->kind != MD_GPU_STREAM_GRAPHICS) return md_mtl_fail("md_gpu_render_begin: stream '%s' is not a GRAPHICS stream", s->label);
    if (s->in_pass) return md_mtl_fail("md_gpu_render_begin: stream '%s' is already inside a render pass", s->label);
    if (!md_mtl_check_no_upload(s, "md_gpu_render_begin")) return false;
    if (desc->color_count > MD_GPU_MAX_COLOR_TARGETS) {
        return md_mtl_fail("md_gpu_render_begin: %u colour attachments, at most %u", desc->color_count, MD_GPU_MAX_COLOR_TARGETS);
    }
    if (desc->color_count == 0 && !desc->depth.texture) return md_mtl_fail("md_gpu_render_begin: the pass has no attachments");

    uint32_t w = 0, h = 0;
    for (uint32_t i = 0; i < desc->color_count; ++i) {
        const md_gpu_color_attachment_t* a = &desc->color[i];
        if (!md_mtl_check_attachment(s, a->texture, a->mip, a->layer, false, i, &w, &h)) return false;
    }
    if (desc->depth.texture &&
        !md_mtl_check_attachment(s, desc->depth.texture, desc->depth.mip, desc->depth.layer, true, 0, &w, &h)) return false;

    if (!md_mtl_stream_ensure_cmd(s)) return false;
    md_mtl_close_encoder(s, true);
    @autoreleasepool {
        MTLRenderPassDescriptor* rp = [MTLRenderPassDescriptor renderPassDescriptor];
        for (uint32_t i = 0; i < desc->color_count; ++i) {
            const md_gpu_color_attachment_t* a = &desc->color[i];
            MTLRenderPassColorAttachmentDescriptor* ca = rp.colorAttachments[i];
            ca.texture     = a->texture->texture;
            ca.level       = a->mip;
            ca.slice       = a->layer;
            ca.loadAction  = md_mtl_load_action(a->load);
            ca.storeAction = a->store == MD_GPU_STORE_DISCARD ? MTLStoreActionDontCare : MTLStoreActionStore;
            if (md_mtl_format_is_uint(a->texture->desc.format)) {
                ca.clearColor = MTLClearColorMake(a->clear.u32[0], a->clear.u32[1], a->clear.u32[2], a->clear.u32[3]);
            } else {
                ca.clearColor = MTLClearColorMake(a->clear.f32[0], a->clear.f32[1], a->clear.f32[2], a->clear.f32[3]);
            }
        }
        if (desc->depth.texture) {
            const md_gpu_depth_attachment_t* a = &desc->depth;
            MTLRenderPassDepthAttachmentDescriptor* da = rp.depthAttachment;
            da.texture     = a->texture->texture;
            da.level       = a->mip;
            da.slice       = a->layer;
            da.loadAction  = md_mtl_load_action(a->load);
            da.storeAction = a->store == MD_GPU_STORE_DISCARD ? MTLStoreActionDontCare : MTLStoreActionStore;
            da.clearDepth  = a->clear_depth;
        }
        id<MTLRenderCommandEncoder> enc = [s->cmd renderCommandEncoderWithDescriptor:rp];
        if (!enc) return md_mtl_fail("renderCommandEncoderWithDescriptor failed on stream '%s'", s->label);
        MD_MTL_RETAIN(enc);
        if (desc->label) enc.label = [NSString stringWithUTF8String:desc->label];
        if (s->fence_valid) [enc waitForFence:s->fence beforeStages:MTLRenderStageVertex | MTLRenderStageFragment];
        s->render_enc = enc;
        s->res_gen    = 0;
    }

    s->in_pass          = true;
    s->has_work         = true;
    s->pass_width       = w;
    s->pass_height      = h;
    s->pass_color_count = desc->color_count;
    for (uint32_t i = 0; i < MD_GPU_MAX_COLOR_TARGETS; ++i) {
        s->pass_color[i] = i < desc->color_count ? desc->color[i].texture->desc.format : MD_GPU_FORMAT_INVALID;
    }
    s->pass_depth     = desc->depth.texture ? desc->depth.texture->desc.format : MD_GPU_FORMAT_INVALID;
    s->bound_pipeline = NULL;

    md_gpu_draw_state_t zero;
    memset(&zero, 0, sizeof(zero));
    md_mtl_apply_draw_state(s, &zero, true);
    const md_gpu_viewport_t vp = {0, 0, (float)w, (float)h, 0, 0};
    md_mtl_apply_viewport(s, &vp);
    const md_gpu_rect_t sc = {0, 0, w, h};
    md_mtl_apply_scissor(s, &sc);
    return true;
}

bool md_gpu_render_end(md_gpu_stream_t s) {
    if (!s) return md_mtl_fail("md_gpu_render_end: null stream");
    if (!s->in_pass) return md_mtl_fail("md_gpu_render_end: stream '%s' has no open render pass", s->label);
    md_mtl_close_encoder(s, true);
    s->in_pass        = false;
    s->bound_pipeline = NULL;
    return true;
}

static void md_mtl_abandon_pass(md_gpu_stream_t s) {
    if (s->in_pass) md_gpu_render_end(s);
}

void md_gpu_set_draw_state(md_gpu_stream_t s, const md_gpu_draw_state_t* state) {
    if (!s) return;
    if (!s->in_pass) { md_mtl_fail("md_gpu_set_draw_state: stream '%s' has no open render pass", s->label); return; }
    md_gpu_draw_state_t n;
    memset(&n, 0, sizeof(n));
    if (state) n = *state;
    if ((unsigned)n.depth_compare > (unsigned)MD_GPU_COMPARE_GREATER || (unsigned)n.cull > (unsigned)MD_GPU_CULL_FRONT) {
        md_mtl_fail("md_gpu_set_draw_state: invalid depth_compare or cull");
        return;
    }
    md_mtl_apply_draw_state(s, &n, false);
}

void md_gpu_set_viewport(md_gpu_stream_t s, const md_gpu_viewport_t* viewport) {
    if (!s) return;
    if (!s->in_pass) { md_mtl_fail("md_gpu_set_viewport: stream '%s' has no open render pass", s->label); return; }
    md_gpu_viewport_t v = {0, 0, (float)s->pass_width, (float)s->pass_height, 0, 0};
    if (viewport) v = *viewport;
    if (!(v.width > 0.0f) || !(v.height > 0.0f)) { md_mtl_fail("md_gpu_set_viewport: width and height must be positive"); return; }
    if (v.min_depth < 0.0f || v.min_depth > 1.0f || v.max_depth < 0.0f || v.max_depth > 1.0f) {
        md_mtl_fail("md_gpu_set_viewport: depth range must lie in [0, 1]");
        return;
    }
    md_mtl_apply_viewport(s, &v);
}

void md_gpu_set_scissor(md_gpu_stream_t s, const md_gpu_rect_t* scissor) {
    if (!s) return;
    if (!s->in_pass) { md_mtl_fail("md_gpu_set_scissor: stream '%s' has no open render pass", s->label); return; }
    md_gpu_rect_t r = {0, 0, s->pass_width, s->pass_height};
    if (scissor) r = *scissor;
    md_mtl_apply_scissor(s, &r);
}

/* ---- Draws ------------------------------------------------------------------ */

static bool md_mtl_draw_check(md_gpu_stream_t s, md_gpu_pipeline_t p, const void* args, size_t args_size, const char* what) {
    if (!s || !p) return md_mtl_fail("%s: null stream or pipeline", what);
    if (p->device != s->device) return md_mtl_fail("%s: pipeline '%s' belongs to another device", what, p->label);
    if (!s->in_pass) return md_mtl_fail("%s: stream '%s' has no open render pass", what, s->label);
    if (p->color_count != s->pass_color_count) {
        return md_mtl_fail("%s: pipeline '%s' has %u colour targets but the pass has %u attachments",
                           what, p->label, p->color_count, s->pass_color_count);
    }
    for (uint32_t i = 0; i < p->color_count; ++i) {
        if (p->color[i] != s->pass_color[i]) {
            return md_mtl_fail("%s: pipeline '%s' colour target %u is %s but the pass attachment is %s",
                               what, p->label, i, md_mtl_format_info(p->color[i]).name, md_mtl_format_info(s->pass_color[i]).name);
        }
    }
    if (p->depth != s->pass_depth) {
        return md_mtl_fail("%s: pipeline '%s' depth format is %s but the pass has %s", what, p->label,
                           p->depth ? md_mtl_format_info(p->depth).name : "none",
                           s->pass_depth ? md_mtl_format_info(s->pass_depth).name : "none");
    }
    if (p->args_size != 0 && args_size != p->args_size) {
        return md_mtl_fail("%s: pipeline '%s' expects a %u-byte argument struct but %zu bytes were passed",
                           what, p->label, p->args_size, args_size);
    }
    if (args_size > 0 && !args) return md_mtl_fail("%s: null args with non-zero size", what);
    return true;
}

/* Arguments into the arena, the root cell pointing at them, and both bound
   to the stages that read a root. See md_mtl_launch_common for why the root
   is a cell of its own. */
static bool md_mtl_draw_bind(md_gpu_stream_t s, md_gpu_pipeline_t p, const void* args, size_t args_size) {
    id<MTLRenderCommandEncoder> enc = s->render_enc;
    md_mtl_declare_residency_render(s, enc);
    if (s->bound_pipeline != p) {
        [enc setRenderPipelineState:p->pso];
        s->bound_pipeline = p;
    }
    if (args_size > 0) {
        uint64_t arg_addr, root_addr, root_off;
        void* arg_host; void* root_host;
        id<MTLBuffer> root_buf = nil;
        if (!md_mtl_arena_alloc(s, args_size, MD_MTL_ARG_ALIGN, &arg_addr, &arg_host, NULL, NULL)) return false;
        memcpy(arg_host, args, args_size);
        if (!md_mtl_arena_alloc(s, sizeof(uint64_t), MD_MTL_ROOT_ALIGN, &root_addr, &root_host, &root_buf, &root_off)) return false;
        memcpy(root_host, &arg_addr, sizeof(arg_addr));
        if (p->vs_root != UINT32_MAX) [enc setVertexBuffer:root_buf offset:root_off atIndex:p->vs_root];
        if (p->fs_root != UINT32_MAX) [enc setFragmentBuffer:root_buf offset:root_off atIndex:p->fs_root];
    }
    return true;
}

bool md_gpu_draw(md_gpu_stream_t s, md_gpu_pipeline_t p, uint32_t vertex_count, uint32_t instance_count,
                 const void* args, size_t args_size) {
    if (!md_mtl_draw_check(s, p, args, args_size, "md_gpu_draw")) return false;
    if (vertex_count == 0 || instance_count == 0) return true;
    if (!md_mtl_draw_bind(s, p, args, args_size)) return false;
    [s->render_enc drawPrimitives:p->primitive vertexStart:0 vertexCount:vertex_count instanceCount:instance_count];
    return true;
}

static bool md_mtl_resolve_indices(md_gpu_stream_t s, md_gpu_addr_t indices, md_gpu_index_type_t type, uint64_t count,
                                   md_mtl_span_t* out, MTLIndexType* out_type, const char* what) {
    if (type != MD_GPU_INDEX_U32 && type != MD_GPU_INDEX_U16) return md_mtl_fail("%s: invalid index type %d", what, (int)type);
    const uint64_t isize = type == MD_GPU_INDEX_U16 ? 2 : 4;
    if (!indices) return md_mtl_fail("%s: null index address", what);
    if (indices % isize != 0) return md_mtl_fail("%s: index address must be %llu-byte aligned", what, (unsigned long long)isize);
    *out_type = type == MD_GPU_INDEX_U16 ? MTLIndexTypeUInt16 : MTLIndexTypeUInt32;
    return md_mtl_resolve(s->device, indices, count * isize, out, what);
}

bool md_gpu_draw_indexed(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                         uint32_t index_count, uint32_t instance_count, const void* args, size_t args_size) {
    const char* what = "md_gpu_draw_indexed";
    if (!md_mtl_draw_check(s, p, args, args_size, what)) return false;
    if (index_count == 0 || instance_count == 0) return true;
    md_mtl_span_t ib;
    MTLIndexType it;
    if (!md_mtl_resolve_indices(s, indices, index_type, index_count, &ib, &it, what)) return false;
    if (!md_mtl_draw_bind(s, p, args, args_size)) return false;
    [s->render_enc drawIndexedPrimitives:p->primitive indexCount:index_count indexType:it
                             indexBuffer:ib.buffer indexBufferOffset:ib.offset instanceCount:instance_count];
    return true;
}

static bool md_mtl_resolve_cmds(md_gpu_stream_t s, md_gpu_addr_t cmds, uint32_t count, uint32_t stride,
                                md_mtl_span_t* out, const char* what) {
    if (!cmds) return md_mtl_fail("%s: null command address", what);
    if (cmds % 4 != 0) return md_mtl_fail("%s: command address must be 4-byte aligned", what);
    return md_mtl_resolve(s->device, cmds, (uint64_t)count * stride, out, what);
}

/* Metal 3 draws one indirect command per call; the layouts match
   MTLDrawPrimitivesIndirectArguments / MTLDrawIndexedPrimitivesIndirectArguments
   exactly, so each call reads a record in place. */
bool md_gpu_draw_indirect(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t cmds, uint32_t count,
                          const void* args, size_t args_size) {
    const char* what = "md_gpu_draw_indirect";
    if (!md_mtl_draw_check(s, p, args, args_size, what)) return false;
    if (count == 0) return true;
    const uint32_t stride = (uint32_t)sizeof(md_gpu_draw_cmd_t);
    md_mtl_span_t cb;
    if (!md_mtl_resolve_cmds(s, cmds, count, stride, &cb, what)) return false;
    if (!md_mtl_draw_bind(s, p, args, args_size)) return false;
    for (uint32_t i = 0; i < count; ++i) {
        [s->render_enc drawPrimitives:p->primitive indirectBuffer:cb.buffer indirectBufferOffset:cb.offset + (uint64_t)i * stride];
    }
    return true;
}

bool md_gpu_draw_indexed_indirect(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                                  md_gpu_addr_t cmds, uint32_t count, const void* args, size_t args_size) {
    const char* what = "md_gpu_draw_indexed_indirect";
    if (!md_mtl_draw_check(s, p, args, args_size, what)) return false;
    if (count == 0) return true;
    const uint32_t stride = (uint32_t)sizeof(md_gpu_draw_indexed_cmd_t);
    md_mtl_span_t cb, ib;
    MTLIndexType it;
    if (!md_mtl_resolve_cmds(s, cmds, count, stride, &cb, what)) return false;
    if (!md_mtl_resolve_indices(s, indices, index_type, 1, &ib, &it, what)) return false;
    if (!md_mtl_draw_bind(s, p, args, args_size)) return false;
    for (uint32_t i = 0; i < count; ++i) {
        [s->render_enc drawIndexedPrimitives:p->primitive indexType:it indexBuffer:ib.buffer indexBufferOffset:ib.offset
                              indirectBuffer:cb.buffer indirectBufferOffset:cb.offset + (uint64_t)i * stride];
    }
    return true;
}

/* ---- Texture to texture -------------------------------------------------------- */

bool md_gpu_copy_texture(md_gpu_stream_t s, md_gpu_texture_t dst, const md_gpu_tex_region_t* dst_region,
                         md_gpu_texture_t src, const md_gpu_tex_region_t* src_region) {
    if (!s || !dst || !src) return md_mtl_fail("md_gpu_copy_texture: null argument");
    if (dst->device != s->device || src->device != s->device) return md_mtl_fail("md_gpu_copy_texture: texture from another device");
    if (dst->desc.format != src->desc.format) {
        return md_mtl_fail("md_gpu_copy_texture: formats differ ('%s' is %s, '%s' is %s); a converting copy is a draw or a kernel",
                           src->label, src->fi.name, dst->label, dst->fi.name);
    }
    if ((dst->desc.type == MD_GPU_TEX_3D) != (src->desc.type == MD_GPU_TEX_3D)) {
        return md_mtl_fail("md_gpu_copy_texture: cannot copy between a 3D texture and a 2D one ('%s' -> '%s')", src->label, dst->label);
    }
    if (!md_mtl_check_no_upload(s, "md_gpu_copy_texture")) return false;
    md_mtl_copy_region_t sr, dr;
    if (!md_mtl_resolve_region(src, src_region, &sr, "md_gpu_copy_texture (src)")) return false;
    md_gpu_tex_region_t d;
    memset(&d, 0, sizeof(d));
    if (dst_region) d = *dst_region;
    d.extent[0] = (uint32_t)sr.size.width;
    d.extent[1] = (uint32_t)sr.size.height;
    d.extent[2] = src->desc.type == MD_GPU_TEX_3D ? (uint32_t)sr.size.depth : sr.layer_count;
    if (dst->desc.type == MD_GPU_TEX_2D && d.extent[2] != 1) {
        return md_mtl_fail("md_gpu_copy_texture: %u layers of '%s' do not fit 2D texture '%s'", d.extent[2], src->label, dst->label);
    }
    if (!md_mtl_resolve_region(dst, &d, &dr, "md_gpu_copy_texture (dst)")) return false;

    id<MTLBlitCommandEncoder> enc = md_mtl_blit_encoder(s);
    if (!enc) return false;
    for (uint32_t l = 0; l < sr.layer_count; ++l) {
        [enc copyFromTexture:src->texture sourceSlice:sr.first_layer + l sourceLevel:sr.mip
                sourceOrigin:sr.origin sourceSize:sr.size
                   toTexture:dst->texture destinationSlice:dr.first_layer + l destinationLevel:dr.mip
           destinationOrigin:dr.origin];
    }
    md_mtl_did_op(s);
    return true;
}

/* =========================================================================
   12. Presentation
   =========================================================================

   A surface is a CAMetalLayer. Acquire takes the next drawable (blocking,
   as nextDrawable does, with its one-second timeout) and wraps its texture;
   present encodes presentDrawable: into the stream's command buffer and
   commits it. The wrapper goes through deferred destruction like any
   texture, so other streams that read it stay safe. Mailbox does not exist
   in Core Animation and falls back to VSYNC, as documented. */

typedef void (*md_mtl_msg_vb_t)(id, SEL, BOOL);

md_gpu_surface_t md_gpu_surface_create(md_gpu_device_t dev, const md_gpu_surface_desc_t* desc) {
    if (!dev || !desc) { md_mtl_fail("md_gpu_surface_create: null argument"); return NULL; }
    const char* label = desc->label ? desc->label : "surface";
    if (desc->system != MD_GPU_WINDOW_COCOA) {
        md_mtl_fail("surface '%s': the Metal backend presents to COCOA windows only (window system %d)", label, (int)desc->system);
        return NULL;
    }
    if (!desc->window) { md_mtl_fail("surface '%s': no window handle", label); return NULL; }
    if (!dev->can_present) { md_mtl_fail("surface '%s': CAMetalLayer is unavailable", label); return NULL; }
    md_gpu_format_t fmt = desc->format ? desc->format : MD_GPU_FORMAT_BGRA8_UNORM;
    if (fmt != MD_GPU_FORMAT_BGRA8_UNORM && fmt != MD_GPU_FORMAT_BGRA8_SRGB && fmt != MD_GPU_FORMAT_RGBA16_FLOAT &&
        fmt != MD_GPU_FORMAT_RGB10A2_UNORM) {
        md_mtl_fail("surface '%s': %s is not supported (CAMetalLayer offers BGRA8_UNORM, BGRA8_SRGB, RGBA16_FLOAT, RGB10A2_UNORM)",
                    label, md_mtl_format_info(fmt).name);
        return NULL;
    }
    if (desc->usage & ~(MD_GPU_TEX_SAMPLED | MD_GPU_TEX_STORAGE | MD_GPU_TEX_RENDER_TARGET)) {
        md_mtl_fail("surface '%s': invalid usage bits", label);
        return NULL;
    }
    if ((desc->usage & MD_GPU_TEX_STORAGE) && md_mtl_format_info(fmt).srgb) {
        md_mtl_fail("surface '%s': STORAGE usage is not supported for %s", label, md_mtl_format_info(fmt).name);
        return NULL;
    }

    md_gpu_surface_t sf = (md_gpu_surface_t)md_alloc(dev->alloc, sizeof(md_gpu_surface));
    if (!sf) { md_mtl_fail("out of memory"); return NULL; }
    memset(sf, 0, sizeof(*sf));
    sf->device = dev;
    sf->desc   = *desc;
    sf->desc.format = fmt;
    sf->desc.usage &= ~MD_GPU_TEX_RENDER_TARGET;
    snprintf(sf->label, sizeof(sf->label), "%s", label);
    sf->desc.label  = sf->label;
    sf->want_width  = desc->width;
    sf->want_height = desc->height;

    @autoreleasepool {
        id obj = (__bridge id)desc->window;
        CAMetalLayer* layer = nil;
        if ([obj isKindOfClass:[CAMetalLayer class]]) {
            layer = (CAMetalLayer*)obj;
        } else if ([obj respondsToSelector:@selector(setWantsLayer:)] && [obj respondsToSelector:@selector(setLayer:)]) {
            /* An NSView: give it a layer of ours. AppKit requires the main
               thread here, which the header states. */
            layer = [CAMetalLayer layer];
            ((md_mtl_msg_vb_t)objc_msgSend)(obj, @selector(setWantsLayer:), YES);
            md_mtl_msg1(obj, @selector(setLayer:), layer);
        }
        if (!layer) {
            md_free(dev->alloc, sf, sizeof(*sf));
            md_mtl_fail("surface '%s': window is neither an NSView nor a CAMetalLayer", label);
            return NULL;
        }
        MD_MTL_RETAIN(layer);
        sf->layer = layer;
        layer.device          = dev->device;
        layer.pixelFormat     = md_mtl_format_info(fmt).fmt;
        /* Copies from the drawable (screenshots, readback) are legal on every
           backend, so the drawable is never framebuffer-only. */
        layer.framebufferOnly = NO;
#if defined(TARGET_OS_OSX) && TARGET_OS_OSX
        layer.displaySyncEnabled = desc->present_mode == MD_GPU_PRESENT_IMMEDIATE ? NO : YES;
#endif
        if (sf->want_width && sf->want_height) {
            layer.drawableSize = CGSizeMake(sf->want_width, sf->want_height);
            sf->width  = sf->want_width;
            sf->height = sf->want_height;
        }
    }

    md_mutex_lock(&dev->device_mutex);
    md_gpu_surface_t* slot = (md_gpu_surface_t*)md_mtl_vec_push(&dev->surfaces, dev->alloc);
    if (slot) *slot = sf;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) {
        MD_MTL_RELEASE(sf->layer);
        md_free(dev->alloc, sf, sizeof(*sf));
        md_mtl_fail("out of memory");
        return NULL;
    }
    return sf;
}

/* Hand the wrapper of an acquired drawable to deferred destruction. Caller
   holds device_mutex. */
static void md_mtl_surface_release_locked(md_gpu_surface_t sf) {
    if (sf->texture) {
        md_mtl_retire_locked(sf->device, MD_MTL_RETIRE_TEXTURE, sf->texture);
        sf->texture = NULL;
    }
    if (sf->drawable) {
        MD_MTL_RELEASE(sf->drawable);
        sf->drawable = nil;
    }
    sf->acquired_stream = NULL;
}

void md_gpu_surface_destroy(md_gpu_surface_t sf) {
    if (!sf) return;
    md_gpu_device_t dev = sf->device;
    md_mutex_lock(&dev->device_mutex);
    md_mtl_vec_remove_ptr(&dev->surfaces, sf);
    md_mtl_surface_release_locked(sf);
    md_mutex_unlock(&dev->device_mutex);
    /* The layer is retained by its view and by any command buffer presenting
       from it; ours can go now. */
    MD_MTL_RELEASE(sf->layer);
    md_free(dev->alloc, sf, sizeof(*sf));
}

void md_gpu_surface_resize(md_gpu_surface_t sf, uint32_t width, uint32_t height) {
    if (!sf) return;
    md_mutex_lock(&sf->device->device_mutex);
    sf->want_width  = width;
    sf->want_height = height;
    md_mutex_unlock(&sf->device->device_mutex);
}

md_gpu_texture_t md_gpu_surface_acquire(md_gpu_stream_t s, md_gpu_surface_t sf) {
    if (!s || !sf) { md_mtl_fail("md_gpu_surface_acquire: null argument"); return NULL; }
    md_gpu_device_t dev = s->device;
    if (sf->device != dev) { md_mtl_fail("md_gpu_surface_acquire: surface '%s' belongs to another device", sf->label); return NULL; }
    if (s->kind != MD_GPU_STREAM_GRAPHICS) { md_mtl_fail("md_gpu_surface_acquire: stream '%s' is not a GRAPHICS stream", s->label); return NULL; }
    if (!md_mtl_not_in_pass(s, "md_gpu_surface_acquire") || !md_mtl_check_no_upload(s, "md_gpu_surface_acquire")) return NULL;
    if (sf->drawable) { md_mtl_fail("md_gpu_surface_acquire: surface '%s' already has an acquired image; present it first", sf->label); return NULL; }

    md_mutex_lock(&dev->device_mutex);
    const uint32_t w = sf->want_width, h = sf->want_height;
    md_mutex_unlock(&dev->device_mutex);
    if (w == 0 || h == 0) return NULL;     /* minimised: skip the frame */

    md_gpu_texture_t t = NULL;
    @autoreleasepool {
        if (w != sf->width || h != sf->height) {
            sf->layer.drawableSize = CGSizeMake(w, h);
            sf->width  = w;
            sf->height = h;
        }
        id<CAMetalDrawable> drawable = [sf->layer nextDrawable];
        if (!drawable) { md_mtl_fail("md_gpu_surface_acquire: surface '%s': nextDrawable timed out", sf->label); return NULL; }
        id<MTLTexture> tex = drawable.texture;

        t = (md_gpu_texture_t)md_alloc(dev->alloc, sizeof(md_gpu_texture));
        if (!t) { md_mtl_fail("out of memory"); return NULL; }
        memset(t, 0, sizeof(*t));
        MD_MTL_RETAIN(drawable);
        MD_MTL_RETAIN(tex);
        t->device   = dev;
        t->texture  = tex;
        t->external = true;
        t->fi       = md_mtl_format_info(sf->desc.format);
        t->desc.type            = MD_GPU_TEX_2D;
        t->desc.format          = sf->desc.format;
        t->desc.usage           = MD_GPU_TEX_RENDER_TARGET | sf->desc.usage;
        t->desc.width           = (uint32_t)tex.width;
        t->desc.height          = (uint32_t)tex.height;
        t->desc.depth_or_layers = 1;
        t->desc.mip_levels      = 1;
        snprintf(t->label, sizeof(t->label), "%s[%llu]", sf->label, (unsigned long long)sf->frame++);
        t->desc.label = t->label;
        if (sf->desc.usage & MD_GPU_TEX_SAMPLED) t->sampled_handle = (uint64_t)tex.gpuResourceID._impl;
        if (sf->desc.usage & MD_GPU_TEX_STORAGE) {
            t->storage_handles = (uint64_t*)md_alloc(dev->alloc, sizeof(uint64_t));
            if (t->storage_handles) t->storage_handles[0] = (uint64_t)tex.gpuResourceID._impl;
        }
        md_mutex_lock(&dev->device_mutex);
        if (sf->desc.usage & (MD_GPU_TEX_SAMPLED | MD_GPU_TEX_STORAGE)) {
            /* Reached by resource id from shaders, so it must be resident. */
            md_mtl_make_resident_locked(dev, tex);
            t->resident = true;
        }
        sf->drawable        = drawable;
        sf->texture         = t;
        sf->acquired_stream = s;
        md_mutex_unlock(&dev->device_mutex);
    }
    return t;
}

bool md_gpu_surface_present(md_gpu_stream_t s, md_gpu_surface_t sf) {
    if (!s || !sf) return md_mtl_fail("md_gpu_surface_present: null argument");
    if (!sf->drawable) return md_mtl_fail("md_gpu_surface_present: surface '%s' has no acquired image", sf->label);
    if (sf->acquired_stream != s) return md_mtl_fail("md_gpu_surface_present: surface '%s' was acquired on another stream", sf->label);
    if (!md_mtl_not_in_pass(s, "md_gpu_surface_present") || !md_mtl_check_no_upload(s, "md_gpu_surface_present")) return false;

    bool ok = md_mtl_stream_ensure_cmd(s);
    if (ok) {
        md_mtl_close_encoder(s, true);
        [s->cmd presentDrawable:sf->drawable];
        s->has_work = true;
        ok = md_mtl_stream_submit(s);
    }
    md_mutex_lock(&s->device->device_mutex);
    md_mtl_surface_release_locked(sf);
    md_mutex_unlock(&s->device->device_mutex);
    return ok;
}

/* =========================================================================
   13. Host callbacks and polling
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
    if (!md_mtl_not_in_pass(s, "md_gpu_launch_host_fn")) return false;
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
    md_mtl_process_pending_frees_locked(dev);
    md_mutex_unlock(&dev->device_mutex);
    return fired;
}

/* =========================================================================
   14. Device
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

/* Every Metal device. Call inside an @autoreleasepool. macOS can have several
   (Intel Macs with a discrete GPU, eGPUs); elsewhere only the default one. */
static NSArray* md_mtl_all_devices(void) {
#if defined(TARGET_OS_OSX) && TARGET_OS_OSX
#  if __has_feature(objc_arc)
    return MTLCopyAllDevices();
#  else
    return [MTLCopyAllDevices() autorelease];
#  endif
#else
    id<MTLDevice> d = MTLCreateSystemDefaultDevice();
    if (!d) return @[];
    NSArray* a = @[d];
#  if !__has_feature(objc_arc)
    [d release];
#  endif
    return a;
#endif
}

/* Metal has no PCI ids; the name tells the vendor. */
static uint32_t md_mtl_vendor_id(const char* name) {
    char low[256];
    size_t i = 0;
    for (; name[i] && i + 1 < sizeof(low); ++i) low[i] = (char)tolower((unsigned char)name[i]);
    low[i] = 0;
    if (strstr(low, "apple"))                           return 0x106B;
    if (strstr(low, "amd") || strstr(low, "radeon"))    return 0x1002;
    if (strstr(low, "intel"))                           return 0x8086;
    if (strstr(low, "nvidia") || strstr(low, "geforce")) return 0x10DE;
    return 0;
}

static md_gpu_device_type_t md_mtl_device_type(id<MTLDevice> d) {
    if ([d hasUnifiedMemory]) return MD_GPU_DEVICE_TYPE_INTEGRATED;
#if defined(TARGET_OS_OSX) && TARGET_OS_OSX
    if ([d isLowPower]) return MD_GPU_DEVICE_TYPE_INTEGRATED;
#endif
    return MD_GPU_DEVICE_TYPE_DISCRETE;
}

static void md_mtl_adapter_info(id<MTLDevice> d, md_gpu_adapter_info_t* out) {
    memset(out, 0, sizeof(*out));
    snprintf(out->name, sizeof(out->name), "%s", [[d name] UTF8String]);
    snprintf(out->driver, sizeof(out->driver), "Metal, %s",
             [[[NSProcessInfo processInfo] operatingSystemVersionString] UTF8String]);
    out->vendor_id = md_mtl_vendor_id(out->name);
    out->type      = md_mtl_device_type(d);
    out->usable    = true;
}

uint32_t md_gpu_enumerate_adapters(md_gpu_adapter_info_t* out, uint32_t max) {
    uint32_t n = 0;
    @autoreleasepool {
        NSArray* all = md_mtl_all_devices();
        n = (uint32_t)[all count];
        for (uint32_t i = 0; i < n && i < max; ++i) {
            id<MTLDevice> d = all[i];
            md_mtl_adapter_info(d, &out[i]);
        }
    }
    return n;
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
    md_mtl_vec_init(&dev->registry, sizeof(md_mtl_range_t*));
    md_mtl_vec_init(&dev->pending_frees, sizeof(md_mtl_pending_free_t));
    md_mtl_vec_init(&dev->textures, sizeof(md_gpu_texture_t));
    for (uint32_t k = 0; k < MD_GPU_MEM_KIND_COUNT; ++k) {
        md_tlsf_init(&dev->heaps[k].tlsf, alloc, MD_MTL_HEAP_ALIGN);
        md_mtl_vec_init(&dev->heaps[k].chunks, sizeof(md_mtl_chunk_t*));
    }
    dev->heap_cache_limit = (desc && desc->heap_cache_limit) ? desc->heap_cache_limit : MD_MTL_HEAP_CACHE_DEFAULT;
    md_mtl_vec_init(&dev->kernels,  sizeof(md_gpu_kernel_t));
    md_mtl_vec_init(&dev->pipelines, sizeof(md_gpu_pipeline_t));
    md_mtl_vec_init(&dev->surfaces, sizeof(md_gpu_surface_t));
    md_mtl_vec_init(&dev->streams,  sizeof(md_gpu_stream_t));
    md_mtl_vec_init(&dev->hostfns,  sizeof(md_mtl_hostfn_t));
    md_mtl_vec_init(&dev->retires,  sizeof(md_mtl_retire_t));
    md_mtl_vec_init(&dev->live_res, sizeof(id));
    md_mutex_init(&dev->queue_mutex);
    md_mutex_init(&dev->device_mutex);

    @autoreleasepool {
        /* Adapter: the system default unless an adapter or a preference is
           given (md_gpu_device_desc_t.adapter, MD_GPU_DEVICE). */
        NSArray* all = md_mtl_all_devices();
        enum { MD_MTL_MAX_ADAPTERS = 16 };
        md_gpu_adapter_info_t infos[MD_MTL_MAX_ADAPTERS];
        const uint32_t n_all = (uint32_t)MIN([all count], (NSUInteger)MD_MTL_MAX_ADAPTERS);
        for (uint32_t i = 0; i < n_all; ++i) {
            id<MTLDevice> d = all[i];
            md_mtl_adapter_info(d, &infos[i]);
        }
        char sel[128];
        const bool explicit_pick = md_gpu_sel_get(desc, sel, sizeof(sel)) ||
                                   (desc && desc->preference != MD_GPU_DEVICE_PREFER_DEFAULT);
        id<MTLDevice> mtl = nil;
        if (explicit_pick) {
            char why[2048];
            const int pick = md_gpu_sel_pick(infos, n_all, desc, why, sizeof(why));
            if (pick < 0) {
                md_mtl_fail("md_gpu_device_create: %s", why);
            } else {
                mtl = all[(NSUInteger)pick];
#if !__has_feature(objc_arc)
                [mtl retain];   /* owned like the result of MTLCreateSystemDefaultDevice */
#endif
                dev->adapter_index = (uint32_t)pick;
            }
        } else {
            mtl = MTLCreateSystemDefaultDevice();
            if (!mtl) md_mtl_fail("MTLCreateSystemDefaultDevice returned nil");
            for (uint32_t i = 0; mtl && i < n_all; ++i) {
                id<MTLDevice> d = all[i];
                if ([d registryID] == [mtl registryID]) dev->adapter_index = i;
            }
        }
        if (mtl) {
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
    dev->default_graphics = md_mtl_stream_create_internal(dev, MD_GPU_STREAM_GRAPHICS, "default graphics", true);
    if (!dev->default_compute || !dev->default_transfer || !dev->default_graphics) goto fail;
    dev->can_present = NSClassFromString(@"CAMetalLayer") != Nil;
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
    info->supports_graphics        = true;
    info->supports_present         = dev->can_present;
    snprintf(info->name, sizeof(info->name), "%s", [[dev->device name] UTF8String]);
    info->vendor_id                = md_mtl_vendor_id(info->name);
    info->type                     = md_mtl_device_type(dev->device);
    info->adapter_index            = dev->adapter_index;
    snprintf(info->driver, sizeof(info->driver), "Metal, %s",
             [[[NSProcessInfo processInfo] operatingSystemVersionString] UTF8String]);
    return true;
}

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;

    /* Flush and idle every stream, then fire what is pending. */
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
        s->upload_open = false;
        md_mtl_abandon_pass(s);
        md_mtl_stream_submit(s);
    }
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t s = MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i);
        if (s->submitted_value > 0) md_mtl_event_wait(s->timeline, s->submitted_value);
    }
    md_gpu_device_poll(dev);

    /* Everything the caller did not destroy, the device does. */
    while (dev->textures.count > 0) md_gpu_texture_destroy(MD_MTL_VEC_AT(dev->textures, md_gpu_texture_t, 0));
    while (dev->kernels.count > 0) md_gpu_kernel_destroy(MD_MTL_VEC_AT(dev->kernels, md_gpu_kernel_t, 0));
    while (dev->pipelines.count > 0) md_gpu_pipeline_destroy(MD_MTL_VEC_AT(dev->pipelines, md_gpu_pipeline_t, 0));
    while (dev->surfaces.count > 0) md_gpu_surface_destroy(MD_MTL_VEC_AT(dev->surfaces, md_gpu_surface_t, 0));
    if (dev->make_grid_kernel) md_mtl_kernel_free(dev, dev->make_grid_kernel);
    if (dev->byte_op_kernel)   md_mtl_kernel_free(dev, dev->byte_op_kernel);

    md_mutex_lock(&dev->device_mutex);
    md_mtl_process_retires_locked(dev, true);
    md_mutex_unlock(&dev->device_mutex);

    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_mtl_stream_free(MD_MTL_VEC_AT(dev->streams, md_gpu_stream_t, i));
    }
    md_mtl_vec_free(&dev->streams, alloc);
    md_mutex_lock(&dev->device_mutex);
    md_mtl_heaps_free_locked(dev);
    md_mutex_unlock(&dev->device_mutex);
    md_mtl_vec_free(&dev->textures, alloc);
    md_mtl_vec_free(&dev->kernels, alloc);
    md_mtl_vec_free(&dev->pipelines, alloc);
    md_mtl_vec_free(&dev->surfaces, alloc);
    for (int c = 0; c < 8; ++c) for (int w = 0; w < 2; ++w) MD_MTL_RELEASE(dev->depth_states[c][w]);

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

/*
md_gpu.h

A thin layer over Vulkan and Metal (and, as an alternative compute-only build,
CUDA). One backend per build. The compute side is CUDA-shaped; the raster side,
added later, follows the "no graphics API" model: pointers, one root argument
pointer, a bindless heap, coarse stage barriers and no resource state.

The model, in full
------------------
  * Work issued into a stream executes in issue order. Streams are unordered
    with respect to each other unless joined with md_gpu_stream_wait().

  * Ordering inside a stream is IMPLICIT by default: every operation observes
    all writes made by the operations before it in that stream, exactly as in
    CUDA. A stream may switch to EXPLICIT ordering, in which md_gpu inserts
    nothing and the caller places md_gpu_barrier(producer_stages,
    consumer_stages). There are no resource lists, no layouts and no resource
    state in either mode.

  * Device memory is a 64-bit GPU address, md_gpu_addr_t. Host-visible
    allocations also hand back a CPU pointer. The two are distinct values of
    distinct types and neither converts to the other implicitly -- which is
    also why copies name their direction instead of inferring it.

  * A texture is an object. Shaders never see it; they see a storage handle, a
    sampled handle or a sampler handle, each placed in the argument struct and
    received as a Slang DescriptorHandle<T>.

  * Every call that is ordered against the GPU takes the stream as its first
    argument. Nothing blocks the calling thread unless it says so:
    md_gpu_stream_sync, md_gpu_sync_wait, md_gpu_stream_destroy (its own work
    only) and md_gpu_device_destroy. In particular, creating or destroying
    textures, pools and kernels never waits for work in flight -- a compute
    job spanning many frames stalls nothing but its own stream.

Correspondence with CUDA
------------------------
    cudaStreamCreate            md_gpu_stream_create
    cudaMallocFromPoolAsync     md_gpu_malloc
    cudaFreeAsync               md_gpu_free
    cudaMemcpyAsync (H2D)       md_gpu_upload / md_gpu_upload_begin+end
    cudaMemcpyAsync (D2D)       md_gpu_copy
    cudaMemcpyAsync (D2H)       md_gpu_copy into MD_GPU_MEM_HOST_READ memory
    cudaMemsetAsync             md_gpu_memset
    kernel<<<g,b,0,s>>>(args)   md_gpu_launch
    cudaEventRecord             md_gpu_stream_record
    cudaStreamWaitEvent         md_gpu_stream_wait
    cudaStreamSynchronize       md_gpu_stream_sync
    cudaEventQuery              md_gpu_sync_is_complete
    cudaLaunchHostFunc          md_gpu_launch_host_fn
    cudaSurfaceObject_t         md_gpu_storage_tex_t
    cudaTextureObject_t         md_gpu_sampled_tex_t

Kernel arguments
----------------
A kernel receives one pointer to a caller-defined argument struct. The backend
copies the struct into device memory and passes its address in an 8-byte push
constant, so there is no size limit and no portability cliff:

    struct Args { uint n; uint _pad; float* dst; };
    MD_KERNEL_ARGS(Args);                   // from md_gpu.slang

    [shader("compute")][numthreads(64,1,1)]
    void main(uint3 tid : SV_DispatchThreadID) {
        Args a = MD_ARGS;
        ...
    }

compile_gpu_shaders() emits, per entry point, a function returning a ready
md_gpu_kernel_desc_t (code, group size, argument-struct size), so a call site
never repeats [numthreads] by hand.

Threading
---------
  * A stream is used by one thread at a time. Different streams may be used
    concurrently from different threads.
  * Pool, allocation, texture, sampler and kernel creation/destruction are
    thread-safe.
  * Host callbacks run inside md_gpu_device_poll(), on the thread calling it.
*/

#ifndef MD_GPU_H
#define MD_GPU_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__cplusplus)
#  define MD_GPU_STATIC_ASSERT(c, m) static_assert(c, m)
#else
#  define MD_GPU_STATIC_ASSERT(c, m) _Static_assert(c, m)
#endif

MD_GPU_STATIC_ASSERT(sizeof(void*) == 8, "md_gpu requires a 64-bit target");

/* =========================================================================
   Handles and value types
   ========================================================================= */

typedef struct md_gpu_device*  md_gpu_device_t;
typedef struct md_gpu_stream*  md_gpu_stream_t;
typedef struct md_gpu_pool*    md_gpu_pool_t;
typedef struct md_gpu_texture* md_gpu_texture_t;   /* identity; host side only */
typedef struct md_gpu_kernel*  md_gpu_kernel_t;

/* A GPU virtual address. Byte arithmetic is valid (`base + offsetof(T, f)`),
   and it is the type of every pointer field in a C argument-struct mirror, so
   no casts are needed there. Zero is null. Not dereferenceable on the host. */
typedef uint64_t md_gpu_addr_t;

/* An allocation. `cpu` is non-NULL only for memory from an MD_GPU_MEM_HOST_*
   pool, and then addresses the same bytes as `gpu`. */
typedef struct md_gpu_mem_t {
    md_gpu_addr_t gpu;
    void*         cpu;
} md_gpu_mem_t;

/* Shader-visible handles. Distinct C types, so that a storage handle cannot be
   assigned to a sampled field by accident. Each is 8 bytes, 8-aligned, and is
   received in Slang as:

       md_gpu_storage_tex_t   DescriptorHandle<RWTexture2D / RWTexture3D / RWTexture2DArray<T>>
       md_gpu_sampled_tex_t   DescriptorHandle<Texture2D / Texture3D / Texture2DArray<T>>
       md_gpu_sampler_t       DescriptorHandle<SamplerState>

   A zero handle is null. Slang type-checks the handle against the resource
   type, so using a 3D storage handle where a 2D sampled texture is expected is
   a compile error rather than silent corruption. */
typedef struct md_gpu_storage_tex_t { uint64_t handle; } md_gpu_storage_tex_t;
typedef struct md_gpu_sampled_tex_t { uint64_t handle; } md_gpu_sampled_tex_t;
typedef struct md_gpu_sampler_t     { uint64_t handle; } md_gpu_sampler_t;

/* A point on a stream's timeline. Value type: copy it, store it, pass it by
   value. A zero-initialised sync is the "none" sync -- waiting on it is a
   no-op, so an optional dependency needs no special case. */
typedef struct md_gpu_sync_t {
    md_gpu_stream_t stream;
    uint64_t        value;
} md_gpu_sync_t;

static inline md_gpu_sync_t md_gpu_sync_none(void) {
    md_gpu_sync_t s;
    s.stream = 0;
    s.value  = 0;
    return s;
}

static inline bool md_gpu_sync_is_valid(md_gpu_sync_t s) {
    return s.stream != 0 && s.value != 0;
}

/* =========================================================================
   Vector types for argument structs
   =========================================================================

   A kernel's argument struct is mirrored in C, and the two shader backends do
   not lay vectors out identically. Slang emits SPIR-V with scalar layout -- a
   vector's alignment is its component's -- while MSL gives vec2 8/8 and vec4
   16/16. Crucially the *sizes* agree; only the alignment differs. So a single
   C struct is correct on both backends exactly when every vector member sits
   at an offset satisfying the stricter (Metal) alignment:

       md_gpu_*2, handle types     offset must be a multiple of 8
       md_gpu_*4, md_gpu_float4x4  offset must be a multiple of 16

   That is the entire ABI rule. It is a property of the struct, not of the
   backend, which is why the types below carry no #ifdef and this header needs
   to know nothing about which backend it was built against. The types declare
   the alignment they require, so C places them correctly on its own; the
   shader side is held to the same rule by tools/check_gpu_arg_layout.py, which
   compiles every kernel for both targets and rejects any argument struct whose
   two layouts disagree. A struct that violates the rule is a build error, not
   a silently wrong number several launches downstream.

   NO THREE-COMPONENT VECTORS IN ARGUMENT STRUCTS
   ----------------------------------------------
   uint3 is 12 bytes in SPIR-V and 16 in MSL, and no amount of padding closes
   that gap: MSL's slack sits *inside* the vector, SPIR-V's sits after it, so
   every following member shifts by 4 on one target only. Landing the vector on
   a 16-byte boundary fixes where it starts and nothing else. There is
   deliberately no md_gpu_uint3 -- use a 4-vector and ignore .w, or three
   scalars, whichever reads better:

       struct Args {                    // Slang
           uint4  dim;                  // xyz used
           float  scale;
           uint   _pad;
           float* dst;
       };

       typedef struct {                 // C, correct on every backend
           md_gpu_uint4 dim;
           float        scale;
           uint32_t     _pad;
           md_gpu_addr_t dst;
       } my_args_t;

   Inside shader code -- locals, groupshared, arithmetic -- 3-vectors are fine
   and cost nothing; `a.dim.xyz` is the usual way to read one back out. The
   rule constrains the argument struct only.

   Scalars and 8-byte device addresses place themselves and need no help.
   Texture and sampler handles are the one non-vector member the rule covers:
   Slang lowers DescriptorHandle<T> to two 32-bit words, so SPIR-V aligns it to
   4 while MSL aligns it to 8. A handle preceded by an odd number of 32-bit
   scalars therefore lands 4 bytes apart on the two targets -- pad to an 8-byte
   offset, exactly as for a 2-vector.
   (The alignment specifier sits on the first member only -- applying it to a
   multi-declarator line would align every declarator and inflate the struct.) */

#if defined(__cplusplus)
#  define MD_GPU_ALIGNAS(n) alignas(n)
#elif defined(_MSC_VER)
   /* _Alignas needs /std:c11 or later in MSVC's C mode; __declspec(align)
      works regardless of the language level. */
#  define MD_GPU_ALIGNAS(n) __declspec(align(n))
#else
#  define MD_GPU_ALIGNAS(n) _Alignas(n)
#endif

#define MD_GPU_DEFINE_VEC(name, T)                                             \
    typedef struct { MD_GPU_ALIGNAS(8)  T x; T y; }       name##2;             \
    typedef struct { MD_GPU_ALIGNAS(16) T x; T y, z, w; } name##4

MD_GPU_DEFINE_VEC(md_gpu_uint,  uint32_t);
MD_GPU_DEFINE_VEC(md_gpu_int,   int32_t);
MD_GPU_DEFINE_VEC(md_gpu_float, float);

/* Slang lowers a float4x4 to four 16-byte columns: 64 bytes on both targets.
   Same rule as a 4-vector -- store column-major, matching Slang. */
typedef struct {
    MD_GPU_ALIGNAS(16) float m[16];
} md_gpu_float4x4;

#undef MD_GPU_DEFINE_VEC
#undef MD_GPU_ALIGNAS

/* A toolchain that quietly drops the alignment specifier would reintroduce
   exactly the bug these types exist to prevent, so say so at compile time.
   offsetof rather than _Alignof: MSVC's C mode did not accept _Alignof until
   recently, and the probe works everywhere. */
struct md_gpu_align_probe2_ { char c; md_gpu_uint2    v; };
struct md_gpu_align_probe4_ { char c; md_gpu_uint4    v; };
struct md_gpu_align_probem_ { char c; md_gpu_float4x4 v; };

MD_GPU_STATIC_ASSERT(sizeof(md_gpu_uint2)    ==  8, "md_gpu_*2 must be 8 bytes");
MD_GPU_STATIC_ASSERT(sizeof(md_gpu_uint4)    == 16, "md_gpu_*4 must be 16 bytes");
MD_GPU_STATIC_ASSERT(sizeof(md_gpu_float4x4) == 64, "md_gpu_float4x4 must be 64 bytes");
MD_GPU_STATIC_ASSERT(offsetof(struct md_gpu_align_probe2_, v) ==  8,
                     "md_gpu_*2 must be 8-byte aligned");
MD_GPU_STATIC_ASSERT(offsetof(struct md_gpu_align_probe4_, v) == 16,
                     "md_gpu_*4 must be 16-byte aligned");
MD_GPU_STATIC_ASSERT(offsetof(struct md_gpu_align_probem_, v) == 16,
                     "md_gpu_float4x4 must be 16-byte aligned");


/* Handles are plain 8-byte values in argument structs. */
MD_GPU_STATIC_ASSERT(sizeof(md_gpu_storage_tex_t) == 8, "handle must be 8 bytes");
MD_GPU_STATIC_ASSERT(sizeof(md_gpu_sampled_tex_t) == 8, "handle must be 8 bytes");
MD_GPU_STATIC_ASSERT(sizeof(md_gpu_sampler_t)     == 8, "handle must be 8 bytes");
#undef MD_GPU_STATIC_ASSERT

/* Launch geometry, in thread groups. */
typedef struct md_gpu_grid_t {
    uint32_t x, y, z;
} md_gpu_grid_t;

static inline md_gpu_grid_t md_gpu_grid(uint32_t x, uint32_t y, uint32_t z) {
    md_gpu_grid_t g;
    g.x = x; g.y = y; g.z = z;
    return g;
}

/* =========================================================================
   Errors
   ========================================================================= */

/* Human-readable description of the most recent failure on the calling
   thread, or NULL. Owned by md_gpu; valid until the next failing call on the
   same thread. */
const char* md_gpu_last_error(void);

/* =========================================================================
   Device
   ========================================================================= */

typedef struct md_gpu_device_desc_t {
    /* Allocator for host-side allocations. NULL selects the default heap
       allocator. Must outlive the device. */
    struct md_allocator_i* alloc;

    /* Request backend validation (Vulkan validation layers). On Metal it must
       be enabled from the environment; see md_gpu_metal.m. */
    bool enable_validation;

    const char* label;
} md_gpu_device_desc_t;

typedef struct md_gpu_device_info_t {
    /* False implies UMA / integrated. An allocation hint only. */
    bool     is_discrete;
    uint32_t max_threads_per_group;
    uint32_t preferred_group_multiple;   /* warp / SIMD width */
    char     name[256];
} md_gpu_device_info_t;

/* `desc` may be NULL for defaults. Returns NULL on failure; see
   md_gpu_last_error(). */
md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc);

/* Waits for every stream to go idle, then destroys the device and everything
   created from it: streams, pools, allocations, textures, samplers, kernels. */
void md_gpu_device_destroy(md_gpu_device_t device);

bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* out_info);

/* Fires host callbacks whose sync point has completed and releases objects
   whose deferred destruction has become safe. Callbacks run on the calling
   thread and nowhere else. Returns the number of callbacks fired.

   Call once per frame. Cheap when there is nothing to do. */
uint32_t md_gpu_device_poll(md_gpu_device_t device);

/* =========================================================================
   Streams
   ========================================================================= */

typedef enum md_gpu_stream_kind_t {
    MD_GPU_STREAM_COMPUTE,    /* async compute; work may span many frames   */
    MD_GPU_STREAM_TRANSFER,   /* prefers a dedicated DMA engine if present;
                                 copies and fills only, no kernel launches  */
} md_gpu_stream_kind_t;

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t device, md_gpu_stream_kind_t kind, const char* label);

/* Waits for this stream's own work to complete, then destroys it. Other
   streams are not waited for. */
void md_gpu_stream_destroy(md_gpu_stream_t stream);

/* Device-owned default streams. Always valid; never destroyed by the caller. */
md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t device, md_gpu_stream_kind_t kind);

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t stream);

/* cudaEventRecord: submit whatever is pending and return the sync point that
   is signalled when it completes. With nothing pending this returns the sync
   of the most recent submission (still a correct "everything so far" point),
   or the none sync if the stream has never submitted. */
md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t stream);

/* cudaStreamWaitEvent: work issued into `stream` after this call waits for
   `sync`. Work already issued is unaffected. A none sync, a sync from `stream`
   itself and an already completed sync are no-ops. */
void md_gpu_stream_wait(md_gpu_stream_t stream, md_gpu_sync_t sync);

/* Submit pending work without blocking. */
void md_gpu_stream_flush(md_gpu_stream_t stream);

/* cudaStreamSynchronize: submit pending work and block until it completes. */
void md_gpu_stream_sync(md_gpu_stream_t stream);

bool md_gpu_sync_is_complete(md_gpu_sync_t sync);
void md_gpu_sync_wait(md_gpu_sync_t sync);

/* ---- Ordering ---------------------------------------------------------------

   IMPLICIT (the default): each operation is ordered after everything before it
   in the stream.

   EXPLICIT: md_gpu inserts nothing between operations. The caller states the
   producer and consumer stages with md_gpu_barrier -- no resource lists, no
   layouts. Switching back to IMPLICIT orders the next operation after
   everything before it, so an explicit region can never leak unordered work
   past its end.

       md_gpu_stream_set_ordering(s, MD_GPU_ORDER_EXPLICIT);
       md_gpu_launch(s, k_a, ...);                 // independent of k_b,
       md_gpu_launch(s, k_b, ...);                 //   may overlap
       md_gpu_barrier(s, MD_GPU_STAGE_COMPUTE, MD_GPU_STAGE_COMPUTE);
       md_gpu_launch(s, k_consume_both, ...);
       md_gpu_stream_set_ordering(s, MD_GPU_ORDER_IMPLICIT);

   Vulkan: one global VkMemoryBarrier2 with the given stage masks. Metal 3:
   serial encoders already order everything, so md_gpu_barrier is a no-op and
   EXPLICIT is merely slower than possible, never incorrect. Code written for
   EXPLICIT is correct on every backend. */

typedef enum md_gpu_ordering_t {
    MD_GPU_ORDER_IMPLICIT,
    MD_GPU_ORDER_EXPLICIT,
} md_gpu_ordering_t;

typedef uint32_t md_gpu_stage_flags_t;
enum {
    MD_GPU_STAGE_TRANSFER = 1u << 0,   /* copy, upload, memset, texture copies   */
    MD_GPU_STAGE_COMPUTE  = 1u << 1,   /* kernel launches                        */
    MD_GPU_STAGE_INDIRECT = 1u << 2,   /* consumer side: reading indirect grids  */
    /* Raster stages arrive with the raster API. */
    MD_GPU_STAGE_ALL      = 0xFFFFFFFFu,
};

void md_gpu_stream_set_ordering(md_gpu_stream_t stream, md_gpu_ordering_t ordering);
md_gpu_ordering_t md_gpu_stream_ordering(md_gpu_stream_t stream);

/* Everything `producers` wrote before this point is visible to `consumers`
   after it. Valid in either mode; redundant, but harmless, in IMPLICIT. */
void md_gpu_barrier(md_gpu_stream_t stream, md_gpu_stage_flags_t producers, md_gpu_stage_flags_t consumers);

/* =========================================================================
   Memory
   ========================================================================= */

typedef enum md_gpu_mem_kind_t {
    MD_GPU_MEM_DEVICE,       /* device-local; no CPU pointer                  */
    MD_GPU_MEM_HOST_WRITE,   /* CPU-written, write-combined; device-local too
                                where the platform allows (UMA / ReBAR):
                                uploads and per-frame data                   */
    MD_GPU_MEM_HOST_READ,    /* CPU-cached: readback destinations            */
} md_gpu_mem_kind_t;

/* A pool is the space allocations are drawn from. It serves exactly one kind
   of memory, and it groups lifetimes: destroying or resetting it releases
   everything drawn from it, textures included. */
typedef struct md_gpu_pool_desc_t {
    md_gpu_mem_kind_t kind;
    /* Bytes a pool keeps cached after md_gpu_free for reuse without a new
       device allocation. 0 means no limit: cached memory is only returned to
       the driver by md_gpu_pool_trim or md_gpu_pool_destroy. */
    uint64_t          cache_limit;
    const char*       label;
} md_gpu_pool_desc_t;

md_gpu_pool_t     md_gpu_pool_create(md_gpu_device_t device, const md_gpu_pool_desc_t* desc);
md_gpu_mem_kind_t md_gpu_pool_kind(md_gpu_pool_t pool);

/* Never blocks. Every allocation and texture drawn from the pool is released
   once every stream has completed the work issued before this call. Work
   issued afterwards that still references the pool is a caller error. */
void md_gpu_pool_destroy(md_gpu_pool_t pool);

/* Free everything the pool has handed out, in one call, without returning its
   memory to the driver -- the CPU arena-reset pattern. Allocations become
   reusable at this point in `stream`, exactly like md_gpu_free; textures are
   released like md_gpu_texture_destroy. Every address and texture previously
   obtained from the pool dangles once this returns. */
void md_gpu_pool_reset(md_gpu_stream_t stream, md_gpu_pool_t pool);

/* Release cached (free and idle) memory down to `keep_bytes`. */
void md_gpu_pool_trim(md_gpu_pool_t pool, uint64_t keep_bytes);

typedef struct md_gpu_pool_stats_t {
    uint64_t bytes_in_use;      /* handed out right now                       */
    uint64_t bytes_reserved;    /* committed by the pool, in use or cached    */
    uint64_t bytes_cached;      /* reserved - in_use                          */
    uint64_t bytes_peak_in_use; /* high-water mark, for sizing                */
    uint32_t blocks_in_use;
    uint32_t blocks_cached;
    uint64_t alloc_count;       /* md_gpu_malloc calls served                 */
    uint64_t reuse_count;       /* of those, served from cache. A ratio near
                                   1 means the pool is doing its job          */
} md_gpu_pool_stats_t;

void md_gpu_pool_stats(md_gpu_pool_t pool, md_gpu_pool_stats_t* out_stats);

/* cudaMallocFromPoolAsync. The allocation is usable by work issued into
   `stream` after this call. `.gpu == 0` on failure. Never waits on another
   stream: a cached block freed elsewhere is reused only once its free point
   has completed. */
md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t stream, md_gpu_pool_t pool, size_t size);

/* cudaFreeAsync. The memory returns to its pool at this point in `stream`, so
   later work in the same stream may reuse it with no synchronisation. A zero
   address is a no-op; `stream` is required. */
void md_gpu_free(md_gpu_stream_t stream, md_gpu_addr_t addr);

/* ---- Copies ------------------------------------------------------------------
   The direction is in the name and in the types; nothing is inferred from
   address values. To read results on the host, copy into MD_GPU_MEM_HOST_READ
   memory and read its `.cpu` pointer once the copy's sync has completed --
   there is deliberately no copy into an arbitrary host pointer, which would
   need hidden staging and a write to host memory at an unspecified later
   poll. */

/* Device to device. Both ranges must lie within live allocations. */
bool md_gpu_copy(md_gpu_stream_t stream, md_gpu_addr_t dst, md_gpu_addr_t src, size_t size);

/* Fill `size` bytes with a repeating byte value. Unaligned heads and tails are
   handled. */
bool md_gpu_memset(md_gpu_stream_t stream, md_gpu_addr_t dst, uint8_t value, size_t size);

/* Host to device. `src` is consumed before this returns (staged when needed),
   so it may be reused immediately. */
bool md_gpu_upload(md_gpu_stream_t stream, md_gpu_addr_t dst, const void* src, size_t size);

/* Zero-copy upload: reserve `size` bytes and build the payload in place,
   avoiding an intermediate buffer plus memcpy. Returns a host pointer that
   lands directly in `dst` when that is safe, and in staging otherwise.

       float* p = md_gpu_upload_begin(s, coeff, n * sizeof(float));
       if (!p) return false;
       pack_coefficients(p, ...);
       md_gpu_upload_end(s);

   Returns NULL on failure, in which case upload_end must not be called.
   At most one upload may be open per stream. */
void* md_gpu_upload_begin(md_gpu_stream_t stream, md_gpu_addr_t dst, size_t size);
bool  md_gpu_upload_end(md_gpu_stream_t stream);

/* =========================================================================
   Textures
   ========================================================================= */

typedef enum md_gpu_tex_type_t {
    MD_GPU_TEX_TYPE_INVALID = 0,   /* a zero-initialised desc is an error */
    MD_GPU_TEX_2D,
    MD_GPU_TEX_2D_ARRAY,
    MD_GPU_TEX_3D,                 /* stays 3D even with depth 1 */
} md_gpu_tex_type_t;

typedef enum md_gpu_format_t {
    MD_GPU_FORMAT_INVALID = 0,
    /* 8-bit normalised */
    MD_GPU_FORMAT_R8_UNORM,
    MD_GPU_FORMAT_RG8_UNORM,
    MD_GPU_FORMAT_RGBA8_UNORM,
    MD_GPU_FORMAT_RGBA8_SRGB,
    MD_GPU_FORMAT_BGRA8_UNORM,
    MD_GPU_FORMAT_BGRA8_SRGB,
    /* 16-bit float */
    MD_GPU_FORMAT_R16_FLOAT,
    MD_GPU_FORMAT_RG16_FLOAT,
    MD_GPU_FORMAT_RGBA16_FLOAT,
    /* 32-bit */
    MD_GPU_FORMAT_R32_FLOAT,
    MD_GPU_FORMAT_RG32_FLOAT,
    MD_GPU_FORMAT_RGBA32_FLOAT,
    MD_GPU_FORMAT_R32_UINT,
    MD_GPU_FORMAT_RG32_UINT,
    MD_GPU_FORMAT_RGBA32_UINT,
    /* packed */
    MD_GPU_FORMAT_RG11B10_FLOAT,
    MD_GPU_FORMAT_RGB10A2_UNORM,
    /* depth */
    MD_GPU_FORMAT_D32_FLOAT,
    MD_GPU_FORMAT_D32_FLOAT_S8_UINT,
    MD_GPU_FORMAT_COUNT,
} md_gpu_format_t;

/* Bytes per texel as laid out in buffers by the texture copy calls. For
   MD_GPU_FORMAT_D32_FLOAT_S8_UINT that is the depth plane only (4 bytes). */
uint32_t md_gpu_format_texel_size(md_gpu_format_t format);

typedef uint32_t md_gpu_tex_usage_t;
enum {
    MD_GPU_TEX_STORAGE       = 1u << 0,  /* shader read/write, random access */
    MD_GPU_TEX_SAMPLED       = 1u << 1,  /* shader sampled read              */
    MD_GPU_TEX_RENDER_TARGET = 1u << 2,  /* colour or depth attachment (by
                                            format); takes no heap slot      */
};

typedef struct md_gpu_texture_desc_t {
    md_gpu_tex_type_t  type;
    md_gpu_format_t    format;
    md_gpu_tex_usage_t usage;            /* any non-empty combination. A
                                            format/usage pair the device
                                            cannot do fails at creation,
                                            naming both                       */
    uint32_t           width;
    uint32_t           height;
    uint32_t           depth_or_layers;  /* 3D: depth. 2D_ARRAY: layers.
                                            2D: must be 0 or 1                */
    uint32_t           mip_levels;       /* 0 = 1                            */
    const char*        label;
} md_gpu_texture_desc_t;

/* Stream-ordered creation, exactly like md_gpu_malloc: the texture is usable
   by work issued into `stream` after this call (other streams join with
   md_gpu_stream_wait). Never blocks. `pool` must be an MD_GPU_MEM_DEVICE pool;
   it owns the texture's lifetime. Returns NULL on failure. */
md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t stream, md_gpu_pool_t pool,
                                       const md_gpu_texture_desc_t* desc);

/* Deferred and non-blocking: the texture and its handles are released once
   every stream has completed the work issued before this call. */
void md_gpu_texture_destroy(md_gpu_texture_t tex);

/* The description the texture was created with, normalised (mip_levels >= 1,
   depth_or_layers >= 1). Valid for the texture's lifetime. */
const md_gpu_texture_desc_t* md_gpu_texture_desc(md_gpu_texture_t tex);

/* Shader handles, created with the texture and valid for its lifetime. Null if
   the texture lacks the usage or `mip` is out of range. A storage handle
   addresses a single mip level; the sampled handle covers all of them. */
md_gpu_storage_tex_t md_gpu_texture_storage(md_gpu_texture_t tex, uint32_t mip);
md_gpu_sampled_tex_t md_gpu_texture_sampled(md_gpu_texture_t tex);

/* ---- Samplers -----------------------------------------------------------------
   Samplers are immutable values: the same desc returns the same handle, the
   device owns them, and there is nothing to destroy. A zero-initialised desc
   is nearest filtering with clamp-to-edge addressing. */

typedef enum md_gpu_filter_t {
    MD_GPU_FILTER_NEAREST,
    MD_GPU_FILTER_LINEAR,
} md_gpu_filter_t;

typedef enum md_gpu_address_mode_t {
    MD_GPU_ADDRESS_CLAMP_TO_EDGE,
    MD_GPU_ADDRESS_REPEAT,
    MD_GPU_ADDRESS_MIRRORED_REPEAT,
} md_gpu_address_mode_t;

typedef struct md_gpu_sampler_desc_t {
    md_gpu_filter_t       min_filter, mag_filter, mip_filter;
    md_gpu_address_mode_t address_u, address_v, address_w;
} md_gpu_sampler_desc_t;

md_gpu_sampler_t md_gpu_sampler(md_gpu_device_t device, const md_gpu_sampler_desc_t* desc);

/* ---- Texture copies -----------------------------------------------------------
   A region is in texels of one mip level. For MD_GPU_TEX_2D_ARRAY the z axis
   addresses layers. A zero extent component means "to the end along that
   axis", so a zero-initialised region (or NULL) is the whole of mip 0. Buffer
   data is tightly packed, md_gpu_format_texel_size() bytes per texel; byte
   counts are derived from the region and checked against the allocation. */

typedef struct md_gpu_tex_region_t {
    uint32_t offset[3];
    uint32_t extent[3];
    uint32_t mip;
} md_gpu_tex_region_t;

/* Bytes a copy of `region` of `tex` moves (0 if the region is invalid). */
size_t md_gpu_texture_region_size(md_gpu_texture_t tex, const md_gpu_tex_region_t* region);

bool md_gpu_copy_to_texture(md_gpu_stream_t stream, md_gpu_texture_t dst,
                            const md_gpu_tex_region_t* region, md_gpu_addr_t src);
bool md_gpu_copy_from_texture(md_gpu_stream_t stream, md_gpu_addr_t dst,
                              md_gpu_texture_t src, const md_gpu_tex_region_t* region);

/* Host to texture. `src` is consumed before return; `size` must equal the
   region's byte size. */
bool md_gpu_upload_texture(md_gpu_stream_t stream, md_gpu_texture_t dst,
                           const md_gpu_tex_region_t* region, const void* src, size_t size);

/* =========================================================================
   Kernels
   ========================================================================= */

typedef struct md_gpu_kernel_desc_t {
    /* SPIR-V on Vulkan. On Metal, either a compiled metallib or Metal Shading
       Language source text -- the backend identifies which from the bytes. */
    const void* code;
    size_t      code_size;
    const char* entry_point;   /* NULL = "main" */
    const char* label;

    /* Threads per group, i.e. [numthreads]. Required on every backend: zero
       is an error. On Vulkan it is also checked against the SPIR-V. */
    uint32_t    group_size[3];

    /* sizeof the argument struct. A launch passing a different size fails.
       0 = unchecked. */
    uint32_t    args_size;
} md_gpu_kernel_desc_t;

/* Normally fed straight from the generated descriptor:

       md_gpu_kernel_desc_t d = md_shader_topo_critical_points_main_kernel();
       md_gpu_kernel_t k = md_gpu_kernel_create(dev, &d);                  */
md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t device, const md_gpu_kernel_desc_t* desc);

/* Deferred and non-blocking, like md_gpu_texture_destroy. */
void md_gpu_kernel_destroy(md_gpu_kernel_t kernel);

typedef struct md_gpu_kernel_info_t {
    uint32_t group_size[3];
    uint32_t args_size;
    uint32_t max_threads_per_group;
    uint32_t preferred_group_multiple;
} md_gpu_kernel_info_t;

bool md_gpu_kernel_info(md_gpu_kernel_t kernel, md_gpu_kernel_info_t* out_info);

/* The grid of groups covering nx * ny * nz threads with this kernel's group
   size: ceil(n / group_size) per axis. */
md_gpu_grid_t md_gpu_grid_for(md_gpu_kernel_t kernel, uint32_t nx, uint32_t ny, uint32_t nz);

/* kernel<<<grid, block, 0, stream>>>(args).

   `args` is copied immediately; the caller may reuse or free the memory as
   soon as this returns. There is no practical size limit. An empty grid is a
   no-op. Kernels cannot be launched into an MD_GPU_STREAM_TRANSFER stream. */
bool md_gpu_launch(md_gpu_stream_t stream, md_gpu_kernel_t kernel, md_gpu_grid_t grid,
                   const void* args, size_t args_size);

/* The grid is read from device memory: 3 consecutive uint32 at `grid`. */
bool md_gpu_launch_indirect(md_gpu_stream_t stream, md_gpu_kernel_t kernel, md_gpu_addr_t grid,
                            const void* args, size_t args_size);

/* Pass an argument struct by value; its size comes from sizeof. */
#define MD_GPU_LAUNCH(stream, kernel, grid, args) \
    md_gpu_launch((stream), (kernel), (grid), &(args), sizeof(args))

/* Built-in helper: reads one uint32 thread count at `count` and writes the
   indirect grid covering it with `kernel`'s group size -- { ceil(count /
   group_size.x), 1, 1 } -- to `out_grid` as 3 uint32. Turns a device-side
   count into an indirect launch without a readback. */
bool md_gpu_make_grid(md_gpu_stream_t stream, md_gpu_addr_t out_grid,
                      md_gpu_addr_t count, md_gpu_kernel_t kernel);

/* =========================================================================
   Host-side ordering
   ========================================================================= */

typedef void (*md_gpu_host_fn)(void* user);

/* cudaLaunchHostFunc: `fn` runs after all work issued into `stream` before
   this call has completed. It runs inside md_gpu_device_poll(), on whichever
   thread calls that -- so the caller picks the thread and there is no locking
   in user code. */
bool md_gpu_launch_host_fn(md_gpu_stream_t stream, md_gpu_host_fn fn, void* user);

/* Same, keyed to an explicit sync point. A none sync fires on the next poll. */
bool md_gpu_sync_on_complete(md_gpu_device_t device, md_gpu_sync_t sync, md_gpu_host_fn fn, void* user);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

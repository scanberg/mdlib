/*
md_gpu.h

A thin layer over Vulkan and Metal (and, as an alternative compute-only build,
CUDA). One backend per build. The compute side is CUDA-shaped, and the raster
side is the same model applied to draws, after Aaltonen's "No Graphics API":
pointers, one root argument pointer, a bindless heap, coarse stage barriers,
no vertex formats and no resource state.

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
    textures and kernels, or freeing memory, never waits for work in flight -- a compute
    job spanning many frames stalls nothing but its own stream.

Correspondence with CUDA
------------------------
    cudaStreamCreate            md_gpu_stream_create
    cudaMallocAsync             md_gpu_malloc
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
    MD_SHADER_ARGS(Args);                   // from md_gpu.slang

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
  * Allocation (md_gpu_malloc / md_gpu_free), texture, sampler and kernel
    creation/destruction are thread-safe. Temp scopes belong to their stream
    and follow the stream's one-thread rule.
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
typedef struct md_gpu_texture* md_gpu_texture_t;   /* identity; host side only */
typedef struct md_gpu_kernel*  md_gpu_kernel_t;

/* A GPU virtual address. Byte arithmetic is valid (`base + offsetof(T, f)`),
   and it is the type of every pointer field in a C argument-struct mirror, so
   no casts are needed there. Zero is null. Not dereferenceable on the host. */
typedef uint64_t md_gpu_addr_t;

/* An allocation. `cpu` is non-NULL only for MD_GPU_MEM_HOST_* memory, and
   then addresses the same bytes as `gpu`. */
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

    /* Bytes of empty heap chunks each memory kind keeps for reuse rather than
       returning them to the driver. 0 selects 256 MiB. */
    uint64_t heap_cache_limit;

    const char* label;
} md_gpu_device_desc_t;

typedef struct md_gpu_device_info_t {
    /* False implies UMA / integrated. An allocation hint only. */
    bool     is_discrete;
    uint32_t max_threads_per_group;
    uint32_t preferred_group_multiple;   /* warp / SIMD width */
    bool     supports_graphics;          /* GRAPHICS streams and rendering */
    bool     supports_present;           /* surfaces can be created        */
    char     name[256];
} md_gpu_device_info_t;

/* `desc` may be NULL for defaults. Returns NULL on failure; see
   md_gpu_last_error(). */
md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc);

/* Waits for every stream to go idle, then destroys the device and everything
   created from it: streams, allocations, textures, samplers, kernels. */
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
    MD_GPU_STREAM_GRAPHICS,   /* render passes and presentation, plus all a
                                 COMPUTE stream does                         */
} md_gpu_stream_kind_t;

/* Why a third kind: Vulkan puts COMPUTE streams on an async compute engine
   where one exists, and that engine cannot rasterise. Keeping the kinds apart
   keeps long-running compute off the queue that presents, while a GRAPHICS
   stream can still launch per-frame kernels without a cross-queue hop. Metal
   queues are universal, so there the kind only gates validation. GRAPHICS
   streams exist when md_gpu_device_info_t.supports_graphics is set. */

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
   layouts. Both edges of an explicit region are ordered: its first operation
   runs after everything before it, and switching back to IMPLICIT orders the
   next operation after everything in it, so an explicit region can never leak
   unordered work past either end.

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
    MD_GPU_STAGE_INDIRECT = 1u << 2,   /* consumer side: indirect grids and draws */
    MD_GPU_STAGE_VERTEX     = 1u << 3, /* index fetch and vertex shading          */
    MD_GPU_STAGE_FRAGMENT   = 1u << 4, /* fragment shading                        */
    MD_GPU_STAGE_ATTACHMENT = 1u << 5, /* depth test, colour output, attachment
                                          load and store                          */
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
    MD_GPU_MEM_KIND_COUNT,
} md_gpu_mem_kind_t;

/* Two lifetimes, two pairs of calls:

     persistent   md_gpu_malloc / md_gpu_free. Freed one at a time; drawn from
                  a device-wide heap per memory kind.
     temporary    md_gpu_temp_alloc between md_gpu_temp_begin and
                  md_gpu_temp_end. Bump-allocated from the stream's own arena
                  and released together when the scope ends.

   Both are stream-ordered, like everything else. Memory is usable by work
   issued into the stream after it was allocated, and a free (or scope end)
   takes effect when the GPU reaches that point in the stream. The CPU never
   waits: memory freed at a point the GPU has not yet reached is reused only
   when that is safe (see md_gpu_free).

   Every allocation is 256-byte aligned. Allocations made before the device
   is destroyed are released with it. */

/* cudaMallocAsync. `.gpu == 0` on failure; `.cpu` is non-NULL for the HOST_*
   kinds. The heap carves allocations out of large chunks, so small
   allocations are cheap and do not count against the driver's allocation
   limit. */
md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t stream, md_gpu_mem_kind_t kind, size_t size);

/* cudaFreeAsync. The memory is released at this point in `stream`. Later
   MD_GPU_MEM_DEVICE allocations on the same stream may reuse it at once,
   because stream order already puts their work after the free (in EXPLICIT
   mode md_gpu inserts a barrier when it does this). Host-visible memory, which
   the CPU writes the moment it is handed out, and every allocation on other
   streams reuse it only once the GPU has passed the free. `addr` must be the
   start of a live allocation. A zero address is a no-op; `stream` is
   required. */
void md_gpu_free(md_gpu_stream_t stream, md_gpu_addr_t addr);

/* ---- Temporary memory -------------------------------------------------------
   The GPU form of a CPU temp arena. Allocation bumps a pointer; there is no
   per-allocation free. The free is md_gpu_temp_end, and like md_gpu_free it
   is recorded on the stream: everything allocated since the matching begin
   is reclaimed once the GPU passes that point. The CPU never waits for it.
   Allocations after the end go into other memory until then, so double or
   triple buffering falls out without the caller counting frames.

       md_gpu_temp_t frame = md_gpu_temp_begin(gfx);
       md_gpu_mem_t v = md_gpu_temp_alloc(gfx, MD_GPU_MEM_HOST_WRITE, bytes);
       memcpy(v.cpu, verts, bytes);
       ... work that reads v.gpu ...
       md_gpu_temp_end(gfx, frame);

   Scopes nest, and must end in reverse order of beginning, exactly as CPU temp
   arenas do. That lets a library take temp memory inside its caller's scope
   without touching the caller's allocations. md_gpu_temp_alloc outside any
   scope fails.

   Rules:
     * Temp memory belongs to its stream. Another stream that reads it must be
       joined before the scope ends, e.g.
       md_gpu_stream_wait(s, md_gpu_stream_record(other)), as for md_gpu_free.
     * Only MD_GPU_MEM_DEVICE and MD_GPU_MEM_HOST_WRITE. Readback memory is
       read by host callbacks, which run at md_gpu_device_poll -- possibly
       after a later scope has already reused the memory -- so readbacks use
       md_gpu_malloc / md_gpu_free.
     * Scratch whose lifetime is not nested in anything on the stream (say, a
       compute job spanning many frames) belongs in md_gpu_malloc.
     * Bounds are checked per arena chunk, not per temp allocation. */

typedef struct md_gpu_temp_t {
    md_gpu_stream_t stream;
    uint32_t        depth;
} md_gpu_temp_t;

md_gpu_temp_t md_gpu_temp_begin(md_gpu_stream_t stream);
md_gpu_mem_t  md_gpu_temp_alloc(md_gpu_stream_t stream, md_gpu_mem_kind_t kind, size_t size);
void          md_gpu_temp_end(md_gpu_stream_t stream, md_gpu_temp_t scope);

/* ---- Statistics -------------------------------------------------------------- */

typedef struct md_gpu_memory_stats_t {
    uint64_t bytes_in_use;       /* live md_gpu_malloc allocations              */
    uint64_t bytes_peak_in_use;  /* high-water mark of bytes_in_use             */
    uint64_t bytes_reserved;     /* heap chunks held from the driver, in use
                                    or not                                      */
    uint64_t bytes_temp;         /* held by temp arenas, over all streams       */
    uint64_t bytes_textures;     /* textures; MD_GPU_MEM_DEVICE only            */
    uint32_t allocations;        /* live md_gpu_malloc allocations              */
    uint32_t chunks;             /* heap chunks                                 */
} md_gpu_memory_stats_t;

bool md_gpu_memory_stats(md_gpu_device_t device, md_gpu_mem_kind_t kind, md_gpu_memory_stats_t* out_stats);

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
                                            format); 2D and 2D_ARRAY only.
                                            Combines with SAMPLED, and with
                                            STORAGE for colour formats the
                                            device can write                 */
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
   md_gpu_stream_wait). Never blocks. Textures are device-local. Returns NULL
   on failure. */
md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t stream, const md_gpu_texture_desc_t* desc);

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

/* Texture to texture: same format, same extent, no scaling (a scaled copy is
   a draw). The extent is taken from `src_region`; `dst_region` supplies the
   destination offset, mip and layer. */
bool md_gpu_copy_texture(md_gpu_stream_t stream,
                         md_gpu_texture_t dst, const md_gpu_tex_region_t* dst_region,
                         md_gpu_texture_t src, const md_gpu_tex_region_t* src_region);

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
   Rendering
   =========================================================================

   Rasterisation on the same terms as compute:

     * A draw works like a launch. It names its pipeline and takes a copied
       argument struct that both shader stages read through the same root
       pointer (MD_ARGS). There are no vertex buffers or input layouts: the
       vertex shader reads what it needs through pointers in the struct.
       Index and indirect-command buffers are GPU addresses.

     * A pipeline holds what Vulkan and Metal both compile into it: shaders,
       topology, attachment formats, blend and write masks. The state both
       APIs make dynamic -- depth test and write, culling, winding, depth
       bias, blend constant, viewport, scissor -- is set inside the pass.

     * A render pass is a scope on a GRAPHICS stream, not an object:
       md_gpu_render_begin names the attachments and their load and store
       actions, md_gpu_render_end closes it. Attachments are ordinary
       textures with MD_GPU_TEX_RENDER_TARGET usage. There are no image
       layouts anywhere in the API.

     * Ordering is the stream's. In IMPLICIT mode a pass is ordered after
       everything before it and before everything after it. In EXPLICIT mode
       the caller places md_gpu_barrier with the raster stages.

     * Presentation goes through a surface made from native window handles.
       Acquire returns a texture; present runs on the stream.

   CONVENTIONS

     Clip space as in Metal and Direct3D: +Y up, depth 0..1. Framebuffer and
     texture space: origin at the top-left texel. Viewport and scissor are
     in framebuffer pixels. The Vulkan backend flips its viewport so that
     the same shader and matrices produce the same image on every backend,
     and a rendered texture reads back the right way up.

     Front faces are counter-clockwise in clip space unless the draw state
     says otherwise, on every backend.

     In a direct draw, SV_VertexID runs 0..vertex_count-1 (for an indexed
     draw, it is the index value) and SV_InstanceID runs
     0..instance_count-1. Direct draws have no base vertex or base instance:
     with vertex pulling those are pointer arithmetic in the argument struct.
     Indirect commands do carry them, because their layout is fixed by the
     hardware, and there the backends disagree: Metal's ids include the
     bases, Slang's SPIR-V subtracts them. Shaders drawn indirectly read ids
     through md_draw_vertex() / md_draw_instance() from md_gpu.slang, which
     return the draw-relative id everywhere. SV_StartInstanceLocation is
     portable and is the way to give each indirect draw its own record: a
     culling kernel writes draw i with first_instance = i.

     The all-ones index (0xFFFF / 0xFFFFFFFF) restarts TRIANGLE_STRIP and
     LINE_STRIP and must not appear with list topologies. Lines are one
     pixel wide. A POINTS pipeline's vertex shader must write the point
     size (an SV_PointSize output).

     Per-frame data -- constants too large for the argument struct, vertices
     built on the CPU -- belongs in a temp scope on the graphics stream (see
     Memory): begin it at the start of the frame, end it after present. */

#define MD_GPU_MAX_COLOR_TARGETS 8

/* ---- Pipelines --------------------------------------------------------------
   Creation is synchronous and may take a while (driver compilation), so
   create pipelines up front and keep them. Destruction is deferred and
   non-blocking, like kernels. Blend is part of the pipeline because Metal
   compiles it into the fragment shader and Vulkan only makes it dynamic
   through an optional extension. */

typedef struct md_gpu_pipeline* md_gpu_pipeline_t;

/* One shader entry point. compile_gpu_shaders generates these as
   <namespace>_<stem>_<entry>_shader() for VERTEX and FRAGMENT entries; they
   are not written by hand. `args_size` is read from the compiled shader and
   checked against every draw. */
typedef struct md_gpu_shader_t {
    const void* code;             /* SPIR-V, or metallib / MSL source       */
    size_t      code_size;
    const char* entry_point;
    uint32_t    args_size;        /* 0: the entry point takes no arguments  */
} md_gpu_shader_t;

typedef enum md_gpu_topology_t {
    MD_GPU_TOPOLOGY_TRIANGLES = 0,
    MD_GPU_TOPOLOGY_TRIANGLE_STRIP,
    MD_GPU_TOPOLOGY_LINES,
    MD_GPU_TOPOLOGY_LINE_STRIP,
    MD_GPU_TOPOLOGY_POINTS,
} md_gpu_topology_t;

typedef enum md_gpu_blend_factor_t {
    MD_GPU_BLEND_ZERO = 0,
    MD_GPU_BLEND_ONE,
    MD_GPU_BLEND_SRC_COLOR,
    MD_GPU_BLEND_ONE_MINUS_SRC_COLOR,
    MD_GPU_BLEND_SRC_ALPHA,
    MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA,
    MD_GPU_BLEND_DST_COLOR,
    MD_GPU_BLEND_ONE_MINUS_DST_COLOR,
    MD_GPU_BLEND_DST_ALPHA,
    MD_GPU_BLEND_ONE_MINUS_DST_ALPHA,
    MD_GPU_BLEND_CONSTANT,             /* md_gpu_draw_state_t.blend_constant */
    MD_GPU_BLEND_ONE_MINUS_CONSTANT,
    MD_GPU_BLEND_SRC_ALPHA_SATURATE,
} md_gpu_blend_factor_t;

typedef enum md_gpu_blend_op_t {
    MD_GPU_BLEND_OP_ADD = 0,
    MD_GPU_BLEND_OP_SUBTRACT,
    MD_GPU_BLEND_OP_REVERSE_SUBTRACT,
    MD_GPU_BLEND_OP_MIN,
    MD_GPU_BLEND_OP_MAX,
} md_gpu_blend_op_t;

typedef uint32_t md_gpu_color_mask_t;
enum {
    MD_GPU_COLOR_R   = 1u << 0,
    MD_GPU_COLOR_G   = 1u << 1,
    MD_GPU_COLOR_B   = 1u << 2,
    MD_GPU_COLOR_A   = 1u << 3,
    MD_GPU_COLOR_ALL = 0xFu,
};

/* Zero-initialised: blending off. With `enable`, the factors are used as
   given. Common setups:

       premultiplied  src ONE,       dst ONE_MINUS_SRC_ALPHA  (colour and alpha)
       straight       src SRC_ALPHA, dst ONE_MINUS_SRC_ALPHA  (colour),
                      src ONE,       dst ONE_MINUS_SRC_ALPHA  (alpha)
       additive       src ONE,       dst ONE

   Blending an integer target is an error at pipeline creation. */
typedef struct md_gpu_blend_t {
    bool                  enable;
    md_gpu_blend_factor_t src_color, dst_color;
    md_gpu_blend_op_t     color_op;
    md_gpu_blend_factor_t src_alpha, dst_alpha;
    md_gpu_blend_op_t     alpha_op;
} md_gpu_blend_t;

typedef struct md_gpu_color_target_t {
    md_gpu_format_t     format;
    md_gpu_blend_t      blend;
    md_gpu_color_mask_t write_disable;   /* channels NOT written; 0 = all written */
} md_gpu_color_target_t;

typedef struct md_gpu_pipeline_desc_t {
    md_gpu_shader_t       vertex;
    md_gpu_shader_t       fragment;      /* code == NULL: no fragment stage
                                            (depth-only passes)                */
    md_gpu_topology_t     topology;
    md_gpu_color_target_t color[MD_GPU_MAX_COLOR_TARGETS];
    uint32_t              color_count;
    md_gpu_format_t       depth_format;  /* MD_GPU_FORMAT_INVALID: no depth
                                            attachment                         */
    const char*           label;
} md_gpu_pipeline_desc_t;

/* Fails if the stages disagree on args_size (they share one argument
   struct), or if a colour format cannot be rendered to or blended as asked. */
md_gpu_pipeline_t md_gpu_pipeline_create(md_gpu_device_t device, const md_gpu_pipeline_desc_t* desc);
void              md_gpu_pipeline_destroy(md_gpu_pipeline_t pipeline);

/* ---- Render passes ------------------------------------------------------------
   A pass is a scope on a GRAPHICS stream:

       md_gpu_render_begin(gfx, &pass);
           ... md_gpu_set_* and md_gpu_draw* ...
       md_gpu_render_end(gfx);

   Inside it, anything else that records GPU work on that stream fails:
   launches, copies, uploads, barriers, texture creation, record, flush,
   sync, md_gpu_stream_wait, surface calls, and a nested begin. md_gpu_malloc, md_gpu_free and
   temp scopes stay legal because they record nothing.

   Every attachment must be the same size at the mip it names; that is the
   render area. At begin, viewport and scissor cover the render area and the
   draw state is the zero-initialised md_gpu_draw_state_t.

   Draws in a pass are ordered for their attachments (rasterisation order:
   depth testing and blending see earlier draws). A draw's storage writes
   are not visible to later draws in the same pass; end the pass and begin
   another with MD_GPU_LOAD for that. Sampling a texture that the same pass
   has bound as an attachment is undefined. */

typedef enum md_gpu_load_t {
    MD_GPU_LOAD = 0,           /* keep the existing contents                   */
    MD_GPU_LOAD_CLEAR,         /* fill with the attachment's clear value       */
    MD_GPU_LOAD_DONT_CARE,     /* undefined; cheapest, especially on tilers    */
} md_gpu_load_t;

typedef enum md_gpu_store_t {
    MD_GPU_STORE = 0,          /* keep what the pass wrote                     */
    MD_GPU_STORE_DISCARD,      /* undefined afterwards: transient depth, say   */
} md_gpu_store_t;

/* Read by the attachment's format: f32 for float, unorm and sRGB, u32 for
   the UINT formats. */
typedef union md_gpu_clear_color_t {
    float    f32[4];
    uint32_t u32[4];
} md_gpu_clear_color_t;

typedef struct md_gpu_color_attachment_t {
    md_gpu_texture_t     texture;
    uint32_t             mip;
    uint32_t             layer;       /* MD_GPU_TEX_2D_ARRAY only             */
    md_gpu_load_t        load;
    md_gpu_store_t       store;
    md_gpu_clear_color_t clear;
} md_gpu_color_attachment_t;

typedef struct md_gpu_depth_attachment_t {
    md_gpu_texture_t texture;         /* NULL: no depth attachment            */
    uint32_t         mip;
    uint32_t         layer;
    md_gpu_load_t    load;
    md_gpu_store_t   store;
    float            clear_depth;
} md_gpu_depth_attachment_t;

typedef struct md_gpu_render_desc_t {
    md_gpu_color_attachment_t color[MD_GPU_MAX_COLOR_TARGETS];
    uint32_t                  color_count;
    md_gpu_depth_attachment_t depth;
    const char*               label;  /* debug marker around the pass         */
} md_gpu_render_desc_t;

bool md_gpu_render_begin(md_gpu_stream_t stream, const md_gpu_render_desc_t* desc);
bool md_gpu_render_end(md_gpu_stream_t stream);

/* ---- Dynamic state ----------------------------------------------------------
   Set inside a pass; kept until changed or the pass ends. Zero means off in
   every field, and a zero-initialised struct (or NULL) is the state at
   md_gpu_render_begin. The backend compares against what is set and emits
   only the difference. */

typedef enum md_gpu_compare_t {
    MD_GPU_COMPARE_ALWAYS = 0,        /* = no depth test */
    MD_GPU_COMPARE_NEVER,
    MD_GPU_COMPARE_LESS,
    MD_GPU_COMPARE_LESS_EQUAL,
    MD_GPU_COMPARE_EQUAL,
    MD_GPU_COMPARE_NOT_EQUAL,
    MD_GPU_COMPARE_GREATER_EQUAL,
    MD_GPU_COMPARE_GREATER,
} md_gpu_compare_t;

typedef enum md_gpu_cull_t {
    MD_GPU_CULL_NONE = 0,
    MD_GPU_CULL_BACK,
    MD_GPU_CULL_FRONT,
} md_gpu_cull_t;

typedef struct md_gpu_draw_state_t {
    /* Depth. depth_write with ALWAYS writes unconditionally. Depth written
       by the fragment shader (SV_Depth*) replaces the interpolated value for
       both test and write; impostors should prefer SV_DepthGreaterEqual
       (SV_DepthLessEqual under reverse-Z), which keeps early rejection. */
    md_gpu_compare_t depth_compare;
    bool             depth_write;

    md_gpu_cull_t    cull;
    bool             front_clockwise;     /* default: counter-clockwise      */

    float            depth_bias;          /* constant, in depth units        */
    float            depth_bias_slope;    /* times the polygon's depth slope */
    float            depth_bias_clamp;    /* 0 = unclamped                   */

    float            blend_constant[4];   /* for MD_GPU_BLEND_CONSTANT       */
} md_gpu_draw_state_t;

typedef struct md_gpu_viewport_t {
    float x, y, width, height;            /* framebuffer pixels              */
    float min_depth, max_depth;           /* both 0: 0..1                    */
} md_gpu_viewport_t;

typedef struct md_gpu_rect_t {
    uint32_t x, y, width, height;         /* framebuffer pixels              */
} md_gpu_rect_t;

/* Errors outside a pass. NULL restores the value the pass began with. */
void md_gpu_set_draw_state(md_gpu_stream_t stream, const md_gpu_draw_state_t* state);
void md_gpu_set_viewport(md_gpu_stream_t stream, const md_gpu_viewport_t* viewport);
void md_gpu_set_scissor(md_gpu_stream_t stream, const md_gpu_rect_t* scissor);

/* ---- Draws --------------------------------------------------------------------
   Shaped like md_gpu_launch: a pipeline, an amount of work, and an argument
   struct that is copied at the call and read by both stages. A zero count is
   a no-op. Every draw checks that it is inside a pass on this stream, that
   the pipeline's attachment formats equal the pass's, in order, and the
   argument size -- a format mismatch is undefined behaviour in Vulkan and
   Metal alike, so the check stays in release builds. */

typedef enum md_gpu_index_type_t {
    MD_GPU_INDEX_U32 = 0,
    MD_GPU_INDEX_U16,
} md_gpu_index_type_t;

bool md_gpu_draw(md_gpu_stream_t stream, md_gpu_pipeline_t pipeline,
                 uint32_t vertex_count, uint32_t instance_count,
                 const void* args, size_t args_size);

/* `indices` addresses the first index and must be aligned to the index size;
   offsetting into an index buffer is address arithmetic. */
bool md_gpu_draw_indexed(md_gpu_stream_t stream, md_gpu_pipeline_t pipeline,
                         md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                         uint32_t index_count, uint32_t instance_count,
                         const void* args, size_t args_size);

/* Indirect commands, in the layout Vulkan and Metal both read directly, so a
   kernel writes them and nothing converts them. Tightly packed arrays,
   4-byte aligned. A culling kernel rejects a draw by writing an
   instance_count of 0. */
typedef struct md_gpu_draw_cmd_t {
    uint32_t vertex_count;
    uint32_t instance_count;
    uint32_t first_vertex;
    uint32_t first_instance;
} md_gpu_draw_cmd_t;

typedef struct md_gpu_draw_indexed_cmd_t {
    uint32_t index_count;
    uint32_t instance_count;
    uint32_t first_index;          /* in indices, from `indices` */
    int32_t  vertex_offset;
    uint32_t first_instance;
} md_gpu_draw_indexed_cmd_t;

/* Draws `count` consecutive commands starting at `cmds`. Producers are
   ordered by the stream, or by md_gpu_barrier(s, ..., MD_GPU_STAGE_INDIRECT). */
bool md_gpu_draw_indirect(md_gpu_stream_t stream, md_gpu_pipeline_t pipeline,
                          md_gpu_addr_t cmds, uint32_t count,
                          const void* args, size_t args_size);

bool md_gpu_draw_indexed_indirect(md_gpu_stream_t stream, md_gpu_pipeline_t pipeline,
                                  md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                                  md_gpu_addr_t cmds, uint32_t count,
                                  const void* args, size_t args_size);

/* Pass an argument struct by value; its size comes from sizeof. */
#define MD_GPU_DRAW(stream, pipeline, vertex_count, instance_count, args) \
    md_gpu_draw((stream), (pipeline), (vertex_count), (instance_count), &(args), sizeof(args))

#define MD_GPU_DRAW_INDEXED(stream, pipeline, indices, type, index_count, instance_count, args) \
    md_gpu_draw_indexed((stream), (pipeline), (indices), (type), (index_count), (instance_count), &(args), sizeof(args))

/* ---- Presentation -------------------------------------------------------------
   A surface is the swapchain of one window, made from native handles so that
   md_gpu does not depend on a windowing library. From GLFW:
   glfwGetWin32Window, glfwGetCocoaWindow + contentView, glfwGetX11Window +
   glfwGetX11Display, glfwGetWaylandWindow + glfwGetWaylandDisplay.

       md_gpu_texture_t back = md_gpu_surface_acquire(gfx, surface);
       if (back) {
           ... passes that render into `back` ...
           md_gpu_surface_present(gfx, surface);
       }

   Acquire waits for the presentation engine to hand back an image; with
   VSYNC that is what paces the frame loop. Apart from the sync calls and
   teardown it is the only call in md_gpu that blocks. It returns NULL when
   the drawable size is zero (a minimised window: skip the frame) or on
   failure (md_gpu_last_error says which). Out-of-date swapchains are rebuilt
   inside acquire at the size last given to md_gpu_surface_resize.

   The acquired texture is 2D, in the surface's format, with RENDER_TARGET
   usage plus what the desc asked for. It is valid from acquire to present,
   has undefined contents at acquire (load CLEAR or DONT_CARE), and is never
   destroyed by the caller. Other streams that touch it must be joined with
   md_gpu_stream_record / md_gpu_stream_wait.

   Present submits the stream, as md_gpu_stream_record does, and queues the
   image for display once that work completes. There is one present per
   acquire, on the stream that acquired.

   Frames in flight need no API: the argument arena and temp scopes recycle
   themselves, and anything else is bounded with syncs:

       md_gpu_sync_wait(frame_done[frame % 2]);         // at most 2 ahead
       ... frame ...
       md_gpu_surface_present(gfx, surface);
       frame_done[frame % 2] = md_gpu_stream_record(gfx); */

typedef struct md_gpu_surface* md_gpu_surface_t;

typedef enum md_gpu_window_system_t {
    MD_GPU_WINDOW_INVALID = 0,
    MD_GPU_WINDOW_WIN32,      /* window: HWND.  display: HINSTANCE or NULL            */
    MD_GPU_WINDOW_COCOA,      /* window: NSView* or CAMetalLayer*                     */
    MD_GPU_WINDOW_X11,        /* window: Window (the XID, cast).  display: Display*   */
    MD_GPU_WINDOW_WAYLAND,    /* window: wl_surface*.  display: wl_display*           */
    MD_GPU_WINDOW_HEADLESS,   /* no window; images go nowhere. For tests and CI      */
} md_gpu_window_system_t;

typedef enum md_gpu_present_mode_t {
    MD_GPU_PRESENT_VSYNC = 0, /* always available                                     */
    MD_GPU_PRESENT_MAILBOX,   /* low latency, no tearing. Falls back to VSYNC         */
    MD_GPU_PRESENT_IMMEDIATE, /* no wait, may tear. Falls back to MAILBOX, then VSYNC */
} md_gpu_present_mode_t;

typedef struct md_gpu_surface_desc_t {
    md_gpu_window_system_t system;
    void*                  window;
    void*                  display;
    uint32_t               width, height;  /* drawable size in pixels            */
    md_gpu_format_t        format;         /* INVALID = BGRA8_UNORM; BGRA8_SRGB and
                                              RGBA16_FLOAT where supported. No
                                              silent substitute                  */
    md_gpu_tex_usage_t     usage;          /* beyond RENDER_TARGET: SAMPLED and/or
                                              STORAGE (compose with a kernel).
                                              Fails if unsupported               */
    md_gpu_present_mode_t  present_mode;
    const char*            label;
} md_gpu_surface_desc_t;

/* COCOA with an NSView must be called on the main thread (it attaches a
   CAMetalLayer); passing a CAMetalLayer lifts that. */
md_gpu_surface_t md_gpu_surface_create(md_gpu_device_t device, const md_gpu_surface_desc_t* desc);

/* Deferred and non-blocking. The window must outlive the surface's last present. */
void md_gpu_surface_destroy(md_gpu_surface_t surface);

/* From the window's framebuffer-size callback; applied at the next acquire. */
void md_gpu_surface_resize(md_gpu_surface_t surface, uint32_t width, uint32_t height);

md_gpu_texture_t md_gpu_surface_acquire(md_gpu_stream_t stream, md_gpu_surface_t surface);
bool             md_gpu_surface_present(md_gpu_stream_t stream, md_gpu_surface_t surface);

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

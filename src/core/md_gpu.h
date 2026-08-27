/*
md_gpu.h

A CUDA-shaped compute API over Vulkan and Metal.

The model, in full:

    Work issued into a stream executes in issue order, and every operation
    observes all writes made by the operations before it in that stream.
    Work in different streams is unordered unless joined by a sync point.

That is the entire dependency model. There are no barriers, no resource state,
no usage declarations and no descriptor sets in this API, because there is
nothing for the caller to get wrong. Concurrency is expressed by using more
streams, exactly as in CUDA.

Correspondence with CUDA
------------------------
    cudaStreamCreate            md_gpu_stream_create
    cudaMallocAsync             md_gpu_malloc
    cudaFreeAsync               md_gpu_free
    cudaMemcpyAsync             md_gpu_memcpy_async     (direction inferred)
    cudaMemsetAsync             md_gpu_memset_async
    kernel<<<g,b,0,s>>>(args)   md_gpu_launch
    cudaEventRecord             md_gpu_stream_record
    cudaStreamWaitEvent         md_gpu_stream_wait
    cudaStreamSynchronize       md_gpu_stream_sync
    cudaEventQuery              md_gpu_sync_is_complete
    cudaLaunchHostFunc          md_gpu_launch_host_fn
    cudaStreamBeginCapture      md_gpu_capture_begin
    cudaStreamEndCapture        md_gpu_capture_end      (also instantiates)
    cudaGraphLaunch             md_gpu_graph_launch
    cudaCreateTextureObject     md_gpu_tex_create

Device memory is a pointer
--------------------------
md_gpu_malloc returns a real address in a flat 64-bit space. Pointer
arithmetic works and means what it means:

    float* v = md_gpu_malloc(device_pool, n * 4 * sizeof(float), s);
    float* w = v + n;   // valid; the second half of the allocation

It is *not* dereferenceable on the host unless allocated with MD_GPU_MEM_HOST_*
and obtained through md_gpu_host_ptr(). This is exactly cudaMalloc's contract.

Kernel arguments
----------------
A kernel receives one pointer to a caller-defined argument struct. The backend
copies the struct into device memory and passes its address in an 8-byte push
constant, so there is no size limit and no portability cliff. Shaders declare:

    struct Args { uint n; float* dst; };
    struct Root { Args* args; };
    [[vk::push_constant]] ConstantBuffer<Root> root;

    [shader("compute")][numthreads(64,1,1)]
    void main(uint3 tid : SV_DispatchThreadID) {
        Args a = *root.args;
        ...
    }

See md_gpu.slang for the bindless texture declarations to include.

Threading
---------
  - A given md_gpu_stream_t must be used by one thread at a time. Different
    streams may be used concurrently from different threads.
  - md_gpu_malloc / md_gpu_free / texture and kernel creation are thread-safe.
  - md_gpu_device_poll should be called from the thread that wants host
    callbacks to run on it (typically the main/UI thread).
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

/* Device pointers are real 64-bit addresses. */
MD_GPU_STATIC_ASSERT(sizeof(void*) == 8, "md_gpu requires a 64-bit target");

/* =========================================================================
   Handles and value types
   ========================================================================= */

typedef struct md_gpu_device* md_gpu_device_t;
typedef struct md_gpu_stream* md_gpu_stream_t;
typedef struct md_gpu_pool*   md_gpu_pool_t;
typedef struct md_gpu_kernel* md_gpu_kernel_t;
typedef struct md_gpu_graph*  md_gpu_graph_t;

/* Device memory. Arithmetic is valid; host dereference is not (see above). */
typedef void* md_gpu_ptr_t;

/* A texture or sampler handle. Put it in an argument struct as-is; the shader
   receives it as a Slang DescriptorHandle<T>:

       struct Args {
           uint4 dim;                                  // xyz used, w spare
           DescriptorHandle<RWTexture3D<float>> vol;   // <- this field
           float* dst;
       };

   The value is a bindless heap slot: 8 bytes, but a pair of 32-bit words in
   SPIR-V against a single 8-aligned value in MSL, so like a 2-vector it has to
   sit at an offset that is a multiple of 8 -- see the ABI rule below. Zero is
   the null handle. There is no binding table, no class,
   and nothing to declare in the shader -- adding a new texture type costs no
   host code at all. Slang type-checks the handle against the resource type,
   so using a 3D storage handle where a 2D sampled texture is expected is a
   compile error rather than silent corruption. */
typedef uint64_t md_gpu_tex_t;
typedef uint64_t md_gpu_sampler_t;

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

       md_gpu_*2, md_gpu_tex_t     offset must be a multiple of 8
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
           uint64_t     dst;
       } my_args_t;

   Inside shader code -- locals, groupshared, arithmetic -- 3-vectors are fine
   and cost nothing; `a.dim.xyz` is the usual way to read one back out. The
   rule constrains the argument struct only.

   Scalars and 8-byte device pointers place themselves and need no help.
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

/* Nothing below needs it, and this header leaves no macros behind. */
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

/* ceil(count / local) groups along x. */
static inline md_gpu_grid_t md_gpu_grid_1d(uint64_t count, uint32_t local) {
    md_gpu_grid_t g;
    g.x = (uint32_t)((count + local - 1) / local);
    g.y = 1; g.z = 1;
    return g;
}

/* =========================================================================
   Errors
   ========================================================================= */

/* Human-readable description of the most recent failure on the calling
   thread, or NULL. Owned by md_gpu; valid until the next failing call on the
   same thread. Thread-safe. */
const char* md_gpu_last_error(void);

/* =========================================================================
   Device
   ========================================================================= */

typedef struct md_gpu_device_desc_t {
    /* Allocator for host-side allocations. NULL selects the default heap
       allocator. Must outlive the device. */
    struct md_allocator_i* alloc;

    /* Request backend validation (Vulkan validation layers / Metal API
       validation). Ignored if unavailable. */
    bool enable_validation;

    /* Optional debug label. */
    const char* label;
} md_gpu_device_desc_t;

typedef struct md_gpu_device_info_t {
    /* False implies UMA / integrated: MD_GPU_MEM_HOST_WRITE memory is
       directly device-accessible at no cost. An allocation hint only --
       md_gpu_memcpy_async picks the right path either way. */
    bool     is_discrete;
    uint32_t max_threads_per_group;
    uint32_t preferred_group_multiple;   /* warp / SIMD width */
    uint64_t timestamp_period_ns_num;    /* ticks -> ns, numerator   */
    uint64_t timestamp_period_ns_den;    /* ticks -> ns, denominator */
    char     name[256];
} md_gpu_device_info_t;

/* `desc` may be NULL for defaults. Returns NULL on failure; see
   md_gpu_last_error(). */
md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc);

/* Waits for all streams to go idle, then destroys the device and everything
   created from it. */
void md_gpu_device_destroy(md_gpu_device_t device);

bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* out_info);

/* Fires host callbacks whose sync point has completed, retires freed memory
   and recycles internal transient storage. Callbacks run on the calling
   thread and nowhere else. Returns the number of callbacks fired.

   Call once per frame. Cheap when there is nothing to do. */
uint32_t md_gpu_device_poll(md_gpu_device_t device);

/* =========================================================================
   Streams
   ========================================================================= */

typedef enum md_gpu_stream_kind_t {
    MD_GPU_STREAM_COMPUTE,
    MD_GPU_STREAM_TRANSFER,   /* prefers a dedicated DMA engine if present */
} md_gpu_stream_kind_t;

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t device, md_gpu_stream_kind_t kind, const char* label);

/* Waits for the stream to go idle, then destroys it. */
void md_gpu_stream_destroy(md_gpu_stream_t stream);

/* Device-owned default streams. Always valid; never destroyed by the caller. */
md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t device, md_gpu_stream_kind_t kind);

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t stream);

/* cudaEventRecord: submit whatever is pending and return the sync point that
   is signalled when it completes. Returns the none sync if nothing has been
   issued since the last record. */
md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t stream);

/* cudaStreamWaitEvent: work issued into `stream` after this call waits for
   `sync`. Work already issued is unaffected. A none sync is a no-op, and a
   sync from `stream` itself is a no-op. */
void md_gpu_stream_wait(md_gpu_stream_t stream, md_gpu_sync_t sync);

/* Submit pending work without blocking. */
void md_gpu_stream_flush(md_gpu_stream_t stream);

/* cudaStreamSynchronize: submit pending work and block until it completes. */
void md_gpu_stream_sync(md_gpu_stream_t stream);

bool md_gpu_sync_is_complete(md_gpu_sync_t sync);
void md_gpu_sync_wait(md_gpu_sync_t sync);

/* =========================================================================
   Memory
   ========================================================================= */

typedef uint32_t md_gpu_mem_flags_t;
enum {
    MD_GPU_MEM_DEVICE     = 0,       /* device-local, not host-visible      */
    MD_GPU_MEM_HOST_WRITE = 1 << 0,  /* host-writable, persistently mapped  */
    MD_GPU_MEM_HOST_READ  = 1 << 1,  /* host-readable; readback destination */
};

/* A pool is the space allocations are drawn from, and it serves exactly one
   kind of memory. That is not incidental: a VkDeviceMemory block has one
   memory type and an MTLHeap has one storage mode, so a pool maps onto them
   1:1 only if its kind is fixed at creation. It also means the grouping is
   real -- a device pool and a readback pool are different objects rather than
   two flavours multiplexed inside one.

   There is deliberately no implicit default pool. Every allocation names the
   space it comes from. This differs from cudaMallocAsync, which falls back to
   a per-device default pool; the explicitness is worth the extra line.

       md_gpu_pool_desc_t pd = {0};
       pd.flags = MD_GPU_MEM_DEVICE;
       pd.label = "md_topo scratch";
       md_gpu_pool_t pool = md_gpu_pool_create(dev, &pd);
*/
typedef struct md_gpu_pool_desc_t {
    md_gpu_mem_flags_t flags;             /* the one kind of memory served    */
    uint64_t           block_size;        /* suballocation granularity hint;
                                             0 picks a backend default        */
    uint64_t           release_threshold; /* bytes kept cached by trim; 0 is
                                             release-eagerly                  */
    const char*        label;
} md_gpu_pool_desc_t;

md_gpu_pool_t      md_gpu_pool_create(md_gpu_device_t device, const md_gpu_pool_desc_t* desc);
void               md_gpu_pool_destroy(md_gpu_pool_t pool);
md_gpu_mem_flags_t md_gpu_pool_flags(md_gpu_pool_t pool);

/* Release cached blocks down to `keep_bytes`. */
void md_gpu_pool_trim(md_gpu_pool_t pool, uint64_t keep_bytes);

/* Free everything the pool has handed out, in one call, without destroying the
   pool or returning its memory to the driver -- the CPU arena-reset pattern.
   The next round of allocations is then served entirely from cache.

   Stream-ordered exactly like md_gpu_free: the blocks become reusable at this
   point in `stream`, so work already in flight is undisturbed. Use
   md_gpu_pool_trim afterwards to hand the memory back.

   Every pointer previously obtained from this pool is dangling once this
   returns. That is the point of the call, but it does mean a pool being reset
   wholesale should not be shared with code that outlives the reset. */
void md_gpu_pool_reset(md_gpu_pool_t pool, md_gpu_stream_t stream);

typedef struct md_gpu_pool_stats_t {
    uint64_t bytes_in_use;      /* handed out right now                      */
    uint64_t bytes_reserved;    /* committed by the pool, in use or cached   */
    uint64_t bytes_cached;      /* reserved - in_use; reusable without a new
                                   device allocation                         */
    uint64_t bytes_peak_in_use; /* high-water mark, for sizing               */
    uint32_t blocks_in_use;
    uint32_t blocks_cached;
    uint64_t alloc_count;       /* md_gpu_malloc calls served               */
    uint64_t reuse_count;       /* of those, served from cache. A ratio near
                                   1 means the pool is doing its job         */
} md_gpu_pool_stats_t;

void md_gpu_pool_stats(md_gpu_pool_t pool, md_gpu_pool_stats_t* out_stats);

/* cudaMallocAsync, drawn from `pool`. The memory kind comes from the pool.
   The allocation is usable by work issued into `stream` after this point.
   Returns NULL on failure. */
md_gpu_ptr_t md_gpu_malloc(md_gpu_pool_t pool, size_t size, md_gpu_stream_t stream);

/* cudaFreeAsync. Always legal, never blocks. The memory returns to the pool
   at this point in `stream`, so later work in the same stream may reuse it
   with no synchronisation. Passing NULL is a no-op. */
void md_gpu_free(md_gpu_ptr_t ptr, md_gpu_stream_t stream);

/* Host-side view of MD_GPU_MEM_HOST_* memory, valid until the pointer is
   freed. NULL for device-local memory. */
void* md_gpu_host_ptr(md_gpu_ptr_t ptr);

/* Size of the allocation containing `ptr`, and the base of that allocation.
   Both return 0 / NULL if the pointer is not a live device allocation. */
size_t       md_gpu_ptr_size(md_gpu_ptr_t ptr);
md_gpu_ptr_t md_gpu_ptr_base(md_gpu_ptr_t ptr);

/* cudaMemcpyAsync with cudaMemcpyDefault: each pointer is looked up in the
   live-allocation map, so host-to-device, device-to-host and device-to-device
   are all this one call. Host staging is internal. */
bool md_gpu_memcpy_async(void* dst, const void* src, size_t size, md_gpu_stream_t stream);

/* Fills `size` bytes with a repeating byte value. offset and size are rounded
   to 4-byte alignment internally; sub-word tails are handled correctly. */
bool md_gpu_memset_async(md_gpu_ptr_t dst, uint8_t value, size_t size, md_gpu_stream_t stream);

/* Zero-copy upload: reserve `size` bytes and build the payload in place,
   avoiding an intermediate buffer plus memcpy. Returns a host pointer that
   lands directly in `dst` when that is safe, and in staging otherwise.

       float* p = md_gpu_upload_begin(s, coeff, n * sizeof(float));
       if (!p) return false;
       pack_coefficients(p, ...);
       md_gpu_upload_end(s);

   Returns NULL on failure, in which case upload_end must not be called.
   Only one upload may be open per stream at a time. */
void* md_gpu_upload_begin(md_gpu_stream_t stream, md_gpu_ptr_t dst, size_t size);
bool  md_gpu_upload_end(md_gpu_stream_t stream);

/* =========================================================================
   Textures
   ========================================================================= */

typedef enum md_gpu_format_t {
    MD_GPU_FORMAT_R32_FLOAT,
    MD_GPU_FORMAT_R32_UINT,
    MD_GPU_FORMAT_RGBA32_FLOAT,
    MD_GPU_FORMAT_RGBA8_UNORM,
    MD_GPU_FORMAT_COUNT,
} md_gpu_format_t;

typedef uint32_t md_gpu_tex_flags_t;
enum {
    MD_GPU_TEX_STORAGE = 1 << 0,   /* shader read/write, random access */
    MD_GPU_TEX_SAMPLED = 1 << 1,   /* shader sampled read              */
};

typedef struct md_gpu_tex_desc_t {
    uint32_t           width, height, depth;  /* depth 0 or 1 => 2D */
    md_gpu_format_t    format;
    md_gpu_tex_flags_t flags;
    const char*        label;
} md_gpu_tex_desc_t;

/* A subregion in texels. A zero extent component means "to the end along that
   axis", so a zero-initialised region covers the whole texture. */
typedef struct md_gpu_tex_region_t {
    uint32_t offset[3];
    uint32_t extent[3];
} md_gpu_tex_region_t;

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
    md_gpu_filter_t       min_filter, mag_filter;
    md_gpu_address_mode_t address_u, address_v, address_w;
    const char*           label;
} md_gpu_sampler_desc_t;

md_gpu_tex_t md_gpu_tex_create(md_gpu_device_t device, const md_gpu_tex_desc_t* desc);

/* Always legal with work in flight; retired once every stream has passed the
   last sync value that touched it. `stream` may be NULL for "any". */
void md_gpu_tex_destroy(md_gpu_tex_t tex, md_gpu_stream_t stream);

bool md_gpu_tex_desc(md_gpu_tex_t tex, md_gpu_tex_desc_t* out_desc);

/* Storage and sampled images are separate descriptor arrays, so a texture
   created with both STORAGE and SAMPLED occupies two heap slots.
   md_gpu_tex_create returns the storage handle; this returns the sampled one,
   for use in a DescriptorHandle<Texture...> field. For single-usage textures
   it returns the texture's own handle.

   Only the handle from md_gpu_tex_create identifies the texture to
   md_gpu_tex_destroy, md_gpu_tex_desc and the texture copy calls. */
md_gpu_tex_t md_gpu_tex_sampled(md_gpu_tex_t tex);

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t device, const md_gpu_sampler_desc_t* desc);
void             md_gpu_sampler_destroy(md_gpu_sampler_t sampler);

/* `region` may be NULL for the whole texture. `src` is a host pointer for
   *to_tex and a host or device pointer for *from_tex. */
bool md_gpu_memcpy_to_tex_async(md_gpu_tex_t dst, const md_gpu_tex_region_t* region,
                                const void* src, size_t size, md_gpu_stream_t stream);
bool md_gpu_memcpy_from_tex_async(void* dst, md_gpu_tex_t src, const md_gpu_tex_region_t* region,
                                  size_t size, md_gpu_stream_t stream);

/* =========================================================================
   Kernels
   ========================================================================= */

typedef struct md_gpu_kernel_desc_t {
    /* SPIR-V on Vulkan, a metallib on Metal. */
    const void* code;
    size_t      code_size;
    const char* entry_point;   /* NULL = "main" */
    const char* label;

    /* Threads per group. Required on Metal, which cannot recover it from the
       library; on Vulkan it is taken from the SPIR-V when left {0,0,0}. */
    uint32_t group_size[3];
} md_gpu_kernel_desc_t;

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t device, const md_gpu_kernel_desc_t* desc);
void            md_gpu_kernel_destroy(md_gpu_kernel_t kernel);

typedef struct md_gpu_kernel_info_t {
    uint32_t max_threads_per_group;
    uint32_t preferred_group_multiple;
    uint32_t group_size[3];
} md_gpu_kernel_info_t;

bool md_gpu_kernel_info(md_gpu_kernel_t kernel, md_gpu_kernel_info_t* out_info);

/* kernel<<<grid, block, 0, stream>>>(args).

   `args` is copied immediately; the caller may reuse or free the memory as
   soon as this returns. There is no practical size limit.

   Launches into one stream are strictly ordered and never overlap, exactly as
   in CUDA -- there is no flag to opt out. Concurrency comes from using more
   streams, which md_gpu spreads across the device's hardware queues. That is
   a deliberate omission: an "these two do not alias" flag is a promise the API
   cannot check, and getting it wrong is silent corruption. */
bool md_gpu_launch(md_gpu_stream_t stream, md_gpu_kernel_t kernel, md_gpu_grid_t grid,
                    const void* args, size_t args_size);

/* The grid is read from device memory: 3 consecutive uint32 at `grid_ptr`. */
bool md_gpu_launch_indirect(md_gpu_stream_t stream, md_gpu_kernel_t kernel, md_gpu_ptr_t grid_ptr,
                             const void* args, size_t args_size);

/* Built-in helper: writes { ceil(*count/local[0]), ceil(1/local[1]),
   ceil(1/local[2]) } to `out_grid` as 3 uint32, reading a single uint32
   element count from device memory. Turns a device-side count into an
   indirect grid without every caller writing the same three-line shader. */
bool md_gpu_make_grid(md_gpu_stream_t stream, md_gpu_ptr_t out_grid,
                      const md_gpu_ptr_t count, const uint32_t local[3]);

/* =========================================================================
   Graphs
   ========================================================================= */

/* Stream capture. Everything issued into `stream` between begin and end is
   recorded instead of executed. Because argument blocks live in device
   memory, relaunching with new parameters is a struct write -- the graph is
   never re-recorded.

   Capture is also how a worker thread records work for later launch. */
bool           md_gpu_capture_begin(md_gpu_stream_t stream, const char* label);
md_gpu_graph_t md_gpu_capture_end(md_gpu_stream_t stream);
bool           md_gpu_is_capturing(md_gpu_stream_t stream);

/* The ordinal the next md_gpu_launch into this stream will be given. Use it
   to remember which launch's arguments you want to patch later. */
uint32_t md_gpu_capture_next_index(md_gpu_stream_t stream);

void     md_gpu_graph_destroy(md_gpu_graph_t graph);
uint32_t md_gpu_graph_launch_count(md_gpu_graph_t graph);

/* Host pointer to launch `index`'s argument block. Write to it, relaunch, and
   no re-recording happens. NULL if the index is out of range. */
void* md_gpu_graph_args(md_gpu_graph_t graph, uint32_t index);

/* Replay. The graph may be launched into any stream of the same kind it was
   captured on. */
bool md_gpu_graph_launch(md_gpu_graph_t graph, md_gpu_stream_t stream);

/* =========================================================================
   Host-side ordering
   ========================================================================= */

typedef void (*md_gpu_host_fn)(void* user);

/* cudaLaunchHostFunc: `fn` runs after all work issued into `stream` before
   this call has completed. It runs inside md_gpu_device_poll(), on whichever
   thread calls that -- so the caller picks the thread and there is no locking
   in user code. */
bool md_gpu_launch_host_fn(md_gpu_stream_t stream, md_gpu_host_fn fn, void* user);

/* Same, keyed to an explicit sync point rather than the stream's current
   position. */
bool md_gpu_sync_on_complete(md_gpu_sync_t sync, md_gpu_host_fn fn, void* user);

/* =========================================================================
   Escape hatch
   ========================================================================= */

/* Program order is conservative by design. If profiling ever proves a
   specific barrier over-strict, disable automatic ordering for a region and
   place barriers by hand. Expect never to use this.

   Vulkan only. On Metal, program order is a property of the command encoder
   rather than something the backend emits, so there is no per-region barrier to
   relax: both calls are no-ops there and ordering is unaffected. Code using the
   hatch stays correct on both, it simply does not speed up on Metal. */
void md_gpu_stream_set_auto_order(md_gpu_stream_t stream, bool enabled);
void md_gpu_stream_barrier(md_gpu_stream_t stream);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

/*
md_gpu.h

Minimal, compute-focused GPU abstraction designed for:
  - Vulkan
  - Metal

Design goals:
  - Explicit resource ownership
  - Explicit synchronization (fences)
  - Async compute across frames
  - Copy-based readback / upload paths
  - C-style, flags-based, opaque handles
*/

#ifndef MD_GPU_H
#define MD_GPU_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// =============================
// Opaque handles
// =============================

typedef struct md_gpu_device*           md_gpu_device_t;
typedef struct md_gpu_queue*            md_gpu_queue_t;
typedef struct md_gpu_command_buffer*   md_gpu_command_buffer_t;
typedef struct md_gpu_pass*             md_gpu_pass_t;
typedef struct md_gpu_compute_pipeline* md_gpu_compute_pipeline_t;
typedef struct md_gpu_buffer*           md_gpu_buffer_t;
typedef struct md_gpu_image*            md_gpu_image_t;
typedef struct md_gpu_sampler*          md_gpu_sampler_t;
typedef struct md_gpu_fence*            md_gpu_fence_t;

typedef struct md_gpu_pass_id_t {
    uint64_t value;
} md_gpu_pass_id_t;

// =============================
// Flags & enums
// =============================

typedef uint32_t md_gpu_buffer_flags_t;
enum {
    MD_GPU_BUFFER_NONE        = 0,
    MD_GPU_BUFFER_CPU_VISIBLE = 1 << 0,
};

// Pipeline stages for barrier synchronization, used with md_gpu_cmd_barrier*.
// Specifying explicit src/dst stages lets the backend avoid a full pipeline stall
// when only a subset of stages are involved.
// On Metal, stage hints are accepted but ignored (resource-scoped barriers are already maximally granular).
// On Vulkan, they map to VkPipelineStageFlags2.
typedef uint32_t md_gpu_barrier_stage_t;
enum {
    MD_GPU_BARRIER_STAGE_COMPUTE  = 1 << 0,  /* compute shader dispatch */
    MD_GPU_BARRIER_STAGE_TRANSFER = 1 << 1,  /* copy / fill operations  */
};

typedef uint32_t md_gpu_image_flags_t;
enum {
    MD_GPU_IMAGE_NONE          = 0,
    MD_GPU_IMAGE_STORAGE       = 1 << 0, /* shader read/write (random access); transfer src+dst set implicitly */
    MD_GPU_IMAGE_SAMPLED       = 1 << 1, /* shader sampled read; transfer src+dst set implicitly */
    MD_GPU_IMAGE_RENDER_TARGET = 1 << 2, /* color attachment; transfer NOT set implicitly (may be tile-memory resident) */
};

typedef enum md_gpu_image_format_t {
    MD_GPU_IMAGE_FORMAT_R8_UINT,
    MD_GPU_IMAGE_FORMAT_R16_UINT,
    MD_GPU_IMAGE_FORMAT_R32_UINT,
    MD_GPU_IMAGE_FORMAT_R32_FLOAT,
    MD_GPU_IMAGE_FORMAT_RGBA8_UNORM,
    MD_GPU_IMAGE_FORMAT_RGBA16_FLOAT,
    MD_GPU_IMAGE_FORMAT_RGBA32_FLOAT,
} md_gpu_image_format_t;

enum {
    MD_GPU_MAX_BIND_SLOTS = 16,
    MD_GPU_PUSH_CONSTANTS_SLOT = (MD_GPU_MAX_BIND_SLOTS - 1),
    MD_GPU_MAX_PUSH_CONSTANTS = 256,
    MD_GPU_STORAGE_IMAGE_BINDING_BASE = MD_GPU_MAX_BIND_SLOTS,
    MD_GPU_SAMPLED_IMAGE_BINDING_BASE = 2 * MD_GPU_MAX_BIND_SLOTS,
    MD_GPU_SAMPLER_BINDING_BASE = 3 * MD_GPU_MAX_BIND_SLOTS,
};

typedef enum md_gpu_filter_t {
    MD_GPU_FILTER_NEAREST,
    MD_GPU_FILTER_LINEAR,
} md_gpu_filter_t;

typedef enum md_gpu_address_mode_t {
    MD_GPU_ADDRESS_MODE_CLAMP_TO_EDGE,
    MD_GPU_ADDRESS_MODE_REPEAT,
    MD_GPU_ADDRESS_MODE_MIRRORED_REPEAT,
} md_gpu_address_mode_t;

// =============================
// Descriptors
// =============================

typedef struct md_gpu_buffer_desc_t {
    size_t size;
    md_gpu_buffer_flags_t flags;
} md_gpu_buffer_desc_t;

typedef struct md_gpu_image_desc_t {
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    md_gpu_image_format_t format;
    md_gpu_image_flags_t flags;
} md_gpu_image_desc_t;

typedef struct md_gpu_sampler_desc_t {
    md_gpu_filter_t min_filter;
    md_gpu_filter_t mag_filter;
    md_gpu_address_mode_t address_u;
    md_gpu_address_mode_t address_v;
    md_gpu_address_mode_t address_w;
} md_gpu_sampler_desc_t;

// A subregion of a 3-D image in texel coordinates.
typedef struct md_gpu_image_region_t {
    uint32_t offset[3]; /* origin in texels (x, y, z) */
    uint32_t extent[3]; /* size   in texels (w, h, d) */
} md_gpu_image_region_t;

typedef struct md_gpu_compute_pipeline_desc_t {
    const void* shader_bytes;
    size_t      shader_byte_size;
    uint32_t threadgroup_size[3]; /* {0,0,0} = auto */
} md_gpu_compute_pipeline_desc_t;

typedef uint32_t md_gpu_pass_resource_usage_flags_t;
enum {
    MD_GPU_PASS_RESOURCE_READ         = 1 << 0,
    MD_GPU_PASS_RESOURCE_WRITE        = 1 << 1,
    MD_GPU_PASS_RESOURCE_TRANSFER_SRC = 1 << 2,
    MD_GPU_PASS_RESOURCE_TRANSFER_DST = 1 << 3,
};

typedef struct md_gpu_pass_buffer_usage_t {
    md_gpu_buffer_t buffer;
    size_t offset;
    size_t size; /* 0 = full remaining buffer */
    md_gpu_pass_resource_usage_flags_t usage;
} md_gpu_pass_buffer_usage_t;

typedef struct md_gpu_pass_image_usage_t {
    md_gpu_image_t image;
    md_gpu_image_region_t region; /* zero extent = full image */
    md_gpu_pass_resource_usage_flags_t usage;
} md_gpu_pass_image_usage_t;

typedef struct md_gpu_pass_desc_t {
    const char* label;

    const md_gpu_pass_buffer_usage_t* buffers;
    uint32_t buffer_count;

    const md_gpu_pass_image_usage_t* images;
    uint32_t image_count;
} md_gpu_pass_desc_t;


// =============================
// Device info / hints
// =============================
// Hints about the physical device, intended to help callers make allocation decisions.
// Fields may be left zero-initialized on backends that cannot determine them.
typedef struct md_gpu_device_info_t {
    // True for dedicated discrete GPUs (e.g. NVIDIA, AMD dGPU).
    // When false, assume UMA / integrated — buffers created with MD_GPU_BUFFER_CPU_VISIBLE
    // are directly accessible and staging copies can be skipped.
    bool is_discrete;

    // Human-readable device name for logging / diagnostics.
    char name[256];
} md_gpu_device_info_t;

// Populate *info with hints for the given device.
// Returns false if device is NULL, true otherwise.
// Any field that cannot be determined is left at its zero value.
bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* info);

// =============================
// Device / queue
// =============================

md_gpu_device_t md_gpu_device_create(void);
void            md_gpu_device_destroy(md_gpu_device_t device);

// Returns a compute queue for use by the calling thread.
// The implementation may return a unique queue per thread if available, otherwise a shared one.
// The handle is owned by the device — no destroy required; the caller may cache and reuse it.
// Note: submission may be internally serialized if the underlying queue is shared.
md_gpu_queue_t  md_gpu_queue_acquire(md_gpu_device_t device);

// =============================
// Buffers
// =============================

md_gpu_buffer_t md_gpu_buffer_create(md_gpu_device_t device, const md_gpu_buffer_desc_t* desc);

void md_gpu_buffer_destroy(md_gpu_buffer_t buffer);

// Returns the persistent CPU-accessible pointer for a CPU_VISIBLE buffer.
// The pointer is valid from creation until md_gpu_buffer_destroy.
// Returns NULL if the buffer was not created with MD_GPU_BUFFER_CPU_VISIBLE.
void*                 md_gpu_buffer_cpu_ptr(md_gpu_buffer_t buffer);
uint64_t              md_gpu_buffer_address(md_gpu_buffer_t buffer);
md_gpu_buffer_flags_t md_gpu_buffer_flags(md_gpu_buffer_t buffer);
size_t                md_gpu_buffer_size(md_gpu_buffer_t buffer);

// =============================
// Images (volumes, storage/sampled images)
// =============================

md_gpu_image_t md_gpu_image_create(md_gpu_device_t device, const md_gpu_image_desc_t* desc);
void md_gpu_image_destroy(md_gpu_image_t image);

bool md_gpu_image_desc_extract(md_gpu_image_t image, md_gpu_image_desc_t* out_desc);

// =============================
// Samplers (filtering/addressing for sampled images)
// =============================

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t device, const md_gpu_sampler_desc_t* desc);
void md_gpu_sampler_destroy(md_gpu_sampler_t sampler);

// =============================
// Pipelines
// =============================

md_gpu_compute_pipeline_t md_gpu_compute_pipeline_create(md_gpu_device_t device, const md_gpu_compute_pipeline_desc_t* desc);
void md_gpu_compute_pipeline_destroy(md_gpu_compute_pipeline_t pipeline);

// =============================
// Passes
// =============================
// A pass is a transient recording scope. The caller records work between begin/end,
// then receives a durable completion id from md_gpu_pass_end().
// v0 backends implement passes on top of the existing transient command buffer path;
// future backends may batch, reorder, or map pass ids to timeline/shared events.
// Queue selection is backend-owned; the v0 implementation records on the primary queue.

md_gpu_pass_t md_gpu_pass_begin(md_gpu_device_t device, const md_gpu_pass_desc_t* desc);
md_gpu_device_t md_gpu_pass_device(md_gpu_pass_t pass);

// Ends recording, submits the pass, and returns its completion id.
// A zero id means submission failed.
md_gpu_pass_id_t md_gpu_pass_end(md_gpu_pass_t pass);

bool md_gpu_pass_is_complete(md_gpu_device_t device, md_gpu_pass_id_t id);
void md_gpu_pass_wait(md_gpu_device_t device, md_gpu_pass_id_t id);

// Inserts a dependency on a previously ended pass. The initial compatibility
// implementation may resolve this as a host wait; timeline/event backends can
// lower it to a GPU-side stream/queue wait.
void md_gpu_pass_wait_pass(md_gpu_pass_t pass, md_gpu_pass_id_t id);

void md_gpu_pass_bind_compute_pipeline(md_gpu_pass_t pass, md_gpu_compute_pipeline_t pipeline);
void md_gpu_pass_bind_buffer(md_gpu_pass_t pass, uint32_t slot, md_gpu_buffer_t buffer);
void md_gpu_pass_bind_buffer_range(md_gpu_pass_t pass, uint32_t slot, md_gpu_buffer_t buffer, size_t offset, size_t size);
void md_gpu_pass_bind_image(md_gpu_pass_t pass, uint32_t slot, md_gpu_image_t image);
void md_gpu_pass_bind_sampled_image(md_gpu_pass_t pass, uint32_t slot, md_gpu_image_t image);
void md_gpu_pass_bind_sampler(md_gpu_pass_t pass, uint32_t slot, md_gpu_sampler_t sampler);
void md_gpu_pass_push_constants(md_gpu_pass_t pass, const void* data, size_t size);
void md_gpu_pass_dispatch(md_gpu_pass_t pass, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z);

void md_gpu_pass_barrier(md_gpu_pass_t pass, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);
void md_gpu_pass_barrier_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer,
                                md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);
void md_gpu_pass_barrier_image(md_gpu_pass_t pass, md_gpu_image_t image,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

void md_gpu_pass_push_debug_group(md_gpu_pass_t pass, const char* label);
void md_gpu_pass_pop_debug_group(md_gpu_pass_t pass);

void md_gpu_pass_copy_buffer(
    md_gpu_pass_t pass,
    md_gpu_buffer_t src,
    md_gpu_buffer_t dst,
    size_t size,
    size_t src_offset,
    size_t dst_offset);

void md_gpu_pass_copy_image_region_to_buffer(
    md_gpu_pass_t pass,
    md_gpu_image_t src_image,
    md_gpu_image_region_t src_region,
    md_gpu_buffer_t dst_buffer,
    size_t dst_offset);

void md_gpu_pass_copy_buffer_to_image(
    md_gpu_pass_t pass,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image);

void md_gpu_pass_copy_buffer_to_image_region(
    md_gpu_pass_t pass,
    md_gpu_buffer_t src_buffer,
    size_t src_offset,
    md_gpu_image_t dst_image,
    md_gpu_image_region_t dst_region);

void md_gpu_pass_fill_buffer(
    md_gpu_pass_t pass,
    md_gpu_buffer_t buffer,
    size_t offset,
    size_t size,
    uint8_t value);

// =============================
// Command buffers
// =============================
// Binding model:
//   Binding slots share a single 0..N-1 index namespace for buffers, storage images, sampled images, and samplers.
//   Push constants occupy a reserved slot (MD_GPU_PUSH_CONSTANTS_SLOT) and must not be used for explicit binds.
//   Vulkan convention:
//     buffers        at descriptor set 0, binding = slot;
//     storage images at binding = slot + MD_GPU_STORAGE_IMAGE_BINDING_BASE;
//     sampled images at binding = slot + MD_GPU_SAMPLED_IMAGE_BINDING_BASE;
//     samplers       at binding = slot + MD_GPU_SAMPLER_BINDING_BASE.
//   This avoids descriptor-type collisions since Vulkan cannot share a binding number between buffers and images,
//   while Metal has separate hardware index spaces for buffers, textures, and samplers.
//   Preferred pattern: pass GPU buffer device addresses via push constants using md_gpu_buffer_address()
//   and rely on the bindless address path; explicit slot binds remain available for images and legacy shaders.

// Acquire a transient command buffer from the queue's pool and begin recording.
// md_gpu_cmd_* calls encode directly into the backend command buffer.
// The command buffer is automatically recycled by md_gpu_queue_submit() and must not be used after that call.
md_gpu_command_buffer_t md_gpu_command_buffer_acquire(md_gpu_queue_t queue);

// Returns the device that owns the queue this command buffer was acquired from.
md_gpu_device_t md_gpu_command_buffer_device(md_gpu_command_buffer_t cmd);

void md_gpu_cmd_bind_compute_pipeline(md_gpu_command_buffer_t cmd, md_gpu_compute_pipeline_t pipeline);
void md_gpu_cmd_bind_buffer(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer);
void md_gpu_cmd_bind_buffer_range(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer, size_t offset, size_t size);
void md_gpu_cmd_bind_image(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_image_t image);
void md_gpu_cmd_bind_sampled_image(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_image_t image);
void md_gpu_cmd_bind_sampler(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_sampler_t sampler);
void md_gpu_cmd_push_constants(md_gpu_command_buffer_t cmd, const void* data, size_t size);
// Dispatch compute work. group_count_* map directly to Vulkan workgroup counts / Metal threadgroup counts.
void md_gpu_cmd_dispatch(md_gpu_command_buffer_t cmd, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z);

// Barrier model:
//   md_gpu_cmd_barrier(): conservative global memory barrier covering all resources in the command buffer.
//   md_gpu_cmd_barrier_buffer/image(): resource-scoped barrier. Guarantees all prior GPU writes to that
//   resource are visible to all subsequent GPU accesses of that resource within this command buffer.
//   Stage hints (src_stage / dst_stage) allow backends to narrow the stall scope; on Metal they are ignored.

void md_gpu_cmd_barrier(md_gpu_command_buffer_t cmd, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

// These will be removed once a fully bindless model is adopted.
void md_gpu_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

void md_gpu_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image,
                              md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

// Debug markers — no-op on release builds or backends without debug support.
// Groups may be nested. On Metal they map to pushDebugGroup/popDebugGroup on the active encoder;
// on Vulkan they use VK_EXT_debug_utils when available.
void md_gpu_cmd_push_debug_group(md_gpu_command_buffer_t cmd, const char* label);
void md_gpu_cmd_pop_debug_group(md_gpu_command_buffer_t cmd);

// =============================
// Copy / readback
// =============================

void md_gpu_cmd_copy_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t src,
    md_gpu_buffer_t dst,
    size_t size,
    size_t src_offset,
    size_t dst_offset);

// Copy a subregion of src_image into dst_buffer at dst_offset bytes.
// Texels are written tightly packed (row-major, no padding).
// Pass a zero-initialized region (all extents zero) to copy the full image.
void md_gpu_cmd_copy_image_region_to_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_image_t src_image,
    md_gpu_image_region_t src_region,
    md_gpu_buffer_t dst_buffer,
    size_t dst_offset);

// Copy the full src_buffer into dst_image.
void md_gpu_cmd_copy_buffer_to_image(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image);

// Copy src_buffer at src_offset bytes into a subregion of dst_image.
// Source data must be tightly packed (row-major, no padding).
// Pass a zero-initialized region (all extents zero) to copy into the full image.
void md_gpu_cmd_copy_buffer_to_image_region(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t src_buffer,
    size_t src_offset,
    md_gpu_image_t dst_image,
    md_gpu_image_region_t dst_region);

// Byte-fill: writes value repeatedly across the specified byte range.
void md_gpu_cmd_fill_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t buffer,
    size_t offset,
    size_t size,
    uint8_t value);

// =============================
// Submission & synchronization
// =============================

// Create a fence for async GPU tracking. The caller owns it and must call md_gpu_fence_destroy when done.
md_gpu_fence_t md_gpu_fence_create(md_gpu_device_t device);

// Finalize and submit a transient command buffer for GPU execution.
// If fence is NULL, blocks until GPU completion.
// If fence is non-NULL it must have been created with md_gpu_fence_create; the GPU signals it on completion.
// The command buffer is recycled by this call and must not be used afterwards.
bool md_gpu_queue_submit(md_gpu_queue_t queue, md_gpu_command_buffer_t cmd, md_gpu_fence_t fence);

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence);
void md_gpu_fence_wait(md_gpu_fence_t fence);
void md_gpu_fence_destroy(md_gpu_fence_t fence);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

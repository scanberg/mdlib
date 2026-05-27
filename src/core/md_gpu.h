/*
md_gpu.h

Minimal, compute-focused GPU abstraction designed for:
  - Vulkan
  - Metal

Design goals:
  - Explicit resource ownership
    - Explicit synchronization (event ids)
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
typedef struct md_gpu_cmd*              md_gpu_cmd_t;
typedef struct md_gpu_compute_pipeline* md_gpu_compute_pipeline_t;
typedef struct md_gpu_buffer*           md_gpu_buffer_t;
typedef struct md_gpu_image*            md_gpu_image_t;
typedef struct md_gpu_sampler*          md_gpu_sampler_t;

typedef struct md_gpu_event_t {
    uint64_t value;
} md_gpu_event_t;

// =============================
// Flags & enums
// =============================

typedef uint32_t md_gpu_buffer_flags_t;
enum {
    MD_GPU_BUFFER_NONE        = 0,
    MD_GPU_BUFFER_CPU_VISIBLE = 1 << 0,
};

// Pipeline stages for barrier synchronization, used with md_gpu_barrier*.
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

typedef uint32_t gpu_usage_flags_t;
enum {
    GPU_USAGE_READ  = 1 << 0,
    GPU_USAGE_WRITE = 1 << 1,
};

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
// Device
// =============================

md_gpu_device_t md_gpu_device_create(void);
void            md_gpu_device_destroy(md_gpu_device_t device);

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
// Command buffers
// =============================
// A command buffer is a transient recording scope. The caller records work,
// submits it, then receives a durable completion event id.

md_gpu_cmd_t md_gpu_cmd_begin(md_gpu_device_t device, const char* label);

// Ends recording, submits the command buffer, and returns its completion event.
// A zero id means submission failed.
md_gpu_event_t md_gpu_cmd_submit(md_gpu_cmd_t cmd);

bool md_gpu_event_is_complete(md_gpu_device_t device, md_gpu_event_t event);
void md_gpu_event_wait(md_gpu_device_t device, md_gpu_event_t event);

// Inserts an ordering dependency on a previously submitted event.
void md_gpu_cmd_event_wait(md_gpu_cmd_t cmd, md_gpu_event_t event);

// Declare shader-visible resource usage within a recording scope.
// These declarations are merged into cmd state and consumed by later dispatches.
// Copy/fill commands do not require md_gpu_cmd_use_*.
void md_gpu_cmd_use_buffer(md_gpu_cmd_t cmd, md_gpu_buffer_t buffer, gpu_usage_flags_t usage);
void md_gpu_cmd_use_image(md_gpu_cmd_t cmd, md_gpu_image_t image, gpu_usage_flags_t usage);
void md_gpu_cmd_use_sampled_resource(md_gpu_cmd_t cmd, md_gpu_image_t image, md_gpu_sampler_t sampler, gpu_usage_flags_t usage);

void md_gpu_cmd_bind_compute_pipeline(md_gpu_cmd_t cmd, md_gpu_compute_pipeline_t pipeline);
void md_gpu_cmd_dispatch(md_gpu_cmd_t cmd, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z, const void* root_args, size_t root_args_size);

void md_gpu_cmd_barrier(md_gpu_cmd_t cmd, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

void md_gpu_cmd_push_debug_group(md_gpu_cmd_t cmd, const char* label);
void md_gpu_cmd_pop_debug_group(md_gpu_cmd_t cmd);

void md_gpu_cmd_copy_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src,
    md_gpu_buffer_t dst,
    size_t size,
    size_t src_offset,
    size_t dst_offset);

void md_gpu_cmd_copy_image_region_to_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_image_t src_image,
    md_gpu_image_region_t src_region,
    md_gpu_buffer_t dst_buffer,
    size_t dst_offset);

void md_gpu_cmd_copy_buffer_to_image(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image);

void md_gpu_cmd_copy_buffer_to_image_region(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src_buffer,
    size_t src_offset,
    md_gpu_image_t dst_image,
    md_gpu_image_region_t dst_region);

void md_gpu_cmd_fill_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t buffer,
    size_t offset,
    size_t size,
    uint8_t value);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

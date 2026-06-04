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
typedef struct md_gpu_queue*            md_gpu_queue_t;
typedef struct md_gpu_cmd*              md_gpu_cmd_t;
typedef struct md_gpu_compute_pipeline* md_gpu_compute_pipeline_t;
typedef struct md_gpu_buffer*           md_gpu_buffer_t;
typedef struct md_gpu_image*            md_gpu_image_t;
typedef struct md_gpu_sampler*          md_gpu_sampler_t;
typedef struct md_gpu_readback*         md_gpu_readback_t;

typedef struct md_gpu_event_t {
    md_gpu_queue_t queue;
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

typedef uint32_t md_gpu_usage_flags_t;
enum {
    MD_GPU_USAGE_READ  = 1 << 0,
    MD_GPU_USAGE_WRITE = 1 << 1,
};

typedef enum md_gpu_resource_kind_t {
    /* Buffer usage declaration for hazard tracking; binding is via root_args/device addresses. */
    MD_GPU_RESOURCE_BUFFER_USAGE,
    MD_GPU_RESOURCE_STORAGE_IMAGE,
    MD_GPU_RESOURCE_SAMPLED_IMAGE,
    MD_GPU_RESOURCE_SAMPLER,
} md_gpu_resource_kind_t;

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

typedef struct md_gpu_buffer_image_copy_t {
    md_gpu_image_region_t image_region;
    size_t buffer_offset;
    size_t bytes_per_row;   /* 0 = tightly packed */
    size_t bytes_per_image; /* 0 = tightly packed */
} md_gpu_buffer_image_copy_t;

typedef struct md_gpu_resource_t {
    md_gpu_resource_kind_t kind;
    md_gpu_usage_flags_t usage;
    uint32_t set;
    uint32_t binding;
    md_gpu_buffer_t buffer;
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
} md_gpu_resource_t;

typedef struct md_gpu_buffer_resource_t {
    md_gpu_buffer_t buffer;
    uint64_t offset;
    md_gpu_usage_flags_t usage;
} md_gpu_buffer_resource_t;

typedef struct md_gpu_resource_binding_t {
    md_gpu_resource_kind_t kind;
    uint32_t set;
    uint32_t binding;
    uint32_t backend_binding;
} md_gpu_resource_binding_t;

typedef struct md_gpu_compute_pipeline_desc_t {
    const void* shader_bytes;
    size_t      shader_byte_size;
    const char* entry_point; /* NULL = "main" */
    const char* label;       /* optional debug label */
    uint32_t threadgroup_size[3]; /* {0,0,0} = auto */
    const md_gpu_resource_binding_t* resource_bindings;
    uint32_t resource_binding_count;
} md_gpu_compute_pipeline_desc_t;

// A single compute dispatch description.
// Buffers participate in usage registration even when their device addresses are
// passed via root_args. Non-buffer resources use logical descriptor set + binding.
typedef struct md_gpu_compute_dispatch_t {
    md_gpu_compute_pipeline_t pipeline;
    const md_gpu_resource_t* resources;
    uint32_t resource_count;
    uint32_t group_count[3];
    const void* root_args;
    size_t root_args_size;
} md_gpu_compute_dispatch_t;

typedef struct md_gpu_queue_wait_t {
    md_gpu_event_t event;
    md_gpu_barrier_stage_t dst_stage;
} md_gpu_queue_wait_t;

typedef struct md_gpu_queue_submit_desc_t {
    const md_gpu_cmd_t* cmds;
    size_t cmd_count;
    const md_gpu_queue_wait_t* waits;
    size_t wait_count;
    const char* label;
} md_gpu_queue_submit_desc_t;

typedef struct md_gpu_transient_t {
    md_gpu_buffer_t buffer;
    void*           cpu_ptr;
    size_t          offset;
    size_t          size;
} md_gpu_transient_t;

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

static inline bool md_gpu_event_is_valid(md_gpu_event_t event) {
    return event.value != 0 && event.queue != NULL;
}

// Populate *info with hints for the given device.
// Returns false if device is NULL, true otherwise.
// Any field that cannot be determined is left at its zero value.
bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* info);

// =============================
// Device
// =============================

md_gpu_device_t md_gpu_device_create(void);
void            md_gpu_device_destroy(md_gpu_device_t device);

// Device-owned logical queues. These handles always exist; backends may map
// several logical queues to the same physical queue when dedicated hardware queues are unavailable.
// The graphics queue is the primary queue for rendering and general-purpose compute.
// The transfer queue is an async compute queue optimized for copy/fill operations, if supported by the device.
md_gpu_queue_t md_gpu_queue_graphics(md_gpu_device_t device);
md_gpu_queue_t md_gpu_queue_compute(md_gpu_device_t device);
md_gpu_queue_t md_gpu_queue_transfer(md_gpu_device_t device);

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
// Command buffers and queue submission
// =============================
// A command buffer is a transient recording scope created from a queue.
// A handle returned from md_gpu_cmd_begin must be closed by either
// md_gpu_queue_submit, md_gpu_queue_submit_one, or md_gpu_cmd_discard.
// For now, md_gpu_queue_submit must be called from the main thread. Command
// recording may happen on worker threads as long as each command is recorded by
// only one thread and ownership is handed to the main thread before submission.

md_gpu_cmd_t md_gpu_cmd_begin(md_gpu_queue_t queue, const char* label);
bool md_gpu_cmd_end(md_gpu_cmd_t cmd);
void md_gpu_cmd_discard(md_gpu_cmd_t cmd);

// Allocates transient memory for data uploads. The returned buffer and CPU pointer are valid until the command buffer is ended or discarded.
md_gpu_transient_t md_gpu_cmd_temp_alloc(md_gpu_cmd_t cmd, size_t size);

md_gpu_event_t md_gpu_queue_submit(md_gpu_queue_t queue, const md_gpu_queue_submit_desc_t* desc);

static inline md_gpu_event_t md_gpu_queue_submit_one(md_gpu_queue_t queue, md_gpu_cmd_t cmd) {
    const md_gpu_queue_submit_desc_t desc = { .cmds = &cmd, .cmd_count = 1 };
    return md_gpu_queue_submit(queue, &desc);
}

static inline md_gpu_event_t md_gpu_queue_submit_one_after(md_gpu_queue_t queue, md_gpu_cmd_t cmd, md_gpu_event_t wait, md_gpu_barrier_stage_t dst_stage) {
    const md_gpu_queue_wait_t waits[] = {{ .event = wait, .dst_stage = dst_stage }};
    const md_gpu_queue_submit_desc_t desc = { .cmds = &cmd, .cmd_count = 1, .waits = waits, .wait_count = 1 };
    return md_gpu_queue_submit(queue, &desc);
}

// CPU poll and wait for event
bool md_gpu_event_is_complete(md_gpu_event_t event);
void md_gpu_event_wait(md_gpu_event_t event);

// Records one compute dispatch.
// Resources are consumed by this dispatch only and cleared after it is encoded.
// Copy/fill commands do not use md_gpu_compute_dispatch_t.
bool md_gpu_cmd_dispatch(md_gpu_cmd_t cmd, const md_gpu_compute_dispatch_t* dispatch);

bool md_gpu_cmd_barrier(md_gpu_cmd_t cmd, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage);

void md_gpu_cmd_push_debug_group(md_gpu_cmd_t cmd, const char* label);
void md_gpu_cmd_pop_debug_group(md_gpu_cmd_t cmd);

bool md_gpu_cmd_copy_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src,
    md_gpu_buffer_t dst,
    size_t size,
    size_t src_offset,
    size_t dst_offset);

bool md_gpu_cmd_copy_image_to_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_image_t src_image,
    md_gpu_buffer_t dst_buffer,
    const md_gpu_buffer_image_copy_t* copy);

bool md_gpu_cmd_copy_buffer_to_image_layout(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image,
    const md_gpu_buffer_image_copy_t* copy);

static inline bool md_gpu_cmd_copy_image_region_to_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_image_t src_image,
    md_gpu_image_region_t src_region,
    md_gpu_buffer_t dst_buffer,
    size_t dst_offset)
{
    md_gpu_buffer_image_copy_t copy = { .image_region = src_region, .buffer_offset = dst_offset };
    return md_gpu_cmd_copy_image_to_buffer(cmd, src_image, dst_buffer, &copy);
}

static inline bool md_gpu_cmd_copy_buffer_to_image(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image)
{
    return md_gpu_cmd_copy_buffer_to_image_layout(cmd, src_buffer, dst_image, NULL);
}

static inline bool md_gpu_cmd_copy_buffer_to_image_region(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t src_buffer,
    size_t src_offset,
    md_gpu_image_t dst_image,
    md_gpu_image_region_t dst_region)
{
    md_gpu_buffer_image_copy_t copy = { .image_region = dst_region, .buffer_offset = src_offset };
    return md_gpu_cmd_copy_buffer_to_image_layout(cmd, src_buffer, dst_image, &copy);
}

bool md_gpu_cmd_fill_buffer(
    md_gpu_cmd_t cmd,
    md_gpu_buffer_t buffer,
    size_t offset,
    size_t size,
    uint8_t value);

// =============================
// Readback (GPU to CPU copy)
// =============================
// Readback operations are issued as async transfer queue submissions that return an md_gpu_readback_t handle.
// The buffer and CPU pointer in the returned handle are valid until the readback is complete and the handle is destroyed.
// The caller must poll or wait for completion before reading from the CPU pointer.
md_gpu_readback_t md_gpu_readback_buffer(md_gpu_buffer_t src_buffer, size_t src_offset, size_t size, md_gpu_event_t after);
md_gpu_readback_t md_gpu_readback_image(md_gpu_image_t src_image, md_gpu_image_region_t src_region, md_gpu_event_t after);

bool md_gpu_readback_is_complete(md_gpu_readback_t readback);
void md_gpu_readback_wait(md_gpu_readback_t readback);

void* md_gpu_readback_cpu_ptr(md_gpu_readback_t readback);
void md_gpu_readback_destroy(md_gpu_readback_t readback);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

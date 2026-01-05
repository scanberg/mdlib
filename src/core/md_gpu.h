/*
md_gpu.h

```
Minimal, compute-focused GPU abstraction designed for:
  - Vulkan
  - Metal

Design goals:
  - Explicit resource ownership
  - Explicit synchronization (fences)
  - Async compute across frames
  - Copy-based readback / upload paths
  - C-style, flags-based, opaque handles
```

*/

#ifndef MD_GPU_H
#define MD_GPU_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =============================
Opaque handles
============================= */

typedef struct md_gpu_device*           md_gpu_device_t;
typedef struct md_gpu_queue*            md_gpu_queue_t;
typedef struct md_gpu_command_buffer*   md_gpu_command_buffer_t;
typedef struct md_gpu_compute_pipeline* md_gpu_compute_pipeline_t;
typedef struct md_gpu_buffer*           md_gpu_buffer_t;
typedef struct md_gpu_image*            md_gpu_image_t;
typedef struct md_gpu_fence*            md_gpu_fence_t;

/* =============================
Flags & enums
============================= */

typedef uint32_t md_gpu_buffer_flags_t;
enum {
    MD_GPU_BUFFER_NONE        = 0,
    MD_GPU_BUFFER_CPU_VISIBLE = 1 << 0,
};

typedef uint32_t md_gpu_image_flags_t;
enum {
    MD_GPU_IMAGE_NONE         = 0,
    MD_GPU_IMAGE_STORAGE     = 1 << 0, /* shader read/write */
    MD_GPU_IMAGE_CPU_VISIBLE = 1 << 1, /* readback via copy */
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

/* =============================
Descriptors
============================= */

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

typedef struct md_gpu_shader_blob_t {
    const void* data;
    size_t size;
} md_gpu_shader_blob_t;

typedef struct md_gpu_compute_pipeline_desc_t {
    md_gpu_shader_blob_t shader;
    uint32_t threadgroup_size[3]; /* {0,0,0} = auto */
} md_gpu_compute_pipeline_desc_t;


/* =============================
Device / queue
============================= */

md_gpu_device_t md_gpu_create_device(void);
void            md_gpu_destroy_device(md_gpu_device_t device);

/* Returns a compute queue intended for use by the calling thread.
    The implementation may return a unique queue if available, otherwise a shared queue.
    The returned handle is owned by the device (no destroy required).
    The caller may cache and reuse the returned queue.
    Note: submission may be internally serialized if the queue is shared. */
md_gpu_queue_t  md_gpu_acquire_compute_queue(md_gpu_device_t device);

/* =============================
Buffers
============================= */

md_gpu_buffer_t md_gpu_create_buffer(md_gpu_device_t device, const md_gpu_buffer_desc_t* desc);

void md_gpu_destroy_buffer(md_gpu_buffer_t buffer);

void* md_gpu_map_buffer(md_gpu_buffer_t buffer);
void  md_gpu_unmap_buffer(md_gpu_buffer_t buffer);

/* =============================
Images (volumes, storage images)
============================= */

md_gpu_image_t md_gpu_create_image(md_gpu_device_t device, const md_gpu_image_desc_t* desc);

void md_gpu_destroy_image(md_gpu_image_t image);

/* =============================
Pipelines
============================= */

md_gpu_compute_pipeline_t md_gpu_create_compute_pipeline(md_gpu_device_t device, const md_gpu_compute_pipeline_desc_t* desc);
void md_gpu_destroy_compute_pipeline(md_gpu_compute_pipeline_t pipeline);

/* =============================
Command buffers
============================= */
/* Binding model notes:
     - Binding slots are a simple 0..N-1 namespace.
     - Push constants are implemented via a reserved buffer binding slot on backends
         that share the same index space for buffers and "setBytes" (e.g. Metal).
     - Vulkan binding convention (recommended):
         - Buffers are declared/bound at descriptor set 0, binding = slot.
         - Images are declared/bound at descriptor set 0, binding = slot + MD_GPU_MAX_BIND_SLOTS.
       This avoids descriptor-type collisions since Vulkan cannot have a buffer and an image at the
       same descriptor binding number, while Metal has separate index spaces for buffers/textures.
     - Do not bind user buffers/images to the reserved slot. */
enum {
    MD_GPU_MAX_BIND_SLOTS = 16,
    MD_GPU_PUSH_CONSTANTS_SLOT = (MD_GPU_MAX_BIND_SLOTS - 1),
    MD_GPU_MAX_PUSH_CONSTANTS = 256,
};

/* Command buffers are pooled and owned by the device.
    Acquire returns a command buffer ready for recording.
    The command buffer is automatically recycled by md_gpu_queue_submit(). */
md_gpu_command_buffer_t md_gpu_acquire_command_buffer(md_gpu_device_t device);

void md_gpu_cmd_bind_compute_pipeline(md_gpu_command_buffer_t cmd, md_gpu_compute_pipeline_t pipeline);
void md_gpu_cmd_bind_buffer(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer);
void md_gpu_cmd_bind_buffer_range(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer, size_t offset, size_t size);
void md_gpu_cmd_bind_image(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_image_t image);
void md_gpu_cmd_push_constants(md_gpu_command_buffer_t cmd, const void* data, size_t size);
void md_gpu_cmd_dispatch(md_gpu_command_buffer_t cmd, uint32_t total_x, uint32_t total_y, uint32_t total_z);

/* Barrier model:
     - md_gpu_cmd_barrier(): conservative global barrier for the command buffer.
     - md_gpu_cmd_barrier_buffer/image(): resource-scoped barrier.
         Guarantees all prior GPU writes to that resource in this command buffer are
         visible to all subsequent GPU reads/writes of that resource in this command buffer. */
void md_gpu_cmd_barrier(md_gpu_command_buffer_t cmd);
void md_gpu_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer);
void md_gpu_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image);

/* =============================
Copy / readback
============================= */

void md_gpu_cmd_copy_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t src,
    md_gpu_buffer_t dst,
    size_t size,
    size_t src_offset,
    size_t dst_offset);

void md_gpu_cmd_copy_image_to_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_image_t src_image,
    md_gpu_buffer_t dst_buffer);

void md_gpu_cmd_copy_buffer_to_image(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t src_buffer,
    md_gpu_image_t dst_image);

void md_gpu_cmd_fill_buffer(
    md_gpu_command_buffer_t cmd,
    md_gpu_buffer_t buffer,
    size_t offset,
    size_t size,
    uint32_t value);

/* =============================
Submission & synchronization
============================= */

md_gpu_fence_t md_gpu_queue_submit(md_gpu_queue_t queue, md_gpu_command_buffer_t cmd);

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence);
void md_gpu_fence_wait(md_gpu_fence_t fence);
void md_gpu_destroy_fence(md_gpu_fence_t fence);

#ifdef __cplusplus
}
#endif

#endif /* MD_GPU_H */

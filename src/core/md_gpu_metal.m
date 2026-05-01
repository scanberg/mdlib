// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <core/md_common.h> // Assert
#include <stdlib.h>         // malloc, free

/* Limits are defined in md_gpu.h */

// Metal-specific: maximum number of recorded commands per command buffer
#define MD_GPU_METAL_MAX_COMMANDS 1024
#define MD_GPU_MAX_COMMANDS MD_GPU_METAL_MAX_COMMANDS

/* =============================
   Internal structs
   ============================= */

struct md_gpu_queue {
    struct md_gpu_device* dev;
    struct md_gpu_command_buffer* cmd_free_list;
};

struct md_gpu_device {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;

    struct md_gpu_queue compute_queue;
};

struct md_gpu_buffer {
    id<MTLBuffer> buffer;
    void* cpu_ptr;  /* non-NULL for CPU_VISIBLE buffers; valid until destroy */
    size_t size;
    md_gpu_buffer_flags_t flags;
};

struct md_gpu_image {
    id<MTLTexture> texture;
    md_gpu_image_desc_t desc;
};

struct md_gpu_compute_pipeline {
    id<MTLComputePipelineState> pso;
    uint32_t tg_size[3];
};

typedef enum md_gpu_cmd_type_t {
    CMD_DISPATCH,
    CMD_COPY_BUFFER,
    CMD_COPY_IMAGE_TO_BUFFER,
    CMD_COPY_IMAGE_REGION_TO_BUFFER,
    CMD_COPY_BUFFER_TO_IMAGE,
    CMD_COPY_BUFFER_TO_IMAGE_REGION,
    CMD_FILL_BUFFER,
    CMD_BARRIER,
    CMD_BARRIER_BUFFER,
    CMD_BARRIER_IMAGE
} md_gpu_cmd_type_t;

struct md_gpu_cmd_copy_buffer {
    struct md_gpu_buffer* src;
    struct md_gpu_buffer* dst;
    size_t size;
    size_t src_offset;
    size_t dst_offset;
};

struct md_gpu_cmd_copy_image_to_buffer {
    struct md_gpu_image* image;
    struct md_gpu_buffer* buffer;
};

struct md_gpu_cmd_copy_image_region_to_buffer {
    struct md_gpu_image*  image;
    struct md_gpu_buffer* buffer;
    md_gpu_image_region_t region;
    size_t dst_offset;
};

struct md_gpu_cmd_copy_buffer_to_image {
    struct md_gpu_buffer* buffer;
    struct md_gpu_image* image;
};

struct md_gpu_cmd_copy_buffer_to_image_region {
    struct md_gpu_buffer* buffer;
    struct md_gpu_image*  image;
    md_gpu_image_region_t region;
    size_t src_offset;
};

struct md_gpu_cmd_fill_buffer {
    struct md_gpu_buffer* buffer;
    size_t offset;
    size_t size;
    uint8_t value;
};

struct md_gpu_cmd_barrier_buffer {
    struct md_gpu_buffer* buffer;
};

struct md_gpu_cmd_barrier_image {
    struct md_gpu_image* image;
};

struct md_gpu_cmd_dispatch {
    uint32_t group_count[3];
    uint8_t  push_constants[MD_GPU_MAX_PUSH_CONSTANTS];
    size_t   push_constant_size;
    /* Full binding snapshot captured at record time — not read from live cmd->bound_* at submit */
    struct md_gpu_compute_pipeline* pipeline;
    struct md_gpu_buffer*           buffers[MD_GPU_MAX_BIND_SLOTS];
    size_t                          buffer_offsets[MD_GPU_MAX_BIND_SLOTS];
    struct md_gpu_image*            images[MD_GPU_MAX_BIND_SLOTS];
};

struct md_gpu_recorded_cmd {
    md_gpu_cmd_type_t type;
    union {
        struct md_gpu_cmd_dispatch dispatch;
        struct md_gpu_cmd_copy_buffer copy_buf;
        struct md_gpu_cmd_copy_image_to_buffer        copy_img_buf;
        struct md_gpu_cmd_copy_image_region_to_buffer  copy_img_region_buf;
        struct md_gpu_cmd_copy_buffer_to_image         copy_buf_img;
        struct md_gpu_cmd_copy_buffer_to_image_region  copy_buf_img_region;
        struct md_gpu_cmd_fill_buffer fill_buf;
        struct md_gpu_cmd_barrier_buffer barrier_buf;
        struct md_gpu_cmd_barrier_image  barrier_img;
    } u;
};

struct md_gpu_command_buffer {
    struct md_gpu_device* dev;

    struct md_gpu_compute_pipeline* bound_pipeline;
    struct md_gpu_buffer* bound_buffers[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_offsets[MD_GPU_MAX_BIND_SLOTS];
    struct md_gpu_image*  bound_images[MD_GPU_MAX_BIND_SLOTS];
    /* push_constants/push_constant_size: staging area, snapshotted into each dispatch record */
    uint8_t push_constants[MD_GPU_MAX_PUSH_CONSTANTS];
    size_t  push_constant_size;

    struct md_gpu_recorded_cmd commands[MD_GPU_METAL_MAX_COMMANDS];
    uint32_t command_count;

    struct md_gpu_command_buffer* next_free;
};

struct md_gpu_fence {
    id<MTLCommandBuffer> cmd;
};

/* =============================
   Helpers
   ============================= */

static MTLPixelFormat to_mtl_format(md_gpu_image_format_t fmt) {
    switch (fmt) {
        case MD_GPU_IMAGE_FORMAT_R8_UINT:      return MTLPixelFormatR8Uint;
        case MD_GPU_IMAGE_FORMAT_R16_UINT:     return MTLPixelFormatR16Uint;
        case MD_GPU_IMAGE_FORMAT_R32_UINT:     return MTLPixelFormatR32Uint;
        case MD_GPU_IMAGE_FORMAT_R32_FLOAT:    return MTLPixelFormatR32Float;
        case MD_GPU_IMAGE_FORMAT_RGBA8_UNORM:  return MTLPixelFormatRGBA8Unorm;
        case MD_GPU_IMAGE_FORMAT_RGBA16_FLOAT: return MTLPixelFormatRGBA16Float;
        case MD_GPU_IMAGE_FORMAT_RGBA32_FLOAT: return MTLPixelFormatRGBA32Float;
        default: return MTLPixelFormatInvalid;
    }
}

static inline uint32_t md_gpu_format_bytes_per_pixel(md_gpu_image_format_t fmt) {
    switch (fmt) {
        case MD_GPU_IMAGE_FORMAT_R8_UINT:       return 1;
        case MD_GPU_IMAGE_FORMAT_R16_UINT:      return 2;
        case MD_GPU_IMAGE_FORMAT_R32_UINT:      return 4;
        case MD_GPU_IMAGE_FORMAT_R32_FLOAT:     return 4;
        case MD_GPU_IMAGE_FORMAT_RGBA8_UNORM:   return 4;
        case MD_GPU_IMAGE_FORMAT_RGBA16_FLOAT:  return 8;
        case MD_GPU_IMAGE_FORMAT_RGBA32_FLOAT:  return 16;
        default: return 0;
    }
}

static inline void cmd_reset(struct md_gpu_command_buffer* cmd) {
    cmd->bound_pipeline = NULL;
    MEMSET(cmd->bound_buffers, 0, sizeof(cmd->bound_buffers));
    MEMSET(cmd->bound_buffer_offsets, 0, sizeof(cmd->bound_buffer_offsets));
    MEMSET(cmd->bound_images,  0, sizeof(cmd->bound_images));
    cmd->push_constant_size = 0;
    cmd->command_count = 0;
    cmd->next_free = NULL;
}

/* =============================
   Device / queue
   ============================= */

md_gpu_device_t md_gpu_device_create(void) {
    struct md_gpu_device* dev = (struct md_gpu_device*)calloc(1, sizeof(struct md_gpu_device));
    dev->device = MTLCreateSystemDefaultDevice();
    ASSERT(dev->device);
    dev->queue = [dev->device newCommandQueue];
    dev->compute_queue.dev = dev;
    /* cmd_free_list is zero-initialized by calloc */
    return dev;
}

void md_gpu_device_destroy(md_gpu_device_t device) {
    struct md_gpu_command_buffer* it = device->compute_queue.cmd_free_list;
    while (it) {
        struct md_gpu_command_buffer* next = it->next_free;
        free(it);
        it = next;
    }
    device->compute_queue.cmd_free_list = NULL;

    [device->queue release];
    [device->device release];
    free(device);
}

md_gpu_queue_t md_gpu_queue_acquire(md_gpu_device_t device) {
    return &device->compute_queue;
}

/* =============================
   Buffers
   ============================= */

md_gpu_buffer_t md_gpu_buffer_create(md_gpu_device_t device,
                                     const md_gpu_buffer_desc_t* desc) {
    struct md_gpu_buffer* buf = (struct md_gpu_buffer*)calloc(1, sizeof(struct md_gpu_buffer));
    buf->size = desc->size;
    buf->flags = desc->flags;

    MTLResourceOptions opts = MTLResourceStorageModePrivate;
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE)
        opts = MTLResourceStorageModeShared;

    buf->buffer = [device->device newBufferWithLength:desc->size options:opts];
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE)
        buf->cpu_ptr = [buf->buffer contents];
    return buf;
}

void md_gpu_buffer_destroy(md_gpu_buffer_t buffer) {
    [buffer->buffer release];
    free(buffer);
}

void* md_gpu_buffer_cpu_ptr(md_gpu_buffer_t buffer) {
    ASSERT(buffer->flags & MD_GPU_BUFFER_CPU_VISIBLE);
    return buffer->cpu_ptr;
}

uint64_t md_gpu_buffer_address(md_gpu_buffer_t buffer) {
    return (uint64_t)[buffer->buffer gpuAddress];
}

md_gpu_buffer_flags_t md_gpu_buffer_flags(md_gpu_buffer_t buffer) {
    if (!buffer) return MD_GPU_BUFFER_NONE;
    return buffer->flags;
}

size_t md_gpu_buffer_size(md_gpu_buffer_t buffer) {
    if (!buffer) return 0;
    return buffer->size;
}

/* =============================
   Images
   ============================= */

md_gpu_image_t md_gpu_image_create(md_gpu_device_t device,
                                   const md_gpu_image_desc_t* desc) {
    struct md_gpu_image* img = (struct md_gpu_image*)calloc(1, sizeof(struct md_gpu_image));
    img->desc = *desc;

    MTLTextureDescriptor* td = [[MTLTextureDescriptor alloc] init];
    td.textureType = desc->depth > 1 ? MTLTextureType3D : MTLTextureType2D;
    td.width  = desc->width;
    td.height = desc->height;
    td.depth  = desc->depth;
    td.pixelFormat = to_mtl_format(desc->format);
    td.usage = MTLTextureUsageUnknown;
    if (desc->flags & MD_GPU_IMAGE_STORAGE) {
        td.usage |= MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;
    }
    if (desc->flags & MD_GPU_IMAGE_SAMPLED) {
        td.usage |= MTLTextureUsageShaderRead;
    }
    if (desc->flags & MD_GPU_IMAGE_RENDER_TARGET) {
        td.usage |= MTLTextureUsageRenderTarget;
    }

    img->texture = [device->device newTextureWithDescriptor:td];
    [td release];
    return img;
}

void md_gpu_image_destroy(md_gpu_image_t image) {
    [image->texture release];
    free(image);
}

/* =============================
   Pipelines
   ============================= */

md_gpu_compute_pipeline_t md_gpu_compute_pipeline_create(
    md_gpu_device_t device,
    const md_gpu_compute_pipeline_desc_t* desc)
{
    NSError* err = nil;

    dispatch_data_t data = dispatch_data_create(
        desc->shader_bytes,
        desc->shader_byte_size,
        dispatch_get_main_queue(),
        DISPATCH_DATA_DESTRUCTOR_DEFAULT);

    id<MTLLibrary> lib =
        [device->device newLibraryWithData:data error:&err];
    dispatch_release(data);
    ASSERT(lib && !err);

    // Slang renames 'main' to 'main_0' for Metal (Metal reserves 'main').
    // spirv-cross used 'main0'; try both for compatibility.
    id<MTLFunction> fn = [lib newFunctionWithName:@"main_0"];
    if (!fn) fn = [lib newFunctionWithName:@"main0"];
    ASSERT(fn);

    struct md_gpu_compute_pipeline* p =
        (struct md_gpu_compute_pipeline*)calloc(1, sizeof(struct md_gpu_compute_pipeline));

    p->pso = [device->device newComputePipelineStateWithFunction:fn error:&err];
    ASSERT(p->pso && !err);

    NSUInteger max_threads = [p->pso maxTotalThreadsPerThreadgroup];
    NSUInteger exec_width  = [p->pso threadExecutionWidth];

    if (desc->threadgroup_size[0] == 0 &&
        desc->threadgroup_size[1] == 0 &&
        desc->threadgroup_size[2] == 0) {
        // Auto-configure threadgroup size
        p->tg_size[0] = (uint32_t)exec_width;
        p->tg_size[1] = (uint32_t)(max_threads / exec_width);
        if (p->tg_size[1] == 0) p->tg_size[1] = 1;
        p->tg_size[2] = 1;
    } else {
        MEMCPY(p->tg_size, desc->threadgroup_size, sizeof(uint32_t) * 3);
    }

    [fn release];
    [lib release];

    return p;
}

void md_gpu_compute_pipeline_destroy(md_gpu_compute_pipeline_t pipeline) {
    [pipeline->pso release];
    free(pipeline);
}

/* =============================
   Command buffers
   ============================= */

md_gpu_command_buffer_t md_gpu_command_buffer_acquire(md_gpu_queue_t queue) {
    struct md_gpu_command_buffer* cmd = queue->cmd_free_list;
    if (cmd) {
        queue->cmd_free_list = cmd->next_free;
    } else {
        cmd = (struct md_gpu_command_buffer*)calloc(1, sizeof(struct md_gpu_command_buffer));
        cmd->dev = queue->dev;
    }

    cmd_reset(cmd);
    return cmd;
}

void md_gpu_cmd_bind_compute_pipeline(md_gpu_command_buffer_t cmd,
                                     md_gpu_compute_pipeline_t pipeline) {
    cmd->bound_pipeline = pipeline;
}

void md_gpu_cmd_bind_buffer(md_gpu_command_buffer_t cmd,
                            uint32_t slot,
                            md_gpu_buffer_t buffer) {
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    ASSERT(slot != MD_GPU_PUSH_CONSTANTS_SLOT);
    cmd->bound_buffers[slot] = buffer;
    cmd->bound_buffer_offsets[slot] = 0;
}

void md_gpu_cmd_bind_buffer_range(md_gpu_command_buffer_t cmd,
                                  uint32_t slot,
                                  md_gpu_buffer_t buffer,
                                  size_t offset,
                                  size_t size) {
    (void)size;
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    ASSERT(slot != MD_GPU_PUSH_CONSTANTS_SLOT);
    ASSERT(buffer);
    ASSERT(offset <= buffer->size);
    ASSERT(offset + size <= buffer->size);
    cmd->bound_buffers[slot] = buffer;
    cmd->bound_buffer_offsets[slot] = offset;
}

void md_gpu_cmd_bind_image(md_gpu_command_buffer_t cmd,
                           uint32_t slot,
                           md_gpu_image_t image) {
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    ASSERT(slot != MD_GPU_PUSH_CONSTANTS_SLOT);
    cmd->bound_images[slot] = image;
}

void md_gpu_cmd_push_constants(md_gpu_command_buffer_t cmd,
                               const void* data,
                               size_t size) {
    ASSERT(size <= MD_GPU_MAX_PUSH_CONSTANTS);
    MEMCPY(cmd->push_constants, data, size);
    cmd->push_constant_size = size;
}

void md_gpu_cmd_dispatch(md_gpu_command_buffer_t cmd,
                         uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_DISPATCH;
    rc->u.dispatch.group_count[0] = group_count_x;
    rc->u.dispatch.group_count[1] = group_count_y;
    rc->u.dispatch.group_count[2] = group_count_z;
    /* Snapshot push constants at record time so each dispatch is independent */
    rc->u.dispatch.push_constant_size = cmd->push_constant_size;
    if (cmd->push_constant_size > 0)
        MEMCPY(rc->u.dispatch.push_constants, cmd->push_constants, cmd->push_constant_size);
    /* Snapshot pipeline and bindings — the live cmd->bound_* state must not be read at submit time */
    rc->u.dispatch.pipeline = cmd->bound_pipeline;
    MEMCPY(rc->u.dispatch.buffers,        cmd->bound_buffers,        sizeof(cmd->bound_buffers));
    MEMCPY(rc->u.dispatch.buffer_offsets, cmd->bound_buffer_offsets, sizeof(cmd->bound_buffer_offsets));
    MEMCPY(rc->u.dispatch.images,         cmd->bound_images,         sizeof(cmd->bound_images));
}

void md_gpu_cmd_barrier(md_gpu_command_buffer_t cmd) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_BARRIER;
}

void md_gpu_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    ASSERT(buffer);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_BARRIER_BUFFER;
    rc->u.barrier_buf.buffer = buffer;
}

void md_gpu_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    ASSERT(image);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_BARRIER_IMAGE;
    rc->u.barrier_img.image = image;
}

/* _ex variants: stage hints accepted but ignored — Metal resource-scoped barriers
   are already as granular as the API allows. Forwards to the non-_ex path. */
void md_gpu_cmd_barrier_buffer_ex(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer,
                                   md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    md_gpu_cmd_barrier_buffer(cmd, buffer);
}

void md_gpu_cmd_barrier_image_ex(md_gpu_command_buffer_t cmd, md_gpu_image_t image,
                                  md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    md_gpu_cmd_barrier_image(cmd, image);
}

void md_gpu_cmd_copy_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t src,
                            md_gpu_buffer_t dst,
                            size_t size,
                            size_t src_offset,
                            size_t dst_offset) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_COPY_BUFFER;
    rc->u.copy_buf.src = src;
    rc->u.copy_buf.dst = dst;
    rc->u.copy_buf.size = size;
    rc->u.copy_buf.src_offset = src_offset;
    rc->u.copy_buf.dst_offset = dst_offset;
}

void md_gpu_cmd_copy_image_to_buffer(md_gpu_command_buffer_t cmd,
                                     md_gpu_image_t src_image,
                                     md_gpu_buffer_t dst_buffer) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_COPY_IMAGE_TO_BUFFER;
    rc->u.copy_img_buf.image = src_image;
    rc->u.copy_img_buf.buffer = dst_buffer;
}

void md_gpu_cmd_copy_buffer_to_image(md_gpu_command_buffer_t cmd,
                                     md_gpu_buffer_t src_buffer,
                                     md_gpu_image_t dst_image) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_COPY_BUFFER_TO_IMAGE;
    rc->u.copy_buf_img.buffer = src_buffer;
    rc->u.copy_buf_img.image = dst_image;
}

void md_gpu_cmd_copy_image_region_to_buffer(md_gpu_command_buffer_t cmd,
                                            md_gpu_image_t src_image,
                                            md_gpu_image_region_t src_region,
                                            md_gpu_buffer_t dst_buffer,
                                            size_t dst_offset) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_COPY_IMAGE_REGION_TO_BUFFER;
    rc->u.copy_img_region_buf.image      = src_image;
    rc->u.copy_img_region_buf.buffer     = dst_buffer;
    rc->u.copy_img_region_buf.region     = src_region;
    rc->u.copy_img_region_buf.dst_offset = dst_offset;
}

void md_gpu_cmd_copy_buffer_to_image_region(md_gpu_command_buffer_t cmd,
                                            md_gpu_buffer_t src_buffer,
                                            size_t src_offset,
                                            md_gpu_image_t dst_image,
                                            md_gpu_image_region_t dst_region) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_COPY_BUFFER_TO_IMAGE_REGION;
    rc->u.copy_buf_img_region.buffer     = src_buffer;
    rc->u.copy_buf_img_region.image      = dst_image;
    rc->u.copy_buf_img_region.region     = dst_region;
    rc->u.copy_buf_img_region.src_offset = src_offset;
}

void md_gpu_cmd_fill_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t buffer,
                            size_t offset,
                            size_t size,
                            uint8_t value) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_FILL_BUFFER;
    rc->u.fill_buf.buffer = buffer;
    rc->u.fill_buf.offset = offset;
    rc->u.fill_buf.size = size;
    rc->u.fill_buf.value = value;
}

/* =============================
   Submission & fences
   ============================= */

md_gpu_fence_t md_gpu_queue_submit(md_gpu_queue_t queue, md_gpu_command_buffer_t cmd) {

    id<MTLCommandBuffer> mtl_cmd = [queue->dev->queue commandBuffer];
    id<MTLComputeCommandEncoder> compute_enc = nil;

    for (uint32_t i = 0; i < cmd->command_count; ++i) {
        struct md_gpu_recorded_cmd* c = &cmd->commands[i];

        switch (c->type) {
            case CMD_DISPATCH: {
                if (!compute_enc)
                    compute_enc = [mtl_cmd computeCommandEncoder];

                [compute_enc setComputePipelineState:c->u.dispatch.pipeline->pso];

                for (uint32_t s = 0; s < MD_GPU_MAX_BIND_SLOTS; ++s) {
                    if (c->u.dispatch.buffers[s])
                        [compute_enc setBuffer:c->u.dispatch.buffers[s]->buffer offset:c->u.dispatch.buffer_offsets[s] atIndex:s + 1];
                    if (c->u.dispatch.images[s])
                        [compute_enc setTexture:c->u.dispatch.images[s]->texture atIndex:s];
                }

                if (c->u.dispatch.push_constant_size > 0)
                    [compute_enc setBytes:c->u.dispatch.push_constants
                                  length:c->u.dispatch.push_constant_size
                                 atIndex:0];

                struct md_gpu_compute_pipeline* p = c->u.dispatch.pipeline;
                MTLSize tg     = MTLSizeMake(p->tg_size[0], p->tg_size[1], p->tg_size[2]);
                MTLSize groups = MTLSizeMake(c->u.dispatch.group_count[0],
                                             c->u.dispatch.group_count[1],
                                             c->u.dispatch.group_count[2]);
                [compute_enc dispatchThreadgroups:groups threadsPerThreadgroup:tg];
                break;
            }

            case CMD_BARRIER:
                if (compute_enc)
                    [compute_enc memoryBarrierWithScope:(MTLBarrierScopeBuffers | MTLBarrierScopeTextures)];
                break;

            case CMD_BARRIER_BUFFER: {
                if (compute_enc) {
                    id<MTLResource> res = c->u.barrier_buf.buffer->buffer;
                    [compute_enc memoryBarrierWithResources:&res count:1];
                }
                break;
            }

            case CMD_BARRIER_IMAGE: {
                if (compute_enc) {
                    id<MTLResource> res = c->u.barrier_img.image->texture;
                    [compute_enc memoryBarrierWithResources:&res count:1];
                }
                break;
            }

            case CMD_COPY_BUFFER: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                [blit copyFromBuffer:c->u.copy_buf.src->buffer
                             sourceOffset:c->u.copy_buf.src_offset
                                 toBuffer:c->u.copy_buf.dst->buffer
                        destinationOffset:c->u.copy_buf.dst_offset
                                     size:c->u.copy_buf.size];
                [blit endEncoding];
                break;
            }

            case CMD_COPY_IMAGE_TO_BUFFER: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                MTLSize sz = MTLSizeMake(c->u.copy_img_buf.image->desc.width,
                                         c->u.copy_img_buf.image->desc.height,
                                         c->u.copy_img_buf.image->desc.depth);
                                const uint32_t bpp = md_gpu_format_bytes_per_pixel(c->u.copy_img_buf.image->desc.format);
                                ASSERT(bpp > 0);
                                const size_t bytes_per_row   = (size_t)sz.width * (size_t)bpp;
                                const size_t bytes_per_image = bytes_per_row * (size_t)sz.height;
                [blit copyFromTexture:c->u.copy_img_buf.image->texture
                           sourceSlice:0
                           sourceLevel:0
                          sourceOrigin:MTLOriginMake(0,0,0)
                            sourceSize:sz
                              toBuffer:c->u.copy_img_buf.buffer->buffer
                     destinationOffset:0
                                destinationBytesPerRow:bytes_per_row
                            destinationBytesPerImage:bytes_per_image];
                [blit endEncoding];
                break;
            }

            case CMD_COPY_BUFFER_TO_IMAGE: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                MTLSize sz = MTLSizeMake(c->u.copy_buf_img.image->desc.width,
                                         c->u.copy_buf_img.image->desc.height,
                                         c->u.copy_buf_img.image->desc.depth);
                const uint32_t bpp = md_gpu_format_bytes_per_pixel(c->u.copy_buf_img.image->desc.format);
                ASSERT(bpp > 0);
                const size_t bytes_per_row   = (size_t)sz.width * (size_t)bpp;
                const size_t bytes_per_image = bytes_per_row * (size_t)sz.height;
                [blit copyFromBuffer:c->u.copy_buf_img.buffer->buffer
                        sourceOffset:0
                   sourceBytesPerRow:bytes_per_row
                 sourceBytesPerImage:bytes_per_image
                          sourceSize:sz
                           toTexture:c->u.copy_buf_img.image->texture
                    destinationSlice:0
                    destinationLevel:0
                   destinationOrigin:MTLOriginMake(0,0,0)];
                [blit endEncoding];
                break;
            }

            case CMD_COPY_IMAGE_REGION_TO_BUFFER: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                const md_gpu_image_region_t* r = &c->u.copy_img_region_buf.region;
                MTLOrigin origin = MTLOriginMake(r->x, r->y, r->z);
                MTLSize   sz     = MTLSizeMake(r->width, r->height, r->depth);
                const uint32_t bpp = md_gpu_format_bytes_per_pixel(c->u.copy_img_region_buf.image->desc.format);
                ASSERT(bpp > 0);
                const size_t bytes_per_row   = (size_t)r->width * (size_t)bpp;
                const size_t bytes_per_image = bytes_per_row * (size_t)r->height;
                [blit copyFromTexture:c->u.copy_img_region_buf.image->texture
                           sourceSlice:0
                           sourceLevel:0
                          sourceOrigin:origin
                            sourceSize:sz
                              toBuffer:c->u.copy_img_region_buf.buffer->buffer
                     destinationOffset:c->u.copy_img_region_buf.dst_offset
                destinationBytesPerRow:bytes_per_row
              destinationBytesPerImage:bytes_per_image];
                [blit endEncoding];
                break;
            }

            case CMD_COPY_BUFFER_TO_IMAGE_REGION: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                const md_gpu_image_region_t* r = &c->u.copy_buf_img_region.region;
                MTLOrigin origin = MTLOriginMake(r->x, r->y, r->z);
                MTLSize   sz     = MTLSizeMake(r->width, r->height, r->depth);
                const uint32_t bpp = md_gpu_format_bytes_per_pixel(c->u.copy_buf_img_region.image->desc.format);
                ASSERT(bpp > 0);
                const size_t bytes_per_row   = (size_t)r->width * (size_t)bpp;
                const size_t bytes_per_image = bytes_per_row * (size_t)r->height;
                [blit copyFromBuffer:c->u.copy_buf_img_region.buffer->buffer
                        sourceOffset:c->u.copy_buf_img_region.src_offset
                   sourceBytesPerRow:bytes_per_row
                 sourceBytesPerImage:bytes_per_image
                          sourceSize:sz
                           toTexture:c->u.copy_buf_img_region.image->texture
                    destinationSlice:0
                    destinationLevel:0
                   destinationOrigin:origin];
                [blit endEncoding];
                break;
            }

            case CMD_FILL_BUFFER: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                [blit fillBuffer:c->u.fill_buf.buffer->buffer
                           range:NSMakeRange(c->u.fill_buf.offset, c->u.fill_buf.size)
                           value:c->u.fill_buf.value];
                [blit endEncoding];
                break;
            }
        }
    }

    if (compute_enc)
        [compute_enc endEncoding];

    [mtl_cmd commit];

    // Safe to recycle immediately: command translation/encoding is complete at commit.
    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;

    struct md_gpu_fence* fence = (struct md_gpu_fence*)calloc(1, sizeof(struct md_gpu_fence));
    fence->cmd = [mtl_cmd retain];

    return fence;
}

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence) {
    MTLCommandBufferStatus s = fence->cmd.status;
    return s == MTLCommandBufferStatusCompleted || s == MTLCommandBufferStatusError;
}

void md_gpu_fence_wait(md_gpu_fence_t fence) {
    [fence->cmd waitUntilCompleted];
}

void md_gpu_fence_destroy(md_gpu_fence_t fence) {
    [fence->cmd release];
    free(fence);
}

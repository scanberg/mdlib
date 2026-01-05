// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <core/md_common.h> // Assert
#include <stdlib.h>         // malloc, free
#include <pthread.h>        // pthread_mutex_t

/* Limits are defined in md_gpu.h */

// Metal-specific: maximum number of recorded commands per command buffer
#define MD_GPU_METAL_MAX_COMMANDS 1024

/* =============================
   Internal structs
   ============================= */

struct md_gpu_queue {
    struct md_gpu_device* dev;
    pthread_mutex_t submit_mutex;
};

struct md_gpu_device {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;

    struct md_gpu_queue compute_queue;

    pthread_mutex_t cmd_pool_mutex;
    struct md_gpu_command_buffer* cmd_free_list;
};

struct md_gpu_buffer {
    id<MTLBuffer> buffer;
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

struct md_gpu_cmd_fill_buffer {
    struct md_gpu_buffer* buffer;
    size_t offset;
    size_t size;
    uint32_t value;
};

struct md_gpu_cmd_barrier_buffer {
    struct md_gpu_buffer* buffer;
};

struct md_gpu_cmd_barrier_image {
    struct md_gpu_image* image;
};

struct md_gpu_recorded_cmd {
    md_gpu_cmd_type_t type;
    union {
        uint32_t dispatch_size[3];
        struct md_gpu_cmd_copy_buffer copy_buf;
        struct md_gpu_cmd_copy_image_to_buffer copy_img_buf;
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

md_gpu_device_t md_gpu_create_device(void) {
    struct md_gpu_device* dev = (struct md_gpu_device*)calloc(1, sizeof(struct md_gpu_device));
    dev->device = MTLCreateSystemDefaultDevice();
    ASSERT(dev->device);
    dev->queue = [dev->device newCommandQueue];

    dev->compute_queue.dev = dev;
    pthread_mutex_init(&dev->compute_queue.submit_mutex, NULL);

    pthread_mutex_init(&dev->cmd_pool_mutex, NULL);
    dev->cmd_free_list = NULL;
    return dev;
}

void md_gpu_destroy_device(md_gpu_device_t device) {
    pthread_mutex_destroy(&device->compute_queue.submit_mutex);

    pthread_mutex_lock(&device->cmd_pool_mutex);
    struct md_gpu_command_buffer* it = device->cmd_free_list;
    while (it) {
        struct md_gpu_command_buffer* next = it->next_free;
        free(it);
        it = next;
    }
    device->cmd_free_list = NULL;
    pthread_mutex_unlock(&device->cmd_pool_mutex);
    pthread_mutex_destroy(&device->cmd_pool_mutex);

    [device->queue release];
    [device->device release];
    free(device);
}

md_gpu_queue_t md_gpu_acquire_compute_queue(md_gpu_device_t device) {
    return &device->compute_queue;
}

/* =============================
   Buffers
   ============================= */

md_gpu_buffer_t md_gpu_create_buffer(md_gpu_device_t device,
                                     const md_gpu_buffer_desc_t* desc) {
    struct md_gpu_buffer* buf = (struct md_gpu_buffer*)calloc(1, sizeof(struct md_gpu_buffer));
    buf->size = desc->size;
    buf->flags = desc->flags;

    MTLResourceOptions opts = MTLResourceStorageModePrivate;
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE)
        opts = MTLResourceStorageModeShared;

    buf->buffer = [device->device newBufferWithLength:desc->size options:opts];
    return buf;
}

void md_gpu_destroy_buffer(md_gpu_buffer_t buffer) {
    [buffer->buffer release];
    free(buffer);
}

void* md_gpu_map_buffer(md_gpu_buffer_t buffer) {
    ASSERT(buffer->flags & MD_GPU_BUFFER_CPU_VISIBLE);
    return [buffer->buffer contents];
}

void md_gpu_unmap_buffer(md_gpu_buffer_t buffer) {
    (void)buffer;
}

/* =============================
   Images
   ============================= */

md_gpu_image_t md_gpu_create_image(md_gpu_device_t device,
                                   const md_gpu_image_desc_t* desc) {
    struct md_gpu_image* img = (struct md_gpu_image*)calloc(1, sizeof(struct md_gpu_image));
    img->desc = *desc;

    MTLTextureDescriptor* td = [[MTLTextureDescriptor alloc] init];
    td.textureType = desc->depth > 1 ? MTLTextureType3D : MTLTextureType2D;
    td.width  = desc->width;
    td.height = desc->height;
    td.depth  = desc->depth;
    td.pixelFormat = to_mtl_format(desc->format);
    td.usage = MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;

    img->texture = [device->device newTextureWithDescriptor:td];
    [td release];
    return img;
}

void md_gpu_destroy_image(md_gpu_image_t image) {
    [image->texture release];
    free(image);
}

/* =============================
   Pipelines
   ============================= */

md_gpu_compute_pipeline_t md_gpu_create_compute_pipeline(
    md_gpu_device_t device,
    const md_gpu_compute_pipeline_desc_t* desc)
{
    NSError* err = nil;

    NSData* data = [NSData dataWithBytes:desc->shader.data
                                  length:desc->shader.size];

    id<MTLLibrary> lib =
        [device->device newLibraryWithData:data error:&err];
    ASSERT(lib && !err);

    id<MTLFunction> fn = [lib newFunctionWithName:@"main0"];
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

void md_gpu_destroy_compute_pipeline(md_gpu_compute_pipeline_t pipeline) {
    [pipeline->pso release];
    free(pipeline);
}

/* =============================
   Command buffers
   ============================= */

md_gpu_command_buffer_t md_gpu_acquire_command_buffer(md_gpu_device_t device) {
    pthread_mutex_lock(&device->cmd_pool_mutex);
    struct md_gpu_command_buffer* cmd = device->cmd_free_list;
    if (cmd) {
        device->cmd_free_list = cmd->next_free;
    }
    pthread_mutex_unlock(&device->cmd_pool_mutex);

    if (!cmd) {
        cmd = (struct md_gpu_command_buffer*)calloc(1, sizeof(struct md_gpu_command_buffer));
        cmd->dev = device;
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
                         uint32_t x, uint32_t y, uint32_t z) {
    ASSERT(cmd->command_count < MD_GPU_MAX_COMMANDS);
    struct md_gpu_recorded_cmd* rc = &cmd->commands[cmd->command_count++];
    rc->type = CMD_DISPATCH;
    rc->u.dispatch_size[0] = x;
    rc->u.dispatch_size[1] = y;
    rc->u.dispatch_size[2] = z; 
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

void md_gpu_cmd_fill_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t buffer,
                            size_t offset,
                            size_t size,
                            uint32_t value) {
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

    pthread_mutex_lock(&queue->submit_mutex);

    id<MTLCommandBuffer> mtl_cmd = [queue->dev->queue commandBuffer];
    id<MTLComputeCommandEncoder> compute_enc = nil;

    for (uint32_t i = 0; i < cmd->command_count; ++i) {
        struct md_gpu_recorded_cmd* c = &cmd->commands[i];

        switch (c->type) {
            case CMD_DISPATCH: {
                if (!compute_enc)
                    compute_enc = [mtl_cmd computeCommandEncoder];

                [compute_enc setComputePipelineState:cmd->bound_pipeline->pso];

                for (uint32_t s = 0; s < MD_GPU_MAX_BIND_SLOTS; ++s) {
                    if (cmd->bound_buffers[s])
                        [compute_enc setBuffer:cmd->bound_buffers[s]->buffer offset:cmd->bound_buffer_offsets[s] atIndex:s];
                    if (cmd->bound_images[s])
                        [compute_enc setTexture:cmd->bound_images[s]->texture atIndex:s];
                }

                if (cmd->push_constant_size > 0)
                    [compute_enc setBytes:cmd->push_constants length:cmd->push_constant_size atIndex:MD_GPU_PUSH_CONSTANTS_SLOT];

                struct md_gpu_compute_pipeline* p = cmd->bound_pipeline;

                uint32_t gx = c->u.dispatch_size[0];
                uint32_t gy = c->u.dispatch_size[1];
                uint32_t gz = c->u.dispatch_size[2];
                MTLSize tg = MTLSizeMake(p->tg_size[0], p->tg_size[1], p->tg_size[2]);
                MTLSize grid = MTLSizeMake(
                    ((gx + tg.width  - 1) / tg.width)  * tg.width,
                    ((gy + tg.height - 1) / tg.height) * tg.height,
                    ((gz + tg.depth  - 1) / tg.depth)  * tg.depth
                );
                [compute_enc dispatchThreads:grid threadsPerThreadgroup:tg];
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

            case CMD_FILL_BUFFER: {
                if (compute_enc) {
                    [compute_enc endEncoding];
                    compute_enc = nil;
                }
                id<MTLBlitCommandEncoder> blit = [mtl_cmd blitCommandEncoder];
                uint8_t byte_val = (uint8_t)(c->u.fill_buf.value & 0xFF);
                [blit fillBuffer:c->u.fill_buf.buffer->buffer
                           range:NSMakeRange(c->u.fill_buf.offset, c->u.fill_buf.size)
                           value:byte_val];
                [blit endEncoding];
                break;
            }
        }
    }

    if (compute_enc)
        [compute_enc endEncoding];

    [mtl_cmd commit];
    pthread_mutex_unlock(&queue->submit_mutex);

    // Safe to recycle immediately: command translation/encoding is complete at commit.
    pthread_mutex_lock(&queue->dev->cmd_pool_mutex);
    cmd->next_free = queue->dev->cmd_free_list;
    queue->dev->cmd_free_list = cmd;
    pthread_mutex_unlock(&queue->dev->cmd_pool_mutex);

    struct md_gpu_fence* fence = (struct md_gpu_fence*)calloc(1, sizeof(struct md_gpu_fence));
    fence->cmd = [mtl_cmd retain];

    return fence;
}

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence) {
    return fence->cmd.status == MTLCommandBufferStatusCompleted;
}

void md_gpu_fence_wait(md_gpu_fence_t fence) {
    [fence->cmd waitUntilCompleted];
}

void md_gpu_destroy_fence(md_gpu_fence_t fence) {
    [fence->cmd release];
    free(fence);
}

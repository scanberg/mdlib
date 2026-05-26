// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <core/md_common.h> // Assert
#include <stdlib.h>         // malloc, free

/* Limits are defined in md_gpu.h */

#define MD_GPU_PASS_HISTORY_SIZE 256

struct md_gpu_pass_record {
    uint64_t id;
    id<MTLCommandBuffer> cmd;
};

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

    uint64_t next_pass_id;
    uint32_t pass_record_head;
    struct md_gpu_pass_record pass_records[MD_GPU_PASS_HISTORY_SIZE];
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

struct md_gpu_sampler {
    id<MTLSamplerState> sampler;
};

struct md_gpu_compute_pipeline {
    id<MTLComputePipelineState> pso;
    uint32_t tg_size[3];
};

struct md_gpu_command_buffer {
    struct md_gpu_device* dev;

    /* Live Metal encoding state */
    id<MTLCommandBuffer>         mtl_cmd;
    id<MTLComputeCommandEncoder> compute_enc;
    id<MTLBlitCommandEncoder>    blit_enc;

    /* Live bind state (updated by md_gpu_cmd_bind_* and md_gpu_cmd_push_constants) */
    struct md_gpu_compute_pipeline* bound_pipeline;
    struct md_gpu_buffer* bound_buffers[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_offsets[MD_GPU_MAX_BIND_SLOTS];
    struct md_gpu_image*  bound_images[MD_GPU_MAX_BIND_SLOTS];
    struct md_gpu_image*  bound_sampled_images[MD_GPU_MAX_BIND_SLOTS];
    struct md_gpu_sampler* bound_samplers[MD_GPU_MAX_BIND_SLOTS];
    uint8_t push_constants[MD_GPU_MAX_PUSH_CONSTANTS];
    size_t  push_constant_size;

    struct md_gpu_command_buffer* next_free;
};

struct md_gpu_fence {
    id<MTLCommandBuffer> cmd;
};

struct md_gpu_pass {
    struct md_gpu_queue* queue;
    struct md_gpu_command_buffer* cmd;
    bool pushed_debug_group;
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

static MTLSamplerMinMagFilter to_mtl_filter(md_gpu_filter_t filter) {
    switch (filter) {
        case MD_GPU_FILTER_NEAREST: return MTLSamplerMinMagFilterNearest;
        case MD_GPU_FILTER_LINEAR:  return MTLSamplerMinMagFilterLinear;
        default: return MTLSamplerMinMagFilterNearest;
    }
}

static MTLSamplerAddressMode to_mtl_address_mode(md_gpu_address_mode_t mode) {
    switch (mode) {
        case MD_GPU_ADDRESS_MODE_CLAMP_TO_EDGE:   return MTLSamplerAddressModeClampToEdge;
        case MD_GPU_ADDRESS_MODE_REPEAT:          return MTLSamplerAddressModeRepeat;
        case MD_GPU_ADDRESS_MODE_MIRRORED_REPEAT: return MTLSamplerAddressModeMirrorRepeat;
        default: return MTLSamplerAddressModeClampToEdge;
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
    if (cmd->compute_enc) {
        [cmd->compute_enc endEncoding];
        [cmd->compute_enc release];
        cmd->compute_enc = nil;
    }
    if (cmd->blit_enc) {
        [cmd->blit_enc endEncoding];
        [cmd->blit_enc release];
        cmd->blit_enc = nil;
    }
    if (cmd->mtl_cmd) {
        [cmd->mtl_cmd release];
        cmd->mtl_cmd = nil;
    }
    cmd->bound_pipeline = NULL;
    MEMSET(cmd->bound_buffers, 0, sizeof(cmd->bound_buffers));
    MEMSET(cmd->bound_buffer_offsets, 0, sizeof(cmd->bound_buffer_offsets));
    MEMSET(cmd->bound_images,  0, sizeof(cmd->bound_images));
    MEMSET(cmd->bound_sampled_images, 0, sizeof(cmd->bound_sampled_images));
    MEMSET(cmd->bound_samplers, 0, sizeof(cmd->bound_samplers));
    cmd->push_constant_size = 0;
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
    for (uint32_t i = 0; i < MD_GPU_PASS_HISTORY_SIZE; ++i) {
        if (device->pass_records[i].cmd) {
            [device->pass_records[i].cmd waitUntilCompleted];
            [device->pass_records[i].cmd release];
            device->pass_records[i].cmd = nil;
        }
        device->pass_records[i].id = 0;
    }

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

bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* info) {
    if (!device || !info) return false;

    MEMSET(info, 0, sizeof(*info));

    // Metal targets Apple Silicon and Intel iGPUs — always UMA, never discrete.
    info->is_discrete = false;

    NSString* name = [device->device name];
    if (name) {
        [name getCString:info->name maxLength:sizeof(info->name) encoding:NSUTF8StringEncoding];
    }

    return true;
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
    if (!buffer) return;
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
    if (!image) return;
    [image->texture release];
    free(image);
}

bool md_gpu_image_desc_extract(md_gpu_image_t image, md_gpu_image_desc_t* out_desc) {
    if (!image || !out_desc) return false;
    *out_desc = image->desc;
    return true;
}

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t device,
                                       const md_gpu_sampler_desc_t* desc) {
    if (!device || !desc) return NULL;
    struct md_gpu_sampler* sampler = (struct md_gpu_sampler*)calloc(1, sizeof(struct md_gpu_sampler));
    if (!sampler) return NULL;

    MTLSamplerDescriptor* sd = [[MTLSamplerDescriptor alloc] init];
    sd.minFilter = to_mtl_filter(desc->min_filter);
    sd.magFilter = to_mtl_filter(desc->mag_filter);
    sd.mipFilter = MTLSamplerMipFilterNotMipmapped;
    sd.sAddressMode = to_mtl_address_mode(desc->address_u);
    sd.tAddressMode = to_mtl_address_mode(desc->address_v);
    sd.rAddressMode = to_mtl_address_mode(desc->address_w);
    sd.normalizedCoordinates = YES;

    sampler->sampler = [device->device newSamplerStateWithDescriptor:sd];
    [sd release];
    if (!sampler->sampler) {
        free(sampler);
        return NULL;
    }

    return sampler;
}

void md_gpu_sampler_destroy(md_gpu_sampler_t sampler) {
    if (!sampler) return;
    [sampler->sampler release];
    free(sampler);
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
    if (!pipeline) return;
    [pipeline->pso release];
    free(pipeline);
}

/* =============================
   Encoder helpers
   ============================= */

static inline id<MTLComputeCommandEncoder> ensure_compute_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->blit_enc) {
        [cmd->blit_enc endEncoding];
        [cmd->blit_enc release];
        cmd->blit_enc = nil;
    }
    if (!cmd->compute_enc)
        cmd->compute_enc = [[cmd->mtl_cmd computeCommandEncoder] retain];
    return cmd->compute_enc;
}

static inline void end_compute_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->compute_enc) {
        [cmd->compute_enc endEncoding];
        [cmd->compute_enc release];
        cmd->compute_enc = nil;
    }
}

static inline id<MTLBlitCommandEncoder> ensure_blit_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->compute_enc) {
        [cmd->compute_enc endEncoding];
        [cmd->compute_enc release];
        cmd->compute_enc = nil;
    }
    if (!cmd->blit_enc)
        cmd->blit_enc = [[cmd->mtl_cmd blitCommandEncoder] retain];
    return cmd->blit_enc;
}

static inline void end_blit_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->blit_enc) {
        [cmd->blit_enc endEncoding];
        [cmd->blit_enc release];
        cmd->blit_enc = nil;
    }
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
    cmd->mtl_cmd = [[queue->dev->queue commandBuffer] retain];
    return cmd;
}

md_gpu_device_t md_gpu_command_buffer_device(md_gpu_command_buffer_t cmd) {
    return cmd ? cmd->dev : NULL;
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

void md_gpu_cmd_bind_sampled_image(md_gpu_command_buffer_t cmd,
                                   uint32_t slot,
                                   md_gpu_image_t image) {
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    ASSERT(slot != MD_GPU_PUSH_CONSTANTS_SLOT);
    cmd->bound_sampled_images[slot] = image;
}

void md_gpu_cmd_bind_sampler(md_gpu_command_buffer_t cmd,
                             uint32_t slot,
                             md_gpu_sampler_t sampler) {
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    ASSERT(slot != MD_GPU_PUSH_CONSTANTS_SLOT);
    cmd->bound_samplers[slot] = sampler;
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
    ASSERT(cmd->bound_pipeline);
    id<MTLComputeCommandEncoder> enc = ensure_compute_enc(cmd);

    [enc setComputePipelineState:cmd->bound_pipeline->pso];

    for (uint32_t s = 0; s < MD_GPU_MAX_BIND_SLOTS; ++s) {
        if (cmd->bound_buffers[s])
            [enc setBuffer:cmd->bound_buffers[s]->buffer offset:cmd->bound_buffer_offsets[s] atIndex:s + 1];
        if (cmd->bound_images[s])
            [enc setTexture:cmd->bound_images[s]->texture atIndex:s];
        if (cmd->bound_sampled_images[s])
            [enc setTexture:cmd->bound_sampled_images[s]->texture atIndex:s];
        if (cmd->bound_samplers[s])
            [enc setSamplerState:cmd->bound_samplers[s]->sampler atIndex:s];
    }

    if (cmd->push_constant_size > 0)
        [enc setBytes:cmd->push_constants length:cmd->push_constant_size atIndex:0];

    MTLSize tg     = MTLSizeMake(cmd->bound_pipeline->tg_size[0], cmd->bound_pipeline->tg_size[1], cmd->bound_pipeline->tg_size[2]);
    MTLSize groups = MTLSizeMake(group_count_x, group_count_y, group_count_z);
    [enc dispatchThreadgroups:groups threadsPerThreadgroup:tg];
}

void md_gpu_cmd_barrier(md_gpu_command_buffer_t cmd, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    if (cmd->compute_enc)
        [cmd->compute_enc memoryBarrierWithScope:(MTLBarrierScopeBuffers | MTLBarrierScopeTextures)];
}

void md_gpu_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    ASSERT(buffer);
    if (cmd->compute_enc) {
        id<MTLResource> res = buffer->buffer;
        [cmd->compute_enc memoryBarrierWithResources:&res count:1];
    }
}

void md_gpu_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image,
                              md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    ASSERT(image);
    if (cmd->compute_enc) {
        id<MTLResource> res = image->texture;
        [cmd->compute_enc memoryBarrierWithResources:&res count:1];
    }
}

void md_gpu_cmd_copy_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t src,
                            md_gpu_buffer_t dst,
                            size_t size,
                            size_t src_offset,
                            size_t dst_offset) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit copyFromBuffer:src->buffer sourceOffset:src_offset toBuffer:dst->buffer destinationOffset:dst_offset size:size];
}

void md_gpu_cmd_copy_buffer_to_image(md_gpu_command_buffer_t cmd,
                                     md_gpu_buffer_t src_buffer,
                                     md_gpu_image_t dst_image) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    MTLSize sz = MTLSizeMake(dst_image->desc.width, dst_image->desc.height, dst_image->desc.depth);
    const uint32_t bpp = md_gpu_format_bytes_per_pixel(dst_image->desc.format);
    ASSERT(bpp > 0);
    const size_t bytes_per_row   = (size_t)sz.width * (size_t)bpp;
    const size_t bytes_per_image = bytes_per_row * (size_t)sz.height;
    [blit copyFromBuffer:src_buffer->buffer
            sourceOffset:0
       sourceBytesPerRow:bytes_per_row
     sourceBytesPerImage:bytes_per_image
              sourceSize:sz
               toTexture:dst_image->texture
        destinationSlice:0
        destinationLevel:0
       destinationOrigin:MTLOriginMake(0,0,0)];
}

void md_gpu_cmd_copy_image_region_to_buffer(md_gpu_command_buffer_t cmd,
                                            md_gpu_image_t src_image,
                                            md_gpu_image_region_t src_region,
                                            md_gpu_buffer_t dst_buffer,
                                            size_t dst_offset) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    MTLOrigin origin = MTLOriginMake(src_region.offset[0], src_region.offset[1], src_region.offset[2]);
    uint32_t ew = src_region.extent[0] ? src_region.extent[0] : src_image->desc.width;
    uint32_t eh = src_region.extent[1] ? src_region.extent[1] : src_image->desc.height;
    uint32_t ed = src_region.extent[2] ? src_region.extent[2] : src_image->desc.depth;
    MTLSize sz = MTLSizeMake(ew, eh, ed);
    const uint32_t bpp = md_gpu_format_bytes_per_pixel(src_image->desc.format);
    ASSERT(bpp > 0);
    const size_t bytes_per_row   = (size_t)ew * (size_t)bpp;
    const size_t bytes_per_image = bytes_per_row * (size_t)eh;
    [blit copyFromTexture:src_image->texture
               sourceSlice:0
               sourceLevel:0
              sourceOrigin:origin
                sourceSize:sz
                  toBuffer:dst_buffer->buffer
         destinationOffset:dst_offset
    destinationBytesPerRow:bytes_per_row
  destinationBytesPerImage:bytes_per_image];
}

void md_gpu_cmd_copy_buffer_to_image_region(md_gpu_command_buffer_t cmd,
                                            md_gpu_buffer_t src_buffer,
                                            size_t src_offset,
                                            md_gpu_image_t dst_image,
                                            md_gpu_image_region_t dst_region) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    uint32_t ew = dst_region.extent[0] ? dst_region.extent[0] : dst_image->desc.width;
    uint32_t eh = dst_region.extent[1] ? dst_region.extent[1] : dst_image->desc.height;
    uint32_t ed = dst_region.extent[2] ? dst_region.extent[2] : dst_image->desc.depth;
    MTLOrigin origin = MTLOriginMake(dst_region.offset[0], dst_region.offset[1], dst_region.offset[2]);
    MTLSize   sz     = MTLSizeMake(ew, eh, ed);
    const uint32_t bpp = md_gpu_format_bytes_per_pixel(dst_image->desc.format);
    ASSERT(bpp > 0);
    const size_t bytes_per_row   = (size_t)ew * (size_t)bpp;
    const size_t bytes_per_image = bytes_per_row * (size_t)eh;
    [blit copyFromBuffer:src_buffer->buffer
            sourceOffset:src_offset
       sourceBytesPerRow:bytes_per_row
     sourceBytesPerImage:bytes_per_image
              sourceSize:sz
               toTexture:dst_image->texture
        destinationSlice:0
        destinationLevel:0
       destinationOrigin:origin];
}

void md_gpu_cmd_fill_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t buffer,
                            size_t offset,
                            size_t size,
                            uint8_t value) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit fillBuffer:buffer->buffer range:NSMakeRange(offset, size) value:value];
}

void md_gpu_cmd_push_debug_group(md_gpu_command_buffer_t cmd, const char* label) {
    if (!cmd || !label) return;
    NSString* ns_label = [NSString stringWithUTF8String:label];
    if (cmd->compute_enc)
        [cmd->compute_enc pushDebugGroup:ns_label];
    else if (cmd->blit_enc)
        [cmd->blit_enc pushDebugGroup:ns_label];
    else
        [cmd->mtl_cmd pushDebugGroup:ns_label];
}

void md_gpu_cmd_pop_debug_group(md_gpu_command_buffer_t cmd) {
    if (!cmd) return;
    if (cmd->compute_enc)
        [cmd->compute_enc popDebugGroup];
    else if (cmd->blit_enc)
        [cmd->blit_enc popDebugGroup];
    else
        [cmd->mtl_cmd popDebugGroup];
}

bool md_gpu_queue_submit(md_gpu_queue_t queue, md_gpu_command_buffer_t cmd, md_gpu_fence_t fence) {
    ASSERT(cmd && cmd->mtl_cmd);

    end_compute_enc(cmd);
    end_blit_enc(cmd);

    [cmd->mtl_cmd commit];

    if (fence) {
        fence->cmd = [cmd->mtl_cmd retain];
    } else {
        [cmd->mtl_cmd waitUntilCompleted];
    }

    [cmd->mtl_cmd release];
    cmd->mtl_cmd = nil;

    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;

    return true;
}

md_gpu_fence_t md_gpu_fence_create(md_gpu_device_t device) {
    (void)device;
    struct md_gpu_fence* f = (struct md_gpu_fence*)calloc(1, sizeof(struct md_gpu_fence));
    return f; // cmd is set by md_gpu_queue_submit
}

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence) {
    MTLCommandBufferStatus s = fence->cmd.status;
    return s == MTLCommandBufferStatusCompleted || s == MTLCommandBufferStatusError;
}

void md_gpu_fence_wait(md_gpu_fence_t fence) {
    [fence->cmd waitUntilCompleted];
}

void md_gpu_fence_destroy(md_gpu_fence_t fence) {
    if (!fence) return;
    if (fence->cmd) [fence->cmd release];
    free(fence);
}

static id<MTLCommandBuffer> md_mtl_queue_submit_pass(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd) {
    ASSERT(queue && cmd && cmd->mtl_cmd);

    end_compute_enc(cmd);
    end_blit_enc(cmd);

    [cmd->mtl_cmd commit];
    id<MTLCommandBuffer> submitted = [cmd->mtl_cmd retain];

    [cmd->mtl_cmd release];
    cmd->mtl_cmd = nil;

    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;

    return submitted;
}

static struct md_gpu_pass_record* md_mtl_pass_record_find(struct md_gpu_device* dev, md_gpu_pass_id_t id) {
    if (!dev || id.value == 0) return NULL;
    for (uint32_t i = 0; i < MD_GPU_PASS_HISTORY_SIZE; ++i) {
        if (dev->pass_records[i].id == id.value) return &dev->pass_records[i];
    }
    return NULL;
}

md_gpu_pass_t md_gpu_pass_begin(md_gpu_device_t device, const md_gpu_pass_desc_t* desc) {
    if (!device) return NULL;
    struct md_gpu_device* dev = device;
    struct md_gpu_queue* queue = &dev->compute_queue;

    struct md_gpu_pass* pass = (struct md_gpu_pass*)calloc(1, sizeof(struct md_gpu_pass));
    if (!pass) return NULL;

    pass->queue = queue;
    pass->cmd = md_gpu_command_buffer_acquire((md_gpu_queue_t)queue);
    if (!pass->cmd) {
        free(pass);
        return NULL;
    }

    if (desc && desc->label) {
        md_gpu_cmd_push_debug_group(pass->cmd, desc->label);
        pass->pushed_debug_group = true;
    }
    return pass;
}

md_gpu_device_t md_gpu_pass_device(md_gpu_pass_t pass) {
    if (!pass) return NULL;
    return pass->queue->dev;
}

md_gpu_pass_id_t md_gpu_pass_end(md_gpu_pass_t pass) {
    md_gpu_pass_id_t id = {0};
    if (!pass) return id;

    struct md_gpu_queue* queue = pass->queue;
    struct md_gpu_device* dev = queue ? queue->dev : NULL;
    if (!queue || !dev || !pass->cmd) {
        free(pass);
        return id;
    }

    if (pass->pushed_debug_group) {
        md_gpu_cmd_pop_debug_group(pass->cmd);
    }

    struct md_gpu_pass_record* record = &dev->pass_records[dev->pass_record_head++ % MD_GPU_PASS_HISTORY_SIZE];
    if (record->cmd) {
        [record->cmd waitUntilCompleted];
        [record->cmd release];
    }
    record->id = 0;
    record->cmd = nil;

    id.value = ++dev->next_pass_id;
    if (id.value == 0) id.value = ++dev->next_pass_id;

    id<MTLCommandBuffer> submitted = md_mtl_queue_submit_pass(queue, pass->cmd);
    if (submitted) {
        record->id = id.value;
        record->cmd = submitted;
    } else {
        id.value = 0;
    }

    free(pass);
    return id;
}

bool md_gpu_pass_is_complete(md_gpu_device_t device, md_gpu_pass_id_t id) {
    if (!device || id.value == 0) return true;
    struct md_gpu_pass_record* record = md_mtl_pass_record_find(device, id);
    if (!record || !record->cmd) return true;
    MTLCommandBufferStatus s = record->cmd.status;
    return s == MTLCommandBufferStatusCompleted || s == MTLCommandBufferStatusError;
}

void md_gpu_pass_wait(md_gpu_device_t device, md_gpu_pass_id_t id) {
    if (!device || id.value == 0) return;
    struct md_gpu_pass_record* record = md_mtl_pass_record_find(device, id);
    if (!record || !record->cmd) return;

    [record->cmd waitUntilCompleted];
    [record->cmd release];
    record->cmd = nil;
}

void md_gpu_pass_wait_pass(md_gpu_pass_t pass, md_gpu_pass_id_t id) {
    (void)pass;
    (void)id;
}

void md_gpu_pass_bind_compute_pipeline(md_gpu_pass_t pass, md_gpu_compute_pipeline_t pipeline) {
    if (!pass) return;
    md_gpu_cmd_bind_compute_pipeline(pass->cmd, pipeline);
}

void md_gpu_pass_bind_buffer(md_gpu_pass_t pass, uint32_t slot, md_gpu_buffer_t buffer) {
    if (!pass) return;
    md_gpu_cmd_bind_buffer(pass->cmd, slot, buffer);
}

void md_gpu_pass_bind_buffer_range(md_gpu_pass_t pass, uint32_t slot, md_gpu_buffer_t buffer, size_t offset, size_t size) {
    if (!pass) return;
    md_gpu_cmd_bind_buffer_range(pass->cmd, slot, buffer, offset, size);
}

void md_gpu_pass_bind_image(md_gpu_pass_t pass, uint32_t slot, md_gpu_image_t image) {
    if (!pass) return;
    md_gpu_cmd_bind_image(pass->cmd, slot, image);
}

void md_gpu_pass_bind_sampled_image(md_gpu_pass_t pass, uint32_t slot, md_gpu_image_t image) {
    if (!pass) return;
    md_gpu_cmd_bind_sampled_image(pass->cmd, slot, image);
}

void md_gpu_pass_bind_sampler(md_gpu_pass_t pass, uint32_t slot, md_gpu_sampler_t sampler) {
    if (!pass) return;
    md_gpu_cmd_bind_sampler(pass->cmd, slot, sampler);
}

void md_gpu_pass_push_constants(md_gpu_pass_t pass, const void* data, size_t size) {
    if (!pass) return;
    md_gpu_cmd_push_constants(pass->cmd, data, size);
}

void md_gpu_pass_dispatch(md_gpu_pass_t pass, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z) {
    if (!pass) return;
    md_gpu_cmd_dispatch(pass->cmd, group_count_x, group_count_y, group_count_z);
}

void md_gpu_pass_barrier(md_gpu_pass_t pass, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass) return;
    md_gpu_cmd_barrier(pass->cmd, src_stage, dst_stage);
}

void md_gpu_pass_barrier_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer,
                                md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass) return;
    md_gpu_cmd_barrier_buffer(pass->cmd, buffer, src_stage, dst_stage);
}

void md_gpu_pass_barrier_image(md_gpu_pass_t pass, md_gpu_image_t image,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass) return;
    md_gpu_cmd_barrier_image(pass->cmd, image, src_stage, dst_stage);
}

void md_gpu_pass_push_debug_group(md_gpu_pass_t pass, const char* label) {
    if (!pass) return;
    md_gpu_cmd_push_debug_group(pass->cmd, label);
}

void md_gpu_pass_pop_debug_group(md_gpu_pass_t pass) {
    if (!pass) return;
    md_gpu_cmd_pop_debug_group(pass->cmd);
}

void md_gpu_pass_copy_buffer(md_gpu_pass_t pass, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!pass) return;
    md_gpu_cmd_copy_buffer(pass->cmd, src, dst, size, src_offset, dst_offset);
}

void md_gpu_pass_copy_image_region_to_buffer(md_gpu_pass_t pass,
                                             md_gpu_image_t src_image,
                                             md_gpu_image_region_t src_region,
                                             md_gpu_buffer_t dst_buffer,
                                             size_t dst_offset) {
    if (!pass) return;
    md_gpu_cmd_copy_image_region_to_buffer(pass->cmd, src_image, src_region, dst_buffer, dst_offset);
}

void md_gpu_pass_copy_buffer_to_image(md_gpu_pass_t pass, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image) {
    if (!pass) return;
    md_gpu_cmd_copy_buffer_to_image(pass->cmd, src_buffer, dst_image);
}

void md_gpu_pass_copy_buffer_to_image_region(md_gpu_pass_t pass,
                                             md_gpu_buffer_t src_buffer,
                                             size_t src_offset,
                                             md_gpu_image_t dst_image,
                                             md_gpu_image_region_t dst_region) {
    if (!pass) return;
    md_gpu_cmd_copy_buffer_to_image_region(pass->cmd, src_buffer, src_offset, dst_image, dst_region);
}

void md_gpu_pass_fill_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer, size_t offset, size_t size, uint8_t value) {
    if (!pass) return;
    md_gpu_cmd_fill_buffer(pass->cmd, buffer, offset, size, value);
}

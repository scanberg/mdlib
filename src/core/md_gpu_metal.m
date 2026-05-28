// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <core/md_common.h> // Assert
#include <stdlib.h>         // malloc, free

/* Limits are defined in md_gpu.h */

#define MD_GPU_EVENT_HISTORY_SIZE 256
#define MD_GPU_MAX_BIND_SLOTS 16
#define MD_GPU_MAX_PUSH_CONSTANTS 256

typedef struct md_gpu_queue* md_gpu_queue_t;
typedef struct md_gpu_command_buffer* md_gpu_command_buffer_t;

typedef struct md_gpu_cmd_buffer_usage_t {
    md_gpu_buffer_t buffer;
    gpu_usage_flags_t usage;
} md_gpu_cmd_buffer_usage_t;

typedef struct md_gpu_cmd_image_usage_t {
    md_gpu_image_t image;
    gpu_usage_flags_t usage;
    md_gpu_image_region_t region;
} md_gpu_cmd_image_usage_t;

typedef struct md_gpu_cmd_sampled_resource_usage_t {
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
    gpu_usage_flags_t usage;
} md_gpu_cmd_sampled_resource_usage_t;

struct md_gpu_event_record {
    uint64_t value;
    id<MTLCommandBuffer> mtl_cmd;
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

    uint64_t next_event_id;
    uint32_t event_record_head;
    struct md_gpu_event_record event_records[MD_GPU_EVENT_HISTORY_SIZE];
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

    md_gpu_cmd_buffer_usage_t cmd_buffers[MD_GPU_MAX_BIND_SLOTS];
    uint32_t cmd_buffer_count;
    md_gpu_cmd_image_usage_t cmd_images[MD_GPU_MAX_BIND_SLOTS];
    uint32_t cmd_image_count;
    md_gpu_cmd_sampled_resource_usage_t cmd_sampled_resources[MD_GPU_MAX_BIND_SLOTS];
    uint32_t cmd_sampled_resource_count;

    struct md_gpu_command_buffer* next_free;
};

struct md_gpu_cmd {
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

static bool mtl_cmd_merge_buffer_usage(struct md_gpu_command_buffer* cmd, md_gpu_buffer_t buffer, gpu_usage_flags_t usage) {
    if (!cmd || !buffer || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->cmd_buffer_count; ++i) {
        if (cmd->cmd_buffers[i].buffer == buffer) {
            cmd->cmd_buffers[i].usage |= usage;
            return true;
        }
    }

    if (cmd->cmd_buffer_count >= MD_GPU_MAX_BIND_SLOTS) return false;
    cmd->cmd_buffers[cmd->cmd_buffer_count].buffer = buffer;
    cmd->cmd_buffers[cmd->cmd_buffer_count].usage = usage;
    cmd->cmd_buffer_count++;
    return true;
}

static bool mtl_cmd_merge_image_usage(struct md_gpu_command_buffer* cmd, md_gpu_image_t image, gpu_usage_flags_t usage) {
    if (!cmd || !image || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->cmd_image_count; ++i) {
        if (cmd->cmd_images[i].image == image) {
            cmd->cmd_images[i].usage |= usage;
            return true;
        }
    }

    if (cmd->cmd_image_count >= MD_GPU_MAX_BIND_SLOTS) return false;
    cmd->cmd_images[cmd->cmd_image_count].image = image;
    cmd->cmd_images[cmd->cmd_image_count].usage = usage;
    MEMSET(&cmd->cmd_images[cmd->cmd_image_count].region, 0, sizeof(cmd->cmd_images[cmd->cmd_image_count].region));
    cmd->cmd_image_count++;
    return true;
}

static bool mtl_cmd_merge_sampled_resource_usage(struct md_gpu_command_buffer* cmd, md_gpu_image_t image, md_gpu_sampler_t sampler, gpu_usage_flags_t usage) {
    if (!cmd || !image || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->cmd_sampled_resource_count; ++i) {
        if (cmd->cmd_sampled_resources[i].image == image && cmd->cmd_sampled_resources[i].sampler == sampler) {
            cmd->cmd_sampled_resources[i].usage |= usage;
            return true;
        }
    }

    if (cmd->cmd_sampled_resource_count >= MD_GPU_MAX_BIND_SLOTS) return false;
    cmd->cmd_sampled_resources[cmd->cmd_sampled_resource_count].image = image;
    cmd->cmd_sampled_resources[cmd->cmd_sampled_resource_count].sampler = sampler;
    cmd->cmd_sampled_resources[cmd->cmd_sampled_resource_count].usage = usage;
    cmd->cmd_sampled_resource_count++;
    return true;
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
    cmd->cmd_buffer_count = 0;
    cmd->cmd_image_count = 0;
    cmd->cmd_sampled_resource_count = 0;
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
    for (uint32_t i = 0; i < MD_GPU_EVENT_HISTORY_SIZE; ++i) {
        if (device->event_records[i].mtl_cmd) {
            [device->event_records[i].mtl_cmd waitUntilCompleted];
            [device->event_records[i].mtl_cmd release];
            device->event_records[i].mtl_cmd = nil;
        }
        device->event_records[i].value = 0;
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

static inline MTLResourceUsage mtl_resource_usage(gpu_usage_flags_t usage) {
    MTLResourceUsage mtl_usage = 0;
    if (usage & GPU_USAGE_READ) {
        mtl_usage |= MTLResourceUsageRead;
    }
    if (usage & GPU_USAGE_WRITE) {
        mtl_usage |= MTLResourceUsageWrite;
    }
    return mtl_usage ? mtl_usage : MTLResourceUsageRead;
}

static inline void use_cmd_resources(struct md_gpu_command_buffer* cmd, id<MTLComputeCommandEncoder> enc) {
    for (uint32_t i = 0; i < cmd->cmd_buffer_count; ++i) {
        md_gpu_cmd_buffer_usage_t usage = cmd->cmd_buffers[i];
        if (usage.buffer) {
            [enc useResource:usage.buffer->buffer usage:mtl_resource_usage(usage.usage)];
        }
    }
    for (uint32_t i = 0; i < cmd->cmd_image_count; ++i) {
        md_gpu_cmd_image_usage_t usage = cmd->cmd_images[i];
        if (usage.image) {
            [enc useResource:usage.image->texture usage:mtl_resource_usage(usage.usage)];
        }
    }
    for (uint32_t i = 0; i < cmd->cmd_sampled_resource_count; ++i) {
        md_gpu_cmd_sampled_resource_usage_t usage = cmd->cmd_sampled_resources[i];
        if (usage.image) {
            [enc useResource:usage.image->texture usage:MTLResourceUsageRead];
        }
    }
}

static inline void apply_cmd_declared_bindings(struct md_gpu_command_buffer* cmd) {
    MEMSET(cmd->bound_images, 0, sizeof(cmd->bound_images));
    MEMSET(cmd->bound_sampled_images, 0, sizeof(cmd->bound_sampled_images));
    MEMSET(cmd->bound_samplers, 0, sizeof(cmd->bound_samplers));

    const uint32_t image_count = cmd->cmd_image_count < MD_GPU_MAX_BIND_SLOTS ? cmd->cmd_image_count : MD_GPU_MAX_BIND_SLOTS;
    for (uint32_t i = 0; i < image_count; ++i) {
        cmd->bound_images[i] = cmd->cmd_images[i].image;
    }

    const uint32_t sampled_resource_count = cmd->cmd_sampled_resource_count < MD_GPU_MAX_BIND_SLOTS ? cmd->cmd_sampled_resource_count : MD_GPU_MAX_BIND_SLOTS;
    for (uint32_t i = 0; i < sampled_resource_count; ++i) {
        cmd->bound_sampled_images[i] = cmd->cmd_sampled_resources[i].image;
        cmd->bound_samplers[i] = cmd->cmd_sampled_resources[i].sampler;
    }
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

static md_gpu_command_buffer_t md_gpu_command_buffer_acquire(md_gpu_queue_t queue) {
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

static void mtl_cmd_bind_compute_pipeline(md_gpu_command_buffer_t cmd,
                                     md_gpu_compute_pipeline_t pipeline) {
    cmd->bound_pipeline = pipeline;
}

static void mtl_cmd_dispatch(md_gpu_command_buffer_t cmd,
                         uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z) {
    ASSERT(cmd->bound_pipeline);
    id<MTLComputeCommandEncoder> enc = ensure_compute_enc(cmd);

    [enc setComputePipelineState:cmd->bound_pipeline->pso];
    apply_cmd_declared_bindings(cmd);
    use_cmd_resources(cmd, enc);

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

    /* Reset per-dispatch resource declarations so the next md_gpu_cmd_use_* calls
       start from slot 0. Resource hazard tracking for useResource: is re-issued per dispatch. */
    cmd->cmd_buffer_count = 0;
    cmd->cmd_image_count = 0;
    cmd->cmd_sampled_resource_count = 0;
}

static void mtl_cmd_barrier(md_gpu_command_buffer_t cmd, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    if (cmd->compute_enc)
        [cmd->compute_enc memoryBarrierWithScope:(MTLBarrierScopeBuffers | MTLBarrierScopeTextures)];
}

static void mtl_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    ASSERT(buffer);
    if (cmd->compute_enc) {
        id<MTLResource> res = buffer->buffer;
        [cmd->compute_enc memoryBarrierWithResources:&res count:1];
    }
}

static void mtl_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image,
                              md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    (void)src_stage; (void)dst_stage;
    ASSERT(image);
    if (cmd->compute_enc) {
        id<MTLResource> res = image->texture;
        [cmd->compute_enc memoryBarrierWithResources:&res count:1];
    }
}

static void mtl_cmd_copy_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t src,
                            md_gpu_buffer_t dst,
                            size_t size,
                            size_t src_offset,
                            size_t dst_offset) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit copyFromBuffer:src->buffer sourceOffset:src_offset toBuffer:dst->buffer destinationOffset:dst_offset size:size];
}

static void mtl_cmd_copy_buffer_to_image(md_gpu_command_buffer_t cmd,
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

static void mtl_cmd_copy_image_region_to_buffer(md_gpu_command_buffer_t cmd,
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

static void mtl_cmd_copy_buffer_to_image_region(md_gpu_command_buffer_t cmd,
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

static void mtl_cmd_fill_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t buffer,
                            size_t offset,
                            size_t size,
                            uint8_t value) {
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit fillBuffer:buffer->buffer range:NSMakeRange(offset, size) value:value];
}

static void mtl_cmd_push_debug_group(md_gpu_command_buffer_t cmd, const char* label) {
    if (!cmd || !label) return;
    NSString* ns_label = [NSString stringWithUTF8String:label];
    if (cmd->compute_enc)
        [cmd->compute_enc pushDebugGroup:ns_label];
    else if (cmd->blit_enc)
        [cmd->blit_enc pushDebugGroup:ns_label];
    else
        [cmd->mtl_cmd pushDebugGroup:ns_label];
}

static void mtl_cmd_pop_debug_group(md_gpu_command_buffer_t cmd) {
    if (!cmd) return;
    if (cmd->compute_enc)
        [cmd->compute_enc popDebugGroup];
    else if (cmd->blit_enc)
        [cmd->blit_enc popDebugGroup];
    else
        [cmd->mtl_cmd popDebugGroup];
}

static id<MTLCommandBuffer> md_mtl_cmd_submit(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd) {
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

static struct md_gpu_event_record* md_mtl_event_record_find(struct md_gpu_device* dev, md_gpu_event_t event) {
    if (!dev || event.value == 0) return NULL;
    for (uint32_t i = 0; i < MD_GPU_EVENT_HISTORY_SIZE; ++i) {
        if (dev->event_records[i].value == event.value) return &dev->event_records[i];
    }
    return NULL;
}

md_gpu_cmd_t md_gpu_cmd_begin(md_gpu_device_t device, const char* label) {
    if (!device) return NULL;
    struct md_gpu_device* dev = device;
    struct md_gpu_queue* queue = &dev->compute_queue;

    struct md_gpu_cmd* c = (struct md_gpu_cmd*)calloc(1, sizeof(struct md_gpu_cmd));
    if (!c) return NULL;

    c->queue = queue;
    c->cmd = md_gpu_command_buffer_acquire((md_gpu_queue_t)queue);
    if (!c->cmd) {
        free(c);
        return NULL;
    }

    if (label) {
        mtl_cmd_push_debug_group(c->cmd, label);
        c->pushed_debug_group = true;
    }
    return (md_gpu_cmd_t)c;
}

md_gpu_event_t md_gpu_cmd_submit(md_gpu_cmd_t cmd_handle) {
    md_gpu_event_t event = {0};
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!c) return event;

    struct md_gpu_queue* queue = c->queue;
    struct md_gpu_device* dev = queue ? queue->dev : NULL;
    if (!queue || !dev || !c->cmd) {
        free(c);
        return event;
    }

    if (c->pushed_debug_group) {
        mtl_cmd_pop_debug_group(c->cmd);
    }

    struct md_gpu_event_record* record = &dev->event_records[dev->event_record_head++ % MD_GPU_EVENT_HISTORY_SIZE];
    if (record->mtl_cmd) {
        [record->mtl_cmd waitUntilCompleted];
        [record->mtl_cmd release];
    }
    record->value = 0;
    record->mtl_cmd = nil;

    event.value = ++dev->next_event_id;
    if (event.value == 0) event.value = ++dev->next_event_id;

    id<MTLCommandBuffer> submitted = md_mtl_cmd_submit(queue, c->cmd);
    if (submitted) {
        record->value = event.value;
        record->mtl_cmd = submitted;
    } else {
        event.value = 0;
    }

    free(c);
    return event;
}

bool md_gpu_event_is_complete(md_gpu_device_t device, md_gpu_event_t event) {
    if (!device || event.value == 0) return true;
    struct md_gpu_event_record* record = md_mtl_event_record_find(device, event);
    if (!record || !record->mtl_cmd) return true;
    MTLCommandBufferStatus s = record->mtl_cmd.status;
    return s == MTLCommandBufferStatusCompleted || s == MTLCommandBufferStatusError;
}

void md_gpu_event_wait(md_gpu_device_t device, md_gpu_event_t event) {
    if (!device || event.value == 0) return;
    struct md_gpu_event_record* record = md_mtl_event_record_find(device, event);
    if (!record || !record->mtl_cmd) return;

    [record->mtl_cmd waitUntilCompleted];
    [record->mtl_cmd release];
    record->mtl_cmd = nil;
    record->value = 0;
}

void md_gpu_cmd_event_wait(md_gpu_cmd_t cmd, md_gpu_event_t event) {
    /* Metal: single queue, submission order is preserved — no explicit fence needed. */
    (void)cmd;
    (void)event;
}

void md_gpu_cmd_use_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t buffer, gpu_usage_flags_t usage) {
    if (!cmd_handle || !buffer || usage == 0) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_merge_buffer_usage(c->cmd, buffer, usage);
}

void md_gpu_cmd_use_image(md_gpu_cmd_t cmd_handle, md_gpu_image_t image, gpu_usage_flags_t usage) {
    if (!cmd_handle || !image || usage == 0) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_merge_image_usage(c->cmd, image, usage);
}

void md_gpu_cmd_use_sampled_resource(md_gpu_cmd_t cmd_handle, md_gpu_image_t image, md_gpu_sampler_t sampler, gpu_usage_flags_t usage) {
    if (!cmd_handle || !image) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_merge_sampled_resource_usage(c->cmd, image, sampler, usage ? usage : GPU_USAGE_READ);
}

void md_gpu_cmd_bind_compute_pipeline(md_gpu_cmd_t cmd_handle, md_gpu_compute_pipeline_t pipeline) {
    if (!cmd_handle || !pipeline) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_bind_compute_pipeline(c->cmd, pipeline);
}

void md_gpu_cmd_dispatch(md_gpu_cmd_t cmd_handle, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z, const void* root_args, size_t root_args_size) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (root_args && root_args_size) {
        ASSERT(root_args_size <= MD_GPU_MAX_PUSH_CONSTANTS);
        MEMCPY(c->cmd->push_constants, root_args, root_args_size);
        c->cmd->push_constant_size = root_args_size;
    }
    mtl_cmd_dispatch(c->cmd, group_count_x, group_count_y, group_count_z);
}

void md_gpu_cmd_barrier(md_gpu_cmd_t cmd_handle, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_barrier(c->cmd, src_stage, dst_stage);
}

void md_gpu_cmd_push_debug_group(md_gpu_cmd_t cmd_handle, const char* label) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_push_debug_group(c->cmd, label);
}

void md_gpu_cmd_pop_debug_group(md_gpu_cmd_t cmd_handle) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_pop_debug_group(c->cmd);
}

void md_gpu_cmd_copy_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_copy_buffer(c->cmd, src, dst, size, src_offset, dst_offset);
}

void md_gpu_cmd_copy_image_region_to_buffer(md_gpu_cmd_t cmd_handle,
                                            md_gpu_image_t src_image,
                                            md_gpu_image_region_t src_region,
                                            md_gpu_buffer_t dst_buffer,
                                            size_t dst_offset) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_copy_image_region_to_buffer(c->cmd, src_image, src_region, dst_buffer, dst_offset);
}

void md_gpu_cmd_copy_buffer_to_image(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_copy_buffer_to_image(c->cmd, src_buffer, dst_image);
}

void md_gpu_cmd_copy_buffer_to_image_region(md_gpu_cmd_t cmd_handle,
                                            md_gpu_buffer_t src_buffer,
                                            size_t src_offset,
                                            md_gpu_image_t dst_image,
                                            md_gpu_image_region_t dst_region) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_copy_buffer_to_image_region(c->cmd, src_buffer, src_offset, dst_image, dst_region);
}

void md_gpu_cmd_fill_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t buffer, size_t offset, size_t size, uint8_t value) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    mtl_cmd_fill_buffer(c->cmd, buffer, offset, size, value);
}

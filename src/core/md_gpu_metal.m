// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <dispatch/dispatch.h>

#include <core/md_common.h> // Assert
#include <stdlib.h>         // malloc, free

/* Limits are defined in md_gpu.h */

#define MD_GPU_MAX_BIND_SLOTS 16
#define MD_GPU_MAX_PUSH_CONSTANTS 256
#define MD_GPU_MAX_RESOURCE_BINDINGS (MD_GPU_MAX_BIND_SLOTS * 3)

typedef struct md_gpu_queue* md_gpu_queue_t;
typedef struct md_gpu_command_buffer* md_gpu_command_buffer_t;

typedef struct md_gpu_cmd_buffer_usage_t {
    md_gpu_buffer_t buffer;
    gpu_usage_flags_t usage;
} md_gpu_cmd_buffer_usage_t;

typedef struct md_gpu_cmd_resource_binding_t {
    md_gpu_resource_kind_t kind;
    gpu_usage_flags_t usage;
    uint32_t set;
    uint32_t binding;
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
} md_gpu_cmd_resource_binding_t;

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
    id<MTLSharedEvent> event_timeline;

    struct md_gpu_queue compute_queue;

    uint64_t next_event_id;
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
    md_gpu_resource_binding_t resource_bindings[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint32_t resource_binding_count;
};

struct md_gpu_command_buffer {
    struct md_gpu_device* dev;

    /* Live Metal encoding state */
    id<MTLCommandBuffer>         mtl_cmd;
    id<MTLComputeCommandEncoder> compute_enc;
    id<MTLBlitCommandEncoder>    blit_enc;

    /* Live dispatch state */
    struct md_gpu_compute_pipeline* bound_pipeline;
    struct md_gpu_compute_pipeline* active_pipeline;
    id<MTLTexture> active_textures[MD_GPU_MAX_RESOURCE_BINDINGS];
    id<MTLSamplerState> active_samplers[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint8_t push_constants[MD_GPU_MAX_PUSH_CONSTANTS];
    size_t  push_constant_size;

    md_gpu_cmd_buffer_usage_t cmd_buffers[MD_GPU_MAX_BIND_SLOTS];
    uint32_t cmd_buffer_count;
    md_gpu_cmd_resource_binding_t cmd_resources[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint32_t cmd_resource_count;

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

static inline void mtl_cmd_reset_active_resource_slots(struct md_gpu_command_buffer* cmd) {
    MEMSET(cmd->active_textures, 0, sizeof(cmd->active_textures));
    MEMSET(cmd->active_samplers, 0, sizeof(cmd->active_samplers));
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

static bool mtl_cmd_merge_resource_binding(struct md_gpu_command_buffer* cmd, const md_gpu_resource_t* resource) {
    if (!cmd || !resource) return false;

    if (resource->kind == MD_GPU_RESOURCE_STORAGE_IMAGE ||
        resource->kind == MD_GPU_RESOURCE_SAMPLED_IMAGE) {
        if (!resource->image) return false;
    } else if (resource->kind == MD_GPU_RESOURCE_SAMPLER) {
        if (!resource->sampler) return false;
    } else {
        return false;
    }

    gpu_usage_flags_t usage = resource->usage;
    if (resource->kind != MD_GPU_RESOURCE_SAMPLER && usage == 0) {
        usage = GPU_USAGE_READ;
    }

    for (uint32_t i = 0; i < cmd->cmd_resource_count; ++i) {
        md_gpu_cmd_resource_binding_t* existing = &cmd->cmd_resources[i];
        if (existing->kind != resource->kind || existing->set != resource->set || existing->binding != resource->binding) {
            continue;
        }

        if ((resource->kind == MD_GPU_RESOURCE_STORAGE_IMAGE ||
             resource->kind == MD_GPU_RESOURCE_SAMPLED_IMAGE) && existing->image != resource->image) {
            return false;
        }

        if (resource->kind == MD_GPU_RESOURCE_SAMPLER && existing->sampler != resource->sampler) {
            return false;
        }

        existing->usage |= usage;
        return true;
    }

    if (cmd->cmd_resource_count >= MD_GPU_MAX_RESOURCE_BINDINGS) return false;

    md_gpu_cmd_resource_binding_t* binding = &cmd->cmd_resources[cmd->cmd_resource_count++];
    binding->kind = resource->kind;
    binding->usage = usage;
    binding->set = resource->set;
    binding->binding = resource->binding;
    binding->image = resource->image;
    binding->sampler = resource->sampler;
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
    cmd->cmd_resource_count = 0;
    cmd->bound_pipeline = NULL;
    cmd->active_pipeline = NULL;
    mtl_cmd_reset_active_resource_slots(cmd);
    cmd->push_constant_size = 0;
    cmd->next_free = NULL;
}

/* =============================
   Device / queue
   ============================= */

md_gpu_device_t md_gpu_device_create(void) {
    struct md_gpu_device* dev = (struct md_gpu_device*)calloc(1, sizeof(struct md_gpu_device));
    if (!dev) return NULL;

    dev->device = MTLCreateSystemDefaultDevice();
    if (!dev->device) {
        free(dev);
        return NULL;
    }

    dev->queue = [dev->device newCommandQueue];
    dev->event_timeline = [dev->device newSharedEvent];

    ASSERT(dev->queue);
    ASSERT(dev->event_timeline);
    if (!dev->queue || !dev->event_timeline) {
        [dev->event_timeline release];
        [dev->queue release];
        [dev->device release];
        free(dev);
        return NULL;
    }

    dev->compute_queue.dev = dev;
    /* cmd_free_list is zero-initialized by calloc */
    return dev;
}

void md_gpu_device_destroy(md_gpu_device_t device) {
    if (!device) return;

    if (device->next_event_id > 0 && [device->event_timeline signaledValue] < device->next_event_id) {
        md_gpu_event_t event = { .value = device->next_event_id };
        md_gpu_event_wait(device, event);
    }

    struct md_gpu_command_buffer* it = device->compute_queue.cmd_free_list;
    while (it) {
        struct md_gpu_command_buffer* next = it->next_free;
        free(it);
        it = next;
    }
    device->compute_queue.cmd_free_list = NULL;

    [device->event_timeline release];
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

    ASSERT(desc->resource_binding_count <= MD_GPU_MAX_RESOURCE_BINDINGS);
    if (desc->resource_bindings && desc->resource_binding_count > 0) {
        p->resource_binding_count = desc->resource_binding_count;
        MEMCPY(p->resource_bindings, desc->resource_bindings, sizeof(md_gpu_resource_binding_t) * desc->resource_binding_count);
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
    for (uint32_t i = 0; i < cmd->cmd_resource_count; ++i) {
        md_gpu_cmd_resource_binding_t usage = cmd->cmd_resources[i];
        if (!usage.image) continue;

        switch (usage.kind) {
            case MD_GPU_RESOURCE_STORAGE_IMAGE:
                [enc useResource:usage.image->texture usage:mtl_resource_usage(usage.usage)];
                break;
            case MD_GPU_RESOURCE_SAMPLED_IMAGE:
                [enc useResource:usage.image->texture usage:MTLResourceUsageRead];
                break;
            default:
                break;
        }
    }
}

static const md_gpu_resource_binding_t* mtl_pipeline_find_resource_binding(const struct md_gpu_compute_pipeline* pipeline, md_gpu_resource_kind_t kind, uint32_t set, uint32_t binding) {
    if (!pipeline) return NULL;

    for (uint32_t i = 0; i < pipeline->resource_binding_count; ++i) {
        const md_gpu_resource_binding_t* shader_binding = &pipeline->resource_bindings[i];
        if (shader_binding->kind == kind && shader_binding->set == set && shader_binding->binding == binding) {
            return shader_binding;
        }
    }

    return NULL;
}

static bool bind_cmd_resources(struct md_gpu_command_buffer* cmd, id<MTLComputeCommandEncoder> enc) {
    if (!cmd || !cmd->bound_pipeline || !enc) return false;

    id<MTLTexture> desired_textures[MD_GPU_MAX_RESOURCE_BINDINGS] = {0};
    id<MTLSamplerState> desired_samplers[MD_GPU_MAX_RESOURCE_BINDINGS] = {0};

    for (uint32_t i = 0; i < cmd->cmd_resource_count; ++i) {
        md_gpu_cmd_resource_binding_t binding = cmd->cmd_resources[i];
        const md_gpu_resource_binding_t* shader_binding = mtl_pipeline_find_resource_binding(cmd->bound_pipeline, binding.kind, binding.set, binding.binding);
        if (!shader_binding) return false;
        if (shader_binding->backend_binding >= MD_GPU_MAX_RESOURCE_BINDINGS) return false;

        switch (binding.kind) {
            case MD_GPU_RESOURCE_STORAGE_IMAGE:
            case MD_GPU_RESOURCE_SAMPLED_IMAGE:
                if (desired_textures[shader_binding->backend_binding] && desired_textures[shader_binding->backend_binding] != binding.image->texture) return false;
                desired_textures[shader_binding->backend_binding] = binding.image->texture;
                break;
            case MD_GPU_RESOURCE_SAMPLER:
                if (desired_samplers[shader_binding->backend_binding] && desired_samplers[shader_binding->backend_binding] != binding.sampler->sampler) return false;
                desired_samplers[shader_binding->backend_binding] = binding.sampler->sampler;
                break;
            default:
                return false;
        }
    }

    for (uint32_t i = 0; i < MD_GPU_MAX_RESOURCE_BINDINGS; ++i) {
        if (cmd->active_textures[i] != desired_textures[i]) {
            [enc setTexture:desired_textures[i] atIndex:i];
            cmd->active_textures[i] = desired_textures[i];
        }
        if (cmd->active_samplers[i] != desired_samplers[i]) {
            [enc setSamplerState:desired_samplers[i] atIndex:i];
            cmd->active_samplers[i] = desired_samplers[i];
        }
    }

    return true;
}

static inline void end_compute_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->compute_enc) {
        [cmd->compute_enc endEncoding];
        [cmd->compute_enc release];
        cmd->compute_enc = nil;
        cmd->active_pipeline = NULL;
        mtl_cmd_reset_active_resource_slots(cmd);
    }
}

static inline id<MTLBlitCommandEncoder> ensure_blit_enc(struct md_gpu_command_buffer* cmd) {
    if (cmd->compute_enc) {
        [cmd->compute_enc endEncoding];
        [cmd->compute_enc release];
        cmd->compute_enc = nil;
        cmd->active_pipeline = NULL;
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

    if (cmd->active_pipeline != cmd->bound_pipeline) {
        [enc setComputePipelineState:cmd->bound_pipeline->pso];
        cmd->active_pipeline = cmd->bound_pipeline;
    }
    ASSERT(bind_cmd_resources(cmd, enc));
    use_cmd_resources(cmd, enc);

    if (cmd->push_constant_size > 0)
        [enc setBytes:cmd->push_constants length:cmd->push_constant_size atIndex:0];

    MTLSize tg     = MTLSizeMake(cmd->bound_pipeline->tg_size[0], cmd->bound_pipeline->tg_size[1], cmd->bound_pipeline->tg_size[2]);
    MTLSize groups = MTLSizeMake(group_count_x, group_count_y, group_count_z);
    [enc dispatchThreadgroups:groups threadsPerThreadgroup:tg];

        /* Reset dispatch-local resource declarations so the next compute dispatch starts
            from a clean slate. Resource hazard tracking for useResource is re-issued per
            dispatch. */
    cmd->cmd_buffer_count = 0;
     cmd->cmd_resource_count = 0;
     cmd->push_constant_size = 0;
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

static bool md_mtl_cmd_submit(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd, uint64_t signal_value) {
    ASSERT(queue && cmd && cmd->mtl_cmd);
    ASSERT(signal_value != 0);

    end_compute_enc(cmd);
    end_blit_enc(cmd);

    [cmd->mtl_cmd encodeSignalEvent:queue->dev->event_timeline value:signal_value];

    [cmd->mtl_cmd commit];

    [cmd->mtl_cmd release];
    cmd->mtl_cmd = nil;

    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;

    return true;
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

    event.value = ++dev->next_event_id;
    if (event.value == 0) event.value = ++dev->next_event_id;

    if (!md_mtl_cmd_submit(queue, c->cmd, event.value)) {
        event.value = 0;
    }

    free(c);
    return event;
}

bool md_gpu_event_is_complete(md_gpu_device_t device, md_gpu_event_t event) {
    if (!device || event.value == 0) return true;
    return [device->event_timeline signaledValue] >= event.value;
}

void md_gpu_event_wait(md_gpu_device_t device, md_gpu_event_t event) {
    if (!device || event.value == 0) return;
    if ([device->event_timeline signaledValue] >= event.value) return;

    dispatch_semaphore_t semaphore = dispatch_semaphore_create(0);
    ASSERT(semaphore);
    if (!semaphore) return;

    MTLSharedEventListener* listener = [[MTLSharedEventListener alloc] initWithDispatchQueue:dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0)];
    ASSERT(listener);
    if (!listener) {
#if !OS_OBJECT_USE_OBJC
        dispatch_release(semaphore);
#endif
        return;
    }

    [device->event_timeline notifyListener:listener atValue:event.value block:^(id<MTLSharedEvent> shared_event, uint64_t value) {
        (void)shared_event;
        (void)value;
        dispatch_semaphore_signal(semaphore);
    }];

    dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
    [listener release];

#if !OS_OBJECT_USE_OBJC
    dispatch_release(semaphore);
#endif
}

void md_gpu_cmd_event_wait(md_gpu_cmd_t cmd_handle, md_gpu_event_t event) {
    if (!cmd_handle || event.value == 0) return;

    struct md_gpu_cmd* cmd = (struct md_gpu_cmd*)cmd_handle;
    if (!cmd->cmd || !cmd->cmd->mtl_cmd) return;

    end_compute_enc(cmd->cmd);
    end_blit_enc(cmd->cmd);
    [cmd->cmd->mtl_cmd encodeWaitForEvent:cmd->cmd->dev->event_timeline value:event.value];
}

static bool mtl_cmd_declare_resources(struct md_gpu_command_buffer* cmd, const md_gpu_resource_t* resources, uint32_t resource_count) {
    if (!cmd) return false;
    cmd->cmd_buffer_count = 0;
    cmd->cmd_resource_count = 0;

    if (!resources || resource_count == 0) return true;

    for (uint32_t i = 0; i < resource_count; ++i) {
        const md_gpu_resource_t* resource = &resources[i];
        switch (resource->kind) {
            case MD_GPU_RESOURCE_BUFFER:
                if (!mtl_cmd_merge_buffer_usage(cmd, resource->buffer, resource->usage)) return false;
                break;
            case MD_GPU_RESOURCE_STORAGE_IMAGE:
            case MD_GPU_RESOURCE_SAMPLED_IMAGE:
            case MD_GPU_RESOURCE_SAMPLER:
                if (!mtl_cmd_merge_resource_binding(cmd, resource)) return false;
                break;
            default:
                return false;
        }
    }

    return true;
}

bool md_gpu_cmd_dispatch(md_gpu_cmd_t cmd_handle, const md_gpu_compute_dispatch_t* dispatch) {
    if (!cmd_handle || !dispatch || !dispatch->pipeline) return false;
    if (dispatch->resource_count > 0 && !dispatch->resources) return false;
    if (dispatch->group_count[0] == 0 || dispatch->group_count[1] == 0 || dispatch->group_count[2] == 0) return false;
    if (dispatch->root_args_size > 0 && !dispatch->root_args) return false;
    if (dispatch->root_args_size > MD_GPU_MAX_PUSH_CONSTANTS) {
        ASSERT(false && "Root argument payload exceeds push constant limit");
        return false;
    }

    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!mtl_cmd_declare_resources(c->cmd, dispatch->resources, dispatch->resource_count)) return false;

    mtl_cmd_bind_compute_pipeline(c->cmd, dispatch->pipeline);
    c->cmd->push_constant_size = 0;
    if (dispatch->root_args && dispatch->root_args_size) {
        MEMCPY(c->cmd->push_constants, dispatch->root_args, dispatch->root_args_size);
        c->cmd->push_constant_size = dispatch->root_args_size;
    }

    mtl_cmd_dispatch(c->cmd, dispatch->group_count[0], dispatch->group_count[1], dispatch->group_count[2]);
    return true;
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

// md_gpu_metal.mm
// Metal backend implementation for md_gpu.h

#include "md_gpu.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <dispatch/dispatch.h>

#include <core/md_common.h> // Assert
#include <core/md_os.h>
#include <stdlib.h>         // malloc, free
#include <string.h>         // strcmp

/* Limits are defined in md_gpu.h */

#define MD_GPU_MAX_BIND_SLOTS 16
#define MD_GPU_MAX_PUSH_CONSTANTS 256
#define MD_GPU_MAX_RESOURCE_BINDINGS (MD_GPU_MAX_BIND_SLOTS * 3)
#define MD_GPU_MTL_TRANSIENT_PAGE_SIZE (4u * 1024u * 1024u)
typedef struct md_gpu_command_buffer* md_gpu_command_buffer_t;

typedef struct md_gpu_cmd_buffer_usage_t {
    md_gpu_buffer_t buffer;
    md_gpu_usage_flags_t usage;
} md_gpu_cmd_buffer_usage_t;

typedef struct md_gpu_cmd_resource_binding_t {
    md_gpu_resource_kind_t kind;
    md_gpu_usage_flags_t usage;
    uint32_t set;
    uint32_t binding;
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
} md_gpu_cmd_resource_binding_t;

typedef struct md_mtl_transient_page_t md_mtl_transient_page_t;

typedef enum md_mtl_transient_page_list_state_t {
    MD_MTL_TRANSIENT_PAGE_LIST_NONE = 0,
    MD_MTL_TRANSIENT_PAGE_LIST_AVAILABLE,
    MD_MTL_TRANSIENT_PAGE_LIST_RETIRED,
} md_mtl_transient_page_list_state_t;

typedef struct md_mtl_transient_page_pool_t {
    md_mutex_t mutex;
    md_mtl_transient_page_t* available_pages;
    md_mtl_transient_page_t* retired_pages;
} md_mtl_transient_page_pool_t;

struct md_mtl_transient_page_t {
    md_gpu_buffer_t buffer;
    void* cpu_ptr;
    size_t capacity;
    size_t offset;
    uint64_t retire_value;
    bool dedicated;
    md_mtl_transient_page_list_state_t list_state;
    md_mtl_transient_page_t* next;
    md_mtl_transient_page_t* cmd_next;
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
    id<MTLSharedEvent> event_timeline;

    struct md_gpu_queue compute_queue;
    md_mtl_transient_page_pool_t transient_page_pool;

    uint64_t next_event_id;
};

struct md_gpu_buffer {
    struct md_gpu_device* device;
    id<MTLBuffer> buffer;
    void* cpu_ptr;  /* non-NULL for CPU_VISIBLE buffers; valid until destroy */
    size_t size;
    md_gpu_buffer_flags_t flags;
};

struct md_gpu_image {
    struct md_gpu_device* device;
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
    md_mtl_transient_page_t* current_transient_page;
    md_mtl_transient_page_t* transient_pages;

    struct md_gpu_command_buffer* next_free;
};

struct md_gpu_readback {
    struct md_gpu_device* dev;
    md_mtl_transient_page_t* page;
    md_gpu_event_t event;
    size_t size;
};

struct md_gpu_cmd {
    struct md_gpu_queue* queue;
    struct md_gpu_command_buffer* cmd;
    enum {
        MD_MTL_CMD_STATE_ACQUIRED = 0,
        MD_MTL_CMD_STATE_RECORDING,
        MD_MTL_CMD_STATE_ENDED,
    } state;
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

static bool mtl_transient_pool_init(struct md_gpu_device* dev) {
    if (!dev) return false;
    MEMSET(&dev->transient_page_pool, 0, sizeof(dev->transient_page_pool));
    return md_mutex_init(&dev->transient_page_pool.mutex);
}

static void mtl_transient_page_destroy(md_mtl_transient_page_t* page) {
    if (!page) return;
    if (page->buffer) {
        md_gpu_buffer_destroy(page->buffer);
        page->buffer = NULL;
    }
    free(page);
}

static void mtl_transient_page_release(struct md_gpu_device* dev, md_mtl_transient_page_t* page) {
    if (!dev || !page) return;

    page->offset = 0;
    page->retire_value = 0;
    page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_NONE;
    page->next = NULL;
    page->cmd_next = NULL;

    if (page->dedicated) {
        mtl_transient_page_destroy(page);
        return;
    }

    md_mtl_transient_page_pool_t* pool = &dev->transient_page_pool;
    md_mutex_lock(&pool->mutex);
    page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_AVAILABLE;
    page->next = pool->available_pages;
    pool->available_pages = page;
    md_mutex_unlock(&pool->mutex);
}

static void mtl_transient_page_enqueue_retired(struct md_gpu_device* dev, md_mtl_transient_page_t* page) {
    if (!dev || !page) return;

    md_mtl_transient_page_pool_t* pool = &dev->transient_page_pool;
    md_mutex_lock(&pool->mutex);
    if (page->list_state != MD_MTL_TRANSIENT_PAGE_LIST_RETIRED) {
        page->next = pool->retired_pages;
        pool->retired_pages = page;
        page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_RETIRED;
    }
    md_mutex_unlock(&pool->mutex);
}

static void mtl_transient_pool_retire(struct md_gpu_device* dev) {
    if (!dev || !dev->event_timeline) return;

    const uint64_t completed = [dev->event_timeline signaledValue];
    md_mtl_transient_page_pool_t* pool = &dev->transient_page_pool;
    md_mutex_lock(&pool->mutex);

    md_mtl_transient_page_t** pp = &pool->retired_pages;
    while (*pp) {
        md_mtl_transient_page_t* page = *pp;
        if (page->retire_value == 0 || page->retire_value > completed) {
            pp = &page->next;
            continue;
        }

        *pp = page->next;
        page->next = NULL;
        page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_NONE;
        md_mutex_unlock(&pool->mutex);
        mtl_transient_page_release(dev, page);
        md_mutex_lock(&pool->mutex);
        pp = &pool->retired_pages;
    }

    md_mutex_unlock(&pool->mutex);
}

static void mtl_transient_pool_free(struct md_gpu_device* dev) {
    if (!dev) return;

    md_mtl_transient_page_pool_t* pool = &dev->transient_page_pool;
    md_mutex_lock(&pool->mutex);
    md_mtl_transient_page_t* available = pool->available_pages;
    md_mtl_transient_page_t* retired = pool->retired_pages;
    pool->available_pages = NULL;
    pool->retired_pages = NULL;
    md_mutex_unlock(&pool->mutex);

    while (available) {
        md_mtl_transient_page_t* next = available->next;
        mtl_transient_page_destroy(available);
        available = next;
    }
    while (retired) {
        md_mtl_transient_page_t* next = retired->next;
        mtl_transient_page_destroy(retired);
        retired = next;
    }

    md_mutex_destroy(&pool->mutex);
    MEMSET(pool, 0, sizeof(*pool));
}

static md_mtl_transient_page_t* mtl_transient_page_acquire(struct md_gpu_device* dev, size_t min_capacity, bool dedicated) {
    if (!dev || min_capacity == 0) return NULL;

    const size_t page_capacity = dedicated ? min_capacity : MAX((size_t)MD_GPU_MTL_TRANSIENT_PAGE_SIZE, min_capacity);
    md_mtl_transient_page_pool_t* pool = &dev->transient_page_pool;

    if (!dedicated) {
        md_mutex_lock(&pool->mutex);
        md_mtl_transient_page_t** pp = &pool->available_pages;
        while (*pp) {
            md_mtl_transient_page_t* page = *pp;
            if (page->capacity >= min_capacity) {
                *pp = page->next;
                page->next = NULL;
                page->cmd_next = NULL;
                page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_NONE;
                md_mutex_unlock(&pool->mutex);
                return page;
            }
            pp = &page->next;
        }
        md_mutex_unlock(&pool->mutex);
    }

    md_mtl_transient_page_t* page = (md_mtl_transient_page_t*)calloc(1, sizeof(*page));
    if (!page) return NULL;

    page->buffer = md_gpu_buffer_create((md_gpu_device_t)dev, &(md_gpu_buffer_desc_t){
        .size = page_capacity,
        .flags = MD_GPU_BUFFER_CPU_VISIBLE,
    });
    if (!page->buffer) {
        free(page);
        return NULL;
    }

    page->cpu_ptr = md_gpu_buffer_cpu_ptr(page->buffer);
    page->capacity = page_capacity;
    page->dedicated = dedicated;
    page->list_state = MD_MTL_TRANSIENT_PAGE_LIST_NONE;
    return page;
}

static bool mtl_cmd_attach_transient_page(struct md_gpu_command_buffer* cmd, md_mtl_transient_page_t* page) {
    if (!cmd || !page) return false;

    for (md_mtl_transient_page_t* it = cmd->transient_pages; it; it = it->cmd_next) {
        if (it == page) return true;
    }

    page->cmd_next = cmd->transient_pages;
    cmd->transient_pages = page;
    return true;
}

static void mtl_cmd_release_transient_pages(struct md_gpu_command_buffer* cmd, bool submitted) {
    if (!cmd || !cmd->dev) return;

    md_mtl_transient_page_t* page = cmd->transient_pages;
    while (page) {
        md_mtl_transient_page_t* next = page->cmd_next;
        page->cmd_next = NULL;

        if (!submitted) {
            mtl_transient_page_release(cmd->dev, page);
        }

        page = next;
    }

    cmd->current_transient_page = NULL;
    cmd->transient_pages = NULL;
}

static inline void mtl_cmd_reset_active_resource_slots(struct md_gpu_command_buffer* cmd) {
    MEMSET(cmd->active_textures, 0, sizeof(cmd->active_textures));
    MEMSET(cmd->active_samplers, 0, sizeof(cmd->active_samplers));
}

static bool mtl_cmd_merge_buffer_usage(struct md_gpu_command_buffer* cmd, md_gpu_buffer_t buffer, md_gpu_usage_flags_t usage) {
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

    md_gpu_usage_flags_t usage = resource->usage;
    if (resource->kind != MD_GPU_RESOURCE_SAMPLER && usage == 0) {
        usage = MD_GPU_USAGE_READ;
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
    cmd->current_transient_page = NULL;
    cmd->transient_pages = NULL;
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
    if (!mtl_transient_pool_init(dev)) {
        [dev->event_timeline release];
        [dev->queue release];
        [dev->device release];
        free(dev);
        return NULL;
    }
    /* cmd_free_list is zero-initialized by calloc */
    return dev;
}

void md_gpu_device_destroy(md_gpu_device_t device) {
    if (!device) return;

    if (device->next_event_id > 0 && [device->event_timeline signaledValue] < device->next_event_id) {
        md_gpu_event_t event = { .queue = (md_gpu_queue_t)&device->compute_queue, .value = device->next_event_id };
        md_gpu_event_wait(event);
    }

    struct md_gpu_command_buffer* it = device->compute_queue.cmd_free_list;
    while (it) {
        struct md_gpu_command_buffer* next = it->next_free;
        free(it);
        it = next;
    }
    device->compute_queue.cmd_free_list = NULL;

    mtl_transient_pool_retire(device);
    mtl_transient_pool_free(device);

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

md_gpu_queue_t md_gpu_queue_graphics(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&device->compute_queue;
}

md_gpu_queue_t md_gpu_queue_compute(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&device->compute_queue;
}

md_gpu_queue_t md_gpu_queue_transfer(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&device->compute_queue;
}

/* =============================
   Buffers
   ============================= */

md_gpu_buffer_t md_gpu_buffer_create(md_gpu_device_t device,
                                     const md_gpu_buffer_desc_t* desc) {
    struct md_gpu_buffer* buf = (struct md_gpu_buffer*)calloc(1, sizeof(struct md_gpu_buffer));
    buf->device = device;
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
    img->device = device;
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
    if (!device || !desc || !desc->shader_bytes || desc->shader_byte_size == 0) return NULL;

    NSError* err = nil;
    const char* entry_point = desc->entry_point ? desc->entry_point : "main";

    dispatch_data_t data = dispatch_data_create(
        desc->shader_bytes,
        desc->shader_byte_size,
        dispatch_get_main_queue(),
        DISPATCH_DATA_DESTRUCTOR_DEFAULT);

    id<MTLLibrary> lib =
        [device->device newLibraryWithData:data error:&err];
    dispatch_release(data);
    ASSERT(lib && !err);
    if (!lib || err) return NULL;

    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:entry_point]];
    if (!fn && strcmp(entry_point, "main") == 0) {
        // Slang renames 'main' to 'main_0' for Metal (Metal reserves 'main').
        // spirv-cross used 'main0'; try both for compatibility.
        fn = [lib newFunctionWithName:@"main_0"];
        if (!fn) fn = [lib newFunctionWithName:@"main0"];
    }
    ASSERT(fn);
    if (!fn) {
        [lib release];
        return NULL;
    }

    struct md_gpu_compute_pipeline* p =
        (struct md_gpu_compute_pipeline*)calloc(1, sizeof(struct md_gpu_compute_pipeline));

    MTLComputePipelineDescriptor* pipeline_desc = [[MTLComputePipelineDescriptor alloc] init];
    pipeline_desc.computeFunction = fn;
    if (desc->label) {
        pipeline_desc.label = [NSString stringWithUTF8String:desc->label];
    }

    p->pso = [device->device newComputePipelineStateWithDescriptor:pipeline_desc options:0 reflection:nil error:&err];
    [pipeline_desc release];
    ASSERT(p->pso && !err);
    if (!p->pso || err) {
        [fn release];
        [lib release];
        free(p);
        return NULL;
    }

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

static inline MTLResourceUsage mtl_resource_usage(md_gpu_usage_flags_t usage) {
    MTLResourceUsage mtl_usage = 0;
    if (usage & MD_GPU_USAGE_READ) {
        mtl_usage |= MTLResourceUsageRead;
    }
    if (usage & MD_GPU_USAGE_WRITE) {
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
        if (!cmd) return NULL;
        cmd->dev = queue->dev;
    }

    cmd_reset(cmd);
    cmd->mtl_cmd = [[queue->dev->queue commandBuffer] retain];
    return cmd;
}

static inline bool md_mtl_cmd_is_recording(const struct md_gpu_cmd* cmd) {
    return cmd && cmd->cmd && cmd->state == MD_MTL_CMD_STATE_RECORDING;
}

static inline void md_mtl_cmd_discard(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd) {
    if (!queue || !cmd) return;
    mtl_cmd_release_transient_pages(cmd, false);
    cmd_reset(cmd);
    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;
}

static md_gpu_transient_t mtl_cmd_temp_alloc(struct md_gpu_command_buffer* cmd, size_t size, size_t alignment) {
    md_gpu_transient_t result = {0};
    if (!cmd || !cmd->dev || size == 0) return result;
    if (alignment == 0) alignment = 256;

    md_mtl_transient_page_t* page = cmd->current_transient_page;
    const bool dedicated = (size > ((size_t)MD_GPU_MTL_TRANSIENT_PAGE_SIZE / 2));
    size_t aligned_offset = SIZE_MAX;

    if (!dedicated && page) {
        aligned_offset = (page->offset + alignment - 1) & ~(alignment - 1);
        if (aligned_offset + size > page->capacity) {
            cmd->current_transient_page = NULL;
            page = NULL;
        }
    }

    if (!page) {
        page = mtl_transient_page_acquire(cmd->dev, size, dedicated);
        if (!page) return result;
        if (!dedicated) {
            cmd->current_transient_page = page;
        }
        aligned_offset = 0;
    }

    if (aligned_offset == SIZE_MAX) {
        aligned_offset = (page->offset + alignment - 1) & ~(alignment - 1);
    }
    if (aligned_offset + size > page->capacity) return result;
    if (!mtl_cmd_attach_transient_page(cmd, page)) return result;

    page->offset = aligned_offset + size;
    result.buffer = page->buffer;
    result.cpu_ptr = page->cpu_ptr ? (uint8_t*)page->cpu_ptr + aligned_offset : NULL;
    result.offset = aligned_offset;
    result.size = size;
    return result;
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

static bool mtl_cmd_copy_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t src,
                            md_gpu_buffer_t dst,
                            size_t size,
                            size_t src_offset,
                            size_t dst_offset) {
    if (!cmd || !src || !dst || size == 0) return false;
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit copyFromBuffer:src->buffer sourceOffset:src_offset toBuffer:dst->buffer destinationOffset:dst_offset size:size];
    return true;
}

static bool mtl_image_copy_layout(MTLOrigin* origin, MTLSize* size, size_t* buffer_offset, size_t* bytes_per_row, size_t* bytes_per_image, const md_gpu_image_desc_t* desc, const md_gpu_buffer_image_copy_t* copy) {
    if (!origin || !size || !buffer_offset || !bytes_per_row || !bytes_per_image || !desc) return false;

    const uint32_t bpp = md_gpu_format_bytes_per_pixel(desc->format);
    if (bpp == 0) return false;

    const md_gpu_image_region_t image_region = copy ? copy->image_region : (md_gpu_image_region_t){0};
    const uint32_t ox = image_region.offset[0];
    const uint32_t oy = image_region.offset[1];
    const uint32_t oz = image_region.offset[2];
    if (ox >= desc->width || oy >= desc->height || oz >= desc->depth) return false;

    const uint32_t ew = image_region.extent[0] ? image_region.extent[0] : desc->width  - ox;
    const uint32_t eh = image_region.extent[1] ? image_region.extent[1] : desc->height - oy;
    const uint32_t ed = image_region.extent[2] ? image_region.extent[2] : desc->depth  - oz;
    if (ew == 0 || eh == 0 || ed == 0) return false;
    if (ew > desc->width - ox || eh > desc->height - oy || ed > desc->depth - oz) return false;

    const size_t tight_row = (size_t)ew * (size_t)bpp;
    const size_t row_pitch = (copy && copy->bytes_per_row) ? copy->bytes_per_row : tight_row;
    if (row_pitch < tight_row) return false;

    const size_t tight_image = row_pitch * (size_t)eh;
    const size_t image_pitch = (copy && copy->bytes_per_image) ? copy->bytes_per_image : tight_image;
    if (image_pitch < tight_image) return false;

    *origin = MTLOriginMake(ox, oy, oz);
    *size = MTLSizeMake(ew, eh, ed);
    *buffer_offset = copy ? copy->buffer_offset : 0;
    *bytes_per_row = row_pitch;
    *bytes_per_image = image_pitch;
    return true;
}

static bool mtl_cmd_copy_buffer_to_image_layout(md_gpu_command_buffer_t cmd,
                                            md_gpu_buffer_t src_buffer,
                                            md_gpu_image_t dst_image,
                                            const md_gpu_buffer_image_copy_t* copy) {
    if (!cmd || !src_buffer || !dst_image) return false;

    MTLOrigin origin;
    MTLSize sz;
    size_t src_offset = 0;
    size_t bytes_per_row = 0;
    size_t bytes_per_image = 0;
    if (!mtl_image_copy_layout(&origin, &sz, &src_offset, &bytes_per_row, &bytes_per_image, &dst_image->desc, copy)) return false;

    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit copyFromBuffer:src_buffer->buffer
            sourceOffset:src_offset
       sourceBytesPerRow:bytes_per_row
     sourceBytesPerImage:bytes_per_image
              sourceSize:sz
               toTexture:dst_image->texture
        destinationSlice:0
        destinationLevel:0
       destinationOrigin:origin];
        return true;
}

static bool mtl_cmd_copy_image_to_buffer(md_gpu_command_buffer_t cmd,
                                                                                        md_gpu_image_t src_image,
                                                                                        md_gpu_buffer_t dst_buffer,
                                                                                        const md_gpu_buffer_image_copy_t* copy) {
        if (!cmd || !src_image || !dst_buffer) return false;

        MTLOrigin origin;
        MTLSize sz;
        size_t dst_offset = 0;
        size_t bytes_per_row = 0;
        size_t bytes_per_image = 0;
        if (!mtl_image_copy_layout(&origin, &sz, &dst_offset, &bytes_per_row, &bytes_per_image, &src_image->desc, copy)) return false;

        id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
        [blit copyFromTexture:src_image->texture
                             sourceSlice:0
                             sourceLevel:0
                            sourceOrigin:origin
                                sourceSize:sz
                                    toBuffer:dst_buffer->buffer
                 destinationOffset:dst_offset
        destinationBytesPerRow:bytes_per_row
    destinationBytesPerImage:bytes_per_image];
        return true;
}

static bool mtl_cmd_fill_buffer(md_gpu_command_buffer_t cmd,
                            md_gpu_buffer_t buffer,
                            size_t offset,
                            size_t size,
                            uint8_t value) {
        if (!cmd || !buffer || size == 0) return false;
    id<MTLBlitCommandEncoder> blit = ensure_blit_enc(cmd);
    [blit fillBuffer:buffer->buffer range:NSMakeRange(offset, size) value:value];
        return true;
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

static void md_mtl_cmd_encode_waits(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd, const md_gpu_event_t* waits, size_t wait_count) {
    if (!queue || !cmd || !cmd->mtl_cmd || !waits) return;

    for (size_t i = 0; i < wait_count; ++i) {
        const md_gpu_event_t event = waits[i];
        if (!md_gpu_event_is_valid(event)) continue;

        struct md_gpu_queue* wait_queue = (struct md_gpu_queue*)event.queue;
        if (!wait_queue) wait_queue = queue;
        if (!wait_queue->dev || !wait_queue->dev->event_timeline) continue;

        [cmd->mtl_cmd encodeWaitForEvent:wait_queue->dev->event_timeline value:event.value];
    }
}

static bool md_mtl_cmd_commit(struct md_gpu_queue* queue, struct md_gpu_command_buffer* cmd, uint64_t signal_value) {
    ASSERT(queue && cmd && cmd->mtl_cmd);

    if (signal_value) {
        [cmd->mtl_cmd encodeSignalEvent:queue->dev->event_timeline value:signal_value];
    }

    [cmd->mtl_cmd commit];

    [cmd->mtl_cmd release];
    cmd->mtl_cmd = nil;
    cmd->current_transient_page = NULL;
    cmd->transient_pages = NULL;

    cmd->next_free = queue->cmd_free_list;
    queue->cmd_free_list = cmd;

    return true;
}

static bool md_mtl_cmd_begin_acquired(md_gpu_cmd_t cmd_handle, const char* label) {
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!c || !c->queue || !c->cmd || c->state != MD_MTL_CMD_STATE_ACQUIRED) return false;

    if (label) {
        mtl_cmd_push_debug_group(c->cmd, label);
        c->pushed_debug_group = true;
    }

    c->state = MD_MTL_CMD_STATE_RECORDING;
    return true;
}

md_gpu_cmd_t md_gpu_cmd_begin(md_gpu_queue_t queue_handle, const char* label) {
    if (!queue_handle) return NULL;

    struct md_gpu_cmd* c = (struct md_gpu_cmd*)calloc(1, sizeof(struct md_gpu_cmd));
    if (!c) return NULL;

    c->queue = (struct md_gpu_queue*)queue_handle;
    c->cmd = md_gpu_command_buffer_acquire(queue_handle);
    if (!c->cmd) {
        free(c);
        return NULL;
    }
    c->state = MD_MTL_CMD_STATE_ACQUIRED;

    if (!md_mtl_cmd_begin_acquired((md_gpu_cmd_t)c, label)) {
        md_gpu_cmd_discard((md_gpu_cmd_t)c);
        return NULL;
    }

    return (md_gpu_cmd_t)c;
}

bool md_gpu_cmd_end(md_gpu_cmd_t cmd_handle) {
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;

    end_compute_enc(c->cmd);
    end_blit_enc(c->cmd);

    if (c->pushed_debug_group) {
        [c->cmd->mtl_cmd popDebugGroup];
        c->pushed_debug_group = false;
    }

    c->state = MD_MTL_CMD_STATE_ENDED;
    return true;
}

void md_gpu_cmd_discard(md_gpu_cmd_t cmd_handle) {
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!c) return;

    if (c->queue && c->cmd) {
        md_mtl_cmd_discard(c->queue, c->cmd);
        c->cmd = NULL;
    }

    free(c);
}

md_gpu_transient_t md_gpu_cmd_temp_alloc(md_gpu_cmd_t cmd_handle, size_t size) {
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return (md_gpu_transient_t){0};
    return mtl_cmd_temp_alloc(c->cmd, size, 256);
}

md_gpu_event_t md_gpu_queue_submit(md_gpu_queue_t queue_handle, const md_gpu_queue_submit_desc_t* desc) {
    md_gpu_event_t event = {0};
    if (!queue_handle || !desc || !desc->cmds || desc->cmd_count == 0) return event;
    if (desc->wait_count > 0 && !desc->waits) return event;

    struct md_gpu_queue* queue = (struct md_gpu_queue*)queue_handle;
    struct md_gpu_device* dev = queue ? queue->dev : NULL;
    if (!queue || !dev) {
        return event;
    }

    const md_gpu_cmd_t* cmds = desc->cmds;
    const size_t cmd_count = desc->cmd_count;

    for (size_t i = 0; i < cmd_count; ++i) {
        struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmds[i];
        if (!c || c->queue != queue || !c->cmd || c->state != MD_MTL_CMD_STATE_ENDED) {
            return event;
        }
    }

    event.queue = (md_gpu_queue_t)queue;
    event.value = ++dev->next_event_id;
    if (event.value == 0) event.value = ++dev->next_event_id;

    for (size_t i = 0; i < cmd_count; ++i) {
        struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmds[i];
        md_mtl_transient_page_t* page = c->cmd->transient_pages;
        while (page) {
            md_mtl_transient_page_t* next = page->cmd_next;
            page->cmd_next = NULL;
            page->retire_value = event.value;
            mtl_transient_page_enqueue_retired(dev, page);
            page = next;
        }
        c->cmd->current_transient_page = NULL;
        c->cmd->transient_pages = NULL;
        if (i == 0) md_mtl_cmd_encode_waits(queue, c->cmd, desc->waits, desc->wait_count);
        const uint64_t signal_value = (i + 1 == cmd_count) ? event.value : 0;
        if (!md_mtl_cmd_commit(queue, c->cmd, signal_value)) {
            event.queue = NULL;
            event.value = 0;
            break;
        }
        c->cmd = NULL;
        free(c);
    }

    return event;
}

bool md_gpu_event_is_complete(md_gpu_event_t event) {
    if (event.value == 0) return true;
    struct md_gpu_queue* queue = (struct md_gpu_queue*)event.queue;
    md_gpu_device_t event_device = queue ? queue->dev : NULL;
    if (!event_device) return true;
    mtl_transient_pool_retire(event_device);
    return [event_device->event_timeline signaledValue] >= event.value;
}

void md_gpu_event_wait(md_gpu_event_t event) {
    if (event.value == 0) return;
    struct md_gpu_queue* queue = (struct md_gpu_queue*)event.queue;
    md_gpu_device_t event_device = queue ? queue->dev : NULL;
    if (!event_device) return;
    if ([event_device->event_timeline signaledValue] >= event.value) {
        mtl_transient_pool_retire(event_device);
        return;
    }

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

    [event_device->event_timeline notifyListener:listener atValue:event.value block:^(id<MTLSharedEvent> shared_event, uint64_t value) {
        (void)shared_event;
        (void)value;
        dispatch_semaphore_signal(semaphore);
    }];

    dispatch_semaphore_wait(semaphore, DISPATCH_TIME_FOREVER);
    [listener release];
    mtl_transient_pool_retire(event_device);

#if !OS_OBJECT_USE_OBJC
    dispatch_release(semaphore);
#endif
}

static bool mtl_cmd_declare_resources(struct md_gpu_command_buffer* cmd, const md_gpu_resource_t* resources, uint32_t resource_count) {
    if (!cmd) return false;
    cmd->cmd_buffer_count = 0;
    cmd->cmd_resource_count = 0;

    if (!resources || resource_count == 0) return true;

    for (uint32_t i = 0; i < resource_count; ++i) {
        const md_gpu_resource_t* resource = &resources[i];
        switch (resource->kind) {
            case MD_GPU_RESOURCE_BUFFER_USAGE:
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
    if (!md_mtl_cmd_is_recording(c)) return false;
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

bool md_gpu_cmd_barrier(md_gpu_cmd_t cmd_handle, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!cmd_handle) return false;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;
    mtl_cmd_barrier(c->cmd, src_stage, dst_stage);
    return true;
}

void md_gpu_cmd_push_debug_group(md_gpu_cmd_t cmd_handle, const char* label) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return;
    mtl_cmd_push_debug_group(c->cmd, label);
}

void md_gpu_cmd_pop_debug_group(md_gpu_cmd_t cmd_handle) {
    if (!cmd_handle) return;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return;
    mtl_cmd_pop_debug_group(c->cmd);
}

bool md_gpu_cmd_copy_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!cmd_handle) return false;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;
    return mtl_cmd_copy_buffer(c->cmd, src, dst, size, src_offset, dst_offset);
}

bool md_gpu_cmd_copy_image_to_buffer(md_gpu_cmd_t cmd_handle, md_gpu_image_t src_image, md_gpu_buffer_t dst_buffer, const md_gpu_buffer_image_copy_t* copy) {
    if (!cmd_handle) return false;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;
    return mtl_cmd_copy_image_to_buffer(c->cmd, src_image, dst_buffer, copy);
}

bool md_gpu_cmd_copy_buffer_to_image_layout(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image, const md_gpu_buffer_image_copy_t* copy) {
    if (!cmd_handle) return false;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;
    return mtl_cmd_copy_buffer_to_image_layout(c->cmd, src_buffer, dst_image, copy);
}

bool md_gpu_cmd_fill_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t buffer, size_t offset, size_t size, uint8_t value) {
    if (!cmd_handle) return false;
    struct md_gpu_cmd* c = (struct md_gpu_cmd*)cmd_handle;
    if (!md_mtl_cmd_is_recording(c)) return false;
    return mtl_cmd_fill_buffer(c->cmd, buffer, offset, size, value);
}

md_gpu_readback_t md_gpu_readback_buffer(md_gpu_buffer_t src_buffer, size_t src_offset, size_t size, md_gpu_event_t after) {
    if (!src_buffer || size == 0) return NULL;

    struct md_gpu_buffer* src = (struct md_gpu_buffer*)src_buffer;
    struct md_gpu_device* dev = src->device;
    if (!dev) return NULL;
    if (src_offset > src->size || size > src->size - src_offset) return NULL;

    md_mtl_transient_page_t* page = mtl_transient_page_acquire(dev, size, true);
    if (!page) return NULL;

    struct md_gpu_readback* readback = (struct md_gpu_readback*)calloc(1, sizeof(*readback));
    if (!readback) {
        mtl_transient_page_release(dev, page);
        return NULL;
    }

    md_gpu_queue_t queue = md_gpu_queue_transfer((md_gpu_device_t)dev);
    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "readback_buffer");
    if (!cmd) goto fail;
    if (!md_gpu_cmd_copy_buffer(cmd, src_buffer, page->buffer, size, src_offset, 0)) goto fail;
    if (!md_gpu_cmd_end(cmd)) goto fail;

    md_gpu_event_t event = md_gpu_event_is_valid(after)
        ? md_gpu_queue_submit_one_after(queue, cmd, after, MD_GPU_BARRIER_STAGE_TRANSFER)
        : md_gpu_queue_submit_one(queue, cmd);
    if (!md_gpu_event_is_valid(event)) goto fail;

    readback->dev = dev;
    readback->page = page;
    readback->event = event;
    readback->size = size;
    return (md_gpu_readback_t)readback;

fail:
    if (cmd) md_gpu_cmd_discard(cmd);
    free(readback);
    mtl_transient_page_release(dev, page);
    return NULL;
}

md_gpu_readback_t md_gpu_readback_image(md_gpu_image_t src_image, md_gpu_image_region_t src_region, md_gpu_event_t after) {
    if (!src_image) return NULL;

    struct md_gpu_image* image = (struct md_gpu_image*)src_image;
    struct md_gpu_device* dev = image->device;
    if (!dev) return NULL;

    const md_gpu_image_desc_t desc = image->desc;
    const uint32_t bpp = md_gpu_format_bytes_per_pixel(desc.format);
    if (bpp == 0) return NULL;

    const uint32_t ox = src_region.offset[0];
    const uint32_t oy = src_region.offset[1];
    const uint32_t oz = src_region.offset[2];
    if (ox >= desc.width || oy >= desc.height || oz >= desc.depth) return NULL;

    const uint32_t ex = src_region.extent[0] ? src_region.extent[0] : (desc.width - ox);
    const uint32_t ey = src_region.extent[1] ? src_region.extent[1] : (desc.height - oy);
    const uint32_t ez = src_region.extent[2] ? src_region.extent[2] : (desc.depth - oz);
    if (ex == 0 || ey == 0 || ez == 0) return NULL;
    if (ex > desc.width - ox || ey > desc.height - oy || ez > desc.depth - oz) return NULL;

    const size_t row_bytes = (size_t)ex * (size_t)bpp;
    if (bpp != 0 && row_bytes / (size_t)bpp != (size_t)ex) return NULL;
    if (ey != 0 && row_bytes > SIZE_MAX / (size_t)ey) return NULL;
    const size_t image_bytes = row_bytes * (size_t)ey;
    if (ez != 0 && image_bytes > SIZE_MAX / (size_t)ez) return NULL;
    const size_t total_bytes = image_bytes * (size_t)ez;
    if (total_bytes == 0) return NULL;

    md_mtl_transient_page_t* page = mtl_transient_page_acquire(dev, total_bytes, true);
    if (!page) return NULL;

    struct md_gpu_readback* readback = (struct md_gpu_readback*)calloc(1, sizeof(*readback));
    if (!readback) {
        mtl_transient_page_release(dev, page);
        return NULL;
    }

    md_gpu_buffer_image_copy_t copy = {
        .image_region = src_region,
        .buffer_offset = 0,
        .bytes_per_row = 0,
        .bytes_per_image = 0,
    };

    md_gpu_queue_t queue = md_gpu_queue_transfer((md_gpu_device_t)dev);
    md_gpu_cmd_t cmd = md_gpu_cmd_begin(queue, "readback_image");
    if (!cmd) goto fail;
    if (!md_gpu_cmd_copy_image_to_buffer(cmd, src_image, page->buffer, &copy)) goto fail;
    if (!md_gpu_cmd_end(cmd)) goto fail;

    md_gpu_event_t event = md_gpu_event_is_valid(after)
        ? md_gpu_queue_submit_one_after(queue, cmd, after, MD_GPU_BARRIER_STAGE_TRANSFER)
        : md_gpu_queue_submit_one(queue, cmd);
    if (!md_gpu_event_is_valid(event)) goto fail;

    readback->dev = dev;
    readback->page = page;
    readback->event = event;
    readback->size = total_bytes;
    return (md_gpu_readback_t)readback;

fail:
    if (cmd) md_gpu_cmd_discard(cmd);
    free(readback);
    mtl_transient_page_release(dev, page);
    return NULL;
}

bool md_gpu_readback_is_complete(md_gpu_readback_t readback) {
    if (!readback) return true;
    return md_gpu_event_is_complete(readback->event);
}

void md_gpu_readback_wait(md_gpu_readback_t readback) {
    if (!readback) return;
    md_gpu_event_wait(readback->event);
}

void* md_gpu_readback_cpu_ptr(md_gpu_readback_t readback) {
    if (!readback || !readback->page) return NULL;
    return readback->page->cpu_ptr;
}

void md_gpu_readback_destroy(md_gpu_readback_t readback) {
    if (!readback) return;
    md_gpu_readback_wait(readback);
    if (readback->dev && readback->page) {
        mtl_transient_page_release(readback->dev, readback->page);
        readback->page = NULL;
    }
    free(readback);
}

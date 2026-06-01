#include "md_gpu.h"

#include <core/md_common.h>
#include <core/md_intrinsics.h>
#include <core/md_log.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

// We avoid including platform-specific Vulkan headers via GLFW here.
// This is a compute-only backend for now, so no surface extensions are required.
#define VOLK_IMPLEMENTATION
#include <volk.h>

// Vulkan backend (C).
// Compute-only for now; designed to be extended with graphics later.

// Command buffer pool size. Must be <= 64 to fit the availability bitmask.
#define MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE 64
STATIC_ASSERT(MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE <= 64, "Vulkan command buffer pool must fit in a uint64_t bitmask");

// Maximum number of descriptor set slots available per command buffer recording.
#define MD_GPU_VK_MAX_DISPATCHES_PER_CMD 32
#define MD_GPU_MAX_BIND_SLOTS 16
#define MD_GPU_MAX_RESOURCE_BINDINGS (MD_GPU_MAX_BIND_SLOTS * 3)
#define MD_GPU_VK_MAX_TRACKED_BUFFERS 128
#define MD_GPU_VK_MAX_TRACKED_IMAGES 128
#define MD_GPU_VK_ALL_CMD_SLOTS UINT64_MAX

// Descriptor set indices (one set per resource category)
enum {
    MD_GPU_DESC_SET_STORAGE_IMAGES,
    MD_GPU_DESC_SET_SAMPLED_IMAGES,
    MD_GPU_DESC_SET_SAMPLERS,
    MD_GPU_DESC_SET_COUNT
};

#define MD_GPU_MAX_PUSH_CONSTANTS 256

typedef struct md_gpu_command_buffer* md_gpu_command_buffer_t;

typedef struct md_vk_resource_sync_state_t {
    VkPipelineStageFlags2 stage;
    VkAccessFlags2 access;
    bool valid;
} md_vk_resource_sync_state_t;

typedef struct md_gpu_cmd_resource_binding_t {
    md_gpu_resource_kind_t kind;
    md_gpu_usage_flags_t usage;
    uint32_t set;
    uint32_t binding;
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
} md_gpu_cmd_resource_binding_t;

typedef struct md_vk_cmd_buffer_state_t {
    struct md_gpu_buffer* buffer;
    md_vk_resource_sync_state_t initial_sync;
    md_vk_resource_sync_state_t current_sync;
} md_vk_cmd_buffer_state_t;

typedef struct md_vk_cmd_image_state_t {
    struct md_gpu_image* image;
    VkImageLayout initial_layout;
    VkImageLayout current_layout;
    md_vk_resource_sync_state_t initial_sync;
    md_vk_resource_sync_state_t current_sync;
} md_vk_cmd_image_state_t;

typedef enum md_vk_cmd_state_t {
    MD_VK_CMD_STATE_AVAILABLE = 0,
    MD_VK_CMD_STATE_ACQUIRED,
    MD_VK_CMD_STATE_RECORDING,
    MD_VK_CMD_STATE_ENDED,
    MD_VK_CMD_STATE_SUBMITTED,
} md_vk_cmd_state_t;

typedef struct md_vk_recording_context md_vk_recording_context;

typedef struct md_gpu_command_buffer {
    struct md_gpu_queue*  queue;
    md_vk_recording_context* context;

    VkCommandBuffer vk_cmd;
    VkDescriptorSet (*dispatch_sets)[MD_GPU_DESC_SET_COUNT]; /* pointer into recorder-local descriptor sets */
    uint32_t dispatch_count; /* slots consumed from the pool this recording */
    uint64_t retire_value;
    md_vk_cmd_state_t state;
    bool pushed_debug_group;

    md_gpu_compute_pipeline_t bound_pipeline;
    md_vk_cmd_buffer_state_t tracked_buffers[MD_GPU_VK_MAX_TRACKED_BUFFERS];
    md_vk_cmd_image_state_t tracked_images[MD_GPU_VK_MAX_TRACKED_IMAGES];
    uint32_t tracked_buffer_count;
    uint32_t tracked_image_count;
} md_gpu_command_buffer;

struct md_vk_recording_context {
    struct md_gpu_queue* queue;
    VkCommandPool vk_command_pool;
    VkDescriptorPool descriptor_pool;
    VkDescriptorSet dispatch_sets[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD][MD_GPU_DESC_SET_COUNT];
    md_gpu_command_buffer cmd_wrappers[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    uint64_t available_cmd_slots;
    md_vk_recording_context* next;
    md_vk_recording_context* tls_next;
};

typedef struct md_gpu_queue {
    struct md_gpu_device* dev;
    VkQueue               vk_queue;
    uint32_t              family_index;
    VkSemaphore           event_timeline;
    uint64_t              next_event_id;

    _Atomic(md_vk_recording_context*) recording_contexts;
} md_gpu_queue;

typedef struct md_gpu_device {
    VkInstance instance;
    VkPhysicalDevice physical_device;
    VkDevice device;

    uint32_t compute_queue_family;

    VkDescriptorSetLayout descriptor_set_layouts[MD_GPU_DESC_SET_COUNT];
    VkPipelineLayout pipeline_layout;

    // Embedded primary compute queue used by cmd submission.
    struct md_gpu_queue compute_queue;
} md_gpu_device;

typedef struct md_gpu_buffer {
    md_gpu_device_t device;
    VkBuffer buffer;
    VkDeviceMemory memory;
    VkDeviceSize size;
    md_gpu_buffer_flags_t flags;
    void* cpu_ptr;  /* non-NULL for CPU_VISIBLE buffers; valid until destroy */
    md_vk_resource_sync_state_t sync;
} md_gpu_buffer;

typedef struct md_gpu_image {
    md_gpu_device_t device;
    VkImage image;
    VkImageView view;
    VkDeviceMemory memory;
    uint32_t width;
    uint32_t height;
    uint32_t depth;
    VkFormat format;
    md_gpu_image_flags_t flags;
    VkImageLayout layout;
    md_vk_resource_sync_state_t sync;
} md_gpu_image;

typedef struct md_gpu_sampler {
    md_gpu_device_t device;
    VkSampler sampler;
} md_gpu_sampler;

typedef struct md_gpu_compute_pipeline {
    md_gpu_device_t device;
    VkPipeline pipeline;
    VkShaderModule shader_module;
    md_gpu_resource_binding_t resource_bindings[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint32_t resource_binding_count;
} md_gpu_compute_pipeline;

static bool md_vk_check(VkResult res, const char* what) {
    if (res == VK_SUCCESS) return true;
    MD_LOG_ERROR("Vulkan error %d: %s", (int)res, what ? what : "(unknown)");
    return false;
}

static VkAccessFlags2 md_vk_access_from_gpu_usage(md_gpu_usage_flags_t usage) {
    VkAccessFlags2 access = 0;
    if (usage & MD_GPU_USAGE_READ)  access |= VK_ACCESS_2_SHADER_READ_BIT;
    if (usage & MD_GPU_USAGE_WRITE) access |= VK_ACCESS_2_SHADER_WRITE_BIT;
    return access;
}

static bool md_vk_access_has_write(VkAccessFlags2 access) {
    return (access & (VK_ACCESS_2_SHADER_WRITE_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT)) != 0;
}

static bool md_vk_sync_state_equal(md_vk_resource_sync_state_t a, md_vk_resource_sync_state_t b) {
    return a.stage == b.stage && a.access == b.access && a.valid == b.valid;
}

static md_vk_cmd_buffer_state_t* md_vk_cmd_track_buffer(md_gpu_command_buffer* cmd, md_gpu_buffer* buf) {
    if (!cmd || !buf) return NULL;

    for (uint32_t i = 0; i < cmd->tracked_buffer_count; ++i) {
        md_vk_cmd_buffer_state_t* state = &cmd->tracked_buffers[i];
        if (state->buffer == buf) return state;
    }

    if (cmd->tracked_buffer_count >= MD_GPU_VK_MAX_TRACKED_BUFFERS) {
        ASSERT(false && "Too many tracked buffers in Vulkan command buffer");
        return NULL;
    }

    md_vk_cmd_buffer_state_t* state = &cmd->tracked_buffers[cmd->tracked_buffer_count++];
    state->buffer = buf;
    state->initial_sync = buf->sync;
    state->current_sync = buf->sync;
    return state;
}

static md_vk_cmd_image_state_t* md_vk_cmd_track_image(md_gpu_command_buffer* cmd, md_gpu_image* img) {
    if (!cmd || !img) return NULL;

    for (uint32_t i = 0; i < cmd->tracked_image_count; ++i) {
        md_vk_cmd_image_state_t* state = &cmd->tracked_images[i];
        if (state->image == img) return state;
    }

    if (cmd->tracked_image_count >= MD_GPU_VK_MAX_TRACKED_IMAGES) {
        ASSERT(false && "Too many tracked images in Vulkan command buffer");
        return NULL;
    }

    md_vk_cmd_image_state_t* state = &cmd->tracked_images[cmd->tracked_image_count++];
    state->image = img;
    state->initial_layout = img->layout;
    state->current_layout = img->layout;
    state->initial_sync = img->sync;
    state->current_sync = img->sync;
    return state;
}

static VkImageLayout md_vk_cmd_image_layout(const md_gpu_command_buffer* cmd, const md_gpu_image* img) {
    if (!img) return VK_IMAGE_LAYOUT_UNDEFINED;
    if (cmd) {
        for (uint32_t i = 0; i < cmd->tracked_image_count; ++i) {
            const md_vk_cmd_image_state_t* state = &cmd->tracked_images[i];
            if (state->image == img) return state->current_layout;
        }
    }

    return img->layout;
}

static const md_vk_resource_sync_state_t* md_vk_batch_find_prior_buffer_sync(const md_gpu_cmd_t* cmds, size_t count, md_gpu_buffer* buf) {
    for (size_t i = count; i-- > 0;) {
        const md_gpu_command_buffer* cmd = (const md_gpu_command_buffer*)cmds[i];
        if (!cmd) continue;
        for (uint32_t j = 0; j < cmd->tracked_buffer_count; ++j) {
            const md_vk_cmd_buffer_state_t* state = &cmd->tracked_buffers[j];
            if (state->buffer == buf) return &state->current_sync;
        }
    }
    return NULL;
}

static const md_vk_cmd_image_state_t* md_vk_batch_find_prior_image_state(const md_gpu_cmd_t* cmds, size_t count, md_gpu_image* img) {
    for (size_t i = count; i-- > 0;) {
        const md_gpu_command_buffer* cmd = (const md_gpu_command_buffer*)cmds[i];
        if (!cmd) continue;
        for (uint32_t j = 0; j < cmd->tracked_image_count; ++j) {
            const md_vk_cmd_image_state_t* state = &cmd->tracked_images[j];
            if (state->image == img) return state;
        }
    }
    return NULL;
}

static bool md_vk_cmd_validate_initial_state(const md_gpu_cmd_t* cmds, size_t cmd_index, const md_gpu_command_buffer* cmd) {
    if (!cmd) return false;

    for (uint32_t i = 0; i < cmd->tracked_buffer_count; ++i) {
        const md_vk_cmd_buffer_state_t* state = &cmd->tracked_buffers[i];
        const md_vk_resource_sync_state_t* expected = md_vk_batch_find_prior_buffer_sync(cmds, cmd_index, state->buffer);
        md_vk_resource_sync_state_t current = expected ? *expected : state->buffer->sync;
        if (!md_vk_sync_state_equal(state->initial_sync, current)) {
            MD_LOG_ERROR("Vulkan: command buffer recorded against stale buffer state; concurrent writable sharing is not supported yet");
            return false;
        }
    }

    for (uint32_t i = 0; i < cmd->tracked_image_count; ++i) {
        const md_vk_cmd_image_state_t* state = &cmd->tracked_images[i];
        const md_vk_cmd_image_state_t* expected = md_vk_batch_find_prior_image_state(cmds, cmd_index, state->image);
        const VkImageLayout layout = expected ? expected->current_layout : state->image->layout;
        const md_vk_resource_sync_state_t sync = expected ? expected->current_sync : state->image->sync;
        if (state->initial_layout != layout || !md_vk_sync_state_equal(state->initial_sync, sync)) {
            MD_LOG_ERROR("Vulkan: command buffer recorded against stale image state; concurrent writable sharing is not supported yet");
            return false;
        }
    }

    return true;
}

static void md_vk_cmd_apply_final_state(const md_gpu_command_buffer* cmd) {
    if (!cmd) return;

    for (uint32_t i = 0; i < cmd->tracked_buffer_count; ++i) {
        const md_vk_cmd_buffer_state_t* state = &cmd->tracked_buffers[i];
        state->buffer->sync = state->current_sync;
    }

    for (uint32_t i = 0; i < cmd->tracked_image_count; ++i) {
        const md_vk_cmd_image_state_t* state = &cmd->tracked_images[i];
        state->image->layout = state->current_layout;
        state->image->sync = state->current_sync;
    }
}

static void md_vk_ensure_buffer_sync(md_gpu_command_buffer* cmd, md_gpu_buffer* buf, VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    if (!cmd || !buf || dst_stage == 0 || dst_access == 0) return;

    md_vk_cmd_buffer_state_t* state = md_vk_cmd_track_buffer(cmd, buf);
    if (!state) return;

    const bool need_barrier =
        state->current_sync.valid &&
        (md_vk_access_has_write(state->current_sync.access) || md_vk_access_has_write(dst_access));

    if (need_barrier) {
        VkBufferMemoryBarrier2 barrier = {0};
        barrier.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER_2;
        barrier.srcStageMask = state->current_sync.stage;
        barrier.srcAccessMask = state->current_sync.access;
        barrier.dstStageMask = dst_stage;
        barrier.dstAccessMask = dst_access;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.buffer = buf->buffer;
        barrier.offset = 0;
        barrier.size = VK_WHOLE_SIZE;

        VkDependencyInfo dep = {0};
        dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        dep.bufferMemoryBarrierCount = 1;
        dep.pBufferMemoryBarriers = &barrier;
        vkCmdPipelineBarrier2(cmd->vk_cmd, &dep);
    }

    state->current_sync.stage = dst_stage;
    state->current_sync.access = dst_access;
    state->current_sync.valid = true;
}

static void md_vk_ensure_image_sync(md_gpu_command_buffer* cmd, md_gpu_image* img, VkImageLayout new_layout, VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    if (!cmd || !img || dst_stage == 0 || dst_access == 0) return;

    md_vk_cmd_image_state_t* state = md_vk_cmd_track_image(cmd, img);
    if (!state) return;

    const bool layout_change = (state->current_layout != new_layout);
    const bool hazard = state->current_sync.valid && (md_vk_access_has_write(state->current_sync.access) || md_vk_access_has_write(dst_access));
    const bool need_barrier = layout_change || hazard;

    if (need_barrier) {
        VkImageMemoryBarrier2 barrier = {0};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
        barrier.srcStageMask = state->current_sync.valid ? state->current_sync.stage : VK_PIPELINE_STAGE_2_TOP_OF_PIPE_BIT;
        barrier.srcAccessMask = state->current_sync.valid ? state->current_sync.access : 0;
        barrier.dstStageMask = dst_stage;
        barrier.dstAccessMask = dst_access;
        barrier.oldLayout = state->current_layout;
        barrier.newLayout = new_layout;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = img->image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;

        VkDependencyInfo dep = {0};
        dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
        dep.imageMemoryBarrierCount = 1;
        dep.pImageMemoryBarriers = &barrier;
        vkCmdPipelineBarrier2(cmd->vk_cmd, &dep);
    }

    state->current_layout = new_layout;
    state->current_sync.stage = dst_stage;
    state->current_sync.access = dst_access;
    state->current_sync.valid = true;
}

typedef struct {
    md_gpu_cmd_resource_binding_t bindings[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint32_t count;
} md_vk_dispatch_resources_t;

static bool md_vk_cmd_merge_resource_binding(md_vk_dispatch_resources_t* res, const md_gpu_resource_t* resource) {
    if (!res || !resource) return false;

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

    for (uint32_t i = 0; i < res->count; ++i) {
        md_gpu_cmd_resource_binding_t* existing = &res->bindings[i];
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

    if (res->count >= MD_GPU_MAX_RESOURCE_BINDINGS) return false;

    md_gpu_cmd_resource_binding_t* binding = &res->bindings[res->count++];
    binding->kind = resource->kind;
    binding->usage = usage;
    binding->set = resource->set;
    binding->binding = resource->binding;
    binding->image = resource->image;
    binding->sampler = resource->sampler;
    return true;
}

static void md_vk_reset_cmd_state(md_gpu_command_buffer* cmd) {
    cmd->dispatch_count = 0;
    cmd->retire_value = 0;
    cmd->state = MD_VK_CMD_STATE_AVAILABLE;
    cmd->pushed_debug_group = false;
    cmd->bound_pipeline = NULL;
    cmd->tracked_buffer_count = 0;
    cmd->tracked_image_count = 0;
}

static uint32_t md_vk_cmd_slot_index(const md_vk_recording_context* ctx, const md_gpu_command_buffer* cmd) {
    ASSERT(ctx);
    ASSERT(cmd);
    ASSERT(cmd >= ctx->cmd_wrappers);
    ASSERT(cmd < ctx->cmd_wrappers + MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE);
    return (uint32_t)(cmd - ctx->cmd_wrappers);
}

static uint64_t md_vk_cmd_slot_mask(uint32_t idx) {
    ASSERT(idx < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE);
    return UINT64_C(1) << idx;
}

static void md_vk_release_cmd_slot(md_gpu_command_buffer* cmd) {
    ASSERT(cmd);
    ASSERT(cmd->context);

    md_vk_recording_context* ctx = cmd->context;
    const uint32_t idx = md_vk_cmd_slot_index(ctx, cmd);
    const uint64_t slot_mask = md_vk_cmd_slot_mask(idx);
    ASSERT((ctx->available_cmd_slots & slot_mask) == 0);

    md_vk_reset_cmd_state(cmd);
    ctx->available_cmd_slots |= slot_mask;
}

static bool md_vk_has_extension(uint32_t count, const VkExtensionProperties* props, const char* name) {
    for (uint32_t i = 0; i < count; ++i) {
        if (strcmp(props[i].extensionName, name) == 0) return true;
    }
    return false;
}

static uint32_t md_vk_find_memory_type(VkPhysicalDevice phys, uint32_t type_bits, VkMemoryPropertyFlags required, VkMemoryPropertyFlags preferred) {
    VkPhysicalDeviceMemoryProperties mem_props;
    vkGetPhysicalDeviceMemoryProperties(phys, &mem_props);

    uint32_t best = UINT32_MAX;
    for (uint32_t i = 0; i < mem_props.memoryTypeCount; ++i) {
        if ((type_bits & (1u << i)) == 0) continue;
        VkMemoryPropertyFlags flags = mem_props.memoryTypes[i].propertyFlags;
        if ((flags & required) != required) continue;

        if (preferred == 0) return i;
        if ((flags & preferred) == preferred) return i;
        if (best == UINT32_MAX) best = i;
    }
    return best;
}

static VkFormat md_vk_format_from_md(md_gpu_image_format_t fmt) {
    switch (fmt) {
        case MD_GPU_IMAGE_FORMAT_R8_UINT:       return VK_FORMAT_R8_UINT;
        case MD_GPU_IMAGE_FORMAT_R16_UINT:      return VK_FORMAT_R16_UINT;
        case MD_GPU_IMAGE_FORMAT_R32_UINT:      return VK_FORMAT_R32_UINT;
        case MD_GPU_IMAGE_FORMAT_R32_FLOAT:     return VK_FORMAT_R32_SFLOAT;
        case MD_GPU_IMAGE_FORMAT_RGBA8_UNORM:   return VK_FORMAT_R8G8B8A8_UNORM;
        case MD_GPU_IMAGE_FORMAT_RGBA16_FLOAT:  return VK_FORMAT_R16G16B16A16_SFLOAT;
        case MD_GPU_IMAGE_FORMAT_RGBA32_FLOAT:  return VK_FORMAT_R32G32B32A32_SFLOAT;
        default: return VK_FORMAT_UNDEFINED;
    }
}

static md_gpu_image_format_t md_vk_format_to_md(VkFormat fmt) {
    switch (fmt) {
        case VK_FORMAT_R8_UINT:       return MD_GPU_IMAGE_FORMAT_R8_UINT;
        case VK_FORMAT_R16_UINT:      return MD_GPU_IMAGE_FORMAT_R16_UINT;
        case VK_FORMAT_R32_UINT:      return MD_GPU_IMAGE_FORMAT_R32_UINT;
        case VK_FORMAT_R32_SFLOAT:    return MD_GPU_IMAGE_FORMAT_R32_FLOAT;
        case VK_FORMAT_R8G8B8A8_UNORM: return MD_GPU_IMAGE_FORMAT_RGBA8_UNORM;
        case VK_FORMAT_R16G16B16A16_SFLOAT: return MD_GPU_IMAGE_FORMAT_RGBA16_FLOAT;
        case VK_FORMAT_R32G32B32A32_SFLOAT: return MD_GPU_IMAGE_FORMAT_RGBA32_FLOAT;
        default: return MD_GPU_IMAGE_FORMAT_R8_UINT;
    }
}

static uint32_t md_vk_format_bytes_per_pixel(VkFormat fmt) {
    switch (fmt) {
        case VK_FORMAT_R8_UINT: return 1;
        case VK_FORMAT_R16_UINT: return 2;
        case VK_FORMAT_R32_UINT:
        case VK_FORMAT_R32_SFLOAT:
        case VK_FORMAT_R8G8B8A8_UNORM: return 4;
        case VK_FORMAT_R16G16B16A16_SFLOAT: return 8;
        case VK_FORMAT_R32G32B32A32_SFLOAT: return 16;
        default: return 0;
    }
}

static VkFilter md_vk_filter_from_md(md_gpu_filter_t filter) {
    switch (filter) {
        case MD_GPU_FILTER_NEAREST: return VK_FILTER_NEAREST;
        case MD_GPU_FILTER_LINEAR:  return VK_FILTER_LINEAR;
        default: return VK_FILTER_NEAREST;
    }
}

static VkSamplerAddressMode md_vk_address_mode_from_md(md_gpu_address_mode_t mode) {
    switch (mode) {
        case MD_GPU_ADDRESS_MODE_CLAMP_TO_EDGE:    return VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
        case MD_GPU_ADDRESS_MODE_REPEAT:           return VK_SAMPLER_ADDRESS_MODE_REPEAT;
        case MD_GPU_ADDRESS_MODE_MIRRORED_REPEAT:  return VK_SAMPLER_ADDRESS_MODE_MIRRORED_REPEAT;
        default: return VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE;
    }
}

static bool md_vk_pick_physical_device(md_gpu_device* out_dev) {
    uint32_t phys_count = 0;
    if (!md_vk_check(vkEnumeratePhysicalDevices(out_dev->instance, &phys_count, NULL), "vkEnumeratePhysicalDevices(count)")) return false;
    if (phys_count == 0) {
        MD_LOG_ERROR("Vulkan: no physical devices found");
        return false;
    }

    VkPhysicalDevice* phys = (VkPhysicalDevice*)malloc(sizeof(VkPhysicalDevice) * phys_count);
    if (!phys) return false;
    bool ok = md_vk_check(vkEnumeratePhysicalDevices(out_dev->instance, &phys_count, phys), "vkEnumeratePhysicalDevices(list)");
    if (!ok) {
        free(phys);
        return false;
    }

    // Prefer discrete GPUs, but accept anything with a compute queue.
    VkPhysicalDevice best = VK_NULL_HANDLE;
    uint32_t best_compute_family = UINT32_MAX;
    int best_score = -1;

    for (uint32_t p = 0; p < phys_count; ++p) {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(phys[p], &props);

        uint32_t qcount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(phys[p], &qcount, NULL);
        if (qcount == 0) continue;
        VkQueueFamilyProperties* qprops = (VkQueueFamilyProperties*)malloc(sizeof(VkQueueFamilyProperties) * qcount);
        if (!qprops) continue;
        vkGetPhysicalDeviceQueueFamilyProperties(phys[p], &qcount, qprops);

        uint32_t compute_family = UINT32_MAX;
        for (uint32_t qi = 0; qi < qcount; ++qi) {
            if ((qprops[qi].queueFlags & VK_QUEUE_COMPUTE_BIT) != 0) {
                compute_family = qi;
                // Prefer a compute-only queue if present.
                if ((qprops[qi].queueFlags & VK_QUEUE_GRAPHICS_BIT) == 0) break;
            }
        }
        free(qprops);
        if (compute_family == UINT32_MAX) continue;

        int score = 0;
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) score += 100;
        if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU) score += 50;
        score += (int)props.limits.maxComputeSharedMemorySize / 1024;

        if (score > best_score) {
            best_score = score;
            best = phys[p];
            best_compute_family = compute_family;
        }
    }

    free(phys);

    if (best == VK_NULL_HANDLE) {
        MD_LOG_ERROR("Vulkan: no suitable physical device with compute queue found");
        return false;
    }

    out_dev->physical_device = best;
    out_dev->compute_queue_family = best_compute_family;
    return true;
}

static bool md_vk_create_instance(md_gpu_device* out_dev) {
    // For compute-only we don't need platform surface extensions.
    const char* instance_exts[8];
    uint32_t instance_ext_count = 0;

    // Enable portability enumeration if supported (helps MoltenVK, harmless elsewhere).
    uint32_t avail_count = 0;
    vkEnumerateInstanceExtensionProperties(NULL, &avail_count, NULL);
    VkExtensionProperties* avail = NULL;
    if (avail_count) {
        avail = (VkExtensionProperties*)malloc(sizeof(VkExtensionProperties) * avail_count);
        if (avail) vkEnumerateInstanceExtensionProperties(NULL, &avail_count, avail);
    }

    VkInstanceCreateFlags flags = 0;
    if (avail && md_vk_has_extension(avail_count, avail, VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME)) {
        instance_exts[instance_ext_count++] = VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME;
        flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
    }

    if (avail && md_vk_has_extension(avail_count, avail, VK_EXT_DEBUG_UTILS_EXTENSION_NAME)) {
        instance_exts[instance_ext_count++] = VK_EXT_DEBUG_UTILS_EXTENSION_NAME;
    }

    free(avail);

    VkApplicationInfo app = {0};
    app.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app.pApplicationName = "md_gpu";
    app.applicationVersion = VK_MAKE_VERSION(0, 1, 0);
    app.pEngineName = "md_gpu";
    app.engineVersion = VK_MAKE_VERSION(0, 1, 0);
    // Request a modern desktop baseline. Driver may return a lower supported version.
    app.apiVersion = VK_API_VERSION_1_3;

    VkInstanceCreateInfo ci = {0};
    ci.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    ci.flags = flags;
    ci.pApplicationInfo = &app;
    ci.enabledExtensionCount = instance_ext_count;
    ci.ppEnabledExtensionNames = (instance_ext_count ? instance_exts : NULL);

    return md_vk_check(vkCreateInstance(&ci, NULL, &out_dev->instance), "vkCreateInstance");
}

static void md_vk_destroy_recording_context(md_gpu_device* dev, md_vk_recording_context* ctx) {
    if (!dev || !ctx) return;
    if (ctx->vk_command_pool) vkDestroyCommandPool(dev->device, ctx->vk_command_pool, NULL);
    if (ctx->descriptor_pool) vkDestroyDescriptorPool(dev->device, ctx->descriptor_pool, NULL);
    free(ctx);
}

static THREAD_LOCAL md_vk_recording_context* md_vk_tls_recording_contexts;

static void md_vk_queue_register_recording_context(md_gpu_queue* queue, md_vk_recording_context* ctx) {
    md_vk_recording_context* head = atomic_load_explicit(&queue->recording_contexts, memory_order_acquire);
    do {
        ctx->next = head;
    } while (!atomic_compare_exchange_weak_explicit(&queue->recording_contexts, &head, ctx, memory_order_release, memory_order_acquire));
}

static md_vk_recording_context* md_vk_create_recording_context(md_gpu_queue* queue) {
    if (!queue || !queue->dev) return NULL;
    md_gpu_device* dev = queue->dev;

    md_vk_recording_context* ctx = (md_vk_recording_context*)calloc(1, sizeof(md_vk_recording_context));
    if (!ctx) return NULL;

    ctx->queue = queue;
    ctx->available_cmd_slots = MD_GPU_VK_ALL_CMD_SLOTS;

    VkCommandPoolCreateInfo pool_ci = {0};
    pool_ci.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    pool_ci.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pool_ci.queueFamilyIndex = queue->family_index;
    if (!md_vk_check(vkCreateCommandPool(dev->device, &pool_ci, NULL, &ctx->vk_command_pool), "vkCreateCommandPool(recording context)")) {
        md_vk_destroy_recording_context(dev, ctx);
        return NULL;
    }

    VkCommandBuffer vk_cmds[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    VkCommandBufferAllocateInfo cb_ai = {0};
    cb_ai.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    cb_ai.commandPool = ctx->vk_command_pool;
    cb_ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cb_ai.commandBufferCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;
    if (!md_vk_check(vkAllocateCommandBuffers(dev->device, &cb_ai, vk_cmds), "vkAllocateCommandBuffers(recording context)")) {
        md_vk_destroy_recording_context(dev, ctx);
        return NULL;
    }

    VkDescriptorPoolSize pool_sizes[MD_GPU_DESC_SET_COUNT];
    pool_sizes[0].type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    pool_sizes[0].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[1].type = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
    pool_sizes[1].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[2].type = VK_DESCRIPTOR_TYPE_SAMPLER;
    pool_sizes[2].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;

    VkDescriptorPoolCreateInfo dp_ci = {0};
    dp_ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    dp_ci.maxSets = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_DESC_SET_COUNT;
    dp_ci.poolSizeCount = MD_GPU_DESC_SET_COUNT;
    dp_ci.pPoolSizes = pool_sizes;
    if (!md_vk_check(vkCreateDescriptorPool(dev->device, &dp_ci, NULL, &ctx->descriptor_pool), "vkCreateDescriptorPool(recording context)")) {
        md_vk_destroy_recording_context(dev, ctx);
        return NULL;
    }

    const uint32_t dispatches_total = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD;
    const uint32_t total_sets = dispatches_total * MD_GPU_DESC_SET_COUNT;
    VkDescriptorSetLayout* set_layouts = (VkDescriptorSetLayout*)malloc(total_sets * sizeof(VkDescriptorSetLayout));
    VkDescriptorSet* all_sets = (VkDescriptorSet*)malloc(total_sets * sizeof(VkDescriptorSet));
    if (!set_layouts || !all_sets) {
        free(set_layouts);
        free(all_sets);
        md_vk_destroy_recording_context(dev, ctx);
        return NULL;
    }
    for (uint32_t i = 0; i < dispatches_total; ++i) {
        for (uint32_t s = 0; s < MD_GPU_DESC_SET_COUNT; ++s) {
            set_layouts[i * MD_GPU_DESC_SET_COUNT + s] = dev->descriptor_set_layouts[s];
        }
    }

    VkDescriptorSetAllocateInfo ds_ai = {0};
    ds_ai.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    ds_ai.descriptorPool = ctx->descriptor_pool;
    ds_ai.descriptorSetCount = total_sets;
    ds_ai.pSetLayouts = set_layouts;
    const bool sets_ok = md_vk_check(vkAllocateDescriptorSets(dev->device, &ds_ai, all_sets), "vkAllocateDescriptorSets(recording context)");
    if (sets_ok) {
        for (uint32_t i = 0; i < dispatches_total; ++i) {
            for (uint32_t s = 0; s < MD_GPU_DESC_SET_COUNT; ++s) {
                ctx->dispatch_sets[i][s] = all_sets[i * MD_GPU_DESC_SET_COUNT + s];
            }
        }
    }
    free(set_layouts);
    free(all_sets);
    if (!sets_ok) {
        md_vk_destroy_recording_context(dev, ctx);
        return NULL;
    }

    for (uint32_t i = 0; i < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE; ++i) {
        ctx->cmd_wrappers[i].queue = queue;
        ctx->cmd_wrappers[i].context = ctx;
        ctx->cmd_wrappers[i].vk_cmd = vk_cmds[i];
        ctx->cmd_wrappers[i].dispatch_sets = &ctx->dispatch_sets[i * MD_GPU_VK_MAX_DISPATCHES_PER_CMD];
    }

    return ctx;
}

static md_vk_recording_context* md_vk_queue_context_for_current_thread(md_gpu_queue* queue) {
    if (!queue) return NULL;

    for (md_vk_recording_context* ctx = md_vk_tls_recording_contexts; ctx; ctx = ctx->tls_next) {
        if (ctx->queue == queue) return ctx;
    }

    md_vk_recording_context* ctx = md_vk_create_recording_context(queue);
    if (!ctx) return NULL;
    ctx->tls_next = md_vk_tls_recording_contexts;
    md_vk_tls_recording_contexts = ctx;
    md_vk_queue_register_recording_context(queue, ctx);
    return ctx;
}

static bool md_vk_create_device(md_gpu_device* out_dev) {
    float prio = 1.0f;
    VkDeviceQueueCreateInfo qci = {0};
    qci.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    qci.queueFamilyIndex = out_dev->compute_queue_family;
    qci.queueCount = 1;
    qci.pQueuePriorities = &prio;

    VkPhysicalDeviceVulkan12Features f12 = {0};
    f12.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES;
    f12.bufferDeviceAddress = VK_TRUE;
    f12.timelineSemaphore = VK_TRUE;
    f12.descriptorBindingPartiallyBound = VK_TRUE;

    VkPhysicalDeviceVulkan13Features f13 = {0};
    f13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
    f13.pNext = &f12;
    f13.synchronization2 = VK_TRUE;

    VkPhysicalDeviceFeatures2 feats2 = {0};
    feats2.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2;
    feats2.pNext = &f13;

    VkDeviceCreateInfo dci = {0};
    dci.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos = &qci;
    dci.pNext = &feats2;

    if (!md_vk_check(vkCreateDevice(out_dev->physical_device, &dci, NULL, &out_dev->device), "vkCreateDevice")) return false;

    md_gpu_queue* q = &out_dev->compute_queue;
    vkGetDeviceQueue(out_dev->device, out_dev->compute_queue_family, 0, &q->vk_queue);
    q->dev          = out_dev;
    q->family_index = out_dev->compute_queue_family;
    if (q->vk_queue == VK_NULL_HANDLE) return false;

    // Create three descriptor set layouts, one per resource category:
    //   set 0 = storage images, set 1 = sampled images, set 2 = samplers
    {
        static const VkDescriptorType set_types[MD_GPU_DESC_SET_COUNT] = {
            VK_DESCRIPTOR_TYPE_STORAGE_IMAGE,
            VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE,
            VK_DESCRIPTOR_TYPE_SAMPLER,
        };
        VkDescriptorSetLayoutBinding slot_bindings[MD_GPU_MAX_BIND_SLOTS];
        VkDescriptorBindingFlags binding_flags[MD_GPU_MAX_BIND_SLOTS];
        for (uint32_t s = 0; s < MD_GPU_DESC_SET_COUNT; ++s) {
            for (uint32_t i = 0; i < MD_GPU_MAX_BIND_SLOTS; ++i) {
                slot_bindings[i].binding            = i;
                slot_bindings[i].descriptorType     = set_types[s];
                slot_bindings[i].descriptorCount    = 1;
                slot_bindings[i].stageFlags         = VK_SHADER_STAGE_COMPUTE_BIT;
                slot_bindings[i].pImmutableSamplers = NULL;
                binding_flags[i]                   = VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT;
            }
            VkDescriptorSetLayoutBindingFlagsCreateInfo flags_ci = {0};
            flags_ci.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO;
            flags_ci.bindingCount  = MD_GPU_MAX_BIND_SLOTS;
            flags_ci.pBindingFlags = binding_flags;
            VkDescriptorSetLayoutCreateInfo dsl_ci = {0};
            dsl_ci.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
            dsl_ci.pNext        = &flags_ci;
            dsl_ci.bindingCount = MD_GPU_MAX_BIND_SLOTS;
            dsl_ci.pBindings    = slot_bindings;
            if (!md_vk_check(vkCreateDescriptorSetLayout(out_dev->device, &dsl_ci, NULL, &out_dev->descriptor_set_layouts[s]), "vkCreateDescriptorSetLayout")) return false;
        }
    }

    // Create pipeline layout with all three sets and push constants
    VkPushConstantRange push_range = {0};
    push_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    push_range.offset = 0;
    push_range.size = MD_GPU_MAX_PUSH_CONSTANTS;

    VkPipelineLayoutCreateInfo pl_ci = {0};
    pl_ci.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    pl_ci.setLayoutCount = MD_GPU_DESC_SET_COUNT;
    pl_ci.pSetLayouts = out_dev->descriptor_set_layouts;
    pl_ci.pushConstantRangeCount = 1;
    pl_ci.pPushConstantRanges = &push_range;
    if (!md_vk_check(vkCreatePipelineLayout(out_dev->device, &pl_ci, NULL, &out_dev->pipeline_layout), "vkCreatePipelineLayout")) return false;

    VkSemaphoreTypeCreateInfo timeline_ci = {0};
    timeline_ci.sType = VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO;
    timeline_ci.semaphoreType = VK_SEMAPHORE_TYPE_TIMELINE;
    timeline_ci.initialValue = 0;

    VkSemaphoreCreateInfo sem_ci = {0};
    sem_ci.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
    sem_ci.pNext = &timeline_ci;
    if (!md_vk_check(vkCreateSemaphore(out_dev->device, &sem_ci, NULL, &q->event_timeline), "vkCreateSemaphore(event timeline)")) return false;

    return true;
}

md_gpu_device_t md_gpu_device_create(void) {
    if (volkInitialize() != VK_SUCCESS) {
        return NULL;
    }

    md_gpu_device* dev = (md_gpu_device*)calloc(1, sizeof(md_gpu_device));
    if (!dev) return NULL;

    if (!md_vk_create_instance(dev)) {
        free(dev);
        return NULL;
    }

    volkLoadInstance(dev->instance);

    if (!md_vk_pick_physical_device(dev)) {
        vkDestroyInstance(dev->instance, NULL);
        free(dev);
        return NULL;
    }

    if (!md_vk_create_device(dev)) {
        vkDestroyInstance(dev->instance, NULL);
        free(dev);
        return NULL;
    }

    volkLoadDevice(dev->device);

    return (md_gpu_device_t)dev;
}

void md_gpu_device_destroy(md_gpu_device_t device) {
    if (!device) return;
    md_gpu_device* dev = (md_gpu_device*)device;
    if (dev->device) {
        vkDeviceWaitIdle(dev->device);
        md_vk_recording_context* ctx = atomic_load_explicit(&dev->compute_queue.recording_contexts, memory_order_acquire);
        while (ctx) {
            md_vk_recording_context* next = ctx->next;
            md_vk_destroy_recording_context(dev, ctx);
            ctx = next;
        }
        if (dev->compute_queue.event_timeline) vkDestroySemaphore(dev->device, dev->compute_queue.event_timeline, NULL);
        if (dev->pipeline_layout) vkDestroyPipelineLayout(dev->device, dev->pipeline_layout, NULL);
        for (uint32_t s = 0; s < MD_GPU_DESC_SET_COUNT; ++s) {
            if (dev->descriptor_set_layouts[s]) vkDestroyDescriptorSetLayout(dev->device, dev->descriptor_set_layouts[s], NULL);
        }
        vkDestroyDevice(dev->device, NULL);
    }
    if (dev->instance) {
        vkDestroyInstance(dev->instance, NULL);
    }
    free(dev);
}

md_gpu_queue_t md_gpu_queue_graphics(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&((md_gpu_device*)device)->compute_queue;
}

md_gpu_queue_t md_gpu_queue_compute(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&((md_gpu_device*)device)->compute_queue;
}

md_gpu_queue_t md_gpu_queue_transfer(md_gpu_device_t device) {
    if (!device) return NULL;
    return (md_gpu_queue_t)&((md_gpu_device*)device)->compute_queue;
}

bool md_gpu_device_info(md_gpu_device_t device, md_gpu_device_info_t* info) {
    if (!device || !info) return false;
    md_gpu_device* dev = (md_gpu_device*)device;

    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(dev->physical_device, &props);
    size_t name_len = strlen(props.deviceName);
    if (name_len >= sizeof(info->name)) name_len = sizeof(info->name) - 1;
    memcpy(info->name, props.deviceName, name_len);
    info->name[name_len] = '\0';

    // Discrete GPUs have dedicated VRAM separate from host memory and require staging.
    // Integrated GPUs (UMA) share memory and can skip staging.
    // deviceType is the authoritative Vulkan signal for this distinction.
    info->is_discrete = (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU);

    return true;
}

md_gpu_buffer_t md_gpu_buffer_create(md_gpu_device_t device, const md_gpu_buffer_desc_t* desc) {
    if (!device || !desc || desc->size == 0) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;

    md_gpu_buffer* buf = (md_gpu_buffer*)calloc(1, sizeof(md_gpu_buffer));
    if (!buf) return NULL;
    buf->device = device;
    buf->size = (VkDeviceSize)desc->size;
    buf->flags = desc->flags;

    // Conservative usage: allow compute shaders and transfer ops.
    VkBufferCreateInfo bci = {0};
    bci.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bci.size = buf->size;
    bci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    bci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    if (!md_vk_check(vkCreateBuffer(dev->device, &bci, NULL, &buf->buffer), "vkCreateBuffer")) {
        free(buf);
        return NULL;
    }

    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(dev->device, buf->buffer, &req);

    VkMemoryPropertyFlags required = 0;
    VkMemoryPropertyFlags preferred = 0;
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE) {
        required |= VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    } else {
        preferred |= VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    }

    uint32_t mem_type = md_vk_find_memory_type(dev->physical_device, req.memoryTypeBits, required, preferred);
    if (mem_type == UINT32_MAX) {
        MD_LOG_ERROR("Vulkan: failed to find memory type for buffer");
        vkDestroyBuffer(dev->device, buf->buffer, NULL);
        free(buf);
        return NULL;
    }

    VkMemoryAllocateFlagsInfo alloc_flags = {0};
    alloc_flags.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO;
    alloc_flags.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;

    VkMemoryAllocateInfo mai = {0};
    mai.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    mai.pNext = &alloc_flags;
    mai.allocationSize = req.size;
    mai.memoryTypeIndex = mem_type;

    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &buf->memory), "vkAllocateMemory")) {
        vkDestroyBuffer(dev->device, buf->buffer, NULL);
        free(buf);
        return NULL;
    }

    if (!md_vk_check(vkBindBufferMemory(dev->device, buf->buffer, buf->memory, 0), "vkBindBufferMemory")) {
        vkFreeMemory(dev->device, buf->memory, NULL);
        vkDestroyBuffer(dev->device, buf->buffer, NULL);
        free(buf);
        return NULL;
    }

    // Persistent map for CPU-visible buffers.
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE) {
        if (!md_vk_check(vkMapMemory(dev->device, buf->memory, 0, VK_WHOLE_SIZE, 0, &buf->cpu_ptr), "vkMapMemory")) {
            vkFreeMemory(dev->device, buf->memory, NULL);
            vkDestroyBuffer(dev->device, buf->buffer, NULL);
            free(buf);
            return NULL;
        }
    }

    return (md_gpu_buffer_t)buf;
}

void md_gpu_buffer_destroy(md_gpu_buffer_t buffer) {
    if (!buffer) return;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    md_gpu_device* dev = (md_gpu_device*)buf->device;
    if (!dev) {
        free(buf);
        return;
    }

    if (buf->cpu_ptr) {
        vkUnmapMemory(dev->device, buf->memory);
        buf->cpu_ptr = NULL;
    }
    if (buf->buffer) {
        vkDestroyBuffer(dev->device, buf->buffer, NULL);
        buf->buffer = VK_NULL_HANDLE;
    }
    if (buf->memory) {
        vkFreeMemory(dev->device, buf->memory, NULL);
        buf->memory = VK_NULL_HANDLE;
    }
    free(buf);
}

void* md_gpu_buffer_cpu_ptr(md_gpu_buffer_t buffer) {
    if (!buffer) return NULL;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    if (!(buf->flags & MD_GPU_BUFFER_CPU_VISIBLE)) {
        MD_LOG_ERROR("md_gpu_buffer_cpu_ptr called on non-CPU-visible buffer");
        return NULL;
    }
    return buf->cpu_ptr;
}

uint64_t md_gpu_buffer_address(md_gpu_buffer_t buffer) {
    if (!buffer) return 0;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    md_gpu_device* dev = (md_gpu_device*)buf->device;

    VkBufferDeviceAddressInfo info = {0};
    info.sType = VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO;
    info.buffer = buf->buffer;
    return (uint64_t)vkGetBufferDeviceAddress(dev->device, &info);
}

md_gpu_buffer_flags_t md_gpu_buffer_flags(md_gpu_buffer_t buffer) {
    if (!buffer) return MD_GPU_BUFFER_NONE;
    return ((md_gpu_buffer*)buffer)->flags;
}

size_t md_gpu_buffer_size(md_gpu_buffer_t buffer) {
    if (!buffer) return 0;
    return (size_t)((md_gpu_buffer*)buffer)->size;
}

md_gpu_image_t md_gpu_image_create(md_gpu_device_t device, const md_gpu_image_desc_t* desc) {
    if (!device || !desc) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;

    md_gpu_image* img = (md_gpu_image*)calloc(1, sizeof(md_gpu_image));
    if (!img) return NULL;
    img->device = device;
    img->width = desc->width;
    img->height = desc->height;
    img->depth = desc->depth;
    img->format = md_vk_format_from_md(desc->format);
    img->flags = desc->flags;
    img->layout = VK_IMAGE_LAYOUT_UNDEFINED;

    VkImageCreateInfo ici = {0};
    ici.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    ici.imageType = (desc->depth > 1) ? VK_IMAGE_TYPE_3D : (desc->height > 1) ? VK_IMAGE_TYPE_2D : VK_IMAGE_TYPE_1D;
    ici.format = img->format;
    ici.extent.width = desc->width;
    ici.extent.height = desc->height;
    ici.extent.depth = desc->depth;
    ici.mipLevels = 1;
    ici.arrayLayers = 1;
    ici.samples = VK_SAMPLE_COUNT_1_BIT;
    ici.tiling = VK_IMAGE_TILING_OPTIMAL;
    ici.usage = 0;
    if (desc->flags & MD_GPU_IMAGE_STORAGE) {
        ici.usage |= VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    }
    if (desc->flags & MD_GPU_IMAGE_SAMPLED) {
        ici.usage |= VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    }
    if (desc->flags & MD_GPU_IMAGE_RENDER_TARGET) {
        ici.usage |= VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    }
    ici.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    ici.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    if (!md_vk_check(vkCreateImage(dev->device, &ici, NULL, &img->image), "vkCreateImage")) {
        free(img);
        return NULL;
    }

    VkMemoryRequirements req;
    vkGetImageMemoryRequirements(dev->device, img->image, &req);

    uint32_t mem_type = md_vk_find_memory_type(dev->physical_device, req.memoryTypeBits, 0, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    if (mem_type == UINT32_MAX) {
        MD_LOG_ERROR("Vulkan: failed to find memory type for image");
        vkDestroyImage(dev->device, img->image, NULL);
        free(img);
        return NULL;
    }

    VkMemoryAllocateInfo mai = {0};
    mai.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    mai.allocationSize = req.size;
    mai.memoryTypeIndex = mem_type;

    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &img->memory), "vkAllocateMemory(image)")) {
        vkDestroyImage(dev->device, img->image, NULL);
        free(img);
        return NULL;
    }

    if (!md_vk_check(vkBindImageMemory(dev->device, img->image, img->memory, 0), "vkBindImageMemory")) {
        vkFreeMemory(dev->device, img->memory, NULL);
        vkDestroyImage(dev->device, img->image, NULL);
        free(img);
        return NULL;
    }

    // Create image view
    VkImageViewCreateInfo ivci = {0};
    ivci.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    ivci.image = img->image;
    ivci.viewType = (desc->depth > 1) ? VK_IMAGE_VIEW_TYPE_3D : (desc->height > 1) ? VK_IMAGE_VIEW_TYPE_2D : VK_IMAGE_VIEW_TYPE_1D;
    ivci.format = img->format;
    ivci.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    ivci.subresourceRange.baseMipLevel = 0;
    ivci.subresourceRange.levelCount = 1;
    ivci.subresourceRange.baseArrayLayer = 0;
    ivci.subresourceRange.layerCount = 1;

    if (!md_vk_check(vkCreateImageView(dev->device, &ivci, NULL, &img->view), "vkCreateImageView")) {
        vkFreeMemory(dev->device, img->memory, NULL);
        vkDestroyImage(dev->device, img->image, NULL);
        free(img);
        return NULL;
    }

    return (md_gpu_image_t)img;
}

void md_gpu_image_destroy(md_gpu_image_t image) {
    if (!image) return;
    md_gpu_image* img = (md_gpu_image*)image;
    md_gpu_device* dev = (md_gpu_device*)img->device;
    if (!dev) {
        free(img);
        return;
    }

    if (img->view) {
        vkDestroyImageView(dev->device, img->view, NULL);
        img->view = VK_NULL_HANDLE;
    }
    if (img->image) {
        vkDestroyImage(dev->device, img->image, NULL);
        img->image = VK_NULL_HANDLE;
    }
    if (img->memory) {
        vkFreeMemory(dev->device, img->memory, NULL);
        img->memory = VK_NULL_HANDLE;
    }
    free(img);
}

bool md_gpu_image_desc_extract(md_gpu_image_t image, md_gpu_image_desc_t* out_desc) {
    if (!image || !out_desc) return false;
    md_gpu_image* img = (md_gpu_image*)image;

    out_desc->width = img->width;
    out_desc->height = img->height;
    out_desc->depth = img->depth;
    out_desc->format = md_vk_format_to_md(img->format);
    out_desc->flags = img->flags;

    return true;
}

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t device, const md_gpu_sampler_desc_t* desc) {
    if (!device || !desc) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;

    md_gpu_sampler* sampler = (md_gpu_sampler*)calloc(1, sizeof(md_gpu_sampler));
    if (!sampler) return NULL;
    sampler->device = device;

    VkSamplerCreateInfo sci = {0};
    sci.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
    sci.magFilter = md_vk_filter_from_md(desc->mag_filter);
    sci.minFilter = md_vk_filter_from_md(desc->min_filter);
    sci.mipmapMode = VK_SAMPLER_MIPMAP_MODE_NEAREST;
    sci.addressModeU = md_vk_address_mode_from_md(desc->address_u);
    sci.addressModeV = md_vk_address_mode_from_md(desc->address_v);
    sci.addressModeW = md_vk_address_mode_from_md(desc->address_w);
    sci.mipLodBias = 0.0f;
    sci.anisotropyEnable = VK_FALSE;
    sci.maxAnisotropy = 1.0f;
    sci.compareEnable = VK_FALSE;
    sci.compareOp = VK_COMPARE_OP_ALWAYS;
    sci.minLod = 0.0f;
    sci.maxLod = 0.0f;
    sci.borderColor = VK_BORDER_COLOR_FLOAT_TRANSPARENT_BLACK;
    sci.unnormalizedCoordinates = VK_FALSE;

    if (!md_vk_check(vkCreateSampler(dev->device, &sci, NULL, &sampler->sampler), "vkCreateSampler")) {
        free(sampler);
        return NULL;
    }

    return (md_gpu_sampler_t)sampler;
}

void md_gpu_sampler_destroy(md_gpu_sampler_t sampler) {
    if (!sampler) return;
    md_gpu_sampler* smp = (md_gpu_sampler*)sampler;
    md_gpu_device* dev = (md_gpu_device*)smp->device;
    if (dev && smp->sampler) {
        vkDestroySampler(dev->device, smp->sampler, NULL);
    }
    free(smp);
}

md_gpu_compute_pipeline_t md_gpu_compute_pipeline_create(md_gpu_device_t device, const md_gpu_compute_pipeline_desc_t* desc) {
    if (!device || !desc || !desc->shader_bytes || desc->shader_byte_size == 0) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;
    const char* entry_point = desc->entry_point ? desc->entry_point : "main";

    md_gpu_compute_pipeline* pipeline = (md_gpu_compute_pipeline*)calloc(1, sizeof(md_gpu_compute_pipeline));
    if (!pipeline) return NULL;
    pipeline->device = device;

    // Create shader module
    VkShaderModuleCreateInfo sm_ci = {0};
    sm_ci.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    sm_ci.codeSize = desc->shader_byte_size;
    sm_ci.pCode = (const uint32_t*)desc->shader_bytes;

    if (!md_vk_check(vkCreateShaderModule(dev->device, &sm_ci, NULL, &pipeline->shader_module), "vkCreateShaderModule")) {
        free(pipeline);
        return NULL;
    }

    VkPipelineShaderStageCreateInfo stage = {0};
    stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    stage.module = pipeline->shader_module;
    stage.pName = entry_point;

    VkComputePipelineCreateInfo cp_ci = {0};
    cp_ci.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    cp_ci.stage = stage;
    cp_ci.layout = dev->pipeline_layout;

    if (!md_vk_check(vkCreateComputePipelines(dev->device, VK_NULL_HANDLE, 1, &cp_ci, NULL, &pipeline->pipeline), "vkCreateComputePipelines")) {
        vkDestroyShaderModule(dev->device, pipeline->shader_module, NULL);
        free(pipeline);
        return NULL;
    }

    ASSERT(desc->resource_binding_count <= MD_GPU_MAX_RESOURCE_BINDINGS);
    if (desc->resource_bindings && desc->resource_binding_count > 0) {
        pipeline->resource_binding_count = desc->resource_binding_count;
        memcpy(pipeline->resource_bindings, desc->resource_bindings, sizeof(md_gpu_resource_binding_t) * desc->resource_binding_count);
    }

    if (desc->label && vkSetDebugUtilsObjectNameEXT) {
        VkDebugUtilsObjectNameInfoEXT name_info = {0};
        name_info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_OBJECT_NAME_INFO_EXT;
        name_info.pObjectName = desc->label;
        name_info.objectType = VK_OBJECT_TYPE_SHADER_MODULE;
        name_info.objectHandle = (uint64_t)pipeline->shader_module;
        vkSetDebugUtilsObjectNameEXT(dev->device, &name_info);
        name_info.objectType = VK_OBJECT_TYPE_PIPELINE;
        name_info.objectHandle = (uint64_t)pipeline->pipeline;
        vkSetDebugUtilsObjectNameEXT(dev->device, &name_info);
    }

    return (md_gpu_compute_pipeline_t)pipeline;
}

void md_gpu_compute_pipeline_destroy(md_gpu_compute_pipeline_t pipeline) {
    if (!pipeline) return;
    md_gpu_compute_pipeline* p = (md_gpu_compute_pipeline*)pipeline;
    md_gpu_device* dev = (md_gpu_device*)p->device;
    if (!dev) {
        free(p);
        return;
    }

    if (p->pipeline) {
        vkDestroyPipeline(dev->device, p->pipeline, NULL);
        p->pipeline = VK_NULL_HANDLE;
    }
    if (p->shader_module) {
        vkDestroyShaderModule(dev->device, p->shader_module, NULL);
        p->shader_module = VK_NULL_HANDLE;
    }
    free(p);
}

static uint64_t md_vk_queue_timeline_value(md_gpu_queue* queue) {
    if (!queue || !queue->dev || !queue->event_timeline) return UINT64_MAX;
    uint64_t value = 0;
    if (!md_vk_check(vkGetSemaphoreCounterValue(queue->dev->device, queue->event_timeline, &value), "vkGetSemaphoreCounterValue(event timeline)")) {
        return 0;
    }
    return value;
}

static bool md_vk_wait_queue_timeline(md_gpu_queue* queue, uint64_t value) {
    if (!queue || !queue->dev || !queue->event_timeline || value == 0) return true;

    VkSemaphoreWaitInfo wait_info = {0};
    wait_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO;
    wait_info.semaphoreCount = 1;
    wait_info.pSemaphores = &queue->event_timeline;
    wait_info.pValues = &value;
    return md_vk_check(vkWaitSemaphores(queue->dev->device, &wait_info, UINT64_MAX), "vkWaitSemaphores(event timeline)");
}

static void md_vk_reclaim_context_completed_cmds(md_vk_recording_context* ctx, uint64_t completed) {
    if (!ctx) return;
    uint64_t busy_slots = MD_GPU_VK_ALL_CMD_SLOTS & ~ctx->available_cmd_slots;
    while (busy_slots) {
        const uint32_t i = (uint32_t)ctz64(busy_slots);
        busy_slots &= busy_slots - 1;

        md_gpu_command_buffer* cmd = &ctx->cmd_wrappers[i];
        if (cmd->retire_value == 0 || cmd->retire_value > completed) continue;

        md_vk_release_cmd_slot(cmd);
    }
}

static void md_vk_reclaim_completed_cmds(md_gpu_queue* queue) {
    if (!queue || !queue->dev || !queue->event_timeline) return;
    const uint64_t completed = md_vk_queue_timeline_value(queue);
    md_vk_recording_context* head = atomic_load_explicit(&queue->recording_contexts, memory_order_acquire);
    for (md_vk_recording_context* ctx = head; ctx; ctx = ctx->next) {
        md_vk_reclaim_context_completed_cmds(ctx, completed);
    }
}

static md_gpu_command_buffer_t md_gpu_command_buffer_acquire(md_gpu_queue_t queue) {
    if (!queue) return NULL;
    md_gpu_queue* q = (md_gpu_queue*)queue;

    md_vk_recording_context* ctx = md_vk_queue_context_for_current_thread(q);
    if (!ctx) return NULL;

    if (ctx->available_cmd_slots == 0 && q->event_timeline) {
        md_vk_reclaim_context_completed_cmds(ctx, md_vk_queue_timeline_value(q));
    }

    if (ctx->available_cmd_slots == 0) {
        MD_LOG_ERROR("md_gpu_command_buffer_acquire: no free command buffers");
        return NULL;
    }

    const uint32_t idx = (uint32_t)ctz64(ctx->available_cmd_slots);
    ctx->available_cmd_slots &= ~md_vk_cmd_slot_mask(idx);

    md_gpu_command_buffer* cmd = &ctx->cmd_wrappers[idx];
    /* Reset transient recording state; queue, context, vk_cmd, and dispatch_sets pointer are permanent. */
    md_vk_reset_cmd_state(cmd);
    cmd->state = MD_VK_CMD_STATE_ACQUIRED;

    return (md_gpu_command_buffer_t)cmd;
}

static const md_gpu_resource_binding_t* md_vk_pipeline_find_resource_binding(const md_gpu_compute_pipeline* pipeline, md_gpu_resource_kind_t kind, uint32_t set, uint32_t binding) {
    if (!pipeline) return NULL;

    for (uint32_t i = 0; i < pipeline->resource_binding_count; ++i) {
        const md_gpu_resource_binding_t* shader_binding = &pipeline->resource_bindings[i];
        if (shader_binding->kind == kind && shader_binding->set == set && shader_binding->binding == binding) {
            return shader_binding;
        }
    }

    return NULL;
}

// Helper: flush declared dispatch resources into the per-category descriptor sets.
static bool md_vk_flush_descriptor_sets(md_gpu_command_buffer* cmd, md_gpu_device* dev, const md_gpu_compute_pipeline* pipeline, const md_vk_dispatch_resources_t* res, VkDescriptorSet sets[MD_GPU_DESC_SET_COUNT]) {
    VkWriteDescriptorSet writes[MD_GPU_MAX_RESOURCE_BINDINGS];
    VkDescriptorImageInfo descriptor_infos[MD_GPU_MAX_RESOURCE_BINDINGS];
    uint32_t write_count = 0;

    for (uint32_t i = 0; i < res->count; ++i) {
        md_gpu_cmd_resource_binding_t binding = res->bindings[i];
        const md_gpu_resource_binding_t* shader_binding = md_vk_pipeline_find_resource_binding(pipeline, binding.kind, binding.set, binding.binding);
        if (!shader_binding || shader_binding->set >= MD_GPU_DESC_SET_COUNT) return false;

        VkDescriptorImageInfo* info = &descriptor_infos[write_count];
        VkWriteDescriptorSet* w = &writes[write_count++];
        w->sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        w->pNext = NULL;
        w->dstSet = sets[shader_binding->set];
        w->dstBinding = shader_binding->backend_binding;
        w->dstArrayElement = 0;
        w->descriptorCount = 1;
        w->pBufferInfo = NULL;
        w->pImageInfo = info;
        w->pTexelBufferView = NULL;

        switch (binding.kind) {
            case MD_GPU_RESOURCE_STORAGE_IMAGE: {
                md_gpu_image* img = (md_gpu_image*)binding.image;
                info->sampler = VK_NULL_HANDLE;
                info->imageView = img->view;
                info->imageLayout = md_vk_cmd_image_layout(cmd, img);
                w->descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            } break;
            case MD_GPU_RESOURCE_SAMPLED_IMAGE: {
                md_gpu_image* img = (md_gpu_image*)binding.image;
                info->sampler = VK_NULL_HANDLE;
                info->imageView = img->view;
                info->imageLayout = md_vk_cmd_image_layout(cmd, img);
                w->descriptorType = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
            } break;
            case MD_GPU_RESOURCE_SAMPLER: {
                md_gpu_sampler* sampler = (md_gpu_sampler*)binding.sampler;
                info->sampler = sampler->sampler;
                info->imageView = VK_NULL_HANDLE;
                info->imageLayout = VK_IMAGE_LAYOUT_UNDEFINED;
                w->descriptorType = VK_DESCRIPTOR_TYPE_SAMPLER;
            } break;
            default:
                return false;
        }
    }

    if (write_count > 0) {
        vkUpdateDescriptorSets(dev->device, write_count, writes, 0, NULL);
    }

    return true;
}

static VkPipelineStageFlags2 md_vk_stage_flags(md_gpu_barrier_stage_t stages) {
    VkPipelineStageFlags2 flags = 0;
    if (stages & MD_GPU_BARRIER_STAGE_COMPUTE)  flags |= VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT;
    if (stages & MD_GPU_BARRIER_STAGE_TRANSFER) flags |= VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
    if (flags == 0) flags = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT; // safe fallback
    return flags;
}

static VkAccessFlags2 md_vk_access_flags_read(md_gpu_barrier_stage_t stages) {
    VkAccessFlags2 flags = 0;
    if (stages & MD_GPU_BARRIER_STAGE_COMPUTE)  flags |= VK_ACCESS_2_SHADER_READ_BIT;
    if (stages & MD_GPU_BARRIER_STAGE_TRANSFER) flags |= VK_ACCESS_2_TRANSFER_READ_BIT;
    if (flags == 0) flags = VK_ACCESS_2_MEMORY_READ_BIT; // safe fallback
    return flags;
}

static VkAccessFlags2 md_vk_access_flags_write(md_gpu_barrier_stage_t stages) {
    VkAccessFlags2 flags = 0;
    if (stages & MD_GPU_BARRIER_STAGE_COMPUTE)  flags |= VK_ACCESS_2_SHADER_WRITE_BIT;
    if (stages & MD_GPU_BARRIER_STAGE_TRANSFER) flags |= VK_ACCESS_2_TRANSFER_WRITE_BIT;
    if (flags == 0) flags = VK_ACCESS_2_MEMORY_WRITE_BIT; // safe fallback
    return flags;
}

static bool md_vk_cmd_is_recording(const md_gpu_command_buffer* cmd) {
    return cmd && cmd->state == MD_VK_CMD_STATE_RECORDING;
}

static bool md_vk_cmd_begin_acquired(md_gpu_cmd_t cmd_handle, const char* label) {
    if (!cmd_handle) return false;
    md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)cmd_handle;
    if (!cmd->queue || cmd->state != MD_VK_CMD_STATE_ACQUIRED) return false;

    VkCommandBufferBeginInfo begin_info = {0};
    begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    if (!md_vk_check(vkBeginCommandBuffer(cmd->vk_cmd, &begin_info), "vkBeginCommandBuffer")) {
        return false;
    }

    cmd->state = MD_VK_CMD_STATE_RECORDING;

    if (label) {
        VkDebugUtilsLabelEXT info = {0};
        info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT;
        info.pLabelName = label;
        if (vkCmdBeginDebugUtilsLabelEXT) vkCmdBeginDebugUtilsLabelEXT(cmd->vk_cmd, &info);
        cmd->pushed_debug_group = true;
    }

    return true;
}

md_gpu_cmd_t md_gpu_cmd_begin(md_gpu_queue_t queue_handle, const char* label) {
    md_gpu_cmd_t cmd = (md_gpu_cmd_t)md_gpu_command_buffer_acquire(queue_handle);
    if (!cmd) return NULL;
    if (!md_vk_cmd_begin_acquired(cmd, label)) {
        md_gpu_cmd_discard(cmd);
        return NULL;
    }
    return cmd;
}

bool md_gpu_cmd_end(md_gpu_cmd_t cmd_handle) {
    if (!cmd_handle) return false;
    md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)cmd_handle;
    if (cmd->state != MD_VK_CMD_STATE_RECORDING) return false;

    if (cmd->pushed_debug_group && vkCmdEndDebugUtilsLabelEXT) {
        vkCmdEndDebugUtilsLabelEXT(cmd->vk_cmd);
    }
    cmd->pushed_debug_group = false;

    if (!md_vk_check(vkEndCommandBuffer(cmd->vk_cmd), "vkEndCommandBuffer(cmd)")) {
        return false;
    }

    cmd->state = MD_VK_CMD_STATE_ENDED;
    return true;
}

void md_gpu_cmd_discard(md_gpu_cmd_t cmd_handle) {
    if (!cmd_handle) return;
    md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)cmd_handle;
    if (!cmd->queue) return;
    if (cmd->state == MD_VK_CMD_STATE_AVAILABLE || cmd->state == MD_VK_CMD_STATE_SUBMITTED) return;

    md_vk_release_cmd_slot(cmd);
}

md_gpu_event_t md_gpu_queue_submit(md_gpu_queue_t queue_handle, const md_gpu_queue_submit_desc_t* desc) {
    md_gpu_event_t event = {0};
    if (!queue_handle || !desc || !desc->cmds || desc->cmd_count == 0 || desc->cmd_count > UINT32_MAX || desc->wait_count > UINT32_MAX) return event;
    if (desc->wait_count > 0 && !desc->waits) return event;

    md_gpu_queue* queue = (md_gpu_queue*)queue_handle;
    const md_gpu_cmd_t* cmds = desc->cmds;
    const size_t cmd_count = desc->cmd_count;
    VkCommandBufferSubmitInfo* cmd_infos = (VkCommandBufferSubmitInfo*)malloc(sizeof(VkCommandBufferSubmitInfo) * cmd_count);
    if (!cmd_infos) return event;

    for (size_t i = 0; i < cmd_count; ++i) {
        md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)cmds[i];
        if (!cmd || cmd->queue != queue || cmd->state != MD_VK_CMD_STATE_ENDED) {
            free(cmd_infos);
            return event;
        }
        if (!md_vk_cmd_validate_initial_state(cmds, i, cmd)) {
            free(cmd_infos);
            return event;
        }

        VkCommandBufferSubmitInfo* info = &cmd_infos[i];
        info->sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO;
        info->pNext = NULL;
        info->commandBuffer = cmd->vk_cmd;
        info->deviceMask = 0;
    }

    md_vk_reclaim_completed_cmds(queue);

    VkSemaphoreSubmitInfo* wait_infos = NULL;
    if (desc->wait_count > 0) {
        wait_infos = (VkSemaphoreSubmitInfo*)malloc(sizeof(VkSemaphoreSubmitInfo) * desc->wait_count);
        if (!wait_infos) {
            free(cmd_infos);
            return event;
        }
    }

    uint32_t wait_count = 0;
    for (size_t i = 0; i < desc->wait_count; ++i) {
        const md_gpu_queue_wait_t* pending = &desc->waits[i];
        if (pending->event.value == 0) continue;

        md_gpu_queue* wait_queue = (md_gpu_queue*)pending->event.queue;
        if (!wait_queue) wait_queue = queue;
        if (!wait_queue->event_timeline) continue;

        VkSemaphoreSubmitInfo* wait_info = &wait_infos[wait_count++];
        wait_info->sType = VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO;
        wait_info->pNext = NULL;
        wait_info->semaphore = wait_queue->event_timeline;
        wait_info->value = pending->event.value;
        wait_info->stageMask = md_vk_stage_flags(pending->dst_stage);
        wait_info->deviceIndex = 0;
    }

    event.queue = (md_gpu_queue_t)queue;
    event.value = ++queue->next_event_id;
    if (event.value == 0) event.value = ++queue->next_event_id;

    VkSemaphoreSubmitInfo signal_info = {0};
    signal_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO;
    signal_info.semaphore = queue->event_timeline;
    signal_info.value = event.value;
    signal_info.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;

    VkSubmitInfo2 submit = {0};
    submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO_2;
    submit.waitSemaphoreInfoCount = wait_count;
    submit.pWaitSemaphoreInfos = wait_count ? wait_infos : NULL;
    submit.commandBufferInfoCount = (uint32_t)cmd_count;
    submit.pCommandBufferInfos = cmd_infos;
    submit.signalSemaphoreInfoCount = 1;
    submit.pSignalSemaphoreInfos = &signal_info;

    if (!md_vk_check(vkQueueSubmit2(queue->vk_queue, 1, &submit, VK_NULL_HANDLE), "vkQueueSubmit2(event timeline)")) {
        free(wait_infos);
        free(cmd_infos);
        event.queue = NULL;
        event.value = 0;
        return event;
    }

    for (size_t i = 0; i < cmd_count; ++i) {
        md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)cmds[i];
        md_vk_cmd_apply_final_state(cmd);
        cmd->retire_value = event.value;
        cmd->state = MD_VK_CMD_STATE_SUBMITTED;
    }

    free(wait_infos);
    free(cmd_infos);
    return event;
}

bool md_gpu_event_is_complete(md_gpu_event_t event) {
    if (event.value == 0) return true;
    md_gpu_queue* queue = (md_gpu_queue*)event.queue;
    if (!queue) return true;
    return md_vk_queue_timeline_value(queue) >= event.value;
}

void md_gpu_event_wait(md_gpu_event_t event) {
    if (event.value == 0) return;
    md_gpu_queue* queue = (md_gpu_queue*)event.queue;
    if (!queue) return;
    if (!md_vk_wait_queue_timeline(queue, event.value)) {
        MD_LOG_ERROR("md_gpu_event_wait: failed to wait for event timeline");
    }
}

static bool md_vk_cmd_declare_resources(md_gpu_command_buffer* cmd, md_vk_dispatch_resources_t* res, const md_gpu_resource_t* resources, uint32_t resource_count) {
    if (!cmd || !res) return false;

    if (!resources || resource_count == 0) return true;

    for (uint32_t i = 0; i < resource_count; ++i) {
        const md_gpu_resource_t* resource = &resources[i];
        switch (resource->kind) {
            case MD_GPU_RESOURCE_BUFFER_USAGE: {
                if (!resource->buffer || resource->usage == 0) return false;
                md_vk_ensure_buffer_sync(cmd, (md_gpu_buffer*)resource->buffer, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, md_vk_access_from_gpu_usage(resource->usage));
            } break;
            case MD_GPU_RESOURCE_STORAGE_IMAGE: {
                const md_gpu_usage_flags_t usage = resource->usage ? resource->usage : MD_GPU_USAGE_READ;
                if (!md_vk_cmd_merge_resource_binding(res, resource)) return false;
                md_vk_ensure_image_sync(cmd, (md_gpu_image*)resource->image, VK_IMAGE_LAYOUT_GENERAL, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, md_vk_access_from_gpu_usage(usage));
            } break;
            case MD_GPU_RESOURCE_SAMPLED_IMAGE: {
                const md_gpu_usage_flags_t usage = resource->usage ? resource->usage : MD_GPU_USAGE_READ;
                if (!md_vk_cmd_merge_resource_binding(res, resource)) return false;
                md_vk_ensure_image_sync(cmd, (md_gpu_image*)resource->image, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL, VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, md_vk_access_from_gpu_usage(usage));
            } break;
            case MD_GPU_RESOURCE_SAMPLER:
                if (!md_vk_cmd_merge_resource_binding(res, resource)) return false;
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

    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    md_gpu_queue* q = c->queue;
    md_gpu_compute_pipeline* pipeline = (md_gpu_compute_pipeline*)dispatch->pipeline;

    if (!md_vk_cmd_is_recording(c)) return false;
    ASSERT(c->state == MD_VK_CMD_STATE_RECORDING);
    ASSERT(dispatch->pipeline && "Pipeline must be supplied for dispatch");

    if (dispatch->root_args_size > MD_GPU_MAX_PUSH_CONSTANTS) {
        ASSERT(false && "Root argument payload exceeds push constant limit");
        return false;
    }

    const bool pipeline_changed = (c->bound_pipeline != dispatch->pipeline);
    c->bound_pipeline = dispatch->pipeline;

    md_vk_dispatch_resources_t res = {0};
    if (!md_vk_cmd_declare_resources(c, &res, dispatch->resources, dispatch->resource_count)) return false;

    if (pipeline_changed) {
        vkCmdBindPipeline(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline->pipeline);
    }

    if (dispatch->root_args_size > 0) {
        vkCmdPushConstants(c->vk_cmd, q->dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, (uint32_t)dispatch->root_args_size, dispatch->root_args);
    }

    ASSERT(c->dispatch_count < MD_GPU_VK_MAX_DISPATCHES_PER_CMD && "Too many dispatches per command buffer");
    if (c->dispatch_count >= MD_GPU_VK_MAX_DISPATCHES_PER_CMD) {
        return false;
    }
    const uint32_t dispatch_slot = c->dispatch_count++;
    VkDescriptorSet* descriptor_sets = c->dispatch_sets[dispatch_slot];
    if (!md_vk_flush_descriptor_sets(c, q->dev, pipeline, &res, descriptor_sets)) return false;

    vkCmdBindDescriptorSets(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, q->dev->pipeline_layout, 0, MD_GPU_DESC_SET_COUNT, descriptor_sets, 0, NULL);

    vkCmdDispatch(c->vk_cmd, dispatch->group_count[0], dispatch->group_count[1], dispatch->group_count[2]);

    return true;
}

bool md_gpu_cmd_barrier(md_gpu_cmd_t cmd_handle, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!cmd_handle) return false;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return false;

    VkMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER_2;
    barrier.srcStageMask = md_vk_stage_flags(src_stage);
    barrier.srcAccessMask = md_vk_access_flags_write(src_stage);
    barrier.dstStageMask = md_vk_stage_flags(dst_stage);
    barrier.dstAccessMask = md_vk_access_flags_read(dst_stage);

    VkDependencyInfo dep = {0};
    dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    dep.memoryBarrierCount = 1;
    dep.pMemoryBarriers = &barrier;
    vkCmdPipelineBarrier2(c->vk_cmd, &dep);
    return true;
}

void md_gpu_cmd_push_debug_group(md_gpu_cmd_t cmd_handle, const char* label) {
    if (!cmd_handle || !label) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return;
    VkDebugUtilsLabelEXT info = {0};
    info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT;
    info.pLabelName = label;
    if (vkCmdBeginDebugUtilsLabelEXT) vkCmdBeginDebugUtilsLabelEXT(c->vk_cmd, &info);
}

void md_gpu_cmd_pop_debug_group(md_gpu_cmd_t cmd_handle) {
    if (!cmd_handle) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return;
    if (vkCmdEndDebugUtilsLabelEXT) vkCmdEndDebugUtilsLabelEXT(c->vk_cmd);
}

static bool md_vk_buffer_image_copy_region(VkBufferImageCopy* region, const md_gpu_image* image, const md_gpu_buffer_image_copy_t* copy) {
    if (!region || !image) return false;

    const uint32_t bpp = md_vk_format_bytes_per_pixel(image->format);
    if (bpp == 0) return false;

    const md_gpu_image_region_t image_region = copy ? copy->image_region : (md_gpu_image_region_t){0};
    const uint32_t ox = image_region.offset[0];
    const uint32_t oy = image_region.offset[1];
    const uint32_t oz = image_region.offset[2];
    if (ox >= image->width || oy >= image->height || oz >= image->depth) return false;

    const uint32_t ew = image_region.extent[0] ? image_region.extent[0] : image->width  - ox;
    const uint32_t eh = image_region.extent[1] ? image_region.extent[1] : image->height - oy;
    const uint32_t ed = image_region.extent[2] ? image_region.extent[2] : image->depth  - oz;
    if (ew == 0 || eh == 0 || ed == 0) return false;
    if (ew > image->width - ox || eh > image->height - oy || ed > image->depth - oz) return false;

    const size_t tight_row = (size_t)ew * (size_t)bpp;
    const size_t bytes_per_row = (copy && copy->bytes_per_row) ? copy->bytes_per_row : tight_row;
    if (bytes_per_row < tight_row || (bytes_per_row % bpp) != 0) return false;

    uint32_t row_length = 0;
    if (bytes_per_row != tight_row) {
        const size_t texel_row_length = bytes_per_row / bpp;
        if (texel_row_length > UINT32_MAX) return false;
        row_length = (uint32_t)texel_row_length;
    }

    uint32_t image_height = 0;
    if (copy && copy->bytes_per_image) {
        if (copy->bytes_per_image < bytes_per_row * (size_t)eh || (copy->bytes_per_image % bytes_per_row) != 0) return false;
        const size_t texel_image_height = copy->bytes_per_image / bytes_per_row;
        if (texel_image_height != eh) {
            if (texel_image_height > UINT32_MAX) return false;
            image_height = (uint32_t)texel_image_height;
        }
    }

    region->bufferOffset = copy ? copy->buffer_offset : 0;
    region->bufferRowLength = row_length;
    region->bufferImageHeight = image_height;
    region->imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region->imageSubresource.mipLevel = 0;
    region->imageSubresource.baseArrayLayer = 0;
    region->imageSubresource.layerCount = 1;
    region->imageOffset.x = (int32_t)ox;
    region->imageOffset.y = (int32_t)oy;
    region->imageOffset.z = (int32_t)oz;
    region->imageExtent.width = ew;
    region->imageExtent.height = eh;
    region->imageExtent.depth = ed;
    return true;
}

bool md_gpu_cmd_copy_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!cmd_handle || !src || !dst || size == 0) return false;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return false;
    md_gpu_buffer* src_buf = (md_gpu_buffer*)src;
    md_gpu_buffer* dst_buf = (md_gpu_buffer*)dst;

    md_vk_ensure_buffer_sync(c, src_buf, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
    md_vk_ensure_buffer_sync(c, dst_buf, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    VkBufferCopy region = {0};
    region.srcOffset = src_offset;
    region.dstOffset = dst_offset;
    region.size = size;
    vkCmdCopyBuffer(c->vk_cmd, src_buf->buffer, dst_buf->buffer, 1, &region);
    return true;
}

bool md_gpu_cmd_copy_image_to_buffer(md_gpu_cmd_t cmd_handle, md_gpu_image_t src_image, md_gpu_buffer_t dst_buffer, const md_gpu_buffer_image_copy_t* copy) {
    if (!cmd_handle || !src_image || !dst_buffer) return false;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return false;
    md_gpu_image* src = (md_gpu_image*)src_image;
    md_gpu_buffer* dst = (md_gpu_buffer*)dst_buffer;

    VkBufferImageCopy region = {0};
    if (!md_vk_buffer_image_copy_region(&region, src, copy)) return false;

    md_vk_ensure_image_sync(c, src, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
    md_vk_ensure_buffer_sync(c, dst, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    vkCmdCopyImageToBuffer(c->vk_cmd, src->image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, dst->buffer, 1, &region);
    return true;
}

bool md_gpu_cmd_copy_buffer_to_image_layout(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image, const md_gpu_buffer_image_copy_t* copy) {
    if (!cmd_handle || !src_buffer || !dst_image) return false;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return false;
    md_gpu_buffer* src = (md_gpu_buffer*)src_buffer;
    md_gpu_image* dst = (md_gpu_image*)dst_image;

    VkBufferImageCopy region = {0};
    if (!md_vk_buffer_image_copy_region(&region, dst, copy)) return false;

    md_vk_ensure_buffer_sync(c, src, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);
    md_vk_ensure_image_sync(c, dst, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    vkCmdCopyBufferToImage(c->vk_cmd, src->buffer, dst->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);
    return true;
}

bool md_gpu_cmd_fill_buffer(md_gpu_cmd_t cmd_handle, md_gpu_buffer_t buffer, size_t offset, size_t size, uint8_t value) {
    if (!cmd_handle || !buffer) return false;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd_handle;
    if (!md_vk_cmd_is_recording(c)) return false;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;

    ASSERT((offset % 4 == 0) && "vkCmdFillBuffer offset must be 4-byte aligned");
    ASSERT((size == VK_WHOLE_SIZE || size % 4 == 0) && "vkCmdFillBuffer size must be 4-byte aligned or VK_WHOLE_SIZE");
    if ((offset % 4) != 0 || (size != VK_WHOLE_SIZE && (size % 4) != 0)) return false;

    md_vk_ensure_buffer_sync(c, buf, VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    uint32_t pattern = (value << 24) | (value << 16) | (value << 8) | value;
    vkCmdFillBuffer(c->vk_cmd, buf->buffer, offset, size, pattern);
    return true;
}

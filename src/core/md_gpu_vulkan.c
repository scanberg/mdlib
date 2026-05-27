#include "md_gpu.h"

#include <core/md_common.h>
#include <core/md_log.h>
#include <stdlib.h>
#include <string.h>

// We avoid including platform-specific Vulkan headers via GLFW here.
// This is a compute-only backend for now, so no surface extensions are required.
#define VOLK_IMPLEMENTATION
#include <volk.h>

// Vulkan backend (C).
// Compute-only for now; designed to be extended with graphics later.

// Command buffer pool size (not to be confused with Metal's MD_GPU_MAX_COMMANDS,
// which is the max number of recorded commands per command buffer)
#define MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE 64
// Maximum number of dispatch calls recorded per command buffer.
// One fresh descriptor set is allocated per dispatch, so the descriptor pool is
// sized as POOL_SIZE * MAX_DISPATCHES_PER_CMD.
#define MD_GPU_VK_MAX_DISPATCHES_PER_CMD 32
#define MD_GPU_PASS_HISTORY_SIZE 256

typedef struct md_gpu_queue* md_gpu_queue_t;
typedef struct md_gpu_command_buffer* md_gpu_command_buffer_t;

typedef struct md_gpu_pass_buffer_usage_t {
    md_gpu_buffer_t buffer;
    gpu_usage_flags_t usage;
} md_gpu_pass_buffer_usage_t;

typedef struct md_gpu_pass_image_usage_t {
    md_gpu_image_t image;
    gpu_usage_flags_t usage;
    md_gpu_image_region_t region;
} md_gpu_pass_image_usage_t;

typedef struct md_gpu_pass_sampled_resource_usage_t {
    md_gpu_image_t image;
    md_gpu_sampler_t sampler;
    gpu_usage_flags_t usage;
} md_gpu_pass_sampled_resource_usage_t;

typedef struct md_gpu_pass_record {
    uint64_t id;
    struct md_gpu_queue* queue;
    struct md_gpu_command_buffer* cmd;
} md_gpu_pass_record;

typedef struct md_gpu_queue {
    struct md_gpu_device* dev;
    VkQueue               vk_queue;
    uint32_t              family_index;
    VkCommandPool         command_pool;
    // Command buffer pool - wrappers allocated once per queue, recycled by index
    struct md_gpu_command_buffer* cmd_wrappers;
    uint32_t              free_list[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    uint32_t              free_count;
} md_gpu_queue;

typedef struct md_gpu_device {
    VkInstance instance;
    VkPhysicalDevice physical_device;
    VkDevice device;

    uint32_t compute_queue_family;

    VkDescriptorPool descriptor_pool;
    VkDescriptorSetLayout descriptor_set_layout;
    VkPipelineLayout pipeline_layout;
    VkSemaphore pass_timeline;

    // Embedded primary compute queue used by pass submission.
    struct md_gpu_queue compute_queue;

    uint64_t next_pass_id;
    uint32_t pass_record_head;
    md_gpu_pass_record pass_records[MD_GPU_PASS_HISTORY_SIZE];
} md_gpu_device;

static void md_vk_reclaim_completed_passes(md_gpu_device* dev);

typedef struct md_gpu_buffer {
    md_gpu_device_t device;
    VkBuffer buffer;
    VkDeviceMemory memory;
    VkDeviceSize size;
    md_gpu_buffer_flags_t flags;
    void* cpu_ptr;  /* non-NULL for CPU_VISIBLE buffers; valid until destroy */
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
    VkImageLayout layout; // current layout for implicit transitions
} md_gpu_image;

typedef struct md_gpu_sampler {
    md_gpu_device_t device;
    VkSampler sampler;
} md_gpu_sampler;

typedef struct md_gpu_compute_pipeline {
    md_gpu_device_t device;
    VkPipeline pipeline;
    VkShaderModule shader_module;
} md_gpu_compute_pipeline;

typedef struct md_gpu_command_buffer {
    struct md_gpu_device* dev;
    struct md_gpu_queue*  queue;
    VkCommandBuffer vk_cmd;
    /* Per-dispatch descriptor sets, allocated at dispatch-record time and freed on recycle */
    VkDescriptorSet dispatch_sets[MD_GPU_VK_MAX_DISPATCHES_PER_CMD];
    uint32_t dispatch_count;
    uint32_t index; /* position in queue->cmd_wrappers[], used to return slot to free list */
    bool is_recording;

    // Live bound state consumed by dispatch
    md_gpu_compute_pipeline_t bound_pipeline;
    md_gpu_buffer_t bound_buffers[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_offsets[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_sizes[MD_GPU_MAX_BIND_SLOTS];
    md_gpu_image_t bound_images[MD_GPU_MAX_BIND_SLOTS];
    md_gpu_image_t bound_sampled_images[MD_GPU_MAX_BIND_SLOTS];
    md_gpu_sampler_t bound_samplers[MD_GPU_MAX_BIND_SLOTS];
    uint8_t push_constants[MD_GPU_MAX_PUSH_CONSTANTS];
    size_t push_constant_size;

    md_gpu_pass_buffer_usage_t* pass_buffers;
    uint32_t pass_buffer_count;
    md_gpu_pass_image_usage_t* pass_images;
    uint32_t pass_image_count;
    md_gpu_pass_sampled_resource_usage_t* pass_sampled_resources;
    uint32_t pass_sampled_resource_count;
} md_gpu_command_buffer;

typedef struct md_gpu_pass {
    struct md_gpu_queue* queue;
    struct md_gpu_command_buffer* cmd;
    bool pushed_debug_group;
} md_gpu_pass;

static bool md_vk_check(VkResult res, const char* what) {
    if (res == VK_SUCCESS) return true;
    MD_LOG_ERROR("Vulkan error %d: %s", (int)res, what ? what : "(unknown)");
    return false;
}

static bool md_vk_pass_merge_buffer_usage(md_gpu_command_buffer* cmd, md_gpu_buffer_t buffer, gpu_usage_flags_t usage) {
    if (!cmd || !buffer || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->pass_buffer_count; ++i) {
        if (cmd->pass_buffers[i].buffer == buffer) {
            cmd->pass_buffers[i].usage |= usage;
            return true;
        }
    }

    const uint32_t new_count = cmd->pass_buffer_count + 1;
    md_gpu_pass_buffer_usage_t* items = (md_gpu_pass_buffer_usage_t*)realloc(cmd->pass_buffers, sizeof(md_gpu_pass_buffer_usage_t) * new_count);
    if (!items) return false;

    cmd->pass_buffers = items;
    cmd->pass_buffers[cmd->pass_buffer_count].buffer = buffer;
    cmd->pass_buffers[cmd->pass_buffer_count].usage = usage;
    cmd->pass_buffer_count = new_count;
    return true;
}

static bool md_vk_pass_merge_image_usage(md_gpu_command_buffer* cmd, md_gpu_image_t image, gpu_usage_flags_t usage) {
    if (!cmd || !image || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->pass_image_count; ++i) {
        if (cmd->pass_images[i].image == image) {
            cmd->pass_images[i].usage |= usage;
            return true;
        }
    }

    const uint32_t new_count = cmd->pass_image_count + 1;
    md_gpu_pass_image_usage_t* items = (md_gpu_pass_image_usage_t*)realloc(cmd->pass_images, sizeof(md_gpu_pass_image_usage_t) * new_count);
    if (!items) return false;

    cmd->pass_images = items;
    cmd->pass_images[cmd->pass_image_count].image = image;
    cmd->pass_images[cmd->pass_image_count].usage = usage;
    memset(&cmd->pass_images[cmd->pass_image_count].region, 0, sizeof(cmd->pass_images[cmd->pass_image_count].region));
    cmd->pass_image_count = new_count;
    return true;
}

static bool md_vk_pass_merge_sampled_resource_usage(md_gpu_command_buffer* cmd, md_gpu_image_t image, md_gpu_sampler_t sampler, gpu_usage_flags_t usage) {
    if (!cmd || !image || usage == 0) return true;

    for (uint32_t i = 0; i < cmd->pass_sampled_resource_count; ++i) {
        if (cmd->pass_sampled_resources[i].image == image && cmd->pass_sampled_resources[i].sampler == sampler) {
            cmd->pass_sampled_resources[i].usage |= usage;
            return true;
        }
    }

    const uint32_t new_count = cmd->pass_sampled_resource_count + 1;
    md_gpu_pass_sampled_resource_usage_t* items = (md_gpu_pass_sampled_resource_usage_t*)realloc(cmd->pass_sampled_resources, sizeof(md_gpu_pass_sampled_resource_usage_t) * new_count);
    if (!items) return false;

    cmd->pass_sampled_resources = items;
    cmd->pass_sampled_resources[cmd->pass_sampled_resource_count].image = image;
    cmd->pass_sampled_resources[cmd->pass_sampled_resource_count].sampler = sampler;
    cmd->pass_sampled_resources[cmd->pass_sampled_resource_count].usage = usage;
    cmd->pass_sampled_resource_count = new_count;
    return true;
}

static void md_vk_reset_cmd_state(md_gpu_command_buffer* cmd, md_gpu_queue* queue) {
    free(cmd->pass_buffers);
    free(cmd->pass_images);
    free(cmd->pass_sampled_resources);
    cmd->pass_buffers = NULL;
    cmd->pass_images = NULL;
    cmd->pass_sampled_resources = NULL;
    cmd->pass_buffer_count = 0;
    cmd->pass_image_count = 0;
    cmd->pass_sampled_resource_count = 0;

    cmd->queue = queue;
    cmd->dispatch_count = 0;
    cmd->is_recording = false;
    cmd->bound_pipeline = NULL;
    memset(cmd->bound_buffers, 0, sizeof(cmd->bound_buffers));
    memset(cmd->bound_buffer_offsets, 0, sizeof(cmd->bound_buffer_offsets));
    memset(cmd->bound_buffer_sizes, 0, sizeof(cmd->bound_buffer_sizes));
    memset(cmd->bound_images, 0, sizeof(cmd->bound_images));
    memset(cmd->bound_sampled_images, 0, sizeof(cmd->bound_sampled_images));
    memset(cmd->bound_samplers, 0, sizeof(cmd->bound_samplers));
    memset(cmd->push_constants, 0, sizeof(cmd->push_constants));
    cmd->push_constant_size = 0;
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

static bool md_vk_create_device(md_gpu_device* out_dev) {
    // Require sync2 for sane barrier implementation.
    // Use VK_KHR_synchronization2 extension if Vulkan < 1.3 (but we request 1.3 above).
    uint32_t ext_count = 0;
    if (!md_vk_check(vkEnumerateDeviceExtensionProperties(out_dev->physical_device, NULL, &ext_count, NULL), "vkEnumerateDeviceExtensionProperties(count)")) return false;
    VkExtensionProperties* exts = NULL;
    if (ext_count) {
        exts = (VkExtensionProperties*)malloc(sizeof(VkExtensionProperties) * ext_count);
        if (!exts) return false;
        if (!md_vk_check(vkEnumerateDeviceExtensionProperties(out_dev->physical_device, NULL, &ext_count, exts), "vkEnumerateDeviceExtensionProperties(list)")) {
            free(exts);
            return false;
        }
    }

    const bool has_sync2 = (exts && md_vk_has_extension(ext_count, exts, VK_KHR_SYNCHRONIZATION_2_EXTENSION_NAME));
    free(exts);

    float prio = 1.0f;
    VkDeviceQueueCreateInfo qci = {0};
    qci.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    qci.queueFamilyIndex = out_dev->compute_queue_family;
    qci.queueCount = 1;
    qci.pQueuePriorities = &prio;

    const char* dev_exts[8];
    uint32_t dev_ext_count = 0;
    if (has_sync2) {
        dev_exts[dev_ext_count++] = VK_KHR_SYNCHRONIZATION_2_EXTENSION_NAME;
    }

    VkPhysicalDeviceVulkan12Features f12 = {0};
    f12.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES;
    f12.bufferDeviceAddress = VK_TRUE;
    f12.timelineSemaphore = VK_TRUE;

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
    dci.enabledExtensionCount = dev_ext_count;
    dci.ppEnabledExtensionNames = (dev_ext_count ? dev_exts : NULL);
    dci.pNext = &feats2;

    if (!md_vk_check(vkCreateDevice(out_dev->physical_device, &dci, NULL, &out_dev->device), "vkCreateDevice")) return false;

    md_gpu_queue* q = &out_dev->compute_queue;
    vkGetDeviceQueue(out_dev->device, out_dev->compute_queue_family, 0, &q->vk_queue);
    q->dev          = out_dev;
    q->family_index = out_dev->compute_queue_family;
    if (q->vk_queue == VK_NULL_HANDLE) return false;

    // Create command pool owned by the compute queue
    VkCommandPoolCreateInfo pool_ci = {0};
    pool_ci.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    pool_ci.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pool_ci.queueFamilyIndex = out_dev->compute_queue_family;
    if (!md_vk_check(vkCreateCommandPool(out_dev->device, &pool_ci, NULL, &q->command_pool), "vkCreateCommandPool")) return false;

    // Allocate command buffer wrappers owned by the queue
    q->cmd_wrappers = (struct md_gpu_command_buffer*)calloc(MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE, sizeof(struct md_gpu_command_buffer));
    if (!q->cmd_wrappers) return false;

    VkCommandBuffer vk_cmds[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    VkCommandBufferAllocateInfo cb_ai = {0};
    cb_ai.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    cb_ai.commandPool = q->command_pool;
    cb_ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cb_ai.commandBufferCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;
    if (!md_vk_check(vkAllocateCommandBuffers(out_dev->device, &cb_ai, vk_cmds), "vkAllocateCommandBuffers")) return false;

    // Wire each wrapper to its VkCommandBuffer slot
    for (uint32_t i = 0; i < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE; ++i) {
        q->cmd_wrappers[i].dev   = out_dev;
        q->cmd_wrappers[i].queue = q;
        q->cmd_wrappers[i].vk_cmd = vk_cmds[i];
        q->cmd_wrappers[i].index  = i;
        q->free_list[i] = i;
    }
    q->free_count = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;

    // Create descriptor set layout: buffers, storage images, sampled images, and samplers occupy separate Vulkan binding ranges.
    VkDescriptorSetLayoutBinding bindings[4 * MD_GPU_MAX_BIND_SLOTS];
    for (uint32_t i = 0; i < MD_GPU_MAX_BIND_SLOTS; ++i) {
        bindings[i].binding = i;
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[i].pImmutableSamplers = NULL;

        bindings[MD_GPU_STORAGE_IMAGE_BINDING_BASE + i].binding = MD_GPU_STORAGE_IMAGE_BINDING_BASE + i;
        bindings[MD_GPU_STORAGE_IMAGE_BINDING_BASE + i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        bindings[MD_GPU_STORAGE_IMAGE_BINDING_BASE + i].descriptorCount = 1;
        bindings[MD_GPU_STORAGE_IMAGE_BINDING_BASE + i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[MD_GPU_STORAGE_IMAGE_BINDING_BASE + i].pImmutableSamplers = NULL;

        bindings[MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i].binding = MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i;
        bindings[MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i].descriptorType = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
        bindings[MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i].descriptorCount = 1;
        bindings[MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i].pImmutableSamplers = NULL;

        bindings[MD_GPU_SAMPLER_BINDING_BASE + i].binding = MD_GPU_SAMPLER_BINDING_BASE + i;
        bindings[MD_GPU_SAMPLER_BINDING_BASE + i].descriptorType = VK_DESCRIPTOR_TYPE_SAMPLER;
        bindings[MD_GPU_SAMPLER_BINDING_BASE + i].descriptorCount = 1;
        bindings[MD_GPU_SAMPLER_BINDING_BASE + i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[MD_GPU_SAMPLER_BINDING_BASE + i].pImmutableSamplers = NULL;
    }

    VkDescriptorSetLayoutCreateInfo dsl_ci = {0};
    dsl_ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    dsl_ci.bindingCount = 4 * MD_GPU_MAX_BIND_SLOTS;
    dsl_ci.pBindings = bindings;
    if (!md_vk_check(vkCreateDescriptorSetLayout(out_dev->device, &dsl_ci, NULL, &out_dev->descriptor_set_layout), "vkCreateDescriptorSetLayout")) return false;

    // Create pipeline layout with push constants
    VkPushConstantRange push_range = {0};
    push_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    push_range.offset = 0;
    push_range.size = MD_GPU_MAX_PUSH_CONSTANTS;

    VkPipelineLayoutCreateInfo pl_ci = {0};
    pl_ci.sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    pl_ci.setLayoutCount = 1;
    pl_ci.pSetLayouts = &out_dev->descriptor_set_layout;
    pl_ci.pushConstantRangeCount = 1;
    pl_ci.pPushConstantRanges = &push_range;
    if (!md_vk_check(vkCreatePipelineLayout(out_dev->device, &pl_ci, NULL, &out_dev->pipeline_layout), "vkCreatePipelineLayout")) return false;

    // Create descriptor pool
    VkDescriptorPoolSize pool_sizes[4];
    pool_sizes[0].type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    pool_sizes[0].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[1].type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    pool_sizes[1].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[2].type = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
    pool_sizes[2].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[3].type = VK_DESCRIPTOR_TYPE_SAMPLER;
    pool_sizes[3].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD * MD_GPU_MAX_BIND_SLOTS;

    VkDescriptorPoolCreateInfo dp_ci = {0};
    dp_ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    dp_ci.flags = 0; /* no FREE_DESCRIPTOR_SET_BIT; sets are pre-allocated and reused */
    dp_ci.maxSets = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD;
    dp_ci.poolSizeCount = 4;
    dp_ci.pPoolSizes = pool_sizes;
    if (!md_vk_check(vkCreateDescriptorPool(out_dev->device, &dp_ci, NULL, &out_dev->descriptor_pool), "vkCreateDescriptorPool")) return false;

    VkSemaphoreTypeCreateInfo timeline_ci = {0};
    timeline_ci.sType = VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO;
    timeline_ci.semaphoreType = VK_SEMAPHORE_TYPE_TIMELINE;
    timeline_ci.initialValue = 0;

    VkSemaphoreCreateInfo sem_ci = {0};
    sem_ci.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
    sem_ci.pNext = &timeline_ci;
    if (!md_vk_check(vkCreateSemaphore(out_dev->device, &sem_ci, NULL, &out_dev->pass_timeline), "vkCreateSemaphore(pass timeline)")) return false;

    // Pre-allocate all descriptor sets for every command buffer slot so they can be
    // reused without per-submission alloc/free (vkUpdateDescriptorSets overwrites in-place).
    {
        const uint32_t total_sets = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_VK_MAX_DISPATCHES_PER_CMD;
        VkDescriptorSetLayout* set_layouts = (VkDescriptorSetLayout*)malloc(total_sets * sizeof(VkDescriptorSetLayout));
        VkDescriptorSet*       all_sets    = (VkDescriptorSet*)      malloc(total_sets * sizeof(VkDescriptorSet));
        if (!set_layouts || !all_sets) { free(set_layouts); free(all_sets); return false; }
        for (uint32_t i = 0; i < total_sets; ++i)
            set_layouts[i] = out_dev->descriptor_set_layout;
        VkDescriptorSetAllocateInfo all_ds_ai = {0};
        all_ds_ai.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
        all_ds_ai.descriptorPool = out_dev->descriptor_pool;
        all_ds_ai.descriptorSetCount = total_sets;
        all_ds_ai.pSetLayouts = set_layouts;
        bool ok = md_vk_check(vkAllocateDescriptorSets(out_dev->device, &all_ds_ai, all_sets), "vkAllocateDescriptorSets (pre-alloc)");
        if (ok) {
            for (uint32_t i = 0; i < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE; ++i) {
                for (uint32_t j = 0; j < MD_GPU_VK_MAX_DISPATCHES_PER_CMD; ++j)
                    q->cmd_wrappers[i].dispatch_sets[j] = all_sets[i * MD_GPU_VK_MAX_DISPATCHES_PER_CMD + j];
            }
        }
        free(set_layouts);
        free(all_sets);
        if (!ok) return false;
    }

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
        for (uint32_t i = 0; i < MD_GPU_PASS_HISTORY_SIZE; ++i) {
            dev->pass_records[i].id = 0;
            dev->pass_records[i].queue = NULL;
            dev->pass_records[i].cmd = NULL;
        }
        if (dev->pass_timeline) vkDestroySemaphore(dev->device, dev->pass_timeline, NULL);
        if (dev->descriptor_pool) vkDestroyDescriptorPool(dev->device, dev->descriptor_pool, NULL);
        if (dev->pipeline_layout) vkDestroyPipelineLayout(dev->device, dev->pipeline_layout, NULL);
        if (dev->descriptor_set_layout) vkDestroyDescriptorSetLayout(dev->device, dev->descriptor_set_layout, NULL);
        if (dev->compute_queue.command_pool) vkDestroyCommandPool(dev->device, dev->compute_queue.command_pool, NULL);
        vkDestroyDevice(dev->device, NULL);
    }
    if (dev->instance) {
        vkDestroyInstance(dev->instance, NULL);
    }
    free(dev->compute_queue.cmd_wrappers);
    free(dev);
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
    stage.pName = "main";

    VkComputePipelineCreateInfo cp_ci = {0};
    cp_ci.sType = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    cp_ci.stage = stage;
    cp_ci.layout = dev->pipeline_layout;

    if (!md_vk_check(vkCreateComputePipelines(dev->device, VK_NULL_HANDLE, 1, &cp_ci, NULL, &pipeline->pipeline), "vkCreateComputePipelines")) {
        vkDestroyShaderModule(dev->device, pipeline->shader_module, NULL);
        free(pipeline);
        return NULL;
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

static uint64_t md_vk_pass_timeline_value(md_gpu_device* dev) {
    if (!dev || !dev->pass_timeline) return UINT64_MAX;
    uint64_t value = 0;
    if (!md_vk_check(vkGetSemaphoreCounterValue(dev->device, dev->pass_timeline, &value), "vkGetSemaphoreCounterValue(pass timeline)")) {
        return 0;
    }
    return value;
}

static bool md_vk_wait_pass_timeline(md_gpu_device* dev, uint64_t value) {
    if (!dev || !dev->pass_timeline || value == 0) return true;

    VkSemaphoreWaitInfo wait_info = {0};
    wait_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO;
    wait_info.semaphoreCount = 1;
    wait_info.pSemaphores = &dev->pass_timeline;
    wait_info.pValues = &value;
    return md_vk_check(vkWaitSemaphores(dev->device, &wait_info, UINT64_MAX), "vkWaitSemaphores(pass timeline)");
}

static void md_vk_reclaim_completed_passes(md_gpu_device* dev) {
    if (!dev || !dev->pass_timeline) return;
    uint64_t completed = md_vk_pass_timeline_value(dev);
    for (uint32_t i = 0; i < MD_GPU_PASS_HISTORY_SIZE; ++i) {
        md_gpu_pass_record* record = &dev->pass_records[i];
        if (record->id == 0 || record->id > completed || !record->queue || !record->cmd) continue;
        md_vk_reset_cmd_state(record->cmd, record->queue);
        ASSERT(record->queue->free_count < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE);
        record->queue->free_list[record->queue->free_count++] = record->cmd->index;
        record->id = 0;
        record->queue = NULL;
        record->cmd = NULL;
    }
}

static md_gpu_command_buffer_t md_gpu_command_buffer_acquire(md_gpu_queue_t queue) {
    if (!queue) return NULL;
    md_gpu_queue* q = (md_gpu_queue*)queue;

    md_vk_reclaim_completed_passes(q->dev);

    if (q->free_count == 0) {
        MD_LOG_ERROR("md_gpu_command_buffer_acquire: no free command buffers");
        return NULL;
    }

    uint32_t idx = q->free_list[--q->free_count];

    md_gpu_command_buffer* cmd = &q->cmd_wrappers[idx];
    /* Reset transient recording state; dev, vk_cmd, dispatch_sets, and index are permanent */
    md_vk_reset_cmd_state(cmd, q);

    // Begin command buffer
    VkCommandBufferBeginInfo begin_info = {0};
    begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    if (!md_vk_check(vkBeginCommandBuffer(cmd->vk_cmd, &begin_info), "vkBeginCommandBuffer")) {
        q->free_list[q->free_count++] = idx;
        return NULL;
    }

    cmd->is_recording = true;

    return (md_gpu_command_buffer_t)cmd;
}

static void md_vk_transition_image(VkCommandBuffer vk_cmd, md_gpu_image* img, VkImageLayout new_layout, VkPipelineStageFlags2 src_stage, VkAccessFlags2 src_access, VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    if (img->layout == new_layout) return;

    VkImageMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
    barrier.srcStageMask = src_stage;
    barrier.srcAccessMask = src_access;
    barrier.dstStageMask = dst_stage;
    barrier.dstAccessMask = dst_access;
    barrier.oldLayout = img->layout;
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

    vkCmdPipelineBarrier2(vk_cmd, &dep);
    img->layout = new_layout;
}

static void md_vk_apply_pass_declared_bindings(md_gpu_command_buffer* cmd) {
    memset(cmd->bound_buffers, 0, sizeof(cmd->bound_buffers));
    memset(cmd->bound_buffer_offsets, 0, sizeof(cmd->bound_buffer_offsets));
    memset(cmd->bound_buffer_sizes, 0, sizeof(cmd->bound_buffer_sizes));
    memset(cmd->bound_images, 0, sizeof(cmd->bound_images));
    memset(cmd->bound_sampled_images, 0, sizeof(cmd->bound_sampled_images));
    memset(cmd->bound_samplers, 0, sizeof(cmd->bound_samplers));

    const uint32_t buffer_count = cmd->pass_buffer_count < MD_GPU_MAX_BIND_SLOTS ? cmd->pass_buffer_count : MD_GPU_MAX_BIND_SLOTS;
    for (uint32_t i = 0; i < buffer_count; ++i) {
        cmd->bound_buffers[i] = cmd->pass_buffers[i].buffer;
    }

    const uint32_t image_count = cmd->pass_image_count < MD_GPU_MAX_BIND_SLOTS ? cmd->pass_image_count : MD_GPU_MAX_BIND_SLOTS;
    for (uint32_t i = 0; i < image_count; ++i) {
        cmd->bound_images[i] = cmd->pass_images[i].image;
    }

    const uint32_t sampled_resource_count = cmd->pass_sampled_resource_count < MD_GPU_MAX_BIND_SLOTS ? cmd->pass_sampled_resource_count : MD_GPU_MAX_BIND_SLOTS;
    for (uint32_t i = 0; i < sampled_resource_count; ++i) {
        cmd->bound_sampled_images[i] = cmd->pass_sampled_resources[i].image;
        cmd->bound_samplers[i] = cmd->pass_sampled_resources[i].sampler;
    }
}

// Helper: flush current bound state into the given descriptor set
static void md_vk_flush_descriptor_set(md_gpu_command_buffer* cmd, VkDescriptorSet ds) {
    md_gpu_device* dev = cmd->dev;

    VkWriteDescriptorSet writes[4 * MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorBufferInfo buffer_infos[MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorImageInfo storage_image_infos[MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorImageInfo sampled_image_infos[MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorImageInfo sampler_infos[MD_GPU_MAX_BIND_SLOTS];
    uint32_t write_count = 0;

    for (uint32_t i = 0; i < MD_GPU_MAX_BIND_SLOTS; ++i) {
        if (cmd->bound_buffers[i]) {
            md_gpu_buffer* buf = (md_gpu_buffer*)cmd->bound_buffers[i];
            buffer_infos[i].buffer = buf->buffer;
            buffer_infos[i].offset = cmd->bound_buffer_offsets[i];
            buffer_infos[i].range = (cmd->bound_buffer_sizes[i] > 0) ? cmd->bound_buffer_sizes[i] : VK_WHOLE_SIZE;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = ds;
            writes[write_count].dstBinding = i;
            writes[write_count].dstArrayElement = 0;
            writes[write_count].descriptorCount = 1;
            writes[write_count].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            writes[write_count].pBufferInfo = &buffer_infos[i];
            writes[write_count].pImageInfo = NULL;
            writes[write_count].pTexelBufferView = NULL;
            write_count++;
        }

        if (cmd->bound_images[i]) {
            md_gpu_image* img = (md_gpu_image*)cmd->bound_images[i];

            // Transition to GENERAL for storage use
            md_vk_transition_image(cmd->vk_cmd, img, VK_IMAGE_LAYOUT_GENERAL,
                                    VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                                    VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT);

            storage_image_infos[i].sampler = VK_NULL_HANDLE;
            storage_image_infos[i].imageView = img->view;
            storage_image_infos[i].imageLayout = VK_IMAGE_LAYOUT_GENERAL;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = ds;
            writes[write_count].dstBinding = MD_GPU_STORAGE_IMAGE_BINDING_BASE + i;
            writes[write_count].dstArrayElement = 0;
            writes[write_count].descriptorCount = 1;
            writes[write_count].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            writes[write_count].pBufferInfo = NULL;
            writes[write_count].pImageInfo = &storage_image_infos[i];
            writes[write_count].pTexelBufferView = NULL;
            write_count++;
        }

        if (cmd->bound_sampled_images[i]) {
            md_gpu_image* img = (md_gpu_image*)cmd->bound_sampled_images[i];

            md_vk_transition_image(cmd->vk_cmd, img, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
                                    VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                                    VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT, VK_ACCESS_2_SHADER_READ_BIT);

            sampled_image_infos[i].sampler = VK_NULL_HANDLE;
            sampled_image_infos[i].imageView = img->view;
            sampled_image_infos[i].imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = ds;
            writes[write_count].dstBinding = MD_GPU_SAMPLED_IMAGE_BINDING_BASE + i;
            writes[write_count].dstArrayElement = 0;
            writes[write_count].descriptorCount = 1;
            writes[write_count].descriptorType = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
            writes[write_count].pBufferInfo = NULL;
            writes[write_count].pImageInfo = &sampled_image_infos[i];
            writes[write_count].pTexelBufferView = NULL;
            write_count++;
        }

        if (cmd->bound_samplers[i]) {
            md_gpu_sampler* sampler = (md_gpu_sampler*)cmd->bound_samplers[i];

            sampler_infos[i].sampler = sampler->sampler;
            sampler_infos[i].imageView = VK_NULL_HANDLE;
            sampler_infos[i].imageLayout = VK_IMAGE_LAYOUT_UNDEFINED;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = ds;
            writes[write_count].dstBinding = MD_GPU_SAMPLER_BINDING_BASE + i;
            writes[write_count].dstArrayElement = 0;
            writes[write_count].descriptorCount = 1;
            writes[write_count].descriptorType = VK_DESCRIPTOR_TYPE_SAMPLER;
            writes[write_count].pBufferInfo = NULL;
            writes[write_count].pImageInfo = &sampler_infos[i];
            writes[write_count].pTexelBufferView = NULL;
            write_count++;
        }
    }

    if (write_count > 0) {
        vkUpdateDescriptorSets(dev->device, write_count, writes, 0, NULL);
    }
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

static bool md_vk_queue_submit_pass(md_gpu_queue* q, md_gpu_command_buffer* c, uint64_t signal_value) {
    if (!q || !c || signal_value == 0) return false;
    md_gpu_device* dev = q->dev;
    ASSERT(c->queue == q);
    ASSERT(c->is_recording);

    if (!md_vk_check(vkEndCommandBuffer(c->vk_cmd), "vkEndCommandBuffer(pass)")) {
        return false;
    }
    c->is_recording = false;

    VkCommandBufferSubmitInfo cmd_info = {0};
    cmd_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO;
    cmd_info.commandBuffer = c->vk_cmd;

    VkSemaphoreSubmitInfo signal_info = {0};
    signal_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO;
    signal_info.semaphore = dev->pass_timeline;
    signal_info.value = signal_value;
    signal_info.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;

    VkSubmitInfo2 submit = {0};
    submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO_2;
    submit.commandBufferInfoCount = 1;
    submit.pCommandBufferInfos = &cmd_info;
    submit.signalSemaphoreInfoCount = 1;
    submit.pSignalSemaphoreInfos = &signal_info;

    if (!md_vk_check(vkQueueSubmit2(q->vk_queue, 1, &submit, VK_NULL_HANDLE), "vkQueueSubmit2(pass timeline)")) {
        md_vk_reset_cmd_state(c, q);
        ASSERT(q->free_count < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE);
        q->free_list[q->free_count++] = c->index;
        return false;
    }

    return true;
}

md_gpu_pass_t md_gpu_begin(md_gpu_device_t device, const char* label) {
    if (!device) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;
    md_gpu_queue* queue = &dev->compute_queue;

    md_gpu_pass* pass = (md_gpu_pass*)calloc(1, sizeof(md_gpu_pass));
    if (!pass) return NULL;

    pass->queue = queue;
    pass->cmd = (md_gpu_command_buffer*)md_gpu_command_buffer_acquire((md_gpu_queue_t)queue);
    if (!pass->cmd) {
        free(pass);
        return NULL;
    }

    if (label) {
        VkDebugUtilsLabelEXT info = {0};
        info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT;
        info.pLabelName = label;
        if (vkCmdBeginDebugUtilsLabelEXT) vkCmdBeginDebugUtilsLabelEXT(pass->cmd->vk_cmd, &info);
        pass->pushed_debug_group = true;
    }
    return (md_gpu_pass_t)pass;
}

md_gpu_pass_id_t md_gpu_end(md_gpu_pass_t pass) {
    md_gpu_pass_id_t id = {0};
    if (!pass) return id;

    md_gpu_pass* p = (md_gpu_pass*)pass;
    md_gpu_queue* queue = p->queue;
    md_gpu_device* dev = queue ? queue->dev : NULL;
    if (!queue || !dev || !p->cmd) {
        free(p);
        return id;
    }

    if (p->pushed_debug_group && vkCmdEndDebugUtilsLabelEXT) {
        vkCmdEndDebugUtilsLabelEXT(p->cmd->vk_cmd);
    }

    md_vk_reclaim_completed_passes(dev);

    md_gpu_pass_record* record = &dev->pass_records[dev->pass_record_head++ % MD_GPU_PASS_HISTORY_SIZE];
    if (record->id != 0) {
        md_vk_wait_pass_timeline(dev, record->id);
        md_vk_reclaim_completed_passes(dev);
    }
    record->id = 0;
    record->queue = NULL;
    record->cmd = NULL;

    id.value = ++dev->next_pass_id;
    if (id.value == 0) id.value = ++dev->next_pass_id;

    bool ok = md_vk_queue_submit_pass(queue, p->cmd, id.value);
    if (ok) {
        record->id = id.value;
        record->queue = queue;
        record->cmd = p->cmd;
    } else {
        id.value = 0;
    }

    free(p);
    return id;
}

bool md_gpu_is_complete(md_gpu_device_t device, md_gpu_pass_id_t id) {
    if (!device || id.value == 0) return true;
    md_gpu_device* dev = (md_gpu_device*)device;
    md_vk_reclaim_completed_passes(dev);
    return md_vk_pass_timeline_value(dev) >= id.value;
}

void md_gpu_wait(md_gpu_device_t device, md_gpu_pass_id_t id) {
    if (!device || id.value == 0) return;
    md_gpu_device* dev = (md_gpu_device*)device;
    if (md_vk_wait_pass_timeline(dev, id.value)) {
        md_vk_reclaim_completed_passes(dev);
    }
}

void md_gpu_wait_for(md_gpu_pass_t pass, md_gpu_pass_id_t id) {
    (void)pass;
    (void)id;
}

void md_gpu_use_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer, gpu_usage_flags_t usage) {
    if (!pass || !buffer || usage == 0) return;
    md_vk_pass_merge_buffer_usage(((md_gpu_pass*)pass)->cmd, buffer, usage);
}

void md_gpu_use_image(md_gpu_pass_t pass, md_gpu_image_t image, gpu_usage_flags_t usage) {
    if (!pass || !image || usage == 0) return;
    md_vk_pass_merge_image_usage(((md_gpu_pass*)pass)->cmd, image, usage);
}

void md_gpu_use_sampled_resource(md_gpu_pass_t pass, md_gpu_image_t image, md_gpu_sampler_t sampler, gpu_usage_flags_t usage) {
    if (!pass || !image) return;
    md_vk_pass_merge_sampled_resource_usage(((md_gpu_pass*)pass)->cmd, image, sampler, usage ? usage : GPU_USAGE_READ);
}

void md_gpu_bind_compute_pipeline(md_gpu_pass_t pass, md_gpu_compute_pipeline_t pipeline) {
    if (!pass || !pipeline) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    ASSERT(c->is_recording);
    c->bound_pipeline = pipeline;
}

void md_gpu_dispatch(md_gpu_pass_t pass, uint32_t group_count_x, uint32_t group_count_y, uint32_t group_count_z, const void* root_args, size_t root_args_size) {
    if (!pass) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_device* dev = c->dev;
    md_gpu_compute_pipeline* pipeline = (md_gpu_compute_pipeline*)c->bound_pipeline;

    ASSERT(c->is_recording);
    ASSERT(c->bound_pipeline && "Pipeline must be bound before dispatch");
    ASSERT(c->dispatch_count < MD_GPU_VK_MAX_DISPATCHES_PER_CMD && "Too many dispatches per command buffer");

    if (root_args && root_args_size) {
        ASSERT(root_args_size <= MD_GPU_MAX_PUSH_CONSTANTS);
        memcpy(c->push_constants, root_args, root_args_size);
        c->push_constant_size = root_args_size;
    }

    md_vk_apply_pass_declared_bindings(c);

    vkCmdBindPipeline(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline->pipeline);

    if (c->push_constant_size > 0) {
        vkCmdPushConstants(c->vk_cmd, dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, (uint32_t)c->push_constant_size, c->push_constants);
    }

    VkDescriptorSet ds = c->dispatch_sets[c->dispatch_count];
    md_vk_flush_descriptor_set(c, ds);
    vkCmdBindDescriptorSets(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, dev->pipeline_layout, 0, 1, &ds, 0, NULL);
    vkCmdDispatch(c->vk_cmd, group_count_x, group_count_y, group_count_z);

    c->dispatch_count++;
}

void md_gpu_barrier(md_gpu_pass_t pass, md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;

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
}

void md_gpu_barrier_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer,
                                md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass || !buffer) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;

    VkBufferMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER_2;
    barrier.srcStageMask = md_vk_stage_flags(src_stage);
    barrier.srcAccessMask = md_vk_access_flags_write(src_stage);
    barrier.dstStageMask = md_vk_stage_flags(dst_stage);
    barrier.dstAccessMask = md_vk_access_flags_read(dst_stage);
    barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.buffer = buf->buffer;
    barrier.offset = 0;
    barrier.size = VK_WHOLE_SIZE;

    VkDependencyInfo dep = {0};
    dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    dep.bufferMemoryBarrierCount = 1;
    dep.pBufferMemoryBarriers = &barrier;
    vkCmdPipelineBarrier2(c->vk_cmd, &dep);
}

void md_gpu_barrier_image(md_gpu_pass_t pass, md_gpu_image_t image,
                               md_gpu_barrier_stage_t src_stage, md_gpu_barrier_stage_t dst_stage) {
    if (!pass || !image) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_image* img = (md_gpu_image*)image;

    VkImageMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
    barrier.srcStageMask = md_vk_stage_flags(src_stage);
    barrier.srcAccessMask = md_vk_access_flags_write(src_stage);
    barrier.dstStageMask = md_vk_stage_flags(dst_stage);
    barrier.dstAccessMask = md_vk_access_flags_read(dst_stage);
    barrier.oldLayout = img->layout;
    barrier.newLayout = img->layout;
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
    vkCmdPipelineBarrier2(c->vk_cmd, &dep);
}

void md_gpu_push_debug_group(md_gpu_pass_t pass, const char* label) {
    if (!pass || !label) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    VkDebugUtilsLabelEXT info = {0};
    info.sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT;
    info.pLabelName = label;
    if (vkCmdBeginDebugUtilsLabelEXT) vkCmdBeginDebugUtilsLabelEXT(c->vk_cmd, &info);
}

void md_gpu_pop_debug_group(md_gpu_pass_t pass) {
    if (!pass) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    if (vkCmdEndDebugUtilsLabelEXT) vkCmdEndDebugUtilsLabelEXT(c->vk_cmd);
}

void md_gpu_copy_buffer(md_gpu_pass_t pass, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!pass || !src || !dst) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_buffer* src_buf = (md_gpu_buffer*)src;
    md_gpu_buffer* dst_buf = (md_gpu_buffer*)dst;

    VkBufferCopy region = {0};
    region.srcOffset = src_offset;
    region.dstOffset = dst_offset;
    region.size = size;
    vkCmdCopyBuffer(c->vk_cmd, src_buf->buffer, dst_buf->buffer, 1, &region);
}

void md_gpu_copy_image_region_to_buffer(md_gpu_pass_t pass,
                                             md_gpu_image_t src_image,
                                             md_gpu_image_region_t src_region,
                                             md_gpu_buffer_t dst_buffer,
                                             size_t dst_offset) {
    if (!pass || !src_image || !dst_buffer) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_image* src = (md_gpu_image*)src_image;
    md_gpu_buffer* dst = (md_gpu_buffer*)dst_buffer;

    md_vk_transition_image(c->vk_cmd, src, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
                           VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);

    VkBufferImageCopy region = {0};
    region.bufferOffset = dst_offset;
    region.bufferRowLength = 0;
    region.bufferImageHeight = 0;
    region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region.imageSubresource.mipLevel = 0;
    region.imageSubresource.baseArrayLayer = 0;
    region.imageSubresource.layerCount = 1;
    region.imageOffset.x = (int32_t)src_region.offset[0];
    region.imageOffset.y = (int32_t)src_region.offset[1];
    region.imageOffset.z = (int32_t)src_region.offset[2];
    region.imageExtent.width = src_region.extent[0] ? src_region.extent[0] : src->width;
    region.imageExtent.height = src_region.extent[1] ? src_region.extent[1] : src->height;
    region.imageExtent.depth = src_region.extent[2] ? src_region.extent[2] : src->depth;
    vkCmdCopyImageToBuffer(c->vk_cmd, src->image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, dst->buffer, 1, &region);
}

void md_gpu_copy_buffer_to_image(md_gpu_pass_t pass, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image) {
    if (!pass || !src_buffer || !dst_image) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_buffer* src = (md_gpu_buffer*)src_buffer;
    md_gpu_image* dst = (md_gpu_image*)dst_image;

    md_vk_transition_image(c->vk_cmd, dst, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                           VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                           VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    VkBufferImageCopy region = {0};
    region.bufferOffset = 0;
    region.bufferRowLength = 0;
    region.bufferImageHeight = 0;
    region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region.imageSubresource.mipLevel = 0;
    region.imageSubresource.baseArrayLayer = 0;
    region.imageSubresource.layerCount = 1;
    region.imageOffset.x = 0;
    region.imageOffset.y = 0;
    region.imageOffset.z = 0;
    region.imageExtent.width = dst->width;
    region.imageExtent.height = dst->height;
    region.imageExtent.depth = dst->depth;
    vkCmdCopyBufferToImage(c->vk_cmd, src->buffer, dst->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);
}

void md_gpu_copy_buffer_to_image_region(md_gpu_pass_t pass,
                                             md_gpu_buffer_t src_buffer,
                                             size_t src_offset,
                                             md_gpu_image_t dst_image,
                                             md_gpu_image_region_t dst_region) {
    if (!pass || !src_buffer || !dst_image) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_buffer* src = (md_gpu_buffer*)src_buffer;
    md_gpu_image* dst = (md_gpu_image*)dst_image;

    md_vk_transition_image(c->vk_cmd, dst, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                           VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                           VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    VkBufferImageCopy region = {0};
    region.bufferOffset = src_offset;
    region.bufferRowLength = 0;
    region.bufferImageHeight = 0;
    region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region.imageSubresource.mipLevel = 0;
    region.imageSubresource.baseArrayLayer = 0;
    region.imageSubresource.layerCount = 1;
    region.imageOffset.x = (int32_t)dst_region.offset[0];
    region.imageOffset.y = (int32_t)dst_region.offset[1];
    region.imageOffset.z = (int32_t)dst_region.offset[2];
    region.imageExtent.width = dst_region.extent[0] ? dst_region.extent[0] : dst->width;
    region.imageExtent.height = dst_region.extent[1] ? dst_region.extent[1] : dst->height;
    region.imageExtent.depth = dst_region.extent[2] ? dst_region.extent[2] : dst->depth;
    vkCmdCopyBufferToImage(c->vk_cmd, src->buffer, dst->image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);
}

void md_gpu_fill_buffer(md_gpu_pass_t pass, md_gpu_buffer_t buffer, size_t offset, size_t size, uint8_t value) {
    if (!pass || !buffer) return;
    md_gpu_command_buffer* c = ((md_gpu_pass*)pass)->cmd;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;

    ASSERT((offset % 4 == 0) && "vkCmdFillBuffer offset must be 4-byte aligned");
    ASSERT((size == VK_WHOLE_SIZE || size % 4 == 0) && "vkCmdFillBuffer size must be 4-byte aligned or VK_WHOLE_SIZE");

    uint32_t pattern = (value << 24) | (value << 16) | (value << 8) | value;
    vkCmdFillBuffer(c->vk_cmd, buf->buffer, offset, size, pattern);
}

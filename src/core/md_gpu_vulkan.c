#include "md_gpu.h"

#include <core/md_common.h>
#include <core/md_log.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#if !defined(MD_GPU_VULKAN)
// Allow forcing this compilation unit in non-Vulkan builds, but it should only
// be compiled/linked when Vulkan is available.
#define MD_GPU_VULKAN 1
#endif

// We avoid including platform-specific Vulkan headers via GLFW here.
// This is a compute-only backend for now, so no surface extensions are required.
#include <vulkan/vulkan.h>

#ifndef MD_ENABLE_GPU
#define MD_ENABLE_GPU 0
#endif

// Vulkan backend (C).
// Compute-only for now; designed to be extended with graphics later.

// Command buffer pool size (not to be confused with Metal's MD_GPU_MAX_COMMANDS,
// which is the max number of recorded commands per command buffer)
#define MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE 64

typedef struct md_gpu_device {
    VkInstance instance;
    VkPhysicalDevice physical_device;
    VkDevice device;

    uint32_t compute_queue_family;
    VkQueue compute_queue;

    VkCommandPool command_pool;
    VkDescriptorPool descriptor_pool;
    VkDescriptorSetLayout descriptor_set_layout;
    VkPipelineLayout pipeline_layout;

    // Command buffer pool
    VkCommandBuffer command_buffers[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    uint32_t command_buffer_free_list[MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE];
    uint32_t command_buffer_free_count;
    pthread_mutex_t cmd_pool_mutex;
} md_gpu_device;

typedef struct md_gpu_queue {
    md_gpu_device_t device;
} md_gpu_queue;

typedef struct md_gpu_buffer {
    md_gpu_device_t device;
    VkBuffer buffer;
    VkDeviceMemory memory;
    VkDeviceSize size;
    md_gpu_buffer_flags_t flags;
    void* mapped;
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

typedef struct md_gpu_compute_pipeline {
    md_gpu_device_t device;
    VkPipeline pipeline;
    VkShaderModule shader_module;
} md_gpu_compute_pipeline;

typedef struct md_gpu_command_buffer {
    md_gpu_device_t device;
    VkCommandBuffer vk_cmd;
    VkDescriptorSet descriptor_set;
    uint32_t index; // for recycling

    // Bound state
    md_gpu_compute_pipeline_t bound_pipeline;
    md_gpu_buffer_t bound_buffers[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_offsets[MD_GPU_MAX_BIND_SLOTS];
    size_t bound_buffer_sizes[MD_GPU_MAX_BIND_SLOTS];
    md_gpu_image_t bound_images[MD_GPU_MAX_BIND_SLOTS];
    bool descriptor_set_dirty;
} md_gpu_command_buffer;

typedef struct md_gpu_fence {
    md_gpu_device_t device;
    VkFence vk_fence;
    uint32_t cmd_buffer_index;
    VkDescriptorSet descriptor_set;
} md_gpu_fence;

static bool md_vk_check(VkResult res, const char* what) {
    if (res == VK_SUCCESS) return true;
    MD_LOG_ERROR("Vulkan error %d: %s", (int)res, what ? what : "(unknown)");
    return false;
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

    VkPhysicalDeviceVulkan13Features f13 = {0};
    f13.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES;
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

    vkGetDeviceQueue(out_dev->device, out_dev->compute_queue_family, 0, &out_dev->compute_queue);
    if (out_dev->compute_queue == VK_NULL_HANDLE) return false;

    // Create command pool
    VkCommandPoolCreateInfo pool_ci = {0};
    pool_ci.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    pool_ci.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pool_ci.queueFamilyIndex = out_dev->compute_queue_family;
    if (!md_vk_check(vkCreateCommandPool(out_dev->device, &pool_ci, NULL, &out_dev->command_pool), "vkCreateCommandPool")) return false;

    // Allocate command buffers
    VkCommandBufferAllocateInfo cb_ai = {0};
    cb_ai.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    cb_ai.commandPool = out_dev->command_pool;
    cb_ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cb_ai.commandBufferCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;
    if (!md_vk_check(vkAllocateCommandBuffers(out_dev->device, &cb_ai, out_dev->command_buffers), "vkAllocateCommandBuffers")) return false;

    // Initialize free list
    for (uint32_t i = 0; i < MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE; ++i) {
        out_dev->command_buffer_free_list[i] = i;
    }
    out_dev->command_buffer_free_count = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;

    // Create descriptor set layout: buffers at binding 0..(MAX_BIND_SLOTS-1), images at MAX_BIND_SLOTS..(2*MAX_BIND_SLOTS-1)
    VkDescriptorSetLayoutBinding bindings[2 * MD_GPU_MAX_BIND_SLOTS];
    for (uint32_t i = 0; i < MD_GPU_MAX_BIND_SLOTS; ++i) {
        bindings[i].binding = i;
        bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[i].pImmutableSamplers = NULL;

        bindings[MD_GPU_MAX_BIND_SLOTS + i].binding = MD_GPU_MAX_BIND_SLOTS + i;
        bindings[MD_GPU_MAX_BIND_SLOTS + i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
        bindings[MD_GPU_MAX_BIND_SLOTS + i].descriptorCount = 1;
        bindings[MD_GPU_MAX_BIND_SLOTS + i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[MD_GPU_MAX_BIND_SLOTS + i].pImmutableSamplers = NULL;
    }

    VkDescriptorSetLayoutCreateInfo dsl_ci = {0};
    dsl_ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    dsl_ci.bindingCount = 2 * MD_GPU_MAX_BIND_SLOTS;
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
    VkDescriptorPoolSize pool_sizes[2];
    pool_sizes[0].type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    pool_sizes[0].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_MAX_BIND_SLOTS;
    pool_sizes[1].type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    pool_sizes[1].descriptorCount = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE * MD_GPU_MAX_BIND_SLOTS;

    VkDescriptorPoolCreateInfo dp_ci = {0};
    dp_ci.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    dp_ci.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    dp_ci.maxSets = MD_GPU_VK_COMMAND_BUFFER_POOL_SIZE;
    dp_ci.poolSizeCount = 2;
    dp_ci.pPoolSizes = pool_sizes;
    if (!md_vk_check(vkCreateDescriptorPool(out_dev->device, &dp_ci, NULL, &out_dev->descriptor_pool), "vkCreateDescriptorPool")) return false;

    return true;
}

md_gpu_device_t md_gpu_create_device(void) {
    md_gpu_device* dev = (md_gpu_device*)calloc(1, sizeof(md_gpu_device));
    if (!dev) return NULL;

    if (!md_vk_create_instance(dev)) {
        free(dev);
        return NULL;
    }

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

    return (md_gpu_device_t)dev;
}

void md_gpu_destroy_device(md_gpu_device_t device) {
    if (!device) return;
    md_gpu_device* dev = (md_gpu_device*)device;
    if (dev->device) {
        vkDeviceWaitIdle(dev->device);
        if (dev->descriptor_pool) vkDestroyDescriptorPool(dev->device, dev->descriptor_pool, NULL);
        if (dev->pipeline_layout) vkDestroyPipelineLayout(dev->device, dev->pipeline_layout, NULL);
        if (dev->descriptor_set_layout) vkDestroyDescriptorSetLayout(dev->device, dev->descriptor_set_layout, NULL);
        if (dev->command_pool) vkDestroyCommandPool(dev->device, dev->command_pool, NULL);
        vkDestroyDevice(dev->device, NULL);
    }
    if (dev->instance) {
        vkDestroyInstance(dev->instance, NULL);
    }
    pthread_mutex_destroy(&dev->cmd_pool_mutex);
    free(dev);
}

md_gpu_queue_t md_gpu_acquire_compute_queue(md_gpu_device_t device) {
    if (!device) return NULL;
    md_gpu_queue* q = (md_gpu_queue*)calloc(1, sizeof(md_gpu_queue));
    if (!q) return NULL;
    q->device = device;
    return (md_gpu_queue_t)q;
}

md_gpu_buffer_t md_gpu_create_buffer(md_gpu_device_t device, const md_gpu_buffer_desc_t* desc) {
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
    bci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT;
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
        required |= VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT;
        preferred |= VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
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

    VkMemoryAllocateInfo mai = {0};
    mai.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
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

    // Optional persistent map for CPU-visible buffers.
    if (desc->flags & MD_GPU_BUFFER_CPU_VISIBLE) {
        void* mapped = NULL;
        VkResult res = vkMapMemory(dev->device, buf->memory, 0, VK_WHOLE_SIZE, 0, &mapped);
        if (res == VK_SUCCESS) {
            buf->mapped = mapped;
        } else {
            // Mapping is optional; allow later map call to retry.
            (void)res;
        }
    }

    return (md_gpu_buffer_t)buf;
}

void md_gpu_destroy_buffer(md_gpu_buffer_t buffer) {
    if (!buffer) return;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    md_gpu_device* dev = (md_gpu_device*)buf->device;
    if (!dev) {
        free(buf);
        return;
    }

    if (buf->mapped) {
        vkUnmapMemory(dev->device, buf->memory);
        buf->mapped = NULL;
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

void* md_gpu_map_buffer(md_gpu_buffer_t buffer) {
    if (!buffer) return NULL;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    if (!(buf->flags & MD_GPU_BUFFER_CPU_VISIBLE)) {
        MD_LOG_ERROR("md_gpu_map_buffer called on non-CPU-visible buffer");
        return NULL;
    }
    if (buf->mapped) return buf->mapped;

    md_gpu_device* dev = (md_gpu_device*)buf->device;
    void* mapped = NULL;
    if (!md_vk_check(vkMapMemory(dev->device, buf->memory, 0, VK_WHOLE_SIZE, 0, &mapped), "vkMapMemory")) {
        return NULL;
    }
    buf->mapped = mapped;
    return mapped;
}

void md_gpu_unmap_buffer(md_gpu_buffer_t buffer) {
    if (!buffer) return;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;
    if (!(buf->flags & MD_GPU_BUFFER_CPU_VISIBLE)) return;
    if (!buf->mapped) return;

    // We keep buffers persistently mapped for simplicity; unmap is a no-op.
    // If you later prefer explicit unmap, switch this to vkUnmapMemory.
}

md_gpu_image_t md_gpu_create_image(md_gpu_device_t device, const md_gpu_image_desc_t* desc) {
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
    ici.usage = VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    if (desc->flags & MD_GPU_IMAGE_STORAGE) {
        ici.usage |= VK_IMAGE_USAGE_STORAGE_BIT;
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

void md_gpu_destroy_image(md_gpu_image_t image) {
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

md_gpu_compute_pipeline_t md_gpu_create_compute_pipeline(md_gpu_device_t device, const md_gpu_compute_pipeline_desc_t* desc) {
    if (!device || !desc || !desc->shader.data || desc->shader.size == 0) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;

    md_gpu_compute_pipeline* pipeline = (md_gpu_compute_pipeline*)calloc(1, sizeof(md_gpu_compute_pipeline));
    if (!pipeline) return NULL;
    pipeline->device = device;

    // Create shader module
    VkShaderModuleCreateInfo sm_ci = {0};
    sm_ci.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    sm_ci.codeSize = desc->shader.size;
    sm_ci.pCode = (const uint32_t*)desc->shader.data;

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

void md_gpu_destroy_compute_pipeline(md_gpu_compute_pipeline_t pipeline) {
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

md_gpu_command_buffer_t md_gpu_acquire_command_buffer(md_gpu_device_t device) {
    if (!device) return NULL;
    md_gpu_device* dev = (md_gpu_device*)device;
    
    pthread_mutex_lock(&dev->cmd_pool_mutex);
    
    if (dev->command_buffer_free_count == 0) {
        pthread_mutex_unlock(&dev->cmd_pool_mutex);
        MD_LOG_ERROR("md_gpu_acquire_command_buffer: no free command buffers");
        return NULL;
    }

    uint32_t idx = dev->command_buffer_free_list[--dev->command_buffer_free_count];
    pthread_mutex_unlock(&dev->cmd_pool_mutex);

    md_gpu_command_buffer* cmd = (md_gpu_command_buffer*)calloc(1, sizeof(md_gpu_command_buffer));
    if (!cmd) {
        pthread_mutex_lock(&dev->cmd_pool_mutex);
        dev->command_buffer_free_list[dev->command_buffer_free_count++] = idx;
        pthread_mutex_unlock(&dev->cmd_pool_mutex);
        return NULL;
    }

    cmd->device = device;
    cmd->vk_cmd = dev->command_buffers[idx];
    cmd->index = idx;

    // Allocate descriptor set
    VkDescriptorSetAllocateInfo dsa_i = {0};
    dsa_i.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    dsa_i.descriptorPool = dev->descriptor_pool;
    dsa_i.descriptorSetCount = 1;
    dsa_i.pSetLayouts = &dev->descriptor_set_layout;

    if (!md_vk_check(vkAllocateDescriptorSets(dev->device, &dsa_i, &cmd->descriptor_set), "vkAllocateDescriptorSets")) {
        pthread_mutex_lock(&dev->cmd_pool_mutex);
        dev->command_buffer_free_list[dev->command_buffer_free_count++] = idx;
        pthread_mutex_unlock(&dev->cmd_pool_mutex);
        free(cmd);
        return NULL;
    }

    // Begin command buffer
    VkCommandBufferBeginInfo begin_info = {0};
    begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    if (!md_vk_check(vkBeginCommandBuffer(cmd->vk_cmd, &begin_info), "vkBeginCommandBuffer")) {
        vkFreeDescriptorSets(dev->device, dev->descriptor_pool, 1, &cmd->descriptor_set);
        pthread_mutex_lock(&dev->cmd_pool_mutex);
        dev->command_buffer_free_list[dev->command_buffer_free_count++] = idx;
        pthread_mutex_unlock(&dev->cmd_pool_mutex);
        free(cmd);
        return NULL;
    }

    return (md_gpu_command_buffer_t)cmd;
}

// Helper: transition image layout if needed
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

// Helper: flush descriptor updates
static void md_vk_flush_descriptor_set(md_gpu_command_buffer* cmd) {
    if (!cmd->descriptor_set_dirty) return;
    md_gpu_device* dev = (md_gpu_device*)cmd->device;

    VkWriteDescriptorSet writes[2 * MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorBufferInfo buffer_infos[MD_GPU_MAX_BIND_SLOTS];
    VkDescriptorImageInfo image_infos[MD_GPU_MAX_BIND_SLOTS];
    uint32_t write_count = 0;

    for (uint32_t i = 0; i < MD_GPU_MAX_BIND_SLOTS; ++i) {
        if (cmd->bound_buffers[i]) {
            md_gpu_buffer* buf = (md_gpu_buffer*)cmd->bound_buffers[i];
            buffer_infos[i].buffer = buf->buffer;
            buffer_infos[i].offset = cmd->bound_buffer_offsets[i];
            buffer_infos[i].range = (cmd->bound_buffer_sizes[i] > 0) ? cmd->bound_buffer_sizes[i] : VK_WHOLE_SIZE;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = cmd->descriptor_set;
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

            image_infos[i].sampler = VK_NULL_HANDLE;
            image_infos[i].imageView = img->view;
            image_infos[i].imageLayout = VK_IMAGE_LAYOUT_GENERAL;

            writes[write_count].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            writes[write_count].pNext = NULL;
            writes[write_count].dstSet = cmd->descriptor_set;
            writes[write_count].dstBinding = MD_GPU_MAX_BIND_SLOTS + i;
            writes[write_count].dstArrayElement = 0;
            writes[write_count].descriptorCount = 1;
            writes[write_count].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
            writes[write_count].pBufferInfo = NULL;
            writes[write_count].pImageInfo = &image_infos[i];
            writes[write_count].pTexelBufferView = NULL;
            write_count++;
        }
    }

    if (write_count > 0) {
        vkUpdateDescriptorSets(dev->device, write_count, writes, 0, NULL);
    }

    cmd->descriptor_set_dirty = false;
}

void md_gpu_cmd_bind_compute_pipeline(md_gpu_command_buffer_t cmd, md_gpu_compute_pipeline_t pipeline) {
    if (!cmd || !pipeline) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_compute_pipeline* p = (md_gpu_compute_pipeline*)pipeline;
    c->bound_pipeline = pipeline;
    vkCmdBindPipeline(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, p->pipeline);
}

void md_gpu_cmd_bind_buffer(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer) {
    if (!cmd) return;
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    c->bound_buffers[slot] = buffer;
    c->bound_buffer_offsets[slot] = 0;
    c->bound_buffer_sizes[slot] = 0;
    c->descriptor_set_dirty = true;
}

void md_gpu_cmd_bind_buffer_range(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_buffer_t buffer, size_t offset, size_t size) {
    if (!cmd) return;
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    c->bound_buffers[slot] = buffer;
    c->bound_buffer_offsets[slot] = offset;
    c->bound_buffer_sizes[slot] = size;
    c->descriptor_set_dirty = true;
}

void md_gpu_cmd_bind_image(md_gpu_command_buffer_t cmd, uint32_t slot, md_gpu_image_t image) {
    if (!cmd) return;
    ASSERT(slot < MD_GPU_MAX_BIND_SLOTS);
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    c->bound_images[slot] = image;
    c->descriptor_set_dirty = true;
}

void md_gpu_cmd_push_constants(md_gpu_command_buffer_t cmd, const void* data, size_t size) {
    if (!cmd || !data || size == 0) return;
    ASSERT(size <= MD_GPU_MAX_PUSH_CONSTANTS);
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_device* dev = (md_gpu_device*)c->device;
    vkCmdPushConstants(c->vk_cmd, dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, (uint32_t)size, data);
}

void md_gpu_cmd_dispatch(md_gpu_command_buffer_t cmd, uint32_t total_x, uint32_t total_y, uint32_t total_z) {
    if (!cmd) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_device* dev = (md_gpu_device*)c->device;

    ASSERT(c->bound_pipeline && "Pipeline must be bound before dispatch");

    // Flush descriptor set if dirty (includes implicit image transitions)
    md_vk_flush_descriptor_set(c);

    // Bind descriptor set
    vkCmdBindDescriptorSets(c->vk_cmd, VK_PIPELINE_BIND_POINT_COMPUTE, dev->pipeline_layout, 0, 1, &c->descriptor_set, 0, NULL);

    vkCmdDispatch(c->vk_cmd, total_x, total_y, total_z);
}

void md_gpu_cmd_barrier(md_gpu_command_buffer_t cmd) {
    if (!cmd) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;

    VkMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER_2;
    barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.srcAccessMask = VK_ACCESS_2_MEMORY_WRITE_BIT;
    barrier.dstStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;

    VkDependencyInfo dep = {0};
    dep.sType = VK_STRUCTURE_TYPE_DEPENDENCY_INFO;
    dep.memoryBarrierCount = 1;
    dep.pMemoryBarriers = &barrier;

    vkCmdPipelineBarrier2(c->vk_cmd, &dep);
}

void md_gpu_cmd_barrier_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer) {
    if (!cmd || !buffer) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;

    VkBufferMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER_2;
    barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.srcAccessMask = VK_ACCESS_2_MEMORY_WRITE_BIT;
    barrier.dstStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
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

void md_gpu_cmd_barrier_image(md_gpu_command_buffer_t cmd, md_gpu_image_t image) {
    if (!cmd || !image) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_image* img = (md_gpu_image*)image;

    VkImageMemoryBarrier2 barrier = {0};
    barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2;
    barrier.srcStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.srcAccessMask = VK_ACCESS_2_MEMORY_WRITE_BIT;
    barrier.dstStageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    barrier.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
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

void md_gpu_cmd_copy_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t src, md_gpu_buffer_t dst, size_t size, size_t src_offset, size_t dst_offset) {
    if (!cmd || !src || !dst) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_buffer* src_buf = (md_gpu_buffer*)src;
    md_gpu_buffer* dst_buf = (md_gpu_buffer*)dst;

    VkBufferCopy region = {0};
    region.srcOffset = src_offset;
    region.dstOffset = dst_offset;
    region.size = size;

    vkCmdCopyBuffer(c->vk_cmd, src_buf->buffer, dst_buf->buffer, 1, &region);
}

void md_gpu_cmd_copy_image_to_buffer(md_gpu_command_buffer_t cmd, md_gpu_image_t src_image, md_gpu_buffer_t dst_buffer) {
    if (!cmd || !src_image || !dst_buffer) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_image* src = (md_gpu_image*)src_image;
    md_gpu_buffer* dst = (md_gpu_buffer*)dst_buffer;

    // Transition to TRANSFER_SRC_OPTIMAL
    md_vk_transition_image(c->vk_cmd, src, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                            VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
                            VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_READ_BIT);

    VkBufferImageCopy region = {0};
    region.bufferOffset = 0;
    region.bufferRowLength = 0; // tight packing
    region.bufferImageHeight = 0;
    region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region.imageSubresource.mipLevel = 0;
    region.imageSubresource.baseArrayLayer = 0;
    region.imageSubresource.layerCount = 1;
    region.imageOffset.x = 0;
    region.imageOffset.y = 0;
    region.imageOffset.z = 0;
    region.imageExtent.width = src->width;
    region.imageExtent.height = src->height;
    region.imageExtent.depth = src->depth;

    vkCmdCopyImageToBuffer(c->vk_cmd, src->image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, dst->buffer, 1, &region);
}

void md_gpu_cmd_copy_buffer_to_image(md_gpu_command_buffer_t cmd, md_gpu_buffer_t src_buffer, md_gpu_image_t dst_image) {
    if (!cmd || !src_buffer || !dst_image) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_buffer* src = (md_gpu_buffer*)src_buffer;
    md_gpu_image* dst = (md_gpu_image*)dst_image;

    // Transition to TRANSFER_DST_OPTIMAL
    md_vk_transition_image(c->vk_cmd, dst, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                            VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                            VK_PIPELINE_STAGE_2_COPY_BIT, VK_ACCESS_2_TRANSFER_WRITE_BIT);

    VkBufferImageCopy region = {0};
    region.bufferOffset = 0;
    region.bufferRowLength = 0; // tight packing
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

void md_gpu_cmd_fill_buffer(md_gpu_command_buffer_t cmd, md_gpu_buffer_t buffer, size_t offset, size_t size, uint32_t value) {
    if (!cmd || !buffer) return;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_buffer* buf = (md_gpu_buffer*)buffer;

    ASSERT((offset % 4 == 0) && "vkCmdFillBuffer offset must be 4-byte aligned");
    ASSERT((size == VK_WHOLE_SIZE || size % 4 == 0) && "vkCmdFillBuffer size must be 4-byte aligned or VK_WHOLE_SIZE");

    vkCmdFillBuffer(c->vk_cmd, buf->buffer, offset, size, value);
}

md_gpu_fence_t md_gpu_queue_submit(md_gpu_queue_t queue, md_gpu_command_buffer_t cmd) {
    if (!queue || !cmd) return NULL;
    md_gpu_queue* q = (md_gpu_queue*)queue;
    md_gpu_command_buffer* c = (md_gpu_command_buffer*)cmd;
    md_gpu_device* dev = (md_gpu_device*)q->device;

    // End command buffer
    if (!md_vk_check(vkEndCommandBuffer(c->vk_cmd), "vkEndCommandBuffer")) {
        return NULL;
    }

    // Create fence
    md_gpu_fence* fence = (md_gpu_fence*)calloc(1, sizeof(md_gpu_fence));
    if (!fence) return NULL;
    fence->device = q->device;

    VkFenceCreateInfo fence_ci = {0};
    fence_ci.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    if (!md_vk_check(vkCreateFence(dev->device, &fence_ci, NULL, &fence->vk_fence), "vkCreateFence")) {
        free(fence);
        return NULL;
    }

    // Submit
    VkCommandBufferSubmitInfo cmd_info = {0};
    cmd_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO;
    cmd_info.commandBuffer = c->vk_cmd;

    VkSubmitInfo2 submit = {0};
    submit.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO_2;
    submit.commandBufferInfoCount = 1;
    submit.pCommandBufferInfos = &cmd_info;

    if (!md_vk_check(vkQueueSubmit2(dev->compute_queue, 1, &submit, fence->vk_fence), "vkQueueSubmit2")) {
        vkDestroyFence(dev->device, fence->vk_fence, NULL);
        free(fence);
        return NULL;
    }

    // Store command buffer resources in fence for later recycling (after GPU completes)
    fence->cmd_buffer_index = c->index;
    fence->descriptor_set = c->descriptor_set;
    free(c);

    return (md_gpu_fence_t)fence;
}

bool md_gpu_fence_is_signaled(md_gpu_fence_t fence) {
    if (!fence) return true;
    md_gpu_fence* f = (md_gpu_fence*)fence;
    md_gpu_device* dev = (md_gpu_device*)f->device;
    VkResult res = vkGetFenceStatus(dev->device, f->vk_fence);
    return (res == VK_SUCCESS);
}

void md_gpu_fence_wait(md_gpu_fence_t fence) {
    if (!fence) return;
    md_gpu_fence* f = (md_gpu_fence*)fence;
    md_gpu_device* dev = (md_gpu_device*)f->device;
    vkWaitForFences(dev->device, 1, &f->vk_fence, VK_TRUE, UINT64_MAX);
}

void md_gpu_destroy_fence(md_gpu_fence_t fence) {
    if (!fence) return;
    md_gpu_fence* f = (md_gpu_fence*)fence;
    md_gpu_device* dev = (md_gpu_device*)f->device;
    if (!dev) {
        free(f);
        return;
    }
    if (f->vk_fence) {
        // Wait for GPU to finish using the command buffer before recycling resources
        vkWaitForFences(dev->device, 1, &f->vk_fence, VK_TRUE, UINT64_MAX);
        
        // Now safe to recycle command buffer and descriptor set
        vkFreeDescriptorSets(dev->device, dev->descriptor_pool, 1, &f->descriptor_set);
        
        pthread_mutex_lock(&dev->cmd_pool_mutex);
        dev->command_buffer_free_list[dev->command_buffer_free_count++] = f->cmd_buffer_index;
        pthread_mutex_unlock(&dev->cmd_pool_mutex);
        
        vkDestroyFence(dev->device, f->vk_fence, NULL);
        f->vk_fence = VK_NULL_HANDLE;
    }
    free(f);
}

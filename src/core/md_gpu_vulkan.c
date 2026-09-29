/*
md_gpu_vulkan.c -- Vulkan backend for md_gpu.h

Structure of this file:

    1.  Configuration, error handling, small utilities
    2.  Types
    3.  Allocation registry (device address -> allocation lookup)
    4.  Raw buffers and transient arenas (argument blocks and staging)
    5.  Deferred destruction
    6.  Device creation
    7.  Bindless descriptor set
    8.  Streams, command buffers, submission, ordering
    9.  Memory: pools, malloc/free, copies, uploads
    10. Textures and samplers
    11. Kernels and launches
    12. Host callbacks and polling
    13. Device destruction

The dependency model is program order within a stream, implemented as a single
global VkMemoryBarrier2 between consecutive operations in a command buffer
(IMPLICIT ordering), or as caller-placed stage barriers (EXPLICIT ordering).
There is no per-resource state tracking anywhere in this file, and every image
lives in VK_IMAGE_LAYOUT_GENERAL for its entire life.

Nothing here blocks the calling thread except md_gpu_stream_sync,
md_gpu_sync_wait, md_gpu_stream_destroy (on its own stream) and device
creation/destruction. Destroying textures, kernels and pools is deferred: the
object goes onto a retire list stamped with every stream's current position
and is released by md_gpu_device_poll once all of those have passed. That is
what lets a compute job run for many frames without anything else stalling
behind it.

Resources are created VK_SHARING_MODE_CONCURRENT across the compute and
transfer families whenever those differ, so a buffer or image written on one
and read on the other needs no queue-family ownership transfer.
*/

#include "md_gpu.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VOLK_IMPLEMENTATION
#include <volk.h>

#include "md_gpu_builtin_spv.inl"

/* =========================================================================
   1. Configuration, error handling, utilities
   ========================================================================= */

/* ---- Bindless heap -------------------------------------------------------
   Slang's DescriptorHandle<T> lowers to a uint2 whose .x component indexes a
   descriptor heap. With -bindless-space-index N, Slang emits one runtime array
   per resource type a shader actually uses, in set N.

   Which binding each type lands on is Slang's BindlessDescriptorOptions. The
   default, VkMutable, aliases every type onto one binding and therefore needs
   VK_DESCRIPTOR_TYPE_MUTABLE_EXT. md_gpu.slang overrides it to `None`, which
   gives one binding per descriptor type -- no aliasing across types, no
   extension, and it runs on hardware older than Turing. Sampled images of
   different dimensionality still share binding 2, which is ordinary
   descriptor indexing: they are the same VkDescriptorType.

   Keep MD_VK_BINDLESS_SPACE in sync with the -bindless-space-index flag in
   cmake/CompileGpuShaders.cmake, and this table in sync with the comment in
   src/shaders/md_gpu.slang. */
#define MD_VK_BINDLESS_SPACE         0u
#define MD_VK_BINDING_SAMPLER        0u   /* VK_DESCRIPTOR_TYPE_SAMPLER       */
#define MD_VK_BINDING_SAMPLED_IMAGE  2u   /* VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE */
#define MD_VK_BINDING_STORAGE_IMAGE  3u   /* VK_DESCRIPTOR_TYPE_STORAGE_IMAGE */
#define MD_VK_BINDING_COUNT          3u   /* how many of the above we declare */

#define MD_VK_MAX_QUEUES_PER_FAMILY  8u
#define MD_VK_MAX_TEXTURE_SLOTS   4096u   /* sampled + storage views share it */
#define MD_VK_MAX_SAMPLERS         256u
#define MD_VK_ARENA_PAGE_SIZE   (256u * 1024u)
#define MD_VK_ARG_ALIGN            64u
#define MD_VK_ERROR_BUF           512u

#if defined(_MSC_VER)
#define MD_VK_THREAD_LOCAL __declspec(thread)
#else
#define MD_VK_THREAD_LOCAL __thread
#endif

static MD_VK_THREAD_LOCAL char md_vk_error_buf[MD_VK_ERROR_BUF];
static MD_VK_THREAD_LOCAL bool md_vk_has_error;

static bool md_vk_fail(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(md_vk_error_buf, sizeof(md_vk_error_buf), fmt, ap);
    va_end(ap);
    md_vk_has_error = true;
    MD_LOG_ERROR("md_gpu: %s", md_vk_error_buf);
    return false;
}

const char* md_gpu_last_error(void) {
    return md_vk_has_error ? md_vk_error_buf : NULL;
}

static const char* md_vk_result_str(VkResult r) {
    switch (r) {
    case VK_ERROR_OUT_OF_HOST_MEMORY:      return "VK_ERROR_OUT_OF_HOST_MEMORY";
    case VK_ERROR_OUT_OF_DEVICE_MEMORY:    return "VK_ERROR_OUT_OF_DEVICE_MEMORY";
    case VK_ERROR_INITIALIZATION_FAILED:   return "VK_ERROR_INITIALIZATION_FAILED";
    case VK_ERROR_DEVICE_LOST:             return "VK_ERROR_DEVICE_LOST";
    case VK_ERROR_MEMORY_MAP_FAILED:       return "VK_ERROR_MEMORY_MAP_FAILED";
    case VK_ERROR_LAYER_NOT_PRESENT:       return "VK_ERROR_LAYER_NOT_PRESENT";
    case VK_ERROR_EXTENSION_NOT_PRESENT:   return "VK_ERROR_EXTENSION_NOT_PRESENT";
    case VK_ERROR_FEATURE_NOT_PRESENT:     return "VK_ERROR_FEATURE_NOT_PRESENT";
    case VK_ERROR_INCOMPATIBLE_DRIVER:     return "VK_ERROR_INCOMPATIBLE_DRIVER";
    case VK_ERROR_TOO_MANY_OBJECTS:        return "VK_ERROR_TOO_MANY_OBJECTS";
    case VK_ERROR_FORMAT_NOT_SUPPORTED:    return "VK_ERROR_FORMAT_NOT_SUPPORTED";
    case VK_ERROR_FRAGMENTED_POOL:         return "VK_ERROR_FRAGMENTED_POOL";
    case VK_ERROR_OUT_OF_POOL_MEMORY:      return "VK_ERROR_OUT_OF_POOL_MEMORY";
    case VK_ERROR_INVALID_EXTERNAL_HANDLE: return "VK_ERROR_INVALID_EXTERNAL_HANDLE";
    case VK_ERROR_FRAGMENTATION:           return "VK_ERROR_FRAGMENTATION";
    case VK_ERROR_UNKNOWN:                 return "VK_ERROR_UNKNOWN";
    default:                               return "VkResult";
    }
}

static bool md_vk_check(VkResult r, const char* what) {
    if (r == VK_SUCCESS) return true;
    return md_vk_fail("%s failed: %s (%d)", what, md_vk_result_str(r), (int)r);
}

static inline uint64_t md_vk_align_up(uint64_t v, uint64_t a) {
    return (v + a - 1) & ~(a - 1);
}

static inline uint64_t md_vk_next_pow2(uint64_t v) {
    if (v < 256) return 256;
    v--;
    v |= v >> 1;  v |= v >> 2;  v |= v >> 4;
    v |= v >> 8;  v |= v >> 16; v |= v >> 32;
    return v + 1;
}

/* Minimal growable array of pointers / structs. */
typedef struct md_vk_vec_t {
    void*  data;
    size_t count;
    size_t capacity;
    size_t stride;
} md_vk_vec_t;

static void md_vk_vec_init(md_vk_vec_t* v, size_t stride) {
    v->data = NULL; v->count = 0; v->capacity = 0; v->stride = stride;
}

static bool md_vk_vec_reserve(md_vk_vec_t* v, struct md_allocator_i* alloc, size_t n) {
    if (n <= v->capacity) return true;
    size_t cap = v->capacity ? v->capacity * 2 : 16;
    while (cap < n) cap *= 2;
    void* mem = md_alloc(alloc, cap * v->stride);
    if (!mem) return false;
    if (v->data) {
        memcpy(mem, v->data, v->count * v->stride);
        md_free(alloc, v->data, v->capacity * v->stride);
    }
    v->data = mem;
    v->capacity = cap;
    return true;
}

static void* md_vk_vec_push(md_vk_vec_t* v, struct md_allocator_i* alloc) {
    if (!md_vk_vec_reserve(v, alloc, v->count + 1)) return NULL;
    void* slot = (char*)v->data + v->count * v->stride;
    memset(slot, 0, v->stride);
    v->count++;
    return slot;
}

/* Remove element i, preserving order. */
static void md_vk_vec_remove(md_vk_vec_t* v, size_t i) {
    char* base = (char*)v->data;
    memmove(base + i * v->stride, base + (i + 1) * v->stride, (v->count - i - 1) * v->stride);
    v->count--;
}

/* Remove the first element equal to the pointer `p` (for vectors of pointers). */
static void md_vk_vec_remove_ptr(md_vk_vec_t* v, const void* p) {
    void** arr = (void**)v->data;
    for (size_t i = 0; i < v->count; ++i) {
        if (arr[i] == p) { md_vk_vec_remove(v, i); return; }
    }
}

static void md_vk_vec_free(md_vk_vec_t* v, struct md_allocator_i* alloc) {
    if (v->data) md_free(alloc, v->data, v->capacity * v->stride);
    v->data = NULL; v->count = 0; v->capacity = 0;
}

#define MD_VK_VEC_AT(v, type, i) (((type*)(v).data)[i])

/* =========================================================================
   2. Types
   ========================================================================= */

typedef struct md_vk_block_t {
    VkBuffer          buffer;
    VkDeviceMemory    memory;
    uint64_t          address;      /* device address of byte 0            */
    void*             host;         /* mapped pointer, or NULL             */
    uint64_t          capacity;     /* actual allocated size               */
    uint64_t          size;         /* size requested by the current user  */
    md_gpu_mem_kind_t kind;
    md_gpu_pool_t     pool;
    bool              in_use;
    /* Stream-ordered free bookkeeping: the block becomes reusable at
       free_value on free_stream, or immediately for that same stream. */
    md_gpu_stream_t   free_stream;
    uint64_t          free_value;
} md_vk_block_t;

typedef struct md_gpu_pool {
    md_gpu_device_t   device;
    md_gpu_mem_kind_t kind;          /* the one kind this pool serves */
    uint64_t          cache_limit;   /* 0 = unlimited                 */
    md_vk_vec_t       blocks;        /* md_vk_block_t*                */
    md_vk_vec_t       textures;      /* md_gpu_texture_t (live)       */
    uint64_t          in_use_bytes;
    uint64_t          reserved_bytes;
    uint64_t          peak_in_use_bytes;
    uint64_t          alloc_count;
    uint64_t          reuse_count;
    char              label[64];
} md_gpu_pool;

/* One page of host-visible, device-addressable transient memory. */
typedef struct md_vk_page_t {
    VkBuffer       buffer;
    VkDeviceMemory memory;
    uint64_t       address;
    uint8_t*       host;
    uint64_t       capacity;
    uint64_t       cursor;
    uint64_t       retire_value;   /* stream timeline value that frees it */
} md_vk_page_t;

typedef struct md_vk_arena_t {
    md_vk_vec_t pages;     /* md_vk_page_t* */
    size_t      current;   /* index of the page being filled */
} md_vk_arena_t;

typedef struct md_vk_cmd_t {
    VkCommandBuffer cmd;
    uint64_t        value;      /* timeline value signalled by its submit */
    bool            pending;
} md_vk_cmd_t;

typedef struct md_vk_wait_t {
    md_gpu_stream_t stream;
    uint64_t        value;
} md_vk_wait_t;

typedef struct md_gpu_stream {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    uint32_t             family;
    bool                 can_compute;      /* false on a transfer-only family */
    VkQueue              queue;
    VkSemaphore          timeline;
    uint64_t             next_value;       /* value the next submit signals */
    uint64_t             submitted_value;  /* last value submitted          */

    VkCommandPool        cmd_pool;
    md_vk_vec_t          cmds;             /* md_vk_cmd_t */
    VkCommandBuffer      open;             /* currently recording, or NULL  */
    bool                 has_work;
    bool                 needs_barrier;
    md_gpu_ordering_t    ordering;

    md_vk_vec_t          waits;            /* md_vk_wait_t, pending for the next submit */

    md_vk_arena_t        arena;

    /* Open upload, if any. */
    bool                 upload_open;
    bool                 upload_direct;
    md_gpu_addr_t        upload_dst;
    uint64_t             upload_src_addr;
    size_t               upload_size;

    bool                 is_default;
    char                 label[64];
} md_gpu_stream;

typedef struct md_vk_format_info_t {
    VkFormat           format;
    uint32_t           bytes;          /* texel size in buffer copies         */
    VkImageAspectFlags view_aspect;    /* aspect of sampled/storage views     */
    VkImageAspectFlags copy_aspect;    /* aspect moved by buffer<->image copy */
    bool               depth;
    const char*        name;
} md_vk_format_info_t;

typedef struct md_gpu_texture {
    md_gpu_device_t       device;
    md_gpu_pool_t         pool;
    VkImage               image;
    VkDeviceMemory        memory;
    uint64_t              bytes;           /* memory size, for pool stats    */
    VkImageView           sampled_view;    /* all mips, or VK_NULL_HANDLE    */
    uint32_t              sampled_slot;    /* heap slot, or 0                */
    VkImageView*          storage_views;   /* one per mip, or NULL           */
    uint32_t*             storage_slots;   /* one per mip, or NULL           */
    md_vk_format_info_t   fi;
    md_gpu_texture_desc_t desc;            /* normalised                     */
    char                  label[64];
} md_gpu_texture;

typedef struct md_vk_sampler_entry_t {
    md_gpu_sampler_desc_t desc;
    VkSampler             sampler;
    uint32_t              slot;
} md_vk_sampler_entry_t;

typedef struct md_gpu_kernel {
    md_gpu_device_t device;
    VkShaderModule  module;
    VkPipeline      pipeline;
    uint32_t        group_size[3];
    uint32_t        args_size;
    char            label[64];
} md_gpu_kernel;

typedef struct md_vk_hostfn_t {
    md_gpu_sync_t  sync;
    md_gpu_host_fn fn;
    void*          user;
} md_vk_hostfn_t;

typedef enum md_vk_retire_kind_t {
    MD_VK_RETIRE_BLOCK,
    MD_VK_RETIRE_TEXTURE,
    MD_VK_RETIRE_KERNEL,
} md_vk_retire_kind_t;

/* An object whose destruction waits for every stream to pass the point at
   which it was destroyed. Owns `waits`. */
typedef struct md_vk_retire_t {
    md_vk_retire_kind_t kind;
    void*               object;     /* md_vk_block_t* / md_gpu_texture_t / md_gpu_kernel_t */
    md_vk_wait_t*       waits;
    uint32_t            wait_count;
    uint32_t            wait_capacity;   /* what md_free must be told */
} md_vk_retire_t;

/* Everything the backend asks of a device, in one place, so that a device we
   cannot use is rejected by name instead of failing later inside
   vkCreateDevice with a bare VkResult. The optional flags are enabled only
   when the driver reports them. */
typedef struct md_vk_dev_caps_t {
    bool maintenance4;
    bool update_unused_while_pending;
    bool nonuniform_storage_image;
    bool nonuniform_sampled_image;
    bool dynamic_storage_image;
    bool dynamic_sampled_image;
    bool shader_int64;
} md_vk_dev_caps_t;

typedef struct md_gpu_device {
    struct md_allocator_i* alloc;
    VkInstance             instance;
    VkPhysicalDevice       phys;
    VkDevice               device;
    VkDebugUtilsMessengerEXT messenger;
    VkPhysicalDeviceMemoryProperties mem_props;
    VkPhysicalDeviceProperties       props;
    uint32_t                         subgroup_size;

    uint32_t compute_family, transfer_family;
    bool     transfer_can_compute;
    /* Families a resource must be shared across (CONCURRENT when 2). */
    uint32_t share_families[2];
    uint32_t share_family_count;

    /* Streams are spread round-robin over these. A device that exposes several
       compute queues (typical: 8 on NVIDIA) then gets real concurrency from
       "use another stream" rather than just permission to overlap. */
    VkQueue  compute_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t compute_queue_count;
    VkQueue  transfer_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t transfer_queue_count;
    uint32_t next_compute_queue, next_transfer_queue;
    md_mutex_t queue_mutex;      /* streams may share a VkQueue */
    md_mutex_t device_mutex;     /* pools, registry, heap, retire and callback lists */

    /* Bindless */
    VkDescriptorSetLayout set_layout;
    VkDescriptorPool      desc_pool;
    VkDescriptorSet       desc_set;
    VkPipelineLayout      pipeline_layout;

    /* A 1x1x1 placeholder. Freed descriptor slots are pointed at it so that no
       descriptor ever references a destroyed view. */
    VkImage         dummy_image;
    VkImageView     dummy_view;
    VkDeviceMemory  dummy_mem;
    VkSampler       dummy_sampler;

    /* Texture heap slots; slot 0 is never handed out so zero stays null. */
    uint32_t        tex_free[MD_VK_MAX_TEXTURE_SLOTS];
    uint32_t        tex_free_count;

    md_vk_sampler_entry_t samplers[MD_VK_MAX_SAMPLERS];
    uint32_t              sampler_count;   /* slots 1..sampler_count used */

    /* Address -> allocation registry, sorted ascending by address. */
    md_vk_vec_t     registry;   /* md_vk_block_t* */

    md_vk_vec_t     pools;      /* md_gpu_pool_t   */
    md_vk_vec_t     kernels;    /* md_gpu_kernel_t */
    md_vk_vec_t     streams;    /* md_gpu_stream_t */
    md_gpu_stream_t default_compute;
    md_gpu_stream_t default_transfer;

    md_vk_vec_t     hostfns;    /* md_vk_hostfn_t */
    md_vk_vec_t     retires;    /* md_vk_retire_t */

    md_gpu_kernel_t make_grid_kernel;

    md_vk_dev_caps_t caps;

    bool            is_discrete;
    bool            validation;
} md_gpu_device;

/* Forward declarations. */
static bool     md_vk_stream_ensure_cmd(md_gpu_stream_t s);
static bool     md_vk_stream_submit(md_gpu_stream_t s);
static uint64_t md_vk_stream_completed(md_gpu_stream_t s);
static bool     md_vk_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host,
                                  VkBuffer* out_buffer, uint64_t* out_offset);
static void     md_vk_texture_free(md_gpu_device_t dev, md_gpu_texture_t t);
static void     md_vk_kernel_free(md_gpu_device_t dev, md_gpu_kernel_t k);
static void     md_vk_block_free(md_gpu_device_t dev, md_vk_block_t* b);

/* =========================================================================
   3. Allocation registry
   ========================================================================= */

/* Binary search for the block whose [address, address+capacity) contains a.
   Caller holds device_mutex. */
static md_vk_block_t* md_vk_registry_find_locked(md_gpu_device_t dev, uint64_t address) {
    size_t lo = 0, hi = dev->registry.count;
    md_vk_block_t** arr = (md_vk_block_t**)dev->registry.data;
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        md_vk_block_t* b = arr[mid];
        if (address < b->address) {
            hi = mid;
        } else if (address >= b->address + b->capacity) {
            lo = mid + 1;
        } else {
            return b;
        }
    }
    return NULL;
}

static bool md_vk_registry_insert_locked(md_gpu_device_t dev, md_vk_block_t* blk) {
    if (!md_vk_vec_reserve(&dev->registry, dev->alloc, dev->registry.count + 1)) return false;
    md_vk_block_t** arr = (md_vk_block_t**)dev->registry.data;
    size_t i = dev->registry.count;
    while (i > 0 && arr[i - 1]->address > blk->address) {
        arr[i] = arr[i - 1];
        --i;
    }
    arr[i] = blk;
    dev->registry.count++;
    return true;
}

static void md_vk_registry_remove_locked(md_gpu_device_t dev, md_vk_block_t* blk) {
    md_vk_vec_remove_ptr(&dev->registry, blk);
}

/* Resolve [addr, addr + size) to a live allocation. Takes the device lock for
   the lookup: the registry is mutated by malloc on other threads, so an
   unlocked binary search is a data race. The returned block stays valid for
   as long as the caller's use of the memory is legal, i.e. until it frees it.
   `what` names the calling function for the error message. */
static md_vk_block_t* md_vk_resolve(md_gpu_device_t dev, md_gpu_addr_t addr, uint64_t size,
                                    uint64_t* out_offset, const char* what) {
    md_mutex_lock(&dev->device_mutex);
    md_vk_block_t* b = md_vk_registry_find_locked(dev, addr);
    md_mutex_unlock(&dev->device_mutex);
    if (!b || !b->in_use) {
        md_vk_fail("%s: 0x%llx is not a live md_gpu allocation", what, (unsigned long long)addr);
        return NULL;
    }
    uint64_t off = addr - b->address;
    if (size > b->size || off > b->size - size) {
        md_vk_fail("%s: range [0x%llx, +%llu) overruns its %llu-byte allocation",
                   what, (unsigned long long)addr, (unsigned long long)size, (unsigned long long)b->size);
        return NULL;
    }
    if (out_offset) *out_offset = off;
    return b;
}

/* =========================================================================
   4. Raw buffers and transient arenas
   ========================================================================= */

static uint32_t md_vk_find_memory_type(md_gpu_device_t dev, uint32_t type_bits, VkMemoryPropertyFlags required, VkMemoryPropertyFlags preferred) {
    uint32_t best = UINT32_MAX;
    for (uint32_t i = 0; i < dev->mem_props.memoryTypeCount; ++i) {
        if (!(type_bits & (1u << i))) continue;
        VkMemoryPropertyFlags f = dev->mem_props.memoryTypes[i].propertyFlags;
        if ((f & required) != required) continue;
        if (preferred && (f & preferred) == preferred) return i;
        if (best == UINT32_MAX) best = i;
    }
    return best;
}

static void md_vk_set_sharing(md_gpu_device_t dev, VkSharingMode* mode, uint32_t* count, const uint32_t** families) {
    if (dev->share_family_count > 1) {
        *mode     = VK_SHARING_MODE_CONCURRENT;
        *count    = dev->share_family_count;
        *families = dev->share_families;
    } else {
        *mode     = VK_SHARING_MODE_EXCLUSIVE;
        *count    = 0;
        *families = NULL;
    }
}

static bool md_vk_create_raw_buffer(md_gpu_device_t dev, uint64_t size, md_gpu_mem_kind_t kind,
                                    VkBuffer* out_buf, VkDeviceMemory* out_mem, uint64_t* out_addr, void** out_host) {
    VkBufferCreateInfo bci = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bci.size  = size;
    bci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
              | VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT
              | VK_BUFFER_USAGE_TRANSFER_SRC_BIT
              | VK_BUFFER_USAGE_TRANSFER_DST_BIT
              | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT
              | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    md_vk_set_sharing(dev, &bci.sharingMode, &bci.queueFamilyIndexCount, &bci.pQueueFamilyIndices);

    VkBuffer buf = VK_NULL_HANDLE;
    if (!md_vk_check(vkCreateBuffer(dev->device, &bci, NULL, &buf), "vkCreateBuffer")) return false;

    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(dev->device, buf, &req);

    const bool host_visible = kind != MD_GPU_MEM_DEVICE;
    VkMemoryPropertyFlags required = 0, preferred = 0;
    if (host_visible) {
        required = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        if (kind == MD_GPU_MEM_HOST_READ) preferred = required | VK_MEMORY_PROPERTY_HOST_CACHED_BIT;
        else                              preferred = required | VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    } else {
        required = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    }

    uint32_t type = md_vk_find_memory_type(dev, req.memoryTypeBits, required, preferred);
    if (type == UINT32_MAX && !host_visible) {
        /* UMA parts may expose no purely device-local type. */
        type = md_vk_find_memory_type(dev, req.memoryTypeBits, 0, 0);
    }
    if (type == UINT32_MAX) {
        vkDestroyBuffer(dev->device, buf, NULL);
        return md_vk_fail("no suitable memory type for %llu bytes (kind %d)", (unsigned long long)size, (int)kind);
    }

    VkMemoryAllocateFlagsInfo fi = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO};
    fi.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;

    VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.pNext           = &fi;
    mai.allocationSize  = req.size;
    mai.memoryTypeIndex = type;

    VkDeviceMemory mem = VK_NULL_HANDLE;
    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &mem), "vkAllocateMemory")) {
        vkDestroyBuffer(dev->device, buf, NULL);
        return false;
    }
    if (!md_vk_check(vkBindBufferMemory(dev->device, buf, mem, 0), "vkBindBufferMemory")) {
        vkFreeMemory(dev->device, mem, NULL);
        vkDestroyBuffer(dev->device, buf, NULL);
        return false;
    }

    void* host = NULL;
    if (host_visible) {
        if (!md_vk_check(vkMapMemory(dev->device, mem, 0, VK_WHOLE_SIZE, 0, &host), "vkMapMemory")) {
            vkFreeMemory(dev->device, mem, NULL);
            vkDestroyBuffer(dev->device, buf, NULL);
            return false;
        }
    }

    VkBufferDeviceAddressInfo bdai = {VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO};
    bdai.buffer = buf;
    uint64_t addr = vkGetBufferDeviceAddress(dev->device, &bdai);

    *out_buf = buf; *out_mem = mem; *out_addr = addr; *out_host = host;
    return true;
}

static void md_vk_destroy_raw_buffer(md_gpu_device_t dev, VkBuffer buf, VkDeviceMemory mem) {
    if (buf) vkDestroyBuffer(dev->device, buf, NULL);
    if (mem) vkFreeMemory(dev->device, mem, NULL);
}

static md_vk_page_t* md_vk_page_create(md_gpu_device_t dev, uint64_t size) {
    md_vk_page_t* p = (md_vk_page_t*)md_alloc(dev->alloc, sizeof(md_vk_page_t));
    if (!p) return NULL;
    memset(p, 0, sizeof(*p));
    void* host = NULL;
    if (!md_vk_create_raw_buffer(dev, size, MD_GPU_MEM_HOST_WRITE, &p->buffer, &p->memory, &p->address, &host)) {
        md_free(dev->alloc, p, sizeof(md_vk_page_t));
        return NULL;
    }
    p->host     = (uint8_t*)host;
    p->capacity = size;
    return p;
}

static void md_vk_page_destroy(md_gpu_device_t dev, md_vk_page_t* p) {
    md_vk_destroy_raw_buffer(dev, p->buffer, p->memory);
    md_free(dev->alloc, p, sizeof(md_vk_page_t));
}

static void md_vk_page_take(md_vk_page_t* p, uint64_t need, uint64_t* out_addr, void** out_host,
                            VkBuffer* out_buffer, uint64_t* out_offset) {
    *out_addr = p->address + p->cursor;
    if (out_host)   *out_host   = p->host + p->cursor;
    if (out_buffer) *out_buffer = p->buffer;
    if (out_offset) *out_offset = p->cursor;
    p->cursor += need;
    /* The page now holds unsubmitted data. Drop any stamp from an earlier,
       completed submission, or a later allocation would see the page as
       retired and recycle it (cursor back to 0) under this data before it is
       even submitted. md_vk_arena_retire restamps it at the next submit. */
    p->retire_value = 0;
}

/* Reserve `size` bytes of transient, device-addressable, host-writable memory
   from the stream's arena. Valid until the submission that consumes it
   completes. Returns the page's VkBuffer and offset too, so that staging never
   has to be looked up again by address. */
static bool md_vk_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host,
                              VkBuffer* out_buffer, uint64_t* out_offset) {
    md_gpu_device_t dev = s->device;
    md_vk_arena_t*  a   = &s->arena;
    uint64_t need = md_vk_align_up(size, MD_VK_ARG_ALIGN);

    /* Try the current page. */
    if (a->pages.count > 0) {
        md_vk_page_t* p = MD_VK_VEC_AT(a->pages, md_vk_page_t*, a->current);
        if (p->cursor + need <= p->capacity) {
            md_vk_page_take(p, need, out_addr, out_host, out_buffer, out_offset);
            return true;
        }
    }

    /* Look for a retired page with room. */
    uint64_t done = md_vk_stream_completed(s);
    for (size_t i = 0; i < a->pages.count; ++i) {
        md_vk_page_t* p = MD_VK_VEC_AT(a->pages, md_vk_page_t*, i);
        if (p->retire_value != 0 && p->retire_value <= done) {
            p->cursor = 0;
            p->retire_value = 0;
        }
        if (p->retire_value == 0 && p->cursor + need <= p->capacity) {
            a->current = i;
            md_vk_page_take(p, need, out_addr, out_host, out_buffer, out_offset);
            return true;
        }
    }

    /* Allocate a new page. */
    uint64_t page_size = MD_VK_ARENA_PAGE_SIZE;
    while (page_size < need) page_size *= 2;
    md_vk_page_t* p = md_vk_page_create(dev, page_size);
    if (!p) return md_vk_fail("failed to allocate a %llu byte transient page", (unsigned long long)page_size);
    md_vk_page_t** slot = (md_vk_page_t**)md_vk_vec_push(&a->pages, dev->alloc);
    if (!slot) { md_vk_page_destroy(dev, p); return md_vk_fail("out of memory"); }
    *slot = p;
    a->current = a->pages.count - 1;
    md_vk_page_take(p, need, out_addr, out_host, out_buffer, out_offset);
    return true;
}

/* Stamp every page holding data with the value that releases it.

   This must overwrite an existing stamp, not skip it. A page can be appended
   to across several submissions -- the fast path in md_vk_arena_alloc carves
   more out of the current page after it has already been submitted, which is
   safe in itself. But if the stamp were left at the *first* submission's
   value, the page would be recycled as soon as that one completed, while a
   later submission was still reading its argument blocks. Always taking the
   newest value keeps the page alive until everything that touched it has
   finished. */
static void md_vk_arena_retire(md_gpu_stream_t s, uint64_t value) {
    md_vk_arena_t* a = &s->arena;
    for (size_t i = 0; i < a->pages.count; ++i) {
        md_vk_page_t* p = MD_VK_VEC_AT(a->pages, md_vk_page_t*, i);
        if (p->cursor > 0) p->retire_value = value;
    }
}

static void md_vk_arena_free(md_gpu_device_t dev, md_vk_arena_t* a) {
    for (size_t i = 0; i < a->pages.count; ++i) {
        md_vk_page_destroy(dev, MD_VK_VEC_AT(a->pages, md_vk_page_t*, i));
    }
    md_vk_vec_free(&a->pages, dev->alloc);
}

/* =========================================================================
   5. Deferred destruction
   ========================================================================= */

/* Where a stream is right now: the value that will be signalled once all work
   issued into it so far has completed, or 0 if there is nothing to wait for. */
static uint64_t md_vk_stream_position(md_gpu_stream_t s) {
    return s->has_work ? s->next_value : s->submitted_value;
}

/* Queue `object` for destruction once every stream has passed its current
   position. Caller holds device_mutex. If the wait list cannot be allocated,
   falls back to idling the device -- the one case in which this blocks,
   because freeing under the GPU is not an option. */
static void md_vk_retire_locked(md_gpu_device_t dev, md_vk_retire_kind_t kind, void* object) {
    md_vk_retire_t r = {0};
    r.kind   = kind;
    r.object = object;

    uint32_t n = (uint32_t)dev->streams.count;
    if (n > 0) {
        r.waits = (md_vk_wait_t*)md_alloc(dev->alloc, n * sizeof(md_vk_wait_t));
    }
    if (n > 0 && r.waits) {
        md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
        for (uint32_t i = 0; i < n; ++i) {
            uint64_t v = md_vk_stream_position(arr[i]);
            if (v > 0 && md_vk_stream_completed(arr[i]) < v) {
                r.waits[r.wait_count].stream = arr[i];
                r.waits[r.wait_count].value  = v;
                r.wait_count++;
            }
        }
    }

    r.wait_capacity = r.waits ? n : 0;

    md_vk_retire_t* slot = NULL;
    if (n == 0 || r.waits) slot = (md_vk_retire_t*)md_vk_vec_push(&dev->retires, dev->alloc);
    if (!slot) {
        md_vk_fail("out of memory recording a deferred destruction; waiting for the device");
        if (r.waits) md_free(dev->alloc, r.waits, n * sizeof(md_vk_wait_t));
        vkDeviceWaitIdle(dev->device);
        switch (kind) {
        case MD_VK_RETIRE_BLOCK:   md_vk_block_free(dev, (md_vk_block_t*)object);     break;
        case MD_VK_RETIRE_TEXTURE: md_vk_texture_free(dev, (md_gpu_texture_t)object); break;
        case MD_VK_RETIRE_KERNEL:  md_vk_kernel_free(dev, (md_gpu_kernel_t)object);   break;
        }
        return;
    }
    *slot = r;
}

/* Release every retire entry whose waits have all completed (all of them when
   `force`). Caller holds device_mutex. */
static void md_vk_process_retires_locked(md_gpu_device_t dev, bool force) {
    for (size_t i = 0; i < dev->retires.count;) {
        md_vk_retire_t* r = &MD_VK_VEC_AT(dev->retires, md_vk_retire_t, i);
        bool done = true;
        for (uint32_t w = 0; w < r->wait_count && done && !force; ++w) {
            if (md_vk_stream_completed(r->waits[w].stream) < r->waits[w].value) done = false;
        }
        if (!done) { ++i; continue; }
        md_vk_retire_t e = *r;
        md_vk_vec_remove(&dev->retires, i);
        switch (e.kind) {
        case MD_VK_RETIRE_BLOCK:   md_vk_block_free(dev, (md_vk_block_t*)e.object);     break;
        case MD_VK_RETIRE_TEXTURE: md_vk_texture_free(dev, (md_gpu_texture_t)e.object); break;
        case MD_VK_RETIRE_KERNEL:  md_vk_kernel_free(dev, (md_gpu_kernel_t)e.object);   break;
        }
        if (e.waits) md_free(dev->alloc, e.waits, e.wait_capacity * sizeof(md_vk_wait_t));
    }
}

/* A stream is going away after having been synchronised: every wait on it is
   satisfied, so drop the references rather than leave them dangling. Also
   makes blocks it freed reusable outright -- their free_value would otherwise
   never be reached. Caller holds device_mutex. */
static void md_vk_forget_stream_locked(md_gpu_device_t dev, md_gpu_stream_t s) {
    for (size_t i = 0; i < dev->pools.count; ++i) {
        md_gpu_pool_t pool = MD_VK_VEC_AT(dev->pools, md_gpu_pool_t, i);
        for (size_t j = 0; j < pool->blocks.count; ++j) {
            md_vk_block_t* b = MD_VK_VEC_AT(pool->blocks, md_vk_block_t*, j);
            if (b->free_stream == s) { b->free_stream = NULL; b->free_value = 0; }
        }
    }
    /* Clearing the sync makes a pending callback unconditionally ready, so it
       still fires on the next md_gpu_device_poll, on the polling thread. */
    for (size_t i = 0; i < dev->hostfns.count; ++i) {
        md_vk_hostfn_t* h = &MD_VK_VEC_AT(dev->hostfns, md_vk_hostfn_t, i);
        if (h->sync.stream == s) h->sync = md_gpu_sync_none();
    }
    for (size_t i = 0; i < dev->retires.count; ++i) {
        md_vk_retire_t* r = &MD_VK_VEC_AT(dev->retires, md_vk_retire_t, i);
        for (uint32_t w = 0; w < r->wait_count;) {
            if (r->waits[w].stream == s) r->waits[w] = r->waits[--r->wait_count];
            else ++w;
        }
    }
    for (size_t i = 0; i < dev->streams.count; ++i) {
        md_gpu_stream_t o = MD_VK_VEC_AT(dev->streams, md_gpu_stream_t, i);
        for (size_t w = 0; w < o->waits.count;) {
            if (MD_VK_VEC_AT(o->waits, md_vk_wait_t, w).stream == s) md_vk_vec_remove(&o->waits, w);
            else ++w;
        }
    }
}


/* =========================================================================
   6. Device creation
   ========================================================================= */

static VKAPI_ATTR VkBool32 VKAPI_CALL md_vk_debug_cb(
    VkDebugUtilsMessageSeverityFlagBitsEXT severity,
    VkDebugUtilsMessageTypeFlagsEXT types,
    const VkDebugUtilsMessengerCallbackDataEXT* data,
    void* user)
{
    (void)types; (void)user;
    if (severity & VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT) {
        MD_LOG_ERROR("md_gpu validation: %s", data->pMessage);
    } else if (severity & VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT) {
        MD_LOG_DEBUG("md_gpu validation: %s", data->pMessage);
    }
    return VK_FALSE;
}

static bool md_vk_layer_available(const char* name) {
    uint32_t n = 0;
    vkEnumerateInstanceLayerProperties(&n, NULL);
    if (n == 0) return false;
    VkLayerProperties* props = (VkLayerProperties*)malloc(n * sizeof(VkLayerProperties));
    if (!props) return false;
    vkEnumerateInstanceLayerProperties(&n, props);
    bool found = false;
    for (uint32_t i = 0; i < n && !found; ++i) {
        if (strcmp(props[i].layerName, name) == 0) found = true;
    }
    free(props);
    return found;
}

static bool md_vk_ext_available(const char* name) {
    uint32_t n = 0;
    vkEnumerateInstanceExtensionProperties(NULL, &n, NULL);
    if (n == 0) return false;
    VkExtensionProperties* props = (VkExtensionProperties*)malloc(n * sizeof(VkExtensionProperties));
    if (!props) return false;
    vkEnumerateInstanceExtensionProperties(NULL, &n, props);
    bool found = false;
    for (uint32_t i = 0; i < n && !found; ++i) {
        if (strcmp(props[i].extensionName, name) == 0) found = true;
    }
    free(props);
    return found;
}


/* On failure *out_missing names the first unmet requirement. It is always a
   string literal, so the caller may keep it. */
static bool md_vk_probe_device(VkPhysicalDevice pd, struct md_allocator_i* alloc,
                               md_vk_dev_caps_t* out_caps, const char** out_missing)
{
    memset(out_caps, 0, sizeof(*out_caps));
    *out_missing = NULL;

    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(pd, &props);
    /* We call the core 1.3 synchronization2 and dynamic-rendering entry points
       directly rather than their KHR aliases. */
    if (props.apiVersion < VK_API_VERSION_1_3) { *out_missing = "Vulkan 1.3"; return false; }

    /* No device extensions are required: everything md_gpu uses is core 1.3. */
    (void)alloc;

    VkPhysicalDeviceVulkan13Features f13 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES};
    VkPhysicalDeviceVulkan12Features f12 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES, &f13};
    VkPhysicalDeviceFeatures2        f2  = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, &f12};
    vkGetPhysicalDeviceFeatures2(pd, &f2);

#define MD_VK_REQUIRE(cond, name) do { if (!(cond)) { *out_missing = (name); return false; } } while (0)
    MD_VK_REQUIRE(f12.bufferDeviceAddress,   "bufferDeviceAddress");
    MD_VK_REQUIRE(f12.timelineSemaphore,     "timelineSemaphore");
    MD_VK_REQUIRE(f12.descriptorIndexing,    "descriptorIndexing");
    MD_VK_REQUIRE(f12.runtimeDescriptorArray,"runtimeDescriptorArray");
    MD_VK_REQUIRE(f12.descriptorBindingPartiallyBound, "descriptorBindingPartiallyBound");
    /* One per binding in the bindless set, which is UPDATE_AFTER_BIND. */
    MD_VK_REQUIRE(f12.descriptorBindingStorageImageUpdateAfterBind, "descriptorBindingStorageImageUpdateAfterBind");
    MD_VK_REQUIRE(f12.descriptorBindingSampledImageUpdateAfterBind, "descriptorBindingSampledImageUpdateAfterBind");
    /* Slang emits scalar layout for the pointer-reached argument struct. */
    MD_VK_REQUIRE(f12.scalarBlockLayout,     "scalarBlockLayout");
    MD_VK_REQUIRE(f13.synchronization2,      "synchronization2");
    /* The heap arrays carry no format qualifier. */
    MD_VK_REQUIRE(f2.features.shaderStorageImageReadWithoutFormat,  "shaderStorageImageReadWithoutFormat");
    MD_VK_REQUIRE(f2.features.shaderStorageImageWriteWithoutFormat, "shaderStorageImageWriteWithoutFormat");
#undef MD_VK_REQUIRE

    out_caps->maintenance4                = f13.maintenance4;
    out_caps->update_unused_while_pending = f12.descriptorBindingUpdateUnusedWhilePending;
    out_caps->nonuniform_storage_image    = f12.shaderStorageImageArrayNonUniformIndexing;
    out_caps->nonuniform_sampled_image    = f12.shaderSampledImageArrayNonUniformIndexing;
    out_caps->dynamic_storage_image       = f2.features.shaderStorageImageArrayDynamicIndexing;
    out_caps->dynamic_sampled_image       = f2.features.shaderSampledImageArrayDynamicIndexing;
    out_caps->shader_int64                = f2.features.shaderInt64;
    return true;
}

static bool md_vk_create_bindless(md_gpu_device_t dev);
static bool md_vk_immediate(md_gpu_device_t dev, void (*record)(VkCommandBuffer, void*), void* user);
static void md_vk_record_to_general(VkCommandBuffer cmd, void* user);
static bool md_vk_create_dummies(md_gpu_device_t dev);
static md_gpu_stream_t md_vk_stream_create_internal(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label, bool is_default);
static bool md_vk_create_builtin_kernels(md_gpu_device_t dev);

md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc) {
    md_vk_has_error = false;

    struct md_allocator_i* alloc = (desc && desc->alloc) ? desc->alloc : md_get_heap_allocator();

    if (volkInitialize() != VK_SUCCESS) {
        md_vk_fail("volkInitialize failed — no Vulkan loader present");
        return NULL;
    }

    md_gpu_device_t dev = (md_gpu_device_t)md_alloc(alloc, sizeof(md_gpu_device));
    if (!dev) { md_vk_fail("out of memory"); return NULL; }
    memset(dev, 0, sizeof(*dev));
    dev->alloc = alloc;
    dev->validation = desc && desc->enable_validation;

    md_vk_vec_init(&dev->registry, sizeof(md_vk_block_t*));
    md_vk_vec_init(&dev->pools,    sizeof(md_gpu_pool_t));
    md_vk_vec_init(&dev->kernels,  sizeof(md_gpu_kernel_t));
    md_vk_vec_init(&dev->streams,  sizeof(md_gpu_stream_t));
    md_vk_vec_init(&dev->hostfns,  sizeof(md_vk_hostfn_t));
    md_vk_vec_init(&dev->retires,  sizeof(md_vk_retire_t));

    /* ---- instance ---- */
    VkApplicationInfo ai = {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = (desc && desc->label) ? desc->label : "mdlib";
    ai.apiVersion       = VK_API_VERSION_1_3;

    const char* layers[4];     uint32_t layer_count = 0;
    const char* exts[4];       uint32_t ext_count   = 0;

    bool want_debug = dev->validation
        && md_vk_layer_available("VK_LAYER_KHRONOS_validation")
        && md_vk_ext_available(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    if (want_debug) {
        layers[layer_count++] = "VK_LAYER_KHRONOS_validation";
        exts[ext_count++]     = VK_EXT_DEBUG_UTILS_EXTENSION_NAME;
    } else if (dev->validation) {
        MD_LOG_DEBUG("md_gpu: validation requested but unavailable");
    }

    VkInstanceCreateInfo ici = {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ici.pApplicationInfo        = &ai;
    ici.enabledLayerCount       = layer_count;
    ici.ppEnabledLayerNames     = layers;
    ici.enabledExtensionCount   = ext_count;
    ici.ppEnabledExtensionNames = exts;

    if (vkCreateInstance(&ici, NULL, &dev->instance) != VK_SUCCESS) {
        /* Retry without validation. */
        ici.enabledLayerCount = 0;
        ici.enabledExtensionCount = 0;
        want_debug = false;
        if (!md_vk_check(vkCreateInstance(&ici, NULL, &dev->instance), "vkCreateInstance")) {
            md_free(alloc, dev, sizeof(*dev));
            return NULL;
        }
    }
    volkLoadInstance(dev->instance);

    if (want_debug && vkCreateDebugUtilsMessengerEXT) {
        VkDebugUtilsMessengerCreateInfoEXT dci = {VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT};
        dci.messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT;
        dci.messageType     = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT
                            | VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT
                            | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        dci.pfnUserCallback = md_vk_debug_cb;
        vkCreateDebugUtilsMessengerEXT(dev->instance, &dci, NULL, &dev->messenger);
    }

    /* ---- physical device ---- */
    uint32_t pd_count = 0;
    vkEnumeratePhysicalDevices(dev->instance, &pd_count, NULL);
    if (pd_count == 0) {
        md_vk_fail("no Vulkan physical devices");
        goto fail_instance;
    }
    VkPhysicalDevice* pds = (VkPhysicalDevice*)md_alloc(alloc, pd_count * sizeof(VkPhysicalDevice));
    vkEnumeratePhysicalDevices(dev->instance, &pd_count, pds);

    VkPhysicalDevice chosen = VK_NULL_HANDLE;
    md_vk_dev_caps_t caps = {0};
    int best_score = -1;
    /* Kept so that a total failure can name what the best candidate lacked
       rather than saying only that nothing matched. */
    const char* first_missing = NULL;
    char first_missing_dev[256] = {0};

    for (uint32_t i = 0; i < pd_count; ++i) {
        VkPhysicalDeviceProperties p;
        vkGetPhysicalDeviceProperties(pds[i], &p);

        md_vk_dev_caps_t c;
        const char* missing = NULL;
        if (!md_vk_probe_device(pds[i], alloc, &c, &missing)) {
            MD_LOG_DEBUG("md_gpu: skipping '%s' — no %s", p.deviceName, missing ? missing : "?");
            if (!first_missing) {
                first_missing = missing;
                snprintf(first_missing_dev, sizeof(first_missing_dev), "%s", p.deviceName);
            }
            continue;
        }

        int score = 0;
        if (p.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU) score += 100;
        else if (p.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU) score += 50;
        else score += 10;
        if (score > best_score) { best_score = score; chosen = pds[i]; caps = c; }
    }
    md_free(alloc, pds, pd_count * sizeof(VkPhysicalDevice));

    if (!chosen) {
        if (first_missing) {
            md_vk_fail("no usable Vulkan device: '%s' lacks %s", first_missing_dev, first_missing);
        } else {
            md_vk_fail("no usable Vulkan device");
        }
        goto fail_instance;
    }
    dev->phys = chosen;
    vkGetPhysicalDeviceProperties(dev->phys, &dev->props);
    vkGetPhysicalDeviceMemoryProperties(dev->phys, &dev->mem_props);
    dev->is_discrete = dev->props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU;
    {
        VkPhysicalDeviceSubgroupProperties sgp = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
        VkPhysicalDeviceProperties2 p2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2, &sgp};
        vkGetPhysicalDeviceProperties2(dev->phys, &p2);
        dev->subgroup_size = sgp.subgroupSize ? sgp.subgroupSize : 32;
    }

    /* ---- queue families ---- */
    uint32_t qf_count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(dev->phys, &qf_count, NULL);
    VkQueueFamilyProperties* qfs = (VkQueueFamilyProperties*)md_alloc(alloc, qf_count * sizeof(VkQueueFamilyProperties));
    vkGetPhysicalDeviceQueueFamilyProperties(dev->phys, &qf_count, qfs);

    dev->compute_family = UINT32_MAX;
    dev->transfer_family = UINT32_MAX;
    /* Prefer a compute family without graphics (async compute engine). */
    for (uint32_t i = 0; i < qf_count; ++i) {
        if ((qfs[i].queueFlags & VK_QUEUE_COMPUTE_BIT) && !(qfs[i].queueFlags & VK_QUEUE_GRAPHICS_BIT)) {
            dev->compute_family = i; break;
        }
    }
    if (dev->compute_family == UINT32_MAX) {
        for (uint32_t i = 0; i < qf_count; ++i) {
            if (qfs[i].queueFlags & VK_QUEUE_COMPUTE_BIT) { dev->compute_family = i; break; }
        }
    }
    /* Prefer a pure transfer family (DMA engine). */
    for (uint32_t i = 0; i < qf_count; ++i) {
        if ((qfs[i].queueFlags & VK_QUEUE_TRANSFER_BIT) &&
            !(qfs[i].queueFlags & (VK_QUEUE_COMPUTE_BIT | VK_QUEUE_GRAPHICS_BIT))) {
            dev->transfer_family = i; break;
        }
    }
    if (dev->transfer_family == UINT32_MAX) dev->transfer_family = dev->compute_family;
    dev->transfer_can_compute = dev->transfer_family < qf_count &&
                                (qfs[dev->transfer_family].queueFlags & VK_QUEUE_COMPUTE_BIT) != 0;

    if (dev->compute_family != UINT32_MAX) {
        dev->compute_queue_count = qfs[dev->compute_family].queueCount;
        if (dev->compute_queue_count > MD_VK_MAX_QUEUES_PER_FAMILY) dev->compute_queue_count = MD_VK_MAX_QUEUES_PER_FAMILY;
        if (dev->compute_queue_count == 0) dev->compute_queue_count = 1;
        dev->transfer_queue_count = qfs[dev->transfer_family].queueCount;
        if (dev->transfer_queue_count > MD_VK_MAX_QUEUES_PER_FAMILY) dev->transfer_queue_count = MD_VK_MAX_QUEUES_PER_FAMILY;
        if (dev->transfer_queue_count == 0) dev->transfer_queue_count = 1;
    }

    if (dev->compute_family == UINT32_MAX) {
        md_free(alloc, qfs, qf_count * sizeof(VkQueueFamilyProperties));
        md_vk_fail("no compute-capable queue family");
        goto fail_instance;
    }
    md_free(alloc, qfs, qf_count * sizeof(VkQueueFamilyProperties));

    /* ---- logical device ---- */
    static const float prios[MD_VK_MAX_QUEUES_PER_FAMILY] = {1,1,1,1,1,1,1,1};
    VkDeviceQueueCreateInfo qci[2];
    uint32_t qci_count = 0;
    qci[qci_count] = (VkDeviceQueueCreateInfo){VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qci[qci_count].queueFamilyIndex = dev->compute_family;
    qci[qci_count].queueCount       = dev->compute_queue_count;
    qci[qci_count].pQueuePriorities = prios;
    qci_count++;
    if (dev->transfer_family != dev->compute_family) {
        qci[qci_count] = (VkDeviceQueueCreateInfo){VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        qci[qci_count].queueFamilyIndex = dev->transfer_family;
        qci[qci_count].queueCount       = dev->transfer_queue_count;
        qci[qci_count].pQueuePriorities = prios;
        qci_count++;
    }

    /* Everything below was confirmed present by md_vk_probe_device. Enabling a
       feature the driver does not report is VK_ERROR_FEATURE_NOT_PRESENT, so
       the optional ones are gated on the probe's answer. */
    VkPhysicalDeviceVulkan13Features f13 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES};
    f13.synchronization2 = VK_TRUE;
    f13.maintenance4     = caps.maintenance4 ? VK_TRUE : VK_FALSE;

    VkPhysicalDeviceVulkan12Features f12 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    f12.pNext = &f13;
    f12.bufferDeviceAddress                          = VK_TRUE;
    f12.timelineSemaphore                            = VK_TRUE;
    f12.descriptorIndexing                           = VK_TRUE;
    f12.runtimeDescriptorArray                       = VK_TRUE;
    f12.descriptorBindingPartiallyBound              = VK_TRUE;
    f12.descriptorBindingStorageImageUpdateAfterBind = VK_TRUE;
    f12.descriptorBindingSampledImageUpdateAfterBind = VK_TRUE;
    f12.scalarBlockLayout                            = VK_TRUE;
    f12.descriptorBindingUpdateUnusedWhilePending = caps.update_unused_while_pending ? VK_TRUE : VK_FALSE;
    f12.shaderStorageImageArrayNonUniformIndexing = caps.nonuniform_storage_image    ? VK_TRUE : VK_FALSE;
    f12.shaderSampledImageArrayNonUniformIndexing = caps.nonuniform_sampled_image    ? VK_TRUE : VK_FALSE;

    VkPhysicalDeviceFeatures2 f2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    f2.pNext = &f12;
    /* Slang's heap arrays carry no format qualifier. */
    f2.features.shaderStorageImageReadWithoutFormat    = VK_TRUE;
    f2.features.shaderStorageImageWriteWithoutFormat   = VK_TRUE;
    f2.features.shaderStorageImageArrayDynamicIndexing = caps.dynamic_storage_image ? VK_TRUE : VK_FALSE;
    f2.features.shaderSampledImageArrayDynamicIndexing = caps.dynamic_sampled_image ? VK_TRUE : VK_FALSE;
    f2.features.shaderInt64                            = caps.shader_int64          ? VK_TRUE : VK_FALSE;

    VkDeviceCreateInfo dci = {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.pNext                   = &f2;
    dci.queueCreateInfoCount     = qci_count;
    dci.pQueueCreateInfos        = qci;

    if (!md_vk_check(vkCreateDevice(dev->phys, &dci, NULL, &dev->device), "vkCreateDevice")) goto fail_instance;
    dev->caps = caps;
    volkLoadDevice(dev->device);

    for (uint32_t i = 0; i < dev->compute_queue_count; ++i) {
        vkGetDeviceQueue(dev->device, dev->compute_family, i, &dev->compute_queues[i]);
    }
    if (dev->transfer_family != dev->compute_family) {
        for (uint32_t i = 0; i < dev->transfer_queue_count; ++i) {
            vkGetDeviceQueue(dev->device, dev->transfer_family, i, &dev->transfer_queues[i]);
        }
    } else {
        dev->transfer_queue_count = dev->compute_queue_count;
        for (uint32_t i = 0; i < dev->compute_queue_count; ++i) {
            dev->transfer_queues[i] = dev->compute_queues[i];
        }
    }

    md_mutex_init(&dev->queue_mutex);
    md_mutex_init(&dev->device_mutex);

    dev->share_families[0]  = dev->compute_family;
    dev->share_family_count = 1;
    if (dev->transfer_family != dev->compute_family) {
        dev->share_families[1]  = dev->transfer_family;
        dev->share_family_count = 2;
    }

    /* Heap free list, allocated high-to-low so slots are handed out from 1.
       Slot 0 is never handed out so that a zero handle stays the null handle. */
    for (uint32_t i = 0; i < MD_VK_MAX_TEXTURE_SLOTS - 1; ++i) dev->tex_free[i] = MD_VK_MAX_TEXTURE_SLOTS - 1 - i;
    dev->tex_free_count = MD_VK_MAX_TEXTURE_SLOTS - 1;

    if (!md_vk_create_bindless(dev)) goto fail_device;
    if (!md_vk_create_dummies(dev))  goto fail_device;

    dev->default_compute  = md_vk_stream_create_internal(dev, MD_GPU_STREAM_COMPUTE,  "default compute",  true);
    dev->default_transfer = md_vk_stream_create_internal(dev, MD_GPU_STREAM_TRANSFER, "default transfer", true);
    if (!dev->default_compute || !dev->default_transfer) goto fail_device;

    if (!md_vk_create_builtin_kernels(dev)) goto fail_device;

    MD_LOG_DEBUG("md_gpu: device '%s' (compute family %u x%u queues, transfer family %u x%u queues)",
                 dev->props.deviceName, dev->compute_family, dev->compute_queue_count,
                 dev->transfer_family, dev->transfer_queue_count);
    return dev;

fail_device:
    md_gpu_device_destroy(dev);
    return NULL;

fail_instance:
    if (dev->messenger && vkDestroyDebugUtilsMessengerEXT) vkDestroyDebugUtilsMessengerEXT(dev->instance, dev->messenger, NULL);
    if (dev->instance) vkDestroyInstance(dev->instance, NULL);
    md_free(alloc, dev, sizeof(*dev));
    return NULL;
}

bool md_gpu_device_info(md_gpu_device_t dev, md_gpu_device_info_t* info) {
    if (!dev || !info) return false;
    memset(info, 0, sizeof(*info));
    info->is_discrete              = dev->is_discrete;
    info->max_threads_per_group    = dev->props.limits.maxComputeWorkGroupInvocations;
    info->preferred_group_multiple = dev->subgroup_size;
    snprintf(info->name, sizeof(info->name), "%s", dev->props.deviceName);
    return true;
}

/* =========================================================================
   7. Bindless descriptor set
   ========================================================================= */

/* A single 1x1x1 placeholder, valid as both a storage and a sampled image.
   Freed heap slots are pointed at it so that no descriptor ever references a
   destroyed view. It is never accessed by any shader. */
static bool md_vk_create_dummies(md_gpu_device_t dev) {
    VkImageCreateInfo ici = {VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    ici.imageType     = VK_IMAGE_TYPE_3D;
    ici.format        = VK_FORMAT_R32_SFLOAT;
    ici.extent.width  = 1;
    ici.extent.height = 1;
    ici.extent.depth  = 1;
    ici.mipLevels     = 1;
    ici.arrayLayers   = 1;
    ici.samples       = VK_SAMPLE_COUNT_1_BIT;
    ici.tiling        = VK_IMAGE_TILING_OPTIMAL;
    ici.usage         = VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
    ici.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    if (!md_vk_check(vkCreateImage(dev->device, &ici, NULL, &dev->dummy_image), "vkCreateImage (dummy)")) return false;

    VkMemoryRequirements req;
    vkGetImageMemoryRequirements(dev->device, dev->dummy_image, &req);
    uint32_t type = md_vk_find_memory_type(dev, req.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0);
    if (type == UINT32_MAX) type = md_vk_find_memory_type(dev, req.memoryTypeBits, 0, 0);

    VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.allocationSize  = req.size;
    mai.memoryTypeIndex = type;
    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &dev->dummy_mem), "vkAllocateMemory (dummy)")) return false;
    vkBindImageMemory(dev->device, dev->dummy_image, dev->dummy_mem, 0);

    VkImageViewCreateInfo vci = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    vci.image    = dev->dummy_image;
    vci.viewType = VK_IMAGE_VIEW_TYPE_3D;
    vci.format   = VK_FORMAT_R32_SFLOAT;
    vci.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    vci.subresourceRange.levelCount = 1;
    vci.subresourceRange.layerCount = 1;
    if (!md_vk_check(vkCreateImageView(dev->device, &vci, NULL, &dev->dummy_view), "vkCreateImageView (dummy)")) return false;

    if (!md_vk_immediate(dev, md_vk_record_to_general, &dev->dummy_image)) return false;

    VkSamplerCreateInfo sci = {VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    sci.maxLod = VK_LOD_CLAMP_NONE;
    return md_vk_check(vkCreateSampler(dev->device, &sci, NULL, &dev->dummy_sampler), "vkCreateSampler (dummy)");
}

/* The one place descriptor type maps to binding. */
static uint32_t md_vk_binding_for(VkDescriptorType type) {
    switch (type) {
    case VK_DESCRIPTOR_TYPE_SAMPLER:       return MD_VK_BINDING_SAMPLER;
    case VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE: return MD_VK_BINDING_SAMPLED_IMAGE;
    case VK_DESCRIPTOR_TYPE_STORAGE_IMAGE: return MD_VK_BINDING_STORAGE_IMAGE;
    default:                               return UINT32_MAX;
    }
}

/* Point a freed heap slot at the placeholder. */
static void md_vk_clear_slot(md_gpu_device_t dev, uint32_t index, VkDescriptorType type) {
    VkDescriptorImageInfo dii = {VK_NULL_HANDLE, VK_NULL_HANDLE, VK_IMAGE_LAYOUT_GENERAL};
    uint32_t binding = md_vk_binding_for(type);
    if (binding == UINT32_MAX) return;
    if (type == VK_DESCRIPTOR_TYPE_SAMPLER) {
        dii.sampler     = dev->dummy_sampler;
        dii.imageLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    } else {
        dii.imageView = dev->dummy_view;
    }
    VkWriteDescriptorSet w = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    w.dstSet          = dev->desc_set;
    w.dstBinding      = binding;
    w.dstArrayElement = index;
    w.descriptorCount = 1;
    w.descriptorType  = type;
    w.pImageInfo      = &dii;
    vkUpdateDescriptorSets(dev->device, 1, &w, 0, NULL);
}

static bool md_vk_create_bindless(md_gpu_device_t dev) {
    /* One binding per descriptor type, at the indices Slang's `None` bindless
       preset uses. Declared bindings need not be contiguous, so binding 1
       (combined image sampler) is simply absent -- md_gpu passes textures and
       samplers separately and never populates it. */
    VkDescriptorSetLayoutBinding bindings[MD_VK_BINDING_COUNT] = {0};
    VkDescriptorBindingFlags     bflags[MD_VK_BINDING_COUNT];

    bindings[0].binding         = MD_VK_BINDING_SAMPLER;
    bindings[0].descriptorType  = VK_DESCRIPTOR_TYPE_SAMPLER;
    bindings[0].descriptorCount = MD_VK_MAX_SAMPLERS;
    bindings[1].binding         = MD_VK_BINDING_SAMPLED_IMAGE;
    bindings[1].descriptorType  = VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE;
    bindings[1].descriptorCount = MD_VK_MAX_TEXTURE_SLOTS;
    bindings[2].binding         = MD_VK_BINDING_STORAGE_IMAGE;
    bindings[2].descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings[2].descriptorCount = MD_VK_MAX_TEXTURE_SLOTS;

    for (uint32_t i = 0; i < MD_VK_BINDING_COUNT; ++i) {
        bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bflags[i] = VK_DESCRIPTOR_BINDING_PARTIALLY_BOUND_BIT
                  | VK_DESCRIPTOR_BINDING_UPDATE_AFTER_BIND_BIT;
        if (dev->caps.update_unused_while_pending) {
            bflags[i] |= VK_DESCRIPTOR_BINDING_UPDATE_UNUSED_WHILE_PENDING_BIT;
        }
    }

    VkDescriptorSetLayoutBindingFlagsCreateInfo bfci = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_BINDING_FLAGS_CREATE_INFO};
    bfci.bindingCount  = MD_VK_BINDING_COUNT;
    bfci.pBindingFlags = bflags;

    VkDescriptorSetLayoutCreateInfo lci = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    lci.pNext        = &bfci;
    lci.bindingCount = MD_VK_BINDING_COUNT;
    lci.pBindings    = bindings;
    lci.flags        = VK_DESCRIPTOR_SET_LAYOUT_CREATE_UPDATE_AFTER_BIND_POOL_BIT;

    if (!md_vk_check(vkCreateDescriptorSetLayout(dev->device, &lci, NULL, &dev->set_layout), "vkCreateDescriptorSetLayout")) return false;

    VkDescriptorPoolSize sizes[3];
    sizes[0] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_SAMPLER,       MD_VK_MAX_SAMPLERS};
    sizes[1] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, MD_VK_MAX_TEXTURE_SLOTS};
    sizes[2] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, MD_VK_MAX_TEXTURE_SLOTS};

    VkDescriptorPoolCreateInfo pci = {VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    pci.flags         = VK_DESCRIPTOR_POOL_CREATE_UPDATE_AFTER_BIND_BIT;
    pci.maxSets       = 1;
    pci.poolSizeCount = 3;
    pci.pPoolSizes    = sizes;
    if (!md_vk_check(vkCreateDescriptorPool(dev->device, &pci, NULL, &dev->desc_pool), "vkCreateDescriptorPool")) return false;

    VkDescriptorSetAllocateInfo dai = {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    dai.descriptorPool     = dev->desc_pool;
    dai.descriptorSetCount = 1;
    dai.pSetLayouts        = &dev->set_layout;
    if (!md_vk_check(vkAllocateDescriptorSets(dev->device, &dai, &dev->desc_set), "vkAllocateDescriptorSets")) return false;

    VkPushConstantRange pcr = {VK_SHADER_STAGE_COMPUTE_BIT, 0, 8};
    VkPipelineLayoutCreateInfo plci = {VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    plci.setLayoutCount         = 1;
    plci.pSetLayouts            = &dev->set_layout;
    plci.pushConstantRangeCount = 1;
    plci.pPushConstantRanges    = &pcr;
    if (!md_vk_check(vkCreatePipelineLayout(dev->device, &plci, NULL, &dev->pipeline_layout), "vkCreatePipelineLayout")) return false;

    return true;
}


/* One-shot command submission with a CPU wait. Used only during device
   creation, for the placeholder image, when nothing can be in flight. Every
   other layout transition is recorded into a stream. */
static bool md_vk_immediate(md_gpu_device_t dev, void (*record)(VkCommandBuffer, void*), void* user) {
    VkCommandPoolCreateInfo cpci = {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    cpci.flags            = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT;
    cpci.queueFamilyIndex = dev->compute_family;
    VkCommandPool pool;
    if (!md_vk_check(vkCreateCommandPool(dev->device, &cpci, NULL, &pool), "vkCreateCommandPool")) return false;

    VkCommandBufferAllocateInfo cbai = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cbai.commandPool        = pool;
    cbai.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cbai.commandBufferCount = 1;
    VkCommandBuffer cmd;
    if (!md_vk_check(vkAllocateCommandBuffers(dev->device, &cbai, &cmd), "vkAllocateCommandBuffers")) {
        vkDestroyCommandPool(dev->device, pool, NULL);
        return false;
    }

    VkCommandBufferBeginInfo bi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(cmd, &bi);
    record(cmd, user);
    vkEndCommandBuffer(cmd);

    VkCommandBufferSubmitInfo cbsi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO};
    cbsi.commandBuffer = cmd;
    VkSubmitInfo2 si = {VK_STRUCTURE_TYPE_SUBMIT_INFO_2};
    si.commandBufferInfoCount = 1;
    si.pCommandBufferInfos    = &cbsi;

    VkFenceCreateInfo fci = {VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    VkFence fence;
    vkCreateFence(dev->device, &fci, NULL, &fence);

    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueueSubmit2(dev->compute_queues[0], 1, &si, fence);
    md_mutex_unlock(&dev->queue_mutex);

    if (r == VK_SUCCESS) vkWaitForFences(dev->device, 1, &fence, VK_TRUE, UINT64_MAX);
    vkDestroyFence(dev->device, fence, NULL);
    vkDestroyCommandPool(dev->device, pool, NULL);
    return md_vk_check(r, "vkQueueSubmit2 (immediate)");
}

static void md_vk_record_to_general(VkCommandBuffer cmd, void* user) {
    VkImage image = *(VkImage*)user;
    VkImageMemoryBarrier2 b = {VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2};
    b.srcStageMask  = VK_PIPELINE_STAGE_2_TOP_OF_PIPE_BIT;
    b.srcAccessMask = 0;
    b.dstStageMask  = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    b.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
    b.oldLayout     = VK_IMAGE_LAYOUT_UNDEFINED;
    b.newLayout     = VK_IMAGE_LAYOUT_GENERAL;
    b.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.image = image;
    b.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    b.subresourceRange.levelCount = 1;
    b.subresourceRange.layerCount = 1;

    VkDependencyInfo di = {VK_STRUCTURE_TYPE_DEPENDENCY_INFO};
    di.imageMemoryBarrierCount = 1;
    di.pImageMemoryBarriers    = &b;
    vkCmdPipelineBarrier2(cmd, &di);
}

/* Take a slot from the heap free list. Slot 0 is reserved as the null handle. */
static uint32_t md_vk_alloc_slot(md_gpu_device_t dev) {
    if (dev->tex_free_count == 0) return UINT32_MAX;
    return dev->tex_free[--dev->tex_free_count];
}

static void md_vk_write_image_slot(md_gpu_device_t dev, uint32_t slot, VkImageView view, VkDescriptorType type) {
    uint32_t binding = md_vk_binding_for(type);
    if (binding == UINT32_MAX) return;
    VkDescriptorImageInfo dii = {VK_NULL_HANDLE, view, VK_IMAGE_LAYOUT_GENERAL};
    VkWriteDescriptorSet w = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    w.dstSet          = dev->desc_set;
    w.dstBinding      = binding;
    w.dstArrayElement = slot;
    w.descriptorCount = 1;
    w.descriptorType  = type;
    w.pImageInfo      = &dii;
    vkUpdateDescriptorSets(dev->device, 1, &w, 0, NULL);
}


/* =========================================================================
   8. Streams
   ========================================================================= */

static uint64_t md_vk_stream_completed(md_gpu_stream_t s) {
    uint64_t v = 0;
    vkGetSemaphoreCounterValue(s->device->device, s->timeline, &v);
    return v;
}

static md_gpu_stream_t md_vk_stream_create_internal(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label, bool is_default) {
    md_gpu_stream_t s = (md_gpu_stream_t)md_alloc(dev->alloc, sizeof(md_gpu_stream));
    if (!s) { md_vk_fail("out of memory"); return NULL; }
    memset(s, 0, sizeof(*s));
    s->device     = dev;
    s->kind       = kind;
    s->ordering   = MD_GPU_ORDER_IMPLICIT;
    s->next_value = 1;
    s->is_default = is_default;
    snprintf(s->label, sizeof(s->label), "%s", label ? label : "stream");

    md_mutex_lock(&dev->device_mutex);
    if (kind == MD_GPU_STREAM_TRANSFER) {
        s->family      = dev->transfer_family;
        s->can_compute = dev->transfer_can_compute;
        s->queue       = dev->transfer_queues[dev->next_transfer_queue % dev->transfer_queue_count];
        dev->next_transfer_queue++;
    } else {
        s->family      = dev->compute_family;
        s->can_compute = true;
        s->queue       = dev->compute_queues[dev->next_compute_queue % dev->compute_queue_count];
        dev->next_compute_queue++;
    }
    md_mutex_unlock(&dev->device_mutex);

    md_vk_vec_init(&s->cmds, sizeof(md_vk_cmd_t));
    md_vk_vec_init(&s->waits, sizeof(md_vk_wait_t));
    md_vk_vec_init(&s->arena.pages, sizeof(md_vk_page_t*));

    VkSemaphoreTypeCreateInfo stci = {VK_STRUCTURE_TYPE_SEMAPHORE_TYPE_CREATE_INFO};
    stci.semaphoreType = VK_SEMAPHORE_TYPE_TIMELINE;
    stci.initialValue  = 0;
    VkSemaphoreCreateInfo sci = {VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    sci.pNext = &stci;
    if (!md_vk_check(vkCreateSemaphore(dev->device, &sci, NULL, &s->timeline), "vkCreateSemaphore")) {
        md_free(dev->alloc, s, sizeof(*s));
        return NULL;
    }

    VkCommandPoolCreateInfo cpci = {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    cpci.flags            = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    cpci.queueFamilyIndex = s->family;
    if (!md_vk_check(vkCreateCommandPool(dev->device, &cpci, NULL, &s->cmd_pool), "vkCreateCommandPool")) {
        vkDestroySemaphore(dev->device, s->timeline, NULL);
        md_free(dev->alloc, s, sizeof(*s));
        return NULL;
    }

    md_mutex_lock(&dev->device_mutex);
    md_gpu_stream_t* slot = (md_gpu_stream_t*)md_vk_vec_push(&dev->streams, dev->alloc);
    if (slot) *slot = s;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) {
        vkDestroyCommandPool(dev->device, s->cmd_pool, NULL);
        vkDestroySemaphore(dev->device, s->timeline, NULL);
        md_free(dev->alloc, s, sizeof(*s));
        md_vk_fail("out of memory");
        return NULL;
    }
    return s;
}

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label) {
    if (!dev) { md_vk_fail("md_gpu_stream_create: null device"); return NULL; }
    return md_vk_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    return kind == MD_GPU_STREAM_TRANSFER ? dev->default_transfer : dev->default_compute;
}

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t s) { return s ? s->device : NULL; }

/* Acquire a command buffer, recycling any whose submission has completed. */
static bool md_vk_stream_ensure_cmd(md_gpu_stream_t s) {
    if (s->open) return true;
    md_gpu_device_t dev = s->device;

    uint64_t done = md_vk_stream_completed(s);
    md_vk_cmd_t* found = NULL;
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &MD_VK_VEC_AT(s->cmds, md_vk_cmd_t, i);
        if (c->pending && c->value <= done) c->pending = false;
        if (!c->pending && !found) found = c;
    }

    if (!found) {
        VkCommandBufferAllocateInfo cbai = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
        cbai.commandPool        = s->cmd_pool;
        cbai.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
        cbai.commandBufferCount = 1;
        VkCommandBuffer cb = VK_NULL_HANDLE;
        if (!md_vk_check(vkAllocateCommandBuffers(dev->device, &cbai, &cb), "vkAllocateCommandBuffers")) return false;
        md_vk_cmd_t* slot = (md_vk_cmd_t*)md_vk_vec_push(&s->cmds, dev->alloc);
        if (!slot) { vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &cb); return md_vk_fail("out of memory"); }
        slot->cmd = cb;
        found = slot;
    }

    vkResetCommandBuffer(found->cmd, 0);
    VkCommandBufferBeginInfo bi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    if (!md_vk_check(vkBeginCommandBuffer(found->cmd, &bi), "vkBeginCommandBuffer")) return false;

    s->open     = found->cmd;
    s->has_work = false;
    /* If this stream has already submitted, the first operation recorded here
       still needs a dependency on that earlier work. Batches submitted to one
       queue may overlap -- submission order alone creates no memory
       dependency. A barrier's first synchronisation scope covers everything
       earlier in submission order, which spans command buffers and batches, so
       one barrier at the top of the buffer closes it. (In EXPLICIT mode the
       caller's own barriers carry the same scope across submissions.) */
    if (s->submitted_value > 0) s->needs_barrier = true;
    return true;
}

/* The full-strength masks one queue family can express. A transfer-only
   family rejects compute and indirect stages, so they must be dropped there,
   not merely left unused. */
static void md_vk_stream_full_masks(md_gpu_stream_t s, VkPipelineStageFlags2* src_stage, VkAccessFlags2* src_access,
                                    VkPipelineStageFlags2* dst_stage, VkAccessFlags2* dst_access) {
    *src_stage  = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
    *src_access = VK_ACCESS_2_TRANSFER_WRITE_BIT;
    *dst_stage  = VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
    *dst_access = VK_ACCESS_2_TRANSFER_READ_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT;
    if (s->can_compute) {
        *src_stage  |= VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT;
        *src_access |= VK_ACCESS_2_SHADER_WRITE_BIT;
        *dst_stage  |= VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_2_DRAW_INDIRECT_BIT;
        *dst_access |= VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT
                     | VK_ACCESS_2_INDIRECT_COMMAND_READ_BIT;
    }
}

static void md_vk_emit_barrier(md_gpu_stream_t s, VkPipelineStageFlags2 src_stage, VkAccessFlags2 src_access,
                               VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    VkMemoryBarrier2 mb = {VK_STRUCTURE_TYPE_MEMORY_BARRIER_2};
    mb.srcStageMask  = src_stage;
    mb.srcAccessMask = src_access;
    mb.dstStageMask  = dst_stage;
    mb.dstAccessMask = dst_access;
    VkDependencyInfo di = {VK_STRUCTURE_TYPE_DEPENDENCY_INFO};
    di.memoryBarrierCount = 1;
    di.pMemoryBarriers    = &mb;
    vkCmdPipelineBarrier2(s->open, &di);
}

/* The command buffer an operation records into, with IMPLICIT ordering
   applied: one global barrier before the operation if anything precedes it.
   Returns VK_NULL_HANDLE on failure. */
static VkCommandBuffer md_vk_begin_op(md_gpu_stream_t s) {
    if (s->upload_open) { md_vk_fail("stream '%s' has an open upload; call md_gpu_upload_end first", s->label); return VK_NULL_HANDLE; }
    if (!md_vk_stream_ensure_cmd(s)) return VK_NULL_HANDLE;
    if (s->needs_barrier && s->ordering == MD_GPU_ORDER_IMPLICIT) {
        VkPipelineStageFlags2 ss, ds; VkAccessFlags2 sa, da;
        md_vk_stream_full_masks(s, &ss, &sa, &ds, &da);
        md_vk_emit_barrier(s, ss, sa, ds, da);
        s->needs_barrier = false;
    }
    return s->open;
}

/* Mark that an operation was recorded. */
static void md_vk_end_op(md_gpu_stream_t s) {
    s->needs_barrier = true;
    s->has_work      = true;
}

void md_gpu_stream_set_ordering(md_gpu_stream_t s, md_gpu_ordering_t ordering) {
    if (!s) return;
    if (ordering == MD_GPU_ORDER_IMPLICIT && s->ordering != MD_GPU_ORDER_IMPLICIT) {
        /* Order the next operation after everything in the explicit region. */
        s->needs_barrier = true;
    }
    s->ordering = ordering;
}

md_gpu_ordering_t md_gpu_stream_ordering(md_gpu_stream_t s) {
    return s ? s->ordering : MD_GPU_ORDER_IMPLICIT;
}

static void md_vk_stage_masks(md_gpu_stream_t s, md_gpu_stage_flags_t stages, bool producer,
                              VkPipelineStageFlags2* out_stage, VkAccessFlags2* out_access) {
    VkPipelineStageFlags2 st = 0;
    VkAccessFlags2        ac = 0;
    if (stages == MD_GPU_STAGE_ALL) {
        VkPipelineStageFlags2 ss, ds; VkAccessFlags2 sa, da;
        md_vk_stream_full_masks(s, &ss, &sa, &ds, &da);
        *out_stage  = producer ? ss : ds;
        *out_access = producer ? sa : da;
        return;
    }
    if (stages & MD_GPU_STAGE_TRANSFER) {
        st |= VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
        ac |= producer ? VK_ACCESS_2_TRANSFER_WRITE_BIT
                       : (VK_ACCESS_2_TRANSFER_READ_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT);
    }
    if ((stages & MD_GPU_STAGE_COMPUTE) && s->can_compute) {
        st |= VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT;
        ac |= producer ? VK_ACCESS_2_SHADER_WRITE_BIT
                       : (VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT);
    }
    if ((stages & MD_GPU_STAGE_INDIRECT) && s->can_compute && !producer) {
        st |= VK_PIPELINE_STAGE_2_DRAW_INDIRECT_BIT;
        ac |= VK_ACCESS_2_INDIRECT_COMMAND_READ_BIT;
    }
    *out_stage  = st;
    *out_access = ac;
}

void md_gpu_barrier(md_gpu_stream_t s, md_gpu_stage_flags_t producers, md_gpu_stage_flags_t consumers) {
    if (!s) return;
    if (s->upload_open) { md_vk_fail("md_gpu_barrier: stream '%s' has an open upload", s->label); return; }
    VkPipelineStageFlags2 ss, ds; VkAccessFlags2 sa, da;
    md_vk_stage_masks(s, producers, true,  &ss, &sa);
    md_vk_stage_masks(s, consumers, false, &ds, &da);
    if (!ss || !ds) return;
    if (!md_vk_stream_ensure_cmd(s)) return;
    md_vk_emit_barrier(s, ss, sa, ds, da);
    /* A barrier is not an operation: it does not by itself need ordering, but
       it must be submitted along with whatever follows it. */
}

static bool md_vk_stream_submit(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    if (!s->open || !s->has_work) {
        /* Nothing recorded. Pending waits stay queued for the next submit. */
        return true;
    }

    if (!md_vk_check(vkEndCommandBuffer(s->open), "vkEndCommandBuffer")) return false;

    uint64_t signal_value = s->next_value;

    VkCommandBufferSubmitInfo cbsi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO};
    cbsi.commandBuffer = s->open;

    VkSemaphoreSubmitInfo  wait_stack[8];
    VkSemaphoreSubmitInfo* waits = wait_stack;
    const uint32_t wait_count = (uint32_t)s->waits.count;
    if (wait_count > 8) {
        waits = (VkSemaphoreSubmitInfo*)md_alloc(dev->alloc, wait_count * sizeof(VkSemaphoreSubmitInfo));
        if (!waits) return md_vk_fail("out of memory");
    }
    for (uint32_t i = 0; i < wait_count; ++i) {
        md_vk_wait_t w = MD_VK_VEC_AT(s->waits, md_vk_wait_t, i);
        waits[i] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        waits[i].semaphore = w.stream->timeline;
        waits[i].value     = w.value;
        waits[i].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }

    VkSemaphoreSubmitInfo signal = {VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
    signal.semaphore = s->timeline;
    signal.value     = signal_value;
    signal.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;

    VkSubmitInfo2 si = {VK_STRUCTURE_TYPE_SUBMIT_INFO_2};
    si.waitSemaphoreInfoCount   = wait_count;
    si.pWaitSemaphoreInfos      = waits;
    si.commandBufferInfoCount   = 1;
    si.pCommandBufferInfos      = &cbsi;
    si.signalSemaphoreInfoCount = 1;
    si.pSignalSemaphoreInfos    = &signal;

    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueueSubmit2(s->queue, 1, &si, VK_NULL_HANDLE);
    md_mutex_unlock(&dev->queue_mutex);
    if (waits != wait_stack) md_free(dev->alloc, waits, wait_count * sizeof(VkSemaphoreSubmitInfo));
    if (!md_vk_check(r, "vkQueueSubmit2")) return false;

    /* Tag the command buffer and arena pages with the value that frees them. */
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &MD_VK_VEC_AT(s->cmds, md_vk_cmd_t, i);
        if (c->cmd == s->open) { c->value = signal_value; c->pending = true; break; }
    }
    md_vk_arena_retire(s, signal_value);

    s->submitted_value = signal_value;
    s->next_value      = signal_value + 1;
    s->waits.count     = 0;
    s->open            = VK_NULL_HANDLE;
    s->has_work        = false;
    /* md_vk_stream_ensure_cmd re-arms needs_barrier for the next buffer. */
    s->needs_barrier   = false;
    return true;
}

md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t s) {
    md_gpu_sync_t out = md_gpu_sync_none();
    if (!s) return out;
    if (s->upload_open) { md_vk_fail("md_gpu_stream_record: stream '%s' has an open upload", s->label); return out; }
    md_vk_stream_submit(s);
    if (s->submitted_value == 0) return out;
    out.stream = s;
    out.value  = s->submitted_value;
    return out;
}

void md_gpu_stream_wait(md_gpu_stream_t s, md_gpu_sync_t sync) {
    if (!s || !md_gpu_sync_is_valid(sync)) return;
    if (sync.stream == s) return;                    /* already ordered */
    if (sync.stream->device != s->device) { md_vk_fail("md_gpu_stream_wait: sync from another device"); return; }
    if (md_gpu_sync_is_complete(sync)) return;       /* nothing to wait for */

    /* Work already issued must not be retroactively delayed: close it first. */
    if (s->has_work) md_vk_stream_submit(s);

    /* Collapse duplicate waits on the same timeline, keeping the larger value. */
    for (size_t i = 0; i < s->waits.count; ++i) {
        md_vk_wait_t* w = &MD_VK_VEC_AT(s->waits, md_vk_wait_t, i);
        if (w->stream == sync.stream) {
            if (sync.value > w->value) w->value = sync.value;
            return;
        }
    }
    md_vk_wait_t* w = (md_vk_wait_t*)md_vk_vec_push(&s->waits, s->device->alloc);
    if (!w) { md_vk_fail("out of memory queueing a stream wait"); return; }
    w->stream = sync.stream;
    w->value  = sync.value;
}

void md_gpu_stream_flush(md_gpu_stream_t s) {
    if (!s) return;
    md_vk_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    md_vk_stream_submit(s);
    if (s->submitted_value > 0) {
        VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
        wi.semaphoreCount = 1;
        wi.pSemaphores    = &s->timeline;
        wi.pValues        = &s->submitted_value;
        vkWaitSemaphores(s->device->device, &wi, UINT64_MAX);
    }
}

bool md_gpu_sync_is_complete(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    return md_vk_stream_completed(sync.stream) >= sync.value;
}

void md_gpu_sync_wait(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return;
    VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
    wi.semaphoreCount = 1;
    wi.pSemaphores    = &sync.stream->timeline;
    wi.pValues        = &sync.value;
    vkWaitSemaphores(sync.stream->device->device, &wi, UINT64_MAX);
}

/* Free a stream's own objects. The stream must be idle. */
static void md_vk_stream_free(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    if (s->open) {
        vkEndCommandBuffer(s->open);
        s->open = VK_NULL_HANDLE;
    }
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &MD_VK_VEC_AT(s->cmds, md_vk_cmd_t, i);
        vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &c->cmd);
    }
    md_vk_vec_free(&s->cmds, dev->alloc);
    md_vk_vec_free(&s->waits, dev->alloc);
    md_vk_arena_free(dev, &s->arena);
    vkDestroyCommandPool(dev->device, s->cmd_pool, NULL);
    vkDestroySemaphore(dev->device, s->timeline, NULL);
    md_free(dev->alloc, s, sizeof(*s));
}

void md_gpu_stream_destroy(md_gpu_stream_t s) {
    if (!s || s->is_default) return;
    md_gpu_device_t dev = s->device;

    /* Only this stream's own work. Anything else waiting on it is thereby
       satisfied, which is what makes forgetting it below safe. */
    s->upload_open = false;
    md_gpu_stream_sync(s);

    md_mutex_lock(&dev->device_mutex);
    md_vk_vec_remove_ptr(&dev->streams, s);
    md_vk_forget_stream_locked(dev, s);
    md_mutex_unlock(&dev->device_mutex);

    md_vk_stream_free(s);
}

/* =========================================================================
   9. Memory
   ========================================================================= */

md_gpu_pool_t md_gpu_pool_create(md_gpu_device_t dev, const md_gpu_pool_desc_t* desc) {
    if (!dev || !desc) { md_vk_fail("md_gpu_pool_create: null argument"); return NULL; }
    if ((unsigned)desc->kind > (unsigned)MD_GPU_MEM_HOST_READ) { md_vk_fail("md_gpu_pool_create: invalid memory kind %d", (int)desc->kind); return NULL; }
    md_gpu_pool_t p = (md_gpu_pool_t)md_alloc(dev->alloc, sizeof(md_gpu_pool));
    if (!p) { md_vk_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device      = dev;
    p->kind        = desc->kind;
    p->cache_limit = desc->cache_limit;
    snprintf(p->label, sizeof(p->label), "%s", desc->label ? desc->label : "pool");
    md_vk_vec_init(&p->blocks,   sizeof(md_vk_block_t*));
    md_vk_vec_init(&p->textures, sizeof(md_gpu_texture_t));

    md_mutex_lock(&dev->device_mutex);
    md_gpu_pool_t* slot = (md_gpu_pool_t*)md_vk_vec_push(&dev->pools, dev->alloc);
    if (slot) *slot = p;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { md_free(dev->alloc, p, sizeof(*p)); md_vk_fail("out of memory"); return NULL; }
    return p;
}

md_gpu_mem_kind_t md_gpu_pool_kind(md_gpu_pool_t p) { return p ? p->kind : MD_GPU_MEM_DEVICE; }

void md_gpu_pool_stats(md_gpu_pool_t p, md_gpu_pool_stats_t* out) {
    if (!p || !out) return;
    memset(out, 0, sizeof(*out));
    md_mutex_lock(&p->device->device_mutex);
    out->bytes_in_use      = p->in_use_bytes;
    out->bytes_reserved    = p->reserved_bytes;
    out->bytes_cached      = p->reserved_bytes - p->in_use_bytes;
    out->bytes_peak_in_use = p->peak_in_use_bytes;
    out->alloc_count       = p->alloc_count;
    out->reuse_count       = p->reuse_count;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (MD_VK_VEC_AT(p->blocks, md_vk_block_t*, i)->in_use) out->blocks_in_use++;
        else                                                   out->blocks_cached++;
    }
    md_mutex_unlock(&p->device->device_mutex);
}

/* Mark a block free at the current point in `stream`. Caller holds device_mutex. */
static void md_vk_block_release(md_vk_block_t* b, md_gpu_stream_t stream) {
    b->in_use      = false;
    b->free_stream = stream;
    b->free_value  = md_vk_stream_position(stream);
    if (b->free_value == 0 || md_vk_stream_completed(stream) >= b->free_value) {
        b->free_stream = NULL;
        b->free_value  = 0;
    }
    b->pool->in_use_bytes -= b->capacity;
}

static bool md_vk_block_idle(md_vk_block_t* b) {
    return !b->in_use && (b->free_stream == NULL || md_vk_stream_completed(b->free_stream) >= b->free_value);
}

/* Destroy the block's device objects and the block itself. It must already be
   out of the registry and out of its pool. */
static void md_vk_block_free(md_gpu_device_t dev, md_vk_block_t* b) {
    md_vk_destroy_raw_buffer(dev, b->buffer, b->memory);
    md_free(dev->alloc, b, sizeof(*b));
}

/* Release idle cached blocks while more than `keep_bytes` are cached.
   Caller holds device_mutex. */
static void md_vk_pool_trim_locked(md_gpu_pool_t p, uint64_t keep_bytes) {
    md_gpu_device_t dev = p->device;
    for (size_t i = 0; i < p->blocks.count && p->reserved_bytes - p->in_use_bytes > keep_bytes;) {
        md_vk_block_t* b = MD_VK_VEC_AT(p->blocks, md_vk_block_t*, i);
        if (md_vk_block_idle(b)) {
            md_vk_registry_remove_locked(dev, b);
            p->reserved_bytes -= b->capacity;
            md_vk_vec_remove(&p->blocks, i);
            md_vk_block_free(dev, b);
            continue;
        }
        ++i;
    }
}

void md_gpu_pool_trim(md_gpu_pool_t p, uint64_t keep_bytes) {
    if (!p) return;
    md_mutex_lock(&p->device->device_mutex);
    md_vk_pool_trim_locked(p, keep_bytes);
    md_mutex_unlock(&p->device->device_mutex);
}

/* Hand a live texture to the retire list. Caller holds device_mutex and has
   already removed it from its pool's list. */
static void md_vk_texture_retire_locked(md_gpu_device_t dev, md_gpu_texture_t t) {
    if (t->pool) {
        t->pool->in_use_bytes   -= t->bytes;
        t->pool->reserved_bytes -= t->bytes;
        t->pool = NULL;
    }
    md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
}

void md_gpu_pool_reset(md_gpu_stream_t stream, md_gpu_pool_t p) {
    if (!p || !stream) { md_vk_fail("md_gpu_pool_reset: null argument"); return; }
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_vk_block_t* b = MD_VK_VEC_AT(p->blocks, md_vk_block_t*, i);
        if (b->in_use) md_vk_block_release(b, stream);
    }
    while (p->textures.count > 0) {
        md_gpu_texture_t t = MD_VK_VEC_AT(p->textures, md_gpu_texture_t, p->textures.count - 1);
        p->textures.count--;
        md_vk_texture_retire_locked(dev, t);
    }
    md_mutex_unlock(&dev->device_mutex);
}

void md_gpu_pool_destroy(md_gpu_pool_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    /* Nothing is freed here: every block and texture goes onto the retire list
       stamped with every stream's current position, and md_gpu_device_poll
       releases it once all of those have passed. A block that is idle already
       is still routed through the list -- it costs nothing and keeps one path. */
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_vk_block_t* b = MD_VK_VEC_AT(p->blocks, md_vk_block_t*, i);
        md_vk_registry_remove_locked(dev, b);
        b->pool = NULL;
        md_vk_retire_locked(dev, MD_VK_RETIRE_BLOCK, b);
    }
    md_vk_vec_free(&p->blocks, dev->alloc);
    while (p->textures.count > 0) {
        md_gpu_texture_t t = MD_VK_VEC_AT(p->textures, md_gpu_texture_t, p->textures.count - 1);
        p->textures.count--;
        md_vk_texture_retire_locked(dev, t);
    }
    md_vk_vec_free(&p->textures, dev->alloc);
    md_vk_vec_remove_ptr(&dev->pools, p);
    /* Idle objects can go straight away. */
    md_vk_process_retires_locked(dev, false);
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, p, sizeof(*p));
}

md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t stream, md_gpu_pool_t p, size_t size) {
    md_gpu_mem_t out = {0, NULL};
    if (!stream || !p) { md_vk_fail("md_gpu_malloc: null stream or pool"); return out; }
    if (size == 0) return out;
    md_gpu_device_t dev = p->device;
    if (stream->device != dev) { md_vk_fail("md_gpu_malloc: stream and pool belong to different devices"); return out; }
    uint64_t capacity = md_vk_next_pow2(size);

    md_mutex_lock(&dev->device_mutex);

    /* Reuse: big enough, and either freed on this stream (program order makes
       it safe immediately) or its free point has completed. Never waits. */
    md_vk_block_t* best = NULL;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_vk_block_t* b = MD_VK_VEC_AT(p->blocks, md_vk_block_t*, i);
        if (b->in_use || b->capacity < size) continue;
        bool safe = (b->free_stream == NULL)
                 || (b->free_stream == stream)
                 || (md_vk_stream_completed(b->free_stream) >= b->free_value);
        if (!safe) continue;
        if (!best || b->capacity < best->capacity) best = b;
    }

    if (best) {
        best->in_use      = true;
        best->size        = size;
        best->free_stream = NULL;
        best->free_value  = 0;
        p->in_use_bytes  += best->capacity;
        if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
        p->alloc_count++;
        p->reuse_count++;
        md_mutex_unlock(&dev->device_mutex);
        out.gpu = best->address;
        out.cpu = best->host;
        return out;
    }

    md_vk_block_t* b = (md_vk_block_t*)md_alloc(dev->alloc, sizeof(md_vk_block_t));
    if (!b) { md_mutex_unlock(&dev->device_mutex); md_vk_fail("out of memory"); return out; }
    memset(b, 0, sizeof(*b));
    if (!md_vk_create_raw_buffer(dev, capacity, p->kind, &b->buffer, &b->memory, &b->address, &b->host)) {
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        return out;
    }
    b->capacity = capacity;
    b->size     = size;
    b->kind     = p->kind;
    b->pool     = p;
    b->in_use   = true;

    md_vk_block_t** slot = (md_vk_block_t**)md_vk_vec_push(&p->blocks, dev->alloc);
    if (!slot || !md_vk_registry_insert_locked(dev, b)) {
        if (slot) p->blocks.count--;
        md_vk_block_free(dev, b);
        md_mutex_unlock(&dev->device_mutex);
        md_vk_fail("out of memory");
        return out;
    }
    *slot = b;
    p->reserved_bytes += capacity;
    p->in_use_bytes   += capacity;
    if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
    p->alloc_count++;

    md_mutex_unlock(&dev->device_mutex);
    out.gpu = b->address;
    out.cpu = b->host;
    return out;
}

void md_gpu_free(md_gpu_stream_t stream, md_gpu_addr_t addr) {
    if (!addr) return;
    if (!stream) { md_vk_fail("md_gpu_free: a stream is required"); return; }
    md_gpu_device_t dev = stream->device;

    md_mutex_lock(&dev->device_mutex);
    md_vk_block_t* b = md_vk_registry_find_locked(dev, addr);
    if (!b || !b->in_use || b->address != addr) {
        md_mutex_unlock(&dev->device_mutex);
        md_vk_fail("md_gpu_free: 0x%llx is not the start of a live allocation", (unsigned long long)addr);
        return;
    }
    md_vk_block_release(b, stream);
    md_gpu_pool_t p = b->pool;
    if (p->cache_limit != 0) md_vk_pool_trim_locked(p, p->cache_limit);
    md_mutex_unlock(&dev->device_mutex);
}

/* ---- Copies ------------------------------------------------------------------ */

static bool md_vk_record_buffer_copy(md_gpu_stream_t s, VkBuffer src, uint64_t src_off,
                                     VkBuffer dst, uint64_t dst_off, uint64_t size) {
    VkCommandBuffer cmd = md_vk_begin_op(s);
    if (!cmd) return false;
    VkBufferCopy region = {src_off, dst_off, size};
    vkCmdCopyBuffer(cmd, src, dst, 1, &region);
    md_vk_end_op(s);
    return true;
}

bool md_gpu_copy(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_addr_t src, size_t size) {
    if (!s) return md_vk_fail("md_gpu_copy: null stream");
    if (size == 0) return true;
    uint64_t doff, soff;
    md_vk_block_t* d  = md_vk_resolve(s->device, dst, size, &doff, "md_gpu_copy (dst)");
    if (!d) return false;
    md_vk_block_t* sb = md_vk_resolve(s->device, src, size, &soff, "md_gpu_copy (src)");
    if (!sb) return false;
    return md_vk_record_buffer_copy(s, sb->buffer, soff, d->buffer, doff, size);
}

/* True when a host write into `b` cannot race this stream's GPU work: the
   stream has nothing recorded and everything it submitted has completed. */
/* Idle means a host write now cannot race anything this stream is ordered
   after: nothing recorded, nothing in flight, and every cross-stream wait
   queued by md_gpu_stream_wait already satisfied. */
static bool md_vk_stream_idle(md_gpu_stream_t s) {
    if (s->has_work || md_vk_stream_completed(s) < s->submitted_value) return false;
    for (size_t i = 0; i < s->waits.count; ++i) {
        md_vk_wait_t w = MD_VK_VEC_AT(s->waits, md_vk_wait_t, i);
        if (md_vk_stream_completed(w.stream) < w.value) return false;
    }
    return true;
}

bool md_gpu_upload(md_gpu_stream_t s, md_gpu_addr_t dst, const void* src, size_t size) {
    if (!s) return md_vk_fail("md_gpu_upload: null stream");
    if (size == 0) return true;
    if (!src) return md_vk_fail("md_gpu_upload: null source");
    if (s->upload_open) return md_vk_fail("md_gpu_upload: stream '%s' has an open upload", s->label);
    uint64_t doff;
    md_vk_block_t* d = md_vk_resolve(s->device, dst, size, &doff, "md_gpu_upload");
    if (!d) return false;

    /* Fast path: host-writable destination and nothing in flight here. */
    if (d->host && md_vk_stream_idle(s)) {
        memcpy((uint8_t*)d->host + doff, src, size);
        return true;
    }
    uint64_t addr; void* host; VkBuffer buf; uint64_t off;
    if (!md_vk_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    return md_vk_record_buffer_copy(s, buf, off, d->buffer, doff, size);
}

bool md_gpu_memset(md_gpu_stream_t s, md_gpu_addr_t dst, uint8_t value, size_t size) {
    if (!s) return md_vk_fail("md_gpu_memset: null stream");
    if (size == 0) return true;
    uint64_t off;
    md_vk_block_t* b = md_vk_resolve(s->device, dst, size, &off, "md_gpu_memset");
    if (!b) return false;

    const uint64_t begin = off;
    const uint64_t end   = off + size;
    uint64_t abeg = md_vk_align_up(begin, 4);
    uint64_t aend = end & ~3ull;

    if (aend > abeg) {
        VkCommandBuffer cmd = md_vk_begin_op(s);
        if (!cmd) return false;
        vkCmdFillBuffer(cmd, b->buffer, abeg, aend - abeg, ((uint32_t)value) * 0x01010101u);
        md_vk_end_op(s);
    }

    /* Unaligned head and tail go through a small staged copy. */
    uint64_t head_off = begin, head = 0, tail = 0;
    if (aend > abeg) {
        head = abeg - begin;
        tail = end - aend;
    } else {
        head = size;   /* the whole range is smaller than one aligned word */
    }
    uint64_t pieces[2][2] = {{head_off, head}, {aend, tail}};
    for (int i = 0; i < 2; ++i) {
        if (pieces[i][1] == 0) continue;
        uint64_t addr; void* host; VkBuffer buf; uint64_t soff;
        if (!md_vk_arena_alloc(s, (size_t)pieces[i][1], &addr, &host, &buf, &soff)) return false;
        memset(host, value, (size_t)pieces[i][1]);
        if (!md_vk_record_buffer_copy(s, buf, soff, b->buffer, pieces[i][0], pieces[i][1])) return false;
    }
    return true;
}

void* md_gpu_upload_begin(md_gpu_stream_t s, md_gpu_addr_t dst, size_t size) {
    if (!s || !dst || size == 0) { md_vk_fail("md_gpu_upload_begin: null argument"); return NULL; }
    if (s->upload_open) { md_vk_fail("an upload is already open on stream '%s'", s->label); return NULL; }
    uint64_t doff;
    md_vk_block_t* b = md_vk_resolve(s->device, dst, size, &doff, "md_gpu_upload_begin");
    if (!b) return NULL;

    /* Write straight into the destination when that cannot race the GPU. */
    if (b->host && md_vk_stream_idle(s)) {
        s->upload_open   = true;
        s->upload_direct = true;
        s->upload_dst    = dst;
        s->upload_size   = size;
        return (uint8_t*)b->host + doff;
    }

    uint64_t addr; void* host;
    if (!md_vk_arena_alloc(s, size, &addr, &host, NULL, NULL)) return NULL;
    s->upload_open     = true;
    s->upload_direct   = false;
    s->upload_dst      = dst;
    s->upload_src_addr = addr;
    s->upload_size     = size;
    return host;
}

bool md_gpu_upload_end(md_gpu_stream_t s) {
    if (!s || !s->upload_open) return md_vk_fail("md_gpu_upload_end: no upload is open");
    s->upload_open = false;
    if (s->upload_direct) return true;

    uint64_t doff;
    md_vk_block_t* d = md_vk_resolve(s->device, s->upload_dst, s->upload_size, &doff, "md_gpu_upload_end");
    if (!d) return false;
    /* Find the staging page again; it is one of this stream's arena pages. */
    for (size_t i = 0; i < s->arena.pages.count; ++i) {
        md_vk_page_t* p = MD_VK_VEC_AT(s->arena.pages, md_vk_page_t*, i);
        if (s->upload_src_addr >= p->address && s->upload_src_addr < p->address + p->capacity) {
            return md_vk_record_buffer_copy(s, p->buffer, s->upload_src_addr - p->address,
                                            d->buffer, doff, s->upload_size);
        }
    }
    return md_vk_fail("md_gpu_upload_end: staging page not found");
}

/* =========================================================================
   10. Textures and samplers
   ========================================================================= */

static md_vk_format_info_t md_vk_format_info(md_gpu_format_t f) {
    const VkImageAspectFlags C = VK_IMAGE_ASPECT_COLOR_BIT;
    const VkImageAspectFlags D = VK_IMAGE_ASPECT_DEPTH_BIT;
    md_vk_format_info_t i = {VK_FORMAT_UNDEFINED, 0, C, C, false, "invalid"};
#define MD_VK_FMT(e, vk, b, va, ca, dep) case e: i = (md_vk_format_info_t){vk, b, va, ca, dep, #e}; break
    switch (f) {
    MD_VK_FMT(MD_GPU_FORMAT_R8_UNORM,          VK_FORMAT_R8_UNORM,                 1, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RG8_UNORM,         VK_FORMAT_R8G8_UNORM,               2, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGBA8_UNORM,       VK_FORMAT_R8G8B8A8_UNORM,           4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGBA8_SRGB,        VK_FORMAT_R8G8B8A8_SRGB,            4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_BGRA8_UNORM,       VK_FORMAT_B8G8R8A8_UNORM,           4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_BGRA8_SRGB,        VK_FORMAT_B8G8R8A8_SRGB,            4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_R16_FLOAT,         VK_FORMAT_R16_SFLOAT,               2, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RG16_FLOAT,        VK_FORMAT_R16G16_SFLOAT,            4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGBA16_FLOAT,      VK_FORMAT_R16G16B16A16_SFLOAT,      8, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_R32_FLOAT,         VK_FORMAT_R32_SFLOAT,               4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RG32_FLOAT,        VK_FORMAT_R32G32_SFLOAT,            8, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGBA32_FLOAT,      VK_FORMAT_R32G32B32A32_SFLOAT,     16, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_R32_UINT,          VK_FORMAT_R32_UINT,                 4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RG32_UINT,         VK_FORMAT_R32G32_UINT,              8, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGBA32_UINT,       VK_FORMAT_R32G32B32A32_UINT,       16, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RG11B10_FLOAT,     VK_FORMAT_B10G11R11_UFLOAT_PACK32,  4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_RGB10A2_UNORM,     VK_FORMAT_A2B10G10R10_UNORM_PACK32, 4, C, C, false);
    MD_VK_FMT(MD_GPU_FORMAT_D32_FLOAT,         VK_FORMAT_D32_SFLOAT,               4, D, D, true);
    MD_VK_FMT(MD_GPU_FORMAT_D32_FLOAT_S8_UINT, VK_FORMAT_D32_SFLOAT_S8_UINT,       4, D, D, true);
    default: break;
    }
#undef MD_VK_FMT
    return i;
}

uint32_t md_gpu_format_texel_size(md_gpu_format_t format) {
    return md_vk_format_info(format).bytes;
}

static const char* md_vk_tex_type_name(md_gpu_tex_type_t t) {
    switch (t) {
    case MD_GPU_TEX_2D:       return "2D";
    case MD_GPU_TEX_2D_ARRAY: return "2D_ARRAY";
    case MD_GPU_TEX_3D:       return "3D";
    default:                  return "INVALID";
    }
}

static VkImageViewType md_vk_view_type(md_gpu_tex_type_t t) {
    switch (t) {
    case MD_GPU_TEX_2D_ARRAY: return VK_IMAGE_VIEW_TYPE_2D_ARRAY;
    case MD_GPU_TEX_3D:       return VK_IMAGE_VIEW_TYPE_3D;
    default:                  return VK_IMAGE_VIEW_TYPE_2D;
    }
}

/* Take a slot from the heap free list. Caller holds device_mutex. */
static uint32_t md_vk_alloc_slot_locked(md_gpu_device_t dev) {
    if (dev->tex_free_count == 0) return 0;
    return dev->tex_free[--dev->tex_free_count];
}

static void md_vk_free_slot_locked(md_gpu_device_t dev, uint32_t slot, VkDescriptorType type) {
    if (slot == 0) return;
    md_vk_clear_slot(dev, slot, type);
    if (dev->tex_free_count < MD_VK_MAX_TEXTURE_SLOTS) dev->tex_free[dev->tex_free_count++] = slot;
}

/* Destroy a texture's device objects, return its heap slots and free it. It
   must no longer be referenced by any work in flight. */
static void md_vk_texture_free(md_gpu_device_t dev, md_gpu_texture_t t) {
    const uint32_t mips = t->desc.mip_levels;
    if (t->storage_slots) {
        for (uint32_t m = 0; m < mips; ++m) md_vk_free_slot_locked(dev, t->storage_slots[m], VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
        md_free(dev->alloc, t->storage_slots, mips * sizeof(uint32_t));
    }
    md_vk_free_slot_locked(dev, t->sampled_slot, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE);
    if (t->storage_views) {
        for (uint32_t m = 0; m < mips; ++m) if (t->storage_views[m]) vkDestroyImageView(dev->device, t->storage_views[m], NULL);
        md_free(dev->alloc, t->storage_views, mips * sizeof(VkImageView));
    }
    if (t->sampled_view) vkDestroyImageView(dev->device, t->sampled_view, NULL);
    if (t->image)  vkDestroyImage(dev->device, t->image, NULL);
    if (t->memory) vkFreeMemory(dev->device, t->memory, NULL);
    md_free(dev->alloc, t, sizeof(*t));
}

static VkImageView md_vk_create_view(md_gpu_device_t dev, md_gpu_texture_t t, uint32_t base_mip, uint32_t mip_count) {
    VkImageViewCreateInfo vci = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    vci.image    = t->image;
    vci.viewType = md_vk_view_type(t->desc.type);
    vci.format   = t->fi.format;
    vci.subresourceRange.aspectMask     = t->fi.view_aspect;
    vci.subresourceRange.baseMipLevel   = base_mip;
    vci.subresourceRange.levelCount     = mip_count;
    vci.subresourceRange.baseArrayLayer = 0;
    vci.subresourceRange.layerCount     = t->desc.type == MD_GPU_TEX_2D_ARRAY ? t->desc.depth_or_layers : 1;
    VkImageView view = VK_NULL_HANDLE;
    if (!md_vk_check(vkCreateImageView(dev->device, &vci, NULL, &view), "vkCreateImageView")) return VK_NULL_HANDLE;
    return view;
}

md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t s, md_gpu_pool_t pool, const md_gpu_texture_desc_t* desc) {
    if (!s || !pool || !desc) { md_vk_fail("md_gpu_texture_create: null argument"); return NULL; }
    md_gpu_device_t dev = s->device;
    if (pool->device != dev)             { md_vk_fail("md_gpu_texture_create: stream and pool belong to different devices"); return NULL; }
    if (pool->kind != MD_GPU_MEM_DEVICE) { md_vk_fail("md_gpu_texture_create: pool '%s' is not an MD_GPU_MEM_DEVICE pool", pool->label); return NULL; }
    if (s->upload_open)                  { md_vk_fail("md_gpu_texture_create: stream '%s' has an open upload", s->label); return NULL; }

    md_gpu_texture_desc_t d = *desc;
    const char* label = d.label ? d.label : "texture";
    md_vk_format_info_t fi = md_vk_format_info(d.format);
    if (fi.format == VK_FORMAT_UNDEFINED) { md_vk_fail("texture '%s': invalid format %d", label, (int)d.format); return NULL; }
    if (d.type != MD_GPU_TEX_2D && d.type != MD_GPU_TEX_2D_ARRAY && d.type != MD_GPU_TEX_3D) {
        md_vk_fail("texture '%s': invalid type %d (a zero-initialised desc has no type)", label, (int)d.type);
        return NULL;
    }
    if (!(d.usage & (MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED | MD_GPU_TEX_RENDER_TARGET))) {
        md_vk_fail("texture '%s': usage is empty", label);
        return NULL;
    }
    if (d.width == 0 || d.height == 0) { md_vk_fail("texture '%s': zero width or height", label); return NULL; }
    if (d.type == MD_GPU_TEX_2D && d.depth_or_layers > 1) {
        md_vk_fail("texture '%s': a 2D texture has depth_or_layers %u; use MD_GPU_TEX_3D or MD_GPU_TEX_2D_ARRAY",
                   label, d.depth_or_layers);
        return NULL;
    }
    if (d.depth_or_layers == 0) d.depth_or_layers = 1;
    if (d.mip_levels == 0) d.mip_levels = 1;

    VkImageCreateInfo ici = {VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    ici.imageType     = d.type == MD_GPU_TEX_3D ? VK_IMAGE_TYPE_3D : VK_IMAGE_TYPE_2D;
    ici.format        = fi.format;
    ici.extent.width  = d.width;
    ici.extent.height = d.height;
    ici.extent.depth  = d.type == MD_GPU_TEX_3D ? d.depth_or_layers : 1;
    ici.mipLevels     = d.mip_levels;
    ici.arrayLayers   = d.type == MD_GPU_TEX_2D_ARRAY ? d.depth_or_layers : 1;
    ici.samples       = VK_SAMPLE_COUNT_1_BIT;
    ici.tiling        = VK_IMAGE_TILING_OPTIMAL;
    ici.usage         = VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    if (d.usage & MD_GPU_TEX_STORAGE) ici.usage |= VK_IMAGE_USAGE_STORAGE_BIT;
    if (d.usage & MD_GPU_TEX_SAMPLED) ici.usage |= VK_IMAGE_USAGE_SAMPLED_BIT;
    if (d.usage & MD_GPU_TEX_RENDER_TARGET) {
        ici.usage |= fi.depth ? VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT : VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    }
    md_vk_set_sharing(dev, &ici.sharingMode, &ici.queueFamilyIndexCount, &ici.pQueueFamilyIndices);
    ici.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    /* Reject what the device cannot do, by name, before anything is created. */
    {
        VkImageFormatProperties ifp;
        VkResult r = vkGetPhysicalDeviceImageFormatProperties(dev->phys, ici.format, ici.imageType, ici.tiling, ici.usage, 0, &ifp);
        if (r != VK_SUCCESS) {
            md_vk_fail("texture '%s': the device does not support %s %s with usage%s%s%s",
                       label, md_vk_tex_type_name(d.type), fi.name,
                       (d.usage & MD_GPU_TEX_STORAGE) ? " STORAGE" : "",
                       (d.usage & MD_GPU_TEX_SAMPLED) ? " SAMPLED" : "",
                       (d.usage & MD_GPU_TEX_RENDER_TARGET) ? " RENDER_TARGET" : "");
            return NULL;
        }
        if (d.width > ifp.maxExtent.width || d.height > ifp.maxExtent.height || ici.extent.depth > ifp.maxExtent.depth ||
            d.mip_levels > ifp.maxMipLevels || ici.arrayLayers > ifp.maxArrayLayers) {
            md_vk_fail("texture '%s': %ux%ux%u with %u mips and %u layers exceeds the device limits for %s",
                       label, d.width, d.height, ici.extent.depth, d.mip_levels, ici.arrayLayers, fi.name);
            return NULL;
        }
    }

    md_gpu_texture_t t = (md_gpu_texture_t)md_alloc(dev->alloc, sizeof(md_gpu_texture));
    if (!t) { md_vk_fail("out of memory"); return NULL; }
    memset(t, 0, sizeof(*t));
    t->device = dev;
    t->fi     = fi;
    t->desc   = d;
    snprintf(t->label, sizeof(t->label), "%s", label);
    t->desc.label = t->label;

    if (!md_vk_check(vkCreateImage(dev->device, &ici, NULL, &t->image), "vkCreateImage")) goto fail;

    {
        VkMemoryRequirements req;
        vkGetImageMemoryRequirements(dev->device, t->image, &req);
        uint32_t type = md_vk_find_memory_type(dev, req.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0);
        if (type == UINT32_MAX) type = md_vk_find_memory_type(dev, req.memoryTypeBits, 0, 0);
        VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        mai.allocationSize  = req.size;
        mai.memoryTypeIndex = type;
        if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &t->memory), "vkAllocateMemory (texture)")) goto fail;
        if (!md_vk_check(vkBindImageMemory(dev->device, t->image, t->memory, 0), "vkBindImageMemory")) goto fail;
        t->bytes = req.size;
    }

    if (d.usage & MD_GPU_TEX_SAMPLED) {
        t->sampled_view = md_vk_create_view(dev, t, 0, d.mip_levels);
        if (!t->sampled_view) goto fail;
    }
    if (d.usage & MD_GPU_TEX_STORAGE) {
        t->storage_views = (VkImageView*)md_alloc(dev->alloc, d.mip_levels * sizeof(VkImageView));
        t->storage_slots = (uint32_t*)md_alloc(dev->alloc, d.mip_levels * sizeof(uint32_t));
        if (!t->storage_views || !t->storage_slots) { md_vk_fail("out of memory"); goto fail; }
        memset(t->storage_views, 0, d.mip_levels * sizeof(VkImageView));
        memset(t->storage_slots, 0, d.mip_levels * sizeof(uint32_t));
        for (uint32_t m = 0; m < d.mip_levels; ++m) {
            t->storage_views[m] = md_vk_create_view(dev, t, m, 1);
            if (!t->storage_views[m]) goto fail;
        }
    }

    /* Record the one layout transition this image ever gets into the stream,
       so creation is stream-ordered like md_gpu_malloc and never waits for
       anything. The barrier's second scope is ALL_COMMANDS, so it orders
       itself before every later operation in the stream in either ordering
       mode. */
    {
        if (!md_vk_stream_ensure_cmd(s)) goto fail;
        VkImageMemoryBarrier2 b = {VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2};
        b.srcStageMask  = VK_PIPELINE_STAGE_2_NONE;
        b.srcAccessMask = 0;
        b.dstStageMask  = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
        b.dstAccessMask = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
        b.oldLayout     = VK_IMAGE_LAYOUT_UNDEFINED;
        b.newLayout     = VK_IMAGE_LAYOUT_GENERAL;
        b.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        b.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        b.image = t->image;
        b.subresourceRange.aspectMask = fi.depth
            ? (VK_IMAGE_ASPECT_DEPTH_BIT | (d.format == MD_GPU_FORMAT_D32_FLOAT_S8_UINT ? VK_IMAGE_ASPECT_STENCIL_BIT : 0))
            : VK_IMAGE_ASPECT_COLOR_BIT;
        b.subresourceRange.levelCount = VK_REMAINING_MIP_LEVELS;
        b.subresourceRange.layerCount = VK_REMAINING_ARRAY_LAYERS;
        VkDependencyInfo di = {VK_STRUCTURE_TYPE_DEPENDENCY_INFO};
        di.imageMemoryBarrierCount = 1;
        di.pImageMemoryBarriers    = &b;
        vkCmdPipelineBarrier2(s->open, &di);
        s->has_work = true;
    }

    md_mutex_lock(&dev->device_mutex);
    {
        uint32_t needed = ((d.usage & MD_GPU_TEX_SAMPLED) ? 1u : 0u) + ((d.usage & MD_GPU_TEX_STORAGE) ? d.mip_levels : 0u);
        if (dev->tex_free_count < needed) {
            md_mutex_unlock(&dev->device_mutex);
            md_vk_fail("texture '%s': out of bindless heap slots (max %u)", label, MD_VK_MAX_TEXTURE_SLOTS - 1);
            /* The transition is already recorded; destroy through the retire
               path so the image outlives that command buffer. */
            md_mutex_lock(&dev->device_mutex);
            md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
            md_mutex_unlock(&dev->device_mutex);
            return NULL;
        }
        if (d.usage & MD_GPU_TEX_SAMPLED) {
            t->sampled_slot = md_vk_alloc_slot_locked(dev);
            md_vk_write_image_slot(dev, t->sampled_slot, t->sampled_view, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE);
        }
        if (d.usage & MD_GPU_TEX_STORAGE) {
            for (uint32_t m = 0; m < d.mip_levels; ++m) {
                t->storage_slots[m] = md_vk_alloc_slot_locked(dev);
                md_vk_write_image_slot(dev, t->storage_slots[m], t->storage_views[m], VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
            }
        }
        md_gpu_texture_t* slot = (md_gpu_texture_t*)md_vk_vec_push(&pool->textures, dev->alloc);
        if (!slot) {
            md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
            md_mutex_unlock(&dev->device_mutex);
            md_vk_fail("out of memory");
            return NULL;
        }
        *slot = t;
        t->pool = pool;
        pool->in_use_bytes   += t->bytes;
        pool->reserved_bytes += t->bytes;
        if (pool->in_use_bytes > pool->peak_in_use_bytes) pool->peak_in_use_bytes = pool->in_use_bytes;
    }
    md_mutex_unlock(&dev->device_mutex);
    return t;

fail:
    /* Nothing referencing the image has been recorded yet. */
    md_mutex_lock(&dev->device_mutex);
    md_vk_texture_free(dev, t);
    md_mutex_unlock(&dev->device_mutex);
    return NULL;
}

void md_gpu_texture_destroy(md_gpu_texture_t t) {
    if (!t) return;
    md_gpu_device_t dev = t->device;
    md_mutex_lock(&dev->device_mutex);
    if (t->pool) md_vk_vec_remove_ptr(&t->pool->textures, t);
    md_vk_texture_retire_locked(dev, t);
    md_mutex_unlock(&dev->device_mutex);
}

const md_gpu_texture_desc_t* md_gpu_texture_desc(md_gpu_texture_t t) {
    return t ? &t->desc : NULL;
}

md_gpu_storage_tex_t md_gpu_texture_storage(md_gpu_texture_t t, uint32_t mip) {
    md_gpu_storage_tex_t h = {0};
    if (t && t->storage_slots && mip < t->desc.mip_levels) h.handle = t->storage_slots[mip];
    return h;
}

md_gpu_sampled_tex_t md_gpu_texture_sampled(md_gpu_texture_t t) {
    md_gpu_sampled_tex_t h = {0};
    if (t) h.handle = t->sampled_slot;
    return h;
}

/* ---- Samplers ------------------------------------------------------------------ */

static bool md_vk_sampler_desc_eq(const md_gpu_sampler_desc_t* a, const md_gpu_sampler_desc_t* b) {
    return a->min_filter == b->min_filter && a->mag_filter == b->mag_filter && a->mip_filter == b->mip_filter &&
           a->address_u == b->address_u && a->address_v == b->address_v && a->address_w == b->address_w;
}

md_gpu_sampler_t md_gpu_sampler(md_gpu_device_t dev, const md_gpu_sampler_desc_t* desc) {
    md_gpu_sampler_t h = {0};
    if (!dev) { md_vk_fail("md_gpu_sampler: null device"); return h; }
    md_gpu_sampler_desc_t d;
    memset(&d, 0, sizeof(d));
    if (desc) d = *desc;

    md_mutex_lock(&dev->device_mutex);
    for (uint32_t i = 0; i < dev->sampler_count; ++i) {
        if (md_vk_sampler_desc_eq(&dev->samplers[i].desc, &d)) {
            h.handle = dev->samplers[i].slot;
            md_mutex_unlock(&dev->device_mutex);
            return h;
        }
    }
    if (dev->sampler_count + 1 >= MD_VK_MAX_SAMPLERS) {
        md_mutex_unlock(&dev->device_mutex);
        md_vk_fail("md_gpu_sampler: more than %u distinct samplers", MD_VK_MAX_SAMPLERS - 1);
        return h;
    }

    static const VkSamplerAddressMode modes[] = {
        VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE,
        VK_SAMPLER_ADDRESS_MODE_REPEAT,
        VK_SAMPLER_ADDRESS_MODE_MIRRORED_REPEAT,
    };
    VkSamplerCreateInfo sci = {VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    sci.minFilter    = d.min_filter == MD_GPU_FILTER_LINEAR ? VK_FILTER_LINEAR : VK_FILTER_NEAREST;
    sci.magFilter    = d.mag_filter == MD_GPU_FILTER_LINEAR ? VK_FILTER_LINEAR : VK_FILTER_NEAREST;
    sci.mipmapMode   = d.mip_filter == MD_GPU_FILTER_LINEAR ? VK_SAMPLER_MIPMAP_MODE_LINEAR : VK_SAMPLER_MIPMAP_MODE_NEAREST;
    sci.addressModeU = modes[(unsigned)d.address_u % 3u];
    sci.addressModeV = modes[(unsigned)d.address_v % 3u];
    sci.addressModeW = modes[(unsigned)d.address_w % 3u];
    sci.maxLod       = VK_LOD_CLAMP_NONE;
    sci.borderColor  = VK_BORDER_COLOR_FLOAT_TRANSPARENT_BLACK;

    VkSampler sampler;
    if (!md_vk_check(vkCreateSampler(dev->device, &sci, NULL, &sampler), "vkCreateSampler")) {
        md_mutex_unlock(&dev->device_mutex);
        return h;
    }
    /* Slot 0 stays the null handle; entry i lives in slot i + 1. */
    const uint32_t slot = dev->sampler_count + 1;
    md_vk_sampler_entry_t* e = &dev->samplers[dev->sampler_count++];
    e->desc    = d;
    e->sampler = sampler;
    e->slot    = slot;

    VkDescriptorImageInfo dii = {sampler, VK_NULL_HANDLE, VK_IMAGE_LAYOUT_UNDEFINED};
    VkWriteDescriptorSet w = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    w.dstSet          = dev->desc_set;
    w.dstBinding      = MD_VK_BINDING_SAMPLER;
    w.dstArrayElement = slot;
    w.descriptorCount = 1;
    w.descriptorType  = VK_DESCRIPTOR_TYPE_SAMPLER;
    w.pImageInfo      = &dii;
    vkUpdateDescriptorSets(dev->device, 1, &w, 0, NULL);
    md_mutex_unlock(&dev->device_mutex);

    h.handle = slot;
    return h;
}

/* ---- Texture copies ------------------------------------------------------------ */

typedef struct md_vk_copy_region_t {
    VkBufferImageCopy copy;
    uint64_t          bytes;
} md_vk_copy_region_t;

/* Resolve a region against a texture: defaults, bounds, and the equivalent
   VkBufferImageCopy (buffer offset left 0). */
static bool md_vk_resolve_region(md_gpu_texture_t t, const md_gpu_tex_region_t* r, md_vk_copy_region_t* out, const char* what) {
    md_gpu_tex_region_t z;
    memset(&z, 0, sizeof(z));
    if (!r) r = &z;
    const md_gpu_texture_desc_t* d = &t->desc;
    if (r->mip >= d->mip_levels) return md_vk_fail("%s: mip %u out of range (texture '%s' has %u)", what, r->mip, t->label, d->mip_levels);

    uint32_t dim[3];
    dim[0] = d->width  >> r->mip; if (!dim[0]) dim[0] = 1;
    dim[1] = d->height >> r->mip; if (!dim[1]) dim[1] = 1;
    if (d->type == MD_GPU_TEX_3D) { dim[2] = d->depth_or_layers >> r->mip; if (!dim[2]) dim[2] = 1; }
    else                          { dim[2] = d->depth_or_layers; }

    uint32_t ext[3];
    for (int i = 0; i < 3; ++i) {
        if (r->offset[i] >= dim[i]) {
            return md_vk_fail("%s: offset[%d] = %u is outside texture '%s' (extent %u at mip %u)",
                              what, i, r->offset[i], t->label, dim[i], r->mip);
        }
        ext[i] = r->extent[i] ? r->extent[i] : dim[i] - r->offset[i];
        if (ext[i] > dim[i] - r->offset[i]) {
            return md_vk_fail("%s: region [%u, +%u) on axis %d overruns texture '%s' (extent %u at mip %u)",
                              what, r->offset[i], ext[i], i, t->label, dim[i], r->mip);
        }
    }

    memset(out, 0, sizeof(*out));
    VkBufferImageCopy* c = &out->copy;
    c->imageSubresource.aspectMask = t->fi.copy_aspect;
    c->imageSubresource.mipLevel   = r->mip;
    c->imageOffset.x      = (int32_t)r->offset[0];
    c->imageOffset.y      = (int32_t)r->offset[1];
    c->imageExtent.width  = ext[0];
    c->imageExtent.height = ext[1];
    if (d->type == MD_GPU_TEX_2D_ARRAY) {
        c->imageSubresource.baseArrayLayer = r->offset[2];
        c->imageSubresource.layerCount     = ext[2];
        c->imageOffset.z                   = 0;
        c->imageExtent.depth               = 1;
    } else {
        c->imageSubresource.baseArrayLayer = 0;
        c->imageSubresource.layerCount     = 1;
        c->imageOffset.z                   = (int32_t)r->offset[2];
        c->imageExtent.depth               = ext[2];
    }
    out->bytes = (uint64_t)ext[0] * ext[1] * ext[2] * t->fi.bytes;
    return true;
}

size_t md_gpu_texture_region_size(md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!t) return 0;
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_texture_region_size")) return 0;
    return (size_t)cr.bytes;
}

static bool md_vk_check_buffer_offset(md_gpu_texture_t t, uint64_t off, const char* what) {
    /* Vulkan requires a multiple of the texel size, and of 4 for depth. */
    const uint64_t a = t->fi.depth ? 4 : t->fi.bytes;
    if (off % a != 0) {
        return md_vk_fail("%s: buffer address must be %llu-byte aligned for %s", what, (unsigned long long)a, t->fi.name);
    }
    return true;
}

static bool md_vk_record_texture_copy(md_gpu_stream_t s, md_gpu_texture_t t, VkBuffer buffer,
                                      const VkBufferImageCopy* copy, bool to_texture) {
    VkCommandBuffer cmd = md_vk_begin_op(s);
    if (!cmd) return false;
    if (to_texture) vkCmdCopyBufferToImage(cmd, buffer, t->image, VK_IMAGE_LAYOUT_GENERAL, 1, copy);
    else            vkCmdCopyImageToBuffer(cmd, t->image, VK_IMAGE_LAYOUT_GENERAL, buffer, 1, copy);
    md_vk_end_op(s);
    return true;
}

bool md_gpu_copy_to_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, md_gpu_addr_t src) {
    if (!s || !t) return md_vk_fail("md_gpu_copy_to_texture: null argument");
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_copy_to_texture")) return false;
    uint64_t off;
    md_vk_block_t* b = md_vk_resolve(s->device, src, cr.bytes, &off, "md_gpu_copy_to_texture");
    if (!b) return false;
    if (!md_vk_check_buffer_offset(t, off, "md_gpu_copy_to_texture")) return false;
    cr.copy.bufferOffset = off;
    return md_vk_record_texture_copy(s, t, b->buffer, &cr.copy, true);
}

bool md_gpu_copy_from_texture(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!s || !t) return md_vk_fail("md_gpu_copy_from_texture: null argument");
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_copy_from_texture")) return false;
    uint64_t off;
    md_vk_block_t* b = md_vk_resolve(s->device, dst, cr.bytes, &off, "md_gpu_copy_from_texture");
    if (!b) return false;
    if (!md_vk_check_buffer_offset(t, off, "md_gpu_copy_from_texture")) return false;
    cr.copy.bufferOffset = off;
    return md_vk_record_texture_copy(s, t, b->buffer, &cr.copy, false);
}

bool md_gpu_upload_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, const void* src, size_t size) {
    if (!s || !t || !src) return md_vk_fail("md_gpu_upload_texture: null argument");
    if (s->upload_open) return md_vk_fail("md_gpu_upload_texture: stream '%s' has an open upload", s->label);
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_upload_texture")) return false;
    if ((uint64_t)size != cr.bytes) {
        return md_vk_fail("md_gpu_upload_texture: region of texture '%s' is %llu bytes but %zu were given",
                          t->label, (unsigned long long)cr.bytes, size);
    }
    uint64_t addr; void* host; VkBuffer buf; uint64_t off;
    if (!md_vk_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    cr.copy.bufferOffset = off;
    return md_vk_record_texture_copy(s, t, buf, &cr.copy, true);
}

/* =========================================================================
   11. Kernels and launches
   ========================================================================= */

/* The local size recorded in SPIR-V: OpExecutionMode <entry> LocalSize x y z. */
static bool md_vk_spirv_local_size(const uint32_t* words, size_t word_count, uint32_t out[3]) {
    if (word_count < 5 || words[0] != 0x07230203u) return false;
    size_t i = 5;
    while (i < word_count) {
        uint32_t op    = words[i] & 0xFFFFu;
        uint32_t count = words[i] >> 16;
        if (count == 0 || i + count > word_count) break;
        if (op == 16 /* OpExecutionMode */ && count >= 6 && words[i + 2] == 17 /* LocalSize */) {
            out[0] = words[i + 3];
            out[1] = words[i + 4];
            out[2] = words[i + 5];
            return true;
        }
        i += count;
    }
    return false;
}

static void md_vk_kernel_free(md_gpu_device_t dev, md_gpu_kernel_t k) {
    if (k->pipeline) vkDestroyPipeline(dev->device, k->pipeline, NULL);
    if (k->module)   vkDestroyShaderModule(dev->device, k->module, NULL);
    md_free(dev->alloc, k, sizeof(*k));
}

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t dev, const md_gpu_kernel_desc_t* desc) {
    if (!dev || !desc || !desc->code || desc->code_size == 0) {
        md_vk_fail("md_gpu_kernel_create: missing code");
        return NULL;
    }
    const char* label = desc->label ? desc->label : "kernel";
    if (desc->code_size % 4 != 0) { md_vk_fail("kernel '%s': SPIR-V size must be a multiple of 4", label); return NULL; }
    if (desc->group_size[0] == 0 || desc->group_size[1] == 0 || desc->group_size[2] == 0) {
        md_vk_fail("kernel '%s': group_size is required ({%u,%u,%u} given); use the generated kernel descriptor",
                   label, desc->group_size[0], desc->group_size[1], desc->group_size[2]);
        return NULL;
    }
    uint32_t spv[3];
    if (md_vk_spirv_local_size((const uint32_t*)desc->code, desc->code_size / 4, spv) &&
        (spv[0] != desc->group_size[0] || spv[1] != desc->group_size[1] || spv[2] != desc->group_size[2])) {
        md_vk_fail("kernel '%s': group_size {%u,%u,%u} does not match [numthreads(%u,%u,%u)] in the shader",
                   label, desc->group_size[0], desc->group_size[1], desc->group_size[2], spv[0], spv[1], spv[2]);
        return NULL;
    }
    const uint64_t threads = (uint64_t)desc->group_size[0] * desc->group_size[1] * desc->group_size[2];
    if (threads > dev->props.limits.maxComputeWorkGroupInvocations) {
        md_vk_fail("kernel '%s': %llu threads per group exceeds the device limit of %u",
                   label, (unsigned long long)threads, dev->props.limits.maxComputeWorkGroupInvocations);
        return NULL;
    }

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) { md_vk_fail("out of memory"); return NULL; }
    memset(k, 0, sizeof(*k));
    k->device       = dev;
    k->group_size[0] = desc->group_size[0];
    k->group_size[1] = desc->group_size[1];
    k->group_size[2] = desc->group_size[2];
    k->args_size    = desc->args_size;
    snprintf(k->label, sizeof(k->label), "%s", label);

    VkShaderModuleCreateInfo smci = {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    smci.codeSize = desc->code_size;
    smci.pCode    = (const uint32_t*)desc->code;
    if (!md_vk_check(vkCreateShaderModule(dev->device, &smci, NULL, &k->module), "vkCreateShaderModule")) {
        md_vk_kernel_free(dev, k);
        return NULL;
    }

    VkComputePipelineCreateInfo cpci = {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    cpci.stage.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    cpci.stage.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    cpci.stage.module = k->module;
    cpci.stage.pName  = desc->entry_point ? desc->entry_point : "main";
    cpci.layout       = dev->pipeline_layout;
    if (!md_vk_check(vkCreateComputePipelines(dev->device, VK_NULL_HANDLE, 1, &cpci, NULL, &k->pipeline), "vkCreateComputePipelines")) {
        md_vk_kernel_free(dev, k);
        return NULL;
    }

    md_mutex_lock(&dev->device_mutex);
    md_gpu_kernel_t* slot = (md_gpu_kernel_t*)md_vk_vec_push(&dev->kernels, dev->alloc);
    if (slot) *slot = k;
    md_mutex_unlock(&dev->device_mutex);
    if (!slot) { md_vk_kernel_free(dev, k); md_vk_fail("out of memory"); return NULL; }
    return k;
}

void md_gpu_kernel_destroy(md_gpu_kernel_t k) {
    if (!k) return;
    md_gpu_device_t dev = k->device;
    md_mutex_lock(&dev->device_mutex);
    md_vk_vec_remove_ptr(&dev->kernels, k);
    md_vk_retire_locked(dev, MD_VK_RETIRE_KERNEL, k);
    md_mutex_unlock(&dev->device_mutex);
}

bool md_gpu_kernel_info(md_gpu_kernel_t k, md_gpu_kernel_info_t* info) {
    if (!k || !info) return false;
    memset(info, 0, sizeof(*info));
    info->group_size[0]            = k->group_size[0];
    info->group_size[1]            = k->group_size[1];
    info->group_size[2]            = k->group_size[2];
    info->args_size                = k->args_size;
    info->max_threads_per_group    = k->device->props.limits.maxComputeWorkGroupInvocations;
    info->preferred_group_multiple = k->device->subgroup_size;
    return true;
}

md_gpu_grid_t md_gpu_grid_for(md_gpu_kernel_t k, uint32_t nx, uint32_t ny, uint32_t nz) {
    md_gpu_grid_t g = {0, 0, 0};
    if (!k) return g;
    g.x = (uint32_t)(((uint64_t)nx + k->group_size[0] - 1) / k->group_size[0]);
    g.y = (uint32_t)(((uint64_t)ny + k->group_size[1] - 1) / k->group_size[1]);
    g.z = (uint32_t)(((uint64_t)nz + k->group_size[2] - 1) / k->group_size[2]);
    return g;
}

static VkCommandBuffer md_vk_launch_common(md_gpu_stream_t s, md_gpu_kernel_t k, const void* args, size_t args_size) {
    md_gpu_device_t dev = s->device;
    if (k->device != dev) { md_vk_fail("kernel '%s' belongs to a different device", k->label); return VK_NULL_HANDLE; }
    if (s->kind == MD_GPU_STREAM_TRANSFER) {
        md_vk_fail("kernel '%s' launched into '%s', a transfer stream; kernels run on compute streams", k->label, s->label);
        return VK_NULL_HANDLE;
    }
    if (k->args_size != 0 && args_size != k->args_size) {
        md_vk_fail("kernel '%s' expects a %u-byte argument struct but %zu bytes were passed", k->label, k->args_size, args_size);
        return VK_NULL_HANDLE;
    }
    if (args_size > 0 && !args) { md_vk_fail("kernel '%s': null args with non-zero size", k->label); return VK_NULL_HANDLE; }

    uint64_t arg_addr = 0;
    if (args_size > 0) {
        void* host;
        if (!md_vk_arena_alloc(s, args_size, &arg_addr, &host, NULL, NULL)) return VK_NULL_HANDLE;
        memcpy(host, args, args_size);
    }

    VkCommandBuffer cmd = md_vk_begin_op(s);
    if (!cmd) return VK_NULL_HANDLE;
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, k->pipeline);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, dev->pipeline_layout, 0, 1, &dev->desc_set, 0, NULL);
    vkCmdPushConstants(cmd, dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, 8, &arg_addr);
    return cmd;
}

bool md_gpu_launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid, const void* args, size_t args_size) {
    if (!s || !k) return md_vk_fail("md_gpu_launch: null stream or kernel");
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;   /* empty launch */
    const VkPhysicalDeviceLimits* lim = &s->device->props.limits;
    if (grid.x > lim->maxComputeWorkGroupCount[0] || grid.y > lim->maxComputeWorkGroupCount[1] ||
        grid.z > lim->maxComputeWorkGroupCount[2]) {
        return md_vk_fail("kernel '%s': grid {%u,%u,%u} exceeds the device limit {%u,%u,%u}", k->label,
                          grid.x, grid.y, grid.z, lim->maxComputeWorkGroupCount[0],
                          lim->maxComputeWorkGroupCount[1], lim->maxComputeWorkGroupCount[2]);
    }
    VkCommandBuffer cmd = md_vk_launch_common(s, k, args, args_size);
    if (!cmd) return false;
    vkCmdDispatch(cmd, grid.x, grid.y, grid.z);
    md_vk_end_op(s);
    return true;
}

bool md_gpu_launch_indirect(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_addr_t grid, const void* args, size_t args_size) {
    if (!s || !k || !grid) return md_vk_fail("md_gpu_launch_indirect: null argument");
    uint64_t off;
    md_vk_block_t* b = md_vk_resolve(s->device, grid, 3 * sizeof(uint32_t), &off, "md_gpu_launch_indirect");
    if (!b) return false;
    if (off % 4 != 0) return md_vk_fail("md_gpu_launch_indirect: grid address must be 4-byte aligned");
    VkCommandBuffer cmd = md_vk_launch_common(s, k, args, args_size);
    if (!cmd) return false;
    vkCmdDispatchIndirect(cmd, b->buffer, off);
    md_vk_end_op(s);
    return true;
}

/* Mirrors Args in src/shaders/gpu/md_gpu_make_grid.slang. */
typedef struct md_vk_make_grid_args_t {
    md_gpu_addr_t count;
    md_gpu_addr_t out_grid;
    md_gpu_uint4  local;     /* xyz = threads per group, w unused */
} md_vk_make_grid_args_t;

static bool md_vk_create_builtin_kernels(md_gpu_device_t dev) {
    md_gpu_kernel_desc_t d;
    memset(&d, 0, sizeof(d));
    d.code          = md_gpu_make_grid_spv;
    d.code_size     = sizeof(md_gpu_make_grid_spv);
    d.label         = "md_gpu make_grid";
    d.group_size[0] = 1;
    d.group_size[1] = 1;
    d.group_size[2] = 1;
    d.args_size     = (uint32_t)sizeof(md_vk_make_grid_args_t);
    dev->make_grid_kernel = md_gpu_kernel_create(dev, &d);
    if (!dev->make_grid_kernel) return false;
    /* Internal: owned by the device, not listed with the user's kernels. */
    md_mutex_lock(&dev->device_mutex);
    md_vk_vec_remove_ptr(&dev->kernels, dev->make_grid_kernel);
    md_mutex_unlock(&dev->device_mutex);
    return true;
}

bool md_gpu_make_grid(md_gpu_stream_t s, md_gpu_addr_t out_grid, md_gpu_addr_t count, md_gpu_kernel_t k) {
    if (!s || !out_grid || !count || !k) return md_vk_fail("md_gpu_make_grid: null argument");
    if (!md_vk_resolve(s->device, out_grid, 3 * sizeof(uint32_t), NULL, "md_gpu_make_grid (out_grid)")) return false;
    if (!md_vk_resolve(s->device, count, sizeof(uint32_t), NULL, "md_gpu_make_grid (count)")) return false;
    md_vk_make_grid_args_t a;
    memset(&a, 0, sizeof(a));
    a.count    = count;
    a.out_grid = out_grid;
    a.local.x  = k->group_size[0];
    a.local.y  = k->group_size[1];
    a.local.z  = k->group_size[2];
    return md_gpu_launch(s, s->device->make_grid_kernel, md_gpu_grid(1, 1, 1), &a, sizeof(a));
}

/* =========================================================================
   12. Host callbacks and polling
   ========================================================================= */

bool md_gpu_sync_on_complete(md_gpu_device_t dev, md_gpu_sync_t sync, md_gpu_host_fn fn, void* user) {
    if (!dev || !fn) return md_vk_fail("md_gpu_sync_on_complete: null argument");
    if (md_gpu_sync_is_valid(sync) && sync.stream->device != dev) return md_vk_fail("md_gpu_sync_on_complete: sync from another device");
    md_mutex_lock(&dev->device_mutex);
    md_vk_hostfn_t* h = (md_vk_hostfn_t*)md_vk_vec_push(&dev->hostfns, dev->alloc);
    if (h) { h->sync = sync; h->fn = fn; h->user = user; }
    md_mutex_unlock(&dev->device_mutex);
    return h ? true : md_vk_fail("out of memory");
}

bool md_gpu_launch_host_fn(md_gpu_stream_t s, md_gpu_host_fn fn, void* user) {
    if (!s || !fn) return md_vk_fail("md_gpu_launch_host_fn: null argument");
    md_gpu_sync_t sync = md_gpu_stream_record(s);
    return md_gpu_sync_on_complete(s->device, sync, fn, user);
}

/* Timeline values sampled once per poll pass.

   Callbacks must fire in the order they were registered. Deciding readiness by
   re-reading the semaphore for every entry breaks that: two callbacks
   registered against the *same* sync value can resolve differently, because
   the timeline may advance between the two reads, and then the later one
   overtakes the earlier. Snapshotting makes readiness a property of the pass
   rather than of the instant. Unbounded: every stream gets an entry. */
typedef struct md_vk_snapshot_t {
    md_gpu_stream_t stream;
    uint64_t        completed;
} md_vk_snapshot_t;

static bool md_vk_snapshot_complete(md_gpu_device_t dev, md_vk_vec_t* snap, md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    for (size_t i = 0; i < snap->count; ++i) {
        md_vk_snapshot_t* e = &MD_VK_VEC_AT(*snap, md_vk_snapshot_t, i);
        if (e->stream == sync.stream) return e->completed >= sync.value;
    }
    uint64_t completed = md_vk_stream_completed(sync.stream);
    md_vk_snapshot_t* e = (md_vk_snapshot_t*)md_vk_vec_push(snap, dev->alloc);
    if (e) { e->stream = sync.stream; e->completed = completed; }
    /* If the push failed, readiness falls back to this instant's reading; the
       ordering guarantee then degrades, but nothing fires early. */
    return completed >= sync.value;
}

uint32_t md_gpu_device_poll(md_gpu_device_t dev) {
    if (!dev) return 0;
    uint32_t fired = 0;
    md_vk_vec_t snap;
    md_vk_vec_init(&snap, sizeof(md_vk_snapshot_t));

    /* Callbacks run outside the lock, one at a time, in registration order.
       Each may register new callbacks or free memory; those are seen by the
       loop but judged against the same snapshot. */
    for (;;) {
        md_vk_hostfn_t ready;
        bool have = false;
        md_mutex_lock(&dev->device_mutex);
        for (size_t i = 0; i < dev->hostfns.count; ++i) {
            md_vk_hostfn_t* h = &MD_VK_VEC_AT(dev->hostfns, md_vk_hostfn_t, i);
            if (!md_vk_snapshot_complete(dev, &snap, h->sync)) continue;
            ready = *h;
            md_vk_vec_remove(&dev->hostfns, i);
            have = true;
            break;
        }
        md_mutex_unlock(&dev->device_mutex);
        if (!have) break;
        ready.fn(ready.user);
        fired++;
    }
    md_vk_vec_free(&snap, dev->alloc);

    md_mutex_lock(&dev->device_mutex);
    md_vk_process_retires_locked(dev, false);
    for (size_t i = 0; i < dev->pools.count; ++i) {
        md_gpu_pool_t p = MD_VK_VEC_AT(dev->pools, md_gpu_pool_t, i);
        if (p->cache_limit != 0) md_vk_pool_trim_locked(p, p->cache_limit);
    }
    md_mutex_unlock(&dev->device_mutex);
    return fired;
}

/* =========================================================================
   13. Device destruction
   ========================================================================= */

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;
    if (dev->device) {
        /* Flush and idle every stream, then fire what is pending. */
        for (size_t i = 0; i < dev->streams.count; ++i) {
            md_gpu_stream_t s = MD_VK_VEC_AT(dev->streams, md_gpu_stream_t, i);
            s->upload_open = false;
            md_vk_stream_submit(s);
        }
        vkDeviceWaitIdle(dev->device);
        md_gpu_device_poll(dev);

        /* Everything the caller did not destroy, the device does. */
        while (dev->pools.count > 0) md_gpu_pool_destroy(MD_VK_VEC_AT(dev->pools, md_gpu_pool_t, 0));
        while (dev->kernels.count > 0) md_gpu_kernel_destroy(MD_VK_VEC_AT(dev->kernels, md_gpu_kernel_t, 0));
        if (dev->make_grid_kernel) md_vk_kernel_free(dev, dev->make_grid_kernel);

        md_mutex_lock(&dev->device_mutex);
        md_vk_process_retires_locked(dev, true);
        md_mutex_unlock(&dev->device_mutex);

        for (size_t i = 0; i < dev->streams.count; ++i) {
            md_vk_stream_free(MD_VK_VEC_AT(dev->streams, md_gpu_stream_t, i));
        }
        md_vk_vec_free(&dev->streams, alloc);
        md_vk_vec_free(&dev->pools, alloc);
        md_vk_vec_free(&dev->kernels, alloc);

        for (uint32_t i = 0; i < dev->sampler_count; ++i) vkDestroySampler(dev->device, dev->samplers[i].sampler, NULL);
        if (dev->dummy_view)    vkDestroyImageView(dev->device, dev->dummy_view, NULL);
        if (dev->dummy_image)   vkDestroyImage(dev->device, dev->dummy_image, NULL);
        if (dev->dummy_mem)     vkFreeMemory(dev->device, dev->dummy_mem, NULL);
        if (dev->dummy_sampler) vkDestroySampler(dev->device, dev->dummy_sampler, NULL);
        md_vk_vec_free(&dev->hostfns,  alloc);
        md_vk_vec_free(&dev->retires,  alloc);
        md_vk_vec_free(&dev->registry, alloc);
        if (dev->pipeline_layout) vkDestroyPipelineLayout(dev->device, dev->pipeline_layout, NULL);
        if (dev->desc_pool)       vkDestroyDescriptorPool(dev->device, dev->desc_pool, NULL);
        if (dev->set_layout)      vkDestroyDescriptorSetLayout(dev->device, dev->set_layout, NULL);
        md_mutex_destroy(&dev->queue_mutex);
        md_mutex_destroy(&dev->device_mutex);
        vkDestroyDevice(dev->device, NULL);
    }
    if (dev->messenger && vkDestroyDebugUtilsMessengerEXT) {
        vkDestroyDebugUtilsMessengerEXT(dev->instance, dev->messenger, NULL);
    }
    if (dev->instance) vkDestroyInstance(dev->instance, NULL);
    md_free(alloc, dev, sizeof(*dev));
}

/*
md_gpu_vulkan.c — Vulkan backend for md_gpu.h

Structure of this file:

    1.  Configuration, error handling, small utilities
    2.  Types
    3.  Allocation registry (device address -> allocation lookup)
    4.  Transient arenas (argument blocks and staging)
    5.  Device creation
    6.  Bindless descriptor set
    7.  Streams, command buffers, submission, program-order barriers
    8.  Memory: pools, malloc/free, memcpy, memset, uploads
    9.  Textures and samplers
    10. Kernels and launches
    11. Graphs
    12. Host callbacks and polling

The dependency model is program order within a stream, implemented as a single
global VkMemoryBarrier2 between consecutive operations in a command buffer.
There is no per-resource state tracking anywhere in this file, and every image
lives in VK_IMAGE_LAYOUT_GENERAL for its entire life.
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

   Slang's `None` assignment is fixed; md_gpu declares only the bindings it
   can actually populate and leaves gaps for the rest. Adding, say, storage
   buffers to the heap later means declaring binding 7 here and nothing else.
   Keep MD_VK_BINDLESS_SPACE in sync with the -bindless-space-index flag in
   cmake/CompileGpuShaders.cmake, and this table in sync with the comment in
   src/shaders/md_gpu.slang. */
#define MD_VK_BINDLESS_SPACE         0u
#define MD_VK_BINDING_SAMPLER        0u   /* VK_DESCRIPTOR_TYPE_SAMPLER       */
#define MD_VK_BINDING_SAMPLED_IMAGE  2u   /* VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE */
#define MD_VK_BINDING_STORAGE_IMAGE  3u   /* VK_DESCRIPTOR_TYPE_STORAGE_IMAGE */
#define MD_VK_BINDING_COUNT          3u   /* how many of the above we declare */

#define MD_VK_MAX_QUEUES_PER_FAMILY  8u
#define MD_VK_MAX_TEXTURES        4096u
#define MD_VK_MAX_SAMPLERS         256u
#define MD_VK_ARENA_PAGE_SIZE   (256u * 1024u)
#define MD_VK_GRAPH_ARG_PAGE     (64u * 1024u)
#define MD_VK_ARG_ALIGN            64u
#define MD_VK_MAX_PENDING_WAITS     16u
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

static void md_vk_vec_free(md_vk_vec_t* v, struct md_allocator_i* alloc) {
    if (v->data) md_free(alloc, v->data, v->capacity * v->stride);
    v->data = NULL; v->count = 0; v->capacity = 0;
}

#define MD_VK_VEC_AT(v, type, i) (((type*)(v).data)[i])

/* =========================================================================
   2. Types
   ========================================================================= */

typedef struct md_vk_block_t {
    VkBuffer           buffer;
    VkDeviceMemory     memory;
    uint64_t           address;      /* device address of byte 0 */
    void*              host;         /* mapped pointer, or NULL  */
    uint64_t           capacity;     /* actual allocated size    */
    uint64_t           size;         /* size requested by the current user */
    md_gpu_mem_flags_t flags;
    md_gpu_pool_t      pool;
    bool               in_use;
    /* Stream-ordered free bookkeeping: the block becomes reusable at
       free_value on free_stream, or immediately for that same stream. */
    md_gpu_stream_t    free_stream;
    uint64_t           free_value;
} md_vk_block_t;

typedef struct md_gpu_pool {
    md_gpu_device_t device;
    md_gpu_mem_flags_t flags;        /* the one kind this pool serves */
    uint64_t        block_size;
    uint64_t        release_threshold;
    md_vk_vec_t     blocks;          /* md_vk_block_t*  */
    uint64_t        in_use_bytes;
    uint64_t        reserved_bytes;
    uint64_t        peak_in_use_bytes;
    uint64_t        alloc_count;
    uint64_t        reuse_count;
    char            label[64];
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

typedef struct md_gpu_graph {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    VkCommandPool        pool;
    VkCommandBuffer      cmd;
    md_vk_vec_t          arg_pages;    /* md_vk_page_t* — owned, persistent */
    md_vk_vec_t          launches;     /* md_vk_graph_launch_t */
    char                 label[64];
} md_gpu_graph;

typedef struct md_vk_graph_launch_t {
    size_t   page_index;
    uint64_t offset;
    size_t   size;
} md_vk_graph_launch_t;

typedef struct md_gpu_stream {
    md_gpu_device_t      device;
    md_gpu_stream_kind_t kind;
    uint32_t             family;
    VkQueue              queue;
    VkSemaphore          timeline;
    uint64_t             next_value;       /* value the next submit signals */
    uint64_t             submitted_value;  /* last value submitted          */

    VkCommandPool        cmd_pool;
    md_vk_vec_t          cmds;             /* md_vk_cmd_t */
    VkCommandBuffer      open;             /* currently recording, or NULL  */
    bool                 has_work;
    bool                 needs_barrier;
    bool                 auto_order;
    VkPipeline           bound_pipeline;

    VkSemaphore          wait_sems[MD_VK_MAX_PENDING_WAITS];
    uint64_t             wait_vals[MD_VK_MAX_PENDING_WAITS];
    uint32_t             wait_count;

    md_vk_arena_t        arena;

    /* Open upload, if any. */
    bool                 upload_open;
    bool                 upload_direct;
    md_gpu_ptr_t         upload_dst;
    uint64_t             upload_src_addr;
    size_t               upload_size;

    md_gpu_graph_t       capture;

    bool                 is_default;
    char                 label[64];
} md_gpu_stream;

typedef struct md_vk_tex_t {
    VkImage           image;
    VkImageView       view;
    VkDeviceMemory    memory;
    md_gpu_tex_desc_t desc;
    uint32_t          storage_slot;   /* heap slot, or UINT32_MAX */
    uint32_t          sampled_slot;   /* heap slot, or UINT32_MAX */
    uint32_t          generation;
    bool              live;
} md_vk_tex_t;

typedef struct md_vk_sampler_t {
    VkSampler sampler;
    uint32_t  generation;
    bool      live;
} md_vk_sampler_t;

typedef struct md_gpu_kernel {
    md_gpu_device_t device;
    VkShaderModule  module;
    VkPipeline      pipeline;
    uint32_t        group_size[3];
    char            label[64];
} md_gpu_kernel;

typedef struct md_vk_hostfn_t {
    md_gpu_sync_t  sync;
    md_gpu_host_fn fn;
    void*          user;
    bool           internal;
} md_vk_hostfn_t;

typedef struct md_vk_retire_t {
    VkImage        image;
    VkImageView    view;
    VkDeviceMemory memory;
    VkBuffer       buffer;
    uint32_t       storage_slot;   /* UINT32_MAX for none */
    uint32_t       sampled_slot;   /* UINT32_MAX for none */
    md_gpu_sync_t  waits[8];
    uint32_t       wait_count;
} md_vk_retire_t;

/* Everything the backend asks of a device, in one place, so that a device we
   cannot use is rejected by name instead of failing later inside
   vkCreateDevice with a bare VkResult. The optional flags are enabled only
   when the driver reports them, which is why they are carried out of the
   probe rather than re-queried. */
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

    uint32_t compute_family, transfer_family;
    /* Streams are spread round-robin over these. A device that exposes several
       compute queues (typical: 8 on NVIDIA) then gets real concurrency from
       "use another stream" rather than just permission to overlap. */
    VkQueue  compute_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t compute_queue_count;
    VkQueue  transfer_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t transfer_queue_count;
    uint32_t next_compute_queue, next_transfer_queue;
    md_mutex_t queue_mutex;      /* streams may share a VkQueue */
    md_mutex_t device_mutex;     /* pools, registry, textures, retire lists */

    /* Bindless */
    VkDescriptorSetLayout set_layout;
    VkDescriptorPool      desc_pool;
    VkDescriptorSet       desc_set;
    VkPipelineLayout      pipeline_layout;

    /* One 1x1(x1) placeholder per image class. Freed descriptor slots are
       pointed at these so that no descriptor ever references a destroyed
       view, which matters for bindless sets that outlive their contents. */
    VkImage         dummy_image;
    VkImageView     dummy_view;
    VkDeviceMemory  dummy_mem;
    VkSampler       dummy_sampler;

    md_vk_tex_t     textures[MD_VK_MAX_TEXTURES];
    uint32_t        tex_free[MD_VK_MAX_TEXTURES];
    uint32_t        tex_free_count;
    md_vk_sampler_t samplers[MD_VK_MAX_SAMPLERS];
    uint32_t        sampler_free[MD_VK_MAX_SAMPLERS];
    uint32_t        sampler_free_count;

    /* Address -> allocation registry, sorted ascending by address. */
    md_vk_vec_t     registry;   /* md_vk_block_t* */

    md_vk_vec_t     pools;      /* md_gpu_pool_t  */

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

/* A texture or sampler handle is exactly what Slang's DescriptorHandle<T>
   lowers to: a uint2 whose .x component indexes the heap. So the handle *is*
   the heap slot, zero-extended. Slot 0 is never handed out, which lets zero
   remain the null handle. Nothing else may live in the low 32 bits. */
static inline md_gpu_tex_t md_vk_tex_pack(uint32_t slot) { return (uint64_t)slot; }
static inline uint32_t md_vk_tex_slot(uint64_t h) { return (uint32_t)h; }
static inline bool     md_vk_tex_valid(uint64_t h)      { return md_vk_tex_slot(h) != 0 && (h >> 32) == 0; }

/* Process-wide device list, so that a bare device pointer or texture handle
   can be resolved back to the device that owns it. */
static md_gpu_device_t md_vk_devices[8];
static uint32_t        md_vk_device_count;

/* Forward declarations. */
static bool  md_vk_stream_ensure_cmd(md_gpu_stream_t s);
static void  md_vk_stream_barrier_if_needed(md_gpu_stream_t s);
static bool  md_vk_stream_submit(md_gpu_stream_t s);
static uint64_t md_vk_stream_completed(md_gpu_stream_t s);
static md_vk_block_t* md_vk_registry_find(md_gpu_device_t dev, uint64_t address);
static uint32_t md_vk_poll_internal(md_gpu_device_t dev, bool fire_user);
static bool  md_vk_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host);
typedef struct md_vk_tex_t md_vk_tex_t;
static md_gpu_device_t md_vk_tex_device(md_gpu_tex_t tex, md_vk_tex_t** out);
static bool  md_vk_immediate(md_gpu_device_t dev, void (*record)(VkCommandBuffer, void*), void* user);
static void  md_vk_record_to_general(VkCommandBuffer cmd, void* user);
static uint32_t md_vk_find_memory_type(md_gpu_device_t dev, uint32_t type_bits, VkMemoryPropertyFlags required, VkMemoryPropertyFlags preferred);

/* =========================================================================
   3. Allocation registry
   ========================================================================= */

/* Binary search for the block whose [address, address+capacity) contains a. */
static md_vk_block_t* md_vk_registry_find(md_gpu_device_t dev, uint64_t address) {
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

static bool md_vk_registry_insert(md_gpu_device_t dev, md_vk_block_t* blk) {
    if (!md_vk_vec_reserve(&dev->registry, dev->alloc, dev->registry.count + 1)) return false;
    md_vk_block_t** arr = (md_vk_block_t**)dev->registry.data;
    size_t i = dev->registry.count;
    while (i > 0 && arr[i - 1]->address > blk->address) {
        arr[i] = arr[i - 1];
        i--;
    }
    arr[i] = blk;
    dev->registry.count++;
    return true;
}

static void md_vk_registry_remove(md_gpu_device_t dev, md_vk_block_t* blk) {
    md_vk_block_t** arr = (md_vk_block_t**)dev->registry.data;
    for (size_t i = 0; i < dev->registry.count; ++i) {
        if (arr[i] == blk) {
            memmove(&arr[i], &arr[i + 1], (dev->registry.count - i - 1) * sizeof(*arr));
            dev->registry.count--;
            return;
        }
    }
}

/* =========================================================================
   4. Raw buffer helpers and transient arenas
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

static bool md_vk_create_raw_buffer(md_gpu_device_t dev, uint64_t size, md_gpu_mem_flags_t flags,
                                    VkBuffer* out_buf, VkDeviceMemory* out_mem, uint64_t* out_addr, void** out_host) {
    VkBufferCreateInfo bci = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bci.size  = size;
    bci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
              | VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT
              | VK_BUFFER_USAGE_TRANSFER_SRC_BIT
              | VK_BUFFER_USAGE_TRANSFER_DST_BIT
              | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT
              | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    bci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    VkBuffer buf = VK_NULL_HANDLE;
    if (!md_vk_check(vkCreateBuffer(dev->device, &bci, NULL, &buf), "vkCreateBuffer")) return false;

    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(dev->device, buf, &req);

    VkMemoryPropertyFlags required = 0, preferred = 0;
    if (flags & (MD_GPU_MEM_HOST_WRITE | MD_GPU_MEM_HOST_READ)) {
        required = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        if (flags & MD_GPU_MEM_HOST_READ) preferred = required | VK_MEMORY_PROPERTY_HOST_CACHED_BIT;
        else                              preferred = required | VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    } else {
        required = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    }

    uint32_t type = md_vk_find_memory_type(dev, req.memoryTypeBits, required, preferred);
    if (type == UINT32_MAX && required == VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT) {
        /* UMA parts may expose no purely device-local type. */
        type = md_vk_find_memory_type(dev, req.memoryTypeBits, 0, 0);
    }
    if (type == UINT32_MAX) {
        vkDestroyBuffer(dev->device, buf, NULL);
        return md_vk_fail("no suitable memory type for %llu bytes (flags 0x%x)", (unsigned long long)size, flags);
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
    if (flags & (MD_GPU_MEM_HOST_WRITE | MD_GPU_MEM_HOST_READ)) {
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

/* Reserve `size` bytes of transient, device-addressable, host-writable memory
   from the stream's arena. Valid until the submission that consumes it
   completes. */
static bool md_vk_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host) {
    md_gpu_device_t dev = s->device;
    md_vk_arena_t*  a   = &s->arena;
    uint64_t need = md_vk_align_up(size, MD_VK_ARG_ALIGN);

    /* Try the current page. */
    if (a->pages.count > 0) {
        md_vk_page_t* p = MD_VK_VEC_AT(a->pages, md_vk_page_t*, a->current);
        if (p->cursor + need <= p->capacity) {
            *out_addr = p->address + p->cursor;
            *out_host = p->host + p->cursor;
            p->cursor += need;
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
            *out_addr = p->address + p->cursor;
            *out_host = p->host + p->cursor;
            p->cursor += need;
            return true;
        }
    }

    /* Allocate a new page. */
    uint64_t page_size = MD_VK_ARENA_PAGE_SIZE;
    while (page_size < need) page_size *= 2;
    md_vk_page_t* p = md_vk_page_create(dev, page_size);
    if (!p) return md_vk_fail("failed to allocate a %llu byte transient page", (unsigned long long)page_size);
    md_vk_page_t** slot = (md_vk_page_t**)md_vk_vec_push(&a->pages, dev->alloc);
    if (!slot) { md_vk_page_destroy(dev, p); return false; }
    *slot = p;
    a->current = a->pages.count - 1;
    *out_addr = p->address;
    *out_host = p->host;
    p->cursor = need;
    return true;
}

/* Tag every in-flight arena page with the value that releases it. */
/* Stamp every page holding data with the value that releases it.

   This must overwrite an existing stamp, not skip it. A page can be appended
   to across several submissions -- the fast path in md_vk_arena_alloc happily
   carves more out of the current page after it has already been submitted,
   which is safe in itself. But if the stamp were left at the *first*
   submission's value, the page would be recycled as soon as that one
   completed, while a later submission was still reading its argument blocks
   out of the same page. Always taking the newest value keeps the page alive
   until everything that touched it has finished. */
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
   5. Device creation
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

    /* Register before anything that may need an address -> device lookup. */
    if (md_vk_device_count >= (uint32_t)(sizeof(md_vk_devices) / sizeof(md_vk_devices[0]))) {
        md_vk_fail("too many simultaneous md_gpu devices");
        goto fail_device;
    }
    md_vk_devices[md_vk_device_count++] = dev;

    /* Texture / sampler free lists, allocated high-to-low so indices start at 0. */
    /* Slot 0 is never handed out so that a zero handle stays the null handle. */
    for (uint32_t i = 0; i < MD_VK_MAX_TEXTURES - 1; ++i) dev->tex_free[i] = MD_VK_MAX_TEXTURES - 1 - i;
    dev->tex_free_count = MD_VK_MAX_TEXTURES - 1;
    for (uint32_t i = 0; i < MD_VK_MAX_SAMPLERS - 1; ++i) dev->sampler_free[i] = MD_VK_MAX_SAMPLERS - 1 - i;
    dev->sampler_free_count = MD_VK_MAX_SAMPLERS - 1;

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
    info->preferred_group_multiple = 32;
    info->timestamp_period_ns_num  = (uint64_t)(dev->props.limits.timestampPeriod * 1000.0f);
    info->timestamp_period_ns_den  = 1000;
    snprintf(info->name, sizeof(info->name), "%s", dev->props.deviceName);
    return true;
}

/* =========================================================================
   6. Bindless descriptor set
   ========================================================================= */

/* Class id -> descriptor type. Must match md_gpu.slang. */
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
    bindings[1].descriptorCount = MD_VK_MAX_TEXTURES;
    bindings[2].binding         = MD_VK_BINDING_STORAGE_IMAGE;
    bindings[2].descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings[2].descriptorCount = MD_VK_MAX_TEXTURES;

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
    sizes[1] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, MD_VK_MAX_TEXTURES};
    sizes[2] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, MD_VK_MAX_TEXTURES};

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

/* =========================================================================
   7. Streams
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
    s->auto_order = true;
    s->next_value = 1;
    s->is_default = is_default;
    snprintf(s->label, sizeof(s->label), "%s", label ? label : "stream");

    md_mutex_lock(&dev->device_mutex);
    if (kind == MD_GPU_STREAM_TRANSFER) {
        s->family = dev->transfer_family;
        s->queue  = dev->transfer_queues[dev->next_transfer_queue % dev->transfer_queue_count];
        dev->next_transfer_queue++;
    } else {
        s->family = dev->compute_family;
        s->queue  = dev->compute_queues[dev->next_compute_queue % dev->compute_queue_count];
        dev->next_compute_queue++;
    }
    md_mutex_unlock(&dev->device_mutex);

    md_vk_vec_init(&s->cmds, sizeof(md_vk_cmd_t));
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

    return s;
}

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label) {
    if (!dev) return NULL;
    return md_vk_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    return kind == MD_GPU_STREAM_TRANSFER ? dev->default_transfer : dev->default_compute;
}

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t s) { return s ? s->device : NULL; }

void md_gpu_stream_set_auto_order(md_gpu_stream_t s, bool enabled) { if (s) s->auto_order = enabled; }

/* Acquire a command buffer, recycling any whose submission has completed. */
static bool md_vk_stream_ensure_cmd(md_gpu_stream_t s) {
    if (s->open) return true;
    md_gpu_device_t dev = s->device;

    uint64_t done = md_vk_stream_completed(s);
    md_vk_cmd_t* found = NULL;
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &((md_vk_cmd_t*)s->cmds.data)[i];
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
        if (!slot) { vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &cb); return false; }
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
       one barrier at the top of the buffer closes it. */
    s->needs_barrier  = (s->submitted_value > 0);
    s->bound_pipeline = VK_NULL_HANDLE;
    return true;
}

/* The command buffer that operations should record into: the graph's during
   capture, the stream's otherwise. */
static VkCommandBuffer md_vk_target_cmd(md_gpu_stream_t s) {
    if (s->capture) return s->capture->cmd;
    if (!md_vk_stream_ensure_cmd(s)) return VK_NULL_HANDLE;
    return s->open;
}

/* The single global barrier that implements program order. */
static void md_vk_stream_barrier_if_needed(md_gpu_stream_t s) {
    if (!s->needs_barrier) return;
    if (!s->auto_order) return;

    VkCommandBuffer cmd = s->capture ? s->capture->cmd : s->open;
    if (!cmd) return;

    VkMemoryBarrier2 mb = {VK_STRUCTURE_TYPE_MEMORY_BARRIER_2};
    mb.srcStageMask  = VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT;
    mb.srcAccessMask = VK_ACCESS_2_SHADER_WRITE_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT;
    mb.dstStageMask  = VK_PIPELINE_STAGE_2_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_2_ALL_TRANSFER_BIT
                     | VK_PIPELINE_STAGE_2_DRAW_INDIRECT_BIT;
    mb.dstAccessMask = VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT
                     | VK_ACCESS_2_TRANSFER_READ_BIT | VK_ACCESS_2_TRANSFER_WRITE_BIT
                     | VK_ACCESS_2_INDIRECT_COMMAND_READ_BIT;

    VkDependencyInfo di = {VK_STRUCTURE_TYPE_DEPENDENCY_INFO};
    di.memoryBarrierCount = 1;
    di.pMemoryBarriers    = &mb;
    vkCmdPipelineBarrier2(cmd, &di);

    s->needs_barrier = false;
}

void md_gpu_stream_barrier(md_gpu_stream_t s) {
    if (!s) return;
    bool saved_auto = s->auto_order;
    s->auto_order    = true;
    s->needs_barrier = true;
    md_vk_stream_barrier_if_needed(s);
    s->auto_order    = saved_auto;
}

/* Mark that an operation was recorded. */
static void md_vk_stream_did_op(md_gpu_stream_t s) {
    s->needs_barrier = true;
    if (!s->capture) s->has_work = true;
}

static bool md_vk_stream_submit(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    if (!s->open || !s->has_work) {
        /* Nothing recorded. Still apply any pending waits at the next submit. */
        return true;
    }

    if (!md_vk_check(vkEndCommandBuffer(s->open), "vkEndCommandBuffer")) return false;

    uint64_t signal_value = s->next_value;

    VkCommandBufferSubmitInfo cbsi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO};
    cbsi.commandBuffer = s->open;

    VkSemaphoreSubmitInfo waits[MD_VK_MAX_PENDING_WAITS];
    for (uint32_t i = 0; i < s->wait_count; ++i) {
        waits[i] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        waits[i].semaphore = s->wait_sems[i];
        waits[i].value     = s->wait_vals[i];
        waits[i].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }

    VkSemaphoreSubmitInfo signal = {VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
    signal.semaphore = s->timeline;
    signal.value     = signal_value;
    signal.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;

    VkSubmitInfo2 si = {VK_STRUCTURE_TYPE_SUBMIT_INFO_2};
    si.waitSemaphoreInfoCount   = s->wait_count;
    si.pWaitSemaphoreInfos      = waits;
    si.commandBufferInfoCount   = 1;
    si.pCommandBufferInfos      = &cbsi;
    si.signalSemaphoreInfoCount = 1;
    si.pSignalSemaphoreInfos    = &signal;

    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueueSubmit2(s->queue, 1, &si, VK_NULL_HANDLE);
    md_mutex_unlock(&dev->queue_mutex);
    if (!md_vk_check(r, "vkQueueSubmit2")) return false;

    /* Tag the command buffer and arena pages with the value that frees them. */
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &((md_vk_cmd_t*)s->cmds.data)[i];
        if (c->cmd == s->open) { c->value = signal_value; c->pending = true; break; }
    }
    md_vk_arena_retire(s, signal_value);

    s->submitted_value = signal_value;
    s->next_value      = signal_value + 1;
    s->wait_count      = 0;
    s->open            = VK_NULL_HANDLE;
    s->has_work        = false;
    /* md_vk_stream_ensure_cmd re-arms this for the next command buffer. */
    s->needs_barrier   = false;
    s->bound_pipeline  = VK_NULL_HANDLE;
    return true;
}

md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t s) {
    md_gpu_sync_t out = md_gpu_sync_none();
    if (!s) return out;
    if (s->capture) { md_vk_fail("md_gpu_stream_record is not allowed during capture"); return out; }
    md_vk_stream_submit(s);
    if (s->submitted_value == 0) return out;
    out.stream = s;
    out.value  = s->submitted_value;
    return out;
}

void md_gpu_stream_wait(md_gpu_stream_t s, md_gpu_sync_t sync) {
    if (!s || !md_gpu_sync_is_valid(sync)) return;
    if (sync.stream == s) return;                    /* already ordered */
    if (md_gpu_sync_is_complete(sync)) return;       /* nothing to wait for */

    /* Work already issued must not be retroactively delayed: close it first. */
    if (s->has_work) md_vk_stream_submit(s);

    if (s->wait_count >= MD_VK_MAX_PENDING_WAITS) {
        md_vk_fail("more than %u pending waits on stream '%s'", MD_VK_MAX_PENDING_WAITS, s->label);
        return;
    }
    /* Collapse duplicate waits on the same semaphore, keeping the larger value. */
    for (uint32_t i = 0; i < s->wait_count; ++i) {
        if (s->wait_sems[i] == sync.stream->timeline) {
            if (sync.value > s->wait_vals[i]) s->wait_vals[i] = sync.value;
            return;
        }
    }
    s->wait_sems[s->wait_count] = sync.stream->timeline;
    s->wait_vals[s->wait_count] = sync.value;
    s->wait_count++;
}

void md_gpu_stream_flush(md_gpu_stream_t s) {
    if (!s || s->capture) return;
    md_vk_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    if (s->capture) { md_vk_fail("md_gpu_stream_sync is not allowed during capture"); return; }
    md_vk_stream_submit(s);
    if (s->submitted_value > 0) {
        VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
        wi.semaphoreCount = 1;
        wi.pSemaphores    = &s->timeline;
        wi.pValues        = &s->submitted_value;
        vkWaitSemaphores(s->device->device, &wi, UINT64_MAX);
    }
    /* Complete any internal work (device-to-host copies) that is now ready. */
    md_vk_poll_internal(s->device, false);
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

static void md_vk_stream_destroy_internal(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    md_gpu_stream_sync(s);
    if (s->open) {
        vkEndCommandBuffer(s->open);
        s->open = VK_NULL_HANDLE;
    }
    for (size_t i = 0; i < s->cmds.count; ++i) {
        md_vk_cmd_t* c = &((md_vk_cmd_t*)s->cmds.data)[i];
        vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &c->cmd);
    }
    md_vk_vec_free(&s->cmds, dev->alloc);
    md_vk_arena_free(dev, &s->arena);
    vkDestroyCommandPool(dev->device, s->cmd_pool, NULL);
    vkDestroySemaphore(dev->device, s->timeline, NULL);
    md_free(dev->alloc, s, sizeof(*s));
}

void md_gpu_stream_destroy(md_gpu_stream_t s) {
    if (!s || s->is_default) return;
    md_gpu_device_t dev = s->device;
    md_mutex_lock(&dev->device_mutex);
    md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
    for (size_t i = 0; i < dev->streams.count; ++i) {
        if (arr[i] == s) {
            memmove(&arr[i], &arr[i + 1], (dev->streams.count - i - 1) * sizeof(*arr));
            dev->streams.count--;
            break;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    md_vk_stream_destroy_internal(s);
}

/* =========================================================================
   8. Memory
   ========================================================================= */

md_gpu_pool_t md_gpu_pool_create(md_gpu_device_t dev, const md_gpu_pool_desc_t* desc) {
    if (!dev || !desc) return NULL;
    const uint64_t release_threshold = desc->release_threshold;
    const char* label = desc->label;
    md_gpu_pool_t p = (md_gpu_pool_t)md_alloc(dev->alloc, sizeof(md_gpu_pool));
    if (!p) { md_vk_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device            = dev;
    p->release_threshold = release_threshold;
    p->flags             = desc->flags;
    p->block_size        = desc->block_size;
    snprintf(p->label, sizeof(p->label), "%s", label ? label : "pool");
    md_vk_vec_init(&p->blocks, sizeof(md_vk_block_t*));

    md_mutex_lock(&dev->device_mutex);
    md_gpu_pool_t* slot = (md_gpu_pool_t*)md_vk_vec_push(&dev->pools, dev->alloc);
    if (slot) *slot = p;
    md_mutex_unlock(&dev->device_mutex);
    return p;
}

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
    md_vk_block_t** arr = (md_vk_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (arr[i]->in_use) out->blocks_in_use++;
        else                out->blocks_cached++;
    }
    md_mutex_unlock(&p->device->device_mutex);
}

/* Mark a block free at the current point in `stream`. Shared by
   md_gpu_free and md_gpu_pool_reset. Caller holds device_mutex. */
static void md_vk_block_release(md_vk_block_t* b, md_gpu_stream_t stream) {
    b->in_use      = false;
    b->free_stream = stream;
    b->free_value  = stream->has_work ? stream->next_value : stream->submitted_value;
    if (b->free_value == 0) b->free_stream = NULL;
    b->pool->in_use_bytes -= b->capacity;
}

void md_gpu_pool_reset(md_gpu_pool_t p, md_gpu_stream_t stream) {
    if (!p || !stream) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    md_vk_block_t** arr = (md_vk_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        if (arr[i]->in_use) md_vk_block_release(arr[i], stream);
    }
    md_mutex_unlock(&dev->device_mutex);
}

/* Caller holds device_mutex. */
static void md_vk_block_destroy(md_gpu_pool_t p, md_vk_block_t* b) {
    md_gpu_device_t dev = p->device;
    md_vk_registry_remove(dev, b);
    md_vk_destroy_raw_buffer(dev, b->buffer, b->memory);
    p->reserved_bytes -= b->capacity;
    md_free(dev->alloc, b, sizeof(*b));
}

static void md_vk_pool_trim_locked(md_gpu_pool_t p, uint64_t keep_bytes) {
    md_vk_block_t** arr = (md_vk_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count;) {
        md_vk_block_t* b = arr[i];
        bool free_and_idle = !b->in_use &&
            (b->free_stream == NULL || md_vk_stream_completed(b->free_stream) >= b->free_value);
        if (free_and_idle && p->reserved_bytes - p->in_use_bytes > keep_bytes) {
            md_vk_block_destroy(p, b);
            memmove(&arr[i], &arr[i + 1], (p->blocks.count - i - 1) * sizeof(*arr));
            p->blocks.count--;
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

void md_gpu_pool_destroy(md_gpu_pool_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    vkDeviceWaitIdle(dev->device);
    md_mutex_lock(&dev->device_mutex);
    md_vk_block_t** arr = (md_vk_block_t**)p->blocks.data;
    for (size_t i = 0; i < p->blocks.count; ++i) md_vk_block_destroy(p, arr[i]);
    md_vk_vec_free(&p->blocks, dev->alloc);

    md_gpu_pool_t* parr = (md_gpu_pool_t*)dev->pools.data;
    for (size_t i = 0; i < dev->pools.count; ++i) {
        if (parr[i] == p) {
            memmove(&parr[i], &parr[i + 1], (dev->pools.count - i - 1) * sizeof(*parr));
            dev->pools.count--;
            break;
        }
    }
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, p, sizeof(*p));
}

md_gpu_ptr_t md_gpu_malloc(md_gpu_pool_t p, size_t size, md_gpu_stream_t stream) {
    if (!p || size == 0) return NULL;
    const md_gpu_mem_flags_t flags = p->flags;
    md_gpu_device_t dev = p->device;
    uint64_t capacity = md_vk_next_pow2(size);

    md_mutex_lock(&dev->device_mutex);

    /* Reuse: same flags, big enough, and either freed on this stream (program
       order makes it safe immediately) or its free point has completed. */
    md_vk_block_t** arr = (md_vk_block_t**)p->blocks.data;
    md_vk_block_t*  best = NULL;
    for (size_t i = 0; i < p->blocks.count; ++i) {
        md_vk_block_t* b = arr[i];
        if (b->in_use || b->flags != flags || b->capacity < size) continue;
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
        return (md_gpu_ptr_t)(uintptr_t)best->address;
    }

    md_vk_block_t* b = (md_vk_block_t*)md_alloc(dev->alloc, sizeof(md_vk_block_t));
    if (!b) { md_mutex_unlock(&dev->device_mutex); md_vk_fail("out of memory"); return NULL; }
    memset(b, 0, sizeof(*b));
    if (!md_vk_create_raw_buffer(dev, capacity, flags, &b->buffer, &b->memory, &b->address, &b->host)) {
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        return NULL;
    }
    b->capacity = capacity;
    b->size     = size;
    b->flags    = flags;
    b->pool     = p;
    b->in_use   = true;

    md_vk_block_t** slot = (md_vk_block_t**)md_vk_vec_push(&p->blocks, dev->alloc);
    if (!slot || !md_vk_registry_insert(dev, b)) {
        md_vk_destroy_raw_buffer(dev, b->buffer, b->memory);
        md_free(dev->alloc, b, sizeof(*b));
        md_mutex_unlock(&dev->device_mutex);
        md_vk_fail("out of memory");
        return NULL;
    }
    *slot = b;
    p->reserved_bytes += capacity;
    p->in_use_bytes   += capacity;
    if (p->in_use_bytes > p->peak_in_use_bytes) p->peak_in_use_bytes = p->in_use_bytes;
    p->alloc_count++;

    md_mutex_unlock(&dev->device_mutex);
    return (md_gpu_ptr_t)(uintptr_t)b->address;
}

void md_gpu_free(md_gpu_ptr_t ptr, md_gpu_stream_t stream) {
    if (!ptr) return;
    /* Any live stream tells us the device; prefer the one given. */
    if (!stream) return;
    md_gpu_device_t dev = stream->device;

    md_mutex_lock(&dev->device_mutex);
    md_vk_block_t* b = md_vk_registry_find(dev, (uint64_t)(uintptr_t)ptr);
    if (!b || !b->in_use) {
        md_mutex_unlock(&dev->device_mutex);
        return;
    }
    /* The block becomes reusable once everything issued into `stream` up to
       now has completed. next_value is what the next submit will signal. */
    md_vk_block_release(b, stream);

    if (b->pool->release_threshold == 0) md_vk_pool_trim_locked(b->pool, 0);
    md_mutex_unlock(&dev->device_mutex);
}

/* Resolve a bare device pointer without the caller supplying a device. */
static md_vk_block_t* md_vk_find_anywhere(uint64_t addr, md_gpu_device_t* out_dev) {
    for (uint32_t i = 0; i < md_vk_device_count; ++i) {
        md_gpu_device_t d = md_vk_devices[i];
        if (!d) continue;
        md_vk_block_t* b = md_vk_registry_find(d, addr);
        if (b) { if (out_dev) *out_dev = d; return b; }
    }
    return NULL;
}

void* md_gpu_host_ptr(md_gpu_ptr_t ptr) {
    if (!ptr) return NULL;
    uint64_t addr = (uint64_t)(uintptr_t)ptr;
    md_vk_block_t* b = md_vk_find_anywhere(addr, NULL);
    if (!b || !b->host) return NULL;
    return (uint8_t*)b->host + (addr - b->address);
}

size_t md_gpu_ptr_size(md_gpu_ptr_t ptr) {
    if (!ptr) return 0;
    uint64_t addr = (uint64_t)(uintptr_t)ptr;
    md_vk_block_t* b = md_vk_find_anywhere(addr, NULL);
    if (!b || !b->in_use) return 0;
    uint64_t off = addr - b->address;
    return (size_t)(b->size > off ? b->size - off : 0);
}

md_gpu_ptr_t md_gpu_ptr_base(md_gpu_ptr_t ptr) {
    if (!ptr) return NULL;
    md_vk_block_t* b = md_vk_find_anywhere((uint64_t)(uintptr_t)ptr, NULL);
    return b ? (md_gpu_ptr_t)(uintptr_t)b->address : NULL;
}

/* Record a buffer-to-buffer copy in device address space. */
static bool md_vk_record_copy(md_gpu_stream_t s, md_vk_block_t* src, uint64_t src_off,
                              md_vk_block_t* dst, uint64_t dst_off, uint64_t size) {
    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);
    VkBufferCopy region = {src_off, dst_off, size};
    vkCmdCopyBuffer(cmd, src->buffer, dst->buffer, 1, &region);
    md_vk_stream_did_op(s);
    return true;
}

/* Copy from a transient arena address (staging) into a destination block. */
static bool md_vk_record_copy_from_arena(md_gpu_stream_t s, uint64_t src_addr,
                                         md_vk_block_t* dst, uint64_t dst_off, uint64_t size) {
    md_gpu_device_t dev = s->device;
    md_vk_page_t* page = NULL;
    uint64_t page_off = 0;
    md_vk_arena_t* a = s->capture ? NULL : &s->arena;
    md_vk_vec_t* pages = a ? &a->pages : &s->capture->arg_pages;
    for (size_t i = 0; i < pages->count; ++i) {
        md_vk_page_t* p = ((md_vk_page_t**)pages->data)[i];
        if (src_addr >= p->address && src_addr < p->address + p->capacity) {
            page = p; page_off = src_addr - p->address; break;
        }
    }
    if (!page) return md_vk_fail("staging address not found in arena");
    (void)dev;

    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);
    VkBufferCopy region = {page_off, dst_off, size};
    vkCmdCopyBuffer(cmd, page->buffer, dst->buffer, 1, &region);
    md_vk_stream_did_op(s);
    return true;
}

typedef struct md_vk_d2h_t {
    md_gpu_device_t dev;
    void*           dst;
    const void*     staging_host;
    size_t          size;
} md_vk_d2h_t;

static void md_vk_d2h_finish(void* user) {
    md_vk_d2h_t* c = (md_vk_d2h_t*)user;
    memcpy(c->dst, c->staging_host, c->size);
    md_free(c->dev->alloc, c, sizeof(*c));
}

bool md_gpu_memcpy_async(void* dst, const void* src, size_t size, md_gpu_stream_t s) {
    if (!s) return md_vk_fail("md_gpu_memcpy_async requires a stream");
    if (size == 0) return true;
    md_gpu_device_t dev = s->device;

    md_vk_block_t* dblk = md_vk_registry_find(dev, (uint64_t)(uintptr_t)dst);
    md_vk_block_t* sblk = md_vk_registry_find(dev, (uint64_t)(uintptr_t)src);

    if (dblk && sblk) {
        return md_vk_record_copy(s, sblk, (uint64_t)(uintptr_t)src - sblk->address,
                                    dblk, (uint64_t)(uintptr_t)dst - dblk->address, size);
    }

    if (dblk && !sblk) {
        /* Host to device. */
        uint64_t dst_off = (uint64_t)(uintptr_t)dst - dblk->address;
        /* Fast path: host-writable destination with nothing in flight. */
        if (dblk->host && !s->capture && !s->has_work &&
            md_vk_stream_completed(s) >= s->submitted_value) {
            memcpy((uint8_t*)dblk->host + dst_off, src, size);
            return true;
        }
        uint64_t addr; void* host;
        if (!md_vk_arena_alloc(s, size, &addr, &host)) return false;
        memcpy(host, src, size);
        return md_vk_record_copy_from_arena(s, addr, dblk, dst_off, size);
    }

    if (!dblk && sblk) {
        /* Device to host, asynchronous: stage, then memcpy on completion. */
        uint64_t src_off = (uint64_t)(uintptr_t)src - sblk->address;
        if (sblk->host && !s->has_work && md_vk_stream_completed(s) >= s->submitted_value) {
            memcpy(dst, (uint8_t*)sblk->host + src_off, size);
            return true;
        }
        uint64_t addr; void* host;
        if (!md_vk_arena_alloc(s, size, &addr, &host)) return false;

        /* Find the staging page's VkBuffer and record the copy into it. */
        md_vk_page_t* page = NULL; uint64_t page_off = 0;
        md_vk_vec_t* pages = s->capture ? &s->capture->arg_pages : &s->arena.pages;
        for (size_t i = 0; i < pages->count; ++i) {
            md_vk_page_t* p = ((md_vk_page_t**)pages->data)[i];
            if (addr >= p->address && addr < p->address + p->capacity) { page = p; page_off = addr - p->address; break; }
        }
        if (!page) return md_vk_fail("staging address not found");

        VkCommandBuffer cmd = md_vk_target_cmd(s);
        if (!cmd) return false;
        md_vk_stream_barrier_if_needed(s);
        VkBufferCopy region = {src_off, page_off, size};
        vkCmdCopyBuffer(cmd, sblk->buffer, page->buffer, 1, &region);
        md_vk_stream_did_op(s);

        md_vk_d2h_t* c = (md_vk_d2h_t*)md_alloc(dev->alloc, sizeof(md_vk_d2h_t));
        if (!c) return md_vk_fail("out of memory");
        c->dev = dev; c->dst = dst; c->staging_host = host; c->size = size;

        md_gpu_sync_t sync = md_gpu_stream_record(s);
        md_mutex_lock(&dev->device_mutex);
        md_vk_hostfn_t* h = (md_vk_hostfn_t*)md_vk_vec_push(&dev->hostfns, dev->alloc);
        if (h) { h->sync = sync; h->fn = md_vk_d2h_finish; h->user = c; h->internal = true; }
        md_mutex_unlock(&dev->device_mutex);
        return h != NULL;
    }

    /* Host to host. */
    memcpy(dst, src, size);
    return true;
}

bool md_gpu_memset_async(md_gpu_ptr_t dst, uint8_t value, size_t size, md_gpu_stream_t s) {
    if (!s) return md_vk_fail("md_gpu_memset_async requires a stream");
    if (!dst || size == 0) return true;
    md_gpu_device_t dev = s->device;
    md_vk_block_t* b = md_vk_registry_find(dev, (uint64_t)(uintptr_t)dst);
    if (!b) return md_vk_fail("md_gpu_memset_async: not a device pointer");

    uint64_t off   = (uint64_t)(uintptr_t)dst - b->address;
    uint64_t begin = off;
    uint64_t end   = off + size;
    uint64_t abeg  = md_vk_align_up(begin, 4);
    uint64_t aend  = end & ~3ull;

    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);

    if (aend > abeg) {
        uint32_t word = ((uint32_t)value) * 0x01010101u;
        vkCmdFillBuffer(cmd, b->buffer, abeg, aend - abeg, word);
        md_vk_stream_did_op(s);
    }

    /* Unaligned head and tail go through a small staged copy. */
    uint64_t head = (abeg > begin) ? (abeg - begin) : 0;
    uint64_t tail = (end > aend && aend >= abeg) ? (end - aend) : 0;
    if (aend <= abeg) { head = size; tail = 0; abeg = begin; }

    if (head) {
        uint64_t addr; void* host;
        if (!md_vk_arena_alloc(s, (size_t)head, &addr, &host)) return false;
        memset(host, value, (size_t)head);
        if (!md_vk_record_copy_from_arena(s, addr, b, begin, head)) return false;
    }
    if (tail) {
        uint64_t addr; void* host;
        if (!md_vk_arena_alloc(s, (size_t)tail, &addr, &host)) return false;
        memset(host, value, (size_t)tail);
        if (!md_vk_record_copy_from_arena(s, addr, b, aend, tail)) return false;
    }
    return true;
}

void* md_gpu_upload_begin(md_gpu_stream_t s, md_gpu_ptr_t dst, size_t size) {
    if (!s || !dst || size == 0) return NULL;
    if (s->upload_open) { md_vk_fail("an upload is already open on stream '%s'", s->label); return NULL; }
    md_gpu_device_t dev = s->device;
    md_vk_block_t* b = md_vk_registry_find(dev, (uint64_t)(uintptr_t)dst);
    if (!b) { md_vk_fail("md_gpu_upload_begin: not a device pointer"); return NULL; }

    uint64_t dst_off = (uint64_t)(uintptr_t)dst - b->address;

    /* Write straight into the destination when that cannot race the GPU. */
    if (b->host && !s->capture && !s->has_work && md_vk_stream_completed(s) >= s->submitted_value) {
        s->upload_open   = true;
        s->upload_direct = true;
        s->upload_dst    = dst;
        s->upload_size   = size;
        return (uint8_t*)b->host + dst_off;
    }

    uint64_t addr; void* host;
    if (!md_vk_arena_alloc(s, size, &addr, &host)) return NULL;
    s->upload_open     = true;
    s->upload_direct   = false;
    s->upload_dst      = dst;
    s->upload_src_addr = addr;
    s->upload_size     = size;
    return host;
}

bool md_gpu_upload_end(md_gpu_stream_t s) {
    if (!s || !s->upload_open) return md_vk_fail("no upload is open");
    s->upload_open = false;
    if (s->upload_direct) return true;
    md_vk_block_t* b = md_vk_registry_find(s->device, (uint64_t)(uintptr_t)s->upload_dst);
    if (!b) return md_vk_fail("upload destination is no longer live");
    uint64_t dst_off = (uint64_t)(uintptr_t)s->upload_dst - b->address;
    return md_vk_record_copy_from_arena(s, s->upload_src_addr, b, dst_off, s->upload_size);
}

/* =========================================================================
   9. Textures and samplers
   ========================================================================= */

typedef struct md_vk_format_info_t {
    VkFormat format;
    uint32_t bytes;
} md_vk_format_info_t;

static md_vk_format_info_t md_vk_format_info(md_gpu_format_t f) {
    md_vk_format_info_t i = {VK_FORMAT_UNDEFINED, 0};
    switch (f) {
    case MD_GPU_FORMAT_R32_FLOAT:    i.format = VK_FORMAT_R32_SFLOAT;          i.bytes = 4;  break;
    case MD_GPU_FORMAT_R32_UINT:     i.format = VK_FORMAT_R32_UINT;            i.bytes = 4;  break;
    case MD_GPU_FORMAT_RGBA32_FLOAT: i.format = VK_FORMAT_R32G32B32A32_SFLOAT; i.bytes = 16; break;
    case MD_GPU_FORMAT_RGBA8_UNORM:  i.format = VK_FORMAT_R8G8B8A8_UNORM;      i.bytes = 4;  break;
    default: break;
    }
    return i;
}

/* One-shot command submission, used only for the initial layout transition. */
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

md_gpu_tex_t md_gpu_tex_create(md_gpu_device_t dev, const md_gpu_tex_desc_t* desc) {
    if (!dev || !desc) return 0;
    md_vk_format_info_t fi = md_vk_format_info(desc->format);
    if (fi.format == VK_FORMAT_UNDEFINED) { md_vk_fail("unsupported texture format %d", (int)desc->format); return 0; }
    if (!(desc->flags & (MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED))) {
        md_vk_fail("texture needs at least one of MD_GPU_TEX_STORAGE / MD_GPU_TEX_SAMPLED");
        return 0;
    }

    bool is_3d = desc->depth > 1;

    VkImageCreateInfo ici = {VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    ici.imageType     = is_3d ? VK_IMAGE_TYPE_3D : VK_IMAGE_TYPE_2D;
    ici.format        = fi.format;
    ici.extent.width  = desc->width  ? desc->width  : 1;
    ici.extent.height = desc->height ? desc->height : 1;
    ici.extent.depth  = is_3d ? desc->depth : 1;
    ici.mipLevels     = 1;
    ici.arrayLayers   = 1;
    ici.samples       = VK_SAMPLE_COUNT_1_BIT;
    ici.tiling        = VK_IMAGE_TILING_OPTIMAL;
    ici.usage         = VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    if (desc->flags & MD_GPU_TEX_STORAGE) ici.usage |= VK_IMAGE_USAGE_STORAGE_BIT;
    if (desc->flags & MD_GPU_TEX_SAMPLED) ici.usage |= VK_IMAGE_USAGE_SAMPLED_BIT;
    ici.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    VkImage image;
    if (!md_vk_check(vkCreateImage(dev->device, &ici, NULL, &image), "vkCreateImage")) return 0;

    VkMemoryRequirements req;
    vkGetImageMemoryRequirements(dev->device, image, &req);
    uint32_t type = md_vk_find_memory_type(dev, req.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0);
    if (type == UINT32_MAX) type = md_vk_find_memory_type(dev, req.memoryTypeBits, 0, 0);

    VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.allocationSize  = req.size;
    mai.memoryTypeIndex = type;
    VkDeviceMemory mem;
    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, &mem), "vkAllocateMemory (image)")) {
        vkDestroyImage(dev->device, image, NULL);
        return 0;
    }
    vkBindImageMemory(dev->device, image, mem, 0);

    VkImageViewCreateInfo vci = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    vci.image    = image;
    vci.viewType = is_3d ? VK_IMAGE_VIEW_TYPE_3D : VK_IMAGE_VIEW_TYPE_2D;
    vci.format   = fi.format;
    vci.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    vci.subresourceRange.levelCount = 1;
    vci.subresourceRange.layerCount = 1;
    VkImageView view;
    if (!md_vk_check(vkCreateImageView(dev->device, &vci, NULL, &view), "vkCreateImageView")) {
        vkFreeMemory(dev->device, mem, NULL);
        vkDestroyImage(dev->device, image, NULL);
        return 0;
    }

    if (!md_vk_immediate(dev, md_vk_record_to_general, &image)) {
        vkDestroyImageView(dev->device, view, NULL);
        vkFreeMemory(dev->device, mem, NULL);
        vkDestroyImage(dev->device, image, NULL);
        return 0;
    }

    md_mutex_lock(&dev->device_mutex);

    /* Storage and sampled images live in separate descriptor arrays, and slots
       come from one shared free list, so a texture usable both as RWTexture
       and as Texture takes two slots. Sharing the free list costs an unused
       entry in the other array but keeps a slot index unique across the whole
       heap, so a handle still identifies exactly one texture. */
    uint32_t storage_slot = UINT32_MAX, sampled_slot = UINT32_MAX;
    if (desc->flags & MD_GPU_TEX_STORAGE) storage_slot = md_vk_alloc_slot(dev);
    if (desc->flags & MD_GPU_TEX_SAMPLED) sampled_slot = md_vk_alloc_slot(dev);
    bool ok = ((desc->flags & MD_GPU_TEX_STORAGE) == 0 || storage_slot != UINT32_MAX)
           && ((desc->flags & MD_GPU_TEX_SAMPLED) == 0 || sampled_slot != UINT32_MAX);
    if (!ok) {
        if (storage_slot != UINT32_MAX) dev->tex_free[dev->tex_free_count++] = storage_slot;
        if (sampled_slot != UINT32_MAX) dev->tex_free[dev->tex_free_count++] = sampled_slot;
        md_mutex_unlock(&dev->device_mutex);
        vkDestroyImageView(dev->device, view, NULL);
        vkFreeMemory(dev->device, mem, NULL);
        vkDestroyImage(dev->device, image, NULL);
        md_vk_fail("out of bindless heap slots (max %u)", MD_VK_MAX_TEXTURES);
        return 0;
    }

    uint32_t primary = (storage_slot != UINT32_MAX) ? storage_slot : sampled_slot;
    md_vk_tex_t* t = &dev->textures[primary];
    t->image = image; t->view = view; t->memory = mem;
    t->desc  = *desc;
    t->storage_slot = storage_slot;
    t->sampled_slot = sampled_slot;
    t->live  = true;
    t->generation++;

    if (storage_slot != UINT32_MAX) md_vk_write_image_slot(dev, storage_slot, view, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    if (sampled_slot != UINT32_MAX) md_vk_write_image_slot(dev, sampled_slot, view, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE);
    md_mutex_unlock(&dev->device_mutex);

    return md_vk_tex_pack(primary);
}

/* The sampled-view handle of a texture created with both usages. Returns the
   texture's own handle when it has only one usage. */
md_gpu_tex_t md_gpu_tex_sampled(md_gpu_tex_t tex) {
    md_vk_tex_t* t = NULL;
    if (!md_vk_tex_device(tex, &t) || !t) return 0;
    return t->sampled_slot != UINT32_MAX ? md_vk_tex_pack(t->sampled_slot) : tex;
}

bool md_gpu_tex_desc(md_gpu_tex_t tex, md_gpu_tex_desc_t* out) {
    md_vk_tex_t* t = NULL;
    if (!out || !md_vk_tex_device(tex, &t) || !t) return false;
    *out = t->desc;
    return true;
}

static md_gpu_device_t md_vk_tex_device(md_gpu_tex_t tex, md_vk_tex_t** out) {
    if (!md_vk_tex_valid(tex)) return NULL;
    uint32_t slot = md_vk_tex_slot(tex);
    if (slot >= MD_VK_MAX_TEXTURES) return NULL;
    for (uint32_t i = 0; i < md_vk_device_count; ++i) {
        md_gpu_device_t d = md_vk_devices[i];
        if (!d) continue;
        md_vk_tex_t* t = &d->textures[slot];
        if (t->live) { if (out) *out = t; return d; }
    }
    return NULL;
}

void md_gpu_tex_destroy(md_gpu_tex_t tex, md_gpu_stream_t stream) {
    md_vk_tex_t* t = NULL;
    md_gpu_device_t dev = md_vk_tex_device(tex, &t);
    if (!dev || !t) return;
    (void)stream;

    md_mutex_lock(&dev->device_mutex);
    md_vk_retire_t* r = (md_vk_retire_t*)md_vk_vec_push(&dev->retires, dev->alloc);
    if (r) {
        r->image     = t->image;
        r->view      = t->view;
        r->memory    = t->memory;
        r->storage_slot = t->storage_slot;
        r->sampled_slot = t->sampled_slot;
        r->wait_count   = 0;
        md_gpu_stream_t* arr = (md_gpu_stream_t*)dev->streams.data;
        for (size_t i = 0; i < dev->streams.count && r->wait_count < 8; ++i) {
            md_gpu_stream_t s = arr[i];
            uint64_t v = s->has_work ? s->next_value : s->submitted_value;
            if (v > 0) {
                r->waits[r->wait_count].stream = s;
                r->waits[r->wait_count].value  = v;
                r->wait_count++;
            }
        }
    }
    /* Ownership has moved to the retire list; clear the slot so that device
       teardown does not destroy these objects a second time. */
    t->live   = false;
    t->image  = VK_NULL_HANDLE;
    t->view   = VK_NULL_HANDLE;
    t->memory = VK_NULL_HANDLE;
    md_mutex_unlock(&dev->device_mutex);
}

md_gpu_sampler_t md_gpu_sampler_create(md_gpu_device_t dev, const md_gpu_sampler_desc_t* desc) {
    if (!dev) return 0;
    md_gpu_sampler_desc_t d = {0};
    if (desc) d = *desc;

    static const VkSamplerAddressMode modes[] = {
        VK_SAMPLER_ADDRESS_MODE_CLAMP_TO_EDGE,
        VK_SAMPLER_ADDRESS_MODE_REPEAT,
        VK_SAMPLER_ADDRESS_MODE_MIRRORED_REPEAT,
    };

    VkSamplerCreateInfo sci = {VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO};
    sci.minFilter    = d.min_filter == MD_GPU_FILTER_LINEAR ? VK_FILTER_LINEAR : VK_FILTER_NEAREST;
    sci.magFilter    = d.mag_filter == MD_GPU_FILTER_LINEAR ? VK_FILTER_LINEAR : VK_FILTER_NEAREST;
    sci.addressModeU = modes[d.address_u % 3];
    sci.addressModeV = modes[d.address_v % 3];
    sci.addressModeW = modes[d.address_w % 3];
    sci.maxLod       = VK_LOD_CLAMP_NONE;
    sci.borderColor  = VK_BORDER_COLOR_FLOAT_TRANSPARENT_BLACK;

    VkSampler sampler;
    if (!md_vk_check(vkCreateSampler(dev->device, &sci, NULL, &sampler), "vkCreateSampler")) return 0;

    md_mutex_lock(&dev->device_mutex);
    if (dev->sampler_free_count == 0) {
        md_mutex_unlock(&dev->device_mutex);
        vkDestroySampler(dev->device, sampler, NULL);
        md_vk_fail("out of sampler slots");
        return 0;
    }
    uint32_t index = dev->sampler_free[--dev->sampler_free_count];
    dev->samplers[index].sampler = sampler;
    dev->samplers[index].live    = true;
    dev->samplers[index].generation++;

    VkDescriptorImageInfo dii = {sampler, VK_NULL_HANDLE, VK_IMAGE_LAYOUT_UNDEFINED};
    VkWriteDescriptorSet w = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
    w.dstSet          = dev->desc_set;
    w.dstBinding      = md_vk_binding_for(VK_DESCRIPTOR_TYPE_SAMPLER);
    w.dstArrayElement = index;
    w.descriptorCount = 1;
    w.descriptorType  = VK_DESCRIPTOR_TYPE_SAMPLER;
    w.pImageInfo      = &dii;
    vkUpdateDescriptorSets(dev->device, 1, &w, 0, NULL);
    md_mutex_unlock(&dev->device_mutex);

    return md_vk_tex_pack(index);
}

void md_gpu_sampler_destroy(md_gpu_sampler_t sampler) {
    if (!md_vk_tex_valid(sampler)) return;
    uint32_t idx = md_vk_tex_slot(sampler);
    if (idx >= MD_VK_MAX_SAMPLERS) return;
    for (uint32_t i = 0; i < md_vk_device_count; ++i) {
        md_gpu_device_t d = md_vk_devices[i];
        if (!d) continue;
        md_vk_sampler_t* s = &d->samplers[idx];
        if (!s->live) continue;
        vkDeviceWaitIdle(d->device);
        md_mutex_lock(&d->device_mutex);
        vkDestroySampler(d->device, s->sampler, NULL);
        s->live = false;
        d->sampler_free[d->sampler_free_count++] = idx;
        md_mutex_unlock(&d->device_mutex);
        return;
    }
}

static void md_vk_resolve_region(const md_gpu_tex_desc_t* d, const md_gpu_tex_region_t* r, VkOffset3D* off, VkExtent3D* ext) {
    uint32_t w = d->width ? d->width : 1;
    uint32_t h = d->height ? d->height : 1;
    uint32_t z = d->depth ? d->depth : 1;
    if (r) {
        off->x = (int32_t)r->offset[0];
        off->y = (int32_t)r->offset[1];
        off->z = (int32_t)r->offset[2];
        ext->width  = r->extent[0] ? r->extent[0] : w - r->offset[0];
        ext->height = r->extent[1] ? r->extent[1] : h - r->offset[1];
        ext->depth  = r->extent[2] ? r->extent[2] : z - r->offset[2];
    } else {
        off->x = off->y = off->z = 0;
        ext->width = w; ext->height = h; ext->depth = z;
    }
}

bool md_gpu_memcpy_to_tex_async(md_gpu_tex_t tex, const md_gpu_tex_region_t* region,
                                const void* src, size_t size, md_gpu_stream_t s) {
    if (!s) return md_vk_fail("md_gpu_memcpy_to_tex_async requires a stream");
    md_vk_tex_t* t = NULL;
    md_gpu_device_t dev = md_vk_tex_device(tex, &t);
    if (!dev || !t) return md_vk_fail("invalid texture handle");

    VkOffset3D off; VkExtent3D ext;
    md_vk_resolve_region(&t->desc, region, &off, &ext);

    md_vk_block_t* sblk = md_vk_registry_find(dev, (uint64_t)(uintptr_t)src);
    VkBuffer src_buffer; uint64_t src_off;
    if (sblk) {
        src_buffer = sblk->buffer;
        src_off    = (uint64_t)(uintptr_t)src - sblk->address;
    } else {
        uint64_t addr; void* host;
        if (!md_vk_arena_alloc(s, size, &addr, &host)) return false;
        memcpy(host, src, size);
        md_vk_page_t* page = NULL;
        md_vk_vec_t* pages = s->capture ? &s->capture->arg_pages : &s->arena.pages;
        for (size_t i = 0; i < pages->count; ++i) {
            md_vk_page_t* p = ((md_vk_page_t**)pages->data)[i];
            if (addr >= p->address && addr < p->address + p->capacity) { page = p; src_off = addr - p->address; break; }
        }
        if (!page) return md_vk_fail("staging address not found");
        src_buffer = page->buffer;
    }

    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);

    VkBufferImageCopy copy = {0};
    copy.bufferOffset = src_off;
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageOffset = off;
    copy.imageExtent = ext;
    vkCmdCopyBufferToImage(cmd, src_buffer, t->image, VK_IMAGE_LAYOUT_GENERAL, 1, &copy);
    md_vk_stream_did_op(s);
    return true;
}

bool md_gpu_memcpy_from_tex_async(void* dst, md_gpu_tex_t tex, const md_gpu_tex_region_t* region,
                                  size_t size, md_gpu_stream_t s) {
    if (!s) return md_vk_fail("md_gpu_memcpy_from_tex_async requires a stream");
    md_vk_tex_t* t = NULL;
    md_gpu_device_t dev = md_vk_tex_device(tex, &t);
    if (!dev || !t) return md_vk_fail("invalid texture handle");

    VkOffset3D off; VkExtent3D ext;
    md_vk_resolve_region(&t->desc, region, &off, &ext);

    md_vk_block_t* dblk = md_vk_registry_find(dev, (uint64_t)(uintptr_t)dst);

    VkBuffer dst_buffer; uint64_t dst_off = 0; void* staging_host = NULL;
    if (dblk) {
        dst_buffer = dblk->buffer;
        dst_off    = (uint64_t)(uintptr_t)dst - dblk->address;
    } else {
        uint64_t addr;
        if (!md_vk_arena_alloc(s, size, &addr, &staging_host)) return false;
        md_vk_page_t* page = NULL;
        md_vk_vec_t* pages = s->capture ? &s->capture->arg_pages : &s->arena.pages;
        for (size_t i = 0; i < pages->count; ++i) {
            md_vk_page_t* p = ((md_vk_page_t**)pages->data)[i];
            if (addr >= p->address && addr < p->address + p->capacity) { page = p; dst_off = addr - p->address; break; }
        }
        if (!page) return md_vk_fail("staging address not found");
        dst_buffer = page->buffer;
    }

    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);

    VkBufferImageCopy copy = {0};
    copy.bufferOffset = dst_off;
    copy.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    copy.imageSubresource.layerCount = 1;
    copy.imageOffset = off;
    copy.imageExtent = ext;
    vkCmdCopyImageToBuffer(cmd, t->image, VK_IMAGE_LAYOUT_GENERAL, dst_buffer, 1, &copy);
    md_vk_stream_did_op(s);

    if (!dblk) {
        md_vk_d2h_t* c = (md_vk_d2h_t*)md_alloc(dev->alloc, sizeof(md_vk_d2h_t));
        if (!c) return md_vk_fail("out of memory");
        c->dev = dev; c->dst = dst; c->staging_host = staging_host; c->size = size;
        md_gpu_sync_t sync = md_gpu_stream_record(s);
        md_mutex_lock(&dev->device_mutex);
        md_vk_hostfn_t* h = (md_vk_hostfn_t*)md_vk_vec_push(&dev->hostfns, dev->alloc);
        if (h) { h->sync = sync; h->fn = md_vk_d2h_finish; h->user = c; h->internal = true; }
        md_mutex_unlock(&dev->device_mutex);
        return h != NULL;
    }
    return true;
}

/* =========================================================================
   10. Kernels and launches
   ========================================================================= */

/* Recover the local size from SPIR-V: OpExecutionMode <entry> LocalSize x y z. */
static bool md_vk_spirv_local_size(const uint32_t* words, size_t word_count, uint32_t out[3]) {
    if (word_count < 5 || words[0] != 0x07230203u) return false;
    size_t i = 5;
    while (i < word_count) {
        uint32_t op    = words[i] & 0xFFFFu;
        uint32_t count = words[i] >> 16;
        if (count == 0 || i + count > word_count) break;
        if (op == 16 /* OpExecutionMode */ && count >= 6) {
            if (words[i + 2] == 17 /* LocalSize */) {
                out[0] = words[i + 3];
                out[1] = words[i + 4];
                out[2] = words[i + 5];
                return true;
            }
        }
        i += count;
    }
    return false;
}

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t dev, const md_gpu_kernel_desc_t* desc) {
    if (!dev || !desc || !desc->code || desc->code_size == 0) {
        md_vk_fail("md_gpu_kernel_create: missing code");
        return NULL;
    }
    if (desc->code_size % 4 != 0) { md_vk_fail("SPIR-V size must be a multiple of 4"); return NULL; }

    VkShaderModuleCreateInfo smci = {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    smci.codeSize = desc->code_size;
    smci.pCode    = (const uint32_t*)desc->code;
    VkShaderModule module;
    if (!md_vk_check(vkCreateShaderModule(dev->device, &smci, NULL, &module), "vkCreateShaderModule")) return NULL;

    VkComputePipelineCreateInfo cpci = {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    cpci.stage.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    cpci.stage.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    cpci.stage.module = module;
    cpci.stage.pName  = desc->entry_point ? desc->entry_point : "main";
    cpci.layout       = dev->pipeline_layout;

    VkPipeline pipeline;
    if (!md_vk_check(vkCreateComputePipelines(dev->device, VK_NULL_HANDLE, 1, &cpci, NULL, &pipeline), "vkCreateComputePipelines")) {
        vkDestroyShaderModule(dev->device, module, NULL);
        return NULL;
    }

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) {
        vkDestroyPipeline(dev->device, pipeline, NULL);
        vkDestroyShaderModule(dev->device, module, NULL);
        md_vk_fail("out of memory");
        return NULL;
    }
    memset(k, 0, sizeof(*k));
    k->device   = dev;
    k->module   = module;
    k->pipeline = pipeline;
    snprintf(k->label, sizeof(k->label), "%s", desc->label ? desc->label : "kernel");

    if (!md_vk_spirv_local_size((const uint32_t*)desc->code, desc->code_size / 4, k->group_size)) {
        k->group_size[0] = desc->group_size[0] ? desc->group_size[0] : 1;
        k->group_size[1] = desc->group_size[1] ? desc->group_size[1] : 1;
        k->group_size[2] = desc->group_size[2] ? desc->group_size[2] : 1;
    }
    return k;
}

void md_gpu_kernel_destroy(md_gpu_kernel_t k) {
    if (!k) return;
    md_gpu_device_t dev = k->device;
    vkDeviceWaitIdle(dev->device);
    vkDestroyPipeline(dev->device, k->pipeline, NULL);
    vkDestroyShaderModule(dev->device, k->module, NULL);
    md_free(dev->alloc, k, sizeof(*k));
}

bool md_gpu_kernel_info(md_gpu_kernel_t k, md_gpu_kernel_info_t* info) {
    if (!k || !info) return false;
    info->max_threads_per_group    = k->device->props.limits.maxComputeWorkGroupInvocations;
    info->preferred_group_multiple = 32;
    info->group_size[0] = k->group_size[0];
    info->group_size[1] = k->group_size[1];
    info->group_size[2] = k->group_size[2];
    return true;
}

/* Place the argument block and return its device address. During capture the
   block is owned by the graph so that it survives and can be patched. */
static bool md_vk_place_args(md_gpu_stream_t s, const void* args, size_t size, uint64_t* out_addr) {
    if (size == 0) { *out_addr = 0; return true; }

    if (s->capture) {
        md_gpu_graph_t g = s->capture;
        md_gpu_device_t dev = s->device;
        uint64_t need = md_vk_align_up(size, MD_VK_ARG_ALIGN);

        md_vk_page_t* page = NULL;
        size_t page_index = 0;
        for (size_t i = 0; i < g->arg_pages.count; ++i) {
            md_vk_page_t* p = ((md_vk_page_t**)g->arg_pages.data)[i];
            if (p->cursor + need <= p->capacity) { page = p; page_index = i; break; }
        }
        if (!page) {
            uint64_t page_size = MD_VK_GRAPH_ARG_PAGE;
            while (page_size < need) page_size *= 2;
            page = md_vk_page_create(dev, page_size);
            if (!page) return md_vk_fail("failed to allocate graph argument page");
            md_vk_page_t** slot = (md_vk_page_t**)md_vk_vec_push(&g->arg_pages, dev->alloc);
            if (!slot) { md_vk_page_destroy(dev, page); return false; }
            *slot = page;
            page_index = g->arg_pages.count - 1;
        }

        uint64_t offset = page->cursor;
        memcpy(page->host + offset, args, size);
        page->cursor += need;

        md_vk_graph_launch_t* l = (md_vk_graph_launch_t*)md_vk_vec_push(&g->launches, dev->alloc);
        if (!l) return false;
        l->page_index = page_index;
        l->offset     = offset;
        l->size       = size;

        *out_addr = page->address + offset;
        return true;
    }

    void* host;
    if (!md_vk_arena_alloc(s, size, out_addr, &host)) return false;
    memcpy(host, args, size);
    return true;
}

static bool md_vk_launch_common(md_gpu_stream_t s, md_gpu_kernel_t k, const void* args, size_t args_size,
                                VkCommandBuffer* out_cmd) {
    md_gpu_device_t dev = s->device;
    if (k->device != dev) return md_vk_fail("kernel belongs to a different device");

    uint64_t arg_addr = 0;
    if (!md_vk_place_args(s, args, args_size, &arg_addr)) return false;

    VkCommandBuffer cmd = md_vk_target_cmd(s);
    if (!cmd) return false;
    md_vk_stream_barrier_if_needed(s);

    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, k->pipeline);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, dev->pipeline_layout, 0, 1, &dev->desc_set, 0, NULL);
    vkCmdPushConstants(cmd, dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, 8, &arg_addr);

    *out_cmd = cmd;
    return true;
}

bool md_gpu_launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid,
                    const void* args, size_t args_size) {
    if (!s || !k) return md_vk_fail("md_gpu_launch: null stream or kernel");
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;   /* empty launch */

    VkCommandBuffer cmd;
    if (!md_vk_launch_common(s, k, args, args_size, &cmd)) return false;
    vkCmdDispatch(cmd, grid.x, grid.y, grid.z);
    md_vk_stream_did_op(s);
    return true;
}

bool md_gpu_launch_indirect(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_ptr_t grid_ptr,
                             const void* args, size_t args_size) {
    if (!s || !k || !grid_ptr) return md_vk_fail("md_gpu_launch_indirect: null argument");
    md_vk_block_t* b = md_vk_registry_find(s->device, (uint64_t)(uintptr_t)grid_ptr);
    if (!b) return md_vk_fail("md_gpu_launch_indirect: grid_ptr is not a device pointer");

    VkCommandBuffer cmd;
    if (!md_vk_launch_common(s, k, args, args_size, &cmd)) return false;
    vkCmdDispatchIndirect(cmd, b->buffer, (uint64_t)(uintptr_t)grid_ptr - b->address);
    md_vk_stream_did_op(s);
    return true;
}

typedef struct md_vk_make_grid_args_t {
    uint64_t count;
    uint64_t out_grid;
    uint32_t local[3];
    uint32_t _pad;
} md_vk_make_grid_args_t;

static bool md_vk_create_builtin_kernels(md_gpu_device_t dev) {
    md_gpu_kernel_desc_t d = {0};
    d.code      = md_gpu_make_grid_spv;
    d.code_size = sizeof(md_gpu_make_grid_spv);
    d.label     = "md_gpu make_grid";
    dev->make_grid_kernel = md_gpu_kernel_create(dev, &d);
    return dev->make_grid_kernel != NULL;
}

bool md_gpu_make_grid(md_gpu_stream_t s, md_gpu_ptr_t out_grid, const md_gpu_ptr_t count, const uint32_t local[3]) {
    if (!s || !out_grid || !count) return md_vk_fail("md_gpu_make_grid: null argument");
    md_vk_make_grid_args_t a = {0};
    a.count    = (uint64_t)(uintptr_t)count;
    a.out_grid = (uint64_t)(uintptr_t)out_grid;
    a.local[0] = local && local[0] ? local[0] : 1;
    a.local[1] = local && local[1] ? local[1] : 1;
    a.local[2] = local && local[2] ? local[2] : 1;
    return md_gpu_launch(s, s->device->make_grid_kernel, md_gpu_grid(1, 1, 1), &a, sizeof(a));
}

/* =========================================================================
   11. Graphs
   ========================================================================= */

bool md_gpu_capture_begin(md_gpu_stream_t s, const char* label) {
    if (!s) return md_vk_fail("md_gpu_capture_begin: null stream");
    if (s->capture) return md_vk_fail("stream '%s' is already capturing", s->label);
    md_gpu_device_t dev = s->device;

    /* Close anything already issued so it is not swallowed by the capture. */
    md_vk_stream_submit(s);

    md_gpu_graph_t g = (md_gpu_graph_t)md_alloc(dev->alloc, sizeof(md_gpu_graph));
    if (!g) return md_vk_fail("out of memory");
    memset(g, 0, sizeof(*g));
    g->device = dev;
    g->kind   = s->kind;
    snprintf(g->label, sizeof(g->label), "%s", label ? label : "graph");
    md_vk_vec_init(&g->arg_pages, sizeof(md_vk_page_t*));
    md_vk_vec_init(&g->launches,  sizeof(md_vk_graph_launch_t));

    VkCommandPoolCreateInfo cpci = {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    cpci.queueFamilyIndex = s->family;
    if (!md_vk_check(vkCreateCommandPool(dev->device, &cpci, NULL, &g->pool), "vkCreateCommandPool (graph)")) {
        md_free(dev->alloc, g, sizeof(*g));
        return false;
    }
    VkCommandBufferAllocateInfo cbai = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cbai.commandPool        = g->pool;
    cbai.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cbai.commandBufferCount = 1;
    if (!md_vk_check(vkAllocateCommandBuffers(dev->device, &cbai, &g->cmd), "vkAllocateCommandBuffers (graph)")) {
        vkDestroyCommandPool(dev->device, g->pool, NULL);
        md_free(dev->alloc, g, sizeof(*g));
        return false;
    }

    /* Recorded once, launched many times, possibly while a previous launch is
       still in flight. */
    VkCommandBufferBeginInfo bi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_SIMULTANEOUS_USE_BIT;
    if (!md_vk_check(vkBeginCommandBuffer(g->cmd, &bi), "vkBeginCommandBuffer (graph)")) {
        vkDestroyCommandPool(dev->device, g->pool, NULL);
        md_free(dev->alloc, g, sizeof(*g));
        return false;
    }

    s->capture = g;
    /* A graph is replayed back to back, so its first operation must depend on
       whatever ran before this launch -- including the previous replay of the
       same graph. Recording the barrier into the graph makes every launch safe
       regardless of what preceded it. */
    s->needs_barrier = true;
    return true;
}

bool md_gpu_is_capturing(md_gpu_stream_t s) { return s && s->capture != NULL; }

uint32_t md_gpu_capture_next_index(md_gpu_stream_t s) {
    if (!s || !s->capture) return 0;
    return (uint32_t)s->capture->launches.count;
}

md_gpu_graph_t md_gpu_capture_end(md_gpu_stream_t s) {
    if (!s || !s->capture) { md_vk_fail("stream is not capturing"); return NULL; }
    md_gpu_graph_t g = s->capture;
    s->capture = NULL;
    s->needs_barrier = false;
    if (!md_vk_check(vkEndCommandBuffer(g->cmd), "vkEndCommandBuffer (graph)")) {
        md_gpu_graph_destroy(g);
        return NULL;
    }
    return g;
}

uint32_t md_gpu_graph_launch_count(md_gpu_graph_t g) {
    return g ? (uint32_t)g->launches.count : 0;
}

void* md_gpu_graph_args(md_gpu_graph_t g, uint32_t index) {
    if (!g || index >= g->launches.count) return NULL;
    md_vk_graph_launch_t* l = &((md_vk_graph_launch_t*)g->launches.data)[index];
    md_vk_page_t* p = ((md_vk_page_t**)g->arg_pages.data)[l->page_index];
    return p->host + l->offset;
}

bool md_gpu_graph_launch(md_gpu_graph_t g, md_gpu_stream_t s) {
    if (!g || !s) return md_vk_fail("md_gpu_graph_launch: null argument");
    if (s->capture) return md_vk_fail("cannot launch a graph while capturing");
    if (g->kind != s->kind) return md_vk_fail("graph was captured on a different stream kind");
    md_gpu_device_t dev = s->device;

    /* Flush anything pending so the graph is ordered after it. */
    md_vk_stream_submit(s);

    uint64_t signal_value = s->next_value;

    VkCommandBufferSubmitInfo cbsi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO};
    cbsi.commandBuffer = g->cmd;

    VkSemaphoreSubmitInfo waits[MD_VK_MAX_PENDING_WAITS];
    for (uint32_t i = 0; i < s->wait_count; ++i) {
        waits[i] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        waits[i].semaphore = s->wait_sems[i];
        waits[i].value     = s->wait_vals[i];
        waits[i].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }

    VkSemaphoreSubmitInfo signal = {VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
    signal.semaphore = s->timeline;
    signal.value     = signal_value;
    signal.stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;

    VkSubmitInfo2 si = {VK_STRUCTURE_TYPE_SUBMIT_INFO_2};
    si.waitSemaphoreInfoCount   = s->wait_count;
    si.pWaitSemaphoreInfos      = waits;
    si.commandBufferInfoCount   = 1;
    si.pCommandBufferInfos      = &cbsi;
    si.signalSemaphoreInfoCount = 1;
    si.pSignalSemaphoreInfos    = &signal;

    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueueSubmit2(s->queue, 1, &si, VK_NULL_HANDLE);
    md_mutex_unlock(&dev->queue_mutex);
    if (!md_vk_check(r, "vkQueueSubmit2 (graph)")) return false;

    s->submitted_value = signal_value;
    s->next_value      = signal_value + 1;
    s->wait_count      = 0;
    return true;
}

void md_gpu_graph_destroy(md_gpu_graph_t g) {
    if (!g) return;
    md_gpu_device_t dev = g->device;
    vkDeviceWaitIdle(dev->device);
    for (size_t i = 0; i < g->arg_pages.count; ++i) {
        md_vk_page_destroy(dev, ((md_vk_page_t**)g->arg_pages.data)[i]);
    }
    md_vk_vec_free(&g->arg_pages, dev->alloc);
    md_vk_vec_free(&g->launches,  dev->alloc);
    if (g->pool) vkDestroyCommandPool(dev->device, g->pool, NULL);
    md_free(dev->alloc, g, sizeof(*g));
}

/* =========================================================================
   12. Host callbacks and polling
   ========================================================================= */

bool md_gpu_launch_host_fn(md_gpu_stream_t s, md_gpu_host_fn fn, void* user) {
    if (!s || !fn) return md_vk_fail("md_gpu_launch_host_fn: null argument");
    if (s->capture) return md_vk_fail("host functions cannot be captured into a graph");
    md_gpu_sync_t sync = md_gpu_stream_record(s);
    return md_gpu_sync_on_complete(sync, fn, user);
}

bool md_gpu_sync_on_complete(md_gpu_sync_t sync, md_gpu_host_fn fn, void* user) {
    if (!fn) return false;
    md_gpu_device_t dev = NULL;
    if (md_gpu_sync_is_valid(sync)) dev = sync.stream->device;
    else if (md_vk_device_count > 0)  dev = md_vk_devices[0];
    if (!dev) return md_vk_fail("md_gpu_sync_on_complete: no device");

    md_mutex_lock(&dev->device_mutex);
    md_vk_hostfn_t* h = (md_vk_hostfn_t*)md_vk_vec_push(&dev->hostfns, dev->alloc);
    if (h) { h->sync = sync; h->fn = fn; h->user = user; h->internal = false; }
    md_mutex_unlock(&dev->device_mutex);
    return h != NULL;
}

/* Timeline values sampled once per poll pass.

   Callbacks must fire in the order they were registered. Deciding readiness by
   re-reading the semaphore for every entry breaks that: two callbacks
   registered against the *same* sync value can resolve differently, because
   the timeline may advance between the two reads, and then the later one
   overtakes the earlier. A device-to-host copy followed by
   md_gpu_launch_host_fn is exactly that pattern, and the callback would
   observe memory the copy had not written yet. Snapshotting makes readiness a
   property of the pass rather than of the instant. */
typedef struct md_vk_timeline_snapshot_t {
    md_gpu_stream_t stream[16];
    uint64_t         completed[16];
    uint32_t         count;
} md_vk_timeline_snapshot_t;

static bool md_vk_snapshot_complete(md_vk_timeline_snapshot_t* snap, md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    for (uint32_t i = 0; i < snap->count; ++i) {
        if (snap->stream[i] == sync.stream) return snap->completed[i] >= sync.value;
    }
    uint64_t completed = md_vk_stream_completed(sync.stream);
    if (snap->count < (uint32_t)(sizeof(snap->stream) / sizeof(snap->stream[0]))) {
        snap->stream[snap->count]    = sync.stream;
        snap->completed[snap->count] = completed;
        snap->count++;
    }
    return completed >= sync.value;
}

/* Fire completed callbacks and reclaim retired objects.
   `fire_user` selects whether user callbacks run; internal ones always do. */
static uint32_t md_vk_poll_internal(md_gpu_device_t dev, bool fire_user) {
    uint32_t fired = 0;
    md_vk_timeline_snapshot_t snap = {0};

    for (;;) {
        md_vk_hostfn_t ready;
        bool have = false;

        md_mutex_lock(&dev->device_mutex);
        md_vk_hostfn_t* arr = (md_vk_hostfn_t*)dev->hostfns.data;
        for (size_t i = 0; i < dev->hostfns.count; ++i) {
            if (!fire_user && !arr[i].internal) continue;
            if (!md_vk_snapshot_complete(&snap, arr[i].sync)) continue;
            ready = arr[i];
            memmove(&arr[i], &arr[i + 1], (dev->hostfns.count - i - 1) * sizeof(*arr));
            dev->hostfns.count--;
            have = true;
            break;
        }
        md_mutex_unlock(&dev->device_mutex);

        if (!have) break;
        ready.fn(ready.user);       /* outside the lock */
        fired++;
    }

    /* Retire textures whose last users have completed. */
    md_mutex_lock(&dev->device_mutex);
    md_vk_retire_t* rr = (md_vk_retire_t*)dev->retires.data;
    for (size_t i = 0; i < dev->retires.count;) {
        md_vk_retire_t* r = &rr[i];
        bool done = true;
        for (uint32_t w = 0; w < r->wait_count && done; ++w) {
            if (!md_gpu_sync_is_complete(r->waits[w])) done = false;
        }
        if (!done) { ++i; continue; }
        /* Drop the descriptors first: nothing may reference a destroyed view. */
        if (r->storage_slot != UINT32_MAX) md_vk_clear_slot(dev, r->storage_slot, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
        if (r->sampled_slot != UINT32_MAX) md_vk_clear_slot(dev, r->sampled_slot, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE);
        if (r->view)   vkDestroyImageView(dev->device, r->view, NULL);
        if (r->image)  vkDestroyImage(dev->device, r->image, NULL);
        if (r->buffer) vkDestroyBuffer(dev->device, r->buffer, NULL);
        if (r->memory) vkFreeMemory(dev->device, r->memory, NULL);
        if (r->storage_slot != UINT32_MAX && dev->tex_free_count < MD_VK_MAX_TEXTURES) dev->tex_free[dev->tex_free_count++] = r->storage_slot;
        if (r->sampled_slot != UINT32_MAX && dev->tex_free_count < MD_VK_MAX_TEXTURES) dev->tex_free[dev->tex_free_count++] = r->sampled_slot;
        memmove(&rr[i], &rr[i + 1], (dev->retires.count - i - 1) * sizeof(*rr));
        dev->retires.count--;
    }
    md_mutex_unlock(&dev->device_mutex);

    return fired;
}

uint32_t md_gpu_device_poll(md_gpu_device_t dev) {
    if (!dev) return 0;
    return md_vk_poll_internal(dev, true);
}

/* =========================================================================
   Device destruction
   ========================================================================= */

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;

    if (dev->device) {
        /* Flush and idle every stream. */
        for (size_t i = 0; i < dev->streams.count; ++i) {
            md_gpu_stream_t s = ((md_gpu_stream_t*)dev->streams.data)[i];
            if (s->capture) { s->capture = NULL; }
            md_vk_stream_submit(s);
        }
        vkDeviceWaitIdle(dev->device);
        md_vk_poll_internal(dev, true);

        if (dev->make_grid_kernel) md_gpu_kernel_destroy(dev->make_grid_kernel);

        for (size_t i = 0; i < dev->streams.count; ++i) {
            md_vk_stream_destroy_internal(((md_gpu_stream_t*)dev->streams.data)[i]);
        }
        dev->streams.count = 0;
        md_vk_vec_free(&dev->streams, alloc);

        while (dev->pools.count > 0) {
            md_gpu_pool_destroy(((md_gpu_pool_t*)dev->pools.data)[0]);
        }
        md_vk_vec_free(&dev->pools, alloc);

        for (uint32_t i = 0; i < MD_VK_MAX_TEXTURES; ++i) {
            md_vk_tex_t* t = &dev->textures[i];
            if (!t->image) continue;
            if (t->view)   vkDestroyImageView(dev->device, t->view, NULL);
            vkDestroyImage(dev->device, t->image, NULL);
            if (t->memory) vkFreeMemory(dev->device, t->memory, NULL);
            memset(t, 0, sizeof(*t));
        }
        for (uint32_t i = 0; i < MD_VK_MAX_SAMPLERS; ++i) {
            if (dev->samplers[i].live) vkDestroySampler(dev->device, dev->samplers[i].sampler, NULL);
        }
        if (dev->dummy_view)  vkDestroyImageView(dev->device, dev->dummy_view, NULL);
        if (dev->dummy_image) vkDestroyImage(dev->device, dev->dummy_image, NULL);
        if (dev->dummy_mem)   vkFreeMemory(dev->device, dev->dummy_mem, NULL);
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

    for (uint32_t i = 0; i < md_vk_device_count; ++i) {
        if (md_vk_devices[i] == dev) {
            md_vk_devices[i] = md_vk_devices[md_vk_device_count - 1];
            md_vk_device_count--;
            break;
        }
    }

    md_free(alloc, dev, sizeof(*dev));
}

md_gpu_mem_flags_t md_gpu_pool_flags(md_gpu_pool_t p) { return p ? p->flags : 0; }

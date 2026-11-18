/*
md_gpu_vulkan.c -- Vulkan backend for md_gpu.h

Structure of this file:

    1.  Configuration, error handling, small utilities
    2.  Types
    3.  Address lookup (device address -> buffer)
    4.  Raw buffers and transient arenas (argument blocks and staging)
    5.  Deferred destruction
    6.  Device creation
    7.  Bindless descriptor set
    8.  Streams, command buffers, submission, ordering
    9.  Memory: chunks for md_gpu_heap, malloc/free, temp scopes, copies, uploads
    10. Textures and samplers
    11. Kernels and launches
    12. Rendering: pipelines, passes, dynamic state, draws
    13. Presentation: surfaces and swapchains
    14. Host callbacks and polling
    15. Device destruction

The dependency model is program order within a stream, implemented as a single
global VkMemoryBarrier2 between consecutive operations (IMPLICIT ordering) --
from the stages the operations before it wrote in to the stages the next one
reads in, see md_vk_begin_op -- or as caller-placed stage barriers (EXPLICIT
ordering).
There is no per-resource state tracking anywhere in this file, and every image
lives in VK_IMAGE_LAYOUT_GENERAL for its entire life -- render targets
included, since dynamic rendering takes GENERAL attachments. Swapchain images
are the one exception, and only at the edges: they enter GENERAL at acquire
and leave it for PRESENT_SRC at present.

Nothing here blocks the calling thread except md_gpu_stream_sync,
md_gpu_sync_wait, md_gpu_stream_destroy (on its own stream) and device
creation/destruction. Destroying textures and kernels is deferred: the
object goes onto a retire list stamped with every stream's current position
and is released by md_gpu_device_poll once all of those have passed. That is
what lets a compute job run for many frames without anything else stalling
behind it.

Resources are created VK_SHARING_MODE_CONCURRENT across the compute,
transfer and graphics families whenever those differ, so a buffer or image
written on one and read on another needs no queue-family ownership transfer.
*/

#include "md_gpu.h"
#include "md_gpu_heap.h"
#include "md_gpu_error.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_format.h>
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
#include "md_gpu_select.inl"

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
/* Submissions a stream keeps the GPU times of (md_gpu_sync_gpu_seconds): two
   timestamps each, in a ring of query pairs */
#define MD_VK_TS_RING                256u
#define MD_VK_DEFAULT_TEXTURE_SLOTS 4096u /* sampled + storage views share them */
#define MD_VK_MAX_SAMPLERS         256u
#define MD_VK_ARENA_PAGE_SIZE   (256u * 1024u)
/* A transient allocation larger than this (an upload's staging, mostly) gets a page of its own, of its
   exact size, in host memory rather than in the host visible device local heap, and the page is freed
   once the submission that reads it completes. Doubling the shared pages up to it instead would take a
   512 MB page for a 300 MB upload, from a heap that is 256 MB on a discrete GPU without resizable BAR,
   and keep it for the stream's lifetime. */
#define MD_VK_ARENA_DEDICATED   (16ull * 1024u * 1024u)
#define MD_VK_ARG_ALIGN            64u
#define MD_VK_DEBUG_GUARD        (64ull << 10)     /* each side of a debug memory allocation */



/* The instance whose entry points volk holds while a device using it is alive
   (VK_NULL_HANDLE otherwise), so md_gpu_enumerate_adapters can put them back. */
static VkInstance md_vk_live_instance = VK_NULL_HANDLE;

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
    return md_gpu_fail("%s failed: %s (%d)", what, md_vk_result_str(r), (int)r);
}

/* =========================================================================
   2. Types
   ========================================================================= */

/* One VkBuffer with its memory: a chunk of a heap or of a temp arena
   (md_gpu_heap.h), or one debug memory allocation. */
typedef struct md_vk_chunk_t {
    md_gpu_chunk_t    base;          /* first: md_gpu_heap hands out &base      */
    VkBuffer          buffer;
    VkDeviceMemory    memory;
    uint32_t          memory_type;   /* of `memory`                             */
    uint64_t          memory_offset; /* of the buffer in `memory`               */
} md_vk_chunk_t;

/* A resolved address: where to point a command at. */
typedef struct md_vk_span_t {
    VkBuffer       buffer;
    uint64_t       offset;           /* within `buffer`                         */
    uint8_t*       host;             /* host pointer of the address, or NULL    */
    md_vk_chunk_t* chunk;            /* the chunk `buffer` is                   */
} md_vk_span_t;

/* One page of host-visible, device-addressable transient memory. */
typedef struct md_vk_page_t {
    VkBuffer       buffer;
    VkDeviceMemory memory;
    uint64_t       address;
    uint8_t*       host;
    uint64_t       capacity;
    uint64_t       cursor;
    uint64_t       retire_value;   /* stream timeline value that frees it */
    bool           dedicated;      /* one large allocation's own page: destroyed once retired, never refilled */
} md_vk_page_t;

typedef struct md_vk_arena_t {
    MD_ARRAY_RAW_TYPE(md_vk_page_t*) pages;
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
    bool                 can_graphics;     /* GRAPHICS streams                */
    VkQueue              queue;
    VkSemaphore          timeline;
    uint64_t             next_value;       /* value the next submit signals */
    uint64_t             submitted_value;  /* last value submitted          */

    VkCommandPool        cmd_pool;
    MD_ARRAY_RAW_TYPE(md_vk_cmd_t) cmds;
    VkCommandBuffer      open;             /* currently recording, or NULL  */
    bool                 has_work;
    /* Stages of the operations recorded since md_gpu last placed a barrier,
       across submissions: what the next barrier waits for. */
    md_gpu_stage_flags_t written;
    bool                 order_next;       /* order the next op after `written`
                                              in either ordering mode (the first
                                              op of an EXPLICIT region)        */
    bool                 cb_head;          /* nothing recorded yet in `open`   */
    bool                 force_barrier;    /* full barrier before the next op, in
                                              either ordering mode (memory reuse) */
    md_gpu_ordering_t    ordering;
    /* Bindings in the open command buffer, set once and kept. */
    struct md_gpu_kernel* bound_kernel;
    bool                 compute_set_bound;
    bool                 graphics_set_bound;

    MD_ARRAY_RAW_TYPE(md_vk_wait_t) waits; /* pending for the next submit */

    md_vk_arena_t        arena;

    md_gpu_heap_temps_t  temps;            /* md_gpu_temp_* scopes */

    /* Binary semaphores for the next submit: swapchain acquire and present. */
    VkSemaphore          bin_waits[4];
    uint32_t             bin_wait_count;
    VkSemaphore          bin_signals[4];
    uint32_t             bin_signal_count;

    /* The open render pass, if any. */
    bool                 in_pass;
    bool                 pass_labelled;
    uint32_t             pass_color_count;
    md_gpu_format_t      pass_color[MD_GPU_MAX_COLOR_TARGETS];
    md_gpu_format_t      pass_depth;
    uint32_t             pass_width, pass_height;
    struct md_gpu_pipeline* bound_pipeline;
    md_gpu_draw_state_t  draw_state;       /* as currently set in the command buffer */

    /* Open upload, if any. */
    bool                 upload_open;
    bool                 upload_direct;
    md_gpu_addr_t        upload_dst;
    VkBuffer             upload_src;           /* staging: an arena page        */
    uint64_t             upload_src_offset;
    size_t               upload_size;

    /* The GPU time of the last MD_VK_TS_RING submissions: a pair of timestamps
       per submission, written at the top of its command buffer and at its end,
       in slot (value % MD_VK_TS_RING); ts_value names the submission a slot
       holds. VK_NULL_HANDLE where the queue cannot time. */
    VkQueryPool          ts_pool;
    uint64_t             ts_mask;          /* the valid bits of a timestamp */
    uint64_t             ts_value[MD_VK_TS_RING];

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
    VkImage               image;
    VkDeviceMemory        memory;
    uint64_t              bytes;           /* memory size, for stats         */
    VkImageView           sampled_view;    /* all mips, or VK_NULL_HANDLE    */
    uint32_t              sampled_slot;    /* heap slot, or 0                */
    VkImageView*          storage_views;   /* one per mip, or NULL           */
    uint32_t*             storage_slots;   /* one per mip, or NULL           */
    md_vk_format_info_t   fi;
    md_gpu_texture_desc_t desc;            /* normalised                     */
    struct md_vk_attach_view_t* attach_views;   /* render-target views, made
                                                   on first use, one per
                                                   (mip, layer)              */
    uint32_t              attach_view_count, attach_view_cap;
    bool                  external;        /* a swapchain image: md_gpu owns
                                              neither the image nor memory */
    char                  label[64];
} md_gpu_texture;

typedef struct md_vk_attach_view_t {
    uint32_t    mip, layer;
    VkImageView view;
} md_vk_attach_view_t;

typedef struct md_gpu_pipeline {
    md_gpu_device_t   device;
    VkPipeline        pipeline;
    uint32_t          color_count;
    md_gpu_format_t   color[MD_GPU_MAX_COLOR_TARGETS];
    md_gpu_format_t   depth;
    uint32_t          args_size;
    char              label[64];
} md_gpu_pipeline;

/* A swapchain and what hangs off it. Retired as one object: on a rebuild the
   old one goes, on surface destruction the last one takes the VkSurfaceKHR
   and the acquire semaphores with it. */
typedef struct md_vk_swapchain_t {
    VkSwapchainKHR    swapchain;
    VkQueue           present_queue;   /* idled before destruction, or NULL */
    uint32_t          image_count;
    md_gpu_texture_t* textures;        /* external wrappers, one per image  */
    VkSemaphore*      present_sems;    /* one per image                     */
    uint32_t          width, height;
    /* Set only on the final retirement of a surface. */
    VkSurfaceKHR      surface;
    VkSemaphore*      extra_sems;
    uint32_t          extra_sem_count;
} md_vk_swapchain_t;

#define MD_VK_ACQUIRE_SEMS 8u

typedef struct md_gpu_surface {
    md_gpu_device_t        device;
    VkSurfaceKHR           surface;
    md_gpu_surface_desc_t  desc;
    char                   label[64];
    VkSurfaceFormatKHR     vk_format;
    VkPresentModeKHR       present_mode;
    VkImageUsageFlags      image_usage;

    uint32_t               want_width, want_height;
    bool                   dirty;          /* rebuild at the next acquire     */
    md_vk_swapchain_t*     sc;

    /* Acquire semaphores, used round robin. Slot i may be reused once the
       submission that waited on it -- at or before (stream, value) -- has
       completed. */
    VkSemaphore            acquire_sems[MD_VK_ACQUIRE_SEMS];
    md_gpu_stream_t        acquire_stream[MD_VK_ACQUIRE_SEMS];
    uint64_t               acquire_value[MD_VK_ACQUIRE_SEMS];
    uint32_t               acquire_next;

    bool                   acquired;
    uint32_t               image_index;
    uint32_t               acquire_slot;
    md_gpu_stream_t        acquired_stream;
    VkQueue                present_queue;  /* last queue presented on         */
} md_gpu_surface;

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
    char            stats[256];     /* the driver's statistics of the compiled pipeline, if reported */
} md_gpu_kernel;

typedef struct md_vk_hostfn_t {
    md_gpu_sync_t  sync;
    md_gpu_host_fn fn;
    void*          user;
} md_vk_hostfn_t;

typedef enum md_vk_retire_kind_t {
    MD_VK_RETIRE_TEXTURE,
    MD_VK_RETIRE_KERNEL,
    MD_VK_RETIRE_PIPELINE,
    MD_VK_RETIRE_SWAPCHAIN,
} md_vk_retire_kind_t;

/* An object whose destruction waits for every stream to pass the point at
   which it was destroyed. Owns `waits`. */
typedef struct md_vk_retire_t {
    md_vk_retire_kind_t kind;
    void*               object;     /* md_gpu_texture_t / md_gpu_kernel_t */
    md_vk_wait_t*       waits;
    uint32_t            wait_count;
    uint32_t            wait_capacity;   /* what md_free must be told */
} md_vk_retire_t;

/* Everything the backend asks of a device, in one place, so that a device we
   cannot use is rejected by name instead of failing later inside
   vkCreateDevice with a bare VkResult. The optional flags are enabled only
   when the driver reports them. */
typedef struct md_vk_dev_caps_t {
    bool graphics;        /* everything rendering needs, below */
    bool depth_bias_clamp;
    bool maintenance4;
    bool update_unused_while_pending;
    bool nonuniform_storage_image;
    bool nonuniform_sampled_image;
    bool dynamic_storage_image;
    bool dynamic_sampled_image;
    bool shader_int64;
    bool storage_read_without_format;    /* device-wide; otherwise per format */
    bool storage_write_without_format;
    /* Above the floor, enabled where present and reported in md_gpu_device_info_t. */
    bool float16;                        /* shaderFloat16 + storageBuffer16BitAccess */
    bool float64;
    bool int64_atomics;                  /* buffer; shared ones too where present */
    bool int64_shared_atomics;
    bool float_atomics;                  /* VK_EXT_shader_atomic_float           */
    VkPhysicalDeviceShaderAtomicFloatFeaturesEXT  atomic_float;
    bool float_atomics2;                 /* VK_EXT_shader_atomic_float2          */
    VkPhysicalDeviceShaderAtomicFloat2FeaturesEXT atomic_float2;
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
    VkSubgroupFeatureFlags           subgroup_ops;     /* in compute kernels */
    char                             driver_desc[256];

    uint32_t compute_family, transfer_family, graphics_family;   /* graphics: UINT32_MAX if none */
    bool     transfer_can_compute;
    /* Timestamps: the valid bits of each family's (0: the family cannot time)
       and the nanoseconds of a tick */
    uint32_t ts_bits_compute, ts_bits_transfer, ts_bits_graphics;
    float    ts_period;
    /* Families a resource must be shared across (CONCURRENT when > 1). */
    uint32_t share_families[3];
    uint32_t share_family_count;

    /* Streams are spread round-robin over these. A device that exposes several
       compute queues (typical: 8 on NVIDIA) then gets real concurrency from
       "use another stream" rather than just permission to overlap. */
    VkQueue  compute_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t compute_queue_count;
    VkQueue  transfer_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t transfer_queue_count;
    VkQueue  graphics_queues[MD_VK_MAX_QUEUES_PER_FAMILY];
    uint32_t graphics_queue_count;
    uint32_t next_compute_queue, next_transfer_queue, next_graphics_queue;
    md_mutex_t queue_mutex;      /* streams may share a VkQueue */
    md_mutex_t device_mutex;     /* heaps, registry, textures, retire and callback lists */

    /* Bindless */
    VkDescriptorSetLayout set_layout;
    VkDescriptorPool      desc_pool;
    VkDescriptorSet       desc_set;
    VkPipelineLayout      pipeline_layout;         /* compute: push constant for COMPUTE           */
    VkPipelineLayout      raster_layout;           /* draws: push constant for VERTEX | FRAGMENT   */

    /* A 1x1x1 placeholder. Freed descriptor slots are pointed at it so that no
       descriptor ever references a destroyed view. */
    VkImage         dummy_image;
    VkImageView     dummy_view;
    VkDeviceMemory  dummy_mem;
    VkSampler       dummy_sampler;

    /* Texture heap slots; slot 0 is never handed out so zero stays null. */
    uint32_t*       tex_free;              /* tex_slot_count - 1 entries        */
    uint32_t        tex_free_count;
    uint32_t        tex_slot_count;        /* slots in each image binding       */
    uint32_t        device_buffer_types;   /* memory types MD_GPU_MEM_DEVICE chunks may use */

    md_vk_sampler_entry_t samplers[MD_VK_MAX_SAMPLERS];
    uint32_t              sampler_count;   /* slots 1..sampler_count used */

    md_gpu_heap_t   heap;       /* md_gpu_malloc, md_gpu_free, temp scopes, address lookup */
    MD_ARRAY_RAW_TYPE(md_gpu_texture_t) textures; /* live */
    uint64_t        texture_bytes;
    MD_ARRAY_RAW_TYPE(md_gpu_kernel_t) kernels;
    MD_ARRAY_RAW_TYPE(md_gpu_pipeline_t) pipelines;
    MD_ARRAY_RAW_TYPE(md_gpu_surface_t) surfaces;
    MD_ARRAY_RAW_TYPE(md_gpu_stream_t) streams;
    md_gpu_stream_t default_compute;
    md_gpu_stream_t default_transfer;
    md_gpu_stream_t default_graphics;   /* NULL without graphics support */

    MD_ARRAY_RAW_TYPE(md_vk_hostfn_t) hostfns;
    MD_ARRAY_RAW_TYPE(md_vk_retire_t) retires;

    md_gpu_kernel_t make_grid_kernel;

    md_vk_dev_caps_t caps;

    bool            is_discrete;
    uint32_t        adapter_index;
    uint64_t        warned_storage_read;   /* md_gpu_format_t bits already warned about */
    bool            validation;
    bool            debug_utils;        /* VK_EXT_debug_utils enabled: object names and labels for GPU tools */
    bool            supports_graphics;
    bool            supports_present;   /* surface + swapchain extensions enabled */
    bool            depth_bias_clamp;   /* feature; otherwise the clamp must be 0 */
    bool            pipeline_stats;     /* VK_KHR_pipeline_executable_properties: kernels capture statistics */
    bool            has_xlib_surface, has_win32_surface, has_wayland_surface, has_headless_surface;
} md_gpu_device;

/* Forward declarations. */
static bool     md_vk_stream_ensure_cmd(md_gpu_stream_t s);
static bool     md_vk_stream_submit(md_gpu_stream_t s);
static uint64_t md_vk_stream_completed(md_gpu_stream_t s);
static bool     md_vk_arena_alloc(md_gpu_stream_t s, size_t size, uint64_t* out_addr, void** out_host,
                                  VkBuffer* out_buffer, uint64_t* out_offset);
static void     md_vk_texture_free(md_gpu_device_t dev, md_gpu_texture_t t);
static void     md_vk_kernel_free(md_gpu_device_t dev, md_gpu_kernel_t k);
static void     md_vk_pipeline_free(md_gpu_device_t dev, md_gpu_pipeline_t p);
static void     md_vk_swapchain_free(md_gpu_device_t dev, md_vk_swapchain_t* sc);
static void     md_vk_heap_init(md_gpu_device_t dev, uint64_t cache_limit, bool debug);
static void     md_vk_abandon_pass(md_gpu_stream_t s);

/* =========================================================================
   3. Address lookup
   ========================================================================= */

/* Resolve [addr, addr + size) to the buffer it lives in (md_gpu_heap.h keeps
   the registry). `what` names the calling function for the error message. */
static bool md_vk_resolve(md_gpu_device_t dev, md_gpu_addr_t addr, uint64_t size,
                          md_vk_span_t* out, const char* what) {
    md_gpu_heap_span_t sp;
    if (!md_gpu_heap_resolve(&dev->heap, addr, size, &sp, what)) return false;
    if (out) {
        out->chunk  = (md_vk_chunk_t*)sp.chunk;
        out->buffer = out->chunk->buffer;
        out->offset = sp.offset;
        out->host   = sp.host;
    }
    return true;
}

/* =========================================================================
   4. Raw buffers and transient arenas
   ========================================================================= */

/* The memory types to allocate from, in the order to try them: those with every preferred flag and
   none of the avoided ones, then those with every preferred flag, then the rest with the required
   flags. Each is tried in turn when allocating from the one before fails - a heap can be full (the
   host visible device local one of a discrete GPU without resizable BAR is 256 MB, shared with the
   driver) while another one serves. Returns the count. */
static uint32_t md_vk_memory_types(md_gpu_device_t dev, uint32_t type_bits, VkMemoryPropertyFlags required, VkMemoryPropertyFlags preferred,
                                   VkMemoryPropertyFlags avoid, uint32_t out[VK_MAX_MEMORY_TYPES]) {
    uint32_t n = 0, taken = 0;
    for (int pass = 0; pass < 3; ++pass) {
        for (uint32_t i = 0; i < dev->mem_props.memoryTypeCount; ++i) {
            if (!(type_bits & (1u << i)) || (taken & (1u << i))) continue;
            const VkMemoryPropertyFlags f = dev->mem_props.memoryTypes[i].propertyFlags;
            if ((f & required) != required) continue;
            if (pass < 2 && (f & preferred) != preferred) continue;
            if (pass < 1 && (f & avoid)) continue;
            out[n++] = i;
            taken |= 1u << i;
        }
    }
    return n;
}

/* Device local memory of its own for `image`, bound. */
static bool md_vk_image_memory(md_gpu_device_t dev, VkImage image, VkDeviceMemory* out_mem, uint64_t* out_size) {
    VkMemoryRequirements req;
    vkGetImageMemoryRequirements(dev->device, image, &req);
    uint32_t types[VK_MAX_MEMORY_TYPES];
    uint32_t n = md_vk_memory_types(dev, req.memoryTypeBits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0, 0, types);
    if (n == 0) n = md_vk_memory_types(dev, req.memoryTypeBits, 0, 0, 0, types);   /* UMA */
    if (n == 0) return md_gpu_fail("no memory type for an image");
    VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.allocationSize  = req.size;
    mai.memoryTypeIndex = types[0];
    if (!md_vk_check(vkAllocateMemory(dev->device, &mai, NULL, out_mem), "vkAllocateMemory (image)")) return false;
    if (out_size) *out_size = req.size;
    return md_vk_check(vkBindImageMemory(dev->device, image, *out_mem, 0), "vkBindImageMemory");
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

/* staging: host memory, not device local, preferred for a host visible kind (what a copy reads once) */
/* Every md_gpu buffer is made alike, so all of them can use the same memory
   types (memoryTypeBits depends only on usage, flags and sharing). */
static VkBufferCreateInfo md_vk_buffer_info(md_gpu_device_t dev, uint64_t size) {
    VkBufferCreateInfo bci = {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bci.size  = size;
    bci.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT
              | VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT
              | VK_BUFFER_USAGE_TRANSFER_SRC_BIT
              | VK_BUFFER_USAGE_TRANSFER_DST_BIT
              | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT
              | VK_BUFFER_USAGE_INDEX_BUFFER_BIT
              | VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    md_vk_set_sharing(dev, &bci.sharingMode, &bci.queueFamilyIndexCount, &bci.pQueueFamilyIndices);
    return bci;
}

/* The memory types MD_GPU_MEM_DEVICE chunks may use. Where some device local
   type also takes optimally tiled images, only those, so that textures can be
   placed in DEVICE memory (md_gpu_texture_create_at). */
static void md_vk_init_device_types(md_gpu_device_t dev) {
    VkBufferCreateInfo bci = md_vk_buffer_info(dev, 65536);
    VkDeviceBufferMemoryRequirements dbmr = {VK_STRUCTURE_TYPE_DEVICE_BUFFER_MEMORY_REQUIREMENTS};
    dbmr.pCreateInfo = &bci;
    VkMemoryRequirements2 bmr = {VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2};
    vkGetDeviceBufferMemoryRequirements(dev->device, &dbmr, &bmr);

    VkImageCreateInfo ici = {VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    ici.imageType   = VK_IMAGE_TYPE_2D;
    ici.format      = VK_FORMAT_R8G8B8A8_UNORM;
    ici.extent      = (VkExtent3D){64, 64, 1};
    ici.mipLevels   = 1;
    ici.arrayLayers = 1;
    ici.samples     = VK_SAMPLE_COUNT_1_BIT;
    ici.tiling      = VK_IMAGE_TILING_OPTIMAL;
    ici.usage       = VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT |
                      VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT;
    md_vk_set_sharing(dev, &ici.sharingMode, &ici.queueFamilyIndexCount, &ici.pQueueFamilyIndices);
    VkDeviceImageMemoryRequirements dimr = {VK_STRUCTURE_TYPE_DEVICE_IMAGE_MEMORY_REQUIREMENTS};
    dimr.pCreateInfo = &ici;
    VkMemoryRequirements2 imr = {VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2};
    vkGetDeviceImageMemoryRequirements(dev->device, &dimr, &imr);

    const uint32_t buf  = bmr.memoryRequirements.memoryTypeBits;
    const uint32_t both = buf & imr.memoryRequirements.memoryTypeBits;
    uint32_t types[VK_MAX_MEMORY_TYPES];
    dev->device_buffer_types = md_vk_memory_types(dev, both, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 0, 0, types) ? both : buf;
}

/* A VkBuffer with VkDeviceMemory of its own. */
typedef struct md_vk_raw_t {
    VkBuffer       buffer;
    VkDeviceMemory memory;
    uint32_t       memory_type;
    uint64_t       memory_offset;   /* of the buffer in `memory`  */
    uint64_t       address;
    uint8_t*       host;            /* of the buffer, or NULL     */
} md_vk_raw_t;

enum {
    MD_VK_RAW_STAGING = 1u << 0,    /* host memory, not device local, preferred for a host visible kind
                                       (what a copy reads once) */
    MD_VK_RAW_GUARDED = 1u << 1,    /* MD_VK_DEBUG_GUARD bytes of memory no buffer covers on each side */
};

static bool md_vk_raw_create(md_gpu_device_t dev, uint64_t size, md_gpu_mem_kind_t kind, uint32_t flags, md_vk_raw_t* out) {
    const bool staging = (flags & MD_VK_RAW_STAGING) != 0;
    VkBufferCreateInfo bci = md_vk_buffer_info(dev, size);
    VkBuffer buf = VK_NULL_HANDLE;
    if (!md_vk_check(vkCreateBuffer(dev->device, &bci, NULL, &buf), "vkCreateBuffer")) return false;

    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(dev->device, buf, &req);
    if (kind == MD_GPU_MEM_DEVICE && (req.memoryTypeBits & dev->device_buffer_types)) {
        req.memoryTypeBits &= dev->device_buffer_types;
    }

    const bool host_visible = kind != MD_GPU_MEM_DEVICE;
    VkMemoryPropertyFlags required = 0, preferred = 0, avoid = 0;
    if (host_visible) {
        required = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        if (kind == MD_GPU_MEM_HOST_READ) preferred = required | VK_MEMORY_PROPERTY_HOST_CACHED_BIT;
        else if (staging)                 { preferred = required; avoid = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT; }
        else                              preferred = required | VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    } else {
        required = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
    }

    uint32_t types[VK_MAX_MEMORY_TYPES];
    uint32_t num_types = md_vk_memory_types(dev, req.memoryTypeBits, required, preferred, avoid, types);
    if (num_types == 0 && !host_visible) {
        /* UMA parts may expose no purely device-local type. */
        num_types = md_vk_memory_types(dev, req.memoryTypeBits, 0, 0, 0, types);
    }
    if (num_types == 0) {
        vkDestroyBuffer(dev->device, buf, NULL);
        return md_gpu_fail("no suitable memory type for %llu bytes (kind %d)", (unsigned long long)size, (int)kind);
    }

    VkMemoryAllocateFlagsInfo fi = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO};
    fi.flags = VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT;

    /* Guarded: the buffer sits between two stretches of its memory that no
       buffer covers, so that GPU-assisted validation, which checks pointer
       accesses against buffer ranges, sees an access past either end. */
    const uint64_t guard = (flags & MD_VK_RAW_GUARDED) ? MD_ALIGN_TO(MD_VK_DEBUG_GUARD, req.alignment) : 0;
    VkMemoryAllocateInfo mai = {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    mai.pNext           = &fi;
    mai.allocationSize  = guard ? guard + MD_ALIGN_TO(req.size, req.alignment) + guard : req.size;

    /* the first type whose heap has room */
    VkDeviceMemory mem = VK_NULL_HANDLE;
    VkResult res = VK_ERROR_OUT_OF_DEVICE_MEMORY;
    for (uint32_t t = 0; t < num_types && mem == VK_NULL_HANDLE; ++t) {
        mai.memoryTypeIndex = types[t];
        out->memory_type    = types[t];
        res = vkAllocateMemory(dev->device, &mai, NULL, &mem);
        if (res != VK_SUCCESS) mem = VK_NULL_HANDLE;
        if (res != VK_SUCCESS && res != VK_ERROR_OUT_OF_DEVICE_MEMORY && res != VK_ERROR_OUT_OF_HOST_MEMORY) break;
    }
    if (mem == VK_NULL_HANDLE) {
        vkDestroyBuffer(dev->device, buf, NULL);
        md_vk_check(res, "vkAllocateMemory");
        return md_gpu_fail("could not allocate %llu bytes (kind %d) from any of %u memory types", (unsigned long long)size, (int)kind, num_types);
    }
    if (!md_vk_check(vkBindBufferMemory(dev->device, buf, mem, guard), "vkBindBufferMemory")) {
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
    out->buffer        = buf;
    out->memory        = mem;
    out->memory_offset = guard;
    out->address       = vkGetBufferDeviceAddress(dev->device, &bdai);
    out->host          = host ? (uint8_t*)host + guard : NULL;
    return true;
}

static void md_vk_destroy_raw_buffer(md_gpu_device_t dev, VkBuffer buf, VkDeviceMemory mem) {
    if (buf) vkDestroyBuffer(dev->device, buf, NULL);
    if (mem) vkFreeMemory(dev->device, mem, NULL);
}

static md_vk_page_t* md_vk_page_create(md_gpu_device_t dev, uint64_t size, bool dedicated) {
    md_vk_page_t* p = (md_vk_page_t*)md_alloc(dev->alloc, sizeof(md_vk_page_t));
    if (!p) return NULL;
    memset(p, 0, sizeof(*p));
    p->dedicated = dedicated;
    md_vk_raw_t raw;
    if (!md_vk_raw_create(dev, size, MD_GPU_MEM_HOST_WRITE, dedicated ? MD_VK_RAW_STAGING : 0, &raw)) {
        md_free(dev->alloc, p, sizeof(md_vk_page_t));
        return NULL;
    }
    p->buffer   = raw.buffer;
    p->memory   = raw.memory;
    p->address  = raw.address;
    p->host     = raw.host;
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
    uint64_t need = MD_ALIGN_TO(size, MD_VK_ARG_ALIGN);

    /* Try the current page. */
    if (a->current < a->pages.size) {
        md_vk_page_t* p = a->pages.data[a->current];
        if (!p->dedicated && p->cursor + need <= p->capacity) {
            md_vk_page_take(p, need, out_addr, out_host, out_buffer, out_offset);
            return true;
        }
    }

    /* Look for a retired page with room; a retired dedicated page is freed. */
    uint64_t done = md_vk_stream_completed(s);
    for (size_t i = 0; i < a->pages.size;) {
        md_vk_page_t* p = a->pages.data[i];
        const bool retired = p->retire_value != 0 && p->retire_value <= done;
        if (p->dedicated) {
            if (retired) {
                md_vk_page_destroy(dev, p);
                MD_ARRAY_REMOVE(a->pages, i);
                if (a->current > i || a->current >= a->pages.size) a->current = a->current > 0 ? a->current - 1 : 0;
            } else {
                ++i;
            }
            continue;
        }
        if (retired) {
            p->cursor = 0;
            p->retire_value = 0;
        }
        if (p->retire_value == 0 && p->cursor + need <= p->capacity && need <= MD_VK_ARENA_DEDICATED) {
            a->current = i;
            md_vk_page_take(p, need, out_addr, out_host, out_buffer, out_offset);
            return true;
        }
        ++i;
    }

    /* Allocate a new page: shared, or the allocation's own. */
    const bool dedicated = need > MD_VK_ARENA_DEDICATED;
    uint64_t page_size = MD_VK_ARENA_PAGE_SIZE;
    if (dedicated) page_size = MD_ALIGN_TO(need, 65536);
    else while (page_size < need) page_size *= 2;
    md_vk_page_t* p = md_vk_page_create(dev, page_size, dedicated);
    if (!p) return md_gpu_fail("failed to allocate a %llu byte transient page", (unsigned long long)page_size);
    if (!MD_ARRAY_TRY_PUSH(a->pages, p)) { md_vk_page_destroy(dev, p); return md_gpu_fail("out of memory"); }
    if (!dedicated) a->current = a->pages.size - 1;
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
    for (size_t i = 0; i < a->pages.size; ++i) {
        md_vk_page_t* p = a->pages.data[i];
        if (p->cursor > 0) p->retire_value = value;
    }
}

static void md_vk_arena_free(md_gpu_device_t dev, md_vk_arena_t* a) {
    for (size_t i = 0; i < a->pages.size; ++i) {
        md_vk_page_destroy(dev, a->pages.data[i]);
    }
    MD_ARRAY_FREE(a->pages);
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

    uint32_t n = (uint32_t)dev->streams.size;
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

    if (!((n == 0 || r.waits) && MD_ARRAY_TRY_PUSH(dev->retires, r))) {
        md_gpu_fail("out of memory recording a deferred destruction; waiting for the device");
        if (r.waits) md_free(dev->alloc, r.waits, n * sizeof(md_vk_wait_t));
        vkDeviceWaitIdle(dev->device);
        switch (kind) {
        case MD_VK_RETIRE_TEXTURE: md_vk_texture_free(dev, (md_gpu_texture_t)object); break;
        case MD_VK_RETIRE_KERNEL:  md_vk_kernel_free(dev, (md_gpu_kernel_t)object);   break;
        case MD_VK_RETIRE_PIPELINE:  md_vk_pipeline_free(dev, (md_gpu_pipeline_t)object); break;
        case MD_VK_RETIRE_SWAPCHAIN: md_vk_swapchain_free(dev, (md_vk_swapchain_t*)object); break;
        }
        return;
    }
}

/* Release every retire entry whose waits have all completed (all of them when
   `force`). Caller holds device_mutex. */
static void md_vk_process_retires_locked(md_gpu_device_t dev, bool force) {
    for (size_t i = 0; i < dev->retires.size;) {
        md_vk_retire_t* r = &dev->retires.data[i];
        bool done = true;
        for (uint32_t w = 0; w < r->wait_count && done && !force; ++w) {
            if (md_vk_stream_completed(r->waits[w].stream) < r->waits[w].value) done = false;
        }
        if (!done) { ++i; continue; }
        md_vk_retire_t e = *r;
        MD_ARRAY_REMOVE(dev->retires, i);
        switch (e.kind) {
        case MD_VK_RETIRE_TEXTURE: md_vk_texture_free(dev, (md_gpu_texture_t)e.object); break;
        case MD_VK_RETIRE_KERNEL:  md_vk_kernel_free(dev, (md_gpu_kernel_t)e.object);   break;
        case MD_VK_RETIRE_PIPELINE:  md_vk_pipeline_free(dev, (md_gpu_pipeline_t)e.object); break;
        case MD_VK_RETIRE_SWAPCHAIN: md_vk_swapchain_free(dev, (md_vk_swapchain_t*)e.object); break;
        }
        if (e.waits) md_free(dev->alloc, e.waits, e.wait_capacity * sizeof(md_vk_wait_t));
    }
}

/* A stream is going away after having been synchronised: every wait on it is
   satisfied, so drop the references rather than leave them dangling. Also
   releases memory it freed -- the free points would otherwise never be
   reached. Caller holds device_mutex. */
static void md_vk_forget_stream_locked(md_gpu_device_t dev, md_gpu_stream_t s) {
    md_gpu_heap_forget_stream_locked(&dev->heap, s);
    /* Clearing the sync makes a pending callback unconditionally ready, so it
       still fires on the next md_gpu_device_poll, on the polling thread. */
    for (size_t i = 0; i < dev->hostfns.size; ++i) {
        md_vk_hostfn_t* h = &dev->hostfns.data[i];
        if (h->sync.stream == s) h->sync = md_gpu_sync_none();
    }
    for (size_t i = 0; i < dev->retires.size; ++i) {
        md_vk_retire_t* r = &dev->retires.data[i];
        for (uint32_t w = 0; w < r->wait_count;) {
            if (r->waits[w].stream == s) r->waits[w] = r->waits[--r->wait_count];
            else ++w;
        }
    }
    for (size_t i = 0; i < dev->surfaces.size; ++i) {
        md_gpu_surface_t sf = dev->surfaces.data[i];
        for (uint32_t k = 0; k < MD_VK_ACQUIRE_SEMS; ++k) {
            if (sf->acquire_stream[k] == s) { sf->acquire_stream[k] = NULL; sf->acquire_value[k] = 0; }
        }
        if (sf->acquired_stream == s) { sf->acquired = false; sf->acquired_stream = NULL; }
    }
    for (size_t i = 0; i < dev->streams.size; ++i) {
        md_gpu_stream_t o = dev->streams.data[i];
        for (size_t w = 0; w < o->waits.size;) {
            if (o->waits.data[w].stream == s) MD_ARRAY_REMOVE(o->waits, w);
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

static bool md_vk_device_ext_available(VkPhysicalDevice pd, const char* name) {
    uint32_t n = 0;
    vkEnumerateDeviceExtensionProperties(pd, NULL, &n, NULL);
    if (n == 0) return false;
    VkExtensionProperties* props = (VkExtensionProperties*)malloc(n * sizeof(VkExtensionProperties));
    if (!props) return false;
    vkEnumerateDeviceExtensionProperties(pd, NULL, &n, props);
    bool found = false;
    for (uint32_t i = 0; i < n && !found; ++i) {
        if (strcmp(props[i].extensionName, name) == 0) found = true;
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


/* On failure, `missing` lists every unmet requirement (comma separated), so a
   device can be judged in one go rather than one feature at a time. */
static bool md_vk_probe_device(VkPhysicalDevice pd, struct md_allocator_i* alloc,
                               md_vk_dev_caps_t* out_caps, char* missing, size_t missing_cap)
{
    memset(out_caps, 0, sizeof(*out_caps));
    missing[0] = 0;
    size_t missing_len = 0;
    bool ok = true;
#define MD_VK_MISSING(name) do {                                                              \
        ok = false;                                                                            \
        if (missing_len + 1 < missing_cap) {                                                   \
            int w_ = md_snprintf(missing + missing_len, missing_cap - missing_len, "%s%s",        \
                              missing_len ? ", " : "", (name));                                \
            if (w_ > 0) missing_len = MD_MIN(missing_len + (size_t)w_, missing_cap - 1);           \
        }                                                                                      \
    } while (0)
#define MD_VK_REQUIRE(cond, name) do { if (!(cond)) MD_VK_MISSING(name); } while (0)

    VkPhysicalDeviceProperties props;
    vkGetPhysicalDeviceProperties(pd, &props);
    /* We call the core 1.3 synchronization2 and dynamic-rendering entry points
       directly rather than their KHR aliases. */
    if (props.apiVersion < VK_API_VERSION_1_3) {
        char v[64];
        md_snprintf(v, sizeof(v), "Vulkan 1.3 (driver has %u.%u.%u)", VK_API_VERSION_MAJOR(props.apiVersion),
                 VK_API_VERSION_MINOR(props.apiVersion), VK_API_VERSION_PATCH(props.apiVersion));
        MD_VK_MISSING(v);
    }

    /* No device extensions are required: everything md_gpu uses is core 1.3. */
    (void)alloc;

    /* The 1.2 / 1.3 feature structs may only be chained on devices of that version. */
    VkPhysicalDeviceVulkan13Features f13 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES};
    VkPhysicalDeviceVulkan12Features f12 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES};
    VkPhysicalDeviceVulkan11Features f11 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES};
    VkPhysicalDeviceFeatures2        f2  = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    const bool has12 = props.apiVersion >= VK_API_VERSION_1_2;
    const bool has13 = props.apiVersion >= VK_API_VERSION_1_3;
    if (has12) { f2.pNext = &f11; f11.pNext = &f12; if (has13) f12.pNext = &f13; }
    vkGetPhysicalDeviceFeatures2(pd, &f2);

    if (has12) {
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
    }
    if (has13) {
        MD_VK_REQUIRE(f13.synchronization2,      "synchronization2");
    }
    /* The heap arrays carry no format qualifier. The device-wide features say
       "every format in the storage-without-format list"; from Vulkan 1.3 on the
       SPIR-V capabilities are allowed without them and support is per format
       (VK_FORMAT_FEATURE_2_STORAGE_{READ,WRITE}_WITHOUT_FORMAT_BIT), which
       md_gpu_texture_create checks. Older Intel drivers (Gen9 on Windows) lack
       the device-wide read feature but can read R32 formats. */
    out_caps->storage_read_without_format  = f2.features.shaderStorageImageReadWithoutFormat;
    out_caps->storage_write_without_format = f2.features.shaderStorageImageWriteWithoutFormat;
#undef MD_VK_REQUIRE
#undef MD_VK_MISSING
    if (!ok) return false;

    /* Rendering: optional as a whole. Vertex pulling needs draw parameters
       (Slang's SV_VertexID subtracts the base vertex), multi-draw indirect
       with first_instance carries per-draw records, and per-target blend
       needs independentBlend. Dynamic rendering and the dynamic depth/cull
       state are core 1.3 and need no feature bit beyond dynamicRendering. */
    {
        VkPhysicalDeviceVulkan11Features f11 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES};
        VkPhysicalDeviceFeatures2        g2  = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, &f11};
        vkGetPhysicalDeviceFeatures2(pd, &g2);
        uint32_t qn = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qn, NULL);
        VkQueueFamilyProperties qp[16];
        if (qn > 16) qn = 16;
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qn, qp);
        bool has_gfx = false;
        for (uint32_t i = 0; i < qn; ++i) if (qp[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) has_gfx = true;
        out_caps->graphics = has_gfx && f13.dynamicRendering && f11.shaderDrawParameters &&
                             g2.features.multiDrawIndirect && g2.features.drawIndirectFirstInstance &&
                             g2.features.independentBlend;
        out_caps->depth_bias_clamp = g2.features.depthBiasClamp;
    }

    out_caps->maintenance4                = f13.maintenance4;
    out_caps->update_unused_while_pending = f12.descriptorBindingUpdateUnusedWhilePending;
    out_caps->nonuniform_storage_image    = f12.shaderStorageImageArrayNonUniformIndexing;
    out_caps->nonuniform_sampled_image    = f12.shaderSampledImageArrayNonUniformIndexing;
    out_caps->dynamic_storage_image       = f2.features.shaderStorageImageArrayDynamicIndexing;
    out_caps->dynamic_sampled_image       = f2.features.shaderSampledImageArrayDynamicIndexing;
    out_caps->shader_int64                = f2.features.shaderInt64;

    /* Above the floor: optional, enabled where present. */
    out_caps->float16              = f12.shaderFloat16 && f11.storageBuffer16BitAccess;
    out_caps->float64              = f2.features.shaderFloat64;
    out_caps->int64_atomics        = f2.features.shaderInt64 && f12.shaderBufferInt64Atomics;
    out_caps->int64_shared_atomics = out_caps->int64_atomics && f12.shaderSharedInt64Atomics;
    if (md_vk_device_ext_available(pd, VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME)) {
        VkPhysicalDeviceShaderAtomicFloatFeaturesEXT af = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT};
        VkPhysicalDeviceFeatures2 q = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, &af};
        vkGetPhysicalDeviceFeatures2(pd, &q);
        af.pNext = NULL;
        out_caps->atomic_float  = af;
        out_caps->float_atomics = af.shaderBufferFloat32Atomics || af.shaderBufferFloat32AtomicAdd ||
                                  af.shaderSharedFloat32Atomics || af.shaderSharedFloat32AtomicAdd;
    }
    if (md_vk_device_ext_available(pd, VK_EXT_SHADER_ATOMIC_FLOAT_2_EXTENSION_NAME)) {
        VkPhysicalDeviceShaderAtomicFloat2FeaturesEXT af = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_2_FEATURES_EXT};
        VkPhysicalDeviceFeatures2 q = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, &af};
        vkGetPhysicalDeviceFeatures2(pd, &q);
        af.pNext = NULL;
        out_caps->atomic_float2  = af;
        out_caps->float_atomics2 = af.shaderBufferFloat32AtomicMinMax || af.shaderSharedFloat32AtomicMinMax;
    }
    return true;
}

static bool md_vk_create_bindless(md_gpu_device_t dev);
static bool md_vk_immediate(md_gpu_device_t dev, void (*record)(VkCommandBuffer, void*), void* user);
static void md_vk_record_to_general(VkCommandBuffer cmd, void* user);
static bool md_vk_create_dummies(md_gpu_device_t dev);
static md_gpu_stream_t md_vk_stream_create_internal(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label, bool is_default);
static bool md_vk_create_builtin_kernels(md_gpu_device_t dev);

typedef struct md_vk_adapter_t {
    VkPhysicalDevice      pd;
    md_vk_dev_caps_t      caps;
    md_gpu_adapter_info_t info;
} md_vk_adapter_t;

static md_gpu_device_type_t md_vk_device_type(VkPhysicalDeviceType t) {
    switch (t) {
    case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:   return MD_GPU_DEVICE_TYPE_DISCRETE;
    case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return MD_GPU_DEVICE_TYPE_INTEGRATED;
    case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:    return MD_GPU_DEVICE_TYPE_VIRTUAL;
    case VK_PHYSICAL_DEVICE_TYPE_CPU:            return MD_GPU_DEVICE_TYPE_CPU;
    default:                                     return MD_GPU_DEVICE_TYPE_OTHER;
    }
}

static void md_vk_driver_string(VkPhysicalDevice pd, const VkPhysicalDeviceProperties* p, char* out, size_t cap) {
    if (p->apiVersion >= VK_API_VERSION_1_2) {
        VkPhysicalDeviceDriverProperties drv = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES};
        VkPhysicalDeviceProperties2 p2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2, &drv};
        vkGetPhysicalDeviceProperties2(pd, &p2);
        md_snprintf(out, cap, "%s %s (Vulkan %u.%u.%u)", drv.driverName, drv.driverInfo,
                 VK_API_VERSION_MAJOR(p->apiVersion), VK_API_VERSION_MINOR(p->apiVersion), VK_API_VERSION_PATCH(p->apiVersion));
    } else {
        md_snprintf(out, cap, "driver version 0x%08X (Vulkan %u.%u.%u)", p->driverVersion,
                 VK_API_VERSION_MAJOR(p->apiVersion), VK_API_VERSION_MINOR(p->apiVersion), VK_API_VERSION_PATCH(p->apiVersion));
    }
}

/* Every physical device of `instance`, probed. The caller frees *out with
   md_free(alloc, *out, count * sizeof(md_vk_adapter_t)). Needs the instance's
   functions loaded through volk. */
static uint32_t md_vk_collect_adapters(VkInstance instance, struct md_allocator_i* alloc, md_vk_adapter_t** out) {
    *out = NULL;
    uint32_t n = 0;
    if (vkEnumeratePhysicalDevices(instance, &n, NULL) != VK_SUCCESS || n == 0) return 0;
    VkPhysicalDevice* pds = (VkPhysicalDevice*)md_alloc(alloc, n * sizeof(VkPhysicalDevice));
    md_vk_adapter_t*  ad  = (md_vk_adapter_t*)md_alloc(alloc, n * sizeof(md_vk_adapter_t));
    if (!pds || !ad) {
        if (pds) md_free(alloc, pds, n * sizeof(VkPhysicalDevice));
        if (ad)  md_free(alloc, ad,  n * sizeof(md_vk_adapter_t));
        return 0;
    }
    vkEnumeratePhysicalDevices(instance, &n, pds);
    memset(ad, 0, n * sizeof(md_vk_adapter_t));
    for (uint32_t i = 0; i < n; ++i) {
        VkPhysicalDeviceProperties p;
        vkGetPhysicalDeviceProperties(pds[i], &p);
        md_gpu_adapter_info_t* info = &ad[i].info;
        ad[i].pd = pds[i];
        md_snprintf(info->name, sizeof(info->name), "%s", p.deviceName);
        info->vendor_id = p.vendorID;
        info->device_id = p.deviceID;
        info->type      = md_vk_device_type(p.deviceType);
        md_vk_driver_string(pds[i], &p, info->driver, sizeof(info->driver));
        info->usable = md_vk_probe_device(pds[i], alloc, &ad[i].caps, info->missing, sizeof(info->missing));
        if (!info->usable && !info->missing[0]) md_snprintf(info->missing, sizeof(info->missing), "?");
    }
    md_free(alloc, pds, n * sizeof(VkPhysicalDevice));
    *out = ad;
    return n;
}

uint32_t md_gpu_enumerate_adapters(md_gpu_adapter_info_t* out, uint32_t max) {
    md_gpu_clear_error();
    if (volkInitialize() != VK_SUCCESS) {
        md_gpu_fail("volkInitialize failed — no Vulkan loader present");
        return 0;
    }
    VkApplicationInfo ai = {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = "mdlib";
    ai.apiVersion       = VK_API_VERSION_1_3;
    VkInstanceCreateInfo ici = {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ici.pApplicationInfo = &ai;
    VkInstance instance = VK_NULL_HANDLE;
    if (!md_vk_check(vkCreateInstance(&ici, NULL, &instance), "vkCreateInstance")) return 0;

    /* volk keeps one set of instance-level entry points; put back those of a
       live device once done. */
    const VkInstance prev = md_vk_live_instance;
    volkLoadInstanceOnly(instance);

    struct md_allocator_i* alloc = md_get_heap_allocator();
    md_vk_adapter_t* ad = NULL;
    const uint32_t n = md_vk_collect_adapters(instance, alloc, &ad);
    for (uint32_t i = 0; i < n && i < max; ++i) out[i] = ad[i].info;
    if (ad) md_free(alloc, ad, n * sizeof(md_vk_adapter_t));

    vkDestroyInstance(instance, NULL);
    if (prev) volkLoadInstanceOnly(prev);
    return n;
}

md_gpu_device_t md_gpu_device_create(const md_gpu_device_desc_t* desc) {
    md_gpu_clear_error();

    struct md_allocator_i* alloc = (desc && desc->alloc) ? desc->alloc : md_get_heap_allocator();

    if (volkInitialize() != VK_SUCCESS) {
        md_gpu_fail("volkInitialize failed — no Vulkan loader present");
        return NULL;
    }

    md_gpu_device_t dev = (md_gpu_device_t)md_alloc(alloc, sizeof(md_gpu_device));
    if (!dev) { md_gpu_fail("out of memory"); return NULL; }
    memset(dev, 0, sizeof(*dev));
    dev->alloc = alloc;
    dev->validation = desc && desc->enable_validation;
    {
        const char* env = md_os_env_get("MD_GPU_DEBUG_MEMORY");
        const bool debug_memory = (desc && desc->debug_memory) || (env && env[0] && strcmp(env, "0") != 0);
        md_vk_heap_init(dev, desc ? desc->heap_cache_limit : 0, debug_memory);
    }

    MD_ARRAY_INIT(dev->textures, alloc);
    MD_ARRAY_INIT(dev->kernels, alloc);
    MD_ARRAY_INIT(dev->pipelines, alloc);
    MD_ARRAY_INIT(dev->surfaces, alloc);
    MD_ARRAY_INIT(dev->streams, alloc);
    MD_ARRAY_INIT(dev->hostfns, alloc);
    MD_ARRAY_INIT(dev->retires, alloc);

    /* ---- instance ---- */
    VkApplicationInfo ai = {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = (desc && desc->label) ? desc->label : "mdlib";
    ai.apiVersion       = VK_API_VERSION_1_3;

    const char* layers[4];     uint32_t layer_count = 0;
    const char* exts[8];       uint32_t ext_count   = 0;

    /* Presentation: the surface extensions that exist. None of them is
       required; without VK_KHR_surface md_gpu simply cannot present. */
    if (md_vk_ext_available(VK_KHR_SURFACE_EXTENSION_NAME)) {
        exts[ext_count++] = VK_KHR_SURFACE_EXTENSION_NAME;
        if (md_vk_ext_available("VK_KHR_win32_surface"))   { exts[ext_count++] = "VK_KHR_win32_surface";   dev->has_win32_surface   = true; }
        if (md_vk_ext_available("VK_KHR_xlib_surface"))    { exts[ext_count++] = "VK_KHR_xlib_surface";    dev->has_xlib_surface    = true; }
        if (md_vk_ext_available("VK_KHR_wayland_surface")) { exts[ext_count++] = "VK_KHR_wayland_surface"; dev->has_wayland_surface = true; }
        if (md_vk_ext_available(VK_EXT_HEADLESS_SURFACE_EXTENSION_NAME)) {
            exts[ext_count++] = VK_EXT_HEADLESS_SURFACE_EXTENSION_NAME;
            dev->has_headless_surface = true;
        }
    }
    const uint32_t surface_ext_count = ext_count;

    /* VK_EXT_debug_utils whenever the loader offers it, validation or not: kernels and pipelines get
       their labels as object names and every dispatch a command-buffer label, which is what profilers
       and capture tools (Nsight Systems / Graphics, RenderDoc) show. Nothing listens otherwise. */
    const bool has_debug_utils = md_vk_ext_available(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
    if (has_debug_utils) exts[ext_count++] = VK_EXT_DEBUG_UTILS_EXTENSION_NAME;
    const uint32_t base_ext_count = ext_count;

    bool want_debug = dev->validation
        && md_vk_layer_available("VK_LAYER_KHRONOS_validation")
        && has_debug_utils;
    if (want_debug) {
        layers[layer_count++] = "VK_LAYER_KHRONOS_validation";
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
        /* Retry without validation, then without debug utils too. */
        ici.enabledLayerCount = 0;
        want_debug = false;
        if (vkCreateInstance(&ici, NULL, &dev->instance) != VK_SUCCESS) {
            ici.enabledExtensionCount = surface_ext_count;
            if (!md_vk_check(vkCreateInstance(&ici, NULL, &dev->instance), "vkCreateInstance")) {
                md_free(alloc, dev, sizeof(*dev));
                return NULL;
            }
        }
    }
    dev->debug_utils = has_debug_utils && ici.enabledExtensionCount == base_ext_count;
    volkLoadInstance(dev->instance);
    MD_LOG_DEBUG("md_gpu: VK_EXT_debug_utils %s", dev->debug_utils ? "enabled: kernels and dispatches are named for GPU tools" : "unavailable");
    md_vk_live_instance = dev->instance;

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
    md_vk_adapter_t* adapters = NULL;
    const uint32_t adapter_count = md_vk_collect_adapters(dev->instance, alloc, &adapters);
    if (adapter_count == 0) {
        md_gpu_fail("no Vulkan physical devices");
        goto fail_instance;
    }
    md_gpu_adapter_info_t* infos = (md_gpu_adapter_info_t*)md_alloc(alloc, adapter_count * sizeof(md_gpu_adapter_info_t));
    for (uint32_t i = 0; i < adapter_count; ++i) {
        infos[i] = adapters[i].info;
        MD_LOG_DEBUG("md_gpu: adapter [%u] '%s' (%s)%s%s", i, infos[i].name, md_gpu_sel_type_str(infos[i].type),
                     infos[i].usable ? "" : " — lacks ", infos[i].usable ? "" : infos[i].missing);
    }
    char why[2048];
    const int pick = md_gpu_sel_pick(infos, adapter_count, desc, why, sizeof(why));
    md_free(alloc, infos, adapter_count * sizeof(md_gpu_adapter_info_t));

    VkPhysicalDevice chosen = VK_NULL_HANDLE;
    md_vk_dev_caps_t caps = {0};
    if (pick >= 0) {
        chosen = adapters[pick].pd;
        caps   = adapters[pick].caps;
        dev->adapter_index = (uint32_t)pick;
    }
    md_free(alloc, adapters, adapter_count * sizeof(md_vk_adapter_t));
    if (!chosen) {
        md_gpu_fail("md_gpu_device_create: %s", why);
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
        dev->subgroup_ops  = (sgp.supportedStages & VK_SHADER_STAGE_COMPUTE_BIT) ? sgp.supportedOperations : 0;
        md_vk_driver_string(dev->phys, &dev->props, dev->driver_desc, sizeof(dev->driver_desc));
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
    /* The universal family, for GRAPHICS streams. */
    dev->graphics_family = UINT32_MAX;
    if (caps.graphics) {
        for (uint32_t i = 0; i < qf_count; ++i) {
            if (qfs[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) { dev->graphics_family = i; break; }
        }
    }
    if (dev->graphics_family != UINT32_MAX) {
        dev->graphics_queue_count = qfs[dev->graphics_family].queueCount;
        if (dev->graphics_queue_count > MD_VK_MAX_QUEUES_PER_FAMILY) dev->graphics_queue_count = MD_VK_MAX_QUEUES_PER_FAMILY;
        if (dev->graphics_queue_count == 0) dev->graphics_queue_count = 1;
    }
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
        md_gpu_fail("no compute-capable queue family");
        goto fail_instance;
    }
    /* What each family's timestamps hold: the GPU time of a submission
       (md_gpu_sync_gpu_seconds), where the family and the device can time */
    dev->ts_period = dev->props.limits.timestampPeriod;
    if (dev->ts_period > 0.0f) {
        dev->ts_bits_compute  = qfs[dev->compute_family].timestampValidBits;
        dev->ts_bits_transfer = dev->transfer_family < qf_count ? qfs[dev->transfer_family].timestampValidBits : 0;
        dev->ts_bits_graphics = dev->graphics_family < qf_count ? qfs[dev->graphics_family].timestampValidBits : 0;
    }
    md_free(alloc, qfs, qf_count * sizeof(VkQueueFamilyProperties));

    /* ---- logical device ---- */
    static const float prios[MD_VK_MAX_QUEUES_PER_FAMILY] = {1,1,1,1,1,1,1,1};
    VkDeviceQueueCreateInfo qci[3];
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
    if (dev->graphics_family != UINT32_MAX && dev->graphics_family != dev->compute_family &&
        dev->graphics_family != dev->transfer_family) {
        qci[qci_count] = (VkDeviceQueueCreateInfo){VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        qci[qci_count].queueFamilyIndex = dev->graphics_family;
        qci[qci_count].queueCount       = dev->graphics_queue_count;
        qci[qci_count].pQueuePriorities = prios;
        qci_count++;
    }

    /* Everything below was confirmed present by md_vk_probe_device. Enabling a
       feature the driver does not report is VK_ERROR_FEATURE_NOT_PRESENT, so
       the optional ones are gated on the probe's answer. */
    VkPhysicalDeviceVulkan13Features f13 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_3_FEATURES};
    f13.synchronization2 = VK_TRUE;
    f13.maintenance4     = caps.maintenance4 ? VK_TRUE : VK_FALSE;
    f13.dynamicRendering = caps.graphics ? VK_TRUE : VK_FALSE;

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
    f12.shaderFloat16                             = caps.float16                     ? VK_TRUE : VK_FALSE;
    f12.shaderBufferInt64Atomics                  = caps.int64_atomics               ? VK_TRUE : VK_FALSE;
    f12.shaderSharedInt64Atomics                  = caps.int64_shared_atomics        ? VK_TRUE : VK_FALSE;

    VkPhysicalDeviceVulkan11Features f11 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES};
    f11.pNext = &f12;
    f11.shaderDrawParameters     = caps.graphics ? VK_TRUE : VK_FALSE;
    f11.storageBuffer16BitAccess = caps.float16  ? VK_TRUE : VK_FALSE;

    VkPhysicalDeviceFeatures2 f2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    f2.pNext = &f11;
    f2.features.multiDrawIndirect         = caps.graphics ? VK_TRUE : VK_FALSE;
    f2.features.drawIndirectFirstInstance = caps.graphics ? VK_TRUE : VK_FALSE;
    f2.features.independentBlend          = caps.graphics ? VK_TRUE : VK_FALSE;
    f2.features.depthBiasClamp            = (caps.graphics && caps.depth_bias_clamp) ? VK_TRUE : VK_FALSE;
    f2.features.largePoints               = VK_FALSE;
    /* Slang's heap arrays carry no format qualifier. */
    f2.features.shaderStorageImageReadWithoutFormat    = caps.storage_read_without_format  ? VK_TRUE : VK_FALSE;
    f2.features.shaderStorageImageWriteWithoutFormat   = caps.storage_write_without_format ? VK_TRUE : VK_FALSE;
    f2.features.shaderStorageImageArrayDynamicIndexing = caps.dynamic_storage_image ? VK_TRUE : VK_FALSE;
    f2.features.shaderSampledImageArrayDynamicIndexing = caps.dynamic_sampled_image ? VK_TRUE : VK_FALSE;
    f2.features.shaderInt64                            = caps.shader_int64          ? VK_TRUE : VK_FALSE;
    f2.features.shaderFloat64                          = caps.float64               ? VK_TRUE : VK_FALSE;

    /* Swapchains, when the instance can make surfaces and the device can
       present through a graphics queue. */
    const char* dev_exts[4];
    uint32_t    dev_ext_count = 0;
    /* Float atomics, exactly as the probe found them. */
    if (caps.float_atomics) {
        dev_exts[dev_ext_count++] = VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME;
        caps.atomic_float.pNext = f13.pNext;
        f13.pNext = &caps.atomic_float;
    }
    if (caps.float_atomics2) {
        dev_exts[dev_ext_count++] = VK_EXT_SHADER_ATOMIC_FLOAT_2_EXTENSION_NAME;
        caps.atomic_float2.pNext = f13.pNext;
        f13.pNext = &caps.atomic_float2;
    }
    if (caps.graphics && surface_ext_count > 0 && md_vk_device_ext_available(dev->phys, VK_KHR_SWAPCHAIN_EXTENSION_NAME)) {
        dev_exts[dev_ext_count++] = VK_KHR_SWAPCHAIN_EXTENSION_NAME;
        dev->supports_present = true;
    }
    /* The compiled kernels' statistics (registers, shared memory, ...), for profiling, where the
       driver reports them: md_gpu_kernel_info's stats. */
    VkPhysicalDevicePipelineExecutablePropertiesFeaturesKHR fpe = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PIPELINE_EXECUTABLE_PROPERTIES_FEATURES_KHR};
    if (md_vk_device_ext_available(dev->phys, VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME)) {
        VkPhysicalDeviceFeatures2 q2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
        q2.pNext = &fpe;
        vkGetPhysicalDeviceFeatures2(dev->phys, &q2);
        if (fpe.pipelineExecutableInfo) {
            fpe.pNext = f13.pNext;
            f13.pNext = &fpe;
            dev_exts[dev_ext_count++] = VK_KHR_PIPELINE_EXECUTABLE_PROPERTIES_EXTENSION_NAME;
            dev->pipeline_stats = true;
        }
    }

    VkDeviceCreateInfo dci = {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.pNext                   = &f2;
    dci.queueCreateInfoCount     = qci_count;
    dci.pQueueCreateInfos        = qci;
    dci.enabledExtensionCount    = dev_ext_count;
    dci.ppEnabledExtensionNames  = dev_exts;

    if (!md_vk_check(vkCreateDevice(dev->phys, &dci, NULL, &dev->device), "vkCreateDevice")) goto fail_instance;
    caps.atomic_float.pNext  = NULL;
    caps.atomic_float2.pNext = NULL;
    dev->caps = caps;
    dev->supports_graphics = dev->graphics_family != UINT32_MAX;
    dev->depth_bias_clamp  = caps.depth_bias_clamp;
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

    if (dev->graphics_family != UINT32_MAX) {
        for (uint32_t i = 0; i < dev->graphics_queue_count; ++i) {
            vkGetDeviceQueue(dev->device, dev->graphics_family, i, &dev->graphics_queues[i]);
        }
    }

    md_mutex_init(&dev->queue_mutex);
    md_mutex_init(&dev->device_mutex);

    dev->share_families[0]  = dev->compute_family;
    dev->share_family_count = 1;
    if (dev->transfer_family != dev->compute_family) {
        dev->share_families[dev->share_family_count++] = dev->transfer_family;
    }
    if (dev->graphics_family != UINT32_MAX && dev->graphics_family != dev->compute_family &&
        dev->graphics_family != dev->transfer_family) {
        dev->share_families[dev->share_family_count++] = dev->graphics_family;
    }

    /* Heap slots: as many as asked for, within what the device can bind
       update-after-bind in one stage. The free list runs high-to-low so slots
       are handed out from 1; slot 0 is never handed out so that a zero handle
       stays the null handle. */
    {
        VkPhysicalDeviceVulkan12Properties p12 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_PROPERTIES};
        VkPhysicalDeviceProperties2 p2 = {VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2, &p12};
        vkGetPhysicalDeviceProperties2(dev->phys, &p2);
        uint32_t want = (desc && desc->max_textures) ? desc->max_textures + 1 : MD_VK_DEFAULT_TEXTURE_SLOTS;
        uint32_t lim  = p12.maxPerStageDescriptorUpdateAfterBindSampledImages;
        if (p12.maxPerStageDescriptorUpdateAfterBindStorageImages < lim) lim = p12.maxPerStageDescriptorUpdateAfterBindStorageImages;
        if (p12.maxDescriptorSetUpdateAfterBindSampledImages < lim) lim = p12.maxDescriptorSetUpdateAfterBindSampledImages;
        if (p12.maxDescriptorSetUpdateAfterBindStorageImages < lim) lim = p12.maxDescriptorSetUpdateAfterBindStorageImages;
        const uint64_t res = p12.maxPerStageUpdateAfterBindResources;
        if (res > MD_VK_MAX_SAMPLERS && (res - MD_VK_MAX_SAMPLERS) / 2 < lim) lim = (uint32_t)((res - MD_VK_MAX_SAMPLERS) / 2);
        if (want > lim) {
            MD_LOG_INFO("md_gpu: %u texture slots requested; this device binds at most %u", want - 1, lim - 1);
            want = lim;
        }
        if (want < 2) { md_gpu_fail("md_gpu: the device has no bindless texture slots"); goto fail_device; }
        dev->tex_slot_count = want;
        dev->tex_free = (uint32_t*)md_alloc(alloc, (want - 1) * sizeof(uint32_t));
        if (!dev->tex_free) { md_gpu_fail("out of memory"); goto fail_device; }
        for (uint32_t i = 0; i < want - 1; ++i) dev->tex_free[i] = want - 1 - i;
        dev->tex_free_count = want - 1;
    }
    md_vk_init_device_types(dev);

    if (!md_vk_create_bindless(dev)) goto fail_device;
    if (!md_vk_create_dummies(dev))  goto fail_device;

    dev->default_compute  = md_vk_stream_create_internal(dev, MD_GPU_STREAM_COMPUTE,  "default compute",  true);
    dev->default_transfer = md_vk_stream_create_internal(dev, MD_GPU_STREAM_TRANSFER, "default transfer", true);
    if (!dev->default_compute || !dev->default_transfer) goto fail_device;
    if (dev->supports_graphics) {
        dev->default_graphics = md_vk_stream_create_internal(dev, MD_GPU_STREAM_GRAPHICS, "default graphics", true);
        if (!dev->default_graphics) goto fail_device;
    }

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
    if (dev->instance && dev->instance == md_vk_live_instance) md_vk_live_instance = VK_NULL_HANDLE;
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
    info->supports_graphics        = dev->supports_graphics;
    info->supports_present         = dev->supports_present;
    info->vendor_id                = dev->props.vendorID;
    info->device_id                = dev->props.deviceID;
    info->type                     = md_vk_device_type(dev->props.deviceType);
    info->adapter_index            = dev->adapter_index;
    md_snprintf(info->driver, sizeof(info->driver), "%s", dev->driver_desc);
    md_snprintf(info->name, sizeof(info->name), "%s", dev->props.deviceName);

    const VkSubgroupFeatureFlags arith = VK_SUBGROUP_FEATURE_BASIC_BIT | VK_SUBGROUP_FEATURE_ARITHMETIC_BIT |
                                         VK_SUBGROUP_FEATURE_BALLOT_BIT;
    info->max_shared_memory    = dev->props.limits.maxComputeSharedMemorySize;
    info->max_textures         = dev->tex_slot_count - 1;
    info->subgroup_arithmetic  = (dev->subgroup_ops & arith) == arith;
    info->subgroup_shuffle     = (dev->subgroup_ops & VK_SUBGROUP_FEATURE_SHUFFLE_BIT) != 0;
    info->float16              = dev->caps.float16;
    info->float64              = dev->caps.float64;
    info->int64                = dev->caps.shader_int64;
    info->atomic_int64         = dev->caps.int64_atomics;
    info->atomic_float_add     = dev->caps.atomic_float.shaderBufferFloat32AtomicAdd;
    info->atomic_float_min_max = dev->caps.atomic_float2.shaderBufferFloat32AtomicMinMax;
    info->timestamps           = dev->ts_bits_compute > 0;
    info->debug_memory         = dev->heap.debug;
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

    if (!md_vk_image_memory(dev, dev->dummy_image, &dev->dummy_mem, NULL)) return false;

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
    bindings[1].descriptorCount = dev->tex_slot_count;
    bindings[2].binding         = MD_VK_BINDING_STORAGE_IMAGE;
    bindings[2].descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings[2].descriptorCount = dev->tex_slot_count;

    for (uint32_t i = 0; i < MD_VK_BINDING_COUNT; ++i) {
        bindings[i].stageFlags = VK_SHADER_STAGE_ALL;
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
    sizes[1] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, dev->tex_slot_count};
    sizes[2] = (VkDescriptorPoolSize){VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, dev->tex_slot_count};

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

    /* Draws: the same bindless set, and the same 8-byte root pointer seen by
       both raster stages. */
    VkPushConstantRange rpcr = {VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, 8};
    plci.pPushConstantRanges = &rpcr;
    if (!md_vk_check(vkCreatePipelineLayout(dev->device, &plci, NULL, &dev->raster_layout), "vkCreatePipelineLayout (raster)")) return false;

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

/* A layout transition of every subresource of `image`, with the given scopes.
   The only image barrier md_gpu records: images live in GENERAL, and change
   layout only when they start (UNDEFINED, at creation and discards) and
   around a present. */
static void md_vk_cmd_image_barrier(VkCommandBuffer cmd, VkImage image, VkImageAspectFlags aspect,
                                    VkImageLayout from, VkImageLayout to,
                                    VkPipelineStageFlags2 src_stage, VkAccessFlags2 src_access,
                                    VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    VkImageMemoryBarrier2 b = {VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER_2};
    b.srcStageMask  = src_stage;
    b.srcAccessMask = src_access;
    b.dstStageMask  = dst_stage;
    b.dstAccessMask = dst_access;
    b.oldLayout     = from;
    b.newLayout     = to;
    b.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    b.image = image;
    b.subresourceRange.aspectMask = aspect;
    b.subresourceRange.levelCount = VK_REMAINING_MIP_LEVELS;
    b.subresourceRange.layerCount = VK_REMAINING_ARRAY_LAYERS;
    VkDependencyInfo di = {VK_STRUCTURE_TYPE_DEPENDENCY_INFO};
    di.imageMemoryBarrierCount = 1;
    di.pImageMemoryBarriers    = &b;
    vkCmdPipelineBarrier2(cmd, &di);
}

/* The same, into a stream: a barrier, not an operation. */
static bool md_vk_record_image_barrier(md_gpu_stream_t s, VkImage image, VkImageAspectFlags aspect,
                                       VkImageLayout from, VkImageLayout to,
                                       VkPipelineStageFlags2 src_stage, VkAccessFlags2 src_access,
                                       VkPipelineStageFlags2 dst_stage, VkAccessFlags2 dst_access) {
    if (!md_vk_stream_ensure_cmd(s)) return false;
    md_vk_cmd_image_barrier(s->open, image, aspect, from, to, src_stage, src_access, dst_stage, dst_access);
    s->has_work = true;
    return true;
}

static void md_vk_record_to_general(VkCommandBuffer cmd, void* user) {
    md_vk_cmd_image_barrier(cmd, *(VkImage*)user, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL,
                            VK_PIPELINE_STAGE_2_NONE, 0,
                            VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT);
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
    if (!s) { md_gpu_fail("out of memory"); return NULL; }
    memset(s, 0, sizeof(*s));
    s->device     = dev;
    s->kind       = kind;
    s->ordering   = MD_GPU_ORDER_IMPLICIT;
    s->next_value = 1;
    s->is_default = is_default;
    md_snprintf(s->label, sizeof(s->label), "%s", label ? label : "stream");

    md_mutex_lock(&dev->device_mutex);
    if (kind == MD_GPU_STREAM_TRANSFER) {
        s->family      = dev->transfer_family;
        s->can_compute = dev->transfer_can_compute;
        s->queue       = dev->transfer_queues[dev->next_transfer_queue % dev->transfer_queue_count];
        dev->next_transfer_queue++;
    } else if (kind == MD_GPU_STREAM_GRAPHICS) {
        s->family       = dev->graphics_family;
        s->can_compute  = true;
        s->can_graphics = true;
        s->queue        = dev->graphics_queues[dev->next_graphics_queue % dev->graphics_queue_count];
        dev->next_graphics_queue++;
    } else {
        s->family      = dev->compute_family;
        s->can_compute = true;
        s->queue       = dev->compute_queues[dev->next_compute_queue % dev->compute_queue_count];
        dev->next_compute_queue++;
    }
    md_mutex_unlock(&dev->device_mutex);

    MD_ARRAY_INIT(s->cmds, dev->alloc);
    MD_ARRAY_INIT(s->waits, dev->alloc);
    MD_ARRAY_INIT(s->arena.pages, dev->alloc);
    md_gpu_heap_temps_init(&dev->heap, &s->temps);

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

    /* The submissions' GPU times. A query is reset in the command buffer that
       writes it, which a transfer-only family cannot do: those are not timed.
       Without the pool the stream works as before, untimed. */
    {
        const uint32_t bits = kind == MD_GPU_STREAM_TRANSFER ? dev->ts_bits_transfer
                            : kind == MD_GPU_STREAM_GRAPHICS ? dev->ts_bits_graphics : dev->ts_bits_compute;
        if (bits > 0 && s->can_compute) {
            VkQueryPoolCreateInfo qpci = {VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
            qpci.queryType  = VK_QUERY_TYPE_TIMESTAMP;
            qpci.queryCount = 2 * MD_VK_TS_RING;
            if (vkCreateQueryPool(dev->device, &qpci, NULL, &s->ts_pool) == VK_SUCCESS) {
                s->ts_mask = bits >= 64 ? ~0ull : ((1ull << bits) - 1ull);
            } else {
                s->ts_pool = VK_NULL_HANDLE;
            }
        }
    }

    md_mutex_lock(&dev->device_mutex);
    const bool listed = MD_ARRAY_TRY_PUSH(dev->streams, s);
    md_mutex_unlock(&dev->device_mutex);
    if (!listed) {
        vkDestroyCommandPool(dev->device, s->cmd_pool, NULL);
        vkDestroySemaphore(dev->device, s->timeline, NULL);
        md_free(dev->alloc, s, sizeof(*s));
        md_gpu_fail("out of memory");
        return NULL;
    }
    return s;
}

md_gpu_stream_t md_gpu_stream_create(md_gpu_device_t dev, md_gpu_stream_kind_t kind, const char* label) {
    if (!dev) { md_gpu_fail("md_gpu_stream_create: null device"); return NULL; }
    if (kind != MD_GPU_STREAM_COMPUTE && kind != MD_GPU_STREAM_TRANSFER && kind != MD_GPU_STREAM_GRAPHICS) {
        md_gpu_fail("md_gpu_stream_create: invalid stream kind %d", (int)kind);
        return NULL;
    }
    if (kind == MD_GPU_STREAM_GRAPHICS && !dev->supports_graphics) {
        md_gpu_fail("md_gpu_stream_create: the device has no graphics support (md_gpu_device_info_t.supports_graphics)");
        return NULL;
    }
    return md_vk_stream_create_internal(dev, kind, label, false);
}

md_gpu_stream_t md_gpu_stream_default(md_gpu_device_t dev, md_gpu_stream_kind_t kind) {
    if (!dev) return NULL;
    switch (kind) {
    case MD_GPU_STREAM_TRANSFER: return dev->default_transfer;
    case MD_GPU_STREAM_GRAPHICS:
        if (!dev->default_graphics) md_gpu_fail("md_gpu_stream_default: the device has no graphics support");
        return dev->default_graphics;
    default:                     return dev->default_compute;
    }
}

md_gpu_device_t md_gpu_stream_device(md_gpu_stream_t s) { return s ? s->device : NULL; }

/* Acquire a command buffer, recycling any whose submission has completed. */
static bool md_vk_stream_ensure_cmd(md_gpu_stream_t s) {
    if (s->open) return true;
    md_gpu_device_t dev = s->device;

    uint64_t done = md_vk_stream_completed(s);
    md_vk_cmd_t* found = NULL;
    for (size_t i = 0; i < s->cmds.size; ++i) {
        md_vk_cmd_t* c = &s->cmds.data[i];
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
        const md_vk_cmd_t c = {cb, 0, false};
        if (!MD_ARRAY_TRY_PUSH(s->cmds, c)) { vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &cb); return md_gpu_fail("out of memory"); }
        found = &MD_ARRAY_BACK(s->cmds);
    }

    vkResetCommandBuffer(found->cmd, 0);
    VkCommandBufferBeginInfo bi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    if (!md_vk_check(vkBeginCommandBuffer(found->cmd, &bi), "vkBeginCommandBuffer")) return false;

    s->open               = found->cmd;
    s->has_work           = false;
    s->cb_head            = true;
    s->bound_kernel       = NULL;
    s->compute_set_bound  = false;
    s->graphics_set_bound = false;
    /* The GPU time of the submission this buffer becomes, which signals
       next_value: from here to its end (md_vk_stream_submit) */
    if (s->ts_pool) {
        const uint32_t q = 2u * (uint32_t)(s->next_value % MD_VK_TS_RING);
        vkCmdResetQueryPool(s->open, s->ts_pool, q, 2);
        vkCmdWriteTimestamp2(s->open, VK_PIPELINE_STAGE_2_TOP_OF_PIPE_BIT, s->ts_pool, q);
    }
    /* The first operation recorded here still needs a dependency on earlier
       work: batches submitted to one queue may overlap -- submission order
       alone creates no memory dependency. A barrier's first synchronisation
       scope covers everything earlier in submission order, which spans
       command buffers and batches, so md_vk_begin_op gives the buffer's first
       operation a barrier from every stage (cb_head). That also covers other
       streams' work on the same VkQueue, which matters when memory one of
       them freed is reused by another. (In EXPLICIT mode the caller's own
       barriers carry the same scope across submissions.) */
    return true;
}

/* The full-strength masks one queue family can express. A transfer-only
   family rejects compute and indirect stages, so they must be dropped there,
   not merely left unused. */
static void md_vk_stream_full_masks(md_gpu_stream_t s, VkPipelineStageFlags2* src_stage, VkAccessFlags2* src_access,
                                    VkPipelineStageFlags2* dst_stage, VkAccessFlags2* dst_access) {
    if (s->can_graphics) {
        /* The universal queue: everything, attachments included. A global
           memory barrier covers images as well as buffers, and every image
           stays in GENERAL, so this is all a pass needs on either side. */
        *src_stage  = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
        *src_access = VK_ACCESS_2_MEMORY_WRITE_BIT;
        *dst_stage  = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
        *dst_access = VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT;
        return;
    }
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

/* Operations that record work outside a render pass fail inside one. */
static bool md_vk_not_in_pass(md_gpu_stream_t s, const char* what) {
    if (!s->in_pass) return true;
    return md_gpu_fail("%s: stream '%s' is inside a render pass; call md_gpu_render_end first", what, s->label);
}

static void md_vk_stage_masks(md_gpu_stream_t s, md_gpu_stage_flags_t stages, bool producer,
                              VkPipelineStageFlags2* out_stage, VkAccessFlags2* out_access);

/* The command buffer an operation records into, with ordering applied. IMPLICIT
   ordering is one global memory barrier before the operation, from the stages
   the operations since the last such barrier wrote in, to the stages this one
   reads in (`consumers`). Chained, these order every operation after all the
   ones before it: each barrier's first scope includes the previous barrier's
   second. Nothing waits for more than it must -- a kernel after a kernel does
   not flush render targets. Returns VK_NULL_HANDLE on failure. */
static VkCommandBuffer md_vk_begin_op(md_gpu_stream_t s, md_gpu_stage_flags_t consumers) {
    if (s->upload_open) { md_gpu_fail("stream '%s' has an open upload; call md_gpu_upload_end first", s->label); return VK_NULL_HANDLE; }
    if (!md_vk_not_in_pass(s, "operation")) return VK_NULL_HANDLE;
    if (!md_vk_stream_ensure_cmd(s)) return VK_NULL_HANDLE;
    VkPipelineStageFlags2 ss = 0, ds = 0; VkAccessFlags2 sa = 0, da = 0;
    if (s->force_barrier) {
        md_vk_stream_full_masks(s, &ss, &sa, &ds, &da);
    } else if (s->ordering == MD_GPU_ORDER_IMPLICIT || s->order_next) {
        /* At the head of a command buffer: after everything earlier on the
           queue (see md_vk_stream_ensure_cmd). Within it: after what the
           operations since the last barrier wrote. */
        md_gpu_stage_flags_t from = s->cb_head ? MD_GPU_STAGE_ALL : s->written;
        if (from) {
            md_vk_stage_masks(s, from,      true,  &ss, &sa);
            md_vk_stage_masks(s, consumers, false, &ds, &da);
        }
    }
    if (ss && ds) {
        md_vk_emit_barrier(s, ss, sa, ds, da);
        s->written = 0;
    }
    s->force_barrier = false;
    s->order_next    = false;
    s->cb_head       = false;
    return s->open;
}

/* Mark that an operation was recorded, writing in `producers`. */
static void md_vk_end_op(md_gpu_stream_t s, md_gpu_stage_flags_t producers) {
    s->written |= producers;
    s->has_work = true;
}

/* Both edges of an EXPLICIT region are ordered: its first operation after
   everything before it, and, back in IMPLICIT, the next operation after
   everything in it -- which `written` already covers, as it is only cleared
   by a barrier md_gpu places. */
void md_gpu_stream_set_ordering(md_gpu_stream_t s, md_gpu_ordering_t ordering) {
    if (!s) return;
    if (ordering == MD_GPU_ORDER_EXPLICIT && s->ordering == MD_GPU_ORDER_IMPLICIT) s->order_next = true;
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
    /* Raster stages exist only on the universal queue; elsewhere they are
       dropped, like compute on a transfer-only queue. */
    if ((stages & MD_GPU_STAGE_VERTEX) && s->can_graphics) {
        if (producer) {
            st |= VK_PIPELINE_STAGE_2_VERTEX_SHADER_BIT;
            ac |= VK_ACCESS_2_SHADER_WRITE_BIT;
        } else {
            st |= VK_PIPELINE_STAGE_2_INDEX_INPUT_BIT | VK_PIPELINE_STAGE_2_VERTEX_SHADER_BIT;
            ac |= VK_ACCESS_2_INDEX_READ_BIT | VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT;
        }
    }
    if ((stages & MD_GPU_STAGE_FRAGMENT) && s->can_graphics) {
        st |= VK_PIPELINE_STAGE_2_FRAGMENT_SHADER_BIT;
        ac |= producer ? VK_ACCESS_2_SHADER_WRITE_BIT
                       : (VK_ACCESS_2_SHADER_READ_BIT | VK_ACCESS_2_SHADER_WRITE_BIT);
    }
    if ((stages & MD_GPU_STAGE_ATTACHMENT) && s->can_graphics) {
        st |= VK_PIPELINE_STAGE_2_EARLY_FRAGMENT_TESTS_BIT | VK_PIPELINE_STAGE_2_LATE_FRAGMENT_TESTS_BIT
            | VK_PIPELINE_STAGE_2_COLOR_ATTACHMENT_OUTPUT_BIT;
        ac |= VK_ACCESS_2_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;
        if (!producer) ac |= VK_ACCESS_2_COLOR_ATTACHMENT_READ_BIT | VK_ACCESS_2_DEPTH_STENCIL_ATTACHMENT_READ_BIT;
    }
    *out_stage  = st;
    *out_access = ac;
}

void md_gpu_barrier(md_gpu_stream_t s, md_gpu_stage_flags_t producers, md_gpu_stage_flags_t consumers) {
    if (!s) return;
    if (s->upload_open) { md_gpu_fail("md_gpu_barrier: stream '%s' has an open upload", s->label); return; }
    if (!md_vk_not_in_pass(s, "md_gpu_barrier")) return;
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

    uint64_t signal_value = s->next_value;
    if (s->ts_pool) {
        const uint32_t slot = (uint32_t)(signal_value % MD_VK_TS_RING);
        vkCmdWriteTimestamp2(s->open, VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, s->ts_pool, 2u * slot + 1u);
        s->ts_value[slot] = signal_value;
    }

    if (!md_vk_check(vkEndCommandBuffer(s->open), "vkEndCommandBuffer")) return false;

    VkCommandBufferSubmitInfo cbsi = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_SUBMIT_INFO};
    cbsi.commandBuffer = s->open;

    VkSemaphoreSubmitInfo  wait_stack[8];
    VkSemaphoreSubmitInfo* waits = wait_stack;
    const uint32_t tl_wait_count = (uint32_t)s->waits.size;
    const uint32_t wait_count    = tl_wait_count + s->bin_wait_count;
    if (wait_count > 8) {
        waits = (VkSemaphoreSubmitInfo*)md_alloc(dev->alloc, wait_count * sizeof(VkSemaphoreSubmitInfo));
        if (!waits) return md_gpu_fail("out of memory");
    }
    for (uint32_t i = 0; i < tl_wait_count; ++i) {
        md_vk_wait_t w = s->waits.data[i];
        waits[i] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        waits[i].semaphore = w.stream->timeline;
        waits[i].value     = w.value;
        waits[i].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }
    for (uint32_t i = 0; i < s->bin_wait_count; ++i) {
        VkSemaphoreSubmitInfo* w = &waits[tl_wait_count + i];
        *w = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        w->semaphore = s->bin_waits[i];
        w->stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }

    VkSemaphoreSubmitInfo signals[1 + 4];
    signals[0] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
    signals[0].semaphore = s->timeline;
    signals[0].value     = signal_value;
    signals[0].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    for (uint32_t i = 0; i < s->bin_signal_count; ++i) {
        signals[1 + i] = (VkSemaphoreSubmitInfo){VK_STRUCTURE_TYPE_SEMAPHORE_SUBMIT_INFO};
        signals[1 + i].semaphore = s->bin_signals[i];
        signals[1 + i].stageMask = VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT;
    }

    VkSubmitInfo2 si = {VK_STRUCTURE_TYPE_SUBMIT_INFO_2};
    si.waitSemaphoreInfoCount   = wait_count;
    si.pWaitSemaphoreInfos      = waits;
    si.commandBufferInfoCount   = 1;
    si.pCommandBufferInfos      = &cbsi;
    si.signalSemaphoreInfoCount = 1 + s->bin_signal_count;
    si.pSignalSemaphoreInfos    = signals;

    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueueSubmit2(s->queue, 1, &si, VK_NULL_HANDLE);
    md_mutex_unlock(&dev->queue_mutex);
    if (waits != wait_stack) md_free(dev->alloc, waits, wait_count * sizeof(VkSemaphoreSubmitInfo));
    if (!md_vk_check(r, "vkQueueSubmit2")) return false;

    /* Tag the command buffer and arena pages with the value that frees them. */
    for (size_t i = 0; i < s->cmds.size; ++i) {
        md_vk_cmd_t* c = &s->cmds.data[i];
        if (c->cmd == s->open) { c->value = signal_value; c->pending = true; break; }
    }
    md_vk_arena_retire(s, signal_value);

    s->submitted_value = signal_value;
    s->next_value      = signal_value + 1;
    s->waits.size     = 0;
    s->bin_wait_count  = 0;
    s->bin_signal_count = 0;
    s->open            = VK_NULL_HANDLE;
    s->has_work        = false;
    return true;
}

md_gpu_sync_t md_gpu_stream_record(md_gpu_stream_t s) {
    md_gpu_sync_t out = md_gpu_sync_none();
    if (!s) return out;
    if (s->upload_open) { md_gpu_fail("md_gpu_stream_record: stream '%s' has an open upload", s->label); return out; }
    if (!md_vk_not_in_pass(s, "md_gpu_stream_record")) return out;
    md_vk_stream_submit(s);
    if (s->submitted_value == 0) return out;
    out.stream = s;
    out.value  = s->submitted_value;
    return out;
}

void md_gpu_stream_wait(md_gpu_stream_t s, md_gpu_sync_t sync) {
    if (!s || !md_gpu_sync_is_valid(sync)) return;
    if (sync.stream == s) return;                    /* already ordered */
    if (sync.stream->device != s->device) { md_gpu_fail("md_gpu_stream_wait: sync from another device"); return; }
    if (!md_vk_not_in_pass(s, "md_gpu_stream_wait")) return;
    if (md_gpu_sync_is_complete(sync)) return;       /* nothing to wait for */

    /* Work already issued must not be retroactively delayed: close it first. */
    if (s->has_work) md_vk_stream_submit(s);

    /* Collapse duplicate waits on the same timeline, keeping the larger value. */
    for (size_t i = 0; i < s->waits.size; ++i) {
        md_vk_wait_t* w = &s->waits.data[i];
        if (w->stream == sync.stream) {
            if (sync.value > w->value) w->value = sync.value;
            return;
        }
    }
    const md_vk_wait_t w = {sync.stream, sync.value};
    if (!MD_ARRAY_TRY_PUSH(s->waits, w)) md_gpu_fail("out of memory queueing a stream wait");
}

void md_gpu_stream_flush(md_gpu_stream_t s) {
    if (!s) return;
    if (!md_vk_not_in_pass(s, "md_gpu_stream_flush")) return;
    md_vk_stream_submit(s);
}

void md_gpu_stream_sync(md_gpu_stream_t s) {
    if (!s) return;
    if (!md_vk_not_in_pass(s, "md_gpu_stream_sync")) return;
    md_vk_stream_submit(s);
    if (s->submitted_value > 0) {
        VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
        wi.semaphoreCount = 1;
        wi.pSemaphores    = &s->timeline;
        wi.pValues        = &s->submitted_value;
        md_vk_check(vkWaitSemaphores(s->device->device, &wi, UINT64_MAX), "vkWaitSemaphores");   /* a lost device shows here */
    }
}

bool md_gpu_sync_is_complete(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    return md_vk_stream_completed(sync.stream) >= sync.value;
}

bool md_gpu_sync_gpu_seconds(md_gpu_sync_t sync, double* seconds) {
    if (!md_gpu_sync_is_valid(sync)) return false;
    md_gpu_stream_t s = sync.stream;
    if (!s->ts_pool) return false;
    const uint32_t slot = (uint32_t)(sync.value % MD_VK_TS_RING);
    if (s->ts_value[slot] != sync.value) return false;      /* an older submission's: its slot was reused */
    if (md_vk_stream_completed(s) < sync.value) return false;
    uint64_t ts[2] = {0, 0};
    const VkResult r = vkGetQueryPoolResults(s->device->device, s->ts_pool, 2u * slot, 2, sizeof(ts), ts, sizeof(uint64_t), VK_QUERY_RESULT_64_BIT);
    if (r != VK_SUCCESS) return false;
    const uint64_t ticks = ((ts[1] & s->ts_mask) - (ts[0] & s->ts_mask)) & s->ts_mask;
    if (seconds) *seconds = (double)ticks * (double)s->device->ts_period * 1.0e-9;
    return true;
}

void md_gpu_sync_wait(md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return;
    VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
    wi.semaphoreCount = 1;
    wi.pSemaphores    = &sync.stream->timeline;
    wi.pValues        = &sync.value;
    md_vk_check(vkWaitSemaphores(sync.stream->device->device, &wi, UINT64_MAX), "vkWaitSemaphores");
}

/* Free a stream's own objects. The stream must be idle. */
static void md_vk_stream_free(md_gpu_stream_t s) {
    md_gpu_device_t dev = s->device;
    if (s->open) {
        vkEndCommandBuffer(s->open);
        s->open = VK_NULL_HANDLE;
    }
    for (size_t i = 0; i < s->cmds.size; ++i) {
        md_vk_cmd_t* c = &s->cmds.data[i];
        vkFreeCommandBuffers(dev->device, s->cmd_pool, 1, &c->cmd);
    }
    MD_ARRAY_FREE(s->cmds);
    MD_ARRAY_FREE(s->waits);
    md_vk_arena_free(dev, &s->arena);
    md_gpu_heap_temps_free(&dev->heap, &s->temps);
    if (s->ts_pool) vkDestroyQueryPool(dev->device, s->ts_pool, NULL);
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
    md_vk_abandon_pass(s);
    md_gpu_stream_sync(s);

    md_mutex_lock(&dev->device_mutex);
    MD_ARRAY_REMOVE_ITEM(dev->streams, s);
    md_vk_forget_stream_locked(dev, s);
    md_mutex_unlock(&dev->device_mutex);

    md_vk_stream_free(s);
}

/* =========================================================================
   9. Memory
   =========================================================================

   md_gpu_malloc, md_gpu_free and the temp scopes are md_gpu_heap.c's, shared
   with the Metal backend. What is Vulkan's is a chunk: one VkBuffer with
   VkDeviceMemory of its own, every buffer usage md_gpu needs and a device
   address -- guarded in debug memory, so that GPU-assisted validation sees an
   access past either end. Reusing freed memory early needs a barrier, which
   cannot be placed inside a render pass. */

static md_gpu_chunk_t* md_vk_chunk_create(void* user, md_gpu_mem_kind_t kind, uint64_t size) {
    md_gpu_device_t dev = (md_gpu_device_t)user;
    md_vk_chunk_t* c = (md_vk_chunk_t*)md_alloc(dev->alloc, sizeof(md_vk_chunk_t));
    if (!c) { md_gpu_fail("out of memory"); return NULL; }
    memset(c, 0, sizeof(*c));
    md_vk_raw_t raw;
    if (!md_vk_raw_create(dev, size, kind, dev->heap.debug ? MD_VK_RAW_GUARDED : 0, &raw)) {
        md_free(dev->alloc, c, sizeof(*c));
        return NULL;
    }
    c->buffer        = raw.buffer;
    c->memory        = raw.memory;
    c->memory_type   = raw.memory_type;
    c->memory_offset = raw.memory_offset;
    c->base.address  = raw.address;
    c->base.host     = raw.host;
    c->base.size     = size;
    c->base.kind     = kind;
    return &c->base;
}

static void md_vk_chunk_destroy(void* user, md_gpu_chunk_t* chunk) {
    md_gpu_device_t dev = (md_gpu_device_t)user;
    md_vk_chunk_t* c = (md_vk_chunk_t*)chunk;
    md_vk_destroy_raw_buffer(dev, c->buffer, c->memory);
    md_free(dev->alloc, c, sizeof(*c));
}

static bool md_vk_reuse_early(md_gpu_stream_t s) {
    if (s->in_pass) return false;
    s->force_barrier = true;
    return true;
}

static void md_vk_heap_init(md_gpu_device_t dev, uint64_t cache_limit, bool debug) {
    const md_gpu_heap_backend_t backend = {
        dev, md_vk_chunk_create, md_vk_chunk_destroy,
        md_vk_stream_position, md_vk_stream_completed, md_vk_reuse_early,
    };
    md_gpu_heap_init(&dev->heap, dev->alloc, &dev->device_mutex, &backend, cache_limit, debug);
}

md_gpu_mem_t md_gpu_malloc(md_gpu_stream_t s, md_gpu_mem_kind_t kind, size_t size) {
    const md_gpu_mem_t none = {0, NULL};
    if (!s) { md_gpu_fail("md_gpu_malloc: null stream"); return none; }
    return md_gpu_heap_malloc(&s->device->heap, s, kind, size);
}

void md_gpu_free(md_gpu_stream_t s, md_gpu_addr_t addr) {
    if (!addr) return;
    if (!s) { md_gpu_fail("md_gpu_free: a stream is required"); return; }
    md_gpu_heap_free(&s->device->heap, s, addr);
}

md_gpu_temp_t md_gpu_temp_begin(md_gpu_stream_t s) {
    const md_gpu_temp_t none = {NULL, 0};
    if (!s) { md_gpu_fail("md_gpu_temp_begin: null stream"); return none; }
    return md_gpu_heap_temp_begin(&s->temps, s);
}

md_gpu_mem_t md_gpu_temp_alloc(md_gpu_temp_t scope, md_gpu_mem_kind_t kind, size_t size) {
    const md_gpu_mem_t none = {0, NULL};
    md_gpu_stream_t s = scope.stream;
    if (!s) { md_gpu_fail("md_gpu_temp_alloc: null scope"); return none; }
    return md_gpu_heap_temp_alloc(&s->device->heap, &s->temps, scope, kind, size, s->label);
}

void md_gpu_temp_end(md_gpu_temp_t scope) {
    md_gpu_stream_t s = scope.stream;
    if (!s) { md_gpu_fail("md_gpu_temp_end: null scope"); return; }
    md_gpu_heap_temp_end(&s->device->heap, &s->temps, scope);
}

bool md_gpu_memory_stats(md_gpu_device_t dev, md_gpu_mem_kind_t kind, md_gpu_memory_stats_t* out) {
    if (!dev || !out || (unsigned)kind >= (unsigned)MD_GPU_MEM_KIND_COUNT) return md_gpu_fail("md_gpu_memory_stats: invalid argument");
    memset(out, 0, sizeof(*out));
    md_mutex_lock(&dev->device_mutex);
    md_gpu_heap_stats_locked(&dev->heap, kind, out);
    out->bytes_textures = kind == MD_GPU_MEM_DEVICE ? dev->texture_bytes : 0;
    md_mutex_unlock(&dev->device_mutex);
    return true;
}

/* ---- Copies ------------------------------------------------------------------ */

static bool md_vk_record_buffer_copy(md_gpu_stream_t s, VkBuffer src, uint64_t src_off,
                                     VkBuffer dst, uint64_t dst_off, uint64_t size) {
    VkCommandBuffer cmd = md_vk_begin_op(s, MD_GPU_STAGE_TRANSFER);
    if (!cmd) return false;
    VkBufferCopy region = {src_off, dst_off, size};
    vkCmdCopyBuffer(cmd, src, dst, 1, &region);
    md_vk_end_op(s, MD_GPU_STAGE_TRANSFER);
    return true;
}

bool md_gpu_copy(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_addr_t src, size_t size) {
    if (!s) return md_gpu_fail("md_gpu_copy: null stream");
    if (size == 0) return true;
    md_vk_span_t d, sp;
    if (!md_vk_resolve(s->device, dst, size, &d,  "md_gpu_copy (dst)")) return false;
    if (!md_vk_resolve(s->device, src, size, &sp, "md_gpu_copy (src)")) return false;
    return md_vk_record_buffer_copy(s, sp.buffer, sp.offset, d.buffer, d.offset, size);
}

/* Idle means a host write now cannot race anything this stream is ordered
   after: nothing recorded, nothing in flight, and every cross-stream wait
   queued by md_gpu_stream_wait already satisfied. */
static bool md_vk_stream_idle(md_gpu_stream_t s) {
    if (s->has_work || md_vk_stream_completed(s) < s->submitted_value) return false;
    for (size_t i = 0; i < s->waits.size; ++i) {
        md_vk_wait_t w = s->waits.data[i];
        if (md_vk_stream_completed(w.stream) < w.value) return false;
    }
    return true;
}

bool md_gpu_upload(md_gpu_stream_t s, md_gpu_addr_t dst, const void* src, size_t size) {
    if (!s) return md_gpu_fail("md_gpu_upload: null stream");
    if (size == 0) return true;
    if (!src) return md_gpu_fail("md_gpu_upload: null source");
    if (s->upload_open) return md_gpu_fail("md_gpu_upload: stream '%s' has an open upload", s->label);
    if (!md_vk_not_in_pass(s, "md_gpu_upload")) return false;
    md_vk_span_t d;
    if (!md_vk_resolve(s->device, dst, size, &d, "md_gpu_upload")) return false;

    /* Fast path: host-writable destination and nothing in flight here. */
    if (d.host && md_vk_stream_idle(s)) {
        memcpy(d.host, src, size);
        return true;
    }
    uint64_t addr; void* host; VkBuffer buf; uint64_t off;
    if (!md_vk_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    return md_vk_record_buffer_copy(s, buf, off, d.buffer, d.offset, size);
}

bool md_gpu_memset(md_gpu_stream_t s, md_gpu_addr_t dst, uint8_t value, size_t size) {
    if (!s) return md_gpu_fail("md_gpu_memset: null stream");
    if (!md_vk_not_in_pass(s, "md_gpu_memset")) return false;
    if (size == 0) return true;
    md_vk_span_t b;
    if (!md_vk_resolve(s->device, dst, size, &b, "md_gpu_memset")) return false;

    const uint64_t begin = b.offset;
    const uint64_t end   = b.offset + size;
    uint64_t abeg = MD_ALIGN_TO(begin, 4);
    uint64_t aend = end & ~3ull;

    if (aend > abeg) {
        VkCommandBuffer cmd = md_vk_begin_op(s, MD_GPU_STAGE_TRANSFER);
        if (!cmd) return false;
        vkCmdFillBuffer(cmd, b.buffer, abeg, aend - abeg, ((uint32_t)value) * 0x01010101u);
        md_vk_end_op(s, MD_GPU_STAGE_TRANSFER);
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
        if (!md_vk_record_buffer_copy(s, buf, soff, b.buffer, pieces[i][0], pieces[i][1])) return false;
    }
    return true;
}

void* md_gpu_upload_begin(md_gpu_stream_t s, md_gpu_addr_t dst, size_t size) {
    if (!s || !dst || size == 0) { md_gpu_fail("md_gpu_upload_begin: null argument"); return NULL; }
    if (s->upload_open) { md_gpu_fail("an upload is already open on stream '%s'", s->label); return NULL; }
    if (!md_vk_not_in_pass(s, "md_gpu_upload_begin")) return NULL;
    md_vk_span_t b;
    if (!md_vk_resolve(s->device, dst, size, &b, "md_gpu_upload_begin")) return NULL;

    /* Write straight into the destination when that cannot race the GPU. */
    if (b.host && md_vk_stream_idle(s)) {
        s->upload_open   = true;
        s->upload_direct = true;
        s->upload_dst    = dst;
        s->upload_size   = size;
        return b.host;
    }

    uint64_t addr; void* host;
    if (!md_vk_arena_alloc(s, size, &addr, &host, &s->upload_src, &s->upload_src_offset)) return NULL;
    s->upload_open   = true;
    s->upload_direct = false;
    s->upload_dst    = dst;
    s->upload_size   = size;
    return host;
}

bool md_gpu_upload_end(md_gpu_stream_t s) {
    if (!s || !s->upload_open) return md_gpu_fail("md_gpu_upload_end: no upload is open");
    s->upload_open = false;
    if (s->upload_direct) return true;

    md_vk_span_t d;
    if (!md_vk_resolve(s->device, s->upload_dst, s->upload_size, &d, "md_gpu_upload_end")) return false;
    return md_vk_record_buffer_copy(s, s->upload_src, s->upload_src_offset, d.buffer, d.offset, s->upload_size);
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
    if (dev->tex_free_count < dev->tex_slot_count - 1) dev->tex_free[dev->tex_free_count++] = slot;
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
    for (uint32_t i = 0; i < t->attach_view_count; ++i) vkDestroyImageView(dev->device, t->attach_views[i].view, NULL);
    if (t->attach_views) md_free(dev->alloc, t->attach_views, t->attach_view_cap * sizeof(md_vk_attach_view_t));
    if (!t->external) {
        if (t->image)  vkDestroyImage(dev->device, t->image, NULL);
        if (t->memory) vkFreeMemory(dev->device, t->memory, NULL);
    }
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

/* Validate and normalise a texture description and fill in the VkImageCreateInfo
   for it. Fails, by name, on anything the device cannot make. */
static bool md_vk_texture_image_info(md_gpu_device_t dev, const md_gpu_texture_desc_t* desc, const char* what,
                                     md_gpu_texture_desc_t* out_d, md_vk_format_info_t* out_fi, VkImageCreateInfo* ici) {
    md_gpu_texture_desc_t d = *desc;
    const char* label = d.label ? d.label : "texture";
    md_vk_format_info_t fi = md_vk_format_info(d.format);
    if (fi.format == VK_FORMAT_UNDEFINED) return md_gpu_fail("%s '%s': invalid format %d", what, label, (int)d.format);
    if (d.type != MD_GPU_TEX_2D && d.type != MD_GPU_TEX_2D_ARRAY && d.type != MD_GPU_TEX_3D) {
        return md_gpu_fail("%s '%s': invalid type %d (a zero-initialised desc has no type)", what, label, (int)d.type);
    }
    if (!(d.usage & (MD_GPU_TEX_STORAGE | MD_GPU_TEX_SAMPLED | MD_GPU_TEX_RENDER_TARGET))) {
        return md_gpu_fail("%s '%s': usage is empty", what, label);
    }
    if (d.width == 0 || d.height == 0) return md_gpu_fail("%s '%s': zero width or height", what, label);
    if ((d.usage & MD_GPU_TEX_RENDER_TARGET) && d.type == MD_GPU_TEX_3D) {
        return md_gpu_fail("%s '%s': RENDER_TARGET needs a 2D or 2D_ARRAY texture", what, label);
    }
    if ((d.usage & MD_GPU_TEX_RENDER_TARGET) && !dev->supports_graphics) {
        return md_gpu_fail("%s '%s': RENDER_TARGET usage on a device without graphics support", what, label);
    }
    if (d.type == MD_GPU_TEX_2D && d.depth_or_layers > 1) {
        return md_gpu_fail("%s '%s': a 2D texture has depth_or_layers %u; use MD_GPU_TEX_3D or MD_GPU_TEX_2D_ARRAY",
                          what, label, d.depth_or_layers);
    }
    if (d.depth_or_layers == 0) d.depth_or_layers = 1;
    if (d.mip_levels == 0) d.mip_levels = 1;

    *ici = (VkImageCreateInfo){VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    ici->imageType     = d.type == MD_GPU_TEX_3D ? VK_IMAGE_TYPE_3D : VK_IMAGE_TYPE_2D;
    ici->format        = fi.format;
    ici->extent.width  = d.width;
    ici->extent.height = d.height;
    ici->extent.depth  = d.type == MD_GPU_TEX_3D ? d.depth_or_layers : 1;
    ici->mipLevels     = d.mip_levels;
    ici->arrayLayers   = d.type == MD_GPU_TEX_2D_ARRAY ? d.depth_or_layers : 1;
    ici->samples       = VK_SAMPLE_COUNT_1_BIT;
    ici->tiling        = VK_IMAGE_TILING_OPTIMAL;
    ici->usage         = VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT;
    if (d.usage & MD_GPU_TEX_STORAGE) ici->usage |= VK_IMAGE_USAGE_STORAGE_BIT;
    if (d.usage & MD_GPU_TEX_SAMPLED) ici->usage |= VK_IMAGE_USAGE_SAMPLED_BIT;
    if (d.usage & MD_GPU_TEX_RENDER_TARGET) {
        ici->usage |= fi.depth ? VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT : VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    }
    md_vk_set_sharing(dev, &ici->sharingMode, &ici->queueFamilyIndexCount, &ici->pQueueFamilyIndices);
    ici->initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    VkImageFormatProperties ifp;
    if (vkGetPhysicalDeviceImageFormatProperties(dev->phys, ici->format, ici->imageType, ici->tiling, ici->usage, 0, &ifp) != VK_SUCCESS) {
        return md_gpu_fail("%s '%s': the device does not support %s %s with usage%s%s%s",
                          what, label, md_vk_tex_type_name(d.type), fi.name,
                          (d.usage & MD_GPU_TEX_STORAGE) ? " STORAGE" : "",
                          (d.usage & MD_GPU_TEX_SAMPLED) ? " SAMPLED" : "",
                          (d.usage & MD_GPU_TEX_RENDER_TARGET) ? " RENDER_TARGET" : "");
    }
    /* Storage access goes through descriptors without a format qualifier.
       Where the device-wide features are missing, it is up to the format:
       writes are essential; a format that can only be written still works for
       kernels that do not read it, so that is a warning, once per format. */
    if ((d.usage & MD_GPU_TEX_STORAGE) &&
        !(dev->caps.storage_read_without_format && dev->caps.storage_write_without_format)) {
        VkFormatProperties3 fp3 = {VK_STRUCTURE_TYPE_FORMAT_PROPERTIES_3};
        VkFormatProperties2 fp2 = {VK_STRUCTURE_TYPE_FORMAT_PROPERTIES_2, &fp3};
        vkGetPhysicalDeviceFormatProperties2(dev->phys, ici->format, &fp2);
        const VkFormatFeatureFlags2 ff = fp3.optimalTilingFeatures;
        if (!(ff & VK_FORMAT_FEATURE_2_STORAGE_WRITE_WITHOUT_FORMAT_BIT)) {
            return md_gpu_fail("%s '%s': the device cannot write %s storage images without a format qualifier", what, label, fi.name);
        }
        if (!(ff & VK_FORMAT_FEATURE_2_STORAGE_READ_WITHOUT_FORMAT_BIT) && (uint32_t)d.format < 64u &&
            !(dev->warned_storage_read & (1ull << (uint32_t)d.format))) {
            dev->warned_storage_read |= 1ull << (uint32_t)d.format;
            MD_LOG_INFO("md_gpu: this device cannot read %s storage images without a format qualifier; "
                        "kernels may write '%s' but reads of it through a storage handle are undefined", fi.name, label);
        }
    }
    if (d.width > ifp.maxExtent.width || d.height > ifp.maxExtent.height || ici->extent.depth > ifp.maxExtent.depth ||
        d.mip_levels > ifp.maxMipLevels || ici->arrayLayers > ifp.maxArrayLayers) {
        return md_gpu_fail("%s '%s': %ux%ux%u with %u mips and %u layers exceeds the device limits for %s",
                          what, label, d.width, d.height, ici->extent.depth, d.mip_levels, ici->arrayLayers, fi.name);
    }
    *out_d  = d;
    *out_fi = fi;
    return true;
}

/* Where a texture placed in md_gpu memory goes: its own alignment, raised to
   bufferImageGranularity so that it shares no granularity page with the
   linear resources around it, and its size rounded to whole pages. */
typedef struct md_vk_placement_t {
    uint64_t align;
    uint64_t total;     /* the texture's pages, plus room to align any start address:
                           what md_gpu_texture_placement_size reports */
} md_vk_placement_t;

static void md_vk_placement(md_gpu_device_t dev, const VkMemoryRequirements* req, md_vk_placement_t* out) {
    const uint64_t gran = dev->props.limits.bufferImageGranularity > 1 ? dev->props.limits.bufferImageGranularity : 1;
    out->align = req->alignment > gran ? req->alignment : gran;
    out->total = MD_ALIGN_TO(req->size, gran) + out->align - 1;
}

uint64_t md_gpu_texture_placement_size(md_gpu_device_t dev, const md_gpu_texture_desc_t* desc) {
    if (!dev || !desc) { md_gpu_fail("md_gpu_texture_placement_size: null argument"); return 0; }
    md_gpu_texture_desc_t d; md_vk_format_info_t fi; VkImageCreateInfo ici;
    if (!md_vk_texture_image_info(dev, desc, "md_gpu_texture_placement_size: texture", &d, &fi, &ici)) return 0;
    VkDeviceImageMemoryRequirements dimr = {VK_STRUCTURE_TYPE_DEVICE_IMAGE_MEMORY_REQUIREMENTS};
    dimr.pCreateInfo = &ici;
    VkMemoryRequirements2 mr = {VK_STRUCTURE_TYPE_MEMORY_REQUIREMENTS_2};
    vkGetDeviceImageMemoryRequirements(dev->device, &dimr, &mr);
    if (!(mr.memoryRequirements.memoryTypeBits & dev->device_buffer_types)) {
        md_gpu_fail("md_gpu_texture_placement_size: this device cannot place %s textures in MD_GPU_MEM_DEVICE memory", fi.name);
        return 0;
    }
    md_vk_placement_t pl;
    md_vk_placement(dev, &mr.memoryRequirements, &pl);
    return pl.total;
}

/* Record a layout transition of every subresource of `t` from UNDEFINED to
   GENERAL. A barrier, not an operation: it orders itself after everything
   before it and before everything after it, in either ordering mode. */
static bool md_vk_texture_init(md_gpu_stream_t s, md_gpu_texture_t t) {
    const VkImageAspectFlags aspect = t->fi.depth
        ? (VK_IMAGE_ASPECT_DEPTH_BIT | (t->desc.format == MD_GPU_FORMAT_D32_FLOAT_S8_UINT ? VK_IMAGE_ASPECT_STENCIL_BIT : 0))
        : VK_IMAGE_ASPECT_COLOR_BIT;
    return md_vk_record_image_barrier(s, t->image, aspect, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL,
                                      VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
                                      VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT);
}

/* Both texture creations: in memory of its own (`at` 0), or placed at `at`. */
static md_gpu_texture_t md_vk_texture_create(md_gpu_stream_t s, const md_gpu_texture_desc_t* desc, md_gpu_addr_t at, const char* what) {
    if (!s || !desc) { md_gpu_fail("%s: null argument", what); return NULL; }
    md_gpu_device_t dev = s->device;
    if (s->upload_open) { md_gpu_fail("%s: stream '%s' has an open upload", what, s->label); return NULL; }
    if (!md_vk_not_in_pass(s, what)) return NULL;

    md_gpu_texture_desc_t d; md_vk_format_info_t fi; VkImageCreateInfo ici;
    if (!md_vk_texture_image_info(dev, desc, "texture", &d, &fi, &ici)) return NULL;
    const char* label = d.label ? d.label : "texture";

    md_gpu_texture_t t = (md_gpu_texture_t)md_alloc(dev->alloc, sizeof(md_gpu_texture));
    if (!t) { md_gpu_fail("out of memory"); return NULL; }
    memset(t, 0, sizeof(*t));
    t->device = dev;
    t->fi     = fi;
    t->desc   = d;
    md_snprintf(t->label, sizeof(t->label), "%s", label);
    t->desc.label = t->label;

    if (!md_vk_check(vkCreateImage(dev->device, &ici, NULL, &t->image), "vkCreateImage")) goto fail;
    if (!at) {
        if (!md_vk_image_memory(dev, t->image, &t->memory, &t->bytes)) goto fail;
    } else {
        /* Placed: the texture's pages lie inside [at, at + placement size) of
           a live DEVICE allocation, in its chunk's VkDeviceMemory. */
        VkMemoryRequirements req;
        vkGetImageMemoryRequirements(dev->device, t->image, &req);
        md_vk_placement_t pl;
        md_vk_placement(dev, &req, &pl);
        md_vk_span_t sp;
        if (!md_vk_resolve(dev, at, pl.total, &sp, what)) goto fail;
        if (sp.chunk->base.kind != MD_GPU_MEM_DEVICE) {
            md_gpu_fail("%s '%s': textures are placed in MD_GPU_MEM_DEVICE memory", what, label);
            goto fail;
        }
        if (!(req.memoryTypeBits & (1u << sp.chunk->memory_type))) {
            md_gpu_fail("%s '%s': this device cannot place %s textures in MD_GPU_MEM_DEVICE memory", what, label, fi.name);
            goto fail;
        }
        const uint64_t off = MD_ALIGN_TO(sp.chunk->memory_offset + sp.offset, pl.align);
        if (!md_vk_check(vkBindImageMemory(dev->device, t->image, sp.chunk->memory, off), "vkBindImageMemory (placed)")) goto fail;
        /* The memory is the caller's: t->memory stays null, nothing frees it. */
    }

    if (d.usage & MD_GPU_TEX_SAMPLED) {
        t->sampled_view = md_vk_create_view(dev, t, 0, d.mip_levels);
        if (!t->sampled_view) goto fail;
    }
    if (d.usage & MD_GPU_TEX_STORAGE) {
        t->storage_views = (VkImageView*)md_alloc(dev->alloc, d.mip_levels * sizeof(VkImageView));
        t->storage_slots = (uint32_t*)md_alloc(dev->alloc, d.mip_levels * sizeof(uint32_t));
        if (!t->storage_views || !t->storage_slots) { md_gpu_fail("out of memory"); goto fail; }
        memset(t->storage_views, 0, d.mip_levels * sizeof(VkImageView));
        memset(t->storage_slots, 0, d.mip_levels * sizeof(uint32_t));
        for (uint32_t m = 0; m < d.mip_levels; ++m) {
            t->storage_views[m] = md_vk_create_view(dev, t, m, 1);
            if (!t->storage_views[m]) goto fail;
        }
    }

    /* The one layout transition this image ever gets (but for discards),
       recorded into the stream, so creation is stream-ordered like
       md_gpu_malloc and never waits for anything. */
    if (!md_vk_texture_init(s, t)) goto fail;

    md_mutex_lock(&dev->device_mutex);
    {
        uint32_t needed = ((d.usage & MD_GPU_TEX_SAMPLED) ? 1u : 0u) + ((d.usage & MD_GPU_TEX_STORAGE) ? d.mip_levels : 0u);
        if (dev->tex_free_count < needed) {
            /* The transition is already recorded; destroy through the retire
               path so the image outlives that command buffer. */
            md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
            md_mutex_unlock(&dev->device_mutex);
            md_gpu_fail("texture '%s': out of bindless heap slots (%u; md_gpu_device_desc_t.max_textures)", label, dev->tex_slot_count - 1);
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
        if (!MD_ARRAY_TRY_PUSH(dev->textures, t)) {
            md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
            md_mutex_unlock(&dev->device_mutex);
            md_gpu_fail("out of memory");
            return NULL;
        }
        dev->texture_bytes += t->bytes;
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

md_gpu_texture_t md_gpu_texture_create(md_gpu_stream_t s, const md_gpu_texture_desc_t* desc) {
    return md_vk_texture_create(s, desc, 0, "md_gpu_texture_create");
}

md_gpu_texture_t md_gpu_texture_create_at(md_gpu_stream_t s, const md_gpu_texture_desc_t* desc, md_gpu_addr_t memory) {
    if (!memory) { md_gpu_fail("md_gpu_texture_create_at: null memory address"); return NULL; }
    return md_vk_texture_create(s, desc, memory, "md_gpu_texture_create_at");
}

bool md_gpu_texture_discard(md_gpu_stream_t s, md_gpu_texture_t t) {
    if (!s || !t) return md_gpu_fail("md_gpu_texture_discard: null argument");
    if (t->device != s->device) return md_gpu_fail("md_gpu_texture_discard: texture '%s' belongs to another device", t->label);
    if (t->external) return md_gpu_fail("md_gpu_texture_discard: '%s' belongs to a surface", t->label);
    if (s->upload_open) return md_gpu_fail("md_gpu_texture_discard: stream '%s' has an open upload", s->label);
    if (!md_vk_not_in_pass(s, "md_gpu_texture_discard")) return false;
    return md_vk_texture_init(s, t);
}

void md_gpu_texture_destroy(md_gpu_texture_t t) {
    if (!t) return;
    if (t->external) { md_gpu_fail("md_gpu_texture_destroy: '%s' belongs to a surface and is not destroyed by the caller", t->label); return; }
    md_gpu_device_t dev = t->device;
    md_mutex_lock(&dev->device_mutex);
    MD_ARRAY_REMOVE_ITEM(dev->textures, t);
    dev->texture_bytes -= t->bytes;
    md_vk_retire_locked(dev, MD_VK_RETIRE_TEXTURE, t);
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
    if (!dev) { md_gpu_fail("md_gpu_sampler: null device"); return h; }
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
        md_gpu_fail("md_gpu_sampler: more than %u distinct samplers", MD_VK_MAX_SAMPLERS - 1);
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
    if (r->mip >= d->mip_levels) return md_gpu_fail("%s: mip %u out of range (texture '%s' has %u)", what, r->mip, t->label, d->mip_levels);

    uint32_t dim[3];
    dim[0] = d->width  >> r->mip; if (!dim[0]) dim[0] = 1;
    dim[1] = d->height >> r->mip; if (!dim[1]) dim[1] = 1;
    if (d->type == MD_GPU_TEX_3D) { dim[2] = d->depth_or_layers >> r->mip; if (!dim[2]) dim[2] = 1; }
    else                          { dim[2] = d->depth_or_layers; }

    uint32_t ext[3];
    for (int i = 0; i < 3; ++i) {
        if (r->offset[i] >= dim[i]) {
            return md_gpu_fail("%s: offset[%d] = %u is outside texture '%s' (extent %u at mip %u)",
                              what, i, r->offset[i], t->label, dim[i], r->mip);
        }
        ext[i] = r->extent[i] ? r->extent[i] : dim[i] - r->offset[i];
        if (ext[i] > dim[i] - r->offset[i]) {
            return md_gpu_fail("%s: region [%u, +%u) on axis %d overruns texture '%s' (extent %u at mip %u)",
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
        return md_gpu_fail("%s: buffer address must be %llu-byte aligned for %s", what, (unsigned long long)a, t->fi.name);
    }
    return true;
}

static bool md_vk_record_texture_copy(md_gpu_stream_t s, md_gpu_texture_t t, VkBuffer buffer,
                                      const VkBufferImageCopy* copy, bool to_texture) {
    VkCommandBuffer cmd = md_vk_begin_op(s, MD_GPU_STAGE_TRANSFER);
    if (!cmd) return false;
    if (to_texture) vkCmdCopyBufferToImage(cmd, buffer, t->image, VK_IMAGE_LAYOUT_GENERAL, 1, copy);
    else            vkCmdCopyImageToBuffer(cmd, t->image, VK_IMAGE_LAYOUT_GENERAL, buffer, 1, copy);
    md_vk_end_op(s, MD_GPU_STAGE_TRANSFER);
    return true;
}

bool md_gpu_copy_to_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, md_gpu_addr_t src) {
    if (!s || !t) return md_gpu_fail("md_gpu_copy_to_texture: null argument");
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_copy_to_texture")) return false;
    md_vk_span_t b;
    if (!md_vk_resolve(s->device, src, cr.bytes, &b, "md_gpu_copy_to_texture")) return false;
    if (!md_vk_check_buffer_offset(t, b.offset, "md_gpu_copy_to_texture")) return false;
    cr.copy.bufferOffset = b.offset;
    return md_vk_record_texture_copy(s, t, b.buffer, &cr.copy, true);
}

bool md_gpu_copy_from_texture(md_gpu_stream_t s, md_gpu_addr_t dst, md_gpu_texture_t t, const md_gpu_tex_region_t* region) {
    if (!s || !t) return md_gpu_fail("md_gpu_copy_from_texture: null argument");
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_copy_from_texture")) return false;
    md_vk_span_t b;
    if (!md_vk_resolve(s->device, dst, cr.bytes, &b, "md_gpu_copy_from_texture")) return false;
    if (!md_vk_check_buffer_offset(t, b.offset, "md_gpu_copy_from_texture")) return false;
    cr.copy.bufferOffset = b.offset;
    return md_vk_record_texture_copy(s, t, b.buffer, &cr.copy, false);
}

bool md_gpu_upload_texture(md_gpu_stream_t s, md_gpu_texture_t t, const md_gpu_tex_region_t* region, const void* src, size_t size) {
    if (!s || !t || !src) return md_gpu_fail("md_gpu_upload_texture: null argument");
    if (s->upload_open) return md_gpu_fail("md_gpu_upload_texture: stream '%s' has an open upload", s->label);
    if (!md_vk_not_in_pass(s, "md_gpu_upload_texture")) return false;
    md_vk_copy_region_t cr;
    if (!md_vk_resolve_region(t, region, &cr, "md_gpu_upload_texture")) return false;
    if ((uint64_t)size != cr.bytes) {
        return md_gpu_fail("md_gpu_upload_texture: region of texture '%s' is %llu bytes but %zu were given",
                          t->label, (unsigned long long)cr.bytes, size);
    }
    uint64_t addr; void* host; VkBuffer buf; uint64_t off;
    if (!md_vk_arena_alloc(s, size, &addr, &host, &buf, &off)) return false;
    memcpy(host, src, size);
    cr.copy.bufferOffset = off;
    return md_vk_record_texture_copy(s, t, buf, &cr.copy, true);
}

bool md_gpu_copy_texture(md_gpu_stream_t s, md_gpu_texture_t dst, const md_gpu_tex_region_t* dst_region,
                         md_gpu_texture_t src, const md_gpu_tex_region_t* src_region) {
    if (!s || !dst || !src) return md_gpu_fail("md_gpu_copy_texture: null argument");
    if (dst->device != s->device || src->device != s->device) return md_gpu_fail("md_gpu_copy_texture: texture from another device");
    if (dst->desc.format != src->desc.format) {
        return md_gpu_fail("md_gpu_copy_texture: formats differ ('%s' is %s, '%s' is %s); a converting copy is a draw or a kernel",
                          src->label, src->fi.name, dst->label, dst->fi.name);
    }
    if ((dst->desc.type == MD_GPU_TEX_3D) != (src->desc.type == MD_GPU_TEX_3D)) {
        return md_gpu_fail("md_gpu_copy_texture: cannot copy between a 3D texture and a 2D one ('%s' -> '%s')", src->label, dst->label);
    }
    md_vk_copy_region_t sr, dr;
    if (!md_vk_resolve_region(src, src_region, &sr, "md_gpu_copy_texture (src)")) return false;
    /* The destination takes the source's extent. On axis 2 that is layers
       for arrays (and plain 2D, one layer) and depth for 3D. */
    md_gpu_tex_region_t d;
    memset(&d, 0, sizeof(d));
    if (dst_region) d = *dst_region;
    d.extent[0] = sr.copy.imageExtent.width;
    d.extent[1] = sr.copy.imageExtent.height;
    d.extent[2] = src->desc.type == MD_GPU_TEX_3D ? sr.copy.imageExtent.depth : sr.copy.imageSubresource.layerCount;
    if (dst->desc.type == MD_GPU_TEX_2D && d.extent[2] != 1) {
        return md_gpu_fail("md_gpu_copy_texture: %u layers of '%s' do not fit 2D texture '%s'", d.extent[2], src->label, dst->label);
    }
    if (!md_vk_resolve_region(dst, &d, &dr, "md_gpu_copy_texture (dst)")) return false;
    if (!md_vk_not_in_pass(s, "md_gpu_copy_texture")) return false;

    VkImageCopy ic;
    memset(&ic, 0, sizeof(ic));
    ic.srcSubresource = sr.copy.imageSubresource;
    ic.srcOffset      = sr.copy.imageOffset;
    ic.dstSubresource = dr.copy.imageSubresource;
    ic.dstOffset      = dr.copy.imageOffset;
    ic.extent         = sr.copy.imageExtent;

    VkCommandBuffer cmd = md_vk_begin_op(s, MD_GPU_STAGE_TRANSFER);
    if (!cmd) return false;
    vkCmdCopyImage(cmd, src->image, VK_IMAGE_LAYOUT_GENERAL, dst->image, VK_IMAGE_LAYOUT_GENERAL, 1, &ic);
    md_vk_end_op(s, MD_GPU_STAGE_TRANSFER);
    return true;
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

/* The driver's statistics of a kernel's pipeline as "name value, ..." (its first executable). */
static void md_vk_kernel_stats(md_gpu_device_t dev, md_gpu_kernel_t k) {
    VkPipelineInfoKHR pi = {VK_STRUCTURE_TYPE_PIPELINE_INFO_KHR};
    pi.pipeline = k->pipeline;
    uint32_t nexe = 0;
    if (vkGetPipelineExecutablePropertiesKHR(dev->device, &pi, &nexe, NULL) != VK_SUCCESS || nexe == 0) return;
    VkPipelineExecutableInfoKHR ei = {VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_INFO_KHR};
    ei.pipeline = k->pipeline;
    ei.executableIndex = 0;
    uint32_t n = 0;
    if (vkGetPipelineExecutableStatisticsKHR(dev->device, &ei, &n, NULL) != VK_SUCCESS || n == 0) return;
    if (n > 32) n = 32;
    VkPipelineExecutableStatisticKHR st[32];
    for (uint32_t i = 0; i < n; ++i) st[i] = (VkPipelineExecutableStatisticKHR){VK_STRUCTURE_TYPE_PIPELINE_EXECUTABLE_STATISTIC_KHR};
    VkResult r = vkGetPipelineExecutableStatisticsKHR(dev->device, &ei, &n, st);
    if (r != VK_SUCCESS && r != VK_INCOMPLETE) return;
    size_t len = 0;
    for (uint32_t i = 0; i < n && len + 1 < sizeof(k->stats); ++i) {
        char v[32];
        switch (st[i].format) {
            case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_BOOL32_KHR:  md_snprintf(v, sizeof(v), "%s", st[i].value.b32 ? "yes" : "no"); break;
            case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_INT64_KHR:   md_snprintf(v, sizeof(v), "%lld", (long long)st[i].value.i64); break;
            case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_UINT64_KHR:  md_snprintf(v, sizeof(v), "%llu", (unsigned long long)st[i].value.u64); break;
            case VK_PIPELINE_EXECUTABLE_STATISTIC_FORMAT_FLOAT64_KHR: md_snprintf(v, sizeof(v), "%g", st[i].value.f64); break;
            default: v[0] = 0; break;
        }
        len += (size_t)md_snprintf(k->stats + len, sizeof(k->stats) - len, "%s%s %s", i ? ", " : "", st[i].name, v);
    }
    if (len >= sizeof(k->stats)) k->stats[sizeof(k->stats) - 1] = 0;
}

static bool md_vk_debug_labels(md_gpu_device_t dev);
static void md_vk_set_name(md_gpu_device_t dev, VkObjectType type, uint64_t handle, const char* name);

/* A constants block as specialization info: word i is constant_id i. */
typedef struct md_vk_spec_t {
    VkSpecializationInfo     info;
    VkSpecializationMapEntry entries[MD_GPU_MAX_CONSTANTS];
} md_vk_spec_t;

/* NULL when there are no constants, the info otherwise; fails on a malformed
   block. */
static bool md_vk_spec_init(md_vk_spec_t* sp, const void* data, uint32_t size, const char* label,
                            const VkSpecializationInfo** out) {
    *out = NULL;
    if (size == 0) return true;
    if (!data) return md_gpu_fail("'%s': %u bytes of constants but no data", label, size);
    if (size % 4 != 0 || size / 4 > MD_GPU_MAX_CONSTANTS) {
        return md_gpu_fail("'%s': constants_size %u is not a multiple of 4 up to %u", label, size, 4 * MD_GPU_MAX_CONSTANTS);
    }
    const uint32_t n = size / 4;
    for (uint32_t i = 0; i < n; ++i) {
        sp->entries[i].constantID = i;
        sp->entries[i].offset     = 4 * i;
        sp->entries[i].size       = 4;
    }
    sp->info.mapEntryCount = n;
    sp->info.pMapEntries   = sp->entries;
    sp->info.dataSize      = size;
    sp->info.pData         = data;
    *out = &sp->info;
    return true;
}

md_gpu_kernel_t md_gpu_kernel_create(md_gpu_device_t dev, const md_gpu_kernel_desc_t* desc) {
    if (!dev || !desc || !desc->code || desc->code_size == 0) {
        md_gpu_fail("md_gpu_kernel_create: missing code");
        return NULL;
    }
    const char* label = desc->label ? desc->label : "kernel";
    if (desc->code_size % 4 != 0) { md_gpu_fail("kernel '%s': SPIR-V size must be a multiple of 4", label); return NULL; }
    if (desc->group_size[0] == 0 || desc->group_size[1] == 0 || desc->group_size[2] == 0) {
        md_gpu_fail("kernel '%s': group_size is required ({%u,%u,%u} given); use the generated kernel descriptor",
                   label, desc->group_size[0], desc->group_size[1], desc->group_size[2]);
        return NULL;
    }
    uint32_t spv[3];
    if (md_vk_spirv_local_size((const uint32_t*)desc->code, desc->code_size / 4, spv) &&
        (spv[0] != desc->group_size[0] || spv[1] != desc->group_size[1] || spv[2] != desc->group_size[2])) {
        md_gpu_fail("kernel '%s': group_size {%u,%u,%u} does not match [numthreads(%u,%u,%u)] in the shader",
                   label, desc->group_size[0], desc->group_size[1], desc->group_size[2], spv[0], spv[1], spv[2]);
        return NULL;
    }
    const uint64_t threads = (uint64_t)desc->group_size[0] * desc->group_size[1] * desc->group_size[2];
    if (threads > dev->props.limits.maxComputeWorkGroupInvocations) {
        md_gpu_fail("kernel '%s': %llu threads per group exceeds the device limit of %u",
                   label, (unsigned long long)threads, dev->props.limits.maxComputeWorkGroupInvocations);
        return NULL;
    }
    md_vk_spec_t spec;
    const VkSpecializationInfo* spec_info = NULL;
    if (!md_vk_spec_init(&spec, desc->constants, desc->constants_size, label, &spec_info)) return NULL;

    md_gpu_kernel_t k = (md_gpu_kernel_t)md_alloc(dev->alloc, sizeof(md_gpu_kernel));
    if (!k) { md_gpu_fail("out of memory"); return NULL; }
    memset(k, 0, sizeof(*k));
    k->device       = dev;
    k->group_size[0] = desc->group_size[0];
    k->group_size[1] = desc->group_size[1];
    k->group_size[2] = desc->group_size[2];
    k->args_size    = desc->args_size;
    md_snprintf(k->label, sizeof(k->label), "%s", label);

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
    cpci.stage.pSpecializationInfo = spec_info;
    cpci.layout       = dev->pipeline_layout;
    if (dev->pipeline_stats) cpci.flags |= VK_PIPELINE_CREATE_CAPTURE_STATISTICS_BIT_KHR;
    if (!md_vk_check(vkCreateComputePipelines(dev->device, VK_NULL_HANDLE, 1, &cpci, NULL, &k->pipeline), "vkCreateComputePipelines")) {
        md_vk_kernel_free(dev, k);
        return NULL;
    }
    /* Some drivers (Intel Gen9 on Windows) report VK_SUCCESS but hand back no
       pipeline when their shader compiler rejects the module; binding that
       later faults inside the driver, so treat it as the failure it is. */
    if (k->pipeline == VK_NULL_HANDLE) {
        md_gpu_fail("kernel '%s': vkCreateComputePipelines returned VK_SUCCESS but no pipeline "
                   "(the driver could not compile the shader)", label);
        md_vk_kernel_free(dev, k);
        return NULL;
    }
    if (dev->pipeline_stats) md_vk_kernel_stats(dev, k);
    md_vk_set_name(dev, VK_OBJECT_TYPE_SHADER_MODULE, (uint64_t)k->module, k->label);
    md_vk_set_name(dev, VK_OBJECT_TYPE_PIPELINE, (uint64_t)k->pipeline, k->label);

    md_mutex_lock(&dev->device_mutex);
    const bool listed = MD_ARRAY_TRY_PUSH(dev->kernels, k);
    md_mutex_unlock(&dev->device_mutex);
    if (!listed) { md_vk_kernel_free(dev, k); md_gpu_fail("out of memory"); return NULL; }
    return k;
}

void md_gpu_kernel_destroy(md_gpu_kernel_t k) {
    if (!k) return;
    md_gpu_device_t dev = k->device;
    md_mutex_lock(&dev->device_mutex);
    MD_ARRAY_REMOVE_ITEM(dev->kernels, k);
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
    md_snprintf(info->stats, sizeof(info->stats), "%s", k->stats);
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

/* Each dispatch inside a command-buffer label with its kernel's name, so a profiler's timeline and
   per-dispatch tables say which kernel ran. */
static void md_vk_dispatch_label_begin(md_gpu_device_t dev, VkCommandBuffer cmd, md_gpu_kernel_t k) {
    if (!md_vk_debug_labels(dev)) return;
    VkDebugUtilsLabelEXT l = {VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT};
    l.pLabelName = k->label;
    vkCmdBeginDebugUtilsLabelEXT(cmd, &l);
}
static void md_vk_dispatch_label_end(md_gpu_device_t dev, VkCommandBuffer cmd) {
    if (md_vk_debug_labels(dev)) vkCmdEndDebugUtilsLabelEXT(cmd);
}

/* ---- Argument blocks ----------------------------------------------------------
   A launch or draw passes one 8-byte address in its push constant: of a copy
   of the caller's struct in the stream's arena, or (the _at forms) of a block
   the caller placed. */

static bool md_vk_args_copy(md_gpu_stream_t s, uint32_t expected, const void* args, size_t args_size,
                            uint64_t* out_addr, const char* what, const char* name) {
    *out_addr = 0;
    if (expected != 0 && args_size != expected) {
        return md_gpu_fail("%s: '%s' expects a %u-byte argument struct but %zu bytes were passed", what, name, expected, args_size);
    }
    if (args_size == 0) return true;
    if (!args) return md_gpu_fail("%s: '%s': null args with non-zero size", what, name);
    void* host;
    if (!md_vk_arena_alloc(s, args_size, out_addr, &host, NULL, NULL)) return false;
    memcpy(host, args, args_size);
    return true;
}

/* The alignment an argument struct of `size` bytes has in C, and so in an
   array of them: the largest power of two that divides its size, at most 16
   (md_gpu_*4, md_gpu_float4x4), at least 4. */
static uint64_t md_vk_args_align(uint32_t size) {
    uint64_t a = size ? (uint64_t)(size & (~size + 1u)) : 4u;
    return a < 4 ? 4 : (a > 16 ? 16 : a);
}

static bool md_vk_args_at(md_gpu_stream_t s, uint32_t expected, md_gpu_addr_t args, const char* what) {
    if (!args) return md_gpu_fail("%s: null argument address", what);
    const uint64_t align = md_vk_args_align(expected);
    if (args % align != 0) {
        return md_gpu_fail("%s: argument address 0x%llx is not %llu-byte aligned, as its %u-byte struct is",
                          what, (unsigned long long)args, (unsigned long long)align, expected);
    }
    return md_vk_resolve(s->device, args, expected ? expected : 1, NULL, what);
}

/* ---- Launches ------------------------------------------------------------------ */

static bool md_vk_launch_check(md_gpu_stream_t s, md_gpu_kernel_t k, const char* what) {
    if (k->device != s->device) return md_gpu_fail("%s: kernel '%s' belongs to a different device", what, k->label);
    if (s->kind == MD_GPU_STREAM_TRANSFER) {
        return md_gpu_fail("%s: kernel '%s' launched into '%s', a transfer stream; kernels run on compute streams", what, k->label, s->label);
    }
    return true;
}

static bool md_vk_grid_check(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid) {
    const VkPhysicalDeviceLimits* lim = &s->device->props.limits;
    if (grid.x > lim->maxComputeWorkGroupCount[0] || grid.y > lim->maxComputeWorkGroupCount[1] ||
        grid.z > lim->maxComputeWorkGroupCount[2]) {
        return md_gpu_fail("kernel '%s': grid {%u,%u,%u} exceeds the device limit {%u,%u,%u}", k->label,
                          grid.x, grid.y, grid.z, lim->maxComputeWorkGroupCount[0],
                          lim->maxComputeWorkGroupCount[1], lim->maxComputeWorkGroupCount[2]);
    }
    return true;
}

/* Record one dispatch: direct over `grid`, or indirect from `indirect`. */
static bool md_vk_dispatch(md_gpu_stream_t s, md_gpu_kernel_t k, uint64_t args_addr, md_gpu_grid_t grid,
                           const md_vk_span_t* indirect) {
    md_gpu_device_t dev = s->device;
    VkCommandBuffer cmd = md_vk_begin_op(s, indirect ? (MD_GPU_STAGE_COMPUTE | MD_GPU_STAGE_INDIRECT) : MD_GPU_STAGE_COMPUTE);
    if (!cmd) return false;
    /* Every kernel shares one pipeline layout, so the set stays bound across
       pipelines; both are bound once per command buffer and kept. */
    if (s->bound_kernel != k) {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, k->pipeline);
        s->bound_kernel = k;
    }
    if (!s->compute_set_bound) {
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, dev->pipeline_layout, 0, 1, &dev->desc_set, 0, NULL);
        s->compute_set_bound = true;
    }
    vkCmdPushConstants(cmd, dev->pipeline_layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, 8, &args_addr);
    md_vk_dispatch_label_begin(dev, cmd, k);
    if (indirect) vkCmdDispatchIndirect(cmd, indirect->buffer, indirect->offset);
    else          vkCmdDispatch(cmd, grid.x, grid.y, grid.z);
    md_vk_dispatch_label_end(dev, cmd);
    md_vk_end_op(s, MD_GPU_STAGE_COMPUTE);
    return true;
}

static bool md_vk_resolve_grid(md_gpu_stream_t s, md_gpu_addr_t grid, md_vk_span_t* out, const char* what) {
    if (!md_vk_resolve(s->device, grid, 3 * sizeof(uint32_t), out, what)) return false;
    if (out->offset % 4 != 0) return md_gpu_fail("%s: grid address must be 4-byte aligned", what);
    return true;
}

bool md_gpu_launch(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid, const void* args, size_t args_size) {
    const char* what = "md_gpu_launch";
    if (!s || !k) return md_gpu_fail("%s: null stream or kernel", what);
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;   /* empty launch */
    if (!md_vk_launch_check(s, k, what) || !md_vk_grid_check(s, k, grid)) return false;
    uint64_t a;
    if (!md_vk_args_copy(s, k->args_size, args, args_size, &a, what, k->label)) return false;
    return md_vk_dispatch(s, k, a, grid, NULL);
}

bool md_gpu_launch_at(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_grid_t grid, md_gpu_addr_t args) {
    const char* what = "md_gpu_launch_at";
    if (!s || !k) return md_gpu_fail("%s: null stream or kernel", what);
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return true;
    if (!md_vk_launch_check(s, k, what) || !md_vk_grid_check(s, k, grid)) return false;
    if (!md_vk_args_at(s, k->args_size, args, what)) return false;
    return md_vk_dispatch(s, k, args, grid, NULL);
}

bool md_gpu_launch_indirect(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_addr_t grid, const void* args, size_t args_size) {
    const char* what = "md_gpu_launch_indirect";
    if (!s || !k || !grid) return md_gpu_fail("%s: null argument", what);
    if (!md_vk_launch_check(s, k, what)) return false;
    md_vk_span_t b;
    if (!md_vk_resolve_grid(s, grid, &b, what)) return false;
    uint64_t a;
    if (!md_vk_args_copy(s, k->args_size, args, args_size, &a, what, k->label)) return false;
    return md_vk_dispatch(s, k, a, md_gpu_grid(0, 0, 0), &b);
}

bool md_gpu_launch_indirect_at(md_gpu_stream_t s, md_gpu_kernel_t k, md_gpu_addr_t grid, md_gpu_addr_t args) {
    const char* what = "md_gpu_launch_indirect_at";
    if (!s || !k || !grid) return md_gpu_fail("%s: null argument", what);
    if (!md_vk_launch_check(s, k, what)) return false;
    md_vk_span_t b;
    if (!md_vk_resolve_grid(s, grid, &b, what)) return false;
    if (!md_vk_args_at(s, k->args_size, args, what)) return false;
    return md_vk_dispatch(s, k, args, md_gpu_grid(0, 0, 0), &b);
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
    MD_ARRAY_REMOVE_ITEM(dev->kernels, dev->make_grid_kernel);
    md_mutex_unlock(&dev->device_mutex);
    return true;
}

bool md_gpu_make_grid(md_gpu_stream_t s, md_gpu_addr_t out_grid, md_gpu_addr_t count, md_gpu_kernel_t k) {
    if (!s || !out_grid || !count || !k) return md_gpu_fail("md_gpu_make_grid: null argument");
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
   12. Rendering
   =========================================================================

   Dynamic rendering (core 1.3) with every attachment in GENERAL, so a pass
   needs no layout transitions and no render-pass objects. The pass is
   ordered like any other operation: IMPLICIT mode puts one global memory
   barrier before vkCmdBeginRendering, and md_vk_end_op at render_end arms
   the one after it. Inside the pass nothing but draws and dynamic state is
   recorded, because a barrier inside dynamic rendering is only legal as a
   declared self-dependency.

   Pipelines bake what both APIs bake (shaders, topology, formats, blend) and
   leave dynamic what core 1.3 and Metal both make dynamic. Clip space is
   flipped to +Y up with a negative viewport height (core since 1.1); that
   also flips the winding seen in framebuffer space, which is why
   counter-clockwise in clip space maps to VK_FRONT_FACE_COUNTER_CLOCKWISE. */

static bool md_vk_format_is_uint(md_gpu_format_t f) {
    return f == MD_GPU_FORMAT_R32_UINT || f == MD_GPU_FORMAT_RG32_UINT || f == MD_GPU_FORMAT_RGBA32_UINT;
}

static bool md_vk_debug_labels(md_gpu_device_t dev) {
    return dev->debug_utils && vkCmdBeginDebugUtilsLabelEXT && vkCmdEndDebugUtilsLabelEXT;
}

/* The label of a kernel or pipeline as its Vulkan object name, for GPU tools. */
static void md_vk_set_name(md_gpu_device_t dev, VkObjectType type, uint64_t handle, const char* name) {
    if (!dev->debug_utils || !vkSetDebugUtilsObjectNameEXT || !handle || !name || !name[0]) return;
    VkDebugUtilsObjectNameInfoEXT ni = {VK_STRUCTURE_TYPE_DEBUG_UTILS_OBJECT_NAME_INFO_EXT};
    ni.objectType   = type;
    ni.objectHandle = handle;
    ni.pObjectName  = name;
    vkSetDebugUtilsObjectNameEXT(dev->device, &ni);
}

/* ---- Pipelines ------------------------------------------------------------ */

static VkBlendFactor md_vk_blend_factor(md_gpu_blend_factor_t f) {
    switch (f) {
    case MD_GPU_BLEND_ZERO:                 return VK_BLEND_FACTOR_ZERO;
    case MD_GPU_BLEND_ONE:                  return VK_BLEND_FACTOR_ONE;
    case MD_GPU_BLEND_SRC_COLOR:            return VK_BLEND_FACTOR_SRC_COLOR;
    case MD_GPU_BLEND_ONE_MINUS_SRC_COLOR:  return VK_BLEND_FACTOR_ONE_MINUS_SRC_COLOR;
    case MD_GPU_BLEND_SRC_ALPHA:            return VK_BLEND_FACTOR_SRC_ALPHA;
    case MD_GPU_BLEND_ONE_MINUS_SRC_ALPHA:  return VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
    case MD_GPU_BLEND_DST_COLOR:            return VK_BLEND_FACTOR_DST_COLOR;
    case MD_GPU_BLEND_ONE_MINUS_DST_COLOR:  return VK_BLEND_FACTOR_ONE_MINUS_DST_COLOR;
    case MD_GPU_BLEND_DST_ALPHA:            return VK_BLEND_FACTOR_DST_ALPHA;
    case MD_GPU_BLEND_ONE_MINUS_DST_ALPHA:  return VK_BLEND_FACTOR_ONE_MINUS_DST_ALPHA;
    case MD_GPU_BLEND_CONSTANT:             return VK_BLEND_FACTOR_CONSTANT_COLOR;
    case MD_GPU_BLEND_ONE_MINUS_CONSTANT:   return VK_BLEND_FACTOR_ONE_MINUS_CONSTANT_COLOR;
    case MD_GPU_BLEND_SRC_ALPHA_SATURATE:   return VK_BLEND_FACTOR_SRC_ALPHA_SATURATE;
    default:                                return VK_BLEND_FACTOR_MAX_ENUM;
    }
}

static VkBlendOp md_vk_blend_op(md_gpu_blend_op_t op) {
    switch (op) {
    case MD_GPU_BLEND_OP_ADD:              return VK_BLEND_OP_ADD;
    case MD_GPU_BLEND_OP_SUBTRACT:         return VK_BLEND_OP_SUBTRACT;
    case MD_GPU_BLEND_OP_REVERSE_SUBTRACT: return VK_BLEND_OP_REVERSE_SUBTRACT;
    case MD_GPU_BLEND_OP_MIN:              return VK_BLEND_OP_MIN;
    case MD_GPU_BLEND_OP_MAX:              return VK_BLEND_OP_MAX;
    default:                               return VK_BLEND_OP_MAX_ENUM;
    }
}

static VkPrimitiveTopology md_vk_topology(md_gpu_topology_t t) {
    switch (t) {
    case MD_GPU_TOPOLOGY_TRIANGLES:      return VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;
    case MD_GPU_TOPOLOGY_TRIANGLE_STRIP: return VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP;
    case MD_GPU_TOPOLOGY_LINES:          return VK_PRIMITIVE_TOPOLOGY_LINE_LIST;
    case MD_GPU_TOPOLOGY_LINE_STRIP:     return VK_PRIMITIVE_TOPOLOGY_LINE_STRIP;
    case MD_GPU_TOPOLOGY_POINTS:         return VK_PRIMITIVE_TOPOLOGY_POINT_LIST;
    default:                             return VK_PRIMITIVE_TOPOLOGY_MAX_ENUM;
    }
}

static void md_vk_pipeline_free(md_gpu_device_t dev, md_gpu_pipeline_t p) {
    if (p->pipeline) vkDestroyPipeline(dev->device, p->pipeline, NULL);
    md_free(dev->alloc, p, sizeof(*p));
}

static bool md_vk_check_shader(const md_gpu_shader_t* sh, const char* stage, const char* label) {
    if (sh->code_size % 4 != 0) return md_gpu_fail("pipeline '%s': %s shader SPIR-V size must be a multiple of 4", label, stage);
    const uint32_t* w = (const uint32_t*)sh->code;
    if (sh->code_size < 20 || w[0] != 0x07230203u) return md_gpu_fail("pipeline '%s': %s shader is not SPIR-V", label, stage);
    return true;
}

md_gpu_pipeline_t md_gpu_pipeline_create(md_gpu_device_t dev, const md_gpu_pipeline_desc_t* desc) {
    if (!dev || !desc) { md_gpu_fail("md_gpu_pipeline_create: null argument"); return NULL; }
    const char* label = desc->label ? desc->label : "pipeline";
    if (!dev->supports_graphics) { md_gpu_fail("pipeline '%s': the device has no graphics support", label); return NULL; }
    if (!desc->vertex.code || desc->vertex.code_size == 0) { md_gpu_fail("pipeline '%s': no vertex shader", label); return NULL; }
    if (!md_vk_check_shader(&desc->vertex, "vertex", label)) return NULL;
    const bool has_fs = desc->fragment.code != NULL;
    if (has_fs && !md_vk_check_shader(&desc->fragment, "fragment", label)) return NULL;
    if (desc->color_count > MD_GPU_MAX_COLOR_TARGETS) {
        md_gpu_fail("pipeline '%s': %u colour targets, at most %u", label, desc->color_count, MD_GPU_MAX_COLOR_TARGETS);
        return NULL;
    }
    if (!has_fs && desc->color_count > 0) {
        md_gpu_fail("pipeline '%s': colour targets without a fragment shader", label);
        return NULL;
    }
    const VkPrimitiveTopology topo = md_vk_topology(desc->topology);
    if (topo == VK_PRIMITIVE_TOPOLOGY_MAX_ENUM) { md_gpu_fail("pipeline '%s': invalid topology %d", label, (int)desc->topology); return NULL; }
    const uint32_t vs_args = desc->vertex.args_size;
    const uint32_t fs_args = has_fs ? desc->fragment.args_size : 0;
    if (vs_args && fs_args && vs_args != fs_args) {
        md_gpu_fail("pipeline '%s': the vertex shader takes a %u-byte argument struct and the fragment shader %u bytes; "
                   "both stages read the same one", label, vs_args, fs_args);
        return NULL;
    }

    VkFormat color_formats[MD_GPU_MAX_COLOR_TARGETS];
    VkPipelineColorBlendAttachmentState blend[MD_GPU_MAX_COLOR_TARGETS];
    memset(blend, 0, sizeof(blend));
    for (uint32_t i = 0; i < desc->color_count; ++i) {
        const md_gpu_color_target_t* ct = &desc->color[i];
        md_vk_format_info_t fi = md_vk_format_info(ct->format);
        if (fi.format == VK_FORMAT_UNDEFINED || fi.depth) {
            md_gpu_fail("pipeline '%s': colour target %u has %s, which is not a colour format", label, i,
                       fi.format == VK_FORMAT_UNDEFINED ? "no format" : fi.name);
            return NULL;
        }
        VkFormatProperties fp;
        vkGetPhysicalDeviceFormatProperties(dev->phys, fi.format, &fp);
        if (!(fp.optimalTilingFeatures & VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BIT)) {
            md_gpu_fail("pipeline '%s': the device cannot render to %s (colour target %u)", label, fi.name, i);
            return NULL;
        }
        const md_gpu_blend_t* b = &ct->blend;
        if (b->enable) {
            if (md_vk_format_is_uint(ct->format)) {
                md_gpu_fail("pipeline '%s': colour target %u is %s, an integer format, and cannot blend", label, i, fi.name);
                return NULL;
            }
            if (!(fp.optimalTilingFeatures & VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BLEND_BIT)) {
                md_gpu_fail("pipeline '%s': the device cannot blend %s (colour target %u)", label, fi.name, i);
                return NULL;
            }
            blend[i].blendEnable         = VK_TRUE;
            blend[i].srcColorBlendFactor = md_vk_blend_factor(b->src_color);
            blend[i].dstColorBlendFactor = md_vk_blend_factor(b->dst_color);
            blend[i].colorBlendOp        = md_vk_blend_op(b->color_op);
            blend[i].srcAlphaBlendFactor = md_vk_blend_factor(b->src_alpha);
            blend[i].dstAlphaBlendFactor = md_vk_blend_factor(b->dst_alpha);
            blend[i].alphaBlendOp        = md_vk_blend_op(b->alpha_op);
            if (blend[i].srcColorBlendFactor == VK_BLEND_FACTOR_MAX_ENUM || blend[i].dstColorBlendFactor == VK_BLEND_FACTOR_MAX_ENUM ||
                blend[i].srcAlphaBlendFactor == VK_BLEND_FACTOR_MAX_ENUM || blend[i].dstAlphaBlendFactor == VK_BLEND_FACTOR_MAX_ENUM ||
                blend[i].colorBlendOp == VK_BLEND_OP_MAX_ENUM || blend[i].alphaBlendOp == VK_BLEND_OP_MAX_ENUM) {
                md_gpu_fail("pipeline '%s': invalid blend factor or op on colour target %u", label, i);
                return NULL;
            }
        }
        blend[i].colorWriteMask = (VkColorComponentFlags)(MD_GPU_COLOR_ALL & ~ct->write_disable);
        color_formats[i] = fi.format;
    }

    VkFormat depth_format = VK_FORMAT_UNDEFINED;
    if (desc->depth_format != MD_GPU_FORMAT_INVALID) {
        md_vk_format_info_t fi = md_vk_format_info(desc->depth_format);
        if (!fi.depth) {
            md_gpu_fail("pipeline '%s': depth_format %s is not a depth format", label,
                       fi.format == VK_FORMAT_UNDEFINED ? "(invalid)" : fi.name);
            return NULL;
        }
        VkFormatProperties fp;
        vkGetPhysicalDeviceFormatProperties(dev->phys, fi.format, &fp);
        if (!(fp.optimalTilingFeatures & VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT)) {
            md_gpu_fail("pipeline '%s': the device cannot use %s as a depth attachment", label, fi.name);
            return NULL;
        }
        depth_format = fi.format;
    }

    md_vk_spec_t spec;
    const VkSpecializationInfo* spec_info = NULL;
    if (!md_vk_spec_init(&spec, desc->constants, desc->constants_size, label, &spec_info)) return NULL;

    md_gpu_pipeline_t p = (md_gpu_pipeline_t)md_alloc(dev->alloc, sizeof(md_gpu_pipeline));
    if (!p) { md_gpu_fail("out of memory"); return NULL; }
    memset(p, 0, sizeof(*p));
    p->device      = dev;
    p->color_count = desc->color_count;
    for (uint32_t i = 0; i < desc->color_count; ++i) p->color[i] = desc->color[i].format;
    p->depth       = desc->depth_format;
    p->args_size   = vs_args ? vs_args : fs_args;
    md_snprintf(p->label, sizeof(p->label), "%s", label);

    VkShaderModule modules[2] = {VK_NULL_HANDLE, VK_NULL_HANDLE};
    VkPipelineShaderStageCreateInfo stages[2];
    memset(stages, 0, sizeof(stages));
    const md_gpu_shader_t* sh[2] = {&desc->vertex, &desc->fragment};
    const VkShaderStageFlagBits bits[2] = {VK_SHADER_STAGE_VERTEX_BIT, VK_SHADER_STAGE_FRAGMENT_BIT};
    const uint32_t stage_count = has_fs ? 2u : 1u;
    bool ok = true;
    for (uint32_t i = 0; i < stage_count && ok; ++i) {
        VkShaderModuleCreateInfo smci = {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        smci.codeSize = sh[i]->code_size;
        smci.pCode    = (const uint32_t*)sh[i]->code;
        ok = md_vk_check(vkCreateShaderModule(dev->device, &smci, NULL, &modules[i]), "vkCreateShaderModule");
        stages[i].sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        stages[i].stage  = bits[i];
        stages[i].module = modules[i];
        /* A single-entry SPIR-V module from compile_gpu_shaders names its
           entry "main"; the generated descriptors leave entry_point NULL. */
        stages[i].pName  = sh[i]->entry_point ? sh[i]->entry_point : "main";
        stages[i].pSpecializationInfo = spec_info;
    }

    if (ok) {
        VkPipelineRenderingCreateInfo prci = {VK_STRUCTURE_TYPE_PIPELINE_RENDERING_CREATE_INFO};
        prci.colorAttachmentCount    = desc->color_count;
        prci.pColorAttachmentFormats = color_formats;
        prci.depthAttachmentFormat   = depth_format;
        prci.stencilAttachmentFormat = VK_FORMAT_UNDEFINED;

        VkPipelineVertexInputStateCreateInfo vi = {VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};

        VkPipelineInputAssemblyStateCreateInfo ia = {VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        ia.topology               = topo;
        ia.primitiveRestartEnable = (topo == VK_PRIMITIVE_TOPOLOGY_TRIANGLE_STRIP || topo == VK_PRIMITIVE_TOPOLOGY_LINE_STRIP)
                                  ? VK_TRUE : VK_FALSE;

        VkPipelineViewportStateCreateInfo vps = {VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
        vps.viewportCount = 1;
        vps.scissorCount  = 1;

        VkPipelineRasterizationStateCreateInfo rs = {VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        rs.polygonMode = VK_POLYGON_MODE_FILL;
        rs.lineWidth   = 1.0f;

        VkPipelineMultisampleStateCreateInfo ms = {VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
        ms.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

        VkPipelineDepthStencilStateCreateInfo ds = {VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};

        VkPipelineColorBlendStateCreateInfo cb = {VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
        cb.attachmentCount = desc->color_count;
        cb.pAttachments    = blend;

        static const VkDynamicState dyn[] = {
            VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR,
            VK_DYNAMIC_STATE_CULL_MODE, VK_DYNAMIC_STATE_FRONT_FACE,
            VK_DYNAMIC_STATE_DEPTH_TEST_ENABLE, VK_DYNAMIC_STATE_DEPTH_WRITE_ENABLE, VK_DYNAMIC_STATE_DEPTH_COMPARE_OP,
            VK_DYNAMIC_STATE_DEPTH_BIAS_ENABLE, VK_DYNAMIC_STATE_DEPTH_BIAS,
            VK_DYNAMIC_STATE_BLEND_CONSTANTS,
        };
        VkPipelineDynamicStateCreateInfo dy = {VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
        dy.dynamicStateCount = (uint32_t)(sizeof(dyn) / sizeof(dyn[0]));
        dy.pDynamicStates    = dyn;

        VkGraphicsPipelineCreateInfo gpci = {VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        gpci.pNext               = &prci;
        gpci.stageCount          = stage_count;
        gpci.pStages             = stages;
        gpci.pVertexInputState   = &vi;
        gpci.pInputAssemblyState = &ia;
        gpci.pViewportState      = &vps;
        gpci.pRasterizationState = &rs;
        gpci.pMultisampleState   = &ms;
        gpci.pDepthStencilState  = &ds;
        gpci.pColorBlendState    = &cb;
        gpci.pDynamicState       = &dy;
        gpci.layout              = dev->raster_layout;
        ok = md_vk_check(vkCreateGraphicsPipelines(dev->device, VK_NULL_HANDLE, 1, &gpci, NULL, &p->pipeline),
                         "vkCreateGraphicsPipelines");
    }
    for (uint32_t i = 0; i < 2; ++i) if (modules[i]) vkDestroyShaderModule(dev->device, modules[i], NULL);
    if (!ok) { md_vk_pipeline_free(dev, p); return NULL; }
    md_vk_set_name(dev, VK_OBJECT_TYPE_PIPELINE, (uint64_t)p->pipeline, p->label);

    md_mutex_lock(&dev->device_mutex);
    const bool listed = MD_ARRAY_TRY_PUSH(dev->pipelines, p);
    md_mutex_unlock(&dev->device_mutex);
    if (!listed) { md_vk_pipeline_free(dev, p); md_gpu_fail("out of memory"); return NULL; }
    return p;
}

void md_gpu_pipeline_destroy(md_gpu_pipeline_t p) {
    if (!p) return;
    md_gpu_device_t dev = p->device;
    md_mutex_lock(&dev->device_mutex);
    MD_ARRAY_REMOVE_ITEM(dev->pipelines, p);
    md_vk_retire_locked(dev, MD_VK_RETIRE_PIPELINE, p);
    md_mutex_unlock(&dev->device_mutex);
}

/* ---- Render passes ---------------------------------------------------------- */

/* The view a pass renders through for (mip, layer), made on first use and
   kept with the texture. */
static VkImageView md_vk_attach_view(md_gpu_texture_t t, uint32_t mip, uint32_t layer) {
    md_gpu_device_t dev = t->device;
    VkImageView view = VK_NULL_HANDLE;
    md_mutex_lock(&dev->device_mutex);
    for (uint32_t i = 0; i < t->attach_view_count; ++i) {
        if (t->attach_views[i].mip == mip && t->attach_views[i].layer == layer) {
            view = t->attach_views[i].view;
            md_mutex_unlock(&dev->device_mutex);
            return view;
        }
    }
    if (t->attach_view_count == t->attach_view_cap) {
        const uint32_t cap = t->attach_view_cap ? t->attach_view_cap * 2 : 2;
        md_vk_attach_view_t* arr = (md_vk_attach_view_t*)md_alloc(dev->alloc, cap * sizeof(md_vk_attach_view_t));
        if (!arr) { md_mutex_unlock(&dev->device_mutex); md_gpu_fail("out of memory"); return VK_NULL_HANDLE; }
        if (t->attach_views) {
            memcpy(arr, t->attach_views, t->attach_view_count * sizeof(md_vk_attach_view_t));
            md_free(dev->alloc, t->attach_views, t->attach_view_cap * sizeof(md_vk_attach_view_t));
        }
        t->attach_views    = arr;
        t->attach_view_cap = cap;
    }
    VkImageViewCreateInfo vci = {VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    vci.image    = t->image;
    vci.viewType = VK_IMAGE_VIEW_TYPE_2D;
    vci.format   = t->fi.format;
    vci.subresourceRange.aspectMask     = t->fi.view_aspect;
    vci.subresourceRange.baseMipLevel   = mip;
    vci.subresourceRange.levelCount     = 1;
    vci.subresourceRange.baseArrayLayer = layer;
    vci.subresourceRange.layerCount     = 1;
    if (md_vk_check(vkCreateImageView(dev->device, &vci, NULL, &view), "vkCreateImageView (attachment)")) {
        md_vk_attach_view_t* e = &t->attach_views[t->attach_view_count++];
        e->mip   = mip;
        e->layer = layer;
        e->view  = view;
    }
    md_mutex_unlock(&dev->device_mutex);
    return view;
}

static bool md_vk_check_attachment(md_gpu_stream_t s, md_gpu_texture_t t, uint32_t mip, uint32_t layer,
                                   bool depth, uint32_t index, uint32_t* w, uint32_t* h) {
    char what[32];
    if (depth) md_snprintf(what, sizeof(what), "depth attachment");
    else       md_snprintf(what, sizeof(what), "colour attachment %u", index);
    if (!t) return md_gpu_fail("md_gpu_render_begin: %s has no texture", what);
    if (t->device != s->device) return md_gpu_fail("md_gpu_render_begin: %s '%s' belongs to another device", what, t->label);
    if (!(t->desc.usage & MD_GPU_TEX_RENDER_TARGET)) {
        return md_gpu_fail("md_gpu_render_begin: %s '%s' lacks MD_GPU_TEX_RENDER_TARGET usage", what, t->label);
    }
    if (t->fi.depth != depth) {
        return md_gpu_fail("md_gpu_render_begin: %s '%s' is %s; %s", what, t->label, t->fi.name,
                          depth ? "a depth attachment needs a depth format" : "depth formats go in md_gpu_render_desc_t.depth");
    }
    if (mip >= t->desc.mip_levels) {
        return md_gpu_fail("md_gpu_render_begin: %s '%s': mip %u out of range (%u levels)", what, t->label, mip, t->desc.mip_levels);
    }
    const uint32_t layers = t->desc.type == MD_GPU_TEX_2D_ARRAY ? t->desc.depth_or_layers : 1u;
    if (layer >= layers) {
        return md_gpu_fail("md_gpu_render_begin: %s '%s': layer %u out of range (%u layers)", what, t->label, layer, layers);
    }
    uint32_t mw = t->desc.width  >> mip; if (!mw) mw = 1;
    uint32_t mh = t->desc.height >> mip; if (!mh) mh = 1;
    if (*w == 0) { *w = mw; *h = mh; }
    else if (*w != mw || *h != mh) {
        return md_gpu_fail("md_gpu_render_begin: %s '%s' is %ux%u at mip %u but the pass is %ux%u; every attachment must be the same size",
                          what, t->label, mw, mh, mip, *w, *h);
    }
    return true;
}

static VkAttachmentLoadOp md_vk_load_op(md_gpu_load_t l) {
    switch (l) {
    case MD_GPU_LOAD_CLEAR:     return VK_ATTACHMENT_LOAD_OP_CLEAR;
    case MD_GPU_LOAD_DONT_CARE: return VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    default:                    return VK_ATTACHMENT_LOAD_OP_LOAD;
    }
}

static VkAttachmentStoreOp md_vk_store_op(md_gpu_store_t st) {
    return st == MD_GPU_STORE_DISCARD ? VK_ATTACHMENT_STORE_OP_DONT_CARE : VK_ATTACHMENT_STORE_OP_STORE;
}

static VkCompareOp md_vk_compare_op(md_gpu_compare_t c) {
    switch (c) {
    case MD_GPU_COMPARE_NEVER:         return VK_COMPARE_OP_NEVER;
    case MD_GPU_COMPARE_LESS:          return VK_COMPARE_OP_LESS;
    case MD_GPU_COMPARE_LESS_EQUAL:    return VK_COMPARE_OP_LESS_OR_EQUAL;
    case MD_GPU_COMPARE_EQUAL:         return VK_COMPARE_OP_EQUAL;
    case MD_GPU_COMPARE_NOT_EQUAL:     return VK_COMPARE_OP_NOT_EQUAL;
    case MD_GPU_COMPARE_GREATER_EQUAL: return VK_COMPARE_OP_GREATER_OR_EQUAL;
    case MD_GPU_COMPARE_GREATER:       return VK_COMPARE_OP_GREATER;
    default:                           return VK_COMPARE_OP_ALWAYS;
    }
}

static bool md_vk_depth_test(const md_gpu_draw_state_t* d) {
    /* Vulkan writes depth only with the test enabled; ALWAYS makes the test
       pass, so "write without testing" is test-enabled + ALWAYS. */
    return d->depth_compare != MD_GPU_COMPARE_ALWAYS || d->depth_write;
}

static bool md_vk_depth_bias(const md_gpu_draw_state_t* d) {
    return d->depth_bias != 0.0f || d->depth_bias_slope != 0.0f;
}

/* Emit the difference between the state set in the command buffer and `n`
   (everything when `all`). */
static void md_vk_apply_draw_state(md_gpu_stream_t s, const md_gpu_draw_state_t* n, bool all) {
    VkCommandBuffer cmd = s->open;
    const md_gpu_draw_state_t* o = &s->draw_state;
    if (all || md_vk_depth_test(n) != md_vk_depth_test(o)) vkCmdSetDepthTestEnable(cmd, md_vk_depth_test(n) ? VK_TRUE : VK_FALSE);
    if (all || n->depth_write != o->depth_write)           vkCmdSetDepthWriteEnable(cmd, n->depth_write ? VK_TRUE : VK_FALSE);
    if (all || n->depth_compare != o->depth_compare)       vkCmdSetDepthCompareOp(cmd, md_vk_compare_op(n->depth_compare));
    if (all || n->cull != o->cull) {
        vkCmdSetCullMode(cmd, n->cull == MD_GPU_CULL_BACK  ? VK_CULL_MODE_BACK_BIT :
                              n->cull == MD_GPU_CULL_FRONT ? VK_CULL_MODE_FRONT_BIT : VK_CULL_MODE_NONE);
    }
    if (all || n->front_clockwise != o->front_clockwise) {
        vkCmdSetFrontFace(cmd, n->front_clockwise ? VK_FRONT_FACE_CLOCKWISE : VK_FRONT_FACE_COUNTER_CLOCKWISE);
    }
    if (all || md_vk_depth_bias(n) != md_vk_depth_bias(o)) vkCmdSetDepthBiasEnable(cmd, md_vk_depth_bias(n) ? VK_TRUE : VK_FALSE);
    if (all || n->depth_bias != o->depth_bias || n->depth_bias_slope != o->depth_bias_slope || n->depth_bias_clamp != o->depth_bias_clamp) {
        vkCmdSetDepthBias(cmd, n->depth_bias, n->depth_bias_clamp, n->depth_bias_slope);
    }
    if (all || memcmp(n->blend_constant, o->blend_constant, sizeof(n->blend_constant)) != 0) {
        vkCmdSetBlendConstants(cmd, n->blend_constant);
    }
    s->draw_state = *n;
}

static void md_vk_apply_viewport(md_gpu_stream_t s, const md_gpu_viewport_t* v) {
    VkViewport vp;
    vp.x        = v->x;
    vp.y        = v->y + v->height;      /* flipped: clip +Y up */
    vp.width    = v->width;
    vp.height   = -v->height;
    vp.minDepth = v->min_depth;
    vp.maxDepth = v->max_depth;
    if (v->min_depth == 0.0f && v->max_depth == 0.0f) vp.maxDepth = 1.0f;
    vkCmdSetViewport(s->open, 0, 1, &vp);
}

static void md_vk_apply_scissor(md_gpu_stream_t s, const md_gpu_rect_t* r) {
    /* Clamped to the render area: Metal requires it, and it is what the
       caller means either way. */
    uint32_t x = r->x < s->pass_width  ? r->x : s->pass_width;
    uint32_t y = r->y < s->pass_height ? r->y : s->pass_height;
    uint32_t w = r->width  < s->pass_width  - x ? r->width  : s->pass_width  - x;
    uint32_t h = r->height < s->pass_height - y ? r->height : s->pass_height - y;
    VkRect2D sc = {{(int32_t)x, (int32_t)y}, {w, h}};
    vkCmdSetScissor(s->open, 0, 1, &sc);
}

bool md_gpu_render_begin(md_gpu_stream_t s, const md_gpu_render_desc_t* desc) {
    if (!s || !desc) return md_gpu_fail("md_gpu_render_begin: null argument");
    if (!s->can_graphics) return md_gpu_fail("md_gpu_render_begin: stream '%s' is not a GRAPHICS stream", s->label);
    if (s->in_pass) return md_gpu_fail("md_gpu_render_begin: stream '%s' is already inside a render pass", s->label);
    if (s->upload_open) return md_gpu_fail("md_gpu_render_begin: stream '%s' has an open upload", s->label);
    if (desc->color_count > MD_GPU_MAX_COLOR_TARGETS) {
        return md_gpu_fail("md_gpu_render_begin: %u colour attachments, at most %u", desc->color_count, MD_GPU_MAX_COLOR_TARGETS);
    }
    if (desc->color_count == 0 && !desc->depth.texture) return md_gpu_fail("md_gpu_render_begin: the pass has no attachments");

    uint32_t w = 0, h = 0;
    VkRenderingAttachmentInfo ca[MD_GPU_MAX_COLOR_TARGETS];
    for (uint32_t i = 0; i < desc->color_count; ++i) {
        const md_gpu_color_attachment_t* a = &desc->color[i];
        if (!md_vk_check_attachment(s, a->texture, a->mip, a->layer, false, i, &w, &h)) return false;
        VkImageView view = md_vk_attach_view(a->texture, a->mip, a->layer);
        if (!view) return false;
        ca[i] = (VkRenderingAttachmentInfo){VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO};
        ca[i].imageView   = view;
        ca[i].imageLayout = VK_IMAGE_LAYOUT_GENERAL;
        ca[i].loadOp      = md_vk_load_op(a->load);
        ca[i].storeOp     = md_vk_store_op(a->store);
        memcpy(&ca[i].clearValue.color, &a->clear, sizeof(a->clear));
    }
    VkRenderingAttachmentInfo da = {VK_STRUCTURE_TYPE_RENDERING_ATTACHMENT_INFO};
    if (desc->depth.texture) {
        const md_gpu_depth_attachment_t* a = &desc->depth;
        if (!md_vk_check_attachment(s, a->texture, a->mip, a->layer, true, 0, &w, &h)) return false;
        VkImageView view = md_vk_attach_view(a->texture, a->mip, a->layer);
        if (!view) return false;
        da.imageView   = view;
        da.imageLayout = VK_IMAGE_LAYOUT_GENERAL;
        da.loadOp      = md_vk_load_op(a->load);
        da.storeOp     = md_vk_store_op(a->store);
        da.clearValue.depthStencil.depth = a->clear_depth;
    }

    /* A pass reads in every raster stage: indirect commands, indices, vertex
       and fragment shaders, and its attachments' load actions. */
    VkCommandBuffer cmd = md_vk_begin_op(s, MD_GPU_STAGE_INDIRECT | MD_GPU_STAGE_VERTEX | MD_GPU_STAGE_FRAGMENT | MD_GPU_STAGE_ATTACHMENT);
    if (!cmd) return false;
    md_gpu_device_t dev = s->device;
    s->pass_labelled = desc->label && md_vk_debug_labels(dev);
    if (s->pass_labelled) {
        VkDebugUtilsLabelEXT l = {VK_STRUCTURE_TYPE_DEBUG_UTILS_LABEL_EXT};
        l.pLabelName = desc->label;
        vkCmdBeginDebugUtilsLabelEXT(cmd, &l);
    }

    VkRenderingInfo ri = {VK_STRUCTURE_TYPE_RENDERING_INFO};
    ri.renderArea.extent.width  = w;
    ri.renderArea.extent.height = h;
    ri.layerCount           = 1;
    ri.colorAttachmentCount = desc->color_count;
    ri.pColorAttachments    = ca;
    ri.pDepthAttachment     = desc->depth.texture ? &da : NULL;
    vkCmdBeginRendering(cmd, &ri);
    if (!s->graphics_set_bound) {
        vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, dev->raster_layout, 0, 1, &dev->desc_set, 0, NULL);
        s->graphics_set_bound = true;
    }

    s->in_pass          = true;
    s->has_work         = true;
    s->pass_width       = w;
    s->pass_height      = h;
    s->pass_color_count = desc->color_count;
    for (uint32_t i = 0; i < MD_GPU_MAX_COLOR_TARGETS; ++i) {
        s->pass_color[i] = i < desc->color_count ? desc->color[i].texture->desc.format : MD_GPU_FORMAT_INVALID;
    }
    s->pass_depth       = desc->depth.texture ? desc->depth.texture->desc.format : MD_GPU_FORMAT_INVALID;
    s->bound_pipeline   = NULL;

    md_gpu_draw_state_t zero;
    memset(&zero, 0, sizeof(zero));
    md_vk_apply_draw_state(s, &zero, true);
    md_vk_apply_viewport(s, &(md_gpu_viewport_t){0, 0, (float)w, (float)h, 0, 0});
    md_vk_apply_scissor(s, &(md_gpu_rect_t){0, 0, w, h});
    return true;
}

bool md_gpu_render_end(md_gpu_stream_t s) {
    if (!s) return md_gpu_fail("md_gpu_render_end: null stream");
    if (!s->in_pass) return md_gpu_fail("md_gpu_render_end: stream '%s' has no open render pass", s->label);
    vkCmdEndRendering(s->open);
    if (s->pass_labelled) vkCmdEndDebugUtilsLabelEXT(s->open);
    s->in_pass        = false;
    s->pass_labelled  = false;
    s->bound_pipeline = NULL;
    md_vk_end_op(s, MD_GPU_STAGE_VERTEX | MD_GPU_STAGE_FRAGMENT | MD_GPU_STAGE_ATTACHMENT);
    return true;
}

/* Teardown with a pass still open: close it so the command buffer can end. */
static void md_vk_abandon_pass(md_gpu_stream_t s) {
    if (s->in_pass) md_gpu_render_end(s);
}

void md_gpu_set_draw_state(md_gpu_stream_t s, const md_gpu_draw_state_t* state) {
    if (!s) return;
    if (!s->in_pass) { md_gpu_fail("md_gpu_set_draw_state: stream '%s' has no open render pass", s->label); return; }
    md_gpu_draw_state_t n;
    memset(&n, 0, sizeof(n));
    if (state) n = *state;
    if ((unsigned)n.depth_compare > (unsigned)MD_GPU_COMPARE_GREATER || (unsigned)n.cull > (unsigned)MD_GPU_CULL_FRONT) {
        md_gpu_fail("md_gpu_set_draw_state: invalid depth_compare or cull");
        return;
    }
    if (n.depth_bias_clamp != 0.0f && !s->device->depth_bias_clamp) {
        md_gpu_fail("md_gpu_set_draw_state: the device does not support depth_bias_clamp; leave it 0");
        return;
    }
    md_vk_apply_draw_state(s, &n, false);
}

void md_gpu_set_viewport(md_gpu_stream_t s, const md_gpu_viewport_t* viewport) {
    if (!s) return;
    if (!s->in_pass) { md_gpu_fail("md_gpu_set_viewport: stream '%s' has no open render pass", s->label); return; }
    md_gpu_viewport_t v = {0, 0, (float)s->pass_width, (float)s->pass_height, 0, 0};
    if (viewport) v = *viewport;
    if (!(v.width > 0.0f) || !(v.height > 0.0f)) { md_gpu_fail("md_gpu_set_viewport: width and height must be positive"); return; }
    if (v.min_depth < 0.0f || v.min_depth > 1.0f || v.max_depth < 0.0f || v.max_depth > 1.0f) {
        md_gpu_fail("md_gpu_set_viewport: depth range must lie in [0, 1]");
        return;
    }
    md_vk_apply_viewport(s, &v);
}

void md_gpu_set_scissor(md_gpu_stream_t s, const md_gpu_rect_t* scissor) {
    if (!s) return;
    if (!s->in_pass) { md_gpu_fail("md_gpu_set_scissor: stream '%s' has no open render pass", s->label); return; }
    md_gpu_rect_t r = {0, 0, s->pass_width, s->pass_height};
    if (scissor) r = *scissor;
    md_vk_apply_scissor(s, &r);
}

/* ---- Draws ------------------------------------------------------------------ */

/* Every draw goes through md_vk_draw: one kind of four, its arguments copied
   or by address, validated before anything is recorded. */
typedef enum md_vk_draw_kind_t {
    MD_VK_DRAW,
    MD_VK_DRAW_INDEXED,
    MD_VK_DRAW_INDIRECT,
    MD_VK_DRAW_INDEXED_INDIRECT,
} md_vk_draw_kind_t;

typedef struct md_vk_draw_t {
    md_vk_draw_kind_t   kind;
    uint32_t            count;            /* vertices, indices, or indirect commands */
    uint32_t            instance_count;   /* direct draws                            */
    md_gpu_addr_t       indices;
    md_gpu_index_type_t index_type;
    md_gpu_addr_t       cmds;
    /* The argument struct: copied from `args` (args_size bytes), or at
       `args_at` when by_address. */
    bool                by_address;
    const void*         args;
    size_t              args_size;
    md_gpu_addr_t       args_at;
} md_vk_draw_t;

static const char* md_vk_draw_name(const md_vk_draw_t* d) {
    static const char* names[2][4] = {
        {"md_gpu_draw", "md_gpu_draw_indexed", "md_gpu_draw_indirect", "md_gpu_draw_indexed_indirect"},
        {"md_gpu_draw_at", "md_gpu_draw_indexed_at", "md_gpu_draw_indirect_at", "md_gpu_draw_indexed_indirect_at"},
    };
    return names[d->by_address ? 1 : 0][d->kind];
}

/* The pass and the pipeline agree; nothing is recorded. */
static bool md_vk_draw_check(md_gpu_stream_t s, md_gpu_pipeline_t p, const char* what) {
    if (!s || !p) return md_gpu_fail("%s: null stream or pipeline", what);
    if (p->device != s->device) return md_gpu_fail("%s: pipeline '%s' belongs to another device", what, p->label);
    if (!s->in_pass) return md_gpu_fail("%s: stream '%s' has no open render pass", what, s->label);
    if (p->color_count != s->pass_color_count) {
        return md_gpu_fail("%s: pipeline '%s' has %u colour targets but the pass has %u attachments",
                          what, p->label, p->color_count, s->pass_color_count);
    }
    for (uint32_t i = 0; i < p->color_count; ++i) {
        if (p->color[i] != s->pass_color[i]) {
            return md_gpu_fail("%s: pipeline '%s' colour target %u is %s but the pass attachment is %s",
                              what, p->label, i, md_vk_format_info(p->color[i]).name, md_vk_format_info(s->pass_color[i]).name);
        }
    }
    if (p->depth != s->pass_depth) {
        return md_gpu_fail("%s: pipeline '%s' depth format is %s but the pass has %s", what, p->label,
                          p->depth ? md_vk_format_info(p->depth).name : "none",
                          s->pass_depth ? md_vk_format_info(s->pass_depth).name : "none");
    }
    return true;
}

static bool md_vk_resolve_indices(md_gpu_stream_t s, md_gpu_addr_t indices, md_gpu_index_type_t type, uint64_t count,
                                  md_vk_span_t* out, const char* what) {
    if (type != MD_GPU_INDEX_U32 && type != MD_GPU_INDEX_U16) return md_gpu_fail("%s: invalid index type %d", what, (int)type);
    const uint64_t isize = type == MD_GPU_INDEX_U16 ? 2 : 4;
    if (!indices) return md_gpu_fail("%s: null index address", what);
    if (indices % isize != 0) return md_gpu_fail("%s: index address must be %llu-byte aligned", what, (unsigned long long)isize);
    return md_vk_resolve(s->device, indices, count * isize, out, what);
}

static bool md_vk_draw(md_gpu_stream_t s, md_gpu_pipeline_t p, const md_vk_draw_t* d) {
    const char* what = md_vk_draw_name(d);
    if (!md_vk_draw_check(s, p, what)) return false;
    const bool indirect = d->kind == MD_VK_DRAW_INDIRECT || d->kind == MD_VK_DRAW_INDEXED_INDIRECT;
    const bool indexed  = d->kind == MD_VK_DRAW_INDEXED  || d->kind == MD_VK_DRAW_INDEXED_INDIRECT;
    if (d->count == 0 || (!indirect && d->instance_count == 0)) return true;

    /* Indirect commands are tightly packed in the layout both APIs read. The
       index range of an indirect draw is known only to the GPU; the address
       must at least start a live allocation. */
    const uint32_t stride = d->kind == MD_VK_DRAW_INDIRECT ? (uint32_t)sizeof(md_gpu_draw_cmd_t)
                                                           : (uint32_t)sizeof(md_gpu_draw_indexed_cmd_t);
    md_vk_span_t ib = {0}, cb = {0};
    if (indexed && !md_vk_resolve_indices(s, d->indices, d->index_type, indirect ? 1 : d->count, &ib, what)) return false;
    if (indirect) {
        if (!d->cmds) return md_gpu_fail("%s: null command address", what);
        if (d->cmds % 4 != 0) return md_gpu_fail("%s: command address must be 4-byte aligned", what);
        if (!md_vk_resolve(s->device, d->cmds, (uint64_t)d->count * stride, &cb, what)) return false;
    }

    uint64_t args_addr = d->args_at;
    if (d->by_address) {
        if (!md_vk_args_at(s, p->args_size, d->args_at, what)) return false;
    } else if (!md_vk_args_copy(s, p->args_size, d->args, d->args_size, &args_addr, what, p->label)) {
        return false;
    }

    VkCommandBuffer cmd = s->open;
    if (s->bound_pipeline != p) {
        vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, p->pipeline);
        s->bound_pipeline = p;
    }
    vkCmdPushConstants(cmd, s->device->raster_layout, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, 8, &args_addr);
    if (indexed) vkCmdBindIndexBuffer(cmd, ib.buffer, ib.offset, d->index_type == MD_GPU_INDEX_U16 ? VK_INDEX_TYPE_UINT16 : VK_INDEX_TYPE_UINT32);

    switch (d->kind) {
    case MD_VK_DRAW:         vkCmdDraw(cmd, d->count, d->instance_count, 0, 0); break;
    case MD_VK_DRAW_INDEXED: vkCmdDrawIndexed(cmd, d->count, d->instance_count, 0, 0, 0); break;
    default: {
        const uint32_t max = s->device->props.limits.maxDrawIndirectCount;
        for (uint32_t done = 0; done < d->count;) {
            const uint32_t n = d->count - done < max ? d->count - done : max;
            const uint64_t off = cb.offset + (uint64_t)done * stride;
            if (indexed) vkCmdDrawIndexedIndirect(cmd, cb.buffer, off, n, stride);
            else         vkCmdDrawIndirect(cmd, cb.buffer, off, n, stride);
            done += n;
        }
    } break;
    }
    return true;
}

bool md_gpu_draw(md_gpu_stream_t s, md_gpu_pipeline_t p, uint32_t vertex_count, uint32_t instance_count,
                 const void* args, size_t args_size) {
    md_vk_draw_t d = {MD_VK_DRAW, vertex_count, instance_count};
    d.args = args; d.args_size = args_size;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_at(md_gpu_stream_t s, md_gpu_pipeline_t p, uint32_t vertex_count, uint32_t instance_count, md_gpu_addr_t args) {
    md_vk_draw_t d = {MD_VK_DRAW, vertex_count, instance_count};
    d.by_address = true; d.args_at = args;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indexed(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                         uint32_t index_count, uint32_t instance_count, const void* args, size_t args_size) {
    md_vk_draw_t d = {MD_VK_DRAW_INDEXED, index_count, instance_count, indices, index_type};
    d.args = args; d.args_size = args_size;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indexed_at(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                            uint32_t index_count, uint32_t instance_count, md_gpu_addr_t args) {
    md_vk_draw_t d = {MD_VK_DRAW_INDEXED, index_count, instance_count, indices, index_type};
    d.by_address = true; d.args_at = args;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indirect(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t cmds, uint32_t count,
                          const void* args, size_t args_size) {
    md_vk_draw_t d = {MD_VK_DRAW_INDIRECT, count};
    d.cmds = cmds; d.args = args; d.args_size = args_size;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indirect_at(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t cmds, uint32_t count, md_gpu_addr_t args) {
    md_vk_draw_t d = {MD_VK_DRAW_INDIRECT, count};
    d.cmds = cmds; d.by_address = true; d.args_at = args;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indexed_indirect(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                                  md_gpu_addr_t cmds, uint32_t count, const void* args, size_t args_size) {
    md_vk_draw_t d = {MD_VK_DRAW_INDEXED_INDIRECT, count, 0, indices, index_type, cmds};
    d.args = args; d.args_size = args_size;
    return md_vk_draw(s, p, &d);
}

bool md_gpu_draw_indexed_indirect_at(md_gpu_stream_t s, md_gpu_pipeline_t p, md_gpu_addr_t indices, md_gpu_index_type_t index_type,
                                     md_gpu_addr_t cmds, uint32_t count, md_gpu_addr_t args) {
    md_vk_draw_t d = {MD_VK_DRAW_INDEXED_INDIRECT, count, 0, indices, index_type, cmds};
    d.by_address = true; d.args_at = args;
    return md_vk_draw(s, p, &d);
}

/* =========================================================================
   13. Presentation
   =========================================================================

   One swapchain per surface, rebuilt when acquire or present reports it out
   of date, or when the size changes. Images get the same treatment as every
   other image -- GENERAL for their whole life in md_gpu -- with two extra
   transitions: UNDEFINED -> GENERAL recorded at acquire (contents are
   undefined at acquire anyway) and GENERAL -> PRESENT_SRC recorded at present.

   Semaphores: acquire signals a binary semaphore that the next submission
   of the acquiring stream waits on; present has that stream's submission
   signal a per-image semaphore that vkQueuePresentKHR waits on. Acquire
   semaphores are reused round robin once the submission that consumed them
   has completed on the stream's timeline. A per-image present semaphore is
   free again when its image is next acquired. */

typedef struct md_vk_win32_surface_ci_t {
    VkStructureType sType; const void* pNext; VkFlags flags; void* hinstance; void* hwnd;
} md_vk_win32_surface_ci_t;
typedef struct md_vk_xlib_surface_ci_t {
    VkStructureType sType; const void* pNext; VkFlags flags; void* dpy; unsigned long window;
} md_vk_xlib_surface_ci_t;
typedef struct md_vk_wayland_surface_ci_t {
    VkStructureType sType; const void* pNext; VkFlags flags; void* display; void* surface;
} md_vk_wayland_surface_ci_t;
typedef VkResult (VKAPI_PTR *md_vk_create_surface_fn)(VkInstance, const void*, const VkAllocationCallbacks*, VkSurfaceKHR*);

#define MD_VK_STYPE_WIN32_SURFACE   ((VkStructureType)1000009000)
#define MD_VK_STYPE_XLIB_SURFACE    ((VkStructureType)1000004000)
#define MD_VK_STYPE_WAYLAND_SURFACE ((VkStructureType)1000006000)

/* Every image of `sc`, and the swapchain itself, once nothing uses them. The
   presentation engine's semaphore waits are not on any timeline, so the
   queue that presented is idled first -- this runs on retirement, after a
   rebuild or surface destruction, never per frame. Caller holds device_mutex. */
static void md_vk_swapchain_free(md_gpu_device_t dev, md_vk_swapchain_t* sc) {
    if (sc->present_queue) {
        md_mutex_lock(&dev->queue_mutex);
        vkQueueWaitIdle(sc->present_queue);
        md_mutex_unlock(&dev->queue_mutex);
    }
    for (uint32_t i = 0; i < sc->image_count; ++i) {
        if (sc->textures && sc->textures[i]) md_vk_texture_free(dev, sc->textures[i]);
        if (sc->present_sems && sc->present_sems[i]) vkDestroySemaphore(dev->device, sc->present_sems[i], NULL);
    }
    if (sc->textures)     md_free(dev->alloc, sc->textures, sc->image_count * sizeof(md_gpu_texture_t));
    if (sc->present_sems) md_free(dev->alloc, sc->present_sems, sc->image_count * sizeof(VkSemaphore));
    if (sc->swapchain) vkDestroySwapchainKHR(dev->device, sc->swapchain, NULL);
    for (uint32_t i = 0; i < sc->extra_sem_count; ++i) if (sc->extra_sems[i]) vkDestroySemaphore(dev->device, sc->extra_sems[i], NULL);
    if (sc->extra_sems) md_free(dev->alloc, sc->extra_sems, sc->extra_sem_count * sizeof(VkSemaphore));
    if (sc->surface) vkDestroySurfaceKHR(dev->instance, sc->surface, NULL);
    md_free(dev->alloc, sc, sizeof(*sc));
}

/* A texture wrapping one swapchain image. Caller holds device_mutex. */
static md_gpu_texture_t md_vk_wrap_image_locked(md_gpu_surface_t sf, VkImage image, uint32_t w, uint32_t h, uint32_t index) {
    md_gpu_device_t dev = sf->device;
    md_gpu_texture_t t = (md_gpu_texture_t)md_alloc(dev->alloc, sizeof(md_gpu_texture));
    if (!t) { md_gpu_fail("out of memory"); return NULL; }
    memset(t, 0, sizeof(*t));
    t->device   = dev;
    t->image    = image;
    t->external = true;
    t->fi       = md_vk_format_info(sf->desc.format);
    t->desc.type            = MD_GPU_TEX_2D;
    t->desc.format          = sf->desc.format;
    t->desc.usage           = MD_GPU_TEX_RENDER_TARGET | sf->desc.usage;
    t->desc.width           = w;
    t->desc.height          = h;
    t->desc.depth_or_layers = 1;
    t->desc.mip_levels      = 1;
    md_snprintf(t->label, sizeof(t->label), "%s[%u]", sf->label, index);
    t->desc.label = t->label;

    if (sf->desc.usage & MD_GPU_TEX_SAMPLED) {
        t->sampled_view = md_vk_create_view(dev, t, 0, 1);
        if (!t->sampled_view) goto fail;
        t->sampled_slot = md_vk_alloc_slot_locked(dev);
        if (!t->sampled_slot) { md_gpu_fail("surface '%s': out of bindless heap slots", sf->label); goto fail; }
        md_vk_write_image_slot(dev, t->sampled_slot, t->sampled_view, VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE);
    }
    if (sf->desc.usage & MD_GPU_TEX_STORAGE) {
        t->storage_views = (VkImageView*)md_alloc(dev->alloc, sizeof(VkImageView));
        if (t->storage_views) t->storage_views[0] = VK_NULL_HANDLE;
        t->storage_slots = (uint32_t*)md_alloc(dev->alloc, sizeof(uint32_t));
        if (t->storage_slots) t->storage_slots[0] = 0;
        if (!t->storage_views || !t->storage_slots) { md_gpu_fail("out of memory"); goto fail; }
        t->storage_views[0] = md_vk_create_view(dev, t, 0, 1);
        if (!t->storage_views[0]) goto fail;
        t->storage_slots[0] = md_vk_alloc_slot_locked(dev);
        if (!t->storage_slots[0]) { md_gpu_fail("surface '%s': out of bindless heap slots", sf->label); goto fail; }
        md_vk_write_image_slot(dev, t->storage_slots[0], t->storage_views[0], VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    }
    return t;
fail:
    md_vk_texture_free(dev, t);
    return NULL;
}

/* Build (or rebuild) the swapchain at the current size. *zero_size is set,
   and nothing built, while the drawable has no area. */
static bool md_vk_surface_rebuild(md_gpu_surface_t sf, bool* zero_size) {
    md_gpu_device_t dev = sf->device;
    *zero_size = false;
    VkSurfaceCapabilitiesKHR caps;
    if (!md_vk_check(vkGetPhysicalDeviceSurfaceCapabilitiesKHR(dev->phys, sf->surface, &caps),
                     "vkGetPhysicalDeviceSurfaceCapabilitiesKHR")) return false;

    md_mutex_lock(&dev->device_mutex);
    uint32_t w = sf->want_width, h = sf->want_height;
    md_mutex_unlock(&dev->device_mutex);
    VkExtent2D ext;
    if (caps.currentExtent.width != UINT32_MAX) {
        ext = caps.currentExtent;       /* the window decides */
    } else {
        ext.width  = w < caps.minImageExtent.width  ? caps.minImageExtent.width  : (w > caps.maxImageExtent.width  ? caps.maxImageExtent.width  : w);
        ext.height = h < caps.minImageExtent.height ? caps.minImageExtent.height : (h > caps.maxImageExtent.height ? caps.maxImageExtent.height : h);
        if (w == 0 || h == 0) ext.width = ext.height = 0;
    }
    if (ext.width == 0 || ext.height == 0) { *zero_size = true; return true; }

    uint32_t count = caps.minImageCount + 1;
    if (caps.maxImageCount && count > caps.maxImageCount) count = caps.maxImageCount;

    VkCompositeAlphaFlagBitsKHR alpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    if (!(caps.supportedCompositeAlpha & alpha)) {
        for (uint32_t b = 0; b < 32; ++b) {
            if (caps.supportedCompositeAlpha & (1u << b)) { alpha = (VkCompositeAlphaFlagBitsKHR)(1u << b); break; }
        }
    }

    VkSwapchainCreateInfoKHR sci = {VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sci.surface          = sf->surface;
    sci.minImageCount    = count;
    sci.imageFormat      = sf->vk_format.format;
    sci.imageColorSpace  = sf->vk_format.colorSpace;
    sci.imageExtent      = ext;
    sci.imageArrayLayers = 1;
    sci.imageUsage       = sf->image_usage;
    md_vk_set_sharing(dev, &sci.imageSharingMode, &sci.queueFamilyIndexCount, &sci.pQueueFamilyIndices);
    sci.preTransform     = caps.currentTransform;
    sci.compositeAlpha   = alpha;
    sci.presentMode      = sf->present_mode;
    sci.clipped          = VK_TRUE;
    sci.oldSwapchain     = sf->sc ? sf->sc->swapchain : VK_NULL_HANDLE;

    md_vk_swapchain_t* sc = (md_vk_swapchain_t*)md_alloc(dev->alloc, sizeof(md_vk_swapchain_t));
    if (!sc) return md_gpu_fail("out of memory");
    memset(sc, 0, sizeof(*sc));
    sc->width  = ext.width;
    sc->height = ext.height;
    if (!md_vk_check(vkCreateSwapchainKHR(dev->device, &sci, NULL, &sc->swapchain), "vkCreateSwapchainKHR")) {
        md_free(dev->alloc, sc, sizeof(*sc));
        return false;
    }

    uint32_t n = 0;
    vkGetSwapchainImagesKHR(dev->device, sc->swapchain, &n, NULL);
    VkImage images[16];
    if (n > 16) n = 16;
    VkResult r = vkGetSwapchainImagesKHR(dev->device, sc->swapchain, &n, images);
    bool ok = r == VK_SUCCESS || r == VK_INCOMPLETE;
    if (!ok) md_vk_check(r, "vkGetSwapchainImagesKHR");

    md_mutex_lock(&dev->device_mutex);
    if (ok) {
        sc->textures     = (md_gpu_texture_t*)md_alloc(dev->alloc, n * sizeof(md_gpu_texture_t));
        sc->present_sems = (VkSemaphore*)md_alloc(dev->alloc, n * sizeof(VkSemaphore));
        if (!sc->textures || !sc->present_sems) {
            if (sc->textures)     md_free(dev->alloc, sc->textures, n * sizeof(md_gpu_texture_t));
            if (sc->present_sems) md_free(dev->alloc, sc->present_sems, n * sizeof(VkSemaphore));
            sc->textures = NULL; sc->present_sems = NULL;
            ok = md_gpu_fail("out of memory");
        } else {
            sc->image_count = n;
            memset(sc->textures, 0, n * sizeof(md_gpu_texture_t));
            memset(sc->present_sems, 0, n * sizeof(VkSemaphore));
        }
    }
    for (uint32_t i = 0; i < sc->image_count && ok; ++i) {
        sc->textures[i] = md_vk_wrap_image_locked(sf, images[i], ext.width, ext.height, i);
        VkSemaphoreCreateInfo semci = {VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        ok = sc->textures[i] && md_vk_check(vkCreateSemaphore(dev->device, &semci, NULL, &sc->present_sems[i]), "vkCreateSemaphore");
    }
    if (!ok) {
        md_vk_swapchain_free(dev, sc);
        md_mutex_unlock(&dev->device_mutex);
        return false;
    }
    /* The old swapchain goes once every stream has passed the work issued so
       far; its images may still be in flight. */
    if (sf->sc) {
        if (!sf->sc->present_queue) sf->sc->present_queue = sf->present_queue;
        md_vk_retire_locked(dev, MD_VK_RETIRE_SWAPCHAIN, sf->sc);
    }
    sf->sc    = sc;
    sf->dirty = false;
    md_mutex_unlock(&dev->device_mutex);
    return true;
}

static md_gpu_format_t md_vk_surface_format_from_vk(VkFormat f) {
    for (int i = 1; i < (int)MD_GPU_FORMAT_COUNT; ++i) {
        if (md_vk_format_info((md_gpu_format_t)i).format == f) return (md_gpu_format_t)i;
    }
    return MD_GPU_FORMAT_INVALID;
}

md_gpu_surface_t md_gpu_surface_create(md_gpu_device_t dev, const md_gpu_surface_desc_t* desc) {
    if (!dev || !desc) { md_gpu_fail("md_gpu_surface_create: null argument"); return NULL; }
    const char* label = desc->label ? desc->label : "surface";
    if (!dev->supports_present) {
        md_gpu_fail("surface '%s': the device cannot present (VK_KHR_surface / VK_KHR_swapchain unavailable)", label);
        return NULL;
    }
    if (desc->system != MD_GPU_WINDOW_HEADLESS && !desc->window) {
        md_gpu_fail("surface '%s': no window handle", label);
        return NULL;
    }
    if (desc->usage & ~(MD_GPU_TEX_SAMPLED | MD_GPU_TEX_STORAGE | MD_GPU_TEX_RENDER_TARGET)) {
        md_gpu_fail("surface '%s': invalid usage bits", label);
        return NULL;
    }

    VkSurfaceKHR vs = VK_NULL_HANDLE;
    VkResult r = VK_ERROR_EXTENSION_NOT_PRESENT;
    const char* ext_missing = NULL;
    switch (desc->system) {
    case MD_GPU_WINDOW_WIN32: {
        md_vk_create_surface_fn fn = (md_vk_create_surface_fn)vkGetInstanceProcAddr(dev->instance, "vkCreateWin32SurfaceKHR");
        if (!dev->has_win32_surface || !fn) { ext_missing = "VK_KHR_win32_surface"; break; }
        md_vk_win32_surface_ci_t ci = {MD_VK_STYPE_WIN32_SURFACE, NULL, 0, desc->display, desc->window};
        if (!ci.hinstance) ci.hinstance = md_os_module_handle();
        r = fn(dev->instance, &ci, NULL, &vs);
        break;
    }
    case MD_GPU_WINDOW_X11: {
        md_vk_create_surface_fn fn = (md_vk_create_surface_fn)vkGetInstanceProcAddr(dev->instance, "vkCreateXlibSurfaceKHR");
        if (!dev->has_xlib_surface || !fn) { ext_missing = "VK_KHR_xlib_surface"; break; }
        if (!desc->display) { md_gpu_fail("surface '%s': X11 needs the Display* in `display`", label); return NULL; }
        md_vk_xlib_surface_ci_t ci = {MD_VK_STYPE_XLIB_SURFACE, NULL, 0, desc->display, (unsigned long)(uintptr_t)desc->window};
        r = fn(dev->instance, &ci, NULL, &vs);
        break;
    }
    case MD_GPU_WINDOW_WAYLAND: {
        md_vk_create_surface_fn fn = (md_vk_create_surface_fn)vkGetInstanceProcAddr(dev->instance, "vkCreateWaylandSurfaceKHR");
        if (!dev->has_wayland_surface || !fn) { ext_missing = "VK_KHR_wayland_surface"; break; }
        if (!desc->display) { md_gpu_fail("surface '%s': Wayland needs the wl_display* in `display`", label); return NULL; }
        md_vk_wayland_surface_ci_t ci = {MD_VK_STYPE_WAYLAND_SURFACE, NULL, 0, desc->display, desc->window};
        r = fn(dev->instance, &ci, NULL, &vs);
        break;
    }
    case MD_GPU_WINDOW_HEADLESS: {
        PFN_vkCreateHeadlessSurfaceEXT fn = (PFN_vkCreateHeadlessSurfaceEXT)vkGetInstanceProcAddr(dev->instance, "vkCreateHeadlessSurfaceEXT");
        if (!dev->has_headless_surface || !fn) { ext_missing = VK_EXT_HEADLESS_SURFACE_EXTENSION_NAME; break; }
        VkHeadlessSurfaceCreateInfoEXT ci = {VK_STRUCTURE_TYPE_HEADLESS_SURFACE_CREATE_INFO_EXT};
        r = fn(dev->instance, &ci, NULL, &vs);
        break;
    }
    case MD_GPU_WINDOW_COCOA:
        md_gpu_fail("surface '%s': COCOA windows present through the Metal backend", label);
        return NULL;
    default:
        md_gpu_fail("surface '%s': invalid window system %d", label, (int)desc->system);
        return NULL;
    }
    if (ext_missing) { md_gpu_fail("surface '%s': the Vulkan instance lacks %s", label, ext_missing); return NULL; }
    if (!md_vk_check(r, "vkCreate*SurfaceKHR")) return NULL;

    md_gpu_surface_t sf = (md_gpu_surface_t)md_alloc(dev->alloc, sizeof(md_gpu_surface));
    if (!sf) { vkDestroySurfaceKHR(dev->instance, vs, NULL); md_gpu_fail("out of memory"); return NULL; }
    memset(sf, 0, sizeof(*sf));
    sf->device  = dev;
    sf->surface = vs;
    sf->desc    = *desc;
    md_snprintf(sf->label, sizeof(sf->label), "%s", label);
    sf->desc.label = sf->label;
    if (sf->desc.format == MD_GPU_FORMAT_INVALID) sf->desc.format = MD_GPU_FORMAT_BGRA8_UNORM;
    sf->desc.usage &= ~MD_GPU_TEX_RENDER_TARGET;
    sf->want_width  = desc->width;
    sf->want_height = desc->height;
    sf->dirty       = true;
    const md_vk_format_info_t fi = md_vk_format_info(sf->desc.format);

    VkBool32 supported = VK_FALSE;
    vkGetPhysicalDeviceSurfaceSupportKHR(dev->phys, dev->graphics_family, vs, &supported);
    if (!supported) { md_gpu_fail("surface '%s': the graphics queue cannot present to it", label); goto fail; }

    /* Format: exactly the one asked for, no silent substitute. */
    {
        uint32_t n = 0;
        vkGetPhysicalDeviceSurfaceFormatsKHR(dev->phys, vs, &n, NULL);
        VkSurfaceFormatKHR fmts[64];
        if (n > 64) n = 64;
        vkGetPhysicalDeviceSurfaceFormatsKHR(dev->phys, vs, &n, fmts);
        bool found = false;
        for (uint32_t i = 0; i < n; ++i) {
            if (fmts[i].format != fi.format) continue;
            if (!found || fmts[i].colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR) sf->vk_format = fmts[i];
            found = true;
        }
        if (!found || fi.format == VK_FORMAT_UNDEFINED) {
            char avail[256] = {0};
            size_t len = 0;
            for (uint32_t i = 0; i < n && len < sizeof(avail) - 32; ++i) {
                md_gpu_format_t mf = md_vk_surface_format_from_vk(fmts[i].format);
                if (mf == MD_GPU_FORMAT_INVALID) continue;
                len += (size_t)md_snprintf(avail + len, sizeof(avail) - len, "%s%s", len ? ", " : "", md_vk_format_info(mf).name);
            }
            md_gpu_fail("surface '%s': %s is not supported (the surface offers: %s)", label, fi.name, len ? avail : "nothing md_gpu knows");
            goto fail;
        }
    }

    /* Usage: RENDER_TARGET always; transfers when offered (copies to and
       from the image); SAMPLED and STORAGE when asked, or fail. */
    {
        VkSurfaceCapabilitiesKHR caps;
        if (!md_vk_check(vkGetPhysicalDeviceSurfaceCapabilitiesKHR(dev->phys, vs, &caps), "vkGetPhysicalDeviceSurfaceCapabilitiesKHR")) goto fail;
        VkFormatProperties fp;
        vkGetPhysicalDeviceFormatProperties(dev->phys, fi.format, &fp);
        sf->image_usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
        sf->image_usage |= caps.supportedUsageFlags & (VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT);
        if (!(caps.supportedUsageFlags & VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT)) {
            md_gpu_fail("surface '%s': images cannot be rendered to", label);
            goto fail;
        }
        if (sf->desc.usage & MD_GPU_TEX_SAMPLED) {
            if (!(caps.supportedUsageFlags & VK_IMAGE_USAGE_SAMPLED_BIT) || !(fp.optimalTilingFeatures & VK_FORMAT_FEATURE_SAMPLED_IMAGE_BIT)) {
                md_gpu_fail("surface '%s': SAMPLED usage is not supported for %s", label, fi.name);
                goto fail;
            }
            sf->image_usage |= VK_IMAGE_USAGE_SAMPLED_BIT;
        }
        if (sf->desc.usage & MD_GPU_TEX_STORAGE) {
            if (!(caps.supportedUsageFlags & VK_IMAGE_USAGE_STORAGE_BIT) || !(fp.optimalTilingFeatures & VK_FORMAT_FEATURE_STORAGE_IMAGE_BIT)) {
                md_gpu_fail("surface '%s': STORAGE usage is not supported for %s", label, fi.name);
                goto fail;
            }
            sf->image_usage |= VK_IMAGE_USAGE_STORAGE_BIT;
        }
    }

    /* Present mode, with the documented fallbacks. */
    {
        uint32_t n = 0;
        vkGetPhysicalDeviceSurfacePresentModesKHR(dev->phys, vs, &n, NULL);
        VkPresentModeKHR modes[16];
        if (n > 16) n = 16;
        vkGetPhysicalDeviceSurfacePresentModesKHR(dev->phys, vs, &n, modes);
        bool has_mailbox = false, has_immediate = false;
        for (uint32_t i = 0; i < n; ++i) {
            if (modes[i] == VK_PRESENT_MODE_MAILBOX_KHR)   has_mailbox   = true;
            if (modes[i] == VK_PRESENT_MODE_IMMEDIATE_KHR) has_immediate = true;
        }
        sf->present_mode = VK_PRESENT_MODE_FIFO_KHR;
        if (desc->present_mode == MD_GPU_PRESENT_IMMEDIATE) {
            sf->present_mode = has_immediate ? VK_PRESENT_MODE_IMMEDIATE_KHR : (has_mailbox ? VK_PRESENT_MODE_MAILBOX_KHR : VK_PRESENT_MODE_FIFO_KHR);
        } else if (desc->present_mode == MD_GPU_PRESENT_MAILBOX) {
            sf->present_mode = has_mailbox ? VK_PRESENT_MODE_MAILBOX_KHR : VK_PRESENT_MODE_FIFO_KHR;
        }
    }

    for (uint32_t i = 0; i < MD_VK_ACQUIRE_SEMS; ++i) {
        VkSemaphoreCreateInfo semci = {VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        if (!md_vk_check(vkCreateSemaphore(dev->device, &semci, NULL, &sf->acquire_sems[i]), "vkCreateSemaphore")) goto fail;
    }

    /* Build now when there is a size, so that errors surface here. */
    if (sf->want_width && sf->want_height) {
        bool zero;
        if (!md_vk_surface_rebuild(sf, &zero)) goto fail;
    }

    md_mutex_lock(&dev->device_mutex);
    const bool listed = MD_ARRAY_TRY_PUSH(dev->surfaces, sf);
    md_mutex_unlock(&dev->device_mutex);
    if (!listed) { md_gpu_fail("out of memory"); goto fail; }
    return sf;

fail:
    md_mutex_lock(&dev->device_mutex);
    if (sf->sc) { md_vk_swapchain_free(dev, sf->sc); sf->sc = NULL; }
    md_mutex_unlock(&dev->device_mutex);
    for (uint32_t i = 0; i < MD_VK_ACQUIRE_SEMS; ++i) if (sf->acquire_sems[i]) vkDestroySemaphore(dev->device, sf->acquire_sems[i], NULL);
    vkDestroySurfaceKHR(dev->instance, vs, NULL);
    md_free(dev->alloc, sf, sizeof(*sf));
    return NULL;
}

void md_gpu_surface_destroy(md_gpu_surface_t sf) {
    if (!sf) return;
    md_gpu_device_t dev = sf->device;
    md_mutex_lock(&dev->device_mutex);
    MD_ARRAY_REMOVE_ITEM(dev->surfaces, sf);
    md_vk_swapchain_t* sc = sf->sc;
    if (!sc) {
        sc = (md_vk_swapchain_t*)md_alloc(dev->alloc, sizeof(md_vk_swapchain_t));
        if (sc) memset(sc, 0, sizeof(*sc));
    }
    VkSemaphore* sems = (VkSemaphore*)md_alloc(dev->alloc, MD_VK_ACQUIRE_SEMS * sizeof(VkSemaphore));
    if (!sc || !sems) {
        /* Cannot defer: idle and free now. */
        md_gpu_fail("out of memory destroying surface '%s'; waiting for the device", sf->label);
        vkDeviceWaitIdle(dev->device);
        if (sems) md_free(dev->alloc, sems, MD_VK_ACQUIRE_SEMS * sizeof(VkSemaphore));
        if (sc) { sc->present_queue = NULL; md_vk_swapchain_free(dev, sc); }
        for (uint32_t i = 0; i < MD_VK_ACQUIRE_SEMS; ++i) vkDestroySemaphore(dev->device, sf->acquire_sems[i], NULL);
        vkDestroySurfaceKHR(dev->instance, sf->surface, NULL);
    } else {
        memcpy(sems, sf->acquire_sems, MD_VK_ACQUIRE_SEMS * sizeof(VkSemaphore));
        sc->extra_sems      = sems;
        sc->extra_sem_count = MD_VK_ACQUIRE_SEMS;
        sc->surface         = sf->surface;
        if (!sc->present_queue) sc->present_queue = sf->present_queue;
        md_vk_retire_locked(dev, MD_VK_RETIRE_SWAPCHAIN, sc);
    }
    md_mutex_unlock(&dev->device_mutex);
    md_free(dev->alloc, sf, sizeof(*sf));
}

void md_gpu_surface_resize(md_gpu_surface_t sf, uint32_t width, uint32_t height) {
    if (!sf) return;
    md_gpu_device_t dev = sf->device;
    md_mutex_lock(&dev->device_mutex);
    sf->want_width  = width;
    sf->want_height = height;
    if (!sf->sc || sf->sc->width != width || sf->sc->height != height) sf->dirty = true;
    md_mutex_unlock(&dev->device_mutex);
}

md_gpu_texture_t md_gpu_surface_acquire(md_gpu_stream_t s, md_gpu_surface_t sf) {
    if (!s || !sf) { md_gpu_fail("md_gpu_surface_acquire: null argument"); return NULL; }
    md_gpu_device_t dev = s->device;
    if (sf->device != dev) { md_gpu_fail("md_gpu_surface_acquire: surface '%s' belongs to another device", sf->label); return NULL; }
    if (!s->can_graphics) { md_gpu_fail("md_gpu_surface_acquire: stream '%s' is not a GRAPHICS stream", s->label); return NULL; }
    if (!md_vk_not_in_pass(s, "md_gpu_surface_acquire")) return NULL;
    if (s->upload_open) { md_gpu_fail("md_gpu_surface_acquire: stream '%s' has an open upload", s->label); return NULL; }
    if (sf->acquired) { md_gpu_fail("md_gpu_surface_acquire: surface '%s' already has an acquired image; present it first", sf->label); return NULL; }
    if (s->bin_wait_count >= 4) { md_gpu_fail("md_gpu_surface_acquire: too many acquires pending on stream '%s'", s->label); return NULL; }

    for (int attempt = 0; attempt < 3; ++attempt) {
        if (!sf->sc || sf->dirty) {
            bool zero;
            if (!md_vk_surface_rebuild(sf, &zero)) return NULL;
            if (zero) return NULL;
        }
        const uint32_t slot = sf->acquire_next;
        if (sf->acquire_stream[slot] && sf->acquire_value[slot] &&
            md_vk_stream_completed(sf->acquire_stream[slot]) < sf->acquire_value[slot]) {
            /* The semaphore's last wait has not executed yet; a stream must
               have been left unsubmitted for many frames. Submit it and wait. */
            md_gpu_stream_t o = sf->acquire_stream[slot];
            if (o->submitted_value < sf->acquire_value[slot]) md_vk_stream_submit(o);
            VkSemaphoreWaitInfo wi = {VK_STRUCTURE_TYPE_SEMAPHORE_WAIT_INFO};
            wi.semaphoreCount = 1;
            wi.pSemaphores    = &o->timeline;
            wi.pValues        = &sf->acquire_value[slot];
            vkWaitSemaphores(dev->device, &wi, UINT64_MAX);
        }
        const VkSemaphore sem = sf->acquire_sems[slot];
        uint32_t index = 0;
        VkResult r = vkAcquireNextImageKHR(dev->device, sf->sc->swapchain, UINT64_MAX, sem, VK_NULL_HANDLE, &index);
        if (r == VK_ERROR_OUT_OF_DATE_KHR) { sf->dirty = true; continue; }
        if (r == VK_SUBOPTIMAL_KHR) {
            sf->dirty = true;       /* usable now; rebuilt at the next acquire */
        } else if (r != VK_SUCCESS) {
            md_vk_check(r, "vkAcquireNextImageKHR");
            return NULL;
        }
        sf->acquire_next = (slot + 1) % MD_VK_ACQUIRE_SEMS;

        /* Work already recorded must not wait for the image. */
        if (s->has_work) md_vk_stream_submit(s);
        md_gpu_texture_t tex = sf->sc->textures[index];
        if (!md_vk_record_image_barrier(s, tex->image, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL,
                                        VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0,
                                        VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT,
                                        VK_ACCESS_2_MEMORY_READ_BIT | VK_ACCESS_2_MEMORY_WRITE_BIT)) {
            return NULL;
        }
        /* The next submission waits on the acquire; it signals next_value,
           after which the semaphore is free again. */
        s->bin_waits[s->bin_wait_count++] = sem;
        sf->acquire_stream[slot] = s;
        sf->acquire_value[slot]  = s->next_value;
        sf->acquired        = true;
        sf->image_index     = index;
        sf->acquire_slot    = slot;
        sf->acquired_stream = s;
        return tex;
    }
    md_gpu_fail("md_gpu_surface_acquire: surface '%s' stayed out of date", sf->label);
    return NULL;
}

bool md_gpu_surface_present(md_gpu_stream_t s, md_gpu_surface_t sf) {
    if (!s || !sf) return md_gpu_fail("md_gpu_surface_present: null argument");
    if (!sf->acquired) return md_gpu_fail("md_gpu_surface_present: surface '%s' has no acquired image", sf->label);
    if (sf->acquired_stream != s) return md_gpu_fail("md_gpu_surface_present: surface '%s' was acquired on another stream", sf->label);
    if (!md_vk_not_in_pass(s, "md_gpu_surface_present")) return false;
    if (s->upload_open) return md_gpu_fail("md_gpu_surface_present: stream '%s' has an open upload", s->label);
    if (s->bin_signal_count >= 4) return md_gpu_fail("md_gpu_surface_present: too many presents pending on stream '%s'", s->label);

    md_gpu_device_t dev = s->device;
    md_vk_swapchain_t* sc = sf->sc;
    const uint32_t index = sf->image_index;
    sf->acquired = false;
    sf->acquired_stream = NULL;
    if (!md_vk_record_image_barrier(s, sc->textures[index]->image, VK_IMAGE_ASPECT_COLOR_BIT, VK_IMAGE_LAYOUT_GENERAL, VK_IMAGE_LAYOUT_PRESENT_SRC_KHR,
                                    VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, VK_ACCESS_2_MEMORY_WRITE_BIT,
                                    /* Chains into the semaphore signal below,
                                       whose stage is ALL_COMMANDS, so the
                                       transition is inside what present waits
                                       for. */
                                    VK_PIPELINE_STAGE_2_ALL_COMMANDS_BIT, 0)) {
        return false;
    }
    s->bin_signals[s->bin_signal_count++] = sc->present_sems[index];
    if (!md_vk_stream_submit(s)) return false;

    VkPresentInfoKHR pi = {VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    pi.waitSemaphoreCount = 1;
    pi.pWaitSemaphores    = &sc->present_sems[index];
    pi.swapchainCount     = 1;
    pi.pSwapchains        = &sc->swapchain;
    pi.pImageIndices      = &index;
    md_mutex_lock(&dev->queue_mutex);
    VkResult r = vkQueuePresentKHR(s->queue, &pi);
    md_mutex_unlock(&dev->queue_mutex);
    sf->present_queue = s->queue;
    sc->present_queue = s->queue;
    if (r == VK_SUBOPTIMAL_KHR || r == VK_ERROR_OUT_OF_DATE_KHR) {
        sf->dirty = true;
        return true;
    }
    return md_vk_check(r, "vkQueuePresentKHR");
}

/* =========================================================================
   14. Host callbacks and polling
   ========================================================================= */

bool md_gpu_sync_on_complete(md_gpu_device_t dev, md_gpu_sync_t sync, md_gpu_host_fn fn, void* user) {
    if (!dev || !fn) return md_gpu_fail("md_gpu_sync_on_complete: null argument");
    if (md_gpu_sync_is_valid(sync) && sync.stream->device != dev) return md_gpu_fail("md_gpu_sync_on_complete: sync from another device");
    md_mutex_lock(&dev->device_mutex);
    const md_vk_hostfn_t h = {sync, fn, user};
    const bool listed = MD_ARRAY_TRY_PUSH(dev->hostfns, h);
    md_mutex_unlock(&dev->device_mutex);
    return listed ? true : md_gpu_fail("out of memory");
}

bool md_gpu_launch_host_fn(md_gpu_stream_t s, md_gpu_host_fn fn, void* user) {
    if (!s || !fn) return md_gpu_fail("md_gpu_launch_host_fn: null argument");
    if (!md_vk_not_in_pass(s, "md_gpu_launch_host_fn")) return false;
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
MD_ARRAY_RAW_DEFINE_TYPE(md_vk_snapshot_array_t, md_vk_snapshot_t);

static bool md_vk_snapshot_complete(md_gpu_device_t dev, md_vk_snapshot_array_t* snap, md_gpu_sync_t sync) {
    if (!md_gpu_sync_is_valid(sync)) return true;
    for (size_t i = 0; i < snap->size; ++i) {
        md_vk_snapshot_t* e = &snap->data[i];
        if (e->stream == sync.stream) return e->completed >= sync.value;
    }
    uint64_t completed = md_vk_stream_completed(sync.stream);
    const md_vk_snapshot_t e = {sync.stream, completed};
    (void)MD_ARRAY_TRY_PUSH(*snap, e);
    /* If the push failed, readiness falls back to this instant's reading; the
       ordering guarantee then degrades, but nothing fires early. */
    return completed >= sync.value;
}

uint32_t md_gpu_device_poll(md_gpu_device_t dev) {
    if (!dev) return 0;
    uint32_t fired = 0;
    md_vk_snapshot_array_t snap;
    MD_ARRAY_INIT(snap, dev->alloc);

    /* Callbacks run outside the lock, one at a time, in registration order.
       Each may register new callbacks or free memory; those are seen by the
       loop but judged against the same snapshot. */
    for (;;) {
        md_vk_hostfn_t ready;
        bool have = false;
        md_mutex_lock(&dev->device_mutex);
        for (size_t i = 0; i < dev->hostfns.size; ++i) {
            md_vk_hostfn_t* h = &dev->hostfns.data[i];
            if (!md_vk_snapshot_complete(dev, &snap, h->sync)) continue;
            ready = *h;
            MD_ARRAY_REMOVE(dev->hostfns, i);
            have = true;
            break;
        }
        md_mutex_unlock(&dev->device_mutex);
        if (!have) break;
        ready.fn(ready.user);
        fired++;
    }
    MD_ARRAY_FREE(snap);

    md_mutex_lock(&dev->device_mutex);
    md_vk_process_retires_locked(dev, false);
    md_gpu_heap_poll_locked(&dev->heap);
    md_mutex_unlock(&dev->device_mutex);
    return fired;
}

/* =========================================================================
   15. Device destruction
   ========================================================================= */

void md_gpu_device_destroy(md_gpu_device_t dev) {
    if (!dev) return;
    struct md_allocator_i* alloc = dev->alloc;
    if (dev->device) {
        /* Flush and idle every stream, then fire what is pending. */
        for (size_t i = 0; i < dev->streams.size; ++i) {
            md_gpu_stream_t s = dev->streams.data[i];
            s->upload_open = false;
            md_vk_abandon_pass(s);
            md_vk_stream_submit(s);
        }
        vkDeviceWaitIdle(dev->device);
        md_gpu_device_poll(dev);

        /* Everything the caller did not destroy, the device does. */
        while (dev->textures.size > 0) md_gpu_texture_destroy(dev->textures.data[0]);
        while (dev->kernels.size > 0) md_gpu_kernel_destroy(dev->kernels.data[0]);
        while (dev->pipelines.size > 0) md_gpu_pipeline_destroy(dev->pipelines.data[0]);
        while (dev->surfaces.size > 0) md_gpu_surface_destroy(dev->surfaces.data[0]);
        if (dev->make_grid_kernel) md_vk_kernel_free(dev, dev->make_grid_kernel);

        md_mutex_lock(&dev->device_mutex);
        md_vk_process_retires_locked(dev, true);
        md_mutex_unlock(&dev->device_mutex);

        for (size_t i = 0; i < dev->streams.size; ++i) {
            md_vk_stream_free(dev->streams.data[i]);
        }
        MD_ARRAY_FREE(dev->streams);
        md_gpu_heap_destroy(&dev->heap);
        MD_ARRAY_FREE(dev->textures);
        MD_ARRAY_FREE(dev->kernels);
        MD_ARRAY_FREE(dev->pipelines);
        MD_ARRAY_FREE(dev->surfaces);

        for (uint32_t i = 0; i < dev->sampler_count; ++i) vkDestroySampler(dev->device, dev->samplers[i].sampler, NULL);
        if (dev->dummy_view)    vkDestroyImageView(dev->device, dev->dummy_view, NULL);
        if (dev->dummy_image)   vkDestroyImage(dev->device, dev->dummy_image, NULL);
        if (dev->dummy_mem)     vkFreeMemory(dev->device, dev->dummy_mem, NULL);
        if (dev->dummy_sampler) vkDestroySampler(dev->device, dev->dummy_sampler, NULL);
        MD_ARRAY_FREE(dev->hostfns);
        MD_ARRAY_FREE(dev->retires);
        if (dev->pipeline_layout) vkDestroyPipelineLayout(dev->device, dev->pipeline_layout, NULL);
        if (dev->raster_layout)   vkDestroyPipelineLayout(dev->device, dev->raster_layout, NULL);
        if (dev->desc_pool)       vkDestroyDescriptorPool(dev->device, dev->desc_pool, NULL);
        if (dev->set_layout)      vkDestroyDescriptorSetLayout(dev->device, dev->set_layout, NULL);
        if (dev->tex_free)        md_free(alloc, dev->tex_free, (dev->tex_slot_count - 1) * sizeof(uint32_t));
        md_mutex_destroy(&dev->queue_mutex);
        md_mutex_destroy(&dev->device_mutex);
        vkDestroyDevice(dev->device, NULL);
    }
    if (dev->messenger && vkDestroyDebugUtilsMessengerEXT) {
        vkDestroyDebugUtilsMessengerEXT(dev->instance, dev->messenger, NULL);
    }
    if (dev->instance && dev->instance == md_vk_live_instance) md_vk_live_instance = VK_NULL_HANDLE;
    if (dev->instance) vkDestroyInstance(dev->instance, NULL);
    md_free(alloc, dev, sizeof(*dev));
}

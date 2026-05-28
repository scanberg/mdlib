#include <md_topo.h>

#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_grid.h>
#include <core/md_array.h>

#include <float.h>

#if DEBUG
#include <core/md_hash.h>
#include <stdlib.h>

// Comparison function to sort uint32 indices (for quicksort C)
static int uint_compare(const void* a, const void* b) {
    uint32_t val_a = *(const uint32_t*)a;
    uint32_t val_b = *(const uint32_t*)b;
    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0;
}
#endif

static inline void index_to_world_matrix(float out_mat[4][4], const md_grid_t* grid) {
    out_mat[0][0] = grid->orientation.elem[0][0] * grid->spacing.elem[0];
    out_mat[0][1] = grid->orientation.elem[0][1] * grid->spacing.elem[0];
    out_mat[0][2] = grid->orientation.elem[0][2] * grid->spacing.elem[0];
    out_mat[0][3] = 0.0f;
    out_mat[1][0] = grid->orientation.elem[1][0] * grid->spacing.elem[1];
    out_mat[1][1] = grid->orientation.elem[1][1] * grid->spacing.elem[1];
    out_mat[1][2] = grid->orientation.elem[1][2] * grid->spacing.elem[1];
    out_mat[1][3] = 0.0f;
    out_mat[2][0] = grid->orientation.elem[2][0] * grid->spacing.elem[2];
    out_mat[2][1] = grid->orientation.elem[2][1] * grid->spacing.elem[2];
    out_mat[2][2] = grid->orientation.elem[2][2] * grid->spacing.elem[2];
    out_mat[2][3] = 0.0f;
    out_mat[3][0] = grid->origin.elem[0];
    out_mat[3][1] = grid->origin.elem[1];
    out_mat[3][2] = grid->origin.elem[2];
    out_mat[3][3] = 1.0f;

    // Incorporate a half voxel offset to move to voxel centers
    out_mat[3][0] += 0.5f * (out_mat[0][0] + out_mat[1][0] + out_mat[2][0]);
    out_mat[3][1] += 0.5f * (out_mat[0][1] + out_mat[1][1] + out_mat[2][1]);
    out_mat[3][2] += 0.5f * (out_mat[0][2] + out_mat[1][2] + out_mat[2][2]);
}

#if MD_ENABLE_GPU

#include <core/md_gpu.h>
#include <topo_gpu_shaders.inl>
#include <topo_gpu_shaders_reflection.inl>

// Pipeline cache
static md_gpu_compute_pipeline_t pip_bidirectional_manifold    = NULL;
static md_gpu_compute_pipeline_t pip_path_compression          = NULL;
static md_gpu_compute_pipeline_t pip_critical_points           = NULL;
static md_gpu_compute_pipeline_t pip_critical_point_compaction = NULL;
static md_gpu_compute_pipeline_t pip_vertex_edge_extraction    = NULL;

static md_gpu_compute_pipeline_t ensure_pipeline(md_gpu_device_t device, md_gpu_compute_pipeline_t* pipeline, const void* blob_start, size_t blob_size, const char* name, uint32_t wg_x, uint32_t wg_y, uint32_t wg_z, const md_gpu_resource_binding_t* resource_bindings, uint32_t resource_binding_count) {
    if (*pipeline == NULL) {
        md_gpu_compute_pipeline_desc_t desc = {
            .shader_bytes     = blob_start,
            .shader_byte_size = blob_size,
            .threadgroup_size = { wg_x, wg_y, wg_z },
            .resource_bindings = resource_bindings,
            .resource_binding_count = resource_binding_count,
        };
        *pipeline = md_gpu_compute_pipeline_create(device, &desc);
        if (*pipeline == NULL) {
            MD_LOG_ERROR("Failed to create compute pipeline: %s", name);
        }
    }
    return *pipeline;
}

void md_topo_gpu_initialize(md_gpu_device_t device) {
    if (!ensure_pipeline(device, &pip_bidirectional_manifold,    topo_bidirectional_manifold_start,    topo_bidirectional_manifold_size(),    "bidirectional_manifold",    topo_bidirectional_manifold_thread_group_size_x,    topo_bidirectional_manifold_thread_group_size_y,    topo_bidirectional_manifold_thread_group_size_z,    topo_bidirectional_manifold_pipeline_resource_bindings, topo_bidirectional_manifold_pipeline_resource_binding_count) ||
        !ensure_pipeline(device, &pip_path_compression,          topo_path_compression_start,          topo_path_compression_size(),          "path_compression",          topo_path_compression_thread_group_size_x,          topo_path_compression_thread_group_size_y,          topo_path_compression_thread_group_size_z,          NULL, topo_path_compression_pipeline_resource_binding_count) ||
        !ensure_pipeline(device, &pip_critical_points,           topo_critical_points_start,           topo_critical_points_size(),           "critical_points",           topo_critical_points_thread_group_size_x,           topo_critical_points_thread_group_size_y,           topo_critical_points_thread_group_size_z,           topo_critical_points_pipeline_resource_bindings, topo_critical_points_pipeline_resource_binding_count) ||
        !ensure_pipeline(device, &pip_critical_point_compaction, topo_critical_point_compaction_start, topo_critical_point_compaction_size(), "critical_point_compaction", topo_critical_point_compaction_thread_group_size_x, topo_critical_point_compaction_thread_group_size_y, topo_critical_point_compaction_thread_group_size_z, NULL, topo_critical_point_compaction_pipeline_resource_binding_count) ||
        !ensure_pipeline(device, &pip_vertex_edge_extraction,    topo_vertex_edge_extraction_start,    topo_vertex_edge_extraction_size(),    "vertex_edge_extraction",    topo_vertex_edge_extraction_thread_group_size_x,    topo_vertex_edge_extraction_thread_group_size_y,    topo_vertex_edge_extraction_thread_group_size_z,    topo_vertex_edge_extraction_pipeline_resource_bindings, topo_vertex_edge_extraction_pipeline_resource_binding_count))
    {
        MD_LOG_ERROR("md_topo_gpu_initialize: failed to create compute pipelines");
    }
}

void md_topo_gpu_shutdown(void) {
    if (pip_bidirectional_manifold)    md_gpu_compute_pipeline_destroy(pip_bidirectional_manifold);
    if (pip_path_compression)          md_gpu_compute_pipeline_destroy(pip_path_compression);
    if (pip_critical_points)           md_gpu_compute_pipeline_destroy(pip_critical_points);
    if (pip_critical_point_compaction) md_gpu_compute_pipeline_destroy(pip_critical_point_compaction);
    if (pip_vertex_edge_extraction)    md_gpu_compute_pipeline_destroy(pip_vertex_edge_extraction);
    pip_bidirectional_manifold = NULL;
    pip_path_compression = NULL;
    pip_critical_points = NULL;
    pip_critical_point_compaction = NULL;
    pip_vertex_edge_extraction = NULL;
}

// Meta buffer: bound as a single SSBO to all compute shaders.
typedef struct {
    uint32_t vertex_count;  // total CP count (critical_points output)
    uint32_t edge_count;    // actual edge count (extraction output)
    uint32_t changed_read;  // path-compression convergence flag, read by shader
    uint32_t changed_write; // path-compression convergence flag, written by shader
    uint32_t counter;       // compaction write cursor
} topo_meta_t;

// Worst-case capacity ratios:
//   vert_cap = num_points / TOPO_VERT_RATIO  (1 CP per 8 voxels is very generous)
//   edge_cap = vert_cap * TOPO_EDGE_RATIO
#define TOPO_VERT_RATIO  8
#define TOPO_EDGE_RATIO  4
#define TOPO_VERT_CAP_MIN 64

struct md_topo_gpu_context {
    md_gpu_device_t device;
    uint32_t        num_points;
    uint32_t        dim[3];
    uint32_t        vert_cap;      // pre-allocated worst-case vertex capacity
    uint32_t        edge_cap;      // pre-allocated worst-case edge capacity

    // Persistent GPU buffers — scratch (each u32[num_points], device-local)
    md_gpu_buffer_t ascending_buf;     // steepest-ascent parent indices
    md_gpu_buffer_t descending_buf;    // steepest-descent parent indices
    md_gpu_buffer_t voxel_types_buf;   // per-voxel critical-point type scratch
    md_gpu_buffer_t voxel_to_vert_buf; // voxel index → vertex index map
    md_gpu_buffer_t meta_buf;          // sizeof(topo_meta_t) bytes, device-local
    
    // Persistent GPU buffers — results (pre-allocated to worst-case capacity)
    md_gpu_buffer_t indices_buf;   // u32[vert_cap]
    md_gpu_buffer_t vert_buf;      // float4[vert_cap]
    md_gpu_buffer_t type_buf;      // u32[vert_cap]: per-vertex CP type
    md_gpu_buffer_t edge_buf;      // md_topo_edge_t[edge_cap]
    
    // These are optional staging buffers
    md_gpu_buffer_t staging_buf;   // sizeof(topo_meta_t) bytes, CPU-visible
    md_gpu_buffer_t staging_verts; // float4[vert_cap], CPU-visible
    md_gpu_buffer_t staging_types; // u32[vert_cap], CPU-visible
    md_gpu_buffer_t staging_edges; // md_topo_edge_t[edge_cap], CPU-visible
};

md_topo_gpu_context_t* md_topo_gpu_context_create(md_gpu_device_t device, uint32_t dim_x, uint32_t dim_y, uint32_t dim_z) {
    if (!device || !dim_x || !dim_y || !dim_z) {
        MD_LOG_ERROR("md_topo_gpu_context_create: invalid arguments");
        return NULL;
    }

    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;

    md_gpu_device_info_t info = {0};
    md_gpu_device_info(device, &info);

    md_gpu_buffer_flags_t buf_flags = MD_GPU_BUFFER_NONE;
    if (!info.is_discrete) {
        // For integrated GPUs, we can skip staging buffers and write results directly to CPU-visible buffers
        buf_flags |= MD_GPU_BUFFER_CPU_VISIBLE;
    }

    ctx->device     = device;
    ctx->num_points = dim_x * dim_y * dim_z;
    ctx->dim[0]     = dim_x;
    ctx->dim[1]     = dim_y;
    ctx->dim[2]     = dim_z;

    uint32_t vert_cap = MAX(ctx->num_points / TOPO_VERT_RATIO, TOPO_VERT_CAP_MIN);
    uint32_t edge_cap = vert_cap * TOPO_EDGE_RATIO;
    ctx->vert_cap = vert_cap;
    ctx->edge_cap = edge_cap;

    const size_t voxel_buf_size = ctx->num_points * sizeof(uint32_t);
    ctx->ascending_buf      = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = voxel_buf_size,     .flags = buf_flags });
    ctx->descending_buf     = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = voxel_buf_size,     .flags = buf_flags });
    ctx->voxel_types_buf    = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = voxel_buf_size,     .flags = buf_flags });
    ctx->voxel_to_vert_buf  = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = voxel_buf_size,     .flags = buf_flags });
    ctx->meta_buf           = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = sizeof(topo_meta_t), .flags = buf_flags });
    ctx->indices_buf        = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * sizeof(uint32_t),        .flags = buf_flags });
    
    ctx->vert_buf           = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * 4 * sizeof(float),       .flags = buf_flags });
    ctx->type_buf           = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * sizeof(uint32_t),        .flags = buf_flags });
    ctx->edge_buf           = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = edge_cap * sizeof(md_topo_edge_t),  .flags = buf_flags });
    
    if (info.is_discrete) {
        // For discrete GPUs, we need staging buffers for readback
        ctx->staging_buf    = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = sizeof(topo_meta_t),                .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        ctx->staging_verts  = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * 4 * sizeof(float),       .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        ctx->staging_types  = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * sizeof(uint32_t),        .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        ctx->staging_edges  = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = edge_cap * 2 * sizeof(uint32_t),    .flags = MD_GPU_BUFFER_CPU_VISIBLE });
    }

    bool success = ctx->ascending_buf && ctx->descending_buf && ctx->voxel_types_buf && ctx->voxel_to_vert_buf && ctx->meta_buf && ctx->indices_buf && ctx->vert_buf && ctx->type_buf && ctx->edge_buf;
    if (info.is_discrete) {
        success = success && ctx->staging_buf && ctx->staging_verts && ctx->staging_types && ctx->staging_edges;
    }

    if (!success) {
        MD_LOG_ERROR("md_topo_gpu_context_create: failed to allocate GPU buffers");
        md_topo_gpu_context_destroy((md_topo_gpu_context_t*)ctx);
        return NULL;
    }

    MD_LOG_DEBUG("md_topo_gpu_context_create: %ux%ux%u (%u voxels), vert_cap=%u edge_cap=%u",
                 dim_x, dim_y, dim_z, ctx->num_points, vert_cap, edge_cap);
    return (md_topo_gpu_context_t*)ctx;
}

void md_topo_gpu_context_destroy(md_topo_gpu_context_t* context) {
    if (!context) return;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;
    if (ctx->ascending_buf)     md_gpu_buffer_destroy(ctx->ascending_buf);
    if (ctx->descending_buf)    md_gpu_buffer_destroy(ctx->descending_buf);
    if (ctx->voxel_types_buf)   md_gpu_buffer_destroy(ctx->voxel_types_buf);
    if (ctx->voxel_to_vert_buf) md_gpu_buffer_destroy(ctx->voxel_to_vert_buf);
    if (ctx->meta_buf)          md_gpu_buffer_destroy(ctx->meta_buf);
    if (ctx->staging_buf)   md_gpu_buffer_destroy(ctx->staging_buf);
    if (ctx->indices_buf)   md_gpu_buffer_destroy(ctx->indices_buf);
    if (ctx->vert_buf)      md_gpu_buffer_destroy(ctx->vert_buf);
    if (ctx->type_buf)      md_gpu_buffer_destroy(ctx->type_buf);
    if (ctx->edge_buf)      md_gpu_buffer_destroy(ctx->edge_buf);
    if (ctx->staging_types) md_gpu_buffer_destroy(ctx->staging_types);
    if (ctx->staging_verts) md_gpu_buffer_destroy(ctx->staging_verts);
    if (ctx->staging_edges) md_gpu_buffer_destroy(ctx->staging_edges);
    free(ctx);
}

void md_topo_gpu_record(md_gpu_cmd_t cmd, md_topo_gpu_context_t* context, md_gpu_image_t volume, const md_grid_t* grid, float scalar_threshold) {
    if (!cmd || !context || !volume || !grid) {
        MD_LOG_ERROR("md_topo_gpu_record: invalid input");
        return;
    }
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    if (!pip_bidirectional_manifold || !pip_path_compression || !pip_critical_points || !pip_critical_point_compaction || !pip_vertex_edge_extraction) {
        MD_LOG_ERROR("md_topo_gpu_record: compute pipelines not initialized");
        return;
    }

#define TOPO_FILL_COMMON_ARGS(root_args) \
    do { \
        index_to_world_matrix((root_args).index_to_world, grid); \
        (root_args).dims[0] = grid->dim[0]; \
        (root_args).dims[1] = grid->dim[1]; \
        (root_args).dims[2] = grid->dim[2]; \
        (root_args).dims[3] = 0; \
        (root_args).scalar_threshold = scalar_threshold; \
    } while (0)

    const uint32_t wg[3]         = { DIV_UP(grid->dim[0], 8), DIV_UP(grid->dim[1], 8), DIV_UP(grid->dim[2], 8) };
    const size_t   voxel_buf_size = ctx->num_points * sizeof(uint32_t);

    // Per-call resets (voxel_to_vert_idx = -1, count = 0, counter = 0, edge_count = 0, vertex_types = 0)
    md_gpu_cmd_push_debug_group(cmd, "Per-call resets");
    md_gpu_cmd_fill_buffer(cmd, ctx->voxel_to_vert_buf, 0, voxel_buf_size, 0xFF); // voxel_to_vert_idx = -1
    md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf, 0, sizeof(topo_meta_t), 0);
    md_gpu_cmd_fill_buffer(cmd, ctx->type_buf, 0, ctx->vert_cap * sizeof(uint32_t), 0); // vertex_types = 0
    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);
    md_gpu_cmd_pop_debug_group(cmd);

    // Step 1: Bidirectional manifold
    md_gpu_cmd_push_debug_group(cmd, "Bidirectional manifold");
    topo_bidirectional_manifold_dispatch_t bidirectional_dispatch = topo_bidirectional_manifold_dispatch_init();
    TOPO_FILL_COMMON_ARGS(bidirectional_dispatch.args);
    bidirectional_dispatch.resources.ascending = (md_gpu_buffer_resource_t){ .buffer = ctx->ascending_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    bidirectional_dispatch.resources.descending = (md_gpu_buffer_resource_t){ .buffer = ctx->descending_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    bidirectional_dispatch.resources.volumeTex.image = volume;
    bidirectional_dispatch.resources.volumeTex.usage = GPU_USAGE_READ;
    bidirectional_dispatch.group_count[0] = wg[0];
    bidirectional_dispatch.group_count[1] = wg[1];
    bidirectional_dispatch.group_count[2] = wg[2];
    topo_bidirectional_manifold_cmd_dispatch(cmd, pip_bidirectional_manifold, &bidirectional_dispatch);
    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);
    md_gpu_cmd_pop_debug_group(cmd);

    // Step 2: Path compression (iterative, GPU-side early-exit via changed flag)
    // Each dispatch does one grandparent pointer-jump. A path of length L needs
    // ceil(log2(L)) iterations; 2 * ceil(log2(max_dim)) is a safe upper bound.
    uint32_t num_iterations = 0;
    {
        uint32_t max_dim = (uint32_t)MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
        while (max_dim > (1U << num_iterations)) num_iterations++;
        num_iterations *= 2;
    }

    md_gpu_cmd_push_debug_group(cmd, "Path compression");
    md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf, offsetof(topo_meta_t, changed_read), 4, 0xFF);
    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);

    topo_path_compression_dispatch_t path_dispatch = topo_path_compression_dispatch_init();
    TOPO_FILL_COMMON_ARGS(path_dispatch.args);
    path_dispatch.resources.ascending = (md_gpu_buffer_resource_t){ .buffer = ctx->ascending_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    path_dispatch.resources.descending = (md_gpu_buffer_resource_t){ .buffer = ctx->descending_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    path_dispatch.resources.meta = (md_gpu_buffer_resource_t){ .buffer = ctx->meta_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    path_dispatch.group_count[0] = wg[0];
    path_dispatch.group_count[1] = wg[1];
    path_dispatch.group_count[2] = wg[2];
    for (uint32_t i = 0; i < num_iterations; i++) {
        topo_path_compression_cmd_dispatch(cmd, pip_path_compression, &path_dispatch);
        md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf, ctx->meta_buf, 4, offsetof(topo_meta_t, changed_write), offsetof(topo_meta_t, changed_read));
        md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf, offsetof(topo_meta_t, changed_write), 4, 0);
        md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);
    }
    md_gpu_cmd_pop_debug_group(cmd);

    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);

    // Step 3: Critical-point detection → meta_buf.count
    md_gpu_cmd_push_debug_group(cmd, "Critical-point detection");
    topo_critical_points_dispatch_t critical_dispatch = topo_critical_points_dispatch_init();
    TOPO_FILL_COMMON_ARGS(critical_dispatch.args);
    critical_dispatch.resources.ascending = (md_gpu_buffer_resource_t){ .buffer = ctx->ascending_buf, .offset = 0, .usage = GPU_USAGE_READ };
    critical_dispatch.resources.descending = (md_gpu_buffer_resource_t){ .buffer = ctx->descending_buf, .offset = 0, .usage = GPU_USAGE_READ };
    critical_dispatch.resources.types = (md_gpu_buffer_resource_t){ .buffer = ctx->voxel_types_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    critical_dispatch.resources.meta = (md_gpu_buffer_resource_t){ .buffer = ctx->meta_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    critical_dispatch.resources.volumeTex.image = volume;
    critical_dispatch.resources.volumeTex.usage = GPU_USAGE_READ;
    critical_dispatch.group_count[0] = wg[0];
    critical_dispatch.group_count[1] = wg[1];
    critical_dispatch.group_count[2] = wg[2];
    topo_critical_points_cmd_dispatch(cmd, pip_critical_points, &critical_dispatch);
    md_gpu_cmd_pop_debug_group(cmd);

    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);

    // Step 4: Critical-point compaction (single atomic counter)
    md_gpu_cmd_push_debug_group(cmd, "Critical-point compaction");
    topo_critical_point_compaction_dispatch_t compaction_dispatch = topo_critical_point_compaction_dispatch_init();
    TOPO_FILL_COMMON_ARGS(compaction_dispatch.args);
    compaction_dispatch.resources.types = (md_gpu_buffer_resource_t){ .buffer = ctx->voxel_types_buf, .offset = 0, .usage = GPU_USAGE_READ };
    compaction_dispatch.resources.cp_indices = (md_gpu_buffer_resource_t){ .buffer = ctx->indices_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    compaction_dispatch.resources.meta = (md_gpu_buffer_resource_t){ .buffer = ctx->meta_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    compaction_dispatch.resources.voxel_to_vertex_idx = (md_gpu_buffer_resource_t){ .buffer = ctx->voxel_to_vert_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    compaction_dispatch.resources.vertex_types = (md_gpu_buffer_resource_t){ .buffer = ctx->type_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    compaction_dispatch.group_count[0] = wg[0];
    compaction_dispatch.group_count[1] = wg[1];
    compaction_dispatch.group_count[2] = wg[2];
    topo_critical_point_compaction_cmd_dispatch(cmd, pip_critical_point_compaction, &compaction_dispatch);
    md_gpu_cmd_pop_debug_group(cmd);

    md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);

    // Step 5: Vertex + edge extraction
    // Dispatched at worst-case capacity; unused slots have vertex_types == 0 → shader returns early.
    md_gpu_cmd_push_debug_group(cmd, "Vertex + edge extraction");
    topo_vertex_edge_extraction_dispatch_t extraction_dispatch = topo_vertex_edge_extraction_dispatch_init();
    TOPO_FILL_COMMON_ARGS(extraction_dispatch.args);
    extraction_dispatch.resources.cp_indices = (md_gpu_buffer_resource_t){ .buffer = ctx->indices_buf, .offset = 0, .usage = GPU_USAGE_READ };
    extraction_dispatch.resources.vertex_types = (md_gpu_buffer_resource_t){ .buffer = ctx->type_buf, .offset = 0, .usage = GPU_USAGE_READ };
    extraction_dispatch.resources.vertex_data = (md_gpu_buffer_resource_t){ .buffer = ctx->vert_buf, .offset = 0, .usage = GPU_USAGE_WRITE };
    extraction_dispatch.resources.edges = (md_gpu_buffer_resource_t){ .buffer = ctx->edge_buf, .offset = 0, .usage = GPU_USAGE_WRITE };
    extraction_dispatch.resources.ascending = (md_gpu_buffer_resource_t){ .buffer = ctx->ascending_buf, .offset = 0, .usage = GPU_USAGE_READ };
    extraction_dispatch.resources.descending = (md_gpu_buffer_resource_t){ .buffer = ctx->descending_buf, .offset = 0, .usage = GPU_USAGE_READ };
    extraction_dispatch.resources.voxel_to_vertex_idx = (md_gpu_buffer_resource_t){ .buffer = ctx->voxel_to_vert_buf, .offset = 0, .usage = GPU_USAGE_READ };
    extraction_dispatch.resources.meta = (md_gpu_buffer_resource_t){ .buffer = ctx->meta_buf, .offset = 0, .usage = GPU_USAGE_READ | GPU_USAGE_WRITE };
    extraction_dispatch.resources.volumeTex.image = volume;
    extraction_dispatch.resources.volumeTex.usage = GPU_USAGE_READ;
    extraction_dispatch.group_count[0] = DIV_UP(ctx->vert_cap, 64);
    extraction_dispatch.group_count[1] = 1;
    extraction_dispatch.group_count[2] = 1;
    topo_vertex_edge_extraction_cmd_dispatch(cmd, pip_vertex_edge_extraction, &extraction_dispatch);
    md_gpu_cmd_pop_debug_group(cmd);

#undef TOPO_FILL_COMMON_ARGS

    // Copy all results to CPU-visible staging in one batch, only if we have a staging buffers (discrete GPU).
    if (ctx->staging_buf) {
        md_gpu_cmd_barrier(cmd, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_TRANSFER);
        md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf, ctx->staging_buf, sizeof(topo_meta_t), 0, 0);
        md_gpu_cmd_copy_buffer(cmd, ctx->vert_buf, ctx->staging_verts, ctx->vert_cap * 4 * sizeof(float),       0, 0);
        md_gpu_cmd_copy_buffer(cmd, ctx->type_buf, ctx->staging_types, ctx->vert_cap * 1 * sizeof(uint32_t),    0, 0);
        md_gpu_cmd_copy_buffer(cmd, ctx->edge_buf, ctx->staging_edges, ctx->edge_cap * 2 * sizeof(uint32_t),    0, 0);
    }
}

bool md_topo_gpu_context_extract(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context) {
    if (!context || !out_graph) return false;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    ASSERT(out_graph->alloc);

    const topo_meta_t* meta = (const topo_meta_t*)md_gpu_buffer_cpu_ptr(ctx->staging_buf ? ctx->staging_buf : ctx->meta_buf);

    const uint32_t num_vertices = meta->vertex_count;
    const uint32_t num_edges    = meta->edge_count;

    if (meta->changed_read != 0) {
        MD_LOG_ERROR("md_topo_gpu_context_extract: path compression did not fully converge — results may be approximate");
    }

    MD_LOG_DEBUG("Topology: %u vertices, %u edges", num_vertices, num_edges);

    if (num_vertices == 0) return false;

    md_allocator_i* alloc = out_graph->alloc;
    MEMSET(out_graph, 0, sizeof(*out_graph));
    out_graph->alloc       = alloc;
    out_graph->num_vertices = num_vertices;
    out_graph->num_edges    = num_edges;

    out_graph->vertices = (md_topo_vert_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
    out_graph->types    = (md_topo_critical_point_type_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_critical_point_type_t));

    const float* vp = (const float*)md_gpu_buffer_cpu_ptr(ctx->staging_verts ? ctx->staging_verts : ctx->vert_buf);
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->vertices[i].x     = vp[i * 4 + 0];
        out_graph->vertices[i].y     = vp[i * 4 + 1];
        out_graph->vertices[i].z     = vp[i * 4 + 2];
        out_graph->vertices[i].value = vp[i * 4 + 3];
    }

    // Read vertex types directly from the staging buffer
    const uint32_t* tp = (const uint32_t*)md_gpu_buffer_cpu_ptr(ctx->staging_types ? ctx->staging_types : ctx->type_buf);
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->types[i] = (md_topo_critical_point_type_t)tp[i];
    }

    if (num_edges > 0) {
        out_graph->edges = (md_topo_edge_t*)md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        MEMCPY(out_graph->edges, md_gpu_buffer_cpu_ptr(ctx->staging_edges ? ctx->staging_edges : ctx->edge_buf), num_edges * sizeof(md_topo_edge_t));
    }

    return true;
}

#elif !MD_PLATFORM_OSX

#include <core/md_gl_util.h>
#include <topo_shaders.inl>
#include <GL/gl3w.h>

static GLuint create_compute_program(const char* source, size_t length) {
    GLuint program = 0;
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    if (md_gl_shader_compile(shader, (str_t){(const char*)source, length}, 0, 0)) {
        GLuint prog = glCreateProgram();
        if (md_gl_program_attach_and_link(prog, &shader, 1)) {
			program = prog;
        }
    }
    glDeleteShader(shader);
    return program;
}

// Shader program cache
static GLuint get_bidirectional_manifold_program(void) {
    static GLuint prog = 0;
    if (prog == 0) {
        prog = create_compute_program((const char*)bidirectional_manifold_comp, bidirectional_manifold_comp_size);
        if (prog == 0) {
            MD_LOG_ERROR("Failed to create bidirectional_manifold compute program");
        }
    }
    return prog;
}

static GLuint get_path_compression_program(void) {
    static GLuint prog = 0;
    if (prog == 0) {
        prog = create_compute_program((const char*)path_compression_comp, path_compression_comp_size);
        if (prog == 0) {
            MD_LOG_ERROR("Failed to create path_compression compute program");
        }
    }
    return prog;
}

static GLuint get_critical_points_program(void) {
    static GLuint prog = 0;
    if (prog == 0) {
        prog = create_compute_program((const char*)critical_points_comp, critical_points_comp_size);
        if (prog == 0) {
            MD_LOG_ERROR("Failed to create critical_points compute program");
        }
    }
    return prog;
}

static GLuint get_critical_point_compaction_program(void) {
	static GLuint prog = 0;
    if (prog == 0) {
        prog = create_compute_program((const char*)critical_point_compaction_comp, critical_point_compaction_comp_size);
        if (prog == 0) {
            MD_LOG_ERROR("Failed to create critical_point_compaction compute program");
        }
    }
    return prog;
}

static GLuint get_vertex_edge_extraction_program(void) {
    static GLuint prog = 0;
    if (prog == 0) {
        prog = create_compute_program((const char*)vertex_edge_extraction_comp, vertex_edge_extraction_comp_size);
        if (prog == 0) {
            MD_LOG_ERROR("Failed to create vertex_edge_extraction compute program");
        }
    }
    return prog;
}

// Helper to create a buffer
static GLuint create_buffer(size_t size, const void* data, GLenum usage) {
    GLuint buffer = 0;
    glGenBuffers(1, &buffer);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, size, data, usage);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    return buffer;
}

static void delete_buffer(GLuint buffer) {
    if (buffer) {
        glDeleteBuffers(1, &buffer);
    }
}

bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const md_grid_t* grid, float scalar_threshold) {
    if (!out_graph || vol_tex == 0 || !grid) {
        MD_LOG_ERROR("Invalid input: out_graph=%p, vol_tex=%u, grid=%p", (void*)out_graph, vol_tex, (void*)grid);
        return false;
    }

	md_gl_debug_push("Compute Extremum Graph");
    
    // Use heap allocator if none specified
    md_allocator_i* alloc = out_graph->alloc ? out_graph->alloc : md_get_heap_allocator();
    
    const uint32_t num_points = (uint32_t)(grid->dim[0] * grid->dim[1] * grid->dim[2]);
    const uint32_t workgroup_size = 8;
    const uint32_t num_workgroups[3] = {
        (grid->dim[0] + workgroup_size - 1) / workgroup_size,
        (grid->dim[1] + workgroup_size - 1) / workgroup_size,
        (grid->dim[2] + workgroup_size - 1) / workgroup_size
    };
    
    // Create UBO which is shared across all shaders
    struct {
        float index_to_world[4][4]; // mat4 in column-major order
        uint32_t dims[3];
        float scalar_threshold;
    } ubo_data;

    index_to_world_matrix(ubo_data.index_to_world, grid);
    ubo_data.dims[0] = grid->dim[0];
    ubo_data.dims[1] = grid->dim[1];
    ubo_data.dims[2] = grid->dim[2];
    // Use a permissive threshold by default to avoid dropping valid low-amplitude features
    ubo_data.scalar_threshold = scalar_threshold; // Set to >0.0 to filter noise if desired
    
    GLuint ubo_buf = create_buffer(sizeof(ubo_data), &ubo_data, GL_STATIC_DRAW);
    if (!ubo_buf) {
        return false;
    }
    
    bool success = false;
    GLuint ascending_buf = 0;
    GLuint descending_buf = 0;
    GLuint types_buf = 0;
    GLuint counts_buf = 0;
    GLuint indices_buf = 0;
    GLuint counter_buf = 0;
    
    // === Step 1: Compute bidirectional manifolds (steepest ascent/descent) ===
    GLuint manifold_prog = get_bidirectional_manifold_program();
    if (!manifold_prog) goto cleanup;

	md_gl_debug_push("Bidirectional Manifolds");
    
    ascending_buf  = create_buffer(num_points * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    descending_buf = create_buffer(num_points * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    
    glUseProgram(manifold_prog);
    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R32F);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, descending_buf);
    
    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	md_gl_debug_pop(); // Bidirectional Manifolds

#if DEBUG
    // Download ascending_buffer into local uint32_t array for debugging
    {
        uint32_t* asc_data = (uint32_t*)md_alloc(alloc, num_points * sizeof(uint32_t));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ascending_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_points * sizeof(uint32_t), asc_data);

        size_t count = 0;
        for (size_t i = 0; i < num_points; ++i) {
            if (asc_data[i] == i) {
                // This voxel is a minima (self-pointing in ascending manifold)
                count++;
            }
        }
        MD_LOG_INFO("Number of minima: %zu", count);
        md_free(alloc, asc_data, num_points * sizeof(uint32_t));
    }
#endif
    
    // === Step 2: Path compression (iteratively) ===
    GLuint compression_prog = get_path_compression_program();
    if (!compression_prog) goto cleanup;

    uint32_t num_iterations = 0;
    uint32_t max_dim = MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
    // Log2 ceiling
    while (max_dim > (1U << num_iterations)) {
        num_iterations++;
    }
    num_iterations += 2; // A couple of extra iterations to be safe

	md_gl_debug_push("Path Compression");

    GLuint changed_flag_buf = create_buffer(sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    
    glUseProgram(compression_prog);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, descending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, changed_flag_buf);
    
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, changed_flag_buf);
    for (size_t i = 0; i < num_iterations; i++) {
        glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	delete_buffer(changed_flag_buf);

	md_gl_debug_pop(); // Path Compression

#if DEBUG
    {  
        md_temp_t temp_scope = md_temp_begin();
        md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
        float*    vol_data  = (float*)md_alloc(temp_alloc, num_points * sizeof(float));
        uint32_t* asc_data  = (uint32_t*)md_alloc(temp_alloc, num_points * sizeof(uint32_t));
        uint32_t* desc_data = (uint32_t*)md_alloc(temp_alloc, num_points * sizeof(uint32_t));
        glBindTexture(GL_TEXTURE_3D, vol_tex);
        glGetTexImage(GL_TEXTURE_3D, 0, GL_RED, GL_FLOAT, vol_data);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ascending_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_points * sizeof(uint32_t), asc_data);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, descending_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_points * sizeof(uint32_t), desc_data);
        uint64_t vol_hash  = md_hash64(vol_data,  num_points * sizeof(float), 0);
        uint64_t asc_hash  = md_hash64(asc_data,  num_points * sizeof(uint32_t), 0);
        uint64_t desc_hash = md_hash64(desc_data, num_points * sizeof(uint32_t), 0);
        MD_LOG_INFO("Volume              hash: 0x%016llX", (unsigned long long)vol_hash);
        MD_LOG_INFO("Ascending  manifold hash: 0x%016llX", (unsigned long long)asc_hash);
        MD_LOG_INFO("Descending manifold hash: 0x%016llX", (unsigned long long)desc_hash);
        md_temp_end(temp_scope);
    }
#endif
    
    // === Step 3: Identify critical points ===
    GLuint critical_prog = get_critical_points_program();
    if (!critical_prog) goto cleanup;

	md_gl_debug_push("Critical Points");
    
    types_buf = create_buffer(num_points * sizeof(int), NULL, GL_DYNAMIC_COPY);
    
    // Counts buffer: single total critical-point count
    uint32_t counts_init = 0;
    counts_buf = create_buffer(sizeof(counts_init), &counts_init, GL_DYNAMIC_COPY);
    
    glUseProgram(critical_prog);
    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R32F);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, descending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, counts_buf);
    
    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	md_gl_debug_pop(); // Critical Points
    
    // Read back single count
    uint32_t num_vertices = 0;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, counts_buf);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(uint32_t), &num_vertices);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    
    uint32_t num_edges = 8 * num_vertices; // conservative estimate; updated after extraction
    
    MD_LOG_DEBUG("Topology: %u critical points", num_vertices);
    
    // === Step 4: Compact critical point indices into ordered array ===
    GLuint compaction_prog = get_critical_point_compaction_program();
    if (!compaction_prog) goto cleanup;

	md_gl_debug_push("Critical Point Compaction");

    // Allocate indices buffer for all critical points
    if (num_vertices > 0) {
        indices_buf = create_buffer(num_vertices * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        if (!indices_buf) {
            MD_LOG_ERROR("Failed to allocate indices buffer");
            goto cleanup;
        }
    } else {
        indices_buf = 0;
    }
    
    // Single atomic write cursor, starting at 0
    uint32_t counter_init = 0;
    counter_buf = create_buffer(sizeof(counter_init), &counter_init, GL_DYNAMIC_COPY);

    // Per-vertex type buffer (written by compaction, read by extraction)
    GLuint type_buf = create_buffer(num_vertices * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    
    glUseProgram(compaction_prog);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, indices_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, counter_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, type_buf);

    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	md_gl_debug_pop(); // Critical Point Compaction

#if 0
    // Print out all of the critical point indices found (sorted)
    #if DEBUG
    {
        md_temp_t temp_scope = md_temp_begin();
        md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
        if (num_vertices > 0) {
            uint32_t* data = (uint32_t*)md_alloc(temp_alloc, num_vertices * sizeof(uint32_t));
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, indices_buf);
            glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_vertices * sizeof(uint32_t), data);
            qsort(data, num_vertices, sizeof(uint32_t), uint_compare);
            printf("Maxima indices:");
            for (uint32_t i = 0; i < num_vertices; i++) {
                printf("  %u", data[i]);
            }
            printf("\n");
        }
        md_temp_end(temp_scope);
    }
    #endif
#endif
    
    // === Step 5: Extract vertices and edges using GPU shader ===
    GLuint extraction_prog = get_vertex_edge_extraction_program();
    if (!extraction_prog) goto cleanup;

	md_gl_debug_push("Graph Extraction");
    
    // Create output buffers for vertices and edges
    GLuint vert_buf = create_buffer(num_vertices * sizeof(md_topo_vert_t), NULL, GL_DYNAMIC_COPY);
    GLuint edge_buf = create_buffer(num_edges    * sizeof(md_topo_edge_t), NULL, GL_DYNAMIC_COPY);
    
    uint32_t edge_count_init = 0;
    GLuint edge_count_buf = create_buffer(sizeof(uint32_t), &edge_count_init, GL_DYNAMIC_COPY);
    
    glUseProgram(extraction_prog);
    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R32F);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, indices_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, type_buf);  // per-vertex type
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, descending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, vert_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, edge_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, edge_count_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, types_buf);         // voxel -> vertex index
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, counts_buf);        // num_vertices
    
    uint32_t num_extraction_workgroups = (num_vertices + 63) / 64;
    glDispatchCompute(num_extraction_workgroups, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	md_gl_debug_pop(); // Graph Extraction

    md_gl_debug_push("Readback Results");

    // Read back edge count
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, edge_count_buf);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(uint32_t), &num_edges);
    
    // Read back vertex data
    md_topo_vert_t* vertices = NULL;
    if (num_vertices > 0) {
        vertices = md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, vert_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_vertices * sizeof(md_topo_vert_t), vertices);
    }

    // Read back vertex types
    md_topo_critical_point_type_t* types_out = NULL;
    if (num_vertices > 0) {
        types_out = md_alloc(alloc, num_vertices * sizeof(md_topo_critical_point_type_t));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, type_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_vertices * sizeof(uint32_t), types_out);
    }

    // Read back edges
    md_topo_edge_t* edges = NULL;
    if (num_edges > 0) {
        edges = md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, edge_buf);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, num_edges * sizeof(md_topo_edge_t), edges);
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    md_gl_debug_pop(); // Readback Results
    
    delete_buffer(indices_buf);
    delete_buffer(type_buf);
    delete_buffer(vert_buf);
    delete_buffer(edge_buf);
    delete_buffer(edge_count_buf);
    
    // Fill output structure
    MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
    out_graph->num_vertices = num_vertices;
    out_graph->vertices     = vertices;
    out_graph->types        = types_out;
    out_graph->num_edges    = num_edges;
    out_graph->edges        = edges;
    out_graph->alloc        = alloc;
    
    success = true;
cleanup:
    delete_buffer(ubo_buf);
    delete_buffer(ascending_buf);
    delete_buffer(descending_buf);
    delete_buffer(types_buf);
    delete_buffer(counts_buf);
    delete_buffer(counter_buf);
    
    md_gl_debug_pop();
    return success;
}

#else

// macOS stub (no GL, no md_gpu)
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const md_grid_t* grid, float scalar_threshold) {
    (void)out_graph;
    (void)vol_tex;
    (void)grid;
    (void)scalar_threshold;
    MD_LOG_ERROR("Topology GPU computation not available (enable MD_ENABLE_GPU)");
    return false;
}

#endif

void md_topo_simplify(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* in_graph,
    float threshold, bool do_prune_duplicate_saddles)
{
    ASSERT(out_graph);
    ASSERT(in_graph);

    md_allocator_i* alloc = out_graph->alloc ? out_graph->alloc : md_get_heap_allocator();
    md_topo_extremum_graph_free(out_graph);
    out_graph->alloc = alloc;

    if (!in_graph->vertices || in_graph->num_vertices == 0) return;

    md_allocator_i* conflicts[] = { alloc };
    md_temp_t temp_scope = md_temp_begin_avoid(conflicts, ARRAY_SIZE(conflicts));
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    // Working copies (may grow during pruning)
    md_array(int)             vertex_type = md_array_create(int,            in_graph->num_vertices, temp_alloc);
    md_array(md_topo_vert_t)  vertex      = md_array_create(md_topo_vert_t, in_graph->num_vertices, temp_alloc);
    md_array(md_topo_edge_t)  edge        = md_array_create(md_topo_edge_t, in_graph->num_edges,    temp_alloc);
    md_array(md_array(int))   vertex_adj  = md_array_create(md_array(int),  in_graph->num_vertices, temp_alloc);

    MEMCPY(vertex_type, in_graph->types,    in_graph->num_vertices * sizeof(int));
    MEMCPY(vertex,      in_graph->vertices, in_graph->num_vertices * sizeof(md_topo_vert_t));
    MEMSET(vertex_adj,  0,                  in_graph->num_vertices * sizeof(md_array(int)));
    MEMCPY(edge,        in_graph->edges,    in_graph->num_edges    * sizeof(md_topo_edge_t));

    // Build adjacency
    for (size_t i = 0; i < md_array_size(edge); ++i) {
        md_topo_edge_t e = edge[i];
        md_array_push(vertex_adj[e.from], (int)e.to,   temp_alloc);
        md_array_push(vertex_adj[e.to],   (int)e.from, temp_alloc);
    }

    // Kill vertices below threshold
    if (threshold > 0.0f) {
        for (size_t i = 0; i < md_array_size(vertex); ++i) {
            if (vertex[i].value < threshold) {
                vertex_type[i] = 0;
            }
        }
    }

    // Prune duplicate saddles between maxima pairs
    if (do_prune_duplicate_saddles) {
        // Split multi-connected saddles (3+ adjacencies) into per-pair saddles
        for (size_t i = 0; i < md_array_size(vertex); ++i) {
            if (vertex_type[i] != MD_TOPO_SPLIT_SADDLE) continue;
            size_t num_adj = md_array_size(vertex_adj[i]);
            if (num_adj <= 2) continue;
            for (size_t j = 0; j < num_adj - 1; ++j) {
                for (size_t k = j + 1; k < num_adj; ++k) {
                    int max_a = vertex_adj[i][j];
                    int max_b = vertex_adj[i][k];
                    md_topo_vert_t new_saddle = vertex[i];
                    md_array_push(vertex,      new_saddle,           temp_alloc);
                    md_array_push(vertex_adj,  NULL,                 temp_alloc);
                    md_array_push(vertex_type, MD_TOPO_SPLIT_SADDLE, temp_alloc);
                    int ns = (int)(md_array_size(vertex) - 1);
                    md_array_push(vertex_adj[ns], max_a, temp_alloc);
                    md_array_push(vertex_adj[ns], max_b, temp_alloc);
                    md_topo_edge_t ea = { (uint32_t)max_a, (uint32_t)ns };
                    md_topo_edge_t eb = { (uint32_t)max_b, (uint32_t)ns };
                    md_array_push(edge, ea, temp_alloc);
                    md_array_push(edge, eb, temp_alloc);
                }
            }
            vertex_type[i] = 0;  // Kill original multi-saddle
        }

        // For each maxima pair, keep only the highest-valued connecting saddle
        md_array(int) saddle_list = 0;
        for (size_t i = 0; i < md_array_size(vertex) - 1; ++i) {
            if (vertex_type[i] != MD_TOPO_MAXIMUM) continue;
            for (size_t j = i + 1; j < md_array_size(vertex); ++j) {
                if (vertex_type[j] != MD_TOPO_MAXIMUM) continue;

                md_array_shrink(saddle_list, 0);
                for (size_t k = 0; k < md_array_size(vertex); ++k) {
                    if (vertex_type[k] != MD_TOPO_SPLIT_SADDLE) continue;
                    if (md_array_size(vertex_adj[k]) != 2) continue;
                    bool ci = false, cj = false;
                    for (size_t m = 0; m < 2; ++m) {
                        int v = vertex_adj[k][m];
                        if (v == (int)i) ci = true;
                        if (v == (int)j) cj = true;
                    }
                    if (ci && cj) md_array_push(saddle_list, (int)k, temp_alloc);
                }

                if (md_array_size(saddle_list) > 1) {
                    float best_val = -FLT_MAX;
                    int   best_idx = -1;
                    for (size_t k = 0; k < md_array_size(saddle_list); ++k) {
                        int idx = saddle_list[k];
                        if (vertex[idx].value > best_val) {
                            best_val = vertex[idx].value;
                            best_idx = idx;
                        }
                    }
                    for (size_t k = 0; k < md_array_size(saddle_list); ++k) {
                        int idx = saddle_list[k];
                        if (idx != best_idx) vertex_type[idx] = 0;
                    }
                }
            }
        }
    }

    // Build compact remap: surviving vertices get sequential indices
    uint32_t num_vertices = 0;
    md_array(int) vertex_remap = md_array_create(int, md_array_size(vertex), temp_alloc);
    MEMSET(vertex_remap, -1, md_array_bytes(vertex_remap));
    for (size_t i = 0; i < md_array_size(vertex_type); ++i) {
        if (vertex_type[i] == 0) continue;
        vertex_remap[i] = (int)num_vertices++;
    }

    if (num_vertices == 0) {
        goto done;
    }

    out_graph->num_vertices = num_vertices;
    out_graph->vertices     = (md_topo_vert_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
    out_graph->types        = (md_topo_critical_point_type_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_critical_point_type_t));

    for (size_t i = 0; i < md_array_size(vertex_type); ++i) {
        int idx = vertex_remap[i];
        if (idx == -1) continue;
        out_graph->vertices[idx] = vertex[i];
        out_graph->types[idx]    = (md_topo_critical_point_type_t)vertex_type[i];
    }

    // Count surviving edges
    uint32_t num_edges = 0;
    for (size_t i = 0; i < md_array_size(edge); ++i) {
        if (vertex_remap[edge[i].from] != -1 && vertex_remap[edge[i].to] != -1) num_edges++;
    }
    out_graph->num_edges = num_edges;
    out_graph->edges     = (md_topo_edge_t*)md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));

    uint32_t ec = 0;
    for (size_t i = 0; i < md_array_size(edge); ++i) {
        int from = vertex_remap[edge[i].from];
        int to   = vertex_remap[edge[i].to];
        if (from != -1 && to != -1) {
            out_graph->edges[ec].from = (uint32_t)from;
            out_graph->edges[ec].to   = (uint32_t)to;
            ec++;
        }
    }
done:
    md_temp_end(temp_scope);
}

void md_topo_count_vertex_types(uint32_t out_counts[MD_TOPO_NUM_TYPES], const md_topo_extremum_graph_t* graph) {
    if (!graph || !out_counts) return;
    MEMSET(out_counts, 0, sizeof(uint32_t) * MD_TOPO_NUM_TYPES);
    for (uint32_t i = 0; i < graph->num_vertices; ++i) {
        md_topo_critical_point_type_t type = graph->types[i];
        if (0 <= type && type < MD_TOPO_NUM_TYPES) {
            out_counts[type]++;
        }
    }
}

void md_topo_extremum_graph_free(md_topo_extremum_graph_t* graph) {
    if (graph && graph->alloc) {
        md_allocator_i* alloc = graph->alloc;
        if (graph->vertices) md_free(alloc, graph->vertices, graph->num_vertices * sizeof(md_topo_vert_t));
        if (graph->types)    md_free(alloc, graph->types,    graph->num_vertices * sizeof(md_topo_critical_point_type_t));
        if (graph->edges)    md_free(alloc, graph->edges,    graph->num_edges    * sizeof(md_topo_edge_t));
        MEMSET(graph, 0, sizeof(md_topo_extremum_graph_t));
        graph->alloc = alloc;
    }
}

void md_topo_extremum_graph_copy(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* src_graph) {
    md_allocator_i* alloc = out_graph->alloc ? out_graph->alloc : md_get_heap_allocator();
    MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
    out_graph->alloc        = alloc;
    out_graph->num_vertices = src_graph->num_vertices;
    out_graph->num_edges    = src_graph->num_edges;

    out_graph->vertices = (md_topo_vert_t*)md_alloc(alloc, src_graph->num_vertices * sizeof(md_topo_vert_t));
    out_graph->types    = (md_topo_critical_point_type_t*)md_alloc(alloc, src_graph->num_vertices * sizeof(md_topo_critical_point_type_t));
    out_graph->edges    = (md_topo_edge_t*)md_alloc(alloc, src_graph->num_edges * sizeof(md_topo_edge_t));

    MEMCPY(out_graph->vertices, src_graph->vertices, src_graph->num_vertices * sizeof(md_topo_vert_t));
    MEMCPY(out_graph->types,    src_graph->types,    src_graph->num_vertices * sizeof(md_topo_critical_point_type_t));
    MEMCPY(out_graph->edges,    src_graph->edges,    src_graph->num_edges    * sizeof(md_topo_edge_t));
}

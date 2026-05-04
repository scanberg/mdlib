#include <md_topo.h>

#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_grid.h>
#include <core/md_array.h>

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

// Pipeline cache
static md_gpu_device_t cached_device = NULL;
static md_gpu_compute_pipeline_t pip_bidirectional_manifold    = NULL;
static md_gpu_compute_pipeline_t pip_path_compression          = NULL;
static md_gpu_compute_pipeline_t pip_critical_points           = NULL;
static md_gpu_compute_pipeline_t pip_topo_setup                = NULL;
static md_gpu_compute_pipeline_t pip_critical_point_compaction = NULL;
static md_gpu_compute_pipeline_t pip_vertex_edge_extraction    = NULL;

static md_gpu_compute_pipeline_t ensure_pipeline(md_gpu_device_t device, md_gpu_compute_pipeline_t* cached, const void* blob_start, size_t blob_size, const char* name, uint32_t wg_x, uint32_t wg_y, uint32_t wg_z) {
    if (cached_device != device) {
        // Device changed, invalidate all cached pipelines
        if (pip_bidirectional_manifold)    { md_gpu_compute_pipeline_destroy(pip_bidirectional_manifold);    pip_bidirectional_manifold    = NULL; }
        if (pip_path_compression)          { md_gpu_compute_pipeline_destroy(pip_path_compression);          pip_path_compression          = NULL; }
        if (pip_critical_points)           { md_gpu_compute_pipeline_destroy(pip_critical_points);           pip_critical_points           = NULL; }
        if (pip_topo_setup)                { md_gpu_compute_pipeline_destroy(pip_topo_setup);                pip_topo_setup                = NULL; }
        if (pip_critical_point_compaction) { md_gpu_compute_pipeline_destroy(pip_critical_point_compaction); pip_critical_point_compaction = NULL; }
        if (pip_vertex_edge_extraction)    { md_gpu_compute_pipeline_destroy(pip_vertex_edge_extraction);    pip_vertex_edge_extraction    = NULL; }
        cached_device = device;
    }
    if (*cached == NULL) {
        md_gpu_compute_pipeline_desc_t desc = {
            .shader_bytes     = blob_start,
            .shader_byte_size = blob_size,
            .threadgroup_size = { wg_x, wg_y, wg_z },
        };
        *cached = md_gpu_compute_pipeline_create(device, &desc);
        if (*cached == NULL) {
            MD_LOG_ERROR("Failed to create compute pipeline: %s", name);
        }
    }
    return *cached;
}

// ---------------------------------------------------------------------------
// Scratch buffer layout:
//   ascending          @ [0 * stride .. 1 * stride)
//   descending         @ [1 * stride .. 2 * stride)
//   types              @ [2 * stride .. 3 * stride)
//   voxel_to_vert_idx  @ [3 * stride .. 4 * stride)
// where stride = ALIGN_UP(num_points * sizeof(uint32_t), 256)
//
// Meta buffer layout (device-local, 256-byte aligned slots):
//   [   0..  16) counts[4]         critical-point type counts (critical_points output)
//   [ 256.. 260) changed_read      path-compression flag, read  by shader
//   [ 512.. 516) changed_write     path-compression flag, write by shader
//   [ 768.. 784) counters[4]       write-cursor offsets for compaction (topo_setup output)
//   [1024..1044) type_counts[5]    per-type counts + total for extraction (topo_setup output)
//   [1280..1284) edge_count        actual edge count (extraction output)
//
// Staging buffer layout (CPU-visible, 256-byte aligned slots):
//   [   0..  16) counts[4]         readback of type counts
//   [ 256.. 260) changed           readback of convergence flag
//   [ 512.. 516) edge_count        readback of actual edge count
// ---------------------------------------------------------------------------
#define TOPO_BUF_ALIGN            256

#define TOPO_META_COUNTS_OFF         0
#define TOPO_META_CHANGED_R_OFF    256
#define TOPO_META_CHANGED_W_OFF    512
#define TOPO_META_COUNTERS_OFF     768
#define TOPO_META_TYPE_COUNTS_OFF 1024
#define TOPO_META_EDGE_COUNT_OFF  1280
#define TOPO_META_BUF_SIZE        1536

#define TOPO_STAGING_COUNTS_OFF    0
#define TOPO_STAGING_CHANGED_OFF   256
#define TOPO_STAGING_EDGE_CNT_OFF  512
#define TOPO_STAGING_BUF_SIZE      768

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

    // Persistent GPU buffers — scratch
    md_gpu_buffer_t scratch_buf;   // scratch_stride*4 bytes, device-local (stride=ALIGN_UP(num_points*4,256))
    md_gpu_buffer_t meta_buf;      // TOPO_META_BUF_SIZE bytes, device-local
    md_gpu_buffer_t staging_buf;   // TOPO_STAGING_BUF_SIZE bytes, CPU-visible

    // Persistent GPU buffers — results (pre-allocated to worst-case capacity)
    md_gpu_buffer_t indices_buf;   // u32[vert_cap]
    md_gpu_buffer_t vert_buf;      // float4[vert_cap]
    md_gpu_buffer_t staging_verts; // float4[vert_cap], CPU-visible
    md_gpu_buffer_t edge_buf;      // md_topo_edge_t[edge_cap]
    md_gpu_buffer_t staging_edges; // md_topo_edge_t[edge_cap], CPU-visible
};

md_topo_gpu_context_t* md_topo_gpu_context_create(md_gpu_device_t device, uint32_t dim_x, uint32_t dim_y, uint32_t dim_z) {
    if (!device || !dim_x || !dim_y || !dim_z) {
        MD_LOG_ERROR("md_topo_gpu_context_create: invalid arguments");
        return NULL;
    }

    if (!ensure_pipeline(device, &pip_bidirectional_manifold,    topo_bidirectional_manifold_start,    topo_bidirectional_manifold_size(),    "bidirectional_manifold",    8,  8,  8) ||
        !ensure_pipeline(device, &pip_path_compression,          topo_path_compression_start,          topo_path_compression_size(),          "path_compression",          8,  8,  8) ||
        !ensure_pipeline(device, &pip_critical_points,           topo_critical_points_start,           topo_critical_points_size(),           "critical_points",           8,  8,  8) ||
        !ensure_pipeline(device, &pip_topo_setup,                topo_topo_setup_start,                topo_topo_setup_size(),                "topo_setup",                1,  1,  1) ||
        !ensure_pipeline(device, &pip_critical_point_compaction, topo_critical_point_compaction_start, topo_critical_point_compaction_size(), "critical_point_compaction", 8,  8,  8) ||
        !ensure_pipeline(device, &pip_vertex_edge_extraction,    topo_vertex_edge_extraction_start,    topo_vertex_edge_extraction_size(),    "vertex_edge_extraction",    64, 1,  1))
    {
        MD_LOG_ERROR("md_topo_gpu_context_create: failed to create compute pipelines");
        return NULL;
    }

    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;

    ctx->device     = device;
    ctx->num_points = dim_x * dim_y * dim_z;
    ctx->dim[0]     = dim_x;
    ctx->dim[1]     = dim_y;
    ctx->dim[2]     = dim_z;

    uint32_t vert_cap = ctx->num_points / TOPO_VERT_RATIO;
    if (vert_cap < TOPO_VERT_CAP_MIN) vert_cap = TOPO_VERT_CAP_MIN;
    uint32_t edge_cap = vert_cap * TOPO_EDGE_RATIO;
    ctx->vert_cap = vert_cap;
    ctx->edge_cap = edge_cap;

    size_t scratch_stride = ((size_t)ctx->num_points * sizeof(uint32_t) + (TOPO_BUF_ALIGN - 1)) & ~(size_t)(TOPO_BUF_ALIGN - 1);
    ctx->scratch_buf   = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = scratch_stride * 4 });
    ctx->meta_buf      = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = TOPO_META_BUF_SIZE });
    ctx->staging_buf   = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = TOPO_STAGING_BUF_SIZE, .flags = MD_GPU_BUFFER_CPU_VISIBLE });
    ctx->indices_buf   = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * sizeof(uint32_t) });
    ctx->vert_buf      = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * 4 * sizeof(float) });
    ctx->staging_verts = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = vert_cap * 4 * sizeof(float), .flags = MD_GPU_BUFFER_CPU_VISIBLE });
    ctx->edge_buf      = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = edge_cap * sizeof(md_topo_edge_t) });
    ctx->staging_edges = md_gpu_buffer_create(device, &(md_gpu_buffer_desc_t){ .size = edge_cap * sizeof(md_topo_edge_t), .flags = MD_GPU_BUFFER_CPU_VISIBLE });

    if (!ctx->scratch_buf || !ctx->meta_buf || !ctx->staging_buf ||
        !ctx->indices_buf || !ctx->vert_buf || !ctx->staging_verts ||
        !ctx->edge_buf    || !ctx->staging_edges)
    {
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
    if (ctx->scratch_buf)   md_gpu_buffer_destroy(ctx->scratch_buf);
    if (ctx->meta_buf)      md_gpu_buffer_destroy(ctx->meta_buf);
    if (ctx->staging_buf)   md_gpu_buffer_destroy(ctx->staging_buf);
    if (ctx->indices_buf)   md_gpu_buffer_destroy(ctx->indices_buf);
    if (ctx->vert_buf)      md_gpu_buffer_destroy(ctx->vert_buf);
    if (ctx->staging_verts) md_gpu_buffer_destroy(ctx->staging_verts);
    if (ctx->edge_buf)      md_gpu_buffer_destroy(ctx->edge_buf);
    if (ctx->staging_edges) md_gpu_buffer_destroy(ctx->staging_edges);
    free(ctx);
}

void md_topo_gpu_cmd_record(md_gpu_command_buffer_t cmd, md_topo_gpu_context_t* context, md_gpu_image_t volume, const md_grid_t* grid, float scalar_threshold) {
    if (!cmd || !context || !volume || !grid) return;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    struct {
        float    index_to_world[4][4];
        uint32_t dims[4];
        float    scalar_threshold;
        float    _pad[3];
    } ubo;
    index_to_world_matrix(ubo.index_to_world, grid);
    ubo.dims[0] = grid->dim[0];
    ubo.dims[1] = grid->dim[1];
    ubo.dims[2] = grid->dim[2];
    ubo.dims[3] = 0;
    ubo.scalar_threshold = scalar_threshold;
    ubo._pad[0] = ubo._pad[1] = ubo._pad[2] = 0.0f;

    const uint32_t wg[3]  = { DIV_UP(grid->dim[0], 8), DIV_UP(grid->dim[1], 8), DIV_UP(grid->dim[2], 8) };
    const size_t   stride = ((size_t)ctx->num_points * sizeof(uint32_t) + (TOPO_BUF_ALIGN - 1)) & ~(size_t)(TOPO_BUF_ALIGN - 1);

    // Per-call resets (voxel_to_vert_idx = -1, edge_count = 0, counts = 0)
    md_gpu_cmd_fill_buffer(cmd, ctx->scratch_buf, stride * 3,               stride, 0xFF); // voxel_to_vert_idx = -1
    md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf,    TOPO_META_COUNTS_OFF,     16,     0);    // counts[4] = 0
    md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf,    TOPO_META_EDGE_COUNT_OFF, 4,      0);    // edge_count = 0
    md_gpu_cmd_barrier(cmd);

    // Step 1: Bidirectional manifold
    md_gpu_cmd_bind_compute_pipeline(cmd, pip_bidirectional_manifold);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_image(cmd, 0, volume);
    md_gpu_cmd_bind_buffer_range(cmd, 0, ctx->scratch_buf, 0,      stride);  // ascending
    md_gpu_cmd_bind_buffer_range(cmd, 1, ctx->scratch_buf, stride, stride);  // descending
    md_gpu_cmd_dispatch(cmd, wg[0], wg[1], wg[2]);
    md_gpu_cmd_barrier(cmd);

    // Step 2: Path compression (iterative, GPU-side early-exit via changed flag)
    // Each dispatch does one grandparent pointer-jump. A path of length L needs
    // ceil(log2(L)) iterations; 2 * ceil(log2(max_dim)) is a safe upper bound.
    uint32_t num_iterations = 0;
    {
        uint32_t max_dim = (uint32_t)MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
        while (max_dim > (1U << num_iterations)) num_iterations++;
        num_iterations *= 2;
    }

    md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf, TOPO_META_CHANGED_R_OFF, 4, 1);
    md_gpu_cmd_barrier_buffer_ex(cmd, ctx->meta_buf, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);

    md_gpu_cmd_bind_compute_pipeline(cmd, pip_path_compression);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_buffer_range(cmd, 0, ctx->scratch_buf, 0,      stride);          // ascending
    md_gpu_cmd_bind_buffer_range(cmd, 1, ctx->scratch_buf, stride, stride);          // descending
    md_gpu_cmd_bind_buffer_range(cmd, 2, ctx->meta_buf, TOPO_META_CHANGED_R_OFF, 4); // changed_read
    md_gpu_cmd_bind_buffer_range(cmd, 3, ctx->meta_buf, TOPO_META_CHANGED_W_OFF, 4); // changed_write
    for (uint32_t i = 0; i < num_iterations; i++) {
        md_gpu_cmd_fill_buffer(cmd, ctx->meta_buf, TOPO_META_CHANGED_W_OFF, 4, 0);
        md_gpu_cmd_barrier_buffer_ex(cmd, ctx->meta_buf, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);

        md_gpu_cmd_dispatch(cmd, wg[0], wg[1], wg[2]);

        md_gpu_cmd_barrier_buffer_ex(cmd, ctx->scratch_buf, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);
        md_gpu_cmd_barrier_buffer_ex(cmd, ctx->meta_buf,    MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_TRANSFER);

        // Rotate changed_write → changed_read (non-overlapping self-copy)
        md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf, ctx->meta_buf, 4, TOPO_META_CHANGED_W_OFF, TOPO_META_CHANGED_R_OFF);
        md_gpu_cmd_barrier_buffer_ex(cmd, ctx->meta_buf, MD_GPU_BARRIER_STAGE_TRANSFER, MD_GPU_BARRIER_STAGE_COMPUTE);
    }

    // Step 3: Critical-point detection → meta_buf.counts[4]
    md_gpu_cmd_bind_compute_pipeline(cmd, pip_critical_points);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_image(cmd, 0, volume);
    md_gpu_cmd_bind_buffer_range(cmd, 0, ctx->scratch_buf, 0,          stride);  // ascending
    md_gpu_cmd_bind_buffer_range(cmd, 1, ctx->scratch_buf, stride,     stride);  // descending
    md_gpu_cmd_bind_buffer_range(cmd, 2, ctx->scratch_buf, stride * 2, stride);  // types
    md_gpu_cmd_bind_buffer_range(cmd, 3, ctx->meta_buf, TOPO_META_COUNTS_OFF, 16); // counts
    md_gpu_cmd_dispatch(cmd, wg[0], wg[1], wg[2]);
    md_gpu_cmd_barrier_buffer_ex(cmd, ctx->meta_buf, MD_GPU_BARRIER_STAGE_COMPUTE, MD_GPU_BARRIER_STAGE_COMPUTE);

    // Step 3.5: topo_setup (1,1,1) — prefix-sum: counts → counters + type_counts
    md_gpu_cmd_bind_compute_pipeline(cmd, pip_topo_setup);
    md_gpu_cmd_bind_buffer_range(cmd, 0, ctx->meta_buf, TOPO_META_COUNTS_OFF,       16); // counts[4]
    md_gpu_cmd_bind_buffer_range(cmd, 1, ctx->meta_buf, TOPO_META_COUNTERS_OFF,     16); // counters[4]
    md_gpu_cmd_bind_buffer_range(cmd, 2, ctx->meta_buf, TOPO_META_TYPE_COUNTS_OFF,  20); // type_counts[5]
    md_gpu_cmd_dispatch(cmd, 1, 1, 1);
    md_gpu_cmd_barrier(cmd);

    // Step 4: Critical-point compaction
    md_gpu_cmd_bind_compute_pipeline(cmd, pip_critical_point_compaction);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_buffer_range(cmd, 0, ctx->scratch_buf, stride * 2, stride);           // types
    md_gpu_cmd_bind_buffer(cmd, 1, ctx->indices_buf);                                     // indices (out)
    md_gpu_cmd_bind_buffer_range(cmd, 2, ctx->meta_buf, TOPO_META_COUNTERS_OFF,    16);   // counters
    md_gpu_cmd_bind_buffer_range(cmd, 3, ctx->scratch_buf, stride * 3, stride);           // voxel_to_vert_idx
    md_gpu_cmd_dispatch(cmd, wg[0], wg[1], wg[2]);
    md_gpu_cmd_barrier(cmd);

    // Step 5: Vertex + edge extraction
    // Dispatched at worst-case capacity; shader exits early for indices >= num_vertices.
    md_gpu_cmd_bind_compute_pipeline(cmd, pip_vertex_edge_extraction);
    md_gpu_cmd_push_constants(cmd, &ubo, sizeof(ubo));
    md_gpu_cmd_bind_buffer(cmd, 0, ctx->indices_buf);                                     // cp_indices
    md_gpu_cmd_bind_buffer_range(cmd, 1, ctx->meta_buf, TOPO_META_TYPE_COUNTS_OFF, 20);   // type_counts
    md_gpu_cmd_bind_buffer_range(cmd, 2, ctx->scratch_buf, 0,          stride);           // ascending
    md_gpu_cmd_bind_buffer_range(cmd, 3, ctx->scratch_buf, stride,     stride);           // descending
    md_gpu_cmd_bind_buffer(cmd, 4, ctx->vert_buf);
    md_gpu_cmd_bind_buffer(cmd, 5, ctx->edge_buf);
    md_gpu_cmd_bind_buffer_range(cmd, 6, ctx->meta_buf, TOPO_META_EDGE_COUNT_OFF, 4);     // edge_count
    md_gpu_cmd_bind_buffer_range(cmd, 7, ctx->scratch_buf, stride * 3, stride);           // voxel_to_vert_idx
    md_gpu_cmd_bind_image(cmd, 0, volume);
    md_gpu_cmd_dispatch(cmd, DIV_UP(ctx->vert_cap, 64), 1, 1);
    md_gpu_cmd_barrier(cmd);

    // Copy all results + counts to CPU-visible staging in one batch
    md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf,  ctx->staging_buf,    16, TOPO_META_COUNTS_OFF,     TOPO_STAGING_COUNTS_OFF);
    md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf,  ctx->staging_buf,     4, TOPO_META_CHANGED_R_OFF,  TOPO_STAGING_CHANGED_OFF);
    md_gpu_cmd_copy_buffer(cmd, ctx->meta_buf,  ctx->staging_buf,     4, TOPO_META_EDGE_COUNT_OFF, TOPO_STAGING_EDGE_CNT_OFF);
    md_gpu_cmd_copy_buffer(cmd, ctx->vert_buf,  ctx->staging_verts,   ctx->vert_cap * 4 * sizeof(float),        0, 0);
    md_gpu_cmd_copy_buffer(cmd, ctx->edge_buf,  ctx->staging_edges,   ctx->edge_cap * sizeof(md_topo_edge_t),   0, 0);
}

bool md_topo_gpu_context_get_result(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context) {
    if (!context || !out_graph) return false;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    ASSERT(out_graph->alloc);

    const uint8_t* stg = (const uint8_t*)md_gpu_buffer_cpu_ptr(ctx->staging_buf);

    const uint32_t* counts = (const uint32_t*)(stg + TOPO_STAGING_COUNTS_OFF);
    const uint32_t num_maxima        = counts[0];
    const uint32_t num_split_saddles = counts[1];
    const uint32_t num_minima        = counts[2];
    const uint32_t num_join_saddles  = counts[3];
    const uint32_t num_vertices      = num_maxima + num_split_saddles + num_minima + num_join_saddles;
    const uint32_t num_edges         = *(const uint32_t*)(stg + TOPO_STAGING_EDGE_CNT_OFF);

    if (*(const uint32_t*)(stg + TOPO_STAGING_CHANGED_OFF) != 0) {
        MD_LOG_ERROR("md_topo_gpu_context_get_result: path compression did not fully converge — results may be approximate");
    }

    MD_LOG_DEBUG("Topology: %u maxima, %u split saddles, %u minima, %u join saddles (total %u vertices, %u edges)",
                 num_maxima, num_split_saddles, num_minima, num_join_saddles, num_vertices, num_edges);

    if (num_vertices == 0) return false;

    md_allocator_i* alloc = out_graph->alloc;
    MEMSET(out_graph, 0, sizeof(*out_graph));
    out_graph->alloc       = alloc;
    out_graph->num_vertices = num_vertices;
    out_graph->num_edges    = num_edges;

    out_graph->vertices = (md_topo_vert_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
    out_graph->types    = (md_topo_critical_point_type_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_critical_point_type_t));

    const float* vp = (const float*)md_gpu_buffer_cpu_ptr(ctx->staging_verts);
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->vertices[i].x     = vp[i * 4 + 0];
        out_graph->vertices[i].y     = vp[i * 4 + 1];
        out_graph->vertices[i].z     = vp[i * 4 + 2];
        out_graph->vertices[i].value = vp[i * 4 + 3];
    }

    // Vertices are type-sorted: [maxima, split_saddles, minima, join_saddles]
    uint32_t off = 0;
    for (uint32_t i = 0; i < num_maxima;        i++) out_graph->types[off++] = MD_TOPO_MAXIMUM;
    for (uint32_t i = 0; i < num_split_saddles; i++) out_graph->types[off++] = MD_TOPO_SPLIT_SADDLE;
    for (uint32_t i = 0; i < num_minima;        i++) out_graph->types[off++] = MD_TOPO_MINIMUM;
    for (uint32_t i = 0; i < num_join_saddles;  i++) out_graph->types[off++] = MD_TOPO_JOIN_SADDLE;

    if (num_edges > 0) {
        out_graph->edges = (md_topo_edge_t*)md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        MEMCPY(out_graph->edges, md_gpu_buffer_cpu_ptr(ctx->staging_edges), num_edges * sizeof(md_topo_edge_t));
    }

#if DEBUG
    {
        uint32_t n = num_vertices < 8u ? num_vertices : 8u;
        MD_LOG_DEBUG("[topo] First %u vertices:", n);
        for (uint32_t i = 0; i < n; i++) {
            MD_LOG_DEBUG("  vertex[%u]: type=%u pos=(%.4f, %.4f, %.4f)  val=%.6f",
                i, (uint32_t)out_graph->types[i],
                out_graph->vertices[i].x, out_graph->vertices[i].y,
                out_graph->vertices[i].z, out_graph->vertices[i].value);
        }
    }
#endif

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
        md_allocator_i* temp_alloc = md_get_heap_allocator();
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
        md_free(temp_alloc, vol_data,  num_points * sizeof(float));
        md_free(temp_alloc, asc_data,  num_points * sizeof(uint32_t));
        md_free(temp_alloc, desc_data, num_points * sizeof(uint32_t));
    }
#endif
    
    // === Step 3: Identify critical points ===
    GLuint critical_prog = get_critical_points_program();
    if (!critical_prog) goto cleanup;

	md_gl_debug_push("Critical Points");
    
    types_buf = create_buffer(num_points * sizeof(int), NULL, GL_DYNAMIC_COPY);
    
    // Counts buffer: maximaCount, splitSaddleCount, minimaCount, joinSaddleCount
    uint32_t counts_init[4] = {0};
    counts_buf = create_buffer(sizeof(counts_init), counts_init, GL_DYNAMIC_COPY);
    
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
    
    // Read back counts
    uint32_t counts[4] = {0};
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, counts_buf);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(counts), counts);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    
    uint32_t num_maxima = counts[0];
    uint32_t num_split_saddles = counts[1];
    uint32_t num_minima = counts[2];
    uint32_t num_join_saddles = counts[3];
    uint32_t num_vertices = num_maxima + num_split_saddles + num_minima + num_join_saddles;
    uint32_t num_edges = 8 * (num_split_saddles + num_join_saddles); // We estimate this from split saddles + join saddles for allocation, this will be updated later to the exact count
    
    MD_LOG_DEBUG("Topology: %u maxima, %u split saddles, %u minima, %u join saddles (total: %u)",
                num_maxima, num_split_saddles, num_minima, num_join_saddles, num_vertices);
    
    // === Step 4: Compact critical point indices into ordered array ===
    GLuint compaction_prog = get_critical_point_compaction_program();
    if (!compaction_prog) goto cleanup;

	md_gl_debug_push("Critical Point Compaction");

    // Allocate unified indices buffer (packed: maxima, split saddles, minima, join saddles)
    if (num_vertices > 0) {
        indices_buf = create_buffer(num_vertices * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
        if (!indices_buf) {
            MD_LOG_ERROR("Failed to allocate indices buffer");
            goto cleanup;
        }
    } else {
        indices_buf = 0;
    }
    
    // These are the write offsets into the unified indices buffer
    uint32_t counters_init[4] = {0, num_maxima, num_maxima + num_split_saddles, num_maxima + num_split_saddles + num_minima};
    counter_buf = create_buffer(sizeof(counters_init), counters_init, GL_DYNAMIC_COPY);
    
    glUseProgram(compaction_prog);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, indices_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, counter_buf);

    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	md_gl_debug_pop(); // Critical Point Compaction

#if 0
    // Print out all of the critical point indices found (sorted)
    #if DEBUG
    {
        md_allocator_i* temp_alloc = md_get_heap_allocator();
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
            md_free(temp_alloc, data, num_vertices * sizeof(uint32_t));
        }
    }
    #endif
#endif
    
    // === Step 5: Extract vertices and edges using GPU shader ===
    GLuint extraction_prog = get_vertex_edge_extraction_program();
    if (!extraction_prog) goto cleanup;

	md_gl_debug_push("Graph Extraction");
    
    // Create type counts buffer
    struct {
        uint32_t num_maxima;
        uint32_t num_split_saddles;
        uint32_t num_minima;
        uint32_t num_join_saddles;
        uint32_t num_vertices;
    } type_counts = {num_maxima, num_split_saddles, num_minima, num_join_saddles, num_vertices};
    GLuint type_counts_buf = create_buffer(sizeof(type_counts), &type_counts, GL_STATIC_DRAW);
    
    // Create output buffers for vertices and edges
    GLuint vert_buf = create_buffer(num_vertices * sizeof(md_topo_vert_t), NULL, GL_DYNAMIC_COPY);
    GLuint edge_buf = create_buffer(num_edges    * sizeof(md_topo_edge_t), NULL, GL_DYNAMIC_COPY);
    
    uint32_t edge_count_init = 0;
    GLuint edge_count_buf = create_buffer(sizeof(uint32_t), &edge_count_init, GL_DYNAMIC_COPY);
    
    glUseProgram(extraction_prog);
    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R32F);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, indices_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, type_counts_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, descending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, vert_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, edge_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 7, edge_count_buf);
    // Bind types buffer as voxel_id -> vertex_index mapping (binding = 8)
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, types_buf);
    
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
    delete_buffer(type_counts_buf);
    delete_buffer(vert_buf);
    delete_buffer(edge_buf);
    delete_buffer(edge_count_buf);
    
    // Fill output structure
    MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
    out_graph->num_vertices = num_vertices;
    out_graph->vertices = vertices;
    out_graph->num_maxima = num_maxima;
    out_graph->num_split_saddles = num_split_saddles;
    out_graph->num_minima = num_minima;
    out_graph->num_join_saddles = num_join_saddles;
    out_graph->num_edges = num_edges;
    out_graph->edges = edges;
    out_graph->alloc = alloc;
    
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

    md_allocator_i* temp_alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

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
        md_arena_allocator_destroy(temp_alloc);
        return;
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

    md_arena_allocator_destroy(temp_alloc);
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

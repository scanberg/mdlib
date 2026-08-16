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

/* ---------------------------------------------------------------------------
   Kernels. One per Slang file; md_gpu needs no reflection or binding tables,
   so these are just the compiled blobs.
   --------------------------------------------------------------------------- */
static md_gpu_kernel_t k_bidirectional_manifold    = NULL;
static md_gpu_kernel_t k_path_compression          = NULL;
static md_gpu_kernel_t k_critical_points           = NULL;
static md_gpu_kernel_t k_critical_point_compaction = NULL;
static md_gpu_kernel_t k_vertex_edge_extraction    = NULL;

static md_gpu_kernel_t ensure_kernel(md_gpu_device_t device, md_gpu_kernel_t* slot,
                                     const void* code, size_t code_size, const char* name,
                                     uint32_t gx, uint32_t gy, uint32_t gz) {
    if (*slot) return *slot;
    md_gpu_kernel_desc_t desc = {0};
    desc.code          = code;
    desc.code_size     = code_size;
    desc.label         = name;
    desc.group_size[0] = gx;
    desc.group_size[1] = gy;
    desc.group_size[2] = gz;
    *slot = md_gpu_kernel_create(device, &desc);
    if (!*slot) MD_LOG_ERROR("md_topo: failed to create kernel '%s': %s", name, md_gpu_last_error());
    return *slot;
}

void md_topo_gpu_initialize(md_gpu_device_t device) {
    if (!device) return;
    ensure_kernel(device, &k_bidirectional_manifold,
                  md_shader_bidirectional_manifold_main_start,
                  md_shader_bidirectional_manifold_main_byte_size, "topo bidirectional_manifold", 8, 8, 8);
    ensure_kernel(device, &k_path_compression,
                  md_shader_path_compression_main_start,
                  md_shader_path_compression_main_byte_size, "topo path_compression", 8, 8, 8);
    ensure_kernel(device, &k_critical_points,
                  md_shader_critical_points_main_start,
                  md_shader_critical_points_main_byte_size, "topo critical_points", 8, 8, 8);
    ensure_kernel(device, &k_critical_point_compaction,
                  md_shader_critical_point_compaction_main_start,
                  md_shader_critical_point_compaction_main_byte_size, "topo critical_point_compaction", 8, 8, 8);
    ensure_kernel(device, &k_vertex_edge_extraction,
                  md_shader_vertex_edge_extraction_main_start,
                  md_shader_vertex_edge_extraction_main_byte_size, "topo vertex_edge_extraction", 64, 1, 1);
}

void md_topo_gpu_shutdown(void) {
    if (k_bidirectional_manifold)    { md_gpu_kernel_destroy(k_bidirectional_manifold);    k_bidirectional_manifold = NULL; }
    if (k_path_compression)          { md_gpu_kernel_destroy(k_path_compression);          k_path_compression = NULL; }
    if (k_critical_points)           { md_gpu_kernel_destroy(k_critical_points);           k_critical_points = NULL; }
    if (k_critical_point_compaction) { md_gpu_kernel_destroy(k_critical_point_compaction); k_critical_point_compaction = NULL; }
    if (k_vertex_edge_extraction)    { md_gpu_kernel_destroy(k_vertex_edge_extraction);    k_vertex_edge_extraction = NULL; }
}

// Meta buffer: shared by every kernel in the chain.
typedef struct {
    uint32_t vertex_count;  // total CP count (critical_points output)
    uint32_t edge_count;    // actual edge count (extraction output)
    uint32_t changed_read;  // path-compression convergence flag, read by shader
    uint32_t changed_write; // path-compression convergence flag, written by shader
    uint32_t counter;       // compaction write cursor
} topo_meta_t;

/* Argument structs, mirroring src/shaders/topo/*.slang. The leading float4x4
   is at offset 0, where SPIR-V and MSL agree; md_gpu_float4x4 keeps it that way
   if anything is ever inserted before it.
   tools/check_gpu_arg_layout.py verifies these against the compiled shaders. */
typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    uint64_t        ascending;
    uint64_t        descending;
    md_gpu_tex_t    vol_tex;
} topo_manifold_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    uint64_t        ascending;
    uint64_t        descending;
    uint64_t        meta;
} topo_path_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    uint64_t        ascending;
    uint64_t        descending;
    uint64_t        types;
    uint64_t        meta;
    md_gpu_tex_t    vol_tex;
} topo_critical_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    uint64_t        types;
    uint64_t        cp_indices;
    uint64_t        meta;
    uint64_t        voxel_to_vertex_idx;
    uint64_t        vertex_types;
} topo_compact_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    uint64_t        cp_indices;
    uint64_t        vertex_types;
    uint64_t        vertex_data;
    uint64_t        edges;
    uint64_t        ascending;
    uint64_t        descending;
    uint64_t        voxel_to_vertex_idx;
    uint64_t        meta;
    md_gpu_tex_t    vol_tex;
} topo_extract_args_t;

// Worst-case capacity ratios:
//   vert_cap = num_points / TOPO_VERT_RATIO  (1 CP per 8 voxels is very generous)
//   edge_cap = vert_cap * TOPO_EDGE_RATIO
#define TOPO_VERT_RATIO  8
#define TOPO_EDGE_RATIO  4
#define TOPO_VERT_CAP_MIN 64

struct md_topo_gpu_context {
    md_gpu_device_t device;
    /* One pool per memory kind: a pool maps 1:1 onto a VkDeviceMemory chain
       or an MTLHeap only if its kind is fixed. Both are owned by the context
       and torn down with it. */
    md_gpu_pool_t   dev_pool;      // device-local scratch and results
    md_gpu_pool_t   read_pool;     // host-readable mirrors
    md_gpu_stream_t stream;        // the stream the last record used
    uint32_t        num_points;
    uint32_t        dim[3];
    uint32_t        vert_cap;
    uint32_t        edge_cap;

    // Scratch, device-local. These used to be cached for the lifetime of the
    // context because allocating per call was too expensive; the pool makes
    // that a non-issue, but keeping them avoids re-zeroing bookkeeping.
    uint32_t*       ascending;
    uint32_t*       descending;
    uint32_t*       voxel_types;
    int32_t*        voxel_to_vert;
    topo_meta_t*    meta;
    uint32_t*       grid_args;     // 3 x uint32, indirect dispatch dimensions

    // Results, device-local, sized to worst case.
    uint32_t*       indices;
    float*          verts;         // float4 per vertex
    uint32_t*       types;
    uint32_t*       edges;         // uint2 per edge

    // Host-readable mirrors, filled at the end of md_topo_gpu_record.
    topo_meta_t*    host_meta;
    float*          host_verts;
    uint32_t*       host_types;
    uint32_t*       host_edges;
};

md_topo_gpu_context_t* md_topo_gpu_context_create(md_gpu_device_t device, uint32_t dim_x, uint32_t dim_y, uint32_t dim_z) {
    if (!device) return NULL;

    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)md_alloc(md_get_heap_allocator(), sizeof(struct md_topo_gpu_context));
    if (!ctx) return NULL;
    MEMSET(ctx, 0, sizeof(*ctx));

    ctx->device     = device;
    ctx->num_points = dim_x * dim_y * dim_z;
    ctx->dim[0] = dim_x; ctx->dim[1] = dim_y; ctx->dim[2] = dim_z;

    uint32_t vert_cap = ctx->num_points / TOPO_VERT_RATIO;
    if (vert_cap < TOPO_VERT_CAP_MIN) vert_cap = TOPO_VERT_CAP_MIN;
    ctx->vert_cap = vert_cap;
    ctx->edge_cap = vert_cap * TOPO_EDGE_RATIO;

    md_gpu_pool_desc_t pd = {0};
    pd.flags = MD_GPU_MEM_DEVICE;
    pd.label = "md_topo device";
    ctx->dev_pool = md_gpu_pool_create(device, &pd);

    pd.flags = MD_GPU_MEM_HOST_READ;
    pd.label = "md_topo readback";
    ctx->read_pool = md_gpu_pool_create(device, &pd);

    if (!ctx->dev_pool || !ctx->read_pool) {
        md_topo_gpu_context_destroy((md_topo_gpu_context_t*)ctx);
        return NULL;
    }

    md_gpu_stream_t s = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    const size_t voxel_bytes = (size_t)ctx->num_points * sizeof(uint32_t);

    ctx->ascending     = (uint32_t*)md_gpu_malloc(ctx->dev_pool, voxel_bytes, s);
    ctx->descending    = (uint32_t*)md_gpu_malloc(ctx->dev_pool, voxel_bytes, s);
    ctx->voxel_types   = (uint32_t*)md_gpu_malloc(ctx->dev_pool, voxel_bytes, s);
    ctx->voxel_to_vert = (int32_t*) md_gpu_malloc(ctx->dev_pool, voxel_bytes, s);
    ctx->meta          = (topo_meta_t*)md_gpu_malloc(ctx->dev_pool, sizeof(topo_meta_t), s);
    ctx->grid_args     = (uint32_t*)md_gpu_malloc(ctx->dev_pool, 3 * sizeof(uint32_t), s);

    ctx->indices = (uint32_t*)md_gpu_malloc(ctx->dev_pool, vert_cap * sizeof(uint32_t), s);
    ctx->verts   = (float*)   md_gpu_malloc(ctx->dev_pool, vert_cap * 4 * sizeof(float), s);
    ctx->types   = (uint32_t*)md_gpu_malloc(ctx->dev_pool, vert_cap * sizeof(uint32_t), s);
    ctx->edges   = (uint32_t*)md_gpu_malloc(ctx->dev_pool, ctx->edge_cap * 2 * sizeof(uint32_t), s);

    ctx->host_meta  = (topo_meta_t*)md_gpu_malloc(ctx->read_pool, sizeof(topo_meta_t), s);
    ctx->host_verts = (float*)   md_gpu_malloc(ctx->read_pool, vert_cap * 4 * sizeof(float), s);
    ctx->host_types = (uint32_t*)md_gpu_malloc(ctx->read_pool, vert_cap * sizeof(uint32_t), s);
    ctx->host_edges = (uint32_t*)md_gpu_malloc(ctx->read_pool, ctx->edge_cap * 2 * sizeof(uint32_t), s);

    if (!ctx->ascending || !ctx->descending || !ctx->voxel_types || !ctx->voxel_to_vert ||
        !ctx->meta || !ctx->grid_args || !ctx->indices || !ctx->verts || !ctx->types || !ctx->edges ||
        !ctx->host_meta || !ctx->host_verts || !ctx->host_types || !ctx->host_edges) {
        MD_LOG_ERROR("md_topo_gpu_context_create: allocation failed: %s", md_gpu_last_error());
        md_topo_gpu_context_destroy((md_topo_gpu_context_t*)ctx);
        return NULL;
    }
    return (md_topo_gpu_context_t*)ctx;
}

void md_topo_gpu_context_destroy(md_topo_gpu_context_t* context) {
    if (!context) return;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;
    if (ctx->dev_pool)  md_gpu_pool_destroy(ctx->dev_pool);   // waits for idle
    if (ctx->read_pool) md_gpu_pool_destroy(ctx->read_pool);
    md_free(md_get_heap_allocator(), ctx, sizeof(*ctx));
}

void md_topo_gpu_record(md_gpu_stream_t stream, md_topo_gpu_context_t* context,
                        md_gpu_tex_t volume, const md_grid_t* grid, float scalar_threshold) {
    if (!stream || !context || !grid) return;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    if (!k_bidirectional_manifold || !k_path_compression || !k_critical_points ||
        !k_critical_point_compaction || !k_vertex_edge_extraction) {
        MD_LOG_ERROR("md_topo_gpu_record: kernels not initialized");
        return;
    }
    ctx->stream = stream;

    md_gpu_float4x4 i2w;
    index_to_world_matrix((float(*)[4])i2w.m, grid);
    const md_gpu_uint4 dims = { grid->dim[0], grid->dim[1], grid->dim[2], 0 };

    const md_gpu_grid_t wg = md_gpu_grid(DIV_UP(grid->dim[0], 8), DIV_UP(grid->dim[1], 8), DIV_UP(grid->dim[2], 8));
    const size_t voxel_bytes = (size_t)ctx->num_points * sizeof(uint32_t);

    /* No barriers anywhere below: everything issued into `stream` runs in
       order and observes the previous step's writes. */

    // Per-call resets.
    md_gpu_memset_async(ctx->voxel_to_vert, 0xFF, voxel_bytes, stream);   // = -1
    md_gpu_memset_async(ctx->meta, 0, sizeof(topo_meta_t), stream);
    md_gpu_memset_async(ctx->types, 0, ctx->vert_cap * sizeof(uint32_t), stream);

    // Step 1: bidirectional manifold.
    topo_manifold_args_t ma = {0};
    ma.index_to_world = i2w; ma.dims = dims; ma.scalar_threshold = scalar_threshold;
    ma.ascending  = (uint64_t)(uintptr_t)ctx->ascending;
    ma.descending = (uint64_t)(uintptr_t)ctx->descending;
    ma.vol_tex    = volume;
    md_gpu_launch(stream, k_bidirectional_manifold, wg, &ma, sizeof(ma));

    // Step 2: path compression, iterated with a GPU-side early-exit flag.
    uint32_t iterations = 0;
    {
        uint32_t max_dim = (uint32_t)MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
        while (max_dim > (1U << iterations)) iterations++;
        iterations *= 2;
    }
    md_gpu_memset_async((char*)ctx->meta + offsetof(topo_meta_t, changed_read), 0xFF, 4, stream);

    topo_path_args_t pa = {0};
    pa.index_to_world = i2w; pa.dims = dims; pa.scalar_threshold = scalar_threshold;
    pa.ascending  = (uint64_t)(uintptr_t)ctx->ascending;
    pa.descending = (uint64_t)(uintptr_t)ctx->descending;
    pa.meta       = (uint64_t)(uintptr_t)ctx->meta;
    for (uint32_t i = 0; i < iterations; ++i) {
        md_gpu_launch(stream, k_path_compression, wg, &pa, sizeof(pa));
        md_gpu_memcpy_async((char*)ctx->meta + offsetof(topo_meta_t, changed_read),
                            (char*)ctx->meta + offsetof(topo_meta_t, changed_write), 4, stream);
        md_gpu_memset_async((char*)ctx->meta + offsetof(topo_meta_t, changed_write), 0, 4, stream);
    }

    // Step 3: critical-point detection.
    topo_critical_args_t ca = {0};
    ca.index_to_world = i2w; ca.dims = dims; ca.scalar_threshold = scalar_threshold;
    ca.ascending  = (uint64_t)(uintptr_t)ctx->ascending;
    ca.descending = (uint64_t)(uintptr_t)ctx->descending;
    ca.types      = (uint64_t)(uintptr_t)ctx->voxel_types;
    ca.meta       = (uint64_t)(uintptr_t)ctx->meta;
    ca.vol_tex    = volume;
    md_gpu_launch(stream, k_critical_points, wg, &ca, sizeof(ca));

    // Step 4: compaction.
    topo_compact_args_t ka = {0};
    ka.index_to_world = i2w; ka.dims = dims; ka.scalar_threshold = scalar_threshold;
    ka.types               = (uint64_t)(uintptr_t)ctx->voxel_types;
    ka.cp_indices          = (uint64_t)(uintptr_t)ctx->indices;
    ka.meta                = (uint64_t)(uintptr_t)ctx->meta;
    ka.voxel_to_vertex_idx = (uint64_t)(uintptr_t)ctx->voxel_to_vert;
    ka.vertex_types        = (uint64_t)(uintptr_t)ctx->types;
    md_gpu_launch(stream, k_critical_point_compaction, wg, &ka, sizeof(ka));

    /* Step 5: vertex + edge extraction. This used to be dispatched at
       worst-case capacity with a per-thread early-out; the vertex count that
       step 4 produced now drives an indirect launch, so only the compacted
       vertices are covered and the count never reaches the CPU. */
    topo_extract_args_t ea = {0};
    ea.index_to_world = i2w; ea.dims = dims; ea.scalar_threshold = scalar_threshold;
    ea.cp_indices          = (uint64_t)(uintptr_t)ctx->indices;
    ea.vertex_types        = (uint64_t)(uintptr_t)ctx->types;
    ea.vertex_data         = (uint64_t)(uintptr_t)ctx->verts;
    ea.edges               = (uint64_t)(uintptr_t)ctx->edges;
    ea.ascending           = (uint64_t)(uintptr_t)ctx->ascending;
    ea.descending          = (uint64_t)(uintptr_t)ctx->descending;
    ea.voxel_to_vertex_idx = (uint64_t)(uintptr_t)ctx->voxel_to_vert;
    ea.meta                = (uint64_t)(uintptr_t)ctx->meta;
    ea.vol_tex             = volume;

    const uint32_t local[3] = {64, 1, 1};
    md_gpu_make_grid(stream, ctx->grid_args,
                     (char*)ctx->meta + offsetof(topo_meta_t, vertex_count), local);
    md_gpu_launch_indirect(stream, k_vertex_edge_extraction, ctx->grid_args, &ea, sizeof(ea));

    // Mirror the results where the CPU can read them.
    md_gpu_memcpy_async(ctx->host_meta,  ctx->meta,  sizeof(topo_meta_t), stream);
    md_gpu_memcpy_async(ctx->host_verts, ctx->verts, ctx->vert_cap * 4 * sizeof(float), stream);
    md_gpu_memcpy_async(ctx->host_types, ctx->types, ctx->vert_cap * sizeof(uint32_t), stream);
    md_gpu_memcpy_async(ctx->host_edges, ctx->edges, ctx->edge_cap * 2 * sizeof(uint32_t), stream);
}

bool md_topo_gpu_context_extract(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context) {
    if (!context || !out_graph) return false;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    ASSERT(out_graph->alloc);

    const topo_meta_t* meta = (const topo_meta_t*)md_gpu_host_ptr(ctx->host_meta);
    if (!meta) return false;

    const uint32_t num_vertices = meta->vertex_count;
    const uint32_t num_edges    = meta->edge_count;

    if (meta->changed_read != 0) {
        MD_LOG_ERROR("md_topo_gpu_context_extract: path compression did not fully converge - results may be approximate");
    }

    MD_LOG_DEBUG("Topology: %u vertices, %u edges", num_vertices, num_edges);

    if (num_vertices == 0) return false;

    md_allocator_i* alloc = out_graph->alloc;
    MEMSET(out_graph, 0, sizeof(*out_graph));
    out_graph->alloc        = alloc;
    out_graph->num_vertices = num_vertices;
    out_graph->num_edges    = num_edges;

    out_graph->vertices = (md_topo_vert_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
    out_graph->types    = (md_topo_critical_point_type_t*)md_alloc(alloc, num_vertices * sizeof(md_topo_critical_point_type_t));

    const float* vp = (const float*)md_gpu_host_ptr(ctx->host_verts);
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->vertices[i].x     = vp[i * 4 + 0];
        out_graph->vertices[i].y     = vp[i * 4 + 1];
        out_graph->vertices[i].z     = vp[i * 4 + 2];
        out_graph->vertices[i].value = vp[i * 4 + 3];
    }

    const uint32_t* tp = (const uint32_t*)md_gpu_host_ptr(ctx->host_types);
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->types[i] = (md_topo_critical_point_type_t)tp[i];
    }

    if (num_edges > 0) {
        out_graph->edges = (md_topo_edge_t*)md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        MEMCPY(out_graph->edges, md_gpu_host_ptr(ctx->host_edges), num_edges * sizeof(md_topo_edge_t));
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
        md_temp_scope_t temp_scope = md_temp_begin();
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
        md_temp_scope_t temp_scope = md_temp_begin();
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

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
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

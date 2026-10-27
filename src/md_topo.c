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
   Kernels. One per Slang file. The generated *_kernel() descriptors carry the
   group size and argument-struct size, so nothing here repeats [numthreads].
   --------------------------------------------------------------------------- */
static md_gpu_kernel_t k_bidirectional_manifold    = NULL;
static md_gpu_kernel_t k_path_compression          = NULL;
static md_gpu_kernel_t k_critical_points           = NULL;
static md_gpu_kernel_t k_critical_point_compaction = NULL;
static md_gpu_kernel_t k_vertex_edge_extraction    = NULL;

static void ensure_kernel(md_gpu_device_t device, md_gpu_kernel_t* slot, md_gpu_kernel_desc_t desc) {
    if (*slot) return;
    *slot = md_gpu_kernel_create(device, &desc);
    if (!*slot) MD_LOG_ERROR("md_topo: failed to create kernel '%s': %s", desc.label, md_gpu_last_error());
}

void md_topo_gpu_initialize(md_gpu_device_t device) {
    if (!device) return;
    ensure_kernel(device, &k_bidirectional_manifold,    md_shader_bidirectional_manifold_main_kernel());
    ensure_kernel(device, &k_path_compression,          md_shader_path_compression_main_kernel());
    ensure_kernel(device, &k_critical_points,           md_shader_critical_points_main_kernel());
    ensure_kernel(device, &k_critical_point_compaction, md_shader_critical_point_compaction_main_kernel());
    ensure_kernel(device, &k_vertex_edge_extraction,    md_shader_vertex_edge_extraction_main_kernel());
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

/* Argument structs, mirroring the kernels in src/shaders/topo/. The leading float4x4
   is at offset 0, where SPIR-V and MSL agree; md_gpu_float4x4 keeps it that way
   if anything is ever inserted before it.
   tools/check_gpu_arg_layout.py verifies these against the compiled shaders. */
typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    md_gpu_addr_t   ascending;
    md_gpu_addr_t   descending;
    md_gpu_storage_tex_t vol_tex;
} topo_manifold_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    md_gpu_addr_t   ascending;
    md_gpu_addr_t   descending;
    md_gpu_addr_t   meta;
} topo_path_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    md_gpu_addr_t   ascending;
    md_gpu_addr_t   descending;
    md_gpu_addr_t   types;
    md_gpu_addr_t   meta;
    md_gpu_storage_tex_t vol_tex;
} topo_critical_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    md_gpu_addr_t   types;
    md_gpu_addr_t   cp_indices;
    md_gpu_addr_t   meta;
    md_gpu_addr_t   voxel_to_vertex_idx;
    md_gpu_addr_t   vertex_types;
} topo_compact_args_t;

typedef struct {
    md_gpu_float4x4 index_to_world;
    md_gpu_uint4    dims;
    float           scalar_threshold;
    md_gpu_addr_t   cp_indices;
    md_gpu_addr_t   vertex_types;
    md_gpu_addr_t   vertex_data;
    md_gpu_addr_t   edges;
    md_gpu_addr_t   ascending;
    md_gpu_addr_t   descending;
    md_gpu_addr_t   voxel_to_vertex_idx;
    md_gpu_addr_t   meta;
    md_gpu_storage_tex_t vol_tex;
} topo_extract_args_t;

// Worst-case capacity ratios:
//   vert_cap = num_points / TOPO_VERT_RATIO  (1 CP per 8 voxels is very generous)
//   edge_cap = vert_cap * TOPO_EDGE_RATIO
#define TOPO_VERT_RATIO  8
#define TOPO_EDGE_RATIO  4
#define TOPO_VERT_CAP_MIN 64

struct md_topo_gpu_context {
    md_gpu_device_t device;
    md_gpu_stream_t stream;        // the stream the last record used; frees go here
    uint32_t        num_points;
    uint32_t        dim[3];
    uint32_t        vert_cap;
    uint32_t        edge_cap;

    // Scratch, device-local, allocated once for the lifetime of the context.
    md_gpu_addr_t   ascending;
    md_gpu_addr_t   descending;
    md_gpu_addr_t   voxel_types;
    md_gpu_addr_t   voxel_to_vert;
    md_gpu_addr_t   meta;          // topo_meta_t
    md_gpu_addr_t   grid_args;     // 3 x uint32, indirect dispatch dimensions

    // Results, device-local, sized to worst case.
    md_gpu_addr_t   indices;
    md_gpu_addr_t   verts;         // float4 per vertex
    md_gpu_addr_t   types;
    md_gpu_addr_t   edges;         // uint2 per edge

    // Host-readable mirrors, filled at the end of md_topo_gpu_record.
    md_gpu_mem_t    host_meta;     // topo_meta_t
    md_gpu_mem_t    host_verts;
    md_gpu_mem_t    host_types;
    md_gpu_mem_t    host_edges;
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

    md_gpu_stream_t s = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    ctx->stream = s;
    const size_t voxel_bytes = (size_t)ctx->num_points * sizeof(uint32_t);

    ctx->ascending     = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, voxel_bytes).gpu;
    ctx->descending    = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, voxel_bytes).gpu;
    ctx->voxel_types   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, voxel_bytes).gpu;
    ctx->voxel_to_vert = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, voxel_bytes).gpu;
    ctx->meta          = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, sizeof(topo_meta_t)).gpu;
    ctx->grid_args     = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, 3 * sizeof(uint32_t)).gpu;

    ctx->indices = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, vert_cap * sizeof(uint32_t)).gpu;
    ctx->verts   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, vert_cap * 4 * sizeof(float)).gpu;
    ctx->types   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, vert_cap * sizeof(uint32_t)).gpu;
    ctx->edges   = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, ctx->edge_cap * 2 * sizeof(uint32_t)).gpu;

    ctx->host_meta  = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, sizeof(topo_meta_t));
    ctx->host_verts = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, vert_cap * 4 * sizeof(float));
    ctx->host_types = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, vert_cap * sizeof(uint32_t));
    ctx->host_edges = md_gpu_malloc(s, MD_GPU_MEM_HOST_READ, ctx->edge_cap * 2 * sizeof(uint32_t));

    if (!ctx->ascending || !ctx->descending || !ctx->voxel_types || !ctx->voxel_to_vert ||
        !ctx->meta || !ctx->grid_args || !ctx->indices || !ctx->verts || !ctx->types || !ctx->edges ||
        !ctx->host_meta.cpu || !ctx->host_verts.cpu || !ctx->host_types.cpu || !ctx->host_edges.cpu) {
        MD_LOG_ERROR("md_topo_gpu_context_create: allocation failed: %s", md_gpu_last_error());
        md_topo_gpu_context_destroy((md_topo_gpu_context_t*)ctx);
        return NULL;
    }
    return (md_topo_gpu_context_t*)ctx;
}

void md_topo_gpu_context_destroy(md_topo_gpu_context_t* context) {
    if (!context) return;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;
    // Non-blocking: freed at the end of the stream the context last recorded
    // into, so the memory is reused only once that work has completed.
    md_gpu_stream_t s = ctx->stream;
    if (s) {
        const md_gpu_addr_t all[] = {
            ctx->ascending, ctx->descending, ctx->voxel_types, ctx->voxel_to_vert, ctx->meta, ctx->grid_args,
            ctx->indices, ctx->verts, ctx->types, ctx->edges,
            ctx->host_meta.gpu, ctx->host_verts.gpu, ctx->host_types.gpu, ctx->host_edges.gpu,
        };
        for (size_t i = 0; i < sizeof(all) / sizeof(all[0]); ++i) md_gpu_free(s, all[i]);
    }
    md_free(md_get_heap_allocator(), ctx, sizeof(*ctx));
}

void md_topo_gpu_record(md_gpu_stream_t stream, md_topo_gpu_context_t* context,
                        md_gpu_texture_t volume, const md_grid_t* grid, float scalar_threshold) {
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

    const md_gpu_storage_tex_t vol = md_gpu_texture_storage(volume, 0);
    if (!vol.handle) {
        MD_LOG_ERROR("md_topo_gpu_record: the volume texture needs MD_GPU_TEX_STORAGE usage");
        return;
    }
    const size_t voxel_bytes = (size_t)ctx->num_points * sizeof(uint32_t);

    /* No barriers anywhere below: everything issued into `stream` runs in
       order and observes the previous step's writes. */

    // Per-call resets.
    md_gpu_memset(stream, ctx->voxel_to_vert, 0xFF, voxel_bytes);   // = -1
    md_gpu_memset(stream, ctx->meta, 0, sizeof(topo_meta_t));
    md_gpu_memset(stream, ctx->types, 0, ctx->vert_cap * sizeof(uint32_t));

    // Step 1: bidirectional manifold.
    topo_manifold_args_t ma = {0};
    ma.index_to_world = i2w; ma.dims = dims; ma.scalar_threshold = scalar_threshold;
    ma.ascending  = ctx->ascending;
    ma.descending = ctx->descending;
    ma.vol_tex    = vol;
    md_gpu_launch(stream, k_bidirectional_manifold, md_gpu_grid_for(k_bidirectional_manifold, dims.x, dims.y, dims.z), &ma, sizeof(ma));

    // Step 2: path compression, iterated with a GPU-side early-exit flag.
    uint32_t iterations = 0;
    {
        uint32_t max_dim = (uint32_t)MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
        while (max_dim > (1U << iterations)) iterations++;
        iterations *= 2;
    }
    md_gpu_memset(stream, ctx->meta + offsetof(topo_meta_t, changed_read), 0xFF, 4);

    topo_path_args_t pa = {0};
    pa.index_to_world = i2w; pa.dims = dims; pa.scalar_threshold = scalar_threshold;
    pa.ascending  = ctx->ascending;
    pa.descending = ctx->descending;
    pa.meta       = ctx->meta;
    const md_gpu_grid_t path_grid = md_gpu_grid_for(k_path_compression, dims.x, dims.y, dims.z);
    for (uint32_t i = 0; i < iterations; ++i) {
        md_gpu_launch(stream, k_path_compression, path_grid, &pa, sizeof(pa));
        md_gpu_copy(stream, ctx->meta + offsetof(topo_meta_t, changed_read),
                            ctx->meta + offsetof(topo_meta_t, changed_write), 4);
        md_gpu_memset(stream, ctx->meta + offsetof(topo_meta_t, changed_write), 0, 4);
    }

    // Step 3: critical-point detection.
    topo_critical_args_t ca = {0};
    ca.index_to_world = i2w; ca.dims = dims; ca.scalar_threshold = scalar_threshold;
    ca.ascending  = ctx->ascending;
    ca.descending = ctx->descending;
    ca.types      = ctx->voxel_types;
    ca.meta       = ctx->meta;
    ca.vol_tex    = vol;
    md_gpu_launch(stream, k_critical_points, md_gpu_grid_for(k_critical_points, dims.x, dims.y, dims.z), &ca, sizeof(ca));

    // Step 4: compaction.
    topo_compact_args_t ka = {0};
    ka.index_to_world = i2w; ka.dims = dims; ka.scalar_threshold = scalar_threshold;
    ka.types               = ctx->voxel_types;
    ka.cp_indices          = ctx->indices;
    ka.meta                = ctx->meta;
    ka.voxel_to_vertex_idx = ctx->voxel_to_vert;
    ka.vertex_types        = ctx->types;
    md_gpu_launch(stream, k_critical_point_compaction, md_gpu_grid_for(k_critical_point_compaction, dims.x, dims.y, dims.z), &ka, sizeof(ka));

    /* Step 5: vertex + edge extraction. This used to be dispatched at
       worst-case capacity with a per-thread early-out; the vertex count that
       step 4 produced now drives an indirect launch, so only the compacted
       vertices are covered and the count never reaches the CPU. */
    topo_extract_args_t ea = {0};
    ea.index_to_world = i2w; ea.dims = dims; ea.scalar_threshold = scalar_threshold;
    ea.cp_indices          = ctx->indices;
    ea.vertex_types        = ctx->types;
    ea.vertex_data         = ctx->verts;
    ea.edges               = ctx->edges;
    ea.ascending           = ctx->ascending;
    ea.descending          = ctx->descending;
    ea.voxel_to_vertex_idx = ctx->voxel_to_vert;
    ea.meta                = ctx->meta;
    ea.vol_tex             = vol;

    md_gpu_make_grid(stream, ctx->grid_args, ctx->meta + offsetof(topo_meta_t, vertex_count),
                     k_vertex_edge_extraction);
    md_gpu_launch_indirect(stream, k_vertex_edge_extraction, ctx->grid_args, &ea, sizeof(ea));

    // Mirror the results where the CPU can read them.
    md_gpu_copy(stream, ctx->host_meta.gpu,  ctx->meta,  sizeof(topo_meta_t));
    md_gpu_copy(stream, ctx->host_verts.gpu, ctx->verts, ctx->vert_cap * 4 * sizeof(float));
    md_gpu_copy(stream, ctx->host_types.gpu, ctx->types, ctx->vert_cap * sizeof(uint32_t));
    md_gpu_copy(stream, ctx->host_edges.gpu, ctx->edges, ctx->edge_cap * 2 * sizeof(uint32_t));
}

bool md_topo_gpu_context_extract(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context) {
    if (!context || !out_graph) return false;
    struct md_topo_gpu_context* ctx = (struct md_topo_gpu_context*)context;

    ASSERT(out_graph->alloc);

    const topo_meta_t* meta = (const topo_meta_t*)ctx->host_meta.cpu;
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

    const float* vp = (const float*)ctx->host_verts.cpu;
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->vertices[i].x     = vp[i * 4 + 0];
        out_graph->vertices[i].y     = vp[i * 4 + 1];
        out_graph->vertices[i].z     = vp[i * 4 + 2];
        out_graph->vertices[i].value = vp[i * 4 + 3];
    }

    const uint32_t* tp = (const uint32_t*)ctx->host_types.cpu;
    for (uint32_t i = 0; i < num_vertices; i++) {
        out_graph->types[i] = (md_topo_critical_point_type_t)tp[i];
    }

    if (num_edges > 0) {
        out_graph->edges = (md_topo_edge_t*)md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        MEMCPY(out_graph->edges, ctx->host_edges.cpu, num_edges * sizeof(md_topo_edge_t));
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

// =====================================================================================================
// Certified critical points of a GTO electron density (CPU reference)
// =====================================================================================================
//
// Finds every non-degenerate critical point of
//
//     rho(r) = sum_{mu,nu} D_{mu nu} phi_mu(r) phi_nu(r)
//
// with rho >= rho_min, without sampling the field on a grid. Space is partitioned by an octree of
// cubes and every cube is PROVEN to contain either no critical point, or exactly one (which Newton
// from the cube centre then converges to). Cubes that reach h_min without either proof are reported
// as unresolved (genuinely near-degenerate topology), grouped into clusters, each with the Brouwer
// degree of grad rho on its boundary (0: a cancelling pair or nothing, +-1: at least one CP).
//
// The proofs need rigorous enclosures of grad rho and of the Hessian over a cube. They are built from
// the structure of the basis:
//   * every AO derivative is a product of 1D factors d^a/dx^a [x^i exp(-alpha x^2)] whose supremum over
//     an interval has a closed form (each monomial term peaks at |x| = sqrt(m / 2 alpha)),
//   * each AO is expanded to second order at the cube centre with a rigorous per-AO remainder, and the
//     density matrix only ever multiplies exact centre data, so its (large) core/valence cancellations
//     survive. Only a sixth order term ever sees |D|.
//
// Tests per cube (centre c, half-width h; g, A, T = grad, Hessian, third derivatives of rho at c):
//   exclusion : grad rho(c+d) in g + A d + 1/2 T[d,d] +- r. No zero if, along some direction w,
//               |w.g| > sum_j |(A w)_j| h + sum_k |w_k| (q_k + r_k), or if the Newton point lies far
//               enough outside the cube.
//   Krawczyk  : K = c - A^-1 g + (I - A^-1 [H]) (X - c), [H] = A +- dH. K inside X => exactly one zero.
//               K disjoint from X => no zero. Run on a 1.5x inflated cube so roots on faces are caught;
//               a root is accepted by any cube that (closed, plus 1e-9) contains it and deduplicated.
// Contributions of AOs screened away from a cube are bounded, so screening does not break rigour: those
// of the shells between the cube's local radius and their far radius by their tails at the cube's
// distance, charged to the cube; all the others by one global term. The search domain is derived from
// the same tail bounds (rho <= rho_min outside).
//
// Separatrices (bond paths BCP -> 2 maxima, ring lines RCP -> CCP or infinity) are traced from the
// saddles along their unique eigen-directions with adaptive Dormand-Prince RK45; they become the graph
// edges (from = saddle, to = extremum), exactly the edges the voxel pipeline produces.
//
// Rigour holds up to floating-point rounding (no outward rounding; the enclosures carry far more slack
// than rounding error, and mdlib builds with -ffast-math, so no infinities are used as sentinels).
// Units: Bohr. Multithreaded over cubes (fixed interleaved chunks, results merged in cube order), so the
// output is bit-identical for any thread count. Atomic core zones (a ~20-25% optimisation, see the design
// notes) are not part of this reference.

#include <md_gto.h>
#include <core/md_os.h>
#include <core/md_intrinsics.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define CPG_LMAX    MD_GTO_MAX_ANGULAR_MOMENTUM
#define CPG_KI      (CPG_LMAX + 8)
#define CPG_KA      7               // 1D derivative orders 0..6 (the GPU kernel's separable remainders)
#define CPG_KM      (CPG_KI + CPG_KA + 1)
#define CPG_NMI     35              // multi-indices of order 0..4
#define CPG_NV      20              // centre values: orders 0..3
#define CPG_NDV     20              // D-products: orders 0..3
#define CPG_MAXPRIM 64
#define CPG_MAX_THREADS 64
#define CPG_SEP_EXTRA 2             // exact centre terms per 1D factor beyond each remainder's order (cpg_shell_rem)
#define CPG_SEP_A   (4 + CPG_SEP_EXTRA + 1)  // 1D derivative orders 0 .. 4 + CPG_SEP_EXTRA
#define CPG_CHUNK   32              // cubes per work chunk
#define CPG_DTILE   16              // AOs per tile of D: the pairs that meet (D on the GPU)
#define CPG_RUNGS   21              // screening rungs: tail thresholds tau 10^(m/2), m = 0..20 (cpg_screen_rung)
#define CPG_SCREEN_TAIL 2.5e-4      // the local radius: the Hessian tail a region may carry, relative to rho_min (cpg_run_init)

typedef struct cpg_shell_t {
    double   A[3];
    double   radius;                // far screening radius: all derivatives (order <= 2) below tau beyond it
    double   radius_trace;          // the separatrices' points: the local radius (cpg_run_init)
    double   tw[3];                 // tail weights: beyond a radius where its AOs are below t, the shell adds at
                                    // most t tw[0] to rho, t tw[1] to grad, t tw[2] to the Hessian (cpg_run_init)
    uint32_t prim_offset;
    uint32_t num_prims;
    uint32_t ao_offset;
    int      l;
    int      ncart;
    int      cls;                   // radial class: the rung radii (ctx->rung_r2)
} cpg_shell_t;

typedef struct cpg_cp_t {
    double x[3];
    double rho;
    double ev[3];
    double evec[3][3];              // columns
    double h;                       // half-width of the owning cube
    int    type;
} cpg_cp_t;

typedef struct cpg_box_t {
    double c[3];
    double h;
} cpg_box_t;

typedef struct cpg_ctx_t {
    md_allocator_i* alloc;
    int nao;
    int nshell;
    cpg_shell_t* shell;
    double* alpha;
    double* coeff;
    int    (*ao_ijk)[3];
    double* ao_nrm;
    const double* D;
    double tail_rho, tail_g, tail_H;    // every shell beyond its far radius, at once (the global tails)
    double tau;
    // Screening by region: a shell is local to a cube or batch within rung_loc's radius; between that and
    // its far radius it is left out and its tail charged to the region (cpg_screen_rung); beyond the far
    // radius it is in the global tails. rung_tau[m] = tau 10^(m/2); rung_r2[cls * CPG_RUNGS + m] the squared
    // radius beyond which a shell of class cls is below rung_tau[m] (descending in m; m = 0: the far radius).
    int     rung_loc;
    int     ncls;
    double  rung_tau[CPG_RUNGS];
    double* rung_r2;
    double eps;
    double kappa[CPG_KI][CPG_KA][CPG_KM];
    double kap_sep[CPG_LMAX + 1][CPG_SEP_A][CPG_LMAX + CPG_SEP_A];  // d^a [x^i E] = sum_m kap_sep[i][a][m] alpha^((m-i+a)/2) x^m E
    int    mi[CPG_NMI][3];
    int    mi_idx[5][5][5];
    int    m1[3], m2[3][3], m3[3][3][3], m4[3][3][3][3];
} cpg_ctx_t;

// Per-thread scratch; the context above is shared and read-only once set up.
typedef struct cpg_scratch_t {
    int*    L;          // local AO list
    int*    shells;     // local shell list
    double* V;          // [n][CPG_NV]
    double* Dv;         // [n][CPG_NDV]
    double* E;          // per-AO remainder vectors [n][16]
    double* Q;          // |D| products [n][8]
    uint64_t trace_evals;   // separatrices: point evaluations of this worker, and their local AOs summed
    uint64_t trace_aos;
    uint64_t trace_block_aos;   // the AOs of the cached block per evaluation, summed (cpg_point_trace)
    // separatrices: the local AOs of a box around the path and their block of D, kept while the path stays
    // in the box (cpg_point_trace)
    double  tc_c[3];
    int     tc_n, tc_ns;
    bool    tc_valid;
    bool    tc_direct;  // this box goes through cpg_point (its block would be too large)
    int*    tc_L;       // [N]
    int*    tc_shells;  // [nshell]
    double* tc_B;       // [tc_n * tc_n] row-major: B[i n + j] = D[L[i]][L[j]]
    size_t  tc_cap;     // doubles allocated for tc_B
} cpg_scratch_t;

typedef struct cpg_eval_t {
    double rho, g[3], A[3][3], T[3][3][3];
    double r[3];        // remainder of grad beyond g + A d + 1/2 T[d,d]
    double quadT[3];    // bound of 1/2 T[d,d]
    double dH[3][3];    // Hessian enclosure half-width
    double rho_up;      // upper bound of rho over the cube
} cpg_eval_t;

// ----------------------------------------------------------------------------------------------- small math

static inline double cpg_powi(double x, int n) {
    double r = 1.0;
    while (n > 0) { if (n & 1) r *= x; x *= x; n >>= 1; }
    return r;
}

static inline double cpg_nrm(int i, int j, int k) {
    double d = 1.0;
    for (int n = 2 * i - 1; n > 1; n -= 2) d *= n;
    for (int n = 2 * j - 1; n > 1; n -= 2) d *= n;
    for (int n = 2 * k - 1; n > 1; n -= 2) d *= n;
    return 1.0 / sqrt(d);
}

static inline double cpg_mono_max(int a, int b, int c) {
    const int L = a + b + c;
    if (L == 0) return 1.0;
    double num = 1.0;
    if (a) num *= pow((double)a, a);
    if (b) num *= pow((double)b, b);
    if (c) num *= pow((double)c, c);
    return sqrt(num / pow((double)L, L));
}

// Jacobi eigen decomposition of a symmetric 3x3 matrix; eigenvalues ascending, eigenvectors in columns.
static void cpg_eigen_sym3(double out_val[3], double out_vec[3][3], double M[3][3]) {
    double a[3][3], v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    memcpy(a, M, sizeof(a));
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
        double scale = fabs(a[0][0]) + fabs(a[1][1]) + fabs(a[2][2]);
        if (off <= 1e-300 || off <= 1e-17 * scale) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (a[p][q] == 0.0) continue;
                double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    int idx[3] = {0, 1, 2};
    for (int i = 0; i < 2; ++i) for (int j = i + 1; j < 3; ++j) if (a[idx[j]][idx[j]] < a[idx[i]][idx[i]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    for (int i = 0; i < 3; ++i) {
        out_val[i] = a[idx[i]][idx[i]];
        for (int k = 0; k < 3; ++k) out_vec[k][i] = v[k][idx[i]];
    }
}

// A^-1 via the eigen decomposition. Returns false if A is (numerically) singular.
static bool cpg_inverse_sym3(double out[3][3], const double val[3], double vec[3][3]) {
    const double amax = fmax(fabs(val[0]), fabs(val[2]));
    for (int i = 0; i < 3; ++i) if (!(fabs(val[i]) > 1e-14 * amax) || amax == 0.0) return false;
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) {
        double s = 0.0;
        for (int k = 0; k < 3; ++k) s += vec[r][k] * vec[c][k] / val[k];
        out[r][c] = s;
    }
    return true;
}

// ----------------------------------------------------------------------------------------------- setup

static void cpg_init_tables(cpg_ctx_t* ctx) {
    memset(ctx->kappa, 0, sizeof(ctx->kappa));
    for (int i = 0; i < CPG_KI; ++i) ctx->kappa[i][0][i] = 1.0;
    for (int a = 0; a < CPG_KA - 1; ++a) {
        for (int i = 0; i < CPG_KI - 1 - a; ++i) {
            for (int m = 0; m < CPG_KM; ++m) {
                double v = -2.0 * ctx->kappa[i + 1][a][m];
                if (i >= 1) v += i * ctx->kappa[i - 1][a][m];
                ctx->kappa[i][a + 1][m] = v;
            }
        }
    }
    // the same expansion with the recurrence in a at fixed i (d/dx [x^m E] = m x^(m-1) E - 2 alpha x^(m+1) E),
    // to the higher orders the separable remainders need
    memset(ctx->kap_sep, 0, sizeof(ctx->kap_sep));
    for (int i = 0; i <= CPG_LMAX; ++i) {
        ctx->kap_sep[i][0][i] = 1.0;
        for (int a = 0; a + 1 < CPG_SEP_A; ++a) {
            for (int m = 0; m <= i + a + 1; ++m) {
                double v = m + 1 <= i + a ? (m + 1) * ctx->kap_sep[i][a][m + 1] : 0.0;
                if (m >= 1) v -= 2.0 * ctx->kap_sep[i][a][m - 1];
                ctx->kap_sep[i][a + 1][m] = v;
            }
        }
    }
    int n = 0;
    for (int o = 0; o <= 4; ++o) {
        for (int i = o; i >= 0; --i) {
            for (int j = o - i; j >= 0; --j) {
                const int k = o - i - j;
                ctx->mi[n][0] = i; ctx->mi[n][1] = j; ctx->mi[n][2] = k;
                ctx->mi_idx[i][j][k] = n++;
            }
        }
    }
    for (int a = 0; a < 3; ++a) {
        int e[3] = {0, 0, 0}; e[a]++;
        ctx->m1[a] = ctx->mi_idx[e[0]][e[1]][e[2]];
        for (int b = 0; b < 3; ++b) {
            int f[3] = {e[0], e[1], e[2]}; f[b]++;
            ctx->m2[a][b] = ctx->mi_idx[f[0]][f[1]][f[2]];
            for (int c = 0; c < 3; ++c) {
                int g[3] = {f[0], f[1], f[2]}; g[c]++;
                ctx->m3[a][b][c] = ctx->mi_idx[g[0]][g[1]][g[2]];
                for (int d = 0; d < 3; ++d) {
                    int q[3] = {g[0], g[1], g[2]}; q[d]++;
                    ctx->m4[a][b][c][d] = ctx->mi_idx[q[0]][q[1]][q[2]];
                }
            }
        }
    }
}

#define M1(a)       (ctx->m1[a])
#define M2(a,b)     (ctx->m2[a][b])
#define M3(a,b,c)   (ctx->m3[a][b][c])
#define M4(a,b,c,d) (ctx->m4[a][b][c][d])

// sup over |x - centre| >= d of |d^(p,q,w) [x^i y^j z^k exp(-alpha r^2)]|, by monomial terms.
static double cpg_tail_prim(const cpg_ctx_t* ctx, int i, int j, int k, int p, int q, int w, double alpha, double d) {
    double acc = 0.0;
    for (int m1 = 0; m1 <= i + p; ++m1) {
        const double k1 = ctx->kappa[i][p][m1];
        if (k1 == 0.0) continue;
        for (int m2 = 0; m2 <= j + q; ++m2) {
            const double k2 = ctx->kappa[j][q][m2];
            if (k2 == 0.0) continue;
            for (int m3 = 0; m3 <= k + w; ++m3) {
                const double k3 = ctx->kappa[k][w][m3];
                if (k3 == 0.0) continue;
                const int M = m1 + m2 + m3;
                const int apow = ((m1 - i + p) + (m2 - j + q) + (m3 - k + w)) / 2;
                const double ts = fmax(d, sqrt(M / (2.0 * alpha)));
                acc += fabs(k1 * k2 * k3) * cpg_powi(alpha, apow) * cpg_mono_max(m1, m2, m3) * cpg_powi(ts, M) * exp(-alpha * ts * ts);
            }
        }
    }
    return acc;
}

// Max over the shell's AOs and all multi-indices of order <= max_order of the tail bound at distance d.
static double cpg_tail_shell(const cpg_ctx_t* ctx, const cpg_shell_t* s, int max_order, double d) {
    double best = 0.0;
    for (int ci = 0; ci < s->ncart; ++ci) {
        const int* ijk = ctx->ao_ijk[s->ao_offset + ci];
        const double nrm = ctx->ao_nrm[s->ao_offset + ci];
        for (int o = 0; o <= max_order; ++o) {
            for (int p = o; p >= 0; --p) for (int q = o - p; q >= 0; --q) {
                const int w = o - p - q;
                double v = 0.0;
                for (uint32_t ip = 0; ip < s->num_prims; ++ip) {
                    v += fabs(ctx->coeff[s->prim_offset + ip]) * cpg_tail_prim(ctx, ijk[0], ijk[1], ijk[2], p, q, w, ctx->alpha[s->prim_offset + ip], d);
                }
                best = fmax(best, nrm * v);
            }
        }
    }
    return best;
}

// Same, but per AO and for exactly order o (max over multi-indices of that order).
static double cpg_tail_ao(const cpg_ctx_t* ctx, const cpg_shell_t* s, int ci, int o, double d) {
    const int* ijk = ctx->ao_ijk[s->ao_offset + ci];
    double best = 0.0;
    for (int p = o; p >= 0; --p) for (int q = o - p; q >= 0; --q) {
        const int w = o - p - q;
        double v = 0.0;
        for (uint32_t ip = 0; ip < s->num_prims; ++ip) {
            v += fabs(ctx->coeff[s->prim_offset + ip]) * cpg_tail_prim(ctx, ijk[0], ijk[1], ijk[2], p, q, w, ctx->alpha[s->prim_offset + ip], d);
        }
        best = fmax(best, v);
    }
    return best * ctx->ao_nrm[s->ao_offset + ci];
}

// ----------------------------------------------------------------------------------------------- AO evaluation

// 1D tables g[i][a] = d^a/dx^a [x^i exp(-alpha x^2)] for i <= imax - a.
static void cpg_tables_1d(double g[CPG_LMAX + 5][4], double x, double alpha, int imax, int amax) {
    const double E = exp(-alpha * x * x);
    double xp = 1.0;
    for (int i = 0; i <= imax; ++i) { g[i][0] = xp * E; xp *= x; }
    for (int a = 0; a < amax; ++a) {
        for (int i = 0; i < imax - a; ++i) {
            double v = -2.0 * alpha * g[i + 1][a];
            if (i >= 1) v += i * g[i - 1][a];
            g[i][a + 1] = v;
        }
    }
}

// Values of all multi-index derivatives up to 'deriv' (<= 3) of the shell's AOs at x. out[ci * stride + mi].
static void cpg_shell_eval(const cpg_ctx_t* ctx, const cpg_shell_t* s, const double x[3], int deriv, double* out, int stride) {
    const int nmi = deriv == 0 ? 1 : deriv == 1 ? 4 : deriv == 2 ? 10 : 20;
    const double dx = x[0] - s->A[0], dy = x[1] - s->A[1], dz = x[2] - s->A[2];
    for (int ci = 0; ci < s->ncart; ++ci) for (int m = 0; m < nmi; ++m) out[ci * stride + m] = 0.0;
    double gx[CPG_LMAX + 5][4], gy[CPG_LMAX + 5][4], gz[CPG_LMAX + 5][4];
    for (uint32_t ip = 0; ip < s->num_prims; ++ip) {
        const double al = ctx->alpha[s->prim_offset + ip];
        const double cf = ctx->coeff[s->prim_offset + ip];
        const int imax = s->l + deriv;
        cpg_tables_1d(gx, dx, al, imax, deriv);
        cpg_tables_1d(gy, dy, al, imax, deriv);
        cpg_tables_1d(gz, dz, al, imax, deriv);
        for (int ci = 0; ci < s->ncart; ++ci) {
            const int* ijk = ctx->ao_ijk[s->ao_offset + ci];
            const double c = cf * ctx->ao_nrm[s->ao_offset + ci];
            double* o = out + ci * stride;
            for (int m = 0; m < nmi; ++m) {
                const int* e = ctx->mi[m];
                o[m] += c * gx[ijk[0]][e[0]] * gy[ijk[1]][e[1]] * gz[ijk[2]][e[2]];
            }
        }
    }
}

// 1D factor of a primitive, X_{i,a}(x) = d^a/dx^a [x^i exp(-alpha x^2)], i <= l, on [xc - h, xc + h]:
// centre values for a <= CPG_SEP_A - 2 and interval sups for 3 <= a <= CPG_SEP_A - 1 (the orders the
// remainders below use). A sup adds up the monomial terms of the expansion, each at its own maximum on
// the interval: sup |t|^m exp(-alpha t^2) sits at t = sqrt(m / 2 alpha), where it is (m / 2 alpha)^(m/2)
// exp(-m/2), or at the nearer end of the interval.
static void cpg_sep_1d(const cpg_ctx_t* ctx, int l, double alpha, double xc, double h,
                       double val[CPG_LMAX + 1][CPG_SEP_A], double sup[CPG_LMAX + 1][CPG_SEP_A]) {
    static const double exp_neg_half[CPG_LMAX + CPG_SEP_A] = {
        1.0, 0.60653065971263342, 0.36787944117144233, 0.22313016014842982, 0.13533528323661270,
        0.08208499862389880, 0.04978706836786394, 0.03019738342231850, 0.01831563888873418,
        0.01110899653824231, 0.00673794699908547,
    };
    STATIC_ASSERT(CPG_LMAX + CPG_SEP_A <= 11, "exp(-m/2) table");
    // centre values by the recurrence g[i][a+1] = i g[i-1][a] - 2 alpha g[i+1][a]
    double g[CPG_LMAX + CPG_SEP_A][CPG_SEP_A - 1];
    const int imax = l + CPG_SEP_A - 2;
    const double E = exp(-alpha * xc * xc);
    double xp = 1.0;
    for (int i = 0; i <= imax; ++i) { g[i][0] = xp * E; xp *= xc; }
    for (int a = 0; a < CPG_SEP_A - 2; ++a) {
        for (int i = 0; i < imax - a; ++i) {
            double v = -2.0 * alpha * g[i + 1][a];
            if (i >= 1) v += i * g[i - 1][a];
            g[i][a + 1] = v;
        }
    }
    for (int i = 0; i <= l; ++i) for (int a = 0; a < CPG_SEP_A - 1; ++a) val[i][a] = g[i][a];
    // monomial sups Sm[m] = sup over the interval of |t|^m exp(-alpha t^2)
    const double x0 = xc - h, x1 = xc + h, ax0 = fabs(x0), ax1 = fabs(x1);
    const double tlo = (x0 <= 0.0 && x1 >= 0.0) ? 0.0 : fmin(ax0, ax1);
    const double thi = fmax(ax0, ax1);
    const double elo = exp(-alpha * tlo * tlo), ehi = exp(-alpha * thi * thi);
    const int mmax = l + CPG_SEP_A - 1;
    double Sm[CPG_LMAX + CPG_SEP_A], plo = 1.0, phi = 1.0;
    for (int m = 0; m <= mmax; ++m) {
        const double t2 = m / (2.0 * alpha);
        if (t2 <= tlo * tlo)      Sm[m] = plo * elo;
        else if (t2 >= thi * thi) Sm[m] = phi * ehi;
        else                      Sm[m] = pow(t2, 0.5 * m) * exp_neg_half[m];
        plo *= tlo;
        phi *= thi;
    }
    double ap[CPG_LMAX + CPG_SEP_A];       // alpha^n
    ap[0] = 1.0;
    for (int n = 1; n < CPG_LMAX + CPG_SEP_A; ++n) ap[n] = ap[n - 1] * alpha;
    for (int i = 0; i <= l; ++i) {
        for (int a = 3; a < CPG_SEP_A; ++a) {
            double acc = 0.0;
            for (int m = (i + a) & 1; m <= i + a; m += 2) {      // only the parity of i + a occurs
                const double k = ctx->kap_sep[i][a][m];
                if (k != 0.0) acc += fabs(k) * ap[(m - i + a) / 2] * Sm[m];
            }
            sup[i][a] = acc;
        }
    }
}

// Per-AO remainder bounds over the cube c +- h, in the layout cpg_box_eval uses: E[0] E1, 1 E2, 2 E3,
// 3..5 Ea2, 6..8 Ea3, 9..14 Eab2 (a <= b). E_K of f = d^beta phi bounds |f(c + d) - T_{K-1} f (d)|,
// T_{K-1} the Taylor polynomial of degree K - 1 at the centre.
//
// Each primitive is a product of 1D factors X_k = d^beta_k/dx^beta_k [x^i exp(-alpha x^2)]. Each factor
// is written as its Taylor polynomial at the centre of degree P = K - 1 + CPG_SEP_EXTRA (exact
// coefficients) plus a 1D Lagrange remainder (sup on the interval of the next derivative). Multiplied
// out, the terms of total degree < K are exactly T_{K-1} f, so |f - T_{K-1} f| is at most the sum of
// the absolute values of all the other terms (a remainder counts with degree P + 1). Only the 1D
// remainders use sups; the exact low-order cross terms carry the variation over the cube.
//
// Against a multivariate Lagrange bound (sups of every order-K multi-index derivative over the cube,
// each at its own worst point) this is 3-4x tighter per AO (geometric mean bound / sampled remainder
// 1.6-2.3 instead of 5-8), and a sweep needs 50-60% fewer cubes. Checked against the remainder
// sampled on a 5^3 lattice per cube for 68M (AO, cube, slot) triples (water/cc-pVDZ, acro-xps):
// never below it beyond the sampling's own rounding.
//
// The sum over degrees >= K is formed from suffix sums, without subtracting anything (so it also
// carries over to float with a plain relative error bound): with u_k[m] the absolute degree-m term of
// axis k and S_k[r] = sum over m >= r of u_k[m],
//   sum over m0 + m1 + m2 >= K = sum over m0 of u_0[m0] T12[K - m0],  T12[j] = sum over m1 + m2 >= j,
//   T12[j] = S_1[j] S_2[0] + sum over m1 < j of u_1[m1] S_2[j - m1]   (T12[j <= 0] = S_1[0] S_2[0]).
static void cpg_shell_rem(const cpg_ctx_t* ctx, const cpg_shell_t* s, const double c[3], double h, double* E, int stride) {
    enum { NBP = 6 };   // (beta, P) per axis: (0,2) (0,3) (1,3) (2,3) (0,4) (1,4)
    static const int bp_beta[NBP] = { 0, 0, 1, 2, 0, 1 }, bp_P[NBP] = { 2, 3, 3, 3, 4, 4 };
    // per slot: K and the (beta, P) index per axis
    static const int slot_K[15] = { 1, 2, 3, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2 };
    static const int slot_bp[15][3] = {
        {0,0,0}, {1,1,1}, {4,4,4},          // E1 E2 E3
        {2,1,1}, {1,2,1}, {1,1,2},          // Ea2
        {5,4,4}, {4,5,4}, {4,4,5},          // Ea3
        {3,1,1}, {2,2,1}, {2,1,2}, {1,3,1}, {1,2,2}, {1,1,3},   // Eab2: xx xy xz yy yz zz
    };
    static const double inv_fact[CPG_SEP_A] = { 1.0, 1.0, 1.0 / 2, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720 };
    STATIC_ASSERT(CPG_SEP_EXTRA == 2 && CPG_SEP_A == 7, "the (beta, P) table assumes two extra terms");
    for (int ci = 0; ci < s->ncart; ++ci) for (int t = 0; t < 15; ++t) E[ci * stride + t] = 0.0;
    double val[CPG_LMAX + 1][CPG_SEP_A], sup[CPG_LMAX + 1][CPG_SEP_A], hp[CPG_SEP_A];
    double u[3][CPG_LMAX + 1][NBP][CPG_SEP_A], S[3][CPG_LMAX + 1][NBP][CPG_SEP_A + 1];
    hp[0] = 1.0;
    for (int m = 1; m < CPG_SEP_A; ++m) hp[m] = hp[m - 1] * h;
    for (uint32_t ip = 0; ip < s->num_prims; ++ip) {
        const double al = ctx->alpha[s->prim_offset + ip];
        const double cf = fabs(ctx->coeff[s->prim_offset + ip]);
        for (int k = 0; k < 3; ++k) {
            cpg_sep_1d(ctx, s->l, al, c[k] - s->A[k], h, val, sup);
            for (int i = 0; i <= s->l; ++i) {
                for (int b = 0; b < NBP; ++b) {
                    const int be = bp_beta[b], P = bp_P[b];
                    double* uu = u[k][i][b];
                    double* SS = S[k][i][b];
                    for (int m = 0; m <= P; ++m) uu[m] = fabs(val[i][be + m]) * hp[m] * inv_fact[m];
                    uu[P + 1] = sup[i][be + P + 1] * hp[P + 1] * inv_fact[P + 1];
                    SS[P + 2] = 0.0;
                    for (int m = P + 1; m >= 0; --m) SS[m] = SS[m + 1] + uu[m];
                }
            }
        }
        for (int ci = 0; ci < s->ncart; ++ci) {
            const int* ijk = ctx->ao_ijk[s->ao_offset + ci];
            const double w = cf * ctx->ao_nrm[s->ao_offset + ci];
            for (int t = 0; t < 15; ++t) {
                const int K = slot_K[t];
                const double* u0 = u[0][ijk[0]][slot_bp[t][0]];
                const double* u1 = u[1][ijk[1]][slot_bp[t][1]];
                const double* S0 = S[0][ijk[0]][slot_bp[t][0]];
                const double* S1 = S[1][ijk[1]][slot_bp[t][1]];
                const double* S2 = S[2][ijk[2]][slot_bp[t][2]];
                double T12[4];
                T12[0] = S1[0] * S2[0];
                for (int j = 1; j <= K; ++j) {
                    double v = S1[j] * S2[0];
                    for (int m1 = 0; m1 < j; ++m1) v += u1[m1] * S2[j - m1];
                    T12[j] = v;
                }
                double bound = S0[K] * T12[0];
                for (int m0 = 0; m0 < K; ++m0) bound += u0[m0] * T12[K - m0];
                E[ci * stride + t] += w * bound;
            }
        }
    }
}

static inline double cpg_box_dist2(const double A[3], const double c[3], double h) {
    double d2 = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double t = fabs(A[k] - c[k]) - h;
        if (t > 0.0) d2 += t * t;
    }
    return d2;
}

// A shell at squared distance d2 from a region: -1 local (within the local radius), 0 far (beyond the far
// radius: in the global tails), m >= 1: left out, its AOs and their derivatives (order <= 2) below
// rung_tau[m] over the region. Conservative: a shell is called near or far only if it surely is.
static inline int cpg_screen_rung(const cpg_ctx_t* ctx, const cpg_shell_t* s, double d2) {
    const double* r2 = ctx->rung_r2 + (size_t)s->cls * CPG_RUNGS;
    const double dl = d2 * (1.0 - 1e-12);
    const int loc = ctx->rung_loc;
    if (dl < r2[loc] * (1.0 + 1e-6)) return -1;                     // most shells of a list: one test
    for (int m = 0; m < loc; ++m) if (dl >= r2[m] * (1.0 + 1e-6)) return m;   // the smallest m with r_m <= d
    return loc;
}

// The cube's local AOs (shells within the local radius) and its tails: the global ones plus those of the
// shells between the local and the far radius. Returns the AO count.
static int cpg_local_aos_region(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double c[3], double h, int* out_num_shells, double tail[3]) {
    int n = 0, ns = 0;
    double t[3] = { 0.0, 0.0, 0.0 };
    for (int si = 0; si < ctx->nshell; ++si) {
        const cpg_shell_t* s = &ctx->shell[si];
        const int m = cpg_screen_rung(ctx, s, cpg_box_dist2(s->A, c, h));
        if (m < 0) {
            sc->shells[ns++] = si;
            for (int ci = 0; ci < s->ncart; ++ci) sc->L[n++] = (int)s->ao_offset + ci;
        } else if (m > 0) {
            for (int k = 0; k < 3; ++k) t[k] += ctx->rung_tau[m] * s->tw[k];
        }
    }
    tail[0] = ctx->tail_rho + t[0];
    tail[1] = ctx->tail_g + t[1];
    tail[2] = ctx->tail_H + t[2];
    *out_num_shells = ns;
    return n;
}

// Local AO list of a point or a box around a path: the shells whose far screening radius reaches it (the
// polish: CP positions as if nothing were left out), or for the separatrices their local radius (a path's
// direction needs no more than a cube's tests do). Returns count.
static int cpg_local_aos(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double c[3], double h, int* out_num_shells, bool trace) {
    int n = 0, ns = 0;
    for (int si = 0; si < ctx->nshell; ++si) {
        const cpg_shell_t* s = &ctx->shell[si];
        const double r = trace ? s->radius_trace : s->radius;
        if (cpg_box_dist2(s->A, c, h) < r * r) {
            sc->shells[ns++] = si;
            for (int ci = 0; ci < s->ncart; ++ci) sc->L[n++] = (int)s->ao_offset + ci;
        }
    }
    *out_num_shells = ns;
    return n;
}

// ----------------------------------------------------------------------------------------------- point field

// rho, grad, Hessian at x (deriv 0..2) over all shells within their screening radius (far, or for the
// separatrices local: cpg_local_aos).
static void cpg_point_r(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], int deriv, double* rho, double g[3], double H[3][3], bool trace) {
    int ns = 0;
    const int n = cpg_local_aos(ctx, sc, x, 0.0, &ns, trace);
    sc->trace_aos += (uint64_t)n;   // separatrix statistics (read as differences around each path)
    int k = 0;
    for (int t = 0; t < ns; ++t) {
        const cpg_shell_t* s = &ctx->shell[sc->shells[t]];
        cpg_shell_eval(ctx, s, x, deriv, sc->V + (size_t)k * CPG_NV, CPG_NV);
        k += s->ncart;
    }
    // D V for the channels needed (1 for the gradient, 4 for the Hessian). With one channel a row's sum
    // is a single chain of dependent adds; CPG_PT_ROWS rows at a time keep that many chains in flight
    // (separatrices: ~5x in the matrix form). Each row still sums over j in order: the same result.
    const int ndv = deriv >= 2 ? 4 : 1;
    #define CPG_PT_ROWS 8
    for (int i0 = 0; i0 < n; i0 += CPG_PT_ROWS) {
        const int nr = MIN(CPG_PT_ROWS, n - i0);
        double acc[CPG_PT_ROWS][4] = {{0}};
        const double* row[CPG_PT_ROWS];
        for (int t = 0; t < CPG_PT_ROWS; ++t) row[t] = ctx->D + (size_t)sc->L[i0 + MIN(t, nr - 1)] * ctx->nao;
        if (ndv == 1) {
            for (int j = 0; j < n; ++j) {
                const int lj = sc->L[j];
                const double v0 = sc->V[(size_t)j * CPG_NV];
                for (int t = 0; t < CPG_PT_ROWS; ++t) acc[t][0] += row[t][lj] * v0;
            }
        } else {
            for (int j = 0; j < n; ++j) {
                const int lj = sc->L[j];
                const double* v = sc->V + (size_t)j * CPG_NV;
                for (int t = 0; t < CPG_PT_ROWS; ++t) {
                    const double d = row[t][lj];
                    for (int m = 0; m < 4; ++m) acc[t][m] += d * v[m];
                }
            }
        }
        for (int t = 0; t < nr; ++t) for (int m = 0; m < ndv; ++m) sc->Dv[(size_t)(i0 + t) * CPG_NDV + m] = acc[t][m];
    }
    #undef CPG_PT_ROWS
    double r = 0.0, gg[3] = {0, 0, 0}, HH[3][3] = {{0}};
    for (int i = 0; i < n; ++i) {
        const double* v = sc->V + (size_t)i * CPG_NV;
        const double* dv = sc->Dv + (size_t)i * CPG_NDV;
        r += v[0] * dv[0];
        if (deriv >= 1) for (int a = 0; a < 3; ++a) gg[a] += 2.0 * v[1 + a] * dv[0];
        if (deriv >= 2) {
            for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
                HH[a][b] += 2.0 * (v[M2(a, b)] * dv[0] + v[1 + a] * dv[1 + b]);
            }
        }
    }
    if (rho) *rho = r;
    if (g) for (int a = 0; a < 3; ++a) g[a] = gg[a];
    if (H) for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) H[a][b] = a <= b ? HH[a][b] : HH[b][a];
}

static void cpg_point(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], int deriv, double* rho, double g[3], double H[3][3]) {
    cpg_point_r(ctx, sc, x, deriv, rho, g, H, false);
}

// ----------------------------------------------------------------------------------------------- cube enclosures

static void cpg_box_eval(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double c[3], double h, cpg_eval_t* ev) {
    int ns = 0;
    double tail[3];
    const int n = cpg_local_aos_region(ctx, sc, c, h, &ns, tail);
    memset(ev, 0, sizeof(*ev));
    // centre values (orders 0..3) and per-AO remainder bounds over the cube
    int k = 0;
    for (int t = 0; t < ns; ++t) {
        const cpg_shell_t* s = &ctx->shell[sc->shells[t]];
        cpg_shell_eval(ctx, s, c, 3, sc->V + (size_t)k * CPG_NV, CPG_NV);
        cpg_shell_rem(ctx, s, c, h, sc->E + (size_t)k * 16, 16);
        k += s->ncart;
    }
    // The rows of rho = sum_ij D_ij phi_i phi_j, per local AO: V (Vr), the remainder vector (Er), Dv = D V
    // and Q = |D| E[0..5].
    const int nr = n;
    const double* Vr = sc->V;
    const double* Er = sc->E;
    {
        // D-products Dv = D V (all 20 channels) and the |D| products Q[i*8 + .]: 0 |D|E1, 1 |D|E2,
        // 2 |D|E3, 3..5 |D|Ea2, in one pass over the local rows of D (read in place: no n x n copy,
        // which was the per-thread scratch's N^2 for a large basis)
        for (int i = 0; i < n; ++i) {
            const double* row = ctx->D + (size_t)sc->L[i] * ctx->nao;
            double acc[CPG_NDV] = {0};
            double q[6] = {0};
            for (int j = 0; j < n; ++j) {
                const double d = row[sc->L[j]];
                const double* v = sc->V + (size_t)j * CPG_NV;
                for (int m = 0; m < CPG_NDV; ++m) acc[m] += d * v[m];
                const double ad = fabs(d);
                const double* e = sc->E + (size_t)j * 16;
                q[0] += ad * e[0]; q[1] += ad * e[1]; q[2] += ad * e[2];
                q[3] += ad * e[3]; q[4] += ad * e[4]; q[5] += ad * e[5];
            }
            memcpy(sc->Dv + (size_t)i * CPG_NDV, acc, sizeof(acc));
            memcpy(sc->Q + (size_t)i * 8, q, sizeof(q));
        }
    }
    // DOT[a][b] = V_a . X V_b for a < 20, b < 10
    double DOT[CPG_NV][10];
    memset(DOT, 0, sizeof(DOT));
    for (int i = 0; i < nr; ++i) {
        const double* v = Vr + (size_t)i * CPG_NV;
        const double* dv = sc->Dv + (size_t)i * CPG_NDV;
        for (int a = 0; a < CPG_NV; ++a) {
            const double va = v[a];
            for (int b = 0; b < 10; ++b) DOT[a][b] += va * dv[b];
        }
    }
    ev->rho = DOT[0][0];
    for (int a = 0; a < 3; ++a) ev->g[a] = 2.0 * DOT[M1(a)][0];
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) ev->A[a][b] = 2.0 * (DOT[M2(a, b)][0] + DOT[M1(a)][M1(b)]);
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) for (int cc = 0; cc < 3; ++cc) {
        ev->T[a][b][cc] = 2.0 * (DOT[M3(a, b, cc)][0] + DOT[M2(a, b)][M1(cc)] + DOT[M2(a, cc)][M1(b)] + DOT[M2(b, cc)][M1(a)]);
    }
    // per-AO remainder vectors (cpg_shell_rem above). E[i*16 + .]: 0 E1, 1 E2, 2 E3, 3..5 Ea2, 6..8 Ea3, 9..14 Eab2 (a<=b)
    const double h2 = h * h, h3 = h2 * h;
    int ab_idx[3][3];
    { int t = 0; for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) { ab_idx[a][b] = ab_idx[b][a] = t++; } }
    // row sums of the remainder terms
    double rho_lin = 0, rho_quad = 0;
    double err_g[3] = {0}, err_H[6] = {0};
    for (int i = 0; i < nr; ++i) {
        const double* dv = sc->Dv + (size_t)i * CPG_NDV;
        const double* e = Er + (size_t)i * 16;
        const double* q = sc->Q + (size_t)i * 8;
        double ad[CPG_NDV];
        for (int m = 0; m < CPG_NDV; ++m) ad[m] = fabs(dv[m]);
        double sum1 = 0, sum2 = 0;
        for (int p = 0; p < 3; ++p) { sum1 += ad[M1(p)]; for (int r = 0; r < 3; ++r) sum2 += ad[M2(p, r)]; }
        const double W  = ad[0] + h * sum1 + 0.5 * h2 * sum2;
        const double W1 = ad[0] + h * sum1;
        double Wa[3], W1a[3], W1ab[6];
        for (int a = 0; a < 3; ++a) {
            double s1 = 0, s2 = 0;
            for (int p = 0; p < 3; ++p) { s1 += ad[M2(a, p)]; for (int r = 0; r < 3; ++r) s2 += ad[M3(a, p, r)]; }
            Wa[a]  = ad[M1(a)] + h * s1 + 0.5 * h2 * s2;
            W1a[a] = ad[M1(a)] + h * s1;
        }
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
            double s1 = 0;
            for (int p = 0; p < 3; ++p) s1 += ad[M3(a, b, p)];
            W1ab[ab_idx[a][b]] = ad[M2(a, b)] + h * s1;
        }
        rho_lin  += ad[0] * e[0];
        rho_quad += e[0] * q[0];
        for (int a = 0; a < 3; ++a) err_g[a] += e[6 + a] * W + e[2] * Wa[a] + e[6 + a] * q[2];
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
            const int t = ab_idx[a][b];
            err_H[t] += e[9 + t] * W1 + W1ab[t] * e[1] + e[9 + t] * q[1] + e[3 + a] * W1a[b] + W1a[a] * e[3 + b] + e[3 + a] * q[3 + b];
        }
    }
    // polynomial terms with exact coefficients
    for (int a = 0; a < 3; ++a) {
        double deg3 = 0, deg4 = 0, quadT = 0;
        for (int p = 0; p < 3; ++p) for (int q = 0; q < 3; ++q) {
            quadT += fabs(ev->T[a][p][q]);
            for (int r = 0; r < 3; ++r) {
                deg3 += fabs(DOT[M3(a, p, q)][M1(r)] + DOT[M2(a, p)][M2(q, r)]);
                for (int s = 0; s < 3; ++s) deg4 += fabs(DOT[M3(a, p, q)][M2(r, s)]);
            }
        }
        ev->r[a] = h3 * deg3 + 0.5 * h2 * h2 * deg4 + 2.0 * err_g[a] + tail[1];
        ev->quadT[a] = 0.5 * h2 * quadT;
    }
    for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
        double lin = 0, quad = 0;
        for (int p = 0; p < 3; ++p) {
            lin += fabs(ev->T[a][b][p]);
            for (int q = 0; q < 3; ++q) quad += fabs(DOT[M3(a, b, p)][M1(q)] + DOT[M2(a, p)][M2(b, q)]);
        }
        const double v = h * lin + 2.0 * h2 * quad + 2.0 * err_H[ab_idx[a][b]] + tail[2];
        ev->dH[a][b] = ev->dH[b][a] = v;
    }
    ev->rho_up = ev->rho + 2.0 * rho_lin + rho_quad + tail[0];
}

// ----------------------------------------------------------------------------------------------- tests

// Range of q(d) = c0 + b.d + 1/2 d^T M d over the box |d_k| <= h, exactly: the extremes of a quadratic
// over a box lie at stationary points of its restriction to the relative interior of some face (the
// box itself, a face, an edge or a corner), so all 27 are tried.
static void cpg_quad_box_range(double c0, const double b[3], double M[3][3], double h, double* out_min, double* out_max) {
    double lo = DBL_MAX, hi = -DBL_MAX;
    for (int pat = 0; pat < 27; ++pat) {
        int st[3] = { pat % 3, (pat / 3) % 3, pat / 9 };   // 0 free, 1 -h, 2 +h
        double d[3] = { 0, 0, 0 };
        int F[3], nf = 0;
        for (int k = 0; k < 3; ++k) {
            if (st[k] == 0) F[nf++] = k;
            else d[k] = st[k] == 1 ? -h : h;
        }
        if (nf > 0) {
            // M_FF x = -(b_F + M_F,fixed d_fixed)
            double r[3], m[3][3];
            for (int i = 0; i < nf; ++i) {
                r[i] = -b[F[i]];
                for (int k = 0; k < 3; ++k) if (st[k] != 0) r[i] -= M[F[i]][k] * d[k];
                for (int j = 0; j < nf; ++j) m[i][j] = M[F[i]][F[j]];
            }
            double x[3];
            double scale = 0.0;
            for (int i = 0; i < nf; ++i) for (int j = 0; j < nf; ++j) scale = fmax(scale, fabs(m[i][j]));
            if (scale == 0.0) continue;
            if (nf == 1) {
                if (fabs(m[0][0]) < 1e-14 * scale) continue;
                x[0] = r[0] / m[0][0];
            } else if (nf == 2) {
                const double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
                if (fabs(det) < 1e-14 * scale * scale) continue;
                x[0] = (r[0] * m[1][1] - m[0][1] * r[1]) / det;
                x[1] = (m[0][0] * r[1] - r[0] * m[1][0]) / det;
            } else {
                const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                if (fabs(det) < 1e-14 * scale * scale * scale) continue;
                x[0] = (r[0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (r[1] * m[2][2] - m[1][2] * r[2]) + m[0][2] * (r[1] * m[2][1] - m[1][1] * r[2])) / det;
                x[1] = (m[0][0] * (r[1] * m[2][2] - m[1][2] * r[2]) - r[0] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * r[2] - r[1] * m[2][0])) / det;
                x[2] = (m[0][0] * (m[1][1] * r[2] - r[1] * m[2][1]) - m[0][1] * (m[1][0] * r[2] - r[1] * m[2][0]) + r[0] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])) / det;
            }
            bool inside = true;
            for (int i = 0; i < nf; ++i) { if (!(fabs(x[i]) <= h)) inside = false; d[F[i]] = x[i]; }
            if (!inside) continue;
        }
        double q = c0;
        for (int i = 0; i < 3; ++i) {
            q += b[i] * d[i];
            for (int j = 0; j < 3; ++j) q += 0.5 * M[i][j] * d[i] * d[j];
        }
        lo = fmin(lo, q);
        hi = fmax(hi, q);
    }
    *out_min = lo;
    *out_max = hi;
}

static bool cpg_exclude(const cpg_eval_t* ev, double h, const double val[3], double vec[3][3]) {
    double res[3];
    for (int k = 0; k < 3; ++k) res[k] = ev->r[k] + ev->quadT[k];
    double dirs[7][3] = {{1,0,0},{0,1,0},{0,0,1}};
    int nd = 3;
    const double gn = sqrt(ev->g[0] * ev->g[0] + ev->g[1] * ev->g[1] + ev->g[2] * ev->g[2]);
    if (gn > 0.0) { for (int k = 0; k < 3; ++k) dirs[nd][k] = ev->g[k] / gn; nd++; }
    for (int e = 0; e < 3; ++e) { for (int k = 0; k < 3; ++k) dirs[nd][k] = vec[k][e]; nd++; }
    for (int t = 0; t < nd; ++t) {
        const double* w = dirs[t];
        const double lhs = fabs(w[0] * ev->g[0] + w[1] * ev->g[1] + w[2] * ev->g[2]);
        double rhs = 0.0;
        for (int j = 0; j < 3; ++j) {
            const double Aw = ev->A[j][0] * w[0] + ev->A[j][1] * w[1] + ev->A[j][2] * w[2];
            rhs += fabs(Aw) * h + fabs(w[j]) * res[j];
        }
        if (lhs > rhs) return true;
        // the same direction with the quadratic part of the model minimised exactly over the box rather
        // than bounded term by term (which ignores where on the box each term peaks)
        double b[3], M[3][3], qlo, qhi, Rw = 0.0, scale = lhs;
        for (int j = 0; j < 3; ++j) {
            b[j] = ev->A[j][0] * w[0] + ev->A[j][1] * w[1] + ev->A[j][2] * w[2];
            for (int k = 0; k < 3; ++k) M[j][k] = w[0] * ev->T[0][j][k] + w[1] * ev->T[1][j][k] + w[2] * ev->T[2][j][k];
            Rw += fabs(w[j]) * ev->r[j];
            scale += fabs(b[j]) * h;
            for (int k = 0; k < 3; ++k) scale += 0.5 * fabs(M[j][k]) * h * h;
        }
        cpg_quad_box_range(w[0] * ev->g[0] + w[1] * ev->g[1] + w[2] * ev->g[2], b, M, h, &qlo, &qhi);
        const double margin = Rw + 1e-12 * scale;     // the extremes are found by small solves in double
        if (qlo > margin || qhi < -margin) return true;
    }
    // Newton point far outside the cube
    const double smin = fmin(fabs(val[0]), fmin(fabs(val[1]), fabs(val[2])));
    if (smin > 0.0) {
        double dstar[3] = {0, 0, 0};
        for (int e = 0; e < 3; ++e) {
            const double proj = vec[0][e] * ev->g[0] + vec[1][e] * ev->g[1] + vec[2][e] * ev->g[2];
            for (int k = 0; k < 3; ++k) dstar[k] -= vec[k][e] * proj / val[e];
        }
        double dist2 = 0.0;
        for (int k = 0; k < 3; ++k) { const double t = fabs(dstar[k]) - h; if (t > 0) dist2 += t * t; }
        const double rn = sqrt(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]);
        if (smin * sqrt(dist2) > rn) return true;
    }
    return false;
}

// Krawczyk on X = c +- h with Y = A^-1: returns +1 (unique zero), -1 (no zero), 0 (undecided).
static int cpg_krawczyk(const double g[3], double Y[3][3], double dH[3][3], double h) {
    double m[3], rad[3];
    for (int i = 0; i < 3; ++i) {
        m[i] = -(Y[i][0] * g[0] + Y[i][1] * g[1] + Y[i][2] * g[2]);
        double s = 0.0;
        for (int j = 0; j < 3; ++j) s += fabs(Y[i][j]) * (dH[j][0] + dH[j][1] + dH[j][2]) * h;
        rad[i] = s;
    }
    bool inside = true, disjoint = false;
    for (int i = 0; i < 3; ++i) {
        if (!(fabs(m[i]) + rad[i] < h)) inside = false;
        if (fabs(m[i]) - rad[i] > h) disjoint = true;
    }
    return inside ? 1 : (disjoint ? -1 : 0);
}

// ----------------------------------------------------------------------------------------------- Newton polish

static bool cpg_newton(const cpg_ctx_t* ctx, cpg_scratch_t* sc, double x[3], double* out_rho, double out_val[3], double out_vec[3][3]) {
    double rho, g[3], H[3][3], val[3], vec[3][3], Y[3][3];
    for (int it = 0; it < 60; ++it) {
        cpg_point(ctx, sc, x, 2, &rho, g, H);
        cpg_eigen_sym3(val, vec, H);
        if (!cpg_inverse_sym3(Y, val, vec)) return false;
        double dx[3], n2 = 0.0;
        for (int i = 0; i < 3; ++i) { dx[i] = -(Y[i][0] * g[0] + Y[i][1] * g[1] + Y[i][2] * g[2]); n2 += dx[i] * dx[i]; }
        for (int i = 0; i < 3; ++i) x[i] += dx[i];
        if (sqrt(n2) < 1e-13) break;
        if (it == 59) return false;
    }
    cpg_point(ctx, sc, x, 2, &rho, g, H);
    cpg_eigen_sym3(out_val, out_vec, H);
    *out_rho = rho;
    return true;
}

static int cpg_classify(const double val[3]) {
    int npos = 0;
    for (int i = 0; i < 3; ++i) npos += val[i] > 0.0;
    switch (npos) {
        case 0:  return MD_TOPO_MAXIMUM;
        case 1:  return MD_TOPO_SPLIT_SADDLE;   // (3,-1) bond critical point
        case 2:  return MD_TOPO_JOIN_SADDLE;    // (3,+1) ring critical point
        default: return MD_TOPO_MINIMUM;        // (3,+3) cage critical point
    }
}

// ----------------------------------------------------------------------------------------------- separatrices

// rho and grad rho for the separatrix tracer. The matrix form reads n^2 elements of D per point, and the AOs
// near a point are scattered over D's columns (shells are ordered by l before atom): on a 9636-AO protein
// every evaluation was a gather of ~640k doubles from DRAM (3.7 ms). Here the local AOs of a box
// x +- CPG_TRACE_BOX are taken once with their block of D, contiguous, and used for every point of the path
// inside the box. The box's list holds each point's own list; the AOs of the shells that do not reach the
// point are 0 and skipped. The same sums as cpg_point's up to the order of the terms and D's symmetry (the
// upper triangle is read). A box whose block would exceed CPG_TRACE_BLOCK doubles sends its points to
// cpg_point.
#define CPG_TRACE_BOX   0.5
#define CPG_TRACE_BLOCK (1536u * 1536u)   // 18 MB per worker at most
static void cpg_point_trace(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], double* rho, double g[3], md_allocator_i* heap) {
    bool inside = sc->tc_valid;
    for (int k = 0; k < 3 && inside; ++k) inside = fabs(x[k] - sc->tc_c[k]) <= CPG_TRACE_BOX;
    if (inside && sc->tc_direct) {
        cpg_point_r(ctx, sc, x, 1, rho, g, NULL, true);
        return;
    }
    if (!inside) {
        int ns = 0;
        const int n = cpg_local_aos(ctx, sc, x, CPG_TRACE_BOX, &ns, true);
        MEMCPY(sc->tc_c, x, sizeof(double) * 3);
        sc->tc_valid = true;
        sc->tc_direct = (size_t)n * n > CPG_TRACE_BLOCK;
        if (sc->tc_direct) {
            cpg_point_r(ctx, sc, x, 1, rho, g, NULL, true);
            return;
        }
        if ((size_t)n * n > sc->tc_cap) {
            if (sc->tc_B) md_free(heap, sc->tc_B, sizeof(double) * sc->tc_cap);
            sc->tc_cap = MAX((size_t)n * n, (size_t)256 * 256);
            sc->tc_B = (double*)md_alloc(heap, sizeof(double) * sc->tc_cap);
        }
        MEMCPY(sc->tc_L, sc->L, sizeof(int) * n);
        MEMCPY(sc->tc_shells, sc->shells, sizeof(int) * ns);
        for (int i = 0; i < n; ++i) {     // the upper triangle (read as symmetric below)
            const double* row = ctx->D + (size_t)sc->tc_L[i] * ctx->nao;
            double* b = sc->tc_B + (size_t)i * n;
            for (int j = i; j < n; ++j) b[j] = row[sc->tc_L[j]];
        }
        sc->tc_n = n;
        sc->tc_ns = ns;
    }
    const int n = sc->tc_n;
    // the AOs of the shells that reach x (cpg_local_aos's test) are the active rows; the others are 0
    int* act = sc->L;
    int k = 0, na = 0;
    for (int t = 0; t < sc->tc_ns; ++t) {
        const cpg_shell_t* s = &ctx->shell[sc->tc_shells[t]];
        double* v = sc->V + (size_t)k * CPG_NV;
        if (cpg_box_dist2(s->A, x, 0.0) < s->radius_trace * s->radius_trace) {
            cpg_shell_eval(ctx, s, x, 1, v, CPG_NV);
            for (int c = 0; c < s->ncart; ++c) act[na++] = k + c;
        } else {
            for (int c = 0; c < s->ncart; ++c) for (int m = 0; m < 4; ++m) v[(size_t)c * CPG_NV + m] = 0.0;
        }
        k += s->ncart;
    }
    sc->trace_aos += (uint64_t)na;
    sc->trace_block_aos += (uint64_t)n;
    // y = B phi for the active rows, B read as its upper triangle (D is symmetric): row i adds its j >= i part
    // to y_i and B_ij phi_i to y_j. Memory-bound at this size, so half the bytes is half the time; each
    // active row once, the inactive ones (phi_i = 0) not at all. Then rho = phi . y, grad = 2 dphi . y.
    double* vv = sc->E;      // phi, contiguous
    double* y  = sc->Dv;
    for (int j = 0; j < n; ++j) { vv[j] = sc->V[(size_t)j * CPG_NV]; y[j] = 0.0; }
    for (int q = 0; q < na; ++q) {
        const int i = act[q];
        const double* b = sc->tc_B + (size_t)i * n;
        const double vi = vv[i];
        double sl[8] = { b[i] * vi, 0, 0, 0, 0, 0, 0, 0 };
        int j = i + 1;
        for (; j + 8 <= n; j += 8) {
            for (int u = 0; u < 8; ++u) { const double bj = b[j + u]; sl[u] += bj * vv[j + u]; y[j + u] += bj * vi; }
        }
        for (; j < n; ++j) { sl[0] += b[j] * vv[j]; y[j] += b[j] * vi; }
        y[i] += ((sl[0] + sl[1]) + (sl[2] + sl[3])) + ((sl[4] + sl[5]) + (sl[6] + sl[7]));
    }
    double r = 0.0, gg[3] = { 0, 0, 0 };
    for (int q = 0; q < na; ++q) {
        const int i = act[q];
        const double* v = sc->V + (size_t)i * CPG_NV;
        r += v[0] * y[i];
        for (int a = 0; a < 3; ++a) gg[a] += 2.0 * v[1 + a] * y[i];
    }
    if (rho) *rho = r;
    for (int a = 0; a < 3; ++a) g[a] = gg[a];
}

static void cpg_unit_grad(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], double sign, double out[3], double* rho, double* gnorm) {
    double g[3];
    cpg_point_trace(ctx, sc, x, rho, g, md_get_heap_allocator());
    sc->trace_evals++;
    const double n = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    *gnorm = n;
    for (int k = 0; k < 3; ++k) out[k] = n > 0 ? sign * g[k] / n : 0.0;
}

// Follows sign * grad rho / |grad rho| from x0 + 1e-3 dir. Returns index of the CP of type 'target'
// it ends in, or -1 (density fell below eps, or path length exhausted). Only the end point is kept,
// and it is found within 'hit' of a CP that attracts the path, so the local error per step (1e-6)
// and the step cap (0.25 Bohr) are loose: on every test input the end points are those of 1e-8 and
// 0.05, at 40% fewer evaluations and a third of the longest path's time.
// rho_stop: the path is abandoned once rho passes it (ascending: above it; descending: below it). rho is
// monotone along the flow, so with rho_stop the highest maximum's density (ascending) or the lowest cage
// point's (descending), with a margin, no target can be reached beyond it: on a protein, 439 of 4450 paths
// reached no CP, the longest after 11947 evaluations.
static int cpg_trace(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x0[3], const double dir[3], double sign, const cpg_cp_t* cps, int ncp, int target,
                     double rho_stop) {
    static const double a21 = 1.0/5;
    static const double a31 = 3.0/40, a32 = 9.0/40;
    static const double a41 = 44.0/45, a42 = -56.0/15, a43 = 32.0/9;
    static const double a51 = 19372.0/6561, a52 = -25360.0/2187, a53 = 64448.0/6561, a54 = -212.0/729;
    static const double a61 = 9017.0/3168, a62 = -355.0/33, a63 = 46732.0/5247, a64 = 49.0/176, a65 = -5103.0/18656;
    static const double b1 = 35.0/384, b3 = 500.0/1113, b4 = 125.0/192, b5 = -2187.0/6784, b6 = 11.0/84;
    static const double e1 = 71.0/57600, e3 = -71.0/16695, e4 = 71.0/1920, e5 = -17253.0/339200, e6 = 22.0/525, e7 = -1.0/40;
    double x[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], y[3], rho, gn;
    for (int k = 0; k < 3; ++k) x[k] = x0[k] + 1e-3 * dir[k];
    double step = 1e-3, len = 0.0;
    const double tol = 1e-6, hit = 2e-3;
    cpg_unit_grad(ctx, sc, x, sign, k1, &rho, &gn);
    for (int iter = 0; iter < 200000 && len < 40.0; ++iter) {
        for (int k = 0; k < 3; ++k) y[k] = x[k] + step * (a21 * k1[k]);
        cpg_unit_grad(ctx, sc, y, sign, k2, &rho, &gn);
        for (int k = 0; k < 3; ++k) y[k] = x[k] + step * (a31 * k1[k] + a32 * k2[k]);
        cpg_unit_grad(ctx, sc, y, sign, k3, &rho, &gn);
        for (int k = 0; k < 3; ++k) y[k] = x[k] + step * (a41 * k1[k] + a42 * k2[k] + a43 * k3[k]);
        cpg_unit_grad(ctx, sc, y, sign, k4, &rho, &gn);
        for (int k = 0; k < 3; ++k) y[k] = x[k] + step * (a51 * k1[k] + a52 * k2[k] + a53 * k3[k] + a54 * k4[k]);
        cpg_unit_grad(ctx, sc, y, sign, k5, &rho, &gn);
        for (int k = 0; k < 3; ++k) y[k] = x[k] + step * (a61 * k1[k] + a62 * k2[k] + a63 * k3[k] + a64 * k4[k] + a65 * k5[k]);
        cpg_unit_grad(ctx, sc, y, sign, k6, &rho, &gn);
        double xn[3];
        for (int k = 0; k < 3; ++k) xn[k] = x[k] + step * (b1 * k1[k] + b3 * k3[k] + b4 * k4[k] + b5 * k5[k] + b6 * k6[k]);
        cpg_unit_grad(ctx, sc, xn, sign, k7, &rho, &gn);
        double err = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double e = step * (e1 * k1[k] + e3 * k3[k] + e4 * k4[k] + e5 * k5[k] + e6 * k6[k] + e7 * k7[k]);
            err = fmax(err, fabs(e));
        }
        if (err <= tol || step <= 1e-7) {
            for (int k = 0; k < 3; ++k) { x[k] = xn[k]; k1[k] = k7[k]; }
            len += step;
            for (int i = 0; i < ncp; ++i) {
                if (cps[i].type != target) continue;
                const double dx = x[0] - cps[i].x[0], dy = x[1] - cps[i].x[1], dz = x[2] - cps[i].x[2];
                if (dx * dx + dy * dy + dz * dz < hit * hit) return i;
            }
            if (rho < ctx->eps) return -1;
            if (sign > 0.0 ? rho > rho_stop : rho < rho_stop) return -1;
            if (gn < 1e-14) {   // stalled on a critical point: take the nearest one of the right type
                int best = -1; double bd = 1e-2;
                for (int i = 0; i < ncp; ++i) {
                    if (cps[i].type != target) continue;
                    const double d = sqrt((x[0] - cps[i].x[0]) * (x[0] - cps[i].x[0]) + (x[1] - cps[i].x[1]) * (x[1] - cps[i].x[1]) + (x[2] - cps[i].x[2]) * (x[2] - cps[i].x[2]));
                    if (d < bd) { bd = d; best = i; }
                }
                return best;
            }
        }
        const double fac = err > 0 ? 0.9 * pow(tol / err, 0.2) : 4.0;
        step *= fmin(4.0, fmax(0.2, fac));
        step = fmin(step, 0.25);
        step = fmax(step, 1e-7);
    }
    return -1;
}

// ----------------------------------------------------------------------------------------------- degree

static double cpg_solid_angle(const double a[3], const double b[3], const double c[3]) {
    const double bc[3] = { b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0] };
    const double num = a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2];
    const double den = 1.0 + (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) + (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) + (c[0] * a[0] + c[1] * a[1] + c[2] * a[2]);
    return 2.0 * atan2(num, den);
}

// Brouwer degree of grad rho on the boundary of [lo,hi] (solid angle of grad/|grad| over a triangulated surface).
static double cpg_degree(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double lo[3], const double hi[3], int n, double* out_min_grad) {
    double total = 0.0, gmin = DBL_MAX;
    double* grid = (double*)md_alloc(md_get_heap_allocator(), sizeof(double) * 3 * (n + 1) * (n + 1));
    for (int axis = 0; axis < 3; ++axis) {
        const int a = (axis + 1) % 3, b = (axis + 2) % 3;     // e_a x e_b = e_axis
        for (int side = -1; side <= 1; side += 2) {
            for (int u = 0; u <= n; ++u) for (int v = 0; v <= n; ++v) {
                double x[3];
                x[a] = lo[a] + (hi[a] - lo[a]) * u / n;
                x[b] = lo[b] + (hi[b] - lo[b]) * v / n;
                x[axis] = side > 0 ? hi[axis] : lo[axis];
                double rho, g[3];
                cpg_point(ctx, sc, x, 1, &rho, g, NULL);
                const double gn = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
                gmin = fmin(gmin, gn);
                double* d = grid + 3 * (u * (n + 1) + v);
                for (int k = 0; k < 3; ++k) d[k] = gn > 0 ? g[k] / gn : 0.0;
            }
            double face = 0.0;
            for (int u = 0; u < n; ++u) for (int v = 0; v < n; ++v) {
                const double* p00 = grid + 3 * (u * (n + 1) + v);
                const double* p10 = grid + 3 * ((u + 1) * (n + 1) + v);
                const double* p01 = grid + 3 * (u * (n + 1) + v + 1);
                const double* p11 = grid + 3 * ((u + 1) * (n + 1) + v + 1);
                face += cpg_solid_angle(p00, p10, p11) + cpg_solid_angle(p00, p11, p01);
            }
            total += side > 0 ? face : -face;
        }
    }
    md_free(md_get_heap_allocator(), grid, sizeof(double) * 3 * (n + 1) * (n + 1));
    *out_min_grad = gmin;
    return total / (4.0 * 3.14159265358979323846);
}

// ----------------------------------------------------------------------------------------------- driver

static int cpg_cp_cmp(const void* pa, const void* pb) {
    static const int order[5] = {9, 0, 1, 3, 2};   // MAX, SPLIT(BCP), JOIN(RCP), MIN(CCP)
    const cpg_cp_t* a = (const cpg_cp_t*)pa;
    const cpg_cp_t* b = (const cpg_cp_t*)pb;
    if (order[a->type] != order[b->type]) return order[a->type] < order[b->type] ? -1 : 1;
    if (a->rho != b->rho) return a->rho > b->rho ? -1 : 1;
    for (int k = 0; k < 3; ++k) if (a->x[k] != b->x[k]) return a->x[k] < b->x[k] ? -1 : 1;
    return 0;
}

static int cpg_find(int* parent, int i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
}

// ----------------------------------------------------------------------------------------------- workers

enum { CPG_DISCARD = 0, CPG_SPLIT = 1, CPG_UNRESOLVED = 2, CPG_ROOT = 3 };
#define CPG_KIND(k)            ((k) & 0xFF)
#define CPG_DISCARD_CHILDREN   (CPG_DISCARD | (1 << 16))   // split, but every child is excluded (bit 16: for the statistics)

typedef struct cpg_rootrec_t {
    size_t   box;
    cpg_cp_t cp;
} cpg_rootrec_t;

// Work items handed out one at a time, in order, to whichever worker asks next: for jobs whose items
// differ a lot in cost and are few (separatrices, polish). Every result is written to the item's own
// slot, so which thread takes which item changes nothing.
typedef struct cpg_queue_t {
    md_mutex_t mutex;
    int        next;
    int        count;
} cpg_queue_t;

static int cpg_queue_pop(cpg_queue_t* q) {
    md_mutex_lock(&q->mutex);
    const int i = q->next < q->count ? q->next++ : -1;
    md_mutex_unlock(&q->mutex);
    return i;
}

typedef struct cpg_polish_t {
    int      r;                         // cpg_polish_root's result
    cpg_cp_t cp;
} cpg_polish_t;

typedef struct cpg_worker_t {
    const cpg_ctx_t* ctx;
    cpg_scratch_t sc;
    int tid, nthreads;
    int job;                            // 0: cube sweep, 1: separatrices, 2: polish
    volatile int32_t* cancel;
    cpg_queue_t* queue;                 // jobs 1, 2 and 3
    // sweep
    const cpg_box_t* boxes;
    size_t nboxes;
    uint32_t* kind;                     // one slot per cube, written by exactly one worker
    double h_min;
    md_array(cpg_rootrec_t) roots;      // this worker's certified roots
    uint64_t inflated;
    // separatrices: item t traces path traces[t] % 2 of saddle saddles[traces[t] / 2]
    const cpg_cp_t* cps;
    int ncp;
    const int* saddles;
    const int* traces;
    int* ends;                          // two per saddle
    uint64_t sep_evals, sep_aos, sep_block_aos;   // point evaluations of this worker's paths, their local AOs summed, the block's
    uint32_t sep_max_evals, sep_open;   // the most for one path; paths that reached no CP
    double sep_rho_lo, sep_rho_hi;      // the lowest minimum's density and the highest maximum's (cpg_trace)
    // polish
    const cpg_box_t* pboxes;
    cpg_polish_t* pres;
} cpg_worker_t;

// The root of a cube Krawczyk certified (unique in the 1.5x inflated cube): double Newton from the cube
// centre. Returns 1 if the root belongs to this cube (written to out_cp), 0 if it lies outside it or
// below rho_min, -1 if Newton failed.
static int cpg_polish_root(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const cpg_box_t* box, cpg_cp_t* out_cp) {
    cpg_cp_t cp = {0};
    MEMCPY(cp.x, box->c, sizeof(cp.x));
    if (!cpg_newton(ctx, sc, cp.x, &cp.rho, cp.ev, cp.evec)) return -1;
    // Ownership is decided on the ROOT, not per cube: a CP on a shared face or corner (common by
    // symmetry, e.g. a ring CP at the origin) is polished by several cubes whose Newton results differ in
    // the last bits, so a half-open test per cube can reject it everywhere. Accept it in the closed cube
    // plus a tolerance and deduplicate afterwards. The Krawczyk certificate makes the root unique in the
    // inflated cube, so a duplicate can only be the same CP.
    const double tol = 1e-9;
    bool inside = true;
    for (int k = 0; k < 3; ++k) inside &= (cp.x[k] >= box->c[k] - box->h - tol) && (cp.x[k] <= box->c[k] + box->h + tol);
    if (!inside || cp.rho < ctx->eps) return 0;
    cp.h = box->h;
    cp.type = cpg_classify(cp.ev);
    *out_cp = cp;
    return 1;
}

// Children of a cube that splits which the cube's own expansion already excludes: g + A d + 1/2 T[d,d]
// re-centred at the child's centre, the quadratic term over the child, and the cube's remainder bound,
// which holds on all of the cube; then the exclusion tests (directions, Newton distance). Returns the mask
// of the children that still need an evaluation (bit s: child s, + side in x if s & 1, y if s & 2, z if
// s & 4). On the test systems 27-38% of all children are excluded this way, and no child excluded so was
// left undecided by its own evaluation.
static uint32_t cpg_child_mask(const cpg_eval_t* P, double h) {
    const double hc = 0.5 * h;
    uint32_t mask = 0;
    for (int s = 0; s < 8; ++s) {
        const double d[3] = { (s & 1) ? hc : -hc, (s & 2) ? hc : -hc, (s & 4) ? hc : -hc };
        cpg_eval_t C;
        MEMSET(&C, 0, sizeof(C));
        for (int i = 0; i < 3; ++i) {
            double gi = P->g[i];
            for (int j = 0; j < 3; ++j) {
                gi += P->A[i][j] * d[j];
                for (int l = 0; l < 3; ++l) gi += 0.5 * P->T[i][j][l] * d[j] * d[l];
            }
            C.g[i] = gi;
            for (int j = 0; j < 3; ++j) {
                double a = P->A[i][j];
                for (int l = 0; l < 3; ++l) a += P->T[i][j][l] * d[l];
                C.A[i][j] = a;
            }
            C.quadT[i] = 0.25 * P->quadT[i];      // 1/2 sum |T| hc^2
            C.r[i] = P->r[i];
        }
        MEMCPY(C.T, P->T, sizeof(C.T));           // a quadratic re-centred: the same T (the exact quadratic test uses it)
        double val[3], vec[3][3], M[3][3];
        MEMCPY(M, C.A, sizeof(M));
        cpg_eigen_sym3(val, vec, M);
        if (!cpg_exclude(&C, hc, val, vec)) mask |= 1u << s;
    }
    return mask;
}

// Decides one cube. Writes the polished root when the result is CPG_ROOT. A split comes with the mask of
// the children to evaluate in bits 8..15 (cpg_child_mask); a split with no child left is a discard.
static int cpg_process_box(cpg_worker_t* w, const cpg_box_t* box, cpg_cp_t* out_cp) {
    const cpg_ctx_t* ctx = w->ctx;
    cpg_eval_t ev, evi;
    cpg_box_eval(ctx, &w->sc, box->c, box->h, &ev);
    if (ev.rho_up < ctx->eps) return CPG_DISCARD;
    double val[3], vec[3][3], Y[3][3];
    cpg_eigen_sym3(val, vec, ev.A);
    if (cpg_exclude(&ev, box->h, val, vec)) return CPG_DISCARD;
    const bool invertible = cpg_inverse_sym3(Y, val, vec);
    int kr = invertible ? cpg_krawczyk(ev.g, Y, ev.dH, box->h) : 0;
    if (kr == -1) return CPG_DISCARD;
    if (invertible && kr == 0) {
        const double hi = 1.5 * box->h;
        w->inflated++;
        cpg_box_eval(ctx, &w->sc, box->c, hi, &evi);   // centre data is identical; only the enclosure differs
        kr = cpg_krawczyk(ev.g, Y, evi.dH, hi);
        if (kr == -1) return CPG_DISCARD;
    }
    if (kr == 1) {
        const int r = cpg_polish_root(ctx, &w->sc, box, out_cp);
        if (r == 1) return CPG_ROOT;
        if (r == 0) return CPG_DISCARD;
        // Newton failed despite the certificate (should not happen): treat the cube as undecided
    }
    if (box->h <= w->h_min) return CPG_UNRESOLVED;
    const uint32_t children = cpg_child_mask(&ev, box->h);
    return children ? (int)(CPG_SPLIT | (children << 8)) : CPG_DISCARD_CHILDREN;
}

static void cpg_worker_entry(void* data) {
    cpg_worker_t* w = (cpg_worker_t*)data;
    md_allocator_i* heap = md_get_heap_allocator();
    if (w->job == 0) {
        // fixed interleaved chunks: the result of every cube is independent of which thread computed it
        for (size_t chunk = (size_t)w->tid; chunk * CPG_CHUNK < w->nboxes; chunk += (size_t)w->nthreads) {
            if (w->cancel && *w->cancel) return;
            const size_t end = (chunk + 1) * CPG_CHUNK < w->nboxes ? (chunk + 1) * CPG_CHUNK : w->nboxes;
            for (size_t bi = chunk * CPG_CHUNK; bi < end; ++bi) {
                cpg_rootrec_t rec;
                const int k = cpg_process_box(w, &w->boxes[bi], &rec.cp);
                w->kind[bi] = (uint32_t)k;
                if (CPG_KIND(k) == CPG_ROOT) { rec.box = bi; md_array_push(w->roots, rec, heap); }
            }
        }
    } else if (w->job == 1) {
        for (int t; (t = cpg_queue_pop(w->queue)) >= 0; ) {
            if (w->cancel && *w->cancel) return;
            const int i = w->traces[t] / 2, e = w->traces[t] % 2;
            const cpg_cp_t* cp = &w->cps[w->saddles[i]];
            const bool bcp = cp->type == MD_TOPO_SPLIT_SADDLE;
            const int col = bcp ? 2 : 0;                    // the unique eigen-direction
            const double sign = bcp ? 1.0 : -1.0;           // ascend to maxima / descend to minima
            const int target = bcp ? MD_TOPO_MAXIMUM : MD_TOPO_MINIMUM;
            double dir[3];
            for (int k = 0; k < 3; ++k) dir[k] = (e ? -1.0 : 1.0) * cp->evec[k][col];
            const uint64_t e0 = w->sc.trace_evals, a0 = w->sc.trace_aos, b0 = w->sc.trace_block_aos;
            w->ends[2 * i + e] = cpg_trace(w->ctx, &w->sc, cp->x, dir, sign, w->cps, w->ncp, target, bcp ? w->sep_rho_hi : w->sep_rho_lo);
            const uint64_t ne = w->sc.trace_evals - e0;
            w->sep_evals += ne;
            w->sep_aos += w->sc.trace_aos - a0;
            w->sep_block_aos += w->sc.trace_block_aos - b0;
            if (ne > w->sep_max_evals) w->sep_max_evals = (uint32_t)MIN(ne, UINT32_MAX);
            if (w->ends[2 * i + e] < 0) w->sep_open++;
        }
    } else {
        for (int t; (t = cpg_queue_pop(w->queue)) >= 0; ) {
            if (w->cancel && *w->cancel) return;
            w->pres[t].r = cpg_polish_root(w->ctx, &w->sc, &w->pboxes[t], &w->pres[t].cp);
        }
    }
}

// Runs every worker; worker 0 on the calling thread. A worker whose thread cannot be created runs inline.
static void cpg_run_workers(cpg_worker_t* workers, int n) {
    md_thread_t* th[CPG_MAX_THREADS] = {0};
    for (int t = 1; t < n; ++t) th[t] = md_thread_create(cpg_worker_entry, &workers[t]);
    cpg_worker_entry(&workers[0]);
    for (int t = 1; t < n; ++t) {
        if (th[t]) md_thread_join(th[t]);
        else cpg_worker_entry(&workers[t]);
    }
}

static int cpg_rootrec_cmp(const void* pa, const void* pb) {
    const cpg_rootrec_t* a = (const cpg_rootrec_t*)pa;
    const cpg_rootrec_t* b = (const cpg_rootrec_t*)pb;
    return a->box < b->box ? -1 : (a->box > b->box ? 1 : 0);
}

static void cpg_scratch_alloc(cpg_scratch_t* sc, int N, int nshell, md_allocator_i* heap) {
    sc->L      = (int*)md_alloc(heap, sizeof(int) * N);
    sc->shells = (int*)md_alloc(heap, sizeof(int) * (nshell + 1));
    sc->V      = (double*)md_alloc(heap, sizeof(double) * N * CPG_NV);
    sc->Dv     = (double*)md_alloc(heap, sizeof(double) * N * CPG_NDV);
    sc->E      = (double*)md_alloc(heap, sizeof(double) * N * 16);
    sc->Q      = (double*)md_alloc(heap, sizeof(double) * N * 8);
    sc->tc_L      = (int*)md_alloc(heap, sizeof(int) * N);
    sc->tc_shells = (int*)md_alloc(heap, sizeof(int) * (nshell + 1));
    sc->tc_B      = NULL;
    sc->tc_cap    = 0;
    sc->tc_valid  = false;
}

static void cpg_scratch_free(cpg_scratch_t* sc, int N, int nshell, md_allocator_i* heap) {
    if (sc->tc_B) md_free(heap, sc->tc_B, sizeof(double) * sc->tc_cap);
    md_free(heap, sc->tc_shells, sizeof(int) * (nshell + 1));
    md_free(heap, sc->tc_L, sizeof(int) * N);
    md_free(heap, sc->Q,  sizeof(double) * N * 8);
    md_free(heap, sc->E,  sizeof(double) * N * 16);
    md_free(heap, sc->Dv, sizeof(double) * N * CPG_NDV);
    md_free(heap, sc->V,  sizeof(double) * N * CPG_NV);
    md_free(heap, sc->shells, sizeof(int) * (nshell + 1));
    md_free(heap, sc->L,  sizeof(int) * N);
}


// ----------------------------------------------------------------------------------------------- driver
// A run is set up once (basis, screening, search domain, root lattice), swept, then finished (order,
// separatrices, clusters, graph). The sweep is either the CPU level loop below or the GPU one
// (md_topo_compute_extremum_graph_gto_gpu), which hands the cubes fp32 cannot decide to the CPU loop.
// Both feed the same root and unresolved lists, so everything after the sweep is shared.

typedef struct cpg_run_t {
    const md_topo_gto_desc_t* desc;
    md_allocator_i*           heap;
    cpg_ctx_t*                ctx;
    cpg_worker_t*             workers;
    int                       nthreads;
    double                    h_min;
    bool                      ok;
    md_array(cpg_cp_t)        cps;
    md_array(cpg_box_t)       cur;
    md_array(cpg_box_t)       nxt;
    md_array(cpg_box_t)       unres;
    md_array(uint32_t)        kind;
    md_array(cpg_rootrec_t)   roots;
    uint64_t*                 tile_need;    // ntiles x tile_words: bit J of row I, both orders (cpg_run_init)
    int                       ntiles, tile_words;
    md_topo_gto_info_t        info;
} cpg_run_t;

static bool cpg_run_cancelled(cpg_run_t* R) {
    if (R->desc->cancel && *R->desc->cancel) R->info.cancelled = true;
    return R->info.cancelled;
}

// Accepts a polished root unless it is one already accepted (the same CP reached from another cube).
static void cpg_run_accept(cpg_run_t* R, const cpg_cp_t* cp) {
    for (size_t i = 0; i < md_array_size(R->cps); ++i) {
        const double dx = R->cps[i].x[0] - cp->x[0], dy = R->cps[i].x[1] - cp->x[1], dz = R->cps[i].x[2] - cp->x[2];
        if (dx * dx + dy * dy + dz * dz < 1e-16) return;
    }
    md_array_push(R->cps, *cp, R->heap);
}

// Upper bound of rho at distance >= pad from every centre: T^T |D| T, T_i = cpg_tail_ao(i, order 0, pad).
// With M (K x K, |D| summed over the AO keys, the class's AO offset + Cartesian index), over key pairs;
// without, over AO pairs. T: K or N doubles of scratch.
static double cpg_pad_bound(const cpg_ctx_t* ctx, double pad, const double* M, int K, int ncls, const int* cls_rep,
                           const int* cls_key, double* T) {
    double b = 0.0;
    if (M) {
        for (int c = 0; c < ncls; ++c) {
            const cpg_shell_t* sh = &ctx->shell[cls_rep[c]];
            for (int ci = 0; ci < sh->ncart; ++ci) T[cls_key[c] + ci] = cpg_tail_ao(ctx, sh, ci, 0, pad);
        }
        for (int k = 0; k < K; ++k) {
            if (T[k] == 0.0) continue;
            double s = 0.0;
            for (int l = 0; l < K; ++l) s += M[(size_t)k * K + l] * T[l];
            b += T[k] * s;
        }
    } else {
        const int N = ctx->nao;
        for (int si = 0; si < ctx->nshell; ++si) {
            const cpg_shell_t* sh = &ctx->shell[si];
            for (int ci = 0; ci < sh->ncart; ++ci) T[sh->ao_offset + ci] = cpg_tail_ao(ctx, sh, ci, 0, pad);
        }
        for (int i = 0; i < N; ++i) {
            if (T[i] == 0.0) continue;
            for (int j = 0; j < N; ++j) b += T[i] * fabs(ctx->D[(size_t)i * N + j]) * T[j];
        }
    }
    return b;
}

static void cpg_run_init(cpg_run_t* R, const md_topo_gto_desc_t* desc) {
    MEMSET(R, 0, sizeof(*R));
    R->desc = desc;
    R->heap = md_get_heap_allocator();
    md_allocator_i* heap = R->heap;
    const md_gto_basis_t* basis = desc->basis;

    cpg_ctx_t* ctx = (cpg_ctx_t*)md_alloc(heap, sizeof(cpg_ctx_t));
    MEMSET(ctx, 0, sizeof(*ctx));
    R->ctx = ctx;
    cpg_init_tables(ctx);
    ctx->eps = desc->rho_min > 0.0 ? desc->rho_min : 1.0e-4;
    R->h_min = desc->h_min > 0.0 ? desc->h_min : 1.0e-4;
    const double h_root = desc->h_root > 0.0 ? desc->h_root : 1.0;
    const size_t stride = desc->atom_xyz_stride ? desc->atom_xyz_stride : sizeof(float) * 3;

    int nthreads = (int)desc->num_threads;
    if (nthreads <= 0) {
        md_os_sys_info_t si = {0};
        nthreads = md_os_sys_info_query(&si) && si.num_virtual_cores > 0 ? si.num_virtual_cores : 1;
    }
    if (nthreads > CPG_MAX_THREADS) nthreads = CPG_MAX_THREADS;
    R->nthreads = nthreads;
    R->info.num_threads = (uint32_t)nthreads;

    // --- basis -> shells / AOs (double precision copies)
    ctx->nshell = (int)basis->num_shells;
    ctx->nao = (int)md_gto_basis_num_ao(basis);
    ctx->shell  = (cpg_shell_t*)md_alloc(heap, sizeof(cpg_shell_t) * ctx->nshell);
    ctx->alpha  = (double*)md_alloc(heap, sizeof(double) * basis->num_primitives);
    ctx->coeff  = (double*)md_alloc(heap, sizeof(double) * basis->num_primitives);
    ctx->ao_ijk = (int(*)[3])md_alloc(heap, sizeof(int[3]) * ctx->nao);
    ctx->ao_nrm = (double*)md_alloc(heap, sizeof(double) * ctx->nao);
    for (uint32_t p = 0; p < basis->num_primitives; ++p) { ctx->alpha[p] = basis->alpha[p]; ctx->coeff[p] = basis->coeff[p]; }
    uint32_t ao = 0;
    bool ok = true;
    for (int si = 0; si < ctx->nshell; ++si) {
        const md_gto_shell_t* bs = &basis->shells[si];
        cpg_shell_t* s = &ctx->shell[si];
        if (bs->l > CPG_LMAX || bs->num_primitives > CPG_MAXPRIM) { ok = false; break; }
        const float* xyz = (const float*)((const uint8_t*)desc->atom_xyz + bs->atom_idx * stride);
        s->A[0] = xyz[0]; s->A[1] = xyz[1]; s->A[2] = xyz[2];
        s->l = (int)bs->l;
        s->ncart = (int)md_gto_num_cart_ao(bs->l);
        s->prim_offset = bs->primitive_offset;
        s->num_prims = bs->num_primitives;
        s->ao_offset = ao;
        int ci = 0;
        for (int i = s->l; i >= 0; --i) for (int j = s->l - i; j >= 0; --j) {
            const int k = s->l - i - j;
            ctx->ao_ijk[ao + ci][0] = i; ctx->ao_ijk[ao + ci][1] = j; ctx->ao_ijk[ao + ci][2] = k;
            ctx->ao_nrm[ao + ci] = cpg_nrm(i, j, k);
            ci++;
        }
        ao += s->ncart;
    }
    if (!ok) {
        MD_LOG_ERROR("md_topo_compute_extremum_graph_gto: unsupported shell (l > %d or more than %d primitives)", CPG_LMAX, CPG_MAXPRIM);
    }
    ctx->D = desc->density_matrix;
    const int N = ctx->nao;
    const int NS = ctx->nshell;

    R->workers = (cpg_worker_t*)md_alloc(heap, sizeof(cpg_worker_t) * nthreads);
    MEMSET(R->workers, 0, sizeof(cpg_worker_t) * nthreads);
    for (int t = 0; t < nthreads; ++t) {
        R->workers[t].ctx = ctx;
        R->workers[t].tid = t;
        R->workers[t].nthreads = nthreads;
        R->workers[t].cancel = desc->cancel;
        R->workers[t].h_min = R->h_min;
        cpg_scratch_alloc(&R->workers[t].sc, N, NS, heap);
    }
    R->ok = ok;
    if (!ok) return;

    // --- radial classes: the tail bounds (cpg_tail_shell, cpg_tail_ao) of a shell depend only on its l and
    //     contraction, which a large system repeats for every atom of an element, so they are computed once
    //     per class. AO key: class AO offset + Cartesian index (the AO maxima and domain pad tails).
    int* shell_cls = (int*)md_alloc(heap, sizeof(int) * NS);
    int* cls_rep   = (int*)md_alloc(heap, sizeof(int) * NS);       // a shell of each class
    int* cls_key   = (int*)md_alloc(heap, sizeof(int) * (NS + 1)); // first AO key of each class
    int ncls = 0;
    {
        uint32_t hcap = 16;
        while (hcap < 2u * (uint32_t)NS) hcap *= 2;
        int* htab = (int*)md_alloc(heap, sizeof(int) * hcap);
        for (uint32_t i = 0; i < hcap; ++i) htab[i] = -1;
        for (int si = 0; si < NS; ++si) {
            const cpg_shell_t* sh = &ctx->shell[si];
            const size_t pb = sizeof(double) * sh->num_prims;
            uint64_t hv = 1469598103934665603ull;
            const uint8_t* bytes[2] = { (const uint8_t*)(ctx->alpha + sh->prim_offset), (const uint8_t*)(ctx->coeff + sh->prim_offset) };
            hv = (hv ^ (uint64_t)sh->l) * 1099511628211ull;
            hv = (hv ^ (uint64_t)sh->num_prims) * 1099511628211ull;
            for (int a = 0; a < 2; ++a) for (size_t b = 0; b < pb; ++b) hv = (hv ^ bytes[a][b]) * 1099511628211ull;
            uint32_t slot = (uint32_t)(hv ^ (hv >> 32)) & (hcap - 1);
            int c = -1;
            for (; htab[slot] >= 0; slot = (slot + 1) & (hcap - 1)) {
                const cpg_shell_t* rs = &ctx->shell[cls_rep[htab[slot]]];
                if (rs->l == sh->l && rs->num_prims == sh->num_prims &&
                    memcmp(ctx->alpha + rs->prim_offset, ctx->alpha + sh->prim_offset, pb) == 0 &&
                    memcmp(ctx->coeff + rs->prim_offset, ctx->coeff + sh->prim_offset, pb) == 0) { c = htab[slot]; break; }
            }
            if (c < 0) { c = ncls++; cls_rep[c] = si; htab[slot] = c; }
            shell_cls[si] = c;
        }
        md_free(heap, htab, sizeof(int) * hcap);
    }
    cls_key[0] = 0;
    for (int c = 0; c < ncls; ++c) cls_key[c + 1] = cls_key[c] + ctx->shell[cls_rep[c]].ncart;
    const int K = cls_key[ncls];
    int* ao_key = (int*)md_alloc(heap, sizeof(int) * N);
    for (int si = 0; si < NS; ++si) {
        const cpg_shell_t* sh = &ctx->shell[si];
        for (int ci = 0; ci < sh->ncart; ++ci) ao_key[sh->ao_offset + ci] = cls_key[shell_cls[si]] + ci;
    }

    // --- global AO maxima (order 0..2) and |D|-weighted row sums -> screening threshold tau
    double* phimax = (double*)md_alloc(heap, sizeof(double) * K * 3);
    for (int c = 0; c < ncls; ++c) {
        const cpg_shell_t* sh = &ctx->shell[cls_rep[c]];
        for (int ci = 0; ci < sh->ncart; ++ci) for (int o = 0; o < 3; ++o) phimax[(cls_key[c] + ci) * 3 + o] = cpg_tail_ao(ctx, sh, ci, o, 0.0);
    }
    // Ri[3 i + o] = sum_j |D_ij| phimax_j^(o). A shell s whose AOs (and their derivatives, order <= 2) are
    // below t at a point adds at most 2 t sum_{i in s} (Ri0, Ri0 + Ri1, Ri0 + 2 Ri1 + Ri2) to rho, grad and
    // the Hessian there (the pairs it is in, both orders): its tail weights tw. Every shell at once, each
    // below tau: the global tails.
    double* Ri = (double*)md_alloc(heap, sizeof(double) * 3 * (size_t)N);
    MEMSET(Ri, 0, sizeof(double) * 3 * (size_t)N);
    for (int i = 0; i < N; ++i) {
        const double* Drow = ctx->D + (size_t)i * N;
        double r0 = 0.0, r1 = 0.0, r2 = 0.0;
        for (int j = 0; j < N; ++j) {
            const double d = fabs(Drow[j]);
            const double* pm = phimax + ao_key[j] * 3;
            r0 += d * pm[0]; r1 += d * pm[1]; r2 += d * pm[2];
        }
        Ri[3 * i + 0] = r0; Ri[3 * i + 1] = r1; Ri[3 * i + 2] = r2;
    }
    md_free(heap, phimax, sizeof(double) * K * 3);
    double Rs[3] = {0, 0, 0};
    for (int si = 0; si < NS; ++si) {
        cpg_shell_t* sh = &ctx->shell[si];
        double R0 = 0.0, R1 = 0.0, R2 = 0.0;
        for (int ci = 0; ci < sh->ncart; ++ci) {
            const double* r = Ri + 3 * ((size_t)sh->ao_offset + ci);
            R0 += r[0]; R1 += r[1]; R2 += r[2];
        }
        sh->tw[0] = 2.0 * R0;
        sh->tw[1] = 2.0 * (R0 + R1);
        sh->tw[2] = 2.0 * (R0 + 2.0 * R1 + R2);
        Rs[0] += R0; Rs[1] += R1; Rs[2] += R2;
        sh->cls = shell_cls[si];
    }
    md_free(heap, Ri, sizeof(double) * 3 * (size_t)N);
    const double tail_target = 1.0e-12;
    ctx->tau = tail_target / (2.0 * (Rs[0] + 2.0 * Rs[1] + Rs[2]) + 1e-300);
    ctx->tail_rho = 2.0 * ctx->tau * Rs[0];
    ctx->tail_g   = 2.0 * ctx->tau * (Rs[0] + Rs[1]);
    ctx->tail_H   = 2.0 * ctx->tau * (Rs[0] + 2.0 * Rs[1] + Rs[2]);
    // rung radii per class, by bisection on the (monotone) tail bound; rung 0 is the far radius
    ctx->ncls = ncls;
    ctx->rung_r2 = (double*)md_alloc(heap, sizeof(double) * (size_t)ncls * CPG_RUNGS);
    for (int m = 0; m < CPG_RUNGS; ++m) ctx->rung_tau[m] = ctx->tau * pow(10.0, 0.5 * m);
    for (int c = 0; c < ncls; ++c) {
        const cpg_shell_t* sh = &ctx->shell[cls_rep[c]];
        double lo = 0.0, hi = 1.0;
        while (cpg_tail_shell(ctx, sh, 2, hi) > ctx->tau && hi < 1e3) hi *= 2.0;
        for (int it = 0; it < 50; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (cpg_tail_shell(ctx, sh, 2, mid) > ctx->tau) lo = mid; else hi = mid;
        }
        ctx->rung_r2[(size_t)c * CPG_RUNGS] = hi * hi;
        // the others down a grid of 1/32 Bohr from the far radius: rung m's radius is the last grid point
        // (from above) where the bound is still below its threshold (the bound is monotone in the distance)
        double d = hi;
        for (int m = 1; m < CPG_RUNGS; ++m) {
            while (d > 0.03125 && cpg_tail_shell(ctx, sh, 2, d - 0.03125) <= ctx->rung_tau[m]) d -= 0.03125;
            ctx->rung_r2[(size_t)c * CPG_RUNGS + m] = d * d;
        }
    }
    for (int si = 0; si < NS; ++si) ctx->shell[si].radius = sqrt(ctx->rung_r2[(size_t)ctx->shell[si].cls * CPG_RUNGS]);
    ctx->rung_loc = 0;          // until the domain is known (below)

    // --- shells that meet: s and t are both in the local list of some cube, batch or point only if
    //     |A_s - A_t| <= r_s + r_t + 2 rho, rho = sqrt(3) e the largest half-diagonal of a region screened
    //     (e = 2.5 h_root: a GPU batch, P +- 2.5 hc with hc <= h_root; a CPU cube's inflated enclosure is
    //     1.5 h). Generous by 1e-6 for the GPU's conservative screen. n_bound: no local list holds more
    //     AOs (every shell of a list meets every other one of it).
    //     The tile map: bit J of row I when some AO of tile I and some AO of tile J belong to shells that
    //     meet (both orders), kept for the GPU's D.
    const int ntiles = (N + CPG_DTILE - 1) / CPG_DTILE;
    const int ntw = (ntiles + 63) / 64;
    uint64_t* tile_need = (uint64_t*)md_alloc(heap, sizeof(uint64_t) * (size_t)ntiles * ntw);
    MEMSET(tile_need, 0, sizeof(uint64_t) * (size_t)ntiles * ntw);
    R->tile_need = tile_need;
    R->ntiles = ntiles;
    R->tile_words = ntw;
    uint32_t n_bound = 0;
    {
        const double rho2 = 2.0 * sqrt(3.0) * 2.5 * h_root;
        uint32_t* nmeet = (uint32_t*)md_alloc(heap, sizeof(uint32_t) * NS);
        for (int si = 0; si < NS; ++si) nmeet[si] = (uint32_t)ctx->shell[si].ncart;
        for (int si = 0; si < NS; ++si) {
            const cpg_shell_t* a = &ctx->shell[si];
            const int ia0 = (int)a->ao_offset / CPG_DTILE, ia1 = ((int)a->ao_offset + a->ncart - 1) / CPG_DTILE;
            for (int ti = si; ti < NS; ++ti) {
                const cpg_shell_t* b = &ctx->shell[ti];
                if (ti != si) {
                    const double dx = a->A[0] - b->A[0], dy = a->A[1] - b->A[1], dz = a->A[2] - b->A[2];
                    const double lim = (a->radius + b->radius + rho2) * (1.0 + 1e-6);
                    if (dx * dx + dy * dy + dz * dz > lim * lim) continue;
                    nmeet[si] += (uint32_t)b->ncart;
                    nmeet[ti] += (uint32_t)a->ncart;
                }
                const int ib0 = (int)b->ao_offset / CPG_DTILE, ib1 = ((int)b->ao_offset + b->ncart - 1) / CPG_DTILE;
                for (int I = ia0; I <= ia1; ++I) for (int J = ib0; J <= ib1; ++J) {
                    tile_need[(size_t)I * ntw + J / 64] |= 1ull << (J % 64);
                    tile_need[(size_t)J * ntw + I / 64] |= 1ull << (I % 64);
                }
            }
        }
        for (int si = 0; si < NS; ++si) n_bound = MAX(n_bound, nmeet[si]);
        md_free(heap, nmeet, sizeof(uint32_t) * NS);
    }
    R->info.local_ao_bound = n_bound;

    // --- search domain: outside AABB(atoms) +- pad every point is >= pad from every centre; rho <= T^T |D| T.
    //     T depends only on the AO key, so T^T |D| T = sum over key pairs of T_k T_l M_kl, M_kl the sum of
    //     |D_ij| over AO pairs with those keys: one pass over D, then K^2 per bisection step.
    double amin[3] = { DBL_MAX, DBL_MAX, DBL_MAX }, amax[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
    for (int si = 0; si < NS; ++si) for (int k = 0; k < 3; ++k) { amin[k] = fmin(amin[k], ctx->shell[si].A[k]); amax[k] = fmax(amax[k], ctx->shell[si].A[k]); }
    const bool by_key = (size_t)K * K <= (size_t)N * N && K <= 2048;
    const int NT = by_key ? K : N;
    double* M  = by_key ? (double*)md_alloc(heap, sizeof(double) * (size_t)K * K) : NULL;
    double* T0 = (double*)md_alloc(heap, sizeof(double) * NT);
    if (by_key) {
        MEMSET(M, 0, sizeof(double) * (size_t)K * K);
        for (int i = 0; i < N; ++i) {
            double* Mi = M + (size_t)ao_key[i] * K;
            for (int j = 0; j < N; ++j) Mi[ao_key[j]] += fabs(ctx->D[(size_t)i * N + j]);
        }
    }
    double plo = 0.0, phi = 1.0;
    for (;;) {
        const double b = cpg_pad_bound(ctx, phi, M, K, ncls, cls_rep, cls_key, T0);
        if (b < ctx->eps || phi > 1e3) break;
        plo = phi; phi *= 2.0;
    }
    for (int it = 0; it < 40; ++it) {
        const double mid = 0.5 * (plo + phi);
        const double b = cpg_pad_bound(ctx, mid, M, K, ncls, cls_rep, cls_key, T0);
        if (b < ctx->eps) phi = mid; else plo = mid;
    }
    md_free(heap, T0, sizeof(double) * NT);
    if (M) md_free(heap, M, sizeof(double) * (size_t)K * K);
    md_free(heap, ao_key, sizeof(int) * N);
    md_free(heap, cls_key, sizeof(int) * (NS + 1));
    md_free(heap, cls_rep, sizeof(int) * NS);
    md_free(heap, shell_cls, sizeof(int) * NS);
    const double pad = phi;
    R->info.domain_pad = pad;

    // --- the local radius: the highest rung whose region tails stay small. Sampled at points of a lattice
    //     over the domain (at most 16 per axis): at each, the tails a region there would carry with each
    //     rung as the local radius; the rung is the highest with the Hessian tail below CPG_SCREEN_TAIL
    //     rho_min everywhere sampled (a region adds its extent to the distances, so a batch can carry more).
    //     A choice of cost only: every region's tails are computed and charged whatever the rung. On mol
    //     (240 AOs, CPU) targets of 2.5e-7 .. 2.5e-3 rho_min took 150291 .. 148456 cubes, the same CPs: the
    //     shorter lists' remainder sums are tighter, and the tails stay far below what a decision needs. On
    //     a 9636-AO protein, 2.5e-5 and 2.5e-4 took 510 and 471 local AOs per GPU batch (737 with the far
    //     radius everywhere), the same CPs and edges bit for bit (3080, M4, and the CPU in double).
    {
        const double target = CPG_SCREEN_TAIL * ctx->eps;
        double ext = 0.0;
        for (int k = 0; k < 3; ++k) ext = fmax(ext, amax[k] - amin[k] + 2.0 * pad);
        const double sp = fmax(2.0 * h_root, ext / 16.0);
        int np[3];
        for (int k = 0; k < 3; ++k) np[k] = MAX(1, (int)ceil((amax[k] - amin[k] + 2.0 * pad) / sp));
        double worst[CPG_RUNGS] = {0};
        ctx->rung_loc = CPG_RUNGS - 1;     // cpg_screen_rung: every rung
        for (int z = 0; z < np[2]; ++z) for (int y = 0; y < np[1]; ++y) for (int x = 0; x < np[0]; ++x) {
            const double p[3] = { amin[0] - pad + (x + 0.5) * sp, amin[1] - pad + (y + 0.5) * sp, amin[2] - pad + (z + 0.5) * sp };
            double hist[CPG_RUNGS] = {0};
            for (int si = 0; si < NS; ++si) {
                const cpg_shell_t* sh = &ctx->shell[si];
                const double d2 = cpg_box_dist2(sh->A, p, 0.0);
                if (d2 >= ctx->rung_r2[(size_t)sh->cls * CPG_RUNGS]) continue;     // far
                const int m = cpg_screen_rung(ctx, sh, d2);
                if (m > 0) hist[m] += ctx->rung_tau[m] * sh->tw[2];
            }
            double acc = 0.0;
            for (int m = 1; m < CPG_RUNGS; ++m) { acc += hist[m]; worst[m] = fmax(worst[m], acc); }
        }
        int loc = 0;
        while (loc + 1 < CPG_RUNGS && worst[loc + 1] <= target) loc++;
        ctx->rung_loc = loc;
        R->info.screening_tail_sampled = worst[loc];
    }
    R->info.screening_rung = (uint32_t)ctx->rung_loc;
    for (int si = 0; si < NS; ++si) {
        cpg_shell_t* sh = &ctx->shell[si];
        sh->radius_trace = sqrt(ctx->rung_r2[(size_t)sh->cls * CPG_RUNGS + ctx->rung_loc]);
    }
    {
        // the local AO count at the atoms with the local and with the far radius (a point there): what the
        // rung saves where the work is
        uint64_t nl = 0, nf = 0, nat = 0;
        for (int si = 0; si < NS; ++si) {
            const cpg_shell_t* a = &ctx->shell[si];
            if (si > 0 && a->A[0] == ctx->shell[si - 1].A[0] && a->A[1] == ctx->shell[si - 1].A[1] && a->A[2] == ctx->shell[si - 1].A[2]) continue;
            nat++;
            for (int ti = 0; ti < NS; ++ti) {
                const cpg_shell_t* b = &ctx->shell[ti];
                const double d2 = cpg_box_dist2(b->A, a->A, 0.0);
                if (d2 < b->radius * b->radius) nf += (uint64_t)b->ncart;
                if (cpg_screen_rung(ctx, b, d2) < 0) nl += (uint64_t)b->ncart;
            }
        }
        R->info.screening_aos_local = nat ? (double)nl / (double)nat : 0.0;
        R->info.screening_aos_far = nat ? (double)nf / (double)nat : 0.0;
    }
    {
        double rl = 0.0, rf = 0.0;      // the most diffuse class: its local and far radius
        for (int c = 0; c < ncls; ++c) {
            const double f = ctx->rung_r2[(size_t)c * CPG_RUNGS];
            if (f > rf * rf) { rf = sqrt(f); rl = sqrt(ctx->rung_r2[(size_t)c * CPG_RUNGS + ctx->rung_loc]); }
        }
        R->info.screening_radius_local = rl;
        R->info.screening_radius_far = rf;
    }

    // --- root lattice, anchored to world multiples of 2 h_root
    const double w = 2.0 * h_root;
    int64_t i0[3], cnt[3];
    for (int k = 0; k < 3; ++k) {
        i0[k] = (int64_t)floor((amin[k] - pad) / w);
        const int64_t i1 = (int64_t)ceil((amax[k] + pad) / w);
        cnt[k] = i1 - i0[k];
    }
    for (int64_t z = 0; z < cnt[2]; ++z) for (int64_t y = 0; y < cnt[1]; ++y) for (int64_t x = 0; x < cnt[0]; ++x) {
        cpg_box_t b = { { (i0[0] + x + 0.5) * w, (i0[1] + y + 0.5) * w, (i0[2] + z + 0.5) * w }, h_root };
        md_array_push(R->cur, b, heap);
    }
}

// Level-synchronous sweep of R->cur on the CPU: cubes in parallel, results merged in cube order.
// The cubes need not share a size (the GPU sweep hands over cubes from several levels).
static void cpg_run_sweep_cpu(cpg_run_t* R) {
    md_allocator_i* heap = R->heap;
    cpg_worker_t* workers = R->workers;
    const int nthreads = R->nthreads;
    const md_tick_t t_start = md_tick_now();
    while (md_array_size(R->cur) > 0) {
        if (cpg_run_cancelled(R)) break;
        const size_t nb = md_array_size(R->cur);
        md_array_resize(R->kind, nb, heap);
        R->info.num_levels++;
        R->info.num_box_evals += nb;
        for (int t = 0; t < nthreads; ++t) {
            workers[t].job = 0;
            workers[t].boxes = R->cur;
            workers[t].nboxes = nb;
            workers[t].kind = R->kind;
            md_array_shrink(workers[t].roots, 0);
        }
        cpg_run_workers(workers, nthreads);
        if (cpg_run_cancelled(R)) break;

        // roots, in cube order, deduplicated against everything accepted so far
        md_array_shrink(R->roots, 0);
        for (int t = 0; t < nthreads; ++t) {
            md_array_push_array(R->roots, workers[t].roots, md_array_size(workers[t].roots), heap);
            R->info.num_inflated_evals += workers[t].inflated;
            workers[t].inflated = 0;
        }
        if (md_array_size(R->roots) > 1) qsort(R->roots, md_array_size(R->roots), sizeof(cpg_rootrec_t), cpg_rootrec_cmp);
        for (size_t r = 0; r < md_array_size(R->roots); ++r) cpg_run_accept(R, &R->roots[r].cp);

        md_array_shrink(R->nxt, 0);
        for (size_t bi = 0; bi < nb; ++bi) {
            const uint32_t k = R->kind[bi];
            if (k == CPG_DISCARD_CHILDREN) {
                R->info.num_children_skipped += 8;
            } else if (CPG_KIND(k) == CPG_UNRESOLVED) {
                md_array_push(R->unres, R->cur[bi], heap);
            } else if (CPG_KIND(k) == CPG_SPLIT) {
                const cpg_box_t box = R->cur[bi];
                const double hc = 0.5 * box.h;
                for (int s = 0; s < 8; ++s) {
                    if (!(k & (1u << (8 + s)))) { R->info.num_children_skipped++; continue; }
                    cpg_box_t child = { { box.c[0] + ((s & 1) ? hc : -hc), box.c[1] + ((s & 2) ? hc : -hc), box.c[2] + ((s & 4) ? hc : -hc) }, hc };
                    md_array_push(R->nxt, child, heap);
                }
            }
        }
        cpg_box_t* t = R->cur; R->cur = R->nxt; R->nxt = t;
    }
    const double ms = md_tick_to_milliseconds(md_tick_now() - t_start);
    R->info.ms_sweep_cpu += ms;
    R->info.ms_sweep += ms;
}

static bool cpg_run_finish(cpg_run_t* R, md_topo_extremum_graph_t* out_graph, md_topo_gto_info_t* out_info) {
    if (!R->ok) {
        if (out_info) *out_info = R->info;
        return false;
    }
    md_allocator_i* heap = R->heap;
    cpg_ctx_t* ctx = R->ctx;
    cpg_worker_t* workers = R->workers;
    const int nthreads = R->nthreads;
    const md_topo_gto_desc_t* desc = R->desc;
    md_topo_gto_info_t* info = &R->info;
    cpg_cp_t* cps = R->cps;
    cpg_box_t* unres = R->unres;

    // --- deterministic order
    const int ncp = (int)md_array_size(cps);
    if (ncp > 1) qsort(cps, ncp, sizeof(cpg_cp_t), cpg_cp_cmp);

    // --- separatrices -> edges (from saddle to extremum), traced in parallel
    md_array(md_topo_edge_t) edges = 0;
    md_tick_t t_phase = md_tick_now();
    if (desc->trace_separatrices && !info->cancelled) {
        md_array(int) saddles = 0;
        for (int i = 0; i < ncp; ++i) {
            if (cps[i].type == MD_TOPO_SPLIT_SADDLE || cps[i].type == MD_TOPO_JOIN_SADDLE) md_array_push(saddles, i, heap);
        }
        const int ns = (int)md_array_size(saddles);
        if (ns > 0) {
            int* ends = (int*)md_alloc(heap, sizeof(int) * 2 * ns);
            for (int i = 0; i < 2 * ns; ++i) ends[i] = -1;
            // One work item per path. A path can only end at a CP of its target type, so with none of
            // that type there is nothing to trace (a molecule without cage points: every RCP descent
            // would run out to rho_min, the longest paths of all, and give no edge). Descents first:
            // they are the longer ones, which a queue should start early.
            bool has_max = false, has_min = false;
            for (int i = 0; i < ncp; ++i) {
                has_max |= cps[i].type == MD_TOPO_MAXIMUM;
                has_min |= cps[i].type == MD_TOPO_MINIMUM;
            }
            md_array(int) traces = 0;
            for (int pass = 0; pass < 2; ++pass) {
                for (int i = 0; i < ns; ++i) {
                    const bool bcp = cps[saddles[i]].type == MD_TOPO_SPLIT_SADDLE;
                    if (bcp != (pass == 1) || !(bcp ? has_max : has_min)) continue;
                    md_array_push(traces, 2 * i, heap);
                    md_array_push(traces, 2 * i + 1, heap);
                }
            }
            double rho_lo = DBL_MAX, rho_hi = 0.0;
            for (int i = 0; i < ncp; ++i) {
                if (cps[i].type == MD_TOPO_MINIMUM) rho_lo = fmin(rho_lo, cps[i].rho);
                if (cps[i].type == MD_TOPO_MAXIMUM) rho_hi = fmax(rho_hi, cps[i].rho);
            }
            info->num_separatrix_paths = (uint32_t)md_array_size(traces);
            // the GPU sweep's tables and D are still on the device: the paths go there, in fp32; the CPU
            // takes those the GPU cannot (a local AO list beyond its scratch), or all if it fails
            md_array(int) cpu_traces = traces;
            md_array(int) rest = 0;
            cpg_queue_t queue = { .next = 0, .count = (int)md_array_size(cpu_traces) };
            md_mutex_init(&queue.mutex);       // in place: a CRITICAL_SECTION or pthread mutex must not be copied
            for (int t = 0; t < nthreads; ++t) {
                workers[t].job = 1;
                workers[t].queue = &queue;
                workers[t].cps = cps;
                workers[t].ncp = ncp;
                workers[t].saddles = saddles;
                workers[t].traces = cpu_traces;
                workers[t].ends = ends;
            }
            for (int t = 0; t < nthreads; ++t) {
                workers[t].sep_evals = workers[t].sep_aos = workers[t].sep_block_aos = 0;
                workers[t].sep_max_evals = workers[t].sep_open = 0;
                workers[t].sep_rho_lo = rho_lo * (1.0 - 1e-6);
                workers[t].sep_rho_hi = rho_hi * (1.0 + 1e-6);
            }
            if (queue.count > 0) cpg_run_workers(workers, MIN(nthreads, queue.count));
            md_mutex_destroy(&queue.mutex);
            md_array_free(rest, heap);
            uint64_t sep_aos = 0, sep_block_aos = 0;
            for (int t = 0; t < nthreads; ++t) {
                info->num_separatrix_evals += workers[t].sep_evals;
                sep_aos += workers[t].sep_aos;
                sep_block_aos += workers[t].sep_block_aos;
                info->max_separatrix_evals = MAX(info->max_separatrix_evals, workers[t].sep_max_evals);
                info->num_separatrix_open += workers[t].sep_open;
            }
            const uint64_t cpu_evals = info->num_separatrix_evals - (info->separatrices_on_gpu ? info->num_separatrix_evals_gpu : 0);
            info->separatrix_mean_local_aos = cpu_evals ? (double)sep_aos / (double)cpu_evals : 0.0;
            info->separatrix_mean_block_aos = cpu_evals ? (double)sep_block_aos / (double)cpu_evals : 0.0;
            md_array_free(traces, heap);
            cpg_run_cancelled(R);
            for (int i = 0; i < ns; ++i) {
                for (int e = 0; e < 2; ++e) {
                    const int to = ends[2 * i + e];
                    if (to < 0) continue;
                    if (e == 1 && to == ends[2 * i]) continue;
                    md_topo_edge_t edge = { (uint32_t)saddles[i], (uint32_t)to };
                    md_array_push(edges, edge, heap);
                }
            }
            md_free(heap, ends, sizeof(int) * 2 * ns);
        }
        md_array_free(saddles, heap);
    }

    info->ms_separatrices = md_tick_to_milliseconds(md_tick_now() - t_phase);
    t_phase = md_tick_now();

    // --- unresolved cubes -> clusters with their boundary degree
    const int nu = (int)md_array_size(unres);
    info->num_unresolved_boxes = (uint32_t)nu;
    if (nu > 0 && !info->cancelled) {
        int* parent = (int*)md_alloc(heap, sizeof(int) * nu);
        for (int i = 0; i < nu; ++i) parent[i] = i;
        for (int i = 0; i < nu; ++i) for (int j = i + 1; j < nu; ++j) {
            bool touch = true;
            for (int k = 0; k < 3; ++k) touch &= fabs(unres[i].c[k] - unres[j].c[k]) <= unres[i].h + unres[j].h + 1e-12;
            if (touch) { const int a = cpg_find(parent, i), b = cpg_find(parent, j); if (a != b) parent[b] = a; }
        }
        for (int i = 0; i < nu; ++i) {
            if (cpg_find(parent, i) != i) continue;
            if (info->num_clusters >= MD_TOPO_GTO_MAX_CLUSTERS) { info->clusters_truncated = true; break; }
            md_topo_gto_cluster_t* cl = &info->clusters[info->num_clusters++];
            double lo[3] = { DBL_MAX, DBL_MAX, DBL_MAX }, hi[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX }, hmax = 0.0;
            for (int j = 0; j < nu; ++j) {
                if (cpg_find(parent, j) != i) continue;
                cl->num_boxes++;
                hmax = fmax(hmax, unres[j].h);
                for (int k = 0; k < 3; ++k) { lo[k] = fmin(lo[k], unres[j].c[k] - unres[j].h); hi[k] = fmax(hi[k], unres[j].c[k] + unres[j].h); }
            }
            for (int k = 0; k < 3; ++k) { lo[k] -= hmax; hi[k] += hmax; cl->lo[k] = (float)lo[k]; cl->hi[k] = (float)hi[k]; }
            double gmin = 0.0;
            const double deg = cpg_degree(ctx, &workers[0].sc, lo, hi, 32, &gmin);
            cl->degree = (int32_t)lround(deg);
            cl->degree_residual = (float)fabs(deg - (double)cl->degree);
            cl->min_boundary_grad = (float)gmin;
        }
        md_free(heap, parent, sizeof(int) * nu);
    }
    info->ms_clusters = md_tick_to_milliseconds(md_tick_now() - t_phase);

    // --- output graph (also when cancelled: what was certified so far, flagged incomplete)
    md_allocator_i* galloc = out_graph->alloc ? out_graph->alloc : heap;
    md_topo_extremum_graph_free(out_graph);
    out_graph->alloc = galloc;
    if (ncp > 0) {
        out_graph->num_vertices = (uint32_t)ncp;
        out_graph->vertices = (md_topo_vert_t*)md_alloc(galloc, sizeof(md_topo_vert_t) * ncp);
        out_graph->types = (md_topo_critical_point_type_t*)md_alloc(galloc, sizeof(md_topo_critical_point_type_t) * ncp);
        int count[5] = {0};
        for (int i = 0; i < ncp; ++i) {
            out_graph->vertices[i] = (md_topo_vert_t){ (float)cps[i].x[0], (float)cps[i].x[1], (float)cps[i].x[2], (float)cps[i].rho };
            out_graph->types[i] = (md_topo_critical_point_type_t)cps[i].type;
            count[cps[i].type]++;
        }
        info->poincare_hopf = count[MD_TOPO_MAXIMUM] - count[MD_TOPO_SPLIT_SADDLE] + count[MD_TOPO_JOIN_SADDLE] - count[MD_TOPO_MINIMUM];
        const size_t ne = md_array_size(edges);
        if (ne > 0) {
            out_graph->num_edges = (uint32_t)ne;
            out_graph->edges = (md_topo_edge_t*)md_alloc(galloc, sizeof(md_topo_edge_t) * ne);
            MEMCPY(out_graph->edges, edges, sizeof(md_topo_edge_t) * ne);
        }
    }
    info->complete = (nu == 0) && !info->cancelled;
    md_array_free(edges, heap);
    if (out_info) *out_info = *info;
    return !info->cancelled;
}

static void cpg_run_free(cpg_run_t* R) {
    md_allocator_i* heap = R->heap;
    cpg_ctx_t* ctx = R->ctx;
    const md_gto_basis_t* basis = R->desc->basis;
    const int N = ctx->nao;
    for (int t = 0; t < R->nthreads; ++t) {
        md_array_free(R->workers[t].roots, heap);
        cpg_scratch_free(&R->workers[t].sc, N, ctx->nshell, heap);
    }
    md_free(heap, R->workers, sizeof(cpg_worker_t) * R->nthreads);
    md_array_free(R->roots, heap);
    md_array_free(R->kind, heap);
    md_array_free(R->cps, heap);
    md_array_free(R->cur, heap);
    md_array_free(R->nxt, heap);
    md_array_free(R->unres, heap);
    if (R->tile_need) md_free(heap, R->tile_need, sizeof(uint64_t) * (size_t)R->ntiles * R->tile_words);
    md_free(heap, ctx->ao_nrm, sizeof(double) * N);
    md_free(heap, ctx->ao_ijk, sizeof(int[3]) * N);
    md_free(heap, ctx->coeff, sizeof(double) * basis->num_primitives);
    md_free(heap, ctx->alpha, sizeof(double) * basis->num_primitives);
    md_free(heap, ctx->shell, sizeof(cpg_shell_t) * ctx->nshell);
    if (ctx->rung_r2) md_free(heap, ctx->rung_r2, sizeof(double) * (size_t)ctx->ncls * CPG_RUNGS);
    md_free(heap, ctx, sizeof(cpg_ctx_t));
}

static bool cpg_check_desc(const md_topo_gto_desc_t* desc, const char* fn) {
    if (!desc->basis || !desc->atom_xyz || !desc->density_matrix) {
        MD_LOG_ERROR("%s: basis, atom_xyz and density_matrix are required", fn);
        return false;
    }
    return true;
}

bool md_topo_compute_extremum_graph_gto(md_topo_extremum_graph_t* out_graph, md_topo_gto_info_t* out_info, const md_topo_gto_desc_t* desc) {
    ASSERT(out_graph);
    ASSERT(desc);
    if (!cpg_check_desc(desc, "md_topo_compute_extremum_graph_gto")) return false;
    cpg_run_t R;
    const md_tick_t t_setup = md_tick_now();
    cpg_run_init(&R, desc);
    R.info.ms_setup = md_tick_to_milliseconds(md_tick_now() - t_setup);
    if (R.ok) cpg_run_sweep_cpu(&R);
    const bool res = cpg_run_finish(&R, out_graph, out_info);
    cpg_run_free(&R);
    return res;
}

#undef M1
#undef M2
#undef M3
#undef M4

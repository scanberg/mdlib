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

// The GTO kernel, see the GPU sweep at the end of the file.
static bool topo_gto_gpu_ensure(md_gpu_device_t device);
static void topo_gto_gpu_release(void);

void md_topo_gpu_initialize(md_gpu_device_t device) {
    if (!device) return;
    ensure_kernel(device, &k_bidirectional_manifold,    md_shader_bidirectional_manifold_main_kernel());
    ensure_kernel(device, &k_path_compression,          md_shader_path_compression_main_kernel());
    ensure_kernel(device, &k_critical_points,           md_shader_critical_points_main_kernel());
    ensure_kernel(device, &k_critical_point_compaction, md_shader_critical_point_compaction_main_kernel());
    ensure_kernel(device, &k_vertex_edge_extraction,    md_shader_vertex_edge_extraction_main_kernel());
    topo_gto_gpu_ensure(device);
}

void md_topo_gpu_shutdown(void) {
    topo_gto_gpu_release();
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
// Contributions of AOs screened away from a cube are bounded by a global tail term, so screening does
// not break rigour. The search domain is derived from the same tail bound (rho <= rho_min outside).
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

typedef struct cpg_shell_t {
    double   A[3];
    double   radius;                // screening radius: all derivatives (order <= 2) below tau beyond it
    uint32_t prim_offset;
    uint32_t num_prims;
    uint32_t ao_offset;
    int      l;
    int      ncart;
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
    // D = sum_k fac_l[k] c_k c_k^T, c_k[i] = fac_C[i * fac_r + k] (cpg_factor_density); fac_r = 0: no factors
    int     fac_r;
    double* fac_C;
    double* fac_l;
    md_topo_gto_density_form_t form;
    double tail_rho, tail_g, tail_H;
    double tau;
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
    double* Dl;         // local D block [n*n]
    double* E;          // per-AO remainder vectors [n][16]
    double* Q;          // |D| products [n][8]
    double* Vf;         // factor rows: c_k . V [r][CPG_NV]
    double* Ef;         // factor rows: |c_k| . E [r][16]
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

// Local AO list of the cube (shells whose screening radius reaches it). Returns count.
static int cpg_local_aos(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double c[3], double h, int* out_num_shells) {
    int n = 0, ns = 0;
    for (int si = 0; si < ctx->nshell; ++si) {
        const cpg_shell_t* s = &ctx->shell[si];
        if (cpg_box_dist2(s->A, c, h) < s->radius * s->radius) {
            sc->shells[ns++] = si;
            for (int ci = 0; ci < s->ncart; ++ci) sc->L[n++] = (int)s->ao_offset + ci;
        }
    }
    *out_num_shells = ns;
    return n;
}

// ----------------------------------------------------------------------------------------------- point field

// Whether the cube is evaluated with D's factors (r rows) rather than as a matrix over its n local AOs.
// Cost per cube: factors r n (20 + 15) multiply-adds for the rows, the matrix n^2 (20 + 6) for the D
// products. The factored remainder bounds are looser (sum_k |l_k (c_k . P)| (|c_k| . E) against
// sum_i |(D P)_i| E_i: 5-13% more cubes on the test inputs), so auto wants a clear margin: r <= n / 3.
static bool cpg_use_factors(const cpg_ctx_t* ctx, int n) {
    if (ctx->fac_r == 0 || ctx->form == MD_TOPO_GTO_DENSITY_MATRIX) return false;
    if (ctx->form == MD_TOPO_GTO_DENSITY_FACTORED) return true;
    return 3 * ctx->fac_r <= n;
}

// rho, grad, Hessian at x (deriv 0..2) over all shells within their screening radius.
static void cpg_point(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], int deriv, double* rho, double g[3], double H[3][3]) {
    int ns = 0;
    const int n = cpg_local_aos(ctx, sc, x, 0.0, &ns);
    const int nmi = deriv == 0 ? 1 : deriv == 1 ? 4 : 10;
    int k = 0;
    for (int t = 0; t < ns; ++t) {
        const cpg_shell_t* s = &ctx->shell[sc->shells[t]];
        cpg_shell_eval(ctx, s, x, deriv, sc->V + (size_t)k * CPG_NV, CPG_NV);
        k += s->ncart;
    }
    if (cpg_use_factors(ctx, n)) {
        // the factors (cpg_box_eval): rho = sum_k l_k psi_k^2, psi_k = c_k . phi, r n nmi multiply-adds
        // instead of n^2 ndv; the same function to rounding level (cpg_factor_density)
        const int r = ctx->fac_r;
        MEMSET(sc->Vf, 0, sizeof(double) * r * CPG_NV);
        for (int i = 0; i < n; ++i) {
            const double* ci = ctx->fac_C + (size_t)sc->L[i] * r;
            const double* v = sc->V + (size_t)i * CPG_NV;
            for (int q = 0; q < r; ++q) {
                const double cq = ci[q];
                if (cq == 0.0) continue;
                double* vf = sc->Vf + (size_t)q * CPG_NV;
                for (int m = 0; m < nmi; ++m) vf[m] += cq * v[m];
            }
        }
        double rr = 0.0, gg[3] = {0, 0, 0}, HH[3][3] = {{0}};
        for (int q = 0; q < r; ++q) {
            const double l = ctx->fac_l[q];
            const double* vf = sc->Vf + (size_t)q * CPG_NV;
            rr += l * vf[0] * vf[0];
            if (deriv >= 1) for (int a = 0; a < 3; ++a) gg[a] += 2.0 * l * vf[0] * vf[1 + a];
            if (deriv >= 2) {
                for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) HH[a][b] += 2.0 * l * (vf[M2(a, b)] * vf[0] + vf[1 + a] * vf[1 + b]);
            }
        }
        if (rho) *rho = rr;
        if (g) for (int a = 0; a < 3; ++a) g[a] = gg[a];
        if (H) for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) H[a][b] = a <= b ? HH[a][b] : HH[b][a];
        return;
    }
    const int ndv = deriv >= 2 ? 4 : 1;
    for (int i = 0; i < n; ++i) {
        double acc[4] = {0, 0, 0, 0};
        const double* row = ctx->D + (size_t)sc->L[i] * ctx->nao;
        for (int j = 0; j < n; ++j) {
            const double d = row[sc->L[j]];
            const double* v = sc->V + (size_t)j * CPG_NV;
            for (int m = 0; m < ndv; ++m) acc[m] += d * v[m];
        }
        for (int m = 0; m < ndv; ++m) sc->Dv[(size_t)i * CPG_NDV + m] = acc[m];
    }
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

// ----------------------------------------------------------------------------------------------- cube enclosures

// Returns whether the factored form was used.
static bool cpg_box_eval(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double c[3], double h, cpg_eval_t* ev) {
    int ns = 0;
    const int n = cpg_local_aos(ctx, sc, c, h, &ns);
    memset(ev, 0, sizeof(*ev));
    // centre values (orders 0..3) and per-AO remainder bounds over the cube
    int k = 0;
    for (int t = 0; t < ns; ++t) {
        const cpg_shell_t* s = &ctx->shell[sc->shells[t]];
        cpg_shell_eval(ctx, s, c, 3, sc->V + (size_t)k * CPG_NV, CPG_NV);
        cpg_shell_rem(ctx, s, c, h, sc->E + (size_t)k * 16, 16);
        k += s->ncart;
    }
    // The rows of rho = sum_ab X_ab psi_a psi_b: the local AOs with X = D, or the factors psi_k = c_k . phi
    // with X = diag(l). A factor's Taylor coefficients are c_k . V exactly, and |c_k| . E bounds each of
    // its remainder terms, so everything below holds for either; per row: V (Vr), the remainder vector
    // (Er), Dv = X V and Q = |X| E[0..5].
    const bool fac = cpg_use_factors(ctx, n);
    int nr = n;
    const double* Vr = sc->V;
    const double* Er = sc->E;
    if (fac) {
        const int r = ctx->fac_r;
        nr = r;
        MEMSET(sc->Vf, 0, sizeof(double) * r * CPG_NV);
        MEMSET(sc->Ef, 0, sizeof(double) * r * 16);
        for (int i = 0; i < n; ++i) {
            const double* ci = ctx->fac_C + (size_t)sc->L[i] * r;
            const double* v = sc->V + (size_t)i * CPG_NV;
            const double* e = sc->E + (size_t)i * 16;
            for (int k = 0; k < r; ++k) {
                const double ck = ci[k];
                if (ck == 0.0) continue;
                const double ak = fabs(ck);
                double* vf = sc->Vf + (size_t)k * CPG_NV;
                double* ef = sc->Ef + (size_t)k * 16;
                for (int m = 0; m < CPG_NV; ++m) vf[m] += ck * v[m];
                for (int t = 0; t < 15; ++t) ef[t] += ak * e[t];
            }
        }
        for (int k = 0; k < r; ++k) {
            const double l = ctx->fac_l[k], al = fabs(l);
            const double* vf = sc->Vf + (size_t)k * CPG_NV;
            const double* ef = sc->Ef + (size_t)k * 16;
            for (int m = 0; m < CPG_NDV; ++m) sc->Dv[(size_t)k * CPG_NDV + m] = l * vf[m];
            for (int t = 0; t < 6; ++t) sc->Q[(size_t)k * 8 + t] = al * ef[t];
        }
        Vr = sc->Vf;
        Er = sc->Ef;
    } else {
        // local D block and D-products Dv = D V (all 20 channels)
        double* Dl = sc->Dl;
        for (int i = 0; i < n; ++i) {
            const double* row = ctx->D + (size_t)sc->L[i] * ctx->nao;
            for (int j = 0; j < n; ++j) Dl[(size_t)i * n + j] = row[sc->L[j]];
        }
        for (int i = 0; i < n; ++i) {
            double acc[CPG_NDV] = {0};
            const double* row = Dl + (size_t)i * n;
            for (int j = 0; j < n; ++j) {
                const double d = row[j];
                const double* v = sc->V + (size_t)j * CPG_NV;
                for (int m = 0; m < CPG_NDV; ++m) acc[m] += d * v[m];
            }
            memcpy(sc->Dv + (size_t)i * CPG_NDV, acc, sizeof(acc));
        }
        // |D| products: Q[i*8 + .]: 0 |D|E1, 1 |D|E2, 2 |D|E3, 3..5 |D|Ea2
        for (int i = 0; i < n; ++i) {
            double q[6] = {0};
            const double* row = Dl + (size_t)i * n;
            for (int j = 0; j < n; ++j) {
                const double d = fabs(row[j]);
                const double* e = sc->E + (size_t)j * 16;
                q[0] += d * e[0]; q[1] += d * e[1]; q[2] += d * e[2];
                q[3] += d * e[3]; q[4] += d * e[4]; q[5] += d * e[5];
            }
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
        ev->r[a] = h3 * deg3 + 0.5 * h2 * h2 * deg4 + 2.0 * err_g[a] + ctx->tail_g;
        ev->quadT[a] = 0.5 * h2 * quadT;
    }
    for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
        double lin = 0, quad = 0;
        for (int p = 0; p < 3; ++p) {
            lin += fabs(ev->T[a][b][p]);
            for (int q = 0; q < 3; ++q) quad += fabs(DOT[M3(a, b, p)][M1(q)] + DOT[M2(a, p)][M2(b, q)]);
        }
        const double v = h * lin + 2.0 * h2 * quad + 2.0 * err_H[ab_idx[a][b]] + ctx->tail_H;
        ev->dH[a][b] = ev->dH[b][a] = v;
    }
    ev->rho_up = ev->rho + 2.0 * rho_lin + rho_quad + ctx->tail_rho;
    return fac;
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

static void cpg_unit_grad(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x[3], double sign, double out[3], double* rho, double* gnorm) {
    double g[3];
    cpg_point(ctx, sc, x, 1, rho, g, NULL);
    const double n = sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    *gnorm = n;
    for (int k = 0; k < 3; ++k) out[k] = n > 0 ? sign * g[k] / n : 0.0;
}

// Follows sign * grad rho / |grad rho| from x0 + 1e-3 dir. Returns index of the CP of type 'target'
// it ends in, or -1 (density fell below eps, or path length exhausted). Only the end point is kept,
// and it is found within 'hit' of a CP that attracts the path, so the local error per step (1e-6)
// and the step cap (0.25 Bohr) are loose: on every test input the end points are those of 1e-8 and
// 0.05, at 40% fewer evaluations and a third of the longest path's time.
static int cpg_trace(const cpg_ctx_t* ctx, cpg_scratch_t* sc, const double x0[3], const double dir[3], double sign, const cpg_cp_t* cps, int ncp, int target) {
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
    cpg_queue_t* queue;                 // jobs 1 and 2
    // sweep
    const cpg_box_t* boxes;
    size_t nboxes;
    uint32_t* kind;                     // one slot per cube, written by exactly one worker
    double h_min;
    md_array(cpg_rootrec_t) roots;      // this worker's certified roots
    uint64_t inflated;
    uint64_t factored;                  // cubes evaluated in the factored form
    // separatrices: item t traces path traces[t] % 2 of saddle saddles[traces[t] / 2]
    const cpg_cp_t* cps;
    int ncp;
    const int* saddles;
    const int* traces;
    int* ends;                          // two per saddle
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
    if (cpg_box_eval(ctx, &w->sc, box->c, box->h, &ev)) w->factored++;
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
            w->ends[2 * i + e] = cpg_trace(w->ctx, &w->sc, cp->x, dir, sign, w->cps, w->ncp, target);
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
    sc->Dl     = (double*)md_alloc(heap, sizeof(double) * N * N);
    sc->E      = (double*)md_alloc(heap, sizeof(double) * N * 16);
    sc->Q      = (double*)md_alloc(heap, sizeof(double) * N * 8);
    sc->Vf     = (double*)md_alloc(heap, sizeof(double) * N * CPG_NV);
    sc->Ef     = (double*)md_alloc(heap, sizeof(double) * N * 16);
}

static void cpg_scratch_free(cpg_scratch_t* sc, int N, int nshell, md_allocator_i* heap) {
    md_free(heap, sc->Ef, sizeof(double) * N * 16);
    md_free(heap, sc->Vf, sizeof(double) * N * CPG_NV);
    md_free(heap, sc->Q,  sizeof(double) * N * 8);
    md_free(heap, sc->E,  sizeof(double) * N * 16);
    md_free(heap, sc->Dl, sizeof(double) * N * N);
    md_free(heap, sc->Dv, sizeof(double) * N * CPG_NDV);
    md_free(heap, sc->V,  sizeof(double) * N * CPG_NV);
    md_free(heap, sc->shells, sizeof(int) * (nshell + 1));
    md_free(heap, sc->L,  sizeof(int) * N);
}


// D = sum_k l_k c_k c_k^T with r <= rmax factors, by LDL^T with 1x1 diagonal pivots, the largest
// residual diagonal first (for a positive semidefinite D: pivoted Cholesky, l_k > 0; the factors come
// out more local than eigenvectors, which keeps the factored remainder bounds tighter: 1.2-2.1 against
// 1.7-3.6 for sum |l||c||c|^T / sum |D| on the test inputs). Stops once every residual diagonal is at
// rounding level, then accepts the factors only if the whole residual R = D - C diag(l) C^T is:
//   |R_ij| <= 16 (r + 2) u (|D_ij| + sum_k |l_k c_ik c_jk|),   u = 2^-53,
// i.e. no larger than what the double sums of rho round anyway, so the reference's guarantee (exact
// arithmetic, evaluated in double) is unchanged. A D that is not low rank at that level (a correlated
// density with small occupations, an indefinite one that 1x1 pivots cannot factor) is used as a matrix.
// Returns r, 0 if not factored. Cost: N r^2 / 2 + N^2 r.
static int cpg_factor_density(cpg_ctx_t* ctx, md_allocator_i* heap, int rmax) {
    const int N = ctx->nao;
    const double* D = ctx->D;
    if (N <= 0 || rmax <= 0) return 0;
    double dmax = 0.0;
    for (int i = 0; i < N; ++i) dmax = fmax(dmax, fabs(D[(size_t)i * N + i]));
    if (!(dmax > 0.0)) return 0;
    const double u = 0.5 * DBL_EPSILON;
    double* d   = (double*)md_alloc(heap, sizeof(double) * N);                   // residual diagonal
    double* F   = (double*)md_alloc(heap, sizeof(double) * (size_t)N * rmax);    // c_k[i] at [i * rmax + k]
    double* lam = (double*)md_alloc(heap, sizeof(double) * rmax);
    double* w   = (double*)md_alloc(heap, sizeof(double) * rmax);
    bool*   used = (bool*)md_alloc(heap, sizeof(bool) * N);
    for (int i = 0; i < N; ++i) { d[i] = D[(size_t)i * N + i]; used[i] = false; }
    const double tol = 64.0 * N * u * dmax;
    int r = 0;
    bool ok = true;
    for (;;) {
        int p = -1;
        double best = tol;
        for (int i = 0; i < N; ++i) if (!used[i] && fabs(d[i]) > best) { best = fabs(d[i]); p = i; }
        if (p < 0) break;                                   // the rest is at rounding level
        if (r == rmax) { ok = false; break; }               // rank too high to pay off
        const double* fp = F + (size_t)p * rmax;
        for (int k = 0; k < r; ++k) w[k] = lam[k] * fp[k];
        // the pivot from scratch rather than the running diagonal
        double dp = D[(size_t)p * N + p];
        for (int k = 0; k < r; ++k) dp -= w[k] * fp[k];
        used[p] = true;
        d[p] = 0.0;
        if (!(fabs(dp) > tol)) continue;
        for (int i = 0; i < N; ++i) {
            double* fi = F + (size_t)i * rmax;
            if (i == p) { fi[r] = 1.0; continue; }
            if (used[i]) { fi[r] = 0.0; continue; }
            double v = D[(size_t)i * N + p];
            for (int k = 0; k < r; ++k) v -= w[k] * fi[k];
            fi[r] = v / dp;
            d[i] -= dp * fi[r] * fi[r];
        }
        lam[r] = dp;
        r++;
    }
    // the residual, everywhere (D must also be symmetric to that level)
    if (ok && r > 0) {
        const double fac = 16.0 * (r + 2) * u;
        for (int i = 0; i < N && ok; ++i) {
            const double* fi = F + (size_t)i * rmax;
            for (int k = 0; k < r; ++k) w[k] = lam[k] * fi[k];
            for (int j = i; j < N; ++j) {
                const double* fj = F + (size_t)j * rmax;
                double v = D[(size_t)i * N + j], a = fabs(v);
                for (int k = 0; k < r; ++k) { const double t = w[k] * fj[k]; v -= t; a += fabs(t); }
                if (fabs(v) > fac * a || fabs(D[(size_t)j * N + i] - D[(size_t)i * N + j]) > fac * a) { ok = false; break; }
            }
        }
    }
    if (ok && r > 0) {
        ctx->fac_r = r;
        ctx->fac_l = (double*)md_alloc(heap, sizeof(double) * r);
        ctx->fac_C = (double*)md_alloc(heap, sizeof(double) * (size_t)N * r);
        MEMCPY(ctx->fac_l, lam, sizeof(double) * r);
        for (int i = 0; i < N; ++i) MEMCPY(ctx->fac_C + (size_t)i * r, F + (size_t)i * rmax, sizeof(double) * r);
    }
    md_free(heap, used, sizeof(bool) * N);
    md_free(heap, w, sizeof(double) * rmax);
    md_free(heap, lam, sizeof(double) * rmax);
    md_free(heap, F, sizeof(double) * (size_t)N * rmax);
    md_free(heap, d, sizeof(double) * N);
    return ctx->fac_r;
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
    ctx->form = desc->density_form;
    if (ok && ctx->form != MD_TOPO_GTO_DENSITY_MATRIX) {
        // auto uses factors only where 3 r <= n (cpg_use_factors), so more than N / 3 is never used
        const int rmax = ctx->form == MD_TOPO_GTO_DENSITY_FACTORED ? N : N / 3;
        R->info.density_rank = (uint32_t)cpg_factor_density(ctx, heap, rmax);
    }

    R->workers = (cpg_worker_t*)md_alloc(heap, sizeof(cpg_worker_t) * nthreads);
    MEMSET(R->workers, 0, sizeof(cpg_worker_t) * nthreads);
    for (int t = 0; t < nthreads; ++t) {
        R->workers[t].ctx = ctx;
        R->workers[t].tid = t;
        R->workers[t].nthreads = nthreads;
        R->workers[t].cancel = desc->cancel;
        R->workers[t].h_min = R->h_min;
        cpg_scratch_alloc(&R->workers[t].sc, N, ctx->nshell, heap);
    }
    R->ok = ok;
    if (!ok) return;

    // --- global AO maxima (order 0..2) and |D|-weighted row sums -> screening threshold tau
    double* phimax = (double*)md_alloc(heap, sizeof(double) * N * 3);
    for (int si = 0; si < ctx->nshell; ++si) {
        const cpg_shell_t* s = &ctx->shell[si];
        for (int ci = 0; ci < s->ncart; ++ci) for (int o = 0; o < 3; ++o) phimax[(s->ao_offset + ci) * 3 + o] = cpg_tail_ao(ctx, s, ci, o, 0.0);
    }
    double Rs[3] = {0, 0, 0};
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
        const double d = fabs(ctx->D[(size_t)i * N + j]);
        for (int o = 0; o < 3; ++o) Rs[o] += d * phimax[j * 3 + o];
    }
    md_free(heap, phimax, sizeof(double) * N * 3);
    const double tail_target = 1.0e-12;
    ctx->tau = tail_target / (2.0 * (Rs[0] + 2.0 * Rs[1] + Rs[2]) + 1e-300);
    ctx->tail_rho = 2.0 * ctx->tau * Rs[0];
    ctx->tail_g   = 2.0 * ctx->tau * (Rs[0] + Rs[1]);
    ctx->tail_H   = 2.0 * ctx->tau * (Rs[0] + 2.0 * Rs[1] + Rs[2]);
    for (int si = 0; si < ctx->nshell; ++si) {
        cpg_shell_t* s = &ctx->shell[si];
        double lo = 0.0, hi = 1.0;
        while (cpg_tail_shell(ctx, s, 2, hi) > ctx->tau && hi < 1e3) hi *= 2.0;
        for (int it = 0; it < 50; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (cpg_tail_shell(ctx, s, 2, mid) > ctx->tau) lo = mid; else hi = mid;
        }
        s->radius = hi;
    }

    // --- search domain: outside AABB(atoms) +- pad every point is >= pad from every centre; rho <= T^T |D| T
    double amin[3] = { DBL_MAX, DBL_MAX, DBL_MAX }, amax[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
    for (int si = 0; si < ctx->nshell; ++si) for (int k = 0; k < 3; ++k) { amin[k] = fmin(amin[k], ctx->shell[si].A[k]); amax[k] = fmax(amax[k], ctx->shell[si].A[k]); }
    double* T0 = (double*)md_alloc(heap, sizeof(double) * N);
    double plo = 0.0, phi = 1.0;
    for (;;) {
        const double pad = phi;
        for (int si = 0; si < ctx->nshell; ++si) { const cpg_shell_t* s = &ctx->shell[si]; for (int ci = 0; ci < s->ncart; ++ci) T0[s->ao_offset + ci] = cpg_tail_ao(ctx, s, ci, 0, pad); }
        double b = 0.0;
        for (int i = 0; i < N; ++i) { if (T0[i] == 0.0) continue; for (int j = 0; j < N; ++j) b += T0[i] * fabs(ctx->D[(size_t)i * N + j]) * T0[j]; }
        if (b < ctx->eps || phi > 1e3) break;
        plo = phi; phi *= 2.0;
    }
    for (int it = 0; it < 40; ++it) {
        const double mid = 0.5 * (plo + phi);
        for (int si = 0; si < ctx->nshell; ++si) { const cpg_shell_t* s = &ctx->shell[si]; for (int ci = 0; ci < s->ncart; ++ci) T0[s->ao_offset + ci] = cpg_tail_ao(ctx, s, ci, 0, mid); }
        double b = 0.0;
        for (int i = 0; i < N; ++i) { if (T0[i] == 0.0) continue; for (int j = 0; j < N; ++j) b += T0[i] * fabs(ctx->D[(size_t)i * N + j]) * T0[j]; }
        if (b < ctx->eps) phi = mid; else plo = mid;
    }
    md_free(heap, T0, sizeof(double) * N);
    const double pad = phi;
    R->info.domain_pad = pad;

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
            R->info.num_factored_evals += workers[t].factored;
            workers[t].factored = 0;
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
            cpg_queue_t queue = { .next = 0, .count = (int)md_array_size(traces) };
            md_mutex_init(&queue.mutex);       // in place: a CRITICAL_SECTION or pthread mutex must not be copied
            for (int t = 0; t < nthreads; ++t) {
                workers[t].job = 1;
                workers[t].queue = &queue;
                workers[t].cps = cps;
                workers[t].ncp = ncp;
                workers[t].saddles = saddles;
                workers[t].traces = traces;
                workers[t].ends = ends;
            }
            if (queue.count > 0) cpg_run_workers(workers, MIN(nthreads, queue.count));
            md_mutex_destroy(&queue.mutex);
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
    if (ctx->fac_r) {
        md_free(heap, ctx->fac_C, sizeof(double) * (size_t)N * ctx->fac_r);
        md_free(heap, ctx->fac_l, sizeof(double) * ctx->fac_r);
    }
    md_free(heap, ctx->ao_nrm, sizeof(double) * N);
    md_free(heap, ctx->ao_ijk, sizeof(int[3]) * N);
    md_free(heap, ctx->coeff, sizeof(double) * basis->num_primitives);
    md_free(heap, ctx->alpha, sizeof(double) * basis->num_primitives);
    md_free(heap, ctx->shell, sizeof(cpg_shell_t) * ctx->nshell);
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

#if MD_ENABLE_GPU
// ----------------------------------------------------------------------------------------------- GPU sweep
// The octree levels run on the GPU (src/shaders/topo/topo_gto_cube.slang, fp32 with rigorous rounding
// margins), in batches of the 8 children of each cube that split; the CPU keeps the parts that need
// double precision or are inherently few:
//   * Newton polish and classification of the root cubes the GPU certified (tens to hundreds);
//   * every cube fp32 cannot decide (its gradient lies inside its own rounding noise, or it reached
//     h_min, or the kernels cannot form its coordinates exactly in float): that cube and its subtree
//     go through the CPU level loop above, in double;
//   * the local shell lists (screened per batch against its parent batch's list), separatrices,
//     clusters and the graph (cpg_run_finish).
// The batch lists live on the host and are rebuilt in batch and child order after every level, so the
// result does not depend on scheduling or on how a level is cut into chunks.

typedef struct cpg_gpu_shell_t {
    float    A[3];
    float    r2;            // squared screening radius, rounded up (unused by the kernels)
    uint32_t prim_offset;
    uint32_t num_prims;
    uint32_t ao_offset;
    uint32_t l_ncart;       // l | ncart << 8
} cpg_gpu_shell_t;

// Mirrors Batch in topo_gto_cube.slang.
typedef struct cpg_gpu_batch_t {
    float    P_h[4];        // parent centre, child half-width
    uint32_t info[4];       // child mask, local AOs, first row, first local shell
    uint32_t info2[4];      // number of local shells, first factor row (CPG_GPU_FACT_NONE: matrix form)
} cpg_gpu_batch_t;
#define CPG_GPU_FACT_NONE 0xFFFFFFFFu

// Mirrors RootArgs in topo_gto_cube.slang.
typedef struct cpg_gpu_args_t {
    uint32_t num_batches;
    uint32_t num_shells;
    uint32_t num_ao;
    uint32_t num_row_tiles;
    float    eps;
    float    h_min;
    float    tail_rho;
    float    tail_g;
    float    tail_H;
    uint32_t flags;
    uint32_t fac_r;
    uint32_t num_fbatches;
    md_gpu_addr_t batches;
    md_gpu_addr_t batch_shells;
    md_gpu_addr_t row_ao;
    md_gpu_addr_t row_tiles;
    md_gpu_addr_t shells;
    md_gpu_addr_t alpha;
    md_gpu_addr_t coeff;
    md_gpu_addr_t ao_nrm;
    md_gpu_addr_t ao_ijk;
    md_gpu_addr_t D;
    md_gpu_addr_t phi;
    md_gpu_addr_t dv;
    md_gpu_addr_t ep;
    md_gpu_addr_t outcome;
    md_gpu_addr_t debug;
    md_gpu_addr_t fac_C;
    md_gpu_addr_t fac_l;
    md_gpu_addr_t fphi;
    md_gpu_addr_t fbatches;
} cpg_gpu_args_t;

// Layout constants shared with topo_gto_cube.slang.
#define CPG_GPU_PHI_W          456             // floats per row of phi: 57 channels x 8 children
#define CPG_GPU_DV_W           272             // floats per row of dv: 34 x 8
#define CPG_GPU_EP_W           260             // floats per cube of ep
#define CPG_GPU_TAB_CAP        4096            // floats of 1D tables per chunk of shells (ao_main)
#define CPG_GPU_DEC_WG         64              // decide_main group size

#define CPG_GPU_SCRATCH_BUDGET (192u << 20)    // bytes of phi + dv per chunk
#define CPG_GPU_MAX_BATCHES    16384u          // per chunk (grid dimensions stay below 65535)
#define CPG_GPU_MAX_TILES      65535u
#define CPG_GPU_DISPATCH_MS    20.0            // target duration of one chunk (four launches)
#define CPG_GPU_OUT_KIND(o)    ((o) & 0xFu)
#define CPG_GPU_OUT_INFLATED   16u
#define CPG_GPU_WHY_CHILDREN   11u             // discarded: split, but no child survived the parent's expansion

// The four stages; GEMM and EPI have a second kernel for the batches in the factored form.
enum { CPG_K_AO, CPG_K_GEMM, CPG_K_EPI, CPG_K_DECIDE, CPG_K_FGEMM16, CPG_K_FGEMM32, CPG_K_FGEMM64, CPG_K_EPIF, CPG_K_COUNT };
static md_gpu_kernel_t k_topo_gto[CPG_K_COUNT] = {0};
static md_gpu_device_t k_topo_gto_device = NULL;

#define CPG_GPU_GEMM_ROWS 64             // gemm_main: rows per group tile

// Factored GEMM tilings (topo_gto_cube.slang: FGEMM): factor rows per tile; any r works with any of
// them (ceil(r / rows) tiles per batch), the default takes the smallest that holds r in one tile.
static const struct { uint32_t rows; int kernel; const char* name; } cpg_fgemm_variants[] = {
    { 16, CPG_K_FGEMM16, "f0: 16 factor rows per group, 32x2 threads, 8x3 outputs per thread" },
    { 32, CPG_K_FGEMM32, "f1: 32 factor rows per group, 32x4 threads, 8x3 outputs per thread" },
    { 64, CPG_K_FGEMM64, "f2: 64 factor rows per group, 32x8 threads, 8x3 outputs per thread" },
};
#define CPG_FGEMM_NV ((uint32_t)(sizeof(cpg_fgemm_variants) / sizeof(cpg_fgemm_variants[0])))
uint32_t md_topo_gto_gpu_fgemm_variant_count(void) { return CPG_FGEMM_NV; }
const char* md_topo_gto_gpu_fgemm_variant_name(uint32_t v) { return v < CPG_FGEMM_NV ? cpg_fgemm_variants[v].name : ""; }
uint32_t md_topo_gto_gpu_fgemm_variant_auto(uint32_t rank) { return rank <= 16 ? 0 : (rank <= 32 ? 1 : 2); }

static void topo_gto_gpu_release(void) {
    for (int i = 0; i < CPG_K_COUNT; ++i) {
        if (k_topo_gto[i]) md_gpu_kernel_destroy(k_topo_gto[i]);
        k_topo_gto[i] = NULL;
    }
    k_topo_gto_device = NULL;
}

// Created by md_topo_gpu_initialize, or on first use. Not synchronised: run one GTO search at a time
// per process, or initialise first.
static bool topo_gto_gpu_ensure(md_gpu_device_t device) {
    if (k_topo_gto_device == device && k_topo_gto[CPG_K_DECIDE]) return true;
    topo_gto_gpu_release();
    const md_gpu_kernel_desc_t kd[CPG_K_COUNT] = {
        md_shader_topo_gto_cube_ao_main_kernel(),
        md_shader_topo_gto_cube_gemm_main_kernel(),
        md_shader_topo_gto_cube_epi_main_kernel(),
        md_shader_topo_gto_cube_decide_main_kernel(),
        md_shader_topo_gto_cube_fgemm16_kernel(),
        md_shader_topo_gto_cube_fgemm32_kernel(),
        md_shader_topo_gto_cube_fgemm64_kernel(),
        md_shader_topo_gto_cube_epi_fact_kernel(),
    };
    for (int i = 0; i < CPG_K_COUNT; ++i) {
        k_topo_gto[i] = md_gpu_kernel_create(device, &kd[i]);
        if (!k_topo_gto[i]) {
            MD_LOG_ERROR("md_topo: failed to create kernel '%s': %s", kd[i].label, md_gpu_last_error());
            topo_gto_gpu_release();
            return false;
        }
    }
    k_topo_gto_device = device;
    return true;
}

static float cpg_up_f(double x) {
    float f = (float)x;
    if ((double)f < x) f = nextafterf(f, FLT_MAX);
    return f;
}

// Every buffer the kernels read is padded: a compiler may issue a loop's next load ahead of the exit
// test (llvmpipe does), which must not run off the end of an allocation.
#define CPG_GPU_PAD 4096

static bool cpg_gpu_upload(md_gpu_stream_t s, md_gpu_addr_t* out, const void* src, size_t size) {
    *out = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, size + CPG_GPU_PAD).gpu;
    return *out && (size == 0 || md_gpu_upload(s, *out, src, size));
}

// Grow-only device buffer of 'need' elements of 'elem' bytes.
static bool cpg_gpu_reserve(md_gpu_stream_t s, md_gpu_addr_t* buf, size_t* cap, size_t need, size_t elem) {
    if (*buf && need <= *cap) return true;
    if (*buf) md_gpu_free(s, *buf);
    size_t c = *cap + *cap / 2;
    if (c < need) c = need;
    if (c < 64) c = 64;
    *buf = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, c * elem + CPG_GPU_PAD).gpu;
    *cap = *buf ? c : 0;
    return *buf != 0;
}

// A batch: the 8 children (c = P +- hc per axis, child s on the + side of x if s & 1, y if s & 2,
// z if s & 4) of one cube, those in 'mask' to be evaluated, and its local shells: a range of the
// pool, screened when the batch is made from its parent batch's local shells. A child batch's union
// of inflated cubes, P +- 2.5 hc, lies inside its parent batch's, so screening against the parent's
// list loses nothing.
typedef struct cpg_hbatch_t {
    double   P[3];
    double   hc;
    uint32_t mask;
    uint32_t level;
    uint32_t lsh_cnt;
    uint32_t nao;                       // local AOs; UINT32_MAX: the kernels cannot form it, the CPU takes it
    size_t   lsh_off;                   // first local shell, in the pool's numbering (pool index + pool base)
} cpg_hbatch_t;

static inline bool cpg_is_float(double x) { return (double)(float)x == x; }

// The kernels form the children, their intervals and the inflated intervals in float: P +- j hc / 2,
// j = -5..5, and 1.5 hc. All must be exact.
static bool cpg_gpu_batch_exact(const cpg_hbatch_t* b) {
    if (!cpg_is_float(b->hc) || !cpg_is_float(1.5 * b->hc)) return false;
    for (int k = 0; k < 3; ++k) {
        for (int j = -5; j <= 5; ++j) if (!cpg_is_float(b->P[k] + 0.5 * j * b->hc)) return false;
    }
    return true;
}

// Screens batch b (P, hc, mask and level set) against the candidates pool[cand_off - pool_base ...],
// cand_cnt of them, and appends its local shells to the pool.
static void cpg_gpu_batch_screen(const cpg_ctx_t* ctx, cpg_hbatch_t* b, md_array(uint32_t)* pool, size_t pool_base,
                                 size_t cand_off, uint32_t cand_cnt, md_allocator_i* heap) {
    md_array_ensure(*pool, md_array_size(*pool) + cand_cnt, heap);   // the candidates must not move below
    b->lsh_off = pool_base + md_array_size(*pool);
    b->lsh_cnt = 0;
    if (!cpg_gpu_batch_exact(b)) {
        b->nao = UINT32_MAX;
        return;
    }
    const uint32_t* cand = cand_cnt ? *pool + (cand_off - pool_base) : NULL;
    const double hs = 2.5 * b->hc;
    uint32_t n = 0;
    for (uint32_t c = 0; c < cand_cnt; ++c) {
        const uint32_t s = cand[c];
        const cpg_shell_t* cs = &ctx->shell[s];
        double d2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            const double t = fabs(cs->A[k] - b->P[k]) - hs;
            if (t > 0.0) d2 += t * t;
        }
        // a lower bound of the distance, conservatively: any shell that may reach is kept
        if (d2 * (1.0 - 1e-12) <= cs->radius * cs->radius * (1.0 + 1e-6)) {
            md_array_push_no_grow(*pool, s);
            b->lsh_cnt++;
            n += (uint32_t)cs->ncart;
        }
    }
    b->nao = n;
}

// A chunk: the batches q[q0..q1), 'count' of them for the GPU (the others go to the CPU).
typedef struct cpg_gpu_chunk_t {
    size_t        q0, q1;
    size_t        count;
    bool          full;                 // as many GPU batches as the chunk size allowed
    bool          gpu_idle;             // no other GPU work in flight when it was launched
    int           slot;                 // its readback buffer
    md_gpu_sync_t sync;
    md_tick_t     t_launch;
    double        ms_profiled;          // with profiling: its time, every kernel waited for
} cpg_gpu_chunk_t;

static inline cpg_box_t cpg_child_box(const cpg_hbatch_t* b, int s) {
    cpg_box_t box = { { b->P[0] + ((s & 1) ? b->hc : -b->hc), b->P[1] + ((s & 2) ? b->hc : -b->hc), b->P[2] + ((s & 4) ? b->hc : -b->hc) }, b->hc };
    return box;
}

typedef struct cpg_rootkey_t {
    int64_t  g[3];
    uint32_t bit;
    uint32_t idx;
} cpg_rootkey_t;

static int cpg_rootkey_cmp(const void* pa, const void* pb) {
    const cpg_rootkey_t* a = (const cpg_rootkey_t*)pa;
    const cpg_rootkey_t* b = (const cpg_rootkey_t*)pb;
    for (int k = 2; k >= 0; --k) {
        if (a->g[k] != b->g[k]) return a->g[k] < b->g[k] ? -1 : 1;
    }
    return a->bit < b->bit ? -1 : (a->bit > b->bit ? 1 : 0);
}

static inline int64_t cpg_floor_div2(int64_t i) { return i >= 0 ? i / 2 : -((-i + 1) / 2); }

// Sweeps R->cur on the GPU. On return R->cur holds the cubes left to the CPU (escalated ones, or, if
// the GPU failed part way, everything not yet decided). Returns false if the GPU could not be used.
static bool cpg_run_sweep_gpu(cpg_run_t* R, md_gpu_stream_t stream) {
    md_gpu_device_t dev = md_gpu_stream_device(stream);
    if (!dev || !topo_gto_gpu_ensure(dev)) return false;
    const md_tick_t t_start = md_tick_now();
    const cpg_ctx_t* ctx = R->ctx;
    md_allocator_i* heap = R->heap;
    const int N = ctx->nao;
    const uint32_t gemm_rows = CPG_GPU_GEMM_ROWS;
    const int NS = ctx->nshell;
    int NP = 0;
    for (int s = 0; s < NS; ++s) {
        const cpg_shell_t* cs = &ctx->shell[s];
        NP = MAX(NP, (int)(cs->prim_offset + cs->num_prims));
        if (cs->l > 4 || (size_t)cs->num_prims * 6 * (20 * (cs->l + 1) + 2) > CPG_GPU_TAB_CAP) {
            MD_LOG_INFO("md_topo_compute_extremum_graph_gto_gpu: a shell (l = %d, %d primitives) exceeds the GPU kernel's tables, using the CPU", cs->l, (int)cs->num_prims);
            return false;
        }
    }

    // --- tables in float. The basis and the atom positions are float to begin with (exact); the
    //     tails are rounded up; D and the normalisation are rounded, which the kernels' error model
    //     accounts for.
    cpg_gpu_shell_t* sh = (cpg_gpu_shell_t*)md_alloc(heap, sizeof(cpg_gpu_shell_t) * NS);
    for (int s = 0; s < NS; ++s) {
        const cpg_shell_t* cs = &ctx->shell[s];
        for (int k = 0; k < 3; ++k) sh[s].A[k] = (float)cs->A[k];
        sh[s].r2 = cpg_up_f(cs->radius * cs->radius * (1.0 + 1e-6));
        sh[s].prim_offset = cs->prim_offset;
        sh[s].num_prims = cs->num_prims;
        sh[s].ao_offset = cs->ao_offset;
        sh[s].l_ncart = (uint32_t)cs->l | ((uint32_t)cs->ncart << 8);
    }
    float* fa = (float*)md_alloc(heap, sizeof(float) * MAX(NP, 1));
    float* fc = (float*)md_alloc(heap, sizeof(float) * MAX(NP, 1));
    for (int p = 0; p < NP; ++p) { fa[p] = (float)ctx->alpha[p]; fc[p] = (float)ctx->coeff[p]; }
    float* fn = (float*)md_alloc(heap, sizeof(float) * N);
    uint32_t* ijk = (uint32_t*)md_alloc(heap, sizeof(uint32_t) * N);
    for (int i = 0; i < N; ++i) {
        fn[i] = (float)ctx->ao_nrm[i];
        ijk[i] = (uint32_t)ctx->ao_ijk[i][0] | ((uint32_t)ctx->ao_ijk[i][1] << 8) | ((uint32_t)ctx->ao_ijk[i][2] << 16);
    }
    float* fD = (float*)md_alloc(heap, sizeof(float) * (size_t)N * N);
    for (size_t i = 0; i < (size_t)N * N; ++i) fD[i] = (float)ctx->D[i];

    cpg_gpu_args_t a = {0};
    a.num_shells = (uint32_t)NS;
    a.num_ao = (uint32_t)N;
    a.eps = (float)ctx->eps;            // a threshold, not a rigorous quantity
    a.h_min = (float)R->h_min;
    a.tail_rho = cpg_up_f(ctx->tail_rho);
    a.tail_g = cpg_up_f(ctx->tail_g);
    a.tail_H = cpg_up_f(ctx->tail_H);

    // the factors of D in float, if any batch may use them (the same rule as the CPU's, per batch)
    const uint32_t R_f = ctx->fac_r && ctx->form != MD_TOPO_GTO_DENSITY_MATRIX ? (uint32_t)ctx->fac_r : 0;
    a.fac_r = R_f;
    uint32_t fv = R->desc->gpu_fgemm_variant ? R->desc->gpu_fgemm_variant - 1 : md_topo_gto_gpu_fgemm_variant_auto(R_f);
    if (fv >= CPG_FGEMM_NV) fv = md_topo_gto_gpu_fgemm_variant_auto(R_f);
    const uint32_t fgemm_rows = cpg_fgemm_variants[fv].rows;
    const int k_fgemm = cpg_fgemm_variants[fv].kernel;
    const uint32_t fgemm_tpb = R_f ? (R_f + fgemm_rows - 1) / fgemm_rows : 0;     // row tiles per batch
    float* fC = R_f ? (float*)md_alloc(heap, sizeof(float) * (size_t)N * R_f) : NULL;
    float* fL = R_f ? (float*)md_alloc(heap, sizeof(float) * R_f) : NULL;
    for (size_t i = 0; i < (size_t)N * R_f; ++i) fC[i] = (float)ctx->fac_C[i];
    for (uint32_t k = 0; k < R_f; ++k) fL[k] = (float)ctx->fac_l[k];

    bool ok = cpg_gpu_upload(stream, &a.shells, sh, sizeof(cpg_gpu_shell_t) * NS)
           && cpg_gpu_upload(stream, &a.alpha, fa, sizeof(float) * NP)
           && cpg_gpu_upload(stream, &a.coeff, fc, sizeof(float) * NP)
           && cpg_gpu_upload(stream, &a.ao_nrm, fn, sizeof(float) * N)
           && cpg_gpu_upload(stream, &a.ao_ijk, ijk, sizeof(uint32_t) * N)
           && cpg_gpu_upload(stream, &a.D, fD, sizeof(float) * (size_t)N * N)
           && (!R_f || (cpg_gpu_upload(stream, &a.fac_C, fC, sizeof(float) * (size_t)N * R_f)
                        && cpg_gpu_upload(stream, &a.fac_l, fL, sizeof(float) * R_f)));
    if (R_f) {
        md_free(heap, fL, sizeof(float) * R_f);
        md_free(heap, fC, sizeof(float) * (size_t)N * R_f);
    }
    md_free(heap, fD, sizeof(float) * (size_t)N * N);
    md_free(heap, ijk, sizeof(uint32_t) * N);
    md_free(heap, fn, sizeof(float) * N);
    md_free(heap, fc, sizeof(float) * MAX(NP, 1));
    md_free(heap, fa, sizeof(float) * MAX(NP, 1));
    md_free(heap, sh, sizeof(cpg_gpu_shell_t) * NS);
    if (!ok) MD_LOG_ERROR("md_topo_compute_extremum_graph_gto_gpu: GPU allocation failed: %s", md_gpu_last_error());

    md_array(cpg_box_t)       esc = 0;        // cubes for the CPU, in the order met
    md_array(cpg_box_t)       proots = 0;     // cubes with a certified root, polished after the sweep
    md_array(cpg_hbatch_t)    q = 0;          // batches, first in first out: level by level, each level in
                                              // the order the previous level's outcomes made it
    md_array(uint32_t)        pool = 0;       // local shells of the batches in q, in the same order
    size_t                    pool_base = 0;  // pool numbering of pool[0]
    md_array(cpg_gpu_batch_t) hb = 0;
    md_array(uint32_t)        hsh = 0;        // the chunk's local shells
    md_array(uint32_t)        htile = 0;      // pairs (batch, first row)
    md_array(uint32_t)        hfb = 0;        // the chunk's factored batches

    // --- root lattice -> batches of 8, children in the lattice's 2x2x2 groups; candidates: every shell
    for (int s = 0; s < NS; ++s) md_array_push(pool, (uint32_t)s, heap);
    {
        const size_t nr = md_array_size(R->cur);
        cpg_rootkey_t* keys = (cpg_rootkey_t*)md_alloc(heap, sizeof(cpg_rootkey_t) * MAX(nr, 1));
        size_t nk = 0;
        for (size_t i = 0; i < nr; ++i) {
            const cpg_box_t* b = &R->cur[i];
            cpg_rootkey_t key = { {0, 0, 0}, 0, (uint32_t)i };
            bool fits = true;
            for (int k = 0; k < 3; ++k) {
                const int64_t I = (int64_t)floor(b->c[k] / (2.0 * b->h));
                const int64_t G = cpg_floor_div2(I);
                const uint32_t bit = (uint32_t)(I - 2 * G);
                const double P = ((double)G + 0.5) * 4.0 * b->h;
                fits &= (P + (bit ? b->h : -b->h)) == b->c[k];
                key.g[k] = G;
                key.bit |= bit << k;
            }
            if (fits) keys[nk++] = key;
            else md_array_push(esc, *b, heap);
        }
        if (nk > 1) qsort(keys, nk, sizeof(cpg_rootkey_t), cpg_rootkey_cmp);
        for (size_t i = 0; i < nk; ) {
            const cpg_box_t* b0 = &R->cur[keys[i].idx];
            cpg_hbatch_t hbt = { { 0, 0, 0 }, b0->h, 0, 0, 0, 0, 0 };
            for (int k = 0; k < 3; ++k) hbt.P[k] = ((double)keys[i].g[k] + 0.5) * 4.0 * b0->h;
            size_t j = i;
            while (j < nk && keys[j].g[0] == keys[i].g[0] && keys[j].g[1] == keys[i].g[1] && keys[j].g[2] == keys[i].g[2] && R->cur[keys[j].idx].h == b0->h) {
                hbt.mask |= 1u << keys[j].bit;
                j++;
            }
            cpg_gpu_batch_screen(ctx, &hbt, &pool, pool_base, 0, (uint32_t)NS, heap);
            md_array_push(q, hbt, heap);
            i = j;
        }
        md_free(heap, keys, sizeof(cpg_rootkey_t) * MAX(nr, 1));
        md_array_shrink(R->cur, 0);
    }

    // --- device scratch, grown on demand. One set: the stream orders a chunk's uploads and kernels after
    //     the previous chunk's. Only the readback is double, as the host reads one chunk's outcomes
    //     while the next chunk runs.
    md_gpu_addr_t d_batches = 0, d_pool = 0, d_row_ao = 0, d_tiles = 0, d_phi = 0, d_dv = 0, d_ep = 0, d_out = 0, d_fphi = 0, d_fb = 0;
    size_t c_batches = 0, c_pool = 0, c_row_ao = 0, c_tiles = 0, c_phi = 0, c_dv = 0, c_ep = 0, c_out = 0, c_fphi = 0, c_fb = 0;
    md_gpu_mem_t h_out[2] = { {0}, {0} };
    size_t c_hout[2] = { 0, 0 };
    const size_t float_cap = CPG_GPU_SCRATCH_BUDGET / sizeof(float);   // phi + dv + fphi per chunk
    const bool profile = R->desc->profile_gpu_kernels;   // one chunk at a time, every kernel waited for and timed
    // Two chunks in flight: while the GPU runs one, the host reads the other's outcomes, makes and screens
    // the batches they split into, and builds and launches the next chunk. Chunks follow q, so they run
    // on across level boundaries. The order everything is decided in is q's, whatever the chunks.
    const int depth = profile ? 1 : 2;
    cpg_gpu_chunk_t fl[2];
    int nfl = 0, slot = 0;
    size_t chunk = 16;                    // GPU batches per chunk, adapted to CPG_GPU_DISPATCH_MS
    size_t qd = 0, qh = 0;                // q[qd..qh) launched, not yet decided; q[qh..] not launched
    uint32_t num_levels = 0;
    bool have_prev_done = false;          // when the previous GPU chunk completed, if seen to the tick
    md_tick_t t_prev_done = 0;

    while (ok) {
        if (cpg_run_cancelled(R)) break;

        // --- launch while there is room in flight and batches to launch
        while (ok && nfl < depth && qh < md_array_size(q)) {
            const md_tick_t t_form = md_tick_now();
            md_array_shrink(hb, 0);
            md_array_shrink(hsh, 0);
            md_array_shrink(htile, 0);
            md_array_shrink(hfb, 0);
            const size_t nq = md_array_size(q);
            size_t b = qh, rows = 0, tiles = 0, count = 0, frows = 0;
            for (; b < nq && count < chunk && count < CPG_GPU_MAX_BATCHES; ++b) {
                const cpg_hbatch_t* bt = &q[b];
                if (bt->nao == UINT32_MAX) continue;
                const uint32_t n = bt->nao;
                const bool fac = R_f && cpg_use_factors(ctx, (int)n);
                const size_t nt = fac ? 0 : (n + gemm_rows - 1) / gemm_rows;
                const size_t nfb = md_array_size(hfb) + (fac ? 1 : 0);
                if (count > 0 && ((rows + n) * (CPG_GPU_PHI_W + CPG_GPU_DV_W) + (frows + (fac ? R_f : 0)) * CPG_GPU_PHI_W > float_cap ||
                                  tiles + nt > CPG_GPU_MAX_TILES || nfb * fgemm_tpb * 5 > CPG_GPU_MAX_TILES)) break;
                cpg_gpu_batch_t g = { { (float)bt->P[0], (float)bt->P[1], (float)bt->P[2], (float)bt->hc },
                                      { bt->mask, n, (uint32_t)rows, (uint32_t)md_array_size(hsh) },
                                      { bt->lsh_cnt, fac ? (uint32_t)frows : CPG_GPU_FACT_NONE, 0, 0 } };
                md_array_push(hb, g, heap);
                if (fac) {
                    md_array_push(hfb, (uint32_t)count, heap);
                    frows += R_f;
                    R->info.gpu_gemm_flop += (double)fgemm_tpb * 2.0 * fgemm_rows * 480.0 * (double)((n + 15) / 16 * 16);
                }
                const uint32_t* ls = pool + (bt->lsh_off - pool_base);
                md_array_push_array(hsh, ls, bt->lsh_cnt, heap);   // rows (AO per row) are laid out by ao_main
                for (uint32_t r0 = 0; !fac && r0 < n; r0 += gemm_rows) {
                    md_array_push(htile, (uint32_t)count, heap);
                    md_array_push(htile, r0, heap);
                }
                rows += n;
                tiles += nt;
                count++;
                R->info.gpu_gemm_flop += (double)nt * 2.0 * gemm_rows * 288.0 * (double)((n + 15) / 16 * 16);
            }
            const size_t nfb = md_array_size(hfb);
            cpg_gpu_chunk_t ch = { qh, b, count, count == chunk, true, slot, md_gpu_sync_none(), 0, 0.0 };
            for (int i = 0; i < nfl; ++i) ch.gpu_idle &= fl[i].count == 0;
            if (count > 0) {
                R->info.num_gpu_batches += count;
                R->info.num_gpu_rows += rows;
                R->info.num_gpu_factored_batches += nfb;
                const size_t nsh = md_array_size(hsh);
                ok = cpg_gpu_reserve(stream, &d_batches, &c_batches, count, sizeof(cpg_gpu_batch_t))
                  && cpg_gpu_reserve(stream, &d_pool, &c_pool, MAX(nsh, 1), sizeof(uint32_t))
                  && cpg_gpu_reserve(stream, &d_row_ao, &c_row_ao, MAX(rows, 1), sizeof(uint32_t))
                  && cpg_gpu_reserve(stream, &d_tiles, &c_tiles, MAX(tiles, 1), 2 * sizeof(uint32_t))
                  && cpg_gpu_reserve(stream, &d_phi, &c_phi, MAX(rows, 1), CPG_GPU_PHI_W * sizeof(float))
                  && cpg_gpu_reserve(stream, &d_dv, &c_dv, MAX(rows, 1), CPG_GPU_DV_W * sizeof(float))
                  && cpg_gpu_reserve(stream, &d_ep, &c_ep, 8 * count, CPG_GPU_EP_W * sizeof(float))
                  && cpg_gpu_reserve(stream, &d_out, &c_out, 8 * count, sizeof(uint32_t))
                  && cpg_gpu_reserve(stream, &d_fphi, &c_fphi, MAX(frows, 1), CPG_GPU_PHI_W * sizeof(float))
                  && cpg_gpu_reserve(stream, &d_fb, &c_fb, MAX(nfb, 1), sizeof(uint32_t));
                if (ok && (!h_out[slot].gpu || c_hout[slot] < 8 * count)) {
                    if (h_out[slot].gpu) md_gpu_free(stream, h_out[slot].gpu);
                    c_hout[slot] = MAX(8 * count, c_hout[slot] + c_hout[slot] / 2);
                    h_out[slot] = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, c_hout[slot] * sizeof(uint32_t));
                    ok = h_out[slot].gpu && h_out[slot].cpu;
                }
                const md_tick_t tu = md_tick_now();
                ok = ok && md_gpu_upload(stream, d_batches, hb, count * sizeof(cpg_gpu_batch_t))
                        && (nsh == 0 || md_gpu_upload(stream, d_pool, hsh, nsh * sizeof(uint32_t)))
                        && (tiles == 0 || md_gpu_upload(stream, d_tiles, htile, tiles * 2 * sizeof(uint32_t)))
                        && (nfb == 0 || md_gpu_upload(stream, d_fb, hfb, nfb * sizeof(uint32_t)));
                if (!ok) break;
                a.num_batches = (uint32_t)count;
                a.num_row_tiles = (uint32_t)tiles;
                a.batches = d_batches;
                a.batch_shells = d_pool;
                a.row_ao = d_row_ao;
                a.row_tiles = d_tiles;
                a.phi = d_phi;
                a.dv = d_dv;
                a.ep = d_ep;
                a.outcome = d_out;
                a.fphi = d_fphi;
                a.fbatches = d_fb;
                a.num_fbatches = (uint32_t)nfb;
                const size_t nmb = count - nfb;           // batches in the matrix form
                // in order; a zero grid is skipped
                const struct { md_gpu_kernel_t k; md_gpu_grid_t grid; double* ms; } launch[] = {
                    { k_topo_gto[CPG_K_AO],     md_gpu_grid((uint32_t)count, 1, 1),                                      &R->info.ms_gpu_ao },
                    { k_topo_gto[CPG_K_GEMM],   md_gpu_grid((uint32_t)tiles, 1, 1),                                      &R->info.ms_gpu_gemm },
                    { k_topo_gto[k_fgemm],      md_gpu_grid((uint32_t)(nfb * fgemm_tpb * 5), 1, 1),                      &R->info.ms_gpu_gemm },
                    { k_topo_gto[CPG_K_EPI],    md_gpu_grid(nmb ? (uint32_t)((count + 1) / 2) : 0, 1, 1),                &R->info.ms_gpu_epilogue },  // two batches per group
                    { k_topo_gto[CPG_K_EPIF],   md_gpu_grid(nfb ? (uint32_t)((count + 1) / 2) : 0, 1, 1),                &R->info.ms_gpu_epilogue },
                    { k_topo_gto[CPG_K_DECIDE], md_gpu_grid((uint32_t)((8 * count + CPG_GPU_DEC_WG - 1) / CPG_GPU_DEC_WG), 1, 1), &R->info.ms_gpu_decide },
                };
                for (size_t k = 0; k < sizeof(launch) / sizeof(launch[0]) && ok; ++k) {
                    if (launch[k].grid.x == 0) continue;
                    if (profile) md_gpu_stream_sync(stream);      // uploads and the previous kernel are not this kernel's time
                    const md_tick_t tk = md_tick_now();
                    ok = md_gpu_launch(stream, launch[k].k, launch[k].grid, &a, sizeof(a));
                    if (profile && ok) {
                        md_gpu_stream_sync(stream);
                        *launch[k].ms += md_tick_to_milliseconds(md_tick_now() - tk);
                    }
                    R->info.num_gpu_dispatches++;
                }
                ok = ok && md_gpu_copy(stream, h_out[slot].gpu, d_out, 8 * count * sizeof(uint32_t));
                if (!ok) break;
                ch.sync = md_gpu_stream_record(stream);
                if (profile) {
                    md_gpu_sync_wait(ch.sync);
                    ch.ms_profiled = md_tick_to_milliseconds(md_tick_now() - tu);
                    R->info.ms_sweep_gpu_wait += ch.ms_profiled;
                }
            }
            ch.t_launch = md_tick_now();
            R->info.ms_sweep_host_launch += md_tick_to_milliseconds(ch.t_launch - t_form) - ch.ms_profiled;
            fl[nfl++] = ch;
            slot ^= 1;
            qh = b;
        }
        if (!ok || nfl == 0) break;

        // --- the oldest chunk in flight
        const cpg_gpu_chunk_t c = fl[0];
        if (c.count > 0) {
            const bool was_done = md_gpu_sync_is_complete(c.sync);
            const md_tick_t tw = md_tick_now();
            if (!was_done) md_gpu_sync_wait(c.sync);
            const md_tick_t t1 = md_tick_now();
            if (!profile) R->info.ms_sweep_gpu_wait += md_tick_to_milliseconds(t1 - tw);
            // The chunk size adapts to a target time per chunk on the GPU (TDR on Windows trips at 2 s,
            // and a long chunk starves rendering on the same GPU); results do not depend on it. A chunk's
            // time runs from its launch, or from the previous chunk's completion if it was queued behind
            // it. Had it completed before the host came to wait, that is only an upper bound, enough to
            // grow the chunk but not to shrink it.
            if (c.full) {
                double est;
                bool exact;
                if (profile) {
                    est = c.ms_profiled;
                    exact = true;
                } else {
                    md_tick_t from = c.t_launch;
                    exact = !was_done;
                    if (!c.gpu_idle) {
                        if (have_prev_done && t_prev_done > from) from = t_prev_done;
                        else if (!have_prev_done) exact = false;
                    }
                    est = md_tick_to_milliseconds(t1 - from);
                }
                const double f = est > 0.0 ? CPG_GPU_DISPATCH_MS / est : 2.0;
                if (exact || f > 1.0) chunk = (size_t)CLAMP((double)chunk * CLAMP(f, 0.25, 2.0), 1.0, (double)CPG_GPU_MAX_BATCHES);
            }
            have_prev_done = !was_done && !profile;
            t_prev_done = t1;
        }

        // --- its outcomes in batch and child order: certified roots (polished later), splits, the rest to the CPU
        const uint32_t* out = c.count > 0 ? (const uint32_t*)h_out[c.slot].cpu : NULL;
        const md_tick_t t_out = md_tick_now();
        size_t gi = 0;
        for (size_t bi = c.q0; bi < c.q1; ++bi) {
            const cpg_hbatch_t bt = q[bi];      // a copy: the children made below grow q
            num_levels = MAX(num_levels, bt.level + 1);
            if (bt.nao == UINT32_MAX) {
                for (int s = 0; s < 8; ++s) if (bt.mask & (1u << s)) md_array_push(esc, cpg_child_box(&bt, s), heap);
                continue;
            }
            ASSERT(gi < c.count);
            for (int s = 0; s < 8; ++s) {
                if (!(bt.mask & (1u << s))) continue;
                const cpg_box_t box = cpg_child_box(&bt, s);
                const uint32_t o = out[gi * 8 + s];
                R->info.num_box_evals++;
                R->info.num_gpu_box_evals++;
                if (o & CPG_GPU_OUT_INFLATED) R->info.num_inflated_evals++;
                switch (CPG_GPU_OUT_KIND(o)) {
                case 0:
                    if (((o >> 8) & 0xFFu) == CPG_GPU_WHY_CHILDREN) R->info.num_children_skipped += 8;
                    break;
                case 1:
                    md_array_push(proots, box, heap);
                    break;
                case 2: {
                    // the children the cube's own expansion could not exclude (topo_gto_cube.slang: child_excluded)
                    const uint32_t cm = (o >> 16) & 0xFFu;
                    for (int k = 0; k < 8; ++k) R->info.num_children_skipped += (cm >> k) & 1u ? 0 : 1;
                    if (cm) {
                        cpg_hbatch_t child = { { box.c[0], box.c[1], box.c[2] }, 0.5 * box.h, cm, bt.level + 1, 0, 0, 0 };
                        cpg_gpu_batch_screen(ctx, &child, &pool, pool_base, bt.lsh_off, bt.lsh_cnt, heap);
                        md_array_push(q, child, heap);
                    }
                    break;
                }
                default:
                    md_array_push(esc, box, heap);
                    break;
                }
            }
            gi++;
        }
        R->info.ms_sweep_host_outcomes += md_tick_to_milliseconds(md_tick_now() - t_out);
        qd = c.q1;
        fl[0] = fl[1];
        nfl--;

        // --- drop the decided batches and their shells once they are most of q
        if (qd >= 4096 && 2 * qd >= md_array_size(q)) {
            const size_t nq = md_array_size(q) - qd;
            memmove(q, q + qd, nq * sizeof(cpg_hbatch_t));
            md_array_shrink(q, nq);
            for (int i = 0; i < nfl; ++i) { fl[i].q0 -= qd; fl[i].q1 -= qd; }
            qh -= qd;
            qd = 0;
            const size_t base = nq > 0 ? q[0].lsh_off : pool_base + md_array_size(pool);
            const size_t drop = base - pool_base;
            memmove(pool, pool + drop, (md_array_size(pool) - drop) * sizeof(uint32_t));
            md_array_shrink(pool, md_array_size(pool) - drop);
            pool_base = base;
        }
    }
    if (!ok) MD_LOG_ERROR("md_topo_compute_extremum_graph_gto_gpu: GPU sweep failed, the CPU takes over: %s", md_gpu_last_error());
    if (nfl > 0) md_gpu_stream_sync(stream);     // left with work in flight: let it end before its buffers go
    R->info.num_levels += num_levels;

    // --- the certified roots: Newton in double on every worker, then accepted in the order they were met
    const size_t npr = md_array_size(proots);
    if (npr > 0) {
        const md_tick_t tp = md_tick_now();
        cpg_polish_t* pres = (cpg_polish_t*)md_alloc(heap, sizeof(cpg_polish_t) * npr);
        for (size_t i = 0; i < npr; ++i) pres[i].r = -1;   // left so if cancelled: the CPU takes the cube
        cpg_queue_t queue = { .next = 0, .count = (int)npr };
        md_mutex_init(&queue.mutex);
        for (int t = 0; t < R->nthreads; ++t) {
            R->workers[t].job = 2;
            R->workers[t].queue = &queue;
            R->workers[t].pboxes = proots;
            R->workers[t].pres = pres;
        }
        cpg_run_workers(R->workers, (int)MIN((size_t)R->nthreads, npr));
        md_mutex_destroy(&queue.mutex);
        for (size_t i = 0; i < npr; ++i) {
            if (pres[i].r == 1) cpg_run_accept(R, &pres[i].cp);
            else if (pres[i].r < 0) md_array_push(esc, proots[i], heap);   // Newton failed despite the certificate
        }
        md_free(heap, pres, sizeof(cpg_polish_t) * npr);
        R->info.ms_sweep_polish += md_tick_to_milliseconds(md_tick_now() - tp);
    }

    // whatever is left (escalated, or not reached if the GPU failed or was cancelled) goes to the CPU loop
    R->info.num_escalated_boxes = (uint32_t)md_array_size(esc);
    for (size_t bi = qd; bi < md_array_size(q); ++bi) {
        for (int s = 0; s < 8; ++s) if (q[bi].mask & (1u << s)) md_array_push(esc, cpg_child_box(&q[bi], s), heap);
    }
    md_array_free(R->cur, heap);
    R->cur = esc;

    md_gpu_addr_t bufs[] = { d_batches, d_pool, d_row_ao, d_tiles, d_phi, d_dv, d_ep, d_out, d_fphi, d_fb, h_out[0].gpu, h_out[1].gpu,
                             a.shells, a.alpha, a.coeff, a.ao_nrm, a.ao_ijk, a.D, a.fac_C, a.fac_l };
    for (size_t i = 0; i < sizeof(bufs) / sizeof(bufs[0]); ++i) if (bufs[i]) md_gpu_free(stream, bufs[i]);
    md_array_free(proots, heap);
    md_array_free(q, heap);
    md_array_free(pool, heap);
    md_array_free(hb, heap);
    md_array_free(hsh, heap);
    md_array_free(htile, heap);
    md_array_free(hfb, heap);
    R->info.ms_sweep += md_tick_to_milliseconds(md_tick_now() - t_start);
    return ok;
}

bool md_topo_compute_extremum_graph_gto_gpu(md_topo_extremum_graph_t* out_graph, md_topo_gto_info_t* out_info, const md_topo_gto_desc_t* desc, md_gpu_stream_t stream) {
    ASSERT(out_graph);
    ASSERT(desc);
    if (!cpg_check_desc(desc, "md_topo_compute_extremum_graph_gto_gpu")) return false;
    cpg_run_t R;
    const md_tick_t t_setup = md_tick_now();
    cpg_run_init(&R, desc);
    R.info.ms_setup = md_tick_to_milliseconds(md_tick_now() - t_setup);
    if (R.ok) {
        R.info.used_gpu = stream && cpg_run_sweep_gpu(&R, stream);
        cpg_run_sweep_cpu(&R);          // escalated cubes, or everything if the GPU was unavailable
    }
    const bool res = cpg_run_finish(&R, out_graph, out_info);
    cpg_run_free(&R);
    return res;
}
#endif

#undef M1
#undef M2
#undef M3
#undef M4

#include <md_topo.h>

#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_grid.h>

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
static md_gpu_compute_pipeline_t pip_bidirectional_manifold = NULL;
static md_gpu_compute_pipeline_t pip_path_compression = NULL;
static md_gpu_compute_pipeline_t pip_critical_points = NULL;
static md_gpu_compute_pipeline_t pip_critical_point_compaction = NULL;
static md_gpu_compute_pipeline_t pip_vertex_edge_extraction = NULL;

static md_gpu_compute_pipeline_t ensure_pipeline(md_gpu_device_t device, md_gpu_compute_pipeline_t* cached, const uint32_t* blob_start, size_t blob_size, const char* name, uint32_t wg_x, uint32_t wg_y, uint32_t wg_z) {
    if (cached_device != device) {
        // Device changed, invalidate all cached pipelines
        if (pip_bidirectional_manifold)    { md_gpu_destroy_compute_pipeline(pip_bidirectional_manifold);    pip_bidirectional_manifold = NULL; }
        if (pip_path_compression)          { md_gpu_destroy_compute_pipeline(pip_path_compression);          pip_path_compression = NULL; }
        if (pip_critical_points)           { md_gpu_destroy_compute_pipeline(pip_critical_points);           pip_critical_points = NULL; }
        if (pip_critical_point_compaction) { md_gpu_destroy_compute_pipeline(pip_critical_point_compaction); pip_critical_point_compaction = NULL; }
        if (pip_vertex_edge_extraction)    { md_gpu_destroy_compute_pipeline(pip_vertex_edge_extraction);    pip_vertex_edge_extraction = NULL; }
        cached_device = device;
    }
    if (*cached == NULL) {
        md_gpu_compute_pipeline_desc_t desc = {
            .shader = { .data = blob_start, .size = blob_size },
            .threadgroup_size = { wg_x, wg_y, wg_z },
        };
        *cached = md_gpu_create_compute_pipeline(device, &desc);
        if (*cached == NULL) {
            MD_LOG_ERROR("Failed to create compute pipeline: %s", name);
        }
    }
    return *cached;
}

bool md_topo_compute_extremum_graph_gpu(md_topo_extremum_graph_t* out_graph, md_gpu_device_t device, md_gpu_image_t volume, const md_grid_t* grid, float scalar_threshold) {
    if (!out_graph || !device || !volume || !grid) {
        MD_LOG_ERROR("Invalid input: out_graph=%p, device=%p, volume=%p, grid=%p", (void*)out_graph, (void*)device, (void*)volume, (void*)grid);
        return false;
    }

    md_allocator_i* alloc = out_graph->alloc ? out_graph->alloc : md_get_heap_allocator();

    const uint32_t num_points = (uint32_t)(grid->dim[0] * grid->dim[1] * grid->dim[2]);

    // UBO / push constants (same layout as shader UBO struct)
    struct {
        float index_to_world[4][4];
        uint32_t dims[3];
        float scalar_threshold;
    } ubo_data;
    index_to_world_matrix(ubo_data.index_to_world, grid);
    ubo_data.dims[0] = grid->dim[0];
    ubo_data.dims[1] = grid->dim[1];
    ubo_data.dims[2] = grid->dim[2];
    ubo_data.scalar_threshold = scalar_threshold;

    // Create pipelines (lazy, cached)
    md_gpu_compute_pipeline_t p_manifold   = ensure_pipeline(device, &pip_bidirectional_manifold,    topo_bidirectional_manifold_start,    topo_bidirectional_manifold_size(),    "bidirectional_manifold",    8, 8, 8);
    md_gpu_compute_pipeline_t p_compress   = ensure_pipeline(device, &pip_path_compression,          topo_path_compression_start,          topo_path_compression_size(),          "path_compression",          8, 8, 8);
    md_gpu_compute_pipeline_t p_critical   = ensure_pipeline(device, &pip_critical_points,           topo_critical_points_start,           topo_critical_points_size(),           "critical_points",           8, 8, 8);
    md_gpu_compute_pipeline_t p_compact    = ensure_pipeline(device, &pip_critical_point_compaction, topo_critical_point_compaction_start, topo_critical_point_compaction_size(), "critical_point_compaction", 8, 8, 8);
    md_gpu_compute_pipeline_t p_extract    = ensure_pipeline(device, &pip_vertex_edge_extraction,    topo_vertex_edge_extraction_start,    topo_vertex_edge_extraction_size(),    "vertex_edge_extraction",    64, 1, 1);

    if (!p_manifold || !p_compress || !p_critical || !p_compact || !p_extract) {
        return false;
    }

    md_gpu_queue_t queue = md_gpu_acquire_compute_queue(device);

    // Create GPU buffers
    md_gpu_buffer_t ascending_buf  = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_points * sizeof(uint32_t) });
    md_gpu_buffer_t descending_buf = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_points * sizeof(uint32_t) });
    md_gpu_buffer_t types_buf      = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_points * sizeof(int32_t) });
    md_gpu_buffer_t counts_buf     = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = 4 * sizeof(uint32_t) });
    md_gpu_buffer_t changed_buf    = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = sizeof(uint32_t) });

    bool success = false;
    md_gpu_buffer_t indices_buf    = NULL;
    md_gpu_buffer_t counters_buf   = NULL;
    md_gpu_buffer_t type_counts_buf = NULL;
    md_gpu_buffer_t vert_buf       = NULL;
    md_gpu_buffer_t edge_buf       = NULL;
    md_gpu_buffer_t edge_count_buf = NULL;

    // Staging buffers for readback
    md_gpu_buffer_t staging_counts     = NULL;
    md_gpu_buffer_t staging_edge_count = NULL;
    md_gpu_buffer_t staging_verts      = NULL;
    md_gpu_buffer_t staging_edges      = NULL;

    // === Steps 1-2: Bidirectional manifold + Path compression ===
    {
        md_gpu_command_buffer_t cmd = md_gpu_acquire_command_buffer(device);

        // Step 1: Bidirectional manifold
        md_gpu_cmd_bind_compute_pipeline(cmd, p_manifold);
        md_gpu_cmd_push_constants(cmd, &ubo_data, sizeof(ubo_data));
        md_gpu_cmd_bind_image(cmd, 0, volume);
        md_gpu_cmd_bind_buffer(cmd, 0, ascending_buf);
        md_gpu_cmd_bind_buffer(cmd, 1, descending_buf);
        md_gpu_cmd_dispatch(cmd, grid->dim[0], grid->dim[1], grid->dim[2]);
        md_gpu_cmd_barrier(cmd);

        // Step 2: Path compression (iterative)
        uint32_t num_iterations = 0;
        uint32_t max_dim = grid->dim[0];
        if (grid->dim[1] > max_dim) max_dim = grid->dim[1];
        if (grid->dim[2] > max_dim) max_dim = grid->dim[2];
        while (max_dim > (1U << num_iterations)) {
            num_iterations++;
        }
        num_iterations += 2;

        md_gpu_cmd_bind_compute_pipeline(cmd, p_compress);
        md_gpu_cmd_push_constants(cmd, &ubo_data, sizeof(ubo_data));
        md_gpu_cmd_bind_buffer(cmd, 0, ascending_buf);
        md_gpu_cmd_bind_buffer(cmd, 1, descending_buf);
        md_gpu_cmd_bind_buffer(cmd, 2, changed_buf);
        for (uint32_t i = 0; i < num_iterations; i++) {
            md_gpu_cmd_dispatch(cmd, grid->dim[0], grid->dim[1], grid->dim[2]);
            md_gpu_cmd_barrier(cmd);
        }

        // Step 3: Critical points
        md_gpu_cmd_fill_buffer(cmd, counts_buf, 0, 4 * sizeof(uint32_t), 0);
        md_gpu_cmd_barrier(cmd);

        md_gpu_cmd_bind_compute_pipeline(cmd, p_critical);
        md_gpu_cmd_push_constants(cmd, &ubo_data, sizeof(ubo_data));
        md_gpu_cmd_bind_image(cmd, 0, volume);
        md_gpu_cmd_bind_buffer(cmd, 0, ascending_buf);
        md_gpu_cmd_bind_buffer(cmd, 1, descending_buf);
        md_gpu_cmd_bind_buffer(cmd, 2, types_buf);
        md_gpu_cmd_bind_buffer(cmd, 3, counts_buf);
        md_gpu_cmd_dispatch(cmd, grid->dim[0], grid->dim[1], grid->dim[2]);
        md_gpu_cmd_barrier(cmd);

        // Copy counts to staging for readback
        staging_counts = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = 4 * sizeof(uint32_t), .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        md_gpu_cmd_copy_buffer(cmd, counts_buf, staging_counts, 4 * sizeof(uint32_t), 0, 0);

        md_gpu_fence_t fence = md_gpu_queue_submit(queue, cmd);
        md_gpu_fence_wait(fence);
        md_gpu_destroy_fence(fence);
    }

    // Read back counts
    uint32_t counts[4] = {0};
    {
        void* ptr = md_gpu_map_buffer(staging_counts);
        MEMCPY(counts, ptr, sizeof(counts));
        md_gpu_unmap_buffer(staging_counts);
    }

    uint32_t num_maxima        = counts[0];
    uint32_t num_split_saddles = counts[1];
    uint32_t num_minima        = counts[2];
    uint32_t num_join_saddles  = counts[3];
    uint32_t num_vertices      = num_maxima + num_split_saddles + num_minima + num_join_saddles;
    uint32_t num_edges         = 8 * (num_split_saddles + num_join_saddles);

    MD_LOG_DEBUG("Topology: %u maxima, %u split saddles, %u minima, %u join saddles (total: %u)",
                num_maxima, num_split_saddles, num_minima, num_join_saddles, num_vertices);

    if (num_vertices == 0) {
        MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
        out_graph->alloc = alloc;
        success = true;
        goto cleanup;
    }

    // === Steps 4-5: Compaction + Extraction ===
    indices_buf    = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_vertices * sizeof(uint32_t) });
    counters_buf   = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = 4 * sizeof(uint32_t) });
    type_counts_buf = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = 5 * sizeof(uint32_t) });
    vert_buf       = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_vertices * 4 * sizeof(float) });
    edge_buf       = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_edges * sizeof(md_topo_edge_t) });
    edge_count_buf = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = sizeof(uint32_t) });

    {
        md_gpu_command_buffer_t cmd = md_gpu_acquire_command_buffer(device);

        // Upload counter initial values: write offsets for each type
        // We need to fill counters_buf with {0, num_maxima, num_maxima+num_split_saddles, num_maxima+num_split_saddles+num_minima}
        // Using fill_buffer for each word
        md_gpu_cmd_fill_buffer(cmd, counters_buf, 0 * sizeof(uint32_t), sizeof(uint32_t), 0);
        md_gpu_cmd_fill_buffer(cmd, counters_buf, 1 * sizeof(uint32_t), sizeof(uint32_t), num_maxima);
        md_gpu_cmd_fill_buffer(cmd, counters_buf, 2 * sizeof(uint32_t), sizeof(uint32_t), num_maxima + num_split_saddles);
        md_gpu_cmd_fill_buffer(cmd, counters_buf, 3 * sizeof(uint32_t), sizeof(uint32_t), num_maxima + num_split_saddles + num_minima);
        md_gpu_cmd_barrier(cmd);

        // Step 4: Critical point compaction
        md_gpu_cmd_bind_compute_pipeline(cmd, p_compact);
        md_gpu_cmd_push_constants(cmd, &ubo_data, sizeof(ubo_data));
        md_gpu_cmd_bind_buffer(cmd, 0, types_buf);
        md_gpu_cmd_bind_buffer(cmd, 1, indices_buf);
        md_gpu_cmd_bind_buffer(cmd, 2, counters_buf);
        md_gpu_cmd_dispatch(cmd, grid->dim[0], grid->dim[1], grid->dim[2]);
        md_gpu_cmd_barrier(cmd);

        // Upload type_counts
        md_gpu_cmd_fill_buffer(cmd, type_counts_buf, 0 * sizeof(uint32_t), sizeof(uint32_t), num_maxima);
        md_gpu_cmd_fill_buffer(cmd, type_counts_buf, 1 * sizeof(uint32_t), sizeof(uint32_t), num_split_saddles);
        md_gpu_cmd_fill_buffer(cmd, type_counts_buf, 2 * sizeof(uint32_t), sizeof(uint32_t), num_minima);
        md_gpu_cmd_fill_buffer(cmd, type_counts_buf, 3 * sizeof(uint32_t), sizeof(uint32_t), num_join_saddles);
        md_gpu_cmd_fill_buffer(cmd, type_counts_buf, 4 * sizeof(uint32_t), sizeof(uint32_t), num_vertices);
        md_gpu_cmd_fill_buffer(cmd, edge_count_buf, 0, sizeof(uint32_t), 0);
        md_gpu_cmd_barrier(cmd);

        // Step 5: Vertex + edge extraction
        md_gpu_cmd_bind_compute_pipeline(cmd, p_extract);
        md_gpu_cmd_push_constants(cmd, &ubo_data, sizeof(ubo_data));
        md_gpu_cmd_bind_buffer(cmd, 0, indices_buf);
        md_gpu_cmd_bind_buffer(cmd, 1, type_counts_buf);
        md_gpu_cmd_bind_buffer(cmd, 2, ascending_buf);
        md_gpu_cmd_bind_buffer(cmd, 3, descending_buf);
        md_gpu_cmd_bind_buffer(cmd, 4, vert_buf);
        md_gpu_cmd_bind_buffer(cmd, 5, edge_buf);
        md_gpu_cmd_bind_buffer(cmd, 6, edge_count_buf);
        md_gpu_cmd_bind_buffer(cmd, 7, types_buf);
        md_gpu_cmd_bind_image(cmd, 0, volume);
        md_gpu_cmd_dispatch(cmd, num_vertices, 1, 1);
        md_gpu_cmd_barrier(cmd);

        // Copy results to staging
        staging_edge_count = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = sizeof(uint32_t), .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        staging_verts      = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_vertices * 4 * sizeof(float), .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        staging_edges      = md_gpu_create_buffer(device, &(md_gpu_buffer_desc_t){ .size = num_edges * sizeof(md_topo_edge_t), .flags = MD_GPU_BUFFER_CPU_VISIBLE });
        md_gpu_cmd_copy_buffer(cmd, edge_count_buf, staging_edge_count, sizeof(uint32_t), 0, 0);
        md_gpu_cmd_copy_buffer(cmd, vert_buf, staging_verts, num_vertices * 4 * sizeof(float), 0, 0);
        md_gpu_cmd_copy_buffer(cmd, edge_buf, staging_edges, num_edges * sizeof(md_topo_edge_t), 0, 0);

        md_gpu_fence_t fence = md_gpu_queue_submit(queue, cmd);
        md_gpu_fence_wait(fence);
        md_gpu_destroy_fence(fence);
    }

    // Read back edge count
    {
        void* ptr = md_gpu_map_buffer(staging_edge_count);
        MEMCPY(&num_edges, ptr, sizeof(uint32_t));
        md_gpu_unmap_buffer(staging_edge_count);
    }

    // Read back vertices (shader writes float4 per vertex: xyz = position, w = value)
    md_topo_vert_t* vertices = md_alloc(alloc, num_vertices * sizeof(md_topo_vert_t));
    {
        float* ptr = (float*)md_gpu_map_buffer(staging_verts);
        for (uint32_t i = 0; i < num_vertices; i++) {
            vertices[i].x     = ptr[i * 4 + 0];
            vertices[i].y     = ptr[i * 4 + 1];
            vertices[i].z     = ptr[i * 4 + 2];
            vertices[i].value = ptr[i * 4 + 3];
        }
        md_gpu_unmap_buffer(staging_verts);
    }

    // Read back edges
    md_topo_edge_t* edges = NULL;
    if (num_edges > 0) {
        edges = md_alloc(alloc, num_edges * sizeof(md_topo_edge_t));
        void* ptr = md_gpu_map_buffer(staging_edges);
        MEMCPY(edges, ptr, num_edges * sizeof(md_topo_edge_t));
        md_gpu_unmap_buffer(staging_edges);
    }

    // Fill output structure
    MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
    out_graph->num_vertices     = num_vertices;
    out_graph->vertices         = vertices;
    out_graph->num_maxima       = num_maxima;
    out_graph->num_split_saddles = num_split_saddles;
    out_graph->num_minima       = num_minima;
    out_graph->num_join_saddles = num_join_saddles;
    out_graph->num_edges        = num_edges;
    out_graph->edges            = edges;
    out_graph->alloc            = alloc;
    success = true;

cleanup:
    if (staging_counts)     md_gpu_destroy_buffer(staging_counts);
    if (staging_edge_count) md_gpu_destroy_buffer(staging_edge_count);
    if (staging_verts)      md_gpu_destroy_buffer(staging_verts);
    if (staging_edges)      md_gpu_destroy_buffer(staging_edges);
    if (ascending_buf)      md_gpu_destroy_buffer(ascending_buf);
    if (descending_buf)     md_gpu_destroy_buffer(descending_buf);
    if (types_buf)          md_gpu_destroy_buffer(types_buf);
    if (counts_buf)         md_gpu_destroy_buffer(counts_buf);
    if (changed_buf)        md_gpu_destroy_buffer(changed_buf);
    if (indices_buf)        md_gpu_destroy_buffer(indices_buf);
    if (counters_buf)       md_gpu_destroy_buffer(counters_buf);
    if (type_counts_buf)    md_gpu_destroy_buffer(type_counts_buf);
    if (vert_buf)           md_gpu_destroy_buffer(vert_buf);
    if (edge_buf)           md_gpu_destroy_buffer(edge_buf);
    if (edge_count_buf)     md_gpu_destroy_buffer(edge_count_buf);
    return success;
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

void md_topo_extremum_graph_free(md_topo_extremum_graph_t* graph) {
    if (graph && graph->alloc) {
        md_allocator_i* alloc = graph->alloc;
        if (graph->vertices) {
            md_free(alloc, graph->vertices, graph->num_vertices * sizeof(md_topo_vert_t));
        }
        if (graph->edges) {
            md_free(alloc, graph->edges, graph->num_edges * sizeof(md_topo_edge_t));
        }
        MEMSET(graph, 0, sizeof(md_topo_extremum_graph_t));
        graph->alloc = alloc;
    }
}

void md_topo_extremum_graph_copy(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* src_graph) {
    // Copy structure, assume that the allocator is set within out_graph
    md_allocator_i* alloc = out_graph->alloc ? out_graph->alloc : md_get_heap_allocator();
    MEMSET(out_graph, 0, sizeof(md_topo_extremum_graph_t));
    out_graph->alloc = alloc;
    out_graph->num_vertices = src_graph->num_vertices;
    out_graph->num_maxima = src_graph->num_maxima;
    out_graph->num_split_saddles = src_graph->num_split_saddles;
    out_graph->num_minima = src_graph->num_minima;
    out_graph->num_join_saddles = src_graph->num_join_saddles;
    out_graph->num_edges = src_graph->num_edges;

    out_graph->vertices = md_alloc(alloc, src_graph->num_vertices * sizeof(md_topo_vert_t));
    out_graph->edges    = md_alloc(alloc, src_graph->num_edges    * sizeof(md_topo_edge_t));


    MEMCPY(out_graph->vertices, src_graph->vertices, src_graph->num_vertices * sizeof(md_topo_vert_t));
    MEMCPY(out_graph->edges, src_graph->edges, src_graph->num_edges * sizeof(md_topo_edge_t));
}

#include <md_topo.h>

#include <core/md_platform.h>
#include <core/md_log.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_grid.h>

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

#if !MD_PLATFORM_OSX

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

bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const md_grid_t* grid) {
    if (!out_graph || vol_tex == 0 || !grid) {
        MD_LOG_ERROR("Invalid input: out_graph=%p, vol_tex=%u, grid=%p", (void*)out_graph, vol_tex, (void*)grid);
        return false;
    }

    GLuint query;
    glGenQueries(1, &query);

    // Start timing
    glBeginQuery(GL_TIME_ELAPSED, query);
    
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
    ubo_data.scalar_threshold = 0.0f; // Set to >0.0 to filter noise if desired
    
    GLuint ubo_buf = create_buffer(sizeof(ubo_data), &ubo_data, GL_STATIC_DRAW);
    if (!ubo_buf) {
        return false;
    }
    
    // Auto-compute compression iterations based on grid dimensions
    uint32_t max_dim = MAX(grid->dim[0], MAX(grid->dim[1], grid->dim[2]));
    uint32_t num_compression_iterations = 0;
    uint32_t temp = max_dim;
    while (temp > 1) {
        temp >>= 1;
    }
    num_compression_iterations += 1; // Add safety margin
    
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
    
    ascending_buf  = create_buffer(num_points * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    descending_buf = create_buffer(num_points * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    
    glUseProgram(manifold_prog);
    glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_READ_ONLY, GL_R32F);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, descending_buf);
    
    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    
    // === Step 2: Path compression (iteratively) ===
    GLuint compression_prog = get_path_compression_program();
    if (!compression_prog) goto cleanup;
    
    glUseProgram(compression_prog);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ascending_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, descending_buf);
    
    for (uint32_t i = 0; i < num_compression_iterations; ++i) {
        glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }
    
    // === Step 3: Identify critical points ===
    GLuint critical_prog = get_critical_points_program();
    if (!critical_prog) goto cleanup;
    
    types_buf = create_buffer(num_points * sizeof(int), NULL, GL_DYNAMIC_COPY);
    
    // Counts buffer: maximaCount, splitSaddleCount, minimaCount, joinSaddleCount, edgeCount
    uint32_t counts_init[5] = {0, 0, 0, 0, 0};
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

    {
		// Ensure all writes are done before reading back
        GLsync fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glClientWaitSync(fence, GL_SYNC_FLUSH_COMMANDS_BIT, GL_TIMEOUT_IGNORED);
        glDeleteSync(fence);
    }
    
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
    
    MD_LOG_INFO("Topology: %u maxima, %u split saddles, %u minima, %u join saddles (total: %u)",
                num_maxima, num_split_saddles, num_minima, num_join_saddles, num_vertices);
    
    // === Step 4: Compact critical point indices into ordered array ===
    GLuint compaction_prog = get_critical_point_compaction_program();
    if (!compaction_prog) goto cleanup;
    
    // Create separate buffers for each type
    GLuint maxima_buf = create_buffer(num_maxima * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    GLuint split_saddles_buf = create_buffer(num_split_saddles * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    GLuint minima_buf = create_buffer(num_minima * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    GLuint join_saddles_buf = create_buffer(num_join_saddles * sizeof(uint32_t), NULL, GL_DYNAMIC_COPY);
    
    uint32_t counters_init[4] = {0};
    counter_buf = create_buffer(sizeof(counters_init), counters_init, GL_DYNAMIC_COPY);
    
    glUseProgram(compaction_prog);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, types_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, maxima_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, split_saddles_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, minima_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, join_saddles_buf);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, counter_buf);

    glDispatchCompute(num_workgroups[0], num_workgroups[1], num_workgroups[2]);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    
    // === Step 5: Extract vertices and edges using GPU shader ===
    GLuint extraction_prog = get_vertex_edge_extraction_program();
    if (!extraction_prog) goto cleanup;
    
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

    // Copy indices from separate buffers into unified buffer in correct order
    if (num_maxima > 0 && indices_buf) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, maxima_buf);
        glBindBuffer(GL_COPY_WRITE_BUFFER, indices_buf);
        glCopyBufferSubData(GL_SHADER_STORAGE_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, num_maxima * sizeof(uint32_t));
    }
    uint32_t offset = num_maxima;
    if (num_split_saddles > 0 && indices_buf) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, split_saddles_buf);
        glBindBuffer(GL_COPY_WRITE_BUFFER, indices_buf);
        glCopyBufferSubData(GL_SHADER_STORAGE_BUFFER, GL_COPY_WRITE_BUFFER, 0, offset * sizeof(uint32_t), num_split_saddles * sizeof(uint32_t));
    }
    offset += num_split_saddles;
    if (num_minima > 0 && indices_buf) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, minima_buf);
        glBindBuffer(GL_COPY_WRITE_BUFFER, indices_buf);
        glCopyBufferSubData(GL_SHADER_STORAGE_BUFFER, GL_COPY_WRITE_BUFFER, 0, offset * sizeof(uint32_t), num_minima * sizeof(uint32_t));
    }
    offset += num_minima;
    if (num_join_saddles > 0 && indices_buf) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, join_saddles_buf);
        glBindBuffer(GL_COPY_WRITE_BUFFER, indices_buf);
        glCopyBufferSubData(GL_SHADER_STORAGE_BUFFER, GL_COPY_WRITE_BUFFER, 0, offset * sizeof(uint32_t), num_join_saddles * sizeof(uint32_t));
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    
    delete_buffer(maxima_buf);
    delete_buffer(split_saddles_buf);
    delete_buffer(minima_buf);
    delete_buffer(join_saddles_buf);
    
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

    {
		// Ensure all writes are done before reading back
        GLsync fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glClientWaitSync(fence, GL_SYNC_FLUSH_COMMANDS_BIT, GL_TIMEOUT_IGNORED);
        glDeleteSync(fence);
    }

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
    
    MD_LOG_INFO("Extracted %u edges in Morse-Smale complex", num_edges);
    
    success = true;

    // End timing
    glEndQuery(GL_TIME_ELAPSED);

    // Retrieve the result (blocking until GPU finishes)
    GLuint64 elapsedTime = 0;
    glGetQueryObjectui64v(query, GL_QUERY_RESULT, &elapsedTime); // nanoseconds

	MD_LOG_INFO("Topology computation of [%u,%u,%u] GPU time: %.3f ms", grid->dim[0], grid->dim[1], grid->dim[2], elapsedTime / 1e6);
    
cleanup:
    delete_buffer(ubo_buf);
    delete_buffer(ascending_buf);
    delete_buffer(descending_buf);
    delete_buffer(types_buf);
    delete_buffer(counts_buf);
    delete_buffer(counter_buf);
    
    return success;
}

#else

// macOS stub implementation
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const md_grid_t* grid) {
    (void)out_graph;
    (void)vol_tex;
    (void)grid;
    MD_LOG_ERROR("Topology computation not supported on macOS");
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
    }
}

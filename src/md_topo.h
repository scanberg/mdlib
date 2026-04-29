#pragma once
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

struct md_grid_t;
struct md_allocator_i;

// Topology analysis for scalar fields
// Constructs extremum graphs using Morse-Smale complex decomposition

// Critical point types
typedef enum md_topo_critical_point_type_t {
    MD_TOPO_UNDEFINED = 0,      // Sentinel value (should not appear in extremum graph)
    MD_TOPO_MAXIMUM = 1,
    MD_TOPO_SPLIT_SADDLE = 2,
    MD_TOPO_MINIMUM = 3,
    MD_TOPO_JOIN_SADDLE = 4,
} md_topo_critical_point_type_t;

static const char* md_topo_critical_point_type_str[] = {
    "Undefined",
    "Maximum",
    "Split Saddle",
    "Minimum",
    "Join Saddle"
};

// Edge in the Morse-Smale complex
// Edges connect critical points along integral lines:
// - Maxima -> Split saddles (descending 1-manifolds)
// - Split saddles -> Minima (descending 2-manifolds)
// - Minima -> Join saddles (ascending 1-manifolds)
// - Join saddles -> Maxima (ascending 2-manifolds)
typedef struct md_topo_edge_t {
    uint32_t from;  // Source vertex index
    uint32_t to;    // Target vertex index
} md_topo_edge_t;

typedef struct md_topo_vert_t {
    float x, y, z; // World-space position
    float value;   // Scalar field value at the vertex
} md_topo_vert_t;

// Result structure containing the Morse-Smale complex
// Vertices are stored contiguously: [maxima, split_saddles, minima, join_saddles]
// Edges encode connectivity between critical points
typedef struct md_topo_extremum_graph_t {
    // Vertex data
    uint32_t  num_vertices;       // Total number of critical points
    md_topo_vert_t* vertices;     // Vertex array [num_vertices]
    
    // Count by type (for efficient indexing into vertex arrays)
    uint32_t num_maxima;          // Maxima are at indices [0, num_maxima)
    uint32_t num_split_saddles;   // Split saddles at [num_maxima, num_maxima + num_split_saddles)
    uint32_t num_minima;          // Minima at [num_maxima + num_split_saddles, ...)
    uint32_t num_join_saddles;    // Join saddles at [... + num_minima, num_vertices)
    
    // Edge connectivity (Morse-Smale complex structure)
    uint32_t        num_edges;
    md_topo_edge_t* edges;        // Edge list [num_edges]

    struct md_allocator_i* alloc; // Allocator used for this structure
} md_topo_extremum_graph_t;

#ifdef __cplusplus
extern "C" {
#endif

#if MD_ENABLE_GPU

// Opaque handle for an in-flight GPU topology computation.
// Owns all intermediate GPU buffers and the submission fence.
typedef struct md_topo_gpu_work md_topo_gpu_work_t;

// Asynchronously compute the extremum graph on the GPU.
//   Phase 1 (manifold + critical-point detection) blocks internally to read back
//   vertex counts, which are required to size the phase-2 allocations.
//   Phase 2 (compaction + graph extraction + readback copy) is submitted
//   asynchronously and the function returns immediately.
// Returns NULL on failure. The caller owns the returned handle.
md_topo_gpu_work_t* md_topo_compute_extremum_graph_gpu(
    md_gpu_device_t device,
    md_gpu_image_t  volume,
    const struct md_grid_t* grid,
    float scalar_threshold);

// Returns true if the GPU has finished phase 2 (non-blocking poll).
bool md_topo_gpu_work_is_done(const md_topo_gpu_work_t* work);

// Block until GPU work completes, optionally extract results into out_graph, then free work.
// Pass out_graph = NULL to discard results (cancel-like).
bool md_topo_gpu_work_finish(md_topo_extremum_graph_t* out_graph, md_topo_gpu_work_t* work);

#else
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const struct md_grid_t* grid, float scalar_threshold);
#endif

// Free an extremum graph structure
void md_topo_extremum_graph_free(md_topo_extremum_graph_t* graph);

// Copy an extremum graph structure (deep copy)
void md_topo_extremum_graph_copy(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* src_graph);

// Total counts
static inline size_t md_topo_total_critical_points(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_vertices : 0;
}

static inline size_t md_topo_total_edges(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_edges : 0;
}

// Count by type
static inline size_t md_topo_count_maxima(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_maxima : 0;
}

static inline size_t md_topo_count_split_saddles(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_split_saddles : 0;
}

static inline size_t md_topo_count_minima(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_minima : 0;
}

static inline size_t md_topo_count_join_saddles(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_join_saddles : 0;
}

// Vertex offsets by type
static inline size_t md_topo_offset_maxima(const md_topo_extremum_graph_t* graph) {
    (void)graph;
    return 0;
}

static inline size_t md_topo_offset_split_saddles(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_maxima : 0;
}

static inline size_t md_topo_offset_minima(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_maxima + graph->num_split_saddles : 0;
}

static inline size_t md_topo_offset_join_saddles(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_maxima + graph->num_split_saddles + graph->num_minima : 0;
}

static inline md_topo_critical_point_type_t md_topo_vertex_type(const md_topo_extremum_graph_t* graph, size_t vertex_idx) {
    if (!graph || vertex_idx >= graph->num_vertices) return MD_TOPO_UNDEFINED;
    if (vertex_idx < graph->num_maxima) return MD_TOPO_MAXIMUM;
    if (vertex_idx < graph->num_maxima + graph->num_split_saddles) return MD_TOPO_SPLIT_SADDLE;
    if (vertex_idx < graph->num_maxima + graph->num_split_saddles + graph->num_minima) return MD_TOPO_MINIMUM;
    return MD_TOPO_JOIN_SADDLE;
}

// Extract vertex types into separate array (returns number of types written)
static inline size_t md_topo_extract_vertex_types(int* out_types, size_t capacity, const md_topo_extremum_graph_t* graph) {
    if (!graph || !out_types) return 0;
    size_t total_cp = md_topo_total_critical_points(graph);
    if (capacity < total_cp) return 0;

    size_t offset = 0;
    for (size_t i = 0; i < graph->num_maxima; ++i) {
        out_types[offset++] = MD_TOPO_MAXIMUM;
    }
    for (size_t i = 0; i < graph->num_split_saddles; ++i) {
        out_types[offset++] = MD_TOPO_SPLIT_SADDLE;
    }
    for (size_t i = 0; i < graph->num_minima; ++i) {
        out_types[offset++] = MD_TOPO_MINIMUM;
    }
    for (size_t i = 0; i < graph->num_join_saddles; ++i) {
        out_types[offset++] = MD_TOPO_JOIN_SADDLE;
    }
    return offset;
}

#ifdef __cplusplus
}
#endif

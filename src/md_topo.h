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

#define MD_TOPO_NUM_TYPES 5

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

// Result structure containing the Morse-Smale complex.
// vertices[i] holds position + scalar value; types[i] holds the critical-point type.
// Edge indices (from, to) are indices into the vertices/types arrays.
typedef struct md_topo_extremum_graph_t {
    md_topo_vert_t*                vertices;    // [num_vertices] position + value
    md_topo_critical_point_type_t* types;       // [num_vertices] per-vertex type
    uint32_t                       num_vertices;

    md_topo_edge_t*                edges;       // [num_edges]
    uint32_t                       num_edges;

    struct md_allocator_i*         alloc;
} md_topo_extremum_graph_t;

#ifdef __cplusplus
extern "C" {
#endif

#if MD_ENABLE_GPU

// Persistent context for GPU topology computation.
// Create once per volume size; reuse across multiple computations.
// All scratch and result buffers are pre-allocated at create time.
typedef struct md_topo_gpu_context md_topo_gpu_context_t;

// Allocate a context sized for a dim_x × dim_y × dim_z volume.
// Result buffers are pre-allocated to a worst-case vertex/edge capacity.
// Returns NULL on failure.
md_topo_gpu_context_t* md_topo_gpu_context_create(md_gpu_device_t device, uint32_t dim_x, uint32_t dim_y, uint32_t dim_z);

// Release all GPU resources owned by the context.
void md_topo_gpu_context_destroy(md_topo_gpu_context_t* context);

// Record the full topology pipeline into caller-supplied cmd:
//   bidirectional manifold → path compression → critical-point detection
//   → compaction setup → compaction → vertex/edge extraction → staging copies.
// Submit cmd with NULL fence (implicit sync), then call get_result.
void md_topo_gpu_cmd_record(md_gpu_command_buffer_t cmd, md_topo_gpu_context_t* context,
    md_gpu_image_t volume, const struct md_grid_t* grid, float scalar_threshold);

// Call after cmd has been submitted and completed.
// Reads CPU-visible staging into out_graph (vertices, types, edges).
// Returns false if zero critical points were found (out_graph left unchanged).
bool md_topo_gpu_context_get_result(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context);

#else
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const struct md_grid_t* grid, float scalar_threshold);
#endif

// Free an extremum graph structure
void md_topo_extremum_graph_free(md_topo_extremum_graph_t* graph);

// Copy an extremum graph structure (deep copy)
void md_topo_extremum_graph_copy(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* src_graph);

// Simplify an extremum graph into out_graph.
// - threshold: vertices with value < threshold are killed (pass 0 to skip).
// - prune_duplicate_saddles: for each pair of maxima connected by more than one
//   split-saddle, keep only the highest-value saddle and remove the rest.
// out_graph->alloc must be set before calling; it is freed and rebuilt.
void md_topo_simplify(md_topo_extremum_graph_t* out_graph, const md_topo_extremum_graph_t* in_graph,
    float threshold, bool prune_duplicate_saddles);

void md_topo_count_vertex_types(uint32_t out_counts[MD_TOPO_NUM_TYPES], const md_topo_extremum_graph_t* graph);

static inline size_t md_topo_num_critical_points(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_vertices : 0;
}

static inline size_t md_topo_num_edges(const md_topo_extremum_graph_t* graph) {
    return graph ? graph->num_edges : 0;
}

static inline md_topo_critical_point_type_t md_topo_vertex_type(const md_topo_extremum_graph_t* graph, size_t vertex_idx) {
    if (!graph || !graph->types || vertex_idx >= graph->num_vertices) return MD_TOPO_UNDEFINED;
    return graph->types[vertex_idx];
}

#ifdef __cplusplus
}
#endif

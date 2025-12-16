#pragma once
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

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

// Compute topology on GPU from a 3D scalar field volume texture
// - out_graph: Output extremum graph structure (allocated inside the function)
//      Should have alloc field set to the desired allocator (0 for default heap allocator)
// - vol_tex: The texture handle to the volume (must be a 3D texture with float format)
// - grid: The grid defining the volume dimensions and spacing
// - scalar_threshold: Minimum scalar value to consider for critical points (to filter noise)
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const struct md_grid_t* grid, float scalar_threshold);

// Free an extremum graph structure
void md_topo_extremum_graph_free(md_topo_extremum_graph_t* graph);

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

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

static inline const char* md_topo_critical_point_type_str(int type) {
    static const char* table[] = {
        "Undefined",
        "Maximum",
        "Split Saddle",
        "Minimum",
        "Join Saddle"
    };
    return (type >= 0 && type <= (int)MD_TOPO_JOIN_SADDLE) ? table[type] : table[MD_TOPO_UNDEFINED];
}

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

void md_topo_gpu_initialize(md_gpu_device_t device);
void md_topo_gpu_shutdown(void);

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

// Issue the full topology pipeline into `stream`:
//   bidirectional manifold → path compression → critical-point detection
//   → compaction → vertex/edge extraction → copies into host-readable memory.
// Everything is ordered by the stream, so no barriers or resource declarations
// are needed. `volume` must have MD_GPU_TEX_STORAGE usage; mip 0 is read.
// Wait for the stream (md_gpu_stream_sync, or a recorded sync) before calling
// md_topo_gpu_context_extract.
void md_topo_gpu_record(md_gpu_stream_t stream, md_topo_gpu_context_t* context,
    md_gpu_texture_t volume, const struct md_grid_t* grid, float scalar_threshold);

// Call once the work issued by md_topo_gpu_record has completed.
// Reads the host-readable results into out_graph (vertices, types, edges).
// Returns false if zero critical points were found (out_graph left unchanged).
bool md_topo_gpu_context_extract(md_topo_extremum_graph_t* out_graph, md_topo_gpu_context_t* context);

#else
bool md_topo_compute_extremum_graph_GPU(md_topo_extremum_graph_t* out_graph, uint32_t vol_tex, const struct md_grid_t* grid, float scalar_threshold);
#endif

// ---------------------------------------------------------------------------
// Certified critical points of a GTO electron density (CPU reference)
// ---------------------------------------------------------------------------
// Grid-free alternative to the volume based extraction above. Works directly on
//     rho(r) = sum_{mu,nu} D_{mu nu} phi_mu(r) phi_nu(r)
// and PROVES, cube by cube, that a region holds no critical point or exactly one (Newton then
// polishes it to |dx| < 1e-13 Bohr). Every non-degenerate critical point with rho >= rho_min is found;
// near-degenerate topology that cannot be resolved down to h_min is reported, never guessed.
// Vertices are in Bohr; types follow md_topo_critical_point_type_t (SPLIT_SADDLE = bond CP (3,-1),
// JOIN_SADDLE = ring CP (3,+1), MINIMUM = cage CP (3,+3)). Edges run from a saddle to the extremum its
// separatrix ends in (bond paths for BCPs, ring lines for RCPs), as in the volume pipeline.
// out_graph->alloc must be set (as for md_topo_simplify).
// Multithreaded and deterministic: the output is bit-identical for any thread count. Blocking; run it
// from a worker thread for interactive use and cancel through 'cancel' (returns false, info.cancelled).

struct md_gto_basis_t;

// How the enclosures treat rho (CPU and GPU sweeps, the same rule per cube / batch). The factored form
// writes D = sum_k l_k c_k c_k^T (pivoted LDL^T at setup, residual checked to rounding level) and bounds
// rho = sum_k l_k (c_k . phi)^2: per cube r factor rows instead of n local AO rows, so the D products
// cost r n instead of n^2 (r = rank of D: the occupied orbitals of an SCF density). Its remainder bounds
// are looser (more cubes), so it pays only when r is well below the local AO count.
typedef enum md_topo_gto_density_form_t {
    MD_TOPO_GTO_DENSITY_AUTO = 0,        // factored where it is cheaper, per cube
    MD_TOPO_GTO_DENSITY_MATRIX,          // sum_ij D_ij phi_i phi_j over the local AOs
    MD_TOPO_GTO_DENSITY_FACTORED,        // the factors wherever D factors (else the matrix)
} md_topo_gto_density_form_t;

typedef struct md_topo_gto_desc_t {
    const struct md_gto_basis_t* basis;  // Cartesian basis, md_gto conventions
    const float*  atom_xyz;              // atom positions in Bohr, indexed by shell.atom_idx
    size_t        atom_xyz_stride;       // bytes between positions, 0 = packed float[3]
    const double* density_matrix;        // [num_ao * num_ao] row major, Cartesian AO order of the basis
    double        rho_min;               // only critical points with rho >= rho_min are sought (0 -> 1e-4)
    double        h_min;                 // smallest cube half-width before a cube is reported unresolved (0 -> 1e-4 Bohr)
    double        h_root;                // root cube half-width (0 -> 1 Bohr)
    bool          trace_separatrices;    // trace separatrices to produce the graph edges
    uint32_t      num_threads;           // worker threads including the caller (0 -> all logical cores)
    volatile int32_t* cancel;            // optional: set non-zero from another thread to stop early
    bool          profile_gpu_kernels;   // GPU sweep only: wait for each kernel and time it (ms_gpu_*); slower, for benchmarks
    uint32_t      gpu_gemm_variant;      // GPU sweep only: GEMM tiling, for tuning: 0 picks one for the GPU, v + 1 forces tiling v
    uint32_t      gpu_fgemm_variant;     // GPU sweep only: factored GEMM tiling, likewise (0 picks by the rank of D)
    md_topo_gto_density_form_t density_form;   // see above (0 = auto)
} md_topo_gto_desc_t;

typedef struct md_topo_gto_cluster_t {
    float    lo[3], hi[3];               // bounds (Bohr) of a group of adjacent unresolved cubes
    uint32_t num_boxes;
    int32_t  degree;                     // Brouwer degree of grad rho on the boundary = sum of sign(det H) over the
                                         // CPs inside (max -1, BCP +1, RCP -1, CCP +1). 0: a cancelling pair or nothing
    float    degree_residual;            // distance of the computed degree from the nearest integer
    float    min_boundary_grad;          // smallest |grad rho| seen on the boundary
} md_topo_gto_cluster_t;

#define MD_TOPO_GTO_MAX_CLUSTERS 32

typedef struct md_topo_gto_info_t {
    uint64_t num_box_evals;
    uint64_t num_inflated_evals;
    uint32_t num_levels;
    uint32_t num_threads;
    uint32_t num_unresolved_boxes;
    uint32_t num_clusters;
    bool     clusters_truncated;
    bool     cancelled;                  // stopped through desc->cancel; the graph holds what was certified so far
    bool     complete;                   // true: no unresolved cubes, the graph holds every CP with rho >= rho_min
    int32_t  poincare_hopf;              // n_max - n_bcp + n_rcp - n_ccp (1 for a complete, isolated molecule)
    double   domain_pad;                 // derived padding (Bohr) around the atoms outside which rho < rho_min
    uint32_t density_rank;               // number of factors of D (0: D did not factor at rounding level; matrix form only)
    uint64_t num_factored_evals;         // cube evaluations in the factored form (CPU; part of num_box_evals)
    md_topo_gto_cluster_t clusters[MD_TOPO_GTO_MAX_CLUSTERS];
    // GPU sweep only (md_topo_compute_extremum_graph_gto_gpu)
    bool     used_gpu;                   // the octree sweep ran on the GPU
    uint64_t num_gpu_box_evals;          // cubes decided (or handed over) by the GPU; included in num_box_evals
    uint32_t num_escalated_boxes;        // cubes fp32 could not decide, finished on the CPU with their subtree
    uint32_t num_gpu_dispatches;
    // Wall-clock milliseconds per phase. ms_sweep is the whole octree sweep; the three after it are parts of it.
    double   ms_setup;                   // tables, screening radii, domain
    double   ms_sweep;
    double   ms_sweep_gpu_wait;          // waiting for GPU dispatches and readbacks
    double   ms_sweep_polish;            // Newton polish of the roots the GPU certified (host)
    double   ms_sweep_cpu;               // CPU levels: the escalated subtrees, or the whole sweep without a GPU
    double   ms_separatrices;
    double   ms_clusters;
    // Per GPU kernel, with desc.profile_gpu_kernels (part of ms_sweep_gpu_wait): AO values and remainders,
    // the D products, the dot products and AO sums, the tests.
    double   ms_gpu_ao;
    double   ms_gpu_gemm;
    double   ms_gpu_epilogue;
    double   ms_gpu_decide;
    uint64_t num_children_skipped;       // children of split cubes excluded by their parent's expansion, never evaluated
    uint64_t num_gpu_batches;            // batches of 8 sibling cubes the GPU evaluated
    uint64_t num_gpu_rows;               // their local AOs, summed (rows / batches = mean local AO count)
    uint64_t num_gpu_factored_batches;   // of them, evaluated in the factored form (desc.density_form)
    double   gpu_gemm_flop;              // arithmetic issued by the D-product kernel, padding included
} md_topo_gto_info_t;

bool md_topo_compute_extremum_graph_gto(md_topo_extremum_graph_t* out_graph, md_topo_gto_info_t* out_info, const md_topo_gto_desc_t* desc);

#if MD_ENABLE_GPU
// The same search with the octree sweep on the GPU, in fp32 with rigorous rounding margins (every bound
// the tests use carries the rounding error of its own computation), so the guarantee is unchanged. The
// CPU keeps what needs double precision or is inherently few: Newton polish of the certified roots, the
// cubes fp32 cannot decide (with their subtree, see num_escalated_boxes), separatrices and clusters.
// Blocking: works through 'stream' in chunks, two in flight. If the GPU cannot be used the whole sweep
// runs on the CPU (info.used_gpu false). Same inputs, outputs and determinism as above.
bool md_topo_compute_extremum_graph_gto_gpu(md_topo_extremum_graph_t* out_graph, md_topo_gto_info_t* out_info,
                                            const md_topo_gto_desc_t* desc, md_gpu_stream_t stream);

// The GEMM tilings desc.gpu_gemm_variant can select. Each sums the same products in the same order, so
// results are identical (md_topo_gto_bench --gemm-sweep checks); only speed differs, per GPU.
uint32_t    md_topo_gto_gpu_gemm_variant_count(void);
const char* md_topo_gto_gpu_gemm_variant_name(uint32_t variant);
// The tiling used for 'device' when desc.gpu_gemm_variant is 0 (measured per GPU vendor).
uint32_t    md_topo_gto_gpu_gemm_variant_auto(md_gpu_device_t device);
// The same for the factored form's GEMM (desc.gpu_fgemm_variant), whose default depends on the rank r of D.
uint32_t    md_topo_gto_gpu_fgemm_variant_count(void);
const char* md_topo_gto_gpu_fgemm_variant_name(uint32_t variant);
uint32_t    md_topo_gto_gpu_fgemm_variant_auto(uint32_t rank);
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

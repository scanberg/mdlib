#pragma once

#include <stdint.h>
#include <stddef.h>

struct md_allocator_i;
struct md_coord_stream_t;
struct md_unitcell_t;

// Returned by the nearest queries when no element was found within the supplied maximum distance
#define MD_SPATIAL_ACC_INVALID_IDX 0xFFFFFFFFU

typedef enum {
    MD_SPATIAL_ACC_FLAG_NONE = 0,

    // If this flag is not set, the internal idx is used for the callback, which is a dense range of 0..count-1.
    // If this flag is set, the supplied idx is used for the callback, which can be any arbitrary integer value per point.
    MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX = 0x1,
} md_spatial_acc_flags_t;

typedef struct md_spatial_acc_t {
    size_t num_elems;
    float* elem_x;
    float* elem_y;
    float* elem_z;
    uint32_t* elem_idx;

    // Optional per element radius in cartesian units, stored in the same (cell sorted) order as elem_x/y/z.
    // NULL if no radii were supplied upon initialization, which is equivalent to all radii being zero.
    float* elem_rad;

    size_t num_cells;
    uint32_t* cell_off;

    // Optional largest element radius within each cell (num_cells entries). NULL if no radii were supplied.
    float* cell_rad_max;

    uint64_t  cell_mask[3][16];
    uint32_t  cell_dim[3];
    float     inv_cell_ext[3];

    // Coarse tier over the same fractional frame, used to skip empty regions in unbounded queries.
    // Coarse cell c along an axis covers the fine cell range [c * cell_dim / coarse_dim, (c+1) * cell_dim / coarse_dim),
    // so the coarse grid partitions the period exactly and wraps consistently with the fine grid.
    size_t    num_coarse_cells;
    uint32_t  coarse_dim[3];
    uint32_t* coarse_count;     // Number of elements within each coarse cell
    float*    coarse_rad_max;   // Largest element radius within each coarse cell (NULL if no radii)

    // Largest element radius over all elements (0 if no radii were supplied)
    float max_rad;

    float G00, G11, G22;
    float H01, H02, H12;

    float A[3][3];
    float I[3][3];
    float origin[3];

    uint32_t flags;

    struct md_allocator_i* alloc;
} md_spatial_acc_t;

// Description used to initialize a spatial acceleration structure.
// Zero initialize and fill in the fields of interest, the zero value is a valid default for every optional field.
typedef struct md_spatial_acc_desc_t {
    // Input positions of the points, in cartesian coordinates. Required.
    const struct md_coord_stream_t* coords;

    // Optional per element radius in cartesian units. NULL means all radii are zero.
    // Indexed exactly like the coordinates of the stream, i.e. through the same optional indirection:
    // element i reads radii[coords->idx[i]] when the stream supplies an index array, and radii[i] otherwise.
    const float* radii;

    // Cell extent for the spatial acceleration structure. If 0, it is automatically determined.
    double cell_ext;

    // Unit cell information for periodic boundary conditions. If NULL, no unit cell is used.
    const struct md_unitcell_t* unitcell;

    md_spatial_acc_flags_t flags;
} md_spatial_acc_desc_t;

#ifdef __cplusplus
extern "C" {
#endif

// Callback signatures
// Internally the query procedures batch up results and pass them in larger chunks to the callback for better performance.
// The arrays passed to the callback are padded so vectorized loads can safely be performed without worrying about out-of-bounds access,
// the num_points/num_pairs parameter indicates the actual number of valid points/pairs in the arrays for the callback to process.

// Callback for single point query.
typedef void (*md_spatial_acc_point_callback_t)(const uint32_t* idx, const float* x, const float* y, const float* z, size_t num_points, void* user_param);

// Callback for pairwise interactions
typedef void (*md_spatial_acc_pair_callback_t)(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param);

// Initialize a spatial acceleration structure
// - in_x, in_y, in_z:  Input positions of points. They are expected to be in cartesian coordinates.
// - in_idx (optional): Input indices for points to extract. If NULL, it is assumed to be a dense range of 0..count-1.
// If not NULL, the supplied indices are used for the callback if the flag MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX is set, otherwise the linear idx is used for the callback.
// - count: Number of points OR number of indices if in_idx is not NULL.
// - cell_ext (optional): Cell extent for the spatial acceleration structure. If 0, it is automatically determined.
// - unitcell (optional): Unit cell information for periodic boundary conditions. If NULL, no unit cell is used.
// - flags: (optional) Flags to control the behavior of the spatial acceleration structure. See md_spatial_acc_init_flags_t for details.
void md_spatial_acc_init(md_spatial_acc_t* acc, const struct md_coord_stream_t* coords, double cell_ext, const struct md_unitcell_t* unitcell, md_spatial_acc_flags_t flags);

// Initialize a spatial acceleration structure from a description, which additionally allows per element radii to be supplied.
void md_spatial_acc_init_desc(md_spatial_acc_t* acc, const md_spatial_acc_desc_t* desc);

// Free the data allocated for the spatial acceleration structure. This should be called when the spatial acceleration structure is no longer needed to free the allocated memory.
void md_spatial_acc_free(md_spatial_acc_t* acc);

// --- INTERNAL PAIR TESTS ---

// Perform full N^2 test of points within the spatial acceleration structure for a supplied radius
// It is recommended that the radius <~ cell_ext, then a tight neighbor search is performed
void md_spatial_acc_for_each_internal_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

// Iterate over each point within the spatial acceleration structure within a 1-cell neighborhood (periodic if applicable)
void md_spatial_acc_for_each_internal_pair_in_neighboring_cells(const md_spatial_acc_t* acc, md_spatial_acc_pair_callback_t callback, void* user_param);

// --- EXTERNAL PAIR TESTS ---

// Test external points against internal points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices in the callback
void md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(const md_spatial_acc_t* acc, const struct md_coord_stream_t* ext_coords, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags);

// Perform a spatial query for points within the spatial acceleration structure within a bounding box defined by center and half extent (radius)
void md_spatial_acc_for_each_point_in_aabb(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param);

void md_spatial_acc_for_each_point_in_sphere(const md_spatial_acc_t* acc, const double center[3], double radius, md_spatial_acc_point_callback_t callback, void* user_param);

// --- NEAREST ELEMENT QUERY ---

// Find, for each supplied query point, the element which minimizes the additively weighted distance
//
//     d(p, i) = |p - c_i| - R_i
//
// where c_i is the position of element i and R_i its radius (zero if no radii were supplied upon initialization).
// This is the signed distance to the surface of the sphere of element i, so it is negative for a point inside an element.
// Without radii it degenerates to the plain euclidean nearest neighbor.
//
// - points:   Query points in cartesian coordinates. They do not have to lie within the unit cell.
// - max_dist: Search is bounded to this distance. Query points with no element within it report
//             MD_SPATIAL_ACC_INVALID_IDX and a distance of max_dist. Pass a large value for an unbounded search,
//             the cost of which scales with the distance actually travelled to reach the nearest element.
// - out_idx:  (optional) Element index per query point, following the same convention as the callbacks, i.e. the
//             supplied coord stream index if MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX was set, otherwise the linear index.
// - out_dist: (optional) Additively weighted distance per query point.
//
// The query only reads from acc and keeps its state on the stack, so it is safe to call concurrently from several
// threads on the same structure. It exploits spatial coherence within the supplied batch of points: supplying a
// compact block of points (a voxel block for example) is considerably cheaper than supplying scattered ones.
void md_spatial_acc_query_nearest(const md_spatial_acc_t* acc, const struct md_coord_stream_t* points, double max_dist, uint32_t* out_idx, float* out_dist);

// --- HELPER FUNCTIONS ---

// Helper functions for partial functionality
static inline void md_spatial_acc_cell_range(uint32_t out_cell_range[2], const md_spatial_acc_t* acc, size_t cell_idx) {
	out_cell_range[0] = 0;
	out_cell_range[1] = 0;
    if (cell_idx < acc->num_cells) {
        out_cell_range[0] = acc->cell_off[cell_idx];
        out_cell_range[1] = acc->cell_off[cell_idx + 1];
    }
}

// Linear cell index from cell coordinates. The layout is x major, z minor, matching the internal ordering.
static inline size_t md_spatial_acc_cell_index(const md_spatial_acc_t* acc, uint32_t cx, uint32_t cy, uint32_t cz) {
    return ((size_t)cz * acc->cell_dim[1] + (size_t)cy) * acc->cell_dim[0] + (size_t)cx;
}

// Cell coordinates from a linear cell index
static inline void md_spatial_acc_cell_coord(uint32_t out_coord[3], const md_spatial_acc_t* acc, size_t cell_idx) {
    const size_t c0  = acc->cell_dim[0];
    const size_t c01 = (size_t)acc->cell_dim[0] * acc->cell_dim[1];
    out_coord[0] = (uint32_t)((cell_idx % c01) % c0);
    out_coord[1] = (uint32_t)((cell_idx % c01) / c0);
    out_coord[2] = (uint32_t)(cell_idx / c01);
}

// Perpendicular thickness of a single cell along the supplied axis, in cartesian units.
// This is the distance between the two cell faces which are normal to the reciprocal axis, not the length of the
// cell edge, so it is the conservative measure to use when reasoning about search radii in a triclinic frame.
static inline double md_spatial_acc_cell_extent(const md_spatial_acc_t* acc, int axis) {
    const double d = (double)acc->inv_cell_ext[axis] * (double)acc->cell_dim[axis];
    return d > 0.0 ? 1.0 / d : 0.0;
}

// Cartesian position of the element at the supplied (cell sorted) element index
static inline void md_spatial_acc_elem_pos(float out_pos[3], const md_spatial_acc_t* acc, size_t elem_idx) {
    const float s[3] = { acc->elem_x[elem_idx], acc->elem_y[elem_idx], acc->elem_z[elem_idx] };
    out_pos[0] = acc->A[0][0] * s[0] + acc->A[1][0] * s[1] + acc->A[2][0] * s[2] + acc->origin[0];
    out_pos[1] = acc->A[0][1] * s[0] + acc->A[1][1] * s[1] + acc->A[2][1] * s[2] + acc->origin[1];
    out_pos[2] = acc->A[0][2] * s[0] + acc->A[1][2] * s[1] + acc->A[2][2] * s[2] + acc->origin[2];
}

#ifdef __cplusplus
}
#endif

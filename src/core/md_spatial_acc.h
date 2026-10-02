#pragma once

#include <stdint.h>
#include <stddef.h>

struct md_allocator_i;
struct md_coord_stream_t;
struct md_unitcell_t;

typedef enum {
    MD_SPATIAL_ACC_FLAG_NONE = 0,

    // If this flag is not set, the internal idx is used for the callback, which is a dense range of 0..count-1.
    // If this flag is set, the supplied idx is used for the callback, which can be any arbitrary integer value per point.
    MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX = 0x1,
} md_spatial_acc_flags_t;

#define MD_SPATIAL_ACC_MAX_TIERS 10

typedef struct md_spatial_acc_t {
    size_t num_elems;
    float* elem_x;
    float* elem_y;
    float* elem_z;
    uint32_t* elem_idx;

    size_t num_cells;
    uint32_t* cell_off;

    uint32_t  cell_dim[3];
    float     inv_cell_ext[3];

    float G00, G11, G22;
    float H01, H02, H12;

    float A[3][3];
    float I[3][3];
    float origin[3];

    uint32_t flags;

    // Cell index (see CELL INDEX in md_spatial_acc.c). Only occupied cells are stored: cell_off is indexed by occupied
    // cell and num_cells counts those. The cells are grouped 4x4x4 into the nodes of tier 1, those 4x4x4 into the
    // nodes of tier 2 and so on. Only the top tier (num_tiers) is a dense grid, the tiers below hold occupied nodes
    // only. A node has a 64 bit mask of its occupied children and the index of its first child; a child's index is
    // base + the number of mask bits below it.
    uint32_t  num_tiers;
    uint32_t  top_dim[3];
    uint64_t* top_mask;                                 // [top_dim[0] * top_dim[1] * top_dim[2]]
    uint32_t* top_base;
    uint64_t* tier_mask[MD_SPATIAL_ACC_MAX_TIERS];      // Tiers 1 .. num_tiers - 1, [0] is unused
    uint32_t* tier_base[MD_SPATIAL_ACC_MAX_TIERS];

    struct md_allocator_i* alloc;
} md_spatial_acc_t;

// Description of a spatial acceleration structure, the only way to initialize one. Zero initialize and fill in what
// applies: zero is a valid default for every field but the coordinates, so fields can be added without touching
// callers which do not need them.
typedef struct md_spatial_acc_desc_t {
    // The points, in cartesian coordinates. Required.
    const struct md_coord_stream_t* coords;

    // The largest cutoff the structure will be queried with. The cells follow from it and from how the points are
    // spread: as large as the cutoff where the points are dense, larger (up to 4 times the cutoff) where they are
    // sparse, which keeps the queries busy rather than visiting cells with a point or two. Any query up to this cutoff
    // works; up to twice it still works, at a cost. 0 means 6 A. Point queries (AABB) take any cutoff of the scale of
    // their boxes.
    double cutoff;

    // Periodic boundary conditions. NULL for none.
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

// Initialize a spatial acceleration structure from its description. Builds from scratch: call again to rebuild for
// new coordinates. acc->alloc has to be set.
void md_spatial_acc_init(md_spatial_acc_t* acc, const md_spatial_acc_desc_t* desc);

// Free the data allocated for the spatial acceleration structure. This should be called when the spatial acceleration structure is no longer needed to free the allocated memory.
void md_spatial_acc_free(md_spatial_acc_t* acc);

// --- INTERNAL PAIR TESTS ---

// Every pair of points of the structure within cutoff of each other (minimum image), each pair once.
// The cutoff may reach at most two cells: cutoff <= 2 * cell extent.
void md_spatial_acc_for_each_internal_pair_within_cutoff(const md_spatial_acc_t* acc, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param);

// --- EXTERNAL PAIR TESTS ---

// Test external points against internal points within the spatial acceleration structure for a supplied cutoff
// The external points are not part of the spatial acceleration structure and will be represented in the callback as the 'i' indices and the internal points are the 'j' indices in the callback
void md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(const md_spatial_acc_t* acc, const struct md_coord_stream_t* ext_coords, double cutoff, md_spatial_acc_pair_callback_t callback, void* user_param, md_spatial_acc_flags_t flags);

// Perform a spatial query for points within the spatial acceleration structure within a bounding box defined by center and half extent (radius)
// The coordinates handed to the callback are cartesian, in the periodic images nearest the image of aabb_cen
// which md_spatial_acc_aabb_query_center reports - which is not aabb_cen itself when that lies outside the cell.
void md_spatial_acc_for_each_point_in_aabb(const md_spatial_acc_t* acc, const double aabb_cen[3], const double aabb_rad[3], md_spatial_acc_point_callback_t callback, void* user_param);

// The image of an AABB query centre the query works in: folded into the cell along the periodic axes.
// A caller which relates the returned coordinates to its own centre offsets them by (center - out_center),
// a lattice vector. Deriving it independently (e.g. wrapping with md_util_pbc) can disagree by a whole cell
// for a centre on a cell face.
void md_spatial_acc_aabb_query_center(double out_center[3], const md_spatial_acc_t* acc, const double center[3]);

#ifdef __cplusplus
}
#endif

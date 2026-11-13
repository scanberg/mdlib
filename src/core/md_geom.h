#pragma once

// Coordinates: distances, minimum images and periodic boundaries, centers of mass, bounding boxes, alignment (rmsd)
// and the interpolation of frames.
//
// Where coordinates come with a stride, it is the distance in bytes from one point to the next, and 0 means packed
// floats: three for a position (md_float3_t), one for a weight. Points of md_float4_t are (const float*)xyzw with a
// stride of sizeof(md_float4_t), and their w as weights is (const float*)xyzw + 3 with the same stride.

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_vec.h>

struct md_unitcell_t;

#ifdef __cplusplus
extern "C" {
#endif

// Computes an array of distances between two sets of coordinates in a periodic domain (cell)
// out_dist:  Output array of distances, must have length of (num_a * num_b)
// coord_a:   Array of coordinates (a)
// num_a:     Length of coord_a
// coord_b:   Array of coordinates (b)
// num_b:     Length of coord_b
// cell:      Periodic boundary cell
void md_geom_distance_array(float* out_dist_arr, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const struct md_unitcell_t* cell);

// The minimum (maximum) distance between each group of points of a and the points of b, under the minimum image
// convention of cell (NULL, or a cell without periodic axes: plain euclidean). The minimum image is exact in triclinic
// cells too, also where distances reach past half the cell (see MINIMUM AND MAXIMUM DISTANCE BETWEEN SETS in md_geom.c).
// Group g is coord_a[a_offsets[g] .. a_offsets[g + 1]), a_offsets holds num_groups + 1 entries.
// out_dist:   The distance per group, 0 for an empty group or an empty b
// out_idx_a:  (optional) Per group, the index into coord_a of the nearest (farthest) pair, -1 where there is none
// out_idx_b:  (optional) Per group, the index into coord_b of the pair, -1 where there is none
// Each group starts from the pair of the one before, so groups which follow each other in space are cheaper.
void md_geom_min_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* coord_a, const size_t* a_offsets, size_t num_groups, const md_float3_t* coord_b, size_t num_b, const struct md_unitcell_t* cell);
void md_geom_max_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const md_float3_t* coord_a, const size_t* a_offsets, size_t num_groups, const md_float3_t* coord_b, size_t num_b, const struct md_unitcell_t* cell);

// Every pair whose distance is within tol of the minimum (maximum), in one pass with the same culling: for a caller
// which computes the distances another way, where they may differ from these by some ulps, and wants the extreme by
// its own. Up to cap pairs are written, as indices into coord_a and coord_b, in no particular order. Returns how many
// there are, or more than cap when they do not all fit (the ones written are then only some of them); 0 when a or b
// is empty.
size_t md_geom_min_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const struct md_unitcell_t* cell);
size_t md_geom_max_distance_near(int64_t* out_idx_a, int64_t* out_idx_b, size_t cap, float tol, const md_float3_t* coord_a, size_t num_a, const md_float3_t* coord_b, size_t num_b, const struct md_unitcell_t* cell);

// Reduces separation vectors to their shortest periodic image in cell, in place. Without a cell, or a cell without
// periodic axes, they are left as they are. The image of x relative to a reference r is r + min_image(x - r).
void md_geom_min_image(float* in_out_dx, size_t stride, size_t count, const struct md_unitcell_t* cell);

// Wraps coordinates into the cell, in place: all of them, or those of in_idx (count of them) when given. Without a
// cell, or a cell without periodic axes, they are left as they are.
void md_geom_pbc(float* in_out_xyz, size_t stride, const int32_t* in_idx, size_t count, const struct md_unitcell_t* cell);

// Makes connected atoms whole across the cell, in place: atom in_idx[i] is moved to the minimum image of itself
// relative to atom in_parent[i], in order, so a parent has to come before the atoms reached from it. That is the tree
// of a molecule (md_molecule_t: atom_idx and parent_idx, a root its own parent and left where it is), and
// md_system_unwrap does it for every molecule of a system. Linear, allocation free, and without consulting bonds.
// Without a cell, or a cell without periodic axes, nothing moves.
void md_geom_unwrap(float* in_out_xyz, size_t stride, const int32_t* in_idx, const int32_t* in_parent, size_t count, const struct md_unitcell_t* cell);

// PERIODIC IMAGE SELECTION
//
// These let the estimator choose each point's periodic image, using its own objective as the
// criterion. Prefer them over unwrapping followed by estimation: unwrap has to commit to an image
// assignment before it knows what is being measured, and it commits using topology, which is
// unrelated to the measurement. These do not need topology at all, which is also why they work for
// sparse selections and for points which are not atoms.
//
// @NOTE: md_geom_unwrap remains the right tool for making a whole connected structure whole - for
// rendering or export - where there is no estimator and no reference to align against, and where the
// structure may be more extended than half a cell. That is the only case it is for.

// Places a set of points in mutually consistent images and reports its centre.
//
// Seeded from the circular mean (which requires no image assignment, so it cannot be skewed by the
// images the input arrived in) and then alternated to a fixed point: place every point in the image
// nearest the centre, recompute the centre, repeat. Settles in one or two passes.
//
// The set is left in the image it arrived in: the points are made mutually consistent and the
// cluster as a whole is NOT moved into the reference cell, so the reported centre can be combined
// with coordinates that were not passed in. Call md_geom_pbc afterwards if the cell is wanted.
// A set straddling a boundary has no image to preserve and may come back on either side.
//
// Correct while the set spans less than half a cell. Beyond that no per point image choice can
// represent it and the structure has to be unwrapped topologically instead.
void md_geom_deperiodize_self(md_float4_t* in_out_xyzw, size_t count, const struct md_unitcell_t* cell, md_float3_t* out_com);

// Optimal rigid rotation between a reference set and a target set under periodic boundaries.
//
// Solves for the rotation, the target centre AND each target point's periodic image together:
//     min over R, c, n_k  of  sum_k w_k | R (q_k + A n_k - c) - (p_k - ref_com) |^2
// Given R and c the best n_k is the image nearest the predicted position, per point and in closed
// form; given the n_k the best R and c is ordinary Kabsch. Alternating the two decreases the
// objective monotonically and terminates when no point changes image.
//
// out_rot maps the TARGET frame onto the REFERENCE frame: R * (q - out_com) ~= p - ref_com.
// out_trg_xyzw optionally receives the placed target points (weights preserved); pass NULL if only
// the transform is wanted, and it may alias trg_xyzw.
//
// out_com and the placed points are expressed in the image the TARGET arrived in, not folded into
// the reference cell - see md_geom_deperiodize_self. ref_com is used only as the reference
// set's own centre, so its image does not matter.
//
// Returns the largest per point residual, which is the margin to an ambiguous image choice. Small
// against half a cell means the assignment is nowhere near flipping; approaching it means this
// frame's alignment should not be trusted.
float md_geom_optimal_rotation_pbc(md_float3x3_t* out_rot, md_float3_t* out_com, md_float4_t* out_trg_xyzw,
                                   const md_float4_t* ref_xyzw, md_float3_t ref_com,
                                   const md_float4_t* trg_xyzw, size_t count,
                                   const struct md_unitcell_t* cell, int max_iter, float tol);

// The optimal rotation (Kabsch) of two sets of points already relative to their centres, weighted by w
md_float3x3_t md_geom_optimal_rotation_rel(const md_float4_t* ref_rel_xyzw, const md_float4_t* trg_rel_xyzw, size_t count);

// The axis aligned bounding box of points: all of them, or those of in_idx (count of them) when given. A margin, as
// the radius of the largest atom, is the caller's to add.
void md_geom_aabb(float out_min[3], float out_max[3], const float* in_xyz, size_t stride, const int32_t* in_idx, size_t count);

// The center of mass of points: all of them, or those of in_idx (count of them) when given, by their weights in_w
// (NULL: equal weights). In a periodic cell it is the circular mean, which needs no image of the points to be chosen.
md_float3_t md_geom_com(const float* in_xyz, size_t xyz_stride, const float* in_w, size_t w_stride, const int32_t* in_idx, size_t count, const struct md_unitcell_t* cell);

// The root mean square deviation of two sets of points, the second rotated optimally onto the first about the
// centres in_com, by the mean of the weights (w) of each pair. With in_idx, pair i is in_idx[0][i] and in_idx[1][i].
double md_geom_rmsd(const md_float4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const md_float3_t in_com[2]);

// Computes linear shape descriptor weights (linear, planar, isotropic) from a covariance matrix
md_float3_t md_geom_shape_weights(const md_float3x3_t* covariance_matrix);

// Perform linear interpolation of supplied coordinates
// out_xyz:     Destination, packed
// in_xyz:      Sources [2], packed
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t: interpolation factor (0..1)
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
void md_geom_interpolate_linear(md_float3_t* out_xyz, const md_float3_t* const in_xyz[2], size_t count, const struct md_unitcell_t* cell, float t);

// The weights of four consecutive frames at t (0..1) for a cubic scheme over them, on the interval between the
// two middle ones: out = w[0] in[0] + w[1] in[1] + w[2] in[2] + w[3] in[3]. They sum to one.
// Cardinal spline with tangent scale s (0.5 is Catmull-Rom): C1, and passes through in[1] at t = 0 and in[2] at t = 1.
void md_geom_cubic_spline_weights(float out_w[4], float t, float s);
// Uniform cubic B-spline with the four frames as its control points: C2, and damps frame to frame jitter, but it
// approximates rather than interpolates: at t = 0 it is (in[0] + 4 in[1] + in[2]) / 6, not in[1].
void md_geom_cubic_bspline_weights(float out_w[4], float t);

// Weighted sum of four consecutive frames of coordinates, out = w[0] in[0] + w[1] in[1] + w[2] in[2] + w[3] in[3]:
// any cubic scheme over them at a given t (see the weights above). In a periodic cell in[0] and in[2] are taken at
// the minimum image of in[1], and in[3] at that of in[2].
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
void md_geom_interpolate_cubic_weighted(md_float3_t* out_xyz, const md_float3_t* const in_xyz[4], size_t count, const struct md_unitcell_t* cell, const float w[4]);

#ifdef __cplusplus
}
#endif

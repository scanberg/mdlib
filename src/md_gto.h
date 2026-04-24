#pragma once
#include <stdint.h>
#include <stddef.h>

#include <core/md_grid.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

// Gaussian Type Orbital
// Are evaluated as f(x',y',z') = (x'-x)^i (y'-y)^j (z'-z)^k c exp(-a ((x'-x)^2 + (y'-y)^2 + (z'-z)^2))
// Where x' y' and z' are the observer coordinates we evaluate the function at
typedef struct md_gto_t {
	float x;
	float y;
	float z;
	float coeff;		// Baked coefficient (Including normalization factors)
	float alpha;		// Exponent alpha
	float cutoff;		// Radial cutoff
	uint8_t i, j, k, l;
	uint32_t _pad;
} md_gto_t;

typedef struct md_pgto_t {
	float coeff;
	float alpha;
	float radius;
	uint8_t i, j, k, l;
} md_pgto_t;

// Orbital Data 
typedef struct md_orbital_data_t {
	size_t num_gtos;
	md_gto_t* gtos;

	// Optional parameters to evaluate multiple orbitals in one 'go'
	// This is necessary in case the evaluation mode is psi squared as the squaring has to occur after each orbital
	size_t num_orbs;
	uint32_t* orb_offsets;
	float* orb_scaling;
} md_orbital_data_t;

// ---------------------------------------------------------------------------
// General basis representation
// ---------------------------------------------------------------------------
// A contracted shell: a set of (l+1)(l+2)/2 cartesian angular functions at
// one atom center, sharing num_primitives radial primitives.
// All coefficients are expected to be fully normalized (contraction + angular
// normalization factors baked in).  The struct contains no atom coordinates —
// those are passed as a separate argument to any procedure that needs them,
// keeping this struct static after construction.
// 16 bytes, 4 per cache line, no internal padding.
typedef struct md_gto_shell_t {
    uint32_t  atom_idx;           // index into the caller-supplied atom_xyz array
    uint32_t  primitive_offset;   // first index into basis alpha/coeff arrays
    uint32_t  l;                  // angular momentum (0=s, 1=p, 2=d, ...)
    uint32_t  num_primitives;
} md_gto_shell_t;

// A complete contracted Gaussian basis set for a molecular system.
// Shells and primitives are owned by this struct (single allocation recommended).
// Atom coordinates are NOT stored here; pass them on the side as const float* atom_xyz
// (packed xyz, length num_atoms*3) to any procedure that requires them.
typedef struct md_gto_basis_t {
    uint32_t         num_shells;
    uint32_t         num_primitives;
    md_gto_shell_t*  shells;   // [num_shells]
    float*           alpha;    // [num_primitives]  Gaussian exponents
    float*           coeff;    // [num_primitives]  normalized contraction coefficients
} md_gto_basis_t;

// ---------------------------------------------------------------------------
// Compute-optimized (derived) representation
// ---------------------------------------------------------------------------
// Flat SoA layout ready for CPU/GPU evaluation. Built from md_gto_basis_t via
// md_gto_data_build(). Each CGTO corresponds to one cartesian angular component
// of a shell; primitives are purely radial (ijkl lives at the CGTO level).
typedef struct md_gto_data_t {
	size_t     num_cgtos;
	vec4_t*    cgto_xyzr;
	uint32_t*  cgto_offset;		// offsets into pgtos (contains num_cgtos + 1 entries) such that a 'range' can be represented by cgto_offset[i] -> cgto_offset[i+1]

	size_t     num_pgtos;
	float*     pgto_alpha;
	float*     pgto_coeff;
	float*     pgto_radius;
	uint32_t*  pgto_ijkl;
} md_gto_data_t;

typedef enum {
	MD_GTO_EVAL_MODE_PSI = 0,
	MD_GTO_EVAL_MODE_PSI_SQUARED = 1,
} md_gto_eval_mode_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline uint32_t md_gto_pack_ijkl(int i, int j, int k, int l) {
	uint32_t res = 0;
	res |= ((uint32_t)i) << 0;
	res |= ((uint32_t)j) << 8;
	res |= ((uint32_t)k) << 16;
	res |= ((uint32_t)l) << 24;
	return res;
}

static inline void md_gto_unpack_ijkl(uint32_t packed, int* i, int* j, int* k, int* l) {
	if (i) *i = (packed >> 0) & 0xFF;
	if (j) *j = (packed >> 8) & 0xFF;
	if (k) *k = (packed >> 16) & 0xFF;
	if (l) *l = (packed >> 24) & 0xFF;
}

// Evaluates GTOs over a grid on the GPU and stores the result into a supplied volume
// - vol_tex: The texture handle to the volume
// - vol_dim: The dimensions of the volume
// - vol_step: The voxel spacing in world space length units
// - world_to_model: float[4][4] (col-major) transformation matrix to transform a point in world_space coordinates into the volumes model space (note not texture space, but a space which is rotated and translated such that the axes align with the volume and its origin is placed at (0,0,0))
// - index_to_world: float[4][4] (col-major) transformation matrix to transform a point in the volumes index coordinates [0, dim[ into world space coordinates
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

void md_gto_grid_evaluate_orb_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_orbital_data_t* orb, md_gto_eval_mode_t mode);

// Average local ionization energy
void md_gto_grid_evaluate_ALIE_GPU(uint32_t vol_tex, const md_grid_t* vol_grid, const md_orbital_data_t* orb, md_gto_eval_mode_t mode);

// This is malplaced at the moment, but this is for the moment, the best match in where to place the functionality
// Performs voronoi segmentation of the supplied volume to points with a supplied radius and accumulates the value of each voxel into the corresponding point
// - out_values: Destination array holding the point values that are written to, should have length 'num_points'
// - point_xyzr: Point coordinates + radius, packed xyzrxyzrxyzr
// - num_points: Number of points
// - vol_tex: The texture handle to the volume
// - vol_grid: The grid defining the volume
void md_gto_voronoi_segment_GPU(float* out_values, const float* point_xyzr, size_t num_points, uint32_t vol_tex, const md_grid_t* vol_grid);

// Evaluates GTOs over a grid
// - out_grid_values: The grid to write the evaluated values to, should have length 'grid->dim[0] * grid->dim[1] * grid->dim[2]'
// - grid: The grid to evaluate
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate(float* out_grid_values, const md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

// Evaluates CGTOs defined by gto_data with a given 'density' matrix
// - out_grid_values: The grid to write the evaluated values to, should have length 'grid->dim[0] * grid->dim[1] * grid->dim[2]'
// - grid: The grid defining the samples to evaluate
// - gto_data: The gto data to evaluate
// - matrix: The matrix coefficients (Upper triangular format)
// - matrix_dim: The dimension of the square matrix (should match the number of CGTOs in gto_data)
void md_gto_grid_evaluate_matrix(float* out_grid_values, const md_grid_t* grid, const md_gto_data_t* gto_data, const float* matrix_data, size_t matrix_dim);

// Evaluates CGTOs defined by gto_data with a given 'density' matrix
// - vol_tex: The texture handle to the volume
// - grid: The grid defining the location of samples to evaluate
// - gto_data: The gto data to evaluate
// - matrix_data: The matrix coefficients (Upper triangular format)
// - matrix_dim: The dimension of the square matrix (should match the number of CGTOs in gto_data)
void md_gto_grid_evaluate_matrix_GPU(uint32_t vol_tex, const md_grid_t* grid, const md_gto_data_t* gto_data, const float* matrix_data, size_t matrix_dim, bool include_gradients);

#if MD_ENABLE_GPU
// GPU-resident buffer holding all input data for density evaluation.
// Sized at creation time from gto_data; only the matrix sub-region can be updated in-place.
// The caller must have waited on (or destroyed) the previous dispatch fence before calling
// md_gto_density_buf_update or md_gto_density_buf_destroy.
typedef struct md_gto_density_buf_t {
    md_gpu_buffer_t buffer;  // single backing GPU buffer
    uint32_t num_cgtos;
    uint32_t num_pgtos;
    uint32_t matrix_dim;     // == num_cgtos; stored separately for convenience
} md_gto_density_buf_t;

// Allocate a GPU buffer and upload all GTO basis data + initial density matrix.
// Returns false on allocation failure.
bool md_gto_density_buf_create(md_gto_density_buf_t* buf, md_gpu_device_t device, const md_gto_data_t* gto_data, const float* matrix_data, size_t matrix_dim);

// Re-upload the density matrix into the existing buffer.
// matrix_dim must match the value used at creation time.
// Call only after the previous dispatch fence has been waited on or destroyed.
void md_gto_density_buf_update_matrix(md_gto_density_buf_t* buf, const float* matrix_data, size_t matrix_dim);

// Destroy the backing GPU buffer. Safe to call with a zeroed struct.
void md_gto_density_buf_destroy(md_gto_density_buf_t* buf);

// Dispatch an async electron density evaluation.
// The image must be R32_FLOAT with MD_GPU_IMAGE_STORAGE.
// *out_fence receives a fence the caller must wait on (md_gpu_fence_wait) or poll
// (md_gpu_fence_is_signaled), then destroy with md_gpu_destroy_fence.
bool md_gto_grid_evaluate_density_gpu(md_gpu_device_t device, md_gto_density_buf_t* buf, md_gpu_image_t out_image, const md_grid_t* grid, md_gpu_fence_t* out_fence);
#endif

// Evaluate GTOs over subportion of a grid
// - out_grid_values: The grid to write the evaluated values to, should have length 'grid->dim[0] * grid->dim[1] * grid->dim[2]'
// - grid: The grid to evaluate a subportion of
// - grid_idx_off: Index offset for x,y,z
// - grid_idx_len: Index length for x,y,z
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
// @NOTE: It is strongly recommended that the evaluation occurs over a 8x8x8 blocks as this will get a vectorized fastpath
void md_gto_grid_evaluate_sub(float* out_grid_values, const md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t eval_mode);

// Evaluates GTOs over a set of given XYZ coordinates
// - out_xyz_values: Array of values to write evaluated values to, should have length 'num_xyz'
// - xyz: Pointer to base addr of xyz packed coordinates to be evaluated, should have length 'num_xyz'
// - num_xyz: Number of coordinates to evaluate
// - stride_xyz [OPTIONAL]: Stride in bytes between the given xyz coordinates. A value of zero implies fully packed XYZXYZ... -> (12 bytes)
// - gtos: input gtos to be evaluated for the supplied coordinates
// - num_gtos: Number of gtos
// - eval_mode: GTO evaluation mode
void md_gto_xyz_evaluate(float* out_xyz_values, const float* xyz, size_t num_xyz, size_t stride_xyz, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t eval_mode);

// Evaluates GTOs over space 
//void md_gto_xyz_voronoi_evaluate(float* out_val, const float* xyz, size_t num_xyz, size_t stride_xyz, const md_gto_t* gtos, size_t)

// Compute the cutoff parameter within the supplied GTOs based on the given value
// Typically this could be somewhere around 1.0e-6
// returns the number of gtos left after this, as some may have no radius of influence and will be pruned away
size_t md_gto_cutoff_compute_and_filter(md_gto_t* gtos, size_t num_gtos, double value);

// Calculate the radius of influence for a single GTO
double md_gto_compute_radius_of_influence(int i, int j, int k, double coeff, double alpha, double cutoff);

// Extracts a subset of gtos from an input array which overlap a given axis aligned bounding box with its radii of influence
// - out_gtos: The output array to store the gtos in
// - min_ext: The minimum extent of the axis aligned bounding box
// - max_ext: The maximum extent of the axis aligned bounding box

size_t md_gto_aabb_test(md_gto_t* out_gtos, const float aabb_min[3], const float aabb_max[3], const md_gto_t* in_gtos, size_t num_gtos);

// Extracts a subset of gtos from an input array which overlap a given oriented bounding box with its radii of influence
size_t md_gto_obb_test(md_gto_t* out_gtos, const float obb_center[3], const float obb_half_ext[3], const float obb_orientation[3][3], const md_gto_t* in_gtos, size_t num_gtos);

#ifdef __cplusplus
}
#endif

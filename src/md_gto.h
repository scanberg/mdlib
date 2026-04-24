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
// One contracted shell in radial (angular-momentum) form.
// A shell with angular momentum l contributes (2l+1) spherical AOs to the
// MO coefficient vector, or (l+1)(l+2)/2 Cartesian components internally.
// Cartesian angular indices (i,j,k) are derived from l at evaluation time —
// they are NOT stored here, keeping this struct independent of any particular
// Cartesian ordering convention.
// Atom coordinates are NOT stored — they are passed separately so the basis
// remains valid after geometry updates.
// 16 bytes, fits exactly 4 per cache line, no padding.
typedef struct md_gto_shell_t {
    uint32_t  atom_idx;          // index of the parent atom in the caller-supplied atom_xyz array
    uint32_t  primitive_offset;  // first index into basis->alpha and basis->coeff arrays
    uint32_t  l;                 // angular momentum quantum number (0=s, 1=p, 2=d, 3=f, ...)
    uint32_t  num_primitives;    // number of Gaussian primitives in the contraction
} md_gto_shell_t;

// A complete contracted Gaussian basis set for a molecular system,
// stored in the radial (shell) form — one entry per contracted shell.
//
// Shells are ordered to match the AO index convention of the QM program that
// generated the MO coefficient vectors.  The AO index for shell i is simply
// the sum of (2*shells[j].l + 1) for all j < i (spherical counting); this is
// implicit and does not need to be stored.
//
// Contraction coefficients (coeff) are pure radial — any spherical-to-Cartesian
// transformation factors are applied internally at evaluation time (md_gto_data_build).
//
// Atom coordinates are NOT stored here.  Pass them as a flat float array
// (packed xyz, stride 3, length >= max(atom_idx)+1) to any function that needs them.
typedef struct md_gto_basis_t {
    uint32_t         num_shells;      // total number of contracted shells
    uint32_t         num_primitives;  // total number of Gaussian primitives (sum of num_primitives over all shells)
    md_gto_shell_t*  shells;          // [num_shells]
    float*           alpha;           // [num_primitives]  Gaussian exponents (bohr^-2)
    float*           coeff;           // [num_primitives]  normalized radial contraction coefficients
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

// ---------------------------------------------------------------------------
// Basis-centric evaluation API  (primary interface)
// ---------------------------------------------------------------------------
// All evaluation functions take the static basis topology plus the dynamic
// data that changes at runtime: atom positions and electronic-state data
// (MO coefficient vector or AO density matrix).
//
// atom_xyz : packed xyz in Bohr, stride 3 floats, length >= max(shell.atom_idx)+1
// cutoff   : value threshold for pruning GTOs by radius of influence (e.g. 1e-6);
//            pass 0 to disable pruning.
// Coefficients are double to preserve QM-code precision at the boundary.

// Evaluate a single molecular orbital on a grid (psi or psi^2).
// mo_coeffs: AO coefficient vector, length = total num_cgtos implied by basis.
void md_gto_grid_evaluate_mo(
    float* out, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* mo_coeffs, double cutoff, md_gto_eval_mode_t mode);

// Evaluate sum_i weight_i * psi_i(r)^2 on a grid.
// Used for NTO densities, attachment/detachment densities, and MO-based
// electron density.  Each mo_coeffs[i] is an AO coefficient vector of
// length num_cgtos; squaring and weighted accumulation are done internally.
void md_gto_grid_evaluate_mo_sum(
    float* out, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* const* mo_coeffs, const double* weights, size_t num_mos,
    double cutoff);

// Evaluate full electron density via AO density matrix on a grid.
// density_matrix: full N×N row-major double array, N = number of CGTOs implied by basis.
// The function constructs a compact upper-triangular float representation internally.
void md_gto_grid_evaluate_density(
    float* out, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix);

// Returns the upper bound on the number of Cartesian GTOs produced when
// fully expanding all shells in the basis (before any cutoff filtering).
// Use this to pre-allocate the output buffer for md_gto_expand_with_mo().
size_t md_gto_pgto_count(const md_gto_basis_t* basis);

// Expand a contracted GTO basis with baked MO coefficients into a flat md_gto_t[] array.
// mo_coeffs: one double per spherical AO (length = sum of (2*l+1) over all shells)
// atom_xyz:  packed xyz in Bohr, stride 3 floats, length >= max(shell.atom_idx)+1
// cutoff:    if > 0, applies md_gto_cutoff_compute_and_filter; pass 0 to skip
// out:       caller-allocated array of length >= md_gto_pgto_count(basis)
// Returns:   number of elements written (after pruning, or total if cutoff == 0)
size_t md_gto_expand_with_mo(md_gto_t* out, const md_gto_basis_t* basis,
    const float* atom_xyz, const double* mo_coeffs, double cutoff);

// Basis-centric GL density evaluation.
// Expands the basis, converts the NxN double density matrix to upper-triangular float,
// and dispatches the existing GL density compute shader.
// density_matrix: full N×N row-major doubles, N = number of spherical AOs in basis.
void md_gto_grid_evaluate_density_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix, bool include_gradients);

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
// Evaluates GTOs over a grid
// - out_grid_values: The grid to write the evaluated values to, should have length 'grid->dim[0] * grid->dim[1] * grid->dim[2]'
// - grid: The grid to evaluate
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate(float* out_grid_values, const md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

// Evaluates CGTOs defined by gto_data with a given 'density' matrix
// - vol_tex: The texture handle to the volume
// - grid: The grid defining the location of samples to evaluate
// - gto_data: The gto data to evaluate
// - matrix_data: The matrix coefficients (Upper triangular format)
// - matrix_dim: The dimension of the square matrix (should match the number of CGTOs in gto_data)
void md_gto_grid_evaluate_matrix_GPU(uint32_t vol_tex, const md_grid_t* grid, const md_gto_data_t* gto_data, const float* matrix_data, size_t matrix_dim, bool include_gradients);

#if MD_ENABLE_GPU
// GPU-resident buffer holding all input data for density evaluation.
// Created once from a basis + initial positions + initial density matrix.
// Static basis data (shell topology, primitives) is uploaded at creation and
// never needs to change.  Positions and the density matrix can be updated
// independently and cheaply via the update functions.
// The caller must have waited on (or destroyed) the previous dispatch fence
// before calling any update or destroy function.
typedef struct md_gto_density_buf_t {
    md_gpu_buffer_t buffer;  // single backing GPU buffer
    uint32_t num_cgtos;
    uint32_t num_pgtos;
    uint32_t matrix_dim;     // == num_cgtos
} md_gto_density_buf_t;

// Allocate and populate a GPU buffer from a basis, atom positions and density matrix.
// density_matrix: full N×N row-major doubles; converted to upper-triangular floats internally.
bool md_gto_density_buf_create(md_gto_density_buf_t* buf, md_gpu_device_t device,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* density_matrix);

// Re-upload atom positions (cgto_xyzr sub-region) into the existing buffer.
void md_gto_density_buf_update_positions(md_gto_density_buf_t* buf,
    const md_gto_basis_t* basis, const float* atom_xyz);

// Re-upload the density matrix into the existing buffer.
// density_matrix: full N×N row-major doubles; converted to upper-triangular floats internally.
void md_gto_density_buf_update_matrix(md_gto_density_buf_t* buf,
    const double* density_matrix);

// Destroy the backing GPU buffer. Safe to call on a zeroed struct.
void md_gto_density_buf_destroy(md_gto_density_buf_t* buf);

// Dispatch an async electron density evaluation.
// The image must be R32_FLOAT with MD_GPU_IMAGE_STORAGE.
// *out_fence receives a fence the caller must wait on (md_gpu_fence_wait) or poll
// (md_gpu_fence_is_signaled), then destroy with md_gpu_destroy_fence.
bool md_gto_grid_evaluate_density_gpu(md_gpu_device_t device, md_gto_density_buf_t* buf,
    md_gpu_image_t out_image, const md_grid_t* grid, md_gpu_fence_t* out_fence);
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

#ifdef __cplusplus
}
#endif

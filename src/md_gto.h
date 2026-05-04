#pragma once
#include <stdint.h>
#include <stddef.h>

#include <core/md_grid.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

// Stand Alone Gaussian Type Orbital
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
// One contracted Gaussian shell.
// A shell with angular momentum l contributes (2l+1) contiguous coefficients
// to an MO vector (spherical AOs) and (2l+1) contiguous rows/columns to an
// AO-basis matrix (density, overlap, etc.).
// Cartesian angular indices (i,j,k) are derived from l at evaluation time.
// Atom coordinates are NOT stored — they are passed separately so the basis
// remains valid after geometry updates.
// 16 bytes, fits exactly 4 per cache line, no padding.
typedef struct md_gto_shell_t {
    uint32_t  atom_idx;          // index of the parent atom in the caller-supplied atom_xyz array
    uint32_t  primitive_offset;  // first index into basis->alpha and basis->coeff arrays
    uint32_t  l;                 // angular momentum quantum number (0=s, 1=p, 2=d, 3=f, ...)
    uint32_t  num_primitives;    // number of Gaussian primitives in the contraction
} md_gto_shell_t;

// A complete contracted Gaussian basis set for a molecular system.
// One md_gto_shell_t per contracted shell, ordered as: angular momentum first,
// then atom, then contracted function within that (angl, atom) pair.
// MO coefficient vectors and AO-basis matrices must be in the same shell order.
//
// Contraction coefficients (coeff) are pure radial — spherical-to-Cartesian
// transformation factors are applied internally at evaluation time
// to produce the final Cartesian GTO coefficients used in the evaluation.

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

// Evaluate full electron density via AO density matrix on a grid.
// density_matrix: full symmetric N×N row-major double array, N = number of CGTOs implied by basis.
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

// GPU-accelerated versions of the above evaluation functions.  See md_gto.c for details on the expected data layout and GPU buffer formats.
void md_gto_grid_evaluate_mo_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* mo_coeffs, double cutoff, md_gto_eval_mode_t mode);

// mo_scl is optional and if null is supplied, then it is assumed that all orbitals should be scaled by 1.0 (i.e. no relative scaling between orbitals).
void md_gto_grid_evaluate_multi_mo_GL(uint32_t vol_tex, const md_grid_t* grid,
    const md_gto_basis_t* basis, const float* atom_xyz,
    const double* mo_coeffs[], const double mo_scl[], size_t num_mos, double cutoff, md_gto_eval_mode_t mode);

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

// Average local ionization energy
// Evaluates GTOs over a grid
// - out_grid_values: The grid to write the evaluated values to, should have length 'grid->dim[0] * grid->dim[1] * grid->dim[2]'
// - grid: The grid to evaluate
// - gtos: The gtos to evaluate
// - num_gtos: Number of supplied gtos
// - eval_mode: GTO evaluation mode
void md_gto_grid_evaluate(float* out_grid_values, const md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode);

#if MD_ENABLE_GPU
// ---------------------------------------------------------------------------
// GPU basis buffer  (owns device-local atom_idx, r, off_len, pgto regions)
// Create once per basis set; reuse across all density and MO evaluations.
// ---------------------------------------------------------------------------

typedef struct md_gto_gpu_basis* md_gto_gpu_basis_t;

typedef struct md_gto_gpu_basis_desc_t {
	const md_gto_basis_t* basis;      // basis set topology (required)
	double                cutoff;     // pgto radius pruning threshold; pass 0 to disable
} md_gto_gpu_basis_desc_t;

// Allocate a device-local basis buffer and fully populate it from desc.
// Ready for dispatch immediately after creation; no separate upload step required.
md_gto_gpu_basis_t md_gto_gpu_basis_create(md_gpu_device_t device, const md_gto_gpu_basis_desc_t* desc);

void md_gto_gpu_basis_destroy(md_gto_gpu_basis_t basis_buf);

// Return the underlying GPU basis buffer (CGTO metadata + PGTO data).
md_gpu_buffer_t md_gto_gpu_basis_buffer(md_gto_gpu_basis_t basis_buf);

// Basis metadata queries.
uint32_t md_gto_gpu_basis_num_cgtos(md_gto_gpu_basis_t basis_buf);
uint32_t md_gto_gpu_basis_num_pgtos(md_gto_gpu_basis_t basis_buf);
uint32_t md_gto_gpu_basis_num_atoms(md_gto_gpu_basis_t basis_buf);

// Atom buffer helpers.
// Atom buffer layout expected by the GPU shaders:
//   float4[num_atoms], xyz in Bohr and .w ignored.
size_t md_gto_gpu_atom_buffer_size(size_t num_atoms);
void md_gto_gpu_atom_pack(float* dst_atom_xyzw, const float* atom_xyz, size_t num_atoms);

// ---------------------------------------------------------------------------
// Coefficient buffer helpers
// Coefficients live in a plain md_gpu_buffer_t owned by the caller.
// Tightly packed floats, no padding.
//   density : float[(num_cgtos*(num_cgtos+1))/2]  — packed upper-triangular
//   MO      : float[num_mos * num_cgtos]           — packed row-major
// ---------------------------------------------------------------------------

// Returns the required buffer size in bytes.
size_t md_gto_gpu_coeff_size_density(size_t num_cgtos);
size_t md_gto_gpu_coeff_size_mo(size_t num_mos, size_t num_cgtos);

// Pack helpers — convert double-precision coefficients to the tightly-packed float
// layout expected by the GPU shaders.  dst must point to at least
// md_gto_gpu_coeff_size_*() bytes of writable memory; the caller decides where
// that memory comes from (heap, temp arena, staging bump cpu_ptr, …).
//
//   density: upper-triangular row-major, length = num_cgtos*(num_cgtos+1)/2
//   MO:      row-major [num_mos × num_cgtos]; mo_scales may be NULL (defaults to 1.0)
void md_gto_gpu_coeff_pack_density(float* dst, const double* density_matrix, size_t num_cgtos);
void md_gto_gpu_coeff_pack_mo(float* dst, const double* const* mo_coeffs, const double* mo_scales, size_t num_mos, size_t num_cgtos);

// Upload helpers — record a buffer copy (src_buf[src_offset] -> coeff_buf) plus a
// TRANSFER→COMPUTE barrier into cmd.
// On UMA, pack directly into md_gpu_buffer_cpu_ptr(coeff_buf) and skip these calls.
// On discrete GPU, pack into a staging cpu_ptr (e.g. from md_gpu_staging_bump_alloc),
// then call upload to record the copy from that staging buffer into coeff_buf.
void md_gto_gpu_coeff_upload_density(md_gpu_command_buffer_t cmd, md_gpu_buffer_t coeff_buf, md_gpu_buffer_t src_buf, size_t src_offset, size_t num_cgtos);
void md_gto_gpu_coeff_upload_mo(md_gpu_command_buffer_t cmd, md_gpu_buffer_t coeff_buf, md_gpu_buffer_t src_buf, size_t src_offset, size_t num_mos, size_t num_cgtos);

// ---------------------------------------------------------------------------
// GPU dispatch
// ---------------------------------------------------------------------------

// Record an electron density evaluation dispatch into the caller's command buffer.
// atom_buf must contain packed float4 atom positions (xyz in Bohr).
// coeff_buf must contain packed upper-triangular float coefficients
// (use md_gto_gpu_coeff_upload_density to fill it).
// A TRANSFER→COMPUTE barrier is inserted before the dispatch so uploads recorded
// earlier in the same command buffer are visible to the shader.
void md_gto_grid_evaluate_density_gpu(md_gpu_command_buffer_t cmd,
	md_gto_gpu_basis_t basis_buf, md_gpu_buffer_t atom_buf, md_gpu_buffer_t coeff_buf,
    md_gpu_image_t out_image, const md_grid_t* grid);

// Record an MO evaluation dispatch into the caller's command buffer.
// atom_buf must contain packed float4 atom positions (xyz in Bohr).
// coeff_buf must contain num_mos packed rows of num_cgtos floats
// (use md_gto_gpu_coeff_upload_mo to fill it).
// eval_mode controls whether psi or psi^2 is accumulated per MO row.
void md_gto_grid_evaluate_mo_gpu(md_gpu_command_buffer_t cmd,
	md_gto_gpu_basis_t basis_buf, md_gpu_buffer_t atom_buf, md_gpu_buffer_t coeff_buf, uint32_t num_mos,
    md_gpu_image_t out_image, const md_grid_t* grid, md_gto_eval_mode_t eval_mode);

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

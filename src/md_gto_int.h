#pragma once

// md_gto_int - integrals over the Cartesian GTO basis of md_gto.h
//
// The first and so far only customer is the ELECTROSTATIC POTENTIAL of a charge distribution
// made of an AO density matrix and point charges (normally the nuclei):
//
//     V(C) = sum_k q_k / |C - A_k|  +  s * sum_{mu,nu} D_{mu nu} (mu| 1/|r - C| |nu)
//
// with s the charge carried by one unit of density (-1 for electrons). The points can carry a
// dipole and a quadrupole besides their charge - the permanent multipoles of a polarizable
// embedding's sites - which add, with d = C - A_k and r = |d|,
//
//     mu_k . d / r^3  +  1/2 sum_ab Q_k,ab (3 d_a d_b - r^2 delta_ab) / r^5
//
// Q is the Cartesian second moment, NOT traceless: the Taylor convention of the polarizable
// embedding literature, V = sum_k (-1)^k / k! M^(k) . T^(k), and of md_gto_int_moments_t. Its
// trace has no potential.
//
// WHY THERE ARE NO INTEGRAL MATRICES HERE
// The potential never needs (mu|1/|r-C||nu) as a matrix, only contracted with D. So the density is
// contracted ONCE, up front, into Hermite space (McMurchie-Davidson): every product of two
// primitives is a finite sum of Hermite Gaussians Lambda_tuv(r; p, P) = d^t/dPx^t d^u/dPy^u
// d^v/dPz^v exp(-p |r - P|^2), so the whole density becomes
//
//     rho(r) = sum_g sum_{t+u+v <= L_g} h^g_tuv Lambda_tuv(r; p_g, P_g)
//
// One such GAUSSIAN g exists per unique primitive pair (atom A, exponent a, atom B, exponent b):
// shell pairs and generally contracted shells that share primitives fold into the same one, and
// the two orders of a pair are the same distribution. Its potential is then exact and cheap,
//
//     V_g(C) = sum_tuv h'^g_tuv R_tuv(p_g, P_g - C),     h' = (2 pi / p) h
//
// where R_tuv are the Hermite Coulomb integrals - one Boys function evaluation F_0..F_L and a short
// recursion. Per point that is one evaluation per gaussian instead of one per AO pair and
// primitive pair, and no matrix the size of the basis is ever formed.
//
// The same expansion gives the electric field (R up to L+1) and the multipole moments for free.
//
// THREE TIMESCALES, AS ELSEWHERE IN md_gto
//   basis                      static
//   md_gto_int_charges_t       rebuilt when the geometry OR the density changes (milliseconds)
//   evaluation                 any number of times against one md_gto_int_charges_t
//
// PRECISION
// The CPU path is double throughout and is the REFERENCE. Given double inputs its integrals agree
// with PySCF to 1e-11 au for s through g; through md_gto_basis_t, which stores exponents and
// contraction coefficients as float, expect a few 1e-7 au away from the nuclei (test_gto_int.c).
// The GPU path is float arithmetic per term with an exact fixed-point sum over the terms, ~1e-6 au
// from the CPU path; see md_gto_int_gpu_potential_launch.
//
// UNITS: bohr in, hartree / e (atomic units of potential) out, field in hartree / (e bohr).

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <md_gto.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

// Highest Hermite order a gaussian can have: the product of two shells of the highest angular
// momentum md_gto supports.
#define MD_GTO_INT_MAX_ORDER (2 * MD_GTO_MAX_ANGULAR_MOMENTUM)

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// Number of Hermite coefficients of a gaussian of order L: all (t,u,v) with t+u+v <= L.
static inline uint32_t md_gto_int_num_hermite(uint32_t L) { return ((L + 1) * (L + 2) * (L + 3)) / 6; }

// The far field's opening criterion for anything bigger than a few dozen atoms: see
// md_gto_int_charges_desc_t::far_field_theta.
#define MD_GTO_INT_FAR_FIELD_THETA 0.4

// Which far field this is, for benchmarks to print: bumped with every change to how it is built
// or evaluated.
#define MD_GTO_INT_FAR_FIELD_REVISION "r2 (order max(6, L), p r^2 >= 12, leaves of 32; GPU: a list per group of 64, points in Morton order)"

// A node of the far field's octree (md_gto_int_charges_desc_t::far_field_theta). Its members are
// item[item_begin, item_end): the gaussians and points inside its box. Its multipoles are about
// 'center', and they stand for the members wherever the distance d from the centre satisfies
// d >= radius / theta and d >= reach.
typedef struct md_gto_int_node_t {
    double   center[3];         // bohr: the centre of the members' bounding box
    double   radius;            // bohr: max |x - center| over the members' centres
    double   reach;             // bohr: max |x - center| + the extent of the member, see md_gto_int.c
    uint32_t first_child;       // the children are contiguous
    uint32_t num_children;      // 0 for a leaf
    uint32_t item_begin;
    uint32_t item_end;
} md_gto_int_node_t;

// ---------------------------------------------------------------------------
// CHARGE DISTRIBUTION
// ---------------------------------------------------------------------------
// A charge distribution in the form the evaluators consume: Hermite gaussians plus point charges.
// The CHARGE SIGN IS BAKED IN - an electron density built with density_scale = -1 has negative
// coefficients - so every evaluator simply sums, and a difference density or a density without
// nuclei needs nothing special.
//
// HERMITE COEFFICIENT ORDER within a gaussian of order L: by degree n = t+u+v ascending, and
// within a degree in the same order md_gto.h uses for the Cartesian AOs of a shell (t descending,
// then u descending). So the coefficients of degree <= n are always a prefix, and index 0 is the
// (0,0,0) term. Every coefficient already includes the factor 2*pi/p, i.e. it multiplies R_tuv
// directly.
//
// Gaussians are SORTED BY ORDER: those of order L are [order_offset[L], order_offset[L+1]).
// That is what lets a GPU kernel run branch free per order, and what md_gto_int_gpu_* relies on.
//
// The struct is transparent so that it can be inspected and serialised, but it is produced by
// md_gto_int_charges_init only: the bound and the ordering are invariants an evaluator relies on.
typedef struct md_gto_int_charges_t {
    uint32_t  num_gaussians;
    uint32_t  max_order;                                   // highest L present, 0 if none
    uint32_t  order_offset[MD_GTO_INT_MAX_ORDER + 2];      // see above
    double*   center;            // [num_gaussians * 3]  P, bohr
    double*   exponent;          // [num_gaussians]      p, bohr^-2
    double*   bound;             // [num_gaussians]      max over all space of |V_g|, hartree/e
    uint32_t* coeff_offset;      // [num_gaussians]      first coefficient of gaussian g
    size_t    num_coeffs;
    double*   coeff;             // [num_coeffs]         h'_tuv, see HERMITE COEFFICIENT ORDER

    uint32_t  num_points;
    double*   point_xyz;         // [num_points * 3]     bohr
    double*   point_charge;      // [num_points]         e
    double*   point_dipole;      // [num_points * 3]     e bohr, NULL when no point has one
    double*   point_quadrupole;  // [num_points * 6]     e bohr^2, xx xy xz yy yz zz, NULL when no point has one

    // The Boys function tabulated for the CPU evaluators (see md_gto_int.c), built with the
    // gaussians. NULL is valid - a distribution assembled any other way, or deserialised, is
    // evaluated with the series instead, to the same values, about three times more slowly.
    double*   boys_table;

    // FAR FIELD (md_gto_int_charges_desc_t::far_field_theta). num_nodes 0 when there is none, and
    // every evaluator then sums every gaussian and point at every point.
    double              far_field_theta;
    uint32_t            far_field_order;   // of the multipoles; at least max_order
    uint32_t            num_nodes;
    md_gto_int_node_t*  node;              // [num_nodes], node 0 the root
    double*             node_multipole;    // [num_nodes * num_hermite(far_field_order)], see md_gto_int.c
    uint32_t*           item;              // [num_gaussians + num_points]: a gaussian, or num_gaussians + a point
} md_gto_int_charges_t;

typedef struct md_gto_int_charges_desc_t {
    // GAUSSIAN PART: density_scale * sum_{mu,nu} density_matrix[mu][nu] phi_mu(r) phi_nu(r)
    // Optional: leave basis or density_matrix NULL for point charges only.
    const md_gto_basis_t* basis;
    const float*  atom_xyz;          // bohr, one position per atom the shells index
    size_t        atom_xyz_stride;   // bytes between positions, 0 = packed (12 bytes)
    const double* density_matrix;    // [N][N] row major, N = md_gto_basis_num_ao(basis)
    // Charge of one unit of density_matrix. 0 is taken as -1, an electron density, which is what
    // almost every caller has. Pass +1 for a density that is meant to be positive.
    double        density_scale;

    // POINT PART: sum_k point_charge[k] delta(r - point_xyz[k]), normally the nuclei.
    // Use the charge the electrons were computed against: with an effective core potential that
    // is Z minus the core electrons, not the atomic number. The QM readers publish it as
    // qm/atom/nuclear_charge (md_qm_publish_atoms).
    const float*  point_xyz;         // bohr
    size_t        point_xyz_stride;  // bytes, 0 = packed (12 bytes)
    const double* point_charge;      // e
    size_t        num_points;
    // Optional, at the same points: a dipole and a quadrupole each (see the top of the file for the
    // convention). NULL for none; a point without one has zeros.
    const double* point_dipole;      // [num_points * 3] e bohr
    const double* point_quadrupole;  // [num_points * 6] e bohr^2, xx xy xz yy yz zz, not traceless

    // Screening: a gaussian whose potential cannot exceed this anywhere (hartree/e) is dropped.
    // 0 keeps everything. 1e-8 removes a third or more of the gaussians of anything bigger than a
    // few atoms, for an error below 1e-6 au even on C60; at 1e-6 the error reaches ~1e-4 au there.
    double        threshold;

    // FAR FIELD. Far enough from its centre, a gaussian's potential is that of a point multipole
    // (erf(sqrt(p) r) = 1 to 4e-7 at p r^2 = 12), and so is that of a group of them, expanded about
    // a common centre. With far_field_theta > 0 the gaussians and points go into an octree whose
    // nodes carry the multipoles of their members, to order max(6, max_order), and an evaluator takes
    // a node whole at a point at distance d from its centre once d >= radius / far_field_theta and
    // every member is a multipole there; it evaluates the members of the nodes it cannot take that
    // way exactly. The expansion converges as far_field_theta^(order + 1).
    //
    // 0 (the default) builds nothing, and every gaussian is evaluated at every point: the reference,
    // and the right choice for a molecule of a few dozen atoms. MD_GTO_INT_FAR_FIELD_THETA (0.4), on
    // a 1152 atom protein in def2-SV(P) a few bohr outside it: ~25 times less work per point, within
    // 5e-5 au of the exact sum (the potential there spans ~0.1 au).
    double        far_field_theta;
} md_gto_int_charges_desc_t;

// Builds the distribution. Allocates from alloc; release with md_gto_int_charges_free.
//
// The density matrix need not be symmetric - only its symmetric part has a density, and that is
// what is used - so a transition density goes in as it is. It must be over the Cartesian AOs of
// the basis in the order of md_gto.h, like every AO matrix this library takes.
//
// Returns false on invalid input (inconsistent sizes, an angular momentum above
// MD_GTO_MAX_ANGULAR_MOMENTUM, a shell centred on an atom that has no position).
bool md_gto_int_charges_init(md_gto_int_charges_t* out, const md_gto_int_charges_desc_t* desc, struct md_allocator_i* alloc);

void md_gto_int_charges_free(md_gto_int_charges_t* charges, struct md_allocator_i* alloc);

// Work of evaluating the potential at ONE point, in Hermite recursion steps: what to multiply by a
// voxel count to compare grid resolutions, or to choose between the CPU and the GPU. The unit is
// the one md_gto_int_gpu_potential_desc_t::work_per_dispatch is stated in. Roughly 1-2 floating
// point operations per step.
uint64_t md_gto_int_charges_work_per_point(const md_gto_int_charges_t* charges);

// ---------------------------------------------------------------------------
// MOMENTS
// ---------------------------------------------------------------------------
// Exact moments of the whole distribution about 'origin' (bohr; NULL = the coordinate origin).
// For a neutral molecule the dipole is origin independent and is what a QM program reports as the
// dipole moment, which makes this the cheapest end to end check of a basis, a density and a reader.
typedef struct md_gto_int_moments_t {
    double charge;       // e
    double dipole[3];    // e bohr
    double second[6];    // e bohr^2: xx xy xz yy yz zz, Cartesian (NOT traceless)
} md_gto_int_moments_t;

md_gto_int_moments_t md_gto_int_charges_moments(const md_gto_int_charges_t* charges, const double origin[3]);

// ---------------------------------------------------------------------------
// CPU EVALUATION  (reference, double precision)
// ---------------------------------------------------------------------------
// Single threaded and not vectorised: ~30 ns per gaussian and point for the potential, ~130 with the
// field (26 atoms in def2-SVP, ~7500 gaussians: ~0.2 ms per point). Hand the _sub variant, or slices
// of the points, to worker threads for anything but small systems, or use the GPU path.
// At a point charge itself the potential is infinite; points closer than 1e-10 bohr to one skip
// that charge's term rather than returning inf. Everywhere else the result is exact up to the
// screening threshold the distribution was built with.

// Potential at arbitrary points (bohr).
// out_potential : [num_xyz]
// out_field     : [num_xyz * 3] optional, E = -grad V
// xyz_stride    : bytes between points, 0 = packed (12 bytes)
void md_gto_int_potential_xyz(double* out_potential, double* out_field, const float* xyz, size_t num_xyz, size_t xyz_stride,
                              const md_gto_int_charges_t* charges);

// Potential on the points of a grid, written as float, x fastest, like md_gto_grid_evaluate.
// Point (i,j,k) is origin + orientation * (((i,j,k) + sample_offset) * spacing); sample_offset
// NULL means (0,0,0), (0.5,0.5,0.5) puts the points at voxel centres as the GPU paths in viamd do.
void md_gto_int_potential_grid(float* out_values, const md_grid_t* grid, const float sample_offset[3],
                               const md_gto_int_charges_t* charges);

// The same for the sub-box [idx_off, idx_off + idx_len) of the grid, writing only those values of
// out_values (still indexed over the whole grid). The unit to hand out to worker threads.
void md_gto_int_potential_grid_sub(float* out_values, const md_grid_t* grid, const float sample_offset[3],
                                   const int idx_off[3], const int idx_len[3], const md_gto_int_charges_t* charges);

// ---------------------------------------------------------------------------
// GPU EVALUATION
// ---------------------------------------------------------------------------
#if MD_ENABLE_GPU

// Kernels are created on first use against this device. Call shutdown before destroying it.
void md_gto_int_gpu_initialize(md_gpu_device_t device);
void md_gto_int_gpu_shutdown(void);

// A charge distribution resident on the device, as float. Create once per md_gto_int_charges_t
// and evaluate it as often as needed. Uploaded on `stream`; work issued into the same stream
// afterwards sees it, another stream must md_gpu_stream_wait first.
//
// Point dipoles and quadrupoles are not evaluated on the device yet: a distribution with any is
// refused (NULL, logged) rather than evaluated without them. Use the CPU path for those.
typedef struct md_gto_int_gpu_charges* md_gto_int_gpu_charges_t;

md_gto_int_gpu_charges_t md_gto_int_gpu_charges_create(md_gpu_stream_t stream, const md_gto_int_charges_t* charges);

// Stream ordered: the memory is released at this point in `stream`.
void md_gto_int_gpu_charges_destroy(md_gpu_stream_t stream, md_gto_int_gpu_charges_t charges);

typedef struct md_gto_int_gpu_potential_desc_t {
    md_gto_int_gpu_charges_t charges;
    md_gpu_texture_t out_tex;        // 3D, needs MD_GPU_TEX_STORAGE; mip 0 is written
    const md_grid_t* grid;
    float            sample_offset[3];   // in index units, (0.5,0.5,0.5) = voxel centres
    md_gto_op_t      op;

    // Upper bound on the work of one dispatch, in units of one Hermite recursion step at one
    // voxel. The evaluation is split over as many dispatches as needed to stay under it, so that
    // a large system on a slow GPU cannot trip a driver watchdog (2 s on Windows by default).
    // 0 = default, 2^31: a few milliseconds on a discrete GPU, well under a second on anything.
    uint64_t         work_per_dispatch;

    // Optional, host memory, one byte per voxel of the grid (x fastest): nonzero where the potential
    // is wanted. With a far field the grid is evaluated in bricks of 4x4x4 voxels, and only the
    // bricks holding a wanted voxel are (all of their voxels); out_tex is left as it was elsewhere.
    // Without a far field it is ignored and every voxel is evaluated. NULL: every voxel.
    const uint8_t*   voxel_mask;
} md_gto_int_gpu_potential_desc_t;

// Issue a potential evaluation over the grid into `stream`. Asynchronous: no readbacks, no waits.
// With a far field (md_gto_int_charges_desc_t::far_field_theta) the device walks the same tree with
// the same criterion as the CPU path, for a group of 64 voxels at a time.
// Returns false, logged, when nothing was issued - invalid input, a kernel that could not be made,
// no memory for the accumulator - and out_tex is then left as it was: a caller about to read the
// texture back must not take its contents for the potential.
//
// Every term is evaluated in float and summed in 64 bit fixed point (2^-32 hartree/e resolution),
// so the sum over tens of thousands of gaussians - which cancels against the nuclei to a few
// percent of either - loses nothing in the summation, the result does not depend on how the work
// was split, and no fast-math setting can reassociate it away. Plain float summation would be off
// by ~5e-4 au on C60, a few percent of its surface potential.
bool md_gto_int_gpu_potential_launch(md_gpu_stream_t stream, const md_gto_int_gpu_potential_desc_t* desc);

// The potential at arbitrary points already on the device, the counterpart of md_gto_int_potential_xyz:
// e.g. the vertices of an isosurface mesh (md_gto_iso_gpu_mesh_t::position, xyz_stride 16), evaluated
// where they are with no readback in between.
typedef struct md_gto_int_gpu_potential_xyz_desc_t {
    md_gto_int_gpu_charges_t charges;
    md_gpu_addr_t    xyz;               // positions, bohr, float: point i at xyz + i * xyz_stride
    uint32_t         xyz_stride;        // bytes between points, a multiple of 4; 0 = packed (12)
    uint32_t         num_xyz;
    md_gpu_addr_t    out_potential;     // float [num_xyz], hartree / e
    uint64_t         work_per_dispatch; // as md_gto_int_gpu_potential_desc_t; 0 = default
} md_gto_int_gpu_potential_xyz_desc_t;

// Issue the evaluation into `stream`. Asynchronous: no readbacks, no waits. Returns false, logged, when
// nothing was issued (invalid input, a kernel that could not be made, no memory for the accumulator).
// Precision and summation as md_gto_int_gpu_potential_launch. With a far field the points are first
// sorted along a Morton curve through their box (on the device, 8 bytes of temp memory a point, to the
// next power of two), and 64 consecutive ones then share one walk of the tree against their bounding
// sphere: the order they come in does not matter, to the time or to the result.
bool md_gto_int_gpu_potential_xyz_launch(md_gpu_stream_t stream, const md_gto_int_gpu_potential_xyz_desc_t* desc);

#endif

#ifdef __cplusplus
}
#endif

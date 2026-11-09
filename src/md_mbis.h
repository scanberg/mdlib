#pragma once

// md_mbis - Minimal Basis Iterative Stockholder partitioning of a molecular electron density
//
// Verstraelen, Vandenbrande, Heidar-Zadeh, Vanduyfhuys, Van Speybroeck, Waroquier and Ayers,
// J. Chem. Theory Comput. 12, 3894 (2016).
//
// WHAT IT IS
// A real space partitioning: atom A owns the fraction w_A(r) = rho0_A(r) / rho0(r) of the density
// at every point, rho0 = sum_B rho0_B, where each PRO-ATOM is a few spherical Slater shells,
//
//     rho0_A(r) = sum_i N_Ai / (8 pi sigma_Ai^3) exp(-|r - R_A| / sigma_Ai),
//
// one per row of the periodic table the atom's electrons reach (H, He: 1; Li-Ne: 2; Na-Ar: 3 ...).
// The shells are FITTED to the molecule's own density by the fixed point
//
//     N_Ai     = int rho(r) rho0_Ai(r) / rho0(r) dr
//     sigma_Ai = 1 / (3 N_Ai) int rho(r) rho0_Ai(r) / rho0(r) |r - R_A| dr,
//
// so no reference atom densities are needed. The weights are positive and sum to one everywhere;
// the populations depend on the density alone and converge with the basis set, unlike those of
// Mulliken or Lowdin, which partition the basis functions rather than space.
//
// QUADRATURE
// Becke's multicentre scheme: per atom a Mura-Knowles radial grid times a Gauss-Legendre (cos theta)
// x uniform (phi) product grid on the sphere, the atoms' grids partitioned by Stratmann-Scuseria-
// Frisch cell functions. The density is evaluated ONCE on that grid, from the AO density matrix in
// double with the basis' shells screened per batch of points; the iterations then only touch the
// stored masses (quadrature weight x density) and the pro-atoms. Points of negligible density are
// dropped before the iterations: what they hold is below 1e-10 e an atom.
//
// THE RESULT IS THE PARTITION, NOT JUST THE CHARGES
// The fitted shells define w_A(r) anywhere, analytically: whatever else is to be partitioned the
// same way (an orbital density, a hole or particle density) is integrated against those weights.
//
// UNITS: bohr in and out, electrons, e.

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <core/md_vec_math.h>

struct md_gto_basis_t;
struct md_allocator_i;
struct md_executor_t;

typedef struct md_mbis_desc_t {
    const struct md_gto_basis_t* basis;
    const dvec3_t*  atom_coord;        // bohr, one per atom the basis' shells index (basis->shells[].atom_idx)
    const uint8_t*  atomic_number;     // per atom
    const double*   nuclear_charge;    // per atom, e: the charge the density was solved against (Z less any ECP's
                                       // core electrons). NULL: the atomic numbers. It sets the charge
                                       // (nuclear_charge - population) and, rounded, the number of shells.
    size_t          num_atoms;
    const double*   density_matrix;    // md_gto_basis_num_ao(basis) on a side, row major: the TOTAL electron density

    // Quadrature and convergence, 0 for the defaults (in brackets). The defaults hold the charges to ~2e-5 e of a
    // converged grid (measured on H to S, against a 200 x order 59 grid); the angular order is what decides it.
    uint32_t        radial_points;     // per atom [48]
    uint32_t        angular_order;     // the spherical harmonic degree integrated exactly, odd [35]
    double          threshold;         // largest change in a shell population (e) or relative width per iteration [1e-8]
    uint32_t        max_iterations;    // [1000]

    const struct md_executor_t* executor;   // NULL: the calling thread. The result does not depend on it.
    volatile int32_t* cancel;          // optional: set non-zero from another thread to stop (md_mbis_compute returns false)
    volatile float*   progress;        // optional: 0 to 1, for a progress display. An estimate that never goes back:
                                       // the grid is the first half, the iterations (by their residual) the second
} md_mbis_desc_t;

typedef struct md_mbis_result_t {
    size_t    num_atoms;
    double*   charge;              // [num_atoms] e: nuclear_charge - population
    double*   population;          // [num_atoms] e: the electrons atom A owns, the sum of its shells' N_Ai

    // The pro-atoms: atom A's shells are [shell_offset[A], shell_offset[A + 1]), innermost first.
    uint32_t  num_shells;
    uint32_t* shell_offset;        // [num_atoms + 1]
    double*   shell_population;    // [num_shells] e, N_Ai
    double*   shell_width;         // [num_shells] bohr, sigma_Ai

    double    electron_count;      // the density integrated over the grid: what the populations sum to
    uint32_t  num_points;          // grid points kept (of significant density)
    uint32_t  iterations;
    double    change;              // the last iteration's largest change (see threshold)
    bool      converged;
} md_mbis_result_t;

#ifdef __cplusplus
extern "C" {
#endif

// Partitions the density. Every array of the result is allocated from 'alloc' (md_mbis_result_free). Blocking;
// run it from a worker thread for interactive use. Returns false, with nothing allocated, on invalid input, out
// of memory or cancellation (desc->cancel, or the executor). A run that hit max_iterations returns true with converged false: its populations are the last
// iterate's.
bool md_mbis_compute(md_mbis_result_t* out, const md_mbis_desc_t* desc, struct md_allocator_i* alloc);

void md_mbis_result_free(md_mbis_result_t* result, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif

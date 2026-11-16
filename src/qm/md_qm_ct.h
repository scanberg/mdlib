#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

// Charge transfer analysis of an electronic excitation: where the excitation takes its electron
// from (the hole) and where it puts it (the particle), resolved over groups of atoms.
//
// THE QUANTITY
//
// With T the transition density matrix of an excitation in an ORTHONORMAL AO basis, the share of
// the excitation that moves an electron from group A to group B is
//
//     Omega_AB = sum_{mu in A, nu in B} T_{mu nu}^2
//
// (the charge transfer numbers of Plasser and Lischka, as computed by TheoDORE). Every entry is a sum
// of squares, so it is never negative; Omega sums to the norm of the excitation; the diagonal is
// local excitation and the off-diagonal is charge transfer. The orthonormal basis is LOWDIN's,
// S^(1/2) applied to the AOs: of all orthonormal bases it is the one closest to the AOs themselves,
// so it attributes each coefficient to the atom the AO sits on. (A Mulliken partitioning needs no
// orthogonalization, but it can turn negative and is markedly more basis dependent at atom level.)
//
// HOW IT IS COMPUTED HERE
//
// Not from T itself, which is A x A, but from the state's natural transition orbitals (NTOs),
// T = sum_k sqrt(lambda_k) h_k p_k^T. Taken to the Lowdin basis they stay orthonormal, and
//
//     Omega_AB = sum_kl s_k s_l H^A_kl P^B_kl,   s_k = sqrt(lambda_k),
//     H^A_kl = sum_{mu in A} h~_k[mu] h~_l[mu],  P^B the same over the particles,
//
// which is EXACT for the pairs given - the cross terms k != l, the interference between pairs, are
// all there. With every pair it is exactly the expression above; truncated to the pairs that matter
// it misses only what the dropped pairs carry. Dropping the cross terms instead (one product per
// pair) is NOT a small approximation: when two local excitations mix into delocalized NTOs it invents
// charge transfer that the cross terms cancel.
//
// The work is O(A K^2) for the K x K blocks per group and O(G^2 K^2) for Omega, against the O(A^2)
// that T alone would take to store.
//
// LAYOUT
//
// Row major throughout, orbitals as ROWS: an {M,A} coefficient array holds orbital m at [m * A].
// That is how the attribute table publishes orbitals (orbital/alpha/coefficient) and NTOs.

struct md_executor_t;

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// LOWDIN ORBITALS
// ---------------------------------------------------------------------------

// The orbitals in the Lowdin basis, C~ = S^(1/2) C, from C ALONE.
//
// For orbitals orthonormal in the metric S (C^T S C = I), S^(1/2) C has orthonormal columns and
// equals the orthogonal polar factor of C, C (C^T C)^(-1/2) - so neither S nor S^(1/2) is needed,
// and nothing A x A is formed. This holds equally for coefficients that were converted from a
// spherical basis to the Cartesian one, whose S is singular: the polar factor is then S^(1/2) C for
// the S those coefficients are orthonormal in.
//
// C:   {M,A} with M <= A, the orbitals as rows. Out: {M,A}, the same rows in the Lowdin basis.
// Orbitals that are linearly dependent within round off (C^T C singular) are rejected: false.
//
// Cost: two products of M x M x A and one M x M symmetric eigendecomposition. Do it once per
// geometry; every excited state reuses the result. All of it runs on 'executor' (NULL: the calling
// thread).
bool md_qm_lowdin_orbitals(double* out_C_lowdin, const double* C, size_t num_mo, size_t num_ao,
                           const struct md_executor_t* executor);

// S^(1/2) x for AO vectors x in the span of the orbitals (any combination of them, NTOs included):
// out = C~ C^T S x in column terms. x and out are {count,A}; S is {A,A}. Two vectors' worth of scratch.
void md_qm_lowdin_transform(double* out, const double* x, size_t count, const double* S,
                            const double* C, const double* C_lowdin, size_t num_mo, size_t num_ao);

// The inverse, S^(-1/2) x~ = C C~^T x~: a Lowdin basis vector back to AO coefficients, for instance
// to evaluate a block NTO on a grid. No S needed.
void md_qm_lowdin_inverse_transform(double* out, const double* x_lowdin, size_t count,
                                    const double* C, const double* C_lowdin, size_t num_mo, size_t num_ao);

// ---------------------------------------------------------------------------
// GROUP BLOCKS AND THE TRANSFER MATRIX
// ---------------------------------------------------------------------------

// The K x K block of a set of K Lowdin vectors over each group of AOs:
//     out[g][k][l] = sum_{mu : ao_group[mu] == g} x[k][mu] x[l][mu]
// x is {K,A}; out is {G,K,K}, overwritten. An AO whose group is >= num_groups belongs to none.
//
// Blocks ADD: the block of a union of groups is the sum of their blocks. So compute them once per
// atom, and any grouping of atoms into fragments is a sum, with no second pass over the AOs.
void md_qm_ct_group_blocks(double* out, const double* x, size_t K, size_t num_ao, const uint32_t* ao_group, size_t num_groups);

// The transfer matrix: out[a][b] = Omega_ab = sum_kl s_k s_l H[a][k][l] P[b][k][l], s_k = sqrt(lambda_k).
// hole_blocks {Gh,K,K} and part_blocks {Gp,K,K} from md_qm_ct_group_blocks over the hole and the
// particle NTOs; out {Gh,Gp}, overwritten. lambda[k] < 0 counts as 0.
//
// Row a sums to the hole on a, column b to the particle on b, and everything to the sum of the
// lambdas when the groups cover every AO.
void md_qm_ct_omega(double* out, const double* lambda, size_t K,
                    const double* hole_blocks, size_t num_hole_groups,
                    const double* part_blocks, size_t num_part_groups);

// ---------------------------------------------------------------------------
// BLOCK NTOs
// ---------------------------------------------------------------------------

// What one entry of Omega is made of: the singular value decomposition of the a -> b block of the
// transition density, T~_ab = sum_j sigma_j u_j v_j^T, with u_j confined to group a and v_j to group b.
// The pairs are orbitals in their own right: a hole on the donor, a particle on the acceptor, and
// their weights sigma_j^2 sum to Omega_ab exactly and are never negative. Typically one or two of
// them carry all of a transfer.
//
// hole_block_a and part_block_b are one K x K block each (as md_qm_ct_group_blocks makes them, or
// sums of them). Writes up to K weights, descending, and returns how many are above round off.
// Optional (NULL to skip): out_hole_coef and out_part_coef, {count,K} - pair j's hole is
// sum_k out_hole_coef[j][k] h~_k restricted to group a, its particle sum_k out_part_coef[j][k] p~_k
// restricted to group b, both of unit norm in the Lowdin basis.
size_t md_qm_ct_block_nto(double* out_weight, double* out_hole_coef, double* out_part_coef,
                          const double* lambda, size_t K,
                          const double* hole_block_a, const double* part_block_b);

#ifdef __cplusplus
}
#endif

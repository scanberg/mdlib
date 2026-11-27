#pragma once

#include <md_gto.h>

// ---------------------------------------------------------------------------
// md_gto_dmat_t - a density matrix in 16 x 16 tiles of the AO pairs that meet
// ---------------------------------------------------------------------------
// The form a density matrix takes in this library: only the AO pairs whose shells both reach some point (their
// centres no further apart than their radii at the cutoff, md_gto_basis_shell_radii, plus margin), never the whole
// matrix. No point has both AOs of any other pair above the cutoff, so leaving those out changes nothing the cutoff
// does not already. Kept as the tile pairs (I, J), I <= J, of 16 x 16 entries in an order of the AOs that keeps the
// tiles full: shells by radius class (2 bohr), then by atom (the input order, along a chain), then shell; a shell's
// AOs together. The GPU's md_gto_gpu_dmat_t has this same pattern, its values in float (md_gto_gpu_dmat_upload copies
// a density's tiles to the device).
//
// The pattern depends on where the atoms are. margin (bohr) also keeps the pairs up to that much further apart, so
// that it holds while no atom has moved more than margin / 2 from where it was made (md_gto_dmat_pattern_valid).
//
// Entry e of the kept tiles is the tile (I, cols[e]) of row tile I, row_ptr[I] <= e < row_ptr[I + 1], its columns
// ascending; tiles[256 e + 16 a + c] is D(perm[16 I + a], perm[16 J + c]), in the AO order of the basis, 0 where
// perm is -1 (past the last AO). A tile on the diagonal holds the whole 16 x 16 block.

struct md_allocator_i;
struct md_executor_t;

typedef struct md_gto_dmat_t {
    uint32_t  num_ao;       // Cartesian AOs of the basis
    uint32_t  num_atoms;
    uint32_t  num_shells;
    uint32_t  nt;           // tiles per edge: num_ao / 16 rounded up
    uint64_t  num_tiles;    // tile pairs kept
    double    cutoff, margin;
    float*    xyz;          // [3 num_atoms] the positions the pattern was made for (bohr)
    float*    shell_r;      // [num_shells] the shells' radii at the cutoff
    int32_t*  perm;         // [16 nt] tiled AO -> AO, -1 past num_ao
    uint32_t* tao;          // [num_ao] AO -> tiled AO
    uint32_t* shell_tao;    // [num_shells] each shell's first tiled AO; its others follow it
    uint64_t* row_ptr;      // [nt + 1]
    uint32_t* cols;         // [num_tiles] each entry's column tile (>= its row tile)
    double*   tiles;        // [256 num_tiles], or NULL before the first md_gto_dmat_build
    struct md_allocator_i* alloc;
} md_gto_dmat_t;

typedef struct md_gto_dmat_desc_t {
    const md_gto_basis_t* basis;            // required
    const float*          atom_xyz;         // bohr, one per atom of the basis
    size_t                atom_xyz_stride;  // bytes between positions; 0 for tightly packed (12 bytes)
    double                cutoff;           // > 0: of the evaluations it is read by
    double                margin;           // bohr, >= 0
} md_gto_dmat_desc_t;

#ifdef __cplusplus
extern "C" {
#endif

// Makes the pattern; the tiles are allocated by the first md_gto_dmat_build. Everything is allocated from alloc.
// False (logged) for invalid arguments or out of memory, dm then empty.
bool md_gto_dmat_init(md_gto_dmat_t* dm, const md_gto_dmat_desc_t* desc, struct md_allocator_i* alloc);

void md_gto_dmat_free(md_gto_dmat_t* dm);

// A copy of src's pattern (made for the same basis, positions, cutoff and margin); the tiles are allocated by the
// first md_gto_dmat_build, not copied. From alloc. False (logged) for invalid arguments or out of memory, dm then empty.
bool md_gto_dmat_init_like(md_gto_dmat_t* dm, const md_gto_dmat_t* src, struct md_allocator_i* alloc);

// Whether the pattern holds for these positions: the same number of atoms, none further than margin / 2 from where
// the pattern was made.
bool md_gto_dmat_pattern_valid(const md_gto_dmat_t* dm, const float* atom_xyz, size_t atom_xyz_stride);

// D = sum_k w_k c_k c_k^T on the kept tiles, in double, over num_mo orbitals: orbital k's coefficient of AO i (the
// basis order) at coeff[k * mo_stride + i * ao_stride] (rows: mo_stride num_ao, ao_stride 1; columns, [num_ao][ld]:
// mo_stride 1, ao_stride ld). Weights of either sign (0 skips the orbital), so that a difference of densities is one
// set of orbitals. In chunks of orbitals (its scratch is bounded, not num_ao x num_mo), the row tiles spread over
// executor (NULL: the calling thread); bit for bit the same for any executor. About 256 multiply-adds per tile and
// orbital. False when out of memory or cancelled through the executor (the tiles then incomplete).
bool md_gto_dmat_build(md_gto_dmat_t* dm, const double* coeff, size_t mo_stride, size_t ao_stride, const double* weight, size_t num_mo,
                       const struct md_executor_t* executor);

// The same in parts, for orbitals that are not all at hand at once (read from a file a block at a time): clear, then
// add each block in turn. D += sum_k w_k c_k c_k^T, the arguments as md_gto_dmat_build's; the blocks of a build in its
// order give it bit for bit (every tile entry sums its orbitals in order, whatever the blocks). clear allocates the
// tiles where there are none yet and zeroes them. False (logged) as md_gto_dmat_build.
bool md_gto_dmat_clear(md_gto_dmat_t* dm);
bool md_gto_dmat_add(md_gto_dmat_t* dm, const double* coeff, size_t mo_stride, size_t ao_stride, const double* weight, size_t num_mo,
                     const struct md_executor_t* executor);

// D(i, j) in the basis' AO order: 0 for a pair that is not kept (or before a build).
double md_gto_dmat_at(const md_gto_dmat_t* dm, uint32_t i, uint32_t j);

// The entry of the tile pair (I, J), I <= J: its tile is tiles + 256 e. -1 where it is not kept.
int64_t md_gto_dmat_find(const md_gto_dmat_t* dm, uint32_t I, uint32_t J);

// The block of D between shells sa and sb (na and nb their AOs, the basis' shell indices): out[ca * nb + cb] is
// D(ca of sa, cb of sb), 0 where the pair is not kept. For a consumer that works shell pair by shell pair: a shell's
// AOs are together in the tiled order, so its block lies in at most 2 x 2 tiles. Returns whether any of it is kept.
bool md_gto_dmat_shell_block(double* out, const md_gto_dmat_t* dm, uint32_t sa, uint32_t na, uint32_t sb, uint32_t nb);

// rho = sum_ij D_ij phi_i phi_j, its gradient, and tau = 1/2 sum_ij D_ij grad phi_i . grad phi_j at points, in double,
// from the kept tiles, with the AOs of the shells that reach each point (md_gto_dmat_t's radii). For ELF and LOL
// (md_qm_elf, md_qm_lol: rho, |grad rho|^2 and tau of each spin). xyz: count points, xyz_stride bytes apart (0:
// packed float3), bohr; atom_xyz as for the pattern (where the atoms are now). Any output may be NULL; out_grad holds
// 3 per point. The points spread over executor (NULL: the calling thread). Returns count, 0 on failure.
size_t md_gto_dmat_evaluate(double* out_rho, double* out_grad, double* out_tau, const md_gto_dmat_t* dm, const md_gto_basis_t* basis,
                            const float* atom_xyz, size_t atom_xyz_stride, const float* xyz, size_t count, size_t xyz_stride,
                            const struct md_executor_t* executor);

#ifdef __cplusplus
}
#endif

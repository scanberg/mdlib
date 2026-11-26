// gto_gemm_args.slang
// Argument block shared by every entry point of the two-pass (Phi + GEMM) path: eval_gto_density_gemm.slang (the
// counts, Phi, D_b) and the GEMM passes on the GEMM core (gto_gemm_core.slang): the density's
// (eval_gto_density_core*.slang), ELF's and LOL's (eval_gto_localization_core*.slang), the orbitals'
// (eval_gto_mo_gemm*.slang).
// Mirrored in C by gto_density_gemm_args_t in md_gto.c.

#ifndef GTO_GEMM_ARGS_SLANG
#define GTO_GEMM_ARGS_SLANG

#include "../md_gpu.slang"

#define GTO_SCAN_WG_SIZE 512
#include "gto_common.slang"

#define BLOCK_VOXELS 512

// The shape of a block of voxels (BLOCK_VOXELS of them) is an argument (args.block_shape): 8x8x8, or for a grid of one
// layer (the image of a plane) 32x16x1 - as many voxels as near together, none of the block off the grid. One kernel
// for both: a group's threads in a row (numthreads(BLOCK_TX, 1, 1)), voxel p of a block at block_voxel(p), x
// fastest. Measured against a kernel per shape (groups of 8x8x8 and 32x16x1 threads) on an M4 Pro, an RTX 3080, a
// GTX 1060 and an HD 530 (gemm-bench r6): the density's, ELF's and the orbitals' GEMM passes with their counts and Phi
// within 2.2% on volumes, but the orbitals' GEMM pass on the RTX 3080's protein (1.037, rounds of 0.99-1.13); runs of
// about 1 ms or less (planes, mol's volume on the RTX 3080) within their spread, the median at most 1.15.
#define BLOCK_SHAPE (args.block_shape.xyz)
#define BLOCK_TX BLOCK_VOXELS

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint4    block_dim;         // blocks per axis (xyz)
    uint4    block_shape;       // a block's voxels per axis (xyz), powers of two, and w: log2 x | log2 y << 8
    uint     num_shells;
    uint     num_aos;
    uint     operation;
    uint     num_blocks;        // blocks in this launch (phi/gemm), unused by count
    uint     d_mode;            // 0: D read from its tiles in the GEMM pass, 1: D_b pre-gathered by dgather_main
    uint     num_mos;           // orbital GEMM: number of MO rows in D_matrix (the coefficients)
    uint     mo_mode;           // orbital GEMM: 0 = psi, 1 = psi^2
    uint     has_weights;       // orbital GEMM: mo_weights is valid
    uint     dm_pad;
    uint     dm_nt;             // its tiles per edge
    uint     d_spin_offset;     // floats into D_blocks of the spin being gathered and multiplied (ELF, LOL: two spins)
    uint     spin;              // ELF, LOL: the spin of this GEMM pass (its terms)
    uint     num_spins;         // ELF, LOL: 1, the total density of a closed shell (two equal spins), or 2
    uint     lol;               // ELF, LOL: 0 ELF, 1 LOL; 2 the density and its gradient (out_grad), by the operation
    float    rho_min;           // ELF, LOL: below it (both spins) the value is 0
    float    d_min;             // ELF: added to D (md_qm_elf's d_min)
    uint     terms_stride;      // ELF, LOL: entries of a batch's terms, per spin and term
    uint     terms_valid;       // ELF, LOL finish: 0 for the empty blocks (no terms: written 0)

    uint*    shell_atom_idx;
    float*   shell_r;
    uint4*   shell_info;        // (prim_offset, num_prims, l, ao_offset)
    float2*  shell_prim;        // (coeff, alpha * log2(e))
    float4*  atom_xyz;
    float*   D_matrix;          // density: the tile pool of a md_gto_gpu_dmat (dm_tao, dm_map); orbitals: MO coefficients [num_mos][num_aos]

    uint2*   block_counts;      // count_main output, per grid block: (num shells, num AOs)
    uint4*   block_table;       // 2 per block in batch: (bx, by, bz, n_ao), (phi_offset, ao_offset, D_b offset, -)
    uint*    ao_list;           // batch scratch: global AO index per (block, local AO)
    float*   phi;               // batch scratch
    float*   D_blocks;          // batch scratch (d_mode 1): D_b, n_ao x n_ao per block, lower triangle filled
    float*   mo_weights;        // orbital GEMM: per-MO weight applied after psi / psi^2 (if has_weights)
    uint*    dm_tao;            // density: AO -> tiled index
    uint*    dm_map;            // density: tile pair -> pool entry
    float*   terms;             // ELF, LOL: per spin rho, grad rho (3), tau of the batch's voxels: [spin][5][entry][512]
    DescriptorHandle<RWTexture3D<float>> out_tex;
    DescriptorHandle<RWTexture3D<float4>> out_grad;     // lol 2: the density's gradient (xyz)
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)

uint3 block_voxel(uint p) {
    const uint lx = args.block_shape.w & 0xFFu, ly = (args.block_shape.w >> 8u) & 0xFFu;
    return uint3(p & (args.block_shape.x - 1u), (p >> lx) & (args.block_shape.y - 1u), p >> (lx + ly));
}

// D[i][j] of the density, from its tiles
float gto_d(uint i, uint j) {
    return gto_dmat_at(args.dm_tao, args.dm_map, args.D_matrix, args.dm_nt, i, j);
}

#endif

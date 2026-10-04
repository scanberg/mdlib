// gto_gemm_args.slang
// Argument block shared by every entry point of the two-pass (Phi + GEMM) density
// path: eval_gto_density_gemm.slang and the eval_gto_density_gemm2_*.slang kernels.
// Mirrored in C by gto_density_gemm_args_t in md_gto.c.

#ifndef GTO_GEMM_ARGS_SLANG
#define GTO_GEMM_ARGS_SLANG

#include "../md_gpu.slang"

#define GTO_SCAN_WG_SIZE 512
#include "gto_common.slang"

#define BLOCK_VOXELS 512

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint4    block_dim;         // blocks per axis (xyz)
    uint     num_shells;
    uint     num_aos;
    uint     operation;
    uint     num_blocks;        // blocks in this launch (phi/gemm), unused by count
    uint     d_mode;            // 0: gather D on the fly in gemm_main, 1: D_b pre-gathered by dgather_main
    uint     num_mos;           // orbital GEMM: number of MO rows in D_matrix (the coefficients)
    uint     mo_mode;           // orbital GEMM: 0 = psi, 1 = psi^2
    uint     has_weights;       // orbital GEMM: mo_weights is valid

    uint*    shell_atom_idx;
    float*   shell_r;
    uint4*   shell_info;        // (prim_offset, num_prims, l, ao_offset)
    float2*  shell_prim;        // (coeff, alpha * log2(e))
    float4*  atom_xyz;
    float*   D_matrix;          // density: packed upper triangular D; orbitals: MO coefficients [num_mos][num_aos]

    uint2*   block_counts;      // count_main output, per grid block: (num shells, num AOs)
    uint4*   block_table;       // 2 per block in batch: (bx, by, bz, n_ao), (phi_offset, ao_offset, D_b offset, -)
    uint*    ao_list;           // batch scratch: global AO index per (block, local AO)
    float*   phi;               // batch scratch
    float*   D_blocks;          // batch scratch (d_mode 1): D_b, n_ao x n_ao per block, lower triangle filled
    float*   mo_weights;        // orbital GEMM: per-MO weight applied after psi / psi^2 (if has_weights)
    DescriptorHandle<RWTexture3D<float>> out_tex;
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)

#endif

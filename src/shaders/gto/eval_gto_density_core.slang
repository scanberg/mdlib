// eval_gto_density_core.slang
//
// The density's GEMM pass on the GEMM core (gto_gemm_core.slang): what eval_gto_density_gemm2.slang does in a loop of
// its own, here as the core specialised by the density's operands and epilogue. The same product, in the same order:
//
//   X = Phi_b (GP voxels x n) * (w o D_b) (n x n, symmetric: m <= k),   rho_p = sum_m Phi[p][m] X[p][m]
//
// MD_GTO_GPU_DENSITY_ALGO_GEMM_CORE: for comparing the core with the loop it replaces, kernel for kernel.
//
// Compile-time configuration (the eval_gto_density_core_*.slang wrappers): GP voxels per group, GM columns per tile.

#include "gto_gemm_args.slang"

#ifndef GP
#define GP 128
#endif
#ifndef GM
#define GM 64
#endif

#include "gto_gemm_core.slang"

// A: the block's Phi (AO-major, the group's GP voxels from p0); B: w o D_b, gathered (d_mode 1) or from D
struct DensityOperands : IGemmOperands {
    uint phi_off, ao_off, db_off, n, p0;
    float a(uint row, uint k) {
        return args.phi[phi_off + k * BLOCK_VOXELS + p0 + row];
    }
    float b(uint k, uint m) {
        const float w = (k == m) ? 1.0f : 2.0f;
        if (args.d_mode != 0u) return w * args.D_blocks[db_off + k * n + m];
        return w * gto_d(args.ao_list[ao_off + k], args.ao_list[ao_off + m]);
    }
};

// rho_p += sum_m Phi[p][m] X[p][m] over the thread's rows and columns
struct RhoEpilogue : IGemmEpilogue {
    float rho[RP];
    uint  phi_off, p0;
    [mutating] void take(uint m0, uint n, uint tx, uint ty, float X[RP][RM]) {
        [ForceUnroll] for (int j = 0; j < RM; ++j) {
            uint m = m0 + tx * RM + uint(j);
            if (m < n) {
                [ForceUnroll] for (int i = 0; i < RP; ++i) {
                    rho[i] += args.phi[phi_off + m * BLOCK_VOXELS + p0 + ty * RP + uint(i)] * X[i][j];
                }
            }
        }
    }
};

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint groups_per_block = BLOCK_VOXELS / GP;
    const uint b  = gid.x / groups_per_block;
    const uint pt = gid.x % groups_per_block;
    if (b >= args.num_blocks) return;   // uniform per group

    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint3 blk = e0.xyz;
    const uint  n   = e0.w;

    DensityOperands ops;
    ops.phi_off = e1.x;
    ops.ao_off  = e1.y;
    ops.db_off  = e1.z;
    ops.n       = n;
    ops.p0      = pt * GP;
    RhoEpilogue epi;
    epi.phi_off = e1.x;
    epi.p0      = pt * GP;
    [ForceUnroll] for (int i = 0; i < RP; ++i) epi.rho[i] = 0.0f;

    gemm_core<DensityOperands, RhoEpilogue>(ops, epi, n, true, tid, gtid.x, gtid.y);

    float rho;
    gemm_row_sums(epi.rho, tid, gtid.x, gtid.y, rho);
    if (tid < GP) {
        uint  p   = pt * GP + tid;
        uint3 vox = blk * BLOCK_SHAPE + block_voxel(p);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(rho, prev, args.operation);
        }
    }
}

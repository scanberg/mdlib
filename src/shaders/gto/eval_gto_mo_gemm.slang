// eval_gto_mo_gemm.slang
//
// Many molecular orbitals at once, on the two-pass GEMM path (count + Phi passes of eval_gto_density_gemm.slang,
// then this kernel), on the GEMM core (gto_gemm_core.slang):
//
//   Psi = Phi_b (GP voxels x n AOs) * C_b (n AOs x num_mos orbitals)
//   result_p = sum_m w_m * Psi[p][m]        (mo_mode 0)
//            = sum_m w_m * Psi[p][m]^2      (mo_mode 1)
//
// Every AO is evaluated once per voxel whatever the number of orbitals; the orbitals cost one multiply-add per AO
// each. C_b is read from the coefficient rows by the block's AO list as the core streams B (not symmetric: every k
// of every column).
//
// Compile-time configuration (see the eval_gto_mo_gemm_*.slang wrappers):
//   GP : voxels per group, a multiple of 16 dividing 512
//   GM : orbitals per tile, a multiple of 16

#include "gto_gemm_args.slang"

#ifndef GP
#define GP 128
#endif
#ifndef GM
#define GM 32
#endif

#include "gto_gemm_core.slang"

// A: the block's Phi (AO-major, the group's GP voxels from p0); B: C[m][ao_list[k]]
struct MoOperands : IGemmOperands {
    uint phi_off, ao_off, p0;
    float a(uint row, uint k) {
        return args.phi[phi_off + k * BLOCK_VOXELS + p0 + row];
    }
    float b(uint k, uint m) {
        return args.D_matrix[m * args.num_aos + args.ao_list[ao_off + k]];
    }
};

// Per row: sum_m w_m Psi (mo_mode 0) or w_m Psi^2 (mo_mode 1)
struct MoEpilogue : IGemmEpilogue {
    float res[RP];
    [mutating] void take(uint m0, uint ncols, uint tx, uint ty, float X[RP][RM]) {
        [ForceUnroll] for (int j = 0; j < RM; ++j) {
            const uint m = m0 + tx * RM + uint(j);
            if (m < ncols) {
                float w = 1.0f;
                if (args.has_weights != 0u) w = args.mo_weights[m];
                [ForceUnroll] for (int i = 0; i < RP; ++i) {
                    const float x = X[i][j];
                    res[i] += w * ((args.mo_mode != 0u) ? x * x : x);
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

    MoOperands ops;
    ops.phi_off = e1.x;
    ops.ao_off  = e1.y;
    ops.p0      = pt * GP;
    MoEpilogue epi;
    [ForceUnroll] for (int i = 0; i < RP; ++i) epi.res[i] = 0.0f;

    gemm_core<MoOperands, MoEpilogue>(ops, epi, n, args.num_mos, false, tid, gtid.x, gtid.y);

    float res;
    gemm_row_sums(epi.res, tid, gtid.x, gtid.y, res);
    if (tid < GP) {
        uint  p   = pt * GP + tid;
        uint3 vox = blk * BLOCK_SHAPE + block_voxel(p);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(res, prev, args.operation);
        }
    }
}

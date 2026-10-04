// eval_gto_mo_gemm.slang
//
// Many molecular orbitals at once, on the two-pass GEMM path (count + Phi passes of
// eval_gto_density_gemm.slang, then this kernel):
//
//   Psi = Phi_b (GP voxels x n AOs) * C_b (n AOs x GM orbitals)
//   result_p = sum_m w_m * Psi[p][m]        (mo_mode 0)
//            = sum_m w_m * Psi[p][m]^2      (mo_mode 1)
//
// Every AO is evaluated once per voxel whatever the number of orbitals; the orbitals
// cost one FMA per AO each. C_b is gathered from the coefficient rows by the block's
// AO list while the chunks stream through groupshared memory, with the next chunk
// prefetched into registers, as in eval_gto_density_gemm2.slang.
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

#define GK   16
#define TPX  16
#define TPY  16
#define NT   (TPX * TPY)
#define RP   (GP / TPY)
#define RM   (GM / TPX)
#define PF_P ((GK * GP) / NT)
#define PF_C ((GK * GM) / NT)

groupshared float mo_smem[GK * GP + GK * GM];

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint groups_per_block = BLOCK_VOXELS / GP;
    const uint b  = gid.x / groups_per_block;
    const uint pt = gid.x % groups_per_block;
    if (b >= args.num_blocks) return;   // uniform per group

    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint3 blk     = e0.xyz;
    const uint  n       = e0.w;
    const uint  phi_off = e1.x;
    const uint  ao_off  = e1.y;
    const uint  p0      = pt * GP;
    const uint  tx      = gtid.x;
    const uint  ty      = gtid.y;

    float res_part[RP];
    [ForceUnroll] for (int i = 0; i < RP; ++i) res_part[i] = 0.0f;

    for (uint m0 = 0u; m0 < args.num_mos; m0 += GM) {
        float X[RP][RM];
        [ForceUnroll] for (int i = 0; i < RP; ++i)
            [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] = 0.0f;

        float pf_p[PF_P];
        float pf_c[PF_C];

        #define MO_LOAD_CHUNK(n0_)                                                                  \
        {                                                                                           \
            [ForceUnroll] for (int r = 0; r < PF_P; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint nu  = (n0_) + idx / GP;                                                        \
                float ph = args.phi[phi_off + min(nu, n - 1u) * BLOCK_VOXELS + p0 + idx % GP];      \
                pf_p[r]  = (nu < n) ? ph : 0.0f;                                                    \
            }                                                                                       \
            [ForceUnroll] for (int r = 0; r < PF_C; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint nu  = (n0_) + idx / GM;                                                        \
                uint m   = m0 + idx % GM;                                                           \
                float cv = args.D_matrix[min(m, args.num_mos - 1u) * args.num_aos                   \
                                       + args.ao_list[ao_off + min(nu, n - 1u)]];                    \
                pf_c[r]  = (nu < n && m < args.num_mos) ? cv : 0.0f;                                \
            }                                                                                       \
        }

        if (n > 0u) {
            MO_LOAD_CHUNK(0u);
        }

        for (uint n0 = 0u; n0 < n; n0 += GK) {
            GroupMemoryBarrierWithGroupSync();
            [ForceUnroll] for (int r = 0; r < PF_P; ++r) mo_smem[tid + uint(r) * NT] = pf_p[r];
            [ForceUnroll] for (int r = 0; r < PF_C; ++r) mo_smem[GK * GP + tid + uint(r) * NT] = pf_c[r];
            GroupMemoryBarrierWithGroupSync();

            if (n0 + GK < n) {
                MO_LOAD_CHUNK(n0 + GK);
            }

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) {
                float a[RP];
                float bb[RM];
                [ForceUnroll] for (int i = 0; i < RP; ++i) a[i]  = mo_smem[uint(kk) * GP + ty * RP + uint(i)];
                [ForceUnroll] for (int j = 0; j < RM; ++j) bb[j] = mo_smem[GK * GP + uint(kk) * GM + tx * RM + uint(j)];
                [ForceUnroll] for (int i = 0; i < RP; ++i)
                    [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] += a[i] * bb[j];
            }
        }
        #undef MO_LOAD_CHUNK

        [ForceUnroll] for (int j = 0; j < RM; ++j) {
            const uint m = m0 + tx * RM + uint(j);
            if (m < args.num_mos) {
                float w = 1.0f;
                if (args.has_weights != 0u) w = args.mo_weights[m];
                [ForceUnroll] for (int i = 0; i < RP; ++i) {
                    const float x = X[i][j];
                    res_part[i] += w * ((args.mo_mode != 0u) ? x * x : x);
                }
            }
        }
    }

    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll] for (int i = 0; i < RP; ++i) mo_smem[(ty * RP + uint(i)) * TPX + tx] = res_part[i];
    GroupMemoryBarrierWithGroupSync();

    for (uint q = tid; q < GP; q += NT) {
        float res = 0.0f;
        [ForceUnroll] for (int k = 0; k < TPX; ++k) res += mo_smem[q * TPX + uint(k)];
        uint  p   = p0 + q;
        uint3 vox = blk * uint3(8u, 8u, 8u) + uint3(p & 7u, (p >> 3) & 7u, p >> 6);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(res, prev, args.operation);
        }
    }
}

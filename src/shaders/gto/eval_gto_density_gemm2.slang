// eval_gto_density_gemm2.slang
//
// Configurable fused GEMM + rho pass of the two-pass density path (see
// eval_gto_density_gemm.slang for the Phi pass and the data layout):
//
//   X = Phi_b (GP voxels x n) * (w o D_b) (n x GM),   rho_p += sum_mu Phi[p][mu] X[p][mu]
//
// Compared with gemm_main in eval_gto_density_gemm.slang:
//   * GP voxels per group (64, 128 or 256): every D_b chunk brought into
//     groupshared memory serves GP voxels instead of 64, so the D_b traffic per
//     voxel drops by GP/64, and each thread gets GP/16 x GM/16 accumulators.
//   * GM AOs per output tile (32 or 64): the host sends blocks with few AOs to a
//     32-wide variant, which wastes far less work on padding.
//   * The next chunk of Phi and D_b is loaded into registers while the current one
//     is multiplied, so global-memory latency overlaps the FMAs.
//   * The final cross-thread reduction reuses the Phi buffer, keeping groupshared
//     memory at 4 * 16 * (GP + GM) bytes (<= 20 KB).
//
// Compile-time configuration (see the eval_gto_density_gemm2_*.slang wrappers):
//   GP : voxels per group, a multiple of 16 dividing 512
//   GM : AOs per output tile, a multiple of 16
// The group is always 16 x 16 threads.

#include "gto_gemm_args.slang"

#ifndef GP
#define GP 128
#endif
#ifndef GM
#define GM 64
#endif

#define GK   16                 // nu per chunk
#define TPX  16                 // threads along mu
#define TPY  16                 // threads along p
#define NT   (TPX * TPY)
#define RP   (GP / TPY)         // voxels per thread
#define RM   (GM / TPX)         // mu per thread
#define PF_P ((GK * GP) / NT)   // Phi chunk elements per thread
#define PF_D ((GK * GM) / NT)   // D chunk elements per thread

// [0, GK*GP)            : Phi chunk, [kk][p]
// [GK*GP, GK*(GP+GM))   : w o D chunk, [kk][mm]
// after the tile loop   : [p][TPX] partial rho for the reduction (fits in the Phi part)
groupshared float g2_smem[GK * GP + GK * GM];

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
    const uint  db_off  = e1.z;
    const uint  p0      = pt * GP;
    const uint  tx      = gtid.x;
    const uint  ty      = gtid.y;

    float rho_part[RP];
    [ForceUnroll] for (int i = 0; i < RP; ++i) rho_part[i] = 0.0f;

    for (uint m0 = 0u; m0 < n; m0 += GM) {
        float X[RP][RM];
        [ForceUnroll] for (int i = 0; i < RP; ++i)
            [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] = 0.0f;

        // Registers holding the next chunk.
        float pf_p[PF_P];
        float pf_d[PF_D];

        // Loads chunk n0 into pf_p / pf_d. Consecutive threads read consecutive voxels
        // (Phi) and consecutive mu (D_b), so both streams are coalesced.
        #define G2_LOAD_CHUNK(n0_)                                                                  \
        {                                                                                           \
            [ForceUnroll] for (int r = 0; r < PF_P; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint nu  = (n0_) + idx / GP;                                                        \
                float ph = args.phi[phi_off + min(nu, n - 1u) * BLOCK_VOXELS + p0 + idx % GP];      \
                pf_p[r]  = (nu < n) ? ph : 0.0f;                                                    \
            }                                                                                       \
            [ForceUnroll] for (int r = 0; r < PF_D; ++r) {                                          \
                uint idx = tid + uint(r) * NT;                                                      \
                uint nu  = (n0_) + idx / GM;                                                        \
                uint mu  = m0 + idx % GM;                                                           \
                float v  = 0.0f;                                                                    \
                if (nu < n && mu < n && nu >= mu) {                                                 \
                    float w = (nu == mu) ? 1.0f : 2.0f;                                             \
                    if (args.d_mode != 0u) v = w * args.D_blocks[db_off + nu * n + mu];             \
                    else v = w * args.D_matrix[gto_tri_index(args.ao_list[ao_off + nu], args.ao_list[ao_off + mu], args.num_aos)]; \
                }                                                                                   \
                pf_d[r] = v;                                                                        \
            }                                                                                       \
        }

        G2_LOAD_CHUNK(m0);

        for (uint n0 = m0; n0 < n; n0 += GK) {
            GroupMemoryBarrierWithGroupSync();   // previous chunk consumed
            [ForceUnroll] for (int r = 0; r < PF_P; ++r) g2_smem[tid + uint(r) * NT] = pf_p[r];
            [ForceUnroll] for (int r = 0; r < PF_D; ++r) g2_smem[GK * GP + tid + uint(r) * NT] = pf_d[r];
            GroupMemoryBarrierWithGroupSync();

            if (n0 + GK < n) {
                G2_LOAD_CHUNK(n0 + GK);           // in flight while this chunk is multiplied
            }

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) {
                float a[RP];
                float bb[RM];
                [ForceUnroll] for (int i = 0; i < RP; ++i) a[i]  = g2_smem[uint(kk) * GP + ty * RP + uint(i)];
                [ForceUnroll] for (int j = 0; j < RM; ++j) bb[j] = g2_smem[GK * GP + uint(kk) * GM + tx * RM + uint(j)];
                [ForceUnroll] for (int i = 0; i < RP; ++i)
                    [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] += a[i] * bb[j];
            }
        }
        #undef G2_LOAD_CHUNK

        // Fused epilogue: rho_p += sum_mu Phi[p][mu] * X[p][mu]
        [ForceUnroll] for (int j = 0; j < RM; ++j) {
            uint mu = m0 + tx * RM + uint(j);
            if (mu < n) {
                [ForceUnroll] for (int i = 0; i < RP; ++i) {
                    rho_part[i] += args.phi[phi_off + mu * BLOCK_VOXELS + p0 + ty * RP + uint(i)] * X[i][j];
                }
            }
        }
    }

    // Reduce over the TPX threads that share each voxel.
    GroupMemoryBarrierWithGroupSync();
    [ForceUnroll] for (int i = 0; i < RP; ++i) g2_smem[(ty * RP + uint(i)) * TPX + tx] = rho_part[i];
    GroupMemoryBarrierWithGroupSync();

    for (uint q = tid; q < GP; q += NT) {
        float rho = 0.0f;
        [ForceUnroll] for (int k = 0; k < TPX; ++k) rho += g2_smem[q * TPX + uint(k)];
        uint  p   = p0 + q;
        uint3 vox = blk * uint3(8u, 8u, 8u) + uint3(p & 7u, (p >> 3) & 7u, p >> 6);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(rho, prev, args.operation);
        }
    }
}

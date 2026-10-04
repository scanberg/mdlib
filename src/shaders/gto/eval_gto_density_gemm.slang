// eval_gto_density_gemm.slang
//
// Two-pass electron density evaluation, the way DFT grid integrators do it:
//
//   rho_p = sum_mu Phi[p][mu] * X[p][mu],   X = Phi * D_b
//
// for every 8x8x8 block b of the grid, where Phi (512 voxels x n_b AOs) holds the
// AOs that survive screening for that block and D_b is the matching sub-matrix of D.
//
// Entry points (all share one argument struct):
//   count_main : one group per block of the grid. Counts screened shells / AOs.
//                The host reads the counts back and plans batches of blocks whose
//                Phi fits in the scratch budget.
//   phi_main   : one group per block in a batch. Re-screens (same code, same
//                result), writes the block's AO index list and evaluates Phi into
//                scratch, AO-major: Phi[ao * 512 + voxel]. Each AO is evaluated
//                exactly once per voxel, and the radial part once per shell.
//                A prefix of the AO list wins if the re-screen ever disagrees.
//   dgather_main: (optional) one group per block in a batch, gathers D_b once.
//   gemm_main  : 8 groups per block in a batch, 64 voxels each. Register-tiled
//                FP32 GEMM X = Phi * (w o D_b) over the upper triangle (w = 2 above
//                the diagonal, 1 on it), with the rho epilogue fused so X never
//                touches memory. Also handles empty blocks (n_b = 0).

#include "gto_gemm_args.slang"

// Screening test, shared by count_main and phi_main so both see identical results.
bool shell_passes(uint s_in, float3 aabb_min, float3 aabb_max) {
    // Clamped so that no lane ever forms an out-of-range address, even if the loads of
    // inactive lanes are evaluated (seen with llvmpipe); callers check s < num_shells.
    const uint s = min(s_in, args.num_shells - 1u);
    float r = args.shell_r[s];
    if (r <= 0.0f) return false;
    float3 c = args.atom_xyz[args.shell_atom_idx[s]].xyz;
    float3 m = mul(args.world_to_model, float4(c, 1.0f)).xyz;
    return gto_aabb_dist2(m, aabb_min, aabb_max) < r * r;
}

void block_aabb(uint3 blk, out float3 aabb_min, out float3 aabb_max) {
    aabb_min = float3( blk             * uint3(8u, 8u, 8u)) * args.step.xyz;
    aabb_max = float3((blk + uint3(1u)) * uint3(8u, 8u, 8u)) * args.step.xyz;
}

// ---------------------------------------------------------------------------
// Pass 0: counts
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 8)]
void count_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    float3 aabb_min, aabb_max;
    block_aabb(gid, aabb_min, aabb_max);

    uint2 sum = uint2(0u, 0u);
    for (uint base = 0u; base < args.num_shells; base += 512u) {
        uint s = base + tid;
        if (s < args.num_shells && shell_passes(s, aabb_min, aabb_max)) {
            sum += uint2(1u, gto_num_cart(args.shell_info[s].z));
        }
    }
    uint2 total;
    gto_wg_exclusive_scan(sum, tid, total);
    if (tid == 0u) {
        uint b = (gid.z * args.block_dim.y + gid.y) * args.block_dim.x + gid.x;
        args.block_counts[b] = total;
    }
}

// ---------------------------------------------------------------------------
// Pass 1: AO list + Phi
// ---------------------------------------------------------------------------
// Shells are screened in chunks of 512 in ascending order; each chunk's passing
// shells go to groupshared memory and are evaluated straight away, so there is no
// limit on the number of shells per block and no device-memory round trip inside
// the group.
#define INVALID_SHELL 0xFFFFFFFFu
groupshared uint2 chunk_list[512];      // (global shell, local AO offset) of this chunk's passing shells
groupshared uint  phi_na_written;

[shader("compute")]
[numthreads(8, 8, 8)]
void phi_main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint b = gid.x;
    if (b >= args.num_blocks) return;   // uniform per group

    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint3 blk     = e0.xyz;
    const uint  cap_ao  = e0.w;
    const uint  phi_off = e1.x;
    const uint  ao_off  = e1.y;

    float3 aabb_min, aabb_max;
    block_aabb(blk, aabb_min, aabb_max);

    const uint3  vox   = blk * uint3(8u, 8u, 8u) + gtid;
    const float3 coord = mul(args.index_to_world, float4(float3(vox), 1.0f)).xyz;

    if (tid == 0u) phi_na_written = 0u;

    uint2 carry = uint2(0u, 0u);
    for (uint base = 0u; base < args.num_shells; base += 512u) {
        uint s = base + tid;
        bool pass = false;
        uint ncart = 0u;
        uint4 info = uint4(0u);
        if (s < args.num_shells && shell_passes(s, aabb_min, aabb_max)) {
            pass  = true;
            info  = args.shell_info[min(s, args.num_shells - 1u)];
            ncart = gto_num_cart(info.z);
        }
        uint2 total;
        uint2 excl = gto_wg_exclusive_scan(uint2(pass ? 1u : 0u, ncart), tid, total);
        if (pass) {
            uint aoff = carry.y + excl.y;
            // Never beyond the capacity planned from count_main (whole shells only).
            if (aoff + ncart <= cap_ao) {
                chunk_list[excl.x] = uint2(s, aoff);
                for (uint c = 0u; c < ncart; ++c) {
                    args.ao_list[ao_off + aoff + c] = info.w + c;
                }
                InterlockedMax(phi_na_written, aoff + ncart);
            } else {
                chunk_list[excl.x] = uint2(INVALID_SHELL, 0u);
            }
        }
        GroupMemoryBarrierWithGroupSync();

        for (uint k = 0u; k < total.x; ++k) {
            uint2 sl = chunk_list[k];
            if (sl.x == INVALID_SHELL) continue;
            const uint sh = min(sl.x, args.num_shells - 1u);
            uint4 si = args.shell_info[sh];
            float3 d  = coord - args.atom_xyz[args.shell_atom_idx[sh]].xyz;
            float r2  = dot(d, d);
            float R   = 0.0f;
            for (uint p = 0u; p < si.y; ++p) {
                float2 cp = args.shell_prim[si.x + p];
                R += cp.x * exp2(-cp.y * r2);
            }
            uint nc   = gto_num_cart(si.z);
            uint dst  = phi_off + sl.y * BLOCK_VOXELS + tid;
            for (uint c = 0u; c < nc; ++c) {
                uint ijk = gto_cart_ijk_packed(si.z, c);
                args.phi[dst + c * BLOCK_VOXELS] = gto_cart_norm(ijk) * gto_monomial(d, ijk) * R;
            }
        }
        GroupMemoryBarrierWithGroupSync();   // chunk_list is rewritten by the next chunk
        carry += total;
    }

    // Normally phi_na_written == cap_ao. Should the re-screen ever disagree with
    // count_main, publish what was actually written so gemm_main reads no stale rows.
    if (tid == 0u && phi_na_written != cap_ao) {
        args.block_table[2u * b].w = phi_na_written;
    }
}

// ---------------------------------------------------------------------------
// Pass 1b (d_mode 1): gather D_b once per block
// ---------------------------------------------------------------------------
// D_b[i][j] = D[ao_list[i]][ao_list[j]] for i >= j (the only part gemm_main reads),
// so that gemm_main streams it with coalesced loads instead of gathering from the
// packed global triangle in each of its 8 groups per block.
[shader("compute")]
[numthreads(256, 1, 1)]
void dgather_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    const uint b = gid.x;
    if (b >= args.num_blocks) return;
    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint n      = e0.w;
    const uint ao_off = e1.y;
    const uint db_off = e1.z;
    for (uint i = 0u; i < n; ++i) {             // uniform loop over rows
        const uint gi = args.ao_list[ao_off + i];
        for (uint j = tid; j <= i; j += 256u) {
            args.D_blocks[db_off + i * n + j] = args.D_matrix[gto_tri_index(gi, args.ao_list[ao_off + j], args.num_aos)];
        }
    }
}

// ---------------------------------------------------------------------------
// Pass 2: fused GEMM + rho
// ---------------------------------------------------------------------------
#define GP 64      // voxels per group
#define GM 64      // mu per tile
#define GK 16      // nu per chunk
#define TPX 16     // threads along mu
#define TPY 16     // threads along p
#define RP 4       // voxels per thread
#define RM 4       // mu per thread

groupshared float Ps[GK][GP];       // Phi chunk  [nu][p]
groupshared float Ds[GK][GM];       // w o D chunk [nu][mu]
groupshared uint  mu_gao[GM];
groupshared float red[GP][TPX + 1];

[shader("compute")]
[numthreads(16, 16, 1)]
void gemm_main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint b  = gid.x / (BLOCK_VOXELS / GP);
    const uint pt = gid.x % (BLOCK_VOXELS / GP);
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
        if (tid < GM) {
            mu_gao[tid] = args.ao_list[ao_off + min(m0 + tid, n - 1u)];   // padding columns are zeroed via D
        }

        float X[RP][RM];
        [ForceUnroll] for (int i = 0; i < RP; ++i)
            [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] = 0.0f;

        for (uint n0 = m0; n0 < n; n0 += GK) {
            GroupMemoryBarrierWithGroupSync();   // previous chunk consumed; mu_gao visible
            for (uint idx = tid; idx < GK * GP; idx += 256u) {
                uint kk = idx / GP;
                uint p  = idx % GP;
                uint nu = n0 + kk;
                const float ph = args.phi[phi_off + min(nu, n - 1u) * BLOCK_VOXELS + p0 + p];
                Ps[kk][p] = (nu < n) ? ph : 0.0f;
            }
            for (uint idx = tid; idx < GK * GM; idx += 256u) {
                uint kk = idx / GM;
                uint mm = idx % GM;
                uint nu = n0 + kk;
                uint mu = m0 + mm;
                float v = 0.0f;
                if (nu < n && mu < n && nu >= mu) {
                    float w = (nu == mu) ? 1.0f : 2.0f;
                    if (args.d_mode != 0u) v = w * args.D_blocks[db_off + nu * n + mu];
                    else                   v = w * args.D_matrix[gto_tri_index(args.ao_list[ao_off + nu], mu_gao[mm], args.num_aos)];
                }
                Ds[kk][mm] = v;
            }
            GroupMemoryBarrierWithGroupSync();

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) {
                float a[RP];
                float bb[RM];
                [ForceUnroll] for (int i = 0; i < RP; ++i) a[i]  = Ps[kk][ty * RP + i];
                [ForceUnroll] for (int j = 0; j < RM; ++j) bb[j] = Ds[kk][tx * RM + j];
                [ForceUnroll] for (int i = 0; i < RP; ++i)
                    [ForceUnroll] for (int j = 0; j < RM; ++j) X[i][j] += a[i] * bb[j];
            }
        }

        // Fused epilogue: rho_p += sum_mu Phi[p][mu] * X[p][mu]
        [ForceUnroll] for (int j = 0; j < RM; ++j) {
            uint mu = m0 + tx * RM + j;
            if (mu < n) {
                [ForceUnroll] for (int i = 0; i < RP; ++i) {
                    rho_part[i] += args.phi[phi_off + mu * BLOCK_VOXELS + p0 + ty * RP + i] * X[i][j];
                }
            }
        }
        GroupMemoryBarrierWithGroupSync();   // mu_gao is rewritten by the next tile
    }

    [ForceUnroll] for (int i = 0; i < RP; ++i) red[ty * RP + i][tx] = rho_part[i];
    GroupMemoryBarrierWithGroupSync();

    if (tid < GP) {
        float rho = 0.0f;
        for (uint k = 0u; k < TPX; ++k) rho += red[tid][k];
        uint  p   = p0 + tid;
        uint3 vox = blk * uint3(8u, 8u, 8u) + uint3(p & 7u, (p >> 3) & 7u, p >> 6);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(rho, prev, args.operation);
        }
    }
}

// ---------------------------------------------------------------------------
// Pass 2, matrix-unit version: the same fused GEMM + rho on 8x8 cooperative
// matrices (Slang linalg.CoopMat). On Metal these are simdgroup_matrix<float>
// with FP32 accumulation; the SPIR-V form needs VK_KHR_cooperative_matrix, so
// the host only uses this entry point on the Metal backend.
//
// 128 threads = 4 subgroups of 32 per group. The group computes a 64 (voxels)
// x 64 (mu) tile of X; each subgroup a 32 x 32 quarter as 4 x 4 accumulators.
// nu is streamed in chunks of 32 through groupshared memory, which is reused
// for X^T once the chunk loop is done so the epilogue reads it conflict free.
// ---------------------------------------------------------------------------
typealias GtoMatA = linalg.CoopMat<float, MemoryScope.Subgroup, 8, 8, linalg.CoopMatMatrixUse.MatrixA>;
typealias GtoMatB = linalg.CoopMat<float, MemoryScope.Subgroup, 8, 8, linalg.CoopMatMatrixUse.MatrixB>;
typealias GtoMatC = linalg.CoopMat<float, MemoryScope.Subgroup, 8, 8, linalg.CoopMatMatrixUse.MatrixAccumulator>;

#define SG_P   64          // voxels per group
#define SG_M   64          // mu per tile
#define SG_K   32          // nu per chunk
#define SG_THREADS 128

// [0, SG_K*SG_P)            : Phi chunk, element (p, kk) at kk*SG_P + p    (A, column major)
// [SG_K*SG_P, 2*SG_K*SG_P)  : w o D chunk, element (kk, mm) at kk*SG_M + mm (B, row major)
// after the chunk loop, all : X^T, element (p, mu) at mu*SG_P + p        (C, column major)
groupshared float sg_smem[2 * SG_K * SG_P];
groupshared uint  sg_mu_gao[SG_M];
groupshared float sg_red[2][SG_P];

[shader("compute")]
[numthreads(128, 1, 1)]
void gemm_sg_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    const uint b  = gid.x / (BLOCK_VOXELS / SG_P);
    const uint pt = gid.x % (BLOCK_VOXELS / SG_P);
    if (b >= args.num_blocks) return;   // uniform per group

    uint4 e0 = args.block_table[2u * b + 0u];
    uint4 e1 = args.block_table[2u * b + 1u];
    const uint3 blk     = e0.xyz;
    const uint  n       = e0.w;
    const uint  phi_off = e1.x;
    const uint  ao_off  = e1.y;
    const uint  db_off  = e1.z;
    const uint  p0      = pt * SG_P;

    const uint sg = tid / 32u;          // subgroup: (p half, mu half)
    const uint sp = (sg & 1u) * 32u;
    const uint sm = (sg >> 1) * 32u;

    // Epilogue mapping: voxel p, mu half q.
    const uint ep = tid & (SG_P - 1u);
    const uint eq = tid / SG_P;
    float rho = 0.0f;

    for (uint m0 = 0u; m0 < n; m0 += SG_M) {
        if (tid < SG_M) {
            sg_mu_gao[tid] = args.ao_list[ao_off + min(m0 + tid, n - 1u)];   // padding columns are zeroed via D
        }

        GtoMatC acc[4][4];
        [ForceUnroll] for (int i = 0; i < 4; ++i)
            [ForceUnroll] for (int j = 0; j < 4; ++j) acc[i][j] = GtoMatC(0.0f);

        for (uint n0 = m0; n0 < n; n0 += SG_K) {
            GroupMemoryBarrierWithGroupSync();   // previous chunk / epilogue done with sg_smem; mu_gao visible
            for (uint idx = tid; idx < SG_K * SG_P; idx += SG_THREADS) {
                uint kk = idx / SG_P;
                uint p  = idx % SG_P;
                uint nu = n0 + kk;
                const float ph = args.phi[phi_off + min(nu, n - 1u) * BLOCK_VOXELS + p0 + p];
                sg_smem[idx] = (nu < n) ? ph : 0.0f;
            }
            for (uint idx = tid; idx < SG_K * SG_M; idx += SG_THREADS) {
                uint kk = idx / SG_M;
                uint mm = idx % SG_M;
                uint nu = n0 + kk;
                uint mu = m0 + mm;
                float v = 0.0f;
                if (nu < n && mu < n && nu >= mu) {
                    float w = (nu == mu) ? 1.0f : 2.0f;
                    if (args.d_mode != 0u) v = w * args.D_blocks[db_off + nu * n + mu];
                    else                   v = w * args.D_matrix[gto_tri_index(args.ao_list[ao_off + nu], sg_mu_gao[mm], args.num_aos)];
                }
                sg_smem[SG_K * SG_P + idx] = v;
            }
            GroupMemoryBarrierWithGroupSync();

            [ForceUnroll]
            for (int k = 0; k < SG_K; k += 8) {
                GtoMatA a[4];
                GtoMatB bm[4];
                [ForceUnroll] for (int i = 0; i < 4; ++i)
                    a[i] = GtoMatA.Load<linalg.CoopMatMatrixLayout.ColumnMajor>(sg_smem, uint(k) * SG_P + sp + uint(i) * 8u, SG_P);
                [ForceUnroll] for (int j = 0; j < 4; ++j)
                    bm[j] = GtoMatB.Load<linalg.CoopMatMatrixLayout.RowMajor>(sg_smem, SG_K * SG_P + uint(k) * SG_M + sm + uint(j) * 8u, SG_M);
                [ForceUnroll] for (int i = 0; i < 4; ++i)
                    [ForceUnroll] for (int j = 0; j < 4; ++j)
                        acc[i][j] = linalg.coopMatMulAdd<float, false>(a[i], bm[j], acc[i][j]);
            }
        }

        // X^T into groupshared memory, then rho_p += sum_mu Phi[p][mu] X[p][mu].
        GroupMemoryBarrierWithGroupSync();   // every subgroup is done reading the last chunk
        [ForceUnroll] for (int i = 0; i < 4; ++i)
            [ForceUnroll] for (int j = 0; j < 4; ++j)
                acc[i][j].Store<linalg.CoopMatMatrixLayout.ColumnMajor>(sg_smem, (sm + uint(j) * 8u) * SG_P + sp + uint(i) * 8u, SG_P);
        GroupMemoryBarrierWithGroupSync();

        const uint mu_end = min(n - m0, SG_M);
        for (uint ml = eq * (SG_M / 2u); ml < min(mu_end, (eq + 1u) * (SG_M / 2u)); ++ml) {
            rho += sg_smem[ml * SG_P + ep] * args.phi[phi_off + (m0 + ml) * BLOCK_VOXELS + p0 + ep];
        }
    }

    sg_red[eq][ep] = rho;
    GroupMemoryBarrierWithGroupSync();
    if (tid < SG_P) {
        float r   = sg_red[0][tid] + sg_red[1][tid];
        uint  p   = p0 + tid;
        uint3 vox = blk * uint3(8u, 8u, 8u) + uint3(p & 7u, (p >> 3) & 7u, p >> 6);
        if (all(vox < args.grid_dim.xyz)) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox] : 0.0f;
            args.out_tex[vox] = gto_apply_op(r, prev, args.operation);
        }
    }
}

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
//   dgather_main: one group per block in a batch, gathers D_b once.
//
// The GEMM pass itself, X = Phi * (w o D_b) with the rho epilogue fused, is in
// eval_gto_density_gemm2.slang (one build per tile configuration).

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
    // count_main, publish what was actually written so the GEMM pass reads no stale rows.
    if (tid == 0u && phi_na_written != cap_ao) {
        args.block_table[2u * b].w = phi_na_written;
    }
}

// ---------------------------------------------------------------------------
// Pass 1b (d_mode 1): gather D_b once per block
// ---------------------------------------------------------------------------
// D_b[i][j] = D[ao_list[i]][ao_list[j]] for i >= j (the only part the GEMM pass reads),
// so that the GEMM pass streams it with coalesced loads instead of gathering from the
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

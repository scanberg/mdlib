// eval_gto_density_tiled.slang
//
// Electron density rho(r) = sum_{mu,nu} D_{mu nu} phi_mu(r) phi_nu(r) on a 3D grid,
// one workgroup per 8x8x8 block of voxels.
//
// Differences from eval_gto_density.slang (the per-AO reference kernel):
//   * Shell based. Screening is per contracted shell, and the radial part
//     R(r) = sum_p c_p exp(-a_p r^2) is evaluated once per shell and shared by all
//     (l+1)(l+2)/2 Cartesian components, instead of once per component.
//   * Asymmetric tiles. The row ("mu") tile of TM AOs is evaluated once and kept in
//     registers while the column ("nu") tiles of TN AOs stream past it, so each AO
//     is re-evaluated ~N/(2*TM) times instead of ~N/(2*16).
//   * VPT voxels per thread. Every D-tile load from groupshared memory and every
//     staged primitive is reused for VPT voxels, and screening/tile setup is
//     amortised over more work per thread.
//   * Column tiles are double buffered, so a tile pair costs 2 barriers.
//   * Primitives beyond the groupshared staging capacity are read from device
//     memory instead of being dropped.
//
// Compile-time configuration (see the eval_gto_density_tiled_*.slang wrappers):
//   TM  : AOs per row tile      (registers: VPT * TM)
//   TN  : AOs per column tile   (registers: VPT * TN)
//   VPT : voxels per thread (1, 2 or 4); the workgroup is 8 x 8 x (8/VPT) threads

#include "../md_gpu.slang"

#ifndef TM
#define TM 32
#endif
#ifndef TN
#define TN 16
#endif
#ifndef VPT
#define VPT 1
#endif

#define WG_X 8
#define WG_Y 8
#define WG_Z (8 / VPT)
#define WG_SIZE (WG_X * WG_Y * WG_Z)

#define GTO_SCAN_WG_SIZE WG_SIZE
#include "gto_common.slang"

#define MAX_SS     1024        // screened shells per block kept in groupshared memory
#define MAXP       8           // primitives per shell staged in groupshared memory
#define TMAX       ((TM > TN) ? TM : TN)
#define FIRST_BIT  0x80000000u
#define GLOBAL_BIT 0x80000000u

struct RootArgs {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint4    grid_dim;
    uint     num_shells;
    uint     num_aos;
    uint     operation;
    uint     _pad0;

    uint*    shell_atom_idx;
    float*   shell_r;
    uint4*   shell_info;      // (prim_offset, num_prims, l, ao_offset)
    float2*  shell_prim;      // (coeff, alpha * log2(e))
    float4*  atom_xyz;
    float*   D_matrix;        // packed upper triangular, num_aos x num_aos
    DescriptorHandle<RWTexture3D<float>> out_tex;
};

struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;

#define args (*md_root.args_ptr)

// ---- Screened shells of this block -----------------------------------------
groupshared uint ss_gid [MAX_SS];       // screened slot -> global shell index (ascending)
groupshared uint ss_aoff[MAX_SS + 1];   // screened slot -> first screened AO; [count] = total
groupshared uint ss_count;
groupshared uint ss_naos;

// ---- AO tiles: slot 0 = row tile, slots 1,2 = double-buffered column tiles ---
groupshared uint   tl_meta  [3][TMAX];  // ijk (12 bits) | tile shell (8 bits) << 16 | FIRST_BIT
groupshared float  tl_norm  [3][TMAX];  // Cartesian normalisation, 0 for padding
groupshared uint   tl_gao   [3][TMAX];  // global AO index
groupshared float4 tl_sh_xyz[3][TMAX];  // per tile shell: centre
groupshared uint2  tl_sh_prim[3][TMAX]; // per tile shell: (offset, count | GLOBAL_BIT)
groupshared uint   tl_sh_src[3][TMAX];  // per tile shell: global primitive offset
groupshared float2 tl_prim  [3][TMAX * MAXP];

groupshared float  D_tile[2][TM][TN];   // weighted D sub-block for (row tile, column tile)

void screen_shells(uint3 gid, uint tid) {
    float3 aabb_min = float3( gid             * uint3(8u, 8u, 8u)) * args.step.xyz;
    float3 aabb_max = float3((gid + uint3(1u)) * uint3(8u, 8u, 8u)) * args.step.xyz;

    uint2 carry = uint2(0u, 0u);   // (shells, AOs) before this chunk
    for (uint base = 0u; base < args.num_shells; base += WG_SIZE) {
        uint s_raw = base + tid;
        uint s = min(s_raw, args.num_shells - 1u);   // clamped: see eval_gto_mo_shell.slang
        bool pass = false;
        uint ncart = 0u;
        if (s_raw < args.num_shells) {
            float r = args.shell_r[s];
            if (r > 0.0f) {
                float3 c = args.atom_xyz[args.shell_atom_idx[s]].xyz;
                float3 m = mul(args.world_to_model, float4(c, 1.0f)).xyz;
                pass = gto_aabb_dist2(m, aabb_min, aabb_max) < r * r;
            }
            if (pass) ncart = gto_num_cart(args.shell_info[s].z);
        }
        uint2 total;
        uint2 excl = gto_wg_exclusive_scan(uint2(pass ? 1u : 0u, ncart), tid, total);
        if (pass) {
            uint slot = carry.x + excl.x;
            if (slot < MAX_SS) {
                ss_gid[slot]  = s;
                ss_aoff[slot] = carry.y + excl.y;
            }
        }
        carry += total;
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        uint n = min(carry.x, (uint)MAX_SS);
        uint naos = carry.y;
        if (carry.x > MAX_SS) {
            // Overflow: keep the first MAX_SS shells (whole shells only).
            naos = ss_aoff[MAX_SS - 1] + gto_num_cart(args.shell_info[ss_gid[MAX_SS - 1]].z);
        }
        ss_aoff[n] = naos;
        ss_count   = n;
        ss_naos    = naos;
    }
    GroupMemoryBarrierWithGroupSync();
}

// Screened slot containing screened AO k (k < ss_naos).
uint find_ss(uint k) {
    uint lo = 0u;
    uint hi = ss_count;
    while (hi - lo > 1u) {
        uint mid = (lo + hi) >> 1;
        if (ss_aoff[mid] <= k) lo = mid; else hi = mid;
    }
    return lo;
}

// Phase 1 of a tile load: threads [0, T) fill the per-AO and per-shell headers.
void tile_load_meta(uint slot, uint k0, uint T, uint tid) {
    if (tid < T) {
        uint k    = k0 + tid;
        uint meta = 0u;
        float nrm = 0.0f;
        uint gao  = 0u;
        if (k < ss_naos) {
            uint ss0 = find_ss(k0);
            uint ss  = find_ss(k);
            uint c   = k - ss_aoff[ss];
            uint s   = min(ss_gid[ss], args.num_shells - 1u);   // clamped: see eval_gto_mo_shell.slang
            uint4 info = args.shell_info[s];
            uint ijk = gto_cart_ijk_packed(info.z, c);
            nrm = gto_cart_norm(ijk);
            gao = info.w + c;
            uint ts = ss - ss0;
            bool first = (tid == 0u) || (c == 0u);
            meta = ijk | (ts << 16) | (first ? FIRST_BIT : 0u);
            if (first) {
                uint np = info.y;
                tl_sh_xyz [slot][ts] = float4(args.atom_xyz[args.shell_atom_idx[s]].xyz, 0.0f);
                tl_sh_prim[slot][ts] = (np <= MAXP) ? uint2(ts * MAXP, np) : uint2(info.x, np | GLOBAL_BIT);
                tl_sh_src [slot][ts] = info.x;
            }
        }
        tl_meta[slot][tid] = meta;
        tl_norm[slot][tid] = nrm;
        tl_gao [slot][tid] = gao;
    }
}

// Phase 2 of a tile load: stage primitives of the tile's shells (all threads).
void tile_stage_prims(uint slot, uint T, uint tid) {
    for (uint idx = tid; idx < T * MAXP; idx += WG_SIZE) {
        uint a = idx / MAXP;
        uint p = idx % MAXP;
        uint meta = tl_meta[slot][a];
        if ((meta & FIRST_BIT) != 0u) {
            uint ts = (meta >> 16) & 0xFFu;
            uint2 pr = tl_sh_prim[slot][ts];
            if ((pr.y & GLOBAL_BIT) == 0u && p < pr.y) {
                tl_prim[slot][pr.x + p] = args.shell_prim[tl_sh_src[slot][ts] + p];
            }
        }
    }
}

// Weighted D block: w * D[mu][nu] with w = 2 above the diagonal, 1 on it, 0 below,
// so that summing the full rectangular tiles of the upper triangle gives rho.
void gather_D(uint buf, uint I, uint J, uint N, uint tid) {
    for (uint idx = tid; idx < TM * TN; idx += WG_SIZE) {
        uint m  = idx / TN;
        uint n  = idx % TN;
        uint mu = I * TM + m;
        uint nu = J * TN + n;
        float val = 0.0f;
        if (mu < N && nu < N && nu >= mu) {
            float w = (nu == mu) ? 1.0f : 2.0f;
            val = w * args.D_matrix[gto_tri_index(tl_gao[0][m], tl_gao[1 + buf][n], args.num_aos)];
        }
        D_tile[buf][m][n] = val;
    }
}

void eval_tile<let T : int>(uint slot, float3 coord[VPT], out float phi[VPT][T]) {
    float3 d[VPT];
    float  R[VPT];
    [ForceUnroll] for (int v = 0; v < VPT; ++v) { d[v] = float3(0.0f); R[v] = 0.0f; }

    [ForceUnroll]
    for (int a = 0; a < T; ++a) {
        uint meta = tl_meta[slot][a];
        if ((meta & FIRST_BIT) != 0u) {   // uniform across the workgroup
            uint   ts = (meta >> 16) & 0xFFu;
            float3 c  = tl_sh_xyz[slot][ts].xyz;
            uint2  pr = tl_sh_prim[slot][ts];
            uint   np = pr.y & ~GLOBAL_BIT;
            float  r2[VPT];
            [ForceUnroll] for (int v = 0; v < VPT; ++v) {
                d[v]  = coord[v] - c;
                r2[v] = dot(d[v], d[v]);
                R[v]  = 0.0f;
            }
            if ((pr.y & GLOBAL_BIT) == 0u) {
                for (uint p = 0u; p < np; ++p) {
                    float2 cp = tl_prim[slot][pr.x + p];
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) R[v] += cp.x * exp2(-cp.y * r2[v]);
                }
            } else {
                for (uint p = 0u; p < np; ++p) {
                    float2 cp = args.shell_prim[pr.x + p];
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) R[v] += cp.x * exp2(-cp.y * r2[v]);
                }
            }
        }
        float nrm = tl_norm[slot][a];
        uint  ijk = meta & 0xFFFu;
        [ForceUnroll] for (int v = 0; v < VPT; ++v) {
            phi[v][a] = nrm * gto_monomial(d[v], ijk) * R[v];
        }
    }
}

[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    screen_shells(gid, tid);

    float3 coord[VPT];
    bool   active[VPT];
    uint3  vox[VPT];
    float  acc[VPT];
    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
        vox[v]    = gid * uint3(8u, 8u, 8u) + uint3(gtid.x, gtid.y, gtid.z + uint(v) * WG_Z);
        active[v] = all(vox[v] < args.grid_dim.xyz);
        coord[v]  = mul(args.index_to_world, float4(float3(vox[v]), 1.0f)).xyz;
        acc[v]    = 0.0f;
    }

    const uint N = ss_naos;
    const uint num_row_tiles = (N + TM - 1u) / TM;
    const uint num_col_tiles = (N + TN - 1u) / TN;
    uint buf = 0u;

    for (uint I = 0u; I < num_row_tiles; ++I) {
        tile_load_meta(0u, I * TM, TM, tid);
        GroupMemoryBarrierWithGroupSync();
        tile_stage_prims(0u, TM, tid);
        GroupMemoryBarrierWithGroupSync();

        float phi_m[VPT][TM];
        eval_tile<TM>(0u, coord, phi_m);

        for (uint J = (I * TM) / TN; J < num_col_tiles; ++J) {
            uint slot = 1u + buf;
            tile_load_meta(slot, J * TN, TN, tid);
            GroupMemoryBarrierWithGroupSync();
            tile_stage_prims(slot, TN, tid);
            gather_D(buf, I, J, N, tid);
            GroupMemoryBarrierWithGroupSync();

            float phi_n[VPT][TN];
            eval_tile<TN>(slot, coord, phi_n);

            [ForceUnroll]
            for (int m = 0; m < TM; ++m) {
                float t[VPT];
                [ForceUnroll] for (int v = 0; v < VPT; ++v) t[v] = 0.0f;
                [ForceUnroll]
                for (int n = 0; n < TN; ++n) {
                    float dmn = D_tile[buf][m][n];
                    [ForceUnroll] for (int v = 0; v < VPT; ++v) t[v] += dmn * phi_n[v][n];
                }
                [ForceUnroll] for (int v = 0; v < VPT; ++v) acc[v] += phi_m[v][m] * t[v];
            }
            buf ^= 1u;
        }
        // Slot 0 and both column buffers are rewritten by the next row tile.
        GroupMemoryBarrierWithGroupSync();
    }

    [ForceUnroll] for (int v = 0; v < VPT; ++v) {
        if (active[v]) {
            float prev = gto_op_reads_prev(args.operation) ? args.out_tex[vox[v]] : 0.0f;
            args.out_tex[vox[v]] = gto_apply_op(acc[v], prev, args.operation);
        }
    }
}

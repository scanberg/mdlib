// eval_gto_mo.slang
// Evaluates molecular orbital wavefunction(s) on a 3D grid.
//
//   PSI mode:         result[voxel] = sum_m( sum_c( C[m][c] * phi_c(voxel) ) )
//   PSI_SQUARED mode: result[voxel] = sum_m( (sum_c C[m][c] * phi_c(voxel))^2 )
//
// PSI mode collapses the MO sum into a single effective coefficient per screened CGTO
// (effective_coeffs[c] = sum_m C[m][c]) so only one tile pass is needed regardless of
// how many MO rows are present.
//
// PSI_SQUARED mode loops over MO rows in the outer loop and CGTO tiles in the inner
// loop.  PGTO data is reloaded for each MO row; this is the correct trade-off for the
// typical use case of 1–4 orbitals.  For many-orbital density accumulation, use the
// density shader instead.
//
// Buffer bindings:
//   0  cgto_xyzr    float4[num_cgtos]         (xyz = center, w = cutoff radius)
//   1  cgto_off_len uint2[num_cgtos]           (offset, length into pgto buffer)
//   2  pgto         PGTO[num_pgtos]            (coeff, alpha, ijkl)
//   3  coeffs       float[num_rows * num_cgtos] row-major MO coefficient matrix
//  16  out_vol      RWTexture3D<float>          output grid

struct PGTO {
    float coeff;
    float alpha;
    uint  ijkl;
};

struct UBO {
    float4x4 world_to_model;  // world → model space (rotation + translation, no scale)
    float4x4 index_to_world;  // grid index → world
    float4   step;            // xyz = voxel spacing in model-space units
    uint     num_cgtos;       // total number of CGTOs in the basis
    uint     num_rows;        // number of MO coefficient rows
    uint     eval_mode;       // 0 = PSI, 1 = PSI_SQUARED
    uint     _pad;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]] StructuredBuffer<float4> cgto_xyzr;
[[vk::binding(1)]] StructuredBuffer<uint2>  cgto_off_len;
[[vk::binding(2)]] StructuredBuffer<PGTO>   pgto;
[[vk::binding(3)]] StructuredBuffer<float>  coeffs;   // float[num_rows][num_cgtos], row-major

[[vk::binding(16)]] // image slot 0  (slot + MD_GPU_MAX_BIND_SLOTS)
RWTexture3D<float> out_vol;

uint4 unpack_ijkl(uint packed) {
    return uint4(
        (packed      ) & 0xFFu,
        (packed >>  8) & 0xFFu,
        (packed >> 16) & 0xFFu,
        (packed >> 24) & 0xFFu
    );
}

#define WG_X 8
#define WG_Y 8
#define WG_Z 8
#define WG_SIZE        (WG_X * WG_Y * WG_Z)   // 512
#define TILE_SIZE      16
#define MAX_SCREENED   2048
#define MAX_TILE_PGTOS (TILE_SIZE * 32)         // 512 PGTOs per tile

// Global indices of CGTOs that pass the spatial screening test for this workgroup.
groupshared uint  screened_cgtos[MAX_SCREENED];
groupshared uint  num_screened_cgtos;

// PSI mode only: effective_coeffs[s] = sum_m C[m][ screened_cgtos[s] ]
// Computed cooperatively after screening, before the tile evaluation loop.
groupshared float effective_coeffs[MAX_SCREENED];

// Tile-level CGTO/PGTO data loaded cooperatively each iteration.
groupshared float3 tile_center[TILE_SIZE];
groupshared uint2  tile_pgto_range[TILE_SIZE];   // [begin, end) into pgtos_tile[]
groupshared PGTO   pgtos_tile[MAX_TILE_PGTOS];
groupshared uint   num_pgtos_tile;
groupshared uint2  gs_src_off_len[TILE_SIZE];    // [global_offset, length] per tile slot

// Cooperatively loads one tile's worth of CGTO and PGTO data into groupshared.
// Must be called with a preceding barrier that zeros num_pgtos_tile, and followed
// by a barrier before the results are read.
void fill_cgtos_tile(uint tile_number, uint tid) {
    // Phase 1: first TILE_SIZE threads fetch metadata and reserve pgtos_tile slots.
    if (tid < TILE_SIZE) {
        uint slot = tile_number * TILE_SIZE + tid;
        uint gi   = (slot < num_screened_cgtos) ? screened_cgtos[slot] : 0xFFFFFFFFu;

        float3 center  = float3(0.0f);
        uint2  off_len = uint2(0u, 0u);
        uint   dst     = 0u;

        if (gi != 0xFFFFFFFFu) {
            off_len = cgto_off_len[gi];
            center  = cgto_xyzr[gi].xyz;
            InterlockedAdd(num_pgtos_tile, off_len.y, dst);
        }

        gs_src_off_len[tid]  = off_len;
        tile_center[tid]     = center;
        tile_pgto_range[tid] = uint2(dst, dst + off_len.y);
    }

    // Ensure metadata is visible before all threads start copying PGTOs.
    GroupMemoryBarrierWithGroupSync();

    // Phase 2: all threads cooperatively copy PGTOs from global memory to pgtos_tile[].
    for (uint ci = 0u; ci < TILE_SIZE; ++ci) {
        uint2 src_ol = gs_src_off_len[ci];
        uint  src    = src_ol.x;
        uint  len    = src_ol.y;
        uint  dst    = tile_pgto_range[ci].x;
        for (uint k = tid; k < len; k += WG_SIZE) {
            pgtos_tile[dst + k] = pgto[src + k];
        }
    }
}

// Evaluate the radial+angular phi for every CGTO in the current tile at world position
// coord.  Reads tile_center[], tile_pgto_range[], and pgtos_tile[] from groupshared.
void eval_phi_tile(out float phi_tile[TILE_SIZE], float3 coord) {
    for (uint ci = 0u; ci < TILE_SIZE; ++ci) {
        float3 center = tile_center[ci];
        float3 d      = coord - center;
        float  r2     = dot(d, d);
        float  dx2    = d.x * d.x;
        float  dy2    = d.y * d.y;
        float  dz2    = d.z * d.z;

        float phi  = 0.0f;
        uint2 range = tile_pgto_range[ci];
        for (uint j = range.x; j < range.y; ++j) {
            PGTO  p    = pgtos_tile[j];
            uint4 ijkl = unpack_ijkl(p.ijkl);

            float fx = (ijkl.x == 0u) ? 1.0f :
                       (ijkl.x == 1u) ? d.x  :
                       (ijkl.x == 2u) ? dx2  :
                       (ijkl.x == 3u) ? dx2 * d.x : dx2 * dx2;

            float fy = (ijkl.y == 0u) ? 1.0f :
                       (ijkl.y == 1u) ? d.y  :
                       (ijkl.y == 2u) ? dy2  :
                       (ijkl.y == 3u) ? dy2 * d.y : dy2 * dy2;

            float fz = (ijkl.z == 0u) ? 1.0f :
                       (ijkl.z == 1u) ? d.z  :
                       (ijkl.z == 2u) ? dz2  :
                       (ijkl.z == 3u) ? dz2 * d.z : dz2 * dz2;

            phi += p.coeff * fx * fy * fz * exp(-p.alpha * r2);
        }
        phi_tile[ci] = phi;
    }
}

[numthreads(WG_X, WG_Y, WG_Z)]
void main(
    uint3 gid  : SV_GroupID,
    uint3 gtid : SV_GroupThreadID,
    uint  tid  : SV_GroupIndex,
    uint3 dtid : SV_DispatchThreadID)
{
    // -------------------------------------------------------------------------
    // Phase 0: initialise groupshared scalars
    // -------------------------------------------------------------------------
    if (tid == 0u) {
        num_screened_cgtos = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    // -------------------------------------------------------------------------
    // Phase 1: spatial screening
    // Compute the workgroup AABB in model space and collect CGTOs whose cutoff
    // spheres overlap it.  Identical logic to eval_gto_density.slang.
    // -------------------------------------------------------------------------
    float3 model_aabb_min = float3( gid                 * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1,1,1)) * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;

    for (uint i = tid; i < ubo.num_cgtos; i += WG_SIZE) {
        float4 cgto = cgto_xyzr[i];
        if (cgto.w == 0.0f) continue;

        float3 model_xyz = mul(ubo.world_to_model, float4(cgto.xyz, 1.0f)).xyz;
        float  r         = cgto.w;
        float3 d         = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;

        if (dot(d, d) < r * r) {
            uint idx;
            InterlockedAdd(num_screened_cgtos, 1u, idx);
            if (idx < (uint)MAX_SCREENED) {
                screened_cgtos[idx] = i;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();

    if (tid == 0u && num_screened_cgtos > (uint)MAX_SCREENED) {
        num_screened_cgtos = (uint)MAX_SCREENED;
    }
    GroupMemoryBarrierWithGroupSync();

    // -------------------------------------------------------------------------
    // Early exit when no CGTOs contribute to this tile
    // -------------------------------------------------------------------------
    uint3  vol_dim;
    out_vol.GetDimensions(vol_dim.x, vol_dim.y, vol_dim.z);

    float3 coord  = mul(ubo.index_to_world, float4(float3(dtid), 1.0f)).xyz;
    float  result = 0.0f;

    if (num_screened_cgtos == 0u) {
        if (all(dtid < vol_dim)) out_vol[dtid] = 0.0f;
        return;
    }

    uint num_tiles = (num_screened_cgtos + (uint)TILE_SIZE - 1u) / (uint)TILE_SIZE;

    if (ubo.eval_mode == 0u) {
        // =====================================================================
        // PSI MODE
        // Collapse the MO sum into a single effective coefficient per screened CGTO:
        //   effective_coeffs[s] = sum_m C[m][ screened_cgtos[s] ]
        // This allows a single pass over CGTO tiles independent of num_rows.
        // =====================================================================
        for (uint s = tid; s < num_screened_cgtos; s += WG_SIZE) {
            uint  c   = screened_cgtos[s];
            float sum = 0.0f;
            for (uint m = 0u; m < ubo.num_rows; ++m) {
                sum += coeffs[m * ubo.num_cgtos + c];
            }
            effective_coeffs[s] = sum;
        }
        GroupMemoryBarrierWithGroupSync();

        for (uint tile = 0u; tile < num_tiles; ++tile) {
            if (tid == 0u) num_pgtos_tile = 0u;
            GroupMemoryBarrierWithGroupSync();

            fill_cgtos_tile(tile, tid);
            GroupMemoryBarrierWithGroupSync();

            float phi_tile[TILE_SIZE];
            eval_phi_tile(phi_tile, coord);

            uint tile_base = tile * (uint)TILE_SIZE;
            for (uint ci = 0u; ci < (uint)TILE_SIZE; ++ci) {
                uint global_ci = tile_base + ci;
                if (global_ci >= num_screened_cgtos) break;
                result += phi_tile[ci] * effective_coeffs[global_ci];
            }
        }

    } else {
        // =====================================================================
        // PSI_SQUARED MODE
        // For each MO row m compute psi_m = sum_c C[m][c] * phi_c, then add
        // psi_m^2 to the running sum.  PGTO tile data is reloaded for each MO
        // row.  This is efficient for the typical 1–4 orbital case; for density
        // accumulation over many occupied orbitals, use the density shader.
        // =====================================================================
        for (uint mo = 0u; mo < ubo.num_rows; ++mo) {
            float psi    = 0.0f;
            uint  mo_off = mo * ubo.num_cgtos;

            for (uint tile = 0u; tile < num_tiles; ++tile) {
                if (tid == 0u) num_pgtos_tile = 0u;
                GroupMemoryBarrierWithGroupSync();

                fill_cgtos_tile(tile, tid);
                GroupMemoryBarrierWithGroupSync();

                float phi_tile[TILE_SIZE];
                eval_phi_tile(phi_tile, coord);

                uint tile_base = tile * (uint)TILE_SIZE;
                for (uint ci = 0u; ci < (uint)TILE_SIZE; ++ci) {
                    uint global_ci = tile_base + ci;
                    if (global_ci >= num_screened_cgtos) break;
                    uint cgto_idx = screened_cgtos[global_ci];
                    psi += coeffs[mo_off + cgto_idx] * phi_tile[ci];
                }
            }

            result += psi * psi;
        }
    }

    if (all(dtid < vol_dim)) {
        out_vol[dtid] = result;
    }
}

// eval_gto_mo.slang
// Evaluates molecular orbital electron density on a 3D grid.
//
//   result[voxel] = sum_m( (sum_c C[m][c] * phi_c(voxel))^2 )
//
// Design:
//   1) Screening pass does minimal work: append CGTO index if it overlaps tile AABB.
//   2) For each CGTO tile, all threads cooperatively build a flat PGTO list in
//      groupshared memory.  Each primitive entry is annotated with its parent CGTO
//      tile-slot index for center and coefficient lookup.
//   3) The evaluation loop is a uniform flat iteration over num_pgtos_tile — no CSR
//      indirection, no empty-slot guards, no serial CGTO nesting.
//   4) MO rows are processed in blocks; coefficients are staged in groupshared.

struct PGTO {
    float coeff;
    float alpha;
    uint  ijkl;
};

struct UBO {
    float4x4 world_to_model;  // world -> model space
    float4x4 index_to_world;  // grid index -> world
    float4   step;            // xyz = voxel spacing in model-space units
    uint     num_cgtos;
    uint     num_rows;
    uint     _pad0;
    uint     _pad1;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]] StructuredBuffer<uint>   cgto_atom_idx;
[[vk::binding(1)]] StructuredBuffer<float>  cgto_r;
[[vk::binding(2)]] StructuredBuffer<uint2>  cgto_off_len;
[[vk::binding(3)]] StructuredBuffer<PGTO>   pgto;
[[vk::binding(4)]] StructuredBuffer<float4> atom_xyz;
[[vk::binding(5)]] StructuredBuffer<float>  coeffs;   // float[num_rows][num_cgtos], row-major

[[vk::binding(16)]]  // image slot 0  (slot + MD_GPU_MAX_BIND_SLOTS)
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
#define MAX_TILE_PGTOS (TILE_SIZE * 32)        // 512 PGTOs per tile
#define MO_BLOCK       4
#define INVALID_IDX    0xFFFFFFFFu

// Screened CGTO indices for this workgroup.
groupshared uint screened_cgtos[MAX_SCREENED];
groupshared uint num_screened_cgtos;

// Tile-local CGTO headers.
groupshared uint   tile_cgto_idx[TILE_SIZE];
groupshared float4 tile_cgto_xyz[TILE_SIZE];    // xyz used, w ignored
groupshared uint   tile_src_pgto_beg[TILE_SIZE];
groupshared uint   tile_src_pgto_len[TILE_SIZE];
groupshared uint   tile_dst_pgto_beg[TILE_SIZE + 1];  // CSR prefix, [ci, ci+1)
groupshared uint   num_pgtos_tile;

// Flat PGTO tile pool in SoA layout.
groupshared float pgto_coeff_tile[MAX_TILE_PGTOS];  // 2 KB
groupshared float pgto_alpha_tile[MAX_TILE_PGTOS];  // 2 KB
groupshared uint  pgto_ijkl_tile [MAX_TILE_PGTOS];  // 2 KB
groupshared uint  pgto_slot_tile [MAX_TILE_PGTOS];  // 2 KB — parent CGTO slot (0..TILE_SIZE-1)

// Coefficient tile: coeff_tile[mb][ci] = C[mo_base+mb][tile_cgto_idx[ci]]
groupshared float coeff_tile[MO_BLOCK][TILE_SIZE];

void fill_pgto_tile(uint tile_number, uint tid) {
    // Phase 1: first TILE_SIZE threads gather CGTO headers for this tile.
    if (tid < TILE_SIZE) {
        uint slot = tile_number * TILE_SIZE + tid;

        uint   cgto_idx = INVALID_IDX;
        uint   src_beg  = 0u;
        uint   src_len  = 0u;
        float4 center   = float4(0.0f);

        if (slot < num_screened_cgtos) {
            cgto_idx = screened_cgtos[slot];
            uint2 off_len = cgto_off_len[cgto_idx];
            src_beg = off_len.x;
            src_len = off_len.y;
            center  = atom_xyz[cgto_atom_idx[cgto_idx]];
        }

        tile_cgto_idx[tid]     = cgto_idx;
        tile_src_pgto_beg[tid] = src_beg;
        tile_src_pgto_len[tid] = src_len;
        tile_cgto_xyz[tid]     = center;
    }

    GroupMemoryBarrierWithGroupSync();

    // Phase 2: build CSR destination ranges for the flat PGTO tile pool.
    if (tid == 0u) {
        uint running = 0u;
        for (uint ci = 0u; ci < TILE_SIZE; ++ci) {
            tile_dst_pgto_beg[ci] = running;

            uint len = tile_src_pgto_len[ci];
            uint avail = (running < (uint)MAX_TILE_PGTOS) ? ((uint)MAX_TILE_PGTOS - running) : 0u;
            uint clamped_len = min(len, avail);
            tile_src_pgto_len[ci] = clamped_len;
            running += clamped_len;
        }
        tile_dst_pgto_beg[TILE_SIZE] = running;
        num_pgtos_tile = running;
    }

    GroupMemoryBarrierWithGroupSync();

    // Phase 3: copy PGTO primitives into SoA tile pool.
    for (uint ci = 0u; ci < TILE_SIZE; ++ci) {
        uint src_beg = tile_src_pgto_beg[ci];
        uint len     = tile_src_pgto_len[ci];
        uint dst_beg = tile_dst_pgto_beg[ci];

        for (uint k = tid; k < len; k += WG_SIZE) {
            PGTO p = pgto[src_beg + k];
            uint dst = dst_beg + k;
            pgto_coeff_tile[dst] = p.coeff;
            pgto_alpha_tile[dst] = p.alpha;
            pgto_ijkl_tile[dst]  = p.ijkl;
            pgto_slot_tile[dst]  = ci;   // annotate with parent CGTO slot
        }
    }
}

[numthreads(WG_X, WG_Y, WG_Z)]
void main(
    uint3 gid  : SV_GroupID,
    uint3 gtid : SV_GroupThreadID,
    uint  tid  : SV_GroupIndex,
    uint3 dtid : SV_DispatchThreadID)
{
    if (tid == 0u) {
        num_screened_cgtos = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    // Minimal screening pass: test overlap and append index if passed.
    float3 model_aabb_min = float3( gid                 * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1,1,1)) * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;

    for (uint i = tid; i < ubo.num_cgtos; i += WG_SIZE) {
        float3 center = atom_xyz[cgto_atom_idx[i]].xyz;
        float  r      = cgto_r[i];
        if (r == 0.0f) continue;

        float3 model_xyz = mul(ubo.world_to_model, float4(center, 1.0f)).xyz;
        float3 d = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;

        if (dot(d, d) < r * r) {
            uint idx;
            InterlockedAdd(num_screened_cgtos, 1u, idx);
            if (idx < (uint)MAX_SCREENED) {
                screened_cgtos[idx] = i;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();

    if (tid == 0u) {
        num_screened_cgtos = min(num_screened_cgtos, (uint)MAX_SCREENED);
    }
    GroupMemoryBarrierWithGroupSync();

    uint3 vol_dim;
    out_vol.GetDimensions(vol_dim.x, vol_dim.y, vol_dim.z);

    if (num_screened_cgtos == 0u) {
        if (all(dtid < vol_dim)) {
            out_vol[dtid] = 0.0f;
        }
        return;
    }

    float3 coord = mul(ubo.index_to_world, float4(float3(dtid), 1.0f)).xyz;
    float  result = 0.0f;
    uint   num_tiles = (num_screened_cgtos + (uint)TILE_SIZE - 1u) / (uint)TILE_SIZE;

    // MO-blocked contraction: reuse phi over a block of MO rows.
    for (uint mo_base = 0u; mo_base < ubo.num_rows; mo_base += (uint)MO_BLOCK) {
        uint mo_count = min((uint)MO_BLOCK, ubo.num_rows - mo_base);
        float psi_block[MO_BLOCK] = {0.0f, 0.0f, 0.0f, 0.0f};

        for (uint tile = 0u; tile < num_tiles; ++tile) {
            fill_pgto_tile(tile, tid);
            GroupMemoryBarrierWithGroupSync();

            // Load MO coefficients for this tile's CGTOs into groupshared.
            uint work_count = mo_count * (uint)TILE_SIZE;
            for (uint work = tid; work < work_count; work += WG_SIZE) {
                uint mb       = work / (uint)TILE_SIZE;
                uint ci       = work % (uint)TILE_SIZE;
                uint cgto_idx = tile_cgto_idx[ci];
                float c = (cgto_idx != INVALID_IDX)
                        ? coeffs[(mo_base + mb) * ubo.num_cgtos + cgto_idx]
                        : 0.0f;
                coeff_tile[mb][ci] = c;
            }
            GroupMemoryBarrierWithGroupSync();

            // Flat primitive evaluation loop — uniform iteration, no CSR indirection,
            // no empty-slot guards.  Each primitive carries its parent slot index for
            // center and coefficient lookup.
            uint num_prims = num_pgtos_tile;
            for (uint j = 0u; j < num_prims; ++j) {
                uint   slot = pgto_slot_tile[j];
                float3 d    = coord - tile_cgto_xyz[slot].xyz;
                float  r2   = dot(d, d);
                float  dx2  = d.x * d.x;
                float  dy2  = d.y * d.y;
                float  dz2  = d.z * d.z;

                uint4 ijkl = unpack_ijkl(pgto_ijkl_tile[j]);

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

                float prim_val = pgto_coeff_tile[j] * fx * fy * fz * exp(-pgto_alpha_tile[j] * r2);

                for (uint mb = 0u; mb < mo_count; ++mb) {
                    psi_block[mb] += coeff_tile[mb][slot] * prim_val;
                }
            }
            GroupMemoryBarrierWithGroupSync();
        }

        for (uint mb = 0u; mb < mo_count; ++mb) {
            result += psi_block[mb] * psi_block[mb];
        }
    }

    if (all(dtid < vol_dim)) {
        out_vol[dtid] = result;
    }
}

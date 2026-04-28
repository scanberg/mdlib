// eval_gto_density.slang
// Evaluates electron density from density matrix and contracted GTO basis on a 3D grid.
// Ported from eval_gto_density.comp to Slang for md_gpu backend.

struct PGTO {
    float coeff;
    float alpha;
    uint ijkl;
};

struct CGTO {
    float3 center;
    uint2  pgto_range; // beg + end
};

struct UBO {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint     D_matrix_dim;
    uint     _pad0;
    uint     _pad1;
    uint     _pad2;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
StructuredBuffer<float4> cgto_xyzr;

[[vk::binding(1)]]
StructuredBuffer<uint2> cgto_off_len;

[[vk::binding(2)]]
StructuredBuffer<PGTO> pgto;

[[vk::binding(3)]]
StructuredBuffer<float> D_matrix;

[[vk::binding(16)]] // image slot 0 (slot + MD_GPU_MAX_BIND_SLOTS)
RWTexture3D<float> out_vol;

uint4 unpack_ijkl(uint packed) {
    return uint4(
        (packed >>  0) & 0xFFu,
        (packed >>  8) & 0xFFu,
        (packed >> 16) & 0xFFu,
        (packed >> 24) & 0xFFu
    );
}

#define WG_X 8
#define WG_Y 8
#define WG_Z 8

#define WG_SIZE (WG_X * WG_Y * WG_Z)
#define TILE_SIZE 16
#define TILE_AREA (TILE_SIZE*TILE_SIZE)
#define MAX_SCREENED_CGTOS 2048
#define MAX_TILE_PGTOS (TILE_SIZE * 32)
#define INVALID_CGTO_IDX 0xFFFFFFFFu

groupshared uint  screened_cgtos[MAX_SCREENED_CGTOS];
groupshared uint  num_screened_cgtos;

groupshared float D_tile[TILE_SIZE][TILE_SIZE];
groupshared CGTO  cgtos_tile[TILE_SIZE];
groupshared PGTO  pgtos_tile[MAX_TILE_PGTOS];
groupshared uint  num_pgtos;

// Per-CGTO source ranges written by phase 1 of fill_cgtos_tile,
// read by all threads in phase 2.
groupshared uint2  gs_src_off_len[TILE_SIZE];

// Upper triangular symmetric matrix index
uint get_index(uint i, uint j, uint N) {
    uint row = (i < j) ? i : j;
    uint col = (i < j) ? j : i;
    uint row_offset = row * (2 * N - row + 1) / 2;
    return row_offset + (col - row);
}

void fill_D_tile(uint tile_i, uint tile_j, uint tid) {
    uint baseI = tile_i * TILE_SIZE;
    uint baseJ = tile_j * TILE_SIZE;

    for (uint index = tid; index < TILE_AREA; index += WG_SIZE) {
        uint row = index / TILE_SIZE;
        uint col = index % TILE_SIZE;

        uint idxI = baseI + row;
        uint gi = (idxI < num_screened_cgtos) ? screened_cgtos[idxI] : INVALID_CGTO_IDX;

        uint idxJ = baseJ + col;
        uint gj = (idxJ < num_screened_cgtos) ? screened_cgtos[idxJ] : INVALID_CGTO_IDX;

        float value = 0.0;
        if (gi != INVALID_CGTO_IDX && gj != INVALID_CGTO_IDX) {
            value = D_matrix[get_index(gi, gj, ubo.D_matrix_dim)];
        }
        D_tile[row][col] = value;
    }
}

// Two-phase CGTO/PGTO loader that utilises all WG_SIZE threads.
//
// Phase 1 (threads 0 .. TILE_SIZE-1):
//   Read CGTO metadata, atomically reserve pgtos_tile slots, write
//   cgtos_tile[] headers and gs_src_off_len[].
//
// [barrier]
//
// Phase 2 (all WG_SIZE threads):
//   Cooperatively copy every PGTO belonging to this tile from the
//   global buffers into pgtos_tile using a strided loop.
//   Per-PGTO radius filtering is intentionally omitted here: the
//   CGTO-level radius (cgto_xyzr.w = max PGTO radius) already screens
//   out non-contributing CGTOs during Stage 1, so any PGTO whose
//   parent CGTO survived screening is likely relevant. Tight PGTOs
//   that still evaluate to near-zero contribute negligible error via
//   exp(-alpha * r2) and avoiding the filtering removes a bottleneck.
void fill_cgtos_tile(uint tile_number, uint tid) {
    // --- Phase 1: metadata gather (TILE_SIZE threads) ---
    if (tid < TILE_SIZE) {
        uint slot = tile_number * TILE_SIZE + tid;
        uint global_cgto_idx = (slot < num_screened_cgtos) ? screened_cgtos[slot] : INVALID_CGTO_IDX;

        float3 cgto_center = float3(0.0);
        uint2 off_len = uint2(0, 0);
        uint dst     = 0u;

        if (global_cgto_idx != INVALID_CGTO_IDX) {
            off_len = cgto_off_len[global_cgto_idx];
            cgto_center = cgto_xyzr[global_cgto_idx].xyz;
            InterlockedAdd(num_pgtos, off_len.y, dst);
        }

        gs_src_off_len[tid] = off_len;
        cgtos_tile[tid].center = cgto_center;
        cgtos_tile[tid].pgto_range = uint2(dst, dst + off_len.y);
    }

    // Ensure gs_src_beg/gs_src_len and cgtos_tile headers are visible
    // to all threads before phase 2.
    GroupMemoryBarrierWithGroupSync();

    // --- Phase 2: cooperative PGTO load (all WG_SIZE threads) ---
    for (uint ci = 0; ci < TILE_SIZE; ci++) {
        uint2 off_len = gs_src_off_len[ci];
        uint src = off_len.x;
        uint len = off_len.y;
        uint dst = cgtos_tile[ci].pgto_range.x;
        for (uint k = tid; k < len; k += WG_SIZE) {
            PGTO g = pgto[src + k];
            pgtos_tile[dst + k] = g;
        }
    }
}

void eval_phis(out float out_phi[TILE_SIZE], float3 coord) {
    for (uint i = 0; i < TILE_SIZE; ++i) {
        uint2 pgto_range = cgtos_tile[i].pgto_range;

        float3 center = cgtos_tile[i].center;
        float3 d   = coord - center;
        float r2 = dot(d, d);

        float dx2 = d.x * d.x;
        float dy2 = d.y * d.y;
        float dz2 = d.z * d.z;

        float phi = 0.0;
        for (uint j = pgto_range.x; j < pgto_range.y; ++j) {
            PGTO pgto = pgtos_tile[j];
            uint4 ijkl = unpack_ijkl(pgto.ijkl);

            float fx = (ijkl.x == 0u) ? 1.0 :
                       (ijkl.x == 1u) ? d.x :
                       (ijkl.x == 2u) ? dx2 :
                       (ijkl.x == 3u) ? dx2 * d.x : dx2 * dx2;

            float fy = (ijkl.y == 0u) ? 1.0 :
                       (ijkl.y == 1u) ? d.y :
                       (ijkl.y == 2u) ? dy2 :
                       (ijkl.y == 3u) ? dy2 * d.y : dy2 * dy2;

            float fz = (ijkl.z == 0u) ? 1.0 :
                       (ijkl.z == 1u) ? d.z :
                       (ijkl.z == 2u) ? dz2 :
                       (ijkl.z == 3u) ? dz2 * d.z : dz2 * dz2;

            phi += pgto.coeff * fx * fy * fz * exp(-pgto.alpha * r2);
        }
        out_phi[i] = phi;
    }
}

float symmetric_contraction(float phi[TILE_SIZE]) {
    float result = 0.0;
    for (uint i = 0; i < TILE_SIZE; ++i) {
        float ai = phi[i];
        result += D_tile[i][i] * ai * ai;
        float ai2 = 2.0 * ai;
        for (uint j = i + 1; j < TILE_SIZE; ++j) {
            result += D_tile[i][j] * ai2 * phi[j];
        }
    }
    return result;
}

float contraction(float phi_mu[TILE_SIZE], float phi_nu[TILE_SIZE]) {
    float result = 0.0;
    for (uint i = 0; i < TILE_SIZE; ++i) {
        float ai = 2.0 * phi_mu[i];
        for (uint j = 0; j < TILE_SIZE; ++j) {
            result += D_tile[i][j] * ai * phi_nu[j];
        }
    }
    return result;
}

[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex, uint3 dtid : SV_DispatchThreadID) {
    if (tid == 0) {
        num_screened_cgtos = 0;
    }
    GroupMemoryBarrierWithGroupSync();

    float3 model_aabb_min = float3( gid * uint3(WG_X, WG_Y, WG_Z))              * ubo.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1,1,1)) * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;

    // Stage 1: Screening - Prune CGTOs relevant for this workgroup region
    {
        for (uint i = tid; i < ubo.D_matrix_dim; i += WG_SIZE) {
            float4 cgto = cgto_xyzr[i];
            if (cgto.w == 0.0) continue;

            float3 model_xyz = mul(ubo.world_to_model, float4(cgto.xyz, 1.0)).xyz;
            float r2 = cgto.w * cgto.w;

            float3 d = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;
            float d2 = dot(d, d);

            if (d2 < r2) {
                uint idx;
                InterlockedAdd(num_screened_cgtos, 1u, idx);
                if (idx < MAX_SCREENED_CGTOS) {
                    screened_cgtos[idx] = i;
                }
            }
        }

        GroupMemoryBarrierWithGroupSync();

        if (tid == 0) {
            if (num_screened_cgtos > MAX_SCREENED_CGTOS) {
                num_screened_cgtos = MAX_SCREENED_CGTOS;
            }
        }
        GroupMemoryBarrierWithGroupSync();
    }

    uint3 vol_dim;
    out_vol.GetDimensions(vol_dim.x, vol_dim.y, vol_dim.z);

    float phi_tile_mu[TILE_SIZE];
    float phi_tile_nu[TILE_SIZE];
    float3 coord = mul(ubo.index_to_world, float4(float3(dtid), 1.0)).xyz;
    float rho = 0.0;

    if (num_screened_cgtos == 0u) {
        if (all(dtid < vol_dim)) {
            out_vol[dtid] = 0.0;
        }
        return;
    }

    uint num_tiles = (num_screened_cgtos + TILE_SIZE - 1) / TILE_SIZE;

    for (uint tile_i = 0; tile_i < num_tiles; ++tile_i) {
        // DIAGONAL TILE
        if (tid == 0) {
            num_pgtos = 0;
        }
        GroupMemoryBarrierWithGroupSync();

        fill_D_tile(tile_i, tile_i, tid);
        fill_cgtos_tile(tile_i, tid);

        GroupMemoryBarrierWithGroupSync();

        eval_phis(phi_tile_mu, coord);

        rho += symmetric_contraction(phi_tile_mu);

        for (uint tile_j = tile_i + 1; tile_j < num_tiles; ++tile_j) {
            // OFF DIAGONAL TILE
            if (tid == 0) {
                num_pgtos = 0;
            }
            GroupMemoryBarrierWithGroupSync();

            fill_D_tile(tile_i, tile_j, tid);
            fill_cgtos_tile(tile_j, tid);

            GroupMemoryBarrierWithGroupSync();

            eval_phis(phi_tile_nu, coord);

            rho += contraction(phi_tile_mu, phi_tile_nu);
        }
    }

    if (all(dtid < vol_dim)) {
        out_vol[dtid] = rho;
    }
}

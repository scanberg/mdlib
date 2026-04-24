// eval_gto_density.slang
// Evaluates electron density from density matrix and contracted GTO basis on a 3D grid.
// Ported from eval_gto_density.comp to Slang for md_gpu backend.

struct PGTO {
    float coeff;
    float alpha;
    uint ijkl;
};

struct CGTO {
    float3 coord;
    uint pgto_off;
    uint pgto_len;
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
StructuredBuffer<uint> cgto_offset;

[[vk::binding(2)]]
StructuredBuffer<float> pgto_coeff;

[[vk::binding(3)]]
StructuredBuffer<float> pgto_alpha;

[[vk::binding(4)]]
StructuredBuffer<float> pgto_radius;

[[vk::binding(5)]]
StructuredBuffer<uint> pgto_ijkl;

[[vk::binding(6)]]
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
#define MAX_TILE_PGTOS (TILE_SIZE * 16)
#define INVALID_CGTO_IDX 0xFFFFFFFFu

groupshared uint  screened_cgtos[MAX_SCREENED_CGTOS];
groupshared uint  num_screened_cgtos;

groupshared float D_tile[TILE_SIZE][TILE_SIZE];
groupshared CGTO  cgtos_tile[TILE_SIZE];
groupshared PGTO  pgtos_tile[MAX_TILE_PGTOS];
groupshared uint  num_pgtos;

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
        uint gi = (idxI < MAX_SCREENED_CGTOS) ? screened_cgtos[idxI] : INVALID_CGTO_IDX;

        uint idxJ = baseJ + col;
        uint gj = (idxJ < MAX_SCREENED_CGTOS) ? screened_cgtos[idxJ] : INVALID_CGTO_IDX;

        float value = 0.0;
        if (gi != INVALID_CGTO_IDX && gj != INVALID_CGTO_IDX) {
            value = D_matrix[get_index(gi, gj, ubo.D_matrix_dim)];
        }
        D_tile[row][col] = value;
    }
}

void fill_cgtos_tile(uint tile_number, float3 model_aabb_min, float3 model_aabb_max, uint tid) {
    if (tid >= TILE_SIZE) return;

    uint global_cgto_idx = screened_cgtos[tile_number * TILE_SIZE + tid];

    float3 cgto_center = float3(0.0);
    uint cgto_pgto_off = 0u;
    uint cgto_pgto_len = 0u;

    if (global_cgto_idx != INVALID_CGTO_IDX) {
        uint src_beg = cgto_offset[global_cgto_idx];
        uint src_end = cgto_offset[global_cgto_idx + 1u];
        uint src_len = src_end - src_beg;

        uint dst = 0;
        InterlockedAdd(num_pgtos, src_len, dst);

        cgto_center = cgto_xyzr[global_cgto_idx].xyz;
        cgto_pgto_off = dst;

        float3 model_xyz = mul(ubo.world_to_model, float4(cgto_center, 1.0)).xyz;
        float3 d = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;
        float d2 = dot(d, d);

        for (uint k = 0; k < src_len; ++k) {
            float r = pgto_radius[src_beg + k];
            if (d2 < r * r) {
                PGTO g;
                g.coeff = pgto_coeff[src_beg + k];
                g.alpha = pgto_alpha[src_beg + k];
                g.ijkl  = pgto_ijkl[src_beg + k];
                pgtos_tile[cgto_pgto_off + cgto_pgto_len] = g;
                cgto_pgto_len++;
            }
        }
    }

    cgtos_tile[tid].coord    = cgto_center;
    cgtos_tile[tid].pgto_off = cgto_pgto_off;
    cgtos_tile[tid].pgto_len = cgto_pgto_len;
}

void eval_phis(out float out_phi[TILE_SIZE], float3 coord) {
    [unroll]
    for (uint i = 0; i < TILE_SIZE; ++i) out_phi[i] = 0.0;

    for (uint i = 0; i < TILE_SIZE; ++i) {
        uint pgto_off = cgtos_tile[i].pgto_off;
        uint pgto_len = cgtos_tile[i].pgto_len;

        float3 center = cgtos_tile[i].coord;
        float3 d   = coord - center;
        float r2 = dot(d, d);

        float dx2 = d.x * d.x;
        float dy2 = d.y * d.y;
        float dz2 = d.z * d.z;

        float phi = 0.0;
        for (uint j = pgto_off; j < pgto_off + pgto_len; ++j) {
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

        for (uint i = tid; i < MAX_SCREENED_CGTOS; i += WG_SIZE) {
            if (i >= num_screened_cgtos) {
                screened_cgtos[i] = INVALID_CGTO_IDX;
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
        fill_cgtos_tile(tile_i, model_aabb_min, model_aabb_max, tid);

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
            fill_cgtos_tile(tile_j, model_aabb_min, model_aabb_max, tid);

            GroupMemoryBarrierWithGroupSync();

            eval_phis(phi_tile_nu, coord);

            rho += contraction(phi_tile_mu, phi_tile_nu);
        }
    }

    if (all(dtid < vol_dim)) {
        out_vol[dtid] = rho;
    }
}

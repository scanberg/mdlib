// eval_gto_density.slang
// Evaluates electron density from density matrix and contracted GTO basis on a 3D grid.
// Ported from eval_gto_density.comp to Slang for md_gpu backend.

struct PGTO {
    float coeff;
    float alpha;
    uint ijkl;
};

struct UBO {
    float4x4 world_to_model;
    float4x4 index_to_world;
    float4   step;
    uint     num_cgtos;
    uint     _pad0;
    uint     _pad1;
    uint     _pad2;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
StructuredBuffer<uint> cgto_atom_idx;

[[vk::binding(1)]]
StructuredBuffer<float> cgto_r;

[[vk::binding(2)]]
StructuredBuffer<uint2> cgto_off_len;

[[vk::binding(3)]]
StructuredBuffer<PGTO> pgto;

[[vk::binding(4)]]
StructuredBuffer<float4> atom_xyz;

[[vk::binding(5)]]
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
#define MAX_PGTOS_PER_CGTO 16
#define MAX_TILE_PGTOS (TILE_SIZE * MAX_PGTOS_PER_CGTO)

groupshared uint  screened_cgtos[MAX_SCREENED_CGTOS];
groupshared uint  num_screened_cgtos;

groupshared float  D_tile[TILE_SIZE][TILE_SIZE];
groupshared float3 tile_cgto_xyz[TILE_SIZE];
groupshared uint2  tile_pgto_range[TILE_SIZE];

// PGTO data stored as struct-of-arrays in groupshared:
//  - avoids the AoS struct padding (PGTO is 12 bytes, often padded to 16 in LDS)
//  - lets the compiler issue independent LDS loads per channel and overlap
//    them with FMA in the eval_phis hot loop.
groupshared float pgtos_tile_coeff[MAX_TILE_PGTOS];
groupshared float pgtos_tile_alpha[MAX_TILE_PGTOS];
groupshared uint  pgtos_tile_ijkl [MAX_TILE_PGTOS];
groupshared uint  num_pgtos;

groupshared uint  pgtos_tile_offset[TILE_SIZE];

// Upper triangular symmetric matrix index
uint get_index(uint i, uint j, uint N) {
    uint row = (i < j) ? i : j;
    uint col = (i < j) ? j : i;
    uint row_offset = row * (2 * N - row + 1) / 2;
    return row_offset + (col - row);
}

void fill_D_tile(uint tile_i, uint tile_j, uint num_cgtos, uint tid) {
    uint baseI = tile_i * TILE_SIZE;
    uint baseJ = tile_j * TILE_SIZE;

    for (uint index = tid; index < TILE_AREA; index += WG_SIZE) {
        uint row = index / TILE_SIZE;
        uint col = index % TILE_SIZE;
        uint slot_i = baseI + row;
        uint slot_j = baseJ + col;

        float value = 0.0;
        if (slot_i < num_cgtos && slot_j < num_cgtos) {
            uint gi = screened_cgtos[slot_i];
            uint gj = screened_cgtos[slot_j];
            value = D_matrix[get_index(gi, gj, ubo.num_cgtos)];
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
// Two-phase CGTO/PGTO tile loader.
// Phase 1 (first TILE_SIZE threads): fetch CGTO metadata, reserve pgtos_tile slots.
// Phase 2 (all WG_SIZE threads): cooperatively copy PGTOs from global memory.
void fill_cgtos_tile(uint tile_offset, uint num_cgtos, uint tid) {
    // ---- Phase 1: per-CGTO metadata + slot reservation ------------------
    if (tid < TILE_SIZE) {
        uint slot = tile_offset * TILE_SIZE + tid;

        uint src = 0u;
        uint len = 0u;
        uint dst = 0u;
        float3 center = float3(0.0f);

        if (slot < num_cgtos) {
            uint gi = screened_cgtos[slot];
            uint2 off_len = cgto_off_len[gi];
            src = off_len.x;
            len = off_len.y;
            center = atom_xyz[cgto_atom_idx[gi]].xyz;
            InterlockedAdd(num_pgtos, len, dst);
        }

        pgtos_tile_offset[tid] = src;
        tile_pgto_range[tid] = uint2(dst, dst + len);
        tile_cgto_xyz[tid] = center;
    }

    // Make Phase 1 results visible to all threads before Phase 2 reads them.
    GroupMemoryBarrierWithGroupSync();

    // ---- Phase 2: cooperative PGTO copy using all WG_SIZE threads -------
    // Each CGTO has up to MAX_PGTOS_PER_CGTO PGTOs.
    // i becomes the CGTO index, j the intra-CGTO PGTO index.
    uint i = tid / MAX_PGTOS_PER_CGTO;
    if (i < TILE_SIZE) {
        uint j = tid % MAX_PGTOS_PER_CGTO;
        uint2 r = tile_pgto_range[i];
        if (r.x + j < r.y) {
            uint dst = r.x + j;
            uint src = pgtos_tile_offset[i] + j;
            PGTO p = pgto[src];
            pgtos_tile_coeff[dst] = p.coeff;
            pgtos_tile_alpha[dst] = p.alpha;
            pgtos_tile_ijkl [dst] = p.ijkl;
        }
    }
}
// Pow up to 4 via a small unrolled switch. Compiles to the same code as the
// original ?: ladder but is easier on the optimizer.
float pow_le4(float x, float x2, uint n) {
    switch (n) {
        case 0u: return 1.0f;
        case 1u: return x;
        case 2u: return x2;
        case 3u: return x2 * x;
        default: return x2 * x2;
    }
}

// Evaluate phi for every CGTO slot in the current tile.
// Slots beyond num_screened_cgtos have an empty pgto_range so their phi is 0.
void eval_phis(out float out_phi[TILE_SIZE], float3 coord) {
    //[unroll]
    for (uint i = 0u; i < TILE_SIZE; ++i) {
        float3 d   = coord - tile_cgto_xyz[i];
        float  r2  = dot(d, d);
        float  dx2 = d.x * d.x;
        float  dy2 = d.y * d.y;
        float  dz2 = d.z * d.z;

        uint2 range = tile_pgto_range[i];
        uint  base  = range.x;
        uint  cnt   = range.y - range.x;

        float phi = 0.0f;
        // Bounded, fully-unrolled inner loop over PGTOs (<= MAX_PGTOS_PER_CGTO).
        // The early break lets the compiler keep `cnt` in a scalar register and
        // hoist the per-PGTO loads from groupshared into registers.
        [unroll]
        for (uint j = 0u; j < MAX_PGTOS_PER_CGTO; ++j) {
            if (j >= cnt) break;

            uint  slot  = base + j;
            float coeff = pgtos_tile_coeff[slot];
            float alpha = pgtos_tile_alpha[slot];
            uint4 ijkl  = unpack_ijkl(pgtos_tile_ijkl[slot]);

            float fx = pow_le4(d.x, dx2, ijkl.x);
            float fy = pow_le4(d.y, dy2, ijkl.y);
            float fz = pow_le4(d.z, dz2, ijkl.z);

            phi += coeff * fx * fy * fz * exp(-alpha * r2);
        }
        out_phi[i] = phi;
    }
}

float symmetric_contraction(float phi[TILE_SIZE]) {
    float result = 0.0;
    //[unroll]
    for (uint i = 0; i < TILE_SIZE; ++i) {
        float ai = phi[i];
        result += D_tile[i][i] * ai * ai;
        float ai2 = 2.0 * ai;
        [unroll]
        for (uint j = i + 1; j < TILE_SIZE; ++j) {
            result += D_tile[i][j] * ai2 * phi[j];
        }
    }
    return result;
}

float contraction(float phi_mu[TILE_SIZE], float phi_nu[TILE_SIZE]) {
    float result = 0.0;
    //[unroll]
    for (uint i = 0; i < TILE_SIZE; ++i) {
        float ai = 2.0 * phi_mu[i];
        [unroll]
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

    float3 model_aabb_min = float3( gid * uint3(WG_X, WG_Y, WG_Z))                 * ubo.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1,1,1)) * uint3(WG_X, WG_Y, WG_Z)) * ubo.step.xyz;

    // Stage 1: Screening - Prune CGTOs relevant for this workgroup region.
    // Wave-coalesced allocation: one LDS atomic per wave instead of per thread.
    for (uint i = tid; i < ubo.num_cgtos; i += WG_SIZE) {
        bool passed = false;
        if (i < ubo.num_cgtos) {
            float r = cgto_r[i];
            if (r != 0.0) {
                uint   atom_idx  = cgto_atom_idx[i];
                float3 center    = atom_xyz[atom_idx].xyz;
                float3 model_xyz = mul(ubo.world_to_model, float4(center, 1.0)).xyz;
                float  r2        = r * r;
                float3 d         = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;
                passed = (dot(d, d) < r2);
            }
        }

        uint wave_count  = WaveActiveCountBits(passed);
        uint wave_offset = WavePrefixSum(passed ? 1u : 0u);
        uint wave_base   = 0u;
        if (WaveIsFirstLane() && wave_count > 0u) {
            InterlockedAdd(num_screened_cgtos, wave_count, wave_base);
        }
        wave_base = WaveReadLaneFirst(wave_base);

        if (passed) {
            uint idx = wave_base + wave_offset;
            if (idx < MAX_SCREENED_CGTOS) {
                screened_cgtos[idx] = i;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();

    uint num_cgtos = min(num_screened_cgtos, MAX_SCREENED_CGTOS);

    uint3 vol_dim;
    out_vol.GetDimensions(vol_dim.x, vol_dim.y, vol_dim.z);
    bool active = all(dtid < vol_dim);

    if (num_cgtos == 0u) {
        if (active) {
            out_vol[dtid] = 0.0;
        }
        return;
    }

    float phi_tile_mu[TILE_SIZE];
    float phi_tile_nu[TILE_SIZE];
    float3 coord = mul(ubo.index_to_world, float4(float3(dtid), 1.0)).xyz;
    float rho = 0.0;

    uint num_tiles = (num_cgtos + TILE_SIZE - 1) / TILE_SIZE;

    for (uint tile_i = 0; tile_i < num_tiles; ++tile_i) {
        // DIAGONAL TILE
        if (tid == 0) {
            num_pgtos = 0;
        }
        GroupMemoryBarrierWithGroupSync();

        fill_D_tile(tile_i, tile_i, num_cgtos, tid);
        fill_cgtos_tile(tile_i, num_cgtos, tid);

        GroupMemoryBarrierWithGroupSync();

        // Inactive (out-of-bounds) voxels must still participate in the
        // collaborative fills + barriers above, but skip the per-voxel work.
        if (active) {
            eval_phis(phi_tile_mu, coord);
            rho += symmetric_contraction(phi_tile_mu);
        }

        for (uint tile_j = tile_i + 1; tile_j < num_tiles; ++tile_j) {
            // OFF DIAGONAL TILE
            if (tid == 0) {
                num_pgtos = 0;
            }
            GroupMemoryBarrierWithGroupSync();

            fill_D_tile(tile_i, tile_j, num_cgtos, tid);
            fill_cgtos_tile(tile_j, num_cgtos, tid);

            GroupMemoryBarrierWithGroupSync();

            if (active) {
                eval_phis(phi_tile_nu, coord);
                rho += contraction(phi_tile_mu, phi_tile_nu);
            }
        }
    }

    if (active) {
        out_vol[dtid] = rho;
    }
}

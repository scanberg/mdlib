// eval_gto_mo.slang
// Evaluates molecular orbital wavefunction or electron density on a 3D grid.
//
//   mode == 0 (PSI):         result[voxel] = sum_m( sum_c C[m][c] * phi_c(voxel) )
//   mode == 1 (PSI_SQUARED): result[voxel] = sum_m( (sum_c C[m][c] * phi_c(voxel))^2 )
//
// Design:
//   1) Screening pass does minimal work: append CGTO index if it overlaps tile AABB.
//   2) For each CGTO tile, all threads cooperatively build a flat PGTO list in
//      groupshared memory.  Each primitive entry is annotated with its parent CGTO
//      tile-slot index for center and coefficient lookup.
//   3) Phase 3 of fill_pgto_tile is a flat scatter: all WG_SIZE threads iterate
//      output indices in stride WG_SIZE and look up the owning CGTO slot from
//      the prefix array — 100% utilisation regardless of per-CGTO primitive counts.
//   4) MO rows are processed in blocks; coefficients are staged in groupshared.

#define OP_SET 0
#define OP_ADD 1
#define OP_SUB 2
#define OP_MAX 3
#define OP_MIN 4
#define OP_ACCUM_MASK 0x0000000Fu
#define OP_ABS 0x00000100u

struct PGTO {
    float coeff;
    float alpha;
    uint  ijkl;
};

struct RootArgs {
    float4x4 world_to_model;  // world -> model space
    float4x4 index_to_world;  // grid index -> world
    float4   step;            // xyz = voxel spacing in model-space units
    uint4    grid_dim;
    uint     num_cgtos;
    uint     num_rows;
    uint     mode;            // 0 = PSI (raw wavefunction), 1 = PSI_SQUARED (density)
    uint     operation;       // low bits = set/add/sub/max/min, high bits = transforms

    uint*   cgto_atom_idx;
    float*  cgto_r;
    uint2*  cgto_off_len;
    PGTO*   pgto;
    float4* atom_xyz;
    float*  coeffs;
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

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
groupshared uint   tile_dst_pgto_beg[TILE_SIZE + 1];  // CSR prefix, [ci, ci+1)
groupshared uint   num_pgtos_tile;

// Flat PGTO tile pool in SoA layout.
groupshared float pgto_coeff_tile[MAX_TILE_PGTOS];  // 2 KB
groupshared float pgto_alpha_tile[MAX_TILE_PGTOS];  // 2 KB
groupshared uint  pgto_ijkl_tile [MAX_TILE_PGTOS];  // 2 KB
groupshared uint  pgto_slot_tile [MAX_TILE_PGTOS];  // 2 KB — parent CGTO slot (0..TILE_SIZE-1)

// Coefficient tile: coeff_tile[mb][ci] = C[mo_base+mb][tile_cgto_idx[ci]]
groupshared float coeff_tile[MO_BLOCK][TILE_SIZE];

void fill_pgto_tile(uint tile_number, uint num_cgtos, uint tid) {
    // Phase 1: first TILE_SIZE threads gather CGTO headers for this tile.
    if (tid < TILE_SIZE) {
        uint slot = tile_number * TILE_SIZE + tid;

        uint   cgto_idx = INVALID_IDX;
        uint   src_beg  = 0u;
        uint   src_len  = 0u;
        float4 center   = float4(0.0f);

        if (slot < num_cgtos) {
            cgto_idx = screened_cgtos[slot];
            uint2 off_len = args.cgto_off_len[cgto_idx];
            src_beg = off_len.x;
            src_len = off_len.y;
            center  = args.atom_xyz[args.cgto_atom_idx[cgto_idx]];
        }

        tile_cgto_idx[tid]     = cgto_idx;
        tile_src_pgto_beg[tid] = src_beg;
        tile_dst_pgto_beg[tid] = src_len;   // temporarily store len; Phase 2 will prefix-sum
        tile_cgto_xyz[tid]     = center;
    }

    GroupMemoryBarrierWithGroupSync();

    // Phase 2: build CSR destination ranges from lengths stored in tile_dst_pgto_beg[0..TILE_SIZE-1].
    if (tid == 0u) {
        uint running = 0u;
        for (uint ci = 0u; ci < TILE_SIZE; ++ci) {
            uint len = tile_dst_pgto_beg[ci];   // read the temporarily stored len
            uint avail = (running < (uint)MAX_TILE_PGTOS) ? ((uint)MAX_TILE_PGTOS - running) : 0u;
            uint clamped_len = min(len, avail);
            tile_dst_pgto_beg[ci] = running;    // overwrite with the prefix offset
            running += clamped_len;
        }
        tile_dst_pgto_beg[TILE_SIZE] = running;
        num_pgtos_tile = running;
    }

    GroupMemoryBarrierWithGroupSync();

    // Phase 3: flat scatter — all threads cooperate over the output range.
    // A branchless unrolled scan over the TILE_SIZE+1 prefix entries finds the
    // owning CGTO slot for each output index without any serial outer loop.
    for (uint k = tid; k < num_pgtos_tile; k += WG_SIZE) {
        // Find ci = max{ s : tile_dst_pgto_beg[s] <= k }, i.e. the CGTO slot
        // whose destination range contains k.
        uint ci = 0u;
        [[unroll]]
        for (uint s = 1u; s <= (uint)TILE_SIZE; ++s) {
            ci = (tile_dst_pgto_beg[s] <= k) ? s : ci;
        }
        PGTO p = args.pgto[tile_src_pgto_beg[ci] + (k - tile_dst_pgto_beg[ci])];
        pgto_coeff_tile[k] = p.coeff;
        pgto_alpha_tile[k] = p.alpha;
        pgto_ijkl_tile[k]  = p.ijkl;
        pgto_slot_tile[k]  = ci;
    }
}

[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid: SV_GroupID, uint3 gtid: SV_GroupThreadID, uint tid: SV_GroupIndex, uint3 dtid: SV_DispatchThreadID) {
    if (tid == 0) {
        num_screened_cgtos = 0;
    }
    GroupMemoryBarrierWithGroupSync();

    bool active = all(dtid < args.grid_dim.xyz);

    float3 model_aabb_min = float3(gid * uint3(WG_X, WG_Y, WG_Z)) * args.step.xyz;
    float3 model_aabb_max = float3((gid + uint3(1, 1, 1)) * uint3(WG_X, WG_Y, WG_Z)) * args.step.xyz;

    // Stage 1: Screening - Prune CGTOs relevant for this workgroup region.
    // Wave-coalesced allocation: one LDS atomic per wave instead of per thread.
    for (uint i = tid; i < args.num_cgtos; i += WG_SIZE) {
        bool passed = false;
        if (i < args.num_cgtos) {
            float r = args.cgto_r[i];
            if (r != 0.0) {
                uint atom_idx = args.cgto_atom_idx[i];
                float3 center = args.atom_xyz[atom_idx].xyz;
                float3 model_xyz = mul(args.world_to_model, float4(center, 1.0)).xyz;
                float r2 = r * r;
                float3 d = clamp(model_xyz, model_aabb_min, model_aabb_max) - model_xyz;
                passed = (dot(d, d) < r2);
            }
        }

        uint wave_count = WaveActiveCountBits(passed);
        uint wave_offset = WavePrefixSum(passed ? 1u : 0u);
        uint wave_base = 0u;
        if (WaveIsFirstLane() && wave_count > 0u) {
            InterlockedAdd(num_screened_cgtos, wave_count, wave_base);
        }
        wave_base = WaveReadLaneFirst(wave_base);

        if (passed) {
            uint idx = wave_base + wave_offset;
            if (idx < MAX_SCREENED) {
                screened_cgtos[idx] = i;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();

    float result = 0.0f;

    uint num_cgtos = min(num_screened_cgtos, (uint)MAX_SCREENED);
    if (num_cgtos > 0) {
        float3 coord = mul(args.index_to_world, float4(float3(dtid), 1.0f)).xyz;
        uint   num_tiles = (num_screened_cgtos + (uint)TILE_SIZE - 1u) / (uint)TILE_SIZE;

        // MO-blocked contraction: reuse phi over a block of MO rows.
        for (uint mo_base = 0u; mo_base < args.num_rows; mo_base += (uint)MO_BLOCK) {
            uint mo_count = min((uint)MO_BLOCK, args.num_rows - mo_base);
            float psi_block[MO_BLOCK] = {0.0f, 0.0f, 0.0f, 0.0f};

            for (uint tile = 0u; tile < num_tiles; ++tile) {
                fill_pgto_tile(tile, num_cgtos, tid);
                GroupMemoryBarrierWithGroupSync();

                // Load MO coefficients for this tile's CGTOs into groupshared.
                uint work_count = mo_count * (uint)TILE_SIZE;
                for (uint work = tid; work < work_count; work += WG_SIZE) {
                    uint mb       = work / (uint)TILE_SIZE;
                    uint ci       = work % (uint)TILE_SIZE;
                    uint cgto_idx = tile_cgto_idx[ci];
                    float c = (cgto_idx != INVALID_IDX)
                            ? args.coeffs[(mo_base + mb) * args.num_cgtos + cgto_idx]
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

            if (active) {
                for (uint mb = 0u; mb < mo_count; ++mb) {
                    float psi = psi_block[mb];
                    result += (args.mode != 0u) ? psi * psi : psi;
                }
            }
        }
    }

    if (active) {
        if ((args.operation & OP_ABS) != 0u) {
            result = abs(result);
        }

        switch (args.operation & OP_ACCUM_MASK) {
            case OP_SET: break; // result is already the new value
            case OP_ADD: result += out_vol[dtid]; break;
            case OP_SUB: result = out_vol[dtid] - result; break;
            case OP_MAX: result = max(result, out_vol[dtid]); break;
            case OP_MIN: result = min(result, out_vol[dtid]); break;
            default: break;  // invalid operation code, treat as OP_SET
        }

        out_vol[dtid] = result;
    }
}

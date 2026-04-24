// critical_points.slang
// Classifies each voxel as a critical point type based on manifold analysis.

#define WG_X 8
#define WG_Y 8
#define WG_Z 8
#define WG_SIZE (WG_X * WG_Y * WG_Z)

#define SG_X (WG_X + 2)
#define SG_Y (WG_Y + 2)
#define SG_Z (WG_Z + 2)
#define SG_SIZE (SG_X * SG_Y * SG_Z)

#define REGULAR 0
#define MAXIMUM 1
#define SPLIT_SADDLE 2
#define MINIMUM 3
#define JOIN_SADDLE 4

#define NUM_NEIGHBORS 26

struct UBO {
    float4x4 index_to_world;
    uint3 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(16)]] // image slot 0
RWTexture3D<float> volumeTex;

[[vk::binding(0)]]
RWStructuredBuffer<uint> ascending;

[[vk::binding(1)]]
RWStructuredBuffer<uint> descending;

[[vk::binding(2)]]
RWStructuredBuffer<int> types;

struct CriticalPointCounts {
    uint maximaCount;
    uint splitSaddleCount;
    uint minimaCount;
    uint joinSaddleCount;
};

[[vk::binding(3)]]
RWStructuredBuffer<CriticalPointCounts> counts; // counts[0]

bool is_valid_voxel(int3 ijk) {
    return all(ijk >= int3(0, 0, 0)) && all(ijk < int3(ubo.dims));
}

float get_scalar(int3 ijk) {
    if (!is_valid_voxel(ijk)) return 0;
    return volumeTex[uint3(ijk)];
}

uint get_voxel_id(int3 ijk) {
    if (!is_valid_voxel(ijk)) return 0xFFFFFFFFu;
    return uint(ijk.x + ijk.y * int(ubo.dims.x) + ijk.z * int(ubo.dims.x) * int(ubo.dims.y));
}

bool sos_is_greater(float value_a, uint id_a, float value_b, uint id_b) {
    return (value_a > value_b) || (value_a == value_b && id_a > id_b);
}

static const int3 neighbor_offset[NUM_NEIGHBORS] = {
    int3(-1, -1, -1), int3(0, -1, -1), int3(1, -1, -1),
    int3(-1, 0, -1),  int3(0, 0, -1),  int3(1, 0, -1),
    int3(-1, 1, -1),  int3(0, 1, -1),  int3(1, 1, -1),

    int3(-1, -1, 0),  int3(0, -1, 0),  int3(1, -1, 0),
    int3(-1, 0, 0),                     int3(1, 0, 0),
    int3(-1, 1, 0),   int3(0, 1, 0),   int3(1, 1, 0),

    int3(-1, -1, 1),  int3(0, -1, 1),  int3(1, -1, 1),
    int3(-1, 0, 1),   int3(0, 0, 1),   int3(1, 0, 1),
    int3(-1, 1, 1),   int3(0, 1, 1),   int3(1, 1, 1),
};

// Precomputed adjacency masks for 26 neighbors
static const uint neighbor_adjacency_mask[NUM_NEIGHBORS] = {
    0x0000161A, 0x00003E3D, 0x00002C32,
    0x0000D6D3, 0x0001FFEF, 0x0001AD96,
    0x0000D098, 0x0001F178, 0x0001A0B0,
    0x0036141B, 0x007E3A3F, 0x006C2436,
    0x01B6C6DB, 0x036D8DB6,
    0x01B090D8, 0x03F171F8, 0x0360A1B0,
    0x00341600, 0x007A3E00, 0x00642C00,
    0x01A6D600, 0x03DFFE00, 0x032DAC00,
    0x0130D000, 0x02F1F000, 0x0161A000,
};

int count_link_components(uint link_mask) {
    if (link_mask == 0u) return 0;

    uint mask = link_mask;
    int num_components = 0;
    for (int start = 0; start < NUM_NEIGHBORS; ++start) {
        uint bit = 1u << start;
        if ((mask & bit) == 0u) continue;
        num_components++;

        // Bitwise flood-fill
        uint component = bit;
        uint prev = 0u;
        while (component != prev) {
            prev = component;
            for (int nb = 0; nb < NUM_NEIGHBORS; ++nb) {
                if ((component & (1u << nb)) != 0u) {
                    component |= neighbor_adjacency_mask[nb] & link_mask;
                }
            }
        }
        mask &= ~component;
    }

    return num_components;
}

groupshared uint numMaxima;
groupshared uint numSplitSaddles;
groupshared uint numMinima;
groupshared uint numJoinSaddles;

groupshared float sharedScalars[SG_X][SG_Y][SG_Z];
groupshared uint sharedIds[SG_X][SG_Y][SG_Z];

[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex, uint3 group_id : SV_GroupID, uint3 local_id_3 : SV_GroupThreadID) {
    if (any(gid >= ubo.dims)) {
        return;
    }
    uint voxel_id = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;

    if (gtid == 0) {
        numMaxima = 0u;
        numSplitSaddles = 0u;
        numMinima = 0u;
        numJoinSaddles = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    int3 local_idx = int3(local_id_3) + int3(1, 1, 1);
    int3 dims_i = int3(ubo.dims);

    // Cooperative linear fill of the shared region
    {
        int tid = int(gtid);
        int3 block_origin = int3(group_id) * int3(WG_X, WG_Y, WG_Z) - int3(1, 1, 1);
        for (int s = tid; s < SG_SIZE; s += WG_SIZE) {
            int sx = s % SG_X;
            int sy = (s / SG_X) % SG_Y;
            int sz = s / (SG_X * SG_Y);

            int3 gp = block_origin + int3(sx, sy, sz);
            bool valid = is_valid_voxel(gp);
            gp = clamp(gp, int3(0, 0, 0), dims_i - int3(1, 1, 1));

            sharedScalars[sx][sy][sz] = volumeTex[uint3(gp)];
            int lin = gp.x + gp.y * dims_i.x + gp.z * dims_i.x * dims_i.y;
            sharedIds[sx][sy][sz] = valid ? uint(lin) : 0xFFFFFFFFu;
        }
    }

    GroupMemoryBarrierWithGroupSync();

    float fi = sharedScalars[local_idx.x][local_idx.y][local_idx.z];
    int type = REGULAR;

    if (ubo.scalar_threshold < fi) {
        uint  neighbor_ids[NUM_NEIGHBORS];
        float neighbor_scalars[NUM_NEIGHBORS];
        for (int ni = 0; ni < NUM_NEIGHBORS; ni++) {
            int3 offset = neighbor_offset[ni];
            int3 n = local_idx + offset;
            neighbor_ids[ni] = sharedIds[n.x][n.y][n.z];
            neighbor_scalars[ni] = sharedScalars[n.x][n.y][n.z];
        }

        // ===== Analyze ASCENDING manifold (maxima + split saddles) =====
        uint reachable_maxima[NUM_NEIGHBORS];
        int reachable_maxima_count = 0;
        uint upper_link_mask = 0u;

        for (int ni = 0; ni < NUM_NEIGHBORS; ni++) {
            uint n_id = neighbor_ids[ni];
            if (n_id == 0xFFFFFFFFu) continue;
            float n_scalar = neighbor_scalars[ni];
            bool in_upper_link = sos_is_greater(n_scalar, n_id, fi, voxel_id);
            if (in_upper_link) {
                upper_link_mask |= (1u << ni);
                uint max_id = ascending[n_id];
                bool found = false;
                for (int j = 0; j < reachable_maxima_count; j++) {
                    if (reachable_maxima[j] == max_id) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    reachable_maxima[reachable_maxima_count++] = max_id;
                }
            }
        }

        if (reachable_maxima_count == 0) {
            type = MAXIMUM;
            uint _orig;
            InterlockedAdd(numMaxima, 1u, _orig);
        } else if (reachable_maxima_count > 1) {
            int num_components = count_link_components(upper_link_mask);
            if (num_components > 1) {
                type = SPLIT_SADDLE;
                uint _orig;
                InterlockedAdd(numSplitSaddles, 1u, _orig);
            }
        }

        // ===== Analyze DESCENDING manifold (minima + join saddles) =====
        if (type == REGULAR) {
            uint reachable_minima[NUM_NEIGHBORS];
            int reachable_minima_count = 0;
            uint lower_link_mask = 0u;

            for (int ni = 0; ni < NUM_NEIGHBORS; ni++) {
                uint n_id = neighbor_ids[ni];
                if (n_id == 0xFFFFFFFFu) continue;
                float n_scalar = neighbor_scalars[ni];
                bool in_lower_link = sos_is_greater(fi, voxel_id, n_scalar, n_id);
                if (in_lower_link) {
                    lower_link_mask |= (1u << ni);
                    uint min_id = descending[n_id];
                    bool found = false;
                    for (int j = 0; j < reachable_minima_count; j++) {
                        if (reachable_minima[j] == min_id) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        reachable_minima[reachable_minima_count++] = min_id;
                    }
                }
            }

            if (reachable_minima_count == 0) {
                type = MINIMUM;
                uint _orig;
                InterlockedAdd(numMinima, 1u, _orig);
            } else {
                int num_components = count_link_components(lower_link_mask);
                if (num_components > 1) {
                    type = JOIN_SADDLE;
                    uint _orig;
                    InterlockedAdd(numJoinSaddles, 1u, _orig);
                }
            }
        }
    }

    types[voxel_id] = type;

    GroupMemoryBarrierWithGroupSync();
    if (gtid == 0) {
        uint _orig;
        InterlockedAdd(counts[0].maximaCount, numMaxima, _orig);
        InterlockedAdd(counts[0].splitSaddleCount, numSplitSaddles, _orig);
        InterlockedAdd(counts[0].minimaCount, numMinima, _orig);
        InterlockedAdd(counts[0].joinSaddleCount, numJoinSaddles, _orig);
    }
}

// bidirectional_manifold.slang
// Computes steepest ascent/descent neighbors for each voxel.

#define WG_X 8
#define WG_Y 8
#define WG_Z 8
#define WG_SIZE (WG_X * WG_Y * WG_Z)

#define SG_X (WG_X + 2)
#define SG_Y (WG_Y + 2)
#define SG_Z (WG_Z + 2)
#define SG_SIZE (SG_X * SG_Y * SG_Z)

#define NUM_NEIGHBORS 26

struct UBO {
    float4x4 index_to_world;
    uint3 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(16)]] // image slot 0 (binding = slot + MD_GPU_MAX_BIND_SLOTS)
RWTexture3D<float> volumeTex;

[[vk::binding(0)]]
RWStructuredBuffer<uint> ascending;

[[vk::binding(1)]]
RWStructuredBuffer<uint> descending;

// Shared memory for scalar values and voxel IDs in a 10x10x10 region
groupshared float sharedScalars[SG_X][SG_Y][SG_Z];
groupshared uint sharedIds[SG_X][SG_Y][SG_Z];

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

[shader("compute")]
[numthreads(WG_X, WG_Y, WG_Z)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex, uint3 group_id : SV_GroupID, uint3 local_id : SV_GroupThreadID) {
    if (any(gid >= ubo.dims)) {
        return;
    }
    uint i = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;
    int3 dims_i = int3(ubo.dims);

    int3 local_idx = int3(local_id) + int3(1, 1, 1);

    // Cooperative linear fill of the shared region
    int tid = int(gtid);
    int3 block_origin = int3(group_id) * int3(WG_X, WG_Y, WG_Z) - int3(1, 1, 1);
    for (int s = tid; s < SG_SIZE; s += WG_SIZE) {
        int sx = s % SG_X;
        int sy = (s / SG_X) % SG_Y;
        int sz = s / (SG_X * SG_Y);

        int3 gp = block_origin + int3(sx, sy, sz);
        gp = clamp(gp, int3(0, 0, 0), dims_i - int3(1, 1, 1));

        sharedScalars[sx][sy][sz] = volumeTex[uint3(gp)];
        int lin = gp.x + gp.y * dims_i.x + gp.z * dims_i.x * dims_i.y;
        sharedIds[sx][sy][sz] = uint(lin);
    }

    GroupMemoryBarrierWithGroupSync();

    float val = sharedScalars[local_idx.x][local_idx.y][local_idx.z];

    float max_val = val;
    float min_val = val;
    uint max_neighbor = i;
    uint min_neighbor = i;

    for (int j = 0; j < NUM_NEIGHBORS; j++) {
        int3 sidx = local_idx + neighbor_offset[j];
        float s = sharedScalars[sidx.x][sidx.y][sidx.z];
        uint n = sharedIds[sidx.x][sidx.y][sidx.z];

        bool is_max = sos_is_greater(s, n, max_val, max_neighbor);
        bool is_min = sos_is_greater(min_val, min_neighbor, s, n);
        max_val      = is_max ? s : max_val;
        max_neighbor = is_max ? n : max_neighbor;
        min_val      = is_min ? s : min_val;
        min_neighbor = is_min ? n : min_neighbor;
    }

    ascending[i]  = max_neighbor;
    descending[i] = min_neighbor;
}

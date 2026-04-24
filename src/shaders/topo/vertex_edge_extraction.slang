// vertex_edge_extraction.slang
// Extracts vertex data and edges from critical points and manifold connectivity.

#define NUM_NEIGHBORS 26
#define MAX_LOCAL_EDGES 2048

struct UBO {
    float4x4 index_to_world;
    uint3 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
RWStructuredBuffer<uint> cp_indices;

struct TypeCounts {
    uint num_maxima;
    uint num_split_saddles;
    uint num_minima;
    uint num_join_saddles;
    uint num_vertices;
};

[[vk::binding(1)]]
RWStructuredBuffer<TypeCounts> type_counts; // type_counts[0]

[[vk::binding(2)]]
RWStructuredBuffer<uint> ascending;

[[vk::binding(3)]]
RWStructuredBuffer<uint> descending;

[[vk::binding(16)]] // image slot 0 (binding = slot + MD_GPU_MAX_BIND_SLOTS)
RWTexture3D<float> volumeTex;

[[vk::binding(4)]]
RWStructuredBuffer<float4> vertex_data; // xyz = position, w = value

[[vk::binding(5)]]
RWStructuredBuffer<uint2> edges; // x = from vertex index, y = to vertex index

struct EdgeCounter {
    uint edge_count;
};

[[vk::binding(6)]]
RWStructuredBuffer<EdgeCounter> edge_counter; // edge_counter[0]

[[vk::binding(7)]]
RWStructuredBuffer<int> voxel_to_vertex_idx;

groupshared uint local_edge_count;
groupshared uint global_edge_offset;
groupshared uint2 local_edges[MAX_LOCAL_EDGES];

int3 get_voxel_coords(uint id) {
    int x = int(id % ubo.dims.x);
    int y = int((id / ubo.dims.x) % ubo.dims.y);
    int z = int(id / (ubo.dims.x * ubo.dims.y));
    return int3(x, y, z);
}

float3 voxel_to_world(int3 ijk) {
    float4 world_pos = mul(ubo.index_to_world, float4(float3(ijk), 1.0));
    return world_pos.xyz;
}

float get_scalar(int3 ijk) {
    return volumeTex[uint3(ijk)];
}

static const int3 neighbors[26] = {
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
[numthreads(64, 1, 1)]
void main(uint3 gid : SV_DispatchThreadID, uint local_id : SV_GroupIndex) {
    uint vert_idx = gid.x;

    // Initialize shared memory
    if (local_id == 0) {
        local_edge_count = 0;
    }
    GroupMemoryBarrierWithGroupSync();

    uint num_vertices = type_counts[0].num_vertices;
    uint num_maxima = type_counts[0].num_maxima;
    uint num_split_saddles = type_counts[0].num_split_saddles;
    uint num_minima = type_counts[0].num_minima;

    if (vert_idx < num_vertices) {
        // === Extract vertex data ===
        uint voxel_id = cp_indices[vert_idx];
        int3 voxel_coords = get_voxel_coords(voxel_id);
        float3 world_pos = voxel_to_world(voxel_coords);
        float value = get_scalar(voxel_coords);

        vertex_data[vert_idx] = float4(world_pos, value);

        uint split_saddle_beg = num_maxima;
        uint split_saddle_end = num_maxima + num_split_saddles;
        uint join_saddle_beg = num_maxima + num_split_saddles + num_minima;
        uint join_saddle_end = num_vertices;

        bool is_split_saddle = (split_saddle_beg <= vert_idx) && (vert_idx < split_saddle_end);
        bool is_join_saddle = (join_saddle_beg <= vert_idx) && (vert_idx < join_saddle_end);

        if (is_split_saddle || is_join_saddle) {
            // Iterate over 26-neighborhood
            uint unique_targets[NUM_NEIGHBORS];
            uint unique_count = 0;
            for (int n = 0; n < NUM_NEIGHBORS; ++n) {
                int3 nc = voxel_coords + neighbors[n];
                if (any(nc < int3(0, 0, 0)) || any(nc >= int3(ubo.dims))) continue;
                uint nid = uint(nc.x) + uint(nc.y) * ubo.dims.x + uint(nc.z) * ubo.dims.x * ubo.dims.y;

                // Follow manifold from neighbor to its CP
                uint cp_id = is_split_saddle ? ascending[nid] : descending[nid];
                if (cp_id == voxel_id) continue;

                int cp_idx_i = voxel_to_vertex_idx[cp_id];
                if (cp_idx_i == -1) continue;
                uint cp_idx = uint(cp_idx_i);

                // Uniqueness check
                bool found = false;
                for (uint ui = 0; ui < unique_count; ++ui) {
                    if (unique_targets[ui] == cp_idx) { found = true; break; }
                }
                if (!found) {
                    unique_targets[unique_count++] = cp_idx;
                }
            }
            // Emit edges to all unique CPs (shared memory batching)
            uint local_offset;
            InterlockedAdd(local_edge_count, unique_count, local_offset);
            uint write_count = min(unique_count, MAX_LOCAL_EDGES - local_offset);
            for (uint ei = 0; ei < write_count; ++ei) {
                local_edges[local_offset + ei] = uint2(vert_idx, unique_targets[ei]);
            }
        }
    }

    GroupMemoryBarrierWithGroupSync();

    uint write_count = min(local_edge_count, MAX_LOCAL_EDGES);
    // One thread reserves space in global buffer
    if (local_id == 0 && write_count > 0) {
        InterlockedAdd(edge_counter[0].edge_count, write_count, global_edge_offset);
    }
    GroupMemoryBarrierWithGroupSync();

    // Threads cooperatively write local edges to global memory
    for (uint i = local_id; i < write_count; i += 64) {
        edges[global_edge_offset + i] = local_edges[i];
    }
}

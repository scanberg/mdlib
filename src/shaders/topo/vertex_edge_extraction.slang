// vertex_edge_extraction.slang
// Extracts vertex data and edges from critical points and manifold connectivity.

#define NUM_NEIGHBORS 26
#define MAX_LOCAL_EDGES 2048

struct MetaData {
    uint32_t vertex_count;
    uint32_t edge_count;
    uint32_t changed_read;
    uint32_t changed_write;
    uint32_t counter;
};

struct RootArgs {
    float4x4 index_to_world;
    uint4 dims;
    float scalar_threshold;

    uint* cp_indices;
    uint* vertex_types;
    float4* vertex_data;
    uint2* edges;
    uint* ascending;
    uint* descending;
    int* voxel_to_vertex_idx;
    MetaData* meta;
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

[[vk::binding(16)]]  // image slot 0 (binding = slot + MD_GPU_MAX_BIND_SLOTS)
RWTexture3D<float> volumeTex;

groupshared uint local_edge_count;
groupshared uint global_edge_offset;
groupshared uint2 local_edges[MAX_LOCAL_EDGES];

int3 get_voxel_coords(uint id) {
    int x = int(id % args.dims.x);
    int y = int((id / args.dims.x) % args.dims.y);
    int z = int(id / (args.dims.x * args.dims.y));
    return int3(x, y, z);
}

float3 voxel_to_world(int3 ijk) {
    float4 world_pos = mul(args.index_to_world, float4(float3(ijk), 1.0));
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
    int3(-1, 0, 0),                    int3(1, 0, 0),
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

    uint num_vertices = args.meta[0].vertex_count;

    if (vert_idx < num_vertices) {
        // === Extract vertex data ===
        uint voxel_id = args.cp_indices[vert_idx];
        int3 voxel_coords = get_voxel_coords(voxel_id);
        float3 world_pos = voxel_to_world(voxel_coords);
        float value = get_scalar(voxel_coords);

        args.vertex_data[vert_idx] = float4(world_pos, value);

        uint type = args.vertex_types[vert_idx];
        bool is_split_saddle = (type == 2u); // SPLIT_SADDLE
        bool is_join_saddle  = (type == 4u); // JOIN_SADDLE

        if (is_split_saddle || is_join_saddle) {
            // Iterate over 26-neighborhood
            uint unique_targets[NUM_NEIGHBORS];
            uint unique_count = 0;
            for (int n = 0; n < NUM_NEIGHBORS; ++n) {
                int3 nc = voxel_coords + neighbors[n];
                if (any(nc < int3(0, 0, 0)) || any(nc >= int3(args.dims.xyz))) continue;
                uint nid = uint(nc.x) + uint(nc.y) * args.dims.x + uint(nc.z) * args.dims.x * args.dims.y;

                // Follow manifold from neighbor to its CP
                uint cp_id = is_split_saddle ? args.ascending[nid] : args.descending[nid];
                if (cp_id == voxel_id) continue;

                int cp_idx_i = args.voxel_to_vertex_idx[cp_id];
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
        InterlockedAdd(args.meta[0].edge_count, write_count, global_edge_offset);
    }
    GroupMemoryBarrierWithGroupSync();

    // Threads cooperatively write local edges to global memory
    for (uint i = local_id; i < write_count; i += 64) {
        args.edges[global_edge_offset + i] = local_edges[i];
    }
}

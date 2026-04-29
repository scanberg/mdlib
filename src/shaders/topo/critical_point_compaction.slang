// critical_point_compaction.slang
// Compacts critical points into a contiguous array, writing global vertex indices.
// Each critical-point voxel atomically claims a slot in cp_indices and records its
// vertex index in voxel_to_vertex_idx (a separate buffer, initialized to -1).
// Non-critical voxels leave voxel_to_vertex_idx untouched (still -1).

#define REGULAR 0
#define MAXIMUM 1
#define SPLIT_SADDLE 2
#define MINIMUM 3
#define JOIN_SADDLE 4

struct UBO {
    float4x4 index_to_world;
    uint3 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
RWStructuredBuffer<int> types;           // read-only after critical_points pass

[[vk::binding(1)]]
RWStructuredBuffer<uint> cp_indices;     // output: voxel_id at each vertex index

struct Counters {
    uint maximaWritePos;
    uint splitSaddleWritePos;
    uint minimaWritePos;
    uint joinSaddleWritePos;
};

[[vk::binding(2)]]
RWStructuredBuffer<Counters> counters;   // counters[0], initialized with per-type base offsets

[[vk::binding(3)]]
RWStructuredBuffer<int> voxel_to_vertex_idx; // output: vertex index for each CP voxel (-1 = not a CP)

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID) {
    uint num_points = ubo.dims.x * ubo.dims.y * ubo.dims.z;
    uint voxel_id = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;
    if (voxel_id >= num_points) return;

    int type = types[voxel_id];
    if (type == REGULAR) return;

    // Claim a slot in the appropriate section of cp_indices via atomic.
    uint vtx_idx;
    if (type == MAXIMUM) {
        InterlockedAdd(counters[0].maximaWritePos, 1u, vtx_idx);
    } else if (type == SPLIT_SADDLE) {
        InterlockedAdd(counters[0].splitSaddleWritePos, 1u, vtx_idx);
    } else if (type == MINIMUM) {
        InterlockedAdd(counters[0].minimaWritePos, 1u, vtx_idx);
    } else { // JOIN_SADDLE
        InterlockedAdd(counters[0].joinSaddleWritePos, 1u, vtx_idx);
    }

    cp_indices[vtx_idx] = voxel_id;
    voxel_to_vertex_idx[voxel_id] = int(vtx_idx);
}


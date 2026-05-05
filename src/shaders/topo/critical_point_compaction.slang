// critical_point_compaction.slang
// Compacts critical points into a contiguous array, writing global vertex indices.
// Each critical-point voxel atomically claims a slot via a single counter and records
// its vertex index in voxel_to_vertex_idx. Vertex types are written to vertex_types[].
// Non-critical voxels leave voxel_to_vertex_idx untouched (still -1).

#define REGULAR      0
#define MAXIMUM      1
#define SPLIT_SADDLE 2
#define MINIMUM      3
#define JOIN_SADDLE  4

struct UBO {
    float4x4 index_to_world;
    uint4 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
RWStructuredBuffer<int> types;           // read-only after critical_points pass

[[vk::binding(1)]]
RWStructuredBuffer<uint> cp_indices;     // output: voxel_id at each vertex index

[[vk::binding(2)]]
RWStructuredBuffer<uint> counter;        // counter[0]: single atomic write cursor

[[vk::binding(3)]]
RWStructuredBuffer<int> voxel_to_vertex_idx; // output: vertex index for each CP voxel (-1 = not a CP)

[[vk::binding(4)]]
RWStructuredBuffer<uint> vertex_types;   // output: per-vertex critical point type

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID) {
    uint num_points = ubo.dims.x * ubo.dims.y * ubo.dims.z;
    uint voxel_id = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;
    if (voxel_id >= num_points) return;

    int type = types[voxel_id];
    if (type == REGULAR) return;

    uint vtx_idx;
    InterlockedAdd(counter[0], 1u, vtx_idx);

    cp_indices[vtx_idx]           = voxel_id;
    vertex_types[vtx_idx]         = uint(type);
    voxel_to_vertex_idx[voxel_id] = int(vtx_idx);
}


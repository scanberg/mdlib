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

    int* types;
    uint* cp_indices;
    MetaData* meta;
    int* voxel_to_vertex_idx;
    uint* vertex_types;
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID) {
    uint num_points = args.dims.x * args.dims.y * args.dims.z;
    uint voxel_id = gid.x + gid.y * args.dims.x + gid.z * args.dims.x * args.dims.y;
    if (voxel_id >= num_points) return;

    int type = args.types[voxel_id];
    if (type == REGULAR) return;

    uint vtx_idx;
    InterlockedAdd(args.meta[0].counter, 1u, vtx_idx);

    args.cp_indices[vtx_idx]           = voxel_id;
    args.vertex_types[vtx_idx]         = uint(type);
    args.voxel_to_vertex_idx[voxel_id] = int(vtx_idx);
}


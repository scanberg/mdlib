// path_compression.slang
// Pointer-jumping pass for ascending and descending manifolds.

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

    uint* ascending;
    uint* descending;
    MetaData* meta;
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

groupshared uint local_changed;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex) {
    // Early exit: previous iteration changed nothing → already converged.
    if (args.meta[0].changed_read == 0u) {
        return;
    }

    bool active = all(gid < args.dims.xyz);
    uint i = gid.x + gid.y * args.dims.x + gid.z * args.dims.x * args.dims.y;
    uint parent_asc  = active ? args.ascending[i] : i;
    uint parent_desc = active ? args.descending[i] : i;

    if (gtid == 0) {
        local_changed = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    uint change = 0u;

    // Single-step pointer jump for ascending manifold
    if (active && parent_asc != i) {
        uint grandparent_asc = args.ascending[parent_asc];
        if (parent_asc != grandparent_asc) {
            args.ascending[i] = grandparent_asc;
            change = 1u;
        }
    }

    // Single-step pointer jump for descending manifold
    if (active && parent_desc != i) {
        uint grandparent_desc = args.descending[parent_desc];
        if (parent_desc != grandparent_desc) {
            args.descending[i] = grandparent_desc;
            change = 1u;
        }
    }

    InterlockedOr(local_changed, change);

    GroupMemoryBarrierWithGroupSync();

    if (gtid == 0) {
        InterlockedOr(args.meta[0].changed_write, local_changed);
    }
}

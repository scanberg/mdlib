// path_compression.slang
// Pointer-jumping pass for ascending and descending manifolds.

struct UBO {
    float4x4 index_to_world;
    uint4 dims;
    float scalar_threshold;
};

struct MetaData {
    uint32_t vertex_count;
    uint32_t edge_count;
    uint32_t changed_read;
    uint32_t changed_write;
    uint32_t counter;       
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
RWStructuredBuffer<uint> ascending;

[[vk::binding(1)]]
RWStructuredBuffer<uint> descending;

// changed_read[0]: non-zero if the previous iteration made any changes.
// If zero, all paths are already compressed — this dispatch is a no-op.
[[vk::binding(2)]]
RWStructuredBuffer<MetaData> meta;

groupshared uint local_changed;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex) {
    // Early exit: previous iteration changed nothing → already converged.
    if (meta[0].changed_read == 0u) {
        return;
    }

    if (any(gid >= ubo.dims.xyz)) {
        return;
    }

    uint i = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;
    uint parent_asc  = ascending[i];
    uint parent_desc = descending[i];

    if (gtid == 0) {
        local_changed = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    uint change = 0u;

    // Single-step pointer jump for ascending manifold
    if (parent_asc != i) {
        uint grandparent_asc = ascending[parent_asc];
        if (parent_asc != grandparent_asc) {
            ascending[i] = grandparent_asc;
            change = 1u;
        }
    }

    // Single-step pointer jump for descending manifold
    if (parent_desc != i) {
        uint grandparent_desc = descending[parent_desc];
        if (parent_desc != grandparent_desc) {
            descending[i] = grandparent_desc;
            change = 1u;
        }
    }

    InterlockedOr(local_changed, change);

    GroupMemoryBarrierWithGroupSync();

    if (gtid == 0) {
        InterlockedOr(meta[0].changed_write, local_changed);
    }
}

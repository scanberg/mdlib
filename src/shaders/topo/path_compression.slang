// path_compression.slang
// Pointer-jumping pass for ascending and descending manifolds.

struct UBO {
    float4x4 index_to_world;
    uint4 dims;
    float scalar_threshold;
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
RWStructuredBuffer<uint> changed_read;

// changed_write[0]: set non-zero by this iteration if any path step was taken.
// Caller clears this to 0 before each dispatch, then copies it to changed_read
// after the dispatch so the next iteration can check it.
[[vk::binding(3)]]
RWStructuredBuffer<uint> changed_write;

groupshared uint local_changed;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex) {
    // Early exit: previous iteration changed nothing → already converged.
    if (changed_read[0] == 0u) {
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

    uint _orig;
    InterlockedOr(local_changed, change, _orig);

    GroupMemoryBarrierWithGroupSync();

    if (gtid == 0) {
        InterlockedOr(changed_write[0], local_changed, _orig);
    }
}

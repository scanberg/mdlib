// path_compression.slang
// Pointer-jumping pass for ascending and descending manifolds.

struct UBO {
    float4x4 index_to_world;
    uint3 dims;
    float scalar_threshold;
};

[[vk::push_constant]]
ConstantBuffer<UBO> ubo;

[[vk::binding(0)]]
RWStructuredBuffer<uint> ascending;

[[vk::binding(1)]]
RWStructuredBuffer<uint> descending;

[[vk::binding(2)]]
RWStructuredBuffer<uint> changed_flag; // changed_flag[0]

groupshared uint local_changed;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID, uint gtid : SV_GroupIndex) {
    if (any(gid >= ubo.dims)) {
        return;
    }

    uint i = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;

    if (gtid == 0) {
        local_changed = 0u;
    }
    GroupMemoryBarrierWithGroupSync();

    // Single-step pointer jump for ascending manifold
    uint parent_asc = ascending[i];
    if (parent_asc != i) {
        uint grandparent_asc = ascending[parent_asc];
        if (parent_asc != grandparent_asc) {
            ascending[i] = grandparent_asc;
            uint _orig;
            InterlockedOr(local_changed, 1u, _orig);
        }
    }

    // Single-step pointer jump for descending manifold
    uint parent_desc = descending[i];
    if (parent_desc != i) {
        uint grandparent_desc = descending[parent_desc];
        if (parent_desc != grandparent_desc) {
            descending[i] = grandparent_desc;
            uint _orig;
            InterlockedOr(local_changed, 1u, _orig);
        }
    }

    GroupMemoryBarrierWithGroupSync();

    if (gtid == 0) {
        uint _orig;
        InterlockedOr(changed_flag[0], local_changed, _orig);
    }
}

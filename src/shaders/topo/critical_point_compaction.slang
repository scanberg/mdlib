// critical_point_compaction.slang
// Compacts critical points into a contiguous array, writing global vertex indices.

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
RWStructuredBuffer<int> types;

[[vk::binding(1)]]
RWStructuredBuffer<uint> cp_indices;

struct Counters {
    uint maximaWritePos;
    uint splitSaddleWritePos;
    uint minimaWritePos;
    uint joinSaddleWritePos;
};

[[vk::binding(2)]]
RWStructuredBuffer<Counters> counters; // counters[0]

// Shared memory for workgroup-local compaction
groupshared uint localMaxima[64];
groupshared uint localSplitSaddles[64];
groupshared uint localMinima[64];
groupshared uint localJoinSaddles[64];

groupshared uint localMaximaCount;
groupshared uint localSplitSaddleCount;
groupshared uint localMinimaCount;
groupshared uint localJoinSaddleCount;

groupshared uint globalMaximaOffset;
groupshared uint globalSplitSaddleOffset;
groupshared uint globalMinimaOffset;
groupshared uint globalJoinSaddleOffset;

[shader("compute")]
[numthreads(8, 8, 8)]
void main(uint3 gid : SV_DispatchThreadID, uint local_id : SV_GroupIndex) {
    uint voxel_id = gid.x + gid.y * ubo.dims.x + gid.z * ubo.dims.x * ubo.dims.y;

    // Initialize shared counters
    if (local_id == 0) {
        localMaximaCount = 0;
        localSplitSaddleCount = 0;
        localMinimaCount = 0;
        localJoinSaddleCount = 0;
    }
    GroupMemoryBarrierWithGroupSync();

    uint num_points = ubo.dims.x * ubo.dims.y * ubo.dims.z;

    // Step 1: Each thread checks its voxel and compacts locally
    int type = (voxel_id < num_points) ? types[voxel_id] : REGULAR;

    if (type == MAXIMUM) {
        uint pos;
        InterlockedAdd(localMaximaCount, 1u, pos);
        localMaxima[pos] = voxel_id;
    } else if (type == SPLIT_SADDLE) {
        uint pos;
        InterlockedAdd(localSplitSaddleCount, 1u, pos);
        localSplitSaddles[pos] = voxel_id;
    } else if (type == MINIMUM) {
        uint pos;
        InterlockedAdd(localMinimaCount, 1u, pos);
        localMinima[pos] = voxel_id;
    } else if (type == JOIN_SADDLE) {
        uint pos;
        InterlockedAdd(localJoinSaddleCount, 1u, pos);
        localJoinSaddles[pos] = voxel_id;
    }

    GroupMemoryBarrierWithGroupSync();

    // Step 2: One thread per workgroup reserves space in global memory
    if (local_id == 0) {
        if (localMaximaCount > 0)
            InterlockedAdd(counters[0].maximaWritePos, localMaximaCount, globalMaximaOffset);
        if (localSplitSaddleCount > 0)
            InterlockedAdd(counters[0].splitSaddleWritePos, localSplitSaddleCount, globalSplitSaddleOffset);
        if (localMinimaCount > 0)
            InterlockedAdd(counters[0].minimaWritePos, localMinimaCount, globalMinimaOffset);
        if (localJoinSaddleCount > 0)
            InterlockedAdd(counters[0].joinSaddleWritePos, localJoinSaddleCount, globalJoinSaddleOffset);
    }
    GroupMemoryBarrierWithGroupSync();

    // Step 3: Cooperatively write to global memory and overwrite types with vertex indices
    if (local_id < localMaximaCount) {
        uint vtx_idx = globalMaximaOffset + local_id;
        uint vtx_id = localMaxima[local_id];
        cp_indices[vtx_idx] = vtx_id;
        types[vtx_id] = int(vtx_idx);
    }
    if (local_id < localSplitSaddleCount) {
        uint vtx_idx = globalSplitSaddleOffset + local_id;
        uint vtx_id = localSplitSaddles[local_id];
        cp_indices[vtx_idx] = vtx_id;
        types[vtx_id] = int(vtx_idx);
    }
    if (local_id < localMinimaCount) {
        uint vtx_idx = globalMinimaOffset + local_id;
        uint vtx_id = localMinima[local_id];
        cp_indices[vtx_idx] = vtx_id;
        types[vtx_id] = int(vtx_idx);
    }
    if (local_id < localJoinSaddleCount) {
        uint vtx_idx = globalJoinSaddleOffset + local_id;
        uint vtx_id = localJoinSaddles[local_id];
        cp_indices[vtx_idx] = vtx_id;
        types[vtx_id] = int(vtx_idx);
    }
    // For non-critical points, set to -1
    if (type == REGULAR && voxel_id < num_points) {
        types[voxel_id] = -1;
    }
}

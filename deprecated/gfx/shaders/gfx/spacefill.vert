#version 460

#include "common.h"

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = CLUSTER_BOUNDS_BINDING) readonly buffer ClusterBoundsBuffer {
    ClusterBounds cluster_bounds[];
};

layout(binding = CLUSTER_POS_RAD_BINDING) readonly buffer ClusterPosRadBuffer {
    CompressedPosRad cluster_pos_rad[];
};

layout(binding = CLUSTER_INST_DATA_BINDING) readonly buffer VisibleClusterBuffer {
    ClusterInstance cluster_instances[];
};

layout(binding = TRANSFORM_BINDING) readonly buffer TransformBuffer {
    mat4 transforms[];
};

layout(location = 0) out VS_GS {
    flat vec4 view_sphere;
    flat uint index;
} OUT;

void main() {
    uint cluster_atom_idx = gl_VertexID & 0xFFu;
    uint cluster_inst_idx = gl_VertexID >> 8;

    ClusterInstance inst = cluster_instances[cluster_inst_idx];

    uint cluster_range = inst.cluster_range;
    uint cluster_idx   = inst.cluster_idx;
    uint transform_idx = inst.transform_idx;
    float radius_scale = 1.0;

    mat4 model_mat = transforms[transform_idx];
    uint cluster_offset = cluster_range >> 8;
    uint cluster_size   = cluster_range & 0xFF;

    ClusterBounds bounds = cluster_bounds[cluster_idx];
    CompressedPosRad ccoords = cluster_pos_rad[cluster_offset + cluster_atom_idx];

    vec4 ncoords = vec4(unpackUnorm2x16(ccoords.xy), unpackUnorm2x16(ccoords.zr));
    vec4 coords = mix(bounds.cluster_min, bounds.cluster_max, ncoords);

    vec3  c = coords.xyz;
    float r = coords.w * radius_scale;

    vec4 view_pos = ubo.world_to_view * model_mat * vec4(c,1);
    vec4 view_sphere = vec4(view_pos.xyz, r);

	OUT.view_sphere = view_sphere;
	OUT.index = gl_VertexID;
}
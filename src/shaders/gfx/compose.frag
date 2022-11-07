#version 460

#extension GL_ARB_gpu_shader_int64 : enable

#include "common.h"

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = TRANSFORM_BINDING) readonly buffer TransformBuffer {
    mat4 in_transforms[];
};

layout (binding = VIS_BINDING, std430) readonly buffer VisBuffer {
	uint64_t in_vis_buf[];
};

layout (binding = CLUSTER_INST_DATA_BINDING, std430) readonly buffer ClusterInstanceBuffer {
	ClusterInstance in_cluster_inst[];
};

layout (binding = CLUSTER_BOUNDS_BINDING, std430) readonly buffer ClusterBoundsBuffer {
	ClusterBounds in_cluster_bounds[];
};

layout(binding = CLUSTER_POS_RAD_BINDING) readonly buffer ClusterPosRadBuffer {
    CompressedPosRad in_cluster_pos_rad[];
};

layout(location = 0) in vec2 uv;

// These are the outputs of the external FBO, this needs to be configured in some way
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec2 out_normal;
layout(location = 2) out vec2 out_velocity;
layout(location = 3) out vec4 out_index;

uint hash(uint a) {
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
}

vec2 encode_normal (vec3 n) {
   float p = sqrt(n.z * 8 + 8);
   return n.xy / p + 0.5;
}

void main() {
	uint64_t vis_data = in_vis_buf[uint(gl_FragCoord.y) * ubo.render_width + uint(gl_FragCoord.x)];
	float depth  = uintBitsToFloat(uint(vis_data >> 32));
	uint payload = uint(vis_data);
	if (depth == 1) discard;

	uint clust_inst_idx = payload >> 8;
	uint clust_atom_idx = payload & 0xFF;

	ClusterInstance clust_inst = in_cluster_inst[clust_inst_idx];
	uint clust_idx = clust_inst.cluster_idx;
	uint clust_range = clust_inst.cluster_range;
	uint clust_offset = clust_range >> 8;
	uint clust_size	= clust_range & 0xFF;
	uint transform_idx = clust_inst.transform_idx;
	vec4 clust_color = unpackUnorm4x8(clust_inst.color);

	mat4 model_mat = in_transforms[transform_idx];
	ClusterBounds clust_bounds = in_cluster_bounds[clust_idx];
	CompressedPosRad ccoords = in_cluster_pos_rad[clust_offset + clust_atom_idx];

	vec4 ncoords = vec4(unpackUnorm2x16(ccoords.xy), unpackUnorm2x16(ccoords.zr));
	vec4 coords = mix(clust_bounds.cluster_min, clust_bounds.cluster_max, ncoords);

	vec4 view_pos = ubo.world_to_view * model_mat * vec4(coords.xyz, 1);
    vec4 view_sphere = vec4(view_pos.xyz, coords.w);

	vec4 clip_coord = vec4((gl_FragCoord.xy / vec2(ubo.render_width, ubo.render_height)) * 2.0 - 1.0, depth * 2.0 - 1.0, 1);
	vec4 view_coord = ubo.clip_to_view * clip_coord;
	view_coord.xyz /= view_coord.w;

	vec3 view_normal = normalize(view_coord.xyz - view_sphere.xyz);

	uint chash = hash(transform_idx);
	vec3 ccolor = vec3(float(chash & 255), float((chash >> 8) & 255), float((chash >> 16) & 255)) / 255.0;

	gl_FragDepth = depth;
	out_color    = vec4(ccolor, 1);
	out_normal   = encode_normal(view_normal);
	out_velocity = vec2(0,0);
	out_index	 = vec4(0,0,0,0);
}
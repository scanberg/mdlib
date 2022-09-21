#version 460

#include "common.h"

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = POSITION_BINDING, std430) readonly buffer PositionBuffer {
    Position atom_pos[];
};

layout(binding = RADIUS_BINDING) readonly buffer RadiusBuffer {
    float atom_rad[];
};

layout(binding = COLOR_BINDING) readonly buffer ColorBuffer {
    uint colors[];
};

layout(binding = TRANSFORM_BINDING) readonly buffer TransformBuffer {
    mat4 transforms[];
};

layout(binding = DRAW_INDIRECT_BINDING) readonly buffer DrawIndirectBuffer {
    DrawIndirect indirect;
};

layout(binding = DRAW_SPHERE_INDEX_BINDING) readonly buffer SphereIndexBuffer {
    uint indices[];
};

layout(location = 0) out VS_GS {
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat uint color;
    flat uint index;
} OUT;

void main() {
    uint idx = indices[gl_VertexID];
    uint draw_idx  = gl_DrawID;
    uint atom_idx  = indirect.draw_sphere_cmd[draw_idx].atom_offset + idx;
    uint color_idx = indirect.draw_sphere_cmd[draw_idx].color_offset + idx;
    uint transform_idx = indirect.draw_sphere_cmd[draw_idx].transform_idx;
    float radius_scale = indirect.draw_sphere_cmd[draw_idx].radius_scale;

    mat4 model_mat = transforms[transform_idx];
    uint color = colors[color_idx];

    vec4  pos = vec4(atom_pos[atom_idx].x, atom_pos[atom_idx].y, atom_pos[atom_idx].z, 1);
    float rad = atom_rad[atom_idx] * radius_scale;
    vec4 view_pos = ubo.world_to_view * model_mat * pos;
    vec4 view_sphere = vec4(view_pos.xyz, rad);

	OUT.view_sphere = view_sphere;
    OUT.view_velocity = vec3(0,0,0);
	OUT.color = color;
	OUT.index = atom_idx;
}
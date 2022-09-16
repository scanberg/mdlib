#version 460

struct UniformData {
    mat4 world_to_view;
    mat4 world_to_view_normal;
    mat4 world_to_clip;
    mat4 view_to_clip;
    mat4 view_to_world;
    mat4 clip_to_view;
    mat4 prev_world_to_clip;
    mat4 curr_view_to_prev_clip;
    vec4 jitter_uv;
    uint atom_mask;
    uint ortho;
    vec2 fle;
    float radius_scale;
};

struct DrawData {
    uint transform_offset;
    uint color_offset;
};

struct VertexData {
    smooth vec3 view_coord;
    flat   vec4 view_sphere;
    flat   vec3 view_velocity;
    flat   vec4 color;
    flat   uint atom_idx;
    smooth vec2 uv;
};

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = 1) buffer readonly PositionBuffer {
    vec3 atom_pos[];
};

layout(binding = 2) buffer readonly RadiusBuffer {
    float atom_rad[];
}

layout(binding = 3) buffer readonly ColorBuffer {
    uint colors[];
}

layout(binding = 4) buffer readonly TransformBuffer {
    mat4 model_transforms[];
}

layout(binding = 5) buffer readonly DrawDataBuffer {
    DrawData draw_data[];
}

layout(location = 0) out VertexData out_vert;

void main() {
    uint vertex_idx = gl_VertexID;
    uint sphere_idx = gl_InstanceID;
    uint draw_idx = gl_DrawID;

    mat4 mdl_mat = model_mat[draw_data[draw_idx].transform_offset];
    vec4 color   = colors[draw_data[draw_idx].color_offset];

    vec4  pos = vec4(atom_pos[sphere_idx], 1);
    float rad = atom_rad[sphere_idx] * u_radius_scale;
    vec4 view_pos = u_world_to_view * model_mat * pos;
    vec4 view_sphere = vec4(view_pos, rad);
    vec4 view_coord = vec4(view_pos.xyz + vec3((vertex_idx & 1U) ? -rad : +rad, (vertex_idx & 2U) ? -rad + rad, +rad), 1);
    
    out_vert.view_coord  = view_coord
	out_vert.view_sphere = view_sphere;
	out_vert.color = color;
	out_vert.atom_idx = uint(gl_VertexID);

    gl_Position = ubo.view_to_clip * view_coord;
}
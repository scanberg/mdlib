#version 460

struct UniformData {
    mat4 world_to_view;
    mat4 world_to_view_normal;
    mat4 world_to_clip;
    mat4 view_to_clip;
    vec2 fle;
    uint ortho;
    uint atom_mask;
};

struct DrawData {
    uint transform_offset;
    uint color_offset;
};

layout (binding = 0, std140) uniform UboBuffer {
    UniformData ubo;
};

layout(binding = 1) readonly buffer PositionBuffer {
    vec3 atom_pos[];
};

layout(binding = 2) readonly buffer RadiusBuffer {
    float atom_rad[];
};

layout(binding = 3) readonly buffer ColorBuffer {
    uint colors[];
};

layout(binding = 4) readonly buffer TransformBuffer {
    mat4 transforms[];
};

layout(binding = 5) readonly buffer DrawDataBuffer {
    DrawData draw_data[];
};

layout(location = 0) out VertexData {
    smooth vec3 view_coord;
    flat   vec4 view_sphere;
    flat   vec3 view_velocity;
    flat   vec4 color;
    flat   uint atom_idx;
    smooth vec2 uv;
} out_vert;

void main() {
    uint vertex_idx = gl_VertexID;
    uint sphere_idx = gl_InstanceID;
    uint draw_idx = gl_DrawID;

    mat4 model_mat = transforms[draw_data[draw_idx].transform_offset];
    vec4 color   = unpackSnorm4x8(colors[draw_data[draw_idx].color_offset]);

    vec4  pos = vec4(atom_pos[sphere_idx], 1);
    float rad = atom_rad[sphere_idx] * ubo.radius_scale;
    vec4 view_pos = ubo.world_to_view * model_mat * pos;
    vec4 view_sphere = vec4(view_pos.xyz, rad);
    vec4 view_coord = vec4(view_pos.xyz + vec3(((vertex_idx & 1U) != 0U) ? -rad : +rad, ((vertex_idx & 2U) != 0U) ? -rad : +rad, +rad), 1);
    
    out_vert.view_coord  = view_coord.xyz;
	out_vert.view_sphere = view_sphere;
	out_vert.color = color;
	out_vert.atom_idx = uint(gl_VertexID);

    gl_Position = ubo.view_to_clip * view_coord;
}
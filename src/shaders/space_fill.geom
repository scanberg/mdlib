#version 330 core

#ifndef ORTHO
#define ORTHO 0
#endif

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
} in_vert[];

out GS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
    flat vec3 view_velocity;
    flat vec4 color;
    flat uint atom_idx;
    smooth vec2 uv;
} out_frag;

layout (std140) uniform ubo {
    mat4 u_world_to_view;
    mat4 u_world_to_view_normal;
    mat4 u_world_to_clip;
    mat4 u_view_to_clip;
    mat4 u_view_to_world;
    mat4 u_clip_to_view;
    mat4 u_prev_world_to_clip;
    mat4 u_curr_view_to_prev_clip;
    vec4 u_jitter_uv;
    uint u_atom_mask;
    uint _pad0;
    uint _pad1;
    uint _pad2;
    float u_radius_scale;
};

void emit_vertex_ortho(vec2 uv) {
    vec3  pos = in_vert[0].view_sphere.xyz;
    float rad = in_vert[0].view_sphere.w;
    out_frag.uv = uv;
    out_frag.view_coord = vec3(pos + vec3(uv * rad, 0));
    gl_Position = u_view_to_clip * vec4(out_frag.view_coord + vec3(0,0,rad), 1.0);
    EmitVertex();
}

void emit_vertex_persp(vec4 clip_coord) {
    vec4 view_coord = u_clip_to_view * clip_coord;
    out_frag.view_coord = view_coord.xyz / view_coord.w;
    gl_Position = clip_coord;
    EmitVertex();
}

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
                 in vec2 fle,
                 out vec2 axis_a,
                 out vec2 axis_b,
                 out vec2 center) {
    vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
    float z2 = o.z*o.z; 
    float l2 = dot(o,o);
    float c = -r2*(r2-l2)/((l2-z2)*(r2-z2));
    
    // axis
    axis_a = fle*sqrt(c/(r2-z2)) * vec2( o.x,o.y);
    axis_b = fle*sqrt(c/(r2-l2)) * vec2(-o.y,o.x);
    center = -fle*o.z*o.xy/(z2-r2);
}

void main()
{
    if (in_vert[0].view_sphere.w == 0) {
        EndPrimitive();
        return;
    }
    if (in_vert[0].color.a == 0) {
        EndPrimitive();
        return;
    }

    out_frag.view_sphere = in_vert[0].view_sphere;
    out_frag.view_velocity = in_vert[0].view_velocity;
    out_frag.color = in_vert[0].color;
    out_frag.atom_idx = in_vert[0].atom_idx;

#if ORTHO
    emit_vertex_ortho(vec2(-1, -1));
    emit_vertex_ortho(vec2(+1, -1));
    emit_vertex_ortho(vec2(-1, +1));
    emit_vertex_ortho(vec2(+1, +1));
    EndPrimitive();
#else
    // Focal length
    float clip_z = -u_view_to_clip[2][2] - u_view_to_clip[3][2] / (in_vert[0].view_sphere.z + in_vert[0].view_sphere.w);
    vec2 fle = vec2(u_view_to_clip[0][0], u_view_to_clip[1][1]);
    vec2 axis_a;
    vec2 axis_b;
    vec2 center;
    proj_sphere(in_vert[0].view_sphere, fle, axis_a, axis_b, center);

    emit_vertex_persp(vec4(center -axis_a -axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center +axis_a -axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center -axis_a +axis_b, clip_z, 1));
    emit_vertex_persp(vec4(center +axis_a +axis_b, clip_z, 1));
    EndPrimitive();
#endif
}
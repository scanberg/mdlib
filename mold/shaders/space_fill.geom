#version 330 core

#ifndef ATOM_VEL
#define ATOM_VEL 1
#endif

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in VS_GS {
    flat vec4 view_sphere;
#if ATOM_VEL
    flat vec3 view_velocity;
#endif
#if ATOM_COL
    flat vec4 color;
#endif
#if ATOM_IDX
    flat uint atom_idx;
#endif
} in_vert[];

out GS_FS {
    smooth vec3 view_coord;
    flat vec4 view_sphere;
#if ATOM_VEL
    flat vec3 view_velocity;
#endif
#if ATOM_COL
    flat vec4 color;
#endif
#if ATOM_IDX
    flat uint atom_idx;
#endif
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
    
    float u_radius_scale;
};

void emit_vertex(vec4 clip_coord) {
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
    if (in_vert[0].color.a == 0 || in_vert[0].view_sphere.w == 0) { // Zero alpha or zero radius
        EndPrimitive();
        return;
    }

    // Focal length
    vec2 fle = vec2(u_view_to_clip[0][0], u_view_to_clip[1][1]);

    vec2 axis_a;
    vec2 axis_b;
    vec2 center;
    proj_sphere(in_vert[0].view_sphere, fle, axis_a, axis_b, center);
    center += u_jitter_uv.xy * 0.5;

    // Compute clip z (with bias of radius) 
    float clip_z = -u_view_to_clip[2][2] - u_view_to_clip[3][2] / (in_vert[0].view_sphere.z + in_vert[0].view_sphere.w);

    out_frag.view_sphere = in_vert[0].view_sphere;
#if ATOM_VEL
    out_frag.view_velocity = in_vert[0].view_velocity;
#endif
#if ATOM_COL
    out_frag.color = in_vert[0].color;
#endif
#if ATOM_IDX
    out_frag.atom_idx = in_vert[0].atom_idx;
#endif

    emit_vertex(vec4(center -axis_a -axis_b, clip_z, 1));
    emit_vertex(vec4(center +axis_a -axis_b, clip_z, 1));
    emit_vertex(vec4(center -axis_a +axis_b, clip_z, 1));
    emit_vertex(vec4(center +axis_a +axis_b, clip_z, 1));

    EndPrimitive();
}
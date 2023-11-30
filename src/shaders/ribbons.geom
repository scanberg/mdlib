#version 330 core

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
    float u_width_scale;
    float u_height_scale;
};

layout(lines) in;
layout(triangle_strip, max_vertices = 24) out;

in Vertex {
    vec3 control_point;
    vec3 support_vector;
    vec3 support_tangent;
    vec3 view_velocity;
    vec4 color;
    uint picking_idx;
    uint spline_flags;
    uint atom_flags;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec3 view_coord;
    smooth vec3 view_velocity;
    smooth vec3 view_normal;
    flat   uint picking_idx;
} out_frag;

void emit_vertex(in vec4 clip_coord, in vec4 normal, in int idx) {
    vec4 view_coord = u_clip_to_view * clip_coord;
    out_frag.color = in_vert[idx].color;
    out_frag.view_coord = view_coord.xyz / view_coord.w;
    out_frag.view_velocity = in_vert[idx].view_velocity;
    out_frag.view_normal = normal.xyz;
    out_frag.picking_idx = in_vert[idx].picking_idx;
    gl_Position = clip_coord;
    EmitVertex();
}

void main() {
    // We consider the rendered segment to fully belong to index 0
    if ((in_vert[0].atom_flags & u_atom_mask) != u_atom_mask) return;
    if (in_vert[0].color.a == 0.0f) return;

    vec4 p[2];
    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    uint f[2];
    
    float flip = sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    p[0] = vec4(in_vert[0].control_point, 1);
    p[1] = vec4(in_vert[1].control_point, 1);
    x[0] = vec4(in_vert[0].support_vector, 0);
    x[1] = vec4(in_vert[1].support_vector * flip, 0);
    z[0] = vec4(in_vert[0].support_tangent, 0);
    z[1] = vec4(in_vert[1].support_tangent, 0);
    y[0] = vec4(normalize(cross(z[0].xyz, x[0].xyz)), 0); // To maintain right-handedness
    y[1] = vec4(normalize(cross(z[1].xyz, x[1].xyz)), 0);
    f[0] = in_vert[0].spline_flags;
    f[1] = in_vert[1].spline_flags;

#if 0
    // It is possible to fix tesselation direction so it follows the twist of the spline
    // By tesselating with respect to the twist we would avoid some thin long triangles.
    // This implementation is however currently bogus
    //float flip_sign = sign(dot(x[1].xyz, y[0].xyz));
    float flip_sign = -1;
    x[0].xyz *= flip_sign;
    x[1].xyz *= flip_sign;
    y[0].xyz *= flip_sign;
    y[1].xyz *= flip_sign;
    z[0].xyz *= flip_sign;
    z[1].xyz *= flip_sign;
#endif

    vec2 scale = vec2(u_width_scale, u_height_scale);

    mat4 M[2];
    mat4 N[2];
    vec4 clip_coord0[4];
    vec4 clip_coord1[4];

    M[0] = u_world_to_clip * mat4(x[0], y[0], z[0], p[0]);
    M[1] = u_world_to_clip * mat4(x[1], y[1], z[1], p[1]);

    N[0] = u_world_to_view_normal * mat4(x[0], y[0], z[0], p[0]);
    N[1] = u_world_to_view_normal * mat4(x[1], y[1], z[1], p[1]);

    clip_coord0[0] = M[0] * vec4(vec2(-1,-1) * scale, 0, 1);
    clip_coord0[1] = M[0] * vec4(vec2( 1,-1) * scale, 0, 1);
    clip_coord0[2] = M[0] * vec4(vec2(-1, 1) * scale, 0, 1);
    clip_coord0[3] = M[0] * vec4(vec2( 1, 1) * scale, 0, 1);

    clip_coord1[0] = M[1] * vec4(vec2(-1,-1) * scale, 0, 1);
    clip_coord1[1] = M[1] * vec4(vec2( 1,-1) * scale, 0, 1);
    clip_coord1[2] = M[1] * vec4(vec2(-1, 1) * scale, 0, 1);
    clip_coord1[3] = M[1] * vec4(vec2( 1, 1) * scale, 0, 1);
 
    // BOTTOM
    {
        vec4 n0 = N[0] * vec4( 0, -1, 0, 0);
        vec4 n1 = N[1] * vec4( 0, -1, 0, 0);
        emit_vertex(clip_coord0[0], n0, 0);
        emit_vertex(clip_coord0[1], n0, 0);
        emit_vertex(clip_coord1[0], n1, 1);
        emit_vertex(clip_coord1[1], n1, 1);
        EndPrimitive();
    }
    // TOP
    {
        vec4 n0 = N[0] * vec4( 0, 1, 0, 0);
        vec4 n1 = N[1] * vec4( 0, 1, 0, 0);
        emit_vertex(clip_coord0[3], n0, 0);
        emit_vertex(clip_coord0[2], n0, 0);
        emit_vertex(clip_coord1[3], n1, 1);
        emit_vertex(clip_coord1[2], n1, 1);
        EndPrimitive();
    }
    // LEFT
    {
        vec4 n0 = N[0] * vec4( -1, 0, 0, 0);
        vec4 n1 = N[1] * vec4( -1, 0, 0, 0);
        emit_vertex(clip_coord0[2], n0, 0);
        emit_vertex(clip_coord0[0], n0, 0);
        emit_vertex(clip_coord1[2], n1, 1);
        emit_vertex(clip_coord1[0], n1, 1);
        EndPrimitive();
    }
    // RIGHT
    {
        vec4 n0 = N[0] * vec4( 1, 0, 0, 0);
        vec4 n1 = N[1] * vec4( 1, 0, 0, 0);
        emit_vertex(clip_coord0[1], n0, 0);
        emit_vertex(clip_coord0[3], n0, 0);
        emit_vertex(clip_coord1[1], n1, 1);
        emit_vertex(clip_coord1[3], n1, 1);
        EndPrimitive();
    }
    // FRONT
    if ((f[0] & 1U) != 0U) {
        vec4 n = N[0] * vec4(0, 0, -1, 0);
        emit_vertex(clip_coord0[1], n, 0);
        emit_vertex(clip_coord0[0], n, 0);
        emit_vertex(clip_coord0[3], n, 0);
        emit_vertex(clip_coord0[2], n, 0);
        EndPrimitive();
    }
    // BACK
    if ((f[1] & 2U) != 0U) {
        vec4 n = N[1] * vec4(0, 0, 1, 0);
        emit_vertex(clip_coord1[0], n, 1);
        emit_vertex(clip_coord1[1], n, 1);
        emit_vertex(clip_coord1[2], n, 1);
        emit_vertex(clip_coord1[3], n, 1);
        EndPrimitive();
    }
}
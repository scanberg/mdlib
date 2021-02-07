#version 330 core

#ifndef ATOM_COL
#define ATOM_COL 1
#endif

#ifndef ATOM_VEL
#define ATOM_VEL 0
#endif

#ifndef ATOM_IDX
#define ATOM_IDX 1
#endif

#ifndef VIEW_NORM
#define VIEW_NORM 1
#endif

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
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
#if ATOM_COL
    vec4 color;
#endif
#if ATOM_VEL
    vec4 velocity;
#endif
#if ATOM_IDX
    uint picking_idx;
#endif
    uint flags;
} in_vert[];

#if ATOM_COL || ATOM_VEL || ATOM_IDX || VIEW_NORM
out Fragment {
#if ATOM_VEL
    smooth vec4 curr_clip_coord;
    smooth vec4 prev_clip_coord;
#endif
#if ATOM_COL
    smooth vec4 color;
#endif
#if VIEW_NORM
    smooth vec3 view_normal;
#endif
#if ATOM_IDX
    flat   uint picking_idx;
#endif
} out_frag;
#endif

void emit_vertex(in vec4 curr_clip_coord, in vec4 prev_clip_coord, in vec4 normal, in int idx) {
#if ATOM_COL
    out_frag.color = in_vert[idx].color;
#endif
#if ATOM_VEL
    out_frag.curr_clip_coord = curr_clip_coord;
    out_frag.prev_clip_coord = prev_clip_coord;
#endif
#if VIEW_NORM
    out_frag.view_normal = normal.xyz;
#endif
#if ATOM_IDX
    out_frag.picking_idx = in_vert[idx].picking_idx;
#endif
    gl_Position = curr_clip_coord;
    EmitVertex();
}


void main() {
    if ((in_vert[0].flags & u_atom_mask) != u_atom_mask ||
        (in_vert[1].flags & u_atom_mask) != u_atom_mask) {
        return;
    }

    vec4 p[2];
    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    
    p[0] = in_vert[0].control_point;
    p[1] = in_vert[1].control_point;
    x[0] = in_vert[0].support_vector;
    x[1] = in_vert[1].support_vector * sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    z[0] = in_vert[0].support_tangent;
    z[1] = in_vert[1].support_tangent;
    y[0] = vec4(normalize(cross(z[0].xyz, x[0].xyz)), 0); // To maintain right-handedness
    y[1] = vec4(normalize(cross(z[1].xyz, x[1].xyz)), 0);

#if 0
    // It is possible to fix tesselation direction so it follows the twist of the spline
    // By tesselating with respect to the twist we would avoid some thing long triangles.
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
    vec4 v0[4];
    vec4 v1[4];
    vec4 curr_clip0[4];
    vec4 curr_clip1[4];
    vec4 prev_clip0[4];
    vec4 prev_clip1[4];

    M[0] = mat4(x[0], y[0], z[0], p[0]);
    M[1] = mat4(x[1], y[1], z[1], p[1]);

#if VIEW_NORM
    N[0] = u_world_to_view_normal * M[0];
    N[1] = u_world_to_view_normal * M[1];
#endif

    v0[0] = M[0] * vec4(vec2(-1,-1) * scale, 0, 1);
    v0[1] = M[0] * vec4(vec2( 1,-1) * scale, 0, 1);
    v0[2] = M[0] * vec4(vec2(-1, 1) * scale, 0, 1);
    v0[3] = M[0] * vec4(vec2( 1, 1) * scale, 0, 1);

    curr_clip0[0] = u_world_to_clip * v0[0];
    curr_clip0[1] = u_world_to_clip * v0[1];
    curr_clip0[2] = u_world_to_clip * v0[2];
    curr_clip0[3] = u_world_to_clip * v0[3];

#if ATOM_VEL
    prev_clip0[0] = u_prev_world_to_clip * (v0[0] - in_vert[0].velocity);
    prev_clip0[1] = u_prev_world_to_clip * (v0[1] - in_vert[0].velocity);
    prev_clip0[2] = u_prev_world_to_clip * (v0[2] - in_vert[0].velocity);
    prev_clip0[3] = u_prev_world_to_clip * (v0[3] - in_vert[0].velocity);
#endif

    v1[0] = M[1] * vec4(vec2(-1,-1) * scale, 0, 1);
    v1[1] = M[1] * vec4(vec2( 1,-1) * scale, 0, 1);
    v1[2] = M[1] * vec4(vec2(-1, 1) * scale, 0, 1);
    v1[3] = M[1] * vec4(vec2( 1, 1) * scale, 0, 1);

    curr_clip1[0] = u_world_to_clip * v1[0];
    curr_clip1[1] = u_world_to_clip * v1[1];
    curr_clip1[2] = u_world_to_clip * v1[2];
    curr_clip1[3] = u_world_to_clip * v1[3];

#if ATOM_VEL
    prev_clip1[0] = u_prev_world_to_clip * (v1[0] - in_vert[1].velocity);
    prev_clip1[1] = u_prev_world_to_clip * (v1[1] - in_vert[1].velocity);
    prev_clip1[2] = u_prev_world_to_clip * (v1[2] - in_vert[1].velocity);
    prev_clip1[3] = u_prev_world_to_clip * (v1[3] - in_vert[1].velocity);
#endif
 
    // BOTTOM
    {
        vec4 n0 = N[0] * vec4( 0, -1, 0, 0);
        vec4 n1 = N[1] * vec4( 0, -1, 0, 0);
        emit_vertex(curr_clip0[0], prev_clip0[0], n0, 0);
        emit_vertex(curr_clip0[1], prev_clip0[1], n0, 0);
        emit_vertex(curr_clip1[0], prev_clip1[0], n1, 1);
        emit_vertex(curr_clip1[1], prev_clip1[1], n1, 1);
        EndPrimitive();
    }
    // TOP
    {
        vec4 n0 = N[0] * vec4( 0, 1, 0, 0);
        vec4 n1 = N[1] * vec4( 0, 1, 0, 0);
        emit_vertex(curr_clip0[3], prev_clip0[3], n0, 0);
        emit_vertex(curr_clip0[2], prev_clip0[2], n0, 0);
        emit_vertex(curr_clip1[3], prev_clip1[3], n1, 1);
        emit_vertex(curr_clip1[2], prev_clip1[2], n1, 1);
        EndPrimitive();
    }
    // LEFT
    {
        vec4 n0 = N[0] * vec4( -1, 0, 0, 0);
        vec4 n1 = N[1] * vec4( -1, 0, 0, 0);
        emit_vertex(curr_clip0[2], prev_clip0[2], n0, 0);
        emit_vertex(curr_clip0[0], prev_clip0[0], n0, 0);
        emit_vertex(curr_clip1[2], prev_clip1[2], n1, 1);
        emit_vertex(curr_clip1[0], prev_clip1[0], n1, 1);
        EndPrimitive();
    }
    // RIGHT
    {
        vec4 n0 = N[0] * vec4( 1, 0, 0, 0);
        vec4 n1 = N[1] * vec4( 1, 0, 0, 0);
        emit_vertex(curr_clip0[1], prev_clip0[1], n0, 0);
        emit_vertex(curr_clip0[3], prev_clip0[3], n0, 0);
        emit_vertex(curr_clip1[1], prev_clip1[1], n1, 1);
        emit_vertex(curr_clip1[3], prev_clip1[3], n1, 1);
        EndPrimitive();
    }
    // FRONT
    {
        vec4 n = N[0] * vec4(0, 0, -1, 0);
        emit_vertex(curr_clip0[1], prev_clip0[1], n, 0);
        emit_vertex(curr_clip0[0], prev_clip0[0], n, 0);
        emit_vertex(curr_clip1[3], prev_clip1[3], n, 0);
        emit_vertex(curr_clip1[2], prev_clip1[2], n, 0);
        EndPrimitive();
    }
    // BACK
    {
        vec4 n = N[1] * vec4(0, 0, 1, 0);
        emit_vertex(curr_clip0[0], prev_clip0[0], n, 1);
        emit_vertex(curr_clip0[1], prev_clip0[1], n, 1);
        emit_vertex(curr_clip1[2], prev_clip1[2], n, 1);
        emit_vertex(curr_clip1[3], prev_clip1[3], n, 1);
        EndPrimitive();
    }
}
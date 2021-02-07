#version 330 core

#define RES 6
//#define NUM_VERTS 26    // RES * 2 + 2 (Cylinder) + 1 + RES (CAP)
#define NUM_VERTS 26

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
layout(triangle_strip, max_vertices = NUM_VERTS) out;

in Vertex {
    vec4  control_point;
    vec4  support_vector;
    vec4  support_tangent;
#if ATOM_VEL
    vec4  velocity;
#endif
    float segment_t;
    vec3  secondary_structure;
    uint  flags;
#if ATOM_COL
    vec4  color;
#endif
#if ATOM_IDX
    flat uint  picking_idx;
#endif
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
    smooth vec4 view_normal;
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
    out_frag.view_normal = normal;
#endif
#if ATOM_IDX
    out_frag.picking_idx = in_vert[0].picking_idx;
#endif
    gl_Position = curr_clip_coord;
    EmitVertex();
}

void main() {
    // We consider the rendered segment to fully belong to index 0
    if ((in_vert[0].flags & u_atom_mask) != u_atom_mask) return;
#if ATOM_COL
    if (in_vert[0].color.a == 0.0f) return;
#endif

    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    vec4 w[2];
    mat4 M[2];
    mat4 N[2];
    mat4 C[2];
    vec3 ss[2];
    vec2 s[2];
    vec4 v0[RES];
    vec4 v1[RES];
    vec4 p0[RES];
    vec4 p1[RES];   
    vec4 n0[RES];
    vec4 n1[RES];
    
    x[0] = in_vert[0].support_vector;
    x[1] = in_vert[1].support_vector * sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    z[0] = in_vert[0].support_tangent;
    z[1] = in_vert[1].support_tangent;
    y[0] = vec4(normalize(cross(z[0].xyz, x[0].xyz)), 0); // To maintain right-handedness
    y[1] = vec4(normalize(cross(z[1].xyz, x[1].xyz)), 0);
    w[0] = in_vert[0].control_point;
    w[1] = in_vert[1].control_point;

    const float TWO_PI = 2.0 * 3.14159265;

    M[0] = mat4(x[0], y[0], z[0], w[0]);
    M[1] = mat4(x[1], y[1], z[1], w[1]);

#if VIEW_NORM
    N[0] = u_world_to_view_normal * M[0];
    N[1] = u_world_to_view_normal * M[1];
#endif

    C[0] = u_world_to_clip * M[0];
    C[1] = u_world_to_clip * M[1];

    ss[0] = in_vert[0].secondary_structure;
    ss[1] = in_vert[1].secondary_structure;

    // Elipse profile scaling based on secondary structure
    const vec2 coil_scale  = vec2(0.2, 0.2);
    const vec2 helix_scale = vec2(1.2, 0.2);
    const vec2 sheet_scale = vec2(1.5, 0.1);

    vec2 scale = vec2(u_width_scale, u_height_scale);

    s[0] = scale * (ss[0].x * coil_scale + ss[0].y * helix_scale + ss[0].z * sheet_scale);
    s[1] = scale * (ss[1].x * coil_scale + ss[1].y * helix_scale + ss[1].z * sheet_scale);



/*
    if ((in_vert[0].flags & 8U) != 0U || (in_vert[1].flags & 8U) != 0U) {
        s[0] = vec2(1.5, 0.4) * (1.0 - fract(in_vert[0].segment_t));
        s[1] = vec2(1.5, 0.4) * (1.0 - fract(in_vert[1].segment_t));
    }
*/
    for (int i = 0; i < RES; ++i) {
        float t = float(i) / float(RES) * TWO_PI;
        vec2 x = vec2(cos(t), sin(t));
        v0[i] = C[0] * vec4(x * s[0],0,1);
        v1[i] = C[1] * vec4(x * s[1],0,1);
#if ATOM_VEL
        p0[i] = u_prev_world_to_clip * (M[0] * vec4(x * s[0], 0, 1) - in_vert[0].velocity);
        p1[i] = u_prev_world_to_clip * (M[1] * vec4(x * s[1], 0, 1) - in_vert[1].velocity);
#endif
#if VIEW_NORM
        n0[i] = N[0] * vec4(x / s[0],0,0);
        n1[i] = N[1] * vec4(x / s[1],0,0);
#endif
    }
    
    for (int i = 0; i < RES; ++i) {
        emit_vertex(v1[i], p1[i], n1[i], 1);
        emit_vertex(v0[i], p0[i], n0[i], 0);
    }
    emit_vertex(v1[0], p1[0], n1[0], 1);
    emit_vertex(v0[0], p0[0], n0[0], 0);
    EndPrimitive();

    // MAKE CAP
    if ((in_vert[0].flags & 1U) != 0U) { // BEG_CHAIN
        vec4 cc = u_world_to_clip * M[0][3];
        vec4 cp = u_prev_world_to_clip * M[0][3];
        vec4 cn = N[0] * vec4(0,0,-1,0);
        for (int i = 0; i < RES-1; i += 2) {
            emit_vertex(v0[i], p0[i], cn, 0);
            emit_vertex(cc, cp, cn, 0);
            emit_vertex(v0[i+1], p0[i+1], cn, 0);
            emit_vertex(v0[(i+2)%RES], p0[(i+2)%RES], cn, 0);
        }
        EndPrimitive();
    } else if ((in_vert[1].flags & 2U) != 0U) { // END_CHAIN
        vec4 cc = u_world_to_clip * M[1][3];
        vec4 cp = u_prev_world_to_clip * M[1][3];
        vec4 cn = N[1] * vec4(0,0,1,0);
        for (int i = 0; i < RES-1; i += 2) {
            emit_vertex(v1[i], p1[i], cn, 0);
            emit_vertex(v1[i+1], p1[i+1], cn, 0);
            emit_vertex(cc, cp, cn, 0);
            emit_vertex(v1[(i+2)%RES], p1[(i+2)%RES], cn, 0);
        }
        EndPrimitive();
    }
}
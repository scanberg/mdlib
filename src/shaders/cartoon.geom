#version 330 core

#define RES 6
//#define NUM_VERTS 26    // RES * 2 + 2 (Cylinder) + 1 + RES (CAP)
#define NUM_VERTS 26
#define TWO_PI 6.28318530717958647693

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
    vec4 u_scale;
};

layout(lines) in;
layout(triangle_strip, max_vertices = NUM_VERTS) out;

in Vertex {
    vec3  control_point;
    vec3  support_vector;
    vec3  support_tangent;
    vec3  view_velocity;
    vec4  color;
    vec3  secondary_structure;
    float segment_t;
    uint  spline_flags;
    uint  atom_flags;
    uint  picking_idx;
} in_vert[];

out Fragment {
    smooth vec3 view_coord;
    smooth vec3 view_velocity;
    smooth vec4 color;
    smooth vec3 view_normal;
    flat   uint picking_idx;
} out_frag;

void emit_vertex(in vec4 clip_coord, in vec4 normal, in int idx) {
    out_frag.color = in_vert[idx].color;
    //out_frag.color = mix(vec4(1,1,1,1), vec4(1,0,0,1), fract(in_vert[0].segment_t));
    vec4 view_coord = u_clip_to_view * clip_coord;
    out_frag.view_coord     = view_coord.xyz / view_coord.w;
    out_frag.view_velocity  = in_vert[idx].view_velocity;
    out_frag.view_normal = normal.xyz;
    out_frag.picking_idx = in_vert[0].picking_idx;
    gl_Position = clip_coord;
    EmitVertex();
}

void main() {
    // We consider the rendered segment to fully belong to index 0
    if ((in_vert[0].atom_flags & u_atom_mask) != u_atom_mask) return;
    if (in_vert[0].color.a == 0.0f) return;

    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    vec4 w[2];
    mat4 M[2];
    mat4 N[2];
    vec3 ss[2];
    vec2 s[2];
    vec4 v0[RES];
    vec4 v1[RES];
    vec4 n0[RES];
    vec4 n1[RES];
    
    float flip = sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    x[0] = vec4(in_vert[0].support_vector, 0);
    x[1] = vec4(in_vert[1].support_vector * flip, 0);
    z[0] = vec4(in_vert[0].support_tangent, 0);
    z[1] = vec4(in_vert[1].support_tangent, 0);
    y[0] = vec4(normalize(cross(z[0].xyz, x[0].xyz)), 0); // To maintain right-handedness
    y[1] = vec4(normalize(cross(z[1].xyz, x[1].xyz)), 0);
    w[0] = vec4(in_vert[0].control_point, 1);
    w[1] = vec4(in_vert[1].control_point, 1);

    M[0] = u_world_to_clip * mat4(x[0], y[0], z[0], w[0]);
    M[1] = u_world_to_clip * mat4(x[1], y[1], z[1], w[1]);

    N[0] = u_world_to_view_normal * mat4(x[0], y[0], z[0], w[0]);
    N[1] = u_world_to_view_normal * mat4(x[1], y[1], z[1], w[1]);

    ss[0] = in_vert[0].secondary_structure;
    ss[1] = in_vert[1].secondary_structure;

    // Elipse profile scaling based on secondary structure
    vec2 coil_scale  = vec2(0.2, 0.2) * u_scale[0];
    vec2 helix_scale = vec2(1.2, 0.2) * u_scale[1];
    vec2 sheet_scale = vec2(1.5, 0.1) * u_scale[2];

    s[0] = ss[0].x * coil_scale + ss[0].y * helix_scale + ss[0].z * sheet_scale;
    s[1] = ss[1].x * coil_scale + ss[1].y * helix_scale + ss[1].z * sheet_scale;

    for (int i = 0; i < RES; ++i) {
        float t = float(i) / float(RES) * TWO_PI;
        vec2 x = vec2(cos(t), sin(t));
        v0[i] = M[0] * vec4(x * s[0],0,1);
        v1[i] = M[1] * vec4(x * s[1],0,1);

        n0[i] = N[0] * vec4(x / s[0],0,0);
        n1[i] = N[1] * vec4(x / s[1],0,0);
    }
    
    for (int i = 0; i < RES; ++i) {
        emit_vertex(v1[i], n1[i], 1);
        emit_vertex(v0[i], n0[i], 0);
    }
    emit_vertex(v1[0], n1[0], 1);
    emit_vertex(v0[0], n0[0], 0);
    EndPrimitive();

    // MAKE CAP
    if ((in_vert[0].spline_flags & 1U) != 0U) { // BEG_CHAIN
        vec4 cc = M[0][3];
        vec4 cn = N[0] * vec4(0,0,-1,0);
        for (int i = 0; i < RES-1; i += 2) {
            emit_vertex(v0[i], cn, 0);
            emit_vertex(cc, cn, 0);
            emit_vertex(v0[i+1], cn, 0);
            emit_vertex(v0[(i+2)%RES], cn, 0);
        }
        EndPrimitive();
    } 
    if ((in_vert[1].spline_flags & 2U) != 0U) { // END_CHAIN
        vec4 cc = M[1][3];
        vec4 cn = N[1] * vec4(0,0,1,0);
        for (int i = 0; i < RES-1; i += 2) {
            emit_vertex(v1[i], cn, 0);
            emit_vertex(v1[i+1], cn, 0);
            emit_vertex(cc, cn, 0);
            emit_vertex(v1[(i+2)%RES], cn, 0);
        }
        EndPrimitive();
    }
}
#version 330 core
#extension GL_ARB_shading_language_packing : enable

uniform int u_enable_smoothing = 0;
uniform float u_smooth_weight_coil = 0.20;
uniform float u_smooth_weight_helix = 0.10;
uniform float u_smooth_weight_sheet = 0.15;

uniform usamplerBuffer u_buf_control_point_words;                      // Interleaved gl_control_point_t as uint32 words
uniform usamplerBuffer u_buf_smooth_neighbors;                         // uvec4(prev1, next1, prev2, next2)

layout (location = 0) in vec3  in_position;
layout (location = 1) in uint  in_atom_idx;
layout (location = 2) in vec3  in_velocity;
layout (location = 3) in float in_segment_t;
layout (location = 4) in uint  in_secondary_structure_and_flags;
layout (location = 5) in uvec3 in_support_and_tangent_vector;

out vec3  out_position;
out uint  out_atom_idx;
out vec3  out_velocity;
out float out_segment_t;
out uint  out_secondary_structure_and_flags;
out uvec3 out_support_and_tangent_vector;

vec3 safe_normalize(vec3 v, vec3 fallback) {
    float d2 = dot(v, v);
    return d2 > 0.0 ? v * inversesqrt(d2) : fallback;
}

vec3 unpack_support(uint w0, uint w1, uint w2) {
    vec2 xy = unpackSnorm2x16(w0);
    vec2 zx = unpackSnorm2x16(w1);
    return vec3(xy.x, xy.y, zx.x);
}

vec3 unpack_tangent(uint w0, uint w1, uint w2) {
    vec2 zx = unpackSnorm2x16(w1);
    vec2 yz = unpackSnorm2x16(w2);
    return vec3(zx.y, yz.x, yz.y);
}

vec3 load_cp_position(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return vec3(
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 0).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 1).r),
        uintBitsToFloat(texelFetch(u_buf_control_point_words, base + 2).r)
    );
}

uvec3 load_cp_support_tangent(uint cp_idx) {
    int base = int(cp_idx * 12u);
    return uvec3(
        texelFetch(u_buf_control_point_words, base + 9).r,
        texelFetch(u_buf_control_point_words, base + 10).r,
        texelFetch(u_buf_control_point_words, base + 11).r
    );
}

vec3 project_point_on_axis(vec3 p, vec3 c, vec3 a) {
    return c + a * dot(p - c, a);
}

vec3 orient_like(vec3 v, vec3 ref) {
    return dot(v, ref) < 0.0 ? -v : v;
}

void main() {
    out_position = in_position;
    out_atom_idx = in_atom_idx;
    out_velocity = in_velocity;
    out_segment_t = in_segment_t;
    out_secondary_structure_and_flags = in_secondary_structure_and_flags;
    out_support_and_tangent_vector = in_support_and_tangent_vector;

    if (u_enable_smoothing == 0) {
        return;
    }

    uint cp_idx = uint(gl_VertexID);
    uvec4 nb = texelFetch(u_buf_smooth_neighbors, int(cp_idx));
    uint prev_idx = nb.x;
    uint next_idx = nb.y;
    uint prev2_idx = nb.z;
    uint next2_idx = nb.w;

    if (prev_idx == cp_idx || next_idx == cp_idx) {
        return;
    }

    uvec4 nb_prev2 = texelFetch(u_buf_smooth_neighbors, int(prev2_idx));
    uvec4 nb_next2 = texelFetch(u_buf_smooth_neighbors, int(next2_idx));
    uint prev3_idx = nb_prev2.x;
    uint next3_idx = nb_next2.y;
    uvec4 nb_prev3 = texelFetch(u_buf_smooth_neighbors, int(prev3_idx));
    uvec4 nb_next3 = texelFetch(u_buf_smooth_neighbors, int(next3_idx));
    uint prev4_idx = nb_prev3.x;
    uint next4_idx = nb_next3.y;

    vec3 ss = unpackUnorm4x8(in_secondary_structure_and_flags).xyz;
    float coil_w  = ss.x;
    float helix_w = ss.y;
    float sheet_w = ss.z;

    vec3 prev_pos = load_cp_position(prev_idx);
    vec3 next_pos = load_cp_position(next_idx);
    vec3 prev2_pos = load_cp_position(prev2_idx);
    vec3 next2_pos = load_cp_position(next2_idx);
    vec3 prev3_pos = load_cp_position(prev3_idx);
    vec3 next3_pos = load_cp_position(next3_idx);
    vec3 prev4_pos = load_cp_position(prev4_idx);
    vec3 next4_pos = load_cp_position(next4_idx);

    vec3 cur_support = unpack_support(in_support_and_tangent_vector.x, in_support_and_tangent_vector.y, in_support_and_tangent_vector.z);
    vec3 cur_tangent = unpack_tangent(in_support_and_tangent_vector.x, in_support_and_tangent_vector.y, in_support_and_tangent_vector.z);
    cur_tangent = safe_normalize(cur_tangent, safe_normalize(next_pos - prev_pos, vec3(1.0, 0.0, 0.0)));

    uvec3 prev_svt  = load_cp_support_tangent(prev_idx);
    uvec3 next_svt  = load_cp_support_tangent(next_idx);
    uvec3 prev2_svt = load_cp_support_tangent(prev2_idx);
    uvec3 next2_svt = load_cp_support_tangent(next2_idx);
    uvec3 prev3_svt = load_cp_support_tangent(prev3_idx);
    uvec3 next3_svt = load_cp_support_tangent(next3_idx);
    uvec3 prev4_svt = load_cp_support_tangent(prev4_idx);
    uvec3 next4_svt = load_cp_support_tangent(next4_idx);

    vec3 prev_support  = unpack_support(prev_svt.x, prev_svt.y, prev_svt.z);
    vec3 next_support  = unpack_support(next_svt.x, next_svt.y, next_svt.z);
    vec3 prev2_support = unpack_support(prev2_svt.x, prev2_svt.y, prev2_svt.z);
    vec3 next2_support = unpack_support(next2_svt.x, next2_svt.y, next2_svt.z);
    vec3 prev3_support = unpack_support(prev3_svt.x, prev3_svt.y, prev3_svt.z);
    vec3 next3_support = unpack_support(next3_svt.x, next3_svt.y, next3_svt.z);
    vec3 prev4_support = unpack_support(prev4_svt.x, prev4_svt.y, prev4_svt.z);
    vec3 next4_support = unpack_support(next4_svt.x, next4_svt.y, next4_svt.z);

    prev_support  *= sign(dot(prev_support,  cur_support));
    next_support  *= sign(dot(next_support,  cur_support));
    prev2_support *= sign(dot(prev2_support, cur_support));
    next2_support *= sign(dot(next2_support, cur_support));
    prev3_support *= sign(dot(prev3_support, cur_support));
    next3_support *= sign(dot(next3_support, cur_support));
    prev4_support *= sign(dot(prev4_support, cur_support));
    next4_support *= sign(dot(next4_support, cur_support));

    float coil_alpha  = clamp(u_smooth_weight_coil,  0.0, 1.0);
    float helix_alpha = clamp(u_smooth_weight_helix, 0.0, 1.0);
    float sheet_alpha = clamp(u_smooth_weight_sheet, 0.0, 1.0);

    vec3 avg_pos = (prev_pos + in_position + next_pos) * (1.0 / 3.0);
    vec3 avg_pos_long = (prev2_pos + prev_pos + in_position + next_pos + next2_pos) * 0.2;
    vec3 avg_pos_sheet = (prev4_pos + prev3_pos + prev2_pos + prev_pos + in_position + next_pos + next2_pos + next3_pos + next4_pos) * (1.0 / 9.0);

    vec3 pos_coil = in_position;

    vec3 helix_axis = safe_normalize(
        (next3_pos - prev3_pos) * 0.50 +
        (next2_pos - prev2_pos) * 0.75 +
        (next_pos  - prev_pos)  * 1.00,
        safe_normalize(next_pos - prev_pos, cur_tangent)
    );
    helix_axis = orient_like(helix_axis, cur_tangent);
    vec3 helix_center = (prev2_pos + prev_pos + next_pos + next2_pos) * 0.25;
    vec3 proj_cur = project_point_on_axis(in_position, helix_center, helix_axis);
    vec3 proj_tgt = project_point_on_axis(avg_pos_long, helix_center, helix_axis);
    vec3 radial_cur = in_position - proj_cur;
    vec3 pos_helix_target = proj_tgt + radial_cur;
    vec3 pos_helix = mix(in_position, pos_helix_target, helix_alpha);

    vec3 pos_sheet = mix(in_position, avg_pos_sheet, sheet_alpha);

    out_position = pos_coil * coil_w + pos_helix * helix_w + pos_sheet * sheet_w;

    vec3 tan_pos_short = safe_normalize(next_pos - prev_pos, cur_tangent);
    vec3 tan_pos_long  = safe_normalize(next2_pos - prev2_pos, tan_pos_short);
    vec3 tan_pos_longer = safe_normalize(next3_pos - prev3_pos, tan_pos_long);
    vec3 tan_pos_sheet = safe_normalize(
        (next4_pos - prev4_pos) * 0.55 +
        (next3_pos - prev3_pos) * 0.80 +
        (next2_pos - prev2_pos) * 1.00,
        tan_pos_longer
    );
    tan_pos_short = orient_like(tan_pos_short, cur_tangent);
    tan_pos_long  = orient_like(tan_pos_long, tan_pos_short);
    tan_pos_longer = orient_like(tan_pos_longer, tan_pos_long);
    tan_pos_sheet = orient_like(tan_pos_sheet, tan_pos_longer);
    vec3 tan_coil  = safe_normalize(mix(cur_tangent, tan_pos_short, coil_alpha  * 0.40), tan_pos_short);
    vec3 tan_helix = safe_normalize(mix(cur_tangent, helix_axis,      helix_alpha * 0.90), helix_axis);
    vec3 tan_sheet = safe_normalize(mix(cur_tangent, tan_pos_sheet,   sheet_alpha * 0.88), tan_pos_sheet);
    vec3 smooth_tangent = safe_normalize(tan_coil * coil_w + tan_helix * helix_w + tan_sheet * sheet_w, cur_tangent);
    smooth_tangent = orient_like(smooth_tangent, cur_tangent);

    vec3 support_avg_short = safe_normalize(prev_support + cur_support + next_support, cur_support);
    vec3 support_avg_long  = safe_normalize(prev2_support + prev_support + cur_support + next_support + next2_support, support_avg_short);
    vec3 support_avg_longer = safe_normalize(prev3_support + prev2_support + prev_support + cur_support + next_support + next2_support + next3_support, support_avg_long);
    vec3 support_avg_sheet = safe_normalize(prev4_support + prev3_support + prev2_support + prev_support + cur_support + next_support + next2_support + next3_support + next4_support, support_avg_longer);
    vec3 support_coil = safe_normalize(mix(cur_support, support_avg_short, coil_alpha * 0.35), cur_support);

    vec3 transport_prev2 = safe_normalize(prev2_support - smooth_tangent * dot(prev2_support, smooth_tangent), support_avg_long);
    vec3 transport_prev  = safe_normalize(prev_support  - smooth_tangent * dot(prev_support,  smooth_tangent), support_avg_long);
    vec3 transport_next  = safe_normalize(next_support  - smooth_tangent * dot(next_support,  smooth_tangent), support_avg_long);
    vec3 transport_next2 = safe_normalize(next2_support - smooth_tangent * dot(next2_support, smooth_tangent), support_avg_long);
    vec3 helix_transport = safe_normalize(transport_prev2 + transport_prev + cur_support + transport_next + transport_next2, support_avg_long);
    vec3 helix_ref = safe_normalize(cross(smooth_tangent, cross(helix_axis, smooth_tangent)), helix_transport);
    helix_ref = orient_like(helix_ref, cur_support);
    vec3 support_helix = safe_normalize(mix(helix_transport, helix_ref, helix_alpha * 0.35), helix_transport);
    support_helix = orient_like(support_helix, cur_support);

    vec3 c0 = cross(next2_pos - prev2_pos, next4_pos - prev4_pos);
    vec3 c1 = cross(next3_pos - in_position, next_pos - in_position);
    vec3 sheet_plane_n = safe_normalize(c0 + c1, support_avg_sheet);
    sheet_plane_n = orient_like(sheet_plane_n, cur_support);
    vec3 sheet_support = safe_normalize(sheet_plane_n - smooth_tangent * dot(sheet_plane_n, smooth_tangent), support_avg_sheet);
    sheet_support = orient_like(sheet_support, cur_support);
    vec3 support_sheet = safe_normalize(mix(support_avg_sheet, sheet_support, sheet_alpha * 0.93), sheet_support);
    support_sheet = orient_like(support_sheet, cur_support);

    vec3 smooth_support = safe_normalize(support_coil * coil_w + support_helix * helix_w + support_sheet * sheet_w, cur_support);
    smooth_support = safe_normalize(smooth_support - smooth_tangent * dot(smooth_support, smooth_tangent), smooth_support);
    smooth_support = orient_like(smooth_support, cur_support);

    out_support_and_tangent_vector[0] = packSnorm2x16(smooth_support.xy);
    out_support_and_tangent_vector[1] = packSnorm2x16(vec2(smooth_support.z, smooth_tangent.x));
    out_support_and_tangent_vector[2] = packSnorm2x16(smooth_tangent.yz);
}

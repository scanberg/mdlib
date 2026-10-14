#version 330 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_shading_language_packing : enable

// Aromatic rings raycast within the quads of ring.vert: a torus in the plane of the ring.
//
// The ray is intersected analytically (the quartic of Inigo Quilez, torIntersect), which is only accurate in float when
// its coefficients are of order one: they grow with the fourth power of the distance of the origin from the torus, and
// a hit is decided by terms of the order R^2 r^2. So the ray is taken to the frame of the torus scaled by 1 / R (the
// torus is then (1, r / R) whatever the size of the ring, which also gives the fixed threshold of torIntersect the same
// meaning for every ring), and its origin moved to where it enters a sphere about it.
//
// Measured against a double precision reference (rays at rings 2 to 2000 Å away, a quarter of the views near edge on,
// r / R from 0.03 to 0.2), three choices matter in float:
//   - The sphere has the radius 1 + 2 r / R, so the origin stays a tube radius clear of the torus. The bounding sphere
//     itself touches the torus along its outer equator, an origin there makes the constant term of the quartic vanish,
//     and the inversion torIntersect then takes (1 / k0) overflows: about 1% of the rays that hit missed.
//   - torIntersect inverts the quartic below 0.001 rather than 0.01: the inverted quartic has large coefficients and
//     its discriminant cancels in float, so it is best kept to the cases that need it. Never inverting is far worse.
//   - Two steps along the ray by the distance to the torus polish the root, a tenth of the error where it was largest.
// What remains is about one ray in 10^4 at r / R ~ 0.1 (one in 10^3 for the thinnest tubes), almost all of them
// grazing the outline of a ring seen nearly edge on.
//
// The circle is colored by the atoms of the ring as a bond is by its two: the hit lies between two consecutive atoms in
// angle about the normal, as along a bond between its ends (u_mode). It picks as the ring.

#ifndef ORTHO
#define ORTHO 0
#endif

#define MODE_NEAREST 0
#define MODE_SMOOTH  1
#define MODE_UNIFORM 2

#define PI 3.14159265358979

#ifndef GL_ARB_shading_language_packing
vec4 unpackUnorm4x8(uint p) {
    uvec4 iv = uvec4(p & 0xFFU, (p >> 8U) & 0xFFU, (p >> 16U) & 0xFFU, (p >> 24U) & 0xFFU);
    return vec4(iv) * (1.0f / 255.0f);
}
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
    uint u_atom_base_index;
    uint u_bond_base_index;
    uint u_backbone_base_index;

    float u_radius;
    float u_max_d2;
    int   u_mode;
    float u_sharpness;
    uint  u_uniform_color;
    float u_spread;
    uint  u_ring_base_index;
};

in Fragment {
    flat vec3  center;
    flat vec3  normal;
    flat vec3  tangent;
    flat vec2  radii;
    flat vec4  angle[2];
    flat uvec4 color[2];
    flat uvec2 count_picking;
    flat vec3  view_vel;
    smooth vec3 view_pos;
} in_frag;

// The quad of ring.vert is in front of the ring, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

// Ray (ro, rd) against the torus about the z axis of major radius tor.x and tube radius tor.y. rd is unit length.
// Returns t of the nearest hit ahead of ro, or -1 or 1e20 for none. From Inigo Quilez
// (https://iquilezles.org/articles/intersectors), with the argument of acos clamped (rounding takes it past +-1 near
// tangency, and the NaN would turn every root into a miss) and the inversion below 0.001 rather than 0.01 (see above).
float torIntersect(in vec3 ro, in vec3 rd, in vec2 tor) {
    float po = 1.0;
    float Ra2 = tor.x * tor.x;
    float ra2 = tor.y * tor.y;
    float m = dot(ro, ro);
    float n = dot(ro, rd);
    float k = (m + Ra2 - ra2) / 2.0;
    float k3 = n;
    float k2 = n * n - Ra2 * dot(rd.xy, rd.xy) + k;
    float k1 = n * k - Ra2 * dot(rd.xy, ro.xy);
    float k0 = k * k - Ra2 * dot(ro.xy, ro.xy);

    if (abs(k3 * (k3 * k3 - k2) + k1) < 0.001) {
        po = -1.0;
        float tmp = k1; k1 = k3; k3 = tmp;
        k0 = 1.0 / k0;
        k1 = k1 * k0;
        k2 = k2 * k0;
        k3 = k3 * k0;
    }

    float c2 = k2 * 2.0 - 3.0 * k3 * k3;
    float c1 = k3 * (k3 * k3 - k2) + k1;
    float c0 = k3 * (k3 * (c2 + 2.0 * k2) - 8.0 * k1) + 4.0 * k0;
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;
    float Q = c2 * c2 + c0;
    float R = c2 * c2 * c2 - 3.0 * c2 * c0 + c1 * c1;
    float h = R * R - Q * Q * Q;

    if (h >= 0.0) {
        h = sqrt(h);
        float v = sign(R + h) * pow(abs(R + h), 1.0 / 3.0);    // Cube root
        float u = sign(R - h) * pow(abs(R - h), 1.0 / 3.0);    // Cube root
        vec2 s = vec2((v + u) + 4.0 * c2, (v - u) * sqrt(3.0));
        float y = sqrt(0.5 * (length(s) + s.x));
        float x = 0.5 * s.y / y;
        float r = 2.0 * c1 / (x * x + y * y);
        float t1 =  x - r - k3; t1 = (po < 0.0) ? 2.0 / t1 : t1;
        float t2 = -x - r - k3; t2 = (po < 0.0) ? 2.0 / t2 : t2;
        float t = 1e20;
        if (t1 > 0.0) t = t1;
        if (t2 > 0.0) t = min(t, t2);
        return t;
    }

    float sQ = sqrt(Q);
    float w = sQ * cos(acos(clamp(-R / (sQ * Q), -1.0, 1.0)) / 3.0);
    float d2 = -(w + c2); if (d2 < 0.0) return -1.0;
    float d1 = sqrt(d2);
    float h1 = sqrt(w - 2.0 * c2 + c1 / d1);
    float h2 = sqrt(w - 2.0 * c2 - c1 / d1);
    float t1 = -d1 - h1 - k3; t1 = (po < 0.0) ? 2.0 / t1 : t1;
    float t2 = -d1 + h1 - k3; t2 = (po < 0.0) ? 2.0 / t2 : t2;
    float t3 =  d1 - h2 - k3; t3 = (po < 0.0) ? 2.0 / t3 : t3;
    float t4 =  d1 + h2 - k3; t4 = (po < 0.0) ? 2.0 / t4 : t4;
    float t = 1e20;
    if (t1 > 0.0) t = t1;
    if (t2 > 0.0) t = min(t, t2);
    if (t3 > 0.0) t = min(t, t3);
    if (t4 > 0.0) t = min(t, t4);
    return t;
}

// Normal of the torus of torIntersect at pos, unit length
vec3 torNormal(in vec3 pos, in vec2 tor) {
    return normalize(pos * (dot(pos, pos) - tor.y * tor.y - tor.x * tor.x * vec3(1.0, 1.0, -1.0)));
}

float atom_angle(int k) {
    return in_frag.angle[k / 4][k % 4];
}

vec4 atom_color(int k) {
    return unpackUnorm4x8(in_frag.color[k / 4][k % 4]);
}

void main() {
    // The ray of the fragment, from its point on the quad, which is in front of the ring
    vec3 ro = in_frag.view_pos;
#if ORTHO
    vec3 rd = vec3(0, 0, -1);
#else
    vec3 rd = normalize(in_frag.view_pos);
#endif

    // Frame of the torus (local -> view), scaled by 1 / R: the torus is (1, r / R)
    vec3  nrm   = in_frag.normal;
    vec3  tng   = in_frag.tangent;
    mat3  basis = mat3(tng, cross(nrm, tng), nrm);
    float R     = in_frag.radii.x;
    float inv_R = 1.0 / R;
    vec2  tor   = vec2(1.0, in_frag.radii.y * inv_R);

    vec3 lo = (ro - in_frag.center) * basis * inv_R;    // transpose(basis) * v: view -> local
    vec3 ld = rd * basis;                               // Still unit length

    // Re-origin on a sphere about the torus, a tube radius clear of it, so the quartic sees |lo| ~ 1 and k0 away from 0
    float rb = 1.0 + 2.0 * tor.y;
    float b  = dot(lo, ld);
    float h  = b * b - (dot(lo, lo) - rb * rb);
    if (h < 0.0) {
        discard;
    }
    h = sqrt(h);
    if (-b + h < 0.0) {
        // The sphere is behind the origin
        discard;
    }
    lo += ld * max(-b - h, 0.0);

    float t = torIntersect(lo, ld, tor);
    if (t <= 0.0 || t >= 1e19) {
        discard;
    }
    // Polish the root: two steps by the (signed) distance to the torus
    for (int i = 0; i < 2; ++i) {
        vec3 p = lo + ld * t;
        t += length(vec2(length(p.xy) - tor.x, p.z)) - tor.y;
    }

    vec3 lp          = lo + ld * t;
    vec3 view_coord  = in_frag.center + basis * (lp * R);
    vec3 view_normal = basis * torNormal(lp, tor);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    // Between atoms k and k + 1 in angle, at seg_t from k. The angles increase from 0 at the first atom.
    int   n   = int(in_frag.count_picking.x);
    float phi = atan(lp.y, lp.x);
    if (phi < 0.0) phi += 2.0 * PI;
    int k = 0;
    for (int i = 1; i < n; ++i) {
        if (atom_angle(i) <= phi) k = i;
    }
    float a0 = atom_angle(k);
    float a1 = k + 1 < n ? atom_angle(k + 1) : 2.0 * PI;
    float seg_t = a1 > a0 ? clamp((phi - a0) / (a1 - a0), 0.0, 1.0) : 0.0;
    int   k1    = k + 1 < n ? k + 1 : 0;

    vec4 color;
    if (u_mode == MODE_NEAREST) {
        color = atom_color(seg_t < 0.5 ? k : k1);
    } else if (u_mode == MODE_SMOOTH) {
        float s  = u_sharpness * u_sharpness;
        float kk = mix(1.0, 32.0, s);
        float pa = pow(seg_t, kk);
        float pb = pow(1.0 - seg_t, kk);
        color = mix(atom_color(k), atom_color(k1), pa / (pa + pb));
    } else {
        color = unpackUnorm4x8(u_uniform_color);
    }

    write_fragment(view_coord, in_frag.view_vel, view_normal, color, in_frag.count_picking.y);
}

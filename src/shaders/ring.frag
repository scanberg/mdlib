#version 330 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_shading_language_packing : enable

// Aromatic rings raycast within the quads of ring.vert: a torus in the plane of the ring.
//
// The ray is intersected analytically (the quartic of Inigo Quilez, torRoots), which is only accurate in float when
// its coefficients are of order one: they grow with the fourth power of the distance of the origin from the torus, and
// a hit is decided by terms of the order R^2 r^2. So the ray is taken to the frame of the torus scaled by 1 / R (the
// torus is then (1, r / R) whatever the size of the ring, which also gives the fixed threshold of torRoots the same
// meaning for every ring), and its origin moved to where it enters a sphere about it.
//
// Measured against a double precision reference (rays at rings 2 to 2000 Å away, a quarter of the views near edge on,
// r / R from 0.03 to 0.2), three choices matter in float:
//   - The sphere has the radius 1 + 2 r / R, so the origin stays a tube radius clear of the torus. The bounding sphere
//     itself touches the torus along its outer equator, an origin there makes the constant term of the quartic vanish,
//     and the inversion torRoots then takes (1 / k0) overflows: about 1% of the rays that hit missed.
//   - torRoots inverts the quartic below 0.001 rather than 0.01: the inverted quartic has large coefficients and
//     its discriminant cancels in float, so it is best kept to the cases that need it. Never inverting is far worse.
//   - Two steps along the ray by the distance to the torus polish the root, a tenth of the error where it was largest.
// What remains is about one ray in 10^4 at r / R ~ 0.1 (one in 10^3 for the thinnest tubes), almost all of them
// grazing the outline of a ring seen nearly edge on.
//
// The circle is colored by the atoms of the ring as a bond is by its two: the hit lies between two consecutive atoms in
// angle about the normal, as along a bond between its ends (u_mode). It picks as the ring.
//
// A dashed circle (the circle of a ring of a fused aromatic system, ring.vert gives its dashes over each bond) is the
// torus cut into dashes with flat ends: the gaps are centred on the atoms. All crossings of the ray with the torus are
// found (the near and the far side of the circle), and in each stretch of the ray within the tube it either enters
// within a dash or hits the flat end of the next dash along (dashed_tube_hit). The ends are planes through the axis of
// the torus, intersected exactly. In the gaps nothing is drawn, so neither does the circle pick there.

#ifndef ORTHO
#define ORTHO 0
#endif

#define MODE_NEAREST 0
#define MODE_SMOOTH  1
#define MODE_UNIFORM 2

#define PI 3.14159265358979

#define RING_DASH_FILL 0.5  // Of each dash and its gap, the dash (as ring.vert)

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
    flat uint  flags;
    flat uint  dashes;
    smooth vec3 view_pos;
} in_frag;

// The quad of ring.vert is in front of the ring, so the depth only increases
#ifdef GL_ARB_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif

#pragma EXTRA_SRC

// Ray (ro, rd) against the torus about the z axis of major radius tor.x and tube radius tor.y. rd is unit length.
// Writes the t where the ray crosses the surface to t in increasing order (1e20 past them) and returns how many there
// are: 0, 2 (in at t[0], out at t[1]) or 4 (in, out, in, out). From Inigo Quilez's torIntersect
// (https://iquilezles.org/articles/intersectors), which returns the nearest of them ahead of ro, with the argument of
// acos clamped (rounding takes it past +-1 near tangency, and the NaN would turn every root into a miss) and the
// inversion below 0.001 rather than 0.01 (see above). Each pair of roots is taken where its square root is real.
int torRoots(in vec3 ro, in vec3 rd, in vec2 tor, out float t[4]) {
    t = float[4](1e20, 1e20, 1e20, 1e20);
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

    int count = 0;
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
        t[0] = t1;
        t[1] = t2;
        count = 2;
    } else {
        float sQ = sqrt(Q);
        float w = sQ * cos(acos(clamp(-R / (sQ * Q), -1.0, 1.0)) / 3.0);
        float d2 = -(w + c2);
        if (d2 < 0.0) return 0;
        float d1 = sqrt(d2);
        float e1 = w - 2.0 * c2 + c1 / d1;
        float e2 = w - 2.0 * c2 - c1 / d1;
        if (e1 >= 0.0) {
            float h1 = sqrt(e1);
            float t1 = -d1 - h1 - k3; t1 = (po < 0.0) ? 2.0 / t1 : t1;
            float t2 = -d1 + h1 - k3; t2 = (po < 0.0) ? 2.0 / t2 : t2;
            t[0] = t1;
            t[1] = t2;
            count = 2;
        }
        if (e2 >= 0.0) {
            float h2 = sqrt(e2);
            float t3 =  d1 - h2 - k3; t3 = (po < 0.0) ? 2.0 / t3 : t3;
            float t4 =  d1 + h2 - k3; t4 = (po < 0.0) ? 2.0 / t4 : t4;
            t[count + 0] = t3;
            t[count + 1] = t4;
            count += 2;
        }
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3 - i; ++j) {
            if (t[j + 1] < t[j]) {
                float tmp = t[j]; t[j] = t[j + 1]; t[j + 1] = tmp;
            }
        }
    }
    return count;
}

// Normal of the torus of torRoots at pos, unit length
vec3 torNormal(in vec3 pos, in vec2 tor) {
    return normalize(pos * (dot(pos, pos) - tor.y * tor.y - tor.x * tor.x * vec3(1.0, 1.0, -1.0)));
}

float atom_angle(int k) {
    return in_frag.angle[k / 4][k % 4];
}

vec4 atom_color(int k) {
    return unpackUnorm4x8(in_frag.color[k / 4][k % 4]);
}

// The angle about the normal of a point of the frame of the torus, in [0, 2 pi), increasing from 0 at the first atom
float ring_angle(vec2 p) {
    float phi = atan(p.y, p.x);
    return phi < 0.0 ? phi + 2.0 * PI : phi;
}

// Between atoms k and k + 1 (the bond k) at the angle phi, at seg_t from k: linear in the angle between theirs
void ring_segment(float phi, int n, out int k, out float seg_t) {
    k = 0;
    for (int i = 1; i < n; ++i) {
        if (atom_angle(i) <= phi) k = i;
    }
    float a0 = atom_angle(k);
    float a1 = k + 1 < n ? atom_angle(k + 1) : 2.0 * PI;
    seg_t = a1 > a0 ? clamp((phi - a0) / (a1 - a0), 0.0, 1.0) : 0.0;
}

// DASHED CIRCLES: m dashes over each bond, between the angles of its atoms. Each m-th of the bond holds a dash of
// RING_DASH_FILL of it with half of the gap on either side, so the gaps are centred on the atoms.

// The angle where the dash after the gap at seg_t along bond k begins, going around the circle in the direction of s
// (increasing angle where s > 0): its flat end
float next_dash_angle(int k, float seg_t, float m, float s, int n) {
    float g = 0.5 * (1.0 - RING_DASH_FILL);
    float c = floor(seg_t * m);
    float w = seg_t * m - c;
    float e;
    int   kk = k;
    if (s > 0.0) {
        e = (w < g ? c + g : c + 1.0 + g) / m;
        if (e > 1.0) {
            kk = k + 1 < n ? k + 1 : 0;
            e = g / m;
        }
    } else {
        e = (w > 1.0 - g ? c + 1.0 - g : c - g) / m;
        if (e < 0.0) {
            kk = k > 0 ? k - 1 : n - 1;
            e = 1.0 - g / m;
        }
    }
    float a0 = atom_angle(kk);
    float a1 = kk + 1 < n ? atom_angle(kk + 1) : 2.0 * PI;
    return a0 + (a1 - a0) * e;
}

bool in_dash(float seg_t, float m) {
    float g = 0.5 * (1.0 - RING_DASH_FILL);
    float w = fract(seg_t * m);
    return g <= w && w <= 1.0 - g;
}

// The ray (lo, ld) within the tube from t_in to t_out against the dashes, in the frame of the torus. As
// dashed_cylinder_intersect (bond.frag): if it enters within a dash it hits the side of the tube there (side), otherwise
// it enters in a gap and hits the flat end of the next dash in the direction it goes around the circle, if it reaches
// that end before it leaves the tube (the end is a disc of the tube in a plane through the axis of the torus). The angle
// about the axis only ever increases or decreases along a line, so no other dash comes first. Constant time.
bool dashed_tube_hit(vec3 lo, vec3 ld, vec2 tor, float t_in, float t_out, float m, int n, out float t_hit, out vec3 normal, out bool side) {
    t_hit  = -1.0;
    normal = vec3(0.0);
    side   = true;
    int   k;
    float seg_t;
    if (t_in <= 0.0) {
        // The ray starts within the tube: where it leaves it, if that is in a dash
        vec3 p = lo + ld * t_out;
        ring_segment(ring_angle(p.xy), n, k, seg_t);
        if (t_out <= 0.0 || !in_dash(seg_t, m)) return false;
        t_hit = t_out;
        return true;
    }

    vec3 p = lo + ld * t_in;
    ring_segment(ring_angle(p.xy), n, k, seg_t);
    if (in_dash(seg_t, m)) {
        t_hit = t_in;
        return true;
    }

    // In a gap: the rate of the angle along the ray, up to a positive factor
    float s = p.x * ld.y - p.y * ld.x;
    if (s == 0.0) return false;
    float a   = next_dash_angle(k, seg_t, m, s, n);
    vec2  dir = vec2(cos(a), sin(a));
    vec3  nb  = vec3(-dir.y, dir.x, 0.0);   // Normal of the plane of the end, towards increasing angle
    float den = dot(ld, nb);
    if (den == 0.0) return false;
    float tb = -dot(lo, nb) / den;
    vec3  q  = lo + ld * tb;
    float rho = dot(q.xy, dir);
    float dr  = rho - tor.x;
    if (tb < t_in - 1.0e-4 || rho <= 0.0 || dr * dr + q.z * q.z > tor.y * tor.y) return false;
    t_hit  = max(tb, t_in);
    normal = s > 0.0 ? -nb : nb;
    side   = false;
    return true;
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

    float roots[4];
    int   count = torRoots(lo, ld, tor, roots);
    int   n     = int(in_frag.count_picking.x);

    float t = -1.0;
    bool  side = true;
    vec3  end_normal = vec3(0.0);
    if (in_frag.dashes == 0u) {
        // Solid: the nearest crossing ahead
        for (int i = 0; i < 4; ++i) {
            if (i < count && roots[i] > 0.0) {
                t = roots[i];
                break;
            }
        }
    } else {
        // Dashed: through the tube where the ray is within it, the near side of the circle and then the far side
        float m = float(in_frag.dashes);
        for (int i = 0; i < 4; i += 2) {
            if (i + 1 >= count) break;
            if (dashed_tube_hit(lo, ld, tor, roots[i], roots[i + 1], m, n, t, end_normal, side)) break;
            t = -1.0;
        }
    }
    if (t <= 0.0 || t >= 1e19) {
        discard;
    }

    vec3 lp;
    vec3 view_normal;
    if (side) {
        // Polish the root: two steps by the (signed) distance to the torus
        for (int i = 0; i < 2; ++i) {
            vec3 p = lo + ld * t;
            t += length(vec2(length(p.xy) - tor.x, p.z)) - tor.y;
        }
        lp          = lo + ld * t;
        view_normal = basis * torNormal(lp, tor);
    } else {
        // The flat end of a dash, in a plane through the axis
        lp          = lo + ld * t;
        view_normal = basis * end_normal;
    }
    vec3 view_coord = in_frag.center + basis * (lp * R);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    gl_FragDepth = (clip_coord.z / clip_coord.w) * 0.5 + 0.5;

    // Between atoms k and k + 1 in angle, at seg_t from k. The angles increase from 0 at the first atom.
    int   k;
    float seg_t;
    ring_segment(ring_angle(lp.xy), n, k, seg_t);
    int   k1 = k + 1 < n ? k + 1 : 0;

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

    write_fragment(view_coord, in_frag.view_vel, view_normal, color, in_frag.count_picking.y, in_frag.flags);
}

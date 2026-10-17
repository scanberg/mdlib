// gto_iso.slang - PROTOTYPE kernels for gridless GTO isosurface extraction (see test_gto_iso.c).
//
// The pipeline, one launch per step, all sizes driven by device-side counts:
//
//   iso_octree_level   one thread per node: bounds psi over its 8 children (GROUP bound) and
//                      appends the ones that may contain psi = +iso or psi = -iso
//   iso_brick_hash     brick coordinate -> brick index, open addressing
//   iso_corners        psi at the 5^3 lattice corners of every brick (4^3 cells)
//   iso_cells          Surface Nets vertex per cell and sheet: mass point of the edge crossings
//   iso_quads          one quad per crossing lattice edge, from the four cells around it
//   iso_project        Newton steps onto sheet * psi = iso along the analytic gradient
//   iso_relax          tangential Laplacian toward the face-connected neighbours (Jacobi)
//   iso_triangulate    quad -> two triangles along the shorter diagonal
//   iso_grid           writes an indirect grid from a device-side count
//
// Field representation (built on the host from an MO, see iso_gpu_upload):
//   psi = sum_g sum_n x^i y^j z^k R_gn(|d|^2),  R_gn(rho) = sum_q w_gnq exp(-a_gq rho)
// one group g per (atom, l). R_gn is 1D, so its range over a rho interval is exact from
// precomputed critical intervals (rho_lo, rho_hi, g_lo, g_hi) between which it is monotone.

#include "../../src/shaders/md_gpu.slang"

struct GpuGroup {
    float4 c;       // xyz centre, w: r2 beyond which |psi_g| < point evaluation epsilon
    float4 cut;     // x: r2 beyond which |psi_g| < cut.y (dropped from bounds into the slack)
    uint4  info;    // l, ncart, nprim, alpha offset
    uint4  info2;   // weight offset (ncart * nprim, row major), component offset
};

struct IsoArgs {
    float4    origin_h;       // lattice origin xyz, cell edge h
    uint      num_groups;
    uint      level;          // octree: level of the input nodes (root = 0)
    uint      brick_level;    // octree level whose nodes are bricks (edge 4h)
    uint      cap;            // capacity of the list being appended to
    float     iso;
    float     rel_eps;        // relative widening of float bounds
    float     step_max;       // Newton step clamp
    float     lambda;         // relaxation step
    uint      hash_mask;
    uint      newton_iters;
    uint      grid_mul;       // iso_grid: threads per counted item
    uint      _pad0;
    GpuGroup* groups;
    float*    alpha;
    float*    weights;
    uint4*    comps;          // packed ijk, crit offset, crit count, live
    float4*   crit;
    uint4*    nodes_in;       // x, y, z, sheets
    uint4*    nodes_out;
    uint*     count_in;
    uint*     count_out;
    uint*     grid_out;       // iso_grid output
    uint*     hash_keys;
    uint*     hash_vals;
    float*    corners;        // brick * 125
    uint*     cell_vert;      // brick * 128: sheet * 64 + cell -> vertex, or ~0
    float4*   vpos;           // xyz, w: sheet (+1 / -1)
    float4*   vpos_out;       // relaxation target
    float4*   vnrm;           // outward unit normal, w: |f| / |grad f|
    uint4*    vcell;          // brick, cell, sheet bit, 0
    uint4*    quads;
    uint*     tris;
    uint*     counters;       // 0: vertices, 1: quads, 2: open edges, 3: overflow
};
MD_SHADER_ARGS(IsoArgs);

static const uint INVALID = 0xFFFFFFFFu;

// ------------------------------------------------------------------------------------------
// Intervals
// ------------------------------------------------------------------------------------------

float2 iv_mul(float2 a, float2 b) {
    float p0 = a.x * b.x, p1 = a.x * b.y, p2 = a.y * b.x, p3 = a.y * b.y;
    return float2(min(min(p0, p1), min(p2, p3)), max(max(p0, p1), max(p2, p3)));
}

float iv_mag(float2 a) { return max(abs(a.x), abs(a.y)); }

float2 iv_pow(float2 a, uint n) {
    if (n == 0u) return float2(1.0, 1.0);
    if (n == 1u) return a;
    float pl = a.x, ph = a.y;
    for (uint i = 1u; i < n; ++i) { pl *= a.x; ph *= a.y; }
    if ((n & 1u) != 0u || a.x >= 0.0) return float2(min(pl, ph), max(pl, ph));
    if (a.y <= 0.0) return float2(ph, pl);
    return float2(0.0, max(pl, ph));
}

// ------------------------------------------------------------------------------------------
// The GROUP bound: an enclosure of psi over the box [lo, hi], widened for float rounding
// and for the groups dropped by distance.
// ------------------------------------------------------------------------------------------

float2 psi_bound(float3 lo, float3 hi) {
    IsoArgs a = MD_ARGS;
    float2 sum   = float2(0.0, 0.0);
    float  slack = 0.0;
    float  err   = 0.0;
    for (uint gi = 0u; gi < a.num_groups; ++gi) {
        GpuGroup g = a.groups[gi];
        float3 dlo = lo - g.c.xyz;
        float3 dhi = hi - g.c.xyz;
        float3 nearest = max(max(dlo, -dhi), float3(0.0));
        float3 farthest = max(abs(dlo), abs(dhi));
        float rho_lo = dot(nearest, nearest);
        float rho_hi = dot(farthest, farthest);
        if (rho_lo >= g.cut.x) { slack += g.cut.y; continue; }

        uint l = g.info.x, ncart = g.info.y, nprim = g.info.z;
        float eu[32], ev[32];
        for (uint q = 0u; q < nprim; ++q) {
            float al = a.alpha[g.info.w + q];
            eu[q] = exp(-al * rho_lo);
            ev[q] = exp(-al * rho_hi);
        }
        float2 pw[3][5];
        for (uint n = 0u; n <= l; ++n) {
            pw[0][n] = iv_pow(float2(dlo.x, dhi.x), n);
            pw[1][n] = iv_pow(float2(dlo.y, dhi.y), n);
            pw[2][n] = iv_pow(float2(dlo.z, dhi.z), n);
        }
        for (uint n = 0u; n < ncart; ++n) {
            uint4 cp = a.comps[g.info2.y + n];
            if (cp.w == 0u) continue;
            float gu = 0.0, gv = 0.0, absr = 0.0;
            for (uint q = 0u; q < nprim; ++q) {
                float w = a.weights[g.info2.x + n * nprim + q];
                gu += w * eu[q];
                gv += w * ev[q];
                absr += abs(w) * eu[q];
            }
            float2 R = float2(min(gu, gv), max(gu, gv));
            for (uint k = 0u; k < cp.z; ++k) {
                float4 c = a.crit[cp.y + k];
                if (c.y >= rho_lo && c.x <= rho_hi) R = float2(min(R.x, c.z), max(R.y, c.w));
            }
            uint i = cp.x & 0xFu, j = (cp.x >> 4) & 0xFu, kk = (cp.x >> 8) & 0xFu;
            float2 M = iv_mul(iv_mul(pw[0][i], pw[1][j]), pw[2][kk]);
            float2 t = iv_mul(M, R);
            sum += t;
            err += iv_mag(M) * (absr + iv_mag(R));
        }
    }
    float w = slack + a.rel_eps * err;
    return float2(sum.x - w, sum.y + w);
}

// ------------------------------------------------------------------------------------------
// Point evaluation: value and gradient
// ------------------------------------------------------------------------------------------

float psi_eval(float3 p, out float3 grad) {
    IsoArgs a = MD_ARGS;
    float v = 0.0;
    grad = float3(0.0);
    for (uint gi = 0u; gi < a.num_groups; ++gi) {
        GpuGroup g = a.groups[gi];
        float3 d = p - g.c.xyz;
        float rho = dot(d, d);
        if (rho > g.c.w) continue;
        uint l = g.info.x, ncart = g.info.y, nprim = g.info.z;
        float e[32];
        for (uint q = 0u; q < nprim; ++q) e[q] = exp(-a.alpha[g.info.w + q] * rho);
        float px[5], py[5], pz[5];
        px[0] = 1.0; py[0] = 1.0; pz[0] = 1.0;
        for (uint n = 1u; n <= l; ++n) { px[n] = px[n - 1] * d.x; py[n] = py[n - 1] * d.y; pz[n] = pz[n - 1] * d.z; }
        for (uint n = 0u; n < ncart; ++n) {
            uint4 cp = a.comps[g.info2.y + n];
            if (cp.w == 0u) continue;
            float R = 0.0, R1 = 0.0;
            for (uint q = 0u; q < nprim; ++q) {
                float t = a.weights[g.info2.x + n * nprim + q] * e[q];
                R  += t;
                R1 -= a.alpha[g.info.w + q] * t;
            }
            uint i = cp.x & 0xFu, j = (cp.x >> 4) & 0xFu, k = (cp.x >> 8) & 0xFu;
            float M = px[i] * py[j] * pz[k];
            float3 Mg = float3(i > 0u ? float(i) * px[i - 1] * py[j] * pz[k] : 0.0,
                               j > 0u ? float(j) * px[i] * py[j - 1] * pz[k] : 0.0,
                               k > 0u ? float(k) * px[i] * py[j] * pz[k - 1] : 0.0);
            v += M * R;
            grad += Mg * R + M * (2.0 * R1) * d;
        }
    }
    return v;
}

float psi_value(float3 p) {
    float3 g;
    return psi_eval(p, g);
}

// ------------------------------------------------------------------------------------------
// Small helpers
// ------------------------------------------------------------------------------------------

uint brick_key(uint3 b) { return (b.x | (b.y << 10) | (b.z << 20)) + 1u; }

uint hash_u32(uint x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return x;
}

uint brick_lookup(int3 b) {
    IsoArgs a = MD_ARGS;
    if (any(b < int3(0)) || any(b >= int3(1024))) return INVALID;
    uint key = brick_key(uint3(b));
    uint slot = hash_u32(key) & a.hash_mask;
    for (uint probe = 0u; probe <= a.hash_mask; ++probe) {
        uint k = a.hash_keys[slot];
        if (k == key) return a.hash_vals[slot];
        if (k == 0u) return INVALID;
        slot = (slot + 1u) & a.hash_mask;
    }
    return INVALID;
}

// The vertex of global cell c on sheet s, or INVALID.
uint cell_vertex(int3 c, uint s) {
    IsoArgs a = MD_ARGS;
    if (any(c < int3(0))) return INVALID;
    uint b = brick_lookup(c >> 2);
    if (b == INVALID) return INVALID;
    int3 lc = c & 3;
    return a.cell_vert[b * 128u + s * 64u + uint(lc.x + lc.y * 4 + lc.z * 16)];
}

uint4 brick_node(uint b) { return MD_ARGS.nodes_in[b]; }

float corner_value(uint b, int3 lc) {   // lc in [0, 4]^3
    return MD_ARGS.corners[b * 125u + uint(lc.x + lc.y * 5 + lc.z * 25)];
}

float sheet_sign(uint s) { return s == 0u ? 1.0 : -1.0; }

// ------------------------------------------------------------------------------------------
// Kernels
// ------------------------------------------------------------------------------------------

[shader("compute")]
[numthreads(1, 1, 1)]
void iso_grid(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint threads = a.count_in[0] * a.grid_mul;
    a.grid_out[0] = (threads + 63u) / 64u;
    a.grid_out[1] = 1u;
    a.grid_out[2] = 1u;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_octree_level(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= a.count_in[0]) return;
    uint4 node = a.nodes_in[tid.x];
    float h = a.origin_h.w;
    float child_edge = 4.0 * h * float(1u << (a.brick_level - a.level - 1u));
    for (uint c = 0u; c < 8u; ++c) {
        uint3 ci = node.xyz * 2u + uint3(c & 1u, (c >> 1) & 1u, (c >> 2) & 1u);
        float3 lo = a.origin_h.xyz + float3(ci) * child_edge;
        float2 e = psi_bound(lo, lo + child_edge);
        uint sheets = 0u;
        if (e.x <= a.iso && a.iso <= e.y)   sheets |= 1u;
        if (e.x <= -a.iso && -a.iso <= e.y) sheets |= 2u;
        sheets &= node.w;
        if (sheets == 0u) continue;
        uint idx;
        InterlockedAdd(a.count_out[0], 1u, idx);
        if (idx < a.cap) a.nodes_out[idx] = uint4(ci, sheets);
        else InterlockedAdd(a.counters[3], 1u);
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_brick_hash(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= a.count_in[0]) return;
    uint key = brick_key(a.nodes_in[tid.x].xyz);
    uint slot = hash_u32(key) & a.hash_mask;
    for (uint probe = 0u; probe <= a.hash_mask; ++probe) {
        uint prev;
        InterlockedCompareExchange(a.hash_keys[slot], 0u, key, prev);
        if (prev == 0u || prev == key) { a.hash_vals[slot] = tid.x; return; }
        slot = (slot + 1u) & a.hash_mask;
    }
    InterlockedAdd(a.counters[3], 1u);
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_corners(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint b = tid.x >> 7, k = tid.x & 127u;
    if (b >= a.count_in[0] || k >= 125u) return;
    uint3 lc = uint3(k % 5u, (k / 5u) % 5u, k / 25u);
    uint3 gc = a.nodes_in[b].xyz * 4u + lc;     // integer first: shared corners agree bit for bit
    float3 p = a.origin_h.xyz + float3(gc) * a.origin_h.w;
    a.corners[b * 125u + k] = psi_value(p);
}

static const int3 CUBE_CORNER[8] = { int3(0,0,0), int3(1,0,0), int3(0,1,0), int3(1,1,0),
                                     int3(0,0,1), int3(1,0,1), int3(0,1,1), int3(1,1,1) };
static const uint2 CUBE_EDGE[12] = { uint2(0,1), uint2(2,3), uint2(4,5), uint2(6,7),
                                     uint2(0,2), uint2(1,3), uint2(4,6), uint2(5,7),
                                     uint2(0,4), uint2(1,5), uint2(2,6), uint2(3,7) };

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_cells(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint b = tid.x >> 6, cell = tid.x & 63u;
    if (b >= a.count_in[0]) return;
    uint4 node = a.nodes_in[b];
    int3 lc = int3(cell & 3u, (cell >> 2) & 3u, cell >> 4);
    float cv[8];
    for (uint k = 0u; k < 8u; ++k) cv[k] = corner_value(b, lc + CUBE_CORNER[k]);
    for (uint s = 0u; s < 2u; ++s) {
        uint vi = INVALID;
        if ((node.w & (1u << s)) != 0u) {
            float sg = sheet_sign(s);
            float f[8];
            uint mask = 0u;
            for (uint k = 0u; k < 8u; ++k) { f[k] = sg * cv[k] - a.iso; if (f[k] > 0.0) mask |= 1u << k; }
            if (mask != 0u && mask != 0xFFu) {
                float3 acc = float3(0.0);
                float cnt = 0.0;
                for (uint e = 0u; e < 12u; ++e) {
                    uint e0 = CUBE_EDGE[e].x, e1 = CUBE_EDGE[e].y;
                    if (((mask >> e0) & 1u) == ((mask >> e1) & 1u)) continue;
                    float t = f[e0] / (f[e0] - f[e1]);
                    acc += lerp(float3(CUBE_CORNER[e0]), float3(CUBE_CORNER[e1]), t);
                    cnt += 1.0;
                }
                float3 gcell = float3(node.xyz * 4u + uint3(lc));
                float3 p = a.origin_h.xyz + (gcell + acc / cnt) * a.origin_h.w;
                InterlockedAdd(a.counters[0], 1u, vi);
                if (vi < a.cap) {
                    a.vpos[vi]  = float4(p, sg);
                    a.vcell[vi] = uint4(b, cell, s, 0u);
                } else {
                    InterlockedAdd(a.counters[3], 1u);
                    vi = INVALID;
                }
            }
        }
        a.cell_vert[b * 128u + s * 64u + cell] = vi;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_quads(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint b = tid.x >> 6, cell = tid.x & 63u;
    if (b >= a.count_in[0]) return;
    uint4 node = a.nodes_in[b];
    int3 lc = int3(cell & 3u, (cell >> 2) & 3u, cell >> 4);
    int3 gc = int3(node.xyz * 4u) + lc;
    for (uint s = 0u; s < 2u; ++s) {
        uint v0 = a.cell_vert[b * 128u + s * 64u + cell];
        if (v0 == INVALID) continue;
        float sg = sheet_sign(s);
        float f0 = sg * corner_value(b, lc) - a.iso;
        for (uint ax = 0u; ax < 3u; ++ax) {
            int3 o = int3(ax == 0u ? 1 : 0, ax == 1u ? 1 : 0, ax == 2u ? 1 : 0);
            float f1 = sg * corner_value(b, lc + o) - a.iso;
            if ((f0 > 0.0) == (f1 > 0.0)) continue;
            uint u = (ax + 1u) % 3u, w = (ax + 2u) % 3u;
            int3 eu = int3(u == 0u ? 1 : 0, u == 1u ? 1 : 0, u == 2u ? 1 : 0);
            int3 ew = int3(w == 0u ? 1 : 0, w == 1u ? 1 : 0, w == 2u ? 1 : 0);
            uint v1 = cell_vertex(gc - eu, s);
            uint v2 = cell_vertex(gc - eu - ew, s);
            uint v3 = cell_vertex(gc - ew, s);
            if (v1 == INVALID || v2 == INVALID || v3 == INVALID) { InterlockedAdd(a.counters[2], 1u); continue; }
            uint4 q = (f0 <= 0.0) ? uint4(v0, v3, v2, v1) : uint4(v0, v1, v2, v3);
            uint qi;
            InterlockedAdd(a.counters[1], 1u, qi);
            if (qi < a.cap) a.quads[qi] = q;
            else InterlockedAdd(a.counters[3], 1u);
        }
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_project(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= min(a.count_in[0], a.cap)) return;
    float4 pv = a.vpos[tid.x];
    float3 p = pv.xyz;
    float sg = pv.w;
    float3 n = float3(0.0);
    float res = 0.0;
    for (uint it = 0u; it <= a.newton_iters; ++it) {
        float3 g;
        float f = sg * psi_eval(p, g) - a.iso;
        g *= sg;
        float g2 = dot(g, g);
        if (g2 <= 1e-30) break;
        float gl = sqrt(g2);
        n = -g / gl;
        res = abs(f) / gl;
        if (it == a.newton_iters) break;
        float step = f / g2;
        float len = abs(step) * gl;
        if (len > a.step_max) step *= a.step_max / len;
        p -= step * g;
    }
    a.vpos[tid.x] = float4(p, sg);
    a.vnrm[tid.x] = float4(n, res);
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_relax(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= min(a.count_in[0], a.cap)) return;
    uint4 vc = a.vcell[tid.x];
    uint b = vc.x, s = vc.z;
    uint4 node = a.nodes_in[b];
    int3 lc = int3(vc.y & 3u, (vc.y >> 2) & 3u, vc.y >> 4);
    int3 gc = int3(node.xyz * 4u) + lc;
    float sg = sheet_sign(s);
    float4 pv = a.vpos[tid.x];
    float3 acc = float3(0.0);
    float cnt = 0.0;
    // Face neighbours that share a mesh edge: the shared face has a crossing lattice edge.
    for (uint d = 0u; d < 6u; ++d) {
        uint ax = d >> 1;
        int side = (d & 1u) != 0u ? 1 : 0;              // face at +1 or at 0 along ax
        int3 o = int3(ax == 0u ? 1 : 0, ax == 1u ? 1 : 0, ax == 2u ? 1 : 0);
        int3 u = int3(o.y, o.z, o.x), w = int3(o.z, o.x, o.y);
        int3 base = lc + o * side;
        float f00 = sg * corner_value(b, base) - a.iso;
        float f10 = sg * corner_value(b, base + u) - a.iso;
        float f01 = sg * corner_value(b, base + w) - a.iso;
        float f11 = sg * corner_value(b, base + u + w) - a.iso;
        bool cross = ((f00 > 0.0) != (f10 > 0.0)) || ((f00 > 0.0) != (f01 > 0.0)) ||
                     ((f11 > 0.0) != (f10 > 0.0)) || ((f11 > 0.0) != (f01 > 0.0));
        if (!cross) continue;
        uint vn = cell_vertex(gc + o * (side * 2 - 1), s);
        if (vn == INVALID || vn >= a.cap) continue;
        acc += a.vpos[vn].xyz;
        cnt += 1.0;
    }
    float3 p = pv.xyz;
    if (cnt > 0.0) {
        float3 n = a.vnrm[tid.x].xyz;
        float3 dlt = acc / cnt - p;
        p += a.lambda * (dlt - dot(dlt, n) * n);
    }
    a.vpos_out[tid.x] = float4(p, pv.w);
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_triangulate(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= min(a.count_in[0], a.cap)) return;
    uint4 q = a.quads[tid.x];
    float3 p0 = a.vpos[q.x].xyz, p1 = a.vpos[q.y].xyz, p2 = a.vpos[q.z].xyz, p3 = a.vpos[q.w].xyz;
    uint o = tid.x * 6u;
    if (dot(p0 - p2, p0 - p2) <= dot(p1 - p3, p1 - p3)) {
        a.tris[o + 0] = q.x; a.tris[o + 1] = q.y; a.tris[o + 2] = q.z;
        a.tris[o + 3] = q.x; a.tris[o + 4] = q.z; a.tris[o + 5] = q.w;
    } else {
        a.tris[o + 0] = q.x; a.tris[o + 1] = q.y; a.tris[o + 2] = q.w;
        a.tris[o + 3] = q.y; a.tris[o + 4] = q.z; a.tris[o + 5] = q.w;
    }
}

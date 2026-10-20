// gto_iso_gpu.slang - gridless isosurface meshes of GTO fields on md_gpu (see md_gto_iso.h, md_gto_iso_gpu.inl)
//
// A FIELD SET is one basis, a coefficient matrix C (AOs x functions) and up to 16 output fields, each
// a short list of terms over the functions:
//     f_k(r) = sum_t w_t a_t(r) b_t(r)    (b_t = 1 for a linear term: an orbital)
//     a(r)   = sum_i C[i][a] phi_i(r)
// One run extracts every field: the octree and the corner values are shared, so the AOs at a point
// are evaluated once however many fields there are (the 16 orbitals of a grid, or one density).
//
// The pipeline (md_gto_iso_gpu.inl issues it):
//   iso_octree      one group per node: AO ranges and centre values of every shell over the node,
//                   the functions as products with C (interval extension and centred form,
//                   intersected), the terms per field. Keeps the node when some field may cross its
//                   iso value; 2 bits per field (positive, negative sheet). Rigorous in float: every
//                   rounding step is bounded and added.
//   iso_hash        brick coordinate -> brick index
//   iso_corners     one group per brick: the field values at its 5^3 lattice points, a fused
//                   register-tiled GEMM (AO values x C) with the terms applied in the epilogue
//   iso_cells       Surface Nets vertex per cell, field and sheet; a brick's vertices are contiguous
//   iso_quads       one quad per crossing lattice edge
//   iso_eval        value and gradient at the vertices, the same GEMM with derivative channels
//   iso_corners_orbital, iso_eval_orbital
//                   the same for sets of orbitals (one function per field): one thread per point,
//                   every field's function in registers - a GEMM tile would be mostly empty
//   iso_newton      one Newton step (or the final normals) from iso_eval's values; converged
//                   vertices are flagged and skipped by the evaluators until they move again
//   iso_relax       tangential Laplacian toward the face-connected neighbours
//   iso_bucket_*    vertices and triangles regrouped per field and sheet
//   iso_grid        indirect grids from device-side counts

#include "../md_gpu.slang"

struct IsoShell {
    float4 c;       // xyz centre, w: r^2 beyond which every AO of the shell is below the evaluation epsilon
    uint4  info;    // l, ncart, nprim, prim offset
    uint4  info2;   // first AO, crit offset, crit count, -
};

struct IsoField {
    float iso;      // > 0; the sheets are f = +iso and f = -iso
    float slack;    // added to every bound (a density matrix's splitting residual)
    uint  term_off;
    uint  term_count;
};

struct IsoTerm {
    float w;
    uint  a;
    uint  b;        // INVALID: linear term (w a)
    uint  _pad;
};

struct IsoArgs {
    float4    origin_h;       // lattice origin xyz, cell edge h
    uint      num_shells;
    uint      num_ao;
    uint      num_fn;
    uint      fn_stride;      // row stride of C (a multiple of 64)
    uint      num_fields;
    uint      level;          // iso_octree: level of the input nodes (root = 0)
    uint      brick_level;    // the level whose nodes are bricks (edge 4h)
    uint      cap;            // capacity of the list being appended to
    uint      hash_mask;
    uint      mode;           // iso_cells: 0 count, 1 emit. iso_newton: 0 step, 1 normals only
    uint      count;          // host-known item count (vertices, quads)
    uint      grid_mul;       // iso_grid: groups (or threads) per counted item
    uint      grid_div;       // iso_grid: items per group
    float     step_max;       // Newton step clamp
    float     lambda;         // relaxation step
    float     tol;            // Newton: converged below this distance
    IsoShell* shells;
    float2*   prims;          // alpha, contraction coefficient (normalized primitives)
    float4*   crit;           // per shell radial function: rho_lo, rho_hi, R_lo, R_hi
    float*    ao_nrm;         // per AO: Cartesian normalization
    uint*     ao_ijk;         // per AO: i | j << 4 | k << 8
    float*    C;              // [num_ao][fn_stride]
    IsoTerm*  terms;
    IsoField* fields;
    uint4*    nodes_in;       // x, y, z, sheet mask (2 bits per field)
    uint4*    nodes_out;
    uint*     count_in;
    uint*     count_out;
    uint*     grid_out;
    uint4*    bricks;         // x, y, z, sheet mask
    uint*     brick_slot;     // per brick: its first field slot (one slot per field with a sheet bit)
    uint*     hash_keys;
    uint*     hash_vals;
    float*    corners;        // [slot][125]
    uint*     cell_vert;      // [slot][2][64]: vertex index or INVALID
    float4*   vpos;           // xyz, w: sheet sign
    float4*   vpos_out;
    float4*   vnrm;           // outward unit normal, w: |f| / |grad f|
    float4*   vval;           // iso_eval: gradient xyz, value w
    uint4*    vinfo;          // brick, cell, field, sheet
    uint4*    quads;
    uint*     counters;       // 0 vertices, 1 quads, 2 open edges, 3 overflow, 4 -, 5 corner points evaluated
    uint*     bucket;         // [0, 32) vertex counts, [32, 64) quad counts, [64, 96) vertex offsets,
                              // [96, 128) quad offsets, [128, 160) vertex cursors, [160, 192) quad cursors
    uint*     remap;          // old vertex -> bucketed vertex
    uint*     vdone;          // per vertex: converged, no more evaluations until it moves
    float4*   out_pos;
    float4*   out_nrm;
    uint*     out_idx;
    float*    fn_tail;        // per function: bound of what the shells dropped from a box can add
    uint*     plist;          // corners: the lattice points to evaluate, brick << 7 | point
};
MD_SHADER_ARGS(IsoArgs);

static const uint  INVALID = 0xFFFFFFFFu;
static const float U = 1.0 / 4194304.0;     // 2^-22: four float ulps, generous for fast-math backends
static const float EXP_REL = 4.0e-6;        // relative error allowed for exp() on any backend

uint  node_index(uint3 gid) { return gid.y * 65535u + gid.x; }
uint  field_bits(uint mask) { return (mask | (mask >> 1)) & 0x55555555u; }     // bit 2f: field f active
uint  field_slot(uint mask, uint f) { return countbits(field_bits(mask) & ((1u << (2u * f)) - 1u)); }

// ------------------------------------------------------------------------------------------
// Intervals
// ------------------------------------------------------------------------------------------

float2 iv_mul(float2 a, float2 b) {
    float p0 = a.x * b.x, p1 = a.x * b.y, p2 = a.y * b.x, p3 = a.y * b.y;
    return float2(min(min(p0, p1), min(p2, p3)), max(max(p0, p1), max(p2, p3)));
}

float iv_mag(float2 a) { return max(abs(a.x), abs(a.y)); }

float2 iv_sq(float2 a) {
    float l = a.x * a.x, h = a.y * a.y;
    if (a.x >= 0.0) return float2(l, h);
    if (a.y <= 0.0) return float2(h, l);
    return float2(0.0, max(l, h));
}

float2 iv_pow(float2 a, uint n) {
    if (n == 0u) return float2(1.0, 1.0);
    if (n == 1u) return a;
    float pl = a.x, ph = a.y;
    for (uint i = 1u; i < n; ++i) { pl *= a.x; ph *= a.y; }
    if ((n & 1u) != 0u || a.x >= 0.0) return float2(min(pl, ph), max(pl, ph));
    if (a.y <= 0.0) return float2(ph, pl);
    return float2(0.0, max(pl, ph));
}

float ipow(float x, uint n) {
    float r = 1.0;
    for (uint i = 0u; i < n; ++i) r *= x;
    return r;
}

float box_dist2(float3 c, float3 lo, float3 hi) {
    float3 d = max(max(lo - c, c - hi), float3(0.0));
    return dot(d, d);
}

// ------------------------------------------------------------------------------------------
// Group-wide compaction: the exclusive prefix count of a flag over the group (as gto_common's
// scan: SubgroupId, not tid / lanes, since a subgroup may be partly populated on Vulkan)
// ------------------------------------------------------------------------------------------

groupshared uint g_scan_wave[64];
groupshared uint g_scan_total;

uint iso_wave_index(uint tid, uint lanes) {
    __target_switch {
    case spirv: return spirv_asm { result:$$uint = OpLoad builtin(SubgroupId:uint); };
    default:    return tid / lanes;
    }
}

uint iso_num_waves(uint group_size, uint lanes) {
    __target_switch {
    case spirv: return spirv_asm { result:$$uint = OpLoad builtin(NumSubgroups:uint); };
    default:    return (group_size + lanes - 1u) / lanes;
    }
}

// The exclusive prefix sum of v over the group, and the total. Every thread of the group calls it
// (barriers inside, the first after every thread has passed its earlier groupshared reads).
uint iso_group_prefix_sum(uint v, uint tid, uint group_size, out uint total) {
    uint lanes = WaveGetLaneCount();
    uint w = iso_wave_index(tid, lanes);
    uint pre = WavePrefixSum(v);
    uint cnt = WaveActiveSum(v);
    if (WaveIsFirstLane()) g_scan_wave[w] = cnt;
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        uint nw = iso_num_waves(group_size, lanes);
        uint run = 0u;
        for (uint k = 0u; k < nw; ++k) { uint c = g_scan_wave[k]; g_scan_wave[k] = run; run += c; }
        g_scan_total = run;
    }
    GroupMemoryBarrierWithGroupSync();
    uint r = g_scan_wave[w] + pre;
    total = g_scan_total;
    GroupMemoryBarrierWithGroupSync();
    return r;
}

uint iso_group_prefix(bool flag, uint tid, uint group_size, out uint total) {
    return iso_group_prefix_sum(flag ? 1u : 0u, tid, group_size, total);
}

// ------------------------------------------------------------------------------------------
// Indirect grids: grid_mul groups per counted item, or one group per grid_div items, split over y
// beyond 65535 groups (node_index undoes it)
// ------------------------------------------------------------------------------------------

[shader("compute")]
[numthreads(1, 1, 1)]
void iso_grid(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint n = a.count_in[0];
    if (a.cap != 0u) n = min(n, a.cap);
    uint groups = a.grid_div > 1u ? (n + a.grid_div - 1u) / a.grid_div : n * a.grid_mul;
    a.grid_out[0] = min(groups, 65535u);
    a.grid_out[1] = (groups + 65534u) / 65535u;
    a.grid_out[2] = 1u;
}

// ------------------------------------------------------------------------------------------
// Octree: one group of 64 threads per node
// ------------------------------------------------------------------------------------------

#define O_BLOCK 64          // shells per pass (one per thread)
#define O_AO_CAP 960        // AOs of 64 shells (l <= 4: 15 components)
#define O_FN_CAP 1024       // functions
#define O_QMAX   4          // tiles of 64 functions per pass over the shells

groupshared float  o_mid[O_AO_CAP];
groupshared float  o_rad[O_AO_CAP];
groupshared float  o_pc[O_AO_CAP];
groupshared float  o_dev[O_AO_CAP];
groupshared uint   o_ao[O_AO_CAP];
groupshared float  o_acc[5 * 64];
groupshared float2 o_A[O_FN_CAP];
groupshared uint   o_mask;

// AO data of one shell over the box: ranges (mid, rad), values at the centre and the largest
// deviation from them, written at slot 'off'.
void octree_shell(uint s, float3 lo, float3 hi, float3 cen, uint off) {
    IsoArgs a = MD_ARGS;
    IsoShell sh = a.shells[s];
    uint l = sh.info.x, ncart = sh.info.y, nprim = sh.info.z, poff = sh.info.w;
    float3 dlo = lo - sh.c.xyz, dhi = hi - sh.c.xyz;
    float3 nearest = max(max(dlo, -dhi), float3(0.0));
    float3 farthest = max(abs(dlo), abs(dhi));
    // rho range widened for the rounding of the distances (the exponentials are monotone in rho)
    float rlo = dot(nearest, nearest) * (1.0 - 8.0 * U);
    float rhi = dot(farthest, farthest) * (1.0 + 8.0 * U) + 1.0e-30;
    float3 dc = cen - sh.c.xyz;
    float rc = dot(dc, dc);
    float gu = 0.0, gv = 0.0, au = 0.0, Rc = 0.0, ac = 0.0;
    for (uint q = 0u; q < nprim; ++q) {
        float2 pr = a.prims[poff + q];
        float eu = exp(-pr.x * rlo), ev = exp(-pr.x * rhi), ec = exp(-pr.x * rc);
        gu += pr.y * eu;
        gv += pr.y * ev;
        au += abs(pr.y) * eu;
        Rc += pr.y * ec;
        ac += abs(pr.y) * ec;
    }
    float2 R = float2(min(gu, gv), max(gu, gv));
    for (uint k = 0u; k < sh.info2.z; ++k) {
        float4 c = a.crit[sh.info2.y + k];
        if (c.y >= rlo && c.x <= rhi) R = float2(min(R.x, c.z), max(R.y, c.w));
    }
    float eR = au * (EXP_REL + float(nprim + 4u) * U) + 1.0e-35;
    R += float2(-eR, eR);
    float2 pw[3][5];
    for (uint n = 0u; n <= l; ++n) {
        pw[0][n] = iv_pow(float2(dlo.x, dhi.x), n);
        pw[1][n] = iv_pow(float2(dlo.y, dhi.y), n);
        pw[2][n] = iv_pow(float2(dlo.z, dhi.z), n);
    }
    uint ao0 = sh.info2.x;
    for (uint n = 0u; n < ncart; ++n) {
        uint ijk = a.ao_ijk[ao0 + n];
        uint i = ijk & 0xFu, j = (ijk >> 4) & 0xFu, kk = (ijk >> 8) & 0xFu;
        float nrm = a.ao_nrm[ao0 + n];
        float2 M = iv_mul(iv_mul(pw[0][i], pw[1][j]), pw[2][kk]);
        float2 P = iv_mul(M, R) * nrm;
        if (nrm < 0.0) P = P.yx;
        float pm = iv_mag(P);
        P += float2(-8.0 * U * pm, 8.0 * U * pm);
        float pc = nrm * ipow(dc.x, i) * ipow(dc.y, j) * ipow(dc.z, kk) * Rc;
        uint slot = off + n;
        o_mid[slot] = 0.5 * (P.x + P.y);
        o_rad[slot] = 0.5 * (P.y - P.x) * (1.0 + 4.0 * U) + U * abs(P.x + P.y);
        o_pc[slot]  = pc;
        o_dev[slot] = max(max(P.y - pc, pc - P.x), 0.0) * (1.0 + 4.0 * U);
        o_ao[slot]  = ao0 + n;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_octree(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint ni = node_index(gid);
    if (ni >= min(a.count_in[0], a.cap)) return;          // uniform per group
    uint4 node = a.nodes_in[ni];

    float h = a.origin_h.w;
    float edge = 4.0 * h * float(1u << (a.brick_level - a.level));
    float3 lo = a.origin_h.xyz + float3(node.xyz) * edge;
    float3 hi = lo + edge;
    // the lattice points inside are computed as origin + index * h in float: cover their rounding
    float slop = 1.0e-5 * edge + 4.0 * U * (abs(a.origin_h.x) + abs(a.origin_h.y) + abs(a.origin_h.z) + 3.0 * edge * float(node.x + node.y + node.z + 1u));
    lo -= slop;
    hi += slop;
    float3 cen = 0.5 * (lo + hi);

    if (tid == 0u) o_mask = 0u;
    // Threads over the functions: below 64 functions, tn2 lanes per function share its AOs out in
    // slices; from 64 on, a lane per function of each tile of 64. Up to O_QMAX tiles are accumulated
    // in one pass over the shells, so the shells' ranges over the node are made once per pass.
    uint nf = a.num_fn;
    uint ntile = (nf + 63u) / 64u;
    uint tn2 = 1u;
    while (tn2 < min(64u, nf)) tn2 <<= 1u;
    uint fl = tid % tn2, slice = tid / tn2, nslice = 64u / tn2;
    for (uint q0 = 0u; q0 < ntile; q0 += O_QMAX) {
        float m1[O_QMAX], r1[O_QMAX], m2[O_QMAX], r2[O_QMAX], mg[O_QMAX];
        [ForceUnroll] for (int q = 0; q < O_QMAX; ++q) { m1[q] = 0.0; r1[q] = 0.0; m2[q] = 0.0; r2[q] = 0.0; mg[q] = 0.0; }
        for (uint s0 = 0u; s0 < a.num_shells; s0 += O_BLOCK) {
            // the block's shells that reach the node; the others are below the evaluation epsilon
            // over all of it, which the bound covers with fn_tail
            uint s = min(s0 + tid, a.num_shells - 1u);       // clamped: llvmpipe loads for inactive lanes
            float4 sc = a.shells[s].c;
            uint snc = a.shells[s].info.y;
            bool keep = s0 + tid < a.num_shells && box_dist2(sc.xyz, lo, hi) <= sc.w * (1.0 + 16.0 * U);
            uint nao;
            uint off = iso_group_prefix_sum(keep ? snc : 0u, tid, 64u, nao);   // barriers inside: the previous block is consumed
            if (keep) octree_shell(s, lo, hi, cen, off);
            GroupMemoryBarrierWithGroupSync();
            for (uint k = slice; k < nao; k += nslice) {
                uint row = min(o_ao[k], a.num_ao - 1u) * a.fn_stride;
                float mid = o_mid[k], rad = o_rad[k], pc = o_pc[k], dev = o_dev[k];
                float gm = abs(mid) + rad + abs(pc) + dev;
                [ForceUnroll] for (int q = 0; q < O_QMAX; ++q) {
                    if (q0 + uint(q) < ntile) {                  // uniform
                        float c = a.C[row + (q0 + uint(q)) * 64u + fl];   // padded columns are zero
                        float ca = abs(c);
                        m1[q] += c * mid;
                        r1[q] += ca * rad;
                        m2[q] += c * pc;
                        r2[q] += ca * dev;
                        mg[q] += ca * gm;
                    }
                }
            }
        }
        [ForceUnroll] for (int q = 0; q < O_QMAX; ++q) {
            uint t0 = (q0 + uint(q)) * 64u;
            if (t0 < nf) {                                       // uniform
                float sm[5] = {m1[q], r1[q], m2[q], r2[q], mg[q]};
                if (nslice > 1u) {                               // one tile: sum the slices
                    GroupMemoryBarrierWithGroupSync();
                    [ForceUnroll] for (int c = 0; c < 5; ++c) o_acc[c * 64 + int(tid)] = sm[c];
                    GroupMemoryBarrierWithGroupSync();
                    if (slice == 0u) {
                        for (uint sl = 1u; sl < nslice; ++sl) {
                            [ForceUnroll] for (int c = 0; c < 5; ++c) sm[c] += o_acc[uint(c) * 64u + sl * tn2 + fl];
                        }
                    }
                }
                uint col = t0 + fl;
                if (slice == 0u && col < nf) {
                    float e = float(a.num_ao + 8u) * 2.0 * U * sm[4] + a.fn_tail[col] + 1.0e-35;
                    float2 i1 = float2(sm[0] - sm[1] - e, sm[0] + sm[1] + e);
                    float2 i2 = float2(sm[2] - sm[3] - e, sm[2] + sm[3] + e);
                    float2 A = float2(max(i1.x, i2.x), min(i1.y, i2.y));
                    if (A.x > A.y) A = float2(min(i1.x, i2.x), max(i1.y, i2.y));   // only by rounding; keep the hull
                    o_A[col] = A;
                }
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();

    if (tid < a.num_fields) {
        uint fbits = (node.w >> (2u * tid)) & 3u;
        if (fbits != 0u) {
            IsoField fd = a.fields[min(tid, 15u)];
            float2 rho = float2(0.0, 0.0);
            float rmag = 0.0;
            for (uint t = 0u; t < fd.term_count; ++t) {
                IsoTerm tm = a.terms[fd.term_off + t];
                float2 q;
                if (tm.b == INVALID)  q = o_A[tm.a];
                else if (tm.a == tm.b) q = iv_sq(o_A[tm.a]);
                else                  q = iv_mul(o_A[tm.a], o_A[tm.b]);
                float2 wq = q * tm.w;
                if (tm.w < 0.0) wq = wq.yx;
                rho += wq;
                rmag += abs(tm.w) * iv_mag(q);
            }
            float e = fd.slack + float(fd.term_count + 8u) * 2.0 * U * rmag + 1.0e-35;
            rho += float2(-e, e);
            uint bits = 0u;
            if (rho.x <= fd.iso && rho.y >= fd.iso)   bits |= 1u;
            if (rho.x <= -fd.iso && rho.y >= -fd.iso) bits |= 2u;
            bits &= fbits;
            if (bits != 0u) InterlockedOr(o_mask, bits << (2u * tid));
        }
    }
    GroupMemoryBarrierWithGroupSync();
    if (a.level == a.brick_level && a.mode == 1u) {
        // the bricks keep their slots, a rejected one with no sheet bits: siblings stay together
        // (what the point list of the corners relies on to keep its groups compact)
        if (tid == 0u) {
            a.nodes_out[ni] = uint4(node.xyz, o_mask);
            if (ni == 0u) a.count_out[0] = min(a.count_in[0], a.cap);
        }
        return;
    }
    if (tid == 0u && o_mask != 0u) {
        if (a.level < a.brick_level) {
            uint idx;
            InterlockedAdd(a.count_out[0], 8u, idx);
            if (idx + 8u <= a.cap) {
                for (uint c = 0u; c < 8u; ++c) {
                    a.nodes_out[idx + c] = uint4(node.xyz * 2u + uint3(c & 1u, (c >> 1) & 1u, (c >> 2) & 1u), o_mask);
                }
            } else {
                InterlockedAdd(a.counters[3], 1u);
            }
        } else {
            uint idx;
            InterlockedAdd(a.count_out[0], 1u, idx);
            if (idx < a.cap) a.nodes_out[idx] = uint4(node.xyz, o_mask);
            else InterlockedAdd(a.counters[3], 1u);
        }
    }
}

// ------------------------------------------------------------------------------------------
// Bricks: hash and lookups
// ------------------------------------------------------------------------------------------

uint brick_key(uint3 b) { return (b.x | (b.y << 10) | (b.z << 20)) + 1u; }

uint hash_u32(uint x) {
    x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16;
    return x;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_hash(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    if (tid.x >= a.count || a.bricks[tid.x].w == 0u) return;      // a rejected sibling
    uint key = brick_key(a.bricks[tid.x].xyz);
    uint slot = hash_u32(key) & a.hash_mask;
    for (uint probe = 0u; probe <= a.hash_mask; ++probe) {
        uint prev;
        InterlockedCompareExchange(a.hash_keys[slot], 0u, key, prev);
        if (prev == 0u || prev == key) { a.hash_vals[slot] = tid.x; return; }
        slot = (slot + 1u) & a.hash_mask;
    }
    InterlockedAdd(a.counters[3], 1u);
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

// The vertex of global cell c for field f and sheet s, or INVALID.
uint cell_vertex(int3 c, uint f, uint s) {
    IsoArgs a = MD_ARGS;
    if (any(c < int3(0))) return INVALID;
    uint b = brick_lookup(c >> 2);
    // loads from brick 0 when there is none: llvmpipe also loads for inactive lanes
    uint bs = b == INVALID ? 0u : b;
    uint mask = a.bricks[bs].w;
    if (b == INVALID || ((mask >> (2u * f + s)) & 1u) == 0u) return INVALID;
    uint slot = a.brick_slot[bs] + field_slot(mask, f);
    int3 lc = c & 3;
    return a.cell_vert[slot * 128u + s * 64u + uint(lc.x + lc.y * 4 + lc.z * 16)];
}

float corner_value(uint slot, int3 lc) {   // lc in [0, 4]^3
    return MD_ARGS.corners[slot * 125u + uint(lc.x + lc.y * 5 + lc.z * 25)];
}

// Shared corners: a lattice point on a brick's high face (a local coordinate 4) is the low-face point
// of the neighbour across it. When that neighbour is present with every field of this brick, the
// value is taken from it (iso_corners_fill) instead of evaluated twice. Returns the neighbour, or
// INVALID when the point is the brick's own to evaluate; *owner_p: the point's index there.
uint corner_owner(uint4 brick, uint p, out uint owner_p) {
    IsoArgs a = MD_ARGS;
    uint k = min(p, 124u);
    uint3 lc = uint3(k % 5u, (k / 5u) % 5u, k / 25u);
    uint3 dd = uint3(lc.x == 4u ? 1u : 0u, lc.y == 4u ? 1u : 0u, lc.z == 4u ? 1u : 0u);
    uint3 ol = lc - 4u * dd;
    owner_p = ol.x + ol.y * 5u + ol.z * 25u;
    if (all(dd == uint3(0u))) return INVALID;
    uint nb = brick_lookup(int3(brick.xyz + dd));
    uint nmask = a.bricks[nb == INVALID ? 0u : nb].w;      // brick 0 when absent: llvmpipe loads anyway
    uint need = field_bits(brick.w);
    return (nb != INVALID && (field_bits(nmask) & need) == need) ? nb : INVALID;
}

// ------------------------------------------------------------------------------------------
// Values at points: the AO values of a group of points times C, register tiled, terms in the
// epilogue. iso_corners: 128 points (a brick's 125 corners), value only. iso_eval: 64 vertices,
// value and gradient.
//
// Per group: the shells that reach the group's box are listed 256 at a time, packed into chunks
// of at most 16 AOs (whole shells), and for each chunk the AO values (and derivatives) of every
// point and the chunk's rows of C go to groupshared memory; each thread accumulates 8 x 4
// (corners) or 4 x 4 x 4 channels (vertices) products. Functions are taken 64 at a time.
// ------------------------------------------------------------------------------------------

#define E_LIST 256
groupshared uint  e_list[E_LIST];
groupshared uint  e_list_nc[E_LIST];   // the listed shells' AO counts
groupshared uint  e_nlist;
groupshared uint  e_next;
groupshared uint  e_ch_shell[16];   // chunk: shells (global index)
groupshared uint  e_ch_ao[16];      // chunk: per AO: shell slot | component << 8
groupshared uint  e_ch_first[16];   // chunk: per shell, its first AO slot
groupshared uint  e_ch_nsh;
groupshared uint  e_ch_nao;
groupshared float e_c[16 * 64];     // chunk rows of C
groupshared float e_buf[4096];      // AO values (and radial parts); the epilogue's staging
groupshared float e_out[128 * 16];  // corners: per point and field. vertices: per point, 4 channels
groupshared float4 e_p[64];         // vertices: positions
groupshared float4 e_lo;
groupshared float4 e_hi;

// Lists the shells of block s0 that reach [lo, hi] (all threads; returns the count via e_nlist).
void list_shells(uint s0, float3 lo, float3 hi, uint tid) {
    IsoArgs a = MD_ARGS;
    // indices clamped into range here and below: llvmpipe also loads for inactive lanes
    uint s = min(s0 + tid, a.num_shells - 1u);
    float4 c = a.shells[s].c;
    bool reach = s0 + tid < a.num_shells && box_dist2(c.xyz, lo, hi) <= c.w;
    uint nc = a.shells[s].info.y;
    uint total;
    uint at = iso_group_prefix(reach, tid, E_LIST, total);
    if (reach) { e_list[at] = s; e_list_nc[at] = nc; }
    if (tid == 0u) { e_nlist = total; e_next = 0u; }
    GroupMemoryBarrierWithGroupSync();
}

// Packs the next chunk of listed shells (thread 0). e_ch_nao == 0: the list is done.
void next_chunk(uint tid) {
    IsoArgs a = MD_ARGS;
    if (tid == 0u) {
        uint nsh = 0u, nao = 0u;
        while (e_next < e_nlist) {
            uint s = e_list[e_next];
            uint nc = e_list_nc[e_next];
            if (nao + nc > 16u) break;
            for (uint n = 0u; n < nc; ++n) e_ch_ao[nao + n] = nsh | (n << 8);
            e_ch_shell[nsh] = s;
            e_ch_first[nsh] = nao;
            nsh++;
            nao += nc;
            e_next++;
        }
        e_ch_nsh = nsh;
        e_ch_nao = nao;
    }
    GroupMemoryBarrierWithGroupSync();
}

// The chunk's rows of C for functions [f0, f0 + 64).
void load_c(uint f0, uint tid) {
    IsoArgs a = MD_ARGS;
    for (uint it = tid; it < 16u * 64u; it += 256u) {
        uint k = it >> 6, f = it & 63u;
        float v = 0.0;
        if (k < e_ch_nao) {
            uint sh = e_ch_shell[e_ch_ao[k] & 0xFFu];
            uint ao = min(a.shells[min(sh, a.num_shells - 1u)].info2.x + (e_ch_ao[k] >> 8), a.num_ao - 1u);
            if (f0 + f < a.num_fn) v = a.C[ao * a.fn_stride + f0 + f];
        }
        e_c[it] = v;
    }
}

float3 corner_point(uint4 brick, uint p) {
    IsoArgs a = MD_ARGS;
    uint k = min(p, 124u);
    uint3 lc = uint3(k % 5u, (k / 5u) % 5u, k / 25u);
    uint3 gc = brick.xyz * 4u + lc;     // integer first: shared corners agree bit for bit
    return a.origin_h.xyz + float3(gc) * a.origin_h.w;
}

// Adds field f's terms over the functions [f0, f0 + 64) at one point. X(p, fn) reads the staged
// values. Linear terms use the function alone.
float field_terms(uint f, uint f0, uint row) {
    IsoArgs a = MD_ARGS;
    IsoField fd = a.fields[f];
    float v = 0.0;
    for (uint t = 0u; t < fd.term_count; ++t) {
        IsoTerm tm = a.terms[fd.term_off + t];
        if (tm.a < f0 || tm.a >= f0 + 64u) continue;
        float xa = e_buf[row * 64u + tm.a - f0];
        if (tm.b == INVALID) { v += tm.w * xa; continue; }
        float xb = e_buf[row * 64u + tm.b - f0];
        v += tm.w * xa * xb;
    }
    return v;
}

[shader("compute")]
[numthreads(16, 16, 1)]
void iso_corners(uint3 gid : SV_GroupID, uint3 gt : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    uint4 brick = a.bricks[b];
    float h = a.origin_h.w;
    float3 blo = a.origin_h.xyz + float3(brick.xyz * 4u) * h;
    float3 bhi = blo + 4.0 * h;
    uint tx = gt.x, ty = gt.y;

    for (uint it = tid; it < 128u * 16u; it += 256u) e_out[it] = 0.0;

    for (uint f0 = 0u; f0 < a.num_fn; f0 += 64u) {
        float X[8][4];
        for (int i = 0; i < 8; ++i) for (int j = 0; j < 4; ++j) X[i][j] = 0.0;
        for (uint s0 = 0u; s0 < a.num_shells; s0 += E_LIST) {
            GroupMemoryBarrierWithGroupSync();
            list_shells(s0, blo, bhi, tid);
            for (;;) {
                next_chunk(tid);
                uint nsh = e_ch_nsh, nao = e_ch_nao;
                if (nao == 0u) break;                       // uniform: read after the barrier
                // radial parts per (shell, point) into e_buf[2048 + s * 128 + p]
                for (uint it = tid; it < nsh * 128u; it += 256u) {
                    uint s = it >> 7, p = it & 127u;
                    IsoShell sh = a.shells[min(e_ch_shell[s], a.num_shells - 1u)];
                    float3 d = corner_point(brick, p) - sh.c.xyz;
                    float r2 = dot(d, d);
                    float R = 0.0;
                    for (uint q = 0u; q < sh.info.z; ++q) {
                        float2 pr = a.prims[sh.info.w + q];
                        R += pr.y * exp(-pr.x * r2);
                    }
                    e_buf[2048u + it] = R;
                }
                load_c(f0, tid);
                GroupMemoryBarrierWithGroupSync();
                for (uint it = tid; it < nao * 128u; it += 256u) {
                    uint k = it >> 7, p = it & 127u;
                    uint s = e_ch_ao[k] & 0xFFu, n = e_ch_ao[k] >> 8;
                    IsoShell sh = a.shells[min(e_ch_shell[s], a.num_shells - 1u)];
                    uint ao = sh.info2.x + n;
                    uint ijk = a.ao_ijk[ao];
                    float3 d = corner_point(brick, p) - sh.c.xyz;
                    float M = ipow(d.x, ijk & 0xFu) * ipow(d.y, (ijk >> 4) & 0xFu) * ipow(d.z, (ijk >> 8) & 0xFu);
                    e_buf[it] = a.ao_nrm[ao] * M * e_buf[2048u + s * 128u + p];
                }
                GroupMemoryBarrierWithGroupSync();
                for (uint k = 0u; k < nao; ++k) {
                    float pv[8], cv[4];
                    for (int i = 0; i < 8; ++i) pv[i] = e_buf[k * 128u + ty * 8u + uint(i)];
                    for (int j = 0; j < 4; ++j) cv[j] = e_c[k * 64u + tx * 4u + uint(j)];
                    for (int i = 0; i < 8; ++i) for (int j = 0; j < 4; ++j) X[i][j] += pv[i] * cv[j];
                }
                GroupMemoryBarrierWithGroupSync();
            }
        }
        // epilogue: stage 64 points x 64 functions at a time, then the terms per point and field
        for (uint half = 0u; half < 2u; ++half) {
            GroupMemoryBarrierWithGroupSync();
            if ((ty >> 3) == half) {
                for (int i = 0; i < 8; ++i) {
                    uint row = (ty & 7u) * 8u + uint(i);
                    for (int j = 0; j < 4; ++j) e_buf[row * 64u + tx * 4u + uint(j)] = X[i][j];
                }
            }
            GroupMemoryBarrierWithGroupSync();
            for (uint it = tid; it < 64u * a.num_fields; it += 256u) {
                uint row = it & 63u, f = it >> 6;
                e_out[(half * 64u + row) * 16u + f] += field_terms(f, f0, row);
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();
    uint fbits = field_bits(brick.w);
    uint slot0 = a.brick_slot[b];
    for (uint it = tid; it < 125u * a.num_fields; it += 256u) {
        uint p = it % 125u, f = it / 125u;
        if (((fbits >> (2u * f)) & 1u) == 0u) continue;
        uint slot = slot0 + field_slot(brick.w, f);
        a.corners[slot * 125u + p] = e_out[p * 16u + f];
    }
}

// Value and gradient at 64 vertices per group. e_buf holds 4 channels x 16 AOs x 64 points.
[shader("compute")]
[numthreads(16, 16, 1)]
void iso_eval(uint3 gid : SV_GroupID, uint3 gt : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint v0 = node_index(gid) * 64u;
    if (v0 >= a.count) return;                             // uniform per group
    uint tx = gt.x, ty = gt.y;
    if (tid == 0u) {
        uint done = 1u;
        for (uint k = 0u; k < 64u && v0 + k < a.count; ++k) done &= a.vdone[v0 + k] != 0u ? 1u : 0u;   // 1 converged, 2 predicted
        e_nlist = done;
    }
    GroupMemoryBarrierWithGroupSync();
    if (e_nlist != 0u) return;                             // uniform: every vertex here has converged
    if (tid < 64u) e_p[tid] = float4(a.vpos[min(v0 + tid, a.count - 1u)].xyz, 0.0);
    for (uint it = tid; it < 64u * 4u; it += 256u) e_out[it] = 0.0;
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        float3 lo = e_p[0].xyz, hi = e_p[0].xyz;
        for (uint k = 1u; k < 64u; ++k) { lo = min(lo, e_p[k].xyz); hi = max(hi, e_p[k].xyz); }
        e_lo = float4(lo, 0.0); e_hi = float4(hi, 0.0);
    }
    GroupMemoryBarrierWithGroupSync();
    float3 glo = e_lo.xyz, ghi = e_hi.xyz;

    for (uint f0 = 0u; f0 < a.num_fn; f0 += 64u) {
        float X[4][4][4];
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) for (int c = 0; c < 4; ++c) X[i][j][c] = 0.0;
        for (uint s0 = 0u; s0 < a.num_shells; s0 += E_LIST) {
            GroupMemoryBarrierWithGroupSync();
            list_shells(s0, glo, ghi, tid);
            for (;;) {
                next_chunk(tid);
                uint nao = e_ch_nao;
                if (nao == 0u) break;
                load_c(f0, tid);
                // per (shell, vertex): the radial part and its derivative once, then the shell's AOs
                uint nsh = e_ch_nsh;
                for (uint it = tid; it < nsh * 64u; it += 256u) {
                    uint s = it >> 6, p = it & 63u;
                    IsoShell sh = a.shells[min(e_ch_shell[s], a.num_shells - 1u)];
                    float3 d = e_p[p].xyz - sh.c.xyz;
                    float r2 = dot(d, d);
                    float R = 0.0, R1 = 0.0;
                    for (uint q = 0u; q < sh.info.z; ++q) {
                        float2 pr = a.prims[sh.info.w + q];
                        float t = pr.y * exp(-pr.x * r2);
                        R += t;
                        R1 -= pr.x * t;
                    }
                    uint k0 = e_ch_first[s];
                    for (uint n = 0u; n < sh.info.y; ++n) {
                        uint ao = sh.info2.x + n;
                        uint ijk = a.ao_ijk[ao];
                        uint i = ijk & 0xFu, j = (ijk >> 4) & 0xFu, kk = (ijk >> 8) & 0xFu;
                        float xi = ipow(d.x, i), yj = ipow(d.y, j), zk = ipow(d.z, kk);
                        float M = xi * yj * zk;
                        float3 Mg = float3(i > 0u ? float(i) * ipow(d.x, i - 1u) * yj * zk : 0.0,
                                           j > 0u ? float(j) * xi * ipow(d.y, j - 1u) * zk : 0.0,
                                           kk > 0u ? float(kk) * xi * yj * ipow(d.z, kk - 1u) : 0.0);
                        float nrm = a.ao_nrm[ao];
                        float3 g = nrm * (Mg * R + M * 2.0 * R1 * d);
                        uint at = (k0 + n) * 64u + p;
                        e_buf[0u * 1024u + at] = nrm * M * R;
                        e_buf[1u * 1024u + at] = g.x;
                        e_buf[2u * 1024u + at] = g.y;
                        e_buf[3u * 1024u + at] = g.z;
                    }
                }
                GroupMemoryBarrierWithGroupSync();
                for (uint k = 0u; k < nao; ++k) {
                    float cv[4];
                    for (int j = 0; j < 4; ++j) cv[j] = e_c[k * 64u + tx * 4u + uint(j)];
                    for (int i = 0; i < 4; ++i) {
                        uint p = ty * 4u + uint(i);
                        for (int c = 0; c < 4; ++c) {
                            float pv = e_buf[uint(c) * 1024u + k * 64u + p];
                            for (int j = 0; j < 4; ++j) X[i][j][c] += pv * cv[j];
                        }
                    }
                }
                GroupMemoryBarrierWithGroupSync();
            }
        }
        // epilogue: 16 points at a time, staged as [point][channel][64 functions]
        for (uint quarter = 0u; quarter < 4u; ++quarter) {
            GroupMemoryBarrierWithGroupSync();
            if ((ty >> 2) == quarter) {
                for (int i = 0; i < 4; ++i) {
                    uint row = (ty & 3u) * 4u + uint(i);
                    for (int c = 0; c < 4; ++c)
                        for (int j = 0; j < 4; ++j) e_buf[(row * 4u + uint(c)) * 64u + tx * 4u + uint(j)] = X[i][j][c];
                }
            }
            GroupMemoryBarrierWithGroupSync();
            if (tid < 16u) {
                uint row = tid, p = quarter * 16u + row;
                uint f = a.vinfo[min(v0 + p, a.count - 1u)].z;
                IsoField fd = a.fields[f];
                float4 acc = float4(0.0);
                for (uint t = 0u; t < fd.term_count; ++t) {
                    IsoTerm tm = a.terms[fd.term_off + t];
                    if (tm.a < f0 || tm.a >= f0 + 64u) continue;
                    uint ia = tm.a - f0;
                    float xa = e_buf[(row * 4u + 0u) * 64u + ia];
                    float3 ga = float3(e_buf[(row * 4u + 1u) * 64u + ia], e_buf[(row * 4u + 2u) * 64u + ia], e_buf[(row * 4u + 3u) * 64u + ia]);
                    if (tm.b == INVALID) { acc += tm.w * float4(ga, xa); continue; }
                    uint ib = tm.b - f0;
                    float xb = e_buf[(row * 4u + 0u) * 64u + ib];
                    float3 gb = float3(e_buf[(row * 4u + 1u) * 64u + ib], e_buf[(row * 4u + 2u) * 64u + ib], e_buf[(row * 4u + 3u) * 64u + ib]);
                    acc += tm.w * float4(xa * gb + xb * ga, xa * xb);
                }
                e_out[p * 4u + 0u] += acc.x; e_out[p * 4u + 1u] += acc.y;
                e_out[p * 4u + 2u] += acc.z; e_out[p * 4u + 3u] += acc.w;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid < 64u && v0 + tid < a.count) {
        a.vval[v0 + tid] = float4(e_out[tid * 4u + 0u], e_out[tid * 4u + 1u], e_out[tid * 4u + 2u], e_out[tid * 4u + 3u]);
    }
}

// Corners of sets of single-function fields (orbitals): one thread per corner point, every active
// field's function accumulated in registers (the shape of the dense orbital kernel).
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_corners_orbital(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count || tid >= 125u) return;
    uint4 brick = a.bricks[b];
    uint fbits = field_bits(brick.w);
    float3 p = corner_point(brick, tid);
    uint col[16];
    float wgt[16];
    float val[16];
    for (uint f = 0u; f < 16u; ++f) {
        val[f] = 0.0;
        col[f] = 0u;
        wgt[f] = 0.0;
        if (f < a.num_fields && ((fbits >> (2u * f)) & 1u) != 0u) {
            IsoTerm tm = a.terms[a.fields[f].term_off];
            col[f] = tm.a;
            wgt[f] = tm.w;
        }
    }
    for (uint s = 0u; s < a.num_shells; ++s) {
        IsoShell sh = a.shells[s];
        float3 d = p - sh.c.xyz;
        float r2 = dot(d, d);
        if (r2 > sh.c.w) continue;
        float R = 0.0;
        for (uint q = 0u; q < sh.info.z; ++q) {
            float2 pr = a.prims[sh.info.w + q];
            R += pr.y * exp(-pr.x * r2);
        }
        for (uint n = 0u; n < sh.info.y; ++n) {
            uint ao = sh.info2.x + n;
            uint ijk = a.ao_ijk[ao];
            float phi = a.ao_nrm[ao] * ipow(d.x, ijk & 0xFu) * ipow(d.y, (ijk >> 4) & 0xFu) * ipow(d.z, (ijk >> 8) & 0xFu) * R;
            uint row = ao * a.fn_stride;
            for (uint f = 0u; f < 16u; ++f) {
                if (wgt[f] != 0.0) val[f] += a.C[row + col[f]] * phi;
            }
        }
    }
    uint slot = a.brick_slot[b];
    for (uint f = 0u; f < 16u; ++f) {
        if (wgt[f] == 0.0) continue;
        a.corners[slot * 125u + tid] = wgt[f] * val[f];
        slot++;
    }
}

// ------------------------------------------------------------------------------------------
// Corners, register form: one thread per lattice point of the brick (125 of 128), the point's
// values of F function columns accumulated in registers. The shells that reach the brick are
// listed (128 at a time), packed into chunks of at most K_CK AOs, and each chunk's rows of C are
// staged in groupshared memory, read by every thread at once (broadcast). The terms per field
// are applied from a staged copy of the accumulators. F: 1 (one orbital), 8, 16, 32 (all the
// functions of the set) or 64 (tiles of 64 functions, as the pair terms require).
// ------------------------------------------------------------------------------------------

#define K_GROUP 128                      // corners: threads (points) per group
#define K_LIST  256                      // the most shells listed at once (the group size)
#define K_CK    32                       // AOs per chunk (whole shells)
groupshared uint  k_sh[K_LIST];          // listed shells
groupshared uint  k_nc[K_LIST];          // their AO counts
groupshared uint  k_row[K_LIST];         // their first row within their chunk
groupshared uint  k_chunk[K_LIST + 1];   // chunk c: listed shells [k_chunk[c], k_chunk[c + 1])
groupshared uint  k_nchunk;
groupshared uint  k_nrows[K_LIST];       // rows per chunk
groupshared uint  k_row_ao[K_CK];        // the current chunk: global AO of each row
groupshared uint  k_pts[K_GROUP];        // corners: the points evaluated, in thread order

// The chunk rows of C (times the AO normalization) [row][F], and the epilogue's staging [64 points][F]:
// 64 F floats, sized per F so that a narrow kernel keeps its occupancy. Read as float4 where F allows.
interface IKBuf {
    static float  ld(uint i);
    static float4 ld4(uint i4);
    static void   st(uint i, float v);
}
#define K_BUF(NAME, FLOATS)                                                         \
    groupshared float4 NAME##_data[(FLOATS) / 4];                                   \
    struct NAME : IKBuf {                                                           \
        static float  ld(uint i)          { return NAME##_data[i >> 2][i & 3u]; }   \
        static float4 ld4(uint i4)        { return NAME##_data[i4]; }               \
        static void   st(uint i, float v) { NAME##_data[i >> 2][i & 3u] = v; }      \
    };
K_BUF(KBuf1, 64)
K_BUF(KBuf8, 512)
K_BUF(KBuf16, 1024)
K_BUF(KBuf32, 2048)
K_BUF(KBuf64, 4096)

// The Cartesian monomials of degree l at d, in md_gto's order (gto_*_lmn in md_gto.c)
void cart_monomials(uint l, float3 d, out float m[15]) {
    [ForceUnroll] for (int n = 0; n < 15; ++n) m[n] = 0.0;
    float x = d.x, y = d.y, z = d.z;
    if (l == 0u) {
        m[0] = 1.0;
    } else if (l == 1u) {
        m[0] = x; m[1] = y; m[2] = z;
    } else if (l == 2u) {
        m[0] = x * x; m[1] = x * y; m[2] = x * z; m[3] = y * y; m[4] = y * z; m[5] = z * z;
    } else if (l == 3u) {
        float xx = x * x, yy = y * y, zz = z * z;
        m[0] = xx * x; m[1] = xx * y; m[2] = xx * z; m[3] = x * yy; m[4] = x * y * z;
        m[5] = x * zz; m[6] = yy * y; m[7] = yy * z; m[8] = y * zz; m[9] = zz * z;
    } else {
        float xx = x * x, yy = y * y, zz = z * z;
        m[0] = xx * xx; m[1] = xx * x * y; m[2] = xx * x * z; m[3] = xx * yy; m[4] = xx * y * z;
        m[5] = xx * zz; m[6] = x * yy * y; m[7] = x * yy * z; m[8] = x * y * zz; m[9] = x * zz * z;
        m[10] = yy * yy; m[11] = yy * y * z; m[12] = yy * zz; m[13] = y * zz * z; m[14] = zz * zz;
    }
}

// Lists the shells of block [s0, s0 + gs) that reach [lo, hi] and packs them into chunks (all gs
// threads of the group).
void k_list(uint s0, float3 lo, float3 hi, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    // indices clamped into range: some implementations (llvmpipe) load for inactive lanes too
    uint s = min(s0 + tid, a.num_shells - 1u);
    IsoShell sh = a.shells[s];
    bool reach = s0 + tid < a.num_shells && box_dist2(sh.c.xyz, lo, hi) <= sh.c.w;
    uint nc = sh.info.y;
    uint total;
    uint at = iso_group_prefix(reach, tid, gs, total);
    if (reach) { k_sh[at] = s; k_nc[at] = nc; }
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {        // greedy packing over groupshared memory: at most gs steps
        uint c = 0u, rows = 0u;
        k_chunk[0] = 0u;
        for (uint i = 0u; i < total; ++i) {
            uint n = k_nc[i];
            if (rows > 0u && rows + n > K_CK) { k_nrows[c] = rows; c++; k_chunk[c] = i; rows = 0u; }
            k_row[i] = rows;
            rows += n;
        }
        if (total > 0u) { k_nrows[c] = rows; c++; }
        k_chunk[c] = total;
        k_nchunk = c;
    }
    GroupMemoryBarrierWithGroupSync();
}

// Stages chunk c's rows of C times the AO normalization, columns [f0, f0 + F) (all threads; ends
// with a barrier).
void k_load_chunk<let F : int, B : IKBuf>(uint c, uint f0, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
    uint i = min(i0 + tid, i1 - 1u);       // clamped: see k_list
    uint ao0 = a.shells[min(k_sh[i], a.num_shells - 1u)].info2.x;
    if (tid < i1 - i0) {
        for (uint n = 0u; n < k_nc[i]; ++n) k_row_ao[k_row[i] + n] = ao0 + n;
    }
    GroupMemoryBarrierWithGroupSync();
    uint nrows = k_nrows[c];
    for (uint e = tid; e < nrows * uint(F); e += gs) {
        uint row = e / uint(F), col = e % uint(F);
        uint ao = min(k_row_ao[min(row, uint(K_CK) - 1u)], a.num_ao - 1u);
        float v = a.C[ao * a.fn_stride + min(f0 + col, a.fn_stride - 1u)] * a.ao_nrm[ao];
        B::st(e, f0 + col < a.num_fn ? v : 0.0);
    }
    GroupMemoryBarrierWithGroupSync();
}

void corners_reg<let F : int, B : IKBuf>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    uint4 brick = a.bricks[b];
    float h = a.origin_h.w;
    float3 blo = a.origin_h.xyz + float3(brick.xyz * 4u) * h;
    float3 bhi = blo + 4.0 * h;
    uint fbits = field_bits(brick.w);
    // the points this brick evaluates, first in thread order: all 125, or (mode 1, shared corners)
    // those not taken from a neighbour; the threads past them skip the work, whole warps of them free
    bool mine = tid < 125u;
    if (a.mode == 1u && mine) {
        uint op;
        mine = corner_owner(brick, tid, op) == INVALID;
    }
    uint npts;
    uint at = iso_group_prefix(mine, tid, K_GROUP, npts);
    if (mine) k_pts[at] = tid;
    if (tid == 0u) InterlockedAdd(a.counters[5], npts);
    GroupMemoryBarrierWithGroupSync();
    bool active = tid < npts;
    uint pt = active ? k_pts[tid] : 124u;
    float3 p = corner_point(brick, pt);

    float val[16];
    [ForceUnroll] for (int f = 0; f < 16; ++f) val[f] = 0.0;

    for (uint f0 = 0u; f0 < a.num_fn; f0 += uint(F)) {
        float acc[F];
        [ForceUnroll] for (int j = 0; j < F; ++j) acc[j] = 0.0;
        for (uint s0 = 0u; s0 < a.num_shells; s0 += K_GROUP) {
            GroupMemoryBarrierWithGroupSync();             // the previous block's lists consumed
            k_list(s0, blo, bhi, tid, K_GROUP);
            uint nchunk = k_nchunk;
            for (uint c = 0u; c < nchunk; ++c) {
                k_load_chunk<F, B>(c, f0, tid, K_GROUP);
                uint i1 = active ? k_chunk[c + 1] : k_chunk[c];
                for (uint i = k_chunk[c]; i < i1; ++i) {    // uniform among the active threads
                    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
                    float3 d = p - sh.c.xyz;
                    float r2 = dot(d, d);
                    float R = 0.0;
                    for (uint q = 0u; q < sh.info.z; ++q) {
                        float2 pr = a.prims[sh.info.w + q];
                        R += pr.y * exp(-pr.x * r2);
                    }
                    float m[15];
                    cart_monomials(sh.info.x, d, m);
                    uint nc = sh.info.y;
                    uint row0 = k_row[i];
                    if (F == 1) {
                        // one function: the shell's polynomial, then its radial part once
                        float P = 0.0;
                        [ForceUnroll] for (int n = 0; n < 15; ++n) {
                            if (uint(n) < nc) P += m[n] * B::ld(row0 + uint(n));
                        }
                        acc[0] += P * R;
                    } else {
                        [ForceUnroll] for (int n = 0; n < 15; ++n) {
                            if (uint(n) < nc) {                 // uniform
                                float phi = m[n] * R;
                                uint rb4 = (row0 + uint(n)) * uint(F / 4);
                                [ForceUnroll] for (int j4 = 0; j4 < F / 4; ++j4) {
                                    float4 cv = B::ld4(rb4 + uint(j4));
                                    acc[4 * j4 + 0] += phi * cv.x;
                                    acc[4 * j4 + 1] += phi * cv.y;
                                    acc[4 * j4 + 2] += phi * cv.z;
                                    acc[4 * j4 + 3] += phi * cv.w;
                                }
                            }
                        }
                    }
                }
                GroupMemoryBarrierWithGroupSync();         // the chunk's rows consumed
            }
        }
        if (F == 1) {
            // one orbital: its one linear term
            IsoTerm tm = a.terms[a.fields[0].term_off];
            val[0] += tm.w * acc[0];
        } else {
            // the terms per field, from the accumulators staged 64 points at a time
            for (uint half = 0u; half < 2u; ++half) {
                GroupMemoryBarrierWithGroupSync();
                if ((tid >> 6) == half) {
                    uint row = tid & 63u;
                    [ForceUnroll] for (int j = 0; j < F; ++j) B::st(row * uint(F) + uint(j), acc[j]);
                }
                GroupMemoryBarrierWithGroupSync();
                if ((tid >> 6) == half) {
                    uint row = tid & 63u;
                    [ForceUnroll] for (int f = 0; f < 16; ++f) {
                        if (uint(f) < a.num_fields && ((fbits >> (2u * uint(f))) & 1u) != 0u) {
                            IsoField fd = a.fields[f];
                            float v = 0.0;
                            for (uint t = 0u; t < fd.term_count; ++t) {
                                IsoTerm tm = a.terms[fd.term_off + t];
                                if (tm.a < f0 || tm.a >= f0 + uint(F)) continue;
                                float xa = B::ld(row * uint(F) + tm.a - f0);
                                if (tm.b == INVALID) { v += tm.w * xa; continue; }
                                v += tm.w * xa * B::ld(row * uint(F) + tm.b - f0);
                            }
                            val[f] += v;
                        }
                    }
                }
            }
        }
    }
    if (!active) return;
    uint slot = a.brick_slot[b];
    [ForceUnroll] for (int f = 0; f < 16; ++f) {
        if (uint(f) < a.num_fields && ((fbits >> (2u * uint(f))) & 1u) != 0u) {
            a.corners[slot * 125u + pt] = val[f];
            slot++;
        }
    }
}

// Shared corners: the points a brick took from its neighbours, copied from their slots
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_corners_fill(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count || tid >= 125u) return;
    uint4 brick = a.bricks[b];
    uint op;
    uint nb = corner_owner(brick, tid, op);
    if (nb == INVALID) return;
    uint nmask = a.bricks[nb].w;
    uint slot = a.brick_slot[b], nslot0 = a.brick_slot[nb];
    for (uint f = 0u; f < a.num_fields; ++f) {
        if (((field_bits(brick.w) >> (2u * f)) & 1u) == 0u) continue;
        a.corners[slot * 125u + tid] = a.corners[(nslot0 + field_slot(nmask, f)) * 125u + op];
        slot++;
    }
}

// ------------------------------------------------------------------------------------------
// Corners from a point list: the lattice points to evaluate, brick by brick (all 125 of a brick with
// sheets, or with shared corners those not taken from a neighbour), and groups of 256 consecutive
// points - 2 to 4 neighbouring bricks - evaluated together: every thread busy, and the setup of each
// shell shared by 256 points instead of one brick's
// ------------------------------------------------------------------------------------------

groupshared uint p_base;

// One group per brick: appends its points to plist (one atomic per brick, so the list keeps roughly
// the order of the bricks; a brick's points are contiguous). mode 1: shared corners. The count is
// counters[5].
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_points(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    uint4 brick = a.bricks[b];
    bool mine = tid < 125u && brick.w != 0u;
    if (a.mode == 1u && mine) {
        uint op;
        mine = corner_owner(brick, tid, op) == INVALID;
    }
    uint n;
    uint at = iso_group_prefix(mine, tid, 128u, n);
    if (tid == 0u) {
        uint base = 0u;
        if (n > 0u) InterlockedAdd(a.counters[5], n, base);
        p_base = base;
    }
    GroupMemoryBarrierWithGroupSync();
    if (mine) a.plist[p_base + at] = (b << 7) | tid;
}

#define L_GROUP 256
groupshared float4 l_lohi[2 * 64];       // per wave: the points' box

void corners_list<let F : int, B : IKBuf>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    uint total = a.count_in[0];
    uint e0 = node_index(gid) * uint(L_GROUP);
    if (e0 >= total) return;                               // uniform per group
    uint e = e0 + tid;
    bool active = e < total;
    uint ent = a.plist[min(e, total - 1u)];
    uint b = ent >> 7, pt = ent & 127u;
    uint4 brick = a.bricks[b];
    float3 p = corner_point(brick, pt);

    // the group's box: its points'
    {
        uint lanes = WaveGetLaneCount();
        uint w = iso_wave_index(tid, lanes);
        float3 lo = WaveActiveMin(active ? p : float3(3.0e38));
        float3 hi = WaveActiveMax(active ? p : float3(-3.0e38));
        if (WaveIsFirstLane()) { l_lohi[2u * w] = float4(lo, 0.0); l_lohi[2u * w + 1u] = float4(hi, 0.0); }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) {
            uint nw = iso_num_waves(L_GROUP, lanes);
            float3 glo = l_lohi[0].xyz, ghi = l_lohi[1].xyz;
            for (uint k = 1u; k < nw; ++k) { glo = min(glo, l_lohi[2u * k].xyz); ghi = max(ghi, l_lohi[2u * k + 1u].xyz); }
            l_lohi[0] = float4(glo, 0.0); l_lohi[1] = float4(ghi, 0.0);
        }
        GroupMemoryBarrierWithGroupSync();
    }
    float3 glo = l_lohi[0].xyz, ghi = l_lohi[1].xyz;
    uint fbits = field_bits(brick.w);

    float val[16];
    [ForceUnroll] for (int f = 0; f < 16; ++f) val[f] = 0.0;

    for (uint f0 = 0u; f0 < a.num_fn; f0 += uint(F)) {
        float acc[F];
        [ForceUnroll] for (int j = 0; j < F; ++j) acc[j] = 0.0;
        for (uint s0 = 0u; s0 < a.num_shells; s0 += L_GROUP) {
            GroupMemoryBarrierWithGroupSync();             // the previous block's lists consumed
            k_list(s0, glo, ghi, tid, L_GROUP);
            uint nchunk = k_nchunk;
            for (uint c = 0u; c < nchunk; ++c) {
                k_load_chunk<F, B>(c, f0, tid, L_GROUP);
                uint i1 = active ? k_chunk[c + 1] : k_chunk[c];
                for (uint i = k_chunk[c]; i < i1; ++i) {    // uniform among the active threads
                    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
                    float3 d = p - sh.c.xyz;
                    float r2 = dot(d, d);
                    // the group's box spans several bricks: a shell out of reach of the point is skipped,
                    // by whole waves where none of its points is in reach
                    if (r2 > sh.c.w) continue;
                    float R = 0.0;
                    for (uint q = 0u; q < sh.info.z; ++q) {
                        float2 pr = a.prims[sh.info.w + q];
                        R += pr.y * exp(-pr.x * r2);
                    }
                    float m[15];
                    cart_monomials(sh.info.x, d, m);
                    uint nc = sh.info.y;
                    uint row0 = k_row[i];
                    if (F == 1) {
                        float P = 0.0;
                        [ForceUnroll] for (int n = 0; n < 15; ++n) {
                            if (uint(n) < nc) P += m[n] * B::ld(row0 + uint(n));
                        }
                        acc[0] += P * R;
                    } else {
                        [ForceUnroll] for (int n = 0; n < 15; ++n) {
                            if (uint(n) < nc) {
                                float phi = m[n] * R;
                                uint rb4 = (row0 + uint(n)) * uint(F / 4);
                                [ForceUnroll] for (int j4 = 0; j4 < F / 4; ++j4) {
                                    float4 cv = B::ld4(rb4 + uint(j4));
                                    acc[4 * j4 + 0] += phi * cv.x;
                                    acc[4 * j4 + 1] += phi * cv.y;
                                    acc[4 * j4 + 2] += phi * cv.z;
                                    acc[4 * j4 + 3] += phi * cv.w;
                                }
                            }
                        }
                    }
                }
                GroupMemoryBarrierWithGroupSync();         // the chunk's rows consumed
            }
        }
        if (F == 1) {
            IsoTerm tm = a.terms[a.fields[0].term_off];
            val[0] += tm.w * acc[0];
        } else {
            // the terms of each point's fields (its brick's), from the accumulators staged 64 points at a time
            for (uint round = 0u; round < uint(L_GROUP) / 64u; ++round) {
                GroupMemoryBarrierWithGroupSync();
                if ((tid >> 6) == round) {
                    uint row = tid & 63u;
                    [ForceUnroll] for (int j = 0; j < F; ++j) B::st(row * uint(F) + uint(j), acc[j]);
                }
                GroupMemoryBarrierWithGroupSync();
                if ((tid >> 6) == round && active) {
                    uint row = tid & 63u;
                    [ForceUnroll] for (int f = 0; f < 16; ++f) {
                        if (uint(f) < a.num_fields && ((fbits >> (2u * uint(f))) & 1u) != 0u) {
                            IsoField fd = a.fields[f];
                            float v = 0.0;
                            for (uint t = 0u; t < fd.term_count; ++t) {
                                IsoTerm tm = a.terms[fd.term_off + t];
                                if (tm.a < f0 || tm.a >= f0 + uint(F)) continue;
                                float xa = B::ld(row * uint(F) + tm.a - f0);
                                if (tm.b == INVALID) { v += tm.w * xa; continue; }
                                v += tm.w * xa * B::ld(row * uint(F) + tm.b - f0);
                            }
                            val[f] += v;
                        }
                    }
                }
            }
        }
    }
    if (!active) return;
    uint slot = a.brick_slot[b];
    [ForceUnroll] for (int f = 0; f < 16; ++f) {
        if (uint(f) < a.num_fields && ((fbits >> (2u * uint(f))) & 1u) != 0u) {
            a.corners[slot * 125u + pt] = val[f];
            slot++;
        }
    }
}

[shader("compute")] [numthreads(L_GROUP, 1, 1)] void iso_corners_list_f1 (uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_list<1, KBuf1>(gid, tid); }
[shader("compute")] [numthreads(L_GROUP, 1, 1)] void iso_corners_list_f8 (uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_list<8, KBuf8>(gid, tid); }
[shader("compute")] [numthreads(L_GROUP, 1, 1)] void iso_corners_list_f16(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_list<16, KBuf16>(gid, tid); }
[shader("compute")] [numthreads(L_GROUP, 1, 1)] void iso_corners_list_f32(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_list<32, KBuf32>(gid, tid); }
[shader("compute")] [numthreads(L_GROUP, 1, 1)] void iso_corners_list_f64(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_list<64, KBuf64>(gid, tid); }

[shader("compute")] [numthreads(K_GROUP, 1, 1)] void iso_corners_f1 (uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_reg<1, KBuf1>(gid, tid); }
[shader("compute")] [numthreads(K_GROUP, 1, 1)] void iso_corners_f8 (uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_reg<8, KBuf8>(gid, tid); }
[shader("compute")] [numthreads(K_GROUP, 1, 1)] void iso_corners_f16(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_reg<16, KBuf16>(gid, tid); }
[shader("compute")] [numthreads(K_GROUP, 1, 1)] void iso_corners_f32(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_reg<32, KBuf32>(gid, tid); }
[shader("compute")] [numthreads(K_GROUP, 1, 1)] void iso_corners_f64(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { corners_reg<64, KBuf64>(gid, tid); }

// The derivatives of the Cartesian monomials of degree l at d along axis ax (0 x, 1 y, 2 z), in
// md_gto's order: d/dx x^i y^j z^k = i x^(i-1) y^j z^k
void cart_monomials_d(uint l, float3 d, uint ax, out float g[15]) {
    [ForceUnroll] for (int n = 0; n < 15; ++n) g[n] = 0.0;
    float x = d.x, y = d.y, z = d.z;
    if (l == 1u) {
        g[ax] = 1.0;
    } else if (l == 2u) {
        if (ax == 0u)      { g[0] = 2.0 * x; g[1] = y; g[2] = z; }
        else if (ax == 1u) { g[1] = x; g[3] = 2.0 * y; g[4] = z; }
        else               { g[2] = x; g[4] = y; g[5] = 2.0 * z; }
    } else if (l == 3u) {
        if (ax == 0u)      { g[0] = 3.0 * x * x; g[1] = 2.0 * x * y; g[2] = 2.0 * x * z; g[3] = y * y; g[4] = y * z; g[5] = z * z; }
        else if (ax == 1u) { g[1] = x * x; g[3] = 2.0 * x * y; g[4] = x * z; g[6] = 3.0 * y * y; g[7] = 2.0 * y * z; g[8] = z * z; }
        else               { g[2] = x * x; g[4] = x * y; g[5] = 2.0 * x * z; g[7] = y * y; g[8] = 2.0 * y * z; g[9] = 3.0 * z * z; }
    } else if (l == 4u) {
        float xx = x * x, yy = y * y, zz = z * z;
        if (ax == 0u) {
            g[0] = 4.0 * xx * x; g[1] = 3.0 * xx * y; g[2] = 3.0 * xx * z; g[3] = 2.0 * x * yy; g[4] = 2.0 * x * y * z;
            g[5] = 2.0 * x * zz; g[6] = yy * y; g[7] = yy * z; g[8] = y * zz; g[9] = zz * z;
        } else if (ax == 1u) {
            g[1] = xx * x; g[3] = 2.0 * xx * y; g[4] = xx * z; g[6] = 3.0 * x * yy; g[7] = 2.0 * x * y * z;
            g[8] = x * zz; g[10] = 4.0 * yy * y; g[11] = 3.0 * yy * z; g[12] = 2.0 * y * zz; g[13] = zz * z;
        } else {
            g[2] = xx * x; g[4] = xx * y; g[5] = 2.0 * xx * z; g[7] = x * yy; g[8] = 2.0 * x * y * z;
            g[9] = 3.0 * x * zz; g[11] = yy * y; g[12] = 2.0 * yy * z; g[13] = 3.0 * y * zz; g[14] = 4.0 * zz * z;
        }
    }
}

// Value and gradient at 64 vertices per group, register form: 256 threads, thread (channel, vertex)
// with channel = tid / 64 (a warp holds one channel: 0 the value, 1..3 d/dx, d/dy, d/dz) accumulating
// that channel of F function columns. The shells that reach the group's box are listed and staged
// as for the corners; the terms per field are applied by one thread per vertex from the channels
// staged 16 vertices at a time.
#define V_GROUP 256
groupshared float4 v_res[64];            // per vertex: gradient xyz, value w, over the tiles
groupshared float4 v_lohi[2];
groupshared uint   v_done;

void eval_reg<let F : int, B : IKBuf>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    uint v0 = node_index(gid) * 64u;
    if (v0 >= a.count) return;                             // uniform per group
    uint vl = tid & 63u, ch = tid >> 6;
    uint vi = min(v0 + vl, a.count - 1u);
    if (tid == 0u) {
        uint done = 1u;
        for (uint k = 0u; k < 64u && v0 + k < a.count; ++k) done &= a.vdone[v0 + k] != 0u ? 1u : 0u;   // 1 converged, 2 predicted
        v_done = done;
    }
    GroupMemoryBarrierWithGroupSync();
    if (v_done != 0u) return;                              // uniform: every vertex here has converged
    float3 p = a.vpos[vi].xyz;
    if (ch == 0u) { e_p[vl] = float4(p, 0.0); v_res[vl] = float4(0.0); }
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        float3 lo = e_p[0].xyz, hi = e_p[0].xyz;
        for (uint k = 1u; k < 64u; ++k) { lo = min(lo, e_p[k].xyz); hi = max(hi, e_p[k].xyz); }
        v_lohi[0] = float4(lo, 0.0); v_lohi[1] = float4(hi, 0.0);
    }
    GroupMemoryBarrierWithGroupSync();
    float3 glo = v_lohi[0].xyz, ghi = v_lohi[1].xyz;
    uint ax = ch > 0u ? ch - 1u : 0u;

    for (uint f0 = 0u; f0 < a.num_fn; f0 += uint(F)) {
        float acc[F];
        [ForceUnroll] for (int j = 0; j < F; ++j) acc[j] = 0.0;
        for (uint s0 = 0u; s0 < a.num_shells; s0 += V_GROUP) {
            GroupMemoryBarrierWithGroupSync();             // the previous block's lists consumed
            k_list(s0, glo, ghi, tid, V_GROUP);
            uint nchunk = k_nchunk;
            for (uint c = 0u; c < nchunk; ++c) {
                k_load_chunk<F, B>(c, f0, tid, V_GROUP);
                uint i1 = k_chunk[c + 1];
                for (uint i = k_chunk[c]; i < i1; ++i) {    // uniform
                    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
                    float3 d = p - sh.c.xyz;
                    float r2 = dot(d, d);
                    float R = 0.0, R1 = 0.0;
                    for (uint q = 0u; q < sh.info.z; ++q) {
                        float2 pr = a.prims[sh.info.w + q];
                        float t = pr.y * exp(-pr.x * r2);
                        R += t;
                        R1 -= pr.x * t;
                    }
                    float m[15], g[15];
                    cart_monomials(sh.info.x, d, m);
                    cart_monomials_d(sh.info.x, d, ax, g);
                    // channel 0: m R; channel 1 + ax: g R + m 2 d[ax] R1
                    float da = ax == 0u ? d.x : (ax == 1u ? d.y : d.z);
                    float A = ch == 0u ? R : 2.0 * da * R1;
                    float Bg = ch == 0u ? 0.0 : R;
                    uint nc = sh.info.y;
                    uint row0 = k_row[i];
                    [ForceUnroll] for (int n = 0; n < 15; ++n) {
                        if (uint(n) < nc) {                 // uniform
                            float phi = m[n] * A + g[n] * Bg;
                            uint rb4 = (row0 + uint(n)) * uint(F / 4);
                            [ForceUnroll] for (int j4 = 0; j4 < F / 4; ++j4) {
                                float4 cv = B::ld4(rb4 + uint(j4));
                                acc[4 * j4 + 0] += phi * cv.x;
                                acc[4 * j4 + 1] += phi * cv.y;
                                acc[4 * j4 + 2] += phi * cv.z;
                                acc[4 * j4 + 3] += phi * cv.w;
                            }
                        }
                    }
                }
                GroupMemoryBarrierWithGroupSync();         // the chunk's rows consumed
            }
        }
        // the terms of each vertex's field, from the channels staged 16 vertices at a time
        for (uint quarter = 0u; quarter < 4u; ++quarter) {
            GroupMemoryBarrierWithGroupSync();
            if ((vl >> 4) == quarter) {
                uint row = (vl & 15u) * 4u + ch;
                [ForceUnroll] for (int j = 0; j < F; ++j) B::st(row * uint(F) + uint(j), acc[j]);
            }
            GroupMemoryBarrierWithGroupSync();
            if (ch == 0u && (vl >> 4) == quarter) {
                uint rb = (vl & 15u) * 4u * uint(F);
                IsoField fd = a.fields[min(a.vinfo[vi].z, 15u)];
                float4 r = float4(0.0);
                for (uint t = 0u; t < fd.term_count; ++t) {
                    IsoTerm tm = a.terms[fd.term_off + t];
                    if (tm.a < f0 || tm.a >= f0 + uint(F)) continue;
                    uint ia = tm.a - f0;
                    float xa = B::ld(rb + ia);
                    float3 ga = float3(B::ld(rb + uint(F) + ia), B::ld(rb + 2u * uint(F) + ia), B::ld(rb + 3u * uint(F) + ia));
                    if (tm.b == INVALID) { r += tm.w * float4(ga, xa); continue; }
                    uint ib = tm.b - f0;
                    float xb = B::ld(rb + ib);
                    float3 gb = float3(B::ld(rb + uint(F) + ib), B::ld(rb + 2u * uint(F) + ib), B::ld(rb + 3u * uint(F) + ib));
                    r += tm.w * float4(xa * gb + xb * ga, xa * xb);
                }
                v_res[vl] += r;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();
    if (ch == 0u && v0 + vl < a.count) a.vval[v0 + vl] = v_res[vl];
}

[shader("compute")] [numthreads(V_GROUP, 1, 1)] void iso_eval_f8 (uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { eval_reg<8, KBuf8>(gid, tid); }
[shader("compute")] [numthreads(V_GROUP, 1, 1)] void iso_eval_f16(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { eval_reg<16, KBuf16>(gid, tid); }
[shader("compute")] [numthreads(V_GROUP, 1, 1)] void iso_eval_f32(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { eval_reg<32, KBuf32>(gid, tid); }
[shader("compute")] [numthreads(V_GROUP, 1, 1)] void iso_eval_f64(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { eval_reg<64, KBuf64>(gid, tid); }

// Sets of single-function fields (orbitals): one thread per vertex, its field's one column of C.
// The tiled GEMM above would compute every field's function at every vertex.
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_eval_orbital(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint v = tid.x;
    if (v >= a.count || a.vdone[v] != 0u) return;
    float3 p = a.vpos[v].xyz;
    IsoField fd = a.fields[a.vinfo[v].z];
    IsoTerm tm = a.terms[fd.term_off];
    float val = 0.0;
    float3 grad = float3(0.0);
    for (uint s = 0u; s < a.num_shells; ++s) {
        IsoShell sh = a.shells[s];
        float3 d = p - sh.c.xyz;
        float r2 = dot(d, d);
        if (r2 > sh.c.w) continue;
        float R = 0.0, R1 = 0.0;
        for (uint q = 0u; q < sh.info.z; ++q) {
            float2 pr = a.prims[sh.info.w + q];
            float t = pr.y * exp(-pr.x * r2);
            R += t;
            R1 -= pr.x * t;
        }
        for (uint n = 0u; n < sh.info.y; ++n) {
            uint ao = sh.info2.x + n;
            float c = a.C[ao * a.fn_stride + tm.a];
            if (c == 0.0) continue;
            uint ijk = a.ao_ijk[ao];
            uint i = ijk & 0xFu, j = (ijk >> 4) & 0xFu, kk = (ijk >> 8) & 0xFu;
            float xi = ipow(d.x, i), yj = ipow(d.y, j), zk = ipow(d.z, kk);
            float M = xi * yj * zk;
            float3 Mg = float3(i > 0u ? float(i) * ipow(d.x, i - 1u) * yj * zk : 0.0,
                               j > 0u ? float(j) * xi * ipow(d.y, j - 1u) * zk : 0.0,
                               kk > 0u ? float(kk) * xi * yj * ipow(d.z, kk - 1u) : 0.0);
            float cn = c * a.ao_nrm[ao];
            val  += cn * M * R;
            grad += cn * (Mg * R + M * 2.0 * R1 * d);
        }
    }
    a.vval[v] = tm.w * float4(grad, val);
}

// ------------------------------------------------------------------------------------------
// Surface Nets
// ------------------------------------------------------------------------------------------

static const int3 CUBE_CORNER[8] = { int3(0,0,0), int3(1,0,0), int3(0,1,0), int3(1,1,0),
                                     int3(0,0,1), int3(1,0,1), int3(0,1,1), int3(1,1,1) };
static const uint2 CUBE_EDGE[12] = { uint2(0,1), uint2(2,3), uint2(4,5), uint2(6,7),
                                     uint2(0,2), uint2(1,3), uint2(4,6), uint2(5,7),
                                     uint2(0,4), uint2(1,5), uint2(2,6), uint2(3,7) };

float sheet_sign(uint s) { return s == 0u ? 1.0 : -1.0; }

groupshared uint c_count[64];
groupshared uint c_base;

// One group per brick, one thread per cell. mode 0 counts the vertices, mode 1 writes them, a
// brick's vertices contiguous (one atomic per brick).
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_cells(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;
    uint4 brick = a.bricks[b];
    uint cell = tid;
    int3 lc = int3(cell & 3u, (cell >> 2) & 3u, cell >> 4);
    uint slot0 = a.brick_slot[b];

    uint mine = 0u;
    for (uint f = 0u; f < a.num_fields; ++f) {
        uint fb = (brick.w >> (2u * f)) & 3u;
        if (fb == 0u) continue;
        uint slot = slot0 + field_slot(brick.w, f);
        float iso = a.fields[f].iso;
        for (uint s = 0u; s < 2u; ++s) {
            if ((fb & (1u << s)) == 0u) continue;
            float sg = sheet_sign(s);
            uint mask = 0u;
            for (uint k = 0u; k < 8u; ++k) if (sg * corner_value(slot, lc + CUBE_CORNER[k]) - iso > 0.0) mask |= 1u << k;
            if (mask != 0u && mask != 0xFFu) mine++;
        }
    }
    c_count[tid] = mine;
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {
        uint total = 0u;
        for (uint k = 0u; k < 64u; ++k) { uint c = c_count[k]; c_count[k] = total; total += c; }
        uint base = 0u;
        if (total > 0u) InterlockedAdd(a.counters[0], total, base);
        c_base = base;
    }
    GroupMemoryBarrierWithGroupSync();
    if (a.mode == 0u) return;

    uint vi = c_base + c_count[tid];
    for (uint f = 0u; f < a.num_fields; ++f) {
        uint fb = (brick.w >> (2u * f)) & 3u;
        if (fb == 0u) continue;
        uint slot = slot0 + field_slot(brick.w, f);
        float iso = a.fields[f].iso;
        float cv[8];
        for (uint k = 0u; k < 8u; ++k) cv[k] = corner_value(slot, lc + CUBE_CORNER[k]);
        for (uint s = 0u; s < 2u; ++s) {
            uint out_v = INVALID;
            if ((fb & (1u << s)) != 0u) {
                float sg = sheet_sign(s);
                float fv[8];
                uint mask = 0u;
                for (uint k = 0u; k < 8u; ++k) { fv[k] = sg * cv[k] - iso; if (fv[k] > 0.0) mask |= 1u << k; }
                if (mask != 0u && mask != 0xFFu) {
                    float3 acc = float3(0.0);
                    float cnt = 0.0;
                    for (uint e = 0u; e < 12u; ++e) {
                        uint e0 = CUBE_EDGE[e].x, e1 = CUBE_EDGE[e].y;
                        if (((mask >> e0) & 1u) == ((mask >> e1) & 1u)) continue;
                        float t = fv[e0] / (fv[e0] - fv[e1]);
                        acc += lerp(float3(CUBE_CORNER[e0]), float3(CUBE_CORNER[e1]), t);
                        cnt += 1.0;
                    }
                    float3 t = acc / cnt;
                    float3 gcell = float3(brick.xyz * 4u + uint3(lc));
                    float3 p = a.origin_h.xyz + (gcell + t) * a.origin_h.w;
                    // the trilinear interpolant's gradient there: the normal the relaxation keeps
                    // the vertices tangential to before they are projected
                    float3 g = float3(0.0);
                    for (uint k = 0u; k < 8u; ++k) {
                        float3 c = float3(CUBE_CORNER[k]);
                        float3 w = lerp(1.0 - t, t, c);         // per axis weight of corner k
                        float3 sgn = 2.0 * c - 1.0;
                        g += fv[k] * sgn * float3(w.y * w.z, w.x * w.z, w.x * w.y);
                    }
                    float gl = length(g);
                    if (vi < a.cap) {
                        a.vpos[vi]  = float4(p, sg);
                        a.vnrm[vi]  = float4(gl > 0.0 ? -g / gl : float3(0.0), 0.0);
                        a.vinfo[vi] = uint4(b, cell, f, s);
                        out_v = vi;
                    } else {
                        InterlockedAdd(a.counters[3], 1u);
                    }
                    vi++;
                }
            }
            a.cell_vert[slot * 128u + s * 64u + cell] = out_v;
        }
    }
}

// One quad per crossing lattice edge, from the four cells around it (thread per brick cell).
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_quads(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;
    uint4 brick = a.bricks[b];
    uint cell = tid;
    int3 lc = int3(cell & 3u, (cell >> 2) & 3u, cell >> 4);
    int3 gc = int3(brick.xyz * 4u) + lc;
    uint slot0 = a.brick_slot[b];
    for (uint f = 0u; f < a.num_fields; ++f) {
        uint fb = (brick.w >> (2u * f)) & 3u;
        if (fb == 0u) continue;
        uint slot = slot0 + field_slot(brick.w, f);
        float iso = a.fields[f].iso;
        for (uint s = 0u; s < 2u; ++s) {
            uint v0 = a.cell_vert[slot * 128u + s * 64u + cell];
            if (v0 == INVALID) continue;
            float sg = sheet_sign(s);
            float f0 = sg * corner_value(slot, lc) - iso;
            for (uint ax = 0u; ax < 3u; ++ax) {
                int3 o = int3(ax == 0u ? 1 : 0, ax == 1u ? 1 : 0, ax == 2u ? 1 : 0);
                float f1 = sg * corner_value(slot, lc + o) - iso;
                if ((f0 > 0.0) == (f1 > 0.0)) continue;
                uint u = (ax + 1u) % 3u, w = (ax + 2u) % 3u;
                int3 eu = int3(u == 0u ? 1 : 0, u == 1u ? 1 : 0, u == 2u ? 1 : 0);
                int3 ew = int3(w == 0u ? 1 : 0, w == 1u ? 1 : 0, w == 2u ? 1 : 0);
                uint v1 = cell_vertex(gc - eu, f, s);
                uint v2 = cell_vertex(gc - eu - ew, f, s);
                uint v3 = cell_vertex(gc - ew, f, s);
                if (v1 == INVALID || v2 == INVALID || v3 == INVALID) { InterlockedAdd(a.counters[2], 1u); continue; }
                uint4 q = (f0 <= 0.0) ? uint4(v0, v3, v2, v1) : uint4(v0, v1, v2, v3);
                uint qi;
                InterlockedAdd(a.counters[1], 1u, qi);
                if (qi < a.cap) a.quads[qi] = q;
                else InterlockedAdd(a.counters[3], 1u);
            }
        }
    }
}

// ------------------------------------------------------------------------------------------
// Vertices: Newton steps from iso_eval's values, relaxation
// ------------------------------------------------------------------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_newton(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint v = tid.x;
    if (v >= a.count) return;
    bool normals = a.mode == 1u;
    if (!normals && a.vdone[v] != 0u) return;
    float4 pv = a.vpos[v];
    float4 gv = a.vval[v];
    float sg = pv.w;
    float iso = a.fields[a.vinfo[v].z].iso;
    float f = sg * gv.w - iso;
    float3 g = sg * gv.xyz;
    float g2 = dot(g, g);
    if (g2 <= 1.0e-30) {
        if (normals) a.vnrm[v] = float4(0.0, 0.0, 0.0, 1.0e30);
        return;
    }
    float gl = sqrt(g2);
    if (normals) {
        a.vnrm[v] = float4(-g / gl, abs(f) / gl);    // f > 0 inside: outward is -grad f
        return;
    }
    float s = abs(f) / gl;
    if (s < a.tol) { a.vdone[v] = 1u; return; }      // converged: stays put, its values stay valid
    float step = f / g2;
    float len = abs(step) * gl;
    if (len > a.step_max) step *= a.step_max / len;
    a.vpos[v] = float4(pv.xyz - step * g, sg);
    if (a.mode == 2u && len <= a.step_max) {
        // EXPERIMENT: the error after this step predicted instead of checked by another evaluation:
        // about |f''| / (2 |f'|) s^2 along the gradient, with |f''| / |f'| ~ |f'| / |f| for the
        // exponentially decaying fields here. Predicted well below the tolerance: done (its normal and
        // values then those of the evaluation before the step)
        float K = gl / max(abs(sg * gv.w), 1.0e-30);
        if (0.5 * K * s * s < 0.1 * a.tol) a.vdone[v] = 2u;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_relax(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint v = tid.x;
    if (v >= a.count) return;
    uint4 vi = a.vinfo[v];
    uint b = vi.x, f = vi.z, s = vi.w;
    uint4 brick = a.bricks[b];
    uint slot = a.brick_slot[b] + field_slot(brick.w, f);
    int3 lc = int3(vi.y & 3u, (vi.y >> 2) & 3u, vi.y >> 4);
    int3 gc = int3(brick.xyz * 4u) + lc;
    float sg = sheet_sign(s);
    float iso = a.fields[f].iso;
    float4 pv = a.vpos[v];
    float3 acc = float3(0.0);
    float cnt = 0.0;
    // face neighbours that share a mesh edge: the shared face has a crossing lattice edge
    for (uint d = 0u; d < 6u; ++d) {
        uint ax = d >> 1;
        int side = (d & 1u) != 0u ? 1 : 0;
        int3 o = int3(ax == 0u ? 1 : 0, ax == 1u ? 1 : 0, ax == 2u ? 1 : 0);
        int3 u = int3(o.y, o.z, o.x), w = int3(o.z, o.x, o.y);
        int3 base = lc + o * side;
        float f00 = sg * corner_value(slot, base) - iso;
        float f10 = sg * corner_value(slot, base + u) - iso;
        float f01 = sg * corner_value(slot, base + w) - iso;
        float f11 = sg * corner_value(slot, base + u + w) - iso;
        bool cross = ((f00 > 0.0) != (f10 > 0.0)) || ((f00 > 0.0) != (f01 > 0.0)) ||
                     ((f11 > 0.0) != (f10 > 0.0)) || ((f11 > 0.0) != (f01 > 0.0));
        if (!cross) continue;
        uint vn = cell_vertex(gc + o * (side * 2 - 1), f, s);
        float3 pn = a.vpos[min(vn, a.count - 1u)].xyz;      // clamped: see cell_vertex
        if (vn == INVALID || vn >= a.count) continue;
        acc += pn;
        cnt += 1.0;
    }
    float3 p = pv.xyz;
    if (cnt > 0.0) {
        float3 n = a.vnrm[v].xyz;
        float3 dlt = acc / cnt - p;
        p += a.lambda * (dlt - dot(dlt, n) * n);
    }
    a.vpos_out[v] = float4(p, pv.w);
}

// ------------------------------------------------------------------------------------------
// Output: vertices and triangles regrouped per (field, sheet) bucket
// ------------------------------------------------------------------------------------------

uint vertex_bucket(uint v) { uint4 vi = MD_ARGS.vinfo[v]; return vi.z * 2u + vi.w; }

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_bucket_count(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint i = tid.x;
    if (i < a.count) InterlockedAdd(a.bucket[vertex_bucket(i)], 1u);
    if (i < a.cap)   InterlockedAdd(a.bucket[32u + vertex_bucket(a.quads[i].x)], 1u);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void iso_bucket_scan(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint ov = 0u, oq = 0u;
    for (uint k = 0u; k < 32u; ++k) {
        a.bucket[64u + k] = ov;  a.bucket[128u + k] = ov;  ov += a.bucket[k];
        a.bucket[96u + k] = oq;  a.bucket[160u + k] = oq;  oq += a.bucket[32u + k];
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_bucket_vertices(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint v = tid.x;
    if (v >= a.count) return;
    uint nv;
    InterlockedAdd(a.bucket[128u + vertex_bucket(v)], 1u, nv);
    a.remap[v] = nv;
    a.out_pos[nv] = float4(a.vpos[v].xyz, 1.0);
    a.out_nrm[nv] = a.vnrm[v];
}

// Quads to triangles along the shorter diagonal, written into their bucket with bucketed indices.
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_bucket_triangles(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint i = tid.x;
    if (i >= a.cap) return;
    uint4 q = a.quads[i];
    uint nq;
    InterlockedAdd(a.bucket[160u + vertex_bucket(q.x)], 1u, nq);
    uint4 r = uint4(a.remap[q.x], a.remap[q.y], a.remap[q.z], a.remap[q.w]);
    float3 p0 = a.vpos[q.x].xyz, p1 = a.vpos[q.y].xyz, p2 = a.vpos[q.z].xyz, p3 = a.vpos[q.w].xyz;
    uint o = nq * 6u;
    if (dot(p0 - p2, p0 - p2) <= dot(p1 - p3, p1 - p3)) {
        a.out_idx[o + 0] = r.x; a.out_idx[o + 1] = r.y; a.out_idx[o + 2] = r.z;
        a.out_idx[o + 3] = r.x; a.out_idx[o + 4] = r.z; a.out_idx[o + 5] = r.w;
    } else {
        a.out_idx[o + 0] = r.x; a.out_idx[o + 1] = r.y; a.out_idx[o + 2] = r.w;
        a.out_idx[o + 3] = r.y; a.out_idx[o + 4] = r.z; a.out_idx[o + 5] = r.w;
    }
}

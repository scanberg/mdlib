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
//                   rounding step is bounded and added. iso_octree_small: up to 4 functions;
//                   iso_octree_dform: densities as phi^T D phi. The kept nodes' children are placed
//                   in Morton order (iso_scan, iso_octree_emit)
//   iso_brick_slots per brick its corner slots (one per field with a sheet bit), placed by iso_scan
//   iso_hash        brick coordinate -> brick index
//   iso_points      the lattice points to evaluate, brick by brick (shared corners: a point on a brick's
//                   high face is evaluated by the neighbour that owns it, iso_corners_fill copies it)
//   corners         the field values at the listed points, consecutive points per group:
//                   iso_corners_list_f1v4 (one function), iso_gcorners_f* (a register-tiled GEMM, AO
//                   values x C, the terms in the epilogue), iso_gcorners_w*/y256 (canonical sets of more
//                   than 64 functions), iso_dcorners (densities as phi^T D phi)
//   iso_flood_*     APPROXIMATE (densities with maxima at the atoms): seed and flood, so only the bricks a
//                   surface reaches from a seed are evaluated in full
//   iso_cells       Surface Nets vertex per cell, field and sheet; a brick's vertices are contiguous
//   iso_quads       one quad per crossing lattice edge
//   iso_relax       tangential Laplacian toward the face-connected neighbours
//   vertices        value and gradient (iso_geval_*, iso_deval, iso_eval_orbital) and one Newton step
//                   (iso_newton) per iteration; converged vertices are flagged and skipped until they
//                   move again. One orbital: the whole projection in one launch (iso_project_list_f1)
//   iso_bucket_*    vertices and triangles regrouped per field and sheet
//   iso_grid        indirect grids from device-side counts
//   iso_scan        exclusive prefix sums (the ordered lists)
// A composite set (LOL, ELF; "Composite fields" below): field 0 a function of the others, bounded in octree_finish
// from their ranges, its corners (iso_comp_corners) and vertex values (iso_comp_vertices) from theirs.
// Outside the extraction, the vertex evaluators serve the caller's points (md_gto_iso_gpu_evaluate_launch):
//   iso_points_load the points as vertices of one field; iso_points_store their values from vval

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
    uint      mode;           // iso_cells: 0 count, 1 emit. iso_newton: 0 step, 1 normals only. iso_points,
                              // iso_corners_fill: flags (see there)
    uint      count;          // host-known item count (vertices, quads)
    uint      grid_mul;       // iso_grid: groups (or threads) per counted item
    uint      grid_div;       // iso_grid: items per group
    float     step_max;       // Newton step clamp
    float     lambda;         // relaxation step
    float     tol;            // Newton: converged below this distance
    uint      canon;          // every field's terms in the canonical form: wts (the register-form epilogues)
    uint      _pad0;
    uint      dform;          // the density fields as phi^T D phi: dmat, fdmat (corners from the list, vertices)
    uint      num_dmat;
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
    float4*   vval;           // the vertex evaluators: gradient xyz, value w
    uint4*    vinfo;          // brick, cell, field, sheet
    uint4*    quads;
    uint*     counters;       // 0 vertices, 1 quads, 2 open edges, 3 overflow, 4 corner slots, 5 corner points evaluated,
                              // 6 flood: bricks of the next round, 7 dform list overflow, 8 bricks with a sheet bit
    uint*     bucket;         // [0, 32) vertex counts, [32, 64) quad counts, [64, 96) vertex offsets,
                              // [96, 128) quad offsets, [128, 160) vertex cursors, [160, 192) quad cursors
    uint*     remap;          // old vertex -> bucketed vertex
    uint*     vdone;          // per vertex: converged, no more evaluations until it moves
    float4*   out_pos;
    float4*   out_nrm;
    uint*     out_idx;
    float*    fn_tail;        // per function: bound of what the shells dropped from a box can add
    uint*     plist;          // corners: the lattice points to evaluate, brick << 7 | point
    float*    wts;            // canon: per field 3 x fn_stride weights: w a, w a^2, w a b (a even, b = a + 1)
    float*    dmat;           // dform: num_dmat symmetric matrices [num_ao][num_ao], with the AO normalization
    uint*     fdmat;          // dform: per field, its matrix
    float*    drs;            // dform: per matrix and AO mu, sum over nu of |D_mu nu| / |nrm_nu| (the octree's tails)
    float4*   dprop;          // dform: per matrix: (2 eps)^2 sum |D raw|, 4 eps, positive semidefinite (1), slack
    uint*     scan;           // ordered lists: per brick, its count (iso_points, iso_cells), then its offset (iso_scan)
    uint*     bstate;         // flood: per brick, the round it was evaluated in (1, 2, ...), FLOOD_NONE not (yet)
    uint      group_base;     // a launch in chunks: the first group of this one (node_index adds it)
    uint      grid_max;       // iso_grid: at most this many groups from group_base on (0: no chunks)
    // dform: the matrices in tiles (md_gto_gpu_dmat): dmat holds num_dmat pools of dm_pool floats, one pattern for
    // all (dm_tao, dm_map, tiles of 16 per edge, dm_nt per edge), without the AO normalizations (dform_at applies them)
    uint      _pad1;
    uint      dm_nt;
    uint      dm_pool;
    uint      _pad2;
    uint*     dm_tao;         // per AO: its tiled index
    uint*     dm_map;         // [dm_nt][dm_nt], I <= J: the tile pair's pool entry, INVALID where not kept
    float*    dm_rows;        // iso_dmat_rows: per matrix and AO, sum over nu of |D_mu nu|
    // iso_points_load, iso_points_store (md_gto_iso_gpu_evaluate_launch): the caller's points and values
    float*    pts;            // point i at pts[i * pts_stride], float xyz
    float*    pval;           // iso_points_store: the value at point i to pval[i]
    uint      pts_stride;     // floats from one point to the next
    uint      pts_field;      // iso_points_load: the field the points are evaluated for
    uint      _pad3;
    uint      _pad4;
    // A composite set (comp_kind != 0, see "Composite fields"): field 0 a function of fields 1 .. num_fields - 1
    float4*   vsub;           // the vertex values of field f >= 1 at vsub[(f - 1) * comp_stride + v] (gradient xyz, value w)
    uint      comp_kind;      // 0 none, 1 LOL, 2 ELF
    uint      comp_spins;     // 1: a closed shell's totals; 2: alpha and beta
    float     comp_rho_min;
    float     comp_d_min;
    uint      comp_stride;
    uint      vfield;         // the vertex evaluators: 0 each vertex's own field (vinfo.z), else field vfield - 1 for all
    uint      _pad5;
    uint      _pad6;
};
MD_SHADER_ARGS(IsoArgs);

// Element (i, j) of a D-form matrix, with the AO normalizations: D the matrix's pool of tiles, as IsoArgs describes it
float dform_at(float* D, uint* tao, uint* map, uint nt, float* nrm, uint i, uint j) {
    uint ti = tao[i], tj = tao[j];
    uint lo = min(ti, tj), hi = max(ti, tj);
    uint e = map[(lo >> 4) * nt + (hi >> 4)];
    return e == 0xFFFFFFFFu ? 0.0 : D[e * 256u + (lo & 15u) * 16u + (hi & 15u)] * nrm[i] * nrm[j];
}

static const uint  INVALID = 0xFFFFFFFFu;
static const float U = 1.0 / 4194304.0;     // 2^-22: four float ulps, generous for fast-math backends
static const float EXP_REL = 4.0e-6;        // relative error allowed for exp() on any backend

uint  node_index(uint3 gid) { IsoArgs a = MD_ARGS; return gid.y * 65535u + gid.x + a.group_base; }
uint  field_bits(uint mask) { return (mask | (mask >> 1)) & 0x55555555u; }     // bit 2f: field f active
uint  field_slot(uint mask, uint f) { return countbits(field_bits(mask) & ((1u << (2u * f)) - 1u)); }
// The field the vertex evaluators evaluate at vertex v: its own, or every vertex's (vfield, a composite's parts)
uint  vertex_field(uint v) { IsoArgs a = MD_ARGS; return a.vfield != 0u ? a.vfield - 1u : min(a.vinfo[v].z, 15u); }

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
    if (a.grid_max != 0u) groups = groups > a.group_base ? min(groups - a.group_base, a.grid_max) : 0u;   // a chunk
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

// One shell over a box: its radial range R (rounding and the exponentials' error covered), its radial
// value at the box centre, and what its AOs' ranges are made of
struct OShell {
    float2 R;
    float  Rc;
    float3 dlo, dhi, dc;
    uint   ao0, ncart;
};

OShell octree_shell_setup(uint s, float3 lo, float3 hi, float3 cen) {
    IsoArgs a = MD_ARGS;
    IsoShell sh = a.shells[s];
    uint nprim = sh.info.z, poff = sh.info.w;
    OShell o;
    o.dlo = lo - sh.c.xyz;
    o.dhi = hi - sh.c.xyz;
    float3 nearest = max(max(o.dlo, -o.dhi), float3(0.0));
    float3 farthest = max(abs(o.dlo), abs(o.dhi));
    // rho range widened for the rounding of the distances (the exponentials are monotone in rho)
    float rlo = dot(nearest, nearest) * (1.0 - 8.0 * U);
    float rhi = dot(farthest, farthest) * (1.0 + 8.0 * U) + 1.0e-30;
    o.dc = cen - sh.c.xyz;
    float rc = dot(o.dc, o.dc);
    float gu = 0.0, gv = 0.0, au = 0.0, Rc = 0.0;
    for (uint q = 0u; q < nprim; ++q) {
        float2 pr = a.prims[poff + q];
        float eu = exp(-pr.x * rlo), ev = exp(-pr.x * rhi), ec = exp(-pr.x * rc);
        gu += pr.y * eu;
        gv += pr.y * ev;
        au += abs(pr.y) * eu;
        Rc += pr.y * ec;
    }
    float2 R = float2(min(gu, gv), max(gu, gv));
    for (uint k = 0u; k < sh.info2.z; ++k) {
        float4 c = a.crit[sh.info2.y + k];
        if (c.y >= rlo && c.x <= rhi) R = float2(min(R.x, c.z), max(R.y, c.w));
    }
    float eR = au * (EXP_REL + float(nprim + 4u) * U) + 1.0e-35;
    o.R = R + float2(-eR, eR);
    o.Rc = Rc;
    o.ao0 = sh.info2.x;
    o.ncart = sh.info.y;
    return o;
}

// AO n of the shell over the box: its range (mid, rad), its value at the centre and the largest
// deviation from it; the powers taken directly (an array of them indexed by the AO's exponents is
// spilled to local memory)
void octree_ao(OShell o, uint n, out float mid, out float rad, out float pc, out float dev, out uint ao, bool raw) {
    IsoArgs a = MD_ARGS;
    ao = o.ao0 + n;
    uint ijk = a.ao_ijk[ao];
    uint i = ijk & 0xFu, j = (ijk >> 4) & 0xFu, kk = (ijk >> 8) & 0xFu;
    float nrm = raw ? 1.0 : a.ao_nrm[ao];               // raw: without the normalization (dform's D carries it)
    float2 M = iv_mul(iv_mul(iv_pow(float2(o.dlo.x, o.dhi.x), i), iv_pow(float2(o.dlo.y, o.dhi.y), j)), iv_pow(float2(o.dlo.z, o.dhi.z), kk));
    float2 P = iv_mul(M, o.R) * nrm;
    if (nrm < 0.0) P = P.yx;
    float pm = iv_mag(P);
    P += float2(-8.0 * U * pm, 8.0 * U * pm);
    pc = nrm * ipow(o.dc.x, i) * ipow(o.dc.y, j) * ipow(o.dc.z, kk) * o.Rc;
    mid = 0.5 * (P.x + P.y);
    rad = 0.5 * (P.y - P.x) * (1.0 + 4.0 * U) + U * abs(P.x + P.y);
    dev = max(max(P.y - pc, pc - P.x), 0.0) * (1.0 + 4.0 * U);
}

// AO data of one shell over the box, written at slot 'off'
void octree_shell(uint s, float3 lo, float3 hi, float3 cen, uint off) {
    OShell o = octree_shell_setup(s, lo, hi, cen);
    for (uint n = 0u; n < o.ncart; ++n) {
        float mid, rad, pc, dev;
        uint ao;
        octree_ao(o, n, mid, rad, pc, dev, ao, false);
        uint slot = off + n;
        o_mid[slot] = mid;
        o_rad[slot] = rad;
        o_pc[slot]  = pc;
        o_dev[slot] = dev;
        o_ao[slot]  = ao;
    }
}

// The node's box (the rounding of its lattice points covered) and centre
void octree_box(uint4 node, out float3 lo, out float3 hi, out float3 cen) {
    IsoArgs a = MD_ARGS;
    float h = a.origin_h.w;
    float edge = 4.0 * h * float(1u << (a.brick_level - a.level));
    lo = a.origin_h.xyz + float3(node.xyz) * edge;
    hi = lo + edge;
    // the lattice points inside are computed as origin + index * h in float: cover their rounding
    float slop = 1.0e-5 * edge + 4.0 * U * (abs(a.origin_h.x) + abs(a.origin_h.y) + abs(a.origin_h.z) + 3.0 * edge * float(node.x + node.y + node.z + 1u));
    lo -= slop;
    hi += slop;
    cen = 0.5 * (lo + hi);
}

// A function's range over the node from its five sums (centered form and natural form, intersected)
float2 octree_range(float sm[5], uint col) {
    IsoArgs a = MD_ARGS;
    float e = float(a.num_ao + 8u) * 2.0 * U * sm[4] + a.fn_tail[col] + 1.0e-35;
    float2 i1 = float2(sm[0] - sm[1] - e, sm[0] + sm[1] + e);
    float2 i2 = float2(sm[2] - sm[3] - e, sm[2] + sm[3] + e);
    float2 A = float2(max(i1.x, i2.x), min(i1.y, i2.y));
    if (A.x > A.y) A = float2(min(i1.x, i2.x), max(i1.y, i2.y));   // only by rounding; keep the hull
    return A;
}

// ------------------------------------------------------------------------------------------
// Composite fields (comp_kind != 0, md_gto_iso_field_create_lol, _elf): field 0 is a function of the others,
// which are density fields like any: per spin s (comp_spins; 1 a closed shell's totals, taken as two equal
// spins of half of them) fields 1 + s K: rho, tau (K = 2), and for ELF the components of grad rho (K = 5;
// products of the density's functions and tau's, the derivatives). The octree bounds field 0 from the others'
// ranges (octree_finish), its corners and vertex values come from theirs (iso_comp_corners, iso_comp_vertices).
//     LOL = t / (1 + t),  t = D0 / tau,   D0 = c_F sum_s rho_s^(5/3)
//     ELF = 1 / (1 + chi^2),  chi = (d_min + sum_s D_s) / D0,  D_s = max(0, tau_s - |grad rho_s|^2 / (8 rho_s))
// both 0 where the density is below rho_min (md_gto_iso.c's iso_lol, iso_elf_eval, in float).
// ------------------------------------------------------------------------------------------

#define COMP_MAX 10                                  // parts: two spins of ELF's five
static const float COMP_CF = 4.557799872345596;     // 3/10 (6 pi^2)^(2/3)
static const float COMP_MARGIN = 1.0e-5;            // the float evaluation's error in a bound of LOL or ELF (in [0, 1])

uint comp_parts(IsoArgs a) { return a.comp_kind == 2u ? 5u : 2u; }

groupshared float2 o_frng[16];                       // the octree: the parts' ranges over the node

// Field 0's range over the node from its parts' (o_frng): LOL non-decreasing in every rho_s and decreasing in
// every tau_s; ELF by intervals over its formula (rho counted independently where it appears twice)
float2 comp_range(IsoArgs a) {
    uint K = comp_parts(a);
    float rlo[2], rhi[2], tlo[2], thi[2], glo[2], ghi[2];
    [ForceUnroll] for (int s = 0; s < 2; ++s) { rlo[s] = 0.0; rhi[s] = 0.0; tlo[s] = 0.0; thi[s] = 0.0; glo[s] = 0.0; ghi[s] = 0.0; }
    uint ns = min(a.comp_spins, 2u);
    for (uint s = 0u; s < ns; ++s) {
        float2 r = o_frng[1u + s * K], t = o_frng[2u + s * K];
        rlo[s] = max(r.x, 0.0); rhi[s] = max(r.y, 0.0);
        tlo[s] = max(t.x, 0.0); thi[s] = max(t.y, 0.0);
        if (K == 5u) {
            float2 g = float2(0.0, 0.0);
            for (uint x = 0u; x < 3u; ++x) g += iv_sq(o_frng[3u + s * K + x]);
            glo[s] = g.x; ghi[s] = g.y;
        }
    }
    if (ns == 1u) {
        rlo[0] *= 0.5; rhi[0] *= 0.5; tlo[0] *= 0.5; thi[0] *= 0.5; glo[0] *= 0.25; ghi[0] *= 0.25;
        rlo[1] = rlo[0]; rhi[1] = rhi[0]; tlo[1] = tlo[0]; thi[1] = thi[0]; glo[1] = glo[0]; ghi[1] = ghi[0];
    }
    float rs_lo = rlo[0] + rlo[1], rs_hi = rhi[0] + rhi[1];
    float D0_lo = COMP_CF * (pow(rlo[0], 5.0 / 3.0) + pow(rlo[1], 5.0 / 3.0));
    float D0_hi = COMP_CF * (pow(rhi[0], 5.0 / 3.0) + pow(rhi[1], 5.0 / 3.0));
    if (!(rs_hi >= a.comp_rho_min) || !(D0_hi > 0.0)) return float2(0.0, COMP_MARGIN);
    float f_lo, f_hi;
    if (K == 5u) {
        float dmin = max(a.comp_d_min, 0.0);
        float D_lo = dmin, D_hi = dmin;
        [ForceUnroll] for (int s = 0; s < 2; ++s) {
            if (rlo[s] > 0.0) D_lo += max(0.0, tlo[s] - ghi[s] / (8.0 * rlo[s]));
            if (rhi[s] > 0.0) D_hi += max(0.0, thi[s] - glo[s] / (8.0 * rhi[s]));
        }
        float chi_lo = D_lo / D0_hi;
        f_hi = 1.0 / (1.0 + chi_lo * chi_lo);
        f_lo = 0.0;
        if (rs_lo >= a.comp_rho_min && D0_lo > 0.0) { float chi_hi = D_hi / D0_lo; f_lo = 1.0 / (1.0 + chi_hi * chi_hi); }
    } else {
        float ts_lo = tlo[0] + tlo[1], ts_hi = thi[0] + thi[1];
        f_hi = ts_lo > 0.0 ? D0_hi / (ts_lo + D0_hi) : 1.0;            // t / (1 + t) = D0 / (tau + D0)
        f_lo = 0.0;
        if (rs_lo >= a.comp_rho_min) f_lo = ts_hi > 0.0 ? D0_lo / (ts_hi + D0_lo) : 1.0;
    }
    return float2(f_lo - COMP_MARGIN, f_hi + COMP_MARGIN);
}

// Field 0's value and gradient (xyz, w the value) at a point from its parts' (p[i]: field 1 + i, gradient xyz
// and value w; the gradients of ELF's grad rho components are the rows of the density's Hessian)
float4 comp_eval(IsoArgs a, float4 p[COMP_MAX]) {
    uint K = comp_parts(a);
    uint ns = min(a.comp_spins, 2u);
    float rho[2], tau[2];
    float3 gr[2], gt[2], h0[2], h1[2], h2[2];       // ELF: grad rho, grad tau, the Hessian's rows
    [ForceUnroll] for (int s = 0; s < 2; ++s) {
        uint b = uint(s) * K;
        bool on = uint(s) < ns;
        rho[s] = on ? p[b].w : 0.0;      gr[s] = on ? p[b].xyz : float3(0.0);
        tau[s] = on ? p[b + 1u].w : 0.0; gt[s] = on ? p[b + 1u].xyz : float3(0.0);
        if (K == 5u) {
            gr[s] = on ? float3(p[b + 2u].w, p[b + 3u].w, p[b + 4u].w) : float3(0.0);
            h0[s] = on ? p[b + 2u].xyz : float3(0.0);
            h1[s] = on ? p[b + 3u].xyz : float3(0.0);
            h2[s] = on ? p[b + 4u].xyz : float3(0.0);
        } else {
            h0[s] = float3(0.0); h1[s] = float3(0.0); h2[s] = float3(0.0);
        }
    }
    if (ns == 1u) {
        rho[0] *= 0.5; tau[0] *= 0.5; gr[0] *= 0.5; gt[0] *= 0.5; h0[0] *= 0.5; h1[0] *= 0.5; h2[0] *= 0.5;
        rho[1] = rho[0]; tau[1] = tau[0]; gr[1] = gr[0]; gt[1] = gt[0]; h0[1] = h0[0]; h1[1] = h1[0]; h2[1] = h2[0];
    }
    float rho_sum = max(rho[0], 0.0) + max(rho[1], 0.0);
    float D0 = 0.0;
    float3 gD0 = float3(0.0);
    [ForceUnroll] for (int s = 0; s < 2; ++s) {
        if (!(rho[s] > 0.0)) continue;
        D0 += COMP_CF * pow(rho[s], 5.0 / 3.0);
        gD0 += 5.0 / 3.0 * COMP_CF * pow(rho[s], 2.0 / 3.0) * gr[s];
    }
    if (!(rho_sum >= a.comp_rho_min) || !(D0 > 0.0)) return float4(0.0);
    if (K == 5u) {
        float D = max(a.comp_d_min, 0.0);
        float3 gD = float3(0.0);
        [ForceUnroll] for (int s = 0; s < 2; ++s) {
            if (!(rho[s] > 0.0)) continue;
            float g2 = dot(gr[s], gr[s]);
            float Ds = tau[s] - g2 / (8.0 * rho[s]);
            if (Ds <= 0.0) continue;
            D += Ds;
            // d_x |grad rho|^2 = 2 sum_y grad_y rho H[y][x]
            float3 dg2 = 2.0 * (gr[s].x * h0[s] + gr[s].y * h1[s] + gr[s].z * h2[s]);
            gD += gt[s] - dg2 / (8.0 * rho[s]) + g2 * gr[s] / (8.0 * rho[s] * rho[s]);
        }
        float chi = D / D0;
        float q = 1.0 + chi * chi;
        return float4(-2.0 * chi / (q * q) * (gD - chi * gD0) / D0, 1.0 / q);
    }
    float tau_sum = max(tau[0], 0.0) + max(tau[1], 0.0);
    if (!(tau_sum > 0.0)) return float4(0.0, 0.0, 0.0, 1.0);
    // LOL = D0 / (tau + D0): d LOL = (tau dD0 - D0 dtau) / (tau + D0)^2
    float3 gtau = (tau[0] > 0.0 ? gt[0] : float3(0.0)) + (tau[1] > 0.0 ? gt[1] : float3(0.0));
    float den = tau_sum + D0;
    return float4((tau_sum * gD0 - D0 * gtau) / (den * den), D0 / den);
}

// With the functions' ranges in o_A: the fields' ranges, the node's sheet bits and its children (or
// itself as a brick). Every thread of the group calls it (barriers inside). A composite set: field 0's
// bits from its parts' ranges, and wherever it is kept every part's (positive) bit too, so that the bricks
// get their corners (iso_comp_corners clears them again).
void octree_finish(uint4 node, uint ni, uint tid) {
    IsoArgs a = MD_ARGS;
    GroupMemoryBarrierWithGroupSync();
    bool comp = a.comp_kind != 0u;
    if (tid < a.num_fields) {
        uint fbits = (node.w >> (2u * tid)) & 3u;
        if (fbits != 0u || (comp && tid > 0u)) {
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
            if (comp) {
                o_frng[min(tid, 15u)] = rho;
            } else {
                uint bits = 0u;
                if (rho.x <= fd.iso && rho.y >= fd.iso)   bits |= 1u;
                if (rho.x <= -fd.iso && rho.y >= -fd.iso) bits |= 2u;
                bits &= fbits;
                if (bits != 0u) InterlockedOr(o_mask, bits << (2u * tid));
            }
        }
    }
    if (comp) {     // uniform
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u && (node.w & 1u) != 0u) {
            float2 r = comp_range(a);
            float iso = a.fields[0].iso;
            if (r.x <= iso && r.y >= iso) o_mask = 1u | (0x55555554u & ((1u << (2u * min(a.num_fields, 15u))) - 1u));
        }
    }
    octree_emit(node, ni, tid);
}

// A node with sheet bits 'mask': its children (or itself as a brick). One thread. In the order of the
// nodes: the node keeps its sheet bits and the count of what it emits (8 children, or itself at the
// brick level) in scan, whose prefix sum (iso_scan) places them (iso_octree_emit). The nodes of the
// first level come in Morton order, and a node's children in Morton order, so every level, the bricks
// too, is in Morton order: a parent's eight children together, and bricks close in the list close in
// space (what the point list of the corners relies on to keep its groups compact).
void octree_emit_one(uint4 node, uint ni, uint mask) {
    IsoArgs a = MD_ARGS;
    a.nodes_in[ni].w = mask;
    a.scan[ni] = mask == 0u ? 0u : a.level < a.brick_level ? 8u : 1u;
}

// Thread per node of the level: its children (or itself as a brick) at its offset in scan (iso_scan of
// the counts octree_emit_one left); past cap, an overflow (counters[3]: the octree is made again with
// more room)
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_octree_emit(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint ni = node_index(gid) * 64u + tid;
    uint n = min(a.count_in[0], a.cap);
    if (ni >= n) return;
    ni = min(ni, n - 1u);                                  // llvmpipe loads for inactive lanes too
    uint4 node = a.nodes_in[ni];
    if (node.w == 0u) return;
    uint idx = a.scan[ni];
    if (a.level < a.brick_level) {
        if (idx + 8u <= a.cap) {
            for (uint c = 0u; c < 8u; ++c) {
                a.nodes_out[idx + c] = uint4(node.xyz * 2u + uint3(c & 1u, (c >> 1) & 1u, (c >> 2) & 1u), node.w);
            }
        } else {
            InterlockedAdd(a.counters[3], 1u);
        }
    } else {
        if (idx < a.cap) a.nodes_out[idx] = node;
        else InterlockedAdd(a.counters[3], 1u);
    }
}

// Thread per brick: its corner slots (one per field with a sheet bit) into scan, for their prefix sum
// (iso_scan: brick_slot); the bricks with any sheet bit counted in counters[8]
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_brick_slots(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid) * 64u + tid;
    uint n = min(a.count_in[0], a.cap);
    if (b >= n) return;
    b = min(b, n - 1u);                                    // llvmpipe loads for inactive lanes too
    uint m = a.bricks[b].w;
    a.scan[b] = countbits((m | (m >> 1)) & 0x55555555u);
    if (m != 0u) InterlockedAdd(a.counters[8], 1u);
}

// The node's sheet bits in o_mask: its children (or itself as a brick)
void octree_emit(uint4 node, uint ni, uint tid) {
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) octree_emit_one(node, ni, o_mask);
}

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_octree(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint ni = node_index(gid);
    if (ni >= min(a.count_in[0], a.cap)) return;          // uniform per group
    uint4 node = a.nodes_in[ni];
    float3 lo, hi, cen;
    octree_box(node, lo, hi, cen);

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
                if (slice == 0u && col < nf) o_A[col] = octree_range(sm, col);
            }
        }
    }
    octree_finish(node, ni, tid);
}

// The octree for sets of at most O_SMALL functions (orbitals): a thread per shell, its AOs' ranges
// accumulated straight into registers per function, the sums reduced over the group at the end. No
// per-AO arrays in groupshared memory and no barriers per block of shells, so many more groups fit.
#define O_SMALL 4

[shader("compute")]
[numthreads(64, 1, 1)]
void iso_octree_small(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint ni = node_index(gid);
    if (ni >= min(a.count_in[0], a.cap)) return;          // uniform per group
    uint4 node = a.nodes_in[ni];
    float3 lo, hi, cen;
    octree_box(node, lo, hi, cen);
    if (tid == 0u) o_mask = 0u;
    uint nf = min(a.num_fn, uint(O_SMALL));
    float acc[5 * O_SMALL];
    [ForceUnroll] for (int k = 0; k < 5 * O_SMALL; ++k) acc[k] = 0.0;
    for (uint s0 = 0u; s0 < a.num_shells; s0 += 64u) {
        uint s = min(s0 + tid, a.num_shells - 1u);
        float4 sc = a.shells[s].c;
        if (s0 + tid >= a.num_shells || box_dist2(sc.xyz, lo, hi) > sc.w * (1.0 + 16.0 * U)) continue;
        OShell o = octree_shell_setup(s, lo, hi, cen);
        for (uint n = 0u; n < o.ncart; ++n) {
            float mid, rad, pc, dev;
            uint ao;
            octree_ao(o, n, mid, rad, pc, dev, ao, false);
            float gm = abs(mid) + rad + abs(pc) + dev;
            uint row = min(ao, a.num_ao - 1u) * a.fn_stride;
            [ForceUnroll] for (int f = 0; f < O_SMALL; ++f) {
                if (uint(f) < nf) {
                    float c = a.C[row + uint(f)];
                    float ca = abs(c);
                    acc[5 * f + 0] += c * mid;
                    acc[5 * f + 1] += ca * rad;
                    acc[5 * f + 2] += c * pc;
                    acc[5 * f + 3] += ca * dev;
                    acc[5 * f + 4] += ca * gm;
                }
            }
        }
    }
    // the sums over the group: per wave, then over the waves (o_acc as a table of at most 16 waves of
    // at least 4 lanes; a sum in any order is within the bound's rounding allowance)
    uint lanes = WaveGetLaneCount();
    uint w = iso_wave_index(tid, lanes);
    uint nw = min(iso_num_waves(64u, lanes), 16u);
    [ForceUnroll] for (int k = 0; k < 5 * O_SMALL; ++k) {
        float t = WaveActiveSum(acc[k]);
        if (WaveIsFirstLane() && w < 16u) o_acc[w * uint(5 * O_SMALL) + uint(k)] = t;
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid < nf) {
        float sm[5];
        [ForceUnroll] for (int c = 0; c < 5; ++c) {
            float t = 0.0;
            for (uint k = 0u; k < nw; ++k) t += o_acc[k * uint(5 * O_SMALL) + tid * 5u + uint(c)];
            sm[c] = t;
        }
        o_A[tid] = octree_range(sm, tid);
    }
    octree_finish(node, ni, tid);
}

// The octree for densities as phi^T D phi (dform): over the AOs the node's box reaches (the local
// list), with pc the AOs' values at the centre and dev bounds of their deviation from them over the
// box, rho = pc^T D pc + 2 g . delta + delta^T D delta (g = D pc, |delta| <= dev): its range
// rho_c +- (2 sum |g| dev + sum dev |D| dev), intersected with |rho| <= m^T |D| m (m: the AOs' bounds;
// rho >= 0 for a positive semidefinite D, less the slack), plus what the dropped AOs (each below eps
// over the box) can add. No factors, so no limit on the rank. A list that overflows keeps the node.
// Defined after the dform section (it shares the local list): see iso_octree_dform below.

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
        for (uint attempt = 0u; attempt < 64u; ++attempt) {
            InterlockedCompareExchange(a.hash_keys[slot], 0u, key, prev);
            if (prev != 0u) break;
            // Metal's compare-exchange is the weak form and may fail spuriously: the key is in only if
            // the slot holds it now
            uint now;
            InterlockedOr(a.hash_keys[slot], 0u, now);
            if (now != 0u) { prev = now == key ? 0u : now; break; }
        }
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
// of the neighbour across it. When that neighbour is present with every field of this brick (and, in
// a flood, was evaluated by round 'upto'; 0: any), the value is taken from it (iso_corners_fill)
// instead of evaluated twice. Returns the neighbour, or INVALID when the point is the brick's own to
// evaluate; *owner_p: the point's index there.
uint corner_owner_upto(uint4 brick, uint p, out uint owner_p, uint upto) {
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
    bool done = true;
    if (upto != 0u) { uint st = a.bstate[nb == INVALID ? 0u : nb]; done = st <= upto; }   // FLOOD_NONE is above any round
    return (nb != INVALID && (field_bits(nmask) & need) == need && done) ? nb : INVALID;
}

// ------------------------------------------------------------------------------------------
// The lattice point p (0 .. 124) of a brick
float3 corner_point(uint4 brick, uint p) {
    IsoArgs a = MD_ARGS;
    uint k = min(p, 124u);
    uint3 lc = uint3(k % 5u, (k / 5u) % 5u, k / 25u);
    uint3 gc = brick.xyz * 4u + lc;     // integer first: shared corners agree bit for bit
    return a.origin_h.xyz + float3(gc) * a.origin_h.w;
}

// ------------------------------------------------------------------------------------------
// Values at a group of points (corners: consecutive points of the list; vertices: consecutive
// vertices): the shells that reach the group's box are listed (K_LIST at a time) and packed into
// chunks of whole shells, and each chunk's shells, primitives and rows of C (times the AO
// normalization) are staged in groupshared memory, read by every thread at once (broadcast).
// ------------------------------------------------------------------------------------------

#define K_LIST  128                      // shells listed at once (larger groups list K_LIST per pass)
#define K_CK    32                       // AOs per chunk (whole shells)
#define K_PR    256                      // primitives per chunk
groupshared uint  k_sh[K_LIST];          // listed shells
groupshared uint  k_nc[K_LIST];          // their AO counts
groupshared uint  k_np[K_LIST];          // their primitive counts
groupshared uint  k_row[K_LIST];         // their first row within their chunk
groupshared uint  k_poff[K_LIST];        // their first primitive within their chunk
groupshared uint  k_chunk[K_LIST + 1];   // chunk c: listed shells [k_chunk[c], k_chunk[c + 1])
groupshared uint  k_nchunk;
groupshared uint  k_nrows[K_LIST];       // rows per chunk
groupshared uint  k_row_ao[K_CK];        // the current chunk: global AO of each row
groupshared float4 k_cs[K_CK];           // the current chunk's shells: centre, cutoff r^2
groupshared uint  k_ci[K_CK];            // l | ncart << 4 | nprim << 8 | first primitive << 16
groupshared float2 k_pr[K_PR];           // the current chunk's primitives: alpha, coefficient

// A chunk's rows of C (times the AO normalization) [row][F], and the epilogues' staging: one buffer per
// kernel shape, sized for its F so that a narrow kernel keeps its occupancy (a type parameter, so a
// kernel names only its own: the driver allocates every groupshared array a kernel's code names).
// Read as float4 where F allows.
interface IKBuf {
    static float  ld(uint i);
    static float4 ld4(uint i4);
    static void   st(uint i, float v);
}
// Groupshared buffers that threads write a float at a time and read as float4. SPIR-V: float4
// elements (one 16-byte load per float4); Metal: plain floats, since Metal may write one component
// of a threadgroup vector as a read-modify-write of the whole vector, and two threads writing
// different components of the same float4 then lose one of the writes.
#define K_BUF(NAME, FLOATS)                                                         \
    groupshared float4 NAME##_d4[(FLOATS) / 4];                                     \
    groupshared float  NAME##_d1[(FLOATS)];                                         \
    struct NAME : IKBuf {                                                           \
        static float ld(uint i) {                                                   \
            __target_switch { case spirv: return NAME##_d4[i >> 2][i & 3u];         \
                              default:    return NAME##_d1[i]; }                    \
        }                                                                           \
        static float4 ld4(uint i4) {                                                \
            __target_switch { case spirv: return NAME##_d4[i4];                     \
                              default:    return float4(NAME##_d1[4u * i4], NAME##_d1[4u * i4 + 1u], \
                                                        NAME##_d1[4u * i4 + 2u], NAME##_d1[4u * i4 + 3u]); } \
        }                                                                           \
        static void st(uint i, float v) {                                           \
            __target_switch { case spirv: NAME##_d4[i >> 2][i & 3u] = v;            \
                              default:    NAME##_d1[i] = v; }                       \
        }                                                                           \
    };
K_BUF(KBuf1, 64)

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
void k_list(uint s0, float3 lo, float3 hi, uint tid, uint gs, uint ck) {
    IsoArgs a = MD_ARGS;
    // indices clamped into range: some implementations (llvmpipe) load for inactive lanes too
    uint s = min(s0 + tid, a.num_shells - 1u);
    IsoShell sh = a.shells[s];
    bool reach = tid < min(gs, uint(K_LIST)) && s0 + tid < a.num_shells && box_dist2(sh.c.xyz, lo, hi) <= sh.c.w;
    uint nc = sh.info.y;
    uint total;
    uint at = iso_group_prefix(reach, tid, gs, total);
    if (reach) { k_sh[at] = s; k_nc[at] = nc; k_np[at] = sh.info.z; }
    GroupMemoryBarrierWithGroupSync();
    if (tid == 0u) {        // greedy packing over groupshared memory: at most gs steps
        uint c = 0u, rows = 0u, prims = 0u;
        k_chunk[0] = 0u;
        for (uint i = 0u; i < total; ++i) {
            uint n = k_nc[i], np = k_np[i];
            if (rows > 0u && (rows + n > ck || prims + np > K_PR)) { k_nrows[c] = rows; c++; k_chunk[c] = i; rows = 0u; prims = 0u; }
            k_row[i] = rows;
            k_poff[i] = prims;
            rows += n;
            prims += np;
        }
        if (total > 0u) { k_nrows[c] = rows; c++; }
        k_chunk[c] = total;
        k_nchunk = c;
    }
    GroupMemoryBarrierWithGroupSync();
}

// Stages chunk c: its shells (k_cs, k_ci), their primitives (k_pr) and its rows of C times the AO
// normalization, columns [f0, f0 + F) (all threads; ends with a barrier). The evaluation then reads
// groupshared memory only.
void k_load_chunk<let F : int, B : IKBuf>(uint c, uint f0, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
    uint i = min(i0 + tid, i1 - 1u);       // clamped: see k_list
    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
    uint ao0 = sh.info2.x;
    if (tid < i1 - i0) {
        uint np = sh.info.z, po = k_poff[i];
        for (uint n = 0u; n < k_nc[i]; ++n) k_row_ao[k_row[i] + n] = ao0 + n;
        k_cs[tid] = sh.c;
        k_ci[tid] = sh.info.x | (sh.info.y << 4) | (np << 8) | (po << 16);
        for (uint q = 0u; q < np; ++q) k_pr[min(po + q, uint(K_PR) - 1u)] = a.prims[sh.info.w + q];
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

// Chunk shell li (k_load_chunk) at p: d = p - centre, the radial part R (and R1 = dR/d(r^2)), l and
// the AO count. False, with R = 0, beyond the shell's cutoff (every AO below the evaluation epsilon).
bool k_shell(uint li, float3 p, out float3 d, out float R, out uint l, out uint nc) {
    float4 cs = k_cs[li];
    uint ci = k_ci[li];
    d = p - cs.xyz;
    float r2 = dot(d, d);
    l = ci & 15u;
    nc = (ci >> 4) & 15u;
    R = 0.0;
    if (r2 > cs.w) return false;
    uint np = (ci >> 8) & 255u, po = ci >> 16;
    for (uint q = 0u; q < np; ++q) {
        float2 pr = k_pr[po + q];
        R += pr.y * exp(-pr.x * r2);
    }
    return true;
}

bool k_shell_d(uint li, float3 p, out float3 d, out float R, out float R1, out uint l, out uint nc) {
    float4 cs = k_cs[li];
    uint ci = k_ci[li];
    d = p - cs.xyz;
    float r2 = dot(d, d);
    l = ci & 15u;
    nc = (ci >> 4) & 15u;
    R = 0.0;
    R1 = 0.0;
    if (r2 > cs.w) return false;
    uint np = (ci >> 8) & 255u, po = ci >> 16;
    for (uint q = 0u; q < np; ++q) {
        float2 pr = k_pr[po + q];
        float t = pr.y * exp(-pr.x * r2);
        R += t;
        R1 -= pr.x * t;
    }
    return true;
}

// Shared corners: the points a brick took from its neighbours, copied from their slots
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_corners_fill(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count || tid >= 125u) return;
    // flood (mode bit 8): only the round's bricks, from neighbours evaluated by then; coarse (bit 16): the
    // brick corners only
    uint upto = 0u;
    if ((a.mode & 8u) != 0u) { if (a.bstate[b] != a.level) return; upto = a.level; }
    if ((a.mode & 16u) != 0u && !iso_brick_corner(tid)) return;
    uint4 brick = a.bricks[b];
    uint op;
    uint nb = corner_owner_upto(brick, tid, op, upto);
    if (nb == INVALID) return;
    uint nmask = a.bricks[nb].w;
    uint slot = a.brick_slot[b], nslot0 = a.brick_slot[nb];
    for (uint f = 0u; f < a.num_fields; ++f) {
        if (((field_bits(brick.w) >> (2u * f)) & 1u) == 0u) continue;
        a.corners[slot * 125u + tid] = a.corners[(nslot0 + field_slot(nmask, f)) * 125u + op];
        slot++;
    }
}

// A composite set, after iso_corners_fill: field 0's corners from its parts', then the parts' bits cleared
// (their slots stay, unused), so that the cells, quads and vertices are field 0's only. A group per brick.
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_comp_corners(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    uint mask = a.bricks[b].w;
    if ((mask & 3u) != 0u && tid < 125u) {
        uint s0 = a.brick_slot[b];
        float4 p[COMP_MAX];
        [ForceUnroll] for (int i = 0; i < COMP_MAX; ++i) {
            uint f = uint(i) + 1u;
            bool on = f < a.num_fields && ((field_bits(mask) >> (2u * f)) & 1u) != 0u;
            float c = a.corners[(s0 + field_slot(mask, min(f, a.num_fields - 1u))) * 125u + tid];   // clamped: llvmpipe loads anyway
            p[i] = float4(0.0, 0.0, 0.0, on ? c : 0.0);
        }
        a.corners[(s0 + field_slot(mask, 0u)) * 125u + tid] = comp_eval(a, p).w;
    }
    GroupMemoryBarrierWithGroupSync();                     // every thread has its mask
    if (tid == 0u) a.bricks[b].w = mask & 3u;
}

// A composite set's vertices: field 0's value and gradient (vval) from its parts' (vsub), where the vertex is
// still moving (as the evaluators: a converged one keeps its values)
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_comp_vertices(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint v = tid.x;
    if (v >= a.count || a.vdone[v] != 0u) return;
    float4 p[COMP_MAX];
    [ForceUnroll] for (int i = 0; i < COMP_MAX; ++i) {
        bool on = uint(i) + 1u < a.num_fields;
        float4 q = a.vsub[min(uint(i), a.num_fields - 2u) * a.comp_stride + v];   // clamped: llvmpipe loads anyway
        p[i] = on ? q : float4(0.0);
    }
    a.vval[v] = comp_eval(a, p);
}

// ------------------------------------------------------------------------------------------
// Corners from a point list: the lattice points to evaluate, brick by brick (all 125 of a brick with
// sheets, or with shared corners those not taken from a neighbour), and groups of 256 consecutive
// points - 2 to 4 neighbouring bricks - evaluated together: every thread busy, and the setup of each
// shell shared by 256 points instead of one brick's
// ------------------------------------------------------------------------------------------

// Point p of a brick (0 .. 124) is one of its 8 corners
bool iso_brick_corner(uint p) {
    uint k = min(p, 124u);
    uint3 lc = uint3(k % 5u, (k / 5u) % 5u, k / 25u);
    return all(lc == uint3(0u) || lc == uint3(4u));
}

// ------------------------------------------------------------------------------------------
// Seed and flood (approximate: a surface piece no seed reaches is missed). The candidate bricks'
// corners are evaluated first; a brick whose corners straddle an iso value is a seed (round 1). Each
// round's bricks are evaluated in full, and a candidate across a face of one whose 25 face points
// hold a sign change on a face edge joins the next round. A flooded brick has surface cells.
// ------------------------------------------------------------------------------------------

// A brick's state before a round reaches it. The rounds claim bricks with InterlockedMin (the state only
// falls), not a compare-exchange: Metal's is the weak form, which may fail spuriously.
#define FLOOD_NONE 0xFFFFFFFFu

// Sign of the brick's point (x, y, z) for field f, sheet s against the iso value: above
bool flood_above(uint slot, uint x, uint y, uint z, float sg, float iso) {
    IsoArgs a = MD_ARGS;
    return sg * a.corners[slot * 125u + x + 5u * y + 25u * z] - iso > 0.0;
}

// Thread per brick: the seeds (counted in counters[6])
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_flood_seed(uint3 tid3 : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint b = tid3.x;
    if (b >= a.count) return;
    uint4 brick = a.bricks[b];
    uint seed = 0u;
    uint slot = a.brick_slot[b];
    for (uint f = 0u; f < a.num_fields; ++f) {
        uint fb = (brick.w >> (2u * f)) & 3u;
        if (fb == 0u) continue;
        float iso = a.fields[f].iso;
        for (uint sh = 0u; sh < 2u; ++sh) {
            if ((fb & (1u << sh)) == 0u) continue;
            float sg = sh == 0u ? 1.0 : -1.0;
            uint pos = 0u;
            for (uint k = 0u; k < 8u; ++k) pos += flood_above(slot, (k & 1u) * 4u, ((k >> 1) & 1u) * 4u, (k >> 2) * 4u, sg, iso) ? 1u : 0u;
            if (pos != 0u && pos != 8u) seed = 1u;
        }
        slot++;
    }
    a.bstate[b] = seed != 0u ? 1u : FLOOD_NONE;
    if (seed != 0u) InterlockedAdd(a.counters[6], 1u);
}

// Thread per shell: the candidate brick holding the shell's centre (an atom) is a seed too (after
// iso_flood_seed; counted in counters[6]). A total density's pieces hold an atom each (its maxima are
// at the nuclei), so a piece too small for the brick corners to straddle is still reached.
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_flood_atoms(uint3 tid3 : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint s = tid3.x;
    if (s >= a.num_shells) return;
    s = min(s, a.num_shells - 1u);                         // llvmpipe loads for inactive lanes too
    float3 q = (a.shells[s].c.xyz - a.origin_h.xyz) / (4.0 * a.origin_h.w);
    if (any(q < float3(0.0)) || any(q >= float3(1024.0))) return;
    uint nb = brick_lookup(int3(floor(q)));
    if (nb == INVALID || a.bricks[nb].w == 0u) return;
    uint orig;
    InterlockedMin(a.bstate[nb], 1u, orig);
    if (orig == FLOOD_NONE) InterlockedAdd(a.counters[6], 1u);
}

groupshared uint fl_face;

// Group per brick of round a.level: its faces with a sign change on a face edge; the candidates across
// them not yet flooded join round a.level + 1 (counted in counters[6])
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_flood_face(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    if (a.bstate[b] != a.level) return;                    // uniform per group
    uint4 brick = a.bricks[b];
    if (tid == 0u) fl_face = 0u;
    GroupMemoryBarrierWithGroupSync();
    // thread: face tid / 10 (6 faces), its rows of 25 points
    uint face = tid / 10u, part = tid % 10u;
    if (face < 6u) {
        uint ax = face >> 1, side = face & 1u;
        uint u = (ax + 1u) % 3u, v = (ax + 2u) % 3u;
        bool change = false;
        uint slot = a.brick_slot[b];
        for (uint f = 0u; f < a.num_fields; ++f) {
            uint fb = (brick.w >> (2u * f)) & 3u;
            if (fb == 0u) continue;
            float iso = a.fields[f].iso;
            for (uint sh = 0u; sh < 2u; ++sh) {
                if ((fb & (1u << sh)) == 0u) continue;
                float sg = sh == 0u ? 1.0 : -1.0;
                for (uint e = part; e < 25u; e += 10u) {
                    uint i = e % 5u, j = e / 5u;
                    uint3 p = uint3(0u); p[ax] = side * 4u; p[u] = i; p[v] = j;
                    bool s0 = flood_above(slot, p.x, p.y, p.z, sg, iso);
                    if (i < 4u) { uint3 q = p; q[u] += 1u; change = change || flood_above(slot, q.x, q.y, q.z, sg, iso) != s0; }
                    if (j < 4u) { uint3 q = p; q[v] += 1u; change = change || flood_above(slot, q.x, q.y, q.z, sg, iso) != s0; }
                }
            }
            slot++;
        }
        if (change) InterlockedOr(fl_face, 1u << face);
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid < 6u && ((fl_face >> tid) & 1u) != 0u) {
        uint ax = tid >> 1, side = tid & 1u;
        int3 nbc = int3(brick.xyz);
        nbc[ax] += side != 0u ? 1 : -1;
        uint nb = brick_lookup(nbc);
        if (nb != INVALID && a.bricks[nb].w != 0u) {
            uint orig;
            InterlockedMin(a.bstate[nb], a.level + 1u, orig);   // a brick already reached keeps its round
            if (orig == FLOOD_NONE) InterlockedAdd(a.counters[6], 1u);
        }
    }
}

// Thread per brick: the bricks the flood never reached lose their sheet bits (no cells)
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_flood_finish(uint3 tid3 : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint b = tid3.x;
    if (b >= a.count) return;
    if (a.bstate[b] == FLOOD_NONE) a.bricks[b].w = 0u;
}


// One group per brick: its points into plist, in the order of the bricks. mode bit 2: the count, per
// brick into scan (and the total into counters[5]); else the points at the brick's offset in scan
// (iso_scan of the counts). mode bit 1: shared corners.
[shader("compute")]
[numthreads(128, 1, 1)]
void iso_points(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint b = node_index(gid);
    if (b >= a.count) return;                              // uniform per group
    uint4 brick = a.bricks[b];
    bool mine = tid < 125u && brick.w != 0u;
    // flood (mode bit 8): the round's bricks, sharing with neighbours evaluated by then; coarse (bit 16):
    // the brick corners of every brick
    uint upto = 0u;
    if ((a.mode & 8u) != 0u) { mine = mine && a.bstate[b] == a.level; upto = a.level; }
    if ((a.mode & 16u) != 0u) mine = mine && iso_brick_corner(tid);
    if ((a.mode & 1u) != 0u && mine) {
        uint op;
        mine = corner_owner_upto(brick, tid, op, upto) == INVALID;
    }
    uint n;
    uint at = iso_group_prefix(mine, tid, 128u, n);
    if ((a.mode & 2u) != 0u) {
        if (tid == 0u) {
            a.scan[b] = n;
            if (n > 0u) InterlockedAdd(a.counters[5], n);
        }
        return;
    }
    if (mine) a.plist[a.scan[b] + at] = (b << 7) | tid;
}

// ------------------------------------------------------------------------------------------
// Exclusive prefix sum over scan[0, count) in place, the total in scan[count]: one group of
// SC_GROUP threads, SC_PER consecutive entries each per pass
// ------------------------------------------------------------------------------------------

#define SC_GROUP 512
#define SC_PER   4

[shader("compute")]
[numthreads(SC_GROUP, 1, 1)]
void iso_scan(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    // count ~0: as many as the device says (count_in, at most cap), the sum also to count_out
    uint n = a.count;
    if (n == 0xFFFFFFFFu) n = min(a.count_in[0], a.cap);
    uint run = 0u;
    for (uint base = 0u; base < n; base += uint(SC_GROUP * SC_PER)) {
        uint i0 = base + tid * uint(SC_PER);
        uint v[SC_PER];
        uint sum = 0u;
        [ForceUnroll] for (int k = 0; k < SC_PER; ++k) {
            v[k] = i0 + uint(k) < n ? a.scan[i0 + uint(k)] : 0u;
            sum += v[k];
        }
        uint tot;
        uint pre = iso_group_prefix_sum(sum, tid, uint(SC_GROUP), tot);   // barriers inside: earlier reads done
        uint at = run + pre;
        [ForceUnroll] for (int k = 0; k < SC_PER; ++k) {
            if (i0 + uint(k) < n) a.scan[i0 + uint(k)] = at;
            at += v[k];
        }
        run += tot;
    }
    if (tid == 0u) {
        a.scan[n] = run;
        if (a.count == 0xFFFFFFFFu) a.count_out[0] = run;
    }
}

// One function, four points per thread: 64 threads take 256 consecutive points of the list (four
// consecutive per thread), so a shell's staged data, primitives and coefficients are read once for
// four points, as the volume path's kernel shares them over the voxels of a thread.
#define L4_GROUP 64
#define L4_VPT   4
groupshared float4 l_lohi[2 * 32];       // per wave: the points' box

[shader("compute")]
[numthreads(L4_GROUP, 1, 1)]
void iso_corners_list_f1v4(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint total = a.count_in[0];
    uint e0 = node_index(gid) * uint(L4_GROUP * L4_VPT);
    if (e0 >= total) return;                               // uniform per group
    float3 p[L4_VPT];
    uint pt[L4_VPT], b[L4_VPT];
    bool act[L4_VPT];
    float3 lo = float3(3.0e38), hi = float3(-3.0e38);
    [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) {
        uint e = e0 + tid * uint(L4_VPT) + uint(v);
        act[v] = e < total;
        uint ent = a.plist[min(e, total - 1u)];
        b[v] = ent >> 7;
        pt[v] = ent & 127u;
        p[v] = corner_point(a.bricks[b[v]], pt[v]);
        if (act[v]) { lo = min(lo, p[v]); hi = max(hi, p[v]); }
    }
    {
        uint lanes = WaveGetLaneCount();
        uint w = min(iso_wave_index(tid, lanes), 31u);
        float3 wlo = WaveActiveMin(lo), whi = WaveActiveMax(hi);
        if (WaveIsFirstLane()) { l_lohi[2u * w] = float4(wlo, 0.0); l_lohi[2u * w + 1u] = float4(whi, 0.0); }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) {
            uint nw = min(iso_num_waves(L4_GROUP, lanes), 32u);
            float3 glo = l_lohi[0].xyz, ghi = l_lohi[1].xyz;
            for (uint k = 1u; k < nw; ++k) { glo = min(glo, l_lohi[2u * k].xyz); ghi = max(ghi, l_lohi[2u * k + 1u].xyz); }
            l_lohi[0] = float4(glo, 0.0); l_lohi[1] = float4(ghi, 0.0);
        }
        GroupMemoryBarrierWithGroupSync();
    }
    float3 glo = l_lohi[0].xyz, ghi = l_lohi[1].xyz;
    float acc[L4_VPT];
    [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) acc[v] = 0.0;
    for (uint s0 = 0u; s0 < a.num_shells; s0 += uint(L4_GROUP)) {
        GroupMemoryBarrierWithGroupSync();                 // the previous block's lists consumed
        k_list(s0, glo, ghi, tid, L4_GROUP, K_CK);
        uint nchunk = k_nchunk;
        for (uint c = 0u; c < nchunk; ++c) {
            k_load_chunk<1, KBuf1>(c, 0u, tid, L4_GROUP);
            uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
            for (uint i = i0; i < i1; ++i) {               // uniform
                uint li = i - i0;
                float4 cs = k_cs[li];
                uint ci = k_ci[li];
                uint l = ci & 15u, nc = (ci >> 4) & 15u, np = (ci >> 8) & 255u, po = ci >> 16;
                float3 d[L4_VPT];
                float r2[L4_VPT], R[L4_VPT];
                bool any = false;
                [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) {
                    d[v] = p[v] - cs.xyz;
                    r2[v] = dot(d[v], d[v]);
                    R[v] = 0.0;
                    any = any || r2[v] <= cs.w;
                }
                if (!any) continue;                         // none of the thread's points in reach
                for (uint q = 0u; q < np; ++q) {
                    float2 pr = k_pr[po + q];
                    [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) R[v] += pr.y * exp(-pr.x * r2[v]);
                }
                // the coefficients once (rows of KBuf1), the polynomial per point
                float cf[15];
                uint row0 = k_row[i];
                [ForceUnroll] for (int n = 0; n < 15; ++n) cf[n] = uint(n) < nc ? KBuf1::ld(row0 + uint(n)) : 0.0;
                [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) {
                    if (r2[v] > cs.w) continue;             // beyond the cutoff: dropped, as for one point a thread
                    float m[15];
                    cart_monomials(l, d[v], m);
                    float P = 0.0;
                    [ForceUnroll] for (int n = 0; n < 15; ++n) {
                        if (uint(n) < nc) P += m[n] * cf[n];
                    }
                    acc[v] += P * R[v];
                }
            }
            GroupMemoryBarrierWithGroupSync();             // the chunk's rows consumed
        }
    }
    [ForceUnroll] for (int v = 0; v < L4_VPT; ++v) {
        if (!act[v]) continue;
        uint4 brick = a.bricks[b[v]];
        uint fbits = field_bits(brick.w);
        uint s = a.brick_slot[b[v]];
        for (uint f = 0u; f < a.num_fields; ++f) {
            if (((fbits >> (2u * f)) & 1u) == 0u) continue;
            IsoTerm tm = a.terms[a.fields[f].term_off];
            a.corners[s * 125u + pt[v]] = tm.w * acc[v];
            s++;
        }
    }
}

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

// ------------------------------------------------------------------------------------------
// Corners and vertices as a tiled GEMM on the GEMM core (gto_gemm_core.slang, staged): 128 evaluation rows per group
// (corners: 128 consecutive points of the list; vertices: 32 vertices x 4 channels, value and gradient), F function
// columns per tile. Per chunk of at most G_CK AOs (whole shells) every thread puts its row's AO values in the core's
// A chunk ([kk][r], r = tid) and the chunk's rows of C times the AO normalization in its B chunk ([kk][F]); thread
// (pg, cg) = (tid >> 2, tid & 3), the core's (ty, tx), multiplies rows 4 pg .. 4 pg + 3 by the F / 4 columns cg F / 4
// .. (R = 4, M = F / 4). The epilogues reduce the four column groups' partial sums per row.
// ------------------------------------------------------------------------------------------

#define G_GROUP 128
#define G_CK    16
groupshared float  g_red[4 * G_GROUP];   // epilogues: per column group and row
groupshared float4 g_lohi[2 * 32];       // per wave: the rows' box
groupshared uint   g_done;

#include "../gto/gto_gemm_core.slang"
GC_MEM4(IsoGc16, 128, 16)
GC_MEM4(IsoGc32, 128, 32)
GC_MEM4(IsoGc64, 128, 64)

// k_load_chunk with the rows of C into the core's B chunk
void k_load_chunk_gc<let F : int, Mem : IGcMem>(uint c, uint f0, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
    uint i = min(i0 + tid, i1 - 1u);       // clamped: see k_list
    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
    uint ao0 = sh.info2.x;
    if (tid < i1 - i0) {
        uint np = sh.info.z, po = k_poff[i];
        for (uint n = 0u; n < k_nc[i]; ++n) k_row_ao[k_row[i] + n] = ao0 + n;
        k_cs[tid] = sh.c;
        k_ci[tid] = sh.info.x | (sh.info.y << 4) | (np << 8) | (po << 16);
        for (uint q = 0u; q < np; ++q) k_pr[min(po + q, uint(K_PR) - 1u)] = a.prims[sh.info.w + q];
    }
    GroupMemoryBarrierWithGroupSync();
    uint nrows = k_nrows[c];
    for (uint e = tid; e < nrows * uint(F); e += gs) {
        uint row = e / uint(F), col = e % uint(F);
        uint ao = min(k_row_ao[min(row, uint(K_CK) - 1u)], a.num_ao - 1u);
        float v = a.C[ao * a.fn_stride + min(f0 + col, a.fn_stride - 1u)] * a.ao_nrm[ao];
        Mem::st4c(uint(GK * G_GROUP / 4) + (e >> 2), e & 3u, f0 + col < a.num_fn ? v : 0.0);
    }
    GroupMemoryBarrierWithGroupSync();
}

// The chunks: the shells that reach the group's box, K_LIST listed at a time, packed into chunks of whole shells of at
// most G_CK rows; each staged and its AO values made at the thread's row
struct IsoGemmStager<let F : int, Mem : IGcMem, let VERT : int> : IGemmStager {
    static const int kFull = 0;
    IsoArgs a;                                             // the kernel's (not read from MD_ARGS again per chunk)
    float3 glo, ghi, p;
    uint   ch, ax;
    uint   s0, c, nchunk;
    bool   listed;
    [mutating] void begin(uint m0, uint tid) { s0 = 0u; c = 0u; nchunk = 0u; listed = false; }
    [mutating] uint stage(uint m0, uint tid) {
        while (c >= nchunk) {                              // uniform
            if (listed) s0 += uint(K_LIST);
            if (s0 >= a.num_shells) return GC_END;
            k_list(s0, glo, ghi, tid, G_GROUP, G_CK);      // barriers inside
            listed = true;
            nchunk = k_nchunk;
            c = 0u;
        }
        k_load_chunk_gc<F, Mem>(c, m0, tid, G_GROUP);
        uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
        // this row's AO values of the chunk's shells (zero beyond a shell's cutoff)
        for (uint i = i0; i < i1; ++i) {
            float3 d;
            float R, R1 = 0.0;
            uint l, nc;
            if (VERT == 0) k_shell(i - i0, p, d, R, l, nc);
            else           k_shell_d(i - i0, p, d, R, R1, l, nc);
            float m[15];
            cart_monomials(l, d, m);
            float A = R, Bg = 0.0;
            float g[15];
            if (VERT != 0) {
                cart_monomials_d(l, d, ax, g);
                // channel 0: m R; channel 1 + ax: g R + m 2 d[ax] R1
                float da = ax == 0u ? d.x : (ax == 1u ? d.y : d.z);
                A = ch == 0u ? R : 2.0 * da * R1;
                Bg = ch == 0u ? 0.0 : R;
            }
            uint rb = k_row[i] * uint(G_GROUP / 4) + (tid >> 2);   // row k_row[i] at float4 rb, component tid & 3
            [ForceUnroll] for (int n = 0; n < 15; ++n) {
                if (uint(n) < nc) {                         // uniform
                    float phi = m[n] * A;
                    if (VERT != 0) phi += g[n] * Bg;
                    Mem::st4c(rb + uint(n) * uint(G_GROUP / 4), tid & 3u, phi);
                }
            }
        }
        const uint nrows = k_nrows[c];
        c += 1u;
        return nrows;
    }
};

// The tiled GEMM's epilogues: X[v][c] the rows 4 pg + v, columns m0 + cg FC + c
struct IsoGemmEpilogue<let F : int, Mem : IGcMem, let VERT : int> : IGemmEpilogue {
    uint   pg, cg, pt, slot0, mask, v0, vl, ch, stride;
    bool   active;
    float  res;
    float4 res4;
    [mutating] void take<let R : int, let M : int>(uint f0, uint ncols, uint tx, uint ty, float X[R][M]) {
        IsoArgs a = MD_ARGS;
        const int FC = M;
        const uint tid = pg * 4u + cg;
        uint col0 = f0 + cg * uint(FC);
        if (a.canon != 0u) {
            for (uint f = 0u; f < (VERT == 0 ? a.num_fields : 1u); ++f) {
                uint fw = f;
                if (VERT != 0) fw = vertex_field(min(v0 + pg, a.count - 1u));   // the vertex's field
                uint wb = fw * 3u * stride + col0;
                float part[4] = {0.0, 0.0, 0.0, 0.0};
                [ForceUnroll] for (int c = 0; c < FC; ++c) {
                    float wl = a.wts[wb + uint(c)], wd = a.wts[wb + stride + uint(c)];
                    float wp = (c & 1) == 0 ? a.wts[wb + 2u * stride + uint(c)] : 0.0;
                    if (VERT == 0) {
                        [ForceUnroll] for (int v = 0; v < 4; ++v) {
                            float x = X[v][c];
                            part[v] += x * (wl + wd * x);
                            if ((c & 1) == 0) part[v] += wp * x * X[v][c + 1];
                        }
                    } else {
                        float x = X[0][c];
                        part[0] += x * (wl + wd * x);
                        [ForceUnroll] for (int k = 1; k < 4; ++k) part[k] += X[k][c] * (wl + 2.0 * wd * x);
                        if ((c & 1) == 0) {
                            float x1 = X[0][c + 1];
                            part[0] += wp * x * x1;
                            [ForceUnroll] for (int k = 1; k < 4; ++k) part[k] += wp * (x * X[k][c + 1] + x1 * X[k][c]);
                        }
                    }
                }
                GroupMemoryBarrierWithGroupSync();         // g_red free
                [ForceUnroll] for (int v = 0; v < 4; ++v) g_red[cg * uint(G_GROUP) + pg * 4u + uint(v)] = part[v];
                GroupMemoryBarrierWithGroupSync();
                float sum = g_red[tid] + g_red[uint(G_GROUP) + tid] + g_red[2u * uint(G_GROUP) + tid] + g_red[3u * uint(G_GROUP) + tid];
                if (VERT == 0) {
                    if (active && ((field_bits(mask) >> (2u * f)) & 1u) != 0u) {
                        uint ci = (slot0 + field_slot(mask, f)) * 125u + pt;
                        a.corners[ci] = f0 == 0u ? sum : a.corners[ci] + sum;
                    }
                } else {
                    res += sum;
                }
            }
        } else {
            // other terms: the rows staged RR at a time (in the A chunk's memory), then the terms of each row's field(s)
            const uint RR = uint(4 * G_CK * 32) / uint(F);
            for (uint r0 = 0u; r0 < uint(G_GROUP); r0 += RR) {
                GroupMemoryBarrierWithGroupSync();
                [ForceUnroll] for (int v = 0; v < 4; ++v) {
                    uint rr = pg * 4u + uint(v) - r0;
                    if (rr < RR) {
                        [ForceUnroll] for (int c = 0; c < FC; ++c) Mem::st4c(rr * uint(F / 4) + cg * uint(FC / 4) + uint(c / 4), uint(c % 4), X[v][c]);
                    }
                }
                GroupMemoryBarrierWithGroupSync();
                uint rr = tid - r0;
                if (rr >= RR || !active) continue;
                if (VERT == 0) {
                    uint rb = rr * uint(F);
                    for (uint f = 0u; f < a.num_fields; ++f) {
                        if (((field_bits(mask) >> (2u * f)) & 1u) == 0u) continue;
                        IsoField fd = a.fields[f];
                        float v = 0.0;
                        for (uint t = 0u; t < fd.term_count; ++t) {
                            IsoTerm tm = a.terms[fd.term_off + t];
                            if (tm.a < f0 || tm.a >= f0 + uint(F)) continue;
                            float xa = Mem::ld(rb + tm.a - f0);
                            if (tm.b == INVALID) { v += tm.w * xa; continue; }
                            v += tm.w * xa * Mem::ld(rb + tm.b - f0);
                        }
                        uint ci = (slot0 + field_slot(mask, f)) * 125u + pt;
                        a.corners[ci] = f0 == 0u ? v : a.corners[ci] + v;
                    }
                } else if (ch == 0u) {
                    uint rb = rr * uint(F);                // rows rr .. rr + 3: value, d/dx, d/dy, d/dz
                    IsoField fd = a.fields[vertex_field(v0 + vl)];
                    float4 r = float4(0.0);
                    for (uint t = 0u; t < fd.term_count; ++t) {
                        IsoTerm tm = a.terms[fd.term_off + t];
                        if (tm.a < f0 || tm.a >= f0 + uint(F)) continue;
                        uint ia = tm.a - f0;
                        float xa = Mem::ld(rb + ia);
                        float3 ga = float3(Mem::ld(rb + uint(F) + ia), Mem::ld(rb + 2u * uint(F) + ia), Mem::ld(rb + 3u * uint(F) + ia));
                        if (tm.b == INVALID) { r += tm.w * float4(ga, xa); continue; }
                        uint ib = tm.b - f0;
                        float xb = Mem::ld(rb + ib);
                        float3 gb = float3(Mem::ld(rb + uint(F) + ib), Mem::ld(rb + 2u * uint(F) + ib), Mem::ld(rb + 3u * uint(F) + ib));
                        r += tm.w * float4(xa * gb + xb * ga, xa * xb);
                    }
                    res4 += r;
                }
            }
        }
    }
};

void gemm_eval<let F : int, Mem : IGcMem, let VERT : int>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    uint pg = tid >> 2, cg = tid & 3u;
    uint b = 0u, pt = 0u, slot0 = 0u, mask = 0u, v0 = 0u;
    uint vl = tid >> 2, ch = tid & 3u;
    bool active;
    float3 p;
    if (VERT == 0) {
        uint total = a.count_in[0];
        uint e0 = node_index(gid) * uint(G_GROUP);
        if (e0 >= total) return;                           // uniform per group
        uint e = e0 + tid;
        active = e < total;
        uint ent = a.plist[min(e, total - 1u)];
        b = ent >> 7;
        pt = ent & 127u;
        uint4 brick = a.bricks[b];
        mask = brick.w;
        slot0 = a.brick_slot[b];
        p = corner_point(brick, pt);
    } else {
        v0 = node_index(gid) * 32u;
        if (v0 >= a.count) return;                         // uniform per group
        active = v0 + vl < a.count;
        if (tid == 0u) {
            uint done = 1u;
            for (uint k = 0u; k < 32u && v0 + k < a.count; ++k) done &= a.vdone[v0 + k] != 0u ? 1u : 0u;   // 1 converged, 2 predicted
            g_done = done;
        }
        GroupMemoryBarrierWithGroupSync();
        if (g_done != 0u) return;                          // uniform: every vertex here has converged
        p = a.vpos[min(v0 + vl, a.count - 1u)].xyz;
    }
    // the group's box: its rows' points
    {
        uint lanes = WaveGetLaneCount();
        uint w = min(iso_wave_index(tid, lanes), 31u);
        float3 lo = WaveActiveMin(active ? p : float3(3.0e38));
        float3 hi = WaveActiveMax(active ? p : float3(-3.0e38));
        if (WaveIsFirstLane()) { g_lohi[2u * w] = float4(lo, 0.0); g_lohi[2u * w + 1u] = float4(hi, 0.0); }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) {
            uint nw = min(iso_num_waves(G_GROUP, lanes), 32u);
            float3 glo = g_lohi[0].xyz, ghi = g_lohi[1].xyz;
            for (uint k = 1u; k < nw; ++k) { glo = min(glo, g_lohi[2u * k].xyz); ghi = max(ghi, g_lohi[2u * k + 1u].xyz); }
            g_lohi[0] = float4(glo, 0.0);
            g_lohi[1] = float4(ghi, 0.0);
        }
        GroupMemoryBarrierWithGroupSync();
    }
    IsoGemmStager<F, Mem, VERT> st;
    st.a = a;
    st.glo = g_lohi[0].xyz;
    st.ghi = g_lohi[1].xyz;
    st.p   = p;
    st.ch  = ch;
    st.ax  = ch > 0u ? ch - 1u : 0u;
    st.s0 = 0u; st.c = 0u; st.nchunk = 0u; st.listed = false;
    IsoGemmEpilogue<F, Mem, VERT> epi;
    epi.pg = pg; epi.cg = cg; epi.pt = pt; epi.slot0 = slot0; epi.mask = mask; epi.v0 = v0; epi.vl = vl; epi.ch = ch;
    epi.stride = a.fn_stride;
    epi.active = active;
    epi.res = 0.0;
    epi.res4 = float4(0.0);

    gemm_core_staged<G_GROUP, F, 4, 32, Mem, IsoGemmStager<F, Mem, VERT>, IsoGemmEpilogue<F, Mem, VERT> >(st, epi, a.num_fn, tid, cg, pg);

    if (VERT != 0) {
        if (a.canon != 0u) {
            GroupMemoryBarrierWithGroupSync();
            g_red[tid] = epi.res;
            GroupMemoryBarrierWithGroupSync();
            if (tid < 32u && v0 + tid < a.count) a.vval[v0 + tid] = float4(g_red[4u * tid + 1u], g_red[4u * tid + 2u], g_red[4u * tid + 3u], g_red[4u * tid]);
        } else if (ch == 0u && active) {
            a.vval[v0 + vl] = epi.res4;
        }
    }
}

[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_gcorners_f16(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<16, IsoGc16, 0>(gid, tid); }
[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_gcorners_f32(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<32, IsoGc32, 0>(gid, tid); }
[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_gcorners_f64(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<64, IsoGc64, 0>(gid, tid); }
[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_geval_f16(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<16, IsoGc16, 1>(gid, tid); }
[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_geval_f32(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<32, IsoGc32, 1>(gid, tid); }
[shader("compute")] [numthreads(G_GROUP, 1, 1)] void iso_geval_f64(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_eval<64, IsoGc64, 1>(gid, tid); }

// ------------------------------------------------------------------------------------------
// The wide GEMM (canonical terms only), on the GEMM core: 64 rows per group of 256 threads and F = 128, 192 or 256
// function columns per pass, so a set of up to 256 functions lists its shells, stages its chunks and makes its AO
// values once (the tiled GEMM: once per tile of 64). Thread (pg, cg), pg = tid % 16, cg = tid / 16, the core's
// (ty, tx), multiplies rows 4 pg .. 4 pg + 3 by the F / 16 columns cg F / 16 .. (R = 4, M = F / 16); the AO values of
// row 4 pg + cg % 4 are made by the four threads cg / 4 = 0 .. 3, each for every fourth shell of the chunk (the
// shells cover the chunk's rows once each). Rows: corners, 64 consecutive points of the list; vertices, 16 vertices
// x 4 channels. Results: rows 0 .. 63 by threads 0 .. 63.
// ------------------------------------------------------------------------------------------

#define W_GROUP 256
#define W_ROWS  64

// ------------------------------------------------------------------------------------------
// The wide GEMM's second staging (IsoWide2Stager): the group's shells listed once, all of them in rounds of the
// group's threads (up to KW_LIST shells and KW_ROWS AOs a pass), and staged as they are listed: centres, primitives
// (up to KW_PR in groupshared memory, the rest read globally) and each AO's global index. The first (IsoWideStager)
// lists 128 shells of the set a pass and stages each chunk's shells after a barrier, the loads' latency exposed per
// chunk. Chunks are 16 consecutive AOs of the list, a shell across a chunk's edge made in both (each its own rows),
// so every chunk is full but the last; a chunk stages only its rows of C, loaded into registers during the previous
// chunk's AO values and products and stored after them. The products and the epilogue are the first's, the sums in
// the same order. Used for sets of more than 192 functions (iso_gcorners_y256, iso_geval_y256); below, the first
// measured as fast.
// ------------------------------------------------------------------------------------------

#define KW_LIST 128
#define KW_ROWS 512
#define KW_PR   512
#define KW_PMAX 7                           // primitives of a shell kept in groupshared memory (more: read globally)
groupshared float4 kw_cs[KW_LIST];          // listed shell: centre, cutoff r^2
groupshared uint   kw_ci[KW_LIST];          // l | ncart << 4 | nprim << 8
groupshared uint   kw_po[KW_LIST];          // its first primitive in kw_pr, or 1 << 31 | its first in a.prims
groupshared uint   kw_la[KW_LIST + 1];      // its first AO in the list (kw_la[kw_n] = kw_a)
groupshared uint   kw_cf[KW_ROWS / 16 + 1]; // chunk c: the shell holding the list's AO 16 c (then kw_n)
groupshared uint   kw_rowao[KW_ROWS];       // the list's AOs: global index
groupshared float2 kw_pr[KW_PR];
groupshared uint   kw_n, kw_a, kw_p, kw_next;  // shells, AOs and primitives listed; where the next pass starts

// Lists the shells from s_begin on that reach the box, until KW_LIST shells or KW_ROWS AOs (all
// threads; barriers inside). kw_next: the first shell not taken (num_shells: all taken).
void kw_list(uint s_begin, float3 lo, float3 hi, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    if (tid == 0u) { kw_n = 0u; kw_a = 0u; kw_p = 0u; kw_next = a.num_shells; }
    GroupMemoryBarrierWithGroupSync();
    for (uint s0 = s_begin; s0 < a.num_shells; s0 += gs) {
        uint n0 = kw_n, a0 = kw_a, p0 = kw_p;
        uint s = min(s0 + tid, a.num_shells - 1u);        // clamped: see k_list
        IsoShell sh = a.shells[s];
        bool in = s0 + tid < a.num_shells && box_dist2(sh.c.xyz, lo, hi) <= sh.c.w;
        uint nc = sh.info.y, np = sh.info.z;
        uint npl = np <= uint(KW_PMAX) ? np : 0u;
        uint tot;
        // count < 512, AOs < 4096, primitives < 2048 per round of up to 256 shells
        uint pre = iso_group_prefix_sum(in ? (1u | (nc << 9) | (npl << 21)) : 0u, tid, gs, tot);
        uint at = n0 + (pre & 511u), la = a0 + ((pre >> 9) & 4095u), po = p0 + (pre >> 21);
        // taken: a prefix of the round's shells (at and la + nc grow with the shell)
        bool take = in && at < uint(KW_LIST) && la + nc <= uint(KW_ROWS);
        if (in && !take) InterlockedMin(kw_next, s0 + tid);
        if (take) {
            bool loc = npl != 0u && po + np <= uint(KW_PR);
            kw_cs[at] = sh.c;
            kw_ci[at] = sh.info.x | (nc << 4) | (np << 8);
            kw_po[at] = loc ? po : (0x80000000u | sh.info.w);
            kw_la[at] = la;
            for (uint k = 0u; k < nc; ++k) kw_rowao[la + k] = sh.info2.x + k;
            if (loc) for (uint q = 0u; q < np; ++q) kw_pr[po + q] = a.prims[sh.info.w + q];
            for (uint c = (la + 15u) >> 4; 16u * c < la + nc; ++c) kw_cf[c] = at;
            InterlockedMax(kw_n, at + 1u);
            InterlockedMax(kw_a, la + nc);
            if (loc) InterlockedMax(kw_p, po + np);
        }
        GroupMemoryBarrierWithGroupSync();
        if (kw_next < a.num_shells) break;                 // uniform: the list is full
    }
    if (tid == 0u) { kw_la[kw_n] = kw_a; kw_cf[(kw_a + 15u) >> 4] = kw_n; }
    GroupMemoryBarrierWithGroupSync();
}

// Listed shell i at p: d = p - centre, the radial part R (and R1 = dR/d(r^2) with GRAD), l, the AO count
void kw_shell<let GRAD : int>(uint i, float3 p, out float3 d, out float R, out float R1, out uint l, out uint nc) {
    IsoArgs a = MD_ARGS;
    float4 cs = kw_cs[i];
    uint ci = kw_ci[i];
    d = p - cs.xyz;
    float r2 = dot(d, d);
    l = ci & 15u;
    nc = (ci >> 4) & 15u;
    R = 0.0;
    R1 = 0.0;
    if (r2 > cs.w) return;
    uint np = (ci >> 8) & 255u, po = kw_po[i];
    bool glob = (po & 0x80000000u) != 0u;
    po &= 0x7FFFFFFFu;
    for (uint q = 0u; q < np; ++q) {
        float2 pr;     // one or the other read (llvmpipe reads both sides of a select: each index kept in range)
        if (glob) pr = a.prims[po + q];
        else      pr = kw_pr[min(po + q, uint(KW_PR) - 1u)];
        float t = pr.y * exp(-pr.x * r2);
        R += t;
        if (GRAD != 0) R1 -= pr.x * t;
    }
}

GC_MEM4(IsoGw128, 64, 128)
GC_MEM4(IsoGw192, 64, 192)
GC_MEM4(IsoGw256, 64, 256)

// The rows of a wide group: false when the group has nothing to do (uniform)
struct IsoWideRows {
    float3 p, glo, ghi;
    uint   pg, cg, ce, js, b, pt, slot0, mask, v0;
    bool   active, out_active;
};
bool iso_wide_rows<let VERT : int>(uint3 gid, uint tid, out IsoWideRows w) {
    IsoArgs a = MD_ARGS;
    w.pg = tid & 15u; w.cg = tid >> 4;
    w.ce = w.cg & 3u; w.js = w.cg >> 2;
    w.b = 0u; w.pt = 0u; w.slot0 = 0u; w.mask = 0u; w.v0 = 0u;
    w.active = false; w.out_active = false;
    w.p = float3(0.0); w.glo = float3(0.0); w.ghi = float3(0.0);
    if (VERT == 0) {
        uint total = a.count_in[0];
        uint e0 = node_index(gid) * uint(W_ROWS);
        if (e0 >= total) return false;                     // uniform per group
        w.active = e0 + 4u * w.pg + w.ce < total;
        uint ee = a.plist[min(e0 + 4u * w.pg + w.ce, total - 1u)];
        w.p = corner_point(a.bricks[ee >> 7], ee & 127u);
        w.out_active = tid < uint(W_ROWS) && e0 + tid < total;
        uint ent = a.plist[min(e0 + min(tid, uint(W_ROWS) - 1u), total - 1u)];
        w.b = ent >> 7;
        w.pt = ent & 127u;
        w.mask = a.bricks[w.b].w;
        w.slot0 = a.brick_slot[w.b];
    } else {
        w.v0 = node_index(gid) * 16u;
        if (w.v0 >= a.count) return false;                 // uniform per group
        if (tid == 0u) {
            uint done = 1u;
            for (uint k = 0u; k < 16u && w.v0 + k < a.count; ++k) done &= a.vdone[w.v0 + k] != 0u ? 1u : 0u;
            g_done = done;
        }
        GroupMemoryBarrierWithGroupSync();
        if (g_done != 0u) return false;                    // uniform: every vertex here has converged
        w.active = w.v0 + w.pg < a.count;
        w.out_active = tid < uint(W_ROWS) && w.v0 + (tid >> 2) < a.count;
        w.p = a.vpos[min(w.v0 + w.pg, a.count - 1u)].xyz;
    }
    {
        uint lanes = WaveGetLaneCount();
        uint wi = min(iso_wave_index(tid, lanes), 31u);
        float3 lo = WaveActiveMin(w.active ? w.p : float3(3.0e38));
        float3 hi = WaveActiveMax(w.active ? w.p : float3(-3.0e38));
        if (WaveIsFirstLane()) { g_lohi[2u * wi] = float4(lo, 0.0); g_lohi[2u * wi + 1u] = float4(hi, 0.0); }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) {
            uint nw = min(iso_num_waves(W_GROUP, lanes), 32u);
            float3 glo = g_lohi[0].xyz, ghi = g_lohi[1].xyz;
            for (uint k = 1u; k < nw; ++k) { glo = min(glo, g_lohi[2u * k].xyz); ghi = max(ghi, g_lohi[2u * k + 1u].xyz); }
            g_lohi[0] = float4(glo, 0.0);
            g_lohi[1] = float4(ghi, 0.0);
        }
        GroupMemoryBarrierWithGroupSync();
    }
    w.glo = g_lohi[0].xyz;
    w.ghi = g_lohi[1].xyz;
    return true;
}

// A shell's AO values at the thread's point into rows [rb, rb + nc) of the chunk (those in [0, 16)), row r at
// r * 64 + 4 pg + ce of the A chunk; VERT: channel ce (0 the value, 1 + ax the derivative along ax)
void iso_wide_ao<let VERT : int, Mem : IGcMem, let CLIP : int>(int rb, uint l, uint nc, float3 d, float R, float R1, uint pg, uint ce) {
    uint ax = ce > 0u ? ce - 1u : 0u;
    float m[15];
    cart_monomials(l, d, m);
    float Av = R, Bg = 0.0;
    float g[15];
    if (VERT != 0) {
        cart_monomials_d(l, d, ax, g);
        float da = ax == 0u ? d.x : (ax == 1u ? d.y : d.z);
        Av = ce == 0u ? R : 2.0 * da * R1;
        Bg = ce == 0u ? 0.0 : R;
    }
    [ForceUnroll] for (int n = 0; n < 15; ++n) {
        int r = rb + n;
        if (uint(n) < nc && (CLIP == 0 || (r >= 0 && r < 16))) {   // CLIP: a shell may straddle the chunk's edges
            float phi = m[n] * Av;
            if (VERT != 0) phi += g[n] * Bg;
            Mem::st4c(uint(r) * uint(W_ROWS / 4) + pg, ce, phi);
        }
    }
}

// The wide GEMM's first staging: k_list's (K_LIST shells listed at a time, chunks of whole shells of at most G_CK rows)
struct IsoWideStager<let F : int, Mem : IGcMem, let VERT : int> : IGemmStager {
    static const int kFull = 0;
    IsoArgs a;                                             // the kernel's (not read from MD_ARGS again per chunk)
    IsoWideRows w;
    uint s0, c, nchunk;
    bool listed;
    [mutating] void begin(uint m0, uint tid) { s0 = 0u; c = 0u; nchunk = 0u; listed = false; }
    [mutating] uint stage(uint m0, uint tid) {
        while (c >= nchunk) {                              // uniform
            if (listed) s0 += uint(K_LIST);
            if (s0 >= a.num_shells) return GC_END;
            k_list(s0, w.glo, w.ghi, tid, W_GROUP, G_CK);  // barriers inside
            listed = true;
            nchunk = k_nchunk;
            c = 0u;
        }
        k_load_chunk_gcw<F, Mem>(c, m0, tid, W_GROUP);
        uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
        for (uint i = i0 + w.js; i < i1; i += 4u) {
            float3 d;
            float R, R1 = 0.0;
            uint l, nc;
            if (VERT == 0) k_shell(i - i0, w.p, d, R, l, nc);
            else           k_shell_d(i - i0, w.p, d, R, R1, l, nc);
            iso_wide_ao<VERT, Mem, 0>(int(k_row[i]), l, nc, d, R, R1, w.pg, w.ce);
        }
        const uint nrows = k_nrows[c];
        c += 1u;
        return nrows;
    }
};

// The wide GEMM's second staging: the shells listed once a pass (kw_list), chunks of 16 consecutive AOs of the list (a shell across
// a chunk's edge made in both), the next chunk's rows of C loaded into registers while the current one is used
struct IsoWide2Stager<let F : int, Mem : IGcMem, let VERT : int> : IGemmStager {
    static const int kFull = 0;
    IsoArgs a;                                             // the kernel's (not read from MD_ARGS again per chunk)
    IsoWideRows w;
    uint s_begin, A, nsh, nchunk, c;
    bool done;
    float cr[F / 16];                                      // the next chunk's rows of C: e = tid + W_GROUP k
    [mutating] void load_rows(uint ln, uint f0, uint tid) {
        [ForceUnroll] for (int k = 0; k < F / 16; ++k) {
            uint e = tid + uint(W_GROUP * k), row = e / uint(F), col = e % uint(F);
            uint ao = kw_rowao[min(ln + row, max(A, 1u) - 1u)];
            cr[k] = (ln + row < A && f0 + col < a.num_fn) ? a.C[ao * a.fn_stride + min(f0 + col, a.fn_stride - 1u)] * a.ao_nrm[ao] : 0.0;
        }
    }
    [mutating] void begin(uint m0, uint tid) {
        s_begin = 0u; A = 0u; nsh = 0u; nchunk = 0u; c = 0u;
        done = a.num_shells == 0u;
    }
    [mutating] uint stage(uint m0, uint tid) {
        while (c >= nchunk) {                              // uniform
            if (done) return GC_END;
            kw_list(s_begin, w.glo, w.ghi, tid, W_GROUP);  // barriers inside
            s_begin = kw_next;
            done = s_begin >= a.num_shells;
            A = kw_a;
            nsh = kw_n;
            nchunk = (A + 15u) >> 4;
            c = 0u;
            if (nchunk > 0u) load_rows(0u, m0, tid);
        }
        uint la0 = 16u * c;
        uint nrows = min(16u, A - la0);
        // the prefetched rows of chunk c into B, then chunk c + 1's
        [ForceUnroll] for (int k = 0; k < F / 16; ++k) Mem::st4c(uint(GK * W_ROWS / 4) + (tid >> 2) + uint(W_GROUP / 4 * k), tid & 3u, cr[k]);
        load_rows(la0 + 16u, m0, tid);
        // the chunk's AO values: shells [i0, i1), each its rows within the chunk
        uint sn = kw_cf[c + 1u];
        uint i0 = kw_cf[c], i1 = (sn < nsh && kw_la[min(sn, uint(KW_LIST))] < la0 + 16u) ? sn + 1u : sn;
        for (uint i = i0 + w.js; i < i1; i += 4u) {
            float3 d;
            float R, R1;
            uint l, nc;
            kw_shell<VERT>(i, w.p, d, R, R1, l, nc);
            iso_wide_ao<VERT, Mem, 1>(int(kw_la[i]) - int(la0), l, nc, d, R, R1, w.pg, w.ce);
        }
        c += 1u;
        return nrows;
    }
};

// k_load_chunk with the rows of C into the wide core's B chunk
void k_load_chunk_gcw<let F : int, Mem : IGcMem>(uint c, uint f0, uint tid, uint gs) {
    IsoArgs a = MD_ARGS;
    uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
    uint i = min(i0 + tid, i1 - 1u);       // clamped: see k_list
    IsoShell sh = a.shells[min(k_sh[i], a.num_shells - 1u)];
    uint ao0 = sh.info2.x;
    if (tid < i1 - i0) {
        uint np = sh.info.z, po = k_poff[i];
        for (uint n = 0u; n < k_nc[i]; ++n) k_row_ao[k_row[i] + n] = ao0 + n;
        k_cs[tid] = sh.c;
        k_ci[tid] = sh.info.x | (sh.info.y << 4) | (np << 8) | (po << 16);
        for (uint q = 0u; q < np; ++q) k_pr[min(po + q, uint(K_PR) - 1u)] = a.prims[sh.info.w + q];
    }
    GroupMemoryBarrierWithGroupSync();
    uint nrows = k_nrows[c];
    for (uint e = tid; e < nrows * uint(F); e += gs) {
        uint row = e / uint(F), col = e % uint(F);
        uint ao = min(k_row_ao[min(row, uint(K_CK) - 1u)], a.num_ao - 1u);
        float v = a.C[ao * a.fn_stride + min(f0 + col, a.fn_stride - 1u)] * a.ao_nrm[ao];
        Mem::st4c(uint(GK * W_ROWS / 4) + (e >> 2), e & 3u, f0 + col < a.num_fn ? v : 0.0);
    }
    GroupMemoryBarrierWithGroupSync();
}

// The wide GEMM's epilogue (canonical terms): per row the partial sums of the thread's columns, reduced over the 16 column
// groups through the memory's first GK GP floats
struct IsoWideEpilogue<let F : int, Mem : IGcMem, let VERT : int> : IGemmEpilogue {
    IsoWideRows w;
    uint  tid;
    float res;
    [mutating] void take<let R : int, let M : int>(uint f0, uint ncols, uint tx, uint ty, float X[R][M]) {
        IsoArgs a = MD_ARGS;
        const int FC = M;
        uint stride = a.fn_stride;
        uint col0 = f0 + w.cg * uint(FC);
        for (uint f = 0u; f < (VERT == 0 ? a.num_fields : 1u); ++f) {
            uint fw = f;
            if (VERT != 0) fw = vertex_field(min(w.v0 + w.pg, a.count - 1u));   // the vertex's field
            uint wb = fw * 3u * stride;
            float part[4] = {0.0, 0.0, 0.0, 0.0};
            [ForceUnroll] for (int c = 0; c < FC; ++c) {
                // columns past the row of weights (a pass of 256 over a stride of 320, say): zero (their C is zero
                // too), and not read
                uint cc = col0 + uint(c);
                bool in = cc < stride;
                uint wi = wb + min(cc, stride - 1u);
                float wl = in ? a.wts[wi] : 0.0, wd = in ? a.wts[wi + stride] : 0.0;
                float wp = ((c & 1) == 0 && in) ? a.wts[wi + 2u * stride] : 0.0;
                if (VERT == 0) {
                    [ForceUnroll] for (int v = 0; v < 4; ++v) {
                        float x = X[v][c];
                        part[v] += x * (wl + wd * x);
                        if ((c & 1) == 0) part[v] += wp * x * X[v][c + 1];
                    }
                } else {
                    float x = X[0][c];
                    part[0] += x * (wl + wd * x);
                    [ForceUnroll] for (int k = 1; k < 4; ++k) part[k] += X[k][c] * (wl + 2.0 * wd * x);
                    if ((c & 1) == 0) {
                        float x1 = X[0][c + 1];
                        part[0] += wp * x * x1;
                        [ForceUnroll] for (int k = 1; k < 4; ++k) part[k] += wp * (x * X[k][c + 1] + x1 * X[k][c]);
                    }
                }
            }
            GroupMemoryBarrierWithGroupSync();             // the memory free
            [ForceUnroll] for (int v = 0; v < 4; ++v) Mem::st4c(w.cg * uint(W_ROWS / 4) + w.pg, uint(v), part[v]);
            GroupMemoryBarrierWithGroupSync();
            if (tid < uint(W_ROWS)) {
                float sum = 0.0;
                [ForceUnroll] for (int k = 0; k < 16; ++k) sum += Mem::ld(uint(k) * uint(W_ROWS) + tid);
                if (VERT == 0) {
                    if (w.out_active && ((field_bits(w.mask) >> (2u * f)) & 1u) != 0u) {
                        uint ci = (w.slot0 + field_slot(w.mask, f)) * 125u + w.pt;
                        a.corners[ci] = f0 == 0u ? sum : a.corners[ci] + sum;
                    }
                } else {
                    res += sum;
                }
            }
        }
    }
};

void gemm_wide<let F : int, Mem : IGcMem, let VERT : int, let V2 : int>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    IsoWideRows w;
    if (!iso_wide_rows<VERT>(gid, tid, w)) return;
    IsoWideEpilogue<F, Mem, VERT> epi;
    epi.w = w;
    epi.tid = tid;
    epi.res = 0.0;
    if (V2 != 0) {
        IsoWide2Stager<F, Mem, VERT> st;
        st.a = a;
        st.w = w;
        st.s_begin = 0u; st.A = 0u; st.nsh = 0u; st.nchunk = 0u; st.c = 0u; st.done = true;
        [ForceUnroll] for (int k = 0; k < F / 16; ++k) st.cr[k] = 0.0;
        gemm_core_staged<W_ROWS, F, 16, 16, Mem, IsoWide2Stager<F, Mem, VERT>, IsoWideEpilogue<F, Mem, VERT> >(st, epi, a.num_fn, tid, w.cg, w.pg);
    } else {
        IsoWideStager<F, Mem, VERT> st;
        st.a = a;
        st.w = w;
        st.s0 = 0u; st.c = 0u; st.nchunk = 0u; st.listed = false;
        gemm_core_staged<W_ROWS, F, 16, 16, Mem, IsoWideStager<F, Mem, VERT>, IsoWideEpilogue<F, Mem, VERT> >(st, epi, a.num_fn, tid, w.cg, w.pg);
    }
    if (VERT != 0) {
        GroupMemoryBarrierWithGroupSync();
        if (tid < uint(W_ROWS)) Mem::st(tid, epi.res);
        GroupMemoryBarrierWithGroupSync();
        if (tid < 16u && w.v0 + tid < a.count) a.vval[w.v0 + tid] = float4(Mem::ld(4u * tid + 1u), Mem::ld(4u * tid + 2u), Mem::ld(4u * tid + 3u), Mem::ld(4u * tid));
    }
}

[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_gcorners_w128(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<128, IsoGw128, 0, 0>(gid, tid); }
[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_gcorners_w192(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<192, IsoGw192, 0, 0>(gid, tid); }
[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_gcorners_y256(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<256, IsoGw256, 0, 1>(gid, tid); }
[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_geval_w128(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<128, IsoGw128, 1, 0>(gid, tid); }
[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_geval_w192(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<192, IsoGw192, 1, 0>(gid, tid); }
[shader("compute")] [numthreads(W_GROUP, 1, 1)] void iso_geval_y256(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_wide<256, IsoGw256, 1, 1>(gid, tid); }

// ------------------------------------------------------------------------------------------
// Densities as phi^T D phi (dform): the local list, the shells that reach a box in shell order with
// their first local AO (the D form's corners, vertices and octree work over the A local AOs)
// ------------------------------------------------------------------------------------------

#define D_LS 2048                        // shells in the local list (more: corners and vertices split the group, the octree keeps the node)
#define D_CH 16                          // AOs per chunk
groupshared uint  d_ent[D_LS];           // the local list: shells reaching the box, in order: shell | first local AO << 16
groupshared uint  d_nsh;                 // shells listed
groupshared uint  d_atot;                // A, the local AOs
groupshared uint  d_ovf;                 // the list overflowed
groupshared uint  d_any;
groupshared uint  d_rs[16];              // the D form's ranges of rows still to do: first | end << 16
groupshared uint  d_rn;

uint d_sh(uint i)  { return d_ent[i] & 0xFFFFu; }
uint d_off(uint i) { return i < d_nsh ? d_ent[i] >> 16 : d_atot; }

// The local shell holding local AO la (first AOs ascending)
uint d_find(uint la) {
    uint lo = 0u, hi = d_nsh;                           // d_off(lo) <= la < d_off(hi)
    while (hi - lo > 1u) {
        uint mid = (lo + hi) >> 1;
        if ((d_ent[mid] >> 16) <= la) lo = mid; else hi = mid;
    }
    return lo;
}

// The local list for the box, the shells reaching it within r^2 cutoff * reach (all threads of a
// group of G_GROUP, or gs; barriers inside). An overflow sets d_ovf and keeps the first D_LS.
void d_list(float3 lo, float3 hi, uint tid, float reach) { d_list_gs(lo, hi, tid, reach, uint(G_GROUP)); }

void d_list_gs(float3 lo, float3 hi, uint tid, float reach, uint gs) {
    IsoArgs a = MD_ARGS;
    uint nsh = 0u, nao = 0u;
    for (uint s0 = 0u; s0 < a.num_shells; s0 += gs) {
        uint s = min(s0 + tid, a.num_shells - 1u);
        IsoShell sh = a.shells[s];
        bool in = s0 + tid < a.num_shells && box_dist2(sh.c.xyz, lo, hi) <= sh.c.w * reach;
        uint tot;
        uint pre = iso_group_prefix_sum(in ? (1u | (sh.info.y << 16)) : 0u, tid, gs, tot);
        if (in) {
            uint at = nsh + (pre & 0xFFFFu);
            if (at < uint(D_LS)) d_ent[at] = s | ((nao + (pre >> 16)) << 16);
        }
        nsh += tot & 0xFFFFu;
        nao += tot >> 16;
    }
    if (tid == 0u) {
        d_nsh = min(nsh, uint(D_LS));
        d_ovf = nsh > uint(D_LS) ? 1u : 0u;
        d_atot = nsh > uint(D_LS) ? (d_ent[D_LS - 1] >> 16) + a.shells[d_ent[D_LS - 1] & 0xFFFFu].info.y : nao;
    }
    GroupMemoryBarrierWithGroupSync();
}

// ------------------------------------------------------------------------------------------
// The D form (iso_dcorners, iso_deval), on the GEMM core (gemm_core_tri: staged, triangular by chunks): 128 rows per
// group as in the tiled GEMM (corners: consecutive points
// of the list; vertices: 32 vertices x 4 channels), the local list in chunks of 16 AOs, tiles of 16 NW
// columns (NW = 8 waves of 32 threads), per tile T = Phi D' and rho += Phi . T, with
//  - D' upper triangular by chunks of 16 AOs (2 D before the column's chunk, D on its diagonal block,
//    0 after), so the T of a column chunk is complete at its own row chunk, and Phi . T is taken then
//    (the epilogue's diag) from the AO values still in the A chunk: a tile's AO values are not made a second time;
//  - the column chunks as waves (cg = tid / 32): a wave skips the row chunks after its own (zero in
//    D') and a column chunk past the list's end;
//  - each chunk's first shell from a table made once per group (no searches per chunk), and the next
//    chunk's shells staged while the current one is evaluated (two buffers): two barriers per chunk;
//  - the tile's D' loaded into registers before the AO values are made, stored after;
//  - Phi . T summed over the tiles in registers and reduced once per matrix.
// Thread (pg, cg) makes the AO values of row 4 pg + cg % 4 (corners: that point of the group;
// vertices: vertex pg, channel cg % 4), for the chunk's shells j with j % (NW / 4) == cg / 4, and
// threads 0 .. 127 write the results of rows 0 .. 127. The cost per row is about A^2 / 2 for A local
// AOs, whatever the density's rank. Thread (pg, cg) = (tid & 31, tid >> 5) is the core's (ty, tx): rows 4 pg ..
// 4 pg + 3, the 16 columns of wave cg's chunk of a tile of 128 (R = 4, M = 16). A chunk: A [kk][r] the rows' AO
// values, B [kk][128] the chunk's rows of D'.
// ------------------------------------------------------------------------------------------

// (sizes kept so every kernel stays within 32 KB of groupshared memory, the Apple GPUs' limit)
#define DQ_PMAX 8                        // primitives per staged shell in groupshared memory (more: global memory)
#define DQ_CAP  512                      // chunks in the table: at most 8192 local AOs (more: as for D_LS)

groupshared uint   dq_tab[DQ_CAP / 2];   // per chunk q, 16 bits: the local shell holding local AO 16 q
groupshared float4 dq_cs[2 * D_CH];      // [buffer][staged shell]: centre, r^2 cutoff
groupshared uint   dq_ci[2 * D_CH];      // l | ncart << 4 | nprim << 8 | 1 << 31 when its primitives are read globally
groupshared int    dq_rb[2 * D_CH];      // chunk row of its first AO (negative: it began in the chunk before)
groupshared uint   dq_pg[2 * D_CH];      // its first primitive (global)
groupshared float2 dq_pr[2 * D_CH * DQ_PMAX];
groupshared uint   dq_rao[2 * D_CH];     // [buffer][chunk row]: global AO (INVALID past A)
groupshared uint   dq_jao[2 * 128];      // [tile parity][column]: global AO (INVALID past A)

uint dq_first(uint q) { return (dq_tab[q >> 1] >> ((q & 1u) * 16u)) & 0xFFFFu; }

// Chunk q's shells: the first local shell and how many (at most 16)
uint2 dq_span(uint q, uint nq) {
    uint i0 = dq_first(q);
    uint last = d_nsh - 1u;
    if (q + 1u < nq) {
        uint i1 = dq_first(q + 1u);
        last = (d_ent[i1] >> 16) == (q + 1u) * uint(D_CH) ? i1 - 1u : i1;
    }
    return uint2(i0, last - i0 + 1u);
}

// The chunk table from the local list (all gs threads; barriers inside)
void dq_table(uint tid, uint gs, uint nq) {
    for (uint k = tid; k < (nq + 1u) / 2u; k += gs) dq_tab[k] = 0u;
    GroupMemoryBarrierWithGroupSync();
    for (uint i = tid; i < d_nsh; i += gs) {
        uint off = d_ent[i] >> 16, end = d_off(i + 1u);
        uint q = (off + uint(D_CH) - 1u) / uint(D_CH);    // a shell holds at most 15 AOs: one multiple of 16 at most
        if (q * uint(D_CH) < end) InterlockedOr(dq_tab[q >> 1], i << ((q & 1u) * 16u));
    }
    GroupMemoryBarrierWithGroupSync();
}

// Stages chunk q's shells into buffer buf (threads 0 .. 31; no barriers)
void dq_stage(uint q, uint nq, uint buf, uint tid) {
    IsoArgs a = MD_ARGS;
    uint la0 = q * uint(D_CH);
    if (tid < uint(D_CH)) {
        uint2 sp = dq_span(q, nq);
        if (tid < sp.y) {
            uint li = sp.x + tid;
            IsoShell sh = a.shells[d_sh(li)];
            int rb = int(d_ent[li] >> 16) - int(la0);
            uint np = sh.info.z, nc = sh.info.y;
            bool fits = np <= uint(DQ_PMAX);
            uint s = buf * uint(D_CH) + tid;
            dq_cs[s] = sh.c;
            dq_ci[s] = sh.info.x | (nc << 4) | (np << 8) | (fits ? 0u : 0x80000000u);
            dq_rb[s] = rb;
            dq_pg[s] = sh.info.w;
            for (uint k = 0u; k < nc; ++k) {
                int row = rb + int(k);
                if (row >= 0 && row < D_CH) dq_rao[buf * uint(D_CH) + uint(row)] = sh.info2.x + k;
            }
            if (fits) for (uint k = 0u; k < np; ++k) dq_pr[s * uint(DQ_PMAX) + k] = a.prims[sh.info.w + k];
        }
    } else if (tid < 2u * uint(D_CH)) {                    // the rows past the list's end
        uint r = tid - uint(D_CH);
        if (la0 + r >= d_atot) dq_rao[buf * uint(D_CH) + r] = INVALID;
    }
}

// Tile J's nc columns into dq_jao[par] (threads nc .. 2 nc - 1; no barriers)
void dq_cols(uint J, uint par, uint tid, uint nc) {
    if (tid < nc || tid >= 2u * nc) return;
    IsoArgs a = MD_ARGS;
    uint t = tid - nc;
    uint la = J * nc + t;
    uint ao = INVALID;
    if (la < d_atot) {
        uint i = dq_first(la / uint(D_CH));
        while (i + 1u < d_nsh && (d_ent[i + 1u] >> 16) <= la) ++i;
        ao = a.shells[d_sh(i)].info2.x + (la - (d_ent[i] >> 16));
    }
    dq_jao[par * 128u + t] = ao;
}

#define NW 8                             // waves per group

GC_MEM4(IsoGd, 128, 128)

// This row's values of the chunk in buffer buf (or channel ch of their gradient) into the core's A chunk (row r at
// float4 r * 32 + pg, component ce): the chunk's shells j0, j0 + js, ..., zero where beyond their cutoff; the rows
// past the list's end zero when j0 == 0. The chunk's shells cover its rows below A once each.
void dq_phi<Mem : IGcMem>(uint buf, uint ns, uint j0, uint js, uint nvalid, float3 p, uint pg, uint ce, bool deriv, uint ch) {
    IsoArgs a = MD_ARGS;
    if (j0 == 0u) for (uint k = nvalid; k < uint(D_CH); ++k) Mem::st4c(k * 32u + pg, ce, 0.0);
    uint ax = ch > 0u ? ch - 1u : 0u;
    for (uint j = j0; j < ns; j += js) {
        uint s = buf * uint(D_CH) + j;
        float4 cs = dq_cs[s];
        uint ci = dq_ci[s];
        int rb = dq_rb[s];
        uint l = ci & 15u, nc = (ci >> 4) & 15u, np = (ci >> 8) & 255u;
        float3 d = p - cs.xyz;
        float r2 = dot(d, d);
        if (r2 > cs.w) {                                   // beyond the cutoff: its rows zero
            [ForceUnroll] for (int n = 0; n < 15; ++n) {
                int row = rb + n;
                if (uint(n) < nc && row >= 0 && row < D_CH) Mem::st4c(uint(row) * 32u + pg, ce, 0.0);
            }
            continue;
        }
        float m[15];
        float g[15];
        float A = 0.0, Bg = 0.0;
        {
            float R = 0.0, R1 = 0.0;
            if ((ci >> 31) != 0u) {                        // uniform: the same shell for every thread
                uint po = dq_pg[s];
                for (uint q = 0u; q < np; ++q) {
                    float2 pr = a.prims[po + q];
                    float t = pr.y * exp(-pr.x * r2);
                    R += t;
                    if (deriv) R1 -= pr.x * t;
                }
            } else {
                uint po = s * uint(DQ_PMAX);
                for (uint q = 0u; q < np; ++q) {
                    float2 pr = dq_pr[po + q];
                    float t = pr.y * exp(-pr.x * r2);
                    R += t;
                    if (deriv) R1 -= pr.x * t;
                }
            }
            A = R;
            if (deriv) {
                float da = ax == 0u ? d.x : (ax == 1u ? d.y : d.z);
                A = ch == 0u ? R : 2.0 * da * R1;
                Bg = ch == 0u ? 0.0 : R;
            }
        }
        cart_monomials(l, d, m);
        if (deriv) cart_monomials_d(l, d, ax, g);
        [ForceUnroll] for (int n = 0; n < 15; ++n) {
            int row = rb + n;
            if (uint(n) < nc && row >= 0 && row < D_CH) {
                float phi = m[n] * A;
                if (deriv) phi += g[n] * Bg;
                Mem::st4c(uint(row) * 32u + pg, ce, phi);
            }
        }
    }
}

// The chunks of a tile J (m0 = 128 J): q = 0 .. qe - 1, qe = min(NW (J + 1), nq) (the chunks after the tile's last
// column are zero in D'). Each: this thread's D' elements loaded first (rows tid / 128 + 2 i, column tid % 128), the
// next chunk's shells (or the next tile's first and its columns) staged into the other buffers, the AO values made,
// the D' rows stored
struct IsoDformStager<let VERT : int, Mem : IGcMem> : IGemmStager {
    static const int kFull = 0;
    float* D;
    uint   A, nq;
    uint   dm_nt;                                          // dform_at's arguments, the kernel's (not read from MD_ARGS
    uint*  dm_tao;                                         // again per chunk)
    uint*  dm_map;
    float* ao_nrm;
    float3 p;
    uint   pg, cg, ce, col;
    uint   J, q, qe, jb, ca;
    [mutating] void begin(uint m0, uint tid) {
        J  = m0 / (16u * uint(NW));
        q  = 0u;
        qe = min(uint(NW) * (J + 1u), nq);
        jb = (J & 1u) * 128u;
        if (J == 0u) {                                     // the first tile's first chunk and columns (later: staged ahead)
            dq_stage(0u, nq, 0u, tid);
            dq_cols(0u, 0u, tid, 16u * uint(NW));
        }
    }
    [mutating] uint stage(uint m0, uint tid) {
        if (q >= qe) return GC_END;
        const uint TC = 16u * uint(NW);
        if (q == 0u) ca = dq_jao[jb + col];                // after the barrier: the tile's columns staged
        uint buf = q & 1u;
        uint ccol = uint(NW) * J + col / uint(D_CH);
        float fac = q < ccol ? 2.0 : (q == ccol ? 1.0 : 0.0);
        float dv[8];
        [ForceUnroll] for (int i = 0; i < 8; ++i) {
            uint ra = dq_rao[buf * uint(D_CH) + tid / TC + 2u * uint(i)];
            bool ok = fac != 0.0 && ra != INVALID && ca != INVALID;
            float v = ok ? dform_at(D, dm_tao, dm_map, dm_nt, ao_nrm, ra, ca) : 0.0;
            dv[i] = ok ? fac * v : 0.0;
        }
        if (q + 1u < qe) {
            dq_stage(q + 1u, nq, buf ^ 1u, tid);
        } else if ((J + 1u) * TC < A) {
            dq_stage(0u, nq, 0u, tid);
            dq_cols(J + 1u, (J + 1u) & 1u, tid, TC);
        }
        uint la0 = q * uint(D_CH);
        dq_phi<Mem>(buf, dq_span(q, nq).y, cg >> 2, uint(NW) / 4u, min(A - la0, uint(D_CH)), p, pg, ce, VERT != 0, ce);
        [ForceUnroll] for (int i = 0; i < 8; ++i) Mem::st4c(uint(GK * 32) + (tid >> 2) + 64u * uint(i), tid & 3u, dv[i]);
        q += 1u;
        return uint(GK);
    }
};

// Phi . T at the diagonal chunk: the rows' AO values of the chunk (read again from the A chunk) times the column chunk
struct IsoDformEpilogue<let VERT : int, Mem : IGcMem> : IGemmTriEpilogue {
    uint  pg;
    float part[4];
    [mutating] void take<let R : int, let M : int>(uint m0, uint ncols, uint tx, uint ty, float X[R][M]) {}
    [mutating] void diag<let R : int, let M : int>(uint k0, uint tx, uint ty, float X[R][M]) {
        [ForceUnroll] for (int k = 0; k < GK; ++k) {
            float4 ph = Mem::ld4(uint(k) * 32u + pg);      // rows 4 pg .. 4 pg + 3 of row k of the A chunk
            if (VERT == 0) {
                [ForceUnroll] for (int v = 0; v < 4; ++v) part[v] += ph[v] * X[v][k];
            } else {
                float t0 = X[0][k];
                part[0] += ph.x * t0;
                [ForceUnroll] for (int c = 1; c < 4; ++c) part[c] += ph[c] * t0 + ph.x * X[c][k];
            }
        }
    }
};

void gemm_dform<let VERT : int>(uint3 gid, uint tid) {
    IsoArgs a = MD_ARGS;
    const uint GS = 32u * uint(NW);                        // threads
    const uint TC = 16u * uint(NW);                        // columns per tile
    const uint SP = uint(NW) / 4u;                         // threads per row making AO values
    uint pg = tid & 31u, cg = tid >> 5;
    uint ce = cg & 3u;                                     // the row this thread's AO values are for: 4 pg + ce
    uint b = 0u, pt = 0u, slot0 = 0u, mask = 0u, v0 = 0u, vf = 0u;
    bool active, out_active;
    float3 p;
    if (VERT == 0) {
        uint total = a.count_in[0];
        uint e0 = node_index(gid) * uint(G_GROUP);
        if (e0 >= total) return;                           // uniform per group
        active = e0 + 4u * pg + ce < total;
        uint ee = a.plist[min(e0 + 4u * pg + ce, total - 1u)];
        p = corner_point(a.bricks[ee >> 7], ee & 127u);
        out_active = tid < uint(G_GROUP) && e0 + tid < total;
        uint ent = a.plist[min(e0 + min(tid, uint(G_GROUP) - 1u), total - 1u)];
        b = ent >> 7;
        pt = ent & 127u;
        mask = a.bricks[b].w;
        slot0 = a.brick_slot[b];
    } else {
        v0 = node_index(gid) * 32u;
        if (v0 >= a.count) return;                         // uniform per group
        if (tid == 0u) {
            uint done = 1u;
            for (uint k = 0u; k < 32u && v0 + k < a.count; ++k) done &= a.vdone[v0 + k] != 0u ? 1u : 0u;
            g_done = done;
        }
        GroupMemoryBarrierWithGroupSync();
        if (g_done != 0u) return;                          // uniform: every vertex here has converged
        active = v0 + pg < a.count;
        p = a.vpos[min(v0 + pg, a.count - 1u)].xyz;
        uint vo = min(tid, uint(G_GROUP) - 1u) >> 2;
        out_active = tid < uint(G_GROUP) && v0 + vo < a.count;
        vf = vertex_field(min(v0 + vo, a.count - 1u));
    }
    // The group's rows in ranges, each with the local list of its own rows' box: one range, all rows,
    // unless that list is too long (D_LS shells, DQ_CAP chunks), when the range is halved until it is
    // not. Rows are consecutive points of the list (corners) or vertices; a group that straddles a
    // jump of the Morton order has a box across the molecule, its halves (mostly) do not. Every row
    // is evaluated in every range's pass; a pass writes only its range's rows.
    const uint units = VERT == 0 ? uint(G_GROUP) : 32u;    // corners: rows; vertices: 32 vertices of 4 rows
    const uint unit = VERT == 0 ? 4u * pg + ce : pg;       // the unit this thread makes AO values for
    const uint ounit = VERT == 0 ? tid : tid >> 2;         // the unit this thread writes (tid < G_GROUP)
    if (tid == 0u) { d_rs[0] = units << 16; d_rn = 1u; }
    GroupMemoryBarrierWithGroupSync();
    for (;;) {
        uint rn = d_rn;                                    // uniform (after a barrier)
        if (rn == 0u) break;
        uint rr = d_rs[rn - 1u];
        uint rs = rr & 0xFFFFu, re = rr >> 16;
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) d_rn = rn - 1u;
        bool in = active && unit >= rs && unit < re;
        bool out_in = out_active && ounit >= rs && ounit < re;
        {
            uint lanes = WaveGetLaneCount();
            uint w = min(iso_wave_index(tid, lanes), 31u);
            float3 lo = WaveActiveMin(in ? p : float3(3.0e38));
            float3 hi = WaveActiveMax(in ? p : float3(-3.0e38));
            if (WaveIsFirstLane()) { g_lohi[2u * w] = float4(lo, 0.0); g_lohi[2u * w + 1u] = float4(hi, 0.0); }
            GroupMemoryBarrierWithGroupSync();
            if (tid == 0u) {
                uint nw = min(iso_num_waves(GS, lanes), 32u);
                float3 glo = g_lohi[0].xyz, ghi = g_lohi[1].xyz;
                for (uint k = 1u; k < nw; ++k) { glo = min(glo, g_lohi[2u * k].xyz); ghi = max(ghi, g_lohi[2u * k + 1u].xyz); }
                g_lohi[0] = float4(glo, 0.0);
                g_lohi[1] = float4(ghi, 0.0);
            }
            GroupMemoryBarrierWithGroupSync();
        }
        float3 glo = g_lohi[0].xyz, ghi = g_lohi[1].xyz;
        d_list_gs(glo, ghi, tid, 1.0, GS);
        uint A = d_atot;
        uint nq = (A + uint(D_CH) - 1u) / uint(D_CH);
        if (d_ovf != 0u || nq > uint(DQ_CAP)) {            // uniform: too many shells reach the range's box
            if (tid == 0u) {
                if (re - rs > 1u && d_rn + 2u <= 16u) {
                    uint mid = (rs + re) >> 1;
                    d_rs[d_rn] = mid | (re << 16);
                    d_rs[d_rn + 1u] = rs | (mid << 16);
                    d_rn += 2u;
                } else {
                    InterlockedAdd(a.counters[7], 1u);     // one point that many shells reach: the values would be wrong
                }
            }
            GroupMemoryBarrierWithGroupSync();
            continue;
        }
        if (A == 0u) {                                     // uniform: nothing reaches the box
            if (VERT == 0) {
                if (out_in) {
                    for (uint f = 0u; f < a.num_fields; ++f)
                        if (((field_bits(mask) >> (2u * f)) & 1u) != 0u) a.corners[(slot0 + field_slot(mask, f)) * 125u + pt] = 0.0;
                }
            } else if (tid < 32u && v0 + tid < a.count && tid >= rs && tid < re) {
                a.vval[v0 + tid] = float4(0.0);
            }
            GroupMemoryBarrierWithGroupSync();
            continue;
        }
        dq_table(tid, GS, nq);
        uint nj = (A + TC - 1u) / TC;
        uint n = a.num_ao;
        float res = 0.0;                                   // vertices: row tid's channel of its field

        for (uint u = 0u; u < a.num_dmat; ++u) {
            // skipped where no row of the range needs this matrix
            bool need = false;
            if (out_in) {
                if (VERT == 0) {
                    for (uint f = 0u; f < a.num_fields; ++f) need = need || (((field_bits(mask) >> (2u * f)) & 1u) != 0u && a.fdmat[f] == u);
                } else {
                    need = a.fdmat[vf] == u;
                }
            }
            GroupMemoryBarrierWithGroupSync();
            if (tid == 0u) d_any = 0u;
            GroupMemoryBarrierWithGroupSync();
            if (need) d_any = 1u;
            GroupMemoryBarrierWithGroupSync();
            if (d_any == 0u) continue;
            float* D = a.dmat + u * a.dm_pool;
            float part[4] = {0.0, 0.0, 0.0, 0.0};          // rows 4 pg + v: Phi . T over this wave's column chunks

            IsoDformStager<VERT, IsoGd> st;
            st.D = D; st.A = A; st.nq = nq; st.p = p; st.pg = pg; st.cg = cg; st.ce = ce; st.col = tid % TC;
            st.dm_nt = a.dm_nt; st.dm_tao = a.dm_tao; st.dm_map = a.dm_map; st.ao_nrm = a.ao_nrm;
            st.J = 0u; st.q = 0u; st.qe = 0u; st.jb = 0u; st.ca = INVALID;
            IsoDformEpilogue<VERT, IsoGd> epi;
            epi.pg = pg;
            [ForceUnroll] for (int v = 0; v < 4; ++v) epi.part[v] = 0.0;
            gemm_core_tri<128, 128, 8, 32, IsoGd, IsoDformStager<VERT, IsoGd>, IsoDformEpilogue<VERT, IsoGd> >(st, epi, A, tid, cg, pg);
            [ForceUnroll] for (int v = 0; v < 4; ++v) part[v] = epi.part[v];
            // rows 4 pg + v summed over the waves; row tid's sum
            GroupMemoryBarrierWithGroupSync();             // the memory free
            [ForceUnroll] for (int v = 0; v < 4; ++v) IsoGd::st4c(cg * 32u + pg, uint(v), part[v]);
            GroupMemoryBarrierWithGroupSync();
            float r = 0.0;
            if (tid < uint(G_GROUP)) {
                [ForceUnroll] for (int c = 0; c < NW; ++c) r += IsoGd::ld((uint(c) * 32u + (tid >> 2)) * 4u + (tid & 3u));
            }
            if (VERT == 0) {
                if (out_in) {
                    for (uint f = 0u; f < a.num_fields; ++f) {
                        if (((field_bits(mask) >> (2u * f)) & 1u) == 0u || a.fdmat[f] != u) continue;
                        a.corners[(slot0 + field_slot(mask, f)) * 125u + pt] = r;
                    }
                }
            } else if (a.fdmat[vf] == u) {
                res = r;
            }
        }
        if (VERT != 0) {
            GroupMemoryBarrierWithGroupSync();
            if (tid < uint(G_GROUP)) IsoGd::st(tid, res);
            GroupMemoryBarrierWithGroupSync();
            if (tid < 32u && v0 + tid < a.count && tid >= rs && tid < re) {
                float4 t = float4(IsoGd::ld(4u * tid), IsoGd::ld(4u * tid + 1u), IsoGd::ld(4u * tid + 2u), IsoGd::ld(4u * tid + 3u));
                a.vval[v0 + tid] = float4(t.y, t.z, t.w, t.x);
            }
        }
        GroupMemoryBarrierWithGroupSync();
    }
}

[shader("compute")] [numthreads(256, 1, 1)] void iso_dcorners(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_dform<0>(gid, tid); }
[shader("compute")] [numthreads(256, 1, 1)] void iso_deval(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { gemm_dform<1>(gid, tid); }
#undef NW

#define OD_BLK 256                       // local AOs per block (two per thread)
groupshared float4 od_q[OD_BLK];         // a block of local AOs: value at the centre, deviation bound, magnitude bound, AO
groupshared float  od_sum[8 * 32];       // per wave: the node's sums

// Local AO la over the box: its raw value at the centre, deviation and magnitude bounds, global AO
float4 od_ao(uint la, float3 lo, float3 hi, float3 cen) {
    uint ls = d_find(la);
    OShell o = octree_shell_setup(d_sh(ls), lo, hi, cen);
    float mid, rad, pc, dev;
    uint ao;
    octree_ao(o, la - d_off(ls), mid, rad, pc, dev, ao, true);
    return float4(pc, dev, abs(mid) + rad, asfloat(ao));
}

[shader("compute")]
[numthreads(G_GROUP, 1, 1)]
void iso_octree_dform(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint ni = node_index(gid);
    if (ni >= min(a.count_in[0], a.cap)) return;          // uniform per group
    uint4 node = a.nodes_in[ni];
    float3 lo, hi, cen;
    octree_box(node, lo, hi, cen);
    if (tid == 0u) o_mask = 0u;
    d_list(lo, hi, tid, 1.0 + 16.0 * U);
    bool ovf = d_ovf != 0u;
    uint A = d_atot;
    uint nb = (A + uint(OD_BLK) - 1u) / uint(OD_BLK);
    uint n = a.num_ao;
    uint lanes = WaveGetLaneCount();
    uint w = min(iso_wave_index(tid, lanes), 31u);
    uint nw = min(iso_num_waves(G_GROUP, lanes), 32u);

    for (uint u = 0u; u < a.num_dmat; ++u) {
        bool any = false;
        for (uint f = 0u; f < a.num_fields; ++f) any = any || (((node.w >> (2u * f)) & 3u) != 0u && a.fdmat[f] == u);
        if (!any) continue;                                // uniform
        float sm[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // rho_c, sum |pc| ga, sum |g| dev, sum dev ga, sum dev h, m^T |D| m, sum m rs
        if (!ovf) {
            float* D = a.dmat + u * a.dm_pool;
            for (uint pb = 0u; pb < nb; ++pb) {
                uint la0 = pb * uint(OD_BLK) + tid, la1 = la0 + uint(G_GROUP);
                float4 m0 = la0 < A ? od_ao(la0, lo, hi, cen) : float4(0.0);
                float4 m1 = la1 < A ? od_ao(la1, lo, hi, cen) : float4(0.0);
                uint a0 = asuint(m0.w), a1 = asuint(m1.w);
                float g0 = 0.0, ga0 = 0.0, h0 = 0.0, k0 = 0.0, g1 = 0.0, ga1 = 0.0, h1 = 0.0, k1 = 0.0;
                for (uint qb = 0u; qb < nb; ++qb) {
                    GroupMemoryBarrierWithGroupSync();
                    if (qb == pb) {
                        od_q[tid] = m0;
                        od_q[tid + uint(G_GROUP)] = m1;
                    } else {
                        uint lb = qb * uint(OD_BLK) + tid;
                        od_q[tid] = lb < A ? od_ao(lb, lo, hi, cen) : float4(0.0);
                        od_q[tid + uint(G_GROUP)] = lb + uint(G_GROUP) < A ? od_ao(lb + uint(G_GROUP), lo, hi, cen) : float4(0.0);
                    }
                    GroupMemoryBarrierWithGroupSync();
                    uint qn = min(uint(OD_BLK), A - qb * uint(OD_BLK));
                    // D symmetric: row nu, the threads' columns mu. Four rows' loads issued before any is
                    // used (one load per 4 multiply-adds: the loop waits on memory otherwise). Every lane
                    // loads: a thread past A has AO 0 and its sums go unused. Summed in the order of v.
                    uint v = 0u;
                    for (; v + 4u <= qn; v += 4u) {
                        float4 q[4];
                        float d0[4], d1[4];
                        [ForceUnroll] for (int j = 0; j < 4; ++j) {
                            q[j] = od_q[v + uint(j)];
                            uint qa = asuint(q[j].w);
                            d0[j] = dform_at(D, a.dm_tao, a.dm_map, a.dm_nt, a.ao_nrm, qa, a0);
                            d1[j] = dform_at(D, a.dm_tao, a.dm_map, a.dm_nt, a.ao_nrm, qa, a1);
                        }
                        [ForceUnroll] for (int j = 0; j < 4; ++j) {
                            float apc = abs(q[j].x), e0 = abs(d0[j]), e1 = abs(d1[j]);
                            g0 += d0[j] * q[j].x; ga0 += e0 * apc; h0 += e0 * q[j].y; k0 += e0 * q[j].z;
                            g1 += d1[j] * q[j].x; ga1 += e1 * apc; h1 += e1 * q[j].y; k1 += e1 * q[j].z;
                        }
                    }
                    for (; v < qn; ++v) {
                        float4 q = od_q[v];
                        uint qa = asuint(q.w);
                        float apc = abs(q.x);
                        float d = dform_at(D, a.dm_tao, a.dm_map, a.dm_nt, a.ao_nrm, qa, a0), e = abs(d);
                        g0 += d * q.x; ga0 += e * apc; h0 += e * q.y; k0 += e * q.z;
                        d = dform_at(D, a.dm_tao, a.dm_map, a.dm_nt, a.ao_nrm, qa, a1); e = abs(d);
                        g1 += d * q.x; ga1 += e * apc; h1 += e * q.y; k1 += e * q.z;
                    }
                }
                if (la0 < A) {
                    sm[0] += m0.x * g0; sm[1] += abs(m0.x) * ga0; sm[2] += abs(g0) * m0.y; sm[3] += m0.y * ga0;
                    sm[4] += m0.y * h0; sm[5] += m0.z * k0; sm[6] += m0.z * a.drs[u * n + a0];
                }
                if (la1 < A) {
                    sm[0] += m1.x * g1; sm[1] += abs(m1.x) * ga1; sm[2] += abs(g1) * m1.y; sm[3] += m1.y * ga1;
                    sm[4] += m1.y * h1; sm[5] += m1.z * k1; sm[6] += m1.z * a.drs[u * n + a1];
                }
            }
        }
        // the sums over the group (in any order: the allowances below cover it)
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll] for (int k = 0; k < 7; ++k) {
            float t = WaveActiveSum(sm[k]);
            if (WaveIsFirstLane()) od_sum[w * 8u + uint(k)] = t;
        }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) {
            float2 rho;
            if (ovf) {
                rho = float2(-3.0e38, 3.0e38);
            } else {
                float t[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                for (uint k = 0u; k < nw; ++k) [ForceUnroll] for (int c = 0; c < 7; ++c) t[c] += od_sum[k * 8u + uint(c)];
                float4 pr = a.dprop[u];
                float gam = float(A + 8u) * 2.0 * U;
                float tail = pr.x + pr.y * t[6] * (1.0 + gam) + pr.w;
                float W = (2.0 * t[2] + 2.0 * gam * t[3] + t[4]) * (1.0 + gam) + 2.0 * gam * t[1] + tail + 1.0e-35;
                float2 c = float2(t[0] - W, t[0] + W);
                float nh = t[5] * (1.0 + gam) + tail + 1.0e-35;
                float2 nat = float2(pr.z != 0.0 ? -pr.w : -nh, nh);
                rho = float2(max(c.x, nat.x), min(c.y, nat.y));
                if (rho.x > rho.y) rho = float2(min(c.x, nat.x), max(c.y, nat.y));   // only by rounding; keep the hull
            }
            for (uint f = 0u; f < a.num_fields; ++f) {
                uint fbits = (node.w >> (2u * f)) & 3u;
                if (fbits == 0u || a.fdmat[f] != u) continue;
                float iso = a.fields[f].iso;
                uint bits = 0u;
                if (rho.x <= iso && rho.y >= iso)   bits |= 1u;
                if (rho.x <= -iso && rho.y >= -iso) bits |= 2u;
                bits &= fbits;
                if (bits != 0u) InterlockedOr(o_mask, bits << (2u * f));
            }
        }
    }
    octree_emit(node, ni, tid);
}

#define EL_GROUP 128
groupshared float4 el_lohi[2 * 32];

// One orbital (one function, one field): the whole Newton projection of 128 consecutive vertices per
// group (a brick's vertices are contiguous) in one launch, instead of an evaluation and a step
// (iso_newton) per iteration, each a launch over every vertex and a barrier however few vertices still
// move. Per iteration the group lists and stages the shells that reach its vertices' box with their
// coefficients, as for the corners, and evaluates every vertex from groupshared memory; then the vertices
// still moving take their step (as iso_newton); after a.level steps, or once all of the group's vertices
// have converged, the normals (iso_newton, mode 1).
groupshared uint ep_any;

[shader("compute")]
[numthreads(EL_GROUP, 1, 1)]
void iso_project_list_f1(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint v0 = node_index(gid) * uint(EL_GROUP);
    if (v0 >= a.count) return;                             // uniform per group
    bool active = v0 + tid < a.count;
    uint vi = min(v0 + tid, a.count - 1u);
    float4 pv = a.vpos[vi];
    float3 p = pv.xyz;
    float sg = pv.w;
    float iso = a.fields[a.vinfo[vi].z].iso;
    IsoTerm tm = a.terms[a.fields[0].term_off];
    uint done = active ? 0u : 1u;
    float4 gv = float4(0.0);
    uint iters = a.level;
    for (uint it = 0u; it <= iters; ++it) {
        // any of the group's vertices still moving? (uniform)
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0u) ep_any = 0u;
        GroupMemoryBarrierWithGroupSync();
        if (done == 0u) ep_any = 1u;
        GroupMemoryBarrierWithGroupSync();
        if (ep_any == 0u) break;
        // the evaluation: the vertices' box, its shells listed and staged by blocks
        {
            uint lanes = WaveGetLaneCount();
            uint w = min(iso_wave_index(tid, lanes), 31u);
            float3 lo = WaveActiveMin(active ? p : float3(3.0e38));
            float3 hi = WaveActiveMax(active ? p : float3(-3.0e38));
            if (WaveIsFirstLane()) { el_lohi[2u * w] = float4(lo, 0.0); el_lohi[2u * w + 1u] = float4(hi, 0.0); }
            GroupMemoryBarrierWithGroupSync();
            if (tid == 0u) {
                uint nw = min(iso_num_waves(EL_GROUP, lanes), 32u);
                float3 glo = el_lohi[0].xyz, ghi = el_lohi[1].xyz;
                for (uint k = 1u; k < nw; ++k) { glo = min(glo, el_lohi[2u * k].xyz); ghi = max(ghi, el_lohi[2u * k + 1u].xyz); }
                el_lohi[0] = float4(glo, 0.0);
                el_lohi[1] = float4(ghi, 0.0);
            }
            GroupMemoryBarrierWithGroupSync();
        }
        float3 glo = el_lohi[0].xyz, ghi = el_lohi[1].xyz;
        float val = 0.0;
        float3 grad = float3(0.0);
        for (uint s0 = 0u; s0 < a.num_shells; s0 += uint(K_LIST)) {
            GroupMemoryBarrierWithGroupSync();             // the previous block's lists (and the box) consumed
            k_list(s0, glo, ghi, tid, EL_GROUP, K_CK);
            uint nchunk = k_nchunk;
            for (uint c = 0u; c < nchunk; ++c) {
                k_load_chunk<1, KBuf1>(c, 0u, tid, EL_GROUP);
                uint i0 = k_chunk[c], i1 = k_chunk[c + 1];
                for (uint i = i0; i < i1; ++i) {           // uniform
                    float3 d;
                    float R, R1;
                    uint l, nc;
                    if (!k_shell_d(i - i0, p, d, R, R1, l, nc)) continue;
                    float m[15], gx[15], gy[15], gz[15];
                    cart_monomials(l, d, m);
                    cart_monomials_d(l, d, 0u, gx);
                    cart_monomials_d(l, d, 1u, gy);
                    cart_monomials_d(l, d, 2u, gz);
                    uint row0 = k_row[i];
                    float P = 0.0;
                    float3 Pg = float3(0.0);
                    [ForceUnroll] for (int n = 0; n < 15; ++n) {
                        if (uint(n) < nc) {                 // uniform
                            float cv = KBuf1::ld(row0 + uint(n));
                            P += cv * m[n];
                            Pg += cv * float3(gx[n], gy[n], gz[n]);
                        }
                    }
                    val += P * R;
                    grad += Pg * R + (2.0 * P * R1) * d;
                }
                GroupMemoryBarrierWithGroupSync();         // the chunk's rows consumed
            }
        }
        if (active) gv = tm.w * float4(grad, val);
        if (it == iters) break;                            // uniform: the last evaluation, for the normals
        // the step of a vertex still moving (iso_newton)
        if (done == 0u) {
            float f = sg * gv.w - iso;
            float3 g = sg * gv.xyz;
            float g2 = dot(g, g);
            if (g2 > 1.0e-30) {
                float gl = sqrt(g2);
                float sd = abs(f) / gl;
                if (sd < a.tol) {
                    done = 1u;                             // converged: stays put, its values stay valid
                } else {
                    float step = f / g2;
                    float len = abs(step) * gl;
                    if (len > a.step_max) step *= a.step_max / len;
                    p = p - step * g;
                    if (len <= a.step_max) {      // the error after the step predicted (iso_newton)
                        float K = gl / max(abs(sg * gv.w), 1.0e-30);
                        if (0.5 * K * sd * sd < 0.1 * a.tol) done = 2u;
                    }
                }
            }
        }
    }
    if (!active) return;
    // the normals from the last evaluation (iso_newton, mode 1)
    float f = sg * gv.w - iso;
    float3 g = sg * gv.xyz;
    float g2 = dot(g, g);
    float4 nrm = float4(0.0, 0.0, 0.0, 1.0e30);
    if (g2 > 1.0e-30) {
        float gl = sqrt(g2);
        nrm = float4(-g / gl, abs(f) / gl);                // f > 0 inside: outward is -grad f
    }
    a.vpos[vi] = float4(p, sg);
    a.vval[vi] = gv;
    a.vnrm[vi] = nrm;
    a.vdone[vi] = done;
}

// ------------------------------------------------------------------------------------------
// Points: the vertex evaluators (iso_geval_*, iso_deval, iso_eval_orbital) at the caller's points
// (md_gto_iso_gpu_evaluate_launch). A point is a vertex of field pts_field that has not converged.
// ------------------------------------------------------------------------------------------

// The count points from pts as vertices: mode 1 their positions (and vdone cleared), always the field
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_points_load(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint i = tid.x;
    if (i >= a.count) return;
    if ((a.mode & 1u) != 0u) {
        uint o = i * a.pts_stride;
        a.vpos[i] = float4(a.pts[o], a.pts[o + 1u], a.pts[o + 2u], 1.0);
        a.vdone[i] = 0u;
    }
    a.vinfo[i] = uint4(0u, 0u, a.pts_field, 0u);
}

// The values of the count points evaluated (vval: gradient xyz, value w) to pval
[shader("compute")]
[numthreads(64, 1, 1)]
void iso_points_store(uint3 tid : SV_DispatchThreadID) {
    IsoArgs a = MD_ARGS;
    uint i = tid.x;
    if (i >= a.count) return;
    a.pval[i] = a.vval[i].w;
}

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

// One group per brick, one thread per cell. mode 0 counts the vertices (per brick into scan, the total
// into counters[0]), mode 1 writes them at the brick's offset in scan (iso_scan of the counts): a brick's
// vertices contiguous, in the order of the bricks.
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
        if (a.mode == 0u) {
            a.scan[b] = total;
            if (total > 0u) InterlockedAdd(a.counters[0], total);
        } else {
            c_base = a.scan[b];
        }
    }
    GroupMemoryBarrierWithGroupSync();
    if ((a.mode & 1u) == 0u) return;

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
// Vertices: Newton steps from the evaluated values (vval), relaxation
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
    if (len <= a.step_max) {
        // the error after this step predicted instead of checked by another evaluation: about
        // |f''| / (2 |f'|) s^2 along the gradient, with |f''| / |f'| ~ |f'| / |f| for the exponentially
        // decaying fields here. Predicted well below the tolerance: done (its normal and values then
        // those of the evaluation before the step)
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
    uint nq = min(a.counters[1], a.cap);     // the quads made (cap: their room)
    if (i < a.count) InterlockedAdd(a.bucket[vertex_bucket(i)], 1u);
    if (i < nq)      InterlockedAdd(a.bucket[32u + vertex_bucket(a.quads[i].x)], 1u);
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
    if (i >= min(a.counters[1], a.cap)) return;
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

// The D form's matrices in tiles: per matrix (gid.y) and AO mu, sum over nu of |D_mu nu| (dm_rows), and the octree's
// tails drs = |nrm_mu| times that (as the dense form's: sum over nu of |D_mu nu nrm_mu nrm_nu| / |nrm_nu|), raised
// by 1e-3 for the float sum's rounding. A tile pair (I, J), I <= J, holds rows of I, columns of J; a diagonal one
// its upper triangle.
[shader("compute")] [numthreads(64, 1, 1)]
void iso_dmat_rows(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    IsoArgs a = MD_ARGS;
    uint i = gid.x * 64u + tid;
    uint u = gid.y;
    if (i >= a.num_ao) return;
    uint nt = a.dm_nt;
    uint ti = a.dm_tao[i];
    uint I = ti >> 4, r = ti & 15u;
    float* T = a.dmat + u * a.dm_pool;
    float sum = 0.0;
    for (uint J = 0u; J < nt; ++J) {
        uint e = a.dm_map[min(I, J) * nt + max(I, J)];
        if (e == INVALID) continue;
        float* t = T + e * 256u;
        if (J > I) {
            for (uint c = 0u; c < 16u; ++c) sum += abs(t[r * 16u + c]);
        } else if (J < I) {
            for (uint c = 0u; c < 16u; ++c) sum += abs(t[c * 16u + r]);
        } else {
            for (uint c = 0u; c < 16u; ++c) sum += abs(c >= r ? t[r * 16u + c] : t[c * 16u + r]);
        }
    }
    a.dm_rows[u * a.num_ao + i] = sum;
    a.drs[u * a.num_ao + i] = abs(a.ao_nrm[i]) * sum * 1.001;
}

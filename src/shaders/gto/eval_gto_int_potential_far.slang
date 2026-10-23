// eval_gto_int_potential_far.slang
//
// The electrostatic potential of an md_gto_int charge distribution that has a far field
// (md_gto_int_charges_desc_t::far_field_theta), on a grid (md_gto_int_gpu_potential_launch) or at
// points (md_gto_int_gpu_potential_xyz_launch). The tree and the criterion are those of the CPU path
// in md_gto_int.c; the host lays the distribution out for it (md_gto_int_gpu_charges_create).
//
// GROUPS. The points are taken 64 at a time: a brick of 4x4x4 voxels of the grid, or 64 points
// consecutive along a Morton curve through their box. A group shares one walk of the tree, against
// its bounding spheres: a node is taken whole by its multipoles when, for every sphere,
//     |centre - sphere centre| - sphere radius >= the node's far distance
// (max(radius / theta, reach), see md_gto_int.c), a leaf that is not is evaluated exactly, and
// anything else is opened. A brick is one sphere. 64 points are four: split where the three largest
// steps between consecutive points are. Points that fill space only in part (a surface, a band) leave
// a Morton curve for empty space and come back elsewhere, and one sphere around both ends of such a
// jump holds everything between them: on a protein's ESP band a few groups so made had over 20 times
// the median work, and a dispatch waits for its slowest group.
//
// ORDER (points only). The points come in any order: a mesh's vertices are close to one another in
// the order a mesher makes them, but one that hands out its output slots by atomics interleaves the
// workgroups in flight, and 64 points from all over the molecule make a group whose sphere holds all
// of it - exact, and over ten times slower on a protein. So the points are sorted first, by Morton key in their box
// with the index as the tie break, by a bitonic network over global memory:
//   xyz_bounds    the box, by atomic min and max of the coordinates as ordered integers
//   xyz_keys      (key, index) per point, padded to a power of two with keys that sort last
//   xyz_bitonic   one compare and swap stage (j, k) per dispatch, log2(n) (log2(n) + 1) / 2 of them
// The result does not depend on the order the points come in, nor on the timing of the atomics. A
// point's accumulator is at its place in the sorted order; resolve writes it back to its index.
//
// PASSES, each for the grid (grid_*) and for points (xyz_*):
//   lists     one thread per group walks the tree and writes the group's list: the far nodes, and
//             the near leaves. It never runs out of room: when opening a node would leave less room
//             than the walk could still need, the node goes on the list whole, to be evaluated
//             exactly (its members are one run per order), which is slower and still right.
//   far<K>    a workgroup per group, a thread per point: the far nodes' multipoles (order K), and
//             the near points (the nuclei)
//   near<L>   likewise: the near gaussians of order L, so that each kernel is sized for one recursion
//   resolve   converts and writes
// adding to one 64 bit fixed point accumulator per point (gto_int_potential_common.slang). Every
// thread of a workgroup reads the same list in the same order: nothing diverges, and there is no
// groupshared memory and no barrier anywhere. (A first version walked the tree cooperatively,
// 64 nodes a round through groupshared memory between barriers; lavapipe read garbage out of it on
// a large tree, in a way that moved with unrelated changes to the code.)
//
// A far node is evaluated as a gaussian of order K whose R^n_000 are those of 1/r,
//     s_n = (-1)^n (2n-1)!! / r^(2n+1),
// with its multipoles as the coefficients: sum_k M_k d^k (1/|Y|), Y = centre - C.

#include "../md_gpu.slang"
#include "gto_common.slang"
#include "gto_int_hermite.slang"

#define GROUP      64
#define WALK_STACK 256     // a walk's stack, per thread: depth 32 of an octree needs 7 * 32 + 1

#define GTO_INT_BOYS_DT     0.1f
#define GTO_INT_BOYS_INV_DT 10.0f
#define GTO_INT_BOYS_TMAX   30.0f
#define GTO_INT_BOYS_NN     16

struct RootArgs {
    float4x4 index_to_world;   // grid: voxel index (sample offset included) -> world - ref, bohr
    uint4    grid_dim;         // grid: voxels, xyz
    uint4    brick_dim;        // grid: bricks of 4x4x4 voxels, xyz
    float4   ref;              // points: the frame the distribution is stored in, bohr
    uint     list_cap;         // entries of a group's list
    uint     group_first;      // the first group of this dispatch
    uint     group_end;        // lists: one past the last group of this batch
    uint     accum_first;      // the group whose accumulators and list start the buffers
    uint     num_xyz;          // points
    uint     xyz_stride;       // points: floats between points
    uint     num_ranges;       // per node: one range per order 0..max_order, then the points'
    uint     operation;        // resolve, grid: md_gto_op_t
    float    brick_radius;     // grid: half the diagonal between a brick's voxel centres, bohr
    uint     far_nh;           // coefficients per node: num_hermite(K)
    uint     coeff_first;      // near<L>: the first coefficient of order L ...
    uint     order_first;      // ... which is that of gaussian order_first, the first of order L
    uint     sort_j;           // xyz_bitonic: the partner distance
    uint     sort_k;           // xyz_bitonic: the size of the sequences being merged
    uint     num_sorted;       // points: the sort's length, a power of two >= num_xyz
    uint     sort_stride;      // xyz_bounds, xyz_keys, xyz_bitonic: threads in the dispatch
    float4*  gaussians;        // (P.xyz - ref, p), sorted by order, within an order by node
    float*   coeffs;           // h', contiguous per order in the same order
    float*   boys;             // [rows][GTO_INT_BOYS_NN]
    float4*  points;           // (xyz - ref, q), by node
    float4*  node_geom;        // (centre - ref, far distance)
    uint4*   node_link;        // (first child, number of children, 0, 0)
    uint2*   node_range;       // [node * num_ranges + r]: gaussians of order r, then the points
    float*   multipoles;       // [node * far_nh]
    uint*    bricks;           // grid: the bricks to evaluate (x + bx (y + by z))
    float*   xyz;              // points, bohr: xyz[i * xyz_stride + 0..2]
    uint2*   accum;            // [64 per group], from accum_first on
    float*   out;              // points: the potential
    uint*    lists;            // [list_cap per group], from accum_first on; FAR_BIT: by its multipoles
    uint*    counts;           // [per group], from accum_first on
    uint2*   sorted;           // points: [num_sorted] (Morton key, index), sorted
    uint*    bounds;           // points: [6] their box, min xyz then max xyz, as float_key
    DescriptorHandle<RWTexture3D<float>> out_tex;
    uint2    _pad_end;         // the struct is 16 byte aligned: its size, said out loud for the C mirror
};
struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;
#define args (*md_root.args_ptr)

#include "gto_int_potential_common.slang"

#define FAR_BIT 0x80000000u

// The grid and the points are separate entry points (M: GRID or POINTS) rather than one with a
// branch on a mode: a compiler may evaluate both sides of such a branch - llvmpipe does - and the
// side not taken loads through a pointer that is not there. For the same reason the grid always
// has a brick list, all of its bricks when every voxel is wanted.
#define GRID   0
#define POINTS 1

uint3 group_brick(uint group) {
    const uint b = args.bricks[group];
    const uint bx = args.brick_dim.x, by = args.brick_dim.y;
    return uint3(b % bx, (b / bx) % by, b / (bx * by));
}

float3 point_xyz(uint i) {
    const uint o = i * args.xyz_stride;
    return float3(args.xyz[o + 0u], args.xyz[o + 1u], args.xyz[o + 2u]) - args.ref.xyz;
}

// The point at place i of the sorted order
float3 sorted_xyz(uint i) { return point_xyz(args.sorted[i].y); }

struct GroupPoint {
    float3 C;       // the point, in the distribution's frame
    bool   valid;   // a voxel of the grid, or a point of the set
    uint3  vox;     // grid: its voxel
};

GroupPoint group_point<let M : int>(uint group, uint lid) {
    GroupPoint p;
    if (M == GRID) {
        const uint3 local = uint3(lid & 3u, (lid >> 2u) & 3u, lid >> 4u);
        p.vox   = group_brick(group) * 4u + local;
        p.valid = all(p.vox < args.grid_dim.xyz);
        p.C     = mul(args.index_to_world, float4(float3(p.vox), 1.0f)).xyz;
    } else {
        const uint i = group * GROUP + lid;
        p.vox   = uint3(0u, 0u, 0u);
        p.valid = i < args.num_xyz;
        p.C     = sorted_xyz(min(i, args.num_xyz - 1u));
    }
    return p;
}

// ---------------------------------------------------------------------------
// The points' order
// ---------------------------------------------------------------------------
// A float as an unsigned integer with the same order (for the atomics), and back
uint float_key(float f) {
    const uint u = asuint(f);
    return (u & 0x80000000u) != 0u ? ~u : (u | 0x80000000u);
}
float key_float(uint k) {
    return asfloat((k & 0x80000000u) != 0u ? (k & 0x7FFFFFFFu) : ~k);
}

// The ten low bits of v, two zeros after each
uint morton_spread(uint v) {
    v &= 0x3FFu;
    v = (v | (v << 16u)) & 0x030000FFu;
    v = (v | (v << 8u))  & 0x0300F00Fu;
    v = (v | (v << 4u))  & 0x030C30C3u;
    v = (v | (v << 2u))  & 0x09249249u;
    return v;
}

// A thread per 64 points (grid-stride); points that are not finite are left out of the box
[shader("compute")] [numthreads(GROUP, 1, 1)]
void xyz_bounds(uint3 tid : SV_DispatchThreadID) {
    float3 lo = float3(3.0e38f, 3.0e38f, 3.0e38f), hi = -lo;
    for (uint c = tid.x; c * GROUP < args.num_xyz; c += args.sort_stride) {
        const uint last = min(c * GROUP + GROUP, args.num_xyz);
        for (uint i = c * GROUP; i < last; ++i) {
            const float3 x = point_xyz(i);
            if (!all(abs(x) < 1.0e30f)) continue;   // false for NaN as well
            lo = min(lo, x);
            hi = max(hi, x);
        }
    }
    if (lo.x > hi.x) return;   // none here
    InterlockedMin(args.bounds[0], float_key(lo.x));
    InterlockedMin(args.bounds[1], float_key(lo.y));
    InterlockedMin(args.bounds[2], float_key(lo.z));
    InterlockedMax(args.bounds[3], float_key(hi.x));
    InterlockedMax(args.bounds[4], float_key(hi.y));
    InterlockedMax(args.bounds[5], float_key(hi.z));
}

// A thread per place of the sort (grid-stride): a point's key, 10 bits an axis of its box, or one
// past every key for the padding
[shader("compute")] [numthreads(GROUP, 1, 1)]
void xyz_keys(uint3 tid : SV_DispatchThreadID) {
    const float3 lo  = float3(key_float(args.bounds[0]), key_float(args.bounds[1]), key_float(args.bounds[2]));
    const float3 hi  = float3(key_float(args.bounds[3]), key_float(args.bounds[4]), key_float(args.bounds[5]));
    const float3 inv = 1024.0f / max(hi - lo, float3(1.0e-6f, 1.0e-6f, 1.0e-6f));
    for (uint i = tid.x; i < args.num_sorted; i += args.sort_stride) {
        uint key = 0xFFFFFFFFu;
        if (i < args.num_xyz) {
            // whatever a point that is not finite makes of this, morton_spread keeps it to 30 bits
            const uint3 q = uint3(clamp((point_xyz(i) - lo) * inv, float3(0.0f, 0.0f, 0.0f), float3(1023.0f, 1023.0f, 1023.0f)));
            key = morton_spread(q.x) | (morton_spread(q.y) << 1u) | (morton_spread(q.z) << 2u);
        }
        args.sorted[i] = uint2(key, i);
    }
}

// One stage of the bitonic network: a thread per pair (grid-stride), i with bit j clear and its
// partner i + j, in the direction of the k-sequence i is in. Points in one cell of the key and the
// padding tie on the key, the index breaks the ties: no two entries are equal, and the result is
// the same whatever order the threads ran in.
[shader("compute")] [numthreads(GROUP, 1, 1)]
void xyz_bitonic(uint3 tid : SV_DispatchThreadID) {
    const uint j = args.sort_j, k = args.sort_k;
    for (uint t = tid.x; t < args.num_sorted / 2u; t += args.sort_stride) {
        const uint i = ((t & ~(j - 1u)) << 1u) | (t & (j - 1u));
        const uint l = i | j;
        const uint2 a = args.sorted[i];
        const uint2 b = args.sorted[l];
        const bool greater = a.x > b.x || (a.x == b.x && a.y > b.y);
        if (greater == ((i & k) == 0u)) {
            args.sorted[i] = b;
            args.sorted[l] = a;
        }
    }
}

// ---------------------------------------------------------------------------
// lists: a thread per group
// ---------------------------------------------------------------------------
#define GROUP_SPHERES 4

// The sorted points [b, e): the centre of their box and the largest distance from it, with a little
// slack for the rounding of both (the points are float, and so is the centre)
float4 segment_sphere(uint b, uint e) {
    float3 lo = sorted_xyz(b), hi = lo;
    for (uint i = b + 1u; i < e; ++i) {
        const float3 x = sorted_xyz(i);
        lo = min(lo, x);
        hi = max(hi, x);
    }
    const float3 c = 0.5f * (lo + hi);
    float r = 0.0f;
    for (uint i = b; i < e; ++i) r = max(r, length(sorted_xyz(i) - c));
    return float4(c, r * 1.0001f + 1.0e-5f);
}

// A group's spheres: the brick's, or the points' split at the three largest steps between them.
// Unused ones repeat the first.
void group_spheres<let M : int>(uint group, out float4 s[GROUP_SPHERES]) {
    if (M == GRID) {
        s[0] = float4(mul(args.index_to_world, float4(float3(group_brick(group) * 4u) + 1.5f, 1.0f)).xyz, args.brick_radius);
        [ForceUnroll] for (int k = 1; k < GROUP_SPHERES; ++k) s[k] = s[0];
        return;
    }
    const uint first = group * GROUP;
    const uint last  = min(first + GROUP, args.num_xyz);
    // the three largest steps, as the point after each: c0 >= c1 >= c2 by step
    const uint NONE = 0xFFFFFFFFu;
    uint  c0 = NONE, c1 = NONE, c2 = NONE;
    float v0 = -1.0f, v1 = -1.0f, v2 = -1.0f;
    float3 prev = sorted_xyz(first);
    for (uint i = first + 1u; i < last; ++i) {
        const float3 x = sorted_xyz(i);
        const float3 dx = x - prev;
        const float  v = dot(dx, dx);
        prev = x;
        if (v > v0)      { v2 = v1; c2 = c1; v1 = v0; c1 = c0; v0 = v; c0 = i; }
        else if (v > v1) { v2 = v1; c2 = c1; v1 = v; c1 = i; }
        else if (v > v2) { v2 = v; c2 = i; }
    }
    // in order along the points (NONE last)
    uint t;
    if (c0 > c1) { t = c0; c0 = c1; c1 = t; }
    if (c1 > c2) { t = c1; c1 = c2; c2 = t; }
    if (c0 > c1) { t = c0; c0 = c1; c1 = t; }
    const uint e0 = min(c0, last), e1 = min(c1, last), e2 = min(c2, last);
    s[0] = segment_sphere(first, e0);
    s[1] = e1 > e0 ? segment_sphere(e0, e1) : s[0];
    s[2] = e2 > e1 ? segment_sphere(e1, e2) : s[0];
    s[3] = last > e2 ? segment_sphere(e2, last) : s[0];
}

// How far a node's centre is outside the group's spheres: the least over them
float group_distance<let M : int>(float3 c, float4 s[GROUP_SPHERES]) {
    float d = length(c - s[0].xyz) - s[0].w;
    if (M != GRID) {
        [ForceUnroll] for (int k = 1; k < GROUP_SPHERES; ++k) d = min(d, length(c - s[k].xyz) - s[k].w);
    }
    return d;
}

void lists_pass<let M : int>(uint tid) {
    const uint group = args.group_first + tid;
    if (group >= args.group_end) return;
    float4 spheres[GROUP_SPHERES];
    group_spheres<M>(group, spheres);
    const uint   cap    = args.list_cap;
    const uint   base   = (group - args.accum_first) * cap;

    // count + sp <= cap throughout: popping a node and listing it leaves the sum as it was, and a
    // node is only opened when its children fit on the stack and in what is left of the list
    uint stack[WALK_STACK];
    uint sp = 1u, count = 0u;
    stack[0] = 0u;   // the root
    while (sp > 0u) {
        const uint node = stack[--sp];
        const float4 g = args.node_geom[node];
        const float  d = group_distance<M>(g.xyz, spheres);
        if (d > 0.0f && d >= g.w) {
            args.lists[base + count++] = node | FAR_BIT;
            continue;
        }
        const uint4 l = args.node_link[node];
        if (l.y == 0u || sp + l.y > uint(WALK_STACK) || count + sp + l.y > cap) {
            args.lists[base + count++] = node;   // a near leaf, or a node evaluated whole
            continue;
        }
        for (uint c = 0u; c < l.y; ++c) stack[sp++] = l.x + c;
    }
    args.counts[group - args.accum_first] = count;
}

[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_lists(uint3 tid : SV_DispatchThreadID) { lists_pass<GRID>(tid.x); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_lists(uint3 tid : SV_DispatchThreadID) { lists_pass<POINTS>(tid.x); }

// ---------------------------------------------------------------------------
// far<K>: the far nodes and the near points
// ---------------------------------------------------------------------------
float far_node<let K : int>(uint node, float3 C) {
    const float3 Y  = args.node_geom[node].xyz - C;
    const float  r2 = dot(Y, Y);
    const float  ir2 = 1.0f / r2;
    float s[GTO_INT_NS];
    s[0] = sqrt(ir2);
    [ForceUnroll] for (int n = 0; n < K; ++n) s[n + 1] = -float(2 * n + 1) * s[n] * ir2;
    [ForceUnroll] for (int n = K + 1; n < GTO_INT_NS; ++n) s[n] = 0.0f;
    const uint o = node * args.far_nh;
    if (K == 6) return gto_int_hermite_L6(Y.x, Y.y, Y.z, s, args.multipoles, o);
    if (K == 7) return gto_int_hermite_L7(Y.x, Y.y, Y.z, s, args.multipoles, o);
    return gto_int_hermite_L8(Y.x, Y.y, Y.z, s, args.multipoles, o);
}

void far_pass<let K : int, let M : int>(uint group_in_dispatch, uint lid) {
    const uint group = args.group_first + group_in_dispatch;
    const GroupPoint p = group_point<M>(group, lid);
    if (!p.valid) return;
    const uint slot = group - args.accum_first;
    const uint base = slot * args.list_cap;
    const uint n    = args.counts[slot];
    const uint pr   = args.num_ranges - 1u;   // the points' range
    uint2 acc = args.accum[slot * GROUP + lid];
    for (uint e = 0u; e < n; ++e) {
        const uint it = args.lists[base + e];
        if ((it & FAR_BIT) != 0u) {
            acc = gto_int_fixed_add(acc, far_node<K>(it & ~FAR_BIT, p.C));
        } else {
            const uint2 r = args.node_range[it * args.num_ranges + pr];
            for (uint k = r.x; k < r.y; ++k) {
                const float4 q = args.points[k];
                const float3 d = p.C - q.xyz;
                const float r2 = dot(d, d);
                if (r2 > 1.0e-20f) acc = gto_int_fixed_add(acc, q.w * rsqrt(r2));
            }
        }
    }
    args.accum[slot * GROUP + lid] = acc;
}

[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_far6(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<6, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_far7(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<7, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_far8(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<8, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_far6(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<6, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_far7(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<7, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_far8(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { far_pass<8, POINTS>(gid.x, lid); }

// ---------------------------------------------------------------------------
// near<L>: the near gaussians of order L
// ---------------------------------------------------------------------------
void near_pass<let L : int, let M : int>(uint group_in_dispatch, uint lid) {
    const uint group = args.group_first + group_in_dispatch;
    const GroupPoint p = group_point<M>(group, lid);
    if (!p.valid) return;
    const uint slot = group - args.accum_first;
    const uint base = slot * args.list_cap;
    const uint n    = args.counts[slot];
    const uint nh   = uint(((L + 1) * (L + 2) * (L + 3)) / 6);
    uint2 acc = args.accum[slot * GROUP + lid];
    for (uint e = 0u; e < n; ++e) {
        const uint it = args.lists[base + e];
        if ((it & FAR_BIT) != 0u) continue;
        const uint2 r = args.node_range[it * args.num_ranges + uint(L)];
        for (uint g = r.x; g < r.y; ++g) {
            const float4 gs = args.gaussians[g];
            const float3 X = gs.xyz - p.C;
            float s[GTO_INT_NS];
            gto_int_boys_scaled<L>(s, gs.w * dot(X, X), gs.w);
            acc = gto_int_fixed_add(acc, gto_int_hermite<L>(X, s, args.coeff_first + (g - args.order_first) * nh));
        }
    }
    args.accum[slot * GROUP + lid] = acc;
}

[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near0(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<0, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near1(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<1, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near2(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<2, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near3(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<3, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near4(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<4, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near5(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<5, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near6(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<6, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near7(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<7, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_near8(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<8, GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near0(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<0, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near1(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<1, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near2(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<2, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near3(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<3, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near4(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<4, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near5(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<5, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near6(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<6, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near7(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<7, POINTS>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_near8(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { near_pass<8, POINTS>(gid.x, lid); }

// ---------------------------------------------------------------------------
// resolve: convert, and write the voxel (with the md_gto_op_t) or the point
// ---------------------------------------------------------------------------
void resolve_pass<let M : int>(uint group_in_dispatch, uint lid) {
    const uint group = args.group_first + group_in_dispatch;
    const GroupPoint p = group_point<M>(group, lid);
    if (!p.valid) return;
    const float V = gto_int_fixed_to_float(args.accum[(group - args.accum_first) * GROUP + lid]);
    if (M == GRID) {
        const float prev = gto_op_reads_prev(args.operation) ? args.out_tex[p.vox] : 0.0f;
        args.out_tex[p.vox] = gto_apply_op(V, prev, args.operation);
    } else {
        args.out[args.sorted[group * GROUP + lid].y] = V;
    }
}

[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_resolve(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { resolve_pass<GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_resolve(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { resolve_pass<POINTS>(gid.x, lid); }

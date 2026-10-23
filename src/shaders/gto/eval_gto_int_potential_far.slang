// eval_gto_int_potential_far.slang
//
// The electrostatic potential of an md_gto_int charge distribution that has a far field
// (md_gto_int_charges_desc_t::far_field_theta), on a grid (md_gto_int_gpu_potential_launch) or at
// points (md_gto_int_gpu_potential_xyz_launch). The tree and the criterion are those of the CPU path
// in md_gto_int.c; the host lays the distribution out for it (md_gto_int_gpu_charges_create).
//
// GROUPS. A workgroup evaluates a group of 64 points: a brick of 4x4x4 voxels of the grid, or 64
// consecutive points (consecutive points of a mesh are close to one another). The group walks the
// tree from the root together, against its bounding sphere, so every thread of it sees the same
// nodes in the same order and nothing diverges: a node is taken whole by its multipoles when
//     |centre - group centre| - group radius >= the node's far distance
// (max(radius / theta, reach), see md_gto_int.c), a leaf that is not is evaluated exactly, and
// anything else is opened. The walk is breadth first in rounds: the threads take up to 64 nodes off
// a stack in groupshared memory, one each, decide them, and push the children of the ones to open;
// then all of them evaluate the round's list at their own point. Should the stack ever be full, a
// node is evaluated exactly as a whole instead of opened (its members are one run per order), which
// is slower and still right.
//
// PASSES, as in eval_gto_int_potential.slang, so that each kernel is sized for one recursion, each
// for the grid (grid_*) and for points (xyz_*):
//   far<K>    the far nodes' multipoles (order K) and the near points (the nuclei)
//   near<L>   the near gaussians of order L
//   resolve   converts and writes
// all adding to one 64 bit fixed point accumulator per point (see gto_int_potential_common.slang).
// Every pass walks the tree again: a walk is a few thousand node tests per group, 64 at a time,
// against hundreds of thousands of operations per point for the evaluation.
//
// A far node is evaluated as a gaussian of order K whose R^n_000 are those of 1/r,
//     s_n = (-1)^n (2n-1)!! / r^(2n+1),
// with its multipoles as the coefficients: sum_k M_k d^k (1/|Y|), Y = centre - C.

#include "../md_gpu.slang"
#include "gto_common.slang"
#include "gto_int_hermite.slang"

#define GROUP     64
#define STACK_CAP 1024

#define GTO_INT_BOYS_DT     0.1f
#define GTO_INT_BOYS_INV_DT 10.0f
#define GTO_INT_BOYS_TMAX   30.0f
#define GTO_INT_BOYS_NN     16

struct RootArgs {
    float4x4 index_to_world;   // grid: voxel index (sample offset included) -> world - ref, bohr
    uint4    grid_dim;         // grid: voxels, xyz
    uint4    brick_dim;        // grid: bricks of 4x4x4 voxels, xyz
    float4   ref;              // points: the frame the distribution is stored in, bohr
    uint     _pad_mode;
    uint     group_first;      // the first group of this dispatch
    uint     accum_first;      // the group whose 64 accumulators start accum
    uint     num_groups;       // grid: bricks (listed, or all); points: ceil(num_xyz / 64)
    uint     num_xyz;          // points
    uint     xyz_stride;       // points: floats between points
    uint     num_ranges;       // per node: one range per order 0..max_order, then the points'
    uint     operation;        // resolve, grid: md_gto_op_t
    float    brick_radius;     // grid: half the diagonal between a brick's voxel centres, bohr
    uint     far_nh;           // coefficients per node: num_hermite(K)
    uint     coeff_base;       // near<L>: gaussian g's coefficients start at coeff_base + g * num_hermite(L) (mod 2^32)
    uint     _pad1;
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
    DescriptorHandle<RWTexture3D<float>> out_tex;
    uint2    _pad_end;         // the struct is 16 byte aligned: its size, said out loud for the C mirror
};
struct MdGpuRoot { RootArgs* args_ptr; };
[[vk::push_constant]]
ConstantBuffer<MdGpuRoot> md_root;
#define args (*md_root.args_ptr)

#include "gto_int_potential_common.slang"

groupshared uint   s_stack[STACK_CAP];
groupshared uint   s_sp;
groupshared uint   s_list[GROUP];      // this round's nodes to evaluate; FAR_BIT: by its multipoles
groupshared uint   s_nlist;
groupshared float4 s_sphere;           // the group's bounding sphere: centre, radius (no float3 here: see below)

#define FAR_BIT 0x80000000u

// ---------------------------------------------------------------------------
// The group: which points, and its bounding sphere
// ---------------------------------------------------------------------------
struct GroupPoint {
    float3 C;       // the point, in the distribution's frame
    bool   valid;   // a voxel of the grid, or a point of the set
    uint3  vox;     // grid: its voxel
};

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
        const uint o = min(i, args.num_xyz - 1u) * args.xyz_stride;
        p.C = float3(args.xyz[o + 0u], args.xyz[o + 1u], args.xyz[o + 2u]) - args.ref.xyz;
    }
    return p;
}

// The bounding sphere into s_sphere: a brick's is known, a run of points' is the centre
// of their box and the largest distance from it (64 points: thread 0 does it, which is nothing
// beside what the group goes on to do).
void group_bounds<let M : int>(uint group, uint lid, GroupPoint p) {
    if (M == GRID) {
        if (lid == 0u) {
            s_sphere = float4(mul(args.index_to_world, float4(float3(group_brick(group) * 4u) + 1.5f, 1.0f)).xyz, args.brick_radius);
        }
        GroupMemoryBarrierWithGroupSync();
        return;
    }
    // The points go through the stack, which is not in use yet, rather than through an array of
    // their own: every kernel of this file then has the same groupshared footprint. With an array
    // only the point kernels had, lavapipe (Mesa, LLVM 20) read garbage out of the walk's lists in a
    // point kernel first created after the grid kernels had run, and crashed on it; created before
    // them, it did not. Whatever the cause there, one footprint for all of them avoids it.
    s_stack[lid * 4u + 0u] = asuint(p.C.x);
    s_stack[lid * 4u + 1u] = asuint(p.C.y);
    s_stack[lid * 4u + 2u] = asuint(p.C.z);
    s_stack[lid * 4u + 3u] = p.valid ? 1u : 0u;
    GroupMemoryBarrierWithGroupSync();
    if (lid == 0u) {
        float3 lo = p.C, hi = p.C;   // point 0 of a group always exists
        for (uint i = 1u; i < GROUP; ++i) {
            if (s_stack[i * 4u + 3u] != 0u) {
                const float3 x = float3(asfloat(s_stack[i * 4u + 0u]), asfloat(s_stack[i * 4u + 1u]), asfloat(s_stack[i * 4u + 2u]));
                lo = min(lo, x);
                hi = max(hi, x);
            }
        }
        s_sphere = float4(0.5f * (lo + hi), 0.0f);
    }
    GroupMemoryBarrierWithGroupSync();
    s_stack[lid * 4u + 3u] = asuint(p.valid ? length(p.C - s_sphere.xyz) : 0.0f);
    GroupMemoryBarrierWithGroupSync();
    if (lid == 0u) {
        float r = 0.0f;
        for (uint i = 0u; i < GROUP; ++i) r = max(r, asfloat(s_stack[i * 4u + 3u]));
        // the points are float, and so is the centre: a little slack for the rounding of both
        s_sphere.w = r * 1.0001f + 1.0e-5f;
    }
    GroupMemoryBarrierWithGroupSync();
}

// ---------------------------------------------------------------------------
// The walk
// ---------------------------------------------------------------------------
// One round: decide up to 64 nodes off the stack, push the children of those to open, and leave
// the ones to evaluate in s_list. False when the stack is empty (uniform: every thread reads s_sp
// after the same barrier). with_far: whether the caller evaluates far nodes (far<K>) or only near
// ones (near<L>), which do not need the far ones listed.
bool walk_round(uint lid, bool with_far) {
    GroupMemoryBarrierWithGroupSync();
    const uint sp = s_sp;
    if (sp == 0u) return false;
    const uint n = min(sp, uint(GROUP));
    const uint node = (lid < n) ? s_stack[sp - n + lid] : 0u;
    GroupMemoryBarrierWithGroupSync();
    if (lid == 0u) {
        s_sp = sp - n;
        s_nlist = 0u;
    }
    GroupMemoryBarrierWithGroupSync();
    if (lid < n) {
        const float4 g = args.node_geom[node];
        const uint4  l = args.node_link[node];
        const float  d = length(g.xyz - s_sphere.xyz) - s_sphere.w;
        uint o;
        if (d > 0.0f && d >= g.w) {
            if (with_far) {
                InterlockedAdd(s_nlist, 1u, o);
                s_list[o] = node | FAR_BIT;
            }
        } else if (l.y == 0u) {
            InterlockedAdd(s_nlist, 1u, o);
            s_list[o] = node;
        } else {
            uint at;
            InterlockedAdd(s_sp, l.y, at);
            if (at + l.y <= uint(STACK_CAP)) {
                for (uint c = 0u; c < l.y; ++c) s_stack[at + c] = l.x + c;
            } else {
                // full: take this node back off and evaluate the whole of it exactly. A failed
                // reservation only ever overshoots the true top, so no other thread is handed a
                // slot past it while it stands (see the C reference of this argument in md_gto_int.c)
                InterlockedAdd(s_sp, uint(-int(l.y)), at);
                InterlockedAdd(s_nlist, 1u, o);
                s_list[o] = node;
            }
        }
    }
    GroupMemoryBarrierWithGroupSync();
    return true;
}

void walk_begin(uint lid) {
    if (lid == 0u) {
        s_stack[0] = 0u;   // the root
        s_sp = 1u;
    }
}

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
    group_bounds<M>(group, lid, p);
    walk_begin(lid);
    const uint pr = args.num_ranges - 1u;   // the points' range
    uint2 acc = args.accum[(group - args.accum_first) * GROUP + lid];
    while (walk_round(lid, true)) {
        const uint nl = s_nlist;
        for (uint e = 0u; e < nl; ++e) {
            const uint it = s_list[e];
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
    }
    if (p.valid) args.accum[(group - args.accum_first) * GROUP + lid] = acc;
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
    group_bounds<M>(group, lid, p);
    walk_begin(lid);
    const uint nh = uint(((L + 1) * (L + 2) * (L + 3)) / 6);
    uint2 acc = args.accum[(group - args.accum_first) * GROUP + lid];
    while (walk_round(lid, false)) {
        const uint nl = s_nlist;
        for (uint e = 0u; e < nl; ++e) {
            const uint2 r = args.node_range[s_list[e] * args.num_ranges + uint(L)];
            for (uint g = r.x; g < r.y; ++g) {
                const float4 gs = args.gaussians[g];
                const float3 X = gs.xyz - p.C;
                float s[GTO_INT_NS];
                gto_int_boys_scaled<L>(s, gs.w * dot(X, X), gs.w);
                acc = gto_int_fixed_add(acc, gto_int_hermite<L>(X, s, args.coeff_base + g * nh));
            }
        }
    }
    if (p.valid) args.accum[(group - args.accum_first) * GROUP + lid] = acc;
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
        args.out[group * GROUP + lid] = V;
    }
}

[shader("compute")] [numthreads(GROUP, 1, 1)] void grid_resolve(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { resolve_pass<GRID>(gid.x, lid); }
[shader("compute")] [numthreads(GROUP, 1, 1)] void xyz_resolve(uint3 gid : SV_GroupID, uint lid : SV_GroupIndex) { resolve_pass<POINTS>(gid.x, lid); }

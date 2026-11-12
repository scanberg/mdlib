// md_render_mol.slang -- kernels of md_render_mol_update (md_render_mol.c).
//
// Small molecules (solvent, ions) diffuse: clusters made of them once fall
// apart over a trajectory. So their part of the cluster order is rebuilt from
// a uniform grid when the caller asks (MD_RENDER_MOL_UPDATE_REGROUP), in O(n):
//
//     clear    the bounds and the per-cell counts
//     bounds   of the molecules' first atoms (one atomic per wave and axis)
//     dims     the cell edge: about CELL_TARGET atoms per cell at the bounds'
//              mean density, at most `cells` cells (the capacity the plan
//              allocated, so the host never reads back how many there are)
//     count    each molecule's cell, and its place in it: the atomics return
//              the atoms and bonds counted before it. Molecules in file order
//              are often neighbours, so a group adds a run of molecules with
//              the same cell in one atomic.
//     scan     exclusive, of the per-cell counts (blocks, sums, add)
//     write    each molecule's atoms and bonds into the tails of the
//              topology's arrays, at its cell's offset plus its place: the
//              bonds between the points its atoms now are
//     offsets  the cells' cluster offsets
//     mixed    the points of the other bonds touching small molecules
//
// And on every update, whatever else it does:
//
//     gather   the copy that is drawn: per point (cluster order), its atom's
//              position and radius as a float4, and its velocity. The radius
//              is a scattered fetch like the position; when neither the order
//              nor the radii changed it is the old copy's w instead, read in
//              order. A marker says the w are valid: written after a full
//              gather, for this topology.
//
// Every pass reads the molecules in their own order (the order of the small
// tables, close to the file order of the atoms): only the write lands
// scattered. The order within a cell is the atomics', not reproducible.
//
// A cell is one cluster. The grid spans the bounds, so a far outlier makes
// the cells larger everywhere (still correct, only looser clusters).
//
// Layout rules of md_gpu.slang apply to Args: no 3-vectors, pointers at 8.

#include "../md_gpu.slang"

static const float CELL_TARGET = 64.0;

struct Args {
    uint* position;             // float3 per atom, `stride` bytes apart
    uint* velocity;             // packed half3 per atom, `vel_stride` bytes apart
    uint* radius_override;      // float per atom
    uint* points;               // float4 per point: the copy
    uint* vel_out;              // uint2 per point
    uint* params;               // 16 words: ordered bounds, then the grid
    uint* mol_slot;             // uint4 per molecule: cell, its atoms' and bonds' place in it
    uint* cnt;                  // uint4 per cell, [cells + 1]: atoms, bonds
    uint* sums;                 // uint4 per scan block
    uint* types;                // uint16 per atom: particle type
    uint* type_radius;          // float per particle type (the topology's)
    uint* mol_atom_offset;      // [n_mol + 1]
    uint* mol_atoms;            // atom indices
    uint* mol_bond_offset;      // [n_mol + 1]
    uint* mol_bond_local;       // per bond: its atoms' places in the molecule, k0 | k1 << 16
    uint* mol_bond_ids;
    uint* mixed_pos;            // per mixed bond: its place in bond_pairs
    uint* mixed_atoms;          // uint2 per mixed bond
    uint* atom_order;           // the topology's: the tail from static_atoms on
    uint* atom_point;           // the topology's: atom -> point
    uint* atom_cluster_offset;  // the topology's: from static_atom_clusters on
    uint* bond_pairs;           // the topology's (the covalent class first): from static_bonds on
    uint* bond_ids;
    uint* bond_cluster_offset;  // the covalent class's: from static_bond_clusters on
    uint* marker;               // 4 words: whose radii the copy's w holds (see sol_gather)
    uint* frame0;               // frames: the slots combined (planes x, y, z, plane_words apart)
    uint* frame1;
    uint* frame2;
    uint* frame3;
    uint* fpar0;                // and their parameters (see FP_*)
    uint* fpar1;
    uint* fpar2;
    uint* fpar3;
    uint* qsrc;                 // reorder, permute: planes from; regroup: the uploaded frame's
    uint* qdst;                 // planes to
    uint* qpar;                 // regroup: the uploaded frame's parameters
    uint* src_point;            // regroup with resident frames: each small point's old point
    uint  stride;
    uint  n_mol;
    uint  cells;
    uint  n_atoms;              // atoms in small molecules
    uint  static_atoms;
    uint  static_atom_clusters;
    uint  static_bonds;
    uint  static_bond_clusters;
    uint  has_bonds;
    uint  mode;                 // scan: 0 blocks, 1 the block sums, 2 add
    uint  n;                    // scan: elements
    uint  atom_count;           // radii, gather
    uint  type_count;
    uint  vel_stride;
    uint  has_velocity;
    uint  has_radius;           // radius_override given
    uint  mixed_count;
    uint  full;                 // gather: the radii from their source (the order or the radii changed)
    uint  ident0;               // the topology the copy belongs to: its block address and version
    uint  ident1;
    uint  ident2;
    uint  _pad;
    uint  plane_words;          // frames: words from one plane to the next
    uint  frame_count;          // fr_gather: 1, 2 or 4
    uint  quantized;            // regroup: positions from qsrc, not `position`
    uint  has_src;              // sol_write: record src_point
    uint  has_fvel;             // fr_gather: velocities from dw
    uint  _pad2;
    float w0, w1, w2, w3;       // fr_gather: weights of the frames
    float dw0, dw1, dw2, dw3;   // and of their derivative
};
MD_SHADER_ARGS(Args);

// params
static const uint P_MIN   = 0u;    // 3 ordered words
static const uint P_MAX   = 3u;    // 3 ordered words
static const uint P_LO    = 6u;    // 3 floats
static const uint P_INV   = 9u;    // 1 / edge
static const uint P_DIM   = 10u;   // 3 uints
static const uint P_EDGE  = 13u;

static const uint MARKER = 0x6D64726Du;   // the copy's radii are valid

uint grid_threads(uint n) { return 64u * min((n + 63u) / 64u, 65535u); }
uint group_start(uint3 tid, uint gi) { return tid.x - gi; }

// Floats as uints in the same order, for atomic min and max
uint  ordered(float f)  { uint b = asuint(f); return (b & 0x80000000u) != 0u ? ~b : (b | 0x80000000u); }
float unordered(uint u) { return asfloat((u & 0x80000000u) != 0u ? (u & 0x7FFFFFFFu) : ~u); }

// ---- Frames: 16-bit planes, two values to a word ------------------------------
// Parameters of a frame, 32 words: origin, scale, the cell's columns, the
// rows of its inverse (zero for an axis that is not periodic), periodic.

static const uint FP_ORIGIN   = 0u;
static const uint FP_SCALE    = 3u;
static const uint FP_CELL     = 6u;
static const uint FP_INV      = 15u;
static const uint FP_PERIODIC = 24u;

float3 ld3(uint* p, uint o) { return float3(asfloat(p[o]), asfloat(p[o + 1u]), asfloat(p[o + 2u])); }

uint half16(uint* plane, uint i) { return (plane[i >> 1u] >> ((i & 1u) * 16u)) & 0xFFFFu; }

float3 dequant(uint* planes, uint pw, uint* par, uint i) {
    float3 q = float3(float(half16(planes, i)), float(half16(planes + pw, i)), float(half16(planes + 2u * pw, i)));
    return ld3(par, FP_ORIGIN) + q * ld3(par, FP_SCALE);
}

// d at its minimum image, by the cell's lattice: right while the true
// displacement is well inside the cell (between frames it is)
float3 min_image(uint* par, float3 d) {
    if (par[FP_PERIODIC] == 0u) return d;
    float3 n = round(float3(dot(ld3(par, FP_INV), d), dot(ld3(par, FP_INV + 3u), d), dot(ld3(par, FP_INV + 6u), d)));
    return d - ld3(par, FP_CELL) * n.x - ld3(par, FP_CELL + 3u) * n.y - ld3(par, FP_CELL + 6u) * n.z;
}

float3 first_pos(Args a, uint m) {
    uint i = a.mol_atoms[a.mol_atom_offset[m]];
    if (a.quantized != 0u) return dequant(a.qsrc, a.plane_words, a.qpar, i);
    uint* p = a.position + i * (a.stride >> 2);
    return float3(asfloat(p[0]), asfloat(p[1]), asfloat(p[2]));
}

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_clear(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = max(4u * (a.cells + 1u), 16u);
    for (uint i = tid.x; i < n; i += grid_threads(n)) {
        if (i < 4u * (a.cells + 1u)) a.cnt[i] = 0u;
        if (i < 16u) a.params[i] = i < P_MAX ? 0xFFFFFFFFu : 0u;
    }
}

// Every lane takes part in the wave operations (idle ones with the identity)
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_bounds(uint3 tid : SV_DispatchThreadID, uint gi : SV_GroupIndex) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint i0 = group_start(tid, gi); i0 < n; i0 += grid_threads(n)) {
        uint m = i0 + gi;
        bool valid = m < n;
        uint3 lo = uint3(0xFFFFFFFFu), hi = uint3(0u);
        if (valid) {
            float3 p = first_pos(a, m);
            lo = uint3(ordered(p.x), ordered(p.y), ordered(p.z));
            hi = lo;
        }
        uint3 wlo = uint3(WaveActiveMin(lo.x), WaveActiveMin(lo.y), WaveActiveMin(lo.z));
        uint3 whi = uint3(WaveActiveMax(hi.x), WaveActiveMax(hi.y), WaveActiveMax(hi.z));
        if (WaveIsFirstLane()) {
            ((Atomic<uint>*)(a.params + P_MIN + 0u))->min(wlo.x);
            ((Atomic<uint>*)(a.params + P_MIN + 1u))->min(wlo.y);
            ((Atomic<uint>*)(a.params + P_MIN + 2u))->min(wlo.z);
            ((Atomic<uint>*)(a.params + P_MAX + 0u))->max(whi.x);
            ((Atomic<uint>*)(a.params + P_MAX + 1u))->max(whi.y);
            ((Atomic<uint>*)(a.params + P_MAX + 2u))->max(whi.z);
        }
    }
}

[shader("compute")]
[numthreads(1, 1, 1)]
void sol_dims(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    float3 lo  = float3(unordered(a.params[P_MIN]), unordered(a.params[P_MIN + 1u]), unordered(a.params[P_MIN + 2u]));
    float3 hi  = float3(unordered(a.params[P_MAX]), unordered(a.params[P_MAX + 1u]), unordered(a.params[P_MAX + 2u]));
    float3 ext = max(hi - lo, float3(1.0e-3));
    float  edge = pow(CELL_TARGET * ext.x * ext.y * ext.z / float(max(a.n_atoms, 1u)), 1.0 / 3.0);
    edge = max(edge, 1.0e-3);
    uint3 d = uint3(1u);
    for (uint it = 0u; it < 256u; ++it) {
        d = uint3(min(floor(ext / edge), float3(1.0e6))) + 1u;
        if (float(d.x) * float(d.y) * float(d.z) <= float(a.cells) && d.x * d.y * d.z <= a.cells) break;
        edge *= 1.1;
        d = uint3(1u);
    }
    a.params[P_LO + 0u] = asuint(lo.x);
    a.params[P_LO + 1u] = asuint(lo.y);
    a.params[P_LO + 2u] = asuint(lo.z);
    a.params[P_INV]     = asuint(1.0 / edge);
    a.params[P_DIM + 0u] = d.x;
    a.params[P_DIM + 1u] = d.y;
    a.params[P_DIM + 2u] = d.z;
    a.params[P_EDGE]    = asuint(edge);
}

uint cell_of(Args a, float3 p) {
    float3 lo  = float3(asfloat(a.params[P_LO]), asfloat(a.params[P_LO + 1u]), asfloat(a.params[P_LO + 2u]));
    float  inv = asfloat(a.params[P_INV]);
    uint3  d   = uint3(a.params[P_DIM], a.params[P_DIM + 1u], a.params[P_DIM + 2u]);
    float3 q   = clamp(floor((p - lo) * inv), float3(0.0), float3(d - 1u));   // NaN: max() takes the 0
    uint3  c   = uint3(q);
    return min((c.z * d.y + c.y) * d.x + c.x, a.cells - 1u);
}

// One group adds a run of equal cells with one atomic: each lane's place is
// the run's base plus its exclusive prefix within the run. Groupshared, not
// wave operations: the same on every wave size.
groupshared uint  g_key[64];
groupshared uint2 g_inc[64];
groupshared uint2 g_exc[64];
groupshared uint2 g_base[64];

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_count(uint3 tid : SV_DispatchThreadID, uint gi : SV_GroupIndex) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint i0 = group_start(tid, gi); i0 < n; i0 += grid_threads(n)) {
        uint m = i0 + gi;
        bool valid = m < n;
        uint key = 0xFFFFFFFFu;
        uint2 v = uint2(0u);
        if (valid) {
            key = cell_of(a, first_pos(a, m));
            v.x = a.mol_atom_offset[m + 1u] - a.mol_atom_offset[m];
            if (a.has_bonds != 0u) v.y = a.mol_bond_offset[m + 1u] - a.mol_bond_offset[m];
        }
        g_key[gi] = key;
        g_inc[gi] = v;
        GroupMemoryBarrierWithGroupSync();
        for (uint off = 1u; off < 64u; off <<= 1u) {
            uint2 x = gi >= off ? g_inc[gi - off] : uint2(0u);
            GroupMemoryBarrierWithGroupSync();
            g_inc[gi] += x;
            GroupMemoryBarrierWithGroupSync();
        }
        uint2 exc = g_inc[gi] - v;
        g_exc[gi] = exc;
        GroupMemoryBarrierWithGroupSync();
        uint h = gi;                                    // the run's first lane
        while (h > 0u && g_key[h - 1u] == key) h -= 1u;
        if (valid && h == gi) {
            uint e = gi;                                // its last
            while (e + 1u < 64u && g_key[e + 1u] == key) e += 1u;
            uint2 total = g_inc[e] - exc;
            uint2 base;
            base.x = ((Atomic<uint>*)(a.cnt + 4u * key + 0u))->add(total.x);
            base.y = total.y > 0u ? ((Atomic<uint>*)(a.cnt + 4u * key + 1u))->add(total.y) : 0u;
            g_base[gi] = base;
        }
        GroupMemoryBarrierWithGroupSync();
        if (valid) {
            uint2 place = g_base[h] + (exc - g_exc[h]);
            uint* s = a.mol_slot + 4u * m;
            s[0] = key; s[1] = place.x; s[2] = place.y; s[3] = 0u;
        }
        GroupMemoryBarrierWithGroupSync();              // the arrays are free for the next run
    }
}

// ---- Exclusive scan of uint4 elements (xyz; w is left 0) ----------------------
// A block is 1024 elements: 256 threads, 4 each.

static const uint SCAN_THREADS = 256u;
static const uint SCAN_BLOCK   = 1024u;

groupshared uint4 sh_scan[256];

uint4 load4(uint* p, uint e, uint n) {
    if (e >= n) return uint4(0u);
    return uint4(p[4u * e], p[4u * e + 1u], p[4u * e + 2u], 0u);
}

void store4(uint* p, uint e, uint n, uint4 v) {
    if (e >= n) return;
    p[4u * e] = v.x; p[4u * e + 1u] = v.y; p[4u * e + 2u] = v.z; p[4u * e + 3u] = 0u;
}

// Elements [first, first + SCAN_BLOCK) of p (those below n) replaced by carry
// plus their exclusive prefix; returns the block's total. All threads call it.
uint4 scan_block(uint* p, uint first, uint n, uint t, uint4 carry) {
    uint e = first + 4u * t;
    uint4 v0 = load4(p, e, n), v1 = load4(p, e + 1u, n), v2 = load4(p, e + 2u, n), v3 = load4(p, e + 3u, n);
    uint4 s = v0 + v1 + v2 + v3;
    sh_scan[t] = s;
    GroupMemoryBarrierWithGroupSync();
    for (uint off = 1u; off < SCAN_THREADS; off <<= 1u) {
        uint4 x = t >= off ? sh_scan[t - off] : uint4(0u);
        GroupMemoryBarrierWithGroupSync();
        sh_scan[t] += x;
        GroupMemoryBarrierWithGroupSync();
    }
    uint4 excl  = carry + sh_scan[t] - s;
    uint4 total = sh_scan[SCAN_THREADS - 1u];
    GroupMemoryBarrierWithGroupSync();      // sh_scan is free for the next block
    store4(p, e, n, excl);
    store4(p, e + 1u, n, excl + v0);
    store4(p, e + 2u, n, excl + v0 + v1);
    store4(p, e + 3u, n, excl + v0 + v1 + v2);
    return total;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void sol_scan(uint3 gid : SV_GroupID, uint gi : SV_GroupIndex) {
    Args a = MD_ARGS;
    if (a.mode == 0u) {             // each block on its own; its total to sums
        uint4 total = scan_block(a.cnt, gid.x * SCAN_BLOCK, a.n, gi, uint4(0u));
        if (gi == 0u) { uint* s = a.sums + 4u * gid.x; s[0] = total.x; s[1] = total.y; s[2] = total.z; s[3] = 0u; }
    } else if (a.mode == 1u) {      // one group: the block sums, block after block
        uint nb = (a.n + SCAN_BLOCK - 1u) / SCAN_BLOCK;
        uint4 carry = uint4(0u);
        for (uint first = 0u; first < nb; first += SCAN_BLOCK) carry += scan_block(a.sums, first, nb, gi, carry);
    } else {                        // the block's offset into each of its elements
        uint b = gid.x;
        if (b == 0u) return;
        uint4 off = load4(a.sums, b, b + 1u);
        for (uint k = 0u; k < 4u; ++k) {
            uint e = b * SCAN_BLOCK + 4u * gi + k;
            if (e < a.n) store4(a.cnt, e, a.n, load4(a.cnt, e, a.n) + off);
        }
    }
}

// ---- Write, offsets, radii --------------------------------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_write(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint m = tid.x; m < n; m += grid_threads(n)) {
        uint* s = a.mol_slot + 4u * m;
        uint c = s[0];
        uint p0 = a.static_atoms + a.cnt[4u * c] + s[1];      // the molecule's first point
        uint at = p0;
        for (uint k = a.mol_atom_offset[m]; k < a.mol_atom_offset[m + 1u]; ++k) {
            uint atom = a.mol_atoms[k];
            if (a.has_src != 0u) a.src_point[at - a.static_atoms] = a.atom_point[atom];
            a.atom_order[at] = atom;
            a.atom_point[atom] = at;
            at += 1u;
        }
        if (a.has_bonds != 0u) {
            uint bt = a.static_bonds + a.cnt[4u * c + 1u] + s[2];
            for (uint k = a.mol_bond_offset[m]; k < a.mol_bond_offset[m + 1u]; ++k) {
                uint l = a.mol_bond_local[k];
                a.bond_pairs[2u * bt]      = p0 + (l & 0xFFFFu);
                a.bond_pairs[2u * bt + 1u] = p0 + (l >> 16u);
                a.bond_ids[bt]             = a.mol_bond_ids[k];
                bt += 1u;
            }
        }
    }
}

// cells + 1 entries: the last is the totals (the exclusive scan's end)
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_offsets(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.cells + 1u;
    for (uint c = tid.x; c < n; c += grid_threads(n)) {
        a.atom_cluster_offset[a.static_atom_clusters + c] = a.static_atoms + a.cnt[4u * c];
        if (a.has_bonds != 0u) a.bond_cluster_offset[a.static_bond_clusters + c] = a.static_bonds + a.cnt[4u * c + 1u];
    }
}


// After a regroup: the other bonds touching small molecules, between the
// points their atoms now are
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_mixed(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.mixed_count;
    for (uint i = tid.x; i < n; i += grid_threads(n)) {
        uint pos = a.mixed_pos[i];
        a.bond_pairs[2u * pos]      = a.atom_point[a.mixed_atoms[2u * i]];
        a.bond_pairs[2u * pos + 1u] = a.atom_point[a.mixed_atoms[2u * i + 1u]];
    }
}

bool radii_valid(Args a) {
    return a.full == 0u && a.marker[0] == MARKER && a.marker[1] == a.ident0 && a.marker[2] == a.ident1 && a.marker[3] == a.ident2;
}

// An atom's radius from its source: the caller's array, or its type's
uint radius_of(Args a, uint atom) {
    if (a.has_radius != 0u) return a.radius_override[atom];
    uint t = (a.types[atom >> 1u] >> ((atom & 1u) * 16u)) & 0xFFFFu;
    return t < a.type_count ? a.type_radius[t] : 0u;
}

// The copy that is drawn: point p is atom atom_order[p]. Writes in order,
// reads where the atoms are (in runs: the atoms of a residue or a molecule
// are neighbours in both orders).
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_gather(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    bool reuse = radii_valid(a);
    for (uint p = tid.x; p < n; p += grid_threads(n)) {
        uint atom = a.atom_order[p];
        uint* q = a.position + atom * (a.stride >> 2u);
        uint r = reuse ? a.points[4u * p + 3u] : radius_of(a, atom);
        uint4 v = uint4(q[0], q[1], q[2], r);
        ((uint4*)a.points)[p] = v;
        if (a.has_velocity != 0u) {
            uint* w = a.velocity + atom * (a.vel_stride >> 2u);
            ((uint2*)a.vel_out)[p] = uint2(w[0], w[1]);
        }
    }
}

// After a full gather (its own launch: the gather's threads read the marker)
[shader("compute")]
[numthreads(1, 1, 1)]
void sol_mark(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    a.marker[0] = MARKER;
    a.marker[1] = a.ident0;
    a.marker[2] = a.ident1;
    a.marker[3] = a.ident2;
}

// ---- Frames --------------------------------------------------------------------

// The copy that is drawn, from resident frames in point order: every read in
// order (the atom only for radii from their source)
[shader("compute")]
[numthreads(64, 1, 1)]
void fr_gather(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    uint pw = a.plane_words;
    bool reuse = radii_valid(a);
    for (uint p = tid.x; p < n; p += grid_threads(n)) {
        float3 pos, vel = float3(0.0);
        float3 x0 = dequant(a.frame0, pw, a.fpar0, p);
        if (a.frame_count == 1u) {
            pos = x0;
        } else if (a.frame_count == 2u) {
            float3 x1 = x0 + min_image(a.fpar0, dequant(a.frame1, pw, a.fpar1, p) - x0);
            pos = a.w0 * x0 + a.w1 * x1;
            vel = a.dw0 * x0 + a.dw1 * x1;
        } else {
            float3 x1 = dequant(a.frame1, pw, a.fpar1, p);
            x0 = x1 + min_image(a.fpar1, x0 - x1);
            float3 x2 = x1 + min_image(a.fpar1, dequant(a.frame2, pw, a.fpar2, p) - x1);
            float3 x3 = x2 + min_image(a.fpar1, dequant(a.frame3, pw, a.fpar3, p) - x2);
            pos = a.w0 * x0 + a.w1 * x1 + a.w2 * x2 + a.w3 * x3;
            vel = a.dw0 * x0 + a.dw1 * x1 + a.dw2 * x2 + a.dw3 * x3;
        }
        uint r = reuse ? a.points[4u * p + 3u] : radius_of(a, a.atom_order[p]);
        ((uint4*)a.points)[p] = uint4(asuint(pos.x), asuint(pos.y), asuint(pos.z), r);
        if (a.has_fvel != 0u) ((uint2*)a.vel_out)[p] = uint2(f32tof16(vel.x) | (f32tof16(vel.y) << 16u), f32tof16(vel.z));
    }
}

// An uploaded frame (qsrc, atom order) into its slot (qdst, point order).
// A thread per word of the destination: two points, each from its atom.
[shader("compute")]
[numthreads(64, 1, 1)]
void fr_reorder(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count, words = (n + 1u) >> 1u, pw = a.plane_words;
    for (uint w = tid.x; w < words; w += grid_threads(words)) {
        uint i0 = a.atom_order[2u * w];
        uint i1 = 2u * w + 1u < n ? a.atom_order[2u * w + 1u] : i0;
        for (uint k = 0u; k < 3u; ++k) {
            uint* src = a.qsrc + k * pw;
            a.qdst[k * pw + w] = half16(src, i0) | (half16(src, i1) << 16u);
        }
    }
}

// After a regroup: a resident frame (qsrc, the old point order) into the new
// one (qdst). The static points stay; a small one was src_point[p - static].
[shader("compute")]
[numthreads(64, 1, 1)]
void fr_permute(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count, words = (n + 1u) >> 1u, pw = a.plane_words, st = a.static_atoms;
    for (uint w = tid.x; w < words; w += grid_threads(words)) {
        uint p0 = 2u * w, p1 = min(p0 + 1u, n - 1u);
        uint s0 = p0 < st ? p0 : a.src_point[p0 - st];
        uint s1 = p1 < st ? p1 : a.src_point[p1 - st];
        for (uint k = 0u; k < 3u; ++k) {
            uint* src = a.qsrc + k * pw;
            a.qdst[k * pw + w] = half16(src, s0) | (half16(src, s1) << 16u);
        }
    }
}

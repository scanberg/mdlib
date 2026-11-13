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
    uint* marker;               // 5 words: whose radii the copy's w holds, in which grouping (see sol_gather)
    uint* frame0;               // frames: the slots combined
    uint* frame1;
    uint* frame2;
    uint* frame3;
    uint* fpar0;                // and their parameters (see FP_*)
    uint* fpar1;
    uint* fpar2;
    uint* fpar3;
    uint* qsrc;                 // reorder, permute: the frame from; regroup: the uploaded frame
    uint* qdst;                 // the frame to
    uint* qpar;                 // regroup: the uploaded frame's parameters
    uint* src_point;            // regroup with resident frames: each small point's old point
    uint* resident;             // the resident atoms: bits, then each word's rank (has_resident)
    uint* regroups;             // the topology's: word 0 counts the regroups
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
    uint  frame_count;          // fr_gather: 1, 2 or 4
    uint  quantized;            // regroup: positions from qsrc, not `position`
    uint  has_src;              // sol_write: record src_point
    uint  has_fvel;             // fr_gather: velocities from dw
    uint  has_resident;         // 0: every atom is resident, its place in a frame is itself
    uint  resident_words;       // words of the bits (the ranks follow them)
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

// ---- Frames: x, y, z per point in 16 bits each, three words to two points -------
// Parameters of a frame, 32 words: origin, scale, the cell (a along x, b in
// the xy plane: ax, bx, by, cx, cy, cz) and the reciprocals of its diagonal
// (zero for an axis that is not periodic).

static const uint FP_ORIGIN   = 0u;
static const uint FP_SCALE    = 3u;
static const uint FP_CELL     = 6u;
static const uint FP_RCELL    = 12u;

float3 ld3(uint* p, uint o) { return float3(asfloat(p[o]), asfloat(p[o + 1u]), asfloat(p[o + 2u])); }

// Point i's three values: halves 3i .. 3i+2, in two words. Selects, not a
// branch: even and odd points are neighbours, so a branch split every wave
// in two (and a compiler need not turn it back into selects).
uint3 load_q(uint* q, uint i) {
    uint w = (3u * i) >> 1u;
    uint a = q[w], b = q[w + 1u];
    bool odd = (i & 1u) != 0u;
    return uint3(select(odd, a >> 16u, a & 0xFFFFu), select(odd, b & 0xFFFFu, a >> 16u), select(odd, b >> 16u, b & 0xFFFFu));
}

float3 dequant(uint* q, uint* par, uint i) {
    return ld3(par, FP_ORIGIN) + float3(load_q(q, i)) * ld3(par, FP_SCALE);
}

// Points 2j and 2j+1 (if it is one of n) as words 3j .. 3j+2
void store_pair(uint* q, uint j, uint3 a, uint3 b, bool has_b) {
    q[3u * j]      = a.x | (a.y << 16u);
    q[3u * j + 1u] = a.z | (b.x << 16u);
    if (has_b) q[3u * j + 2u] = b.y | (b.z << 16u);
}


// An atom's place in an uploaded frame: its rank among the resident atoms
uint upload_index(Args a, uint atom) {
    if (a.has_resident == 0u) return atom;
    uint w = atom >> 5u;
    return a.resident[a.resident_words + w] + countbits(a.resident[w] & ((1u << (atom & 31u)) - 1u));
}

float3 first_pos(Args a, uint m) {
    uint i = a.mol_atoms[a.mol_atom_offset[m]];
    if (a.quantized != 0u) return dequant(a.qsrc, a.qpar, upload_index(a, i));
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
    if (tid.x == 0u) a.regroups[0] = a.regroups[0] + 1u;  // every copy gathered before is in the old order
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
    return a.full == 0u && a.marker[0] == MARKER && a.marker[1] == a.ident0 && a.marker[2] == a.ident1 && a.marker[3] == a.ident2 &&
           a.marker[4] == a.regroups[0];
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
    a.marker[4] = a.regroups[0];
}

// ---- Frames --------------------------------------------------------------------

// A frame's origin and scale, and the reference cell, in registers: read
// through pointers inside the loop they were reloaded at every point (the
// stores to the copy may alias them, as far as the compiler knows)
struct FrameP { float3 o; float3 s; };
struct CellP  { float ax; float bx; float by; float cx; float cy; float cz; float3 r; };

FrameP frame_p(uint* par) { FrameP f; f.o = ld3(par, FP_ORIGIN); f.s = ld3(par, FP_SCALE); return f; }

CellP cell_p(uint* par) {
    CellP c;
    c.ax = asfloat(par[FP_CELL]);      c.bx = asfloat(par[FP_CELL + 1u]); c.by = asfloat(par[FP_CELL + 2u]);
    c.cx = asfloat(par[FP_CELL + 3u]); c.cy = asfloat(par[FP_CELL + 4u]); c.cz = asfloat(par[FP_CELL + 5u]);
    c.r  = ld3(par, FP_RCELL);
    return c;
}

float3 deq(uint* q, FrameP f, uint i) { return f.o + float3(load_q(q, i)) * f.s; }

// d at its minimum image: c, then b, then a taken off as often as they fit
// (as GROMACS does for a triclinic cell), right while the true displacement
// is well inside the cell (between frames it is). No branch: a reciprocal is
// zero where an axis is not periodic, so its step is too. Nine values, not
// the eighteen of a cell and its inverse: on Pascal those registers pushed
// the two-frame kernel past a tier.
float3 mi(CellP c, float3 d) {
    float n2 = round(d.z * c.r.z);
    d -= n2 * float3(c.cx, c.cy, c.cz);
    float n1 = round(d.y * c.r.y);
    d.xy -= n1 * float2(c.bx, c.by);
    d.x  -= round(d.x * c.r.x) * c.ax;
    return d;
}

// Before the frames' gather, when the copy's radii are not valid: each
// point's radius from its source, into the copy's w. fr_check sizes its
// launch (none when they are), so the gathers need not carry that path.
[shader("compute")]
[numthreads(1, 1, 1)]
void fr_check(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    uint* g = a.marker + 16u;                               // the indirect grid
    g[0] = radii_valid(a) || n == 0u ? 0u : min((n + 63u) / 64u, 65535u);
    g[1] = 1u;
    g[2] = 1u;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void fr_radii(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    for (uint p = tid.x; p < n; p += grid_threads(n)) a.points[4u * p + 3u] = radius_of(a, a.atom_order[p]);
}

// The copy that is drawn, from resident frames in point order: every read in
// order, the radius the copy's own (fr_radii saw to it). A kernel per number
// of frames (COUNT a constant): one kernel for all three kept every frame's
// parameters and the cell live, and on Pascal the registers that cost made
// the nearest frame twice the copy from positions.
[ForceInline]
void fr_gather_body(uint3 tid, const uint COUNT) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    const bool vel_out = COUNT > 1u && a.has_fvel != 0u;
    uint* points = a.points;
    const FrameP f0 = frame_p(a.fpar0);
    FrameP f1 = f0, f2 = f0, f3 = f0;
    CellP  cell;
    if (COUNT >= 2u) { f1 = frame_p(a.fpar1); cell = cell_p(COUNT == 4u ? a.fpar1 : a.fpar0); }
    if (COUNT == 4u) { f2 = frame_p(a.fpar2); f3 = frame_p(a.fpar3); }
    for (uint p = tid.x; p < n; p += grid_threads(n)) {
        float3 pos, vel = float3(0.0);
        float3 x0 = deq(a.frame0, f0, p);
        if (COUNT == 1u) {
            pos = x0;
        } else if (COUNT == 2u) {
            float3 x1 = x0 + mi(cell, deq(a.frame1, f1, p) - x0);
            pos = a.w0 * x0 + a.w1 * x1;
            vel = a.dw0 * x0 + a.dw1 * x1;
        } else {
            float3 x1 = deq(a.frame1, f1, p);
            x0 = x1 + mi(cell, x0 - x1);
            float3 x2 = x1 + mi(cell, deq(a.frame2, f2, p) - x1);
            float3 x3 = x2 + mi(cell, deq(a.frame3, f3, p) - x2);
            pos = a.w0 * x0 + a.w1 * x1 + a.w2 * x2 + a.w3 * x3;
            vel = a.dw0 * x0 + a.dw1 * x1 + a.dw2 * x2 + a.dw3 * x3;
        }
        uint r = points[4u * p + 3u];
        ((uint4*)points)[p] = uint4(asuint(pos.x), asuint(pos.y), asuint(pos.z), r);
        if (vel_out) ((uint2*)a.vel_out)[p] = uint2(f32tof16(vel.x) | (f32tof16(vel.y) << 16u), f32tof16(vel.z));
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void fr_gather1(uint3 tid : SV_DispatchThreadID) { fr_gather_body(tid, 1u); }

[shader("compute")]
[numthreads(64, 1, 1)]
void fr_gather2(uint3 tid : SV_DispatchThreadID) { fr_gather_body(tid, 2u); }

[shader("compute")]
[numthreads(64, 1, 1)]
void fr_gather4(uint3 tid : SV_DispatchThreadID) { fr_gather_body(tid, 4u); }

// An uploaded frame (qsrc: the resident atoms in atom order) into its slot
// (qdst, point order).
// A thread per two points of the destination, each read from its atom: one
// gather of two words per point.
[shader("compute")]
[numthreads(64, 1, 1)]
void fr_reorder(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count, pairs = (n + 1u) >> 1u;
    for (uint j = tid.x; j < pairs; j += grid_threads(pairs)) {
        bool has_b = 2u * j + 1u < n;
        uint3 qa = load_q(a.qsrc, upload_index(a, a.atom_order[2u * j]));
        uint3 qb = has_b ? load_q(a.qsrc, upload_index(a, a.atom_order[2u * j + 1u])) : uint3(0u);
        store_pair(a.qdst, j, qa, qb, has_b);
    }
}

// After a regroup: a resident frame (qsrc, the old point order) into the new
// one (qdst). The static points stay; a small one was src_point[p - static].
[shader("compute")]
[numthreads(64, 1, 1)]
void fr_permute(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count, pairs = (n + 1u) >> 1u, st = a.static_atoms;
    for (uint j = tid.x; j < pairs; j += grid_threads(pairs)) {
        uint p0 = 2u * j, p1 = p0 + 1u;
        bool has_b = p1 < n;
        uint3 qa = load_q(a.qsrc, p0 < st ? p0 : a.src_point[p0 - st]);
        uint3 qb = has_b ? load_q(a.qsrc, p1 < st ? p1 : a.src_point[p1 - st]) : uint3(0u);
        store_pair(a.qdst, j, qa, qb, has_b);
    }
}

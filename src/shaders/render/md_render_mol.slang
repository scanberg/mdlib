// md_render_mol.slang -- kernels of md_render_mol_update (md_render_mol.c).
//
// Small molecules (solvent, ions) diffuse: clusters made of them once fall
// apart over a trajectory. So their part of the cluster order is rebuilt from
// a uniform grid whenever the positions change, in O(n):
//
//     clear    the bounds and the per-cell counts
//     bounds   of the molecules' first atoms (one atomic per wave and axis)
//     dims     the cell edge: about CELL_TARGET atoms per cell at the bounds'
//              mean density, at most `cells` cells (the capacity the plan
//              allocated, so the host never reads back how many there are)
//     count    molecules, atoms and bonds per cell
//     scan     exclusive, of the three counts at once (blocks, sums, add)
//     scatter  each molecule into its cell's list
//     write    one thread per cell: its molecules sorted (the order is then
//              a function of the positions alone), their atoms and bonds into
//              the tails of the topology's arrays, the cell's cluster offsets
//
// A cell is one cluster. The grid spans the bounds, so a far outlier makes
// the cells larger everywhere (still correct, only looser clusters).
//
// Layout rules of md_gpu.slang apply to Args: no 3-vectors, pointers at 8.

#include "../md_gpu.slang"

static const float CELL_TARGET = 48.0;

struct Args {
    uint* position;             // float3 per atom, `stride` bytes apart
    uint* params;               // 16 words: ordered bounds, then the grid
    uint* mol_cell;             // [n_mol]
    uint* cnt;                  // uint4 per cell, [cells + 1]: molecules, atoms, bonds, fill
    uint* sums;                 // uint4 per scan block
    uint* mol_list;             // [n_mol]: molecules by cell
    uint* mol_atom_offset;      // [n_mol + 1]
    uint* mol_atoms;            // atom indices
    uint* mol_bond_offset;      // [n_mol + 1]
    uint* mol_bond_pairs;       // uint2 per bond
    uint* mol_bond_ids;
    uint* atom_order;           // the topology's: the tail from static_atoms on
    uint* atom_cluster_offset;  // the topology's: from static_atom_clusters on
    uint* bond_pairs;           // the topology's covalent class: from static_bonds on
    uint* bond_ids;
    uint* bond_cluster_offset;  // the covalent class's: from static_bond_clusters on
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
    uint  _pad;
};
MD_SHADER_ARGS(Args);

// params
static const uint P_MIN   = 0u;    // 3 ordered words
static const uint P_MAX   = 3u;    // 3 ordered words
static const uint P_LO    = 6u;    // 3 floats
static const uint P_INV   = 9u;    // 1 / edge
static const uint P_DIM   = 10u;   // 3 uints
static const uint P_EDGE  = 13u;

uint grid_threads(uint n) { return 64u * min((n + 63u) / 64u, 65535u); }
uint group_start(uint3 tid, uint gi) { return tid.x - gi; }

// Floats as uints in the same order, for atomic min and max
uint  ordered(float f)  { uint b = asuint(f); return (b & 0x80000000u) != 0u ? ~b : (b | 0x80000000u); }
float unordered(uint u) { return asfloat((u & 0x80000000u) != 0u ? (u & 0x7FFFFFFFu) : ~u); }

float3 first_pos(Args a, uint m) {
    uint i = a.mol_atoms[a.mol_atom_offset[m]];
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

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_count(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint m = tid.x; m < n; m += grid_threads(n)) {
        uint c = cell_of(a, first_pos(a, m));
        a.mol_cell[m] = c;
        ((Atomic<uint>*)(a.cnt + 4u * c + 0u))->add(1u);
        ((Atomic<uint>*)(a.cnt + 4u * c + 1u))->add(a.mol_atom_offset[m + 1u] - a.mol_atom_offset[m]);
        if (a.has_bonds != 0u) ((Atomic<uint>*)(a.cnt + 4u * c + 2u))->add(a.mol_bond_offset[m + 1u] - a.mol_bond_offset[m]);
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

// ---- Scatter and write -----------------------------------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_scatter(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint m = tid.x; m < n; m += grid_threads(n)) {
        uint c = a.mol_cell[m];
        uint k = ((Atomic<uint>*)(a.cnt + 4u * c + 3u))->add(1u);
        a.mol_list[a.cnt[4u * c] + k] = m;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void sol_write(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.cells;
    for (uint c = tid.x; c < n; c += grid_threads(n)) {
        uint mo = a.cnt[4u * c], ao = a.cnt[4u * c + 1u], bo = a.cnt[4u * c + 2u];
        uint nm = a.cnt[4u * (c + 1u)] - mo;
        // Sorted, so that the order is the same for the same positions. A
        // degenerate cell (everything in one place) stays unsorted rather
        // than quadratic.
        if (nm <= 256u) {
            for (uint i = 1u; i < nm; ++i) {
                uint x = a.mol_list[mo + i];
                uint j = i;
                while (j > 0u && a.mol_list[mo + j - 1u] > x) { a.mol_list[mo + j] = a.mol_list[mo + j - 1u]; j -= 1u; }
                a.mol_list[mo + j] = x;
            }
        }
        uint at = a.static_atoms + ao;
        uint bt = a.static_bonds + bo;
        for (uint i = 0u; i < nm; ++i) {
            uint m = a.mol_list[mo + i];
            for (uint k = a.mol_atom_offset[m]; k < a.mol_atom_offset[m + 1u]; ++k) a.atom_order[at++] = a.mol_atoms[k];
            if (a.has_bonds != 0u) {
                for (uint k = a.mol_bond_offset[m]; k < a.mol_bond_offset[m + 1u]; ++k) {
                    a.bond_pairs[2u * bt]      = a.mol_bond_pairs[2u * k];
                    a.bond_pairs[2u * bt + 1u] = a.mol_bond_pairs[2u * k + 1u];
                    a.bond_ids[bt]             = a.mol_bond_ids[k];
                    bt += 1u;
                }
            }
        }
        a.atom_cluster_offset[a.static_atom_clusters + c] = a.static_atoms + ao;
        if (a.has_bonds != 0u) a.bond_cluster_offset[a.static_bond_clusters + c] = a.static_bonds + bo;
        if (c + 1u == n) {
            a.atom_cluster_offset[a.static_atom_clusters + n] = a.static_atoms + a.cnt[4u * n + 1u];
            if (a.has_bonds != 0u) a.bond_cluster_offset[a.static_bond_clusters + n] = a.static_bonds + a.cnt[4u * n + 2u];
        }
    }
}

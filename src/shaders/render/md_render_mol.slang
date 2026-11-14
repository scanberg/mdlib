// md_render_mol.slang -- kernels of md_render_mol_update (md_render_mol.c).
//
// Small molecules (solvent, ions) diffuse: clusters made of them once fall
// apart over a trajectory. So their part of the cluster order is rebuilt from
// a uniform grid when the caller asks (MD_RENDER_MOL_UPDATE_REGROUP), in O(n):
//
//     clear    the bounds and the per-cell counts
//     bounds   of the molecules' first atoms (one atomic per wave and axis),
//              each staged in mol_slot for the count
//     dims     the cell edge: about CELL_TARGET atoms per cell at the bounds'
//              mean density, at most `cells` cells (the capacity the plan
//              allocated, so the host never reads back how many there are)
//     count    each molecule's cell, and its place in it: the atomics return
//              the atoms, bonds and molecules counted before it; a group adds
//              a run of molecules with the same cell in one atomic
//     scan     exclusive, of the per-cell counts (blocks, sums, add)
//     write    each molecule's atoms, their types and its bonds into the
//              tails of the topology's arrays, at its cell's offset plus its
//              place (the bonds between the points its atoms now are), and
//              its record into the next grouping's order
//     offsets  the cells' cluster offsets; the regroups counted
//     mixed    the points of the other bonds touching small molecules
//
// And on every update, whatever else it does:
//
//     gather   the copy that is drawn: per point (cluster order), its atom's
//              position and a word (x, y, z, w), and its velocity. w is the
//              point's type and flags (POINT_TYPE_* below) from the
//              topology's point_type, in order: md_render takes the radius
//              from the types' table (the topology's type_radius), so a
//              change of radii is a change of that table and nothing else.
//              With the caller's per-atom radii instead, w is the radius,
//              fetched (scattered, as the position) on every update.
//
// The small molecules are records (first atom, its bonds' place among the
// small bonds, kind, first point)
// in the last grouping's order, and every regroup pass goes over them in that
// order: after a short move the molecules of neighbouring slots are still in
// the same or neighbouring cells, so the count's atomics come in runs and the
// write's stores land near each other. Only the bounds read scattered (the
// positions are in file order). The first regroup after the plan starts from
// molecule order and scatters once. The order within a cell is the atomics',
// not reproducible.
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
    uint* points;               // float4 per point (x, y, z, w: the type and flags, or a radius): the copy
    uint* vel_out;              // uint2 per point
    uint* params;               // 16 words: ordered bounds, then the grid
    uint* mol_slot;             // uint4 per slot (a molecule, the last grouping's order): its first
                                // atom's position (bounds), then its cell and its atoms', bonds' and
                                // molecules' place in it (count)
    uint* cnt;                  // uint4 per cell, [cells + 1]: atoms, bonds, molecules
    uint* sums;                 // uint4 per scan block
    uint* types;                // uint16 per atom: particle type
    uint* type_radius;          // float per particle type (the topology's)
    uint* kinds;                // uint4 per kind of small molecule: atoms, bonds, first of each below
    uint* kind_atom;            // per atom of a kind: its atom less the molecule's first
    uint* kind_type;            // and its particle type
    uint* kind_bond_local;      // per bond of a kind: its atoms' places in the molecule, k0 | k1 << 16
    uint* ids_old;              // the small bonds' ids of the last grouping (copied aside by the regroup)
    uint* recs0;                // uint4 per small molecule, in cluster order: first atom, its bonds' place
                                // among the small bonds,
    uint* recs1;                // kind (bit 31: has mixed bonds), its first point. Two: one per parity
                                // of regroups[0], the other the regroup's output
    uint* mixed_pos;            // per mixed bond: its place in bond_pairs
    uint* mixed_atoms;          // uint2 per mixed bond
    uint* atom_order;           // the topology's: the tail from static_atoms on
    uint* mixed_point;          // the topology's: atom -> point, for the static atoms and the atoms of
                                // molecules with mixed bonds (sol_mixed's)
    uint* atom_cluster_offset;  // the topology's: from static_atom_clusters on
    uint* bond_pairs;           // the topology's (the covalent class first): from static_bonds on
    uint* bond_ids;
    uint* bond_cluster_offset;  // the covalent class's: from static_bond_clusters on
    uint* point_type;           // the topology's: per point its type and flags (POINT_TYPE_*), in point order
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
    uint* bb_index;             // the backbone (the cartoon): per control point its CA's point | flags
    uint* bb_segment;           // its segment of sys->protein_backbone
    uint* bb_carbonyl;          // uint2: its C's and O's points
    uint* bb_chain;             // each piece's first control point, [bb_chains + 1]
    uint* bb_axis;              // float4 per control point: the ribbon's axis (md_render's)
    uint* bb_flip;              // per control point: the next one's axis turned against it (bit 0);
                                // then the stamp of the update that wrote them
    uint* bb_ss;                // int32 per segment: md_secondary_structure_t (has_ss)
    uint* bb_out;               // float4 per control point: the profile (md_render's)
    uint* na_index;             // the nucleic band: per control point its C4''s point | flags
    uint* na_face;              // the point it faces (UINT32_MAX for none)
    uint* na_axis;              // float4 per control point: the band's axis (md_render's)
    uint* der_offset;           // the derived points: each one's members, [der_count + 1]
    uint* der_member;           // per member its atom
    uint* der_weight;           // float per member
    uint* der_radius;           // float per derived point: its w
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
    uint  atom_count;           // the points (gather)
    uint  type_count;
    uint  vel_stride;
    uint  has_velocity;
    uint  has_radius;           // radius_override given
    uint  mixed_count;
    uint  fr_threads;           // fr_gather: the threads launched (its stride)
    uint  frame_count;          // fr_gather: 1, 2 or 4
    uint  quantized;            // regroup: positions from qsrc, not `position`
    uint  has_src;              // sol_write: record src_point
    uint  has_fvel;             // fr_gather: velocities from dw
    uint  has_resident;         // 0: every atom is resident, its place in a frame is itself
    uint  resident_words;       // words of the bits (the ranks follow them)
    float w0, w1, w2, w3;       // fr_gather: weights of the frames
    float dw0, dw1, dw2, dw3;   // and of their derivative
    uint  bb_count;             // control points
    uint  bb_chains;            // pieces
    uint  bb_history;           // bb_orient: the stamp the memory must carry to be used (0: none)
    uint  has_ss;
    float bb_coil;              // bb_profile: the scales of the shapes
    float bb_helix;
    float bb_sheet;
    uint  dbg;                  // the bench's probes: sol_write leaves out (1) atom_order, (4) the bonds
    uint  na_count;             // control points of the nucleic band
    uint  der_count;            // derived points
    uint  der_cell;             // the cell of their images: 0 none, 1 fpar0's, 2 fpar1's
    uint  _pad2;
};
MD_SHADER_ARGS(Args);

// params
static const uint P_MIN   = 0u;    // 3 ordered words
static const uint P_MAX   = 3u;    // 3 ordered words
static const uint P_LO    = 6u;    // 3 floats
static const uint P_INV   = 9u;    // 1 / edge
static const uint P_DIM   = 10u;   // 3 uints
static const uint P_EDGE  = 13u;

// point_type: the type in the low 16 bits (md_render's radius table index),
// the point's flags in the high 16 (none yet)
static const uint POINT_TYPE_MASK = 0xFFFFu;

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

// Atom i's position: from the uploaded frame (a frame's regroup), or the caller's
float3 atom_pos(Args a, uint i) {
    if (a.quantized != 0u) return dequant(a.qsrc, a.qpar, upload_index(a, i));
    uint* p = a.position + i * (a.stride >> 2);
    return float3(asfloat(p[0]), asfloat(p[1]), asfloat(p[2]));
}

// The small molecules' records: those of the last grouping (by the parity of
// the regroups so far: sol_offsets counts this one after the write), and
// those this regroup writes
uint4* recs_cur(Args a)  { return (uint4*)(((a.regroups[0] & 1u) != 0u) ? a.recs1 : a.recs0); }
uint4* recs_next(Args a) { return (uint4*)(((a.regroups[0] & 1u) != 0u) ? a.recs0 : a.recs1); }
static const uint REC_MIXED = 0x80000000u;

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
        uint j = i0 + gi;
        bool valid = j < n;
        uint3 lo = uint3(0xFFFFFFFFu), hi = uint3(0u);
        if (valid) {
            // the one read that is scattered: atoms are in file order, slots in
            // the last grouping's; the count reads it back in order
            float3 p = atom_pos(a, loadAligned<16>(recs_cur(a) + j).x);
            storeAligned<16>((uint4*)a.mol_slot + j, uint4(asuint(p), 0u));
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
groupshared uint3 g_inc[64];
groupshared uint3 g_exc[64];
groupshared uint3 g_base[64];

// By slot, in the last grouping's order: after a short move most of a run of
// slots is still in one cell, and its molecules count with one atomic.
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_count(uint3 tid : SV_DispatchThreadID, uint gi : SV_GroupIndex) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    for (uint i0 = group_start(tid, gi); i0 < n; i0 += grid_threads(n)) {
        uint j = i0 + gi;
        bool valid = j < n;
        uint key = 0xFFFFFFFFu;
        uint3 v = uint3(0u);
        if (valid) {
            uint4 st = loadAligned<16>((uint4*)a.mol_slot + j);
            key = cell_of(a, asfloat(st.xyz));
            uint4 k = loadAligned<16>((uint4*)a.kinds + (loadAligned<16>(recs_cur(a) + j).z & ~REC_MIXED));
            v = uint3(k.x, a.has_bonds != 0u ? k.y : 0u, 1u);
        }
        g_key[gi] = key;
        g_inc[gi] = v;
        GroupMemoryBarrierWithGroupSync();
        for (uint off = 1u; off < 64u; off <<= 1u) {
            uint3 x = gi >= off ? g_inc[gi - off] : uint3(0u);
            GroupMemoryBarrierWithGroupSync();
            g_inc[gi] += x;
            GroupMemoryBarrierWithGroupSync();
        }
        uint3 exc = g_inc[gi] - v;
        g_exc[gi] = exc;
        GroupMemoryBarrierWithGroupSync();
        uint h = gi;                                    // the run's first lane
        while (h > 0u && g_key[h - 1u] == key) h -= 1u;
        if (valid && h == gi) {
            uint e = gi;                                // its last
            while (e + 1u < 64u && g_key[e + 1u] == key) e += 1u;
            uint3 total = g_inc[e] - exc;
            uint3 base;
            base.x = ((Atomic<uint>*)(a.cnt + 4u * key + 0u))->add(total.x);
            base.y = total.y > 0u ? ((Atomic<uint>*)(a.cnt + 4u * key + 1u))->add(total.y) : 0u;
            base.z = ((Atomic<uint>*)(a.cnt + 4u * key + 2u))->add(total.z);
            g_base[gi] = base;
        }
        GroupMemoryBarrierWithGroupSync();
        if (valid) {
            uint3 place = g_base[h] + (exc - g_exc[h]);
            storeAligned<16>((uint4*)a.mol_slot + j, uint4(key, place));
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

// ---- Write, offsets ---------------------------------------------------------------

// By slot, in the last grouping's order: each molecule's atoms and bonds
// at its cell's offset plus its place, and its record at its new slot. In
// the order of the last grouping the places of neighbouring slots are near
// each other (most molecules stay in or next to their cell), so the writes
// land near each other too; in file order they were scattered over the
// whole arrays (GTX 1060, 10M atoms: 29 ms; 4.2 when the order is the last
// grouping's and the molecules moved up to 12 A since).
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_write(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.n_mol;
    uint4* cur = recs_cur(a);
    uint4* nxt = recs_next(a);
    for (uint j = tid.x; j < n; j += grid_threads(n)) {
        uint4 r = loadAligned<16>(cur + j);
        uint4 st = loadAligned<16>((uint4*)a.mol_slot + j);
        uint c = st.x;
        uint4 k = loadAligned<16>((uint4*)a.kinds + (r.z & ~REC_MIXED));
        uint p0 = a.static_atoms + a.cnt[4u * c] + st.y;      // the molecule's first point
        for (uint i = 0u; i < k.x; ++i) {
            uint atom = r.x + a.kind_atom[k.z + i];
            if ((a.dbg & 1u) == 0u) a.atom_order[p0 + i] = atom;
            a.point_type[p0 + i] = a.kind_type[k.z + i];   // flags: none yet
            if (a.has_src != 0u) a.src_point[p0 + i - a.static_atoms] = r.w + i;   // its point in the last grouping
            if ((r.z & REC_MIXED) != 0u) a.mixed_point[atom] = p0 + i;
        }
        if (a.has_bonds != 0u && (a.dbg & 4u) == 0u) {
            uint bt = a.static_bonds + a.cnt[4u * c + 1u] + st.z;
            for (uint b = 0u; b < k.y; ++b) {
                uint l = a.kind_bond_local[k.w + b];
                a.bond_pairs[2u * (bt + b)]      = p0 + (l & 0xFFFFu);
                a.bond_pairs[2u * (bt + b) + 1u] = p0 + (l >> 16u);
                a.bond_ids[bt + b]               = a.ids_old[r.y + b];
            }
        }
        uint bt_new = a.cnt[4u * c + 1u] + st.z;                 // among the small bonds
        storeAligned<16>(nxt + a.cnt[4u * c + 2u] + st.w, uint4(r.x, bt_new, r.z, p0));
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
    if (tid.x == 0u) a.regroups[0] = a.regroups[0] + 1u;  // the records' buffers swap
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
        a.bond_pairs[2u * pos]      = a.mixed_point[a.mixed_atoms[2u * i]];
        a.bond_pairs[2u * pos + 1u] = a.mixed_point[a.mixed_atoms[2u * i + 1u]];
    }
}


// A point's w: its type and flags, or the caller's radius for its atom
uint point_w(Args a, uint p, uint atom) {
    return a.has_radius != 0u ? a.radius_override[atom] : a.point_type[p];
}

// A point into the copy: one aligned 16-byte store, whole, so a wave writes
// whole lines (the copy's block is 16-byte aligned). Not the store `*p`
// makes: Slang marks it Aligned 4 and NVIDIA splits it in four (with the
// radius read back from the copy to complete it, the two-frame gather was
// 6.3 ms on a GTX 1060). Not 12 bytes of 16 either: 0.7 ms slower than the
// whole store on the M4. w: the point's type and flags (or a radius).
[ForceInline]
void store_point(uint* points, uint p, uint3 v, uint r) {
    storeAligned<16>((uint4*)points + p, uint4(v, r));
}

// The copy that is drawn: point p is atom atom_order[p]. Writes in order,
// reads where the atoms are (in runs: the atoms of a residue or a molecule
// are neighbours in both orders).
[shader("compute")]
[numthreads(64, 1, 1)]
void sol_gather(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    for (uint p = tid.x; p < n; p += grid_threads(n)) {
        uint atom = a.atom_order[p];
        uint* q = a.position + atom * (a.stride >> 2u);
        store_point(a.points, p, uint3(q[0], q[1], q[2]), point_w(a, p, atom));
        if (a.has_velocity != 0u) {
            uint* w = a.velocity + atom * (a.vel_stride >> 2u);
            storeAligned<8>((uint2*)a.vel_out + p, uint2(w[0], w[1]));
        }
    }
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


// The copy that is drawn, from resident frames in point order: every read in
// order, every write too (w: point_type; the caller's radius instead is
// fetched by atom). A
// kernel per number of frames (COUNT a constant): each frame's parameters
// and the cell in registers for only the frames it combines.
[ForceInline]
void fr_gather_body(uint3 tid, const uint COUNT) {
    Args a = MD_ARGS;
    uint n = a.atom_count;
    uint step = a.fr_threads;
    const bool vel_out = COUNT > 1u && a.has_fvel != 0u;
    uint* points = a.points;
    const FrameP f0 = frame_p(a.fpar0);
    FrameP f1 = f0, f2 = f0, f3 = f0;
    CellP  cell;
    if (COUNT >= 2u) { f1 = frame_p(a.fpar1); cell = cell_p(COUNT == 4u ? a.fpar1 : a.fpar0); }
    if (COUNT == 4u) { f2 = frame_p(a.fpar2); f3 = frame_p(a.fpar3); }
    for (uint p = tid.x; p < n; p += step) {
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
        store_point(points, p, asuint(pos), a.has_radius != 0u ? a.radius_override[a.atom_order[p]] : a.point_type[p]);
        if (vel_out) storeAligned<8>((uint2*)a.vel_out + p, uint2(f32tof16(vel.x) | (f32tof16(vel.y) << 16u), f32tof16(vel.z)));
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


// ---- The protein backbone (the cartoon) ------------------------------------------
//
// As md_gl's cartoon (compute_spline_orient.vert, backbone.vert):
//
//     bb_support  per control point, the ribbon's axis: the binormal of the
//                 CA trace there, B = N x T with T = CA[i+1] - CA[i-1] and N =
//                 CA[i-1] - 2 CA[i] + CA[i+1] (the B-spline's first and second
//                 derivatives at the point; the ends borrow their neighbour's
//                 curvature). The CA virtual bond angle stays well below 180
//                 degrees, so N does not vanish; where B is short (|B| < 0.75
//                 A, a straight stretch) the carbonyl (O - C) across T takes
//                 over, blended in from 0.75 down to 0.25.
//     bb_orient   per piece (a thread walks it): whether each axis is turned
//                 against the one before, with hysteresis: the last update's
//                 decision stands until the angle under it passes 90 degrees
//                 by 40 (md_gl's u_flip_threshold). Then the axes written with
//                 their signs, which md_render takes as given.
//     bb_profile  per control point, (a0, a1, r0, r1) from the secondary
//                 structure (see md_render_mol_cartoon_profile).

static const uint  BB_FIRST = 0x80000000u;
static const uint  BB_LAST  = 0x40000000u;
static const uint  BB_POINT = 0x3FFFFFFFu;
static const float BB_FLIP_THRESHOLD = -0.64278760968;   // -sin(40 degrees)

float3 bb_pos(Args a, uint cp) {
    uint p = a.bb_index[cp] & BB_POINT;
    return ld3(a.points, 4u * p);
}

float3 bb_any_perpendicular(float3 t) {
    float3 v = abs(t.x) < 0.9 ? float3(1.0, 0.0, 0.0) : float3(0.0, 1.0, 0.0);
    return normalize(cross(t, v));
}

[shader("compute")]
[numthreads(64, 1, 1)]
void bb_support(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.bb_count;
    for (uint j = tid.x; j < n; j += grid_threads(n)) {
        uint ie = a.bb_index[j];
        bool hp = (ie & BB_FIRST) == 0u, hn = (ie & BB_LAST) == 0u;
        uint jp = j, jn = j;
        if (hp) jp = j - 1u;
        if (hn) jn = j + 1u;
        // the curvature's three points: the ends borrow their neighbour's
        uint i0 = jp, i1 = j, i2 = jn;
        if (!hp) {
            i0 = j; i1 = jn; i2 = jn;
            if ((a.bb_index[jn] & BB_LAST) == 0u) i2 = jn + 1u;
        } else if (!hn) {
            i0 = jp; i1 = jp; i2 = j;
            if ((a.bb_index[jp] & BB_FIRST) == 0u) i0 = jp - 1u;
        }
        float3 T  = bb_pos(a, jn) - bb_pos(a, jp);
        float3 N  = bb_pos(a, i0) - 2.0 * bb_pos(a, i1) + bb_pos(a, i2);
        float  tt = dot(T, T);
        float3 tg = tt > 0.0 ? T * rsqrt(tt) : float3(0.0, 0.0, 1.0);
        float3 B  = cross(N, tg);
        float  lb = length(B);
        float3 s;
        if (lb >= 0.75) {
            s = B / lb;
        } else {
            uint2  co = uint2(a.bb_carbonyl[2u * j], a.bb_carbonyl[2u * j + 1u]);
            float3 oc = ld3(a.points, 4u * co.y) - ld3(a.points, 4u * co.x);
            oc -= tg * dot(oc, tg);
            float  o2 = dot(oc, oc);
            float3 fb = o2 > 1.0e-6 ? oc * rsqrt(o2) : bb_any_perpendicular(tg);
            if (lb <= 0.25) {
                s = fb;
            } else {
                float3 sb = B / lb;
                if (dot(fb, sb) < 0.0) fb = -fb;
                float3 m = lerp(fb, sb, smoothstep(0.25, 0.75, lb));
                s = normalize(m - tg * dot(m, tg));
            }
        }
        uint* o = a.bb_axis + 4u * j;
        o[0] = asuint(s.x); o[1] = asuint(s.y); o[2] = asuint(s.z); o[3] = 0u;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void bb_orient(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.bb_chains;
    // the memory is the last update's when its stamp is the one asked for
    bool history = a.bb_history != 0u && a.bb_flip[a.bb_count] == a.bb_history;
    for (uint c = tid.x; c < n; c += grid_threads(n)) {
        uint   b = a.bb_chain[c], e = a.bb_chain[c + 1u];
        float  sg = 1.0;
        float3 s  = ld3(a.bb_axis, 4u * b);
        for (uint j = b; j < e; ++j) {
            float3 sn = float3(0.0);
            bool   last = j + 1u >= e;
            if (!last) sn = ld3(a.bb_axis, 4u * (j + 1u));
            uint* o = a.bb_axis + 4u * j;
            o[0] = asuint(s.x * sg); o[1] = asuint(s.y * sg); o[2] = asuint(s.z * sg);
            if (last) break;
            float d = dot(s, sn);
            bool flip = d < 0.0;
            if (history) {
                bool  prev  = (a.bb_flip[j] & 1u) != 0u;
                float sigma = prev ? -1.0 : 1.0;
                flip = (sigma * d < BB_FLIP_THRESHOLD) ? !prev : prev;
            }
            a.bb_flip[j] = flip ? 1u : 0u;
            if (flip) sg = -sg;
            s = sn;
        }
    }
}

// A control point's structure as fractions (helix, sheet), as md_gl's
// md_gl_secondary_structure_convert: 3-10, alpha and pi helices are helix,
// sheet and bridge are sheet, the rest coil
float2 bb_fractions(Args a, uint cp) {
    if (a.has_ss == 0u) return float2(0.0);
    int ss = int(a.bb_ss[a.bb_segment[cp]]);
    if (ss >= 4 && ss <= 6) return float2(1.0, 0.0);
    if (ss == 7 || ss == 8) return float2(0.0, 1.0);
    return float2(0.0);
}

// The fractions where the parts of control points l and l + 1 meet: the
// Catmull-Rom spline of the fractions (md_gl's blend) halfway, over l - 1 ..
// l + 2 clamped to the piece. Both parts compute it alike: they meet
// without a step (and so without caps).
float2 bb_joint(Args a, uint l) {
    uint r  = l + 1u;
    uint lp = l, rn = r;
    if ((a.bb_index[l] & BB_FIRST) == 0u) lp = l - 1u;
    if ((a.bb_index[r] & BB_LAST)  == 0u) rn = r + 1u;
    float2 f = (-bb_fractions(a, lp) + 9.0 * bb_fractions(a, l) + 9.0 * bb_fractions(a, r) - bb_fractions(a, rn)) * (1.0 / 16.0);
    return saturate(f);
}

// (half width, radius) of the profile at fractions f: md_gl's stadiums,
// coil (0, 0.2), helix (1.0, 0.2), sheet (1.4, 0.1), each times its scale
float2 bb_shape(Args a, float2 f) {
    float coil = max(1.0 - f.x - f.y, 0.0);
    return coil * float2(0.0, 0.2) * a.bb_coil + f.x * float2(1.0, 0.2) * a.bb_helix + f.y * float2(1.4, 0.1) * a.bb_sheet;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void bb_profile(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.bb_count;
    for (uint j = tid.x; j < n; j += grid_threads(n)) {
        uint   ie = a.bb_index[j];
        float2 f0 = bb_fractions(a, j), f1 = f0;      // the ends: the control point's own
        if ((ie & BB_FIRST) == 0u) f0 = bb_joint(a, j - 1u);
        if ((ie & BB_LAST)  == 0u) f1 = bb_joint(a, j);
        float2 s0 = bb_shape(a, f0), s1 = bb_shape(a, f1);
        uint* o = a.bb_out + 4u * j;
        o[0] = asuint(s0.x); o[1] = asuint(s1.x); o[2] = asuint(s0.y); o[3] = asuint(s1.y);
    }
}

// The nucleic band (viamd's cartoon_glyphs, after Lindow et al.): per
// control point (a C4'), the long axis of its profile across the curve and
// the way to the base, cross(T, F), F from the C4' to the point it faces (the
// glycosidic N, else the C1'). Its flat faces then face the base. F is on the
// same side all along a strand: the signs need no deciding, no memory. T: the
// chord over its neighbours (one-sided at an end), as the curve's direction
// through it (md_render makes the axis orthogonal to the curve itself).
[shader("compute")]
[numthreads(64, 1, 1)]
void na_axis(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.na_count;
    for (uint j = tid.x; j < n; j += grid_threads(n)) {
        uint ie = a.na_index[j];
        uint jp = j, jn = j;
        if ((ie & BB_FIRST) == 0u) jp = j - 1u;
        if ((ie & BB_LAST)  == 0u) jn = j + 1u;
        float3 c  = ld3(a.points, 4u * (ie & BB_POINT));
        float3 T  = ld3(a.points, 4u * (a.na_index[jn] & BB_POINT)) - ld3(a.points, 4u * (a.na_index[jp] & BB_POINT));
        float  tt = dot(T, T);
        float3 tg = tt > 0.0 ? T * rsqrt(tt) : float3(0.0, 0.0, 1.0);
        uint   f  = a.na_face[j];
        float3 F  = ld3(a.points, 4u * (f != 0xFFFFFFFFu ? f : (ie & BB_POINT))) - c;   // no face: none
        float3 s  = cross(tg, F - tg * dot(F, tg));
        float  ls = length(s);
        s = ls > 1.0e-4 ? s / ls : bb_any_perpendicular(tg);
        uint* o = a.na_axis + 4u * j;
        o[0] = asuint(s.x); o[1] = asuint(s.y); o[2] = asuint(s.z); o[3] = 0u;
    }
}

// The derived points (md_render_mol_topo_t): each the weighted sum of its
// members' points, found where they are now (mixed_point: the atoms of a
// small molecule move among the points with every regroup), each at its
// image nearest the first's where the update has a cell (the reference
// frame's, der_cell). Written after the atoms' points, with the radius it
// stands for as its w; its velocity the same sum of theirs.
[shader("compute")]
[numthreads(64, 1, 1)]
void derive(uint3 tid : SV_DispatchThreadID) {
    Args a = MD_ARGS;
    uint n = a.der_count;
    bool vel = a.has_velocity != 0u || a.has_fvel != 0u;
    CellP cell;
    cell.ax = cell.bx = cell.by = cell.cx = cell.cy = cell.cz = 0.0;
    cell.r  = float3(0.0);                          // no cell: mi leaves d as it is
    if (a.der_cell != 0u) cell = cell_p(a.der_cell == 1u ? a.fpar0 : a.fpar1);
    for (uint d = tid.x; d < n; d += grid_threads(n)) {
        uint   o0  = a.der_offset[d], o1 = a.der_offset[d + 1u];
        float3 ref = ld3(a.points, 4u * a.mixed_point[a.der_member[o0]]);
        float3 x = float3(0.0), v = float3(0.0);
        for (uint o = o0; o < o1; ++o) {
            uint  p = a.mixed_point[a.der_member[o]];
            float w = asfloat(a.der_weight[o]);
            x += w * (ref + mi(cell, ld3(a.points, 4u * p) - ref));
            if (vel) {
                uint2 h = loadAligned<8>((uint2*)a.vel_out + p);
                v += w * float3(f16tof32(h.x & 0xFFFFu), f16tof32(h.x >> 16u), f16tof32(h.y & 0xFFFFu));
            }
        }
        uint q = a.atom_count + d;
        store_point(a.points, q, asuint(x), a.der_radius[d]);
        if (vel) storeAligned<8>((uint2*)a.vel_out + q, uint2(f32tof16(v.x) | (f32tof16(v.y) << 16u), f32tof16(v.z)));
    }
}

// gto_dmat.slang — EXPERIMENT (unittest/test_gto_dmat.c): a density matrix in 16 x 16 tiles of the AO
// pairs that meet, built from the occupied orbitals, D_tile(I, J) = P_I^T diag(w) P_J.
//
// P holds the orbitals' coefficients by tile: P[(tile * num_mo + k) * 16 + m] is orbital k's
// coefficient of the tile's AO m (in the tiled AO order). A work item is a row tile I and up to 8 of
// its column tiles (consecutive entries of the row's list); a group of 256 threads does one, 32
// threads per tile pair, each thread a 2 x 4 block of the tile. The orbitals go through groupshared
// memory KC at a time.

#include "../../src/shaders/md_gpu.slang"

struct DmatArgs {
    uint   num_mo;
    uint   num_items;
    uint   items_x;      // groups along x: item = gid.y * items_x + gid.x
    uint   _pad0;
    float* P;
    float* w;
    uint4* items;        // x: row tile, y: first entry, z: entries (1 .. 8)
    uint*  cols;         // column tile of each entry
    float* tiles;        // entry e: tiles[e * 256 + 16 m + c]
};

struct DmatRoot { DmatArgs* args; };
[[vk::push_constant]] ConstantBuffer<DmatRoot> dmat_root;

#define KC 32
groupshared float4 sa[KC * 4];          // [k][16]: the row tile's chunk, weighted
groupshared float4 sb[8 * KC * 4];      // [slot][k][16]: the column tiles' chunks

[shader("compute")]
[numthreads(256, 1, 1)]
void dmat_tiles(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    DmatArgs a = *dmat_root.args;
    uint item = gid.y * a.items_x + gid.x;
    if (item >= a.num_items) return;                       // uniform per group
    uint4 it = a.items[item];
    uint I = it.x, e0 = it.y, cnt = it.z;
    uint r = a.num_mo;
    uint slot = tid >> 5;                                  // the column tile this thread works on
    uint l = tid & 31u;
    uint m0 = (l >> 2) * 2u;                               // rows m0, m0 + 1
    uint cq = l & 3u;                                      // columns 4 cq .. 4 cq + 3
    float4 acc0 = float4(0.0), acc1 = float4(0.0);

    for (uint k0 = 0u; k0 < r; k0 += uint(KC)) {
        uint kn = min(uint(KC), r - k0);
        // the row tile's chunk (KC x 16 = 128 float4, weighted), then the column tiles' (8 x 128)
        if (tid < uint(KC) * 4u) {
            uint k = tid >> 2;
            float4 v = float4(0.0);
            if (k < kn) {
                float* src = a.P + ((I * r + k0 + k) * 16u + (tid & 3u) * 4u);
                v = float4(src[0], src[1], src[2], src[3]) * a.w[k0 + k];
            }
            sa[tid] = v;
        }
        for (uint i = tid; i < 8u * uint(KC) * 4u; i += 256u) {
            uint s = i / (uint(KC) * 4u), rem = i % (uint(KC) * 4u), k = rem >> 2;
            float4 v = float4(0.0);
            if (s < cnt && k < kn) {
                uint J = a.cols[e0 + s];
                float* src = a.P + ((J * r + k0 + k) * 16u + (rem & 3u) * 4u);
                v = float4(src[0], src[1], src[2], src[3]);
            }
            sb[i] = v;
        }
        GroupMemoryBarrierWithGroupSync();
        if (slot < cnt) {                                  // uniform per wave of 32
            for (uint k = 0u; k < kn; ++k) {
                float4 ra = sa[k * 4u + (m0 >> 2)];         // rows m0 .. m0 + 3 of the chunk's row k
                float x0 = ra[m0 & 3u], x1 = ra[(m0 & 3u) + 1u];
                float4 b = sb[(slot * uint(KC) + k) * 4u + cq];
                acc0 += x0 * b;
                acc1 += x1 * b;
            }
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (slot < cnt) {
        float* dst = a.tiles + (e0 + slot) * 256u + m0 * 16u + cq * 4u;
        dst[0] = acc0.x; dst[1] = acc0.y; dst[2] = acc0.z; dst[3] = acc0.w;
        dst[16] = acc1.x; dst[17] = acc1.y; dst[18] = acc1.z; dst[19] = acc1.w;
    }
}

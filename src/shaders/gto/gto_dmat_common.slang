// gto_dmat_common.slang
// Reading a density matrix in tiles (md_gto_gpu_dmat, md_gto.c). Included by gto_common.slang and
// eval_gto_density.slang.

#ifndef GTO_DMAT_COMMON_SLANG
#define GTO_DMAT_COMMON_SLANG

// Element (i, j) of a matrix in tiles (md_gto_gpu_dmat): AO a sits at tiled index tao[a]; the tile pair
// (I, J), I <= J, of tiles of 16 is entry map[I * nt + J] of the pool (16 x 16 floats, row a, column c at
// 16 a + c), or 0xFFFFFFFF where it is not kept (zero: no point has both AOs above the cutoff).
#define GTO_DMAT_NONE 0xFFFFFFFFu
float gto_dmat_at(uint* tao, uint* map, float* tiles, uint nt, uint i, uint j) {
    uint ti = tao[i], tj = tao[j];
    uint lo = min(ti, tj), hi = max(ti, tj);
    uint e = map[(lo >> 4) * nt + (hi >> 4)];
    return e == GTO_DMAT_NONE ? 0.0f : tiles[e * 256u + (lo & 15u) * 16u + (hi & 15u)];
}

#endif

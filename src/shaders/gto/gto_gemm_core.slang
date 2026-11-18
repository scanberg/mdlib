// gto_gemm_core.slang
//
// THE GEMM CORE of the two-pass grid evaluations (count + Phi passes, eval_gto_density_gemm.slang): one main loop,
// specialised at compile time by what it multiplies and what is made of the product.
//
//   X[r][m] = sum_k A[r][k] B[k][m]     r: the GP rows of a group, m < ncols the columns, k < nk (the block's AOs)
//
// A TEMPLATE over the tile and the threads (Slang generics, as CUTLASS's templates): GP rows per group, GM columns per
// tile, TX x TY threads, and the groupshared memory of the chunks (a type: GC_MEM, sized for the kernel's tile, so that
// a kernel with narrow tiles keeps its occupancy, and several tile shapes live in one source). The product goes
// through the columns in tiles of GM, and for each tile through k in chunks of GK held in groupshared memory. Thread
// (tx, ty) holds the RP x RM part of X (RP = GP / TY, RM = GM / TX) of rows ty * RP + i and columns m0 + tx * RM + j.
//
// Two ways to fill the chunks, one multiply-add loop (gc_mac):
//   gemm_core         OPERANDS read element by element (IGemmOperands), the next chunk loaded into registers while the
//                     current one is multiplied: the grid evaluations, whose Phi is in memory
//   gemm_core_staged  a STAGER (IGemmStager) writes each chunk into groupshared memory itself, as many rows of k as it
//                     has (at most GK): the evaluations that make the AO values of a chunk in the group (the isosurface
//                     corners and vertices, gto_iso_gpu.slang)
//   gemm_core_tri     staged, for a quadratic form streamed over a matrix TRIANGULAR BY CHUNKS (B = D' with 2 D before a
//                     column's chunk, D on its diagonal block, 0 after; M = GK, so a thread's columns are one chunk): a
//                     thread multiplies the chunks up to its own and, at its own, the epilogue's diag() takes A . X while
//                     the chunk's A is still in memory (the AO values made once per column tile, with no Phi in memory:
//                     the isosurface D form, gto_iso_gpu.slang)
//
// The chunks are plain floats in groupshared memory, written and read a float at a time (the driver merges a thread's
// consecutive reads). r3 read them as float4 (a float4 array, a component written at a time): on an RTX 3080 the
// orbitals' GEMM lost 13%, ELF's 5%, the density's about even; r4 is r2's plain floats again. After each column tile the EPILOGUE
// takes that part of X; it may read anything (A at the same rows, say) and keeps what it accumulates in its own state.
//
//   OPERANDS   IGemmOperands: A[r][k] and B[k][m] as the loads read them (consecutive threads read consecutive rows of
//              A and consecutive columns of B: both streams coalesced where the operands are laid out so). The core
//              keeps the loads in range (k < nk, m < ncols) and, for a symmetric product, to k >= m only.
//   STAGER     IGemmStager: begin(m0) before a column tile; stage() after a barrier writes the next chunk (A [kk][r] at
//              [0, GK GP), B [kk][mm] at GK GP of the memory) and returns its rows of k, GC_END after the last.
//   EPILOGUE   IGemmEpilogue: the part of X of each column tile (take<R, M>, R = RP and M = RM).
//
// What uses it: the density (eval_gto_density_core.slang: symmetric, ncols = nk), ELF and LOL
// (eval_gto_localization_core.slang: symmetric, four rows a voxel), the orbitals (eval_gto_mo_gemm.slang: B the
// orbitals' coefficients, ncols the orbitals), a density matrix's tiles from orbitals (gto_dmat.slang: P_I^T diag(w) P_J,
// 16 x 128 per group of 16 x 4 threads), the isosurface GEMMs (gto_iso_gpu.slang, staged).
//
// A SYMMETRIC product (B = w o D of a symmetric D, w 1 on the diagonal and 2 off it, k >= m only) sums
// sum_{m <= k} w D[k][m] A[r][k] - half the multiply-adds - and is exact for a quadratic form a^T D a = sum_m a[m] X[m]
// with X of the same a; for a symmetric bilinear form 2 a^T D b = sum_m (a[m] Xb[m] + b[m] Xa[m]) with the X of both.
//
// What it costs: GP x GM x n / 2 (symmetric) multiply-adds per column tile, of which the epilogue's are a few per
// element of X. The specialisations (Slang generics, instantiated per kernel) are each a kernel of their own, with the
// operands' and the epilogue's code inlined into the loop: one source, as many binaries as before.
//
// A kernel declares its memory, GC_MEM(Name, GP, GM), and calls gemm_core<GP, GM, TX, TY, Name, ...>(...) with its
// thread index and the thread's (tx, ty) (tid = ty TX + tx is the core's loads' order; the layout of the group is the
// kernel's). Past the core, Name::ld / Name::st are the kernel's (gemm_row_sums, an epilogue's staging).
//
// r5: the template (r4: the same loop with the tile and the threads as each file's macros, one shape per source).
// r6: the chunks' memory a kernel's choice (GC_MEM, GC_MEM4). r7: GC_MEM4's float4 reads indexed in float4s.

#ifndef GTO_GEMM_CORE_SLANG
#define GTO_GEMM_CORE_SLANG

// For measurements: bumped with every change to the core's loop
#define GTO_GEMM_CORE_REVISION 7

#define GK     16               // k per chunk
#define GC_END 0xFFFFFFFFu      // IGemmStager.stage: no more chunks

// The groupshared memory of a kernel's chunks: [0, GK GP) the A chunk [kk][r], [GK GP, GK (GP + GM)) the B chunk
// [kk][mm]; after the core, the kernel's (the epilogues' reductions). Two kinds, a kernel's choice (measured per kernel,
// RTX 3080):
//   GC_MEM   floats, read a float at a time: the grid kernels (their GEMM pass lost 5-13% on float4 reads, gemm-core r3)
//   GC_MEM4  a thread's rows (columns) read as float4 where R (M) is a multiple of 4: SPIR-V float4 elements (one
//            16-byte read; a float written at a time), Metal floats (Metal may write one
//            component of a threadgroup vector as a read-modify-write of the whole vector, and two threads writing
//            components of one float4 would lose a write; gto_iso_gpu.slang's K_BUF): the isosurface kernels (their
//            tiled GEMM 1.3-1.5x slower on GC_MEM's float reads at 20-30 functions, gemm-core r5)
interface IGcMem {
    static const int kVec;              // 1: ld4 a 16-byte read (GC_MEM4), the multiply-adds read float4
    static float  ld(uint i);
    static float4 ld4(uint i4);
    static void   st(uint i, float v);
};
#define GC_MEM(NAME, GPV, GMV)                                                                      \
    groupshared float NAME##_smem[GK * ((GPV) + (GMV))];                                            \
    struct NAME : IGcMem {                                                                          \
        static const int kVec = 0;                                                                  \
        static float  ld(uint i)          { return NAME##_smem[i]; }                                \
        static float4 ld4(uint i4)        { return float4(NAME##_smem[4u * i4], NAME##_smem[4u * i4 + 1u], \
                                                          NAME##_smem[4u * i4 + 2u], NAME##_smem[4u * i4 + 3u]); } \
        static void   st(uint i, float v) { NAME##_smem[i] = v; }                                   \
    };
#define GC_MEM4(NAME, GPV, GMV)                                                                     \
    groupshared float4 NAME##_d4[GK * ((GPV) + (GMV)) / 4];                                         \
    groupshared float  NAME##_d1[GK * ((GPV) + (GMV))];                                             \
    struct NAME : IGcMem {                                                                          \
        static const int kVec = 1;                                                                  \
        static float ld(uint i) {                                                                   \
            __target_switch { case spirv: return NAME##_d4[i >> 2][i & 3u];                        \
                              default:    return NAME##_d1[i]; }                                    \
        }                                                                                           \
        static float4 ld4(uint i4) {                                                                \
            __target_switch { case spirv: return NAME##_d4[i4];                                     \
                              default:    return float4(NAME##_d1[4u * i4], NAME##_d1[4u * i4 + 1u], \
                                                        NAME##_d1[4u * i4 + 2u], NAME##_d1[4u * i4 + 3u]); } \
        }                                                                                           \
        static void st(uint i, float v) {                                                           \
            __target_switch { case spirv: NAME##_d4[i >> 2][i & 3u] = v;                           \
                              default:    NAME##_d1[i] = v; }                                       \
        }                                                                                           \
    };

interface IGemmOperands {
    // A[row][k], row < GP, k < nk
    float a(uint row, uint k);
    // B[k][m], k < nk, m < ncols (symmetric: k >= m, weighted)
    float b(uint k, uint m);
};

interface IGemmStager {
    // Before the column tile at m0
    [mutating] void begin(uint m0, uint tid);
    // The next chunk into the memory (every thread, after a barrier; barriers of its own allowed): its rows of k (1 .. GK),
    // or GC_END
    [mutating] uint stage(uint m0, uint tid);
};

interface IGemmEpilogue {
    // The part of X of the columns [m0, m0 + GM): X[i][j] at row ty * R + i, column m0 + tx * M + j (only those below
    // ncols are of the product)
    [mutating] void take<let R : int, let M : int>(uint m0, uint ncols, uint tx, uint ty, float X[R][M]);
};

interface IGemmTriEpilogue : IGemmEpilogue {
    // gemm_core_tri: the thread's columns [k0, k0 + GK) complete at the chunk of k0; that chunk's A is in the memory
    // (row kk at kk GP + r). Called by the threads whose chunk it is (not group-uniform: no barriers)
    [mutating] void diag<let R : int, let M : int>(uint k0, uint tx, uint ty, float X[R][M]);
};

// A fence for groupshared memory without an execution barrier (for code that not every thread runs): on Vulkan
// OpMemoryBarrier; on Metal simdgroup_barrier, as threadgroup_barrier would also wait for the other SIMD-groups
void gc_fence() {
    __target_switch {
    case metal: __intrinsic_asm "simdgroup_barrier(mem_flags::mem_threadgroup)";
    default: GroupMemoryBarrier();
    }
}

// Row kk of the chunk in Mem into the thread's X. The float4 reads' indices in float4s (R, M, CP, CM multiples of 4): a
// float index divided by 4 does not fold into the reads' constant offsets (unsigned wrap), r6 an add and a shift per read
[ForceInline]
void gc_mac<let CP : int, let CM : int, let R : int, let M : int, Mem : IGcMem>(inout float X[R][M], uint kk, uint tx, uint ty) {
    float a[R];
    float bb[M];
    if (Mem::kVec != 0 && R % 4 == 0) {
        [ForceUnroll] for (int q = 0; q < R / 4; ++q) {
            const float4 v = Mem::ld4(kk * uint(CP / 4) + ty * uint(R / 4) + uint(q));
            a[4 * q] = v.x; a[4 * q + 1] = v.y; a[4 * q + 2] = v.z; a[4 * q + 3] = v.w;
        }
    } else {
        [ForceUnroll] for (int i = 0; i < R; ++i) a[i] = Mem::ld(kk * uint(CP) + ty * uint(R) + uint(i));
    }
    if (Mem::kVec != 0 && M % 4 == 0) {
        [ForceUnroll] for (int q = 0; q < M / 4; ++q) {
            const float4 v = Mem::ld4(uint(GK * CP / 4) + kk * uint(CM / 4) + tx * uint(M / 4) + uint(q));
            bb[4 * q] = v.x; bb[4 * q + 1] = v.y; bb[4 * q + 2] = v.z; bb[4 * q + 3] = v.w;
        }
    } else {
        [ForceUnroll] for (int j = 0; j < M; ++j) bb[j] = Mem::ld(uint(GK * CP) + kk * uint(CM) + tx * uint(M) + uint(j));
    }
    [ForceUnroll] for (int i = 0; i < R; ++i)
        [ForceUnroll] for (int j = 0; j < M; ++j) X[i][j] += a[i] * bb[j];
}

[ForceInline]
void gemm_core<let CP : int, let CM : int, let CTX : int, let CTY : int, Mem : IGcMem, O : IGemmOperands, E : IGemmEpilogue>(
        O ops, inout E epi, uint nk, uint ncols, bool symmetric, uint tid, uint tx, uint ty) {
    static const int NT   = CTX * CTY;
    static const int R    = CP / CTY;
    static const int M    = CM / CTX;
    static const int PF_A = (GK * CP) / NT;
    static const int PF_B = (GK * CM) / NT;
    for (uint m0 = 0u; m0 < ncols; m0 += uint(CM)) {
        float X[R][M];
        [ForceUnroll] for (int i = 0; i < R; ++i)
            [ForceUnroll] for (int j = 0; j < M; ++j) X[i][j] = 0.0f;

        // the next chunk, in registers
        float pf_a[PF_A];
        float pf_b[PF_B];

        // Loads the chunk of k from k0 into pf_a / pf_b: consecutive threads, consecutive rows (A) and columns (B)
        #define GC_LOAD_CHUNK(k0_)                                                                  \
        {                                                                                           \
            [ForceUnroll] for (int r = 0; r < PF_A; ++r) {                                          \
                uint idx = tid + uint(r) * uint(NT);                                                \
                uint k   = (k0_) + idx / uint(CP);                                                  \
                float v  = ops.a(idx % uint(CP), min(k, nk - 1u));                                  \
                pf_a[r]  = (k < nk) ? v : 0.0f;                                                     \
            }                                                                                       \
            [ForceUnroll] for (int r = 0; r < PF_B; ++r) {                                          \
                uint idx = tid + uint(r) * uint(NT);                                                \
                uint k   = (k0_) + idx / uint(CM);                                                  \
                uint m   = m0 + idx % uint(CM);                                                     \
                float v  = 0.0f;                                                                    \
                if (k < nk && m < ncols && (!symmetric || k >= m)) v = ops.b(k, m);                 \
                pf_b[r] = v;                                                                        \
            }                                                                                       \
        }

        // a symmetric product starts at the tile's own rows of B: the ones above are zero
        const uint k_first = symmetric ? m0 : 0u;
        if (k_first < nk) {                      // nk 0 (an empty block of the orbitals): nothing to load
            GC_LOAD_CHUNK(k_first);
        }

        for (uint k0 = k_first; k0 < nk; k0 += uint(GK)) {
            GroupMemoryBarrierWithGroupSync();   // the previous chunk consumed
            [ForceUnroll] for (int r = 0; r < PF_A; ++r) Mem::st(tid + uint(r) * uint(NT), pf_a[r]);
            [ForceUnroll] for (int r = 0; r < PF_B; ++r) Mem::st(uint(GK * CP) + tid + uint(r) * uint(NT), pf_b[r]);
            GroupMemoryBarrierWithGroupSync();

            if (k0 + uint(GK) < nk) {
                GC_LOAD_CHUNK(k0 + uint(GK));    // in flight while this chunk is multiplied
            }

            [ForceUnroll]
            for (int kk = 0; kk < GK; ++kk) gc_mac<CP, CM, R, M, Mem>(X, uint(kk), tx, ty);
        }
        #undef GC_LOAD_CHUNK

        epi.take<R, M>(m0, ncols, tx, ty, X);
    }
}

[ForceInline]
void gemm_core_staged<let CP : int, let CM : int, let CTX : int, let CTY : int, Mem : IGcMem, S : IGemmStager, E : IGemmEpilogue>(
        inout S st, inout E epi, uint ncols, uint tid, uint tx, uint ty) {
    static const int R = CP / CTY;
    static const int M = CM / CTX;
    for (uint m0 = 0u; m0 < ncols; m0 += uint(CM)) {
        float X[R][M];
        [ForceUnroll] for (int i = 0; i < R; ++i)
            [ForceUnroll] for (int j = 0; j < M; ++j) X[i][j] = 0.0f;
        st.begin(m0, tid);
        for (;;) {
            GroupMemoryBarrierWithGroupSync();   // the previous chunk consumed (and the epilogue's staging)
            const uint kn = st.stage(m0, tid);   // uniform
            if (kn == GC_END) break;
            GroupMemoryBarrierWithGroupSync();   // the chunk staged
            for (uint kk = 0u; kk < kn; ++kk) gc_mac<CP, CM, R, M, Mem>(X, kk, tx, ty);
        }
        epi.take<R, M>(m0, ncols, tx, ty, X);
    }
}

[ForceInline]
void gemm_core_tri<let CP : int, let CM : int, let CTX : int, let CTY : int, Mem : IGcMem, S : IGemmStager, E : IGemmTriEpilogue>(
        inout S st, inout E epi, uint ncols, uint tid, uint tx, uint ty) {
    static const int R = CP / CTY;
    static const int M = CM / CTX;
    // (M == GK: a thread's columns are one chunk of k; every chunk GK rows, zero past the end)
    for (uint m0 = 0u; m0 < ncols; m0 += uint(CM)) {
        float X[R][M];
        [ForceUnroll] for (int i = 0; i < R; ++i)
            [ForceUnroll] for (int j = 0; j < M; ++j) X[i][j] = 0.0f;
        const uint c0 = m0 + tx * uint(M);                 // the thread's first column
        st.begin(m0, tid);
        for (uint k0 = 0u;; k0 += uint(GK)) {
            GroupMemoryBarrierWithGroupSync();             // the previous chunk consumed
            const uint kn = st.stage(m0, tid);             // uniform
            if (kn == GC_END) break;
            GroupMemoryBarrierWithGroupSync();             // the chunk staged
            if (c0 < ncols && k0 <= c0) {                  // the chunks after the thread's own are zero in B
                [ForceUnroll]
                for (int kk = 0; kk < GK; ++kk) gc_mac<CP, CM, R, M, Mem>(X, uint(kk), tx, ty);
                if (k0 == c0) {
                    gc_fence();
                    epi.diag<R, M>(k0, tx, ty, X);
                }
            }
        }
        epi.take<R, M>(m0, ncols, tx, ty, X);
    }
}

// The sum over the CTX threads of each row of one value per row (v[i]: row ty * R + i), into out_row (row r's sum,
// read by thread r < CP after this returns; every thread calls it). Through Mem, so after the core.
[ForceInline]
void gemm_row_sums<let CP : int, let CTX : int, let R : int, Mem : IGcMem>(float v[R], uint tid, uint tx, uint ty, out float out_row) {
    GroupMemoryBarrierWithGroupSync();   // the memory free
    [ForceUnroll] for (int i = 0; i < R; ++i) Mem::st((ty * uint(R) + uint(i)) * uint(CTX) + tx, v[i]);
    GroupMemoryBarrierWithGroupSync();
    out_row = 0.0f;
    if (tid < uint(CP)) {
        [ForceUnroll] for (int k = 0; k < CTX; ++k) out_row += Mem::ld(tid * uint(CTX) + uint(k));
    }
}

#endif

// topo_gto_cube.slang
//
// The octree levels of the certified critical-point search on a GTO electron density
// (md_topo_compute_extremum_graph_gto_gpu), in fp32.
//
// Each cube c +- h is decided as the CPU reference decides it (see md_topo.c): a rigorous upper
// bound of rho below rho_min, exclusion of grad rho along a few directions, the Newton-distance test,
// or Krawczyk on the cube and on the 1.5x inflated cube. What changes in fp32 is that every quantity
// a test uses carries a bound on its own rounding error, and a test decides only with those bounds
// added. The error model assumes only what any conforming GPU guarantees:
//   * + - * correctly rounded; FMA contraction and reassociation only change which a priori bound
//     gamma_k * sum|terms| applies, and the bounds below hold for any order of every sum;
//   * exp(-t) relative error (6 + 8.1 t) u, u = 2^-24, including the rounding of t. The exp below is
//     computed here from + - * only (Cody-Waite + polynomial, ~6u), so this holds on every backend;
//   * denormals may be flushed, covered by an absolute floor per sum.
// Each AO channel carries its value V, its magnitude MV (the same expression on |terms|) and an error
// bound EV; downstream errors follow from |D| * (EV + n u MV), folded to one value per derivative order.
// The bounds and decisions were validated against double evaluations by an fp32 C emulation
// (Claude outputs/cp_prototype/gpu_fp32/cpg32_experiment.c in the viamd tree).
//
// WORK UNIT. A batch is the 8 children of one cube that split (the root lattice is grouped the same
// way, with a mask for the children that exist). The children share one local AO list, screened by
// the host against the union of their inflated cubes (extra AOs only add exact contributions), and
// so one sub-block of D. A chunk of batches runs as four launches, every intermediate in device
// memory (rows = local AOs of all batches in the chunk, columns = channel * 8 + child):
//
//   ao_main     one group per batch. Per local shell, primitive, axis and child side the 1D factors
//               of the centre values and of the remainder sups (cube and inflated cube) are computed
//               once into groupshared memory; they serve every Cartesian component and the 4 children
//               on that side. Then per (AO, child) the products: centre values with magnitude and error
//               shadows, and the remainder vectors.                                        -> phi
//   gemm_main   32-row tiles: D * V and |D| * (error weights, remainder vectors), over the batch's
//               local AOs.                                                                 -> dv
//   epi_main    one group per batch: DOT[a][b] = V_a . D V_b with its error sums (which depend on a
//               only through its derivative order), and the AO sums of the enclosures.     -> ep
//   decide_main one thread per cube: derivatives, enclosures, tests.                       -> outcome
//
// Every sum runs in a fixed order, so a cube's result depends only on its batch, never on the chunk.
//
// Outcome per cube (low 4 bits): 0 discard, 1 root (Krawczyk certified; the CPU polishes it in
// double), 2 split, 3 escalate (fp32 cannot decide: the gradient is inside its rounding noise and the
// cube is too small for a Krawczyk proof at that noise, or h <= h_min; the CPU takes the cube and its
// subtree in double). Bit 4: the inflated enclosure decided or was tried. Bits 8..15: which test.

#include "../md_gpu.slang"

#define LMAX    4
#define KI      (LMAX + 8)
#define KA      5
#define KM      (KI + KA + 1)
#define UF      5.9604644775390625e-8

// phi channels per (row, child); column = channel * 8 + child
#define PC_V    0       // 20 centre values                              (GEMM with D)
#define PC_WV   20      // 4  error weights of the D-products, per order  (GEMM with |D|)
#define PC_E    24      // 6  remainder terms E[0..5]                     (GEMM with |D|)
#define PC_EI   30      // 4  inflated remainder terms E[1], E[3..5]      (GEMM with |D|)
#define PC_MV   34      // 4  magnitudes, per order
#define PC_EV   38      // 4  value errors, per order
#define PC_E2   42      // 9  remainder terms E[6..14]
#define PC_EI2  51      // 6  inflated remainder terms E[9..14]
#define PHI_CH  57
#define PHI_W   (PHI_CH * 8)
// dv channels: 0..19 D V, 20..23 |D| WV, 24..29 |D| E[0..5], 30..33 |D| Ei
#define GEMM_CH 34
#define DV_W    (GEMM_CH * 8)
#define GEMM_ABS_COL (20 * 8)   // first column multiplied by |D| (= 5 groups of 32)

// ep fields per cube, structure of arrays: field f of cube c at ep[f * 8 * num_batches + c]
#define EP_DOT  0       // 200 DOT[a][b]
#define EP_EO   200     // 40  raw error sums per (order of a, b)
#define EP_PART 240     // 11  enclosure AO sums: rho (2), gradient (3), Hessian (6)
#define EP_PI   251     // 6   Hessian AO sums of the inflated cube
#define EP_W    260

#define OUT_DISCARD  0u
#define OUT_ROOT     1u
#define OUT_SPLIT    2u
#define OUT_ESCALATE 3u
#define OUT_INFLATED 16u

#define WHY_RHO      1u
#define WHY_EXCL     2u
#define WHY_NEWT     3u
#define WHY_KOUT     4u
#define WHY_KOUT_I   5u
#define WHY_ROOT     6u
#define WHY_SPLIT    7u
#define WHY_NOISE    8u
#define WHY_HMIN     9u

struct GtoShell {
    float4 A_r2;                // centre (Bohr), squared screening radius (unused here)
    uint4  info;                // prim_offset, num_prims, ao_offset, l | (ncart << 8)
};

struct Batch {
    float4 P_h;                 // parent centre, child half-width h; children at P +- h (all exact)
    uint4  info;                // child mask (bit s: x side s&1, y side s&2, z side s&4), local AOs n,
                                // first row, first local shell in batch_shells
    uint4  info2;               // number of local shells, -, -, -
};

struct RootArgs {
    uint  num_batches;
    uint  num_shells;
    uint  num_ao;
    uint  num_row_tiles;
    float eps;                  // rho_min
    float h_min;
    float tail_rho;
    float tail_g;
    float tail_H;
    uint  flags;                // bit 0: write debug
    uint  _pad0;
    uint  _pad1;
    Batch*    batches;
    uint*     batch_shells;     // global shell index per local shell, ascending per batch
    uint*     row_ao;           // global AO index per row
    uint2*    row_tiles;        // (batch, first local row) per 32-row tile of the GEMM
    GtoShell* shells;
    float*    alpha;
    float*    coeff;
    float*    ao_nrm;
    uint*     ao_ijk;           // i | j << 8 | k << 16
    float*    D;                // num_ao * num_ao
    float*    kappa;            // [KI][KA][KM]
    float*    phi;              // rows * PHI_W
    float*    dv;               // rows * DV_W
    float*    ep;               // EP_W * 8 * num_batches, by field
    uint*     outcome;          // 8 * num_batches
    float*    debug;            // optional, 16 floats per cube
};
MD_SHADER_ARGS(RootArgs);
#define args (*md_root.args)

// ------------------------------------------------------------------------------------ multi-indices
// Order 0..4 enumerated as in md_topo.c: by order, then i descending, then j descending.

int mi_index(int i, int j, int k) {
    const int o = i + j + k;
    const int base = o * (o + 1) * (o + 2) / 6;
    return base + (o - i) * (o - i + 1) / 2 + (o - i - j);
}
int M1(int a) { return 1 + a; }
int M2(int a, int b) { return mi_index(int(a == 0) + int(b == 0), int(a == 1) + int(b == 1), int(a == 2) + int(b == 2)); }
int M3(int a, int b, int c) {
    return mi_index(int(a == 0) + int(b == 0) + int(c == 0), int(a == 1) + int(b == 1) + int(c == 1), int(a == 2) + int(b == 2) + int(c == 2));
}
int M4(int a, int b, int c, int d) {
    return mi_index(int(a == 0) + int(b == 0) + int(c == 0) + int(d == 0), int(a == 1) + int(b == 1) + int(c == 1) + int(d == 1),
                    int(a == 2) + int(b == 2) + int(c == 2) + int(d == 2));
}

static const uint MI[35] = {
    0x000000, 0x000001, 0x000100, 0x010000, 0x000002, 0x000101, 0x010001, 0x000200, 0x010100, 0x020000,
    0x000003, 0x000102, 0x010002, 0x000201, 0x010101, 0x020001, 0x000300, 0x010200, 0x020100, 0x030000,
    0x000004, 0x000103, 0x010003, 0x000202, 0x010102, 0x020002, 0x000301, 0x010201, 0x020101, 0x030001,
    0x000400, 0x010300, 0x020200, 0x030100, 0x040000 };
int3 mi_of(int m) { const uint v = MI[m]; return int3(int(v & 0xFFu), int((v >> 8) & 0xFFu), int((v >> 16) & 0xFFu)); }
int order_of(int m) { return m == 0 ? 0 : (m < 4 ? 1 : (m < 10 ? 2 : 3)); }
int ab_index(int a, int b) { const int lo = min(a, b), hi = max(a, b); return lo == 0 ? hi : (lo == 1 ? 2 + hi : 5); }

float kappa(int i, int a, int m) { return args.kappa[(i * KA + a) * KM + m]; }

float powi(float x, int n) { float r = 1.0; for (int k = 0; k < n; ++k) r *= x; return r; }

// exp(-t) for t >= 0 from + - * only: Cody-Waite reduction and a degree-7 polynomial, relative error
// about 6u (well inside the (6 + 4t) u budget even if a compiler reassociates the reduction).
float exp_neg(float t) {
    if (t > 87.0) return 0.0;
    const float y = -t;
    const float kf = round(y * 1.44269504088896341);
    const float r = (y - kf * 0.693145751953125) - kf * 1.42860682030941723212e-6;
    float p = 1.0 / 5040.0;
    p = p * r + 1.0 / 720.0;
    p = p * r + 1.0 / 120.0;
    p = p * r + 1.0 / 24.0;
    p = p * r + 1.0 / 6.0;
    p = p * r + 0.5;
    p = p * r + 1.0;
    p = p * r + 1.0;
    const int k = int(kf);
    if (k < -126) return 0.0;
    return p * asfloat(uint(k + 127) << 23);
}

uint3 ao_ijk(uint ao) { const uint v = args.ao_ijk[min(ao, args.num_ao - 1u)]; return uint3(v & 0xFFu, (v >> 8) & 0xFFu, (v >> 16) & 0xFFu); }

// ================================================================================== ao_main
// One group per batch. Local shells are taken in chunks whose 1D tables fit TAB_CAP floats (the host
// guarantees that every single shell fits).
//
// Table entry per (shell, primitive p, axis, side), te(l) = 18 (l + 1) + 1 floats, l1 = l + 1:
//   [a l1 + i]             a = 0..3: d^a/dx^a [x^i exp(-al x^2)] at the child centre      (i = 0..l)
//   [4 l1 + a l1 + i]      the same on |terms| (magnitude shadow)
//   [8 l1]                 rounding term of this axis: (2l + 20 + 8.1 t) u
//   [8 l1 + 1 + 5 i + a]   a = 0..4: sup over the child interval of |d^a/dx^a [x^i exp(-al x^2)]|
//   [13 l1 + 1 + 5 i + a]  the same over the inflated child interval
// at offset ((p * 3 + axis) * 2 + side) * te(l) in the shell's block.

#define AO_WG   128
#define TAB_CAP 4096
#define CS_MAX  64

groupshared float gs_tab[TAB_CAP];
groupshared uint  gs_cs_shell[CS_MAX];
groupshared uint  gs_cs_tab[CS_MAX + 1];
groupshared uint  gs_cs_item[CS_MAX + 1];
groupshared uint  gs_cs_prod[CS_MAX + 1];
groupshared uint  gs_cs_row[CS_MAX];
groupshared uint  gs_ncs;

uint te_size(uint l) { return 18u * (l + 1u) + 1u; }

// sups of |d^a/dx^a [x^i exp(-al x^2)]| over [x0, x1], a = 0..4, for i = 0..l, into gs_tab[base + 5 i + a].
// The table of x^m exp(-al x^2) at its maximiser over the interval does not depend on i.
void sup_tab(uint base, uint l, float x0, float x1, float al) {
    x0 -= 4.0 * UF * abs(x0) + 1e-30;    // the computed interval must contain the exact one
    x1 += 4.0 * UF * abs(x1) + 1e-30;
    const float ax0 = abs(x0), ax1 = abs(x1);
    const float tlo = (x0 <= 0.0 && x1 >= 0.0) ? 0.0 : min(ax0, ax1);
    const float thi = max(ax0, ax1);
    float Sm[LMAX + KA], eS[LMAX + KA];
    [ForceUnroll] for (int m = 0; m < LMAX + KA; ++m) {
        Sm[m] = 0.0; eS[m] = 0.0;
        if (m <= int(l) + 4) {
            float ts = sqrt(float(m) / (2.0 * al));
            ts = clamp(ts, tlo, thi);
            const float t = al * ts * ts;
            Sm[m] = powi(ts, m) * exp_neg(t);
            eS[m] = (float(2 * m + 4) + 6.0 + 8.1 * t) * UF;
        }
    }
    for (int i = 0; i <= int(l); ++i) {
        {
            for (int a = 0; a < KA; ++a) {
                float acc = 0.0;
                [ForceUnroll] for (int m = 0; m < LMAX + KA; ++m) {
                    const float kk = m <= i + a ? kappa(i, a, m) : 0.0;
                    if (kk != 0.0) {
                        const int kp = (m - i + a) / 2;
                        acc += abs(kk) * powi(al, kp) * Sm[m] * (1.0 + eS[m] + float(kp + 2) * UF);
                    }
                }
                gs_tab[base + uint(i * KA + a)] = acc * (1.0 + float(i + a + 3) * UF);
            }
        }
    }
}

// One table entry: centre factors at x (child centre - atom), sups over [x0, x1] and [x0i, x1i].
void tab_entry(uint base, uint l, float x, float x0, float x1, float x0i, float x1i, float al) {
    const uint l1 = l + 1u;
    const float t = al * x * x;
    const float E = exp_neg(t);
    float xp[LMAX + 4], axp[LMAX + 4];
    xp[0] = 1.0; axp[0] = 1.0;
    [ForceUnroll] for (int m = 1; m < LMAX + 4; ++m) { xp[m] = xp[m - 1] * x; axp[m] = axp[m - 1] * abs(x); }
    for (int i = 0; i <= int(l); ++i) {
        {
            for (int a = 0; a < 4; ++a) {
                float v = 0.0, mg = 0.0;
                [ForceUnroll] for (int m = 0; m < LMAX + 4; ++m) {
                    const float k = m <= i + a ? kappa(i, a, m) : 0.0;
                    if (k != 0.0) {
                        const float ka = k * powi(al, (m - i + a) / 2);
                        v += ka * xp[m];
                        mg += abs(ka) * axp[m];
                    }
                }
                gs_tab[base + uint(a) * l1 + uint(i)] = v * E;
                gs_tab[base + 4u * l1 + uint(a) * l1 + uint(i)] = mg * E;
            }
        }
    }
    gs_tab[base + 8u * l1] = (float(2u * l + 14u) + 6.0 + 8.1 * t) * UF;
    sup_tab(base + 8u * l1 + 1u, l, x0, x1, al);
    sup_tab(base + 13u * l1 + 1u, l, x0i, x1i, al);
}

// Remainder vector from the sups S[1..34] of the products (cube of half-width h).
// E[0] E1, [1] E2, [2] E3, [3..5] Ea2, [6..8] Ea3, [9..14] Eab2 (a <= b).
void rem_fold(float S[35], float h, uint np, out float E[15]) {
    const float fs = 1.0 + float(np + 8u) * UF;
    float e1 = 0.0, e2 = 0.0, e3 = 0.0;
    float ea2[3] = { 0.0, 0.0, 0.0 }, ea3[3] = { 0.0, 0.0, 0.0 }, eab2[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    [ForceUnroll] for (int p = 0; p < 3; ++p) {
        e1 += S[M1(p)];
        [ForceUnroll] for (int q = 0; q < 3; ++q) { e2 += S[M2(p, q)]; [ForceUnroll] for (int r = 0; r < 3; ++r) e3 += S[M3(p, q, r)]; }
    }
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int p = 0; p < 3; ++p) [ForceUnroll] for (int q = 0; q < 3; ++q) {
        ea2[a] += S[M3(a, p, q)];
        [ForceUnroll] for (int r = 0; r < 3; ++r) ea3[a] += S[M4(a, p, q, r)];
    }
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int b = 0; b < 3; ++b) {
        if (b >= a) { [ForceUnroll] for (int p = 0; p < 3; ++p) [ForceUnroll] for (int q = 0; q < 3; ++q) eab2[ab_index(a, b)] += S[M4(a, b, p, q)]; }
    }
    const float f = fs * (1.0 + 48.0 * UF), h2 = h * h, h3 = h2 * h;
    E[0] = h * e1 * f; E[1] = 0.5 * h2 * e2 * f; E[2] = h3 * (1.0 / 6.0) * e3 * f;
    [ForceUnroll] for (int a = 0; a < 3; ++a) { E[3 + a] = 0.5 * h2 * ea2[a] * f; E[6 + a] = h3 * (1.0 / 6.0) * ea3[a] * f; }
    [ForceUnroll] for (int t = 0; t < 6; ++t) E[9 + t] = 0.5 * h2 * eab2[t] * f;
}

// The 1D factors of one (AO, child) for primitive p: table entries per axis.
struct AxisBases { uint x, y, z; };

// Products of local AO (shell cs of the chunk, component ci) for child k -> phi row 'row'.
void ao_products(uint cs, uint ci, uint k, uint row, float h, uint n) {
    const GtoShell sh = args.shells[min(gs_cs_shell[cs], args.num_shells - 1u)];
    const uint l = sh.info.w & 0xFFu, l1 = l + 1u, np = sh.info.y;
    const uint te = te_size(l);
    const uint ao = min(sh.info.z + min(ci, 14u), args.num_ao - 1u);
    const int3 e = int3(ao_ijk(ao));
    const float nrm = args.ao_nrm[ao];
    const uint sx = k & 1u, sy = (k >> 1) & 1u, sz = (k >> 2) & 1u;
    const uint tab0 = gs_cs_tab[cs];
    float* out = args.phi + row * PHI_W + k;

    // ---- centre values, magnitudes, errors
    float V[20], MV[20], EV[20];
    [ForceUnroll] for (int m = 0; m < 20; ++m) { V[m] = 0.0; MV[m] = 0.0; EV[m] = 0.0; }
    for (uint p = 0; p < np; ++p) {
        const float cfn = args.coeff[sh.info.x + p] * nrm, acfn = abs(cfn);
        const uint bx = tab0 + ((p * 3u + 0u) * 2u + sx) * te;
        const uint by = tab0 + ((p * 3u + 1u) * 2u + sy) * te;
        const uint bz = tab0 + ((p * 3u + 2u) * 2u + sz) * te;
        float gx[4], gy[4], gz[4], Gx[4], Gy[4], Gz[4];
        [ForceUnroll] for (int a = 0; a < 4; ++a) {
            gx[a] = gs_tab[bx + uint(a) * l1 + uint(e.x)]; Gx[a] = gs_tab[bx + (4u + uint(a)) * l1 + uint(e.x)];
            gy[a] = gs_tab[by + uint(a) * l1 + uint(e.y)]; Gy[a] = gs_tab[by + (4u + uint(a)) * l1 + uint(e.y)];
            gz[a] = gs_tab[bz + uint(a) * l1 + uint(e.z)]; Gz[a] = gs_tab[bz + (4u + uint(a)) * l1 + uint(e.z)];
        }
        float eps = 5.0 * UF;
        eps += gs_tab[bx + 8u * l1];
        eps += gs_tab[by + 8u * l1];
        eps += gs_tab[bz + 8u * l1];
        [ForceUnroll] for (int m = 0; m < 20; ++m) {
            const int3 d = mi_of(m);
            const float term = cfn * gx[d.x] * gy[d.y] * gz[d.z];
            const float mag = acfn * Gx[d.x] * Gy[d.y] * Gz[d.z];
            V[m] += term; MV[m] += mag; EV[m] += mag * eps;
        }
    }
    float WVo[4] = { 0.0, 0.0, 0.0, 0.0 }, MVo[4] = { 0.0, 0.0, 0.0, 0.0 }, EVo[4] = { 0.0, 0.0, 0.0, 0.0 };
    [ForceUnroll] for (int m = 0; m < 20; ++m) {
        const float ev = (EV[m] + float(np + 1u) * UF * MV[m]) * 1.001 + 1e-36;
        const int o = order_of(m);
        WVo[o] = max(WVo[o], ev + float(n + 2u) * UF * MV[m]);
        MVo[o] = max(MVo[o], MV[m]);
        EVo[o] = max(EVo[o], ev);
    }
    [ForceUnroll] for (int m = 0; m < 20; ++m) out[(PC_V + m) * 8] = V[m];
    [ForceUnroll] for (int o = 0; o < 4; ++o) { out[(PC_WV + o) * 8] = WVo[o]; out[(PC_MV + o) * 8] = MVo[o]; out[(PC_EV + o) * 8] = EVo[o]; }

    // ---- remainder vectors, cube and inflated cube
    [ForceUnroll] for (int pass = 0; pass < 2; ++pass) {
        float S[35];
        [ForceUnroll] for (int m = 0; m < 35; ++m) S[m] = 0.0;
        const uint boff = 8u * l1 + 1u + uint(pass) * 5u * l1;
        for (uint p = 0; p < np; ++p) {
            const float w = abs(args.coeff[sh.info.x + p]) * nrm;
            const uint bx = tab0 + ((p * 3u + 0u) * 2u + sx) * te + boff;
            const uint by = tab0 + ((p * 3u + 1u) * 2u + sy) * te + boff;
            const uint bz = tab0 + ((p * 3u + 2u) * 2u + sz) * te + boff;
            float Bx[KA], By[KA], Bz[KA];
            [ForceUnroll] for (int a = 0; a < KA; ++a) {
                Bx[a] = gs_tab[bx + uint(e.x) * KA + uint(a)];
                By[a] = gs_tab[by + uint(e.y) * KA + uint(a)];
                Bz[a] = gs_tab[bz + uint(e.z) * KA + uint(a)];
            }
            // the inflated cube only needs orders 2..4 (E2, Ea2, Eab2)
            [ForceUnroll] for (int m = 1; m < 35; ++m) {
                if (pass == 0 || m >= 4) {
                    const int3 d = mi_of(m);
                    S[m] += w * Bx[d.x] * By[d.y] * Bz[d.z];
                }
            }
        }
        float E[15];
        if (pass == 0) {
            rem_fold(S, h, np, E);
            [ForceUnroll] for (int t = 0; t < 6; ++t) out[(PC_E + t) * 8] = E[t];
            [ForceUnroll] for (int t = 0; t < 9; ++t) out[(PC_E2 + t) * 8] = E[6 + t];
        } else {
            rem_fold(S, 1.5 * h, np, E);
            out[(PC_EI + 0) * 8] = E[1];
            [ForceUnroll] for (int a = 0; a < 3; ++a) out[(PC_EI + 1 + a) * 8] = E[3 + a];
            [ForceUnroll] for (int t = 0; t < 6; ++t) out[(PC_EI2 + t) * 8] = E[9 + t];
        }
    }
}

[shader("compute")]
[numthreads(AO_WG, 1, 1)]
void ao_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    const uint b = gid.x;
    if (b >= args.num_batches) return;      // uniform per group
    const Batch bt = args.batches[b];
    const float3 P = bt.P_h.xyz;
    const float h = bt.P_h.w;
    const uint n = bt.info.y, row0 = bt.info.z, sh0 = bt.info.w, nsh = bt.info2.x;

    uint next = 0, rowc = 0;                // thread 0's cursors
    for (;;) {
        // ---- chunk layout (thread 0). Written only here, read by all after the barrier below and
        //      not again after the barrier that ends the iteration.
        if (tid == 0) {
            uint s = next, tab = 0, item = 0, prod = 0, row = rowc, ncs = 0;
            while (s < nsh && ncs < CS_MAX) {
                const uint g = min(args.batch_shells[sh0 + s], args.num_shells - 1u);
                const uint4 inf = args.shells[g].info;
                const uint l = inf.w & 0xFFu, nc = inf.w >> 8, np = inf.y;
                const uint sz = np * 6u * te_size(l);
                if (tab + sz > TAB_CAP) break;
                gs_cs_shell[ncs] = g; gs_cs_tab[ncs] = tab; gs_cs_item[ncs] = item; gs_cs_prod[ncs] = prod; gs_cs_row[ncs] = row;
                tab += sz; item += np * 6u; prod += nc * 8u; row += nc; ncs++; s++;
            }
            gs_cs_tab[ncs] = tab; gs_cs_item[ncs] = item; gs_cs_prod[ncs] = prod;
            gs_ncs = ncs;                   // 0: done (or a shell that does not fit, which the host never sends)
            next = s; rowc = row;
        }
        GroupMemoryBarrierWithGroupSync();
        const uint ncs = gs_ncs;
        if (ncs == 0) break;                // uniform

        // ---- 1D tables: (shell, primitive, axis, side)
        const uint nitem = gs_cs_item[ncs];
        for (uint f = tid; f < nitem; f += AO_WG) {
            uint cs = 0;
            while (cs + 1u < ncs && gs_cs_item[cs + 1u] <= f) cs++;
            const uint r = f - gs_cs_item[cs];
            const uint p = r / 6u, ax = (r % 6u) / 2u, side = r % 2u;
            const GtoShell sh = args.shells[min(gs_cs_shell[cs], args.num_shells - 1u)];
            const uint l = sh.info.w & 0xFFu;
            const float al = args.alpha[sh.info.x + min(p, sh.info.y - 1u)];
            const float A = sh.A_r2[ax];
            const float c = P[ax] + (side != 0u ? h : -h);
            const float hi = 1.5 * h;
            tab_entry(gs_cs_tab[cs] + r * te_size(l), l, c - A, (c - h) - A, (c + h) - A, (c - hi) - A, (c + hi) - A, al);
        }
        GroupMemoryBarrierWithGroupSync();

        // ---- products: (AO, child), children fastest
        const uint nprod = gs_cs_prod[ncs];
        for (uint f = tid; f < nprod; f += AO_WG) {
            uint cs = 0;
            while (cs + 1u < ncs && gs_cs_prod[cs + 1u] <= f) cs++;
            const uint r = f - gs_cs_prod[cs];
            const uint ci = r / 8u, k = r % 8u;
            ao_products(cs, ci, k, row0 + gs_cs_row[cs] + ci, h, n);
        }
        GroupMemoryBarrierWithGroupSync();
    }
}

// ================================================================================== gemm_main
// dv[row][c] = sum_j X[row][j] phi[j][c] over the batch's local AOs j, X = D for the 160 columns of the
// centre values and |D| for the 112 after them. One group per 32-row tile computes all 272 columns, so
// each D tile is gathered once: 32 x 8 threads, thread (tx, ty) owns rows 4 ty .. 4 ty + 3 and columns
// tx + 32 j, j = 0..8 (j < 5 with D, j >= 5 with |D|). The sum over j runs in increasing j, so the
// result does not depend on the tiling.

#define GT  32                  // rows per tile
#define GKC 16                  // j per chunk
#define GCJ 9                   // column groups of 32: 9 * 32 = 288 >= DV_W
#define GEMM_WG 256

groupshared float gs_A[GKC][GT];
groupshared float gs_B[GKC][GCJ * 32];

[shader("compute")]
[numthreads(32, 8, 1)]
void gemm_main(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {
    const uint rt = gid.x;
    if (rt >= args.num_row_tiles) return;   // uniform per group
    const uint2 rtd = args.row_tiles[rt];
    const Batch bt = args.batches[min(rtd.x, args.num_batches - 1u)];
    const uint n = bt.info.y, base = bt.info.z, r0 = rtd.y;
    const uint tx = gtid.x, ty = gtid.y;
    const uint N = args.num_ao;

    float acc[4][GCJ];
    [ForceUnroll] for (int i = 0; i < 4; ++i) [ForceUnroll] for (int j = 0; j < GCJ; ++j) acc[i][j] = 0.0;

    for (uint k0 = 0; k0 < n; k0 += GKC) {
        // D tile: 32 rows x 16 j, 2 per thread, j fastest (consecutive AOs are mostly consecutive in D)
        [ForceUnroll] for (int q = 0; q < 2; ++q) {
            const uint e = tid + uint(GEMM_WG * q);
            const uint r = e / GKC, kk = e % GKC;
            float d = 0.0;
            if (r0 + r < n && k0 + kk < n) {
                const uint ai = min(args.row_ao[base + r0 + r], N - 1u);
                const uint aj = min(args.row_ao[base + k0 + kk], N - 1u);
                d = args.D[ai * N + aj];
            }
            gs_A[kk][r] = d;
        }
        // phi tile: 16 j x 288 columns, 18 per thread, columns fastest (coalesced rows)
        [ForceUnroll] for (int q = 0; q < 18; ++q) {
            const uint e = tid + uint(GEMM_WG * q);
            const uint kk = e / (GCJ * 32), c = e % (GCJ * 32);
            float v = 0.0;
            if (k0 + kk < n && c < DV_W) v = args.phi[(base + k0 + kk) * PHI_W + c];
            gs_B[kk][c] = v;
        }
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll] for (int kk = 0; kk < GKC; ++kk) {
            float a[4], aa[4], bb[GCJ];
            [ForceUnroll] for (int i = 0; i < 4; ++i) { a[i] = gs_A[kk][ty * 4u + uint(i)]; aa[i] = abs(a[i]); }
            [ForceUnroll] for (int j = 0; j < GCJ; ++j) bb[j] = gs_B[kk][tx + 32u * uint(j)];
            [ForceUnroll] for (int i = 0; i < 4; ++i) [ForceUnroll] for (int j = 0; j < GCJ; ++j) acc[i][j] += (j < 5 ? a[i] : aa[i]) * bb[j];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    [ForceUnroll] for (int i = 0; i < 4; ++i) {
        const uint row = r0 + ty * 4u + uint(i);
        if (row >= n) continue;
        [ForceUnroll] for (int j = 0; j < GCJ; ++j) {
            const uint col = tx + 32u * uint(j);
            if (col < DV_W) args.dv[(base + row) * DV_W + col] = acc[i][j];
        }
    }
}

// ================================================================================== epi_main
// One group per batch, 16 threads per child. Rows are staged IT at a time into groupshared memory.
//   DOT[a][b] = sum_i V_i[a] (D V)_i[b], a = 0..19, b = 0..9, and its error sum, which depends on a only
//   through its derivative order o: e[o][b] = sum_i EV_i[o] (|DV_i[b]| + EDV_i[ob]) + MV_i[o] EDV_i[ob]
//                                           + (n + 1) u MV_i[o] |DV_i[b]|.
//   The enclosure's AO sums (cube and inflated cube), one (row, child) per thread of the first 64.

#define EPI_WG 128
#define IT     8
#define T_V    0        // 20
#define T_MV   20       // 4
#define T_EV   24       // 4
#define T_DV   28       // 20
#define T_EDV  48       // 4   final: raw * safe + n 1e-36
#define T_Q    52       // 6   raw
#define T_QI   58       // 4   raw
#define T_E    62       // 15
#define T_EI   77       // 10  inflated E[1], E[3..5], E[9..14]
#define TCH    87

groupshared float gs_t[IT * TCH * 8];

// source of tile channel ch: phi channel (< 64) or 64 + dv channel
static const uint TSRC[TCH] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,              // V
    34, 35, 36, 37,                                                                    // MV
    38, 39, 40, 41,                                                                    // EV
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,    // DV
    84, 85, 86, 87,                                                                    // EDV
    88, 89, 90, 91, 92, 93,                                                            // Q
    94, 95, 96, 97,                                                                    // QI
    24, 25, 26, 27, 28, 29, 42, 43, 44, 45, 46, 47, 48, 49, 50,                        // E[0..14]
    30, 31, 32, 33, 51, 52, 53, 54, 55, 56 };                                          // EI

[shader("compute")]
[numthreads(EPI_WG, 1, 1)]
void epi_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    const uint b = gid.x;
    if (b >= args.num_batches) return;      // uniform per group
    const Batch bt = args.batches[b];
    const float h = bt.P_h.w, hi = 1.5 * h;
    const uint n = bt.info.y, base = bt.info.z;
    const float fn = float(n);
    const float safe = 1.001 + (fn + 64.0) * UF;
    const float h2 = h * h, hi2 = hi * hi;

    // DOT role: child kd, a = sd and sd + 16, error pairs sd, sd + 16, sd + 32
    const uint kd = tid >> 4, sd = tid & 15u;
    float s0[10], s1[10], eo[3];
    [ForceUnroll] for (int j = 0; j < 10; ++j) { s0[j] = 0.0; s1[j] = 0.0; }
    [ForceUnroll] for (int q = 0; q < 3; ++q) eo[q] = 0.0;
    // row role (tid < 64): child kr, row il of each tile
    const uint kr = tid & 7u, ilr = tid >> 3;
    float part[17];
    [ForceUnroll] for (int t = 0; t < 17; ++t) part[t] = 0.0;

    for (uint i0 = 0; i0 < n; i0 += IT) {
        // ---- stage rows i0 .. i0 + IT - 1
        for (uint e = tid; e < IT * TCH * 8u; e += EPI_WG) {
            const uint kk = e & 7u, ch = (e >> 3) % TCH, il = (e >> 3) / TCH;
            float v = 0.0;
            if (i0 + il < n) {
                const uint row = base + i0 + il, src = TSRC[ch];
                if (src < 64u) v = args.phi[row * PHI_W + src * 8u + kk];
                else           v = args.dv[row * DV_W + (src - 64u) * 8u + kk];
                if (ch >= T_EDV && ch < T_EDV + 4) v = v * safe + fn * 1e-36;
            }
            gs_t[e] = v;
        }
        GroupMemoryBarrierWithGroupSync();

        // ---- DOT and its error sums, rows in increasing order
        for (uint il = 0; il < IT; ++il) {
            if (i0 + il >= n) break;
            const uint tb = il * TCH * 8u + kd;
            const float va0 = gs_t[tb + (T_V + sd) * 8u];
            const float va1 = sd < 4u ? gs_t[tb + (T_V + 16u + sd) * 8u] : 0.0;
            [ForceUnroll] for (int j = 0; j < 10; ++j) {
                const float dvb = gs_t[tb + (T_DV + uint(j)) * 8u];
                s0[j] += va0 * dvb;
                s1[j] += va1 * dvb;
            }
            [ForceUnroll] for (int q = 0; q < 3; ++q) {
                const uint idx = sd + 16u * uint(q);
                if (idx < 40u) {
                    const uint o = idx / 10u, bb = idx % 10u, ob = uint(order_of(int(bb)));
                    const float ea = gs_t[tb + (T_EV + o) * 8u], ma = gs_t[tb + (T_MV + o) * 8u];
                    const float dv = gs_t[tb + (T_DV + bb) * 8u], edv = gs_t[tb + (T_EDV + ob) * 8u], adv = abs(dv);
                    eo[q] += ea * (adv + edv) + ma * edv + (fn + 1.0) * UF * ma * adv;
                }
            }
        }

        // ---- enclosure AO sums, one (row, child) per thread
        if (tid < 64u && i0 + ilr < n) {
            const uint tb = ilr * TCH * 8u + kr;
            float ad[20];
            [ForceUnroll] for (int m = 0; m < 20; ++m) ad[m] = abs(gs_t[tb + uint(T_DV + m) * 8u]) + gs_t[tb + uint(T_EDV + order_of(m)) * 8u];
            float q[6], qi[4], e[15], ei[15];
            [ForceUnroll] for (int t = 0; t < 6; ++t) q[t] = gs_t[tb + uint(T_Q + t) * 8u] * safe;
            [ForceUnroll] for (int t = 0; t < 4; ++t) qi[t] = gs_t[tb + uint(T_QI + t) * 8u] * safe;
            [ForceUnroll] for (int t = 0; t < 15; ++t) { e[t] = gs_t[tb + uint(T_E + t) * 8u]; ei[t] = 0.0; }
            ei[1] = gs_t[tb + uint(T_EI + 0) * 8u];
            [ForceUnroll] for (int a = 0; a < 3; ++a) ei[3 + a] = gs_t[tb + uint(T_EI + 1 + a) * 8u];
            [ForceUnroll] for (int t = 0; t < 6; ++t) ei[9 + t] = gs_t[tb + uint(T_EI + 4 + t) * 8u];
            // inflated |D| sums in the cube's slots: qs[1] = E2, qs[3 + a] = Ea2
            float qs[6] = { 0.0, qi[0], 0.0, qi[1], qi[2], qi[3] };

            float sum1 = 0.0, sum2 = 0.0;
            [ForceUnroll] for (int p = 0; p < 3; ++p) { sum1 += ad[M1(p)]; [ForceUnroll] for (int r = 0; r < 3; ++r) sum2 += ad[M2(p, r)]; }
            const float W = ad[0] + h * sum1 + 0.5 * h2 * sum2;
            const float W1 = ad[0] + h * sum1;
            const float W1i = ad[0] + hi * sum1;
            float Wa[3], W1a[3], W1ai[3], W1ab[6], W1abi[6];
            [ForceUnroll] for (int a = 0; a < 3; ++a) {
                float t1 = 0.0, t2 = 0.0;
                [ForceUnroll] for (int p = 0; p < 3; ++p) { t1 += ad[M2(a, p)]; [ForceUnroll] for (int r = 0; r < 3; ++r) t2 += ad[M3(a, p, r)]; }
                Wa[a] = ad[M1(a)] + h * t1 + 0.5 * h2 * t2;
                W1a[a] = ad[M1(a)] + h * t1;
                W1ai[a] = ad[M1(a)] + hi * t1;
            }
            [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int bq = 0; bq < 3; ++bq) {
                if (bq >= a) {
                    float t1 = 0.0;
                    [ForceUnroll] for (int p = 0; p < 3; ++p) t1 += ad[M3(a, bq, p)];
                    W1ab[ab_index(a, bq)] = ad[M2(a, bq)] + h * t1;
                    W1abi[ab_index(a, bq)] = ad[M2(a, bq)] + hi * t1;
                }
            }
            part[0] += ad[0] * e[0];
            part[1] += e[0] * q[0];
            [ForceUnroll] for (int a = 0; a < 3; ++a) part[2 + a] += e[6 + a] * W + e[2] * Wa[a] + e[6 + a] * q[2];
            [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int bq = 0; bq < 3; ++bq) {
                if (bq >= a) {
                    const int t = ab_index(a, bq);
                    part[5 + t] += e[9 + t] * W1 + W1ab[t] * e[1] + e[9 + t] * q[1] + e[3 + a] * W1a[bq] + W1a[a] * e[3 + bq] + e[3 + a] * q[3 + bq];
                    part[11 + t] += ei[9 + t] * W1i + W1abi[t] * ei[1] + ei[9 + t] * qs[1] + ei[3 + a] * W1ai[bq] + W1ai[a] * ei[3 + bq] + ei[3 + a] * qs[3 + bq];
                }
            }
        }
        GroupMemoryBarrierWithGroupSync();
    }

    // ---- outputs
    const uint fs = 8u * args.num_batches;
    float* out = args.ep + b * 8u + kd;
    [ForceUnroll] for (int j = 0; j < 10; ++j) {
        out[(EP_DOT + sd * 10u + uint(j)) * fs] = s0[j];
        if (sd < 4u) out[(EP_DOT + (sd + 16u) * 10u + uint(j)) * fs] = s1[j];
    }
    [ForceUnroll] for (int q = 0; q < 3; ++q) {
        const uint idx = sd + 16u * uint(q);
        if (idx < 40u) out[(EP_EO + idx) * fs] = eo[q];
    }
    // AO sums: the 8 row lanes of each child, in lane order
    if (tid < 64u) {
        [ForceUnroll] for (int t = 0; t < 17; ++t) gs_t[tid * 17u + uint(t)] = part[t];
    }
    GroupMemoryBarrierWithGroupSync();
    if (tid < 8u) {
        float tot[17];
        [ForceUnroll] for (int t = 0; t < 17; ++t) tot[t] = 0.0;
        for (uint lane = 0; lane < 8u; ++lane) {
            const uint src = (lane * 8u + tid) * 17u;
            [ForceUnroll] for (int t = 0; t < 17; ++t) tot[t] += gs_t[src + uint(t)];
        }
        float* o2 = args.ep + b * 8u + tid;
        [ForceUnroll] for (int t = 0; t < 11; ++t) o2[(EP_PART + t) * fs] = tot[t];
        [ForceUnroll] for (int t = 0; t < 6; ++t) o2[(EP_PI + t) * fs] = tot[11 + t];
    }
}

// ================================================================================== decide_main

struct Eval {
    float rho, e_rho;
    float g[3], e_g[3];
    float A[3][3], e_A[3][3];
    float q_up[3], r_up[3], dH_up[3][3], rho_up;
    float dH_i[3][3];           // inflated cube
};

void eigen3(float M[3][3], out float val[3], out float vec[3][3]) {
    float a[3][3] = M;
    float v[3][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };
    for (int sweep = 0; sweep < 30; ++sweep) {
        const float off = abs(a[0][1]) + abs(a[0][2]) + abs(a[1][2]);
        const float scale = abs(a[0][0]) + abs(a[1][1]) + abs(a[2][2]);
        if (off <= 1e-30 || off <= 1e-7 * scale) break;
        [ForceUnroll] for (int p = 0; p < 2; ++p) [ForceUnroll] for (int q = 0; q < 3; ++q) {
            if (q > p && a[p][q] != 0.0) {
                const float theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const float t = (theta >= 0.0 ? 1.0 : -1.0) / (abs(theta) + sqrt(theta * theta + 1.0));
                const float cs = 1.0 / sqrt(t * t + 1.0), sn = t * cs;
                [ForceUnroll] for (int k = 0; k < 3; ++k) { const float x = a[k][p], y = a[k][q]; a[k][p] = cs * x - sn * y; a[k][q] = sn * x + cs * y; }
                [ForceUnroll] for (int k = 0; k < 3; ++k) { const float x = a[p][k], y = a[q][k]; a[p][k] = cs * x - sn * y; a[q][k] = sn * x + cs * y; }
                [ForceUnroll] for (int k = 0; k < 3; ++k) { const float x = v[k][p], y = v[k][q]; v[k][p] = cs * x - sn * y; v[k][q] = sn * x + cs * y; }
            }
        }
    }
    // sort ascending (3 elements, branch-free on the indices)
    float d[3] = { a[0][0], a[1][1], a[2][2] };
    int i0 = 0, i1 = 1, i2 = 2;
    if (d[i1] < d[i0]) { const int t = i0; i0 = i1; i1 = t; }
    if (d[i2] < d[i0]) { const int t = i0; i0 = i2; i2 = t; }
    if (d[i2] < d[i1]) { const int t = i1; i1 = i2; i2 = t; }
    const int idx[3] = { i0, i1, i2 };
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        val[i] = d[idx[i]];
        [ForceUnroll] for (int k = 0; k < 3; ++k) vec[k][i] = v[k][idx[i]];
    }
}

bool excl_dir(Eval ev, float h, float3 w) {
    float wg = 0.0, wgabs = 0.0, ewg = 0.0;
    [ForceUnroll] for (int k = 0; k < 3; ++k) { wg += w[k] * ev.g[k]; wgabs += abs(w[k] * ev.g[k]); ewg += abs(w[k]) * ev.e_g[k]; }
    const float lhs_lo = abs(wg) - ewg - 3.0 * UF * wgabs;
    float rhs = 0.0;
    [ForceUnroll] for (int j = 0; j < 3; ++j) {
        float Aw = 0.0, Awabs = 0.0, eAw = 0.0;
        [ForceUnroll] for (int k = 0; k < 3; ++k) { Aw += ev.A[j][k] * w[k]; Awabs += abs(ev.A[j][k] * w[k]); eAw += ev.e_A[j][k] * abs(w[k]); }
        rhs += (abs(Aw) + eAw + 3.0 * UF * Awabs) * h + abs(w[j]) * (ev.q_up[j] + ev.r_up[j]);
    }
    rhs *= 1.0 + 16.0 * UF;
    return lhs_lo > rhs;
}

struct Lin {
    bool ok;
    float Y[3][3], R[3][3], eR[3][3];
    float rhoR, zeta;
    float m[3], em[3];
};

Lin linear(Eval ev, float val[3], float vec[3][3]) {
    Lin L;
    L.ok = false; L.rhoR = 0.0; L.zeta = 0.0;
    [ForceUnroll] for (int i = 0; i < 3; ++i) { L.m[i] = 0.0; L.em[i] = 0.0; [ForceUnroll] for (int j = 0; j < 3; ++j) { L.Y[i][j] = 0.0; L.R[i][j] = 0.0; L.eR[i][j] = 0.0; } }
    const float amax = max(abs(val[0]), abs(val[2]));
    if (amax == 0.0) return L;
    [ForceUnroll] for (int i = 0; i < 3; ++i) if (!(abs(val[i]) > 1e-6 * amax)) return L;
    [ForceUnroll] for (int r = 0; r < 3; ++r) [ForceUnroll] for (int c = 0; c < 3; ++c) {
        float s = 0.0;
        [ForceUnroll] for (int k = 0; k < 3; ++k) s += vec[r][k] * vec[c][k] / val[k];
        L.Y[r][c] = s;
    }
    float Ynorm = 0.0;
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        float rs = 0.0, rowR = 0.0;
        [ForceUnroll] for (int j = 0; j < 3; ++j) {
            float s = (i == j) ? 1.0 : 0.0, sabs = (i == j) ? 1.0 : 0.0, yea = 0.0;
            [ForceUnroll] for (int k = 0; k < 3; ++k) { s -= L.Y[i][k] * ev.A[k][j]; sabs += abs(L.Y[i][k] * ev.A[k][j]); yea += abs(L.Y[i][k]) * ev.e_A[k][j]; }
            L.R[i][j] = s;
            L.eR[i][j] = (yea + 4.0 * UF * sabs) * 1.001;
            rowR += abs(s) + L.eR[i][j];
            rs += abs(L.Y[i][j]);
        }
        L.rhoR = max(L.rhoR, rowR * (1.0 + 4.0 * UF));
        Ynorm = max(Ynorm, rs);
    }
    if (!(L.rhoR < 0.9)) return L;
    L.zeta = Ynorm * L.rhoR / (1.0 - L.rhoR) * (1.0 + 16.0 * UF);
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        float s = 0.0, sabs = 0.0, e = 0.0;
        [ForceUnroll] for (int k = 0; k < 3; ++k) { s += L.Y[i][k] * ev.g[k]; sabs += abs(L.Y[i][k] * ev.g[k]); e += abs(L.Y[i][k]) * ev.e_g[k]; }
        L.m[i] = -s;
        L.em[i] = (e + 3.0 * UF * sabs) * 1.001;
    }
    L.ok = true;
    return L;
}

bool newton_dist(Eval ev, Lin L, float h) {
    float gabs = 0.0;
    [ForceUnroll] for (int k = 0; k < 3; ++k) gabs = max(gabs, abs(ev.g[k]) + ev.e_g[k]);
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        float rad = L.em[i] + L.zeta * 3.0 * gabs;
        [ForceUnroll] for (int j = 0; j < 3; ++j) rad += (abs(L.Y[i][j]) + L.zeta) * (ev.q_up[j] + ev.r_up[j]);
        rad *= 1.0 + 16.0 * UF;
        if (abs(L.m[i]) - rad > h * (1.0 + 2.0 * UF)) return true;
    }
    return false;
}

// +1 unique zero in c +- hk, -1 none, 0 undecided
int krawczyk(Lin L, float dH[3][3], float hk) {
    bool inside = true, disjoint = false;
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        float rad = 0.0;
        [ForceUnroll] for (int j = 0; j < 3; ++j) {
            float yd = 0.0;
            [ForceUnroll] for (int k = 0; k < 3; ++k) yd += abs(L.Y[i][k]) * dH[k][j];
            rad += abs(L.R[i][j]) + L.eR[i][j] + yd;
        }
        rad *= hk * (1.0 + 16.0 * UF);
        const float lo_out = abs(L.m[i]) - L.em[i] - rad;
        const float hi_in = (abs(L.m[i]) + L.em[i] + rad) * (1.0 + 4.0 * UF);
        if (!(hi_in < hk)) inside = false;
        if (lo_out > hk * (1.0 + 2.0 * UF)) disjoint = true;
    }
    return inside ? 1 : (disjoint ? -1 : 0);
}

#define DEC_WG 64

[shader("compute")]
[numthreads(DEC_WG, 1, 1)]
void decide_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) {
    const uint cube = gid.x * DEC_WG + tid;
    if (cube >= args.num_batches * 8u) return;
    const uint b = cube >> 3, k = cube & 7u;
    const Batch bt = args.batches[b];
    if (((bt.info.x >> k) & 1u) == 0u) { args.outcome[cube] = OUT_DISCARD; return; }
    const float h = bt.P_h.w, hi = 1.5 * h;
    const float h2 = h * h, h3 = h2 * h, hi2 = hi * hi;
    const uint n = bt.info.y;
    const float fn = float(n);
    const float safe = 1.001 + (fn + 64.0) * UF;
    const uint fs = 8u * args.num_batches;
    float* E = args.ep + cube;
    // the 40 error sums once; DOT read where used (consecutive cubes are consecutive in each field)
    float eo[40];
    [ForceUnroll] for (int f = 0; f < 40; ++f) eo[f] = E[(EP_EO + f) * fs] * safe + fn * 1e-36;
    #define DOT(a, bb)  E[(EP_DOT + (a) * 10 + (bb)) * fs]
    #define EDOT(a, bb) eo[order_of(a) * 10 + (bb)]

    // ---- centre derivatives and their errors
    Eval ev;
    ev.rho = DOT(0, 0); ev.e_rho = EDOT(0, 0);
    [ForceUnroll] for (int a = 0; a < 3; ++a) { ev.g[a] = 2.0 * DOT(M1(a), 0); ev.e_g[a] = 2.0 * EDOT(M1(a), 0); }
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int c = 0; c < 3; ++c) {
        const float x = DOT(M2(a, c), 0), y = DOT(M1(a), M1(c));
        ev.A[a][c] = 2.0 * (x + y);
        ev.e_A[a][c] = (2.0 * (EDOT(M2(a, c), 0) + EDOT(M1(a), M1(c))) + 4.0 * UF * (abs(x) + abs(y))) * 1.001;
    }
    // third derivatives enter only through sums of |T| + e_T
    float qT[3], lT[6];
    [ForceUnroll] for (int a = 0; a < 3; ++a) qT[a] = 0.0;
    [ForceUnroll] for (int t = 0; t < 6; ++t) lT[t] = 0.0;
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int c = 0; c < 3; ++c) [ForceUnroll] for (int d = 0; d < 3; ++d) {
        const float t0 = DOT(M3(a, c, d), 0), t1 = DOT(M2(a, c), M1(d)), t2 = DOT(M2(a, d), M1(c)), t3 = DOT(M2(c, d), M1(a));
        const float T = 2.0 * (t0 + t1 + t2 + t3);
        const float eT = (2.0 * (EDOT(M3(a, c, d), 0) + EDOT(M2(a, c), M1(d)) + EDOT(M2(a, d), M1(c)) + EDOT(M2(c, d), M1(a)))
                          + 8.0 * UF * (abs(t0) + abs(t1) + abs(t2) + abs(t3))) * 1.001;
        qT[a] += abs(T) + eT;                                 // sum over (c, d) of |T[a][c][d]|
        if (a <= c) lT[ab_index(a, c)] += abs(T) + eT;        // sum over d, a <= c
    }

    // ---- enclosure over the cube (and the Hessian over the inflated cube)
    [ForceUnroll] for (int a = 0; a < 3; ++a) {
        float d3u = 0.0, d4u = 0.0;
        [ForceUnroll] for (int p = 0; p < 3; ++p) [ForceUnroll] for (int qq = 0; qq < 3; ++qq) [ForceUnroll] for (int r = 0; r < 3; ++r) {
            const float x = DOT(M3(a, p, qq), M1(r)), y = DOT(M2(a, p), M2(qq, r));
            d3u += abs(x + y) + EDOT(M3(a, p, qq), M1(r)) + EDOT(M2(a, p), M2(qq, r)) + 2.0 * UF * (abs(x) + abs(y));
            [ForceUnroll] for (int s = 0; s < 3; ++s) d4u += abs(DOT(M3(a, p, qq), M2(r, s))) + EDOT(M3(a, p, qq), M2(r, s));
        }
        ev.r_up[a] = (h3 * d3u + 0.5 * h2 * h2 * d4u + 2.0 * E[(EP_PART + 2 + a) * fs] * safe + args.tail_g) * safe;
        ev.q_up[a] = 0.5 * h2 * qT[a] * safe;
    }
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int c = 0; c < 3; ++c) {
        if (c < a) continue;
        float qu = 0.0;
        [ForceUnroll] for (int p = 0; p < 3; ++p) [ForceUnroll] for (int qq = 0; qq < 3; ++qq) {
            const float x = DOT(M3(a, c, p), M1(qq)), y = DOT(M2(a, p), M2(c, qq));
            qu += abs(x + y) + EDOT(M3(a, c, p), M1(qq)) + EDOT(M2(a, p), M2(c, qq)) + 2.0 * UF * (abs(x) + abs(y));
        }
        const int t = ab_index(a, c);
        const float lu = lT[t];
        const float vu = (h * lu + 2.0 * h2 * qu + 2.0 * E[(EP_PART + 5 + t) * fs] * safe + args.tail_H) * safe;
        const float vi = (hi * lu + 2.0 * hi2 * qu + 2.0 * E[(EP_PI + t) * fs] * safe + args.tail_H) * safe;
        ev.dH_up[a][c] = vu; ev.dH_up[c][a] = vu;
        ev.dH_i[a][c] = vi; ev.dH_i[c][a] = vi;
    }
    ev.rho_up = (ev.rho + ev.e_rho + 2.0 * E[(EP_PART + 0) * fs] * safe + E[(EP_PART + 1) * fs] * safe + args.tail_rho) * safe;
    #undef DOT
    #undef EDOT

    // ---- tests
    uint outc = OUT_SPLIT, why = WHY_SPLIT;
    bool decided = false, tried_infl = false;
    if (ev.rho_up < args.eps) { outc = OUT_DISCARD; why = WHY_RHO; decided = true; }
    float val[3], vec[3][3];
    eigen3(ev.A, val, vec);
    if (!decided) {
        float3 dirs[7];
        dirs[0] = float3(1, 0, 0); dirs[1] = float3(0, 1, 0); dirs[2] = float3(0, 0, 1);
        int nd = 3;
        const float gn = sqrt(ev.g[0] * ev.g[0] + ev.g[1] * ev.g[1] + ev.g[2] * ev.g[2]);
        if (gn > 0.0) { dirs[nd] = float3(ev.g[0], ev.g[1], ev.g[2]) / gn; nd++; }
        [ForceUnroll] for (int e = 0; e < 3; ++e) { dirs[nd] = float3(vec[0][e], vec[1][e], vec[2][e]); nd++; }
        for (int t = 0; t < nd && !decided; ++t) if (excl_dir(ev, h, dirs[t])) { outc = OUT_DISCARD; why = WHY_EXCL; decided = true; }
    }
    const Lin lin = linear(ev, val, vec);
    if (!decided && lin.ok) {
        if (newton_dist(ev, lin, h)) { outc = OUT_DISCARD; why = WHY_NEWT; decided = true; }
        if (!decided) {
            const int kr = krawczyk(lin, ev.dH_up, h);
            if (kr == -1) { outc = OUT_DISCARD; why = WHY_KOUT; decided = true; }
            else if (kr == 1) { outc = OUT_ROOT; why = WHY_ROOT; decided = true; }
            else {
                tried_infl = true;
                const int ki = krawczyk(lin, ev.dH_i, hi);
                if (ki == -1) { outc = OUT_DISCARD; why = WHY_KOUT_I; decided = true; }
                else if (ki == 1) { outc = OUT_ROOT; why = WHY_ROOT; decided = true; }
            }
        }
    }
    if (!decided) {
        // fp32 noise: the gradient lies inside its own rounding bound and a Krawczyk proof would need a
        // larger cube than this one; refining in fp32 cannot help.
        float gmax = 0.0, egmax = 0.0, emmax = 0.0;
        [ForceUnroll] for (int c = 0; c < 3; ++c) { gmax = max(gmax, abs(ev.g[c])); egmax = max(egmax, ev.e_g[c]); }
        if (lin.ok) [ForceUnroll] for (int c = 0; c < 3; ++c) emmax = max(emmax, lin.em[c]);
        const bool g_in_noise = gmax <= 4.0 * egmax;
        const bool kraw_hopeless = !lin.ok || h <= 2.0 * emmax / (1.0 - lin.rhoR);
        if (g_in_noise && kraw_hopeless) { outc = OUT_ESCALATE; why = WHY_NOISE; }
        else if (h <= args.h_min) { outc = OUT_ESCALATE; why = WHY_HMIN; }
        else { outc = OUT_SPLIT; why = WHY_SPLIT; }
    }
    if (tried_infl) outc |= OUT_INFLATED;
    args.outcome[cube] = outc | (why << 8);
    if ((args.flags & 1u) != 0u) {
        float* dbg = args.debug + cube * 16u;
        dbg[0] = ev.rho; dbg[1] = ev.e_rho; dbg[2] = ev.rho_up;
        dbg[3] = ev.g[0]; dbg[4] = ev.g[1]; dbg[5] = ev.g[2];
        dbg[6] = ev.e_g[0]; dbg[7] = ev.e_g[1]; dbg[8] = ev.e_g[2];
        dbg[9] = ev.A[0][0]; dbg[10] = ev.A[1][1]; dbg[11] = ev.A[2][2];
        dbg[12] = ev.r_up[0]; dbg[13] = ev.dH_up[0][0]; dbg[14] = fn; dbg[15] = ev.q_up[0];
    }
}

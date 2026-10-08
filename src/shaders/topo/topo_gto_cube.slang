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
// so one sub-block of D. A chunk of batches runs as four stages, every intermediate in device
// memory (rows = local AOs of all batches in the chunk, columns = channel * 8 + child):
//
//   ao_main     one group per batch. Per local shell, primitive, axis and child side the 1D factors
//               of the centre values and of the 1D remainder sups (cube and inflated cube) are
//               computed once into groupshared memory; they serve every Cartesian component and the
//               4 children on that side. Then per (AO, child) the products: centre values with
//               magnitude and error shadows, and the separable remainder bounds.           -> phi
//   gemm_main   64-row tiles: D * V and |D| * (error weights, remainder vectors), over the batch's
//               local AOs.                                                                 -> dv
//   epi_main    32 threads per batch: DOT[a][b] = V_a . D V_b with its error sums (which depend on
//               a only through its derivative order), and the AO sums of the enclosures.   -> ep
//   decide_main one thread per cube: derivatives, enclosures, tests.                       -> outcome
//
// FACTORED FORM. With D = sum_k l_k c_k c_k^T (md_topo.c: cpg_factor_density) a batch whose rank r is
// well below its local AO count n (cpg_use_factors) takes fgemm + epi_fact instead of gemm_main +
// epi_main: fgemm transforms its phi rows to r factor rows, c_k . V and |c_k| . (errors, magnitudes,
// remainders), r n 480 multiply-adds instead of n^2 288; epi_fact runs the same sums over the r rows
// with X = diag(l) in place of D. Both kinds of batch share a chunk; each kernel skips the other kind.
//
// Every sum runs in a fixed order, so a cube's result depends only on its batch, never on the chunk.
//
// Outcome per cube (low 4 bits): 0 discard, 1 root (Krawczyk certified; the CPU polishes it in
// double), 2 split, 3 escalate (fp32 cannot decide: the gradient is inside its rounding noise and the
// cube is too small for a Krawczyk proof at that noise, or h <= h_min; the CPU takes the cube and its
// subtree in double). Bit 4: the inflated enclosure decided or was tried. Bits 8..15: which test.
// Bits 16..23, for a split: the children that still need an evaluation (child_excluded below).

#include "../md_gpu.slang"

#define LMAX    4
#define KA      7       // 1D derivative orders 0..6 (centre values to 5, sups 3..6)
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
#define WHY_CHILDREN 11u        // split, but every child is excluded by this cube's own expansion

struct GtoShell {
    float4 A_r2;                // centre (Bohr), squared screening radius (unused here)
    uint4  info;                // prim_offset, num_prims, ao_offset, l | (ncart << 8)
};

struct Batch {
    float4 P_h;                 // parent centre, child half-width h; children at P +- h (all exact)
    uint4  info;                // child mask (bit s: x side s&1, y side s&2, z side s&4), local AOs n,
                                // first row, first local shell in batch_shells
    uint4  info2;               // number of local shells, first factor row in fphi (FACT_NONE: the
                                // batch uses D as a matrix), -, -
};

#define FACT_NONE 0xFFFFFFFFu

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
    uint  fac_r;                // factors of D (0: none)
    uint  num_fbatches;         // batches in the factored form
    Batch*    batches;
    uint*     batch_shells;     // global shell index per local shell, ascending per batch
    uint*     row_ao;           // global AO index per row (written by ao_main)
    uint2*    row_tiles;        // (batch, first local row) per 32-row tile of the GEMM
    GtoShell* shells;
    float*    alpha;
    float*    coeff;
    float*    ao_nrm;
    uint*     ao_ijk;           // i | j << 8 | k << 16
    float*    D;                // num_ao * num_ao
    float*    phi;              // rows * PHI_W
    float*    dv;               // rows * DV_W
    float*    ep;               // EP_W * 8 * num_batches, by field
    uint*     outcome;          // 8 * num_batches
    float*    debug;            // optional, 16 floats per cube
    float*    fac_C;            // D = sum_k l_k c_k c_k^T: c_k[i] at [i * fac_r + k]
    float*    fac_l;            // l_k
    float*    fphi;             // factor rows of the factored batches, phi's layout
    uint*     fbatches;         // the factored batches' indices
};
MD_SHADER_ARGS(RootArgs);
#define args (*md_root.args)

bool batch_factored(Batch bt) { return bt.info2.y != FACT_NONE; }
// Rows of the batch's quadratic form: its local AOs, or the factors.
uint batch_rows(Batch bt) { return batch_factored(bt) ? args.fac_r : bt.info.y; }

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

static const uint MI[35] = {
    0x000000, 0x000001, 0x000100, 0x010000, 0x000002, 0x000101, 0x010001, 0x000200, 0x010100, 0x020000,
    0x000003, 0x000102, 0x010002, 0x000201, 0x010101, 0x020001, 0x000300, 0x010200, 0x020100, 0x030000,
    0x000004, 0x000103, 0x010003, 0x000202, 0x010102, 0x020002, 0x000301, 0x010201, 0x020101, 0x030001,
    0x000400, 0x010300, 0x020200, 0x030100, 0x040000 };
int3 mi_of(int m) { const uint v = MI[m]; return int3(int(v & 0xFFu), int((v >> 8) & 0xFFu), int((v >> 16) & 0xFFu)); }
int order_of(int m) { return m == 0 ? 0 : (m < 4 ? 1 : (m < 10 ? 2 : 3)); }
int ab_index(int a, int b) { const int lo = min(a, b), hi = max(a, b); return lo == 0 ? hi : (lo == 1 ? 2 + hi : 5); }

// kappa[i][a][m]: d^a/dx^a [x^i exp(-al x^2)] = sum over m of kappa al^((m - i + a) / 2) x^m exp(-al x^2)
// (cpg_init_tables in md_topo.c), i <= LMAX, a < KA. Integers; compile-time constants once the loops
// over i, a and m are unrolled, so the zero terms vanish.
static const float KAPPA[LMAX + 1][KA][LMAX + KA] = {
    { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { -2, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 12, 0, -8, 0, 0, 0, 0, 0, 0, 0 },
      { 12, 0, -48, 0, 16, 0, 0, 0, 0, 0, 0 },
      { 0, -120, 0, 160, 0, -32, 0, 0, 0, 0, 0 },
      { -120, 0, 720, 0, -480, 0, 64, 0, 0, 0, 0 } },
    { { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 1, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, -6, 0, 4, 0, 0, 0, 0, 0, 0, 0 },
      { -6, 0, 24, 0, -8, 0, 0, 0, 0, 0, 0 },
      { 0, 60, 0, -80, 0, 16, 0, 0, 0, 0, 0 },
      { 60, 0, -360, 0, 240, 0, -32, 0, 0, 0, 0 },
      { 0, -840, 0, 1680, 0, -672, 0, 64, 0, 0, 0 } },
    { { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0 },
      { 2, 0, -10, 0, 4, 0, 0, 0, 0, 0, 0 },
      { 0, -24, 0, 36, 0, -8, 0, 0, 0, 0, 0 },
      { -24, 0, 156, 0, -112, 0, 16, 0, 0, 0, 0 },
      { 0, 360, 0, -760, 0, 320, 0, -32, 0, 0, 0 },
      { 360, 0, -3000, 0, 3120, 0, -864, 0, 64, 0, 0 } },
    { { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 3, 0, -2, 0, 0, 0, 0, 0, 0 },
      { 0, 6, 0, -14, 0, 4, 0, 0, 0, 0, 0 },
      { 6, 0, -54, 0, 48, 0, -8, 0, 0, 0, 0 },
      { 0, -120, 0, 300, 0, -144, 0, 16, 0, 0, 0 },
      { -120, 0, 1140, 0, -1320, 0, 400, 0, -32, 0, 0 },
      { 0, 2520, 0, -7560, 0, 5040, 0, -1056, 0, 64, 0 } },
    { { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 4, 0, -2, 0, 0, 0, 0, 0 },
      { 0, 0, 12, 0, -18, 0, 4, 0, 0, 0, 0 },
      { 0, 24, 0, -96, 0, 60, 0, -8, 0, 0, 0 },
      { 24, 0, -336, 0, 492, 0, -176, 0, 16, 0, 0 },
      { 0, -720, 0, 2640, 0, -2040, 0, 480, 0, -32, 0 },
      { -720, 0, 9360, 0, -15480, 0, 7440, 0, -1248, 0, 64 } } };

// al^0 .. al^(KA - 1), each by one more multiplication (the roundings of powi)
void al_pows(float al, out float alp[KA]) { alp[0] = 1.0; [ForceUnroll] for (int k = 1; k < KA; ++k) alp[k] = alp[k - 1] * al; }

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

// An upper bound of exp(-t) that never underflows: exp_neg returns 0 beyond t = 87, and a GPU may
// flush what would be denormal. 2^-100 > exp(-87), and products of it with any factor above 2^-26 stay
// normal (the rest is below the absolute floors added downstream).
#define EXP_FLOOR 7.888609052210118e-31
float exp_neg_up(float t) { return t > 87.0 ? EXP_FLOOR : exp_neg(t); }

uint3 ao_ijk(uint ao) { const uint v = args.ao_ijk[min(ao, args.num_ao - 1u)]; return uint3(v & 0xFFu, (v >> 8) & 0xFFu, (v >> 16) & 0xFFu); }

// ================================================================================== ao_main
// One group per batch. Local shells are taken in chunks whose 1D tables fit TAB_CAP floats (the host
// guarantees that every single shell fits).
//
// Table entry per (shell, primitive p, axis, side), te(l) = 20 (l + 1) + 2 floats, l1 = l + 1:
//   [a l1 + i]               a = 0..5: d^a/dx^a [x^i exp(-al x^2)] at the child centre     (i = 0..l)
//   [6 l1 + a l1 + i]        the same on |terms| (magnitude shadow)
//   [12 l1]                  rounding term of this axis for a <= 3: (2l + 20 + 8.1 t) u
//   [12 l1 + 1]              the same for a = 4, 5: (2l + 30 + 8.1 t) u
//                            (both 1 where t > 87: the value is then 0 and its error its whole magnitude,
//                            the magnitude taken with exp_neg_up)
//   [12 l1 + 2 + 4 i + a-3]  a = 3..6: sup over the child interval of |d^a/dx^a [x^i exp(-al x^2)]|
//   [16 l1 + 2 + 4 i + a-3]  the same over the inflated child interval
// at offset ((p * 3 + axis) * 2 + side) * te(l) in the shell's block.

#define AO_WG   128
#define TAB_CAP 4096
#define CS_MAX  32      // shells per chunk; TAB_CAP / 132 (one s primitive) < 32, so never the binding limit

groupshared float gs_tab[TAB_CAP];
groupshared uint  gs_cs_shell[CS_MAX];     // global shell
groupshared uint  gs_cs_pk[CS_MAX];        // l | nc << 8 | np << 16 of the candidates
groupshared uint  gs_cs_tab[CS_MAX + 1];   // table offset; [ncs] the total
groupshared uint  gs_cs_item[CS_MAX + 1];  // first table item (3 per entry: centre, sups over the cube and the inflated cube)
groupshared uint  gs_cs_prod[CS_MAX + 1];  // first product
groupshared uint  gs_cs_row[CS_MAX + 1];   // first row, relative to the batch's
groupshared uint  gs_ncs;

uint te_size(uint l) { return 20u * (l + 1u) + 2u; }

// sups of |d^a/dx^a [x^i exp(-al x^2)]| over [x0, x1], a = 3..6, for i = 0..l, into gs_tab[base + 4 i + a - 3]:
// the monomial terms of the expansion, each at its own maximum on the interval.
// The table of x^m exp(-al x^2) at its maximiser over the interval does not depend on i.
void sup_tab(uint base, uint l, float x0, float x1, float al, float alp[KA]) {
    x0 -= 4.0 * UF * abs(x0) + 1e-30;    // the computed interval must contain the exact one
    x1 += 4.0 * UF * abs(x1) + 1e-30;
    const float ax0 = abs(x0), ax1 = abs(x1);
    const float tlo = (x0 <= 0.0 && x1 >= 0.0) ? 0.0 : min(ax0, ax1);
    const float thi = max(ax0, ax1);
    float Sm[LMAX + KA], eS[LMAX + KA];
    [ForceUnroll] for (int m = 0; m < LMAX + KA; ++m) {
        Sm[m] = 0.0; eS[m] = 0.0;
        if (m <= int(l) + KA - 1) {
            float ts = sqrt(float(m) / (2.0 * al));
            ts = clamp(ts, tlo, thi);
            const float t = al * ts * ts;
            Sm[m] = powi(ts, m) * exp_neg_up(t);
            eS[m] = (float(2 * m + 4) + 6.0 + 8.1 * t) * UF;
        }
    }
    [ForceUnroll] for (int i = 0; i <= LMAX; ++i) {
        if (i <= int(l)) {
            [ForceUnroll] for (int a = 3; a < KA; ++a) {
                float acc = 0.0;
                [ForceUnroll] for (int m = 0; m < LMAX + KA; ++m) {
                    const float kk = m <= i + a ? KAPPA[i][a][m] : 0.0;
                    if (kk != 0.0) {
                        const int kp = max(m - i + a, 0) / 2;
                        acc += abs(kk) * alp[kp] * Sm[m] * (1.0 + eS[m] + float(kp + 2) * UF);
                    }
                }
                gs_tab[base + uint(i * 4 + a - 3)] = acc * (1.0 + float(i + a + 3) * UF);
            }
        }
    }
}

// The centre part of one table entry: the factors at x (child centre - atom) and the rounding terms.
void tab_centre(uint base, uint l, float x, float al) {
    const uint l1 = l + 1u;
    const float t = al * x * x;
    const float E = exp_neg(t), Eup = exp_neg_up(t);
    float xp[LMAX + 6], axp[LMAX + 6];
    xp[0] = 1.0; axp[0] = 1.0;
    [ForceUnroll] for (int m = 1; m < LMAX + 6; ++m) { xp[m] = xp[m - 1] * x; axp[m] = axp[m - 1] * abs(x); }
    float alp[KA];
    al_pows(al, alp);
    [ForceUnroll] for (int i = 0; i <= LMAX; ++i) {
        if (i <= int(l)) {
            [ForceUnroll] for (int a = 0; a < 6; ++a) {
                float v = 0.0, mg = 0.0;
                [ForceUnroll] for (int m = 0; m < LMAX + 6; ++m) {
                    const float k = m <= i + a ? KAPPA[i][a][m] : 0.0;
                    if (k != 0.0) {
                        const float ka = k * alp[max(m - i + a, 0) / 2];
                        v += ka * xp[m];
                        mg += abs(ka) * axp[m];
                    }
                }
                gs_tab[base + uint(a) * l1 + uint(i)] = v * E;
                gs_tab[base + 6u * l1 + uint(a) * l1 + uint(i)] = mg * Eup;
            }
        }
    }
    // beyond exp_neg's range the value is set to 0: its error is then its whole magnitude
    gs_tab[base + 12u * l1] = t > 87.0 ? 1.0 : (float(2u * l + 14u) + 6.0 + 8.1 * t) * UF;
    gs_tab[base + 12u * l1 + 1u] = t > 87.0 ? 1.0 : (float(2u * l + 24u) + 6.0 + 8.1 * t) * UF;
}

// Separable remainder bounds (cpg_shell_rem in md_topo.c has the derivation). For one 1D factor
// d^be/dx^be [x^i exp(-al x^2)] of a primitive, the absolute degree-graded expansion over the cube:
// u[m] = |value of order be + m| h^m / m! for m <= P (rounding included through the magnitude
// shadow), u[P + 1] = sup of order be + P + 1 times h^(P+1) / (P+1)!. Returns its first terms u[0..2]
// and suffix sums S[r] = sum over m >= r of u[m], r = 0..3, all a slot needs.
static const float INV_FACT[7] = { 1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0 };
void axis_terms(uint tb, uint l1, uint i, int be, int P, float hp[7], uint soff, out float u3[3], out float S4[4]) {
    const float e3 = gs_tab[tb + 12u * l1], e5 = gs_tab[tb + 12u * l1 + 1u];
    float uu[6];
    [ForceUnroll] for (int m = 0; m < 6; ++m) {
        uu[m] = 0.0;
        if (m <= P) {
            const int a = be + m;
            const float v = gs_tab[tb + uint(a) * l1 + i], mg = gs_tab[tb + uint(6 + a) * l1 + i];
            uu[m] = (abs(v) + (a <= 3 ? e3 : e5) * mg) * hp[m] * INV_FACT[m];
        } else if (m == P + 1) {
            uu[m] = gs_tab[tb + soff + 4u * i + uint(be + P + 1 - 3)] * hp[m] * INV_FACT[m];
        }
    }
    float S[7];
    S[6] = 0.0;
    [ForceUnroll] for (int m = 5; m >= 0; --m) S[m] = S[m + 1] + uu[m];
    [ForceUnroll] for (int m = 0; m < 3; ++m) u3[m] = uu[m];
    [ForceUnroll] for (int r = 0; r < 4; ++r) S4[r] = S[r];
}

// Sum over m0 + m1 + m2 >= K of u0[m0] u1[m1] u2[m2] from the first terms and suffix sums (K <= 3).
float sep_bound(int K, float u0[3], float S0[4], float u1[3], float S1[4], float S2[4]) {
    float T12[4];
    T12[0] = S1[0] * S2[0];
    [ForceUnroll] for (int j = 1; j < 4; ++j) {
        float v = S1[j] * S2[0];
        [ForceUnroll] for (int m1 = 0; m1 < 3; ++m1) if (m1 < j) v += u1[m1] * S2[j - m1];
        T12[j] = v;
    }
    float b = S0[K] * T12[0];
    [ForceUnroll] for (int m0 = 0; m0 < 3; ++m0) if (m0 < K) b += u0[m0] * T12[K - m0];
    return b;
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
            gx[a] = gs_tab[bx + uint(a) * l1 + uint(e.x)]; Gx[a] = gs_tab[bx + (6u + uint(a)) * l1 + uint(e.x)];
            gy[a] = gs_tab[by + uint(a) * l1 + uint(e.y)]; Gy[a] = gs_tab[by + (6u + uint(a)) * l1 + uint(e.y)];
            gz[a] = gs_tab[bz + uint(a) * l1 + uint(e.z)]; Gz[a] = gs_tab[bz + (6u + uint(a)) * l1 + uint(e.z)];
        }
        float eps = 5.0 * UF;
        eps += gs_tab[bx + 12u * l1];
        eps += gs_tab[by + 12u * l1];
        eps += gs_tab[bz + 12u * l1];
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

    // ---- remainder vectors, cube and inflated cube: separable bounds (cpg_shell_rem in md_topo.c).
    //      E[0] E1, [1] E2, [2] E3, [3..5] Ea2, [6..8] Ea3, [9..14] Eab2 (a <= b); the inflated cube
    //      needs only the K = 2 slots (E2, Ea2, Eab2). Every term is a product of non-negative upper
    //      bounds and every sum adds non-negative terms, so the rounding of the whole chain (at most
    //      ~40 operations, plus the sum over primitives) is covered by one relative factor.
    const float fE = 1.0 + float(2u * np + 96u) * UF;
    [ForceUnroll] for (int pass = 0; pass < 2; ++pass) {
        const float hh = pass == 0 ? h : 1.5 * h;
        float hp[7];
        hp[0] = 1.0;
        [ForceUnroll] for (int m = 1; m < 7; ++m) hp[m] = hp[m - 1] * hh;
        const uint soff = 12u * l1 + 2u + uint(pass) * 4u * l1;
        float E[15];
        [ForceUnroll] for (int t = 0; t < 15; ++t) E[t] = 0.0;
        for (uint p = 0; p < np; ++p) {
            const float w = abs(args.coeff[sh.info.x + p]) * nrm;
            uint tb[3];
            tb[0] = tab0 + ((p * 3u + 0u) * 2u + sx) * te;
            tb[1] = tab0 + ((p * 3u + 1u) * 2u + sy) * te;
            tb[2] = tab0 + ((p * 3u + 2u) * 2u + sz) * te;
            const uint ei[3] = { uint(e.x), uint(e.y), uint(e.z) };
            // K = 2: E2, Ea2, Eab2 from the factors (be, P) = (0,3) (1,3) (2,3) per axis
            {
                float u[3][3][3], S[3][3][4];
                [ForceUnroll] for (int k = 0; k < 3; ++k) [ForceUnroll] for (int be = 0; be < 3; ++be) axis_terms(tb[k], l1, ei[k], be, 3, hp, soff, u[k][be], S[k][be]);
                E[1] += w * sep_bound(2, u[0][0], S[0][0], u[1][0], S[1][0], S[2][0]);
                E[3] += w * sep_bound(2, u[0][1], S[0][1], u[1][0], S[1][0], S[2][0]);
                E[4] += w * sep_bound(2, u[0][0], S[0][0], u[1][1], S[1][1], S[2][0]);
                E[5] += w * sep_bound(2, u[0][0], S[0][0], u[1][0], S[1][0], S[2][1]);
                E[9]  += w * sep_bound(2, u[0][2], S[0][2], u[1][0], S[1][0], S[2][0]);   // xx
                E[10] += w * sep_bound(2, u[0][1], S[0][1], u[1][1], S[1][1], S[2][0]);   // xy
                E[11] += w * sep_bound(2, u[0][1], S[0][1], u[1][0], S[1][0], S[2][1]);   // xz
                E[12] += w * sep_bound(2, u[0][0], S[0][0], u[1][2], S[1][2], S[2][0]);   // yy
                E[13] += w * sep_bound(2, u[0][0], S[0][0], u[1][1], S[1][1], S[2][1]);   // yz
                E[14] += w * sep_bound(2, u[0][0], S[0][0], u[1][0], S[1][0], S[2][2]);   // zz
            }
            if (pass == 0) {
                // K = 1: E1 from (0,2); K = 3: E3, Ea3 from (0,4) (1,4)
                float u[3][3], S[3][4];
                [ForceUnroll] for (int k = 0; k < 3; ++k) axis_terms(tb[k], l1, ei[k], 0, 2, hp, soff, u[k], S[k]);
                E[0] += w * sep_bound(1, u[0], S[0], u[1], S[1], S[2]);
                float v[3][2][3], T[3][2][4];
                [ForceUnroll] for (int k = 0; k < 3; ++k) [ForceUnroll] for (int be = 0; be < 2; ++be) axis_terms(tb[k], l1, ei[k], be, 4, hp, soff, v[k][be], T[k][be]);
                E[2] += w * sep_bound(3, v[0][0], T[0][0], v[1][0], T[1][0], T[2][0]);
                E[6] += w * sep_bound(3, v[0][1], T[0][1], v[1][0], T[1][0], T[2][0]);
                E[7] += w * sep_bound(3, v[0][0], T[0][0], v[1][1], T[1][1], T[2][0]);
                E[8] += w * sep_bound(3, v[0][0], T[0][0], v[1][0], T[1][0], T[2][1]);
            }
        }
        [ForceUnroll] for (int t = 0; t < 15; ++t) E[t] = E[t] * fE + float(np) * 1e-36;   // + what a flush may have lost
        if (pass == 0) {
            [ForceUnroll] for (int t = 0; t < 6; ++t) out[(PC_E + t) * 8] = E[t];
            [ForceUnroll] for (int t = 0; t < 9; ++t) out[(PC_E2 + t) * 8] = E[6 + t];
        } else {
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

    uint next = 0, rowc = 0;                // the next shell and its row (uniform)
    for (;;) {
        // ---- chunk layout: the next CS_MAX shells are read in parallel, then thread 0 takes as many as
        //      fit TAB_CAP. Written only here, read by all after the barrier below and not again after
        //      the barrier that ends the iteration.
        if (tid < CS_MAX && next + tid < nsh) {
            const uint g = min(args.batch_shells[sh0 + next + tid], args.num_shells - 1u);
            const uint4 inf = args.shells[g].info;
            gs_cs_shell[tid] = g;
            gs_cs_pk[tid] = (inf.w & 0xFFFFu) | (min(inf.y, 0xFFFFu) << 16);
        }
        GroupMemoryBarrierWithGroupSync();
        if (tid == 0) {
            const uint cand = min(nsh - next, uint(CS_MAX));
            uint c = 0, tab = 0, item = 0, prod = 0, row = rowc;
            for (; c < cand; ++c) {
                const uint pk = gs_cs_pk[c];
                const uint l = pk & 0xFFu, nc = (pk >> 8) & 0xFFu, np = pk >> 16;
                const uint sz = np * 6u * te_size(l);
                if (tab + sz > TAB_CAP) break;
                gs_cs_tab[c] = tab; gs_cs_item[c] = item; gs_cs_prod[c] = prod; gs_cs_row[c] = row;
                tab += sz; item += np * 18u; prod += nc * 8u; row += nc;
            }
            gs_cs_tab[c] = tab; gs_cs_item[c] = item; gs_cs_prod[c] = prod; gs_cs_row[c] = row;
            gs_ncs = c;                     // 0: done (or a shell that does not fit, which the host never sends)
        }
        GroupMemoryBarrierWithGroupSync();
        const uint ncs = gs_ncs;
        if (ncs == 0) break;                // uniform
        next += ncs;
        rowc = gs_cs_row[ncs];

        // ---- 1D tables: (shell, part, primitive, axis, side), part 0 the centre factors, 1 and 2 the sups
        //      over the cube and the inflated cube (a thread's item index only grows, and with it its
        //      shell: the search resumes where it stopped)
        const uint nitem = gs_cs_item[ncs];
        uint cs = 0;
        for (uint f = tid; f < nitem; f += AO_WG) {
            while (cs + 1u < ncs && gs_cs_item[cs + 1u] <= f) cs++;
            const uint pk = gs_cs_pk[cs];
            const uint l = pk & 0xFFu, np = pk >> 16, l1 = l + 1u;
            const uint r = f - gs_cs_item[cs];
            const uint part = r / (np * 6u), e = r % (np * 6u);
            const uint p = e / 6u, ax = (e % 6u) / 2u, side = e % 2u;
            const GtoShell sh = args.shells[min(gs_cs_shell[cs], args.num_shells - 1u)];
            const float al = args.alpha[sh.info.x + min(p, sh.info.y - 1u)];
            const float A = sh.A_r2[ax];
            const float c = P[ax] + (side != 0u ? h : -h);
            const uint base = gs_cs_tab[cs] + e * te_size(l);
            if (part == 0u) {
                tab_centre(base, l, c - A, al);
            } else {
                float alp[KA];
                al_pows(al, alp);
                const float hs = part == 1u ? h : 1.5 * h;
                sup_tab(base + (part == 1u ? 12u : 16u) * l1 + 2u, l, (c - hs) - A, (c + hs) - A, al, alp);
            }
        }
        GroupMemoryBarrierWithGroupSync();

        // ---- products: (AO, child), children fastest
        const uint nprod = gs_cs_prod[ncs];
        cs = 0;
        for (uint f = tid; f < nprod; f += AO_WG) {
            while (cs + 1u < ncs && gs_cs_prod[cs + 1u] <= f) cs++;
            const uint r = f - gs_cs_prod[cs];
            const uint ci = r / 8u, k = r % 8u;
            const uint row = row0 + gs_cs_row[cs] + ci;
            if (k == 0u) args.row_ao[row] = args.shells[min(gs_cs_shell[cs], args.num_shells - 1u)].info.z + ci;   // for the GEMMs
            ao_products(cs, ci, k, row, h, n);
        }
        GroupMemoryBarrierWithGroupSync();
    }
}

// ================================================================================== gemm_main
// dv[row][c] = sum_j X[row][j] phi[j][c] over the batch's local AOs j, X = D for the 160 columns of the
// centre values and |D| for the 112 after them. One group per 64-row tile computes all 272 columns, so
// each D tile is gathered once: 32 x 8 threads, thread (tx, ty) owns rows 8 ty .. 8 ty + 7 and columns
// tx + 32 j, j = 0..8 (j < 5 with D, j >= 5 with |D|): 17 groupshared reads per 72 multiply-adds. The
// sum over j runs in increasing j, so the result does not depend on the tiling.
// Of seven tilings measured (64/32/16 rows, 128/256 threads, 288 or 96 columns per group, 24-72
// accumulators), this one was within 11% of the fastest on an M4 Pro, an RTX 3080 and a UHD 770 alike.

#define GT  64                  // rows per tile
#define GKC 16                  // j per chunk
#define GCJ 9                   // column groups of 32: 9 * 32 = 288 >= DV_W
#define GR  8                   // rows per thread
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

    float acc[GR][GCJ];
    [ForceUnroll] for (int i = 0; i < GR; ++i) [ForceUnroll] for (int j = 0; j < GCJ; ++j) acc[i][j] = 0.0;

    for (uint k0 = 0; k0 < n; k0 += GKC) {
        // D tile: 32 rows x 16 j, 4 per thread, j fastest (consecutive AOs are mostly consecutive in D)
        [ForceUnroll] for (int q = 0; q < (GT * GKC) / GEMM_WG; ++q) {
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
        // phi tile: 16 j x 288 columns, 36 per thread, columns fastest (coalesced rows)
        [ForceUnroll] for (int q = 0; q < (GKC * GCJ * 32) / GEMM_WG; ++q) {
            const uint e = tid + uint(GEMM_WG * q);
            const uint kk = e / (GCJ * 32), c = e % (GCJ * 32);
            float v = 0.0;
            if (k0 + kk < n && c < DV_W) v = args.phi[(base + k0 + kk) * PHI_W + c];
            gs_B[kk][c] = v;
        }
        GroupMemoryBarrierWithGroupSync();
        [ForceUnroll] for (int kk = 0; kk < GKC; ++kk) {
            float a[GR], aa[GR], bb[GCJ];
            [ForceUnroll] for (int i = 0; i < GR; ++i) { a[i] = gs_A[kk][ty * GR + uint(i)]; aa[i] = abs(a[i]); }
            [ForceUnroll] for (int j = 0; j < GCJ; ++j) bb[j] = gs_B[kk][tx + 32u * uint(j)];
            [ForceUnroll] for (int i = 0; i < GR; ++i) [ForceUnroll] for (int j = 0; j < GCJ; ++j) acc[i][j] += (j < 5 ? a[i] : aa[i]) * bb[j];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    [ForceUnroll] for (int i = 0; i < GR; ++i) {
        const uint row = r0 + ty * GR + uint(i);
        if (row >= n) continue;
        [ForceUnroll] for (int j = 0; j < GCJ; ++j) {
            const uint col = tx + 32u * uint(j);
            if (col < DV_W) args.dv[(base + row) * DV_W + col] = acc[i][j];
        }
    }
}

// For the column-split tilings (FGEMM): bit j set: column group c * GCJ_ + j of block c is one of the
// |X| groups (>= 5).
#define GEMM_ABS_MASK(c, GCJ_) (0xFFFFu << uint(max(5 - (c) * (GCJ_), 0)))
// One GKC step of products over the groupshared tiles, |D| for the columns in MASK.
#define GEMM_TILE_STEP(NAME, GT_, TY_, GCJ_, MASK)                                                     \
        {                                                                                               \
        [ForceUnroll] for (int kk = 0; kk < GKC; ++kk) {                                                \
            float a[(GT_) / (TY_)], aa[(GT_) / (TY_)], bb[GCJ_];                                        \
            [ForceUnroll] for (int i = 0; i < (GT_) / (TY_); ++i) { a[i] = NAME##_A[kk][ty * uint((GT_) / (TY_)) + uint(i)]; aa[i] = abs(a[i]); } \
            [ForceUnroll] for (int j = 0; j < (GCJ_); ++j) bb[j] = NAME##_B[kk][tx + 32u * uint(j)];    \
            [ForceUnroll] for (int i = 0; i < (GT_) / (TY_); ++i) [ForceUnroll] for (int j = 0; j < (GCJ_); ++j) \
                acc[i][j] += ((((MASK) >> uint(j)) & 1u) != 0u ? aa[i] : a[i]) * bb[j];                 \
        }                                                                                               \
        }

// ---- Factored form (md_topo.c: cpg_factor_density). Per factored batch, its factor rows
//   fphi[k][c] = sum_j X_c[k][j] phi[j][c],   X = C~^T for the 20 centre values, |C~|^T for the rest,
// over the batch's local AOs j, C~ the factors in float. Every phi channel is transformed (the 4 value
// errors go unused), so a factor row has phi's layout: centre values c~_k . V, and |c~_k| . (WV, MV, E)
// as upper bounds of its value errors, magnitudes and remainder terms (epi_main). One group per
// (factored batch, GT_-row tile of factors, column block): 5 blocks of 3 x 32 columns (480 >= PHI_W),
// 32 x TY_ threads, GT_ / TY_ rows and 3 columns each; the sum over j in increasing order.
#define FGEMM(NAME, GT_, TY_)                                                                            \
groupshared float NAME##_A[GKC][GT_];                                                                   \
groupshared float NAME##_B[GKC][96];                                                                    \
[shader("compute")]                                                                                     \
[numthreads(32, TY_, 1)]                                                                                \
void NAME(uint3 gid : SV_GroupID, uint3 gtid : SV_GroupThreadID, uint tid : SV_GroupIndex) {            \
    const uint R = args.fac_r;                                                                          \
    const uint tpb = max((R + uint(GT_) - 1u) / uint(GT_), 1u);                                         \
    const uint cb = gid.x % 5u, t = gid.x / 5u;                                                         \
    const uint fb = t / tpb;                                                                            \
    if (fb >= args.num_fbatches) return;                                                                \
    const uint r0 = (t % tpb) * uint(GT_);                                                              \
    const Batch bt = args.batches[min(args.fbatches[fb], args.num_batches - 1u)];                       \
    const uint n = bt.info.y, base = bt.info.z, fbase = bt.info2.y;                                     \
    const uint tx = gtid.x, ty = gtid.y;                                                                \
    const uint N = args.num_ao;                                                                         \
    const uint c0 = cb * 96u;                                                                           \
    float acc[(GT_) / (TY_)][3];                                                                        \
    [ForceUnroll] for (int i = 0; i < (GT_) / (TY_); ++i) [ForceUnroll] for (int j = 0; j < 3; ++j) acc[i][j] = 0.0; \
    for (uint k0 = 0; k0 < n; k0 += GKC) {                                                              \
        /* factor tile: GT_ factors x 16 AOs, factors fastest (consecutive in fac_C) */                 \
        [ForceUnroll] for (int q = 0; q < ((GT_) * GKC) / (32 * (TY_)); ++q) {                          \
            const uint e = tid + uint(32 * (TY_) * q);                                                  \
            const uint r = e % uint(GT_), kk = e / uint(GT_);                                           \
            float c = 0.0;                                                                              \
            if (r0 + r < R && k0 + kk < n) {                                                            \
                const uint aj = min(args.row_ao[base + k0 + kk], N - 1u);                               \
                c = args.fac_C[aj * R + r0 + r];                                                        \
            }                                                                                           \
            NAME##_A[kk][r] = c;                                                                        \
        }                                                                                               \
        [ForceUnroll] for (int q = 0; q < (GKC * 96) / (32 * (TY_)); ++q) {                             \
            const uint e = tid + uint(32 * (TY_) * q);                                                  \
            const uint kk = e / 96u, c = e % 96u;                                                       \
            float v = 0.0;                                                                              \
            if (k0 + kk < n && c0 + c < PHI_W) v = args.phi[(base + k0 + kk) * PHI_W + c0 + c];         \
            NAME##_B[kk][c] = v;                                                                        \
        }                                                                                               \
        GroupMemoryBarrierWithGroupSync();                                                              \
        if (cb == 0u)      GEMM_TILE_STEP(NAME, GT_, TY_, 3, GEMM_ABS_MASK(0, 3))                    \
        else if (cb == 1u) GEMM_TILE_STEP(NAME, GT_, TY_, 3, GEMM_ABS_MASK(1, 3))                    \
        else               GEMM_TILE_STEP(NAME, GT_, TY_, 3, GEMM_ABS_MASK(2, 3))                    \
        GroupMemoryBarrierWithGroupSync();                                                              \
    }                                                                                                   \
    [ForceUnroll] for (int i = 0; i < (GT_) / (TY_); ++i) {                                             \
        const uint row = r0 + ty * uint((GT_) / (TY_)) + uint(i);                                       \
        if (row >= R) continue;                                                                         \
        [ForceUnroll] for (int j = 0; j < 3; ++j) {                                                     \
            const uint col = c0 + tx + 32u * uint(j);                                                   \
            if (col < PHI_W) args.fphi[(fbase + row) * PHI_W + col] = acc[i][j];                        \
        }                                                                                               \
    }                                                                                                   \
}

//    name      GT  TY      factors per tile; threads; 8 x 3 outputs per thread
FGEMM(fgemm16,  16,  2)  // r <= 16;  64
FGEMM(fgemm32,  32,  4)  // r <= 32; 128
FGEMM(fgemm64,  64,  8)  // larger;  256

// ================================================================================== epi_main
// Per child cube, sums over the batch's local AOs (rows) i:
//   DOT[a][b] = sum_i V_i[a] (D V)_i[b], a = 0..19, b = 0..9, and its error sum, which depends on a only
//   through its derivative order o: e[o][b] = sum_i EV_i[o] (|DV_i[b]| + EDV_i[ob]) + MV_i[o] EDV_i[ob]
//                                           + (n + 1) u MV_i[o] |DV_i[b]|;
//   the enclosure's AO sums for the cube and the inflated cube.
// 32 threads per batch, two batches per group; thread (g, k), g = 0..3, k = child: DOT rows a = 5g .. 5g + 4
// and the error sums of order g, every row in increasing order; then the AO sums of rows i = g and g + 4
// (mod 8), each its own lane, the 8 lanes summed in lane order at the end. Rows are read straight from
// phi and dv (consecutive children are consecutive floats); no staging and one barrier.

#define EPI_WG 64

groupshared float gs_part[EPI_WG * 2][17];

// One row of the quadratic form for child k, as the enclosure sums use it: ad = |X V| plus its error
// bound, per channel; q = |X| E[0..5]; qs = |X| (inflated E2, Ea2) in the cube's slots 1, 3..5; the row's
// own remainder terms e (cube) and ei (inflated cube: slots 1, 3..5, 9..14).
struct RowIn {
    float ad[20];
    float q[6], qs[6];
    float e[15], ei[15];
};

// A local AO, X = D: dv and the |D| sums from gemm_main.
RowIn row_matrix(uint row, uint k, float safe, float fn) {
    RowIn ri;
    float* P = args.phi + row * PHI_W + k;
    float* Q = args.dv + row * DV_W + k;
    float edv[4];
    [ForceUnroll] for (int o = 0; o < 4; ++o) edv[o] = Q[(20 + o) * 8] * safe + fn * 1e-36;
    [ForceUnroll] for (int m = 0; m < 20; ++m) ri.ad[m] = abs(Q[m * 8]) + edv[order_of(m)];
    [ForceUnroll] for (int t = 0; t < 6; ++t) { ri.q[t] = Q[(24 + t) * 8] * safe; ri.qs[t] = 0.0; }
    ri.qs[1] = Q[30 * 8] * safe;
    [ForceUnroll] for (int a = 0; a < 3; ++a) ri.qs[3 + a] = Q[(31 + a) * 8] * safe;
    [ForceUnroll] for (int t = 0; t < 6; ++t) ri.e[t] = P[(PC_E + t) * 8];
    [ForceUnroll] for (int t = 0; t < 9; ++t) ri.e[6 + t] = P[(PC_E2 + t) * 8];
    [ForceUnroll] for (int t = 0; t < 15; ++t) ri.ei[t] = 0.0;
    ri.ei[1] = P[PC_EI * 8];
    [ForceUnroll] for (int a = 0; a < 3; ++a) ri.ei[3 + a] = P[(PC_EI + 1 + a) * 8];
    [ForceUnroll] for (int t = 0; t < 6; ++t) ri.ei[9 + t] = P[(PC_EI2 + t) * 8];
    return ri;
}

// Factor row P (fgemm), X = diag(l~), per derivative order o:
//   MV' = |C~| MV, rounded up;                                                magnitudes
//   EV' = |C~| WV + 4u MV'                                                    value errors
//   EDV' = |l~| (EV' + 2u MV')                                                errors of dv = l~ V'
// |C~| WV bounds the rounding of c~ . V and of V itself (WV = EV + (n + 2) u MV, n the GEMM length).
// 4u MV' covers D = C~ diag(l~) C~^T + Delta with |Delta| <= 4u sum_k |l~_k c~_k c~_k^T| (the factors in
// float; the double factorization's residual is far below that): Delta's polynomial part through the
// DOT error sums (MV' EDV' >= 4u |l~| MV'^2), its remainder part through 'safe'. The GEMM sums have
// non-negative terms, so 'safe' (>= 1 + (n + 64) u) rounds them up.
void fac_row_err(float* P, float alam, float safe, float fn, out float mvp[4], out float evp[4], out float edvp[4]) {
    [ForceUnroll] for (int o = 0; o < 4; ++o) {
        mvp[o] = P[(PC_MV + o) * 8] * safe;
        evp[o] = P[(PC_WV + o) * 8] * safe + 4.0 * UF * mvp[o];
        edvp[o] = alam * (evp[o] + 2.0 * UF * mvp[o]) * safe + fn * 1e-36;
    }
}

RowIn row_factor(uint frow, uint f, uint k, float safe, float fn) {
    RowIn ri;
    float* P = args.fphi + frow * PHI_W + k;
    const float lam = args.fac_l[min(f, args.fac_r - 1u)], alam = abs(lam);
    float mvp[4], evp[4], edvp[4];
    fac_row_err(P, alam, safe, fn, mvp, evp, edvp);
    [ForceUnroll] for (int m = 0; m < 20; ++m) ri.ad[m] = abs(lam * P[(PC_V + m) * 8]) + edvp[order_of(m)];
    [ForceUnroll] for (int t = 0; t < 6; ++t) ri.e[t] = P[(PC_E + t) * 8] * safe;
    [ForceUnroll] for (int t = 0; t < 9; ++t) ri.e[6 + t] = P[(PC_E2 + t) * 8] * safe;
    [ForceUnroll] for (int t = 0; t < 15; ++t) ri.ei[t] = 0.0;
    ri.ei[1] = P[PC_EI * 8] * safe;
    [ForceUnroll] for (int a = 0; a < 3; ++a) ri.ei[3 + a] = P[(PC_EI + 1 + a) * 8] * safe;
    [ForceUnroll] for (int t = 0; t < 6; ++t) ri.ei[9 + t] = P[(PC_EI2 + t) * 8] * safe;
    const float aq = alam * (1.0 + 4.0 * UF);
    [ForceUnroll] for (int t = 0; t < 6; ++t) { ri.q[t] = aq * ri.e[t]; ri.qs[t] = 0.0; }
    ri.qs[1] = aq * ri.ei[1];
    [ForceUnroll] for (int a = 0; a < 3; ++a) ri.qs[3 + a] = aq * ri.ei[3 + a];
    return ri;
}

// The enclosure sums of one row added into part[17]: cube [0..10], inflated cube Hessian [11..16].
void epi_row_parts(RowIn ri, float h, float hi, float h2, inout float part[17]) {
    float sum1 = 0.0, sum2 = 0.0;
    [ForceUnroll] for (int p = 0; p < 3; ++p) { sum1 += ri.ad[M1(p)]; [ForceUnroll] for (int r = 0; r < 3; ++r) sum2 += ri.ad[M2(p, r)]; }
    const float W = ri.ad[0] + h * sum1 + 0.5 * h2 * sum2;
    const float W1 = ri.ad[0] + h * sum1;
    const float W1i = ri.ad[0] + hi * sum1;
    float Wa[3], W1a[3], W1ai[3], W1ab[6], W1abi[6];
    [ForceUnroll] for (int a = 0; a < 3; ++a) {
        float t1 = 0.0, t2 = 0.0;
        [ForceUnroll] for (int p = 0; p < 3; ++p) { t1 += ri.ad[M2(a, p)]; [ForceUnroll] for (int r = 0; r < 3; ++r) t2 += ri.ad[M3(a, p, r)]; }
        Wa[a] = ri.ad[M1(a)] + h * t1 + 0.5 * h2 * t2;
        W1a[a] = ri.ad[M1(a)] + h * t1;
        W1ai[a] = ri.ad[M1(a)] + hi * t1;
    }
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int bq = 0; bq < 3; ++bq) {
        if (bq >= a) {
            float t1 = 0.0;
            [ForceUnroll] for (int p = 0; p < 3; ++p) t1 += ri.ad[M3(a, bq, p)];
            W1ab[ab_index(a, bq)] = ri.ad[M2(a, bq)] + h * t1;
            W1abi[ab_index(a, bq)] = ri.ad[M2(a, bq)] + hi * t1;
        }
    }
    part[0] += ri.ad[0] * ri.e[0];
    part[1] += ri.e[0] * ri.q[0];
    [ForceUnroll] for (int a = 0; a < 3; ++a) part[2 + a] += ri.e[6 + a] * W + ri.e[2] * Wa[a] + ri.e[6 + a] * ri.q[2];
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int bq = 0; bq < 3; ++bq) {
        if (bq >= a) {
            const int t = ab_index(a, bq);
            part[5 + t] += ri.e[9 + t] * W1 + W1ab[t] * ri.e[1] + ri.e[9 + t] * ri.q[1] + ri.e[3 + a] * W1a[bq] + W1a[a] * ri.e[3 + bq] + ri.e[3 + a] * ri.q[3 + bq];
            part[11 + t] += ri.ei[9 + t] * W1i + W1abi[t] * ri.ei[1] + ri.ei[9 + t] * ri.qs[1] + ri.ei[3 + a] * W1ai[bq] + W1ai[a] * ri.ei[3 + bq] + ri.ei[3 + a] * ri.qs[3 + bq];
        }
    }
}

// epi_main takes the batches in the matrix form, epi_fact those in the factored form (rows: the factors,
// X = diag(l~)); each skips the other's.
[ForceInline]
void epi_body(uint3 gid, uint tid, bool fact) {
    const uint slot = tid >> 5, lane = tid & 31u;
    const uint g = lane >> 3, k = lane & 7u;
    const uint b = gid.x * 2u + slot;
    const Batch bt = args.batches[min(b, args.num_batches - 1u)];
    const bool live = b < args.num_batches && batch_factored(bt) == fact;   // not uniform per group: no early return
    const float h = bt.P_h.w, hi = 1.5 * h;
    const uint nrows = batch_rows(bt);
    const uint n = live ? nrows : 0u;
    const uint rbase = fact ? bt.info2.y : bt.info.z;
    const float fn = float(max(bt.info.y, nrows));
    const float safe = 1.001 + (fn + 64.0) * UF;
    const float h2 = h * h;

    // ---- DOT rows 5g .. 5g + 4 and the error sums of order g, rows in increasing order
    float s[5][10], eo[10];
    [ForceUnroll] for (int a = 0; a < 5; ++a) [ForceUnroll] for (int j = 0; j < 10; ++j) s[a][j] = 0.0;
    [ForceUnroll] for (int j = 0; j < 10; ++j) eo[j] = 0.0;
    for (uint i = 0; i < n; ++i) {
        float va[5], dvb[10], edv[3], ea, ma;
        if (fact) {
            float* P = args.fphi + (rbase + i) * PHI_W + k;
            const float lam = args.fac_l[min(i, args.fac_r - 1u)];
            float mvp[4], evp[4], edvp[4];
            fac_row_err(P, abs(lam), safe, fn, mvp, evp, edvp);
            [ForceUnroll] for (int j = 0; j < 10; ++j) dvb[j] = lam * P[(PC_V + j) * 8];
            [ForceUnroll] for (int o = 0; o < 3; ++o) edv[o] = edvp[o];
            [ForceUnroll] for (int a = 0; a < 5; ++a) va[a] = P[(PC_V + 5u * g + uint(a)) * 8u];
            ea = evp[g]; ma = mvp[g];
        } else {
            float* P = args.phi + (rbase + i) * PHI_W + k;
            float* Q = args.dv + (rbase + i) * DV_W + k;
            [ForceUnroll] for (int j = 0; j < 10; ++j) dvb[j] = Q[j * 8];
            [ForceUnroll] for (int o = 0; o < 3; ++o) edv[o] = Q[(20 + o) * 8] * safe + fn * 1e-36;
            [ForceUnroll] for (int a = 0; a < 5; ++a) va[a] = P[(PC_V + 5u * g + uint(a)) * 8u];
            ea = P[(PC_EV + g) * 8u]; ma = P[(PC_MV + g) * 8u];
        }
        [ForceUnroll] for (int a = 0; a < 5; ++a) [ForceUnroll] for (int j = 0; j < 10; ++j) s[a][j] += va[a] * dvb[j];
        [ForceUnroll] for (int j = 0; j < 10; ++j) {
            const float adv = abs(dvb[j]), ed = edv[order_of(j)];
            eo[j] += ea * (adv + ed) + ma * ed + (fn + 1.0) * UF * ma * adv;
        }
    }

    // ---- enclosure sums: lanes g and g + 4 of the rows taken 8 at a time
    float pa[17], pb[17];
    [ForceUnroll] for (int t = 0; t < 17; ++t) { pa[t] = 0.0; pb[t] = 0.0; }
    for (uint i0 = 0; i0 < n; i0 += 8u) {
        const uint ia = i0 + g, ib = i0 + g + 4u;
        if (ia < n) epi_row_parts(fact ? row_factor(rbase + ia, ia, k, safe, fn) : row_matrix(rbase + ia, k, safe, fn), h, hi, h2, pa);
        if (ib < n) epi_row_parts(fact ? row_factor(rbase + ib, ib, k, safe, fn) : row_matrix(rbase + ib, k, safe, fn), h, hi, h2, pb);
    }

    // ---- outputs
    const uint fs = 8u * args.num_batches;
    if (live) {
        float* out = args.ep + b * 8u + k;
        [ForceUnroll] for (int a = 0; a < 5; ++a) [ForceUnroll] for (int j = 0; j < 10; ++j) out[(EP_DOT + (5u * g + uint(a)) * 10u + uint(j)) * fs] = s[a][j];
        [ForceUnroll] for (int j = 0; j < 10; ++j) out[(EP_EO + g * 10u + uint(j)) * fs] = eo[j];
    }
    // lane L of child k sits at gs_part[slot * 64 + L * 8 + k]
    [ForceUnroll] for (int t = 0; t < 17; ++t) {
        gs_part[slot * 64u + g * 8u + k][t] = pa[t];
        gs_part[slot * 64u + (g + 4u) * 8u + k][t] = pb[t];
    }
    GroupMemoryBarrierWithGroupSync();
    if (live && g == 0u) {
        float tot[17];
        [ForceUnroll] for (int t = 0; t < 17; ++t) tot[t] = 0.0;
        for (uint L = 0; L < 8u; ++L) {
            [ForceUnroll] for (int t = 0; t < 17; ++t) tot[t] += gs_part[slot * 64u + L * 8u + k][t];
        }
        float* o2 = args.ep + b * 8u + k;
        [ForceUnroll] for (int t = 0; t < 11; ++t) o2[(EP_PART + t) * fs] = tot[t];
        [ForceUnroll] for (int t = 0; t < 6; ++t) o2[(EP_PI + t) * fs] = tot[11 + t];
    }
}

[shader("compute")]
[numthreads(EPI_WG, 1, 1)]
void epi_main(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { epi_body(gid, tid, false); }

[shader("compute")]
[numthreads(EPI_WG, 1, 1)]
void epi_fact(uint3 gid : SV_GroupID, uint tid : SV_GroupIndex) { epi_body(gid, tid, true); }

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

// q(d) = c0 + b.d + 1/2 d^T M d evaluated in float, and a bound of its rounding error.
float quad_eval(float c0, float b[3], float M[3][3], float d[3], out float err) {
    float q = c0, qa = abs(c0);
    [ForceUnroll] for (int j = 0; j < 3; ++j) {
        q += b[j] * d[j]; qa += abs(b[j] * d[j]);
        [ForceUnroll] for (int l = 0; l < 3; ++l) { const float t = 0.5 * M[j][l] * d[j] * d[l]; q += t; qa += abs(t); }
    }
    err = 20.0 * UF * qa + 1e-36;
    return q;
}

// A lower bound of min over |d_k| <= h of q(d) = c0 + b.d + 1/2 d^T M d (M symmetric), in exact
// arithmetic on the given float coefficients. The minimum lies at a stationary point of q restricted
// to the relative interior of one of the box's 27 faces (the box, 6 faces, 12 edges, 8 corners): a
// corner gives its value; a face whose restricted matrix M_FF is not positive definite has its
// minimum on its own boundary (indefinite, or flat directions) and is skipped, as is a positive
// definite one whose minimiser is certainly outside it; otherwise, with x~ the computed stationary
// point, rho the residual (b + M x~)_F and l_min > 0 a lower bound of M_FF's smallest eigenvalue,
//   q >= q(x~) - |rho|^2 / (2 l_min)   on the whole face.
// Returns false when a face's definiteness cannot be decided in float (the test then just fails).
// Division is only assumed accurate to 2.5 ulp (<= 5u relative, the Vulkan guarantee for OpFDiv), and
// no square root is taken: |rho|_2 <= s = |rho|_1 serves in both the distance bound and the value bound.
bool quad_box_min_lb(float c0, float b[3], float M[3][3], float h, out float lb) {
    lb = 3.0e38;
    for (int pat = 0; pat < 27; ++pat) {
        const int st[3] = { pat % 3, (pat / 3) % 3, pat / 9 };    // 0 free, 1 -h, 2 +h
        float d[3];
        int F[3] = { 0, 0, 0 }, nf = 0;
        [ForceUnroll] for (int k = 0; k < 3; ++k) {
            d[k] = st[k] == 1 ? -h : (st[k] == 2 ? h : 0.0);
            if (st[k] == 0) { F[nf] = k; nf++; }
        }
        if (nf > 0) {
            float m[3][3], r[3];
            [ForceUnroll] for (int i = 0; i < 3; ++i) {
                r[i] = 0.0;
                [ForceUnroll] for (int j = 0; j < 3; ++j) m[i][j] = 0.0;
                if (i < nf) {
                    r[i] = -b[F[i]];
                    [ForceUnroll] for (int k = 0; k < 3; ++k) if (st[k] != 0) r[i] -= M[F[i]][k] * d[k];
                    [ForceUnroll] for (int j = 0; j < 3; ++j) if (j < nf) m[i][j] = M[F[i]][F[j]];
                }
            }
            // definiteness of M_FF by its leading minors, each with a bound of its rounding
            float lmin = 0.0, x[3] = { 0.0, 0.0, 0.0 };
            if (nf == 1) {
                if (!(m[0][0] > 0.0)) continue;                     // exact: not positive definite
                lmin = m[0][0];
                x[0] = r[0] / m[0][0];
            } else if (nf == 2) {
                if (!(m[0][0] > 0.0)) continue;
                const float det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
                const float edet = 3.0 * UF * (abs(m[0][0] * m[1][1]) + abs(m[0][1] * m[1][0])) + 1e-36;
                if (det < -edet) continue;
                if (!(det > edet)) return false;
                lmin = (det - edet) / ((m[0][0] + m[1][1]) * (1.0 + 4.0 * UF)) * (1.0 - 16.0 * UF);
                x[0] = (r[0] * m[1][1] - m[0][1] * r[1]) / det;
                x[1] = (m[0][0] * r[1] - r[0] * m[1][0]) / det;
            } else {
                if (!(m[0][0] > 0.0)) continue;
                const float m2 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
                const float em2 = 3.0 * UF * (abs(m[0][0] * m[1][1]) + abs(m[0][1] * m[1][0])) + 1e-36;
                if (m2 < -em2) continue;
                if (!(m2 > em2)) return false;
                const float c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
                const float c01 = m[1][0] * m[2][2] - m[1][2] * m[2][0];
                const float c02 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
                const float det = m[0][0] * c00 - m[0][1] * c01 + m[0][2] * c02;
                const float edet = 8.0 * UF * (abs(m[0][0]) * (abs(m[1][1] * m[2][2]) + abs(m[1][2] * m[2][1]))
                                             + abs(m[0][1]) * (abs(m[1][0] * m[2][2]) + abs(m[1][2] * m[2][0]))
                                             + abs(m[0][2]) * (abs(m[1][0] * m[2][1]) + abs(m[1][1] * m[2][0]))) + 1e-36;
                if (det < -edet) continue;
                if (!(det > edet)) return false;
                const float tr = (m[0][0] + m[1][1] + m[2][2]) * (1.0 + 2.0 * UF);
                lmin = (det - edet) / (0.25 * tr * tr * (1.0 + 4.0 * UF)) * (1.0 - 16.0 * UF);  // l1 l2 l3 = det, l2 l3 <= (tr/2)^2
                const float r0 = r[0], r1 = r[1], r2 = r[2];
                x[0] = (r0 * c00 - m[0][1] * (r1 * m[2][2] - m[1][2] * r2) + m[0][2] * (r1 * m[2][1] - m[1][1] * r2)) / det;
                x[1] = (m[0][0] * (r1 * m[2][2] - m[1][2] * r2) - r0 * c01 + m[0][2] * (m[1][0] * r2 - r1 * m[2][0])) / det;
                x[2] = (m[0][0] * (m[1][1] * r2 - r1 * m[2][1]) - m[0][1] * (m[1][0] * r2 - r1 * m[2][0]) + r0 * c02) / det;
            }
            if (!(lmin >= 1e-30)) return false;                     // also keeps the divisions by lmin normal
            [ForceUnroll] for (int i = 0; i < 3; ++i) if (i < nf) d[F[i]] = x[i];
            // residual of the restricted stationarity at d, with its rounding: s >= |rho|_1 >= |rho|_2
            float sr = 0.0;
            [ForceUnroll] for (int i = 0; i < 3; ++i) {
                if (i < nf) {
                    const int j = F[i];
                    float g = b[j], ga = abs(b[j]);
                    [ForceUnroll] for (int l = 0; l < 3; ++l) { g += M[j][l] * d[l]; ga += abs(M[j][l] * d[l]); }
                    sr += abs(g) + 5.0 * UF * ga + 1e-36;
                }
            }
            const float dist = sr / lmin * (1.0 + 16.0 * UF);                  // >= |x~ - x*|
            bool outside = false;
            [ForceUnroll] for (int i = 0; i < 3; ++i) if (i < nf && abs(x[i]) - dist > h) outside = true;
            if (outside) continue;
            float eq;
            const float q = quad_eval(c0, b, M, d, eq);
            const float v = q - eq - 0.5 * sr * dist;                            // >= |rho|^2 / (2 l_min)
            lb = min(lb, v - 4.0 * UF * abs(v) - 1e-36);
        } else {
            float eq;
            const float q = quad_eval(c0, b, M, d, eq);
            const float v = q - eq;
            lb = min(lb, v - 2.0 * UF * abs(v) - 1e-36);
        }
    }
    return true;
}

// Exclusion along w with the quadratic part of the model minimised over the box rather than bounded
// term by term (excl_dir). w . grad rho(c + d) >= q(d) - Etot on the cube, q built from the computed
// g, A, T; Etot covers their errors, the rounding of q's coefficients and the remainder r_up.
bool excl_quad(Eval ev, float T[3][3][3], float eT[3][3][3], float h, float3 w) {
    float c0 = 0.0, ec0 = 0.0, b[3], M[3][3], Etot = 0.0;
    [ForceUnroll] for (int k = 0; k < 3; ++k) { c0 += w[k] * ev.g[k]; ec0 += abs(w[k]) * ev.e_g[k] + 3.0 * UF * abs(w[k] * ev.g[k]); }
    Etot = ec0;
    [ForceUnroll] for (int j = 0; j < 3; ++j) {
        float bj = 0.0, eb = 0.0;
        [ForceUnroll] for (int k = 0; k < 3; ++k) { bj += w[k] * ev.A[k][j]; eb += abs(w[k]) * ev.e_A[k][j] + 3.0 * UF * abs(w[k] * ev.A[k][j]); }
        b[j] = bj;
        Etot += eb * h;
        [ForceUnroll] for (int l = 0; l < 3; ++l) {
            float mj = 0.0, em = 0.0;
            [ForceUnroll] for (int k = 0; k < 3; ++k) { mj += w[k] * T[k][j][l]; em += abs(w[k]) * eT[k][j][l] + 3.0 * UF * abs(w[k] * T[k][j][l]); }
            M[j][l] = mj;
            Etot += 0.5 * em * h * h;
        }
        Etot += abs(w[j]) * ev.r_up[j];
    }
    Etot *= 1.0 + 16.0 * UF;
    // The corners are candidates of both extremes: if one is within Etot of zero on a side, that
    // side cannot succeed, and the 27-face search is skipped (it is for most cubes and directions).
    bool try_min = true, try_max = true;
    [ForceUnroll] for (int c = 0; c < 8; ++c) {
        float d[3] = { (c & 1) != 0 ? h : -h, (c & 2) != 0 ? h : -h, (c & 4) != 0 ? h : -h };
        float eq;
        const float q = quad_eval(c0, b, M, d, eq);
        if (q <= Etot) try_min = false;
        if (q >= -Etot) try_max = false;
    }
    float lo, lo_neg;
    if (try_min && quad_box_min_lb(c0, b, M, h, lo) && lo > Etot) return true;
    if (!try_max) return false;
    // the maximum: the minimum of -q
    float nb[3], nM[3][3];
    [ForceUnroll] for (int j = 0; j < 3; ++j) { nb[j] = -b[j]; [ForceUnroll] for (int l = 0; l < 3; ++l) nM[j][l] = -M[j][l]; }
    return quad_box_min_lb(-c0, nb, nM, h, lo_neg) && lo_neg > Etot;
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

// Whether child s (centre c + (+-h/2)^3, s & 1: + side in x, s & 2: y, s & 4: z) is excluded by this
// cube's expansion alone: g + A d + 1/2 T[d,d] re-centred at the child's centre (with the error bounds
// of g, A, T carried through, and the rounding of the re-centring), the quadratic term over the child,
// and this cube's remainder bound r_up, which holds on all of the cube. Then the exclusion tests the
// cube itself got: the seven directions and the Newton distance. Skips an evaluation for 27-38% of the
// children on the test systems; no child it excludes was ever undecided by its own evaluation.
bool child_excluded(Eval P, float T[3][3][3], float eT[3][3][3], float qT[3], float h, uint s, float safe) {
    const float hc = 0.5 * h;                   // exact; every product with it below is exact too
    float d[3];
    d[0] = (s & 1u) != 0u ? hc : -hc;
    d[1] = (s & 2u) != 0u ? hc : -hc;
    d[2] = (s & 4u) != 0u ? hc : -hc;
    Eval C = P;
    [ForceUnroll] for (int i = 0; i < 3; ++i) {
        float v = P.g[i], S = abs(P.g[i]), e = P.e_g[i];
        [ForceUnroll] for (int j = 0; j < 3; ++j) {
            const float t = P.A[i][j] * d[j];
            v += t; S += abs(t); e += P.e_A[i][j] * hc;
            [ForceUnroll] for (int l = 0; l < 3; ++l) {
                const float t2 = 0.5 * T[i][j][l] * d[j] * d[l];
                v += t2; S += abs(t2); e += 0.5 * eT[i][j][l] * hc * hc;
            }
        }
        C.g[i] = v;
        C.e_g[i] = e * (1.0 + 16.0 * UF) + 17.0 * UF * S + 13.0e-36;
        [ForceUnroll] for (int j = 0; j < 3; ++j) {
            float a = P.A[i][j], Sa = abs(P.A[i][j]), ea = P.e_A[i][j];
            [ForceUnroll] for (int l = 0; l < 3; ++l) {
                const float t = T[i][j][l] * d[l];
                a += t; Sa += abs(t); ea += eT[i][j][l] * hc;
            }
            C.A[i][j] = a;
            C.e_A[i][j] = ea * (1.0 + 8.0 * UF) + 9.0 * UF * Sa + 4.0e-36;
        }
        C.q_up[i] = 0.5 * hc * hc * qT[i] * safe;
        C.r_up[i] = P.r_up[i];
    }
    float val[3], vec[3][3];
    eigen3(C.A, val, vec);
    float3 dirs[7];
    dirs[0] = float3(1, 0, 0); dirs[1] = float3(0, 1, 0); dirs[2] = float3(0, 0, 1);
    int nd = 3;
    const float gn = sqrt(C.g[0] * C.g[0] + C.g[1] * C.g[1] + C.g[2] * C.g[2]);
    if (gn > 0.0) { dirs[nd] = float3(C.g[0], C.g[1], C.g[2]) / gn; nd++; }
    [ForceUnroll] for (int e = 0; e < 3; ++e) { dirs[nd] = float3(vec[0][e], vec[1][e], vec[2][e]); nd++; }
    for (int t = 0; t < nd; ++t) if (excl_dir(C, hc, dirs[t])) return true;
    const Lin lin = linear(C, val, vec);
    if (lin.ok && newton_dist(C, lin, hc)) return true;
    for (int t = 0; t < nd; ++t) if (excl_quad(C, T, eT, hc, dirs[t])) return true;     // the costly one last
    return false;
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
    const float fn = float(max(n, batch_rows(bt)));       // terms per sum: the GEMM's and the rows'
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
    // third derivatives: kept for the children's test, and their sums of |T| + e_T for the enclosure
    float T[3][3][3], eT[3][3][3], qT[3], lT[6];
    [ForceUnroll] for (int a = 0; a < 3; ++a) qT[a] = 0.0;
    [ForceUnroll] for (int t = 0; t < 6; ++t) lT[t] = 0.0;
    [ForceUnroll] for (int a = 0; a < 3; ++a) [ForceUnroll] for (int c = 0; c < 3; ++c) [ForceUnroll] for (int d = 0; d < 3; ++d) {
        const float t0 = DOT(M3(a, c, d), 0), t1 = DOT(M2(a, c), M1(d)), t2 = DOT(M2(a, d), M1(c)), t3 = DOT(M2(c, d), M1(a));
        T[a][c][d] = 2.0 * (t0 + t1 + t2 + t3);
        eT[a][c][d] = (2.0 * (EDOT(M3(a, c, d), 0) + EDOT(M2(a, c), M1(d)) + EDOT(M2(a, d), M1(c)) + EDOT(M2(c, d), M1(a)))
                       + 8.0 * UF * (abs(t0) + abs(t1) + abs(t2) + abs(t3))) * 1.001;
        qT[a] += abs(T[a][c][d]) + eT[a][c][d];                               // sum over (c, d)
        if (a <= c) lT[ab_index(a, c)] += abs(T[a][c][d]) + eT[a][c][d];      // sum over d, a <= c
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
    float3 dirs[7];
    dirs[0] = float3(1, 0, 0); dirs[1] = float3(0, 1, 0); dirs[2] = float3(0, 0, 1);
    int nd = 3;
    {
        const float gn = sqrt(ev.g[0] * ev.g[0] + ev.g[1] * ev.g[1] + ev.g[2] * ev.g[2]);
        if (gn > 0.0) { dirs[nd] = float3(ev.g[0], ev.g[1], ev.g[2]) / gn; nd++; }
        [ForceUnroll] for (int e = 0; e < 3; ++e) { dirs[nd] = float3(vec[0][e], vec[1][e], vec[2][e]); nd++; }
    }
    if (!decided) {
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
        // the exact quadratic test, the costliest, once the cheap ones have failed
        for (int t = 0; t < nd && !decided; ++t) if (excl_quad(ev, T, eT, h, dirs[t])) { outc = OUT_DISCARD; why = WHY_EXCL; decided = true; }
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
    uint children = 0u;
    if (outc == OUT_SPLIT) {
        for (uint s = 0; s < 8u; ++s) if (!child_excluded(ev, T, eT, qT, h, s, safe)) children |= 1u << s;
        if (children == 0u) { outc = OUT_DISCARD; why = WHY_CHILDREN; }
    }
    if (tried_infl) outc |= OUT_INFLATED;
    args.outcome[cube] = outc | (why << 8) | (children << 16);
    if ((args.flags & 1u) != 0u) {
        float* dbg = args.debug + cube * 16u;
        dbg[0] = ev.rho; dbg[1] = ev.e_rho; dbg[2] = ev.rho_up;
        dbg[3] = ev.g[0]; dbg[4] = ev.g[1]; dbg[5] = ev.g[2];
        dbg[6] = ev.e_g[0]; dbg[7] = ev.e_g[1]; dbg[8] = ev.e_g[2];
        dbg[9] = ev.A[0][0]; dbg[10] = ev.A[1][1]; dbg[11] = ev.A[2][2];
        dbg[12] = ev.r_up[0]; dbg[13] = ev.dH_up[0][0]; dbg[14] = fn; dbg[15] = ev.q_up[0];
    }
}

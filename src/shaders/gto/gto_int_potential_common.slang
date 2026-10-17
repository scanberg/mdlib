// gto_int_potential_common.slang
//
// What the two potential evaluators share (eval_gto_int_potential.slang on a grid,
// eval_gto_int_potential_xyz.slang at points): the fixed point accumulator, the Boys function from
// the table and the Hermite recursion of one order. Included after the evaluator defines `args`
// (with .boys and .coeffs) and the GTO_INT_BOYS_* constants.

// ---------------------------------------------------------------------------
// 64 bit fixed point, 32 fractional bits, as (lo, hi) with hi two's complement.
// ---------------------------------------------------------------------------
uint2 gto_int_fixed_add(uint2 acc, float v) {
    // Magnitude into (integer, fraction); a - floor(a) is exact for a >= 0, and the fraction
    // times 2^32 stays below 2^32 (at most 2^32 - 2^8), so round() cannot overflow.
    const float a  = min(abs(v), 2147483520.0f);
    const float fl = floor(a);
    const uint  ih = uint(fl);
    const uint  il = uint(round((a - fl) * 4294967296.0f));
    uint2 r;
    if (v >= 0.0f) {
        r.x = acc.x + il;
        r.y = acc.y + ih + ((r.x < acc.x) ? 1u : 0u);
    } else {
        r.x = acc.x - il;
        r.y = acc.y - ih - ((acc.x < il) ? 1u : 0u);
    }
    return r;
}

float gto_int_fixed_to_float(uint2 acc) {
    return float(int(acc.y)) + float(acc.x) * 2.3283064365386963e-10f;
}

// ---------------------------------------------------------------------------
// s[n] = (-2p)^n F_n(T), n = 0..L
// ---------------------------------------------------------------------------
// Below TMAX: 5 term Taylor expansion of F_L around the nearest tabulated T (|dT| <= 0.05, error
// below 3e-9 relative), then downward recursion, which is stable. Above: the asymptotic F_0 and
// upward recursion, stable there for every order this kernel handles.
void gto_int_boys_scaled<let L : int>(out float s[GTO_INT_NS], float T, float p) {
    float F[GTO_INT_NS];
    const float e = exp(-T);
    if (T < GTO_INT_BOYS_TMAX) {
        const int   k = int(T * GTO_INT_BOYS_INV_DT + 0.5f);
        const float d = float(k) * GTO_INT_BOYS_DT - T;
        const uint  row = uint(k) * GTO_INT_BOYS_NN + uint(L);
        float f = args.boys[row + 4u] * (1.0f / 24.0f);
        f = args.boys[row + 3u] * (1.0f / 6.0f) + d * f;
        f = args.boys[row + 2u] * 0.5f + d * f;
        f = args.boys[row + 1u] + d * f;
        f = args.boys[row + 0u] + d * f;
        F[L] = f;
        [ForceUnroll] for (int n = L - 1; n >= 0; --n) {
            F[n] = (2.0f * T * F[n + 1] + e) * (1.0f / float(2 * n + 1));
        }
    } else {
        const float inv2T = 0.5f / T;
        F[0] = 0.5f * sqrt(3.14159265358979f / T);
        [ForceUnroll] for (int n = 0; n < L; ++n) {
            F[n + 1] = (float(2 * n + 1) * F[n] - e) * inv2T;
        }
    }
    const float m2p = -2.0f * p;
    float scale = 1.0f;
    [ForceUnroll] for (int n = 0; n <= L; ++n) {
        s[n] = scale * F[n];
        scale *= m2p;
    }
    [ForceUnroll] for (int n = L + 1; n < GTO_INT_NS; ++n) s[n] = 0.0f;
}

float gto_int_hermite<let L : int>(float3 X, float s[GTO_INT_NS], uint o) {
    if (L == 0) return gto_int_hermite_L0(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 1) return gto_int_hermite_L1(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 2) return gto_int_hermite_L2(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 3) return gto_int_hermite_L3(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 4) return gto_int_hermite_L4(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 5) return gto_int_hermite_L5(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 6) return gto_int_hermite_L6(X.x, X.y, X.z, s, args.coeffs, o);
    if (L == 7) return gto_int_hermite_L7(X.x, X.y, X.z, s, args.coeffs, o);
    return gto_int_hermite_L8(X.x, X.y, X.z, s, args.coeffs, o);
}

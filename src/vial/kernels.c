// The interpreter's kernels: elementwise operators with NumPy's broadcasting, conversions, indexing, stacking and the
// kinds. Built with floating-point contraction off: where a multiply and an add fuse is the compiler's decision, made
// once in the IR, so every backend gives the same bits.

// No contraction: where a multiply and an add fuse is the IR's decision (VIAL_IR_FMA). GCC and Clang get
// -ffp-contract=off for this file from CMake; MSVC's /fp:contract is taken back here, which it allows per file.
#if defined(_MSC_VER) && !defined(__clang__)
#pragma fp_contract(off)
#endif

#include "interp.h"
#include "ast.h"
#include "infer_internal.h"

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>

#include <math.h>
#include <stdint.h>
#include <string.h>

// ### CONVERSION ###

static double get_f64(const vial_val_t* a, size_t i) {
    switch (a->dtype) {
    case VIAL_DT_BOOL:    return ((const uint8_t*)a->data)[i];
    case VIAL_DT_INT:     return ((const int32_t*)a->data)[i];
    case VIAL_DT_FLOAT:   return ((const float*)a->data)[i];
    case VIAL_DT_FLOAT64: return ((const double*)a->data)[i];
    default:              return 0.0;
    }
}

vial_val_t vial_val_as(const vial_val_t* a, uint8_t dtype, md_allocator_i* alloc) {
    if (a->dtype == dtype) return *a;
    const bool num_in = a->dtype == VIAL_DT_BOOL || a->dtype == VIAL_DT_INT || a->dtype == VIAL_DT_FLOAT || a->dtype == VIAL_DT_FLOAT64;
    const bool num_out = dtype == VIAL_DT_BOOL || dtype == VIAL_DT_INT || dtype == VIAL_DT_FLOAT || dtype == VIAL_DT_FLOAT64;
    if (!num_in || !num_out) return *a;
    vial_val_t o = vial_val_make(dtype, a->ndim, a->dims, alloc);
    const size_t n = vial_val_count(a);
    switch (dtype) {
    case VIAL_DT_BOOL:    for (size_t i = 0; i < n; ++i) ((uint8_t*)o.data)[i] = get_f64(a, i) != 0.0; break;
    case VIAL_DT_INT:
        // Toward zero, as int() does; a float that does not fit is missing
        for (size_t i = 0; i < n; ++i) {
            const double d = get_f64(a, i);
            if (d >= -2147483648.0 && d < 2147483648.0) ((int32_t*)o.data)[i] = (int32_t)d;
            else vial_val_set_missing(&o, i, alloc);
        }
        break;
    case VIAL_DT_FLOAT:
        if (a->dtype == VIAL_DT_INT) for (size_t i = 0; i < n; ++i) ((float*)o.data)[i] = (float)((const int32_t*)a->data)[i];
        else for (size_t i = 0; i < n; ++i) ((float*)o.data)[i] = (float)get_f64(a, i);
        break;
    case VIAL_DT_FLOAT64: for (size_t i = 0; i < n; ++i) ((double*)o.data)[i] = get_f64(a, i); break;
    }
    if (a->missing) for (size_t i = 0; i < n; ++i) if (vial_val_is_missing(a, i)) vial_val_set_missing(&o, i, alloc);
    return o;
}

bool vial_k_convert(vial_kctx_t* k, const vial_val_t* a, uint8_t dtype, vial_val_t* out) {
    *out = vial_val_as(a, dtype, k->alloc);
    return true;
}

// ### BROADCASTING ###

typedef struct bc_t {
    uint8_t  nd;
    uint32_t dims[VIAL_MAX_DIMS];
    size_t   sa[VIAL_MAX_DIMS], sb[VIAL_MAX_DIMS], sc[VIAL_MAX_DIMS];   // element strides, 0 along broadcast dimensions
    size_t   n;
} bc_t;

static void strides(const vial_val_t* v, uint8_t nd, const uint32_t* dims, size_t* st) {
    size_t s = 1;
    for (int d = nd - 1; d >= 0; --d) {
        const int vd = d - (nd - v->ndim);
        if (vd < 0) { st[d] = 0; continue; }
        st[d] = v->dims[vd] == 1 && dims[d] != 1 ? 0 : s;
        s *= v->dims[vd];
    }
}

// NumPy's rule: aligned on the right, each pair of extents equal or one of them 1
static bool bc_init(vial_kctx_t* k, bc_t* b, const vial_val_t* x, const vial_val_t* y, const vial_val_t* z) {
    memset(b, 0, sizeof(*b));
    const vial_val_t* v[3] = {x, y, z};
    for (int i = 0; i < 3; ++i) if (v[i] && v[i]->ndim > b->nd) b->nd = v[i]->ndim;
    for (uint8_t d = 0; d < b->nd; ++d) b->dims[d] = 1;
    for (int i = 0; i < 3; ++i) {
        if (!v[i]) continue;
        for (uint8_t d = 0; d < v[i]->ndim; ++d) {
            const uint32_t e = v[i]->dims[d];
            uint32_t* o = &b->dims[b->nd - v[i]->ndim + d];
            if (*o == 1) *o = e;
            else if (e != 1 && e != *o) return vial_rt_fail(k, "R0103", "extents %u and %u do not broadcast", *o, e);
        }
    }
    strides(x, b->nd, b->dims, b->sa);
    if (y) strides(y, b->nd, b->dims, b->sb);
    if (z) strides(z, b->nd, b->dims, b->sc);
    b->n = 1;
    for (uint8_t d = 0; d < b->nd; ++d) {
        if (b->dims[d] && b->n > SIZE_MAX / b->dims[d]) { b->n = SIZE_MAX; break; }
        b->n *= b->dims[d];
    }
    return vial_k_room(k, b->n);
}

// Calls body(io, ia, ib, ic) for every element of the broadcast shape, in order
#define BC_FOR(b, ...)                                                                                  \
    do {                                                                                                \
        uint32_t idx_[VIAL_MAX_DIMS] = {0};                                                             \
        size_t ia = 0, ib = 0, ic = 0;                                                                  \
        for (size_t io = 0; io < (b).n; ++io) {                                                         \
            (void)ia; (void)ib; (void)ic;                                                               \
            __VA_ARGS__;                                                                                \
            for (int d_ = (int)(b).nd - 1; d_ >= 0; --d_) {                                             \
                if (++idx_[d_] < (b).dims[d_]) { ia += (b).sa[d_]; ib += (b).sb[d_]; ic += (b).sc[d_]; break; } \
                ia -= (b).sa[d_] * (idx_[d_] - 1); ib -= (b).sb[d_] * (idx_[d_] - 1); ic -= (b).sc[d_] * (idx_[d_] - 1); \
                idx_[d_] = 0;                                                                           \
            }                                                                                           \
        }                                                                                               \
    } while (0)

static void bc_missing(const bc_t* b, const vial_val_t* x, const vial_val_t* y, const vial_val_t* z, vial_val_t* out, md_allocator_i* alloc) {
    if (!(x && x->missing) && !(y && y->missing) && !(z && z->missing)) return;
    BC_FOR(*b, {
        if ((x && vial_val_is_missing(x, ia)) || (y && vial_val_is_missing(y, ib)) || (z && vial_val_is_missing(z, ic))) vial_val_set_missing(out, io, alloc);
    });
}

static bool is_num(uint8_t dt) { return dt == VIAL_DT_BOOL || dt == VIAL_DT_INT || dt == VIAL_DT_FLOAT || dt == VIAL_DT_FLOAT64; }

// ### ARITHMETIC ###

// NumPy's divmod for floats: the floor quotient and a remainder with the divisor's sign
static float floordiv_f(float a, float b, float* mod_out) {
    float mod = fmodf(a, b);
    if (b == 0.0f) { if (mod_out) *mod_out = mod; return a / b; }
    float div = (a - mod) / b;
    if (mod != 0.0f) {
        if ((b < 0.0f) != (mod < 0.0f)) { mod += b; div -= 1.0f; }
    } else {
        mod = copysignf(0.0f, b);
    }
    float fl;
    if (div != 0.0f) { fl = floorf(div); if (div - fl > 0.5f) fl += 1.0f; }
    else fl = copysignf(0.0f, a / b);
    if (mod_out) *mod_out = mod;
    return fl;
}

static double floordiv_d(double a, double b, double* mod_out) {
    double mod = fmod(a, b);
    if (b == 0.0) { if (mod_out) *mod_out = mod; return a / b; }
    double div = (a - mod) / b;
    if (mod != 0.0) {
        if ((b < 0.0) != (mod < 0.0)) { mod += b; div -= 1.0; }
    } else {
        mod = copysign(0.0, b);
    }
    double fl;
    if (div != 0.0) { fl = floor(div); if (div - fl > 0.5) fl += 1.0; }
    else fl = copysign(0.0, a / b);
    if (mod_out) *mod_out = mod;
    return fl;
}

// x ** y in float: computed in double and rounded once, so the result does not depend on a libm's powf, which
// promises less than half an ulp; x ** 2 is x * x, the same bits (a float's square is exact in double)
static inline float pow_f(float a, float b) { return b == 2.0f ? a * a : (float)pow((double)a, (double)b); }
static inline double pow_d(double a, double b) { return b == 2.0 ? a * a : pow(a, b); }

double vial_k_arith(uint16_t op, uint8_t dtype, double a, double b) {
    if (dtype == VIAL_DT_FLOAT) {
        const float x = (float)a, y = (float)b;
        float m;
        switch (op) {
        case VIAL_OP_PLUS:    return x + y;
        case VIAL_OP_MINUS:   return x - y;
        case VIAL_OP_STAR:    return x * y;
        case VIAL_OP_SLASH:   return x / y;
        case VIAL_OP_DSLASH:  return floordiv_f(x, y, NULL);
        case VIAL_OP_PERCENT: floordiv_f(x, y, &m); return m;
        case VIAL_OP_POW:     return pow_f(x, y);
        default:              return NAN;
        }
    }
    double m;
    switch (op) {
    case VIAL_OP_PLUS:    return a + b;
    case VIAL_OP_MINUS:   return a - b;
    case VIAL_OP_STAR:    return a * b;
    case VIAL_OP_SLASH:   return a / b;
    case VIAL_OP_DSLASH:  return floordiv_d(a, b, NULL);
    case VIAL_OP_PERCENT: floordiv_d(a, b, &m); return m;
    case VIAL_OP_POW:     return pow_d(a, b);
    default:              return NAN;
    }
}

bool vial_k_binary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    const uint8_t dt = k->type->dtype;
    if (dt == VIAL_DT_DRAWING) {
        // Drawings combine: one per element of the broadcast shape
        bc_t bc;
        if (!bc_init(k, &bc, a, b, NULL)) return false;
        *out = (vial_val_t){.dtype = VIAL_DT_DRAWING, .ndim = bc.nd};
        memcpy(out->dims, bc.dims, sizeof(bc.dims));
        return true;
    }
    if (!is_num(a->dtype) || !is_num(b->dtype)) return vial_rt_fail(k, "R0103", "arithmetic on what is not a number");
    if (op == VIAL_OP_AT) return vial_k_matmul(k, a, b, out);
    // Computed in the result's dtype, as NumPy does once both sides are promoted
    const uint8_t cd = dt == VIAL_DT_BOOL ? VIAL_DT_INT : dt;
    const vial_val_t x = vial_val_as(a, cd, k->alloc), y = vial_val_as(b, cd, k->alloc);
    bc_t bc;
    if (!bc_init(k, &bc, &x, &y, NULL)) return false;
    vial_val_t o; if (!vial_k_make(k, cd, bc.nd, bc.dims, &o)) return false;
    if (cd == VIAL_DT_FLOAT) {
        const float* A = x.data; const float* B = y.data; float* O = o.data;
        switch (op) {
        case VIAL_OP_PLUS:    BC_FOR(bc, O[io] = A[ia] + B[ib]); break;
        case VIAL_OP_MINUS:   BC_FOR(bc, O[io] = A[ia] - B[ib]); break;
        case VIAL_OP_STAR:    BC_FOR(bc, O[io] = A[ia] * B[ib]); break;
        case VIAL_OP_SLASH:   BC_FOR(bc, O[io] = A[ia] / B[ib]); break;
        case VIAL_OP_DSLASH:  BC_FOR(bc, O[io] = floordiv_f(A[ia], B[ib], NULL)); break;
        case VIAL_OP_PERCENT: BC_FOR(bc, { float m; floordiv_f(A[ia], B[ib], &m); O[io] = m; }); break;
        case VIAL_OP_POW:     BC_FOR(bc, O[io] = pow_f(A[ia], B[ib])); break;
        default: return vial_rt_fail(k, "R0103", "`%s` on floats", vial_op_str((vial_op_t)op));
        }
    } else if (cd == VIAL_DT_FLOAT64) {
        const double* A = x.data; const double* B = y.data; double* O = o.data;
        switch (op) {
        case VIAL_OP_PLUS:    BC_FOR(bc, O[io] = A[ia] + B[ib]); break;
        case VIAL_OP_MINUS:   BC_FOR(bc, O[io] = A[ia] - B[ib]); break;
        case VIAL_OP_STAR:    BC_FOR(bc, O[io] = A[ia] * B[ib]); break;
        case VIAL_OP_SLASH:   BC_FOR(bc, O[io] = A[ia] / B[ib]); break;
        case VIAL_OP_DSLASH:  BC_FOR(bc, O[io] = floordiv_d(A[ia], B[ib], NULL)); break;
        case VIAL_OP_PERCENT: BC_FOR(bc, { double m; floordiv_d(A[ia], B[ib], &m); O[io] = m; }); break;
        case VIAL_OP_POW:     BC_FOR(bc, O[io] = pow_d(A[ia], B[ib])); break;
        default: return vial_rt_fail(k, "R0103", "`%s` on floats", vial_op_str((vial_op_t)op));
        }
    } else {
        // Ints: computed in 64 bits; a result outside 32 is missing (R0102), as is a division by zero
        const int32_t* A = x.data; const int32_t* B = y.data; int32_t* O = o.data;
        BC_FOR(bc, {
            const int64_t p = A[ia], q = B[ib];
            int64_t r = 0;
            bool ok = true;
            switch (op) {
            case VIAL_OP_PLUS:    r = p + q; break;
            case VIAL_OP_MINUS:   r = p - q; break;
            case VIAL_OP_STAR:    r = p * q; break;
            case VIAL_OP_DSLASH:  if (!q) ok = false; else { r = p / q; if ((p % q) && ((p < 0) != (q < 0))) r -= 1; } break;
            case VIAL_OP_PERCENT: if (!q) ok = false; else { r = p % q; if (r && ((r < 0) != (q < 0))) r += q; } break;
            case VIAL_OP_POW:
                if (q < 0) { ok = false; break; }
                r = 1;
                for (int64_t e = 0; e < q && ok; ++e) { r *= p; if (r > INT32_MAX || r < INT32_MIN) ok = false; if (r == 0 || r == 1) break; }
                break;
            case VIAL_OP_AMP:     r = p & q; break;
            case VIAL_OP_PIPE:    r = p | q; break;
            case VIAL_OP_CARET:   r = p ^ q; break;
            case VIAL_OP_LSHIFT:  if (q < 0 || q > 31) ok = false; else r = (int64_t)((uint64_t)p << q); break;
            case VIAL_OP_RSHIFT:  if (q < 0) ok = false; else r = q > 62 ? (p < 0 ? -1 : 0) : p >> q; break;
            default: ok = false; break;
            }
            if (ok && (r > INT32_MAX || r < INT32_MIN)) ok = false;
            if (ok) O[io] = (int32_t)r; else { O[io] = 0; vial_val_set_missing(&o, io, k->alloc); }
        });
    }
    bc_missing(&bc, &x, &y, NULL, &o, k->alloc);
    *out = o;
    return true;
}

bool vial_k_fma(vial_kctx_t* k, const vial_val_t* a, const vial_val_t* b, const vial_val_t* c, vial_val_t* out) {
    const uint8_t cd = k->type->dtype;
    const vial_val_t x = vial_val_as(a, cd, k->alloc), y = vial_val_as(b, cd, k->alloc), z = vial_val_as(c, cd, k->alloc);
    bc_t bc;
    if (!bc_init(k, &bc, &x, &y, &z)) return false;
    vial_val_t o; if (!vial_k_make(k, cd, bc.nd, bc.dims, &o)) return false;
    if (cd == VIAL_DT_FLOAT) { const float *A = x.data, *B = y.data, *C = z.data; float* O = o.data; BC_FOR(bc, O[io] = fmaf(A[ia], B[ib], C[ic])); }
    else if (cd == VIAL_DT_FLOAT64) { const double *A = x.data, *B = y.data, *C = z.data; double* O = o.data; BC_FOR(bc, O[io] = fma(A[ia], B[ib], C[ic])); }
    else return vial_rt_fail(k, "R0103", "fma on ints");
    bc_missing(&bc, &x, &y, &z, &o, k->alloc);
    *out = o;
    return true;
}

bool vial_k_unary(vial_kctx_t* k, uint16_t op, const vial_val_t* a, vial_val_t* out) {
    const uint8_t cd = k->type->dtype;
    const vial_val_t x = vial_val_as(a, cd, k->alloc);
    vial_val_t o; if (!vial_k_make(k, cd, x.ndim, x.dims, &o)) return false;
    const size_t n = vial_val_count(&x);
    if (op == VIAL_OP_PLUS) { *out = x; return true; }
    switch (cd) {
    case VIAL_DT_FLOAT:   for (size_t i = 0; i < n; ++i) ((float*)o.data)[i] = -((const float*)x.data)[i]; break;
    case VIAL_DT_FLOAT64: for (size_t i = 0; i < n; ++i) ((double*)o.data)[i] = -((const double*)x.data)[i]; break;
    case VIAL_DT_INT:
        for (size_t i = 0; i < n; ++i) {
            const int32_t v = ((const int32_t*)x.data)[i];
            if (op == VIAL_OP_TILDE) ((int32_t*)o.data)[i] = ~v;
            else if (v == INT32_MIN) vial_val_set_missing(&o, i, k->alloc);
            else ((int32_t*)o.data)[i] = -v;
        }
        break;
    default: return vial_rt_fail(k, "R0103", "a sign on what is not a number");
    }
    if (x.missing) for (size_t i = 0; i < n; ++i) if (vial_val_is_missing(&x, i)) vial_val_set_missing(&o, i, k->alloc);
    *out = o;
    return true;
}

bool vial_k_not(vial_kctx_t* k, const vial_val_t* a, vial_val_t* out) {
    if (a->dtype != VIAL_DT_BOOL) return vial_rt_fail(k, "R0103", "not of what is not a bool");
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, a->ndim, a->dims, &o)) return false;
    const size_t n = vial_val_count(a);
    for (size_t i = 0; i < n; ++i) ((uint8_t*)o.data)[i] = !((const uint8_t*)a->data)[i];
    o.missing = a->missing;
    *out = o;
    return true;
}

bool vial_k_logic(vial_kctx_t* k, uint16_t kw, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    if (a->dtype != VIAL_DT_BOOL || b->dtype != VIAL_DT_BOOL) return vial_rt_fail(k, "R0103", "and/or of what is not a bool");
    bc_t bc;
    if (!bc_init(k, &bc, a, b, NULL)) return false;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, bc.nd, bc.dims, &o)) return false;
    const uint8_t *A = a->data, *B = b->data;
    uint8_t* O = o.data;
    if (kw == VIAL_KW_AND) BC_FOR(bc, O[io] = A[ia] && B[ib]);
    else BC_FOR(bc, O[io] = A[ia] || B[ib]);
    bc_missing(&bc, a, b, NULL, &o, k->alloc);
    *out = o;
    return true;
}

// ### COMPARISONS ###

static bool sets_equal(const vial_sets_t* a, uint32_t i, const vial_sets_t* b, uint32_t j) {
    const uint32_t n = vial_sets_size(a, i);
    return n == vial_sets_size(b, j) && (!n || !memcmp(vial_sets_items(a, i), vial_sets_items(b, j), n * sizeof(uint32_t)));
}

bool vial_k_compare(vial_kctx_t* k, uint16_t cmp, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    if ((cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) && b->dtype == VIAL_DT_RANGE) {
        // a..b includes both ends
        const double lo = ((const double*)b->data)[0], hi = ((const double*)b->data)[1];
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, a->ndim, a->dims, &o)) return false;
        const size_t n = vial_val_count(a);
        for (size_t i = 0; i < n; ++i) {
            const double v = get_f64(a, i);
            const bool r = v >= lo && v <= hi;
            ((uint8_t*)o.data)[i] = cmp == VIAL_CMP_IN ? r : !r;
        }
        o.missing = a->missing;
        *out = o;
        return true;
    }
    if (cmp == VIAL_CMP_IN || cmp == VIAL_CMP_NOT_IN) {
        // Each of a's values among b's, as NumPy's isin
        if (!is_num(a->dtype) || !is_num(b->dtype)) return vial_rt_fail(k, "R0103", "`in` of what is not a number");
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, a->ndim, a->dims, &o)) return false;
        const size_t n = vial_val_count(a), m = vial_val_count(b);
        for (size_t i = 0; i < n; ++i) {
            const double v = get_f64(a, i);
            bool found = false;
            for (size_t j = 0; j < m && !found; ++j) found = !vial_val_is_missing(b, j) && get_f64(b, j) == v;
            ((uint8_t*)o.data)[i] = cmp == VIAL_CMP_IN ? found : !found;
        }
        o.missing = a->missing;
        *out = o;
        return true;
    }
    bc_t bc;
    if (!bc_init(k, &bc, a, b, NULL)) return false;
    vial_val_t o; if (!vial_k_make(k, VIAL_DT_BOOL, bc.nd, bc.dims, &o)) return false;
    uint8_t* O = o.data;
    if (is_num(a->dtype) && is_num(b->dtype)) {
        // In doubles, which hold every int and float of 32 bits exactly
        const vial_val_t x = vial_val_as(a, VIAL_DT_FLOAT64, k->alloc), y = vial_val_as(b, VIAL_DT_FLOAT64, k->alloc);
        const double *A = x.data, *B = y.data;
        switch (cmp) {
        case VIAL_CMP_LT: BC_FOR(bc, O[io] = A[ia] < B[ib]); break;
        case VIAL_CMP_GT: BC_FOR(bc, O[io] = A[ia] > B[ib]); break;
        case VIAL_CMP_LE: BC_FOR(bc, O[io] = A[ia] <= B[ib]); break;
        case VIAL_CMP_GE: BC_FOR(bc, O[io] = A[ia] >= B[ib]); break;
        case VIAL_CMP_EQ: BC_FOR(bc, O[io] = A[ia] == B[ib]); break;
        case VIAL_CMP_NE: BC_FOR(bc, O[io] = A[ia] != B[ib]); break;
        default: return vial_rt_fail(k, "R0103", "a comparison not evaluated yet");
        }
    } else if (a->dtype == VIAL_DT_STR && b->dtype == VIAL_DT_STR && (cmp == VIAL_CMP_EQ || cmp == VIAL_CMP_NE)) {
        const str_t *A = a->data, *B = b->data;
        BC_FOR(bc, O[io] = str_eq(A[ia], B[ib]) == (cmp == VIAL_CMP_EQ));
    } else if (a->dtype == VIAL_DT_SET && b->dtype == VIAL_DT_SET && (cmp == VIAL_CMP_EQ || cmp == VIAL_CMP_NE)) {
        BC_FOR(bc, O[io] = sets_equal(a->sets, (uint32_t)ia, b->sets, (uint32_t)ib) == (cmp == VIAL_CMP_EQ));
    } else {
        return vial_rt_fail(k, "R0103", "a comparison of what does not compare");
    }
    bc_missing(&bc, a, b, NULL, &o, k->alloc);
    *out = o;
    return true;
}

// ### SETS ###

bool vial_k_setop(vial_kctx_t* k, uint16_t sub, const vial_val_t* a, const vial_val_t* b, vial_val_t* out) {
    if (a->dtype != VIAL_DT_SET || b->dtype != VIAL_DT_SET) return vial_rt_fail(k, "R0103", "set algebra on what is not a set");
    bc_t bc;
    if (!bc_init(k, &bc, a, b, NULL)) return false;
    // Element by element; each result a set of its own, gathered into one array
    vial_sets_t** parts = md_alloc(k->alloc, (bc.n ? bc.n : 1) * sizeof(vial_sets_t*));
    uint32_t total = 0;
    BC_FOR(bc, {
        const uint32_t i = (uint32_t)ia, j = (uint32_t)ib;
        switch (sub) {
        case VIAL_KW_AND:  parts[io] = vial_sets_and(a->sets, i, b->sets, j, k->alloc); break;
        case VIAL_KW_OR:   parts[io] = vial_sets_or(a->sets, i, b->sets, j, k->alloc); break;
        case VIAL_KW_NOT:  parts[io] = vial_sets_not(a->sets, b->sets, j, k->alloc); break;     // [scope, a]
        default: {
            // a and not b
            vial_sets_t* one = vial_sets_pick(a->sets, &i, 1, k->alloc);
            parts[io] = vial_sets_not(one, b->sets, j, k->alloc);
        }
        }
        total += parts[io]->offset[parts[io]->count];
    });
    if (bc.nd == 0) {
        parts[0]->single = true;
        *out = vial_val_sets(parts[0]);
        return true;
    }
    vial_sets_t* s = md_alloc(k->alloc, sizeof(vial_sets_t));
    memset(s, 0, sizeof(*s));
    s->count = (uint32_t)bc.n;
    s->offset = md_alloc(k->alloc, (bc.n + 1) * sizeof(uint32_t));
    s->idx = md_alloc(k->alloc, (total ? total : 1) * sizeof(uint32_t));
    uint32_t at = 0;
    s->offset[0] = 0;
    for (size_t i = 0; i < bc.n; ++i) {
        const uint32_t sz = vial_sets_size(parts[i], 0);
        memcpy(s->idx + at, vial_sets_items(parts[i], 0), sz * sizeof(uint32_t));
        at += sz;
        s->offset[i + 1] = at;
    }
    vial_val_t o = {.dtype = VIAL_DT_SET, .ndim = bc.nd, .sets = s};
    memcpy(o.dims, bc.dims, sizeof(bc.dims));
    *out = o;
    return true;
}

// ### INDEXING ###

// What an index item selects along the dimensions it takes: the offsets it adds, one list per output dimension
typedef struct sel_t {
    uint8_t   src_dims;     // how many dimensions of a it takes
    uint8_t   out_dims;     // 0 or 1
    size_t*   offsets;      // per output position, the element offset it adds
    uint32_t  count;        // output extent (1 when out_dims is 0)
} sel_t;

static int64_t int_at(const vial_val_t* v, size_t i) { return (int64_t)get_f64(v, i); }

static bool index_set_particles(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* it, vial_val_t* out);

bool vial_k_index(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* items, uint32_t num, vial_val_t* out) {
    if (a->dtype == VIAL_DT_SET && a->ndim == 0) {
        if (num != 1) return vial_rt_fail(k, "R0103", "a set takes one position");
        return index_set_particles(k, a, items, out);
    }
    vial_interp_t* in = k->in;
    const vial_val_t* vals = in->vals;
    // Strides of a, in elements
    size_t st[VIAL_MAX_DIMS];
    {
        size_t s = 1;
        for (int d = a->ndim - 1; d >= 0; --d) { st[d] = s; s *= a->dims[d]; }
    }
    // How many dimensions the items take, so ... knows how many it stands for
    uint32_t takes = 0;
    for (uint32_t i = 0; i < num; ++i) {
        if (items[i].kind == VIAL_INDEX_NEWAXIS || items[i].kind == VIAL_INDEX_ELLIPSIS) continue;
        if (items[i].kind == VIAL_INDEX_ARRAY && vals[items[i].a].dtype == VIAL_DT_BOOL) takes += vals[items[i].a].ndim;
        else takes += 1;
    }
    if (takes > a->ndim) return vial_rt_fail(k, "R0103", "%u positions for an array of %u dimensions", takes, a->ndim);
    sel_t sels[VIAL_MAX_DIMS * 2 + 2];
    uint32_t ns = 0;
    uint8_t d = 0;
    bool missing_all = false;
    for (uint32_t i = 0; i <= num; ++i) {
        // After the items, the dimensions left are taken whole
        const vial_index_item_t* it = i < num ? &items[i] : NULL;
        if (!it) {
            while (d < a->ndim) {
                sel_t s = {.src_dims = 1, .out_dims = 1, .count = a->dims[d]};
                s.offsets = md_alloc(k->alloc, (s.count ? s.count : 1) * sizeof(size_t));
                for (uint32_t q = 0; q < s.count; ++q) s.offsets[q] = q * st[d];
                sels[ns++] = s;
                d++;
            }
            break;
        }
        if (it->kind == VIAL_INDEX_ELLIPSIS) {
            const uint32_t skip = a->ndim - takes;
            for (uint32_t q = 0; q < skip; ++q, ++d) {
                sel_t s = {.src_dims = 1, .out_dims = 1, .count = a->dims[d]};
                s.offsets = md_alloc(k->alloc, (s.count ? s.count : 1) * sizeof(size_t));
                for (uint32_t p = 0; p < s.count; ++p) s.offsets[p] = p * st[d];
                sels[ns++] = s;
            }
            continue;
        }
        if (it->kind == VIAL_INDEX_NEWAXIS) {
            sel_t s = {.src_dims = 0, .out_dims = 1, .count = 1};
            s.offsets = md_alloc(k->alloc, sizeof(size_t));
            s.offsets[0] = 0;
            sels[ns++] = s;
            continue;
        }
        const uint32_t n = a->dims[d];
        if (it->kind == VIAL_INDEX_SCALAR) {
            const vial_val_t* iv = &vals[it->a];
            int64_t p = int_at(iv, 0);
            if (p < 0) p += n;
            sel_t s = {.src_dims = 1, .out_dims = 0, .count = 1};
            s.offsets = md_alloc(k->alloc, sizeof(size_t));
            if (vial_val_is_missing(iv, 0) || p < 0 || p >= n) { missing_all = true; s.offsets[0] = 0; }    // R0101: missing
            else s.offsets[0] = (size_t)p * st[d];
            sels[ns++] = s;
            d++;
            continue;
        }
        if (it->kind == VIAL_INDEX_SLICE) {
            int64_t step = it->c ? int_at(&vals[it->c], 0) : 1;
            if (step == 0) return vial_rt_fail(k, "R0103", "a slice's step is 0");
            int64_t lo, hi;
            if (step > 0) {
                lo = it->a ? int_at(&vals[it->a], 0) : 0;
                hi = it->b ? int_at(&vals[it->b], 0) : n;
                if (lo < 0) lo = lo + n < 0 ? 0 : lo + n; else if (lo > n) lo = n;
                if (hi < 0) hi = hi + n < 0 ? 0 : hi + n; else if (hi > n) hi = n;
            } else {
                lo = it->a ? int_at(&vals[it->a], 0) : (int64_t)n - 1;
                hi = it->b ? int_at(&vals[it->b], 0) : -1 - (int64_t)n;
                if (lo < 0) lo = lo + n < -1 ? -1 : lo + n; else if (lo >= n) lo = (int64_t)n - 1;
                if (it->b) { if (hi < 0) hi = hi + n < -1 ? -1 : hi + n; else if (hi >= n) hi = (int64_t)n - 1; }
                else hi = -1;
            }
            uint32_t cnt = 0;
            if (step > 0 && hi > lo) cnt = (uint32_t)((hi - lo + step - 1) / step);
            if (step < 0 && lo > hi) cnt = (uint32_t)((lo - hi - step - 1) / (-step));
            sel_t s = {.src_dims = 1, .out_dims = 1, .count = cnt};
            s.offsets = md_alloc(k->alloc, (cnt ? cnt : 1) * sizeof(size_t));
            for (uint32_t q = 0; q < cnt; ++q) s.offsets[q] = (size_t)(lo + (int64_t)q * step) * st[d];
            sels[ns++] = s;
            d++;
            continue;
        }
        if (it->kind == VIAL_INDEX_SET) {
            // The positions along dimension c whose particle (from b's list) is in a's set, in the dimension's order
            const vial_val_t* sv = &vals[it->a];
            const vial_sets_t* over = vals[it->b].sets;
            const uint32_t pd = (uint32_t)int_at(&vals[it->c], 0);
            if (pd != d) return vial_rt_fail(k, "R0103", "a set picks along a dimension not first");
            const vial_sets_t* both = vial_sets_and(over, 0, sv->sets, 0, k->alloc);
            const uint32_t m = vial_sets_size(both, 0);
            sel_t s = {.src_dims = 1, .out_dims = 1, .count = m};
            s.offsets = md_alloc(k->alloc, (m ? m : 1) * sizeof(size_t));
            const uint32_t* oi = vial_sets_items(over, 0);
            const uint32_t on = vial_sets_size(over, 0);
            const uint32_t* bi = vial_sets_items(both, 0);
            uint32_t q = 0;
            for (uint32_t p = 0; p < on && q < m; ++p) if (oi[p] == bi[q]) s.offsets[q++] = (size_t)p * st[d];
            sels[ns++] = s;
            d++;
            continue;
        }
        // An array: ints pick positions; a bool mask over the dimensions it covers keeps where it is true
        const vial_val_t* iv = &vals[it->a];
        if (iv->dtype == VIAL_DT_BOOL) {
            const size_t mn = vial_val_count(iv);
            uint32_t cnt = 0;
            for (size_t q = 0; q < mn; ++q) cnt += ((const uint8_t*)iv->data)[q] && !vial_val_is_missing(iv, q);
            sel_t s = {.src_dims = iv->ndim, .out_dims = 1, .count = cnt};
            s.offsets = md_alloc(k->alloc, (cnt ? cnt : 1) * sizeof(size_t));
            for (uint8_t q = 0; q < iv->ndim; ++q) {
                if (iv->dims[q] != a->dims[d + q]) return vial_rt_fail(k, "R0103", "a mask of %u for a dimension of %u", iv->dims[q], a->dims[d + q]);
            }
            uint32_t w = 0;
            for (size_t q = 0; q < mn; ++q) {
                if (!((const uint8_t*)iv->data)[q] || vial_val_is_missing(iv, q)) continue;
                size_t rem = q, off = 0;
                for (int e = iv->ndim - 1; e >= 0; --e) { off += (rem % iv->dims[e]) * st[d + e]; rem /= iv->dims[e]; }
                s.offsets[w++] = off;
            }
            sels[ns++] = s;
            d += iv->ndim;
            continue;
        }
        const size_t cnt = vial_val_count(iv);
        sel_t s = {.src_dims = 1, .out_dims = 1, .count = (uint32_t)cnt};
        s.offsets = md_alloc(k->alloc, (cnt ? cnt : 1) * sizeof(size_t));
        for (size_t q = 0; q < cnt; ++q) {
            int64_t p = int_at(iv, q);
            if (p < 0) p += n;
            if (p < 0 || p >= n) return vial_rt_fail(k, "R0101", "position %lld is outside an extent of %u", (long long)int_at(iv, q), n);
            s.offsets[q] = (size_t)p * st[d];
        }
        sels[ns++] = s;
        d++;
    }
    // The result: the selections' output dimensions in order, gathered
    vial_val_t o = {.dtype = a->dtype};
    size_t total = 1;
    for (uint32_t i = 0; i < ns; ++i) {
        if (!sels[i].out_dims) continue;
        if (o.ndim >= VIAL_MAX_DIMS) return vial_rt_fail(k, "R0103", "more than 6 dimensions");
        o.dims[o.ndim++] = sels[i].count;
        if (sels[i].count && total > SIZE_MAX / sels[i].count) total = SIZE_MAX;
        else total *= sels[i].count;
    }
    if (!vial_k_room(k, total)) return false;
    size_t* src = md_alloc(k->alloc, (total ? total : 1) * sizeof(size_t));
    {
        uint32_t idx[VIAL_MAX_DIMS * 2 + 2] = {0};
        for (size_t q = 0; q < total; ++q) {
            size_t off = 0;
            for (uint32_t i = 0; i < ns; ++i) off += sels[i].offsets[sels[i].out_dims ? idx[i] : 0];
            src[q] = off;
            for (int i = (int)ns - 1; i >= 0; --i) {
                if (!sels[i].out_dims) continue;
                if (++idx[i] < sels[i].count) break;
                idx[i] = 0;
            }
        }
    }
    if (a->dtype == VIAL_DT_SET) {
        // A set's particles by position, or sets of an array picked
        uint32_t* pos = md_alloc(k->alloc, (total ? total : 1) * sizeof(uint32_t));
        for (size_t q = 0; q < total; ++q) pos[q] = (uint32_t)src[q];
        if (a->ndim == 0 && a->sets->single) {
            return vial_rt_fail(k, "R0103", "positions in a set are taken on its particles");
        }
        vial_sets_t* s = vial_sets_pick(a->sets, pos, (uint32_t)total, k->alloc);
        if (o.ndim == 0) s->single = true;
        o.sets = s;
        *out = o;
        return true;
    }
    const size_t es = vial_dtype_size(a->dtype);
    if (!es) return vial_rt_fail(k, "R0103", "indexing what holds no values");
    vial_val_t m; if (!vial_k_make(k, a->dtype, o.ndim, o.dims, &m)) return false;
    o.data = m.data;
    // A position past the end of an empty dimension reads nothing: the result is missing throughout
    if (!missing_all && a->data) for (size_t q = 0; q < total; ++q) memcpy((uint8_t*)o.data + q * es, (const uint8_t*)a->data + src[q] * es, es);
    for (size_t q = 0; q < total; ++q) {
        if (missing_all || vial_val_is_missing(a, src[q])) vial_val_set_missing(&o, q, k->alloc);
    }
    *out = o;
    return true;
}

// A set's particles, picked by position: s[i], s[2:5], s[[0, 3]]
static bool index_set_particles(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* it, vial_val_t* out) {
    const vial_val_t* vals = k->in->vals;
    const uint32_t n = vial_sets_size(a->sets, 0);
    const vial_val_t view = {.dtype = VIAL_DT_INT, .ndim = 1, .dims = {n}};
    (void)view;
    // The positions, as the index would give them on [0, n)
    int32_t* iota = md_alloc(k->alloc, (n ? n : 1) * sizeof(int32_t));
    for (uint32_t i = 0; i < n; ++i) iota[i] = (int32_t)i;
    vial_val_t iv = {.dtype = VIAL_DT_INT, .ndim = 1, .dims = {n}, .data = iota};
    vial_val_t pos;
    vial_kctx_t k2 = *k;
    if (!vial_k_index(&k2, &iv, it, 1, &pos)) return false;
    (void)vals;
    const size_t cnt = vial_val_count(&pos);
    uint32_t* p = md_alloc(k->alloc, (cnt ? cnt : 1) * sizeof(uint32_t));
    for (size_t i = 0; i < cnt; ++i) {
        if (vial_val_is_missing(&pos, i)) {
            // A position past the end: an empty set, which has no position (R0101)
            vial_sets_t* e = vial_sets_one(p, 0, k->alloc);
            e->single = true;
            *out = vial_val_sets(e);
            return true;
        }
        p[i] = (uint32_t)((const int32_t*)pos.data)[i];
    }
    vial_sets_t* s = vial_sets_pick_particles(a->sets, 0, p, (uint32_t)cnt, k->alloc);
    s->single = true;
    *out = vial_val_sets(s);
    return true;
}

// ### STACKING ###

bool vial_k_stack(vial_kctx_t* k, const vial_val_t* items, uint32_t n, vial_val_t* out) {
    const uint8_t dt = k->type->dtype;
    vial_val_t o = {.dtype = dt, .ndim = 1, .dims = {n}};
    if (n) {
        if (items[0].ndim + 1 > VIAL_MAX_DIMS) return vial_rt_fail(k, "R0103", "more than 6 dimensions");
        for (uint8_t d = 0; d < items[0].ndim; ++d) o.dims[o.ndim++] = items[0].dims[d];
        for (uint32_t i = 1; i < n; ++i) {
            if (items[i].ndim != items[0].ndim) return vial_rt_fail(k, "R0103", "the items of a list differ in shape");
            for (uint8_t d = 0; d < items[0].ndim; ++d) {
                if (items[i].dims[d] != items[0].dims[d]) return vial_rt_fail(k, "R0103", "the items of a list differ in shape");
            }
        }
    } else {
        // An empty list: the type's shape, with its first extent 0
        o.ndim = k->type->ndim ? k->type->ndim : 1;
        for (uint8_t d = 0; d < o.ndim; ++d) o.dims[d] = d == 0 ? 0 : (k->type->dims[d].kind == VIAL_EXT_STATIC ? k->type->dims[d].n : 0);
    }
    if (dt == VIAL_DT_SET) {
        uint32_t total = 0, count = 0;
        for (uint32_t i = 0; i < n; ++i) { total += items[i].sets->offset[items[i].sets->count]; count += items[i].sets->count; }
        vial_sets_t* s = md_alloc(k->alloc, sizeof(vial_sets_t));
        memset(s, 0, sizeof(*s));
        s->count = count;
        s->offset = md_alloc(k->alloc, (count + 1) * sizeof(uint32_t));
        s->idx = md_alloc(k->alloc, (total ? total : 1) * sizeof(uint32_t));
        uint32_t at = 0, c = 0;
        s->offset[0] = 0;
        for (uint32_t i = 0; i < n; ++i) {
            const vial_sets_t* e = items[i].sets;
            for (uint32_t q = 0; q < e->count; ++q) {
                const uint32_t sz = vial_sets_size(e, q);
                memcpy(s->idx + at, vial_sets_items(e, q), sz * sizeof(uint32_t));
                at += sz;
                s->offset[++c] = at;
            }
        }
        o.sets = s;
        *out = o;
        return true;
    }
    const size_t es = vial_dtype_size(dt);
    if (!es) { *out = o; return true; }
    vial_val_t m; if (!vial_k_make(k, dt, o.ndim, o.dims, &m)) return false;
    o.data = m.data;
    size_t inner = 1;
    for (uint8_t d = 1; d < o.ndim; ++d) inner *= o.dims[d];
    for (uint32_t i = 0; i < n; ++i) {
        const vial_val_t c = vial_val_as(&items[i], dt, k->alloc);
        if (c.data) memcpy((uint8_t*)o.data + i * inner * es, c.data, inner * es);
        for (size_t q = 0; q < inner; ++q) if (vial_val_is_missing(&c, q)) vial_val_set_missing(&o, i * inner + q, k->alloc);
    }
    *out = o;
    return true;
}

// ### KINDS ###

bool vial_k_positions(vial_kctx_t* k, const vial_val_t* s, vial_val_t* out);

bool vial_k_kind(vial_kctx_t* k, uint16_t kind, const vial_val_t* a, vial_val_t* out) {
    if (kind == VIAL_KIND_EACH) {
        // A set's particles, each a set of one
        if (a->dtype != VIAL_DT_SET) return vial_rt_fail(k, "R0103", "each particle of what is not a set");
        const uint32_t n = vial_sets_size(a->sets, 0);
        vial_sets_t* s = md_alloc(k->alloc, sizeof(vial_sets_t));
        memset(s, 0, sizeof(*s));
        s->count = n;
        s->offset = md_alloc(k->alloc, (n + 1) * sizeof(uint32_t));
        s->idx = md_alloc(k->alloc, (n ? n : 1) * sizeof(uint32_t));
        for (uint32_t i = 0; i <= n; ++i) s->offset[i] = i;
        if (n) memcpy(s->idx, vial_sets_items(a->sets, 0), n * sizeof(uint32_t));
        *out = vial_val_sets(s);
        return true;
    }
    if (a->dtype != VIAL_DT_SET) {
        // A point as a group of one: [3] as [1][3]
        vial_val_t o = *a;
        if (a->ndim == 1 && (kind == VX_KIND_POINTS || kind == VX_KIND_GROUP)) { o.ndim = 2; o.dims[0] = 1; o.dims[1] = a->dims[0]; }
        *out = o;
        return true;
    }
    vial_val_t p;
    if (!vial_k_positions(k, a, &p)) return false;
    if (kind == VX_KIND_POINT) {
        // One particle's position; an empty set has none
        vial_val_t o; if (!vial_k_make(k, VIAL_DT_FLOAT, 1, (uint32_t[]){3}, &o)) return false;
        if (p.dims[0] >= 1) memcpy(o.data, p.data, 3 * sizeof(float));
        else for (int i = 0; i < 3; ++i) vial_val_set_missing(&o, i, k->alloc);
        *out = o;
        return true;
    }
    *out = p;
    return true;
}

// Index helper exposed to the built-ins: a single set's particles by an index item
bool vial_k_index_set(vial_kctx_t* k, const vial_val_t* a, const vial_index_item_t* it, vial_val_t* out) {
    return index_set_particles(k, a, it, out);
}

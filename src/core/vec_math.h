#ifndef _VEC_MATH_H_
#define _VEC_MATH_H_

#ifndef VEC_MATH_USE_SSE_H
#define VEC_MATH_USE_SSE_H 1
#endif

#ifdef VEC_MATH_USE_SSE_H
#include <xmmintrin.h>
#endif

typedef union md_vec4 {
    struct {
        float x, y, z, w;
    };
    float elem[4];
#if VEC_MATH_USE_SSE_H
    __m128 mm128;
#endif
} md_vec4;

typedef union md_mat4 {
    float elem[4][4];
    md_vec4 col[4];
} md_mat4;

inline md_vec4 md_vec4_mul(const md_vec4 a, const md_vec4 b) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_mul_ps(a.mm128, b.mm128);
#else
    c.x = a.x * b.x;
    c.y = a.y * b.y;
    c.z = a.z * b.z;
    c.w = a.w * b.w;
#endif
    return c;
}

inline md_vec4 md_vec4_mul_f(const md_vec4 a, const float s) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_mul_ps(a.mm128, _mm_set1_ps(s));
#else
    c.x = a.x * s;
    c.y = a.y * s;
    c.z = a.z * s;
    c.w = a.w * s;
#endif
    return c;
}

inline md_vec4 md_vec4_div(const md_vec4 a, const md_vec4 b) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_div_ps(a.mm128, b.mm128);
#else
    c.x = a.x / b.x;
    c.y = a.y / b.y;
    c.z = a.z / b.z;
    c.w = a.w / b.w;
#endif
    return c;
}

inline md_vec4 vec4_div_f(const md_vec4 a, const float s) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_div_ps(a.mm128, _mm_set1_ps(s));
#else
    c.x = a.x / s;
    c.y = a.y / s;
    c.z = a.z / s;
    c.w = a.w / s;
#endif
    return c;
}

inline md_vec4 md_vec4_add(const md_vec4 a, const md_vec4 b) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_add_ps(a.mm128, b.mm128);
#else
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    c.w = a.w + b.w;
#endif
    return c;
}

inline md_vec4 md_vec4_add_f(const md_vec4 a, const float s) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_add_ps(a.mm128, _mm_set1_ps(s));
#else
    c.x = a.x + s;
    c.y = a.y + s;
    c.z = a.z + s;
    c.w = a.w + s;
#endif
    return c;
}

inline md_vec4 md_vec4_sub(const md_vec4 a, const md_vec4 b) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_sub_ps(a.mm128, b.mm128);
#else
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    c.w = a.w - b.w;
#endif
    return c;
}

inline md_vec4 md_vec4_sub_f(const md_vec4 a, const float s) {
    md_vec4 c;
#if VEC_MATH_USE_SSE_H
    c.mm128 = _mm_sub_ps(a.mm128, _mm_set1_ps(s));
#else
    c.x = a.x - s;
    c.y = a.y - s;
    c.z = a.z - s;
    c.w = a.w - s;
#endif
    return c;
}

#if VEC_MATH_USE_SSE_H
inline __m128 md_linear_combine_sse(__m128 a, md_mat4 B) {
    __m128 res;
    res = _mm_mul_ps(_mm_shuffle_ps(a, a, 0x00), B.col[0].mm128);
    res = _mm_add_ps(res, _mm_mul_ps(_mm_shuffle_ps(a, a, 0x55), B.col[1].mm128));
    res = _mm_add_ps(res, _mm_mul_ps(_mm_shuffle_ps(a, a, 0xaa), B.col[2].mm128));
    res = _mm_add_ps(res, _mm_mul_ps(_mm_shuffle_ps(a, a, 0xff), B.col[3].mm128));
    return res;
}
#endif

inline md_mat4 md_mat4_mul(const md_mat4 A, const md_mat4 B) {
    md_mat4 C;
#if VEC_MATH_USE_SSE_H
    C.col[0].mm128 = md_linear_combine_sse(B.col[0].mm128, A);
    C.col[1].mm128 = md_linear_combine_sse(B.col[1].mm128, A);
    C.col[2].mm128 = md_linear_combine_sse(B.col[2].mm128, A);
    C.col[3].mm128 = md_linear_combine_sse(B.col[3].mm128, A);
#else
#define MULT(col, row) \
    A.elem[0][row] * B.elem[col][0] + A.elem[1][row] * B.elem[col][1] + A.elem[2][row] * B.elem[col][2] + A.elem[3][row] * B.elem[col][3]
    C.elem[0][0] = MULT(0, 0);
    C.elem[0][1] = MULT(0, 1);
    C.elem[0][2] = MULT(0, 2);
    C.elem[0][3] = MULT(0, 3);

    C.elem[1][0] = MULT(1, 0);
    C.elem[1][1] = MULT(1, 1);
    C.elem[1][2] = MULT(1, 2);
    C.elem[1][3] = MULT(1, 3);

    C.elem[2][0] = MULT(2, 0);
    C.elem[2][1] = MULT(2, 1);
    C.elem[2][2] = MULT(2, 2);
    C.elem[2][3] = MULT(2, 3);

    C.elem[3][0] = MULT(3, 0);
    C.elem[3][1] = MULT(3, 1);
    C.elem[3][2] = MULT(3, 2);
    C.elem[3][3] = MULT(3, 3);
#undef MULT
#endif
    return C;
}

inline md_mat4 md_mat4_mul_f(const md_mat4 M, float s) {
    md_mat4 C;
#if VEC_MATH_USE_SSE_H
    C.col[0] = md_vec4_mul_f(M.col[0], s);
    C.col[1] = md_vec4_mul_f(M.col[1], s);
    C.col[2] = md_vec4_mul_f(M.col[2], s);
    C.col[3] = md_vec4_mul_f(M.col[3], s);
#else
    C.elem[0][0] = M.elem[0][0] * s;
    C.elem[0][1] = M.elem[0][1] * s;
    C.elem[0][2] = M.elem[0][2] * s;
    C.elem[0][3] = M.elem[0][3] * s;

    C.elem[1][0] = M.elem[1][0] * s;
    C.elem[1][1] = M.elem[1][1] * s;
    C.elem[1][2] = M.elem[1][2] * s;
    C.elem[1][3] = M.elem[1][3] * s;

    C.elem[2][0] = M.elem[2][0] * s;
    C.elem[2][1] = M.elem[2][1] * s;
    C.elem[2][2] = M.elem[2][2] * s;
    C.elem[2][3] = M.elem[2][3] * s;

    C.elem[3][0] = M.elem[3][0] * s;
    C.elem[3][1] = M.elem[3][1] * s;
    C.elem[3][2] = M.elem[3][2] * s;
    C.elem[3][3] = M.elem[3][3] * s;
#endif
    return C;
}

inline md_mat4 md_mat4_ident() {
    md_mat4 M = {0};
    M.elem[0][0] = 1.0f;
    M.elem[1][1] = 1.0f;
    M.elem[2][2] = 1.0f;
    M.elem[3][3] = 1.0f;
    return M;
}

inline md_mat4 md_mat4_inverse(const md_mat4 M) {
    const float c00 = M.elem[2][2] * M.elem[3][3] - M.elem[3][2] * M.elem[2][3];
    const float c02 = M.elem[1][2] * M.elem[3][3] - M.elem[3][2] * M.elem[1][3];
    const float c03 = M.elem[1][2] * M.elem[2][3] - M.elem[2][2] * M.elem[1][3];

    const float c04 = M.elem[2][1] * M.elem[3][3] - M.elem[3][1] * M.elem[2][3];
    const float c06 = M.elem[1][1] * M.elem[3][3] - M.elem[3][1] * M.elem[1][3];
    const float c07 = M.elem[1][1] * M.elem[2][3] - M.elem[2][1] * M.elem[1][3];

    const float c08 = M.elem[2][1] * M.elem[3][2] - M.elem[3][1] * M.elem[2][2];
    const float c10 = M.elem[1][1] * M.elem[3][2] - M.elem[3][1] * M.elem[1][2];
    const float c11 = M.elem[1][1] * M.elem[2][2] - M.elem[2][1] * M.elem[1][2];

    const float c12 = M.elem[2][0] * M.elem[3][3] - M.elem[3][0] * M.elem[2][3];
    const float c14 = M.elem[1][0] * M.elem[3][3] - M.elem[3][0] * M.elem[1][3];
    const float c15 = M.elem[1][0] * M.elem[2][3] - M.elem[2][0] * M.elem[1][3];

    const float c16 = M.elem[2][0] * M.elem[3][2] - M.elem[3][0] * M.elem[2][2];
    const float c18 = M.elem[1][0] * M.elem[3][2] - M.elem[3][0] * M.elem[1][2];
    const float c19 = M.elem[1][0] * M.elem[2][2] - M.elem[2][0] * M.elem[1][2];

    const float c20 = M.elem[2][0] * M.elem[3][1] - M.elem[3][0] * M.elem[2][1];
    const float c22 = M.elem[1][0] * M.elem[3][1] - M.elem[3][0] * M.elem[1][1];
    const float c23 = M.elem[1][0] * M.elem[2][1] - M.elem[2][0] * M.elem[1][1];

    const md_vec4 f0 = {c00, c00, c02, c03};
    const md_vec4 f1 = {c04, c04, c06, c07};
    const md_vec4 f2 = {c08, c08, c10, c11};
    const md_vec4 f3 = {c12, c12, c14, c15};
    const md_vec4 f4 = {c16, c16, c18, c19};
    const md_vec4 f5 = {c20, c20, c22, c23};

    const md_vec4 v0 = {M.elem[1][0], M.elem[0][0], M.elem[0][0], M.elem[0][0]};
    const md_vec4 v1 = {M.elem[1][1], M.elem[0][1], M.elem[0][1], M.elem[0][1]};
    const md_vec4 v2 = {M.elem[1][2], M.elem[0][2], M.elem[0][2], M.elem[0][2]};
    const md_vec4 v3 = {M.elem[1][3], M.elem[0][3], M.elem[0][3], M.elem[0][3]};

    const md_vec4 i0 = md_vec4_add(md_vec4_sub(md_vec4_mul(v1, f0), md_vec4_mul(v2, f1)), md_vec4_mul(v3, f2));
    const md_vec4 i1 = md_vec4_add(md_vec4_sub(md_vec4_mul(v0, f0), md_vec4_mul(v2, f3)), md_vec4_mul(v3, f4));
    const md_vec4 i2 = md_vec4_add(md_vec4_sub(md_vec4_mul(v0, f1), md_vec4_mul(v1, f3)), md_vec4_mul(v3, f5));
    const md_vec4 i3 = md_vec4_add(md_vec4_sub(md_vec4_mul(v0, f2), md_vec4_mul(v1, f4)), md_vec4_mul(v2, f5));

    const md_vec4 sign_a = {+1, -1, +1, -1};
    const md_vec4 sign_b = {-1, +1, -1, +1};

    md_mat4 I;
    I.col[0] = md_vec4_mul(i0, sign_a);
    I.col[1] = md_vec4_mul(i1, sign_b);
    I.col[2] = md_vec4_mul(i2, sign_a);
    I.col[3] = md_vec4_mul(i3, sign_b);

    const md_vec4 row0 = {I.elem[0][0], I.elem[1][0], I.elem[2][0], I.elem[3][0]};
    const md_vec4 dot0 = md_vec4_mul(M.col[0], row0);

    return md_mat4_mul_f(I, 1.0f / (dot0.x + dot0.y + dot0.z + dot0.w));
}

inline md_mat4 md_mat4_transpose(const md_mat4 M) {
    md_mat4 T;

#if VEC_MATH_USE_SSE_H
    T = M;
    _MM_TRANSPOSE4_PS(T.col[0].mm128, T.col[1].mm128, T.col[2].mm128, T.col[3].mm128);
#else
    T.elem[0][0] = M.elem[0][0];
    T.elem[0][1] = M.elem[1][0];
    T.elem[0][2] = M.elem[2][0];
    T.elem[0][3] = M.elem[3][0];

    T.elem[1][0] = M.elem[0][1];
    T.elem[1][1] = M.elem[1][1];
    T.elem[1][2] = M.elem[2][1];
    T.elem[1][3] = M.elem[3][1];

    T.elem[2][0] = M.elem[0][2];
    T.elem[2][1] = M.elem[1][2];
    T.elem[2][2] = M.elem[2][2];
    T.elem[2][3] = M.elem[3][2];

    T.elem[3][0] = M.elem[0][3];
    T.elem[3][1] = M.elem[1][3];
    T.elem[3][2] = M.elem[2][3];
    T.elem[3][3] = M.elem[3][3];
#endif
    return T;
}


#endif
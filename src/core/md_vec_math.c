#include "md_vec_math.h"

#include "ext/svd3/svd3.h"
#include <math.h>

#define SWAP(x, y) {int t = x; x = y; y = t;}

void mat3_svd(const mat3_t M, mat3_t* U, mat3_t* S, mat3_t* V) {
    mat3_t Mt = mat3_transpose(M);
    // the external svd library uses row major matrix convention...
    svd(Mt.elem, U->elem, S->elem, V->elem);
    *U = mat3_transpose(*U);
    *S = mat3_transpose(*S);
    *V = mat3_transpose(*V);
}

quat_t quat_angle_axis(float angle, vec3_t axis) {
    float half_angle = angle * 0.5f;
    float sin_angle = sinf(half_angle);

    quat_t q;
    q.x = axis.x * sin_angle;
    q.y = axis.y * sin_angle;
    q.z = axis.z * sin_angle;
    q.w = cosf(half_angle);
    return q;
}

quat_t quat_normalize(quat_t q) {
    quat_t r = q;

    float mag2 = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
    if (mag2 != 0.0f && (fabsf(mag2 - 1.0f) > 0.000001)) {
        float mag = 1.0f / sqrtf(mag2);
        r.x = q.x * mag;
        r.y = q.y * mag;
        r.z = q.z * mag;
        r.w = q.w * mag;
    }

    return r;
}

quat_t quat_from_mat3(mat3_t M) {
    quat_t q;
    q.w = sqrtf( MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2]) ) * 0.5f;
    q.x = sqrtf( MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2]) ) * 0.5f;
    q.y = sqrtf( MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2]) ) * 0.5f;
    q.z = sqrtf( MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2]) ) * 0.5f;

    q.x = copysignf( q.x, M.elem[2][1] - M.elem[1][2] );
    q.y = copysignf( q.y, M.elem[0][2] - M.elem[2][0] );
    q.z = copysignf( q.z, M.elem[1][0] - M.elem[0][1] );
    return q;
}

quat_t quat_from_mat4(mat4_t M) {
    quat_t q;
    q.w = sqrtf( MAX(0, 1 + M.elem[0][0] + M.elem[1][1] + M.elem[2][2]) ) * 0.5f;
    q.x = sqrtf( MAX(0, 1 + M.elem[0][0] - M.elem[1][1] - M.elem[2][2]) ) * 0.5f;
    q.y = sqrtf( MAX(0, 1 - M.elem[0][0] + M.elem[1][1] - M.elem[2][2]) ) * 0.5f;
    q.z = sqrtf( MAX(0, 1 - M.elem[0][0] - M.elem[1][1] + M.elem[2][2]) ) * 0.5f;

    q.x = copysignf( q.x, M.elem[2][1] - M.elem[1][2] );
    q.y = copysignf( q.y, M.elem[0][2] - M.elem[2][0] );
    q.z = copysignf( q.z, M.elem[1][0] - M.elem[0][1] );
    return q;
}

void mat3_eigen(mat3_t M, vec3_t vectors[3], float values[3]) {
    mat3_t U, S, V;
    mat3_svd(M, &U, &S, &V);

    const float max_val = MAX(S.elem[0][0], MAX(S.elem[1][1], S.elem[2][2]));
    const float  e_val[] = {S.elem[0][0] / max_val, S.elem[1][1] / max_val, S.elem[2][2] / max_val};
    const vec3_t e_vec[] = {U.col[0], U.col[1], U.col[2]};

    int l[3] = {0, 1, 2};

    if (e_val[l[0]] < e_val[l[1]]) SWAP(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) SWAP(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) SWAP(l[0], l[1]);

    values[0] = e_val[l[0]];
    values[1] = e_val[l[1]];
    values[2] = e_val[l[2]];

    vectors[0] = e_vec[l[0]];
    vectors[1] = e_vec[l[1]];
    vectors[2] = e_vec[l[2]];
}

mat3_t mat3_covariance_matrix( const float* x, const float* y, const float* z, vec3_t com, int64_t count) {
    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const float px = x[i] - com.x;
        const float py = y[i] - com.y;
        const float pz = z[i] - com.z;

        A.elem[0][0] += px * px;
        A.elem[0][1] += px * py;
        A.elem[0][2] += px * pz;
        A.elem[1][0] += py * px;
        A.elem[1][1] += py * py;
        A.elem[1][2] += py * pz;
        A.elem[2][0] += pz * px;
        A.elem[2][1] += pz * py;
        A.elem[2][2] += pz * pz;
    }

    return A;
}

mat3_t mat3_covariance_matrix_vec3(const vec3_t* xyz, vec3_t com, int64_t count) {
    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const float px = xyz[i].x - com.x;
        const float py = xyz[i].y - com.y;
        const float pz = xyz[i].z - com.z;

        A.elem[0][0] += px * px;
        A.elem[0][1] += px * py;
        A.elem[0][2] += px * pz;
        A.elem[1][0] += py * px;
        A.elem[1][1] += py * py;
        A.elem[1][2] += py * pz;
        A.elem[2][0] += pz * px;
        A.elem[2][1] += pz * py;
        A.elem[2][2] += pz * pz;
    }

    return A;
}

mat3_t mat3_weighted_covariance_matrix(const float* x, const float* y, const float* z, const float* weight, vec3_t com, int64_t count) {
    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const float px = x[i] - com.x;
        const float py = y[i] - com.y;
        const float pz = z[i] - com.z;
        const float w = weight[i];

        A.elem[0][0] += w * px * px;
        A.elem[0][1] += w * px * py;
        A.elem[0][2] += w * px * pz;
        A.elem[1][0] += w * py * px;
        A.elem[1][1] += w * py * py;
        A.elem[1][2] += w * py * pz;
        A.elem[2][0] += w * pz * px;
        A.elem[2][1] += w * pz * py;
        A.elem[2][2] += w * pz * pz;
    }

    return A;
}

mat3_t mat3_cross_covariance_matrix(
    const float* x0, const float* y0, const float* z0,
    const float* x1, const float* y1, const float* z1,
    vec3_t com0,
    vec3_t com1,
    int64_t count) {

    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const float px = x0[i] - com0.x;
        const float py = y0[i] - com0.y;
        const float pz = z0[i] - com0.z;

        const float qx = x1[i] - com1.x;
        const float qy = y1[i] - com1.y;
        const float qz = z1[i] - com1.z;

        A.elem[0][0] += px * qx;
        A.elem[0][1] += px * qy;
        A.elem[0][2] += px * qz;
        A.elem[1][0] += py * qx;
        A.elem[1][1] += py * qy;
        A.elem[1][2] += py * qz;
        A.elem[2][0] += pz * qx;
        A.elem[2][1] += pz * qy;
        A.elem[2][2] += pz * qz;
    }

    return A;
}

mat3_t mat3_weighted_cross_covariance_matrix(
    const float* x0, const float* y0, const float* z0,
    const float* x1, const float* y1, const float* z1,
    const float* weight,
    vec3_t com0,
    vec3_t com1,
    int64_t count) {

    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const float px = x0[i] - com0.x;
        const float py = y0[i] - com0.y;
        const float pz = z0[i] - com0.z;

        const float qx = x1[i] - com1.x;
        const float qy = y1[i] - com1.y;
        const float qz = z1[i] - com1.z;

        const float w = weight[i];

        A.elem[0][0] += w * px * qx;
        A.elem[0][1] += w * px * qy;
        A.elem[0][2] += w * px * qz;
        A.elem[1][0] += w * py * qx;
        A.elem[1][1] += w * py * qy;
        A.elem[1][2] += w * py * qz;
        A.elem[2][0] += w * pz * qx;
        A.elem[2][1] += w * pz * qy;
        A.elem[2][2] += w * pz * qz;
    }

    return A;
}

mat3_t mat3_extract_rotation(mat3_t M) {
    mat3_t U, S, V;
    mat3_svd(M, &U, &S, &V);

    mat3_t Ut = mat3_transpose(U);
    float  d = mat3_determinant(mat3_mul(V, Ut));
    mat3_t D = {1, 0, 0, 0, 1, 0, 0, 0, d};
    mat3_t R = mat3_mul(mat3_mul(V, D), Ut);
    return R;
}

mat4_t mat4_inverse(mat4_t M) {
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

    const vec4_t f0 = {c00, c00, c02, c03};
    const vec4_t f1 = {c04, c04, c06, c07};
    const vec4_t f2 = {c08, c08, c10, c11};
    const vec4_t f3 = {c12, c12, c14, c15};
    const vec4_t f4 = {c16, c16, c18, c19};
    const vec4_t f5 = {c20, c20, c22, c23};

    const vec4_t v0 = {M.elem[1][0], M.elem[0][0], M.elem[0][0], M.elem[0][0]};
    const vec4_t v1 = {M.elem[1][1], M.elem[0][1], M.elem[0][1], M.elem[0][1]};
    const vec4_t v2 = {M.elem[1][2], M.elem[0][2], M.elem[0][2], M.elem[0][2]};
    const vec4_t v3 = {M.elem[1][3], M.elem[0][3], M.elem[0][3], M.elem[0][3]};

    const vec4_t i0 = vec4_add(vec4_sub(vec4_mul(v1, f0), vec4_mul(v2, f1)), vec4_mul(v3, f2));
    const vec4_t i1 = vec4_add(vec4_sub(vec4_mul(v0, f0), vec4_mul(v2, f3)), vec4_mul(v3, f4));
    const vec4_t i2 = vec4_add(vec4_sub(vec4_mul(v0, f1), vec4_mul(v1, f3)), vec4_mul(v3, f5));
    const vec4_t i3 = vec4_add(vec4_sub(vec4_mul(v0, f2), vec4_mul(v1, f4)), vec4_mul(v2, f5));

    const vec4_t sign_a = {+1, -1, +1, -1};
    const vec4_t sign_b = {-1, +1, -1, +1};

    mat4_t I = {0};
    I.col[0] = vec4_mul(i0, sign_a);
    I.col[1] = vec4_mul(i1, sign_b);
    I.col[2] = vec4_mul(i2, sign_a);
    I.col[3] = vec4_mul(i3, sign_b);

    const vec4_t row0 = {I.elem[0][0], I.elem[1][0], I.elem[2][0], I.elem[3][0]};
    const vec4_t dot0 = vec4_mul(M.col[0], row0);

    return mat4_mul_f(I, 1.0f / (dot0.x + dot0.y + dot0.z + dot0.w));
}


void vec3_batch_translate(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64_t count, vec3_t translation) {
    int64_t i = 0;

    const int64_t simd_count = ROUND_DOWN(count, md_simd_widthf);
    if (simd_count > 0) {
        md_simd_typef t_x = md_simd_set1f(translation.x);
        md_simd_typef t_y = md_simd_set1f(translation.y);
        md_simd_typef t_z = md_simd_set1f(translation.z);

        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef x = md_simd_loadf(in_out_x + i);
            md_simd_typef y = md_simd_loadf(in_out_y + i);
            md_simd_typef z = md_simd_loadf(in_out_z + i);

            x = md_simd_addf(x, t_x);
            y = md_simd_addf(y, t_y);
            z = md_simd_addf(z, t_z);

            md_simd_storef(in_out_x + i, x);
            md_simd_storef(in_out_y + i, y);
            md_simd_storef(in_out_z + i, z);
        }
    }

    for (; i < count; i++) {
        in_out_x[i] += translation.x;
        in_out_y[i] += translation.y;
        in_out_z[i] += translation.z;
    }
}

void vec3_batch_translate2(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, int64_t count, vec3_t translation) {
    int64_t i = 0;

    const int64_t simd_count = ROUND_DOWN(count, md_simd_widthf);
    if (simd_count > 0) {
        md_simd_typef t_x = md_simd_set1f(translation.x);
        md_simd_typef t_y = md_simd_set1f(translation.y);
        md_simd_typef t_z = md_simd_set1f(translation.z);

        for (; i < simd_count; i += md_simd_widthf) {
            md_simd_typef p_x = md_simd_loadf(in_x + i);
            md_simd_typef p_y = md_simd_loadf(in_y + i);
            md_simd_typef p_z = md_simd_loadf(in_z + i);

            p_x = md_simd_addf(p_x, t_x);
            p_y = md_simd_addf(p_y, t_y);
            p_z = md_simd_addf(p_z, t_z);

            md_simd_storef(out_x + i, p_x);
            md_simd_storef(out_y + i, p_y);
            md_simd_storef(out_z + i, p_z);
        }
    }

    for (; i < count; i++) {
        out_x[i] += translation.x;
        out_y[i] += translation.y;
        out_z[i] += translation.z;
    }
}

void mat4_batch_transform(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, float w_comp, int64_t count, mat4_t M) {
    const md_simd_typef m11 = md_simd_set1f(M.elem[0][0]);
    const md_simd_typef m12 = md_simd_set1f(M.elem[0][1]);
    const md_simd_typef m13 = md_simd_set1f(M.elem[0][2]);

    const md_simd_typef m21 = md_simd_set1f(M.elem[1][0]);
    const md_simd_typef m22 = md_simd_set1f(M.elem[1][1]);
    const md_simd_typef m23 = md_simd_set1f(M.elem[1][2]);

    const md_simd_typef m31 = md_simd_set1f(M.elem[2][0]);
    const md_simd_typef m32 = md_simd_set1f(M.elem[2][1]);
    const md_simd_typef m33 = md_simd_set1f(M.elem[2][2]);

    const md_simd_typef m41 = md_simd_set1f(M.elem[3][0]);
    const md_simd_typef m42 = md_simd_set1f(M.elem[3][1]);
    const md_simd_typef m43 = md_simd_set1f(M.elem[3][2]);

    const md_simd_typef w = md_simd_set1f(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_widthf);
    for (; i < simd_count; i += md_simd_widthf) {
        const md_simd_typef x = md_simd_loadf(in_out_x + i);
        const md_simd_typef y = md_simd_loadf(in_out_y + i);
        const md_simd_typef z = md_simd_loadf(in_out_z + i);

        const md_simd_typef m11x = md_simd_mulf(m11, x);
        const md_simd_typef m21y = md_simd_mulf(m21, y);
        const md_simd_typef m31z = md_simd_mulf(m31, z);
        const md_simd_typef m41w = md_simd_mulf(m41, w);

        const md_simd_typef m12x = md_simd_mulf(m12, x);
        const md_simd_typef m22y = md_simd_mulf(m22, y);
        const md_simd_typef m32z = md_simd_mulf(m32, z);
        const md_simd_typef m42w = md_simd_mulf(m42, w);

        const md_simd_typef m13x = md_simd_mulf(m13, x);
        const md_simd_typef m23y = md_simd_mulf(m23, y);
        const md_simd_typef m33z = md_simd_mulf(m33, z);
        const md_simd_typef m43w = md_simd_mulf(m43, w);

        const md_simd_typef res_x = md_simd_addf(md_simd_addf(m11x, m21y), md_simd_addf(m31z, m41w));
        const md_simd_typef res_y = md_simd_addf(md_simd_addf(m12x, m22y), md_simd_addf(m32z, m42w));
        const md_simd_typef res_z = md_simd_addf(md_simd_addf(m13x, m23y), md_simd_addf(m33z, m43w));

        md_simd_storef(in_out_x + i, res_x);
        md_simd_storef(in_out_y + i, res_y);
        md_simd_storef(in_out_z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_out_x[i];
        const float y = in_out_y[i];
        const float z = in_out_z[i];

        in_out_x[i] = x * M.elem[0][0] + y * M.elem[1][0] + z * M.elem[2][0] + w_comp * M.elem[3][0];
        in_out_y[i] = x * M.elem[0][1] + y * M.elem[1][1] + z * M.elem[2][1] + w_comp * M.elem[3][1];
        in_out_z[i] = x * M.elem[0][2] + y * M.elem[1][2] + z * M.elem[2][2] + w_comp * M.elem[3][2];
    }
}

void mat4_batch_transform2(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, float w_comp, int64_t count, mat4_t M) {
    const md_simd_typef m11 = md_simd_set1f(M.elem[0][0]);
    const md_simd_typef m12 = md_simd_set1f(M.elem[0][1]);
    const md_simd_typef m13 = md_simd_set1f(M.elem[0][2]);

    const md_simd_typef m21 = md_simd_set1f(M.elem[1][0]);
    const md_simd_typef m22 = md_simd_set1f(M.elem[1][1]);
    const md_simd_typef m23 = md_simd_set1f(M.elem[1][2]);

    const md_simd_typef m31 = md_simd_set1f(M.elem[2][0]);
    const md_simd_typef m32 = md_simd_set1f(M.elem[2][1]);
    const md_simd_typef m33 = md_simd_set1f(M.elem[2][2]);

    const md_simd_typef m41 = md_simd_set1f(M.elem[3][0]);
    const md_simd_typef m42 = md_simd_set1f(M.elem[3][1]);
    const md_simd_typef m43 = md_simd_set1f(M.elem[3][2]);

    const md_simd_typef w = md_simd_set1f(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_widthf);
    for (; i < simd_count; i += md_simd_widthf) {
        const md_simd_typef x = md_simd_loadf(in_x + i);
        const md_simd_typef y = md_simd_loadf(in_y + i);
        const md_simd_typef z = md_simd_loadf(in_z + i);

        const md_simd_typef m11x = md_simd_mulf(m11, x);
        const md_simd_typef m21y = md_simd_mulf(m21, y);
        const md_simd_typef m31z = md_simd_mulf(m31, z);
        const md_simd_typef m41w = md_simd_mulf(m41, w);

        const md_simd_typef m12x = md_simd_mulf(m12, x);
        const md_simd_typef m22y = md_simd_mulf(m22, y);
        const md_simd_typef m32z = md_simd_mulf(m32, z);
        const md_simd_typef m42w = md_simd_mulf(m42, w);

        const md_simd_typef m13x = md_simd_mulf(m13, x);
        const md_simd_typef m23y = md_simd_mulf(m23, y);
        const md_simd_typef m33z = md_simd_mulf(m33, z);
        const md_simd_typef m43w = md_simd_mulf(m43, w);

        const md_simd_typef res_x = md_simd_addf(md_simd_addf(m11x, m21y), md_simd_addf(m31z, m41w));
        const md_simd_typef res_y = md_simd_addf(md_simd_addf(m12x, m22y), md_simd_addf(m32z, m42w));
        const md_simd_typef res_z = md_simd_addf(md_simd_addf(m13x, m23y), md_simd_addf(m33z, m43w));

        md_simd_storef(out_x + i, res_x);
        md_simd_storef(out_y + i, res_y);
        md_simd_storef(out_z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_x[i];
        const float y = in_y[i];
        const float z = in_z[i];

        out_x[i] = x * M.elem[0][0] + y * M.elem[1][0] + z * M.elem[2][0] + w_comp * M.elem[3][0];
        out_y[i] = x * M.elem[0][1] + y * M.elem[1][1] + z * M.elem[2][1] + w_comp * M.elem[3][1];
        out_z[i] = x * M.elem[0][2] + y * M.elem[1][2] + z * M.elem[2][2] + w_comp * M.elem[3][2];
    }
}

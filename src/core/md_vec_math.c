#include "md_vec_math.h"

#include "ext/svd3/svd3.h"
#include <math.h>

#define SWAP(x, y) {int t = x; x = y; y = t;}

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
    svd(M.elem, U.elem, S.elem, V.elem);

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

void mat3_svd(const mat3_t M, mat3_t* U, mat3_t* S, mat3_t* V) {
    mat3_t Mt = mat3_transpose(M);
    // the external svd library uses row major matrix convention...
    svd(Mt.elem, U->elem, S->elem, V->elem);
    *U = mat3_transpose(*U);
    *S = mat3_transpose(*S);
    *V = mat3_transpose(*V);
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
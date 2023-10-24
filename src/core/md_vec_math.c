#include <core/md_vec_math.h>

#include <svd3.h>
#include <math.h>

#define SWAP_INT(x, y) {int t = x; x = y; y = t;}

mat3_svd_t mat3_svd(const mat3_t M) {
    mat3_t Mt = mat3_transpose(M);
    // the external svd library uses row major matrix convention...
    mat3_t U,S,V;
    svd((const float(*)[3])Mt.elem, U.elem, S.elem, V.elem);
    
    mat3_svd_t res = {
        .U = mat3_transpose(U),
        .V = mat3_transpose(V),
        .s = {S.elem[0][0], S.elem[1][1], S.elem[2][2]},
    };
    
    return res;
}

mat3_eigen_t mat3_eigen(mat3_t M) {
    mat3_svd_t svd = mat3_svd(M);

    const float max_val = MAX(svd.s.elem[0], MAX(svd.s.elem[1], svd.s.elem[2]));
    const float  e_val[] = {svd.s.elem[0] / max_val, svd.s.elem[1] / max_val, svd.s.elem[2] / max_val};
    const vec3_t e_vec[] = {svd.U.col[0], svd.U.col[1], svd.U.col[2]};

    int l[3] = {0, 1, 2};

    if (e_val[l[0]] < e_val[l[1]]) SWAP_INT(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) SWAP_INT(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) SWAP_INT(l[0], l[1]);

    mat3_eigen_t res;
    res.values = (vec3_t){e_val[l[0]], e_val[l[1]], e_val[l[2]]},
    res.vectors.col[0] = e_vec[l[0]];
    res.vectors.col[1] = e_vec[l[1]];
    res.vectors.col[2] = e_vec[l[2]];
    
    return res;
}

mat3_t mat3_covariance_matrix( const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* indices, vec3_t com, int64_t count) {
    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const int64_t idx = indices ? indices[i] : i;
        const float x = in_x[idx] - com.x;
        const float y = in_y[idx] - com.y;
        const float z = in_z[idx] - com.z;
        const float w = in_w ? in_w[idx] : 1.0f;

        A.elem[0][0] += w * x * x;
        A.elem[0][1] += w * x * y;
        A.elem[0][2] += w * x * z;
        A.elem[1][0] += w * y * x;
        A.elem[1][1] += w * y * y;
        A.elem[1][2] += w * y * z;
        A.elem[2][0] += w * z * x;
        A.elem[2][1] += w * z * y;
        A.elem[2][2] += w * z * z;
    }

    return A;
}

mat3_t mat3_covariance_matrix_vec3(const vec3_t* in_xyz, const float* in_w, const int32_t* indices, vec3_t com, int64_t count) {
    mat3_t A = {0};
    for (int64_t i = 0; i < count; i++) {
        const int64_t idx = indices ? indices[i] : i;
        const float x = in_xyz[idx].x - com.x;
        const float y = in_xyz[idx].y - com.y;
        const float z = in_xyz[idx].z - com.z;
        const float w = in_w ? in_w[idx] : 1.0f;

        A.elem[0][0] += w * x * x;
        A.elem[0][1] += w * x * y;
        A.elem[0][2] += w * x * z;
        A.elem[1][0] += w * y * x;
        A.elem[1][1] += w * y * y;
        A.elem[1][2] += w * y * z;
        A.elem[2][0] += w * z * x;
        A.elem[2][1] += w * z * y;
        A.elem[2][2] += w * z * z;
    }

    return A;
}

mat3_t mat3_cross_covariance_matrix(
    const float* in_x0, const float* in_y0, const float* in_z0,
    const float* in_x1, const float* in_y1, const float* in_z1,
    const float* in_w,
    vec3_t com0,
    vec3_t com1,
    int64_t count) {

    mat3_t A = {0};

    if (in_w) {
        for (int64_t i = 0; i < count; i++) {
            const float px = in_x0[i] - com0.x;
            const float py = in_y0[i] - com0.y;
            const float pz = in_z0[i] - com0.z;

            const float qx = in_x1[i] - com1.x;
            const float qy = in_y1[i] - com1.y;
            const float qz = in_z1[i] - com1.z;

            const float w = in_w[i];

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
    } else {
        for (int64_t i = 0; i < count; i++) {
            const float px = in_x0[i] - com0.x;
            const float py = in_y0[i] - com0.y;
            const float pz = in_z0[i] - com0.z;

            const float qx = in_x1[i] - com1.x;
            const float qy = in_y1[i] - com1.y;
            const float qz = in_z1[i] - com1.z;

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
    }

    return A;
}

mat3_t mat3_extract_rotation(mat3_t M) {
    mat3_svd_t svd = mat3_svd(M);

    mat3_t Ut = mat3_transpose(svd.U);
    float  d = mat3_determinant(mat3_mul(svd.V, Ut));
    mat3_t D = {1, 0, 0, 0, 1, 0, 0, 0, d};
    mat3_t R = mat3_mul(mat3_mul(svd.V, D), Ut);
    return R;
}

mat3_t mat3_optimal_rotation(
    const float* x0, const float* y0, const float* z0,
    const float* x1, const float* y1, const float* z1,
    const float* weight,
    vec3_t com0,
    vec3_t com1,
    int64_t count) {
    
    if (count < 1) {
        return mat3_ident();
    }

    const mat3_t cov_mat = mat3_cross_covariance_matrix(x0, y0, z0, x1, y1, z1, weight, com0, com1, count);
    return mat3_extract_rotation(cov_mat);
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


void vec3_batch_translate_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, int64_t count, vec3_t translation) {
    int64_t i = 0;

    const int64_t simd_count = ROUND_DOWN(count, 8);
    if (simd_count > 0) {
        __m256 t_x = md_mm256_set1_ps(translation.x);
        __m256 t_y = md_mm256_set1_ps(translation.y);
        __m256 t_z = md_mm256_set1_ps(translation.z);

        for (; i < simd_count; i += 8) {
            __m256 x = md_mm256_loadu_ps(in_out_x + i);
            __m256 y = md_mm256_loadu_ps(in_out_y + i);
            __m256 z = md_mm256_loadu_ps(in_out_z + i);

            x = md_mm256_add_ps(x, t_x);
            y = md_mm256_add_ps(y, t_y);
            z = md_mm256_add_ps(z, t_z);

            md_mm256_storeu_ps(in_out_x + i, x);
            md_mm256_storeu_ps(in_out_y + i, y);
            md_mm256_storeu_ps(in_out_z + i, z);
        }
    }

    for (; i < count; i++) {
        in_out_x[i] += translation.x;
        in_out_y[i] += translation.y;
        in_out_z[i] += translation.z;
    }
}

void vec3_batch_translate(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, int64_t count, vec3_t translation) {
    int64_t i = 0;

    const int64_t simd_count = ROUND_DOWN(count, 8);
    if (simd_count > 0) {
        __m256 t_x = md_mm256_set1_ps(translation.x);
        __m256 t_y = md_mm256_set1_ps(translation.y);
        __m256 t_z = md_mm256_set1_ps(translation.z);

        for (; i < simd_count; i += 8) {
            __m256 p_x = md_mm256_loadu_ps(in_x + i);
            __m256 p_y = md_mm256_loadu_ps(in_y + i);
            __m256 p_z = md_mm256_loadu_ps(in_z + i);

            p_x = md_mm256_add_ps(p_x, t_x);
            p_y = md_mm256_add_ps(p_y, t_y);
            p_z = md_mm256_add_ps(p_z, t_z);

            md_mm256_storeu_ps(out_x + i, p_x);
            md_mm256_storeu_ps(out_y + i, p_y);
            md_mm256_storeu_ps(out_z + i, p_z);
        }
    }

    for (; i < count; i++) {
        out_x[i] += translation.x;
        out_y[i] += translation.y;
        out_z[i] += translation.z;
    }
}

void mat4_batch_transform_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, float w_comp, int64_t count, mat4_t M) {
    const __m256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const __m256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const __m256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const __m256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const __m256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const __m256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const __m256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const __m256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const __m256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    const __m256 m41 = md_mm256_set1_ps(M.elem[3][0]);
    const __m256 m42 = md_mm256_set1_ps(M.elem[3][1]);
    const __m256 m43 = md_mm256_set1_ps(M.elem[3][2]);

    const __m256 w = md_mm256_set1_ps(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        const __m256 x = md_mm256_loadu_ps(in_out_x + i);
        const __m256 y = md_mm256_loadu_ps(in_out_y + i);
        const __m256 z = md_mm256_loadu_ps(in_out_z + i);

        const __m256 m11x = md_mm256_mul_ps(m11, x);
        const __m256 m21y = md_mm256_mul_ps(m21, y);
        const __m256 m31z = md_mm256_mul_ps(m31, z);
        const __m256 m41w = md_mm256_mul_ps(m41, w);

        const __m256 m12x = md_mm256_mul_ps(m12, x);
        const __m256 m22y = md_mm256_mul_ps(m22, y);
        const __m256 m32z = md_mm256_mul_ps(m32, z);
        const __m256 m42w = md_mm256_mul_ps(m42, w);

        const __m256 m13x = md_mm256_mul_ps(m13, x);
        const __m256 m23y = md_mm256_mul_ps(m23, y);
        const __m256 m33z = md_mm256_mul_ps(m33, z);
        const __m256 m43w = md_mm256_mul_ps(m43, w);

        const __m256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), md_mm256_add_ps(m31z, m41w));
        const __m256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), md_mm256_add_ps(m32z, m42w));
        const __m256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), md_mm256_add_ps(m33z, m43w));

        md_mm256_storeu_ps(in_out_x + i, res_x);
        md_mm256_storeu_ps(in_out_y + i, res_y);
        md_mm256_storeu_ps(in_out_z + i, res_z);
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

void mat4_batch_transform(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, float w_comp, int64_t count, mat4_t M) {
    const __m256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const __m256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const __m256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const __m256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const __m256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const __m256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const __m256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const __m256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const __m256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    const __m256 m41 = md_mm256_set1_ps(M.elem[3][0]);
    const __m256 m42 = md_mm256_set1_ps(M.elem[3][1]);
    const __m256 m43 = md_mm256_set1_ps(M.elem[3][2]);

    const __m256 w = md_mm256_set1_ps(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        __m256 x = md_mm256_loadu_ps(in_x + i);
        __m256 y = md_mm256_loadu_ps(in_y + i);
        __m256 z = md_mm256_loadu_ps(in_z + i);

        __m256 m11x = md_mm256_mul_ps(m11, x);
        __m256 m21y = md_mm256_mul_ps(m21, y);
        __m256 m31z = md_mm256_mul_ps(m31, z);
        __m256 m41w = md_mm256_mul_ps(m41, w);

        __m256 m12x = md_mm256_mul_ps(m12, x);
        __m256 m22y = md_mm256_mul_ps(m22, y);
        __m256 m32z = md_mm256_mul_ps(m32, z);
        __m256 m42w = md_mm256_mul_ps(m42, w);

        __m256 m13x = md_mm256_mul_ps(m13, x);
        __m256 m23y = md_mm256_mul_ps(m23, y);
        __m256 m33z = md_mm256_mul_ps(m33, z);
        __m256 m43w = md_mm256_mul_ps(m43, w);

        __m256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), md_mm256_add_ps(m31z, m41w));
        __m256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), md_mm256_add_ps(m32z, m42w));
        __m256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), md_mm256_add_ps(m33z, m43w));

        md_mm256_storeu_ps(out_x + i, res_x);
        md_mm256_storeu_ps(out_y + i, res_y);
        md_mm256_storeu_ps(out_z + i, res_z);
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

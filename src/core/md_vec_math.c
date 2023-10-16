#include <core/md_vec_math.h>

#include <svd3.h>
#include <math.h>

#define SWAP(x, y) {int t = x; x = y; y = t;}

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

    if (e_val[l[0]] < e_val[l[1]]) SWAP(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) SWAP(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) SWAP(l[0], l[1]);

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

    const int64_t simd_count = ROUND_DOWN(count, md_simd_width_f32);
    if (simd_count > 0) {
        md_simd_f32_t t_x = md_simd_set1_f32(translation.x);
        md_simd_f32_t t_y = md_simd_set1_f32(translation.y);
        md_simd_f32_t t_z = md_simd_set1_f32(translation.z);

        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t x = md_simd_load_f32(in_out_x + i);
            md_simd_f32_t y = md_simd_load_f32(in_out_y + i);
            md_simd_f32_t z = md_simd_load_f32(in_out_z + i);

            x = md_simd_add(x, t_x);
            y = md_simd_add(y, t_y);
            z = md_simd_add(z, t_z);

            md_simd_store_f32(in_out_x + i, x);
            md_simd_store_f32(in_out_y + i, y);
            md_simd_store_f32(in_out_z + i, z);
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

    const int64_t simd_count = ROUND_DOWN(count, md_simd_width_f32);
    if (simd_count > 0) {
        md_simd_f32_t t_x = md_simd_set1_f32(translation.x);
        md_simd_f32_t t_y = md_simd_set1_f32(translation.y);
        md_simd_f32_t t_z = md_simd_set1_f32(translation.z);

        for (; i < simd_count; i += md_simd_width_f32) {
            md_simd_f32_t p_x = md_simd_load_f32(in_x + i);
            md_simd_f32_t p_y = md_simd_load_f32(in_y + i);
            md_simd_f32_t p_z = md_simd_load_f32(in_z + i);

            p_x = md_simd_add(p_x, t_x);
            p_y = md_simd_add(p_y, t_y);
            p_z = md_simd_add(p_z, t_z);

            md_simd_store(out_x + i, p_x);
            md_simd_store(out_y + i, p_y);
            md_simd_store(out_z + i, p_z);
        }
    }

    for (; i < count; i++) {
        out_x[i] += translation.x;
        out_y[i] += translation.y;
        out_z[i] += translation.z;
    }
}

void mat4_batch_transform_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, float w_comp, int64_t count, mat4_t M) {
    const md_simd_f32_t m11 = md_simd_set1_f32(M.elem[0][0]);
    const md_simd_f32_t m12 = md_simd_set1_f32(M.elem[0][1]);
    const md_simd_f32_t m13 = md_simd_set1_f32(M.elem[0][2]);

    const md_simd_f32_t m21 = md_simd_set1_f32(M.elem[1][0]);
    const md_simd_f32_t m22 = md_simd_set1_f32(M.elem[1][1]);
    const md_simd_f32_t m23 = md_simd_set1_f32(M.elem[1][2]);

    const md_simd_f32_t m31 = md_simd_set1_f32(M.elem[2][0]);
    const md_simd_f32_t m32 = md_simd_set1_f32(M.elem[2][1]);
    const md_simd_f32_t m33 = md_simd_set1_f32(M.elem[2][2]);

    const md_simd_f32_t m41 = md_simd_set1_f32(M.elem[3][0]);
    const md_simd_f32_t m42 = md_simd_set1_f32(M.elem[3][1]);
    const md_simd_f32_t m43 = md_simd_set1_f32(M.elem[3][2]);

    const md_simd_f32_t w = md_simd_set1_f32(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_width_f32);
    for (; i < simd_count; i += md_simd_width_f32) {
        const md_simd_f32_t x = md_simd_load_f32(in_out_x + i);
        const md_simd_f32_t y = md_simd_load_f32(in_out_y + i);
        const md_simd_f32_t z = md_simd_load_f32(in_out_z + i);

        const md_simd_f32_t m11x = md_simd_mul(m11, x);
        const md_simd_f32_t m21y = md_simd_mul(m21, y);
        const md_simd_f32_t m31z = md_simd_mul(m31, z);
        const md_simd_f32_t m41w = md_simd_mul(m41, w);

        const md_simd_f32_t m12x = md_simd_mul(m12, x);
        const md_simd_f32_t m22y = md_simd_mul(m22, y);
        const md_simd_f32_t m32z = md_simd_mul(m32, z);
        const md_simd_f32_t m42w = md_simd_mul(m42, w);

        const md_simd_f32_t m13x = md_simd_mul(m13, x);
        const md_simd_f32_t m23y = md_simd_mul(m23, y);
        const md_simd_f32_t m33z = md_simd_mul(m33, z);
        const md_simd_f32_t m43w = md_simd_mul(m43, w);

        const md_simd_f32_t res_x = md_simd_add(md_simd_add(m11x, m21y), md_simd_add(m31z, m41w));
        const md_simd_f32_t res_y = md_simd_add(md_simd_add(m12x, m22y), md_simd_add(m32z, m42w));
        const md_simd_f32_t res_z = md_simd_add(md_simd_add(m13x, m23y), md_simd_add(m33z, m43w));

        md_simd_store(in_out_x + i, res_x);
        md_simd_store(in_out_y + i, res_y);
        md_simd_store(in_out_z + i, res_z);
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
    const md_simd_f32_t m11 = md_simd_set1_f32(M.elem[0][0]);
    const md_simd_f32_t m12 = md_simd_set1_f32(M.elem[0][1]);
    const md_simd_f32_t m13 = md_simd_set1_f32(M.elem[0][2]);

    const md_simd_f32_t m21 = md_simd_set1_f32(M.elem[1][0]);
    const md_simd_f32_t m22 = md_simd_set1_f32(M.elem[1][1]);
    const md_simd_f32_t m23 = md_simd_set1_f32(M.elem[1][2]);

    const md_simd_f32_t m31 = md_simd_set1_f32(M.elem[2][0]);
    const md_simd_f32_t m32 = md_simd_set1_f32(M.elem[2][1]);
    const md_simd_f32_t m33 = md_simd_set1_f32(M.elem[2][2]);

    const md_simd_f32_t m41 = md_simd_set1_f32(M.elem[3][0]);
    const md_simd_f32_t m42 = md_simd_set1_f32(M.elem[3][1]);
    const md_simd_f32_t m43 = md_simd_set1_f32(M.elem[3][2]);

    const md_simd_f32_t w = md_simd_set1_f32(w_comp);

    int64_t i = 0;
    const int64_t simd_count = ROUND_DOWN(count, md_simd_width_f32);
    for (; i < simd_count; i += md_simd_width_f32) {
        md_simd_f32_t x = md_simd_load_f32(in_x + i);
        md_simd_f32_t y = md_simd_load_f32(in_y + i);
        md_simd_f32_t z = md_simd_load_f32(in_z + i);

        md_simd_f32_t m11x = md_simd_mul(m11, x);
        md_simd_f32_t m21y = md_simd_mul(m21, y);
        md_simd_f32_t m31z = md_simd_mul(m31, z);
        md_simd_f32_t m41w = md_simd_mul(m41, w);

        md_simd_f32_t m12x = md_simd_mul(m12, x);
        md_simd_f32_t m22y = md_simd_mul(m22, y);
        md_simd_f32_t m32z = md_simd_mul(m32, z);
        md_simd_f32_t m42w = md_simd_mul(m42, w);

        md_simd_f32_t m13x = md_simd_mul(m13, x);
        md_simd_f32_t m23y = md_simd_mul(m23, y);
        md_simd_f32_t m33z = md_simd_mul(m33, z);
        md_simd_f32_t m43w = md_simd_mul(m43, w);

        md_simd_f32_t res_x = md_simd_add(md_simd_add(m11x, m21y), md_simd_add(m31z, m41w));
        md_simd_f32_t res_y = md_simd_add(md_simd_add(m12x, m22y), md_simd_add(m32z, m42w));
        md_simd_f32_t res_z = md_simd_add(md_simd_add(m13x, m23y), md_simd_add(m33z, m43w));

        md_simd_store_f32(out_x + i, res_x);
        md_simd_store_f32(out_y + i, res_y);
        md_simd_store_f32(out_z + i, res_z);
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

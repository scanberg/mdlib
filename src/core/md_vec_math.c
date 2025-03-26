#include <core/md_vec_math.h>

#include <svd3.h>

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

mat3_t mat3_covariance_matrix(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, vec3_t in_mean) {   
    // The covariance matrix is symmetric, so we only need to compute the upper triangular part
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        for (size_t i = 0; i < count; i++) {
            const int32_t idx = in_idx[i];
            const float x = in_x[idx] - in_mean.x;
            const float y = in_y[idx] - in_mean.y;
            const float z = in_z[idx] - in_mean.z;
            const float w = in_w ? in_w[idx] : 1.0f;

            A[0][0] += w * x * x;
            A[0][1] += w * x * y;
            A[0][2] += w * x * z;
            A[1][0] += w * y * x;
            A[1][1] += w * y * y;
            A[1][2] += w * y * z;
            A[2][0] += w * z * x;
            A[2][1] += w * z * y;
            A[2][2] += w * z * z;
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; i++) {
            const float x = in_x[i] - in_mean.x;
            const float y = in_y[i] - in_mean.y;
            const float z = in_z[i] - in_mean.z;
            const float w = in_w ? in_w[i] : 1.0f;

            A[0][0] += w * x * x;
            A[0][1] += w * x * y;
            A[0][2] += w * x * z;
            A[1][0] += w * y * x;
            A[1][1] += w * y * y;
            A[1][2] += w * y * z;
            A[2][0] += w * z * x;
            A[2][1] += w * z * y;
            A[2][2] += w * z * z;
            w_sum += w;
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            A[i][j] /= w_sum;
        }
    }

    return (mat3_t) {
        (float)A[0][0], (float)A[0][1], (float)A[0][2],
        (float)A[1][0], (float)A[1][1], (float)A[1][2],
        (float)A[2][0], (float)A[2][1], (float)A[2][2],
    };
}

mat3_t mat3_covariance_matrix_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, vec3_t com) {
    // The covariance matrix is symmetric, so we only need to compute the upper triangular part
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        for (size_t i = 0; i < count; i++) {
            const int32_t idx = in_idx[i];
            const float x = in_xyzw[idx].x - com.x;
            const float y = in_xyzw[idx].y - com.y;
            const float z = in_xyzw[idx].z - com.z;
            const float w = in_xyzw[idx].w;

            A[0][0] += w * x * x;
            A[0][1] += w * x * y;
            A[0][2] += w * x * z;
            A[1][0] += w * y * x;
            A[1][1] += w * y * y;
            A[1][2] += w * y * z;
            A[2][0] += w * z * x;
            A[2][1] += w * z * y;
            A[2][2] += w * z * z;
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; i++) {
            const float x = in_xyzw[i].x - com.x;
            const float y = in_xyzw[i].y - com.y;
            const float z = in_xyzw[i].z - com.z;
            const float w = in_xyzw[i].w;

            A[0][0] += w * x * x;
            A[0][1] += w * x * y;
            A[0][2] += w * x * z;
            A[1][0] += w * y * x;
            A[1][1] += w * y * y;
            A[1][2] += w * y * z;
            A[2][0] += w * z * x;
            A[2][1] += w * z * y;
            A[2][2] += w * z * z;
            w_sum += w;
        }
    }

    for (size_t i = 0; i < 3; i++) {
    	for (size_t j = 0; j < 3; j++) {
            A[i][j] /= w_sum;
        }
    }

    return (mat3_t) {
        (float)A[0][0], (float)A[0][1], (float)A[0][2],
        (float)A[1][0], (float)A[1][1], (float)A[1][2],
        (float)A[2][0], (float)A[2][1], (float)A[2][2],
    };
}

mat3_t mat3_cross_covariance_matrix(const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        ASSERT(in_idx[0]);
        ASSERT(in_idx[1]);
        for (size_t i = 0; i < count; i++) {
            const int32_t i0 = in_idx[0][i];
            const int32_t i1 = in_idx[1][i];
            const float px = in_x[0][i0] - com[0].x;
            const float py = in_y[0][i0] - com[0].y;
            const float pz = in_z[0][i0] - com[0].z;

            const float qx = in_x[1][i1] - com[1].x;
            const float qy = in_y[1][i1] - com[1].y;
            const float qz = in_z[1][i1] - com[1].z;

            const float w = in_w ? (in_w[0][i0] + in_w[1][i1]) * 0.5f : 1.0f;

            A[0][0] += w * px * qx;
            A[0][1] += w * px * qy;
            A[0][2] += w * px * qz;
            A[1][0] += w * py * qx;
            A[1][1] += w * py * qy;
            A[1][2] += w * py * qz;
            A[2][0] += w * pz * qx;
            A[2][1] += w * pz * qy;
            A[2][2] += w * pz * qz;
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; i++) {
            const float px = in_x[0][i] - com[0].x;
            const float py = in_y[0][i] - com[0].y;
            const float pz = in_z[0][i] - com[0].z;

            const float qx = in_x[1][i] - com[1].x;
            const float qy = in_y[1][i] - com[1].y;
            const float qz = in_z[1][i] - com[1].z;

            const float w = in_w ? (in_w[0][i] + in_w[1][i]) * 0.5f : 1.0f;

            A[0][0] += w * px * qx;
            A[0][1] += w * px * qy;
            A[0][2] += w * px * qz;
            A[1][0] += w * py * qx;
            A[1][1] += w * py * qy;
            A[1][2] += w * py * qz;
            A[2][0] += w * pz * qx;
            A[2][1] += w * pz * qy;
            A[2][2] += w * pz * qz;
            w_sum += w;
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            A[i][j] /= w_sum;
        }
    }

    return (mat3_t) {
        (float)A[0][0], (float)A[0][1], (float)A[0][2],
        (float)A[1][0], (float)A[1][1], (float)A[1][2],
        (float)A[2][0], (float)A[2][1], (float)A[2][2],
    };
}

mat3_t mat3_cross_covariance_matrix_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        ASSERT(in_idx[0]);
        ASSERT(in_idx[1]);
        for (size_t i = 0; i < count; ++i) {
            const int32_t i0 = in_idx[0][i];
            const int32_t i1 = in_idx[1][i];
            vec4_t p = vec4_sub(in_xyzw[0][i0], vec4_from_vec3(com[0], 0));
            vec4_t q = vec4_sub(in_xyzw[1][i1], vec4_from_vec3(com[1], 0));

            // The question here is how to combine the weights.
            // For now we just take the average.
            // This should be equivalent to the other case where the weights for both sets are equal.
            const float w = (p.w + q.w) * 0.5f;

            A[0][0] += w * p.x * q.x;
            A[0][1] += w * p.x * q.y;
            A[0][2] += w * p.x * q.z;
            A[1][0] += w * p.y * q.x;
            A[1][1] += w * p.y * q.y;
            A[1][2] += w * p.y * q.z;
            A[2][0] += w * p.z * q.x;
            A[2][1] += w * p.z * q.y;
            A[2][2] += w * p.z * q.z;
            w_sum += w;
        }
    } else {
        for (size_t i = 0; i < count; ++i) {
            vec4_t p = vec4_sub(in_xyzw[0][i], vec4_from_vec3(com[0], 0));
            vec4_t q = vec4_sub(in_xyzw[1][i], vec4_from_vec3(com[1], 0));

            // The question here is how to combine the weights.
            // For now we just take the average.
            // This should be equivalent to the other case where the weights for both sets are equal.
            const float w = (p.w + q.w) * 0.5f;

            A[0][0] += w * p.x * q.x;
            A[0][1] += w * p.x * q.y;
            A[0][2] += w * p.x * q.z;
            A[1][0] += w * p.y * q.x;
            A[1][1] += w * p.y * q.y;
            A[1][2] += w * p.y * q.z;
            A[2][0] += w * p.z * q.x;
            A[2][1] += w * p.z * q.y;
            A[2][2] += w * p.z * q.z;
            w_sum += w;
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            A[i][j] /= w_sum;
        }
    }

    return (mat3_t) {
        (float)A[0][0], (float)A[0][1], (float)A[0][2],
        (float)A[1][0], (float)A[1][1], (float)A[1][2],
        (float)A[2][0], (float)A[2][1], (float)A[2][2],
    };
}

mat3_t mat3_extract_rotation(mat3_t M) {
    mat3_svd_t svd = mat3_svd(M);

    mat3_t Ut = mat3_transpose(svd.U);
    float  d = mat3_determinant(mat3_mul(svd.V, Ut));
    mat3_t D = {1, 0, 0, 0, 1, 0, 0, 0, d};
    mat3_t R = mat3_mul(mat3_mul(svd.V, D), Ut);
    return R;
}

static double highp_dot(vec3_t a, vec3_t b) {
    return (double)a.x * (double)b.x + (double)a.y * (double)b.y + (double)a.z * (double)b.z;
}

static vec3_t highp_normalize(vec3_t v) {
    double len = sqrt(highp_dot(v, v));
    vec3_t result = {
        (float)(v.x / len),
        (float)(v.y / len),
        (float)(v.z / len),
    };
    return result;
}

mat3_t mat3_orthonormalize(mat3_t M) {
    M.col[0] = vec3_normalize(M.col[0]);

    M.col[1] = vec3_sub(M.col[1], vec3_mul_f(M.col[0], vec3_dot(M.col[0], M.col[1])));
    M.col[1] = vec3_normalize(M.col[1]);

    M.col[2] = vec3_sub(M.col[2], vec3_add(vec3_mul_f(M.col[0], vec3_dot(M.col[0], M.col[2])), vec3_mul_f(M.col[1], vec3_dot(M.col[1], M.col[2]))));
    M.col[2] = vec3_normalize(M.col[2]);

    return M;
}

mat3_t mat3_optimal_rotation(const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    if (count < 1) {
        return mat3_ident();
    }

    const mat3_t cov_mat = mat3_cross_covariance_matrix(in_x, in_y, in_z, in_w, in_idx, count, com);
    return mat3_extract_rotation(cov_mat);
}

mat3_t mat3_optimal_rotation_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    if (count < 1) {
		return mat3_ident();
	}

	const mat3_t cov_mat = mat3_cross_covariance_matrix_vec4(in_xyzw, in_idx, count, com);
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

vec3_t mat4_unproject(vec3_t window_coords, mat4_t inv_view_proj_mat, vec4_t viewport) {
    vec4_t tmp = vec4_from_vec3(window_coords, 1.f);
    tmp.x = (tmp.x - viewport.elem[0]) / viewport.elem[2];
    tmp.y = (tmp.y - viewport.elem[1]) / viewport.elem[3];
    tmp = vec4_sub_f(vec4_mul_f(tmp, 2.f), 1.f);

    vec4_t obj = mat4_mul_vec4(inv_view_proj_mat, tmp);
    obj = vec4_div_f(obj, obj.w);

    return vec3_from_vec4(obj);
}

mat4_t mat4_look_at(vec3_t look_from, vec3_t look_at, vec3_t look_up) {
    const vec3_t f = vec3_normalize(vec3_sub(look_at, look_from));
    const vec3_t s = vec3_normalize(vec3_cross(f, look_up));
    const vec3_t u = vec3_cross(s, f);
    const mat4_t M = {
        s.x, u.x, -f.x, 0.0f,
        s.y, u.y, -f.y, 0.0f,
        s.z, u.z, -f.z, 0.0f,
        -vec3_dot(s, look_from), -vec3_dot(u, look_from), vec3_dot(f, look_from), 1.0f,
    };
    return M;
}

mat4_t mat4_ortho(float l, float r, float b, float t, float n, float f) {
    mat4_t M = {0};
    M.elem[0][0] = 2 / (r-l);
    M.elem[1][1] = 2 / (t-b);
    M.elem[2][2] = -2 / (f-n);
    M.elem[3][0] = -(r+l) / (r-l);
    M.elem[3][1] = -(t+b) / (t-b);
    M.elem[3][2] = -(f+n) / (f-n);
    M.elem[3][3] = 1;
    return M;
}

mat4_t mat4_ortho_inv(float l, float r, float b, float t, float n, float f) {
    mat4_t M = {0};
    M.elem[0][0] = (r-l) / 2;
    M.elem[1][1] = (t-b) / 2;
    M.elem[2][2] = (n-f) / 2;
    M.elem[3][0] = (l+r) / 2;
    M.elem[3][1] = (b+t) / 2;
    M.elem[3][2] = -(n+f) / 2;
    M.elem[3][3] = 1;
    return M;
}

mat4_t mat4_ortho_2d(float l, float r, float b, float t) {
    mat4_t M = {0};
    M.elem[0][0] = 2 / (r-l);
    M.elem[1][1] = 2 / (t-b);
    M.elem[2][2] = -1;
    M.elem[3][0] = -(r+l) / (r-l);
    M.elem[3][1] = -(t+b) / (t-b);
    M.elem[3][3] = 1;
    return M;
}

mat4_t mat4_ortho_2d_inv(float l, float r, float b, float t) {
    mat4_t M = {0};
    M.elem[0][0] = (r-l) / 2;
    M.elem[1][1] = (t-b) / 2;
    M.elem[2][2] = -1;
    M.elem[3][0] = (l+r) / 2;
    M.elem[3][1] = (b+t) / 2;
    M.elem[3][3] = 1;
    return M;
}

mat4_t mat4_persp(float fovy, float aspect, float near, float far) {
    const float tan_half_fovy = tanf(fovy * 0.5f);
    mat4_t M = {0};
    M.elem[0][0] = 1.0f / (aspect * tan_half_fovy);
    M.elem[1][1] = 1.0f / (tan_half_fovy);
    M.elem[2][2] = -(far + near) / (far - near);
    M.elem[2][3] = -1;
    M.elem[3][2] = -(2 * far * near) / (far - near);
    return M;
}

mat4_t mat4_persp_inv(float fovy, float aspect, float near, float far) {
    const float tan_half_fovy = tanf(fovy * 0.5f);
    mat4_t M = {0};
    M.elem[0][0] = aspect * tan_half_fovy;
    M.elem[1][1] = tan_half_fovy;
    M.elem[2][3] = (near - far) / (2 * far * near);
    M.elem[3][2] = -1;
    M.elem[3][3] = (near + far) / (2 * far * near);
    return M;
}

mat4_t mat4_frustum(float l, float r, float b, float t, float n, float f) {
    mat4_t M = {0};
    M.elem[0][0] = (2*n) / (r-l);
    M.elem[1][1] = (2*n) / (t-b);
    M.elem[2][0] = (r+l) / (r-l);
    M.elem[2][1] = (t+b) / (t-b);
    M.elem[2][2] = -(f+n) / (f-n);
    M.elem[2][3] = -1;
    M.elem[3][2] = -(2*n*f) / (f-n);
    return M;
}

mat4_t mat4_frustum_inv(float l, float r, float b, float t, float n, float f) {
    mat4_t M = {0};
    M.elem[0][0] = (r-l) / (2*n);
    M.elem[1][1] = (t-b) / (2*n);
    M.elem[2][3] = (n-f) / (2*n*f);
    M.elem[3][0] = (l+r) / (2*n);
    M.elem[3][1] = (b+t) / (2*n);
    M.elem[3][2] = -1;
    M.elem[3][3] = (n+f) / (2*n*f);
    return M;
}

void vec3_batch_translate_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, size_t count, vec3_t translation) {
    size_t i = 0;

    const size_t simd_count = ROUND_DOWN(count, 8);
    if (simd_count > 0) {
        md_256 t_x = md_mm256_set1_ps(translation.x);
        md_256 t_y = md_mm256_set1_ps(translation.y);
        md_256 t_z = md_mm256_set1_ps(translation.z);

        for (; i < simd_count; i += 8) {
            md_256 x = md_mm256_loadu_ps(in_out_x + i);
            md_256 y = md_mm256_loadu_ps(in_out_y + i);
            md_256 z = md_mm256_loadu_ps(in_out_z + i);

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

void vec3_batch_translate(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, size_t count, vec3_t translation) {
    size_t i = 0;

    const size_t simd_count = ROUND_DOWN(count, 8);
    if (simd_count > 0) {
        md_256 t_x = md_mm256_set1_ps(translation.x);
        md_256 t_y = md_mm256_set1_ps(translation.y);
        md_256 t_z = md_mm256_set1_ps(translation.z);

        for (; i < simd_count; i += 8) {
            md_256 p_x = md_mm256_loadu_ps(in_x + i);
            md_256 p_y = md_mm256_loadu_ps(in_y + i);
            md_256 p_z = md_mm256_loadu_ps(in_z + i);

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

void mat3_batch_transform_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, size_t count, mat3_t M) {
    const md_256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const md_256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const md_256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const md_256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const md_256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const md_256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const md_256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const md_256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const md_256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    size_t i = 0;
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        const md_256 x = md_mm256_loadu_ps(in_out_x + i);
        const md_256 y = md_mm256_loadu_ps(in_out_y + i);
        const md_256 z = md_mm256_loadu_ps(in_out_z + i);

        const md_256 m11x = md_mm256_mul_ps(m11, x);
        const md_256 m21y = md_mm256_mul_ps(m21, y);
        const md_256 m31z = md_mm256_mul_ps(m31, z);

        const md_256 m12x = md_mm256_mul_ps(m12, x);
        const md_256 m22y = md_mm256_mul_ps(m22, y);
        const md_256 m32z = md_mm256_mul_ps(m32, z);

        const md_256 m13x = md_mm256_mul_ps(m13, x);
        const md_256 m23y = md_mm256_mul_ps(m23, y);
        const md_256 m33z = md_mm256_mul_ps(m33, z);

        const md_256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), m31z);
        const md_256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), m32z);
        const md_256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), m33z);

        md_mm256_storeu_ps(in_out_x + i, res_x);
        md_mm256_storeu_ps(in_out_y + i, res_y);
        md_mm256_storeu_ps(in_out_z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_out_x[i];
        const float y = in_out_y[i];
        const float z = in_out_z[i];

        in_out_x[i] = x * M.elem[0][0] + y * M.elem[1][0] + z * M.elem[2][0];
        in_out_y[i] = x * M.elem[0][1] + y * M.elem[1][1] + z * M.elem[2][1];
        in_out_z[i] = x * M.elem[0][2] + y * M.elem[1][2] + z * M.elem[2][2];
    }
}

void mat3_batch_transform(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, size_t count, mat3_t M) {
    const md_256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const md_256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const md_256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const md_256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const md_256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const md_256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const md_256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const md_256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const md_256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    size_t i = 0;
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        md_256 x = md_mm256_loadu_ps(in_x + i);
        md_256 y = md_mm256_loadu_ps(in_y + i);
        md_256 z = md_mm256_loadu_ps(in_z + i);

        md_256 m11x = md_mm256_mul_ps(m11, x);
        md_256 m21y = md_mm256_mul_ps(m21, y);
        md_256 m31z = md_mm256_mul_ps(m31, z);

        md_256 m12x = md_mm256_mul_ps(m12, x);
        md_256 m22y = md_mm256_mul_ps(m22, y);
        md_256 m32z = md_mm256_mul_ps(m32, z);

        md_256 m13x = md_mm256_mul_ps(m13, x);
        md_256 m23y = md_mm256_mul_ps(m23, y);
        md_256 m33z = md_mm256_mul_ps(m33, z);

        md_256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), m31z);
        md_256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), m32z);
        md_256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), m33z);

        md_mm256_storeu_ps(out_x + i, res_x);
        md_mm256_storeu_ps(out_y + i, res_y);
        md_mm256_storeu_ps(out_z + i, res_z);
    }

    for (; i < count; i++) {
        const float x = in_x[i];
        const float y = in_y[i];
        const float z = in_z[i];

        out_x[i] = x * M.elem[0][0] + y * M.elem[1][0] + z * M.elem[2][0];
        out_y[i] = x * M.elem[0][1] + y * M.elem[1][1] + z * M.elem[2][1];
        out_z[i] = x * M.elem[0][2] + y * M.elem[1][2] + z * M.elem[2][2];
    }
}

void mat4_batch_transform_inplace(float* RESTRICT in_out_x, float* RESTRICT in_out_y, float* RESTRICT in_out_z, float w_comp, size_t count, mat4_t M) {
    const md_256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const md_256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const md_256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const md_256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const md_256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const md_256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const md_256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const md_256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const md_256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    const md_256 m41 = md_mm256_set1_ps(M.elem[3][0]);
    const md_256 m42 = md_mm256_set1_ps(M.elem[3][1]);
    const md_256 m43 = md_mm256_set1_ps(M.elem[3][2]);

    const md_256 w = md_mm256_set1_ps(w_comp);

    size_t i = 0;
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        const md_256 x = md_mm256_loadu_ps(in_out_x + i);
        const md_256 y = md_mm256_loadu_ps(in_out_y + i);
        const md_256 z = md_mm256_loadu_ps(in_out_z + i);

        const md_256 m11x = md_mm256_mul_ps(m11, x);
        const md_256 m21y = md_mm256_mul_ps(m21, y);
        const md_256 m31z = md_mm256_mul_ps(m31, z);
        const md_256 m41w = md_mm256_mul_ps(m41, w);

        const md_256 m12x = md_mm256_mul_ps(m12, x);
        const md_256 m22y = md_mm256_mul_ps(m22, y);
        const md_256 m32z = md_mm256_mul_ps(m32, z);
        const md_256 m42w = md_mm256_mul_ps(m42, w);

        const md_256 m13x = md_mm256_mul_ps(m13, x);
        const md_256 m23y = md_mm256_mul_ps(m23, y);
        const md_256 m33z = md_mm256_mul_ps(m33, z);
        const md_256 m43w = md_mm256_mul_ps(m43, w);

        const md_256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), md_mm256_add_ps(m31z, m41w));
        const md_256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), md_mm256_add_ps(m32z, m42w));
        const md_256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), md_mm256_add_ps(m33z, m43w));

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

void mat4_batch_transform(float* out_x, float* out_y, float* out_z, const float* in_x, const float* in_y, const float* in_z, float w_comp, size_t count, mat4_t M) {
    const md_256 m11 = md_mm256_set1_ps(M.elem[0][0]);
    const md_256 m12 = md_mm256_set1_ps(M.elem[0][1]);
    const md_256 m13 = md_mm256_set1_ps(M.elem[0][2]);

    const md_256 m21 = md_mm256_set1_ps(M.elem[1][0]);
    const md_256 m22 = md_mm256_set1_ps(M.elem[1][1]);
    const md_256 m23 = md_mm256_set1_ps(M.elem[1][2]);

    const md_256 m31 = md_mm256_set1_ps(M.elem[2][0]);
    const md_256 m32 = md_mm256_set1_ps(M.elem[2][1]);
    const md_256 m33 = md_mm256_set1_ps(M.elem[2][2]);

    const md_256 m41 = md_mm256_set1_ps(M.elem[3][0]);
    const md_256 m42 = md_mm256_set1_ps(M.elem[3][1]);
    const md_256 m43 = md_mm256_set1_ps(M.elem[3][2]);

    const md_256 w = md_mm256_set1_ps(w_comp);

    size_t i = 0;
    const size_t simd_count = ROUND_DOWN(count, 8);
    for (; i < simd_count; i += 8) {
        md_256 x = md_mm256_loadu_ps(in_x + i);
        md_256 y = md_mm256_loadu_ps(in_y + i);
        md_256 z = md_mm256_loadu_ps(in_z + i);

        md_256 m11x = md_mm256_mul_ps(m11, x);
        md_256 m21y = md_mm256_mul_ps(m21, y);
        md_256 m31z = md_mm256_mul_ps(m31, z);
        md_256 m41w = md_mm256_mul_ps(m41, w);

        md_256 m12x = md_mm256_mul_ps(m12, x);
        md_256 m22y = md_mm256_mul_ps(m22, y);
        md_256 m32z = md_mm256_mul_ps(m32, z);
        md_256 m42w = md_mm256_mul_ps(m42, w);

        md_256 m13x = md_mm256_mul_ps(m13, x);
        md_256 m23y = md_mm256_mul_ps(m23, y);
        md_256 m33z = md_mm256_mul_ps(m33, z);
        md_256 m43w = md_mm256_mul_ps(m43, w);

        md_256 res_x = md_mm256_add_ps(md_mm256_add_ps(m11x, m21y), md_mm256_add_ps(m31z, m41w));
        md_256 res_y = md_mm256_add_ps(md_mm256_add_ps(m12x, m22y), md_mm256_add_ps(m32z, m42w));
        md_256 res_z = md_mm256_add_ps(md_mm256_add_ps(m13x, m23y), md_mm256_add_ps(m33z, m43w));

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

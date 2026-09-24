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

mat3_t mat3_covariance_matrix(const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, vec3_t in_mean) {   
    // The covariance matrix is symmetric, so we only need to compute the upper triangular part
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        for (size_t i = 0; i < count; i++) {
            const int32_t idx = in_idx[i];
            const float x = in_xyz[idx].x - in_mean.x;
            const float y = in_xyz[idx].y - in_mean.y;
            const float z = in_xyz[idx].z - in_mean.z;
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
            const float x = in_xyz[i].x - in_mean.x;
            const float y = in_xyz[i].y - in_mean.y;
            const float z = in_xyz[i].z - in_mean.z;
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

mat3_t mat3_cross_covariance_matrix(const vec3_t* const in_xyz[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    double A[3][3] = {0};
    double w_sum = 0.0;

    if (in_idx) {
        ASSERT(in_idx[0]);
        ASSERT(in_idx[1]);
        for (size_t i = 0; i < count; i++) {
            const int32_t i0 = in_idx[0][i];
            const int32_t i1 = in_idx[1][i];
            const float px = in_xyz[0][i0].x - com[0].x;
            const float py = in_xyz[0][i0].y - com[0].y;
            const float pz = in_xyz[0][i0].z - com[0].z;

            const float qx = in_xyz[1][i1].x - com[1].x;
            const float qy = in_xyz[1][i1].y - com[1].y;
            const float qz = in_xyz[1][i1].z - com[1].z;

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
        // A native vector of atoms of each set at a time, split into x, y and z: every entry pairs one axis of
        // the first with one of the second. The weight goes on the first factor.
        size_t i = 0;
        md_xv acc[3][3];
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) acc[r][c] = md_xv_setzero_ps();
        md_xv acc_w = md_xv_setzero_ps();
        const md_xv c0[3] = { md_xv_set1_ps(com[0].x), md_xv_set1_ps(com[0].y), md_xv_set1_ps(com[0].z) };
        const md_xv c1[3] = { md_xv_set1_ps(com[1].x), md_xv_set1_ps(com[1].y), md_xv_set1_ps(com[1].z) };
        for (; i + MD_XV_WIDTH <= count; i += MD_XV_WIDTH) {
            md_xv p[3], q[3];
            md_xv_load_xyz_packed_ps(&p[0], &p[1], &p[2], (const float*)(in_xyz[0] + i));
            md_xv_load_xyz_packed_ps(&q[0], &q[1], &q[2], (const float*)(in_xyz[1] + i));
            md_xv w = md_xv_set1_ps(1.0f);
            if (in_w) {
                w = md_xv_mul_ps(md_xv_add_ps(md_xv_loadu_ps(in_w[0] + i), md_xv_loadu_ps(in_w[1] + i)), md_xv_set1_ps(0.5f));
            }
            acc_w = md_xv_add_ps(acc_w, w);
            for (int k = 0; k < 3; ++k) {
                p[k] = md_xv_mul_ps(md_xv_sub_ps(p[k], c0[k]), w);
                q[k] = md_xv_sub_ps(q[k], c1[k]);
            }
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    acc[r][c] = md_xv_fmadd_ps(p[r], q[c], acc[r][c]);
                }
            }
        }
        for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) A[r][c] += md_xv_reduce_add_ps(acc[r][c]);
        w_sum += md_xv_reduce_add_ps(acc_w);

        for (; i < count; i++) {
            const float px = in_xyz[0][i].x - com[0].x;
            const float py = in_xyz[0][i].y - com[0].y;
            const float pz = in_xyz[0][i].z - com[0].z;

            const float qx = in_xyz[1][i].x - com[1].x;
            const float qy = in_xyz[1][i].y - com[1].y;
            const float qz = in_xyz[1][i].z - com[1].z;

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

mat3_t mat3_cross_covariance_matrix_raw_vec4(const vec4_t* const in_xyzw[2], size_t count) {
    double A[3][3] = {0};
    double w_sum = 0.0;

    for (size_t i = 0; i < count; ++i) {
        vec4_t p = in_xyzw[0][i];
        vec4_t q = in_xyzw[1][i];

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
    const vec4_t com0 = vec4_from_vec3(com[0], 0);
    const vec4_t com1 = vec4_from_vec3(com[1], 0);

    if (in_idx) {
        ASSERT(in_idx[0]);
        ASSERT(in_idx[1]);
        for (size_t i = 0; i < count; ++i) {
            const int32_t i0 = in_idx[0][i];
            const int32_t i1 = in_idx[1][i];
            vec4_t p = vec4_sub(in_xyzw[0][i0], com0);
            vec4_t q = vec4_sub(in_xyzw[1][i1], com1);

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
            vec4_t p = vec4_sub(in_xyzw[0][i], com0);
            vec4_t q = vec4_sub(in_xyzw[1][i], com1);

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
    float  d  = mat3_determinant(mat3_mul(svd.V, Ut));
    float  s  = signf(d);
    mat3_t D  = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, s};
    mat3_t R  = mat3_mul(mat3_mul(svd.V, D), Ut);
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

    M.col[1] = vec3_sub(M.col[1], vec3_mul1(M.col[0], vec3_dot(M.col[0], M.col[1])));
    M.col[1] = vec3_normalize(M.col[1]);

    M.col[2] = vec3_sub(M.col[2], vec3_add(vec3_mul1(M.col[0], vec3_dot(M.col[0], M.col[2])), vec3_mul1(M.col[1], vec3_dot(M.col[1], M.col[2]))));
    M.col[2] = vec3_normalize(M.col[2]);

    return M;
}

mat3_t mat3_optimal_rotation(const vec3_t* const in_xyz[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t com[2]) {
    if (count < 1) {
        return mat3_ident();
    }

    const mat3_t cov_mat = mat3_cross_covariance_matrix(in_xyz, in_w, in_idx, count, com);
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

    return mat4_mul1(I, 1.0f / (dot0.x + dot0.y + dot0.z + dot0.w));
}

vec3_t mat4_unproject(vec3_t window_coords, mat4_t inv_view_proj_mat, vec4_t viewport) {
    vec4_t tmp = vec4_from_vec3(window_coords, 1.f);
    tmp.x = (tmp.x - viewport.elem[0]) / viewport.elem[2];
    tmp.y = (tmp.y - viewport.elem[1]) / viewport.elem[3];
    tmp = vec4_sub1(vec4_mul1(tmp, 2.f), 1.f);

    vec4_t obj = mat4_mul_vec4(inv_view_proj_mat, tmp);
    obj = vec4_div1(obj, obj.w);

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

// NOT 'near' and 'far': windef.h, which windows.h drags in, does '#define near' and
// '#define far' as empty macros for the old segmented-memory keywords, so parameters by
// those names preprocess away and leave a bare 'float,' in the list. It only shows up in
// builds whose include chain reaches that header, so it reads as an architecture problem
// and is not one.
mat4_t mat4_persp(float fovy, float aspect, float z_near, float z_far) {
    const float tan_half_fovy = tanf(fovy * 0.5f);
    mat4_t M = {0};
    M.elem[0][0] = 1.0f / (aspect * tan_half_fovy);
    M.elem[1][1] = 1.0f / (tan_half_fovy);
    M.elem[2][2] = -(z_far + z_near) / (z_far - z_near);
    M.elem[2][3] = -1;
    M.elem[3][2] = -(2 * z_far * z_near) / (z_far - z_near);
    return M;
}

mat4_t mat4_persp_inv(float fovy, float aspect, float z_near, float z_far) {
    const float tan_half_fovy = tanf(fovy * 0.5f);
    mat4_t M = {0};
    M.elem[0][0] = aspect * tan_half_fovy;
    M.elem[1][1] = tan_half_fovy;
    M.elem[2][3] = (z_near - z_far) / (2 * z_far * z_near);
    M.elem[3][2] = -1;
    M.elem[3][3] = (z_near + z_far) / (2 * z_far * z_near);
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

// Translation is the same for every atom, so packed coordinates need no split: 24 floats are eight
// atoms, and the translation repeats with them every three floats.
void vec3_batch_translate(vec3_t* out_xyz, const vec3_t* in_xyz, size_t count, vec3_t t) {
    const float* src = (const float*)in_xyz;
    float* dst = (float*)out_xyz;
    enum { N = 24 / MD_XV_WIDTH };
    float pattern[24];
    for (int k = 0; k < 24; ++k) pattern[k] = t.elem[k % 3];
    md_xv tv[N];
    for (int j = 0; j < N; ++j) tv[j] = md_xv_loadu_ps(pattern + j * MD_XV_WIDTH);

    size_t i = 0;
    for (; i + 8 <= count; i += 8) {
        for (int j = 0; j < N; ++j) {
            const md_xv v = md_xv_loadu_ps(src + i * 3 + j * MD_XV_WIDTH);
            md_xv_storeu_ps(dst + i * 3 + j * MD_XV_WIDTH, md_xv_add_ps(v, tv[j]));
        }
    }
    for (; i < count; i++) {
        out_xyz[i] = vec3_add(in_xyz[i], t);
    }
}

void vec3_batch_translate_inplace(vec3_t* in_out_xyz, size_t count, vec3_t translation) {
    vec3_batch_translate(in_out_xyz, in_out_xyz, count, translation);
}

// A native vector of atoms at a time, split into x, y and z and packed again. In place is allowed:
// each group is read before it is written.
void mat3_batch_transform(vec3_t* out_xyz, const vec3_t* in_xyz, size_t count, mat3_t M) {
    const md_xv m11 = md_xv_set1_ps(M.elem[0][0]);
    const md_xv m12 = md_xv_set1_ps(M.elem[0][1]);
    const md_xv m13 = md_xv_set1_ps(M.elem[0][2]);

    const md_xv m21 = md_xv_set1_ps(M.elem[1][0]);
    const md_xv m22 = md_xv_set1_ps(M.elem[1][1]);
    const md_xv m23 = md_xv_set1_ps(M.elem[1][2]);

    const md_xv m31 = md_xv_set1_ps(M.elem[2][0]);
    const md_xv m32 = md_xv_set1_ps(M.elem[2][1]);
    const md_xv m33 = md_xv_set1_ps(M.elem[2][2]);

    size_t i = 0;
    for (; i + MD_XV_WIDTH <= count; i += MD_XV_WIDTH) {
        md_xv x, y, z;
        md_xv_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

        md_xv rx = md_xv_mul_ps(m11, x);
        md_xv ry = md_xv_mul_ps(m12, x);
        md_xv rz = md_xv_mul_ps(m13, x);

        rx = md_xv_fmadd_ps(m21, y, rx);
        ry = md_xv_fmadd_ps(m22, y, ry);
        rz = md_xv_fmadd_ps(m23, y, rz);

        rx = md_xv_fmadd_ps(m31, z, rx);
        ry = md_xv_fmadd_ps(m32, z, ry);
        rz = md_xv_fmadd_ps(m33, z, rz);

        md_xv_store_xyz_packed_ps((float*)(out_xyz + i), rx, ry, rz);
    }

    for (; i < count; i++) {
        const vec3_t p = in_xyz[i];
        out_xyz[i] = (vec3_t) {
            p.x * M.elem[0][0] + p.y * M.elem[1][0] + p.z * M.elem[2][0],
            p.x * M.elem[0][1] + p.y * M.elem[1][1] + p.z * M.elem[2][1],
            p.x * M.elem[0][2] + p.y * M.elem[1][2] + p.z * M.elem[2][2],
        };
    }
}

void mat3_batch_transform_inplace(vec3_t* in_out_xyz, size_t count, mat3_t M) {
    mat3_batch_transform(in_out_xyz, in_out_xyz, count, M);
}

void mat4_batch_transform(vec3_t* out_xyz, const vec3_t* in_xyz, float w_comp, size_t count, mat4_t M) {
    const md_xv m11 = md_xv_set1_ps(M.elem[0][0]);
    const md_xv m12 = md_xv_set1_ps(M.elem[0][1]);
    const md_xv m13 = md_xv_set1_ps(M.elem[0][2]);

    const md_xv m21 = md_xv_set1_ps(M.elem[1][0]);
    const md_xv m22 = md_xv_set1_ps(M.elem[1][1]);
    const md_xv m23 = md_xv_set1_ps(M.elem[1][2]);

    const md_xv m31 = md_xv_set1_ps(M.elem[2][0]);
    const md_xv m32 = md_xv_set1_ps(M.elem[2][1]);
    const md_xv m33 = md_xv_set1_ps(M.elem[2][2]);

    // The fourth row times w is the same for every atom
    const md_xv w  = md_xv_set1_ps(w_comp);
    const md_xv tx = md_xv_mul_ps(md_xv_set1_ps(M.elem[3][0]), w);
    const md_xv ty = md_xv_mul_ps(md_xv_set1_ps(M.elem[3][1]), w);
    const md_xv tz = md_xv_mul_ps(md_xv_set1_ps(M.elem[3][2]), w);

    size_t i = 0;
    for (; i + MD_XV_WIDTH <= count; i += MD_XV_WIDTH) {
        md_xv x, y, z;
        md_xv_load_xyz_packed_ps(&x, &y, &z, (const float*)(in_xyz + i));

        md_xv rx = md_xv_fmadd_ps(m11, x, tx);
        md_xv ry = md_xv_fmadd_ps(m12, x, ty);
        md_xv rz = md_xv_fmadd_ps(m13, x, tz);

        rx = md_xv_fmadd_ps(m21, y, rx);
        ry = md_xv_fmadd_ps(m22, y, ry);
        rz = md_xv_fmadd_ps(m23, y, rz);

        rx = md_xv_fmadd_ps(m31, z, rx);
        ry = md_xv_fmadd_ps(m32, z, ry);
        rz = md_xv_fmadd_ps(m33, z, rz);

        md_xv_store_xyz_packed_ps((float*)(out_xyz + i), rx, ry, rz);
    }

    for (; i < count; i++) {
        const vec3_t p = in_xyz[i];
        out_xyz[i] = (vec3_t) {
            p.x * M.elem[0][0] + p.y * M.elem[1][0] + p.z * M.elem[2][0] + w_comp * M.elem[3][0],
            p.x * M.elem[0][1] + p.y * M.elem[1][1] + p.z * M.elem[2][1] + w_comp * M.elem[3][1],
            p.x * M.elem[0][2] + p.y * M.elem[1][2] + p.z * M.elem[2][2] + w_comp * M.elem[3][2],
        };
    }
}

void mat4_batch_transform_inplace(vec3_t* in_out_xyz, float w_comp, size_t count, mat4_t M) {
    mat4_batch_transform(in_out_xyz, in_out_xyz, w_comp, count, M);
}

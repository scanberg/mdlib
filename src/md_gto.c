#include <md_gto.h>

#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

#include <stdbool.h>
#include <float.h>

#include <gto_shaders.inl>

#include <GL/gl3w.h>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

static GLuint get_gto_eval_program(void) {
	static GLuint eval_gto_program = 0;
	if (!eval_gto_program) {
		char buf[1024];
		GLint success;
		const char* c_src[] = {(const char*)eval_gto_comp};

		GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
		glShaderSource(shader, ARRAY_SIZE(c_src), c_src, 0);
		glCompileShader(shader);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
		if (!success) {
			int len = 0;
			glGetShaderInfoLog(shader, ARRAY_SIZE(buf), &len, buf);
			MD_LOG_ERROR("Shader compile error:\n%.*s\n", len, buf);
			return 0;
		}

		GLuint program = glCreateProgram();
		glAttachShader(program, shader);
		glLinkProgram(program);
		glGetProgramiv(program, GL_LINK_STATUS, &success);
		if (!success) {
			int len = 0;
			glGetProgramInfoLog(program, ARRAY_SIZE(buf), &len, buf);
			MD_LOG_ERROR("Program link error:\n%.*s\n", len, buf);
			return 0;
		}

		glDeleteShader(shader);
		eval_gto_program = program;
	}
	return eval_gto_program;
}

static GLuint get_gto_eval_ssbo(void) {
	static GLuint ssbo = 0;
	if (!ssbo) {
		glCreateBuffers(1, &ssbo);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
		glBufferData(GL_SHADER_STORAGE_BUFFER, MEGABYTES(4), 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}
	return ssbo;
}

void md_gto_grid_evaluate_GPU(uint32_t vol_tex, const int vol_dim[3], const float vol_step[3], const float* world_to_model, const float* index_to_world, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	ASSERT(world_to_model);
	ASSERT(index_to_world);
	ASSERT(vol_dim);
	ASSERT(vol_step);

	PUSH_GPU_SECTION("EVAL GTO")

	GLuint ssbo = get_gto_eval_ssbo();
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
	glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(md_gto_t) * num_gtos, gtos);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssbo);

	glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);

	GLuint program = get_gto_eval_program();
	glUseProgram(program);

	glUniformMatrix4fv(0, 1, GL_FALSE, world_to_model);
	glUniformMatrix4fv(1, 1, GL_FALSE, index_to_world);
	glUniform3fv(2, 1, (const float*)vol_step);
	glUniform1ui(3, (GLuint)num_gtos);
	glUniform1i(4, (GLint)mode);

	glBindImageTexture(0, vol_tex, 0, GL_TRUE, 0, GL_WRITE_ONLY, GL_R16F);

	int num_groups[3] = {
		DIV_UP(vol_dim[0], 8),
		DIV_UP(vol_dim[1], 8),
		DIV_UP(vol_dim[2], 8),
	};

	glDispatchCompute(num_groups[0], num_groups[1], num_groups[2]);
	glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);

	glUseProgram(0);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	POP_GPU_SECTION()
}

static inline float fast_powf(float base, int exp) {
	float val = 1.0f;
	switch(exp) {
	case 4: val *= base;
	case 3: val *= base;
	case 2: val *= base;
	case 1: val *= base;
	case 0: break;
	}
	return val;
}

static inline double fast_pow(double base, int exp){
	double val = 1.0;
	switch(exp) {
	case 4: val *= base;
	case 3: val *= base;
	case 2: val *= base;
	case 1: val *= base;
	case 0: break;
	}
	return val;
}

static inline md_128 md_mm_fast_pow(md_128 base, md_128i exp) {
	md_128 base2 = md_mm_mul_ps(base,  base);
	md_128 base3 = md_mm_mul_ps(base2, base);
	md_128 base4 = md_mm_mul_ps(base2, base2);

	md_128 mask1 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(1)));
	md_128 mask2 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(2)));
	md_128 mask3 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(3)));
	md_128 mask4 = md_mm_castsi128_ps(md_mm_cmpeq_epi32(exp, md_mm_set1_epi32(4)));

	md_128 res = md_mm_set1_ps(1.0f);
	res = md_mm_blendv_ps(res, base, mask1);
	res = md_mm_blendv_ps(res, base2, mask2);
	res = md_mm_blendv_ps(res, base3, mask3);
	res = md_mm_blendv_ps(res, base4, mask4);
	return res;
}

FORCE_INLINE md_256 md_mm256_fast_pow(md_256 base1, md_256i exp) {
	md_256 base2 = md_mm256_mul_ps(base1, base1);
	md_256 base3 = md_mm256_mul_ps(base2, base1);
	md_256 base4 = md_mm256_mul_ps(base2, base2);

	md_256 mask1 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(1)));
	md_256 mask2 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(2)));
	md_256 mask3 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(3)));
	md_256 mask4 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(4)));

	md_256 res = md_mm256_set1_ps(1.0f);
	res = md_mm256_blendv_ps(res, base1, mask1);
	res = md_mm256_blendv_ps(res, base2, mask2);
	res = md_mm256_blendv_ps(res, base3, mask3);
	res = md_mm256_blendv_ps(res, base4, mask4);
	return res;
}

#ifdef __AVX512F__
static inline __m512 md_mm512_fast_pow(__m512 base, __m512i exp) {
	__m512 base2 = _mm512_mul_ps(base,  base);
	__m512 base3 = _mm512_mul_ps(base2, base);
	__m512 base4 = _mm512_mul_ps(base2, base2);

	__mmask16 mask1 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(1), _MM_CMPINT_EQ);
	__mmask16 mask2 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(2), _MM_CMPINT_EQ);
	__mmask16 mask3 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(3), _MM_CMPINT_EQ);
	__mmask16 mask4 = _mm512_cmp_epi32_mask(exp, _mm512_set1_epi32(4), _MM_CMPINT_EQ);

	__m512 res = _mm512_set1_ps(1.0f);
	res = _mm512_mask_blend_ps(mask1, res, base);
	res = _mm512_mask_blend_ps(mask2, res, base2);
	res = _mm512_mask_blend_ps(mask3, res, base3);
	res = _mm512_mask_blend_ps(mask4, res, base4);
	return res;
}
#endif

static inline void evaluate_grid_ref(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
		const int z_stride = iz * grid_dim[0] * grid_dim[1];
		for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
			const int y_stride = iy * grid_dim[0];
			for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
				const int x_stride = ix;

				float x = grid_origin[0] + ix * grid_step_x[0] + iy * grid_step_y[0] + iz * grid_step_z[0];
				float y = grid_origin[1] + ix * grid_step_x[1] + iy * grid_step_y[1] + iz * grid_step_z[1];
				float z = grid_origin[2] + ix * grid_step_x[2] + iy * grid_step_y[2] + iz * grid_step_z[2];

				double psi = 0.0;
				for (size_t i = 0; i < num_gtos; ++i) {
					float px	= gtos[i].x;
					float py	= gtos[i].y;
					float pz	= gtos[i].z;
					float alpha	= gtos[i].alpha;
					float coeff	= gtos[i].coeff;
					int   pi	= gtos[i].i;
					int   pj	= gtos[i].j;
					int   pk	= gtos[i].k;

					float dx = x - px;
					float dy = y - py;
					float dz = z - pz;
					float d2 = dx * dx + dy * dy + dz * dz;
					float fx = powf(dx, (float)pi);
					float fy = powf(dy, (float)pj);
					float fz = powf(dz, (float)pk);
					float exp_term = (alpha == 0.0f) ? 1.0f : expf(-alpha * d2);
					float powxyz = fx * fy * fz;
					float prod = coeff * powxyz * exp_term;
					psi += prod;
				}

				if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
					psi *= psi;
				}

				int index = x_stride + y_stride + z_stride;
				grid_data[index] = (float)psi;
			}
		}
	}
}

#if 0
// These are old code-paths that were vectorized over the PGTOs
// It is not ideal as the occupation of vector registers is quite poor due to the large variation
// Of how many PGTOs that influence a spatial region

#ifdef __AVX512F__
static inline void evaluate_grid_512(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
	for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
		float z = grid_origin[2] + iz * grid_stepsize[2];
		__m512 vz = _mm512_set1_ps(z);
		int z_stride = iz * grid_dim[0] * grid_dim[1];
		for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
			float y = grid_origin[1] + iy * grid_stepsize[1];
			__m512 vy = _mm512_set1_ps(y);
			int y_stride = iy * grid_dim[0];
			for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
				float x = grid_origin[0] + ix * grid_stepsize[0];
				__m512 vx = _mm512_set1_ps(x);
				int x_stride = ix;

				__m512 vpsi = _mm512_setzero_ps();
				for (size_t i = 0; i < gto->count; i += 16) {
					__m512  px = _mm512_loadu_ps(gto->x + i);
					__m512  py = _mm512_loadu_ps(gto->y + i);
					__m512  pz = _mm512_loadu_ps(gto->z + i);
					__m512  pa = _mm512_loadu_ps(gto->neg_alpha + i);
					__m512  pc = _mm512_loadu_ps(gto->coeff + i);
					__m512i pi = _mm512_loadu_si512(gto->i + i);
					__m512i pj = _mm512_loadu_si512(gto->j + i);
					__m512i pk = _mm512_loadu_si512(gto->k + i);

					__m512 dx = _mm512_sub_ps(vx, px);
					__m512 dy = _mm512_sub_ps(vy, py);
					__m512 dz = _mm512_sub_ps(vz, pz);
					__m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
					__m512 fx = md_mm512_fast_pow(dx, pi);
					__m512 fy = md_mm512_fast_pow(dy, pj);
					__m512 fz = md_mm512_fast_pow(dz, pk);
					__m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));

					__m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);

					if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
						prod = _mm512_mul_ps(prod, prod);
					}

					vpsi = _mm512_add_ps(vpsi, prod);
				}

				float psi = _mm512_reduce_add_ps(vpsi);
				int index = x_stride + y_stride + z_stride;
				grid_data[index] = psi;
			}
		}
	}
}
#endif

static inline void evaluate_grid_256(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
	for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
		float z = grid_origin[2] + iz * grid_stepsize[2];
		md_256 vz = md_mm256_set1_ps(z);
		int z_stride = iz * grid_dim[0] * grid_dim[1];
		for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
			float y = grid_origin[1] + iy * grid_stepsize[1];
			md_256 vy = md_mm256_set1_ps(y);
			int y_stride = iy * grid_dim[0];
			for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
				float x = grid_origin[0] + ix * grid_stepsize[0];
				md_256 vx = md_mm256_set1_ps(x);
				int x_stride = ix;

				md_256 vpsi = md_mm256_setzero_ps();
				for (size_t i = 0; i < gto->count; i += 8) {
					md_256  px = md_mm256_loadu_ps(gto->x + i);
					md_256  py = md_mm256_loadu_ps(gto->y + i);
					md_256  pz = md_mm256_loadu_ps(gto->z + i);
					md_256  pa = md_mm256_loadu_ps(gto->neg_alpha + i);
					md_256  pc = md_mm256_loadu_ps(gto->coeff + i);
					md_256i pi = md_mm256_loadu_si256((const md_256i*)(gto->i + i));
					md_256i pj = md_mm256_loadu_si256((const md_256i*)(gto->j + i));
					md_256i pk = md_mm256_loadu_si256((const md_256i*)(gto->k + i));

					md_256 dx = md_mm256_sub_ps(vx, px);
					md_256 dy = md_mm256_sub_ps(vy, py);
					md_256 dz = md_mm256_sub_ps(vz, pz);
					md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
					md_256 fx = md_mm256_fast_pow(dx, pi);
					md_256 fy = md_mm256_fast_pow(dy, pj);
					md_256 fz = md_mm256_fast_pow(dz, pk);
					md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));

					md_256 prod = md_mm256_mul_ps(md_mm256_mul_ps(md_mm256_mul_ps(pc, fx), md_mm256_mul_ps(fy, fz)), ex);
					if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
						prod = _mm256_mul_ps(prod, prod);
					}

					vpsi = md_mm256_add_ps(vpsi, prod);
				}

				float psi = md_mm256_reduce_add_ps(vpsi);
				int index = x_stride + y_stride + z_stride;
				grid_data[index] = psi;
			}
		}
	}
}

static inline void evaluate_grid_128(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto, md_gto_eval_mode_t mode) {
	for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
		float z = grid_origin[2] + iz * grid_stepsize[2];
		md_128 vz = md_mm_set1_ps(z);
		int z_stride = iz * grid_dim[0] * grid_dim[1];
		for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
			float y = grid_origin[1] + iy * grid_stepsize[1];
			md_128 vy = md_mm_set1_ps(y);
			int y_stride = iy * grid_dim[0];
			for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
				float x = grid_origin[0] + ix * grid_stepsize[0];
				md_128 vx = md_mm_set1_ps(x);
				int x_stride = ix;

				md_128 vpsi = md_mm_setzero_ps();
				for (size_t i = 0; i < gto->count; i += 4) {
					md_128  px = md_mm_loadu_ps(gto->x + i);
					md_128  py = md_mm_loadu_ps(gto->y + i);
					md_128  pz = md_mm_loadu_ps(gto->z + i);
					md_128  pa = md_mm_loadu_ps(gto->neg_alpha + i);
					md_128  pc = md_mm_loadu_ps(gto->coeff + i);
					md_128i pi = md_mm_loadu_si128((const md_128i*)(gto->i + i));
					md_128i pj = md_mm_loadu_si128((const md_128i*)(gto->j + i));
					md_128i pk = md_mm_loadu_si128((const md_128i*)(gto->k + i));

					md_128 dx = md_mm_sub_ps(vx, px);
					md_128 dy = md_mm_sub_ps(vy, py);
					md_128 dz = md_mm_sub_ps(vz, pz);
					md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
					md_128 fx = md_mm_fast_pow(dx, pi);
					md_128 fy = md_mm_fast_pow(dy, pj);
					md_128 fz = md_mm_fast_pow(dz, pk);
					md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));

					md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);
					if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
						prod = _mm_mul_ps(prod, prod);
					}

					vpsi = md_mm_add_ps(vpsi, prod);
				}

				float psi = md_mm_reduce_add_ps(vpsi);
				int index = x_stride + y_stride + z_stride;
				grid_data[index] = psi;
			}
		}
	}
}
#endif

#if defined(__AVX512F__) && defined(__AVX512DQ__)

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_512(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
	const md_256  vxh = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));
	const __m512  vx  = _mm512_insertf32x8(_mm512_castps256_ps512(vxh), vxh, 1);
	const int x_stride = grid_idx_min[0];

	__m512 vpsi[8][4] = {0};

	for (size_t i = 0; i < num_gtos; ++i) {
		const __m512  px = _mm512_set1_ps(gtos[i].x);
		const __m512  py = _mm512_set1_ps(gtos[i].y);
		const __m512  pz = _mm512_set1_ps(gtos[i].z);
		const __m512  pc = _mm512_set1_ps(gtos[i].coeff);
		const __m512  pa = _mm512_set1_ps(-gtos[i].alpha); // Negate alpha here
		const __m512i pi = _mm512_set1_epi32(gtos[i].i);
		const __m512i pj = _mm512_set1_epi32(gtos[i].j);
		const __m512i pk = _mm512_set1_epi32(gtos[i].k);

		for (int iz = 0; iz < 8; ++iz) {
			float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
			__m512 vz = _mm512_set1_ps(z);
			for (int iy = 0; iy < 4; ++iy) {
				float y[2] = {
					grid_origin[1] + (grid_idx_min[1] + iy * 2 + 0) * grid_step[1],
					grid_origin[1] + (grid_idx_min[1] + iy * 2 + 1) * grid_step[1],
				};
				__m512 vy = _mm512_insertf32x8(_mm512_set1_ps(y[0]), _mm256_set1_ps(y[1]), 1);

				__m512 dx = _mm512_sub_ps(vx, px);
				__m512 dy = _mm512_sub_ps(vy, py);
				__m512 dz = _mm512_sub_ps(vz, pz);
				__m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
				__m512 fx = md_mm512_fast_pow(dx, pi);
				__m512 fy = md_mm512_fast_pow(dy, pj);
				__m512 fz = md_mm512_fast_pow(dz, pk);
				__m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));
				__m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);

				vpsi[iz][iy] = _mm512_add_ps(vpsi[iz][iy], prod);
			}
		}
	}

	// Write result block to memory
	for (int iz = 0; iz < 8; ++iz) {
		int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
		for (int iy = 0; iy < 4; ++iy) {
			int y_stride[2] = {
				(grid_idx_min[1] + iy * 2 + 0) * grid_dim[0],
				(grid_idx_min[1] + iy * 2 + 1) * grid_dim[0],
			};
			int index[2] = {
				x_stride + y_stride[0] + z_stride,
				x_stride + y_stride[1] + z_stride,
			};

			md_512 psi = vpsi[iz][iy];
			if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
				psi = _mm512_mul_ps(psi, psi);
			}

			md_256 tpsi[2] = {
				_mm512_castps512_ps256(psi),
				_mm512_extractf32x8_ps(psi, 1),
			};

			md_mm256_storeu_ps(grid_data + index[0], tpsi[0]);
			md_mm256_storeu_ps(grid_data + index[1], tpsi[1]);
		}
	}
}

#endif

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	const md_256i vix  = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
	const int x_stride = grid_idx_min[0];
	const md_256   vx  = md_mm256_fmadd_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step[0]), md_mm256_set1_ps(grid_origin[0]));

	// Operate on local block to avoid cache-line contention across threads
	md_256 vpsi[8][8] = {0};

	for (size_t i = 0; i < num_gtos; ++i) {
		const md_256  px = md_mm256_set1_ps(gtos[i].x);
		const md_256  py = md_mm256_set1_ps(gtos[i].y);
		const md_256  pz = md_mm256_set1_ps(gtos[i].z);
		const md_256  pc = md_mm256_set1_ps(gtos[i].coeff);
		const md_256  pa = md_mm256_set1_ps(-gtos[i].alpha); // Negate alpha here
		const md_256i pi = md_mm256_set1_epi32(gtos[i].i);
		const md_256i pj = md_mm256_set1_epi32(gtos[i].j);
		const md_256i pk = md_mm256_set1_epi32(gtos[i].k);

		for (int iz = 0; iz < 8; ++iz) {
			float z = grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2];
			md_256 vz = md_mm256_set1_ps(z);

			for (int iy = 0; iy < 8; ++iy) {
				float y = grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1];
				md_256 vy = md_mm256_set1_ps(y);

				md_256 dx = md_mm256_sub_ps(vx, px);
				md_256 dy = md_mm256_sub_ps(vy, py);
				md_256 dz = md_mm256_sub_ps(vz, pz);
				md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
				md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));
				md_256 fx = md_mm256_fast_pow(dx, pi);
				md_256 fy = md_mm256_fast_pow(dy, pj);
				md_256 fz = md_mm256_fast_pow(dz, pk);

				md_256 prod_a = md_mm256_mul_ps(pc, fx);
				md_256 prod_b = md_mm256_mul_ps(fy, fz);
				vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(prod_a, prod_b), ex, vpsi[iz][iy]);
			}
		}
	}

	// Write result block to memory
	for (int iz = 0; iz < 8; ++iz) {
		int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
		for (int iy = 0; iy < 8; ++iy) {
			int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
			int index = x_stride + y_stride + z_stride;
			md_256 psi = vpsi[iz][iy];

			if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
				psi = md_mm256_mul_ps(psi, psi);
			}

			md_mm256_storeu_ps(grid_data + index, psi);
		}
	}
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_256(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	const md_256i vix = md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[0]), md_mm256_set_epi32(7,6,5,4,3,2,1,0));
	const int x_stride = grid_idx_min[0];

	const md_256 gsx[3] = {
		md_mm256_add_ps(md_mm256_set1_ps(grid_origin[0]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[0]))),
		md_mm256_add_ps(md_mm256_set1_ps(grid_origin[1]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[1]))),
		md_mm256_add_ps(md_mm256_set1_ps(grid_origin[2]), md_mm256_mul_ps(md_mm256_cvtepi32_ps(vix), md_mm256_set1_ps(grid_step_x[2]))),
	};

	const md_256 gsz[3] = {
		md_mm256_set1_ps(grid_step_z[0]),
		md_mm256_set1_ps(grid_step_z[1]),
		md_mm256_set1_ps(grid_step_z[2]),
	};

	const md_256 gsy[3] = {
		md_mm256_set1_ps(grid_step_y[0]),
		md_mm256_set1_ps(grid_step_y[1]),
		md_mm256_set1_ps(grid_step_y[2]),
	};

	// Operate on local block to avoid cache-line contention across threads
	md_256 vpsi[8][8] = {0};

	for (size_t i = 0; i < num_gtos; ++i) {
		const md_256  px = md_mm256_set1_ps(gtos[i].x);
		const md_256  py = md_mm256_set1_ps(gtos[i].y);
		const md_256  pz = md_mm256_set1_ps(gtos[i].z);
		const md_256  pc = md_mm256_set1_ps(gtos[i].coeff);
		const md_256  pa = md_mm256_set1_ps(-gtos[i].alpha); // Negate alpha here
		const md_256i pi = md_mm256_set1_epi32(gtos[i].i);
		const md_256i pj = md_mm256_set1_epi32(gtos[i].j);
		const md_256i pk = md_mm256_set1_epi32(gtos[i].k);
	
		for (int iz = 0; iz < 8; ++iz) {
			const md_256 tz = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[2]), md_mm256_set1_epi32(iz)));

			const md_256 xz[3] = {
				md_mm256_fmadd_ps(tz, gsz[0], gsx[0]),
				md_mm256_fmadd_ps(tz, gsz[1], gsx[1]),
				md_mm256_fmadd_ps(tz, gsz[2], gsx[2]),
			};

			for (int iy = 0; iy < 8; ++iy) {
				const md_256 ty = md_mm256_cvtepi32_ps(md_mm256_add_epi32(md_mm256_set1_epi32(grid_idx_min[1]), md_mm256_set1_epi32(iy)));

				md_256 vx = md_mm256_fmadd_ps(ty, gsy[0], xz[0]);
				md_256 vy = md_mm256_fmadd_ps(ty, gsy[1], xz[1]);
				md_256 vz = md_mm256_fmadd_ps(ty, gsy[2], xz[2]);

				md_256 dx = md_mm256_sub_ps(vx, px);
				md_256 dy = md_mm256_sub_ps(vy, py);
				md_256 dz = md_mm256_sub_ps(vz, pz);
				md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
				md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));
				md_256 fx = md_mm256_fast_pow(dx, pi);
				md_256 fy = md_mm256_fast_pow(dy, pj);
				md_256 fz = md_mm256_fast_pow(dz, pk);

				md_256 prod_a = md_mm256_mul_ps(pc, fx);
				md_256 prod_b = md_mm256_mul_ps(fy, fz);
				vpsi[iz][iy] = md_mm256_fmadd_ps(md_mm256_mul_ps(prod_a, prod_b), ex, vpsi[iz][iy]);
			}
		}
	}

	// Write result block to memory
	for (int iz = 0; iz < 8; ++iz) {
		int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
		for (int iy = 0; iy < 8; ++iy) {
			int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
			int index = x_stride + y_stride + z_stride;
			md_256 psi = vpsi[iz][iy];

			if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
				psi = md_mm256_mul_ps(psi, psi);
			}

			md_mm256_storeu_ps(grid_data + index, psi);
		}
	}
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_ortho_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	const md_128i vix[2] = {
		md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 0), md_mm_set_epi32(3,2,1,0)),
		md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
	};
	const md_128   vx[2] = {
		md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
		md_mm_fmadd_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step[0]), md_mm_set1_ps(grid_origin[0])),
	};
	const int x_stride[2] = {
		grid_idx_min[0] + 0,
		grid_idx_min[0] + 4,
	};

	// Operate on local block to avoid cache-line contention across threads
	md_128 vpsi[8][8][2] = {0};

	for (size_t i = 0; i < num_gtos; ++i) {
		const md_128  px = md_mm_set1_ps(gtos[i].x);
		const md_128  py = md_mm_set1_ps(gtos[i].y);
		const md_128  pz = md_mm_set1_ps(gtos[i].z);
		const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
		const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
		const md_128i pi = md_mm_set1_epi32(gtos[i].i);
		const md_128i pj = md_mm_set1_epi32(gtos[i].j);
		const md_128i pk = md_mm_set1_epi32(gtos[i].k);

		for (int iz = 0; iz < 8; ++iz) {
			md_128 vz = md_mm_set1_ps(grid_origin[2] + (grid_idx_min[2] + iz) * grid_step[2]);
			for (int iy = 0; iy < 8; ++iy) {
				md_128 vy = md_mm_set1_ps(grid_origin[1] + (grid_idx_min[1] + iy) * grid_step[1]);

				for (int ix = 0; ix < 2; ++ix) {
					md_128 dx = md_mm_sub_ps(vx[ix], px);
					md_128 dy = md_mm_sub_ps(vy,	 py);	 
					md_128 dz = md_mm_sub_ps(vz,	 pz);	 
					md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
					md_128 fx = md_mm_fast_pow(dx, pi);
					md_128 fy = md_mm_fast_pow(dy, pj);
					md_128 fz = md_mm_fast_pow(dz, pk);
					md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));
					md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

					vpsi[iz][iy][ix] = md_mm_add_ps(vpsi[iz][iy][ix], prod);
				}
			}
		}
	}

	// Write result block to memory
	for (int iz = 0; iz < 8; ++iz) {
		int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
		for (int iy = 0; iy < 8; ++iy) {
			int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
			int index[2] = {
				x_stride[0] + y_stride + z_stride,
				x_stride[1] + y_stride + z_stride,
			};

			md_128 psi[2] = {
				vpsi[iz][iy][0],
				vpsi[iz][iy][1],
			};

			if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
				psi[0] = md_mm_mul_ps(psi[0], psi[0]);
				psi[1] = md_mm_mul_ps(psi[1], psi[1]);
			}

			md_mm_storeu_ps(grid_data + index[0], psi[0]);
			md_mm_storeu_ps(grid_data + index[1], psi[1]);
		}
	}
}

// Evaluate 8 voxels per gto
static inline void evaluate_grid_8x8x8_128(float grid_data[], const int grid_idx_min[3], const int grid_dim[3], const float grid_origin[3], const float grid_step_x[3], const float grid_step_y[3], const float grid_step_z[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	const int x_stride[2] = {
		grid_idx_min[0] + 0,
		grid_idx_min[0] + 4,
	};

	const md_128i vix[2] = {
		md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0])    , md_mm_set_epi32(3,2,1,0)),
		md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[0] + 4), md_mm_set_epi32(3,2,1,0)),
	};

	const md_128 gsx[2][3] = {
		{
			md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[0]))),
			md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[1]))),
			md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[0]), md_mm_set1_ps(grid_step_x[2]))),
		},
		{
			md_mm_add_ps(md_mm_set1_ps(grid_origin[0]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[0]))),
			md_mm_add_ps(md_mm_set1_ps(grid_origin[1]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[1]))),
			md_mm_add_ps(md_mm_set1_ps(grid_origin[2]), md_mm_mul_ps(md_mm_cvtepi32_ps(vix[1]), md_mm_set1_ps(grid_step_x[2]))),
		},
	};

	const md_128 gsz[3] = {
		md_mm_set1_ps(grid_step_z[0]),
		md_mm_set1_ps(grid_step_z[1]),
		md_mm_set1_ps(grid_step_z[2]),
	};

	const md_128 gsy[3] = {
		md_mm_set1_ps(grid_step_y[0]),
		md_mm_set1_ps(grid_step_y[1]),
		md_mm_set1_ps(grid_step_y[2]),
	};

	// Operate on local block to avoid cache-line contention across threads
	md_128 vpsi[8][8][2] = {0};

	for (size_t i = 0; i < num_gtos; ++i) {
		const md_128  px = md_mm_set1_ps(gtos[i].x);
		const md_128  py = md_mm_set1_ps(gtos[i].y);
		const md_128  pz = md_mm_set1_ps(gtos[i].z);
		const md_128  pc = md_mm_set1_ps(gtos[i].coeff);
		const md_128  pa = md_mm_set1_ps(-gtos[i].alpha); // Negate alpha here
		const md_128i pi = md_mm_set1_epi32(gtos[i].i);
		const md_128i pj = md_mm_set1_epi32(gtos[i].j);
		const md_128i pk = md_mm_set1_epi32(gtos[i].k);

		for (int iz = 0; iz < 8; ++iz) {
			const md_128 tz = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[2]), md_mm_set1_epi32(iz)));
			const md_128 xz[2][3] = {
				{
					md_mm_fmadd_ps(tz, gsz[0], gsx[0][0]),
					md_mm_fmadd_ps(tz, gsz[1], gsx[0][1]),
					md_mm_fmadd_ps(tz, gsz[2], gsx[0][2]),
				},
				{
					md_mm_fmadd_ps(tz, gsz[0], gsx[1][0]),
					md_mm_fmadd_ps(tz, gsz[1], gsx[1][1]),
					md_mm_fmadd_ps(tz, gsz[2], gsx[1][2]),
				},
			};
			for (int iy = 0; iy < 8; ++iy) {
				const md_128 ty = md_mm_cvtepi32_ps(md_mm_add_epi32(md_mm_set1_epi32(grid_idx_min[1]), md_mm_set1_epi32(iy)));

				const md_128 vx[2] = {
					md_mm_fmadd_ps(ty, gsy[0], xz[0][0]),
					md_mm_fmadd_ps(ty, gsy[0], xz[1][0]),
				};
				const md_128 vy[2] = {
					md_mm_fmadd_ps(ty, gsy[1], xz[0][1]),
					md_mm_fmadd_ps(ty, gsy[1], xz[1][1]),
				};
				const md_128 vz[2] = {
					md_mm_fmadd_ps(ty, gsy[2], xz[0][2]),
					md_mm_fmadd_ps(ty, gsy[2], xz[1][2]),
				};

				for (int ix = 0; ix < 2; ++ix) {
					md_128 dx = md_mm_sub_ps(vx[ix], px);
					md_128 dy = md_mm_sub_ps(vy[ix], py);
					md_128 dz = md_mm_sub_ps(vz[ix], pz);
					md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
					md_128 fx = md_mm_fast_pow(dx, pi);
					md_128 fy = md_mm_fast_pow(dy, pj);
					md_128 fz = md_mm_fast_pow(dz, pk);
					md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));
					md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

					vpsi[iz][iy][ix] = md_mm_add_ps(vpsi[iz][iy][ix], prod);
				}
			}
		}
	}

	// Write result block to memory
	for (int iz = 0; iz < 8; ++iz) {
		int z_stride = (grid_idx_min[2] + iz) * grid_dim[0] * grid_dim[1];
		for (int iy = 0; iy < 8; ++iy) {
			int y_stride = (grid_idx_min[1] + iy) * grid_dim[0];
			int index[2] = {
				x_stride[0] + y_stride + z_stride,
				x_stride[1] + y_stride + z_stride,
			};

			md_128 psi[2] = {
				vpsi[iz][iy][0],
				vpsi[iz][iy][1],
			};

			if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
				psi[0] = md_mm_mul_ps(psi[0], psi[0]);
				psi[1] = md_mm_mul_ps(psi[1], psi[1]);
			}

			md_mm_storeu_ps(grid_data + index[0], psi[0]);
			md_mm_storeu_ps(grid_data + index[1], psi[1]);
		}
	}
}

void md_gto_grid_evaluate_sub(md_grid_t* grid, const int grid_idx_off[3], const int grid_idx_len[3], const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	ASSERT(grid);
	ASSERT(gtos);

	const int* grid_idx_min = grid_idx_off;
	const int  grid_idx_max[3] = {
		grid_idx_off[0] + grid_idx_len[0],
		grid_idx_off[1] + grid_idx_len[1],
		grid_idx_off[2] + grid_idx_len[2],
	};

	//printf("Number of pgtos in volume region: %zu\n", gto.count);

	bool ortho =
		(grid->step_x[1] == 0 && grid->step_x[2] == 0) &&
		(grid->step_y[0] == 0 && grid->step_y[2] == 0) &&
		(grid->step_z[0] == 0 && grid->step_z[1] == 0);

	float ortho_step[3] = {grid->step_x[0], grid->step_y[1], grid->step_z[2]};

	// There are specialized versions for evaluating 8x8x8 subgrids
	// 8x8x8 Is a good chunk size to operate on as it probably fits in L1 Cache together with the GTOs
	// Then we vectorize over the spatial domain rather than the GTOs to get better register occupation
	if (grid_idx_len[0] == 8 && grid_idx_len[1] == 8 && grid_idx_len[2] == 8) {

#if defined(__AVX512F__) && defined(__AVX512DQ__) || defined (__AVX2__)
		// @TODO: Implement real AVX512 path
		if (ortho) {
			evaluate_grid_ortho_8x8x8_256(grid->data, grid_idx_min, grid->dim, grid->origin, ortho_step, gtos, num_gtos, mode);
		} else {
			evaluate_grid_8x8x8_256(grid->data, grid_idx_min, grid->dim, grid->origin, grid->step_x, grid->step_y, grid->step_z, gtos, num_gtos, mode);
		}
#elif defined(__SSE2__)
		if (ortho) {
			evaluate_grid_ortho_8x8x8_128(grid->data, grid_idx_min, grid->dim, grid->origin, ortho_step, gtos, num_gtos, mode);
		} else {
			evaluate_grid_8x8x8_128(grid->data, grid_idx_min, grid->dim, grid->origin, grid->step_x, grid->step_y, grid->step_z, gtos, num_gtos, mode);
		}
#else
		evaluate_grid_ref(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->step_x, grid->step_y, grid->step_z, gtos, num_gtos, mode);
#endif
	} else {
		// Slowpath
		evaluate_grid_ref(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->step_x, grid->step_y, grid->step_z, gtos, num_gtos, mode);
	}
}

void md_gto_grid_evaluate(md_grid_t* grid, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	ASSERT(grid);
	ASSERT(gtos);

	int idx_off[3] = {0};
	int idx_len[3] = {0};

	size_t temp_pos = md_temp_get_pos();
	md_gto_t* sub_gtos = (md_gto_t*)md_temp_push(sizeof(md_gto_t) * num_gtos);

	for (idx_off[2] = 0; idx_off[2] < grid->dim[2]; idx_off[2] += 8) {
		idx_len[2] = MIN(8, grid->dim[2] - idx_off[2]);
		for (idx_off[1] = 0; idx_off[1] < grid->dim[1]; idx_off[1] += 8) {
			idx_len[1] = MIN(8, grid->dim[1] - idx_off[1]);
			for (idx_off[0] = 0; idx_off[0] < grid->dim[0]; idx_off[0] += 8) {
				idx_len[0] = MIN(8, grid->dim[0] - idx_off[0]);

				float aabb_min[3] = {
					grid->origin[0] + idx_off[0] * (grid->step_x[0] + grid->step_y[0] + grid->step_z[0]),
					grid->origin[1] + idx_off[1] * (grid->step_x[1] + grid->step_y[1] + grid->step_z[1]),
					grid->origin[2] + idx_off[2] * (grid->step_x[2] + grid->step_y[2] + grid->step_z[2]),
				};
				float aabb_max[3] = {
					grid->origin[0] + (idx_off[0] + idx_len[0]) * (grid->step_x[0] + grid->step_y[0] + grid->step_z[0]),
					grid->origin[1] + (idx_off[1] + idx_len[1]) * (grid->step_x[1] + grid->step_y[1] + grid->step_z[1]),
					grid->origin[2] + (idx_off[2] + idx_len[2]) * (grid->step_x[2] + grid->step_y[2] + grid->step_z[2]),
				};

				size_t num_sub_gtos = md_gto_aabb_test(sub_gtos, aabb_min, aabb_max, gtos, num_gtos);
				md_gto_grid_evaluate_sub(grid, idx_off, idx_len, sub_gtos, num_sub_gtos, mode);
			}
		}
	}

	md_temp_set_pos_back(temp_pos);
}

// Evaluate GTOs over a set of passed in packed XYZ coordinates with a bytestride
static void evaluate_gtos(float* out_psi, const float* in_xyz, size_t num_xyz, size_t xyz_stride, const md_gto_t* in_gto, size_t num_gtos, md_gto_eval_mode_t mode) {
	for (size_t j = 0; j < num_xyz; ++j) {
		const float* xyz = (const float*)((const char*)in_xyz + j * xyz_stride);
		double x = xyz[0];
		double y = xyz[1];
		double z = xyz[2];

		double psi = 0.0;
		for (size_t i = 0; i < num_gtos; ++i) {
			double cutoff	= in_gto[i].cutoff;
			double rx		= x - in_gto[i].x;
			double ry		= y - in_gto[i].y;
			double rz		= z - in_gto[i].z;
			double r2		= rx * rx + ry * ry + rz * rz;
			if (r2 > cutoff * cutoff) {
				continue;
			}

			double alpha	= in_gto[i].alpha;
			double coeff	= in_gto[i].coeff;
			int   pi		= in_gto[i].i;
			int   pj		= in_gto[i].j;
			int   pk		= in_gto[i].k;

			double fx = pow(rx, pi);
			double fy = pow(ry, pj);
			double fz = pow(rz, pk);
			double powxyz = fx * fy * fz;
			double exp_term = alpha == 0 ? 1.0 : exp(-alpha * r2);

			double prod = coeff * powxyz * exp_term;
			psi += prod;
		}

		if (mode == MD_GTO_EVAL_MODE_PSI_SQUARED) {
			psi = psi * psi;
		}

		out_psi[j] = (float)psi;
	}
}

void md_gto_xyz_evaluate(float* out_psi, const float* in_xyz, size_t num_xyz, size_t stride, const md_gto_t* gtos, size_t num_gtos, md_gto_eval_mode_t mode) {
	if (!out_psi) {
		MD_LOG_ERROR("out_psi array is NULL!");
		return;
	}
	if (!in_xyz) {
		MD_LOG_ERROR("in_xyz base pointer is NULL!");
		return;
	}
	if (!gtos) {
		MD_LOG_ERROR("gtos is NULL!");
		return;
	}
	if (stride != 0 && stride < sizeof(float) * 3) {
		MD_LOG_ERROR("Invalid xyz stride: expected value >= 12 Bytes, got %zu", stride);
		return;
	}

	stride = (stride == 0) ? sizeof(float) * 3 : stride;
	evaluate_gtos(out_psi, in_xyz, num_xyz, stride, gtos, num_gtos, mode);
}

static inline double eval_G(double d, double C, int l, double neg_alpha) {
	return C * fast_pow(d, l) * exp(neg_alpha * d * d);
}

static inline void eval_G_and_G_prime(double* out_G, double* out_G_prime, double d, double C, int l, double neg_alpha) {
	double exp_term = exp(neg_alpha * d * d);
	*out_G		 = C * fast_pow(d, l) * exp_term;
	*out_G_prime = C * fast_pow(d, l-1) * (l + 2 * neg_alpha * d * d) * exp_term;
}

#define PRINT_RESULT 0

static double compute_distance_cutoff(double cutoff_value, int i, int j, int k, int l, double coeff, double alpha) {
	double d = 0.0;

	const double neg_alpha = -alpha;

	// Bake into single constant C
	const double C = fabs(coeff * sqrt((fast_pow(i,i) * fast_pow(j,j) * fast_pow(k,k)) / fast_pow(l,l)));

	// Compute maxima
	const double d_maxima = sqrt(l / (2.0 * fabs(neg_alpha)));

	// Check the contribution at the maxima
	const double y_max = eval_G(d_maxima, C, l, neg_alpha);
	if (y_max < cutoff_value) {
		d = 0.0;
		goto done;
	}

	// If we have an S-type orbital (l == 0) the expression collapses into an expression we can invert and evaluate
	if (l == 0) {
		double y = cutoff_value;
		double a = fabs(coeff) / y;
		double la = log(a);
		d = sqrt(fabs(la) / fabs(neg_alpha));
		goto done;
	}

	// If we end up here we need to perform a numerical search for a d value where the value G(d) < cutoff_value
	// We do not want to overestimate d which will result in a too large radius of influence for the PGTO.
	// And will directly negatively impact performance when evaluating on the grid.
	// Therefore we want to find d where G(d) < cutoff_value but within a tolerance of cutoff_value

	// Search parameters
	const double d_min = d_maxima + 0.001;
	const double d_max = d_maxima + 100.0;
	const double y_tol = cutoff_value * 0.001;
	const double d_tol = 1.0e-9;

	// Initial guess
	// This is the analytical solution for d^2/dx^2 G(d) = 0
	// Which should give us a value which corresponds to the point where we have the maximum negative slope
	d = 0.5 * sqrt(sqrt(neg_alpha*neg_alpha * (8*l + 1)) / (neg_alpha*neg_alpha) + (2*l + 1) / fabs(neg_alpha));

	// Newton-Rhapson iterative search
	for (int iter = 0; iter < 100; ++iter) {
		double y, yp;
		eval_G_and_G_prime(&y, &yp, d, C, l, neg_alpha);

		// Shift function so it intersects the x axis at the point we seek (with a bias towards values less than cutoff value)
		y = y - cutoff_value + y_tol;

		//printf ("d: %.10f, y: %.10f, yp: %.10f\n", d, y, yp);

		if (y < 0 && fabs(y) < y_tol) {
			//printf ("y tolerance met after %i iterations\n", iter);
			break;
		}

		if (fabs(yp) < DBL_EPSILON) {
			//printf ("Denominator is too small!\n");
			break;
		}

		double dn = d - y / yp;
		dn = CLAMP(dn, d_min, d_max);

		if (fabs(dn - d) < d_tol) {
			//printf ("d tolerance met after %i iterations\n", iter);
			break;
		}

		d = dn;
	}

done:
#if PRINT_RESULT
	if (d > 0.0) {
		printf("Cutoff dist and value: %15.5f, %15.12f\n", d, eval_G(d, C, l, neg_alpha));
	}
#endif
	return d;
}

void md_gto_cutoff_compute(md_gto_t* gtos, size_t count, double value) {
	if (value == 0) {
		for (size_t i = 0; i < count; ++i) {
			gtos[i].cutoff = FLT_MAX;
		}
	} else {
		for (size_t i = 0; i < count; ++i) {
			gtos[i].cutoff = (float)compute_distance_cutoff(value, gtos[i].i, gtos[i].j, gtos[i].k, gtos[i].l, gtos[i].coeff, gtos[i].alpha);
		}
	}
}

size_t md_gto_aabb_test(md_gto_t* out_gtos, const float aabb_min[3], const float aabb_max[3], const md_gto_t* in_gtos, size_t in_num_gtos) {
	// Extract a subset of gtos that overlap with the evaluated subportion of the grid
	// @TODO: This can be vectorized as well
	size_t num_gtos = 0;
	for (size_t i = 0; i < in_num_gtos; ++i) {
		float x  = in_gtos[i].x;
		float y  = in_gtos[i].y;
		float z  = in_gtos[i].z;
		float cutoff = in_gtos[i].cutoff;

		float cx = CLAMP(x, aabb_min[0], aabb_max[0]);
		float cy = CLAMP(y, aabb_min[1], aabb_max[1]);
		float cz = CLAMP(z, aabb_min[2], aabb_max[2]);

		float dx = x - cx;
		float dy = y - cy;
		float dz = z - cz;

		float d2 = dx * dx + dy * dy + dz * dz;

		if (d2 > cutoff * cutoff) {
			continue;
		}
		out_gtos[num_gtos++] = in_gtos[i];
	}
	return num_gtos;
}

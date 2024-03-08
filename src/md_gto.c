#include <md_gto.h>

#include <core/md_simd.h>
#include <core/md_allocator.h>
#include <stdbool.h>
#include <float.h>

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

static inline md_256 md_mm256_fast_pow(md_256 base, md_256i exp) {
	md_256 base2 = md_mm256_mul_ps(base,  base);
	md_256 base3 = md_mm256_mul_ps(base2, base);
	md_256 base4 = md_mm256_mul_ps(base2, base2);

	md_256 mask1 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(1)));
	md_256 mask2 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(2)));
	md_256 mask3 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(3)));
	md_256 mask4 = md_mm256_castsi256_ps(md_mm256_cmpeq_epi32(exp, md_mm256_set1_epi32(4)));

	md_256 res = md_mm256_set1_ps(1.0f);
	res = md_mm256_blendv_ps(res, base, mask1);
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

static inline void evaluate_grid_ref(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto) {
	for (int iz = grid_idx_min[2]; iz < grid_idx_max[2]; iz++) {
		float z = grid_origin[2] + iz * grid_stepsize[2];
		int z_stride = iz * grid_dim[0] * grid_dim[1];
		for (int iy = grid_idx_min[1]; iy < grid_idx_max[1]; ++iy) {
			float y = grid_origin[1] + iy * grid_stepsize[1];
			int y_stride = iy * grid_dim[0];
			for (int ix = grid_idx_min[0]; ix < grid_idx_max[0]; ++ix) {
				float x = grid_origin[0] + ix * grid_stepsize[0];
				int x_stride = ix;

				double psi = 0.0;

				for (size_t i = 0; i < gto->count; ++i) {
					float px		= gto->x[i];
					float py		= gto->y[i];
					float pz		= gto->z[i];
					float neg_alpha = gto->neg_alpha[i];
					float coeff		= gto->coeff[i];
					int   pi		= gto->i[i];
					int   pj		= gto->j[i];
					int   pk		= gto->k[i];

					float dx = px - x;
					float dy = py - y;
					float dz = pz - z;
					float d2 = dx * dx + dy * dy + dz * dz;
					float fx = fast_powf(dx, pi);
					float fy = fast_powf(dy, pj);
					float fz = fast_powf(dz, pk);
					float exponentTerm = neg_alpha == 0 ? 1.0f : expf(neg_alpha * d2);
					psi += coeff * fx * fy * fz * exponentTerm;
				}

				int index = x_stride + y_stride + z_stride;
				grid_data[index] = (float)psi;
			}
		}
	}
}

#ifdef __AVX512F__
static inline void evaluate_grid_512(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto) {
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
					__m512  px = _mm512_load_ps(gto->x + i);
					__m512  py = _mm512_load_ps(gto->y + i);
					__m512  pz = _mm512_load_ps(gto->z + i);
					__m512  pa = _mm512_load_ps(gto->neg_alpha + i);
					__m512  pc = _mm512_load_ps(gto->coeff + i);
					__m512i pi = _mm512_load_si512(gto->i + i);
					__m512i pj = _mm512_load_si512(gto->j + i);
					__m512i pk = _mm512_load_si512(gto->k + i);

					__m512 dx = _mm512_sub_ps(px, vx);
					__m512 dy = _mm512_sub_ps(py, vy);
					__m512 dz = _mm512_sub_ps(pz, vz);
					__m512 d2 = _mm512_fmadd_ps(dx, dx, _mm512_fmadd_ps(dy, dy, _mm512_mul_ps(dz, dz)));
					__m512 fx = md_mm512_fast_pow(dx, pi);
					__m512 fy = md_mm512_fast_pow(dy, pj);
					__m512 fz = md_mm512_fast_pow(dz, pk);
					__m512 ex = md_mm512_exp_ps(_mm512_mul_ps(pa, d2));

					__m512 prod = _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(pc, fx), _mm512_mul_ps(fy, fz)), ex);
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

static inline void evaluate_grid_256(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto) {
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
					md_256  px = md_mm256_load_ps(gto->x + i);
					md_256  py = md_mm256_load_ps(gto->y + i);
					md_256  pz = md_mm256_load_ps(gto->z + i);
					md_256  pa = md_mm256_load_ps(gto->neg_alpha + i);
					md_256  pc = md_mm256_load_ps(gto->coeff + i);
					md_256i pi = md_mm256_load_si256((const md_256i*)(gto->i + i));
					md_256i pj = md_mm256_load_si256((const md_256i*)(gto->j + i));
					md_256i pk = md_mm256_load_si256((const md_256i*)(gto->k + i));

					md_256 dx = md_mm256_sub_ps(px, vx);
					md_256 dy = md_mm256_sub_ps(py, vy);
					md_256 dz = md_mm256_sub_ps(pz, vz);
					md_256 d2 = md_mm256_fmadd_ps(dx, dx, md_mm256_fmadd_ps(dy, dy, md_mm256_mul_ps(dz, dz)));
					md_256 fx = md_mm256_fast_pow(dx, pi);
					md_256 fy = md_mm256_fast_pow(dy, pj);
					md_256 fz = md_mm256_fast_pow(dz, pk);
					md_256 ex = md_mm256_exp_ps(md_mm256_mul_ps(pa, d2));

					md_256 prod = md_mm256_mul_ps(md_mm256_mul_ps(md_mm256_mul_ps(pc, fx), md_mm256_mul_ps(fy, fz)), ex);
					vpsi = md_mm256_add_ps(vpsi, prod);
				}

				float psi = md_mm256_reduce_add_ps(vpsi);
				int index = x_stride + y_stride + z_stride;
				grid_data[index] = psi;
			}
		}
	}
}

static inline void evaluate_grid_128(float grid_data[], const int grid_idx_min[3], const int grid_idx_max[3], const int grid_dim[3], const float grid_origin[3], const float grid_stepsize[3], const md_gto_data_t* gto) {
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

				double psi = 0.0;

				md_128 vpsi = md_mm_setzero_ps();
				for (size_t i = 0; i < gto->count; i += 4) {
					md_128  px = md_mm_load_ps(gto->x + i);
					md_128  py = md_mm_load_ps(gto->y + i);
					md_128  pz = md_mm_load_ps(gto->z + i);
					md_128  pa = md_mm_load_ps(gto->neg_alpha + i);
					md_128  pc = md_mm_load_ps(gto->coeff + i);
					md_128i pi = md_mm_load_si128((const md_128i*)(gto->i + i));
					md_128i pj = md_mm_load_si128((const md_128i*)(gto->j + i));
					md_128i pk = md_mm_load_si128((const md_128i*)(gto->k + i));

					md_128 dx = md_mm_sub_ps(px, vx);
					md_128 dy = md_mm_sub_ps(py, vy);
					md_128 dz = md_mm_sub_ps(pz, vz);
					md_128 d2 = md_mm_fmadd_ps(dx, dx, md_mm_fmadd_ps(dy, dy, md_mm_mul_ps(dz, dz)));
					md_128 fx = md_mm_fast_pow(dx, pi);
					md_128 fy = md_mm_fast_pow(dy, pj);
					md_128 fz = md_mm_fast_pow(dz, pk);
					md_128 ex = md_mm_exp_ps(md_mm_mul_ps(pa, d2));

					md_128 prod = md_mm_mul_ps(md_mm_mul_ps(md_mm_mul_ps(pc, fx), md_mm_mul_ps(fy, fz)), ex);

					vpsi = md_mm_add_ps(vpsi, prod);
				}

				psi = md_mm_reduce_add_ps(vpsi);

				int index = x_stride + y_stride + z_stride;
				grid_data[index] = (float)psi;
			}
		}
	}
}

void md_gto_grid_evaluate_sub(md_grid_t* grid, const int grid_idx_min[3], const int grid_idx_max[3], const md_gto_data_t* in_gto) {
	ASSERT(grid);
	ASSERT(in_gto);

	size_t temp_pos = md_temp_get_pos();
	size_t stride = ALIGN_TO(in_gto->count, 16);
	size_t bytes  = stride * (sizeof(float) * 5 + sizeof(int) * 3);
	void* mem = md_temp_push_aligned(bytes, 64);
	md_gto_data_t gto = {
		.count	   = 0,
		.x		   = (float*)mem,
		.y		   = (float*)mem + stride * 1,
		.z		   = (float*)mem + stride * 2,
		.neg_alpha = (float*)mem + stride * 3,
		.coeff	   = (float*)mem + stride * 4,
		.i		   = (int*)  mem + stride * 5,
		.j		   = (int*)  mem + stride * 6,
		.k		   = (int*)  mem + stride * 7,
	};

	for (size_t i = in_gto->count; i < stride; ++i) {
		gto.coeff[i] = 0.0f;
	}

	// Compute the extent of the spatial region which is to be evaluated
	// This will be used to limit the number of pgtos to a set which only has a valid contribution in this region
	float min_box[3] = {
		grid->origin[0] + grid_idx_min[0] * grid->stepsize[0],
		grid->origin[1] + grid_idx_min[1] * grid->stepsize[1],
		grid->origin[2] + grid_idx_min[2] * grid->stepsize[2],
	};
	float max_box[3] = {
		grid->origin[0] + grid_idx_max[0] * grid->stepsize[0],
		grid->origin[1] + grid_idx_max[1] * grid->stepsize[1],
		grid->origin[2] + grid_idx_max[2] * grid->stepsize[2],
	};

	// Extract a subset of gtos that overlap with the evaluated subportion of the grid
	// @TODO: This can be vectorized as well
	for (size_t i = 0; i < in_gto->count; ++i) {
		float x  = in_gto->x[i];
		float y  = in_gto->y[i];
		float z  = in_gto->z[i];
		float cutoff = in_gto->cutoff[i];

		float cx = CLAMP(x, min_box[0], max_box[0]);
		float cy = CLAMP(y, min_box[1], max_box[1]);
		float cz = CLAMP(z, min_box[2], max_box[2]);

		float dx = x - cx;
		float dy = y - cy;
		float dz = z - cz;

		float d2 = dx * dx + dy * dy + dz * dz;
			
		if (d2 > cutoff * cutoff) {
			continue;
		}

		size_t idx = gto.count;
		gto.x[idx]			= x;
		gto.y[idx]			= y;
		gto.z[idx]			= z;
		gto.coeff[idx]		= in_gto->coeff[i];
		gto.neg_alpha[idx]	= in_gto->neg_alpha[i];
		gto.i[idx]			= in_gto->i[i];
		gto.j[idx]			= in_gto->j[i];
		gto.k[idx]			= in_gto->k[i];

		gto.count += 1;
	}

	// This is to allow for zero contributions in the vectorized paths evaluates chunks of 4, 8 or 16 pgtos
	for (size_t i = gto.count; i < ROUND_UP(gto.count, 16); ++i) {
		gto.coeff[i] = 0;
	}

	//printf("Number of pgtos in volume region: %zu\n", gto.count);

#if defined(__AVX512F__)
	evaluate_grid_512(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->stepsize, &gto);
#elif defined(__AVX2__)
	evaluate_grid_256(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->stepsize, &gto);
#elif defined(__SSE2__)
	evaluate_grid_128(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->stepsize, &gto);
#else
	evaluate_grid_ref(grid->data, grid_idx_min, grid_idx_max, grid->dim, grid->origin, grid->stepsize, &gto);
#endif

	md_temp_set_pos_back(temp_pos);
}

void md_gto_grid_evaluate(md_grid_t* grid, const md_gto_data_t* gto) {
	int grid_idx_min[3] = {0,0,0};

#if defined(__AVX512F__)
	evaluate_grid_512(grid->data, grid_idx_min, grid->dim, grid->dim, grid->origin, grid->stepsize, gto);
#elif defined(__AVX2__)
	evaluate_grid_256(grid->data, grid_idx_min, grid->dim, grid->dim, grid->origin, grid->stepsize, gto);
#elif defined(__SSE2__)
	evaluate_grid_128(grid->data, grid_idx_min, grid->dim, grid->dim, grid->origin, grid->stepsize, gto);
#else
	evaluate_grid_ref(grid->data, grid_idx_min, grid->dim, grid->dim, grid->origin, grid->stepsize, gto);
#endif
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

static double compute_distance_cutoff(double cutoff_value, int i, int j, int k, int l, double neg_alpha, double coeff) {
	double d = 0.0;

	// Bake into single constant C
	const double C = fabs(coeff * sqrt((fast_pow(i,i) * fast_pow(j,j) * fast_pow(k,k)) / fast_pow(l,l)));

	// Compute maxima
	const double d_maxima = sqrt(l / (2.0 * fabs(neg_alpha)));

	// Check the contribution at the maxima
	if (eval_G(d_maxima, C, l, neg_alpha) < cutoff_value) {
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

static double compute_distance_cutoff_v2(double cutoff_value, int i, int j, int k, int l, double neg_alpha, double coeff) {
	double min_d = 0.0001;
	double max_d = 100.0;
	double val	 = 0.0;
	bool decrementing = true;

	// Bake into single constant C
	const double C = fabs(coeff * sqrt((fast_pow(i,i) * fast_pow(j,j) * fast_pow(k,k)) / fast_pow(l,l)));

	// Compute maxima
	const double d_maxima = sqrt(l / (2.0 * fabs(neg_alpha)));

	double d = 0.0;
	if (eval_G(d_maxima, C, l, neg_alpha) < cutoff_value) {
		goto done;
	}

	// Initial guess
	d = d_maxima + 0.5;

	while (true) {
		if (d < min_d){
			d = min_d;
			break;
		}
		if (d > max_d){
			d = max_d;
			break;
		}
		val = eval_G(d, C, l, neg_alpha);
		if (decrementing) {
			d = d * 0.5;
			if (val > cutoff_value) {
				decrementing = false;
			}
		} else {
			d = d * 2.0;
			if (val < cutoff_value) {
				d = d * 0.5;
				break;
			}
		}
	}

done:
#if PRINT_RESULT
	if (d > 0.0) {
		printf("Cutoff dist and value: %15.5f, %15.12f\n", d, eval_G(d, C, l, neg_alpha));
	}
#endif
	return d;
}

static double compute_distance_cutoff_v1(double cutoff_value, int i, int j, int k, int l, double neg_alpha, double coeff) {
	(void)l;
	int maxijk = MAX(i, MAX(j,k));
	double d	 = 1.0;
	double min_d = 0.0001;
	double max_d = 100.0;
	double val	 = 0.0;
	bool decrementing = true;

	double f1 = fast_pow(d, maxijk);
	double f3 = f1 * f1 * f1;
	double f  = MAX(f1, f3);
	val = fabs(coeff * f * exp(neg_alpha * d * d));

	if (val < cutoff_value) {
		d = 0.0;
		goto done;
	}

	while (true) {
		if (d < min_d){
			d = min_d;
			break;
		}
		if (d > max_d){
			d = max_d;
			break;
		}
		f1 = fast_pow(d, maxijk);
		f3 = f1 * f1 * f1;
		f  = MAX(f1, f3);
		val = fabs(coeff * f * exp(neg_alpha * d * d));
		if (decrementing) {
			d = d * 0.5;
			if (val > cutoff_value) {
				decrementing = false;
			}
		} else {
			d = d * 2.0;
			if (val < cutoff_value) {
				d = d * 0.5;
				break;
			}
		}
	}

done:

#if PRINT_RESULT
	if (d > 0.0) {
		f1 = fast_pow(d, maxijk);
		f3 = f1 * f1 * f1;
		f  = MAX(f1, f3);
		val = fabs(coeff * f * exp(neg_alpha * d * d));
		printf("Cutoff dist and value: %15.5f, %15.12f\n", d, val);
	}
#endif

	return d;
}

void md_gto_compute_cutoff(md_gto_data_t* gto, double value) {
	for (size_t i = 0; i < gto->count; ++i) {
		gto->cutoff[i] = (float)compute_distance_cutoff(value, gto->i[i], gto->j[i], gto->k[i], gto->l[i], gto->neg_alpha[i], gto->coeff[i]);
	}
}

